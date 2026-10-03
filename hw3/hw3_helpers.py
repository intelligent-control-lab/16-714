"""Provided numerical, geometry, plotting, and output utilities for HW3.

Exercise-specific matrices and control laws are supplied by callbacks from
the problem scripts. Keep this file beside them. States use time-major order.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import osqp
from scipy import sparse
from scipy.optimize import linprog

DT = 0.1
A = np.array([[1.0, DT], [0.0, 1.0]])
B = np.array([[0.5 * DT**2], [DT]])
STATE_BOUND = np.array([5.0, 5.0])
INPUT_BOUND = 1.0


def parser(description, output_dir):
    result = argparse.ArgumentParser(description=description)
    result.add_argument("--output-dir", type=Path, default=output_dir)
    return result


def save_json(path, values):
    Path(path).write_text(json.dumps(values, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def save_csv(path, header, values):
    np.savetxt(path, values, delimiter=",", header=",".join(header), comments="")


def save_figure(fig, directory, name):
    for suffix in ("png", "pdf"):
        fig.savefig(Path(directory) / (name + "." + suffix), dpi=180, bbox_inches="tight")
    plt.close(fig)


def lp_feasible(matrix, bound, tolerance=1e-7):
    """Test feasibility; distinguish infeasibility from solver failure.

    All variable bounds must be included in matrix/bound. SciPy's default
    nonnegative-variable bounds would incorrectly exclude negative controls.
    """
    result = linprog(np.zeros(matrix.shape[1]), A_ub=matrix, b_ub=bound,
                     bounds=(None, None), method="highs")
    if result.status == 2:
        return False
    if not result.success:
        raise RuntimeError(f"LP failed (status {result.status}): {result.message}")
    if np.max(np.einsum("ij,j->i", matrix, result.x) - bound) > tolerance:
        raise RuntimeError("LP returned a point outside the constraint tolerance.")
    return True


def feasibility_grid(constraint_function, dynamics_function, horizon, resolution=0.1):
    if resolution <= 0:
        raise ValueError("Resolution must be positive.")
    count = int(round(10.0 / resolution))
    if count < 1 or not np.isclose(count * resolution, 10.0):
        raise ValueError("Resolution must be positive and divide the 10-unit state interval.")
    grid = np.linspace(-5.0, 5.0, count + 1)
    f, g = dynamics_function(A, B, horizon)
    mask = np.empty((len(grid), len(grid)), dtype=bool)
    # The matrix is constant. Evaluate the student's formula at zero once,
    # then use its bounds at each sampled initial state.
    matrix, _ = constraint_function(np.zeros(2), f, g, STATE_BOUND, INPUT_BOUND)
    for row, velocity in enumerate(grid):
        for column, position in enumerate(grid):
            _, bound = constraint_function(np.array([position, velocity]), f, g,
                                           STATE_BOUND, INPUT_BOUND)
            mask[row, column] = lp_feasible(matrix, bound)
    return grid, mask


def plot_feasible_sets(grid, masks, directory):
    cmap = ListedColormap(["#3d26a8", "#f9fa14"])
    fig, axes = plt.subplots(1, len(masks), figsize=(4.1 * len(masks), 4), squeeze=False)
    step = grid[1] - grid[0]
    for axis, (horizon, mask) in zip(axes[0], masks.items()):
        axis.imshow(mask, origin="lower", extent=[grid[0]-step/2, grid[-1]+step/2]*2,
                    interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
        axis.set(xlabel="Position", ylabel="Velocity", title=f"N = {horizon} ({mask.sum():,} feasible)")
        axis.set_xlim(-5, 5)
        axis.set_ylim(-5, 5)
    fig.suptitle("Question 1.2: feasible initial states (yellow)")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save_figure(fig, directory, "feasible_sets")


def eliminate_scalar(matrix, bound, tolerance=1e-12):
    """Provided Fourier–Motzkin projection of the LAST scalar coordinate.

    Input inequalities are matrix @ [x; u] <= bound. Output inequalities
    describe the x for which at least one u exists. This generic helper does
    not construct a dynamics predecessor or add state/input constraints.
    """
    coefficient = matrix[:, -1]
    positive = np.flatnonzero(coefficient > tolerance)
    negative = np.flatnonzero(coefficient < -tolerance)
    zero = np.flatnonzero(abs(coefficient) <= tolerance)
    rows = list(matrix[zero, :-1])
    rhs = list(bound[zero])
    for upper in positive:
        for lower in negative:
            rows.append(matrix[upper, :-1] / coefficient[upper]
                        - matrix[lower, :-1] / coefficient[lower])
            rhs.append(bound[upper] / coefficient[upper] - bound[lower] / coefficient[lower])
    return np.asarray(rows).reshape(-1, matrix.shape[1]-1), np.asarray(rhs)


def box_halfspaces(bound=STATE_BOUND):
    return np.vstack([np.eye(2), -np.eye(2)]), np.tile(bound, 2)


def polygon_from_halfspaces(matrix, bound, domain=STATE_BOUND, tolerance=1e-10):
    """Clip a bounded 2-D polygon; the stated domain box is always included."""
    px, py = domain
    vertices = np.array([[-px, -py], [px, -py], [px, py], [-px, py]])
    for normal, limit in zip(matrix, bound):
        norm = np.linalg.norm(normal)
        if norm < 1e-13:
            if limit < -tolerance:
                raise ValueError("Empty polytope.")
            continue
        normal, limit = normal / norm, limit / norm
        distance = vertices @ normal - limit
        if np.all(distance <= tolerance):
            continue
        clipped = []
        for i in range(len(vertices)):
            j = (i + 1) % len(vertices)
            inside_i, inside_j = distance[i] <= tolerance, distance[j] <= tolerance
            if inside_i:
                clipped.append(vertices[i])
            if inside_i != inside_j:
                fraction = distance[i] / (distance[i] - distance[j])
                clipped.append(vertices[i] + fraction * (vertices[j] - vertices[i]))
        if len(clipped) < 3:
            raise ValueError("Expected a nonempty, full-dimensional invariant polygon.")
        vertices = np.asarray(clipped)
    # Remove duplicate and collinear vertices introduced by clipping.
    changed = True
    while changed and len(vertices) > 3:
        changed = False
        for i in range(len(vertices)):
            before, after = vertices[i] - vertices[i-1], vertices[(i+1) % len(vertices)] - vertices[i]
            length = np.linalg.norm(before) * np.linalg.norm(after)
            cross = before[0] * after[1] - before[1] * after[0]
            if length < 1e-18 or abs(cross) <= 1e-10 * length:
                vertices = np.delete(vertices, i, axis=0)
                changed = True
                break
    return vertices


def polygon_halfspaces(vertices):
    edges = np.roll(vertices, -1, axis=0) - vertices
    normals = np.column_stack([edges[:, 1], -edges[:, 0]])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    return normals, np.einsum("ij,ij->i", normals, vertices)


def polygon_area(vertices):
    return float(abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1))
                     - np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1))) / 2)


def invariant_iteration(predecessor_function, tolerance=1e-9, max_iterations=200):
    """Provided driver for C_{j+1} = X intersect Pre(C_j), C_0 = X."""
    matrix, bound = box_halfspaces()
    vertices = polygon_from_halfspaces(matrix, bound)
    history = [[0, len(vertices), polygon_area(vertices)]]
    for iteration in range(1, max_iterations + 1):
        candidate_matrix, candidate_bound = predecessor_function(matrix, bound, A, B, STATE_BOUND, INPUT_BOUND)
        new_vertices = polygon_from_halfspaces(candidate_matrix, candidate_bound)
        new_matrix, new_bound = polygon_halfspaces(new_vertices)
        old_in_new = float(np.max(vertices @ new_matrix.T - new_bound))
        new_in_old = float(np.max(new_vertices @ matrix.T - bound))
        history.append([iteration, len(new_vertices), polygon_area(new_vertices)])
        matrix, bound, vertices = new_matrix, new_bound, new_vertices
        if max(old_in_new, new_in_old) <= tolerance:
            return matrix, bound, vertices, np.asarray(history), max(0.0, old_in_new, new_in_old)
    raise RuntimeError("Controlled-predecessor iteration did not reach a fixed point.")


def plot_invariant(vertices, directory):
    fig, axis = plt.subplots(figsize=(6, 5))
    closed = np.vstack([vertices, vertices[0]])
    axis.fill(closed[:, 0], closed[:, 1], color="#dcebc1", label="Maximal controlled invariant set")
    axis.plot(closed[:, 0], closed[:, 1], color="#346320")
    axis.set(xlabel="Position", ylabel="Velocity", xlim=(-5, 5), ylim=(-5, 5),
             title="Question 1.3: persistent feasibility", aspect="equal")
    axis.grid(alpha=0.25)
    axis.legend(loc="upper right", fontsize=8)
    save_figure(fig, directory, "invariant_set")


class QPSolver:
    """Provided OSQP adapter; H and constraint matrix remain fixed across solves."""
    def __init__(self, hessian, matrix):
        self.matrix = np.asarray(matrix)
        self.solver = osqp.OSQP()
        self.solver.setup(P=sparse.csc_matrix((hessian + hessian.T) / 2),
                          q=np.zeros(hessian.shape[0]), A=sparse.csc_matrix(matrix),
                          l=np.full(matrix.shape[0], -np.inf), u=np.zeros(matrix.shape[0]),
                          verbose=False, eps_abs=1e-9, eps_rel=1e-9, max_iter=20000)

    def solve(self, gradient, upper):
        self.solver.update(q=gradient, u=upper)
        result = self.solver.solve()
        if result.info.status_val != 1:
            raise RuntimeError(f"MPC QP failed: {result.info.status}")
        violation = float(max(0.0, np.max(self.matrix @ result.x - upper)))
        if not np.all(np.isfinite(result.x)) or violation > 1e-7:
            raise RuntimeError("MPC solution violates numerical/constraint checks.")
        return result.x.copy(), violation, result.info.iter


def mpc_rollout(matrix_function, online_function, dynamics_function, cost_function,
                *, a, b, q, r, s, initial_state, reference, state_bound,
                input_bound, horizon, steps):
    """Run the supplied two-state, one-input MPC experiment.

    The problem script supplies all parameters and the reference; callbacks
    supply the assessed mathematics. Use a fixed prediction horizon at every
    step, applying only the first control from each newly optimized plan.
    """
    f, g = dynamics_function(a, b, horizon)
    qbar, rbar = cost_function(q, r, s, horizon)
    hessian, gradient_map, constraint_matrix = matrix_function(g, qbar, rbar)
    solver = QPSolver(hessian, constraint_matrix)
    states, controls = np.zeros((steps+1, 2)), np.zeros((steps, 1))
    predicted_states = np.zeros((steps, horizon+1, 2))
    predicted_controls = np.zeros((steps, horizon, 1))
    diagnostics = np.zeros((steps, 2))
    states[0] = initial_state
    for k in range(steps):
        gradient, upper = online_function(states[k], reference[k:k+horizon+1], f,
                                          gradient_map, state_bound, input_bound)
        try:
            plan, violation, iterations = solver.solve(gradient, upper)
        except RuntimeError as error:
            raise RuntimeError(f"At MPC step {k}: {error}") from error
        predicted_states[k] = (f @ states[k] + g @ plan).reshape(horizon+1, 2)
        predicted_controls[k, :, 0] = plan
        controls[k, 0] = plan[0]
        states[k+1] = a @ states[k] + b[:, 0] * controls[k, 0]
        diagnostics[k] = violation, iterations
    return dict(states=states, controls=controls, reference=reference,
                predicted_states=predicted_states, predicted_controls=predicted_controls,
                solver_diagnostics=diagnostics, F=f, G=g, H=hessian)


def plot_mpc(result, directory, input_bound):
    states, reference = result["states"], result["reference"]
    predictions = result["predicted_states"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True)
    for component, axis in enumerate(axes):
        for k, prediction in enumerate(predictions):
            axis.plot(k + np.arange(len(prediction)), prediction[:, component],
                      color="blue", lw=0.7, alpha=0.65, label="Predicted" if k == 0 else None)
        axis.plot(np.arange(len(states)), states[:, component], color="red", lw=2, label="Executed")
        axis.plot(np.arange(len(reference)), reference[:, component], "k--", lw=1.5, label="Reference")
        axis.set(ylabel=("Position", "Velocity")[component])
        axis.grid(alpha=0.2)
    axes[0].set_title("Question 1.4: receding-horizon MPC")
    axes[0].legend()
    axes[-1].set_xlabel("Time step k")
    fig.tight_layout()
    save_figure(fig, directory, "mpc_trajectories")
    fig, axis = plt.subplots(figsize=(9, 3))
    axis.step(np.arange(len(result["controls"])), result["controls"][:, 0], where="post")
    axis.axhline(input_bound, color="gray", ls="--")
    axis.axhline(-input_bound, color="gray", ls="--")
    axis.set(xlabel="Time step k", ylabel="Control u", title="Executed control and input bounds")
    axis.grid(alpha=0.2)
    save_figure(fig, directory, "mpc_controls")
