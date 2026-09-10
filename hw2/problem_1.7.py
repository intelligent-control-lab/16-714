#!/usr/bin/env python3
"""Student scaffold for Question 1.7: one iLQR update.

Implement only delta_stage_step (P,b recursion and control coefficients)
and delta_control_law (the delta-control correction with fixed feedforward
scaling). The model terms, backward driver, nominal-control addition,
nonlinear simulation, co-state evaluation, plotting, and saving are provided.
Use the prescribed step size; no tuning or line search.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# %% Parameters: N transitions, N+1 state/co-state samples.
DT = 0.1
NUM_STATE_SAMPLES = 101
N = NUM_STATE_SAMPLES - 1
STEP_SIZE = 0.1  # Provided numerical setting: fixed feedforward relaxation, no line search.
X0 = np.zeros(4)
X_GOAL = np.array([10.0, 10.0, 0.0, np.pi / 2.0])
Q = np.eye(4)
R = 0.1 * np.eye(2)
S = 10.0 * np.eye(4)
B = np.array([[0.0, 0.0], [0.0, 0.0], [DT, 0.0], [0.0, DT]])
A_GOAL = np.array(
    [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, DT, 0.0],
     [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
)
STATE_NAMES = ("p_x", "p_y", "v", "theta")
CONTROL_NAMES = ("acceleration", "steering_rate")


# %% Provided nonlinear Euler model and the Question 1.6 reference.
# State x=[p_x,p_y,v,theta], control u=[a,omega]; speed is a state.
def nonlinear_step(state: np.ndarray, control: np.ndarray) -> np.ndarray:
    px, py, speed, theta = state
    acceleration, steering_rate = control
    return np.array([
        px + DT * speed * np.cos(theta),
        py + DT * speed * np.sin(theta),
        speed + DT * acceleration,
        theta + DT * steering_rate,
    ])


def state_jacobian(state: np.ndarray) -> np.ndarray:
    speed, theta = state[2], state[3]
    return np.array([
        [1.0, 0.0, DT * np.cos(theta), -DT * speed * np.sin(theta)],
        [0.0, 1.0, DT * np.sin(theta), DT * speed * np.cos(theta)],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])


def problem_1_6_reference() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Question 1.6 nonlinear LQR rollout with approximate co-states."""
    P = np.zeros((N + 1, 4, 4))
    K = np.zeros((N, 2, 4))
    P[N] = S
    for k in range(N - 1, -1, -1):
        K[k] = np.linalg.solve(R + B.T @ P[k + 1] @ B, B.T @ P[k + 1] @ A_GOAL)
        P[k] = Q + A_GOAL.T @ P[k + 1] @ A_GOAL - A_GOAL.T @ P[k + 1] @ B @ K[k]
        P[k] = 0.5 * (P[k] + P[k].T)
    states = np.zeros((N + 1, 4))
    controls = np.zeros((N, 2))
    states[0] = X0
    for k in range(N):
        controls[k] = -K[k] @ (states[k] - X_GOAL)
        states[k + 1] = nonlinear_step(states[k], controls[k])
    costates = np.einsum("kij,kj->ki", P, states - X_GOAL)
    return states, controls, costates


def trajectory_cost(states: np.ndarray, controls: np.ndarray) -> float:
    errors = states[:-1] - X_GOAL
    terminal_error = states[-1] - X_GOAL
    return float(
        0.5 * np.einsum("ni,ij,nj->", errors, Q, errors)
        + 0.5 * np.einsum("ni,ij,nj->", controls, R, controls)
        + 0.5 * terminal_error @ S @ terminal_error
    )


# %% Delta dynamics and affine co-state recursion: lambda_k = P_k delta_x_k + b_k.
def delta_stage_step(
    A_k: np.ndarray, d_k: np.ndarray, q_k: np.ndarray, r_k: np.ndarray,
    P_next: np.ndarray, b_next: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (P_k, b_k, K_k, k_ff_k); shapes (4,4), (4,), (2,4), (2,).

    Use the affine-quadratic recursion. The policy convention is
    delta_u = -K_k delta_x - k_ff_k.
    """
    # TODO(student, Problem 1.7): Compute the Riccati matrix P_k and bias
    # b_k and the two control coefficients K_k and k_ff_k.
    # Keep d_k (the four-state dynamics defect) distinct from k_ff_k.
    # Use np.linalg.solve for the two control-Hessian systems.
    raise NotImplementedError("Implement Problem 1.7 delta_stage_step")


def delta_model_terms(
    nominal_states: np.ndarray, nominal_controls: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return A, d, q, r for the current nominal trajectory.

    Shapes: A (N,4,4), d (N,4), q (N+1,4), r (N,2).
    q has N+1 rows: stage gradients q[0:N] and terminal gradient q[N].
    Delta dynamics: delta_x_next = A_k delta_x + B delta_u + d_k.
    """
    A = np.array([state_jacobian(x) for x in nominal_states[:-1]])
    d = np.array([
        nonlinear_step(nominal_states[k], nominal_controls[k]) - nominal_states[k + 1]
        for k in range(N)
    ])
    # Tracking-cost gradients: q_k=Q(x_bar_k-x_G), q_N=S(x_bar_N-x_G).
    q = np.zeros((N + 1, 4))
    q[:N] = (nominal_states[:-1] - X_GOAL) @ Q.T
    q[N] = S @ (nominal_states[N] - X_GOAL)
    r = nominal_controls @ R.T
    return A, d, q, r


def delta_backward_pass(nominal_states: np.ndarray, nominal_controls: np.ndarray) -> dict:
    """Provided backward driver: P_N=S, b_N=q_N, then k=N-1,...,0."""
    A, d, q, r = delta_model_terms(nominal_states, nominal_controls)
    P = np.zeros((N + 1, 4, 4))
    b = np.zeros((N + 1, 4))
    K = np.zeros((N, 2, 4))
    k_ff = np.zeros((N, 2))
    P[N] = S
    # lambda_N = S delta_x_N + q_N = P_N delta_x_N + b_N.
    b[N] = q[N]
    for k in range(N - 1, -1, -1):
        P[k], b[k], K[k], k_ff[k] = delta_stage_step(
            A[k], d[k], q[k], r[k], P[k + 1], b[k + 1]
        )
    return {"A": A, "B": B.copy(), "d": d, "q": q, "r": r,
            "P": P, "b": b, "K": K, "k_ff": k_ff}


def delta_control_law(
    state: np.ndarray, nominal_state: np.ndarray,
    K_k: np.ndarray, k_ff_k: np.ndarray, step_size: float,
) -> np.ndarray:
    """Return the prescribed delta control -K_k (x-x_bar) - eta k_ff_k.

    state is x_k in the current rollout; nominal_state is the reference x_bar_k.
    Their difference is delta_x_k. The return value is delta_u_k, not u_k.
    Shapes: state and nominal_state (4,), K_k (2,4), k_ff_k and result (2,).
    The passed step_size is eta; scale only the feedforward term.
    rollout_update adds the nominal control and propagates the dynamics.
    """
    # TODO(student, Problem 1.7): Form the state deviation and return the
    # delta-control correction in the question. Use the passed step_size
    # exactly once on k_ff_k; leave K_k unscaled. Do not add a nominal
    # control or simulate the dynamics here; rollout_update handles both.
    raise NotImplementedError("Implement Problem 1.7 delta_control_law")


# %% Provided nonlinear rollout and local co-state evaluation.
def rollout_update(
    nominal_states: np.ndarray, nominal_controls: np.ndarray,
    backward: dict, step_size: float = STEP_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Simulate x_k, u_k while keeping the reference x_bar_k, u_bar_k fixed.

    rollout_states[k] and rollout_controls[k] are x_k and u_k.
    nominal_states[k] and nominal_controls[k] are x_bar_k and u_bar_k.
    Each x_k is known from initialization or the preceding simulation step.
    """
    if not 0.0 < step_size <= 1.0:
        raise ValueError("step_size must lie in (0, 1].")
    rollout_states = np.zeros_like(nominal_states)
    rollout_controls = np.zeros_like(nominal_controls)
    rollout_states[0] = X0
    for k in range(N):
        delta_control = delta_control_law(
            rollout_states[k], nominal_states[k],
            backward["K"][k], backward["k_ff"][k], step_size,
        )
        rollout_controls[k] = nominal_controls[k] + delta_control
        rollout_states[k + 1] = nonlinear_step(rollout_states[k], rollout_controls[k])
    # Local co-state lambda_k = P_k (x_k - x_bar_k) + b_k.
    rollout_costates = (
        np.einsum("kij,kj->ki", backward["P"], rollout_states - nominal_states)
        + backward["b"]
    )
    return rollout_states, rollout_controls, rollout_costates


def nonlinear_costate(states: np.ndarray) -> np.ndarray:
    """Nonlinear adjoint diagnostic: lambda_k=Q(x_k-x_G)+A(x_k).T lambda_next.

    This is distinct from the local approximation P delta_x + b.
    A small stationarity residual is a necessary check, not a proof of global
    optimality.
    """
    costates = np.zeros((N + 1, 4))
    costates[N] = S @ (states[N] - X_GOAL)
    for k in range(N - 1, -1, -1):
        costates[k] = Q @ (states[k] - X_GOAL) + state_jacobian(states[k]).T @ costates[k + 1]
    return costates


def stationarity_residual(controls: np.ndarray, costates: np.ndarray) -> np.ndarray:
    return controls @ R.T + costates[1:] @ B


def save_csv(path: Path, header: tuple[str, ...], values: np.ndarray) -> None:
    np.savetxt(path, values, delimiter=",", header=",".join(header), comments="")


def save_trajectory(output_dir: Path, states: np.ndarray, controls: np.ndarray,
                    costates: np.ndarray, prefix: str = "") -> None:
    for name, names, values in (
        ("states", STATE_NAMES, states), ("controls", CONTROL_NAMES, controls),
        ("costates", tuple("lambda_" + name for name in STATE_NAMES), costates),
    ):
        save_csv(output_dir / (prefix + name + ".csv"), ("k",) + names,
                 np.column_stack((np.arange(len(values)), values)))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-size", type=float, default=STEP_SIZE,
                        help="Constant step throughout the run; 1 uses the full control update.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args(argv)
    if not 0.0 < args.step_size <= 1.0:
        parser.error("--step-size must lie in (0, 1].")
    return args


OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "student_problem_1_7"


def plot_trajectories(states: np.ndarray, controls: np.ndarray, costates: np.ndarray) -> Path:
    output = OUTPUT_DIR / "problem_1_7_next_iteration.png"
    state_steps = np.arange(N + 1)
    control_steps = np.arange(N)
    colors = ("tab:blue", "tab:orange", "#EDB120", "tab:purple")
    fig = plt.figure(figsize=(11.0, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.15))
    axes = (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :]))
    axes[0].plot(control_steps, controls[:, 0], label="acceleration")
    axes[0].plot(control_steps, controls[:, 1], label="steering rate")
    axes[0].set_title("Control")
    for index, (name, color) in enumerate(zip(STATE_NAMES, colors)):
        axes[1].plot(state_steps, costates[:, index], label=name, color=color)
        axes[2].plot(state_steps, states[:, index], label=name, color=color)
    axes[1].set_title("Co-state (local model)")
    axes[2].set_title("State")
    for axis in axes:
        axis.set_xlabel("time step $k$")
        axis.grid(alpha=0.25)
        axis.legend()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_delta(
    nominal_states: np.ndarray,
    nominal_controls: np.ndarray,
    new_states: np.ndarray,
    new_controls: np.ndarray,
) -> Path:
    output = OUTPUT_DIR / "problem_1_7_delta_update.png"
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0), sharex=True, constrained_layout=True)
    state_steps = np.arange(N + 1)
    control_steps = np.arange(N)
    delta_states = new_states - nominal_states
    delta_controls = new_controls - nominal_controls
    for index, name in enumerate(STATE_NAMES):
        axes[0, 0].plot(state_steps, nominal_states[:, index], label=name)
        axes[0, 1].plot(state_steps, new_states[:, index], label=name)
        axes[1, 0].plot(state_steps, delta_states[:, index], label=f"delta {name}")
    axes[1, 1].plot(control_steps, delta_controls[:, 0], label="delta acceleration")
    axes[1, 1].plot(control_steps, delta_controls[:, 1], label="delta steering rate")
    titles = ("Q1.6 nominal state", "next-iteration state", "state correction", "control correction")
    for axis, title in zip(axes.flat, titles):
        axis.set_title(title)
        axis.set_xlabel("time step $k$")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def main(argv=None) -> None:
    global OUTPUT_DIR
    args = parse_args(argv)
    OUTPUT_DIR = args.output_dir.resolve()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    nominal_states, nominal_controls, nominal_costates = problem_1_6_reference()
    backward = delta_backward_pass(nominal_states, nominal_controls)
    states, controls, costates = rollout_update(
        nominal_states, nominal_controls, backward, args.step_size
    )
    nonlinear_costates = nonlinear_costate(states)
    old_cost = trajectory_cost(nominal_states, nominal_controls)
    new_cost = trajectory_cost(states, controls)
    save_trajectory(OUTPUT_DIR, states, controls, costates)
    save_trajectory(OUTPUT_DIR, nominal_states, nominal_controls, nominal_costates, "nominal_")
    save_csv(OUTPUT_DIR / "nonlinear_costates.csv",
             ("k",) + tuple("lambda_" + name for name in STATE_NAMES),
             np.column_stack((np.arange(N + 1), nonlinear_costates)))
    np.savez_compressed(
        OUTPUT_DIR / "delta_backward_pass.npz", **backward,
        nominal_states=nominal_states, nominal_controls=nominal_controls,
        states=states, controls=controls, costates=costates,
        nonlinear_costates=nonlinear_costates,
        delta_states=states - nominal_states, delta_controls=controls - nominal_controls,
        step_size=args.step_size,
    )
    figures = [plot_trajectories(states, controls, costates),
               plot_delta(nominal_states, nominal_controls, states, controls)]
    summary = {
        "problem": "1.7",
        "initialization": "Question 1.6 nonlinear LQR trajectory",
        "method": "affine co-state P,b recursion with a constant step",
        "horizon": N, "state_samples": NUM_STATE_SAMPLES, "dt": DT,
        "step_size": args.step_size,
        "costate_definition": "P_k (x_new_k - x_nominal_k) + b_k; local approximation",
        "nominal_cost": old_cost, "new_cost": new_cost,
        "relative_cost_decrease": (old_cost - new_cost) / max(1.0, abs(old_cost)),
        "nominal_final_state": nominal_states[-1].tolist(),
        "new_final_state": states[-1].tolist(),
        "terminal_error_norm": float(np.linalg.norm(states[-1] - X_GOAL)),
        "maximum_nominal_dynamics_defect": float(np.max(np.abs(backward["d"]))),
        "nonlinear_stationarity_infinity_norm": float(np.max(np.abs(
            stationarity_residual(controls, nonlinear_costates)))),
        "figures": [path.name for path in figures],
    }
    (OUTPUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Question 1.7: constant step {args.step_size:g}, J {old_cost:.9f} -> {new_cost:.9f}")
    print(f"final state = {np.array2string(states[-1], precision=9)}")
    print(f"saved results to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
