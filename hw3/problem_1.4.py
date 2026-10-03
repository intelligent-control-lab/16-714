#!/usr/bin/env python3
"""Student scaffold for 16-714 HW3, Question 1.4: constrained receding-horizon MPC.

Keep hw3_helpers.py beside this script. Simulation, solver calls, plots, and
artifact I/O are provided; the mathematical functions below are the exercise.
Use the QP from Question 1.1 to track the reference specified below. Plot all
predictions and the executed trajectory together, as requested in Question 1.4.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

# Provided sibling import; works when this script is run from any directory.
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import hw3_helpers as helpers

OUTPUT_DIR = HERE / "results" / "student_problem_1_4"


# %% Given parameters and experiment from Question 1.4
DT = 0.1
A = np.array([[1.0, DT], [0.0, 1.0]])
B = np.array([[0.5 * DT**2], [DT]])
X0 = np.array([2.0, 0.0])
Q = np.eye(2)
R = np.array([[0.1]])
S = 10.0 * np.eye(2)
STATE_BOUND = np.array([5.0, 5.0])
INPUT_BOUND = 1.0
HORIZON = 100
EXECUTED_STEPS = HORIZON


# %% Mathematical implementation
def lifted_dynamics(a, b, horizon):
    """Return F, G such that stacked states X = F @ x0 + G @ U.

    Stack x0 through x_N and u0 through u_(N-1), in time order.
    For n states and m inputs, return shapes ((N+1)*n,n) and ((N+1)*n,N*m).
    You may reuse your own implementation from Question 1.2 here.
    """
    # TODO(student, Problem 1.4): Implement the lifted dynamics derived in Question 1.1.
    raise NotImplementedError("Implement Problem 1.4 lifted_dynamics")


def cost_matrices(q, r, s, horizon):
    """Return the stacked state-error and control cost matrices Qbar, Rbar.

    Include the terminal cost from the handout. With n states and m inputs,
    the matrices have shapes ((N+1)*n,(N+1)*n) and (N*m,N*m).
    """
    # TODO(student, Problem 1.4): Implement the cost matrices derived in Question 1.1.
    raise NotImplementedError("Implement Problem 1.4 cost_matrices")


def condensed_qp_matrices(g, qbar, rbar):
    """Return H, C, A_ineq for the fixed part of the condensed QP.

    The solver uses 0.5 U.T H U + gradient.T U and A_ineq U <= bounds.
    C is the constant gradient map passed to online_qp_terms. Choose a
    consistent constraint row order in this function and online_qp_terms.
    """
    # TODO(student, Problem 1.4): Form the condensed quadratic cost and state/input constraint matrix.
    raise NotImplementedError("Implement Problem 1.4 condensed_qp_matrices")


def online_qp_terms(state, reference, f, gradient_map, state_bound, input_bound):
    """Return gradient and inequality bounds for the current MPC step.

    At step k, reference contains r_k through r_(k+N), with shape (N+1,2).
    Flatten in ordinary C order. Bound rows must match condensed_qp_matrices.
    State bounds include the current x_k and the predicted terminal x_(k+N).
    """
    # TODO(student, Problem 1.4): Insert the current state and moving reference window into the QP.
    raise NotImplementedError("Implement Problem 1.4 online_qp_terms")


# %% Provided solver, receding-horizon rollout, plotting, and saving
def solve():
    # Given reference: r_k = (N-k)/N * x0 up to N, then zero through 2N.
    times = np.arange(EXECUTED_STEPS + HORIZON + 1)
    reference = np.maximum(1 - times / HORIZON, 0)[:, None] * X0
    # The driver solves at k=0,...,N-1, keeping a horizon of N each time.
    # It applies only the first control and records every predicted trajectory.
    return helpers.mpc_rollout(
        condensed_qp_matrices, online_qp_terms, lifted_dynamics, cost_matrices,
        a=A, b=B, q=Q, r=R, s=S, initial_state=X0, reference=reference,
        state_bound=STATE_BOUND, input_bound=INPUT_BOUND,
        horizon=HORIZON, steps=EXECUTED_STEPS)


def main():
    args = helpers.parser(__doc__, OUTPUT_DIR).parse_args()
    result = solve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    for key, names in [("states", ("position", "velocity")),
                       ("controls", ("acceleration",)),
                       ("reference", ("position", "velocity"))]:
        values = result[key]
        helpers.save_csv(out / (key + ".csv"), ("k",) + names,
                         np.column_stack([np.arange(len(values)), values]))
    helpers.save_csv(out / "solver_diagnostics.csv", ("k", "constraint_violation", "iterations"),
                     np.column_stack([np.arange(len(result["controls"])), result["solver_diagnostics"]]))
    np.savez_compressed(out / "mpc_trajectories.npz", **result)
    helpers.plot_mpc(result, out, INPUT_BOUND)
    helpers.save_json(out / "summary.json", {
        "problem": "1.4", "method": "receding-horizon condensed QP (OSQP)",
        "horizon": HORIZON, "executed_steps": len(result["controls"]),
        "state_order": ["position", "velocity"],
        "final_state": result["states"][-1].tolist(),
        "minimum_velocity": float(result["states"][:, 1].min()),
        "max_abs_control": float(np.max(abs(result["controls"]))),
        "max_constraint_violation": float(result["solver_diagnostics"][:, 0].max()),
        "terminal_invariant_constraint": False, "figure": "mpc_trajectories.png"})
    print(f"Question 1.4: final state {result['states'][-1]}")
    print(f"Saved results to {out}")


if __name__ == "__main__":
    main()
