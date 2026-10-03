#!/usr/bin/env python3
"""Student scaffold for 16-714 HW3, Question 2.4: infinite-horizon LQR with a state-control cross term.

Keep hw3_helpers.py beside this script. Simulation, solver calls, plots, and
artifact I/O are provided; the mathematical functions below are the exercise.
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

OUTPUT_DIR = HERE / "results" / "student_problem_2_4"


# %% Question 2.4 uses T=1 and R=1, unlike the MPC questions.
A = np.array([[1.0, 1.0], [0.0, 1.0]])
B = np.array([[0.5], [1.0]])
Q = np.eye(2)
R = np.array([[1.0]])
CROSS_TERM = np.array([[0.1], [0.01]])  # N in the written question

# %% Mathematical implementation
def riccati_step(p, a, b, q, r, cross_term):
    """Return one generalized DARE fixed-point update, including the cross term.

    The stage cost is (x.T Q x + 2 x.T N u + u.T R u)/2.
    Use a linear solve, not an explicit inverse; return a symmetric matrix.
    """
    # TODO(student, Problem 2.4): Implement the Riccati equation derived in Question 2.3.
    raise NotImplementedError("Implement Problem 2.4 riccati_step")


def feedback_gain(p, a, b, r, cross_term):
    """Return K for the convention u=-Kx, using the converged P."""
    # TODO(student, Problem 2.4): Implement the control gain derived in Question 2.2.
    raise NotImplementedError("Implement Problem 2.4 feedback_gain")


# %% Provided fixed-point driver
def solve(tolerance=1e-12, max_iterations=10000):
    p = np.zeros_like(Q)
    history = []
    for iteration in range(1, max_iterations + 1):
        updated = riccati_step(p, A, B, Q, R, CROSS_TERM)
        change = float(np.linalg.norm(updated - p, ord=np.inf))
        if not np.all(np.isfinite(updated)):
            raise RuntimeError("Riccati iteration produced nonfinite values.")
        history.append([iteration, change])
        p = updated
        if change < tolerance:
            return p, feedback_gain(p, A, B, R, CROSS_TERM), np.asarray(history)
    raise RuntimeError("Riccati iteration did not converge.")


def main():
    args = helpers.parser(__doc__, OUTPUT_DIR).parse_args()
    p, gain, history = solve()
    poles = np.linalg.eigvals(A - B @ gain)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    helpers.save_csv(out / "P.csv", ("P_column_0", "P_column_1"), p)
    helpers.save_csv(out / "K.csv", ("position_gain", "velocity_gain"), gain)
    helpers.save_csv(out / "iterations.csv", ("iteration", "update_inf_norm"), history)
    np.savez_compressed(out / "lqr_solution.npz", P=p, K=gain, poles=poles,
                        A=A, B=B, Q=Q, R=R, cross_term=CROSS_TERM)
    helpers.save_json(out / "summary.json", {
        "problem": "2.4", "method": "generalized Riccati fixed-point iteration",
        "P": p.tolist(), "K": gain.tolist(), "control_convention": "u = -K x",
        "iterations": len(history),
        "dare_residual_inf": float(np.linalg.norm(riccati_step(p,A,B,Q,R,CROSS_TERM)-p, ord=np.inf)),
        "closed_loop_poles": [{"real": float(z.real), "imag": float(z.imag)} for z in poles],
        "spectral_radius": float(max(abs(poles)))})
    print(f"Question 2.4:\nK = {gain}\nP = {p}")
    print(f"Saved results to {out}")


if __name__ == "__main__":
    main()
