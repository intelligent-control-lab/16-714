#!/usr/bin/env python3
"""Student scaffold for 16-714 HW4, Question 1.5: LQR feedback gain.

Keep hw4_helpers.py and hw4_inputs.py beside this script. Solvers, plotting,
and file output are provided; the functions below implement the exercise.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from scipy.linalg import solve_discrete_are

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import hw4_helpers as helpers

OUTPUT_DIR = HERE / "results" / "student_problem_1_5"


# %% Your implementation
def lqr_feedback(a, b, q, r):
    """Return the state cost matrix P, feedback gain K, and closed-loop matrix A_cl.

    Use the control convention u = -K x."""
    # TODO(student, Problem 1.5): Implement lqr_feedback using your derivation.
    raise NotImplementedError("Implement Problem 1.5 lqr_feedback")


# %% Provided experiment and output

def solve():
    return lqr_feedback(helpers.A, helpers.B, helpers.Q, helpers.R)


def main():
    args = helpers.parser(__doc__, OUTPUT_DIR).parse_args()
    p, k, a_cl = solve()
    out = helpers.output_directory(args.output_dir)
    np.savez_compressed(out / "lqr_gain.npz", P=p, K=k, A_cl=a_cl)
    for name, value in [("P", p), ("K", k), ("A_cl", a_cl)]:
        helpers.save_csv(out / f"{name}.csv", ("position", "velocity"), value)
    helpers.save_json(out / "summary.json", {"problem": "1.5", "P": p.tolist(),
                                            "K": k.tolist(), "A_cl": a_cl.tolist()})
    print(f"K = {k}")
    print(f"Saved Question 1.5 results to {out}")


if __name__ == "__main__":
    main()
