#!/usr/bin/env python3
"""Student scaffold for 16-714 HW4, Question 2.4: steady-state estimation gain.

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

OUTPUT_DIR = HERE / "results" / "student_problem_2_4"


# %% Your implementation
def augmented_model(a, b, c):
    """Return A_e, B_e, C_e, D for the augmented dynamics from Question 2.1.

    Their shapes are (3,3), (3,1), (1,3), and (3,1), respectively."""
    # TODO(student, Problem 2.4): Implement augmented_model using your derivation.
    raise NotImplementedError("Implement Problem 2.4 augmented_model")


def steady_covariance(a, c, d, process_variance, measurement_variance):
    """Return the steady-state prior covariance for the augmented filter.

    Use the covariance equations from Question 2.3."""
    # TODO(student, Problem 2.4): Implement steady_covariance using your derivation.
    raise NotImplementedError("Implement Problem 2.4 steady_covariance")


def kalman_gain(covariance, c, measurement_variance):
    """Return the measurement-update gain as a (3,1) array."""
    # TODO(student, Problem 2.4): Implement kalman_gain using your derivation.
    raise NotImplementedError("Implement Problem 2.4 kalman_gain")


# %% Provided experiment and output

def solve():
    a, b, c, d = augmented_model(helpers.A, helpers.B, helpers.C)
    covariance = steady_covariance(a, c, d, helpers.PROCESS_VARIANCE, helpers.MEASUREMENT_VARIANCE)
    gain = kalman_gain(covariance, c, helpers.MEASUREMENT_VARIANCE)
    return dict(A=a, B=b, C=c, D=d, M=covariance, F=gain)


def main():
    args = helpers.parser(__doc__, OUTPUT_DIR).parse_args()
    model = solve()
    out = helpers.output_directory(args.output_dir)
    np.savez_compressed(out / "kalman_filter.npz", **model)
    helpers.save_csv(out / "M.csv", ("position", "velocity", "disturbance"), model["M"])
    helpers.save_csv(out / "F.csv", ("gain",), model["F"])
    helpers.save_json(out / "summary.json", {"problem": "2.4", "M": model["M"].tolist(),
                                            "F": model["F"].tolist(), "covariance_convention": "prior"})
    print(f"Filter gain = {model['F'].ravel()}")
    print(f"Saved Question 2.4 results to {out}")


if __name__ == "__main__":
    main()
