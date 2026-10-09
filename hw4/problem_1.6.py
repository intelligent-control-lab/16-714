#!/usr/bin/env python3
"""Student scaffold for 16-714 HW4, Question 1.6: lifted ILC model and learning gain.

Keep hw4_helpers.py and hw4_inputs.py beside this script. Solvers, plotting,
and file output are provided; the functions below implement the exercise.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import hw4_helpers as helpers

OUTPUT_DIR = HERE / "results" / "student_problem_1_6"
INPUT_DIR = HERE / "results" / "student_problem_1_5"


# %% Your implementation
def lifted_error_matrix(a_cl, b, horizon):
    """Return the lifted matrix for positive input through B, including the e_0 rows.

    Stack states in time order with adjacent position and velocity components.
    For horizon 100, return shape (202,100), using your model from Question 1.4."""
    # TODO(student, Problem 1.6): Implement lifted_error_matrix using your derivation.
    raise NotImplementedError("Implement Problem 1.6 lifted_error_matrix")


def learning_gain(lifted_matrix):
    """Return the learning matrix L using your design from Question 1.6.

    For this stacking convention, its shape is (100,202)."""
    # TODO(student, Problem 1.6): Implement learning_gain using your derivation.
    raise NotImplementedError("Implement Problem 1.6 learning_gain")


# %% Provided experiment and output

def solve(input_dir=INPUT_DIR):
    previous = helpers.load_results(input_dir, "lqr_gain.npz")
    matrix = lifted_error_matrix(previous["A_cl"], helpers.B, helpers.HORIZON)
    return matrix, learning_gain(matrix)


def main():
    args = helpers.parser(__doc__, OUTPUT_DIR, INPUT_DIR).parse_args()
    matrix, gain = solve(args.input_dir)
    out = helpers.output_directory(args.output_dir)
    np.savez_compressed(out / "learning_gain.npz", G=matrix, L=gain)
    for name, values in [("position_gain", gain[:, 2::2]), ("velocity_gain", gain[:, 3::2])]:
        helpers.save_csv(out / f"{name}.csv", tuple(f"e_{k}" for k in range(1, helpers.HORIZON+1)), values)
    helpers.plot_learning_gain(gain, out)
    helpers.save_json(out / "summary.json", {"problem": "1.6", "lifted_shape": list(matrix.shape),
                                            "gain_shape": list(gain.shape), "plotted_shape": [100, 100]})
    print(f"Saved Question 1.6 results to {out}")


if __name__ == "__main__":
    main()
