#!/usr/bin/env python3
"""Student scaffold for 16-714 HW4, Question 2.2: augmented-state feedback gain.

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

OUTPUT_DIR = HERE / "results" / "student_problem_2_2"
INPUT_DIR = HERE / "results" / "student_problem_1_5"


# %% Your implementation
def cross_cost_block(a, b, p, k, r):
    """Return P12 with shape (2,1), using your block equations from Question 2.2.

    P and K are your saved state cost matrix and gain from Question 1.5."""
    # TODO(student, Problem 2.2): Implement cross_cost_block using your derivation.
    raise NotImplementedError("Implement Problem 2.2 cross_cost_block")


def augmented_gain(a, b, p, p12, r):
    """Return the feedback row for [position, velocity, disturbance], with u = -K_e x_e."""
    # TODO(student, Problem 2.2): Implement augmented_gain using your derivation.
    raise NotImplementedError("Implement Problem 2.2 augmented_gain")


# %% Provided experiment and output

def solve(input_dir=INPUT_DIR):
    earlier = helpers.load_results(input_dir, "lqr_gain.npz")
    p, k = earlier["P"], earlier["K"]
    p12 = cross_cost_block(helpers.A, helpers.B, p, k, helpers.R)
    gain = augmented_gain(helpers.A, helpers.B, p, p12, helpers.R)
    return p, p12, gain


def main():
    args = helpers.parser(__doc__, OUTPUT_DIR, INPUT_DIR).parse_args()
    p, p12, gain = solve(args.input_dir)
    out = helpers.output_directory(args.output_dir)
    np.savez_compressed(out / "control_gain.npz", P11=p, P12=p12, K=gain)
    helpers.save_csv(out / "P11.csv", ("position", "velocity"), p)
    helpers.save_csv(out / "P12.csv", ("P12",), p12)
    helpers.save_csv(out / "K.csv", ("position", "velocity", "disturbance"), gain)
    helpers.save_json(out / "summary.json", {"problem": "2.2", "P11": p.tolist(),
                                            "P12": p12.tolist(), "K": gain.tolist()})
    print(f"Augmented gain = {gain}")
    print(f"Saved Question 2.2 results to {out}")


if __name__ == "__main__":
    main()
