#!/usr/bin/env python3
"""Student scaffold for 16-714 HW4, Question 1.2: tracking error.

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

OUTPUT_DIR = HERE / "results" / "student_problem_1_2"
INPUT_DIR = HERE / "results" / "student_problem_1_1"


# %% Your implementation
def trajectory_error(nominal_states, actual_states):
    """Return nominal minus executed states for all times 0 through 100.

    Each row contains position and velocity. The first row is zero; the
    written stack in Question 1.2 omits that row, while Question 1.4 retains it."""
    # TODO(student, Problem 1.2): Implement trajectory_error using your derivation.
    raise NotImplementedError("Implement Problem 1.2 trajectory_error")


# %% Provided experiment and output

def solve(input_dir=INPUT_DIR):
    trial = helpers.load_results(input_dir, "trial.npz")
    return trajectory_error(trial["nominal_states"], trial["states"])


def main():
    args = helpers.parser(__doc__, OUTPUT_DIR, INPUT_DIR).parse_args()
    error = np.asarray(solve(args.input_dir))
    if error.shape != (helpers.HORIZON+1, 2):
        raise ValueError("Return all 101 error samples with shape (101,2).")
    out = helpers.output_directory(args.output_dir)
    np.savez_compressed(out / "tracking_error.npz", errors=error, lifted_error=error[1:].reshape(-1))
    helpers.save_csv(out / "error.csv", ("k", "position_error", "velocity_error"),
                     np.column_stack([np.arange(len(error)), error]))
    helpers.plot_error(error, out)
    helpers.save_json(out / "summary.json", {"problem": "1.2", "error_l2": float(np.linalg.norm(error)),
                                            "max_abs_error": float(np.max(abs(error)))})
    print(f"Saved Question 1.2 results to {out}")


if __name__ == "__main__":
    main()
