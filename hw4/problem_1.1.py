#!/usr/bin/env python3
"""Student scaffold for 16-714 HW4, Question 1.1: disturbed MPC execution.

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

OUTPUT_DIR = HERE / "results" / "student_problem_1_1"


# %% Your implementation
def disturbed_step(state, control, disturbance, a, b):
    """Return the next physical state for the disturbed plant."""
    # TODO(student, Problem 1.1): Implement disturbed_step using your derivation.
    raise NotImplementedError("Implement Problem 1.1 disturbed_step")


# %% Provided experiment and output

def solve(seed=None):
    return helpers.disturbed_rollout(disturbed_step, helpers.disturbance_sequence(seed))


def main():
    parser = helpers.parser(__doc__, OUTPUT_DIR)
    helpers.add_seed_argument(parser)
    args = parser.parse_args()
    result = solve(args.seed)
    out = helpers.output_directory(args.output_dir)
    helpers.save_trial(result, out)
    helpers.save_json(out / "summary.json", {
        "problem": "1.1", "seed": args.seed, "reference_inputs": args.seed is None,
        "executed_steps": len(result["controls"]),
        "max_abs_control": float(np.max(abs(result["controls"]))),
        "max_nominal_constraint_violation": float(result["solver_diagnostics"][:, 0].max())})
    print(f"Saved Question 1.1 results to {out}")


if __name__ == "__main__":
    main()
