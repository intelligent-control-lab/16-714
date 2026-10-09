#!/usr/bin/env python3
"""Student scaffold for 16-714 HW4, Question 1.7: repeated ILC trials.

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

OUTPUT_DIR = HERE / "results" / "student_problem_1_7"
INPUT_DIR = HERE / "results" / "student_problem_1_6"
TRIAL_DIR = HERE / "results" / "student_problem_1_1"


# %% Your implementation
def disturbed_step(state, control, disturbance, a, b):
    """Reuse your own disturbed-plant implementation from Question 1.1."""
    # TODO(student, Problem 1.7): Implement disturbed_step using your derivation.
    raise NotImplementedError("Implement Problem 1.7 disturbed_step")


def trajectory_error(nominal_states, actual_states):
    """Reuse your own error calculation from Question 1.2, including the zero initial row."""
    # TODO(student, Problem 1.7): Implement trajectory_error using your derivation.
    raise NotImplementedError("Implement Problem 1.7 trajectory_error")


def update_feedforward(feedforward, gain, error):
    """Return the next trial's feedforward sequence using the previous trial's error.

    error has shape (101,2), including the initial row; gain has shape (100,202)."""
    # TODO(student, Problem 1.7): Implement update_feedforward using your derivation.
    raise NotImplementedError("Implement Problem 1.7 update_feedforward")


# %% Provided experiment and output

def main():
    parser = helpers.parser(__doc__, OUTPUT_DIR, INPUT_DIR)
    parser.add_argument("--trial-dir", type=Path, default=TRIAL_DIR,
                        help="Directory containing your Question 1.1 disturbance.")
    parser.add_argument("--trials", type=int, default=5)
    helpers.add_seed_argument(parser)
    args = parser.parse_args()
    gain = helpers.load_results(args.input_dir, "learning_gain.npz")["L"]
    if args.seed is None:
        disturbance = helpers.load_results(args.trial_dir, "trial.npz")["disturbance"]
    else:
        disturbance = helpers.disturbance_sequence(args.seed)
    result = helpers.learning_trials(disturbed_step, trajectory_error, update_feedforward,
                                     gain, disturbance, args.trials)
    out = helpers.output_directory(args.output_dir)
    norms = helpers.save_learning(result, out)
    helpers.save_json(out / "summary.json", {"problem": "1.7", "seed": args.seed,
        "trials": args.trials, "error_l2_by_trial": norms.tolist(),
        "max_abs_input": float(np.max(abs(result["controls"]))),
        "max_active_nominal_constraints": int(result["solver_diagnostics"][:, :, 1].max())})
    print(f"Error norms: {norms}")
    print(f"Saved Question 1.7 results to {out}")


if __name__ == "__main__":
    main()
