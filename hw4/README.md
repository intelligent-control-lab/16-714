# Homework 4: iterative learning control and LQG

This code accompanies the ILC and LQG handout. Use the handout for the
mathematical questions and written discussions. These scripts implement
the numerical experiments with the original parameters and timing conventions.

Keep these files together:

```text
README.md
requirements.txt
hw4_helpers.py
hw4_inputs.py
problem_1.1.py
problem_1.2.py
problem_1.5.py
problem_1.6.py
problem_1.7.py
problem_2.2.py
problem_2.4.py
problem_2.5.py
```

## Setup

Use Python 3.10 or newer with NumPy, SciPy, Matplotlib, and OSQP. HW4
runs headlessly without SPARK, MuJoCo, or a GPU. From the repository root,
install the standalone dependencies with:

```bash
python -m pip install -r hw4/requirements.txt
```

## Complete the exercises

Implement the functions marked `TODO(student, Problem ...)`. The helpers
provide the previous homework's nominal MPC solver, experiment drivers,
plotting, and file output. The exercise-specific equations belong in your
functions.

| Script | Functions to complete | Earlier results used |
|---|---|---|
| `problem_1.1.py` | `disturbed_step` | None |
| `problem_1.2.py` | `trajectory_error` | 1.1 |
| `problem_1.5.py` | `lqr_feedback` | None |
| `problem_1.6.py` | `lifted_error_matrix`, `learning_gain` | 1.5 |
| `problem_1.7.py` | `disturbed_step`, `trajectory_error`, `update_feedforward` | 1.1, 1.6 |
| `problem_2.2.py` | `cross_cost_block`, `augmented_gain` | 1.5 |
| `problem_2.4.py` | `augmented_model`, `steady_covariance`, `kalman_gain` | None |
| `problem_2.5.py` | `initial_distribution`, `feedback_control`, `kalman_step` | 2.2, 2.4 |

Copy **your own completed** plant update and error calculation from Questions
1.1 and 1.2 into Question 1.7. Later scripts load your saved matrices where
needed; they do not supply completed implementations of earlier exercises.
Question 2.2 reuses your state cost matrix and feedback gain from Question 1.5.
Questions 1.3, 1.4, 2.1, and 2.3 require written derivations used in the later
code. Submit the remaining explanations requested in the handout as well.

## Run and save results

From the repository root, complete the TODOs and run in this order:

```bash
python hw4/problem_1.1.py
python hw4/problem_1.2.py
python hw4/problem_1.5.py
python hw4/problem_1.6.py
python hw4/problem_1.7.py
python hw4/problem_2.2.py
python hw4/problem_2.4.py
python hw4/problem_2.5.py
```

An unfinished function raises `NotImplementedError`. Importing a script does
not run the experiment or write results. Scripts also work when invoked by
absolute path from another directory.

Every script accepts `--output-dir PATH`. The default output directory for
Question 1.1 is `results/student_problem_1_1/`, with corresponding names for
the other questions. Scripts with earlier dependencies use these directories
by default. Use `--input-dir PATH` to select earlier results; Question 1.7
also accepts `--trial-dir PATH` for the Question 1.1 run. Question 2.5 accepts
`--control-dir PATH` and `--filter-dir PATH` for Questions 2.2 and 2.4.

| Question | Generated files |
|---|---|
| 1.1 | Executed and predicted trajectories, controls, disturbance, nominal trajectory, solver diagnostics, and PNG/PDF plots |
| 1.2 | Error samples, stacked error, and PNG/PDF plots |
| 1.5 | State cost matrix, feedback gain, and closed-loop matrix |
| 1.6 | Lifted model, learning gain, and two gain heatmaps |
| 1.7 | All trial trajectories, errors, feedforward sequences, convergence table, and PNG/PDF plots |
| 2.2 | Cost blocks and augmented feedback gain |
| 2.4 | Augmented model, prior covariance, and filter gain |
| 2.5 | True and estimated states, controls, measurements, noise samples, and PNG/PDF plots |

Each question saves a JSON summary and NPZ numerical data. CSV files expose
the main tables. Arrays use one row per time sample, except that repeated
ILC trials add a leading trial dimension.

## Inputs and conventions

`hw4_inputs.py` contains the supplied disturbance, the disturbance-free
nominal trajectory from the previous MPC homework, and fixed random draws
for a reproducible LQG experiment. These are given experiment inputs.

The ILC experiment executes 100 steps, starting from `[2, 0]`, with a fixed
100-step MPC prediction horizon. The reference is
`r_k = max(1 - k/100, 0) * [2, 0]`. The state bounds are `[5, 5]`, and the
nominal input bound is `1`. The plant receives the nominal control plus the
learned feedforward and disturbance; the driver does not clip this sum.
In the Question 1.1 plot, predictions extend beyond the 100 executed steps
and the reference continues to time 200.

Errors contain all 101 state samples, including the initial zero row.
Question 1.2 also saves the written stack for times 1 through 100 as
`lifted_error`. The ILC calculation retains the initial row, giving a
`(202,100)` lifted matrix and a `(100,202)` learning matrix. The gain plots
omit the two columns corresponding to the initial error and display two
`(100,100)` panels using the same fixed color limits, `[-10,10]`.

By default, Question 1.7 runs five trials with the disturbance saved by your
Question 1.1 script. Each trial resets the state and reuses the same
disturbance. For the three disturbance sequences requested in the handout,
use the default run plus two additional runs, for example:

```bash
python hw4/problem_1.7.py --seed 1 --output-dir hw4/results/student_problem_1_7_seed1
python hw4/problem_1.7.py --seed 2 --output-dir hw4/results/student_problem_1_7_seed2
```

Questions 1.1, 1.7, and 2.5 accept `--seed INTEGER` to generate a new experiment.
Omitting it replays the provided inputs. The same seed reproduces the same
inputs within this Python workflow; numerical solver tolerances can still
cause small differences between installations.

For LQG, the initial estimate is the prior mean. The driver first uses a
measurement at time 1, then predicts and updates through time 99. Saved true
states and measurements have 101 samples (times 0–100); estimates and controls
have 100 samples (times 0–99). The fixed replay inputs reproduce a reference
run of the original numerical procedure. The historical LQG figure's random
inputs were not saved, so its individual noise trace may differ.

The dots in filenames are intentional; run by file path. Use ordinary
`python` on macOS as well. Generated `hw4/results/` and local `hw4/release/`
packages are ignored by Git. No completed HW4 implementations, instructor
solutions, grading rubrics, or generated HW4 results are distributed.
The nominal HW3 trajectory in `hw4_inputs.py` is a supplied HW4 input.

## Check the installation

Run any script with `--help` to check its imports and available options
without executing an exercise or writing results, for example:

```bash
python hw4/problem_1.1.py --help
python hw4/problem_2.5.py --help
```

If an earlier result is missing, the script identifies the missing file.
Complete and run the prerequisite question, or select its saved results
using the input-directory options described above. No files from another
course checkout are required.
