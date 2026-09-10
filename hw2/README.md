# Homework 2: discrete-time optimal control

This folder contains the student scaffolds for Questions 1.4–1.8. Use the
assignment handout for the mathematical definitions, required plots, and
convergence criterion. Complete only the marked `TODO(student, ...)` blocks;
the remaining code supplies model helpers, numerical safeguards, plotting,
and result saving.

## Setup and commands

HW2 runs headlessly with Python 3.10 or newer. It does not require SPARK,
MuJoCo, a GPU, or another course repository. From this repository's root,
install the numerical dependencies in your active environment:

```bash
python -m pip install -r hw2/requirements.txt
```

After completing the relevant TODOs, run:

```bash
python hw2/problem_1.4.py
python hw2/problem_1.5.py
python hw2/problem_1.6.py
python hw2/problem_1.7.py
python hw2/problem_1.8.py
```

The dots in the filenames are intentional; run them by file path, not with
`python -m hw2.problem_1.8`. Unfinished mathematical blocks raise
`NotImplementedError`, identifying the question and function to complete.
This is expected, not a missing dependency. HW2 does not open a viewer, so
macOS also uses ordinary `python` for these commands.

## Implementation scope

| Question | Marked implementations |
|---|---|
| 1.4 | `backward_step`, `forward_step` |
| 1.5 | `linearize_about_reference`, `riccati_step`, `lqr_forward_step` |
| 1.6 | `nonlinear_lqr_step` |
| 1.7 | `delta_stage_step`, `delta_control_law` |
| 1.8 | Three blocks in `solve_ilqr`: iteration calls, cost decrease and nominal update, termination |

Questions 1.4–1.7 can be run independently. Question 1.8 imports the shared
functions and parameters from the sibling `problem_1.7.py`, including your
completed mathematical functions. Keep these two files together. You do not
need to run Question 1.7 first or generate its output files. Both scripts
initialize from the supplied Question 1.6 reference; iteration 0 is that
reference, and the first Q1.8 update corresponds to Q1.7.

Use the prescribed fixed step size in Questions 1.7 and 1.8. No line search,
step-size tuning, or symbolic-library implementation is required. The
iteration safeguards and diagnostics are provided; follow the handout for
the termination condition. See [the symbolic-computation guide](README_symbolic_solver.md)
for the supplied SymPy and NumPy workflow.

## Outputs and checks

Each script saves CSV trajectories, compressed NumPy arrays, a JSON summary,
and PNG figures under its own directory:
`hw2/results/student_problem_1_4/` through
`hw2/results/student_problem_1_8/`. Output paths are relative to the script
location, not the current working directory. Questions 1.7 and 1.8 also
accept `--output-dir PATH` to select a different destination.

The `states` and `costates` arrays have 101 rows; `controls` has 100 rows.
The state order is `[p_x, p_y, v, theta]`, and the control order is
`[acceleration, steering_rate]`. Use the generated figures for your report
and the numerical files to check array shapes and trajectory consistency.
Generated outputs are ignored by Git.

Run the release smoke checks without completing the TODOs:

```bash
python -B -m unittest unit_tests.test_hw2_release -v
```

These checks validate imports, the supplied helpers, and the scaffold
boundaries. They are release checks, not an answer checker or a grading
suite; the placeholder checks intentionally expect the original TODOs.
No completed homework scripts, written solutions, grading rubrics, or
reference-output files are included in this release.
