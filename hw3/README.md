# Homework 3: constrained MPC and generalized LQR

This code accompanies the HW3 handout on constrained MPC and LQR with a
state-control cross term. Use the handout for the mathematical questions,
parameters, plots, and written discussions. This README describes the Python
workflow, required implementations, and generated outputs.

Keep these files together:

```text
README.md
requirements.txt
hw3_helpers.py
problem_1.2.py
problem_1.3.py
problem_1.4.py
problem_2.4.py
```

## Setup

Use Python 3.10 or newer with NumPy, SciPy, Matplotlib, and OSQP. HW3
runs headlessly and does not require SPARK, MuJoCo, or a GPU. From this
repository root, install the standalone dependencies with:

```bash
python -m pip install -r hw3/requirements.txt
```

## Complete the exercises

Implement the functions marked `TODO(student, Problem ...)` in each problem
script. The helper supplies solver adapters, numerical iteration, generic
polygon operations, simulation, plotting, and file output. It takes the
exercise-specific mathematics from your functions.

| Script | Functions to complete |
|---|---|
| `problem_1.2.py` | `lifted_dynamics`, `feasibility_constraints` |
| `problem_1.3.py` | `predecessor_intersection` |
| `problem_1.4.py` | `lifted_dynamics`, `cost_matrices`, `condensed_qp_matrices`, `online_qp_terms` |
| `problem_2.4.py` | `riccati_step`, `feedback_gain` |

Use your derivation from Question 1.1 for the lifted dynamics, costs, and
constraints in Questions 1.2 and 1.4. You may copy **your own completed**
`lifted_dynamics` function from Question 1.2 into Question 1.4. Both scripts
start with an empty function and run independently.

Question 1.4 lists its given parameters in the script: `X0 = [2, 0]`,
`DT = 0.1`, `Q = I`, `R = 0.1`, `S = 10I`, `HORIZON = 100`, state bounds
`[5, 5]`, and input bound `1`. The provided reference is
`r_k = max(1 - k/N, 0) * X0`. The driver executes 100 MPC steps with a fixed
100-step prediction horizon, using the reference window `r_k, ..., r_(k+N)`
and applying only the first control at each step. The four TODOs implement
the QP from Question 1.1 with the stated box constraints. The provided driver
does not add a terminal invariant-set constraint.

Question 1.3 uses the provided generic scalar-projection and polygon routines.
Question 2.4 implements your control law and Riccati equation from Questions
2.2–2.3. It uses different system and cost parameters from the MPC questions;
use the constants in its own script. Its variable `CROSS_TERM` denotes the
matrix called N in the LQR question.

Questions 1.1, 2.1–2.3, and 2.5 require written work. Submit the plots and
discussions requested in the handout as well. The Python functions
implement the handout's tasks; they are not additional written questions.

## Run and save results

From this repository root, run after completing the relevant TODOs:

```bash
python hw3/problem_1.2.py
python hw3/problem_1.3.py
python hw3/problem_1.4.py
python hw3/problem_2.4.py
```

An unfinished function raises `NotImplementedError` with its name. Importing
a script does not run the experiment or write files. Scripts also work when
invoked by absolute path from another directory.

All scripts accept `--output-dir PATH`. By default they save results in
`results/student_problem_1_2/`, `student_problem_1_3/`, `student_problem_1_4/`,
and `student_problem_2_4/`, respectively, beside the scripts.

Question 1.2 samples states with spacing 0.1 by default. For a quick
development run, use `python hw3/problem_1.2.py --resolution 1`; use the default
spacing for your final plots.

| Question | Generated files |
|---|---|
| 1.2 | Feasibility grids in CSV/NPZ, a summary, and feasible-set PNG/PDF plots |
| 1.3 | Polygon vertices and halfspaces, iteration history, a summary, and invariant-set PNG/PDF plots |
| 1.4 | Executed and predicted trajectories, controls, reference, solver diagnostics, a summary, and PNG/PDF plots |
| 2.4 | Gain and Riccati matrices, iteration history, closed-loop poles, and a summary |

State arrays use `[position, velocity]`, with one row per time step. Stack
states in time order, keeping both components adjacent; function docstrings
specify the required shapes. In Question 1.4, distinguish the reference,
predicted trajectories, and executed trajectory when interpreting the plots.
For Question 1.4, `mpc_trajectories.png` and `.pdf` put all three on the same
figure, with separate position and velocity panels: blue predictions, red
execution, and a black dashed reference. In `mpc_trajectories.npz`,
`predicted_states[k, j]` is the state predicted at time `k+j` by the solve at
time `k`. There are 100 predictions of 101 states each, 101 executed states,
100 executed controls, and reference samples for `k = 0, ..., 200`.

The dots in the filenames are intentional; invoke these scripts by file path.
Use ordinary `python` on macOS as well. Generated `results/` and locally
assembled `release/` directories are ignored by Git.

## Check the installation

From the repository root, check that each script imports and exposes its CLI:

```bash
python hw3/problem_1.2.py --help
python hw3/problem_1.3.py --help
python hw3/problem_1.4.py --help
python hw3/problem_2.4.py --help
```

These commands do not execute the exercises or write results. Running an
unfinished exercise raises `NotImplementedError` identifying the mathematical
function to complete. Completed implementations, written solutions, grading
rubrics, and reference outputs are excluded from this release.
