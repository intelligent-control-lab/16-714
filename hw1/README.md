# Homework 1: unicycle control with bicycle-model mismatch

Homework 1 compares the future predicted by a four-state unicycle model with
the trajectory produced by a six-state dynamic bicycle model. The released
student entry point is `problem.py`. Complete its marked TODOs before running.

The unicycle state and control are

```text
x_u = [p_x, p_y, v, theta]
u_u = [v_dot, theta_dot]
```

The bicycle follows the Lecture 2 and SPARK ordering

```text
x_b = [X, Y, v_x, v_y, r, psi]
u_b = [a_x, delta]
```

Students complete the marked implementations and analysis:

1. Implement the two-gain unicycle controller and plot its trajectory and
   state histories.
2. Implement the provided pure-rolling state/control conversion, then repeat
   the plots with the tire-force dynamic bicycle plant.
3. Run the provided nine-pair gain study and analyze the position and
   orientation model-error plots.

Run the student scaffold after completing its TODOs:

```bash
python -m hw1.problem --model unicycle
python -m hw1.problem --model bicycle
python -m hw1.problem --gain-study
```

Add `--viewer --real-time --show-simulation-info` to either model command to
replay its saved numerical trajectory in the SPARK MuJoCo viewer. Results are
written under `hw1/results/`. The locally assembled distribution package is
written under `hw1/release/`; both directories are reproducible and ignored by
Git.

On macOS, viewer commands use `mjpython` instead of `python`, for example
`mjpython -m hw1.problem --model unicycle --viewer --real-time`.
The launcher is included with MuJoCo. Headless commands keep using `python`.
See the [macOS viewer note](../README.md#run-examples).

The model transitions are provided in `vehicle_models.py`. Plotting, CSV
output, numerical rollout, and viewer replay helpers are isolated in
`unicycle_bicycle_helpers.py`, keeping the assignment entry script below 200
lines.
