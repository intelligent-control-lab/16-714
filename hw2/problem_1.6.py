#!/usr/bin/env python3
"""Student scaffold for 16-714 HW2, Question 1.6.

The Question 1.5 LQR backward pass is provided so students can focus on the
new mathematical step: applying those gains to the nonlinear plant and
evaluating the associated approximate co-state.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# %% Scenario and parameters
DT = 0.1
NUM_STATE_SAMPLES = 101
HORIZON = NUM_STATE_SAMPLES - 1
X0 = np.zeros(4)
X_GOAL = np.array([10.0, 10.0, 0.0, np.pi / 2.0])
Q = np.eye(4)
R = 0.1 * np.eye(2)
S = 10.0 * np.eye(4)
A_GOAL = np.array(
    [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, DT, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
)
B = np.array([[0.0, 0.0], [0.0, 0.0], [DT, 0.0], [0.0, DT]])
OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "student_problem_1_6"
STATE_NAMES = ("p_x", "p_y", "v", "theta")
CONTROL_NAMES = ("acceleration", "steering_rate")


# %% Provided dynamics and Question 1.5 backward pass
def nonlinear_step(state: np.ndarray, control: np.ndarray) -> np.ndarray:
    px, py, speed, theta = state
    acceleration, steering_rate = control
    return np.array(
        [
            px + DT * speed * np.cos(theta),
            py + DT * speed * np.sin(theta),
            speed + DT * acceleration,
            theta + DT * steering_rate,
        ]
    )


def finite_horizon_lqr() -> tuple[np.ndarray, np.ndarray]:
    """Provided result from Question 1.5: return all P_k and K_k."""

    p_matrices = np.zeros((NUM_STATE_SAMPLES, 4, 4))
    gains = np.zeros((HORIZON, 2, 4))
    p_matrices[-1] = S
    for k in range(HORIZON - 1, -1, -1):
        p_next = p_matrices[k + 1]
        control_hessian = R + B.T @ p_next @ B
        gains[k] = np.linalg.solve(control_hessian, B.T @ p_next @ A_GOAL)
        p_matrices[k] = (
            Q + A_GOAL.T @ p_next @ A_GOAL - A_GOAL.T @ p_next @ B @ gains[k]
        )
        p_matrices[k] = 0.5 * (p_matrices[k] + p_matrices[k].T)
    return p_matrices, gains


# %% Student mathematical implementation: LQR on nonlinear dynamics
def nonlinear_lqr_step(
    state: np.ndarray, p_current: np.ndarray, gain: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (u_k, x_{k+1}, lambda_k) using the fixed LQR gain."""

    # TODO(student, Problem 1.6): Evaluate the Question 1.5 LQR control law,
    # apply it to the original nonlinear dynamics, and evaluate lambda_k.
    raise NotImplementedError("Implement Problem 1.6 nonlinear_lqr_step")


def forward_simulation(
    p_matrices: np.ndarray, gains: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    states = np.zeros((NUM_STATE_SAMPLES, 4))
    controls = np.zeros((HORIZON, 2))
    costates = np.zeros((NUM_STATE_SAMPLES, 4))
    states[0] = X0
    for k in range(HORIZON):
        controls[k], states[k + 1], costates[k] = nonlinear_lqr_step(
            states[k], p_matrices[k], gains[k]
        )
    costates[-1] = S @ (states[-1] - X_GOAL)
    return states, controls, costates


# %% Provided verification, saving, and plotting
def trajectory_cost(states: np.ndarray, controls: np.ndarray) -> float:
    errors = states[:-1] - X_GOAL
    terminal_error = states[-1] - X_GOAL
    return float(
        0.5 * np.einsum("ni,ij,nj->", errors, Q, errors)
        + 0.5 * np.einsum("ni,ij,nj->", controls, R, controls)
        + 0.5 * terminal_error @ S @ terminal_error
    )


def save_csv(path: Path, header: tuple[str, ...], values: np.ndarray) -> None:
    np.savetxt(path, values, delimiter=",", header=",".join(header), comments="")


def save_plot(states: np.ndarray, controls: np.ndarray, costates: np.ndarray) -> Path:
    output = OUTPUT_DIR / "problem_1_6_trajectories.png"
    fig = plt.figure(figsize=(11.0, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.15))
    axes = (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :]))
    axes[0].plot(np.arange(HORIZON), controls[:, 0], label="acceleration")
    axes[0].plot(np.arange(HORIZON), controls[:, 1], label="steering rate")
    axes[0].set_title("Control")
    for index, name in enumerate(STATE_NAMES):
        axes[1].plot(np.arange(NUM_STATE_SAMPLES), costates[:, index], label=name)
        axes[2].plot(np.arange(NUM_STATE_SAMPLES), states[:, index], label=name)
    axes[1].set_title("Co-state")
    axes[2].set_title("State")
    for axis in axes:
        axis.set_xlabel("time step $k$")
        axis.grid(alpha=0.25)
        axis.legend()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    p_matrices, gains = finite_horizon_lqr()
    states, controls, costates = forward_simulation(p_matrices, gains)
    total_cost = trajectory_cost(states, controls)

    save_csv(OUTPUT_DIR / "states.csv", ("k",) + STATE_NAMES, np.column_stack((np.arange(NUM_STATE_SAMPLES), states)))
    save_csv(OUTPUT_DIR / "controls.csv", ("k",) + CONTROL_NAMES, np.column_stack((np.arange(HORIZON), controls)))
    save_csv(OUTPUT_DIR / "costates.csv", ("k",) + tuple(f"lambda_{name}" for name in STATE_NAMES), np.column_stack((np.arange(NUM_STATE_SAMPLES), costates)))
    np.savez_compressed(OUTPUT_DIR / "lqr_on_nonlinear.npz", P=p_matrices, K=gains, states=states, controls=controls, costates=costates)
    figure = save_plot(states, controls, costates)
    summary = {"problem": "1.6", "total_cost": total_cost, "final_state": states[-1].tolist(), "figure": figure.name}
    (OUTPUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Problem 1.6 complete: J = {total_cost:.9f}; outputs: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
