#!/usr/bin/env python3
"""Student scaffold for 16-714 HW2, Question 1.5.

Students implement the general linearization, one Riccati backward step, and
one LQR forward step. Plotting and result serialization are provided.
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
OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "student_problem_1_5"
STATE_NAMES = ("p_x", "p_y", "v", "theta")
CONTROL_NAMES = ("acceleration", "steering_rate")


# %% Student mathematical implementation: general linearization
def linearize_about_reference(reference_state: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (A^r, B, c^r) such that f(x,u) ~= c^r + A^r x + B u."""

    # TODO(student, Problem 1.5): Derive and implement the affine
    # linearization about the general state [p_x^r, p_y^r, v^r, theta^r].
    raise NotImplementedError("Implement Problem 1.5 linearize_about_reference")


# %% Student mathematical implementation: backward propagation
def riccati_step(
    p_next: np.ndarray, a: np.ndarray, b: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return (P_k, K_k) for one finite-horizon LQR backward step."""

    # TODO(student, Problem 1.5): Implement the Riccati recursion and LQR
    # gain. Solve the control Hessian system without forming its inverse.
    raise NotImplementedError("Implement Problem 1.5 riccati_step")


def backward_propagation(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p_matrices = np.zeros((NUM_STATE_SAMPLES, 4, 4))
    gains = np.zeros((HORIZON, 2, 4))
    p_matrices[-1] = S
    for k in range(HORIZON - 1, -1, -1):
        p_matrices[k], gains[k] = riccati_step(p_matrices[k + 1], a, b)
    return p_matrices, gains


# %% Student mathematical implementation: forward simulation
def lqr_forward_step(
    state: np.ndarray,
    p_current: np.ndarray,
    gain: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    affine_offset: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (u_k, x_{k+1}, lambda_k) for LQR on the linearized dynamics."""

    # TODO(student, Problem 1.5): Evaluate the optimal LQR law, propagate the
    # affine linear model, and evaluate the co-state.
    raise NotImplementedError("Implement Problem 1.5 lqr_forward_step")


def forward_simulation(
    a: np.ndarray,
    b: np.ndarray,
    affine_offset: np.ndarray,
    p_matrices: np.ndarray,
    gains: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    states = np.zeros((NUM_STATE_SAMPLES, 4))
    controls = np.zeros((HORIZON, 2))
    costates = np.zeros((NUM_STATE_SAMPLES, 4))
    states[0] = X0
    for k in range(HORIZON):
        controls[k], states[k + 1], costates[k] = lqr_forward_step(
            states[k], p_matrices[k], gains[k], a, b, affine_offset
        )
    costates[-1] = S @ (states[-1] - X_GOAL)
    return states, controls, costates


# %% Provided verification, saving, and plotting
def controllability_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.hstack([np.linalg.matrix_power(a, power) @ b for power in range(4)])


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
    output = OUTPUT_DIR / "problem_1_5_trajectories.png"
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
    a, b, affine_offset = linearize_about_reference(X_GOAL)
    p_matrices, gains = backward_propagation(a, b)
    states, controls, costates = forward_simulation(a, b, affine_offset, p_matrices, gains)
    total_cost = trajectory_cost(states, controls)
    rank = int(np.linalg.matrix_rank(controllability_matrix(a, b)))

    save_csv(OUTPUT_DIR / "states.csv", ("k",) + STATE_NAMES, np.column_stack((np.arange(NUM_STATE_SAMPLES), states)))
    save_csv(OUTPUT_DIR / "controls.csv", ("k",) + CONTROL_NAMES, np.column_stack((np.arange(HORIZON), controls)))
    save_csv(OUTPUT_DIR / "costates.csv", ("k",) + tuple(f"lambda_{name}" for name in STATE_NAMES), np.column_stack((np.arange(NUM_STATE_SAMPLES), costates)))
    np.savez_compressed(OUTPUT_DIR / "lqr_backward_recursion.npz", P=p_matrices, K=gains, A=a, B=b, c=affine_offset)
    figure = save_plot(states, controls, costates)
    summary = {"problem": "1.5", "controllability_rank": rank, "total_cost": total_cost, "final_state": states[-1].tolist(), "figure": figure.name}
    (OUTPUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Problem 1.5 complete: rank(C) = {rank}, J = {total_cost:.9f}; outputs: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
