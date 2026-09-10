#!/usr/bin/env python3
"""Student scaffold for 16-714 HW2, Question 1.4.

The SymPy declarations, differentiation, numerical conversion, plotting, and
artifact I/O are provided. Students implement only the mathematical backward
recursion and the forward state/co-state/control update.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


# %% Scenario and parameters
DT = 0.1
NUM_STATE_SAMPLES = 101
HORIZON = NUM_STATE_SAMPLES - 1
X0 = np.zeros(4)
X_GOAL = np.array([10.0, 10.0, 0.0, np.pi / 2.0])
Q = np.eye(4)
R = 0.1 * np.eye(2)
S = 10.0 * np.eye(4)
OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "student_problem_1_4"
STATE_NAMES = ("p_x", "p_y", "v", "theta")
CONTROL_NAMES = ("acceleration", "steering_rate")


# %% Provided symbolic setup
def symbolic_model() -> dict[str, object]:
    """Construct f(x,u), df/dx, and df/du; no student work is required here."""

    px, py, speed, theta = sp.symbols("p_x p_y v theta", real=True)
    acceleration, steering_rate = sp.symbols("a omega", real=True)
    dt = sp.symbols("Delta_t", positive=True, real=True)
    state = sp.Matrix([px, py, speed, theta])
    control = sp.Matrix([acceleration, steering_rate])
    dynamics = sp.Matrix(
        [
            px + dt * speed * sp.cos(theta),
            py + dt * speed * sp.sin(theta),
            speed + dt * acceleration,
            theta + dt * steering_rate,
        ]
    )
    return {
        "state": state,
        "control": control,
        "dt": dt,
        "dynamics": dynamics,
        "A": sp.simplify(dynamics.jacobian(state)),
        "B": sp.simplify(dynamics.jacobian(control)),
    }


def evaluated_jacobians(model: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    substitutions = {model["dt"]: DT}
    substitutions.update(dict(zip(model["state"], X_GOAL)))
    substitutions.update({symbol: 0.0 for symbol in model["control"]})
    a_goal = np.asarray(model["A"].subs(substitutions), dtype=float)
    b = np.asarray(model["B"].subs({model["dt"]: DT}), dtype=float)
    a_goal[np.abs(a_goal) < 1.0e-14] = 0.0
    b[np.abs(b) < 1.0e-14] = 0.0
    return a_goal, b


def uncontrolled_step(state: np.ndarray) -> np.ndarray:
    px, py, speed, theta = state
    return np.array(
        [px + DT * speed * np.cos(theta), py + DT * speed * np.sin(theta), speed, theta]
    )


def nonlinear_step(state: np.ndarray, control: np.ndarray, b: np.ndarray) -> np.ndarray:
    return uncontrolled_step(state) + b @ control


# %% Student mathematical implementation: backward propagation
def backward_step(
    p_next: np.ndarray, a_goal: np.ndarray, b: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return (P_k, G_k) from P_{k+1}.

    G_k parameterizes the nonlinear policy as a function of the
    uncontrolled next-state expression. Derive it by solving the stationarity
    equation and state equation simultaneously. Then retain the first-order
    co-state terms about X_GOAL to obtain P_k.
    """

    # TODO(student, Problem 1.4): Implement the implicit optimal-control solve
    # and first-order co-state/Riccati update. Use np.linalg.solve rather than
    # forming an explicit numerical inverse.
    raise NotImplementedError("Implement Problem 1.4 backward_step")


def backward_propagation(a_goal: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p_matrices = np.zeros((NUM_STATE_SAMPLES, 4, 4))
    policy_coefficients = np.zeros((HORIZON, 2, 4))
    p_matrices[-1] = S
    for k in range(HORIZON - 1, -1, -1):
        p_matrices[k], policy_coefficients[k] = backward_step(
            p_matrices[k + 1], a_goal, b
        )
    return p_matrices, policy_coefficients


# %% Student mathematical implementation: forward simulation
def forward_step(
    state: np.ndarray,
    p_current: np.ndarray,
    policy_coefficient: np.ndarray,
    b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (u_k, x_{k+1}, lambda_k) for the nonlinear-policy forward pass."""

    # TODO(student, Problem 1.4): Evaluate the nonlinear symbolic policy,
    # propagate the nonlinear Euler dynamics, and evaluate lambda_k.
    raise NotImplementedError("Implement Problem 1.4 forward_step")


def forward_simulation(
    p_matrices: np.ndarray, policy_coefficients: np.ndarray, b: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    states = np.zeros((NUM_STATE_SAMPLES, 4))
    controls = np.zeros((HORIZON, 2))
    costates = np.zeros((NUM_STATE_SAMPLES, 4))
    states[0] = X0
    for k in range(HORIZON):
        controls[k], states[k + 1], costates[k] = forward_step(
            states[k], p_matrices[k], policy_coefficients[k], b
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
    output = OUTPUT_DIR / "problem_1_4_trajectories.png"
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
    model = symbolic_model()
    a_goal, b = evaluated_jacobians(model)
    p_matrices, policy_coefficients = backward_propagation(a_goal, b)
    states, controls, costates = forward_simulation(p_matrices, policy_coefficients, b)
    total_cost = trajectory_cost(states, controls)

    save_csv(OUTPUT_DIR / "states.csv", ("k",) + STATE_NAMES, np.column_stack((np.arange(NUM_STATE_SAMPLES), states)))
    save_csv(OUTPUT_DIR / "controls.csv", ("k",) + CONTROL_NAMES, np.column_stack((np.arange(HORIZON), controls)))
    save_csv(OUTPUT_DIR / "costates.csv", ("k",) + tuple(f"lambda_{name}" for name in STATE_NAMES), np.column_stack((np.arange(NUM_STATE_SAMPLES), costates)))
    np.savez_compressed(OUTPUT_DIR / "backward_recursion.npz", P=p_matrices, G=policy_coefficients, A_goal=a_goal, B=b)
    figure = save_plot(states, controls, costates)
    summary = {"problem": "1.4", "total_cost": total_cost, "final_state": states[-1].tolist(), "figure": figure.name}
    (OUTPUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Problem 1.4 complete: J = {total_cost:.9f}; outputs: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
