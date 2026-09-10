#!/usr/bin/env python3
"""Student scaffold for Question 1.8: fixed-step iLQR iterations.

Complete the iteration calls, nominal-trajectory update, and termination
test in solve_ilqr. Pass the current trajectory, backward-pass result, and
fixed step size to the appropriate helpers on every iteration.
First complete delta_stage_step and delta_control_law in problem_1.7.py,
and keep that file beside this one. The supplied loader reuses its common
parameters and helpers: dynamics, Jacobian, initialization, cost, backward
pass, rollout, diagnostics, and trajectory saving. It does not run the
Question 1.7 main function. Only the iteration logic and Question 1.8-specific
plotting/output setup live here; all non-TODO code is provided.
No symbolic-library coding, line search, or step-size tuning is required.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# %% Provided import: reuse your Question 1.7 functions without copying their answers.
def _load_problem_1_7():
    """Load the sibling file by path because its filename contains a dot."""
    path = Path(__file__).resolve().with_name("problem_1.7.py")
    if not path.is_file():
        raise FileNotFoundError(
            f"Question 1.8 requires {path}. Keep your completed problem_1.7.py "
            "in the same directory as problem_1.8.py."
        )
    spec = importlib.util.spec_from_file_location("hw2_problem_1_7", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Question 1.7 from {path}.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_problem_1_7 = _load_problem_1_7()

# Reuse the same function objects, including their Question 1.7 dependencies.
nonlinear_step = _problem_1_7.nonlinear_step
state_jacobian = _problem_1_7.state_jacobian
problem_1_6_reference = _problem_1_7.problem_1_6_reference
trajectory_cost = _problem_1_7.trajectory_cost
delta_stage_step = _problem_1_7.delta_stage_step
delta_model_terms = _problem_1_7.delta_model_terms
delta_backward_pass = _problem_1_7.delta_backward_pass
delta_control_law = _problem_1_7.delta_control_law
rollout_update = _problem_1_7.rollout_update
nonlinear_costate = _problem_1_7.nonlinear_costate
stationarity_residual = _problem_1_7.stationarity_residual
save_csv = _problem_1_7.save_csv
save_trajectory = _problem_1_7.save_trajectory


# %% Share the prescribed parameters with the imported Question 1.7 functions.
# Change common settings in problem_1.7.py, not independently in this file.
DT = _problem_1_7.DT
NUM_STATE_SAMPLES = _problem_1_7.NUM_STATE_SAMPLES
N = _problem_1_7.N
STEP_SIZE = _problem_1_7.STEP_SIZE  # Fixed feedforward relaxation, no line search.
MAX_ITERATIONS = 200  # Provided internal safety cap, not a convergence criterion.
RELATIVE_COST_TOLERANCE = 1.0e-6  # Sole convergence test; trajectory changes are diagnostics.
X0 = _problem_1_7.X0
X_GOAL = _problem_1_7.X_GOAL
Q = _problem_1_7.Q
R = _problem_1_7.R
S = _problem_1_7.S
B = _problem_1_7.B
A_GOAL = _problem_1_7.A_GOAL
STATE_NAMES = _problem_1_7.STATE_NAMES
CONTROL_NAMES = _problem_1_7.CONTROL_NAMES


# %% Question-specific command-line settings and output directory.
def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-size", type=float, default=STEP_SIZE,
                        help="Constant step throughout the run; 1 uses the full control update.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args(argv)
    if not 0.0 < args.step_size <= 1.0:
        parser.error("--step-size must lie in (0, 1].")
    return args


OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "student_problem_1_8"


# %% Complete the iteration and termination logic; bookkeeping is provided.
def solve_ilqr(step_size: float = STEP_SIZE) -> dict:
    """Complete the three TODO blocks in the fixed-step iteration loop.

    Each backward pass must use the latest nominal trajectory. Retain every
    finite update, even if its cost increases. Check the signed cost decrease
    after storing the update, and stop at the first qualifying iteration.
    State/control changes are diagnostics, not stopping conditions.
    Here states/controls hold the reference; rollout_update returns its
    fresh rollout arrays as new_states/new_controls before reference replacement.
    """
    states, controls, costates = problem_1_6_reference()
    state_history = [states.copy()]
    control_history = [controls.copy()]
    costate_history = [costates.copy()]
    backward_history = []
    metrics = [{
        "iteration": 0, "cost": trajectory_cost(states, controls),
        "relative_cost_decrease": 0.0, "state_change": 0.0, "control_change": 0.0,
    }]
    stop_reason = "maximum_iterations"

    for iteration in range(1, MAX_ITERATIONS + 1):
        # TODO(student, Problem 1.8, iteration): Call delta_backward_pass
        # with the current states and controls, saving its dict as backward.
        # Call rollout_update with that same nominal trajectory, backward,
        # and the passed step_size. Unpack new_states, new_controls, and
        # new_costates. Read old_cost from metrics[-1]["cost"] and compute
        # new_cost with trajectory_cost on the updated trajectory.
        raise NotImplementedError("Implement Problem 1.8 backward/forward iteration")

        # Provided finite-value guard and trajectory-change diagnostics.
        if not (np.isfinite(new_cost) and np.all(np.isfinite(new_states))
                and np.all(np.isfinite(new_controls)) and np.all(np.isfinite(new_costates))):
            stop_reason = "nonfinite_rollout"
            break
        state_change = float(np.max(np.abs(new_states - states)))
        control_change = float(np.max(np.abs(new_controls - controls)))

        # TODO(student, Problem 1.8, update): Compute relative_cost_decrease
        # using the signed, normalized decrease in the problem statement.
        # Make new_states and new_controls the next nominal states and
        # controls. Retain them even when new_cost is larger than old_cost.
        raise NotImplementedError("Implement Problem 1.8 cost decrease and nominal update")

        # Provided storage: record the completed update before testing convergence.
        state_history.append(states.copy())
        control_history.append(controls.copy())
        costate_history.append(new_costates.copy())
        backward_history.append(backward)
        metrics.append({
            "iteration": iteration, "cost": new_cost,
            "relative_cost_decrease": relative_cost_decrease,
            "state_change": state_change, "control_change": control_change,
        })
        # TODO(student, Problem 1.8, termination): Test the nonnegative
        # guard and strict upper bound using RELATIVE_COST_TOLERANCE.
        # If satisfied, set stop_reason = "relative_cost_tolerance" and
        # break. Otherwise continue, including after a cost increase.
        raise NotImplementedError("Implement Problem 1.8 termination condition")

    return {
        "states": np.asarray(state_history), "controls": np.asarray(control_history),
        "costates": np.asarray(costate_history), "backward_history": backward_history,
        "metrics": metrics, "stop_reason": stop_reason,
    }


def selected_iterations(number_of_iterates: int) -> list[int]:
    choices = [0, 1, 2, 5, 10, number_of_iterates - 1]
    return sorted(set(index for index in choices if 0 <= index < number_of_iterates))


def iteration_colors(count: int):
    color_map = plt.get_cmap("viridis")
    if count == 1:
        return [color_map(0.5)]
    return [color_map(index / (count - 1)) for index in range(count)]


def plot_state_iterations(state_history: list[np.ndarray], selected: list[int]) -> Path:
    output = OUTPUT_DIR / "problem_1_8_state_iterations.png"
    colors = iteration_colors(len(selected))
    steps = np.arange(N + 1)
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 7.5), constrained_layout=True)
    for iteration, color in zip(selected, colors):
        states = state_history[iteration]
        label = f"iteration {iteration}"
        axes[0, 0].plot(states[:, 0], states[:, 1], color=color, label=label)
        for state_index, axis in enumerate(axes.flat[1:5]):
            axis.plot(steps, states[:, state_index], color=color, label=label)
    axes[0, 0].plot(X_GOAL[0], X_GOAL[1], "r*", markersize=11, label="goal")
    axes[0, 0].set_xlabel("$p_x$")
    axes[0, 0].set_ylabel("$p_y$")
    axes[0, 0].set_title("planar path")
    for axis, name in zip(axes.flat[1:5], STATE_NAMES):
        axis.set_title(name)
        axis.set_xlabel("time step $k$")
    axes[1, 2].axis("off")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    axes[1, 2].legend(handles, labels, loc="center", fontsize=9)
    for axis in axes.flat[:5]:
        axis.grid(alpha=0.25)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_control_iterations(control_history: list[np.ndarray], selected: list[int]) -> Path:
    output = OUTPUT_DIR / "problem_1_8_control_iterations.png"
    colors = iteration_colors(len(selected))
    steps = np.arange(N)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharex=True, constrained_layout=True)
    for iteration, color in zip(selected, colors):
        for control_index, axis in enumerate(axes):
            axis.plot(steps, control_history[iteration][:, control_index], color=color, label=f"iteration {iteration}")
    for axis, name in zip(axes, CONTROL_NAMES):
        axis.set_title(name)
        axis.set_xlabel("time step $k$")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_costate_iterations(costate_history: list[np.ndarray], selected: list[int]) -> Path:
    output = OUTPUT_DIR / "problem_1_8_costate_iterations.png"
    colors = iteration_colors(len(selected))
    steps = np.arange(N + 1)
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True, constrained_layout=True)
    for iteration, color in zip(selected, colors):
        for state_index, axis in enumerate(axes.flat):
            axis.plot(steps, costate_history[iteration][:, state_index], color=color, label=f"iteration {iteration}")
    for axis, name in zip(axes.flat, STATE_NAMES):
        axis.set_title(f"$\\lambda_{{{name}}}$")
        axis.set_xlabel("time step $k$")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_convergence(metrics: list[dict[str, float]]) -> Path:
    output = OUTPUT_DIR / "problem_1_8_convergence.png"
    iterations = np.array([item["iteration"] for item in metrics])
    costs = np.array([item["cost"] for item in metrics])
    relative = np.array([item["relative_cost_decrease"] for item in metrics[1:]])
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    axes[0].plot(iterations, costs)
    axes[0].set_title("objective")
    axes[0].set_ylabel("$J$")
    axes[1].plot(iterations[1:], relative)
    # Keep signed values visible: a log axis cannot show zero or increases.
    if relative.size and np.all(relative > 0.0):
        axes[1].set_yscale("log")
    else:
        axes[1].axhline(0.0, color="gray", linestyle=":", label="zero decrease")
    axes[1].axhline(RELATIVE_COST_TOLERANCE, color="red", linestyle="--", label="tolerance")
    axes[1].set_title("relative objective decrease")
    axes[1].legend()
    for axis in axes:
        axis.set_xlabel("iteration")
        axis.grid(alpha=0.25)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def main(argv=None) -> None:
    global OUTPUT_DIR
    args = parse_args(argv)
    OUTPUT_DIR = args.output_dir.resolve()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    run = solve_ilqr(args.step_size)
    state_history, control_history, costate_history = run["states"], run["controls"], run["costates"]
    metrics, stop_reason = run["metrics"], run["stop_reason"]
    states, controls, costates = state_history[-1], control_history[-1], costate_history[-1]
    nonlinear_costates = nonlinear_costate(states)
    final_backward = delta_backward_pass(states, controls)
    selected = selected_iterations(len(state_history))
    save_trajectory(OUTPUT_DIR, states, controls, costates, "final_")
    save_csv(OUTPUT_DIR / "nonlinear_costates.csv",
             ("k",) + tuple("lambda_" + name for name in STATE_NAMES),
             np.column_stack((np.arange(N + 1), nonlinear_costates)))
    metric_names = ("iteration", "cost", "relative_cost_decrease", "state_change", "control_change")
    save_csv(OUTPUT_DIR / "iteration_metrics.csv", metric_names,
             np.array([[item[name] for name in metric_names] for item in metrics]))
    # Entry j stores the backward model used to produce trajectory j+1.
    backward_arrays = {
        name + "_history": np.asarray([item[name] for item in run["backward_history"]])
        for name in ("A", "d", "q", "r", "P", "b", "K", "k_ff")
    }
    np.savez_compressed(
        OUTPUT_DIR / "all_iterations.npz", states=state_history, controls=control_history,
        costates=costate_history, **backward_arrays, selected_iterations=np.asarray(selected),
        Q=Q, R=R, S=S, B=B, x0=X0, x_goal=X_GOAL, dt=DT, step_size=args.step_size,
    )
    np.savez_compressed(OUTPUT_DIR / "final_backward_pass.npz", **final_backward)
    figures = [
        plot_state_iterations(state_history, selected),
        plot_control_iterations(control_history, selected),
        plot_costate_iterations(costate_history, selected),
        plot_convergence(metrics),
    ]
    summary = {
        "problem": "1.8",
        "initialization": "Question 1.6 nonlinear LQR trajectory",
        "method": "affine co-state P,b recursion with a constant step",
        "horizon": N, "state_samples": NUM_STATE_SAMPLES, "dt": DT,
        "step_size": args.step_size,
        "costate_definition": "P_k (x_new_k - x_nominal_k) + b_k; iteration 0 uses Q1.6",
        "converged": stop_reason == "relative_cost_tolerance",
        "stop_reason": stop_reason, "iterations": len(state_history) - 1,
        "maximum_iterations": MAX_ITERATIONS,
        "relative_cost_tolerance": RELATIVE_COST_TOLERANCE,
        "convergence_criterion": "0 <= (J_old-J_new)/max(1,abs(J_old)) < cost_tol",
        "trajectory_changes_are_diagnostics": True,
        "initial_cost": metrics[0]["cost"], "final_cost": metrics[-1]["cost"],
        "final_relative_cost_decrease": metrics[-1]["relative_cost_decrease"],
        "final_state_change": metrics[-1]["state_change"],
        "final_control_change": metrics[-1]["control_change"],
        "final_state": states[-1].tolist(),
        "terminal_error_norm": float(np.linalg.norm(states[-1] - X_GOAL)),
        "nonlinear_stationarity_infinity_norm": float(np.max(np.abs(
            stationarity_residual(controls, nonlinear_costates)))),
        "local_vs_nonlinear_costate_infinity_norm": float(np.max(np.abs(costates - nonlinear_costates))),
        "selected_plot_iterations": selected,
        "figures": [path.name for path in figures],
    }
    (OUTPUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Question 1.8: constant step {args.step_size:g}, {summary['iterations']} iterations")
    print(f"J {summary['initial_cost']:.9f} -> {summary['final_cost']:.9f}; {stop_reason}")
    print(f"final state = {np.array2string(states[-1], precision=9)}")
    print(f"saved results to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
