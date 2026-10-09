"""Provided HW3 MPC solver, HW4 experiment drivers, plotting, and file output.

Exercise-specific equations are supplied by callbacks in the problem scripts.
State arrays have one time sample per row. The reference inputs are provided
in hw4_inputs.py so all experiments run without external data files.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import osqp
from scipy import sparse

import hw4_inputs as inputs

DT = 0.1
A = np.array([[1.0, DT], [0.0, 1.0]])
B = np.array([[0.5 * DT**2], [DT]])
C = np.array([[1.0, 0.0]])
Q = np.eye(2)
R = np.array([[0.1]])
S = 10.0 * np.eye(2)
X0 = np.array([2.0, 0.0])
HORIZON = 100
STATE_BOUND = np.array([5.0, 5.0])
INPUT_BOUND = 1.0
PROCESS_VARIANCE = 0.0001
MEASUREMENT_VARIANCE = 0.001
INITIAL_STATE_COVARIANCE = 0.001 * np.eye(2)
INITIAL_DISTURBANCE_VARIANCE = 0.0001


def parser(description, output_dir, input_dir=None):
    result = argparse.ArgumentParser(description=description)
    result.add_argument("--output-dir", type=Path, default=output_dir)
    if input_dir is not None:
        result.add_argument("--input-dir", type=Path, default=input_dir,
                            help="Directory containing your earlier results.")
    return result


def add_seed_argument(result):
    result.add_argument("--seed", type=int, default=None,
                        help="Generate a new experiment; omit to replay the reference inputs.")


def output_directory(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def save_csv(path, columns, values):
    np.savetxt(path, values, delimiter=",", header=",".join(columns), comments="")


def save_figure(fig, directory, name):
    for suffix in ("png", "pdf"):
        fig.savefig(Path(directory) / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def load_results(directory, filename):
    path = Path(directory) / filename
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}. Run the earlier question or set its input directory.")
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key].copy() for key in data.files}


def disturbance_sequence(seed=None):
    if seed is None:
        return inputs.REFERENCE_DISTURBANCE.copy()
    return np.random.default_rng(seed).normal(0.0, 0.3, HORIZON)


class NominalMPC:
    """Provided HW3 nominal MPC with a fixed 100-step prediction horizon."""

    def __init__(self):
        n = HORIZON
        powers = [np.linalg.matrix_power(A, k) for k in range(n + 1)]
        self.f = np.vstack(powers)
        self.g = np.zeros((2 * (n + 1), n))
        for k in range(1, n + 1):
            for j in range(k):
                self.g[2*k:2*k+2, j:j+1] = powers[k-j-1] @ B
        qbar = np.kron(np.eye(n + 1), Q)
        qbar[-2:, -2:] = S
        self.gradient_map = self.g.T @ qbar
        hessian = self.gradient_map @ self.g + np.kron(np.eye(n), R)
        self.matrix = np.vstack([self.g, -self.g, np.eye(n), -np.eye(n)])
        self.reference = np.maximum(1 - np.arange(2*n+1) / n, 0)[:, None] * X0
        self.solver = osqp.OSQP()
        self.solver.setup(P=sparse.csc_matrix((hessian + hessian.T) / 2),
                          q=np.zeros(n), A=sparse.csc_matrix(self.matrix),
                          l=np.full(len(self.matrix), -np.inf), u=np.zeros(len(self.matrix)),
                          verbose=False, eps_abs=1e-10, eps_rel=1e-10, max_iter=50000)

    def plan(self, state, k):
        offset = self.f @ state
        reference = self.reference[k:k+HORIZON+1].reshape(-1)
        upper = np.concatenate([np.tile(STATE_BOUND, HORIZON+1) - offset,
                                np.tile(STATE_BOUND, HORIZON+1) + offset,
                                np.full(2*HORIZON, INPUT_BOUND)])
        self.solver.update(q=self.gradient_map @ (offset-reference), u=upper)
        result = self.solver.solve()
        if result.info.status_val != 1 or not np.all(np.isfinite(result.x)):
            raise RuntimeError(f"MPC failed at step {k}: {result.info.status}")
        residual = self.matrix @ result.x - upper
        violation = max(0.0, float(residual.max()))
        if violation > 1e-7:
            raise RuntimeError(f"MPC constraint residual {violation} exceeds tolerance.")
        predicted = (offset + self.g @ result.x).reshape(HORIZON+1, 2)
        active = int(np.count_nonzero(abs(residual) < 1e-6))
        return result.x.copy(), predicted, (violation, active)


def disturbed_rollout(step_function, disturbance, feedforward=None):
    disturbance = np.asarray(disturbance, dtype=float)
    feedforward = np.zeros(HORIZON) if feedforward is None else np.asarray(feedforward, dtype=float)
    if disturbance.shape != (HORIZON,) or feedforward.shape != (HORIZON,):
        raise ValueError("Disturbance and feedforward must each have 100 samples.")
    controller = NominalMPC()
    states = np.zeros((HORIZON+1, 2))
    states[0] = X0
    controls, feedback = np.zeros(HORIZON), np.zeros(HORIZON)
    predictions = np.zeros((HORIZON, HORIZON+1, 2))
    plans, diagnostics = np.zeros((HORIZON, HORIZON)), np.zeros((HORIZON, 2))
    for k in range(HORIZON):
        plans[k], predictions[k], diagnostics[k] = controller.plan(states[k], k)
        feedback[k] = plans[k, 0]
        controls[k] = feedback[k] + feedforward[k]
        states[k+1] = step_function(states[k], controls[k], disturbance[k], A, B)
        if not np.all(np.isfinite(states[k+1])):
            raise ValueError(f"The plant update returned a nonfinite state at step {k+1}.")
    return dict(states=states, controls=controls, feedback=feedback,
                disturbance=disturbance.copy(), feedforward=feedforward.copy(),
                predicted_states=predictions, predicted_controls=plans,
                solver_diagnostics=diagnostics, reference=controller.reference,
                nominal_states=inputs.NOMINAL_STATES.copy())


def save_trial(result, out):
    np.savez_compressed(out / "trial.npz", **result)
    save_csv(out / "states.csv", ("k", "position", "velocity"),
             np.column_stack([np.arange(HORIZON+1), result["states"]]))
    save_csv(out / "controls.csv", ("k", "u", "feedback", "feedforward"),
             np.column_stack([np.arange(HORIZON), result["controls"],
                              result["feedback"], result["feedforward"]]))
    save_csv(out / "disturbance.csv", ("k", "w"),
             np.column_stack([np.arange(HORIZON), result["disturbance"]]))
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for k in range(HORIZON):
        for component in range(2):
            axes[component].plot(k+np.arange(HORIZON+1), result["predicted_states"][k, :, component],
                                 "b", linewidth=.5, label="Predicted" if k == 0 else None)
        axes[2].plot(k+np.arange(HORIZON), result["predicted_controls"][k],
                     "b", linewidth=.5, label="Predicted" if k == 0 else None)
    for component in range(2):
        axes[component].plot(result["states"][:, component], "r", label="Executed")
        axes[component].plot(result["reference"][:, component], "k--", label="Reference")
    axes[2].plot(result["controls"], "r", label="Executed")
    for axis, label in zip(axes, ("Position", "Velocity", "Control")):
        axis.set_ylabel(label)
        axis.legend(fontsize=8)
    axes[-1].set_xlabel("Time step k")
    fig.tight_layout()
    save_figure(fig, out, "disturbed_trajectory")


def plot_error(error, out):
    fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    for component, axis in enumerate(axes):
        axis.plot(error[:, component], "r")
        axis.set_ylabel(("Position error", "Velocity error")[component])
    axes[-1].set_xlabel("Time step k")
    fig.tight_layout()
    save_figure(fig, out, "tracking_error")


def plot_learning_gain(gain, out):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, constrained_layout=True)
    for component, axis in enumerate(axes):
        # The two e_0 columns are retained in the calculation and omitted from the plots.
        artist = axis.imshow(gain[:, 2+component::2], origin="upper", aspect="auto",
                             interpolation="nearest", cmap="viridis", vmin=-10, vmax=10,
                             extent=(.5, HORIZON+.5, HORIZON+.5, .5))
        axis.set(xlabel="Error time k (1 to 100)",
                 title=("Position-error gain", "Velocity-error gain")[component])
        fig.colorbar(artist, ax=axis, label="Gain", shrink=.85)
    axes[0].set_ylabel("Feedforward index (1 to 100)")
    save_figure(fig, out, "learning_gain")


def learning_trials(step_function, error_function, update_function, gain, disturbance, trials):
    if trials < 1:
        raise ValueError("The number of trials must be positive.")
    feedforward = np.zeros(HORIZON)
    runs = []
    for iteration in range(trials):
        result = disturbed_rollout(step_function, disturbance, feedforward)
        error = np.asarray(error_function(result["nominal_states"], result["states"]))
        if error.shape != (HORIZON+1, 2) or not np.all(np.isfinite(error)):
            raise ValueError("Return finite errors for times 0 through 100 with shape (101,2).")
        runs.append(dict(states=result["states"], controls=result["controls"], errors=error,
                         feedforward=feedforward.copy(), feedback=result["feedback"],
                         solver_diagnostics=result["solver_diagnostics"]))
        if iteration + 1 < trials:
            feedforward = np.asarray(update_function(feedforward, gain, error))
            if feedforward.shape != (HORIZON,) or not np.all(np.isfinite(feedforward)):
                raise ValueError("Return 100 finite feedforward values.")
    return dict(disturbance=np.asarray(disturbance), nominal_states=inputs.NOMINAL_STATES.copy(),
                reference=result["reference"],
                **{key: np.asarray([run[key] for run in runs]) for key in runs[0]})


def save_learning(result, out):
    np.savez_compressed(out / "learning_trials.npz", **result)
    count = len(result["states"])
    norms = np.linalg.norm(result["errors"].reshape(count, -1), axis=1)
    save_csv(out / "convergence.csv", ("trial", "error_l2", "max_abs_error", "max_abs_input"),
             np.column_stack([np.arange(1, count+1), norms,
                              np.max(abs(result["errors"]), axis=(1, 2)),
                              np.max(abs(result["controls"]), axis=1)]))
    fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    for iteration in range(count):
        color = str(1-(iteration+1)/count)
        for component, axis in enumerate(axes):
            axis.plot(result["errors"][iteration, :, component], color=color, label=f"Trial {iteration+1}")
    for component, axis in enumerate(axes):
        axis.set_ylabel(("Position error", "Velocity error")[component])
        axis.legend(fontsize=8)
    axes[-1].set_xlabel("Time step k")
    fig.tight_layout()
    save_figure(fig, out, "learning_errors")
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for iteration in range(count):
        color = str(1-(iteration+1)/count)
        for component in range(2):
            axes[component].plot(result["states"][iteration, :, component], color=color,
                                 label=f"Trial {iteration+1}")
        axes[2].plot(result["controls"][iteration], color=color, label=f"Trial {iteration+1}")
    for component in range(2):
        axes[component].plot(result["reference"][:HORIZON+1, component], "k--", label="Reference")
    for axis, label in zip(axes, ("Position", "Velocity", "Control")):
        axis.set_ylabel(label)
        axis.legend(fontsize=8)
    axes[-1].set_xlabel("Time step k")
    fig.tight_layout()
    save_figure(fig, out, "learning_trajectories")
    return norms


def simulate_lqg(distribution_function, control_function, filter_function, control_gain, model, seed):
    a, b, c, d, gain = (model[key] for key in ("A", "B", "C", "D", "F"))
    mean, covariance = distribution_function(X0, INITIAL_STATE_COVARIANCE,
                                              INITIAL_DISTURBANCE_VARIANCE)
    if seed is None:
        initial_draw = inputs.INITIAL_DRAW
        process_noise = np.sqrt(PROCESS_VARIANCE) * inputs.PROCESS_DRAWS
        measurement_noise = np.sqrt(MEASUREMENT_VARIANCE) * inputs.MEASUREMENT_DRAWS
    else:
        rng = np.random.default_rng(seed)
        initial_draw = rng.standard_normal(3)
        process_noise, measurement_noise = np.zeros(HORIZON), np.zeros(HORIZON+1)
        measurement_noise[0] = rng.normal(0.0, np.sqrt(MEASUREMENT_VARIANCE))
        for k in range(HORIZON):
            process_noise[k] = rng.normal(0.0, np.sqrt(PROCESS_VARIANCE))
            measurement_noise[k+1] = rng.normal(0.0, np.sqrt(MEASUREMENT_VARIANCE))
    states, estimates = np.zeros((HORIZON+1, 3)), np.zeros((HORIZON, 3))
    controls, measurements = np.zeros(HORIZON), np.zeros(HORIZON+1)
    states[0] = mean + np.linalg.cholesky(covariance) @ initial_draw
    estimates[0] = mean
    measurements[0] = (c @ states[0]).item() + measurement_noise[0]
    for k in range(HORIZON):
        # The first control uses the prior mean; filtering begins with y_1.
        if k > 0:
            estimates[k] = filter_function(estimates[k-1], controls[k-1], measurements[k],
                                           a, b, c, gain)
        controls[k] = control_function(estimates[k], control_gain)
        states[k+1] = a @ states[k] + b[:, 0]*controls[k] + d[:, 0]*process_noise[k]
        measurements[k+1] = (c @ states[k+1]).item() + measurement_noise[k+1]
    if not all(np.all(np.isfinite(x)) for x in (states, estimates, controls, measurements)):
        raise ValueError("The simulation produced nonfinite values.")
    return dict(states=states, estimates=estimates, controls=controls, measurements=measurements,
                process_noise=process_noise, measurement_noise=measurement_noise,
                initial_mean=mean, initial_covariance=covariance)


def save_lqg(result, out):
    np.savez_compressed(out / "lqg_trajectory.npz", **result)
    save_csv(out / "states.csv", ("k", "position", "velocity", "disturbance", "measurement"),
             np.column_stack([np.arange(HORIZON+1), result["states"], result["measurements"]]))
    save_csv(out / "estimates.csv", ("k", "position", "velocity", "disturbance"),
             np.column_stack([np.arange(HORIZON), result["estimates"]]))
    save_csv(out / "controls.csv", ("k", "u"),
             np.column_stack([np.arange(HORIZON), result["controls"]]))
    fig, axes = plt.subplots(5, 1, figsize=(9, 10), sharex=True)
    for component in range(3):
        axes[component].plot(result["states"][:, component], "r", label="True")
        axes[component].plot(result["estimates"][:, component], "k--", label="Estimated")
        axes[component].legend(fontsize=8)
    axes[2].set_ylim(-.5, .5)
    axes[3].plot(result["controls"], "k")
    axes[4].plot(result["measurements"], "k")
    for axis, label in zip(axes, ("Position", "Velocity", "Disturbance", "Control", "Measurement")):
        axis.set_ylabel(label)
    axes[-1].set_xlabel("Time step k")
    fig.tight_layout()
    save_figure(fig, out, "lqg_trajectory")
