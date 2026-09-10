# 16-714 SPARK Control Examples

This repository contains Python course-control examples built on top of SPARK.
It keeps lecture-level problem setup, plotting, and artifact generation here,
while robot execution and SPARK policy implementations stay in SPARK.

The examples fall into two groups:

- SPARK robot examples: lecture 3 runs on the AgiBot G1 right arm, while
  lectures 7, 9, and 10 run on its mobile-base configurations.
- Local numerical examples: lectures 2, 12, 13, 14, 15, 16, 17, 18, 19,
  21, 22, and 23 use SPARK dynamics, estimation, adaptive control, iterative
  learning control, and model stepping. Course-only rollout state, plotting,
  artifact I/O, and reinforcement-learning experiments live in this repo.

Each lecture is a self-contained Python package: its experiment configuration,
implementation, entry point, and generated results all live in the same lecture
folder. Generated `results/` directories are ignored by git.

## Repository Layout

```text
16-714/
  README.md
  requirements.txt  # course-wide Python dependencies
  shared/
    artifacts.py
    control_pipeline.py
    numerical_rollout.py
    rl/
  lecture2/
    config.py
    main.py
    __main__.py
    results/          # generated and git-ignored
  ...
  lecture7/
    lqr/
      config.py
      main.py
      results/
    ilqr/
      config.py
      main.py
      results/
  lecture10/
    config.py
    main.py
    results/
  lecture21/
    config.py
    main.py
    value_learning.py # when lecture-specific support code is needed
    results/
  ...
  lecture23/
    config.py
    main.py
    results/
  hw*/               # assignment-specific student material
    README.md
    problem*.py
    ...              # supplied helpers and assignment dependencies
    results/         # generated and git-ignored
  scripts/
    generate_results.py
    validate_results.py
  unit_tests/
    README.md
    test_spark_compatibility.py
    test_simulation_info.py
    test_lecture2_models.py
    test_course_architecture.py
    test_lecture_dynamics_configs.py
    test_spark_library_migration.py
```

Lecture-specific code stays inside its lecture package. `shared/`
contains only infrastructure genuinely shared by several lectures: artifact
writing, the trace-collecting SPARK pipeline, course numerical rollouts, and
decomposed RL building blocks. Reusable controllers, estimators, and dynamics
models remain SPARK dependencies rather than being copied into course folders.

The root compatibility scripts and the former `lib/` implementations have been
removed. `scripts/` and `unit_tests/` are tracked repository code because they
define the reproducible generation and compatibility-validation workflow.

## Homework

Each `hwN/` folder contains its student scripts, supplied helpers, and a
README with assignment-specific setup, commands, and implementation tasks.
Complete the marked TODOs before running; an unfinished block raises
`NotImplementedError`. Homework answers and reference outputs are not
distributed here.

## Dependencies

Use a SPARK environment installed from the public upstream
[intelligent-control-lab/spark](https://github.com/intelligent-control-lab/spark). The scripts
expect these imports to work:

```python
import spark_pipeline
import spark_policy
import spark_robot
```

The SPARK checkout must include:

- `AgiBotG1RightArmTeleopPipelineConfig`
- `AgiBotG1TeleopPipelineConfig`
- `AgiBotG1RightArmAgent`
- `AgiBotG1MobileBaseAgent`
- `AgiBotG1RightArmDynamic1Config`
- `AgiBotG1RightArmKinematics`
- `AgiBotG1MobileBaseDynamic1Config`
- `AgiBotG1MobileBaseDynamic2Config`
- `AgiBotG1MobileBaseUnicycleDynamic1Config`
- `AgiBotG1MobileBaseBicycleDynamic2Config`
- `ILQRPolicy`
- `LQRPolicy`
- `FiniteHorizonLQRPolicy`
- `LinearMPCPolicy`
- `InputConstrainedMPCPolicy`
- `ConstrainedMPCPolicy`
- `KalmanFilterEstimator`
- `SteadyStateKalmanFilterEstimator`
- `ExtendedKalmanFilterEstimator`
- `UnscentedKalmanFilterEstimator`
- `RecursiveLeastSquaresEstimator`
- `GradientParameterEstimator`
- `ModelReferenceAdaptiveController`
- `FrequencyDomainILC`
- `TimeDomainILC`

The course-wide Python dependencies are listed in `requirements.txt`,
covering numerical computation, plotting, optimization, and symbolic algebra.
CUDA PyTorch, Isaac, ROS, and hardware SDK packages are not required for the
course examples.

## Install SPARK

Clone public SPARK main and install the MuJoCo development profile first:

```bash
git clone --branch main --single-branch https://github.com/intelligent-control-lab/spark.git
cd spark
./install.sh --name spark_course --profile mujoco --dev
conda activate spark_course
```

The course targets public SPARK `main`. Confirm the checkout with:

```bash
git branch --show-current
git rev-parse HEAD
```

The compatibility suite below is the authority for a newer SPARK revision;
course code must not depend on a private course branch or copied SPARK source.

If you already have a local SPARK checkout, run the same installer from that
checkout instead.

After activating the environment, confirm that SPARK imports correctly:

```bash
python - <<'PY'
import spark_pipeline
import spark_policy
import spark_robot
print("SPARK imports ok")
PY
```

If the imports fail, install the SPARK packages from the SPARK checkout before
running this repository.

From this course repository's root, install the course-wide dependencies:

```bash
cd /path/to/16-714
python -m pip install -r requirements.txt
```

GitHub Actions installs the same public SPARK MuJoCo profile and course
dependencies, then runs the compatibility tests, result validation, and
repository hygiene checks.
No additional repository-access token or private checkout is required, and
checkout credentials are not persisted.

## Repository Model

This repository is an experiment and teaching repository, not a reusable Python
library. Run its module entry points from the repository root and install SPARK
itself into the active environment. A `pip install -e .` workflow is not needed
for the current layout: installing top-level packages named `lecture2`,
`lecture3`, and so on would add course-specific names to the environment without
providing a stable public API.

If reusable course infrastructure emerges later, move that API into a normal
source package such as `src/spark_course/` and add a `pyproject.toml` then. The
lecture folders can remain experiment entry points that depend on that package.

## Run Examples

Run from this repository root:

```bash
cd /path/to/16-714
conda activate spark_course
python -m lecture12
```

Each lecture writes outputs to its own `results/` directory and records its
standalone configuration in `config.json`. For example, Lecture 12 writes to
`lecture12/results/`, while Lecture 7 iLQR writes to
`lecture7/ilqr/results/`. Robot lectures additionally save planned and
executed state/control traces plus a compact `summary.json`.

```bash
python -m lecture7.lqr --show-simulation-info
python -m lecture7.ilqr --show-simulation-info
python -m lecture9 --viewer --show-simulation-info
python -m lecture10 --viewer --show-simulation-info
```

On macOS, commands that open the MuJoCo viewer may need `mjpython` in place of
`python`; the passive viewer requires rendering on the main thread.
`mjpython` is included with MuJoCo and accepts the same arguments. For example,
after activating the environment:

```bash
mjpython -m lecture9 --viewer --show-simulation-info
```

Headless numerical commands can still use `python`. See the
[MuJoCo viewer documentation](https://mujoco.readthedocs.io/en/stable/python.html#passive-viewer).

Run the SPARK/course compatibility suite after installing SPARK or updating the
SPARK checkout:

```bash
python -B -m unittest discover -s unit_tests -v
python -B -m scripts.generate_results --all
python -B -m scripts.validate_results
python scripts/check_repository.py
```

The repository check performs a cache-free syntax pass and rejects generated
metadata, obsolete compatibility scripts, missing homework release sources,
instructor-only homework files, and references to the retired course branch.

## SPARK Robot Examples

These scripts use the Agibot G1 mobile base through SPARK:

```bash
python -m lecture7.lqr --show-simulation-info
python -m lecture7.ilqr --show-simulation-info
python -m lecture9 --viewer --show-simulation-info
python -m lecture10 --viewer --show-simulation-info
```

Course viewer configurations explicitly opt in to SPARK 2.0's simulation
information panel. The `--show-simulation-info` flag is therefore enabled by
default whenever a viewer is active; use `--no-show-simulation-info` to hide
the panel for a particular run.

They use SPARK pipeline configuration and the shared `ControlPipeline` wrapper
in this repository. The wrapper follows the standard SPARK runtime contract:

```python
action, action_info = policy.act(agent_feedback, task_info)
```

Planner-style policies may expose a planned trajectory through:

```python
action_info["policy_plan"] = {
    "states": states,
    "controls": controls,
    "state_t": state_t,
    "control_t": control_t,
    "state_names": state_names,
    "control_names": control_names,
    "solver_info": solver_info,
}
```

`ControlPipeline` saves planned and executed traces when that information is
available.

## Implemented Examples

| Lecture | Topic | Entry Point |
|---|---|---|
| 2 | Vehicle model simulation | `python -m lecture2` |
| 3 | AgiBot G1 joint and Cartesian arm servo | `python -m lecture3` |
| 7 LQR | Infinite-horizon LQR | `python -m lecture7.lqr` |
| 7 iLQR | iLQR on unicycle base model | `python -m lecture7.ilqr` |
| 9 | LQR, finite-horizon LQR, linear MPC | `python -m lecture9` |
| 10 | Linear, input-constrained, and constrained MPC | `python -m lecture10` |
| 12 | Frequency-domain ILC | `python -m lecture12` |
| 13 | Time-domain ILC | `python -m lecture13` |
| 14 | Least-squares conditional estimate | `python -m lecture14` |
| 15 | Linear and steady-state Kalman filters | `python -m lecture15` |
| 16 | EKF and UKF | `python -m lecture16` |
| 17 | RLS and SGD parameter estimation | `python -m lecture17` |
| 18 | Separation principle | `python -m lecture18` |
| 19 | Model reference adaptive control | `python -m lecture19` |
| 21 | MC, SARSA, and Q-learning | `python -m lecture21` |
| 22 | Policy gradient and ILC comparison | `python -m lecture22` |
| 23 | REINFORCE, baseline, and actor critic | `python -m lecture23` |

Homework 1 is released as a three-part scaffold. Its models run independently:

```bash
python -m hw1.problem --model unicycle
python -m hw1.problem --model bicycle
python -m hw1.problem --gain-study
```

Each command saves state/control CSV data, one XY trajectory figure, and one
single-curve figure per state under `hw1/results/<model>/`. Students implement
the controller and unicycle/bicycle conversion, then analyze the provided
gain/error study. See
[`hw1/README.md`](hw1/README.md) for the model conventions, commands, and
deliverables. Add
`--viewer --real-time --show-simulation-info` to either command for visual
replay. The bicycle uses the six-state `(X, Y, v_x, v_y, r, psi)` convention
shared by Lecture 2 and SPARK.

## Notes

- SPARK robot configurations own all dynamics equations. Policies obtain a
  reduced dynamics view from the selected robot configuration instead of
  selecting a model by name.
- SPARK owns reusable estimators, MRAC, ILC, and dynamics models. Course
  numerical state, stopping, and trajectory collection live in
  `shared/numerical_rollout.py`; reinforcement-learning representations and
  update methods remain under `shared/rl`.
- Numerical experiments call the selected SPARK model's public `step` method.
  The course wrapper owns experiment state and time, while the model owns the
  Euler, RK4, ZOH, or native-discrete transition.
- Robot examples depend on SPARK for policy construction, robot configuration,
  agent execution, and simulation.
- Every package-local `results/` directory is generated output and is ignored.
- Homework implementations, constants, plots, and result files stay in their
  corresponding `hwN/` folder. Reusable toolbox improvements are not applied
  to SPARK merely to encode an assignment.

Regenerate deterministic numerical results and validate the compact committed
baselines with:

```bash
python -m scripts.generate_results
python -m scripts.validate_results
```

Use `python -m scripts.generate_results --all` to include the simulator-backed
robot lectures in headless mode.
Use `python -m scripts.generate_results --all --viewer --show-simulation-info`
to replay all viewer-capable lectures with their information panels visible.

Lecture 3 is headless by default for deterministic artifact generation. Use
`python -m lecture3 --viewer --real-time --show-simulation-info` to replay its
AgiBot right-arm experiments in SPARK's MuJoCo viewer.

Lecture 2 also runs headless by default. Use
`python -m lecture2 --viewer --real-time --show-simulation-info` to replay its
integrator, unicycle, and bicycle trajectories with the AgiBot mobile base in
the SPARK viewer.
The viewer overlay follows the active trajectory source, including its robot
configuration, dynamics variant, reduced state/control dimensions, integrator,
sample period, and playback time; the persistent MuJoCo visualization agent is
not presented as the lecture's active dynamics model.
