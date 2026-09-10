# SPARK compatibility checks

Run these tests from the `16-714` repository root whenever the SPARK
checkout changes:

```bash
conda activate spark_course
python -B -m unittest discover -s unit_tests -v
python -B -m scripts.generate_results --all
python -B -m scripts.validate_results
python scripts/check_repository.py
```

The supported SPARK baseline is the public
[intelligent-control-lab/spark](https://github.com/intelligent-control-lab/spark)
`main` branch. Install it with the `mujoco` profile and `--dev`, then install
`sympy` for HW2; the course does not require learned-policy, Isaac, ROS, or
hardware SDK extras. CI runs the full suite without a private repository token.

The suite checks:

- public SPARK imports required by the lectures;
- robot-config-owned dynamics and course-owned numerical model stepping;
- LQR, MPC, and iLQR construction from the expected robot configurations;
- estimator, MRAC, and ILC migration boundaries;
- lecture package/configuration structure;
- deterministic result summaries against compact baselines;
- student-scaffold imports, supplied helper interfaces, and homework release
  boundaries without any dependency on completed homework scripts.

HW2-only checks need just `hw2/requirements.txt`, not SPARK:

```bash
python -B -m unittest unit_tests.test_hw2_release -v
```

These are release checks, not grading tests. They expect the distributed
TODO placeholders to remain unfinished. Homework 1's output smoke tests
use fixed dummy inputs, not the assignment controller. Neither suite embeds
homework answers or loads an instructor repository.

A future SPARK commit is compatible if this suite and the result validation
both pass. Do not add course-specific behavior to SPARK to satisfy this suite.
