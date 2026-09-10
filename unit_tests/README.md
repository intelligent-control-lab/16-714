# SPARK compatibility checks

Run these tests from the `16-714` repository root whenever the SPARK
checkout changes:

```bash
conda activate spark_course
python -m pip install -r requirements.txt
python -B -m unittest discover -s unit_tests -v
python -B -m scripts.generate_results --all
python -B -m scripts.validate_results
python scripts/check_repository.py
```

The supported SPARK baseline is the public
[intelligent-control-lab/spark](https://github.com/intelligent-control-lab/spark)
`main` branch. Install it with the `mujoco` profile and `--dev`, then install
the course dependencies from the root `requirements.txt`. The course does
not require learned-policy, Isaac, ROS, or hardware SDK extras. CI runs the
full suite without a private repository token.

The suite checks:

- public SPARK imports required by the lectures;
- robot-config-owned dynamics and course-owned numerical model stepping;
- LQR, MPC, and iLQR construction from the expected robot configurations;
- estimator, MRAC, and ILC migration boundaries;
- lecture package/configuration structure;
- deterministic result summaries against compact baselines;
- student-scaffold imports, supplied helper interfaces, and homework release
  boundaries without any dependency on completed homework scripts.

For assignment-specific setup and smoke-test commands, see the README in
the corresponding homework folder.

These are release checks, not grading tests. They expect the distributed
TODO placeholders to remain unfinished. Output smoke tests use fixed dummy
inputs, not assignment controllers. The released suites do not embed
homework answers or load an instructor repository.

A future SPARK commit is compatible if this suite and the result validation
both pass. Do not add course-specific behavior to SPARK to satisfy this suite.
