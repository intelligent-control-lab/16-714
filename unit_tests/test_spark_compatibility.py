import importlib
from pathlib import Path
import unittest

import numpy as np

from spark_robot import AgiBotG1MobileBaseUnicycleDynamic1Config


REQUIRED_PUBLIC_EXPORTS = {
    "spark_agent": (
        "DynamicsExecutor",
        "DynamicsModelAgent",
        "AgiBotG1MobileBaseAgent",
        "AgiBotG1RightArmAgent",
    ),
    "spark_env": (
        "SparkEnvironment",
    ),
    "spark_pipeline": (
        "AgiBotG1TeleopPipelineConfig",
        "AgiBotG1RightArmTeleopPipelineConfig",
    ),
    "spark_policy": (
        "ILQRPolicy",
        "LQRPolicy",
        "FiniteHorizonLQRPolicy",
        "LinearMPCPolicy",
        "InputConstrainedMPCPolicy",
        "ConstrainedMPCPolicy",
        "KalmanFilterEstimator",
        "SteadyStateKalmanFilterEstimator",
        "ExtendedKalmanFilterEstimator",
        "UnscentedKalmanFilterEstimator",
        "RecursiveLeastSquaresEstimator",
        "GradientParameterEstimator",
        "ModelReferenceAdaptiveController",
        "FrequencyDomainILC",
        "TimeDomainILC",
    ),
    "spark_robot": (
        "RobotDynamicsModel",
        "AgiBotG1RightArmDynamic1Config",
        "AgiBotG1RightArmKinematics",
        "AgiBotG1MobileBaseDynamic1Config",
        "AgiBotG1MobileBaseDynamic2Config",
        "AgiBotG1MobileBaseUnicycleDynamic1Config",
        "AgiBotG1MobileBaseBicycleDynamic2Config",
        "LinearDiscreteDynamicsConfig",
    ),
}


class SparkCompatibilityTests(unittest.TestCase):
    def test_required_public_exports_are_available(self):
        for module_name, names in REQUIRED_PUBLIC_EXPORTS.items():
            module = importlib.import_module(module_name)
            for name in names:
                with self.subTest(module=module_name, name=name):
                    self.assertTrue(hasattr(module, name))

    def test_robot_config_exposes_executable_dynamics_model(self):
        robot_cfg = AgiBotG1MobileBaseUnicycleDynamic1Config()
        model = robot_cfg.create_dynamics_model()
        self.assertIs(model.robot_cfg, robot_cfg)
        self.assertEqual(model.variant, "unicycle")
        state = model.step(np.zeros(3), np.array([1.0, 0.2]), 0.1, "RK4")
        np.testing.assert_allclose(
            state,
            [0.0999933334666654, 0.0009999666671111, 0.02],
            rtol=1e-8,
            atol=1e-10,
        )

    def test_ci_uses_public_spark_and_unconditional_checks(self):
        workflow = (
            Path(__file__).resolve().parents[1] / ".github/workflows/ci.yml"
        ).read_text(encoding="utf-8")
        self.assertEqual(workflow.count("actions/checkout@v6"), 2)
        self.assertIn("actions/setup-python@v6", workflow)
        self.assertIn("repository: intelligent-control-lab/spark", workflow)
        self.assertIn("ref: main", workflow)
        self.assertIn("--python 3.10 --profile mujoco --dev", workflow)
        self.assertIn("Install course dependencies", workflow)
        self.assertIn("pip install -r requirements.txt", workflow)
        self.assertNotRegex(workflow.lower(), r"hw\d+|homework\s+\d+")
        self.assertIn("unittest discover -s unit_tests -v", workflow)
        self.assertIn("scripts.generate_results --all", workflow)
        self.assertIn("scripts.validate_results", workflow)
        self.assertNotIn("secrets.", workflow)
        self.assertNotIn("github.repository_owner", workflow)
        self.assertNotIn("if:", workflow)
        self.assertNotIn("actions/checkout@v4", workflow)
        self.assertNotIn("actions/setup-python@v5", workflow)

    def test_course_setup_uses_shared_dependency_manifest(self):
        root = Path(__file__).resolve().parents[1]
        readme = (root / "README.md").read_text(encoding="utf-8")
        self.assertIn("pip install -r requirements.txt", readme)
        self.assertNotIn("hw2", readme.lower())
        self.assertNotIn("homework 2", readme.lower())

        def requirements(path):
            return {
                line.strip() for line in path.read_text(encoding="utf-8").splitlines()
                if line.strip() and not line.lstrip().startswith("#")
            }

        course_dependencies = requirements(root / "requirements.txt")
        self.assertTrue({"numpy", "scipy", "matplotlib", "osqp", "sympy"}
                        <= course_dependencies)
        for manifest in root.glob("hw*/requirements.txt"):
            with self.subTest(manifest=manifest.name):
                self.assertTrue(requirements(manifest) <= course_dependencies)


if __name__ == "__main__":
    unittest.main()
