"""Student-release smoke checks; no completed controllers or reference outputs."""

import ast
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np


HW2 = Path(__file__).resolve().parents[1] / "hw2"
QUESTIONS = tuple(f"1.{number}" for number in range(4, 9))


def load_problem(question, directory=HW2):
    path = Path(directory) / f"problem_{question}.py"
    spec = importlib.util.spec_from_file_location(
        "released_problem_" + question.replace(".", "_"), path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class HW2ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modules = {question: load_problem(question) for question in QUESTIONS}

    def test_imports_and_common_array_contract(self):
        for question, module in self.modules.items():
            with self.subTest(question=question):
                self.assertTrue(callable(module.main))
                self.assertEqual(module.NUM_STATE_SAMPLES, 101)
                self.assertEqual(module.X0.shape, (4,))
                self.assertEqual(module.X_GOAL.shape, (4,))
                self.assertEqual(module.Q.shape, (4, 4))
                self.assertEqual(module.R.shape, (2, 2))
                self.assertEqual(module.S.shape, (4, 4))
                self.assertEqual(module.OUTPUT_DIR.parent, HW2 / "results")
                self.assertEqual(
                    module.OUTPUT_DIR.name, "student_problem_" + question.replace(".", "_")
                )

    def test_all_released_placeholders_are_preserved(self):
        expected = {
            "1.4": {"backward_step": 1, "forward_step": 1},
            "1.5": {"linearize_about_reference": 1, "riccati_step": 1,
                    "lqr_forward_step": 1},
            "1.6": {"nonlinear_lqr_step": 1},
            "1.7": {"delta_stage_step": 1, "delta_control_law": 1},
            "1.8": {"solve_ilqr": 3},
        }
        for question in QUESTIONS:
            tree = ast.parse((HW2 / f"problem_{question}.py").read_text())
            actual = {}
            for function in (node for node in tree.body if isinstance(node, ast.FunctionDef)):
                count = sum(
                    isinstance(node, ast.Raise) and isinstance(node.exc, ast.Call)
                    and isinstance(node.exc.func, ast.Name)
                    and node.exc.func.id == "NotImplementedError"
                    for node in ast.walk(function)
                )
                if count:
                    actual[function.name] = count
            with self.subTest(question=question):
                self.assertEqual(actual, expected[question])

    def test_each_entry_point_reports_an_unfinished_math_block(self):
        for question, module in self.modules.items():
            with self.subTest(question=question), tempfile.TemporaryDirectory() as directory:
                original_output = module.OUTPUT_DIR
                try:
                    module.OUTPUT_DIR = Path(directory)
                    with self.assertRaisesRegex(NotImplementedError, "Problem " + question):
                        if question in ("1.7", "1.8"):
                            module.main(["--output-dir", directory])
                        else:
                            module.main()
                finally:
                    module.OUTPUT_DIR = original_output

    def test_supplied_symbolic_conversion_produces_numeric_matrices(self):
        module = self.modules["1.4"]
        symbolic = module.symbolic_model()
        self.assertEqual(symbolic["dynamics"].shape, (4, 1))
        self.assertEqual(symbolic["A"].shape, (4, 4))
        self.assertEqual(symbolic["B"].shape, (4, 2))
        for array in module.evaluated_jacobians(symbolic):
            self.assertTrue(np.issubdtype(array.dtype, np.floating))
            self.assertTrue(np.all(np.isfinite(array)))

    def test_question_1_8_reuses_sibling_helpers_and_parameters(self):
        module = self.modules["1.8"]
        owner = module._problem_1_7
        for name in (
            "nonlinear_step", "state_jacobian", "problem_1_6_reference",
            "trajectory_cost", "delta_stage_step", "delta_model_terms",
            "delta_backward_pass", "delta_control_law", "rollout_update",
            "nonlinear_costate", "stationarity_residual", "save_csv", "save_trajectory",
        ):
            with self.subTest(function=name):
                self.assertIs(getattr(module, name), getattr(owner, name))
        for name in ("DT", "N", "NUM_STATE_SAMPLES", "STEP_SIZE", "X0", "X_GOAL",
                     "Q", "R", "S", "B", "A_GOAL", "STATE_NAMES", "CONTROL_NAMES"):
            self.assertIs(getattr(module, name), getattr(owner, name))
        with self.assertRaisesRegex(NotImplementedError, "Problem 1.7"):
            module.delta_control_law(np.zeros(4), np.zeros(4), np.zeros((2, 4)),
                                     np.zeros(2), module.STEP_SIZE)

    def test_supplied_reference_and_local_model_shapes(self):
        module = self.modules["1.7"]
        states, controls, costates = module.problem_1_6_reference()
        self.assertEqual(states.shape, (module.N + 1, 4))
        self.assertEqual(controls.shape, (module.N, 2))
        self.assertEqual(costates.shape, states.shape)
        np.testing.assert_array_equal(states[0], module.X0)
        self.assertTrue(np.isfinite(module.trajectory_cost(states, controls)))
        model = module.delta_model_terms(states, controls)
        for array, shape in zip(model, ((module.N, 4, 4), (module.N, 4),
                                       (module.N + 1, 4), (module.N, 2))):
            self.assertEqual(array.shape, shape)
            self.assertTrue(np.all(np.isfinite(array)))

    def test_supplied_jacobian_matches_finite_differences(self):
        module = self.modules["1.7"]
        state = np.array([0.2, -0.3, 0.7, 0.4])
        control = np.array([0.1, -0.2])
        delta = 1e-6 * np.eye(4)
        numeric = np.column_stack([
            (module.nonlinear_step(state + direction, control)
             - module.nonlinear_step(state - direction, control)) / (2e-6)
            for direction in delta
        ])
        np.testing.assert_allclose(module.state_jacobian(state), numeric, atol=1e-8)

    def test_sibling_import_is_portable_and_does_not_run_main(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script_dir = root / "homework files"
            script_dir.mkdir()
            other_dir = root / "other working directory"
            other_dir.mkdir()
            for question in ("1.7", "1.8"):
                shutil.copyfile(HW2 / f"problem_{question}.py",
                                script_dir / f"problem_{question}.py")
            for question in ("1.7", "1.8"):
                run = subprocess.run(
                    [sys.executable, "-B", str(script_dir / f"problem_{question}.py"), "--help"],
                    cwd=other_dir, capture_output=True, text=True, timeout=30,
                )
                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertIn("--output-dir", run.stdout)
            self.assertFalse((script_dir / "results").exists())

    def test_missing_sibling_has_a_clear_error(self):
        with tempfile.TemporaryDirectory() as directory:
            shutil.copyfile(HW2 / "problem_1.8.py", Path(directory) / "problem_1.8.py")
            with self.assertRaisesRegex(FileNotFoundError, "completed problem_1.7.py"):
                load_problem("1.8", directory)

    def test_no_instructor_files_or_external_course_imports(self):
        self.assertFalse(list(HW2.glob("solution*.py")))
        self.assertFalse((HW2 / "soln").exists())
        self.assertFalse((HW2 / "grad").exists())
        allowed_modules = {"__future__", "argparse", "importlib", "json", "pathlib",
                           "matplotlib", "numpy", "sympy"}
        for question in QUESTIONS:
            source = (HW2 / f"problem_{question}.py").read_text()
            self.assertNotIn("solution_", source)
            self.assertNotIn("/home/", source)
            for node in ast.walk(ast.parse(source)):
                if isinstance(node, ast.Import):
                    imports = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom):
                    imports = [node.module]
                else:
                    continue
                for name in imports:
                    self.assertIn(name.split(".")[0], allowed_modules)


if __name__ == "__main__":
    unittest.main()
