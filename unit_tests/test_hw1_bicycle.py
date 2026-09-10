"""Release-safe HW1 helper checks; fixed test inputs are not homework answers."""

import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from spark_robot import (
    AgiBotG1MobileBaseBicycleDynamic2Config,
    AgiBotG1MobileBaseDynamic1Config,
)

from hw1.problem import (
    ComparisonPipeline,
    Config,
    run_gain_study,
)
from hw1 import problem
from hw1.unicycle_bicycle_helpers import ModelViewer


class _ViewerAgentSpy:
    def __init__(self):
        self.information = []
        self.reset_count = 0

    def set_viewer_simulation_info(self, information, *, replace=False):
        self.information.append((dict(information), replace))

    def is_running(self):
        return True

    def reset(self, _reset_info):
        self.reset_count += 1

    def begin_render_frame(self):
        pass

    def render_line_segment(self, *_args):
        pass

    def render_sphere(self, *_args):
        pass

    def render(self):
        pass


class Homework1BicycleTests(unittest.TestCase):
    def test_released_scaffold_exposes_expected_api(self):
        for name in ("unicycle_control", "bicycle_initial_state", "bicycle_control",
                     "run_gain_study"):
            self.assertTrue(callable(getattr(problem, name)))
        with self.assertRaisesRegex(NotImplementedError, "Problem 2.1"):
            problem.unicycle_control(np.zeros(4), problem.Config())

    def test_bicycle_uses_course_state_and_control_order(self):
        with tempfile.TemporaryDirectory() as directory:
            pipeline = ComparisonPipeline(Config(duration=0.1), directory)

        self.assertEqual(
            pipeline.bicycle_model.state_names,
            ("X", "Y", "v_x", "v_y", "r", "psi"),
        )
        self.assertEqual(pipeline.bicycle_model.control_names, ("a_x", "delta"))

    def test_local_euler_transition_matches_spark_euler_step(self):
        config = Config(duration=0.1)
        with tempfile.TemporaryDirectory() as directory:
            pipeline = ComparisonPipeline(config, directory)
        state = np.array([0.0, 0.0, 0.5, 0.0, 0.0, 0.0])
        control = np.array([1.0, 0.2])
        spark_model = AgiBotG1MobileBaseBicycleDynamic2Config().create_dynamics_model()

        np.testing.assert_allclose(
            pipeline.bicycle_model.step(state, control, config.dt, "native"),
            spark_model.step(state, control, config.dt, "Euler"),
            atol=1e-12,
        )

    def test_bicycle_conversion_remains_a_mathematical_placeholder(self):
        with self.assertRaisesRegex(NotImplementedError, "Problem 2.2"):
            problem.bicycle_initial_state(np.zeros(4))
        with self.assertRaisesRegex(NotImplementedError, "Problem 2.2"):
            problem.bicycle_control(np.zeros(6), Config())

    def test_stopped_bicycle_does_not_amplify_lateral_roundoff(self):
        config = Config(duration=0.1)
        with tempfile.TemporaryDirectory() as directory:
            model = ComparisonPipeline(config, directory).bicycle_model
        state = np.array([0.0, 0.0, 0.09, 1e-10, 1e-10, 0.0])

        next_state = model.step(state, np.zeros(2), config.dt, "native")

        np.testing.assert_array_equal(next_state[3:5], np.zeros(2))

    def test_simulation_reaches_the_correct_unfinished_controller(self):
        with tempfile.TemporaryDirectory() as directory:
            pipeline = ComparisonPipeline(Config(duration=0.02), directory)
            with self.assertRaisesRegex(NotImplementedError, "Problem 2.1"):
                pipeline.simulate("unicycle")
            with self.assertRaisesRegex(NotImplementedError, "Problem 2.2"):
                pipeline.simulate("bicycle")

    def test_viewer_decimates_playback_and_updates_simulation_time(self):
        viewer = object.__new__(ModelViewer)
        viewer.config = Config(duration=0.04)
        viewer.robot_config = AgiBotG1MobileBaseDynamic1Config()
        viewer.default_q = np.array(
            [viewer.robot_config.DefaultDoFVal[dof] for dof in viewer.robot_config.DoFs]
        )
        viewer.agent = _ViewerAgentSpy()
        with tempfile.TemporaryDirectory() as directory:
            model = ComparisonPipeline(viewer.config, directory).bicycle_model
        states = np.zeros((5, 6))

        viewer.replay(states, model)

        self.assertEqual(viewer.agent.reset_count, 3)
        initial_information, replace = viewer.agent.information[0]
        self.assertTrue(replace)
        self.assertEqual(initial_information["Dynamics model"], "bicycle_dynamic_6_state")
        self.assertEqual(initial_information["Robot config"], "DiscreteTimeDynamicsConfig")
        self.assertEqual(initial_information["Propagation"], "external forward Euler replay")
        self.assertEqual(viewer.agent.information[-1][0]["Simulation time"], "0.040 s")

    def test_models_run_separately_and_save_independent_results(self):
        config = Config(duration=0.02)
        # Exercise provided simulation/I/O with arbitrary constant test inputs.
        with tempfile.TemporaryDirectory() as directory, patch.object(
            problem, "unicycle_control", return_value=np.zeros(2)
        ), patch.object(
            problem, "bicycle_initial_state", return_value=np.zeros(6)
        ), patch.object(problem, "bicycle_control", return_value=np.zeros(2)):
            pipeline = ComparisonPipeline(config, directory)
            unicycle = pipeline.run("unicycle")
            bicycle = pipeline.run("bicycle")

            self.assertEqual(unicycle["states"].shape, (3, 4))
            self.assertEqual(bicycle["states"].shape, (3, 6))
            self.assertEqual(
                (unicycle["output_dir"] / "states.csv").read_text().splitlines()[0],
                "time,p_x,p_y,v,theta",
            )
            self.assertEqual(
                (bicycle["output_dir"] / "states.csv").read_text().splitlines()[0],
                "time,X,Y,v_x,v_y,r,psi",
            )
            self.assertEqual(len(list(unicycle["output_dir"].glob("*.png"))), 5)
            self.assertEqual(len(list(bicycle["output_dir"].glob("*.png"))), 7)

    def test_gain_study_saves_all_controller_configurations(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(
            problem, "unicycle_control", return_value=np.zeros(2)
        ), patch.object(
            problem, "bicycle_initial_state", return_value=np.zeros(6)
        ), patch.object(problem, "bicycle_control", return_value=np.zeros(2)):
            result = run_gain_study(Config(duration=0.02), directory)
            self.assertEqual(len(result["traces"]), 9)
            self.assertTrue((result["output_dir"] / "errors.csv").is_file())
            self.assertTrue((result["output_dir"] / "position_error.png").is_file())
            self.assertTrue((result["output_dir"] / "orientation_error.png").is_file())


if __name__ == "__main__":
    unittest.main()
