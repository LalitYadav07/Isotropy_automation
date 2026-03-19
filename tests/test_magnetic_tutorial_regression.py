import json
import unittest
from pathlib import Path

from scripts.isodistort_distortion_parser import parse_distortion_file


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = REPO_ROOT / "fixtures" / "magnetic_tutorial_regression.json"


class MagneticTutorialRegressionTests(unittest.TestCase):
    maxDiff = None

    @classmethod
    def setUpClass(cls):
        cls.fixtures = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))

    def test_all_magnetic_tutorial_fixtures_match_parser_output(self):
        for fixture in self.fixtures:
            with self.subTest(example=fixture["id"]):
                data = parse_distortion_file(REPO_ROOT / fixture["distortion_file"])
                expectations = fixture["parser_expectations"]

                for label in expectations.get("classification_includes", []):
                    self.assertIn(label, data.classification)

                if "order_parameter_equals" in expectations:
                    self.assertEqual(expectations["order_parameter_equals"], data.order_parameter)

                for snippet in expectations.get("order_parameter_contains", []):
                    self.assertIn(snippet, data.order_parameter)

                for filename in expectations.get("associated_files_include", []):
                    self.assertIn(filename, data.associated_files)

                self.assertEqual(
                    expectations.get("child_basic_structure_atom_count"),
                    None if data.child_basic_structure is None else data.child_basic_structure["atom_count"],
                )

                expected_primary_irreps = expectations.get("primary_irreps")
                if expected_primary_irreps == []:
                    self.assertEqual([], data.primary_irreps)
                elif expected_primary_irreps is not None:
                    self.assertEqual(len(expected_primary_irreps), len(data.primary_irreps))
                    for expected, observed in zip(expected_primary_irreps, data.primary_irreps, strict=True):
                        self.assertEqual(expected["irrep"], observed["irrep"])
                        self.assertEqual(expected["k_vector"], observed["k_vector"])
                        self.assertEqual(
                            expected["magnetic_nonzero"],
                            int(observed["category_summaries"]["magnetic"]["nonzero_count"]),
                        )

                for family, count in expectations.get("mode_nonzero_counts", {}).items():
                    self.assertEqual(count, data.mode_summaries[family]["nonzero_mode_count"])

                for key in expectations.get("internal_checks_true", []):
                    self.assertTrue(data.internal_checks[key], msg=f"expected internal check {key} to be true")

                for key in expectations.get("internal_checks_false", []):
                    self.assertFalse(data.internal_checks[key], msg=f"expected internal check {key} to be false")


if __name__ == "__main__":
    unittest.main()
