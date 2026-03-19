from pathlib import Path
import unittest

from scripts.isodistort_distortion_parser import parse_distortion_file


FIXTURES = Path("isotutorials/exercises-completed")


class DistortionParserMagneticMetadataTests(unittest.TestCase):
    def test_dymn6ge6_includes_structured_magnetic_metadata(self) -> None:
        data = parse_distortion_file(FIXTURES / "dymn6ge6-distortion.txt")

        self.assertEqual(data.magnetic_mode_metadata["mode_count"], 7)
        self.assertTrue(data.magnetic_mode_metadata["has_amp_phase_table"])
        self.assertEqual(len(data.superspace_metadata["incommensurate_irreps"]), 2)

        first_mode = data.magnetic_mode_metadata["modes"][0]
        self.assertEqual(first_mode["label"], "m1")
        self.assertEqual(first_mode["irrep_label"], "mGM2+, mk16t3")
        self.assertEqual(first_mode["irrep_index"], 1)
        self.assertEqual(len(first_mode["amp_phase_rows"]), 13)

    def test_tbmno3_secondary_channels_link_to_magnetic_primary_summary(self) -> None:
        data = parse_distortion_file(FIXTURES / "tbmno3-distortion2.txt")

        self.assertEqual(len(data.magnetic_secondary_linkages), 2)
        first_linkage = data.magnetic_secondary_linkages[0]
        families = {item["family"] for item in first_linkage["secondary_channels"]}
        self.assertIn("displacive", families)
        self.assertTrue(all("linkage_scope" in item for item in first_linkage["secondary_channels"]))


if __name__ == "__main__":
    unittest.main()
