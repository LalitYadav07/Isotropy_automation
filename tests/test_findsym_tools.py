from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from findsym_tools import standardize_cif_with_fallback


STRUCTURAL_CIF = textwrap.dedent(
    """
    # CIF file created by FINDSYM, version 7.1.3
    data_findsym-output
    _audit_creation_method FINDSYM
    _cell_length_a     5.6300000000
    _cell_length_b     5.6300000000
    _cell_length_c     5.6300000000
    _cell_angle_alpha  90.0000000000
    _cell_angle_beta   90.0000000000
    _cell_angle_gamma  90.0000000000
    _symmetry_space_group_name_H-M "F 4/m -3 2/m"
    _symmetry_Int_Tables_number 225
    _space_group.reference_setting '225:-F 4 2 3'
    _space_group.transform_Pp_abc a,b,c;0,0,0
    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_symmetry_multiplicity
    _atom_site_Wyckoff_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    Na1 Na   4 a  0.0000000000  0.0000000000  0.0000000000  1.0000000000
    Cl1 Cl   4 b  0.5000000000  0.5000000000  0.5000000000  1.0000000000
    """
).strip() + "\n"

MAGNETIC_CIF = textwrap.dedent(
    """
    # CIF file created by FINDSYM, version 7.1.3
    data_findsym-output
    _audit_creation_method FINDSYM
    _cell_length_a     3.9810111781
    _cell_length_b     3.9810111781
    _cell_length_c     5.6300000000
    _cell_angle_alpha  90.0000000000
    _cell_angle_beta   90.0000000000
    _cell_angle_gamma  90.0000000000
    _space_group_magn.number_BNS "140.550"
    _space_group_magn.name_UNI "I4/mcm.1'_c[rP4/mmm]"
    _space_group_magn.name_BNS "I_c4/mcm"
    loop_
    _space_group_symop_magn_operation.id
    _space_group_symop_magn_operation.xyz
    1 x,y,z,+1
    2 x,-y,-z+1/2,+1
    loop_
    _space_group_symop_magn_centering.id
    _space_group_symop_magn_centering.xyz
    1 x,y,z,+1
    2 x,y,z+1/2,-1
    3 x+1/2,y+1/2,z+1/2,+1
    4 x+1/2,y+1/2,z,-1
    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_symmetry_multiplicity
    _atom_site_Wyckoff_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    Fe1 Fe   4 a  0.0000000000  0.0000000000  0.0000000000  1.0000000000
    loop_
    _atom_site_moment.label
    _atom_site_moment.crystalaxis_x
    _atom_site_moment.crystalaxis_y
    _atom_site_moment.crystalaxis_z
    Fe1  0.0000000000  0.0000000000  1.0000000000
    """
).strip() + "\n"


class FindsymToolsMagneticTest(unittest.TestCase):
    def test_structural_standardization_reports_structural_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "nacl.cif"
            output = Path(tmpdir) / "nacl_standardized.cif"
            source.write_text(STRUCTURAL_CIF, encoding="utf-8")

            result = standardize_cif_with_fallback(source, output)

        self.assertEqual(result.method, "findsym")
        self.assertEqual(result.standardized_space_group_number, 225)
        self.assertEqual(result.standardized_space_group_symbol, "Fm-3m")
        self.assertIsNone(result.magnetic_space_group_number_bns)
        self.assertIsNone(result.magnetic_space_group_symbol_bns)
        self.assertFalse(result.magnetic_symmetry_evaluated)
        self.assertEqual(result.origin_shift, (0.0, 0.0, 0.0))
        self.assertEqual(result.basis_vectors, ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -0.0, 1.0)))

    def test_magnetic_standardization_reports_bns_group(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "afm.cif"
            output = Path(tmpdir) / "afm_standardized.cif"
            source.write_text(MAGNETIC_CIF, encoding="utf-8")

            result = standardize_cif_with_fallback(source, output)

        self.assertEqual(result.method, "findsym")
        self.assertEqual(result.standardized_space_group_number, 221)
        self.assertEqual(result.standardized_space_group_symbol, "Pm-3m")
        self.assertEqual(result.magnetic_space_group_number_bns, "140.550")
        self.assertEqual(result.magnetic_space_group_symbol_bns, "I_c4/mcm")
        self.assertTrue(result.magnetic_symmetry_evaluated)
        self.assertEqual(result.origin_shift, (0.0, 0.0, 0.5))
        self.assertEqual(result.basis_vectors, ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -0.0, 1.0)))


if __name__ == "__main__":
    unittest.main()
