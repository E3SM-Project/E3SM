from typing import TYPE_CHECKING
from unittest import TestCase
from unittest.mock import Mock

import numpy as np

from acme_diags.derivations.acme import adjust_prs_val_units, determine_cloud_level

if TYPE_CHECKING:
    from cdms2.axis import FileAxis


class TestAdjustPrsValUnits(TestCase):
    def setUp(self):
        self.mock_file_axis: "FileAxis" = Mock()
        self.mock_file_axis.getData.return_value = np.arange(0, 1002)

        self.prs_val = 1
        self.prs_val0 = 10

    def test_without_swapping_units(self):
        actual = adjust_prs_val_units(self.mock_file_axis, self.prs_val, None)
        expected = 1

        self.assertEqual(actual, expected)

    def test_swap_units_without_multiplier_if_no_matched_id(self):
        self.mock_file_axis.id = "no_match"

        actual = adjust_prs_val_units(self.mock_file_axis, self.prs_val, self.prs_val0)
        expected = 10
        self.assertEqual(actual, expected)

    def test_swap_units_without_multiplier_if_max_less_than_1000(
        self,
    ):
        self.mock_file_axis.id = "cosp_prs"
        self.mock_file_axis.getData.return_value = np.arange(0, 1)

        actual = adjust_prs_val_units(self.mock_file_axis, self.prs_val, self.prs_val0)
        expected = 10
        self.assertEqual(actual, expected)

    def test_swap_units_and_apply_multiplier_for_all_matched_ids(self):
        adjust_ids = {"cosp_prs": 100, "cosp_htmsir": 1000}

        for id, multiplier in adjust_ids.items():
            self.mock_file_axis.id = id

            actual = adjust_prs_val_units(
                self.mock_file_axis, self.prs_val, self.prs_val0
            )
            expected = self.prs_val0 * multiplier
            self.assertEqual(actual, expected)


class TestDetermineCloudLevel(TestCase):
    def setUp(self):
        self.low_bnds = (1, 2)
        self.high_bdns = (3, 4)

    def test_returns_middle_cloud_fraction(self):
        actual = determine_cloud_level(1, 3, self.low_bnds, self.high_bdns)
        expected = "middle cloud fraction"
        self.assertEqual(actual, expected)

    def test_returns_high_cloud_fraction(self):
        actual = determine_cloud_level(1, 1, self.low_bnds, self.high_bdns)
        expected = "high cloud fraction"
        self.assertEqual(actual, expected)

    def test_returns_low_cloud_fraction(self):
        actual = determine_cloud_level(3, 3, self.low_bnds, self.high_bdns)
        expected = "low cloud fraction"
        self.assertEqual(actual, expected)

    def test_returns_total_cloud_fraction(self):
        actual = determine_cloud_level(0, 5, self.low_bnds, self.high_bdns)
        expected = "total cloud fraction"
        self.assertEqual(actual, expected)
