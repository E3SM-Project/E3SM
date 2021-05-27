import shutil
import tempfile
from typing import TYPE_CHECKING, Dict, Union
from unittest import TestCase
from unittest.mock import MagicMock

import numpy as np
from netCDF4 import Dataset

from acme_diags.driver.tc_analysis_driver import (
    BASIN_DICT,
    _calc_mean_ace,
    _calc_num_storms_and_max_len,
    _calc_seasonal_cycle,
    _calc_ts_intensity_dist,
    _get_mon_wind,
    _get_monthmc_yearic,
    _get_vars_from_te_stitch,
)

if TYPE_CHECKING:
    from cdms2.dataset import DatasetVariable


class TestRunDiags(TestCase):
    # TODO: Add tests
    pass


class TestGenerateTCMetricsFromTEStitchFile(TestCase):
    # TODO: Add tests
    pass


class TestCalcNumStormAndMaxLen(TestCase):
    def test_correct_output(self):
        lines = ["s", "a", "a", "s"]

        expected_num_storms = 2
        expected_max_len = 2
        num_storms, max_len = _calc_num_storms_and_max_len(lines)

        self.assertEqual(num_storms, expected_num_storms)
        self.assertEqual(max_len, expected_max_len)


class TestGetVarsFromTEStitch(TestCase):
    def test_correct_output(self):
        lines = [
            "storm\tfoo\t90\t180\tfoo\t1\t1\t1",
            "test\tfoo\t90\t180\tfoo\t1\t1\t1",
        ]
        max_len = 2
        num_storms = 2

        expected: Dict[str, Union[np.ndarray, int]] = {
            "longmc": np.array([[90, np.nan]]),
            "latmc": np.array([[180, np.nan]]),
            "vsmc": np.array([[1.94, np.nan]]),
            "yearmc": np.array([[1, np.nan]]),
            "monthmc": np.array([[1, np.nan]]),
            "year_start": 90,
            "year_end": 90,
            "num_years": 1,
        }
        result = _get_vars_from_te_stitch(lines, max_len, num_storms)

        np.array_equal(result["longmc"], expected["longmc"])
        np.array_equal(result["latmc"], expected["latmc"])
        np.array_equal(result["vsmc"], expected["vsmc"])
        np.array_equal(result["yearmc"], expected["yearmc"])
        np.array_equal(result["monthmc"], expected["monthmc"])
        self.assertEqual(result["year_start"], expected["year_start"])
        self.assertEqual(result["year_end"], expected["year_end"])
        self.assertEqual(result["num_years"], expected["num_years"])


class TestDeriveMetricsPerBasin(TestCase):
    def test_correct_output(self):
        te_stitch_vars: Dict[str, Union[np.ndarray, int]] = {  # noqa
            "longmc": np.array([[90, np.nan]]),
            "latmc": np.array([[180, np.nan]]),
            "vsmc": np.array([[1.94, np.nan]]),
            "yearmc": np.array([[1, np.nan]]),
            "monthmc": np.array([[1, np.nan]]),
            "year_start": 90,
            "year_end": 90,
            "num_years": 1,
        }
        num_storms = 2  # noqa
        basin_info = BASIN_DICT["NA"]  # noqa

        # FIXME: Figure out how to add return value when calling a DatasetVariable
        # For example, ocnfranc[0,0]
        ocnfranc: "DatasetVariable" = MagicMock()
        ocnfranc.getLatitude.return_value = np.array([-90, -45, 0, 45, 90])
        ocnfranc.getLongitude.return_value = np.array([-180, -90, 0, 90, 180])

        # TODO: Add an expected value
        expected = None  # noqa
        # result = _derive_metrics_per_basin(num_storms, te_stitch_vars, ocnfranc, basin_info)  # noqa

        # self.assertEqual(result, expected)


class TestGenerateTCMetricsFromObsFiles(TestCase):
    # TODO: Add tests
    def setUp(self):
        # Variables
        # wmo_wind = None
        # time = np.ma.array(
        #     [[-3000.01, -5000.02], [50000.01, 50000.02]],
        #     mask=[[False, False], [True, False]],
        # )

        self.test_dir = tempfile.mkdtemp()
        nc = Dataset("IBTrACS.NA.v04r00.nc", "w", format="NETCDF4")
        nc.createDimension("time", None)
        nc.createVariable("time", "f8", ("time",))

    def tearDown(self):
        shutil.rmtree(self.test_dir)


class TestGetMonthsMCYearsIC(TestCase):
    def test_extracts_vars(self):
        time = np.ma.array(
            [[-3000.01, -5000.02], [50000.01, 50000.02]],
            mask=[[False, False], [True, False]],
        )

        expected_monthmc = np.array([8, 3])
        expected_yearic = np.array([1850, 1995])
        monthmc, yearic = _get_monthmc_yearic(time)

        np.array_equal(monthmc, expected_monthmc)
        np.array_equal(yearic, expected_yearic)


class TestGetMonWnd(TestCase):
    def test_extracts_vars(self):
        monthmc = np.array([[8, 8], [3, 3]])
        yearic = np.array([1850, 1995])
        vsmc = np.ma.array(
            [[80, 0], [50, 30]],
            mask=[[False, True], [False, False]],
        )

        # The values for 1850 should be ignored since it is outside of the obs year range
        expected_mon = [3]
        expected_wind = [50]
        mon, wind = _get_mon_wind(vsmc, monthmc, yearic, num_rows=vsmc.shape[0])

        self.assertEqual(mon, expected_mon)
        self.assertEqual(wind, expected_wind)


class TestCalcMeanAce(TestCase):
    def test_calc_mean_ace(self):
        vsmc = np.ma.array(
            [[80, 0], [50, 30]],
            mask=[[False, True], [False, False]],
        )
        yearic = np.array([1850, 1995])
        num_rows = vsmc.shape[0]

        expected = 0.00625
        result = _calc_mean_ace(vsmc, yearic, num_rows)
        self.assertAlmostEqual(result, expected)


class TestCalcTCIntensity(TestCase):
    def setUp(self):
        self.num_categories = 6

    def test_distribution_properly_allocates_counts(self):
        wind_speeds = [35, 64, 83, 96, 114, 137]
        expected = np.ones(self.num_categories)

        result = _calc_ts_intensity_dist(wind_speeds)
        np.array_equal(result, expected)

    def test_distribution_ignores_speeds_below_min_threshold(self):
        wind_speeds = [34]
        expected = np.zeros(len(wind_speeds))

        result = _calc_ts_intensity_dist(wind_speeds)
        np.array_equal(result, expected)


class TestCalcSeasonalCycle(TestCase):
    def test_calculates_seasonal_cycle(self):
        mon = [1, 2, 3, 3, 5]
        expected = np.array(
            [0.2, 0.2, 0.4, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        result = _calc_seasonal_cycle(mon)

        np.array_equal(result, expected)
