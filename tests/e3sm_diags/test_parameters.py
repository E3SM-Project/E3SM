import numpy as np
import pytest

from e3sm_diags.parameter.annual_cycle_zonal_mean_parameter import ACzonalmeanParameter
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.lat_lon_land_parameter import LatLonLandParameter
from e3sm_diags.parameter.lat_lon_river_parameter import LatLonRiverParameter
from e3sm_diags.parameter.meridional_mean_2d_parameter import MeridionalMean2dParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.parameter.time_series_parameter import TimeSeriesParameter
from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from e3sm_diags.parameter.zonal_mean_2d_stratosphere_parameter import (
    ZonalMean2dStratosphereParameter,
)


class TestCoreParameter:
    def test__init__(self):
        CoreParameter()

    def test__add__copies_attributes_from_other_object(self):
        param1 = CoreParameter()
        param2 = CoreParameter()

        # Add custom attributes to the second object
        param2.test_start_yr = 2000  # type: ignore
        param2.test_end_yr = 2001  # type: ignore

        new_param = param2 + param1

        assert new_param.test_start_yr == 2000
        assert new_param.test_end_yr == 2001

    def test_check_values_does_not_raise_error_if_required_args_are_set(self):
        param = CoreParameter()
        param.reference_data_path = "path"
        param.test_data_path = "path"
        param.results_dir = "path"

        param.check_values()

    def test_check_values_raises_error_if_required_args_are_not_set(self):
        param = CoreParameter()

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_raises_error_if_ref_timeseries_input_and_no_ref_start_and_end_year_set(
        self,
    ):
        param = CoreParameter()
        param.reference_data_path = "path"
        param.test_data_path = "path"
        param.results_dir = "path"

        param.ref_timeseries_input = True

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_raises_error_if_test_timeseries_input_and_no_test_start_and_end_year_set(
        self,
    ):
        param = CoreParameter()
        param.reference_data_path = "path"
        param.test_data_path = "path"
        param.results_dir = "path"

        param.test_timeseries_input = True

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_returns_parameter_with_results(self):
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]

        results = parameter._run_diag()
        expected = [parameter]

        # NOTE: We are only testing that the function returns a list of
        # parameter objects, not the results themselves. There are integration
        # tests validates the results.
        assert results == expected

    def test_logs_error_if_driver_module_for_set_not_found(self, caplog):
        parameter = CoreParameter()
        parameter.sets = ["invalid_set"]

        parameter._run_diag()

        assert (
            "ModuleNotFoundError: No module named 'e3sm_diags.driver.invalid_set_driver'"
            in caplog.text
        )

    @pytest.mark.xfail
    def test_logs_exception_if_driver_run_diag_function_fails(self, caplog):
        # TODO: Need to implement this test by raising an exception through
        # the driver's `run_diag` function
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]

        # Make this attribute an invalid value to test exception is thrown.
        parameter.seasons = None  # type: ignore

        # NOTE: Comment out temporarily to avoid polluting the test output log.
        # parameter._run_diag()

        assert "TypeError: 'NoneType' object is not iterable" in caplog.text


def test_ac_zonal_mean_parameter():
    param = ACzonalmeanParameter()

    assert "seasons" not in param.granulate


def test_area_mean_time_series_parameter(caplog):
    param = AreaMeanTimeSeriesParameter()

    assert "regions" not in param.granulate
    assert "seasons" not in param.granulate

    # Raises error if no start year or end year set.
    with pytest.raises(RuntimeError):
        param.check_values()

    assert "You have no value for ref_names. Calculate test data only" in caplog.text


def test_arms_diags_parameter():
    param = ARMDiagsParameter()

    assert "seasons" not in param.granulate


def test_diurnal_cycle_parameter():
    DiurnalCycleParameter()


def test_enso_diags_parameter():
    param1 = EnsoDiagsParameter()

    assert "seasons" not in param1.granulate

    # Raises error if nino region is not valid.
    param1.nino_region = "INVALID_NINO"
    with pytest.raises(RuntimeError):
        param1.check_values()

    param2 = EnsoDiagsParameter()

    # Raises error invalid plot type
    param2.plot_type = "invalid_plot_Type"
    with pytest.raises(RuntimeError):
        param2.check_values()


def test_lat_lon_land_parameter():
    LatLonLandParameter()


def test_lat_lon_river_parameter():
    LatLonRiverParameter()


class TestMeridionalMean2dParameter:
    def test__init__(self):
        param = MeridionalMean2dParameter()

        assert "plevs" not in param.granulate
        np.testing.assert_allclose(param.plevs, np.logspace(2.0, 3.0, num=17))

    def test_check_values_raises_error_if_plevs_is_not_a_list(self):
        param = MeridionalMean2dParameter()
        param.plevs = "not_an_array"  # type: ignore

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_raises_error_if_plevs_are_not_monotonically_increasing_or_decreasing(
        self,
    ):
        param = MeridionalMean2dParameter()
        param.plevs = [1, 0, 2]

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_raises_error_if_less_than_2_plevs_are_set(self):
        param = MeridionalMean2dParameter()
        param.plevs = [1]

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_forces_monotonically_increasing_list(self):
        param = MeridionalMean2dParameter()
        param.plevs = [3, 2, 1]

        param.check_values()

        assert param.plevs == [1, 2, 3]


def test_qbo_parameter():
    param = QboParameter()

    assert "seasons" not in param.granulate


def test_streamflow_parameter():
    param = StreamflowParameter()

    assert "seasons" not in param.granulate


def test_tc_analysis_parameter():
    param = TCAnalysisParameter()

    assert "seasons" not in param.granulate


class TestTimeSeriesParameter:
    def test__init__(self):
        TimeSeriesParameter()

    @pytest.mark.xfail
    def test_check_values_raises_error_if_start_year_set_with_no_end_year(self):
        assert 0

    @pytest.mark.xfail
    def test_check_values_raises_error_if_end_year_set_with_no_start_year(self):
        assert 0


class TestZonalMean2dParameter:
    def test__init__(self):
        param = ZonalMean2dParameter()

        assert "plevs" not in param.granulate
        np.testing.assert_allclose(param.plevs, np.linspace(50, 1000, 20))

    def test_check_values_raises_error_if_plevs_is_not_a_list(self):
        param = ZonalMean2dParameter()
        param.plevs = "not_an_array"  # type: ignore

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_raises_error_if_plevs_are_not_monotonically_increasing_or_decreasing(
        self,
    ):
        param = ZonalMean2dParameter()
        param.plevs = [1, 0, 2]

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_raises_error_if_less_than_2_plevs_are_set(self):
        param = ZonalMean2dParameter()
        param.plevs = [1]

        with pytest.raises(RuntimeError):
            param.check_values()

    def test_check_values_forces_monotonically_increasing_list(self):
        param = ZonalMean2dParameter()
        param.plevs = [3, 2, 1]

        param.check_values()

        assert param.plevs == [1, 2, 3]


def test_zonal_mean_2d_stratosphere_parameter():
    param = ZonalMean2dStratosphereParameter()
    np.testing.assert_allclose(param.plevs, np.logspace(0, 2.0, num=10))
