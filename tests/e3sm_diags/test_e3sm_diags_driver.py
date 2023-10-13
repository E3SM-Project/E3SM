import pytest

from e3sm_diags.e3sm_diags_driver import _run_serially, _run_with_dask
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger("e3sm_diags.e3sm_diags_driver", propagate=True)


class TestRunDiag:
    def test_run_diag_serially_returns_parameters_with_results(self):
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]

        results = _run_serially([parameter])
        expected = [parameter]

        # NOTE: We are only testing that the function returns a list of
        # parameter objects, not the results themselves. There are integration
        # tests validates the results.
        assert results == expected

    def test_run_diag_with_dask_returns_parameters_with_results(self):
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]

        results = _run_with_dask([parameter])

        expected_parameter = CoreParameter()
        expected_parameter.sets = ["lat_lon"]
        expected_parameter.current_set = "lat_lon"
        expected_parameter.test_name_yrs = ""
        expected_parameter.ref_name_yrs = ""
        expected_parameter.model_only = False
        expected = [expected_parameter]

        # NOTE: We are only testing that the function returns a list of
        # parameter objects, not the results themselves. There are integration
        # tests validates the results.
        assert results[0].__dict__ == expected[0].__dict__

    def test_run_diag_with_dask_raises_error_if_num_workers_attr_not_set(
        self,
    ):
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]
        del parameter.num_workers

        with pytest.raises(ValueError):
            _run_with_dask([parameter])
