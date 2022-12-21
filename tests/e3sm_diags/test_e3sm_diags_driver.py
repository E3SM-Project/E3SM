import pytest

from e3sm_diags.e3sm_diags_driver import _run_serially, _run_with_dask, run_diag
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger("e3sm_diags.e3sm_diags_driver", propagate=True)


class TestRunDiag:
    def test_returns_parameter_with_results(self):
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]

        results = run_diag(parameter)
        expected = [parameter]

        # NOTE: We are only testing that the function returns a list of
        # parameter objects, not the results themselves. There are integration
        # tests validates the results.
        assert results == expected

    def test_logs_error_if_driver_module_for_set_not_found(self, caplog):
        parameter = CoreParameter()
        parameter.sets = ["invalid_set"]

        run_diag(parameter)

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
        # run_diag(parameter)

        assert "TypeError: 'NoneType' object is not iterable" in caplog.text

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
