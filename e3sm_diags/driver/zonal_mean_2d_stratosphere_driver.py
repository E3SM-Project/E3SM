from e3sm_diags.driver.zonal_mean_2d_driver import create_metrics as base_create_metrics
from e3sm_diags.driver.zonal_mean_2d_driver import run_diag as base_run_diag
from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from e3sm_diags.parameter.zonal_mean_2d_stratosphere_parameter import (
    ZonalMean2dStratosphereParameter,
)


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    return base_create_metrics(ref, test, ref_regrid, test_regrid, diff)


def run_diag(
    parameter: ZonalMean2dStratosphereParameter,
) -> ZonalMean2dParameter:
    return base_run_diag(
        parameter, default_plevs=ZonalMean2dStratosphereParameter().plevs
    )
