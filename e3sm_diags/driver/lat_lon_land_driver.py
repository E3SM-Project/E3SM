from __future__ import annotations

from typing import TYPE_CHECKING

from e3sm_diags.driver.lat_lon_driver import (
    create_and_save_data_and_metrics as base_create_and_save_data_and_metrics,
)
from e3sm_diags.driver.lat_lon_driver import create_metrics as base_create_metrics
from e3sm_diags.driver.lat_lon_driver import run_diag as base_run_diag

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter
    from e3sm_diags.parameter.lat_lon_land_parameter import LatLonLandParameter


def create_and_save_data_and_metrics(parameter, test, ref):
    return base_create_and_save_data_and_metrics(parameter, test, ref)


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    return base_create_metrics(ref, test, ref_regrid, test_regrid, diff)


def run_diag(parameter: LatLonLandParameter) -> CoreParameter:
    return base_run_diag(parameter)
