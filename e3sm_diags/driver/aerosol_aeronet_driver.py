from __future__ import annotations

import os
from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd
from scipy import interpolate

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.logger import custom_logger
from e3sm_diags.plot.cartopy import aerosol_aeronet_plot

if TYPE_CHECKING:
    from cdms2.tvariable import TransientVariable

    from e3sm_diags.parameter.core_parameter import CoreParameter


logger = custom_logger(__name__)

# This aerosol diagnostics scripts based on AERONET sites data was originally developed by Feng Yan and adapted and integrated in e3sm_diags by Jill Zhang.
# Years include 2006–2015 average climatology for observation according to Feng et al. 2022:doi:10.1002/essoar.10510950.1, and Golaz et al. 2022 E3SMv2 paper.


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Runs the aerosol aeronet diagnostic.

    :param parameter: Parameters for the run
    :type parameter: CoreParameter
    :raises ValueError: Invalid run type
    :return: Parameters for the run
    :rtype: CoreParameter
    """
    variables = parameter.variables
    run_type = parameter.run_type
    seasons = parameter.seasons

    for season in seasons:
        test_data = utils.dataset.Dataset(parameter, test=True)
        parameter.test_name_yrs = utils.general.get_name_and_yrs(
            parameter, test_data, season
        )
        parameter.ref_name_yrs = "AERONET (2006-2015)"

        for var in variables:
            logger.info("Variable: {}".format(var))
            parameter.var_id = var

            test = test_data.get_climo_variable(var, season)
            test_site = interpolate_model_output_to_obs_sites(test, var)

        if run_type == "model_vs_model":
            ref_data = utils.dataset.Dataset(parameter, ref=True)
            parameter.ref_name_yrs = utils.general.get_name_and_yrs(
                parameter, ref_data, season
            )
            ref = ref_data.get_climo_variable(var, season)
            ref_site = interpolate_model_output_to_obs_sites(ref, var)

        elif run_type == "model_vs_obs":
            ref_site = interpolate_model_output_to_obs_sites(None, var)
        else:
            raise ValueError("Invalid run_type={}".format(run_type))

        parameter.output_file = (
            f"{parameter.ref_name}-{parameter.var_id}-{season}-global"
        )
        aerosol_aeronet_plot.plot(test, test_site, ref_site, parameter)

    return parameter


def interpolate_model_output_to_obs_sites(
    var: Optional[TransientVariable], var_id: str
):
    """Interpolate model outputs (on regular lat lon grids) to observational sites

    :param var: Input model variable, var_id: name of the variable
    :type var: TransientVariable or NoneType, var_id: str
    :raises IOError: Invalid variable input
    :return: interpolated values over all observational sites
    :rtype: 1-D numpy.array

    """
    logger.info(
        "Interpolate model outputs (on regular lat lon grids) to observational sites"
    )
    if var_id == "AODABS":
        aeronet_file = os.path.join(
            e3sm_diags.INSTALL_PATH, "aerosol_aeronet/aaod550_AERONET_2006-2015.txt"
        )
        var_header = "aaod"
    elif var_id == "AODVIS":
        aeronet_file = os.path.join(
            e3sm_diags.INSTALL_PATH, "aerosol_aeronet/aod550_AERONET_2006-2015.txt"
        )
        var_header = "aod"
    else:
        raise IOError("Invalid variable input.")

    data_obs = pd.read_csv(aeronet_file, dtype=object, sep=",")

    lonloc = np.array(data_obs["lon"].astype(float))
    latloc = np.array(data_obs["lat"].astype(float))
    obsloc = np.array(data_obs[var_header].astype(float))
    # sitename = np.array(data_obs["site"].astype(str))
    nsite = len(obsloc)

    # express lonloc from 0 to 360
    lonloc[lonloc < 0.0] = lonloc[lonloc < 0.0] + 360.0

    if var is not None:
        f_intp = interpolate.RectBivariateSpline(
            var.getLatitude()[:], var.getLongitude()[:], var
        )
        var_intp = np.zeros(nsite)
        for i in range(nsite):
            var_intp[i] = f_intp(latloc[i], lonloc[i])

        return var_intp
    return obsloc
