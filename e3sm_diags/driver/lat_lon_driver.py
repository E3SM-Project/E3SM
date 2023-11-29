from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING

import cdms2

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics import corr, mean, rmse, std
from e3sm_diags.plot import plot

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def create_and_save_data_and_metrics(parameter, mv1_domain, mv2_domain):
    if not parameter.model_only:
        # Regrid towards the lower resolution of the two
        # variables for calculating the difference.
        mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
            mv1_domain,
            mv2_domain,
            parameter.regrid_tool,
            parameter.regrid_method,
        )

        diff = mv1_reg - mv2_reg
    else:
        mv2_domain = None
        mv2_reg = None
        mv1_reg = mv1_domain
        diff = None

    metrics_dict = create_metrics(mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)

    # Saving the metrics as a json.
    metrics_dict["unit"] = mv1_domain.units

    fnm = os.path.join(
        utils.general.get_output_dir(parameter.current_set, parameter),
        parameter.output_file + ".json",
    )
    with open(fnm, "w") as outfile:
        json.dump(metrics_dict, outfile)

    logger.info(f"Metrics saved in {fnm}")

    plot(
        parameter.current_set,
        mv2_domain,
        mv1_domain,
        diff,
        metrics_dict,
        parameter,
    )
    utils.general.save_ncfiles(
        parameter.current_set,
        mv1_domain,
        mv2_domain,
        diff,
        parameter,
    )


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    # For input None, metrics are instantiated to 999.999.
    # Apply float() to make sure the elements in metrics_dict are JSON serializable, i.e. np.float64 type is JSON serializable, but not np.float32.
    missing_value = 999.999
    metrics_dict = {}
    metrics_dict["ref"] = {
        "min": float(ref.min()) if ref is not None else missing_value,
        "max": float(ref.max()) if ref is not None else missing_value,
        "mean": float(mean(ref)) if ref is not None else missing_value,
    }
    metrics_dict["ref_regrid"] = {
        "min": float(ref_regrid.min()) if ref_regrid is not None else missing_value,
        "max": float(ref_regrid.max()) if ref_regrid is not None else missing_value,
        "mean": float(mean(ref_regrid)) if ref_regrid is not None else missing_value,
        "std": float(std(ref_regrid)) if ref_regrid is not None else missing_value,
    }
    metrics_dict["test"] = {
        "min": float(test.min()),
        "max": float(test.max()),
        "mean": float(mean(test)),
    }
    metrics_dict["test_regrid"] = {
        "min": float(test_regrid.min()),
        "max": float(test_regrid.max()),
        "mean": float(mean(test_regrid)),
        "std": float(std(test_regrid)),
    }
    metrics_dict["diff"] = {
        "min": float(diff.min()) if diff is not None else missing_value,
        "max": float(diff.max()) if diff is not None else missing_value,
        "mean": float(mean(diff)) if diff is not None else missing_value,
    }
    metrics_dict["misc"] = {
        "rmse": float(rmse(test_regrid, ref_regrid))
        if ref_regrid is not None
        else missing_value,
        "corr": float(corr(test_regrid, ref_regrid))
        if ref_regrid is not None
        else missing_value,
    }
    return metrics_dict


def run_diag(parameter: CoreParameter) -> CoreParameter:  # noqa: C901
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)

    for season in seasons:
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = utils.general.get_name_and_yrs(
            parameter, test_data, season
        )
        parameter.ref_name_yrs = utils.general.get_name_and_yrs(
            parameter, ref_data, season
        )

        # Get land/ocean fraction for masking.
        try:
            land_frac = test_data.get_climo_variable("LANDFRAC", season)
            ocean_frac = test_data.get_climo_variable("OCNFRAC", season)
        except Exception:
            mask_path = os.path.join(
                e3sm_diags.INSTALL_PATH, "acme_ne30_ocean_land_mask.nc"
            )
            with cdms2.open(mask_path) as f:
                land_frac = f("LANDFRAC")
                ocean_frac = f("OCNFRAC")

        parameter.model_only = False
        for var in variables:
            logger.info("Variable: {}".format(var))
            parameter.var_id = var

            mv1 = test_data.get_climo_variable(var, season)
            try:
                mv2 = ref_data.get_climo_variable(var, season)
            except (RuntimeError, IOError):
                mv2 = mv1
                logger.info("Can not process reference data, analyse test data only")

                parameter.model_only = True

            parameter.viewer_descr[var] = (
                mv1.long_name
                if hasattr(mv1, "long_name")
                else "No long_name attr in test data."
            )

            # For variables with a z-axis.
            if mv1.getLevel() and mv2.getLevel():
                plev = parameter.plevs
                logger.info("Selected pressure level: {}".format(plev))

                mv1_p = utils.general.convert_to_pressure_levels(
                    mv1, plev, test_data, var, season
                )
                mv2_p = utils.general.convert_to_pressure_levels(
                    mv2, plev, ref_data, var, season
                )

                # Select plev.
                for ilev in range(len(plev)):
                    mv1 = mv1_p[ilev,]
                    mv2 = mv2_p[ilev,]

                    for region in regions:
                        parameter.var_region = region
                        logger.info(f"Selected regions: {region}")
                        mv1_domain = utils.general.select_region(
                            region, mv1, land_frac, ocean_frac, parameter
                        )
                        mv2_domain = utils.general.select_region(
                            region, mv2, land_frac, ocean_frac, parameter
                        )

                        parameter.output_file = "-".join(
                            [
                                ref_name,
                                var,
                                str(int(plev[ilev])),
                                season,
                                region,
                            ]
                        )
                        parameter.main_title = str(
                            " ".join(
                                [
                                    var,
                                    str(int(plev[ilev])),
                                    "mb",
                                    season,
                                    region,
                                ]
                            )
                        )

                        create_and_save_data_and_metrics(
                            parameter, mv1_domain, mv2_domain
                        )

            # For variables without a z-axis.
            elif mv1.getLevel() is None and mv2.getLevel() is None:
                for region in regions:
                    parameter.var_region = region

                    logger.info(f"Selected region: {region}")
                    mv1_domain = utils.general.select_region(
                        region, mv1, land_frac, ocean_frac, parameter
                    )
                    mv2_domain = utils.general.select_region(
                        region, mv2, land_frac, ocean_frac, parameter
                    )

                    parameter.output_file = "-".join([ref_name, var, season, region])
                    parameter.main_title = str(" ".join([var, season, region]))

                    create_and_save_data_and_metrics(parameter, mv1_domain, mv2_domain)

            else:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting."
                )

    return parameter
