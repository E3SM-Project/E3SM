from __future__ import print_function

import os

import cdms2

import acme_diags
from acme_diags.driver import utils
from acme_diags.plot import plot


def run_diag(parameter):
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
                acme_diags.INSTALL_PATH, "acme_ne30_ocean_land_mask.nc"
            )
            with cdms2.open(mask_path) as f:
                land_frac = f("LANDFRAC")
                ocean_frac = f("OCNFRAC")

        for var in variables:
            print("Variable: {}".format(var))
            # test = test_data.get_timeseries_variable(var)
            # ref = ref_data.get_timeseries_variable(var)
            test = test_data.get_climo_variable(var, season)
            ref = ref_data.get_climo_variable(var, season)

            parameter.var_id = var
            parameter.viewer_descr[var] = (
                test.long_name
                if hasattr(test, "long_name")
                else "No long_name attr in test data."
            )

            for region in regions:
                # print("Selected region: {}".format(region))

                test_domain = utils.general.select_region(
                    region, test, land_frac, ocean_frac, parameter
                )
                ref_domain = utils.general.select_region(
                    region, ref, land_frac, ocean_frac, parameter
                )

                parameter.output_file = "-".join([ref_name, var, season, region])
                parameter.main_title = str(
                    " ".join([var, "Diurnal Cycle ", season, region])
                )

                (
                    test_cmean,
                    test_amplitude,
                    test_maxtime,
                ) = utils.diurnal_cycle.composite_diurnal_cycle(test_domain, season)
                (
                    ref_cmean,
                    ref_amplitude,
                    ref_maxtime,
                ) = utils.diurnal_cycle.composite_diurnal_cycle(ref_domain, season)
                parameter.var_region = region
                plot(
                    parameter.current_set,
                    test_maxtime,
                    test_amplitude,
                    ref_maxtime,
                    ref_amplitude,
                    parameter,
                )
                utils.general.save_ncfiles(
                    parameter.current_set,
                    test_cmean,
                    ref_cmean,
                    None,
                    parameter,
                )
                utils.general.save_ncfiles(
                    parameter.current_set,
                    test_amplitude,
                    ref_amplitude,
                    None,
                    parameter,
                )
                utils.general.save_ncfiles(
                    parameter.current_set,
                    test_maxtime,
                    ref_maxtime,
                    None,
                    parameter,
                )

    return parameter
