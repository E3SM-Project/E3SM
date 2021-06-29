from __future__ import print_function

from typing import TYPE_CHECKING, Any, Dict, List

import cdutil
import MV2
from cdms2 import createAxis

from acme_diags.driver.utils.dataset import Dataset
from acme_diags.driver.utils.general import (
    get_name_and_yrs,
    regrid_to_lower_res,
    save_ncfiles,
    select_region_lat_lon,
)
from acme_diags.plot import plot

if TYPE_CHECKING:
    from cdms2.axis import TransientAxis
    from cdms2.tvariable import TransientVariable

    from acme_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: "CoreParameter"):
    """Runs the annual cycle zonal mean diagnostic.

    :param parameter: Parameters for the run
    :type parameter: CoreParameter
    :return: Parameters for the run
    :rtype: CoreParameter
    """
    variables: List[str] = parameter.variables
    ref_name = getattr(parameter, "ref_name", "")

    test_data = Dataset(parameter, test=True)
    ref_data = Dataset(parameter, ref=True)

    parameter.test_name_yrs = get_name_and_yrs(parameter, test_data, "01")
    parameter.ref_name_yrs = get_name_and_yrs(parameter, ref_data, "01")

    for var in variables:
        test_ac = _create_annual_cycle(test_data, var)
        ref_ac = _create_annual_cycle(ref_data, var)

        test_ac_reg, ref_ac_reg = regrid_to_lower_res(
            test_ac,
            ref_ac,
            parameter.regrid_tool,
            parameter.regrid_method,
        )

        test_ac_zonal_mean = cdutil.averager(test_ac, axis="x", weights="generate")
        test_ac_reg_zonal_mean = cdutil.averager(
            test_ac_reg, axis="x", weights="generate"
        )

        if (
            parameter.ref_name == "OMI-MLS"
        ):  # SCO from OMI-MLS only available as (time, lat)
            test_ac_reg_zonal_mean = select_region_lat_lon(
                "60S60N", test_ac_reg_zonal_mean, parameter
            )
            test_ac_zonal_mean = select_region_lat_lon(
                "60S60N", test_ac_zonal_mean, parameter
            )
            if var == "SCO":
                ref_ac_zonal_mean = ref_ac
                ref_ac_reg_zonal_mean = ref_ac_reg
            else:
                ref_ac_zonal_mean = cdutil.averager(
                    ref_ac, axis="x", weights="generate"
                )
                ref_ac_reg_zonal_mean = cdutil.averager(
                    ref_ac_reg, axis="x", weights="generate"
                )

        else:
            ref_ac_zonal_mean = cdutil.averager(ref_ac, axis="x", weights="generate")
            ref_ac_reg_zonal_mean = cdutil.averager(
                ref_ac_reg, axis="x", weights="generate"
            )

        # if var == 'SCO' and parameter.ref_name=='OMI-MLS':  # SCO from OMI-MLS only available as (time, lat)
        #    ref_ac_zonal_mean = ref_ac
        #    ref_ac_reg_zonal_mean = ref_ac_reg

        #    test_ac_reg_zonal_mean = select_region_lat_lon("60S60N", test_ac_reg_zonal_mean, parameter)
        #    test_ac_zonal_mean = select_region_lat_lon("60S60N", test_ac_zonal_mean, parameter)
        # else:
        #    ref_ac_zonal_mean = cdutil.averager(ref_ac, axis="x", weights="generate")
        #    ref_ac_reg_zonal_mean = cdutil.averager(
        #        ref_ac_reg, axis="x", weights="generate"
        #    )

        diff_ac = test_ac_reg_zonal_mean - ref_ac_reg_zonal_mean
        diff_ac.setAxis(1, test_ac_reg_zonal_mean.getAxis(1))
        diff_ac.setAxis(0, test_ac_reg_zonal_mean.getAxis(0))

        parameter.var_id = var
        parameter.output_file = "-".join([ref_name, var, "Annual-Cycle"])
        parameter.main_title = str(" ".join([var, "Zonel Mean Annual Cycle"]))

        parameter.viewer_descr[var] = (
            test_ac.long_name
            if hasattr(test_ac, "long_name")
            else "No long_name attr in test data."
        )

        metrics_dict: Dict[str, Any] = {}

        plot(
            parameter.current_set,
            ref_ac_zonal_mean,
            test_ac_zonal_mean,
            diff_ac,
            metrics_dict,
            parameter,
        )
        save_ncfiles(
            parameter.current_set,
            ref_ac_zonal_mean,
            test_ac_zonal_mean,
            diff_ac,
            parameter,
        )

    return parameter


def _create_annual_cycle(dataset: Dataset, variable: str) -> "TransientVariable":
    """Creates the annual climatology cycle for a dataset variable.

    :param dataset: Dataset
    :type dataset: Dataset
    :param variable: Dataset variable
    :type variable: str
    :return: Variable's annual climatology cycle
    :rtype: tvariable.TransientVariable
    """
    months = range(1, 13)
    month_list = [f"{x:02}" for x in list(months)]

    for index, month in enumerate(month_list):
        var = dataset.get_climo_variable(variable, month)
        if month == "01":
            var_ann_cycle: "TransientVariable" = MV2.zeros([12] + list(var.shape))
            var_ann_cycle.id = var.id
            var_ann_cycle.long_name = var.long_name
            var_ann_cycle.units = var.units

            time: "TransientAxis" = createAxis(months)
            time.id = "time"

            var_ann_cycle.setAxis(0, time)
            time.designateTime()

            for iax in list(range(len(var.shape))):
                var_ann_cycle.setAxis(1 + iax, var.getAxis(iax))
            # var_ann_cycle.setAxis(2, var.getAxis(1))

            var_ann_cycle[0] = var
        else:
            var_ann_cycle[index] = var

    return var_ann_cycle
