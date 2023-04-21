from __future__ import annotations

import csv
import json
import os
from typing import TYPE_CHECKING  # , Optional

import e3sm_diags
from e3sm_diags.driver import utils

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter

import cdutil
import numpy

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)

# This aerosol budget set is requested by the E3SM Aerosol Working Group. The script is integrated in e3sm_diags by Jill Zhang, with input from Kai Zhang, Taufiq Hassan, Xue Zheng, Ziming Ke, Susannah Burrows, and Naser Mahfouz.


def global_integral(var, area_m2):
    """Compute global integral of 2 dimentional properties"""
    return numpy.sum(numpy.sum(abs(var) * area_m2, axis=0), axis=0)


def calc_column_integral(data, aerosol, season):
    """Calculate column integrated mass"""
    mass = data.get_climo_variable(f"Mass_{aerosol}", season)
    hyai, hybi, ps = data.get_extra_variables_only(
        f"Mass_{aerosol}", season, extra_vars=["hyai", "hybi", "PS"]
    )

    p0 = 100000.0  # Pa
    ps = ps  # Pa
    pressure_levs = cdutil.vertical.reconstructPressureFromHybrid(ps, hyai, hybi, p0)

    # (72,lat,lon)
    delta_p = numpy.diff(pressure_levs, axis=0)
    mass_3d = mass * delta_p / 9.8  # mass density * mass air   kg/m2
    burden = numpy.nansum(mass_3d, axis=0)  # kg/m2
    return burden


def generate_metrics_dic(data, aerosol, season):
    metrics_dict = {}
    wetdep = data.get_climo_variable(f"{aerosol}_SFWET", season)
    drydep = data.get_climo_variable(f"{aerosol}_DDF", season)
    srfemis = data.get_climo_variable(f"SF{aerosol}", season)
    if aerosol in ["bc", "pom", "so4"]:
        elvemis = data.get_climo_variable(f"{aerosol}_CLXF", season)
    area = data.get_extra_variables_only(f"{aerosol}_DDF", season, extra_vars=["area"])
    area_m2 = area * REARTH**2

    burden = calc_column_integral(data, aerosol, season)
    burden_total = global_integral(burden, area_m2) * 1e-9  # kg to Tg
    sink = global_integral((drydep - wetdep), area_m2) * UNITS_CONV
    drydep = global_integral(drydep, area_m2) * UNITS_CONV
    wetdep = global_integral(wetdep, area_m2) * UNITS_CONV
    srfemis = global_integral(srfemis, area_m2) * UNITS_CONV
    if aerosol in ["bc", "pom", "so4"]:
        elvemis = global_integral(elvemis, area_m2) * UNITS_CONV
    else:
        elvemis = 0.0
    metrics_dict = {
        "Surface Emission (Tg/yr)": f"{srfemis:.2f}",
        "Elevated Emission (Tg/yr)": f"{elvemis:.2f}",
        "Sink (Tg/yr)": f"{sink:.2f}",
        "Dry Deposition (Tg/yr)": f"{drydep:.2f}",
        "Wet Deposition (Tg/yr)": f"{wetdep:.2f}",
        "Burden (Tg)": f"{burden_total:.2f}",
        "Lifetime (Days)": f"{burden_total/sink*365:.2f}",
    }
    return metrics_dict


REARTH = 6.37122e6  # km
UNITS_CONV = 86400.0 * 365.0 * 1e-9  # kg/s to Tg/yr

# species = ["bc", "dst", "mom", "ncl", "pom", "so4", "soa"]
SPECIES_NAMES = {
    "bc": "Black Carbon",
    "dst": "Dust",
    "mom": "Marine Organic Matter",
    "ncl": "Sea Salt",
    "pom": "Primary Organic Matter",
    "so4": "Sulfate",
    "soa": "Secondary Organic Aerosol",
}
MISSING_VALUE = 999.999


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Runs the aerosol aeronet diagnostic.

    :param parameter: Parameters for the run
    :type parameter: CoreParameter
    :raises ValueError: Invalid run type
    :return: Parameters for the run
    :rtype: CoreParameter
    """
    variables = parameter.variables[0].split(", ")
    run_type = parameter.run_type
    seasons = parameter.seasons

    for season in seasons:
        metrics_dict_test = {}
        metrics_dict_ref = {}
        test_data = utils.dataset.Dataset(parameter, test=True)
        parameter.test_name_yrs = utils.general.get_name_and_yrs(
            parameter, test_data, season
        )
        parameter.ref_name_yrs = "Aerosol Global Benchmarks (Present Day)"

        for aerosol in variables:
            logger.info("Variable: {}".format(aerosol))
            metrics_dict_test[aerosol] = generate_metrics_dic(
                test_data, aerosol, season
            )

        if run_type == "model_vs_model":
            ref_data = utils.dataset.Dataset(parameter, ref=True)
            parameter.ref_name_yrs = utils.general.get_name_and_yrs(
                parameter, ref_data, season
            )
            for aerosol in variables:
                metrics_dict_ref[aerosol] = generate_metrics_dic(
                    ref_data, aerosol, season
                )

        elif run_type == "model_vs_obs":
            parameter.ref_name = parameter.ref_name_yrs
            ref_data_path = os.path.join(
                e3sm_diags.INSTALL_PATH,
                "control_runs",
                "aerosol_global_metrics_benchmarks.json",
            )

            with open(ref_data_path, "r") as myfile:
                ref_file = myfile.read()

            metrics_ref = json.loads(ref_file)

            for aerosol in variables:
                metrics_dict_ref[aerosol] = metrics_ref[aerosol]
        else:
            raise ValueError("Invalid run_type={}".format(run_type))

        parameter.output_file = f"{parameter.test_name}-{season}-budget-table"
        fnm = os.path.join(
            utils.general.get_output_dir(parameter.current_set, parameter),
            parameter.output_file + ".csv",
        )

        with open(fnm, "w") as table_csv:
            writer = csv.writer(
                table_csv,
                delimiter=",",
                quotechar="'",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator="\n",
            )
            writer.writerow([f"Test: {parameter.test_name_yrs}"])
            writer.writerow([f"Ref: {parameter.ref_name_yrs}"])
            writer.writerow(
                [
                    " ",
                    "Test",
                    "Ref",
                ]
            )
            for key, values in metrics_dict_test.items():
                writer.writerow([SPECIES_NAMES[key]])
                for value in values:
                    line = []
                    line.append(value)
                    line.append(values[value])
                    line.append(metrics_dict_ref[key][value])
                    writer.writerows([line])
                writer.writerows([""])

    logger.info(f"Metrics saved in {fnm}")

    return parameter
