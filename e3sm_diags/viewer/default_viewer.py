"""
The viewer for the current plotsets supported by
E3SM Diagnostics as of v1.7.0.
"""

import collections
import json
import os
from collections import OrderedDict
from typing import Dict

import numpy

from e3sm_diags.logger import custom_logger
from e3sm_diags.parser import SET_TO_PARSER
from e3sm_diags.viewer.core_viewer import OutputViewer

from . import lat_lon_viewer, utils

logger = custom_logger(__name__)
# A dictionary of the sets to a better name which
# is displayed in the viewer.
# These are all of the sets that this viewer supports.
SET_TO_NAME = {
    "zonal_mean_xy": "Zonal mean line plots",
    "lat_lon": "Latitude-Longitude contour maps",
    "lat_lon_land": "Latitude-Longitude contour maps (land variables)",
    "lat_lon_river": "Latitude-Longitude contour maps (river variables)",
    "polar": "Polar contour maps",
    "cosp_histogram": "CloudTopHeight-Tau joint histograms",
    "diurnal_cycle": "Diurnal cycle phase maps",
    "arm_diags": "Diagnostics at ARM stations",
    "aerosol_aeronet": "Aerosol Diags at AERONET sites",
    "aerosol_budget": "Aerosol Budget Tables",
}

# The ordering of the columns in the viewer.
# These 'seasons' can be months as well, if a user
# wants that:
#     SEASONS = ['01', '02', ..., '12']
SEASONS = [
    "ANN",
    "DJF",
    "MAM",
    "JJA",
    "SON",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
]


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for a certain set of diagnostics,
    create a single page.

    Return the title and url for this page.
    """
    # Information for each row in the viewer.
    # The order is:
    # {
    #  case_id: {
    #       row_name: {
    #             season: filename
    #       }
    #   }
    # }
    ROW_INFO: OrderedDict[
        str, Dict[str, Dict[str, Dict[str, Dict[str, str]]]]
    ] = collections.OrderedDict()
    # A similar dict, but for creating the lat-lon tables.
    LAT_LON_TABLE_INFO: OrderedDict[
        str, Dict[str, Dict[str, Dict[str, str]]]
    ] = collections.OrderedDict()

    # Since we're only having one set for each
    # create_viewer() call, this works.
    set_name = parameters[0].sets[0]

    viewer = OutputViewer(path=root_dir)
    cols = ["Description"] + seasons_used(parameters)
    viewer.add_page(SET_TO_NAME[set_name], short_name=set_name, columns=cols)

    for parameter in parameters:
        results_dir = parameter.results_dir

        # Filenames are like one in the two formats:
        #   ref_name-variable-season-region
        #   ref_name-variable-plev'mb'-season-region
        ref_name = getattr(parameter, "ref_name", "")
        for var in parameter.variables:
            for season in parameter.seasons:
                for region in parameter.regions:
                    # Since some parameters have plevs, there might be
                    # more than one row_name, filename pair.
                    # Hence, this is why we have a list.
                    row_name_and_filename = []

                    if parameter.plevs == []:  # 2d variables.
                        if parameter.run_type == "model_vs_model":
                            row_name = "{} {}".format(var, region)
                        else:
                            row_name = "{} {} {}".format(var, region, ref_name)
                        fnm = "{}-{}-{}-{}".format(ref_name, var, season, region)
                        row_name_and_filename.append((row_name, fnm))
                    else:  # 3d variables.
                        for plev in parameter.plevs:
                            if parameter.run_type == "model_vs_model":
                                row_name = "{}-{} {}".format(
                                    var, str(int(plev)) + "mb", region
                                )
                            else:
                                row_name = "{}-{} {} {}".format(
                                    var,
                                    str(int(plev)) + "mb",
                                    region,
                                    ref_name,
                                )
                            fnm = "{}-{}-{}-{}-{}".format(
                                ref_name, var, int(plev), season, region
                            )
                            row_name_and_filename.append((row_name, fnm))

                    if set_name in ["lat_lon", "lat_lon_land", "lat_lon_river"]:
                        metrics_path = os.path.join(
                            results_dir,
                            "{}".format(set_name),
                            parameter.case_id,
                            fnm,
                        )
                        if os.path.exists(metrics_path + ".json"):
                            _add_to_lat_lon_metrics_table(
                                LAT_LON_TABLE_INFO,
                                metrics_path,
                                season,
                                row_name,
                            )
                        else:
                            logger.warning(
                                (
                                    "JSON does not exist: {}".format(
                                        metrics_path + ".json"
                                    )
                                )
                            )
                            continue
                    for row_name, fnm in row_name_and_filename:
                        if set_name not in ROW_INFO:
                            ROW_INFO[set_name] = collections.OrderedDict()
                        if parameter.case_id not in ROW_INFO[set_name]:
                            ROW_INFO[set_name][
                                parameter.case_id
                            ] = collections.OrderedDict()
                        if row_name not in ROW_INFO[set_name][parameter.case_id]:
                            ROW_INFO[set_name][parameter.case_id][
                                row_name
                            ] = collections.OrderedDict()
                            ROW_INFO[set_name][parameter.case_id][row_name][
                                "descr"
                            ] = _get_description(var, parameter)
                        # Each season has a image_path and metadata linked to it, thus we use a dict.
                        ROW_INFO[set_name][parameter.case_id][row_name][season] = {}
                        # Format the filename to support relative paths.
                        ROW_INFO[set_name][parameter.case_id][row_name][season][
                            "image_path"
                        ] = os.path.join(
                            "..", "{}".format(set_name), parameter.case_id, fnm
                        )
                        ROW_INFO[set_name][parameter.case_id][row_name][season][
                            "metadata"
                        ] = create_metadata(parameter)

    save_netcdf = parameters[0].save_netcdf
    ext = parameters[0].output_format[0]
    _add_information_to_viewer(viewer, ROW_INFO, set_name, save_netcdf, ext)

    viewer.set_page(SET_TO_NAME[set_name])
    # Return the name and the url of the current page of the viewer.
    name, url = SET_TO_NAME[set_name], viewer.generate_page()

    # Edit the final URL to add in the header.
    utils.add_header(root_dir, os.path.join(root_dir, url), parameters)
    utils.h1_to_h3(os.path.join(root_dir, url))

    if set_name == "lat_lon":
        table_tuple = lat_lon_viewer.generate_lat_lon_metrics_table(
            LAT_LON_TABLE_INFO, SEASONS, viewer, root_dir, parameters
        )
        taylor_diag_tuple = lat_lon_viewer.generate_lat_lon_taylor_diag(
            LAT_LON_TABLE_INFO, SEASONS, viewer, root_dir, parameters
        )
        if parameter.run_type == "model_vs_obs":
            cmip6_comparison_tuple = lat_lon_viewer.generate_lat_lon_cmip6_comparison(
                LAT_LON_TABLE_INFO, SEASONS, viewer, root_dir, parameters
            )
            return [(name, url), table_tuple, taylor_diag_tuple, cmip6_comparison_tuple]

        print((name, url), table_tuple, taylor_diag_tuple)
        return [(name, url), table_tuple, taylor_diag_tuple]

    if set_name == "lat_lon_land" or set_name == "lat_lon_river":
        table_tuple = lat_lon_viewer.generate_lat_lon_metrics_table(
            LAT_LON_TABLE_INFO, SEASONS, viewer, root_dir, parameters
        )
        return [(name, url), table_tuple]

    return (name, url)


def seasons_used(parameters):
    """
    Get a list of the seasons used for this set of parameters.
    """
    seasons_used = set([s for p in parameters for s in p.seasons])
    # Make sure this list is ordered by SEASONS.
    return [season for season in SEASONS if season in seasons_used]


def _get_description(var, parameters):
    """
    Get the description for the variable from the parameters.
    This was manually added in the driver to the viewer_descr parameter.
    """
    if hasattr(parameters, "viewer_descr") and var in parameters.viewer_descr:
        return parameters.viewer_descr[var]
    return var


def _add_to_lat_lon_metrics_table(lat_lon_table_info, metrics_path, season, row_name):
    """
    Add the metrics for the current season and
    row_name to the lat-lon table.
    """
    with open(metrics_path + ".json") as json_file:
        metrics_dict = json.load(json_file)

        if season not in lat_lon_table_info:
            lat_lon_table_info[season] = collections.OrderedDict()
        if row_name not in lat_lon_table_info[season]:
            lat_lon_table_info[season][row_name] = collections.OrderedDict()
        lat_lon_table_info[season][row_name]["metrics"] = metrics_dict


def create_metadata(parameter):
    """
    From a set of parameters, extract the metadata.
    """
    metadata = collections.OrderedDict()
    msg = "Use this command to recreate this image:"
    metadata[msg] = ""

    set_name = parameter.sets[0]
    parser = SET_TO_PARSER[set_name]()
    cmd = "e3sm_diags {} --no_viewer ".format(set_name)

    args = parser.view_args()
    supported_cmd_args = list(args.__dict__.keys())

    if "other_parameters" in supported_cmd_args:
        supported_cmd_args.remove("other_parameters")

    if "parameters" in supported_cmd_args:
        supported_cmd_args.remove("parameters")

    for param_name in parameter.__dict__:
        param = parameter.__dict__[param_name]
        # We don't want to include blank values.
        if (isinstance(param, numpy.ndarray) and not param.all()) or (
            not isinstance(param, numpy.ndarray) and not param
        ):
            continue

        if param_name in supported_cmd_args:
            if (
                isinstance(param, list)
                or isinstance(param, tuple)
                or isinstance(param, numpy.ndarray)
            ):
                # ex: --diff_levels -7, -6, -5, -4
                cmd += "--{} ".format(param_name)
                for p in param:
                    if isinstance(p, str) and p.isdigit():
                        cmd += " {} ".format(str(p))
                    else:
                        cmd += " '{}' ".format(str(p))

            elif isinstance(param, bool):
                # ex: --multiprocessing
                # note there's no value after the parameter, it's just a flag
                if param:  # command is True, so add --command to set it to True
                    cmd += "--{} ".format(param_name)

            elif isinstance(param, str) and param.isdigit():
                cmd += "--{} {} ".format(param_name, param)
            else:
                cmd += "--{} '{}' ".format(param_name, param)

    metadata[msg] = cmd

    return metadata


def _add_information_to_viewer(viewer, row_info, set_name, save_netcdf, ext="png"):
    """
    Based on the information in row_info,
    add the rows/cols to the viewer object.
    """
    viewer.set_page(SET_TO_NAME[set_name])

    for group in row_info[set_name]:
        try:
            viewer.set_group(group)
        except RuntimeError:
            viewer.add_group(group)

        for row_name in row_info[set_name][group]:
            try:
                viewer.set_row(row_name)
            except RuntimeError:
                # A row of row_name wasn't in the viewer, so add it.
                viewer.add_row(row_name)
                # Adding the description for the row.
                viewer.add_col(row_info[set_name][group][row_name]["descr"])

            # [1:] is to ignore 'Description' col.
            for col_season in viewer.page.columns[1:]:
                if col_season not in row_info[set_name][group][row_name]:
                    viewer.add_col("-----", is_file=True, title="-----")
                else:
                    metadata = row_info[set_name][group][row_name][col_season][
                        "metadata"
                    ]
                    fnm = row_info[set_name][group][row_name][col_season]["image_path"]
                    formatted_files = []
                    if save_netcdf:
                        nc_files = [
                            fnm + nc_ext
                            for nc_ext in ["_test.nc", "_ref.nc", "_diff.nc"]
                        ]
                        formatted_files = [{"url": f, "title": f} for f in nc_files]
                    viewer.add_col(
                        fnm + "." + ext,
                        is_file=True,
                        title=col_season,
                        other_files=formatted_files,
                        meta=metadata,
                    )
