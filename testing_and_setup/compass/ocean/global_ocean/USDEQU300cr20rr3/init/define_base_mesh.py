#!/usr/bin/env python

import numpy as np
import jigsaw_to_MPAS.coastal_tools as ct


def cellWidthVsLatLon():
    km = 1000.0

    params = ct.default_params

    params["mesh_type"] = "QU"
    params["dx_max_global"] = 300.0 * km
    ct.Delaware_Region["include"][0] = np.array([-77, -69.8, 37.25, 40.25])
    params["region_box"] = ct.Delaware_Region
    params["plot_box"] = ct.Western_Atlantic
    params["dx_min_coastal"] = 20.0 * km
    params["trans_width"] = 600.0 * km
    params["trans_start"] = 100.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    ct.Delaware_Bay["include"][0] = np.array(
        [-75.61903, -74.7, 38.5, 40.312747])
    params["region_box"] = ct.Delaware_Bay
    params["plot_box"] = ct.Delaware
    params["restrict_box"] = ct.Delaware_restrict
    params["dx_min_coastal"] = 3.0 * km
    params["trans_width"] = 100.0 * km
    params["trans_start"] = 17.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    return cell_width / 1000, lon, lat
