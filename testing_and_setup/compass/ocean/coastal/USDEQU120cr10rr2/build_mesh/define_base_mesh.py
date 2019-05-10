#!/usr/bin/env python
'''
name: define_base_mesh
authors: Steven Brus, Phillip J. Wolfram

This function specifies the resolution for a coastal refined mesh for Delaware Bay.
It contains the following resolution resgions:
  1) a QU 120km global background resolution
  2) 10km refinement region from the coast to past the shelf-break from North Carolina to New Hampshire
  3) 5km refinement region from the coast up to the shelf-break from Virginia to past Long Island, New York
  4) 2km refinement region inside Delaware Bay

'''
import numpy as np
import jigsaw_to_MPAS.coastal_tools as ct


def cellWidthVsLatLon():
    km = 1000.0

    params = ct.default_params

    print("****QU120 background mesh and 10km refinement from NC to NH****")
    params["mesh_type"] = "QU"
    params["dx_max_global"] = 120.0 * km
    params["region_box"] = ct.Delaware_Bay
    params["plot_box"] = ct.Western_Atlantic
    params["dx_min_coastal"] = 10.0 * km
    params["trans_width"] = 600.0 * km
    params["trans_start"] = 400.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    print("****5km refinement along coast from VA to NY****")
    params["region_box"] = ct.Delaware_Region
    params["plot_box"] = ct.Delaware
    params["dx_min_coastal"] = 5.0 * km
    params["trans_width"] = 175.0 * km
    params["trans_start"] = 75.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    print("****2km refinement inside Delaware Bay****")
    params["region_box"] = ct.Delaware_Bay
    params["plot_box"] = ct.Delaware
    params["restrict_box"] = ct.Delaware_restrict
    params["dx_min_coastal"] = 2.0 * km
    params["trans_width"] = 100.0 * km
    params["trans_start"] = 17.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    return cell_width / 1000, lon, lat
