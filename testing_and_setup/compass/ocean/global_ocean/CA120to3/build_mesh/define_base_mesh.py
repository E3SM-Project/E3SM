#!/usr/bin/env python
'''
name: define_base_mesh
authors: Phillip J. Wolfram

This function specifies the resolution for a coastal refined mesh for the CA coast from SF to LA for
Chris Jeffrey and Mark Galassi.
It contains the following resolution resgions:
  1) a QU 120km global background resolution
  2) 3km refinement region along the CA coast from SF to LA, with 30km transition region

'''
import numpy as np
import mpas_tools.mesh.creation.coastal_tools as ct


def cellWidthVsLatLon():
    km = 1000.0

    params = ct.default_params

    SFtoLA        = {"include": [np.array([-124.0, -117.5, 34.2, 38.0])],  # SF to LA
                     "exclude": [np.array([-122.1, -120.8, 37.7, 39.2])]}  # SF Bay Delta

    WestCoast = np.array([-136.0, -102.0, 22.0, 51])

    print("****QU120 background mesh and 300m refinement from SF to LA****")
    params["mesh_type"] = "QU"
    params["dx_max_global"] = 120.0 * km
    params["region_box"] = SFtoLA
    params["plot_box"] = WestCoast
    params["dx_min_coastal"] = 3.0 * km
    params["trans_width"] = 100.0 * km
    params["trans_start"] = 30.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    return cell_width / 1000, lon, lat
