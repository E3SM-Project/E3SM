#!/usr/bin/env python
import mpas_tools.ocean.coastal_tools as ct
from mpas_tools.ocean import build_spherical_mesh


def cellWidthVsLatLon():
    """
    This function creates cell width array for this mesh on a regular
    latitude-longitude grid.

    Returns
    -------
       cellWidth : ndarray
            m x n array, entries are desired cell width in km

       lat : ndarray
            latitude, vector of length m, with entries between -90 and 90,
            degrees

       lon : ndarray
            longitude, vector of length n, with entries between -180 and 180,
            degrees
    """
    km = 1000.0

    params = ct.default_params

    print("****QU 120 background mesh and enhanced Atlantic (30km)****")
    params["mesh_type"] = "QU"
    params["dx_max_global"] = 120.0 * km
    params["region_box"] = ct.Atlantic
    params["restrict_box"] = ct.Atlantic_restrict
    params["plot_box"] = ct.Western_Atlantic
    params["dx_min_coastal"] = 30.0 * km
    params["trans_width"] = 5000.0 * km
    params["trans_start"] = 500.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    print("****Northeast refinement (10km)***")
    params["region_box"] = ct.Delaware_Bay
    params["plot_box"] = ct.Western_Atlantic
    params["dx_min_coastal"] = 10.0 * km
    params["trans_width"] = 600.0 * km
    params["trans_start"] = 400.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    print("****Delaware regional refinement (5km)****")
    params["region_box"] = ct.Delaware_Region
    params["plot_box"] = ct.Delaware
    params["dx_min_coastal"] = 5.0 * km
    params["trans_width"] = 175.0 * km
    params["trans_start"] = 75.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    print("****Delaware Bay high-resolution (2km)****")
    params["region_box"] = ct.Delaware_Bay
    params["plot_box"] = ct.Delaware
    params["restrict_box"] = ct.Delaware_restrict
    params["dx_min_coastal"] = 2.0 * km
    params["trans_width"] = 100.0 * km
    params["trans_start"] = 17.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    return cell_width / 1000, lon, lat


def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


if __name__ == '__main__':
    main()
