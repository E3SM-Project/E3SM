#!/usr/bin/env python
import mpas_tools.ocean.coastal_tools as ct
from mpas_tools.ocean import build_spherical_mesh


def cellWidthVsLatLon():
    """
    This function specifies the resolution for a coastal refined mesh for
    Delaware Bay.  It creates cell width array for this mesh on a regular
    latitude-longitude grid.

    It contains the following resolution regions:
      1) a QU 120km global background resolution
      2) 10km refinement region from the coast to past the shelf-break from
         North Carolina to New Hampshire
      3) 5km refinement region from the coast up to the shelf-break from
         Virginia to past Long Island, New York
      4) 2km refinement region inside Delaware Bay

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
    # authors: Steven Brus, Phillip J. Wolfram
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


def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc',
                         do_inject_bathymetry=True, preserve_floodplain=True,
                         floodplain_elevation=20.0)


if __name__ == '__main__':
    main()
