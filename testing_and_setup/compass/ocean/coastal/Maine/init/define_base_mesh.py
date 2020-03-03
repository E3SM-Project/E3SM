#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import jigsaw_to_MPAS.coastal_tools as ct


def cellWidthVsLatLon():

    km = 1000.0

    params = ct.default_params

    # Use 1.0 degree for fast scoping, and 0.1 for the final creation:
    params["ddeg"] = 0.2
    #params["ddeg"] = 1.0
    #params["ddeg"] = 0.5

    params["dx_min_coastal"] = 12.0*km
    params["trans_start"] = 400.0*km
    params["trans_width"] = 600.0*km
    params["n_longest"] = 20
    params["mesh_type"] = "EC"

    params["plot_box"] = ct.Entire_Globe
    params["plot_option"] = False

    print("***Gulf Coast***")
    params["region_box"] = ct.US_Gulf_Coast
    params["restrict_box"] = ct.Gulf_restrict
    params["trans_start"] = 400.0*km
    params["trans_width"] = 600.0*km
    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    print("***East Coast***")
    params["region_box"] = ct.US_East_Coast
    params["restrict_box"] = ct.East_Coast_restrict
    params["trans_start"] = 400.0*km
    params["trans_width"] = 600.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Gulf Stream extension***")
    params["restrict_box"] = ct.Empty
    params["trans_width"] = 600.0*km
    params["region_box"] = {"include":[np.array([-78,-70,40,42])],
                            "exclude":[]}
    params["trans_start"] = 700.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)
    params["region_box"] = {"include":[np.array([-74,-50,42,55])],
                            "exclude":[]}
    params["trans_start"] = 1100.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Gulf Coast spot***")
    params["restrict_box"] = ct.Empty
    params["trans_width"] = 100.0*km
    params["region_box"] = {"include":[np.array([-98,-94,28,30])],
                            "exclude":[]}
    params["trans_start"] = 400.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)
    
    print("***Greenland***")
    params["region_box"] = ct.Greenland
    params["restrict_box"] = ct.Empty
    params["trans_width"] = 600.0*km
    params["trans_start"] = 275.0*km
    cell_width, lon, lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Greenland-Iceland smooth kink***")
    params["region_box"] = {"include":[np.array([-45,-40,59,62])],
                            "exclude":[]}
    params["trans_start"] = 700.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print( "***Hudson Bay***")
    params["region_box"] = {"include":[np.array([-100,-70,50,65])],
                            "exclude":[]}
    params["restrict_box"] = {"include":[np.array([-100,-70,50,65])],
                            "exclude":[]}
    params["trans_start"] = 600.0*km
    params["trans_width"] = 100.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print( "***Northern Canada***")
    params["region_box"] = {"include":[np.array([-168,-59,67,85])],
                            "exclude":[]}
    params["restrict_box"] = ct.Empty
    params["trans_start"] = 275.0*km
    params["trans_width"] = 600.0*km
    params["n_longest"] = 20
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Alaska***")
    params["region_box"] = ct.Alaska
    params["trans_start"] = 400.0*km
    params["trans_width"] = 600.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print( "***Crook at Alaska/Canada***")
    params["region_box"] = {"include":[np.array([-144,-127,68.5,71])],
                            "exclude":[]}
    params["restrict_box"] = ct.Empty
    params["trans_start"] = 550.0*km
    params["trans_width"] = 600.0*km
    params["n_longest"] = 20
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Caribbean***")
    params["region_box"] = ct.Caribbean
    params["restrict_box"] = ct.Caribbean_restrict
    params["trans_width"] = 400.0*km
    params["trans_start"] = 300.0*km
    params["n_longest"] = 50
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***West Coast***")
    params["region_box"] = ct.US_West_Coast
    params["restrict_box"] = ct.Empty
    params["trans_width"] = 600.0*km
    params["trans_start"] = 400.0*km
    params["n_longest"] = 10
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Aleutian Islands (West)***")
    params["region_box"] = ct.Aleutian_Islands_W
    params["n_longest"] = 100
    params["trans_start"] = 200.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("****Aleutian Islands (East)***")
    params["region_box"] = ct.Aleutian_Islands_E
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Bering Sea (East)****")
    params["region_box"] = ct.Bering_Sea_E
    params["trans_start"] = 400.0*km
    params["trans_width"] = 600.0*km
    params["n_longest"] = 10
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Bering Sea (West)***")
    params["region_box"] = ct.Bering_Sea_W
    params["restrict_box"] = ct.Bering_Sea_restrict
    params["trans_start"] = 450.0*km
    params["trans_width"] = 600.0*km
    params["n_longest"] = 10
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Newfoundland***")
    Newfoundland = {"include":[np.array([-65.0,-50.0,44.0,60.0]),
                               np.array([-65.5,-64.5,61.0,62.0])],
                    "exclude":[]}
    params["region_box"] = Newfoundland
    params["restrict_box"] = ct.Empty
    params["trans_width"] = 600.0*km
    params["trans_start"] = 400.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Labrador Sea***")
    params["region_box"] = ct.Empty
    params["restrict_box"] = ct.Empty
    params["point_list"] = [np.array([-50.0,55.0])]
    params["trans_width"] = 600.0*km
    params["trans_start"] = 400.0*km
    cell_width, lon, lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print("***Central America (West Coast)***")
    Central_America = {"include":[np.array([[-110.26,20.69],
                                            [-87.84, 8.94 ],
                                            [-84.55, 12.03],
                                            [-104.26,23.11]]),
                                  np.array([[-88.02, 10.47],
                                            [-81.53, 6.14],
                                            [-81.45, 8.07],
                                            [-84.80, 11.51]]),
                                  np.array([[-81.92, 7.76],
                                            [-76.84, 4.51],
                                            [-77.41, 8.22],
                                            [-79.23, 9.28]])],
                       "exclude":[]}
    params["region_box"] = Central_America
    params["restrict_box"] = ct.Empty
    params["point_list"] = None
    params["trans_width"] = 600.0*km
    params["trans_start"] = 400.0*km
    cell_width, lon, lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    # coast of Maine
    params["dx_min_coastal"] = 6.0*km
    params["trans_start"] = 500.0*km
    params["trans_width"] = 300.0*km
    params["n_longest"] = 20
    params["region_box"] = {"include":[np.array([-70,-68,43,45])],
                            "exclude":[]}
    params["restrict_box"] = {"include":[],
                            "exclude":[np.array([-73,-60,46,50])]}
    params["point_list"] = [[-68.8,43.8]]
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    # coast of Maine
    params["dx_min_coastal"] = 2.0*km
    params["trans_start"] = 150.0*km
    params["trans_width"] = 300.0*km
    params["n_longest"] = 20
    params["region_box"] = {"include":[np.array([-70,-68,43,45])],
                            "exclude":[]}
    params["restrict_box"] = {"include":[],
                            "exclude":[np.array([-73,-50,46,50])]}
    params["point_list"] = [[-68.8,43.8]]
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    return cell_width / km, lon, lat
