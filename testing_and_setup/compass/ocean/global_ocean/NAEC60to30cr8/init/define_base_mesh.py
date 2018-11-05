#!/usr/bin/env python

import numpy as np
import coastal_tools as ct


def cellWidthVsLatLon():

    km = 1000.0

    params = ct.default_params
    params["dx_min_coastal"] = 8.0*km
    params["trans_width"] = 600.0*km
    params["trans_start"] = 400.0*km
    params["mesh_type"] = "EC"

    params["plot_box"] = np.array([-80,-30,30,90]) 

    print "***Gulf Coast***"
    params["region_box"] = ct.US_Gulf_Coast
    params["restrict_box"] = ct.Gulf_restrict
    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    print "***East Coast***"
    params["region_box"] = ct.US_East_Coast
    params["restrict_box"] = ct.East_Coast_restrict
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Caribbean***"
    params["region_box"] = ct.Caribbean
    params["restrict_box"] = ct.Caribbean_restrict
    params["trans_width"] = 400.0*km
    params["trans_start"] = 300.0*km
    params["n_longest"] = 50
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***West Coast***"
    params["region_box"] = ct.US_West_Coast
    params["restrict_box"] = ct.Empty
    params["trans_width"] = 600.0*km
    params["trans_start"] = 400.0*km
    params["n_longest"] = 10
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Alaska***"
    params["region_box"] = ct.Alaska
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Aleutian Islands (West)***"
    params["region_box"] = ct.Aleutian_Islands_W
    params["lon_wrap"]  = True
    params["n_longest"] = 100
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "****Aleutian Islands (East)***"
    params["region_box"] = ct.Aleutian_Islands_E
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Bering Sea (East)****"
    params["region_box"] = ct.Bering_Sea_E
    params["trans_start"] = 400.0*km
    params["trans_width"] = 600.0*km
    params["n_longest"] = 10
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Bering Sea (West)***"
    params["region_box"] = ct.Bering_Sea_W
    params["restrict_box"] = ct.Bering_Sea_restrict
    params["trans_start"] = 450.0*km
    params["trans_width"] = 600.0*km
    params["n_longest"] = 10
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Hawaii***"
    params["region_box"] = ct.Hawaii
    params["restrict_box"] = ct.Empty
    params["trans_start"] = 100*km
    params["trans_width"] = 200*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)
 
    print "***Newfoundland***" 
    Newfoundland = {"include":[np.array([-65.0,-50.0,44.0,60.0]),
                               np.array([-65.5,-64.5,61.0,62.0])],
                    "exclude":[]}
    params["region_box"] = Newfoundland 
    params["restrict_box"] = ct.Empty
    params["trans_width"] = 600.0*km
    params["trans_start"] = 400.0*km
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Greenland***"
    Greenland = {"include":[np.array([-49,-45,60,61])],
                    "exclude":[]}
    params["region_box"] = Greenland 
    params["restrict_box"] = ct.Empty
    params["trans_start"] = 900.0*km
    cell_width, lon, lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)
    
    return cell_width / km, lon, lat
