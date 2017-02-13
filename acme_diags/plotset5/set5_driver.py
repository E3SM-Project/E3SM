#!/usr/bin/env python
import json
import copy
import numpy
import cdutil
import genutil
import cdms2
import MV2
import acme_diags.acme_parser
import plot_set_5
import glob
import os
import fnmatch
import re


def make_parameters(orginal_parameter):
    f_data = open('set5_diags.json').read()
    json_file = json.loads(f_data)

    parameters = []
    for key in json_file:
        for single_run in json_file[key]:
            p = copy.deepcopy(orginal_parameter)
            for attr_name in single_run:
                setattr(p, attr_name, single_run[attr_name])
            parameters.append(p)
    return parameters

parser = acme_diags.acme_parser.ACMEParser()
orginal_parameter = parser.get_parameter()
parameters = make_parameters(orginal_parameter)


for parameter in parameters:

    reference_data_path = parameter.reference_data_path
#    reference_data_set = parameter.reference_data_set # observation

    test_data_path = parameter.test_data_path
#    test_data_set = parameter.test_data_set # model

    var = parameter.variables
    season = parameter.season
    ref_name = parameter.ref_name
    test_name = parameter.test_name


    for it in season:
        test_files = glob.glob(os.path.join(test_data_path,'*'+test_name+'*.nc'))
        for filename in fnmatch.filter(test_files, '*'+it+'*'):
            print filename
            filename1 = filename
 
        ref_files = glob.glob(os.path.join(reference_data_path,'*'+ref_name+'*.nc'))
        for filename in fnmatch.filter(ref_files, '*'+it+'*'):
            print filename
            filename2 = filename
    
#        filename1=test_data_path + test_data_set
#        filename2=reference_data_path + reference_data_set
        print filename1, filename2
        if var == 'PRECT':
            f_obs = cdms2.open(filename2)
            f_mod = cdms2.open(filename1)
            if ref_name == 'TRMM': #TRMM region following AMWG  latitude=(-38,38)
                obs_pr = f_obs(var, latitude=(-38,38))
                mod_pr = (f_mod('PRECC', latitude=(-38,38)) + f_mod('PRECL', latitude=(-38,38)))*3600.0*24.0*1000.0
            else:
                obs_pr = f_obs(var)
                mod_pr = (f_mod('PRECC') + f_mod('PRECL'))*3600.0*24.0*1000.0
            mod_pr.units = 'mm/day'
            parameter.output_file = '_'.join([ref_name,it])
            print parameter.output_file
    
        elif var == 'T':
    
            plev = parameter.levels[0]
            print 'selected pressure level', plev
    
            for filename in [filename1,filename2]:
                f_in = cdms2.open(filename)
                mv = f_in[var] # Square brackets for metadata preview
                mv_plv = mv.getLevel()
    
                if mv_plv.long_name.lower().find('hybrid') != -1: # var(time,lev,lon,lat) convert from hybrid level to pressure
                    mv = f_in (var) # Parentheses actually load the transient variable
                    hyam = f_in('hyam')
                    hybm = f_in('hybm')
                    ps = f_in('PS')/100.    #convert unit from 'Pa' to mb
                    p0 = 1000. #mb
                    levels_orig = cdutil.vertical.reconstructPressureFromHybrid(ps,hyam,hybm,p0)
                    levels_orig.units = 'mb'
                    mv_p=cdutil.vertical.logLinearInterpolation(mv, levels_orig, plev)
    
                elif mv_plv.long_name.lower().find('pressure') != -1 or mv_plv.long_name.lower().find('isobaric') != -1: # levels are presure levels
                    # if desired plev exists, read from input file, otherwise, interpolate
                    try:
                        plev_ind = mv_plv[:].tolist().index(plev)
                        mv_p = f_in(var,level= plev)
    
                    except Exception as e:
                        print str(plev)+ ' is not standard pressure level, execute vertical interpolation'
                        mv = f_in (var)
                        #Construct pressure level for interpolation
                        levels_orig = MV2.array(mv_plv[:])
                        levels_orig.setAxis(0,mv_plv)
                        mv,levels_orig = genutil.grower(mv,levels_orig) # grow 1d levels_orig to mv dimention
                        #levels_orig.info()
    
                        mv_p=cdutil.vertical.logLinearInterpolation(mv[:,::-1,], levels_orig[:,::-1,], plev)   #logLinearInterpolation only takes positive down plevel: "I :      interpolation field (usually Pressure or depth) from TOP (level 0) to BOTTOM (last level), i.e P value going up with each level"
    
                else:
                    print( 'Vertical level is neither hybrid nor pressure. Abort')
                    quit()
                print filename
                parameter.output_file = '_'.join([ref_name,it, str(plev)+'mb'])
                print parameter.output_file
    
                if filename == filename1:
                    mod_pr = mv_p
    
                if filename == filename2:
                    obs_pr = mv_p
    
        axes1 = mod_pr.getAxisList()
        axes2 = obs_pr.getAxisList()
    
    
        # For plotting, original grid is plotted for model observation, differece plot is regridded to coaser grid. Need if statement to evaluate grid size. aminusb_2ax from uvcmetrics takes care of this,which also considers complex corner cases.
        if len(axes1[1]) <= len(axes2[1]): # use nlat to decide data resolution, higher number means higher data resolution. For the difference plot, regrid toward lower resolution
            model_grid = mod_pr.getGrid()
            mod_pr_reg = mod_pr
            obs_pr_reg = obs_pr.regrid(model_grid, regridTool=parameter.regrid_tool, regridMethod=parameter.regrid_method)
        else:
            obs_grid = obs_pr.getGrid()
            obs_pr_reg = obs_pr
            mod_pr_reg = mod_pr.regrid(obs_grid, regridTool=parameter.regrid_tool, regridMethod=parameter.regrid_method)
    
        #if var == 'T' and plev == 850:
        #    levels = [230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300]
        #    diff_levels = [-8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 8]
        #    parameter.reference_levels = levels
        #    parameter.test_levels = levels
        #    parameter.diff_levels = diff_levels
        #elif var == 'T' and plev == 200:
        #    levels = [190, 193, 196, 199, 202, 205, 208, 211, 214, 217, 220, 223, 226, 229, 232]
        #    diff_levels = [-10, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10]
        #    parameter.reference_levels = levels
        #    parameter.test_levels = levels
        #    parameter.diff_levels = diff_levels
    
    ###not working
    #    if var == 'T':
    #        #parameter.main_title = ' '.join([var, str(plev), 'mb', season])
    #        parameter.main_title = ' '.join([var, season])
    #    else:
    #        parameter.main_title = ' '.join([var, season])
    
        plot_set_5.plot(obs_pr, mod_pr, obs_pr_reg, mod_pr_reg, parameter)
