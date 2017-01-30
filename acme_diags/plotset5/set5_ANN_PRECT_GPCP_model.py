#!/usr/bin/env python
import numpy
import cdutil
import cdms2
import acme_diags.acme_parser
import plot_set_5
    

parser = acme_diags.acme_parser.ACMEParser()
parameter = parser.get_parameter()

reference_data_path = parameter.reference_data_path
reference_data_set = parameter.reference_data_set # observation
test_data_path = parameter.test_data_path
test_data_set = parameter.test_data_set # model

var = parameter.variables
season = parameter.season

filename1=test_data_path + test_data_set
filename2=reference_data_path + reference_data_set

if var == 'PRECT':

    f_obs = cdms2.open(reference_data_path + reference_data_set)
    f_mod = cdms2.open(test_data_path + test_data_set)
    obs_pr = f_obs(var, longitude=(-180, 540))
    mod_pr = (f_mod('PRECC', longitude=(-180, 540)) + f_mod('PRECL', longitude=(-180, 540)))*3600.0*24.0*1000.0
    mod_pr.units = 'mm/day'
    obs_name = 'GPCP (yrs1979-2009)'

elif var == 'T':

#    plev = 850
    plev = parameter.plev

    for filename in [filename1,filename2]:
        f_in = cdms2.open(filename)
        mv = f_in(var)
        mv_plv = mv.getLevel()

        if mv_plv.long_name.lower().find('hybrid') != -1: # var(time,lev,lon,lat) convert from hybrid level to pressure
            hyam = f_in('hyam')
            hybm = f_in('hybm')
            ps = f_in('PS')/100.    #convert unit from 'Pa' to mb
            p0 = 1000. #mb
            levels_orig = cdutil.vertical.reconstructPressureFromHybrid(ps,hyam,hybm,p0)
            levels_orig.units = 'mb'
            mv_p=cdutil.vertical.logLinearInterpolation(mv, levels_orig, plev)

        elif mv_plv.long_name.lower().find('pressure') != -1 or mv_plv.long_name.lower().find('isobaric') != -1: # levels are presure levels
            #Construct pressure level for interpolation
            levels_orig = mv.clone()
            levels_orig.units = 'mb'
            levels_orig.id = 'pressure'

            for ilev in range(len(mv_plv[:])):
                levels_orig.data[:,ilev,:,:]=mv_plv[ilev]

            mv_p=cdutil.vertical.logLinearInterpolation(mv[:,::-1,:,:], levels_orig[:,::-1,:,:], plev)   #logLinearInterpolation only takes positive down plevel: "I :      interpolation field (usually Pressure or depth) from TOP (level 0) to BOTTOM (last level), i.e P value going up with each level"

        else:
            print( 'Vertical level is neither hybrid nor pressure. Abort')
            quit()
        print filename

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

if var == 'T':
    parameter.main_title = ' '.join([var, str(plev), 'mb', season])

plot_set_5.plot(obs_pr, mod_pr, obs_pr_reg, mod_pr_reg, parameter)
