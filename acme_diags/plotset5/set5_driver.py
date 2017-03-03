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
from acme_diags.derivations import acme
from acme_diags.derivations.default_regions import regions_specs


def make_parameters(orginal_parameter):
    #f_data = open('set5_diags_default.json').read()
    #f_data = open('set5_diags.json').read()
#    f_data = open('set5_diags_MERRA.json').read()
    #f_data = open('set5_diags_HADISST.json').read()
    #f_data = open('set5_diags_CRU.json').read()
    #f_data = open('set5_diags_LEGATES.json').read()
    f_data = open('set5_diags_WILLMOTT.json').read()
    #f_data = open('set5_diags_NVA.json').read()
    json_file = json.loads(f_data)

    parameters = []
    for key in json_file:
        for single_run in json_file[key]:
            p = copy.deepcopy(orginal_parameter)
            for attr_name in single_run:
                setattr(p, attr_name, single_run[attr_name])
            parameters.append(p)
    return parameters


def regrid_to_lower_res(mv1, mv2, regrid_tool, regrid_method):
    """regrid transient variable toward lower resolution of two variables"""

    axes1 = mv1.getAxisList()
    axes2 = mv2.getAxisList()

    if len(axes1[1]) <= len(axes2[1]): # use nlat to decide data resolution, higher number means higher data resolution. For the difference plot, regrid toward lower resolution
        mv_grid = mv1.getGrid()
        mv1_reg = mv1
        mv2_reg = mv2.regrid(mv_grid, regridTool=regrid_tool, regridMethod=regrid_method)
    else:
        mv_grid = mv2.getGrid()
        mv2_reg = mv2
        mv1_reg = mv1.regrid(mv_grid, regridTool=regrid_tool, regridMethod=regrid_method)
    return mv1_reg, mv2_reg



parser = acme_diags.acme_parser.ACMEParser()
orginal_parameter = parser.get_parameter()
parameters = make_parameters(orginal_parameter)

for parameter in parameters:

    reference_data_path = parameter.reference_data_path
    test_data_path = parameter.test_data_path

    var = parameter.variables
    seasons = parameter.season
    ref_name = parameter.ref_name
    test_name = parameter.test_name
    regions = parameter.region
#    domain = regions_specs[region]

    for season in seasons:
        test_files = glob.glob(os.path.join(test_data_path,'*'+test_name+'*.nc'))
        for filename in fnmatch.filter(test_files, '*'+season+'*'):
            print filename
            filename1 = filename

        ref_files = glob.glob(os.path.join(reference_data_path,'*'+ref_name+'*.nc'))
        for filename in fnmatch.filter(ref_files, '*'+season+'*'):
            print filename
            filename2 = filename

        #print filename1, filename2
        f_mod = cdms2.open(filename1)
        f_obs = cdms2.open(filename2)

        # domain can pass in process_derived_var after Charles fix cdutil.domain'unit problem
        #mv1 = acme.process_derived_var(var, acme.derived_variables, f_mod, domain)
        #mv2 = acme.process_derived_var(var, acme.derived_variables, f_obs, domain)
        mv1 = acme.process_derived_var(var, acme.derived_variables, f_mod)
        mv2 = acme.process_derived_var(var, acme.derived_variables, f_obs)
        #try:
        #    mv1 = f_mod(var)
        #except:
        #    mv1 = acme.process_derived_var(var, acme.derived_variables, f_mod)
        #try:
        #    mv2 = f_obs(var)
        #except:
        #    mv2 = acme.process_derived_var(var, acme.derived_variables, f_obs)
        # special case, cdms didn't properly convert mask with fill value -999.0, filed issue with denise
        if ref_name == 'WILLMOTT':
            #mv2=MV2.masked_where(mv2==mv2.fill_value,mv2)
            mv2=MV2.masked_where(mv2==-999.,mv2) # this is cdms2 for bad mask, denise's fix should fix
            #following should move to derived variable
            days_season = {'ANN':365,'DJF':90,'MAM':92,'JJA':92,'SON':91}
            mv1 = mv1 * days_season[season] * 0.1 #following AMWG approximate way to convert to seasonal cumulative precipitation, need to have solution in derived variable, unit convert from mm/day to cm
            mv1.units = 'cm'
            mv2.units = 'cm'
            
            

         
        #for variables without z axis:
        if mv1.getLevel()==None and mv2.getLevel()==None:
            if len(regions) == 0:
                regions = ['global']
    
            for region in regions: 
                print region
                domain = None
                if region != 'global':
                    domain = regions_specs[region]['domain']
                    # below 7 lines are temporary solution for the cdutil error
                print mv1.units
                mv1_domain = mv1(domain)
                mv2_domain = mv2(domain)
                mv1_domain.units = mv1.units
                mv2_domain.units = mv1.units
                parameter.output_file = '_'.join([ref_name,var, season,region])
                parameter.main_title = str(' '.join([var, season]))
        
                #else:
                #    
                #    parameter.output_file = '_'.join([ref_name, season])
                #    parameter.main_title = str(' '.join([var, season]))

                #parameter.output_file = '_'.join([ref_name, season,region])
                #parameter.main_title = str(' '.join([var, season]))
        
    
                #regrid towards lower resolution of two variables for calculating difference
                mv1_reg, mv2_reg = regrid_to_lower_res(mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                #if var is 'SST' or var is 'TREFHT_LAND': #special case
                
                if var == 'TREFHT_LAND'or var == 'SST': #use "==" instead of "is"
                    if ref_name == 'WILLMOTT':
                       mv2_reg=MV2.masked_where(mv2_reg==mv2_reg.fill_value,mv2_reg)
                       print ref_name
                        
                        #if mv.mask is False: 
                        #    mv = MV2.masked_less_equal(mv, mv._FillValue)
                        #    print "*************",mv.count()
                    land_mask = MV2.logical_or(mv1_reg.mask, mv2_reg.mask)
                    mv1_reg = MV2.masked_where(land_mask, mv1_reg)
                    mv2_reg = MV2.masked_where(land_mask, mv2_reg)
           
                plot_set_5.plot(mv2_domain, mv1_domain, mv2_reg, mv1_reg, parameter)
                #plot_set_5.plot(mv2_reg, mv1_reg, mv2_reg, mv1_reg, parameter)
 
    
        #elif mv1.rank() == 4 and mv2.rank() == 4: #for variables with z axis:
        elif mv1.getLevel() and mv2.getLevel(): #for variables with z axis:
            plev = parameter.levels
            print 'selected pressure level', plev

            for filename in [filename1, filename2]:
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
                    mv = f_in (var)
                    #Construct pressure level for interpolation
                    levels_orig = MV2.array(mv_plv[:])
                    levels_orig.setAxis(0,mv_plv)
                    mv,levels_orig = genutil.grower(mv,levels_orig) # grow 1d levels_orig to mv dimention
                    #levels_orig.info()
                    mv_p=cdutil.vertical.logLinearInterpolation(mv[:,::-1,], levels_orig[:,::-1,], plev)   #logLinearInterpolation only takes positive down plevel: "I :      interpolation field (usually Pressure or depth) from TOP (level 0) to BOTTOM (last level), i.e P value going up with each level"

                else:
                    raise RuntimeError("Vertical level is neither hybrid nor pressure. Abort")

                if filename == filename1:
                    mv1_p = mv_p

                if filename == filename2:
                    mv2_p = mv_p

            for ilev in range(len(plev)):
                mv1 = mv1_p[:,ilev,]
                mv2 = mv2_p[:,ilev,]

                if len(regions) == 0:
                    regions = ['global']
    
                for region in regions: 
                    print region
                    domain = None
                    if region != 'global':
                        domain = regions_specs[region]['domain']
                        # below 7 lines are temporary solution for the cdutil error
                    mv1_domain = mv1(domain)
                    mv2_domain = mv2(domain)
                    mv1_domain.units = mv1.units
                    mv2_domain.units = mv2.units

                    parameter.output_file = '_'.join([ref_name,var,str(int(plev[ilev])),season,region,var, str(int(plev[ilev])), 'mb'])
                    parameter.main_title = str(' '.join([var, str(int(plev[ilev])), 'mb', season]))

                    # Regrid towards lower resolution of two variables for calculating difference
                    mv1_reg, mv2_reg = regrid_to_lower_res(mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                    # Plotting
                    plot_set_5.plot(mv2_domain, mv1_domain, mv2_reg, mv1_reg, parameter)
        
        else:
            raise RuntimeError("Dimensions of two variables are difference. Abort")
