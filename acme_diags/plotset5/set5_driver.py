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


def make_parameters(orginal_parameter):
    #f_data = open('set5_diags.json').read()
    #f_data = open('set5_diags_MERRA.json').read()
    f_data = open('set5_diags_HADISST.json').read()
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

# Below three functions are copied from uvcmetrics.
def mask_by( var, maskvar, lo=None, hi=None ):
    """masks a variable var to be missing except where maskvar>=lo and maskvar<=hi.  That is, the missing-data mask is True where maskvar<lo or maskvar>hi or where it was True on input. For lo and hi, None means to omit the constrint, i.e. lo=-infinity or hi=infinity. var is changed and returned; we don't make a new variable.
var and maskvar: dimensioned the same variables.  
lo and hi: scalars.
    """
    if lo is None and hi is None:
        return var
    if lo is None and hi is not None:
        maskvarmask = maskvar>hi
    elif lo is not None and hi is None:
        maskvarmask = maskvar<lo
    else:
        maskvarmask = (maskvar<lo) | (maskvar>hi)
    if var.mask is False:
        newmask = maskvarmask
    else:
        newmask = var.mask | maskvarmask
    var.mask = newmask
    return var

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
#    reference_data_set = parameter.reference_data_set # observation

    test_data_path = parameter.test_data_path
#    test_data_set = parameter.test_data_set # model
#    filename1=test_data_path + test_data_set
#    filename2=reference_data_path + reference_data_set

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
        for filename in fnmatch.filter(ref_files, '*'+ref_name+'_'+it+'*'):
            print filename
            filename2 = filename
    

        if var == 'PRECT':
            f_mod = cdms2.open(filename1)
            f_obs = cdms2.open(filename2)
            if ref_name == 'TRMM': #TRMM region following AMWG  latitude=(-38,38)
                mv1 = (f_mod('PRECC', latitude=(-38,38)) + f_mod('PRECL', latitude=(-38,38)))*3600.0*24.0*1000.0
                mv2 = f_obs(var, latitude=(-38,38))
            else:
                mv1 = (f_mod('PRECC') + f_mod('PRECL'))*3600.0*24.0*1000.0
                mv2 = f_obs(var)
            mv1.units = 'mm/day'
        elif var == 'PREH2O':
            f_in = cdms2.open(filename1)
            mv1 = f_in('TMQ')
            f_in = cdms2.open(filename2)
            mv2 = f_in(var)
        elif var == 'SST':
            f_in = cdms2.open(filename1)
            TS = f_in('TS')
            TS = TS -273.15
            OCNFRAC = f_in('OCNFRAC')
            mv1 = mask_by(TS, OCNFRAC, lo = 0.9) #following AMWG
            f_in = cdms2.open(filename2)
            mv2 = f_in(var)

            mv1_reg, mv2_reg = regrid_to_lower_res(mv1, mv2, parameter.regrid_tool, parameter.regrid_method)
            mv1_reg=mv1_reg(squeeze=1)
            #create common mask for SST
            land_mask = MV2.logical_or(mv1_reg.mask,mv2_reg.mask)

            # Note, NCAR plots native mv1, mv2, but mean is from newly masker mv1 and mv2.while we plot everything already masked for now.Need to modify plotting routine for accomodating printing different sets of means to plots.
            mv1_new = MV2.masked_where(land_mask,mv1_reg)
            mv2_new = MV2.masked_where(land_mask,mv2_reg)

            plot_set_5.plot(mv2, mv1, mv2_new, mv1_new, parameter)
        

        else:
            f_in = cdms2.open(filename1)
            mv1 = f_in(var)
            f_in = cdms2.open(filename2)
            mv2 = f_in(var)
        parameter.output_file = '_'.join([ref_name,it])
        parameter.main_title = str(' '.join([var,it]))
        print parameter.output_file

        #for variables without z axis:
        if mv1.getLevel()==None and mv2.getLevel()==None:
            #regrid towards lower resolution of two variables for calculating difference
            if var != 'SST':
                mv1_reg, mv2_reg = regrid_to_lower_res(mv1, mv2, parameter.regrid_tool, parameter.regrid_method)
                plot_set_5.plot(mv2, mv1, mv2_reg, mv1_reg, parameter)

        
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
                    print( 'Vertical level is neither hybrid nor pressure. Abort')
                    quit()
           
                if filename == filename1:
                    mv1_p = mv_p
        
                if filename == filename2:
                    mv2_p = mv_p

            for ilev in range(len(plev)):
                mv1 = mv1_p[:,ilev,]
                mv2 = mv2_p[:,ilev,]
    
                parameter.main_title = str(' '.join([var, str(int(plev[ilev])), 'mb', it]))
                parameter.output_file = str('_'.join([ref_name,it, str(int(plev[ilev]))]))

                # Regrid towards lower resolution of two variables for calculating difference
                mv1_reg, mv2_reg = regrid_to_lower_res(mv1, mv2, parameter.regrid_tool, parameter.regrid_method)
        
                # Ploting
                plot_set_5.plot(mv2, mv1, mv2_reg, mv1_reg, parameter)
        else:
            print "Dimensions of two variables are difference. Abort"
