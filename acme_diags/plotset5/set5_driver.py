#!/usr/bin/env python
import json
import copy
import glob
import os
import fnmatch
import numpy
import cdutil
import genutil
import cdms2
import MV2
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_parameter import ACMEParameter
import acme_diags.plotting.set5.plot
from acme_diags.derivations import acme
from acme_diags.derivations.default_regions import regions_specs
from acme_diags.viewer import OutputViewer


def make_parameters(orginal_parameter):
    """ Create multiple parameters given a list of 
    parameters in a json and an original parameter """
    #f_data = open('set5_diags.json').read()
    f_data = open('set5_diags_GPCP.json').read()
    #f_data = open('set5_diags_MERRA.json').read()
    #f_data = open('set5_diags_HADISST.json').read()
    #f_data = open('set5_diags_NVA.json').read()
    json_file = json.loads(f_data)

    parameters = []
    for key in json_file:
        for single_run in json_file[key]:
            #p = copy.deepcopy(orginal_parameter)
            p = ACMEParameter()
            for attr_name in single_run:
                setattr(p, attr_name, single_run[attr_name])
            # Add attributes of original_parameter to p
            print orginal_parameter.__dict__
            p.__dict__.update(orginal_parameter.__dict__)
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

def add_page_and_top_row(viewer, parameters):
    """ Setup for OutputViewer """
    col_labels = ['Description']
    seasons = []
    for p in parameters:
        for s in p.season:
            if s not in seasons:
                seasons.append(s)
    for s in seasons:
        col_labels.append(s)

    viewer.add_page("Set 5", col_labels)


parser = ACMEParser()
original_parameter = parser.get_parameter(default_vars=False)
parameters = make_parameters(original_parameter)
viewer = OutputViewer(path=parameters[0].case_id, index_name='index name')
add_page_and_top_row(viewer, parameters)

for parameter in parameters:
    viewer.add_group(parameter.reference_name)

    reference_data_path = parameter.reference_data_path
    test_data_path = parameter.test_data_path

    var = parameter.variables
    seasons = parameter.season
    ref_name = parameter.ref_name
    test_name = parameter.test_name
    regions = parameter.region
    #domain = regions_specs[region]
    for season in seasons:
        test_files = glob.glob(os.path.join(test_data_path,'*'+test_name+'*.nc'))
        for filename in fnmatch.filter(test_files, '*'+season+'*'):
            print filename
            filename1 = filename

        ref_files = glob.glob(os.path.join(reference_data_path,'*'+ref_name+'*.nc'))
        for filename in fnmatch.filter(ref_files, '*'+season+'*'):
            print filename
            filename2 = filename

        f_mod = cdms2.open(filename1)
        f_obs = cdms2.open(filename2)

        # domain can pass in process_derived_var after Charles fix cdutil.domain'unit problem
        #mv1 = acme.process_derived_var(var, acme.derived_variables, f_mod, domain)
        #mv2 = acme.process_derived_var(var, acme.derived_variables, f_obs, domain)
        try:
            mv1 = f_mod(var)
        except:
            mv1 = acme.process_derived_var(var, acme.derived_variables, f_mod)
        try:
            mv2 = f_obs(var)
        except:
            mv2 = acme.process_derived_var(var, acme.derived_variables, f_obs)
        print regions,"regions"

        # Temporary fix to bypass bug (Zeshawn or Jill will fix)
        mv1.units = parameter.test_units
        mv2.units = parameter.test_units
         
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
                mv1_domain = mv1(domain)
                mv2_domain = mv2(domain)
                mv1_domain.units = mv1.units
                mv2_domain.units = mv2.units
                parameter.output_file = '_'.join([ref_name, season,region])
                parameter.main_title = str(' '.join([var, season]))

                #else:
                #
                #    parameter.output_file = '_'.join([ref_name, season])
                #    parameter.main_title = str(' '.join([var, season]))

                #parameter.output_file = '_'.join([ref_name, season,region])
                #parameter.main_title = str(' '.join([var, season]))


                #regrid towards lower resolution of two variables for calculating difference
                mv1_reg, mv2_reg = regrid_to_lower_res(mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                if var is 'SST': #special case

                    land_mask = MV2.logical_or(mv1_reg.mask, mv2_reg.mask)
                    mv1_reg = MV2.masked_where(land_mask, mv1_reg)
                    mv2_reg = MV2.masked_where(land_mask, mv2_reg)

                print 'hi1'
                acme_diags.plotting.set5.plot.plot(mv2_domain, mv1_domain, mv2_reg, mv1_reg, parameter)
                viewer.add_row('%s %s' % (var, region), 'Description for %s' % var, file_name=parameter.case_id + '/' + parameter.output_file)


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

                    parameter.output_file = '_'.join([ref_name, season, region, var, str(int(plev[ilev])), 'mb'])
                    parameter.main_title = str(' '.join([var, str(int(plev[ilev])), 'mb', season]))

                    # Regrid towards lower resolution of two variables for calculating difference
                    mv1_reg, mv2_reg = regrid_to_lower_res(mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                    # Plotting
                    print 'hi2'
                    acme_diags.plotting.set5.plot.plot(mv2_domain, mv1_domain, mv2_reg, mv1_reg, parameter)
                    viewer.add_row('%s %s %s' % (var, plev[ilev], region), 'Description for %s' % var, file_name=parameter.case_id + '/' + parameter.output_file)

        else:
            raise RuntimeError("Dimensions of two variables are difference. Abort")

viewer.generate_viewer()
