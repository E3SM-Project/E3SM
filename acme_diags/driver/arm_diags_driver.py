import os
import collections
import cdms2
import cdutil
import acme_diags
from acme_diags.driver import utils
#from acme_diags.metrics import mean
import numpy as np
import json
from acme_diags.plot.cartopy import arm_diags_plot
import acme_diags.derivations.acme

RefsTestMetrics = collections.namedtuple('RefsTestMetrics', ['refs', 'test', 'metrics', 'misc'])

def get_vars_funcs_for_derived_var(data_file, var):
    vars_to_func_dict =  acme_diags.derivations.acme.derived_variables[var]
    vars_in_file = set(data_file.variables)
    possible_vars = list(vars_to_func_dict.keys())  # ex: [('pr',), ('PRECC', 'PRECL')]

    for list_of_vars in possible_vars:
        if vars_in_file.issuperset(list_of_vars):
            # All of the variables (list_of_vars) are in data_file.
            # Return the corresponding dict.
            return {list_of_vars: vars_to_func_dict[list_of_vars]}

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def create_metrics(test, ref):
    """
    For this plotset, calculate the mean of the
    reference data and return a dict of that.
    """
    return {'test_mean': float(test.mean()),
            'ref_mean': float(ref.mean()),
            'test_std': float(test.std(ddof=1)),
            'ref_std': float(ref.std(ddof=1)),
            'rmse': float(rmse(test, ref)),
            'corr': float(np.corrcoef(test, ref)[0,1])
            }

def run_diag_diurnal_cycle(parameter):
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    run_type = parameter.run_type
    ref_path = parameter.reference_data_path

    seasons = ['DJF', 'MAM', 'JJA', 'SON']

    for region in regions:
        print("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            print("Season: {}".format(season))
            for var in variables:
                print('Variable: {}'.format(var))
                test_data = utils.dataset.Dataset(parameter, test=True)
                test = test_data.get_timeseries_variable(var, single_point = True)
                test.lat = test_data.get_static_variable('lat', var)
                test.lon = test_data.get_static_variable('lon', var)
                test_diurnal, lst = utils.diurnal_cycle.composite_diurnal_cycle(test, season, fft = False)

                parameter.viewer_descr[var] = getattr(test, 'long_name', var)
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)
                parameter.var_name = getattr(test, 'long_name', var)
                parameter.var_units = getattr(test, 'units', var)

                refs = []

                if 'armdiags' in ref_name:
                    if region != 'sgp': 
                        msg = 'Diurnal cycle of {} at Site: {} is not supported yet'.format(region, var)
                        raise RuntimeError(msg)
                    else:
                        ref_file_name = 'sgparmdiagsmondiurnalC1.c1.nc'
 
                        ref_file = os.path.join(ref_path,ref_file_name)
                        ref_data = cdms2.open(ref_file)
                        
                        if var == 'PRECT':
                            ref = ref_data('pr')*3600.*24 #Converting mm/second to mm/day"
                            ref.lat = test.lat
                            ref.lon = test.lon
                            ref_diurnal, lst = utils.diurnal_cycle.composite_diurnal_cycle(ref, season, fft = False)
                            ref.long_name = ref.standard_name
                            ref = ref_diurnal

                else:
                    ref_data = utils.dataset.Dataset(parameter, ref=True)
                    ref = ref_data.get_timeseries_variable(var, single_point = True)
                    ref.lat = test_data.get_static_variable('lat', var)
                    ref.lon = test_data.get_static_variable('lon', var)
                    ref_diurnal, lst = utils.diurnal_cycle.composite_diurnal_cycle(ref, season, fft = False)
                    ref = ref_diurnal


                refs.append(ref)

                metrics_dict = {}
                result = RefsTestMetrics(test=test_diurnal, refs=refs, metrics=None,misc = lst)
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict['unit'] = test.units
                parameter.output_file = '-'.join(
                            [ref_name, var, season, region])
                fnm = os.path.join(utils.general.get_output_dir(
                    parameter.current_set, parameter), parameter.output_file + '.json')
                with open(fnm, 'w') as outfile:
                     json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                # When running in a container, the paths are modified.
                fnm = os.path.join(utils.general.get_output_dir(parameter.current_set,
                    parameter, ignore_container=True), parameter.output_file + '.json')
                print('Metrics saved in: ' + fnm)
   
                arm_diags_plot.plot_diurnal_cycle(var, vars_to_data[season], parameter)


    return parameter


def run_diag_diurnal_cycle_zt(parameter):
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    run_type = parameter.run_type
    ref_path = parameter.reference_data_path

    seasons = ['ANNUALCYCLE']
    plevs = np.linspace(100,1000,37)

    for region in regions:
        print("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            print("Season: {}".format(season))
            for var in variables:
                print('Variable: {}'.format(var))
                test_data = utils.dataset.Dataset(parameter, test=True)
                test = test_data.get_timeseries_variable(var, single_point = True)
                test.lat = test_data.get_static_variable('lat', var)
                test.lon = test_data.get_static_variable('lon', var)
                if test.getLevel():
                    test_p = utils.general.convert_to_pressure_levels(test, plevs, test_data, var, season)
                    test_diurnal, lst = utils.diurnal_cycle.composite_diurnal_cycle(test_p, season, fft=False)

                parameter.viewer_descr[var] = getattr(test, 'long_name', var)
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)
                parameter.var_name = getattr(test, 'long_name', var)
                parameter.var_units = getattr(test, 'units', var)

                refs = []

                if 'armdiags' in ref_name:
                    if region == 'sgp': 
                        ref_file_name = 'sgparmdiagsmondiurnalclimC1.c1.nc'
                    elif region == 'nsa':
                        ref_file_name = region[:3]+'armdiagsmondiurnalclim' + 'C1.c1.nc'
                    else:
                        ref_file_name = region[:3]+'armdiagsmondiurnalclim' + region[3:5].upper()+'.c1.nc'
 
                    ref_file = os.path.join(ref_path,ref_file_name)
                    ref_data = cdms2.open(ref_file)
                    if var == 'CLOUD':
                        ref_var = ref_data('cl_p')
                        ref_var.long_name = "Cloud Fraction"
                        ref = ref_var
                        ref = np.reshape(ref,(12, 24, ref.shape[1]))
                        ref.ref_name = ref_name
                        ref.lat = test.lat
                        ref.lon = test.lon

                else:
                    ref_data = utils.dataset.Dataset(parameter, ref=True)
                    ref = ref_data.get_timeseries_variable(var, single_point = True)
                    ref.lat = ref_data.get_static_variable('lat', var)
                    ref.lon = ref_data.get_static_variable('lon', var)
                    if ref.getLevel():
                        ref_p = utils.general.convert_to_pressure_levels(ref, plevs, ref_data, var, season)
                        ref_diurnal, lst = utils.diurnal_cycle.composite_diurnal_cycle(ref_p, season, fft=False)
                    ref = ref_diurnal

                refs.append(ref)

                metrics_dict = {}
                result = RefsTestMetrics(test=test_diurnal, refs=refs, metrics=None,misc = lst)
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict['unit'] = test.units
                parameter.output_file = '-'.join(
                            [ref_name, var, season, region])
                fnm = os.path.join(utils.general.get_output_dir(
                    parameter.current_set, parameter), parameter.output_file + '.json')
                with open(fnm, 'w') as outfile:
                     json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                # When running in a container, the paths are modified.
                fnm = os.path.join(utils.general.get_output_dir(parameter.current_set,
                    parameter, ignore_container=True), parameter.output_file + '.json')
                print('Metrics saved in: ' + fnm)
   
            if season == 'ANNUALCYCLE':
                arm_diags_plot.plot_diurnal_cycle_zt(var, vars_to_data[season], parameter)

    return parameter


def run_diag_annual_cycle(parameter):
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    run_type = parameter.run_type
    ref_path = parameter.reference_data_path

    seasons = ['ANNUALCYCLE']
    plevs = np.linspace(100,1000,37)

    for region in regions:
        # The regions that are supported are in acme_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        print("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            print("Season: {}".format(season))
            for var in variables:
                print('Variable: {}'.format(var))
                test_data = utils.dataset.Dataset(parameter, test=True)
                test = test_data.get_climo_variable(var, season)
                if test.getLevel():
                    test_p = utils.general.convert_to_pressure_levels(test, plevs, test_data, var, season)
                    test = utils.climo.climo(test_p, season)

                parameter.viewer_descr[var] = getattr(test, 'long_name', var)
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)
                parameter.var_name = getattr(test, 'long_name', var)
                parameter.var_units = getattr(test, 'units', var)

                refs = []

                if 'armdiags' in ref_name:
                    if region == 'sgp': 
                        ref_file = os.path.join(ref_path,'sgparmdiagsmonC1.c1.nc')
                    elif region == 'nsa':
                        ref_file = os.path.join(ref_path,region[:3]+'armdiagsmonclim' + 'C1.c1.nc')
                    else:
                        ref_file = os.path.join(ref_path,region[:3]+'armdiagsmonclim' + region[3:5].upper()+'.c1.nc')
                    ref_data = cdms2.open(ref_file)
                    vars_funcs = get_vars_funcs_for_derived_var(ref_data, var)
                    target_var = list(vars_funcs.keys())[0][0]
                    ref_var = ref_data(target_var)
                    ref_var.long_name = ref_var.standard_name 
                    ref = vars_funcs[(target_var,)](utils.climo.climo(ref_var,season))

                else:
                    ref_data = utils.dataset.Dataset(parameter, ref=True)
                    ref = ref_data.get_climo_variable(var, season)
                    if ref.getLevel():
                        ref_p = utils.general.convert_to_pressure_levels(ref, plevs, ref_data, var, season)
                        ref = utils.climo.climo(ref_p, season)
                ref_domain = utils.general.select_point(region, ref)
                ref.ref_name = ref_name
                refs.append(ref_domain)

                test_domain = utils.general.select_point(region, test)

                metrics_dict = create_metrics(test_domain,ref_domain)

                result = RefsTestMetrics(test=test_domain, refs=refs, metrics=metrics_dict, misc = None)
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict['unit'] = test.units
                parameter.output_file = '-'.join(
                            [ref_name, var, season, region])
                fnm = os.path.join(utils.general.get_output_dir(
                    parameter.current_set, parameter), parameter.output_file + '.json')
                with open(fnm, 'w') as outfile:
                     json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                # When running in a container, the paths are modified.
                fnm = os.path.join(utils.general.get_output_dir(parameter.current_set,
                    parameter, ignore_container=True), parameter.output_file + '.json')
                print('Metrics saved in: ' + fnm)
   
            if season == 'ANNUALCYCLE':
                arm_diags_plot.plot_annual_cycle(var, vars_to_data[season], parameter)


    return parameter

def run_diag_convection_onset(parameter):
    regions = parameter.regions
    ref_name = parameter.ref_name
    run_type = parameter.run_type
    ref_path = parameter.reference_data_path
     #Read in observation data
    variables =['PRECT', 'TMQ']

    for region in regions:
        # The regions that are supported are in acme_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        print("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        test_data = utils.dataset.Dataset(parameter, test=True)

        test_pr = test_data.get_timeseries_variable('PRECT',single_point = True)/24.0
        test_prw = test_data.get_timeseries_variable('TMQ',single_point = True)

        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)

        if 'armdiags' in ref_name:
            if region == 'sgp':
                ref_file_name = 'sgparmdiags1hrC1.c1.nc'
            else:
                ref_file_name = region[:3]+'armdiags1hr' + region[3:5].upper()+'.c1.nc'
            ref_file = os.path.join(ref_path, ref_file_name)
            ref_data = cdms2.open(ref_file)
            ref_pr = ref_data('pr')        #mm/hr
            ref_pr[ref_pr<-900] = np.nan
            ref_prw = ref_data('prw')      #mm 
            ref_prw[ref_prw<-900] = np.nan
        else:
            ref_data = utils.dataset.Dataset(parameter, ref=True)
            ref_pr = test_data.get_timeseries_variable('PRECT', single_point = True)/24.0
            ref_prw = test_data.get_timeseries_variable('TMQ', single_point = True)
        parameter.output_file = '-'.join(
                            [ref_name, 'convection-onset', region]) 
       
        arm_diags_plot.plot_convection_onset_statistics(test_pr, test_prw, ref_pr, ref_prw,parameter, region)

    return parameter


def run_diag(parameter):

    if parameter.diags_set == 'annual_cycle':
        return run_diag_annual_cycle(parameter)
    elif parameter.diags_set == 'diurnal_cycle':
        return run_diag_diurnal_cycle(parameter)
    elif parameter.diags_set == 'diurnal_cycle_zt':
        return run_diag_diurnal_cycle_zt(parameter)
    elif parameter.diags_set == 'pdf_daily':
        return run_diag_pdf_daily(parameter)
    elif parameter.diags_set == 'convection_onset':
        return run_diag_convection_onset(parameter)
    else:
        raise Exception('Invalid diags_set={}'.format(parameter.diags_set))
