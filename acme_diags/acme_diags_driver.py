import os
import sys
import json
import datetime
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_parameter import ACMEParameter
from acme_diags.acme_viewer import create_viewer

def _get_default_diags(set_num):
    """Get the data from the json corresponding to set_num"""
    folder = 'set{}'.format(set_num)
    fnm = 'set{}_diags_AMWG_default.json'.format(set_num)
    pth = os.path.join(sys.prefix, 'share', 'acme_diags', folder, fnm)
    with open(pth) as json_file:
        json_data = json.loads(json_file.read())
    return json_data

def make_parameters(original_parameter, vars_to_ignore=[]):
    """Create multiple parameters given a list of
    parameters in a json and an original parameter"""

    if hasattr(original_parameter, 'custom_diags'):
        with open(original_parameter.custom_diags) as json_file:
            json_data = json.loads(json_file.read())
    else:
        json_data = {'': []}  # the first key doesn't hold any value, so it's ''
        for set_num in original_parameter.sets:
            default_set_runs = _get_default_diags(set_num)
            for _, set_runs in default_set_runs.iteritems():
                for single_run in set_runs:
                    json_data[''].append(single_run)
    parameters = []
    for key in json_data:
        for single_run in json_data[key]:
            p = ACMEParameter()
            for attr_name in single_run:
                setattr(p, attr_name, single_run[attr_name])

            # Add attributes of original_parameter to p
            for var in original_parameter.__dict__:
                if var not in vars_to_ignore:
                    p.__dict__[var] = original_parameter.__dict__[var]
            p.check_values()
            parameters.append(p)
    return parameters

if __name__ == '__main__':
    parser = ACMEParser()
    original_parameter = parser.get_parameter(default_vars=False)
    if not hasattr(original_parameter, 'results_dir'):
        dt = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        original_parameter.results_dir = '{}-{}'.format('acme_diags_results', dt)

    # if user wants all of the default diags for the sets (ex: sets=[5, 7]), then
    # don't overwrite the sets keyword in the default json files.
    ignore_vars = [] if hasattr(original_parameter, 'custom_diags') else ['sets']
    parameters = make_parameters(original_parameter, vars_to_ignore=ignore_vars)
    
    for parameter in parameters:
        for pset in parameter.sets:
            pset = str(pset)

            parameter.current_set = pset
            if pset == '3' or pset == 'zonal_mean_xy':
                from acme_diags.driver.zonal_mean_xy_driver import run_diag
            elif pset == '4' or pset == 'zonal_mean_2d':
                from acme_diags.driver.zonal_mean_2d_driver import run_diag
            elif pset == '5' or pset == 'lat_lon':
                from acme_diags.driver.lat_lon_driver import run_diag
            elif pset == '7' or pset == 'polar':
                from acme_diags.driver.polar_driver import run_diag
            elif pset == '13' or pset == 'cosp_histogram':
                from acme_diags.driver.cosp_histogram_driver import run_diag

            else:
                print('Plot set {} is not supported yet. Please give us time.'.format(pset))
                continue
            print('Starting to run ACME diags')
            run_diag(parameter)
    
    pth = os.path.join(original_parameter.results_dir, 'viewer')
    if not os.path.exists(pth):
        os.makedirs(pth)

    create_viewer(pth, parameters, parameters[0].output_format[0])
    