#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import json
import getpass
import datetime
import importlib
from dask.distributed import Client
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_parameter import ACMEParameter
from acme_diags.acme_viewer import create_viewer
from acme_diags.driver.utils import get_set_name

def _get_default_diags(set_num):
    """Get the data from the json corresponding to set_num"""
    set_num = get_set_name(set_num)
    folder = '{}'.format(set_num)
    fnm = '{}_AMWG_default.json'.format(set_num)
    pth = os.path.join(sys.prefix, 'share', 'acme_diags', folder, fnm)
    if not os.path.exists(pth):
        raise RuntimeError("Plotting via set '{}' not supported".format(set_num))
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
            p.user = getpass.getuser()
            parameters.append(p)
    return parameters


def run_diag(parameters):
    """For a single set of parameters, run the corresponding diags."""
    for pset in parameters.sets:
        set_name = get_set_name(pset)

        parameters.current_set = set_name
        mod_str = 'acme_diags.driver.{}_driver'.format(set_name)
        try:
            module = importlib.import_module(mod_str)
        except ImportError as e:
            print('Set {} is not supported yet. Please give us time.'.format(set_name))
            continue
        print('Starting to run ACME Diagnostics.')
        module.run_diag(parameters)

            
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
   
    if not os.path.exists(original_parameter.results_dir):
        os.makedirs(original_parameter.results_dir, 0775)

    if not parameters[0].distributed:
        for p in parameters:
            run_diag(p)
    else:
        client = Client(parameters[0].scheduler_addr)
        results = client.map(run_diag, parameters)
        client.gather(results)
        client.close()

    pth = os.path.join(original_parameter.results_dir, 'viewer')
    if not os.path.exists(pth):
        os.makedirs(pth)

    create_viewer(pth, parameters, parameters[0].output_format[0])
    
