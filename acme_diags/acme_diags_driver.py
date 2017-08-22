#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import json
import getpass
import datetime
import importlib
import cdp.cdp_run
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_parameter import ACMEParameter
from acme_diags.acme_viewer import create_viewer
from acme_diags.driver.utils import get_set_name

def _get_default_diags(set_num, dataset):
    """Returns the path from the json corresponding to set_num"""
    set_num = get_set_name(set_num)
    folder = '{}'.format(set_num)
    fnm = '{}_{}.json'.format(set_num, dataset)
    pth = os.path.join(sys.prefix, 'share', 'acme_diags', folder, fnm)
    print('Using {} for {}'.format(pth, set_num))
    if not os.path.exists(pth):
        raise RuntimeError("Plotting via set '{}' not supported, file {} not installed".format(set_num, fnm))
    return pth

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
    args = parser.view_args()

    if args.parameters and not args.other_parameters:  # -p only
        original_parameter = parser.get_orig_parameters(default_vars=False)

        if not hasattr(original_parameter, 'sets'):
            original_parameter.sets = ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']
            
        # load the default jsons
        default_jsons_paths = []

        for set_num in original_parameter.sets:
            datasets = ['ACME']
            if hasattr(original_parameter, 'datasets'):
                datasets = original_parameter.datasets
            for ds in datasets:
                default_jsons_paths.append(_get_default_diags(set_num, ds))
        other_parameters = parser.get_other_parameters(files_to_open=default_jsons_paths, check_values=False)
        # Don't put the sets from the Python parameters to the default.
        # Ex. if sets=[5, 7] in the Python parameters, don't change sets in the default jsons like lat_lon_AMWG_default.json
        vars_to_ignore = ['sets']
        parameters = parser.get_parameters(orig_parameters=original_parameter, other_parameters=other_parameters, vars_to_ignore=vars_to_ignore)

    elif not args.parameters and args.other_parameters:  # -d only
        other_parameters = parser.get_other_parameters(check_values=True)
        parameters = parser.get_parameters(other_parameters=other_parameters)

    elif args.parameters and args.other_parameters:  # -p and -d
        original_parameter = parser.get_orig_parameters(default_vars=False)
        other_parameters = parser.get_other_parameters(check_values=False)
        parameters = parser.get_parameters(orig_parameters=original_parameter, other_parameters=other_parameters)

    else:
        raise RuntimeError('You tried running the diags without -p and/or -d')

    dt = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    for p in parameters:
        # needed for distributed running
        # chown of all generated files to the user who ran the diags
        p.user = getpass.getuser()

        if not hasattr(p, 'results_dir'):
            p.results_dir = '{}-{}'.format('acme_diags_results', dt)


    if not os.path.exists(parameters[0].results_dir):
        os.makedirs(parameters[0].results_dir, 0775)
    if parameters[0].multiprocessing:
        cdp.cdp_run.multiprocess(run_diag, parameters)
    elif parameters[0].distributed:
        cdp.cdp_run.distribute(run_diag, parameters)
    else:
        cdp.cdp_run.serial(run_diag, parameters)

    pth = os.path.join(parameters[0].results_dir, 'viewer')
    if not os.path.exists(pth):
        os.makedirs(pth)

    create_viewer(pth, parameters, parameters[0].output_format[0])
