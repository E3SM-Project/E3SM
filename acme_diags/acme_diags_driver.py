#!/usr/bin/env python
from __future__ import print_function

import os
# Must be done before any CDAT library is called
if 'UVCDAT_ANONYMOUS_LOG' not in os.environ:
    os.environ['UVCDAT_ANONYMOUS_LOG'] = 'no'

import sys
import getpass
import datetime
import importlib
import traceback
import cdp.cdp_run
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_viewer import create_viewer
from acme_diags.driver.utils import get_set_name, SET_NAMES


def _get_default_diags(set_num, dataset, run_type):
    """Returns the path for the default diags corresponding to set_num"""
    set_num = get_set_name(set_num)

    if dataset and run_type in ['model_vs_obs', 'obs_vs_model']:  # either 'ACME' or 'AMWG', use the jsons
        fnm = '{}_{}.json'.format(set_num, dataset)
    else:  # use the cfgs
        fnm = '{}_{}.cfg'.format(set_num, run_type)

    folder = '{}'.format(set_num)
    pth = os.path.join(sys.prefix, 'share', 'acme_diags', folder, fnm)
    print('Using {} for {}'.format(pth, set_num))
    if not os.path.exists(pth):
        raise RuntimeError(
            "Plotting via set '{}' not supported, file {} not installed".format(set_num, fnm))
    return pth


def _collaspse_results(parameters):
    """When using cdp_run, parameters is a list of lists: [[Parameters], ...].
       Make this just a list: [Parameters]."""
    output_parameters = []

    for p1 in parameters:
        if isinstance(p1, list):
            for p2 in p1:
                output_parameters.append(p2)
        else:
            output_parameters.append(p1)

    return output_parameters


def run_diag(parameters):
    """For a single set of parameters, run the corresponding diags."""
    results = []
    for pset in parameters.sets:
        set_name = get_set_name(pset)

        parameters.current_set = set_name
        mod_str = 'acme_diags.driver.{}_driver'.format(set_name)
        try:
            module = importlib.import_module(mod_str)
            print('Starting to run ACME Diagnostics.')
            single_result = module.run_diag(parameters)
            print('')
            results.append(single_result)
        except Exception as e:
            print('Error in {}'.format(mod_str))
            traceback.print_exc()
            if parameters.debug:
                sys.exit()

    return results


def main():
    parser = ACMEParser()
    args = parser.view_args()

    if args.parameters and not args.other_parameters:  # -p only
        cmdline_parameter = parser.get_cmdline_parameters(default_vars=False, cmd_default_vars=False)
        # If just a -p with no command line parameters, check the py for errors
        if not cmdline_parameter:
            original_parameter = parser.get_orig_parameters(default_vars=False, cmd_default_vars=False)
        # Otherwise don't check for errors, the command line args might have some missing ones
        else:
            original_parameter = parser.get_orig_parameters(default_vars=False, cmd_default_vars=False, check_values=False)

        if not hasattr(original_parameter, 'sets'):
            original_parameter.sets = SET_NAMES

        # load the default jsons
        default_diags_paths = []

        for set_num in original_parameter.sets:
            run_type = 'model_vs_obs'
            dataset = ''
            if hasattr(original_parameter, 'dataset'):
                dataset = original_parameter.dataset
            if hasattr(original_parameter, 'run_type'):
                run_type = original_parameter.run_type
            default_diags_paths.append(_get_default_diags(set_num, dataset, run_type))
        other_parameters = parser.get_other_parameters(
            files_to_open=default_diags_paths, check_values=False, cmd_default_vars=False)
        # Don't put the sets from the Python parameters to each of the parameters.
        # Ex. if sets=[5, 7] in the Python parameters, don't change sets in the
        # default jsons like lat_lon_ACME_default.json
        vars_to_ignore = ['sets']
        parameters = parser.get_parameters(cmdline_parameters=cmdline_parameter,
            orig_parameters=original_parameter, other_parameters=other_parameters, vars_to_ignore=vars_to_ignore)

    elif not args.parameters and args.other_parameters:  # -d only
        cmdline_parameter = parser.get_cmdline_parameters(default_vars=False, cmd_default_vars=False)
        other_parameters = parser.get_other_parameters(check_values=True, cmd_default_vars=False)
        parameters = parser.get_parameters(cmdline_parameters=cmdline_parameter, other_parameters=other_parameters)

    elif args.parameters and args.other_parameters:  # -p and -d
        cmdline_parameter = parser.get_cmdline_parameters(default_vars=False, cmd_default_vars=False)
        original_parameter = parser.get_orig_parameters(default_vars=False, cmd_default_vars=False)
        other_parameters = parser.get_other_parameters(check_values=False, cmd_default_vars=False)
        parameters = parser.get_parameters(cmdline_parameters=cmdline_parameter,
            orig_parameters=original_parameter, other_parameters=other_parameters)

    else:  # command line args without -p or -d
        parameters = parser.get_parameters(cmd_default_vars=False)

    dt = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    for p in parameters:
        # needed for distributed running
        # chown of all generated files to the user who ran the diags
        p.user = getpass.getuser()

        if not hasattr(p, 'results_dir'):
            p.results_dir = '{}-{}'.format('acme_diags_results', dt)

    if not os.path.exists(parameters[0].results_dir):
        os.makedirs(parameters[0].results_dir, 0o775)
    if parameters[0].multiprocessing:
        parameters = cdp.cdp_run.multiprocess(run_diag, parameters)
    elif parameters[0].distributed:
        parameters = cdp.cdp_run.distribute(run_diag, parameters)
    else:
        parameters = cdp.cdp_run.serial(run_diag, parameters)

    parameters = _collaspse_results(parameters)

    if not parameters:
        print('There was not a single valid diagnostics run, no viewer created.')

    elif parameters[0].no_viewer:
        print('Viewer not created because the no_viewer parameter is True.')

    else:
        pth = os.path.join(parameters[0].results_dir, 'viewer')
        if not os.path.exists(pth):
            os.makedirs(pth)
        create_viewer(pth, parameters, parameters[0].output_format[0])


if __name__ == '__main__':
    main()

