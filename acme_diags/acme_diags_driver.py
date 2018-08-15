#!/usr/bin/env python
from __future__ import print_function

import os
# Must be done before any CDAT library is called.
os.environ['UVCDAT_ANONYMOUS_LOG'] = 'no'
# Needed for when using hdf5 >= 1.10.0,
# without this, errors are thrown on Edison compute nodes.
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
# Used by numpy, causes too many threads to spawn otherwise.
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import sys
import getpass
import datetime
import importlib
import traceback
import subprocess
import cdp.cdp_run
import acme_diags
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_viewer import create_viewer
from acme_diags.driver import utils


def _get_default_diags(set_num, dataset, run_type):
    """
    Returns the path for the default diags corresponding to set_num.
    """
    set_num = utils.general.get_set_name(set_num)

    if dataset and run_type == 'model_vs_obs':  # either 'ACME' or 'AMWG', use the jsons
        fnm = '{}_{}.json'.format(set_num, dataset)
    else:  # use the cfgs
        fnm = '{}_{}.cfg'.format(set_num, run_type)

    folder = '{}'.format(set_num)
    pth = os.path.join(acme_diags.INSTALL_PATH, folder, fnm)
    print('Using {} for {}'.format(pth, set_num))
    if not os.path.exists(pth):
        raise RuntimeError(
            "Plotting via set '{}' not supported, file {} not installed".format(set_num, fnm))
    return pth


def _collapse_results(parameters):
    """
    When using cdp_run, parameters is a list of lists: [[Parameters], ...].
    Make this just a list: [Parameters, ...].
    """
    output_parameters = []

    for p1 in parameters:
        if isinstance(p1, list):
            for p2 in p1:
                output_parameters.append(p2)
        else:
            output_parameters.append(p1)

    return output_parameters


def _save_env_yml(results_dir):
    """
    Save the yml to recreate the environment in results_dir.
    """
    cmd = 'conda env export'
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, err = p.communicate()

    if err:
        print('Error when creating env yml file:')
        print(err)

    else:
        fnm = os.path.join(results_dir, 'environment.yml')
        with open(fnm, 'w') as f:
            f.write(output.decode('utf-8'))

        print('Saved environment yml file to: {}'.format(fnm))


def _save_parameter_files(results_dir, parser):
    """
    Save the command line arguments used, and any py or cfg files.
    """
    cmd_used = ' '.join(sys.argv)
    fnm = os.path.join(results_dir, 'cmd_used.txt')
    with open(fnm, 'w') as f:
        f.write(cmd_used)
    print('Saved command used to: {}'.format(fnm))

    args = parser.view_args()

    if hasattr(args, 'parameters') and args.parameters:
        fnm = args.parameters
        if not os.path.isfile(fnm):
            print('File does not exist: {}'.format(fnm))
        else:
            with open(fnm, 'r') as f:
                contents = ''.join(f.readlines())
            # Remove any path, just keep the filename
            new_fnm = fnm.split('/')[-1]
            new_fnm = os.path.join(results_dir, new_fnm)
            with open(new_fnm, 'w') as f:
                f.write(contents)
            print('Saved py file to: {}'.format(new_fnm))

    if hasattr(args, 'other_parameters') and args.other_parameters:
        fnm = args.other_parameters[0]
        if not os.path.isfile(fnm):
            print('File does not exist: {}'.format(fnm))
        else:
            with open(fnm, 'r') as f:
                contents = ''.join(f.readlines())
            # Remove any path, just keep the filename
            new_fnm = fnm.split('/')[-1]
            new_fnm = os.path.join(results_dir, new_fnm)
            with open(new_fnm, 'w') as f:
                f.write(contents)
            print('Saved cfg file to: {}'.format(new_fnm))


def save_provenance(results_dir, parser):
    """
    Store the provenance in results_dir.
    """
    results_dir = os.path.join(results_dir, 'prov')
    if not os.path.exists(results_dir):
        os.makedirs(results_dir, 0o775)

    try:
        _save_env_yml(results_dir)
    except:
        traceback.print_exc()

    _save_parameter_files(results_dir, parser)


def run_diag(parameters):
    """
    For a single set of parameters, run the corresponding diags.
    """
    results = []
    for pset in parameters.sets:
        set_name = utils.general.get_set_name(pset)

        parameters.current_set = set_name
        mod_str = 'acme_diags.driver.{}_driver'.format(set_name)
        try:
            module = importlib.import_module(mod_str)
            single_result = module.run_diag(parameters)
            print('')
            results.append(single_result)
        except:
            print('Error in {}'.format(mod_str))
            traceback.print_exc()
            if parameters.debug:
                sys.exit()

    return results


def main():
    parser = ACMEParser()
    args = parser.view_args()

    # There weren't any arguments defined
    if not any(getattr(args, arg) for arg in vars(args)):
        parser.print_help()
        sys.exit()

    if args.parameters and not args.other_parameters:  # -p only
        cmdline_parameter = parser.get_cmdline_parameters()
        # If just a -p with no command line parameters, check the py for errors
        # Otherwise don't check for errors, the command line args might have some missing ones
        check_values = True if not cmdline_parameter else False
        original_parameter = parser.get_orig_parameters(check_values)

        # Special case, until selector parameter is added.
        if hasattr(cmdline_parameter, 'sets'):
            original_parameter.sets = cmdline_parameter.sets
        elif not hasattr(original_parameter, 'sets'):
            original_parameter.sets = utils.general.SET_NAMES

        # load the default cfg files
        default_diags_paths = []
        for set_num in original_parameter.sets:
            run_type = 'model_vs_obs'
            dataset = ''
            if hasattr(original_parameter, 'dataset'):
                dataset = original_parameter.dataset
            if hasattr(original_parameter, 'run_type'):
                run_type = original_parameter.run_type
            default_diags_paths.append(_get_default_diags(set_num, dataset, run_type))

        other_parameters = parser.get_other_parameters(files_to_open=default_diags_paths)

        # Don't put the sets from the Python parameters to each of the parameters.
        # Ex. if sets=[5, 7] in the Python parameters, don't change sets in the
        # default jsons like lat_lon_ACME_default.json
        vars_to_ignore = ['sets']
        parameters = parser.get_parameters(cmdline_parameters=cmdline_parameter,
            orig_parameters=original_parameter, other_parameters=other_parameters,
            vars_to_ignore=vars_to_ignore, cmd_default_vars=False)

    elif not args.parameters and args.other_parameters:  # -d only
        cmdline_parameter = parser.get_cmdline_parameters()
        other_parameters = parser.get_other_parameters()
        parameters = parser.get_parameters(cmdline_parameters=cmdline_parameter, 
            other_parameters=other_parameters, cmd_default_vars=False)

    else:  # -p and -d, or just command line arguments
        parameters = parser.get_parameters(cmd_default_vars=False)

    dt = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    for p in parameters:
        # needed for distributed running
        # chown of all generated files to the user who ran the diags
        p.user = getpass.getuser()

        if not hasattr(p, 'results_dir'):
            p.results_dir = '{}-{}'.format('e3sm_diags_results', dt)

    if not os.path.exists(parameters[0].results_dir):
        os.makedirs(parameters[0].results_dir, 0o775)
    if not parameters[0].no_viewer:  # Only save provenance for full runs.
        save_provenance(parameters[0].results_dir, parser)

    if parameters[0].multiprocessing:
        parameters = cdp.cdp_run.multiprocess(run_diag, parameters)
    elif parameters[0].distributed:
        parameters = cdp.cdp_run.distribute(run_diag, parameters)
    else:
        parameters = cdp.cdp_run.serial(run_diag, parameters)

    parameters = _collapse_results(parameters)

    if not parameters:
        print('There was not a single valid diagnostics run, no viewer created.')
    else:        
        if parameters[0].no_viewer:
            print('Viewer not created because the no_viewer parameter is True.')
        else:
            pth = os.path.join(parameters[0].results_dir, 'viewer')
            if not os.path.exists(pth):
                os.makedirs(pth)
            create_viewer(pth, parameters, parameters[0].output_format[0])


if __name__ == '__main__':
    main()
