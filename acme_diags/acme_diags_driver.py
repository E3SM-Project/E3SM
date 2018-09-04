#!/usr/bin/env python
from __future__ import print_function

import os
# Must be done before any CDAT library is called.
if 'UVCDAT_ANONYMOUS_LOG' not in os.environ:
    os.environ['UVCDAT_ANONYMOUS_LOG'] = 'no'
# Needed for when using hdf5 >= 1.10.0,
# without this, errors are thrown on Edison compute nodes.
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

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
from acme_diags.driver.utils import get_set_name, SET_NAMES


def _get_default_diags(set_name, run_type):
    """
    Returns the path for the default diags for plotset set_name.
    These are different depending on the run_type.
    """
    set_name = get_set_name(set_name)

    folder = '{}'.format(set_name)
    fnm = '{}_{}.cfg'.format(set_name, run_type)
    pth = os.path.join(acme_diags.INSTALL_PATH, folder, fnm)

    print('Using {} for {}.'.format(pth, set_name))
    if not os.path.exists(pth):
        raise RuntimeError(
            "Plotting via set '{}' not supported, file {} not installed".format(set_name, fnm))
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


def get_parameters(parser=ACMEParser()):
    """
    Get the parameters from the parser.
    """
    args = parser.view_args()

    # There weren't any arguments defined
    if not any(getattr(args, arg) for arg in vars(args)):
        parser.print_help()
        sys.exit()

    if args.parameters and not args.other_parameters:  # -p only
        cmdline_parameter = parser.get_cmdline_parameters()
        # If just a -p with no command line parameters, check the py for errors.
        # Otherwise don't check for errors, the command line args might have some missing ones.
        check_values = True if not cmdline_parameter else False
        # default_vars needs to be True in original_parameter b/c 
        original_parameter = parser.get_orig_parameters(check_values=check_values)
        #if not hasattr(original_parameter, 'sets'):
        #    # Since sets is a selector, it must be in there.

        # WITHOUT THIS, THERE ARE A BUNCH OF RESULTS.
        # It's need to get the default selector and everything.
        # But when this is done, it even gets the default EVERYTHING, case_id, etc.
        # And all of them overwrite all of the stuff in the CFG.
        # THIS CANNOT HAPPEN.
        # parser.add_default_values(original_parameter, default_vars=True)
        # We need some way to get the default values in original_parameter before it gets run through select().

        # Maybe the selector should use something other than the *py to look for the values??
        # Maybe create an original_param with all of the defaults with user options on top of it.
        # And just use this to select the values.

        # Load the default cfg files.
        default_diags_paths = []
        # TODO: Fix how the run_type is retrieved.
        run_type = 'model_vs_obs' # original_parameter.run_type
        for set_name in SET_NAMES:
            default_diags_paths.append(_get_default_diags(set_name, run_type))

        other_parameters = parser.get_other_parameters(files_to_open=default_diags_paths)

        # NEW INFO:
        # Need to add selectors in the original_parameter before combine_parameter()
        # Because it's used as vars_to_ignore.
        '''
        parameters = parser.get_parameters(cmdline_parameters=cmdline_parameter,
            orig_parameters=original_parameter, other_parameters=other_parameters)
        '''
        parameters = parser.get_parameters(other_parameters=other_parameters, cmd_default_vars=False)

    elif not args.parameters and args.other_parameters:  # -d only
        cmdline_parameter = parser.get_cmdline_parameters()
        other_parameters = parser.get_other_parameters()
        parameters = parser.get_parameters(cmdline_parameters=cmdline_parameter, 
            other_parameters=other_parameters, cmd_default_vars=False)
    else:  # -p and -d, or just command line arguments.
        parameters = parser.get_parameters(cmd_default_vars=False)

    return parameters


def run_diag(parameters):
    """
    For a single set of parameters, run the corresponding diags.
    """
    results = []
    for pset in parameters.sets:
        set_name = get_set_name(pset)

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
    parameters = get_parameters(parser)

    print(len(parameters))
    for p in parameters:
        print(p.selectors, end=' ')
        # print(p.granulate, end=' ')
        print(p.sets, end=' ')
        print(p.variables, end=' ')
        print(p.seasons, end=' ')
        #print(p.case_id)
        print(p.ref_name)

    quit()
    dt = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    for p in parameters:
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
