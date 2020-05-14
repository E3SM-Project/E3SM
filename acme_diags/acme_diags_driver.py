#!/usr/bin/env python
from __future__ import print_function

import os
# Must be done before any CDAT library is called.
os.environ['UVCDAT_ANONYMOUS_LOG'] = 'no'
os.environ['CDAT_ANONYMOUS_LOG'] = 'no'
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
import cdms2.tvariable 
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.parser import SET_TO_PARSER
from acme_diags.parser.core_parser import CoreParser
from acme_diags.viewer.main import create_viewer
from acme_diags.driver import utils
from acme_diags import container


# turn off MPI in cdms2 -- not currently supported by e3sm_diags
cdms2.tvariable.HAVE_MPI = False


def get_default_diags_path(set_name, run_type, print_path=True):
    """
    Returns the path for the default diags for plotset set_name.
    These are different depending on the run_type.
    """
    folder = '{}'.format(set_name)
    fnm = '{}_{}.cfg'.format(set_name, run_type)
    pth = os.path.join(acme_diags.INSTALL_PATH, folder, fnm)

    if print_path:
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
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
        if container.is_container():
            f.write('# e3sm_diags was ran in a container.\n')
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
            # Remove any path, just keep the filename.
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
            # Remove any path, just keep the filename.
            new_fnm = fnm.split('/')[-1]
            new_fnm = os.path.join(results_dir, new_fnm)
            with open(new_fnm, 'w') as f:
                f.write(contents)
            print('Saved cfg file to: {}'.format(new_fnm))


def _save_python_script(results_dir, parser):
    """
    When using a Python script to run the
    diags via the API, dump a copy of the script.
    """
    args = parser.view_args()
    # If running the legacy way, there's
    # nothing to be saved.
    if args.parameters:
        return

    # Get the last argument that has .py in it.
    # The reason we're getting the last .py file
    # is for this case when running in a container:
    #     python e3sm_diag_container.py run_diags.py
    py_files = [f for f in sys.argv if f.endswith('.py')]
    # User didn't pass in a Python file, so they maybe ran:
    #    e3sm_diags -d diags.cfg
    if not py_files:
        return

    fnm = py_files[-1]

    if not os.path.isfile(fnm):
        print('File does not exist: {}'.format(fnm))
        return

    with open(fnm, 'r') as f:
        contents = ''.join(f.readlines())
    # Remove any path, just keep the filename.
    new_fnm = fnm.split('/')[-1]
    new_fnm = os.path.join(results_dir, new_fnm)
    with open(new_fnm, 'w') as f:
        f.write(contents)
    print('Saved Python script to: {}'.format(new_fnm))


def save_provenance(results_dir, parser):
    """
    Store the provenance in results_dir.
    """
    results_dir = os.path.join(results_dir, 'prov')
    if not os.path.exists(results_dir):
        os.makedirs(results_dir, 0o755)

    # Create a PHP file to list the contents of the prov dir.
    php_path = os.path.join(results_dir, 'index.php')
    with open(php_path, 'w') as f:
        contents = """
        <?php
        # Taken from:
        # https://stackoverflow.com/questions/3785055/how-can-i-create-a-simple-index-html-file-which-lists-all-files-directories
        $path = ".";
        $dh = opendir($path);
        $i=1;
        while (($file = readdir($dh)) !== false) {
            if($file != "." && $file != ".." && $file != "index.php" && $file != ".htaccess" && $file != "error_log" && $file != "cgi-bin") {
                echo "<a href='$path/$file'>$file</a><br /><br />";
                $i++;
            }
        }
        closedir($dh);
        ?>
        """
        f.write(contents)
    try:
        _save_env_yml(results_dir)
    except:
        traceback.print_exc()

    _save_parameter_files(results_dir, parser)

    _save_python_script(results_dir, parser)


def get_parameters(parser=CoreParser()):
    """
    Get the parameters from the parser.
    """
    # A separate parser to just get the args used.
    # The reason it's a separate object than `parser`
    # is so we can parse the known args.
    parser_for_args = CoreParser()
    # The unknown args are _.
    # These are any set-specific args that aren't needed
    # for now, we just want to know what args are used.
    args, _ = parser_for_args.parse_known_args()

    # Below is the legacy way to run this software, pre v2.0.0.
    # There weren't any arguments defined.
    if not any(getattr(args, arg) for arg in vars(args)):
        parser.print_help()
        sys.exit()

    # For when a user runs the software with commands like:
    #    e3sm_diags lat_lon [the other parameters]
    # This use-case is usually ran when the provenance
    # command is copied and pasted from the viewers.
    if args.set_name in SET_TO_PARSER:
       parser = SET_TO_PARSER[args.set_name]()
       parameters = parser.get_parameters(cmd_default_vars=False, argparse_vals_only=False) 

    # The below two clauses are for the legacy way to
    # run this software, pre v2.0.0.
    # Ex: e3sm_diags -p params.py -d diags.cfg
    elif args.parameters and not args.other_parameters:  # -p only
        original_parameter = parser.get_orig_parameters(argparse_vals_only=False)

        # Load the default cfg files.
        run_type = getattr(original_parameter, 'run_type', 'model_vs_obs')
        default_diags_paths = [get_default_diags_path(set_name, run_type) for set_name in CoreParameter().sets]

        other_parameters = parser.get_other_parameters(files_to_open=default_diags_paths, argparse_vals_only=False)

        parameters = parser.get_parameters(orig_parameters=original_parameter,
            other_parameters=other_parameters, cmd_default_vars=False, argparse_vals_only=False)

    else:
        parameters = parser.get_parameters(cmd_default_vars=False, argparse_vals_only=False)

    parser.check_values_of_params(parameters)

    if not parameters:
        msg = 'No parameters were able to be created. Please check your .py '
        msg += 'file, and any .cfg files or command line args you\'re using.'
        raise RuntimeError(msg)

    return parameters


def run_diag(parameters):
    """
    For a single set of parameters, run the corresponding diags.
    """
    results = []
    for set_name in parameters.sets:

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


def main(parameters=[]):
    parser = CoreParser()
    if not parameters:
        parameters = get_parameters(parser)

    # each case id (aka, variable) has a set of parameters specified.
    # print out the parameters for each variable for trouble shootting
    #for p in parameters:
    #    attrs = vars(p)
    #    print (', '.join("%s: %s" % item for item in attrs.items()))

    if not os.path.exists(parameters[0].results_dir):
        os.makedirs(parameters[0].results_dir, 0o755)
    if not parameters[0].no_viewer:  # Only save provenance for full runs.
        save_provenance(parameters[0].results_dir, parser)

    if container.is_container():
        print('Running e3sm_diags in a container.')
        # Modify the parmeters so that it runs in
        # the container as if it usually runs.
        for p in parameters:
            container.containerize_parameter(p)

    if parameters[0].multiprocessing:
        parameters = cdp.cdp_run.multiprocess(run_diag, parameters)
    elif parameters[0].distributed:
        parameters = cdp.cdp_run.distribute(run_diag, parameters)
    else:
        parameters = cdp.cdp_run.serial(run_diag, parameters)

    parameters = _collapse_results(parameters)

    if container.is_container():
        for p in parameters:
            container.decontainerize_parameter(p)

    if not parameters:
        print('There was not a single valid diagnostics run, no viewer created.')
    else:
        # If you get `AttributeError: 'NoneType' object has no attribute 'no_viewer'` on this line
        # then `run_diag` likely returns `None`.
        if parameters[0].no_viewer:
            print('Viewer not created because the no_viewer parameter is True.')
        else:
            path = os.path.join(parameters[0].results_dir, 'viewer')
            if not os.path.exists(path):
                os.makedirs(path)
            
            index_path = create_viewer(path, parameters)
            print('Viewer HTML generated at {}'.format(index_path))

if __name__ == '__main__':
    main()
