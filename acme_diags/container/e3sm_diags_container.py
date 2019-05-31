#!/usr/bin/env python
"""
A standalone script that prepares a set of input to e3sm_diags
to be ran as a container, and then runs e3sm_diags as a container.
"""
import os
import sys
import importlib
import argparse
import subprocess

# Change these commands if needed.
SHIFTER_COMMAND = 'shifter --volume=$REFERENCE_DATA_PATH:/reference_data_path'
SHIFTER_COMMAND += ' --volume=$TEST_DATA_PATH:/test_data_path --volume=$RESULTS_DIR:/results_dir'
SHIFTER_COMMAND += ' --image=docker:e3sm/e3sm_diags:{}'
# Shifter doesn't use the entrypoint defined in the Dockerfile, so we need to specify what command to use.
SHIFTER_COMMAND += ' -- e3sm_diags'

DOCKER_COMMAND = 'docker run --mount type=bind,source=$REFERENCE_DATA_PATH,target=/reference_data_path' 
DOCKER_COMMAND += ' --mount type=bind,source=$TEST_DATA_PATH,target=/test_data_path'
DOCKER_COMMAND += ' --mount type=bind,source=$RESULTS_DIR,target=/results_dir'
# Docker needs the cwd mounted as well, otherwise the input parameter files will not be found.
DOCKER_COMMAND += ' --mount type=bind,source="$(pwd)",target=/e3sm_diags_container_cwd'
DOCKER_COMMAND += ' e3sm/e3sm_diags:{}'

SINGULARITY_COMMAND = 'singularity run -B $REFERENCE_DATA_PATH:/reference_data_path,'
SINGULARITY_COMMAND += '$TEST_DATA_PATH:/test_data_path,$RESULTS_DIR:/results_dir '
# It seems like Singularity doesn't automatically look in $SINGULARITY_PULLFOLDER.
SINGULARITY_COMMAND += '$SINGULARITY_PULLFOLDER/e3sm_diags-{}.simg'

def run_cmd(cmd):
    """
    Given a command, run it.
    """
    print('Using the command: {}'.format(cmd))
    # p = subprocess.Popen(cmd, shell=True)
    p = subprocess.Popen(cmd, shell=True).wait()
    # This doesn't work: p = subprocess.Popen(cmd.split(), shell=True)


def run_container(args):
    e3sm_diags_args = get_user_args_for_e3sm_diags()

    # Append the e3sm_diags arguments to the container command.
    if args.shifter:
        cmd = SHIFTER_COMMAND.format(args.container_version)
        cmd += ' ' + ' '.join(e3sm_diags_args)
        run_cmd(cmd)
    elif args.singularity:
        cmd = SINGULARITY_COMMAND.format(args.container_version)
        cmd += ' ' + ' '.join(e3sm_diags_args)
        run_cmd(cmd)
    elif args.docker:
        cmd = DOCKER_COMMAND.format(args.container_version)
        cmd += ' ' + ' '.join(e3sm_diags_args)
        run_cmd(cmd)
    else:
        msg = 'Invalid container runtime option. Please choose '
        msg += 'from "--shifter", "--singularity", or "--docker".'
        raise RuntimeError(msg)


def get_parameter_from_file(path, param):
    """
    From the Python parameters file located at path,
    get the value of the param attribute from it.
    """
    if not os.path.exists(path) or not os.path.isfile(path):
        raise IOError('The parameter file passed in does not exist.')

    path_to_module = os.path.split(path)[0]
    module_name = os.path.split(path)[1]
    if '.' in module_name:
        module_name = module_name.split('.')[0]

    sys.path.insert(0, path_to_module)
    module = importlib.import_module(module_name)

    if not hasattr(module, param):
        msg = 'Error: {} was not defined in the command'.format(param)
        msg += ' line nor in the Python file.'
        raise AttributeError(msg)
    else:
        return getattr(module, param)
    

def set_env_vars(args):
    """
    From the user's input, set the right environmental
    variables to run the container.
    """
    # Get the parameters from the command line.
    reference_data_path = args.reference_data_path
    test_data_path = args.test_data_path
    results_dir = args.results_dir

    # If they are empty, try to get them from the Python file.
    param_file = args.parameters
    if not reference_data_path:
        reference_data_path = get_parameter_from_file(param_file, 'reference_data_path')
    if not test_data_path:
        test_data_path = get_parameter_from_file(param_file, 'test_data_path')
    if not results_dir:
        results_dir = get_parameter_from_file(param_file, 'results_dir')

    # Need to make sure paths are valid before actually setting the environmental variables.
    if not os.path.exists(reference_data_path):
        msg = '{} does not exist.'.format(reference_data_path)
        raise IOError(msg)
    if not os.path.exists(test_data_path):
        msg = '{} does not exist.'.format(test_data_path)
        raise IOError(msg)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir, 0o755)
    
    # Make the paths absolute.
    # Docker needs this.
    reference_data_path = os.path.abspath(reference_data_path)
    test_data_path = os.path.abspath(test_data_path)
    results_dir = os.path.abspath(results_dir)

    # Then set them as the environmental variables.
    os.environ['REFERENCE_DATA_PATH'] = reference_data_path
    os.environ['TEST_DATA_PATH'] = test_data_path
    os.environ['RESULTS_DIR'] = results_dir


def get_user_args_for_e3sm_diags():
    """
    Extract the correct passed in arguments from this script that are needed for e3sm_diags.
    """
    params_to_ignore = ['e3sm_diags_container', 'python', 'e3sm_diags_container.py',
                        '--shifter', '--singularity', '--docker', '--container_version']

    e3sm_diags_args = []

    i = 0
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg not in params_to_ignore:
            # Commands must be in '' so it works correctly. 
            e3sm_diags_args.append("'{}'".format(arg))
            i += 1
        elif arg == '--container_version':
            # Skip the 'container_version' arg, as well as the version specified.
            # That's why we skip by two.
            i += 2
        else:
            i += 1

    return e3sm_diags_args


# For this preprocessing step, these are the only parameters we care about.
parser = argparse.ArgumentParser()
parser.add_argument(
    '-p', '--parameters',
    dest='parameters',
    default='',
    help='Path to the user-defined parameter file.',
    required=False)
parser.add_argument(
    '--reference_data_path',
    dest='reference_data_path',
    help='Path for the reference data.',
    required=False)
parser.add_argument(
    '--test_data_path',
    dest='test_data_path',
    help='Path for the test data.',
    required=False)
parser.add_argument(
    '--results_dir',
    dest='results_dir',
    help='Path of where to save the results.',
    required=False)
parser.add_argument(
    '--container_version',
    dest='container_version',
    help='Which version of the container to use.',
    default='latest',
    required=False)
# A separate group of arguments for the container runtimes.
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    '--shifter',
    dest='shifter',
    help='Run the diags in a container using shifter.',
    action='store_const',
    const=True,
    required=False)
group.add_argument(
    '--singularity',
    dest='singularity',
    help='Run the diags in a container using singularity.',
    action='store_const',
    const=True,
    required=False)
group.add_argument(
    '--docker',
    dest='docker',
    help='Run the diags in a container using docker.',
    action='store_const',
    const=True,
    required=False)

args, unknown = parser.parse_known_args()
set_env_vars(args)
run_container(args)
