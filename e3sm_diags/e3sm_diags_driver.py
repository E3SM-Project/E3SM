#!/usr/bin/env python
# The above line is needed for `test_all_sets.test_all_sets_mpl`.
# Otherwise, OSError: [Errno 8] Exec format error: 'e3sm_diags_driver.py'.
import os
import subprocess
import sys
import traceback
from typing import Dict, List, Tuple

import dask
import dask.bag as db

import e3sm_diags
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parser import SET_TO_PARSER
from e3sm_diags.parser.core_parser import CoreParser
from e3sm_diags.viewer.main import create_viewer

logger = custom_logger(__name__)


def get_default_diags_path(set_name, run_type, print_path=True):
    """
    Returns the path for the default diags for plotset set_name.
    These are different depending on the run_type.
    """
    folder = "{}".format(set_name)
    fnm = "{}_{}.cfg".format(set_name, run_type)
    pth = os.path.join(e3sm_diags.INSTALL_PATH, folder, fnm)

    if print_path:
        logger.info("Using {} for {}.".format(pth, set_name))
    if not os.path.exists(pth):
        raise RuntimeError(
            "Plotting via set '{}' not supported, file {} not installed".format(
                set_name, fnm
            )
        )
    return pth


def _save_env_yml(results_dir):
    """
    Save the yml to recreate the environment in results_dir.
    """
    cmd = "conda env export"
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()

    if err:
        logger.exception("Error when creating env yml file: ")
        logger.exception(err)
    else:
        fnm = os.path.join(results_dir, "environment.yml")
        with open(fnm, "w") as f:
            f.write(output.decode("utf-8"))
        logger.info("Saved environment yml file to: {}".format(fnm))


def _save_parameter_files(results_dir, parser):
    """
    Save the command line arguments used, and any py or cfg files.
    """
    cmd_used = " ".join(sys.argv)
    fnm = os.path.join(results_dir, "cmd_used.txt")
    with open(fnm, "w") as f:
        f.write(cmd_used)
    logger.info("Saved command used to: {}".format(fnm))

    args = parser.view_args()

    if hasattr(args, "parameters") and args.parameters:
        fnm = args.parameters
        if not os.path.isfile(fnm):
            logger.warning("File does not exist: {}".format(fnm))
        else:
            with open(fnm, "r") as f:
                contents = "".join(f.readlines())
            # Remove any path, just keep the filename.
            new_fnm = fnm.split("/")[-1]
            new_fnm = os.path.join(results_dir, new_fnm)
            with open(new_fnm, "w") as f:
                f.write(contents)
            logger.info("Saved py file to: {}".format(new_fnm))

    if hasattr(args, "other_parameters") and args.other_parameters:
        fnm = args.other_parameters[0]
        if not os.path.isfile(fnm):
            logger.warning("File does not exist: {}".format(fnm))
        else:
            with open(fnm, "r") as f:
                contents = "".join(f.readlines())
            # Remove any path, just keep the filename.
            new_fnm = fnm.split("/")[-1]
            new_fnm = os.path.join(results_dir, new_fnm)
            with open(new_fnm, "w") as f:
                f.write(contents)
            logger.info("Saved cfg file to: {}".format(new_fnm))


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
    py_files = [f for f in sys.argv if f.endswith(".py")]
    # User didn't pass in a Python file, so they maybe ran:
    #    e3sm_diags -d diags.cfg
    if not py_files:
        return

    fnm = py_files[-1]

    if not os.path.isfile(fnm):
        logger.warning("File does not exist: {}".format(fnm))
        return

    with open(fnm, "r") as f:
        contents = "".join(f.readlines())
    # Remove any path, just keep the filename.
    new_fnm = fnm.split("/")[-1]
    new_fnm = os.path.join(results_dir, new_fnm)
    with open(new_fnm, "w") as f:
        f.write(contents)
    logger.info("Saved Python script to: {}".format(new_fnm))


def save_provenance(results_dir, parser):
    """
    Store the provenance in results_dir.
    """
    results_dir = os.path.join(results_dir, "prov")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir, 0o755)

    # Create a PHP file to list the contents of the prov dir.
    php_path = os.path.join(results_dir, "index.php")
    with open(php_path, "w") as f:
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
    except Exception:
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
    args, _ = parser_for_args.parser.parse_known_args()

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
        parameters = parser.get_parameters(
            cmd_default_vars=False, argparse_vals_only=False
        )

    # The below two clauses are for the legacy way to
    # run this software, pre v2.0.0.
    # Ex: e3sm_diags -p params.py -d diags.cfg
    elif args.parameters and not args.other_parameters:  # -p only
        original_parameter = parser.get_orig_parameters(argparse_vals_only=False)

        # Load the default cfg files.
        run_type = getattr(original_parameter, "run_type", "model_vs_obs")
        default_diags_paths = [
            get_default_diags_path(set_name, run_type)
            for set_name in CoreParameter().sets
        ]

        other_parameters = parser.get_cfg_parameters(
            files_to_open=default_diags_paths, argparse_vals_only=False
        )

        parameters = parser.get_parameters(
            orig_parameters=original_parameter,
            other_parameters=other_parameters,
            cmd_default_vars=False,
            argparse_vals_only=False,
        )

    else:
        parameters = parser.get_parameters(
            cmd_default_vars=False, argparse_vals_only=False
        )

    parser.check_values_of_params(parameters)

    if not parameters:
        msg = "No parameters were able to be created. Please check your .py "
        msg += "file, and any .cfg files or command line args you're using."
        raise RuntimeError(msg)

    return parameters


def create_parameter_dict(parameters):
    d: Dict[type, int] = dict()
    for parameter in parameters:
        t = type(parameter)
        if t in d.keys():
            d[t] += 1
        else:
            d[t] = 1
    return d


def _run_serially(parameters: List[CoreParameter]) -> List[CoreParameter]:
    """Run diagnostics with the parameters serially.

    Parameters
    ----------
    parameters : List[CoreParameter]
        The list of CoreParameter objects to run diagnostics on.

    Returns
    -------
    List[CoreParameter]
        The list of CoreParameter objects with results from the diagnostic run.
    """
    # A nested list of lists, where a sub-list represents the results of
    # the sets related to the CoreParameter object.
    nested_results: List[List[CoreParameter]] = []

    for parameter in parameters:
        nested_results.append(parameter._run_diag())

    # `results` becomes a list of lists of parameters so it needs to be
    # collapsed a level.
    collapsed_results = _collapse_results(nested_results)

    return collapsed_results


def _run_with_dask(parameters: List[CoreParameter]) -> List[CoreParameter]:
    """Run diagnostics with the parameters in parallel using Dask.

    This function passes ``run_diag`` to ``dask.bag.map``, which gets executed
    in parallel with ``.compute``.

    The first CoreParameter object's `num_workers` attribute is used to set
    the number of workers for ``.compute``.

    Parameters
    ----------
    parameters : List[CoreParameter]
        The list of CoreParameter objects to run diagnostics on.

    Returns
    -------
    List[CoreParameter]
        The list of CoreParameter objects with results from the diagnostic run.

    Notes
    -----
    https://docs.dask.org/en/stable/generated/dask.bag.map.html
    https://docs.dask.org/en/stable/generated/dask.dataframe.DataFrame.compute.html
    """
    bag = db.from_sequence(parameters)
    config = {"scheduler": "processes", "multiprocessing.context": "fork"}

    num_workers = getattr(parameters[0], "num_workers", None)
    if num_workers is None:
        raise ValueError(
            "The `num_workers` attribute is required for multiprocessing but it is not "
            "defined on the CoreParameter object. Set this attribute and try running "
            "again."
        )

    with dask.config.set(config):
        results = bag.map(CoreParameter._run_diag).compute(num_workers=num_workers)

    # `results` becomes a list of lists of parameters so it needs to be
    # collapsed a level.
    collapsed_results = _collapse_results(results)

    return collapsed_results


def _collapse_results(parameters: List[List[CoreParameter]]) -> List[CoreParameter]:
    """Collapses the results of diagnostic runs by one list level.

    Parameters
    ----------
    parameters : List[List[CoreParameter]]
        A list of lists of CoreParameter objects with results from the
        diagnostic run.

    Returns
    -------
    List[CoreParameter]
        A list of CoreParameter objects with results from the diagnostic run.
    """
    output_parameters = []

    for p1 in parameters:
        if isinstance(p1, list):
            for p2 in p1:
                output_parameters.append(p2)
        else:
            output_parameters.append(p1)

    return output_parameters


def main(parameters=[]) -> List[CoreParameter]:
    # Get the diagnostic run parameters
    # ---------------------------------
    parser = CoreParser()

    # If no parameters are passed, use the parser args as defaults. Otherwise,
    # create the dictionary of expected parameters.
    if len(parameters) == 0:
        parameters = get_parameters(parser)

    expected_parameters = create_parameter_dict(parameters)

    if not os.path.exists(parameters[0].results_dir):
        os.makedirs(parameters[0].results_dir, 0o755)
    if not parameters[0].no_viewer:  # Only save provenance for full runs.
        save_provenance(parameters[0].results_dir, parser)

    # Perform the diagnostic run
    # --------------------------
    if parameters[0].multiprocessing:
        parameters_results = _run_with_dask(parameters)
    else:
        parameters_results = _run_serially(parameters)

    # Generate the viewer outputs using results
    # -----------------------------------------
    if not parameters_results:
        logger.warning(
            "There was not a single valid diagnostics run, no viewer created."
        )
    else:
        # If you get `AttributeError: 'NoneType' object has no attribute
        # `'no_viewer'` on this line then `run_diag` likely returns `None`.
        if parameters_results[0].no_viewer:
            logger.info("Viewer not created because the no_viewer parameter is True.")
        else:
            path = os.path.join(parameters_results[0].results_dir, "viewer")
            if not os.path.exists(path):
                os.makedirs(path)

            index_path = create_viewer(path, parameters_results)
            logger.info("Viewer HTML generated at {}".format(index_path))

    # Validate actual and expected parameters are aligned
    # ---------------------------------------------------
    actual_parameters = create_parameter_dict(parameters_results)
    if parameters_results[0].fail_on_incomplete and (
        actual_parameters != expected_parameters
    ):
        d: Dict[type, Tuple[int, int]] = dict()

        # Loop through all expected parameter types.
        for t in expected_parameters.keys():
            d[t] = (actual_parameters[t], expected_parameters[t])

        message = (
            "Not all parameters completed successfully. Check output above for "
            "errors/exceptions. The following dictionary maps parameter types to their "
            f"actual and expected numbers: {d}"
        )
        raise Exception(message)

    return parameters_results


if __name__ == "__main__":
    main()
