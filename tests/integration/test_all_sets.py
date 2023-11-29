import ast
import configparser
import os
import re
import shutil
from typing import List

import pytest

from e3sm_diags.parameter import SET_TO_PARAMETERS
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner
from tests.integration.config import TEST_DATA_DIR
from tests.integration.utils import run_cmd_and_pipe_stderr

# The path to the integration test data, which needs to be downloaded
# prior to running this test file.
MODULE_PATH = os.path.dirname(__file__)
TEST_DATA_PATH = os.path.join(MODULE_PATH, TEST_DATA_DIR)

# The path to the integration test diagnostics .cfg file.
CFG_PATH = os.path.join(MODULE_PATH, "all_sets_modified.cfg")
CFG_PATH = os.path.abspath(CFG_PATH)


class TestAllSets:
    def test_all_sets(self):
        expected_num_diags = 12

        # *_data_path needs to be added b/c the tests runs the diags from a different location
        cmd = (
            f"e3sm_diags_driver.py -d {CFG_PATH} "
            f"--reference_data_path {TEST_DATA_PATH} "
            f"--test_data_path {TEST_DATA_PATH}"
        )

        stderr = run_cmd_and_pipe_stderr(cmd)

        # count the number of pngs in viewer_dir
        results_dir = self._get_results_dir(stderr)
        count = self._count_images(results_dir)

        # -1 is needed because of the E3SM logo in the viewer html
        assert count - 1 == expected_num_diags

        shutil.rmtree(results_dir)  # remove all generated results from the diags

    def _get_results_dir(self, output: List[str]):
        """Given output from e3sm_diags_driver, extract the path to results_dir."""
        for line in output:
            match = re.search("Viewer HTML generated at (.*)viewer.*.html", line)
            if match:
                results_dir = match.group(1)
                return results_dir

        raise RuntimeError("No viewer directory listed in output: {}".format(output))

    def _count_images(self, directory: str):
        """Count the number of images of type file_type in directory"""
        count = 0

        for _, __, files in os.walk(directory):
            for f in files:
                if f.endswith("png"):
                    count += 1

        return count

    @pytest.mark.xfail
    def test_all_sets_directly(self):
        # TODO: This test is meant to replace `test_all_sets`. It should create
        # CoreParameter objects per diagnostic set defined in
        # `all_sets_modified.cfg`. These CoreParameter objects should then be
        # passed to `runner.run_diags()`. The benefit with this approach is
        # that we don't need to run the command using subprocess, which means
        # immediate unit testing feedback rather than waiting to pipe the
        # complete stderr. We can also step through the code for debugging using
        # an interactive console such as VS Code's Python debugger.
        params = self._convert_cfg_to_param_objs()

        for param in params:
            runner.run_diags([param])

    def _convert_cfg_to_param_objs(self) -> List[CoreParameter]:
        """Convert diagnostic cfg entries to parameter objects.

        NOTE: ast.literal_eval is not considered "safe" on untrusted data.
        The reason why it is used is because `configparser.ConfigParser`
        doesn't work well with parsing Python types from strings in
        `.cfg` files, resulting in things such as nested strings or string
        representation of lists. Since we are only calling literal_eval on
        `.cfg` files hosted in this repo, there is minimal risk here.

        Returns
        -------
        List[CoreParameter]
            A list of CoreParameter objects, one for each diagnotic set.
        """
        config = configparser.ConfigParser()
        config.read(CFG_PATH)
        params = []

        for set_name in config.sections():
            param = SET_TO_PARAMETERS[set_name]()

            for option in config.options(set_name):
                val = config.get(set_name, option)
                val = ast.literal_eval(val)

                setattr(param, option, val)

            param.reference_data_path = TEST_DATA_PATH
            param.test_data_path = TEST_DATA_PATH

            params.append(param)

        return params
