import os
import sys

# Sets path for "import CIME", etc
CIMEROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "cime"))
sys.path.insert(0, CIMEROOT)
# Sets path for "import provenance"
sys.path.insert(1, os.getcwd())

import pytest

from CIME import utils
from CIME.tests import scripts_regression_tests

os.environ["SRCROOT"] = os.path.join(os.getcwd(), "..", "..")
os.environ["CIME_GLOBAL_WALLTIME"] = "0:05:00"


def pytest_addoption(parser):
    # set addoption as add_argument to use common argument setup
    # pytest's addoption has same signature as add_argument
    setattr(parser, "add_argument", parser.addoption)

    scripts_regression_tests.setup_arguments(parser)

    # verbose and debug flags already exist
    parser.addoption("--silent", action="store_true", help="Disable all logging")


def pytest_configure(config):
    kwargs = vars(config.option)

    utils.configure_logging(kwargs["verbose"], kwargs["debug"], kwargs["silent"])

    scripts_regression_tests.configure_tests(**kwargs)
