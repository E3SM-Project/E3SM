"""
Common functions used by cime python scripts
"""
import logging
import sys
import os

_MODEL = None
_CIMEROOT = None

def expect(condition, error_msg):
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.

    >>> expect(True, "error1")
    >>> expect(False, "error2")
    Traceback (most recent call last):
        ...
    SystemExit: FAIL: error2
    """
    if (not condition):
        raise SystemExit('ERROR: '+error_msg)
    

def get_python_libs_location_within_cime():
    """
    From within CIME, return subdirectory of python libraries
    """
    return os.path.join("utils", "python")

def get_cime_root():
    """
    Return the absolute path to the root of CIME that contains this script

    >>> os.path.isdir(os.path.join(get_cime_root(), get_acme_scripts_location_within_cime()))
    True
    """
    global _CIMEROOT
    if (_CIMEROOT is None):
        try:
            _CIMEROOT = os.environ["CIMEROOT"]
        except KeyError:
            script_absdir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
            assert script_absdir.endswith(get_python_libs_location_within_cime()), script_absdir
            _CIMEROOT = os.path.abspath(os.path.join(script_absdir,"..",".."))
    logging.info( "CIMEROOT is " + _CIMEROOT)
    return _CIMEROOT


def set_model(model):
    global _MODEL
    _MODEL = model


def get_model():
    global _MODEL
    if (_MODEL is None):
        try:
            _MODEL = os.environ["CIME_MODEL"]
        except KeyError:
            modelroot = os.path.join(get_cime_root(), "cime_config")
            models = os.listdir(modelroot)
            msg = "Environment variable CIME_MODEL must be set to one of: "
            msg += ", ".join([model for model in models if os.path.isdir(os.path.join(modelroot,model)) and model != "xml_schemas"])
            expect(False, msg)

    return _MODEL


