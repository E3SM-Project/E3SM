"""
Creates a test suite-specific cs.status file from a template
"""

from CIME.XML.standard_module_setup import *
import CIME.utils
import os
import stat

def create_cs_status(test_root, test_id, extra_args='', filename=None):
    """Create a test suite-specific cs.status file from the template

    Arguments:
    test_root (string): path to test root; the file will be put here. If
        this directory doesn't exist, it is created.
    test_id (string): test id for this test suite. This can contain
        shell wildcards if you want this one cs.status file to work
        across multiple test suites. However, be careful not to make
        this too general: for example, ending this with '*' will pick up
        the *.ref1 directories for ERI and other tests, which is NOT
        what you want.
    extra_args (string): extra arguments to the cs.status command
        (If there are multiple arguments, these should be in a space-delimited string.)
    filename (string): name of the generated cs.status file. If not
        given, this will be built from the test_id.
    """
    python_libs_root = CIME.utils.get_python_libs_root()
    cime_root = CIME.utils.get_cime_root()
    template_file = os.path.join(python_libs_root, "cs.status.template")
    template = open(template_file, "r").read()
    template = template.replace("<PATH>",
                                os.path.join(cime_root,"scripts","Tools")).replace\
                                ("<EXTRA_ARGS>", extra_args).replace\
                                ("<TESTID>", test_id).replace\
                                ("<TESTROOT>", test_root)
    if not os.path.exists(test_root):
        os.makedirs(test_root)
    if filename is None:
        filename = "cs.status.{}".format(test_id)
    cs_status_file = os.path.join(test_root, filename)
    with open(cs_status_file, "w") as fd:
        fd.write(template)
    os.chmod(cs_status_file, os.stat(cs_status_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)
