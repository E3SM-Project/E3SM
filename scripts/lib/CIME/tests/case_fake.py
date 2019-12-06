"""
This module contains a fake implementation of the Case class that can be used
for testing the tests.
"""

import os
from copy import deepcopy

class CaseFake(object):
    def __init__(self, case_root, create_case_root=True):
        """
        Initialize a new case object for the given case_root directory.

        Args:
            case_root (str): path to CASEROOT
            create_case_root (bool): If True, creates the directory given by case_root
        """
        self.vars = dict()
        if create_case_root:
            os.makedirs(case_root)
        self.set_value('CASEROOT', case_root)
        casename = os.path.basename(case_root)
        # Typically, CIME_OUTPUT_ROOT is independent of the case. Here,
        # we nest it under CASEROOT so that (1) tests don't interfere
        # with each other; (2) a cleanup that removes CASEROOT will also
        # remove CIME_OUTPUT_ROOT.
        self.set_value('CIME_OUTPUT_ROOT',
                       os.path.join(case_root, 'CIME_OUTPUT_ROOT'))
        self.set_value('CASE', casename)
        self.set_value('CASEBASEID', casename)
        self.set_value('RUN_TYPE', 'startup')
        self.set_exeroot()
        self.set_rundir()

    def get_value(self, item):
        """
        Get the value of the given item

        Returns None if item isn't set for this case

        Args:
            item (str): variable of interest
        """
        return self.vars.get(item)

    def set_value(self, item, value):
        """
        Set the value of the given item to the given value

        Args:
            item (str): variable of interest
            value (any type): new value for item
        """
        self.vars[item] = value

    def copy(self, newcasename, newcaseroot):
        """
        Create and return a copy of self, but with CASE and CASEBASEID set to newcasename,
        CASEROOT set to newcaseroot, and RUNDIR set appropriately.

        Args:
            newcasename (str): new value for CASE
            newcaseroot (str): new value for CASEROOT
        """
        newcase = deepcopy(self)
        newcase.set_value('CASE', newcasename)
        newcase.set_value('CASEBASEID', newcasename)
        newcase.set_value('CASEROOT', newcaseroot)
        newcase.set_exeroot()
        newcase.set_rundir()

        return newcase

    def create_clone(self, newcase, keepexe=False, mach_dir=None, project=None,
                     cime_output_root=None, exeroot=None, rundir=None):
        # Need to disable unused-argument checking: keepexe is needed to match
        # the interface of Case, but is not used in this fake implementation
        #
        # pylint: disable=unused-argument
        """
        Create a clone of the current case. Also creates the CASEROOT directory
        for the clone case (given by newcase).

        Args:
            newcase (str): full path to the new case. This directory should not
                already exist; it will be created
            keepexe (bool, optional): Ignored
            mach_dir (str, optional): Ignored
            project (str, optional): Ignored
            cime_output_root (str, optional): New CIME_OUTPUT_ROOT for the clone
            exeroot (str, optional): New EXEROOT for the clone
            rundir (str, optional): New RUNDIR for the clone

        Returns the clone case object
        """
        newcaseroot = os.path.abspath(newcase)
        newcasename = os.path.basename(newcase)
        os.makedirs(newcaseroot)
        clone = self.copy(newcasename = newcasename, newcaseroot = newcaseroot)
        if cime_output_root is not None:
            clone.set_value('CIME_OUTPUT_ROOT', cime_output_root)
        if exeroot is not None:
            clone.set_value('EXEROOT', exeroot)
        if rundir is not None:
            clone.set_value('RUNDIR', rundir)

        return clone

    def flush(self):
        pass

    def make_rundir(self):
        """
        Make directory given by RUNDIR
        """
        os.makedirs(self.get_value('RUNDIR'))

    def set_exeroot(self):
        """
        Assumes CASEROOT is already set; sets an appropriate EXEROOT
        (nested inside CASEROOT)
        """
        self.set_value('EXEROOT', os.path.join(self.get_value('CASEROOT'), 'bld'))

    def set_rundir(self):
        """
        Assumes CASEROOT is already set; sets an appropriate RUNDIR (nested
        inside CASEROOT)
        """
        self.set_value('RUNDIR', os.path.join(self.get_value('CASEROOT'), 'run'))

    def case_setup(self, clean=False, test_mode=False, reset=False):
        pass

    def load_env(self, reset=False):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *_):
        pass
