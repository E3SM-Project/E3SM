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
        self.set_value('CASE', casename)
        self.set_value('CASEBASEID', casename)
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
        newcase.set_rundir()

        return newcase

    def create_clone(self, newcase, keepexe=False):
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

        Returns the clone case object
        """
        newcaseroot = os.path.abspath(newcase)
        newcasename = os.path.basename(newcase)
        os.makedirs(newcaseroot)
        clone = self.copy(newcasename = newcasename, newcaseroot = newcaseroot)

        return clone

    def flush(self):
        pass

    def make_rundir(self):
        """
        Make directory given by RUNDIR
        """
        os.makedirs(self.get_value('RUNDIR'))

    def set_rundir(self):
        """
        Assumes CASEROOT is already set; sets an appropriate RUNDIR (nested
        inside CASEROOT)
        """
        self.set_value('RUNDIR', os.path.join(self.get_value('CASEROOT'), 'run'))

