"""
Contains the definition of a class to hold information on expected failures for a single test
"""

from CIME.XML.standard_module_setup import *

EXPECTED_FAILURE_COMMENT = "(EXPECTED FAILURE)"
UNEXPECTED_FAILURE_COMMENT_START = "(UNEXPECTED" # There will be some additional text after this, before the end parentheses

class ExpectedFails(object):

    def __init__(self):
        """Initialize an empty ExpectedFails object"""
        self._fails = {}

    def __eq__(self, rhs):
        expect(isinstance(rhs, ExpectedFails), "Wrong type")
        return self._fails == rhs._fails # pylint: disable=protected-access

    def __ne__(self, rhs):
        result = self.__eq__(rhs)
        return not result

    def __repr__(self):
        return repr(self._fails)

    def add_failure(self, phase, expected_status):
        """Add an expected failure to the list"""
        expect(phase not in self._fails, "Phase {} already present in list".format(phase))
        self._fails[phase] = expected_status

    def expected_fails_comment(self, phase, status):
        """Returns a string giving the expected fails comment for this phase and status"""
        if phase not in self._fails:
            return ''

        if self._fails[phase] == status:
            return EXPECTED_FAILURE_COMMENT
        else:
            return "{}: expected {})".format(UNEXPECTED_FAILURE_COMMENT_START,
                                             self._fails[phase])
