"""
Interface to the env_test.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from env_base import EnvBase

class SMS(SystemTestsCommon):
    def __init__(self):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCommon.__init__(self)

    def run(self):
        SystemTestsCommon.run(self)

    def report(self):
        SystemTestsCommon.report(self)
