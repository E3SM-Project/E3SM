"""
This is the "testscripts" module.
"""
import unittest

#from standard_script_setup import *

    
# Standard setup
import sys, os, logging, doctest, argparse, logging.config
import __main__ as main

# Set environment and import cime modules
_CIMEROOT = os.environ.get("CIMEROOT")
if(_CIMEROOT is None):
    _CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..")
    os.environ["CIMEROOT"] = _CIMEROOT
_LIB_DIR = os.path.join(_CIMEROOT, "utils", "python")
sys.path.append(_LIB_DIR)


import CIME.utils
CIME.utils.check_minimum_python_version(2, 7)
CIME.utils.stop_buffering_output()
# set up logging to file
logging.basicConfig(level=logging.INFO)



from CIME.utils import expect, run_cmd
from CIME.case import Case

import argparse, doctest, shutil, glob
import textwrap



class TestCaseModule(unittest.TestCase):
    
    def test_init_case(self):
        
        caseroots = [  os.getcwd() , "path/to/test/data" , None , '' ]
        for d in caseroots :
            logging.info("Testing creating case with caseroot: %s" ,d)
            case = None
            try:
                case    = Case(d)
            except:
                logging.info("Can't create case with %s" , d)
                
            self.assertIsInstance(case,Case)
    
    
    
    
    



if __name__ == "__main__":
    unittest.main()