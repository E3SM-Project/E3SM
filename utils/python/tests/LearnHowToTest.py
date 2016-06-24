import unittest, sys, os, shutil 
import logging
import subprocess
from subprocess import call



LIB_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(LIB_DIR)
from CIME.utils import run_cmd
import CIME.utils, update_acme_tests, wait_for_tests
import CIME.system_test
from  CIME.system_test import SystemTest
from  CIME.XML.machines import Machines
from  CIME.XML.files import Files
from  CIME.case import Case

SCRIPT_DIR  = CIME.utils.get_scripts_root()
TOOLS_DIR   = os.path.join(SCRIPT_DIR,"Tools")
DART_CONFIG = "DartConfiguration.tcl"
DART_BACKUP = "temp_dart_backup"
MACHINE     = Machines()

# Global settings
SCRIPT_DIR  = CIME.utils.get_scripts_root()

logging.basicConfig( stream=sys.stderr )
logger = logging.getLogger( "LHTT" )
logger.setLevel( logging.DEBUG )
    

def fun(x):
    return x + 1
    
def debug_test() :
    logger.warning("Searching for xmlquery")
    TOOLS_DIR = SCRIPT_DIR + "/Tools"
    xmlquery = TOOLS_DIR + "/xmlquery"
    
    opt = "--listall --valonly"
    out = subprocess.check_output([xmlquery , "--listall" , "--valonly"])
    logger.warning(out)
        

class XMLQUERY(unittest.TestCase):
    
    def setUp(self):
        # Create case directory
        self._testroot = "/tmp/" # MACHINE.get_value("CESMSCRATCHROOT")
        self._testdirs = []
        self._do_teardown = []
        
        logger.debug(self)
        
        testdir = os.path.join(self._testroot, 'scripts_regression_tests.testscripts')
        machine = 'edison'
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        self._testdirs.append(testdir)
        cmd = "%s/create_newcase --case %s --compset X --res f19_g16 --mach %s " % (SCRIPT_DIR, testdir, machine)
        stat, output, errput = run_cmd(cmd, ok_to_fail=True, from_dir=SCRIPT_DIR)
        self.assertEqual(stat, 0, msg="COMMAND '%s' SHOULD HAVE WORKED\noutput:\n%s\n\nerrput:\n%s" % (cmd, output, errput))
        # cmd = "./case.setup "
#         stat, output, errput = run_cmd(cmd, ok_to_fail=True, from_dir=testdir)
#         self.assertEqual(stat, 0, msg="COMMAND '%s' from case directory '%s' SHOULD HAVE WORKED\noutput:\n%s\n\nerrput:\n%s" % (cmd, testdir, output, errput))
#         cmd = "./case.build"
#         stat, output, errput = run_cmd(cmd, ok_to_fail=True, from_dir=testdir)
#         self.assertEqual(stat, 0, msg="COMMAND '%s' from case directory '%s' SHOULD HAVE WORKED\noutput:\n%s\n\nerrput:\n%s" % (cmd, testdir,output, errput))
        self._do_teardown.append(testdir)
        
        
    def tearDown(self):
        pass
        
        
    def test_XMLQUERY(self):
     
        print(SCRIPT_DIR)
        logger.warning(SCRIPT_DIR)
        
        TOOLS_DIR = SCRIPT_DIR + "/Tools"
        xmlquery = TOOLS_DIR + "/xmlquery"
        testdir  = self._testdirs[0]
        
        #call(["ls", "-l"])
        #c=subprocess.check_call([xmlquery])
        
        # Check for environment 
        self.assertTrue(os.path.isdir(SCRIPT_DIR))
        self.assertTrue(os.path.isdir(TOOLS_DIR))
        self.assertTrue(os.path.isfile(xmlquery) )
        
        # Test command line options
        
        options = [ 
            [] , 
            ['--listall'] , 
            ['--listall' , '--valonly'] , 
            #['-value' , ' JOB_QUEUE'] ,
            #['--value' , 'JOB_QUEUE'] ,
            #['-value' , '-attribute' , 'JOB_QUEUE'] ,
            #['-caseroot' , '/Volumes/Blues/ACME/cime/scripts/mycase3/' , '-valonly' ,
            # '-subgroup' , 'case.run' , '-value' ,  'JOB_QUEUE' ] ,
            ['-caseroot' , testdir , '-valonly' ,
             '-subgroup' , 'case.run' , '-value' ,  'JOB_QUEUE' ] ,
        ]
        
        for opt in options :
            
            out = ''
            call = [xmlquery] + opt
            out = subprocess.check_output(call)
            logger.warning( "Testing %s %s" , xmlquery , " ".join(opt))
            logger.debug("Output for %s : %s" , " ".join(opt) , out)
            
            self.assertTrue(len(out) , msg="output exists for:\n " + xmlquery + " " + " ".join(opt) + '\n' + out + 'END' )
            
        # self.assertTrue(0 , msg='check output')
        
        #self.assertFalse(os.path.isdir(SCRIPT_DIR))
        
        
                 




if __name__ == '__main__':
   
    
    #print(fun(10))
    #debug_test()
    unittest.main()        