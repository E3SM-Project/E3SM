#!/usr/bin/env python

"""
The build and test structure (BATS) module. BATS is primarily intended to allow
users and developers of CISM to quickly generate a set of regression tests for
use with the Land Ice Validation and Verification toolkit (LIVVkit). 
"""

import os
import sys
import time
import argparse
import subprocess

from util import paths
from util import dicts
from util import runnit

# get the tests to run! See util/dicts.py for details.
test_dict = dicts.test_dict

# setup our input argument parser
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

def unsigned_int(x):
    """
    Small helper function so argparse will understand unsigned integers.
    """
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

parser.add_argument('-p','--platform', default='linux',  
        help="Your computer platform.")
parser.add_argument('-c','--compiler', default='gnu', 
        help="The compiler.\n")

parser.add_argument('-i','--cism-dir', default=os.pardir+os.sep+os.pardir,
        help="Location of the CISM source code.")
parser.add_argument('-b','--build-dir',default=os.getcwd()+os.sep+'build',
        help="Location to build CISM.")
parser.add_argument('-j', type=unsigned_int, default='8',
        help="Number of processors to use when making CISM.")
parser.add_argument('-s','--skip-build', action='store_true', 
        help="Skip build and look for the CISM driver in BUILD_DIR.")

parser.add_argument('-o','--out-dir', default='reg_test',
        help="The location of the directoy to output the test data.")
parser.add_argument('-f','--force', action='store_true',
        help="Supress any warning about possibly overwriting test data.")
parser.add_argument('--timing', action='store_true',
        help="Run the timing test. This is needed for creating a new benchmark dataset. This is needed for creating a new benchmark dataset. Selecting this option will force the --performance option to true as well.")

#NOTE: These last two are just for personal desktop/laptop type machines. 
#      --performance is always turned on for HPC systems, and tests are
#      only run on personal machines.
parser.add_argument('--performance', action='store_true',
        help="Run the performance tests. Will be forced on for HPC platforms or if the --timing option is specified.")
parser.add_argument('--sleep', type=unsigned_int, default=30,
        help="Number of seconds to sleep between checks on running tests. Only used on personal machines, when HPC platforms are detected, BATS does not run the tests and just creates the batch jobs.")

#TODO: Use these options to turn on and off dycores/trilinos in build (is there already),
#      And then update the [dycore, which_ho_sparse, which_ho_nonlinear] options
#      in the config files (still todo -- option will need to be added to the 
#      test files run*.py.
#           --dycore option is used in util.paths
#           --library option is used here and in util.paths
#           Uncomment them when turning back on.
#parser.add_argument('--dycore', 
#        help="The dycore used. May be a '-' separated list for special cases. \n"
#            +"   Example: bisicles \n"
#            +"   Special: bisicles-petsc-python")
#parser.add_argument('--library',
#        help="Use an external solver library like trilinos.")


# sets up the default options in args. Will be overwritten when running in
# script mode, and can be overwritten by calling the parse function.
global args
args = parser.parse_args([])

def parse(arg_list):
    global args
    args = parser.parse_args(arg_list)
    
def main():
    global args
    # used to modify timing file names
    args.tmod = None 

    # setup the needed paths
    args = paths.make_absolute(args)
    cmake_dir, cmake_file = paths.cmake(args)
    paths.mkdir_p(args.build_dir)
    cism_driver = args.build_dir+os.sep+'cism_driver'+os.sep+'cism_driver'

    # always run performance tests on HPC systems.
    if args.platform.lower() in dicts.hpc_dict.keys():
        isHPC = True
        args.performance = True
    else:
        isHPC = False

    # always run performance tests if timing runs selected.
    if args.timing == True:
        args.performance = True

    if not args.skip_build:
        print("\nPreparing to build CISM")
        print(  "=======================")
          
        print("Build options:")
        print("   Platform: "+args.platform)
        print("   Compiler: "+args.compiler)
        print("-----------------------")
        print("cmake directory: "+cmake_dir)
        print("cmake file: "+cmake_file)
        
        print("\nBuilding CISM")
        print(  "=============\n")

        #TODO: turn on args.library option.
        #if args.library and args.library.lower() == 'trilinos':
        #    trilinos_string = "CISM_USE_TRILINOS=ON"
        #else:
        #    trilinos_string = "CISM_USE_TRILINOS=OFF"

        prep_commands = ["cd "+args.build_dir,
                         #TODO: turn on args.library option.
                         #"export "+trilinos_string,
                         "source "+cmake_dir+os.sep+cmake_file+" "+args.cism_dir,
                         "make -j "+str(args.j),
                         "exit"]

        #print(str.join("; ",prep_commands))
        process = subprocess.check_call(str.join("; ",prep_commands),executable='/bin/bash',shell=True)

        print("\nCISM built!")
    
    print("\nSetting up regression tests directory")
    print(  "=====================================\n")
    data_dir = paths.mkdir_test(args, test_dict)
   
    print(  "   Copying CMake cache into regression test directory.")
    
    mod_list = paths.file_modifier_list(args)

    cache_mod = ""
    if mod_list:
        cache_mod = "-"+str.join("-", mod_list)
    
    cache_root = "CMakeCache"
    cache_ext = ".txt"
    cache_file = args.build_dir+os.sep+cache_root+cache_ext
    cache_new = data_dir+os.sep+cache_root+cache_mod+cache_ext
    
    subprocess.check_call("cp "+cache_file+" "+cache_new, shell=True)
    

    if isHPC:
        print("\nPreparing HPC batch jobs")
        print(  "========================\n")
        
        runnit.hpc(args, cism_driver, data_dir, test_dict)
        
        print("\nDone! You can now submit the job scripts.")
        print(  "=========================================")
    
    else:
        print("\nRunning regression tests")
        print(  "========================\n")
        
        runnit.personal(args, cism_driver, data_dir, test_dict)

        if args.timing:
            print("\nRe-running regression tests for timing data.")
            print(  "This is going to take a while. A long while.")
            print(  "============================================\n")
            
            for rnd in range(10):
                print("\nTiming round: "+str(rnd+1)+" of 10")
                args.tmod = 't'+str(rnd)
                runnit.personal(args, cism_driver, data_dir, test_dict)

            # turn timing modifier off now that we're done.
            args.tmod = None
            
    
        print("\nAll regression tests finished.")
        print(  "==============================")

if __name__=='__main__':
    args = parser.parse_args()
    main()
