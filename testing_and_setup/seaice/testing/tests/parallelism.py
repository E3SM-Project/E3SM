#!/usr/bin/env python

import os, shutil
from compare_mpas_files import compare_files
from testing_utils import *

#-------------------------------------------------------------------------

def parallelism(mpasDevelopmentDir, domainsDir, domain, configuration, options, check):

    # find available directory name
    iTest = 1
    dirExists = True
    while (dirExists):
        testDir = "parallelism_%i.%s.%s" %(iTest,configuration,domain)
        iTest = iTest + 1
        dirExists = os.path.isdir(testDir)

    # make a test directory
    create_test_directory(testDir)

    title = "Test: Parallelism, Configuration: %s, Domain: %s" %(configuration,domain)

    print_colour(title, "title")

    logfile = open("log_test.txt","w")
    logfile.write(title)

    multipleBlocks = False
    if ("multipleBlocks" in options.keys() and options["multipleBlocks"] == "True"):
        multipleBlocks = True

    print "multipleBlocks: ", multipleBlocks
    logfile.write("multipleBlocks: %s" %(multipleBlocks))

    # development run
    nProcs = 16

    nmlChanges = {"seaice_model": {"config_run_duration":'24:00:00'}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":"24:00:00"}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("development1", mpasDevelopmentDir, domainsDir, domain, configuration, nmlChanges, streamChanges, nProcs, logfile) != 0):
        run_failed("parallelism")
        os.chdir("..")
        return 1

    # base run
    nProcs = 32

    if (not multipleBlocks):
        nmlChanges = {"seaice_model": {"config_run_duration":'24:00:00'}}
    else:
        nmlChanges = {"seaice_model": {"config_run_duration":'24:00:00'},
                     "decomposition": {"config_block_decomp_file_prefix":'graphs/graph.info.eq.part.',
                                       "config_number_of_blocks": 96,
                                       "config_explicit_proc_decomp": True,
                                       "config_proc_decomp_file_prefix":'graphs/graph.info.eq_block.part.'}}

    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":"24:00:00"}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("development2", mpasDevelopmentDir, domainsDir, domain, configuration, nmlChanges, streamChanges, nProcs, logfile) != 0):
        run_failed("parallelism")
        os.chdir("..")
        return 1


    # compare
    restart_file = "restart.2000-01-02_00.00.00.nc"

    file1 = "./development1/restarts/%s" %(restart_file)
    file2 = "./development2/restarts/%s" %(restart_file)

    ignoreVarname = ["cellsOnCell","verticesOnCell","edgesOnEdge","edgesOnCell"]
    if (check):
        ignoreVarname.append("testArrayReproducibility")
        ignoreVarname.append("testArrayRestartability")

    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile,ignoreVarname)

    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "parallelism")

    logfile.close()

    os.chdir("..")

    return failed

#-------------------------------------------------------------------------
