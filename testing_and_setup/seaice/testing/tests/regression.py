#!/usr/bin/env python

import os, shutil
from compare_mpas_files import compare_files
from testing_utils import *

#-------------------------------------------------------------------------

def regression(mpasDevelopmentDir, mpasBaseDir, domainsDir, domain, configuration, options, check):

    # find available directory name
    iTest = 1
    dirExists = True
    while (dirExists):
        testDir = "regression_%i.%s.%s" %(iTest,configuration,domain)
        iTest = iTest + 1
        dirExists = os.path.isdir(testDir)

    # make a test directory
    create_test_directory(testDir)

    title = "Test: Regression, Configuration: %s, Domain: %s" %(configuration,domain)

    print_colour(title, "title")

    logfile = open("log_test.txt","w")
    logfile.write(title)

    # development run
    nProcs = 16

    nmlChanges = {"seaice_model": {"config_run_duration":'24:00:00'}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":"24:00:00"}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("development", mpasDevelopmentDir, domainsDir, domain, configuration, nmlChanges, streamChanges, nProcs, logfile) != 0):
        run_failed("regression")
        os.chdir("..")
        return 1

    # base run
    nProcs = 16

    nmlChanges = {"seaice_model": {"config_run_duration":'24:00:00'}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":"24:00:00"}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("base", mpasBaseDir, domainsDir, domain, configuration, nmlChanges, streamChanges, nProcs, logfile) != 0):
        run_failed("regression")
        os.chdir("..")
        return 1


    # compare
    restart_file = "restart.2000-01-02_00.00.00.nc"

    file1 = "./development/restarts/%s" %(restart_file)
    file2 = "./base/restarts/%s" %(restart_file)

    ignoreVarname = None
    if (check):
        ignoreVarname = ["testArrayParallelism","testArrayRestartability"]

    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile)

    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "regression")

    logfile.close()

    os.chdir("..")

    return failed

#-------------------------------------------------------------------------
