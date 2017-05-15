#!/usr/bin/env python

import os, shutil, sys

#-------------------------------------------------------------------------

def colour_init():

    try:
        import colorama
        colorama.init()
    except ImportError:
        pass

#-------------------------------------------------------------------------

def print_colour(message, messageType):

    try:
        import colorama

        if (messageType == "fail"):

            print(colorama.Fore.RED + colorama.Style.BRIGHT + \
                "FAIL!:" + \
                colorama.Style.RESET_ALL + colorama.Fore.RED + \
                " " + message + \
                colorama.Style.RESET_ALL)

        elif (messageType == "pass"):

            print(colorama.Fore.GREEN + colorama.Style.BRIGHT + \
                "PASS:" + \
                colorama.Style.RESET_ALL + colorama.Fore.GREEN + \
                " " + message + \
                colorama.Style.RESET_ALL)

        elif (messageType == "title"):

            print
            print(colorama.Style.BRIGHT + colorama.Fore.MAGENTA + message + colorama.Style.RESET_ALL)
            underline = "=" * len(message)
            print(colorama.Style.BRIGHT + colorama.Fore.MAGENTA + underline + colorama.Style.RESET_ALL)

        elif (messageType == "green"):

            print(colorama.Fore.GREEN + message + colorama.Style.RESET_ALL)

        elif (messageType == "red"):

            print(colorama.Fore.RED + message + colorama.Style.RESET_ALL)

    except ImportError:

        if (messageType == "fail"):

            print("FAIL!: " + message)

        elif (messageType == "pass"):

            print("PASS: " + message)

        elif (messageType == "title"):

            print
            print(message)
            underline = "=" * len(message)
            print(underline)

        else:

            print(message)

#-------------------------------------------------------------------------

def test_summary(nErrorsNonArray, nErrorsArray, logfile, testname):

    if (nErrorsNonArray == 0):
        print "No non-contents errors"
        logfile.write("No non-contents errors\n")
    else:
        print_colour("%i non-content errors" %(nErrorsNonArray), "red")
        logfile.write("%i non-content errors\n" %(nErrorsNonArray))

    if (nErrorsArray == 0):
        print_colour("Test %s passed" %(testname), "pass")
        logfile.write("Test %s passed\n" %(testname))
        failed = 0
    else:
        print_colour("Test %s failed" %(testname), "fail")
        logfile.write("Test %s failed\n" %(testname))
        failed = 1

    return failed

#-------------------------------------------------------------------------

def final_summary(nTests, nFails):

    nPasses = nTests - nFails

    try:
        import colorama

        if (nFails == 0):

            print("%i tests (" %(nTests) + \
                colorama.Fore.GREEN + colorama.Style.BRIGHT + \
                "%i passed" %(nPasses) + \
                colorama.Style.RESET_ALL + \
                ")")

        elif (nPasses == 0):

            print("%i tests (" %(nTests) + \
                colorama.Fore.RED + colorama.Style.BRIGHT + \
                "%i failed" %(nFails) + \
                colorama.Style.RESET_ALL + \
                ")")

        else:

            print("%i tests (" %(nTests) + \
                colorama.Fore.GREEN + colorama.Style.BRIGHT + \
                "%i passed" %(nPasses) + \
                colorama.Style.RESET_ALL + \
                ", " + \
                colorama.Fore.RED + colorama.Style.BRIGHT + \
                "%i failed" %(nFails) + \
                colorama.Style.RESET_ALL + \
                ")")

    except ImportError:

        if (nFails == 0):
            print("%i tests (%i passed)" %(nTests,nPasses))
        elif (nPasses == 0):
            print("%i tests (%i failed)" %(nTests,nFails))
        else:
            print("%i tests (%i passed, %i failed)" %(nTests, nPasses, nFails))

#-------------------------------------------------------------------------
# namelist manupulation
#-------------------------------------------------------------------------

def create_new_namelist(filenameIn, filenameOut, nmlPatch):

    try:
        import f90nml
    except ImportError:
        print "Module f90nml needed and not available"
        sys.exit()

    f90nml.patch(filenameIn, nmlPatch, filenameOut)

#-------------------------------------------------------------------------
# streams manupulation
#-------------------------------------------------------------------------

def create_new_streams(filenameIn, filenameOut, changes):

    try:
        import xml.etree.ElementTree as ET
    except ImportError:
        print "Module xml.etree.ElementTree needed and not available"
        sys.exit()

    tree = ET.parse(filenameIn)
    streams = tree.getroot()

    for stream in streams:
        for change in changes:
            if (stream.get('name') == change["streamName"]):
                stream.set(change["attributeName"],change["newValue"])

    tree.write(filenameOut)

#-------------------------------------------------------------------------
# run the model
#-------------------------------------------------------------------------

def create_test_directory(directory):

    if (os.path.exists(directory)):
        shutil.rmtree(directory)
    os.mkdir(directory)
    os.chdir(directory)

#-------------------------------------------------------------------------

def run_model(runName, mpasDir, domainsDir, domain, configuration, nmlChanges, streamChanges, nProcs, logfile):

    # create development directory
    os.mkdir(runName)
    os.chdir(runName)

    # sym link executable
    os.symlink(mpasDir + "/cice_model", "cice_model")

    # get domain
    cmd = domainsDir + "/" + domain + "/get_domain.py"
    os.system(cmd)

    # create namelist file
    create_new_namelist(mpasDir+"/testing_and_setup/seaice/configurations/"+configuration+"/namelist.cice", "namelist.cice", nmlChanges)

    # create streams file
    create_new_streams(mpasDir+"/testing_and_setup/seaice/configurations/"+configuration+"/streams.cice", "streams.cice", streamChanges)

    # run the model
    returnCode = execute_model(nProcs, logfile)

    # up one level
    os.chdir("..")

    return returnCode

#-------------------------------------------------------------------------

def restart_model(runName, nmlChanges, streamChanges, nProcs, logfile):

    # move to run dir
    os.chdir(runName)

    # update namelist
    os.rename("namelist.cice","namelist.cice.prev")
    os.rename("streams.cice","streams.cice.prev")

    # update namelist
    create_new_namelist("namelist.cice.prev", "namelist.cice", nmlChanges)

    # create streams file
    create_new_streams("streams.cice.prev", "streams.cice", streamChanges)

    # run the model
    returnCode = execute_model(nProcs, logfile)

    # up one level
    os.chdir("..")

    return returnCode

#-------------------------------------------------------------------------

def execute_model(nProcs, logfile):

    import subprocess

    # execute mpirun -np np cice_model
    process = subprocess.Popen(["mpirun", "-np", "%i" %(nProcs), "cice_model"], stdout=logfile, stderr=logfile)
    returnCode = process.wait()
    logfile.flush()
    logfile.write("Return code: %i\n" %(returnCode))
    logfile.flush()

    return returnCode

#-------------------------------------------------------------------------
