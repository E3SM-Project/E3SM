"""
long term archiving
"""

import shutil, glob, re, os

from CIME.XML.standard_module_setup import *
from CIME.utils                     import expect, does_file_have_string, append_status, run_cmd
from distutils.spawn import find_executable

import time

logger = logging.getLogger(__name__)

###############################################################################
def _copy_dirs_hsi(dout_s_root, dout_l_msroot, dout_l_hpss_accnt, dout_l_dryrun, 
                   dout_l_delete):
###############################################################################

    logger.debug('In copy_dirs_hsi...')

    success = False

    hsi = find_executable("hsi")
    expect(hsi == None,"lt_archive: asked for copy_dirs_hsi - but hsi not found, check path")

    # check the local delete optin
    saveFlag="-PR"
    if dout_l_delete:
        saveFlag="-PRd"

    # send files to HPSS
    # Note - setting chmod +t applies the sticky bit permission mask
    os.chdir(dout_s_root)
    hsiArgs='mkdir -p' + dout_l_msroot + ' ; chmod +t ' + dout_l_msroot + ' ; cd ' + dout_l_msroot + ' ; put ' + saveFlag + ' *'

    logger.debug(hsiArgs)

    if not dout_l_hpss_accnt.startswith('0000'):
        if dout_l_dryrun:
           print ('hsi -a ' + dout_l_hpss_accnt + ' ' + hsiArg)
        else:
	   stat, msg, error = run_cmd('hsi -a ' + dout_l_hpss_accnt + ' ' + hsiArgs)
    else:
        if dout_l_dryrun:
           print ('hsi ' + hsiArg)
        else:
	   stat, msg, error = run_cmd('hsi ' + hsiArgs)
    
    return stat, msg


###############################################################################
def _copy_files(dout_s_root, dout_l_msroot, dout_l_hpss_accnt, dout_l_dryrun,
                dout_l_delete):
###############################################################################

    logger.debug('In copy_files...')
    success = True
    msg = ""
    
    return success, msg


###############################################################################
def _copy_dirs_ssh(dout_s_root, dout_l_msroot, dout_l_ssh_loc, dout_l_dryrun,
                   dout_l_delete):
###############################################################################

    logger.debug('In copy_dirs_ssh...')
    success = True
    msg = ""
    
    return success, msg


###############################################################################
def _copy_dirs_local(dout_s_root, dout_l_msroot, dout_l_arc_root, dout_l_dryrun,
                     dout_l_delete):
###############################################################################

    logger.debug('In copy_dirs_local...')
    success = True
    msg = ""
    
    return success, msg


###############################################################################
def case_lt_archive(case):
###############################################################################
    """"
    perform long term archiving of DOUT_S_ROOT files to HPSS
    """
    caseroot = case.get_value("CASEROOT")

    # max number of threads needed by scripts
    os.environ["maxthrds"] = 1

    # document start
    append_status("lt_archive starting",caseroot=caseroot,sfile="CaseStatus")

    # determine status of run and short term archiving
    # TODO want to reverse search from the bottom of these files
    # add more to this logic for edge cases
    runComplete = does_file_have_string(os.path.join(caseroot, "CaseStatus"),
                                        "run SUCCESSFUL")
    staComplete = does_file_have_string(os.path.join(caseroot, "CaseStatus"),
                                        "st_archive_complete")

    # get env variables and call the different archive methods based on mode
    if runComplete and staComplete:
        dout_s_root = case.get_value("DOUT_S_ROOT")
        dout_l_msroot = case.get_value("DOUT_L_MSROOT")
        dout_l_hpss_accnt = case.get_value("DOUT_L_HPSS_ACCNT")
        dout_l_mode = case.get_value("DOUT_L_MODE")
        dout_l_dryrun = case.get_value("DOUT_L_DRYRUN")
        dout_l_delete = case.get_value("DOUT_L_DELETE_LOC_FILES")
        dout_l_ssh_loc = case.get_value("DOUT_L_SSH_LOC")
        dout_l_arc_root = case.get_value("DOUT_L_ARC_ROOT")
        lid = time.strftime("%y%m%d-%H%M%S")

        msg = ""
        # perform archiving based on the mode requested
        if dout_l_mode == "copy_dirs_hsi":
           (success, msg) = _copy_dirs_hsi(dout_s_root, dout_l_msroot, dout_l_hpss_accnt, 
                                           dout_l_dryrun, dout_l_delete)
        elif dout_l_mode == "copy_files":
           (success, msg) = _copy_files(dout_s_root, dout_l_msroot, dout_l_hpss_accnt, 
                                        dout_l_dryrun, dout_l_delete)
        elif dout_l_mode == "copy_dirs_ssh":
           (success, msg) = _copy_dirs_ssh(dout_s_root, dout_l_msroot, dout_l_ssh_loc, 
                                           dout_l_dryrun, dout_l_delete)
        elif dout_l_mode == "copy_dirs_local":
           (success, msg) = _copy_dirs_local(dout_s_root, dout_l_msroot, dout_l_arc_loc, 
                                             dout_l_dryrun, dout_l_delete)
        else:
            expect(False,
                   "lt_archive: unrecognized DOUT_L_MODE '"+dout_l_mode+"'."
                   "Unable to perform long term archive...")

        if not success:
            expect(False,
                   "lt_archive: "+msg+
                   "Unable to perform long term archive...")

    else:
        expect(False,
               "lt_archive: run or st_archive is not yet complete or was not successful."
               "Unable to perform long term archive...")

    # document completion
    append_status("lt_archive completed" ,caseroot=caseroot, sfile="CaseStatus")

    return True
