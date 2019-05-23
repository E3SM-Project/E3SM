#!/usr/bin/env python

import os, shutil, sys

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *

logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")
#    cimeroot  = case.get_value("CIMEROOT")
    comp_wav  = case.get_value("COMP_WAV")

    atm_cpl_dt = coupling_times["atm_cpl_dt"]
    ocn_cpl_dt = coupling_times["ocn_cpl_dt"]

    outfile   = open(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), "w")

    if comp_wav == 'swav':
        outfile.write ("runSeq::                               \n")
        outfile.write ("@" + str(ocn_cpl_dt) + "               \n")
        outfile.write ("  MED med_phases_prep_ocn_accum_avg    \n")
        outfile.write ("  MED -> OCN :remapMethod=redist       \n")
        if (atm_cpl_dt < ocn_cpl_dt):
            outfile.write ("@" + str(atm_cpl_dt) + "           \n")
        outfile.write ("    MED med_phases_prep_ocn_map        \n")
        outfile.write ("    MED med_phases_aofluxes_run        \n")
        outfile.write ("    MED med_phases_prep_ocn_merge      \n")
        outfile.write ("    MED med_phases_prep_ocn_accum_fast \n")
        outfile.write ("    MED med_phases_ocnalb_run          \n")
        outfile.write ("    MED med_phases_prep_ice            \n")
        outfile.write ("    MED -> ICE :remapMethod=redist     \n")
        outfile.write ("    ICE                                \n")
        outfile.write ("    ROF                                \n")
        outfile.write ("    ATM                                \n")
        outfile.write ("    ICE -> MED :remapMethod=redist     \n")
        outfile.write ("    MED med_fraction_set               \n")
        outfile.write ("    ROF -> MED :remapMethod=redist     \n")
        outfile.write ("    ATM -> MED :remapMethod=redist     \n")
        outfile.write ("    MED med_phases_history_write       \n")
        outfile.write ("    MED med_phases_profile             \n")
        if (atm_cpl_dt < ocn_cpl_dt):
            outfile.write ("  @                                \n")
        outfile.write ("  OCN                                  \n")
        outfile.write ("  OCN -> MED :remapMethod=redist       \n")
        outfile.write ("  MED med_phases_restart_write         \n")
        outfile.write ("@                                      \n")
        outfile.write ("::                                     \n")

    elif comp_wav == 'ww' or comp_wav == "dwav":
        outfile.write ("runSeq::                             \n")
        outfile.write ("@" + str(atm_cpl_dt) + "             \n") # Assume that atm_cpl_dt >= ocn_cpl_dt
        outfile.write ("  MED med_phases_prep_ocn_map        \n") # map to ocean (including wav)
        outfile.write ("  MED med_phases_aofluxes_run        \n") # run atm/ocn flux calculation
        outfile.write ("  MED med_phases_prep_ocn_merge      \n")
        outfile.write ("  MED med_phases_prep_ocn_accum_fast \n")
        outfile.write ("  MED med_phases_prep_ocn_accum_avg  \n")
        outfile.write ("  MED med_phases_ocnalb_run          \n")
        outfile.write ("  MED -> OCN :remapMethod=redist     \n")
        outfile.write ("  MED med_phases_prep_ice            \n")
        outfile.write ("  MED -> ICE :remapMethod=redist     \n")
        outfile.write ("  MED med_phases_prep_wav            \n")
        outfile.write ("  MED -> WAV :remapMethod=redist     \n")
        outfile.write ("  ICE                                \n")
        outfile.write ("  ROF                                \n")
        outfile.write ("  WAV                                \n")
        outfile.write ("  ATM                                \n")
        outfile.write ("  ICE -> MED :remapMethod=redist     \n")
        outfile.write ("  MED med_fraction_set               \n")
        outfile.write ("  ROF -> MED :remapMethod=redist     \n")
        outfile.write ("  WAV -> MED :remapMethod=redist     \n")
        outfile.write ("  ATM -> MED :remapMethod=redist     \n")
        outfile.write ("  @ocn_cpl_dt   #ocean coupling step \n")
        outfile.write ("    OCN                              \n")
        outfile.write ("  @                                  \n")
        outfile.write ("  OCN -> MED :remapMethod=redist     \n")
        outfile.write ("  MED med_phases_restart_write       \n")
        outfile.write ("  MED med_phases_history_write       \n")
        outfile.write ("  MED med_phases_profile             \n")
        outfile.write ("@                                    \n")
        outfile.write ("::                                   \n")

    outfile.close()
    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
