#!/usr/bin/env python

import os, shutil, sys

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *

logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")

    atm_cpl_dt = coupling_times["atm_cpl_dt"]
    ocn_cpl_dt = coupling_times["ocn_cpl_dt"]

    outfile   = open(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), "w")

    if case.get_value("CONTINUE_RUN") or case.get_value("MEDIATOR_READ_RESTART"):
        logger.info("NUOPC run sequence: warm start (concurrent)")

        outfile.write ("runSeq::                                \n")
        outfile.write ("@" + str(ocn_cpl_dt) + "                \n")
        outfile.write ("   MED med_phases_prep_ocn_accum_avg    \n")
        outfile.write ("   MED -> OCN :remapMethod=redist       \n")
        outfile.write ("   OCN                                  \n")
        outfile.write ("   @" + str(atm_cpl_dt) + "             \n")
        outfile.write ("     MED med_phases_prep_atm            \n")
        outfile.write ("     MED med_phases_prep_ice            \n")
        outfile.write ("     MED -> ATM :remapMethod=redist     \n")
        outfile.write ("     MED -> ICE :remapMethod=redist     \n")
        outfile.write ("     ATM                                \n")
        outfile.write ("     ICE                                \n")
        outfile.write ("     ATM -> MED :remapMethod=redist     \n")
        outfile.write ("     ICE -> MED :remapMethod=redist     \n")
        outfile.write ("     MED med_fraction_set               \n")
        outfile.write ("     MED med_phases_prep_ocn_map        \n")
        outfile.write ("     MED med_phases_aofluxes_run        \n")
        outfile.write ("     MED med_phases_prep_ocn_merge      \n")
        outfile.write ("     MED med_phases_prep_ocn_accum_fast \n")
        outfile.write ("     MED med_phases_history_write       \n")
        outfile.write ("     MED med_phases_profile             \n")
        outfile.write ("   @                                    \n")
        outfile.write ("   OCN -> MED :remapMethod=redist       \n")
        outfile.write ("   MED med_phases_restart_write         \n")
        outfile.write ("@                                       \n")
        outfile.write ("::                                      \n")

    else:
        logger.info("NUOPC run sequence: cold start (sequential)")
        outfile.write ("runSeq::                                \n")
        outfile.write ("@" + str(ocn_cpl_dt) + "                \n")
        outfile.write ("   @" + str(atm_cpl_dt) + "             \n")
        outfile.write ("     MED med_phases_prep_atm            \n")
        outfile.write ("     MED -> ATM :remapMethod=redist     \n")
        outfile.write ("     ATM                                \n")
        outfile.write ("     ATM -> MED :remapMethod=redist     \n")
        outfile.write ("     MED med_phases_prep_ice            \n")
        outfile.write ("     MED -> ICE :remapMethod=redist     \n")
        outfile.write ("     ICE                                \n")
        outfile.write ("     ICE -> MED :remapMethod=redist     \n")
        outfile.write ("     MED med_fraction_set               \n")
        outfile.write ("     MED med_phases_prep_ocn_map        \n")
        outfile.write ("     MED med_phases_aofluxes_run        \n")
        outfile.write ("     MED med_phases_prep_ocn_merge      \n")
        outfile.write ("     MED med_phases_prep_ocn_accum_fast \n")
        outfile.write ("     MED med_phases_history_write       \n")
        outfile.write ("     MED med_phases_profile             \n")
        outfile.write ("   @                                    \n")
        outfile.write ("   MED med_phases_prep_ocn_accum_avg    \n")
        outfile.write ("   MED -> OCN :remapMethod=redist       \n")
        outfile.write ("   OCN                                  \n")
        outfile.write ("   OCN -> MED :remapMethod=redist       \n")
        outfile.write ("   MED med_phases_restart_write         \n")
        outfile.write ("@                                       \n")
        outfile.write ("::                                      \n")

    outfile.close()
    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
