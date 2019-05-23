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
#    comp_atm  = case.get_value("COMP_ATM")
    comp_glc  = case.get_value("COMP_GLC")
#    comp_lnd  = case.get_value("COMP_LND")
    comp_rof  = case.get_value("COMP_ROF")

    outfile   = open(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), "w")

#    glc_cpl_dt = coupling_times["glc_cpl_dt"]
    rof_cpl_dt = coupling_times["rof_cpl_dt"]
    atm_cpl_dt = coupling_times["atm_cpl_dt"]

    if (comp_rof == 'srof' and comp_glc == "sglc"):

        outfile.write ("runSeq::                                \n" )
        outfile.write ("@" + str(atm_cpl_dt) + "                \n" )
        outfile.write ("  MED med_phases_prep_lnd               \n" )
        outfile.write ("  MED -> LND :remapMethod=redist        \n" )
        outfile.write ("  LND                                   \n" )
        outfile.write ("  LND -> MED :remapMethod=redist        \n" )
        outfile.write ("  MED med_fraction_set                  \n" )
        outfile.write ("  ATM                                   \n" )
        outfile.write ("  ATM -> MED :remapMethod=redist        \n" )
        outfile.write ("  MED med_phases_history_write          \n" )
        outfile.write ("  MED med_phases_profile                \n" )
        outfile.write ("  MED med_phases_restart_write          \n" )
        outfile.write ("@                                       \n" )
        outfile.write ("::                                      \n" )

    elif ((comp_rof == 'mosart' or comp_rof == 'rtm') and comp_glc == "sglc"):

        outfile.write ("runSeq::                                  \n" )
        outfile.write ("@" + str(rof_cpl_dt) + "                  \n" )
        outfile.write ("  MED med_phases_prep_rof_avg             \n" )
        outfile.write ("  MED -> ROF :remapMethod=redist          \n" )
        outfile.write ("  ROF                                     \n" )
        if atm_cpl_dt < rof_cpl_dt:
            outfile.write ("@" + str(atm_cpl_dt) + "              \n" )
        outfile.write ("  MED med_phases_prep_lnd                 \n" )
        outfile.write ("  MED -> LND :remapMethod=redist          \n" )
        outfile.write ("  LND                                     \n" )
        outfile.write ("  LND -> MED :remapMethod=redist          \n" )
        outfile.write ("  MED med_phases_prep_rof_accum_fast      \n" )
        outfile.write ("  ATM                                     \n" )
        outfile.write ("  ATM -> MED :remapMethod=redist          \n" )
        outfile.write ("  MED med_phases_profile                  \n" )
        if atm_cpl_dt < rof_cpl_dt:
            outfile.write ("@                                     \n" )
        outfile.write ("  ROF -> MED :remapMethod=redist          \n" )
        outfile.write ("  MED med_phases_history_write            \n" )
        outfile.write ("  MED med_phases_restart_write            \n" )
        outfile.write ("@                                         \n" )
        outfile.write ("::                                        \n" )

    elif ((comp_rof == 'mosart' or comp_rof == 'rtm') and comp_glc == "cism"):

        # TODO: need to add CISM below
        outfile.write ("runSeq::                                  \n" )
        outfile.write ("@" + str(rof_cpl_dt) + "                  \n" )
        outfile.write ("  MED med_phases_prep_rof_avg             \n" )
        outfile.write ("  MED -> ROF :remapMethod=redist          \n" )
        outfile.write ("  ROF                                     \n" )
        if atm_cpl_dt < rof_cpl_dt:
            outfile.write ("@" + str(atm_cpl_dt) + "              \n" )
        outfile.write ("  @atm_cpl_dt  # atmosphere coupling step \n" )
        outfile.write ("  MED med_phases_prep_lnd                 \n" )
        outfile.write ("  MED -> LND :remapMethod=redist          \n" )
        outfile.write ("  LND                                     \n" )
        outfile.write ("  LND -> MED :remapMethod=redist          \n" )
        outfile.write ("  MED med_phases_prep_rof_accum_fast      \n" )
        outfile.write ("  ATM                                     \n" )
        outfile.write ("  ATM -> MED :remapMethod=redist          \n" )
        outfile.write ("  MED med_phases_profile                  \n" )
        if atm_cpl_dt < rof_cpl_dt:
            outfile.write ("@                                     \n" )
        outfile.write ("  ROF -> MED :remapMethod=redist          \n" )
        outfile.write ("  MED med_phases_history_write            \n" )
        outfile.write ("  MED med_phases_restart_write            \n" )
        outfile.write ("@                                         \n" )
        outfile.write ("::                                        \n" )

    outfile.close()
    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
