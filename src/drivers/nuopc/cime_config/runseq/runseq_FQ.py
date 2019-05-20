#!/usr/bin/env python

import os, shutil, sys, glob, itertools

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *
from CIME.case import Case
from CIME.utils import expect
from textwrap import dedent

logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")
    cimeroot  = case.get_value("CIMEROOT")
    comp_atm  = case.get_value("COMP_ATM")
    comp_ice  = case.get_value("COMP_ICE")
    comp_ocn  = case.get_value("COMP_OCN")
    comp_lnd  = case.get_value("COMP_LND")
    comp_rof  = case.get_value("COMP_ROF")
    comp_glc  = case.get_value("COMP_GLC")
    docn_mode = case.get_value("DOCN_MODE")

    output    = os.path.join(caseroot, "nuopc.runseq")
    outfile   = open(output, "w")

    # determine if prognostic ocn
    prognostic_ocn = ('som' in docn_mode or 'interannual' in docn_mode)

    # TODO: need to add cism

    outfile.write ("runSeq::                                  \n")
    if comp_rof != xrof:
        outfile.write ("@rof_cpl_dt                           \n")
        outfile.write ("   MED med_phases_prep_rof_avg        \n")
        outfile.write ("   MED -> ROF :remapMethod=redist     \n")
    if (atm_cpl_dt < rof_cpl_dt):
        outfile.write ("@atm_cpl_dt                           \n")
    if prognostic_ocn:
        outfile.write ("   MED med_phases_prep_ocn_accum_avg  \n")
        outfile.write ("   MED -> OCN :remapMethod=redist     \n")
    if (comp_lnd != "xlnd"):
        outfile.write ("   MED med_phases_prep_lnd            \n")
        outfile.write ("   MED -> LND :remapMethod=redist     \n")
        outfile.write ("   LND                                \n")
        outfile.write ("   LND -> MED :remapMethod=redist     \n")
    if (comp_ice != "xice"):
        outfile.write ("   MED med_phases_prep_ice            \n")
        outfile.write ("   MED -> ICE :remapMethod=redist     \n")
        outfile.write ("   ICE                                \n")
        outfile.write ("   ICE -> MED :remapMethod=redist     \n")
    outfile.write ("   OCN                                    \n")
    outfile.write ("   OCN -> MED :remapMethod=redist         \n")
    outfile.write ("   MED med_fraction_set                   \n")
    outfile.write ("   MED med_phases_aofluxes_run            \n")
    if prognostic_ocn:
        outfile.write ("   MED med_phases_prep_ocn_map        \n")
        outfile.write ("   MED med_phases_prep_ocn_merge      \n")
        outfile.write ("   MED med_phases_prep_ocn_accum_fast \n")
    outfile.write ("   MED med_phases_ocnalb_run              \n")
    outfile.write ("   MED med_phases_prep_atm                \n")
    outfile.write ("   MED -> ATM :remapMethod=redist         \n")
    outfile.write ("   ATM                                    \n")
    outfile.write ("   ATM -> MED :remapMethod=redist         \n")
    if (atm_cpl_dt < rof_cpl_dt):
        outfile.write ("   MED med_phases_history_write       \n")
    if comp_rof != "xrof":
        outfile.write ("@                                     \n")
    outfile.write ("   MED med_phases_restart_write           \n")
    outfile.write ("   MED med_phases_history_write           \n")
    outfile.write ("   MED med_phases_profile                 \n")
    outfile.write ("@                                         \n")
    outfile.write ("::                                        \n")

    outfile.close()
    shutil.copy(output, rundir)
