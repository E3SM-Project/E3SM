"""
Perturbation Growth New (PGN) - The CESM/ACME model's
multi-instance capability is used to conduct an ensemble
of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.

"""

from __future__ import division

import os
import re
import json
import shutil
import logging

from collections import OrderedDict

import CIME.test_status
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.build import post_build
from CIME.utils import expect

import evv4esm  # pylint: disable=import-error
from evv4esm.__main__ import main as evv  # pylint: disable=import-error

evv_lib_dir = os.path.abspath(os.path.dirname(evv4esm.__file__))

# Logic for PGN ensemble runs:
# ----------------------------
# We need two inputs:
# A. Number of inic cond files
# B. perturbations (e.g. currently we have: without prt, pos prt and neg prt)

# Based off the above, we compute number of instances to have (A * B)
# Build phase               : Change the user_nl* files to add perturbations and
#                             other flags
# Run phase                 : Rename history files
# Baselines generation phase: Compute cloud, store in netcdf file, copy netcdf
#                             file to baseline folder
# Baseline Comparison phase : Compare against baseline cloud to know pass/fail,
#                             and plot (optional)


# TODO:
# 1. The loop to change user_nl* files is run twice as we call build again after
#    changing ninst, which changes ntasks
# 2. Do we want to remove user_nl_cam, user_nl_clm etc. files as they are not
#    used in the simulation?
# 3. change code so that pergro_ptend_names.txt is not generated if it is
#    already there or only one instance writes this file...
# 4. Plot generation is very basic at this point(no labels, legends etc.),
#    improve it!
# 5. Decision making about PASS/FAIL should have multiple criteria

logger = logging.getLogger(__name__)

# --------------------------------------------------------
# Variables which needs global scope for various functions
# --------------------------------------------------------
# number of initial conditions
NINIT_COND = 6  # 12
# perturbations for runs
PERTURBATIONS = OrderedDict([('woprt', 0.0),
                             ('posprt', 1.0e-14),
                             ('negprt', -1.0e-14),
                             ])
# file name for file containing PGE cloud
FCLD_NC = 'cam.h0.cloud.nc'
# For preparing paths for namelist files for initial condition files
INIT_COND_FILE_TEMPLATE = \
    "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.{}.{}.0002-{:02d}-01-00000.nc"
# For finding the instance files
# FIXME: should 'cam' be 'atm' now?
INSTANCE_FILE_TEMPLATE = '{}.cam_{:04d}.h0.0001-01-01-00000.nc'
# ------------------------------------------------------------
# Some flags for debugging or invoking extra features
# ------------------------------------------------------------
# prints out max rmse diffs if set to True
INDEX = False
# Building the model can take a significant amount of time. Setting fake_bld to
# True can save that time
FAKE_BUILD = False


class PGN(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the PGN test
        """
        super(PGN, self).__init__(self, case)


    def build_phase(self, sharedlib_only=False, model_only=False):

        ninst = NINIT_COND * len(PERTURBATIONS)
        logger.debug('PGN_INFO: number of instance: '+str(ninst))

        # Find number of instance in the default setup
        default_ninst = self._case.get_value("NINST_ATM")

        # Sanity check: see if NINST is same for all model components, otherwise
        # exit with error
        for comp in self._case.get_values("COMP_CLASSES"):
            if comp != 'CPL':
                iinst = self._case.get_value("NINST_{}".format(comp))
                expect(default_ninst == iinst,
                       "ERROR: component {}  NINST({}) is different from  component"
                       " ATM NINST({}).".format(comp, iinst, default_ninst))

        # ------------------------------------------------------
        # Setup multi-instances for model components:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Set the model for multi-instance ONLY if NINST == 1
        # for all model components. This is because, for
        # NINST > 1 (e.g. rebuilding an old case) the following
        # loop will increase the ntasks to a multiple of ninst
        # (requiring a clean build again). We hit this issue if
        # we launch ./case.build in the case directory of PGN
        # test
        # ------------------------------------------------------

        if default_ninst == 1:  # if multi-instance is not already set
            # Only want this to happen once. It will impact the sharedlib build
            # so it has to happen here.
            if not model_only:
                # Lay all of the components out concurrently
                logger.debug("PGN_INFO: Updating NINST for multi-instance in "
                             "env_mach_pes.xml")
                for comp in ['ATM', 'OCN', 'WAV', 'GLC', 'ICE', 'ROF', 'LND']:
                    ntasks = self._case.get_value("NTASKS_{}".format(comp))
                    self._case.set_value("ROOTPE_{}".format(comp), 0)
                    self._case.set_value("NINST_{}".format(comp),  ninst)
                    self._case.set_value("NTASKS_{}".format(comp), ntasks*ninst)

                self._case.set_value("ROOTPE_CPL", 0)
                self._case.set_value("NTASKS_CPL", ntasks*ninst)
                self._case.flush()

                case_setup(self._case, test_mode=False, reset=True)

        if FAKE_BUILD:
            logger.debug("PGN_INFO: FAKE Build")
            if not sharedlib_only:
                post_build(self._case, [])
        else:
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

        logger.debug("PGN_INFO: Updating user_nl_* files")

        csmdata_root = self._case.get_value("DIN_LOC_ROOT")
        csmdata_atm = os.path.join(csmdata_root, "atm/cam/inic/homme/ne4_v1_init")
        csmdata_lnd = os.path.join(csmdata_root, "lnd/clm2/initdata/ne4_v1_init/b58d55680")

        iinst = 1
        for icond in range(1, NINIT_COND+1):
            fatm_in = os.path.join(csmdata_atm, INIT_COND_FILE_TEMPLATE.format('atm', 'i', icond))
            flnd_in =  os.path.join(csmdata_lnd, INIT_COND_FILE_TEMPLATE.format('clm2', 'r', icond))
            for iprt in PERTURBATIONS.values():
                with open('user_nl_cam_{:04d}'.format(iinst), 'w') as atmnlfile, \
                        open('user_nl_clm_{:04d}'.format(iinst), 'w') as lndnlfile:

                    atmnlfile.write("ncdata  = '{}' \n".format(fatm_in))
                    lndnlfile.write("finidat = '{}' \n".format(flnd_in))

                    # atm model output
                    atmnlfile.write("avgflag_pertape = 'I' \n")
                    atmnlfile.write("nhtfrq = 1 \n")
                    atmnlfile.write("mfilt  = 2  \n")
                    atmnlfile.write("ndens  = 1  \n")
                    atmnlfile.write("pergro_mods  = .true. \n")
                    atmnlfile.write("pergro_test_active = .true. \n")

                    # atmnlfile.write("empty_htapes = .true. \n")
                    # atmnlfile.write("fincl1 = 'PS','U','V','T','Q','CLDLIQ',"
                    #                 "'CLDICE','NUMLIQ','NUMICE','num_a1',"
                    #                 "'num_a2','num_a3','LANDFRAC' \n")
                    # atmnlfile.write("phys_debug_lat = 41.3495891345")
                    # atmnlfile.write("phys_debug_lon = 45.0" )

                    if iprt != 0.0:
                        atmnlfile.write("pertlim = {} \n".format(iprt))

                    iinst += 1

        self._case.set_value("STOP_N", "1")
        self._case.set_value("STOP_OPTION", "nsteps")


    # FIXME: Can this be done during __init__? Only should be done once right?
    def get_var_list(self):
        """
        Get variable list for pergro specific output vars
        """
        rundir = self._case.get_value("RUNDIR")
        prg_fname = 'pergro_ptend_names.txt'
        var_file = os.path.join(rundir, prg_fname)
        expect(os.path.isfile(var_file),
               "File {} does not exist in: {}".format(prg_fname, rundir))

        with open(var_file, 'r') as fvar:
            var_list = fvar.readlines()

        return map(str.strip, var_list)


    def _compare_baseline(self):
        """
        Compare baselines in the pergro test sense. That is,
        compare PGE from the test simulation with the baseline 
        cloud
        """
        with self._test_status:
            self._test_status.set_status(CIME.test_status.BASELINE_PHASE,
                                         CIME.test_status.TEST_FAIL_STATUS)

            logger.debug("PGN_INFO:BASELINE COMPARISON STARTS")

            run_dir = self._case.get_value("RUNDIR")
            case_name = self._case.get_value("CASE")
            base_dir = os.path.join(self._case.get_value("BASELINE_ROOT"),
                                    self._case.get_value("BASECMP_CASE"))

            var_list = self.get_var_list()
            path_cld_nc = os.path.join(base_dir, FCLD_NC)

            test_name = "{}".format(case_name.split('.')[-1])
            evv_config = {
                test_name: {
                    "module": os.path.join(evv_lib_dir, "extensions", "pg.py"),
                    "test-case": case_name,
                    "test-name": "Test",
                    "test-dir": run_dir,
                    "ref-name": "Baseline",
                    "ref-dir": base_dir,
                    "variables": var_list,
                    "perturbations": PERTURBATIONS,
                    "pge-cld": FCLD_NC,
                    "ninit": NINIT_COND,
                    "init-file-template": INIT_COND_FILE_TEMPLATE,
                    "instance-file-template": INSTANCE_FILE_TEMPLATE,

                }
            }

            json_file = os.path.join(run_dir, '.'.join([case_name, 'json']))
            with open(json_file, 'w') as config_file:
                json.dump(evv_config, config_file, indent=4)

            evv_out_dir = os.path.join(run_dir, '.'.join([case_name, 'evv']))
            evv(['-e', json_file, '-o', evv_out_dir])

            with open(os.path.join(evv_out_dir, 'index.json'), 'r') as evv_f:
                evv_status = json.load(evv_f)

            for evv_elem in evv_status['Data']['Elements']:
                # FIXME: TableTitle!
                if evv_elem['Type'] == 'ValSummary' \
                        and evv_elem['TableTitle'] == 'Peturbation-Growth':
                    if evv_elem['Data'][test_name]['']['Ensembles'] == 'identical':
                        self._test_status.set_status(CIME.test_status.BASELINE_PHASE,
                                                     CIME.test_status.TEST_PASS_STATUS)
                        break

    # FIXME: REFACTOR
    def run_phase(self):
        logger.debug("PGN_INFO: RUN PHASE")

        self.run_indv()

        # Here were are in case directory, we need to go to the run directory
        # and rename files
        rundir = self._case.get_value("RUNDIR")
        casename = self._case.get_value("CASE")
        logger.debug("PGN_INFO: Case name is:{}".format(casename))

        for icond in range(NINIT_COND):
            for iprt, (prt, prt_name) in enumerate(PERTURBATIONS.items()):

                # ------------------------------------------------------
                # SANITY CHECK - To confirm that history file extension
                # ~~~~~~~~~~~~
                # corresponds to the right perturbation
                # ------------------------------------------------------
                # find corresponding atm_in_*
                iinst = _sub2instance(icond, iprt)
                fatm_in = os.path.join(rundir, 'atm_in_{:04d}'.format(iinst))
                # see if atm_in_* file has pertlim same as PRT[iprt] (sanity check)
                found = False
                prtval = 0.0
                with open(fatm_in) as atmfile:
                    for line in atmfile:
                        if line.find('pertlim') > 0:
                            found = True
                            prtval = float(line.split('=')[1])
                            expect(prtval == prt,
                                   "ERROR: prtval doesn't match, "
                                   "prtval:{}; prt[{}]:{}".format(
                                           prtval, iprt, prt))
                            logger.debug("PGN_INFO:prtval:{}; prt[{}]:{}".format(
                                    prtval, iprt, prt))

                if not found:
                    expect(prtval == prt,
                           "ERROR: default prtval doesn't match, "
                           "prtval:{}; prt[{}]:{}".format(prtval, iprt, prt))
                    logger.debug("PGN_INFO:def prtval:{}; prt[{}]:{}".format(
                            prtval, iprt, prt))

                # ---------------------------------------------------------
                # Rename file
                # ---------------------------------------------------------
                fname = os.path.join(rundir, INSTANCE_FILE_TEMPLATE.format(casename, iinst))
                renamed_fname = re.sub(r'\.nc$', '_{}.nc'.format(prt_name), fname)

                logger.debug("PGN_INFO: fname to rename:{}".format(fname))
                logger.debug("PGN_INFO: Renamed file:{}".format(renamed_fname))
                try:
                    shutil.move(fname, renamed_fname)
                except IOError:
                    expect(os.path.isfile(renamed_fname),
                           "ERROR: File {} does not exist".format(renamed_fname))
                    logger.debug("PGN_INFO: Renamed file already exists:"
                                 "{}".format(renamed_fname))

        # cloud generation

        logger.debug("PGN_INFO: cloud generation-gen base")

        var_list = self.get_var_list()
        len_var_list = len(var_list)
        # nprt = len(PRT)

        # for trusted cloud sims
        logger.debug("PGN_INFO: Computing cloud")

        # cld_res = np.empty([len_var_list])

        # ---------------------------------------------
        # Write netcdf file for cloud in rundir
        # ---------------------------------------------
        os.chdir(rundir)
        fcld, cld_rmse_nc = self.nc_write_handle(FCLD_NC, 'cld_rmse')

        iinst = 0
        for icond in range(NINIT_COND):
            iinst += 1
            ifile_cntl = os.path.join(rundir, '{}_{}.nc'.format(
                    self.get_fname_wo_ext('', casename, iinst), PERTURBATIONS))
            expect(os.path.isfile(ifile_cntl),
                   "ERROR: File {} does not exist".format(ifile_cntl))
            logger.debug("PGN_INFO:CNTL_CLD:{}".format(ifile_cntl))

            for iprt, prt_name in enumerate(PERTURBATIONS):
                if prt_name == 'woprt':
                    continue
                iinst += 1

                ifile_test = os.path.join(rundir, '{}_{}.nc'.format(
                        self.get_fname_wo_ext('', casename, iinst), prt_name))
                expect(os.path.isfile(ifile_test),
                       "ERROR: File {} does not exist".format(ifile_test))
                # NOTE: iprt-1 as we will only get two curves for each inic
                # (wo-pos and wo-neg)
                cld_rmse_nc[icond, iprt-1, 0:len_var_list] = self.rmse_var(
                        ifile_test, ifile_cntl,  var_list, 't_')
                logger.debug("PGN_INFO:Compared to CLD:{}".format(ifile_test))

        fcld.close()

        if self._case.get_value("GENERATE_BASELINE"):

            # baseline directory names
            base_root = self._case.get_value("BASELINE_ROOT")
            base_gen = self._case.get_value("BASEGEN_CASE")

            # baseline directory is:base_root/base_gen
            base_dir = os.path.join(base_root, base_gen)

            # first copy files to the baseline directory
            self._generate_baseline()  # BALLI-CHECK IF THIS IS OKAY

            # copy cloud.nc file to baseline directory
            logger.debug("PGN_INFO:copy:{} to {}".format(FCLD_NC, base_dir))
            shutil.copy(FCLD_NC, base_dir)

        logger.debug("PGN_INFO: RUN PHASE ENDS")

# =====================================================
# Debugging:
# =====================================================

# -----------------------------------------------------
# DEBUG type 1 (DB1): Ensure that model produces BFB instances
# -----------------------------------------------------
# Replace PERTURBATIONS dictionary at the top of the script with the following:
# PERTURBATIONS = OrderedDict([('woprt', 0.0), ('posprt', 0.0)])

# Comment out all namelist changes so that all instances have same namelist values
# For testing effect of namelist changes, uncomment namelist changes such that
# namelists are same for all instances

# Comment out sanity checks as well if needed....


def _instance2sub(instance_number):
    """
    Converts an instance number (ii) to initial condition index (ci) and
    perturbation index (pi)  subscripts

    instances use 1-based indexes and vary according to this function:
        ii = ci * len(PERTURBATIONS) + pi + 1
    where both pi and ci use 0-based indexes.
    """
    perturbation_index = (instance_number - 1) % len(PERTURBATIONS)
    initial_condition = (instance_number - 1 - perturbation_index) // len(PERTURBATIONS)
    return initial_condition, perturbation_index


def _sub2instance(initial_condition, perturbation_index):
    """
    Converts initial condition index (ci) and perturbation index (pi) subscripts
    to an instance number (ii)

    instances use 1-based indexes and vary according to this function:
        ii = ci * len(PERTURBATIONS) + pi + 1
    where both pi and ci use 0-based indexes.
    """
    instance = initial_condition * len(PERTURBATIONS) + perturbation_index + 1
    return instance
