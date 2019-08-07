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

import pandas as pd
import numpy as np
from collections import OrderedDict

import CIME.test_status
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.utils import expect

import evv4esm  # pylint: disable=import-error
from evv4esm.extensions import pg  # pylint: disable=import-error
from evv4esm.__main__ import main as evv  # pylint: disable=import-error

evv_lib_dir = os.path.abspath(os.path.dirname(evv4esm.__file__))

logger = logging.getLogger(__name__)

NUMBER_INITIAL_CONDITIONS = 6
PERTURBATIONS = OrderedDict([('woprt', 0.0),
                             ('posprt', 1.0e-14),
                             ('negprt', -1.0e-14),
                             ])
FCLD_NC = 'cam.h0.cloud.nc'
INIT_COND_FILE_TEMPLATE = \
    "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.{}.{}.0002-{:02d}-01-00000.nc"
# FIXME: should 'cam' be 'atm' now?
INSTANCE_FILE_TEMPLATE = '{}cam_{:04d}.h0.0001-01-01-00000{}.nc'


class PGN(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the PGN test
        """
        super(PGN, self).__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        ninst = NUMBER_INITIAL_CONDITIONS * len(PERTURBATIONS)
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

        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

        logger.debug("PGN_INFO: Updating user_nl_* files")

        csmdata_root = self._case.get_value("DIN_LOC_ROOT")
        csmdata_atm = os.path.join(csmdata_root, "atm/cam/inic/homme/ne4_v1_init")
        csmdata_lnd = os.path.join(csmdata_root, "lnd/clm2/initdata/ne4_v1_init/b58d55680")

        iinst = 1
        for icond in range(1, NUMBER_INITIAL_CONDITIONS + 1):
            fatm_in = os.path.join(csmdata_atm, INIT_COND_FILE_TEMPLATE.format('atm', 'i', icond))
            flnd_in = os.path.join(csmdata_lnd, INIT_COND_FILE_TEMPLATE.format('clm2', 'r', icond))
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

        return list(map(str.strip, var_list))

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
                    "ninit": NUMBER_INITIAL_CONDITIONS,
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
                if evv_elem['Type'] == 'ValSummary' \
                        and evv_elem['TableTitle'] == 'Perturbation growth test':
                    if evv_elem['Data'][test_name]['']['Test status'].lower() == 'pass':
                        self._test_status.set_status(CIME.test_status.BASELINE_PHASE,
                                                     CIME.test_status.TEST_PASS_STATUS)
                        break

    def run_phase(self):
        logger.debug("PGN_INFO: RUN PHASE")

        self.run_indv()

        # Here were are in case directory, we need to go to the run directory
        # and rename files
        rundir = self._case.get_value("RUNDIR")
        casename = self._case.get_value("CASE")
        logger.debug("PGN_INFO: Case name is:{}".format(casename))

        for icond in range(NUMBER_INITIAL_CONDITIONS):
            for iprt, (prt, prt_name) in enumerate(PERTURBATIONS.items()):
                iinst = pg._sub2instance(icond, iprt, len(PERTURBATIONS))
                fname = os.path.join(rundir, INSTANCE_FILE_TEMPLATE.format(casename + '.', iinst, ''))
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

        logger.debug("PGN_INFO: RUN PHASE ENDS")

    def _generate_baseline(self):
        super(PGN, self)._generate_baseline()

        base_root = self._case.get_value("BASELINE_ROOT")
        base_gen = self._case.get_value("BASEGEN_CASE")
        base_dir = os.path.join(base_root, base_gen)

        rundir = self._case.get_value("RUNDIR")
        casename = self._case.get_value("CASE")

        var_list = self.get_var_list()
        nvar = len(var_list)
        nprt = len(PERTURBATIONS)
        rmse_prototype = {}
        for icond in range(NUMBER_INITIAL_CONDITIONS):
            prt_rmse = {}
            for iprt, prt_name in enumerate(PERTURBATIONS):
                if prt_name == 'woprt':
                    continue
                iinst_ctrl = pg._sub2instance(icond, 0, nprt)
                ifile_ctrl = os.path.join(rundir,
                                          INSTANCE_FILE_TEMPLATE.format(casename + '.', iinst_ctrl, '_woprt'))

                iinst_test = pg._sub2instance(icond, iprt, nprt)
                ifile_test = os.path.join(rundir,
                                          INSTANCE_FILE_TEMPLATE.format(casename + '.', iinst_test, '_' + prt_name))

                prt_rmse[prt_name] = pg.variables_rmse(ifile_test, ifile_ctrl, var_list, 't_')
            rmse_prototype[icond] = pd.concat(prt_rmse)
        rmse = pd.concat(rmse_prototype)
        cld_rmse = np.reshape(rmse.RMSE.values, (NUMBER_INITIAL_CONDITIONS, nprt - 1, nvar))

        pg.rmse_writer(os.path.join(rundir, FCLD_NC),
                       cld_rmse, PERTURBATIONS.keys(), var_list, INIT_COND_FILE_TEMPLATE)

        logger.debug("PGN_INFO:copy:{} to {}".format(FCLD_NC, base_dir))
        shutil.copy(os.path.join(rundir, FCLD_NC), base_dir)
