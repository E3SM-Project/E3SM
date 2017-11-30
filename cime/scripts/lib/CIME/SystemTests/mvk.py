"""
Multivariate test for climate reproducibility using the Kolmogrov-Smirnov (K-S)
test and based on The CESM/E3SM model's multi-instance capability is used to
conduct an ensemble of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case_setup import case_setup
from CIME.build import post_build
# import CIME.utils
# from CIME.check_lockedfiles import *

logger = logging.getLogger(__name__)


class MVK(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the MVK test
        """
        SystemTestsCommon.__init__(self, case)


    def build_phase(self, sharedlib_only=False, model_only=False):

        # Build executable with multiple instances
        ninst = 30

        # BSINGH: For debugging only - Building the model can take
        # significant amount of time. Setting fake_bld to True can save
        # that time
        fake_bld = False

        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        if not model_only:
            logging.warn('Starting to build multi-instance exe')
            for comp in ['ATM', 'OCN', 'WAV', 'GLC', 'ICE', 'ROF', 'LND']:
                ntasks = self._case.get_value("NTASKS_{}".format(comp))
                self._case.set_value('ROOTPE_{}'.format(comp), 0)
                self._case.set_value('NINST_{}'.format(comp), ninst)
                self._case.set_value('NTASKS_{}'.format(comp), ntasks*ninst)

            self._case.set_value('ROOTPE_CPL', 0)
            self._case.set_value('NTASKS_CPL', ntasks*ninst)
            self._case.flush()

            case_setup(self._case, test_mode=False, reset=True)

        # BSINGH: Faking a bld can save the time code spend in building the model components
        if fake_bld:
            if not sharedlib_only:
                post_build(self._case, [])
        else:
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

        # =================================================================
        # Run-time settings.
        # Do this already in build_phase so that we can check the xml and
        # namelist files before job starts running.
        # =================================================================

        # namelist specifications for each instance
        data_root = self._case.get_value("DIN_LOC_ROOT")
        for iinst in range(1, ninst+1):
            with open('user_nl_cam_{:04d}'.format(iinst), 'w') as nl_atm_file:
                nl_atm_file.write("fincl2 = 'TS:A','TSMN:M','TSMX:X','PRECT:A','PRECTMX:X'")
                nl_atm_file.write('mfilt = 1,30')
                nl_atm_file.write('nhtfrq = 0,-24')

                nl_atm_file.write('cldfrc_rhminl = 0.9100')
                nl_atm_file.write('zmconv_c0_ocn = 0.0110')

                nl_atm_file.write("rsf_file = '{}/atm/waccm/phot/RSF_GT200nm_v3.0_c080416.nc'".format(data_root))
                nl_atm_file.write("bnd_topo = '{}/atm/cam/topo/USGS-gtopo30_1.9x2.5_smooth500-50_ne16np4_c050602.nc'".format(data_root))

                nl_atm_file.write('phys_loadbalance = 0')
                nl_atm_file.write('new_random = .true.')
                nl_atm_file.write('pertlim = 1.0e-10')
                nl_atm_file.write('seed_custom = {}'.format(iinst))

            with open('user_nl_clm_{:04d}'.format(iinst), 'w') as nl_clm_file:
                nl_clm_file.write("finidat = '/lustre/atlas/proj-shared/cli106/salil/archive/A_F.1850c5.ne16_g37.l2/rest/0100-11-01-00000/A_F.1850c5.ne16_g37.l2.clm2.r.0100-11-01-00000.nc'")


# ------
# Notes:
# ------
# Follow setup of cases in: /autofs/nccs-svm1_home1/salil/acme_cases/ne16_F1850C5_master/
# User_nl_* changes:
#   CAM:
#       fincl2 = 'TS:A','TSMN:M','TSMX:X','PRECT:A','PRECTMX:X'
#       mfilt = 1,30
#       nhtfrq = 0,-24
#
#       cldfrc_rhminl = 0.9100
#       zmconv_c0_ocn = 0.0110
#
#       rsf_file = '/lustre/atlas1/cli900/world-shared/cesm/inputdata/atm/waccm/phot/RSF_GT200nm_v3.0_c080416.nc'
#       bnd_topo = '/lustre/atlas1/cli900/world-shared/cesm/inputdata/atm/cam/topo/USGS-gtopo30_1.9x2.5_smooth500-50_ne16np4_c050602.nc'
#
#       phys_loadbalance = 0
#       new_random = .true.
#       pertlim = 1.e-10
#       seed_custom = 2
#
#   CLM:
#       finidat = '/lustre/atlas/proj-shared/cli106/salil/archive/A_F.1850c5.ne16_g37.l2/rest/0100-11-01-00000/A_F.1850c5.ne16_g37.l2.clm2.r.0100-11-01-00000.nc'

# Original case setup:
# ./create_newcase -mach titan -compiler pgi -res ne16_g37 -compset F1850C5 -case /ccs/home/salil/acme_cases/F1850C5_08162017_master -project cli106
