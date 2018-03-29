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
        # The default NTASKS seem to be much lower than we'd need, so scale it
        scale_ntasks = 1

        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        if not model_only:
            logging.warn('Starting to build multi-instance exe')
            for comp in self._case.get_values("COMP_CLASSES"):
                self._case.set_value('ROOTPE_{}'.format(comp), 0)
                self._case.set_value('NTHRDS_{}'.format(comp), 1)

                ntasks = self._case.get_value("NTASKS_{}".format(comp))
                print('MVK:NTASKS_{}: {}'.format(comp, ntasks))
                self._case.set_value('NTASKS_{}'.format(comp), ntasks*scale_ntasks*ninst)
                if comp != 'CPL':
                    self._case.set_value('NINST_{}'.format(comp), ninst)

            self._case.flush()

            case_setup(self._case, test_mode=False, reset=True)


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
                nl_atm_file.write("fincl2 = 'TS:A','TSMN:M','TSMX:X','PRECT:A','PRECTMX:X'\n")
                nl_atm_file.write('mfilt = 1,30\n')
                nl_atm_file.write('nhtfrq = 0,-24\n')

                # nl_atm_file.write('cldfrc_rhminl = 0.9100\n')
                nl_atm_file.write('zmconv_c0_ocn = 0.0074D0\n')
                nl_atm_file.write('dtime = 5400\n')

                # nl_atm_file.write("rsf_file = '{}/atm/waccm/phot/RSF_GT200nm_v3.0_c080416.nc'\n".format(data_root))
                # nl_atm_file.write("bnd_topo = '{}/atm/cam/topo/USGS-gtopo30_1.9x2.5_smooth500-50_ne16np4_c050602.nc'\n".format(data_root))

                # nl_atm_file.write('phys_loadbalance = 0\n')
                nl_atm_file.write('new_random = .true.\n')
                nl_atm_file.write('pertlim = 1.0e-10\n')
                nl_atm_file.write('seed_custom = {}\n'.format(iinst))

            # with open('user_nl_clm_{:04d}'.format(iinst), 'w') as nl_clm_file:
            #     # nl_clm_file.write("finidat = '{}/lnd/clm2/initdata_map/clmi.I1850CLM45.ne30_g16.1155ba0.clm2.r.nc'\n".format(data_root))
            #     nl_clm_file.write("finidat = '{}/lnd/clm2/initdata_map/clmi.I1850CLM45.ne16_qu240.d02c96d.clm2.r.0021-01-01-00000.nc'\n".format(data_root))


    def run_phase(self):

        self.run_indv()
