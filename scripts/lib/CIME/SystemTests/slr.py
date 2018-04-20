"""
Solution reproducibility test - common parts shared by
different test methods. The CESM/ACME model's
multi-instance capability is used to conduct an ensemble
of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.

Different solution reproducibility test methods use
different namelist settings and postprocessing.
We will decide later whether separate test types
will be created for those different methods
(pergro, time step convergence, statistical methods)
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.build import post_build
#import CIME.utils
#from CIME.check_lockedfiles import *

logger = logging.getLogger(__name__)

class SLR(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SLR test
        """
        SystemTestsCommon.__init__(self, case)

    #=================================================================
    # Compile model with multiple instances
    #=================================================================
    def build_phase(self, sharedlib_only=False, model_only=False):
        #total number of initial conditions
        n_inic_cond = 12

        #For each initial condition, we need to run 3 simulations:
        #1. Default
        #2. With Positive Perturbation
        #3. With Negative Perturb

        #ninst = n_inic_cond * 3
        ninst =  3

        #BSINGH: For debugging only - Building the model can take
        #significant amount of time. Setting fake_bld to True can save
        #that time
        fake_bld = False

        # Build exectuable with multiple instances
        # Lay all of the components out concurrently

        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        if not model_only:
            logging.warn("Starting to build multi-instance exe")
            for comp in ['ATM','OCN','WAV','GLC','ICE','ROF','LND']:
                ntasks = self._case.get_value("NTASKS_%s"%comp)
                self._case.set_value("ROOTPE_%s"%comp, 0)
                self._case.set_value("NINST_%s"%comp,  ninst)
                self._case.set_value("NTASKS_%s"%comp, ntasks*ninst)

            self._case.set_value("ROOTPE_CPL",0)
            self._case.set_value("NTASKS_CPL",ntasks*ninst)
            self._case.flush()

            case_setup(self._case, test_mode=False, reset=True)

        #BSINGH: Faking a bld can save the time code spend in building the model components
        if fake_bld:
            if (not sharedlib_only):
                post_build(self._case, [])
        else:
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)


        #=================================================================
        # Run-time settings.
        # Do this already in build_phase so that we can check the xml and
        # namelist files before job starts running.
        #=================================================================

        #=================================================================
        #Settings common to all instances
        #=================================================================

        #Coupler settings which are applicable to ALL the instances (as there is only one coupler for all instances)
        #self._case.set_value("STOP_N",     "1")
        #self._case.set_value("STOP_OPTION","nsteps")

        #Prepare paths for namelist file

        # generate paths/file names to get files to set initial conditons
        #csmdata_root = self._case.get_value("DIN_LOC_ROOT")
        #csmdata_atm  = csmdata_root+"atm/cam/inic/homme/"
        #csmdata_lnd  = csmdata_root+"lnd/clm2/initdata/"
        #file_pref_atm = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.cam.i.0002-"
        #file_pref_lnd = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.clm2.r.0002-"

        #file_suf_atm = "-01-00000.nc"
        #file_suf_lnd = "-01-00000.nc"

        #def write_nl(iinst,prt=0):

        #     with open('user_nl_cam_'+str(iinst).zfill(4), 'w') as atmnlfile, \
        #          open('user_nl_clm_'+str(iinst).zfill(4), 'w') as lndnlfile:

        #           atm/lnd intitial conditions

        #          inst_label_2digits = str(iinst).zfill(2)

        #          atmnlfile.write("ncdata  = '"+ csmdata_atm + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
        #          lndnlfile.write("finidat = '"+ csmdata_lnd + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")

        #           initial condition files to use for atm and land
        #          atmnlfile.write("ncdata  = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
        #          lndnlfile.write("finidat = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")

        #           atm model output
        #          atmnlfile.write("avgflag_pertape = 'I' \n")
        #          atmnlfile.write("nhtfrq = 1 \n")
        #          atmnlfile.write("mfilt  = 1  \n")
        #          atmnlfile.write("ndens  = 1  \n")
        #          atmnlfile.write("empty_htapes = .true. \n")
        #          atmnlfile.write("fincl1 = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC' \n")
        #          atmnlfile.write("pergro_mods  = .true. \n")
        #          atmnlfile.write("pergro_test_active = .true. \n")

        #          if(prt <> 0):
        #              atmnlfile.write(" pertlim = "+str(prt)+" \n")

        # for iinst in range(1, ninst+1,3):
        #     write_nl(iinst,prt=0)


    #=================================================================
    # run_phase
    #=================================================================
    def run_phase(self):

        self.run_indv()
