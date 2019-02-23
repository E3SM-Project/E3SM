"""
Solution reproducibility test based on time-step convergence
The CESM/ACME model's
multi-instance capability is used to conduct an ensemble
of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.build import post_build
from CIME.hist_utils import rename_all_hist_files
#import CIME.utils
#from CIME.check_lockedfiles import *

logger = logging.getLogger(__name__)

class TSC(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the TSC test
        """
        SystemTestsCommon.__init__(self, case)

    #=================================================================
    # Compile model with multiple instances
    #=================================================================
    def build_phase(self, sharedlib_only=False, model_only=False):

        ninst = 12

        #BSINGH: For debugging only - Building the model can take
        #significant amount of time. Setting fake_bld to True can save
        #that time
        fake_bld = False

        # Build exectuable with multiple instances
        # (in the development phase we use 3)
        # Lay all of the components out concurrently

        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        if not model_only:
            logging.warning("Starting to build multi-instance exe")
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
    # Conduct one multi-instance run with specified time step size.
    # The complete run_phase (defined below) will contain two calls
    # of this function using different dtime.
    #=================================================================
    def _run_with_specified_dtime(self,dtime=2):

        # dtime is model time step size in seconds. Default is set to 2.

        # Change the coupling frequency accordingly
        ncpl  = 86400/dtime
        self._case.set_value("ATM_NCPL",ncpl)

        # Simulation length = 10 minutes
        nsteps = 600/dtime
        self._case.set_value("STOP_N",     nsteps)
        self._case.set_value("STOP_OPTION","nsteps")

        # Write output every 10 seconds
        nhtfrq = 10/dtime

        # generate paths/file names for initial conditons
        csmdata_root = self._case.get_value("DIN_LOC_ROOT")
        csmdata_atm  = csmdata_root+"/atm/cam/inic/homme/ne4_v1_init/"
        csmdata_lnd  = csmdata_root+"/lnd/clm2/initdata/ne4_v1_init/b58d55680/"
        file_pref_atm = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.cam.i.0002-"
        file_pref_lnd = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.clm2.r.0002-"

        file_suf_atm = "-01-00000.nc"
        file_suf_lnd = "-01-00000.nc"

        # user_nl_... files are modified for every instance.
        # The number instances specified here has to be consistent with the
        # value set in build_phase.
        ninst = 12
        for iinst in range(1, ninst+1):

            with open('user_nl_cam_'+str(iinst).zfill(4), 'w') as atmnlfile, \
                 open('user_nl_clm_'+str(iinst).zfill(4), 'w') as lndnlfile:

                # atm/lnd intitial conditions

                inst_label_2digits = str(iinst).zfill(2)

                atmnlfile.write("ncdata  = '"+ csmdata_atm + "/" + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
                lndnlfile.write("finidat = '"+ csmdata_lnd + "/" + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")

                # for initial testing on constance@pnnl
                #atmnlfile.write("ncdata  = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
                #lndnlfile.write("finidat = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")

                # time step sizes

                atmnlfile.write("dtime = "+str(dtime)+" \n")
                atmnlfile.write("iradsw = 2 \n")
                atmnlfile.write("iradlw = 2 \n")

                lndnlfile.write("dtime = "+str(dtime)+" \n")

                # atm model output

                atmnlfile.write("avgflag_pertape = 'I' \n")
                atmnlfile.write("nhtfrq = "+str(nhtfrq)+" \n")
                atmnlfile.write("mfilt  = 1  \n")
                atmnlfile.write("ndens  = 1  \n")
                atmnlfile.write("empty_htapes = .true. \n")
                atmnlfile.write("fincl1 = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC' \n")

        # Force rebuild namelists
        self._skip_pnl = False

        # Namelist settings done. Now run the model.
        self.run_indv()

        # Append "DT000x" to the history file names
        rename_all_hist_files( self._case, suffix="DT"+str(dtime).zfill(4) )


    #=================================================================
    # run_phase
    #=================================================================
    def run_phase(self):

        # Run simulations with 2-s time step
        self._run_with_specified_dtime(dtime=2)

        # Run simulations with 1-s time step. This is NOT needed
        # for testing new code or new computing environment.
        if self._case.get_value("GENERATE_BASELINE"):
            self._run_with_specified_dtime(dtime=1)
