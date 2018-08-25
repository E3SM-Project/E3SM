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

import os
import glob
import shutil
import numpy as np

from shutil import copyfile
from netCDF4 import Dataset
from math import sqrt
from sklearn.metrics import mean_squared_error

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.build import post_build
from CIME.hist_utils import rename_all_hist_files
from CIME.utils import expect

#Logic for SLR ensemble runs:
# We need two inputs:
# 1. Number of inic cond files
# 2. perturbations (e.g. currently we have no prt, pos prt and neg prt)
# Based of the above, we compute number of instances to have
# We change the user_nl* files accordingly (during build time) 
# and rename history files (after run).


#TODO: 
#1. The loop to change user_nl* files is run twice as we call build again after changing ninst, which changes ntasks
#2. Do we want to remove user_nl_cam, user_nl_clm etc. files as they are not used in the simulation?
#3. The loop for modifying namelists and renaming history files should be EXACTLY same. Try devising a function
#   with these loops and calling another function within the loops for either modifying namelists (for build pahse)
#   or renaming history files (run phase) so that these loops stay the same for both modifying namelists
#   and renaming history files.



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

        #-----------------------------------------------------
        #Compute number of instances:
        #------------------------------------------------------

        #number of inital conditions
        ninit_cond = 2 #12
        
        #perturbations for runs        
        prt        = [0.0, 1.0e-14, -1.0e-14]

        #------------------------------------------------------
        #Number of instances:
        #--------------------
        #Comput it from number of initial conditions to use and 
        #number of perturbation ensemble members to use
        #------------------------------------------------------
        nprt       = len(prt)
        ninst      = ninit_cond * nprt

        logging.warn('SLR_INFO: number of instance: '+str(ninst))
        
        #------------------------------------------------------
        #Fake Build:
        #-----------
        #(for debugging only) Building the model can take 
        #significant amount of time. Setting fake_bld to True 
        #can savethat time
        #------------------------------------------------------
        fake_bld = False

        #Find number of instance in the default setup
        default_ninst = self._case.get_value("NINST_ATM")
        
        #Sanity check: see if NINST is same for all model components, otherwise exit with error
        for comp in ['OCN','WAV','GLC','ICE','ROF','LND']:
            iinst = self._case.get_value("NINST_%s"%comp)
            expect(default_ninst == iinst, "ERROR: component "+str(comp)+" NINST("+str(iinst)+")"
                   " is different from component ATM NINST("+str(default_ninst)+")")

        #------------------------------------------------------
        #Setup multi-instances for model components:
        #-------------------------------------------
        #Set the model for multi-instance ONLY if NINST == 1 
        #for all model components. This is because, for 
        #NINST > 1 (e.g. rebuilding an old case) the following 
        #loop will increase the ntasks to a multiple of ninst 
        #(requiring a clean build again). We hit this issue if 
        #we launch ./case.build in the case directory of SLR 
        #test
        #------------------------------------------------------

        if(default_ninst == 1): #if multi-instance is not already set
            # Only want this to happen once. It will impact the sharedlib build
            # so it has to happen here.
            if not model_only:
                # Lay all of the components out concurrently
                logging.warn("SLR_INFO: Updating NINST for multi-instance in env_mach_pes.xml")
                for comp in ['ATM','OCN','WAV','GLC','ICE','ROF','LND']:
                    ntasks = self._case.get_value("NTASKS_%s"%comp)
                    self._case.set_value("ROOTPE_%s"%comp, 0)
                    self._case.set_value("NINST_%s"%comp,  ninst)
                    self._case.set_value("NTASKS_%s"%comp, ntasks*ninst)

                self._case.set_value("ROOTPE_CPL",0)
                self._case.set_value("NTASKS_CPL",ntasks*ninst)
                self._case.flush()

                case_setup(self._case, test_mode=False, reset=True)

        #Faking a bld can save the time code spend in building the model components
        if fake_bld:
            logging.warn("SLR_INFO: FAKE Build")
            if (not sharedlib_only):
                post_build(self._case, [])
        else:
            # Build exectuable with multiple instances
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)


        #----------------------------------------------------------------
        # Run-time settings:
        # Do this already in build_phase so that we can check the xml and
        # namelist files before job starts running.
        #----------------------------------------------------------------

        #Prepare paths for namelist file

        # generate paths/file names to get files to set initial conditons
        csmdata_root = self._case.get_value("DIN_LOC_ROOT")
        csmdata_atm  = csmdata_root+"atm/cam/inic/homme/"
        csmdata_lnd  = csmdata_root+"lnd/clm2/initdata/"
        file_pref_atm = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.cam.i.0002-"
        file_pref_lnd = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.clm2.r.0002-"

        file_suf_atm = "-01-00000.nc"
        file_suf_lnd = "-01-00000.nc"

        def write_nl():
            logging.warn("SLR_INFO: Updating user_nl_* files")
            iinst = 1
            for icond in range(ninit_cond):
                for iprt in prt:
                    with open('user_nl_cam_'+str(iinst).zfill(4), 'w') as atmnlfile, \
                            open('user_nl_clm_'+str(iinst).zfill(4), 'w') as lndnlfile:
                    
                        #atm/lnd intitial conditions                   
                        inst_label_2digits = str(iinst).zfill(2)
                    
                        #initial condition files to use for atm and land
                        atmnlfile.write("ncdata  = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
                        lndnlfile.write("finidat = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")
                        
                        #uncomment the following when there files are on SVN server
                        #atmnlfile.write("ncdata  = '"+ csmdata_atm + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
                        #lndnlfile.write("finidat = '"+ csmdata_lnd + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")

                    
                        #atm model output
                        atmnlfile.write("avgflag_pertape = 'I' \n")
                        atmnlfile.write("nhtfrq = 1 \n")
                        atmnlfile.write("mfilt  = 1  \n")
                        atmnlfile.write("ndens  = 1  \n")
                        #atmnlfile.write("empty_htapes = .true. \n")
                        #atmnlfile.write("fincl1 = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC' \n")
                        atmnlfile.write("pergro_mods  = .true. \n")
                        atmnlfile.write("pergro_test_active = .true. \n")
                
                        if(iprt != 0.0):
                            atmnlfile.write("pertlim = "+str(iprt)+" \n")
                        
                        iinst = iinst + 1

        write_nl()
        #=================================================================
        #Settings common to all instances
        #=================================================================

        #Coupler settings which are applicable to ALL the instances (as there is only one coupler for all instances)
        self._case.set_value("STOP_N",     "1")
        self._case.set_value("STOP_OPTION","nsteps")

    #=================================================================
    # run_phase
    #=================================================================
    def run_phase(self):
        print('SLR_INFO: RUN PHASE')
        
        self.run_indv()
        #Here were are in case directory, we need to go to the run directory and rename files
        cwd = os.getcwd()
        
        
