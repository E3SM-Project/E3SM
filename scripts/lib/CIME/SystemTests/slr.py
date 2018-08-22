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



logger = logging.getLogger(__name__)

class SLR(SystemTestsCommon):

    #Number of instances are total number of initial conditions we have
    global ninst
    ninst = 6

    def __init__(self, case):
        """
        initialize an object interface to the SLR test
        """
        SystemTestsCommon.__init__(self, case)

    #=================================================================
    # Compile model with multiple instances
    #=================================================================
    def build_phase(self, sharedlib_only=False, model_only=False):

        #BSINGH: For debugging only - Building the model can take
        #significant amount of time. Setting fake_bld to True can save
        #that time
        fake_bld = False

        #Find number of instance
        default_ninst = self._case.get_value("NINST_ATM")
        
        #Sanity check: see if NINST is same for all model components, otherwise exit with error
        for comp in ['OCN','WAV','GLC','ICE','ROF','LND']:
            iinst = self._case.get_value("NINST_%s"%comp)
            expect(default_ninst == iinst, "ERROR: component "+str(comp)+" NINST("+str(iinst)+")"
                   " is different from componant ATM NINST("+str(default_ninst)+")")

        #Set the model for multi-instance ONLY if NINST == 1 for all model components.
        #This is because, for NINST > 1 (e.g. rebuilding an old case) the following loop
        #will increase the ntasks to a multiple of ninst (requiring a clean build again).

        if(default_ninst == 1):
            # Only want this to happen once. It will impact the sharedlib build
            # so it has to happen here.
            if not model_only:
                # Lay all of the components out concurrently
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
            # Build exectuable with multiple instances
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
        self._case.set_value("STOP_N",     "1")
        self._case.set_value("STOP_OPTION","nsteps")

    #=================================================================
    # run_phase
    #=================================================================
    def run_phase(self):
        
        #Prepare paths for namelist file

        # generate paths/file names to get files to set initial conditons
        csmdata_root = self._case.get_value("DIN_LOC_ROOT")
        csmdata_atm  = csmdata_root+"atm/cam/inic/homme/"
        csmdata_lnd  = csmdata_root+"lnd/clm2/initdata/"
        file_pref_atm = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.cam.i.0002-"
        file_pref_lnd = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.clm2.r.0002-"

        file_suf_atm = "-01-00000.nc"
        file_suf_lnd = "-01-00000.nc"

        def copy_nl(directory):
            if os.path.exists(directory):
                shutil.rmtree(directory)
                os.makedirs(directory)
            else:
                os.makedirs(directory)
                
            user_nl_files = glob.glob('user_nl_*')
            for unl_file in user_nl_files:
                shutil.copy(unl_file,directory)


        def write_nl(directory,prt=0):
            for iinst in range(1, ninst+1):
                with open('user_nl_cam_'+str(iinst).zfill(4), 'w') as atmnlfile, \
                        open('user_nl_clm_'+str(iinst).zfill(4), 'w') as lndnlfile:
                    
                    #atm/lnd intitial conditions                   
                    inst_label_2digits = str(iinst).zfill(2)
                    
                    atmnlfile.write("ncdata  = '"+ csmdata_atm + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
                    lndnlfile.write("finidat = '"+ csmdata_lnd + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")
                    
                    #initial condition files to use for atm and land
                    atmnlfile.write("ncdata  = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_atm + inst_label_2digits + file_suf_atm+"' \n")
                    lndnlfile.write("finidat = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_lnd + inst_label_2digits + file_suf_lnd+"' \n")
                    
                    #atm model output
                    atmnlfile.write("avgflag_pertape = 'I' \n")
                    atmnlfile.write("nhtfrq = 1 \n")
                    atmnlfile.write("mfilt  = 1  \n")
                    atmnlfile.write("ndens  = 1  \n")
                    #atmnlfile.write("empty_htapes = .true. \n")
                    #atmnlfile.write("fincl1 = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC' \n")
                    atmnlfile.write("pergro_mods  = .true. \n")
                    atmnlfile.write("pergro_test_active = .true. \n")
                
                    if(prt != 0):
                        atmnlfile.write("pertlim = "+str(prt)+" \n")

            #copy namelist files for keeping record
            copy_nl(directory)

        #We need to lauch 3 sets:
        #1. Default
        #2. With Positive Perturbation
        #3. With Negative Perturb

        #=================================================================
        #namelist changes for the first set:
        #=================================================================

        write_nl("default_user_nl")
        self.run_indv()
        rename_all_hist_files( self._case, suffix="def" )

        write_nl("pos_user_nl",prt=1.e-14)
        self.run_indv()
        rename_all_hist_files( self._case, suffix="pos" )

        write_nl("neg_user_nl",prt=-1.e-14)
        self.run_indv()        
        rename_all_hist_files( self._case, suffix="neg" )

        
        def rmse_diff_var(ifile_test, ifile_cntl,  var_list, var_suffix):
            print( os.getcwd())
            print( ifile_test)
            print(ifile_cntl)
            ftest = Dataset(ifile_test, mode='r')
            fcntl = Dataset(ifile_cntl, mode='r')
            
            diff = np.zeros(shape=(len(var_list)))#
            i = 0
            
            is_se = (len(ftest.variables[var_suffix+var_list[i]].dimensions)==3) # see if it is SE grid
            nz = 1
                            
            for ivar in var_list:
                var = var_suffix+ivar
                if var in ftest.variables:
                    vtest    = ftest.variables[var.strip()][0,...] # first dimention is time (=0)
                    vcntl    = fcntl.variables[var.strip()][0,...] # first dimention is time (=0)
                    
                    #reshape for RMSE
                    if(is_se ):
                        nx, ny = vtest.shape #shape will be same for both arrays
                    else:
                        nx, ny, nz = vtest.shape #shape will be same for both arrays
                    diff[i]  = sqrt(mean_squared_error(vtest.reshape((nx,ny*nz)), vcntl.reshape((nx,ny*nz))))
                    print(str(var)+":"+str(diff[i]))
                    i += 1
                    vtest = None
                    vcntl = None
                    var   = None
                else:
                    print (var+' not in file')

            ftest.close()
            fcntl.close()
            return diff

        #post processing
        #cd into run directory
        #os.chdir(self._case.get_value("RUNDIR"))
        
        #with open('pergro_ptend_names.txt', 'r') as fvar:
        #    var_list = fvar.readlines()
        #fvar.close()
        #var_list = map(str.strip,var_list)

        
        ##def_cam_hist = glob.glob("*.cam*.h0.*.def")
        ##pos_cam_hist = glob.glob("*.cam*.h0.*.pos")

        ##def_cam_hist = glob.glob("*.cam_0001*.h0.*-00000*.def")
        ##pos_cam_hist = glob.glob("*.cam_0001*.h0.*-00000*.pos")

        #ts = "07200"
        #for iinst in range(1, 3): #ninst+1):
        #    cam_inst_str = "cam_"+str(iinst).zfill(4)
        #    def_cam_hist = glob.glob("*"+cam_inst_str+"*.h0.*-"+ts+"*.def")
        #    pos_cam_hist = glob.glob("*"+cam_inst_str+"*.h0.*-"+ts+"*.pos")
            
        #    rmse = rmse_diff_var(str(pos_cam_hist[0]), str(def_cam_hist[0]), var_list, 't_')

            
        
        
