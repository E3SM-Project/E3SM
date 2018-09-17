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
# In baseline generation step, we compute cloud
# In comparison step, we compare against baseline cloud to know pass/fail


#TODO: 
#1. The loop to change user_nl* files is run twice as we call build again after changing ninst, which changes ntasks
#2. Do we want to remove user_nl_cam, user_nl_clm etc. files as they are not used in the simulation?
#3. The loop for modifying namelists and renaming history files should be EXACTLY same. Try devising a function
#   with these loops and calling another function within the loops for either modifying namelists (for build pahse)
#   or renaming history files (run phase) so that these loops stay the same for both modifying namelists
#   and renaming history files.
#4. For generating baselines: rename files and generate cloud. Store that cloud so that we can re-use it when we run it for comparison
#5. For comparison: Compute rmse of current netcdf files w.r.t baselines and give pass/fail
#6. change code so that pergro_ptend_names.txt is not generated if it is already there or only one instance writes this file...





logger = logging.getLogger(__name__)

#number of inital conditions
ninit_cond = 1 #12
        
#perturbations for runs
prt        = [0.0, 1.0e-14, -1.0e-14]
prtstr     = ['woprt','posprt','negprt']


index = False #DEBUGGING: print out max rmse diffs if set to True

class SLR(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SLR test
        """

        #perturbation values and number of perturbation strings should be same
        expect(len(prt)== len(prtstr),"Number of perturbation values ("+str(len(prt))+") are NOT " 
               "equal to number of perturbation strings("+ str(len(prtstr))+")")
        SystemTestsCommon.__init__(self, case)

    #=================================================================
    # Compile model with multiple instances
    #=================================================================
    def build_phase(self, sharedlib_only=False, model_only=False):

        #------------------------------------------------------
        #Compute number of instances:
        #------------------------------------------------------

        #------------------------------------------------------
        #Number of instances:
        #~~~~~~~~~~~~~~~~~~~
        #Comput it from number of initial conditions to use and 
        #number of perturbation ensemble members to use
        #------------------------------------------------------
        nprt       = len(prt)
        ninst      = ninit_cond * nprt

        logging.warn('SLR_INFO: number of instance: '+str(ninst))
        
        #------------------------------------------------------
        #Fake Build:
        #~~~~~~~~~~
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
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        #~~~~~~~~~~~~~~~~~~
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

        
        logging.warn("SLR_INFO: Updating user_nl_* files")
        iinst = 1
        for icond in range(ninit_cond):
            icond_label_2digits = str(icond+1).zfill(2)
            for iprt in prt:
                with open('user_nl_cam_'+str(iinst).zfill(4), 'w') as atmnlfile, \
                        open('user_nl_clm_'+str(iinst).zfill(4), 'w') as lndnlfile:
                    
                    #atm/lnd intitial conditions                   
                    
                    #initial condition files to use for atm and land
                    atmnlfile.write("ncdata  = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_atm + icond_label_2digits + file_suf_atm+"' \n")
                    lndnlfile.write("finidat = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + file_pref_lnd + icond_label_2digits + file_suf_lnd+"' \n")
                    
                    #uncomment the following when there files are on SVN server
                    #atmnlfile.write("ncdata  = '"+ csmdata_atm + file_pref_atm + icond_label_2digits + file_suf_atm+"' \n")
                    #lndnlfile.write("finidat = '"+ csmdata_lnd + file_pref_lnd + icond_label_2digits + file_suf_lnd+"' \n")
                    
                    
                    #atm model output
                    atmnlfile.write("avgflag_pertape = 'I' \n")
                    atmnlfile.write("nhtfrq = 1 \n")
                    atmnlfile.write("mfilt  = 2  \n")
                    atmnlfile.write("ndens  = 1  \n")
                    #atmnlfile.write("empty_htapes = .true. \n")
                    #atmnlfile.write("fincl1 = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC' \n")
                    atmnlfile.write("pergro_mods  = .true. \n")
                    atmnlfile.write("pergro_test_active = .true. \n")
                    #atmnlfile.write("phys_debug_lat = 41.3495891345")
                    #atmnlfile.write("phys_debug_lon = 45.0" )
                
                    if(iprt != 0.0):
                        atmnlfile.write("pertlim = "+str(iprt)+" \n")
                        
                    iinst += 1

        #=================================================================
        #Settings common to all instances
        #=================================================================

        #Coupler settings which are applicable to ALL the instances (as there is only one coupler for all instances)
        self._case.set_value("STOP_N",     "1")
        self._case.set_value("STOP_OPTION","nsteps")

    def get_fname_wo_ext(self,rundir,casename,iinst):
        #form file name withOUT extension
        tstmp    = '00000'
        date_str = '.h0.0001-01-01-'
        cam_str  = 'cam_'
        if(casename != "" ):
            cam_str = '.'+cam_str
        return os.path.join(rundir,casename+cam_str+str(iinst).zfill(4)+date_str+tstmp)

    def get_var_list(self):
        rundir    = self._case.get_value("RUNDIR")
        casename  = self._case.get_value("CASE")    
        prg_fname = 'pergro_ptend_names.txt'
        var_file = os.path.join(rundir,prg_fname)
        expect(os.path.isfile(var_file),"File "+prg_fname+" does not exist in: "+rundir)

        with open(var_file, 'r') as fvar:
            var_list = fvar.readlines()
            fvar.close()
            
        return  map(str.strip,var_list)
        
    def rmse_var(self,ifile_test, ifile_cntl,  var_list, var_suffix ):
        from netCDF4 import Dataset
        from math import sqrt
        from sklearn.metrics import mean_squared_error
        import numpy as np
        
        
        #------------------ARGS-------------------------
        #ifile_test: path of test file
        #ifile_cntl: path of control file
        #var_list  : List of all variables
        #var_suffix: Suffix for var_list (e.g. T_, S_ QV_ etc.)
        #-----------------------------------------------

        # See if the files exists or not....
        expect(os.path.isfile(ifile_test), "ERROR: Test file "+ifile_test+" does not exist")
        expect(os.path.isfile(ifile_cntl), "ERROR: CNTL file "+ifile_cntl+" does not exist")
        
        ftest = Dataset(ifile_test,  mode='r')
        fcntl = Dataset(ifile_cntl, mode='r')

        #if max RMSE/DIFF is to be printed, extract lat lons from a file
        if(index):
            lat = ftest.variables['lat']
            lon = ftest.variables['lon']
        
        rmse = np.zeros(shape=(len(var_list)))
        icntvar = 0
        
        is_se = (len(ftest.variables[var_suffix+var_list[icntvar]].dimensions)==3) # see if it is SE grid
        
        for ivar in var_list:
            var = var_suffix+ivar
            if var in ftest.variables:
                vtest    = ftest.variables[var.strip()][0,...] # first dimention is time (=0)
                vcntl    = fcntl.variables[var.strip()][0,...] # first dimention is time (=0)

                #reshape for RMSE
                if(is_se ):
                    nx, ny = vtest.shape #shape will be same for both arrays
                    nz = 1
                else:
                    nx, ny, nz = vtest.shape #shape will be same for both arrays

                rmse[icntvar]  = sqrt(mean_squared_error(vtest.reshape((nx,ny*nz)), vcntl.reshape((nx,ny*nz))))

                if(index):
                    diff_arr = abs(vtest[...] - vcntl[...])
                    max_diff = np.amax(diff_arr)
                    ind_max = np.unravel_index(diff_arr.argmax(),diff_arr.shape)
                    print(ind_max)
                    print_str = var+' '+str(max_diff)+' '+str(vtest[ind_max])+' ' +str(vcntl[ind_max])+' '+str(lat[ind_max[1]])+' '+ str(lon[ind_max[1]])+' '+str(ind_max[0])
                    print("{}".format(print_str))


                #normalize by mean values of the field in the control case
                rmse[icntvar] = rmse[icntvar]/np.mean(vcntl)
                icntvar += 1
                vtest = None
                vcntl = None
                var   = None
        ftest.close()
        fcntl.close()
        return rmse

    def _compare_baseline(self):
        
        print('SLR_INFO:BASELINE COMPARISON STARTS')


        rundir    = self._case.get_value("RUNDIR")
        casename  = self._case.get_value("CASE") 

        var_list = self.get_var_list()
        len_var_list = len(var_list)

        #baseline directory names
        base_root = self._case.get_value("BASELINE_ROOT")
        base_comp = self._case.get_value("BASECMP_CASE")        

        #baseline directory is:base_root/base_comp
        base_dir = os.path.join(base_root,base_comp)

        #for test cases (new environment etc.)
        print('SLR_INFO: Test case comparison...')
        tst_res = np.empty([ninit_cond,len(prt),len_var_list])
        iinst = 0
        for icond in range(ninit_cond):
            iprt = 0
            for aprt in prtstr:
                iinst += 1;
                ifile_cntl = os.path.join(base_dir,self.get_fname_wo_ext('','',iinst)+'_'+aprt+'.nc')
                expect(os.path.isfile(ifile_cntl), "ERROR: File "+ifile_cntl+" does not exist")
                print('SLR_INFO:CNTL_TST:'+ifile_cntl)
                ifile_test = os.path.join(rundir,self.get_fname_wo_ext(rundir,casename,iinst)+'_'+aprt+'.nc')
                expect(os.path.isfile(ifile_test), "ERROR: File "+ifile_test+" does not exist")
                tst_res[icond,iprt,0:len_var_list] = self.rmse_var(ifile_test, ifile_cntl,  var_list, 't_' )
                print('SLR_INFO:Compared to TEST_TST:'+ifile_test)
                print(tst_res[icond,iprt,0:len_var_list])
                iprt += 1

        print('SLR_INFO: POST PROCESSING PHASE ENDS')
        

    #=================================================================
    # run_phase
    #=================================================================
    def run_phase(self):
        print('SLR_INFO: RUN PHASE')
        
        self.run_indv()

        #cwd = os.getcwd()
        #Here were are in case directory, we need to go to the run directory and rename files
        rundir   = self._case.get_value("RUNDIR")
        casename = self._case.get_value("CASE") 
        print('SLR_INFO: Case name is:'+casename )

        iinst = 1
        for icond in range(ninit_cond):
            for iprt in range(len(prt)):
                
                #------------------------------------------------------
                #SANITY CHECK - To confirm that history file extension
                #~~~~~~~~~~~~
                # corresponds to the right perturbation
                #------------------------------------------------------
                #find corresponding atm_in_*
                fatm_in = os.path.join(rundir,'atm_in_'+str(iinst).zfill(4))
                #see if atm_in_* file has pertlim same as prt[iprt] (sanity check)
                found  = False
                prtval = 0.0
                with open(fatm_in) as atmfile:
                    for line in atmfile:
                        if(line.find('pertlim') > 0):
                            found = True
                            prtval = float(line.split('=')[1])
                            expect(prtval == prt[iprt], "ERROR: prtval doesn't match, prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(prt[iprt]))
                            print("SLR_INFO:prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(prt[iprt]))
                
                if(not found):
                    expect(prtval == prt[iprt], "ERROR: default prtval doesn't match, prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(prt[iprt]))
                    print("SLR_INFO:def prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(prt[iprt]))
                  
                #---------------------------------------------------------
                #Rename file
                #---------------------------------------------------------

                #form file name
                fname = self.get_fname_wo_ext(rundir,casename,iinst) #NOTE:without extension
                print('SLR_INFO: fname to rename:'+fname )
                
                renamed_fname = fname +'_'+prtstr[iprt]+'.nc'
                fname_w_ext   = fname + '.nc' #add extention

                #see if file exists
                if (os.path.isfile(fname_w_ext)):
                    #rename file
                    shutil.move(fname_w_ext,renamed_fname)
                    print('SLR_INFO: Renamed file:'+renamed_fname)
                else:
                    expect(os.path.isfile(renamed_fname), "ERROR: File "+renamed_fname+" does not exist")
                    print('SLR_INFO: Renamed file already exists:'+renamed_fname)
                
                iinst += 1
        self._generate_baseline()

        if(self._case.get_value("GENERATE_BASELINE")):
            print('SLR_INFO: cloud generation-gen base' )


            var_list = self.get_var_list()
            len_var_list = len(var_list)

            #baseline directory names
            base_root = self._case.get_value("BASELINE_ROOT")
            base_gen  = self._case.get_value("BASEGEN_CASE")

            #baseline directory is:base_root/base_gen
            base_dir = os.path.join(base_root,base_gen)

            #for trusted cloud sims
            print('SLR_INFO: Computing cloud')
            cld_res = np.empty([ninit_cond,len(prt)-1,len_var_list])
            iinst = 0
            for icond in range(ninit_cond):

                iinst += 1;  iprt = 0
                ifile_cntl = os.path.join(base_dir,self.get_fname_wo_ext('','',iinst)+'_'+prtstr[iprt]+'.nc')
                expect(os.path.isfile(ifile_cntl), "ERROR: File "+ifile_cntl+" does not exist")
                print('SLR_INFO:CNTL_CLD:'+ifile_cntl)            
            
                for aprt in prtstr[1:]:
                    iinst += 1;
                    ifile_test = os.path.join(base_dir,self.get_fname_wo_ext('','',iinst)+'_'+aprt+'.nc')
                    expect(os.path.isfile(ifile_test), "ERROR: File "+ifile_test+" does not exist")
                    cld_res[icond,iprt,0:len_var_list] = self.rmse_var(ifile_test, ifile_cntl,  var_list, 't_' )
                    print('SLR_INFO:Compared to CLD:'+ifile_test)
                    print(cld_res[icond,iprt,0:len_var_list])

        #store results of cld_res in a text file and use it when comparing baselines...
        print('SLR_INFO: RUN PHASE ENDS' )


#=====================================================
#Debugging:
#=====================================================

#-----------------------------------------------------
#DEBUG type 1 (DB1): Ensure that model produces BFB instances
#-----------------------------------------------------
##Replace the following lines in the script with the following updated values
#prt        = [0.0, 0.0] #DB1
#prtstr     = ['woprt','posprt']

## Comment out all namelist changes so that all instances have same namelist values
## For testing effect of namelist changes, uncomment namlist changes such that
## namelists are same for all instances

## Comment out sanity checks as well if needed....
