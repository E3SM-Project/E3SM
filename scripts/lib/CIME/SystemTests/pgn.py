"""
Perturbation Growth New (PGN) - The CESM/ACME model's
multi-instance capability is used to conduct an ensemble
of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.

"""

import os
import glob
import shutil
import numpy as np
import scipy.stats as stats


from netCDF4 import Dataset
from math import sqrt
from sklearn.metrics import mean_squared_error

from CIME.test_status import *
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.build import post_build
from CIME.hist_utils import rename_all_hist_files
from CIME.utils import expect

#Logic for PGN ensemble runs:
# We need two inputs:
# A. Number of inic cond files
# B. perturbations (e.g. currently we have: without prt, pos prt and neg prt)

# Based off the above, we compute number of instances to have (A * B)
# Build phase               : Change the user_nl* files to add perturbations and other flags 
# Run pahse                 : Rename history files
# Baselines generation phase: Compute cloud, store in netcdf file, copy netcdf file to baseline folder
# Baseline Comparison phase : Compare against baseline cloud to know pass/fail, and plot (optional)


#TODO: 
#1. The loop to change user_nl* files is run twice as we call build again after changing ninst, which changes ntasks
#2. Do we want to remove user_nl_cam, user_nl_clm etc. files as they are not used in the simulation?
#3. change code so that pergro_ptend_names.txt is not generated if it is already there or only one instance writes this file...
#4. Plot generation is very basic at this point(no lables, legends etc.), improve it! 
#5. Decision making about PASS/FAIL should have multiple criteria

logger = logging.getLogger(__name__)

#--------------------------------------------------------
#Variables which needs global scope for various functions
#--------------------------------------------------------

#number of inital conditions
NINIT_COND = 6 #12
        
#perturbations for runs
PRT        = [0.0, 1.0e-14, -1.0e-14]
PRTSTR     = ['woprt','posprt','negprt']

#file name for file containing PGE cloud
FCLD_NC  = 'cam.h0.cloud.nc'

#For preparing paths for namelist files for initial condition files
FILE_PREF_ATM = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.cam.i.0002-"
FILE_PREF_LND = "SMS_Ly5.ne4_ne4.FC5AV1C-04P2.eos_intel.ne45y.clm2.r.0002-"

FILE_SUF_ATM = "-01-00000.nc"
FILE_SUF_LND = "-01-00000.nc"

#------------------------------------------------------------
# Some flags for debugging or invoking extra features
#------------------------------------------------------------

#DEBUGGING: prints out max rmse diffs if set to True
INDEX = False 

# generate a plot
DOPLOT = False  #It impacts performance!

class PGN(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the PGN test
        """

        #perturbation values and number of perturbation strings should be same
        expect(len(PRT)== len(PRTSTR),"Number of perturbation values ("+str(len(PRT))+") are NOT " 
               "equal to number of perturbation strings("+ str(len(PRTSTR))+")")
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
        nprt       = len(PRT)
        ninst      = NINIT_COND * nprt

        logger.debug('PGN_INFO: number of instance: '+str(ninst))
        
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
        #we launch ./case.build in the case directory of PGN 
        #test
        #------------------------------------------------------

        if(default_ninst == 1): #if multi-instance is not already set
            # Only want this to happen once. It will impact the sharedlib build
            # so it has to happen here.
            if not model_only:
                # Lay all of the components out concurrently
                logger.debug("PGN_INFO: Updating NINST for multi-instance in env_mach_pes.xml")
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
            logger.debug("PGN_INFO: FAKE Build")
            if (not sharedlib_only):
                post_build(self._case, [])
        else:
            # Build exectuable with multiple instances
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)


        #----------------------------------------------------------------
        # Namelist settings:
        #~~~~~~~~~~~~~~~~~~
        # Do this already in build_phase so that we can check the xml and
        # namelist files before job starts running.
        #----------------------------------------------------------------
        logger.debug("PGN_INFO: Updating user_nl_* files")

        csmdata_root = self._case.get_value("DIN_LOC_ROOT")
        csmdata_atm  = csmdata_root+"/atm/cam/inic/homme/ne4_v1_init/"
        csmdata_lnd  = csmdata_root+"/lnd/clm2/initdata/ne4_v1_init/b58d55680/"

        iinst = 1
        for icond in range(NINIT_COND):
            icond_label_2digits = str(icond+1).zfill(2)
            fatm_in = FILE_PREF_ATM + icond_label_2digits + FILE_SUF_ATM
            flnd_in = FILE_PREF_LND + icond_label_2digits + FILE_SUF_LND
            for iprt in PRT:
                with open('user_nl_cam_'+str(iinst).zfill(4), 'w') as atmnlfile, \
                        open('user_nl_clm_'+str(iinst).zfill(4), 'w') as lndnlfile:
                    
                    #atm/lnd intitial conditions                   
                    
                    #initial condition files to use for atm and land
                    #atmnlfile.write("ncdata  = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + fatm_in+"' \n")
                    #lndnlfile.write("finidat = '"+ "/pic/projects/uq_climate/wanh895/acme_input/ne4_v1_init/" + flnd_in+"' \n")
                    
                    #uncomment the following when there files are on SVN server
                    atmnlfile.write("ncdata  = '"+ csmdata_atm + "/" + fatm_in+"' \n")
                    lndnlfile.write("finidat = '"+ csmdata_lnd + "/" + flnd_in+"' \n")
                    
                    
                    #atm model output
                    atmnlfile.write("avgflag_pertape = 'I' \n")
                    atmnlfile.write("nhtfrq = 1 \n")
                    atmnlfile.write("mfilt  = 2  \n")
                    atmnlfile.write("ndens  = 1  \n")
                    atmnlfile.write("pergro_mods  = .true. \n")
                    atmnlfile.write("pergro_test_active = .true. \n")

                    #atmnlfile.write("empty_htapes = .true. \n")
                    #atmnlfile.write("fincl1 = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC' \n")
                    #atmnlfile.write("phys_debug_lat = 41.3495891345")
                    #atmnlfile.write("phys_debug_lon = 45.0" )
                
                    if(iprt != 0.0):
                        atmnlfile.write("pertlim = "+str(iprt)+" \n")
                        
                    iinst += 1

        #--------------------------------
        #Settings common to all instances
        #--------------------------------

        #Coupler settings which are applicable to ALL the instances (as there is only one coupler for all instances)
        self._case.set_value("STOP_N",     "1")
        self._case.set_value("STOP_OPTION","nsteps")


    #===========================================================
    # Some user defined functions to avoid repeated calculations
    #===========================================================

    def get_fname_wo_ext(self,rundir,casename,iinst):
        """
        construct file name given the input
        """
        #form file name withOUT extension
        tstmp    = '00000'
        date_str = '.h0.0001-01-01-'
        cam_str  = 'cam_'
        if(casename != "" ):
            cam_str = '.'+cam_str
        return os.path.join(rundir,casename+cam_str+str(iinst).zfill(4)+date_str+tstmp)

    def get_var_list(self):
        """
        Get variable list for pergro specific output vars
        """
        rundir    = self._case.get_value("RUNDIR")
        casename  = self._case.get_value("CASE")    
        prg_fname = 'pergro_ptend_names.txt'
        var_file = os.path.join(rundir,prg_fname)
        expect(os.path.isfile(var_file),"File "+prg_fname+" does not exist in: "+rundir)

        with open(var_file, 'r') as fvar:
            var_list = fvar.readlines()
                    
        return  map(str.strip,var_list)

    def nc_write_handle(self,fname_nc,rmse_nc_var):
        """
        Opens and write netcdf file for PGE curves
        This function is here purely to avoid duplicate 
        codes so that it is easy to maintain code longterm        
        """

        fhndl = Dataset(fname_nc,'w', format='NETCDF4')

        var_list     = self.get_var_list()
        len_var_list = len(var_list)
        nprt         = len(PRT)

        #create dims
        fhndl.createDimension('ninit', NINIT_COND)
        fhndl.createDimension('nprt', nprt)
        fhndl.createDimension('nprt_m1', nprt-1)
        fhndl.createDimension('nvars', len_var_list)


        #create variables in the file
        init_cond_nc = fhndl.createVariable('init_cond_fnames', 'S100', 'ninit')
        prt_nc       = fhndl.createVariable('perturb_strings', 'S10', 'nprt')
        variables_nc = fhndl.createVariable('perturb_vnames', 'S20', 'nvars')
        rmse_nc      = fhndl.createVariable(rmse_nc_var, 'f8', ('ninit', 'nprt_m1', 'nvars'))


        #assign variables for writing to the netcdf file
        for iprt in range(nprt):
            prt_nc[iprt]       = PRTSTR[iprt]

        for ivar in range(len_var_list):
            variables_nc[ivar] = var_list[ivar]

        for icond in range(NINIT_COND):
            icond_label_2digits = str(icond+1).zfill(2)
            init_cond_nc[icond] =  FILE_PREF_ATM + icond_label_2digits + FILE_SUF_ATM

        return fhndl, rmse_nc
        
    def rmse_var(self,ifile_test, ifile_cntl,  var_list, var_suffix ):
        """
        Compute RMSE difference between ifile_test and ifile_cntl for
        variables listed in var_list
        """
        
        #------------------ARGS-------------------------
        #ifile_test: path of test file
        #ifile_cntl: path of control file
        #var_list  : List of all variables
        #var_suffix: Suffix for var_list (e.g. t_, t_ qv_ etc.)
        #-----------------------------------------------

        # See if the files exists or not....
        expect(os.path.isfile(ifile_test), "ERROR: Test file "+ifile_test+" does not exist")
        expect(os.path.isfile(ifile_cntl), "ERROR: CNTL file "+ifile_cntl+" does not exist")
        
        ftest = Dataset(ifile_test,  mode='r')
        fcntl = Dataset(ifile_cntl, mode='r')

        #if max RMSE/DIFF is to be printed, extract lat lons from a file
        if(INDEX):
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

                if(INDEX):
                    #NEED to work on formatting this....
                    diff_arr = abs(vtest[...] - vcntl[...])
                    max_diff = np.amax(diff_arr)
                    ind_max = np.unravel_index(diff_arr.argmax(),diff_arr.shape)
                    print(ind_max)
                    print_str = var+' '+str(max_diff)+' '+str(vtest[ind_max])+' ' +str(vcntl[ind_max])+' '+str(lat[ind_max[1]])+' '+ str(lon[ind_max[1]])+' '+str(ind_max[0])
                    print("{}".format(print_str))
                    #print('{0:15}{1:15}{2:15}\n'.format(name, a, b))#TRY THIS!!


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
        
        #import time
        #t0 = time.time()
        """
        Compare baselines in the pergro test sense. That is,
        compare PGE from the test simulation with the baseline 
        cloud
        """
        logger.debug("PGN_INFO:BASELINE COMPARISON STARTS")

        rundir    = self._case.get_value("RUNDIR")
        casename  = self._case.get_value("CASE") 

        var_list = self.get_var_list()
        len_var_list = len(var_list)
        nprt         = len(PRT)

        #baseline directory names
        base_root = self._case.get_value("BASELINE_ROOT")
        base_comp = self._case.get_value("BASECMP_CASE")        

        #baseline directory is:base_root/base_comp
        base_dir = os.path.join(base_root,base_comp)

        #for test cases (new environment etc.)        
        logger.debug("PGN_INFO: Test case comparison...")

        #---------------------------------------------
        # Write netcdf file for comparison
        #---------------------------------------------
        fcomp_nc = 'comp_cld.nc'
        fcomp_cld, comp_rmse_nc  = self.nc_write_handle(fcomp_nc,'comp_rmse')

        iinst = 0
        for icond in range(NINIT_COND):
            iinst += 1;  iprt = 0
            ifile_cntl = os.path.join(base_dir,self.get_fname_wo_ext('','',iinst)+'_woprt.nc')
            expect(os.path.isfile(ifile_cntl), "ERROR: File "+ifile_cntl+" does not exist")
            logger.debug("PGN_INFO:CNTL_TST:"+ifile_cntl)

            for aprt in PRTSTR[1:]:
                iinst += 1;                
                ifile_test = os.path.join(rundir,self.get_fname_wo_ext(rundir,casename,iinst)+'_'+aprt+'.nc')
                expect(os.path.isfile(ifile_test), "ERROR: File "+ifile_test+" does not exist")
                comp_rmse_nc[icond,iprt,0:len_var_list] = self.rmse_var(ifile_test, ifile_cntl,  var_list, 't_' )
                logger.debug("PGN_INFO:Compared to TEST_TST:"+ifile_test)

                iprt += 1
        
        #--------------------------------------------
        #Student's t-test
        #-------------------------------------------

        #Extract last element of each PGE curve of the cloud
        path_cld_nc  = os.path.join(base_dir,FCLD_NC)
        expect(os.path.isfile(path_cld_nc), "ERROR: "+FCLD_NC+"  does not exist at:"+path_cld_nc)

        fcld         = Dataset(path_cld_nc,'r', format='NETCDF4')

        #Get dimentions of cloud PGE curves
        ninic_cld = fcld.variables['cld_rmse'].shape[0]
        nprt_cld  = fcld.variables['cld_rmse'].shape[1]
        nvar_cld  = fcld.variables['cld_rmse'].shape[2]

        expect(ninic_cld == NINIT_COND, "BASELINE COMPARISON: Number of initial conditions should be same:" \
               "inic cld:"+str(ninic_cld)+"  inic comp:"+str(NINIT_COND))

        expect(nprt_cld == nprt-1, "BASELINE COMPARISON: Number of perturbations should be same:" \
               "prt cld:"+str(nprt_cld)+"  prt comp:"+str(nprt - 1))

        expect(nvar_cld == len_var_list, "BASELINE COMPARISON: Number of phys update variables should be same:" \
               "nvar cld:"+str(nvar_cld)+"  nvar comp:"+str(len_var_list))


        pgecld =  fcld.variables['cld_rmse']

        if (DOPLOT):
            import pylab as pl            
            #generate plot
            ax = pl.gca()
            for icond in range(NINIT_COND):
                for iprt in range(nprt-1):
                    ax.semilogy(pgecld[icond,iprt,:],color='b')
                    ax.semilogy(comp_rmse_nc[icond,iprt,:],color='r')
            pl.savefig('plot_comp.png')

        pge_ends_cld  = fcld.variables['cld_rmse'][:,:,len_var_list-1]
        pge_ends_comp = comp_rmse_nc[:,0:nprt-1,len_var_list-1]

        #run the t-test
        pge_ends_cld = pge_ends_cld.flatten()
        pge_ends_comp = pge_ends_comp.flatten()        

        t_stat, p_val = stats.ttest_ind(pge_ends_cld, pge_ends_comp)

        t_stat = abs(t_stat) #only value of t_stat matters, not the sign

        print("PGN_INFO: T value:"+str(t_stat))
        print("PGN_INFO: P value:"+str(p_val))

        logger.warn(" T value:"+str(t_stat))
        logger.warn(" P value:"+str(p_val))
        
        
        #There should be multiple criteria to determin pass/fail
        #This part of decision making should be more polished

        if (t_stat > 3):
            with self._test_status:
                self._test_status.set_status(BASELINE_PHASE, TEST_FAIL_STATUS)
            print('PGN_INFO:TEST FAIL')
        else:
            with self._test_status:
                self._test_status.set_status(BASELINE_PHASE, TEST_PASS_STATUS)
            print('PGN_INFO:TEST PASS')
                    
        #Close this file after you access "comp_rmse_nc" variable
        fcomp_cld.close()
        logger.debug("PGN_INFO: POST PROCESSING PHASE ENDS")
        
        #t1 = time.time()
        #print('time: '+str(t1-t0))
    #=================================================================
    # run_phase
    #=================================================================
    def run_phase(self):
        logger.debug("PGN_INFO: RUN PHASE")
        
        self.run_indv()

        #cwd = os.getcwd()

        #Here were are in case directory, we need to go to the run directory and rename files
        rundir   = self._case.get_value("RUNDIR")
        casename = self._case.get_value("CASE") 
        logger.debug("PGN_INFO: Case name is:"+casename )

        iinst = 1
        for icond in range(NINIT_COND):
            for iprt in range(len(PRT)):
                
                #------------------------------------------------------
                #SANITY CHECK - To confirm that history file extension
                #~~~~~~~~~~~~
                # corresponds to the right perturbation
                #------------------------------------------------------
                #find corresponding atm_in_*
                fatm_in = os.path.join(rundir,'atm_in_'+str(iinst).zfill(4))
                #see if atm_in_* file has pertlim same as PRT[iprt] (sanity check)
                found  = False
                prtval = 0.0
                with open(fatm_in) as atmfile:
                    for line in atmfile:
                        if(line.find('pertlim') > 0):
                            found = True
                            prtval = float(line.split('=')[1])
                            expect(prtval == PRT[iprt], "ERROR: prtval doesn't match, prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(PRT[iprt]))
                            logger.debug("PGN_INFO:prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(PRT[iprt]))
                
                if(not found):
                    expect(prtval == PRT[iprt], "ERROR: default prtval doesn't match, prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(PRT[iprt]))
                    logger.debug("PGN_INFO:def prtval:"+str(prtval)+"; prt["+str(iprt)+"]:"+str(PRT[iprt]))
                  
                #---------------------------------------------------------
                #Rename file
                #---------------------------------------------------------

                #form file name
                fname = self.get_fname_wo_ext(rundir,casename,iinst) #NOTE:without extension
                logger.debug("PGN_INFO: fname to rename:"+fname )
                
                renamed_fname = fname +'_'+PRTSTR[iprt]+'.nc'
                fname_w_ext   = fname + '.nc' #add extention

                #see if file exists
                if (os.path.isfile(fname_w_ext)):
                    #rename file
                    shutil.move(fname_w_ext,renamed_fname)
                    logger.debug("PGN_INFO: Renamed file:"+renamed_fname)
                else:
                    expect(os.path.isfile(renamed_fname), "ERROR: File "+renamed_fname+" does not exist")
                    logger.debug("PGN_INFO: Renamed file already exists:"+renamed_fname)
                
                iinst += 1

        #cloud generation

        logger.debug("PGN_INFO: cloud generation-gen base")
            
        var_list     = self.get_var_list()
        len_var_list = len(var_list)
        nprt         = len(PRT)
        
        #for trusted cloud sims
        logger.debug("PGN_INFO: Computing cloud")
            
        cld_res = np.empty([len_var_list])
            
        #---------------------------------------------
        # Write netcdf file for cloud in rundir
        #---------------------------------------------
        os.chdir(rundir)
        fcld, cld_rmse_nc = self.nc_write_handle(FCLD_NC,'cld_rmse')

        iinst = 0
        for icond in range(NINIT_COND):
            iinst += 1;  iprt = 0
            ifile_cntl = os.path.join(rundir,self.get_fname_wo_ext('',casename,iinst)+'_'+PRTSTR[iprt]+'.nc')
            expect(os.path.isfile(ifile_cntl), "ERROR: File "+ifile_cntl+" does not exist")
            logger.debug("PGN_INFO:CNTL_CLD:"+ifile_cntl)            
            
            for aprt in PRTSTR[1:]:
                iinst += 1;
                iprt = PRTSTR.index(aprt) - 1 # subtracted 1 as we will only get two curves for each inic (wo-pos and wo-neg)
                ifile_test = os.path.join(rundir,self.get_fname_wo_ext('',casename,iinst)+'_'+aprt+'.nc')
                expect(os.path.isfile(ifile_test), "ERROR: File "+ifile_test+" does not exist")
                cld_rmse_nc[icond,iprt,0:len_var_list] = self.rmse_var(ifile_test, ifile_cntl,  var_list, 't_' )
                logger.debug("PGN_INFO:Compared to CLD:"+ifile_test)

        fcld.close()

        if(self._case.get_value("GENERATE_BASELINE")):
            
            #baseline directory names
            base_root = self._case.get_value("BASELINE_ROOT")
            base_gen  = self._case.get_value("BASEGEN_CASE")
        
            #baseline directory is:base_root/base_gen
            base_dir = os.path.join(base_root,base_gen)

            #first copy files to the baseline directory
            self._generate_baseline()#BALLI-CHECK IF THIS IS OKAY

            #copy cloud.nc file to baseline directory
            logger.debug("PGN_INFO:copy:"+FCLD_NC+" to "+base_dir)
            shutil.copy(FCLD_NC,base_dir)

        logger.debug("PGN_INFO: RUN PHASE ENDS" )


#=====================================================
#Debugging:
#=====================================================

#-----------------------------------------------------
#DEBUG type 1 (DB1): Ensure that model produces BFB instances
#-----------------------------------------------------
##Replace the following lines in the script with the following updated values
#PRT        = [0.0, 0.0] #DB1
#PRTSTR     = ['woprt','posprt']

## Comment out all namelist changes so that all instances have same namelist values
## For testing effect of namelist changes, uncomment namlist changes such that
## namelists are same for all instances

## Comment out sanity checks as well if needed....
