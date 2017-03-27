#!/usr/bin/env python

import getpass, os, sys, csv, math
from optparse import OptionParser
import subprocess
import numpy

parser = OptionParser();

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--runroot", dest="runroot", default="", \
                  help="Directory where the run would be created")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup",default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--machine", dest="machine", default = 'oic2', \
                  help = "machine to use (default = oic2)")
parser.add_option("--compiler", dest="compiler", default = 'gnu', \
                  help = "compiler to use (pgi*, gnu)")
parser.add_option("--mpilib", dest="mpilib", default="mpi-serial", \
                      help = "mpi library (openmpi*, mpich, ibm, mpi-serial)")
parser.add_option("--csmdir", dest="csmdir", default='../../../../..', \
                  help = "base CESM directory (default = ../../..)")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='', \
                  help = "input data directory for CESM (required)")
parser.add_option("--nopointdata", action="store_true", \
                  dest="nopointdata", help="Do NOT make point data (use data already created)", \
                  default=False)
parser.add_option("--noad", action="store_true", dest="noad", default=False, \
                  help='Do not perform ad spinup simulation')
parser.add_option("--notrans", action="store_true", dest="notrans", default=False, \
	          help='Do not perform transient simulation (spinup only)')
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--nyears_final_spinup", dest="nyears_final_spinup", default='200', \
                  help="base no. of years for final spinup")
parser.add_option("--clean_build", action="store_true", default=False, \
                  help="Perform a clean build")
parser.add_option("--hist_mfilt_trans", dest="hist_mfilt", default="365", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_nhtfrq_trans", dest="hist_nhtfrq", default="-24", \
                  help = 'output file timestep (transient only)')
parser.add_option("--hist_mfilt_spinup", dest="hist_mfilt_spinup", default="-999", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_nhtfrq_spinup", dest="hist_nhtfrq_spinup", default="-999", \
                  help = 'output file timestep (transient only)')
parser.add_option("--parm_file", dest="parm_file", default="", \
                  help = 'parameter file to use')
parser.add_option("--nopftdyn", action="store_true", dest="nopftdyn", \
                      default=False, help='Do not use dynamic PFT file')
parser.add_option("--nofire", dest="nofire", default=False, \
                  action="store_true", help='Turn off fire algorithms')
parser.add_option("--regional", action="store_true", \
                   dest="regional", default=False, \
                   help="Flag for regional run (2x2 or greater)")
parser.add_option("--np", dest="np", default=1, \
                  help = 'number of processors')
parser.add_option("--tstep", dest="tstep", default=0.5, \
                  help = 'CLM timestep (hours)')
parser.add_option("--co2_file", dest="co2_file", default="fco2_datm_1765-2007_c100614.nc", \
                  help = 'CLM timestep (hours)')
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=250, \
                  help = 'number of years to run ad_spinup')
parser.add_option("--nyears_transient", dest="nyears_transient", default=-1, \
                  help = 'number of years to run transient')
parser.add_option("--metdir", dest="metdir", default="none", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--C13", dest="C13", default=False, \
                      help = 'Switch to turn on C13', action="store_true")
parser.add_option("--C14", dest="C14", default=False, \
                  help = 'Use C14 as C13 (no decay)', action="store_true")
parser.add_option("--ninst", dest="ninst", default=1, \
                      help = 'number of land model instances')
parser.add_option("--BGC", action="store_true", dest="bgc", default=False, \
                  help = 'Use BGC compset')
parser.add_option("--harvmod", action="store_true", dest='harvmod', default=False, \
                    help="turn on harvest modification:  All harvest at first timestep")
parser.add_option("--bulk_denitrif", dest="bulk_denitrif", default=False, \
                  help = 'To turn on BGC nitrification-denitrification', action="store_true")
parser.add_option("--no_dynroot", dest="no_dynroot", default=False, \
                  help = 'Turn off dynamic root distribution', action="store_true")
parser.add_option("--vertsoilc", dest="vsoilc", default=False, \
                  help = 'To turn on CN with multiple soil layers, excluding CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me', action="store_true")
parser.add_option("--makemetdata", action="store_true", dest="makemet", default=False, \
                    help="generate site meteorology")
parser.add_option("--xpts", dest="xpts", default=1, \
                      help = 'for regional runs: xpts')
parser.add_option("--ypts", dest="ypts", default=1, \
                      help = 'for regional runs: ypts')
parser.add_option("--cruncep", dest="cruncep", default=False, \
                  action="store_true", help = 'Use CRU-NCEP meteorology')
parser.add_option("--surfdata_grid", dest="surfdata_grid", default=False, \
                  help = 'Use gridded surface data instead of site data', action="store_true")
parser.add_option("--siteparms",dest = "siteparms", default=False, \
                  action="store_true", help = 'Use default PFT parameters')
parser.add_option("--cpl_bypass", dest = "cpl_bypass", default=False, \
                  help = "Bypass coupler", action="store_true")
parser.add_option("--spinup_vars", dest = "spinup_vars", default=False, help = "limit output variables for spinup", action="store_true")
parser.add_option("--cn_only", dest="cn_only", default=False, help='Carbon/Nitrogen only (saturated P)', action ="store_true")
parser.add_option("--ensemble_file", dest="ensemble_file", default='', \
                  help = 'Parameter sample file to generate ensemble')
parser.add_option("--parm_list", dest="parm_list", default='parm_list', \
                  help = 'File containing list of parameters to vary')
parser.add_option("--mc_ensemble", dest="mc_ensemble", default=-1, \
                  help = 'Monte Carlo ensemble (argument is # of simulations)')
parser.add_option("--ng", dest="ng", default=256, \
                  help = 'number of groups to run in ensemble mode')
#Changed by Ming for mesabi
parser.add_option("--archiveroot", dest="archiveroot", default='', \
                  help = "archive root directory only for mesabi")
#Added by Kirk to include the modified parameter file
parser.add_option("--mod_parm_file", dest="mod_parm_file", default='', \
                  help = "adding the path to the modified parameter file")
parser.add_option("--walltime", dest="walltime", default=6, \
                  help = "desired walltime for each job (hours)")

(options, args) = parser.parse_args()

#------------ define function for pbs submission
def submit(fname, submit_type='qsub', job_depend=''):
    if (job_depend != ''):
        os.system(submit_type+' -W depend=afterok:'+job_depend+' '+fname+' > temp/jobinfo')
    else:
        os.system(submit_type+' '+fname+' > temp/jobinfo')
    myinput = open('temp/jobinfo')
    for s in myinput:
        thisjob = (s.split('.')[0]).strip('\n')
    myinput.close()
    os.system('rm temp/jobinfo')
    return thisjob

#----------------------------------------------------------
if (options.ccsm_input != ''):
    ccsm_input = options.ccsm_input
elif (options.machine == 'titan' or options.machine == 'eos'):
    ccsm_input = '/lustre/atlas/world-shared/cli900/cesm/inputdata'
elif (options.machine == 'cades'):
    ccsm_input = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/ACME_inputdata/'
elif (options.machine == 'edison' or options.machine == 'cori'):
    ccsm_input = '/project/projectdirs/acme/inputdata'

if (options.compiler != ''):
    if (options.machine == 'titan'):
        options.compiler = 'pgi'
    if (options.machine == 'eos' or options.machine == 'edison'):
        options.compiler = 'intel'
    if (options.machine == 'cades'):
        options.compiler = 'gnu'
    

mycaseid   = options.mycaseid
srcmods    = options.srcmods_loc

#get start and year of input meteorology from site data file
PTCLMfiledir = ccsm_input+'/lnd/clm2/PTCLM'
fname = PTCLMfiledir+'/'+options.sitegroup+'_sitedata.txt'
AFdatareader = csv.reader(open(fname, "rt"))

translen = int(options.nyears_transient)

if (options.bgc):
    mybgc='BGC'
else:
    mybgc='CN'

csmdir = os.getcwd()+'/'+options.csmdir

#case run and case root directories
if (options.runroot == '' or (os.path.exists(options.runroot) == False)):
    myuser = getpass.getuser()
    if (options.machine == 'titan' or options.machine == 'eos'):
        myinput = open('/ccs/home/'+myuser+'/.cesm_proj','r')
        for s in myinput:
	   myproject=s[:-1]
        runroot='/lustre/atlas/scratch/'+myuser+'/'+myproject
    elif (options.machine == 'cades'):
        runroot='/lustre/or-hydra/cades-ccsi/scratch/'+myuser
    else:
        runroot = csmdir+'/run'
else:
    runroot = os.path.abspath(options.runroot)
if (options.caseroot == '' or (os.path.exists(options.caseroot) == False)):
    caseroot = os.path.abspath(csmdir+'/cime/scripts')
else:
    caseroot = os.path.abspath(options.caseroot)


isfirstsite = True
#create ensemble file if requested (so that all cases use the same)
if (int(options.mc_ensemble) != -1):
    if (not(os.path.isfile(options.parm_list))):
	print('parm_list file does not exist')
        sys.exit()
    else:
        param_names=[]
        param_min=[]
        param_max=[]
        input = open(options.parm_list,'r')
        for s in input:
	    if (s):
                param_names.append(s.split()[0])
                if (int(options.mc_ensemble) > 0):
                    if (len(s.split()) == 3):
                        param_min.append(float(s.split()[1]))
                        param_max.append(float(s.split()[2]))
                    else:
                        param_min.append(float(s.split()[2]))
                        param_max.append(float(s.split()[3]))
        input.close() 
        n_parameters = len(param_names)
    nsamples = int(options.mc_ensemble)
    samples=numpy.zeros((n_parameters,nsamples), dtype=numpy.float)
    for i in range(0,nsamples):
        for j in range(0,n_parameters):
            samples[j][i] = param_min[j]+(param_max[j]-param_min[j])*numpy.random.rand(1)
    numpy.savetxt('mcsamples_'+options.mycaseid+'_'+str(options.mc_ensemble)+'.txt', \
                  numpy.transpose(samples))
    options.ensemble_file = 'mcsamples_'+options.mycaseid+'_'+str(options.mc_ensemble)+'.txt'

for row in AFdatareader:
    if row[0] == options.site or (options.site == 'all' and row[0] !='site_code' \
                                      and row[0] != ''):
        site      = row[0]
        if (options.cruncep):
                startyear = 1901
                endyear = 1920
        else:
            startyear = int(row[6])
            endyear   = int(row[7])

        site_endyear = int(row[7])
        ncycle   = endyear-startyear+1   #number of years in met cycle
        ny_ad = options.ny_ad
	ny_fin = options.nyears_final_spinup
        if (int(options.ny_ad) % ncycle != 0):
          #AD spinup and final spinup lengths must be multiples of met data cyle.
          ny_ad = str(int(ny_ad) + ncycle - (int(ny_ad) % ncycle))
        if (int(options.nyears_final_spinup) % ncycle !=0):
          ny_fin = str(int(ny_fin) + ncycle - (int(ny_fin) % ncycle))

        if (translen == -1):
          translen = endyear-1850+1        #length of transient run
	  if (options.cpl_bypass and options.cruncep):
 	    translen = min(site_endyear,2010)-1850+1

        #use site parameter file if it exists
        if (options.siteparms):
            if (os.path.exists(PTCLMfiledir+'/parms_'+site)):
                print ('Using parameter file PTCLM_Files/parms_'+site)
                options.parm_file = PTCLMfiledir+'/parms_'+site
            else:
                options.parm_file = ''

        fsplen = int(ny_fin)
 
        #get align_year
        year_align = (endyear-1850+1) % ncycle

        #print year_align, fsplen
        basecmd = 'python pointCLM.py --site '+site+' --ccsm_input '+ \
            os.path.abspath(ccsm_input)+' --rmold --no_submit --sitegroup ' + \
            options.sitegroup+' --xpts '+str(options.xpts)+' --ypts '+str(options.ypts)+ \
            ' --machine '+options.machine
        if (srcmods != ''):
            srcmods    = os.path.abspath(srcmods)
            basecmd = basecmd+' --srcmods_loc '+srcmods
        if (mycaseid != ''):
            basecmd = basecmd+' --caseidprefix '+mycaseid
        if (options.parm_file != ''):
            basecmd = basecmd+' --parm_file '+options.parm_file
        if (options.clean_build):
            basecmd = basecmd+' --clean_build '
        if (options.regional):
            basecmd = basecmd+' --regional --xpts 2 --ypts 1 '
        if (options.metdir !='none'):
            basecmd = basecmd+' --metdir '+options.metdir
        if (options.C13):
            basecmd = basecmd+' --C13 '
        if (options.C13):
            basecmd = basecmd+' --C14 '
        if (options.ninst > 1):
            basecmd = basecmd+' --ninst '+str(options.ninst)
        if (options.nofire):
            basecmd = basecmd+' --nofire'
        if (options.harvmod):
            basecmd = basecmd+' --harvmod'
        if (options.nopftdyn):
            basecmd = basecmd+' --nopftdyn'
        if (options.no_dynroot):
            basecmd = basecmd+' --no_dynroot'
        if (options.bulk_denitrif):
            basecmd = basecmd+' --bulk_denitrif'
        if (options.vsoilc):
            basecmd = basecmd+' --vertsoilc'
        if (options.centbgc):
            basecmd = basecmd+' --centbgc'
        if (options.cn_only):
            basecmd = basecmd+' --cn_only'
        if (options.CH4):
            basecmd = basecmd+' --CH4'
        if (options.cruncep):
            basecmd = basecmd+' --cruncep'
        if (options.nopointdata):
            basecmd = basecmd+' --nopointdata'
        if (options.surfdata_grid):
            basecmd = basecmd+' --surfdata_grid'
        if (options.ensemble_file != ''):   
            basecmd = basecmd+' --ensemble_file '+options.ensemble_file
            basecmd = basecmd+' --parm_list '+options.parm_list
        if (options.archiveroot !=''):
            basecmd = basecmd+' --archiveroot '+options.archiveroot
        if (options.mod_parm_file !=''):
            basecmd = basecmd+' --mod_parm_file '+options.mod_parm_file
        basecmd = basecmd + ' --ng '+str(options.ng)
        basecmd = basecmd + ' --np '+str(options.np)
        basecmd = basecmd + ' --tstep '+str(options.tstep)
        basecmd = basecmd + ' --co2_file '+options.co2_file
        basecmd = basecmd + ' --compiler '+options.compiler
        basecmd = basecmd + ' --mpilib '+options.mpilib
        basecmd = basecmd+' --caseroot '+caseroot
        basecmd = basecmd+' --runroot '+runroot
        basecmd = basecmd+' --walltime '+str(options.walltime)

#----------------------- build commands for runCLM.py -----------------------------

        #AD spinup
        cmd_adsp = basecmd+' --ad_spinup --nyears_ad_spinup '+ \
            str(ny_ad)+' --align_year '+str(year_align+1)
        if (int(options.hist_mfilt_spinup) == -999):
            cmd_adsp = cmd_adsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
            str((endyear-startyear+1)*8760)
        else:
            cmd_adsp = cmd_adsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
                   +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)
        if (not isfirstsite):
            cmd_adsp = cmd_adsp+' --exeroot '+ad_exeroot+' --no_build'
        if (options.cpl_bypass):
            cmd_adsp = cmd_adsp+' --compset I1850CLM45CB'+mybgc
            ad_case = site+'_I1850CLM45CB'+mybgc
        else:
            cmd_adsp = cmd_adsp+' --compset I1850CLM45'+mybgc
            ad_case = site+'_I1850CLM45'+mybgc
        if (options.noad == False):
	    ad_case = ad_case+'_ad_spinup'
        if (options.makemet):
            cmd_adsp = cmd_adsp+' --makemetdat'
        if (options.spinup_vars):
            cmd_adsp = cmd_adsp+' --spinup_vars'
        if (mycaseid != ''):
            ad_case = mycaseid+'_'+ad_case
        if (isfirstsite):
            ad_exeroot = os.path.abspath(runroot+'/'+ad_case+'/bld')

        #final spinup
        if mycaseid !='':
            basecase=mycaseid+'_'+site
            if (options.cpl_bypass):
                basecase = basecase+'_I1850CLM45CB'+mybgc
            else: 
                basecase = basecase+'_I1850CLM45'+mybgc
        else:
            if (options.cpl_bypass):
                basecase=site+'_I1850CLM45CB'+mybgc
            else:
                basecase=site+'_I1850CLM45'+mybgc
        if (options.noad):
            cmd_fnsp = basecmd+' --run_units nyears --run_n '+str(fsplen)+' --align_year '+ \
                       str(year_align+1)+' --coldstart'
            if (not isfirstsite):
                cmd_fnsp = cmd_fnsp+' --exeroot '+ad_exeroot+' --no_build'
        else:
            cmd_fnsp = basecmd+' --finidat_case '+basecase+'_ad_spinup '+ \
                       '--finidat_year '+str(int(ny_ad)+1)+' --run_units nyears --run_n '+ \
                       str(fsplen)+' --align_year '+str(year_align+1)+' --no_build' + \
                       ' --exeroot '+ad_exeroot
        if (int(options.hist_mfilt_spinup) == -999):
            cmd_fnsp = cmd_fnsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
            str((endyear-startyear+1)*8760)
        else:
            cmd_fnsp = cmd_fnsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
                   +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)
        if (options.cpl_bypass):
            cmd_fnsp = cmd_fnsp+' --compset I1850CLM45CB'+mybgc
        else:
            cmd_fnsp = cmd_fnsp+' --compset I1850CLM45'+mybgc
        if (options.spinup_vars):
		cmd_fnsp = cmd_fnsp+' --spinup_vars'
        if (options.ensemble_file != '' and options.notrans):	
                cmd_fnsp = cmd_fnsp + ' --postproc_file postproc_vars'

        #transient
        cmd_trns = basecmd+' --finidat_case '+basecase+ \
            ' --finidat_year '+str(fsplen+1)+' --run_units nyears' \
            +' --run_n '+str(translen)+' --align_year '+ \
            str(year_align+1850)+' --hist_nhtfrq '+ \
            options.hist_nhtfrq+' --hist_mfilt '+options.hist_mfilt+' --no_build' + \
            ' --exeroot '+ad_exeroot
        if (options.cpl_bypass):
            cmd_trns = cmd_trns+' --compset I20TRCLM45CB'+mybgc
        else:
            cmd_trns = cmd_trns+' --compset I20TRCLM45'+mybgc
        if (options.spinup_vars):
            cmd_trns = cmd_trns + ' --spinup_vars'
        if (options.ensemble_file != ''):  #Transient post-processing
            cmd_trns = cmd_trns + ' --postproc_file postproc_vars'
        #transient phase 2 (CRU-NCEP only, without coupler bypass)
        if (options.cruncep and not options.cpl_bypass):
            basecase=basecase.replace('1850','20TR')+'_phase1'
            thistranslen = site_endyear - 1921 + 1
            cmd_trns2 = basecmd+' --trans2 --finidat_case '+basecase+ \
                ' --finidat_year 1921 --run_units nyears --branch ' \
                +' --run_n '+str(thistranslen)+' --align_year 1921'+ \
                ' --hist_nhtfrq '+options.hist_nhtfrq+' --hist_mfilt '+ \
                options.hist_mfilt+' --no_build'+' --exeroot '+ad_exeroot + \
                ' --compset I20TRCLM45'+mybgc
            print(cmd_trns2)

#---------------------------------------------------------------------------------

        #set site environment variable
        os.environ['SITE']=site
        basecase = site
        if (mycaseid != ''):
                basecase = mycaseid+'_'+site
        os.system('mkdir -p temp')

        #build cases
        print('\nSetting up ad_spinup case\n')
        if (options.noad == False):
            print cmd_adsp
            os.system(cmd_adsp)
        print('\nSetting up final spinup case\n')
        os.system(cmd_fnsp)
        if (options.notrans == False):
            print('\nSetting up transient case\n')
            print cmd_trns
            os.system(cmd_trns)
            if (options.cruncep and not options.cpl_bypass):
                 print('\nSetting up transient case phase 2\n')
                 os.system(cmd_trns2)
        
        output = open('./temp/site_fullrun.pbs','w')

        mysubmit_type = 'qsub'
        if (options.machine == 'cori' or options.machine == 'edison'):
            mysubmit_type = 'sbatch'
        #Create a .PBS site fullrun script to launch the full job (all 3 cases)
        if (options.cpl_bypass):
          input = open(caseroot+'/'+basecase+"_I1850CLM45CB"+mybgc+'/case.run')
        else:
          input = open(caseroot+'/'+basecase+"_I1850CLM45"+mybgc+'/case.run')
        for s in input:
            if ("perl" in s or "python" in s):
                output.write("#!/bin/csh -f\n")
                if (mysubmit_type == 'qsub'):
                    output.write('#PBS -l walltime='+str(options.walltime)+':00:00\n')
                else:
                    output.write('#SBATCH --time='+str(options.walltime)+':00:00\n')
                    if ('edison' in options.machine):
                        output.write('#SBATCH --partition=regular\n')
	    elif ("#!" in s or "#PBS" in s or "#SBATCH" in s):
                output.write(s)
        input.close()
        output.write("\n")
        
        if (options.machine == 'eos'):
            output.write('source $MODULESHOME/init/csh\n')
            output.write('module load nco\n')
            output.write('module unload python\n')
            output.write('module load python/2.7.5\n')
            output.write('module unload PrgEnv-intel\n')
            output.write('module load PrgEnv-gnu\n')
            output.write('module load python_numpy\n')
            output.write('module load python_scipy\n')
            output.write('module load python_mpi4py/2.0.0\n')
            output.write('module unload PrgEnv-gnu\n')
            output.write('module load PrgEnv-intel\n')
        if (options.machine == 'titan'):
            output.write('source $MODULESHOME/init/csh\n')
            output.write('module load nco\n')
            output.write('module load python\n')
            output.write('module load python_numpy/1.9.2\n')
            output.write('module load python_scipy/0.15.1\n')
            output.write('module load python_mpi4py/2.0.0\n')

        if (options.cpl_bypass):
          modelst = 'CLM45CB'+mybgc
        else:
          modelst = 'CLM45'+mybgc

        basecase = site
        if (mycaseid != ''):
                basecase = mycaseid+'_'+site
        if (options.noad == False):
            output.write("cd "+caseroot+'/'+basecase+"_I1850"+modelst+"_ad_spinup/\n")
            output.write("./case.submit --no-batch\n")
            output.write("cd "+os.path.abspath(".")+'\n')
            if (options.bgc):
                output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+ \
                             '/'+ad_case+'/run/ --casename '+ ad_case+' --restart_year '+ \
                             str(int(ny_ad)+1)+' --BGC\n')
            else:
                output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+ \
                             '/'+ad_case+'/run/ --casename '+ad_case+' --restart_year '+ \
                             str(int(ny_ad)+1)+'\n')
        output.write("cd "+caseroot+'/'+basecase+"_I1850"+modelst+"\n")
        output.write('./case.submit --no-batch\n')	
        if (options.notrans == False):
            output.write("cd "+caseroot+'/'+basecase+"_I20TR"+modelst+"\n")
            output.write('./case.submit --no-batch\n')
        output.close()

        #if ensemble simulations requested, submit jobs created by pointclm.py in correct order
        if (options.ensemble_file != ''):
            cases=[]
            #build list of cases for fullrun
            if (options.noad == False):
                cases.append(basecase+'_I1850'+modelst+'_ad_spinup')
            cases.append(basecase+'_I1850'+modelst)
            if (options.notrans == False):
                cases.append(basecase+'_I20TR'+modelst)
            job_depend_run=''    
            for thiscase in cases:
                job_depend_run = submit('temp/ensemble_run_'+thiscase+'.pbs',job_depend= \
                                        job_depend_run, submit_type=mysubmit_type)
        else:  #submit single job
            job_fullrun = submit('temp/site_fullrun.pbs', submit_type=mysubmit_type)
        isfirstsite = False
