#!/usr/bin/env python

import getpass, os, sys, csv, math
from optparse import OptionParser
import subprocess
import numpy
import re

parser = OptionParser();

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='', \
                  help = "input data directory for CESM (required)")
parser.add_option("--clean_build", action="store_true", default=False, \
                  help="Perform a clean build")
parser.add_option("--csmdir", dest="csmdir", default='../../../../..', \
                  help = "base CESM directory (default = ../../..)")
parser.add_option("--compiler", dest="compiler", default = 'gnu', \
                  help = "compiler to use (pgi*, gnu)")
parser.add_option("--diags", dest="diags", default=False, \
                 action="store_true", help="Write special outputs for diagnostics")
parser.add_option("--debug", dest="debug", default=False, \
                 action="store_true", help='Use debug queue and options')
parser.add_option("--hist_mfilt_trans", dest="hist_mfilt", default="365", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_mfilt_spinup", dest="hist_mfilt_spinup", default="-999", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_nhtfrq_spinup", dest="hist_nhtfrq_spinup", default="-999", \
                  help = 'output file timestep (transient only)')
parser.add_option("--hist_nhtfrq_trans", dest="hist_nhtfrq", default="-24", \
                  help = 'output file timestep (transient only)')
parser.add_option("--humhol", dest="humhol", default=False, \
                  help = 'Use hummock/hollow microtopography', action="store_true")
parser.add_option("--marsh", dest="marsh", default=False, \
                  help = 'Use marsh hydrology/elevation', action="store_true")
parser.add_option("--machine", dest="machine", default = 'oic2', \
                  help = "machine to use (default = oic2)")
parser.add_option("--mpilib", dest="mpilib", default="mpi-serial", \
                      help = "mpi library (openmpi*, mpich, ibm, mpi-serial)")
parser.add_option("--noad", action="store_true", dest="noad", default=False, \
                  help='Do not perform ad spinup simulation')
parser.add_option("--nofire", dest="nofire", default=False, \
                  action="store_true", help='Turn off fire algorithms')
parser.add_option("--nopointdata", action="store_true", \
                  dest="nopointdata", help="Do NOT make point data (use data already created)", \
                  default=False)
parser.add_option("--notrans", action="store_true", dest="notrans", default=False, \
                  help='Do not perform transient simulation (spinup only)')
parser.add_option("--nyears_final_spinup", dest="nyears_final_spinup", default='200', \
                  help="base no. of years for final spinup")
parser.add_option("--runroot", dest="runroot", default="", \
                  help="Directory where the run would be created")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup",default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--parm_file", dest="parm_file", default="", \
                  help = 'parameter file to use')
parser.add_option("--parm_file_P", dest="parm_file_P", default="", \
                  help = 'parameter file to use')
parser.add_option("--parm_vals", dest="parm_vals", default="", \
                  help = 'User specified parameter values')
parser.add_option("--postproc_file", dest="postproc_file", default="postproc_vars", \
                  help = 'File for ensemble post processing')
parser.add_option("--nopftdyn", action="store_true", dest="nopftdyn", \
                      default=False, help='Do not use dynamic PFT file')
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
parser.add_option("--ECA", action="store_true", dest="eca", default=False, \
                  help = 'Use ECA compset')
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
parser.add_option("--cruncep", dest="cruncep", default=False, \
                  action="store_true", help = 'Use CRU-NCEP meteorology')
parser.add_option("--gswp3", dest="gswp3", default=False, \
                  action="store_true", help = 'Use GSWP3 meteorology')
parser.add_option("--surfdata_grid", dest="surfdata_grid", default=False, \
                  help = 'Use gridded surface data instead of site data', action="store_true")
parser.add_option("--siteparms",dest = "siteparms", default=False, \
                  action="store_true", help = 'Use default PFT parameters')
parser.add_option("--cpl_bypass", dest = "cpl_bypass", default=False, \
                  help = "Bypass coupler", action="store_true")
parser.add_option("--spinup_vars", dest = "spinup_vars", default=False, \
                  help = "limit output variables for spinup", action="store_true")
parser.add_option("--trans_varlist", dest = "trans_varlist", default='', help = "Transient outputs")
parser.add_option("--cn_only", dest="cn_only", default=False, \
                  help='Carbon/Nitrogen only (saturated P)', action ="store_true")
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
parser.add_option("--mod_parm_file_P", dest="mod_parm_file_P", default='', \
                  help = "adding the path to the modified parameter file")
parser.add_option("--walltime", dest="walltime", default=6, \
                  help = "desired walltime for each job (hours)")

(options, args) = parser.parse_args()

#------------ define function for pbs submission

def submit(fname, project='', submit_type='qsub', job_depend=''):
    job_depend_flag = ' -W depend=afterok:'
    if (submit_type == 'sbatch' and project != ''):
        submit_type = submit_type+' --account '+project
    if ('sbatch' in submit_type):
	job_depend_flag = ' --dependency=afterok:'
    if (job_depend != ''):
        os.system(submit_type+job_depend_flag+job_depend+' '+fname+' > temp/jobinfo')
    else:
        os.system(submit_type+' '+fname+' > temp/jobinfo')
    myinput = open('temp/jobinfo')
    for s in myinput:
        thisjob = re.search('[0-9]+', s).group(0)
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
elif (options.machine == 'edison' or 'cori' in options.machine):
    ccsm_input = '/project/projectdirs/acme/inputdata'

if (options.compiler != ''):
    if (options.machine == 'titan'):
        options.compiler = 'pgi'
    if (options.machine == 'eos' or options.machine == 'edison' or 'cori' in options.machine):
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

csmdir = os.getcwd()+'/'+options.csmdir

#case run and case root directories
myproject=''
if (options.runroot == '' or (os.path.exists(options.runroot) == False)):
    myuser = getpass.getuser()
    if (options.machine == 'titan' or options.machine == 'eos'):
        myinput = open('/ccs/home/'+myuser+'/.cesm_proj','r')
        for s in myinput:
	   myproject=s[:-1]
        runroot='/lustre/atlas/scratch/'+myuser+'/'+myproject
    elif (options.machine == 'cades'):
        runroot='/lustre/or-hydra/cades-ccsi/scratch/'+myuser
    elif ('cori' in options.machine):
        runroot='/global/cscratch1/sd/'+myuser
        myinput = open(os.environ.get('HOME')+'/.cesm_proj','r')
        for s in myinput:
           myproject=s[:-1] 
        print('Project = '+myproject)
    else:
        runroot = csmdir+'/run'
else:
    runroot = os.path.abspath(options.runroot)
if (options.caseroot == '' or (os.path.exists(options.caseroot) == False)):
    caseroot = os.path.abspath(csmdir+'/cime/scripts')
else:
    caseroot = os.path.abspath(options.caseroot)


sitenum=0
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

mysites = options.site.split(',')
timezone = 0
for row in AFdatareader:
    if (row[0] in mysites) or ('all' in mysites and row[0] !='site_code' \
                                      and row[0] != ''):
        site      = row[0]
        if (sitenum == 0):
            firstsite=site
        site_lat  = row[4]
        site_lon  = row[3]
        if (options.cruncep or options.gswp3):
                startyear = 1901
                endyear = 1920
        else:
            startyear = int(row[6])
            endyear   = int(row[7])
        #if (options.diags): 
            timezone = int(row[9])

        site_endyear = int(row[7])
        ncycle   = endyear-startyear+1   #number of years in met cycle
        ny_ad = options.ny_ad
	ny_fin = options.nyears_final_spinup
        if (int(options.ny_ad) % ncycle != 0):
          #AD spinup and final spinup lengths must be multiples of met data cyle.
          ny_ad = str(int(ny_ad) + ncycle - (int(ny_ad) % ncycle))
        if (int(options.nyears_final_spinup) % ncycle !=0):
          ny_fin = str(int(ny_fin) + ncycle - (int(ny_fin) % ncycle))

        if (options.nyears_transient == -1):
          translen = endyear-1850+1        #length of transient run
	  if (options.cpl_bypass and (options.cruncep or options.gswp3)):
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
            options.sitegroup+' --machine '+options.machine
        if (srcmods != ''):
            srcmods    = os.path.abspath(srcmods)
            basecmd = basecmd+' --srcmods_loc '+srcmods
        if (mycaseid != ''):
            basecmd = basecmd+' --caseidprefix '+mycaseid
        if (options.parm_file != ''):
            basecmd = basecmd+' --parm_file '+options.parm_file
        if (options.parm_file_P != ''):
            basecmd = basecmd+' --parm_file_P '+options.parm_file_P
        if (options.parm_vals != ''):
            basecmd = basecmd+' --parm_vals '+options.parm_vals
        if (options.clean_build):
            basecmd = basecmd+' --clean_build '
        if (options.metdir !='none'):
            basecmd = basecmd+' --metdir '+options.metdir
        if (options.C13):
            basecmd = basecmd+' --C13 '
        if (options.C14):
            basecmd = basecmd+' --C14 '
        if (options.debug):
            basecmd = basecmd+' --debug'
        if (options.ninst > 1):
            basecmd = basecmd+' --ninst '+str(options.ninst)
        if (options.nofire):
            basecmd = basecmd+' --nofire'
        if (options.harvmod):
            basecmd = basecmd+' --harvmod'
        if (options.humhol):
            basecmd = basecmd+' --humhol'
        if (options.marsh):
            basecmd = basecmd+' --marsh'
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
        if (options.gswp3):
            basecmd = basecmd+' --gswp3'
        if (options.surfdata_grid):
            basecmd = basecmd+' --surfdata_grid'
        if (options.ensemble_file != ''):   
            basecmd = basecmd+' --ensemble_file '+options.ensemble_file
            basecmd = basecmd+' --parm_list '+options.parm_list
        if (options.archiveroot !=''):
            basecmd = basecmd+' --archiveroot '+options.archiveroot
        if (options.mod_parm_file !=''):
            basecmd = basecmd+' --mod_parm_file '+options.mod_parm_file
        if (options.mod_parm_file_P !=''):
            basecmd = basecmd+' --mod_parm_file_P '+options.mod_parm_file_P
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

        #ECA or CTC
        if (options.cn_only):
            nutrients = 'CN'
        else: 
            nutrients = 'CNP'
        if (options.centbgc):
            decomp_model = 'CNT'
        else:
            decomp_model = 'CTC'
        if (options.eca):
	    mymodel = nutrients+'ECA'+decomp_model
        else:
            mymodel = nutrients+'RD'+decomp_model
        mymodel_adsp = mymodel.replace('CNP','CN')

        #AD spinup
        cmd_adsp = basecmd+' --ad_spinup --nyears_ad_spinup '+ \
            str(ny_ad)+' --align_year '+str(year_align+1)
        if (int(options.hist_mfilt_spinup) == -999):
            cmd_adsp = cmd_adsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
            str((endyear-startyear+1)*8760)
        else:
            cmd_adsp = cmd_adsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
                   +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)
        if (sitenum > 0):
            cmd_adsp = cmd_adsp+' --exeroot '+ad_exeroot+' --no_build'
        if (options.cpl_bypass):
            cmd_adsp = cmd_adsp+' --compset ICB1850'+mymodel_adsp+'BC'
            ad_case = site+'_ICB1850'+mymodel_adsp+'BC'
        else:
            cmd_adsp = cmd_adsp+' --compset I1850'+mymodel_adsp+'BC'
            ad_case = site+'_I1850'+mymodel_adsp+'BC'
        if (options.noad == False):
	    ad_case = ad_case+'_ad_spinup'
        if (options.makemet):
            cmd_adsp = cmd_adsp+' --makemetdat'
        if (options.spinup_vars):
            cmd_adsp = cmd_adsp+' --spinup_vars'
        if (mycaseid != ''):
            ad_case = mycaseid+'_'+ad_case
        if (sitenum == 0):
            ad_exeroot = os.path.abspath(runroot+'/'+ad_case+'/bld')

        #final spinup
        if mycaseid !='':
            basecase=mycaseid+'_'+site
            if (options.cpl_bypass):
                basecase = basecase+'_ICB1850'+mymodel+'BC'
            else: 
                basecase = basecase+'_I1850'+mymodel+'BC'
        else:
            if (options.cpl_bypass):
                basecase=site+'_ICB1850'+mymodel+'BC'
            else:
                basecase=site+'_I1850'+mymodel+'BC'
        if (options.noad):
            cmd_fnsp = basecmd+' --run_units nyears --run_n '+str(fsplen)+' --align_year '+ \
                       str(year_align+1)+' --coldstart'
            if (sitenum > 0):
                cmd_fnsp = cmd_fnsp+' --exeroot '+ad_exeroot+' --no_build'
        else:
            cmd_fnsp = basecmd+' --finidat_case '+ad_case+ \
                       ' --finidat_year '+str(int(ny_ad)+1)+' --run_units nyears --run_n '+ \
                       str(fsplen)+' --align_year '+str(year_align+1)+' --no_build' + \
                       ' --exeroot '+ad_exeroot+' --nopointdata'
        if (int(options.hist_mfilt_spinup) == -999):
            cmd_fnsp = cmd_fnsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
            str((endyear-startyear+1)*8760)
        else:
            cmd_fnsp = cmd_fnsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
                   +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)
        if (options.cpl_bypass):
            cmd_fnsp = cmd_fnsp+' --compset ICB1850'+mymodel+'BC'
        else:
            cmd_fnsp = cmd_fnsp+' --compset I1850'+mymodel+'BC'
        if (options.spinup_vars):
		cmd_fnsp = cmd_fnsp+' --spinup_vars'
        if (options.ensemble_file != '' and options.notrans):	
                cmd_fnsp = cmd_fnsp + ' --postproc_file '+options.postproc_file

        #transient
        cmd_trns = basecmd+' --finidat_case '+basecase+ \
            ' --finidat_year '+str(fsplen+1)+' --run_units nyears' \
            +' --run_n '+str(translen)+' --align_year '+ \
            str(year_align+1850)+' --hist_nhtfrq '+ \
            options.hist_nhtfrq+' --hist_mfilt '+options.hist_mfilt+' --no_build' + \
            ' --exeroot '+ad_exeroot+' --nopointdata'
        if (options.cpl_bypass):
            cmd_trns = cmd_trns+' --compset ICB20TR'+mymodel+'BC'
        else:
            cmd_trns = cmd_trns+' --compset I20TR'+mymodel+'BC'
        if (options.spinup_vars):
            cmd_trns = cmd_trns + ' --spinup_vars'
        if (options.trans_varlist != ''):
            cmd_trns = cmd_trns + ' --trans_varlist '+options.trans_varlist
        if (options.ensemble_file != ''):  #Transient post-processing
            cmd_trns = cmd_trns + ' --postproc_file '+options.postproc_file
        if (options.diags):
            cmd_trns = cmd_trns + ' --diags'
        if (not options.nofire):
            #Turn wildfire off in transient simulations (disturbances are known)
            cmd_trns = cmd_trns + ' --nofire'
        #transient phase 2 (CRU-NCEP only, without coupler bypass)
        if ((options.cruncep or options.gswp3) and not options.cpl_bypass):
            basecase=basecase.replace('1850','20TR')+'_phase1'
            thistranslen = site_endyear - 1921 + 1
            cmd_trns2 = basecmd+' --trans2 --finidat_case '+basecase+ \
                ' --finidat_year 1921 --run_units nyears --branch ' \
                +' --run_n '+str(thistranslen)+' --align_year 1921'+ \
                ' --hist_nhtfrq '+options.hist_nhtfrq+' --hist_mfilt '+ \
                options.hist_mfilt+' --no_build'+' --exeroot '+ad_exeroot + \
                ' --compset I20TR'+mymodel+'BC --nopointdata'
            print(cmd_trns2)

#---------------------------------------------------------------------------------

        #set site environment variable
        os.environ['SITE']=site
        basecase = site
        if (mycaseid != ''):
                basecase = mycaseid+'_'+site
        os.system('mkdir -p temp')

        #If not the first site, create point data here
        if ((sitenum > 0) and not options.nopointdata):
                print 'Creating point data for '+site
                ptcmd = 'python makepointdata.py --caseroot '+caseroot+ \
                        ' --site '+site+' --sitegroup '+options.sitegroup+ \
                        ' --csmdir '+csmdir+' --ccsm_input '+ccsm_input
                os.system(ptcmd)

        #Build Cases
        print('\nSetting up ad_spinup case\n')
        if (options.noad == False):
            if (sitenum == 0):
                ad_case_firstsite = ad_case
                os.system(cmd_adsp)
            else:
                ptcmd = 'python case_copy.py --runroot '+runroot+' --case_copy '+ \
                        ad_case_firstsite+' --site_orig '+firstsite +\
                        ' --site_new '+site+' --nyears '+str(ny_ad)+' --spin_cycle ' \
                        +str(endyear-startyear+1)
                os.system(ptcmd)

        print('\nSetting up final spinup case\n')
        if (sitenum == 0):
            fin_case_firstsite = ad_case_firstsite.replace('_ad_spinup','')
            if (nutrients == 'CNP'):
                fin_case_firstsite = fin_case_firstsite.replace('1850CN','1850CNP')
            os.system(cmd_fnsp)
        else:
            ptcmd = 'python case_copy.py --runroot '+runroot+' --case_copy '+ \
                    fin_case_firstsite+' --site_orig '+firstsite +\
                    ' --site_new '+site+' --nyears '+str(ny_fin)+' --finidat_year ' \
                    +str(int(ny_ad)+1)+' --spin_cycle '+str(endyear-startyear+1)
            os.system(ptcmd)
        if (options.notrans == False):
            print('\nSetting up transient case\n')
            if (sitenum == 0):
                tr_case_firstsite = fin_case_firstsite.replace('1850','20TR')
                os.system(cmd_trns)
            else:
                 ptcmd = 'python case_copy.py --runroot '+runroot+' --case_copy '+ \
                        tr_case_firstsite+' --site_orig '+firstsite +\
                        ' --site_new '+site+' --finidat_year '+str(int(ny_fin)+1)+ \
                        ' --nyears '+str(translen)
                 os.system(ptcmd)
            if ((options.cruncep or options.gswp3) and not options.cpl_bypass):
                 print('\nSetting up transient case phase 2\n')
                 os.system(cmd_trns2)


                 
        #Create .pbs scripts
        case_list = []
        if (options.noad == False):
            case_list.append('ad_spinup')
            case_list.append('iniadjust')
        case_list.append('fn_spinup')
        case_list.append('spinup_diags')
        if (options.notrans == False):
            case_list.append('transient')
            case_list.append('trans_diags')
        
        for c in case_list:
            
            mysubmit_type = 'qsub'
            groupnum = sitenum/32
            if ('cori' in options.machine or options.machine == 'edison'):
                mysubmit_type = 'sbatch'
            if ((sitenum % 32) == 0):
                input = open(caseroot+'/'+ad_case_firstsite+'/.case.run')
                output = open('./temp/'+c+'_group'+str(groupnum)+'.pbs','w')
                for s in input:
                    if ("perl" in s or "python" in s):
                        output.write("#!/bin/csh -f\n")
                        timestr=str(int(float(options.walltime)))+':'+str(int((float(options.walltime)- \
                                     int(float(options.walltime)))*60))+':00'
                        if (mysubmit_type == 'qsub'):
                            output.write('#PBS -l walltime='+timestr+'\n')
                        else:
                            output.write('#SBATCH --time='+timestr+'\n')
                            if ('edison' in options.machine or 'cori' in options.machine):
                                if (options.debug):
                                    output.write('#SBATCH --partition=debug\n')
                                else:
                                    output.write('#SBATCH --partition=regular\n')
                    elif ("#" in s and "ppn" in s):
                        if ('cades' in options.machine):
                            #if ('diags' in c or 'iniadjust' in c):
                            #    output.write("#PBS -l nodes=1:ppn=1\n")
                            #else:
                            output.write("#PBS -l nodes=1:ppn=32\n")
                        else:
                            output.write("#PBS -l nodes=1\n")
                    elif ("#!" in s or "#PBS" in s or "#SBATCH" in s):
                        output.write(s.replace(firstsite,site))
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
                if (options.machine == 'edison' or 'cori' in options.machine):
                     output.write('module unload python\n')
                     output.write('module unload scipy\n')
                     output.write('module unload numpy\n')
                     output.write('module load python/2.7-anaconda\n')
                     output.write('module load nco\n')     
            else:
                output = open('./temp/'+c+'_group'+str(groupnum)+'.pbs','a')   
                
            modelst = 'I1850'+mymodel+'BC'
            if (options.cpl_bypass):
                modelst = 'ICB1850'+mymodel+'BC'

            basecase = site
            if (mycaseid != ''):
                basecase = mycaseid+'_'+site
            if (sitenum == 0 and 'ad_spinup' in c):
                output.write("cd "+caseroot+'/'+basecase+"_"+modelst.replace('CNP','CN')+"_ad_spinup/\n")
                output.write("./case.submit --no-batch &\n")
            elif ('ad_spinup' in c):
                output.write("cd "+runroot+'/'+basecase+"_"+modelst.replace('CNP','CN')+"_ad_spinup/run\n")
                output.write(runroot+'/'+ad_case_firstsite+'/bld/acme.exe &\n')
            if ('iniadjust' in c):
                output.write("cd "+os.path.abspath(".")+'\n')
                if (options.centbgc):
                    output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+ \
                                 '/'+ad_case+'/run/ --casename '+ ad_case+' --restart_year '+ \
                                 str(int(ny_ad)+1)+' --BGC\n')
                else:
                    output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+ \
                                 '/'+ad_case+'/run/ --casename '+ad_case+' --restart_year '+ \
                                 str(int(ny_ad)+1)+'\n')
            if ('spinup_diags' in c):
                 if (options.cpl_bypass):
                     mycompset = 'ICB1850'+mymodel+'BC'
                 else:
                     mycompset = 'I1850'+mymodel+'BC'
                 output.write("cd "+os.path.abspath(".")+'\n')
                 output.write("python plotcase.py --site "+site+" --spinup --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE --ylog --csmdir "+os.path.abspath(runroot)+ \
                               " --pdf\n")
                 output.write("python plotcase.py --site "+site+" --spinup --compset "+mycompset \
                               +" --case "+mycaseid+" --vars TLAI,NPP,GPP,TOTVEGC,TOTSOMC " \
	                       +" --csmdir "+os.path.abspath(runroot)+" --pdf\n")
                 output.write("scp -r ./plots/"+mycaseid+" acme-webserver.ornl.gov:~/www/single_point/plots\n")
            if (sitenum == 0 and 'fn_spinup' in c):
                output.write("cd "+caseroot+'/'+basecase+"_"+modelst+"\n")
                output.write('./case.submit --no-batch &\n')
            elif ('fn_spinup' in c):
                output.write("cd "+runroot+'/'+basecase+"_"+modelst+"/run\n")
                output.write(runroot+'/'+ad_case_firstsite+'/bld/acme.exe &\n')
            if (sitenum == 0 and 'transient' in c):
                output.write("cd "+caseroot+'/'+basecase+"_"+modelst.replace('1850','20TR')+"\n")
                output.write('./case.submit --no-batch &\n')
            elif ('transient' in c):
                output.write("cd "+runroot+'/'+basecase+"_"+modelst.replace('1850','20TR')+"/run\n")
                output.write(runroot+'/'+ad_case_firstsite+'/bld/acme.exe &\n')
            if ('trans_diags' in c):
                 if (options.cpl_bypass):
                     mycompset = 'ICB20TR'+mymodel+'BC'
                 else:
                     mycompset = 'I20TR'+mymodel+'BC'
                 output2 = open('./temp/transdiag_'+site+'.csh','w')
                 #Average seasonal cycle
                 output.write("cd "+os.path.abspath(".")+'\n')
                 output2.write("cd "+os.path.abspath(".")+'\n')

                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE,GPP,EFLX_LH_TOT,FSH,FPG,FPG_P,"+ \
                               "FPI,FPI_P,NPP,QOVER --csmdir "+os.path.abspath(runroot)+ \
                               " --obs --seasonal --pdf --yend "+str(site_endyear)+"\n")
                 #Diurnal cycle (JF)
                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE,GPP,EFLX_LH_TOT,FSH --csmdir "+ \
                              os.path.abspath(runroot)+' --h1 --timezone '+str(timezone)+ \
                               " --obs --diurnal --dstart 1 --dend 59 --pdf --yend "+str(site_endyear)+"\n")
                 #Diurnal cycle (MAM)
                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE,GPP,EFLX_LH_TOT,FSH --csmdir "+ \
                              os.path.abspath(runroot)+' --h1 --timezone '+str(timezone)+ \
                               " --obs --diurnal --dstart 60 --dend 151 --pdf --yend "+str(site_endyear)+"\n")
                 #Diurnal cycle (JJA)
                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE,GPP,EFLX_LH_TOT,FSH --csmdir "+ \
                              os.path.abspath(runroot)+' --h1 --timezone '+str(timezone)+ \
                               " --obs --diurnal --dstart 152 --dend 243 --pdf --yend "+str(site_endyear)+"\n")
                 #Diurnal cycle (SON)
                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE,GPP,EFLX_LH_TOT,FSH --csmdir "+ \
                              os.path.abspath(runroot)+' --h1 --timezone '+str(timezone)+ \
                               " --obs --diurnal --dstart 244 --dend 334 --pdf --yend "+str(site_endyear)+"\n")
                 #Interannual variability
                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE,GPP,EFLX_LH_TOT,FSH --csmdir "+ \
                               os.path.abspath(runroot)+" --obs --h4 --pdf\n")
                 #monthly data
                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars NEE,GPP,EFLX_LH_TOT,FSH,FPG,FPG_P,"+ \
                               "FPI,FPI_P,NPP,QOVER --csmdir "+os.path.abspath(runroot)+ \
                               " --obs --pdf --yend "+str(site_endyear)+"\n")
                 #1850-present
                 output2.write("python plotcase.py --site "+site+" --compset "+mycompset \
                               +" --case "+mycaseid+" --vars TOTLITC,CWDC,TOTVEGC,TOTSOMC --csmdir " \
                               +os.path.abspath(runroot)+" --ystart 1850 --yend "+str(site_endyear)+" --h4 --pdf\n")
                 output2.close()
                 os.system('chmod u+x '+'./temp/transdiag_'+site+'.csh')
                 output.write('./temp/transdiag_'+site+'.csh &\n')
            output.write('sleep 1\n')
            output.close()

        #if ensemble simulations requested, submit jobs created by pointclm.py in correct order
        if (options.ensemble_file != ''):
            cases=[]
            #build list of cases for fullrun
            if (options.noad == False):
                cases.append(basecase+'_'+modelst.replace('CNP','CN')+'_ad_spinup')
            cases.append(basecase+'_'+modelst)
            if (options.notrans == False):
                cases.append(basecase+'_'+modelst.replace('1850','20TR'))
            job_depend_run=''    
            for thiscase in cases:
                job_depend_run = submit('temp/ensemble_run_'+thiscase+'.pbs',job_depend= \
                                        job_depend_run, project=myproject, submit_type=mysubmit_type)
        #else:  #submit single job
        #    job_fullrun = submit('temp/site_fullrun.pbs', submit_type=mysubmit_type)
        sitenum=sitenum+1

#Submit PBS scripts for multi-site simualtions on 1 node
if (options.ensemble_file == ''):
    for g in range(0,groupnum+1):
        job_depend_run=''
        for thiscase in case_list:
            output = open('./temp/'+thiscase+'_group'+str(g)+'.pbs','a')
            output.write('wait\n')
            if ('trans_diags' in thiscase):
                output.write("scp -r ./plots/"+mycaseid+" acme-webserver.ornl.gov:~/www/single_point/plots\n")
            output.close()
            job_depend_run = submit('temp/'+thiscase+'_group'+str(g)+'.pbs',job_depend= \
                                    job_depend_run, project=myproject, submit_type=mysubmit_type)
