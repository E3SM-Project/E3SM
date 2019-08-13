#!/usr/bin/env python

import os, sys, csv, math
from optparse import OptionParser
import subprocess

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
                  default='../../../../../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--clm40", dest="clm40", action="store_true", \
                  default=False, help="Use CLM 4.0 code")
parser.add_option("--nopointdata", action="store_true", \
                  dest="nopointdata", help="Do NOT make point data (use data already created)", \
                  default=False)
parser.add_option("--noad", action="store_true", dest="noad", default=False, \
                  help='Do not perform ad spinup simulation')
parser.add_option("--notrans", action="store_true", dest="notrans", default=False, \
	          help='Do not perform transient simulation (spinup only)')
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--nyears_final_spinup", dest="nyears_final_spinup", default='1000', \
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
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=600, \
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
                  help = "Bypass couplwer", action="store_true")
parser.add_option("--spinup_vars", dest = "spinup_vars", default=False, help = "limit output variables for spinup", action="store_true")
parser.add_option("--cn_only", dest="cn_only", default=False, help='Carbon/Nitrogen only (saturated P)', action ="store_true")
parser.add_option("--ensemble_file", dest="ensemble_file", default='', \
                  help = 'Parameter sample file to generate ensemble')
parser.add_option("--mc_ensemble", dest="mc_ensemble", default=-1, \
                  help = 'Monte Carlo ensemble (argument is # of simulations)')
parser.add_option("--ng", dest="ng", default=64, \
                  help = 'number of groups to run in ensemble mode')
#Changed by Ming for mesabi
parser.add_option("--archiveroot", dest="archiveroot", default='', \
                  help = "archive root directory only for mesabi")
#Added by Kirk to include the modified parameter file
parser.add_option("--mod_parm_file", dest="mod_parm_file", default='', \
                  help = "adding the path to the modified parameter file")

(options, args) = parser.parse_args()

#------------ define function for pbs submission
def submit(fname, submit_type='qsub', job_depend=''):
    if (job_depend != ''):
        subprocess.call(submit_type+' -W depend=afterok:'+job_depend+' '+fname+' > temp/jobinfo', shell=True)
    else:
        subprocess.call(submit_type+' '+fname+' > temp/jobinfo', shell=True)
    myinput = open('temp/jobinfo')
    for s in myinput:
        thisjob = s.split('.')[0]
    myinput.close()
    os.system('rm temp/jobinfo')
    return thisjob
#----------------------------------------------------------

ccsm_input = options.ccsm_input
mycaseid   = options.mycaseid
srcmods    = options.srcmods_loc

#get start and year of input meteorology from site data file
PTCLMfiledir = options.ccsm_input+'/lnd/clm2/PTCLM'
fname = PTCLMfiledir+'/'+options.sitegroup+'_sitedata.txt'
AFdatareader = csv.reader(open(fname, "rt"))

translen = int(options.nyears_transient)

if (options.bgc):
    mybgc='BGC'
else:
    mybgc='CN'

csmdir = os.getcwd()+'/'+options.csmdir

#case run root directory
if (options.runroot == '' or (os.path.exists(options.runroot) == False)):
    runroot = csmdir+'/run'
else:
    runroot = os.path.abspath(options.runroot)

isfirstsite = True
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
        if (int(options.ny_ad) % ncycle != 0):
          #AD spinup and final spinup lengths must be multiples of met data cyle.
          options.ny_ad = str(int(options.ny_ad) + ncycle - (int(options.ny_ad) % ncycle))
        if (int(options.nyears_final_spinup) % ncycle !=0):
          options.nyears_final_spinup = str(int(options.nyears_final_spinup) + ncycle - \
                (int(options.nyears_final_spinup) % ncycle))

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

        fsplen = int(options.nyears_final_spinup)
 
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
        if (options.caseroot != ''):
            basecmd = basecmd+' --caseroot '+options.caseroot
        if (options.runroot != ''):
            basecmd = basecmd+' --runroot '+options.runroot
        if (options.ensemble_file != ''):
            basecmd = basecmd+' --ensemble_file '+options.ensemble_file
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

#----------------------- build commands for runCLM.py -----------------------------



        #AD spinup
        cmd_adsp = basecmd+' --ad_spinup --nyears_ad_spinup '+ \
            str(options.ny_ad)+' --align_year '+str(year_align+1)
        if (int(options.hist_mfilt_spinup) == -999):
            cmd_adsp = cmd_adsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
            str((endyear-startyear+1)*8760)
        else:
            cmd_adsp = cmd_adsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
                   +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)
        if (not isfirstsite):
	  cmd_adsp = cmd_adsp+' --exeroot '+ad_exeroot+' --no_build'
        if (options.clm40):
            cmd_adsp = cmd_adsp+' --compset I1850'+mybgc
            ad_case = site+'_I1850'+mybgc+'_ad_spinup'
        else:
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
        if (options.clm40):
            cmd_exsp = basecmd+' --exit_spinup --compset I1850'+mybgc+ \
                ' --finidat_case '+ad_case+' --finidat_year '+str(options.ny_ad)+ \
                ' --run_units nyears --run_n 1 --nyears_ad_spinup '+str(options.ny_ad)
	    if (options.spinup_vars):
               cmd_exsp = cmd_exsp + ' --spinup_vars'
        #final spinup
        if mycaseid !='':
            basecase=mycaseid+'_'+site
            if (options.clm40):
                basecase = basecase+'_I1850CN'
            else:
                if (options.cpl_bypass):
                  basecase = basecase+'_I1850CLM45CB'+mybgc
                else: 
                  basecase = basecase+'_I1850CLM45'+mybgc
        else:
            if (options.clm40):
                basecase = site+'_I1850CN'
            else:
                if (options.cpl_bypass):
                  basecase=site+'_I1850CLM45CB'+mybgc
                else:
                  basecase=site+'_I1850CLM45'+mybgc
        if (options.clm40):
           cmd_fnsp = basecmd+' --finidat_case '+basecase+'_exit_spinup '+ \
                '--finidat_year '+str(int(options.ny_ad)+2)+' --run_units nyears --run_n '+ \
                str(fsplen)
        else:
            if (options.noad):
                cmd_fnsp = basecmd+' --run_units nyears --run_n '+str(fsplen)+' --align_year '+ \
                    str(year_align+1)+' --coldstart'
                if (not isfirstsite):
                    cmd_fnsp = cmd_fnsp+' --exeroot '+ad_exeroot+' --no_build'
            else:
                cmd_fnsp = basecmd+' --finidat_case '+basecase+'_ad_spinup '+ \
                    '--finidat_year '+str(int(options.ny_ad)+1)+' --run_units nyears --run_n '+ \
                    str(fsplen)+' --align_year '+str(year_align+1)+' --no_build' + \
                    ' --exeroot '+ad_exeroot
        if (int(options.hist_mfilt_spinup) == -999):
            cmd_fnsp = cmd_fnsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
            str((endyear-startyear+1)*8760)
        else:
            cmd_fnsp = cmd_fnsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
                   +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)
        if (options.clm40):
            cmd_fnsp = cmd_fnsp+' --compset I1850CN'
        elif (options.cpl_bypass):
            cmd_fnsp = cmd_fnsp+' --compset I1850CLM45CB'+mybgc
        if (options.spinup_vars):
		cmd_fnsp = cmd_fnsp+' --spinup_vars'
        #transient
        cmd_trns = basecmd+' --finidat_case '+basecase+ \
            ' --finidat_year '+str(fsplen+1)+' --run_units nyears' \
            +' --run_n '+str(translen)+' --align_year '+ \
            str(year_align+1850)+' --hist_nhtfrq '+ \
            options.hist_nhtfrq+' --hist_mfilt '+options.hist_mfilt+' --no_build' + \
            ' --exeroot '+ad_exeroot
        if (options.clm40):
             cmd_trns = cmd_trns+' --compset I20TRCN'
        elif (options.cpl_bypass):
            cmd_trns = cmd_trns+' --compset I20TRCLM45CB'+mybgc
        else:
            cmd_trns = cmd_trns+' --compset I20TRCLM45'+mybgc
        if (options.spinup_vars):
               cmd_trns = cmd_trns + ' --spinup_vars'
        #transient phase 2 (CRU-NCEP only)
        if (options.cruncep and not options.cpl_bypass):
            basecase=basecase.replace('1850','20TR')+'_phase1'
            thistranslen = site_endyear - 1921 + 1
            cmd_trns2 = basecmd+' --trans2 --finidat_case '+basecase+ \
                ' --finidat_year 1921 --run_units nyears --branch ' \
                +' --run_n '+str(thistranslen)+' --align_year 1921'+ \
                ' --hist_nhtfrq '+options.hist_nhtfrq+' --hist_mfilt '+ \
                options.hist_mfilt+' --no_build'+' --exeroot '+ad_exeroot
            if (options.clm40):
                cmd_trns2 = cmd_trns2+' --compset I20TRCN'
            else:
                cmd_trns2 = cmd_trns2+' --compset I20TRCLM45'+mybgc
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
            os.system(cmd_adsp)
            if (options.clm40):
                print("\nSetting up exit_spinup case\n")
                os.system(cmd_exsp)
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

        #Create a .PBS site fullrun script to launch the full job (all 3 cases)
        if (options.cpl_bypass):
          input = open(csmdir+'/cime/scripts/'+basecase+"_I1850CLM45CB"+mybgc+ \
                  '/'+basecase+"_I1850CLM45CB"+mybgc+'.run')
        else:
          input = open(csmdir+'/cime/scripts/'+basecase+"_I1850CLM45"+mybgc+ \
                           '/'+basecase+"_I1850CLM45"+mybgc+'.run')
        for s in input:
            if ("perl" in s):
                output.write("#!/bin/csh -f\n")
            elif ("#PBS" in s or "#!" in s):
                output.write(s.replace('24:00','72:00'))
            elif ("#SBATCH" in s or "#!" in s):
                output.write(s)
        input.close()
        output.write("\n")
        
        if (options.cpl_bypass):
          modelst = 'CLM45CB'+mybgc
        else:
          modelst = 'CLM45'+mybgc

        basecase = site
        if (mycaseid != ''):
                basecase = mycaseid+'_'+site
        if (options.noad == False):
            output.write("cd "+os.path.abspath(csmdir+"/cime/scripts/"+basecase+"_I1850"+modelst+"_ad_spinup\n"))
            output.write("./"+basecase+"_I1850"+modelst+"_ad_spinup.run\n")
            output.write("cd "+os.path.abspath(".")+'\n')
            if (options.bgc):
                output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+'/'+ad_case+ \
                                 '/run/ --casename '+ ad_case+' --restart_year '+str(int(options.ny_ad)+1)+ \
                                 ' --BGC\n')
            else:
                output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+'/'+ad_case+ \
                                 '/run/ --casename '+ad_case+' --restart_year '+str(int(options.ny_ad)+1)+'\n')
        output.write("cd "+os.path.abspath(csmdir+"/cime/scripts/"+basecase+"_I1850"+modelst+"\n"))
        output.write("./"+basecase+"_I1850"+modelst+".run\n")
        if (options.notrans == False):
            output.write("cd "+os.path.abspath(csmdir+"/cime/scripts/"+basecase+"_I20TR"+modelst+"\n"))
            output.write("./"+basecase+"_I20TR"+modelst+".run\n")
        output.close()


        #if ensemble simulations requested, submit jobs created by pointclm.py in correct order
        if (options.ensemble_file != ''):
            nsamples = 0
            myinput = open(options.csmdir+'/components/clm/tools/clm4_5/pointclm/'+options.ensemble_file)
            #count number of samples
            for s in myinput:
                nsamples = nsamples+1
            myinput.close()
            n_qsub_files = int(math.ceil(nsamples*1.0/int(options.ng)))
            cases=[]
            #build list of cases for fullrun
            if (options.noad == False):
                cases.append(basecase+'_I1850'+modelst+'_ad_spinup')
            cases.append(basecase+'_I1850'+modelst)
            if (options.notrans == False):
                cases.append(basecase+'_I20TR'+modelst)
                
            job_depend_copy = ''
            job_depend_run = ''
            for thiscase in cases:
                for i in range(0,n_qsub_files):
                    if (i == 0):
                        job_depend_copy = submit('temp/ensemble_copy_'+thiscase+'_'+str(1000+i)[1:]+'.pbs', \
                                                     job_depend=job_depend_run)
                        job_depend_run = submit('temp/ensemble_run_'+thiscase+'_'+str(1000+i)[1:]+'.pbs', \
                                                    job_depend=job_depend_copy)
                    else:
                        job_depend_copy = submit('temp/ensemble_copy_'+thiscase+'_'+str(1000+i)[1:]+'.pbs', \
                                                     job_depend=job_depend_copy)
                        job_depend_run = submit('temp/ensemble_run_'+thiscase+'_'+str(1000+i)[1:]+'.pbs', \
                                                    job_depend=job_depend_run)         
        else:  #submit single job
            if ('edison' in options.machine):
                job_fullrun = submit('temp/site_fullrun.pbs', submit_type='sbatch')
            else:
                job_fullrun = submit('temp/site_fullrun.pbs')
        isfirstsite = False
