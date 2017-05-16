#!/usr/bin/env python

import os, sys, csv
from optparse import OptionParser

parser = OptionParser();

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--machine", dest="machine", default = 'oic2', \
                  help = "machine to use (default = oic2)")
parser.add_option("--compiler", dest="compiler", default = 'gnu', \
                  help = "compiler to use (pgi*, gnu)")
parser.add_option("--mpilib", dest="mpilib", default="mpi-serial", \
                      help = "mpi library (openmpi*, mpich, ibm, mpi-serial)")
parser.add_option("--csmdir", dest="csmdir", default='../../..', \
                  help = "base CESM directory (default = ../../..)")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--clm40", dest="clm40", action="store_true", \
                  default=False, help="Use CLM 4.0 code")
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--nyears_final_spinup", dest="nyears_final_spinup", default='1000', \
                  help="base no. of years for final spinup")
parser.add_option("--clean_build", action="store_true", default=False, \
                  help="Perform a clean build")
parser.add_option("--hist_mfilt", dest="hist_mfilt", default="365", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_nhtfrq", dest="hist_nhtfrq", default="-24", \
                  help = 'output file timestep (transient only)')
parser.add_option("--parm_file", dest="parm_file", default="", \
                  help = 'parameter file to use')
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
parser.add_option("--harvmod", action="store_true", dest='harvmod', default=False, \
                    help="turn on harvest modification:  All harvest at first timestep")
parser.add_option("--bulk_denitrif", dest="nitrif", default=False, \
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
parser.add_option("--refcase", dest="refcase" , default='none', \
                  help = 'Use already compiled CLM case')
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

(options, args) = parser.parse_args()

ccsm_input = options.ccsm_input
mycaseid   = options.mycaseid
srcmods    = options.srcmods_loc

#get start and year of input meteorology from site data file
PTCLMfiledir = options.ccsm_input+'/lnd/clm2/PTCLM'
fname = PTCLMfiledir+'/'+options.sitegroup+'_sitedata.txt'
AFdatareader = csv.reader(open(fname, "rt"))

translen = int(options.nyears_transient)
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
	  if (options.cpl_bypass):
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
        if (options.refcase != 'none'):
            basecmd = basecmd+' --refcase '+options.refcase
        if (options.nofire):
            basecmd = basecmd+' --nofire'
        if (options.harvmod):
            basecmd = basecmd+' --harvmod'
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
        if (options.MICROBE):
	    basecmd = basecmd+' --MICROBE'
        if (options.cruncep):
            basecmd = basecmd+' --cruncep'
        if (options.surfdata_grid):
            basecmd = basecmd+' --surfdata_grid'
        basecmd = basecmd + ' --np '+str(options.np)
        basecmd = basecmd + ' --tstep '+str(options.tstep)
        basecmd = basecmd + ' --co2_file '+options.co2_file
        basecmd = basecmd + ' --compiler '+options.compiler
        basecmd = basecmd + ' --mpilib '+options.mpilib

#----------------------- build commands for runCLM.py -----------------------------

        #AD spinup
        cmd_adsp = basecmd+' --ad_spinup --nyears_ad_spinup '+ \
            str(options.ny_ad)+' --hist_mfilt 1 --hist_nhtfrq -'+ \
            str((endyear-startyear+1)*8760)+' --align_year '+str(year_align+1)
        if (options.clm40):
            cmd_adsp = cmd_adsp+' --compset I1850CN'
        elif (options.cpl_bypass):
            cmd_adsp = cmd_adsp+' --compset I1850CLM45CBCN'
        if (options.makemet):
            cmd_adsp = cmd_adsp+' --makemetdat'
        if (options.spinup_vars):
            cmd_adsp = cmd_adsp+' --spinup_vars'
        #exit spinup (CLM 4.0 only)
        if (options.clm40):
            ad_case = site+'_I1850CN_ad_spinup'
            if (mycaseid != ''):
                ad_case = mycaseid+'_'+ad_case
            cmd_exsp = basecmd+' --exit_spinup --compset I1850CN '+ \
                '--finidat_case '+ad_case+' --finidat_year '+str(options.ny_ad)+ \
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
                  basecase = basecase+'_I1850CLM45CBCN'
                else: 
                  basecase = basecase+'_I1850CLM45CN'
        else:
            if (options.clm40):
                basecase = site+'_I1850CN'
            else:
                if (options.cpl_bypass):
                  basecase=site+'_I1850CLM45CBCN'
                else:
                  basecase=site+'_I1850CLM45CN'
        if (options.clm40):
           cmd_fnsp = basecmd+' --finidat_case '+basecase+'_exit_spinup '+ \
                '--finidat_year '+str(int(options.ny_ad)+2)+' --run_units nyears --run_n '+ \
                str(fsplen)+' --hist_mfilt 1 --hist_nhtfrq -'+ \
                str((endyear-startyear+1)*8760)
        else:
            cmd_fnsp = basecmd+' --finidat_case '+basecase+'_ad_spinup '+ \
                '--finidat_year '+str(int(options.ny_ad)+1)+' --run_units nyears --run_n '+ \
                str(fsplen)+' --hist_mfilt 1 --hist_nhtfrq -'+ \
                str((endyear-startyear+1)*8760)+' --align_year '+str(year_align+1)
        if (options.clm40):
            cmd_fnsp = cmd_fnsp+' --compset I1850CN'
        elif (options.cpl_bypass):
            cmd_fnsp = cmd_fnsp+' --compset I1850CLM45CBCN'
        if (options.spinup_vars):
		cmd_fnsp = cmd_fnsp+' --spinup_vars'
        #transient
        cmd_trns = basecmd+' --finidat_case '+basecase+ \
            ' --finidat_year '+str(fsplen+1)+' --run_units nyears' \
            +' --run_n '+str(translen)+' --align_year '+ \
            str(year_align+1850)+' --hist_nhtfrq '+ \
            options.hist_nhtfrq+' --hist_mfilt '+options.hist_mfilt
        if (options.clm40):
             cmd_trns = cmd_trns+' --compset I20TRCN'
        elif (options.cpl_bypass):
            cmd_trns = cmd_trns+' --compset I20TRCLM45CBCN'
        else:
            cmd_trns = cmd_trns+' --compset I20TRCLM45CN'
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
                options.hist_mfilt
            if (options.clm40):
                cmd_trns2 = cmd_trns2+' --compset I20TRCN'
            else:
                cmd_trns2 = cmd_trns2+' --compset I20TRCLM45CN'
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
        os.system(cmd_adsp)
        if (options.clm40):
            print("\nSetting up exit_spinup case\n")
            os.system(cmd_exsp)
        print('\nSetting up final spinup case\n')
        os.system(cmd_fnsp)
        print('\nSetting up transient case\n')
        os.system(cmd_trns)
        if (options.cruncep and not options.cpl_bypass):
             print('\nSetting up transient case phase 2\n')
             os.system(cmd_trns2)
        
        output = open('./temp/site_fullrun.pbs','w')
        #make site-specific pbs script

        if (options.cpl_bypass):
          input = open('../../scripts/'+basecase+"_I1850CLM45CBCN"+'/'+basecase+"_I1850CLM45CBCN"+'.run')
        else:
          input = open('../../scripts/'+basecase+"_I1850CLM45CN"+'/'+basecase+"_I1850CLM45CN"+'.run')
        for s in input:
            if ("perl" in s):
                output.write("#!/bin/csh -f\n")
            elif ("#PBS" in s or "#!" in s):
                output.write(s)
        input.close()
        output.write("\n")
        
        if (options.cpl_bypass):
          modelst = 'CLM45CBCN'
        else:
          modelst = 'CLM45CN'

        basecase = site
        if (mycaseid != ''):
                basecase = mycaseid+'_'+site
        output.write("cd "+os.path.abspath("../../scripts/"+basecase+"_I1850"+modelst+"_ad_spinup\n"))
        output.write("./"+basecase+"_I1850"+modelst+"_ad_spinup.run\n")
        output.write("cd "+os.path.abspath("../../scripts/"+basecase+"_I1850"+modelst+"\n"))
        output.write("./"+basecase+"_I1850"+modelst+".run\n")
        output.write("cd "+os.path.abspath("../../scripts/"+basecase+"_I20TR"+modelst+"\n"))
        output.write("./"+basecase+"_I20TR"+modelst+".run\n")
        output.close()


        #submit
        os.system('qsub temp/site_fullrun.pbs')

