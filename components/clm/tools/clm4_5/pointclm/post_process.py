#!/usr/bin/env python

import netcdf_functions as nffun
import os, sys, csv, time, math, numpy
from optparse import OptionParser
#from Numeric import *

parser = OptionParser()

parser.add_option('--UQ_rundir', dest='UQ_rundir', default="../../../../../run/UQ", \
                      help="location of UQ run directory")
parser.add_option('--case', dest='case', default='' , \
                      help="full case name")
parser.add_option('--ninst', dest='ninst', default=1, \
                      help="number of land model instances")
parser.add_option('--start_year', dest="start_year", default=1, \
                      help ="start year for post-processing")
parser.add_option('--end_year', dest='end_year', default=1, \
                      help='end year for post processing')
parser.add_option('--start_day', dest="start_day", default=1, \
                      help ="start year for post-processing")
parser.add_option('--end_day', dest='end_day', default=365, \
                      help='end year for post processing')
parser.add_option('--start_month', dest="start_month", default=1, \
                      help ="start year for post-processing")
parser.add_option('--end_month', dest='end_month', default=12, \
                      help='end year for post processing')
parser.add_option('--vars', dest='vars', default="", \
                      help='variables to post process')
parser.add_option('--ensnum_start', dest='ensnum_start', default=1, \
                      help='Beginning ensemble member number to start')
parser.add_option('--ens_size', dest="ens_size", default=-1, \
                      help='Ensemble size to process')
parser.add_option('--avpd', dest="avpd", default=1, help='Averaging period')

(options, args) = parser.parse_args()

casedir = options.UQ_rundir+'/'+options.case
if (options.ninst == 1):
    input_fname = casedir+'/g00001/lnd_in'
else:
    input_fname = casedir+'/g00001/lnd_in_0001'

#get timestep info from lnd_in file
if (os.path.isfile(input_fname)):
    myinput = open(input_fname,'r')
    for s in myinput:
        if ('hist_mfilt' in s):
            mfilt = s.split()
            h0_mfilt = int(mfilt[2].split(",")[0])
        if ('hist_nhtfrq' in s):
            nhtfrq = s.split()
            h0_nhtfrq = int(nhtfrq[2].split(",")[0])
n_years = int(options.end_year)-int(options.start_year)+1

if (h0_nhtfrq != 0 and (h0_nhtfrq*h0_mfilt % 8760) != 0):
    print ('Must have either monthly or annual history files.')
    sys.exit()

if (int(h0_nhtfrq) == 0):
    n_files = n_years*(options.end_month-options.start_month+1)
else:
    n_files = max((-8760*n_years)/(int(h0_nhtfrq)*int(h0_mfilt)),1)
var_list = options.vars.split(",")

os.system('mkdir -p ./processed/'+options.case)
if (int(options.ens_size) == -1):
    print 'Ensemble size not specified'
    options.ens_size = len([f for f in os.listdir(casedir) if ((len(f) == 6) and 'g' in f)])
    print 'Processing all '+str(options.ens_size)+' ensemble members'

for v in var_list:
   print ('Processing '+v)
   fname_proc = './processed/'+options.case+'/'+v+'_'+str(options.start_year)+'-'+ \
                 options.end_year
   if (options.start_day != 1 or options.end_day != 365):
       fname_proc = fname_proc+'_d'+str(1000+int(options.start_day))[1:]+'-'+ \
                                    str(1000+int(options.end_day))[1:]+'_processed.txt'
   elif (options.start_month != 1 or options.end_month != 12):
       fname_proc = fname_proc+'_m'+str(100+int(options.start_month))[1:]+'-'+ \
                                    str(100+int(options.end_month))[1:]+'_processed.txt'
   else:
       fname_proc = fname_proc+'_processed.txt'
   output = open(fname_proc,'w')
   for n in range(0, int(options.ens_size)):
        sys.stdout.write(str(int((n*100)/int(options.ens_size)))+ ' % complete\r')
        ensdirnum = (int(options.ensnum_start)+n-1)/int(options.ninst)+1
        ensstr = casedir+'/g'+str(100000+ensdirnum)[1:]

        myvar=[]
        for i in range(0,n_files):
            yst = str(10000+int(options.start_year)+(n_years/n_files)*i)[1:]
            if (int(h0_nhtfrq) == 0):
                mst = str(100+options.start_month + (i%n_years))[1:]
                fname = casedir+'/'+ensstr+'/'+options.case+'.clm2.h0.'+yst+'-'+mst+'.nc'
            else:
                fname = casedir+'/'+ensstr+'/'+options.case+'.clm2.h0.'+yst+'-01-01-00000.nc'
            if (os.path.isfile(fname)):
                var = nffun.getvar(fname, v)
                #convert carbon fluxes to gC/m2/day
                if (v == 'GPP' or v == 'NEE' or v == 'NPP' or v == 'ER'):
                    var = var*24*3600
                if (len(var) == 365):
                    for d in range(int(options.start_day)-1,int(options.end_day)):
                        myvar.append(var[d][0])
                else:
                    myvar.append(var)
            else:
                print('Error: '+fname+' does not exist')
                sys.exit()
        myvarst=''
        avpd = int(options.avpd)
        for i in range(0,len(myvar)/avpd):
            vartemp=0.0
            for j in range(0,avpd):
                vartemp = vartemp+myvar[i*avpd+j]/avpd 
            myvarst  = myvarst+(str(float(vartemp)))+' '
        output.write(myvarst+'\n')
   output.close()
