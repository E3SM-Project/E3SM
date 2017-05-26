#!/usr/bin/env python

import os, sys, csv, time, math, numpy
from optparse import OptionParser
#from Numeric import *

parser = OptionParser()

parser.add_option('--UQ_rundir', dest='UQ_rundir', default="", \
                      help="location of UQ run directory")
parser.add_option('--case', dest='case', default='' , \
                      help="full case name")
parser.add_option('--ninst', dest='ninst', default=1, \
                      help="number of land model instances")
parser.add_option('--start_year', dest="start_year", default=1, \
                      help ="start year for post-processing")
parser.add_option('--end_year', dest='end_year', default=1, \
                      help='end year for post processing')
parser.add_option('--vars', dest='vars', default="", \
                      help='variables to post process')
parser.add_option('--ensnum_start', dest='ensnum_start', default=1, \
                      help='Beginning ensemble member number to start')
parser.add_option('--ens_size', dest="ens_size", default=1, \
                      help='Ensemble size to process')
parser.add_option('--nav', dest='nav', default = 1, \
                      help='number of points to average over')

(options, args) = parser.parse_args()

#================= netcdf manipulation functions ===============================#
def getvar(fname, varname):
    usescipy = False
    try:
    	import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"r")
        var = nffile.variables[varname]
        varvals = var[:].copy()    #works for vector only?
        nffile.close()
    else:
    	nffile = netcdf.NetCDFFile(fname,"r")
    	var = nffile.variables[varname]
    	varvals = var.getValue()
    	nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"a")
        var = nffile.variables[varname]
        var[:] = varvals[:]
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"a")
        var = nffile.variables[varname]
        var.assignValue(varvals)
        nffile.close()
    ierr = 0
    return ierr


os.chdir(os.path.abspath(options.UQ_rundir+'/'+options.case))
if (options.ninst == 1):
    input_fname = './g00001/lnd_in'
else:
    input_fname = './g00001/lnd_in_0001'

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
    n_files = n_years*12
else:
    n_files = max((-8760*n_years)/(int(h0_nhtfrq)*int(h0_mfilt)),1)

var_list = options.vars.split(",")

for v in var_list:
   output = open(v+'_processed.txt','w')
   for n in range(0, int(options.ens_size)):
        ensdirnum = (int(options.ensnum_start)+n-1)/int(options.ninst)+1
        ensstr = './g'+str(100000+ensdirnum)[1:]

        myvar=[]
        for i in range(0,n_files):
            if (int(h0_nhtfrq) == 0):
                fff
            else:
                yst = str(10000+int(options.start_year)+(n_years/n_files)*i)[1:]
                fname = ensstr+'/'+options.case+'.clm2.h0.'+yst+'-01-01-00000.nc'
                if (os.path.isfile(fname)):
                    var = getvar(fname, v)
                    #convert carbon fluxes to gC/m2/yr
                    if (v == 'GPP' or v == 'NEE' or v == 'NPP' or v == 'ER'):
                        var = var*24*3600*365
                    myvar.append(var)
                else:
                    print('Warning: '+fname+' does not exist.  Setting to 0')
                    myvar.append(0)
        myvarst=''
        vartemp=0.0
        for i in range(0,len(myvar)):
            vartemp = vartemp+myvar[i]/int(options.nav)
            if ((i+1) % int(options.nav) == 0):
                myvarst  = myvarst+(str(float(vartemp)))+' '
                vartemp=0
        output.write(myvarst+'\n')
   output.close()
