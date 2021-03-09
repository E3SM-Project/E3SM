#!/usr/bin/python

import os, sys, csv, glob
import numpy, scipy
from scipy.io import netcdf
from optparse import OptionParser
import matplotlib as mpl

def getvar(fname, varname, npf, index, scale_factor):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"r",mmap=False)
        var = nffile.variables[varname]
        varvals = var[0:npf,index].copy() * scale_factor    #works for vector only?
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"r")
        var = nffile.variables[varname]
        varvals = var.getValue()[0:npf,index] * scale_factor
        nffile.close()
    return varvals


parser = OptionParser()
parser.add_option("--ad_Pinit", dest="ad_Pinit", default=False, action="store_true",\
                  help="Initialize AD spinup with P pools and use CNP mode")
parser.add_option("--csmdir", dest="mycsmdir", default='', \
                  help = 'Base CESM directory (default = ..)')
parser.add_option("--cases", dest="mycase", default='', \
                  help = "name of case id prefixs to plot (comma delmited)")
parser.add_option("--compset", dest="compset", default="I20TRCLM45CN", \
                  help = "Compset to plot")
parser.add_option("--titles", dest="titles", default='', \
                  help = "titles of case to plot (for legend)")
parser.add_option("--obs", action="store_true", default=False, \
                  help = "plot observations", dest="myobs")
parser.add_option("--sites", dest="site", default="none", \
                  help = 'site (to plot observations)')
parser.add_option("--timezone", dest="timezone", default=0, \
                  help = 'time zone (relative to UTC')
parser.add_option("--varfile", dest="myvarfile", default='varfile', \
                  help = 'file containing list of variables to plot')
parser.add_option("--vars", dest="myvar", default='', \
                  help="variable to plot (overrides varfile, " \
                  +"sends plot to screen")
parser.add_option("--avpd", dest="myavpd", default=1, \
                  help = 'averaging period in # of output timesteps' \
                  +' (default = 1)')
parser.add_option("--hist_mfilt", dest="myhist_mfilt", default=-999, \
                  help = 'beginning model year to plot')
parser.add_option("--hist_nhtfrq", dest="myhist_nhtfrq", default=-999, \
                  help = 'beginning model year to plot')
parser.add_option("--ystart", dest="myystart", default=1, \
                  help = 'beginning model year to plot')
parser.add_option("--yend", dest="myyend", default=9999, \
                  help = 'final model year to plot')
parser.add_option("--ystart_obs", dest="ystart_obs", default=0, \
                  help = 'beginning model year to plot')
parser.add_option("--yend_obs", dest="yend_obs", default=0, \
                  help = 'final model year to plot')
parser.add_option("--diurnal", dest="mydiurnal", default=False, \
                  action="store_true", help = 'plot diurnal cycle')
parser.add_option("--dstart", dest="dstart", default=1, \
                  help = 'beginning model DOY to plot (for diruanl average)')
parser.add_option("--dend", dest="dend", default=365, \
                  help = 'final model DOY to plot (for diurnal average)')
parser.add_option("--seasonal", dest="myseasonal", default=False, \
                  action="store_true", help = 'plot seasonal cycle')
parser.add_option("--h1", dest="h1", default=False, \
                  action="store_true", help = 'Use h1 history files')
parser.add_option("--h2", dest="h2", default=False, \
                  action="store_true", help = 'Use h2 history files')
parser.add_option("--h3", dest="h3", default=False, \
                  action="store_true", help = 'Use h3 history files')
parser.add_option("--h4", dest="h4", default=False, \
                  action="store_true", help = 'Use h4 history files')
parser.add_option("--index", dest="index", help = 'index (site or pft)', \
                   default=0)
parser.add_option("--spinup", dest="spinup", help = 'plot Ad and final spinup', \
                   default=False, action="store_true")
parser.add_option("--scale_factor", dest="scale_factor", help = 'scale factor', \
                   default=-999)
parser.add_option("--ylog", dest="ylog", help="log scale for Y axis", \
                  action="store_true", default=False)
parser.add_option("--pdf", dest="pdf", help="save plot to pdf", \
                  action="store_true", default=False)
parser.add_option("--noplot", dest="noplot", help="Do not make plots", \
                  action="store_true", default=False)
(options,args) = parser.parse_args()
               

cesmdir=os.path.abspath(options.mycsmdir)                 

if (options.pdf):
    mpl.use('Agg')	
import matplotlib.pyplot as plt

mycases = options.mycase.split(',')
mysites = options.site.split(',')
mycompsets = options.compset.split(',')

ncases = 1
if (len(mycases) > 1):
  ncases = len(mycases)
  mysites=[]
  for c in range(0,ncases):
    mysites.append(options.site)
    #if (len(mycompsets) == 1):
    mycompsets.append(options.compset)
  mytitles = mycases
elif (len(mysites) > 1):
  ncases = len(mysites)
  mycases=[]
  mycompsets=[]
  for c in range(0,ncases):
    mycases.append(options.mycase)
    mycompsets.append(options.compset)
  mytitles = mysites
elif (len(mycompsets) > 1):
  ncases = len(mycompsets)
  mycases=[]
  mysites=[]
  for c in range(0,ncases):
    mycases.append(options.mycase)
    mysites.append(options.site)  
  mytitles = mycompsets
else:
  mytitles=[]
  mytitles.append(mysites[0])

if (options.titles != ''):
  mytitles = options.titles.split(',')

obs     = options.myobs
myobsdir = '/lustre/or-hydra/cades-ccsi/scratch/dmricciuto/fluxnet'

#get list of variables from varfile
myvars=[]

if (options.myvar == ''):
    if os.path.isfile('./'+options.myvarfile):
        input = open('./'+options.myvarfile)
        for s in input:
            myvars.append(s.strip())
    else:
        print('Error:  invalid varfile')
        sys.exit()
    terminal = 'postscript'
else:
    terminal=''
    myvars = options.myvar.split(',')
    

avpd      = int(options.myavpd)        # desired averaging period in output timestep
ystart    = int(options.myystart)      # beginning year to plot/average
yend      = int(options.myyend)        # final year to plot/average 
yend_all  = yend                       # keep track of last year for which datasets exist
mylat_vals =[]
mylon_vals= []

avtype = 'default'
if (options.mydiurnal):
    avtype = 'diurnal'
    avpd=1
if (options.myseasonal):
    avtype = 'seasonal'

#------------------------------------------------------------------------------

#site = options.site
#compset = options.compset

#dirs=[]
nvar = len(myvars)    
x_toplot    = numpy.zeros([ncases, 2000000], numpy.float)
data_toplot = numpy.zeros([ncases, nvar, 2000000], numpy.float)
obs_toplot  = numpy.zeros([ncases, nvar, 2000000], numpy.float)+numpy.NaN
err_toplot  = numpy.zeros([ncases, nvar, 2000000], numpy.float)+numpy.NaN
snum        = numpy.zeros([ncases], numpy.int)

for c in range(0,ncases):
    if (mycases[c] == ''):
        mydir = cesmdir+'/'+mysites[c]+'_'+mycompsets[c]+'/run/'
    else:
        mydir = cesmdir+'/'+mycases[c]+'_'+mysites[c]+'_'+mycompsets[c]+'/run/'
    print('Processing '+mydir)

    #query lnd_in file for output file information
    if ((options.myhist_mfilt == -999 or options.myhist_nhtfrq == -999)):
        #print('Obtaining output resolution information from lnd_in')
        input = open(mydir+"/lnd_in")
        npf=-999
        tstep=-999
        input = open(mydir+"/lnd_in")
        for s in input:
	    if ('hist_mfilt' in s):
	        mfiltinfo = s.split()[2]
                npf = int(mfiltinfo.split(',')[0])
                if (options.h1): 
                    npf = int(mfiltinfo.split(',')[1])
                if (options.h2):
                    npf = int(mfiltinfo.split(',')[2])
                if (options.h3):
                    npf = int(mfiltinfo.split(',')[3])
                if (options.h4):
                    npf = int(mfiltinfo.split(',')[4])
            if ('hist_nhtfrq' in s):
                nhtfrqinfo = s.split()[2]
                tstep = int(nhtfrqinfo.split(',')[0])
                if (options.h1):
                    tstep = int(nhtfrqinfo.split(',')[1])
                if (options.h2):
                    tstep = int(nhtfrqinfo.split(',')[2])
                if (options.h3):
                    tstep = int(nhtfrqinfo.split(',')[3])
                if (options.h4):
                    tstep = int(nhtfrqinfo.split(',')[4])
        input.close()
    else:
        npf   = int(options.myhist_mfilt)
        tstep = int(options.myhist_nhtfrq)
   
    if (npf == -999 or tstep == -999):
        print('Unable to obtain output file information from lnd_in.  Exiting')
        sys.exit()

    yststr=str(100000+ystart)
    #determine type of file to plot
    if (options.h4):
        hst = 'h4'
    elif (options.h3):
        hst = 'h3'
    elif (options.h2):
        hst = 'h2'
    elif (options.h1):
        hst = 'h1'
    else:
        hst = 'h0'

    if (tstep == 0):
        ftype = 'default'
        mytstep = 'monthly'
        npy=12
    else:
        ftype = 'custom'
        if (abs(npf) == 8760):
            mytstep = 'halfhourly'
            npy=8760
        elif (abs(npf) == 365):
            mytstep = 'daily'
            npy=365
        elif (abs(npf) == 1):
            mytstep = 'annual'
            npy=1
        nhtot=-1*tstep*npf
        nypf = max(1, nhtot/8760)

    #initialize data arrays
    mydata      = numpy.zeros([nvar,2000000], numpy.float)
    myobs       = numpy.zeros([nvar,2000000], numpy.float)+numpy.NaN
    myerr       = numpy.zeros([nvar,2000000], numpy.float)+numpy.NaN
    x           = numpy.zeros([2000000], numpy.float)
    nsteps=0
 
    if (c == 0):   
        var_units=[]
        var_long_names=[]
        myscalefactors=[]

    #Get observations
    if (obs):
        myobsfiles = os.listdir(myobsdir+'/'+mytstep+'/')
        for f in myobsfiles:
            if mysites[c] in f and '.csv' in f:
                myobsfile = myobsdir+'/'+mytstep+'/'+f
        avpd_obs = 1
        if (mytstep == 'halfhourly' and '_HH_' in myobsfile):
            avpd_obs = 2 
        if (os.path.exists(myobsfile) and ystart < 1900):
            print ('Getting start and end year information from observation file')
            thisrow=0
            myobs_input = open(myobsfile)
            for s in myobs_input:
                if thisrow == 1:
                    ystart = int(s[0:4])
                elif (thisrow > 1):
                    if (int(options.myyend) > 9000):
                        yend = min(int(s[0:4]), int(options.myyend))
                thisrow=thisrow+1
            myobs_input.close
        for v in range(0,nvar):
            if (os.path.exists(myobsfile)):
                myobs_in = open(myobsfile)
                thisrow=0
                thisob=0
                for s in myobs_in:
                    if (thisrow == 0):
                        header = s.split(',')
                    else:
                        myvals = s.split(',')
                        thiscol=0
                        if int(myvals[0][0:4]) >= ystart and int(myvals[0][0:4]) <= yend:
                            if (thisob == 0):
                                thisob = (int(myvals[0][0:4])-ystart)*npy*avpd_obs
                            if (thisob % avpd_obs ==  0):
                                myobs[v,thisob/avpd_obs] = 0.0
                                myerr[v,thisob/avpd_obs] = 0.0
                            for h in header:
                                if (h.strip() == 'NEE_CUT_REF' and 'NEE' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] + float(myvals[thiscol])/avpd_obs
                                elif (h.strip () == 'NEE_CUT_REF_JOINTUNC' and \
                                      'NEE' in myvars[v]):
                                    myerr[v,thisob/avpd_obs] = myerr[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'GPP_NT_CUT_REF' and 'GPP' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'GPP_NT_CUT_SE' and 'GPP' in myvars[v]):
                                    myerr[v,thisob/avpd_obs] = myerr[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'RECO_NT_CUT_REF' and 'ER' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'RECO_NET_CUT_SE' and 'ER' in myvars[v]):
                                    myerr[v,thisob/avpd_obs] = myerr[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'LE_F_MDS' and 'EFLX_LH_TOT' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'LE_RANDUNC' and 'EFLX_LH_TOT' in myvars[v]):
                                    myerr[v,thisob/avpd_obs] = myerr[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'H_F_MDS' and 'FSH' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'H_RANDUNC' and 'FSH' in myvars[v]):
                                    myerr[v,thisob/avpd_obs] = myerr[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'TA_F_MDS' and 'TBOT' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'SWIN_F_MDS' and 'FSDS' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'WS_F' and 'WIND' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                elif (h.strip() == 'P_F' and 'RAIN' in myvars[v]):
                                    myobs[v,thisob/avpd_obs] = myobs[v,thisob/avpd_obs] +float(myvals[thiscol])/avpd_obs
                                if myobs[v,thisob/avpd_obs] < -4000:
                                    myobs[v,thisob/avpd_obs] = numpy.NaN
                                thiscol=thiscol+1
                            thisob=thisob+1
                    thisrow = thisrow+1
                myobs_in.close()
            else:
                print 'Error reading observations for '+mysites[c]

    #read monthly .nc files (default output)
    if (ftype == 'default'):
        for v in range(0,nvar):
            nsteps = 0
            for y in range(ystart,yend+1):
                yst=str(10000+y)[1:5]
                for m in range(0,12):
                    mst=str(101+m)[1:3]
                    myfile = os.path.abspath(mydir+'/'+mycases[c]+'_'+mysites[c]+'_'+mycompsets[c]+ \
                                             ".clm2."+hst+"."+yst+"-"+mst+".nc")
                    #get units/long names from first file
                    if (os.path.exists(myfile)):
                        if (y == ystart and m == 0 and c == 0):
                            nffile = netcdf.netcdf_file(myfile,"r")
                            varout=nffile.variables[myvars[v]]
                            var_long_names.append(varout.long_name)
                            nffile.close()
                            if (float(options.scale_factor) < -900):
                                if ('gC/m^2/s' in varout.units):
                                    myscalefactors.append(3600*24)
                                    var_units.append('g.C/m2/day')
                                else:
                                    myscalefactors.append(1.0)
                                    var_units.append(varout.units.replace('^',''))
                            else:
                                myscalefactors.append(float(options.scale_factor))
                                var_units.append(varout.units.replace('^',''))
                        
                        if (y == ystart and m == 0 and v == 0):      # get lat/lon info
                            nffile = netcdf.netcdf_file(myfile,"r")
                            mylat_vals.append(nffile.variables['lat'][0])
                            mylon_vals.append(nffile.variables['lon'][0])
                            nffile.close()
                        x[nsteps] = y+m/12.0
                        myvar_temp = getvar(myfile, myvars[v],npf,int(options.index), \
                                            myscalefactors[v])
                        mydata[v,nsteps] = myvar_temp
                        if (myvars[v] == 'RAIN'):
                            myvar_temp2 = getvar(myfile,'SNOW',npf,int(options.index), \
                                            myscalefactors[v])
                            mydata[v,nsteps] = mydata[v,nsteps]+myvar_temp2
                    else:
                        if (v == 0 and m == 0):
                           print 'Warning: '+myfile+' does not exist'
                        x[nsteps] = y+m/12.0
                        mydata[v,nsteps] = numpy.NaN
                        if (y-1 < yend_all):
                            yend_all = y-1
                    nsteps = nsteps+1

    #read annual .nc files
    if (ftype == 'custom'):
        for v in range(0,nvar):
            nsteps=0
            nfiles = (yend-ystart)/nypf
            nc=1
            starti=0
            ylast=0
            if (options.spinup):
                nc=2
            if (npf == 1):
                starti = 1
            for n in range(0,nc):
                if ((options.spinup)and n== 0):
                    if (mycases[c] == ''):
                        if (options.ad_Pinit):
                            mydir = cesmdir+'/'+mysites[c]+'_'+mycompsets[c]+'_ad_spinup/run/'
                        else:
                            mydir = cesmdir+'/'+mysites[c]+'_'+mycompsets[c].replace('CNP','CN')+ \
	                          '_ad_spinup/run/'
                    else:
                        if (options.ad_Pinit):
                            mydir = cesmdir+'/'+mycases[c]+'_'+mysites[c]+'_'+ \
                                    mycompsets[c]+'_ad_spinup/run/'
                            thiscompset = mycompsets[c]+'_ad_spinup'
                        else:
                            mydir = cesmdir+'/'+mycases[c]+'_'+mysites[c]+'_'+ \
                                    mycompsets[c].replace('CNP','CN')+'_ad_spinup/run/'
                            thiscompset = mycompsets[c].replace('CNP','CN')+'_ad_spinup'
                else:
                    if (mycases[c] == ''):
                        mydir = cesmdir+'/'+mysites[c]+'_'+mycompsets[c]+'/run/'
                    else:
                        mydir = cesmdir+'/'+mycases[c]+'_'+mysites[c]+'_'+ \
                                mycompsets[c]+'/run/'
                    thiscompset = mycompsets[c]
                for y in range(starti,nfiles+1):   #skip first file in spinup
                    yst=str(10000+ystart+(y*nypf))[1:5]
                    if (mycases[c].strip() == ''):
                        myfile = os.path.abspath(mydir+'/'+mycases[c]+'_'+thiscompset+".clm2."+hst+ \
                                                 "."+yst+"-01-01-00000.nc")
                    else:
                        myfile = os.path.abspath(mydir+'/'+mycases[c]+"_"+mysites[c]+'_'+thiscompset+ \
                                                 ".clm2."+hst+"."+yst+"-01-01-00000.nc")
                    if (os.path.exists(myfile)):
                        if (n == 0):
                            ylast = y
                        if (y == starti and n == 0 and c == 0):
                            nffile = netcdf.netcdf_file(myfile,"r")
                            varout=nffile.variables[myvars[v]]
                            var_long_names.append(varout.long_name)
                            if (float(options.scale_factor) < -900):
                                if ('gC/m^2/s' in varout.units):
                                    if (npf >= 365):
                                        myscalefactors.append(3600*24)
                                        var_units.append('g.C/m2/day')
                                    else:
                                        myscalefactors.append(3600*24*365)
                                        var_units.append('g.C/m2/yr')
                                else:
                                    myscalefactors.append(1.0)
                                    var_units.append(varout.units.replace('^',''))
                            else:
                                 myscalefactors.append(float(options.scale_factor))
                                 var_units.append(varout.units.replace('^',''))
                            nffile.close()
                        if (y == starti and n == 0 and v == 0):      # get lat/lon info
                            nffile = netcdf.netcdf_file(myfile,"r")
                            mylat_vals.append(nffile.variables['lat'][0])
                            mylon_vals.append(nffile.variables['lon'][0])
                            nffile.close()

                        myvar_temp = getvar(myfile,myvars[v],npf,int(options.index), \
                                            myscalefactors[v])
                        if (myvars[v] == 'RAIN'):
                            myvar_temp2 = getvar(myfile,'SNOW',npf,int(options.index), \
                                            myscalefactors[v])

                        if len(myvar_temp) == npf:
                            for i in range(0,npf):
                                myind = ylast*n*npf+y*npf+i
                                x[nsteps] = ystart+(ylast*n*nypf+y*nypf) + nypf*(i*1.0-0.5)/npf
                                mydata[v,nsteps] = myvar_temp[i]
                                if (myvars[v] == 'RAIN'):    #add snow for total precip
                                    mydata[v,nsteps] = mydata[v,nsteps]+myvar_temp2[i]
                                nsteps=nsteps+1
                        else:
                            for i in range(0,npf):
                                myind=ylast*n*npf+(y-1)*npf+i
                                x[myind] = ystart+(ylast*n*nypf+y*nypf) + nypf*(i*1.0-0.5)/npf
                                mydata[v,myind] = numpy.NaN
                                nsteps=nsteps+1 
  		    else:	
                         if (v == 0):
                             print('Warning: '+myfile+' does not exist')
                         if (y-1 < yend_all):
                             yend_all = y-1
                         for i in range(0,npf):
                             if (n == nc-1):
                                 myind=ylast*n*npf+y*npf+i
                                 x[myind] = ystart+(ylast*n*nypf+y*nypf) + nypf*(i*1.0-0.5)/npf
                                 mydata[v,myind] = numpy.NaN
                                 nsteps=nsteps+1

    #perform averaging and write output files 
    if (avtype == 'default'):
        for v in range(0,nvar):
            snum[c] = 0
            for s in range(0,int(nsteps/avpd)): 
                x_toplot[c, snum[c]]       = sum(x[s*avpd:(s+1)*avpd])/avpd
                data_toplot[c, v, snum[c]] = sum(mydata[v,s*avpd:(s+1)*avpd])/avpd
                obs_toplot[c, v, snum[c]]  = sum(myobs[v,s*avpd:(s+1)*avpd])/avpd
                if (min(myerr[v,s*avpd:(s+1)*avpd]) < -9000):
                    err_toplot[c,v,snum[c]] = 0
                else:
                    err_toplot[c, v, snum[c]]  = sum(myerr[v,s*avpd:(s+1)*avpd])/avpd
                snum[c] = snum[c]+1
          

    #diurnal average (must have hourly output)
    if (avtype == 'diurnal'):
        snum[c]=24
        for v in range(0,nvar):
            mysum = numpy.zeros(snum[c], numpy.float)
            mysum_obs = numpy.zeros(snum[c], numpy.float)
            myct = numpy.zeros(snum[c],numpy.float)
            myct_obs = numpy.zeros(snum[c],numpy.float)
            for y in range(0,(yend_all-ystart+1)):
                for d in range (int(options.dstart),int(options.dend)):
                    for s in range(0,snum[c]):        
                        h=s
                        if (h >= 24):
                            h=h-24
                        mysum[s] = mysum[s]+mydata[v,y*8760+(d-1)*24+h-int(options.timezone)+1]
                        myct[s] = myct[s]+1
                        if (myobs[v,y*8760+(d-1)*24+h] > -900):
                            mysum_obs[s] = mysum_obs[s]+myobs[v,y*8760+(d-1)*24+h]
                            myct_obs[s] = myct_obs[s]+1
            for s in range(0,snum[c]):
                if (myct_obs[s] > 0):
                    mysum_obs[s] = mysum_obs[s]/myct_obs[s]
                else:
                    mysum_obs[s] = numpy.NaN
                x_toplot[c,s] = s+1
                obs_toplot[c, v, s] = mysum_obs[s]
                data_toplot[c, v, s] = mysum[s]/myct[s]

    #seasonal average (assumes default monthly output)
    if (avtype == 'seasonal'):
        for v in range(0,nvar):
            snum[c] = 12
            mysum=numpy.zeros(snum[c], numpy.float)
            mysum_obs = numpy.zeros(snum[c], numpy.float)
            mycount_obs = numpy.zeros(snum[c], numpy.int)
            for y in range(0,(yend_all-ystart+1)):
                for s in range(0,snum[c]):
                    mysum[s]=mysum[s]+mydata[v,(y*12+s)]/float(yend_all-ystart+1)
                    if (myobs[v,(y*12+s)] > -900):
                        mysum_obs[s]=mysum_obs[s]+myobs[v,(y*12+s)]
                        mycount_obs[s] = mycount_obs[s]+1
            for s in range(0,snum[c]):
                if (mycount_obs[s] > 0):
                    mysum_obs[s] = mysum_obs[s]/mycount_obs[s]
                else:
                    mysum_obs[s] = numpy.NaN
                x_toplot[c,s] = s+1.5
                obs_toplot[c,v,s]  = mysum_obs[s]
                data_toplot[c,v,s] = mysum[s]

                
#diagnostics, outputs and plots

if (options.spinup):
    analysis_type = 'spinup'
elif (options.mydiurnal):
    analysis_type = 'diurnalcycle_'+str(options.dstart)+'_'+str(options.dend)
elif (options.myseasonal):
    analysis_type = 'seasonalcycle'
elif (mytstep == 'halfhourly'):
    analysis_type = 'hourly'
else:
    analysis_type=mytstep

rmse = numpy.zeros([len(myvars),ncases],numpy.float)
bias = numpy.zeros([len(myvars),ncases],numpy.float)
corr = numpy.zeros([len(myvars),ncases],numpy.float)

for v in range(0,len(myvars)):
    if (not options.noplot):
        fig = plt.figure()
        ax = plt.subplot(111)
    colors=['b','g','r','c','m','y','k','b','g','r','c','m','y','k','b','g','r','c','m','y','k']
    styles=['-','-','-','-','-','-','-','--','--','--','--','--','--','--','-.','-.','-.','-.','-.','-.','-.']
    for c in range(0,ncases):
        #Output data in netcdf format
        if (c == 0):
            if (v == 0):
                ftype_suffix=['model','obs']
                os.system('mkdir -p ./plots/'+mycases[0]+'/'+analysis_type)
                for ftype in range(0,2):
                    outdata = netcdf.netcdf_file('./plots/'+mycases[0]+'/'+analysis_type+'/'+mycases[0]+"_"+mysites[0]+'_'+ftype_suffix[ftype]+".nc","w",mmap=False)
                    outdata.createDimension('time',snum[c])
                    #outdata.createDimension('lat',ncases)
                    #outdata.createDimension('lon',ncases)
                    outdata.createDimension('gridcell',ncases)
                    outdata.createDimension('strlen',6)
                    mylat = outdata.createVariable('lat','f',('gridcell',))
                    mylat.long_name='coordinate latitude'
                    mylat.units='degrees_north'
                    mylon = outdata.createVariable('lon','f',('gridcell',))
                    mylon.long_name='coordinate longitude'
                    mylon.units='degrees_east'
                    mytime = outdata.createVariable('time','f',('time',))
                    mytime.long_name='time'
                    mytime.units='days since '+str(ystart)+'-01-01 00:00:00'
                    mytime.calendar='noleap'
                    mytime[:] = (x_toplot[0,0:snum[c]]-ystart)*365
                    #myname = outdata.createVariable('site_name','c',('lat','lon','strlen'))
                    myname = outdata.createVariable('site_name','c',('gridcell','strlen'))
                    myname[:,:] = ''   #changed for gridcell
                    outdata.close()
        for ftype in range(0,2):
            outdata = netcdf.netcdf_file('./plots/'+mycases[0]+'/'+analysis_type+'/'+mycases[0]+"_"+mysites[0]+'_'+ftype_suffix[ftype]+".nc","a",mmap=False)
            if (c == 0):
                #myvar = outdata.createVariable(myvars[v],'f',('time','lat','lon'))
                myvar = outdata.createVariable(myvars[v],'f',('time','gridcell'))
                myvar.units=var_units[v]
                myvar.missing_value=1e36
                myvar[:,:]=myvar.missing_value   #changed for gridcell
            else:
                myvar=outdata.variables[myvars[v]]
            scalefac = 1.0
            if (var_units[v] == 'g.C/m2/day'):
                myvar.units = 'kg.C/m2/s'
                scalefac = 1.0 / (3600*24*1000.0)
            if (ftype == 0):
                myvar[:,c] = data_toplot[c,v,0:snum[c]]*scalefac  #changed for gridcell
            if (ftype == 1):
                myvar[:,c] = obs_toplot[c,v,0:snum[c]]*scalefac   #changed for gridcell
            if (v == 0):
                myname = outdata.variables['site_name']
                myname[c,0:6] = str(mysites[c])[0:6]   #changed for gridcell
                mylat = outdata.variables['lat']
                mylat[c] = mylat_vals[c]
                mylon = outdata.variables['lon']
                mylon[c] = mylon_vals[c]
            outdata.close()
        
        #----------------------- plotting ---------------------------------
        if (options.noplot == False):
            gind=[]
            for i in range(0,snum[c]):
                if (obs_toplot[c,v,i] < -900):
                    obs_toplot[c,v,i] = numpy.nan
                else:
                    gind.append(i)
                    rmse[v,c] = rmse[v,c] + (data_toplot[c,v,i]-obs_toplot[c,v,i])**2.0
                    bias[v,c] = bias[v,c] + (data_toplot[c,v,i]-obs_toplot[c,v,i])
            rmse[v,c] = (rmse[v,c]/len(gind))**0.5
            bias[v,c] = bias[v,c]/len(gind)
            corr[v,c] = numpy.corrcoef(data_toplot[c,v,gind],obs_toplot[c,v,gind])[1,0]

            if (options.ylog):
                ax.plot(x_toplot[c, 0:snum[c]], abs(data_toplot[c,v,0:snum[c]]), label=mytitles[c], color=colors[c], \
                  linestyle=styles[c], linewidth=3)
            else:
                ax.plot(x_toplot[c, 0:snum[c]], (data_toplot[c,v,0:snum[c]]), label=mytitles[c], color=colors[c], \
	          linestyle=styles[c], linewidth=3)
                if (options.myobs and c == 0):
                    ax.errorbar(x_toplot[c, 0:snum[c]], obs_toplot[c,v,0:snum[c]], yerr = err_toplot[c,v,0:snum[c]], \
                                color=colors[c], fmt='o')

    if (options.noplot == False):
        if (avtype == 'seasonal'):
            plt.xlabel('Model Month')
        elif (avtype == 'diurnal'):
            plt.xlabel('Model Hour (LST)')
        else:
            plt.xlabel('Model Year')
 
        plt.ylabel(myvars[v]+' ('+var_units[v]+')')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 10})
        plt.title(var_long_names[v]+' at '+mysites[0])
        if (options.ylog):
            plt.yscale('log')

        os.system('mkdir -p ./plots/'+mycases[0]+'/'+analysis_type)
        if (options.pdf):
            fig_filename = './plots/'+mycases[0]+'/'+analysis_type+'/'+mysites[0]+'_'+myvars[v]
            fig_filename = fig_filename+'_'+analysis_type
            fig.savefig(fig_filename+'.pdf')

if (not options.pdf and not options.noplot):
    plt.show()
