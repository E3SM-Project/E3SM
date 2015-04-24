#!/usr/bin/python

import os, sys, csv, glob
import Scientific.IO.NetCDF
from Scientific.IO import NetCDF
import numpy
import pp
from optparse import OptionParser

def getvar(fname, var, npf):
    nffile = Scientific.IO.NetCDF.NetCDFFile(fname,"r")
    varout=nffile.variables[var]
    thisval = varout.getValue()[0:npf,0]
    nffile.close()
    return thisval


parser = OptionParser()
parser.add_option("--csmdir", dest="mycsmdir", default='../../..', \
                  help = 'Base CESM directory (default = ../../../)')
parser.add_option("--case", dest="mycase", default='', \
                  help = "name of case id prefix to plot")
parser.add_option("--case2", dest="mycase2", default='', \
                  help = "name of second case to plot")
parser.add_option("--case3", dest="mycase3", default='', \
                  help = "name of third case to plot")
parser.add_option("--case4", dest="mycase4", default='', \
                  help = "name of second case to plot")
parser.add_option("--case5", dest="mycase5", default='', \
                  help = "name of third case to plot")
#parser.add_option("--compset", dest="compset", default="I20TRCLM45CN", \
#                  help = "Compset to plot")
parser.add_option("--title", dest="mytitle", default='', \
                  help = "title of case to plot (for legend)")
parser.add_option("--title2", dest="mytitle2", default='', \
                  help = "title of second case to plot (for legend)")
parser.add_option("--title3", dest="mytitle3", default='', \
                  help = "title of third case to plot (for legend)")
parser.add_option("--title4", dest="mytitle4", default='', \
                  help = "title of second case to plot (for legend)")
parser.add_option("--title5", dest="mytitle5", default='', \
                  help = "title of third case to plot (for legend)")
parser.add_option("--obs", action="store_true", default=False, \
                  help = "plot observations", dest="myobs")
#parser.add_option("--site", dest="site", default="none", \
#                  help = 'site (to plot observations)')
parser.add_option("--varfile", dest="myvarfile", default='varfile', \
                  help = 'file containing list of variables to plot')
parser.add_option("--var", dest="myvar", default='', \
                  help="variable to plot (overrides varfile, " \
                  +"sends plot to screen")
parser.add_option("--avpd", dest="myavpd", default=1, \
                  help = 'averaging period in # of output timesteps' \
                  +' (default = 1)')
parser.add_option("--hist_mfilt", dest="myhist_mfilt", default=1, \
                  help = 'beginning model year to plot')
parser.add_option("--hist_nhtfrq", dest="myhist_nhtfrq", default=0, \
                  help = 'beginning model year to plot')
parser.add_option("--ystart", dest="myystart", default=1, \
                  help = 'beginning model year to plot')
parser.add_option("--yend", dest="myyend", default=1, \
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

(options,args) = parser.parse_args()
                  
#set up parallel
ppservers = ("*",)
#ppservers=()
job_server = pp.Server(ppservers=ppservers, secret="mypswd")
print "Starting pp with", job_server.get_ncpus(), "workers"

cesmdir=os.path.abspath(options.mycsmdir)                 
if (options.mycase == ''): # or os.path.exists(options.mycase) == False):
    print('Error: invalid CESM root directory')
    sys.exit()
else:
    mycase1 = options.mycase
    if (options.mytitle == ''):
        mytitle1=mycase1
    else:
        mytitle1=options.mytitle

mycase2 = options.mycase2
if (options.mytitle2 == ''):
    mytitle2=mycase2
else:
    mytitle2=options.mytitle2
mycase3 = options.mycase3
if (options.mytitle3 == ''):
    mytitle3=mycase3
else:
    mytitle3=options.mytitle3
mycase4 = options.mycase4
if (options.mytitle4 == ''):
    mytitle4=mycase4
else:
    mytitle4=options.mytitle4
mycase5 = options.mycase5
if (options.mytitle5 == ''):
    mytitle5=mycase5
else:
    mytitle5=options.mytitle5

obs     = options.myobs

#get list of variables from varfile
myvars=[]

if (options.myvar == ''):
    if os.path.isfile('./PTCLM_files/'+options.myvarfile):
        input = open('./PTCLM_files/'+options.myvarfile)
        for s in input:
            myvars.append(s.strip())
    else:
        print('Error:  invalid varfile')
        sys.exit()
    terminal = 'postscript'
else:
    terminal=''
    myvars.append(options.myvar)
    

avpd      = int(options.myavpd)        # desired averaging period in output timestep
ystart    = int(options.myystart)      # beginning year to plot/average
yend      = int(options.myyend)        # final year to plot/average 

avtype = 'default'
if (options.mydiurnal):
    avtype = 'diurnal'
    avpd=1
if (options.myseasonal):
    avtype = 'seasonal'

#------------------------------------------------------------------------------

#get # of cases to plot
ncases=1
#ncases=2    #TEMP
if mycase2 != '':
    ncases=2
    if mycase3 !='':
        ncases=3
        if mycase4 != '':
            ncases=4
            if mycase5 != '':
                ncases=5
if (obs):
    ncases=ncases+1

#site = options.site
#compset = options.compset
for c in range(0,ncases):
    if (obs and c == ncases-1):   #observations
        diro='/home/zdr/models/LoTEC/assim/observations'
        os.chdir(diro)
    elif (c == 0):
        mycase = mycase1
        dir1=cesmdir+'/run/'+mycase1+'/run/' #+'_'+site+'_'+compset+'/run/'
        os.chdir(dir1)
    elif (c == 1):
        mycase = mycase2
        dir2=cesmdir+'/run/'+mycase2+'/run/' #+'_'+site+'_'+compset+'/run/'
        #dir2=cesmdir+'/run/'+site+'_'+compset+'/run'  #TEMP
        os.chdir(dir2)
    elif (c == 2):
        mycase = mycase3
        dir3=cesmdir+'/run/'+mycase3+'/run/' #+'_'+site+'_'+compset+'/run/'
        os.chdir(dir3)
    elif (c == 3):
        mycase = mycase4
        dir4=cesmdir+'/run/'+mycase4+'/run/' #+'_'+site+'_'+compset+'/run/'
        os.chdir(dir4)
    elif (c == 4):
        mycase = mycase5
        dir5=cesmdir+'/run/'+mycase5+'/run/' #+'_'+site+'_'+compset+'/run/'
        os.chdir(dir5)
 
    #query lnd_in file for output file information
    os.system('pwd')
    if ((options.myhist_mfilt == -999 or options.myhist_nhtfrq == -999) or (obs and  c  < ncases-1)):
        print('Obtaining output resolution information from lnd_in')
        input = open("lnd_in")
        hash1 = False
        hash2 = False
        for s in input:
            s2=s.strip()
            if (s2.startswith('hist_fincl2')):
                hash1 = True
            if (s2.startswith('hist_fincl3')):
                hash2 = True
        input.close()
        npf=-999
        tstep=-999
        input = open("lnd_in")
        print(hash1, hash2)
        for s in input:
            if (s[0:11].strip() == 'hist_mfilt'):
                if (hash2):
                    s1,s2,val,val1,val2=s.split()
                elif (hash1):
                    s1,s2,val,val1 = s.split()
                else:
                    s1,s2,val=s.split() 
                
                if (options.h2):
                    if (hash2):
                        npf = int(val2)
                    else:
                        print('h2 files requested but do not exist')
                elif (options.h1):
                    if (hash1 and hash2):
                        npf = int(val1[:-1])
                    elif (hash1):
                        npf = int(val1)
                    else:
                        print('h1 files requested but do not exist')
                else:
                    if (hash1 or hash2):
                        npf = int(val[:-1])
                    else:
                        npf = int(val)
                print(npf)
                
            if (s[0:12].strip() == 'hist_nhtfrq'):
                if (hash2):
                    s1,s2,val,val1,val2=s.split()
                    print(val, val1, val2)
                elif (hash1):
                    s1,s2,val,val1 = s.split()
                else:
                    s1,s2,val=s.split() 
                
                if (options.h2):
                    if (hash2):
                        tstep= int(val2)
                elif (options.h1):
                    if (hash1 and hash2):
                        tstep = int(val1[:-1])
                    elif (hash1):
                        tstep = int(val1)
                else:
                    if (hash1 or hash2):
                        tstep = int(val[:-1])
                    else:
                        tstep = int(val)
                print(tstep)
        input.close()
    elif (c == ncases-1 and obs):
        npf   = 8760
        tstep = -1
    else:
        npf   = int(options.myhist_mfilt)
        tstep = int(options.myhist_nhtfrq)
   
    print npf, tstep 
    if (npf == -999 or tstep == -999):
        print('Unable to obtain output file information from lnd_in.  Exiting')
        sys.exit()

    yststr=str(100000+ystart)
    #determine type of file to plot
    if (tstep == 0):
        ftype = 'default'
        if (options.h2):
            hst = 'h2'
        elif (options.h1):
            hst = 'h1'
        else:
            hst = 'h0'
        testfile = mycase #+'_'+site+'_'+compset+'.clm2.'+hst+'.'+yststr[2:6]+'-01.nc'
    else:
        ftype = 'custom'
        nhtot=-1*tstep*npf
        if (nhtot != 8760):
            print('Only annual or default (monthly) files are supported')
            sys.exit()
        if (options.h2):
            hst='h2'
        elif (options.h1):
            hst='h1'
        else:
            hst='h0'
        #testfile = mycase+'_'+site+'_'+compset+'.clm2.'+hst+'.'+yststr[2:6]+'-01-01-00000.nc'
        testfile = mycase+'.clm2.'+hst+'.'+yststr[2:6]+'-01-01-00000.nc'
        

    if (obs and c == ncases-1):
        ftype='obs'
        if (options.site == 'none'):
            print 'Error:  No site specified.  Cannot load observations'
        testfile = options.site+'obs.nc'

    #check for output files (if not here, change to archive directory)
    print(testfile)
    os.system('pwd')
    if (os.path.isfile(testfile) == False):
        print('Output not in run directory.  Switching to archive directory')
        archdir=cesmdir+'/run/archive/'+mycase+'/lnd/hist' #+'_'+site+'_'+compset+'/lnd/hist'
        print(archdir)
        if (os.path.exists(archdir) == False):
            print('Archive directory does not exist.  Exiting')
            sys.exit()
        else:
            os.chdir(archdir)
            if (c == 0):
                dir1=archdir
            if (c == 1):
                dir2=archdir
            if (c == 2):
                dir3=archdir
            #if (ftype == 'obs'):
            #    diro=archdir
            if (os.path.isfile(testfile) == False):
                print('Output not found.  Exiting')
                sys.exit()
                

    #initialize data arrays
    nvar=len(myvars)
    mydata=numpy.zeros([nvar,2000000], numpy.float)
    x=numpy.zeros([2000000], numpy.float)
    nsteps=0
 
    if (c == 0):   
        var_units=[]
        var_long_names=[]

    #read monthly .nc files (default output)
    if (ftype == 'default'):
        jobs=[]
        for y in range(ystart,yend+1):
            yst=str(10000+y)[1:5]
            for m in range(0,12):
                mst=str(101+m)[1:3]
                #myfile = os.path.abspath('./'+mycase+'_'+site+'_'+compset+".clm2."+hst+"."+yst+"-"+mst+".nc")
                myfile = os.path.abspath('./'+mycase+".clm2."+hst+"."+yst+"-"+mst+".nc")
                nstepslast=nsteps
                for v in range(0,nvar):
                    #get units/long names from first file
                    if (y == ystart and m == 0 and c == 0):
                        nffile = Scientific.IO.NetCDF.NetCDFFile(myfile,"r")
                        varout=nffile.variables[myvars[v]]
                        var_long_names.append(varout.long_name)
                        var_units.append(varout.units)
                        nffile.close()
                    nsteps=nstepslast
                    x[nsteps] = y+m/12.0
                    nsteps=nsteps+1
                    jobs.append(job_server.submit(getvar,(myfile,myvars[v],npf),(), \
                                    ("Scientific.IO.NetCDF",)))
                    nsteps = 0
                    for job in jobs:
                        mydata[v,nsteps] = job()
            #print(v, nsteps, mydata[v,nsteps])
#                    if (myvars[v] == 'NEE' or myvars[v] == 'GPP'):
#                        mydata[v,nsteps]=mydata[v,nsteps]*1e6/12.0
                        nsteps = nsteps + 1


    #read annual .nc files
    if (ftype == 'custom'):
        for v in range(0,nvar):
            jobs = []
            for y in range(ystart,yend+1):
                yst=str(10000+y)[1:5]
                #myfile = os.path.abspath('./'+mycase+"_"+site+'_'+compset+".clm2."+hst+"."+yst+\
                #                             "-01-01-00000.nc")
                myfile = os.path.abspath('./'+mycase+".clm2."+hst+"."+yst+\
                                             "-01-01-00000.nc")
                if (y == ystart and c == 0):
                    nffile = Scientific.IO.NetCDF.NetCDFFile(myfile,"r")
                    varout=nffile.variables[myvars[v]]
                    var_long_names.append(varout.long_name)
                    var_units.append(varout.units)
                    nffile.close()
                jobs.append(job_server.submit(getvar,(myfile,myvars[v],npf,),(), \
                                                  ("Scientific.IO.NetCDF",)))
            nsteps = 0
            for job in jobs:
                for i in range(0,npf):
                    mydata[v,nsteps] = job()[i]
                    x[nsteps] = ystart+(nsteps*1.0)/npf
                    nsteps=nsteps+1
        
    #read obervation file, assumes it is in case directory (years must match!)
    #will work for NEE only!
    print ftype
    if (ftype == 'obs'):
        npd=24                #number of time steps per day
        lst=-5 #local standard time (diff from UTC)
        npf=(yend-ystart+1)*365*npd
        nffile = NetCDF.NetCDFFile(options.site+'obs.nc',"r")
        nsteps=-(lst-1)
        nstepslast=nsteps
        for v in range(0,nvar):
            if (myvars[v] == 'NEE' or myvars[v] == 'GPP'):
                nsteps=nstepslast
                varout=nffile.variables[myvars[v]+'_filled']
                print(npf)
                indata = varout.getValue()[0:npf]*12/1e6
                x[0:9]=ystart
                for i in range(0,npf*24/npd):
                    x[nsteps] = ystart+(i*1.0)/8760
                    if (npd == 24):
                        mydata[v,max(nsteps,0)] = float(indata[i])
                    if (npd == 48):
                        mydata[v,max(nsteps,0)] = (float(indata[i*2])+ \
                                                       float(indata[i*2+1]))/2
                    nsteps=nsteps+1
        nffile.close()
        nsteps=nsteps+lst-1
    
    #perform averaging and write output files for gnuplot
    if (avtype == 'default'):
        print(nsteps, avpd)
        for v in range(0,nvar):
            output=open(myvars[v]+'_'+str(c)+".txt","w")
            for s in range(0,int(nsteps/avpd)):        
                output.write(str(sum(x[s*avpd:(s+1)*avpd])/avpd)+ " " + \
                             str(sum(mydata[v,s*avpd:(s+1)*avpd])/avpd)+"\n")
            output.close()

    #diurnal average (must have hourly output)
    if (avtype == 'diurnal'):
        for v in range(0,nvar):
            output=open(myvars[v]+'_'+str(c)+".txt","w")
            mysum = numpy.zeros(36, numpy.float)
            myct = numpy.zeros(36,numpy.float)
            for y in range(0,(yend-ystart+1)):
                for d in range (int(options.dstart),int(options.dend)):
                    for s in range(0,36):        
                        h=s
                        if (h >= 24):
                            h=h-24
                        mysum[s] = mysum[s]+mydata[v,y*8760+(d-1)*24+h]/((yend-ystart+1)* \
                                                 (int(options.dend)-int(options.dstart)+1))
            for s in range(0,36):
                output.write(str(s+0.5)+ " " + \
                                 str(mysum[s])+"\n")
            output.close()
        
    #seasonal average (assumes default monthly output)
    if (avtype == 'seasonal'):
        for v in range(0,nvar):
            output=open(myvars[v]+'_'+str(c)+".txt","w")
            np = 12
            mysum=numpy.zeros(np, numpy.float)
            
            for y in range(0,(yend-ystart+1)):
                for s in range(0,np):
                    print(y, s, v, mydata[v,y*12+s])
                    mysum[s]=mysum[s]+mydata[v,(y*12+s)]/(yend-ystart+1)
        
            for s in range(0,np):
                output.write(str(s+0.5)+" "+str(mysum[s])+"\n")
            output.close()
                
                
#create gnuplot script
output=open('plotmyvar.p','w')
#output.write('set terminal '+terminal+' enhanced color\n')


if terminal != '':
    output.write('set terminal '+terminal+'\n') #+' enhanced color\n')
    #plotfile = mycase1+'_'+site+'_'+compset+'_plots.ps'
    plotfile =  mycase1+'_plots.ps'
    if (avtype == 'seasonal'):
        #plotfile = mycase1+'_'+site+'_'+compset+'_seasonal_plots.ps'
        plotfile = mycase1+'_seasonal_plots.ps'
    elif (avtype == 'diurnal'):
        #plotfile = mycase1+'_'+site+'_'+compset+'_diurnal_plots.ps'
        plotfile = mycase1+'_diurnal_plots.ps'
    output.write('set output "'+plotfile+'"\n')
output.write('set style line 1 lt 1 lw 3 lc rgb "red"\n')
output.write('set style line 2 lt 1 lw 3 lc rgb "blue"\n')
output.write('set style line 3 lt 1 lw 3 lc rgb "black"\n')
output.write('set style line 4 lt 1 lw 3 lc rgb "cyan"\n')
output.write('set style line 5 lt 1 lw 3 lc rgb "green"\n')

for v in range(0,nvar):
    print var_units
    cmd  = 'plot "'+dir1+'/'+myvars[v]+'_0.txt" with lines linestyle 1 title "' \
           +mytitle1+'" ' 
    if (avtype == 'seasonal'):
        cmd_xlab = 'set xlabel "model month" '
    elif (avtype == 'diurnal'):
        cmd_xlab = 'set xlabel "model hour (UTC)" '
    else:
        cmd_xlab = 'set xlabel "model year" '
    cmd_ylab = 'set ylabel "'+myvars[v]+' ('+var_units[v]+')" '
    cmd_titl = 'set title  "'+var_long_names[v] #+' at '+site+'" '
    if ((ncases >= 2 and obs == False) or (ncases >= 3 and obs == True)):
        cmd = cmd+', "'+dir2+'/'+myvars[v]+'_1.txt" with lines linestyle 2 title "' \
              +mytitle2+'" '
    if ((ncases >= 3 and obs == False) or (ncases >= 4 and obs == True)):
        cmd = cmd+', "'+dir3+'/'+myvars[v]+'_2.txt" with lines linestyle 3 title "' \
              +mytitle3+'" '
    if ((ncases >= 4 and obs == False) or (ncases >= 5 and obs == True)):
        cmd = cmd+', "'+dir4+'/'+myvars[v]+'_3.txt" with lines linestyle 4 title "' \
              +mytitle4+'" '
    if ((ncases >= 5 and obs == False) or (ncases >= 6 and obs == True)):
        cmd = cmd+', "'+dir5+'/'+myvars[v]+'_4.txt" with lines linestyle 5 title "' \
              +mytitle5+'" '
    if (obs and myvars[v] == 'NEE' or myvars[v] == 'GPP'):
        cmd = cmd+', "'+diro+'/'+myvars[v]+'_'+str(ncases-1)+'.txt" with lines linestyle '+str(ncases)+ \
            ' title "observed"'
    output.write(cmd_xlab+'\n')
    output.write(cmd_ylab+'\n')
    output.write(cmd_titl+'\n')
    #output.write('set xrange [1997:2002]\n')
    output.write(cmd+'\n')
    #output.write('set xrange [1997:2002]\n')
    output.write('replot\n')
output.close()

#execute gnuplot script
os.system('gnuplot -persist plotmyvar.p')
os.system('mkdir '+cesmdir+'/scripts/plots/'+mycase1) #+'_'+site+'_'+compset)
if (terminal == 'postscript'):
    os.system('ps2pdf '+plotfile+' '+plotfile[:-4]+'.pdf')
    os.system('mv '+plotfile[:-4]+'.pdf '+cesmdir+'/scripts/plots/'+ \
                  mycase1) #+'_'+site+'_'+compset)
    os.system('mv '+plotfile+' '+cesmdir+'/scripts/plots/'+ \
                  mycase1) #+'_'+site+'_'+compset) 
#os.system('ps2pdf '+cesmdir+'/scripts/plots/'+mycase1+'/'+mycase1+'_plots.ps')
