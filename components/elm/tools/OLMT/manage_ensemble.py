#!/usr/bin/env python3
import sys,os, time
import numpy as np
import netcdf4_functions as nffun
import subprocess
from mpi4py import MPI
from optparse import OptionParser

#MPI python code used to manage the ensemble simulations 
#  and perform post-processing of model output.
#  DMRicciuto 7/14/2016

parser = OptionParser()

parser.add_option("--runroot", dest="runroot", default="../../run", \
                  help="Directory where the run would be created")
parser.add_option("--exeroot", dest="exeroot", default="../../run", \
                  help="Directory where the executable would be created")
parser.add_option("--n_ensemble", dest="n", default=0, \
                  help="Number of ensemble members")
parser.add_option("--case", dest="casename", default="", \
                  help="Name of case")
parser.add_option("--ens_file", dest="ens_file", default="", \
                  help="Name of samples file")
parser.add_option("--mc_ensemble", dest="mc_ensemble", default=0, \
                  help = 'Create monte carlo ensemble')
parser.add_option("--microbe", dest="microbe", default = False, action="store_true", \
                  help = 'CNP mode - initialize P pools')
parser.add_option("--postproc_file", dest="postproc_file", default="", \
                  help="Location of post_processing info")
parser.add_option("--postproc_only", dest="postproc_only", default=False, \
                  action="store_true", help='Only do post-processing')
parser.add_option("--parm_list", dest="parm_list", default='parm_list', \
                  help = 'File containing list of parameters to vary')
parser.add_option("--cnp", dest="cnp", default = False, action="store_true", \
                  help = 'CNP mode - initialize P pools')
parser.add_option("--site", dest="site", default='parm_list', \
                  help = 'Site name')
#parser.add_option("--calc_costfuction")
(options, args) = parser.parse_args()

options.n = int(options.n)
#Get number of samples from ensemble file
if (os.path.isfile(options.ens_file)):
    if (options.n == 0):
        #get # of lines
        myinput=open(options.ens_file)
        for s in myinput:
            options.n = options.n+1
        myinput.close()
else:
    if (not options.postproc_only):
      if (options.mc_ensemble > 0):
        nsamples = int(options.mc_ensemble)
        options.ens_file = 'mcsamples_'+options.casename+'.txt'
        options.n = int(options.mc_ensemble)
        print('Creating Monte Carlo ensemble with '+str(nsamples)+' members')
      else:
        print('ensemble file does not exist.  Exiting')
        sys.exit()
    else:
      print('Ensemble file not provided')
      print('Getting parameter information from output files')

#Define function to perform ensemble member post-processing
def postproc(myvars, myyear_start, myyear_end, myday_start, myday_end, myavg, \
             myfactor, myoffset, mypft, thisjob, runroot, case, pnames, ppfts, data, parms):
    rundir = options.runroot+'/UQ/'+case+'/g'+str(100000+thisjob)[1:]+'/'
    index=0
    ierr = 0
    thiscol = 0
    for v in myvars:
        ndays_total = 0
        output = []
        n_years = myyear_end[index]-myyear_start[index]+1
        for y in range(myyear_start[index],myyear_end[index]+1):
            if (mypft[index] == -1):
              fname = rundir+case+'.clm2.h0.'+str(10000+y)[1:]+'-01-01-00000.nc'
              myindex = 0
              hol_add = 1
            else:
              fname = rundir+case+'.clm2.h1.'+str(10000+y)[1:]+'-01-01-00000.nc'
              myindex = mypft[index]
              hol_add = 17
            if (os.path.exists(fname)):
              mydata = nffun.getvar(fname,v)   
            else:
              mydata = np.zeros([365,1], np.float)+np.NaN
            #get output and average over days/years
            n_days = myday_end[index]-myday_start[index]+1
            ndays_total = ndays_total + n_days
            #get number of timesteps per output file
            
            if (('20TR' in case or (not '1850' in case)) and (not 'ED' in case)):     #Transient assume daily ouput
                for d in range(myday_start[index]-1,myday_end[index]):
                    if ('US-SPR' in case):
                      output.append(0.25*(mydata[d][myindex+hol_add]*myfactor[index] \
                             +myoffset[index]) + 0.75*(mydata[d][myindex]*myfactor[index] \
                             +myoffset[index]))
                    else:
                      output.append(mydata[d][myindex]*myfactor[index] \
                             +myoffset[index])
            else:                    #Assume annual output (ignore days)
               for d in range(myday_start[index]-1,myday_end[index]):    #28-38 was myindex
                 if ('SCPF' in v):
                   output.append(sum(mydata[0,28:38])/10.0*myfactor[index]+myoffset[index])
                 else:
                   output.append(mydata[0,myindex]*myfactor[index]+myoffset[index])
        for i in range(0,ndays_total/myavg[index]):
 	    data[thiscol] = sum(output[(i*myavg[index]):((i+1)*myavg[index])])/myavg[index]
            thiscol=thiscol+1
        index=index+1
    if (options.microbe):
      pfname = rundir+'microbepar_in'
      pnum=0
      for p in pnames:
        myinput = open(pfname, 'r')
        for s in myinput:
          if (p == s.split()[0]):
            parms[pnum] = s.split()[1]
        myinput.close()
        pnum=pnum+1
    else:
      pfname = rundir+'clm_params_'+str(100000+thisjob)[1:]+'.nc'
      fpfname = rundir+'fates_params_'+str(100000+thisjob)[1:]+'.nc'
      sfname = rundir+'surfdata_'+str(100000+thisjob)[1:]+'.nc'
      pnum=0
      for p in pnames:
         if (p == 'lai'):     #Surface data file
           mydata = nffun.getvar(sfname,'MONTHLY_LAI')
           parms[pnum] = mydata[0,0,0,0]
         elif (p == 'co2'):   #CO2 value from namelist
           lnd_infile = open(rundir+'lnd_in','r')
           for s in lnd_infile:
             if ('co2_ppm' in s):
               ppmv = float(s.split()[2])
           parms[pnum] = ppmv
           lnd_infile.close()
         elif ('fates' in p):   #fates parameter file
           mydata = nffun.getvar(fpfname,p) 
           if (int(ppfts[pnum]) >= 0):
             if (p == 'fates_prt_nitr_stoich_p1'):
               parms[pnum] = mydata[int(ppfts[pnum]) % 6,int(ppfts[pnum])/6]
             else:
               parms[pnum] = mydata[int(ppfts[pnum])]
           else:
             parms[pnum] = mydata[0]
         else:                #Regular parameter file
           mydata = nffun.getvar(pfname,p) 
           if (int(ppfts[pnum]) >= 0):
             parms[pnum] = mydata[int(ppfts[pnum])]
           else:
             parms[pnum] = mydata[0]
         pnum=pnum+1

    return ierr
            

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

workdir = os.getcwd()

#get postproc info
do_postproc=False
if (os.path.isfile(options.postproc_file)):
    do_postproc=True
    myvars=[]
    myyear_start=[]
    myyear_end=[]
    myday_start=[]
    myday_end=[]
    myavg_pd=[]
    myfactor=[]
    myoffset=[]
    mypft=[]
    time.sleep(rank)
    postproc_input = open(options.postproc_file,'r')
    data_cols = 0
    for s in postproc_input:
        if (s[0:1] != '#'):
            myvars.append(s.split()[0])
            myyear_start.append(int(s.split()[1]))
            myyear_end.append(int(s.split()[2]))
            myday_start.append(int(s.split()[3]))
            myday_end.append(int(s.split()[4]))
            myavg_pd.append(int(s.split()[5]))
            myfactor.append(float(s.split()[6]))
            myoffset.append(float(s.split()[7]))
            if (len(s.split()) == 9):
              mypft.append(int(s.split()[8]))
            else:
              mypft.append(-1)
            days_total = (int(s.split()[2]) - int(s.split()[1])+1)*(int(s.split()[4]) - int(s.split()[3])+1)        
            data_cols = data_cols + days_total / int(s.split()[5])
    if (rank == 0):
        data = np.zeros([data_cols,options.n], np.float)-999
    data_row = np.zeros([data_cols], np.float)-999
    postproc_input.close()
    #get the parameter names
    pnames=[]
    ppfts=[]
    pmin=[]
    pmax=[]
    pfile = open(options.parm_list,'r')
    nparms = 0
    for s in pfile:
      pnames.append(s.split()[0])
      ppfts.append(s.split()[1])
      pmin.append(s.split()[2])
      pmax.append(s.split()[3])
      nparms = nparms+1
    pfile.close()
    parm_row = np.zeros([nparms], np.float)-999
    if (rank == 0):
      parms = np.zeros([nparms, options.n], np.float)-999
      

if (rank == 0):
    n_done=0
    if (options.mc_ensemble > 0):
      #Create a parameter samples file
      #get the parameter names
      pnames=[]
      ppfts=[]
      pmin=[]
      pmax=[]
      pfile = open(options.parm_list,'r')
      nparms = 0
      for s in pfile:
        pnames.append(s.split()[0])
        ppfts.append(s.split()[1])
        pmin.append(float(s.split()[2]))
        pmax.append(float(s.split()[3]))
        nparms = nparms+1
      pfile.close()
      nsamples = int(options.mc_ensemble)
      samples=np.zeros((nparms,nsamples), dtype=np.float)
      for i in range(0,nsamples):
          for j in range(0,nparms):
              samples[j][i] = pmin[j]+(pmax[j]-pmin[j])*np.random.rand(1)
      np.savetxt('mcsamples_'+options.casename+'.txt', np.transpose(samples))
    #send first np-1 jobs where np is number of processes
    for n_job in range(1,size):
        comm.send(n_job, dest=n_job, tag=1)
        comm.send(0,     dest=n_job, tag=2)
        if (options.postproc_only == False):
            time.sleep(0.2)
    #Assign rest of jobs on demand
    for n_job in range(size,options.n+1):
        process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
        thisjob = comm.recv(source=process, tag=4)
        if (do_postproc):
            data_row = comm.recv(source=process, tag=5)
            data[:,thisjob-1] = data_row
            parm_row = comm.recv(source=process, tag=6)
            parms[:,thisjob-1] = parm_row
        print 'Received', thisjob
        n_done = n_done+1
        comm.send(n_job, dest=process, tag=1)
        comm.send(0,     dest=process, tag=2)
    #receive remaining messages and finalize
    while (n_done < options.n):
        process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
        thisjob = comm.recv(source=process, tag=4)
        if (do_postproc):
            data_row = comm.recv(source=process, tag=5)
            data[:,thisjob-1] = data_row
            parm_row = comm.recv(source=process, tag=6)
            parms[:,thisjob-1] = parm_row
        print 'Received', thisjob
        n_done = n_done+1
        comm.send(-1, dest=process, tag=1)
        comm.send(-1, dest=process, tag=2)
    
    if (do_postproc):
        data_out = data.transpose()
        parm_out = parms.transpose()
        good=[]
        for i in range(0,options.n):
          #only save valid runs (no NaNs)
          if not np.isnan(sum(data_out[i,:])):
            good.append(i)
        data_out = data_out[good,:]
        parm_out = parm_out[good,:]
        np.savetxt(options.casename+'_postprocessed.txt', data_out)
        #UQ-ready outputs (80% of data for traning, 20% for validation)
        UQ_output = 'UQ_output/'+options.casename
        os.system('mkdir -p '+UQ_output)
        np.savetxt(UQ_output+'/ytrain.dat', data_out[0:int(len(good)*0.8),:])
        np.savetxt(UQ_output+'/yval.dat',   data_out[int(len(good)*0.8):,:])
        np.savetxt(UQ_output+'/ptrain.dat', parm_out[0:int(len(good)*0.8),:])
        np.savetxt(UQ_output+'/pval.dat', parm_out[int(len(good)*0.8):,:])
        myoutput = open(UQ_output+'/pnames.txt', 'w')
        eden_header=''
        for p in pnames:
          myoutput.write(p+'\n')
          eden_header=eden_header+p+','
        myoutput.close()
        myoutput = open(UQ_output+'/outnames.txt', 'w')
        for v in myvars:
          myoutput.write(v+'\n')
          eden_header=eden_header+v+','
        myoutput.close()
        myoutput = open(UQ_output+'/param_range.txt', 'w')
        for p in range(0,len(pmin)):
          myoutput.write(pmin[p]+' '+pmax[p]+'\n')
        myoutput.close()
        print np.hstack((parm_out,data_out))
        np.savetxt(UQ_output+'/foreden.csv', np.hstack((parm_out,data_out)), delimiter=',', header=eden_header[:-1])

    MPI.Finalize()

#Slave
else:
    status=0
    while status == 0:
        myjob = comm.recv(source=0, tag=1)
        status = comm.recv(source=0, tag=2) 

        if (status == 0):
            if (options.postproc_only == False):
                os.chdir(workdir)
                cnp = 'False'
                if (options.cnp):
                    cnp='True'
                #Python script to set up the ensemble run directory and manipulate parameters
                os.system('python ensemble_copy.py --case '+options.casename+' --runroot '+ \
                      options.runroot +' --ens_num '+str(myjob)+' --ens_file '+options.ens_file+ \
                      ' --parm_list '+options.parm_list+' --cnp '+cnp+' --site '+options.site)
                jobst = str(100000+int(myjob))
                rundir = options.runroot+'/UQ/'+options.casename+'/g'+jobst[1:]+'/'
                os.chdir(rundir)
                #Run the executable
                exedir = options.exeroot
                if os.path.isfile(exedir+'/acme.exe'):
                   os.system(exedir+'/acme.exe > acme_log.txt')
                elif os.path.isfile(exedir+'/e3sm.exe'):
                   os.system(exedir+'/e3sm.exe > e3sm_log.txt')
                elif os.path.isfile(exedir+'/cesm.exe'):
                   os.system(exedir+'/cesm.exe > cesm_log.txt')
            if (do_postproc):
                ierr = postproc(myvars, myyear_start, myyear_end, myday_start, \
                         myday_end, myavg_pd, myfactor, myoffset, mypft, myjob, \
                         options.runroot, options.casename, pnames, ppfts, data_row, parm_row)
                comm.send(rank, dest=0, tag=3)
                comm.send(myjob, dest=0, tag=4)
                comm.send(data_row, dest=0, tag=5)
                comm.send(parm_row, dest=0, tag=6)
            else:
                comm.send(rank,  dest=0, tag=3)
                comm.send(myjob, dest=0, tag=4)

    print rank, ' complete'
    MPI.Finalize()
