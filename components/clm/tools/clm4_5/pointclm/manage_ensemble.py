#!/usr/bin/env python
import sys,os, time
import numpy as np
import netcdf_functions as nffun
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
        for s in input:
            options.n = options.n+1
        myinput.close()
else:
    print('ensemble file does not exist.  Exiting')
    sys.exit()

#Define function to perform ensemble member post-processing
def postproc(myvars, myyear_start, myyear_end, myday_start, myday_end, \
             myavg, myfactor, myoffset, thisjob, runroot, case, data, outnum):
    rundir = options.runroot+'/UQ/'+case+'/g'+str(100000+thisjob)[1:]+'/'
    index=0
    outnum = 0
    ierr = 0
    for v in myvars:
        ndays_total = 0
        output = []
        n_years = myyear_end[index]-myyear_start[index]+1
        for y in range(myyear_start[index],myyear_end[index]+1):
            fname = rundir+case+'.clm2.h0.'+str(10000+y)[1:]+'-01-01-00000.nc'
            mydata = nffun.getvar(fname,v)   
            #get output and average over days/years
            n_days = myday_end[index]-myday_start[index]+1
            ndays_total = ndays_total + n_days
            if ('20TR' in case):     #Transient assume daily ouput
                for d in range(myday_start[index]-1,myday_end[index]):
                    output.append(mydata[d][0]*myfactor[index] \
                             +myoffset[index])
            else:                    #Assume annual output (ignore days)
               for d in range(myday_start[index]-1,myday_end[index]):
                    output.append(mydata*myfactor[index]+myoffset[index])
        for i in range(0,ndays_total/myavg[index]):
            print v, index, i
	    data[outnum,thisjob-1] = sum(output[(i*myavg[index]):((i+1)*myavg[index])])/myavg[index]
            outnum=outnum+1
        index = index+1
    return ierr
            

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

workdir = os.getcwd()

#Master
if (rank == 0):
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
                days_total = (int(s.split()[2]) - int(s.split()[1])+1)*(int(s.split()[4]) - int(s.split()[3])+1)        
                data_cols = data_cols + days_total / int(s.split()[5])
                print data_cols
        data=np.zeros([data_cols,options.n], np.float)-999
        postproc_input.close()

    n_done=0
    #send first np-1 jobs where np is number of processes
    for n_job in range(1,size):
        comm.send(n_job, dest=n_job, tag=1)
        comm.send(0,     dest=n_job, tag=2)
        time.sleep(1)
    #Assign rest of jobs on demand
    for n_job in range(size,options.n+1):
        process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
        thisjob = comm.recv(source=process, tag=4)
        n_done = n_done+1
        outnum=0
        if (do_postproc):
            ierr = postproc(myvars, myyear_start, myyear_end, myday_start, \
                     myday_end, myavg_pd, myfactor, myoffset, thisjob, \
                     options.runroot, options.casename, data, outnum)
        print outnum
        comm.send(n_job, dest=process, tag=1)
        comm.send(0,     dest=process, tag=2)
    #receive remaining messages and finalize
    while (n_done < options.n):
        process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
        thisjob = comm.recv(source=process, tag=4)
        n_done = n_done+1
        outnum = 0
        if (do_postproc):
            ierr = postproc(myvars, myyear_start, myyear_end, myday_start, \
                            myday_end, myavg_pd, myfactor, myoffset, thisjob, \
                            options.runroot, options.casename, data, outnum)
            #Update post-processed output file
        comm.send(-1, dest=process, tag=1)
        comm.send(-1, dest=process, tag=2)
    
    if (do_postproc):
        data_out = data.transpose()
        np.savetxt(options.casename+'_postprocessed.txt', data_out)
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
                os.system(options.exeroot+'/acme.exe > acme.log')
            comm.send(rank,  dest=0, tag=3)
            comm.send(myjob, dest=0, tag=4)
    print rank, ' complete'
    MPI.Finalize()
