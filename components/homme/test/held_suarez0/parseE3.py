#!/usr/bin/env python3
from string import *
import os, getopt, sys, builtins
import numpy as np
import matplotlib.pyplot as plt

class search_failed(Exception):
    def __init__(self,args=None):
        self.args=args
class eof(Exception):
    def __init__(self,args=None):
        self.args=args




#
# look for key1.  if we find key2 instead, raise error
#
def lookfor1(fid,key1,key2="",allow_eof=0):
    line=fid.readline()
    while line:
        pos=line.find(key1)
        if (-1 != pos ):
            sline = line[pos+len(key1):-1]
            return sline
        if (len(key2)>0):
           pos=line.find(key2)
           if (-1 != pos ):
               #print ("error looking for: "+key1)
               #raise (search_failed,"run not complete, found: "+key2)
               sline="0 0 0 0 0 0 0 0 0 0"
               return sline
        line=fid.readline()

    if (allow_eof==1):
        raise eof("EOF")
    raise search_failed("Search failed:  string='"+key1+"'")
    return 0

showplot=1
if len(os.sys.argv) >= 2:
   if os.sys.argv[1]=="-noX11":
       showplot=0

mumin=[]
mumax=[]
Tmin=[]
Tmax=[]
KE=[]
IE=[]
PE=[]
KEdiss=[]
IEdiss=[]
PEdiss=[]
IPEdiss=[]
DELE=[]
time=[]
try:
    startstr="number of MPI processes:"
    str = lookfor1(sys.stdin,startstr,"",0)
    str=str.split()
    ncpu=int(str[0])

    str = lookfor1(sys.stdin,"theta_hydrostatic_mode ","",0)
    str=str.split()
    hydrostatic_mode  = (str[1]=="T")

    str = lookfor1(sys.stdin,"tstep ","",0)
    str=str.split()
    tstep=float(str[1])
    print('NCPU = %i tstep=%f NH=%i' % (ncpu, tstep, not hydrostatic_mode))
    while 1:
        # nstep=           3  time=  1.041666666666667E-002  [day]
        str = lookfor1(sys.stdin,"nstep=","",1)
        str=str.split()
        n=int(str[0])
     
        time.extend([n*tstep/(24*3600)])

        # look for mu. return all zeros if we find dz(m) first:
        str = lookfor1(sys.stdin,"mu    =","dz(m) =",0)
        str=str.split()
        mumin.extend([float(str[0])])
        mumax.extend([float(str[3])])

        # look for TBOT. return all zeros if we find ps first:
        str = lookfor1(sys.stdin,"TBOT=","ps=",0)
        str=str.split()
        Tmin.extend([float(str[0])])
        Tmax.extend([float(str[1])])

        # KE,d/dt,diss:
        # IE,d/dt,diss:
        # PE,d/dt,diss:
        str = lookfor1(sys.stdin,"KE,d/dt","",0)
        str=str.split()
        KE.extend([float(str[1])])
        if (len(str) >= 4):
           KEdiss.extend([float(str[3])])

        str = lookfor1(sys.stdin,"IE,d/dt","",0)
        str=str.split()
        IE.extend([float(str[1])])
        if (len(str) >= 4):
           IEdiss.extend([float(str[3])])

        str = lookfor1(sys.stdin,"PE,d/dt","",0)
        str=str.split()
        PE.extend([float(str[1])])
        if (len(str) >= 4):
           PEdiss.extend([float(str[3])])

        if ( hydrostatic_mode ):
           str = lookfor1(sys.stdin,"I+P,d/dt","",0)
           str=str.split()
           IPEdiss.extend([float(str[3])])

        str = lookfor1(sys.stdin," E,d/dt","",0)
        str=str.split()
        DELE.extend([float(str[2])])

        #print 'parsed nstep=%i' % ( n)

        
except search_failed as e:
    print(''.join(e.args))
    sys.exit(1)
    
except eof as e:
    print('eof reached. analysing data.')

    print('plotting energy...')
    KE2= np.array(KE)
    nlen=KE2.size
    print('data parsed size=%i' % (nlen))
    time=time[0:nlen]
    KE2=KE2*1e3
    PE2= np.array(PE)
    IE2= np.array(IE)
    plt.figure()
    plt.plot(time,KE2,label='KE*1e3')
    plt.plot(time,PE2,label='PE')
    plt.plot(time,IE2,label='IE')
    #plt.axis([0, 500, 0, 2.5e9])
    plt.grid(True)
    plt.legend()
    plt.savefig("HS-E.png")

    print('plotting dissipation rates...'    )
    print('data parsed size=%i' % (len(IPEdiss)))
    plt.figure()
    if (len(IPEdiss)>0):
       print('plotting Hydrostatic I+P,d/dt,diss data')
       plt.plot(time,IPEdiss,label='IE+PE dissipation')

    if (len(IEdiss)>0):
       print('plotting NH data, sum of IEdiss and PEdiss')
       IPEdiss = np.array(PEdiss) + np.array(IEdiss)       
       plt.plot(time,IPEdiss,label='IE+PE dissipation')

    if (len(KEdiss)>0): 
       plt.plot(time,KEdiss,label='KE dissipation')

    plt.plot(time,DELE,label='TOT E dissipation')

    plt.axis([0, 700, -.5, .5])
    #plt.axis([0, 500, -.1, .1])
    #plt.axis([1600, 1700, -.4, .2])
    plt.grid(True)
    plt.legend()
    plt.savefig("HS-diss.png")

    if ( not hydrostatic_mode ):
        plt.figure()
        print ('plotting mu...std min,max=%f %f' % (np.std(mumin),np.std(mumax)))
        print ('avg min,max=%f %f' % (sum(mumin)/len(mumin),sum(mumax)/len(mumax)))
        legend1=("min avg: %.3f std: %.4f" % (sum(mumin)/len(mumin),np.std(mumin)) )
        plt.plot(time,mumin,label=legend1)
        legend1=("max avg: %.3f std: %.4f" % (sum(mumax)/len(mumax),np.std(mumax)) )
        plt.plot(time,mumax,label=legend1)
        plt.axis([min(time), max(min(time)+200,max(time)+10), -1.0,4])
        plt.grid(True)
        plt.legend()
        plt.savefig("mu.png")

    plt.figure()
    print ('plotting TBOT..std min,max=%f %f' % (np.std(Tmin),np.std(Tmax)))
    print ('min,max=%f %f' % (min(Tmin),max(Tmax)))
    legend1=("min: %.3f std: %.4f" % (min(Tmin),np.std(Tmin)) )
    plt.plot(time,Tmin,label=legend1)
    legend1=("max: %.3f std: %.4f" % (max(Tmax),np.std(Tmax)) )
    plt.plot(time,Tmax,label=legend1)
    plt.axis([min(time), max(min(time)+200,max(time)+10), 0,400])
    plt.grid(True)
    plt.legend()
    plt.savefig("TBOT.png")
   
    if showplot:
        plt.show()
    sys.exit(0)

