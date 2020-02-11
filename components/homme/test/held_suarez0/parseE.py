#!/usr/bin/env python
from string import *
import os, commands, getopt, sys, exceptions
import numpy as np
import matplotlib.pyplot as plt

class search_failed(exceptions.Exception):
    def __init__(self,args=None):
        self.args=args
class eof(exceptions.Exception):
    def __init__(self,args=None):
        self.args=args



#
# look for key1.  if we find key2 instead, raise error
#
def lookfor1(fid,key1,key2="",allow_eof=0):
    line=fid.readline()
    while line:
        pos=find(line,key1)
        if (-1 <> pos ):
            sline = line[pos+len(key1):-1]
            return sline
        if (len(key2)>0):
           pos=find(line,key2)
           if (-1 <> pos ):
               print "error looking for: "+key1
               raise search_failed,"run not complete, found: "+key2
        line=fid.readline()

    if (allow_eof==1):
        raise eof,"EOF"
#    raise search_failed,"Search failed  string='"+key1+"'"
    raise search_failed,"Search failed  string="+key1
    return 0

    

mumin=[]
mumax=[]
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
    str=split(str)
    ncpu=atoi(str[0])

    str = lookfor1(sys.stdin,"theta_hydrostatic_mode ","",0)
    str=split(str)
    hydrostatic_mode  = (str[1]=="T")

    str = lookfor1(sys.stdin,"tstep ","",0)
    str=split(str)
    tstep=atof(str[1])
    print 'NCPU = %i tstep=%f' % (ncpu, tstep)
    while 1:
        # nstep=           3  time=  1.041666666666667E-002  [day]
        str = lookfor1(sys.stdin,"nstep=","",1)
        str=split(str)
        n=atoi(str[0])
     
        time.extend([n*tstep/(24*3600)])

        if ( ~hydrostatic_mode ):
            str = lookfor1(sys.stdin,"mu    =","",0)
            str=split(str)
            mumin.extend([atof(str[0])])
            mumax.extend([atof(str[3])])

        # KE,d/dt,diss:
        # IE,d/dt,diss:
        # PE,d/dt,diss:
        str = lookfor1(sys.stdin,"KE,d/dt","",0)
        str=split(str)
        KE.extend([atof(str[1])])
        if (len(str) >= 4):
           KEdiss.extend([atof(str[3])])

        str = lookfor1(sys.stdin,"IE,d/dt","",0)
        str=split(str)
        IE.extend([atof(str[1])])
        if (len(str) >= 4):
           IEdiss.extend([atof(str[3])])

        str = lookfor1(sys.stdin,"PE,d/dt","",0)
        str=split(str)
        PE.extend([atof(str[1])])
        if (len(str) >= 4):
           PEdiss.extend([atof(str[3])])

        if ( hydrostatic_mode ):
           str = lookfor1(sys.stdin,"I+P,d/dt","",0)
           str=split(str)
           IPEdiss.extend([atof(str[3])])

        str = lookfor1(sys.stdin," E,d/dt","",0)
        str=split(str)
        DELE.extend([atof(str[2])])

        #print 'parsed nstep=%i' % ( n)

        
except search_failed,e:
    print "search failed"
    print "".join(e)
    sys.exit(1)
    
except eof,e:
    print 'plotting energy...'
    KE2= np.array(KE)
    nlen=KE2.size
    print 'data parsed size=%i' % (nlen)
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

    print 'plotting dissipation rates...'    
    print 'data parsed size=%i' % (len(IPEdiss))
    plt.figure()
    if (len(IPEdiss)>0):
       print 'plotting Hydrostatic I+P,d/dt,diss data'
       plt.plot(time,IPEdiss,label='IE+PE dissipation')

    if (len(IEdiss)>0):
       print 'plotting NH data, sum of IEdiss and PEdiss'
       IPEdiss = np.array(PEdiss) + np.array(IEdiss)       
       plt.plot(time,IPEdiss,label='IE+PE dissipation')

    if (len(KEdiss)>0): 
       plt.plot(time,KEdiss,label='KE dissipation')

    plt.plot(time,DELE,label='TOT E dissipation')

    plt.axis([0, 500, -.5, .1])
    #plt.axis([0, 500, -.1, .1])
    #plt.axis([1600, 1700, -.4, .2])
    plt.grid(True)
    plt.legend()
    plt.savefig("HS-diss.png")

    if ( ~hydrostatic_mode ):
        plt.figure()
        print ('plotting mu...std min,max=%f %f' % (np.std(mumin),np.std(mumax)))
        legend1=("min avg: %.3f std: %.4f" % (sum(mumin)/len(mumin),np.std(mumin)) )
        plt.plot(time,mumin,label=legend1)
        legend1=("min avg: %.3f std: %.4f" % (sum(mumax)/len(mumax),np.std(mumax)) )
        plt.plot(time,mumax,label=legend1)
        plt.axis([min(time), max(min(time)+200,max(time)), -0.,2.0])
        plt.grid(True)
        plt.legend()
        plt.savefig("mu.png")
    

    plt.show()
    sys.exit(0)

