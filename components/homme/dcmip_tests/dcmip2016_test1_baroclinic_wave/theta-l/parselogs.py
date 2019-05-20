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

    


KE=[]
IE=[]
PE=[]
KEdiss=[]
IEdiss=[]
PEdiss=[]
IPEdiss=[]
TOTE=[]
DELE=[]
time=[]
try:
    startstr="number of MPI processes:"
    str = lookfor1(sys.stdin,startstr,"",0)
    str=split(str)
    ncpu=atoi(str[0])

    str = lookfor1(sys.stdin,"nmax ","",0)
    str=split(str)
    nmax=atof(str[1])

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

        # KE,d/dt,diss:
        # IE,d/dt,diss:
        # PE,d/dt,diss:
        str = lookfor1(sys.stdin,"KE,d/dt","",1)
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
        TOTE.extend([atof(str[2])])

        str = lookfor1(sys.stdin,"E-E0","",0)
        str=split(str)
        DELE.extend([atof(str[1])])

        #print 'parsed nstep=%i' % ( n)

        
except search_failed,e:
    print "search failed"
    print "".join(e)
    sys.exit(1)
    
except eof,e:
    KE2= np.array(KE)
    nlen=KE2.size
    print 'data parsed size=%i nmax=%i' % (nlen,nmax)

    # data from 2nd diags outptu is skipped since E-E0 is missing
    # nstep>=1 (3rd output) is correct

    nstep=1
    time=(24*3600*15 + tstep*(nstep+2))/3600./24.
    print '%i %.4f %e %e %e %e' % (tstep,time, KEdiss[nstep],IEdiss[nstep],PEdiss[nstep],DELE[nstep] )

    nstep=nlen-1
    time=(24*3600*15 + tstep*(nstep+2))/3600./24.
    print '%i %.4f %e %e %e %e' % (tstep,time, KEdiss[nstep],IEdiss[nstep],PEdiss[nstep],DELE[nstep] )

    sys.exit(0)

