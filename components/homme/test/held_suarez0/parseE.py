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
time=[]
try:
    startstr="number of MPI processes:"
    str = lookfor1(sys.stdin,startstr,"",0)
    str=split(str)
    ncpu=atoi(str[0])

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

        str = lookfor1(sys.stdin,"IE,d/dt","",0)
        str=split(str)
        IE.extend([atof(str[1])])

        str = lookfor1(sys.stdin,"PE,d/dt","",0)
        str=split(str)
        PE.extend([atof(str[1])])
        #print 'parsed nstep=%i' % ( n)

        
except search_failed,e:
    print "search failed"
    print "".join(e)
    sys.exit(1)
    
except eof,e:
    print 'plotting...'
    KE2= np.array(KE)
    nlen=KE2.size
    print 'data parsed size=%i' % (nlen)
    time=time[0:nlen]
    #KE2=KE2/KE2[0]
    KE2=KE2*1e3
    PE2= np.array(PE)
    #PE2=PE2/PE2[0]
    IE2= np.array(IE)
    #IE2=IE2/IE2[0]
    p1,=plt.plot(time,KE2,label='KE*1e3')
    p2,=plt.plot(time,PE2,label='PE')
    p2,=plt.plot(time,IE2,label='IE')
    plt.axis([0, 500, 0, 2.5e9])
    plt.legend()
    plt.show()

    sys.exit(0)

