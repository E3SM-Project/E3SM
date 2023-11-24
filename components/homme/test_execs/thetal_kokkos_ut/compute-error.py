#!/usr/bin/env python3

import os, sys, re, math
#from array import array
import numpy as np

basef="baseline-output"
outf="output"

patt="dp["
s1="NP"
s2="NELEM"
s3="NUM_PNHYS_LEVELS"

nelem1=0
np1=0
nlev1=0

dparray=[]

def read_dims(s1, s2, s3, ffile):

 np=0
 nelem=0
 nlev=0
 with open(ffile, 'r') as f:

 #read dimensions
  for ln in f:
     if s1 in ln: 
        #print(ln)
        np = int(ln.split()[2])
        #print(np)
  f.seek(0)
  for ln in f:
     if s2 in ln:   
        #print(ln)
        nelem = int(ln.split()[2])
  f.seek(0)
  for ln in f:
     if s3 in ln:   
        #print(ln)
        nlev = int(ln.split()[2])
  f.seek(0)

  print("dimensions: nelem="+str(nelem) + " np=" + str(np) + " nlev="  + str(nlev))
  return (nelem,np,nlev)


def read_field(patt,ffile):

 dparray = [];

 with open(ffile, 'r') as f:
  for ln in f:
     if patt in ln:
        #print(ln)
        dpval = float(ln.split()[1])
        dparray.extend([dpval])
        #print (dpval)
  f.seek(0)
  return dparray


##################  main

[nelem1,np1,nlev1]=read_dims(s1,s2,s3,basef)
[nelem2,np2,nlev2]=read_dims(s1,s2,s3,outf)

if( nelem1 != nelem2 ):
  sys.exit("nelem dim does not match")
if( np1 != np2 ):
  sys.exit("np dim does not match")
if( nlev1 != nlev2 ):
  sys.exit("nlev dim does not match")

dp1 = read_field(patt,basef)
dp2 = read_field(patt,outf)

dd1=np.array(dp1)
dd2=np.array(dp2)

res=math.sqrt(sum( (dd1-dd2)*(dd1-dd2) )) / math.sqrt(sum( dd1*dd1 ))

print (res)


