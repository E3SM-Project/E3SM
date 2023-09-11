#!/usr/bin/env python
from string import *
import os, commands, getopt, sys, exceptions
import numpy as np
import matplotlib.pyplot as plt


#call is as
#python parse-5files.py ../../output-cpstarnh-l128-dt0.0001 ../../output-cpstarhy-l128-dt0.0001 ../../output-cv-l128-dt0.0001 ../../output-eamcpstar-l128-dt0.0001 ../../output-eamcpdry-l128-dt0.0001 "const pressure" "const pressure HY" "const volume" "EAM cpstar" "EAM cpdry"

#first and last plots -- total precip and cl energy are set, the rest is not


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
        #raise eof,"EOF"
        return "EOF"
    raise search_failed,"Search failed  string="+key1
    return 0


############################## plot ###############################
### assume all 3 arrays have the same length, indarray too
### start and end are indices, not lengths
def plot4arrays(arr1, arr2, arr3, arr4, arr1name, arr2name, arr3name, arr4name, plotname, FS=10, ttitle="none"):

  s1=arr1.size
  s2=arr2.size
  s3=arr3.size
  s4=arr4.size

  #if( s1>=startt and s1>=endd ): #DO NOT CHECK SIZES
  plt.figure()
  b1=arr1
  b2=arr2
  b3=arr3
  b4=arr4

  legend0=(arr1name+": avg %g , std: %g" % (sum(b1)/len(b1),np.std(b1)) )
  legend1=(arr2name+": avg %g , std: %g" % (sum(b2)/len(b2),np.std(b2)) )
  legend2=(arr3name+": avg %g , std: %g" % (sum(b3)/len(b3),np.std(b3)) )
  legend3=(arr4name+": avg %g , std: %g" % (sum(b4)/len(b4),np.std(b4)) )
  #print(log(abs(b1)))

  plt.semilogy(abs(b1)*100.0, alpha=0.5,  color='green',label=legend0)
  print("aaa")
  plt.semilogy(abs(b2)*100.0, alpha=0.5,  color='blue',label=legend1)
  print("222")
  plt.semilogy(abs(b3)*100.0, alpha=0.5,  color='red',label=legend2)
  print("333")
  plt.semilogy(abs(b4)*100.0, alpha=0.5,  color='black',label=legend3)
  print("444")

  plt.legend(fontsize=FS)
  plt.title(ttitle,fontsize=FS)
  plt.tick_params(labelsize=FS)
  plt.grid(True)

  plt.savefig(plotname+".pdf")

  return;


############### START MAIN


# p Mass [kg/m2/sec]   798.000000000000       9.506296313586071E-004
# p En [W/m2/sec]   798.000000000000        1289.41111562664
# p cl En [W/m2/sec]   798.000000000000        1511.59617682332
# leak En [W/m2/sec]   798.000000000000       3.227660888293464E-009
# leak M [ks/m2/sec]   798.000000000000       2.142152774586057E-014

latvap=2.501e6
latice=3.337e5
Tforcl=290.0
cl=4188.0


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
print 'file1 is', str(sys.argv[1])
print 'file2 is', str(sys.argv[2])
print 'file3 is', str(sys.argv[3])
print 'file4 is', str(sys.argv[4])
print 'file5 is', str(sys.argv[5])
print 'LABEL1: 6 is', str(sys.argv[6])
print 'LABEL2: 7 is', str(sys.argv[7])
print 'LABEL3: 8 is', str(sys.argv[8])
print 'LABEL4: 9 is', str(sys.argv[9])
print 'LABEL5: 10 is', str(sys.argv[10])
#print 'NAME OF OUTPUT PLOT: 11 is', str(sys.argv[11])
#print 'PATTERN: 12 is', str(sys.argv[12])
#will be used to decide what to plot
#print 'what: 13 is', str(sys.argv[13])

label1=str(sys.argv[6])
label2=str(sys.argv[7])
label3=str(sys.argv[8])
label4=str(sys.argv[9])
label5=str(sys.argv[10])
#plotname=str(sys.argv[11])
#pattern=str(sys.argv[12])
#ttitle=str(sys.argv[13])

stri1="p Mass"
stri2="p En"
stri3="p cl En"
stri4="leak En"

pmass1=[]
pmass2=[]
pmass3=[]
pmass4=[]
pmass5=[]
nstep1=[]
nstep2=[]
nstep3=[]
nstep4=[]
nstep5=[]

pen1=[]
pen2=[]
pen3=[]
pen4=[]
pen5=[]

pclen1=[]
pclen2=[]
pclen3=[]
pclen4=[]
pclen5=[]

pleake1=[]
pleake2=[]
pleake3=[]
pleake4=[]
pleake5=[]


try:
   fff1=open(str(sys.argv[1]),"r")
   fff2=open(str(sys.argv[2]),"r")
   fff3=open(str(sys.argv[3]),"r")
   fff4=open(str(sys.argv[4]),"r")
   fff5=open(str(sys.argv[5]),"r")
finally:

   ff=fff1 #1st file

   noteof=1
   while noteof :
      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        nstep1.extend([atof(strr[1])])
        pmass1.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri2,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pen1.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri3,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pclen1.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri4,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pleake1.extend([atof(strr[2])])
      else :
        noteof=0
   print "finished readfing the 1st file"

   ff=fff2 #file

   noteof=1
   while noteof :
      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        nstep2.extend([atof(strr[1])])
        pmass2.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri2,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pen2.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri3,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pclen2.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri4,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pleake2.extend([atof(strr[2])])
      else :
        noteof=0
   print "finished readfing the file 2"

   ff=fff3 #file

   noteof=1
   while noteof :
      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        nstep3.extend([atof(strr[1])])
        pmass3.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri2,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pen3.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri3,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pclen3.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri4,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pleake3.extend([atof(strr[2])])
      else :
        noteof=0
   print "finished readfing the file 3"


   ff=fff4 #file

   noteof=1
   while noteof :
      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        nstep4.extend([atof(strr[1])])
        pmass4.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri2,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pen4.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri3,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pclen4.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri4,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pleake4.extend([atof(strr[2])])
      else :
        noteof=0
   print "finished readfing the file 4"

   ff=fff5 #file

   noteof=1
   while noteof :
      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        nstep5.extend([atof(strr[1])])
        pmass5.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri2,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pen5.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri3,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pclen5.extend([atof(strr[2])])
      else :
        noteof=0
   ff.seek(0)
   noteof=1
   while noteof :
      strr=lookfor1(ff,stri4,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        pleake5.extend([atof(strr[2])])
      else :
        noteof=0
   print "finished readfing the file 5"




   ####################################################

   #f1=(np.array(te_before_dme1)-np.array(te_after_dycore1))/dtime
   ns1=np.array(nstep1)
   ns2=np.array(nstep2)
   ns3=np.array(nstep3)
   ns4=np.array(nstep4)
   ns5=np.array(nstep5)

   nn1=ns1.size
   print "size1="+str(nn1)

   apmass1=np.array(pmass1)
   apen1=np.array(pen1)
   apclen1=np.array(pclen1)
   apleake1=np.array(pleake1)

   apmass2=np.array(pmass2)
   apen2=np.array(pen2)
   apclen2=np.array(pclen2)
   apleake2=np.array(pleake2)

   apmass3=np.array(pmass3)
   apen3=np.array(pen3)
   apclen3=np.array(pclen3)
   apleake3=np.array(pleake3)

   apmass4=np.array(pmass4)
   apen4=np.array(pen4)
   apclen4=np.array(pclen4)
   apleake4=np.array(pleake4)

   apmass5=np.array(pmass5)
   apen5=np.array(pen5)
   apclen5=np.array(pclen5)
   apleake5=np.array(pleake5)

   print "arr 1",apmass1[0:400]


   #print "First log, fixer: avg and std",sum(fixer1)/len(fixer1),np.std(fixer1) 
   #print "second log, fixer: avg and std",sum(fixer2)/len(fixer2),np.std(fixer2) 


######################################################################################

#   plot4arrays(arr1, arr2, arr3, arr4, label1, label2, label3, label4, "RR-"+plotname, 10, "energy imbalance")
#   plot4arrays(aww1, aww2, aww3, aww4, label1, label2, label3, label4, "WATER-"+plotname, 10, "water imbalance")

#        !above neither mass not energy have factor 1/gravit
#        prect_mass_forint(i,j,ie) = mass_prect/gravit/dt
#        prect_en_forint(i,j,ie) = energy_prect/gravit/dt
#        cl_en_forint(i,j,ie) = encl/gravit/dt
#        !make sure sign of dE matches sigh of energy_prect
#        en_leak_forint(i,j,ie) = (en1glob_cp - en2glob_cp - energy_prect)/gravit/dt
#        !make sure sign of dM matches sigh of mass_prect
#        mass_leak_forint(i,j,ie) = (mass1global - mass2global - mass_prect)/gravit/dt

#to get EN1-EN2, do pleake + pen, (energy_leak)+en_precip


   plotname="mass-for5"
   plt.figure()

   bblue='#377eb8'
   oorange='#ff7f00'
   ggreen='#4daf4a'
   ppink='#f781bf'
   bbrown='#a65628'
   ppurple='#984ea3'
   ggray='#999999'
   rred='#e41a1c'
   yyellow='#dede00'

   #legend0=(arr1name+": avg %g , std: %g" % (sum(b1)/len(b1),np.std(b1)) )
   #legend1=(arr2name+": avg %g , std: %g" % (sum(b2)/len(b2),np.std(b2)) )
   #legend2=(arr3name+": avg %g , std: %g" % (sum(b3)/len(b3),np.std(b3)) )
   #legend3=(arr4name+": avg %g , std: %g" % (sum(b4)/len(b4),np.std(b4)) )
   #print(log(abs(b1)))

   plt.plot(ns1,apmass1,  color=rred,  label=label1, alpha=1, linewidth=3)
   plt.plot(ns1,apmass2,  color='black',label=label2, alpha=1, linewidth=1)
   plt.plot(ns1,apmass3,  color=ggreen,label=label3, alpha=1, linewidth=2)
   plt.plot(ns1,apmass4,  color=yyellow,label=label4, alpha=1, linewidth=3)
   plt.plot(ns1,apmass5,  color=bblue, label=label5, alpha=1, linewidth=1)

   FS=14
   plt.legend(fontsize=FS)
   plt.title("precipitation, kg/m^2/sec",fontsize=FS)
   
   plt.tick_params(labelsize=FS)
   plt.grid(True)

   plt.xlabel("simulated time, sec", fontsize=FS)
   #plt.ylabel("", fontsize=FS)
   #plt.ylim([-0.0005,0.0105])

   plt.savefig(plotname+".pdf")
   plt.savefig(plotname+".png")


   #print *, 'CL no L, (encl-enprect)/encl rel', (en1glob_cp-en2glob_cp-latice*mass_prect)/(encl-latice*mass_prect)
   rr=300
   plotname="comparecl-for5"
   plt.figure()

   #this does not work exactly the same way for 3 HY updates
   me1 = (apen1[rr:nn1-1] - latice*apmass1[rr:nn1-1])/(apclen1[rr:nn1-1] - latice*apmass1[rr:nn1-1])
   me3 = (apen3[rr:nn1-1] - latice*apmass3[rr:nn1-1])/(apclen3[rr:nn1-1] - latice*apmass3[rr:nn1-1])

   plt.plot(ns1[rr:nn1-1],me1,  color='green',label=label1)
   plt.plot(ns1[rr:nn1-1],me3,  color='brown',label=label3)

   FS=14
   plt.legend(fontsize=FS)
   plt.title("deviation from liquid precip. energy",fontsize=FS)

   plt.tick_params(labelsize=FS)
   plt.grid(True)

   plt.xlabel("simulated time, sec", fontsize=FS)
   #plt.ylabel("", fontsize=FS)
   #plt.ylim([-0.0005,0.0105])

   plt.savefig(plotname+".pdf")
   plt.savefig(plotname+".png")

   #################333

   #this plot needs removal of L terms

   #from the code, each of 5 updates compure encl with latice*mass term

   #also, from the code, each of 5 updates computes energy_precip with latice*mass part, too
   #thus we just remove latice*mass part from encl (which is then IEFLX in the paper)
   #and from E_before - E_after
   #in other words, consider exact cases where en_precip = energy_before-energy_after=cpTerm + latice*mass
   #since we want to compare only cpTerms, we need to substact from en_before-en_after latice term

   plotname="encl-diagn-for5"
   plt.figure()

   clcolor=ppurple
   #plt.plot(ns1,apclen1 - latice*apmass1, color=clcolor)
   #plt.plot(ns1,apclen2 - latice*apmass2, color=clcolor)
   #plt.plot(ns1,apclen3 - latice*apmass3, color=clcolor)
   #plt.plot(ns1,apclen4 - latice*apmass4, color=clcolor)
   #plt.plot(ns1,apclen5 - latice*apmass5, color=clcolor, label='E precipitation')

   plt.plot(ns1,Tforcl*cl*apmass1, color=clcolor)
   plt.plot(ns1,Tforcl*cl*apmass2, color=clcolor)
   plt.plot(ns1,Tforcl*cl*apmass3, color=clcolor)
   plt.plot(ns1,Tforcl*cl*apmass4, color=clcolor)
   plt.plot(ns1,Tforcl*cl*apmass5, color=clcolor, label='E precipitation')

   #apleake = en_leak = E_before - E_after - E_precip
   #apen = E_precip, so together
   #apleake + apen = E_before - E_after
   plt.plot(ns1,apleake1+apen1 - latice*apmass1,  color=rred,      label=label1, linewidth=3)
   plt.plot(ns1,apleake2+apen2 - latice*apmass2,  color='black',    label=label2, linewidth=1)
   plt.plot(ns1,apleake3+apen3 - latice*apmass3,  color=ggreen,label=label3, linewidth=2)
   plt.plot(ns1,apleake4+apen4 - latice*apmass4,  color=yyellow,   label=label4, linewidth=3)
   plt.plot(ns1,apleake5+apen5 - latice*apmass5,  color=bblue,     label=label5, linewidth=1)

   #dE for cpstarhy option
   #need relative error
   #print(apleake2+apen2)

   FS=14
   plt.legend(fontsize=FS)
   plt.title("comparing energy change with IEFLX, W/m^2",fontsize=FS)

   plt.tick_params(labelsize=FS)
   plt.grid(True)

   plt.xlabel("simulated time, sec", fontsize=FS)
   #plt.ylabel("", fontsize=FS)
   #plt.ylim([-0.0005,0.0105])

   plt.savefig(plotname+".pdf")
   plt.savefig(plotname+".png")



