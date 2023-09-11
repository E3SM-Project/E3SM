#!/usr/bin/env python
from string import *
import os, commands, getopt, sys, exceptions
import numpy as np
import matplotlib.pyplot as plt


#call is as
#python parse-precip.py ../output-cpstarnh-l128-dt0.2 ../output-cpstarnh-l128-dt0.1 ../output-cpstarnh-l128-dt0.05 ../output-cpstarnh-l128-dt0.02 ../output-cpstarnh-l128-dt0.01  ../output-cpstarnh-l128-dt0.005    dt0.2 dt0.1 dt0.05 dt0.02 dt0.01 dt0.005 NAME PATTERN TITLE


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
  s3=arr4.size

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


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
print 'file1 is', str(sys.argv[1])
print 'file2 is', str(sys.argv[2])
print 'file3 is', str(sys.argv[3])
print 'file4 is', str(sys.argv[4])
print 'file5 is', str(sys.argv[5])
print 'file6 is', str(sys.argv[6])
print 'LABEL1: 7 is', str(sys.argv[7])
print 'LABEL2: 8 is', str(sys.argv[8])
print 'LABEL3: 9 is', str(sys.argv[9])
print 'LABEL4: 10 is', str(sys.argv[10])
print 'LABEL5: 11 is', str(sys.argv[11])
print 'LABEL6: 12 is', str(sys.argv[12])
print 'NAME OF OUTPUT PLOT: 13 is', str(sys.argv[13])
print 'PATTERN: 14 is', str(sys.argv[14])
print 'title: 15 is', str(sys.argv[15])

label1=str(sys.argv[7])
label2=str(sys.argv[8])
label3=str(sys.argv[9])
label4=str(sys.argv[10])
label5=str(sys.argv[11])
label6=str(sys.argv[12])
plotname=str(sys.argv[13])
pattern=str(sys.argv[14])
ttitle=str(sys.argv[15])


##############example from the file
#  ptt=   556.000000000000       max_prcl =   5187.12019299356
#  pts=   556.000000000000       sum_prcl =   821764080.619827
#  cls=   556.000000000000       ecl_diff =   2.15049951069576



stri1=pattern
 
dwater1=[]
dwater2=[]
dwater3=[]
dwater4=[]
dwater5=[]
dwater6=[]
nstep1=[]
nstep2=[]
nstep3=[]
nstep4=[]
nstep5=[]
nstep6=[]

try:
   kkk=open(str(sys.argv[1]),"r")
   lll=open(str(sys.argv[2]),"r")
   mmm=open(str(sys.argv[3]),"r")
   nnn=open(str(sys.argv[4]),"r")
   nnn2=open(str(sys.argv[5]),"r")
   nnn3=open(str(sys.argv[6]),"r")
finally:

   ff=kkk #1st file

   #look for water
   noteof=1
   while noteof :

      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        #print(strr[0])  #0
    	#print(strr[1]) #double
    	#print(strr[2]) #double
    	#print(strr[3]) #double
        nstep1.extend([atof(strr[1])])
        dwater1.extend([atof(strr[2])])
      else :
        noteof=0
   print "finished readfing the 1st file"

#end of the 1st file, read the second:
   ff=lll #1st file

   #look for water
   noteof=1
   while noteof :

      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        #print(str[0]) #0
        #print(str[1]) #double
        nstep2.extend([atof(strr[1])])
        dwater2.extend([atof(strr[2])])
      else :
        noteof=0

   print "finished readfing the 2nd file"

   #3rd file
   ff=mmm #1st file

   #look for water
   noteof=1
   while noteof :

      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        #print(str[0]) #0
        #print(str[1]) #double
        nstep3.extend([atof(strr[1])])
        dwater3.extend([atof(strr[2])])
      else :
        noteof=0

   print "finished readfing the 3rd file"

   #4 file
   ff=nnn #1st file

   #look for water
   noteof=1
   while noteof :

      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        #print(str[0]) #0
        #print(str[1]) #double
        nstep4.extend([atof(strr[1])])
        dwater4.extend([atof(strr[2])])
      else :
        noteof=0

   print "finished readfing the 4th file"

   #5 file
   ff=nnn2 #1st file

   #look for water
   noteof=1
   while noteof :

      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        #print(str[0]) #0
        #print(str[1]) #double
        nstep5.extend([atof(strr[1])])
        dwater5.extend([atof(strr[2])])
      else :
        noteof=0

   print "finished readfing the 5th file"

   #6 file
   ff=nnn3 #1st file

   #look for water
   noteof=1
   while noteof :

      strr=lookfor1(ff,stri1,"",1)
      if (strr <> "EOF") :
        strr=split(strr)
        #print(str[0]) #0
        #print(str[1]) #double
        nstep6.extend([atof(strr[1])])
        dwater6.extend([atof(strr[2])])
      else :
        noteof=0

   print "finished readfing the 6th file"






   ####################################################

   #f1=(np.array(te_before_dme1)-np.array(te_after_dycore1))/dtime
   ns1=np.array(nstep1)
   nn1=ns1.size
   print "size1="+str(nn1)

   ns2=np.array(nstep2)
   nn2=ns2.size
   print "size2="+str(nn2)
   ns3=np.array(nstep3)
   nn3=ns3.size
   print "size3="+str(nn3)
   ns4=np.array(nstep4)
   nn4=ns4.size
   print "size4="+str(nn4)
   ns5=np.array(nstep5)
   nn5=ns5.size
   print "size5="+str(nn5)

   ns6=np.array(nstep6)
   nn6=ns6.size
   print "size6="+str(nn6)


   aw1=np.array(dwater1)
   aw2=np.array(dwater2)
   aw3=np.array(dwater3)
   aw4=np.array(dwater4)
   aw5=np.array(dwater5)
   aw6=np.array(dwater6)

   print "rr1",aw1[0:400]


   #print "First log, fixer: avg and std",sum(fixer1)/len(fixer1),np.std(fixer1) 
   #print "second log, fixer: avg and std",sum(fixer2)/len(fixer2),np.std(fixer2) 


######################################################################################

#   plot4arrays(arr1, arr2, arr3, arr4, label1, label2, label3, label4, "RR-"+plotname, 10, "energy imbalance")
#   plot4arrays(aww1, aww2, aww3, aww4, label1, label2, label3, label4, "WATER-"+plotname, 10, "water imbalance")


   plt.figure()

   #legend0=(arr1name+": avg %g , std: %g" % (sum(b1)/len(b1),np.std(b1)) )
   #legend1=(arr2name+": avg %g , std: %g" % (sum(b2)/len(b2),np.std(b2)) )
   #legend2=(arr3name+": avg %g , std: %g" % (sum(b3)/len(b3),np.std(b3)) )
   #legend3=(arr4name+": avg %g , std: %g" % (sum(b4)/len(b4),np.std(b4)) )
   #print(log(abs(b1)))

   bblue='#377eb8'
   oorange='#ff7f00'
   ggreen='#4daf4a'
   ppink='#f781bf'
   bbrown='#a65628'
   ppurple='#984ea3'
   ggray='#999999'
   rred='#e41a1c'
   yyellow='#dede00'

   plt.plot(ns1,aw1,  color=ggreen,label=label1)
   plt.plot(ns2,aw2,  color=bblue,label=label2)
   plt.plot(ns3,aw3,  color=bbrown,label=label3)
   plt.plot(ns4,aw4,  color=rred,label=label4)
   plt.plot(ns5,aw5,  color=oorange,label=label5)
   plt.plot(ns6,aw6,  color='black',label=label6)

   FS=14
   plt.legend(fontsize=FS)

   if ( ttitle == "precl" ) :
     plt.title("precipitation, kg/m^2/sec",fontsize=FS)
   if( ttitle=="encl" ) :
     plt.title("accum. diff. of precipitation energy",fontsize=FS)
   
   plt.tick_params(labelsize=FS)
   plt.grid(True)

   plt.xlabel("simulated time, sec", fontsize=FS)
   #plt.ylabel("", fontsize=FS)
   plt.ylim([-0.0005,0.0105])

   plt.savefig(plotname+".pdf")
   plt.savefig(plotname+".png")



#   print(abs(arr1))
#   print(abs(arr2))
#   print(abs(arr3))
#   print(abs(arr4))
   #print(aww1/arr1)



