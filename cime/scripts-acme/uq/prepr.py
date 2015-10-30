#!/usr/bin/env python

import csv
import sys
try:
    import numpy as np
except ImportError:
    print "Numpy was not found. "


filename=sys.argv[1]
nsam=int(sys.argv[2])
remove_pars=sys.argv[3:]

print "Reading parameter file       : ", filename


pnames_clm_file='pnames_clm.dat'
input_clm_file='input_clm.dat'
domain_x_file='pdomain_x.dat'
pnames_x_file='pnames_x.dat'
input_x_file='input_x.dat'

f_clm = open(pnames_clm_file,'w')
f_x = open(pnames_x_file,'w')

pdomain_x_list=[]
remove_par_ind=[]
with open(filename, 'rU') as csvfile:
    linereader = csv.reader(csvfile, delimiter=',')
    next(linereader, None)
    parnum=1
    for row in linereader:
        pmin=float(row[0])
        pmax=float(row[1])
        if (pmin>=pmax):
            print "Error in reading input domain bounds."
            sys.exit()
        use=int(row[3])
        name=row[4]
        if (use==1):
            f_clm.write(name+'\n')
            if (name not in remove_pars):
                pdomain_x_list.append([pmin,pmax])
                f_x.write(name+'\n')
            else:
                remove_par_ind.append(parnum)
            parnum+=1
f_clm.close()
f_x.close()

pdomain_x=np.array(pdomain_x_list)
np.savetxt(domain_x_file,pdomain_x)
print "Written parameter names file (CLM)   : ", pnames_clm_file
print "Written domain file          (comp.) : ", domain_x_file
print "Written parameter names file (comp.) : ", pnames_x_file

npar=len(pdomain_x_list)
input_x=np.random.rand(nsam,npar)
print "Generated %d uniform random i.i.d. samples for %d parameters" % (nsam,npar)
np.savetxt(input_x_file,2.*input_x-1.)
print "Written computational inputs to ", input_x_file

input_clm=np.dot(input_x, np.diag(pdomain_x[:,1]-pdomain_x[:,0]))+np.dot(np.ones((nsam,npar)),np.diag(pdomain_x[:,0]))

for ind in remove_par_ind:
    tmp=np.insert(input_clm, ind-1, 1-input_clm[:,ind-2]-input_clm[:,ind-3], axis=1)
    input_clm=tmp.copy()

np.savetxt(input_clm_file,input_clm)
print "Written CLM inputs to           ", input_clm_file

#print remove_par_ind

