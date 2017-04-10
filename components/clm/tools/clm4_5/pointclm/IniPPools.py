#!/usr/bin/env python

from optparse import OptionParser
import numpy as np
import os, shutil
from netcdf_functions import putvar
import glob

parser = OptionParser()
parser.add_option("--casesite", dest="casesite", default="", \
                  help="case plus site name?")
parser.add_option("--casename", dest="casename", default="", \
                  help="case name?")
parser.add_option("--diricase", dest="diricase", default="", \
                  help="case input directory")
parser.add_option("--dirocase", dest="dirocase", default="", \
                  help="case output directory")
parser.add_option("--sitephos", dest="sitephos", default="", \
                  help="site P pool data")
parser.add_option('--restart_year', dest='restart_year', default="", \
                    help='Year of restart file to modify')

(options, args) = parser.parse_args()

# the directory of input
casesite = options.casesite
casename = options.casename
diricase = options.diricase
if (options.dirocase == ''):
    dirocase = options.diricase
else:
    dirocase = options.dirocase
sitephos = options.sitephos



if (options.restart_year == ''):
        #if restart_year not provided, take the last existing restart file
        restart_file = glob.glob(diricase+'/'+casename+'.clm2.r.*.nc')
        restart_file_last = restart_file[-1]
        year = int(restart_file_last[-19:-15])
else:
        year = int(options.restart_year)


#site = casesite[len(casesite)-6:len(casesite)]

site = casesite

solutionP = {}
labileP   = {}
secondP   = {}
occlP     = {}
primP     = {}

with open(sitephos) as f:
     lines = f.readlines()
     contents = [x.rstrip('\n') for x in lines]
     for content in contents[1:]:
         print not content.strip()
         if content.strip():
            sl = content.split(',')
            sitename = sl[0]
            solutionP[sitename]=float(sl[1])
            labileP  [sitename]=float(sl[2])
            secondP  [sitename]=float(sl[3])
            occlP    [sitename]=float(sl[4])
            primP    [sitename]=float(sl[5])
            print sitename, solutionP[sitename], labileP  [sitename], secondP  [sitename], occlP    [sitename], primP    [sitename]
f.close()

if site in solutionP.keys():
   fileinp = diricase + '/' + casename + ".clm2.r."+ str(10000+year)[1:] + "-01-01-00000.nc"
   fileout = dirocase + '/' + casename + ".clm2.r."+ str(10000+year)[1:] + "-01-01-00000.nc"

   if diricase != dirocase:
      shutil.copy(fileinp, dirocase)
   shutil.copy(fileinp, fileinp+".b")

   putvar(fileout, 'solutionp_vr', solutionP[site])
   putvar(fileout, 'labilep_vr'  , labileP  [site])
   putvar(fileout, 'secondp_vr'  , secondP  [site])
   putvar(fileout, 'occlp_vr'    , occlP    [site])
   putvar(fileout, 'primp_vr'    , primP    [site])

      
