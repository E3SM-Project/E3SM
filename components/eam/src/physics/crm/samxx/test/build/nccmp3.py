import netCDF4
import sys
import numpy as np

#######################################################################################
#######################################################################################
##
## nccmp3.py: A simple python-based NetCDF 3-way comparison tool. The purpose of this
## is to show whether files 2 and 3 have more of a difference than files 1 and 2. The
## specific purpose is to compare refactored differences against presumed bit-level
## differences (like -O0 and -O3 compiler flags). This prints the relative 2-norm of
## the absolute differences between files 1 & 2 and files 2 & 3, as well as the ratio
## of the relative 2-norms between the 2-3 comparison and the 1-2 comparison.
##
## python nccmp.py file1.nc file2.nc file3.nc
##
#######################################################################################
#######################################################################################

#Complain if there aren't two arguments
if (len(sys.argv) < 4) :
  print("Usage: python nccmp.py file1.nc file2.nc")
  sys.exit(1)

#Open the files
nc1 = netCDF4.Dataset(sys.argv[1])
nc2 = netCDF4.Dataset(sys.argv[2])
nc3 = netCDF4.Dataset(sys.argv[3])

#Print column headers
print("Var Name".ljust(20)+":  "+"|1-2|".ljust(20)+"  ,  "+"|2-3|".ljust(20)+"  ,  "+"|2-3|/|1-2|")

#Loop through all variables
for v in nc1.variables.keys() :
  #Only compare floats
  if (nc2.variables[v].dtype == np.float64 or nc2.variables[v].dtype == np.float32) :
    #Grab the variables
    a1 = nc1.variables[v][:]
    a2 = nc2.variables[v][:]
    a3 = nc3.variables[v][:]
    #Compute the absolute difference vectors
    adiff12 = abs(a2-a1)
    adiff23 = abs(a3-a2)

    #Compute relative 2-norm between files 1 & 2 and files 2 & 3
    norm12 = np.sum( adiff12**2 )
    norm23 = np.sum( adiff23**2 )
    #Assume file 1 is "truth" for the normalization
    norm_denom = np.sum( a1**2 )
    #Only normalize if this denominator is != 0
    if (norm_denom != 0) :
      norm12 = norm12 / norm_denom
      norm23 = norm23 / norm_denom

    #Compute the ratio between the 2-3 norm and the 1-2 norm
    normRatio = norm23
    #If the denom is != 0, then go ahead and compute the ratio
    if (norm12 != 0) :
      normRatio = norm23 / norm12
    else :
      #If they're both zero, then just give a ratio of "1", showing they are the same
      if (norm23 == 0) :
        normRatio = 1
      #If the 1-2 norm is zero but the 2-3 norm is not, give a very large number so the user is informed
      else :
        normRatio = 1e50

    #Only print ratios that are > 2, meaning 2-3 diff is >2x more than the 1-2 diff.
    #In the future, this should be added as a command line parameter for the user to choose.
    if (normRatio > 2) :
      print(v.ljust(20)+":  %20.10e  ,  %20.10e  ,  %20.10e"%(norm12,norm23,norm23/norm12))
