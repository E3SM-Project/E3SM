import netCDF4
import sys
import numpy as np

#######################################################################################
#######################################################################################
##
## nccmp.py: A simple python-based NetCDF comparison tool that assumes the two files
## being compared have the same variables and dimensions. It computes the relative
## 2-norm, relative infinity-norm, and maximum value of the absolute difference b/t
## the two files for every variables that is of float type.
##
## python nccmp.py file1.nc file2.nc
##
#######################################################################################
#######################################################################################

#Complain if there aren't two arguments
if (len(sys.argv) < 3) :
  print("Usage: python nccmp.py file1.nc file2.nc")
  sys.exit(1)

#Open the two files
nc1 = netCDF4.Dataset(sys.argv[1])
nc2 = netCDF4.Dataset(sys.argv[2])

#Print column headers
print("Var Name".ljust(20)+":  "+"rel 2-norm".ljust(20)+"  ,  "+"rel inf-norm".ljust(20)+"  ,  "+"max abs")

#Loop through all variables
for v in nc1.variables.keys() :
  #Only compare floats
  if (nc2.variables[v].dtype == np.float64 or nc2.variables[v].dtype == np.float32) :
    #Grab the variables
    a1 = nc1.variables[v][:]
    a2 = nc2.variables[v][:]
    #Compute the absolute difference vector
    adiff = abs(a2-a1)

    #Compute the 2-norm and a normalization term
    norm2 = np.sum( adiff**2 )
    norm2_denom = np.sum( a1**2 )
    #Only apply the normalizaiton if it's non-zero
    if (norm2_denom != 0) :
      norm2 = norm2 / norm2_denom

    #Compute the inf-norm and a normalization term
    normi = np.amax( adiff )
    normi_denom = np.amax(a1) - np.amin(a1)
    #Only apply the normalizaiton if it's non-zero
    if (normi_denom != 0) :
      normi = normi / normi_denom

    #Compute the maximum absolute difference
    maxabs = np.amax( adiff )

    #Print to terminal
    print(v.ljust(20)+":  %20.10e  ,  %20.10e  ,  %20.10e"%(norm2,normi,maxabs))
