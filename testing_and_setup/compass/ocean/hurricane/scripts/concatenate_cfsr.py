# Author: Steven Brus
# Date: August 2019
# Description: This script concatenates a series of CFSR data files into a single file

import subprocess
import os
import glob

pwd = os.getcwd()+'/'
 
filetypes = ['prmsl','wnd10m']

for ft in filetypes:

  files = sorted(glob.glob(pwd+ft+'*.nc'))
  print(files) 

  filenames = [f.split('/')[-1] for f in files]
  print(filenames)

  command = ['ncrcat']
  command.extend(filenames)
  command.append(ft+'.nc')
  print(command)
  subprocess.call(command)
