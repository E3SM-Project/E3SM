#!/usr/bin/env python
# For python 2.6-2.7
from __future__ import print_function
# For python2.5
# from __future__ import with_statement

import re
import sys

includeFile="generated.inc"
cmakeIncludeFile="cmakeGenerated.inc"

generatedFiles = []
with open(includeFile,'r') as f:
   l0 = f.readline()
   lines = f.readlines()
   generatedFiles = [(re.match(r'(.*)\+\= (.*)',line)).group(2) for line in lines]
#   for line in lines:
#      obj = re.match(r'(.*)\+\= (.*)',line)
#   print (' : "'+obj.group()+'"')
#   print ('1: "'+obj.group(2)+'"')
# print (generatedFiles)
#+++ print (''+' '.join(generatedFiles))
# print ('\n'+'\n'.join(generatedFiles))

with open(cmakeIncludeFile,'w') as f:
    for i in generatedFiles:
        f.write('list(APPEND srcs '+i+' )\n')


# print (cmakeIncludeFile)
sys.stdout.write(cmakeIncludeFile)






    
    
