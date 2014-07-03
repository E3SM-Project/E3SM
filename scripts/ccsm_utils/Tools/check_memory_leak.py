#! /usr/bin/env python 
#
# script name: check_memory_leak.py 
# Usage: check_memory_leak.py filename filename_baseline (optional)  (where filename* is a coupler log, e.g. cpl.log*,  produced by cesm) 
#
# purpose:
# Detect possible memory leaks by comparing memory highwater marks recorded in the coupler 
# log (cpl.log*) for consecutive model days. If memory usage continually increases by more than 2MB/day, 
# consider it an indication of a memory leak and set test status to FAIL.  
# If test is being compared to to baseline results ( e.g. -compare cesm1_1_X ) then compare "pes max memory" 
# high water mark with current run and determine whether is has increased by more than 5MB since baseline.   
# If current run "pes max memory" highwater value is more than 5MB greater than baseline cpl.log "pes max memory"   
# highwater value consider it a failure and set test status to FAIL.
  
import sys
import re
memleak = 0 
memIncreaseFromBaseline = 0
memHighWater = []
pesmemHighWater = []
modelDate = [] 

if len(sys.argv) < 2: 
    print 'Usage: memory_leak_check cpl.log.name baseline.cpl.log'
    sys.exit(1) 

# read cpl.log*  into a variable 
fileContents  = file(sys.argv[1]).read() 

# open cpl.log* file
infile = open( sys.argv[1] )

# find pes max memory highwater in MB in cpl.log* and baseline cpl.log
if len(sys.argv) == 3: 
   fileContentsBaseline = file(sys.argv[2]).read()
   regex = re.compile(r'pes max memory highwater.*MB.(.{8})')
   matchObj1 = re.search( regex, str(fileContentsBaseline))
   matchObj2 = re.search( regex, str(fileContents))
   if matchObj1: 
       pesMaxMemHighWaterBaseline = matchObj1.group(1)
   if matchObj2: 
       pesMaxMemHighWater = matchObj2.group(1)


# compare pes memory high water mark with previous baseline if >5MB increase FAIL  
if len(sys.argv) == 3:
   if (float(pesMaxMemHighWater) - float(pesMaxMemHighWaterBaseline) > 5 ): 
      print "memory highwater mark per tasks has increased by more than 5MB since last baseline"  
      memIncreaseFromBaseline = 1 

# Find lines containing memory highwater mark values per model day in MBs, store model day and 
# corresponding memory highwater mark.  Find maximum memory high water mark per pe.    
for line in infile:
    matchObj3 = re.search( r'model date = (.{8}).*memory =(.*)MB.*highwater', line)
    if matchObj3: 
       modelDate.append( matchObj3.group(1) )
       memHighWater.append( matchObj3.group(2) )

# for consecutive model days determine whether memory usage has increased by more than 2MB
# every day 
for i in range ( len(memHighWater) -1 ): 
    if ( (int(modelDate[i+1]) - int(modelDate[i]) )  == 1):        
       diff = float(memHighWater[i+1]) - float(memHighWater[i]) 
       if (diff  >= 2): 
           memleak = 1  
       else: 
           memleak = 0
           break

if (memleak == 1): 
    print "memory leak. memory increased by 2 MB/day or more" 
if ( (memIncreaseFromBaseline == 1) | (memleak ==1) ): 
    print "FAIL"
elif ( memleak == 0): 
    print "not a memory leak" 
    print "PASS"                   
infile.close() 
