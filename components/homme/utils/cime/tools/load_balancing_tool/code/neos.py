#!/usr/bin/env python
#########################################

# NeosClient.py
#########################################
import sys
import xmlrpclib
import time

NEOS_HOST="neos-server.org"
NEOS_PORT=3332

if len(sys.argv) != 2:
  sys.stderr.write("Usage: NeosClient <xmlfilename | queue> ")
  sys.exit(1)

neos=xmlrpclib.Server("http://%s:%d" % (NEOS_HOST, NEOS_PORT))

if sys.argv[1] == "queue":
  #Print NEOS job queue
  msg = neos.printQueue()
  sys.stdout.write(msg)
else:
  #Read XML file
  xmlfile = open(sys.argv[1],"r")
  xml=""
  buffer=1

  while buffer:
    buffer =  xmlfile.read()
    xml+= buffer
  xmlfile.close()

  (jobNumber,password) = neos.submitJob(xml)
  sys.stdout.write("JobNumber = %d " % jobNumber)

  offset=0

  status=""
  #Print out partial job output while job is running
  while status != "Done":
    (msg,offset) = neos.getIntermediateResults(jobNumber,password,offset)
    sys.stdout.write(msg.data)
    status = neos.getJobStatus(jobNumber, password)

  #Print out the final result
  msg = neos.getFinalResults(jobNumber, password).data

  sys.stdout.write(msg)

