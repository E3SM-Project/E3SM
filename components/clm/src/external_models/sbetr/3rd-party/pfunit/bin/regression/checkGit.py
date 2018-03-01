#!/usr/bin/python
import commands
import sys
import os
import subprocess

# This script manages the pFUnit regresion tests submission by checking
# whether selected branches have been updated or not.

def syscmd(cmd, stdout=subprocess.PIPE):
     p = subprocess.Popen( cmd , shell = True, \
         stdout = stdout, stderr = subprocess.PIPE )
     (out, err) = p.communicate()
     return out

def main():
   branches = {'master'      : 'NO', \
               'release-3.0' : 'NO', \
               'development' : 'NO', \
               'pfunit_2.1.0': 'NO'};

   print "Check for repository updates...";
   # We compare the SHA keys between reference repository, assumed
   # to be in the users's $HOME, and pFUnit's central repository.
   home = os.environ['HOME'];
   for branch, dotest in branches.iteritems():
     os.chdir(home + '/pFUnit.git/' + branch);
     befKey = syscmd('git rev-parse HEAD')
     pull = syscmd('git pull') # update repository for "next" time
     aftKey = syscmd('git rev-parse HEAD')
     if  (str(befKey) != str(aftKey)):
        print " -- Changes detected in branch " + branch;
        branches[branch] = 'YES';
     else:
        print " -- NO changes detected in branch " + branch;

   for branch, dotest in branches.iteritems():
     if (dotest == 'YES'):
        print " -- Running tests for branch " + branch;
        with open(home + '/bin/' + branch + '.log', 'w') as log:
           rc = syscmd(home + '/bin/mainRegress.sh ' + branch, \
                       stdout=log)
     else:
        print " -- Nothing to do for branch " + branch;
      
if __name__ == "__main__":
    main()

# TODO: Need to parallelize branch tests. Right now branch tests are run
#       sequentially.
