#!/usr/bin/python

import sys
import os
import re
import inspect

## Regular expression for object file line
objline = re.compile(r"^([A-Za-z0-9_-]+)[.]o[ ]+:([^:]+)$")
modfile = re.compile(r"^([A-Za-z0-9_-]+)[.]mod$")

def Usage ():
  thisFile = os.path.realpath(inspect.getfile(inspect.currentframe()))
  print """
Usage: %s <Depends file>

  Looks for circular dependencies in CESM Depends file

  """%(os.path.basename(thisFile))
# End Usage

def keylen(key):
  return len(key[1])
# End def

def findDependencies(dependDict, callTree, goodMods, badMods):
  obj = callTree[-1]
  badLen = len(badMods)
  if ((obj in dependDict) and (obj not in goodMods)):
    for dep in dependDict[obj]:
      if (dep in callTree):
        # The circle might not start at the beginning of the call tree
        cds = 0
        while ((cds < len(callTree)) and (callTree[cds] != dep)):
          cds = cds + 1
        # End while
        print "Circular dependency: "," <== ".join(callTree[cds:] + [ dep, ])
        badMods.append(dep)
      elif ((dep not in goodMods) and (dep not in badMods)):
        newtree = list(callTree)
        newtree.append(dep)
        goodMods, badMods = findDependencies(dependDict, newtree, goodMods, badMods)
      # No else (do nothing if we have dealt with this dependency before
      # End if
    # End for
    if (len(badMods) == badLen):
      # obj is a good mod
      goodMods.append(obj)
    # End if
  # End if
  return goodMods, badMods
# End def

def findCircularDep(filename):
  depends = {}
  goodMods = list()
  badMods = list()
  # Read in all the object file dependencies
  with open(filename) as file:
    for line in file:
      lmatch = objline.match(line.strip())
      if (lmatch is not None):
        obj = lmatch.group(1)
        deps = [ x.split(".")[0].lower() for x in lmatch.group(2).strip().split(" ") if x.split(".")[1] == "mod" ]
        depends[obj.strip().lower()] = deps
      # End if
    # End for
  # End with
  # Now, iterate through object files looking for circle (do short lists first)
  for obj in [ x[0] for x in sorted(depends.items(), key=keylen) ]:
    goodMods, badMods = findDependencies(depends, ( obj, ), goodMods, badMods)
  # End for
  return (len(badMods) == 0)
# End def

def main(filename):
  DepFile = os.path.abspath(filename)
  if (not os.path.exists(DepFile)):
    print "ERROR: File '%s', does not exist"%filename
    return 1
  # End if
  cleanTree = findCircularDep(DepFile)
  if (cleanTree):
    print "No circular dependencies found"
    return 0
  else:
    return -1
  # End if
# End def

if __name__ == "__main__":
  if len(sys.argv) == 2:
    retcode = main(sys.argv[1])
    exit(retcode)
  else:
    Usage()
    sys.exit(1)
  # End if
# End if
