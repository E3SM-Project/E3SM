#!/usr/bin/python

import sys
import os
import re
import inspect

## Regular expression for object file line
objline = re.compile(r"^([A-Za-z0-9_-]+)[.]o[ ]+:([^:]+)$")
modline = re.compile(r"^([A-Za-z0-9_-]+)[.]mod[ ]+:[ ]+([A-Za-z0-9_-]+)[.]o$")
modfile = re.compile(r"^([A-Za-z0-9_-]+)[.]mod$")
help_re = re.compile(r"^(-h)|(--help)$")

def Usage ():
    thisFile = os.path.realpath(inspect.getfile(inspect.currentframe()))
    print("""
Usage: {} <Depends file>

  Looks for circular dependencies in CESM Depends file

  """.format(os.path.basename(thisFile)))

def keylen(key):
    return len(key[1])

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

                print("Circular dependency: "," <== ".join(callTree[cds:] + [ dep, ]))
                badMods.append(dep)
            elif ((dep not in goodMods) and (dep not in badMods)):
                newtree = list(callTree)
                newtree.append(dep)
                goodMods, badMods = findDependencies(dependDict, newtree, goodMods, badMods)
            # No else (do nothing if we have dealt with this dependency before
        if (len(badMods) == badLen):
            # obj is a good mod
            goodMods.append(obj)
    return goodMods, badMods

def findCircularDep(filename):
    depends = {}
    modfiles = {}
    goodMods = list()
    badMods = list()
    # Read in all the object file dependencies
    with open(filename) as file:
        for line in file:
            lmatch = objline.match(line.strip())
            if lmatch is not None:
                obj = lmatch.group(1)
                deps = [ x.split(".")[0].lower() for x in lmatch.group(2).strip().split(" ") if "." in x and x.split(".")[1] == "mod" ]
                depends[obj.strip()] = deps

            lmatch = modline.match(line.strip())
            if lmatch is not None:
                mod = lmatch.group(1)
                obj = lmatch.group(2)
                modfiles[mod.strip().lower()] = obj.strip()

    # Fix up the depends to use object names instead of module names
    for key in depends:
        for loc in range(len(depends[key])):
            mfile = depends[key][loc]
            if mfile in modfiles:
                depends[key][loc] = modfiles[mfile]

    # Now, iterate through object files looking for circle (short lists first)
    for obj in [ x[0] for x in sorted(depends.items(), key=keylen) ]:
        goodMods, badMods = findDependencies(depends, ( obj, ), goodMods, badMods)
    return (len(badMods) == 0)

def main(filename):
    DepFile = os.path.abspath(filename)
    if (not os.path.exists(DepFile)):
        print("ERROR: File '{0}', does not exist".format(filename))
        return 1
    cleanTree = findCircularDep(DepFile)
    if (cleanTree):
        print("No circular dependencies found")
        return 0
    else:
        return -1

if __name__ == "__main__":
    if len(sys.argv) == 2:
        if help_re.match(sys.argv[1]) is not None:
            Usage()
            sys.exit(0)
        else:
            retcode = main(sys.argv[1])
            exit(retcode)
    else:
        Usage()
        sys.exit(1)
