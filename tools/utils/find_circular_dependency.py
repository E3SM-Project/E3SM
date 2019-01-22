#!/usr/bin/python
"""Look for and output any circular dependencies in a CESM Depends file"""

import sys
import os
import re
import inspect

## Regular expression for object file line
_OBJLINE = re.compile(r"^([A-Za-z0-9_-]+)[.]o[ ]+:([^:]+)$")
_MODLINE = re.compile(r"^([A-Za-z0-9_-]+)[.]mod[ ]+:[ ]+([A-Za-z0-9_-]+)[.]o$")
_HELP_RE = re.compile(r"^(-h)|(--help)$")

def usage():
    """Looks for circular dependencies in CESM Depends file"""
    this_file = os.path.realpath(inspect.getfile(inspect.currentframe()))
    print("""
Usage: {} <Depends file>

  {}

    """.format(os.path.basename(this_file), usage.__doc__))

def keylen(key):
    """Return the length of the input, <key>"""
    return len(key[1])

def find_dependencies(depend_dict, call_tree, good_mods, bad_mods):
    """Find circular dependencies in depend_dict"""
    obj = call_tree[-1]
    bad_len = len(bad_mods)
    if (obj in depend_dict) and (obj not in good_mods):
        for dep in depend_dict[obj]:
            if dep in call_tree:
                # The circle might not start at the beginning of the call tree
                cds = 0
                while (cds < len(call_tree)) and (call_tree[cds] != dep):
                    cds = cds + 1

                print("Circular dependency: {}".format(" <== ".join(call_tree[cds:] + [dep, ])))
                bad_mods.append(dep)
            elif (dep not in good_mods) and (dep not in bad_mods):
                newtree = list(call_tree)
                newtree.append(dep)
                good_mods, bad_mods = find_dependencies(depend_dict, newtree, good_mods, bad_mods)
            # No else (do nothing if we have dealt with this dependency before
        if len(bad_mods) == bad_len:
            # obj is a good mod
            good_mods.append(obj)
    return good_mods, bad_mods

def find_circular_dep(filename):
    """Parse <filename>, then check for circular dependencies"""
    depends = {}
    modfiles = {}
    good_mods = list()
    bad_mods = list()
    # Read in all the object file dependencies
    with open(filename) as fh1:
        for line in fh1:
            lmatch = _OBJLINE.match(line.strip())
            if lmatch is not None:
                obj = lmatch.group(1)
                deps = [x.split(".")[0].lower() for x in
                        lmatch.group(2).strip().split(" ") if "." in x and x.split(".")[1] == "mod"]
                depends[obj.strip()] = deps

            lmatch = _MODLINE.match(line.strip())
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
    for obj in [x[0] for x in sorted(depends.items(), key=keylen)]:
        good_mods, bad_mods = find_dependencies(depends, (obj, ), good_mods, bad_mods)
    return len(bad_mods) == 0

def main(filename):
    """Check <filename>, then examine it for circular dependencies"""
    dep_file = os.path.abspath(filename)
    if not os.path.exists(dep_file):
        print("ERROR: File '{0}', does not exist".format(filename))
        return 1
    clean_tree = find_circular_dep(dep_file)
    if clean_tree:
        print("No circular dependencies found")
        return 0

    return -1

if __name__ == "__main__":
    if len(sys.argv) == 2:
        if _HELP_RE.match(sys.argv[1]) is not None:
            usage()
            sys.exit(0)
        else:
            exit(main(sys.argv[1]))
    else:
        usage()
        sys.exit(1)
