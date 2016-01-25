#!/usr/bin/env python
import sys, os
cimeroot = os.getenv("CIMEROOT",\
        os.path.dirname(os.path.abspath(os.path.join(__file__,"..",".."))))
model = os.getenv("CIME_MODEL")
print cimeroot
libdir = os.path.join(cimeroot,"utils","python")
print libdir
sys.path.append(libdir)
from  CIME.XML.Machines import Machines
from CIME.XML.Files import Files

files = Files()

machinefile = files.get_resolved_value(files.get_value("MACHINES_SPEC_FILE"))
print machinefile




