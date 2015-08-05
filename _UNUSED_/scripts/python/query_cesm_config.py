#!/usr/bin/env python
""" Query the CESM config_compilers file for compiler settings for cesm supported machines

Author: Jim Edwards <jedwards@ucar.edu>

"""

from __future__ import print_function
import sys
if sys.hexversion < 0x02070000:
    print(70 * "*")
    print("ERROR: {0} requires python >= 2.7.x. ".format(sys.argv[0]))
    print("It appears that you are running python {0}".format(
        ".".join(str(x) for x in sys.version_info[0:3])))
    print(70 * "*")
    sys.exit(1)

import os
import traceback


lib_path = os.path.join('scripts/python/contrib/unit_testing')
sys.path.append(lib_path)
from machine_setup import MachineCompilerSettings
from machine_setup import get_machine_name



# -------------------------------------------------------------------------------
#
# User input
#
# -------------------------------------------------------------------------------

# def commandline_options():
#     """Process the command line arguments.

#     """
#     parser = argparse.ArgumentParser(
#         description='FIXME: python program template.')



#     options = parser.parse_args()
#     return options


# # -------------------------------------------------------------------------------
# #
# # main
# #
# # -------------------------------------------------------------------------------

# def main(options):
#     machinefilename = options.xmlpath[0]+"/config_machines.xml"
#     compilerfilename = options.xmlpath[0]+"/config_compilers.xml"
#     if options.mach:
#         machine_name = options.mach[0] 
#     else:
#         machine_name = get_machine_name()

#     compiler_name = options.compiler[0]

#     print("Reading machines xml file : {0}".format(machinefilename))

#     xmlfile = os.path.abspath(machinefilename)
#     if not os.path.isfile(xmlfile):
#         raise RuntimeError("Could not find machines file: {0}".format(xmlfile))

# #    print("searching machine = {0} for field {1}".format(machine_name, field))
    



#     return 0

    


# if __name__ == "__main__":
#     options = commandline_options()
#     try:
#         status = main(options)
#         sys.exit(status)
#     except Exception as error:
#         print(str(error))
#         if options.backtrace:
#             traceback.print_exc()
#         sys.exit(1)
