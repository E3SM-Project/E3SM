#!/usr/bin/env python2
'''
This script provides a interface between the CESM CASE environment 
and the Python package for time slice-to-series operation.
It is called from the run script and resides in the CASEROOT/Tools directory.

__________________________
Created on May 21, 2014
Update Sept 4, 2014 - only monthly stream specifiers currently working due
to a bug in the reshaper module. 

@author: CSEG <cseg@cgd.ucar.edu>
'''

# set the PYTHONPATH env var to point to the correct python locations in the 
# CASEROOT/Tools directory

import sys
import os
import glob
import re
import string
import pprint
import xml.etree.ElementTree as ET
from mpi4py import MPI

# append the CASEROOT/Tools location to the PYTHONPATH
sys.path.append('../Tools/pythonlib')

# import the cesm environment module
import cesmEnvLib

def readArchiveXML(cesmEnv):
    '''
    returns a fully defined list of reshaper specifications
    '''
    from pyreshaper import specification

    specifiers = list()
    xml_tree = ET.ElementTree()
# check if the env_archive.xml file exists
    if ( not os.path.isfile('../env_archive.xml') ):
        err_msg = "cesm_tseries_generator.py ERROR: env_archive.xml does not exists."
        raise OSError(err_msg)
    else:
# parse the xml
        xml_tree.parse('../env_archive.xml')

# loop through all the comp_archive_spec elements to find the tseries related elements
        for comp_archive_spec in xml_tree.findall("components/comp_archive_spec"):
            comp = comp_archive_spec.get("name")
            rootdir = comp_archive_spec.find("rootdir").text
            multi_instance = comp_archive_spec.find("multi_instance").text

# for now, set instance value to empty string implying only 1 instance
            instance = ""

# loop through all the files/file_spec elements
            for file_spec in comp_archive_spec.findall("files/file_extension"):
                file_extension = file_spec.get("suffix")
                subdir = file_spec.find("subdir").text

# check if tseries_create is an element for this file_spec
                if file_spec.find("tseries_create") is not None:
                    tseries_create = file_spec.find("tseries_create").text

# check if the tseries_create element is set to TRUE            
                    if tseries_create.upper() in ["T","TRUE"]:

# check if tseries_format is an element for this file_spec and if it is valid
                        if file_spec.find("tseries_output_format") is not None:
                            tseries_output_format = file_spec.find("tseries_output_format").text
                            if tseries_output_format not in ["netcdf","netcdf4","netcdf4c"]:
                                err_msg = "cesm_tseries_generator.py error: tseries_output_format invalid for data stream {0}.*.{1}".format(comp,file_extension)
                                raise TypeError(err_msg)
                        else:
                            err_msg = "cesm_tseries_generator.py error: tseries_output_format undefined for data stream {0}.*.{1}".format(comp,file_extension)
                            raise TypeError(err_msg)


# check if the tseries_output_subdir is specified and create the tseries_output_dir
                        if file_spec.find("tseries_output_subdir") is not None:
                            tseries_output_subdir = file_spec.find("tseries_output_subdir").text
                            tseries_output_dir = '/'.join( [cesmEnv["DOUT_S_ROOT"], rootdir,tseries_output_subdir] )
                            if not os.path.exists(tseries_output_dir):
                                os.makedirs(tseries_output_dir)
                        else:
                            err_msg = "cesm_tseries_generator.py error: tseries_output_subdir undefined for data stream {0}.*.{1}".format(comp,file_extension)
                            raise TypeError(err_msg)
                        

# check if tseries_tper is specified and is valid 
                        if file_spec.find("tseries_tper") is not None:
                            tseries_tper = file_spec.find("tseries_tper").text
                            if tseries_tper not in ["yearly","monthly","weekly","daily","hourly6","hourly3","hourly1","min30"]:
                                err_msg = "cesm_tseries_generator.py error: tseries_tper invalid for data stream {0}.*.{1}".format(comp,file_extension)
                                raise TypeError(err_msg)
                        else:
                            err_msg = "cesm_tseries_generator.py error: tseries_tper undefined for data stream {0}.*.{1}".format(comp,file_extension)
                            raise TypeError(err_msg)


# load the tseries_time_variant_variables into a list
                        if comp_archive_spec.find("tseries_time_variant_variables") is not None:
                            variable_list = list()
                            for variable in comp_archive_spec.findall("tseries_time_variant_variables/variable"):
                                variable_list.append(variable.text)


# get a list of all the input files for this stream from the archive location
                        history_files = list()
                        in_file_path = '/'.join( [cesmEnv["DOUT_S_ROOT"],rootdir,subdir] )                        
                        all_in_files = os.listdir(in_file_path)

# check that there are actually a list of history files to work with
                        for in_file in all_in_files:
                            if re.search(file_extension, in_file):
                                history_files.append(in_file_path+"/"+in_file)

# sort the list of input history files in order to get the output suffix from the first and last file
                        if len(history_files) > 0:
                            history_files.sort()

                            start_file = history_files[0]
                            start_file_parts = list()
                            start_file_parts = start_file.split( "." )
                            start_file_time = start_file_parts[-2]

                            last_file = history_files[-1]
                            last_file_parts = list()
                            last_file_parts = last_file.split( "." )
                            last_file_time = last_file_parts[-2]

# get the actual component name from the history file - will also need to deal with the instance numbers based on the comp_name
                            comp_name = last_file_parts[-4]
                            stream = last_file_parts[-3]

# check for pop.h nday1 and nyear1 history streams
                            if last_file_parts[-3] in ["nday1","nyear1"]:
                                comp_name = last_file_parts[-5]
                                stream = last_file_parts[-4]+"."+last_file_parts[-3]

# create the tseries output prefix needs to end with a "."
                            tseries_output_prefix = tseries_output_dir+"/"+cesmEnv["CASE"]+"."+comp_name+"."+stream+"."

# format the time series variable output suffix based on the tseries_tper setting suffix needs to start with a "."
                            if tseries_tper == "yearly":
                                tseries_output_suffix = "."+start_file_time+"-"+last_file_time+".nc"
                            elif tseries_tper == "monthly":
                                start_time_parts = start_file_time.split( "-" )
                                last_time_parts = last_file_time.split( "-" )
                                tseries_output_suffix = "."+start_time_parts[0]+start_time_parts[1]+"-"+last_time_parts[0]+last_time_parts[1]+".nc"
                            elif tseries_tper in ["weekly","daily","hourly6","hourly3","hourly1","min30"]:
                                start_time_parts = start_file_time.split( "-" )
                                last_time_parts = last_file_time.split( "-" )
                                tseries_output_suffix = "."+start_time_parts[0]+start_time_parts[1]+start_time_parts[2]+"-"+last_time_parts[0]+last_time_parts[1]+last_time_parts[2]+".nc"

# START HERE... need to create specifiers based on the tseries_filecat_years spec

# get a reshpaer specification object
                            spec = specification.create_specifier()

# populate the spec object with data for this history stream
                            spec.input_file_list = history_files
                            spec.netcdf_format = tseries_output_format
                            spec.output_file_prefix = tseries_output_prefix
                            spec.output_file_suffix = tseries_output_suffix
                            spec.time_variant_metadata = variable_list

                            dbg = list()
                            pp = pprint.PrettyPrinter(indent=5)
                            dbg = [comp_name, spec.input_file_list, spec.netcdf_format, spec.output_file_prefix, spec.output_file_suffix, spec.time_variant_metadata]
                            pp.pprint(dbg)
                            
# append this spec to the list of specifiers
                            specifiers.append(spec)

    return specifiers

#==========================
# main
#==========================

# cesmEnv is a dictionary with the case environment 
# cesmEnv["id"] = "value" parsed from the CASEROOT/env_*.xml files
cesmEnv = cesmEnvLib.readCesmXML()

# append the path to the pyreshaper package
pyReshaperPath = cesmEnv["CESMDATAROOT"] + "/tools/pyReshaper/lib/python2.7/site-packages"

sys.path.append(pyReshaperPath)

# load the pyreshaper modules
from pyreshaper import specification, reshaper

specifiers = list()
# loading the specifiers from the env_archive.xml  only needs to run on the master task (rank=0) 
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    specifiers = readArchiveXML(cesmEnv)
comm.Barrier()

# specifiers is a list of pyreshaper specification objects ready to pass to the reshaper
specifiers = comm.bcast(specifiers, root=0)

# create the PyReshaper object - uncomment when multiple specifiers is allowed
reshpr = reshaper.create_reshaper(specifiers, serial=False, verbosity=2)

# Set the output limit - the limit on the number of time-series files per processor to write.
# A limit of 0 means write all output files.
reshpr.output_limit = 0

# Run the conversion (slice-to-series) process 
reshpr.convert()

# Print timing diagnostics
reshpr.print_diagnostics()

# check if DOUT_S_SAVE_HISTORY_FILES is true or false

# Print a success statement to stdout
