#!/usr/bin/env python2
'''
This module provides functions for parsing all the $CASEROOT/env_*.xml files
except the env_archive.xml file which is parsed in the cesm_tseries_generate.py
This script lives in the CASEROOT/Tools directory after create_newcase is run.
__________________________
Created on Apr 30, 2014
Updated Sept 4, 2014 - make sure execute permission is allowed

@author: NCAR - CSEG
'''
import os
import re
import xml.etree.ElementTree as ET

re_val = re.compile(r'\$(\{([A-Za-z0-9_]+)\}|[A-Za-z0-9_]+)')

def expand (val, src):
  '''
  recursively traverse the environment dictionary to expand all environment variables
  '''
  return re_val.sub(lambda m: expand(src.get(m.group(1), ''), src), val)

def readCesmXML():
  '''
  for the <entry id=value> xml files
  returns a dictionary output["id"]="value"
  '''
  output = os.environ.copy()
  env_file_list = ["../env_case.xml","../env_run.xml","../env_build.xml","../env_mach_pes.xml"]
  for env_file in env_file_list:
    if os.path.isfile(env_file):
      xml_tree = ET.ElementTree()
      xml_tree.parse(env_file)
      for entry_tag in xml_tree.findall("entry"):
        output[entry_tag.get("id")] = entry_tag.get("value")
    else:
      err_msg = "cesmEnvLib.py ERROR: {0} does not exist or cannot be read".format(env_file)
      raise TypeError(err_msg)

  for k,v in output.iteritems():
    output[k] = expand(v, output)

  return output


