import logging
import os
import sys
#import libxml2
#import xml.etree.ElementTree as ET
#import lxml
import lxml.etree as ET
import re
LIB_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(LIB_DIR)
