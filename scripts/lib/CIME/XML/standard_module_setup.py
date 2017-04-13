
# pragma pylint: disable=unused-import

import logging
import os
import sys
import xml.etree.ElementTree as ET
import re
LIB_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(LIB_DIR)
from CIME.utils import expect, run_cmd, run_cmd_no_fail, get_cime_root
