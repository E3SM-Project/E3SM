#!/usr/bin/env python

import sys
import os
import re
import subprocess
import collections
import bit4bit_check
import glob


data_file_path = '/tmp/work/ab3/higher-order/reg_test/confined-shelf/data/'
bench_file_path = '/tmp/work/ab3/higher-order/reg_test/confined-shelf/bench/'

flag = bit4bit_check.bit4bit(data_file_path,bench_file_path)

dict = {}
dict['diagnostic'] = flag

if 0 in dict:
    print dict['diagnostic']




