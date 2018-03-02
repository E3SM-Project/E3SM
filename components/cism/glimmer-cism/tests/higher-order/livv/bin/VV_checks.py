#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import glob

#bit4bit check
def bit4bit(model_file_path,bench_file_path):
    a = []
    b = []
    flag = 0

    model_file_list = glob.glob(model_file_path + '/*.nc')

    bench_file_list = glob.glob(bench_file_path + '/*.nc')

    for model_file in model_file_list:
            md5sum_model = "md5sum " + model_file
            p = subprocess.Popen(md5sum_model, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            output = p.communicate()[0]
            a.append(output[0:32])

    for bench_file in bench_file_list:
            md5sum_bench = "md5sum " + bench_file
            q = subprocess.Popen(md5sum_bench, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            output = q.communicate()[0]
            b.append(output[0:32])

#match up the md5sums in a and b, tell if bit-for-bit for each specific case
    for b_file in b:
        if b_file not in a:
            i = b.index(b_file)
            flag = 1
            break

    return flag

#ncdump_model = "ncdump -c " + model_file
#ncdump_bench = "ncdump -c " + bench_file

#p = subprocess.Popen(ncdump_model, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
#output_model = p.communicate()

#q = subprocess.Popen(ncdump_bench, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
#output_bench = q.communicate()

#if output_model == output_bench:
#       print "bit-for-bit"
#else:
#       print "not bit-for-bit: in terms of ncdump -c"

##need to write results to GIS-con-diag.html
#make red if not bit-for-bit



#test failure check

def failcheck(job_path, path):
    failedt = 0
    input = open(job_path + path, 'r')
    while 1:
        line = input.readline()
        if line == '':
            return
        if line == "-- Final Status Test Results --\n":
            line = input.readline()
            for letter in line:
                if letter == 'F':
                    return 1
                else:
                    return 0



def emptycheck(checkpath):
        noplot = 0
        file = 'temp.txt'
        comline = 'ncdump -c ' + checkpath + '> temp.txt'
        subprocess.call(comline, shell = True)

        input = open('temp.txt', 'r')
        line = ''

        for line in input:
                if line.find('// (0 currently)') != -1:
                        noplot = 1

        return noplot
