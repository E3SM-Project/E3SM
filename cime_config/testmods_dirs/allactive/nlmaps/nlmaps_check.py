#!/usr/bin/env python3

"""
POSTRUN_SCRIPT that searches the cpl.log files for the string 'ALARM' that
is output if nlmaps_verbosity > 0 and there is a failure to preserve properties
sufficiently accurately.
"""

import os, sys, re, glob

def findall_multiline(pattern, fn):
    """
    Apply re.findall to the text in the file FN with PATTERN as the group.
    """
    with open(fn,'r') as f:
        txt = f.read()
    return re.findall('(?:' + pattern + ').*', txt, flags=re.MULTILINE)

def get_run_dir(case_dir):
    """
    Given the case directory, look in atm_modelio.nml for the run directory.
    """
    filename = case_dir + os.path.sep + 'CaseDocs' + os.path.sep + 'atm_modelio.nml'
    ln = findall_multiline('diro = ', filename)[0]
    return ln.split()[2].split('"')[1]

def get_cpl_log(run_dir):
    """
    Given the run directory, return the most recent cpl.log* file.
    """
    filenames = glob.glob(run_dir + os.path.sep + 'cpl.log.*')
    if len(filenames) == 0: return None
    filenames.sort()
    return filenames[-1]

def uncompress(filename):
    """
    Run gunzip on FILENAME if FILENAME ends in '.gz'; otherwise, return without
    doing anything
    """
    if '.gz' in filename:
        os.system('gunzip {}'.format(filename))
        return filename[:-3]
    return filename

def main(case_dir):
    """
    Search the cpl.log files for the string 'ALARM' output if nlmaps_verbosity >
    0 and there is a failure to preserve properties sufficiently accurately.
    """
    run_dir = get_run_dir(case_dir)
    cpl_fn = get_cpl_log(run_dir)
    if cpl_fn is None:
        print('Did not find any cpl.log files.')
        return False
    cpl_fn = uncompress(cpl_fn)
    print(f'Using log file {cpl_fn}')
    alarms = findall_multiline('nlmap>.*ALARM', cpl_fn)
    ok = len(alarms) == 0
    if not ok:
        for ln in alarms:
            print(ln)
    return ok

good = main(sys.argv[1])

if good:
    print('PASS')
    sys.exit(0)
else:
    print('FAIL')
    sys.exit(1) # non-0 exit will make test RUN phase fail, as desired
