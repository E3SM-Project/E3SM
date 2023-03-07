#!/usr/bin/env python3

import os, sys, re, glob

def readall(fn):
    with open(fn,'r') as f:
        txt = f.read()
    return txt

def greptxt(pattern, txt):
    return re.findall('(?:' + pattern + ').*', txt, flags=re.MULTILINE)

def grep(pattern, fn):
    txt = readall(fn)
    return greptxt(pattern, txt)

def get_run_dir(case_dir):
    filename = case_dir + os.path.sep + 'CaseDocs' + os.path.sep + 'atm_modelio.nml'
    ln = grep('diro = ', filename)[0]
    return ln.split()[2].split('"')[1]

def get_cpl_log(run_dir):
    filenames = glob.glob(run_dir + os.path.sep + 'cpl.log.*')
    if len(filenames) == 0: return None
    filenames.sort()
    return filenames[-1]

def uncompress(filename):
    if '.gz' in filename:
        os.system('gunzip {}'.format(filename))
        return filename[:-3]
    return filename

def main(case_dir):
    run_dir = get_run_dir(case_dir)
    cpl_fn = get_cpl_log(run_dir)
    if cpl_fn is None: return False
    cpl_fn = uncompress(cpl_fn)
    print('Using log file {}'.format(cpl_fn))
    alarms = grep('nlmap>.*ALARM', cpl_fn)
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
