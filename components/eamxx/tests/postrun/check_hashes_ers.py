#!/usr/bin/env python3

import os, sys, re, glob, subprocess

def readall(fn):
    with open(fn,'r') as f:
        txt = f.read()
    return txt

def greptxt(pattern, txt):
    return re.findall('(?:' + pattern + ').*', txt, flags=re.MULTILINE)

def grep(pattern, fn):
    txt = readall(fn)
    return greptxt(pattern, txt)

def get_log_glob_from_atm_modelio(case_dir):
    filename = case_dir + os.path.sep + 'CaseDocs' + os.path.sep + 'atm_modelio.nml'
    ln = grep('diro = ', filename)[0]
    run_dir = ln.split()[2].split('"')[1]
    ln = grep('logfile = ', filename)[0]
    atm_log_fn = ln.split()[2].split('"')[1]
    id = atm_log_fn.split('.')[2]
    return run_dir + os.path.sep + 'e3sm.log.' + id + '*'

def get_hash_lines(fn):
    rlns = subprocess.run(['zgrep', 'exxhash', fn], capture_output=True)
    rlns = rlns.stdout.decode().split('\n')
    lns = []
    if len(rlns) == 0: return lns
    for rln in rlns:
        pos = rln.find('exxhash')
        lns.append(rln[pos:])
    return lns

def parse_time(hash_ln):
    return hash_ln.split()[1:3]

def all_equal(t1, t2):
    if len(t1) != len(t2): return False
    for i in range(len(t1)):
        if t1[i] != t2[i]: return False
    return True

def find_first_index_at_time(lns, time):
    for i, ln in enumerate(lns):
        t = parse_time(ln)
        if all_equal(time, t): return i
    return None

def diff(l1, l2):
    diffs = []
    for i in range(len(l1)):
        if l1[i] != l2[i]:
            diffs.append((l1[i], l2[i]))
    return diffs
    
def main(case_dir):
    # Look for the two e3sm.log files.
    glob_pat = get_log_glob_from_atm_modelio(case_dir)
    e3sm_fns = glob.glob(glob_pat)
    if len(e3sm_fns) == 0:
        print('Could not find e3sm.log files with glob string {}'.format(glob_pat))
        return False
    e3sm_fns.sort()
    if len(e3sm_fns) == 1:
        # This is the first run. Exit and wait for the second
        # run. (POSTRUN_SCRIPT is called after each of the two runs.)
        print('Exiting on first run.')
        return True
    print('Diffing base {} and restart {}'.format(e3sm_fns[0], e3sm_fns[1]))

    # Because of the prefixed 1: and 2: on some systems, we can't just use
    # zdiff.
    lns = []
    for f in e3sm_fns:
        lns.append(get_hash_lines(f))
    time = parse_time(lns[1][0])
    time_idx = find_first_index_at_time(lns[0], time)
    if time_idx is None:
        print('Could not find a start time.')
        return False
    lns[0] = lns[0][time_idx:]
    if len(lns[0]) != len(lns[1]):
        print('Number of hash lines starting at restart time do not agree.')
        return False
    diffs = diff(lns[0], lns[1])

    # Flushed prints to e3sm.log can sometimes conflict with other
    # output. Permit up to 'thr' diffs so we don't fail due to badly printed
    # lines. This isn't a big loss in checking because an ERS_Ln22 second run
    # writes > 1000 hash lines, and a true loss of BFBness is nearly certain to
    # propagate to a large number of subsequent hashes.
    thr = 5
    if len(lns[0]) < 100: thr = 0

    ok = True
    if len(diffs) > thr:
        print('DIFF')
        print(diffs[-10:])
        ok = False
    else:
        print('OK')
        
    return ok

case_dir = sys.argv[1]
sys.exit(0 if main(case_dir) else 1)
