#!/usr/bin/env python

import os, sys, re

def readall(fn):
    with open(fn,'r') as f:
        txt = f.read()
    return txt

def greptxt(pattern, txt):
    return re.findall('(?:' + pattern + ').*', txt, flags=re.MULTILINE)

def grep(pattern, fn):
    txt = readall(fn)
    return greptxt(pattern, txt)

def read_atm_modelio(case_dir):
    filename = case_dir + os.path.sep + 'CaseDocs' + os.path.sep + 'atm_modelio.nml'
    ln = grep('diro = ', filename)[0]
    return ln.split()[2].split('"')[1]

def get_atm_log(run_dir):
    filenames = os.listdir(run_dir)
    atm_fn = None
    for f in filenames:
        if 'atm.log' in f:
            atm_fn = f
            break
    return run_dir + os.path.sep + atm_fn

def uncompress(filename):
    if '.gz' in filename:
        os.system('gunzip {}'.format(filename))
        return filename[:-3]
    return filename

def parse_tracer_index(atm_log_fn, tracer_name):
    state = 0
    with open(atm_log_fn, 'r') as f:
        for ln in f:
            if state == 0:
                if 'Advected constituent list:' in ln: state = 1
            elif state == 1:
                if tracer_name in ln:
                    toks = ln.split()
                    tracer_idx = int(toks[0])
                    break
    return tracer_idx

def parsetime(ln):
    toks = ln.split()
    val = float(toks[3])
    unit = toks[4]
    if unit == '[s]': val = val/86400.0
    return val

def parseqmass(ln):
    toks = ln.split()
    return float(toks[-1])
    
def gather_mass_data(atm_log_fn, tracer_idx):
    d = {'day': [], 'dryM': [], 'qmass': []}
    qline = 'qv({:3d})='.format(tracer_idx)
    with open(atm_log_fn, 'r') as f:
        for ln in f:
            if 'nstep=' in ln: t = parsetime(ln)
            elif qline in ln: qmass = parseqmass(ln)
            elif 'dry M' in ln:
                dryM = float(ln.split()[3])
                d['day'].append(t)
                d['qmass'].append(qmass)
                d['dryM'].append(dryM)
    return d

def conservative(label, time, mass, tol_per_year, verbose):
    t_delta = (time[-1] - time[0])/365.0
    mass_rel = abs(mass[-1] - mass[0])/abs(mass[0])
    thr = tol_per_year*t_delta
    if (verbose):
        print('{:20s}: mass rel err {:1.3e} tol: {:1.3e}'.format(label, mass_rel, thr))
    return mass_rel <= thr

tracer_name = 'CO2_FFF'

case_dir = sys.argv[1]
run_dir = read_atm_modelio(case_dir)
atm_fn = get_atm_log(run_dir)
atm_fn = uncompress(atm_fn)
print('Using log file {}'.format(atm_fn))
tracer_idx = parse_tracer_index(atm_fn, tracer_name)
print('Tracer {} has index {}'.format(tracer_name, tracer_idx))
d = gather_mass_data(atm_fn, tracer_idx)
good = (conservative('dry M', d['day'], d['dryM'], 1e-11, True) and
        conservative('tracer {}'.format(tracer_idx), d['day'], d['qmass'], 2e-13, True))

if good:
    print('PASS')
    sys.exit(0)
else:
    print('FAIL')
    sys.exit(1) # non-0 exit will make test RUN phase fail, as desired
