#!/usr/bin/env python3

import os, sys, re
#import numpy as np

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

def gather_energy_data(atm_log_fn):
    #process only water, energy leaks, ignore statistics on cflx, shflx
    dgmean1 = {'nstep': [], 'te': []}
    dgmean2 = {'nstep': [], 'dt': [], 'tw': [], 'wdiff': [], 'ediff': []}

    #file lines
    # nstep, te        2   0.25846082320679531E+10   0.25845989127171683E+10  -0.12885108527156013E-03   0.98516311352621502E+05
    # n, dt, W tot mass [kg/m2]        1   0.72000000000000000E+04   0.23608972524228417E+02
    # n, W flux, dWater [kg/m2]        1   0.69364453203771062E-01   0.69458133498910668E-01
    # n, W flux-dWater [kg/m2]        1  -0.93680295139605962E-04
    # n, W cflx*dt loss [kg/m2]        1  -0.83649589965878964E-04
    # n, E d(TE)/dt, RR [W/m2]        1   0.14448525118182053E+02   0.14477581732178772E+02
    # n, E difference [W/m2]        1  -0.29056613996718994E-01
    # n, E shf loss [W/m2]        1   0.29056614514536561E-01

    l1 = "nstep, te "
    l2 = "n, dt, W tot mass [kg/m2]    "
    l3 = "n, W flux-dWater [kg/m2]     "
    l4 = "n, E diff"
    with open(atm_log_fn, 'r') as f:
        for ln in f:
            if l1 in ln: 
                nstep = float(ln.split()[2])
                te = float(ln.split()[3])
                dgmean1['nstep'].append(nstep)
                dgmean1['te'].append(te)
                #print (nstep)
                #print (rr)
        #reset
        f.seek(0)

        for ln in f:
            if l2 in ln: 
                nstep = float(ln.split()[6])
                dt = float(ln.split()[7])
                tw = float(ln.split()[8])
                dgmean2['nstep'].append(nstep)
                dgmean2['dt'].append(dt)
                dgmean2['tw'].append(tw)
        #reset
        f.seek(0)

        for ln in f:
            if l3 in ln:
                diff = float(ln.split()[5])
                dgmean2['wdiff'].append(diff)
        #reset
        f.seek(0)

        for ln in f:
            if l4 in ln:
                diff = float(ln.split()[5])
                dgmean2['ediff'].append(diff)
        #reset
        f.seek(0)

    return dgmean1, dgmean2

def conservativeW(start_time, nstep, tw, wdiff, tol, verbose):
   
    short_tw=tw[start_time:]
    short_wdiff=wdiff[start_time:]
    res = []
    for i in range(len(short_wdiff)):
       res.append(short_wdiff[i] / short_tw[i])
    short_wdiff = res

    aver = sum(short_wdiff)/len(short_wdiff)
    short_wdiff_abs = [abs(ele) for ele in short_wdiff]

    #stdd = np.std(np.array(short_wdiff))
    if (verbose):
        print('Water leak rel:', short_wdiff)
        print('starting index is {}'.format(start_time))
        print('number of samples: {}'.format(len(short_wdiff)))
        print('Water leak average: {:1.3e}'.format(aver))
        print('abs(Water leak): min {:1.3e}, max {:1.3e}'.format(min(short_wdiff_abs),max(short_wdiff_abs)))
        #maxabs = max(abs(rr))
        maxx = abs(aver)
    return maxx <= tol


def conservativeE(s1, s2, te, nstep, ediff, tol, verbose):

    #TE will be used for rel computations
    #it is not essential that TE and ediff are outputed for lagged time steps, TE is practically the same
    #for all time steps

    short_te=te[s1:]
    short_ediff=ediff[s2:]
    res = []
    for i in range(len(short_ediff)):
       res.append(short_ediff[i] / short_te[i])
    short_ediff = res 

    aver = sum(short_ediff)/len(short_ediff)
    short_ediff_abs = [abs(ele) for ele in short_ediff]

    aver = sum(short_ediff)/len(short_ediff)
    #stdd = np.std(np.array(short_ediff))
    if (verbose):
        print('Energy leak rel:', short_ediff)
        print('starting indices are {},{}'.format(s1,s2))
        print('number of samples: {}'.format(len(short_ediff)))
        print('Ediff average: {:1.3e}'.format(aver))
        print('abs(Ediff): min {:1.3e}, max {:1.3e}'.format(min(short_ediff_abs),max(short_ediff_abs)))
        #maxabs = max(abs(rr))
        maxx = abs(aver)
    return maxx <= tol



#############################################################3


print('1st sys arg {}'.format(sys.argv[1]))

case_dir = sys.argv[1]

#uncomment in the repo!!!!!!
run_dir = read_atm_modelio(case_dir)

print('run dir is {}'.format(run_dir))

atm_fn = get_atm_log(run_dir)
#debug in run folder
#atm_fn = get_atm_log(".")

atm_fn = uncompress(atm_fn)
print('Using log file {}'.format(atm_fn))
[d1,d2] = gather_energy_data(atm_fn)

goodw = (conservativeW(0,d2['nstep'],d2['tw'],d2['wdiff'],2e-7,True))
#d1 and d2 arrays are of different size
#d1 array has values for nstep 0, 1, 2, 3
#d2 array has values for nstep 2, 3 (but they are named as 1 and 2, they are shifted by 1 wrt te arrays)
goode = (conservativeE(2,0,d1['te'],d2['nstep'],d2['ediff'],2e-12,True))

if (goodw and goode):
    print('PASS')
    sys.exit(0)
else:
    print('FAIL')
    sys.exit(1) # non-0 exit will make test RUN phase fail, as desired
