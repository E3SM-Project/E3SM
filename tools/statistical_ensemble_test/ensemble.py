#!/usr/bin/python
import os, sys, getopt
import random
import utility

#==============================================================================
# set up and submit 12-month (original) or 9-time step (uf) run.  then create 
# clones for a complete ensemble or a set of (3) test cases
#==============================================================================

#generate <num_pick> positive random integers in [0, end-1]
#can't have any duplicates
def random_pick(num_pick, end):
    ar = range(0, end)
    rand_list = random.sample(ar, num_pick)
    #for i in rand_list:
    #    print i
    return rand_list

#get the pertlim corressponding to the random int
def get_pertlim_uf(rand_num):
    i = rand_num
    if i == 0:
       ptlim = 0
    else:
       j = 2*((i - 1)/100) + 101
       k = (i - 1)%100
       if i%2 != 0:
          ll = j + (k/2)*18
          ippt = '{0:03d}'.format(ll)
          ptlim = "0."+ippt+"d-13"
       else:
          ll = j + ((k-1)/2)*18
          ippt = '{0:03d}'.format(ll)
          ptlim = "-0."+ippt+"d-13"
    return ptlim 

#I'm here
def create_cases(case_path, run_type, case_name):
    os.chdir(os.path.dirname(case_path))
    thisdir=os.getcwd()
    singlerun=thisdir+"/single_run.py"
    if not os.path.exists(singlerun):
       print "ERROR: cannot find script to produce single run in "+thisdir
    
    ret=1
    if runtype == 'verify':
       firstpertlim = get_pertlim_uf(30)
       command = 'python '+singlerun+' --pertlim '+firstpertlim
       print command
       #ret=os.system(command)
    else:
       print command
       command = 'python '+singlerun
       #ret=os.system(command)

    if ret != 0:
       print "Exit: "+ str(ret)
       sys.exit(0)

    ######### END OF BUILDING ROOT CASE, NOW CLONING

    case_root=os.path.dirname(case_name)
    case_pfx=os.path.basename(case_name)
    for i in range(1,clonecount):
        iens='{0:03d}'.format(i)
        if runtype == 'validation':
           pertlim=get_pertlim_uf(rand_ints(i))
        else:
           pertlim=get_pertlim_uf(i)
        case1_name=case_pfx+"."+iens
        case1=case_root+"/"+case1_name

        os.chdir(scripts_root)
        print "=== SCRIPTS_ROOT ==="
        print scripts_root
        command='scripts_root/create_clone --keepexe --case case1 --clone case'
        ret=os.system(command)


        # Get value for EXEROOT from $CASE
        # Note return string is "EXEROOT = $EXEROOT"
        if test_suite == TRUE:
           os.chdir(scripts_root+"/"+case)
        else:
           os.chdir(case)


def main(argv):

    caller = 'ensemble.py'

    #default is verification mode (3 runs)
    run_type = 'verify'
    clone_count = 2

    # Pull in and analyze the command line arguements
    s='case= mach= project= compiler= compset= res= uf nb ns ensemble= verbose silent test multi-driver pecount= nist= mpilib= pesfile= gridfile= srcroot= output-root= script-root= queue= user-modes-dir= input-dir= pertlim= h'
    optkeys=s.split()
    try: 
        opts, args = getopt.getopt(argv,"hf:",optkeys)
    except getopt.GetoptError:
        disp_usage(caller)
        sys.exit(2)

    for opt,arg in opts:
        if opt == '-h':
            disp_usage(caller)
            sys.exit()

    #opts_dict and defaults    
    opts_dict={}
    opts_dict['res']='f19_f19'
    opts_dict['compset']='F2000'
    opts_dict['walltime']='04:30'
    opts_dict['pertlim']= 0
    opts_dict['nb'] = False
    opts_dict['ns'] = False
    opts_dict['uf'] = False
    opts_dict['ensemble'] = 0

    opts_dict=utility.getopt_parseconfig(opts, optkeys, caller, opts_dict)

    uf = opts_dict['uf']

    #check for run_type change
    ens_size = opts_dict['ensemble'] 
    if  ens_size > 0:
        run_type = 'ensemble'
        clone_count = ens_size - 1
        if ens_size > 999:
            print 'Error: cannot have an ensemble size greater than 999.'
            sys.exit()

    #generate 3 random pertlims for verify
    if run_type == 'verify':
        if uf:
            end_range = 350
        else:
            end_range = 150
        rand_ints = random_pick(3, end_range)
       
        
    #create cases

    
    
