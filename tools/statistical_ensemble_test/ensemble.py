#!/usr/bin/python
import os, sys, getopt
import random
#import utility
import single_run

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


def main(argv):

    caller = 'ensemble.py'

    opts_dict, case_flags  = process_args_dict(caller, argv)

    #default is verification mode (3 runs)
    run_type = 'verify'
    clone_count = 2

    uf = opts_dict['uf']

    #check for run_type change (i.e., if doing ensemble instead of verify)
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
        
    #now create cases
    thisdir = os.getcwd()
   
    #create first case - then clone
    if runtype == 'verify':
        opts_dict['pertlim'] = get_pertlim_uf(rand_ints[0])
    else: #full ensemble
        opts_dict['pertlim'] = "0"
   
    #first case
    single_case(opts_dict, case_flags)

    #now clone FIX THIS
#    case_root=os.path.dirname(case_name)
#    case_pfx=os.path.basename(case_name)
#    for i in range(1,clonecount):
#        iens='{0:03d}'.format(i)
#        if runtype == 'validation':
#           pertlim=get_pertlim_uf(rand_ints(i))
#        else:
#           pertlim=get_pertlim_uf(i)
#        case1_name=case_pfx+"."+iens
#        case1=case_root+"/"+case1_name

#        os.chdir(scripts_root)
#        print "=== SCRIPTS_ROOT ==="
#        print scripts_root
#        command='scripts_root/create_clone --keepexe --case case1 --clone case'
#        ret=os.system(command)


        # Get value for EXEROOT from $CASE
        # Note return string is "EXEROOT = $EXEROOT"
        if test_suite == TRUE:
           os.chdir(scripts_root+"/"+case)
        else:
           os.chdir(case)


if __name__ == "__main__":
    main(sys.argv[1:])
