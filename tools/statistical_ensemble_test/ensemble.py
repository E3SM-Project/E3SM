#!/usr/bin/python
import os, sys, getopt
import random
from single_run import process_args_dict, single_case

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
    print('STATUS: creating first case ...')

    #create first case - then clone
    if run_type == 'verify':
        opts_dict['pertlim'] = get_pertlim_uf(rand_ints[0])
    else: #full ensemble
        opts_dict['pertlim'] = "0"
   
    #first case
    single_case(opts_dict, case_flags)

    #now clone
    print('STATUS: cloning additional cases ...')

    #scripts dir
    root_dir = os.path.basename(__file__)
    ret = os.chdir(root_dir)
    ret = os.chdir('../../scripts')
    scripts_dir = os.getcwd()

    #we know case name ends in '.000' (already checked)
    clone_case = opts_dict['case']
    case_pfx = clone_case[:-4]

    for i in range(1, clone_count + 1): #1: clone_count
        if runtype == 'verify':
           pertlim = get_pertlim_uf(rand_ints[i])
        else: #full ensemble
           pertlim = get_pertlim_uf(i)

        iens = '{0:03d}'.format(i)
        new_case = case_pfx + "." + iens
        full_new_case = scripts_dir +"/" + clone_name

        os.chdir(scripts_dir)
        print ("STATUS: creating cloned case: " + clone_name)

        command = scripts_dir + "/create_clone --keepexe --case " + full_new_case + " --clone" + clone_case
        ret = os.system(command)

        print ("STATUS: running setup for cloned case: " + clone_name)
        os.chdir(full_new_case)
        command = './case_setup'
        ret = os.system(command)

        #modify the perlim in the file
        if runtype == 'verify': #remove old pertlim first
            f = open("user_nl_cam.nl","r+")
            all_lines = f.readlines()
            f.seek(0)
            for line in all_lines:
                if line.find("pertlim") == -1:
                    f.write(line)
            f.truncate()
            f.close()

        #now append new pertlim
        with open("user_nl_cam", "a") as f:
            text = "pertlim = " + opts_dict['pertlim']
            f.write(text)

        #preview namelists
        command = './preview_namelists'
        ret = os.system(command)
    
        #submit?
        if opts_dict["ns"] == False:
            command = './case.submit'
            ret = os.system(command)

    #Final output
    if runtype == "verify":
        print ("STATUS: ---VERIFICATION CASES COMPLETE---")
        print ("Set up three cases using the following pertlim values:")
        print get_pertlim_uf(rand_ints[0]) + '   ' + get_pertlim_uf(rand_ints[1]) + "   " + get_pertlim_uf(rand_ints[2])
    else:
       print ("STATUS: --ENSEMBLE CASES COMPLETE---")




if __name__ == "__main__":
    main(sys.argv[1:])
