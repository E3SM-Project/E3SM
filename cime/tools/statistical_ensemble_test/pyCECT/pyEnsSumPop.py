#!/usr/bin/env python
import ConfigParser
import sys, getopt, os 
import numpy as np 
import Nio 
import time
import re
from asaptools.partition import EqualStride, Duplicate
import asaptools.simplecomm as simplecomm 
import pyEnsLib
#import pdb

def main(argv):
    print 'Running pyEnsSumPop!'

    # Get command line stuff and store in a dictionary
    s = 'nyear= nmonth= npert= tag= res= mach= compset= sumfile= indir= tslice= verbose jsonfile= mpi_enable nrand= rand seq= jsondir='
    optkeys = s.split()
    try: 
        opts, args = getopt.getopt(argv, "h", optkeys)
    except getopt.GetoptError:
        pyEnsLib.EnsSumPop_usage()
        sys.exit(2)

    # Put command line options in a dictionary - also set defaults
    opts_dict={}

    # Defaults
    opts_dict['tag'] = 'cesm2_0_0'
    opts_dict['compset'] = 'G'
    opts_dict['mach'] = 'cheyenne'
    opts_dict['tslice'] = 0 
    opts_dict['nyear'] = 1
    opts_dict['nmonth'] = 12
    opts_dict['npert'] = 40
    opts_dict['nbin'] = 40
    opts_dict['minrange'] = 0.0
    opts_dict['maxrange'] = 4.0
    opts_dict['res'] = 'T62_g17'
    opts_dict['sumfile'] = 'ens.pop.summary.nc'
    opts_dict['indir'] = './'
    opts_dict['jsonfile'] = ''
    opts_dict['verbose'] = True
    opts_dict['mpi_enable'] = False
    opts_dict['zscoreonly'] = True
    opts_dict['popens'] = True
    opts_dict['nrand'] = 40 
    opts_dict['rand'] = False
    opts_dict['seq'] = 0 
    opts_dict['jsondir'] = './' 

    # This creates the dictionary of input arguments 
    #print "before parseconfig"
    opts_dict = pyEnsLib.getopt_parseconfig(opts,optkeys,'ESP',opts_dict)

    verbose = opts_dict['verbose']
    nbin = opts_dict['nbin']

    if verbose:
        print "opts_dict = "
        print opts_dict
       
    # Now find file names in indir
    input_dir = opts_dict['indir']

    # Create a mpi simplecomm object
    if opts_dict['mpi_enable']:
        me=simplecomm.create_comm()
    else:
        me=simplecomm.create_comm(not opts_dict['mpi_enable'])
    if opts_dict['jsonfile']:
        # Read in the included var list
        Var2d,Var3d=pyEnsLib.read_jsonlist(opts_dict['jsonfile'],'ESP')
        str_size=0
        for str in Var3d:
            if str_size < len(str):
               str_size=len(str)
        for str in Var2d:
            if str_size < len(str):
               str_size=len(str)


    in_files=[]
    if(os.path.exists(input_dir)):
        # Pick up the 'nrand' random number of input files to generate summary files
        if opts_dict['rand']:
           in_files=pyEnsLib.Random_pickup_pop(input_dir,opts_dict,opts_dict['nrand'])
        else:    
           # Get the list of files
           in_files_temp = os.listdir(input_dir)
           in_files=sorted(in_files_temp)
        num_files = len(in_files)
#        if (verbose == True):
#            print in_files

    else:
        print 'ERROR: Input directory: ',input_dir,' not found'
        sys.exit(2)

    # Create a mpi simplecomm object
    if opts_dict['mpi_enable']:
        me=simplecomm.create_comm()
    else:
        me=simplecomm.create_comm(not opts_dict['mpi_enable'])
    #Partition the input file list 
    in_file_list=me.partition(in_files,func=EqualStride(),involved=True)

    
    # Open the files in the input directory
    o_files=[]
    for onefile in in_file_list:
        if (os.path.isfile(input_dir+'/' + onefile)):
            o_files.append(Nio.open_file(input_dir+'/' + onefile,"r"))
        else:
            print "ERROR: Could not locate file: "+ input_dir + onefile 
            sys.exit() 


    #print in_file_list

    # Store dimensions of the input fields
    if (verbose == True):
        print "Getting spatial dimensions"
    nlev = -1
    nlat = -1
    nlon = -1

    # Look at first file and get dims
    input_dims = o_files[0].dimensions
    ndims = len(input_dims)

    # Make sure all files have the same dimensions
    if (verbose == True):
        print "Checking dimensions ..."
    for key in input_dims:
        if key == "z_t":
            nlev = input_dims["z_t"]
        elif key == "nlon":
            nlon = input_dims["nlon"]
        elif key == "nlat":
            nlat = input_dims["nlat"]

    for count, this_file in enumerate(o_files):
        input_dims = this_file.dimensions     
        if ( nlev != int(input_dims["z_t"]) or ( nlat != int(input_dims["nlat"]))\
              or ( nlon != int(input_dims["nlon"]))):
            print "Dimension mismatch between ", in_file_list[0], 'and', in_file_list[count], '!!!'
            sys.exit() 


    # Create new summary ensemble file
    this_sumfile = opts_dict["sumfile"]

    if verbose:
       print "Creating ", this_sumfile, "  ..."
    if (me.get_rank() == 0 ):
       if os.path.exists(this_sumfile):
           os.unlink(this_sumfile)
       opt =Nio.options()
       opt.PreFill = False
       opt.Format = 'NetCDF4Classic'

       nc_sumfile = Nio.open_file(this_sumfile, 'w', options=opt)

       # Set dimensions
       if (verbose == True):
           print "Setting dimensions ....."
       nc_sumfile.create_dimension('nlat', nlat)
       nc_sumfile.create_dimension('nlon', nlon)
       nc_sumfile.create_dimension('nlev', nlev)
       nc_sumfile.create_dimension('time',None)
       nc_sumfile.create_dimension('ens_size', opts_dict['npert'])
       nc_sumfile.create_dimension('nbin', opts_dict['nbin'])
       nc_sumfile.create_dimension('nvars', len(Var3d) + len(Var2d))
       nc_sumfile.create_dimension('nvars3d', len(Var3d))
       nc_sumfile.create_dimension('nvars2d', len(Var2d))
       nc_sumfile.create_dimension('str_size', str_size)

       # Set global attributes
       now = time.strftime("%c")
       if (verbose == True):
           print "Setting global attributes ....."
       setattr(nc_sumfile, 'creation_date',now)
       setattr(nc_sumfile, 'title', 'POP verification ensemble summary file')
       setattr(nc_sumfile, 'tag', opts_dict["tag"]) 
       setattr(nc_sumfile, 'compset', opts_dict["compset"]) 
       setattr(nc_sumfile, 'resolution', opts_dict["res"]) 
       setattr(nc_sumfile, 'machine', opts_dict["mach"]) 

       # Create variables
       if (verbose == True):
           print "Creating variables ....."
       v_lev = nc_sumfile.create_variable("lev", 'f', ('nlev',))
       v_vars = nc_sumfile.create_variable("vars", 'S1', ('nvars', 'str_size'))
       v_var3d = nc_sumfile.create_variable("var3d", 'S1', ('nvars3d', 'str_size'))
       v_var2d = nc_sumfile.create_variable("var2d", 'S1', ('nvars2d', 'str_size'))
       v_time = nc_sumfile.create_variable("time",'d',('time',))
       v_ens_avg3d = nc_sumfile.create_variable("ens_avg3d", 'f', ('time','nvars3d', 'nlev', 'nlat', 'nlon'))
       v_ens_stddev3d = nc_sumfile.create_variable("ens_stddev3d", 'f', ('time','nvars3d', 'nlev', 'nlat', 'nlon'))
       v_ens_avg2d = nc_sumfile.create_variable("ens_avg2d", 'f', ('time','nvars2d', 'nlat', 'nlon'))
       v_ens_stddev2d = nc_sumfile.create_variable("ens_stddev2d", 'f', ('time','nvars2d', 'nlat', 'nlon'))

       v_RMSZ = nc_sumfile.create_variable("RMSZ", 'f', ('time','nvars', 'ens_size','nbin'))
       if not opts_dict['zscoreonly']:
          v_gm = nc_sumfile.create_variable("global_mean", 'f', ('time','nvars', 'ens_size'))



       # Assign vars, var3d and var2d
       if (verbose == True):
           print "Assigning vars, var3d, and var2d ....."

       eq_all_var_names =[]
       eq_d3_var_names = []
       eq_d2_var_names = []
       all_var_names = list(Var3d)
       all_var_names += Var2d
       l_eq = len(all_var_names)
       for i in range(l_eq):
           tt = list(all_var_names[i])
           l_tt = len(tt)
           if (l_tt < str_size):
               extra = list(' ')*(str_size - l_tt)
               tt.extend(extra)
           eq_all_var_names.append(tt)

       l_eq = len(Var3d)
       for i in range(l_eq):
           tt = list(Var3d[i])
           l_tt = len(tt)
           if (l_tt < str_size):
               extra = list(' ')*(str_size - l_tt)
               tt.extend(extra)
           eq_d3_var_names.append(tt)

       l_eq = len(Var2d)
       for i in range(l_eq):
           tt = list(Var2d[i])
           l_tt = len(tt)
           if (l_tt < str_size):
               extra = list(' ')*(str_size - l_tt)
               tt.extend(extra)
           eq_d2_var_names.append(tt)

       v_vars[:] = eq_all_var_names[:]
       v_var3d[:] = eq_d3_var_names[:]
       v_var2d[:] = eq_d2_var_names[:]

       # Time-invarient metadata
       if (verbose == True):
           print "Assigning time invariant metadata ....."
       vars_dict = o_files[0].variables
       lev_data = vars_dict["z_t"]
       v_lev = lev_data
       
    # Time-varient metadata
    if verbose:
       print "Assigning time variant metadata ....."
    vars_dict = o_files[0].variables
    time_value = vars_dict['time']
    time_array = np.array([time_value])
    time_array = pyEnsLib.gather_npArray_pop(time_array,me,(me.get_size(),))
    if me.get_rank() == 0:
       v_time[:]=time_array[:]

    #Assign zero values to first time slice of RMSZ and avg and stddev for 2d & 3d 
    #in case of a calculation problem before finishing
    e_size = opts_dict['npert']
    b_size =  opts_dict['nbin']
    z_ens_avg3d=np.zeros((len(Var3d),nlev,nlat,nlon),dtype=np.float32)
    z_ens_stddev3d=np.zeros((len(Var3d),nlev,nlat,nlon),dtype=np.float32)
    z_ens_avg2d=np.zeros((len(Var2d),nlat,nlon),dtype=np.float32)
    z_ens_stddev2d=np.zeros((len(Var2d),nlat,nlon),dtype=np.float32)
    z_RMSZ = np.zeros(((len(Var3d)+len(Var2d)),e_size,b_size), dtype=np.float32)
    if me.get_rank() == 0 :
        v_RMSZ[0,:,:,:]=z_RMSZ[:,:,:]
        v_ens_avg3d[0,:,:,:,:]=z_ens_avg3d[:,:,:,:]
        v_ens_stddev3d[0,:,:,:,:]=z_ens_stddev3d[:,:,:,:]
        v_ens_avg2d[0,:,:,:]=z_ens_avg2d[:,:,:]
        v_ens_stddev2d[0,:,:,:]=z_ens_stddev2d[:,:,:]

    # Calculate global mean, average, standard deviation and rmse 
    if verbose:
       if not opts_dict['zscoreonly']: 
           print "Calculating global means ....."
    is_SE = False
    tslice=0
    if not opts_dict['zscoreonly']:
       gm3d,gm2d = pyEnsLib.generate_global_mean_for_summary(o_files,Var3d,Var2d, is_SE,False,opts_dict)
    if verbose:
        if not opts_dict['zscoreonly']:
            print "Finish calculating global means ....."

    # Calculate RMSZ scores  
    if (verbose == True):
       print "Calculating RMSZ scores ....."
    zscore3d,zscore2d,ens_avg3d,ens_stddev3d,ens_avg2d,ens_stddev2d,temp1,temp2=pyEnsLib.calc_rmsz(o_files,Var3d,Var2d,is_SE,opts_dict)    

    if (verbose == True):
        print "Finished with RMSZ scores ....."

    # Collect from all processors
    if opts_dict['mpi_enable'] :
        # Gather the 3d variable results from all processors to the master processor
        # Gather global means 3d results
        if not opts_dict['zscoreonly']:
           gmall=np.concatenate((gm3d,gm2d),axis=0)
           #print "before gather, gmall.shape=",gmall.shape
           gmall=pyEnsLib.gather_npArray_pop(gmall,me,(me.get_size(),len(Var3d)+len(Var2d),len(o_files)))
        zmall=np.concatenate((zscore3d,zscore2d),axis=0)
        zmall=pyEnsLib.gather_npArray_pop(zmall,me,(me.get_size(),len(Var3d)+len(Var2d),len(o_files),nbin))
        #print 'zmall=',zmall
        
        #print "after gather, gmall.shape=",gmall.shape
        ens_avg3d=pyEnsLib.gather_npArray_pop(ens_avg3d,me,(me.get_size(),len(Var3d),nlev,(nlat),nlon))
        ens_avg2d=pyEnsLib.gather_npArray_pop(ens_avg2d,me,(me.get_size(),len(Var2d),(nlat),nlon))
        ens_stddev3d=pyEnsLib.gather_npArray_pop(ens_stddev3d,me,(me.get_size(),len(Var3d),nlev,(nlat),nlon))
        ens_stddev2d=pyEnsLib.gather_npArray_pop(ens_stddev2d,me,(me.get_size(),len(Var2d),(nlat),nlon))

    # Assign to file:
    if me.get_rank() == 0 :
        #Zscoreall=np.concatenate((zscore3d,zscore2d),axis=0)
        v_RMSZ[:,:,:,:]=zmall[:,:,:,:]
        v_ens_avg3d[:,:,:,:,:]=ens_avg3d[:,:,:,:,:]
        v_ens_stddev3d[:,:,:,:,:]=ens_stddev3d[:,:,:,:,:]
        v_ens_avg2d[:,:,:,:]=ens_avg2d[:,:,:,:]
        v_ens_stddev2d[:,:,:,:]=ens_stddev2d[:,:,:,:]
        if not opts_dict['zscoreonly']:
           v_gm[:,:,:]=gmall[:,:,:]
        print "All done"

if __name__ == "__main__":
    main(sys.argv[1:])
