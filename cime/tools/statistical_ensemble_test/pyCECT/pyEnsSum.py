#!/usr/bin/env python
import ConfigParser
import sys, getopt, os 
import numpy as np 
import Nio 
import time
import re
from asaptools.partition import EqualStride, Duplicate,EqualLength
import asaptools.simplecomm as simplecomm 
import pyEnsLib

#This routine creates a summary file from an ensemble of CAM
#output files

def main(argv):

    print 'Running pyEnsSum!'

    # Get command line stuff and store in a dictionary
    s = 'tag= compset= esize= tslice= res= sumfile= indir= sumfiledir= mach= verbose jsonfile= mpi_enable maxnorm gmonly popens cumul regx= startMon= endMon= fIndex='
    optkeys = s.split()
    try: 
        opts, args = getopt.getopt(argv, "h", optkeys)
    except getopt.GetoptError:
        pyEnsLib.EnsSum_usage()
        sys.exit(2)

    # Put command line options in a dictionary - also set defaults
    opts_dict={}
    
    # Defaults
    opts_dict['tag'] = ''
    opts_dict['compset'] = ''
    opts_dict['mach'] = ''
    opts_dict['esize'] = 151
    opts_dict['tslice'] = 0
    opts_dict['res'] = ''
    opts_dict['sumfile'] = 'ens.summary.nc'
    opts_dict['indir'] = './'
    opts_dict['sumfiledir'] = './'
    opts_dict['jsonfile'] = ''
    opts_dict['verbose'] = True
    opts_dict['mpi_enable'] = False
    opts_dict['maxnorm'] = False
    opts_dict['gmonly'] = False
    opts_dict['popens'] = False
    opts_dict['cumul'] = False
    opts_dict['regx'] = 'test'
    opts_dict['startMon'] = 1
    opts_dict['endMon'] = 1
    opts_dict['fIndex'] = 151

    # This creates the dictionary of input arguments 
    opts_dict = pyEnsLib.getopt_parseconfig(opts,optkeys,'ES',opts_dict)

    verbose = opts_dict['verbose']

    st = opts_dict['esize']
    esize = int(st)

    if (verbose == True):
        print opts_dict
        print 'Ensemble size for summary = ', esize

    if not (opts_dict['tag'] and opts_dict['compset'] and opts_dict['mach'] or opts_dict['res']):
       print 'Please specify --tag, --compset, --mach and --res options'
       sys.exit()
       
    # Now find file names in indir
    input_dir = opts_dict['indir']
    # The var list that will be excluded
    ex_varlist=[]

    # Create a mpi simplecomm object
    if opts_dict['mpi_enable']:
        me=simplecomm.create_comm()
    else:
        me=simplecomm.create_comm(not opts_dict['mpi_enable'])


    if me.get_rank() == 0:
	if opts_dict['jsonfile']:
	    # Read in the excluded var list
	    ex_varlist=pyEnsLib.read_jsonlist(opts_dict['jsonfile'],'ES')

    # Broadcast the excluded var list to each processor
    if opts_dict['mpi_enable']:
	ex_varlist=me.partition(ex_varlist,func=Duplicate(),involved=True)
        
    in_files=[]
    if(os.path.exists(input_dir)):
        # Get the list of files
        in_files_temp = os.listdir(input_dir)
        in_files=sorted(in_files_temp)
        #print in_files
        # Make sure we have enough
        num_files = len(in_files)
        if (verbose == True):
            print 'Number of files in input directory = ', num_files
        if (num_files < esize):
            print 'Number of files in input directory (',num_files,\
                ') is less than specified ensemble size of ', esize
            sys.exit(2)
        if (num_files > esize):
            print 'NOTE: Number of files in ', input_dir, \
                'is greater than specified ensemble size of ', esize ,\
                '\nwill just use the first ',  esize, 'files'
    else:
        print 'Input directory: ',input_dir,' not found'
        sys.exit(2)

    if opts_dict['cumul']:
        if opts_dict['regx']:
           in_files_list=get_cumul_filelist(opts_dict,opts_dict['indir'],opts_dict['regx'])
        in_files=me.partition(in_files_list,func=EqualLength(),involved=True)
        if me.get_rank()==0:
           print 'in_files=',in_files

    # Open the files in the input directory
    o_files=[]
    for onefile in in_files[0:esize]:
        if (os.path.isfile(input_dir+'/' + onefile)):
            o_files.append(Nio.open_file(input_dir+'/' + onefile,"r"))
        else:
            print "COULD NOT LOCATE FILE "+ input_dir + onefile + "! EXITING...."
            sys.exit() 

    # Store dimensions of the input fields
    if (verbose == True):
        print "Getting spatial dimensions"
    nlev = -1
    ncol = -1
    nlat = -1
    nlon = -1
    lonkey=''
    latkey=''
    # Look at first file and get dims
    input_dims = o_files[0].dimensions
    ndims = len(input_dims)

    for key in input_dims:
        if key == "lev":
            nlev = input_dims["lev"]
        elif key == "ncol":
            ncol = input_dims["ncol"]
        elif (key == "nlon") or (key =="lon"):
            nlon = input_dims[key]
            lonkey=key
        elif (key == "nlat") or (key == "lat"):
            nlat = input_dims[key]
            latkey=key
        
    if (nlev == -1) : 
        print "COULD NOT LOCATE valid dimension lev => EXITING...."
        sys.exit() 

    if (( ncol == -1) and ((nlat == -1) or (nlon == -1))):
        print "Need either lat/lon or ncol  => EXITING...."
        sys.exit()            

    # Check if this is SE or FV data
    if (ncol != -1):
        is_SE = True 
    else:
        is_SE = False    

    # Make sure all files have the same dimensions
    if (verbose == True):
        print "Checking dimensions across files...."
        print 'lev = ', nlev
        if (is_SE == True):
            print 'ncol = ', ncol
        else:
            print 'nlat = ', nlat
            print 'nlon = ', nlon

    for count, this_file in enumerate(o_files):
        input_dims = this_file.dimensions     
        if (is_SE == True):
            if ( nlev != int(input_dims["lev"]) or ( ncol != int(input_dims["ncol"]))):
                print "Dimension mismatch between ", in_files[0], 'and', in_files[0], '!!!'
                sys.exit() 
        else:
            if ( nlev != int(input_dims["lev"]) or ( nlat != int(input_dims[latkey]))\
                  or ( nlon != int(input_dims[lonkey]))): 
                print "Dimension mismatch between ", in_files[0], 'and', in_files[0], '!!!'
                sys.exit() 

    # Get 2d vars, 3d vars and all vars (For now include all variables) 
    vars_dict = o_files[0].variables
    # Remove the excluded variables (specified in json file) from variable dictionary
    if ex_varlist:
	for i in ex_varlist:
            if i in vars_dict:
	       del vars_dict[i]
    num_vars = len(vars_dict)
    if (verbose == True):
        print 'Number of variables (including metadata) found =  ', num_vars
    str_size = 0

    d2_var_names = []
    d3_var_names = []
    num_2d = 0
    num_3d = 0

    # Which are 2d, which are 3d and max str_size 
    for k,v in vars_dict.iteritems():  
        var = k
        vd = v.dimensions # all the variable's dimensions (names)
        vr = v.rank # num dimension
        vs = v.shape # dim values
        is_2d = False
        is_3d = False
        if (is_SE == True): # (time, lev, ncol) or (time, ncol)
	    if ((vr == 2) and (vs[1] == ncol)):  
		is_2d = True 
		num_2d += 1
	    elif ((vr == 3) and (vs[2] == ncol and vs[1] == nlev )):  
		is_3d = True 
		num_3d += 1
        else: # (time, lev, nlon, nlon) or (time, nlat, nlon)
            if ((vr == 3) and (vs[1] == nlat and vs[2] == nlon)):  
                is_2d = True 
                num_2d += 1
            elif ((vr == 4) and (vs[2] == nlat and vs[3] == nlon and vs[1] == nlev )):  
                is_3d = True 
                num_3d += 1
        if (is_3d == True) :
            str_size = max(str_size, len(k))
            d3_var_names.append(k)
        elif  (is_2d == True):    
            str_size = max(str_size, len(k))
            d2_var_names.append(k)


    # Now sort these and combine (this sorts caps first, then lower case - 
    # which is what we want)
    d2_var_names.sort()       
    d3_var_names.sort()


    # All vars is 3d vars first (sorted), the 2d vars
    all_var_names = list(d3_var_names)
    all_var_names += d2_var_names
    n_all_var_names = len(all_var_names)

    if (verbose == True):
        print 'num vars = ', n_all_var_names, '(3d = ', num_3d, ' and 2d = ', num_2d, ")"

    # Create new summary ensemble file
    this_sumfile = opts_dict["sumfile"]

    if (verbose == True):
        print "Creating ", this_sumfile, "  ..."
    if(me.get_rank() ==0 | opts_dict["popens"]):
	if os.path.exists(this_sumfile):
	    os.unlink(this_sumfile)

	opt = Nio.options()
	opt.PreFill = False
	opt.Format = 'NetCDF4Classic'
	nc_sumfile = Nio.open_file(this_sumfile, 'w', options=opt)

	# Set dimensions
	if (verbose == True):
	    print "Setting dimensions ....."
	if (is_SE == True):
	    nc_sumfile.create_dimension('ncol', ncol)
	else:
	    nc_sumfile.create_dimension('nlat', nlat)
	    nc_sumfile.create_dimension('nlon', nlon)
	nc_sumfile.create_dimension('nlev', nlev)
	nc_sumfile.create_dimension('ens_size', esize)
	nc_sumfile.create_dimension('nvars', num_3d + num_2d)
	nc_sumfile.create_dimension('nvars3d', num_3d)
	nc_sumfile.create_dimension('nvars2d', num_2d)
	nc_sumfile.create_dimension('str_size', str_size)

	# Set global attributes
	now = time.strftime("%c")
	if (verbose == True):
	    print "Setting global attributes ....."
	setattr(nc_sumfile, 'creation_date',now)
	setattr(nc_sumfile, 'title', 'CAM verification ensemble summary file')
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
        if not opts_dict['gmonly']:
	    if (is_SE == True):
		v_ens_avg3d = nc_sumfile.create_variable("ens_avg3d", 'f', ('nvars3d', 'nlev', 'ncol'))
		v_ens_stddev3d = nc_sumfile.create_variable("ens_stddev3d", 'f', ('nvars3d', 'nlev', 'ncol'))
		v_ens_avg2d = nc_sumfile.create_variable("ens_avg2d", 'f', ('nvars2d', 'ncol'))
		v_ens_stddev2d = nc_sumfile.create_variable("ens_stddev2d", 'f', ('nvars2d', 'ncol'))
	    else:
		v_ens_avg3d = nc_sumfile.create_variable("ens_avg3d", 'f', ('nvars3d', 'nlev', 'nlat', 'nlon'))
		v_ens_stddev3d = nc_sumfile.create_variable("ens_stddev3d", 'f', ('nvars3d', 'nlev', 'nlat', 'nlon'))
		v_ens_avg2d = nc_sumfile.create_variable("ens_avg2d", 'f', ('nvars2d', 'nlat', 'nlon'))
		v_ens_stddev2d = nc_sumfile.create_variable("ens_stddev2d", 'f', ('nvars2d', 'nlat', 'nlon'))

	    v_RMSZ = nc_sumfile.create_variable("RMSZ", 'f', ('nvars', 'ens_size'))
	v_gm = nc_sumfile.create_variable("global_mean", 'f', ('nvars', 'ens_size'))
	v_loadings_gm = nc_sumfile.create_variable('loadings_gm','f',('nvars','nvars'))
	v_mu_gm = nc_sumfile.create_variable('mu_gm','f',('nvars',))
	v_sigma_gm = nc_sumfile.create_variable('sigma_gm','f',('nvars',))
	v_sigma_scores_gm = nc_sumfile.create_variable('sigma_scores_gm','f',('nvars',))


	# Assign vars, var3d and var2d
	if (verbose == True):
	    print "Assigning vars, var3d, and var2d ....."

	eq_all_var_names =[]
	eq_d3_var_names = []
	eq_d2_var_names = []

	l_eq = len(all_var_names)
	for i in range(l_eq):
	    tt = list(all_var_names[i])
	    l_tt = len(tt)
	    if (l_tt < str_size):
		extra = list(' ')*(str_size - l_tt)
		tt.extend(extra)
	    eq_all_var_names.append(tt)

	l_eq = len(d3_var_names)
	for i in range(l_eq):
	    tt = list(d3_var_names[i])
	    l_tt = len(tt)
	    if (l_tt < str_size):
		extra = list(' ')*(str_size - l_tt)
		tt.extend(extra)
	    eq_d3_var_names.append(tt)

	l_eq = len(d2_var_names)
	for i in range(l_eq):
	    tt = list(d2_var_names[i])
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
	lev_data = vars_dict["lev"]
	v_lev = lev_data

    # Form ensembles, each missing one member; compute RMSZs and global means
    #for each variable, we also do max norm also (currently done in pyStats)
    tslice = opts_dict['tslice']

    if not opts_dict['cumul']:
        # Partition the var list
        var3_list_loc=me.partition(d3_var_names,func=EqualStride(),involved=True)
        var2_list_loc=me.partition(d2_var_names,func=EqualStride(),involved=True)
    else:
        var3_list_loc=d3_var_names
        var2_list_loc=d2_var_names

    # Calculate global means #
    if (verbose == True):
        print "Calculating global means ....."
    if not opts_dict['cumul']:
        gm3d,gm2d = pyEnsLib.generate_global_mean_for_summary(o_files,var3_list_loc,var2_list_loc , is_SE, False,opts_dict)      
    if (verbose == True):
        print "Finish calculating global means ....."

    # Calculate RMSZ scores  
    if (verbose == True):
        print "Calculating RMSZ scores ....."
    if (not opts_dict['gmonly']) | (opts_dict['cumul']):
        zscore3d,zscore2d,ens_avg3d,ens_stddev3d,ens_avg2d,ens_stddev2d,temp1,temp2=pyEnsLib.calc_rmsz(o_files,var3_list_loc,var2_list_loc,is_SE,opts_dict)    

    # Calculate max norm ensemble
    if opts_dict['maxnorm']:
	if (verbose == True):
	    print "Calculating max norm of ensembles ....."
	pyEnsLib.calculate_maxnormens(opts_dict,var3_list_loc)
	pyEnsLib.calculate_maxnormens(opts_dict,var2_list_loc)

    if opts_dict['mpi_enable'] & ( not opts_dict['popens']):

        if not opts_dict['cumul']:
	    # Gather the 3d variable results from all processors to the master processor
	    slice_index=get_stride_list(len(d3_var_names),me)
	 
	    # Gather global means 3d results
	    gm3d=gather_npArray(gm3d,me,slice_index,(len(d3_var_names),len(o_files)))

	    if not opts_dict['gmonly']:
		# Gather zscore3d results
		zscore3d=gather_npArray(zscore3d,me,slice_index,(len(d3_var_names),len(o_files)))

		# Gather ens_avg3d and ens_stddev3d results
		shape_tuple3d=get_shape(ens_avg3d.shape,len(d3_var_names),me.get_rank())
		ens_avg3d=gather_npArray(ens_avg3d,me,slice_index,shape_tuple3d) 
		ens_stddev3d=gather_npArray(ens_stddev3d,me,slice_index,shape_tuple3d) 

	    # Gather 2d variable results from all processors to the master processor
	    slice_index=get_stride_list(len(d2_var_names),me)

	    # Gather global means 2d results
	    gm2d=gather_npArray(gm2d,me,slice_index,(len(d2_var_names),len(o_files)))

	    if not opts_dict['gmonly']:
		# Gather zscore2d results
		zscore2d=gather_npArray(zscore2d,me,slice_index,(len(d2_var_names),len(o_files)))

		# Gather ens_avg3d and ens_stddev2d results
		shape_tuple2d=get_shape(ens_avg2d.shape,len(d2_var_names),me.get_rank())
		ens_avg2d=gather_npArray(ens_avg2d,me,slice_index,shape_tuple2d) 
		ens_stddev2d=gather_npArray(ens_stddev2d,me,slice_index,shape_tuple2d) 

        else:
	    gmall=np.concatenate((temp1,temp2),axis=0)
            gmall=pyEnsLib.gather_npArray_pop(gmall,me,(me.get_size(),len(d3_var_names)+len(d2_var_names)))
    # Assign to file:
    if me.get_rank() == 0 | opts_dict['popens'] :
        if not opts_dict['cumul']:
	    gmall=np.concatenate((gm3d,gm2d),axis=0)
	    if not opts_dict['gmonly']:
		Zscoreall=np.concatenate((zscore3d,zscore2d),axis=0)
		v_RMSZ[:,:]=Zscoreall[:,:]
	    if not opts_dict['gmonly']:
		if (is_SE == True):
		    v_ens_avg3d[:,:,:]=ens_avg3d[:,:,:]
		    v_ens_stddev3d[:,:,:]=ens_stddev3d[:,:,:]
		    v_ens_avg2d[:,:]=ens_avg2d[:,:]
		    v_ens_stddev2d[:,:]=ens_stddev2d[:,:]
		else:
		    v_ens_avg3d[:,:,:,:]=ens_avg3d[:,:,:,:]
		    v_ens_stddev3d[:,:,:,:]=ens_stddev3d[:,:,:,:]
		    v_ens_avg2d[:,:,:]=ens_avg2d[:,:,:]
		    v_ens_stddev2d[:,:,:]=ens_stddev2d[:,:,:]
        else:
            gmall_temp=np.transpose(gmall[:,:])
            gmall=gmall_temp
	mu_gm,sigma_gm,standardized_global_mean,loadings_gm,scores_gm=pyEnsLib.pre_PCA(gmall)
	v_gm[:,:]=gmall[:,:]
	v_mu_gm[:]=mu_gm[:]
	v_sigma_gm[:]=sigma_gm[:].astype(np.float32)
	v_loadings_gm[:,:]=loadings_gm[:,:]
	v_sigma_scores_gm[:]=scores_gm[:]
                
	print "All Done"

def get_cumul_filelist(opts_dict,indir,regx):
   if not opts_dict['indir']:
      print 'input dir is not specified'
      sys.exit(2)
   #regx='(pgi(.)*-(01|02))'
   regx_list=["mon","gnu","pgi"]
   all_files=[]
   for prefix in regx_list: 
       for i in range(opts_dict['fIndex'],opts_dict['fIndex']+opts_dict['esize']/3):
	   for j in range(opts_dict['startMon'],opts_dict['endMon']+1):
	       mon_str=str(j).zfill(2)
	       regx='(^'+prefix+'(.)*'+str(i)+'(.)*-('+mon_str+'))'
	       print 'regx=',regx
	       res=[f for f in os.listdir(indir) if re.search(regx,f)]
	       in_files=sorted(res)
	       all_files.extend(in_files)
   print "all_files=",all_files
   #in_files=res
   return all_files
   
   
      
   

#
# Get the shape of all variable list in tuple for all processor
# 
def get_shape(shape_tuple,shape1,rank):
    lst=list(shape_tuple)
    lst[0]=shape1
    shape_tuple=tuple(lst)
    return shape_tuple
 
#
# Get the mpi partition list for each processor
#
def get_stride_list(len_of_list,me):
    slice_index=[]
    for i in range(me.get_size()):
	index_arr=np.arange(len_of_list)
	slice_index.append(index_arr[i::me.get_size()])
    return slice_index

# 
# Gather arrays from each processor by the var_list to the master processor and make it an array
#
def gather_npArray(npArray,me,slice_index,array_shape):
    the_array=np.zeros(array_shape,dtype=np.float32)
    if me.get_rank()==0:
	k=0
	for j in slice_index[me.get_rank()]:
	     the_array[j,:]=npArray[k,:]
	     k=k+1
    for i in range(1,me.get_size()):
	if me.get_rank() == 0:
	    rank,npArray=me.collect()
	    k=0
	    for j in slice_index[rank]:
		the_array[j,:]=npArray[k,:]
		k=k+1
    if me.get_rank() != 0: 
	message={"from_rank":me.get_rank(),"shape":npArray.shape}
	me.collect(npArray)
    me.sync()
    return the_array
        
if __name__ == "__main__":
    main(sys.argv[1:])
