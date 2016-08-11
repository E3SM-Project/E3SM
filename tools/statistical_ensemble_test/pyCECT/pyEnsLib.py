#!/usr/bin/env python
import ConfigParser
import sys, getopt, os 
import numpy as np 
import Nio 
import time
import re
import json
import random

#
# Create RMSZ zscores for 
#
def calc_rmsz(o_files,openfile,var_name3d,var_name2d,tslice,is_SE,verbose):

    input_dims = openfile.dimensions
    nlev=input_dims["lev"]
    # Create array variables based on is_SE
    if (is_SE == True):
      ncol=input_dims["ncol"]
      npts2d=ncol
      npts3d=nlev*ncol
      output3d = np.zeros((len(o_files),nlev,ncol),dtype=np.float32)
      output2d = np.zeros((len(o_files),ncol),dtype=np.float32)
      ens_avg3d=np.zeros((len(var_name3d),nlev,ncol),dtype=np.float32)
      ens_stddev3d=np.zeros((len(var_name3d),nlev,ncol),dtype=np.float32)
      ens_avg2d=np.zeros((len(var_name2d),ncol),dtype=np.float32)
      ens_stddev2d=np.zeros((len(var_name2d),ncol),dtype=np.float32)
    else:
      nlon = input_dims["nlon"]
      nlat = input_dims["nlat"]
      npts2d=nlat*nlon
      npts3d=nlev*nlat*nlon
      output3d = np.zeros((len(o_files),nlev,nlat,nlon),dtype=np.float32)
      output2d = np.zeros((len(o_files),nlat,nlon),dtype=np.float32)
      ens_avg3d=np.zeros((len(var_name3d),nlev,nlat,nlon),dtype=np.float32)
      ens_stddev3d=np.zeros((len(var_name3d),nlev,nlat,nlon),dtype=np.float32)
      ens_avg2d=np.zeros((len(var_name2d),nlat,nlon),dtype=np.float32)
      ens_stddev2d=np.zeros((len(var_name2d),nlat,nlon),dtype=np.float32)
    Zscore3d = np.zeros((len(var_name3d),len(o_files)),dtype=np.float32) 
    Zscore2d = np.zeros((len(var_name2d),len(o_files)),dtype=np.float32) 
    avg3d={}
    stddev3d={}
    avg2d={}
    stddev2d={}
    indices = np.arange(0,len(o_files),1)
   
    
    for vcount,vname in enumerate(var_name3d):
      #Read in vname's data of all files
      for fcount, this_file in enumerate(o_files):
        data=this_file.variables[vname]
        if (is_SE == True):
          output3d[fcount,:,:]=data[tslice,:,:]
      
        else:
          output3d[fcount,:,:,:]=data[tslice,:,:,:]
      #Generate avg, stddev and zscore for 3d variable
      for fcount,this_file in enumerate(o_files):

        new_index=np.where(indices!=fcount)
        ensemble3d = output3d[new_index] 
        avg3d=np.average(ensemble3d,axis=0)
        stddev3d=np.std(ensemble3d,axis=0,dtype=np.float64) 

        flag3d = False
        count3d = 0
        count3d,ret_val=calc_Z(output3d[fcount].astype(np.float64),avg3d.astype(np.float64),stddev3d.astype(np.float64),count3d,flag3d)
        Zscore=np.sum(np.square(ret_val))

        if (count3d < npts3d):
          Zscore3d[vcount,fcount]=np.sqrt(Zscore/(npts3d-count3d))
        else:
          print "WARNING: no variance in "+vname
      #Generate ens_avg and esn_stddev to store in the ensemble summary file
      ens_avg3d[vcount]=np.average(output3d,axis=0)
      ens_stddev3d[vcount]=np.std(output3d.astype(np.float64),axis=0,dtype=np.float64)

    for vcount,vname in enumerate(var_name2d):
      #Read in vname's data of all files
      for fcount, this_file in enumerate(o_files):
        data=this_file.variables[vname]
        if (is_SE == True):
          output2d[fcount,:]=data[tslice,:]
      
        else:
          output2d[fcount,:,:]=data[tslice,:,:]
      #Generate avg, stddev and zscore for 3d variable
      for fcount,this_file in enumerate(o_files):

        new_index=np.where(indices!=fcount)
        ensemble2d = output2d[new_index] 
        avg2d=np.average(ensemble2d,axis=0)
        stddev2d=np.std(ensemble2d,axis=0,dtype=np.float64) 

        flag2d = False
        count2d = 0
        count2d,ret_val=calc_Z(output2d[fcount].astype(np.float64),avg2d.astype(np.float64),stddev2d.astype(np.float64),count2d,flag2d)
        Zscore=np.sum(np.square(ret_val))

        if (count2d < npts2d):
          Zscore2d[vcount,fcount]=np.sqrt(Zscore/(npts2d-count2d))
        else:
          print "WARNING: no variance in "+vname
      #Generate ens_avg and esn_stddev to store in the ensemble summary file
      ens_avg2d[vcount]=np.average(output2d,axis=0)
      ens_stddev2d[vcount]=np.std(output2d.astype(np.float64),axis=0,dtype=np.float64)
     
    return Zscore3d,Zscore2d,ens_avg3d,ens_stddev3d,ens_avg2d,ens_stddev2d

#
# Calculate rmsz score by compare the run file with the ensemble summary file 
#
def calculate_raw_score(k,v,npts3d,npts2d,ens_avg,ens_stddev,is_SE):
  count=0
  Zscore=0
  has_zscore=True
  if k in ens_avg:
    if is_SE:
	if ens_avg[k].ndim == 1:
	  npts=npts2d
	else:
	  npts=npts3d
    else:
	if ens_avg[k].ndim == 2:
	  npts=npts2d
	else:
	  npts=npts3d
        
    count,return_val=calc_Z(v,ens_avg[k].astype(np.float64),ens_stddev[k].astype(np.float64),count,False) 
    Zscore=np.sum(np.square(return_val))
    if npts == count: 
      Zscore=0
    else:
      Zscore=np.sqrt(Zscore/(npts-count))
  else:
    has_zscore=False
   
  return Zscore,has_zscore

#
# Create some variables and call a function to calculate PCA
#
def pre_PCA(gm):
    gm_len=gm.shape
    nvar=gm_len[0]
    nfile=gm_len[1]
    mu_gm=np.average(gm,axis=1)
    sigma_gm=np.std(gm,axis=1,dtype=np.float32)
    standardized_global_mean=np.zeros(gm.shape,dtype=np.float32)
    scores_gm=np.zeros(gm.shape,dtype=np.float32)

    for var in range(nvar):
      for file in range(nfile):
        standardized_global_mean[var,file]=(gm[var,file]-mu_gm[var])/sigma_gm[var]
    loadings_gm=princomp(standardized_global_mean)

    #now do coord transformation on the standardized meana to get the scores
    scores_gm=np.dot(loadings_gm.T,standardized_global_mean)
    sigma_scores_gm =np.std(scores_gm,axis=1,dtype=np.float64)
    return mu_gm,sigma_gm,standardized_global_mean,loadings_gm.astype(np.float32),sigma_scores_gm.astype(np.float32)

#
# Performs principal components analysis  (PCA) on the p-by-n data matrix A
# rows of A correspond to (p) variables AND cols of A correspond to the (n) tests 
# assume already standardized
#
# Returns the loadings: p-by-p matrix, each column containing coefficients 
# for one principal component. 
#
def princomp(standardized_global_mean):
  # find covariance matrix (will be pxp)
  co_mat= np.cov(standardized_global_mean)    
  # Calculate evals and evecs of covariance matrix (evecs are also pxp)
  [evals, evecs] = np.linalg.eig(co_mat)
  # Above may not be sorted - sort largest first
  new_index = np.argsort(evals)[::-1]
  evecs = evecs[:,new_index]
  evals = evals[new_index]

  return evecs

#          
# Calculate (val-avg)/stddev and exclude zero value
#          
def calc_Z(val,avg,stddev,count,flag):
  return_val=np.empty(val.shape,dtype=np.float32,order='C')
  tol =1e-12   
  if stddev[(stddev > tol)].size ==0:
    if flag: 
      print "WARNING: ALL standard dev = 0"
      flag = False
    count =count + stddev[(stddev <= tol)].size
    return_val = 0.
  else:        
    if stddev[(stddev <= tol)].size > 0:
      if flag:
        print "WARNING: some standard dev = 0"
        flag =False
      count =count + stddev[(stddev <= tol)].size
      return_val[np.where(stddev <= tol)]=0.
      return_val[np.where(stddev > tol)]= (val[np.where(stddev> tol)]-avg[np.where(stddev> tol)])/stddev[np.where(stddev>tol)]
    else:
      return_val=(val-avg)/stddev
  return count,return_val

#
# Read a json file for the excluded list of variables
#
def read_jsonlist(metajson):
  fd=open(metajson)
  metainfo = json.load(fd)
  varList = metainfo['ExcludedVar']
  return varList


# 
# Calculate Normalized RMSE metric
#
def calc_nrmse(orig_array,comp_array):
  
  orig_size=orig_array.size
  sumsqr=np.sum(np.square(orig_array.astype(np.float64)-comp_array.astype(np.float64)))
  rng=np.max(orig_array)-np.min(orig_array)
  if abs(rng) < 1e-18:
    rmse=0.0
  else:
    rmse=np.sqrt(sumsqr/orig_size)/rng

  return rmse

#
# Calculate weighted global mean for one level of CAM output
#
def area_avg(data, weight, is_SE):

    #TO DO: tke into account missing values
    if (is_SE == True):
        a = np.average(data, weights=weight)
    else: #FV
        #a = wgt_areaave(data, weight, 1.0, 1)
        #weights are for lat 
        a_lat = np.average(data,axis=0, weights=weight)
        a = np.average(a_lat)
    return a

#
# Open input files,compute area_wgts, and then loop through all files to call calc_global_means_for_onefile
#
def generate_global_mean_for_summary(o_files,var_name3d,var_name2d,tslice,is_SE,verbose):

    #openfile - should have already been opened by Nio.open_file()
    n3d = len(var_name3d)
    n2d = len(var_name2d)
    tot = n3d + n2d

    gm3d = np.zeros((n3d,len(o_files)),dtype=np.float32)
    gm2d = np.zeros((n2d,len(o_files)),dtype=np.float32)

    #get dimensions and compute area_wgts
    input_dims = o_files[0].dimensions
    nlev = input_dims["lev"]

   
    if (is_SE == True):
        ncol = input_dims["ncol"]
        output3d = np.zeros((nlev, ncol))
        output2d = np.zeros(ncol)
        area_wgt = np.zeros(ncol)
        area = o_files[0].variables["area"]
        area_wgt[:] = area[:]
        total = np.sum(area_wgt)
        area_wgt[:] /= total
    else:
        nlon = input_dims["lon"]
        nlat = input_dims["lat"]
        output3d = np.zeros((nlev, nlat, nlon))
        output2d = np.zeros((nlat, nlon))
        area_wgt = np.zeros(nlat) #note gaues weights are length nlat
        gw = o_files[0].variables["gw"]
        area_wgt[:] = gw[:]
    
    #loop through the input file list to calculate global means
    #var_name3d=[]
    for fcount,fname in enumerate(o_files):
        gm3d[:,fcount],gm2d[:,fcount]=calc_global_mean_for_onefile(fname,area_wgt,var_name3d,var_name2d,output3d,output2d,tslice,is_SE,nlev,verbose)
        # Generate global mean for pepsi challenge data timeseries daily files, they all are 2d variables
        #var_name2d=[]
        #for k,v in fname.variables.iteritems():
        #  if v.typecode() == 'f': 
        #    var_name2d.append(k)
        #    fout = open(k+"_32.txt","w")
        #  if k == 'time':
        #    ntslice=v[:]
        #for i in np.nditer(ntslice):
        #  print i
        #  temp1,temp2=calc_global_mean_for_onefile(fname,area_wgt,var_name3d,var_name2d,output3d,output2d,int(i),is_SE,nlev,verbose)
        #  fout.write(str(temp2[0])+'\n')

    return gm3d,gm2d

#
# Calculate global means for one input file
#
def calc_global_mean_for_onefile(fname, area_wgt,var_name3d, var_name2d,output3d,output2d, tslice, is_SE, nlev,verbose):
    
    ##openfile - should have already been opened by Nio.open_file()
    n3d = len(var_name3d)
    n2d = len(var_name2d)
    #tot = n3d + n2d

    gm3d = np.zeros((n3d),dtype=np.float32)
    gm2d = np.zeros((n2d),dtype=np.float32)

    #calculate global mean for each 3D variable 
    for count, vname in enumerate(var_name3d):
	#if (verbose == True):
	#    print "calculating GM for variable ", vname
	gm_lev = np.zeros(nlev)
	data = fname.variables[vname]
	if (is_SE == True):
	    output3d[:,:] = data[tslice,:,:] 
	    for k in range(nlev):
		gm_lev[k] = area_avg(output3d[k,:], area_wgt, is_SE)
	else:
	    output3d[:,:,:] = data[tslice,:,:,:] 
	    for k in range(nlev):
		gm_lev[k] = area_avg(output3d[k,:,:], area_wgt, is_SE)
	#note: averaging over levels should probably be pressure-weighted(TO DO)        
	gm3d[count] = np.mean(gm_lev)         

	
    #calculate global mean for each 2D variable 
    for count, vname in enumerate(var_name2d):
	#if (verbose == True):
	#    print "calculating GM for variable ", vname
	data = fname.variables[vname]
	if (is_SE == True):
	    output2d[:] = data[tslice,:] 
	    gm2d_mean = area_avg(output2d[:], area_wgt, is_SE)
	else:
	    output2d[:,:] = data[tslice,:,:] 
	    gm2d_mean = area_avg(output2d[:,:], area_wgt, is_SE)
	gm2d[count]=gm2d_mean

    return gm3d,gm2d        

#
# Read variable values from ensemble summary file
#
def read_ensemble_summary(ens_file):
  if(os.path.isfile(ens_file)):
     fens = Nio.open_file(ens_file,"r")
  else:
     print 'file ens summary: ',ens_file,' Not found'
     sys.exit(2)

  ens_avg={}
  ens_stddev={}
  ens_var_name=[]
  ens_rmsz={}
  ens_gm={}
  #mu_gm={}
  #sigma_gm={}
  #loadings_gm={}
  #sigma_scores_gm={}

  # Retrieve the variable list from ensemble file
  for k,v in fens.variables.iteritems():
    if k== 'vars':
      for i in v[0:len(v)]:
        l=0
        for j in i:
          if j: 
            l=l+1
        ens_var_name.append(i[0:l].tostring().strip())
    elif k== 'var3d':
      num_var3d=len(v)
    elif k== 'var2d':
      num_var2d=len(v)
  for k,v in fens.variables.iteritems():
    # Retrieve the ens_avg3d or ens_avg2d array
    if k == 'ens_avg3d' or k=='ens_avg2d':
      if k== 'ens_avg2d':
        m=num_var3d
      else:
        m=0
      if v:
        for i in v[0:len(v)]:
          temp_name=ens_var_name[m]
          ens_avg[temp_name] = i
          m=m+1
    
    # Retrieve the ens_stddev3d or ens_stddev2d  array
    elif k == 'ens_stddev3d' or k == 'ens_stddev2d':
      if k== 'ens_stddev2d':
        m=num_var3d
      else:
        m=0
      if v:
        for i in v[0:len(v)]:
          temp_name=ens_var_name[m]
          ens_stddev[temp_name] = i
          m=m+1
    # Retrieve the RMSZ score array
    elif k == 'RMSZ':
      m=0
      for i in v[0:len(v)]:
        temp_name=ens_var_name[m]
        ens_rmsz[temp_name]=i
        m=m+1
    elif k == 'global_mean':
      m=0
      for i in v[0:len(v)]:
        temp_name=ens_var_name[m]
        ens_gm[temp_name]=i
        m=m+1
     
    elif k == 'mu_gm':
      mu_gm=np.zeros((num_var3d+num_var2d),dtype=np.float32)
      mu_gm[:]=v[:]
    elif k == 'sigma_gm':
      sigma_gm=np.zeros((num_var3d+num_var2d),dtype=np.float32)
      sigma_gm[:]=v[:]
    elif k == 'loadings_gm':
      loadings_gm=np.zeros((num_var3d+num_var2d,num_var3d+num_var2d),dtype=np.float32)
      loadings_gm[:,:]=v[:,:]
    elif k == 'sigma_scores_gm':
      sigma_scores_gm=np.zeros((num_var3d+num_var2d),dtype=np.float32)
      sigma_scores_gm[:]=v[:]
   
  return ens_var_name,ens_avg,ens_stddev,ens_rmsz,ens_gm,num_var3d,mu_gm,sigma_gm,loadings_gm,sigma_scores_gm


#
# Get the ncol and nlev value from cam run file
#
def get_ncol_nlev(frun):
  input_dims=frun.dimensions
  ncol = -1
  nlev = -1
  ilev = -1
  nlat = -1
  nlon = -1
  icol = -1
  ilat = -1
  ilon = -1
  for k,v in input_dims.iteritems():
    if k == 'lev':
      nlev = v
    if k == 'ncol':
      ncol = v
    if k == 'lat':
      nlat = v
      
    if k == 'lon':
      nlon = v
    
  if ncol == -1 :
    one_spatial_dim = False
  else:
    one_spatial_dim = True
  #if nlev == -1 or ncol == -1:
  #  print "Error: cannot find ncol or nlev dimension in "+run_file
  #  sys.exit(2) 

  if one_spatial_dim:
    npts3d=float(nlev*ncol)
    npts2d=float(ncol)
  else:
    npts3d=float(nlev*nlat*nlon)
    npts2d=float(nlat*nlon)
    
  return npts3d,npts2d,one_spatial_dim


#
# Calculate max norm ensemble value for each variable base on all ensemble run files
# the inputdir should only have all ensemble run files
#
def calculate_maxnormens(opts_dict,var_list):
  ifiles=[]
  Maxnormens={}
  threshold=1e-12
  # input file directory
  inputdir=opts_dict['indir']
  
  # the timeslice that we want to process
  tstart=opts_dict['tslice']
  
  # open all files
  for frun_file in os.listdir(inputdir):
    if (os.path.isfile(inputdir+frun_file)):
      ifiles.append(Nio.open_file(inputdir+frun_file,"r"))
    else:
      print "COULD NOT LOCATE FILE "+inputdir+frun_file+" EXISTING"
      sys.exit() 
  comparision={}
  # loop through each variable
  for k in var_list:
    output=[]
    # read all data of variable k from all files
    for f in ifiles:
      v=f.variables
      output.append(v[k][tstart])
    max_val=0
    # open an output file
    outmaxnormens=k+"_ens_maxnorm.txt"
    fout=open(outmaxnormens,"w")
    Maxnormens[k]=[]
   
    # calculate E(i=0:n)(maxnormens[i][x])=max(comparision[i]-E(x=0:n)(output[x])) 
    for n in range(len(ifiles)):
      Maxnormens[k].append(0)
      comparision[k]=ifiles[n].variables[k][tstart]
      for m in range(len(ifiles)):
        max_val=np.max(np.abs(comparision[k]-output[m]))
        if Maxnormens[k][n] < max_val:
          Maxnormens[k][n]=max_val
      range_max=np.max((comparision[k]))
      range_min=np.min((comparision[k]))
      if range_max-range_min < threshold:
	Maxnormens[k][n]=0.
      else:
	Maxnormens[k][n]=Maxnormens[k][n]/(range_max-range_min)
      fout.write(str(Maxnormens[k][n])+'\n')
    strtmp = k + ' : '  + 'ensmax min max' + ' : ' + '{0:9.2e}'.format(min(Maxnormens[k]))+' '+'{0:9.2e}'.format(max(Maxnormens[k]))
    print strtmp
    fout.close()

#
# Parse options from command line or from config file
#
def getopt_parseconfig(opts,optkeys,caller,opts_dict):
  # integer 
  integer   = '-[0-9]+'
  int_p=re.compile(integer)
  # scientific notation
  flt  = '-*[0-9]+\.[0-9]+'
  flt_p=re.compile(flt)

  #opts_dict={}
  #put all optkeys in a dictionary and set default value
  #for key in optkeys:
  #  if key.find('=') != -1:
  #    if key.find('tstart') !=-1 or key.find('tend') != -1:
  #      opts_dict[key[0:key.find('=')]]=0
  #    else:
  #      opts_dict[key[0:key.find('=')]]=''
  #  else:
  #    opts_dict[key]=False

  for opt,arg in opts:
    if opt =='-h' and caller=='CECT':
      CECT_usage()
      sys.exit()
    elif opt == '-h' and caller == 'Ec':
      EnsSum_usage()
      sys.exit()
    elif opt == '-f':
      opts_dict['orig']=arg
    elif opt == '-m':
      opts_dict['reqmeth']=arg
    #parse config file
    elif opt in ("--config"):
      configfile=arg
      config=ConfigParser.ConfigParser()
      config.read(configfile)
      for sec in config.sections():
        for name,value in config.items(sec):
          if sec== 'bool_arg' or sec == 'metrics':
            opts_dict[name]=config.getboolean(sec,name)
          elif sec == 'int_arg':
            opts_dict[name]=config.getint(sec,name)
          elif sec == 'float_arg':
              opts_dict[name]=config.getfloat(sec,name)
          else:
            opts_dict[name]=value
    #parse command line options which might replace the settings in the config file 
    else:
      for k in optkeys:
	if k.find("=") != -1:
	  keyword=k[0:k.find('=')]
	  #if opt.find(keyword)!= -1:
          if opt == '--'+keyword:
            if arg.isdigit():
	      opts_dict[keyword]=int(arg)
            else:
	      if flt_p.match(arg)  :
                opts_dict[keyword]=float(arg)
              elif int_p.match(arg) :
                print arg
                opts_dict[keyword]=int(arg)
              else: 
	        opts_dict[keyword]=arg
	else:
	  #if opt.find(k) != -1:
	  if opt == '--'+k:
	    opts_dict[k]=True
  return opts_dict 
   
#
# Figure out the scores of the 3 new runs, standardized global means, then multiple by the loadings_gm
#
def standardized(gm,mu_gm,sigma_gm,loadings_gm):
    nvar=gm.shape[0]
    nfile=gm.shape[1] 
    standardized_mean=np.zeros(gm.shape,dtype=np.float32)
    for var in range(nvar):
      for file in range(nfile):
        standardized_mean[var,file]=(gm[var,file]-mu_gm[var])/sigma_gm[var]
    new_scores=np.dot(loadings_gm.T,standardized_mean)

    return new_scores

#
# Insert rmsz scores, global mean of new run to the dictionary results
#
def addresults(results,key,value,var,thefile):
    if var in results:
       temp = results[var]
       if key in temp:
          temp2 = temp[key]
          if thefile in temp2:
             temp3 = results[var][key][thefile]
          else:
             temp3={}
       else:
          temp[key]={}
          temp2={}
          temp3={}
       temp3=value
       temp2[thefile]=temp3
       temp[key]=temp2
       results[var]=temp
    else:
       results[var]={}
       results[var][key]={}
       results[var][key][thefile]=value
       
          
    return results

#
# Print out rmsz score failure, global mean failure summary
#
def printsummary(results,key,name,namerange,thefilecount,variables,label):
  thefile='f'+str(thefilecount)
  for k,v in results.iteritems():
    if 'status' in v:
      temp0 = v['status']
      if key in temp0: 
	 if thefile in temp0[key]:
	    temp = temp0[key][thefile]
	    if temp < 1:
               print ' '
	       print k+'       ('+'{0:9.2e}'.format(v[name][thefile])+' outside of ['+'{0:9.2e}'.format(variables[k][namerange][0])+' '+'{0:9.2e}'.format(variables[k][namerange][1])+'])'

#
# Insert the range of rmsz score and global mean of the ensemble summary file to the dictionary variables
#
def addvariables(variables,var,vrange,thearray):
  if var in variables:
     variables[var][vrange]=(np.min(thearray),np.max(thearray))
  else:
     variables[var]={}
     variables[var][vrange]=(np.min(thearray),np.max(thearray))

  return variables

# 
# Evaluate if the new run rmsz score/global mean in the range of rmsz scores/global mean of the ensemble summary
#
def evaluatestatus(name,rangename,variables,key,results,thefile):
  totalcount=0
  for k,v in results.iteritems():
     if name in v and rangename in variables[k]:
       temp0=results[k]
       xrange = variables[k][rangename]
       if v[name][thefile] > xrange[1] or v[name][thefile] < xrange[0]:
         val=0
       else:
         val=1
       if 'status' in temp0:
         temp=temp0['status']
         if key in temp:
            temp2 = temp[key]
         else:
            temp[key] = temp2 = {}
            
         if val == 0:
            totalcount = totalcount+1
         temp2[thefile]=val
         temp[key]=temp2
         results[k]['status']==temp
       else:
         temp0['status']={}
         temp0['status'][key]={}
         temp0['status'][key][thefile]=val
         if val == 0:
            totalcount = totalcount+1

  return totalcount
             
#
# Evaluate if the new run PCA scores pass or fail by comparing with the PCA scores of the ensemble summary
#
def comparePCAscores(ifiles,new_scores,sigma_scores_gm,opts_dict):

   comp_array=np.zeros(new_scores.shape,dtype=np.int32)
   sum=np.zeros(new_scores.shape[0],dtype=np.int32)
   eachruncount=np.zeros(new_scores.shape[1],dtype=np.int32)
   totalcount=0
   sum_index=[]
   print '*********************************************** '
   print 'PCA Test Results'
   print '*********************************************** '
  
   #Test to check if new_scores out of range of sigMul*sigma_scores_gm
   #for i in range(new_scores.shape[0]):
   for i in range(opts_dict['nPC']):
      for j in range(new_scores.shape[1]):
         if abs(new_scores[i][j]) > opts_dict['sigMul'] * (sigma_scores_gm[i]):
            comp_array[i][j] = 1
            eachruncount[j]=eachruncount[j]+1
         #Only check the first nPC number of scores, and sum comp_array together
         sum[i]=sum[i]+comp_array[i][j]   
       

   if len(ifiles) >= opts_dict['minRunFail']:
     num_run_less = False
   else:
     num_run_less = True
   #Check to see if sum is larger than min_run_fail, if so save the index of the sum
   for i in range(opts_dict['nPC']):
      if sum[i] >= opts_dict['minRunFail']:
        totalcount=totalcount+1
        sum_index.append(i+1)

   false_positive=check_falsepositive(opts_dict,sum_index)
   
   #If the length of sum_index is larger than min_PC_fail, the three runs failed
   if len(sum_index) >= opts_dict['minPCFail']:
      decision='FAILED'
   else:
      decision='PASSED'
   if num_run_less == False:
     print ' '
     print "Summary: "+str(totalcount)+" PC scores failed at least "+str(opts_dict['minRunFail'])+" runs: ",sum_index 
     print ' '
     print 'These runs '+decision+' according to our testing criterion.'
     if decision == 'FAILED' and false_positive != 1.0:
       print 'The probability of this test failing although everything functions correctly (false positive) is '+'{0:5.2f}'.format(false_positive*100)+'%.'
     print ' '
     print ' '
   else:
     print ' '
     print 'The number of run files is less than minRunFail (=2), so we cannot determin an overall pass or fail.'
     print ' ' 

   #Record the histogram of comp_array which value is one by the PCA scores
   for i in range(opts_dict['nPC']):
       index_list=[]
       for j in range(comp_array.shape[1]):
           if comp_array[i][j] == 1:
              index_list.append(j+1)
       if len(index_list) > 0:
          print "PC "+str(i+1)+": failed "+str(len(index_list))+" runs ",index_list
   print ' '

   #Record the index of comp_array which value is one
   for j in range(comp_array.shape[1]):
       index_list=[]
       for i in range(opts_dict['nPC']):
           if comp_array[i][j] == 1:
              index_list.append(i+1)
       print "Run "+str(j+1)+": "+str(eachruncount[j])+" PC scores failed ",index_list

   print ' '
#
# Command options for pyCECT.py
#
def CECT_usage():
    print '\n Compare test runs to an ensemble summary file. \n'
    print '  ----------------------------'
    print '   Args for pyTest :'
    print '  ----------------------------'
    print '   pyTest.py'
    print '   -h                      : prints out this usage message'
    print '   --verbose               : prints out in verbose mode (turned off by default)'
    print '   --printVarTest          : print out variable comparisons to RMSZ and global means (turned off by default)'
    print '   --sumfile  <ifile>      : the ensemble summary file (generated by pyEnsSum.py)'
    print '   --indir    <path>       : directory containing the input run files (at least 3 files)'
    print '   --timeslice <num>       : which time slice to use from input run files (default = 1)'
    print '   --nPC <num>             : number of principal components (PCs) to check (default = 50, but can\'t be greater than the number of variables)'
    print '   --sigMul   <num>        : number of standard deviations away from the mean defining the "acceptance region" (default = 2)'
    print '   --minPCFail <num>       : minimum number of PCs that must fail the specified number of runs for a FAILURE (default = 3)'
    print '   --minRunFail <num>      : minimum number of runs that <minPCfail> PCs must fail for a FAILURE (default = 2)'
    print '   --numRunFile <num>      : total number of runs to include in test (default = 3)'

    print 'Version 1.0.0'

#
# Command options for pyEnsSum.py
#
def EnsSum_usage(): 
    print '\n Creates the summary file for an ensemble of CAM data. \n'
    print '  ------------------------' 
    print '   Args for pyEnsSum : '
    print '  ------------------------' 
    print '   pyStats.py'
    print '   -h                   : prints out this usage message'
    print '   --verbose            : prints out in verbose mode (turned off by default)'
    print '   --sumfile    <ofile> : the output summary data file (default = ens.summary.nc)'
    print '   --indir      <path>  : directory containing all of the ensemble runs (default = ./)'
    print '   --esize  <num>       : Number of ensemble members (default = 151)'
    print '   --tag <name>         : Tag name used in metadata (default = cesm1_2_0)'
    print '   --compset <name>     : Compset used in metadata (default = FC5)'
    print '   --res <name>         : Resolution (used in metadata), (default = ne30_ne30)'
    print '   --tslice <num>       : the index into the time dimension, (default = 1)'
    print '   --mach <num>         : Machine name used in the metadata, (default = yellowstone)'
    print '   --jsonfile <fname>   : Jsonfile to provide that a list of variables that will be excluded  (no default)'
    print '   --mpi_enable         : Enable mpi mode if True'
    print '   --maxnorm            : Enable to generate max norm ensemble files'
    print '   --gmonly             : Only generate global_mean and PCA loadings (omit RMSZ information)'
    print '   '
    print 'Version 1.0.0'


#
# Random pick up three files out of a lot files
#
def Random_pickup(ifiles,opts_dict):
    if len(ifiles) > opts_dict['minRunFail']:
      random_index=random.sample(range(0,len(ifiles)),opts_dict['numRunFile'])
    else:
      random_index=range(len(ifiles))
    new_ifiles=[]
    for i in random_index:
       new_ifiles.append(ifiles[i])
       print ifiles[i]


    return new_ifiles

#
# Check the false positive rate
#
def check_falsepositive(opts_dict,sum_index):

    fp=np.zeros((opts_dict['nPC'],),dtype=np.float32)
    fp[0]=0.30305
    fp[1]=0.05069
    fp[2]=0.005745
    fp[3]=0.000435
    fp[4]=5.0e-05
    nPC = 50
    sigMul = 2
    minPCFail = 3
    minRunFail = 2
    numRunFile = 3

    if (nPC == opts_dict['nPC']) and (sigMul == opts_dict['sigMul']) and (minPCFail == opts_dict['minPCFail']) and (minRunFail == opts_dict['minRunFail']) and (numRunFile == opts_dict['numRunFile']):
       false_positive=fp[len(sum_index)-1]
    else:
       false_positive=1.0

    return false_positive
