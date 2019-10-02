#!/usr/bin/env cdat
"""
This script takes an observational dataset that contains monthly data 
and creates the appropriate climo files. Aspects of this code have been taken 
from Susannah Burrows's ACME_obs_for_diagnostics_metadata_template.sh
See: https://gist.github.com/susburrows/fb0fe95ad2a3348d8ce1#file-acme_obs_for_diagnostics_metadata_template-sh

IMPORTANT: The script assumes that the observational dataset already 
contains all the necessary attributes (correct units, standard_name, etc).

WHAT THIS DOES: The script creates 17 climo files: one for each month of the 
year (12), four seasonal datasets (4), and one annual mean dataset.

For DJF, the dataset only calculates means over contiguous DJF seasons, so the 
first JF and last D are not included in the mean. - This should be reflected in 
the climatology_bnds attribute.

In order to run this script, UVCDAT (2.2.0 or later) and NCO (4.5.0 or later) 
need to be installed on your system. 

Author: Chris Terai (terai1@llnl.gov) 


Modified by Jill Zhang (zhang40@llnl.gov) 08-27-2018
"""

# Upload the following libraries to call in this script
import argparse,datetime,gc,re,sys,time
import cdat_info,cdtime,code,inspect,os,re,string,sys,pytz
import cdms2 as cdm
import MV2 as MV     #functions for dealing with masked values.
import cdutil as cdu
import glob
import os
from socket import gethostname
#from string import replace
import numpy
#import metrics.computation.reductions as reductions       # primarily for reconcile units
from subprocess import call  #for calling NCO functions

# Set nc classic as outputs
cdm.setCompressionWarnings(0) ; # Suppress warnings
cdm.setNetcdfShuffleFlag(0)
cdm.setNetcdfDeflateFlag(1) ; # was 0 130717
cdm.setNetcdfDeflateLevelFlag(9) ; # was 0 130717
cdm.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated


#^^^^^^^^^^^^^^^   MODIFY FOLLOWING THINGS DEPENDING ON THE DATASET ^^^^^^^^^^^^^^^
#userinitials="CRT"
userinitials="CZ"


#input_hostpath="/work/terai1/ACME/OAFLUX/"               #location of input data
#output_hostpath="/work/terai1/ACME/OAFLUX/" #location to save climatology files
input_hostpath="/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/OAFlux/"               #location of input data
output_hostpath="".join([input_hostpath,'climatology/']) #location to save climatology files
input_filename="OAFlux_MonthAvg_197901-201312.nc"                           #orig file name

input_hostANDfilename="".join([input_hostpath,'original_data/', input_filename])        #string of location+filename of input file

rootname='OAFlux' #name of the climo files
startyear="1980"  #start year of the dataset
startmonth='01'   #start month of the dataset, typically 01
endyear="2005"    #end year of the dataset
endmonth='12'     #end month of the dataset, typically 12
#^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^



#for running climatology.py
fisc=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC','MAM','JJA','SON','DJF','YEAR']
#for naming files
fisc2=['01','02','03','04','05','06','07','08','09','10','11','12','MAM','JJA','SON','DJF','ANN']
#fisc=['JAN']
clim_edges=["-1-1","-2-1","-3-1","-4-1","-5-1","-6-1","-7-1","-8-1","-9-1","-10-1","-11-1","-12-1"]
clim_edge_values=numpy.zeros((17,2),dtype=numpy.int)
#JAN
clim_edge_values[0,0]=0 
clim_edge_values[0,1]=1
#FEB
clim_edge_values[1,0]=1
clim_edge_values[1,1]=2
#MAR
clim_edge_values[2,0]=2 
clim_edge_values[2,1]=3
#APR
clim_edge_values[3,0]=3 
clim_edge_values[3,1]=4
#MAY
clim_edge_values[4,0]=4 
clim_edge_values[4,1]=5
#JUN
clim_edge_values[5,0]=5 
clim_edge_values[5,1]=6
#JUL
clim_edge_values[6,0]=6 
clim_edge_values[6,1]=7
#AUG
clim_edge_values[7,0]=7 
clim_edge_values[7,1]=8
#SEP
clim_edge_values[8,0]=8 
clim_edge_values[8,1]=9
#OCT
clim_edge_values[9,0]=9 
clim_edge_values[9,1]=10
#NOV
clim_edge_values[10,0]=10
clim_edge_values[10,1]=11
#DEC
clim_edge_values[11,0]=11 
clim_edge_values[11,1]=0
#MAM
clim_edge_values[12,0]=2 
clim_edge_values[12,1]=5
#JJA
clim_edge_values[13,0]=5
clim_edge_values[13,1]=8
#SON
clim_edge_values[14,0]=8
clim_edge_values[14,1]=11
#DJF
clim_edge_values[15,0]=11
clim_edge_values[15,1]=2
#YEAR
clim_edge_values[16,0]=0
clim_edge_values[16,1]=0

f_in=cdm.open(input_hostANDfilename)
varlist=f_in.listvariables()

print(varlist)
print(f_in.listdimension())

ignore_variables=['bounds_lat','bounds_time','bounds_lon','bnds_lat','bnds_time','bnds_lon','date','lon_bnds','lat_bnds','time_bnds']  #variables to ignore when making climatology means

f_index=0
#loop through the different seasons/months to make climatologies
for fi in fisc:
    print(fi)
    fi2=fisc2[f_index]
    print(fi2)
    #output_filename="".join([rootname,'_',fi2,'_',startyear,startmonth,'-',endyear,endmonth,'_climo.nc'])  #old option of including start and end months
    output_filename="".join([rootname,'_',fi2,'_climo.nc'])
    output_hostANDfilename="".join([output_hostpath, output_filename])
    f_out=cdm.open(output_hostANDfilename,'w')
    inputtime=f_in.getAxis('time')
    for ji in varlist:
        if ji not in ignore_variables:
            print(ji)
            dattable=f_in(ji)
            print("evaluate")
            print(".".join(['cdu',fi,'climatology(dattable, criteriaarg = [0.99,None])']))
            outputdattable=eval(".".join(['cdu',fi,'climatology(dattable, criteriaarg = [0.99,None])'])) #only full DJF ranges are included in the calculation
            outputdattable.cell_methods="time: mean over years"
            
            #specify the time axis long_name and the units - uses the startyear-1-1 as the starting point
            time_axis=outputdattable.getAxis(0)
            time_axis.climatology="climatology_bnds"
            time_axis.long_name="time"
            time_axis.units="".join(["days since ",startyear,"-1-1"])
            outputdattable.setAxis(0,time_axis)
            
            
            
            var_missing_value=outputdattable.missing_value  #obtain missing value for use later in appending _FillValue
#            #Write out the attributes of the variable (dattable) onto the climatological-mean variable (outputdattable)
#            att_keys = dattable.attributes.keys()
#            att_dic = {}
#            dattable.attributes.keys()
#            for i in range(len(att_keys)):
#                print(att_keys[i])
#                att_dic[i]=att_keys[i],dattable.attributes[att_keys[i]]
#                to_out = att_dic[i]
#                setattr(outputdattable,to_out[0],to_out[1])
#            time_axis=outputdattable.getAxis(0)
#            time_axis.climatology="climatology_bnds"
#            outputdattable.setAxis(0,time_axis)
            
            # Make changes to variable ID, standard_name, and units to make them fit cf-conventions
            if ji == 'evapr':
                outputdattable=outputdattable/365.25/24./3600.*10.  #conversion from cm/yr to mm s-1
                outputdattable.id='QFLX'
                #outputdattable,outputdattable=reductions.reconcile_units(outputdattable,outputdattable,preferred_units='mm s-1')
                outputdattable.units='kg m-2 s-1'
                outputdattable.standard_name='water_evaporation_flux'
                #print ji,' modified'
            if ji == 'err':
                outputdattable=outputdattable/365.25/24./3600.*10.  #conversion from cm/yr to mm s-1
                outputdattable.id='QFLX_ERR'
                #outputdattable,outputdattable=reductions.reconcile_units(outputdattable,outputdattable,preferred_units='mm s-1')
                outputdattable.units='kg m-2 s-1'
                outputdattable.standard_name='water_evaporation_flux_error'
                #print ji,' modified'
            if ji == 'shtfl':
                outputdattable.id='SHFLX'
                outputdattable.units='W m-2'
                outputdattable.standard_name='surface_upward_sensible_heat_flux'
                #print ji,' modified'
            if ji == 'lhtfl':
                outputdattable.id='LHFLX'
                outputdattable.units='W m-2'
                outputdattable.standard_name='surface_upward_latent_heat_flux'
                #print ji,' modified'
            f_out.write(outputdattable)
    
    time_axis=f_out.getAxis('time')
    bnds_axis=f_out.getAxis('bound')
    climatology=None
    climo_start="".join([startyear,clim_edges[clim_edge_values[f_index,0]]])
    climo_end="".join([endyear,clim_edges[clim_edge_values[f_index,1]]])
    #climatology_bnds_value={climo_start, climo_end}
    
    start_date_str=clim_edges[clim_edge_values[f_index,0]]
    end_date_str=clim_edges[clim_edge_values[f_index,1]]
    nothing,start_month,startdate=start_date_str.split("-")
    nothing,end_month,enddate=end_date_str.split("-")
    cstart=cdtime.comptime(int(startyear),int(start_month),int(startdate))
    cend=cdtime.comptime(int(endyear),int(end_month),int(enddate))
    if f_index==11:
        cend=cdtime.comptime(int(endyear)+1,int(end_month),int(enddate))
    if f_index==16:
        cend=cdtime.comptime(int(endyear)+1,int(end_month),int(enddate))
    rstart=cstart.torel("".join(["days since ",startyear,"-1-1"]))  #convert startdate string into days since
    rend=cend.torel("".join(["days since ",startyear,"-1-1"]))  #convert enddate string into days since
    climatology_bnds_value=numpy.zeros([1,2])
    #daysSinceStart,day_str,since_str,date_str=rstart.split(" ")
    #daysSinceEnd,day_str,since_str,date_str=rend.split(" ")
    climatology_bnds_value[0,0]=int(rstart.value)
    climatology_bnds_value[0,1]=int(rend.value)
    climatology_bnds=cdm.createVariable(climatology_bnds_value,id="climatology_bnds",fill_value=None)
    climatology_bnds.units="".join(["days since ",startyear,"-1-1"])
    climatology_bnds.setAxis(0,time_axis)
    climatology_bnds.setAxis(1,bnds_axis)
    f_out.write(climatology_bnds)
    
#    att_keys = f_in.attributes.keys()
#    att_dic = {}
#    for i in range(len(att_keys)):
#        att_dic[i]=att_keys[i],f_in.attributes[att_keys[i]]
#        to_out = att_dic[i]
#        setattr(f_out,to_out[0],to_out[1])
    setattr(f_out,'season',fi)
    year_str="".join([startyear,'-',endyear])
    setattr(f_out,'yrs_averaged',year_str)
    setattr(f_out,'data_name','OAFlux')
    local                       = pytz.timezone("America/Los_Angeles")
    time_now                    = datetime.datetime.now();
    local_time_now              = time_now.replace(tzinfo = local)
    utc_time_now                = local_time_now.astimezone(pytz.utc)
    time_format                 = utc_time_now.strftime("%d-%b-%Y %H:%M:%S %p")
    f_out.history="".join([userinitials," [",time_format,"]: changed variable id, units, and standard name to fit cf-compliance and created time-mean climatology file using climatology.py"])
    f_out.contact="For original data: Dr. Lisan Yu (lyu@whoi.edu) and Dr. Robert A. Weller (rweller@whoi.edu) // For processed data: Chris Terai (terai1@llnl.gov)"
    f_out.references="Yu, L., X. Jin, and R. A. Weller, 2008: Multidecade Global Flux Datasets from the Objectively Analyzed Air-sea Fluxes (OAFlux) Project: Latent and sensible heat fluxes, ocean evaporation, and related surface meteorological variables. Woods Hole Oceanographic Institution, OAFlux Project Technical Report. OA-2008-01, 64pp. Woods Hole. Massachusetts"
    f_out.title="OAFlux LHF, SHF, and Evaporation rate 1980-01 to 2004-12"
    f_out.institution="Woods Hole Oceanographic Institution, Woods Hole, Massachusetts"
    f_index=f_index+1
    f_out.close()
    
    #introduce _FillValue attribute using NCO tools and remove missing_value from other variables
    #print "      adding _FillValue attribute"
    #print "".join(["       ncatted -a _FillValue,,o,d,",str(var_missing_value)," ",output_hostANDfilename])
    #call("".join(["ncatted -a _FillValue,,o,d,",str(var_missing_value)," ",output_hostANDfilename]), shell=True)
    #print "      removing _FillValue and missing_value from climatology_bnds attribute"
    call("".join(["ncatted -a missing_value,climatology_bnds,d,, ",output_hostANDfilename]), shell=True)
    call("".join(["ncatted -a _FillValue,climatology_bnds,d,, ",output_hostANDfilename]), shell=True)
    call("".join(["ncatted -a _FillValue,'^bounds',d,, ",output_hostANDfilename]), shell=True)
    call("".join(["ncatted -a _FillValue,'lat',d,, ",output_hostANDfilename]), shell=True)
    call("".join(["ncatted -a _FillValue,'lon',d,, ",output_hostANDfilename]), shell=True)
    call("".join(["ncatted -a _FillValue,'time',d,, ",output_hostANDfilename]), shell=True)
    
    #Add git hash to climo file
#    call("".join(["~/ACME/PreAndPostProcessingScripts/utils/add_git_hash_to_netcdf_metadata ",output_hostANDfilename]), shell=True)
