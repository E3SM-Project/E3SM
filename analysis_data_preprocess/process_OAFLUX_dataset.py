
#!/usr/bin/env cdat
"""
Based on Chris Terai's original script for GPCP
Modified for OAFlux by Jill Zhang (zhang40@llnl.gov) 2018-09-15
"""
import argparse,datetime,gc,re,sys,time
import cdat_info,cdtime,code,inspect,os,re,string,sys,pytz
import cdms2 as cdm
import MV2 as MV #stuff for dealing with masked values.
import cdutil as cdu
import glob
import os
from socket import gethostname
from string import replace
import numpy
from subprocess import call

# Set nc classic as outputs
cdm.setCompressionWarnings(0) ; # Suppress warnings
cdm.setNetcdfShuffleFlag(0)
cdm.setNetcdfDeflateFlag(1) ; # was 0 130717
cdm.setNetcdfDeflateLevelFlag(9) ; # was 0 130717
cdm.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated

userinitials = 'CZ'

#=============================================
# I/O location and names
#input_hostpath='/work/terai1/ACME/GPCP/'    #!!! change directory to wherever you have downloaded GPCP data
input_hostpath='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/OAFlux/original_data/'    #!!! change directory to wherever you have downloaded GPCP data
input_filename='OAFlux_MonthAvg_197901-201312.nc'
input_hostANDfilename="".join([input_hostpath, input_filename])
#output_hostpath='/work/terai1/ACME/GPCP/'    #!!! change directory to wherever you want to save processed data
output_hostpath='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/OAFlux/time_series/'    #!!! change directory to wherever you want to save processed data
output_filename='197901_201312.nc'

#==============================================
# Git has script location
#GIT_HASH = "~/ACME/PreAndPostProcessingScripts/utils/add_git_hash_to_netcdf_metadata "
#==============================================


# Input data
print input_hostANDfilename
f_in=cdm.open(input_hostANDfilename)


varlist=f_in.listvariables()

print(varlist)
print(f_in.listdimension())

ignore_variables=['bounds_lat','bounds_time','bounds_lon','bnds_lat','bnds_time','bnds_lon','date','lon_bnds','lat_bnds','time_bnds','err']  #variables to ignore when making climatology means

for ji in varlist:
    if ji not in ignore_variables:
    # Make changes to variable ID, standard_name, and units to make them fit cf-conventions
        if ji == 'evapr':
            outputdattable = f_in(ji)
            outputdattable=outputdattable/365.25/24./3600.*10.  #conversion from cm/yr to mm s-1
            outputdattable.id='QFLX'
            #outputdattable,outputdattable=reductions.reconcile_units(outputdattable,outputdattable,preferred_units='mm s-1')
            outputdattable.units='kg m-2 s-1'
            outputdattable.standard_name='water_evaporation_flux'
            print ji,' modified'
#        if ji == 'err':
#            outputdattable = f_in(ji)
#            outputdattable=outputdattable/365.25/24./3600.*10.  #conversion from cm/yr to mm s-1
#            outputdattable.id='QFLX_ERR'
#            #outputdattable,outputdattable=reductions.reconcile_units(outputdattable,outputdattable,preferred_units='mm s-1')
#            outputdattable.units='kg m-2 s-1'
#            outputdattable.standard_name='water_evaporation_flux_error'
#            #print ji,' modified'
        if ji == 'shtfl':
            outputdattable = f_in(ji)
            outputdattable.id='SHFLX'
            outputdattable.units='W m-2'
            outputdattable.standard_name='surface_upward_sensible_heat_flux'
            #print ji,' modified'
        if ji == 'lhtfl':
            outputdattable = f_in(ji)
            outputdattable.id='LHFLX'
            outputdattable.units='W m-2'
            outputdattable.standard_name='surface_upward_latent_heat_flux'
            #print ji,' modified'
        print ji
        f_out=cdm.open(output_hostpath+outputdattable.id+'_'+output_filename,'w')
        f_out.write(outputdattable)

        att_keys = f_in.attributes.keys()
        att_dic = {}
        for i in range(len(att_keys)):
            att_dic[i]=att_keys[i],f_in.attributes[att_keys[i]]
            to_out = att_dic[i]
            setattr(f_out,to_out[0],to_out[1])
        local                       = pytz.timezone("America/Los_Angeles")
        time_now                    = datetime.datetime.now();
        local_time_now              = time_now.replace(tzinfo = local)
        utc_time_now                = local_time_now.astimezone(pytz.utc)
        time_format                 = utc_time_now.strftime("%d-%b-%Y %H:%M:%S %p")
        f_out.history="".join([userinitials," [",time_format,"]: changed variable id, units, and standard name to fit cf-compliance y"])
        f_out.contact="For original data: Dr. Lisan Yu (lyu@whoi.edu) and Dr. Robert A. Weller (rweller@whoi.edu) // For processed data: Chris Terai (terai1@llnl.gov)"
        f_out.references="Yu, L., X. Jin, and R. A. Weller, 2008: Multidecade Global Flux Datasets from the Objectively Analyzed Air-sea Fluxes (OAFlux) Project: Latent and sensible heat fluxes, ocean evaporation, and related surface meteorological variables. Woods Hole Oceanographic Institution, OAFlux Project Technical Report. OA-2008-01, 64pp. Woods Hole. Massachusetts"
        f_out.title="OAFlux LHF, SHF, and Evaporation rate 1980-01 to 2004-12"
        f_out.institution="Woods Hole Oceanographic Institution, Woods Hole, Massachusetts"
        f_out.close()



#Add git hash to obs
#call("".join([GIT_HASH,outfile]), shell=True) #the directories may differ depending on where you have saved PreAndPostProcessingScripts
