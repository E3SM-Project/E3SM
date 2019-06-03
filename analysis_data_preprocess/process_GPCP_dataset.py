
#!/usr/bin/env cdat
"""
This script takes the GPCP dataset with monthly data on it from
"http://www.esrl.noaa.gov/psd/data/gridded/data.gpcp.html

This script then creates a netcdf file with the same data but that also fulfills
the cf-compliance and provenance requirements of ACME.
UVCDAT needs to be downloaded onto the system to run this.
Author: Chris Terai (terai1@llnl.gov) 2015-08-03
Modified for GPCP v2.3 by Jill Zhang (zhang40@llnl.gov) 2018-08-15
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


#=============================================
# I/O location and names
#input_hostpath='/work/terai1/ACME/GPCP/'    #!!! change directory to wherever you have downloaded GPCP data
data_name = 'GPCP'
input_hostpath='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/'+data_name+'/original_data/'    #!!! change directory to wherever you have downloaded GPCP data
input_filename='precip.mon.mean.nc'
input_hostANDfilename="".join([input_hostpath, input_filename])
#output_hostpath='/work/terai1/ACME/GPCP/'    #!!! change directory to wherever you want to save processed data
output_hostpath='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/'+data_name+'/time_series/'    #!!! change directory to wherever you want to save processed data
if data_name == 'GPCP_v2.2':
    output_filename = 'PRECT_197901_201412.nc'
elif data_name == 'GPCP':
    output_filename='PRECT_197901_201712.nc'
#==============================================
# Git has script location
#GIT_HASH = "~/ACME/PreAndPostProcessingScripts/utils/add_git_hash_to_netcdf_metadata "
#==============================================


# Input data
f_in=cdm.open(input_hostANDfilename)
varlist=['precip']
varidname=['PRECT']
varstandardname=['precipitation_flux']
data=[]
locs=[]


filecount = 0
fisc=varlist
f_index=0
for fi in fisc:
    print "##########################"
    print "".join(["** Processing: ",fi," **"])
    #GRAB DATA - make sure it works and print the shape of data
    try:
        if data_name == 'GPCP_v2.2':
            dattable=f_in(fi,time=("1978-12-31 23:59:59","2015-01-01 00:00:00"))
        elif data_name == 'GPCP':
            dattable=f_in(fi,time=("1978-12-31 23:59:59","2018-01-01 00:00:00"))
        print dattable.shape
    except:
        print "Couldn't access data for ",fi
        continue #go on to next var in loop.
    dattable.standard_name=varstandardname[f_index]
    dattable.id=varidname[f_index]
    dattable.units='mm/d'
    dattable.long_name='Average Monthly Rate of Precipitation'
    var_missing_value=dattable.missing_value
    
    if data_name =='GPCP_v2.2':
        #right side up the globe
        lat=dattable.getAxis(1)
        bnd=lat.getBounds()[::-1]
        lat[:]=lat[::-1]
        lat.setBounds(bnd)
        dattable=dattable[:,::-1,:]
        dattable.setAxis(1,lat)
    if filecount==0:
        outfile="".join([output_hostpath,output_filename])
        f_out=cdm.open(outfile,'w')
        
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
        f_out.history="".join([" CRT[",time_format,"]: added axes bounds to data using UVCDAT, added standard_name to variable, and appended provenance information:",f_out.history])
        f_out.data_provider = "NOAA Earth System Research Laboratory // NASA/Goddard Space Flight Center"
        f_out.Conventions = "CF-1.7"
#        f_out.references = "".join([f_out.references," // Adler, R.F., G.J. Huffman, A. Chang, R. Ferraro, P. Xie, J. Janowiak, B. Rudolf, U. Schneider, S. Curtis, D. Bolvin, A. Gruber, J. Susskind, and P. Arkin, 2003: The Version 2 Global Precipitation Climatology Project (GPCP) Monthly Precipitation Analysis (1979-Present). J. Hydrometeor., 4,1147-1167."])
        f_out.references = " Adler, R.F., G.J. Huffman, A. Chang, R. Ferraro, P. Xie, J. Janowiak, B. Rudolf, U. Schneider, S. Curtis, D. Bolvin, A. Gruber, J. Susskind, and P. Arkin, 2003: The Version 2 Global Precipitation Climatology Project (GPCP) Monthly Precipitation Analysis (1979-Present). J. Hydrometeor., 4,1147-1167."
        f_out.institution     = "NOAA Earth System Research Laboratory // Data processed at: Program for Climate Model Diagnosis and Intercomparison (LLNL)"
        f_out.contact    = "Physical Sciences Division: Data Management, NOAA/ESRL/PSD, esrl.psd.data@noaa.gov // Processed by: Chris Terai; terai1@llnl.gov; +1 925 422 8830"
        f_out.host            = "".join([gethostname(),'; UVCDAT version: ',".".join(["%s" % el for el in cdat_info.version()]),
                                           '; Python version: ',replace(replace(sys.version,'\n','; '),') ;',');')])
        
        
    
    f_out.write(dattable) ; 
    print "".join(["** Finished processing: ",fi," **"])
    filecount = filecount + 1; filecount_s = '%06d' % filecount
    f_index=f_index + 1

f_out.close()
print 'got data from '+str(filecount)+' variables.'
print "      adding _FillValue attribute"
call("".join(["ncatted -a _FillValue,PRECT,o,f,",str(var_missing_value)," ",outfile]), shell=True)
call("".join(["ncatted -a _FillValue,'^bounds',d,, ",outfile]), shell=True)

#Add git hash to obs
#call("".join([GIT_HASH,outfile]), shell=True) #the directories may differ depending on where you have saved PreAndPostProcessingScripts
