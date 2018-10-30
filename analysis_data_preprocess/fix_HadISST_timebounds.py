import cdms2
import cdutil

#This script is to fix overlapping timebounds from original data

# Set nc classic as outputs
cdms2.setCompressionWarnings(0) ; # Suppress warnings
cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(1) ; # was 0 130717
cdms2.setNetcdfDeflateLevelFlag(9) ; # was 0 130717
cdms2.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated

filepath = '/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/HadISST/time_series/'
filename = 'sst_187001_201612_orig.nc'
fin = cdms2.open(filepath +filename)
var = fin('sst')
cdutil.setTimeBoundsMonthly(var)
outfilename = 'sst_187001_201612.nc'
fout = cdms2.open(filepath + outfilename,'w')
fout.write(var)

att_keys = fin.attributes.keys()
att_dic = {}
for i in range(len(att_keys)):
    att_dic[i]=att_keys[i],fin.attributes[att_keys[i]]
    to_out = att_dic[i]
    
    setattr(fout,to_out[0],to_out[1])
print fout.attributes
fout.close()
