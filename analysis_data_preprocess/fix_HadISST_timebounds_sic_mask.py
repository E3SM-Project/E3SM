import cdms2
import cdutil
import MV2
#Data downloaded from here https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
#This script is to fix overlapping timebounds from original data
#also add ice mask

# Set nc classic as outputs
cdms2.setCompressionWarnings(0) ; # Suppress warnings
cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(1) ; # was 0 130717
cdms2.setNetcdfDeflateLevelFlag(9) ; # was 0 130717
cdms2.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated

filepath = '/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/HadISST/original_data/'
filename1 = 'HadISST_ice.nc'
filename2 = 'HadISST_sst.nc'
fin1 = cdms2.open(filepath +filename1)
fin2 = cdms2.open(filepath +filename2)
ice = fin1('sic')
sst = fin2('sst')

fout = cdms2.open(filepath + 'HadISST_sst_ice_masked.nc','w')
sst_masked = MV2.masked_where(ice>0,sst,copy=True)
sst_masked.id = 'sst'
cdutil.setTimeBoundsMonthly(sst_masked)
#reverse latitude so that latitude in ascending
sst_masked = sst_masked[:,::-1,:]
fout.write(sst_masked)

att_keys = fin2.attributes.keys()
att_dic = {}
for i in range(len(att_keys)):
    att_dic[i]=att_keys[i],fin2.attributes[att_keys[i]]
    to_out = att_dic[i]
    
    setattr(fout,to_out[0],to_out[1])
print fout.attributes
fout.close()
