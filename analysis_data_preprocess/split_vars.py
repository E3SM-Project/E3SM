import cdms2
import pytz,time,datetime

# Set nc classic as outputs
cdms2.setCompressionWarnings(0) ; # Suppress warnings
cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(1) ; # was 0 130717
cdms2.setNetcdfDeflateLevelFlag(9) ; # was 0 130717
cdms2.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated

#apply to 4.0_toa, 2.8_toa, 4.0_surface, 2.8_toa
datapath = '/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/CERES-EBAF/4.0_surface/time_series/'
data = 'ceres_ebaf_surface_v4.0'
years = '200101_201512'
filename = datapath + data + '_' + years +'.nc'
print filename  
fin = cdms2.open(filename)
vars = fin.listvariable()
print "variable list: ", vars

userinitials = "C.Zhang, zhang40@llnl.gov"
for var in vars:
    print var
    fout = cdms2.open(datapath+var+'_'+years+'.nc', 'w')
    fout.write(fin(var))
    local                       = pytz.timezone("America/Los_Angeles")
    time_now                    = datetime.datetime.now();
    local_time_now              = time_now.replace(tzinfo = local)
    utc_time_now                = local_time_now.astimezone(pytz.utc)
    time_format                 = utc_time_now.strftime("%d-%b-%Y %H:%M:%S %p")
    fout.history="".join([userinitials," [",time_format,"]: split variables from multi-var files"])
    fout.title = fin.title   
    try:
        fout.Version = fin.Version   
    except: 
        fout.Version = fin.version
    fout.institution = fin.institution   
    fout.close()
    

