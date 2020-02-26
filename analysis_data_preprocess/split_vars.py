import cdms2
import pytz,time,datetime

# Set nc classic as outputs
cdms2.setCompressionWarnings(0) ; # Suppress warnings
cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(1) ; # was 0 130717
cdms2.setNetcdfDeflateLevelFlag(9) ; # was 0 130717
cdms2.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated

#apply to 4.0_toa, 2.8_toa, 4.0_surface, 2.8_toa, 4.1_toa, 4.1_surface
dataset = '4.1_surface'
data = 'ceres_ebaf_surface_v4.1'
userinitials = "C.Zhang, zhang40@llnl.gov"
if dataset in ['4.1_toa','4.1_surface']:
    years = '200101_201812'
else:
    years = '200101_201512'
datapath = '/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/CERES-EBAF/'+dataset+'/time_series/'
filename = datapath + data + '_' + years +'.nc'
print(filename ) 
fin = cdms2.open(filename)
vars = fin.listvariable()
print("variable list: ", vars)
for var in vars:
    variable = fin(var)
    if(var.find('_clr_t')>=0):
       new_var = var.replace('_clr_t','_clr')
       variable.id = new_var
       var = new_var
       print(new_var)

    fout = cdms2.open(datapath+var+'_'+years+'.nc', 'w')
    fout.write(variable)
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
    

