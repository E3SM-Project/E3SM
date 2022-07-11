import cdms2
from pathlib import Path
import MV2

datapath='/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/TRMM-3B43v-7_3hr/'
for path in Path(datapath).rglob('*.nc'):
    print(path.name)
    print(path)
    filename = path.name.split("/")[-1]
    print(filename)

    #filename ='TRMM-3B43v-7_3hr_ANN_199801_201312_climo.nc'
    f_in = cdms2.open(datapath+filename)
    var0 = f_in('pr')
    var = f_in('pr')
    lat = var0.getLatitude()
    lon = var0.getLongitude()
    nlat = len(lat)
    var.id = 'pr'
    # Fix units problem for obs4mip
    #var = var / 3600.0 * 1000 / 1000.0

    # Manually shift time by 1.5 hour to work around an ncclimo problem:
    # The shift may be due to the assumption that ncclimo makes that times are encoded using the CESM/E3SM convention where times refer to the end of the timestep
    print('var.time',var.getTime()[:])
    for ilat in range(nlat):
        print,'ilat= ',ilat
        var[0,ilat,:] = var0[7,ilat,:]
        var[1:8,ilat,:] = var0[0:7,ilat,:]
    newtime=cdms2.createAxis(MV2.zeros(8))
    newtime.id="time" #name of dimension
    newtime.designateTime()  # tell cdms to add attributes that make it time
    newtime.units = var0.getTime().units
    newtime[:] = [x-0.0625 for x in var0.getTime()[:]]
    var.setAxis(0,newtime)
    print(newtime)
    print(var.getTime())


    f_out = cdms2.open(filename,'w')
    f_out.write(var)
    att_keys = list(f_in.attributes.keys())
    att_dic = {}
    for i in range(len(att_keys)):
        att_dic[i]=att_keys[i],f_in.attributes[att_keys[i]]
        to_out = att_dic[i]
        setattr(f_out,to_out[0],to_out[1])
    f_in.close()
    f_out.close()
