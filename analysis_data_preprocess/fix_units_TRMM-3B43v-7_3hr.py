from pathlib import Path

import cdms2

datapath = (
    "/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/TRMM/climatology_diurnal_cycle/"
)
for path in Path(datapath).rglob("*.nc"):
    print(path.name)
    print(path)
    filename = path.name.split("/")[-1]
    print(filename)

    # filename ='TRMM-3B43v-7_3hr_ANN_199801_201312_climo.nc'
    f_in = cdms2.open(datapath + filename)
    var = f_in("pr")
    var = var / 3600.0 * 1000 / 1000.0
    var.id = "pr"
    f_out = cdms2.open(datapath + "units_fix_" + filename, "w")
    f_out.write(var)

    att_keys = list(f_in.attributes.keys())
    att_dic = {}
    for i in range(len(att_keys)):
        att_dic[i] = att_keys[i], f_in.attributes[att_keys[i]]
        to_out = att_dic[i]
        setattr(f_out, to_out[0], to_out[1])
    print(var.mean())
    f_in.close()
    f_out.close()
