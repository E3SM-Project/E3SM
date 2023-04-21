import e3sm_diags
from e3sm_diags.driver import utils
import cdms2
import cdutil
import numpy
import os
import json

import csv

from e3sm_diags.parameter.core_parameter import CoreParameter


def global_integral(var, area_m2):
    """ Compute global integral of 2 dimentional properties"""
    return numpy.sum(numpy.sum(abs(var)*area_m2,axis = 0), axis=0)

def calc_column_integral(data, aerosol, season):
    """ Calculate column integrated mass """
    mass = data.get_climo_variable(f'Mass_{aerosol}', season)
    hyai, hybi, ps = data.get_extra_variables_only(
         f'Mass_{aerosol}', season, extra_vars=["hyai", "hybi", "PS"]
     )

    p0 = 100000.0  # Pa
    ps = ps   # Pa
    pressure_levs = cdutil.vertical.reconstructPressureFromHybrid(ps, hyai, hybi, p0)

    #(72,lat,lon)
    delta_p = numpy.diff(pressure_levs,axis = 0)
    mass_3d = mass*delta_p/9.8 #mass density * mass air   kg/m2
    burden = numpy.nansum(mass_3d,axis = 0)   #kg/m2
    return burden
    
def generate_metrics_dic(data, aerosol, season):
    metrics_dict = {}
    wetdep = data.get_climo_variable(f'{aerosol}_SFWET', season)
    drydep = data.get_climo_variable(f'{aerosol}_DDF', season)
    srfemis = data.get_climo_variable(f'SF{aerosol}', season)
    area = data.get_extra_variables_only(
                f'{aerosol}_DDF', season, extra_vars=["area"]
            )
    area_m2 = area * REARTH**2

    burden = calc_column_integral(data, aerosol, season)
    burden_total= global_integral(burden, area_m2)*1e-9 # kg to Tg
    print(f'{aerosol} Burden (Tg): ',f'{burden_total:.3f}')
    sink = global_integral((drydep-wetdep),area_m2)*UNITS_CONV
    drydep = global_integral(drydep,area_m2)*UNITS_CONV
    wetdep = global_integral(wetdep,area_m2)*UNITS_CONV
    srfemis = global_integral(srfemis,area_m2)*UNITS_CONV
    print(f'{aerosol} Sink (Tg/year): ',f'{sink:.3f}')
    print(f'{aerosol} Lifetime (days): ',f'{burden_total/sink*365:.3f}')
    metrics_dict = {
    "Surface Emission (Tg/yr)": f'{srfemis:.3f}',
    "Sink (Tg/yr)": f'{sink:.3f}',
    "Dry Deposition (Tg/yr)": f'{drydep:.3f}',
    "Wet Deposition (Tg/yr)": f'{wetdep:.3f}',
    "Burden (Tg)": f'{burden_total:.3f}',
    "Lifetime (Days)": f'{burden_total/sink*365:.3f}',
    }
    return metrics_dict

param = CoreParameter()
param.test_name = 'v2.LR.historical_0101'
param.test_name = 'F2010.PD.NGD_v3atm.0096484.compy'
param.test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
param.test_data_path = '/compyfs/mahf708/E3SMv3_dev/F2010.PD.NGD_v3atm.0096484.compy/post/atm/180x360_aave/clim/10yr'
test_data = utils.dataset.Dataset(param, test=True)

#rearth = 6.37122e6 #km
#UNITS_CONV = 86400.0*365.0*1e-9 # kg/s to Tg/yr
REARTH = 6.37122e6 #km
UNITS_CONV = 86400.0*365.0*1e-9 # kg/s to Tg/yr
# TODO:
# Convert so4 unit to TgS
#mwso4 = 115.0
#mws = 32.066
#UNITS_CONV_S = UNITS_CONV/mwso4*mws # kg/s to TgS/yr


species = ["bc", "dst", "mom", "ncl","pom","so4","soa"]
SPECIES_NAMES = {"bc": "Black Carbon",
    "dst": "Dust",
    "mom": "Marine Organic Matter",
    "ncl": "Sea Salt",
    "pom": "Primary Organic Matter",
    "so4": "Sulfate",
    "soa": "Secondary Organic Aerosol"}
MISSING_VALUE = 999.999
metrics_dict = {}
metrics_dict_ref = {}

seasons = ["ANN"]

ref_data_path = os.path.join(
        e3sm_diags.INSTALL_PATH,
        "control_runs",
        "aerosol_global_metrics_benchmarks.json",
    )

with open(ref_data_path, 'r') as myfile:
    ref_file=myfile.read()

metrics_ref = json.loads(ref_file)

for season in seasons:
    for aerosol in species:
        print(f'Aerosol species: {aerosol}')
        metrics_dict[aerosol] = generate_metrics_dic(test_data, aerosol, season)
        metrics_dict_ref[aerosol] = metrics_ref[aerosol]
        #metrics_dict_ref[aerosol] = {
        #    "Surface Emission (Tg/yr)": f'{MISSING_VALUE:.3f}',
        #    "Sink (Tg/yr)": f'{MISSING_VALUE:.3f}',
        #    "Dry Deposition (Tg/yr)": f'{MISSING_VALUE:.3f}',
        #    "Wet Deposition (Tg/yr)": f'{MISSING_VALUE:.3f}',
        #    "Burden (Tg)": f'{MISSING_VALUE:.3f}',
        #    "Lifetime (Days)": f'{MISSING_VALUE:.3f}',
        #    }
    
        with open(f'aerosol_table_{season}.csv', "w") as table_csv:
            writer = csv.writer(
                table_csv,
                delimiter=",",
                quotechar="'",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator='\n',
            )
            #writer.writerow([" ", "test","ref",])
            for key, values in metrics_dict.items():
                writer.writerow([SPECIES_NAMES[key]])
                print('key',key, values)
                for value in values:
                    print(value)
                    line = []
                    line.append(value)
                    line.append(values[value])
                    line.append(metrics_dict_ref[key][value])
                    print(line, 'line')
                    writer.writerows([line])
                writer.writerows([""])




