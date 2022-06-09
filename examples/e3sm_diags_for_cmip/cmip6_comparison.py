# Script to perform comparison with CMIP6 models

# On NERSC
"""
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
"""

import glob
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.transforms as mtransforms
import numpy as np
import numpy.ma as ma
import os

# --- Function to read E3SM Diags metrics for CMIP6 models ---
def read_e3sm_diags_metrics(path, variables, seasons, names=None):

  # List of available models
  models = []
  paths = []
  dirs = sorted( glob.glob(path + os.path.sep) )
  for d in dirs:
      tmp = d.split(os.path.sep)
      model = tmp[-6]
      paths.append(d)
      models.append(model)
  if names:
      models = names

  # Array to hold data
  nmodels = len(models)
  nvariables = len(variables)
  nseasons = len(seasons)
  data = ma.array(np.zeros((nmodels,nvariables,nseasons)), mask=True)

  # Fill data
  for imodel in range(nmodels):
    for iseason in range(nseasons):
      # Open metrics file
      fname = paths[imodel]+'/%s_metrics_table.csv' % (seasons[iseason])
      with open(fname, 'r') as f:
          content = f.readlines()
          for ivariable in range(nvariables):
              # Skip of model has been flagged for this variable
              if models[imodel] in variables[ivariable]['exclude']:
                print("Excluding: %s, %s, %s" % (models[imodel],variables[ivariable]['name'],seasons[iseason]) )
                continue
              lines = [ l for l in content if l.startswith(variables[ivariable]['id']) ]
              if len(lines) > 1:
                raise "Found unexpected multiple entries"
              elif len(lines) == 1:
                rmse = lines[0].split(',')[-2]
                if rmse.upper() == 'NAN':
                  print("NAN: %s, %s, %s" % (models[imodel],variables[ivariable]['name'],seasons[iseason]) )
                else:
                  data[imodel,ivariable,iseason] = float(rmse)
                  #print(float(rmse),models[imodel])
              else:
                 print("Missing: %s, %s, %s" % (models[imodel],variables[ivariable]['name'],seasons[iseason]) )

  # Dictionary to hold data
  d = {}
  d['data'] = data.copy()
  d['models'] = models.copy()
  d['variables'] = variables.copy()
  d['seasons'] = seasons.copy()

  return d

# --- Function to save data into csv file ---
def write_csv(path, metrics):

  import csv

  with open(path, mode='w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    # Header: comes from the first set of metrics
    header1 = ["Model", ]
    header2 = ["", ]
    header3 = ["", ]
    for v in metrics[0]['variables']:
      for s in metrics[0]['seasons']:
          header1.append( f"{v['name']}" )
          header2.append( f"{v['units']}" )
          header3.append( f"{s}" )
    writer.writerow(header1)
    writer.writerow(header2)
    writer.writerow(header3)

    # Data: loop over all metrics
    for metric in metrics:
      for imodel in range(len(metric['models'])):
        line = [metric['models'][imodel], ]
        for ivariable in range(len(metric['variables'])):
          for iseason in range(len(metric['seasons'])):
            line.append( f"{metric['data'][imodel,ivariable,iseason]}" )
        writer.writerow(line)
  
  return

# --- Main ---

# Variables
variables = \
[

  {'name':'Net TOA',
   'units':'W m$^{-2}$',
   'id':'RESTOM global ceres_ebaf_toa_v4.1',
   'exclude':()},

  {'name':'SW CRE',
   'units':'W m$^{-2}$',
   'id':'SWCF global ceres_ebaf_toa_v4.1',
   'exclude':()},

  {'name':'LW CRE',
   'units':'W m$^{-2}$',
   'id':'LWCF global ceres_ebaf_toa_v4.1',
   'exclude':()},

  {'name':'prec',
   'units':'mm day$^{-1}$',
   'id':'PRECT global GPCP_v2.3',
   'exclude':('CIESM',)},

  {'name':'tas land',
   'units':'K',
   'id':'TREFHT land ERA5',
   'exclude':()},

  {'name':'SLP',
   'units':'hPa',
   'id':'PSL global ERA5',
   'exclude':()},

  {'name':'u-200',
   'units':'m s$^{-1}$',
   'id':'U-200mb global ERA5',
   'exclude':()},

  {'name':'u-850',
   'units':'m s$^{-1}$',
   'id':'U-850mb global ERA5',
   'exclude':()},

  {'name':'Zg-500',
   'units':'hm',
   'id':'Z3-500mb global ERA5',
   'exclude':('KIOST-ESM',)},

]

# Seasons
seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']

# Read CMIP6 data
#path = '/global/cfs/cdirs/e3sm/www/e3sm_diags_for_cmip/*/historical/r1i1p1f1/viewer/table-data'
#path = '/global/cfs/cdirs/e3sm/www/chengzhu/e3sm_diags_for_cmip/CMIP6_20211216_v2_paper_linear/*/historical/r1i1p1f1/viewer/table-data'
#path = '/global/cfs/cdirs/e3sm/www/chengzhu/e3sm_diags_for_cmip/CMIP6_20220324_v2_paper_linear/*/historical/r1i1p1f1/viewer/table-data'
path = '/global/cfs/cdirs/e3sm/www/CMIP6_comparison_1985-2014_E3SMv2_golaz_etal_2022/*/historical/r1i1p1f1/viewer/table-data'
cmip6 = read_e3sm_diags_metrics(path, variables, seasons)

# Read E3SMv2 (coupled)
path = '/global/cfs/cdirs/e3sm/www/golaz/E3SMv2/v2.LR.historical_0101/e3sm_diags/180x360_aave_cmip6/model_vs_obs_1985-2014/viewer/table-data'
E3SMv2 = read_e3sm_diags_metrics(path, variables, seasons, names=['E3SMv2',])

# Read E3SMv2 (composite base)
path = '/global/cfs/cdirs/e3sm/www/golaz/E3SMv2/v2.LR.hist-composite-all/e3sm_diags/180x360_aave_cmip6/model_vs_obs_1985-2014/viewer/table-data'
E3SMv2_composite_base = read_e3sm_diags_metrics(path, variables, seasons, names=['E3SMv2 (composite base)',])

# Read E3SMv2 (composite best)
path = '/global/cfs/cdirs/e3sm/www/golaz/E3SMv2/v2.LR.hist-composite-073_017/e3sm_diags/180x360_aave_cmip6/model_vs_obs_1985-2014/viewer/table-data'
E3SMv2_composite_best = read_e3sm_diags_metrics(path, variables, seasons, names=['E3SMv2 (composite best)',])

# Save to csv
write_csv('cmip6.csv', [cmip6, E3SMv2, E3SMv2_composite_base, E3SMv2_composite_best])

# -----------------------------------------------------------------------------
# Create plot: first only with CMIP6, E3SMv1 and v2
fig = plt.figure(figsize=[12,9])
nsx = 4
nsy = 3

nmodels = len(cmip6['models'])
nvariables = len(variables)
nseasons = len(seasons)
for ivariable in range(nvariables):

  # CMIP6 data for box and whiskers
  data = []
  labels = []
  for iseason in range(nseasons):
     # Identify model with lowest RMSE
     ibest = ma.argmin( cmip6['data'][:,ivariable,iseason].compressed() )
     print("Best model %s %s %s" % (variables[ivariable]['name'],seasons[iseason],cmip6['models'][ibest]))
     # Remove missing data using 'compressed()' function
     data.append( cmip6['data'][:,ivariable,iseason].compressed() )
     labels.append(seasons[iseason])
  cmip6_stats = cbook.boxplot_stats(data,whis=[0,100],labels=labels)

  # Plot panel
  ax = plt.subplot(nsy, nsx, ivariable+int(ivariable/3)+1)
  ax.set_box_aspect(1)

  # CMIP6 ensemble
  ax.bxp(cmip6_stats)

  # E3SMv1
  x = np.arange(nseasons)+0.8
  iE3SMv1 = cmip6['models'].index('E3SM-1-0')
  ax.scatter(x,cmip6['data'][iE3SMv1,ivariable,:],color='b',marker='>',label="E3SMv1 (0101)")

  # E3SMv2 (coupled)
  x = np.arange(nseasons)+1.2
  ax.scatter(x,E3SMv2['data'][0,ivariable,:],color='r',marker='<',label="E3SMv2 (0101)")

  ## E3SMv2 (amip)
  #x = np.arange(nseasons)+1.4
  #ax.scatter(x,E3SMv2_amip['data'][0,ivariable,:],color='r',marker="^")

  ## E3SMv2 (composite base)
  #x = np.arange(nseasons)+0.75
  #ax.scatter(x,E3SMv2_composite_base['data'][0,ivariable,:],color='r',marker="*",label="Composite base")

  ## E3SMv2 (composite best)
  #x = np.arange(nseasons)+0.75
  #ax.scatter(x,E3SMv2_composite_best['data'][0,ivariable,:],color='GoldenRod',marker="*",label='Composite best')

  # Customize plot
  ax.set_title('('+chr(97+ivariable)+')', loc="left")
  ax.set_title(variables[ivariable]['name']+' ('+variables[ivariable]['units']+')', loc="right")
  ax.set_xlim([0.4,nseasons+0.9])

fig.subplots_adjust(wspace=0.3,hspace=0.3)

# Legend base on last subplot
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc=(0.76,0.8))

fig.savefig("cmip6.pdf",bbox_inches='tight')

# -----------------------------------------------------------------------------
# Create plot: second CMIP6, E3SMv2 and composite configurations
fig = plt.figure(figsize=[12,9])
nsx = 4
nsy = 3

nmodels = len(cmip6['models'])
nvariables = len(variables)
nseasons = len(seasons)
for ivariable in range(nvariables):

  # CMIP6 data for box and whiskers
  data = []
  labels = []
  for iseason in range(nseasons):
     # Identify model with lowest RMSE
     ibest = ma.argmin( cmip6['data'][:,ivariable,iseason].compressed() )
     print("Best model %s %s %s" % (variables[ivariable]['name'],seasons[iseason],cmip6['models'][ibest]))
     # Remove missing data using 'compressed()' function
     data.append( cmip6['data'][:,ivariable,iseason].compressed() )
     labels.append(seasons[iseason])
  cmip6_stats = cbook.boxplot_stats(data,whis=[0,100],labels=labels)

  # Plot panel
  ax = plt.subplot(nsy, nsx, ivariable+int(ivariable/3)+1)
  ax.set_box_aspect(1)

  # CMIP6 ensemble
  ax.bxp(cmip6_stats)

  # E3SMv2 (coupled)
  x = np.arange(nseasons)+1.2
  ax.scatter(x,E3SMv2['data'][0,ivariable,:],color='r',marker='<',label="E3SMv2 (0101)")

  # E3SMv2 (composite base)
  x = np.arange(nseasons)+0.8
  ax.scatter(x,E3SMv2_composite_base['data'][0,ivariable,:],color='r',marker="*",label="Composite base")

  # E3SMv2 (composite best)
  x = np.arange(nseasons)+0.8
  ax.scatter(x,E3SMv2_composite_best['data'][0,ivariable,:],color='GoldenRod',marker="*",label='Composite best')

  # Customize plot
  ax.set_title('('+chr(97+ivariable)+')', loc="left")
  ax.set_title(variables[ivariable]['name']+' ('+variables[ivariable]['units']+')', loc="right")
  ax.set_xlim([0.4,nseasons+0.9])

fig.subplots_adjust(wspace=0.3,hspace=0.3)

# Legend base on last subplot
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc=(0.76,0.8))

fig.savefig("cmip6_composite.pdf",bbox_inches='tight')

