#!/usr/bin/env python3
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.transforms as mtransforms

# This script is to compare rmse from e3sm_diags tables and those from cmip6 models in a pre-compiled csv file. 
# The script was written by Chris Golaz (Golaz et al. 2022 E3SM v2 overview paper) and adapted by Jill Zhang for implementation in e3sm_diags.  

# --- Function to read CMIP6 model metrics  ---
def read_cmip6_metrics_from_csv(path, variables, seasons): 

    models=[]
    
    with open (path, 'r') as fin:
        # skip 3 header lines and last 2 E3SMv2 composites 
        rmse = fin.readlines()[3:-2]
        nmodels = len(rmse)
        nvariables = len(variables)
        nseasons = len(seasons)
        print(nmodels, nvariables,nseasons)
        data = ma.array(np.zeros((nmodels,nvariables,nseasons)), mask=True)
        for imodel,line in enumerate(rmse):
            models.append(line.split(',')[0])
            model_line = line
            print(model_line)
            for ivariable in range(nvariables):
                rmse_seasons = model_line.split(',')[1+ivariable*5:6+ivariable*5]
                print(rmse_seasons)
                if any(x=="--" for x in rmse_seasons):
                    print(model_line.split(',')[1+ivariable*5:6+ivariable*5])
                    pass
                    print(type(model_line.split(',')[1+ivariable*5:6+ivariable*5][0]))
                else:
                    print('second', model_line.split(',')[1+ivariable*5:6+ivariable*5])
                    data[imodel,ivariable,:] = [float(x) for x in model_line.split(',')[1+ivariable*5:6+ivariable*5]]
                print(data[imodel,ivariable,:])

    # Dictionary to hold data
    d = {}
    d['data'] = data.copy()
    d['models'] = models.copy()
    d['variables'] = variables.copy()
    d['seasons'] = seasons.copy()
    
    return d

# --- Function to read E3SM Diags metrics  ---
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
nseasons = len(seasons)


cmip6 = read_cmip6_metrics_from_csv('cmip6_seasonal_rmse.csv', variables, seasons)

# Create plot: first only with CMIP6, E3SMv1 and v2
fig = plt.figure(figsize=[12,9])
nsx = 4
nsy = 3
for ivariable in range(len(variables)):
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
  iE3SMv2 = cmip6['models'].index('E3SMv2')
  ax.scatter(x,cmip6['data'][iE3SMv2,ivariable,:],color='r',marker='<',label="E3SMv2 (0101)")

  # Customize plot
  ax.set_title('('+chr(97+ivariable)+')', loc="left")
  ax.set_title(variables[ivariable]['name']+' ('+variables[ivariable]['units']+')', loc="right")
  ax.set_xlim([0.4,nseasons+0.9])

fig.subplots_adjust(wspace=0.3,hspace=0.3)

# Legend base on last subplot
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc=(0.76,0.8))

fig.savefig("cmip6.pdf",bbox_inches='tight')

    
    
    
    
    
