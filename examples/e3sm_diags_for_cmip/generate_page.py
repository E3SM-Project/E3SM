# Script to generate html summary page of E3SM Diags applied to CMIP6 models

import collections
import jinja2
import glob
import os.path
import pandas as pd
import pprint
import shlex
import subprocess
import utils

def tree(): return collections.defaultdict(tree)

def table_elements(search):

    # Performance data to retrieve
    variables = [ \
       "PRECT global GPCP_v2.3",
       "RESTOM global ceres_ebaf_toa_v4.1",
       "FSNTOA global ceres_ebaf_toa_v4.1", "FLUT global ceres_ebaf_toa_v4.1",
       "SWCF global ceres_ebaf_toa_v4.1", "LWCF global ceres_ebaf_toa_v4.1",
       "U-850mb global ERA5", "U-200mb global ERA5",
       "T-850mb global ERA5", "T-200mb global ERA5",
       "Z3-500mb global ERA5"
    ]
    metrics = ["Unit", "RMSE", "Mean_Bias", "Correlation"]
    seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

    # Translation lexicon
    lex = { \
    "PRECT global GPCP_v2.3":"Precipitation (vs GPCP v2.3)",
    "RESTOM global ceres_ebaf_toa_v4.1":"TOA net (vs CERES-EBAAF Ed4.1)",
    "FSNTOA global ceres_ebaf_toa_v4.1":"TOA SW (vs CERES-EBAAF Ed4.1)",
    "FLUT global ceres_ebaf_toa_v4.1":"TOA LW (vs CERES-EBAAF Ed4.1)",
    "SWCF global ceres_ebaf_toa_v4.1":"TOA SWCRE (vs CERES-EBAAF Ed4.1)",
    "LWCF global ceres_ebaf_toa_v4.1":"TOA LWCRE (vs CERES-EBAAF Ed4.1)",
    "U-850mb global ERA5":"u 850 hPa (vs ERA5)",
    "U-200mb global ERA5":"u 200 hPa (vs ERA5)",
    "T-850mb global ERA5":"t 850 hPa (vs ERA5)",
    "T-200mb global ERA5":"t 200 hPa (vs ERA5)",
    "Z3-500mb global ERA5":"Geo-Z 500 hPA (vs ERA5)",
    "Unit":"unit",
    "RMSE":"rmse",
    "Mean_Bias":"bias",
    "Correlation":"corr",
    }

    # Mapping between variables and viewer relative path
    relpath = { \
    "PRECT global GPCP_v2.3":"lat_lon/gpcp_v23/prect-global-gpcp_v23",
    "RESTOM global ceres_ebaf_toa_v4.1":"lat_lon/ceres-ebaf-toa-v41/restom-global-ceres_ebaf_toa_v41",
    "FSNTOA global ceres_ebaf_toa_v4.1":"lat_lon/ceres-ebaf-toa-v41/fsntoa-global-ceres_ebaf_toa_v41",
    "FLUT global ceres_ebaf_toa_v4.1":"lat_lon/ceres-ebaf-toa-v41/flut-global-ceres_ebaf_toa_v41",
    "SWCF global ceres_ebaf_toa_v4.1":"lat_lon/ceres-ebaf-toa-v41/swcf-global-ceres_ebaf_toa_v41",
    "LWCF global ceres_ebaf_toa_v4.1":"lat_lon/ceres-ebaf-toa-v41/lwcf-global-ceres_ebaf_toa_v41",
    "U-850mb global ERA5":"lat_lon/era5/u-850mb-global-era5",
    "U-200mb global ERA5":"lat_lon/era5/u-200mb-global-era5",
    "T-850mb global ERA5":"lat_lon/era5/t-850mb-global-era5",
    "T-200mb global ERA5":"lat_lon/era5/t-200mb-global-era5",
    "Z3-500mb global ERA5":"lat_lon/era5/z3-500mb-global-era5",
    }

    # Loop over all simulations to gather data
    table = []
    simulations = glob.glob(search)
    n_simulations = len(simulations)
    for simulation in simulations:
        print('\n=== %s ===' % (simulation))

        # Split source path, remove empty strings
        p = list(filter(None, simulation.split('/')))

        # Extract relevant data
        c = {}
        c['xmls'] = simulation
        c['realization'] = p[-1]
        c['experiment'] = p[-2]
        c['model'] = p[-3]
        c['institution'] = p[-4]
        #c['www'] = "/var/www/acme/acme-diags/zhang40/CMIP6/%s/%s/%s" \
        #c['www'] = "/var/www/acme/acme-diags/e3sm_diags_for_cmip/%s/%s/%s" \
        c['www'] = "/var/www/acme/acme-diags/zhang40/CMIP6_20240109_1985-2014/%s/%s/%s" \
                   % (c['model'],c['experiment'],c['realization'])

        print(c['www'])
        d = tree()
        try:
            for season in seasons:

              csv = "%s/viewer/table-data/%s_metrics_table.csv" % (c['www'],season)
              #print(csv)
              df = pd.read_csv(csv, index_col="Variables")
              for variable in variables:
                for metric in metrics:
                  if variable in df.index:
                    val = df.loc[variable,metric]
                  else:
                    val = None
                  #print(variable,metric,season,val)
                  d[lex[variable]][lex[metric]][season] = val

            c['data'] = d
            table.append(c)
        except OSError as err:
            print("OS error: {0}".format(err))
            n_simulations = n_simulations - 1
            print(n_simulations, ' number of simulations')
            pass
        
            

    """
    # Print for verification
    pp = pprint.PrettyPrinter(indent=2)
    for s in table:
        print(s['xmls'])
        pp.pprint(s['data'])
    """
    
    # Output portion to html page.
    # Data should be separated into header and content
    # - header should be a list
    # - content should be a 2d list (nested list)
    # Once that done, should be straightforward to fill the table
    # within the Jinja2 template.


    # List of fields
    n = len(seasons)
    fields = {}
    fields['name'] = []
    fields['visible'] = []
    for i,variable in enumerate(variables):
      print(i,variable)
      fields['name'].append(lex[variable])
      visible = [0,1] + list(range(n*i+2,n*(i+1)+2))
      fields['visible'].append([0,1] + list(range(n*i+2,n*(i+1)+2)))
    # Initial field visibility
    fields['initial'] = list(range(n*1+2,n*(len(variables))+2))

    # Table header
    header = ["Model", "Realization"]
    for variable in variables:
      header.extend(["ANN", "DJF", "MAM", "JJA", "SON"])

    # Table content 
    ncols = len(header)
    #nrows = len(simulations)
    print(n_simulations, ' number of simulations')
    nrows = n_simulations
    content = [[None] * ncols for i in range(nrows)]
    for i in range(nrows):
        content[i][0] = table[i]['model']
        content[i][1] = \
         '<a href="%s/%s/%s/viewer/">%s</a>' \
          % (table[i]['model'],table[i]['experiment'],table[i]['realization'],table[i]['realization'])
        j = 2
        for variable in variables:
          for metric in ['rmse',]:
              for season in seasons:
                val = table[i]['data'][lex[variable]][metric][season]
                if val == None:
                  val = 999.999
                  content[i][j] = "%.3f" % val
                else:
                  content[i][j] = \
                   '<a href="%s/%s/%s/viewer/%s/%s.html">%.3f</a>' \
                    % ( table[i]['model'],table[i]['experiment'],table[i]['realization'], \
                        relpath[variable],season.lower(),val )
                j += 1


    return fields, header, content


# Gether table elements
c = {'fields':[], 'header':[], 'content':[] }

# amip simulations
#fields, header, content = table_elements('/var/www/acme/acme-diags/zhang40/CMIP6/*/amip/r1i1p1f1/')
#fields, header, content = table_elements('/var/www/acme/acme-diags/e3sm_diags_for_cmip/*/amip/*/')
fields, header, content = table_elements('/var/www/acme/acme-diags/zhang40/CMIP6_20240109_1985-2014/*/amip/*/')
c['fields'].append(fields)
c['header'].append(header)
c['content'].append(content)

# historical simulations
#fields, header, content = table_elements('/var/www/acme/acme-diags/zhang40/CMIP6/*/historical/r1i1p1f1/')
#fields, header, content = table_elements('/var/www/acme/acme-diags/e3sm_diags_for_cmip/*/historical/*/')
fields, header, content = table_elements('/var/www/acme/acme-diags/zhang40/CMIP6_20240109_1985-2014/*/historical/*/')
c['fields'].append(fields)
c['header'].append(header)
c['content'].append(content)

# Initialize jinja2 template engine
templateLoader = jinja2.FileSystemLoader(searchpath="./templates")
templateEnv = jinja2.Environment( loader=templateLoader )
template = templateEnv.get_template( 'cmip6_template.html' )

# Instantiate page
path = os.path.join('/var/www/acme/acme-diags/zhang40/CMIP6_20240109_1985-2014','index.html')
with open(path, 'w') as f:
    f.write(template.render( **c ))

