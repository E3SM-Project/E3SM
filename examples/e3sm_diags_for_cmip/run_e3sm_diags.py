# Script to generate e3sm_diags batch script for CMIP6 model output

import jinja2
import glob
import os.path
import shlex
import subprocess
import utils

# Location of xml files
input = '/home/zhang40/e3sm_diags_for_CMIP6/CMIP6_20240109_wE3SM/CMIP/*/*/amip/r1i1p1f1/'
#input = '/home/zhang40/e3sm_diags_for_CMIP6/CMIP6_20240109/CMIP/*/*/historical/r1i1p1f1/'
print(input)

# Initialize jinja2 template engine
templateLoader = jinja2.FileSystemLoader(searchpath='.')
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template('e3sm_diags.bash')

# Create script subdirectory
if not os.path.exists('scripts'):
    os.makedirs('scripts')
os.chdir('scripts')

# Loop over all simulations
simulations = glob.glob(input)
for simulation in simulations:
    print('\n=== %s ===' % (simulation))

    # Split source path, remove empty strings
    p = list(filter(None, simulation.split('/')))

    # Extract relevant data
    c = {}
    c['simulation'] = simulation
    c['realization'] = p[-1]
    c['experiment'] = p[-2]
    c['model'] = p[-3]
    c['institution'] = p[-4]
   
    # Create script
    #scriptFile = 'e3sm_diags_%s.bash' % (c['model'])
    scriptFile = 'e3sm_diags_{}_{}_{}.bash'.format(c['model'], c['experiment'], c['realization'])
    with open(scriptFile, 'w') as f:
        f.write(template.render( **c ))

    # Submit script
    utils.submitScript(scriptFile)

