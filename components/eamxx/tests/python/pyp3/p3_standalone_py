#!/usr/bin/env python3

import sys

# Add path to scream libs
sys.path.append('@SCREAM_BASE_DIR@/scripts')

# Add path to pyeamxx libs
sys.path.append('@CMAKE_BINARY_DIR@/src/python')

# Without these, and manual init/finalize, on my laptop I get
# obscure MPI errors at exit.
import mpi4py
mpi4py.rc.initialize = False  # do not initialize MPI automatically
mpi4py.rc.finalize = False    # do not finalize MPI automatically

from mpi4py import MPI
import pyeamxx
from pathlib import Path

from utils import ensure_yaml
ensure_yaml()
import yaml

#########################################
def main ():
#########################################

    # Get timestepping params
    with open('input.yaml','r') as fd:
        yaml_input = yaml.load(fd,Loader=yaml.SafeLoader)
    nsteps = yaml_input['time_stepping']['number_of_steps']
    dt     = yaml_input['time_stepping']['time_step']
    t0_str = yaml_input['time_stepping']['run_t0']

    ic_file = Path('@SCREAM_DATA_DIR@/init/screami_unit_tests_ne2np4L72_20220822.nc')

    # Create the grid
    ncols = 218
    nlevs = 72
    pyeamxx.create_grids_manager(ncols,nlevs,str(ic_file))

    p3 = pyeamxx.AtmProc(yaml_input['atmosphere_processes']['p3'],'p3')
    params = p3.get_params()
    old = params.get_dbl('max_total_ni')
    print (f"max_total_ni: {params.get_dbl('max_total_ni')}")
    params.set("max_total_ni",1000000.0)
    print (f"max_total_ni: {params.get_dbl('max_total_ni')}")
    params.set("max_total_ni",old)
    print (f"max_total_ni: {params.get_dbl('max_total_ni')}")

    missing = p3.read_ic(str(ic_file))
    if len(missing)>0:
        print (f"WARNING! The following input fields were not found in the IC file, and must be manually initialized: {missing}")
    p3.initialize(t0_str)
    p3.setup_output("output_py.yaml")
    
    # Time looop
    for n in range(0,nsteps):
        p3.run(dt)

####################################
if  __name__  == "__main__":
    # This level of indirection ensures all pybind structs are destroyed
    # before we finalize eamxx (and hence kokkos)
    MPI.Init()
    pyeamxx.init()
    main ()
    pyeamxx.finalize()
    MPI.Finalize()
