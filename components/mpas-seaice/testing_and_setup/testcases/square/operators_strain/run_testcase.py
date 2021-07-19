from create_grids import create_grids
from create_ics import create_ics
from run_model import run_model
from average_variational_strains import average_variational_strains
from strain_map import strain_map
from strain_scaling import strain_scaling
from error_analysis_strain import error_analysis_strain

create_grids()

create_ics()

run_model()

average_variational_strains()

strain_map()

strain_scaling()

error_analysis_strain()
