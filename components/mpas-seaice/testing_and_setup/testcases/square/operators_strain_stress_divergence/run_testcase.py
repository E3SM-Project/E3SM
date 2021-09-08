from create_grids import create_grids
from create_ics import create_ics
from run_model import run_model
from strain_stress_divergence_map import strain_stress_divergence_map
from strain_stress_divergence_scaling import strain_stress_divergence_scaling

create_grids()

create_ics()

run_model()

strain_stress_divergence_map()

strain_stress_divergence_scaling()
