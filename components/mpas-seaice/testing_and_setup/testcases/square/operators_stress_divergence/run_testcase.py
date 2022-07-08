from create_grids import create_grids
from create_ics import create_ics
from run_model import run_model
from stress_divergence_scaling import stress_divergence_scaling
from error_analysis_stress_divergence import error_analysis_stress_divergence

create_grids()

create_ics()

run_model()

stress_divergence_scaling()

error_analysis_stress_divergence()
