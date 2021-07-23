from get_testcase_data import get_testcase_data
from create_ic import create_ic
from run_model import run_model
from average_variational_strains import average_variational_strains
from strain_map import strain_map
from strain_scaling import strain_scaling

#-------------------------------------------------------------------------------

def run_strain_testcase():

    get_testcase_data()

    create_ic()

    run_model()

    average_variational_strains()

    strain_map()

    strain_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_strain_testcase()
