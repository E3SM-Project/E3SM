from create_grids import create_grids
from create_ics import create_ics
from run_model import run_model
from set_difference_fields import set_difference_fields
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_square_quadhex_testcase():

    create_grids()

    create_ics()

    run_model()

    set_difference_fields()

    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_square_quadhex_testcase()
