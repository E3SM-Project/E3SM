from create_grids import create_grids
from create_ics import create_ics
from run_model import run_model

#-------------------------------------------------------------------------------

def run_square_quadhex_testcase():

    create_grids()

    create_ics()

    run_model()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_square_quadhex_testcase()
