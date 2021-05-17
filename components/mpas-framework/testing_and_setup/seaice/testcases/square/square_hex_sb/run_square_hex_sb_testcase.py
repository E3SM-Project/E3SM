from create_grids import create_grids
from create_boundary_info import create_boundary_info
from create_ic import create_ic
from run_model import run_model
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_square_hex_sb_testcase():

    create_grids()

    create_boundary_info()

    create_ic()

    run_model()

    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_square_hex_sb_testcase()
