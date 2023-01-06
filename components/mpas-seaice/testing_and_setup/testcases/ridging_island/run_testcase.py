from create_grids import create_grids
from run_model import run_model
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_ridging_island_testcase():

    # create grids
    create_grids()

    # run the model
    run_model()

    # plot test case
    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_ridging_island_testcase()
