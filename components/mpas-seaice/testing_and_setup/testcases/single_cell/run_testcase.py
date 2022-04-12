import os
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_testcase():

    # copy namelist and streams file
    os.system("cp ../../configurations/standard_physics_single_cell/namelist.seaice .")
    os.system("cp ../../configurations/standard_physics_single_cell/streams.seaice .")

    # forcing
    os.system("python ../../testing/DATA/domain_sc_71.35_-156.5/get_domain.py")
    
    # run MPAS-Seaice
    os.system("../../../seaice_model")

    # plot output
    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
