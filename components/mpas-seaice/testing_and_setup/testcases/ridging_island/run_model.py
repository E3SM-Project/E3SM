import os

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    islandTypes = ["no_island","island"]

    for islandType in islandTypes:

        print("Island type: ", islandType)

        os.system("rm grid.nc")
        os.system("rm ic.nc")
        os.system("rm special_boundaries.nc")
        os.system("ln -s grid_%s.nc grid.nc" %(islandType))
        os.system("ln -s ic_%s.nc ic.nc" %(islandType))
        os.system("ln -s special_boundaries_%s.nc special_boundaries.nc" %(islandType))

        os.system("rm -rf namelist.seaice streams.seaice output_%s" %(islandType))
        os.system("ln -s namelist.seaice.ridging_island namelist.seaice")
        os.system("ln -s streams.seaice.ridging_island streams.seaice")

        os.system("%s ../../../../seaice_model" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND))

        os.system("mv output output_%s" %(islandType))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
