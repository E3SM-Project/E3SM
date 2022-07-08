import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    advectionMethods = ["IR","upwind"]
    #advectionMethods = ["IR"]

    icTypes = ["cosine_bell","slotted_cylinder"]
    #icTypes = ["cosine_bell"]

    gridSizes = [2562, 10242, 40962, 163842]
    #gridSizes = [2562]

    for advectionMethod in advectionMethods:

        print("Advection method: ", advectionMethod)

        for icType in icTypes:

            print("  IC type: ", icType)

            for gridSize in gridSizes:

                print("    Gridsize: ", gridSize)

                os.system("rm grid.nc ic.nc namelist.seaice")
                os.system("ln -s x1.%i.grid.nc grid.nc" %(gridSize))
                os.system("ln -s ic_%s_%i.nc ic.nc" %(icType, gridSize))

                if (advectionMethod == "IR"):
                    nmlPatch = {"advection": {"config_advection_type":"incremental_remap"}}
                elif (advectionMethod == "upwind"):
                    nmlPatch = {"advection": {"config_advection_type":"upwind"}}

                f90nml.patch("namelist.seaice.advection", nmlPatch, "namelist.seaice.%s" %(advectionMethod))

                os.system("rm -rf output_%s_%s_%i" %(advectionMethod, icType, gridSize))
                os.system("ln -s namelist.seaice.%s namelist.seaice" %(advectionMethod))

                os.system("%s ../../../../seaice_model" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND))

                os.system("mv output output_%s_%s_%i" %(advectionMethod, icType, gridSize))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
