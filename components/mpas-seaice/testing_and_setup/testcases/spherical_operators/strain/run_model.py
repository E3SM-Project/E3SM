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

    gridSizes = [2562, 10242, 40962, 163842]
    #gridSizes = [10242]

    operatorMethods = ["wachspress","pwl","weak","weakwachs"]
    #operatorMethods = ["pwl"]

    for operatorMethod in operatorMethods:

        print("Operator Method: ", operatorMethod)

        for gridSize in gridSizes:

            print("  Gridsize: ", gridSize)

            os.system("rm grid.nc ic.nc")
            os.system("ln -s x1.%i.grid.nc grid.nc" %(gridSize))
            os.system("ln -s ic_%i.nc ic.nc" %(gridSize))

            if (operatorMethod == "wachspress"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"wachspress"}}
            elif (operatorMethod == "pwl"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"pwl"}}
            elif (operatorMethod == "weak"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                "config_stress_divergence_scheme":"weak"}}
            elif (operatorMethod == "weakwachs"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"wachspress"}}

            f90nml.patch("namelist.seaice.strain", nmlPatch, "namelist.seaice.%s.%i" %(operatorMethod, gridSize))

            os.system("rm -rf namelist.seaice streams.seaice output_%s_%i" %(operatorMethod, gridSize))
            os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(operatorMethod, gridSize))
            os.system("ln -s streams.seaice.strain streams.seaice")

            os.system("%s ../../../../../seaice_model" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND))

            os.system("mv output output_%s_%i" %(operatorMethod, gridSize))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
