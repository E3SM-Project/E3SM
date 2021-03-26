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

    operatorMethods = ["wachspress","pwl","weak"]

    gridTypes = ["hex"]

    for gridType in gridTypes:

        print("Grid type: ", gridType)

        for operatorMethod in operatorMethods:

            print("  Operator Method: ", operatorMethod)

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

            f90nml.patch("namelist.seaice.strain", nmlPatch, "namelist.seaice.%s" %(operatorMethod))

            os.system("rm -rf namelist.seaice streams.seaice output_%s_%s" %(gridType, operatorMethod))
            os.system("ln -s namelist.seaice.%s namelist.seaice" %(operatorMethod))
            os.system("ln -s streams.seaice.strain streams.seaice")

            os.system("%s ../../../../seaice_model" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND))

            os.system("mv output output_%s_%s" %(gridType, operatorMethod))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
