import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model():

    operatorMethods = ["wachspress","pwl","weak"]
    subcycleNumbers = [120,240,480,960,1920,3840,7680]

    for operatorMethod in operatorMethods:

        print("Operator Method: ", operatorMethod)

        for subcycleNumber in subcycleNumbers:

            print("  Subcycle number: ", subcycleNumber)

            if (operatorMethod == "wachspress"):
                nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"wachspress",
                                                "config_elastic_subcycle_number":subcycleNumber}}
            elif (operatorMethod == "pwl"):
                nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"pwl",
                                                "config_elastic_subcycle_number":subcycleNumber}}
            elif (operatorMethod == "weak"):
                nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"weak",
                                                "config_elastic_subcycle_number":subcycleNumber}}

            f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.%s.%i" %(operatorMethod, subcycleNumber))

            os.system("rm -rf namelist.seaice streams.seaice output_%s_%i" %(operatorMethod, subcycleNumber))
            os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(operatorMethod, subcycleNumber))
            os.system("ln -s streams.seaice.square streams.seaice")

            os.system("../../../../../seaice_model")

            os.system("mv output output_hex_%s_%i" %(operatorMethod, subcycleNumber))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
