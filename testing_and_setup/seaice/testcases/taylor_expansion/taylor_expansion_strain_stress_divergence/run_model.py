import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise
import sys

#-------------------------------------------------------------------------------

def run_model():

    gridTypes = ["hex"]
    operatorMethods = ["wachspress","pwl","weakwachs","weakpwl"]
    velocityTypes = \
        ["xlinear",   "ylinear",   "xylinear",\
         "xquadratic","yquadratic","xyquadratic",\
         "xcubic",    "ycubic",    "xycubic",\
         "xquartic",  "yquartic",  "xyquartic"]
    operatorSupermethod = "var"

    for gridType in gridTypes:

        print("Grid type: ", gridType)

        for operatorMethod in operatorMethods:

            print("  Operator Method: ", operatorMethod)

            for velocityType in velocityTypes:

                print("    velocityType: ", velocityType)

                os.system("rm grid.nc")
                os.system("rm ic.nc")
                os.system("ln -s grid_%s_%s.nc grid.nc" %(operatorSupermethod, gridType))
                os.system("ln -s ic_%s_%s_%s.nc ic.nc" %(operatorSupermethod, gridType, velocityType))

                if (operatorMethod == "wachspress"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                    "config_stress_divergence_scheme":"variational",
                                                    "config_variational_basis":"wachspress",
                                                    "config_variational_denominator_type":"original"}}
                elif (operatorMethod == "pwl"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                    "config_stress_divergence_scheme":"variational",
                                                    "config_variational_basis":"pwl",
                                                    "config_variational_denominator_type":"original"}}
                elif (operatorMethod == "weak"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                    "config_stress_divergence_scheme":"weak",
                                                    "config_variational_denominator_type":"original"}}
                elif (operatorMethod == "weakwachs"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                    "config_stress_divergence_scheme":"variational",
                                                    "config_variational_basis":"wachspress",
                                                    "config_variational_denominator_type":"original"}}
                elif (operatorMethod == "weakpwl"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                    "config_stress_divergence_scheme":"variational",
                                                    "config_variational_basis":"pwl",
                                                    "config_variational_denominator_type":"original"}}

                f90nml.patch("namelist.seaice.strain_stress_divergence", nmlPatch, "namelist.seaice.%s" %(operatorMethod))

                os.system("rm -rf namelist.seaice streams.seaice output_%s_%s_%s" %(operatorMethod, gridType, velocityType))
                os.system("ln -s namelist.seaice.%s namelist.seaice" %(operatorMethod))
                os.system("ln -s streams.seaice.strain_stress_divergence streams.seaice")

                os.system("../../../../../seaice_model")

                os.system("mv output output_%s_%s_%s" %(operatorMethod, gridType, velocityType))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
