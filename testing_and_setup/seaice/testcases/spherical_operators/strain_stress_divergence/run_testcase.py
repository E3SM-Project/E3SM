import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

gridSizes = [2562, 10242, 40962, 163842]

operatorMethods = ["wachspress","pwl","weak"]

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

        f90nml.patch("namelist.seaice.strain_stress_divergence", nmlPatch, "namelist.seaice.%s.%i" %(operatorMethod, gridSize))

        os.system("rm -rf namelist.seaice streams.seaice output_%s_%i" %(operatorMethod, gridSize))
        os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(operatorMethod, gridSize))
        os.system("ln -s streams.seaice.strain_stress_divergence streams.seaice")

        os.system("../../../../../seaice_model")

        os.system("mv output output_%s_%i" %(operatorMethod, gridSize))
