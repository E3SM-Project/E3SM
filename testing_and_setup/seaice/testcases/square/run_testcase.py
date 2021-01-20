import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

operatorMethods = ["wachspress","pwl","weak"]

gridTypes = ["hex","quad"]
#gridTypes = ["quad"]

grids = {"hex" :["0082x0094",
                 "0164x0188",
                 "0328x0376",
                 "0656x0752"],
         "quad":["0080x0080",
                 "0160x0160",
                 "0320x0320",
                 "0640x0640"]}
#grids = {"hex" :["0082x0094"],
#         "quad":["0080x0080"]}


#subcycleNumbers = [120,240,480,960,1920,3840,7680,15360,30720]
subcycleNumbers = [120,240,480,960,1920,3840,7680]
#subcycleNumbers = [120,240]

for gridType in gridTypes:

    print("Grid type: ", gridType)

    for operatorMethod in operatorMethods:

        print("  Operator Method: ", operatorMethod)

        for grid in grids[gridType]:

            print("    Grid: ", grid)

            os.system("rm grid.nc")
            os.system("rm ic.nc")
            os.system("ln -s grid_%s_%s.nc grid.nc" %(gridType,grid))
            os.system("ln -s ic_%s_%s.nc ic.nc" %(gridType,grid))

            for subcycleNumber in subcycleNumbers:

                print("      Subcycle number: ", subcycleNumber)

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

                os.system("rm -rf namelist.seaice streams.seaice output_%s_%s_%s_%i" %(gridType, operatorMethod, grid, subcycleNumber))
                os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(operatorMethod, subcycleNumber))
                os.system("ln -s streams.seaice.square streams.seaice")

                os.system("../../../../seaice_model")

                os.system("mv output output_%s_%s_%s_%i" %(gridType, operatorMethod, grid, subcycleNumber))
