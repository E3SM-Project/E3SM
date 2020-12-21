import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

# quad
print("Quad")
os.system("rm grid.nc")
os.system("ln -s square_mesh_80x80_culled.nc grid.nc")
os.system("rm ic.nc")
os.system("ln -s ic_quad.nc ic.nc")

# Wachspress
print("  Wachspress")
nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"variational",
                                "config_variational_basis":"wachspress"}}

f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.wachspress")

os.system("rm -rf namelist.seaice streams.seaice output_quad_wachspress")
os.system("ln -s namelist.seaice.wachspress namelist.seaice")
os.system("ln -s streams.seaice.square streams.seaice")

os.system("../../../../seaice_model")

os.system("mv output output_quad_wachspress")

# PWL
print("  PWL")
nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"variational",
                                "config_variational_basis":"pwl"}}

f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.pwl")

os.system("rm -rf namelist.seaice streams.seaice output_quad_pwl")
os.system("ln -s namelist.seaice.pwl namelist.seaice")
os.system("ln -s streams.seaice.square streams.seaice")

os.system("../../../../seaice_model")

os.system("mv output output_quad_pwl")

# Weak
print("  Weak")
nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"weak"}}

f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.weak")

os.system("rm -rf namelist.seaice streams.seaice output_quad_weak")
os.system("ln -s namelist.seaice.weak namelist.seaice")
os.system("ln -s streams.seaice.square streams.seaice")

os.system("../../../../seaice_model")

os.system("mv output output_quad_weak")

# hex
print("Hex")
os.system("rm grid.nc")
os.system("ln -s square_mesh_82x94_culled.nc grid.nc")
os.system("rm ic.nc")
os.system("ln -s ic_hex.nc ic.nc")

# Wachspress
print("  Wachspress")
nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"variational",
                                "config_variational_basis":"wachspress"}}

f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.wachspress")

os.system("rm -rf namelist.seaice streams.seaice output_hex_wachspress")
os.system("ln -s namelist.seaice.wachspress namelist.seaice")
os.system("ln -s streams.seaice.square streams.seaice")

os.system("../../../../seaice_model")

os.system("mv output output_hex_wachspress")

# PWL
print("  PWL")
nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"variational",
                                "config_variational_basis":"pwl"}}

f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.pwl")

os.system("rm -rf namelist.seaice streams.seaice output_hex_pwl")
os.system("ln -s namelist.seaice.pwl namelist.seaice")
os.system("ln -s streams.seaice.square streams.seaice")

os.system("../../../../seaice_model")

os.system("mv output output_hex_pwl")

# Weak
print("  Weak")
nmlPatch = {"velocity_solver": {"config_stress_divergence_scheme":"weak"}}

f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.weak")

os.system("rm -rf namelist.seaice streams.seaice output_hex_weak")
os.system("ln -s namelist.seaice.weak namelist.seaice")
os.system("ln -s streams.seaice.square streams.seaice")

os.system("../../../../seaice_model")

os.system("mv output output_hex_weak")
