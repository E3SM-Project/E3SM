import subprocess

filenames = ["ic_hex.nc",
             "ic_quad.nc",
             "square_mesh_80x80_culled.nc",
             "square_mesh_82x94_culled.nc"]

dirName = "https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/testcases/square/"

for filename in filenames:

    args = ["wget", dirName+filename]

    process = subprocess.Popen(args, stdout=subprocess.PIPE)

    while process.poll() is None:
        line = process.stdout.readline()
