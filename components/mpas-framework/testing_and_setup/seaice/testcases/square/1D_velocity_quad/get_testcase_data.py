import subprocess

#-------------------------------------------------------------------------------

def get_testcase_data():

    filenames = ["square_mesh_80x80_culled.nc"]

    dirName = "https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/testcases/square/"

    for filename in filenames:

        args = ["wget", dirName+filename]

        process = subprocess.Popen(args, stdout=subprocess.PIPE)

        while process.poll() is None:
            line = process.stdout.readline()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    get_testcase_data()
