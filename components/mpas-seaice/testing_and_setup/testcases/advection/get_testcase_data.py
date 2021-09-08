import subprocess

#-------------------------------------------------------------------------------

def get_testcase_data():

    filenames = ["x1.10242.grid.nc",
                 "x1.163842.grid.nc",
                 "x1.2562.grid.nc",
                 "x1.40962.grid.nc"]

    dirName = "https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/testcases/strain_stress_divergence/"

    for filename in filenames:

        args = ["wget", dirName+filename]

        process = subprocess.Popen(args, stdout=subprocess.PIPE)

        while process.poll() is None:
            line = process.stdout.readline()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    get_testcase_data()
