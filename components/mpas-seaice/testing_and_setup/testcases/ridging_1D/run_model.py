import os

#-------------------------------------------------------------------------------

def run_model():

    os.system("rm -rf namelist.seaice streams.seaice")
    os.system("ln -s namelist.seaice.ridging_1D namelist.seaice")
    os.system("ln -s streams.seaice.ridging_1D streams.seaice")

    os.system("../../../seaice_model")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
