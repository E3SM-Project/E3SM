# C++ FFT Testing workflow

The FFT tests are not really "tests", but rather a collection of simple scripts to serve as sanity checks on how the various FFT routines in the CRM code(s) behave. This include simple forward and backwards transforms and print statements to verify the order of FFT weights.

The following commands can be used to test the YAKL C++ FFTs:
```bash
################################################################
################################################################

# set up environment
module purge
module load cmake xl/16.1.1-5 spectrum-mpi/10.3.1.2-20200121 hsi/5.0.2.p5 xalt/1.2.1 lsf-tools/2.0 darshan-runtime/3.1.7 DefApps

export YAKL_HOME="`pwd`/../../../../../../externals/YAKL"

# clean the previous build
./clean_build.sh

# configure and build
./cmakescript.sh

# execute
./test_fft_cpp.exe


```
