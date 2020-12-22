# Fortran FFT Testing workflow

The FFT tests are not really "tests", but rather a collection of simple scripts to serve as sanity checks on how the various FFT routines in the CRM code(s) behave. This include simple forward and backwards transforms and print statements to verify the order of FFT weights.

The following commands can be used to test the Fortran FFTs:
```bash
################################################################
################################################################

# clean previous build
rm -f test_fft.exe fftpack51d.mod

# build
gfortran ../sam/fftpack5_1d.F90 test_fft.F90 -o test_fft.exe 

#execute
./test_fft.exe


```
