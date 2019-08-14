#!/bin/sh
# This is a test script for PIO.
# Ed Hartnett

# Stop execution of script if error is returned.
set -e

# Stop loop if ctrl-c is pressed.
trap exit INT TERM

printf 'running Fortran tests for PIO netCDF integration...\n'
PIO_TESTS='tst_c_pio ftst_pio_orig ftst_pio'

success1=true
for TEST in $PIO_TESTS
do
    success1=false
    echo "running ${TEST}"
    mpiexec -n 4 ./${TEST} && success1=true
    if test $success1 = false; then
        break
    fi
done

# Did we succeed?
if test x$success1 = xtrue; then
    exit 0
fi
exit 1
