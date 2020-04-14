#!/bin/sh

# This is a test script for PIO. It runs performance tests for the
# netCDF intergration of PIO.

# Ed Hartnett

# Stop execution of script if error is returned.
set -e

# Stop loop if ctrl-c is pressed.
trap exit INT TERM

printf 'running PIO performance tests...\n'

PIO_TESTS='tst_ncint_perf tst_ncint_async_perf'

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
