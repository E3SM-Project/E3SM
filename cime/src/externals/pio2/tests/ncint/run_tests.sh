#!/bin/sh
# This is a test script for PIO.
# Ed Hartnett

# Stop execution of script if error is returned.
set -e

# Stop loop if ctrl-c is pressed.
trap exit INT TERM

printf 'running PIO tests...\n'

PIO_TESTS='tst_pio_udf tst_pio_async tst_async_multi'

success1=true
success2=true
for TEST in $PIO_TESTS
do
    success1=false
    echo "running ${TEST}"
    mpiexec -n 4 ./${TEST} && success1=true
    if test $success1 = false; then
        break
    fi
done

# PIO_TESTS_8='test_async_multi2  test_async_manyproc'

# for TEST in $PIO_TESTS_8
# do
#     success2=false
#     echo "running ${TEST}"
#     mpiexec -n 8 ./${TEST} && success2=true
#     if test $success2 = false; then
#         break
#     fi
# done

# Did we succeed?
if test x$success1 = xtrue -a x$success2 = xtrue; then
    exit 0
fi
exit 1
