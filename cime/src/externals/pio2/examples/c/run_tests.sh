#!/bin/sh
# This is a test script for PIO examples.
# Ed Hartnett 5/7/18

# Stop execution of script if error is returned.
set -e

# Stop loop if ctrl-c is pressed.
trap exit INT TERM

printf 'running PIO examples...\n'

#PIO_EXAMPLES='examplePio'
PIO_EXAMPLES='example1 examplePio'
PIO_EXAMPLES_16='darray_no_async'

success1=true
for EXAMPLE in $PIO_EXAMPLES
do
    success1=false
    echo "running ${EXAMPLE}"
    mpiexec -n 4 ./${EXAMPLE} && success1=true
    if test $success1 = false; then
        break
    fi
done
success2=true
for EXAMPLE in $PIO_EXAMPLES_16
do
    success2=false
    echo "running ${EXAMPLE}"
    mpiexec -n 16 ./${EXAMPLE} && success2=true
    if test $success2 = false; then
        break
    fi
done

# Did we succeed?
if test x$success1 = xtrue -a  x$success2 = xtrue; then
    exit 0
fi
exit 1
