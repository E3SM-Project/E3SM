# This is a test script for PIO examples.
# Ed Hartnett 5/7/18

# Stop execution of script if error is returned.
set -e

# Stop loop if ctrl-c is pressed.
trap exit SIGINT SIGTERM

printf 'running PIO examples...\n'

#PIO_EXAMPLES='examplePio'
PIO_EXAMPLES='example1 examplePio darray_no_async'

success1=true
for EXAMPLE in $PIO_EXAMPLES
do
    success1=false
    echo "running ${EXAMPLE}"
    mpiexec -n 4 ./${EXAMPLE} && success1=true || break
done

# Did we succeed?
if test x$success1 = xtrue; then
    exit 0
fi
exit 1
