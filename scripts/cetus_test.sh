#!/bin/csh -f
setenv LOCARGS "--block ${COBALT_PARTNAME}"
ctest --verbose
mail -s'cetus pio2 tests' jedwards@ucar.edu -a pio2tests.out
