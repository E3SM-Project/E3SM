#!/bin/csh -f
setenv LOCARGS "--block ${COBALT_PARTNAME}"
ctest --verbose
grep 'tests passed' ./pio2build.out | mail -s'cetus pio2 tests' jedwards@ucar.edu
