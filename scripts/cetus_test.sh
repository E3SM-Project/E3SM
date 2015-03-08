#!/bin/csh -f
setenv LOCARGS "--block ${COBALT_PARTNAME}"
ctest --verbose
