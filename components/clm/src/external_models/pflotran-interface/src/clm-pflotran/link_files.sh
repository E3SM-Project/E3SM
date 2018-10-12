#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    ln -sf ${file} ${file##*/}
done

if [ -f pflotran.F90 ]; then \
   unlink pflotran.F90
fi

if [ -f pflotran_derivative.F90 ]; then \
   unlink pflotran_derivative.F90
fi

if [ -f pflotran_provenance.F90 ]; then \
   unlink pflotran_provenance.F90
fi

