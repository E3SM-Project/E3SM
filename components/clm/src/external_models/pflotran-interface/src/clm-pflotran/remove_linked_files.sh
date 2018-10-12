#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    rm -f ${file##*/}
done

