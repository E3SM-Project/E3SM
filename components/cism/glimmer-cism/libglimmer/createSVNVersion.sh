#!/bin/bash

# create a sed script containing svn version

OUT=svnversion.sed

echo "s|svn_vers_string|$(svnversion $1)|" > $OUT.tmp

if [ -f $OUT ]; then
  if [ -n "$(diff $OUT $OUT.tmp)" ]; then
     cp $OUT.tmp $OUT
  fi
else
  cp $OUT.tmp $OUT
fi

rm -f $OUT.tmp
