#!/bin/bash

SCRIPT=`cut -f 1 -d" " args.in`
OUT=`cut -f 2 -d" " args.in`
ARGUM=`cut -f 3- -d" " args.in`

THIS=`basename $PWD`
echo "Running $SCRIPT $ARGUM in $THIS"

#echo $(< args.in)
../$SCRIPT $ARGUM > $OUT
