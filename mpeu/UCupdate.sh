#!/bin/sh
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
#BOP -------------------------------------------------------------------
#
# !ROUTINE: UCupdate.sh - Find and copy F90 modules to a directory
#
# !DESCRIPTION:
#
#   This version handles modules with upper-case names.
#
# !INTERFACE:
#
#   [env M=<sfx>] sh UCupdate.sh <dir> <mod.o>+
#
# !REVISION HISTORY:
# 	04Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
#		- initial prototype/prolog/code
#EOP ___________________________________________________________________

M=${M-:mod}
c=`basename $0 .sh`

if [ $# -le 1 ]; then
  echo "Usage: [env M=<sfx>] $c <dir> <modfile>+"
  exit 1
fi

D=$1
shift 1

if [ ! -d $D ]; then
  mkdir $D
  if [ $? != 0 ]; then echo "$c: cann't mkdir, $D" 1>&2; exit 2; fi
fi

remove(){
  echo "$c: removing the target, $D" 1>&2
  rm -fr ${D}
}

trash="${D}"
trap "remove; trap ''  1; exit 2" 1
trap "remove; trap ''  2; exit 2" 2
trap "remove; trap ''  3; exit 2" 3
trap "remove; trap '' 15; exit 1" 15

status=0
for o in $@; do
  case $o in
  *.o)
    F=`basename $o .o | tr 'a-z' 'A-Z'`.$M
    ;;
  *)
    F="$o"
    ;;
  esac

  if [ ! -r $F ]; then echo "$c: not found, $F" 1>&2; status=1; fi

	# If two files are the same, do nothing

  if cmp -s $F $D/$F; then
    echo "\c"
  else
    rm -f $D/$F
    cp -p $F $D/
  fi
done
exit $status
