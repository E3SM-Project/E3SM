#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "CAM_compare.sh: incorrect number of input arguments" 
    exit 1
fi

echo "CAM_compare.sh: comparing $1 " 
echo "                with      $2" 

##note syntax here as stderr and stdout from cprnc command go 
##to separate places!
#
# The -d option to cprnc is required for the comparison of sat_hist
# files, it should not affect other files.
#
${CPRNC_EXE} -d ncol:1:-1 $1 $2 2>&1 > cprnc.out
rc=$?
if [ $rc -ne 0 ]; then
    echo "CAM_compare.sh: error doing comparison, cprnc error= $rc" 
    exit 2
fi

if grep -c "the two files seem to be IDENTICAL" cprnc.out > /dev/null; then
    echo "CAM_compare.sh: files are b4b" 
elif grep -c "the two files seem to be DIFFERENT" cprnc.out > /dev/null; then
    echo "CAM_compare.sh: files are NOT b4b...the following fields had diffs" 
    result=`perl -e 'while (my $ll = <>) \
	{ if ($ll =~ /RMS\s+(\S+)\s+(\S+)/) \
	{ print "$1 " }}' cprnc.out`
    echo $result
    exit 3
else
    echo "CAM_compare.sh: unable to determine whether files are identical"
    exit 4
fi

exit 0
