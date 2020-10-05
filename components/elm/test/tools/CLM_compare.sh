#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "CLM_compare.sh: incorrect number of input arguments" 
    exit 1
fi

echo "CLM_compare.sh: comparing $1 " 
echo "                with      $2" 

##note syntax here as stderr and stdout from cprnc command go 
##to separate places!
${CPRNC_EXE} $1 $2 2>&1 > cprnc.out
rc=$?
if [ $rc -ne 0 ]; then
    echo "CLM_compare.sh: error doing comparison, cprnc error= $rc" 
    exit 2
fi

result_old=`perl -e 'while (my $ll = <>) \
    { if ($ll =~ /(\d+)[^0-9]+compared[^0-9]+(\d+)/) \
    { print "PASS" if $1>0 && $2==0 }}' cprnc.out`
if grep -c "the two files seem to be IDENTICAL" cprnc.out > /dev/null; then
   result=PASS
elif grep -c "the two files seem to be DIFFERENT" cprnc.out > /dev/null; then
   result=FAIL
else
   result=$result_old
fi

if [ "$result" = "PASS" ]; then
    echo "CLM_compare.sh: files are b4b" 
else
    echo "CLM_compare.sh: files are NOT b4b" 
    exit 3
fi

exit 0
