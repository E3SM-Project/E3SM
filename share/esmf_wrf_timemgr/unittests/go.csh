#!/bin/csh

rm -f ./test
gmake
rm -f ./test.out
./test >& test.out

tail -5 test.out
set nd = `diff test.out.base test.out | wc -l`

echo "diffs vs baseline = $nd"



