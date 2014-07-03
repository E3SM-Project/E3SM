#!/bin/csh -f

#----------------------------------------

if !( $#argv == 2 ) then
  echo "hist_compare.csh requires two arguments"
  echo "FAIL"
  exit 0
endif

set file1 = $argv[1]
set file2 = $argv[2]

set do_cprnc = 0
if ($?CCSM_CPRNC) then
  if (-e $CCSM_CPRNC) then
    set do_cprnc = 1
  endif
endif

if (${do_cprnc} == 1) then

$CCSM_CPRNC $file1 $file2 >! cprnc.out
tail -13 cprnc.out
set diff_test = `grep "diff_test" cprnc.out | grep IDENTICAL | wc -l`
if (${diff_test} > 0) then
  echo "PASS"
else
  set diff_test = `grep -a "diff_test" cprnc.out | grep IDENTICAL | wc -l`
  if (${diff_test} > 0) then
    echo "PASS"
  else
    echo "FAIL"
  endif
endif

else

set ndir = hist_compare_tmpdir
if (-e $ndir) rm -r -f $ndir
mkdir $ndir
if !(-e $ndir) then
  echo "$ndir not created"
  echo "FAIL"
  exit 0
endif
cd $ndir

ncdump $file1 | grep -v netcdf | grep -v history | split -50000
foreach sfile (x[a-z][a-z])
  mv $sfile $sfile.base
end
echo "ncdump1 done"

ncdump $file2 | grep -v netcdf | grep -v history | split -50000
echo "ncdump2 done"

echo "comparing split files x[a-z][a-z]"
foreach sfile (x[a-z][a-z])
  set sdiff = `diff ${sfile}.base ${sfile} | head -5`
  if !("$sdiff" == "") then
    echo "$sfile"
    echo "$sdiff"
    echo "FAIL"
    cd ../
    rm -r -f $ndir
    exit 0
  endif
end

echo "Files appear identical"
echo "PASS"

cd ..
rm -r -f $ndir

endif

exit 0

