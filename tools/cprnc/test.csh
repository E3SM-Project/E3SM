#!/bin/csh -f

gmake

set exefile = ./cprnc

#-- testdata is on bluevista /ptmp
set histdir = /ptmp/tcraig/cprnc.testdata
set fh1 = $histdir/fm.clmhalf.fmi.base.clm2.h0.1998-01-02-00000.nc
set fh2 = $histdir/fm.clmhalf.fmi.clm2.h0.1998-01-02-00000.nc
set fh3 = $histdir/fm.clmhalf.fmi.clm2.h0.1998-01-01-43200.nc
set fr1 = $histdir/fm.clmhalf.fmi.base.clm2.r.1998-01-02-00000.nc
set fr2 = $histdir/fm.clmhalf.fmi.clm2.r.1998-01-02-00000.nc
set fts = $histdir/fm.clmhalf.fmi.clm2.h0.2ts.nc

rm -f hist.out
echo "create hist.out, compare two identical clm history files"
$exefile $fh1 $fh2 >&! hist.out

rm -f difh.out
echo "create difh.out, compare two different clm history files"
$exefile -m $fh3 $fh2 >&! difh.out

rm -f rest.out
echo "create rest.out, compare two identical clm restart files"
$exefile $fr1 $fr2 >&! rest.out

rm -f bigd.out
echo "create bigd.out, compare a history and restart file"
$exefile $fh2 $fr2 >&! bigd.out

rm -f f2ts.out
echo "create f2ts.out, analyze a two timestep history file"
$exefile $fts >& f2ts.out

echo "==============================="
echo " "
tail -10 hist.out
echo " "
echo "==============================="
echo " "
tail -10 difh.out
echo " "
echo "==============================="
echo " "
tail -10 rest.out
echo " "
echo "==============================="
echo " "
tail -10 bigd.out
echo " "
echo "==============================="
echo " "
tail -10 f2ts.out

