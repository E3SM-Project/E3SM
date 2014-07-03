#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT
touch cice_perf.out

./xmlchange -file env_build.xml -id CICE_AUTO_DECOMP  -val false

# set precheck to 0 to build cesm
# set precheck to 1 to diagnose decomp configurations
set precheck = 0

set cnt = 0
foreach dlist (1 2)

if (${ICE_GRID} =~ gx1* ) then
  if ($dlist == 1) then
    set decvals = (blkrobin blkcart roundrobin spacecurve)
    set bxvals = ( 2  4  5  8  10  16  20  32)
    set byvals = ( 2  4  6  8  12  16  24  32)
  endif
  if ($dlist == 2) then
    set decvals = (slenderX1 slenderX2)
    set bxvals = (  1   2  4  5  8 10 16 20 32 64)
    set byvals = (384 192 96 48 24 12)
  endif
else if (${ICE_GRID} =~ gx3* ) then
  if ($dlist == 1) then
    set decvals = (blkrobin blkcart roundrobin spacecurve)
    set bxvals = ( 5  10  20)
    set byvals = ( 4  29)
  endif
  if ($dlist == 2) then
    set decvals = (slenderX1 slenderX2)
    set bxvals = (  1  2  4  5 10  20  25)
    set byvals = (116 58 29)
  endif
else if (${ICE_GRID} =~ tx1* ) then
  if ($dlist == 1) then
    set decvals = (blkrobin blkcart roundrobin spacecurve)
    set bxvals = ( 5  8  10  12  15)
    set byvals = ( 6  8  10  12  15)
  endif
  if ($dlist == 2) then
    set decvals = (slenderX1 slenderX2)
    set bxvals = (  1  2  3  4  5  8  9  10  12  15  18  20  24  30)
    set byvals = (240 120 60 30)
  endif
else if (${ICE_GRID} =~ tx0.1* ) then
  if ($dlist == 1) then
    set decvals = (blkrobin blkcart roundrobin spacecurve)
    set bxvals = ( 5  8  10  12  15  20  24)
    set byvals = ( 6  8  10  12  15  20  24)
  endif
  if ($dlist == 2) then
    set decvals = (slenderX1 slenderX2)
    set bxvals = (  1  2  3  4  5  8  9  10  12  15  18  20  24  30  36  40  48  50  60 )
    set byvals = (2400 1200 600 300)
  endif
else
  echo "ICE GRID not supported by default testing in ICP test"
  exit -9
endif
 
foreach by ($byvals)
foreach bx ($bxvals)
  @ mxt  = (${ICE_NX} * ${ICE_NY}) / ($bx * $by * ${NTASKS_ICE})
  @ mxtr = (${ICE_NX} * ${ICE_NY}) % ($bx * $by * ${NTASKS_ICE})
  if ($mxtr != 0) then
     @ mxt = $mxt + 1
  endif
  ./xmlchange -file env_build.xml -id CICE_BLCKX   -val $bx
  ./xmlchange -file env_build.xml -id CICE_BLCKY   -val $by
  ./xmlchange -file env_build.xml -id CICE_MXBLCKS -val $mxt
  rm LockedFiles/env_build*  >& /dev/null

  foreach decomp ($decvals)
    set runit = 1
    #--- only check non padded block sizes
    if (${ICE_NX} % $bx != 0) set runit = 0
    if (${ICE_NY} % $by != 0) set runit = 0
    #--- spacecurve must have 2, 3, and 5 only in block factors
    if ($decomp == spacecurve) then
       @ blkx = ${ICE_NX} / $bx
       foreach divid (2 3 5)
          set done = 0
          while ($done == 0)
             if ($blkx % $divid == 0) then
                @ blkx = $blkx / $divid
             else
                set done = 1
             endif
          end
       end
       if ($blkx != 1) set runit = 0
       @ blky = ${ICE_NY} / $by
       foreach divid (2 3 5)
          set done = 0
          while ($done == 0)
             if ($blky % $divid == 0) then
                @ blky = $blky / $divid
             else
                set done = 1
             endif
          end
       end
       if ($blky != 1) set runit = 0
    endif
    #--- blkcart must have multiple of 4 blocks per task
    if ($decomp == blkcart) then
       @ blk4 = (${ICE_NX} * ${ICE_NY}) % ($bx * $by * 4 * ${NTASKS_ICE})
       if ($blk4 != 0) then
          set runit=0
       endif
    endif
    #--- slenderX1 is constrained by certain blocksizes
    if ($decomp == slenderX1) then
       @ blkx = (${ICE_NX}) % ($bx * ${NTASKS_ICE})
       @ blky = (${ICE_NY}) % ($by)
       if ($blkx != 0 || $blky != 0) then
          set runit=0
       endif
    endif
    #--- slenderX2 is constrained by certain blocksizes
    if ($decomp == slenderX2) then
       @ blkx = (${ICE_NX} * 2) % ($bx * ${NTASKS_ICE})
       @ blky = (${ICE_NY}) % ($by * 2)
       if ($blkx != 0 || $blky != 0) then
          set runit=0
       endif
    endif

    @ ice_pes = ${NTASKS_ICE} * ${NTHRDS_ICE}
    if ($runit == 1) then
      ./xmlchange -file env_build.xml -id CICE_DECOMPTYPE -val $decomp
      ./xmlchange -file env_build.xml -id CICE_DECOMPSETTING -val null
      if ($decomp == slenderX1 || $decomp == slenderX2 || $decomp == square-pop) then
        ./xmlchange -file env_build.xml -id CICE_DECOMPTYPE -val cartesian
        ./xmlchange -file env_build.xml -id CICE_DECOMPSETTING -val $decomp
      endif
      cp env_build.xml LockedFiles/env_build.xml.locked

      if ($precheck == 0) then
        ./*.build

        @ cnt = $cnt + 1
        set acnt = $cnt
        if ($cnt < 1000) set acnt = "0${cnt}"
        if ($cnt < 100 ) set acnt = "00${cnt}"
        if ($cnt < 10  ) set acnt = "000${cnt}"

        cp -f env_build.xml env_build.xml.${acnt}
        cp -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.${acnt}
      else
        source ./Tools/ccsm_getenv || exit -2
        echo "precheck... $ICE_GRID ${ice_pes} $NTASKS_ICE $NTHRDS_ICE $CICE_BLCKX $CICE_BLCKY $CICE_MXBLCKS $CICE_DECOMPTYPE $CICE_DECOMPSETTING"
      endif

    else
      source ./Tools/ccsm_getenv || exit -2
      if ($precheck == 0) then
        echo "  skip..... $ICE_GRID ${ice_pes} $NTASKS_ICE $NTHRDS_ICE $bx $by $mxt $decomp" >> cice_perf.out
      else
        echo "  skip..... $ICE_GRID ${ice_pes} $NTASKS_ICE $NTHRDS_ICE $bx $by $mxt $decomp" 
      endif
    endif

  end   # decomp
end   # bxvals
end   # byvals
end   # list

