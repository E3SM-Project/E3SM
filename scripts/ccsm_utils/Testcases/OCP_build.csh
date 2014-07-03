#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT
touch pop_perf.out

./xmlchange -file env_build.xml -id POP_AUTO_DECOMP  -val false

# set precheck to 0 to build cesm
# set precheck to 1 to diagnose decomp configurations
set precheck = 0

set decvals = (cartesian blockone spacecurve)
@ ocn_pes = ${NTASKS_OCN} * ${NTHRDS_OCN}
set cnt = 0

set by = 0
while (${by} < ${OCN_NY})
  @ by = $by + 1
  #--- only check non padded block sizes
  if (${OCN_NY} % $by == 0) then

set bx = 0
while (${bx} < ${OCN_NX})
  @ bx = $bx + 1
  #--- only check non padded block sizes
  if (${OCN_NX} % $bx == 0) then

  set runit = 1
  #--- make sure number of blocks is close to number of procs
  @ nblocks = (${OCN_NX} * ${OCN_NY}) / (${bx} * ${by})
  @ nbmin = 1 * ${ocn_pes} / 3
  @ nbmax = 5 * ${ocn_pes}
  if (${nblocks} < ${nbmin}) set runit = 0
  if (${nblocks} > ${nbmax}) set runit = 0
  #--- skip high aspect ratio blocks
  @ rata = $bx / $by
  @ ratb = $by / $bx
  if (${rata} > 3) set runit = 0
  if (${ratb} > 3) set runit = 0

  if ($runit == 1) then

  @ mxt  = (${OCN_NX} * ${OCN_NY}) / ($bx * $by * ${NTASKS_OCN})
  @ mxtr = (${OCN_NX} * ${OCN_NY}) % ($bx * $by * ${NTASKS_OCN})
  if ($mxtr != 0) then
     @ mxt = $mxt + 1
  endif
  ./xmlchange -file env_build.xml -id POP_BLCKX   -val $bx
  ./xmlchange -file env_build.xml -id POP_BLCKY   -val $by
  ./xmlchange -file env_build.xml -id POP_MXBLCKS -val $mxt
  rm LockedFiles/env_build*  >& /dev/null

  foreach decomp ($decvals)
    set runit = 1
    #--- cartesian must divide grid nice and even
    if ($decomp == cartesian) then
       if ($mxtr != 0) set runit = 0
    endif
    #--- spacecurve must have 2, 3, and 5 only in block factors
    if ($decomp == spacecurve) then
       @ blkx = ${OCN_NX} / $bx
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
       @ blky = ${OCN_NY} / $by
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

    if ($runit == 1) then
      ./xmlchange -file env_build.xml -id POP_DECOMPTYPE -val $decomp
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
#        source ./Tools/ccsm_getenv || exit -2
        echo "precheck... $OCN_GRID ${ocn_pes} $NTASKS_OCN $NTHRDS_OCN $bx $by $mxt $decomp"
      endif

    else
#      source ./Tools/ccsm_getenv || exit -2
      if ($precheck == 0) then
        echo "  skip..... $OCN_GRID ${ocn_pes} $NTASKS_OCN $NTHRDS_OCN $bx $by $mxt $decomp" >> pop_perf.out
      else
        echo "  skip..... $OCN_GRID ${ocn_pes} $NTASKS_OCN $NTHRDS_OCN $bx $by $mxt $decomp" 
      endif
    endif

  end   # decomp
  endif # runit
endif # bx size
end   # bx
endif # by size
end   # by

