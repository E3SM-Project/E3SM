#! /bin/csh -fv

if !(-d $OBJROOT/ocn/obj   ) mkdir -p $OBJROOT/ocn/obj    || exit 2
if !(-d $OBJROOT/ocn/source) mkdir -p $OBJROOT/ocn/source || exit 3 
if !(-d $OBJROOT/ocn/input ) mkdir -p $OBJROOT/ocn/input  || exit 4

set my_path = $CASEROOT/SourceMods/src.pop2

echo -----------------------------------------------------------------
echo  Copy the necessary files into $OBJROOT/ocn/source
echo -----------------------------------------------------------------

cd $OBJROOT/ocn/source

cp -fp $CODEROOT/ocn/pop2/source/*.F90  .
cp -fp $CODEROOT/ocn/pop2/mpi/*.F90  .
cp -fp $CODEROOT/ocn/pop2/drivers/cpl/*.F90 .

# Two files require special attention because they get renamed when copied
# from either SourceMods/src.pop2 or input_templates/:
# ${OCN_GRID}_domain_size.F90       -> domain_size.F90
# ${OCN_GRID}_POP_DomainSizeMod.F90 -> POP_DomainSizeMod.F90
# (The latter is "needed during LANL merge transition")
#
# For these files:
# 1) Make sure SourceMods does not contain copies of the same file both
#    with and without the ${OCN_GRID}_ preface.
if ((-f ${my_path}/domain_size.F90) && \
    (-f ${my_path}/${OCN_GRID}_domain_size.F90)) then
   echo "ERROR: you can not have both domain_size.F90 and " \
        "${OCN_GRID}_domain_size.F90 in SourceMods/src.pop2/"
   exit -1
endif

if ((-f ${my_path}/POP_DomainSizeMod.F90) && \
    (-f ${my_path}/${OCN_GRID}_POP_DomainSizeMod.F90)) then
   echo "ERROR: you can not have both POP_DomainSizeMod.F90 and " \
        "${OCN_GRID}_POP_DomainSizeMod.F90 in SourceMods/src.pop2/"
   exit -1
endif

set POP2_FOUND_d_s = 0
set POP2_FOUND_POP_DSM = 0
# 1) Copy (with name-change) from input templates/ if they exist
if (-f  $CODEROOT/ocn/pop2/input_templates/${OCN_GRID}_domain_size.F90 ) then
   cp -fp $CODEROOT/ocn/pop2/input_templates/${OCN_GRID}_domain_size.F90 domain_size.F90
   set POP2_FOUND_d_s = 1
endif
if (-f $CODEROOT/ocn/pop2/input_templates/${OCN_GRID}_POP_DomainSizeMod.F90) then 
   cp -fp $CODEROOT/ocn/pop2/input_templates/${OCN_GRID}_POP_DomainSizeMod.F90 POP_DomainSizeMod.F90
   set POP2_FOUND_POP_DSM = 1
endif

# 2) Copy everything from SourceMods over
if (-d $my_path ) cp -fp $my_path/*.F90 .
# If domain_size.F90 or POP_DomainSizeMod.F90 exist, they should overwrite 
# anything copied from from input_templates/
if (-f  ${my_path}/domain_size.F90 ) then
   set POP2_FOUND_d_s = 1
endif
if (-f ${my_path}/POP_DomainSizeMod.F90) then 
   set POP2_FOUND_POP_DSM = 1
endif

# 3) If SourceMods/ contains ${OCN_GRID}_domain_size.F90 or
#    ${OCN_GRID}_POP_DomainSizeMod.F90, those files will now exist in the
#    current directory and need to be renamed
if (-f  ${OCN_GRID}_domain_size.F90 ) then
   mv -f ${OCN_GRID}_domain_size.F90 domain_size.F90
   set POP2_FOUND_d_s = 1
endif
if (-f ${OCN_GRID}_POP_DomainSizeMod.F90) then
   mv -f ${OCN_GRID}_POP_DomainSizeMod.F90 POP_DomainSizeMod.F90
   set POP2_FOUND_POP_DSM = 1
endif

# 4) Make sure both domain_size.F90 and POP_DomainSizeMod.F90 exist for the
#    specified grid
if (${POP2_FOUND_d_s} == 0) then
  echo "ERROR: you need either ${OCN_GRID}_domain_size.F90 or domain_size.F90"
  exit -1
endif
if (${POP2_FOUND_POP_DSM} == 0) then
  echo "ERROR: you need either ${OCN_GRID}_POP_DomainSizeMod.F90 or " \
       "POP_DomainSizeMod.F90"
  exit -1
endif

echo -------------------------------------------------------------------------
echo  Checking for any auxilliary ocean-component models before building pop2
echo  For now - only looking for moby component
echo -------------------------------------------------------------------------

setenv USE_OCN_MOBY FALSE
foreach comp (`echo $OCN_TRACER_MODULES`)
  if ($comp == moby) then
    setenv USE_OCN_MOBY TRUE
  endif
end

if ($USE_OCN_MOBY == TRUE ) then

  echo -------------------------------------------------------------------------
  echo  Building moby
  echo -------------------------------------------------------------------------

  cd $CASEBUILD
  if (-f moby.buildexe.csh) then
    ./moby.buildexe.csh
    if ( $status != 0 ) then
      echo ERROR: moby.buildexe.csh failed
      exit 18
    endif
  else 
    cp $CODEROOT/ocn/pop2/aux/moby/pop2/bld/moby.buildexe.csh .
    ./moby.buildexe.csh
    if ( $status != 0 ) then
      echo ERROR: moby.buildexe.csh failed
      exit 18
    endif
  endif
endif

echo -----------------------------------------------------------------
echo  Build pop2 library
echo -----------------------------------------------------------------

cd ${OBJROOT}/ocn/obj 

cat >! Filepath <<EOF
$OBJROOT/ocn/source
EOF

@ NT = 2
foreach module ( `echo $OCN_TRACER_MODULES` )  
   if ($module =~ "iage"   ) @ NT = $NT +  1
   if ($module =~ "cfc"    ) @ NT = $NT +  2   
   if ($module =~ "ecosys" ) @ NT = $NT + 27
   if ($module == moby     ) then
      if (-e $my_path/${OCN_GRID}_data.ptracers) then
        set dir = $my_path
      else if (-e $CODEROOT/ocn/pop2/aux/moby/darwin/input/${OCN_GRID}_data.ptracers) then
        set dir = $CODEROOT/ocn/pop2/aux/moby/darwin/input
      else
        exit 31
      endif
      set nt_moby = `grep PTRACERS_numInUse $dir/${OCN_GRID}_data.ptracers | cut -f 2 -d = | cut -f 1 -d","`
      if ($status != 0) exit 32
      @ NT = $NT + $nt_moby
    endif
end

set cppdefs = "-DCCSMCOUPLED -DBLCKX=$POP_BLCKX -DBLCKY=$POP_BLCKY -DMXBLCKS=$POP_MXBLCKS -DNT=$NT"
if ($OCN_ICE_FORCING == 'inactive' ) set cppdefs = "$cppdefs -DZERO_SEA_ICE_REF_SAL"
if ($OCN_GRID =~ "tx0.1*"          ) set cppdefs = "$cppdefs -D_HIRES";
if ($OCN_ICE_FORCING == 'inactive' ) set cppdefs = "$cppdefs -DZERO_SEA_ICE_REF_SAL"

cat >! $OBJROOT/ocn/obj/POP2_cppdefs.new <<EOF
$cppdefs
EOF

#  recompile if 2d decomp or NT is changed
set recompile = FALSE
if (-e $OBJROOT/ocn/obj/POP2_cppdefs) then
    diff $OBJROOT/ocn/obj/POP2_cppdefs.new $OBJROOT/ocn/obj/POP2_cppdefs || set recompile = TRUE
    echo "recompile is $recompile"
    if ($recompile == 'TRUE') then
      touch `grep -wl BLCKX   $OBJROOT/ocn/source/*`  # force recompile
      touch `grep -wl BLCKY   $OBJROOT/ocn/source/*`  # force recompile
      touch `grep -wl MXBLCKS $OBJROOT/ocn/source/*`  # force recompile
      touch `grep -wl NT      $OBJROOT/ocn/source/*`  # force recompile
    endif  
endif
cp -f $OBJROOT/ocn/obj/POP2_cppdefs.new $OBJROOT/ocn/obj/POP2_cppdefs
set pop2defs = "`cat $OBJROOT/ocn/obj/POP2_cppdefs`"

gmake complib -j $GMAKE_J MODEL=pop2 COMPLIB=$LIBROOT/libocn.a USER_CPPDEFS="$pop2defs" -f $CASETOOLS/Makefile || exit 2

echo " "
echo ----------------------------------------------------------------------------
echo  Note that f90 files may not exist on all machines
echo ----------------------------------------------------------------------------

set f90_dir = $OBJROOT/ocn/source/f90
if !(-d  $f90_dir ) mkdir -p $f90_dir
mv -f *.f90 $f90_dir

if !(-f $LIBROOT/libocn.a) then
  echo "ERROR: pop2 library not available"
  exit -1
endif

echo " "
echo -------------------------------------------------------------------------
echo  Successful completion of the pop2 executable building process
echo -------------------------------------------------------------------------
