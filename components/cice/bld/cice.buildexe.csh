#! /bin/csh -f 

set objdir = $OBJROOT/ice/obj

cd $objdir

cp $CASEBUILD/ciceconf/Filepath ./tmp_filepath 
if (-f Filepath) then
  cmp -s tmp_filepath Filepath || mv -f tmp_filepath Filepath 
else
  mv -f tmp_filepath Filepath 
endif

# Check for recompile if BLCKX, BLCKY or MXBLCKS changes
#-------------------------------------------------------
set recompile = FALSE
echo ${CICE_BLCKX} ${CICE_BLCKY} ${CICE_MXBLCKS} > $objdir/iceres.new
cmp -s $objdir/iceres.new $objdir/iceres.old || set recompile = TRUE
if ($recompile == 'TRUE') then
  cat $objdir/iceres.old
  cat $objdir/iceres.new
  echo "blckx,blcky,mxblcks has changed, removing objdir and cice, preparing for new compile"
  rm -f $objdir/*.o
  rm -f $objdir/*.f
  rm -f $objdir/*.f90
endif

# Build the library
#-------------------------------------------------------
set cicedefs = "`cat $CASEBUILD/ciceconf/CICE_cppdefs`"
set cicedefs = "$cicedefs -DBLCKX=$CICE_BLCKX -DBLCKY=$CICE_BLCKY -DMXBLCKS=$CICE_MXBLCKS"        
gmake complib -j $GMAKE_J MODEL=cice COMPLIB=$LIBROOT/libice.a MACFILE=$CASEROOT/Macros.$MACH USER_CPPDEFS="$cicedefs" -f $CASETOOLS/Makefile || exit 2

mv $objdir/iceres.new $objdir/iceres.old

wait 
exit 0


