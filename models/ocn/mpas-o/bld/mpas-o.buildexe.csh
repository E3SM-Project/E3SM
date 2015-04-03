#! /bin/csh -fv

if !(-d $OBJROOT/ocn/obj   ) mkdir -p $OBJROOT/ocn/obj    || exit 2
if !(-d $OBJROOT/ocn/source) mkdir -p $OBJROOT/ocn/source || exit 3 
if !(-d $OBJROOT/ocn/input ) mkdir -p $OBJROOT/ocn/input  || exit 4

set my_path = $CASEROOT/SourceMods/src.mpas-o

echo -----------------------------------------------------------------
echo  Copy the necessary files into $OBJROOT/ocn/source
echo -----------------------------------------------------------------

cd $OBJROOT/ocn/source

cp -fpR $CODEROOT/ocn/mpas-o/model/src/* .
cp -fpR $CODEROOT/ocn/mpas-o/driver ocean_cesm_driver

if ( $?MACH ) then
	if ( "X$MACH" == "Xedison" ) then
		make all CORE=ocean ESM=ACME DRIVER=ocean_cesm_driver TOOL_TARGET_ARCH="-target-cpu=sandybridge" GEN_F90=true ROOT_DIR=`pwd` || exit 5
	else
		make all CORE=ocean ESM=ACME DRIVER=ocean_cesm_driver GEN_F90=true ROOT_DIR=`pwd` || exit 5
	endif
else
	exit 5
endif

## COPY ALL MODULE FILES TO THE OCEAN OBJ DIRECTORY ##
find . -name "*.mod" -exec cp -p {} $OBJROOT/ocn/obj/. \;

## COPY LIBOCEAN TO LIBOCN IN LIBROOT ##
cp -p libocn.a ${LIBROOT}/libocn.a

