#! /bin/csh -fv

if !(-d $OBJROOT/glc/obj ) mkdir -p $OBJROOT/glc/obj || exit 2
if !(-d $OBJROOT/glc/source) mkdir -p $OBJROOT/glc/source || exit 3
if !(-d $OBJROOT/glc/input ) mkdir -p $OBJROOT/glc/input || exit 4

set my_path = $CASEROOT/SourceMods/src.mpas-li

echo -----------------------------------------------------------------
echo Copy the necessary files into $OBJROOT/glc/source
echo -----------------------------------------------------------------

cd $OBJROOT/glc/source

cp -fpR $CODEROOT/glc/mpas-li/model/src/* .
cp -fpR $CODEROOT/glc/mpas-li/driver glc_cesm_driver

if ( $?MACH ) then
        if ( "X$MACH" == "Xedison" ) then
                make all CORE=landice ESM=ACME DRIVER=glc_cesm_driver TOOL_TARGET_ARCH="-target-cpu=sandybridge" GEN_F90=true ROOT_DIR=`pwd` || exit 5
        else
                make all CORE=landice ESM=ACME DRIVER=glc_cesm_driver GEN_F90=true ROOT_DIR=`pwd` || exit 5
        endif
else
        exit 5
endif


## COPY ALL MODULE FILES TO THE GLC OBJ DIRECTORY ##
find . -name "*.mod" -exec cp -p {} $OBJROOT/glc/obj/. \;

## COPY LIBGLC TO LIBROOT ##
cp -p libglc.a ${LIBROOT}/libglc.a
