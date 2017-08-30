#! /bin/bash

use_pFUnit=1  # yes
use_pFUnit=0  # no

this_dir=`pwd`

# Path to ACME code

ACME_src_root=${this_dir%components*}
echo ACME src root is $ACME_src_root


# Create a clean src/ directory for hosting soft links to ACME's Fortran code.
# When not using pFUnit, src/ also contains links to the test code files.

path_to_pf_src=`pwd`/src/
echo Src dir for test compilation is $path_to_pf_src

if [  -d $path_to_pf_src ]; then
  rm -rf $path_to_pf_src
fi
mkdir    $path_to_pf_src
cd       $path_to_pf_src

# Create soft links

ln -s $ACME_src_root/cime/src/share/util/shr_kind_mod.F90
ln -s $ACME_src_root/cime/src/share/util/shr_const_mod.F90

ln -s $ACME_src_root/components/cam/src/control/cam_logfile.F90

ln -s $ACME_src_root/components/cam/src/physics/cam/ppgrid.F90
ln -s $ACME_src_root/components/cam/src/physics/cam/constituents.F90
ln -s $ACME_src_root/components/cam/src/physics/cam/glb_verif_smry.F90


# if not using pFUnit, also soft link files in tests/ to src/ so that we have 
# everything in src/ for compilation.
if [ $use_pFUnit -eq 0 ]; then
   ln -s $path_to_pf_src/../tests/* $path_to_pf_src
fi
