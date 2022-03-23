#!/bin/bash
#
# Configure and run all DCMIP tests and then collect PDF output
# for use in a PDF HOMME evaluation document.
#
# Run the script from the git clone to be tested.  The script location is used to determine
# the source code directory
#
# Run this script no options to see the syntax
#
# options:
#   clean:     remove CMake files (to reconfigure)
#   configure: configure with cmake
#   submit:    compile each executable, submit jobscript
#   status:    see if PDF files for all tests are present
#   latex:     collect all PDF files, create PDF of results
#
# Currently running in Anvil
#
# Requirements:
#    standalone HOMME supported platform
#    jobscript for each CASE that runs all tests and produces plots
#    jobscript needs to produce plots with the exact names given below
#    NCL, pdflatex 
#
MACH=anvil.cmake            # HOMME's machine file for the given machine
JOBSCRIPT=jobscript-snl.sh  # name of the jobscript to be submitted in every test

# list of DCMIP tests to run
declare -a cases
cases[0]="dcmip_tests/dcmip2012_test2.0_steady_state_with_orography/theta-l"
cases[1]="dcmip_tests/dcmip2012_test2.1_nh_mountain_waves_no_shear/theta-l"
cases[2]="dcmip_tests/dcmip2012_test2.2_nh_mountain_waves_with_shear/theta-l" 
cases[3]="dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta-l"
cases[4]="dcmip_tests/dcmip2012_test4.1_baroclinic_instability/theta-l"
cases[5]="dcmip_tests/dcmip2016_test1_baroclinic_wave/theta-l"
cases[6]="dcmip_tests/dcmip2016_test2_tropical_cyclone/theta-l"
cases[7]="dcmip_tests/dcmip2016_test3_supercell/theta-l"

# list of output files for each DCMIP test to collect:
declare -a outputs
outputs[0]="hydro_test2_0_u.pdf hydro_test2_0_u_z.pdf nonhydro_test2_0_u.pdf nonhydro_test2_0_u_t6.00.pdf"
outputs[1]="hydro_test2_1_T_t10.pdf nonhydro_test2_1_T_t10.pdf"
outputs[2]="hydro_test2_2_T_t10.pdf nonhydro_test2_2_T_t10.pdf"
outputs[3]="hydro_test31_omega.pdf hydro_test3_theta_diff_last.pdf nonhydro_test31_omega.pdf nonhydro_test3_theta_diff_last.pdf"
outputs[4]="hydro-X1-ps.pdf hydro-X1-zeta.pdf nonhydro-X1000-ps.pdf nonhydro-X1000-zeta.pdf \
nonhydro-X100-ps.pdf nonhydro-X100-zeta.pdf \
nonhydro-X10-ps.pdf nonhydro-X10-zeta.pdf \
nonhydro-X1-ps.pdf nonhydro-X1-zeta.pdf \
"
outputs[5]="r400_PS.pdf r400_T850.pdf r100-dry_PS.pdf r100-dry_T850.pdf \
r100-h_PS.pdf r100-h_T850.pdf r100_PS.pdf r100_T850.pdf r50_PS.pdf r50_T850.pdf"
outputs[6]="r400_psmap.pdf r100_psmap.pdf r50_psmap.pdf r50-h_psmap.pdf"
outputs[7]="5km_xsec_r400.pdf  5km_xsec_r200.pdf  5km_xsec_r100.pdf  5km_xsec_r50.pdf"


args=("$@")
if [ "$#" -lt "2" ]; then
    echo "To configure: "
    echo "./homme.sh /path/to/output/dir {clean,configure}"
    echo ""
    echo "To run all tests, check status or collect plots: "
    echo "./homme.sh /path/to/output/dir {submit,status,latex}"
    echo ""
    echo "As above, but for a single single test (0-7):"
    echo "./homme.sh /path/to/output/dir {submit,status,latex} testno"
    exit 1
fi

if ! which ncl >&/dev/null > /dev/null ; then
   echo "'which ncl' failed. NCL is required to produce PDFs. Exiting."
   exit 1
fi


WDIR=${args[0]}
opt=${args[1]}

if [ "$#" -eq "2" ]; then
    len=${#cases[@]}
    case_start=0
    case_stop=len-1
else
    case_start=${args[2]}
    case_stop=${args[2]}
fi

SRCDIR="`dirname \"$0\"`"              # relative
SRCDIR="`( cd \"$SRCDIR/..\" && pwd )`"  # absolutized and normalized
echo SOURCE CODE: $SRCDIR
echo OUTPUT DIR:  $WDIR
if [ -z "$SRCDIR" ] ; then
  echo ERROR: for some reason cannot find src directory. Exiting.
  exit 1  # fail
fi

cd $SRCDIR
if [ ! -d dcmip_tests ]; then 
    echo Error: SRCDIR is missing dcmip tests case directory
    exit 1
fi
mkdir -p $WDIR
cd $WDIR
if [ -f README.cmake ] || [ -d cime ]; then 
    # prevent user from accidently installing in the git clone
    echo Error: /path/to/output/dir appears to be the source code directory. Exiting.
    exit 1
fi


if [ "$opt" = "clean" ]; then
  cd $WDIR
  \rm -Rf CMakeFiles CMakeCache.txt
fi

if [ "$opt" = "configure" ]; then
  cd $WDIR
  echo running cmake:
  cmake -C $SRCDIR/cmake/machineFiles/$MACH \
  -DUSE_NUM_PROCS=24 -DPREQX_PLEV=26 -DPREQX_NP=4 -DPREQX_USE_ENERGY=TRUE  \
   $SRCDIR
  if [ $? -ne 0 ]; then exit $?; fi
  make -j4 theta-l  # make sure theta model compiles
  if [ $? -ne 0 ]; then exit $?; fi
fi


if [ "$opt" = "submit" ]; then
    if [ ! -x src/theta-l/theta-l ]; then
        echo Error: theta-l not built during configure step. Exiting.
        exit 1
    fi
   for (( i=$case_start; i<=$case_stop; i++ )); do 
       echo BUILD: ${cases[i]}
       cd $WDIR/${cases[i]}
       if [ $? -ne 0 ]; then exit $?; fi
       make install
       ./build.sh
       if [ $? -ne 0 ]; then exit $?; fi

       # remove any existing outputs:
       IFS=" " read -a ARRAY <<< "${outputs[i]}"
       for file in "${ARRAY[@]}"
       do
           \rm -f $file
       done

       echo SUBMIT: ${cases[i]}
       cd $WDIR/${cases[i]}
       sbatch $JOBSCRIPT
   done
fi

missing_count=0
if [ "$opt" = "status" ] || [ "$opt" = "latex" ]; then
   for (( i=$case_start; i<=$case_stop; i++ )); do 
       echo CASE $i: ${cases[i]}
       cd $WDIR/${cases[i]}
       if [ $? -ne 0 ]; then exit $?; fi

       IFS=" " read -a ARRAY <<< "${outputs[i]}"
       for file in "${ARRAY[@]}"
       do
           if [ -f "$file" ]; then
              echo "  found: $file"
           else
              (( missing_count++ ))
              echo "  missing(${missing_count}): $file"
           fi
       done
   done
   echo missing output count: $missing_count
fi


if [ "$opt" = "latex" ]; then
    if [ $missing_count -ne 0 ]; then 
        echo ERROR: some output files are missing.  Exiting.
        exit 1; 
    fi   
    mkdir -p $WDIR/latex
    for (( i=$case_start; i<=$case_stop; i++ )); do 
        echo CASE $i: ${cases[i]}
        cd $WDIR/${cases[i]}
        
        IFS=" " read -a ARRAY <<< "${outputs[i]}"
        for file in "${ARRAY[@]}"
        do
            cp $file $WDIR/latex
        done
    done
fi





exit 0
