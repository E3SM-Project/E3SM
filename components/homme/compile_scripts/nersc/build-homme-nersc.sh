#!/bin/bash
# Compiles standalone HOMME for various (but not all) test cases.
# Must be adjusted for each user, and depends on machine.
# Currently set up for account jlturner on Perlmutter.
# Author: Jason Torchinsky

# Parse command-line input
while getopts t: flag
do
	case "${flag}" in
		t) test=${OPTARG};;
	esac
done

# Declare development branch for E3SM
# Useful for developing multiple branches simultaneously
branch=jasonltorchinsky/scripts/current

# Declare directory paths
# Current directory to return to
curr=$(pwd)
# E3SM source code directory
e3sm=$(pwd)/../../../..
# HOMME subdirectory
homme=$e3sm/components/homme
# Scratch space to contain E3SM compiled files
scratch=$PSCRATCH/E3SM
# Working directory to compile the development branch of HOMME to
wdir=$scratch/$branch/homme
# Directory to machine file, must be changed for each machine
mach=$homme/cmake/machineFiles/perlmutter-nocuda-gnu.cmake

# Load the necessary modules
echo "-- WARNING: Ensure that necessary modules are lodaded by running"
echo "            'source get_case_env.sh'"

# Selects which test to compile
case "${test}" in
	dcmip2012_1_1) 
		testname='DCMIP 2012 Test 1.1 - 3D Deformational Flow'
		testdir=$wdir/dcmip_tests/dcmip2012_test1.1_3d_deformational_flow/preqx
		mode='preqx'
		;;
	dcmip2012_1_2) 
		testname='DCMIP 2012 Test 1.2 - Hadley-Like Meridional Circulation'
		testdir=$wdir/dcmip_tests/dcmip2012_test1.2_hadley_meridional_circulation/preqx
		mode='preqx'
		;;
	dcmip2012_1_3) 
		testname='DCMIP 2012 Test 1.3 - Thin Clouds Over Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test1.3_thin_clouds_over_orography/preqx
		mode='preqx'
		;;
	dcmip2012_2_0) 
		testname='DCMIP 2012 Test 2.0 - Steady State with Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test2.0_steady_state_with_orography/theta-l
		mode='theta-l'
		;;
	dcmip2012_2_1) 
		testname='DCMIP 2012 Test 2.1 - Non-Sheared Background Flow with Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test2.1_nh_mountain_waves_no_shear/theta-l
		mode='theta-l'
		;;
	dcmip2012_2_2) 
		testname='DCMIP 2012 Test 2.2 - Sheared Background Flow with Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test2.2_nh_mountain_waves_with_shear/theta-l
		mode='theta-l'
		;;
	dcmip2012_3) 
		testname='DCMIP 2012 Test 3 - Non-Hydrostatic Gravity Waves'
		testdir=$wdir/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta-l
		mode='theta-l'
		;;
	dcmip2012_4_1) 
		testname='DCMIP 2012 Test 4.1 - Baroclinic Instability'
		testdir=$wdir/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/theta-l
		mode='theta-l'
		;;
	dcmip2016_1) 
		testname='DCMIP 2016 Test 1 - Moist Baroclinic Wave'
		testdir=$wdir/dcmip_tests/dcmip2016_test1_baroclinic_wave/theta-l
		mode='theta-l'
		;;
	dcmip2016_2) 
		testname='DCMIP 2016 Test 2 - Tropical Cyclone'
		testdir=$wdir/dcmip_tests/dcmip2016_test2_tropical_cyclone/theta-l
		mode='theta-l'
		;;
	dcmip2016_3) 
		testname='DCMIP 2016 Test 3 - Supercell'
		testdir=$wdir/dcmip_tests/dcmip2016_test3_supercell/theta-l
		mode='theta-l'
		;;
	*)
		echo "-- Unable to parse test number, or test number is unsupported. Aborting..."
                exit 2
		;;
esac

# Ensure directory paths exists
mkdir -p $wdir
mkdir -p $testdir

echo "-- Setting up test: ${testname}..."

export e3sm
export branch
export homme
export scratch
export wdir
export mach
export testdir
export mode

# Compile HOMME
echo "-- Compiling HOMME..."

cd $wdir
cmake -C $mach $homme
make -j4 $mode

# Compile Test
cd $testdir
make install
./build.sh

# Return to directory this script was called from
cd $curr
