#!/bin/sh 
#
# test_driver.sh:  driver script for the offline testing of CLM of tools
#
# interactive usage on all machines:
#
# env CLM_SOFF=FALSE ./test_driver.sh -i
#
# valid arguments: 
# -i    interactive usage
# -d    debug usage -- display tests that will run -- but do NOT actually execute them
# -f    force batch submission (avoids user prompt)
# -h    displays this help message
#
#
# **pass environment variables by preceding above commands 
#   with 'env var1=setting var2=setting '
# **more details in the CLM testing user's guide, accessible 
#   from the CLM developers web page


#will attach timestamp onto end of script name to prevent overwriting
cur_time=`date '+%H:%M:%S'`

hostname=`hostname`
case $hostname in

    ##yellowstone
    ys* )
    submit_script="test_driver_yellowstone${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

interactive="YES"
input_file="tests_pretag_yellowstone_nompi"
c_threads=16

source /glade/apps/opt/lmod/lmod/init/sh

module purge
module load ncarenv/1.0
module load ncarbinlibs/1.0

module load intel/12.1.5
module load ncarcompilers/1.0
module load netcdf/4.2 
module load pnetcdf/1.3.0

module load nco
module load ncl


##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=\$c_threads
fi

export CESM_MACH="yellowstone"
export CESM_COMP="intel"

export NETCDF_DIR=\$NETCDF
export INC_NETCDF=\$NETCDF/include
export LIB_NETCDF=\$NETCDF/lib
export MAKE_CMD="gmake -j "
export CFG_STRING=""
export TOOLS_MAKE_STRING="USER_FC=ifort USER_LINKER=ifort "
export MACH_WORKSPACE="/glade/scratch"
CPRNC_EXE="$CESMDATAROOT/tools/cprnc/cprnc"
dataroot="$CESMDATAROOT"
export TOOLSLIBS=""
export TOOLS_CONF_STRING=""


echo_arg=""

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ## goldbach
    goldbach* ) 
    submit_script="test_driver_goldbach_${cur_time}.sh"
    export PATH=/cluster/torque/bin:${PATH}

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
#PBS -q long
# Number of nodes (CHANGE THIS if needed)
#PBS -l nodes=1:ppn=8
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive="NO"
    input_file="tests_posttag_goldbach"
else
    interactive="YES"
    input_file="tests_posttag_goldbach_nompi"
fi

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=2
fi
export CLM_RESTART_THREADS=1

##mpi tasks
export CLM_TASKS=8
export CLM_RESTART_TASKS=7

export PGI=/usr/local/pgi
export INTEL=/usr/local/intel-cluster
export P4_GLOBMEMSIZE=500000000


export CESM_MACH="goldbach"

if [ "\$CLM_FC" = "PGI" ]; then
   export NETCDF_PATH=/usr/local/netcdf-pgi
   export LD_LIBRARY_PATH=\${PGI}/linux86/lib:/cluster/torque/lib:\${LD_LIBRARY_PATH}
   export LIBRARY_PATH=\${PGI}/linux86/lib:/cluster/torque/lib:\${LIBRARY_PATH}
   export PATH=\${PGI}/linux86/bin:\${PATH}
   export CESM_COMP="pgi"
   export TOOLS_MAKE_STRING=""
   export TOOLS_CONF_STRING=""
   export CFG_STRING=""
else
   export OMPI_INTEL_PATH=/cluster/openmpi-qlc-intel
    # Runtime environment variables
   export PATH=\${OMPI_INTEL_PATH}/bin:\${INTEL}/bin:\${PATH}

   export NETCDF_PATH=/usr/local/netcdf-intel-cluster
   export LD_LIBRARY_PATH=\$INTEL/composer_xe_2011_sp1.6.233/compiler/lib/intel64:\$OMPI_INTEL_PATH/lib64:/usr/mpi/intel/openmpi-1.4.3-qlc/lib64:/usr/local/intel-cluster-2013.4.183/compiler/lib/intel64
   export LIBRARY_PATH=\$INTEL/composer_xe_2011_sp1.6.233/compiler/lib/intel64:\$OMPI_INTEL_PATH/lib64:/usr/mpi/intel/openmpi-1.4.3-qlc/lib64:/usr/local/intel-cluster-2013.4.183/compiler/lib/intel64
   export CESM_COMP="intel"
   export TOOLS_MAKE_STRING="USER_FC=ifort USER_CC=icc "
   export TOOLS_CONF_STRING=""
   export CFG_STRING=""
fi

export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${NETCDF_PATH}/lib
export LIBRARY_PATH=\${LIBRARY_PATH}:\${NETCDF_PATH}/lib
export NETCDF_DIR=\$NETCDF_PATH
export INC_NETCDF=\${NETCDF_PATH}/include
export LIB_NETCDF=\${NETCDF_PATH}/lib
export MAKE_CMD="gmake -j 5"   ##using hyper-threading on goldbach
export MACH_WORKSPACE="/scratch/cluster"
export CPRNC_EXE=/fs/cgd/csm/tools/cprnc_64/cprnc
export DATM_QIAN_DATA_DIR="/project/tss/atm_forcing.datm7.Qian.T62.c080727"
dataroot="/fs/cgd/csm"
export TOOLSSLIBS=""
echo_arg="-e"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##yong
    yong* | vpn* )
    submit_script="test_driver_yong_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

interactive="YES"

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=2
fi
export CLM_RESTART_THREADS=1

##mpi tasks
export CLM_TASKS=2
export CLM_RESTART_TASKS=1

export CESM_MACH="generic_Darwin"
if [ "\$CLM_FC" = "PGI" ]; then
   export CESM_COMP="pgi"
   export NETCDF_PATH=/usr/local/netcdf-3.6.3-pgi-10.9
   export CFG_STRING=""
   export TOOLS_MAKE_STRING=""
else
   export CESM_COMP="intel"
   export NETCDF_PATH=/usr/local/netcdf-3.6.3-intel-11.1
   export MPICH_PATH=/usr/local/mpich2-1.3.1-intel-11.1
   export PATH="\$MPICH_PATH/bin:$PATH"
   export CFG_STRING=""
   export TOOLS_MAKE_STRING="USER_FC=ifort USER_LINKER=ifort USER_CC=icc "
   export DYLD_LIBRARY_PATH=/opt/intel/Compiler/11.1/067/lib
fi
export NETCDF_DIR=\$NETCDF_PATH
export INC_NETCDF=\$NETCDF_PATH/include
export LIB_NETCDF=\$NETCDF_PATH/lib
export MAKE_CMD="make -j 4"
export MACH_WORKSPACE="/glade/scratch"
export CPRNC_EXE=$HOME/bin/newcprnc
export DATM_QIAN_DATA_DIR="/fis/cgd/cseg/csm/inputdata/atm/datm7/atm_forcing.datm7.Qian.T62.c080727"
export ESMFBIN_PATH=\
"/usr/local/esmf_5_2_0/DEFAULTINSTALLDIR/bin/binO/Darwin.intel.64.mpiuni.default"
dataroot="/fis/cgd/cseg/csm"
echo_arg=""
input_file="tests_posttag_yong"
export TOOLSLIBS=""
export TOOLS_CONF_STRING="-scratchroot \$MACH_WORKSPACE/$USER -max_tasks_per_node 2 -din_loc_root \$dataroot/inputdata"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;
    * )
    echo "Only setup to work on: yellowstone, goldbach, and yong"
    exit
 

esac

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat >> ./${submit_script} << EOF

if [ -n "\${CLM_JOBID}" ]; then
    export JOBID=\${CLM_JOBID}
fi
##check if interactive job

if [ "\$interactive" = "YES" ]; then

    if [ -z "\${JOBID}" ]; then
       export JOBID=\$\$
    fi
    echo "test_driver.sh: interactive run - setting JOBID to \$JOBID"
    if [ \$0 = "test_driver.sh" ]; then
	initdir="."
    else
	initdir=\${0%/*}
    fi
else
    echo "ERROR: you *always* need to use the interactive option (-i)"
    echo "       currently doesn't work without it"
    exit 3
fi

##establish script dir and clm_root
if [ -f \${initdir}/test_driver.sh ]; then
    export CLM_SCRIPTDIR=\`cd \${initdir}; pwd \`
    export CLM_ROOT=\`cd \${CLM_SCRIPTDIR}/../../../../..; pwd \`
else
    if [ -n "\${CLM_ROOT}" ] && [ -f \${CLM_ROOT}/models/lnd/clm*/test/tools/test_driver.sh ]; then
	export CLM_SCRIPTDIR=\`cd \${CLM_ROOT}/models/lnd/clm*/test/tools; pwd \`
    else
	echo "ERROR: unable to determine script directory "
	echo "       if initiating batch job from directory other than the one containing test_driver.sh, "
	echo "       you must set the environment variable CLM_ROOT to the full path of directory containing "
        echo "       <models>. "
	exit 3
    fi
fi

##output files
clm_log=\${initdir}/td.\${JOBID}.log
if [ -f \$clm_log ]; then
    rm \$clm_log
fi
clm_status=\${initdir}/td.\${JOBID}.status
if [ -f \$clm_status ]; then
    rm \$clm_status
fi

##setup test work directory
if [ -z "\$CLM_TESTDIR" ]; then
    export CLM_TESTDIR=\${MACH_WORKSPACE}/\$LOGNAME/clmTests/test-driver.\${JOBID}
    if [ -d \$CLM_TESTDIR ] && [ \$CLM_RETAIN_FILES != "TRUE" ]; then
        rm -r \$CLM_TESTDIR
    fi
fi
if [ ! -d \$CLM_TESTDIR ]; then
    mkdir -p \$CLM_TESTDIR
    if [ \$? -ne 0 ]; then
	echo "ERROR: unable to create work directory \$CLM_TESTDIR"
	exit 4
    fi
fi

## MCT and PIO build directorys
export MCT_LIBDIR=\$CLM_TESTDIR/mct
export PIO_LIBDIR=\$CLM_TESTDIR/pio

##set our own environment vars
export CSMDATA=\${dataroot}/inputdata
export DIN_LOC_ROOT=\${CSMDATA}
export MPI_TYPE_MAX=100000

##process other env vars possibly coming in
if [ -z "\$CLM_RETAIN_FILES" ]; then
    export CLM_RETAIN_FILES=FALSE
fi
if [ -n "\${CLM_INPUT_TESTS}" ]; then
    input_file=\$CLM_INPUT_TESTS
else
    input_file=\${CLM_SCRIPTDIR}/\${input_file}
fi
if [ ! -f \${input_file} ]; then
    echo "ERROR: unable to locate input file \${input_file}"
    exit 5
fi

if [ \$interactive = "YES" ]; then
    echo "reading tests from \${input_file}"
else
    echo "reading tests from \${input_file}" >> \${clm_log}
fi

num_tests=\`wc -w < \${input_file}\`
echo "STATUS OF CLM TESTING UNDER JOB \${JOBID};  scheduled to run \$num_tests tests from:" >> \${clm_status}
echo "\$input_file" >> \${clm_status}
echo "" >> \${clm_status}
echo " on machine: $hostname" >> \${clm_status}
if [ -n "${BL_ROOT}" ]; then
   echo "tests of baseline will use source code from:" >> \${clm_status}
   echo "\$BL_ROOT" >> \${clm_status}
fi
if [ \$interactive = "NO" ]; then
    echo "see \${clm_log} for more detailed output" >> \${clm_status}
fi
echo "" >> \${clm_status}

test_list=""
while read input_line; do
    test_list="\${test_list}\${input_line} "
done < \${input_file}


##initialize flags, counter
skipped_tests="NO"
pending_tests="NO"
count=0

##loop through the tests of input file
for test_id in \${test_list}; do
    count=\`expr \$count + 1\`
    while [ \${#count} -lt 3 ]; do
        count="0\${count}"
    done

    master_line=\`grep \$test_id \${CLM_SCRIPTDIR}/input_tests_master\`
    status_out=""
    for arg in \${master_line}; do
        status_out="\${status_out}\${arg} "
    done

    if [ -z "\$status_out" ]; then
	echo "No test matches \$test_id in \${CLM_SCRIPTDIR}/input_tests_master"
        exit 3
    fi

    test_cmd=\${status_out#* }

    status_out="\${count} \${status_out}"

    if [ \$interactive = "YES" ]; then
        echo ""
        echo "***********************************************************************************"
        echo "\${status_out}"
        echo "***********************************************************************************"
    else
        echo "" >> \${clm_log}
        echo "***********************************************************************************"\
            >> \${clm_log}
        echo "\$status_out" >> \${clm_log}
        echo "***********************************************************************************"\
            >> \${clm_log}
    fi

    if [ \${#status_out} -gt 94 ]; then
        status_out=\`echo "\${status_out}" | cut -c1-100\`
    fi
    while [ \${#status_out} -lt 97 ]; do
        status_out="\${status_out}."
    done

    echo \$echo_arg "\$status_out\c" >> \${clm_status}

    if [   \$interactive = "YES" ]; then
        \${CLM_SCRIPTDIR}/\${test_cmd}
        rc=\$?
    else
        \${CLM_SCRIPTDIR}/\${test_cmd} >> \${clm_log} 2>&1
        rc=\$?
    fi
    if [ \$rc -eq 0 ]; then
        echo "PASS" >> \${clm_status}
    elif [ \$rc -eq 255 ]; then
        echo "SKIPPED*" >> \${clm_status}
        skipped_tests="YES"
    elif [ \$rc -eq 254 ]; then
        echo "PENDING**" >> \${clm_status}
        pending_tests="YES"
    else
        echo " rc=\$rc FAIL" >> \${clm_status}
	if [ \$interactive = "YES" ]; then
	    if [ "\$CLM_SOFF" != "FALSE" ]; then
		echo "stopping on first failure"
		echo "stopping on first failure" >> \${clm_status}
		exit 6
	    fi
	else
	    if [ "\$CLM_SOFF" = "TRUE" ]; then
		echo "stopping on first failure" >> \${clm_status}
		echo "stopping on first failure" >> \${clm_log}
		exit 6
	    fi
	fi
    fi
done

echo "end of input" >> \${clm_status}
if [ \$interactive = "YES" ]; then
    echo "end of input"
else
    echo "end of input" >> \${clm_log}
fi

if [ \$skipped_tests = "YES" ]; then
    echo "*  please verify that any skipped tests are not required of your clm commit" >> \${clm_status}
fi
if [ \$pending_tests = "YES" ]; then
    echo "** tests that are pending must be checked manually for a successful completion" >> \${clm_status}
    if [ \$interactive = "NO" ]; then
	echo "   see the test's output in \${clm_log} " >> \${clm_status}
	echo "   for the location of test results" >> \${clm_status}
    fi
fi

if [ "\$interactive" = "YES" ]; then
   passInt="test_driver.sh-i"
else
   passInt="test_driver.sh"
fi

../../bld/unit_testers/xFail/wrapClmTests.pl -statusFile "\${clm_status}" -numberOfTests "\${num_tests}" -callingScript "\${passInt}"

exit 0

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^


chmod a+x $submit_script
if [ ! -z "$CLM_RETAIN_FILES" ]; then
   export CLM_RETAIN_FILES="FALSE"
fi
arg1=${1##*-}
case $arg1 in
    [iI]* )
    debug="NO"
    interactive="YES"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
    ./${submit_script}
    exit 0
    ;;

    [cC]* )
    debug="NO"
    interactive="YES"
    compile_only="YES"
    export debug
    export CLM_RETAIN_FILES="TRUE"
    export interactive
    export compile_only
    export CLM_RETAIN_FILES="TRUE"
    ./${submit_script}
    exit 0
    ;;

    [dD]* )
    debug="YES"
    interactive="YES"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
    ./${submit_script}
    exit 0
    ;;

    [fF]* )
    debug="NO"
    interactive="NO"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
    ;;

    "" )
    echo ""
    echo "**********************"
    echo "$submit_script has been created and will be submitted to the batch queue..."
    echo "(ret) to continue, (a) to abort"
    read ans
    case $ans in
	[aA]* ) 
	echo "aborting...type ./test_driver.sh -h for help message"
	exit 0
	;;
    esac
    debug="NO"
    interactive="NO"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
    ;;

    * )
    echo ""
    echo "**********************"
    echo "usage on yellowstone, goldbach, and yongi: "
    echo "./test_driver.sh -i"
    echo ""
    echo "valid arguments: "
    echo "-i    interactive usage"
    echo "-c    compile-only usage (run configure and compile do not run clm)"
    echo "-d    debug-only  usage (run configure and build-namelist do NOT compile or run clm)"
    echo "-f    force batch submission (avoids user prompt)"
    echo "-h    displays this help message"
    echo ""
    echo "**pass environment variables by preceding above commands "
    echo "  with 'env var1=setting var2=setting '"
    echo ""
    echo "**********************"
    exit 0
    ;;
esac

echo "submitting..."
case $hostname in
    #default
    * )
    echo "no submission capability on this machine use the interactive option: -i"
    exit 0
    ;;

esac
exit 0
