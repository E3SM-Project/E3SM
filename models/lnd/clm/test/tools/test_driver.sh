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

    ##bluefire
    be* )
    submit_script="test_driver_bluefire${cur_time}.sh"

    if [ -z "$CLM_ACCOUNT" ]; then
	export CLM_ACCOUNT=`grep -i "^${LOGNAME}:" /etc/project.ncar | cut -f 1 -d "," | cut -f 2 -d ":" `
	if [ -z "${CLM_ACCOUNT}" ]; then
	    echo "ERROR: unable to locate an account number to charge for this job under user: $LOGNAME"
	    exit 2
	fi
    fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

#BSUB -a poe                      # use LSF poe elim
#BSUB -x                          # exclusive use of node (not_shared)
#BSUB -n 192                      # total tasks needed
#BSUB -R "span[ptile=32]"         # max number of tasks (MPI) per node
#BSUB -o test_dr.o%J              # output filename
#BSUB -e test_dr.o%J              # error filename
#BSUB -q regular                  # queue
#BSUB -W 6:00                     
#BSUB -P $CLM_ACCOUNT
#BSUB -J clmtest

if [ -n "\$LSB_JOBID" ]; then   #batch job
    export JOBID=\${LSB_JOBID}
    initdir=\${LS_SUBCWD}
    interactive="NO"
    input_file="tests_pretag_bluefire"
    c_threads=2
    r_threads=4
else
    interactive="YES"
    export LSB_MCPU_HOSTS="\$hostname 8"
    input_file="tests_pretag_bluefire_nompi"
    c_threads=13
    r_threads=25
fi

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=\$c_threads
fi
export CLM_RESTART_THREADS=\$r_threads

##mpi tasks
export CLM_TASKS=96
export CLM_RESTART_TASKS=46

export OBJECT_MODE=64
export XLSMPOPTS="stack=256000000"
export OMP_DYNAMIC=FALSE
export AIXTHREAD_SCOPE=S
export MALLOCMULTIHEAP=TRUE
export MP_LABELIO=yes

# MPI Environment
export MP_RC_USE_LMC=yes
export LAPI_DEBUG_RC_WAIT_ON_QP_SETUP=yes
export MP_INFOLEVEL=2
export MP_EUIDEVICE=sn_all
export MP_SHARED_MEMORY=yes
export LAPI_USE_SHM=yes
export MP_EUILIB=us
# commenting out the following line because we believe it will be better to use 
# the defaults, which change with processor count
#export MP_EAGER_LIMIT=32k
export MP_BULK_MIN_MSG_SIZE=64k
export MP_POLLING_INTERVAL=20000000
export MEMORY_AFFINITY=MCM
export LAPI_DEBUG_ENABLE_AFFINITY=YES
export LAPI_DEBUG_BINDPROC_AFFINITY=YES
export MP_SYNC_QP=YES
export MP_RFIFO_SIZE=16777216
export MP_SHM_ATTACH_THRESH=500000
export MP_EUIDEVELOP=min
export MP_USE_BULK_XFER=yes
export MP_BUFFER_MEM=64M

export MP_RC_MAX_QP=8192
export LAPI_DEBUG_RC_DREG_THRESHOLD=1000000
export LAPI_DEBUG_QP_NOTIFICATION=no
export LAPI_DEBUG_RC_INIT_SETUP=no

. /contrib/Modules/3.2.6/init/sh
module load netcdf/4.1.3_seq
if [ "\$CLM_FC" = "GENIBM" ]; then
  export CESM_MACH="generic_AIX"
else
  export CESM_MACH="bluefire"
fi
export CESM_COMP="ibm"

export NETCDF_DIR=\$NETCDF
export INC_NETCDF=\$NETCDF/include
export LIB_NETCDF=\$NETCDF/lib
export MAKE_CMD="gmake -j "
export CFG_STRING=""
export TOOLS_MAKE_STRING=""
export MACH_WORKSPACE="/glade/scratch"
CPRNC_EXE="/glade/proj3/cseg/tools/cprnc/cprnc"
newcprnc="\$MACH_WORKSPACE/\$LOGIN/newcprnc"
/bin/cp -fp \$CPRNC_EXE \$newcprnc
export CPRNC_EXE="\$newcprnc"
export DATM_QIAN_DATA_DIR="/glade/proj2/cgd/tss/atm_forcing.datm7.Qian.T62.c080727"
dataroot="/glade/proj3/cseg"
export TOOLSLIBS=""
export TOOLS_CONF_STRING=""


echo_arg=""

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##mirage
    mirage* )
    submit_script="test_driver_mirage_${cur_time}.sh"

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
export CLM_TASKS=4
export CLM_RESTART_TASKS=3

export NETCDF_PATH=/fs/local/apps/netcdf-4.1.3_intel
export NETCDF_DIR=\$NETCDF_PATH
export INC_NETCDF=\$NETCDF_PATH/include
export LIB_NETCDF=\$NETCDF_PATH/lib
export intel=/fs/local
export PATH=\${intel}/bin:\${PATH}
export MAKE_CMD="gmake -j5 "
export CESM_MACH="generic_LINUX"
export CESM_COMP="intel"
export CFG_STRING=""
export TOOLS_MAKE_STRING="USER_FC=ifort USER_LINKER=ifort "
export MACH_WORKSPACE="/glade/scratch"
export CPRNC_EXE=/glade/home/erik/bin/cprnc
export DATM_QIAN_DATA_DIR="/glade/proj2/cgd/tss/atm_forcing.datm7.Qian.T62.c080727"
dataroot="/glade/proj3/cseg"
echo_arg="-e"
export LD_LIBRARY_PATH="\$LIB_NETCDF:\${LD_LIBRARY_PATH}"
export SLIBS="-L\$NETCDF_DIR/lib -lnetcdff -lnetcdf"
export TOOLSLIBS=\$SLIBS
export TOOLS_CONF_STRING="-scratchroot \$MACH_WORKSPACE/$USER -max_tasks_per_node 2 -din_loc_root \$dataroot/inputdata"
input_file="tests_posttag_mirage"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ## frankfurt
    frankfurt* ) 
    submit_script="test_driver_frankfurt_${cur_time}.sh"
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
    input_file="tests_posttag_frankfurt"
else
    interactive="YES"
    input_file="tests_posttag_frankfurt_nompi"
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


export CESM_MACH="frankfurt"

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
   export LD_LIBRARY_PATH=\$INTEL/composer_xe_2011_sp1.6.233/compiler/lib/intel64:\$OMPI_INTEL_PATH/lib64:/usr/mpi/intel/openmpi-1.4.3-qlc/lib64
   export LIBRARY_PATH=\$INTEL/composer_xe_2011_sp1.6.233/compiler/lib/intel64:\$OMPI_INTEL_PATH/lib64:/usr/mpi/intel/openmpi-1.4.3-qlc/lib64
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
export MAKE_CMD="gmake -j 5"   ##using hyper-threading on frankfurt
export MACH_WORKSPACE="/scratch/cluster"
export CPRNC_EXE=/fs/cgd/csm/tools/cprnc_64/cprnc
export DATM_QIAN_DATA_DIR="/project/tss/atm_forcing.datm7.Qian.T62.c080727"
dataroot="/fs/cgd/csm"
export TOOLSSLIBS=""
echo_arg="-e"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ## lynx
    lynx* | l0*) 
    submit_script="test_driver_lynx_${cur_time}.sh"

    shell=bash

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/$shell
#

# Name of the queue (CHANGE THIS if needed)
#PBS -q regular
# 
# Number of nodes (CHANGE THIS if needed)
#PBS -l mppwidth=48
#PBS -l walltime=06:00:00
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# Use bourne shell
#PBS -S /bin/$shell
# End of options

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive="NO"
    input_file="tests_posttag_lynx"
else
    interactive="YES"
    input_file="tests_posttag_lynx_nompi"
fi

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=6
fi
export CLM_RESTART_THREADS=3

##mpi tasks
export CLM_TASKS=8
export CLM_RESTART_TASKS=7

#-------------------------------------------------------------------------------
# Runtime environment variables (from Machines_120309)
#-------------------------------------------------------------------------------

# fix for file system problem with empty namelists, june 2011
export DVS_MAXNODES=1

export MPICH_MAX_SHORT_MSG_SIZE=8000 # default is 128000 bytes
export MPICH_PTL_UNEX_EVENTS=960000  # default is  90000 (unexpected recv queue size)
export MPICH_MSGS_PER_PROC=160000    # default is  32768
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
limit coredumpsize unlimited

# The environment variables below produce corefiles and maybe (?) should be
# moved to DEBUG mode at some point
export MPICH_DBMASK=0x200
export MPSTKZ=64M
export OMP_STACKSIZE=64M

#-------------------------------------------------------------------------------
# Modules
#-------------------------------------------------------------------------------

. /opt/modules/default/init/$shell

module rm PrgEnv-intel
module rm PrgEnv-pgi
module rm PrgEnv-cray
module rm PrgEnv-gnu
module rm PrgEnv-pathscale
module rm intel
module rm pgi
module rm cray
module rm pathscale
module rm netcdf
module load subversion


export CESM_MACH="lynx"
if [ "\$CLM_FC" = "PATH" ]; then
   module load PrgEnv-pathscale
   module load netcdf
   module load torque
   module load netcdf
   export CESM_COMP="pathscale"
   export NETCDF_DIR=\$CRAY_NETCDF_DIR/netcdf-pathscale
   export TOOLS_MAKE_STRING="USER_FC=ftn USER_CC=cc "
   export MAKE_CMD="gmake -j 12"   ##using hyper-threading on lynx
elif [ "\$CLM_FC" = "INTEL" ]; then
   module load PrgEnv-intel
   module load netcdf
   module switch intel intel/12.1.0
   module switch xt-mpt    xt-mpt/5.1.4
   module switch xt-libsci xt-libsci/10.5.02
   module load netcdf
   export NETCDF_DIR=\$CRAY_NETCDF_DIR/netcdf-intel
   export CESM_COMP="intel"
   export TOOLS_MAKE_STRING="USER_FC=ftn USER_CC=cc "
   export MAKE_CMD="gmake -j 12"   ##using hyper-threading on lynx
else
   module load PrgEnv-pgi
   module switch pgi       pgi/11.10.0        
   module switch xt-mpt    xt-mpt/4.0.3     
   module switch xt-libsci xt-libsci/10.4.3 
   module load PGI/netcdf4/4.1.3_seq
   export NETCDF_DIR=\$NETCDF
   export CESM_COMP="pgi"
   export TOOLS_MAKE_STRING=" "
   export MAKE_CMD="gmake -j 12"   ##using hyper-threading on lynx
fi

module load esmf/5.2.0rp1
module load subversion

export CFG_STRING=""

export INC_NETCDF=\${NETCDF_DIR}/include
export LIB_NETCDF=\${NETCDF_DIR}/lib
export INC_MPI=""
export LIB_MPI=""
export MACH_WORKSPACE="/glade/scratch"
export CPRNC_EXE=/glade/proj3/cseg/tools/cprnc.lynx/cprnc
export DATM_QIAN_DATA_DIR="/glade/proj2/cgd/tss/atm_forcing.datm7.Qian.T62.c080727"
dataroot="/glade/proj3/cseg"
export TOOLSSLIBS=""
export TOOLS_CONF_STRING="-scratchroot \$MACH_WORKSPACE/$USER -max_tasks_per_node 2 -din_loc_root \$dataroot/inputdata"
echo_arg="-e"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##jaguarpf
    jaguarpf* ) 
    submit_script="test_driver_jaguarpf_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
# #PBS -q batch
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=12:00:00,size=1584
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Use sh
#PBS -S /bin/sh
# Export all Environment variables
#PBS -V
# filesystems to associate with job
#PBS -l gres=widow1%widow2%widow3
#PBS -A CLI017dev
# End of options

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

echo_arg="-e"
if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive="NO"
    input_file="tests_pretag_jaguarpf"
else
    interactive="YES"
    input_file="tests_pretag_jaguarpf_nompi"
    if [ "\$compile_only" = "YES" ]; then
       input_file="tests_pretag_jaguarpf"
    fi
fi


##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=6
fi
export CLM_RESTART_THREADS=12

##mpi tasks
export CLM_TASKS=264
export CLM_RESTART_TASKS=129

source /opt/modules/default/init/sh

module switch pgi       pgi/11.10.0        
module switch xt-mpich2    xt-mpich2/5.4.1
module switch xt-libsci xt-libsci/11.0.04.4
module swap xt-asyncpe xt-asyncpe/5.05
module load szip/2.1
module load hdf5/1.8.7
module load netcdf/4.1.3
module load parallel-netcdf/1.2.0
module load subversion
export CESM_MACH="titan"
export CESM_COMP="pgi"
module load   ncl
module load esmf/5.2.0r_with-lapack+netcdf_O

# The environment variables below produce corefiles and maybe (?) should be
# moved to DEBUG mode at some point
export MPICH_DBMASK=0x200
#limit coredumpsize unlimited

# The environment variable below increase the stack size, which is necessary for
# CICE to run threaded on this machine. 
export MPSTKZ=64M
export OMP_STACKSIZE=64M


export LIB_NETCDF=\${NETCDF_DIR}/lib
export INC_NETCDF=\${NETCDF_DIR}/include
export MOD_NETCDF=\${NETCDF_DIR}/include
export INC_PNETCDF=\${PNETCDF_DIR}/include
export LIB_PNETCDF=\${PNETCDF_DIR}/lib
export CFG_STRING=""
export TOOLS_MAKE_STRING="USER_FC=ftn USER_CC=cc "
export MAKE_CMD="gmake -j 25 "
export MACH_WORKSPACE="/tmp/work"
export CPRNC_EXE=/tmp/proj/ccsm/tools/ccsm_cprnc/cprnc
export DATM_QIAN_DATA_DIR="/tmp/proj/ccsm/inputdata/atm/datm7/atm_forcing.datm7.Qian.T62.c080727"
export ESMFBIN_PATH=$ESMF_BIN
dataroot="/tmp/proj/ccsm"
export TOOLSLIBS=""
export TOOLS_CONF_STRING=""
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
    echo "Only setup to work on: yellowstone, bluefire, mirage, lynx, jaguarpf, frankfurt, and yong"
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
    echo "usage on yellowstone, bluefire, frankfurt, lynx, mirage, titan: "
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
