#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "CAM_runcmnd.sh: incorrect number of input arguments"
    exit 1
fi

if [ ! -f ${CAM_SCRIPTDIR}/config_files/$1 ]; then
    echo "CAM_runcmnd.sh: configure options file ${CAM_SCRIPTDIR}/config_files/$1 not found"
    exit 2
fi

hostname=`hostname`
case $hostname in

    ##firefly
    ff* )
    ##search config options file for parallelization info; default on aix is hybrid
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""                                   
	else
            ##open-mp only
#	    cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
	    cmnd="env LSB_PJL_TASK_GEOMETRY="\{\(0\)\}" OMP_NUM_THREADS=${CAM_THREADS} "
	fi
    else
	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    CAM_TASKS=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $2 ) ))
	fi

	num_nodes=`echo $LSB_MCPU_HOSTS | wc -w`
	num_nodes=`expr $num_nodes / 2`
	tpn=`expr $CAM_TASKS / $num_nodes `
	proc=0
	geo_string="\{"
	count1=$num_nodes
	while [ "$count1" != "0" ]; do
	    geo_string="${geo_string}\("
	    count2=$tpn
	    while [ "$count2" != "0" ]; do
		if [ "$count2" != "$tpn" ]; then
		    geo_string="${geo_string}\,"
		fi
		geo_string="${geo_string}$proc"
		proc=`expr $proc + 1`
		count2=`expr $count2 - 1`
	    done
	    geo_string="${geo_string}\)"
	    count1=`expr $count1 - 1`
	done
	geo_string="${geo_string}\}"

	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} mpirun.lsf hybrid_launch "
	else
            ##hybrid
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} OMP_NUM_THREADS=${CAM_THREADS} mpirun.lsf hybrid_launch "
	fi
    fi ;;


    ##yellowstone
    ye* | ys* | ca* )
    ##search config options file for parallelization info; default on aix is hybrid
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""                                   
	else
            ##open-mp only
#	    cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
	    cmnd="env LSB_PJL_TASK_GEOMETRY="\{\(0\)\}" OMP_NUM_THREADS=${CAM_THREADS} "
	fi
    else
	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
            CAM_TASKS=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $2 ) ))
	fi

	num_nodes=`echo $LSB_MCPU_HOSTS | wc -w`
	num_nodes=`expr $num_nodes / 2`
	tpn=`expr $CAM_TASKS / $num_nodes `
	proc=0
	geo_string="\{"
	count1=$num_nodes
	while [ "$count1" != "0" ]; do
	    geo_string="${geo_string}\("
	    count2=$tpn
	    while [ "$count2" != "0" ]; do
		if [ "$count2" != "$tpn" ]; then
		    geo_string="${geo_string}\,"
		fi
		geo_string="${geo_string}$proc"
		proc=`expr $proc + 1`
		count2=`expr $count2 - 1`
	    done
	    geo_string="${geo_string}\)"
	    count1=`expr $count1 - 1`
	done
	geo_string="${geo_string}\}"

	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="env TARGET_CPU_LIST="-1" LSB_PJL_TASK_GEOMETRY=${geo_string} mpirun.lsf "
	else
            ##hybrid
	    cmnd="env TARGET_CPU_LIST="-1" LSB_PJL_TASK_GEOMETRY=${geo_string} OMP_NUM_THREADS=${CAM_THREADS} mpirun.lsf " 
	fi
    fi ;;

    ##goldbach
    go* | g[[:digit:]]* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	elif grep -ic SMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
	else
            ##serial
	    cmnd=""
	fi
    else
	if grep -ic '\-smp' ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
            cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
        else
            ##mpi only
            cmnd=""
            CAM_TASKS=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $2 ) ))
        fi
        cmnd="${cmnd}mpiexec -n ${CAM_TASKS} "
    fi ;;

   ##edison
    ed* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
       if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
           cmnd=""
       elif grep -ic SMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
           cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
       else
            ##serial
           cmnd=""
       fi
    else
       if grep -ic '\-smp' ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
            ppn=24
            tpn=`expr $ppn / $CAM_THREADS `
            cmnd="env OMP_NUM_THREADS=${CAM_THREADS} aprun -N $tpn -n ${CAM_TASKS} -d ${CAM_THREADS} "
        else
            ##mpi only
            CAM_TASKS=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $2 ) ))
            cmnd="aprun -n ${CAM_TASKS} "
        fi
    fi ;;

    ##jaguarpf
    ja* | ti* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	elif grep -ic SMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
	else
            ##serial
	    cmnd=""
	fi
    else
	if grep -ic '\-smp' ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
            ppn=16
            tpn=`expr $ppn / $CAM_THREADS `
            cmnd="env OMP_NUM_THREADS=${CAM_THREADS} aprun -N $tpn -n ${CAM_TASKS} -d ${CAM_THREADS} "
        else
            ##mpi only
            CAM_TASKS=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $2 ) ))
            cmnd="aprun -n ${CAM_TASKS} "
        fi
    fi ;;

    ##hopper
    ho* | ni* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
       if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
           cmnd=""
       elif grep -ic SMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
           cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
       else
            ##serial
           cmnd=""
       fi
    else
       if grep -ic '\-smp' ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
            ppn=24
            tpn=`expr $ppn / $CAM_THREADS `
            cmnd="env OMP_NUM_THREADS=${CAM_THREADS} aprun -N $tpn -n ${CAM_TASKS} -d ${CAM_THREADS} "
        else
            ##mpi only
            CAM_TASKS=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $2 ) ))
            cmnd="aprun -n ${CAM_TASKS} "
        fi
    fi ;;

    ##lynx
    ly* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
       if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
           cmnd=""
       elif grep -ic SMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
           cmnd="env OMP_NUM_THREADS=${CAM_THREADS} "
       else
            ##serial
           cmnd=""
       fi
    else
       if grep -ic '\-smp' ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
            tpn=`expr 12 / $CAM_THREADS `
            cmnd="env OMP_NUM_THREADS=${CAM_THREADS} aprun -N $tpn -n ${CAM_TASKS} -d ${CAM_THREADS} "
        else
            ##mpi only
            CAM_TASKS=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $2 ) ))
            cmnd="aprun -n ${CAM_TASKS} "
        fi
    fi ;;

    * ) 
    echo "CAM_runcmnd.sh: unable to construct run command for unsupported machine $hostname "
    exit 3;;
esac

#store command in temporary file for calling script to access
echo ${cmnd} > cam_run_command.txt
exit 0
