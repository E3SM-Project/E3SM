#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "CAM_decomp.sh: incorrect number of input arguments"
    exit 1
fi

if [ ! -f ${CAM_SCRIPTDIR}/config_files/$1 ]; then
    echo "CAM_decomp.sh: configure options file ${CAM_SCRIPTDIR}/config_files/$1 not found"
    exit 2
fi

if [ ! -f ${CAM_SCRIPTDIR}/nl_files/$2 ]; then
    echo "CAM_decomp.sh: namelist options file ${CAM_SCRIPTDIR}/nl_files/$2 not found"
    exit 5
fi

##search config options file for parallelization info
if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
    str=""                                   
else
    if grep -ic npr_yz ${CAM_SCRIPTDIR}/nl_files/$2 > /dev/null; then
        #do not override what's already in test files
	str=""
    else

	if grep -ic 'dyn fv' ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then

	    if grep -ic 4x5 ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		lats=46
	    elif grep -ic 2.5x3.33 ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		lats=72
	    elif grep -ic 10x15 ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		lats=19
	    elif grep -ic 1.9x2.5 ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		lats=96
	    elif grep -ic 0.9x1.25 ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		lats=192
	    else
		echo "CAM_decomp.sh: fv resolution not supported "
		exit 3
	    fi

	    if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
                ##mpi only
                ntasks=$(( $CAM_TASKS * $CAM_THREADS / ( $min_cpus_per_task * $3 ) ))
	    else
                ##hybrid
		ntasks=$CAM_TASKS
	    fi

	    lats_per_task=`expr $lats / $ntasks`

	    if [ $lats_per_task -lt 3 ]; then

		if grep -ic 'chem waccm' ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		    levs=66
		elif grep -ic cam4 ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		    levs=26
		elif grep -ic cam3 ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
		    levs=26
		else
		    levs=30
		fi

		z=0
		done="FALSE"
		while [ "$done" = "FALSE" ]; do
		    z=`expr $z + 1`
		    y=`expr $ntasks / $z`
		    prod=`expr $y "*" $z`
		    lats_per_task=`expr $lats / $y`
		    if [ $prod -eq $ntasks ] && [ $lats_per_task -ge 3 ]; then
			done="TRUE"
		    fi
		    if [ $z -gt $levs ]; then
			echo "CAM_decomp.sh: unable to construct fv decomposition string "
			exit 3
		    fi
		done
		str="npr_yz=$y,$z,$z,$y"

	    else
                #fv 2d decomp not required
		str=""
	    fi
	else
	    str=""
	fi
    fi
fi 

#store command in temporary file for calling script to access
echo ${str} > cam_decomp_string.txt
exit 0
