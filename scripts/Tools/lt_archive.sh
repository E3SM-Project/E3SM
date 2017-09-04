#!/bin/bash

msls () {

    rd=$1
    ssh_loc=$2
    if [ "${ssh_loc}" != "" ]; then
	ssh -q ${ssh_loc} "ls -l ${rd}"
    fi
}

msmkdir () {

    rd=$1
    ssh_loc=$2
    if  [ `which hsi | wc -w` == 1 ]; then
	echo "msmkdir: hsi 'mkdir -p ${rd}'"
	hsi -q "mkdir -p ${rd}"
	if [ $? -eq 0 ]; then
	    return 0
	else
	    echo "mksmkdir: error"
	    return $?
	fi
    else
        echo "msmkdir: ssh -q ${ssh_loc} mkdir -p ${rd} ":
        ssh -q ${ssh_loc} "mkdir -p ${rd}"
        sleep 2
	return 0
    fi
}

msfsize() {
    # function to get the size of a file
    # from its long listing
    if [ $# -lt 5 ]; then
        echo "0"
    else
	echo $5
    fi
}

msget () {

    #------------------------------------------------------------------
    # Copy files from the local mass store
    # "Usage msget mssdir/file2 locdir/file1"
    #    rdf = remote dir/filename        ldf = local  dir/filename
    #    rd  = remote dir                 rf  = remote filename
    #    ld  = local  dir                 lf  = local  filename
    # Split inputs into r(remote) and l(local) d(directories) and f(files)
    # If the local filename is empty, set it to the remote filename
    # If the local filename doesn't exist, exit
    #------------------------------------------------------------------
    rdf=$1; rd=`dirname ${rdf}`;  rf=`basename ${rdf}`
    ldf=$2; ld=`dirname ${ldf}`;  lf=`basename ${ldf}`
    if [ "${lf}" == '' ]; then
	lf=${rf}
    fi
    if [ `which hsi | wc -w` == 1 ];  then
	hsi -q "cd ${rd} ; get ${ldf} : ${rf}" >& /dev/null
	return $?
    fi
}

mscpdir () {
    #------------------------------------------------------------------
    # Copy entire directory to the local mass store -- for NAS pleiades lou
    #------------------------------------------------------------------

    ldr=$1
    rdr=$2
    ssh_loc=$3

    # copy dir tree by first ssh to ssh_loc.
    # ssh to lou (lfe) and copy from /nobackup disk system to the lou /ur/$USER
    #  -- /nobackup/ and /ur/ are both mounted on lou
    myld=`pwd`

    if [ `ssh -q ${ssh_loc} "which cxfscp" | wc -w` == 1 ]; then
       echo ssh -q ${ssh_loc} "cxfscp -r ${ldr} ${rdr}"
       ssh -q ${ssh_loc} "cxfscp -r ${ldr} ${rdr}"
    elif [ `ssh -q ${ssh_loc} "which shiftc" | wc -w` == 1 ]; then
      echo "mscpdir: ssh -q ${ssh_loc} shiftc --wait --hosts 8 -r ${ldr} ${rdr}"
      ssh -q ${ssh_loc} "shiftc --wait --hosts 8 -r ${ldr} ${rdr}"
    else
      echo "mscpdir: ssh -q ${ssh_loc} /bin/cp -r -q ${ldr} ${rdr}"
      ssh -q ${ssh_loc} "/bin/cp -r ${ldr} ${rdr}"
    fi

    return 0
}

msput() {

    #------------------------------------------------------------------
    # Copy files to the local mass store
    #    rdf = remote dir/filename    #    ldf = local  dir/filename
    #    rd  = remote dir             #    rf  = remote filename
    #    ld  = local  dir             #    lf  = local  filename
    # Split inputs into r(remote) and l(local) d(directories) and f(files)
    # If the remote file is empty, set it to the local filename
    # Then execute site dependent mass store write
    #------------------------------------------------------------------

    ldf=$1; ld=`dirname ${ldf}`; lf=`basename ${ldf}`
    rdf=$2; rd=`dirname ${rdf}`; rf=`basename ${rdf}`
    ssh_loc=$3

    if [ "${rf}" == "" ]; then
	rf=$lf
    fi
    if [ `which hsi | wc -w` == 1 ]; then
	opts=" "
	if ! [[ "$DOUT_L_HPSS_ACCNT" =~ "0000*" ]]; then
	    opts=" -a ${DOUT_L_HPSS_ACCNT} "
	fi
	# note that the -d flag will delete the local copy
	echo "msput: hsi ${opts} 'cd ${rd} ; put -d ${ldf} : ${rf}'"
	hsi ${opts} -q "cd ${rd} ; put -d ${ldf} : ${rf} ; chmod +r ${rf}"
	return $?
    fi
}

#***********************************************************************
# Long term archiving functionality
#***********************************************************************

# Assume that have access to the following environment variables
#   $DOUT_S_ROOT, $DOUT_L_MSROOT, $DOUT_L_HPSS_ACCNT
#   Above name for $MACH is there just for brief backwards compatibility

mode="copy_dirs_hsi"
ssh_loc="unknown"
arc_root="unknown"
rm_loc_files="TRUE" # default is to remove files from short-term archive

while [ $# -gt 0 ]; do
   case $1 in
       -m|--mode )
	   mode=$2
	   echo " mode is $2"
	   shift
	   ;;
       --ssh_loc )
	   ssh_loc=$2
	   shift
	   ;;
       --arc_root )
	   arc_root=$2
	   shift
	   ;;
       --rm_loc_files )
	   rm_loc_files=$2
	   shift
	   ;;
       * )
   esac
   shift
done

found=0
for name in copy_files copy_dirs_hsi copy_dirs_ssh copy_dirs_local ; do
    if [ "$name" == "$mode" ] ; then
	found=1
	break
    fi
done
if [ $found -ne 1 ] ; then
    echo "$current value of mode $model not supported"
    exit 1
fi

#----------------------------------------------------------------------

if [ "$mode" == "copy_dirs_hsi" ]; then

    if_hsi=`which hsi | wc -w`
    if  [ $if_hsi != 1 ] ; then
	echo "lt_archive: asked for copy_dirs_hsi - but hsi not found"
	echo "lt_archive: check path"
	exit -1
    fi

    date

    # send files to HPSS
    saveFlag="-PR"
    cd $DOUT_S_ROOT
    hsiArgs="mkdir -p $DOUT_L_MSROOT ; chmod +t $DOUT_L_MSROOT ; cd $DOUT_L_MSROOT ; put $saveFlag *"
    # echo $hsiArgs
    if ! [[ "$DOUT_L_HPSS_ACCNT" =~ "0000*" ]]; then
	hsi -a $DOUT_L_HPSS_ACCNT "$hsiArgs"
    else
	hsi "$hsiArgs"
    fi

    date
fi

#----------------------------------------------------------------------

if [ "$mode" == "copy_files" ]; then

    #------------------------------------------------------------------
    # Copy files and dir structure from short term archiving
    # Assume there are up to two levels of dirs below $DOUT_S_ROOT
    #   $DOUT_S_ROOT/$dirl1/$dirl2
    #   dirl1 =>  normallly [atm,lnd,ocn,ice,cpl,glc,rof,wav,rest]
    #   dirl2 =>  normally  [init,hist,logs,date(for rest)]
    #------------------------------------------------------------------
    cd $DOUT_S_ROOT
    msmkdir $DOUT_L_MSROOT
    for dirl1 in */ ; do
	cd ${DOUT_S_ROOT}/${dirl1}
	msmkdir ${DOUT_L_MSROOT}/${dirl1}
	for dirl2 in */ ; do
	    cd ${DOUT_S_ROOT}/${dirl1}/${dirl2}
	    msmkdir ${DOUT_L_MSROOT}/${dirl1}/${dirl2}
	    for file in * ; do
		if [ -f ${file} ]; then
                    # first remove any local file with name checkmssfile
		    if [ -e checkmssfile ]; then
			rm -f checkmssfile
		    fi
                    # try to copy file from mass store into local checkmssfile
		    msget ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file} checkmssfile
                    # compare local file and remote file, either remove local file
                    # OR write local file to mass store based on cmp return status
		    cmp -s ${file} checkmssfile
		    if [ $? == 0 ]; then
			echo "l_archive.sh rm ${file}"
			rm -f $file
		    else
			echo "l_archive.sh: msput ${file} ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file}"
			msput ${file} ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file}
		    fi
		fi
	    done # for file
	done # for dirl2
    done # for dirl

fi # if copy_files

#----------------------------------------------------------------------

if [ "$mode" == "copy_dirs_ssh" ]; then

    echo "START copy_dirs_ssh :"`date`
    cd $DOUT_S_ROOT

    msmkdir ${DOUT_L_MSROOT} $ssh_loc

    for dirl1 in */ ; do
        cd $DOUT_S_ROOT/${dirl1}
        mscpdir ${DOUT_S_ROOT}/${dirl1} ${DOUT_L_MSROOT} $ssh_loc
        echo "time : "`date`
        for dirl2 in */ ; do
            cd ${DOUT_S_ROOT}/${dirl1}/${dirl2}
            for file in `ls -1`; do
                if [ -f ${file} ]; then
                    echo "local file: $file ...  long-term archive file: ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file}"
                    lta_listing=`msls ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file} $ssh_loc`
                    loc_listing=`ls -l ${file}`
                    lta_size=`msfsize $lta_listing`
                    loc_size=`msfsize $loc_listing`
                    if [ $loc_size -gt 0 ] && [ $loc_size -eq $lta_size ]; then
                        echo "local file and long-term archive file are same size"
                        if [ $rm_loc_files == 'TRUE' ]; then
                            echo "rm -f ${file}"
                            rm -f ${file}
                        fi
                    else
                        echo "*** local file and long-term archive file are NOT the same size... ${file} will remain on local disk ***"
                    fi
                fi
            done # for file
        done # for dirl2
    done # dirl1

    echo "DONE copy_dirs_ssh :"`date`
fi  # if copy_dirs_ssh

#----------------------------------------------------------------------

if [ "$mode" == "copy_dirs_local" ]; then

    cd $DOUT_S_ROOT

    mkdir -p ${arc_root}/${DOUT_L_MSROOT}

    for dirl1 in */ ; do
	cd $DOUT_S_ROOT/${dirl1}
	cp -r ${DOUT_S_ROOT}/${dirl1} ${arc_root}/${DOUT_L_MSROOT}
	for dirl2 in */ ; do
	    cd ${DOUT_S_ROOT}/${dirl1}/${dirl2}
            for file in `ls -1`; do
		if [ -f ${file} ]; then
		    echo "local file: $file ...  long-term archive file: ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file}"
		    lta_listing=`ls -l ${arc_root}/${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file}`
		    loc_listing=`ls -l ${file}`
		    lta_size=`msfsize $lta_listing`
		    loc_size=`msfsize $loc_listing`
		    if [ $loc_size -gt 0 ] && [ $loc_size -eq $lta_size ]; then
			echo "local file and long-term archive file are same size"
			echo rm -f ${file}
			rm -f ${file}
		    else
			echo "local file and long-term archive file are NOT the same size... ${file} will remain on local disk"
		    fi
		fi
	    done # for file
	done # for dirl2
    done # dirl1

fi  # if copy_dirs


