#!/bin/bash 

msls () {

    rd=$1
    ssh_loc=$2
    scp_loc=$3
    if [ "${ssh_loc}" != "" ] && [ "${scp_loc}" != "" ]; then
	ssh -q ${ssh_loc} "ssh -q ${scp_loc} ls -l ${rd}"
    fi
}

msmkdir () {

    rd=$1
    ssh_loc=$2
    scp_loc=$3
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
        echo "msmkdir: ssh -q ${ssh_loc}  ssh -q ${scp_loc} mkdir -p ${rd} ":
        ssh -q ${ssh_loc} "ssh -q ${scp_loc} mkdir -p ${rd}"
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
    # Copy entire directory to the local mass store
    #------------------------------------------------------------------

    ldr=$1
    rdr=$2
    ssh_loc=$3
    scp_loc=$4

    # ssh/scp to ssh_loc by first ssh to ssh_loc. 
    myld=`pwd`

    if  [ `ssh -q ${ssh_loc} "which bbscp" | wc -w` == 1 ]; then
      echo "mscpdir: ssh -q ${ssh_loc}  bbscp -z ${ldr} ${scp_loc}:${rdr}"
      ssh -q ${ssh_loc} "bbscp -z ${ldr} ${scp_loc}:${rdr}"
    else
      echo "mscpdir: ssh -q ${ssh_loc} /usr/bin/scp -r -q ${ldr} ${scp_loc}:${rdr}"
      ssh -q ${ssh_loc} "/usr/bin/scp -r -q ${ldr} ${scp_loc}:${rdr}"
    fi

    sleep 2

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
    scp_loc=$4
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
    if [ "${ssh_loc}" != "" ] && [ "${scp_loc}" != "" ]; then
        ssh -q ${ssh_loc} "scp -q ${ldf} ${scp_loc}:${rdf}"
	sleep 2
    fi
}

#***********************************************************************
# Long term archiving functionality
#***********************************************************************

# Assume that have access to the following environment variables
#   $DOUT_S_ROOT, $DOUT_L_MSROOT, $DOUT_L_HPSS_ACCNT
#   Above name for $MACH is there just for brief backwards compatibility

mode="unknown"
ssh_loc="unknown"
scp_loc="unknown"

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
       --scp_loc )
	   scp_loc=$2
	   shift
	   ;;
       * )
   esac
   shift 
done

found=0
for name in copy_files copy_dirs_hsi copy_dirs_sshscp ; do
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

    # Long-term archiver for HPSS (Trey White, December 6, 2011)
    date

    if [ ! $?DOUT_L_HPSS_ACCNT ]; then
	DOUT_L_HPSS_ACCNT=0
    fi

    # send files to HPSS and delete upon success

    cd $DOUT_S_ROOT
    if [ $DOUT_L_HPSS_ACCNT -gt 0 ]; then
	hsi -a $DOUT_L_HPSS_ACCNT "mkdir -p $DOUT_L_MSROOT ; chmod +t $DOUT_L_MSROOT ; cd $DOUT_L_MSROOT ; put -dPR *"
    else
	hsi "mkdir -p $DOUT_L_MSROOT ; chmod +t $DOUT_L_MSROOT ; cd $DOUT_L_MSROOT ; put -dPR *"
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

if [ "$mode" == "copy_dirs_sshscp" ]; then 

    echo "time : "`date`
    cd $DOUT_S_ROOT

    msmkdir ${DOUT_L_MSROOT} $ssh_loc $scp_loc
    echo "time : "`date`

    for dirl1 in */ ; do 
	cd $DOUT_S_ROOT/${dirl1}
	mscpdir ${DOUT_S_ROOT}/${dirl1} ${DOUT_L_MSROOT} $ssh_loc $scp_loc
	echo "time : "`date`
	for dirl2 in */ ; do
	    cd ${DOUT_S_ROOT}/${dirl1}/${dirl2}
            for file in `ls -1`; do
		if [ -f ${file} ]; then
		    echo "local file: $file ...  long-term archive file: ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file}"
		    lta_listing=`msls ${DOUT_L_MSROOT}/${dirl1}/${dirl2}/${file} $ssh_loc $scp_loc`
		    echo "time : "`date`
		    loc_listing=`ls -l ${file}`
		    lta_size=`msfsize $lta_listing`
		    loc_size=`msfsize $loc_listing`
		    if [ $loc_size -gt 0 ] && [ $loc_size -eq $lta_size ]; then
			echo "local file and long-term archive file are same size"
			echo rm -f ${file}
			rm -f ${file}
		    else
			echo "local file and long-term archive file are NOT the same size... ${file} will remain on local disk"
			#exit -1 #??? ask francis if this is right 
                        # Not sure what to do here... maybe make the log entry and carry on...
		    fi
		fi
	    done # for file
	done # for dirl2
    done # dirl1

fi  # if copy_dirs



