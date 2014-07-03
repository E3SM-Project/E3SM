#!/bin/sh

#short-term archive script - move model output out of run directory
#to free disc space for next invocation of model run
#must be executed from run directory

#function dispose:
#moves output files to specified area of st archive and will
#process interim files along the way: 
#arg1 => interim files flag
#arg2 => destination
#remaining args => actual files to be processed 
dispose() {
    if [ "$1" == "ifiles_y" ] && [ "$DOUT_S_SAVE_INT_REST_FILES" != "TRUE" ]; then
	shift
	shift
	rm $*       2> /dev/null
    else
	shift
	dest=$1
	mkdir -p $dest
	shift
	mv $* $dest 2> /dev/null
    fi
}

#function get_inst_suffix:
#Gets a string corresponding to the current instance index of a given
#component, with a leading '.'; this can be appended to a file name.
#In the general case, the returned string is something like ".###",
#but if there is only a single instance of the given component, then
#the returned string is empty, because in that case the file names
#don't have any instance number (or the associated extra '_').
#
#arg1 => instance index for a given component
#arg2 => number of instances of this component
#
#result is returned in $inst_suffix
get_inst_suffix() {
#    echo "get_inst_suffixa $1 $2 ${inst_suffix}"
    if [ $2 -eq 1 ]; then   # only one instance of this component
	inst_suffix=""
    else                    # multiple instances of this component
        inst_suffix=`printf _%04d $1`
    fi
#    echo "get_inst_suffixb $1 $2 ${inst_suffix}"
}

echo ""
echo "st_archive.sh: start of short-term archiving"

#validate required env var settings
if [ -z "$DOUT_S_ROOT" ]; then
    echo "st_archive.sh: error, environment variable DOUT_S_ROOT is required "
    echo "               for root location of short-term archive"
    echo "st_archive.sh: exiting"
    exit 1
fi

sta=${DOUT_S_ROOT}/.sta-$$-`date +%Y%m%d%H%M%S%N`
mkdir -p ${sta} 2> /dev/null

if [ $? -ne 0 ]; then
    echo "st_archive.sh: error, unable to create short-term archive directory"
    echo "st_archive.sh: exiting"
    exit 1
fi
mv ${DOUT_S_ROOT}/* ${sta}

if [ -z "$DOUT_S_SAVE_INT_REST_FILES" ]; then
    echo "st_archive.sh: warning, environment variable DOUT_S_SAVE_INT_REST_FILES is not "
    echo "               set - using "FALSE" as default for saving interim restart files"
    export DOUT_S_SAVE_INT_REST_FILES=FALSE
fi

if [ "$DOUT_S_SAVE_INT_REST_FILES" == "FALSE" ]; then
    echo "st_archive.sh: restart files from end of run will be saved, "
    echo "               interim restart files will be deleted"
fi

## Check component instance counts

if [ -z "$NINST_ATM" ]; then
    echo "st_archive.sh: warning, NINST_ATM not set -- using 1 instance"
    export NINST_ATM=1
fi

if [ -z "$NINST_LND" ]; then
    echo "st_archive.sh: warning, NINST_LND not set -- using 1 instance"
    export NINST_LND=1
fi

if [ -z "$NINST_ICE" ]; then
    echo "st_archive.sh: warning, NINST_ICE not set -- using 1 instance"
    export NINST_ICE=1
fi

if [ -z "$NINST_OCN" ]; then
    echo "st_archive.sh: warning, NINST_OCN not set -- using 1 instance"
    export NINST_OCN=1
fi

if [ -z "$NINST_GLC" ]; then
    echo "st_archive.sh: warning, NINST_GLC not set -- using 1 instance"
    export NINST_GLC=1
fi

if [ -z "$NINST_WAV" ]; then
    echo "st_archive.sh: warning, NINST_WAV not set -- using 1 instance"
    export NINST_WAV=1
fi

#create directory for restart files
set ${CASE}.cpl.r.*
cplfile=`ls -rt $* 2> /dev/null | tail -1`
dname=`echo $cplfile | sed "s/\.nc//; s/^.*\.r\.//;"`
if [ -d ${sta}/rest/${dname} ]; then
    rm -rf ${sta}/rest/${dname}
fi
mkdir -p ${sta}/rest/${dname}
if [ $? -ne 0 ]; then
    echo "st_archive.sh: error, unable to create rest directory"
    echo "st_archive.sh: exiting"
    exit 1
fi

#populate temp directory with pointer files
set rpointer.*
if [ $# -le 0 ]; then
    echo "st_archive.sh: error, script should be invoked from run directory..."
    echo "               expecting restart pointer files of the form 'rpointer.<component>'"
    echo "               but did not find any: $*"
    echo "st_archive.sh: exiting"
    exit 1
fi
mv $* ${sta}/rest/${dname}

set cpl.log.*;                                                                                                        dispose ifiles_n ${sta}/cpl/logs $*
set cesm*.log.*;                                                                                                      dispose ifiles_n ${sta}/cpl/logs $*
set ${CASE}.cpl.r.*;         latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/cpl/rest $*
set ${CASE}.cpl.h* ;                                                                                                  dispose ifiles_n ${sta}/cpl/hist $*


# DART assimilation-related files
set dart_log.*;                                                                                                       dispose ifiles_n ${sta}/dart/logs $*
set True_State.*.nc;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*
set Prior_Diag.*.nc;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*
set Posterior_Diag.*.nc;                                                                                              dispose ifiles_n ${sta}/dart/hist $*
set obs_seq.*.out;                                                                                                    dispose ifiles_n ${sta}/dart/hist $*
set obs_seq.*.final;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*
set pr*inflate_restart*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null;    dispose ifiles_n ${sta}/dart/rest $*
set po*inflate_restart*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null;    dispose ifiles_n ${sta}/dart/rest $*


IDX=1
while [ $IDX -le $NINST_ATM ]
do
    get_inst_suffix $IDX $NINST_ATM
    set atm${inst_suffix}.log.*;                                                                                                   dispose ifiles_n ${sta}/atm/logs $*
    set ${CASE}.cam*${inst_suffix}.r.*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.rs.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.ra.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.rh0.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.rh1.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.rh2.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.rh3.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.rh4.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.rh5.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam*${inst_suffix}.h0.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam*${inst_suffix}.h1.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam*${inst_suffix}.h2.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam*${inst_suffix}.h3.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam*${inst_suffix}.h4.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam*${inst_suffix}.h5.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam*${inst_suffix}.hs.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam*${inst_suffix}.i.*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/init $*

    set ${CASE}.datm${inst_suffix}.r.* ;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.datm${inst_suffix}.rs* ;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.datm${inst_suffix}.h.* ;                                                                                           dispose ifiles_n ${sta}/atm/hist $*

    set ${CASE}.wrf.r01.*;                latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.wrf.r02.*;                latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.wrf.r03.*;                latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.wrf.h01.*;                                                                                                         dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h02.*;                                                                                                         dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h03.*;                                                                                                         dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h1aux01.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h1aux02.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h1aux03.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h2aux01.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h2aux02.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h2aux03.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h3aux01.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h3aux02.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.wrf.h3aux03.*;                                                                                                     dispose ifiles_n ${sta}/atm/hist $*
 
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_LND ]
do
    get_inst_suffix $IDX $NINST_LND
    set lnd${inst_suffix}.log.*;                                                                                                   dispose ifiles_n ${sta}/lnd/logs $*
    set ${CASE}.clm?${inst_suffix}.r.*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh0.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh1.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh2.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh3.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh4.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh5.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.h0.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h1.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h2.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h3.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h4.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h5.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.hv.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.i.*;                                                                                            dispose ifiles_y ${sta}/lnd/init $*

    set ${CASE}.dlnd${inst_suffix}.r.* ;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.dlnd${inst_suffix}.rs* ;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.dlnd${inst_suffix}.h.* ;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_ROF ]
do
    get_inst_suffix $IDX $NINST_ROF
    set rof${inst_suffix}.log.*;                                                                                                  dispose ifiles_n ${sta}/rof/logs $*
    set ${CASE}.rtm${inst_suffix}.r.*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/rof/rest $*
    set ${CASE}.rtm${inst_suffix}.rh0.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/rof/rest $*
    set ${CASE}.rtm${inst_suffix}.rh1.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/rof/rest $*
    set ${CASE}.rtm${inst_suffix}.rh2.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/rof/rest $*
    set ${CASE}.rtm${inst_suffix}.rh3.*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/rof/rest $*
    set ${CASE}.rtm${inst_suffix}.h0.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/rof/hist $*
    set ${CASE}.rtm${inst_suffix}.h1.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/rof/hist $*
    set ${CASE}.rtm${inst_suffix}.h2.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/rof/hist $*
    set ${CASE}.rtm${inst_suffix}.h3.*;  latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/rof/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_ICE ]
do
    get_inst_suffix $IDX $NINST_ICE
    set ice${inst_suffix}.log.*;                                                                                                      dispose ifiles_n ${sta}/ice/logs $*
    set ${CASE}.cice${inst_suffix}.r.[0-9]*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.volpn*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.dEdd*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.age*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.aero*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.h*;                                                                                                dispose ifiles_n ${sta}/ice/hist $*
    set ${CASE}.cice${inst_suffix}.i.*;                                                                                               dispose ifiles_y ${sta}/ice/init $*

    set ${CASE}.dice${inst_suffix}.r.* ;     latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.dice${inst_suffix}.rs* ;     latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.dice${inst_suffix}.h.* ;                                                                                              dispose ifiles_n ${sta}/ice/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_OCN ]
do
    get_inst_suffix $IDX $NINST_OCN
    set ocn${inst_suffix}.log.*;                                                                                                               dispose ifiles_n ${sta}/ocn/logs $*
    set ${CASE}.pop${inst_suffix}.r.*.hdr;            latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.r.*0000;            latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.r.*0000.nc;         latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.ecosys.*.hdr;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.ecosys.*0000;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.ecosys.*0000.nc; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.*.hdr;           latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.*0000;           latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.*0000.nc;        latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.ro.*;               latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.d?*;                                                                                                         dispose ifiles_n ${sta}/ocn/hist $*
    set ${CASE}.pop${inst_suffix}.h*;                                                                                                          dispose ifiles_n ${sta}/ocn/hist $*

    set ${CASE}.docn${inst_suffix}.r.* ;              latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.docn${inst_suffix}.rs* ;              latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.docn${inst_suffix}.h.* ;                                                                                                       dispose ifiles_n ${sta}/ocn/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_GLC ]
do
    get_inst_suffix $IDX $NINST_GLC
    set glc${inst_suffix}.log.*;                                                                                                      dispose ifiles_n ${sta}/glc/logs $*
    set ${CASE}.cism${inst_suffix}.r.[0-9]*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.volpn*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.dEdd*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.age*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.aero*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.h*;                                                                                                dispose ifiles_n ${sta}/glc/hist $*
    set ${CASE}.cism${inst_suffix}.i.*;                                                                                               dispose ifiles_y ${sta}/glc/init $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_WAV ]
do
    get_inst_suffix $IDX $NINST_WAV
    set wav${inst_suffix}.log.*;                                                                                                     dispose ifiles_n ${sta}/wav/logs $*
    set ${CASE}.ww3${inst_suffix}.r.[0-9]*; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/wav/rest $*
    set ${CASE}.ww3${inst_suffix}.h*;                                                                                                dispose ifiles_n ${sta}/wav/hist $*
    set ${CASE}.ww3${inst_suffix}.i.*;                                                                                               dispose ifiles_y ${sta}/wav/init $*
    set ${CASE}.dwav${inst_suffix}.r.* ;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/wav/rest $*
    set ${CASE}.dwav${inst_suffix}.rs* ;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/wav/rest $*
    set ${CASE}.dwav${inst_suffix}.h.* ;                                                                                             dispose ifiles_n ${sta}/wav/hist $*
    IDX=`expr $IDX + 1`
done

#copy back the required files for next restart
cp ${sta}/rest/${dname}/* .

if mv ${sta}/* ${DOUT_S_ROOT}; then
    rm -fr ${sta}
else
    echo "st_archive.sh: error, final move command unsuccessful"
    echo "               some short-term archive data may be in ${sta}"
    echo "st_archive.sh: exiting"
    exit 1
fi

echo "st_archive.sh: short-term archiving completed successfully"

exit 0
