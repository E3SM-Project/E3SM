#!/bin/bash


function command_yellow
{
    echo -e "\e[33m\$ $1\e[0m"
    eval $1
    return $?
}

function command_yellow_rc
{
    command_yellow "$1"
    rc=$?
    if [ $rc -ne 0 ]; then
        exit $rc
    fi
}

command_yellow_rc "mkdir -p $HOME/projects/e3sm/cesm-inputdata"
command_yellow_rc 'cd cime/scripts'
command_yellow_rc './create_newcase --case master.A_WCYCL1850.ne4_oQU240.baseline --compset A_WCYCL1850 --res ne4_oQU240'
command_yellow_rc 'cd master.A_WCYCL1850.ne4_oQU240.baseline'
command_yellow_rc './case.setup'
command_yellow './case.build'
rc=$?
if [ $rc -ne 0 ]; then
    log=`ls ${HOME}/projects/e3sm/scratch/master.A_WCYCL1850.ne4_oQU240.baseline/bld/pio.bldlog*`
    command_yellow "cat $log"
    exit $rc
fi
