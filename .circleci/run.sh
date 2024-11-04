#!/bin/bash

case="master.WCYCL1850.ne4_oQU240.baseline"
compset="WCYCL1850"
res="ne4_oQU240"

log_dir="${HOME}/projects/e3sm/scratch/$case"

shopt -s extglob
shopt -s nullglob

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

function print_bldlog
{
    bldlog=$1
    bldlog_path=(${log_dir}/bld/$bldlog.bldlog.+([[:digit:]])-+([[:digit:]]))
    if [ -n "${bldlog_path}" ]; then
        command_yellow "cat ${bldlog_path}"
    fi
}

function print_runlog
{
    runlog=$1
    runlog_path=(${log_dir}/run/$runlog.log.+([[:digit:]])-+([[:digit:]]))
    if [ -n "${runlog_path}" ]; then
        command_yellow "cat ${runlog_path}"
        return 1
    else
        return 0
    fi
}

command_yellow_rc "mkdir -p $HOME/projects/e3sm/cesm-inputdata"
command_yellow_rc "cd cime/scripts"
command_yellow_rc "./create_newcase --case $case --compset $compset --res $res"
command_yellow_rc "cd $case"
command_yellow_rc "./xmlchange GMAKE_J=2"
command_yellow_rc "./case.setup"
command_yellow "./case.build -v"
rc=$?
if [ $rc -ne 0 ]; then
    print_bldlog "e3sm"
    print_bldlog "csm_share"
    print_bldlog "spio"
    print_bldlog "mct"
    print_bldlog "gptl"
    exit $rc
fi

# Unfortunately, case.submit always return 0, even if it fails.
# To find out if it succeeded, we check if the run log has been gzipped.
#command_yellow "./case.submit"
#print_runlog "e3sm"
#exit $?
