#!/bin/bash

case="master.A_WCYCL1850.ne4_oQU240.baseline"
compset="A_WCYCL1850"
res="ne4_oQU240"

bldlog_dir="${HOME}/projects/e3sm/scratch/$case/bld"

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
    bldlog_path=(${bldlog_dir}/$bldlog.bldlog.+([[:digit:]])-+([[:digit:]]))
    if [ -n "${bldlog_path}" ]; then
        command_yellow "cat ${bldlog_path}"
    fi
}

command_yellow_rc "mkdir -p $HOME/projects/e3sm/cesm-inputdata"
command_yellow_rc "cd cime/scripts"
command_yellow_rc "./create_newcase --case $case --compset $compset --res $res"
command_yellow_rc "cd $case"
command_yellow_rc "./case.setup"
command_yellow "./case.build"
rc=$?
if [ $rc -ne 0 ]; then
    print_bldlog "csm_share"
    print_bldlog "pio"
    print_bldlog "mct"
    print_bldlog "gptl"
    print_bldlog "e3sm"
    exit $rc
fi
