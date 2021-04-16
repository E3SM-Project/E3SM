#!/bin/bash -xe

SCRIPT_DIR=$( cd "$( dirname "$0" )" && pwd )
DATE_STAMP=$(date "+%Y-%m-%d_%H%M%S")

$SCRIPT_DIR/jenkins_common_impl.sh >& JENKINS_$DATE_STAMP
