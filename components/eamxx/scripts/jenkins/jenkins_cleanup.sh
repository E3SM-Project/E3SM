#!/bin/bash -e

export JENKINS_SCRIPT_DIR=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)
DATE_STAMP=$(date "+%Y-%m-%d_%H%M%S")

$JENKINS_SCRIPT_DIR/jenkins_cleanup_impl.sh >& ${WORKSPACE}/${BUILD_ID}/JENKINS_CLEAN_$DATE_STAMP
