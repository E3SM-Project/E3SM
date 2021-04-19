#!/bin/bash -xe

SCRIPT_DIR=$( cd "$( dirname "$0" )" && pwd )
DATE_STAMP=$(date "+%Y-%m-%d_%H%M%S")

set -o pipefail
$SCRIPT_DIR/jenkins_common_impl.sh 2>&1 | tee JENKINS_$DATE_STAMP
