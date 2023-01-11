#!/bin/bash -xe

cd $JENKINS_SCRIPT_DIR/../..

source scripts/jenkins/${NODE_NAME}_setup

./scripts/scripts-ctest-driver -s -m $SCREAM_MACHINE
