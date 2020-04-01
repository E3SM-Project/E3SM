#!/bin/bash -xe

source ./scream/components/scream/scripts/jenkins/${NODE_NAME}_setup

git config --global user.email "jenkins@ignore.com"
git config --global user.name "Jenkins Jenkins"

if [ -n "$PULLREQUESTNUM" ]; then
    cd ./scream
    git checkout -b pr/$PULLREQUESTNUM
    git submodule update --recursive
    cd -
fi

./scream/components/scream/scripts/gather-all-data './scripts/test-all-scream $compiler -p -m $machine' -l -m $SCREAM_MACHINE
