#!/bin/bash -xe

source ./scream/components/scream/scripts/jenkins/${NODE_NAME}_setup

git config --global user.email "jenkins@ignore.com"
git config --global user.name "Jenkins Jenkins"

SUBMIT="--submit"
if [ -n "$PULLREQUESTNUM" ]; then
    SUBMIT="" # We don't submit AT runs
fi

./scream/components/scream/scripts/gather-all-data "./scripts/test-all-scream \$compiler -p -i -m \$machine $SUBMIT" -l -m $SCREAM_MACHINE
