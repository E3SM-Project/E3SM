#!/bin/bash -xe

# Adjust this number to keep more/less builds
NUM_KEEP=30
KEEP_LAST=${BUILD_ID}
KEEP_FIRST=$((${BUILD_ID}-${NUM_KEEP}))
KEEP="$(seq ${KEEP_FIRST} 1 ${KEEP_LAST})"

cd ${WORKSPACE}
echo "Keeping only builds in range [${KEEP_FIRST} ${KEEP_LAST}]."
REMOVE_THESE="$(ls -1 | grep -vF "${KEEP}")"
echo "Purging old builds: ${REMOVE_THESE}."

/bin/rm -rf $REMOVE_THESE
