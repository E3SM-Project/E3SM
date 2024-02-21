#!/bin/bash -xe

echo "RUNNING CLEANUP FOR WORKSPACE: ${WORKSPACE}, BUILD_ID: ${BUILD_ID}"

cd ${WORKSPACE}

NUM_KEEP=12 # Adjust this number to keep more/fewer builds
KEEP_LAST=${BUILD_ID}
KEEP_FIRST=$((${BUILD_ID}-${NUM_KEEP}))
KEEP="$(seq ${KEEP_FIRST} 1 ${KEEP_LAST})"

echo "Keeping only builds in range [${KEEP_FIRST} ${KEEP_LAST}]."
REMOVE_THESE="$(ls -1 | grep -vF "${KEEP}")"
echo "Purging old builds: ${REMOVE_THESE}."

/bin/rm -rf $REMOVE_THESE

# Now clean up the scratch area
if [[ "$NODE_NAME" == "mappy" ]]; then
    # Ensure we have a newer python
    source scripts/jenkins/${NODE_NAME}_setup

    $JENKINS_SCRIPT_DIR/scratch_cleanup.py
fi
