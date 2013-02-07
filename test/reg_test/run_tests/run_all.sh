#!/bin/bash

# We'll use functions from the test utilities script
source ./test-utilities.sh

# First read the list of test files from test-list.in
readTestList

# Set some directory paths
setTestDirs

# Submit the tests and keep track of the job ids
submitTestsToLSF

# Print a summary of the submissions
printSubmissionSummary

# Wait for the jobs to run through the queue
queueWait

# Diff the output files with those saved in the repo
diffStdOut

# Exit happily
exit 0
