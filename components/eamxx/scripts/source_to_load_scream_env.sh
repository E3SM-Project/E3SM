#
# A script for loading the scream-approved env for your machine.
# Requires SCREAM_MACHINE to be set
#

DIR=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

load_scream () { eval $($DIR/scream-env-cmd $1); }

if [ -z "$SCREAM_MACHINE" ]; then
    echo "Must set SCREAM_MACHINE"
else
    load_scream $SCREAM_MACHINE
fi
