#!/bin/sh
# Test for bad svn:mergeinfo
# Ensures that the top directory has mergeinfo, and nothing else.
# This also verifies that CAM externals do not have mergeinfo below the
# root directory (except HOMME).

# Return codes in use:
# 1: Extra mergeinfo
# 2: Error from running an external command
# 4: No mergeinfo on top level

# Utility to check return code.
# Give it the code and an error message, and it will print stuff and exit.
check_code () {
    if [ "$1" -ne 0 ]; then
        echo "Error: return code from command was $1"
        echo "$2"
        exit 2
    fi
}

# Little utility for finding absolute path to a directory.
get_dir_abspath () {
    echo $(cd $1 && pwd)
}

cam_top_dir=$(get_dir_abspath ${CAM_SCRIPTDIR}/../..)
cesm_top_dir=$(get_dir_abspath $cam_top_dir/../../..)

rc=0

# Check to make sure that the top level directory *does* have mergeinfo.
top_dir_mergeinfo=$(svn pg svn:mergeinfo "$cesm_top_dir")

check_code "$?" "Problem running svn pg on the CAM root directory."

if [ "${#top_dir_mergeinfo}" -lt 1 ]; then
    cat <<EOF
The top directory is missing svn:mergeinfo, which should be restored. Try:
1) Reverting the local change that removed it.
2) Copying from the last revision of this branch that had mergeinfo.
3) Copying mergeinfo from the root directory of the CAM trunk.

EOF
    rc=4
fi

# Check list of mergeinfo files.
handle_mergeinfo () {
    if [ "${#1}" -gt 0 ]; then
        cat <<EOF
Search discovered that the following files have svn:mergeinfo, which should
be removed:
$1

EOF
        rc=1
    fi
}

# Check for mergeinfo in CAM; use "-s" to ignore top level.
mergeinfo_files=$(${CAM_SCRIPTDIR}/find_mergeinfo.sh -s "$cesm_top_dir")

check_code "$?" "Problem running find_mergeinfo.sh on CAM."

handle_mergeinfo "$mergeinfo_files"

# Check for mergeinfo in CAM externals.
# Leave out the HOMME external for now, since CAM devs don't manage it.
# The grep uses extended regex to rule out matches to either absolute or
# relative path.
mergeinfo_files=$(${CAM_SCRIPTDIR}/find_mergeinfo.sh -rs "$cam_top_dir" | \
    grep -Ev "($cam_top_dir|${cam_top_dir#${PWD}/})/src/dynamics/se/share")

# Not checking code because grep will normally return 1.

handle_mergeinfo "$mergeinfo_files"

if [ "$rc" -eq 1 ]; then
    cat <<EOF

You can delete svn:mergeinfo by running the following commands at the top
level of the CAM source tree (or the CAM external with mergeinfo):
> svn pd -R svn:mergeinfo
> svn revert .

If you have changed the externals, they will be reverted, so you will have
to set the svn:externals property again after using these commands.

EOF
fi

if [ "$rc" -eq 0 ]; then
    echo "No problems found in svn:mergeinfo in this working copy."
fi

exit $rc
