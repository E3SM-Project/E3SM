#!/usr/bin/env bash

usage () {
    cat <<EOF
find_mergeinfo.sh [options] directory

Recursively print all files and directories with svn:mergeinfo.

OPTIONS

-d :: Only print directories with svn:mergeinfo.
-h :: Print this help message.
-r :: Recursively descend into svn:externals as well as the directories in
      the current working copy.
-s :: Subdirectories only: skip the top level directory and all top level
      svn:external directories.

LIMITATIONS

Currently, this script only works with old-style external specifications,
where the local directory is on the left, and the repository URL is on the
right.

svn:externals detection is messed up on directories with paths containing a
'!' character.

This is written for Bash 4, and won't work with POSIX shell or earlier Bash
versions.

EOF
}

new_args=$(getopt -o dhrs -- "$@")
if (( $? != 0 )); then
    usage >&2
    exit 1
fi
eval set -- "$new_args"

while (( $# > 0 )); do
    case "$1" in
        -d)
            check_dir=1
            ;;
        -h)
            usage
            exit
            ;;
        -r)
            recursive=1
            ;;
        -s)
            sub_only=1
            ;;
        --)
            shift
            break
            ;;
    esac
    shift
done

if (( $# != 1 )); then
    echo "You must provide exactly one argument." >&2
    usage >&2
    exit 1
fi

find_external_dirs () {
    local more_dirs
    local more_more_dirs
    # Regex below only valid with GNU sed. Use "!" as delimiter because it
    # is unlikely to appear in the path.
    get_external_regex='s!\(\S\+\)\b\s\+.*!'"$1"'/\1!'
    more_dirs=$(svn pg svn:externals "$1" 2>/dev/null | \
        sed "$get_external_regex")
    for dir in $more_dirs; do
        more_more_dirs=$(find_external_dirs "$dir")
    done
    echo "$more_dirs $more_more_dirs"
}

# Get absolute path (to make -s option work robustly).
svn_dirs=`cd $1 && pwd`

if [[ "$recursive" ]]; then
    svn_dirs+=" $(find_external_dirs $svn_dirs)"
fi

find_dir_mergeinfo () {
# Parse recursive listing of svn:mergeinfo, returning all directories
# that have this property.
# Should use pl instead of pg?
# Subversion 1.6 always returns absolute paths from pg, but Subversion 1.7
# can return relative paths. Hack around this by removing "${PWD}/" if it
# is present from both $top_dir and $tmp_path, so that hopefully they match
# up.
    top_dir="${1#${PWD}/}"
    while read line; do
        if [[ "${line}" =~ (.*)\ -\ .* ]]; then
            tmp_path="${BASH_REMATCH[1]#${PWD}/}"
            if [[ ! "$check_dir" || -d $tmp_path ]]; then
                if [[ ! "$sub_only" || "$tmp_path" != "$top_dir" ]]; then
                    echo "$tmp_path"
                fi
            fi
        fi
    done < <(svn pg -R svn:mergeinfo "$top_dir")
}

for dir in $svn_dirs; do
    find_dir_mergeinfo "$dir"
done
