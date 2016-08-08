#!/bin/sh
#
#!
# @file support/dox-filter.sh
#
# @brief Convince Doxygen to parse files other than source code
# as part of documentation set.
#
# This script is invoked from with @b doxygen via @b support/bootjvm.dox
# as the @b INPUT_FILTER value.
#
# For shell scripts, strip non-comments and convert shell comments
# beginning in column 1 from (<b>^#</b>) to (<b> *</b>).  These appear
# to Doxygen to be 'C' style intermediate comments (of stylistic
# interest).
#
# Convert comment with an explanation point (<b>^#!</b>) character
# to look like a Doxygen start-of-documentation tag (<b>/</b><b>*!</b>),
# a special form of a 'C' style open comment (of syntactic interest).
#
# Convert comment with a slash (<b>^#/</b>) character to look like
# a Doxygen end-of-documentation tag (<b>*</b><b>/</b>), a special
# form of a 'C' style close comment (of syntactic interest).
#
# Other selected files may be converted also, but only the header areas
# need to begin with a script-style comment.  For a simple example,
# please refer to @link ./LICENSE LICENSE@endlink.
#
#
# @todo  HARMONY-6-support-dox-filter.sh-1 A Windows .BAT version of this
#        script needs to be written
#
#
# @section Control
#
# \$URL$
#
# \$Id$
#
# Copyright 2005 The Apache Software Foundation
# or its licensors, as applicable.
#
# Licensed under the Apache License, Version 2.0 ("the License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied.
#
# See the License for the specific language governing permissions
# and limitations under the License.
#
# @version \$LastChangedRevision$
#
# @date \$LastChangedDate$
#
# @author \$LastChangedBy$
#
#         Original code contributed by Daniel Lydick on 09/28/2005.
#
# @section Reference
#
#/ /* 
# (Use  #! and #/ with dox-filter.sh to fool Doxygen into
# parsing this non-source text file for the documentation set.
# Use the above open comment to force termination of parsing
# since it is not a Doxygen-style 'C' comment.)
#
#
###################################################################
#
# Look for all INPUT= files with '.sh' extension and convert them
# to look like a block of 'C' code comments so Doxygen can parse
# out their documentation tags.
#

# Magic, but _very_ vanilla Unix, method to extract file extension
FILENAME=`expr "/${1:-.}" : \
           '\(.*[^/]\)/*$' : \
           '.*/\(..*\)' : \
           "\\(.*\\)sh\$" `

FILEEXT=`expr "/${1:-.}" : \
           '\(.*[^/]\)/*$' : \
           '.*/\(..*\)' : \
           "$FILENAME\\(.*\\)\$" `

# Process '*.sh' and specific files with filter to convert their
# comments per above

convertit=0
if test "sh" = "$FILEEXT"
then
    convertit=1
else
    case $FILENAME in
        AUTHORS | INSTALL | LICENSE | README) convertit=2;;
        RELEASE_LEVEL)                        convertit=3;;
        Makefile | MakeSetup | MakeRules)     convertit=4;;
        *);;
    esac
fi

case $convertit in
    0)
       cat $1;;

    1) # Read file, strip /bin/sh line, convert comments
       cat $1 | \
       (read line1; cat -) | \
       grep "^#" | \
       sed 's,^#!,/*!,;s,^#/, */,;s,^#, *,'
       ;;

    2 | 3 | 4)
       # Read file, convert whole file of which only
       # header areas need to contain comments.
       cat $1 | \
       sed 's,^#!,/*!,;s,^#/, */,;s,^#, *,'
esac
exit 0
#
# EOF
