#!/bin/bash

## Argument check
if [ $# -ne 1 ]; then
    printf "Syntax: cleandir.x <CLEANMODE>\n \
       cleandir.x 1   : deletes auxiliary data files generated\n \
       cleandir.x 2   : deletes log files\n \
       cleandir.x 3   : deletes task files and folders\n"
    exit
fi

CLEANMODE=$1

echo "Cleaning .pyc files"
\rm -f *.pyc

if [ "$CLEANMODE" -ge "1" ]; then
	echo "Cleaning data files generated along the way"
	\rm -f sigma2* *.dat
fi
if [ "$CLEANMODE" -ge "2" ]; then
	echo "Cleaning .log files"
	\rm -f *.log
fi
if [ "$CLEANMODE" -ge "3" ]; then
	echo "Cleaning task files and folders"
	\rm -rf tasks task_*
fi



#\rm `find . ! -type d | egrep -v "CVS|ibcs.x|cleandir.x"`

