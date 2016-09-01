#! /bin/bash


while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        start )
	    CASE_DIR="$2"
	    shift # past argument
	    cd $CASE_DIR
	    ./case.setup
	    ./case.build
	    ./case.submit
	    ;;
	setup )
	    CASE_DIR="$2"
	    shift # past argument
	    cd $CASE_DIR
	    ./case.setup
	    ;;
	*)
            # unknown option
	;;
    esac
	    echo shift # past argument or value
done
