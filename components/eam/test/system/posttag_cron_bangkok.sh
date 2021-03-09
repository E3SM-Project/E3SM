#!/bin/sh 
# script currently set to run as nightly cron job on bangkok

tag_to_grab=`svn cat \
    https://svn-ccsm-models.cgd.ucar.edu/cam1/trunk/models/atm/cam/doc/ChangeLog | \
    grep 'Tag name:' | head -1`
tag_to_grab=${tag_to_grab##*:}
tag_to_grab=`echo ${tag_to_grab}`

if [ -z ${tag_to_grab} ]; then
    echo "ERROR getting most recent tag name  ...repository offline?...someone forget a tagname?"
    exit 1
fi

collections=/fs/cgd/csm/models/atm/cam

if [ ! -d ${collections}/${tag_to_grab} ]; then
    cd ${collections}
    svn export https://svn-ccsm-models.cgd.ucar.edu/cam1/trunk_tags/${tag_to_grab}
    chmod -R a-w ${tag_to_grab}

    tag_for_baseline=`svn cat \
	https://svn-ccsm-models.cgd.ucar.edu/cam1/trunk/models/atm/cam/doc/ChangeLog | \
	grep 'Tag name:' | head -2 | tail -1`
    tag_for_baseline=${tag_for_baseline##* }
    tag_for_baseline=`echo ${tag_for_baseline}`

    if [ ! -d ${collections}/${tag_to_grab} ] || \
	[ ! -d ${collections}/${tag_for_baseline} ]; then
	echo "ERROR collections not up to date"
	exit 1
    fi

    #prepare results directory for new tag
    results_dir=${collections}/test_results/${tag_to_grab}_bangkok
    if [ -d $results_dir ]; then
	rm -rf $results_dir
    fi
    mkdir $results_dir
    cd $results_dir

    export PATH=$PATH:/usr/local/torque/bin:/usr/local/bin
    env BL_ROOT=${collections}/${tag_for_baseline} \
	CAM_ROOT=${collections}/${tag_to_grab} \
	CAM_FC=PGI \
	CAM_INPUT_TESTS=${collections}/${tag_to_grab}/models/atm/cam/test/system/tests_posttag_bangkok \
	${collections}/${tag_to_grab}/models/atm/cam/test/system/test_driver.sh -f
else
    echo "collections/testing up to date"
fi

exit 0
