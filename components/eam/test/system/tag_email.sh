#!/bin/sh 
# script currently set to run every 10min as cron job on poorman

admin=mvr@cgd.ucar.edu

cd /fs/cgd/csm/models/atm/cam/doc

svn diff -r HEAD --diff-cmd /usr/bin/diff -x "-b" ChangeLog > ./body
rc=$?
if [ $rc -eq 0 ]; then
    bytes=`wc -c < ./body`
    if [ "$bytes" -ne 0 ]; then
	subject=`grep 'Tag name:' ./body | head -1`
	subject=${subject##*:}
	subject=`echo ${subject}`
	svn list https://svn-ccsm-models.cgd.ucar.edu/cam1/trunk_tags/${subject} > /dev/null
	rc=$?
	if [ $rc -eq 0 ]; then
	    echo ""                                                            > ./message
	    echo "*** RESPONSES TO THIS EMAIL WILL NOT BE READ ***"           >> ./message
	    echo ""                                                           >> ./message
	    cat ./body                                                        >> ./message
	    env NAME='CAM Gatekeeper' mail -s ${subject} cam-dev@cgd.ucar.edu  < ./message

	    rc=$?
	    if [ $rc -eq 0 ]; then
		svn update
		rm /home/cam/.vacation.db
	        cp ChangeLog /web/public_html/cam/versions/.
	        cp ./body /web/public_html/cam/cam_checkins/Documentation/${subject}.diffs
	    else
		echo "Error from mail= $rc"
		echo "Error from mail= $rc" \
		    | env NAME='CAM Gatekeeper' mail -s "cron script error" $admin
	    fi
	    rm ./message
	else
	    echo "Error from svn list= $rc ...ChangeLog committed but $subject not tagged?"
	    echo "Error from svn list= $rc ...ChangeLog committed but $subject not tagged?" \
		| env NAME='CAM Gatekeeper' mail -s "cron script error" $admin
	fi
    else
	echo "ChangeLog unmodified"
    fi
else
    echo "Error from svn diff= $rc ...repository offline?"
    echo "Error from svn diff= $rc ...repository offline?" \
	| env NAME='CAM Gatekeeper' mail -s "cron script error" $admin
fi
rm ./body

exit 0
