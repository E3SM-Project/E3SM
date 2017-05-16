
#!/bin/sh
# script to archive case directory to subversion repository
#
# aliceb - 11/1/11
# commented out prompts for case metadata.
# Case metadata is now being entered and stored in
# mysql run db at http://csegweb.cgd.ucar.edu.
#
# aliceb - 2013/5/30
# adding the archival of user_nl_xxx files to the archive list

if [ $# -gt 0 ]; then
    echo ""
    echo "NAME"
    echo "archive_metadata.sh:  archive information from case directory"
    echo "  to the ccsm run repository - will import the information if"
    echo "  the case is new to the repository, or commit to the trunk"
    echo "  of an existing case"
    echo ""
    echo "OPTIONS"
    echo "no options currently enabled"
    echo ""
    echo "EXAMPLES"
    echo "./archive_metadata.sh"
    echo "./archive_metadata.sh -help"
    exit 0
fi

while [ -z "$EDITOR" ] || [ ! -x $EDITOR ]; do
    echo "A text editor is required for executing this script and the environment"
    echo "variable EDITOR was unset or invalid;  Enter full path for EDITOR (default is vi)"
    read EDITOR
    if [ -z "$EDITOR" ]; then
	EDITOR=`which vi`
    fi
done

#get case name
string_to_parse=`ls *.build`
casename=${string_to_parse%\.*}

#generate date strings for use later
cur_time=`date '+%H:%M:%S'`
cur_date=`date`

svnrepo="https://svn-ccsm-rundb.cgd.ucar.edu"
svnstr="${svnrepo}/${casename}"

no_svn="FALSE"
#look for subversion client
svn --version > /dev/null
if [ $? -ne 0 ]; then
    no_svn="TRUE"
    echo ""
    echo "Error: the subversion client svn was not found"
else
    echo "your subversion password may be requested..."
    svn list $svnrepo > /dev/null
    if [ $? -ne 0 ]; then
	no_svn="TRUE"
	echo ""
	echo "Error: subversion username/password was invalid"
    fi
fi

if [ "$no_svn" == "TRUE" ]; then
    echo "Access to the subversion run database repository is required for this utility"
    echo ""
    exit 99
fi

#look for existence of this case's trunk in repository to see if case should be imported
#or committed to an existing trunk
svn list $svnstr/trunk > /dev/null
if [ $? -ne 0 ]; then

    #no info for this case, but it may be reserved in repository
    svn list $svnstr > /dev/null
    if [ $? -eq 0 ]; then
	#get contact for reserved casename
	line=`svn log $svnrepo/$casename | grep 'reserved by:'`
	reserver=${line##*:}
	reserver=`echo ${reserver}`

	entry="NOT valid"
	while [ "$entry" != "valid" ]; do

	    echo "The casename $casename has been reserved by $reserver..."
            echo "enter 'cont' if you acknowledge using this reservation and wish to continue or 'quit' to exit"
	    read ans
	    if [ "$ans" == 'cont' ]; then
		entry="valid"
	    elif [ "$ans" == 'quit' ]; then
		exit 99
	    fi
	done
    fi

    #prompt for metadata
#
# aliceb - metadata all stored in the rundb at http://csegweb.cgd.ucar.edu/cgi-bin/index.cgi
#
#    entry="NOT valid"
#    while [ "$entry" != "valid" ]; do
#	echo "new case metadata to be archived - please enter a title for $casename - limit to 80 chars"
#	read title
#	if [ ${#title} -le 80 ]; then
#	    entry="valid"
#	else
#	    echo "too many characters - try again"
#	fi
#    done

#    entry="NOT valid"
#    while [ "$entry" != "valid" ]; do
#	echo "please enter an email address for the person acting as a contact for this model run (user@domain)"
#	read contact
#	entry=`echo $contact | perl -e 'if ($ll = <> =~ /^[^\s@]+@[^\s@]+$/) \
#	    { print "valid" }'`
#    done

#    entry="NOT valid"
#    while [ "$entry" != "valid" ]; do
#	echo "enter the responsible ccsm working group"
#	echo "use [AMWG  LMWG  OMWG  PCWG  PWG  BGCWG  CVWG  CWG  SEWG  OTHER]"
#	read wg
#	case $wg in
#	    AMWG | LMWG | OMWG | PCWG | PWG | BGCWG | CVWG | CWG | SEWG | OTHER )
#	    entry="valid"
#	    ;;
#	esac
#    done

#    entry="NOT valid"
#    while [ "$entry" != "valid" ]; do
#	echo "enter the date of the first run submission (yyyy-mm-dd)"
#	read submit_dt
#	entry=`echo $submit_dt | perl -e 'if ($ll = <> =~ /^\d\d\d\d-(\d\d)-(\d\d)$/) \
#	    { print "valid" if $1>0 && $1<13 && $2>0 && $2<32 }'`
#    done

    #get run startdate for documentation to changelog
#    entry="NOT valid"
#    while [ "$entry" != "valid" ]; do
#	echo "enter the model run start date yyyy-mm-dd"
#	read effective_model_dt
#	entry=`echo $effective_model_dt | perl -e 'if ($ll = <> =~ /^\d\d\d\d-(\d\d)-(\d\d)$/) \
#	    { print "valid" if $1>0 && $1<13 && $2>0 && $2<32 }'`
#    done

#    echo "...................      DESCRIPTION FOR $casename      ..................." > description_for_${casename}
#    $EDITOR description_for_${casename}
    echo "thank you - will now attempt import"

    #dump metadata to a file
#    cat description_for_${casename} | perl -e 'while (my $ll = <>) { print STDOUT "# $ll" }' \
#                                                                              > env_metadata
#    rm description_for_${casename}*
#    echo " "                                                                 >> env_metadata
#    echo "setenv    CASE_TITLE            \"$title\""                        >> env_metadata
#    echo "setenv    CASE_CONTACT          \"$contact\""                      >> env_metadata
#    echo "setenv    CASE_WORKING_GROUP    \"$wg\""                           >> env_metadata
#    echo "setenv    CASE_1ST_SUBMIT_DT    \"$submit_dt\""                    >> env_metadata

    #initial tag name, will increment from here on
    casetag=${casename}_0001

    echo "************************************************************************" >  ChangeLog
    echo "Tag Name:                              $casetag"                          >> ChangeLog
    echo "Date:                                  $cur_date"                         >> ChangeLog
    echo "Model date where changes take effect:  $effective_model_dt"               >> ChangeLog
    echo " "                                                                        >> ChangeLog
    echo "...................      EXPLANATION OF CHANGES      ..................." >> ChangeLog
    echo "initial import of case information"                                       >> ChangeLog
    echo " "                                                                        >> ChangeLog

    #make temp directory with unique name for staging files for import to svn
    tempdir="tempdir_$cur_time"
    mkdir -p $tempdir/trunk
    mkdir -p $tempdir/trunk_tags


    #set list of files/directories to be archived by default
    archive_list=""
    archive_list="${archive_list} Buildconf"
    archive_list="${archive_list} CaseDocs"
    archive_list="${archive_list} LockedFiles"
    archive_list="${archive_list} README.science_support"
    archive_list="${archive_list} README.case"
    archive_list="${archive_list} SourceMods"
    archive_list="${archive_list} Tools"
    archive_list="${archive_list} archive_metadata.sh"
    archive_list="${archive_list} case.build"
    archive_list="${archive_list} case.run"
    archive_list="${archive_list} check_input_data"
    archive_list="${archive_list} case.setup"
    archive_list="${archive_list} create_production_test"
    archive_list="${archive_list} env_build.xml"
    archive_list="${archive_list} env_case.xml"
    archive_list="${archive_list} env_derived"
    archive_list="${archive_list} env_mach_pes.xml"
    archive_list="${archive_list} env_mach_specific"
    archive_list="${archive_list} env_run.xml"
    # archive_list="${archive_list} logs"
    archive_list="${archive_list} timing"
    archive_list="${archive_list} xmlchange"
    # aliceb archive_list="${archive_list} env_metadata"
    archive_list="${archive_list} ChangeLog"
    # adding user_nl_xxx files
    archive_list="${archive_list} user_nl_cam"
    archive_list="${archive_list} user_nl_cice"
    archive_list="${archive_list} user_nl_clm"
    archive_list="${archive_list} user_nl_cpl"
    archive_list="${archive_list} user_nl_pop"
    archive_list="${archive_list} user_nl_rtm"
    archive_list="${archive_list} user_nl_mosart"

    for item in ${archive_list}; do
	if [ -e $item ]; then
    	  cp -rp $item $tempdir/trunk/.
	  if [ $? -ne 0 ]; then
	    echo "ERROR: preparing '$item' for archival but file or directory had problems"
	    rm -rf $tempdir
	    exit 1
	  fi
         fi
    done

    archive_o_files="FALSE"
    echo "stderr/stdout files are not required for archival - archive anyway? (y/n - default is no)"
    read ans
    case $ans in
	Y* | y* )
        archive_o_files="TRUE";;
    esac

    file_list=`ls -A`
    for item in ${file_list}; do
	if [ $item != $tempdir ]; then
	    if [ ! -e ${tempdir}/trunk/${item} ]; then
                case $item in
                    poe* | ${casename}.o* | ${casename}_la.o* )
                    if [ $archive_o_files == "TRUE" ]; then
			cp -rp $item $tempdir/trunk/.
			archive_list="$archive_list $item"
		    fi;;

                    * )
                    echo "'$item' is not required for archival - archive anyway? (y/n - default is no)"
		    read ans
		    case $ans in
			Y* | y* )
                        cp -rp $item $tempdir/trunk/.
		        archive_list="$archive_list $item";;
                    esac;;
		esac
            fi
        fi
    done

    cd $tempdir
    svn import . $svnstr/ -m "initial import of case info for $casename"
    rc=$?

    cd ..
    rm -rf $tempdir

    if [ $rc -ne 0 ]; then
	echo "ERROR: trouble importing $casename to repository "
	exit 1
    else
	#convert current sandbox to a subversion sandbox
	for item in ${archive_list}; do
	    rm -rf $item
	done
	svn co $svnstr/trunk . > /dev/null

	echo "successfully imported case information for $casename"
    fi

else
    #check to see if we have an svn sandbox
    svn status > /dev/null
    if [ $? -ne 0 ]; then
	echo " "
	echo "ERROR: $casename exists in the case repository, but the current "
	echo "directory is not a working copy of the development trunk from "
	echo "the repository.  You must first check out the latest version of "
	echo "the trunk and then add your changes before trying to archive: "
	echo " "
	echo "example> svn co https://svn-ccsm-rundb.cgd.ucar.edu/<casename>/trunk  <casename>"
	echo " "
	exit 99
    fi

    #check to see if we are referencing the trunk before an attempt to commit
    svn info | grep "^URL: $svnstr/trunk\$" > /dev/null
    if [ $? -ne 0 ]; then
	echo " "
	echo "ERROR: the current directory is not a working copy of the "
	echo "development trunk for $casename.  You must first check out the "
	echo "latest version of the trunk and then add your changes before "
	echo "trying to archive.  You are currently referencing: "
	echo " "
	svn info | grep URL
	echo " "
	exit 99
    fi

    #check to see if the trunk has changed in the repository since being checked out
    result=`svn status -u | perl -e 'while (my $ll = <>) \
	{ if ($ll =~ /^.......\*/) \
	{ print "update_needed" }}'`
    if [ -n "$result" ]; then
	echo " "
	echo "ERROR: the current directory is now out-of-date with the HEAD "
	echo "revision of the development trunk for $casename in the case "
	echo "repository.  Please resolve before trying to archive"
	echo " "
	exit 99
    fi

    #ask user if unrecognized files should be archived
    result=`svn status | perl -e 'while (my $ll = <>) \
	{ if ($ll =~ /^\?\s+([^\s]+)/) \
	{ print "$1 " }}'`
    for item in ${result}; do
	echo "'$item' has not previously been archived - archive now? (y/n - default is no)"
	read ans
	case $ans in
	    Y* | y* )
	    svn add $item;;
	esac
    done

    #ask user if missing files should be deleted from archive
    result=`svn status | perl -e 'while (my $ll = <>) \
	{ if ($ll =~ /^\!\s+([^\s]+)/) \
	{ print "$1 " }}'`
    for item in ${result}; do
	echo "$item has previously been archived but is now missing - delete from archive? (y/n - default is no)"
	read ans
	case $ans in
	    Y* | y* )
	    svn delete $item;;
	esac
    done

    #check for messed up sandbox - missing files, files in conflict
    result=`svn status | perl -e 'while (my $ll = <>) \
	{ if ($ll =~ /^[C\!]/) \
	{ print "conflict " }}'`
    if [ -n "$result" ]; then
	echo " "
	echo "ERROR: archival of this directory tree cannot continue due "
	echo "to a conflict or missing file.  Please review this output "
	echo "of 'svn status' to identify and correct the problem before "
	echo "attemping to archive: "
	echo " "
	svn status
	echo " "
	#revert files scheduled for addition/deletion?
	exit 99
    fi

    #check to see if there are any changes to commit
    svn status | perl -e 'while (my $ll = <>) { if ($ll =~ /^([MAD].+)$/) { print STDOUT "$1\n" }}' \
                                                                                    > stat_$cur_time
    bytes=`wc -c < ./stat_$cur_time`
    if [ "$bytes" -eq 0 ]; then
	echo " "
	echo "ERROR: there are no detected modifications to the files "
	echo "or directories currently being archived. "
	echo " "
	rm stat_$cur_time*
	exit 99
    fi

    #generate name for new tag - just incrementing for now
    result=`svn list $svnstr/trunk_tags | tail -1`
    result=${result##*_}
    result=${result%/}
    result=`expr $result + 1`
    if [ $result -lt 10 ]; then
	result="000${result}"
    elif [ $result -lt 100 ]; then
	result="00${result}"
    elif [ $result -lt 1000 ]; then
	result="0${result}"
    elif [ $result -gt 9999 ]; then
	echo "ERROR: maximum number of tags exceeded for $casename"
	exit 99
    fi
    casetag=${casename}_${result}

    #prompt for effective_model_date and explanation of mods
    entry="NOT valid"
    while [ "$entry" != "valid" ]; do
	echo "enter the model date where changes take effect yyyy-mm-dd"
	read effective_model_dt
	entry=`echo $effective_model_dt | perl -e 'if ($ll = <> =~ /^\d\d\d\d-(\d\d)-(\d\d)$/) \
	    { print "valid" if $1>0 && $1<13 && $2>0 && $2<32 }'`
    done

    echo "...................      EXPLANATION OF CHANGES      ..................." > xoc_$cur_time
    $EDITOR xoc_$cur_time

    #dump documentation to temporary file for eventual inclusion in ChangeLog
    echo "************************************************************************" >  log_$cur_time
    echo "Tag Name:                              $casetag"                          >> log_$cur_time
    echo "Date:                                  $cur_date"                         >> log_$cur_time
    echo "Model date where changes take effect:  $effective_model_dt"               >> log_$cur_time
    echo " "                                                                        >> log_$cur_time
    cat stat_$cur_time                                                              >> log_$cur_time
    rm stat_$cur_time*
    echo " "                                                                        >> log_$cur_time
    cat xoc_$cur_time                                                               >> log_$cur_time
    rm xoc_$cur_time*
    echo " "                                                                        >> log_$cur_time
    cat ChangeLog                                                                   >> log_$cur_time
    mv log_$cur_time ChangeLog

    svn commit . -m "committing mods to trunk of $casename"
    if [ $? -ne 0 ]; then
	echo "ERROR: problem commiting changes to the repository"
	svn revert ChangeLog
	#revert files scheduled for addition/deletion?
	exit 99
    else
	echo "successfully committed changes to case information for $casename"
    fi

    #won't be used since this won't be for initial tag, but need something for contact in svn tag message
    contact=$LOGNAME
fi

svn copy $svnstr/trunk $svnstr/trunk_tags/${casetag} -m "archive_metadata.sh: tag $casetag created for $casename by: $contact"
if [ $? -ne 0 ]; then
    echo "ERROR: trouble tagging $svnstr/trunk"
    exit 1
else
    echo "successfully tagged the repository trunk of $casename"
fi

exit 0
