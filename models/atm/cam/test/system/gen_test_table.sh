#!/bin/sh 
#

# this script, when executed in the directory containing the test-driver 
# scripts (~models/atm/cam/test/system) will loop through the default test 
# lists for pre and post tag testing of cam and create an html file 
# (test_table.html) with the specifics of each test detailed


echo "<html>" > ./test_table.html
echo "<head>" >> ./test_table.html
echo "<TITLE>CAM Testing Information Page</TITLE>" >> ./test_table.html
echo "</head>" >> ./test_table.html
echo "<body  BGCOLOR=\"666699\" TEXT=\"black\" LINK=\"white\" VLINK=\"FF9933\">" >> ./test_table.html

#########################################################################################
for input_file in `ls tests_*` ; do
    echo "<TABLE border=2 width=750>" >> ./test_table.html
    if [ $input_file = "tests_posttag_bangkok" ] || [ $input_file = "tests_pretag_bangkok" ]; then
	echo "<CAPTION>${input_file}*</CAPTION>" >> ./test_table.html
    else
	echo "<CAPTION>$input_file</CAPTION>" >> ./test_table.html
    fi
    echo "<TR>" >> ./test_table.html
    echo "<TH>test# </TH>" >> ./test_table.html
    echo "<TH>testid </TH>" >> ./test_table.html
    echo "<TH>test script </TH>" >> ./test_table.html
    echo "<TH>arg1 </TH>" >> ./test_table.html
    echo "<TH>arg2 </TH>" >> ./test_table.html
    echo "<TH>arg3 </TH>" >> ./test_table.html
    echo "<TH>arg4 </TH>" >> ./test_table.html
    echo "<TH>arg5 </TH>" >> ./test_table.html
    echo "</TR>" >> ./test_table.html

    test_list=""
    while read input_line; do
	test_list="${test_list}${input_line} "
    done < ./${input_file}

    count=0
    ##loop through the tests of input file
    for test_id in ${test_list}; do
	echo "<TR>" >> ./test_table.html
	count=`expr $count + 1`
	while [ ${#count} -lt 3 ]; do
		count="0${count}"
	done
	echo "<TD> $count </TD>" >> ./test_table.html

	master_line=`grep $test_id ./input_tests_master`
	for arg in ${master_line}; do
	    if [ -f ./nl_files/$arg ]; then
		echo "<TD><A HREF=\"./nl_files/$arg\">$arg </A></TD>" >> ./test_table.html
	    elif [ -f ./config_files/$arg ]; then
		echo "<TD><A HREF=\"./config_files/$arg\">$arg </A></TD>" >> ./test_table.html
	    else
		echo "<TD>$arg </TD>" >> ./test_table.html
	    fi
	done
	echo "</TR>" >> ./test_table.html
    done
    echo "</TABLE>" >> ./test_table.html
    echo "<pre>" >> ./test_table.html
    echo " " >> ./test_table.html
    echo "</pre>" >> ./test_table.html
done
echo "<pre>" >> ./test_table.html
echo "* post-tag testing on bangkok is done with pgi compilers, pre-tag with lahey " >> ./test_table.html
echo "</pre>" >> ./test_table.html

exit 0
