#!/bin/sh 
#

# this script, when executed in the directory containing the test-driver 
# scripts (~/test/system) will loop through the default test 
# lists for pre and post tag testing of clm and create an html file 
# (test_table.html) with the specifics of each test detailed

outfile="./test_table.html"

echo '<?xml version="1.0" encoding="UTF-8"?>' > $outfile
echo '<!DOCTYPE html>' >> $outfile
echo '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">' >> $outfile
echo '<head>' >> $outfile
echo '<title>CLM Testing Information Page</title>' >> $outfile
echo '</head>' >> $outfile
echo '<body  bgcolor="666699" text="black" link="white" vlink="FF9933">' >> $outfile

#########################################################################################
for input_file in `ls tests_*` ; do
    echo '<table border="2" width="750">' >> $outfile
    echo "<caption>$input_file</caption>" >> $outfile
    echo "<tr>" >> $outfile
    echo "<th>test# </th>" >> $outfile
    echo "<th>testid </th>" >> $outfile
    echo "<th>test script </th>" >> $outfile
    echo "<th>arg1 </th>" >> $outfile
    echo "<th>arg2 </th>" >> $outfile
    echo "<th>arg3 </th>" >> $outfile
    echo "</tr>" >> $outfile

    test_list=""
    while read input_line; do
	test_list="${test_list}${input_line} "
    done < ./${input_file}

    count=0
    ##loop through the tests of input file
    for test_id in ${test_list}; do
	echo "<tr>" >> $outfile
	count=`expr $count + 1`
	while [ ${#count} -lt 3 ]; do
		count="0${count}"
	done
	echo "<td> $count </td>" >> $outfile

	master_line=`grep $test_id ./input_tests_master`
        dir=""
	for arg in ${master_line}; do
            arg1=${arg%^*}
            arg2=${arg#*^}
            if [ -d ../../tools/$arg ]; then
                dir=$arg
	    elif [ -f ./nl_files/$arg ]; then
		echo "<td><a href=\"./nl_files/$arg\">$arg </a></td>" >> $outfile
	    elif [ -f ./config_files/$arg ]; then
		echo "<td><a href=\"./config_files/$arg\">$arg </a></td>" >> $outfile
	    elif [ -f ./nl_files/$arg1 ] && [ -f ./nl_files/$arg2 ]; then
		echo  "<td><a href=\"./nl_files/$arg1\">$arg1</a>^" \
                        "<a href=\"./nl_files/$arg2\">$arg2</a></td>" >> $outfile
	    elif [ -f ./nl_files/$arg1 ] && [ -f ./config_files/$arg2 ]; then
		echo "<td><a href=\"./nl_files/$arg1\">$arg1</a>^" \
                        "<a href=\"./config_files/$arg2\">$arg2</a></td>" >> $outfile
	    elif [ -f ../../tools/$dir/$dir.$arg ]; then
		echo "<td><a href=\"../../tools/$dir/$dir.$arg\">$arg </a></td>" >> $outfile
	    else
		echo "<td>$arg </td>" >> $outfile
	    fi
	done
	echo '</tr>' >> $outfile
    done
    echo '</table>' >> $outfile
    echo '<pre>' >> $outfile
    echo ' ' >> $outfile
    echo '</pre>' >> $outfile
done
echo '</body>' >> $outfile
echo '</html>' >> $outfile

exit 0
