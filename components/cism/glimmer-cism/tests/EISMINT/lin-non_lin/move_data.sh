# simple shell script to move data file if it exists
infile=$1
outfile=$2
if [ -f $infile ] 
then
 mv $infile $outfile
fi
