#!/bin/bash

plid=$1
echo -n"" > tasks

MAXSITE=`echo $2`
OUTLIST=`echo $3`


for outname in $OUTLIST; do
for ((siteId=1;siteId<=$MAXSITE;siteId++)); do
#for siteId in {1..$MAXSITE}; do

echo "postp.py out $plid $outname $siteId" >> tasks

done
done
