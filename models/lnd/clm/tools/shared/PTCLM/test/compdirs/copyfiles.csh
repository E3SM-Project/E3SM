#!/bin/csh -f
#
# This script is for initially populating a list of new sites
#
set sites = ( TS-Ts1/ TS-Ts2/ TS-Ts3/ US-CHATS/ CA-Let/ CA-Man/ US-WCr/ CA-Ca1/ CA-Obs/ CA-Ojp/ US-Dk2/ US-Dk3/ US-Me4/ CA-Qfo/ US-Me/ US-MOz/ BE-Vie/ BR-Sa3/ DE_Tha/ ES-ES1/ FL-Hyy/ FL-Kaa/ IT-Col/ IT-Cpz/ US-FPe/ US-NR1/ US-0Brw/ US-ARM/ US-Var/ US-Bo1/ US-Ho1/ US-MMS/ BR-Sa1/ CA-Oas/ US-Ne3/ ) 
foreach site( $sites )
   echo $site
   set dir = `ls -1d *$site`
   echo "dir = $dir"
   set copydir=`ls -qd ../testing_dir/*/1x1pt_$site`
   echo "copydir = $copydir"
   \cp $copydir/README.PTCLM    $dir
   \cp $copydir/*.log           $dir
   rm $copydir/surf*.log
   \cp $copydir/user_nl_clm     $dir
   \cp $copydir/xmlchange_cmnds $dir
end
