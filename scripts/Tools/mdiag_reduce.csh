#!/bin/csh
mdiag -n --xml | perl -pe 's/\<\/node\>/\<\/node\>\n/g' | perl mdiag_reduce.pl

# perl mdiag_reduce.pl
#  # define variables
#  my @mdiag_line;
#  my $mdiag_field;
#  my $job_id;
#  my $nid_id;
#  my $nid_state;
#
#  # extract useful information
#  while ( <STDIN> ){
#    @mdiag_line = split(/ /,$_);
#    $job_id = -1;
#    $nid_id = -1;
#    foreach $mdiag_field (@mdiag_line){
#      if (m/JOBLIST=\"(\d+)\"/){
#        $job_id = $1;
#      }
#      if (m/NODEID=\"(\d+)\"/){
#        $nid_id = $1;
#      }
#      if (m/NODESTATE=\"(\S+)\"/){
#        $nid_state = $1;
#      }
#    }
#    if ($nid_id != -1){
#      print "$nid_id $job_id $nid_state\n";
#    }
#  }

