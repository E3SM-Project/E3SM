#!/usr/local/bin/perl

#  Perl 5 required.
require 5.000;

#  Strict syntax checking.
use strict;

#  Modules.
use Cwd;

=head1 NAME

mdiag_reduce

=head1 SYNOPSIS

extracts node id, node status, and application assignment
for output of 'mdiag -n --xml' after it has been split into
separate 'node' records

=head1 DESCRIPTION

read from standard input the result of applying
''s/\<\/node\>/\<\/node\>\n/g' applied to the output of
'mdiag -n --xml'

=head1 AUTHOR

Pat Worley E<lt>F<worleyph@ornl.gov>E<gt>

=cut

###############################################################################
#
#  Miscellaneous global flags and variables.
#
###############################################################################

###############################################################################
#
#  This is the call to main.
#
###############################################################################

&main();

exit 0;

###############################################################################
#
#  This is main.
#
###############################################################################

sub main {

  # define variables
  my @mdiag_line;
  my $mdiag_field;
  my $job_id;
  my $nid_id;
  my $nid_state;
  my $nid_features;

  # extract useful information
  while ( <STDIN> ){
    @mdiag_line = split(/ /,$_);
    $job_id = -1;
    $nid_id = -1;
    foreach $mdiag_field (@mdiag_line){
      if (m/JOBLIST=\"(\d+)\"/){
        $job_id = $1;
      }
      if (m/NODEID=\"(\d+)\"/){
        $nid_id = $1;
      }
      if (m/NODESTATE=\"(\S+)\"/){
        $nid_state = $1;
      }
      if (m/FEATURES=\"(\S+)\"/){
        $nid_features = $1;
      }
    }
    if ($nid_id != -1){
      print "$nid_id $job_id $nid_state $nid_features\n";
    }
  }

  exit ;

}

1;

__END__
