###############################################################################
#
# Module: NMLTest::CompFiles
#
# Created by Erik Kluzek NCAR
#
# This is a tester built on top of Test::More to compare namelist files
# (or really any ASCII text files). There is a mechanism for telling the
# test object that you should (or should NOT) expect the comparison to be
# exact or not.
# 
###############################################################################

package NMLTest::CompFiles;
use strict;
use Test::More;
use IO::File;

=head1 NAME

NMLTest::CompFiles - A comparision tester for namelist (or ASCII text) files

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

sub new {
   my $self  = {};
   my $class = shift;
   my $dir   = shift;
   my @files = @_;

   my $nm = ref($self)."\:\:new";

   my %diffref        = { };
   $self->{'diffref'} = \%diffref;
   if ( ! -d "$dir" ) {
     die "ERROR::($nm) Input directory ($dir) does NOT exist!\n";
   }
   $self->{'dir'}     = $dir;
   $self->{'files'}   = \@_;
   bless( $self, $class );
}

sub checkfilesexist {
#
# Check that files exist
#
   my $self = shift;
   my $type = shift;
   my $mode = shift;
   my $nm = ref($self)."\:\:checkfilesexist";

   my $filesref = $self->{'files'};
   my $confdir  = $self->{'dir'};
   foreach my $file ( @$filesref ) {
      my $exists = ( -f "$confdir/$file" );
      ok( $exists, "$file file exists" );
      if ( $exists ) {
         $self->dodiffonfile(      $file, $type, $mode );
      } else {
         $self->doNOTdodiffonfile( $file, $type, $mode );
      }
   }
}

sub comparefiles {
#
# Compare the resultant files to the default versions
#
   my $self      = shift;
   my $type      = shift;
   my $comp_mode = shift;
   my $compdir   = shift;
   my $nm = ref($self)."\:\:comparefiles";

   $type      =~ s/ /+/g;
   $comp_mode =~ s/ /+/g;
   my $confdir  = $self->{'dir'};
   my $diffref   = $self->{'diffref'};
   if ( ! defined($type)    ) {
      $type = "default";
   }
   my $compare = "compare to previous tag";
   if ( ! defined($compdir) ) {
      $compdir = ".";
      $compare = undef;
   }
   if ( ! -d "$compdir" ) {
      die "ERROR($nm):: Compare directory $compdir does NOT exist!\n";
   }
   print "Compare files for $type type MODE=$comp_mode $compare\n";
   my $diffstat;
   my %diffhas = %$diffref;
   my $same = "file the same as expected";
   my $diff = "file different as expected";
   my $filesref = $self->{'files'};
   foreach my $file ( @$filesref ) {
      if ( ! -f "$compdir/${file}.$comp_mode.${type}" ) {
         print "WARNING($nm):: File $compdir/${file}.$comp_mode.${type} does NOT exist!\n";
         fail( "compare file $file DNE for $comp_mode and $type" );
      } else {
         if ( ! exists($diffhas{$comp_mode}{$type}{$file}) ) {
            die "ERROR($nm):: difference is NOT setup for $comp_mode ${type} $file!\n";
         }
         system( "diff $confdir/${file}     $compdir/${file}.$comp_mode.${type} > /dev/null" );
         $diffstat = $?;
         if ( $diffhas{$comp_mode}{$type}{$file} ) {
            ok( ! $diffstat, "$file $same for $comp_mode" );
         } else {
            ok( $diffstat,   "$file different as expected for $comp_mode" );
         }
      }
   }

}

sub copyfiles {
#
# Copy the namelist files to default names for comparisions
#
   my $self    = shift;
   my $type    = shift;
   my $mode    = shift;
   my $nm = ref($self)."\:\:copyfiles";

   $type =~ s/ /+/g;
   $mode =~ s/ /+/g;
   my $diffref   = $self->{'diffref'};
   my $filesref = $self->{'files'};
   my $confdir  = $self->{'dir'};
   foreach my $file ( @$filesref ) {
      system( "/bin/cp  $confdir/$file     ${file}.${mode}.${type}"     );
      $$diffref{${mode}}{${type}}{$file} = 1;
   }
   print "$type namelists for $mode\n";
   foreach my $file ( @$filesref ) {
      system( "/bin/cat $file.${mode}.${type}"     );
   }
}


sub shownmldiff {
#
# Show the differences in the namelists
#
   my $self      = shift;
   my $type      = shift;
   my $comp_mode = shift;
   my $nm = ref($self)."\:\:shownmldiff";

   $type      =~ s/ /+/g;
   $comp_mode =~ s/ /+/g;
   my $filesref = $self->{'files'};
   my $confdir  = $self->{'dir'};
   foreach my $file ( @$filesref ) {
      my $file1 = "$confdir/$file";
      if ( ! -f "$file1" ) {
         print "$file1 does NOT exist\n";
         return;
      }
      my $file2 = "${file}.${comp_mode}.${type}";
      if ( ! -f "$file2" ) {
         print "$file2 does NOT exist\n";
         return;
      }
      print "Diff in in $file to $type $comp_mode version";
      system( "diff $file1 $file2" );
   }

}



sub dodiffonfile {
#
# Set it so that it does do a difference on the given input file
#
  my $self    = shift;
  my $file    = shift;
  my $type    = shift;
  my $mode    = shift;
  my $nm = ref($self)."\:\:dodiffonfile";

  $type =~ s/ /+/g;
  $mode =~ s/ /+/g;
  my $diffref   = $self->{'diffref'};
  if ( ! defined($type) ) {
     $type = "default";
  }
  $$diffref{$mode}{$type}{$file} = 1;
}


sub doNOTdodiffonfile {
#
# Set it so that it does NOT do a difference on the given input file
#
  my $self    = shift;
  my $file    = shift;
  my $type    = shift;
  my $mode    = shift;
  my $nm = ref($self)."\:\:doNOTdodiffonfile";

  $type =~ s/ /+/g;
  $mode =~ s/ /+/g;
  my $diffref   = $self->{'diffref'};
  if ( ! defined($type) ) {
     $type = "default";
  }
  $$diffref{$mode}{$type}{$file} = 0;
}

#-----------------------------------------------------------------------------------------------

1   # to make use or require happy
