package File::DirCompare;

use 5.005;
use strict;

use File::Basename;
use File::Spec::Functions;
use File::Compare ();
use File::Glob qw(bsd_glob);
use Carp;

use vars qw($VERSION);

$VERSION = '0.7';

# ----------------------------------------------------------------------------
# Private methods

sub _dir_compare
{
  my $self = shift;
  my ($dir1, $dir2, $sub, $opts) = @_;

  # Glob $dir1 and $dir2
  my (%d1, %d2);
  $d1{basename $_} = 1 foreach bsd_glob(catfile($dir1, ".*"));
  $d1{basename $_} = 1 foreach bsd_glob(catfile($dir1, "*"));
  $d2{basename $_} = 1 foreach bsd_glob(catfile($dir2, ".*"));
  $d2{basename $_} = 1 foreach bsd_glob(catfile($dir2, "*"));

  # Prune dot dirs
  delete $d1{curdir()} if $d1{curdir()};
  delete $d1{updir()}  if $d1{updir()};
  delete $d2{curdir()} if $d2{curdir()};
  delete $d2{updir()}  if $d2{updir()};

  # Setup cmp and matches subs
  my $cmp = $opts->{cmp} && ref $opts->{cmp} eq 'CODE' ? $opts->{cmp} : \&File::Compare::compare;
  my $matches = $opts->{matches} if $opts->{matches} && ref $opts->{matches} eq 'CODE';

  # Iterate over sorted and uniquified file list
  my %u;
  for my $f (map { $u{$_}++ == 0 ? $_ : () } sort(keys(%d1), keys(%d2))) {
    my $f1 = catfile($dir1, $f);
    my $f2 = catfile($dir2, $f);
    # Only in $dir1
    if (! $d2{$f}) {
      $sub->($f1, undef) unless $opts->{ignore_unique};
    }
    # Only in $dir2
    elsif (! $d1{$f}) {
      $sub->(undef, $f2) unless $opts->{ignore_unique};
    }
    # Item exists in both directories
    else {
      # Both symlinks
      if (-l $f1 && -l $f2) {
        my $t1 = readlink $f1 or croak "Cannot read symlink $f1: $!";
        my $t2 = readlink $f2 or croak "Cannot read symlink $f2: $!";
        $sub->($f1, $f2) if $t1 ne $t2;
      }
      # One symlink (i.e. different)
      elsif (-l $f1 || -l $f2) {
        $sub->($f1, $f2);
      }
      # Both directories
      elsif (-d $f1 && -d $f2) {
        $self->_dir_compare($f1, $f2, $sub, $opts);
      }
      # One directory (i.e. different)
      elsif (-d $f1 || -d $f2) {
        $sub->($f1, $f2);
      }
      # Both files - check if different
      else {
        if ($opts->{ignore_cmp}) {
          $sub->($f1, $f2);
        }
        elsif ($cmp->($f1, $f2) != 0) {
          $sub->($f1, $f2);
        }
        elsif ($matches) {
          $matches->($f1, $f2);
        }
      }
    }
  }
}

# ----------------------------------------------------------------------------
# Public methods

sub compare
{
  my $self = shift;
  my ($dir1, $dir2, $sub, $opts) = @_;

  croak "Not a directory: '$dir1'" unless -d $dir1;
  croak "Not a directory: '$dir2'" unless -d $dir2;
  croak "Not a subroutine: '$sub'" unless ref $sub eq 'CODE';
  croak "Not a hashref: '$opts'" if $opts && ref $opts ne 'HASH';

  $self = $self->new unless ref $self;
  $self->_dir_compare(@_);
}

# ----------------------------------------------------------------------------
# Constructors

sub new { bless {}, shift }

# ----------------------------------------------------------------------------

1;

__END__

=head1 NAME

File::DirCompare - Perl module to compare two directories using
callbacks.


=head1 SYNOPSIS

  use File::DirCompare;

  # Simple diff -r --brief replacement
  use File::Basename;
  File::DirCompare->compare($dir1, $dir2, sub {
    my ($a, $b) = @_;
    if (! $b) {
      printf "Only in %s: %s\n", dirname($a), basename($a);
    } elsif (! $a) {
      printf "Only in %s: %s\n", dirname($b), basename($b);
    } else {
      print "Files $a and $b differ\n";
    }
  });

  # Version-control like Deleted/Added/Modified listing
  my (@listing, @modified);     # use closure to collect results
  File::DirCompare->compare('old_tree', 'new_tree', sub {
    my ($a, $b) = @_;
    if (! $b) {
      push @listing, "D   $a";
    } elsif (! $a) {
      push @listing, "A   $b";
    } else {
      if (-f $a && -f $b) {
        push @listing, "M   $b";
        push @modified, $b;
      } else {
        # One file, one directory - treat as delete + add
        push @listing, "D   $a";
        push @listing, "A   $b";
      }
    }
  });


=head1 DESCRIPTION

File::DirCompare is a perl module to compare two directories using
a callback, invoked for all files that are 'different' between the
two directories, and for any files that exist only in one or other
directory ('unique' files).

File::DirCompare has a single public compare() method, with the
following signature:

  File::DirCompare->compare($dir1, $dir2, $sub, $opts);

The first three arguments are required - $dir1 and $dir2 are paths
to the two directories to be compared, and $sub is the subroutine
reference called for all unique or different files. $opts is an
optional hashref of options - see L<OPTIONS> below.

The provided subroutine is called for all unique files, and for
every pair of 'different' files encountered, with the following
signature:

  $sub->($file1, $file2)

where $file1 and $file2 are the paths to the two files. For 'unique'
files i.e. where a file exists in only one directory, the subroutine
is called with the other argument 'undef' i.e. for:

  $sub->($file1, undef)
  $sub->(undef, $file2)

the first indicates $file1 exists only in the first directory given
($dir1), and the second indicates $file2 exists only in the second
directory given ($dir2).

=head2 OPTIONS

The following optional arguments are supported, passed in using a
hash reference after the three required arguments to compare() e.g.

  File::DirCompare->compare($dir1, $dir2, $sub, {
    cmp             => $cmp_sub,
    ignore_cmp      => 1,
    ignore_unique   => 1,
    matches         => $matches_sub,
  });

=over 4

=item cmp

By default, two files are regarded as different if their contents do
not match (tested with File::Compare::compare). That default behaviour
can be overridden by providing a 'cmp' subroutine to do the file
comparison, returning zero if the two files are equal, and non-zero
if not.

E.g. to compare using modification times instead of file contents:

  File::DirCompare->compare($dir1, $dir2, $sub, {
    cmp => sub { -M $_[0] <=> -M $_[1] },
  });

=item ignore_cmp

If you want to see I<all> corresponding files, not just 'different'
ones, set the 'ignore_cmp' flag to tell File::DirCompare to skip its
file comparison checks i.e.

  File::DirCompare->compare($dir1, $dir2, $sub,
    { ignore_cmp => 1 });

=item ignore_unique

If you want to ignore files that only exist in one of the two
directories, set the 'ignore_unique' flag i.e.

  File::DirCompare->compare($dir1, $dir2, $sub,
    { ignore_unique => 1 });

=item matches

Subroutine to be called for file pairs that I<match>, with the
following signature:

  $sub->($file1, $file2)

These pairs are ordinarily ignored (unless C<ignore_cmp> is set).

=back

=head1 SEE ALSO

File::Dircmp, which provides similar functionality (and whose
directory walking code I've adapted for this module), but a simpler
reporting-only interface, something like the first example in the
SYNOPSIS above.

=head1 AUTHOR AND CREDITS

Gavin Carr <gavin@openfusion.com.au>

Thanks to Robin Barker for a bug report and fix for glob problems
with whitespace.

=head1 COPYRIGHT AND LICENSE

Copyright 2006-2012 by Gavin Carr E<lt>gavin@openfusion.com.auE<gt>.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
