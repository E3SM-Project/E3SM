package ProjectTools;
my $pkg_nm = 'ProjectTools';

# Provides tools for obtaining and working with the PROJECT variable, which
# determines the project / account number to use in job submission scripts
# and/or directory paths on some machines.

# Public routines:

# find_project()
#   Tries to find a project / account number to use

# check_project_required_but_unset(project, cfg_ref)
#   Checks if PROJECT is required for this machine, but has not been set

# set_project(project, cfg_ref)
#   Sets the PROJECT xml variable in cfg_ref

use strict;
require ConfigCase;

use constant _PROJECT_UNSET => "";

# ------------------------------------------------------------------------
# Public routines
# ------------------------------------------------------------------------

sub find_project
{
    # Tries to find a project/account number to use. If none can be found,
    # returns _PROJECT_UNSET.

   my $project = _PROJECT_UNSET;

   if (defined $ENV{'PROJECT'}) {
      $project = $ENV{'PROJECT'};
   }
   elsif (defined $ENV{'ACCOUNT'}) {
      $project = $ENV{'ACCOUNT'};
   }
   else {
      # Loop over possible files that can contain a project number.
      # (We should eventually remove .ccsm_proj, but it is kept now for
      # backwards compatibility.)
      my @proj_filenames = ($ENV{'HOME'} . "/.cesm_proj", 
                            $ENV{'HOME'} . "/.ccsm_proj");

      foreach my $proj_filename (@proj_filenames) {
         if (-f $proj_filename) {
            $project = _read_project_from_file($proj_filename);
            last if ($project ne _PROJECT_UNSET);
         }
      }
   }

   return $project;
}


sub check_project_required_but_unset
{
   # Checks if PROJECT is required for this machine, but has not been set

   my ($project, $cfg_ref) = @_;

   if (($project eq _PROJECT_UNSET) && ($cfg_ref->get('PROJECT_REQUIRED') eq "TRUE")) {
      die <<"EOF";
** ERROR **
   A project must be specified for this machine, in order to set paths and/or an
   account number in job submission scripts.

   You can specify a project number in one of the following ways:
   (1) via the -project argument to create_newcase
   (2) via the \$PROJECT or \$ACCOUNT environment variables
   (3) via a file named .cesm_proj in your home directory, whose first line gives a project number
EOF
   }

}


sub set_project
{
   # Sets the PROJECT xml variable in cfg_ref. However, if $project hasn't been
   # set, then this does NOT set the xml variable, instead keeping it at its
   # default value.

   my ($project, $cfg_ref) = @_;

   if ($project ne _PROJECT_UNSET) {
      $cfg_ref->set('PROJECT', "$project");
   }
}

# ------------------------------------------------------------------------
# Private routines
# ------------------------------------------------------------------------

sub _read_project_from_file 
{
   # Try to read a project number from a file with the given file name. Return
   # the project number if found. If we can't find a project number in that
   # file, return _PROJECT_UNSET.
   
   my ($proj_filename) = @_;
   
   my $project = _PROJECT_UNSET;

   # read first line
   open my $proj_file, '<', $proj_filename;
   my $firstline = <$proj_file>;
   close $proj_file;

   # check the first line for something that looks like a project number
   if ($firstline =~ /^\s*(\S+)/) {
      $project = $1;
   } else {
      print "WARNING: $proj_filename found, but I cannot find a project number on its first line\n";
   }

   return $project;
}
   


