package UserModsTools;
my $pck_nm = 'UserModsTools';

# Provides tools for applying modifications via a user_mods directory, a.k.a. a
# testmods directory in the case of automated tests. These modifications can
# take the form of user_nl files, xmlchange commands (or other commands), and
# SourceMods. They are applied recursively, so that one user_mods directory can
# pull in changes from another.

# Public routines:

# apply_mods(user_mods_dir, caseroot, print_level, is_test)
#   Applies the mods in user_mods_dir (and any other included mods directories,
#   recursively) to the case given by caseroot

# Note that entities prefixed with an underscore (_) should be treated as
# private to this module, and are not meant to be accessed by other code.

use strict;
use Cwd;
use File::Spec;

# ------------------------------------------------------------------------
# Public routines
# ------------------------------------------------------------------------

sub apply_mods {
   # Apply mods from a user_mods directory. 
   #
   # Usage: apply_mods($user_mods_dir, $caseroot, $print_level, $is_test)

   my ($user_mods_dir, $caseroot, $print_level, $is_test) = @_;
   my $shell_commands_filename = 'shell_commands';

   _apply_mods_recursively($user_mods_dir, $caseroot, $print_level, $is_test,
                          $shell_commands_filename);

   _apply_shell_commands($caseroot, $shell_commands_filename, $print_level);
}

# ------------------------------------------------------------------------
# Private routines
# ------------------------------------------------------------------------

sub _apply_mods_recursively {
   # Apply mods from a user_mods directory, recursively descending into any
   # included mods directories. Note that this does not execute the xmlchange
   # commands, but instead relies on the apply_mods wrapper to do that.
   #
   # Usage: _apply_mods_recursively($user_mods_dir, $caseroot, $print_level, $is_test, $shell_commands_filename)
   #
   # shell_commands_filename gives the name of the file in caseroot that will
   # contain the xmlchange commands (and possibly other shell commands)

   my ($user_mods_dir, $caseroot, $print_level, $is_test, $shell_commands_filename) = @_;

   # cd to the mods directory; this facilitates processing the includes
   # directories recursively, allowing for relative paths
   my $cwd = getcwd();
   chdir "$user_mods_dir";

   # ------------------------------------------------------------------------
   # Process mods listed in include_user_mods file recursively
   #
   # Each line in this file gives a directory (with relative or absolute path)
   # of additional mods to apply before applying the mods in the current directory.
   # ------------------------------------------------------------------------

   my @include_dirs = _get_include_dir_list('include_user_mods');
   foreach my $include_dir (@include_dirs) {
      if (! -d $include_dir) {
         die "ERROR: Cannot find desired user_mods_dir to include: $include_dir";
      }

      _apply_mods_recursively($include_dir, $caseroot, $print_level, $is_test, 
                             $shell_commands_filename);
   }

   # ------------------------------------------------------------------------
   # Process these mods
   #
   # The processing of user_nl files and shell_commands files is done by
   # appending to the corresponding file in the case directory, so we take the
   # union of all settings. Note that, since these mods are applied after
   # everything in include_user_mods, they will override any settings set
   # through those includes, if the same namelist item, xml variable, or
   # SourceMod file was set in multiple places.
   # ------------------------------------------------------------------------
   
   _apply_mods_from_current_dir($caseroot, $print_level, $is_test, $shell_commands_filename);

   chdir "$cwd";
}   

#-------------------------------------------------------------------------------

sub _get_include_dir_list {
   # Get list of directories to include in applying user mods, by reading a file
   # that lists one directory path per line. Directory paths can be relative to
   # the current directory. Blank lines are skipped.
   #
   # Returns an empty list if $includes_file doesn't exist.
   #
   # Usage: @include_dirs = _get_include_dir_list($includes_file)
   
   my ($includes_file) = @_;
   
   my @include_dirs;
   if (-e $includes_file) {
      open(my $includes_fh, "<", $includes_file);
      while(my $include_dir = <$includes_fh>) {
         # trim trailing whitespace
         chomp $include_dir;
         $include_dir =~ s/\s+$//;

         # skip blank lines (which will include lines that originally had some
         # whitespace, because we have removed all trailing whitespace above)
         next if ($include_dir =~ m/^$/);
         
         # convert relative to absolute paths; this isn't strictly necessary,
         # but makes diagnostic messages more clear
         my $include_dir_abs = File::Spec->rel2abs($include_dir);

         push (@include_dirs, $include_dir_abs);
      }
      close($includes_fh);
   }

   return @include_dirs;

}

#-------------------------------------------------------------------------------

sub _apply_mods_from_current_dir {
   # Apply mods from the current working directory. This routine does NOT handle
   # the recursive application of mods - rather, it just does the application of
   # mods that are contained directly here. Note that it does not actually
   # execute the commands in the shell_commands file - rather, it appends
   # these commands into a file in the case directory, which can be executed later.
   #
   # Usage: _apply_mods_from_current_dir($caseroot, $print_level, $is_test, $shell_commands_filename)
   #
   # shell_commands_filename gives the name of the file in caseroot that will
   # contain the xmlchange commands (and possibly other shell commands)

   my ($caseroot, $print_level, $is_test, $shell_commands_filename) = @_;

   my $cwd = getcwd();
   print "Applying mods from: $cwd\n";

   my @user_nl_files = glob("user_nl*");

   # allow xmlchange_cmnds as well as shell_commands; the latter is
   # preferred, but the former is supported for backwards compatibility
   my @shell_commands_files = ("xmlchange_cmnds", "shell_commands");

   my @source_mods_dir = glob("SourceMods/*");

   if (! @user_nl_files && ! @shell_commands_files && ! @source_mods_dir) {
      die "ERROR: There are no user_nl files, shell_commands file, or SourceMods directories";
   } else {
      # Append the user_nl_* files from this user_mods directory to the end of
      # the corresponding files in the case directory
      foreach my $file (@user_nl_files) {
         if (-f "$file") {
            if ($print_level>=2) {
               print "Append to $file.\n";
            }
            my $append = "cat $file >> $caseroot/$file";
            system($append) == 0 or die "ERROR: $append failed: $?\n";
         }
      }

      # Append to the shell_commands file
      my $num_commands_files_found = 0;
      foreach my $file (@shell_commands_files) {
         if (-f "$file") {
            $num_commands_files_found += 1;
            # In principle, you could have both xmlchange_cmnds and shell_commands,
            # but that could lead to unexpected results, so don't allow that
            if ($num_commands_files_found > 1) {
               die "ERROR: Multiple shell_commands/xmlchange_cmnds files found in the user_mods_dir; please use only one";
            }

            if ($print_level>=2) {
               print "Append to $shell_commands_filename.\n";
            }
            my $append = "cat $file >> $caseroot/$shell_commands_filename";
            system($append) == 0 or die "ERROR: $append failed: $?\n";
         }
      }
        

      # Copy any files in SourceMods over provided the same directory exists in the case
      my $num_sourcemods_found = 0;
      foreach my $dir (@source_mods_dir) {
         if ($print_level>=2) {
            print "Copy $dir over.\n";
         }
         if ( $dir =~ m/(SourceMods\/.+)/ ) {
            my $mydir = "$caseroot/$1";
            if (-d "$dir" && -d "$mydir" ) {
               foreach my $file ( glob("$dir/*") ) {
                  if ($print_level>=2) {
                     print "Copy $file over.\n";
                  }
                  my $copy = "cp $file $mydir";
                  system($copy) == 0 or die "ERROR: $copy failed: $?\n";
                  $num_sourcemods_found += 1;
               }
            }
         }
      }

      if ($num_sourcemods_found > 0 && $is_test) {
         # Don't allow SourceMods in testmods, because (1) you shouldn't be
         # using SourceMods for automated tests, and (2) this likely won't
         # work correctly if multiple tests are sharing a single build.
         die "SourceMods are not allowed in the testmods directory for automated tests";
      }
   }
}

#-------------------------------------------------------------------------------

sub _apply_shell_commands {
   # Applies commands in a shell_commands file
   #
   # Usage: _apply_shell_commands($caseroot, $shell_commands_filename, $print_level)

   my ($caseroot, $shell_commands_filename, $print_level) = @_;

   chmod 0777, "$caseroot/$shell_commands_filename";
   my $cwd = getcwd(); # current working directory
   chdir "$caseroot";
   if ($print_level>=2) {
      print "Execute $shell_commands_filename.\n";
   }
   system ("./$shell_commands_filename");
   chdir "$cwd";
}
   


