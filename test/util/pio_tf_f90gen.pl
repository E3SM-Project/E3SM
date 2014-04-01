#!/usr/bin/env perl

use strict;
use warnings;
use FileHandle;
use File::stat;
use File::Basename;
use Getopt::Long;

my($verbose);
my($nargs);
my($template_fname, $output_fname);
# If the template source has PIO_TF_AUTO_TEST* and the test
# does not have a test driver we need to create it
# Auto test functions do not have any arguments so our life is easy
my($template_has_auto_funcs, $template_has_test_driver);
# Indicates whether auto funcs have already been inserted
my($template_auto_funcs_inserted);
my($annotate_source);
# The current test case number is modified by multiple funcs
# * its logically easier to keep it a global
my($cur_test_case_num);

# Initializing global vars
$verbose = 0;
$nargs = 0;
$template_fname = '';
$output_fname = '';
$template_has_auto_funcs = 0;
$template_auto_funcs_inserted = 0;
$template_has_test_driver = 0;
$annotate_source = 0;
$cur_test_case_num = 1;

# This function transforms the template source one line at a time
# * Does the transformation of the template source
sub transform_src
{
  my($ref_auto_funcs_list, $in_line, $template_fname, $template_line_no, $ref_modif) = @_;
  my($out_line);
  if($verbose){ print "Transforming: $in_line\n"; }
  $_ = $in_line;
  # By default assume that source is transformed/modified
  # reset otherwise
  $$ref_modif = 1;
  if(/^(\s*)(PIO_TF_TEST_SUB_BEGIN)(.*)$/s){
    $out_line = $1 . "SUBROUTINE" . $3 . "\n";
    $out_line = $out_line . $1 . "#ifndef NO_MPIMOD\n";
    $out_line = $out_line . $1 . "  USE mpi\n";
    $out_line = $out_line . $1 . "#else\n";
    $out_line = $out_line . $1 . "  include 'mpif.h'\n";
    $out_line = $out_line . $1 . "#endif\n";
    $out_line = $out_line . $1 . "  USE pio_tutil";
  }
  elsif(/^(\s*)(PIO_TF_TEST_SUB_END)(.*)$/s){
    $out_line = $1 . "END SUBROUTINE" . $3;
  }
  elsif(/^(\s*)(PIO_TF_AUTO_TEST_SUB_BEGIN)(\s*)(.*)$/s){
    $out_line = $1 . "SUBROUTINE" . $3 . $4 . "\n";
    $out_line = $out_line . $1 . "#ifndef NO_MPIMOD\n";
    $out_line = $out_line . $1 . "  USE mpi\n";
    $out_line = $out_line . $1 . "#else\n";
    $out_line = $out_line . $1 . "  include 'mpif.h'\n";
    $out_line = $out_line . $1 . "#endif\n";
    $out_line = $out_line . $1 . "  USE pio_tutil";
    $template_has_auto_funcs = 1;

    # Add auto test func into the auto func list
    push @{$ref_auto_funcs_list}, $4;
  }
  elsif(/^(\s*)(PIO_TF_AUTO_TEST_SUB_END)(.*)$/s){
    $out_line = $1 . "END SUBROUTINE" . $3;
  }
  elsif(/^(\s*)(PIO_TF_TEST_DRIVER_BEGIN)(.*)$/s){
    $out_line = $1 . "SUBROUTINE PIO_TF_Test_driver_" . $3 . "\n";
    $out_line = $out_line . $1 . "#ifndef NO_MPIMOD\n";
    $out_line = $out_line . $1 . "  USE mpi\n";
    $out_line = $out_line . $1 . "#else\n";
    $out_line = $out_line . $1 . "  include 'mpif.h'\n";
    $out_line = $out_line . $1 . "#endif\n";
    $out_line = $out_line . $1 . "  USE pio_tutil";
    $template_has_test_driver = 1;
  }
  elsif(/^(\s*)(PIO_TF_TEST_DRIVER_END)(.*)$/s){
    # Make sure to add any auto test functions before the end
    $out_line = $1 . "\n";
    if($template_auto_funcs_inserted == 0){
    foreach(@{$ref_auto_funcs_list}){
      $out_line = $out_line . $1 . "  pio_tf_retval_utest_ = 0\n";
      $out_line = $out_line . $1 . "  IF (pio_tf_world_rank_ == 0) THEN\n";
      $out_line = $out_line . $1 . "    PRINT *, \"PIO_TF: Starting $_\"\n";
      $out_line = $out_line . $1 . "  END IF\n";
      $out_line = $out_line . $1 . "  CALL $_()\n";
      $out_line = $out_line . $1 . "  IF (pio_tf_retval_utest_ /= 0) THEN\n";
      $out_line = $out_line . $1 . "    pio_tf_nerrs_total_ = pio_tf_nerrs_total_ + 1\n";
      $out_line = $out_line . $1 . "  END IF\n";
      $out_line = $out_line . $1 . "  IF (pio_tf_world_rank_ == 0) THEN\n";
      $out_line = $out_line . $1 . "    IF (pio_tf_retval_utest_ == 0) THEN\n";
      $out_line = $out_line . $1 . "      WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF:&\n";
      $out_line = $out_line . $1 . "        Test $cur_test_case_num:\",&\n";
      $out_line = $out_line . $1 . "        \"$_\",&\n";
      $out_line = $out_line . $1 . "        \"---------\",\"PASSED\"\n";
      $out_line = $out_line . $1 . "    ELSE\n";
      $out_line = $out_line . $1 . "      WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF: &\n";
      $out_line = $out_line . $1 . "        Test $cur_test_case_num:\",&\n";
      $out_line = $out_line . $1 . "        \"$_\",&\n";
      $out_line = $out_line . $1 . "        \"---------\",\"FAILED\"\n";
      $out_line = $out_line . $1 . "    END IF\n";
      $out_line = $out_line . $1 . "  END IF\n";
      $cur_test_case_num += 1;
    }
    }
    $out_line = $out_line . $1 . "END SUBROUTINE PIO_TF_Test_driver_" . $3;
  }
  elsif(/^(\s*)(PIO_TF_AUTO_TESTS_RUN)\((.*)\)\s*$/s){
    # Make sure to add any auto test functions before the end
    $out_line = $1 . "\n";
    $out_line = $out_line . $1 . "  IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . $1 . "    PRINT *, \"PIO_TF: Running AUTO tests: \", $3\n";
    $out_line = $out_line . $1 . "  END IF\n";
    foreach(@{$ref_auto_funcs_list}){
      $out_line = $out_line . $1 . "  pio_tf_retval_utest_ = 0\n";
      $out_line = $out_line . $1 . "  pio_tf_tmp_log_str_=\"$_\"//\"(\"//$3//\")\"\n";
      $out_line = $out_line . $1 . "  IF (pio_tf_world_rank_ == 0) THEN\n";
      $out_line = $out_line . $1 . "    PRINT *, \"PIO_TF: Starting $_\"\n";
      $out_line = $out_line . $1 . "  END IF\n";
      $out_line = $out_line . $1 . "  CALL $_()\n";
      $out_line = $out_line . $1 . "  IF (pio_tf_retval_utest_ /= 0) THEN\n";
      $out_line = $out_line . $1 . "    pio_tf_nerrs_total_ = pio_tf_nerrs_total_ + 1\n";
      $out_line = $out_line . $1 . "  END IF\n";
      $out_line = $out_line . $1 . "  IF (pio_tf_world_rank_ == 0) THEN\n";
      $out_line = $out_line . $1 . "    IF (pio_tf_retval_utest_ == 0) THEN\n";
      $out_line = $out_line . $1 . "      WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF:&\n";
      $out_line = $out_line . $1 . "        Test $cur_test_case_num:\",&\n";
      $out_line = $out_line . $1 . "        pio_tf_tmp_log_str_,&\n";
      $out_line = $out_line . $1 . "        \"---------\",\"PASSED\"\n";
      $out_line = $out_line . $1 . "    ELSE\n";
      $out_line = $out_line . $1 . "      WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF: &\n";
      $out_line = $out_line . $1 . "        Test $cur_test_case_num:\",&\n";
      $out_line = $out_line . $1 . "        pio_tf_tmp_log_str_,&\n";
      $out_line = $out_line . $1 . "        \"---------\",\"FAILED\"\n";
      $out_line = $out_line . $1 . "    END IF\n";
      $out_line = $out_line . $1 . "  END IF\n";
      $cur_test_case_num += 1;
    }
    $template_auto_funcs_inserted = 1;
  }
  elsif(/^(\s*)PIO_TF_TEST_RUN\(\s*([^(]+)\((.*)\),(.*)\)\s*$/s){
    my ($test_name) = $2;
    my ($test_args) = $3;
    # $4 is $comment

    $out_line = $1 . "\n";
    $out_line = $out_line . $1 . "pio_tf_retval_utest_ = 0\n";
    $out_line = $out_line . $1 . "pio_tf_tmp_log_str_=\"$test_name\"//\"(\"//$4//\")\"\n";
    $out_line = $out_line . $1 . "IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . $1 . "  PRINT *, \"PIO_TF: Starting \",&\n";
    $out_line = $out_line . $1 . "    pio_tf_tmp_log_str_\n";
    $out_line = $out_line . $1 . "END IF\n";
    $out_line = $out_line . $1 . "CALL $test_name($test_args)\n";
    $out_line = $out_line . $1 . "IF (pio_tf_retval_utest_ /= 0) THEN\n";
    $out_line = $out_line . $1 . "  pio_tf_nerrs_total_ = pio_tf_nerrs_total_ + 1\n";
    $out_line = $out_line . $1 . "END IF\n";
    $out_line = $out_line . $1 . "IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . $1 . "  IF (pio_tf_retval_utest_ == 0) THEN\n";
    $out_line = $out_line . $1 . "    WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF: &\n";
    $out_line = $out_line . $1 . "      Test $cur_test_case_num:\",&\n";
    $out_line = $out_line . $1 . "      pio_tf_tmp_log_str_,&\n";
    $out_line = $out_line . $1 . "      \"--------\", \"PASSED\"\n";
    $out_line = $out_line . $1 . "  ELSE\n";
    $out_line = $out_line . $1 . "    WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF: &\n";
    $out_line = $out_line . $1 . "      Test $cur_test_case_num:\",&\n";
    $out_line = $out_line . $1 . "      pio_tf_tmp_log_str_,&\n";
    $out_line = $out_line . $1 . "      \"--------\", \"FAILED\"\n";
    $out_line = $out_line . $1 . "  END IF\n";
    $out_line = $out_line . $1 . "END IF";
    $cur_test_case_num += 1;
}
  elsif(/^(\s*)PIO_TF_PASSERT\((.+),([^)]+)\)(.*)$/s){
    $out_line = $1 . $4 . "\n";
    $out_line = $out_line . $1 . "IF (.NOT. (PIO_TF_Passert_($2))) THEN\n";
    $out_line = $out_line . $1 . "  pio_tf_retval_utest_ = -1\n";
    $out_line = $out_line . $1 . "  IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . $1 . "    PRINT *, \"PIO_TF: Assertion failed :\",&\n";
    $out_line = $out_line . $1 . "      " . $3 . ",&\n";
    $out_line = $out_line . $1 . "       \":\", __FILE__, \":\", __LINE__,&\n";
    $out_line = $out_line . $1 . "      " . "\"($template_fname:$template_line_no)\"\n";
    $out_line = $out_line . $1 . "  END IF\n";
    $out_line = $out_line . $1 . "  RETURN\n";
    $out_line = $out_line . $1 . "END IF";
  }
  elsif(/^(\s*)PIO_TF_ASSERT\((.+),([^)]+)\)(.*)$/s){
    $out_line = $1 . $4 . "\n";
    $out_line = $out_line . $1 . "IF (.NOT. ($2)) THEN\n";
    $out_line = $out_line . $1 . "  pio_tf_retval_utest_ = -1\n";
    $out_line = $out_line . $1 . "  PRINT *, \"PIO_TF: Assertion failed :\",&\n";
    $out_line = $out_line . $1 . "    " . $3 . ",&\n";
    $out_line = $out_line . $1 . "     \":\", __FILE__, \":\", __LINE__,&\n";
    $out_line = $out_line . $1 . "    " . "\"($template_fname:$template_line_no)\"\n";
    $out_line = $out_line . $1 . "  RETURN\n";
    $out_line = $out_line . $1 . "END IF";
  }
  elsif(/^(\s*)PIO_TF_CHECK_ERR\(([^,]+),(.+)\)(\s*)$/s){
    $out_line = $1 . $4 . "\n";
    $out_line = $out_line . $1 . "IF (.NOT. (PIO_TF_Passert_(($2) == PIO_NOERR))) THEN\n";
    $out_line = $out_line . $1 . "  pio_tf_retval_utest_ = -1\n";
    $out_line = $out_line . $1 . "  IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . $1 . "    PRINT *, \"PIO_TF: PIO Function failed:\",&\n";
    $out_line = $out_line . $1 . "      " . $3 . ",&\n";
    $out_line = $out_line . $1 . "      \":\", __FILE__, \":\", __LINE__,&\n";
    $out_line = $out_line . $1 . "      " . "\"($template_fname:$template_line_no)\"\n";
    $out_line = $out_line . $1 . "  END IF\n";
    $out_line = $out_line . $1 . "  RETURN\n";
    $out_line = $out_line . $1 . "END IF";
  }
  elsif(/^(\s*)PIO_TF_ERROR\(([^)]+)\)(.*)$/s){
    $out_line = $1 . $3 . "\n";
    $out_line = $out_line . $1 . "pio_tf_retval_utest_ = -1\n";
    $out_line = $out_line . $1 . "IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . $1 . "  PRINT *, \"PIO_TF: Fatal Error:\",&\n";
    $out_line = $out_line . $1 . "    " . $2 . ",&\n";
    $out_line = $out_line . $1 . "     \":\", __FILE__, \":\", __LINE__,&\n";
    $out_line = $out_line . $1 . "    " . "\"($template_fname:$template_line_no)\"\n";
    $out_line = $out_line . $1 . "END IF\n";
    $out_line = $out_line . $1 . "RETURN";
  }
  elsif(/^(\s*)PIO_TF_LOG\(([0-9]+),([^)]+)\)(.*)$/s){
    $out_line = "\n";
    $out_line = $out_line . $1 . "IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . $1 . "  IF (pio_tf_log_level_ >= " . $2 . ") THEN\n";
    $out_line = $out_line . $1 . "    WRITE(*,\"(A)\",ADVANCE=\"NO\") \"PIO_TF: \"\n";
    $out_line = $out_line . $1 . "    WRITE(*," . $3 . ") " . $4 . "\n";
    $out_line = $out_line . $1 . "  END IF\n";
    $out_line = $out_line . $1 . "END IF";
  }
  else{
    $out_line = $_;
    $$ref_modif = 0;
  }
  return $out_line;
}

# Returns the default test main code
sub get_default_test_main
{
  my($out_line);
  $out_line = "\n\n";
  $out_line = $out_line . "  PROGRAM PIO_TF_Test_main_\n";
  $out_line = $out_line . "#ifndef NO_MPIMOD\n";
  $out_line = $out_line . "    USE mpi\n";
  $out_line = $out_line . "#else\n";
  $out_line = $out_line . "    include 'mpif.h'\n";
  $out_line = $out_line . "#endif\n";
  $out_line = $out_line . "    USE pio_tutil\n";
  $out_line = $out_line . "    IMPLICIT NONE\n";
  $out_line = $out_line . "    INTEGER ierr\n";
  $out_line = $out_line . "\n";
  $out_line = $out_line . "    pio_tf_nerrs_total_=0\n";
  $out_line = $out_line . "    pio_tf_retval_utest_=0\n";
  $out_line = $out_line . "    CALL MPI_Init(ierr)\n";
  $out_line = $out_line . "    CALL PIO_TF_Init_()\n";
  $out_line = $out_line . "    CALL PIO_TF_Test_driver_()\n";
  $out_line = $out_line . "    CALL PIO_TF_Finalize_()\n";
  $out_line = $out_line . "    IF (pio_tf_world_rank_ == 0) THEN\n";
  $out_line = $out_line . "      IF (pio_tf_nerrs_total_ == 0) THEN\n";
  $out_line = $out_line . "        IF (pio_tf_retval_utest_ == 0) THEN\n";
  $out_line = $out_line . "          WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF: \",&\n";
  $out_line = $out_line . "           \"All tests\", \"---------\", \"PASSED\"\n";
  $out_line = $out_line . "        ELSE\n";
  $out_line = $out_line . "          WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF: \",&\n";
  $out_line = $out_line . "           \"Test driver\", \"---------\", \"FAILED\"\n";
  $out_line = $out_line . "        END IF\n";
  $out_line = $out_line . "      ELSE\n";
  $out_line = $out_line . "        WRITE(*,PIO_TF_TEST_RES_FMT2) \"PIO_TF:[\",&\n";
  $out_line = $out_line . "          pio_tf_nerrs_total_,\"] Tests\",&\n";
  $out_line = $out_line . "          \"----- FAILED\"\n";
  $out_line = $out_line . "      END IF\n";
  $out_line = $out_line . "    END IF\n";
  $out_line = $out_line . "    CALL MPI_Finalize(ierr)\n";
  $out_line = $out_line . "    IF (pio_tf_nerrs_total_ /= 0) THEN\n";
  $out_line = $out_line . "      STOP PIO_TF_ERR\n";
  $out_line = $out_line . "    END IF\n";
  $out_line = $out_line . "  END PROGRAM\n";

  return $out_line;
}

# Returns the default test driver
sub get_default_test_driver
{
  my($ref_auto_funcs_list) = @_;
  my($out_line);

  $out_line = "\n\n";
  $out_line = $out_line . "  SUBROUTINE PIO_TF_Test_driver_\n";
  $out_line = $out_line . "#ifndef NO_MPIMOD\n";
  $out_line = $out_line . "    USE mpi\n";
  $out_line = $out_line . "#else\n";
  $out_line = $out_line . "    include 'mpif.h'\n";
  $out_line = $out_line . "#endif\n";
  $out_line = $out_line . "    USE pio_tutil\n";
  $out_line = $out_line . "    IMPLICIT NONE\n";
  if($template_auto_funcs_inserted == 1){
    print "Error parsing template file, auto tests can only be inserted (PIO_TF_AUTO_TESTS_RUN) in a test driver\n";
    exit;
  }
  foreach(@{$ref_auto_funcs_list}){
    $out_line = $out_line . "    pio_tf_retval_utest_ = 0\n";
    $out_line = $out_line . "    IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . "      PRINT *, \"PIO_TF: Starting $_\"\n";
    $out_line = $out_line . "    END IF\n";
    $out_line = $out_line . "    CALL $_()\n";
    $out_line = $out_line . "    IF (pio_tf_retval_utest_ /= 0) THEN\n";
    $out_line = $out_line . "      pio_tf_nerrs_total_ = pio_tf_nerrs_total_ + 1\n";
    $out_line = $out_line . "    END IF\n";
    $out_line = $out_line . "    IF (pio_tf_world_rank_ == 0) THEN\n";
    $out_line = $out_line . "      IF (pio_tf_retval_utest_ == 0) THEN\n";
    $out_line = $out_line . "        WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF:&\n";
    $out_line = $out_line . "          Test $cur_test_case_num:\",&\n";
    $out_line = $out_line . "          \"$_\",\"-----------\", \"PASSED\"\n";
    $out_line = $out_line . "      ELSE\n";
    $out_line = $out_line . "        WRITE(*,PIO_TF_TEST_RES_FMT) \"PIO_TF:&\n";
    $out_line = $out_line . "          Test $cur_test_case_num:\",&\n";
    $out_line = $out_line . "          \"$_\",\"-----------\", \"FAILED\"\n";
    $out_line = $out_line . "      END IF\n";
    $out_line = $out_line . "    END IF\n";
    $cur_test_case_num += 1;
  }
  $out_line = $out_line . "  END SUBROUTINE PIO_TF_Test_driver_\n";
  return $out_line;
}

# Header is just a disclaimer for the generated code
sub get_header
{
  my ($template_fname) = @_;
  my ($out_line);
  $out_line = "! DON'T MODIFY THIS FILE, ALL YOUR CHANGES WILL BE LOST\n";
  $out_line = $out_line . "! This file is generated by $0\n";
  $out_line = $out_line . "! from $template_fname\n\n";

  return $out_line;
}

# The footer always contains the default test main code
# The footer can contain the default test driver code is none is specified
# - The default test driver code will contain all the auto test subs
# If a test driver code is specified the list of auto test funcs has already
# been appended the driver
sub get_footer
{
  my($ref_auto_funcs_list) = @_;
  my($out_line);
  if($template_has_test_driver == 0){
    # Add default test driver
    $out_line = &get_default_test_driver($ref_auto_funcs_list);
  }

  $out_line = $out_line . &get_default_test_main();
  return $out_line;
}

# This function processes the template file
# * Reads the template file and writes the transformed output
# * Also does some basic processing of inputs
sub process_template_file
{
  my($ifname, $ofname) = @_;
  # The var below keeps track of input file line number
  my($ifline_num);
  my($ifline);
  my($header,$footer);
  my(@auto_funcs_list);
  my($TEMPLATE_FILE);
  my($base_file_name);
  my($in_line) = "";
  my($orig_line) = "";

  $base_file_name = basename($ifname);

  # Read template file
  open(TEMPLATE_FILE, "< $ifname") or die "Cannot open input template file: $ifname";
  open(OUTPUT_FILE, "> $ofname") or die "Cannot open output file: $ofname";

  $header = get_header($ifname);
  print OUTPUT_FILE $header;

  $ifline_num = 1;
  $in_line = "";
  # Unmodified (orig) line
  $orig_line = "";
  while(<TEMPLATE_FILE>){
    my ($out_line);
  
    $orig_line = $orig_line . $_;
    chomp;
    $_ = $in_line . $_;
    if(/^s*[!#].*$/s){
      # Line starting with comment & preproc directives - print it as it is
      if($verbose){ print "Reading comment/preproc-directive: $_\n"; }
      $out_line = $_;
    }
    else{
      # Strip off comments at the end of the line
      if(/^([^!]*)!(.*)$/s){
        if($verbose){ print "Stripping comment: $2\n"; }
        $_ = $1; 
      }
      if(/^\s*PROGRAM.*/si){
        # Report error and exit if the template has a PROGRAM keyword
        print "Error: A PROGRAM keyword in statements is invalid in a template file, use PIO_TF_TEST_DRIVER instead\n";
        exit;
      }
      elsif(/^\s*$/s){
        # Print blank lines as "\n" - keep the template formatting as it is
        $out_line = "\n";
      }
      elsif(/^(.*)&$/s){
        # Fortran allows multi-line stmts with a & at the end
        # Collect all these parts of stmts into a single line before processing
        $in_line = $in_line . $1;
        $ifline_num += 1;
        next;
      }
      else{
        my ($is_transformed) = 0;
        $out_line = &transform_src(\@auto_funcs_list, $_, $base_file_name,
                        $ifline_num, \$is_transformed);
        if($is_transformed == 0){
          # If the line is not transformed, undo all changes done and
          # output the original line[s] read from the file
          chomp($orig_line);
          $out_line = $orig_line;
        }
        if($annotate_source){
          # Annotate source with template file name and line number
          # Add it after the line
          $out_line = $out_line . "   ! $base_file_name:$ifline_num";
        }
      }
    }
    print OUTPUT_FILE "$out_line\n";
    $in_line = "";
    $orig_line = "";
    $ifline_num += 1;
  }

  $footer = &get_footer(\@auto_funcs_list);
  print OUTPUT_FILE $footer;
}

sub print_usage_and_exit()
{
  print "\n\nUsage: ./pio_tf_f90gen.pl --annotate-source --out=<OUTPUT_FILE_NAME> <INPUT_TEMPLATE_FILE_NAME> \n\n";
  print "eg: ./pio_tf_f90gen.pl --annotate-source --out=pio_init_finalize.F90 pio_init_finalize.F90.in\n";
  exit;
}

# Main program

# Read input args
GetOptions(
  # Annotate generated source with template line numbers etc
  "annotate-source"		=>	\$annotate_source,
  "out=s"				      =>	\$output_fname,
  "verbose"           =>  \$verbose
);

$nargs = @ARGV;
if($nargs == 0){
  &print_usage_and_exit();
}

if($output_fname eq ""){
  &print_usage_and_exit();
}

$template_fname = shift;
if(!stat($template_fname)){
  print "Error: Cannot find input template file\n";
  &print_usage_and_exit();
}

if($verbose){ print "Reading input args complete\n" }
&process_template_file($template_fname, $output_fname);

