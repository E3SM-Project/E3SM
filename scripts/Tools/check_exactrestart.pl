#!/usr/bin/env perl
$ENV{SHELL}="/bin/csh";  # Some of the comands assume csh syntax
#### compare the two log files with check_exactrestart.pl  $file1 $file2 

if($#ARGV != 1){
   print "Usage: check_exactrestart.pl  file1 file2\n";
   print "ERROR\n";
   exit;
}
## compare the two log files
   $f1=shift(@ARGV);
   $f2=shift(@ARGV);
   if($f1 eq $f2){
      print "Error: $f1 and $f2 are the same file.\n";
      print "FAIL \n";
      exit(0);
   }
   if($f2=~"\.gz"||$f1=~"\.gz") {
      print "Error: Files still gzipped\n";
      print "ERROR\n";
      exit(0);
   }
   compare($f1,$f2);

   exit(0);

###end of the main program

# This subroutine searches the cpl.log file for the first model date record
# and returns the date and the offset from the top of the file to that point.
sub findmodeldate{
    my($file, @date) = @_;

    open(F,$file) or die "Could not open $file";
    my @firstdate;
    while(<F>){
	if(defined @date){
	    if(/model date =\s+$date[1]\s+$date[2] wall/){
		$firstdate[0] = tell F;
		last;
	    }
	}else{
	    if(/model date =\s+(\d+)\s+(\d+) wall/){
		$firstdate[0] = tell F;
		$firstdate[1] = $1;
		$firstdate[2] = $2;
		last;
	    }
	}
    }
    close(F);
    return @firstdate;
}

# Given two cpl log files this subroutine aligns file pointers on the first matching
# date record.

sub alignmodeldate{
    my ($f1, $f2) = @_;
    my @date1 = findmodeldate($f1);
    my @date2 = findmodeldate($f2);

    my @align;
    if($date1[0] < $date2[0]) {
	push(@align,findmodeldate($f1,@date2));
	$align[1]=$date2[0];
    }else{
	$align[0] = $date1[0];
	push(@align, findmodeldate($f2,@date1));
    }
    return @align;

}
    

#sub compare
#
#  This subroutine compares two cpl log files by aligning the files
# on the first common model date field and then comparing diagnostics
# from that point forward in the file.   If any diag fields are found to be different
# the next 10 diag fields of the file are printed to stdout and the test fails. 
# The test also fails if no diag fields are found after that aligned point.
#
sub compare{
  local($f1,$f2)=@_;

  my @align = alignmodeldate($f1, $f2);

  open(F1,$f1) or die "Could not open $f1";
  open(F2,$f2) or die "Could not open $f2";

  seek(F1,$align[0]-100,0);
  seek(F2,$align[1]-100,0);

  my $line1;
  my @line2 = <F2>;

  my $line2;
  my $curdate;
  my $curtime;
  my $good_cnt=0;
  my $bad_cnt=0;
  foreach $line1 (<F1>){
      chomp $line1;
      if($line1 =~ /model date =\s+(\d+)\s+(\d+) wall/){
	  $curdate = $1;
	  $curtime = $2;
      }
      next unless ($line1 =~ /comm_diag/);
      $line2 = shift @line2;
      chomp $line2;
      while ($line2 !~ /comm_diag/) {
        #print "skip line2 $line2\n";
        $line2 = shift @line2;
        chomp $line2;
      }
      if($line1 eq $line2) {
	  $good_cnt++;
      }
      else {
	  print "Difference found beginning at $curdate $curtime :\n" if($bad_cnt==0);
	  print "< $line1\n> $line2\n";
	  $bad_cnt++;
      }
       last if($bad_cnt > 10);
  }

  if($bad_cnt>0){
      print "FAIL \n";
  }
  elsif($good_cnt==0 && !($ENV{CASEBASEID} =~ /\.S\./)){
      print "ERROR: No lines compared\n";
      print "FAIL \n";
  }
  else{
      print "PASS \n";
  }


}##end of compare
