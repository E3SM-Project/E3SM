#!/usr/bin/env perl

use strict;
use Cwd;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);
use English;
use IO::File;
use Fcntl;

my %opts = ( file => undef);
GetOptions( "file=s" => \$opts{'file'}) or die "need to specify file as input";
my $file = $opts{'file'};

#----------------------------------------------------------------------
# Add input/inout/output comments to pointer declarations
#----------------------------------------------------------------------

my $fh = new IO::File;
$fh->open("<$file") or die "** can't open file: $file\n";

my $fhout = new IO::File;
$fhout->open(">$file.temp1") or die "** can't open file: $file.temp1\n";

my $select_case = 0;
my $type;
while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /.*pointers to .*implicit in /) {
	$type = 'Input: '; 
    }
    if ($line =~ /.*pointers to .*implicit in\/out /) {
	$type = 'InOut: '; 
    }
    if ($line =~ /.*pointers to .*implicit out /) {
	$type = 'Output:'; 
    }

    if ($line =~ /(.+)(.*,.*pointer.*::.*)(\(:.*\))(.*\!)(.*)/) {
	my $irl = $1;
	$irl = clean($irl);
	print $fhout "$1 $2 $3 $4 $type [$irl $3] $5 \n";
    } elsif ($line =~ /(.+)(.*,.*pointer.*::.*)(\(:.*\))( *$)/) {
	my $irl = $1;
	$irl = clean($irl);
	print $fhout "$1 $2 $3 \! $type [$irl $3] $4\n";
    } elsif ($line =~ /(A|a)ssign local/) {
	# do nothing 
    } else {
	print $fhout "$line\n";
    }
}
$fhout->close();
$fh->close();

#----------------------------------------------------------------------
# Add input/inout/output comments to current pointer associations (=>)
# Output is $file.temp2
#----------------------------------------------------------------------

my  $fh = new IO::File;
$fh->open("<$file.temp1") or die "** can't open file: $file.temp1\n";

my $fhout = new IO::File;
$fhout->open(">$file.temp2") or die "** can't open file: $file.temp2\n";

my $tell = 0;
my $tellm1 = 0;
my $end_file = 0;
my $is_module = 0;
my $module_name; 
while(1) {
    seek $fh, $tell, 0;
    my %subst = ();
    while (my $line = <$fh>) {
	if ($line =~ /(^module)(.+)/) {
	    $is_module = 1;
	    $module_name = $2;
	    print "module name is $module_name\n";
	}
	if ($line =~ /(.+pointer.*::)(.*)(\!.*\n)/) {
	    my $var = $2; 
	    my $comment = $3; 
	    $comment = clean($comment);
	    $var =~ s/\(.*\)//g; 
	    $var = clean($var);
	    $subst{$var} = $comment;
	}
	if ($line =~ /end .*subroutine/) {
	    $tell = tell($fh);
	    last;
	}
    }
    
    seek $fh, $tellm1, 0;
    my @vars = keys %subst;
    my $select_case = 0;
    while (my $line = <$fh>) {
	if ($line =~ /(.+pointer.*::)(.*)(\!.*\n)/) {
	    my $newline = backwards_compatibility ($line, $module_name); 
	    if ($newline) {
		print $fhout $newline;
	    } else {
		print $fhout "$line";
	    }
	} elsif ($line =~ /select case/) {
	    $select_case = 1;
	    print $fhout "$line";
	} elsif ($line =~ /end select/) {
	    $select_case = 0;
	    print $fhout "$line";

	} elsif ($line =~ /(.+)(=\>)(.+%.+)/) {
	    if (!$select_case) {
		my $varmatch = $1;
		$varmatch = clean($varmatch);
		foreach my $var (sort @vars) {
		    $var = clean($var);
		    if ($varmatch eq $var) {
			my $add = $subst{$var};
			$line = $varmatch . " => $3 , & $add";
			$line =~ s/\n//g;
			$line =~ /(.*)(=\>)([^,]*)(\,.+)/;
			printf $fhout ("   %-35s %-3s %-45s %-80s\n",$1,$2,$3,$4); 
		    }
		}
	    } else {
		print $fhout "$line";
	    }		    
	} else {
	    print $fhout "$line";
	}
	if ($line =~ /end .*subroutine/) {
	    if ($is_module) {
		$tellm1 = tell($fh);
	    } else {
		$end_file = 1;
	    }		
	    last;
	}
	if ($line =~ /end .*module/) {
	    $end_file = 1;
	    last;
	}
    }
    if ($end_file) {
	last;
    }
}

$fh->close();
$fhout->close();

#-------------------------------------------------------------------
# Read in file - fill in $line[] array - but skip lines that are
# pointer declarations - after "pointers to implicit" and before
# "LOCAL VARIABLE"
# Add associate(& and ) statements around subroutine declarations
# Output file is $file.temp3
#-------------------------------------------------------------------

my  $fh = new IO::File;
$fh->open("<$file.temp2") or die "** can't open file: $file.temp2\n";

my $fhout = new IO::File;
$fhout->open(">$file.temp3") or die "** can't open file: $file.temp3\n";

my @lines;
my $n = 0;
my $old_pointer_block = 0;
seek $fh, 0, 0;
while (my $line = <$fh>) {
    if ($line =~ /.*pointers to.+implicit in /) {
	$old_pointer_block = 1;
    } 
    if ($line =~ /.*Assign.*pointer.*/) {
	# do nothing
    } elsif ($line =~ /.*local.*pointers.*to.*/) {
	# do nothing
    } elsif ($line =~ /OTHER LOCAL/) {
	# do nothing
    } elsif ($line =~ /^#endif/) {
	if ($lines[$n-1] =~ /(.+)(\=\>)(.+%.+)/) {
	    # do nothing
	} else {
	    $lines[$n] = $line;
	    $n++;
	}
    } elsif ($old_pointer_block) {
	# Do not add lines in pointer declaration to $line[] array
	if ($line !~ /.*pointer.*::.*/) {
	    $lines[$n] = $line;
	    $n++;
	}
    } else {
	$lines[$n] = $line;
	$n++;
    }
    if ($line =~ /.*LOCAL +VARIABLE.*/) {
	$old_pointer_block = 0;
    } 
}
my $nsize = scalar(@lines);

my $associate = 0;
$n = 0;
while ($n <= $nsize) {
    my $line = $lines[$n];
    if ($line =~ /(.+)(\=\>)(.+%.+)/) {
	if ($lines[$n-1] !~ /.+\=\>.+/) {
	    $line = "   associate(& \n$line";
	    $associate = 1;
	}
	if (($lines[$n+1] !~ /.+\=\>.+/) &&
	    ($lines[$n+2] !~ /.+\=\>.+/) &&
	    ($lines[$n+3] !~ /.+\=\>.+/) &&
	    ($lines[$n+4] !~ /.+\=\>.+/)) {
	    $line =~  s/(.+), & \!/$1  & \!/g;
	    $line = $line . "   )\n"; 
	}  
    } elsif ($line =~ /end *subroutine/ && ($associate)) {
	$line = "    end associate \n $line";
	$associate = 0;
    }
    $lines[$n] = $line; 
    $n = $n + 1; 
}

my $nsize = scalar(@lines);
$n = 0;
while ($n <= $nsize) {
    print $fhout $lines[$n];
    $n++;
}
$fhout->close();

#----------------------------------------------------------------------
# Remove non associate statements from associate block declarations
#----------------------------------------------------------------------

my  $fh = new IO::File;
$fh->open("<$file.temp3") or die "** can't open file: $file.temp3\n";

my $fhout = new IO::File;
$fhout->open(">$file.temp4") or die "** can't open file: $file.temp4\n";

my $associate = 0;
while (my $line = <$fh>) {
    if ($line =~ /.+associate\(&/) {
	if (!$associate) {
	    print $fhout $line;
	}
	$associate = 1;
    }
    if ($associate) {
	if ($line =~ /.+=>.+/) {
	    print $fhout $line;
	}
	if ($line =~ /^  +\)$/) {
	    print $fhout $line;
	    $associate = 0;
	}
    } else {
	print $fhout $line;
    }
}
$fhout->close();

#----------------------------------------------------------------------

my $sysmod = "/bin/rm $file.temp1"; 
system($sysmod) == 0 or die "ERROR: $sysmod failed: $?\n";
my $sysmod = "/bin/rm $file.temp2"; 
system($sysmod) == 0 or die "ERROR: $sysmod failed: $?\n";
my $sysmod = "/bin/rm $file.temp3"; 
system($sysmod) == 0 or die "ERROR: $sysmod failed: $?\n";
my $sysmod = "/bin/mv $file.temp4 $file"; 
system($sysmod) == 0 or die "ERROR: $sysmod failed: $?\n";

#----------------------------------------------------------------------

sub clean
{
    my ($name) = @_;
    $name =~ s/^\s+//; # strip any leading whitespace 
    $name =~ s/\s+$//; # strip any trailing whitespace
    return ($name);
}

sub backwards_compatibility 
{
    my ($line, $module_name) = @_; 

    my @matches;
    $module_name = clean($module_name);

    if ($module_name eq "ch4Mod") { 
	my @temp= qw(grnd_ch4_cond  grnd_ch4_cond_col ch4_prod_depth o2_decomp_depth  co2_decomp_depth  conc_o2 conc_o2
                     conc_ch4 ch4_oxid_depth o2_oxid_depth co2_oxid_depth o2_decomp_depth ch4_aere_depth ch4_tran_depth
                     o2_aere_depth co2_aere_depth conc_ch4 conc_o2 ch4_oxid_depth ch4_prod_depth ch4_aere_depth ch4_oxid_depth
                     ch4_ebul_depth ch4_ebul_total conc_ch4 ch4_prod_depth o2_oxid_depth conc_ch4 conc_o2 
                     ch4_oxid_depth ch4_aere_depth ch4_surf_aere ch4_surf_ebul ch4_ebul_depth ch4_ebul_total ch4_surf_diff
                     o2_decomp_depth o2_aere_depth co2_decomp_depth o2stress ch4stress);
	push(@matches, @temp);
    }

    if ($module_name eq 'CanopyFluxesMod') {
	my @temp = qw(lai_z par_z vcmaxcint psn_z lmr_z rs_z ci_z psn psn_wc psn_wj psn_wp lmr rs alphapsn alphapsn par_z);
	push(@matches, @temp);
    }


    foreach my $match (@matches) {
	$match = clean($match); 
	if ($line =~ /(.+)($match)(.+)(!.*$)/) {
	    return "$1$2$3 ! needed for backwards compatiblity\n";
	}
	if ($line =~ /(.+)($match)/) {
	    return "$1$2 ! needed for backwards compatiblity\n";
	}
    }

    #---------------------------------------------------------

    @matches = ();
    if ($module_name eq 'MEGANFactorsMod') {
	@matches = qw(eff comp_factors_table hash_table_indices);
    }
    if ($module_name eq 'CNAllocationMod') {
	@matches = qw(arepr aroot col_plant_ndemand residual_plant_ndemand);
    }
    if ($module_name eq 'CNPhenologyMod') {
	 @matches = qw(inhemi);
    }
    if ($module_name eq 'CNDecompCascadeMod') {
	@matches = qw(fr);
    }
    if ($module_name eq 'CNRestMod') {
	@matches = qw(iptemp ptr1d data_r1);
    }
    if ($module_name eq 'CNDVMod') {
	@matches = qw(rbuf2dg);
    }
    if ($module_name eq 'CNDecompCascadeMod_CENTURY') {
	@matches = qw(fr);
    }
    if ($module_name eq 'CNSummaryMod') {
	@matches = qw(gpp col_gppar col_ar rr col_rr npp col_npp vegfire col_vegfire
                      wood_harvestc col_wood_harvestc totvegc col_totvegc totpftc col_totpftc
                      pft_fire_closs col_pft_fire_closs litfall col_litfall 
                      hrv_xsmrpool_to_atm col_hrv_xsmrpool_to_atm);
    }
    if ($module_name eq 'CNSoilLittVertTranspMod') {
	@matches = qw(spinup_factor is_cwd altmax altmax_lastyear 
                      som_adv_coef som_diffus_coef conc_ptr source 
                      trcr_tendency_ptr);
    }
    if ($module_name eq 'clm_varcon') {
	@matches = qw(zlak dzlak zsoi dzsoi zisoi dzsoi_decomp nlvic)
    }
    foreach my $match (@matches) {
	$match = clean($match); 
	if ($line =~ /(.+)(pointer)(.*)($match)/) {
	    return $line;
	}
    }
}
