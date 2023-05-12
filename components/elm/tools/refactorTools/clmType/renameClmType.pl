#!/usr/bin/env perl

# 
# script that will flatten pft type
#

use strict;
use Cwd;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);
use English;
use IO::File;

my %opts = ( file => undef);
GetOptions( "file=s" => \$opts{'file'}) or die "need to specify file as input";
my $file = $opts{'file'};

my  $fh = new IO::File;
$fh->open("<$file") or die "** can't open file: $file\n";

my $fhout = new IO::File;
$fhout->open(">$file.temp1") or die "** can't open file: $file.temp1\n";

my %subst = ();
while (my $line = <$fh>) {
	
    if ($line =~ /p%active/) {
	$line =~ s/p%active/pft%active/g;
    }
    if ($line =~ /pactive.*[^:]/) {
	$line =~ s/ active/pft%active/;
    }
    if ($line =~ /c%active/) {
	$line =~ s/c%active/col%active/g;
    }
    if ($line =~ /cactive/) {
	$line =~ s/\sactive/col%active/g;
    }

    # --------------------------------------------------------
    if ($line =~ /\bp\b *=\>(.+)/) {
	$line =~ s/$1/ pft/;
    }
    if ($line =~ /\bc\b *=\>(.+)/) {
	$line =~ s/$1/ col/;
    }
    if ($line =~ /\bl\b *=\>(.+)/) {
	$line =~ s/$1/ lun/;
    }
    if ($line =~ /\bg\b *=\>(.+)/) {
	$line =~ s/$1/ grc/;
    }

    # --------------------------------------------------------
    # substitute pptr => clm3%g%l%c%p with pptr => pft
    # --------------------------------------------------------
    if ($line =~ / pptr *=\>(.+)/) {
	$line =~ s/$1/ pft/;
    }
    if ($line =~ / cptr *=\>(.+)/) {
	$line =~ s/$1/ col/;
    }
    if ($line =~ / lptr *=\>(.+)/) {
	$line =~ s/$1/ lun/;
    }
    if ($line =~ / gptr *=\>(.+)/) {
	$line =~ s/$1/ grc/;
    }

    # --------------------------------------------------------
    # special clm3% substitution before general clm3% subst
    # --------------------------------------------------------
    $line =~ s/clm3%g%l%itype/lun%itype/g;
    $line =~ s/clm3%g%l%c%itype/col%itype/g;

    # -----------------------------------------
    # substitute clm3%g%l%c%p%TYPE% to TYPE%
    # -----------------------------------------
    if ($line =~ /clm3%g%l%c%p%.+%/) {
	$line =~ s/clm3%g%l%c%p%//g;
    }	
    if ($line =~ /clm3%g%l%c%.+%/) {
	$line =~ s/clm3%g%l%c%//g;
    }	
    if ($line =~ /clm3%g%l%.+%/) {
	$line =~ s/clm3%g%l%//g;
    }	
    if ($line =~ /clm3%g%.+%/) {
	$line =~ s/clm3%g%//g;
    }

    # -----------------------------------------
    # substitute clm3%g%l%c%p%TYPE to pft%TYPE
    # -----------------------------------------
    if ($line =~ /clm3%g%l%c%p%.+[^%]/) {
	$line =~ s/clm3%g%l%c%p%/ pft%/g;
    }	
    if ($line =~ /clm3%g%l%c%.+[^%]/) {
	$line =~ s/clm3%g%l%c%/ col%/g;
    }	
    if ($line =~ /clm3%g%l%.+[^%]/) {
	$line =~ s/clm3%g%l%/ lun%/g;
    }	
    if ($line =~ /clm3%g%.+[^%]/) {
	$line =~ s/clm3%g%/ grc%/g;
    }

    # -----------------------------------------
    if ($line =~ /\b(g%).+%/) {
	$line =~ s/$1//g;
    }
    if ($line =~ /\b(l%).+%/) {
	$line =~ s/$1//g;
    }
    if ($line =~ /\b(c%).+%/) {
	$line =~ s/$1//g;
    }
    if ($line =~ /\b(p%).+%/) {
	$line =~ s/$1//g;
    }

    # -----------------------------------------
    $line =~ s/ p%(.+)/pft%$1/g;
    $line =~ s/ c%(.+)/col%$1/g;
    $line =~ s/ l%(.+)/lun%$1/g;
    $line =~ s/ g%(.+)/grc%$1/g;

    # --------------------------------------------------------
    # substitution for (p% => (pft%
    # --------------------------------------------------------
    $line =~ s/([^a-zA-Z0-9_])p%/$1pft%/; 
    $line =~ s/([^a-zA-Z0-9_])c%/$1col%/; 
    $line =~ s/([^a-zA-Z0-9_])l%/$1lun%/; 
    $line =~ s/([^a-zA-Z0-9_])g%/$1grc%/; 

    # --------------------------------------------------------
    # another special substitution
    # --------------------------------------------------------
    if ($line =~ /call *CNSet/){
	$line =~ s/pft%/ /;
	$line =~ s/col%/ /;
	$line =~ s/lun%/ /;
	$line =~ s/grc%/ /;
    }	

   # -----------------------------------------
    # substitute ccf%pcf_a to pcf_a example
    # -----------------------------------------
    $line =~ s/cc13f%pcf_a/pc13f_a/g;
    $line =~ s/cc14f%pcf_a/pc14f_a/g;

    $line =~ s/..f%p(.)f_a/p$1f_a/g;
    $line =~ s/..s%p(.)s_a/p$1s_a/g;

    # ----------------------------------------
    $line =~ s/pptr%(p[a-zA-Z0-9_]+%)/$1/g;
    $line =~ s/cptr%(c[a-zA-Z0-9_]+%)/$1/g;
    $line =~ s/lptr%(l[a-zA-Z0-9_]+%)/$1/g;
    $line =~ s/gptr%(g[a-zA-Z0-9_]+%)/$1/g;

    print $fhout "$line"; 
}

$fh->close();
$fhout->close();

my $sysmod = "/bin/mv $file.temp1 $file"; 
system($sysmod) == 0 or die "ERROR: $sysmod failed: $?\n";

# =========================================================================

my  $fh = new IO::File;
$fh->open("<$file") or die "** can't open file: $file\n";

my $fhout = new IO::File;
$fhout->open(">$file.temp1") or die "** can't open file: $file.temp1\n";

my $lmax1 = 0;
my $lmax2 = 0;
my $lmax3 = 3;
while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /(.*)(=\>)([^,]*)([, ] [\& ].*\!.+)/ ) {
	my $var1 = $1;
	my $var2 = $3;
	my $var3 = $4;
	$var1 =~ s/\s+$//g;
	$var2 =~ s/\s+$//g;
	$var3 =~ s/\s+$//g;
	$var1 =~ s/^\s+//; #remove leading spaces
	$var2 =~ s/^\s+//; #remove leading spaces
	$var3 =~ s/^\s+//; #remove leading spaces
	my $lvar1 = length($var1);
	my $lvar2 = length($var2);
	my $lvar3 = length($var3);
	$lmax1 = ($lvar1 > $lmax1) ? $lvar1 : $lmax1;
	$lmax2 = ($lvar2 > $lmax2) ? $lvar2 : $lmax2;
	$lmax3 = ($lvar3 > $lmax3) ? $lvar3 : $lmax3;
    }
}
my $max = $lmax1 + $lmax2 + $lmax3; 
print " lmax1,lmax2,lmax3 are $lmax1,$lmax2,$lmax3 \n";    

if ($lmax1 > 0 && $lmax2 > 0) {
    seek $fh, 0, 0;
    while (my $line = <$fh>) {
	chomp $line;
	if ($line =~ /(.*)(=\>)([^,]*)([, ] \& \!.+)/) {
	    my $var1 = $1;
	    my $arrow = $2;
	    my $var2 = $3;
	    my $var3 = $4;
	    $var1 =~ s/\s+$//g;
	    $var2 =~ s/\s+$//g;
	    $var3 =~ s/\s+$//g;
	    $var1 =~ s/^\s+//; #remove leading spaces
	    $var2 =~ s/^\s+//; #remove leading spaces
	    $var3 =~ s/^\s+//; #remove leading spaces
	    my $format = "   %-${lmax1}s %-3s %-${lmax2}s %-${lmax3}s\n";
	    printf $fhout ($format,$var1,$arrow,$var2,$var3); 
	} else {
	    print $fhout "$line\n";
	}
    }
    my $sysmod = "/bin/mv $file.temp1 $file"; 
    system($sysmod) == 0 or die "ERROR: $sysmod failed: $?\n";
}

$fh->close();
$fhout->close();
