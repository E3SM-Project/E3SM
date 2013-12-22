#!/usr/bin/perl 
use strict;

my $netcdf_inc = "$ENV{NETCDF}/include/netcdf.h";
my $pnetcdf_inc = "$ENV{PNETCDF}/include/pnetcdf.h";

my $functions;

open(PF,$pnetcdf_inc) or die "Could not open $pnetcdf_inc";
my @file = <PF>;
my $func="";
foreach my $line (@file){
  chomp($line);
  next if($line =~ /^\s*\/\*/);
  next if($line =~ /^\s*[#\*]/);
  $line =~ s/\s+/ /g; 
  if($line =~ / ncmpi_(.*)(\(.*)$/) {
    $func = $1;

    $functions->{$func}{pnetcdf} = $2;
  }elsif($func ne ""){
    $functions->{$func}{pnetcdf} .= $line;
  }
  $func="" if($line =~ /\;\s*$/);
}
close(PF);

open(NF,$netcdf_inc) or die "Could not open $netcdf_inc";
my @file = <NF>;

my $func="";
foreach my $line (@file){
  chomp($line);
  next if($line =~ /^\s*\/\*/);
  next if($line =~ /^\s*[#\*]/); 
  $line =~ s/\s+/ /g;
  if($line =~ /^\s*nc_([^_].*)(\(.*)$/) {
    $func = $1;
    $functions->{$func}{netcdf} = $2;
  }elsif($func ne ""){
    $functions->{$func}{netcdf} .= $line;
  }
  $func="" if($line =~ /\;\s*$/);
}
close(NF);

my $func="";
open(PIO,"nf_mod.F90") or die "Could not open file nf_mod.F90";
my @nfmod = <PIO>;
close(PIO);
foreach my $line (@nfmod){
    if ($line =~ /\s*integer\s+function\s+pio_([^\(]+)\(/i){
	$func = $1;
	$func =~ tr/A-Z/a-z/;
	$functions->{$func}{pio} = $line;
    }
    elsif ($line =~ /\s*integer\s+function\s+([^\(]+)\(/i){
	$func = $1;
	$func =~ tr/A-Z/a-z/;
	$functions->{$func}{pio} = $line;
    }

}



open(F,">pio_nc.c");


print "$functions->{inq}{pnetcdf}\n";
print "$functions->{inq}{netcdf}\n";
#print "$functions->{inquire}{pio}\n";
  open(T,"pio_c_template.c") or die "Could not open file pio_c_template.c";
  my @template = <T>;
  close(T);

foreach my $func (keys %{$functions}){
  next if($func=~/get_v/ or $func=~/put_v/);

  my @tmptmp = @template;

  if(defined ($functions->{$func}{pnetcdf}) && defined($functions->{$func}{netcdf})){
#      && defined($functions->{$func}{pio})){

    foreach my $line (@tmptmp){
      if($line =~/function/){
      if($line =~ s/PIO_function/PIO_$func/){
#	$line =~ s/\(\)/ $functions->{$func}{pio} / ; 
      }
      if($line =~ s/ncmpi_function/ncmpi_$func/){
	$line =~ s/\(\)/ $functions->{$func}{pnetcdf} / ;
      }
      if($line =~ s/nc_function/nc_$func/){
	$line =~ s/\(\)/ $functions->{$func}{netcdf} / ; 
      }
      $line =~ s/int ncid/file->fh/;
#      $line =~ s/int varid/vdesc->varid/;
      $line =~ s/const char \*//g;
      $line =~ s/char \*//g;
      $line =~ s/nc_type \*//g;
      $line =~ s/nc_type //g;
      $line =~ s/MPI_Offset \*//g;
      $line =~ s/const int \*//g;
      $line =~ s/int \*//g;
      $line =~ s/ int //g;
      $line =~ s/size_t \*//g;
    }
 
     print F $line;
  }

#    print "$func  $functions->{$func}{pnetcdf} $functions->{$func}{netcdf}\n";
  }
}
close(F);
