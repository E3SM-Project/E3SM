#!/usr/bin/perl 
use strict;

my $netcdf_inc = "$ENV{NETCDF}/include/netcdf.h";
my $pnetcdf_inc = "$ENV{PNETCDF}/include/pnetcdf.h";

my $functions;

my $typemap = {text=>"MPI_CHAR", short=>"MPI_SHORT", uchar=>"MPI_UNSIGNED_CHAR",
                             schar =>"MPI_CHAR", float=>"MPI_FLOAT", ushort=>"MPI_UNSIGNED_SHORT",
                             ulonglong=>"MPI_UNSIGNED_LONG_LONG",longlong=>"MPI_LONG_LONG",
                             double=>"MPI_DOUBLE",uint=>"MPI_UNSIGNED",int=>"MPI_INT",long=>"MPI_LONG"};


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

open(F,">pio_put_nc.c");
print F "#include <pio.h>
#include <pio_internal.h>\n\n";

  open(T,"pio_c_put_template.c") or die "Could not open file pio_c_put_template.c";
  my @template = <T>;
  close(T);

foreach my $func (keys %{$functions}){
  next unless($func=~/put_v/);

  my @tmptmp = @template;

  if($functions->{$func}{pnetcdf} =~ /ngattsp/){
      $functions->{$func}{netcdf} =~ s/nattsp/ngattsp/g;
  }

  if(defined ($functions->{$func}{pnetcdf}) && defined($functions->{$func}{netcdf})){
#      && defined($functions->{$func}{pio})){
      my $bfunc = "b".$func;
      my $pnetfunc;
      if(defined $functions->{$bfunc}{pnetcdf}){
	  $pnetfunc = $functions->{$bfunc}{pnetcdf};
      }else{
	  $pnetfunc = $functions->{$func}{pnetcdf};
	  undef $bfunc;
      }

      my $buf=0;
      foreach my $line (@tmptmp){
	  if($line =~ /msg = 0/){
	      my $msg = "PIO_MSG_".$func;
	      $msg =~ tr/a-z/A-Z/;
	      $line =~ s/msg = 0/msg = $msg/;
#	      print "  $msg,\n";
	  }
	  if($line =~/function/){
	      if($line =~ s/PIO_function/PIOc_$func/){
		  $line =~ s/\(\)/ $pnetfunc / ; 
#	      $line =~ s/int ncid/file_desc_t *file/;
		  $line =~ s/MPI_Offset/PIO_Offset/g;
		  $line =~ s/\;//;

# PIO handles requests internally, remove it from the argument list

		  $line =~ s/, int \*request//;

		  
	      }else{
		  my $args;
		  if(defined $bfunc){
		      if($line =~ s/ncmpi_function/ncmpi_$bfunc/){
			  $args = $pnetfunc ;
			  $args =~ s/request/&request/;
		      }
		  }else{
		      if($line =~ s/ncmpi_function/ncmpi_$func/){
			  $args = $pnetfunc ;
		      }
		  }
		  if($line =~ s/nc_function/nc_$func/){
		      $args = $functions->{$func}{netcdf}  ; 
		      if($pnetfunc =~ /void \*buf/  ){
			  $args =~ s/\*op/ buf/g;
			  $args =~ s/\*ip/ buf/g;
		      }
		      if($pnetfunc =~ /\*ip/  ){
			  $args =~ s/\*op/ ip/g;
		      }
		  }
		  $args =~ s/int ncid/file->fh/;
#      $args =~ s/int varid/vdesc->varid/;
		  $args =~ s/const //g;
		  $args =~ s/unsigned //g;
		  $args =~ s/signed //g;
		  $args =~ s/char //g;
		  $args =~ s/nc_type //g;
		  $args =~ s/MPI_Offset //g;
		  $args =~ s/MPI_Datatype //g;
		  $args =~ s/int //g;
		  $args =~ s/short //g;
		  $args =~ s/long //g;
		  $args =~ s/float //g;
		  $args =~ s/double //g;
		  $args =~ s/void//g;
		  $args =~ s/MPI_Info //g;
		  $args =~ s/,\s+\*/, /g;
		  $args =~ s/\[\]//g;
		  
		  $args =~ s/size_t \*indexp/\(size_t \*\) index/g;
		  $args =~ s/size_t \*startp/\(size_t \*\) start/g;
		  $args =~ s/size_t \*countp/\(size_t \*\) count/g;
		  $args =~ s/ptrdiff_t \*stridep/\(ptrdiff_t \*\) stride/g;
		  $args =~ s/ptrdiff_t \*\s*imapp/\(ptrdiff_t \*\) imap/g;


#		  $args =~ s/ptrdiff_t \*/\(ptrdiff_t \*\)/g;
#		  if($args =~ s/size_t \*/\(size_t \*\)/g){
#		  }else{
#		      $args =~ s/size_t /\(size_t\)/g;
#		  }
		  
		  $line =~ s/\(\)/$args/;
	      }
	  }
#	  if($line =~ /chkerr =/ && $func =~ /def_var/){
#	      print F "  file->varlist[*varidp].record=-1;\n";
#	      print F "  file->varlist[*varidp].buffer=NULL;\n";
#	  }
	  print F $line;	  

      }
    

      print F "\n";
  }
}
close(F);


open(F,">pio_get_nc.c");
print F "#include <pio.h>
#include <pio_internal.h>\n\n";
  open(T,"pio_c_get_template.c") or die "Could not open file pio_c_get_template.c";
  my @template = <T>;
  close(T);

foreach my $func (keys %{$functions}){
  next unless($func=~/get_var/);

  my $bufcount;
  my $buftype;
  if($func =~ /_([^_]+)\s*$/ || $func =~ /_([^_]+)_all/){
      my $nctype = $1;
#      print "nctype $nctype\n";
      $buftype = $typemap->{$nctype};
  }
  
  if($func =~ /var1/){
      $bufcount = 1;
  }
#  print "$func $buftype\n";

  my @tmptmp = @template;

  if($functions->{$func}{pnetcdf} =~ /ngattsp/){
      $functions->{$func}{netcdf} =~ s/nattsp/ngattsp/g;
  }

  if(defined ($functions->{$func}{pnetcdf}) && defined($functions->{$func}{netcdf})){
#      && defined($functions->{$func}{pio})){
      my $buf=0;
      foreach my $line (@tmptmp){
	  if($line =~ /PIO_NOERR/){
	      if( $functions->{$func}{pnetcdf} =~ /buf/){
		  print F "  ibufcnt = bufcount;\n";
		  print F "  ibuftype = buftype;\n";
	      }else{
		  print F "  ibuftype = $buftype;\n";
		  if(defined($bufcount)){
		      print F "  ibufcnt = $bufcount;\n";
		  }elsif($func =~ /vara/){
		      print F "  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);\n";
		      print F "  ibufcnt = 1;\n";
		      print F "  for(int i=0;i<ndims;i++){\n";
		      print F "    ibufcnt *= count[i];\n";
		      print F "  }\n";
		  }elsif($func =~ /vars/ or $func =~ /varm/){
		      print F "  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);\n";
		      print F "  ibufcnt = 1;\n";
		      print F "  for(int i=0;i<ndims;i++){\n";
		      print F "    ibufcnt *= count[i]/stride[i];\n";
		      print F "  }\n";
		  }elsif($func =~ /var_?/){
		      print F "  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);\n";
		      print F "  int dimid[ndims];\n";
		      print F "  PIO_Offset dimsize;\n";
		      print F "  ibufcnt = 1;\n";
		      print F "  PIOc_inq_vardimid(file->fh, varid, dimid);\n";
		      print F "  for(int i=0;i<ndims;i++){\n";
		      print F "    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);\n";
		      print F "    ibufcnt *= dimsize;\n";
		      print F "  }\n";
		  }
	      }


	  }
	  if($line =~ /msg = 0/){
	      my $msg = "PIO_MSG_".$func;
	      $msg =~ tr/a-z/A-Z/;
	      $line =~ s/msg = 0/msg = $msg/;
#	      print "  $msg,\n";
	  }
	  if($line =~/function/){
	      if($line =~ s/PIO_function/PIOc_$func/){
		  $line =~ s/\(\)/ $functions->{$func}{pnetcdf} / ; 
#	      $line =~ s/int ncid/file_desc_t *file/;
		  $line =~ s/MPI_Offset/PIO_Offset/g;
		  $line =~ s/\;//;
		  $line =~ s/\*ip/\*buf/;
	      }else{
		  my $args;
		  my $tall = $func."_all";
		  if(defined($functions->{$tall}{pnetcdf})){
		      if($line =~ s/ncmpi_function/ncmpi_$tall/){
			  $args = $functions->{$tall}{pnetcdf};
		      }
		  }else{
		      if($line =~ s/ncmpi_function/ncmpi_$func/){
			  $args = $functions->{$func}{pnetcdf} ;
		      }
		  }
		  if($line =~ s/nc_function/nc_$func/){
		      $args = $functions->{$func}{netcdf}  ; 
		  }
		  $args =~ s/int ncid/file->fh/;

		  if($functions->{$func}{pnetcdf} =~ /void \*buf/  ){
		      $args =~ s/\*op/ buf/g;
		      $args =~ s/\*ip/ buf/g;
		  }
		  if($functions->{$func}{pnetcdf} =~ /\*ip/  ){
		      $args =~ s/\*ip/ buf/g;
		  }
		  if($functions->{$func}{pnetcdf} =~ /\*op/  ){
		      $args =~ s/\*op/ buf/g;
		  }
#      $args =~ s/int varid/vdesc->varid/;
		  $args =~ s/const //g;
		  $args =~ s/unsigned //g;
		  $args =~ s/signed //g;
		  $args =~ s/char //g;
		  $args =~ s/nc_type //g;
		  $args =~ s/MPI_Offset //g;
		  $args =~ s/MPI_Datatype //g;
		  $args =~ s/int //g;
		  $args =~ s/short //g;
		  $args =~ s/long //g;
		  $args =~ s/float //g;
		  $args =~ s/double //g;
		  $args =~ s/void//g;
		  $args =~ s/MPI_Info //g;
		  $args =~ s/,\s+\*/, /g;
		  $args =~ s/\[\]//g;
		  
		  $args =~ s/size_t \*indexp/\(size_t \*\) index/g;
		  $args =~ s/size_t \*startp/\(size_t \*\) start/g;
		  $args =~ s/size_t \*countp/\(size_t \*\) count/g;
		  $args =~ s/ptrdiff_t \*stridep/\(ptrdiff_t \*\) stride/g;
		  $args =~ s/ptrdiff_t \*\s*imapp/\(ptrdiff_t \*\) imap/g;


#		  $args =~ s/ptrdiff_t \*/\(ptrdiff_t \*\)/g;
#		  if($args =~ s/size_t \*/\(size_t \*\)/g){
#		  }else{
#		      $args =~ s/size_t /\(size_t\)/g;
#		  }
		  
		  $line =~ s/\(\)/$args/;
	      }
	  }
#	  if($line =~ /chkerr =/ && $func =~ /def_var/){
#	      print F "  file->varlist[*varidp].record=-1;\n";
#	      print F "  file->varlist[*varidp].buffer=NULL;\n";
#	  }
	  print F $line;


      }
    
#    print "$func  $functions->{$func}{pnetcdf} $functions->{$func}{netcdf}\n";
      print F "\n";
  }
}
close(F);
