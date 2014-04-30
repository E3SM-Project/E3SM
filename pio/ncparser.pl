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

open(F,">pio_nc.c");
print F 
"/**
* \@file   pio_nc.c
* \@Author Jim Edwards (jedwards\@ucar.edu)
* \@date     Feburary 2014 
* \@brief    PIO interfaces to netcdf support functions
*
* This file provides an interface to the netcdf support functions.  It calls the underlying netcdf or pnetcdf or netcdf4 
* functions from the appropriate subset of mpi tasks.  
*/
#include <pio.h>
#include <pio_internal.h>\n\n";

  open(T,"pio_c_template.c") or die "Could not open file pio_c_template.c";
  my @template = <T>;
  close(T);

my $typemap = {text=>"MPI_CHAR", short=>"MPI_SHORT", uchar=>"MPI_UNSIGNED_CHAR",
                             schar =>"MPI_CHAR", float=>"MPI_FLOAT", ushort=>"MPI_UNSIGNED_SHORT",
                             ulonglong=>"MPI_UNSIGNED_LONG_LONG",longlong=>"MPI_LONG_LONG",ubyte=>"MPI_BYTE",
                             double=>"MPI_DOUBLE",uint=>"MPI_UNSIGNED",int=>"MPI_INT",long=>"MPI_LONG"};



foreach my $func (keys %{$functions}){
  next if($func=~/get_v/ or $func=~/put_v/);
  next if($func=~/strerror/);
  next if($func=~/delete/);
  next if($func=~/open/);
  next if($func=~/close/);
  next if($func=~/create/);
  next if($func=~/default_format/);
  next if($func=~/copy_att/);
  next if($func=~/abort/);
  next if($func=~/string/);
  next if($functions->{$func}{pnetcdf} =~ /\(void\)/);

  my @tmptmp = @template;

  if($functions->{$func}{pnetcdf} =~ /ngattsp/){
      $functions->{$func}{netcdf} =~ s/nattsp/ngattsp/g;
  }

  if(defined ($functions->{$func}{pnetcdf}) && defined($functions->{$func}{netcdf})){
#      && defined($functions->{$func}{pio})){

      foreach my $line (@tmptmp){
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

		  print F 
"/** 
* \@name    PIOc_$func
*/\n";

	      }else{
		  my $args;
		  if($line =~ s/ncmpi_function/ncmpi_$func/){
		      $args = $functions->{$func}{pnetcdf} ;
		  }
		  if($line =~ s/nc_function/nc_$func/){
		      $args = $functions->{$func}{netcdf}  ; 
		  }
		  $args =~ s/int ncid/file->fh/;
#      $args =~ s/int varid/vdesc->varid/;
		  $args =~ s/const //g;
		  $args =~ s/unsigned //g;
		  $args =~ s/signed //g;
		  $args =~ s/char //g;
		  $args =~ s/nc_type //g;
		  $args =~ s/MPI_Offset //g;
		  $args =~ s/int //g;
		  $args =~ s/short //g;
		  $args =~ s/long //g;
		  $args =~ s/float //g;
		  $args =~ s/double //g;
		  $args =~ s/void//g;
		  $args =~ s/MPI_Info //g;
		  $args =~ s/,\s+\*/, /g;
		  if($args =~ s/size_t \*/\(size_t \*\)/g){
		  }else{
		      $args =~ s/size_t /\(size_t\)/g;
		  }
		  
		  $line =~ s/\(\)/$args/;
	      }
	  }
#	  if($line =~ /chkerr =/ && $func =~ /def_var/){
#	      print F "  file->varlist[*varidp].record=-1;\n";
#	      print F "  file->varlist[*varidp].buffer=NULL;\n";
#	  }
	  print F $line;
	  if($func =~ /sync/ && $line =~   /case PIO_IOTYPE_PNETCDF/){
	      printf F "      flush_output_buffer(file);\n";
	  }
	  if($line =~ /msg =/){
	      if($func =~ /inq_varndims/){
		  print F "  if(file->varlist[varid].ndims > 0){\n";
		  print F "    (*ndimsp) = file->varlist[varid].ndims;\n";
		  print F "    return PIO_NOERR;\n  }\n";
#	      }elsif($func =~ /inq_dimlen/){
#		  print F "  printf(\"dimid %d\\n\",dimid);\n";

	      }
	  }

	  if($line =~ /check_netcdf/){

	      if($func =~ /inq_varid/){
		  print F "    mpierr = MPI_Bcast(varidp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_ndims/){
		  print F "    mpierr = MPI_Bcast(ndimsp , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";	
	      }elsif($func =~ /def_var/){
		  print F "    mpierr = MPI_Bcast(varidp , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";	
	      }elsif($func =~ /def_dim/){
		  print F "    mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";	
	      }elsif($func =~ /inq_format/){
		  print F "    mpierr = MPI_Bcast(formatp , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";	
	      }elsif($func =~ /inq_natts/){
		  print F "    mpierr = MPI_Bcast(ngattsp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_varnatts/){
		  print F "    mpierr = MPI_Bcast(nattsp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_nvars/){
		  print F "    mpierr = MPI_Bcast(nvarsp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_varndims/){
		  print F "    mpierr = MPI_Bcast(ndimsp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "    file->varlist[varid].ndims = (*ndimsp);\n";
	      }elsif($func =~ /inq_dimlen/ || $func =~ /inq_attlen/){
		  print F "    mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_varid/){
		  print F "    mpierr = MPI_Bcast(varidp , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_vartype/ || $func =~ /inq_atttype/){
		  print F "    mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_attid/ || $func =~ /inq_dimid/){
		  print F "    mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_att$/){
		  print F "    if(xtypep != NULL) mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "    if(lenp != NULL) mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_var$/){
		  print F "    if(xtypep != NULL) mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "    if(ndimsp != NULL){ mpierr = MPI_Bcast(ndimsp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);\n";
		  print F "      file->varlist[varid].ndims = (*ndimsp);}\n";
		  print F "      if(nattsp != NULL) mpierr = MPI_Bcast(nattsp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "    if(name != NULL){ \n";
                  print F "       char tname[PIO_MAX_NAME];\n";
		  print F  "      if(ios->iomaster)\n";
   		  print F "	       strcpy(tname, name);\n";
		  print F "      mpierr = MPI_Bcast(tname , PIO_MAX_NAME, MPI_CHAR, ios->ioroot, ios->my_comm);\n";
		  print F "     strcpy(name,tname);\n";
	          print F "  }\n";
		  print F "    if(dimidsp != NULL) {int ndims;\n";
		  print F "      PIOc_inq_varndims(file->fh, varid, \&ndims);\n";
		  print F "      mpierr = MPI_Bcast(dimidsp , ndims, MPI_INT, ios->ioroot, ios->my_comm);\n ";
		  print F "    }\n";


		  
	      }elsif($func =~ /inq_dim$/){
		  print F  "    if(name != NULL){ char tname[PIO_MAX_NAME];\n";
		  print F  "      if(ios->iomaster)\n";
   		  print F "	       strcpy(tname, name);\n";
		  print F "      mpierr = MPI_Bcast(tname , PIO_MAX_NAME, MPI_CHAR, ios->ioroot, ios->my_comm);\n";
		  print F "      strcpy(name,tname); }\n";
		  print F "      if(lenp != NULL) mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_dimname/ || $func =~ /inq_varname/ || $func =~ /inq_attname/){
		  print F  "    { char tname[PIO_MAX_NAME];\n";
		  print F  "        if(ios->iomaster)\n";
   		  print F "	       strcpy(tname, name);\n";
		  print F "        mpierr = MPI_Bcast(tname , PIO_MAX_NAME, MPI_CHAR, ios->ioroot, ios->my_comm);\n";
		  print F "      strcpy(name,tname); }\n";
	      }elsif($func =~ /inq_vardimid/){
		  print F "    {int ndims;\n";
		  print F "    PIOc_inq_varndims(file->fh, varid, \&ndims);\n";
		  print F "    mpierr = MPI_Bcast(dimidsp , ndims, MPI_INT, ios->ioroot, ios->my_comm);\n ";
		  print F "    }\n";
	      }elsif($func =~ /get_att_(\w+)/){
		  my $atype = $1;
		  print F "      if(ierr == PIO_NOERR){\n";
		  print F "        PIO_Offset attlen;\n";		  
		  print F "        PIOc_inq_attlen(file->fh, varid, name, \&attlen);\n";
		  print F "        mpierr = MPI_Bcast(ip , (int) attlen, $typemap->{$atype}, ios->ioroot, ios->my_comm);\n ";
		  print F "     }\n";
	      }elsif($func =~/inq$/){
		  print F "      if(ndimsp != NULL)\n";
		  print F "        mpierr = MPI_Bcast(ndimsp, 1, MPI_INT, ios->ioroot, ios->my_comm);\n ";
		  print F "      if(nvarsp != NULL)\n";
		  print F "        mpierr = MPI_Bcast(nvarsp, 1, MPI_INT, ios->ioroot, ios->my_comm);\n ";
		  print F "      if(ngattsp != NULL)\n";
		  print F "        mpierr = MPI_Bcast(ngattsp, 1, MPI_INT, ios->ioroot, ios->my_comm);\n ";
		  print F "      if(unlimdimidp != NULL)\n";
		  print F "        mpierr = MPI_Bcast(unlimdimidp, 1, MPI_INT, ios->ioroot, ios->my_comm);\n ";

	      }

	  }
      }

#    print "$func  $functions->{$func}{pnetcdf} $functions->{$func}{netcdf}\n";
      print F "\n";
  }
} 
close(F);
