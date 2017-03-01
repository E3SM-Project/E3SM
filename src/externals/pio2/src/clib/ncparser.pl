#!/usr/bin/perl 
use strict;

my $netcdf_incdir = $ARGV[0];
my $pnetcdf_incdir = $ARGV[1];

my $netcdf_inc = "$netcdf_incdir/netcdf.h";
my $pnetcdf_inc = "$pnetcdf_incdir/pnetcdf.h";

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
* \@file   
* \@author Jim Edwards (jedwards\@ucar.edu)
* \@date     Feburary 2014 
* \@brief    PIO interfaces to [NetCDF](http://www.unidata.ucar.edu/software/netcdf/docs/modules.html) support functions
* \@details
*  This file provides an interface to the [NetCDF](http://www.unidata.ucar.edu/software/netcdf/docs/modules.html) support functions.
*  Each subroutine calls the underlying netcdf or pnetcdf or netcdf4 functions from 
*  the appropriate subset of mpi tasks (io_comm). Each routine must be called 
*  collectively from union_comm.
*  
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

my $urls    = {files=>"http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html",
               dims=>"http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html",
               vars=>"http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html",
               atts=>"http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html",
               groups=>"http://www.unidata.ucar.edu/software/netcdf/docs/group__groups.html",
               utypes=>"http://www.unidata.ucar.edu/software/netcdf/docs/group__user__types.html",
               other=>"http://www.unidata.ucar.edu/software/netcdf/docs/modules.html"};
my $currurl = $urls->{other};


foreach my $func (keys %{$functions}){
  next if($func=~/get_v/ or $func=~/put_v/);
  next if($func=~/strerror/);
  next if($func=~/delete/);
  next if($func=~/sync/);
  next if($func=~/open/);
  next if($func=~/close/);
  next if($func=~/create/);
  next if($func=~/default_format/);
  next if($func=~/copy_att/);
  next if($func=~/abort/);
  next if($func=~/def_var_fill/);
  next if($func=~/string/);
  next if($func=~/t_att$/); # skip void versions of get and put att
  next if($functions->{$func}{pnetcdf} =~ /\(void\)/);

  my @tmptmp = @template;

  if($functions->{$func}{pnetcdf} =~ /ngattsp/){
      $functions->{$func}{netcdf} =~ s/nattsp/ngattsp/g;
  }

  if($functions->{$func}{pnetcdf} =~ /fill_value/){
      $functions->{$func}{netcdf} =~ s/fill_valuep/fill_value/g;
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

		  if($func=~/.*var/){
		      $currurl = $urls->{vars};
		  } elsif ($func=~/.*att/) {
		      $currurl = $urls->{atts};
		  } elsif ($func=~/.*dim/){
		      $currurl = $urls->{dims};
		  } elsif ($func=~/.*grp/){
		      $currurl = $urls->{groups};
		  } elsif ($func=~/.*def/ || $func=~/.*fill/ || $func=~/.*inq/) {
		      $currurl = $urls->{files};
		  } else {
		      $currurl = $urls->{other}; }

		  print F 
"/** 
 * \@ingroup PIOc_$func
 * The PIO-C interface for the NetCDF function nc_$func.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * ".$currurl."
 *
 * \@param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().\n";
		  if($functions->{$func}{pnetcdf} =~ /varid/){
		      print F " * \@param varid the variable ID.\n";
		  }

		  if($functions->{$func}{pnetcdf} =~ /int attnum/){
		      print F " * \@param attnum the attribute ID.\n";
		  }

		  if($functions->{$func}{pnetcdf} =~ /const char *name/){
		      print F " * \@param name the attribute name.\n";
		  }

		  if($functions->{$func}{pnetcdf} =~ /\*xtypep/){
		      print F " * \@param xtypep a pointer that will get the type of the attribute.\n";
		  }
		  if($functions->{$func}{pnetcdf} =~ /\*idp/){
		      print F " * \@param idp a pointer that will get the id of the variable or attribute.\n";
		  }
		  if($functions->{$func}{pnetcdf} =~ /\*lenp/){
		      print F " * \@param lenp a pointer that will get the number of values \n";
		  }
		  if($functions->{$func}{pnetcdf} =~ /\*varidp/){
		      print F " * \@param varidp a pointer that will get the variable id \n";
		  }
		  if($functions->{$func}{pnetcdf} =~ /\*formatp/){
		      print F " * \@param formatp a pointer that will get the file format \n";
		  }
		  if($functions->{$func}{pnetcdf} =~ /\*nattsp/){
		      print F " * \@param nattsp a pointer that will get the number of attributes \n";
		  }
		  
print F 
" * \@return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
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
	  if($line =~ /check_netcdf/){
	      print F "   if(ierr != PIO_NOERR){\n";
	      if($func =~ /get_att/){
		  print F "    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));\n";
		  print F "    sprintf(errstr,\"name %s in file %s\",name,__FILE__);\n";
	      }else{
		  print F "    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));\n";
		  print F "    sprintf(errstr,\"in file %s\",__FILE__);\n";
	      }	      
              print F "  }\n";
	  }


	  print F $line;
	  if($func =~ /sync/ && $line =~   /case PIO_IOTYPE_PNETCDF/){
	      printf F "      flush_output_buffer(file);\n";
	  }
	  if($func =~ /def_var/ && $line =~ /NETCDF4C/){
	      printf F "      if(ios->io_rank==0){\n";
	      printf F "        ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);\n";
	      printf F "        if(ierr == PIO_NOERR){\n";
	      printf F "          ierr = nc_def_var_deflate(file->fh, *varidp, 0,1,1);\n";
	      printf F "        }\n";  
	      printf F "      }\n";  
	      printf F "      break;\n";  
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
	      }elsif($func =~ /inq_unlimdim/){
		  print F "    mpierr = MPI_Bcast(unlimdimidp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_ndims/){
		  print F "    mpierr = MPI_Bcast(ndimsp , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";	
	      }elsif($func =~ /var_fill/){
                  print F "    mpierr = MPI_Bcast(fill_value, 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
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
		  print F "  if(xtypep != NULL) mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "  if(lenp != NULL) mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_var$/){
		  print F "    if(xtypep != NULL) mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "    if(ndimsp != NULL){ mpierr = MPI_Bcast(ndimsp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);\n";
		  print F "      file->varlist[varid].ndims = (*ndimsp);}\n";
		  print F "      if(nattsp != NULL) mpierr = MPI_Bcast(nattsp,1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "    if(name != NULL){ \n";
		  print F "      int slen;\n";
		  print F "      if(ios->iomaster)\n";
		  print F "        slen = (int) strlen(name) + 1;\n";
		  print F "      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);\n    }\n";

		  print F "    if(dimidsp != NULL) {int ndims;\n";
		  print F "      PIOc_inq_varndims(file->fh, varid, \&ndims);\n";
		  print F "      mpierr = MPI_Bcast(dimidsp , ndims, MPI_INT, ios->ioroot, ios->my_comm);\n ";
		  print F "    }\n";


		  
	      }elsif($func =~ /inq_dim$/){
		  print F  "    if(name != NULL){ \n";
		  print F "      int slen;\n";
		  print F "      if(ios->iomaster)\n";
		  print F "        slen = (int) strlen(name) + 1;\n";
		  print F "      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);\n    }\n";
		  print F "      if(lenp != NULL) mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);\n";
	      }elsif($func =~ /inq_dimname/ || $func =~ /inq_varname/ || $func =~ /inq_attname/){
		  print F "    if(name != NULL){\n";
		  print F "      int slen;\n";
		  print F "      if(ios->iomaster)\n";
		  print F "        slen = (int) strlen(name) + 1;\n";
		  print F "      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);\n";
		  print F "      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);\n    }\n";
	      }elsif($func =~ /inq_vardimid/){
		  print F "    if(ierr==PIO_NOERR){\n";
                  print F "      int ndims;\n";
		  print F "      PIOc_inq_varndims(file->fh, varid, \&ndims);\n";
		  print F "      mpierr = MPI_Bcast(dimidsp , ndims, MPI_INT, ios->ioroot, ios->my_comm);\n ";
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
