#include <pio.h>
#include <pio_internal.h>

int PIOc_put_vars_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_uchar(file->fh, varid, start, count, stride, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_schar (int ncid, int varid, const PIO_Offset index[], signed char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_schar(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_schar(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_schar(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_ulonglong(file->fh, varid, start, count, stride, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned char *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_uchar(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], signed char *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_schar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_schar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_schar(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], short *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_short(file->fh, varid, start, count, stride, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_double (int ncid, int varid, double *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_double(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_double(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_double(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], double *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_double(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_double(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_double(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_int (int ncid, int varid, int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_int(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_int(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_int(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_ushort (int ncid, int varid, unsigned short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_ushort(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_ushort(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_ushort(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_ushort(file->fh, varid, start, count, stride, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], char *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_text(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_text(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_text(file->fh, varid, start, count, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_ulonglong(file->fh, varid, start, count, stride, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], int *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_int(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_int(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_int(file->fh, varid, start, count, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm(file->fh, varid, start, count, stride, imap, buf, bufcount, buftype); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_float (int ncid, int varid, const PIO_Offset index[], float *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_float(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_float(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_float(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_short (int ncid, int varid, const PIO_Offset index[], short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_short(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_short(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_short(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], int *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_int(file->fh, varid, start, count, stride, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_uint(file->fh, varid, start, count, stride, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_text (int ncid, int varid, char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_text(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_text(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_text(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], double *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_double(file->fh, varid, start, count, stride, imap, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_uchar(file->fh, varid, start, count, stride, imap, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_ushort (int ncid, int varid, const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_ushort(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_ushort(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_ushort(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], signed char *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_schar(file->fh, varid, start, count, stride, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_ushort(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_longlong (int ncid, int varid, const PIO_Offset index[], const long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_longlong(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_longlong(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_longlong(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_uchar(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const short *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_short(file->fh, varid, start, count, stride, imap, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_long (int ncid, int varid, const PIO_Offset index[], const long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_long(file->fh, varid, (size_t *) index,  ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_long(file->fh, varid, (size_t *) index,  ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_long(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_long(file->fh, varid, start, count, stride, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_ushort (int ncid, int varid, const PIO_Offset index[], unsigned short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_ushort(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_ushort(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_ushort(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_short (int ncid, int varid, const short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_short(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_short(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_short(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_int(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_int(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_int(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_float (int ncid, int varid, float *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_float(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_float(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_float(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_ushort (int ncid, int varid, const PIO_Offset index[], const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_ushort(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_ushort(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_ushort(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_text(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_text(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_text(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const char *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_text(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_text(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_text(file->fh, varid, start, count, stride, imap, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_uchar(file->fh, varid, start, count, stride, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var (int ncid, int varid, void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var(file->fh, varid,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var(file->fh, varid,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var(file->fh, varid, buf, bufcount, buftype);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_ushort(file->fh, varid, start, count, stride, imap, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_longlong (int ncid, int varid, const PIO_Offset index[], long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_longlong(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_longlong(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_longlong(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_ushort(file->fh, varid, start, count, stride, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_long (int ncid, int varid, long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_long(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_long(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_long(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_double (int ncid, int varid, const PIO_Offset index[], double *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_double(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_double(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_double(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_ulonglong (int ncid, int varid, const unsigned long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_ulonglong(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_ulonglong(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_ulonglong(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_int (int ncid, int varid, const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_int(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_int(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_int(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_uint(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_uint(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_uint(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_longlong (int ncid, int varid, const long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_longlong(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_longlong(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_longlong(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_longlong(file->fh, varid, start, count, stride, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_schar (int ncid, int varid, const signed char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_schar(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_schar(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_schar(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_uint (int ncid, int varid, const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_uint(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_uint(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_uint(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var (int ncid, int varid, const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var(file->fh, varid,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var(file->fh, varid,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var(file->fh, varid, buf, bufcount, buftype);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_ushort(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_longlong (int ncid, int varid, long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_longlong(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_longlong(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_longlong(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], short *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_short(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_short(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_short(file->fh, varid, start, count, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const short *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_short(file->fh, varid, start, count, stride, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_uint(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_uint(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_uint(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const signed char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_schar(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_schar(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_schar(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_ulonglong(file->fh, varid, start, count, stride, imap, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_uchar (int ncid, int varid, const PIO_Offset index[], const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_uchar(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_uchar(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_uchar(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_int(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_int(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_int(file->fh, varid, start, count, stride, imap, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const signed char *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_schar(file->fh, varid, start, count, stride, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_long(file->fh, varid, (size_t *) start, (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_long(file->fh, varid, (size_t *) start, (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_long(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1 (int ncid, int varid, const PIO_Offset index[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1(file->fh, varid, (size_t *) index,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1(file->fh, varid, (size_t *) index,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1(file->fh, varid, index, buf, bufcount, buftype);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_int (int ncid, int varid, const PIO_Offset index[], int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_int(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_int(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_int(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_ulonglong (int ncid, int varid, const PIO_Offset index[], unsigned long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_ulonglong(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_ulonglong(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_ulonglong(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_uchar (int ncid, int varid, unsigned char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_uchar(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_uchar(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_uchar(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const float *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_float(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_float(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_float(file->fh, varid, start, count, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_uchar(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], float *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_float(file->fh, varid, start, count, stride, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_float (int ncid, int varid, const PIO_Offset index[], const float *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_float(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_float(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_float(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const float *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_float(file->fh, varid,(size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_float(file->fh, varid,(size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_float(file->fh, varid, start, count, stride, imap, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_text (int ncid, int varid, const PIO_Offset index[], const char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_text(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_text(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_text(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const char *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_text(file->fh, varid, start, count, stride, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_long(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_long(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_long(file->fh, varid, start, count, stride, imap, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_long(file->fh, varid, start, count, stride, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const double *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_double(file->fh, varid, start, count, stride, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1 (int ncid, int varid, const PIO_Offset index[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1(file->fh, varid, (size_t *) index,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1(file->fh, varid, (size_t *) index,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1(file->fh, varid, index, buf, bufcount, buftype);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_uint (int ncid, int varid, unsigned int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_uint(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_uint(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_uint(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const long long *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_longlong(file->fh, varid, start, count, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara(file->fh, varid, (size_t *) start,  (size_t *) count,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara(file->fh, varid, (size_t *) start,  (size_t *) count,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara(file->fh, varid, start, count, buf, bufcount, buftype);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_double (int ncid, int varid, const double *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_double(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_double(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_double(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], signed char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_schar(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_schar(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_schar(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_float (int ncid, int varid, const float *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_float(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_float(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_float(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_uint (int ncid, int varid, const PIO_Offset index[], unsigned int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_uint(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_uint(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_uint(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_uint(file->fh, varid, start, count, stride, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_ulonglong (int ncid, int varid, const PIO_Offset index[], const unsigned long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_ulonglong(file->fh, varid, (size_t *) index,   ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_ulonglong(file->fh, varid, (size_t *) index,   ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_ulonglong(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_uint(file->fh, varid, start, count, stride, imap, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_uint (int ncid, int varid, const PIO_Offset index[], const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_uint(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_uint(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_uint(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_int (int ncid, int varid, const PIO_Offset index[], const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_int(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_int(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_int(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], float *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_float(file->fh, varid, (size_t *) start, (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_float(file->fh, varid, (size_t *) start, (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_float(file->fh, varid, start, count, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], char *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_text(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_text(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_text(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const float *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_float(file->fh, varid, start, count, stride, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_text (int ncid, int varid, const PIO_Offset index[], char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_text(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_text(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_text(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_short(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_short(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_short(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_schar (int ncid, int varid, const PIO_Offset index[], const signed char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_schar(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_schar(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_schar(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_ulonglong(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const double *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_double(file->fh, varid, start, count, stride, imap, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], int *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_int(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_int(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_int(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara(file->fh, varid, (size_t *) start,  (size_t *) count,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara(file->fh, varid, (size_t *) start,  (size_t *) count,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara(file->fh, varid, start, count, buf, bufcount, buftype);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_long(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_long(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_long(file->fh, varid, start, count, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned int *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_uint(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm(file->fh, varid, start, count, stride, imap, buf, bufcount, buftype); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_double (int ncid, int varid, const PIO_Offset index[], const double *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_double(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_double(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_double(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], double *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_double(file->fh, varid, start, count, stride, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_longlong(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_ulonglong (int ncid, int varid, unsigned long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_ulonglong(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_ulonglong(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_ulonglong(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const signed char *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_schar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_schar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_schar(file->fh, varid, start, count, stride, imap, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vara_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARA_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vara_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vara_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vara_ulonglong(file->fh, varid, start, count, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_short (int ncid, int varid, short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_short(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_short(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_short(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], float *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_float(file->fh, varid,(size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_float(file->fh, varid,(size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_float(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_text (int ncid, int varid, const char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_text(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_text(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_text(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const int *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_int(file->fh, varid, start, count, stride, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_long (int ncid, int varid, const PIO_Offset index[], long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_long(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_long(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_long(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_long(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_long(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_long(file->fh, varid, start, count, stride, imap, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned short *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_ushort(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var1_short (int ncid, int varid, const PIO_Offset index[], const short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR1_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_short(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_short(file->fh, varid, (size_t *) index, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var1_short(file->fh, varid, index, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const long long *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars_longlong(file->fh, varid, start, count, stride, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_longlong(file->fh, varid, start, count, stride, imap, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], char *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars_text(file->fh, varid, start, count, stride, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vara_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const double *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_double(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_double(file->fh, varid, (size_t *) start,  (size_t *) count, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vara_double(file->fh, varid, start, count, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_vars (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARS;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_vars(file->fh, varid, start, count, stride, buf, bufcount, buftype); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_uchar (int ncid, int varid, const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_uchar(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_uchar(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_uchar(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var1_uchar (int ncid, int varid, const PIO_Offset index[], unsigned char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR1_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var1_uchar(file->fh, varid, (size_t *) index, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var1_uchar(file->fh, varid, (size_t *) index, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var1_uchar(file->fh, varid, index, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_var_long (int ncid, int varid, const long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VAR_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var_long(file->fh, varid, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var_long(file->fh, varid, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_var_long(file->fh, varid, op);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_vars (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARS;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_vars(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_vars(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_vars(file->fh, varid, start, count, stride, buf, bufcount, buftype); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], short *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_short(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned long long *ip)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_varm_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_varm_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_varm_ulonglong(file->fh, varid, start, count, stride, imap, ip); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_put_varm_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long long *op)  
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_put_varm_longlong(file->fh, varid, start, count, stride, imap, op); ;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_var_schar (int ncid, int varid, signed char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VAR_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_get_var_schar(file->fh, varid, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_var_schar(file->fh, varid, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_begin_indep_data(file->fh);
      if(ios->io_rank==0 &&  (ierr == PIO_NOERR || ierr == PIO_EINDEP)){
	ierr = ncmpi_get_var_schar(file->fh, varid, ip);;
      }
      if(ierr == PIO_NOERR){
	ierr = ncmpi_end_indep_data(file->fh);
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

