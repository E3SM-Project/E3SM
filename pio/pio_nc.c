#include <mpi.h>
#include <pio.h>
#ifdef _NETCDF
#include <netcdf.h>
#endif
#ifdef _PNETCDF
#include <pnetcdf.h>
#endif
int PIO_inq_att (struct file_desc_t *file, int varid, const char *name, nc_type *xtypep, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_att(file->fh, varid, name, xtypep, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_att(file->fh, varid, name, xtypep, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_att(file->fh, varid, name, xtypep, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_format (struct file_desc_t *file, int *formatp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_format(file->fh, formatp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_format(file->fh, formatp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_format(file->fh, formatp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_varid (struct file_desc_t *file, const char *name, int *varidp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varid(file->fh, name, varidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varid(file->fh, name, varidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varid(file->fh, name, varidp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_varnatts (struct file_desc_t *file, int varid, int *nattsp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varnatts(file->fh, varid, nattsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varnatts(file->fh, varid, nattsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varnatts(file->fh, varid, nattsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_def_var (struct file_desc_t *file, const char *name, nc_type xtype,  int ndims, const int *dimidsp, int *varidp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_def_var(file->fh, name, xtype,  ndims, dimidsp, varidp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_var (struct file_desc_t *file, int varid, char *name, nc_type *xtypep, int *ndimsp, int *dimidsp, int *nattsp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_varname (struct file_desc_t *file, int varid, char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varname(file->fh, varid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varname(file->fh, varid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varname(file->fh, varid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_double (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const double *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_double(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_double(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_double(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_int (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_int(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_int(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_int(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_rename_att (struct file_desc_t *file, int varid, const char *name, const char *newname) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_rename_att(file->fh, varid, name, newname);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_att(file->fh, varid, name, newname);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_att(file->fh, varid, name, newname);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_del_att (struct file_desc_t *file, int varid, const char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_del_att(file->fh, varid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_del_att(file->fh, varid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_del_att(file->fh, varid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_natts (struct file_desc_t *file, int *ngattsp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_natts(file->fh, ngattsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_natts(file->fh, ngattsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_natts(file->fh, ngattsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq (struct file_desc_t *file, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp)  
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq(file->fh, ndimsp, nvarsp, ngattsp, unlimdimidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq(file->fh, ndimsp, nvarsp, ngattsp, unlimdimidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq(file->fh, ndimsp, nvarsp, ngattsp, unlimdimidp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_text (struct file_desc_t *file, int varid, const char *name, char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_text(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_text(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_text(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_short (struct file_desc_t *file, int varid, const char *name, short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_short(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_short(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_short(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_long (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_long(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_long(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_long(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_redef (struct file_desc_t *file) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_redef(file->fh);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_redef(file->fh);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_redef(file->fh);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_set_fill (struct file_desc_t *file, int fillmode, int *old_modep) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_set_fill(file->fh, fillmode, old_modep);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_set_fill(file->fh, fillmode, old_modep);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_set_fill(file->fh, fillmode, old_modep);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_enddef (struct file_desc_t *file) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_enddef(file->fh);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_enddef(file->fh);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_enddef(file->fh);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_rename_var (struct file_desc_t *file, int varid, const char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_rename_var(file->fh, varid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_var(file->fh, varid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_var(file->fh, varid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_short (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_short(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_short(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_short(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_text (struct file_desc_t *file, int varid, const char *name, PIO_Offset len, const char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_text(file->fh, varid, name, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_text(file->fh, varid, name, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_text(file->fh, varid, name, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_attname (struct file_desc_t *file, int varid, int attnum, char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_attname(file->fh, varid, attnum, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attname(file->fh, varid, attnum, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attname(file->fh, varid, attnum, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_ulonglong (struct file_desc_t *file, int varid, const char *name, unsigned long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_ulonglong(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_ulonglong(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_ulonglong(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_ushort (struct file_desc_t *file, int varid, const char *name, unsigned short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_ushort(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_ushort(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_ushort(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_ulonglong (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_ulonglong(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_ulonglong(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_ulonglong(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_dimlen (struct file_desc_t *file, int dimid, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dimlen(file->fh, dimid, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimlen(file->fh, dimid, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimlen(file->fh, dimid, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_uint (struct file_desc_t *file, int varid, const char *name, unsigned int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_uint(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_uint(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_uint(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_longlong (struct file_desc_t *file, int varid, const char *name, long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_longlong(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_longlong(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_longlong(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_schar (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const signed char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_schar(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_schar(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_schar(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_float (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const float *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_float(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_float(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_float(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_nvars (struct file_desc_t *file, int *nvarsp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_nvars(file->fh, nvarsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_nvars(file->fh, nvarsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_nvars(file->fh, nvarsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_rename_dim (struct file_desc_t *file, int dimid, const char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_rename_dim(file->fh, dimid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_dim(file->fh, dimid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_dim(file->fh, dimid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_varndims (struct file_desc_t *file, int varid, int *ndimsp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varndims(file->fh, varid, ndimsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varndims(file->fh, varid, ndimsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varndims(file->fh, varid, ndimsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_long (struct file_desc_t *file, int varid, const char *name, long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_long(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_long(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_long(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_dim (struct file_desc_t *file, int dimid, char *name, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dim(file->fh, dimid, name, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dim(file->fh, dimid, name, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dim(file->fh, dimid, name, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_dimid (struct file_desc_t *file, const char *name, int *idp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dimid(file->fh, name, idp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimid(file->fh, name, idp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimid(file->fh, name, idp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_unlimdim (struct file_desc_t *file, int *unlimdimidp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_unlimdim(file->fh, unlimdimidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_unlimdim(file->fh, unlimdimidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_unlimdim(file->fh, unlimdimidp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_vardimid (struct file_desc_t *file, int varid, int *dimidsp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_vardimid(file->fh, varid, dimidsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_vardimid(file->fh, varid, dimidsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_vardimid(file->fh, varid, dimidsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_attlen (struct file_desc_t *file, int varid, const char *name, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_attlen(file->fh, varid, name, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attlen(file->fh, varid, name, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attlen(file->fh, varid, name, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_dimname (struct file_desc_t *file, int dimid, char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dimname(file->fh, dimid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimname(file->fh, dimid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimname(file->fh, dimid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_ushort (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_ushort(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_ushort(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_ushort(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_float (struct file_desc_t *file, int varid, const char *name, float *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_float(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_float(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_float(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_sync (struct file_desc_t *file) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_sync(file->fh);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_sync(file->fh);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_sync(file->fh);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_longlong (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_longlong(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_longlong(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_longlong(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_uint (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_uint(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_uint(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_uint(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_schar (struct file_desc_t *file, int varid, const char *name, signed char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_schar(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_schar(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_schar(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_attid (struct file_desc_t *file, int varid, const char *name, int *idp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_attid(file->fh, varid, name, idp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attid(file->fh, varid, name, idp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attid(file->fh, varid, name, idp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_def_dim (struct file_desc_t *file, const char *name, PIO_Offset len, int *idp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_def_dim(file->fh, name, (size_t)len, idp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_def_dim(file->fh, name, (size_t)len, idp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_def_dim(file->fh, name, len, idp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_ndims (struct file_desc_t *file, int *ndimsp) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_ndims(file->fh, ndimsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_ndims(file->fh, ndimsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_ndims(file->fh, ndimsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_vartype (struct file_desc_t *file, int varid, nc_type *xtypep) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_vartype(file->fh, varid, xtypep);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_vartype(file->fh, varid, xtypep);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_vartype(file->fh, varid, xtypep);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_int (struct file_desc_t *file, int varid, const char *name, int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_int(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_int(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_int(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_double (struct file_desc_t *file, int varid, const char *name, double *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_double(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_double(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_double(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_inq_atttype (struct file_desc_t *file, int varid, const char *name, nc_type *xtypep) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_atttype(file->fh, varid, name, xtypep);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_atttype(file->fh, varid, name, xtypep);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_atttype(file->fh, varid, name, xtypep);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_put_att_uchar (struct file_desc_t *file, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_uchar(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_uchar(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_uchar(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
int PIO_get_att_uchar (struct file_desc_t *file, int varid, const char *name, unsigned char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = mpi_send(msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = mpi_bcast(file->fh,1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_uchar(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_uchar(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_uchar(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
