int PIO_inq_att(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_att (file->fh,varid, name, xtypep, lenp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_att (file->fh,varid, name, xtypep, lenp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_att (file->fh,varid, name, xtypep, lenp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_format(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_format (file->fh, formatp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_format (file->fh, formatp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_format (file->fh, formatp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_varid(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varid (file->fh, name, varidp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varid (file->fh, name, varidp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varid (file->fh, name, varidp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_varnatts(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varnatts (file->fh,varid, nattsp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varnatts (file->fh,varid, nattsp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varnatts (file->fh,varid, nattsp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_def_var(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_def_var (file->fh, name, xtype,ndims,  dimidsp, varidp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_def_var (file->fh, name, xtype,ndims,  dimidsp, varidp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_def_var (file->fh, name, xtype, ndims, dimidsp, varidp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_var(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_var (file->fh,varid, name, xtypep,  ndimsp, dimidsp, nattsp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_var (file->fh,varid, name, xtypep,  ndimsp, dimidsp, nattsp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_var (file->fh,varid, name, xtypep, ndimsp, dimidsp, nattsp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_varname(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varname (file->fh,varid, name); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varname (file->fh,varid, name); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varname (file->fh,varid, name); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_double(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_double (file->fh,varid, name, xtype, size_t len, const double *op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_double (file->fh,varid, name, xtype, size_t len, const double *op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_double (file->fh,varid, name, xtype, MPI_Offset len, const double *op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_int(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_int (file->fh,varid, name, xtype, size_t len, op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_int (file->fh,varid, name, xtype, size_t len, op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_int (file->fh,varid, name, xtype, MPI_Offset len, op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_rename_att(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_rename_att (file->fh,varid, name, newname); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_att (file->fh,varid, name, newname); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_att (file->fh,varid, name, newname); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_del_att(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_del_att (file->fh,varid, name); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_del_att (file->fh,varid, name); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_del_att (file->fh,varid, name); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_natts(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_natts (file->fh, nattsp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_natts (file->fh, nattsp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_natts (file->fh, ngattsp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq (file->fh, ndimsp, nvarsp, nattsp, unlimdimidp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq (file->fh, ndimsp, nvarsp, nattsp, unlimdimidp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq (file->fh, ndimsp, nvarsp, ngattsp, unlimdimidp);  ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_text(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_text (file->fh,varid, name, ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_text (file->fh,varid, name, ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_text (file->fh,varid, name, ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_short(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_short (file->fh,varid, name, short *ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_short (file->fh,varid, name, short *ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_short (file->fh,varid, name, short *ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_libvers(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_libvers (void); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_libvers (void); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_libvers (void); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_long(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_long (file->fh,varid, name, xtype, size_t len, const long *op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_long (file->fh,varid, name, xtype, size_t len, const long *op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_long (file->fh,varid, name, xtype, MPI_Offset len, const long *op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_redef(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_redef (file->fh); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_redef (file->fh); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_redef (file->fh); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_strerror(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_strerror (int ncerr); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_strerror (int ncerr); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_strerror (int err); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_set_fill(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_set_fill (file->fh,fillmode, old_modep); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_set_fill (file->fh,fillmode, old_modep); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_set_fill (file->fh,fillmode, old_modep); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_delete(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_delete (path); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_delete (path); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_delete (filename, MPI_Info info); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_enddef(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_enddef (file->fh); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_enddef (file->fh); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_enddef (file->fh); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_rename_var(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_rename_var (file->fh,varid, name); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_var (file->fh,varid, name); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_var (file->fh,varid, name); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_short(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_short (file->fh,varid, name, xtype, size_t len, const short *op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_short (file->fh,varid, name, xtype, size_t len, const short *op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_short (file->fh,varid, name, xtype, MPI_Offset len, const short *op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_text(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_text (file->fh,varid, name, size_t len, op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_text (file->fh,varid, name, size_t len, op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_text (file->fh,varid, name, MPI_Offset len, op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_attname(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_attname (file->fh,varid,attnum, name); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attname (file->fh,varid,attnum, name); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attname (file->fh,varid,attnum, name); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_ulonglong(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_ulonglong (file->fh,varid, name,  unsigned long long *ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_ulonglong (file->fh,varid, name,  unsigned long long *ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_ulonglong (file->fh,varid, name, unsigned long long *ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_ushort(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_ushort (file->fh,varid, name, unsigned short *ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_ushort (file->fh,varid, name, unsigned short *ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_ushort (file->fh,varid, name, unsigned short *ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_ulonglong(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_ulonglong (file->fh,varid, name, xtype, size_t len, const unsigned long long *op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_ulonglong (file->fh,varid, name, xtype, size_t len, const unsigned long long *op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_ulonglong (file->fh,varid, name, xtype, MPI_Offset len, const unsigned long long *op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_dimlen(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dimlen (file->fh,dimid, lenp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimlen (file->fh,dimid, lenp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimlen (file->fh,dimid, lenp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_uint(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_uint (file->fh,varid, name, unsigned ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_uint (file->fh,varid, name, unsigned ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_uint (file->fh,varid, name, unsigned ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_longlong(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_longlong (file->fh,varid, name, long long *ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_longlong (file->fh,varid, name, long long *ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_longlong (file->fh,varid, name, long long *ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_schar(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_schar (file->fh,varid, name, xtype, size_t len, const signed op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_schar (file->fh,varid, name, xtype, size_t len, const signed op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_schar (file->fh,varid, name, xtype, MPI_Offset len, const signed op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_float(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_float (file->fh,varid, name, xtype, size_t len, const float *op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_float (file->fh,varid, name, xtype, size_t len, const float *op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_float (file->fh,varid, name, xtype, MPI_Offset len, const float *op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_nvars(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_nvars (file->fh, nvarsp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_nvars (file->fh, nvarsp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_nvars (file->fh, nvarsp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_rename_dim(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_rename_dim (file->fh,dimid, name); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_dim (file->fh,dimid, name); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_dim (file->fh,dimid, name); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_varndims(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_varndims (file->fh,varid, ndimsp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varndims (file->fh,varid, ndimsp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varndims (file->fh,varid, ndimsp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_long(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_long (file->fh,varid, name, long *ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_long (file->fh,varid, name, long *ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_long (file->fh,varid, name, long *ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_abort(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_abort (file->fh); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_abort (file->fh); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_abort (file->fh); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_dim(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dim (file->fh,dimid, name, lenp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dim (file->fh,dimid, name, lenp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dim (file->fh,dimid, name, lenp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_dimid(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dimid (file->fh, name, idp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimid (file->fh, name, idp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimid (file->fh, name, idp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_create(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_create (path,cmode, ncidp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_create (path,cmode, ncidp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_create (MPI_Comm comm, path,cmode, MPI_Info info, ncidp);  ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_set_default_format(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_set_default_format (int format, old_formatp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_set_default_format (int format, old_formatp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_set_default_format (int format, old_formatp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_unlimdim(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_unlimdim (file->fh, unlimdimidp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_unlimdim (file->fh, unlimdimidp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_unlimdim (file->fh, unlimdimidp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_vardimid(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_vardimid (file->fh,varid, dimidsp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_vardimid (file->fh,varid, dimidsp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_vardimid (file->fh,varid, dimidsp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_attlen(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_attlen (file->fh,varid, name, lenp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attlen (file->fh,varid, name, lenp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attlen (file->fh,varid, name, lenp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_dimname(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_dimname (file->fh,dimid, name); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimname (file->fh,dimid, name); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimname (file->fh,dimid, name); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_ushort(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_ushort (file->fh,varid, name, xtype, size_t len, const unsigned short *op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_ushort (file->fh,varid, name, xtype, size_t len, const unsigned short *op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_ushort (file->fh,varid, name, xtype, MPI_Offset len, const unsigned short *op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_float(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_float (file->fh,varid, name, float *ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_float (file->fh,varid, name, float *ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_float (file->fh,varid, name, float *ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_open(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_open (path,mode, ncidp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_open (path,mode, ncidp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_open (MPI_Comm comm, path,omode, MPI_Info info, ncidp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_sync(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_sync (file->fh); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_sync (file->fh); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_sync (file->fh); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_copy_att(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_copy_att (file->fh_in,varid_in, name,ncid_out,varid_out); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_copy_att (file->fh_in,varid_in, name,ncid_out,varid_out); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_copy_att (file->fh_in,varid_in, name,ncid_out,varid_out); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_longlong(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_longlong (file->fh,varid, name, xtype, size_t len, const long long *op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_longlong (file->fh,varid, name, xtype, size_t len, const long long *op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_longlong (file->fh,varid, name, xtype, MPI_Offset len, const long long *op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_uint(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_uint (file->fh,varid, name, xtype, size_t len, const unsigned op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_uint (file->fh,varid, name, xtype, size_t len, const unsigned op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_uint (file->fh,varid, name, xtype, MPI_Offset len, const unsigned op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_schar(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_schar (file->fh,varid, name, signed ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_schar (file->fh,varid, name, signed ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_schar (file->fh,varid, name, signed ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_attid(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_attid (file->fh,varid, name, idp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attid (file->fh,varid, name, idp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attid (file->fh,varid, name, idp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_close(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_close (file->fh); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_close (file->fh); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_close (file->fh); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_def_dim(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_def_dim (file->fh, name, size_t len, idp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_def_dim (file->fh, name, size_t len, idp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_def_dim (file->fh, name, MPI_Offset len, idp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_ndims(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_ndims (file->fh, ndimsp); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_ndims (file->fh, ndimsp); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_ndims (file->fh, ndimsp); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_vartype(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_vartype (file->fh,varid, xtypep); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_vartype (file->fh,varid, xtypep); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_vartype (file->fh,varid, xtypep); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_int(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_int (file->fh,varid, name, ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_int (file->fh,varid, name, ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_int (file->fh,varid, name, ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_double(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_double (file->fh,varid, name, double *ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_double (file->fh,varid, name, double *ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_double (file->fh,varid, name, double *ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_inq_atttype(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_inq_atttype (file->fh,varid, name, xtypep); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_atttype (file->fh,varid, name, xtypep); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_atttype (file->fh,varid, name, xtypep); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_put_att_uchar(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_put_att_uchar (file->fh,varid, name, xtype, size_t len, const unsigned op); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_uchar (file->fh,varid, name, xtype, size_t len, const unsigned op); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_uchar (file->fh,varid, name, xtype, MPI_Offset len, const unsigned op); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
int PIO_get_att_uchar(struct file_desc_t *file)
{
  int ierr;
  int msg;
  int mpierr;
  struct iosystem_desc_t *ios;
  ierr = PIO_NOERR;

  ios = file->iosystem;

  if(ios->async_interface && ! ios->ioproc){
    if(iosys->comp_rank==0) 
      mpierr = mpi_send(msg, 1 mpi_integer, ios%ioroot, 1, ios%union_comm);
    mpierr = mpi_bcast(file->fh,1, mpi_integer, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_att_uchar (file->fh,varid, name, unsigned ip); ;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_uchar (file->fh,varid, name, unsigned ip); ;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_uchar (file->fh,varid, name, unsigned ip); ;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
