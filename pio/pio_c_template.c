int PIO_function(struct file_desc_t *file)
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
      ierr = nc_function();
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_function();
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_function();
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


}
