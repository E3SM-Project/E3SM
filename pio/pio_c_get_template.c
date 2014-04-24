int PIO_function()
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  size_t ibufcnt;


  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = 0;
  ierr = PIO_NOERR;

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
      ierr = nc_function();
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_function();
	if(! (ios->async_interface) &&  (ios->num_comptasks == ios->num_iotasks)){
	  MPI_Bcast(buf, ibufcnt, ibuftype, 0, ios->io_comm);
	}

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

  if(ios->async_interface ||  (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}
