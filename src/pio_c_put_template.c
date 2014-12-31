///
/// PIO interface to nc_function
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIO_function()
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  int request;
  PIO_Offset usage;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = 0;

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
      //      ncmpi_begin_indep_data(file->fh);
      usage = 0;
      if(ios->io_rank==file->indep_rank){
	ierr = ncmpi_function();
	pio_push_request(file, request);
	ierr = ncmpi_inq_buffer_usage(ncid, &usage);
	//	printf("%s %d %d\n",__FILE__,__LINE__,usage);
      }
      // ncmpi_end_indep_data(file->fh);
      MPI_Bcast(&usage, 1,  MPI_LONG_LONG, file->indep_rank, ios->io_comm);
      file->indep_rank = (file->indep_rank + 1) % ios->num_iotasks;
      if(usage >= 0.8*PIO_BUFFER_SIZE_LIMIT){
	flush_output_buffer(file);
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
