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
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request,
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_function();
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
