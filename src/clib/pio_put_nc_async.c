#include <pio.h>
#include <pio_internal.h>

/**
 * Internal PIO function which provides a type-neutral interface to
 * nc_put_vars
 *
 * This routine is called collectively by all tasks in the
 * communicator ios.union_comm.
 *
 * @param ncid identifies the netCDF file
 * @param varid the variable ID number

 * @param start[] an array of start indicies (must have same number of
 * entries as variable has dimensions). If NULL, indices of 0 will be
 * used.
 *
 * @param count[] an array of counts (must have same number of entries
 * as variable has dimensions). If NULL, counts matching the size of
 * the variable will be used.
 *
 * @param stride[] an array of strides (must have same number of
 * entries as variable has dimensions). If NULL, strides of 1 will be
 * used.
 * 
 * @param xtype the netCDF type of the data being passed in buf. Data
 * will be automatically covnerted from this type to the type of the
 * variable being written to.
 *
 * @param buf pointer to the data to be written.
 *
 * @return PIO_NOERR on success, error code otherwise.
 */
int PIOc_put_vars_tc(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		     const PIO_Offset *stride, nc_type xtype, const void *buf) 
{
    iosystem_desc_t *ios;  /** Pointer to io system information. */
    file_desc_t *file;     /** Pointer to file information. */
    int ierr = PIO_NOERR;  /** Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    int ndims; /** The number of dimensions in the variable. */
    int *dimids; /** The IDs of the dimensions for this variable. */
    PIO_Offset typelen; /** Size (in bytes) of the data type of data in buf. */
    size_t num_elem = 1; /** Number of data elements in the buffer. */
    var_desc_t *vdesc;
    PIO_Offset usage;
    int *request;

    LOG((1, "PIOc_put_vars_tc ncid = %d varid = %d start = %d count = %d "
	 "stride = %d xtype = %d", ncid, start, count, stride, xtype));

    /* User must provide some data. */
    if (!buf)
	return PIO_EINVAL;

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* Run these on all tasks if async is not in use, but only on
     * non-IO tasks if async is in use. */
    if (!ios->async_interface || !ios->ioproc)
    {
	/* Get the length of the data type. */
	if ((ierr = PIOc_inq_type(ncid, xtype, NULL, &typelen)))
	    return check_netcdf(file, ierr, __FILE__, __LINE__);

	/* Get the number of dims for this var. */
	if ((ierr = PIOc_inq_varndims(ncid, varid, &ndims)))
	    return check_netcdf(file, ierr, __FILE__, __LINE__);

	PIO_Offset dimlen[ndims];

	/* If no count array was passed, we need to know the dimlens
	 * so we can calculate how many data elements are in the
	 * buf. */
	if (!count)
	{
	    int dimid[ndims];
	    
	    /* Get the dimids for this var. */
	    if ((ierr = PIOc_inq_vardimid(ncid, varid, dimid)))
		return check_netcdf(file, ierr, __FILE__, __LINE__);
	    
	    /* Get the length of each dimension. */
	    for (int vd = 0; vd < ndims; vd++)
		if ((ierr = PIOc_inq_dimlen(ncid, dimid[vd], &dimlen[vd])))
		    return check_netcdf(file, ierr, __FILE__, __LINE__);
	}

	/* Figure out the real start, count, and stride arrays. (The
	 * user may have passed in NULLs.) */
	PIO_Offset rstart[ndims], rcount[ndims], rstride[ndims];
	for (int vd = 0; vd < ndims; vd++)
	{
	    rstart[vd] = start ? start[vd] : 0;
	    rcount[vd] = count ? count[vd] : dimlen[vd];
	    rstride[vd] = stride ? stride[vd] : 1;
	}

	/* How many elements in buf? */
	for (int vd = 0; vd < ndims; vd++)
	    num_elem *= (rcount[vd] - rstart[vd])/rstride[vd];
	LOG((2, "PIOc_put_vars_tc num_elem = %d", num_elem));
    }

    sleep(2);
    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_PUT_VARS;
	    char start_present = start ? true : false;
	    char count_present = count ? true : false;
	    char stride_present = stride ? true : false;

	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

	    if (!mpierr)
		mpierr = MPI_Bcast(&ncid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&start_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr && start_present)
		mpierr = MPI_Bcast(&start, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);		
	    if (!mpierr)
		mpierr = MPI_Bcast(&count_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr && count_present)
		mpierr = MPI_Bcast(&count, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);		
	    if (!mpierr)
		mpierr = MPI_Bcast(&stride_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr && stride_present)
		mpierr = MPI_Bcast(&stride, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);		
	    if (!mpierr)
		mpierr = MPI_Bcast(&xtype, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&num_elem, 1, MPI_OFFSET, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, ios->compmaster, ios->intercomm);
	    LOG((2, "PIOc_put_vars_tc ncid = %d varid = %d ndims = %d start_present = %d "
		 "count_present = %d stride_present = %d xtype = %d num_elem = %d", ncid, varid,
		 ndims, start_present, count_present, stride_present, xtype, num_elem));
	    
	    /* Send the data. */
	    if (!mpierr)
		mpierr = MPI_Bcast(buf, num_elem * typelen, MPI_BYTE, ios->compmaster,
				   ios->intercomm);
	}

	sleep(2);
	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    return check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);

	/* /\* Broadcast values currently only known on computation tasks to IO tasks. *\/ */
	/* if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->comproot, ios->my_comm))) */
	/*     check_mpi(file, mpierr, __FILE__, __LINE__); */
	/* if ((mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, ios->comproot, ios->my_comm))) */
	/*     check_mpi(file, mpierr, __FILE__, __LINE__); */
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
	if (file->iotype == PIO_IOTYPE_PNETCDF)
	{
	    vdesc = file->varlist + varid;
	    if (vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0)
		vdesc->request = realloc(vdesc->request, 
					 sizeof(int) * (vdesc->nreqs + PIO_REQUEST_ALLOC_CHUNK));
	    request = vdesc->request+vdesc->nreqs;

	    if(ios->io_rank == 0)
		switch(xtype)
		{
		case NC_BYTE:
		    ierr = ncmpi_bput_vars_schar(ncid, varid, start, count, stride, buf, request);
		    break;
		case NC_CHAR:
		    ierr = ncmpi_bput_vars_text(ncid, varid, start, count, stride, buf, request);
		    break;
		case NC_SHORT:
		    ierr = ncmpi_bput_vars_short(ncid, varid, start, count, stride, buf, request);
		    break;
		case NC_INT:
		    ierr = ncmpi_bput_vars_int(ncid, varid, start, count, stride, buf, request);
		    break;
		case NC_FLOAT:
		    ierr = ncmpi_bput_vars_float(ncid, varid, start, count, stride, buf, request);
		    break;
		case NC_DOUBLE:
		    ierr = ncmpi_bput_vars_double(ncid, varid, start, count, stride, buf, request);
		    break;
		case NC_INT64:
		    ierr = ncmpi_bput_vars_longlong(ncid, varid, start, count, stride, buf, request);
		    break;
		default:
		    LOG((0, "Unknown type for pnetcdf file! xtype = %d", xtype));
		    
		}
	    else
		*request = PIO_REQ_NULL;
	    
	    vdesc->nreqs++;
	    flush_output_buffer(file, false, 0);
	}	    
#endif /* _PNETCDF */
#ifdef _NETCDF
	if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
	    switch(xtype)
	    {
	    case NC_BYTE:
		ierr = nc_put_vars_schar(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_CHAR:
		ierr = nc_put_vars_schar(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_SHORT:
		ierr = nc_put_vars_short(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_INT:
		ierr = nc_put_vars_int(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_FLOAT:
		ierr = nc_put_vars_float(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_DOUBLE:
		ierr = nc_put_vars_double(ncid, varid, (size_t *)start, (size_t *)count,
					  (ptrdiff_t *)stride, buf);
		break;
#ifdef _NETCDF4		
	    case NC_UBYTE:
		ierr = nc_put_vars_uchar(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_USHORT:
		ierr = nc_put_vars_ushort(ncid, varid, (size_t *)start, (size_t *)count,
					  (ptrdiff_t *)stride, buf);
		break;
	    case NC_UINT:
		ierr = nc_put_vars_uint(ncid, varid, (size_t *)start, (size_t *)count,
					(ptrdiff_t *)stride, buf);
		break;
	    case NC_INT64:
		ierr = nc_put_vars_longlong(ncid, varid, (size_t *)start, (size_t *)count,
					    (ptrdiff_t *)stride, buf);
		break;
	    case NC_UINT64:
		ierr = nc_put_vars_ulonglong(ncid, varid, (size_t *)start, (size_t *)count,
					     (ptrdiff_t *)stride, buf);
		break;
	    /* case NC_STRING: */
	    /* 	ierr = nc_put_vars_string(ncid, varid, (size_t *)start, (size_t *)count, */
	    /* 				  (ptrdiff_t *)stride, (void *)buf); */
	    /* 	break; */
	    default:
		ierr = nc_put_vars(ncid, varid, (size_t *)start, (size_t *)count,
				   (ptrdiff_t *)stride, buf);
#endif /* _NETCDF4 */		
	    }
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	return check_mpi(file, mpierr, __FILE__, __LINE__);
    if (ierr)
	return check_netcdf(file, ierr, __FILE__, __LINE__);

    return ierr;
}

/**
 * PIO interface to nc_put_vars
 *
 * This routine is called collectively by all tasks in the
 * communicator ios.union_comm.

 * Refer to the <A
 * HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html">
 * netcdf documentation. </A> */
int PIOc_put_vars(int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[],
		  const PIO_Offset stride[], const void *buf, PIO_Offset bufcount,
		  MPI_Datatype buftype) 
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
  msg = PIO_MSG_PUT_VARS;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,   buf);;
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
	ierr = ncmpi_bput_vars(file->fh, varid, start, count, stride, buf, bufcount, buftype, request);;
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

int PIOc_put_vars_uchar(int ncid, int varid, const PIO_Offset start[],
			const PIO_Offset count[], const PIO_Offset stride[],
			const unsigned char *op) 
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_uchar(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_vars_ushort
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned short *op) 
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
      ierr = nc_put_vars_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
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
	ierr = ncmpi_bput_vars_ushort(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_vars_ulonglong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned long long *op) 
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
      ierr = nc_put_vars_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
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
	ierr = ncmpi_bput_vars_ulonglong(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_varm
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
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
      ierr = nc_put_varm(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
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
	ierr = ncmpi_bput_varm(file->fh, varid, start, count, stride, imap, buf, bufcount, buftype, request);;
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

///
/// PIO interface to nc_put_vars_uint
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned int *op) 
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
      ierr = nc_put_vars_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
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
	ierr = ncmpi_bput_vars_uint(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_varm_uchar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned char *op) 
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
      ierr = nc_put_varm_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_uchar(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_var_ushort
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_ushort (int ncid, int varid, const unsigned short *op) 
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_ushort(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var1_longlong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_longlong (int ncid, int varid, const PIO_Offset index[], const long long *op) 
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_longlong(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_vara_uchar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned char *op) 
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
      ierr = nc_put_vara_uchar(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_uchar(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_uchar(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_varm_short
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const short *op) 
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
      ierr = nc_put_varm_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_short(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_var1_long
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_long (int ncid, int varid, const PIO_Offset index[], const long *ip) 
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_long(file->fh, varid, index, ip, request);;
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

///
/// PIO interface to nc_put_vars_long
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const long *op) 
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_long(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_var_short
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_short (int ncid, int varid, const short *op) 
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_short(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_vara_int
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_int(int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[],
		      const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;
  int size;
  int ret;
  int ndims;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARA_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&varid, 1, MPI_INT,  ios->compmaster, ios->intercomm);
    if ((ret = PIOc_inq_varndims(ncid, varid, &ndims)))
	return ret;
    mpierr = MPI_Bcast((void *)start, ndims, MPI_INT, ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast((void *)count, ndims, MPI_INT, ios->compmaster, ios->intercomm);
    for (int d = 0, size = 1; d < ndims; d++)
	size *= count[d];
    mpierr = MPI_Bcast((void *)op, size, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_int(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_int(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_int(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_var1_ushort
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_ushort (int ncid, int varid, const PIO_Offset index[], const unsigned short *op) 
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
  msg = PIO_MSG_PUT_VAR1_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_ushort(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_vara_text
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const char *op) 
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
  msg = PIO_MSG_PUT_VARA_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_text(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_text(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_text(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_varm_text
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const char *op) 
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
  msg = PIO_MSG_PUT_VARM_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_text(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_varm_ushort
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned short *op) 
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
  msg = PIO_MSG_PUT_VARM_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_ushort(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_var_ulonglong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_ulonglong (int ncid, int varid, const unsigned long long *op) 
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
  msg = PIO_MSG_PUT_VAR_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_ulonglong(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var_int
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_int (int ncid, int varid, const int *op) 
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
  msg = PIO_MSG_PUT_VAR_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_int(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var_longlong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_longlong (int ncid, int varid, const long long *op) 
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
  msg = PIO_MSG_PUT_VAR_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_longlong(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var_schar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_schar (int ncid, int varid, const signed char *op) 
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
  msg = PIO_MSG_PUT_VAR_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_schar(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var_uint
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_uint (int ncid, int varid, const unsigned int *op) 
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
  msg = PIO_MSG_PUT_VAR_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_uint(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var (int ncid, int varid, const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
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
  msg = PIO_MSG_PUT_VAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var(file->fh, varid, buf, bufcount, buftype, request);;
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

///
/// PIO interface to nc_put_vara_ushort
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned short *op) 
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
  msg = PIO_MSG_PUT_VARA_USHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_ushort(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_ushort(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_ushort(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_vars_short
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const short *op) 
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
  msg = PIO_MSG_PUT_VARS_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_short(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_vara_uint
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned int *op) 
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
  msg = PIO_MSG_PUT_VARA_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_uint(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_uint(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_uint(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_vara_schar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const signed char *op) 
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
  msg = PIO_MSG_PUT_VARA_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_schar(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_schar(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_schar(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_varm_ulonglong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned long long *op) 
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
  msg = PIO_MSG_PUT_VARM_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_ulonglong(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_var1_uchar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_uchar (int ncid, int varid, const PIO_Offset index[], const unsigned char *op) 
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
  msg = PIO_MSG_PUT_VAR1_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_uchar(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_varm_int
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const int *op) 
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
  msg = PIO_MSG_PUT_VARM_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_int(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_vars_schar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const signed char *op) 
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
  msg = PIO_MSG_PUT_VARS_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_schar(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_var1
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1 (int ncid, int varid, const PIO_Offset index[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
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
  msg = PIO_MSG_PUT_VAR1;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1(file->fh, varid, index, buf, bufcount, buftype, request);;
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

///
/// PIO interface to nc_put_vara_float
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const float *op) 
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
  msg = PIO_MSG_PUT_VARA_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vara_float(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_var1_float
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_float (int ncid, int varid, const PIO_Offset index[], const float *op) 
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
  msg = PIO_MSG_PUT_VAR1_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_float(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_varm_float
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const float *op) 
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
  msg = PIO_MSG_PUT_VARM_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_float(file->fh, varid,(size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_float(file->fh, varid,(size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_float(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_var1_text
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_text (int ncid, int varid, const PIO_Offset index[], const char *op) 
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
  msg = PIO_MSG_PUT_VAR1_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_text(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_vars_text
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const char *op) 
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
  msg = PIO_MSG_PUT_VARS_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_text(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_varm_long
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long *op) 
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
  msg = PIO_MSG_PUT_VARM_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_long(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_vars_double
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const double *op) 
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
  msg = PIO_MSG_PUT_VARS_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_double(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_vara_longlong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const long long *op) 
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
  msg = PIO_MSG_PUT_VARA_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_longlong(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_longlong(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_longlong(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_var_double
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_double (int ncid, int varid, const double *op) 
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
  msg = PIO_MSG_PUT_VAR_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_double(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var_float
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_float (int ncid, int varid, const float *op) 
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
  msg = PIO_MSG_PUT_VAR_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_float(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var1_ulonglong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_ulonglong (int ncid, int varid, const PIO_Offset index[], const unsigned long long *op) 
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
  msg = PIO_MSG_PUT_VAR1_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_var1_ulonglong(file->fh, varid, (size_t *) index, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_var1_ulonglong(file->fh, varid, (size_t *) index, op);;
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
	ierr = ncmpi_bput_var1_ulonglong(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_varm_uint
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned int *op) 
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
  msg = PIO_MSG_PUT_VARM_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_uint(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_var1_uint
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_uint (int ncid, int varid, const PIO_Offset index[], const unsigned int *op) 
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
  msg = PIO_MSG_PUT_VAR1_UINT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_uint(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_var1_int
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_int (int ncid, int varid, const PIO_Offset index[], const int *op) 
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
  msg = PIO_MSG_PUT_VAR1_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_int(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_vars_float
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const float *op) 
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
  msg = PIO_MSG_PUT_VARS_FLOAT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_float(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_vara_short
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const short *op) 
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
  msg = PIO_MSG_PUT_VARA_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_short(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_short(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_short(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_var1_schar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_schar (int ncid, int varid, const PIO_Offset index[], const signed char *op) 
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
  msg = PIO_MSG_PUT_VAR1_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_schar(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_vara_ulonglong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned long long *op) 
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
  msg = PIO_MSG_PUT_VARA_ULONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_ulonglong(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_varm_double
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const double *op) 
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
  msg = PIO_MSG_PUT_VARM_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_double(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_vara
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
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
  msg = PIO_MSG_PUT_VARA;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara(file->fh, varid, (size_t *) start, (size_t *) count,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara(file->fh, varid, (size_t *) start, (size_t *) count,   buf);;
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
	ierr = ncmpi_bput_vara(file->fh, varid, start, count, buf, bufcount, buftype, request);;
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

///
/// PIO interface to nc_put_vara_long
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const long *op) 
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
  msg = PIO_MSG_PUT_VARA_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_long(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_long(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_long(file->fh, varid, start, count, op, request);;
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

///
/// PIO interface to nc_put_var1_double
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_double (int ncid, int varid, const PIO_Offset index[], const double *op) 
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
  msg = PIO_MSG_PUT_VAR1_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_double(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_varm_schar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const signed char *op) 
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
  msg = PIO_MSG_PUT_VARM_SCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_schar(file->fh, varid, start, count, stride, imap, op, request);;
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

///
/// PIO interface to nc_put_var_text
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_text (int ncid, int varid, const char *op) 
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
  msg = PIO_MSG_PUT_VAR_TEXT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_text(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_vars_int
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const int *op) 
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
  msg = PIO_MSG_PUT_VARS_INT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_vars_int(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_var1_short
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var1_short (int ncid, int varid, const PIO_Offset index[], const short *op) 
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
  msg = PIO_MSG_PUT_VAR1_SHORT;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var1_short(file->fh, varid, index, op, request);;
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

///
/// PIO interface to nc_put_vars_longlong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vars_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const long long *op) 
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
  msg = PIO_MSG_PUT_VARS_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vars_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vars_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, op);;
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
	ierr = ncmpi_bput_vars_longlong(file->fh, varid, start, count, stride, op, request);;
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

///
/// PIO interface to nc_put_vara_double
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_vara_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const double *op) 
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
  msg = PIO_MSG_PUT_VARA_DOUBLE;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_vara_double(file->fh, varid, (size_t *) start, (size_t *) count, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_vara_double(file->fh, varid, (size_t *) start, (size_t *) count, op);;
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
	ierr = ncmpi_bput_vara_double(file->fh, varid, start, count, op, request);;
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


///
/// PIO interface to nc_put_var_uchar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_uchar (int ncid, int varid, const unsigned char *op) 
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
  msg = PIO_MSG_PUT_VAR_UCHAR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_uchar(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_var_long
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_var_long (int ncid, int varid, const long *op) 
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
  msg = PIO_MSG_PUT_VAR_LONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
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
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_var_long(file->fh, varid, op, request);;
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

///
/// PIO interface to nc_put_varm_longlong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long long *op) 
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
  msg = PIO_MSG_PUT_VARM_LONGLONG;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
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
	ierr = ncmpi_bput_varm_longlong(file->fh, varid, start, count, stride, imap, op, request);;
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

