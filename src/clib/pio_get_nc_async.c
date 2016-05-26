#include <config.h>
#include <pio.h>
#include <pio_internal.h>

int PIOc_get_vars_tc(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		     const PIO_Offset *stride, nc_type xtype, void *buf) 
{
    iosystem_desc_t *ios;  /** Pointer to io system information. */
    file_desc_t *file;     /** Pointer to file information. */
    int ierr = PIO_NOERR;  /** Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    int ndims; /** The number of dimensions in the variable. */
    int *dimids; /** The IDs of the dimensions for this variable. */
    PIO_Offset typelen; /** Size (in bytes) of the data type of data in buf. */
    size_t num_elem = 1; /** Number of data elements in the buffer. */
    int bcast = false;

    LOG((1, "PIOc_get_vars_tc ncid = %d varid = %d start = %d count = %d "
    	 "stride = %d xtype = %d", ncid, varid, start, count, stride, xtype));

    /* User must provide a place to put some data. */
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

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_GET_VARS;
	    char start_present = start ? true : false;
	    char count_present = count ? true : false;
	    char stride_present = stride ? true : false;

	    if(ios->compmaster)
		mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

	    /* Send the function parameters and associated informaiton
	     * to the msg handler. */
	    if (!mpierr)
		mpierr = MPI_Bcast(&ncid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&start_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr && start_present)
		mpierr = MPI_Bcast((PIO_Offset *)start, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&count_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr && count_present)
		mpierr = MPI_Bcast((PIO_Offset *)count, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&stride_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr && stride_present)
		mpierr = MPI_Bcast((PIO_Offset *)stride, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&xtype, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&num_elem, 1, MPI_OFFSET, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, ios->compmaster, ios->intercomm);
	    LOG((2, "PIOc_get_vars_tc ncid = %d varid = %d ndims = %d start_present = %d "
		 "count_present = %d stride_present = %d xtype = %d num_elem = %d", ncid, varid,
		 ndims, start_present, count_present, stride_present, xtype, num_elem));
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    return check_mpi(file, mpierr2, __FILE__, __LINE__);
	check_mpi(file, mpierr, __FILE__, __LINE__);

	/* Broadcast values currently only known on computation tasks to IO tasks. */
	if ((mpierr = MPI_Bcast(&num_elem, 1, MPI_OFFSET, ios->comproot, ios->my_comm)))
	    check_mpi(file, mpierr, __FILE__, __LINE__);
	if ((mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, ios->comproot, ios->my_comm)))
	    check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
	if (file->iotype == PIO_IOTYPE_PNETCDF)
	{
#ifdef PNET_READ_AND_BCAST
	    LOG((1, "PNET_READ_AND_BCAST"));
	    ncmpi_begin_indep_data(file->fh);
	    if (ios->iomaster)
	    {
		switch(xtype)
		{
		case NC_BYTE:
		    ierr = ncmpi_get_vars_schar(ncid, varid, start, count, stride, buf);
		    break;
		case NC_CHAR:
		    ierr = ncmpi_get_vars_text(ncid, varid, start, count, stride, buf);
		    break;
		case NC_SHORT:
		    ierr = ncmpi_get_vars_short(ncid, varid, start, count, stride, buf);
		    break;
		case NC_INT:
		    ierr = ncmpi_get_vars_int(ncid, varid, start, count, stride, buf);
		    break;
		case NC_FLOAT:
		    ierr = ncmpi_get_vars_float(ncid, varid, start, count, stride, buf);
		    break;
		case NC_DOUBLE:
		    ierr = ncmpi_get_vars_double(ncid, varid, start, count, stride, buf);
		    break;
		case NC_INT64:
		    ierr = ncmpi_get_vars_longlong(ncid, varid, start, count, stride, buf);
		    break;
		default:
		    LOG((0, "Unknown type for pnetcdf file! xtype = %d", xtype));
		}
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else /* PNET_READ_AND_BCAST */
	    LOG((1, "not PNET_READ_AND_BCAST"));
	    switch(xtype)
	    {
	    case NC_BYTE:
		ierr = ncmpi_get_vars_schar_all(ncid, varid, start, count, stride, buf);
		break;
	    case NC_CHAR:
		ierr = ncmpi_get_vars_text_all(ncid, varid, start, count, stride, buf);
		break;
	    case NC_SHORT:
		ierr = ncmpi_get_vars_short_all(ncid, varid, start, count, stride, buf);
		break;
	    case NC_INT:
		ierr = ncmpi_get_vars_int_all(ncid, varid, start, count, stride, buf);
		for (int i = 0; i < 4; i++)
		    LOG((2, "((int *)buf)[%d] = %d", i, ((int *)buf)[0]));
		break;
	    case NC_FLOAT:
		ierr = ncmpi_get_vars_float_all(ncid, varid, start, count, stride, buf);
		break;
	    case NC_DOUBLE:
		ierr = ncmpi_get_vars_double_all(ncid, varid, start, count, stride, buf);
		break;
	    case NC_INT64:
		ierr = ncmpi_get_vars_longlong_all(ncid, varid, start, count, stride, buf);
		break;
	    default:
		LOG((0, "Unknown type for pnetcdf file! xtype = %d", xtype));
	    }
#endif /* PNET_READ_AND_BCAST */
	}
#endif /* _PNETCDF */
#ifdef _NETCDF
	if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
	    switch(xtype)
	    {
	    case NC_BYTE:
		ierr = nc_get_vars_schar(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_CHAR:
		ierr = nc_get_vars_schar(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_SHORT:
		ierr = nc_get_vars_short(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_INT:
		ierr = nc_get_vars_int(ncid, varid, (size_t *)start, (size_t *)count,
				       (ptrdiff_t *)stride, buf);
		break;
	    case NC_FLOAT:
		ierr = nc_get_vars_float(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_DOUBLE:
		ierr = nc_get_vars_double(ncid, varid, (size_t *)start, (size_t *)count,
					  (ptrdiff_t *)stride, buf);
		break;
#ifdef _NETCDF4
	    case NC_UBYTE:
		ierr = nc_get_vars_uchar(ncid, varid, (size_t *)start, (size_t *)count,
					 (ptrdiff_t *)stride, buf);
		break;
	    case NC_USHORT:
		ierr = nc_get_vars_ushort(ncid, varid, (size_t *)start, (size_t *)count,
					  (ptrdiff_t *)stride, buf);
		break;
	    case NC_UINT:
		ierr = nc_get_vars_uint(ncid, varid, (size_t *)start, (size_t *)count,
					(ptrdiff_t *)stride, buf);
		break;
	    case NC_INT64:
		ierr = nc_get_vars_longlong(ncid, varid, (size_t *)start, (size_t *)count,
					    (ptrdiff_t *)stride, buf);
		break;
	    case NC_UINT64:
		ierr = nc_get_vars_ulonglong(ncid, varid, (size_t *)start, (size_t *)count,
					     (ptrdiff_t *)stride, buf);
		break;
		/* case NC_STRING: */
		/* 	ierr = nc_get_vars_string(ncid, varid, (size_t *)start, (size_t *)count, */
		/* 				  (ptrdiff_t *)stride, (void *)buf); */
		/* 	break; */
	    default:
		ierr = nc_get_vars(ncid, varid, (size_t *)start, (size_t *)count,
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

    /* Send the data. */
    LOG((2, "PIOc_get_vars_tc bcasting data num_elem = %d typelen = %d", num_elem,
	 typelen));    
    if (!mpierr)
	mpierr = MPI_Bcast((void *)buf, num_elem * typelen, MPI_BYTE, ios->ioroot,
			   ios->my_comm);
    return ierr;
}

int PIOc_get_vars_int(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		      const PIO_Offset *stride, int *buf) 
{
    return PIOc_get_vars_tc(ncid, varid, start, count, stride, NC_INT, buf);
}

int PIOc_get_var1_schar (int ncid, int varid, const PIO_Offset *index, signed char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_SCHAR;
    ibuftype = MPI_CHAR;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_schar(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_schar(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_schar(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_schar_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_ulonglong (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, unsigned long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_ULONGLONG;
    ibuftype = MPI_UNSIGNED_LONG_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_ulonglong(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_ulonglong_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_short (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_SHORT;
    ibuftype = MPI_SHORT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_short(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_short_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_double (int ncid, int varid, double *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_DOUBLE;
    ibuftype = MPI_DOUBLE;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_double(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_double(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_double(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_double_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_double (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, double *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_DOUBLE;
    ibuftype = MPI_DOUBLE;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_double(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_double(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_double(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_double_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_int (int ncid, int varid, int *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_INT;
    ibuftype = MPI_INT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_int(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_int(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_int(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_int_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_ushort (int ncid, int varid, unsigned short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_USHORT;
    ibuftype = MPI_UNSIGNED_SHORT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_ushort(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_ushort(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_ushort(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_ushort_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_text (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_TEXT;
    ibuftype = MPI_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_text(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_text(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_text(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_text_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_int (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, int *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_INT;
    ibuftype = MPI_INT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_int(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_int(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_int(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_int_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_float (int ncid, int varid, const PIO_Offset *index, float *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_FLOAT;
    ibuftype = MPI_FLOAT;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_float(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_float(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_float(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_float_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_short (int ncid, int varid, const PIO_Offset *index, short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_SHORT;
    ibuftype = MPI_SHORT;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_short(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_short(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_short(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_short_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_text (int ncid, int varid, char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_TEXT;
    ibuftype = MPI_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_text(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_text(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_text(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_text_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_schar (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, signed char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_SCHAR;
    ibuftype = MPI_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_schar(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_schar_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_ushort (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, unsigned short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_USHORT;
    ibuftype = MPI_UNSIGNED_SHORT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_ushort(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_ushort(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_ushort(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_ushort_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_ushort (int ncid, int varid, const PIO_Offset *index, unsigned short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_USHORT;
    ibuftype = MPI_UNSIGNED_SHORT;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_ushort(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_ushort(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_ushort(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_ushort_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_float (int ncid, int varid, float *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_FLOAT;
    ibuftype = MPI_FLOAT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_float(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_float(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_float(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_float_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_uchar (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, unsigned char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_UCHAR;
    ibuftype = MPI_UNSIGNED_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_uchar(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_uchar_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var (int ncid, int varid, void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR;
    ibufcnt = bufcount;
    ibuftype = buftype;
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
	    ierr = nc_get_var(file->fh, varid,   buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var(file->fh, varid,   buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var(file->fh, varid, buf, bufcount, buftype);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_all(file->fh, varid, buf, bufcount, buftype);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_longlong (int ncid, int varid, const PIO_Offset *index, long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_LONGLONG;
    ibuftype = MPI_LONG_LONG;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_longlong(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_longlong(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_longlong(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_longlong_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_ushort (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, unsigned short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_USHORT;
    ibuftype = MPI_UNSIGNED_SHORT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_ushort(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_ushort_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_long (int ncid, int varid, long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_LONG;
    ibuftype = MPI_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_long(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_long(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_long(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_long_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_double (int ncid, int varid, const PIO_Offset *index, double *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_DOUBLE;
    ibuftype = MPI_DOUBLE;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_double(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_double(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_double(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_double_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_uint (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, unsigned int *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_UINT;
    ibuftype = MPI_UNSIGNED;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_uint(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_uint(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_uint(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_uint_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_longlong (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_LONGLONG;
    ibuftype = MPI_LONG_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_longlong(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_longlong_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_longlong (int ncid, int varid, long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_LONGLONG;
    ibuftype = MPI_LONG_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_longlong(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_longlong(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_longlong(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_longlong_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_short (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_SHORT;
    ibuftype = MPI_SHORT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_short(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_short(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_short(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_short_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_long (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_LONG;
    ibuftype = MPI_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_long(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_long(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_long(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_long_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_int (int ncid, int varid, const PIO_Offset *index, int *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_INT;
    ibuftype = MPI_INT;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_int(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_int(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_int(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_int_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_ulonglong (int ncid, int varid, const PIO_Offset *index, unsigned long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_ULONGLONG;
    ibuftype = MPI_UNSIGNED_LONG_LONG;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_ulonglong(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_ulonglong(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_ulonglong(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_ulonglong_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_uchar (int ncid, int varid, unsigned char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_UCHAR;
    ibuftype = MPI_UNSIGNED_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_uchar(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_uchar(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_uchar(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_uchar_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_uchar (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, unsigned char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_UCHAR;
    ibuftype = MPI_UNSIGNED_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_uchar(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_uchar(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_uchar(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_uchar_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_float (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, float *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_FLOAT;
    ibuftype = MPI_FLOAT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_float(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_float_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_long (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_LONG;
    ibuftype = MPI_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_long(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_long_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1 (int ncid, int varid, const PIO_Offset *index, void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1;
    ibufcnt = bufcount;
    ibuftype = buftype;
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
	    ierr = nc_get_var1(file->fh, varid, (size_t *) index,   buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1(file->fh, varid, (size_t *) index,   buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1(file->fh, varid, index, buf, bufcount, buftype);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_all(file->fh, varid, index, buf, bufcount, buftype);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_uint (int ncid, int varid, unsigned int *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_UINT;
    ibuftype = MPI_UNSIGNED;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_uint(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_uint(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_uint(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_uint_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA;
    ibufcnt = bufcount;
    ibuftype = buftype;
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
	    ierr = nc_get_vara(file->fh, varid, (size_t *) start, (size_t *) count,   buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara(file->fh, varid, (size_t *) start, (size_t *) count,   buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara(file->fh, varid, start, count, buf, bufcount, buftype);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_all(file->fh, varid, start, count, buf, bufcount, buftype);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_schar (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, signed char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_SCHAR;
    ibuftype = MPI_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_schar(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_schar(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_schar(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_schar_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_uint (int ncid, int varid, const PIO_Offset *index, unsigned int *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_UINT;
    ibuftype = MPI_UNSIGNED;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_uint(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_uint(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_uint(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_uint_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_uint (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, unsigned int *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_UINT;
    ibuftype = MPI_UNSIGNED;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_uint(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_uint_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_float (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, float *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_FLOAT;
    ibuftype = MPI_FLOAT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_float(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_float(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_float(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_float_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_text (int ncid, int varid, const PIO_Offset *index, char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_TEXT;
    ibuftype = MPI_CHAR;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_text(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_text(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_text(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_text_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_double (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, double *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_DOUBLE;
    ibuftype = MPI_DOUBLE;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_double(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_double_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_longlong (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_LONGLONG;
    ibuftype = MPI_LONG_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_longlong(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_longlong(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_longlong(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_longlong_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_ulonglong (int ncid, int varid, unsigned long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_ULONGLONG;
    ibuftype = MPI_UNSIGNED_LONG_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_ulonglong(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_ulonglong(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_ulonglong(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_ulonglong_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vara_ulonglong (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, unsigned long long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARA_ULONGLONG;
    ibuftype = MPI_UNSIGNED_LONG_LONG;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i];
    }
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
	    ierr = nc_get_vara_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vara_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vara_ulonglong(file->fh, varid, start, count,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vara_ulonglong_all(file->fh, varid, start, count,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_short (int ncid, int varid, short *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_SHORT;
    ibuftype = MPI_SHORT;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_short(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_short(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_short(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_short_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_long (int ncid, int varid, const PIO_Offset *index, long *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_LONG;
    ibuftype = MPI_LONG;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_long(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_long(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_long(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_long_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars_text (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS_TEXT;
    ibuftype = MPI_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    ibufcnt = 1;
    for(int i=0;i<ndims;i++){
	ibufcnt *= count[i]/stride[i];
    }
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
	    ierr = nc_get_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars_text(file->fh, varid, start, count, stride,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_text_all(file->fh, varid, start, count, stride,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var1_uchar (int ncid, int varid, const PIO_Offset *index, unsigned char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR1_UCHAR;
    ibuftype = MPI_UNSIGNED_CHAR;
    ibufcnt = 1;
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
	    ierr = nc_get_var1_uchar(file->fh, varid, (size_t *) index,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var1_uchar(file->fh, varid, (size_t *) index,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var1_uchar(file->fh, varid, index,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var1_uchar_all(file->fh, varid, index,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_vars (int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const PIO_Offset *stride, void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VARS;
    ibufcnt = bufcount;
    ibuftype = buftype;
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
	    ierr = nc_get_vars(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,   buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_vars(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,   buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_vars(file->fh, varid, start, count, stride, buf, bufcount, buftype);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_vars_all(file->fh, varid, start, count, stride, buf, bufcount, buftype);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

int PIOc_get_var_schar (int ncid, int varid, signed char *buf) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    MPI_Datatype ibuftype;
    int ndims;
    int ibufcnt;
    bool bcast = false;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_GET_VAR_SCHAR;
    ibuftype = MPI_CHAR;
    ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
    int dimid[ndims];
    PIO_Offset dimsize;
    ibufcnt = 1;
    PIOc_inq_vardimid(file->fh, varid, dimid);
    for(int i=0;i<ndims;i++){
	PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
	ibufcnt *= dimsize;
    }
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
	    ierr = nc_get_var_schar(file->fh, varid,  buf);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    bcast = true;
	    if(ios->iomaster){
		ierr = nc_get_var_schar(file->fh, varid,  buf);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
	    ncmpi_begin_indep_data(file->fh);
	    if(ios->iomaster){
		ierr = ncmpi_get_var_schar(file->fh, varid,  buf);;
	    };
	    ncmpi_end_indep_data(file->fh);
	    bcast=true;
#else
	    ierr = ncmpi_get_var_schar_all(file->fh, varid,  buf);;
#endif
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

    if(ios->async_interface || bcast ||  
       (ios->num_iotasks < ios->num_comptasks)){
	MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

