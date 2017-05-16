/**
 * @file
 * PIO functions to write data. 
 *
 * @author Ed Hartnett
 * @date  2016
 * @see http://code.google.com/p/parallelio/
 */

#include <config.h>
#include <pio.h>
#include <pio_internal.h>

/**
 * Internal PIO function which provides a type-neutral interface to
 * nc_put_vars.
 *
 * Users should not call this function directly. Instead, call one of
 * the derived functions, depending on the type of data you are
 * writing: PIOc_put_vars_text(), PIOc_put_vars_uchar(),
 * PIOc_put_vars_schar(), PIOc_put_vars_ushort(),
 * PIOc_put_vars_short(), PIOc_put_vars_uint(), PIOc_put_vars_int(),
 * PIOc_put_vars_long(), PIOc_put_vars_float(),
 * PIOc_put_vars_longlong(), PIOc_put_vars_double(),
 * PIOc_put_vars_ulonglong().
 *
 * This routine is called collectively by all tasks in the
 * communicator ios.union_comm.
 *
 * @param ncid identifies the netCDF file
 * @param varid the variable ID number
 * @param start an array of start indicies (must have same number of
 * entries as variable has dimensions). If NULL, indices of 0 will be
 * used.
 * @param count an array of counts (must have same number of entries
 * as variable has dimensions). If NULL, counts matching the size of
 * the variable will be used.
 * @param stride an array of strides (must have same number of
 * entries as variable has dimensions). If NULL, strides of 1 will be
 * used.
 * @param xtype the netCDF type of the data being passed in buf. Data
 * will be automatically covnerted from this type to the type of the
 * variable being written to.
 * @param buf pointer to the data to be written.
 *
 * @return PIO_NOERR on success, error code otherwise.
 * @private
 */
int PIOc_put_vars_tc(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		     const PIO_Offset *stride, nc_type xtype, const void *buf)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
    int ndims; /* The number of dimensions in the variable. */
    int *dimids; /* The IDs of the dimensions for this variable. */
    PIO_Offset typelen; /* Size (in bytes) of the data type of data in buf. */
    PIO_Offset num_elem = 1; /* Number of data elements in the buffer. */
    char start_present = start ? true : false; /* Is start non-NULL? */
    char count_present = count ? true : false; /* Is count non-NULL? */
    char stride_present = stride ? true : false; /* Is stride non-NULL? */
    PIO_Offset *rstart, *rcount, *rstride;
    var_desc_t *vdesc;
    PIO_Offset usage;
    int *request;

    LOG((1, "PIOc_put_vars_tc ncid = %d varid = %d start = %d count = %d "
    	 "stride = %d xtype = %d", ncid, varid, start, count, stride, xtype));

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
	/* Get the number of dims for this var. */
	if ((ierr = PIOc_inq_varndims(ncid, varid, &ndims)))
	    return check_netcdf(file, ierr, __FILE__, __LINE__);

	/* Get the length of the data type. */
	if ((ierr = PIOc_inq_type(ncid, xtype, NULL, &typelen)))
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

	/* Allocate memory for these arrays, now that we know ndims. */
	if (!(rstart = malloc(ndims * sizeof(PIO_Offset))))
	    return check_netcdf(file, PIO_ENOMEM, __FILE__, __LINE__);
	if (!(rcount = malloc(ndims * sizeof(PIO_Offset))))
	    return check_netcdf(file, PIO_ENOMEM, __FILE__, __LINE__);
	if (!(rstride = malloc(ndims * sizeof(PIO_Offset))))
	    return check_netcdf(file, PIO_ENOMEM, __FILE__, __LINE__);

	/* Figure out the real start, count, and stride arrays. (The
	 * user may have passed in NULLs.) */
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
	    int msg = PIO_MSG_PUT_VARS;

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
	    LOG((2, "PIOc_put_vars_tc ncid = %d varid = %d ndims = %d start_present = %d "
		 "count_present = %d stride_present = %d xtype = %d num_elem = %d", ncid, varid,
		 ndims, start_present, count_present, stride_present, xtype, num_elem));

	    /* Send the data. */
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)buf, num_elem * typelen, MPI_BYTE, ios->compmaster,
				   ios->intercomm);
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
	    return check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    check_mpi(file, mpierr, __FILE__, __LINE__);
	LOG((2, "PIOc_put_vars_tc checked mpierr = %d", mpierr));

	/* Broadcast values currently only known on computation tasks to IO tasks. */
	LOG((2, "PIOc_put_vars_tc bcast from comproot"));
	if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->comproot, ios->my_comm)))
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
	LOG((2, "PIOc_put_vars_tc complete bcast from comproot ndims = %d", ndims));
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
	if (file->iotype == PIO_IOTYPE_PNETCDF)
	{
	    PIO_Offset *fake_stride;
		
	    if (!stride_present)
	    {
		LOG((2, "stride not present"));
		if (!(fake_stride = malloc(ndims * sizeof(PIO_Offset))))
		    return PIO_ENOMEM;
		for (int d = 0; d < ndims; d++)
		    fake_stride[d] = 1;
	    }
	    else
		fake_stride = (PIO_Offset *)stride;
	    
	    LOG((2, "PIOc_put_vars_tc calling pnetcdf function"));
	    vdesc = file->varlist + varid;
	    if (vdesc->nreqs % PIO_REQUEST_ALLOC_CHUNK == 0)
		vdesc->request = realloc(vdesc->request,
					 sizeof(int) * (vdesc->nreqs + PIO_REQUEST_ALLOC_CHUNK));
	    request = vdesc->request + vdesc->nreqs;
	    LOG((2, "PIOc_put_vars_tc request = %d", vdesc->request));	    

	    /* Only the IO master actually does the call. */
	    if (ios->iomaster)
	    {
/*		LOG((2, "PIOc_put_vars_tc ncid = %d varid = %d start[0] = %d count[0] = %d fake_stride[0] = %d",
		ncid, varid, start[0], count[0], fake_stride[0]));*/
		/* for (int d = 0; d < ndims; d++) */
		/*     LOG((2, "start[%d] = %d count[%d] = %d stride[%d] = %d", d, start[d], d, count[d], d, stride[d])); */
		switch(xtype)
		{
		case NC_BYTE:
		    ierr = ncmpi_bput_vars_schar(ncid, varid, start, count, fake_stride, buf, request);
		    break;
		case NC_CHAR:
		    ierr = ncmpi_bput_vars_text(ncid, varid, start, count, fake_stride, buf, request);
		    break;
		case NC_SHORT:
		    ierr = ncmpi_bput_vars_short(ncid, varid, start, count, fake_stride, buf, request);
		    break;
		case NC_INT:
		    LOG((2, "PIOc_put_vars_tc io_rank 0 doing pnetcdf for int"));
		    ierr = ncmpi_bput_vars_int(ncid, varid, start, count, fake_stride, buf, request);
		    LOG((2, "PIOc_put_vars_tc io_rank 0 done with pnetcdf call for int ierr = %d", ierr));	    		    
		    break;
		case NC_FLOAT:
		    ierr = ncmpi_bput_vars_float(ncid, varid, start, count, fake_stride, buf, request);
		    break;
		case NC_DOUBLE:
		    ierr = ncmpi_bput_vars_double(ncid, varid, start, count, fake_stride, buf, request);
		    break;
		case NC_INT64:
		    ierr = ncmpi_bput_vars_longlong(ncid, varid, start, count, fake_stride, buf, request);
		    break;
		default:
		    LOG((0, "Unknown type for pnetcdf file! xtype = %d", xtype));
		}
		LOG((2, "PIOc_put_vars_tc io_rank 0 done with pnetcdf call"));	    
	    }
	    else
		*request = PIO_REQ_NULL;

	    vdesc->nreqs++;
	    LOG((2, "PIOc_put_vars_tc flushing output buffer"));
	    flush_output_buffer(file, false, 0);
	    LOG((2, "PIOc_put_vars_tc flushed output buffer"));

	    /* Free malloced resources. */
	    if (!stride_present)
		free(fake_stride);
	}
#endif /* _PNETCDF */
#ifdef _NETCDF
	if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
	{
	    LOG((2, "PIOc_put_vars_tc calling netcdf function file->iotype = %d",
		 file->iotype));	    
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
	}
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    LOG((2, "PIOc_put_vars_tc bcasting netcdf return code %d", ierr));
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	return check_mpi(file, mpierr, __FILE__, __LINE__);
    if (ierr)
	return check_netcdf(file, ierr, __FILE__, __LINE__);
    LOG((2, "PIOc_put_vars_tc bcast netcdf return code %d complete", ierr));

    return ierr;
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_text(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		       const PIO_Offset *stride, const char *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_CHAR, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_uchar(int ncid, int varid, const PIO_Offset *start,
			const PIO_Offset *count, const PIO_Offset *stride,
			const unsigned char *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_UBYTE, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_schar(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
			const PIO_Offset *stride, const signed char *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_BYTE, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_ushort(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
			 const PIO_Offset *stride, const unsigned short *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_USHORT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_short(int ncid, int varid, const PIO_Offset *start,
			const PIO_Offset *count, const PIO_Offset *stride, const short *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_SHORT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_uint(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		       const PIO_Offset *stride, const unsigned int *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_UINT, op);
}

/** PIO interface to nc_put_vars_int */
int PIOc_put_vars_int(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		      const PIO_Offset *stride, const int *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_INT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_long(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		       const PIO_Offset *stride, const long *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_INT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_float(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
			const PIO_Offset *stride, const float *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_FLOAT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_longlong(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
			   const PIO_Offset *stride, const long long *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_INT64, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_double(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
			 const PIO_Offset *stride, const double *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_DOUBLE, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vars_ulonglong(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
			    const PIO_Offset *stride, const unsigned long long *op)
{
    return PIOc_put_vars_tc(ncid, varid, start, count, stride, NC_UINT64, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_tc(int ncid, int varid, const PIO_Offset *index, nc_type xtype,
		     const void *op)
{
    int ndims;
    int ierr;

    /* Find the number of dimensions. */
    if ((ierr = PIOc_inq_varndims(ncid, varid, &ndims)))
	return ierr;

    /* Set up count array. */
    PIO_Offset count[ndims];
    for (int c = 0; c < ndims; c++)
	count[c] = 1;

    return PIOc_put_vars_tc(ncid, varid, index, count, NULL, xtype, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_text(int ncid, int varid, const PIO_Offset *index, const char *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_CHAR, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_uchar(int ncid, int varid, const PIO_Offset *index,
			const unsigned char *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_UBYTE, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_schar(int ncid, int varid, const PIO_Offset *index,
			const signed char *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_BYTE, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_ushort(int ncid, int varid, const PIO_Offset *index,
			 const unsigned short *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_USHORT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_short(int ncid, int varid, const PIO_Offset *index,
			const short *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_SHORT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_uint(int ncid, int varid, const PIO_Offset *index,
		       const unsigned int *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_UINT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_int(int ncid, int varid, const PIO_Offset *index, const int *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_INT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_float(int ncid, int varid, const PIO_Offset *index, const float *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_FLOAT, op);    
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_long(int ncid, int varid, const PIO_Offset *index, const long *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_LONG, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_double(int ncid, int varid, const PIO_Offset *index,
			 const double *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_DOUBLE, op);    
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_ulonglong(int ncid, int varid, const PIO_Offset *index,
			    const unsigned long long *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_UINT64, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var1_longlong(int ncid, int varid, const PIO_Offset *index,
			   const long long *op)
{
    return PIOc_put_var1_tc(ncid, varid, index, NC_INT64, op);    
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_text(int ncid, int varid, const PIO_Offset *start,
			const PIO_Offset *count, const char *op)
{
    return PIOc_put_vars_text(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_uchar(int ncid, int varid, const PIO_Offset *start,
			const PIO_Offset *count, const unsigned char *op)
{
    return PIOc_put_vars_uchar(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_schar(int ncid, int varid, const PIO_Offset *start,
			const PIO_Offset *count, const signed char *op)
{
    return PIOc_put_vars_schar(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_ushort(int ncid, int varid, const PIO_Offset *start,
			 const PIO_Offset *count, const unsigned short *op)
{
    return PIOc_put_vars_ushort(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_short(int ncid, int varid, const PIO_Offset *start,
			const PIO_Offset *count, const short *op)
{
    return PIOc_put_vars_short(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_uint(int ncid, int varid, const PIO_Offset *start,
		       const PIO_Offset *count, const unsigned int *op)
{
    return PIOc_put_vars_uint(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_int(int ncid, int varid, const PIO_Offset *start,
		      const PIO_Offset *count, const int *op)
{
    return PIOc_put_vars_int(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_long(int ncid, int varid, const PIO_Offset *start,
		       const PIO_Offset *count, const long *op)
{
    return PIOc_put_vars_long(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_float(int ncid, int varid, const PIO_Offset *start,
			const PIO_Offset *count, const float *op)
{
    return PIOc_put_vars_float(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_ulonglong(int ncid, int varid, const PIO_Offset *start,
			    const PIO_Offset *count, const unsigned long long *op)
{
    return PIOc_put_vars_ulonglong(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_longlong(int ncid, int varid, const PIO_Offset *start,
			   const PIO_Offset *count, const long long *op)
{
    return PIOc_put_vars_longlong(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_vara_double(int ncid, int varid, const PIO_Offset *start,
			 const PIO_Offset *count, const double *op)
{
    return PIOc_put_vars_double(ncid, varid, start, count, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_text(int ncid, int varid, const char *op)
{
    return PIOc_put_vars_text(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_uchar(int ncid, int varid, const unsigned char *op)
{
    return PIOc_put_vars_uchar(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_schar(int ncid, int varid, const signed char *op)
{
    return PIOc_put_vars_schar(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_ushort(int ncid, int varid, const unsigned short *op)
{
    return PIOc_put_vars_tc(ncid, varid, NULL, NULL, NULL, NC_USHORT, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_short(int ncid, int varid, const short *op)
{
    return PIOc_put_vars_short(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_uint(int ncid, int varid, const unsigned int *op)
{
    return PIOc_put_vars_uint(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_int(int ncid, int varid, const int *op)
{
    return PIOc_put_vars_int(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_long(int ncid, int varid, const long *op)
{
    return PIOc_put_vars_long(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_float(int ncid, int varid, const float *op)
{
    return PIOc_put_vars_float(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_ulonglong(int ncid, int varid, const unsigned long long *op)
{
    return PIOc_put_vars_ulonglong(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_longlong(int ncid, int varid, const long long *op)
{
    return PIOc_put_vars_longlong(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var_double(int ncid, int varid, const double *op)
{
    return PIOc_put_vars_double(ncid, varid, NULL, NULL, NULL, op);
}

/** Interface to netCDF data write function. */
int PIOc_put_var(int ncid, int varid, const void *buf, PIO_Offset bufcount,
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

/**
 * PIO interface to nc_put_vars
 *
 * This routine is called collectively by all tasks in the
 * communicator ios.union_comm.

 * Refer to the <A
 * HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html">
 * netcdf documentation. </A> */
int PIOc_put_vars(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
		  const PIO_Offset *stride, const void *buf, PIO_Offset bufcount,
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
		ierr = nc_put_vars(file->fh, varid, (size_t *)start, (size_t *)count,
				   (ptrdiff_t *)stride,   buf);;
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
		ierr = ncmpi_bput_vars(file->fh, varid, start, count, stride, buf,
				       bufcount, buftype, request);;
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

/** Interface to netCDF data write function. */
int PIOc_put_var1(int ncid, int varid, const PIO_Offset *index, const void *buf,
		  PIO_Offset bufcount, MPI_Datatype buftype)
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

/** Interface to netCDF data write function. */
int PIOc_put_vara(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count, const void *buf,
		  PIO_Offset bufcount, MPI_Datatype buftype)
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
