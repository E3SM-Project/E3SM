/** @file
 *
 * Functions to wrap netCDF-4 functions for PIO.
 **/
#include <pio.h>
#include <pio_internal.h>

/**
 * @ingroup PIO_def_var
 * Set deflate (zlib) settings for a variable.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable.
 * @param shuffle non-zero to turn on shuffle filter (can be good for
 * integer data).
 * @param deflate non-zero to turn on zlib compression for this
 * variable.
 * @param deflate_level 1 to 9, with 1 being faster and 9 being more
 * compressed.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_def_var_deflate(int ncid, int varid, int shuffle, int deflate,
			 int deflate_level)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_DEF_VAR_DEFLATE;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    /* Versions of netCDF 4.4 and earlier do not return an
	     * error when attempting to turn on deflation with
	     * parallel I.O. But this is not allowed by HDF5. So
	     * return the correct error code. */
	    ierr = NC_EINVAL;
	    //ierr = nc_def_var_deflate(file->fh, varid, shuffle, deflate, deflate_level);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (!ios->io_rank)
		ierr = nc_def_var_deflate(file->fh, varid, shuffle, deflate, deflate_level);
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* If there is an error, allocate space for the error string. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check the netCDF return code, and broadcast it to all tasks. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    return ierr;
}    

/**
 * @ingroup PIO_inq_var
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * Inquire about deflate (zlib compression) settings for a variable.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param shufflep pointer to an int that will get the status of the
 * shuffle filter.
 * @param deflatep pointer to an int that will be set to non-zero if
 * deflation is in use for this variable.
 * @param deflate_levelp pointer to an int that will get the deflation
 * level (from 1-9) if deflation is in use for this variable.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_inq_var_deflate(int ncid, int varid, int *shufflep,
			 int *deflatep, int *deflate_levelp)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;
    int ret;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_INQ_VAR_DEFLATE;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_var_deflate(file->fh, varid, shufflep, deflatep, deflate_levelp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (ios->io_rank == 0)
		ierr = nc_inq_var_deflate(file->fh, varid, shufflep, deflatep, deflate_levelp);
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* If there is an error, allocate space for the error string. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check the netCDF return code, and broadcast it to all tasks. */
    ierr = check_netcdf(file, ierr, errstr, __LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    /* Broadcast results to all tasks. */
    if (ierr == PIO_NOERR)
    {
	if (shufflep)
	    ierr = MPI_Bcast(shufflep, 1, MPI_INT, ios->ioroot, ios->my_comm);
	if (deflatep)
	    ierr = MPI_Bcast(deflatep, 1, MPI_INT, ios->ioroot, ios->my_comm);
	if (deflate_levelp)
	    ierr = MPI_Bcast(deflate_levelp, 1, MPI_INT, ios->ioroot, ios->my_comm);
    }
    return ierr;
}    

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
 * 
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * Chunksizes have important performance repercussions. NetCDF
 * attempts to choose sensible chunk sizes by default, but for best
 * performance check chunking against access patterns.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 *
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param storage NC_CONTIGUOUS or NC_CHUNKED.
 * @param chunksizep an array of chunksizes. Must have a chunksize for
 * every variable dimension.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_def_var_chunking(int ncid, int varid, int storage,
			  const PIO_Offset *chunksizesp) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_DEF_VAR_CHUNKING;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_def_var_chunking(file->fh, varid, storage, chunksizesp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (ios->io_rank==0)
	    {
		ierr = nc_def_var_chunking(file->fh, varid, storage, chunksizesp);
	    }
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    return PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    return PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Allocate an error string if needed. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check for netCDF error. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    return ierr;
}

/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param storagep pointer to int which will be set to either
 * NC_CONTIGUOUS or NC_CHUNKED.
 * @param chunksizep pointer to memory where chunksizes will be
 * set. There are the same number of chunksizes as there are
 * dimensions.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_inq_var_chunking(int ncid, int varid, int *storagep, PIO_Offset *chunksizesp)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;
    int ndims;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_INQ_VAR_CHUNKING;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_var_chunking(file->fh, varid, storagep, chunksizesp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (ios->io_rank == 0)
	    {
		if ((ierr = nc_inq_var_chunking(file->fh, varid, storagep, chunksizesp)))
		    return ierr;
		if ((ierr = nc_inq_varndims(file->fh, varid, &ndims)))
		    return ierr;
	    }
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    return PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    return PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* If there is an error, allocate space for the error string. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check the netCDF return code, and broadcast it to all tasks. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error stringif it was allocated. */
    if (errstr != NULL)
	free(errstr);

    /* Broadcast results to all tasks. */
    ierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->ioroot, ios->my_comm);
    if (storagep)
	ierr = MPI_Bcast(storagep, 1, MPI_INT, ios->ioroot, ios->my_comm);
    if (chunksizesp)
	ierr = MPI_Bcast(chunksizesp, ndims, MPI_OFFSET, ios->ioroot, ios->my_comm);

    return ierr;
}

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * Chunksizes have important performance repercussions. NetCDF
 * attempts to choose sensible chunk sizes by default, but for best
 * performance check chunking against access patterns.
 *
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param storage NC_CONTIGUOUS or NC_CHUNKED.
 * @param chunksizep an array of chunksizes. Must have a chunksize for
 * every variable dimension.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_def_var_fill(int ncid, int varid, int no_fill, const void *fill_value)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_SET_FILL;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_def_var_fill(file->fh, varid, no_fill, fill_value);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (ios->io_rank==0)
	    {
		ierr = nc_def_var_fill(file->fh, varid, no_fill, fill_value);
	    }
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    return PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    return PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Allocate an error string if needed. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check for netCDF error. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    return ierr;
}

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * Chunksizes have important performance repercussions. NetCDF
 * attempts to choose sensible chunk sizes by default, but for best
 * performance check chunking against access patterns.
 *
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param storage NC_CONTIGUOUS or NC_CHUNKED.
 * @param chunksizep an array of chunksizes. Must have a chunksize for
 * every variable dimension.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_def_var_endian(int ncid, int varid, int endian)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_DEF_VAR_ENDIAN;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_def_var_endian(file->fh, varid, endian);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (ios->io_rank==0)
	    {
		ierr = nc_def_var_endian(file->fh, varid, endian);
	    }
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Allocate an error string if needed. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check for netCDF error. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    return ierr;
}

/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param storagep pointer to int which will be set to either
 * NC_CONTIGUOUS or NC_CHUNKED.
 * @param chunksizep pointer to memory where chunksizes will be
 * set. There are the same number of chunksizes as there are
 * dimensions.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_inq_var_endian(int ncid, int varid, int *endianp)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_INQ_VAR_CHUNKING;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_var_endian(file->fh, varid, endianp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (!ios->io_rank)
		ierr = nc_inq_var_endian(file->fh, varid, endianp);
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* If there is an error, allocate space for the error string. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check the netCDF return code, and broadcast it to all tasks. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    /* Broadcast results to all tasks. */
    if (endianp)
	ierr = MPI_Bcast(endianp, 1, MPI_INT, ios->ioroot, ios->my_comm);

    return ierr;
}    

/**
 * @ingroup PIO_def_var
 *
 * Set chunk cache netCDF files to be opened/created.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * The file chunk cache for HDF5 can be set, and will apply for any
 * files opened or created until the program ends, or the settings are
 * changed again. The cache settings apply only to the open file. They
 * do not persist with the file, and must be set each time the file is
 * opened, before it is opened, if they are to have effect.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * @param iotype the iotype of files to be created or opened.
 * @param size size of file cache.
 * @param nelems number of elements in file cache.
 * @param preemption preemption setting for file cache.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_set_chunk_cache(int iosysid, int iotype, PIO_Offset size,
			 PIO_Offset nelems, float preemption)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;
    msg = PIO_MSG_SET_CHUNK_CACHE;

    ios = pio_get_iosystem_from_id(iosysid);
    if(ios == NULL)
	return PIO_EBADID;

  /* if (ios->async_interface && ! ios->ioproc){ */
    /* 	if (ios->compmaster) */
    /* 	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm); */
    /* 	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm); */
    /* } */

    switch (iotype)
    {
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
	ierr = nc_set_chunk_cache(size, nelems, preemption);
	break;
    case PIO_IOTYPE_NETCDF4C:
	if (!ios->io_rank)
	    ierr = nc_set_chunk_cache(size, nelems, preemption);
	break;
#endif
    case PIO_IOTYPE_NETCDF:
	ierr = PIO_ENOTNC4;
	break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
	ierr = PIO_ENOTNC4;
	break;
#endif
    default:
	ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }

    /* Propogate error code to all processes. */
    MPI_Bcast(&ierr, 1, MPI_INTEGER, ios->ioroot, ios->my_comm);    

    return ierr;
}    

/**
 * @ingroup PIO_def_var
 * Get current file chunk cache settings from HDF5.
 *
 * This function has no effect on netCDF classic files. Calling this
 * function with iotype of PIO_IOTYPE_PNETCDF or PIO_IOTYPE_NETCDF
 * returns an error.
 *
 * The file chunk cache for HDF5 can be set, and will apply for any
 * files opened or created until the program ends, or the settings are
 * changed again. The cache settings apply only to the open file. They
 * do not persist with the file, and must be set each time the file is
 * opened, before it is opened, if they are to have effect.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * Chunksizes have important performance repercussions. NetCDF
 * attempts to choose sensible chunk sizes by default, but for best
 * performance check chunking against access patterns.
 *
 * @param iotype the iotype of files to be created or opened.
 * @param sizep gets the size of file cache.
 * @param nelemsp gets the number of elements in file cache.
 * @param preemptionp gets the preemption setting for file cache.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_get_chunk_cache(int iosysid, int iotype, PIO_Offset *sizep,
			 PIO_Offset *nelemsp, float *preemptionp)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    ios = pio_get_iosystem_from_id(iosysid);
    if(ios == NULL)
	return PIO_EBADID;

    /* Since this is a property of the running HDF5 instance, not the
     * file, it's not clear if this message passing will apply. For
     * now, comment it out. EJH */
    /* msg = PIO_MSG_INQ_VAR_FLETCHER32; */

    /* if (ios->async_interface && ! ios->ioproc){ */
    /* 	if (ios->compmaster)  */
    /* 	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm); */
    /* 	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm); */
    /* } */

    switch (iotype)
    {
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
	ierr = nc_get_chunk_cache((size_t *)sizep, (size_t *)nelemsp, preemptionp);
	break;
    case PIO_IOTYPE_NETCDF4C:
	if (!ios->io_rank)
	    ierr = nc_get_chunk_cache((size_t *)sizep, (size_t *)nelemsp, preemptionp);
	break;
#endif
    case PIO_IOTYPE_NETCDF:
	ierr = PIO_ENOTNC4;
	break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
	ierr = PIO_ENOTNC4;
	break;
#endif
    default:
	ierr = iotype_error(iotype,__FILE__,__LINE__);
    }


    /* Check for netCDF error. */
    MPI_Bcast(&ierr, 1, MPI_INTEGER, ios->ioroot, ios->my_comm);
    if (!ierr)
    {
	if (sizep)
	    if ((ierr = MPI_Bcast(sizep, 1, MPI_OFFSET, ios->ioroot, ios->my_comm)))
		ierr = PIO_EIO;
	if (nelemsp && !ierr)
	    if ((ierr = MPI_Bcast(nelemsp, 1, MPI_OFFSET, ios->ioroot, ios->my_comm)))
		ierr = PIO_EIO;
	if (preemptionp && !ierr)
	    if ((ierr = MPI_Bcast(preemptionp, 1, MPI_FLOAT, ios->ioroot, ios->my_comm)))
		ierr = PIO_EIO;
    }

    return ierr;
}

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * Chunksizes have important performance repercussions. NetCDF
 * attempts to choose sensible chunk sizes by default, but for best
 * performance check chunking against access patterns.
 *
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param storage NC_CONTIGUOUS or NC_CHUNKED.
 * @param chunksizep an array of chunksizes. Must have a chunksize for
 * every variable dimension.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_set_var_chunk_cache(int ncid, int varid, PIO_Offset size, PIO_Offset nelems,
			     float preemption)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_SET_VAR_CHUNK_CACHE;

    if (ios->async_interface && ! ios->ioproc)
    {
	if (ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_set_var_chunk_cache(file->fh, varid, size, nelems, preemption);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (!ios->io_rank)
		ierr = nc_set_var_chunk_cache(file->fh, varid, size, nelems, preemption);
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Allocate an error string if needed. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check for netCDF error. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    return ierr;
}    

/**
 * @ingroup PIO_inq_var
 * Get the variable chunk cache settings.
 *
 * This function only applies to netCDF-4 files. When used with netCDF
 * classic files, the error PIO_ENOTNC4 will be returned.
 *
 * Note that these settings are not part of the data file - they apply
 * only to the open file as long as it is open.
 *
 *  See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param sizep will get the size of the cache in bytes.
 * @param nelemsp will get the number of elements in the cache.
 * @param preemptionp will get the cache preemption value.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_get_var_chunk_cache(int ncid, int varid, PIO_Offset *sizep, PIO_Offset *nelemsp,
			     float *preemptionp)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* Since this is a property of the running HDF5 instance, not the
     * file, it's not clear if this message passing will apply. For
     * now, comment it out. EJH */
    /* msg = PIO_MSG_INQ_VAR_FLETCHER32; */

    /* if (ios->async_interface && ! ios->ioproc){ */
    /* 	if (ios->compmaster)  */
    /* 	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm); */
    /* 	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm); */
    /* } */

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_get_var_chunk_cache(file->fh, varid,  (size_t *)sizep, (size_t *)nelemsp,
					  preemptionp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (ios->io_rank == 0)
		ierr = nc_get_var_chunk_cache(file->fh, varid,  (size_t *)sizep, (size_t *)nelemsp,
					      preemptionp);
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = PIO_ENOTNC4;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* If there is an error, allocate space for the error string. */
    if (ierr != PIO_NOERR)
    {
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }

    /* Check the netCDF return code, and broadcast it to all tasks. */
    ierr = check_netcdf(file, ierr, errstr,__LINE__);

    /* Free the error string if it was allocated. */
    if (errstr != NULL)
	free(errstr);

    /* Broadcast results to all tasks. */
    if (sizep && !ierr)
	ierr = MPI_Bcast(sizep, 1, MPI_OFFSET, ios->ioroot, ios->my_comm);	
    if (nelemsp && !ierr)
    	ierr = MPI_Bcast(nelemsp, 1, MPI_OFFSET, ios->ioroot, ios->my_comm);
    if (preemptionp && !ierr)
    	ierr = MPI_Bcast(preemptionp, 1, MPI_FLOAT, ios->ioroot, ios->my_comm);

    return ierr;
}    

	    
