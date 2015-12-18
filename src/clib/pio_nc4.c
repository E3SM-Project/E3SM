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
	    ierr = nc_def_var_deflate(file->fh, varid, shuffle, deflate, deflate_level);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
		ierr = nc_def_var_deflate(file->fh, varid, shuffle, deflate, deflate_level);
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

    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}    

/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
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
int PIOc_inq_var_deflate(int ncid, int varid, int *shufflep,
			 int *deflatep, int *deflate_levelp);

/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
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
int PIOc_inq_var_szip(int ncid, int varid, int *options_maskp, int *pixels_per_blockp);

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
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
int PIOc_def_var_fletcher32(int ncid, int varid, int fletcher32)
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
    msg = PIO_MSG_DEF_VAR_FLETCHER32;

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
	    ierr = nc_def_var_fletcher32(file->fh, varid, fletcher32);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
		ierr = nc_def_var_fletcher32(file->fh, varid, fletcher32);
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

    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}

/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
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
int PIOc_inq_var_fletcher32(int ncid, int varid, int *fletcher32p);

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
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
			  const size_t *chunksizesp) 
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
	    ierr = nc_def_var_chunking(file->fh, varid, storage, chunksizesp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
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

    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}

/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
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
int PIOc_inq_var_chunking(int ncid, int varid, int *storagep, size_t *chunksizesp)
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
	    ierr = nc_inq_var_chunking(file->fh, varid, storagep, chunksizesp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
		ierr = nc_inq_var_chunking(file->fh, varid, storagep, chunksizesp);
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

    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
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
	    ierr = nc_def_var_fill(file->fh, varid, no_fill, fill_value);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
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

    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}


/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
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
int PIOc_inq_var_fill(int ncid, int varid, int *no_fill, void *fill_valuep);

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
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
	    ierr = nc_def_var_endian(file->fh, varid, endian);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
		ierr = nc_def_var_endian(file->fh, varid, endian);
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

    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}

/**
 * @ingroup PIO_inq_var
 * Inquire about chunksizes for a variable.
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
int PIOc_inq_var_endian(int ncid, int varid, int *endianp);

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
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
int PIOc_set_chunk_cache(int iotype, int io_rank, size_t size, size_t nelems, float preemption)
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

    /* if(ios->async_interface && ! ios->ioproc){ */
    /* 	if(ios->compmaster)  */
    /* 	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm); */
    /* 	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm); */
    /* } */

    switch(iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
	ierr = nc_set_chunk_cache(size, nelems, preemption);
	break;
    case PIO_IOTYPE_NETCDF4C:
	if (!io_rank)
	    ierr = nc_set_chunk_cache(size, nelems, preemption);
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

    /* ierr = check_netcdf(file, ierr, errstr,__LINE__); */
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}    

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
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
int PIOc_get_chunk_cache(size_t *sizep, size_t *nelemsp, float *preemptionp);

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
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
int PIOc_set_var_chunk_cache(int ncid, int varid, size_t size, size_t nelems,
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
	    ierr = nc_set_var_chunk_cache(file->fh, varid, size, nelems, preemption);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
		ierr = nc_set_var_chunk_cache(file->fh, varid, size, nelems, preemption);
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

    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    if(errstr != NULL) free(errstr);
    return ierr;
}    

int PIOc_get_var_chunk_cache(int ncid, int varid, size_t *sizep, size_t *nelemsp,
			     float *preemptionp);
	    
