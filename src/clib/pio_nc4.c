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
	    if (!ios->io_rank)
		ierr = nc_def_var_deflate(file->fh, varid, shuffle, deflate, deflate_level);
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
	    ierr = nc_inq_var_deflate(file->fh, varid, shufflep, deflatep, deflate_levelp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if(ios->io_rank==0){
		ierr = nc_inq_var_deflate(file->fh, varid, shufflep, deflatep, deflate_levelp);
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
 * Inquire about szip settings for a variable.
 *
 * Szip is a read-only compression format in netCDF-4. Only raw HDF5
 * can create szip files, but netcdf-4 can read them.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 * 
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param options_maskp will get the options mask.
 * @param pixels_per_block will get the pixels per block.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_inq_var_szip(int ncid, int varid, int *options_maskp,
		      int *pixels_per_blockp)
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
    msg = PIO_MSG_INQ_VAR_SZIP;

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
	    ierr = nc_inq_var_szip(file->fh, varid, options_maskp, pixels_per_blockp);
		break;
	case PIO_IOTYPE_NETCDF4C:
	    if (!ios->io_rank)
		ierr = nc_inq_var_szip(file->fh, varid, options_maskp, pixels_per_blockp);
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
 * Set the fletcher32 filter for a variable.
 * 
 * The fletcher32 filter may help with compression of integer values.
 *
 * See the <a
 * href="http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html">netCDF
 * variable documentation</a> for details about the operation of this
 * function.
 *
 * @param ncid the ncid of the open file.
 * @param varid the ID of the variable to set chunksizes for.
 * @param fletcher32 non-zero to turn on fletcher32 filter.
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
int PIOc_inq_var_fletcher32(int ncid, int varid, int *fletcher32p)
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
    msg = PIO_MSG_INQ_VAR_FLETCHER32;

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
	    ierr = nc_inq_var_fletcher32(file->fh, varid, fletcher32p);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (!ios->io_rank)
		ierr = nc_inq_var_fletcher32(file->fh, varid, fletcher32p);
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
	    ierr = nc_inq_var_endian(file->fh, varid, endianp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (!ios->io_rank)
		ierr = nc_inq_var_endian(file->fh, varid, endianp);
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
 *
 * Set chunk cache netCDF files to be opened/created.
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
 * @param iotype the iotype of files to be created or opened.
`* @param io_rank the rank of the calling process.
 * @param size size of file cache.
 * @param nelems number of elements in file cache.
 * @param preemption preemption setting for file cache.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_set_chunk_cache(int iotype, int io_rank, size_t size,
			 size_t nelems, float preemption)
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;
    /*msg = PIO_MSG_SET_CHUNK_CACHE;*/

    /* if(ios->async_interface && ! ios->ioproc){ */
    /* 	if(ios->compmaster)  */
    /* 	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm); */
    /* 	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm); */
    /* } */

    switch(iotype)
    {
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
`* @param io_rank the rank of the calling process.
 * @param sizep gets the size of file cache.
 * @param nelemsp gets the number of elements in file cache.
 * @param preemptionp gets the preemption setting for file cache.
 * 
 * @return PIO_NOERR for success, otherwise an error code.
 */
int PIOc_get_chunk_cache(int iotype, int io_rank, size_t *sizep,
			 size_t *nelemsp, float *preemptionp)
{
    int ierr;
    int msg;
    int mpierr;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    /* Since this is a property of the running HDF5 instance, not the
     * file, it's not clear if this message passing will apply. For
     * now, comment it out. EJH */
    /* msg = PIO_MSG_INQ_VAR_FLETCHER32; */

    /* if(ios->async_interface && ! ios->ioproc){ */
    /* 	if(ios->compmaster)  */
    /* 	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm); */
    /* 	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm); */
    /* } */

    switch(iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
	ierr = nc_get_chunk_cache(sizep, nelemsp, preemptionp);
	break;
    case PIO_IOTYPE_NETCDF4C:
	if (!io_rank)
	    ierr = nc_get_chunk_cache(sizep, nelemsp, preemptionp);
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
	ierr = iotype_error(iotype,__FILE__,__LINE__);
    }

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

/**
 * @ingroup PIO_inq_var
 * Get the variable chunk cache settings.
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
int PIOc_get_var_chunk_cache(int ncid, int varid, size_t *sizep, size_t *nelemsp,
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

    /* if(ios->async_interface && ! ios->ioproc){ */
    /* 	if(ios->compmaster)  */
    /* 	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm); */
    /* 	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm); */
    /* } */

    if(ios->ioproc){
	switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_get_var_chunk_cache(file->fh, varid,  sizep, nelemsp, preemptionp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (!ios->io_rank)
		ierr = nc_get_var_chunk_cache(file->fh, varid,  sizep, nelemsp, preemptionp);
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

	    
