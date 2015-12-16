/** @file
 *
 * Functions to wrap netCDF-4 functions for PIO.
 **/
#include <pio.h>
#include <pio_internal.h>

/**
 * @ingroup PIO_def_var
 * Set chunksizes for a variable.
 * 
 * Chunksizes have important performance repercussions. NetCDF
 * attempts to choose sensible chunk sizes by default, but for best
 * performance check chunking against access patterns.
 *
 * See the discussion of chunksizes in this paper from the HDF5
 * developers: <a
 * href="https://www.hdfgroup.org/HDF5/PHDF5/parallelhdf5hints.pdf">Parallel
 * HDF5 Hints</a>
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
int nc_inq_var_chunking(int ncid, int varid, int *storagep, size_t *chunksizesp)
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

	    
