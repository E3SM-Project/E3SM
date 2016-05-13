/**
 * @file   
 * PIO interfaces to
 * [NetCDF](http://www.unidata.ucar.edu/software/netcdf/docs/modules.html)
 * support functions

 *  This file provides an interface to the
 *  [NetCDF](http://www.unidata.ucar.edu/software/netcdf/docs/modules.html)
 *  support functions.  Each subroutine calls the underlying netcdf or
 *  pnetcdf or netcdf4 functions from the appropriate subset of mpi
 *  tasks (io_comm). Each routine must be called collectively from
 *  union_comm.
 *  
 * @author Jim Edwards (jedwards@ucar.edu), Ed Hartnett
 * @date     Feburary 2014, April 2016
 */

#include <config.h>
#ifdef PIO_ENABLE_LOGGING
#include <stdarg.h>
#include <unistd.h>
#endif /* PIO_ENABLE_LOGGING */
#include <pio.h>
#include <pio_internal.h>

#ifdef PIO_ENABLE_LOGGING
int pio_log_level = 0;
int my_rank;
#endif /* PIO_ENABLE_LOGGING */

/** Set the logging level. Set to -1 for nothing, 0 for errors only, 1
 * for important logging, and so on. Log levels below 1 are only
 * printed on the io/component root. If the library is not built with
 * logging, this function does nothing. */
int PIOc_set_log_level(int level)
{
#ifdef PIO_ENABLE_LOGGING
    printf("setting log level to %d\n", level);
    pio_log_level = level;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    return PIO_NOERR;
#endif /* PIO_ENABLE_LOGGING */
}

#ifdef PIO_ENABLE_LOGGING
/** This function prints out a message, if the severity of the message
   is lower than the global pio_log_level. To use it, do something
   like this:
   
   pio_log(0, "this computer will explode in %d seconds", i);

   After the first arg (the severity), use the rest like a normal
   printf statement. Output will appear on stdout.
   This function is heavily based on the function in section 15.5 of
   the C FAQ. 
*/
void 
pio_log(int severity, const char *fmt, ...)
{
   va_list argp;
   int t;

   /* If the severity is greater than the log level, we don't print
      this message. */
   if (severity > pio_log_level)
      return;

   /* If the severity is 0, only print on rank 0. */
   if (severity < 1 && my_rank != 0)
       return;

   /* If the severity is zero, this is an error. Otherwise insert that
      many tabs before the message. */
   if (!severity)
       fprintf(stdout, "ERROR: ");
   for (t = 0; t < severity; t++)
       fprintf(stdout, "\t");

   /* Show the rank. */
   fprintf(stdout, "%d ", my_rank);
   
   /* Print out the variable list of args with vprintf. */
   va_start(argp, fmt);
   vfprintf(stdout, fmt, argp);
   va_end(argp);
   
   /* Put on a final linefeed. */
   fprintf(stdout, "\n");
   fflush(stdout);
}
#endif /* PIO_ENABLE_LOGGING */

/** 
 * @ingroup PIOc_inq
 * The PIO-C interface for the NetCDF function nc_inq.
 *
 * This routine is called collectively by all tasks in the
 * communicator ios.union_comm. For more information on the underlying
 * NetCDF commmand please read about this function in the NetCDF
 * documentation at:
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 *
 * @return PIO_NOERR for success, error code otherwise. See
 * PIOc_Set_File_Error_Handling
 */
int PIOc_inq(int ncid, int *ndimsp, int *nvarsp, int *ngattsp,
	     int *unlimdimidp) 
{
    int msg = PIO_MSG_INQ; /** Message for async notification. */
    iosystem_desc_t *ios;  /** Pointer to io system information. */
    file_desc_t *file;     /** Pointer to file information. */
    int ierr = PIO_NOERR;  /** Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */

    LOG((1, "PIOc_inq ncid = %d", ncid));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    char ndims_present = ndimsp ? true : false;
	    char nvars_present = nvarsp ? true : false;
	    char ngatts_present = ngattsp ? true : false;
	    char unlimdimid_present = unlimdimidp ? true : false;
	    
	    if (ios->compmaster) 
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	    
	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&ndims_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&nvars_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&ngatts_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&unlimdimid_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    LOG((2, "PIOc_inq netcdf Bcast unlimdimid_present = %d", unlimdimid_present));
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if (!file->iosystem->io_rank)
	    {
		/* Should not be necessary to do this - nc_inq should
		 * handle null pointers. This has been reported as a bug
		 * to netCDF developers. */
		int tmp_ndims, tmp_nvars, tmp_ngatts, tmp_unlimdimid;
		ierr = nc_inq(ncid, &tmp_ndims, &tmp_nvars, &tmp_ngatts, &tmp_unlimdimid);
		if (ndimsp)
		    *ndimsp = tmp_ndims;
		if (nvarsp)
		    *nvarsp = tmp_nvars;
		if (ngattsp)
		    *ngattsp = tmp_ngatts;
		if (unlimdimidp)
		    *unlimdimidp = tmp_unlimdimid;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp);
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}

	LOG((2, "PIOc_inq netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	return PIO_EIO;
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (!ierr)
    {
	if (ndimsp)
	    if ((mpierr = MPI_Bcast(ndimsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	
	if (nvarsp)
	    if ((mpierr = MPI_Bcast(nvarsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	
	if (ngattsp)
	    if ((mpierr = MPI_Bcast(ngattsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	
	if (unlimdimidp)
	    if ((mpierr = MPI_Bcast(unlimdimidp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    return ierr;
}

/** 
 * @ingroup PIOc_inq_ndims
 * The PIO-C interface for the NetCDF function nc_inq_ndims.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_ndims (int ncid, int *ndimsp) 
{
    LOG((1, "PIOc_inq_ndims"));
    return PIOc_inq(ncid, ndimsp, NULL, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_nvars
 * The PIO-C interface for the NetCDF function nc_inq_nvars.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_nvars(int ncid, int *nvarsp) 
{
    return PIOc_inq(ncid, NULL, nvarsp, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_natts
 * The PIO-C interface for the NetCDF function nc_inq_natts.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_natts(int ncid, int *ngattsp) 
{
    return PIOc_inq(ncid, NULL, NULL, ngattsp, NULL);
}

/** 
 * @ingroup PIOc_inq_unlimdim
 * The PIO-C interface for the NetCDF function nc_inq_unlimdim.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_unlimdim(int ncid, int *unlimdimidp) 
{
    return PIOc_inq(ncid, NULL, NULL, unlimdimidp, NULL);
}

/** 
 * @ingroup PIOc_typelen
 * The PIO-C interface for the NetCDF function nctypelen.
 */
int PIOc_inq_type(int ncid, nc_type xtype, char *name, PIO_Offset *sizep)
{
    int msg = PIO_MSG_INQ_TYPE; /** Message for async notification. */
    iosystem_desc_t *ios;  /** Pointer to io system information. */
    file_desc_t *file;     /** Pointer to file information. */
    int ierr = PIO_NOERR;  /** Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    int typelen;

    LOG((1, "PIOc_inq_type ncid = %d xtype = %d", ncid, xtype));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    char name_present = name ? true : false;	
	    char size_present = sizep ? true : false;
	    
	    if (ios->compmaster) 
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	    
	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&xtype, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&size_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_type(ncid, xtype, name, (size_t *)sizep);
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if (!file->iosystem->io_rank)
		ierr = nc_inq_type(ncid, xtype, name, (size_t *)sizep);		
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    switch (xtype)
	    {
	    case NC_UBYTE:
	    case NC_BYTE:
	    case NC_CHAR:
		typelen = 1;
		break;
	    case NC_SHORT:
	    case NC_USHORT:
		typelen = 2;
		break;
	    case NC_UINT:
	    case NC_INT:
	    case NC_FLOAT:
		typelen = 4;
		break;
	    case NC_UINT64:
	    case NC_INT64:
	    case NC_DOUBLE:
		typelen = 8;
		break;
	    }
	    
	    if (sizep)
		*sizep = typelen;
	    if (name)
		strcpy(name, "some type");
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}

	LOG((2, "PIOc_inq_type netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	return PIO_EIO;
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (!ierr)
    {
	if (name)
	{ 
	    int slen;
	    if (ios->iomaster)
		slen = strlen(name);
	    if ((mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	    if ((mpierr = MPI_Bcast((void *)name, slen + 1, MPI_CHAR, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	}
	if (sizep)
	    if ((mpierr = MPI_Bcast(sizep , 1, MPI_OFFSET, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);	    
    }

    return ierr;
}

/** 
 * @ingroup PIOc_inq_format
 * The PIO-C interface for the NetCDF function nc_inq_format.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param formatp a pointer that will get the file format 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_format (int ncid, int *formatp) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_INQ_FORMAT;

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
	    ierr = nc_inq_format(file->fh, formatp);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_inq_format(file->fh, formatp);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_format(file->fh, formatp);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(formatp , 1, MPI_INT, ios->ioroot, ios->my_comm);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_inq_dim
 * The PIO-C interface for the NetCDF function nc_inq_dim.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param lenp a pointer that will get the number of values 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_dim(int ncid, int dimid, char *name, PIO_Offset *lenp) 
{
    iosystem_desc_t *ios;
    file_desc_t *file;
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */

    LOG((1, "PIOc_inq_dim"));

    /* Get the file info, based on the ncid. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_INQ_DIM;
	    char name_present = name ? true : false;	
	    char len_present = lenp ? true : false;
	    
	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	    
	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&dimid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    LOG((2, "PIOc_inq netcdf Bcast name_present = %d", name_present));
	    if (!mpierr)
		mpierr = MPI_Bcast(&len_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    LOG((2, "PIOc_inq netcdf Bcast len_present = %d", len_present));
	}
	
	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* Make the call to the netCDF layer. */
    if(ios->ioproc){
	switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_dim(file->fh, dimid, name, (size_t *)lenp);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if (ios->io_rank == 0){
		ierr = nc_inq_dim(file->fh, dimid, name, (size_t *)lenp);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_dim(file->fh, dimid, name, lenp);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	return PIO_EIO;
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (!ierr)
    {
	if (name)
	{ 
	    int slen;
	    if (ios->iomaster)
		slen = strlen(name);
	    if ((mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	    if ((mpierr = MPI_Bcast((void *)name, slen + 1, MPI_CHAR, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	}
	
	if (lenp)
	    if ((mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);	    
    }
    
    return ierr;
}

/** 
 * @ingroup PIOc_inq_dimname
 * The PIO-C interface for the NetCDF function nc_inq_dimname.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_dimname(int ncid, int dimid, char *name) 
{
    return PIOc_inq_dim(ncid, dimid, name, NULL);
}

/** 
 * @ingroup PIOc_inq_dimlen
 * The PIO-C interface for the NetCDF function nc_inq_dimlen.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param lenp a pointer that will get the number of values 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_dimlen(int ncid, int dimid, PIO_Offset *lenp) 
{
    return PIOc_inq_dim(ncid, dimid, NULL, lenp);
}

/** 
 * @ingroup PIOc_inq_dimid
 * The PIO-C interface for the NetCDF function nc_inq_dimid.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param idp a pointer that will get the id of the variable or attribute.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_dimid(int ncid, const char *name, int *idp) 
{
    iosystem_desc_t *ios;
    file_desc_t *file;
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */

    /* Name must be provided. */
    if (!name)
	return PIO_EINVAL;
    
    LOG((1, "PIOc_inq_dimid name = %s", name));

    /* Get the file info, based on the ncid. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If using async, and not an IO task, then send parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_INQ_DIMID;
	    char id_present = idp ? true : false;
	    
	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	    
	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    int namelen = strlen(name);
	    if (!mpierr)
		mpierr = MPI_Bcast(&namelen, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&id_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* IO tasks call the netCDF functions. */
    if (ios->ioproc)
    {
	switch(file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_dimid(file->fh, name, idp);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if (ios->io_rank == 0)
		ierr = nc_inq_dimid(file->fh, name, idp);;
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_dimid(file->fh, name, idp);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	check_mpi(file, mpierr, __FILE__, __LINE__);
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results. */
    if (!ierr)
	if (idp)
	    if ((mpierr = MPI_Bcast(idp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);

    return ierr;
}

/** 
 * @ingroup PIOc_inq_var
 * The PIO-C interface for the NetCDF function nc_inq_var.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param xtypep a pointer that will get the type of the attribute.
 * @param nattsp a pointer that will get the number of attributes 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_var(int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp,
		 int *dimidsp, int *nattsp) 
{
    iosystem_desc_t *ios;
    file_desc_t *file;
    int ndims;    /** The number of dimensions for this variable. */
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */

    LOG((1, "PIOc_inq_var ncid = %d varid = %d", ncid, varid));

    /* Get the file info, based on the ncid. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_INQ_VAR;
	    char name_present = name ? true : false;
	    char xtype_present = xtypep ? true : false;
	    char ndims_present = ndimsp ? true : false;
	    char dimids_present = dimidsp ? true : false;
	    char natts_present = nattsp ? true : false;

	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&xtype_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&ndims_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&dimids_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&natts_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	LOG((2, "PIOc_inq_var name_present = %d xtype_present = %d ndims_present = %d "
	     "dimids_present = %d, natts_present = %d nattsp = %d",
	     name_present, xtype_present, ndims_present, dimids_present, natts_present, nattsp));
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }
	
    /* Call the netCDF layer. */
    if (ios->ioproc)
    {
	switch(file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_varndims(file->fh, varid, &ndims);
	    ierr = nc_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if (ios->io_rank == 0)
	    {
		ierr = nc_inq_varndims(file->fh, varid, &ndims);
		ierr = nc_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_varndims(file->fh, varid, &ndims);
	    ierr = ncmpi_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	check_mpi(file, mpierr, __FILE__, __LINE__);
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast the results for non-null pointers. */
    if (!ierr)
    {
	if (name)
	{ 
	    int slen;
	    if(ios->iomaster)
		slen = strlen(name);
	    if ((mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	    if ((mpierr = MPI_Bcast((void *)name, slen + 1, MPI_CHAR, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	}
	if (xtypep)
	    if ((mpierr = MPI_Bcast(xtypep, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	
	if (ndimsp)
	{
	    if ((mpierr = MPI_Bcast(ndimsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);	    
	    file->varlist[varid].ndims = (*ndimsp);
	}
	if (dimidsp)
	{
	    if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	    if ((mpierr = MPI_Bcast(dimidsp, ndims, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	}
	if (nattsp)
	    if ((mpierr = MPI_Bcast(nattsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    return ierr;
}

/** 
 * @ingroup PIOc_inq_varname
 * The PIO-C interface for the NetCDF function nc_inq_varname.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_varname (int ncid, int varid, char *name) 
{
    return PIOc_inq_var(ncid, varid, name, NULL, NULL, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_vartype
 * The PIO-C interface for the NetCDF function nc_inq_vartype.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param xtypep a pointer that will get the type of the attribute.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_vartype (int ncid, int varid, nc_type *xtypep) 
{
    return PIOc_inq_var(ncid, varid, NULL, xtypep, NULL, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_varndims
 * The PIO-C interface for the NetCDF function nc_inq_varndims.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_varndims (int ncid, int varid, int *ndimsp) 
{
    return PIOc_inq_var(ncid, varid, NULL, NULL, ndimsp, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_vardimid
 * The PIO-C interface for the NetCDF function nc_inq_vardimid.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_vardimid(int ncid, int varid, int *dimidsp) 
{
    return PIOc_inq_var(ncid, varid, NULL, NULL, NULL, dimidsp, NULL);
}

/** 
 * @ingroup PIOc_inq_varnatts
 * The PIO-C interface for the NetCDF function nc_inq_varnatts.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param nattsp a pointer that will get the number of attributes 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_varnatts (int ncid, int varid, int *nattsp) 
{
    return PIOc_inq_var(ncid, varid, NULL, NULL, NULL, NULL, nattsp);
}

/** 
 * @ingroup PIOc_inq_varid
 * The PIO-C interface for the NetCDF function nc_inq_varid.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param varidp a pointer that will get the variable id 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_varid (int ncid, const char *name, int *varidp) 
{
    iosystem_desc_t *ios;
    file_desc_t *file;
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */

    /* Caller must provide name. */
    if (!name)
	return PIO_EINVAL;

    /* Get file info based on ncid. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    LOG((1, "PIOc_inq_varid ncid = %d name = %s", ncid, name));

    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_INQ_VARID;
	    
	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	    
	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    int namelen;
	    namelen = strlen(name);
	    if (!mpierr)
		mpierr = MPI_Bcast(&namelen, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
	switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_varid(file->fh, name, varidp);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_inq_varid(file->fh, name, varidp);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_varid(file->fh, name, varidp);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	check_mpi(file, mpierr, __FILE__, __LINE__);
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (varidp)
	if ((mpierr = MPI_Bcast(varidp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr, __FILE__, __LINE__);

    return ierr;
}

/** 
 * @ingroup PIOc_inq_att
 * The PIO-C interface for the NetCDF function nc_inq_att.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param xtypep a pointer that will get the type of the attribute.
 * @param lenp a pointer that will get the number of values 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_att(int ncid, int varid, const char *name, nc_type *xtypep,
		 PIO_Offset *lenp) 
{
    int msg = PIO_MSG_INQ_ATT;
    iosystem_desc_t *ios;
    file_desc_t *file;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    int ierr = PIO_NOERR;

    /* Caller must provide a name. */
    if (!name)
	return PIO_EINVAL;

    LOG((1, "PIOc_inq_att ncid = %d varid = %d xtpyep = %d lenp = %d",
	 ncid, varid, xtypep, lenp));

    /* Find file based on ncid. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If using async and this is not an IO task, send parameter data
     * over the intercomm. */
    if(ios->async_interface && !ios->ioproc)
    {
	char xtype_present = xtypep ? true : false;
	char len_present = lenp ? true : false;

	if(ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	int namelen = strlen(name);
	mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&xtype_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&len_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
    }

    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_inq_att(file->fh, varid, name, xtypep, (size_t *)lenp);
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if (ios->io_rank == 0)
		ierr = nc_inq_att(file->fh, varid, name, xtypep, (size_t *)lenp);
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_att(file->fh, varid, name, xtypep, lenp);
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Handle MPI errors. */
    if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	check_mpi(file, mpierr2, __FILE__, __LINE__);	    
    check_mpi(file, mpierr, __FILE__, __LINE__);

    /* Broadcast results. */
    if (!ierr)
    {
	if(xtypep)
	    if ((mpierr = MPI_Bcast(xtypep, 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
	if(lenp)
	    if ((mpierr = MPI_Bcast(lenp, 1, MPI_OFFSET, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    return ierr;
}

/** 
 * @ingroup PIOc_inq_attlen
 * The PIO-C interface for the NetCDF function nc_inq_attlen.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param lenp a pointer that will get the number of values 
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_attlen (int ncid, int varid, const char *name, PIO_Offset *lenp) 
{
    return PIOc_inq_att(ncid, varid, name, NULL, lenp);
}

/** 
 * @ingroup PIOc_inq_atttype
 * The PIO-C interface for the NetCDF function nc_inq_atttype.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param xtypep a pointer that will get the type of the attribute.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_atttype(int ncid, int varid, const char *name, nc_type *xtypep) 
{
    return PIOc_inq_att(ncid, varid, name, xtypep, NULL);
}

/** 
 * @ingroup PIOc_inq_attname
 * The PIO-C interface for the NetCDF function nc_inq_attname.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param attnum the attribute ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_attname (int ncid, int varid, int attnum, char *name) 
{
    int ierr = PIO_NOERR;
    int msg = PIO_MSG_INQ_ATTNAME;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    iosystem_desc_t *ios;
    file_desc_t *file;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

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
	    ierr = nc_inq_attname(file->fh, varid, attnum, name);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_inq_attname(file->fh, varid, attnum, name);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_attname(file->fh, varid, attnum, name);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    check_netcdf(file, ierr, __FILE__, __LINE__);
    if (name)
    {
	int slen;
	if(ios->iomaster)
	    slen = (int) strlen(name);
	mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);
	mpierr = MPI_Bcast((void *)name, slen + 1, MPI_CHAR, ios->ioroot, ios->my_comm);
    }

    return ierr;
}

/** 
 * @ingroup PIOc_inq_attid
 * The PIO-C interface for the NetCDF function nc_inq_attid.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param idp a pointer that will get the id of the variable or attribute.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_attid (int ncid, int varid, const char *name, int *idp) 
{
    int ierr = PIO_NOERR;
    int msg = PIO_MSG_INQ_ATTID;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    iosystem_desc_t *ios;
    file_desc_t *file;

    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

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
	    ierr = nc_inq_attid(file->fh, varid, name, idp);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_inq_attid(file->fh, varid, name, idp);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_attid(file->fh, varid, name, idp);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Handle MPI errors. */
    if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	check_mpi(file, mpierr2, __FILE__, __LINE__);	    
    check_mpi(file, mpierr, __FILE__, __LINE__);

    /* Broadcast results. */
    if (!ierr)
    {
	if (idp)
	    if ((mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm)))
		check_mpi(file, mpierr, __FILE__, __LINE__);	
    }

    return ierr;
}

/** 
 * @ingroup PIOc_rename_dim
 * The PIO-C interface for the NetCDF function nc_rename_dim.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_rename_dim (int ncid, int dimid, const char *name) 
{
    int ierr;
    int msg;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_RENAME_DIM;

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
	    ierr = nc_rename_dim(file->fh, dimid, name);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_rename_dim(file->fh, dimid, name);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_rename_dim(file->fh, dimid, name);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_rename_var
 * The PIO-C interface for the NetCDF function nc_rename_var.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_rename_var (int ncid, int varid, const char *name) 
{
    int ierr;
    int msg;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_RENAME_VAR;

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
	    ierr = nc_rename_var(file->fh, varid, name);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_rename_var(file->fh, varid, name);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_rename_var(file->fh, varid, name);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_rename_att
 * The PIO-C interface for the NetCDF function nc_rename_att.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_rename_att (int ncid, int varid, const char *name, const char *newname) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_RENAME_ATT;

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
	    ierr = nc_rename_att(file->fh, varid, name, newname);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_rename_att(file->fh, varid, name, newname);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_rename_att(file->fh, varid, name, newname);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_del_att
 * The PIO-C interface for the NetCDF function nc_del_att.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_del_att (int ncid, int varid, const char *name) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_DEL_ATT;

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
	    ierr = nc_del_att(file->fh, varid, name);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_del_att(file->fh, varid, name);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_del_att(file->fh, varid, name);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_set_fill
 * The PIO-C interface for the NetCDF function nc_set_fill.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_set_fill (int ncid, int fillmode, int *old_modep) 
{
    int ierr;
    int msg;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_SET_FILL;

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
	    ierr = nc_set_fill(file->fh, fillmode, old_modep);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_set_fill(file->fh, fillmode, old_modep);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_set_fill(file->fh, fillmode, old_modep);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_enddef
 * The PIO-C interface for the NetCDF function nc_enddef.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_enddef(int ncid) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_ENDDEF;

    if(ios->async_interface && ! ios->ioproc){
	if(!ios->comp_rank) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	printf("PIOc_enddef file->fh = %d\n", file->fh);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
    }


    if(ios->ioproc){
	switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_enddef(file->fh);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_enddef(file->fh);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_enddef(file->fh);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_redef
 * The PIO-C interface for the NetCDF function nc_redef.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_redef (int ncid) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_REDEF;

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
	    ierr = nc_redef(file->fh);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_redef(file->fh);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_redef(file->fh);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_def_dim
 * The PIO-C interface for the NetCDF function nc_def_dim.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param idp a pointer that will get the id of the variable or attribute.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_def_dim (int ncid, const char *name, PIO_Offset len, int *idp) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;
    int namelen;

    errstr = NULL;
    ierr = PIO_NOERR;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d PIOc_def_dim ncid = %d name = %s len = %d\n", my_rank,
	   ncid, name, len);

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_DEF_DIM;

    if(ios->async_interface && ! ios->ioproc){
	if(ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
	namelen = strlen(name);
	mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&len, 1, MPI_INT,  ios->compmaster, ios->intercomm);
    }

    if(ios->ioproc){
	switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_def_dim(file->fh, name, (size_t)len, idp);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_def_dim(file->fh, name, (size_t)len, idp);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_def_dim(file->fh, name, len, idp);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_def_var
 * The PIO-C interface for the NetCDF function nc_def_var.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @param varidp a pointer that will get the variable id 
 * @return PIO_NOERR for success, error code otherwise. See
 * PIOc_Set_File_Error_Handling
 */
int PIOc_def_var (int ncid, const char *name, nc_type xtype, int ndims,
		  const int *dimidsp, int *varidp) 
{
    iosystem_desc_t *ios;
    file_desc_t *file;
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    int namelen;

    /* User must provide name and storage for varid. */
    if (!name || !varidp)
    {
	check_netcdf(file, PIO_EINVAL, __FILE__, __LINE__);	
	return PIO_EINVAL;
    }

    /* Get the file information. */
    if (!(file = pio_get_file_from_id(ncid)))
    {
	check_netcdf(file, PIO_EBADID, __FILE__, __LINE__);	
	return PIO_EBADID;
    }
    ios = file->iosystem;

    /* If using async, and not an IO task, then send parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_DEF_VAR;
	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&(file->fh), 1, MPI_INT, ios->compmaster, ios->intercomm);
	    namelen = strlen(name);
	    LOG((2, "bcasting namelen = %d name = %s\n", namelen, name));
	    if (!ios->compmaster)
		ios->compmaster = MPI_PROC_NULL;
	    if (!mpierr)
		mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&xtype, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)dimidsp, ndims, MPI_INT, ios->compmaster, ios->intercomm);
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* IO tasks call the netCDF functions. */
    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
	    if (ios->io_rank == 0)
	    {
		ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);
		if (!ierr)
		    ierr = nc_def_var_deflate(file->fh, *varidp, 0,1,1);
	    }
	    break;
#endif
	case PIO_IOTYPE_NETCDF:
	    if (ios->io_rank == 0)
		ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	check_mpi(file, mpierr, __FILE__, __LINE__);
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results. */
    if (!ierr)
	if (varidp)
	    mpierr = MPI_Bcast(varidp , 1, MPI_INT, ios->ioroot, ios->my_comm);

    return ierr;
}

/** 
 * @ingroup PIOc_inq_var_fill
 * The PIO-C interface for the NetCDF function nc_inq_var_fill.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_inq_var_fill (int ncid, int varid, int *no_fill, void *fill_value) 
{
    int ierr;
    int msg;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_INQ_VAR_FILL;

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
	    ierr = nc_inq_var_fill(file->fh, varid, no_fill, fill_value);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_inq_var_fill(file->fh, varid, no_fill, fill_value);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_inq_var_fill(file->fh, varid, no_fill, fill_value);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(fill_value, 1, MPI_INT, ios->ioroot, ios->my_comm);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_get_att
 * The PIO-C interface for the NetCDF function nc_get_att.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise. See
 * PIOc_Set_File_Error_Handling
 */
int PIOc_get_att(int ncid, int varid, const char *name, void *ip) 
{
    iosystem_desc_t *ios;
    file_desc_t *file;
    PIO_Offset attlen;
    nc_type atttype;
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */

    /* User must provide a name and destination pointer. */
    if (!name || !ip)
	return PIO_EINVAL;

    LOG((1, "PIOc_get_att ncid %d varid %d name %s", ncid, varid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_GET_ATT_INT;

	    /* Get the type and length of the attribute. */
	    if ((ierr = PIOc_inq_att(file->fh, varid, name, &atttype, &attlen)))
	    {
		check_netcdf(file, ierr, __FILE__, __LINE__);
		return ierr;
	    }

	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
	    
	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    int namelen = strlen(name);
	    if (!mpierr)
		mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&file->iotype, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&atttype, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&attlen, 1, MPI_OFFSET, ios->compmaster, ios->intercomm);
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }
	
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_get_att(file->fh, varid, name, ip);
	    if (!ierr)
		ierr = nc_inq_att(file->fh, varid, name, &atttype, (size_t *)&attlen);
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if (ios->io_rank == 0)
	    {
		ierr = nc_get_att(file->fh, varid, name, ip);
		if (!ierr)
		    ierr = nc_inq_att(file->fh, varid, name, &atttype, (size_t *)&attlen);
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_get_att(file->fh, varid, name, ip);
	    if (!ierr)
		ierr = ncmpi_inq_att(file->fh, varid, name, &atttype, &attlen);
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }
    LOG((2, "PIOc_get_att called netcdf layer ierr = %d", ierr));

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	return PIO_EIO;
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results to all tasks. */
    if (!ierr)
    {
	LOG((2, "PIOc_get_att broadcasting attlen  = %d", attlen));	
        mpierr = MPI_Bcast(&attlen, 1, MPI_OFFSET, ios->ioroot, ios->my_comm);
	LOG((2, "PIOc_get_att done broadcasting attlen  = %d", attlen));	
	LOG((2, "PIOc_get_att broadcasting att data"));	
        mpierr = MPI_Bcast(ip, (int)attlen, MPI_INT, ios->ioroot, ios->my_comm);
	LOG((2, "PIOc_get_att done broadcasting att data"));	
    }
    return ierr;
}

/** 
 * @ingroup PIOc_put_att
 * The PIO-C interface for the NetCDF function nc_put_att.
 *
 * This routine is called collectively by all tasks in the communicator 
 * ios.union_comm. For more information on the underlying NetCDF commmand
 * please read about this function in the NetCDF documentation at: 
 * http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 *
 * @param ncid the ncid of the open file, obtained from
 * PIOc_openfile() or PIOc_createfile().
 * @param varid the variable ID.
 * @return PIO_NOERR for success, error code otherwise.  See PIOc_Set_File_Error_Handling
 */
int PIOc_put_att(int ncid, int varid, const char *name, nc_type xtype,
		 PIO_Offset len, const void *op) 
{
    iosystem_desc_t *ios;  /** Pointer to io system information. */
    file_desc_t *file;     /** Pointer to file information. */
    int ierr = PIO_NOERR;  /** Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */

    LOG((1, "PIOc_inq_varid ncid = %d name = %s", ncid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
	return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = PIO_MSG_PUT_ATT;

	    if(ios->compmaster) 
		mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    int namelen = strlen(name);
	    if (!mpierr)
		mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&xtype, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast(&len, 1, MPI_OFFSET,  ios->compmaster, ios->intercomm);
	    if (!mpierr)
		mpierr = MPI_Bcast((void *)op, len, MPI_INT,  ios->compmaster, ios->intercomm);
	}

	/* Handle MPI errors. */
	if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	    check_mpi(file, mpierr2, __FILE__, __LINE__);	    
	check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
	switch (file->iotype)
	{
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_put_att(file->fh, varid, name, xtype, (size_t)len, op);
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_put_att(file->fh, varid, name, xtype, (size_t)len, op);
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_put_att(file->fh, varid, name, xtype, len, op);
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
	return PIO_EIO;
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    return ierr;
}

/** 
 * @ingroup PIOc_put_att_short
 * The PIO-C interface for the NetCDF function nc_put_att_short.
 */
int PIOc_put_att_short(int ncid, int varid, const char *name, nc_type xtype,
		       PIO_Offset len, const short *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_get_att_double
 * The PIO-C interface for the NetCDF function nc_get_att_double.
 */
int PIOc_get_att_double(int ncid, int varid, const char *name, double *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_put_att_double
 * The PIO-C interface for the NetCDF function nc_put_att_double.
 */
int PIOc_put_att_double(int ncid, int varid, const char *name, nc_type xtype,
			PIO_Offset len, const double *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_get_att_uchar
 * The PIO-C interface for the NetCDF function nc_get_att_uchar.
 */
int PIOc_get_att_uchar (int ncid, int varid, const char *name, unsigned char *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_put_att_schar
 * The PIO-C interface for the NetCDF function nc_put_att_schar.
 */
int PIOc_put_att_schar(int ncid, int varid, const char *name, nc_type xtype,
		       PIO_Offset len, const signed char *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_get_att_ushort
 * The PIO-C interface for the NetCDF function nc_get_att_ushort.
 */
int PIOc_get_att_ushort (int ncid, int varid, const char *name, unsigned short *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_uint
 * The PIO-C interface for the NetCDF function nc_get_att_uint.
 */
int PIOc_get_att_uint (int ncid, int varid, const char *name, unsigned int *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_long
 * The PIO-C interface for the NetCDF function nc_get_att_long.
 */
int PIOc_get_att_long (int ncid, int varid, const char *name, long *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_ubyte
 * The PIO-C interface for the NetCDF function nc_get_att_ubyte.
 */
int PIOc_get_att_ubyte (int ncid, int varid, const char *name, unsigned char *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_text
 * The PIO-C interface for the NetCDF function nc_get_att_text.
 */
int PIOc_get_att_text (int ncid, int varid, const char *name, char *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_schar
 * The PIO-C interface for the NetCDF function nc_get_att_schar.
 */
int PIOc_get_att_schar (int ncid, int varid, const char *name, signed char *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_ulonglong
 * The PIO-C interface for the NetCDF function nc_get_att_ulonglong.
 */
int PIOc_get_att_ulonglong (int ncid, int varid, const char *name, unsigned long long *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_short
 * The PIO-C interface for the NetCDF function nc_get_att_short.
 */
int PIOc_get_att_short (int ncid, int varid, const char *name, short *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_int
 * The PIO-C interface for the NetCDF function nc_get_att_int.
 */
int PIOc_get_att_int(int ncid, int varid, const char *name, int *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_get_att_longlong
 * The PIO-C interface for the NetCDF function nc_get_att_longlong.
 */
int PIOc_get_att_longlong(int ncid, int varid, const char *name, long long *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_put_att_long
 * The PIO-C interface for the NetCDF function nc_put_att_long.
 */
int PIOc_put_att_long(int ncid, int varid, const char *name, nc_type xtype,
		      PIO_Offset len, const long *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

/** 
 * @ingroup PIOc_get_att_float
 * The PIO-C interface for the NetCDF function nc_get_att_float.
 */
int PIOc_get_att_float (int ncid, int varid, const char *name, float *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
}

/** 
 * @ingroup PIOc_put_att_int
 * The PIO-C interface for the NetCDF function nc_put_att_int.
 */
int PIOc_put_att_int(int ncid, int varid, const char *name, nc_type xtype,
		     PIO_Offset len, const int *op) 
{
    int ierr;
    int msg;
    int mpierr;
    iosystem_desc_t *ios;
    file_desc_t *file;
    char *errstr;
    size_t namelen;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d PIOc_inq_varid ncid = %d name = %s\n", my_rank, ncid, name);

    errstr = NULL;
    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if(file == NULL)
	return PIO_EBADID;
    ios = file->iosystem;
    msg = PIO_MSG_PUT_ATT_INT;

    if(ios->async_interface && ! ios->ioproc){
	if(ios->compmaster) 
	    mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);
	mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
	namelen = strlen(name);
	mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&xtype, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&len, 1, MPI_OFFSET,  ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast((void *)op, len, MPI_INT,  ios->compmaster, ios->intercomm);
    }

    if(ios->ioproc){
	switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
	case PIO_IOTYPE_NETCDF4P:
	    ierr = nc_put_att_int(file->fh, varid, name, xtype, (size_t)len, op);;
	    break;
	case PIO_IOTYPE_NETCDF4C:
#endif
	case PIO_IOTYPE_NETCDF:
	    if(ios->io_rank==0){
		ierr = nc_put_att_int(file->fh, varid, name, xtype, (size_t)len, op);;
	    }
	    break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	    ierr = ncmpi_put_att_int(file->fh, varid, name, xtype, len, op);;
	    break;
#endif
	default:
	    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	}
    }

    if(ierr != PIO_NOERR){
	errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
	sprintf(errstr,"in file %s",__FILE__);
    }
    ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(errstr != NULL) free(errstr);
    return ierr;
}

/** 
 * @ingroup PIOc_put_att_uchar
 * The PIO-C interface for the NetCDF function nc_put_att_uchar.
 */
int PIOc_put_att_uchar(int ncid, int varid, const char *name, nc_type xtype,
		       PIO_Offset len, const unsigned char *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

/** 
 * @ingroup PIOc_put_att_longlong
 * The PIO-C interface for the NetCDF function nc_put_att_longlong.
 */
int PIOc_put_att_longlong(int ncid, int varid, const char *name, nc_type xtype,
			  PIO_Offset len, const long long *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

/** 
 * @ingroup PIOc_put_att_uint
 * The PIO-C interface for the NetCDF function nc_put_att_uint.
 */
int PIOc_put_att_uint(int ncid, int varid, const char *name, nc_type xtype,
		      PIO_Offset len, const unsigned int *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

/** 
 * @ingroup PIOc_put_att_ubyte
 * The PIO-C interface for the NetCDF function nc_put_att_ubyte.
 */
int PIOc_put_att_ubyte(int ncid, int varid, const char *name, nc_type xtype,
		       PIO_Offset len, const unsigned char *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

/** 
 * @ingroup PIOc_put_att_float
 * The PIO-C interface for the NetCDF function nc_put_att_float.
 */
int PIOc_put_att_float(int ncid, int varid, const char *name, nc_type xtype,
		       PIO_Offset len, const float *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

/** 
 * @ingroup PIOc_put_att_ulonglong
 * The PIO-C interface for the NetCDF function nc_put_att_ulonglong.
 */
int PIOc_put_att_ulonglong(int ncid, int varid, const char *name, nc_type xtype,
			   PIO_Offset len, const unsigned long long *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_put_att_ushort
 * The PIO-C interface for the NetCDF function nc_put_att_ushort.
 */
int PIOc_put_att_ushort(int ncid, int varid, const char *name, nc_type xtype,
			PIO_Offset len, const unsigned short *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_put_att_text
 * The PIO-C interface for the NetCDF function nc_put_att_text.
 */
int PIOc_put_att_text(int ncid, int varid, const char *name,
		      PIO_Offset len, const char *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

