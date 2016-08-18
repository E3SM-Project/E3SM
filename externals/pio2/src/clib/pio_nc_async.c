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
#include <pio.h>
#include <pio_internal.h>

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
int PIOc_inq(int ncid, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

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
            int msg = PIO_MSG_INQ; /* Message for async notification. */
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
	    LOG((2, "PIOc_inq ncid = %d ndims_present = %d nvars_present = %d ngatts_present = %d unlimdimid_present = %d",
		 ncid, ndims_present, nvars_present, ngatts_present, unlimdimid_present));
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
	{
	    LOG((2, "PIOc_inq calling ncmpi_inq unlimdimidp = %d", unlimdimidp));
            ierr = ncmpi_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp);
	    LOG((2, "PIOc_inq called ncmpi_inq"));
	    if (unlimdimidp)
		LOG((2, "PIOc_inq returned from ncmpi_inq unlimdimid = %d", *unlimdimidp));		
	}
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype == PIO_IOTYPE_NETCDF && file->do_io)
        {
	    LOG((2, "PIOc_inq calling classic nc_inq"));
            /* Should not be necessary to do this - nc_inq should
             * handle null pointers. This has been reported as a bug
             * to netCDF developers. */
            int tmp_ndims, tmp_nvars, tmp_ngatts, tmp_unlimdimid;
	    LOG((2, "PIOc_inq calling classic nc_inq"));
            ierr = nc_inq(ncid, &tmp_ndims, &tmp_nvars, &tmp_ngatts, &tmp_unlimdimid);
	    LOG((2, "PIOc_inq calling classic nc_inq"));
	    if (unlimdimidp)
		LOG((2, "classic tmp_unlimdimid = %d", tmp_unlimdimid));
            if (ndimsp)
                *ndimsp = tmp_ndims;
            if (nvarsp)
                *nvarsp = tmp_nvars;
            if (ngattsp)
                *ngattsp = tmp_ngatts;
            if (unlimdimidp)
                *unlimdimidp = tmp_unlimdimid;
	    if (unlimdimidp)
		LOG((2, "classic unlimdimid = %d", *unlimdimidp));
        }
	else if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
	{
	    LOG((2, "PIOc_inq calling netcdf-4 nc_inq"));	    
            ierr = nc_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp);
	}
#endif /* _NETCDF */
        LOG((2, "PIOc_inq netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);
    if (ierr)
	return check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (!ierr)
    {
        if (ndimsp)
            if ((mpierr = MPI_Bcast(ndimsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
        
        if (nvarsp)
            if ((mpierr = MPI_Bcast(nvarsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
        
        if (ngattsp)
            if ((mpierr = MPI_Bcast(ngattsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
        
        if (unlimdimidp)
            if ((mpierr = MPI_Bcast(unlimdimidp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    return ierr;
}

/** 
 * @ingroup PIOc_inq_ndims
 * The PIO-C interface for the NetCDF function nc_inq_ndims.
 */
int PIOc_inq_ndims (int ncid, int *ndimsp) 
{
    LOG((1, "PIOc_inq_ndims"));
    return PIOc_inq(ncid, ndimsp, NULL, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_nvars
 * The PIO-C interface for the NetCDF function nc_inq_nvars.
 */
int PIOc_inq_nvars(int ncid, int *nvarsp) 
{
    return PIOc_inq(ncid, NULL, nvarsp, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_natts
 * The PIO-C interface for the NetCDF function nc_inq_natts.
 */
int PIOc_inq_natts(int ncid, int *ngattsp) 
{
    return PIOc_inq(ncid, NULL, NULL, ngattsp, NULL);
}

/** 
 * @ingroup PIOc_inq_unlimdim
 * The PIO-C interface for the NetCDF function nc_inq_unlimdim.
 */
int PIOc_inq_unlimdim(int ncid, int *unlimdimidp) 
{
    LOG((1, "PIOc_inq_unlimdim ncid = %d unlimdimidp = %d", ncid, unlimdimidp));
    return PIOc_inq(ncid, NULL, NULL, NULL, unlimdimidp);
}

/** Internal function to provide inq_type function for pnetcdf. */
int pioc_pnetcdf_inq_type(int ncid, nc_type xtype, char *name,
                          PIO_Offset *sizep)
{
    int typelen;
    char typename[NC_MAX_NAME + 1];
    
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

    /* If pointers were supplied, copy results. */
    if (sizep)
        *sizep = typelen;
    if (name)
        strcpy(name, "some type");

    return PIO_NOERR;
}

/** 
 * @ingroup PIOc_typelen
 * The PIO-C interface for the NetCDF function nctypelen.
 */
int PIOc_inq_type(int ncid, nc_type xtype, char *name, PIO_Offset *sizep)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
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
            int msg = PIO_MSG_INQ_TYPE; /* Message for async notification. */
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
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);            
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = pioc_pnetcdf_inq_type(ncid, xtype, name, sizep);         
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_type(ncid, xtype, name, (size_t *)sizep);
#endif /* _NETCDF */
        LOG((2, "PIOc_inq_type netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);             
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
                return check_mpi(file, mpierr, __FILE__, __LINE__);
	    if (!mpierr)
		if ((mpierr = MPI_Bcast((void *)name, slen + 1, MPI_CHAR, ios->ioroot, ios->my_comm)))
		    return check_mpi(file, mpierr, __FILE__, __LINE__);
        }
        if (sizep)
            if ((mpierr = MPI_Bcast(sizep , 1, MPI_OFFSET, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);         
    }

    return ierr;
}

/** 
 * @ingroup PIOc_inq_format
 * The PIO-C interface for the NetCDF function nc_inq_format.
 */
int PIOc_inq_format (int ncid, int *formatp) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

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
            int msg = PIO_MSG_INQ_FORMAT;
            char format_present = formatp ? true : false;
        
            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&format_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_format(file->fh, formatp);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_format(file->fh, formatp);
#endif /* _NETCDF */
        LOG((2, "PIOc_inq netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);             
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (!ierr)
    {
        if (formatp)
            if ((mpierr = MPI_Bcast(formatp , 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
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
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

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
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);            
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_dim(file->fh, dimid, name, lenp);;     
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_dim(file->fh, dimid, name, (size_t *)lenp);;      
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);             
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
                return check_mpi(file, mpierr, __FILE__, __LINE__);
            if ((mpierr = MPI_Bcast((void *)name, slen + 1, MPI_CHAR, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
        }
        
        if (lenp)
            if ((mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);         
    }
    
    return ierr;
}

/** 
 * @ingroup PIOc_inq_dimname
 * The PIO-C interface for the NetCDF function nc_inq_dimname.
 */
int PIOc_inq_dimname(int ncid, int dimid, char *name) 
{
    return PIOc_inq_dim(ncid, dimid, name, NULL);
}

/** 
 * @ingroup PIOc_inq_dimlen
 * The PIO-C interface for the NetCDF function nc_inq_dimlen.
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
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

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
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);            
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* IO tasks call the netCDF functions. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_dimid(file->fh, name, idp);;
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_dimid(file->fh, name, idp);;
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);             
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results. */
    if (!ierr)
        if (idp)
            if ((mpierr = MPI_Bcast(idp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);

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
    int ndims;    /* The number of dimensions for this variable. */
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

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
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
        
    /* Call the netCDF layer. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
        {
            ierr = ncmpi_inq_varndims(file->fh, varid, &ndims);
            if (!ierr)
                ierr = ncmpi_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
        }
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
        {
            ierr = nc_inq_varndims(file->fh, varid, &ndims);
            if (!ierr)
                ierr = nc_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);
        }           
#endif /* _NETCDF */
    }

    if (ndimsp)
	LOG((2, "PIOc_inq_var ndims = %d ierr = %d", *ndimsp, ierr));

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);             
    if (ierr)
	return check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast the results for non-null pointers. */
    if (!ierr)
    {
        if (name)
        { 
            int slen;
            if(ios->iomaster)
                slen = strlen(name);
            if ((mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
            if ((mpierr = MPI_Bcast((void *)name, slen + 1, MPI_CHAR, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
        }
        if (xtypep)
            if ((mpierr = MPI_Bcast(xtypep, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
        
        if (ndimsp)
        {
	    if (ios->ioroot)
		LOG((2, "PIOc_inq_var about to Bcast ndims = %d ios->ioroot = %d", *ndimsp, ios->ioroot));
            if ((mpierr = MPI_Bcast(ndimsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);         
            file->varlist[varid].ndims = *ndimsp;
	    LOG((2, "PIOc_inq_var Bcast ndims = %d", *ndimsp));
        }
        if (dimidsp)
        {
            if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
            if ((mpierr = MPI_Bcast(dimidsp, ndims, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
        }
        if (nattsp)
            if ((mpierr = MPI_Bcast(nattsp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    return ierr;
}

/** 
 * @ingroup PIOc_inq_varname
 * The PIO-C interface for the NetCDF function nc_inq_varname.
 */
int PIOc_inq_varname (int ncid, int varid, char *name) 
{
    return PIOc_inq_var(ncid, varid, name, NULL, NULL, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_vartype
 * The PIO-C interface for the NetCDF function nc_inq_vartype.
 */
int PIOc_inq_vartype (int ncid, int varid, nc_type *xtypep) 
{
    return PIOc_inq_var(ncid, varid, NULL, xtypep, NULL, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_varndims
 * The PIO-C interface for the NetCDF function nc_inq_varndims.
 */
int PIOc_inq_varndims (int ncid, int varid, int *ndimsp) 
{
    return PIOc_inq_var(ncid, varid, NULL, NULL, ndimsp, NULL, NULL);
}

/** 
 * @ingroup PIOc_inq_vardimid
 * The PIO-C interface for the NetCDF function nc_inq_vardimid.
 */
int PIOc_inq_vardimid(int ncid, int varid, int *dimidsp) 
{
    return PIOc_inq_var(ncid, varid, NULL, NULL, NULL, dimidsp, NULL);
}

/** 
 * @ingroup PIOc_inq_varnatts
 * The PIO-C interface for the NetCDF function nc_inq_varnatts.
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
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    /* Caller must provide name. */
    if (!name || strlen(name) > NC_MAX_NAME)
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
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_varid(file->fh, name, varidp);;        
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_varid(file->fh, name, varidp);
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
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
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
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

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            char xtype_present = xtypep ? true : false;
            char len_present = lenp ? true : false;
            int namelen = strlen(name);
            
            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&xtype_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&len_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_att(file->fh, varid, name, xtypep, lenp);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_att(file->fh, varid, name, xtypep, (size_t *)lenp);
#endif /* _NETCDF */
        LOG((2, "PIOc_inq netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
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
 */
int PIOc_inq_attlen (int ncid, int varid, const char *name, PIO_Offset *lenp) 
{
    return PIOc_inq_att(ncid, varid, name, NULL, lenp);
}

/** 
 * @ingroup PIOc_inq_atttype
 * The PIO-C interface for the NetCDF function nc_inq_atttype.
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
int PIOc_inq_attname(int ncid, int varid, int attnum, char *name) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    LOG((1, "PIOc_inq_attname ncid = %d varid = %d attnum = %d", ncid, varid,
         attnum));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_INQ_ATTNAME;
            char name_present = name ? true : false;

            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
            
            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&attnum, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_attname(file->fh, varid, attnum, name);;
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_attname(file->fh, varid, attnum, name);;
#endif /* _NETCDF */
        LOG((2, "PIOc_inq_attname netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (!ierr)
        if (name)
        {
            int namelen = strlen(name);
            if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                check_mpi(file, mpierr, __FILE__, __LINE__);
            if ((mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->ioroot,
                                    ios->my_comm)))
                check_mpi(file, mpierr, __FILE__, __LINE__);
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
int PIOc_inq_attid(int ncid, int varid, const char *name, int *idp) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    /* User must provide name shorter than NC_MAX_NAME +1. */
    if (!name || strlen(name) > NC_MAX_NAME)
        return PIO_EINVAL;
    
    LOG((1, "PIOc_inq_attid ncid = %d varid = %d name = %s", ncid, varid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_INQ_ATTID;
            int namelen = strlen(name);
            char id_present = idp ? true : false;
            
            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&namelen, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((char *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&id_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_attid(file->fh, varid, name, idp);;
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_attid(file->fh, varid, name, idp);;
#endif /* _NETCDF */
        LOG((2, "PIOc_inq_attname netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results. */
    if (!ierr)
    {
        if (idp)
            if ((mpierr = MPI_Bcast(idp, 1, MPI_INT, ios->ioroot, ios->my_comm)))
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
int PIOc_rename_dim(int ncid, int dimid, const char *name) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    /* User must provide name of correct length. */
    if (!name || strlen(name) > NC_MAX_NAME)
        return PIO_EINVAL;

    LOG((1, "PIOc_rename_dim ncid = %d dimid = %d name = %s", ncid, dimid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_RENAME_DIM; /* Message for async notification. */
            int namelen = strlen(name);

            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&dimid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            LOG((2, "PIOc_rename_dim Bcast file->fh = %d dimid = %d namelen = %d name = %s",
                 file->fh, dimid, namelen, name));
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_rename_dim(file->fh, dimid, name);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_rename_dim(file->fh, dimid, name);;
#endif /* _NETCDF */
        LOG((2, "PIOc_inq netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
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
int PIOc_rename_var(int ncid, int varid, const char *name) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    /* User must provide name of correct length. */
    if (!name || strlen(name) > NC_MAX_NAME)
        return PIO_EINVAL;

    LOG((1, "PIOc_rename_var ncid = %d varid = %d name = %s", ncid, varid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_RENAME_VAR; /* Message for async notification. */
            int namelen = strlen(name);

            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            LOG((2, "PIOc_rename_var Bcast file->fh = %d varid = %d namelen = %d name = %s",
                 file->fh, varid, namelen, name));
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_rename_var(file->fh, varid, name);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_rename_var(file->fh, varid, name);;
#endif /* _NETCDF */
        LOG((2, "PIOc_inq netcdf call returned %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
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
 * @return PIO_NOERR for success, error code otherwise.  See
 * PIOc_Set_File_Error_Handling
 */
int PIOc_rename_att (int ncid, int varid, const char *name,
                     const char *newname) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI functions. */

    /* User must provide names of correct length. */
    if (!name || strlen(name) > NC_MAX_NAME ||
        !newname || strlen(newname) > NC_MAX_NAME)
        return PIO_EINVAL;

    LOG((1, "PIOc_rename_att ncid = %d varid = %d name = %s newname = %s",
         ncid, varid, name, newname));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_RENAME_ATT; /* Message for async notification. */
            int namelen = strlen(name);
            int newnamelen = strlen(newname);

            if (ios->compmaster) 
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((char *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&newnamelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((char *)newname, newnamelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_rename_att(file->fh, varid, name, newname);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_rename_att(file->fh, varid, name, newname);
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);

    LOG((2, "PIOc_rename_att succeeded"));        
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
int PIOc_del_att(int ncid, int varid, const char *name) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI functions. */

    /* User must provide name of correct length. */
    if (!name || strlen(name) > NC_MAX_NAME)
        return PIO_EINVAL;

    LOG((1, "PIOc_del_att ncid = %d varid = %d name = %s", ncid, varid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_DEL_ATT;
            int namelen = strlen(name); /* Length of name string. */

            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&varid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((char *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);
	if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_del_att(file->fh, varid, name);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_del_att(file->fh, varid, name);
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);

    LOG((2, "PIOc_del_att succeeded"));        
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
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI functions. */

    LOG((1, "PIOc_set_fill ncid = %d fillmode = %d old_modep = %d", ncid, fillmode,
         old_modep));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_SET_FILL;

            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
            
            if (!mpierr)
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }


    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_set_fill(file->fh, fillmode, old_modep);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_set_fill(file->fh, fillmode, old_modep);
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);

    LOG((2, "PIOc_set_fill succeeded"));        
    return ierr;
}

/** This is an internal function that handles both PIOc_enddef and
 * PIOc_redef. 
 * @param ncid the ncid of the file to enddef or redef
 * @param is_enddef set to non-zero for enddef, 0 for redef. 
 * @returns PIO_NOERR on success, error code on failure. */
int pioc_change_def(int ncid, int is_enddef)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI functions. */

    LOG((1, "pioc_change_def ncid = %d is_enddef = %d", ncid, is_enddef));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;
    LOG((2, "pioc_change_def found file"));
    
    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = is_enddef ? PIO_MSG_ENDDEF : PIO_MSG_REDEF;
	    LOG((2, "sending message msg = %d", msg));
            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)            
                mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	    LOG((2, "pioc_change_def ncid = %d mpierr = %d", file->fh, mpierr));
        }

        /* Handle MPI errors. */
	LOG((2, "pioc_change_def handling MPI errors"));	
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    LOG((2, "pioc_change_def ios->ioproc = %d", ios->ioproc));    
    if (ios->ioproc)
    {
	LOG((2, "pioc_change_def calling netcdf function"));
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            if (is_enddef)
                ierr = ncmpi_enddef(file->fh);
            else
                ierr = ncmpi_redef(file->fh);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            if (is_enddef)          
                ierr = nc_enddef(file->fh);
            else
                ierr = nc_redef(file->fh);
#endif /* _NETCDF */
	LOG((2, "pioc_change_def ierr = %d", ierr));
    }

    /* Broadcast and check the return code. */
    LOG((2, "pioc_change_def bcasting return code ierr = %d", ierr));
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);            
    if (ierr)
	return check_netcdf(file, ierr, __FILE__, __LINE__);
    LOG((2, "pioc_change_def succeeded"));

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
    return pioc_change_def(ncid, 1);
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
int PIOc_redef(int ncid) 
{
    return pioc_change_def(ncid, 0);    
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
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    /* User must provide name. */
    if (!name || strlen(name) > NC_MAX_NAME)
        return PIO_EINVAL;

    LOG((1, "PIOc_def_dim ncid = %d name = %s len = %d", ncid, name, len));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_DEF_DIM;
            int namelen = strlen(name);

            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);

            if (!mpierr)
                mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&len, 1, MPI_INT,  ios->compmaster, ios->intercomm);
        }


        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_def_dim(file->fh, name, len, idp);;        
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_def_dim(file->fh, name, (size_t)len, idp);            
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr2, __FILE__, __LINE__);           
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if (!ierr)
        if (idp)
            if ((mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm)))
                check_mpi(file, mpierr, __FILE__, __LINE__);            

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
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    /* User must provide name and storage for varid. */
    if (!name || !varidp || strlen(name) > NC_MAX_NAME)
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
            int namelen = strlen(name);
            
            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);
            
            if (!mpierr)
                mpierr = MPI_Bcast(&(ncid), 1, MPI_INT, ios->compmaster, ios->intercomm);
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
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_def_var(ncid, name, xtype, ndims, dimidsp, varidp);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_def_var(ncid, name, xtype, ndims, dimidsp, varidp);
#ifdef _NETCDF4
        /* For netCDF-4 serial files, turn on compression for this variable. */
        if (!ierr && file->iotype == PIO_IOTYPE_NETCDF4C)
            ierr = nc_def_var_deflate(ncid, *varidp, 0, 1, 1);

        /* For netCDF-4 parallel files, set parallel access to collective. */
        if (!ierr && file->iotype == PIO_IOTYPE_NETCDF4P)
            ierr = nc_var_par_access(ncid, *varidp, NC_COLLECTIVE);
#endif /* _NETCDF4 */
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results. */
    if (!ierr)
        if (varidp)
            if ((mpierr = MPI_Bcast(varidp , 1, MPI_INT, ios->ioroot, ios->my_comm)))
                check_mpi(file, mpierr, __FILE__, __LINE__);
    
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
int PIOc_inq_var_fill(int ncid, int varid, int *no_fill, void *fill_valuep) 
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

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
            int msg = PIO_MSG_INQ_VAR_FILL;

            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_inq_var_fill(file->fh, varid, no_fill, fill_valuep);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_inq_var_fill(file->fh, varid, no_fill, fill_valuep);
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);

    /* Broadcast results to all tasks. Ignore NULL parameters. */    
    if (!ierr)
        if (fill_valuep)
            if ((mpierr = MPI_Bcast(fill_valuep, 1, MPI_INT, ios->ioroot, ios->my_comm)))
                check_mpi(file, mpierr, __FILE__, __LINE__);
    
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
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
    PIO_Offset attlen, typelen;
    nc_type atttype;

    /* User must provide a name and destination pointer. */
    if (!name || !ip || strlen(name) > NC_MAX_NAME)
        return PIO_EINVAL;

    LOG((1, "PIOc_get_att ncid %d varid %d name %s", ncid, varid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* Run these on all tasks if async is not in use, but only on
     * non-IO tasks if async is in use. */
    if (!ios->async_interface || !ios->ioproc)
    {
        /* Get the type and length of the attribute. */
        if ((ierr = PIOc_inq_att(file->fh, varid, name, &atttype, &attlen)))
        {
            check_netcdf(file, ierr, __FILE__, __LINE__);
            return ierr;
        }

        /* Get the length (in bytes) of the type. */
        if ((ierr = PIOc_inq_type(file->fh, atttype, NULL, &typelen)))
        {
            check_netcdf(file, ierr, __FILE__, __LINE__);
            return ierr;
        }
    }

    /* If async is in use, and this is not an IO task, bcast the
     * parameters and the attribute and type information we fetched. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_GET_ATT;

            /* Send the message to IO master. */
            if(ios->compmaster) 
                mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

            /* Send the function parameters. */
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
            if (!mpierr)
                mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);
        
        /* Broadcast values currently only known on computation tasks to IO tasks. */
        LOG((2, "PIOc_get_att bcast from comproot = %d attlen = %d typelen = %d", ios->comproot, attlen, typelen));
        if ((mpierr = MPI_Bcast(&attlen, 1, MPI_OFFSET, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr, __FILE__, __LINE__);
        LOG((2, "PIOc_get_att bcast complete attlen = %d typelen = %d", attlen, typelen));
    }
        
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_get_att(file->fh, varid, name, ip);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_get_att(file->fh, varid, name, ip);
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    /* Broadcast results to all tasks. */
    if (!ierr)
    {
        if ((mpierr = MPI_Bcast(ip, (int)attlen * typelen, MPI_BYTE, ios->ioroot,
                                ios->my_comm)))
        {
            check_mpi(file, mpierr, __FILE__, __LINE__);
            return PIO_EIO;
        }
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
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    PIO_Offset typelen; /* Length (in bytes) of the type. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    LOG((1, "PIOc_put_att ncid = %d varid = %d name = %s", ncid, varid, name));

    /* Find the info about this file. */
    if (!(file = pio_get_file_from_id(ncid)))
        return PIO_EBADID;
    ios = file->iosystem;

    /* Run these on all tasks if async is not in use, but only on
     * non-IO tasks if async is in use. */
    if (!ios->async_interface || !ios->ioproc)
    {
        /* Get the length (in bytes) of the type. */
        if ((ierr = PIOc_inq_type(ncid, xtype, NULL, &typelen)))
        {
            check_netcdf(file, ierr, __FILE__, __LINE__);
            return ierr;
        }
        LOG((2, "PIOc_put_att typelen = %d", ncid, typelen));
    }
    
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
                mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET,  ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((void *)op, len * typelen, MPI_BYTE, ios->compmaster,
                                   ios->intercomm);
            LOG((2, "PIOc_put_att finished bcast ncid = %d varid = %d namelen = %d name = %s "
                 "len = %d typelen = %d", file->fh, varid, namelen, name, len, typelen));
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);           
        if (mpierr)
	    return check_mpi(file, mpierr, __FILE__, __LINE__);

        /* Broadcast values currently only known on computation tasks to IO tasks. */
        LOG((2, "PIOc_put_att bcast from comproot = %d typelen = %d", ios->comproot, typelen));
        if ((mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr, __FILE__, __LINE__);
    }
    
    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
            ierr = ncmpi_put_att(file->fh, varid, name, xtype, len, op);
#endif /* _PNETCDF */
#ifdef _NETCDF
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
            ierr = nc_put_att(file->fh, varid, name, xtype, (size_t)len, op);
#endif /* _NETCDF */
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
    {
        check_mpi(file, mpierr, __FILE__, __LINE__);            
        return PIO_EIO;
    }
    check_netcdf(file, ierr, __FILE__, __LINE__);
    
    return ierr;
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
 * @ingroup PIOc_get_att_uchar
 * The PIO-C interface for the NetCDF function nc_get_att_uchar.
 */
int PIOc_get_att_uchar (int ncid, int varid, const char *name, unsigned char *ip) 
{
    return PIOc_get_att(ncid, varid, name, (void *)ip);
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
 * @ingroup PIOc_get_att_float
 * The PIO-C interface for the NetCDF function nc_get_att_float.
 */
int PIOc_get_att_float (int ncid, int varid, const char *name, float *ip) 
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
 * @ingroup PIOc_put_att_long
 * The PIO-C interface for the NetCDF function nc_put_att_long.
 */
int PIOc_put_att_long(int ncid, int varid, const char *name, nc_type xtype,
                      PIO_Offset len, const long *op) 
{
    return PIOc_put_att(ncid, varid, name, NC_CHAR, len, op);
}

/** 
 * @ingroup PIOc_put_att_int
 * The PIO-C interface for the NetCDF function nc_put_att_int.
 */
int PIOc_put_att_int(int ncid, int varid, const char *name, nc_type xtype,
                     PIO_Offset len, const int *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);    
}

/** 
 * @ingroup PIOc_put_att_uchar
 * The PIO-C interface for the NetCDF function nc_put_att_uchar.
 */
int PIOc_put_att_uchar(int ncid, int varid, const char *name, nc_type xtype,
                       PIO_Offset len, const unsigned char *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_put_att_longlong
 * The PIO-C interface for the NetCDF function nc_put_att_longlong.
 */
int PIOc_put_att_longlong(int ncid, int varid, const char *name, nc_type xtype,
                          PIO_Offset len, const long long *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_put_att_uint
 * The PIO-C interface for the NetCDF function nc_put_att_uint.
 */
int PIOc_put_att_uint(int ncid, int varid, const char *name, nc_type xtype,
                      PIO_Offset len, const unsigned int *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_put_att_ubyte
 * The PIO-C interface for the NetCDF function nc_put_att_ubyte.
 */
int PIOc_put_att_ubyte(int ncid, int varid, const char *name, nc_type xtype,
                       PIO_Offset len, const unsigned char *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}

/** 
 * @ingroup PIOc_put_att_float
 * The PIO-C interface for the NetCDF function nc_put_att_float.
 */
int PIOc_put_att_float(int ncid, int varid, const char *name, nc_type xtype,
                       PIO_Offset len, const float *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
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
 * @ingroup PIOc_put_att_double
 * The PIO-C interface for the NetCDF function nc_put_att_double.
 */
int PIOc_put_att_double(int ncid, int varid, const char *name, nc_type xtype,
                        PIO_Offset len, const double *op) 
{
    return PIOc_put_att(ncid, varid, name, xtype, len, op);
}


