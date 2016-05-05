/**
 * @file 
 * @author Ed Hartnett
 * @date  2016
 * @brief PIO async msg handling
 *
 * @see http://code.google.com/p/parallelio/
 */

#include <pio.h>
#include <pio_internal.h>

/** This function is run on the IO tasks to create a netCDF file. */
int create_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len;
    int iotype;
    char *filename;
    int mode;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d create_file_handler comproot = %d\n", my_rank, ios->comproot);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d create_file_handler got parameter len = %d\n", my_rank, len);
    if (!(filename = malloc(len + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&mode, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    printf("%d create_file_handler got parameters len = %d "
    	   "filename = %s iotype = %d mode = %d\n",
    	   my_rank, len, filename, iotype, mode);

    /* Call the create file function. */
    if ((ret = PIOc_createfile(ios->iosysid, &ncid, &iotype, filename, mode)))
	return ret;
    
    /* Free resources. */
    free(filename);
    
    printf("%d create_file_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to close a netCDF file. It is
 * only ever run on the IO tasks. 
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int close_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d close_file_handler\n", my_rank);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d create_file_handler got parameter ncid = %d\n", my_rank, ncid);

    /* Call the close file function. */
    if ((ret = PIOc_closefile(ncid)))
	return ret;
    
    printf("%d close_file_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to inq a netCDF file. It is
 * only ever run on the IO tasks. 
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @param msg the message sent my the comp root task. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int inq_handler(iosystem_desc_t *ios, int msg)
{
    int ncid;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d inq_handler msg = %d\n", my_rank, msg);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;

    /* Call the inq file function. */
    int ndims, nvars, ngatts, unlimdimid;
    int *ndimsp = NULL, *nvarsp = NULL, *ngattsp = NULL, *unlimdimidp = NULL;
    switch (msg)
    {
    case PIO_MSG_INQ:
	ndimsp = &ndims;
	nvarsp = &nvars;
	ngattsp = &ngatts;
	unlimdimidp = &unlimdimid;
	break;
    case PIO_MSG_INQ_NVARS:
	nvarsp = &nvars;
	break;
    case PIO_MSG_INQ_NDIMS:
	ndimsp = &ndims;
	break;
    case PIO_MSG_INQ_NATTS:
	ngattsp = &ngatts;
	break;
    case PIO_MSG_INQ_UNLIMDIM:
	unlimdimidp = &unlimdimid;
	break;
    default:
	return PIO_EINVAL;
    }

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp)))
	return ret;
    
    return PIO_NOERR;
}

/** Do an inq_dim on a netCDF dimension. This function is only run on
 * IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @param msg the message sent my the comp root task. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int inq_dim_handler(iosystem_desc_t *ios, int msg)
{
    int ncid;
    int dimid;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d inq_handler msg = %d\n", my_rank, msg);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&dimid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;

    /* Call the inq_dim function. */
    char *dimnamep = NULL;
    PIO_Offset *dimlenp = NULL;
    char dimname[NC_MAX_NAME + 1];
    PIO_Offset dimlen;
    switch (msg)
    {
    case PIO_MSG_INQ_DIM:
	dimnamep = dimname;
	dimlenp = &dimlen;
	break;
    case PIO_MSG_INQ_DIMLEN:
	dimlenp = &dimlen;
	break;
    case PIO_MSG_INQ_DIMNAME:
	dimnamep = dimname;
	break;
    default:
	return PIO_EINVAL;
    }

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq_dim(ncid, dimid, dimnamep, dimlenp)))
	return ret;
    
    return PIO_NOERR;
}

/** Do an inq_dimid on a netCDF dimension name. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int inq_dimid_handler(iosystem_desc_t *ios)
{
    int ncid;
    int dimid;
    int mpierr;
    int ret;
    int namelen;
    char *name;
    
    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if (!(name = malloc((namelen + 1) * sizeof(char))))
	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;

    /* Call the inq_dimid function. */
    if ((ret = PIOc_inq_dimid(ncid, name, &dimid)))
	return ret;

    /* Free resources. */
    free(name);
    
    return PIO_NOERR;
}

/** Handle attribute operations. This code only runs on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @param msg the message sent my the comp root task. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int att_handler(iosystem_desc_t *ios, int msg)
{
    int ncid;
    int varid;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d inv_var_handler ncid = %d varid = %d\n", my_rank, ncid, varid);

    /* Call the inq_var function. */
    char name[NC_MAX_NAME + 1], *namep;
    nc_type xtype, *xtypep = NULL;
    int *ndimsp = NULL, *dimidsp = NULL, *nattsp = NULL;    
    int ndims, dimids[NC_MAX_DIMS], natts;    
    switch (msg)
    {
    case PIO_MSG_INQ_VAR:
	namep = name;
	xtypep = &xtype;
	ndimsp = &ndims;
	dimidsp = dimids;
	nattsp = &natts;
	break;
    case PIO_MSG_INQ_VARNATTS:
	nattsp = &natts;
	break;
    case PIO_MSG_INQ_VARNAME:
	namep = name;
	break;
    case PIO_MSG_INQ_VARNDIMS:
	ndimsp = &ndims;
	break;
    case PIO_MSG_INQ_VARDIMID:
	dimidsp = dimids;
	break;
    case PIO_MSG_INQ_VARTYPE:
	xtypep = &xtype;
	break;
    default:
	return PIO_EINVAL;
    }

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq_var(ncid, varid, namep, xtypep, ndimsp, dimidsp, nattsp)))
	return ret;
    
    return PIO_NOERR;
}

/** Do an inq_var on a netCDF variable. This function is only run on
 * IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @param msg the message sent my the comp root task. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int inq_var_handler(iosystem_desc_t *ios, int msg)
{
    int ncid;
    int varid;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d inv_var_handler ncid = %d varid = %d\n", my_rank, ncid, varid);

    /* Call the inq_var function. */
    char name[NC_MAX_NAME + 1], *namep;
    nc_type xtype, *xtypep = NULL;
    int *ndimsp = NULL, *dimidsp = NULL, *nattsp = NULL;    
    int ndims, dimids[NC_MAX_DIMS], natts;    
    switch (msg)
    {
    case PIO_MSG_INQ_VAR:
	namep = name;
	xtypep = &xtype;
	ndimsp = &ndims;
	dimidsp = dimids;
	nattsp = &natts;
	break;
    case PIO_MSG_INQ_VARNATTS:
	nattsp = &natts;
	break;
    case PIO_MSG_INQ_VARNAME:
	namep = name;
	break;
    case PIO_MSG_INQ_VARNDIMS:
	ndimsp = &ndims;
	break;
    case PIO_MSG_INQ_VARDIMID:
	dimidsp = dimids;
	break;
    case PIO_MSG_INQ_VARTYPE:
	xtypep = &xtype;
	break;
    default:
	return PIO_EINVAL;
    }

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq_var(ncid, varid, namep, xtypep, ndimsp, dimidsp, nattsp)))
	return ret;
    
    return PIO_NOERR;
}

/** Do an inq_varid on a netCDF variable name. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int inq_varid_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    int ret;
    int namelen;
    char *name;
    
    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if (!(name = malloc((namelen + 1) * sizeof(char))))
	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;

    /* Call the inq_dimid function. */
    if ((ret = PIOc_inq_varid(ncid, name, &varid)))
	return ret;

    /* Free resources. */
    free(name);
    
    return PIO_NOERR;
}

/** This function is run on the IO tasks to sync a netCDF file. 
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int sync_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d sync_file_handler\n", my_rank);

    /* Get the parameters for this function that the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d sync_file_handler got parameter ncid = %d\n", my_rank, ncid);

    /* Call the sync file function. */
    if ((ret = PIOc_sync(ncid)))
	return ret;
    
    printf("%d sync_file_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to enddef a netCDF file. 
 *
 * @param ios pointer to the iosystem_desc_t. 
 * @return PIO_NOERR for success, error code otherwise. 
*/
int enddef_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d enddef_file_handler\n", my_rank);

    /* Get the parameters for this function that the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d enddef_file_handler got parameter ncid = %d\n", my_rank, ncid);

    /* Call the sync file function. */
    if ((ret = PIOc_enddef(ncid)))
	return ret;
    
    printf("%d enddef_file_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to define a netCDF
 *  variable. */
int def_var_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len, namelen;
    int iotype;
    char *name;
    int mode;
    int mpierr;
    int ret;
    int varid;
    nc_type xtype;
    int ndims;
    int *dimids;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d def_var_handler comproot = %d\n", my_rank, ios->comproot);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if (!(name = malloc(namelen + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&xtype, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if (!(dimids = malloc(ndims * sizeof(int))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast(dimids, ndims, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d def_var_handler got parameters namelen = %d "
    	   "name = %s len = %d ncid = %d\n",
    	   my_rank, namelen, name, len, ncid);

    /* Call the create file function. */
    if ((ret = PIOc_def_var(ncid, name, xtype, ndims, dimids, &varid)))
	return ret;
    
    /* Free resources. */
    free(name);
    free(dimids);
    
    printf("%d def_var_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to define a netCDF
 * dimension. */
int def_dim_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len, namelen;
    int iotype;
    char *name;
    int mode;
    int mpierr;
    int ret;
    int dimid;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d def_dim_handler comproot = %d\n", my_rank, ios->comproot);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d def_dim_handler ncid = %d namelen %d\n", my_rank, ncid, namelen);    
    if (!(name = malloc(namelen + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)name, namelen + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d def_dim_handler got parameters namelen = %d "
    	   "name = %s len = %d ncid = %d\n",
    	   my_rank, namelen, name, len, ncid);

    /* Call the create file function. */
    if ((ret = PIOc_def_dim(ncid, name, len, &dimid)))
	return ret;
    
    /* Free resources. */
    free(name);
    
    printf("%d def_dim_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks. It reads or writes an array
 *  of data to a netCDF variable. 
 *
 * @param ios pointer to the iosystem_desc_t data.
 * @param msg the message sent my the comp root task. 
 *
 * @return PIO_NOERR for success, error code otherwise. */
int vara_handler(iosystem_desc_t *ios, int msg)
{
    int ncid;
    int len, namelen;
    int iotype;
    char *name;
    int mode;
    int mpierr;
    int ret;
    int varid;
    nc_type xtype;
    int ndims;
    int *dimids;
    void *data;
    int size_in_bytes;
    PIO_Offset *count, *start;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d def_var_handler comproot = %d\n", my_rank, ios->comproot);

    if (msg == PIO_MSG_PUT_VARA)
    {
	/* Get the parameters for this function that the he comp master
	 * task is broadcasting. */
	if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	if ((ret = PIOc_inq_varndims(ncid, varid, &ndims)))
	    return ret;
	if (!(start = malloc(ndims * sizeof(int))))
	    return PIO_ENOMEM;
	if (!(count = malloc(ndims * sizeof(int))))
	    return PIO_ENOMEM;
	if ((mpierr = MPI_Bcast(start, ndims, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	if ((mpierr = MPI_Bcast(count, ndims, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	int size = 1;
	for (int d = 0; d < ndims; d++)
	    size *= count[d];
	size_in_bytes = size * sizeof(int);
	if (!(data = malloc(size_in_bytes)))
	    return PIO_ENOMEM;	
	if ((mpierr = MPI_Bcast(data, size_in_bytes, MPI_INT, 0,
				ios->intercomm)))
	    return PIO_EIO;
	printf("%d def_var_handler got parameters namelen = %d "
	       "name = %s len = %d ncid = %d\n",
	       my_rank, namelen, name, len, ncid);

	/* Call the create file function. */
	if ((ret = PIOc_put_vara_int(ncid, varid, start, count, data)))
	    return ret;
    
	/* Free resources. */
	free(start);
	free(count);
	free(data);
    }
    else
    {
	
    }

    printf("%d vara_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to open a netCDF file. 
 *
 * @param ios pointer to the iosystem_desc_t data.
 *
 * @return PIO_NOERR for success, error code otherwise. */
int open_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len;
    int iotype;
    char *filename;
    int mode;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d open_file_handler comproot = %d\n", my_rank, ios->comproot);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d open_file_handler got parameter len = %d\n", my_rank, len);
    if (!(filename = malloc(len + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&mode, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    printf("%d open_file_handler got parameters len = %d "
    	   "filename = %s iotype = %d mode = %d\n",
    	   my_rank, len, filename, iotype, mode);

    /* Call the open file function. */
    if ((ret = PIOc_openfile(ios->iosysid, &ncid, &iotype, filename, mode)))
	return ret;
    
    /* Free resources. */
    free(filename);
    
    printf("%d open_file_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to delete a netCDF file. 
 *
 * @param ios pointer to the iosystem_desc_t data.
 *
 * @return PIO_NOERR for success, error code otherwise. */
int delete_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len;
    char *filename;
    int mpierr;
    int ret;
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    printf("%d delete_file_handler comproot = %d\n", my_rank, ios->comproot);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d open_file_handler got parameter len = %d\n", my_rank, len);
    if (!(filename = malloc(len + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    printf("%d delete_file_handler got parameters len = %d filename = %s\n",
	   my_rank, len, filename);

    /* Call the delete file function. */
    if ((ret = PIOc_deletefile(ios->iosysid, filename)))
	return ret;
    
    /* Free resources. */
    free(filename);
    
    printf("%d delete_file_handler succeeded!\n", my_rank);
    return PIO_NOERR;
}

int initdecomp_dof_handler(iosystem_desc_t *ios)
{
    return PIO_NOERR;
}

int writedarray_handler(iosystem_desc_t *ios)
{
    return PIO_NOERR;
}

int readdarray_handler(iosystem_desc_t *ios)
{
    return PIO_NOERR;
}

int seterrorhandling_handler(iosystem_desc_t *ios)
{
    return PIO_NOERR;
}

int var_handler(iosystem_desc_t *ios, int msg)
{
    return PIO_NOERR;
}

int freedecomp_handler(iosystem_desc_t *ios)
{
    return PIO_NOERR;
}

int finalize_handler(iosystem_desc_t *ios)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    printf("%d finalize_handler called\n", my_rank);
    return PIO_NOERR;
}

int pio_callback_handler(iosystem_desc_t *ios, int msg)
{
    return PIO_NOERR;
}

/** This function is called by the IO tasks.  This function will not
 return, unless there is an error. */
int pio_msg_handler(int io_rank, int component_count, iosystem_desc_t *iosys)
{
    iosystem_desc_t *my_iosys;    
    int msg = 0;
    MPI_Request req[component_count];
    MPI_Status status;
    int index;
    int mpierr;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    printf("%d pio_msg_handler called\n", my_rank);
    
    /* Have IO comm rank 0 (the ioroot) register to receive
     * (non-blocking) for a message from each of the comproots. */
    if (!io_rank)
    {
	for (int cmp = 0; cmp < component_count; cmp++)
	{
	    my_iosys = &iosys[cmp];
	    printf("%d about to call MPI_Irecv\n", my_rank);	    
	    mpierr = MPI_Irecv(&msg, 1, MPI_INT, my_iosys->comproot, MPI_ANY_TAG,
			       my_iosys->union_comm, &req[cmp]);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);
	}
    }

    /* If the message is not -1, keep processing messages. */
    while (msg != -1)
    {
	/* Wait until any one of the requests are complete. */
	if (!io_rank)
	{
	    printf("%d about to call MPI_Waitany req[0] = %d MPI_REQUEST_NULL = %d\n",
		   my_rank, req[0], MPI_REQUEST_NULL);	    
	    mpierr = MPI_Waitany(component_count, req, &index, &status);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    printf("%d Waitany returned index = %d req[%d] = %d\n", my_rank,
		   index, index, req[index]);
	}

	/* Broadcast the index of the computational component that
	 * originated the request to the rest of the IO tasks. */
	printf("%d about to call index MPI_Bcast index = %d\n", my_rank, index);	    
	mpierr = MPI_Bcast(&index, 1, MPI_INT, 0, iosys->io_comm);
	CheckMPIReturn(mpierr, __FILE__, __LINE__);
	my_iosys = &iosys[index];
	printf("%d index MPI_Bcast complete index = %d\n", my_rank, index);	    	

	/* Broadcast the msg value to the rest of the IO tasks. */
	printf("%d about to call msg MPI_Bcast\n", my_rank);	    
	mpierr = MPI_Bcast(&msg, 1, MPI_INT, 0, my_iosys->io_comm);
	CheckMPIReturn(mpierr, __FILE__, __LINE__);
	printf("%d msg MPI_Bcast complete msg = %d\n", my_rank, msg);	    	

	/* Handle the message. This code is run on all IO tasks. */
	switch (msg)
	{
	case PIO_MSG_CREATE_FILE:
	    create_file_handler(my_iosys);
	    break;
	case PIO_MSG_SYNC:
	    sync_file_handler(my_iosys);
	    break;
	case PIO_MSG_ENDDEF:
	    enddef_file_handler(my_iosys);
	    break;
	case PIO_MSG_OPEN_FILE:
	    open_file_handler(my_iosys);
	    break;
	case PIO_MSG_CLOSE_FILE:
	    close_file_handler(my_iosys);
	    break;
	case PIO_MSG_DELETE_FILE:
	    delete_file_handler(my_iosys);
	    break;
	case PIO_MSG_DEF_DIM:
	    def_dim_handler(my_iosys);
	    break;
	case PIO_MSG_DEF_VAR:
	    def_var_handler(my_iosys);
	    break;
	case PIO_MSG_INQ:
	case PIO_MSG_INQ_NVARS:
	case PIO_MSG_INQ_NDIMS:
	case PIO_MSG_INQ_NATTS:
	case PIO_MSG_INQ_UNLIMDIM:
	    inq_handler(my_iosys, msg);
	    break;
	case PIO_MSG_INQ_DIM:
	case PIO_MSG_INQ_DIMLEN:
	case PIO_MSG_INQ_DIMNAME:
	    inq_dim_handler(my_iosys, msg);
	    break;
	case PIO_MSG_INQ_DIMID:
	    inq_dimid_handler(my_iosys);
	    break;
	case PIO_MSG_INQ_VAR:
	case PIO_MSG_INQ_VARNATTS:
	case PIO_MSG_INQ_VARNAME:
	case PIO_MSG_INQ_VARNDIMS:
	case PIO_MSG_INQ_VARDIMID:
	case PIO_MSG_INQ_VARTYPE:
	    inq_var_handler(my_iosys, msg);
	    break;
	case PIO_MSG_GET_ATT_INT:
	    att_handler(my_iosys, msg);
	    break;
	case PIO_MSG_PUT_ATT_INT:
	    att_handler(my_iosys, msg);
	    break;
	case PIO_MSG_INQ_VARID:
	    inq_varid_handler(my_iosys);
	    break;
	case PIO_MSG_INITDECOMP_DOF:
	    initdecomp_dof_handler(my_iosys);
	    break;
	case PIO_MSG_WRITEDARRAY:
	    writedarray_handler(my_iosys);
	    break;
	case PIO_MSG_READDARRAY:
	    readdarray_handler(my_iosys);
	    break;
	case PIO_MSG_SETERRORHANDLING:
	    seterrorhandling_handler(my_iosys);
	    break;
	case PIO_MSG_GET_VAR1:
	    var_handler(my_iosys, msg);
	    break;
	case PIO_MSG_PUT_VAR1:
	    var_handler(my_iosys, msg);
	    break;
	case PIO_MSG_PUT_VARA:
	    vara_handler(my_iosys, msg);
	    break;
	case PIO_MSG_FREEDECOMP:
	    freedecomp_handler(my_iosys);
	    break;
	case PIO_MSG_EXIT:
	    finalize_handler(my_iosys);
	    msg = -1;
	    break;
	default:
	    pio_callback_handler(my_iosys, msg);
	}

	/* Unless finalize was called, listen for another msg from the
	 * component whose message we just handled. */
	if (!io_rank && msg != -1)
	{
	    my_iosys = &iosys[index];
	    mpierr = MPI_Irecv(&msg, 1, MPI_INT, my_iosys->comproot, MPI_ANY_TAG, my_iosys->union_comm,
			       &req[index]);
	    printf("%d pio_msg_handler called MPI_Irecv req[%d] = %d\n", my_rank, index, req[index]);	    
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);
	}
	
    }

    return PIO_NOERR;
}

int
pio_iosys_print(int my_rank, iosystem_desc_t *iosys)
{
    printf("%d iosysid: %d\n", my_rank, iosys->iosysid);
    if (iosys->union_comm == MPI_COMM_NULL)
	printf("%d union_comm: MPI_COMM_NULL ", my_rank);
    else
	printf("%d union_comm: %d ", my_rank, iosys->union_comm);
    
    if (iosys->comp_comm == MPI_COMM_NULL)
	printf("comp_comm: MPI_COMM_NULL ");
    else
	printf("comp_comm: %d ", iosys->comp_comm);
    
    if (iosys->io_comm == MPI_COMM_NULL)
	printf("io_comm: MPI_COMM_NULL ");
    else
	printf("io_comm: %d ", iosys->io_comm);
    
    if (iosys->intercomm == MPI_COMM_NULL)
	printf("intercomm: MPI_COMM_NULL\n");
    else
	printf("intercomm: %d\n", iosys->intercomm);
    
    printf("%d num_iotasks=%d num_comptasks=%d union_rank=%d, comp_rank=%d, "
	   "io_rank=%d async_interface=%d\n",
	   my_rank, iosys->num_iotasks, iosys->num_comptasks, iosys->union_rank,
	   iosys->comp_rank, iosys->io_rank, iosys->async_interface);

    printf("%d ioroot=%d comproot=%d iomaster=%d, compmaster=%d\n",
	   my_rank, iosys->ioroot, iosys->comproot, iosys->iomaster,
	   iosys->compmaster);

    printf("%d iotasks:", my_rank);
    for (int i = 0; i < iosys->num_iotasks; i++)
	printf("%d ", iosys->ioranks[i]);
    printf("\n");
    return PIO_NOERR;
}

