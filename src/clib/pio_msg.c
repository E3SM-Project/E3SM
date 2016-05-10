/**
 * @file
 * @author Ed Hartnett
 * @date  2016
 * @brief PIO async msg handling
 *
 * @see http://code.google.com/p/parallelio/
 */

#include <config.h>
#include <pio.h>
#include <pio_internal.h>

#ifdef PIO_ENABLE_LOGGING
extern int my_rank;
extern int pio_log_level;
#endif /* PIO_ENABLE_LOGGING */

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

    LOG((1, "create_file_handler comproot = %d\n", ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    LOG((1, "create_file_handler got parameter len = %d\n", len));
    if (!(filename = malloc(len + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&mode, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    LOG((1, "create_file_handler got parameters len = %d "
    	   "filename = %s iotype = %d mode = %d\n",
	 len, filename, iotype, mode));

    /* Call the create file function. */
    if ((ret = PIOc_createfile(ios->iosysid, &ncid, &iotype, filename, mode)))
	return ret;

    /* Free resources. */
    free(filename);

    LOG((1, "create_file_handler succeeded!\n", my_rank));
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
    LOG((1, "%d close_file_handler\n", my_rank));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    LOG((1, "%d create_file_handler got parameter ncid = %d\n", ncid));

    /* Call the close file function. */
    if ((ret = PIOc_closefile(ncid)))
	return ret;

    LOG((1, "close_file_handler succeeded!\n", my_rank));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to inq a netCDF file. It is
 * only ever run on the IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @return PIO_NOERR for success, error code otherwise.
*/
int inq_handler(iosystem_desc_t *ios)
{
    int ncid;
    int ndims, nvars, ngatts, unlimdimid;
    int *ndimsp = NULL, *nvarsp = NULL, *ngattsp = NULL, *unlimdimidp = NULL;
    char ndims_present, nvars_present, ngatts_present, unlimdimid_present;
    int mpierr;
    int ret;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    LOG((1, "%d inq_handler\n", my_rank));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&ndims_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&nvars_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&ngatts_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&unlimdimid_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    LOG((1, "%d inq_handler ndims_present = %d nvars_present = %d ngatts_present = %d unlimdimid_present = %d\n",
	 ndims_present, nvars_present, ngatts_present, unlimdimid_present));

    /* NULLs passed in to any of the pointers in the original call
     * need to be matched with NULLs here. Assign pointers where
     * non-NULL pointers were passed in. */
    if (ndims_present)
	ndimsp = &ndims;
    if (nvars_present)
	nvarsp = &nvars;
    if (ngatts_present)
	ngattsp = &ngatts;
    if (unlimdimid_present)
	unlimdimidp = &unlimdimid;

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
    char name_present, len_present;
    char *dimnamep = NULL;
    PIO_Offset *dimlenp = NULL;
    char dimname[NC_MAX_NAME + 1];
    PIO_Offset dimlen;

    int mpierr;
    int ret;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    LOG((1, "inq_dim_handler\n", my_rank));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&dimid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&len_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d inq_handler name_present = %d len_present = %d\n",
	   name_present, len_present);

    /* Set the non-null pointers. */
    if (name_present)
	dimnamep = dimname;
    if (len_present)
	dimlenp = &dimlen;

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
    int *dimidp = NULL, dimid;
    int mpierr;
    int id_present;
    int ret;
    int namelen;
    char *name;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    LOG((1, "inq_dimid_handler\n", my_rank));

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
    if ((mpierr = MPI_Bcast(&id_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    LOG((1, "%d inq_dimid_handler ncid = %d namelen = %d name = %s id_present = %d\n",
	 ncid, namelen, name, id_present));

    /* Set non-null pointer. */
    if (id_present)
	dimidp = &dimid;

    /* Call the inq_dimid function. */
    if ((ret = PIOc_inq_dimid(ncid, name, dimidp)))
	return ret;

    /* Free resources. */
    free(name);

    return PIO_NOERR;
}

/** Handle attribute inquiry operations. This code only runs on IO
 * tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @return PIO_NOERR for success, error code otherwise.
*/
int inq_att_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    int ret;
    char *name5;
    int namelen;
    PIO_Offset attlen;
    int *op, *ip;
    nc_type xtype, *xtypep = NULL;
    PIO_Offset len, *lenp = NULL;
    char xtype_present, len_present;

    LOG((1, "inq_att_handler"));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm)))
	return PIO_EIO;
    if (!(name5 = malloc((namelen + 1) * sizeof(char))))
	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)name5, namelen + 1, MPI_CHAR, ios->compmaster,
			    ios->intercomm)))
	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast(&xtype_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&len_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;

    /* Match NULLs in collective function call. */
    if (xtype_present)
	xtypep = &xtype;
    if (len_present)
	lenp = &len;

    /* Call the function to learn about the attribute. */
    if ((ret = PIOc_inq_att(ncid, varid, name5, xtypep, lenp)))
	return ret;

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
    char *name5;
    int namelen;
    PIO_Offset len;
    nc_type xtype;
    int *op, *ip;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    LOG((1, "%d att_handler\n", my_rank));

    if (msg == PIO_MSG_PUT_ATT_INT)
    {
	/* Get the parameters for this function that the the comp master
	 * task is broadcasting. */
	if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	if (!(name5 = malloc((namelen + 1) * sizeof(char))))
	    return PIO_ENOMEM;
	mpierr = MPI_Bcast((void *)name5, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&xtype, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	mpierr = MPI_Bcast(&len, 1, MPI_OFFSET,  ios->compmaster, ios->intercomm);
	if (!(op = malloc(len * sizeof(int))))
	    return PIO_ENOMEM;
	mpierr = MPI_Bcast(op, len, MPI_INT,  ios->compmaster, ios->intercomm);

	/* Call the function to write the attribute. */
	if ((ret = PIOc_put_att_int(ncid, varid, name5, xtype, len, op)))
	    return ret;

	/* Free resources. */
	free(name5);
	free(op);
    }
    else if (msg = PIO_MSG_GET_ATT_INT)
    {
	/* Get the parameters for this function that the the comp master
	 * task is broadcasting. */
	if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
	    return PIO_EIO;
	mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm);
	if (!(name5 = malloc((namelen + 1) * sizeof(char))))
	    return PIO_ENOMEM;
	mpierr = MPI_Bcast((void *)name5, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm);

	/* Allocate space for the attribute data. */
        if ((ret = PIOc_inq_attlen(ncid, varid, name5, &len)))
	    return ret;
	if (!(ip = malloc(len * sizeof(int))))
	    return PIO_ENOMEM;

	/* Call the function to read the attribute. */
	if ((ret = PIOc_get_att_int(ncid, varid, name5, ip)))
	    return ret;

	/* Free resources. */
	free(name5);
	free(ip);
    }

    return PIO_NOERR;
}

/** Do an inq_var on a netCDF variable. This function is only run on
 * IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @return PIO_NOERR for success, error code otherwise.
*/
int inq_var_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    char name_present, xtype_present, ndims_present, dimids_present, natts_present;
    char name[NC_MAX_NAME + 1], *namep;
    nc_type xtype, *xtypep = NULL;
    int *ndimsp = NULL, *dimidsp = NULL, *nattsp = NULL;
    int ndims, dimids[NC_MAX_DIMS], natts;
    int ret;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    LOG((1, "%d inq_var_handler\n", my_rank));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&xtype_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&ndims_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&dimids_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&natts_present, 1, MPI_CHAR, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d inq_var_handler ncid = %d varid = %d name_present = %d xtype_present = %d ndims_present = %d "
	   "dimids_present = %d natts_present = %d\n",
	   my_rank, ncid, varid, name_present, xtype_present, ndims_present, dimids_present, natts_present);

    /* Set the non-NULL pointers. */
    if (name_present)
	namep = name;
    if (xtype_present)
	xtypep = &xtype;
    if (ndims_present)
	ndimsp = &ndims;
    if (dimids_present)
	dimidsp = dimids;
    if (natts_present)
	nattsp = &natts;

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
    LOG((1, "%d sync_file_handler\n", my_rank));

    /* Get the parameters for this function that the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    LOG((1, "%d sync_file_handler got parameter ncid = %d\n", my_rank, ncid));

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
    LOG((1, "%d enddef_file_handler\n", my_rank));

    /* Get the parameters for this function that the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;

    /* Call the sync file function. */
    if ((ret = PIOc_enddef(ncid)))
	return ret;

    LOG((1, "%d enddef_file_handler succeeded!\n", my_rank));
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
    LOG((1, "%d def_var_handler comproot = %d\n", my_rank, ios->comproot));

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
    LOG((1, "%d def_var_handler got parameters namelen = %d "
    	   "name = %s len = %d ncid = %d\n",
	 my_rank, namelen, name, len, ncid));

    /* Call the create file function. */
    if ((ret = PIOc_def_var(ncid, name, xtype, ndims, dimids, &varid)))
	return ret;

    /* Free resources. */
    free(name);
    free(dimids);

    LOG((1, "%d def_var_handler succeeded!\n", my_rank));
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
    LOG((1, "%d def_dim_handler comproot = %d\n", my_rank, ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    printf("%d def_dim_handler ncid = %d\n", my_rank, ncid);
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

    LOG((1, "%d def_dim_handler succeeded!\n", my_rank));
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
    LOG((1, "%d def_var_handler comproot = %d\n", my_rank, ios->comproot));

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

    LOG((1, "%d vara_handler succeeded!\n", my_rank));
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
    LOG((1, "%d open_file_handler comproot = %d\n", my_rank, ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    LOG((2, "open_file_handler got parameter len = %d", len));
    if (!(filename = malloc(len + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    if ((mpierr = MPI_Bcast(&mode, 1, MPI_INT, 0, ios->intercomm)))
    	return PIO_EIO;
    LOG((1, "%d open_file_handler got parameters len = %d "
    	   "filename = %s iotype = %d mode = %d\n",
	 my_rank, len, filename, iotype, mode));

    /* Call the open file function. */
    if ((ret = PIOc_openfile(ios->iosysid, &ncid, &iotype, filename, mode)))
	return ret;

    /* Free resources. */
    free(filename);

    LOG((1, "%d open_file_handler succeeded!\n", my_rank));
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
    LOG((1, "%d delete_file_handler comproot = %d\n", my_rank, ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
	return PIO_EIO;
    if (!(filename = malloc(len + 1 * sizeof(char))))
    	return PIO_ENOMEM;
    if ((mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, 0,
    			    ios->intercomm)))
    	return PIO_EIO;
    LOG((1, "%d delete_file_handler got parameters len = %d filename = %s\n",
	 my_rank, len, filename));

    /* Call the delete file function. */
    if ((ret = PIOc_deletefile(ios->iosysid, filename)))
	return ret;

    /* Free resources. */
    free(filename);

    LOG((1, "%d delete_file_handler succeeded!\n", my_rank));
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
    LOG((1, "%d finalize_handler called\n", my_rank));
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
    int ret;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    LOG((1, "%d pio_msg_handler called\n", my_rank));

    /* Have IO comm rank 0 (the ioroot) register to receive
     * (non-blocking) for a message from each of the comproots. */
    if (!io_rank)
    {
	for (int cmp = 0; cmp < component_count; cmp++)
	{
	    my_iosys = &iosys[cmp];
	    LOG((1, "%d about to call MPI_Irecv\n", my_rank));
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
	    LOG((1, "%d about to call MPI_Waitany req[0] = %d MPI_REQUEST_NULL = %d\n",
		 my_rank, req[0], MPI_REQUEST_NULL));
	    mpierr = MPI_Waitany(component_count, req, &index, &status);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    LOG((3, "Waitany returned index = %d req[%d] = %d",
		 index, index, req[index]));
	}

	/* Broadcast the index of the computational component that
	 * originated the request to the rest of the IO tasks. */
	mpierr = MPI_Bcast(&index, 1, MPI_INT, 0, iosys->io_comm);
	CheckMPIReturn(mpierr, __FILE__, __LINE__);
	my_iosys = &iosys[index];
	LOG((3, "index MPI_Bcast complete index = %d", index));

	/* Broadcast the msg value to the rest of the IO tasks. */
	LOG((3, "about to call msg MPI_Bcast"));
	mpierr = MPI_Bcast(&msg, 1, MPI_INT, 0, my_iosys->io_comm);
	CheckMPIReturn(mpierr, __FILE__, __LINE__);
	LOG((3, "msg MPI_Bcast complete msg = %d", msg));

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
	    inq_handler(my_iosys);
	    break;
	case PIO_MSG_INQ_DIM:
	    inq_dim_handler(my_iosys, msg);
	    break;
	case PIO_MSG_INQ_DIMID:
	    inq_dimid_handler(my_iosys);
	    break;
	case PIO_MSG_INQ_VAR:
	    inq_var_handler(my_iosys);
	    break;
	case PIO_MSG_GET_ATT_INT:
	case PIO_MSG_PUT_ATT_INT:
	    ret = att_handler(my_iosys, msg);
	    break;
	case PIO_MSG_INQ_VARID:
	    inq_varid_handler(my_iosys);
	    break;
	case PIO_MSG_INQ_ATT:
	    inq_att_handler(my_iosys);
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

	/* If an error was returned by the handler, do something! */
	if (ret)
	{
	    LOG((0, "hander returned error code %d", ret));
	    MPI_Finalize();
	}

	/* Unless finalize was called, listen for another msg from the
	 * component whose message we just handled. */
	if (!io_rank && msg != -1)
	{
	    my_iosys = &iosys[index];
	    mpierr = MPI_Irecv(&msg, 1, MPI_INT, my_iosys->comproot, MPI_ANY_TAG, my_iosys->union_comm,
			       &req[index]);
	    LOG((3, "pio_msg_handler called MPI_Irecv req[%d] = %d\n", index, req[index]));
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

/** @ingroup PIO_init
 * Library initialization used when IO tasks are distinct from compute
 * tasks.
 *
 * This is a collective call.  Input parameters are read on
 * comp_rank=0 values on other tasks are ignored.  This variation of
 * PIO_init sets up a distinct set of tasks to handle IO, these tasks
 * do not return from this call.  Instead they go to an internal loop
 * and wait to receive further instructions from the computational
 * tasks.
 *
 * For 4 tasks, to have 2 of them be computational, and 2 of them
 * be IO, I would provide the following:
 *
 * component_count = 1
 *
 * peer_comm = MPI_COMM_WORLD
 *
 * comp_comms = an array with one element, an MPI (intra) communicator
 * that contains the two tasks designated to do computation
 * (processors 0, 1).

 * io_comm = an MPI (intra) communicator with the other two tasks (2,
 * 3).
 *
 * iosysidp = pointer that gets the IO system ID.
 *
 * Fortran function (from PIO1, in piolib_mod.F90) is:
 *
 * subroutine init_intercom(component_count, peer_comm, comp_comms,
 * io_comm, iosystem, rearr_opts)
 *
 * Some notes from Jim:
 *
 * Components and Component Count
 * ------------------------------
 *
 * It's a cesm thing - the cesm model is composed of several component
 * models (atm, ocn, ice, lnd, etc) that may or may not be collocated
 * on mpi tasks.  Since for intercomm the IOCOMM tasks are a subset of
 * the compute tasks for a given component we have a separate iocomm
 * for each model component.  and we call init_inracomm independently
 * for each component.
 *
 * When the IO tasks are independent of any model component then we
 * can have all of the components share one set of iotasks and we call
 * init_intercomm once with the information for all components.
 *
 * Inter vs Intra Communicators
 * ----------------------------
 *
 * ​For an intra you just need to provide the compute comm, pio creates
 * an io comm as a subset of that compute comm.
 *
 * For an inter you need to provide multiple comms - peer comm is the
 * communicator that is going to encompass all of the tasks - usually
 * this will be mpi_comm_world.  Then you need to provide a comm for
 * each component model that will share the io server, then an
 * io_comm.
 *
 * Example of Communicators
 * ------------------------
 *
 * Starting from MPI_COMM_WORLD the calling program will create an
 * IO_COMM and one or more COMP_COMMs, I think an example might be best:
 *
 * Suppose we have 10 tasks and 2 of them will be IO tasks.  Then 0:7
 * are in COMP_COMM and 8:9 are in IO_COMM In this case on tasks 0:7
 * COMP_COMM is defined and IO_COMM is MPI_COMM_NULL and on tasks 8:9
 * IO_COMM is defined and COMP_COMM is MPI_COMM_NULL The communicators
 * to handle communications between COMP_COMM and IO_COMM are defined
 * in init_intercomm and held in a pio internal data structure.
 *
 * Return or Not
 * -------------
 *
 * The io_comm tasks do not return from the init_intercomm routine.
 *
 * Sequence of Events to do Asynch I/O
 * -----------------------------------
 *
 * Here is the sequence of events that needs to occur when an IO
 * operation is called from the collection of compute tasks.  I'm
 * going to use pio_put_var because write_darray has some special
 * characteristics that make it a bit more complicated...
 *
 * Compute tasks call pio_put_var with an integer argument
 *
 * The MPI_Send sends a message from comp_rank=0 to io_rank=0 on
 * union_comm (a comm defined as the union of io and compute tasks)
 * msg is an integer which indicates the function being called, in
 * this case the msg is PIO_MSG_PUT_VAR_INT
 *
 * The iotasks now know what additional arguments they should expect
 * to receive from the compute tasks, in this case a file handle, a
 * variable id, the length of the array and the array itself.
 *
 * The iotasks now have the information they need to complete the
 * operation and they call the pio_put_var routine.  (In pio1 this bit
 * of code is in pio_get_put_callbacks.F90.in)
 *
 * After the netcdf operation is completed (in the case of an inq or
 * get operation) the result is communicated back to the compute
 * tasks.
 *
 *
 * @param component_count The number of computational (ex. model)
 * components to associate with this IO component
 *
 * @param peer_comm The communicator from which all other communicator
 * arguments are derived
 *
 * @param comp_comms An array containing the computational
 * communicator for each of the computational components. The I/O
 * tasks pass MPI_COMM_NULL for this parameter.
 *
`* @param io_comm The io communicator. Processing tasks pass
 * MPI_COMM_NULL for this parameter.
 *
 * @param iosysidp An array of length component_count. It will get the
 * iosysid for each component.
 *
 * @return PIO_NOERR on success, error code otherwise.
 */
int PIOc_Init_Intercomm(int component_count, MPI_Comm peer_comm,
			MPI_Comm *comp_comms, MPI_Comm io_comm, int *iosysidp)
{
    iosystem_desc_t *iosys;
    iosystem_desc_t *my_iosys;
    int ierr = PIO_NOERR;
    int mpierr;
    int my_rank;
    int iam;
    int io_leader, comp_leader;
    int root;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Allocate struct to hold io system info for each component. */
    if (!(iosys = (iosystem_desc_t *) calloc(1, sizeof(iosystem_desc_t) * component_count)))
	ierr = PIO_ENOMEM;

    if (!ierr)
	for (int cmp = 0; cmp < component_count; cmp++)
	{
	    /* These are used when using the intercomm. */
	    int comp_master = MPI_PROC_NULL, io_master = MPI_PROC_NULL;

	    /* Get a pointer to the iosys struct */
	    my_iosys = &iosys[cmp];

	    /* Create an MPI info object. */
	    CheckMPIReturn(MPI_Info_create(&(my_iosys->info)),__FILE__,__LINE__);

	    /* This task is part of the computation communicator. */
	    if (comp_comms[cmp] != MPI_COMM_NULL)
	    {
		/* Copy the computation communicator. */
		mpierr = MPI_Comm_dup(comp_comms[cmp], &my_iosys->comp_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Create an MPI group with the computation tasks. */
		mpierr = MPI_Comm_group(my_iosys->comp_comm, &my_iosys->compgroup);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find out how many tasks are in this communicator. */
		mpierr = MPI_Comm_size(iosys->comp_comm, &my_iosys->num_comptasks);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Set the rank within the comp_comm. */
		mpierr = MPI_Comm_rank(my_iosys->comp_comm, &my_iosys->comp_rank);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the io leader in peer_comm. */
		iam = -1;
		mpierr = MPI_Allreduce(&iam, &io_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the comp leader in peer_comm. */
		if (!my_iosys->comp_rank)
		{
		    mpierr = MPI_Comm_rank(peer_comm, &iam);
		    CheckMPIReturn(mpierr, __FILE__, __LINE__);
		    if (mpierr)
			ierr = PIO_EIO;
		}
		else
		    iam = -1;

		/* Find the lucky comp_leader task. */
		mpierr = MPI_Allreduce(&iam, &comp_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Is this the compmaster? Only if the comp_rank is zero. */
		if (!my_iosys->comp_rank)
		{
		    my_iosys->compmaster = MPI_ROOT;
		    comp_master = MPI_ROOT;
		}
		else
		    my_iosys->compmaster = MPI_PROC_NULL;

		/* Set up the intercomm from the computation side. */
		mpierr = MPI_Intercomm_create(my_iosys->comp_comm, 0, peer_comm,
					      io_leader, cmp, &my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Create the union communicator. */
		mpierr = MPI_Intercomm_merge(my_iosys->intercomm, 0, &my_iosys->union_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;
	    }
	    else
	    {
		my_iosys->comp_comm = MPI_COMM_NULL;
		my_iosys->compgroup = MPI_GROUP_NULL;
		my_iosys->comp_rank = -1;
	    }

	    /* This task is part of the IO communicator, so set up the
	     * IO stuff. */
	    if (io_comm != MPI_COMM_NULL)
	    {
		/* Copy the IO communicator. */
		mpierr = MPI_Comm_dup(io_comm, &my_iosys->io_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Get an MPI group that includes the io tasks. */
		mpierr = MPI_Comm_group(my_iosys->io_comm, &my_iosys->iogroup);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find out how many tasks are in this communicator. */
		mpierr = MPI_Comm_size(iosys->io_comm, &my_iosys->num_iotasks);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Set the rank within the io_comm. */
		mpierr = MPI_Comm_rank(my_iosys->io_comm, &my_iosys->io_rank);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the io leader in peer_comm. */
		if (!my_iosys->io_rank)
		{
		    mpierr = MPI_Comm_rank(peer_comm, &iam);
		    CheckMPIReturn(mpierr, __FILE__, __LINE__);
		    if (mpierr)
			ierr = PIO_EIO;
		}
		else
		    iam = -1;

		/* Find the lucky io_leader task. */
		mpierr = MPI_Allreduce(&iam, &io_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the comp leader in peer_comm. */
		iam = -1;
		mpierr = MPI_Allreduce(&iam, &comp_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* This is an io task. */
		my_iosys->ioproc = true;

		/* Is this the iomaster? Only if the io_rank is zero. */
		if (!my_iosys->io_rank)
		{
		    my_iosys->iomaster = MPI_ROOT;
		    io_master = MPI_ROOT;
		}
		else
		    my_iosys->iomaster = 0;

		/* Set up the intercomm from the I/O side. */
		mpierr = MPI_Intercomm_create(my_iosys->io_comm, 0, peer_comm,
					      comp_leader, cmp, &my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

		/* Create the union communicator. */
		mpierr = MPI_Intercomm_merge(my_iosys->intercomm, 0, &my_iosys->union_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		if (mpierr)
		    ierr = PIO_EIO;

	    }
	    else
	    {
		my_iosys->io_comm = MPI_COMM_NULL;
		my_iosys->iogroup = MPI_GROUP_NULL;
		my_iosys->io_rank = -1;
		my_iosys->ioproc = false;
		my_iosys->iomaster = false;
	    }

	    /* my_comm points to the union communicator for async, and
	     * the comp_comm for non-async. It should not be freed
	     * since it is not a proper copy of the commuicator, just
	     * a copy of the reference to it. */
	    my_iosys->my_comm = my_iosys->union_comm;

	    /* Find rank in union communicator. */
	    mpierr = MPI_Comm_rank(my_iosys->union_comm, &my_iosys->union_rank);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    if (mpierr)
		ierr = PIO_EIO;

	    /* Find the rank of the io leader in the union communicator. */
	    if (!my_iosys->io_rank)
		my_iosys->ioroot = my_iosys->union_rank;
	    else
		my_iosys->ioroot = -1;

	    /* Distribute the answer to all tasks. */
	    mpierr = MPI_Allreduce(&my_iosys->ioroot, &root, 1, MPI_INT, MPI_MAX,
				   my_iosys->union_comm);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    if (mpierr)
		ierr = PIO_EIO;
	    my_iosys->ioroot = root;

	    /* Find the rank of the computation leader in the union
	     * communicator. */
	    if (!my_iosys->comp_rank)
		my_iosys->comproot = my_iosys->union_rank;
	    else
		my_iosys->comproot = -1;

	    /* Distribute the answer to all tasks. */
	    mpierr = MPI_Allreduce(&my_iosys->comproot, &root, 1, MPI_INT, MPI_MAX,
				   my_iosys->union_comm);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    if (mpierr)
		ierr = PIO_EIO;
	    my_iosys->comproot = root;

	    /* Send the number of tasks in the IO and computation
	       communicators to each other over the intercomm. This is
	       a one-to-all bcast from the local task that passes
	       MPI_ROOT as the root (all other local tasks should pass
	       MPI_PROC_NULL as the root). The bcast is recieved by
	       all the members of the leaf group which each pass the
	       rank of the root relative to the root group. */
	    if (io_comm != MPI_COMM_NULL)
	    {
		comp_master = 0;
		mpierr = MPI_Bcast(&my_iosys->num_comptasks, 1, MPI_INT, comp_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		mpierr = MPI_Bcast(&my_iosys->num_iotasks, 1, MPI_INT, io_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    }
	    else
	    {
		io_master = 0;
		mpierr = MPI_Bcast(&my_iosys->num_comptasks, 1, MPI_INT, comp_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
		mpierr = MPI_Bcast(&my_iosys->num_iotasks, 1, MPI_INT, io_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    }

	    /* Allocate an array to hold the ranks of the IO tasks
	     * within the union communicator. */
	    if (!(my_iosys->ioranks = malloc(my_iosys->num_iotasks * sizeof(int))))
		return PIO_ENOMEM;

	    /* Allocate a temp array to help get the IO ranks. */
	    int *tmp_ioranks;
	    if (!(tmp_ioranks = malloc(my_iosys->num_iotasks * sizeof(int))))
		return PIO_ENOMEM;

	    /* Init array, then have IO tasks set their values, then
	     * use allreduce to distribute results to all tasks. */
	    for (int cnt = 0 ; cnt < my_iosys->num_iotasks; cnt++)
		tmp_ioranks[cnt] = -1;
	    if (io_comm != MPI_COMM_NULL)
		tmp_ioranks[my_iosys->io_rank] = my_iosys->union_rank;
	    mpierr = MPI_Allreduce(tmp_ioranks, my_iosys->ioranks, my_iosys->num_iotasks, MPI_INT, MPI_MAX,
				   my_iosys->union_comm);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);

	    /* Free temp array. */
	    free(tmp_ioranks);

	    /* Set the default error handling. */
	    my_iosys->error_handler = PIO_INTERNAL_ERROR;

	    /* We do support asynch interface. */
	    my_iosys->async_interface = true;

	    /* For debug purposes, print the contents of the struct. */
	    /* for (int t = 0; t < my_iosys->num_iotasks + my_iosys->num_comptasks; t++) */
	    /* { */
	    /* 	MPI_Barrier(my_iosys->union_comm); */
	    /* 	if (my_rank == t) */
	    /* 	    pio_iosys_print(my_rank, my_iosys); */
	    /* } */

	    /* Add this id to the list of PIO iosystem ids. */
	    iosysidp[cmp] = pio_add_to_iosystem_list(my_iosys);
	    LOG((2, "added to iosystem_list iosysid = %d", iosysidp[cmp]));

	    /* Now call the function from which the IO tasks will not
	     * return until the PIO_MSG_EXIT message is sent. */
	    if (io_comm != MPI_COMM_NULL)
		if ((ierr = pio_msg_handler(my_iosys->io_rank, component_count, iosys)))
		    return ierr;
	}

    /* If there was an error, make sure all tasks see it. */
    if (ierr)
    {
	/*mpierr = MPI_Bcast(&ierr, 1, MPI_INT, iosys->my_rank, iosys->intercomm);*/
	mpierr = MPI_Bcast(&ierr, 1, MPI_INT, 0, iosys->intercomm);
	CheckMPIReturn(mpierr, __FILE__, __LINE__);
	if (mpierr)
	    ierr = PIO_EIO;
    }

    return ierr;
}
