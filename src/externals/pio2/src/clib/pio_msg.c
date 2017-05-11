/** @file
 *
 * PIO async message handling. This file contains the code which
 * runs on the IO nodes when async is in use. This code waits for
 * messages from the computation nodes, and responds to messages by
 * running the appropriate netCDF function.
 *
 * @author Ed Hartnett
 */

#include <config.h>
#include <pio.h>
#include <pio_internal.h>

#ifdef PIO_ENABLE_LOGGING
extern int my_rank;
extern int pio_log_level;
#endif /* PIO_ENABLE_LOGGING */

/** This function is run on the IO tasks to handle nc_inq_type*()
 * functions.
 *
 * @param ios pointer to the iosystem info.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_type_handler(iosystem_desc_t *ios)
{
    int ncid;
    int xtype;
    char name_present, size_present;
    char *namep = NULL, name[NC_MAX_NAME + 1];
    PIO_Offset *sizep = NULL, size;
    int mpierr;
    int ret;

    LOG((1, "inq_type_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&xtype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&size_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Handle null pointer issues. */
    if (name_present)
        namep = name;
    if (size_present)
        sizep = &size;

    /* Call the function. */
    if ((ret = PIOc_inq_type(ncid, xtype, namep, sizep)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "inq_type_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to find netCDF file
 * format.
 *
 * @param ios pointer to the iosystem info.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_format_handler(iosystem_desc_t *ios)
{
    int ncid;
    int *formatp = NULL, format;
    char format_present;
    int mpierr;
    int ret;

    LOG((1, "inq_format_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&format_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "inq_format_handler got parameters ncid = %d format_present = %d",
         ncid, format_present));

    /* Manage NULL pointers. */
    if (format_present)
        formatp = &format;

    /* Call the function. */
    if ((ret = PIOc_inq_format(ncid, formatp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    if (formatp)
        LOG((2, "inq_format_handler format = %d", *formatp));
    LOG((1, "inq_format_handler succeeded!"));

    return PIO_NOERR;
}

/** This function is run on the IO tasks to set the file fill mode.
 *
 * @param ios pointer to the iosystem info.
 * @returns 0 for success, error code otherwise.
 * @internal
 */
int set_fill_handler(iosystem_desc_t *ios)
{
    int ncid;
    int fillmode;
    int old_modep_present;
    int old_mode, *old_modep = NULL;
    int mpierr;
    int ret;

    LOG((1, "set_fill_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&fillmode, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&old_modep_present, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "set_fill_handler got parameters ncid = %d fillmode = %d old_modep_present = %d",
         ncid, fillmode, old_modep_present));

    /* Manage NULL pointers. */
    if (old_modep_present)
        old_modep = &old_mode;

    /* Call the function. */
    if ((ret = PIOc_set_fill(ncid, fillmode, old_modep)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "set_fill_handler succeeded!"));

    return PIO_NOERR;
}

/** This function is run on the IO tasks to create a netCDF file.
 *
 * @param ios pointer to the iosystem info.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int create_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len;
    int iotype;
    int mode;
    int mpierr;
    int ret;

    LOG((1, "create_file_handler comproot = %d", ios->comproot));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Get space for the filename. */
    char filename[len + 1];

    if ((mpierr = MPI_Bcast(filename, len + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&mode, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "create_file_handler got parameters len = %d filename = %s iotype = %d mode = %d",
         len, filename, iotype, mode));

    /* Call the create file function. */
    if ((ret = PIOc_createfile(ios->iosysid, &ncid, &iotype, filename, mode)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "create_file_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to close a netCDF file. It is
 * only ever run on the IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int close_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int mpierr;
    int ret;

    LOG((1, "close_file_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "create_file_handler got parameter ncid = %d", ncid));

    /* Call the close file function. */
    if ((ret = PIOc_closefile(ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "close_file_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to inq a netCDF file. It is
 * only ever run on the IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_handler(iosystem_desc_t *ios)
{
    int ncid;
    int ndims, nvars, ngatts, unlimdimid;
    int *ndimsp = NULL, *nvarsp = NULL, *ngattsp = NULL, *unlimdimidp = NULL;
    char ndims_present, nvars_present, ngatts_present, unlimdimid_present;
    int mpierr;
    int ret;

    LOG((1, "inq_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ndims_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&nvars_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ngatts_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&unlimdimid_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "inq_handler ndims_present = %d nvars_present = %d ngatts_present = %d unlimdimid_present = %d",
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
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Do an inq_dim on a netCDF dimension. This function is only run on
 * IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
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

    LOG((1, "inq_dim_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&dimid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&len_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "inq_handler name_present = %d len_present = %d", name_present,
         len_present));

    /* Set the non-null pointers. */
    if (name_present)
        dimnamep = dimname;
    if (len_present)
        dimlenp = &dimlen;

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq_dim(ncid, dimid, dimnamep, dimlenp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Do an inq_dimid on a netCDF dimension name. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_dimid_handler(iosystem_desc_t *ios)
{
    int ncid;
    int *dimidp = NULL, dimid;
    int mpierr;
    int id_present;
    int ret;
    int namelen;
    char name[PIO_MAX_NAME + 1];

    LOG((1, "inq_dimid_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&id_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "inq_dimid_handler ncid = %d namelen = %d name = %s id_present = %d",
         ncid, namelen, name, id_present));

    /* Set non-null pointer. */
    if (id_present)
        dimidp = &dimid;

    /* Call the inq_dimid function. */
    if ((ret = PIOc_inq_dimid(ncid, name, dimidp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Handle attribute inquiry operations. This code only runs on IO
 * tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_att_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    int ret;
    char name[PIO_MAX_NAME + 1];
    int namelen;
    nc_type xtype, *xtypep = NULL;
    PIO_Offset len, *lenp = NULL;
    char xtype_present, len_present;

    LOG((1, "inq_att_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, ios->compmaster,
                            ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&xtype_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&len_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Match NULLs in collective function call. */
    if (xtype_present)
        xtypep = &xtype;
    if (len_present)
        lenp = &len;

    /* Call the function to learn about the attribute. */
    if ((ret = PIOc_inq_att(ncid, varid, name, xtypep, lenp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Handle attribute inquiry operations. This code only runs on IO
 * tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_attname_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int attnum;
    char name[NC_MAX_NAME + 1], *namep = NULL;
    char name_present;
    int mpierr;
    int ret;

    LOG((1, "inq_att_name_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&attnum, 1, MPI_INT,  ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "inq_attname_handler got ncid = %d varid = %d attnum = %d name_present = %d",
         ncid, varid, attnum, name_present));

    /* Match NULLs in collective function call. */
    if (name_present)
        namep = name;

    /* Call the function to learn about the attribute. */
    if ((ret = PIOc_inq_attname(ncid, varid, attnum, namep)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Handle attribute inquiry operations. This code only runs on IO
 * tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_attid_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    char name[PIO_MAX_NAME + 1];
    int namelen;
    int id, *idp = NULL;
    char id_present;
    int mpierr;
    int ret;

    LOG((1, "inq_attid_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR,  ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&id_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "inq_attid_handler got ncid = %d varid = %d id_present = %d",
         ncid, varid, id_present));

    /* Match NULLs in collective function call. */
    if (id_present)
        idp = &id;

    /* Call the function to learn about the attribute. */
    if ((ret = PIOc_inq_attid(ncid, varid, name, idp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Handle attribute operations. This code only runs on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int att_put_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    int ret;
    char name[PIO_MAX_NAME + 1];
    int namelen;
    PIO_Offset attlen;  /* Number of elements in att array. */
    nc_type atttype;    /* Type of att in file. */
    PIO_Offset atttype_len; /* Length in bytes of one elementy of type atttype. */
    nc_type memtype;    /* Type of att data in memory. */
    PIO_Offset memtype_len; /* Length of element of memtype. */
    void *op;

    LOG((1, "att_put_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&atttype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&attlen, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&atttype_len, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&memtype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&memtype_len, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Allocate memory for the attribute data. */
    if (!(op = malloc(attlen * memtype_len)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(op, attlen * memtype_len, MPI_BYTE, 0, ios->intercomm)))
    {
        free(op);
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    }
    LOG((1, "att_put_handler ncid = %d varid = %d namelen = %d name = %s"
         "atttype = %d attlen = %d atttype_len = %d memtype = %d memtype_len = 5d",
         ncid, varid, namelen, name, atttype, attlen, atttype_len, memtype, memtype_len));

    /* Call the function to write the attribute. */
    ret = PIOc_put_att_tc(ncid, varid, name, atttype, attlen, memtype, op);

    /* Free resources. */
    free(op);

    /* Did it work? */
    if (ret)
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((2, "att_put_handler complete!"));
    return PIO_NOERR;
}

/** Handle attribute operations. This code only runs on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @param msg the message sent my the comp root task.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int att_get_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    char name[PIO_MAX_NAME + 1];
    int namelen;
    PIO_Offset attlen;
    nc_type atttype;        /* Type of att in file. */
    PIO_Offset atttype_len; /* Length in bytes of an element of attype. */
    nc_type memtype;        /* Type of att in memory. */
    PIO_Offset memtype_len; /* Length in bytes of an element of memype. */
    int *ip;
    int iotype;
    int ret;

    LOG((1, "att_get_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT,  ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, ios->compmaster, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&atttype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&attlen, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&atttype_len, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&memtype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&memtype_len, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "att_get_handler ncid = %d varid = %d namelen = %d name = %s iotype = %d"
         " atttype = %d attlen = %d atttype_len = %d memtype = %d memtype_len = %d",
         ncid, varid, namelen, name, iotype, atttype, attlen, atttype_len, memtype, memtype_len));

    /* Allocate space for the attribute data. */
    if (!(ip = malloc(attlen * memtype_len)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Call the function to read the attribute. */
    ret = PIOc_get_att_tc(ncid, varid, name, memtype, ip);

    /* Free resources. */
    free(ip);

    /* Did it work? */
    if (ret)
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Handle var put operations. This code only runs on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int put_vars_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    PIO_Offset typelen;  /** Length (in bytes) of this type. */
    nc_type xtype;       /** Type of the data being written. */
    char stride_present; /** Zero if user passed a NULL stride. */
    PIO_Offset *startp = NULL;
    PIO_Offset *countp = NULL;
    PIO_Offset *stridep = NULL;
    int ndims;           /** Number of dimensions. */
    void *buf;           /** Buffer for data storage. */
    PIO_Offset num_elem; /** Number of data elements in the buffer. */
    int mpierr;          /** Error code from MPI function calls. */

    LOG((1, "put_vars_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Now we know how big to make these arrays. */
    PIO_Offset start[ndims], count[ndims], stride[ndims];

    if ((mpierr = MPI_Bcast(start, ndims, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "put_vars_handler getting start[0] = %d ndims = %d", start[0], ndims));
    if ((mpierr = MPI_Bcast(count, ndims, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&stride_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (stride_present)
        if ((mpierr = MPI_Bcast(stride, ndims, MPI_OFFSET, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&xtype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&num_elem, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "put_vars_handler ncid = %d varid = %d ndims = %d "
         "stride_present = %d xtype = %d num_elem = %d typelen = %d",
         ncid, varid, ndims, stride_present, xtype, num_elem, typelen));

    for (int d = 0; d < ndims; d++)
    {
        LOG((2, "start[%d] = %d", d, start[d]));
        LOG((2, "count[%d] = %d", d, count[d]));
        if (stride_present)
            LOG((2, "stride[%d] = %d", d, stride[d]));
    }

    /* Allocate room for our data. */
    if (!(buf = malloc(num_elem * typelen)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Get the data. */
    if ((mpierr = MPI_Bcast(buf, num_elem * typelen, MPI_BYTE, 0, ios->intercomm)))
    {
        free(buf);
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    }

    /* for (int e = 0; e < num_elem; e++) */
    /*  LOG((2, "element %d = %d", e, ((int *)buf)[e])); */

    /* Set the non-NULL pointers. */
    startp = start;
    countp = count;
    if (stride_present)
        stridep = stride;

    /* Call the function to write the data. No need to check return
     * values, they are bcast to computation tasks inside function. */
    switch(xtype)
    {
    case NC_BYTE:
        PIOc_put_vars_schar(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_CHAR:
        PIOc_put_vars_text(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_SHORT:
        PIOc_put_vars_short(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_INT:
        PIOc_put_vars_int(ncid, varid, startp, countp, stridep, buf);
        break;
    case PIO_LONG_INTERNAL:
        PIOc_put_vars_long(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_FLOAT:
        PIOc_put_vars_float(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_DOUBLE:
        PIOc_put_vars_double(ncid, varid, startp, countp, stridep, buf);
        break;
#ifdef _NETCDF4
    case NC_UBYTE:
        PIOc_put_vars_uchar(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_USHORT:
        PIOc_put_vars_ushort(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_UINT:
        PIOc_put_vars_uint(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_INT64:
        PIOc_put_vars_longlong(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_UINT64:
        PIOc_put_vars_ulonglong(ncid, varid, startp, countp, stridep, buf);
        break;
        /* case NC_STRING: */
        /*      PIOc_put_vars_string(ncid, varid, startp, countp, */
        /*                                stridep, (void *)buf); */
        /*      break; */
        /*    default:*/
        /* PIOc_put_vars(ncid, varid, startp, countp, */
        /*                   stridep, buf); */
#endif /* _NETCDF4 */
    }

    free(buf);

    return PIO_NOERR;
}

/** Handle var get operations. This code only runs on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int get_vars_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    PIO_Offset typelen; /** Length (in bytes) of this type. */
    nc_type xtype; /** Type of the data being written. */
    char stride_present;
    PIO_Offset *startp = NULL, *countp = NULL, *stridep = NULL;
    int ndims; /** Number of dimensions. */
    void *buf; /** Buffer for data storage. */
    PIO_Offset num_elem; /** Number of data elements in the buffer. */

    LOG((1, "get_vars_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Now we know how big to make these arrays. */
    PIO_Offset start[ndims], count[ndims], stride[ndims];

    if ((mpierr = MPI_Bcast(start, ndims, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "put_vars_handler getting start[0] = %d ndims = %d", start[0], ndims));
    if ((mpierr = MPI_Bcast(count, ndims, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&stride_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (stride_present)
        if ((mpierr = MPI_Bcast(stride, ndims, MPI_OFFSET, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&xtype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&num_elem, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&typelen, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "get_vars_handler ncid = %d varid = %d ndims = %d "
         "stride_present = %d xtype = %d num_elem = %d typelen = %d",
         ncid, varid, ndims, stride_present, xtype, num_elem, typelen));

    for (int d = 0; d < ndims; d++)
    {
        LOG((2, "start[%d] = %d", d, start[d]));
        LOG((2, "count[%d] = %d", d, count[d]));
        if (stride_present)
            LOG((2, "stride[%d] = %d", d, stride[d]));
    }

    /* Allocate room for our data. */
    if (!(buf = malloc(num_elem * typelen)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Set the non-NULL pointers. */
    startp = start;
    countp = count;
    if (stride_present)
        stridep = stride;

    /* Call the function to read the data. */
    switch(xtype)
    {
    case NC_BYTE:
        PIOc_get_vars_schar(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_CHAR:
        PIOc_get_vars_text(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_SHORT:
        PIOc_get_vars_short(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_INT:
        PIOc_get_vars_int(ncid, varid, startp, countp, stridep, buf);
        break;
    case PIO_LONG_INTERNAL:
        PIOc_get_vars_long(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_FLOAT:
        PIOc_get_vars_float(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_DOUBLE:
        PIOc_get_vars_double(ncid, varid, startp, countp, stridep, buf);
        break;
#ifdef _NETCDF4
    case NC_UBYTE:
        PIOc_get_vars_uchar(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_USHORT:
        PIOc_get_vars_ushort(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_UINT:
        PIOc_get_vars_uint(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_INT64:
        PIOc_get_vars_longlong(ncid, varid, startp, countp, stridep, buf);
        break;
    case NC_UINT64:
        PIOc_get_vars_ulonglong(ncid, varid, startp, countp, stridep, buf);
        break;
        /* case NC_STRING: */
        /*      PIOc_get_vars_string(ncid, varid, startp, countp, */
        /*                                stridep, (void *)buf); */
        /*      break; */
        /*    default:*/
        /* PIOc_get_vars(ncid, varid, startp, countp, */
        /*                   stridep, buf); */
#endif /* _NETCDF4 */
    }

    free(buf);
    LOG((1, "get_vars_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * Do an inq_var on a netCDF variable. This function is only run on
 * IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int inq_var_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    char name_present, xtype_present, ndims_present, dimids_present, natts_present;
    char name[NC_MAX_NAME + 1], *namep = NULL;
    nc_type xtype, *xtypep = NULL;
    int *ndimsp = NULL, *dimidsp = NULL, *nattsp = NULL;
    int ndims, dimids[NC_MAX_DIMS], natts;
    int ret;

    LOG((1, "inq_var_handler"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&name_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&xtype_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ndims_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&dimids_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&natts_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2,"inq_var_handler ncid = %d varid = %d name_present = %d xtype_present = %d ndims_present = %d "
         "dimids_present = %d natts_present = %d",
         ncid, varid, name_present, xtype_present, ndims_present, dimids_present, natts_present));

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
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    if (ndims_present)
        LOG((2, "inq_var_handler ndims = %d", ndims));

    return PIO_NOERR;
}

/**
 * Do an inq_var_chunking on a netCDF variable. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int inq_var_chunking_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    char storage_present, chunksizes_present;
    int storage, *storagep = NULL;
    PIO_Offset chunksizes[NC_MAX_DIMS], *chunksizesp = NULL;
    int mpierr;
    int ret;

    assert(ios);
    LOG((1, "inq_var_chunking_handler"));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&storage_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&chunksizes_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2,"inq_var_handler ncid = %d varid = %d storage_present = %d chunksizes_present = %d",
         ncid, varid, storage_present, chunksizes_present));

    /* Set the non-NULL pointers. */
    if (storage_present)
        storagep = &storage;
    if (chunksizes_present)
        chunksizesp = chunksizes;

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq_var_chunking(ncid, varid, storagep, chunksizesp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/**
 * Do an inq_var_fill on a netCDF variable. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int inq_var_fill_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    char fill_mode_present, fill_value_present;
    PIO_Offset type_size;
    int fill_mode, *fill_modep = NULL;
    PIO_Offset *fill_value, *fill_valuep = NULL;
    int mpierr;

    assert(ios);
    LOG((1, "inq_var_fill_handler"));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&type_size, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&fill_mode_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&fill_value_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2,"inq_var_fill_handler ncid = %d varid = %d type_size = %lld, fill_mode_present = %d fill_value_present = %d",
         ncid, varid, type_size, fill_mode_present, fill_value_present));

    /* If we need to, alocate storage for fill value. */
    if (fill_value_present)
        if (!(fill_value = malloc(type_size)))
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Set the non-NULL pointers. */
    if (fill_mode_present)
        fill_modep = &fill_mode;
    if (fill_value_present)
        fill_valuep = fill_value;

    /* Call the inq function to get the values. */
    PIOc_inq_var_fill(ncid, varid, fill_modep, fill_valuep);

    /* Free fill value storage if we allocated some. */
    if (fill_value_present)
        free(fill_value);

    return PIO_NOERR;
}

/**
 * Do an inq_var_endian on a netCDF variable. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int inq_var_endian_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    char endian_present;
    int endian, *endianp = NULL;
    int mpierr;
    int ret;

    assert(ios);
    LOG((1, "inq_var_endian_handler"));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&endian_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2,"inq_var_endian_handler ncid = %d varid = %d endian_present = %d", ncid, varid,
         endian_present));

    /* Set the non-NULL pointers. */
    if (endian_present)
        endianp = &endian;

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq_var_endian(ncid, varid, endianp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/**
 * Do an inq_var_deflate on a netCDF variable. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int inq_var_deflate_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    char shuffle_present;
    char deflate_present;
    char deflate_level_present;
    int shuffle, *shufflep;
    int deflate, *deflatep;
    int deflate_level, *deflate_levelp;
    int mpierr;
    int ret;

    assert(ios);
    LOG((1, "inq_var_deflate_handler"));

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&shuffle_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (shuffle_present && !mpierr)
        if ((mpierr = MPI_Bcast(&shuffle, 1, MPI_INT, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&deflate_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (deflate_present && !mpierr)
        if ((mpierr = MPI_Bcast(&deflate, 1, MPI_INT, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&deflate_level_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (deflate_level_present && !mpierr)
        if ((mpierr = MPI_Bcast(&deflate_level, 1, MPI_INT, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "inq_var_handler ncid = %d varid = %d shuffle_present = %d deflate_present = %d "
         "deflate_level_present = %d", ncid, varid, shuffle_present, deflate_present,
         deflate_level_present));

    /* Set the non-NULL pointers. */
    if (shuffle_present)
        shufflep = &shuffle;
    if (deflate_present)
        deflatep = &deflate;
    if (deflate_level_present)
        deflate_levelp = &deflate_level;

    /* Call the inq function to get the values. */
    if ((ret = PIOc_inq_var_deflate(ncid, varid, shufflep, deflatep, deflate_levelp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** Do an inq_varid on a netCDF variable name. This function is only
 * run on IO tasks.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int inq_varid_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int mpierr;
    int ret;
    int namelen;
    char name[PIO_MAX_NAME + 1];

    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Call the inq_dimid function. */
    if ((ret = PIOc_inq_varid(ncid, name, &varid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/** This function is run on the IO tasks to sync a netCDF file.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int sync_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int mpierr;
    int ret;

    LOG((1, "sync_file_handler"));
    assert(ios);

    /* Get the parameters for this function that the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "sync_file_handler got parameter ncid = %d", ncid));

    /* Call the sync file function. */
    if ((ret = PIOc_sync(ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((2, "sync_file_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to enddef a netCDF file.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int change_def_file_handler(iosystem_desc_t *ios, int msg)
{
    int ncid;
    int mpierr;

    LOG((1, "change_def_file_handler"));
    assert(ios);

    /* Get the parameters for this function that the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Call the function. */
    if (msg == PIO_MSG_ENDDEF)
        PIOc_enddef(ncid);
    else
        PIOc_redef(ncid);

    LOG((1, "change_def_file_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to define a netCDF
 *  variable.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int def_var_handler(iosystem_desc_t *ios)
{
    int ncid;
    int namelen;
    char name[PIO_MAX_NAME + 1];
    int mpierr;
    int ret;
    int varid;
    nc_type xtype;
    int ndims;
    int *dimids;

    LOG((1, "def_var_handler comproot = %d", ios->comproot));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&xtype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (!(dimids = malloc(ndims * sizeof(int))))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(dimids, ndims, MPI_INT, 0, ios->intercomm)))
    {
        free(dimids);
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    }
    LOG((1, "def_var_handler got parameters namelen = %d "
         "name = %s ncid = %d", namelen, name, ncid));

    /* Call the function. */
    if ((ret = PIOc_def_var(ncid, name, xtype, ndims, dimids, &varid)))
    {
        free(dimids);
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    }

    /* Free resources. */
    free(dimids);

    LOG((1, "def_var_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to define chunking for a
 *  netCDF variable.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int def_var_chunking_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int ndims;
    int storage;
    char chunksizes_present;
    PIO_Offset chunksizes[NC_MAX_DIMS], *chunksizesp = NULL;
    int mpierr;
    int ret;

    assert(ios);
    LOG((1, "def_var_chunking_handler comproot = %d", ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&storage, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&chunksizes_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (chunksizes_present)
        if ((mpierr = MPI_Bcast(chunksizes, ndims, MPI_OFFSET, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "def_var_chunking_handler got parameters ncid = %d varid = %d storage = %d "
         "ndims = %d chunksizes_present = %d", ncid, varid, storage, ndims, chunksizes_present));

    /* Set the non-NULL pointers. */
    if (chunksizes_present)
        chunksizesp = chunksizes;

    /* Call the function. */
    if ((ret = PIOc_def_var_chunking(ncid, varid, storage, chunksizesp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "def_var_chunking_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to define fill mode and fill
 * value.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int def_var_fill_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int fill_mode;
    char fill_value_present;
    PIO_Offset type_size;
    PIO_Offset *fill_valuep = NULL;
    int mpierr;

    assert(ios);
    LOG((1, "def_var_fill_handler comproot = %d", ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&fill_mode, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&type_size, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&fill_value_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if (fill_value_present)
    {
        if (!(fill_valuep = malloc(type_size)))
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(fill_valuep, type_size, MPI_CHAR, 0, ios->intercomm)))
        {
            free(fill_valuep);
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
        }
    }
    LOG((1, "def_var_fill_handler got parameters ncid = %d varid = %d fill_mode = %d "
         "type_size = %lld fill_value_present = %d", ncid, varid, fill_mode, type_size, fill_value_present));

    /* Call the function. */
    PIOc_def_var_fill(ncid, varid, fill_mode, fill_valuep);

    /* Free memory allocated for the fill value. */
    if (fill_valuep)
        free(fill_valuep);

    LOG((1, "def_var_fill_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to define endianness for a
 * netCDF variable.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int def_var_endian_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int endian;
    int mpierr;
    int ret;

    assert(ios);
    LOG((1, "def_var_endian_handler comproot = %d", ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&endian, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "def_var_endian_handler got parameters ncid = %d varid = %d endain = %d ",
         ncid, varid, endian));

    /* Call the function. */
    if ((ret = PIOc_def_var_endian(ncid, varid, endian)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "def_var_chunking_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to define deflate settings for
 * a netCDF variable.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int def_var_deflate_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int shuffle;
    int deflate;
    int deflate_level;
    int mpierr;
    int ret;

    assert(ios);
    LOG((1, "def_var_deflate_handler comproot = %d", ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&shuffle, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&deflate, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&deflate_level, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "def_var_deflate_handler got parameters ncid = %d varid = %d shuffle = %d ",
         "deflate = %d deflate_level = %d", ncid, varid, shuffle, deflate, deflate_level));

    /* Call the function. */
    if ((ret = PIOc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "def_var_deflate_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to define chunk cache settings
 * for a netCDF variable.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, error code otherwise.
 */
int set_var_chunk_cache_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    PIO_Offset size;
    PIO_Offset nelems;
    float preemption;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int ret; /* Return code. */

    assert(ios);
    LOG((1, "set_var_chunk_cache_handler comproot = %d", ios->comproot));

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&size, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&nelems, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&preemption, 1, MPI_FLOAT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "set_var_chunk_cache_handler got params ncid = %d varid = %d size = %d "
         "nelems = %d preemption = %g", ncid, varid, size, nelems, preemption));

    /* Call the function. */
    if ((ret = PIOc_set_var_chunk_cache(ncid, varid, size, nelems, preemption)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "def_var_chunk_cache_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to define a netCDF
 * dimension.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int def_dim_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len, namelen;
    char name[PIO_MAX_NAME + 1];
    int mpierr;
    int ret;
    int dimid;

    LOG((1, "def_dim_handler comproot = %d", ios->comproot));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "def_dim_handler got parameters namelen = %d "
         "name = %s len = %d ncid = %d", namelen, name, len, ncid));

    /* Call the function. */
    if ((ret = PIOc_def_dim(ncid, name, len, &dimid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "def_dim_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to rename a netCDF
 * dimension.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int rename_dim_handler(iosystem_desc_t *ios)
{
    int ncid;
    int namelen;
    char name[PIO_MAX_NAME + 1];
    int mpierr;
    int ret;
    int dimid;

    LOG((1, "rename_dim_handler"));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&dimid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "rename_dim_handler got parameters namelen = %d "
         "name = %s ncid = %d dimid = %d", namelen, name, ncid, dimid));

    /* Call the function. */
    if ((ret = PIOc_rename_dim(ncid, dimid, name)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "rename_dim_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to rename a netCDF
 * dimension.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int rename_var_handler(iosystem_desc_t *ios)
{
    int ncid;
    int namelen;
    char name[PIO_MAX_NAME + 1];
    int mpierr;
    int ret;
    int varid;

    LOG((1, "rename_var_handler"));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "rename_var_handler got parameters namelen = %d "
         "name = %s ncid = %d varid = %d", namelen, name, ncid, varid));

    /* Call the function. */
    if ((ret = PIOc_rename_var(ncid, varid, name)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "rename_var_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to rename a netCDF
 * attribute.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int rename_att_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int namelen, newnamelen;
    char name[PIO_MAX_NAME + 1], newname[PIO_MAX_NAME + 1];
    int mpierr;
    int ret;

    LOG((1, "rename_att_handler"));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&newnamelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(newname, newnamelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "rename_att_handler got parameters namelen = %d name = %s ncid = %d varid = %d "
         "newnamelen = %d newname = %s", namelen, name, ncid, varid, newnamelen, newname));

    /* Call the function. */
    if ((ret = PIOc_rename_att(ncid, varid, name, newname)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "rename_att_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to delete a netCDF
 * attribute.
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int delete_att_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    int namelen;
    char name[PIO_MAX_NAME + 1];
    int mpierr;
    int ret;

    LOG((1, "delete_att_handler"));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&namelen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(name, namelen + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "delete_att_handler namelen = %d name = %s ncid = %d varid = %d ",
         namelen, name, ncid, varid));

    /* Call the function. */
    if ((ret = PIOc_del_att(ncid, varid, name)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "delete_att_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to open a netCDF file.
 *
 *
 * @param ios pointer to the iosystem_desc_t.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int open_file_handler(iosystem_desc_t *ios)
{
    int ncid;
    int len;
    int iotype;
    int mode;
    int mpierr;

    LOG((1, "open_file_handler comproot = %d", ios->comproot));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "open_file_handler got parameter len = %d", len));

    /* Get space for the filename. */
    char filename[len + 1];

    if ((mpierr = MPI_Bcast(filename, len + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&mode, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    LOG((2, "open_file_handler got parameters len = %d filename = %s iotype = %d mode = %d",
         len, filename, iotype, mode));

    /* Call the open file function. Errors are handling within
     * function, so return code can be ignored. */
    PIOc_openfile(ios->iosysid, &ncid, &iotype, filename, mode);

    return PIO_NOERR;
}

/** This function is run on the IO tasks to delete a netCDF file.
 *
 * @param ios pointer to the iosystem_desc_t data.
 *
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int delete_file_handler(iosystem_desc_t *ios)
{
    int len;
    int mpierr;
    int ret;

    LOG((1, "delete_file_handler comproot = %d", ios->comproot));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&len, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Get space for the filename. */
    char filename[len + 1];

    if ((mpierr = MPI_Bcast(filename, len + 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "delete_file_handler got parameters len = %d filename = %s",
         len, filename));

    /* Call the delete file function. */
    if ((ret = PIOc_deletefile(ios->iosysid, filename)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "delete_file_handler succeeded!"));
    return PIO_NOERR;
}

/** 
 * This function is run on the IO tasks to initialize a decomposition.
 *
 * @param ios pointer to the iosystem_desc_t data.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 */
int initdecomp_dof_handler(iosystem_desc_t *ios)
{
    int iosysid;
    int basetype;
    int ndims;
    int maplen;
    int ioid;
    char rearranger_present;
    int rearranger;
    int *rearrangerp = NULL;
    char iostart_present;
    PIO_Offset *iostartp = NULL;
    char iocount_present;
    PIO_Offset *iocountp = NULL;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int ret; /* Return code. */

    LOG((1, "initdecomp_dof_handler called"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&iosysid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&basetype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ndims, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Now we know the size of these arrays. */
    int dims[ndims];
    PIO_Offset iostart[ndims];
    PIO_Offset iocount[ndims];

    if ((mpierr = MPI_Bcast(dims, ndims, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&maplen, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    PIO_Offset compmap[maplen];

    if ((mpierr = MPI_Bcast(compmap, maplen, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    if ((mpierr = MPI_Bcast(&rearranger_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    if (rearranger_present)
        if ((mpierr = MPI_Bcast(&rearranger, 1, MPI_INT, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    if ((mpierr = MPI_Bcast(&iostart_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    if (iostart_present)
        if ((mpierr = MPI_Bcast(iostart, ndims, MPI_OFFSET, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    if ((mpierr = MPI_Bcast(&iocount_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    if (iocount_present)
        if ((mpierr = MPI_Bcast(iocount, ndims, MPI_OFFSET, 0, ios->intercomm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    LOG((2, "initdecomp_dof_handler iosysid = %d basetype = %d ndims = %d maplen = %d "
         "rearranger_present = %d iostart_present = %d iocount_present = %d ",
         iosysid, basetype, ndims, maplen, rearranger_present, iostart_present, iocount_present));

    if (rearranger_present)
        rearrangerp = &rearranger;
    if (iostart_present)
        iostartp = iostart;
    if (iocount_present)
        iocountp = iocount;

    /* Call the function. */
    ret = PIOc_InitDecomp(iosysid, basetype, ndims, dims, maplen, compmap, &ioid, rearrangerp,
                          iostartp, iocountp);
    
    LOG((1, "PIOc_InitDecomp returned %d", ret));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to...
 * NOTE: not yet implemented
 *
 * @param ios pointer to the iosystem_desc_t data.
 *
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int writedarray_handler(iosystem_desc_t *ios)
{
    assert(ios);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to...
 * NOTE: not yet implemented
 *
 * @param ios pointer to the iosystem_desc_t data.
 *
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int readdarray_handler(iosystem_desc_t *ios)
{
    assert(ios);
    return PIO_NOERR;
}

/** This function is run on the IO tasks to set the error handler.
 * NOTE: not yet implemented
 *
 * @param ios pointer to the iosystem_desc_t data.
 *
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int seterrorhandling_handler(iosystem_desc_t *ios)
{
    int method;
    int old_method_present;
    int old_method;
    int *old_methodp = NULL;
    int mpierr;
    int ret;

    LOG((1, "seterrorhandling_handler comproot = %d", ios->comproot));
    assert(ios);

    /* Get the parameters for this function that the he comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&method, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&old_method_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    LOG((1, "seterrorhandling_handler got parameters method = %d old_method_present = %d",
         method, old_method_present));

    if (old_method_present)
        old_methodp = &old_method;

    /* Call the function. */
    if ((ret = PIOc_set_iosystem_error_handling(ios->iosysid, method, old_methodp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "seterrorhandling_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to set the chunk cache
 * parameters for netCDF-4.
 *
 * @param ios pointer to the iosystem_desc_t data.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 */
int set_chunk_cache_handler(iosystem_desc_t *ios)
{
    int iosysid;
    int iotype;
    PIO_Offset size;
    PIO_Offset nelems;
    float preemption;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int ret; /* Return code. */

    LOG((1, "set_chunk_cache_handler called"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&iosysid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&size, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&nelems, 1, MPI_OFFSET, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&preemption, 1, MPI_FLOAT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "set_chunk_cache_handler got params iosysid = %d iotype = %d size = %d "
         "nelems = %d preemption = %g", iosysid, iotype, size, nelems, preemption));

    /* Call the function. */
    if ((ret = PIOc_set_chunk_cache(iosysid, iotype, size, nelems, preemption)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "set_chunk_cache_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to get the chunk cache
 * parameters for netCDF-4.
 *
 * @param ios pointer to the iosystem_desc_t data.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 */
int get_chunk_cache_handler(iosystem_desc_t *ios)
{
    int iosysid;
    int iotype;
    char size_present, nelems_present, preemption_present;
    PIO_Offset size, *sizep;
    PIO_Offset nelems, *nelemsp;
    float preemption, *preemptionp;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int ret; /* Return code. */

    LOG((1, "get_chunk_cache_handler called"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&iosysid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&iotype, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&size_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&nelems_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&preemption_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "get_chunk_cache_handler got params iosysid = %d iotype = %d size_present = %d "
         "nelems_present = %d preemption_present = %g", iosysid, iotype, size_present,
         nelems_present, preemption_present));

    /* Set the non-NULL pointers. */
    if (size_present)
        sizep = &size;
    if (nelems_present)
        nelemsp = &nelems;
    if (preemption_present)
        preemptionp = &preemption;

    /* Call the function. */
    if ((ret = PIOc_get_chunk_cache(iosysid, iotype, sizep, nelemsp, preemptionp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "get_chunk_cache_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is run on the IO tasks to get the variable chunk
 * cache parameters for netCDF-4.
 *
 * @param ios pointer to the iosystem_desc_t data.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 */
int get_var_chunk_cache_handler(iosystem_desc_t *ios)
{
    int ncid;
    int varid;
    char size_present, nelems_present, preemption_present;
    PIO_Offset size, *sizep;
    PIO_Offset nelems, *nelemsp;
    float preemption, *preemptionp;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int ret; /* Return code. */

    LOG((1, "get_var_chunk_cache_handler called"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&varid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&size_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&nelems_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&preemption_present, 1, MPI_CHAR, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "get_var_chunk_cache_handler got params ncid = %d varid = %d size_present = %d "
         "nelems_present = %d preemption_present = %g", ncid, varid, size_present,
         nelems_present, preemption_present));

    /* Set the non-NULL pointers. */
    if (size_present)
        sizep = &size;
    if (nelems_present)
        nelemsp = &nelems;
    if (preemption_present)
        preemptionp = &preemption;

    /* Call the function. */
    if ((ret = PIOc_get_var_chunk_cache(ncid, varid, sizep, nelemsp, preemptionp)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "get_var_chunk_cache_handler succeeded!"));
    return PIO_NOERR;
}

/** This function is run on the IO tasks to free the decomp hanlder.
 *
 * @param ios pointer to the iosystem_desc_t data.
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 */
int freedecomp_handler(iosystem_desc_t *ios)
{
    int iosysid;
    int ioid;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int ret; /* Return code. */

    LOG((1, "freedecomp_handler called"));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&iosysid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Bcast(&ioid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "freedecomp_handler iosysid = %d ioid = %d", iosysid, ioid));

    /* Call the function. */
    ret = PIOc_freedecomp(iosysid, ioid);
    
    LOG((1, "PIOc_freedecomp returned %d", ret));
    return PIO_NOERR;
}

/** Handle the finalize call.
 *
 * @param ios pointer to the iosystem info
 * @param index
 * @returns 0 for success, PIO_EIO for MPI Bcast errors, or error code
 * from netCDF base function.
 * @internal
 */
int finalize_handler(iosystem_desc_t *ios, int index)
{
    int iosysid;
    int mpierr;
    int ret;

    LOG((1, "finalize_handler called index = %d", index));
    assert(ios);

    /* Get the parameters for this function that the the comp master
     * task is broadcasting. */
    if ((mpierr = MPI_Bcast(&iosysid, 1, MPI_INT, 0, ios->intercomm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((1, "finalize_handler got parameter iosysid = %d", iosysid));

    /* Call the function. */
    LOG((2, "finalize_handler calling PIOc_finalize for iosysid = %d",
         iosysid));
    if ((ret = PIOc_finalize(iosysid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    LOG((1, "finalize_handler succeeded!"));
    return PIO_NOERR;
}

/**
 * This function is called by the IO tasks.  This function will not
 * return, unless there is an error.
 *
 * @param io_rank
 * @param component_count number of computation components
 * @param iosys pointer to pointer to iosystem info
 * @param io_comm MPI communicator for IO
 * @returns 0 for success, error code otherwise.
 */
int pio_msg_handler2(int io_rank, int component_count, iosystem_desc_t **iosys,
                     MPI_Comm io_comm)
{
    iosystem_desc_t *my_iosys;
    int msg = 0;
    MPI_Request req[component_count];
    MPI_Status status;
    int index;
    int mpierr;
    int ret = PIO_NOERR;
    int open_components = component_count;

    LOG((1, "pio_msg_handler2 called"));
    assert(iosys);

    /* Have IO comm rank 0 (the ioroot) register to receive
     * (non-blocking) for a message from each of the comproots. */
    if (!io_rank)
    {
        for (int cmp = 0; cmp < component_count; cmp++)
        {
            my_iosys = iosys[cmp];
            LOG((1, "about to call MPI_Irecv union_comm = %d", my_iosys->union_comm));
            if ((mpierr = MPI_Irecv(&msg, 1, MPI_INT, my_iosys->comproot, MPI_ANY_TAG,
                                    my_iosys->union_comm, &req[cmp])))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            LOG((1, "MPI_Irecv req[%d] = %d", cmp, req[cmp]));
        }
    }

    /* If the message is not -1, keep processing messages. */
    while (msg != -1)
    {
        LOG((3, "pio_msg_handler2 at top of loop"));

        /* Wait until any one of the requests are complete. Once it
         * returns, the Waitany function automatically sets the
         * appropriate member of the req array to MPI_REQUEST_NULL. */
        if (!io_rank)
        {
            LOG((1, "about to call MPI_Waitany req[0] = %d MPI_REQUEST_NULL = %d",
                 req[0], MPI_REQUEST_NULL));
            for (int c = 0; c < component_count; c++)
                LOG((2, "req[%d] = %d", c, req[c]));
            if ((mpierr = MPI_Waitany(component_count, req, &index, &status)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            LOG((3, "Waitany returned index = %d req[%d] = %d", index, index, req[index]));
        }

        /* Broadcast the index of the computational component that
         * originated the request to the rest of the IO tasks. */
        LOG((3, "About to do Bcast of index = %d io_comm = %d", index, io_comm));
        if ((mpierr = MPI_Bcast(&index, 1, MPI_INT, 0, io_comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        LOG((3, "index MPI_Bcast complete index = %d", index));

        /* Set the correct iosys depending on the index. */
        my_iosys = iosys[index];

        /* Broadcast the msg value to the rest of the IO tasks. */
        LOG((3, "about to call msg MPI_Bcast my_iosys->io_comm = %d", my_iosys->io_comm));
        if ((mpierr = MPI_Bcast(&msg, 1, MPI_INT, 0, my_iosys->io_comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        LOG((1, "pio_msg_handler2 msg MPI_Bcast complete msg = %d", msg));

        /* Handle the message. This code is run on all IO tasks. */
        switch (msg)
        {
        case PIO_MSG_INQ_TYPE:
            inq_type_handler(my_iosys);
            break;
        case PIO_MSG_INQ_FORMAT:
            inq_format_handler(my_iosys);
            break;
        case PIO_MSG_CREATE_FILE:
            create_file_handler(my_iosys);
            LOG((2, "returned from create_file_handler"));
            break;
        case PIO_MSG_SYNC:
            sync_file_handler(my_iosys);
            break;
        case PIO_MSG_ENDDEF:
        case PIO_MSG_REDEF:
            LOG((2, "calling change_def_file_handler"));
            change_def_file_handler(my_iosys, msg);
            LOG((2, "returned from change_def_file_handler"));
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
        case PIO_MSG_RENAME_DIM:
            rename_dim_handler(my_iosys);
            break;
        case PIO_MSG_RENAME_VAR:
            rename_var_handler(my_iosys);
            break;
        case PIO_MSG_RENAME_ATT:
            rename_att_handler(my_iosys);
            break;
        case PIO_MSG_DEL_ATT:
            delete_att_handler(my_iosys);
            break;
        case PIO_MSG_DEF_DIM:
            def_dim_handler(my_iosys);
            break;
        case PIO_MSG_DEF_VAR:
            def_var_handler(my_iosys);
            break;
        case PIO_MSG_DEF_VAR_CHUNKING:
            def_var_chunking_handler(my_iosys);
            break;
        case PIO_MSG_DEF_VAR_FILL:
            def_var_fill_handler(my_iosys);
            break;
        case PIO_MSG_DEF_VAR_ENDIAN:
            def_var_endian_handler(my_iosys);
            break;
        case PIO_MSG_DEF_VAR_DEFLATE:
            def_var_deflate_handler(my_iosys);
            break;
        case PIO_MSG_INQ_VAR_ENDIAN:
            inq_var_endian_handler(my_iosys);
            break;
        case PIO_MSG_SET_VAR_CHUNK_CACHE:
            set_var_chunk_cache_handler(my_iosys);
            break;
        case PIO_MSG_GET_VAR_CHUNK_CACHE:
            get_var_chunk_cache_handler(my_iosys);
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
        case PIO_MSG_INQ_VAR_CHUNKING:
            inq_var_chunking_handler(my_iosys);
            break;
        case PIO_MSG_INQ_VAR_FILL:
            inq_var_fill_handler(my_iosys);
            break;
        case PIO_MSG_INQ_VAR_DEFLATE:
            inq_var_deflate_handler(my_iosys);
            break;
        case PIO_MSG_GET_ATT:
            ret = att_get_handler(my_iosys);
            break;
        case PIO_MSG_PUT_ATT:
            ret = att_put_handler(my_iosys);
            break;
        case PIO_MSG_INQ_VARID:
            inq_varid_handler(my_iosys);
            break;
        case PIO_MSG_INQ_ATT:
            inq_att_handler(my_iosys);
            break;
        case PIO_MSG_INQ_ATTNAME:
            inq_attname_handler(my_iosys);
            break;
        case PIO_MSG_INQ_ATTID:
            inq_attid_handler(my_iosys);
            break;
        case PIO_MSG_GET_VARS:
            get_vars_handler(my_iosys);
            break;
        case PIO_MSG_PUT_VARS:
            put_vars_handler(my_iosys);
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
        case PIO_MSG_SET_CHUNK_CACHE:
            set_chunk_cache_handler(my_iosys);
            break;
        case PIO_MSG_GET_CHUNK_CACHE:
            get_chunk_cache_handler(my_iosys);
            break;
        case PIO_MSG_FREEDECOMP:
            freedecomp_handler(my_iosys);
            break;
        case PIO_MSG_SET_FILL:
            set_fill_handler(my_iosys);
            break;
        case PIO_MSG_EXIT:
            finalize_handler(my_iosys, index);
            msg = -1;
            break;
        default:
            LOG((0, "unknown message received %d", msg));
            return PIO_EINVAL;
        }

        /* If an error was returned by the handler, do nothing! */
        LOG((3, "pio_msg_handler2 checking error ret = %d", ret));

        /* Listen for another msg from the component whose message we
         * just handled. */
        if (!io_rank && msg != -1)
        {
            my_iosys = iosys[index];
            LOG((3, "pio_msg_handler2 about to Irecv index = %d comproot = %d union_comm = %d",
                 index, my_iosys->comproot, my_iosys->union_comm));
            if ((mpierr = MPI_Irecv(&msg, 1, MPI_INT, my_iosys->comproot, MPI_ANY_TAG, my_iosys->union_comm,
                                    &req[index])))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            LOG((3, "pio_msg_handler2 called MPI_Irecv req[%d] = %d", index, req[index]));
        }

        LOG((3, "pio_msg_handler2 done msg = %d open_components = %d",
             msg, open_components));

        /* If there are no more open components, exit. */
        if (msg == -1)
            if (--open_components)
                msg = PIO_MSG_EXIT;
    }

    LOG((3, "returning from pio_msg_handler2"));
    return PIO_NOERR;
}
