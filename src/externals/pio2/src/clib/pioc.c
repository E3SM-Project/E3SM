/**
 * @file
 * Some initialization and support functions.
 * @author Jim Edwards
 * @date  2014
 *
 * @see http://code.google.com/p/parallelio/
 */

#include <config.h>
#include <pio.h>
#include <pio_internal.h>

/** The default error handler used when iosystem cannot be located. */
int default_error_handler = PIO_INTERNAL_ERROR;

/** The target blocksize for each io task when the box rearranger is
 * used (see pio_sc.c). */
extern int blocksize;

/**
 * Check to see if PIO has been initialized.
 *
 * @param iosysid the IO system ID
 * @param active pointer that gets true if IO system is active, false
 * otherwise.
 * @returns 0 on success, error code otherwise
 */
int PIOc_iosystem_is_active(int iosysid, bool *active)
{
    iosystem_desc_t *ios;

    /* Get the ios if there is one. */
    ios = pio_get_iosystem_from_id(iosysid);

    if (active)
    {
        if (!ios || (ios->comp_comm == MPI_COMM_NULL && ios->io_comm == MPI_COMM_NULL))
            *active = false;
        else
            *active = true;
    }

    return PIO_NOERR;
}

/**
 * Check to see if PIO file is open.
 *
 * @param ncid the ncid of an open file
 * @returns 1 if file is open, 0 otherwise.
 */
int PIOc_File_is_Open(int ncid)
{
    file_desc_t *file;

    /* If get file returns non-zero, then this file is not open. */
    if (pio_get_file(ncid, &file))
        return 0;
    else
        return 1;
}

/**
 * Set the error handling method to be used for subsequent pio library
 * calls, returns the previous method setting. Note that this changes
 * error handling for the IO system that was used when this file was
 * opened. Other files opened with the same IO system will also he
 * affected by this call. This function is supported but
 * deprecated. New code should use PIOc_set_iosystem_error_handling().
 * This method has no way to return an error, so any failure will
 * result in MPI_Abort.
 *
 * @param ncid the ncid of an open file
 * @param method the error handling method
 * @returns old error handler
 * @ingroup PIO_error_method
 */
int PIOc_Set_File_Error_Handling(int ncid, int method)
{
    file_desc_t *file;
    int oldmethod;

    /* Get the file info. */
    if (pio_get_file(ncid, &file))
        piodie("Could not find file", __FILE__, __LINE__);

    /* Check that valid error handler was provided. */
    if (method != PIO_INTERNAL_ERROR && method != PIO_BCAST_ERROR &&
        method != PIO_RETURN_ERROR)
        piodie("Invalid error hanlder method", __FILE__, __LINE__);

    /* Get the old method. */
    oldmethod = file->iosystem->error_handler;

    /* Set the error hanlder. */
    file->iosystem->error_handler = method;

    return oldmethod;
}

/**
 * Increment the unlimited dimension of the given variable.
 *
 * @param ncid the ncid of the open file
 * @param varid the variable ID
 * @returns 0 on success, error code otherwise
 */
int PIOc_advanceframe(int ncid, int varid)
{
    file_desc_t *file;
    int ret;

    /* Get the file info. */
    if ((ret = pio_get_file(ncid, &file)))
        return pio_err(NULL, NULL, ret, __FILE__, __LINE__);

    /* Check inputs. */
    if (varid < 0 || varid >= PIO_MAX_VARS)
        return pio_err(NULL, file, PIO_EINVAL, __FILE__, __LINE__);

    file->varlist[varid].record++;

    return PIO_NOERR;
}

/**
 * Set the unlimited dimension of the given variable
 *
 * @param ncid the ncid of the file.
 * @param varid the varid of the variable
 * @param frame the value of the unlimited dimension.  In c 0 for the
 * first record, 1 for the second
 * @return PIO_NOERR for no error, or error code.
 * @ingroup PIO_setframe
 */
int PIOc_setframe(int ncid, int varid, int frame)
{
    file_desc_t *file;
    int ret;

    /* Get file info. */
    if ((ret = pio_get_file(ncid, &file)))
        return pio_err(NULL, NULL, ret, __FILE__, __LINE__);

    /* Check inputs. */
    if (varid < 0 || varid >= PIO_MAX_VARS)
        return pio_err(NULL, file, PIO_EINVAL, __FILE__, __LINE__);

    file->varlist[varid].record = frame;

    return PIO_NOERR;
}

/**
 * Get the number of IO tasks set.
 *
 * @param iosysid the IO system ID
 * @param numiotasks a pointer taht gets the number of IO
 * tasks. Ignored if NULL.
 * @returns 0 on success, error code otherwise
 */
int PIOc_get_numiotasks(int iosysid, int *numiotasks)
{
    iosystem_desc_t *ios;

    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    if (numiotasks)
        *numiotasks = ios->num_iotasks;

    return PIO_NOERR;
}

/**
 * Get the local size of the variable.
 *
 * @param ioid IO descrption ID.
 * @returns the size of the array.
 */
int PIOc_get_local_array_size(int ioid)
{
    io_desc_t *iodesc;

    if (!(iodesc = pio_get_iodesc_from_id(ioid)))
        piodie("Could not get iodesc", __FILE__, __LINE__);

    return iodesc->ndof;
}

/**
 * Set the error handling method used for subsequent calls. This
 * function is deprecated. New code should use
 * PIOc_set_iosystem_error_handling(). This method has no way to
 * return an error, so any failure will result in MPI_Abort.
 *
 * @param iosysid the IO system ID
 * @param method the error handling method
 * @returns old error handler
 * @ingroup PIO_error_method
 */
int PIOc_Set_IOSystem_Error_Handling(int iosysid, int method)
{
    iosystem_desc_t *ios;
    int oldmethod;

    /* Get the iosystem info. */
    if (iosysid != PIO_DEFAULT)
        if (!(ios = pio_get_iosystem_from_id(iosysid)))
            return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Set the error handler. */
    if (PIOc_set_iosystem_error_handling(iosysid, method, &oldmethod))
        piodie("Could not set the IOSystem error hanlder", __FILE__, __LINE__);

    return oldmethod;
}

/**
 * Set the error handling method used for subsequent calls for this IO
 * system.
 *
 * @param iosysid the IO system ID. Passing PIO_DEFAULT instead
 * changes the default error handling for the library.
 * @param method the error handling method
 * @param old_method pointer to int that will get old method. Ignored
 * if NULL.
 * @returns 0 for success, error code otherwise.
 * @ingroup PIO_error_method
 */
int PIOc_set_iosystem_error_handling(int iosysid, int method, int *old_method)
{
    iosystem_desc_t *ios = NULL;
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */

    LOG((1, "PIOc_set_iosystem_error_handling iosysid = %d method = %d", iosysid,
         method));

    /* Find info about this iosystem. */
    if (iosysid != PIO_DEFAULT)
        if (!(ios = pio_get_iosystem_from_id(iosysid)))
            return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Check that valid error handler was provided. */
    if (method != PIO_INTERNAL_ERROR && method != PIO_BCAST_ERROR &&
        method != PIO_RETURN_ERROR)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* If using async, and not an IO task, then send parameters. */
    if (iosysid != PIO_DEFAULT)
        if (ios->async_interface)
        {
            if (!ios->ioproc)
            {
                int msg = PIO_MSG_SETERRORHANDLING;
                char old_method_present = old_method ? true : false;

                if (ios->compmaster == MPI_ROOT)
                    mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

                if (!mpierr)
                    mpierr = MPI_Bcast(&method, 1, MPI_INT, ios->compmaster, ios->intercomm);
                if (!mpierr)
                    mpierr = MPI_Bcast(&old_method_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            }

            /* Handle MPI errors. */
            if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
                check_mpi2(ios, NULL, mpierr2, __FILE__, __LINE__);
            if (mpierr)
                return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
        }

    /* Return the current handler. */
    if (old_method)
        *old_method = iosysid == PIO_DEFAULT ? default_error_handler : ios->error_handler;

    /* Set new error handler. */
    if (iosysid == PIO_DEFAULT)
        default_error_handler = method;
    else
        ios->error_handler = method;

    return PIO_NOERR;
}

/**
 * Initialize the decomposition used with distributed arrays. The
 * decomposition describes how the data will be distributed between
 * tasks.
 *
 * @param iosysid the IO system ID.
 * @param basetype the basic PIO data type used.
 * @param ndims the number of dimensions in the variable, not
 * including the unlimited dimension.
 * @param dims an array of global size of each dimension.
 * @param maplen the local length of the compmap array.
 * @param compmap a 1 based array of offsets into the array record on
 * file. A 0 in this array indicates a value which should not be
 * transfered.
 * @param ioidp pointer that will get the io description ID.
 * @param rearranger pointer to the rearranger to be used for this
 * decomp or NULL to use the default.
 * @param iostart An array of start values for block cyclic
 * decompositions. If NULL ???
 * @param iocount An array of count values for block cyclic
 * decompositions. If NULL ???
 * @returns 0 on success, error code otherwise
 * @ingroup PIO_initdecomp
 */
int PIOc_InitDecomp(int iosysid, int basetype, int ndims, const int *dims, int maplen,
                    const PIO_Offset *compmap, int *ioidp, const int *rearranger,
                    const PIO_Offset *iostart, const PIO_Offset *iocount)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    io_desc_t *iodesc;     /* The IO description. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function calls. */
    int ierr;              /* Return code. */

    LOG((1, "PIOc_InitDecomp iosysid = %d basetype = %d ndims = %d maplen = %d",
         iosysid, basetype, ndims, maplen));

    /* Get IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Caller must provide these. */
    if (!dims || !compmap || !ioidp)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* Check the dim lengths. */
    for (int i = 0; i < ndims; i++)
        if (dims[i] <= 0)
            return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_INITDECOMP_DOF; /* Message for async notification. */
            char rearranger_present = rearranger ? true : false;
            char iostart_present = iostart ? true : false;
            char iocount_present = iocount ? true : false;

            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&iosysid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&basetype, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&ndims, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((int *)dims, ndims, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&maplen, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((PIO_Offset *)compmap, maplen, MPI_OFFSET, ios->compmaster, ios->intercomm);

            if (!mpierr)
                mpierr = MPI_Bcast(&rearranger_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (rearranger_present && !mpierr)
                mpierr = MPI_Bcast((int *)rearranger, 1, MPI_INT, ios->compmaster, ios->intercomm);

            if (!mpierr)
                mpierr = MPI_Bcast(&iostart_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (iostart_present && !mpierr)
                mpierr = MPI_Bcast((PIO_Offset *)iostart, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);

            if (!mpierr)
                mpierr = MPI_Bcast(&iocount_present, 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (iocount_present && !mpierr)
                mpierr = MPI_Bcast((PIO_Offset *)iocount, ndims, MPI_OFFSET, ios->compmaster, ios->intercomm);
            LOG((2, "PIOc_InitDecomp iosysid = %d basetype = %d ndims = %d maplen = %d rearranger_present = %d iostart_present = %d "
                 "iocount_present = %d ", iosysid, basetype, ndims, maplen, rearranger_present, iostart_present, iocount_present));
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi2(ios, NULL, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    }

    /* Allocate space for the iodesc info. */
    if ((ierr = malloc_iodesc(ios, basetype, ndims, &iodesc)))
        return pio_err(ios, NULL, ierr, __FILE__, __LINE__);

    /* Remember the maplen. */
    iodesc->maplen = maplen;

    /* Remember the map. */
    if (!(iodesc->map = malloc(sizeof(PIO_Offset) * maplen)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
    for (int m = 0; m < maplen; m++)
        iodesc->map[m] = compmap[m];

    /* Remember the dim sizes. */
    if (!(iodesc->dimlen = malloc(sizeof(int) * ndims)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
    for (int d = 0; d < ndims; d++)
        iodesc->dimlen[d] = dims[d];

    /* Set the rearranger. */
    if (!rearranger)
        iodesc->rearranger = ios->default_rearranger;
    else
        iodesc->rearranger = *rearranger;
    LOG((2, "iodesc->rearranger = %d", iodesc->rearranger));

    /* Is this the subset rearranger? */
    if (iodesc->rearranger == PIO_REARR_SUBSET)
    {
        LOG((2, "Handling subset rearranger."));
        if (iostart && iocount)
            fprintf(stderr,"%s %s\n","Iostart and iocount arguments to PIOc_InitDecomp",
                    "are incompatable with subset rearrange method and will be ignored");
        iodesc->num_aiotasks = ios->num_iotasks;
        if ((ierr = subset_rearrange_create(ios, maplen, (PIO_Offset *)compmap, dims,
                                            ndims, iodesc)))
            return pio_err(NULL, NULL, ierr, __FILE__, __LINE__);
    }
    else
    {
        LOG((2, "Handling not the subset rearranger."));
        if (ios->ioproc)
        {
            /*  Unless the user specifies the start and count for each
             *  IO task compute it. */
            if (iostart && iocount)
            {
                iodesc->maxiobuflen = 1;
                for (int i = 0; i < ndims; i++)
                {
                    iodesc->firstregion->start[i] = iostart[i];
                    iodesc->firstregion->count[i] = iocount[i];
                    compute_maxIObuffersize(ios->io_comm, iodesc);

                }
                iodesc->num_aiotasks = ios->num_iotasks;
            }
            else
            {
                iodesc->num_aiotasks = CalcStartandCount(basetype, ndims, dims,
                                                         ios->num_iotasks, ios->io_rank,
                                                         iodesc->firstregion->start, iodesc->firstregion->count);
            }
            compute_maxIObuffersize(ios->io_comm, iodesc);
        }

        /* Depending on array size and io-blocksize the actual number
         * of io tasks used may vary. */
        if ((mpierr = MPI_Bcast(&(iodesc->num_aiotasks), 1, MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
        LOG((3, "iodesc->num_aiotasks = %d", iodesc->num_aiotasks));

        /* Compute the communications pattern for this decomposition. */
        if (iodesc->rearranger == PIO_REARR_BOX)
            ierr = box_rearrange_create(ios, maplen, compmap, dims, ndims, iodesc);

        /*
          if (ios->ioproc){
          io_region *ioregion = iodesc->firstregion;
          while(ioregion != NULL){
          for (int i=0;i<ndims;i++)
          printf("%s %d i %d dim %d start %ld count %ld\n",__FILE__,__LINE__,i,dims[i],ioregion->start[i],ioregion->count[i]);
          ioregion = ioregion->next;
          }
          }
        */
    }

    /* Add this IO description to the list. */
    *ioidp = pio_add_to_iodesc_list(iodesc);

    LOG((3, "About to tune rearranger..."));
    performance_tune_rearranger(ios, iodesc);
    LOG((3, "Done with rearranger tune."));

    return PIO_NOERR;
}

/**
 * Initialize the decomposition used with distributed arrays. The
 * decomposition describes how the data will be distributed between
 * tasks.
 *
 * @param iosysid the IO system ID.
 * @param basetype the basic PIO data type used.
 * @param ndims the number of dimensions in the variable, not
 * including the unlimited dimension.
 * @param dims an array of global size of each dimension.
 * @param maplen the local length of the compmap array.
 * @param compmap a 0 based array of offsets into the array record on
 * file. A -1 in this array indicates a value which should not be
 * transfered.
 * @param ioidp pointer that will get the io description ID.
 * @param rearranger pointer to the rearranger to be used for this
 * decomp or NULL to use the default.
 * @param iostart An array of start values for block cyclic
 * decompositions. If NULL ???
 * @param iocount An array of count values for block cyclic
 * decompositions. If NULL ???
 * @returns 0 on success, error code otherwise
 * @ingroup PIO_initdecomp
 */
int PIOc_init_decomp(int iosysid, int basetype, int ndims, const int *dims, int maplen,
                     const PIO_Offset *compmap, int *ioidp, const int *rearranger,
                     const PIO_Offset *iostart, const PIO_Offset *iocount)
{
    PIO_Offset compmap_1_based[maplen];

    LOG((1, "PIOc_init_decomp iosysid = %d basetype = %d ndims = %d maplen = %d",
         iosysid, basetype, ndims, maplen));

    /* Add 1 to all elements in compmap. */
    for (int e = 0; e < maplen; e++)
        compmap_1_based[e] = compmap[e] + 1;

    return PIOc_InitDecomp(iosysid, basetype, ndims, dims, maplen, compmap_1_based,
                           ioidp, rearranger, iostart, iocount);
}

/**
 * This is a simplified initdecomp which can be used if the memory
 * order of the data can be expressed in terms of start and count on
 * the file. In this case we compute the compdof.
 *
 * @param iosysid the IO system ID
 * @param basetype
 * @param ndims the number of dimensions
 * @param dims array of dimensions
 * @param start start array
 * @param count count array
 * @param pointer that gets the IO ID.
 * @returns 0 for success, error code otherwise
 * @ingroup PIO_initdecomp
 */
int PIOc_InitDecomp_bc(int iosysid, int basetype, int ndims, const int *dims,
                       const long int *start, const long int *count, int *ioidp)

{
    iosystem_desc_t *ios;
    int n, i, maplen = 1;
    PIO_Offset prod[ndims], loc[ndims];
    int rearr = PIO_REARR_SUBSET;

    LOG((1, "PIOc_InitDecomp_bc iosysid = %d basetype = %d ndims = %d"));

    /* Get the info about the io system. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Check for required inputs. */
    if (!dims || !start || !count || !ioidp)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* Check that dim, start, and count values are not obviously
     * incorrect. */
    for (int i = 0; i < ndims; i++)
        if (dims[i] <= 0 || start[i] < 0 || count[i] < 0 || (start[i] + count[i]) > dims[i])
            return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* Find the maplen. */
    for (i = 0; i < ndims; i++)
        maplen *= count[i];

    /* Get storage for the compmap. */
    PIO_Offset compmap[maplen];

    /* Find the compmap. */
    prod[ndims - 1] = 1;
    loc[ndims - 1] = 0;
    for (n = ndims - 2; n >= 0; n--)
    {
        prod[n] = prod[n + 1] * dims[n + 1];
        loc[n] = 0;
    }
    for (i = 0; i < maplen; i++)
    {
        compmap[i] = 0;
        for (n = ndims - 1; n >= 0; n--)
            compmap[i] += (start[n] + loc[n]) * prod[n];

        n = ndims - 1;
        loc[n] = (loc[n] + 1) % count[n];
        while (loc[n] == 0 && n > 0)
        {
            n--;
            loc[n] = (loc[n] + 1) % count[n];
        }
    }

    return PIOc_InitDecomp(iosysid, basetype, ndims, dims, maplen, compmap, ioidp,
                           &rearr, NULL, NULL);
}

/**
 * Library initialization used when IO tasks are a subset of compute
 * tasks.
 *
 * This function creates an MPI intracommunicator between a set of IO
 * tasks and one or more sets of computational tasks.
 *
 * The caller must create all comp_comm and the io_comm MPI
 * communicators before calling this function.
 *
 * @param comp_comm the MPI_Comm of the compute tasks
 * @param num_iotasks the number of io tasks to use
 * @param stride the offset between io tasks in the comp_comm
 * @param base the comp_comm index of the first io task
 * @param rearr the rearranger to use by default, this may be
 * overriden in the @ref PIO_initdecomp
 * @param iosysidp index of the defined system descriptor
 * @return 0 on success, otherwise a PIO error code.
 * @ingroup PIO_init
 */
int PIOc_Init_Intracomm(MPI_Comm comp_comm, int num_iotasks, int stride, int base,
                        int rearr, int *iosysidp)
{
    iosystem_desc_t *ios;
    int ustride;
    int num_comptasks; /* The size of the comp_comm. */
    int mpierr;        /* Return value for MPI calls. */
    int ret;           /* Return code for function calls. */

    /* Turn on the logging system. */
    pio_init_logging();

    /* Find the number of computation tasks. */
    if ((mpierr = MPI_Comm_size(comp_comm, &num_comptasks)))
        return check_mpi2(NULL, NULL, mpierr, __FILE__, __LINE__);

    /* Check the inputs. */
    if (!iosysidp || num_iotasks < 1 || num_iotasks * stride > num_comptasks)
        return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((1, "PIOc_Init_Intracomm comp_comm = %d num_iotasks = %d stride = %d base = %d "
         "rearr = %d", comp_comm, num_iotasks, stride, base, rearr));

    /* Allocate memory for the iosystem info. */
    if (!(ios = calloc(1, sizeof(iosystem_desc_t))))
        return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    ios->io_comm = MPI_COMM_NULL;
    ios->intercomm = MPI_COMM_NULL;
    ios->error_handler = default_error_handler;
    ios->default_rearranger = rearr;
    ios->num_iotasks = num_iotasks;
    ios->num_comptasks = num_comptasks;

    /* Initialize the rearranger options. */
    init_rearr_opts(ios);

    /* Copy the computation communicator into union_comm. */
    if ((mpierr = MPI_Comm_dup(comp_comm, &ios->union_comm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Copy the computation communicator into comp_comm. */
    if ((mpierr = MPI_Comm_dup(comp_comm, &ios->comp_comm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "union_comm = %d comp_comm = %d", ios->union_comm, ios->comp_comm));

    ios->my_comm = ios->comp_comm;
    ustride = stride;

    /* Find MPI rank comp_comm communicator. */
    if ((mpierr = MPI_Comm_rank(ios->comp_comm, &ios->comp_rank)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    if (ios->comp_rank == 0)
        ios->compmaster = MPI_ROOT;
    LOG((2, "comp_rank = %d num_comptasks = %d", ios->comp_rank, ios->num_comptasks));

    /* Create an array that holds the ranks of the tasks to be used
     * for IO. NOTE that sizeof(int) should probably be 1, not
     * sizeof(int) ???*/
    if (!(ios->ioranks = calloc(sizeof(int), ios->num_iotasks)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
    for (int i = 0; i < ios->num_iotasks; i++)
    {
        ios->ioranks[i] = (base + i * ustride) % ios->num_comptasks;
        if (ios->ioranks[i] == ios->comp_rank)
            ios->ioproc = true;
    }
    ios->ioroot = ios->ioranks[0];

    for (int i = 0; i < ios->num_iotasks; i++)
        LOG((3, "ios->ioranks[%d] = %d", i, ios->ioranks[i]));

    /* We are not providing an info object. */
    ios->info = MPI_INFO_NULL;

    /* The task that has an iomaster value of MPI_ROOT will be the
     * root of the IO communicator. */
    if (ios->comp_rank == ios->ioranks[0])
        ios->iomaster = MPI_ROOT;

    /* Create a group for the computation tasks. */
    if ((mpierr = MPI_Comm_group(ios->comp_comm, &ios->compgroup)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Create a group for the IO tasks. */
    if ((mpierr = MPI_Group_incl(ios->compgroup, ios->num_iotasks, ios->ioranks,
                                 &ios->iogroup)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Create an MPI communicator for the IO tasks. */
    if ((mpierr = MPI_Comm_create(ios->comp_comm, ios->iogroup, &ios->io_comm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* For the tasks that are doing IO, get their rank within the IO
     * communicator. For some reason when I check the return value of
     * this MPI call, all tests start to fail! */
    if (ios->ioproc)
    {
        if ((mpierr = MPI_Comm_rank(ios->io_comm, &ios->io_rank)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    }
    else
        ios->io_rank = -1;
    LOG((3, "ios->io_comm = %d ios->io_rank = %d", ios->io_comm, ios->io_rank));

    ios->union_rank = ios->comp_rank;

    /* Add this ios struct to the list in the PIO library. */
    *iosysidp = pio_add_to_iosystem_list(ios);

    /* Allocate buffer space for compute nodes. */
    if ((ret = compute_buffer_init(ios)))
        return ret;

    LOG((2, "Init_Intracomm complete iosysid = %d", *iosysidp));

    return PIO_NOERR;
}

/**
 * Interface to call from pio_init from fortran.
 *
 * @param f90_comp_comm
 * @param num_iotasks the number of IO tasks
 * @param stride the stride to use assigning tasks
 * @param base the starting point when assigning tasks
 * @param rearr the rearranger
 * @param rearr_opts the rearranger options
 * @param iosysidp a pointer that gets the IO system ID
 * @returns 0 for success, error code otherwise
 */
int PIOc_Init_Intracomm_from_F90(int f90_comp_comm,
                                 const int num_iotasks, const int stride,
                                 const int base, const int rearr,
                                 rearr_opt_t *rearr_opts, int *iosysidp)
{
    int ret = PIO_NOERR;
    ret = PIOc_Init_Intracomm(MPI_Comm_f2c(f90_comp_comm), num_iotasks,
                              stride, base, rearr,
                              iosysidp);
    if (ret != PIO_NOERR)
    {
        LOG((1, "PIOc_Init_Intracomm failed"));
        return ret;
    }

    if (rearr_opts)
    {
        LOG((1, "Setting rearranger options, iosys=%d", *iosysidp));
        return PIOc_set_rearr_opts(*iosysidp, rearr_opts->comm_type,
                                   rearr_opts->fcd,
                                   rearr_opts->comm_fc_opts_comp2io.enable_hs,
                                   rearr_opts->comm_fc_opts_comp2io.enable_isend,
                                   rearr_opts->comm_fc_opts_comp2io.max_pend_req,
                                   rearr_opts->comm_fc_opts_io2comp.enable_hs,
                                   rearr_opts->comm_fc_opts_io2comp.enable_isend,
                                   rearr_opts->comm_fc_opts_io2comp.max_pend_req);
    }
    return ret;
}

/**
 * Send a hint to the MPI-IO library.
 *
 * @param iosysid the IO system ID
 * @param hint the hint for MPI
 * @param hintval the value of the hint
 * @returns 0 for success, or PIO_BADID if iosysid can't be found.
 */
int PIOc_set_hint(int iosysid, const char *hint, const char *hintval)
{
    iosystem_desc_t *ios;
    int mpierr; /* Return value for MPI calls. */

    /* Get the iosysid. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* User must provide these. */
    if (!hint || !hintval)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((1, "PIOc_set_hint hint = %s hintval = %s", hint, hintval));

    /* Make sure we have an info object. */
    if (ios->info == MPI_INFO_NULL)
        if ((mpierr = MPI_Info_create(&ios->info)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Set the MPI hint. */
    if (ios->ioproc)
        if ((mpierr = MPI_Info_set(ios->info, hint, hintval)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    return PIO_NOERR;
}

/**
 * Clean up internal data structures, free MPI resources, and exit the
 * pio library.
 *
 * @param iosysid: the io system ID provided by PIOc_Init_Intracomm().
 * @returns 0 for success or non-zero for error.
 * @ingroup PIO_finalize
 */
int PIOc_finalize(int iosysid)
{
    iosystem_desc_t *ios;
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
    int ierr = PIO_NOERR;

    LOG((1, "PIOc_finalize iosysid = %d MPI_COMM_NULL = %d", iosysid,
         MPI_COMM_NULL));

    /* Find the IO system information. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* If asynch IO is in use, send the PIO_MSG_EXIT message from the
     * comp master to the IO processes. This may be called by
     * componets for other components iosysid. So don't send unless
     * there is a valid union_comm. */
    if (ios->async_interface && ios->union_comm != MPI_COMM_NULL)
    {
        int msg = PIO_MSG_EXIT;

        LOG((3, "found iosystem info comproot = %d union_comm = %d comp_idx = %d",
             ios->comproot, ios->union_comm, ios->comp_idx));
        if (!ios->ioproc)
        {
            LOG((2, "sending msg = %d ioroot = %d union_comm = %d", msg,
                 ios->ioroot, ios->union_comm));

            /* Send the message to the message handler. */
            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            /* Send the parameters of the function call. */
            if (!mpierr)
                mpierr = MPI_Bcast((int *)&iosysid, 1, MPI_INT, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        LOG((3, "handling async errors mpierr = %d my_comm = %d", mpierr, ios->my_comm));
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi2(ios, NULL, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
        LOG((3, "async errors bcast"));
    }

    /* Free this memory that was allocated in init_intracomm. */
    if (ios->ioranks)
        free(ios->ioranks);
    LOG((3, "Freed ioranks."));

    /* Free the buffer pool. */
    int niosysid;
    if ((ierr = pio_num_iosystem(&niosysid)))
        return pio_err(ios, NULL, ierr, __FILE__, __LINE__);
    LOG((2, "%d iosystems are still open.", niosysid));

    /* Only free the buffer pool if this is the last open iosysid. */
    if (niosysid == 1)
    {
        free_cn_buffer_pool(ios);
        LOG((2, "Freed buffer pool."));
    }

    /* Free the MPI groups. */
    if (ios->compgroup != MPI_GROUP_NULL)
        MPI_Group_free(&ios->compgroup);

    if (ios->iogroup != MPI_GROUP_NULL)
        MPI_Group_free(&(ios->iogroup));

    /* Free the MPI communicators. my_comm is just a copy (but not an
     * MPI copy), so does not have to have an MPI_Comm_free()
     * call. comp_comm and io_comm are MPI duplicates of the comms
     * handed into init_intercomm. So they need to be freed by MPI. */
    if (ios->intercomm != MPI_COMM_NULL)
        MPI_Comm_free(&ios->intercomm);
    if (ios->union_comm != MPI_COMM_NULL)
        MPI_Comm_free(&ios->union_comm);
    if (ios->io_comm != MPI_COMM_NULL)
        MPI_Comm_free(&ios->io_comm);
    if (ios->comp_comm != MPI_COMM_NULL)
        MPI_Comm_free(&ios->comp_comm);
    if (ios->my_comm != MPI_COMM_NULL)
        ios->my_comm = MPI_COMM_NULL;

    /* Free the MPI Info object. */
    if (ios->info != MPI_INFO_NULL)
        MPI_Info_free(&ios->info);

    /* Delete the iosystem_desc_t data associated with this id. */
    LOG((2, "About to delete iosysid %d.", iosysid));
    if ((ierr = pio_delete_iosystem_from_list(iosysid)))
        return pio_err(NULL, NULL, ierr, __FILE__, __LINE__);

    pio_finalize_logging();

    LOG((2, "PIOc_finalize completed successfully"));
    return PIO_NOERR;
}

/**
 * Return a logical indicating whether this task is an IO task.
 *
 * @param iosysid the io system ID
 * @param ioproc a pointer that gets 1 if task is an IO task, 0
 * otherwise. Ignored if NULL.
 * @returns 0 for success, or PIO_BADID if iosysid can't be found.
 */
int PIOc_iam_iotask(int iosysid, bool *ioproc)
{
    iosystem_desc_t *ios;

    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    if (ioproc)
        *ioproc = ios->ioproc;

    return PIO_NOERR;
}

/**
 * Return the rank of this task in the IO communicator or -1 if this
 * task is not in the communicator.
 *
 * @param iosysid the io system ID
 * @param iorank a pointer that gets the io rank, or -1 if task is not
 * in the IO communicator. Ignored if NULL.
 * @returns 0 for success, or PIO_BADID if iosysid can't be found.
 */
int PIOc_iotask_rank(int iosysid, int *iorank)
{
    iosystem_desc_t *ios;

    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    if (iorank)
        *iorank = ios->io_rank;

    return PIO_NOERR;
}

/**
 * Return true if this iotype is supported in the build, 0 otherwise.
 *
 * @param iotype the io type to check
 * @returns 1 if iotype is in build, 0 if not.
 */
int PIOc_iotype_available(int iotype)
{
    switch(iotype)
    {
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
    case PIO_IOTYPE_NETCDF4C:
        return 1;
#endif
    case PIO_IOTYPE_NETCDF:
        return 1;
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        return 1;
        break;
#endif
    default:
        return 0;
    }
}

/**
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
 * @param world the communicator containing all the available tasks.
 *
 * @param num_io_procs the number of processes for the IO component.
 *
 * @param io_proc_list an array of lenth num_io_procs with the
 * processor number for each IO processor. If NULL then the IO
 * processes are assigned starting at processes 0.
 *
 * @param component_count number of computational components
 *
 * @param num_procs_per_comp an array of int, of length
 * component_count, with the number of processors in each computation
 * component.
 *
 * @param proc_list an array of arrays containing the processor
 * numbers for each computation component. If NULL then the
 * computation components are assigned processors sequentially
 * starting with processor num_io_procs.
 *
 * @param user_io_comm pointer to an MPI_Comm. If not NULL, it will
 * get an MPI duplicate of the IO communicator. (It is a full
 * duplicate and later must be freed with MPI_Free() by the caller.)
 *
 * @param user_comp_comm pointer to an array of pointers to MPI_Comm;
 * the array is of length component_count. If not NULL, it will get an
 * MPI duplicate of each computation communicator. (These are full
 * duplicates and each must later be freed with MPI_Free() by the
 * caller.)
 *
 * @param iosysidp pointer to array of length component_count that
 * gets the iosysid for each component.
 *
 * @return PIO_NOERR on success, error code otherwise.
 * @ingroup PIO_init
 */
int PIOc_Init_Async(MPI_Comm world, int num_io_procs, int *io_proc_list,
                    int component_count, int *num_procs_per_comp, int **proc_list,
                    MPI_Comm *user_io_comm, MPI_Comm *user_comp_comm, int *iosysidp)
{
    int my_rank;
    int **my_proc_list;
    int *my_io_proc_list;
    int mpierr;
    int ret;

    /* Check input parameters. */
    if (num_io_procs < 1 || component_count < 1 || !num_procs_per_comp || !iosysidp)
        return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* Temporarily limit to one computational component. */
    if (component_count > 1)
        return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* Turn on the logging system for PIO. */
    pio_init_logging();
    LOG((1, "PIOc_Init_Async component_count = %d", component_count));

    /* If the user did not supply a list of process numbers to use for
     * IO, create it. */
    if (!io_proc_list)
    {
        if (!(my_io_proc_list = malloc(num_io_procs * sizeof(int))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        for (int p = 0; p < num_io_procs; p++)
            my_io_proc_list[p] = p;
    }
    else
        my_io_proc_list = io_proc_list;

    /* If the user did not provide a list of processes for each
     * component, create one. */
    if (!proc_list)
    {
        int last_proc = 0;

        /* Allocate space for array of arrays. */
        if (!(my_proc_list = malloc((component_count + 1) * sizeof(int *))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

        /* Fill the array of arrays. */
        for (int cmp = 0; cmp < component_count + 1; cmp++)
        {
            LOG((3, "calculating processors for component %d", cmp));

            /* Allocate space for each array. */
            if (!(my_proc_list[cmp] = malloc(num_procs_per_comp[cmp] * sizeof(int))))
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

            int proc;
            for (proc = last_proc; proc < num_procs_per_comp[cmp] + last_proc; proc++)
            {
                my_proc_list[cmp][proc - last_proc] = proc;
                LOG((3, "my_proc_list[%d][%d] = %d", cmp, proc - last_proc, proc));
            }
            last_proc = proc;
        }
    }
    else
        my_proc_list = proc_list;

    /* Get rank of this task. */
    if ((ret = MPI_Comm_rank(world, &my_rank)))
        return check_mpi(NULL, ret, __FILE__, __LINE__);

    /* Is this process in the IO component? */
    int pidx;
    for (pidx = 0; pidx < num_procs_per_comp[0]; pidx++)
        if (my_rank == my_proc_list[0][pidx])
            break;
    int in_io = (pidx == num_procs_per_comp[0]) ? 0 : 1;
    LOG((3, "in_io = %d", in_io));

    /* Allocate struct to hold io system info for each computation component. */
    /* Allocate struct to hold io system info for each component. */
    iosystem_desc_t *iosys[component_count], *my_iosys;
    for (int cmp1 = 0; cmp1 < component_count; cmp1++)
        if (!(iosys[cmp1] = (iosystem_desc_t *)calloc(1, sizeof(iosystem_desc_t))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Create group for world. */
    MPI_Group world_group;
    if ((ret = MPI_Comm_group(world, &world_group)))
        return check_mpi(NULL, ret, __FILE__, __LINE__);
    LOG((3, "world group created\n"));

    /* We will create a group for the IO component. */
    MPI_Group io_group;

    /* The shared IO communicator. */
    MPI_Comm io_comm;

    /* Rank of current process in IO communicator. */
    int io_rank = -1;

    /* Set to MPI_ROOT on master process, MPI_PROC_NULL on other
     * processes. */
    int iomaster;

    /* Create a group for the IO component. */
    if ((ret = MPI_Group_incl(world_group, num_io_procs, my_io_proc_list, &io_group)))
        return check_mpi(NULL, ret, __FILE__, __LINE__);
    LOG((3, "created IO group - io_group = %d group empty is %d", io_group, MPI_GROUP_EMPTY));
    for (int p = 0; p < num_io_procs; p++)
        LOG((3, "my_io_proc_list[%d] = %d", p, my_io_proc_list[p]));

    /* There is one shared IO comm. Create it. */
    if ((ret = MPI_Comm_create(world, io_group, &io_comm)))
        return check_mpi(NULL, ret, __FILE__, __LINE__);
    LOG((3, "created io comm io_comm = %d", io_comm));

    /* Does the user want a copy of the IO communicator? */
    if (user_io_comm)
    {
        *user_io_comm = MPI_COMM_NULL;
        if (in_io)
            if ((mpierr = MPI_Comm_dup(io_comm, user_io_comm)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    }

    /* For processes in the IO component, get their rank within the IO
     * communicator. */
    if (in_io)
    {
        LOG((3, "about to get io rank"));
        if ((ret = MPI_Comm_rank(io_comm, &io_rank)))
            return check_mpi(NULL, ret, __FILE__, __LINE__);
        iomaster = !io_rank ? MPI_ROOT : MPI_PROC_NULL;
        LOG((3, "intracomm created for io_comm = %d io_rank = %d IO %s",
             io_comm, io_rank, iomaster == MPI_ROOT ? "MASTER" : "SERVANT"));
    }

    /* We will create a group for each component. */
    MPI_Group group[component_count + 1];

    /* We will also create a group for each component and the IO
     * component processes (i.e. a union of computation and IO
     * processes. */
    MPI_Group union_group[component_count];

    /* For each component, starting with the IO component. */
    for (int cmp = 0; cmp < component_count + 1; cmp++)
    {
        LOG((3, "processing component %d", cmp));

        /* Don't start initing iosys until after IO component. */
        if (cmp)
        {
            /* Get pointer to current iosys. */
            my_iosys = iosys[cmp - 1];

            /* Initialize some values. */
            my_iosys->io_comm = MPI_COMM_NULL;
            my_iosys->comp_comm = MPI_COMM_NULL;
            my_iosys->union_comm = MPI_COMM_NULL;
            my_iosys->intercomm = MPI_COMM_NULL;
            my_iosys->my_comm = MPI_COMM_NULL;
            my_iosys->async_interface = 1;
            my_iosys->error_handler = default_error_handler;
            my_iosys->num_comptasks = num_procs_per_comp[cmp];
            my_iosys->num_iotasks = num_procs_per_comp[0];
            my_iosys->compgroup = MPI_GROUP_NULL;
            my_iosys->iogroup = MPI_GROUP_NULL;

            /* The rank of the computation leader in the union comm. */
            my_iosys->comproot = num_procs_per_comp[0];
            LOG((3, "my_iosys->comproot = %d", my_iosys->comproot));

            /* We are not providing an info object. */
            my_iosys->info = MPI_INFO_NULL;
        }

        /* Create a group for this component. */
        if ((ret = MPI_Group_incl(world_group, num_procs_per_comp[cmp], my_proc_list[cmp],
                                  &group[cmp])))
            return check_mpi(NULL, ret, __FILE__, __LINE__);
        LOG((3, "created component MPI group - group[%d] = %d", cmp, group[cmp]));

        /* For all the computation components (i.e. cmp != 0), create
         * a union group with their processors and the processors of
         * the (shared) IO component. */
        if (cmp)
        {
            /* How many processors in the union comm? */
            int nprocs_union = num_procs_per_comp[0] + num_procs_per_comp[cmp];

            /* This will hold proc numbers from both computation and IO
             * components. */
            int proc_list_union[nprocs_union];

            /* Add proc numbers from IO. */
            for (int p = 0; p < num_procs_per_comp[0]; p++)
                proc_list_union[p] = my_proc_list[0][p];

            /* Add proc numbers from computation component. */
            for (int p = 0; p < num_procs_per_comp[cmp]; p++)
                proc_list_union[p + num_procs_per_comp[0]] = my_proc_list[cmp][p];

            /* Create the union group. */
            if ((ret = MPI_Group_incl(world_group, nprocs_union, proc_list_union,
                                      &union_group[cmp - 1])))
                return check_mpi(NULL, ret, __FILE__, __LINE__);
            LOG((3, "created union MPI_group - union_group[%d] = %d with %d procs", cmp, union_group[cmp-1], nprocs_union));
        }

        /* Remember whether this process is in the IO component. */
        if (cmp)
            my_iosys->ioproc = in_io;

        /* Is this process in this computation component (which is the
         * IO component if cmp == 0)? */
        int in_cmp = 0;
        for (pidx = 0; pidx < num_procs_per_comp[cmp]; pidx++)
            if (my_rank == my_proc_list[cmp][pidx])
                break;
        in_cmp = (pidx == num_procs_per_comp[cmp]) ? 0 : 1;
        LOG((3, "pidx = %d num_procs_per_comp[%d] = %d in_cmp = %d",
             pidx, cmp, num_procs_per_comp[cmp], in_cmp));

        /* Create an intracomm for this component. Only processes in
         * the component need to participate in the intracomm create
         * call. */
        /* Create the intracomm from the group. */
        LOG((3, "creating intracomm cmp = %d from group[%d] = %d", cmp, cmp, group[cmp]));

        /* We handle the IO comm differently (cmp == 0). */
        if (!cmp)
        {
            /* LOG((3, "about to create io comm")); */
            /* if ((ret = MPI_Comm_create_group(world, group[cmp], cmp, &io_comm))) */
            /*     return check_mpi(NULL, ret, __FILE__, __LINE__); */
            /* LOG((3, "about to get io rank"));                 */
            /* if ((ret = MPI_Comm_rank(io_comm, &io_rank))) */
            /*     return check_mpi(NULL, ret, __FILE__, __LINE__); */
            /* iomaster = !io_rank ? MPI_ROOT : MPI_PROC_NULL; */
            /* LOG((3, "intracomm created for cmp = %d io_comm = %d io_rank = %d IO %s", */
            /*      cmp, io_comm, io_rank, iomaster == MPI_ROOT ? "MASTER" : "SERVANT")); */
        }
        else
        {
            if ((ret = MPI_Comm_create(world, group[cmp], &my_iosys->comp_comm)))
                return check_mpi(NULL, ret, __FILE__, __LINE__);

            if (in_cmp)
            {
                /* Does the user want a copy? */
                if (user_comp_comm)
                    if ((mpierr = MPI_Comm_dup(my_iosys->comp_comm, &user_comp_comm[cmp - 1])))
                        return check_mpi(NULL, mpierr, __FILE__, __LINE__);

                /* Get the rank in this comp comm. */
                if ((ret = MPI_Comm_rank(my_iosys->comp_comm, &my_iosys->comp_rank)))
                    return check_mpi(NULL, ret, __FILE__, __LINE__);

                /* Set comp_rank 0 to be the compmaster. It will have
                 * a setting of MPI_ROOT, all other tasks will have a
                 * setting of MPI_PROC_NULL. */
                my_iosys->compmaster = my_iosys->comp_rank ? MPI_PROC_NULL : MPI_ROOT;

                LOG((3, "intracomm created for cmp = %d comp_comm = %d comp_rank = %d comp %s",
                     cmp, my_iosys->comp_comm, my_iosys->comp_rank,
                     my_iosys->compmaster == MPI_ROOT ? "MASTER" : "SERVANT"));
            }
        }


        /* If this is the IO component, make a copy of the IO comm for
         * each computational component. */
        if (in_io)
            if (cmp)
            {
                LOG((3, "making a dup of io_comm = %d io_rank = %d", io_comm, io_rank));
                if ((ret = MPI_Comm_dup(io_comm, &my_iosys->io_comm)))
                    return check_mpi(NULL, ret, __FILE__, __LINE__);
                LOG((3, "dup of io_comm = %d io_rank = %d", my_iosys->io_comm, io_rank));
                my_iosys->iomaster = iomaster;
                my_iosys->io_rank = io_rank;
                my_iosys->ioroot = 0;
                my_iosys->comp_idx = cmp - 1;
            }

        /* All the processes in this component, and the IO component,
         * are part of the union_comm. */
        if (cmp)
        {
            if (in_io || in_cmp)
            {
                LOG((3, "my_iosys->io_comm = %d group = %d", my_iosys->io_comm, union_group[cmp-1]));
                /* Create a group for the union of the IO component
                 * and one of the computation components. */
                if ((ret = MPI_Comm_create(world, union_group[cmp - 1],
                                           &my_iosys->union_comm)))
                    return check_mpi(NULL, ret, __FILE__, __LINE__);

                if ((ret = MPI_Comm_rank(my_iosys->union_comm, &my_iosys->union_rank)))
                    return check_mpi(NULL, ret, __FILE__, __LINE__);

                /* Set my_comm to union_comm for async. */
                my_iosys->my_comm = my_iosys->union_comm;
                LOG((3, "intracomm created for union cmp = %d union_rank = %d union_comm = %d",
                     cmp, my_iosys->union_rank, my_iosys->union_comm));

                if (in_io)
                {
                    LOG((3, "my_iosys->io_comm = %d", my_iosys->io_comm));
                    /* Create the intercomm from IO to computation component. */
                    LOG((3, "about to create intercomm for IO component to cmp = %d "
                         "my_iosys->io_comm = %d", cmp, my_iosys->io_comm));
                    if ((ret = MPI_Intercomm_create(my_iosys->io_comm, 0, my_iosys->union_comm,
                                                    my_proc_list[cmp][0], 0, &my_iosys->intercomm)))
                        return check_mpi(NULL, ret, __FILE__, __LINE__);
                }
                else
                {
                    /* Create the intercomm from computation component to IO component. */
                    LOG((3, "about to create intercomm for cmp = %d my_iosys->comp_comm = %d", cmp,
                         my_iosys->comp_comm));
                    if ((ret = MPI_Intercomm_create(my_iosys->comp_comm, 0, my_iosys->union_comm,
                                                    my_proc_list[0][0], 0, &my_iosys->intercomm)))
                        return check_mpi(NULL, ret, __FILE__, __LINE__);
                }
                LOG((3, "intercomm created for cmp = %d", cmp));
            }

            /* Add this id to the list of PIO iosystem ids. */
            iosysidp[cmp - 1] = pio_add_to_iosystem_list(my_iosys);
            LOG((2, "new iosys ID added to iosystem_list iosysid = %d\n", iosysidp[cmp - 1]));
        }
    }

    /* Now call the function from which the IO tasks will not return
     * until the PIO_MSG_EXIT message is sent. This will handle all
     * components. */
    if (in_io)
    {
        LOG((2, "Starting message handler io_rank = %d component_count = %d",
             io_rank, component_count));
        if ((ret = pio_msg_handler2(io_rank, component_count, iosys, io_comm)))
            return pio_err(NULL, NULL, ret, __FILE__, __LINE__);
        LOG((2, "Returned from pio_msg_handler2() ret = %d", ret));
    }

    /* Free resources if needed. */
    LOG((2, "PIOc_Init_Async starting to free resources"));
    if (!io_proc_list)
        free(my_io_proc_list);

    if (in_io)
        if ((mpierr = MPI_Comm_free(&io_comm)))
            return check_mpi(NULL, ret, __FILE__, __LINE__);

    if (!proc_list)
    {
        for (int cmp = 0; cmp < component_count + 1; cmp++)
            free(my_proc_list[cmp]);
        free(my_proc_list);
    }

    /* Free MPI groups. */
    if ((ret = MPI_Group_free(&io_group)))
        return check_mpi(NULL, ret, __FILE__, __LINE__);

    for (int cmp = 0; cmp < component_count + 1; cmp++)
    {
        if ((ret = MPI_Group_free(&group[cmp])))
            return check_mpi(NULL, ret, __FILE__, __LINE__);
        if (cmp)
            if ((ret = MPI_Group_free(&union_group[cmp - 1])))
                return check_mpi(NULL, ret, __FILE__, __LINE__);
    }

    if ((ret = MPI_Group_free(&world_group)))
        return check_mpi(NULL, ret, __FILE__, __LINE__);

    LOG((2, "successfully done with PIO_Init_Async"));
    return PIO_NOERR;
}

/**
 * Set the target blocksize for the box rearranger.
 *
 * @param newblocksize the new blocksize.
 * @returns 0 for success.
 * @ingroup PIO_set_blocksize
 */
int PIOc_set_blocksize(int newblocksize)
{
    if (newblocksize > 0)
        blocksize = newblocksize;
    return PIO_NOERR;
}
