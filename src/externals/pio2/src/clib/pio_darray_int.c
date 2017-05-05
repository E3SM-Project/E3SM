/** @file
 *
 * Private functions to help read and write distributed arrays in PIO.
 *
 * When arrays are distributed, each processor holds some of the
 * array. Only by combining the distributed arrays from all processor
 * can the full array be obtained.
 *
 * @author Jim Edwards
 */

#include <config.h>
#include <pio.h>
#include <pio_internal.h>

/* 10MB default limit. */
extern PIO_Offset pio_buffer_size_limit;

/* Initial size of compute buffer. */
bufsize pio_cnbuffer_limit = 33554432;

/* Global buffer pool pointer. */
extern void *CN_bpool;

/* Maximum buffer usage. */
extern PIO_Offset maxusage;

/**
 * Initialize the compute buffer to size pio_cnbuffer_limit.
 *
 * This routine initializes the compute buffer pool if the bget memory
 * management is used. If malloc is used (that is, PIO_USE_MALLOC is
 * non zero), this function does nothing.
 *
 * @param ios pointer to the iosystem descriptor which will use the
 * new buffer.
 * @returns 0 for success, error code otherwise.
 */
int compute_buffer_init(iosystem_desc_t *ios)
{
#if !PIO_USE_MALLOC

    if (!CN_bpool)
    {
        if (!(CN_bpool = malloc(pio_cnbuffer_limit)))
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

        bpool(CN_bpool, pio_cnbuffer_limit);
        if (!CN_bpool)
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

        bectl(NULL, malloc, free, pio_cnbuffer_limit);
    }
#endif
    LOG((2, "compute_buffer_init CN_bpool = %d", CN_bpool));

    return PIO_NOERR;
}

/**
 * Write a set of one or more aggregated arrays to output file. This
 * function is only used with parallel-netcdf and netcdf-4 parallel
 * iotypes. Serial io types use pio_write_darray_multi_nc_serial().
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param nvars the number of variables to be written with this
 * decomposition.
 * @param vid: an array of the variable ids to be written.
 * @param iodesc_ndims: the number of dimensions explicitly in the
 * iodesc.
 * @param basetype the basic type of the minimal data unit
 * @param maxregions max number of blocks to be written from
 * this iotask.
 * @param firstregion pointer to the first element of a linked
 * list of region descriptions. May be NULL.
 * @param llen length of the iobuffer on this task for a single
 * field.
 * @param num_aiotasks actual number of iotasks participating.
 * @param iobuf the buffer to be written from this mpi task. May be
 * null. for example we have 8 ionodes and a distributed array with
 * global size 4, then at least 4 nodes will have a null iobuf. In
 * practice the box rearranger trys to have at least blocksize bytes
 * on each io task and so if the total number of bytes to write is
 * less than blocksize*numiotasks then some iotasks will have a NULL
 * iobuf.
 * @param frame the frame or record dimension for each of the nvars
 * variables in iobuf.
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int pio_write_darray_multi_nc(file_desc_t *file, int nvars, const int *vid, int iodesc_ndims,
                              MPI_Datatype basetype, int maxregions, io_region *firstregion,
                              PIO_Offset llen, int num_aiotasks, void *iobuf,
                              const int *frame)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;     /* Pointer to var info struct. */
    int fndims;            /* Number of dims for this var in the file. */
    int dsize;             /* Data size (for one region). */
    int tsize;             /* Size of MPI type. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
    int ierr = PIO_NOERR;

    /* Check inputs. */
    pioassert(file && file->iosystem && vid && vid[0] >= 0 && vid[0] <= PIO_MAX_VARS,
              "invalid input", __FILE__, __LINE__);

    LOG((1, "pio_write_darray_multi_nc nvars = %d iodesc_ndims = %d basetype = %d "
         "maxregions = %d llen = %d num_aiotasks = %d", nvars, iodesc_ndims,
         basetype, maxregions, llen, num_aiotasks));

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:write_darray_multi_nc");
#endif

    /* Get pointer to iosystem. */
    ios = file->iosystem;

    /* Point to var description scruct for first var. */
    vdesc = file->varlist + vid[0];

    /* If async is in use, send message to IO master task. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = 0;
            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->pio_ncid, 1, MPI_INT, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* Find out how many dims this variable has. */
    if ((ierr = PIOc_inq_varndims(file->pio_ncid, vid[0], &fndims)))
        return pio_err(ios, file, ierr, __FILE__, __LINE__);

    /* Find out the size of the MPI type. */
    if ((mpierr = MPI_Type_size(basetype, &tsize)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);
    LOG((2, "fndims = %d tsize = %d", fndims, tsize));

    /* If this is an IO task write the data. */
    if (ios->ioproc)
    {
        io_region *region = firstregion;
        int rrcnt = 0; /* Number of subarray requests (pnetcdf only). */
        void *bufptr;
        size_t start[fndims];
        size_t count[fndims];
        int ndims = iodesc_ndims;
        PIO_Offset *startlist[maxregions]; /* Array of start arrays for ncmpi_iput_varn(). */
        PIO_Offset *countlist[maxregions]; /* Array of count  arrays for ncmpi_iput_varn(). */

        LOG((3, "maxregions = %d", maxregions));

        /* Process each region of data to be written. */
        for (int regioncnt = 0; regioncnt < maxregions; regioncnt++)
        {
            /* Init start/count arrays to zero. */
            for (int i = 0; i < fndims; i++)
            {
                start[i] = 0;
                count[i] = 0;
            }

            if (region)
            {
                if (vdesc->record >= 0)
                {
                    /* This is a record based multidimensional
                     * array. Figure out start/count for all but the
                     * record dimension (dimid 0). */
                    for (int i = fndims - ndims; i < fndims; i++)
                    {
                        start[i] = region->start[i - (fndims - ndims)];
                        count[i] = region->count[i - (fndims - ndims)];
                    }

                    /* Now figure out start/count for record dimension. */
                    if (fndims > 1 && ndims < fndims && count[1] > 0)
                    {
                        count[0] = 1;
                        start[0] = frame[0];
                    }
                    else if (fndims == ndims)
                    {
                        /* ??? */
                        start[0] += vdesc->record;
                    }
                }
                else
                {
                    /* This is a non record variable. */
                    for (int i = 0; i < ndims; i++)
                    {
                        start[i] = region->start[i];
                        count[i] = region->count[i];
                    }
                }
#if PIO_ENABLE_LOGGING
                /* Log arrays for debug purposes. */
                for (int i = 0; i < ndims; i++)
                    LOG((3, "start[%d] = %d count[%d] = %d", i, start[i], i, count[i]));
#endif /* PIO_ENABLE_LOGGING */
            }

            /* IO tasks will run the netCDF/pnetcdf functions to write the data. */
            switch (file->iotype)
            {
#ifdef _NETCDF4
            case PIO_IOTYPE_NETCDF4P:
                /* For each variable to be written. */
                for (int nv = 0; nv < nvars; nv++)
                {
                    /* Set the start of the record dimension. (Hasn't
                     * this already been set above ???) */
                    if (vdesc->record >= 0 && ndims < fndims)
                        start[0] = frame[nv];

                    /* If there is data for this region, get a pointer to it. */
                    if (region)
                        bufptr = (void *)((char *)iobuf + tsize * (nv * llen + region->loffset));

                    /* Ensure collective access. */
                    ierr = nc_var_par_access(file->fh, vid[nv], NC_COLLECTIVE);

                    /* Write the data for this variable. */
                    if (!ierr)
                        ierr = nc_put_vara(file->fh, vid[nv], (size_t *)start, (size_t *)count, bufptr);
                }
                break;
#endif
#ifdef _PNETCDF
            case PIO_IOTYPE_PNETCDF:
                /* Get the total number of data elements we are
                 * writing for this region. */
                dsize = 1;
                for (int i = 0; i < fndims; i++)
                    dsize *= count[i];
                LOG((3, "dsize = %d", dsize));

                /* For pnetcdf's ncmpi_iput_varn() function, we need
                 * to provide arrays of arrays for start/count. */
                if (dsize > 0)
                {
                    /* Allocate storage for start/count arrays for
                     * this region. */
                    if (!(startlist[rrcnt] = calloc(fndims, sizeof(PIO_Offset))))
                        return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);
                    if (!(countlist[rrcnt] = calloc(fndims, sizeof(PIO_Offset))))
                        return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);

                    /* Copy the start/count arrays for this region. */
                    for (int i = 0; i < fndims; i++)
                    {
                        startlist[rrcnt][i] = start[i];
                        countlist[rrcnt][i] = count[i];
                        LOG((3, "startlist[%d][%d] = %d countlist[%d][%d] = %d", rrcnt, i,
                             startlist[rrcnt][i], rrcnt, i, countlist[rrcnt][i]));
                    }
                    rrcnt++;
                }

                /* Do this when we reach the last region. */
                if (regioncnt == maxregions - 1)
                {
                    /* For each variable to be written. */
                    for (int nv = 0; nv < nvars; nv++)
                    {
                        /* Get the var info. */
                        vdesc = file->varlist + vid[nv];

                        /* If this is a record var, set the start for
                         * the record dimension. */
                        if (vdesc->record >= 0 && ndims < fndims)
                            for (int rc = 0; rc < rrcnt; rc++)
                                startlist[rc][0] = frame[nv];

                        /* Get a pointer to the data. */
                        bufptr = (void *)((char *)iobuf + nv * tsize * llen);

                        /* ??? */
                        int reqn = 0;
                        if (vdesc->nreqs % PIO_REQUEST_ALLOC_CHUNK == 0 )
                        {
                            if (!(vdesc->request = realloc(vdesc->request, sizeof(int) *
                                                           (vdesc->nreqs + PIO_REQUEST_ALLOC_CHUNK))))
                                return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);

                            for (int i = vdesc->nreqs; i < vdesc->nreqs + PIO_REQUEST_ALLOC_CHUNK; i++)
                                vdesc->request[i] = NC_REQ_NULL;
                            reqn = vdesc->nreqs;
                        }
                        else
                            while(vdesc->request[reqn] != NC_REQ_NULL)
                                reqn++;

                        /* Write, in non-blocking fashion, a list of subarrays. */
                        LOG((3, "about to call ncmpi_iput_varn() vid[%d] = %d rrcnt = %d, llen = %d",
                             nv, vid[nv], rrcnt, llen));
                        ierr = ncmpi_iput_varn(file->fh, vid[nv], rrcnt, startlist, countlist,
                                               bufptr, llen, basetype, vdesc->request + reqn);

                        /* keeps wait calls in sync */
                        if (vdesc->request[reqn] == NC_REQ_NULL)
                            vdesc->request[reqn] = PIO_REQ_NULL;

                        vdesc->nreqs += reqn + 1;

                    }

                    /* Free resources. */
                    for (int i = 0; i < rrcnt; i++)
                    {
                        free(startlist[i]);
                        free(countlist[i]);
                    }
                }
                break;
#endif
            default:
                return pio_err(ios, file, PIO_EBADIOTYPE, __FILE__, __LINE__);
            }

            /* Go to next region. */
            if (region)
                region = region->next;
        } /* next regioncnt */
    } /* endif (ios->ioproc) */

    /* Check the return code from the netCDF/pnetcdf call. */
    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:write_darray_multi_nc");
#endif

    return ierr;
}

/**
 * Write a set of one or more aggregated arrays to output file in
 * serial mode. This function is called for netCDF classic and
 * netCDF-4 serial iotypes. Parallel iotypes use
 * pio_write_darray_multi_nc().
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to.
 * @param nvars the number of variables to be written with this
 * decomposition.
 * @param vid an array of the variable ids to be written
 * @param iodesc_ndims the number of dimensions explicitly in the
 * iodesc.
 * @param basetype the basic type of the minimal data unit
 * @param maxregions max number of blocks to be written from this
 * iotask.
 * @param firstregion pointer to the first element of a linked
 * list of region descriptions. May be NULL.
 * @param llen length of the iobuffer on this task for a single
 * field.
 * @param num_aiotasks actual number of iotasks participating
 * @param iobuf the buffer to be written from this mpi task. May be
 * null. for example we have 8 ionodes and a distributed array with
 * global size 4, then at least 4 nodes will have a null iobuf. In
 * practice the box rearranger trys to have at least blocksize bytes
 * on each io task and so if the total number of bytes to write is
 * less than blocksize*numiotasks then some iotasks will have a NULL
 * iobuf.
 * @param frame the record dimension for each of the nvars variables
 * in iobuf.
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int pio_write_darray_multi_nc_serial(file_desc_t *file, int nvars, const int *vid, int iodesc_ndims,
                                     MPI_Datatype basetype, int maxregions, io_region *firstregion,
                                     PIO_Offset llen, int num_aiotasks, void *iobuf,
                                     const int *frame)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;     /* Contains info about the variable. */
    int dsize;             /* Size in bytes of one element of basetype. */
    int fndims;            /* Number of dims in the var in the file. */
    int tsize;             /* Size of the MPI type, in bytes. */
    MPI_Status status;     /* Recv status for MPI. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
    int ierr;              /* Return code. */

    /* Check inputs. */
    pioassert(file && file->iosystem && vid && vid[0] >= 0 && vid[0] <= PIO_MAX_VARS,
              "invalid input", __FILE__, __LINE__);

    LOG((1, "pio_write_darray_multi_nc_serial nvars = %d iodesc_ndims = %d basetype = %d "
         "maxregions = %d llen = %d num_aiotasks = %d", nvars, iodesc_ndims,
         basetype, maxregions, llen, num_aiotasks));

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:write_darray_multi_nc_serial");
#endif

    /* Get the iosystem info. */
    ios = file->iosystem;

    /* Get the var info. */
    vdesc = file->varlist + vid[0];

    LOG((2, "vdesc record %d ndims %d nreqs %d ios->async_interface = %d", vdesc->record,
         vdesc->ndims, vdesc->nreqs, ios->async_interface));

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = 0;

            if (ios->comp_rank == 0)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&file->pio_ncid, 1, MPI_INT, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* Get the number of dimensions. */
    if ((ierr = PIOc_inq_varndims(file->pio_ncid, vid[0], &fndims)))
        return pio_err(ios, file, ierr, __FILE__, __LINE__);

    /* Get the size of the type. */
    if ((mpierr = MPI_Type_size(basetype, &tsize)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);
    LOG((3, "fndims = %d tsize = %d", fndims, tsize));

    /* Only IO tasks participate in this code. */
    if (ios->ioproc)
    {
        io_region *region;
        void *bufptr;
        size_t tmp_start[fndims * maxregions]; /* A start array for each region. */
        size_t tmp_count[fndims * maxregions]; /* A count array for each region. */

        LOG((3, "maxregions = %d", maxregions));

        /* Fill the tmp_start and tmp_count arrays, which contain the
         * start and count arrays for all regions. */
        region = firstregion;
        for (int regioncnt = 0; regioncnt < maxregions; regioncnt++)
        {
            /* Initialize the start/count arrays for this region to 0. */
            for (int i = 0; i < fndims; i++)
            {
                tmp_start[i + regioncnt * fndims] = 0;
                tmp_count[i + regioncnt * fndims] = 0;
            }

            if (region)
            {
                if (vdesc->record >= 0)
                {
                    /* This is a record based multidimensional
                     * array. Copy start/count for non-record
                     * dimensions. */
                    for (int i = fndims - iodesc_ndims; i < fndims; i++)
                    {
                        tmp_start[i + regioncnt * fndims] = region->start[i - (fndims - iodesc_ndims)];
                        tmp_count[i + regioncnt * fndims] = region->count[i - (fndims - iodesc_ndims)];
                        LOG((3, "tmp_start[%d] = %d tmp_count[%d] = %d", i + regioncnt * fndims,
                             tmp_start[i + regioncnt * fndims], i + regioncnt * fndims,
                             tmp_count[i + regioncnt * fndims]));
                    }
                }
                else
                {
                    /* This is not a record based multidimensional array. */
                    for (int i = 0; i < iodesc_ndims; i++)
                    {
                        tmp_start[i + regioncnt * fndims] = region->start[i];
                        tmp_count[i + regioncnt * fndims] = region->count[i];
                        LOG((3, "tmp_start[%d] = %d tmp_count[%d] = %d", i + regioncnt * fndims,
                             tmp_start[i + regioncnt * fndims], i + regioncnt * fndims,
                             tmp_count[i + regioncnt * fndims]));
                    }
                }

                /* Move to next region. */
                region = region->next;

            } /* endif region */
        } /* next regioncnt */

        /* Tasks other than 0 will send their data to task 0. */
        if (ios->io_rank > 0)
        {
            /* Do a handshake. */
            if ((mpierr = MPI_Recv(&ierr, 1, MPI_INT, 0, 0, ios->io_comm, &status)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);

            /* Send local length of iobuffer for each field (all
             * fields are the same length). */
            if ((mpierr = MPI_Send((void *)&llen, 1, MPI_OFFSET, 0, ios->io_rank, ios->io_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
            LOG((3, "sent llen = %d", llen));

            /* Send the number of data regions, the start/count for
             * all regions, and the data buffer with all the data. */
            if (llen > 0)
            {
                if ((mpierr = MPI_Send((void *)&maxregions, 1, MPI_INT, 0, ios->io_rank + ios->num_iotasks, ios->io_comm)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(tmp_start, maxregions * fndims, MPI_OFFSET, 0, ios->io_rank + 2 * ios->num_iotasks,
                                       ios->io_comm)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(tmp_count, maxregions * fndims, MPI_OFFSET, 0, ios->io_rank + 3 * ios->num_iotasks,
                                       ios->io_comm)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(iobuf, nvars * llen, basetype, 0, ios->io_rank + 4 * ios->num_iotasks, ios->io_comm)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                LOG((3, "sent data for maxregions = %d", maxregions));
            }
        }
        else
        {
            /* Task 0 will receive data from all other IO tasks. */
            size_t rlen;    /* Length of IO buffer on this task. */
            int rregions;   /* Number of regions in buffer for this task. */
            size_t start[fndims], count[fndims];
            size_t loffset;

            /* Get the size of the MPI data type. */
            if ((mpierr = MPI_Type_size(basetype, &dsize)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
            LOG((3, "dsize = %d", dsize));

            /* For each of the other tasks that are using this task
             * for IO. */
            for (int rtask = 0; rtask < ios->num_iotasks; rtask++)
            {
                /* From the remote tasks, we send information about
                 * the data regions. and also the data. */
                if (rtask)
                {
                    /* handshake - tell the sending task I'm ready */
                    if ((mpierr = MPI_Send(&ierr, 1, MPI_INT, rtask, 0, ios->io_comm)))
                        return check_mpi(file, mpierr, __FILE__, __LINE__);

                    /* Get length of iobuffer for each field on this
                     * task (all fields are the same length). */
                    if ((mpierr = MPI_Recv(&rlen, 1, MPI_OFFSET, rtask, rtask, ios->io_comm, &status)))
                        return check_mpi(file, mpierr, __FILE__, __LINE__);
                    LOG((3, "received rlen = %d", rlen));

                    /* Get the number of regions, the start/count
                     * values for all regions, and the data buffer. */
                    if (rlen > 0)
                    {
                        if ((mpierr = MPI_Recv(&rregions, 1, MPI_INT, rtask, rtask + ios->num_iotasks,
                                               ios->io_comm, &status)))
                            return check_mpi(file, mpierr, __FILE__, __LINE__);
                        if ((mpierr = MPI_Recv(tmp_start, rregions * fndims, MPI_OFFSET, rtask, rtask + 2 * ios->num_iotasks,
                                               ios->io_comm, &status)))
                            return check_mpi(file, mpierr, __FILE__, __LINE__);
                        if ((mpierr = MPI_Recv(tmp_count, rregions * fndims, MPI_OFFSET, rtask, rtask + 3 * ios->num_iotasks,
                                               ios->io_comm, &status)))
                            return check_mpi(file, mpierr, __FILE__, __LINE__);
                        if ((mpierr = MPI_Recv(iobuf, nvars * rlen, basetype, rtask, rtask + 4 * ios->num_iotasks, ios->io_comm,
                                               &status)))
                            return check_mpi(file, mpierr, __FILE__, __LINE__);
                        LOG((3, "received data rregions = %d fndims = %d", rregions, fndims));
                    }
                }
                else /* task 0 */
                {
                    rlen = llen;
                    rregions = maxregions;
                }
                LOG((3, "rtask = %d rlen = %d rregions = %d", rtask, rlen, rregions));

                /* If there is data from this task, write it. */
                if (rlen > 0)
                {
                    loffset = 0;
                    for (int regioncnt = 0; regioncnt < rregions; regioncnt++)
                    {
                        LOG((3, "writing data for region with regioncnt = %d", regioncnt));

                        /* Get the start/count arrays for this region. */
                        for (int i = 0; i < fndims; i++)
                        {
                            start[i] = tmp_start[i + regioncnt * fndims];
                            count[i] = tmp_count[i + regioncnt * fndims];
                            LOG((3, "start[%d] = %d count[%d] = %d", i, start[i], i, count[i]));
                        }

                        /* Process each variable in the buffer. */
                        for (int nv = 0; nv < nvars; nv++)
                        {
                            LOG((3, "writing buffer var %d", nv));

                            /* Get a pointer to the correct part of the buffer. */
                            bufptr = (void *)((char *)iobuf + tsize * (nv * rlen + loffset));

                            /* If this var has an unlimited dim, set
                             * the start on that dim to the frame
                             * value for this variable. */
                            if (vdesc->record >= 0)
                            {
                                if (fndims > 1 && iodesc_ndims < fndims && count[1] > 0)
                                {
                                    count[0] = 1;
                                    start[0] = frame[nv];
                                }
                                else if (fndims == iodesc_ndims)
                                {
                                    start[0] += vdesc->record;
                                }
                            }

                            /* Call the netCDF functions to write the data. */
                            ierr = nc_put_vara(file->fh, vid[nv], start, count, bufptr);

                            if (ierr)
                            {
                                for (int i = 0; i < fndims; i++)
                                    fprintf(stderr, "vid %d dim %d start %ld count %ld \n", vid[nv], i,
                                            start[i], count[i]);
                                return check_netcdf(file, ierr, __FILE__, __LINE__);
                            }
                        } /* next var */

                        /* Calculate the total size. */
                        size_t tsize = 1;
                        for (int i = 0; i < fndims; i++)
                            tsize *= count[i];

                        /* Keep track of where we are in the buffer. */
                        loffset += tsize;

                        LOG((3, " at bottom of loop regioncnt = %d tsize = %d loffset = %d", regioncnt,
                             tsize, loffset));
                    } /* next regioncnt */
                } /* endif (rlen > 0) */
            } /* next rtask */
        }
    }

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:write_darray_multi_nc_serial");
#endif

    return PIO_NOERR;
}

/**
 * Read an array of data from a file to the (parallel) IO library.
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param iodesc a pointer to the defined iodescriptor for the buffer
 * @param vid the variable id to be read
 * @param iobuf the buffer to be read into from this mpi task. May be
 * null. for example we have 8 ionodes and a distributed array with
 * global size 4, then at least 4 nodes will have a null iobuf. In
 * practice the box rearranger trys to have at least blocksize bytes
 * on each io task and so if the total number of bytes to write is
 * less than blocksize*numiotasks then some iotasks will have a NULL
 * iobuf.
 * @return 0 on success, error code otherwise.
 * @ingroup PIO_read_darray
 */
int pio_read_darray_nc(file_desc_t *file, io_desc_t *iodesc, int vid, void *iobuf)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;     /* Information about the variable. */
    int ndims;             /* Number of dims in decomposition. */
    int fndims;            /* Number of dims for this var in file. */
    int mpierr;            /* Return code from MPI functions. */
    int ierr;              /* Return code from netCDF functions. */

    /* Check inputs. */
    pioassert(file && file->iosystem && iodesc && vid <= PIO_MAX_VARS, "invalid input",
              __FILE__, __LINE__);

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:read_darray_nc");
#endif

    /* Get the IO system info. */
    ios = file->iosystem;

    /* Get the variable info. */
    vdesc = file->varlist + vid;

    /* Get the number of dimensions in the decomposition. */
    ndims = iodesc->ndims;

    /* Get the number of dims for this var in the file. */
    if ((ierr = PIOc_inq_varndims(file->pio_ncid, vid, &fndims)))
        return pio_err(ios, file, ierr, __FILE__, __LINE__);

    /* Is this a non-record var? */
    if (fndims == ndims)
        vdesc->record = -1;

    /* IO procs will actially read the data. */
    if (ios->ioproc)
    {
        io_region *region;
        size_t start[fndims];
        size_t count[fndims];
        size_t tmp_bufsize = 1;
        void *bufptr;
        int tsize;
        int rrlen = 0;
        PIO_Offset *startlist[iodesc->maxregions];
        PIO_Offset *countlist[iodesc->maxregions];

        /* buffer is incremented by byte and loffset is in terms of
           the iodessc->basetype so we need to multiply by the size of
           the basetype. */
        region = iodesc->firstregion;

        /* Get the size of the MPI type. */
        if ((mpierr = MPI_Type_size(iodesc->basetype, &tsize)))
            return check_mpi(file, mpierr, __FILE__, __LINE__);

        /* ??? */
        if (fndims > ndims)
        {
            ndims++;
            if (vdesc->record < 0)
                vdesc->record = 0;
        }

        /* For each regions, read the data. */
        for (int regioncnt = 0; regioncnt < iodesc->maxregions; regioncnt++)
        {
            tmp_bufsize = 1;
            if (region == NULL || iodesc->llen == 0)
            {
                /* No data for this region. */
                for (int i = 0; i < fndims; i++)
                {
                    start[i] = 0;
                    count[i] = 0;
                }
                bufptr = NULL;
            }
            else
            {
                /* Get a pointer where we should put the data we read. */
                if (regioncnt == 0 || region == NULL)
                    bufptr = iobuf;
                else
                    bufptr=(void *)((char *)iobuf + tsize * region->loffset);

                LOG((2, "%d %d %d", iodesc->llen - region->loffset, iodesc->llen, region->loffset));

                /* Get the start/count arrays. */
                if (vdesc->record >= 0 && fndims > 1)
                {
                    /* This is a record var. The unlimited dimension
                     * (0) is handled specially. */
                    start[0] = vdesc->record;
                    for (int i = 1; i < ndims; i++)
                    {
                        start[i] = region->start[i-1];
                        count[i] = region->count[i-1];
                    }

                    /* Read one record. */
                    if (count[1] > 0)
                        count[0] = 1;
                }
                else
                {
                    /* Non-time dependent array */
                    for (int i = 0; i < ndims; i++)
                    {
                        start[i] = region->start[i];
                        count[i] = region->count[i];
                    }
                }
            }

            /* Do the read. */
            switch (file->iotype)
            {
#ifdef _NETCDF4
            case PIO_IOTYPE_NETCDF4P:
                ierr = nc_get_vara(file->fh, vid, start, count, bufptr);
                break;
#endif
#ifdef _PNETCDF
            case PIO_IOTYPE_PNETCDF:
            {
                tmp_bufsize = 1;
                for (int j = 0; j < fndims; j++)
                    tmp_bufsize *= count[j];

                if (tmp_bufsize > 0)
                {
                    startlist[rrlen] = bget(fndims * sizeof(PIO_Offset));
                    countlist[rrlen] = bget(fndims * sizeof(PIO_Offset));

                    for (int j = 0; j < fndims; j++)
                    {
                        startlist[rrlen][j] = start[j];
                        countlist[rrlen][j] = count[j];
                    }
                    rrlen++;
                }

                /* Is this is the last region to process? */
                if (regioncnt == iodesc->maxregions - 1)
                {
                    /* Read a list of subarrays. */
                    ierr = ncmpi_get_varn_all(file->fh, vid, rrlen, startlist,
                                              countlist, iobuf, iodesc->llen, iodesc->basetype);

                    /* Release the start and count arrays. */
                    for (int i = 0; i < rrlen; i++)
                    {
                        brel(startlist[i]);
                        brel(countlist[i]);
                    }
                }
            }
            break;
#endif
            default:
                return pio_err(ios, file, PIO_EBADIOTYPE, __FILE__, __LINE__);
            }

            /* Check return code. */
            if (ierr)
                return check_netcdf(file, ierr, __FILE__,__LINE__);

            /* Move to next region. */
            if (region)
                region = region->next;
        } /* next regioncnt */
    }

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:read_darray_nc");
#endif

    return PIO_NOERR;
}

/**
 * Read an array of data from a file to the (serial) IO library. This
 * function is only used with netCDF classic and netCDF-4 serial
 * iotypes.
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param iodesc a pointer to the defined iodescriptor for the buffer
 * @param vid the variable id to be read.
 * @param iobuf the buffer to be written from this mpi task. May be
 * null. for example we have 8 ionodes and a distributed array with
 * global size 4, then at least 4 nodes will have a null iobuf. In
 * practice the box rearranger trys to have at least blocksize bytes
 * on each io task and so if the total number of bytes to write is
 * less than blocksize * numiotasks then some iotasks will have a NULL
 * iobuf.
 * @returns 0 for success, error code otherwise.
 * @ingroup PIO_read_darray
 */
int pio_read_darray_nc_serial(file_desc_t *file, io_desc_t *iodesc, int vid,
                              void *iobuf)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;     /* Information about the variable. */
    int ndims;             /* Number of dims in decomposition. */
    int fndims;            /* Number of dims for this var in file. */
    MPI_Status status;
    int mpierr;  /* Return code from MPI functions. */
    int ierr;

    /* Check inputs. */
    pioassert(file && file->iosystem && iodesc && vid >= 0 && vid <= PIO_MAX_VARS,
              "invalid input", __FILE__, __LINE__);

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:read_darray_nc_serial");
#endif
    ios = file->iosystem;

    /* Get var info for this var. */
    vdesc = file->varlist + vid;

    /* Get the number of dims in our decomposition. */
    ndims = iodesc->ndims;

    /* Get number of dims for this var. */
    if ((ierr = PIOc_inq_varndims(file->pio_ncid, vid, &fndims)))
        return pio_err(ios, file, ierr, __FILE__, __LINE__);

    /* Is this a non-record var? */
    if (fndims == ndims)
        vdesc->record = -1;

    if (ios->ioproc)
    {
        io_region *region;
        size_t start[fndims];
        size_t count[fndims];
        size_t tmp_start[fndims * iodesc->maxregions];
        size_t tmp_count[fndims * iodesc->maxregions];
        size_t tmp_bufsize;
        void *bufptr;
        int tsize;

        /* buffer is incremented by byte and loffset is in terms of
           the iodessc->basetype so we need to multiply by the size of
           the basetype. */
        region = iodesc->firstregion;

        /* Get the size of the MPI type. */
        if ((mpierr = MPI_Type_size(iodesc->basetype, &tsize)))
            return check_mpi(file, mpierr, __FILE__, __LINE__);

        if (fndims > ndims)
        {
            if (vdesc->record < 0)
                vdesc->record = 0;
        }

        /* Put together start/count arrays for all regions. */
        for (int regioncnt = 0; regioncnt < iodesc->maxregions; regioncnt++)
        {
            if (!region || iodesc->llen == 0)
            {
                /* Nothing to write for this region. */
                for (int i = 0; i < fndims; i++)
                {
                    tmp_start[i + regioncnt * fndims] = 0;
                    tmp_count[i + regioncnt * fndims] = 0;
                }
                bufptr = NULL;
            }
            else
            {
                if (vdesc->record >= 0 && fndims > 1)
                {
                    /* This is a record var. Find start for record dims. */
                    tmp_start[regioncnt * fndims] = vdesc->record;

                    /* Find start/count for all non-record dims. */
                    for (int i = 1; i < fndims; i++)
                    {
                        tmp_start[i + regioncnt * fndims] = region->start[i - 1];
                        tmp_count[i + regioncnt * fndims] = region->count[i - 1];
                    }

                    /* Set count for record dimension. */
                    if (tmp_count[1 + regioncnt * fndims] > 0)
                        tmp_count[regioncnt * fndims] = 1;
                }
                else
                {
                    /* Non-time dependent array */
                    for (int i = 0; i < fndims; i++)
                    {
                        tmp_start[i + regioncnt * fndims] = region->start[i];
                        tmp_count[i + regioncnt * fndims] = region->count[i];
                    }
                }
            }

#if PIO_ENABLE_LOGGING
            /* Log arrays for debug purposes. */
            LOG((3, "region = %d", region));
            for (int i = 0; i < fndims; i++)
                LOG((3, "tmp_start[%d] = %d tmp_count[%d] = %d", i + regioncnt * fndims, tmp_start[i + regioncnt * fndims],
                     i + regioncnt * fndims, tmp_count[i + regioncnt * fndims]));
#endif /* PIO_ENABLE_LOGGING */

            /* Move to next region. */
            if (region)
                region = region->next;
        } /* next regioncnt */

        /* IO tasks other than 0 send their starts/counts and data to
         * IO task 0. */
        if (ios->io_rank > 0)
        {
            if ((mpierr = MPI_Send(&iodesc->llen, 1, MPI_OFFSET, 0, ios->io_rank, ios->io_comm)))
                return check_mpi(file, mpierr, __FILE__, __LINE__);
            LOG((3, "sent iodesc->llen = %d", iodesc->llen));

            if (iodesc->llen > 0)
            {
                if ((mpierr = MPI_Send(&(iodesc->maxregions), 1, MPI_INT, 0,
                                       ios->num_iotasks + ios->io_rank, ios->io_comm)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(tmp_count, iodesc->maxregions * fndims, MPI_OFFSET, 0,
                                       2 * ios->num_iotasks + ios->io_rank, ios->io_comm)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(tmp_start, iodesc->maxregions * fndims, MPI_OFFSET, 0,
                                       3 * ios->num_iotasks + ios->io_rank, ios->io_comm)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                LOG((3, "sent iodesc->maxregions = %d tmp_count and tmp_start arrays", iodesc->maxregions));

                if ((mpierr = MPI_Recv(iobuf, iodesc->llen, iodesc->basetype, 0,
                                       4 * ios->num_iotasks + ios->io_rank, ios->io_comm, &status)))
                    return check_mpi(file, mpierr, __FILE__, __LINE__);
                LOG((3, "received %d elements of data", iodesc->llen));
            }
        }
        else if (ios->io_rank == 0)
        {
            /* This is IO task 0. Get starts/counts and data from
             * other IO tasks. */
            int maxregions = 0;
            size_t loffset, regionsize;
            size_t this_start[fndims * iodesc->maxregions];
            size_t this_count[fndims * iodesc->maxregions];

            for (int rtask = 1; rtask <= ios->num_iotasks; rtask++)
            {
                if (rtask < ios->num_iotasks)
                {
                    if ((mpierr = MPI_Recv(&tmp_bufsize, 1, MPI_OFFSET, rtask, rtask, ios->io_comm, &status)))
                        return check_mpi(file, mpierr, __FILE__, __LINE__);
                    LOG((3, "received tmp_bufsize = %d", tmp_bufsize));

                    if (tmp_bufsize > 0)
                    {
                        if ((mpierr = MPI_Recv(&maxregions, 1, MPI_INT, rtask, ios->num_iotasks + rtask,
                                               ios->io_comm, &status)))
                            return check_mpi(file, mpierr, __FILE__, __LINE__);
                        if ((mpierr = MPI_Recv(this_count, maxregions * fndims, MPI_OFFSET, rtask,
                                               2 * ios->num_iotasks + rtask, ios->io_comm, &status)))
                            return check_mpi(file, mpierr, __FILE__, __LINE__);
                        if ((mpierr = MPI_Recv(this_start, maxregions * fndims, MPI_OFFSET, rtask,
                                               3 * ios->num_iotasks + rtask, ios->io_comm, &status)))
                            return check_mpi(file, mpierr, __FILE__, __LINE__);
                        LOG((3, "received maxregions = %d this_count, this_start arrays ", maxregions));
                    }
                }
                else
                {
                    maxregions = iodesc->maxregions;
                    tmp_bufsize = iodesc->llen;
                }
                LOG((3, "maxregions = %d tmp_bufsize = %d", maxregions, tmp_bufsize));

                /* Now get each region of data. */
                loffset = 0;
                for (int regioncnt = 0; regioncnt < maxregions; regioncnt++)
                {
                    /* Get pointer where data should go. */
                    bufptr = (void *)((char *)iobuf + tsize * loffset);
                    regionsize = 1;

                    /* ??? */
                    if (rtask < ios->num_iotasks)
                    {
                        for (int m = 0; m < fndims; m++)
                        {
                            start[m] = this_start[m + regioncnt * fndims];
                            count[m] = this_count[m + regioncnt * fndims];
                            regionsize *= count[m];
                        }
                    }
                    else
                    {
                        for (int m = 0; m < fndims; m++)
                        {
                            start[m] = tmp_start[m + regioncnt * fndims];
                            count[m] = tmp_count[m + regioncnt * fndims];
                            regionsize *= count[m];
                        }
                    }
                    loffset += regionsize;

                    /* Read the data. */
                    ierr = nc_get_vara(file->fh, vid, start, count, bufptr);

                    /* Check error code of netCDF call. */
                    if (ierr)
                        return check_netcdf(file, ierr, __FILE__, __LINE__);
                }

                /* The decomposition may not use all of the active io
                 * tasks. rtask here is the io task rank and
                 * ios->num_iotasks is the number of iotasks actually
                 * used in this decomposition. */
                if (rtask < ios->num_iotasks && tmp_bufsize > 0)
                    if ((mpierr = MPI_Send(iobuf, tmp_bufsize, iodesc->basetype, rtask,
                                           4 * ios->num_iotasks + rtask, ios->io_comm)))
                        return check_mpi(file, mpierr, __FILE__, __LINE__);
            }
        }
    }

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:read_darray_nc_serial");
#endif

    return PIO_NOERR;
}

/**
 * Flush the output buffer. This is only relevant for files opened
 * with pnetcdf.
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param force true to force the flushing of the buffer
 * @param addsize additional size to add to buffer (in bytes)
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int flush_output_buffer(file_desc_t *file, bool force, PIO_Offset addsize)
{
    int mpierr;  /* Return code from MPI functions. */
    int ierr = PIO_NOERR;

#ifdef _PNETCDF
    var_desc_t *vdesc;
    PIO_Offset usage = 0;

    /* Check inputs. */
    pioassert(file, "invalid input", __FILE__, __LINE__);

    /* Find out the buffer usage. */
    if ((ierr = ncmpi_inq_buffer_usage(file->fh, &usage)))
	return pio_err(NULL, file, PIO_EBADID, __FILE__, __LINE__);

    /* If we are not forcing a flush, spread the usage to all IO
     * tasks. */
    if (!force && file->iosystem->io_comm != MPI_COMM_NULL)
    {
        usage += addsize;
        if ((mpierr = MPI_Allreduce(MPI_IN_PLACE, &usage, 1,  MPI_OFFSET,  MPI_MAX,
                                    file->iosystem->io_comm)))
            return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* Keep track of the maximum usage. */
    if (usage > maxusage)
        maxusage = usage;

    /* If the user forces it, or the buffer has exceeded the size
     * limit, then flush to disk. */
    if (force || usage >= pio_buffer_size_limit)
    {
        int rcnt;
        int  maxreq;
        int reqcnt;
        maxreq = 0;
        reqcnt = 0;
        rcnt = 0;
        for (int i = 0; i < PIO_MAX_VARS; i++)
        {
            vdesc = file->varlist + i;
            reqcnt += vdesc->nreqs;
            if (vdesc->nreqs > 0)
                maxreq = i;
        }
        int request[reqcnt];
        int status[reqcnt];

        for (int i = 0; i <= maxreq; i++)
        {
            vdesc = file->varlist + i;
#ifdef MPIO_ONESIDED
            /*onesided optimization requires that all of the requests in a wait_all call represent
              a contiguous block of data in the file */
            if (rcnt > 0 && (prev_record != vdesc->record || vdesc->nreqs==0))
            {
                ierr = ncmpi_wait_all(file->fh, rcnt, request, status);
                rcnt = 0;
            }
            prev_record = vdesc->record;
#endif
            for (reqcnt = 0; reqcnt < vdesc->nreqs; reqcnt++)
                request[rcnt++] = max(vdesc->request[reqcnt], NC_REQ_NULL);

            if (vdesc->request != NULL)
                free(vdesc->request);
            vdesc->request = NULL;
            vdesc->nreqs = 0;

#ifdef FLUSH_EVERY_VAR
            ierr = ncmpi_wait_all(file->fh, rcnt, request, status);
            rcnt = 0;
#endif
        }

        if (rcnt > 0)
            ierr = ncmpi_wait_all(file->fh, rcnt, request, status);

        /* Release resources. */
        for (int i = 0; i < PIO_MAX_VARS; i++)
        {
            vdesc = file->varlist + i;
            if (vdesc->iobuf)
            {
                brel(vdesc->iobuf);
                vdesc->iobuf = NULL;
            }
            if (vdesc->fillbuf)
            {
                brel(vdesc->fillbuf);
                vdesc->fillbuf = NULL;
            }
        }
    }

#endif /* _PNETCDF */
    return ierr;
}

/**
 * Print out info about the buffer for debug purposes. This should
 * only be called when logging is enabled.
 *
 * @param ios pointer to the IO system structure
 * @param collective true if collective report is desired
 * @ingroup PIO_write_darray
 */
void cn_buffer_report(iosystem_desc_t *ios, bool collective)
{
    int mpierr;  /* Return code from MPI functions. */

    LOG((2, "cn_buffer_report ios->iossysid = %d collective = %d CN_bpool = %d",
         ios->iosysid, collective, CN_bpool));
    if (CN_bpool)
    {
        long bget_stats[5];
        long bget_mins[5];
        long bget_maxs[5];

        bstats(bget_stats, bget_stats+1,bget_stats+2,bget_stats+3,bget_stats+4);
        if (collective)
        {
            LOG((3, "cn_buffer_report calling MPI_Reduce ios->comp_comm = %d", ios->comp_comm));
            if ((mpierr = MPI_Reduce(bget_stats, bget_maxs, 5, MPI_LONG, MPI_MAX, 0, ios->comp_comm)))
                check_mpi(NULL, mpierr, __FILE__, __LINE__);
            LOG((3, "cn_buffer_report calling MPI_Reduce"));
            if ((mpierr = MPI_Reduce(bget_stats, bget_mins, 5, MPI_LONG, MPI_MIN, 0, ios->comp_comm)))
                check_mpi(NULL, mpierr, __FILE__, __LINE__);
            if (ios->compmaster == MPI_ROOT)
            {
                LOG((1, "Currently allocated buffer space %ld %ld", bget_mins[0], bget_maxs[0]));
                LOG((1, "Currently available buffer space %ld %ld", bget_mins[1], bget_maxs[1]));
                LOG((1, "Current largest free block %ld %ld", bget_mins[2], bget_maxs[2]));
                LOG((1, "Number of successful bget calls %ld %ld", bget_mins[3], bget_maxs[3]));
                LOG((1, "Number of successful brel calls  %ld %ld", bget_mins[4], bget_maxs[4]));
            }
        }
        else
        {
            LOG((1, "Currently allocated buffer space %ld", bget_stats[0]));
            LOG((1, "Currently available buffer space %ld", bget_stats[1]));
            LOG((1, "Current largest free block %ld", bget_stats[2]));
            LOG((1, "Number of successful bget calls %ld", bget_stats[3]));
            LOG((1, "Number of successful brel calls  %ld", bget_stats[4]));
        }
    }
}

/**
 * Free the buffer pool. If malloc is used (that is, PIO_USE_MALLOC is
 * non zero), this function does nothing.
 *
 * @param ios pointer to the IO system structure.
 * @ingroup PIO_write_darray
 */
void free_cn_buffer_pool(iosystem_desc_t *ios)
{
#if !PIO_USE_MALLOC
    LOG((2, "free_cn_buffer_pool CN_bpool = %d", CN_bpool));
    if (CN_bpool)
    {
        cn_buffer_report(ios, false);
        bpoolrelease(CN_bpool);
        LOG((2, "free_cn_buffer_pool done!"));
        free(CN_bpool);
        CN_bpool = NULL;
    }
#endif /* !PIO_USE_MALLOC */
}

/**
 * Flush the buffer.
 *
 * @param ncid identifies the netCDF file
 * @param wmb May be NULL, in which case function returns.
 * @param flushtodisk
 * @returns 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int flush_buffer(int ncid, wmulti_buffer *wmb, bool flushtodisk)
{
    file_desc_t *file;
    int ret;

    /* Check input. */
    pioassert(wmb, "invalid input", __FILE__, __LINE__);

    /* Get the file info (to get error handler). */
    if ((ret = pio_get_file(ncid, &file)))
        return pio_err(NULL, NULL, ret, __FILE__, __LINE__);

    LOG((1, "flush_buffer ncid = %d flushtodisk = %d", ncid, flushtodisk));

    /* If there are any variables in this buffer... */
    if (wmb->validvars > 0)
    {
        /* Write any data in the buffer. */
        ret = PIOc_write_darray_multi(ncid, wmb->vid,  wmb->ioid, wmb->validvars,
                                      wmb->arraylen, wmb->data, wmb->frame,
                                      wmb->fillvalue, flushtodisk);

        wmb->validvars = 0;

        /* Release the list of variable IDs. */
        brel(wmb->vid);
        wmb->vid = NULL;

        /* Release the data memory. */
        brel(wmb->data);
        wmb->data = NULL;

        /* If there is a fill value, release it. */
        if (wmb->fillvalue)
            brel(wmb->fillvalue);
        wmb->fillvalue = NULL;

        /* Release the record number. */
        if (wmb->frame)
            brel(wmb->frame);
        wmb->frame = NULL;

        if (ret)
            return pio_err(NULL, file, ret, __FILE__, __LINE__);
    }

    return PIO_NOERR;
}

/**
 * Compute the maximum aggregate number of bytes.
 *
 * @param ios the IO system structure
 * @param iodesc a pointer to the defined iodescriptor for the
 * buffer. If NULL, function returns immediately.
 * @returns 0 for success, error code otherwise.
 */
int compute_maxaggregate_bytes(iosystem_desc_t *ios, io_desc_t *iodesc)
{
    int maxbytesoniotask = INT_MAX;
    int maxbytesoncomputetask = INT_MAX;
    int maxbytes;
    int mpierr;  /* Return code from MPI functions. */

    /* Check inputs. */
    pioassert(iodesc, "invalid input", __FILE__, __LINE__);

    LOG((2, "compute_maxaggregate_bytes iodesc->maxiobuflen = %d iodesc->ndof = %d",
         iodesc->maxiobuflen, iodesc->ndof));

    if (ios->ioproc && iodesc->maxiobuflen > 0)
        maxbytesoniotask = pio_buffer_size_limit / iodesc->maxiobuflen;

    if (ios->comp_rank >= 0 && iodesc->ndof > 0)
        maxbytesoncomputetask = pio_cnbuffer_limit / iodesc->ndof;

    maxbytes = min(maxbytesoniotask, maxbytesoncomputetask);
    LOG((2, "compute_maxaggregate_bytes maxbytesoniotask = %d maxbytesoncomputetask = %d",
         maxbytesoniotask, maxbytesoncomputetask));

    if ((mpierr = MPI_Allreduce(MPI_IN_PLACE, &maxbytes, 1, MPI_INT, MPI_MIN,
                                ios->union_comm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    iodesc->maxbytes = maxbytes;

    return PIO_NOERR;
}
