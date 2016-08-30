/** @file
 *
 * This file contains the routines that read and write
 * distributed arrays in PIO.
 *
 * When arrays are distributed, each processor holds some of the
 * array. Only by combining the distributed arrays from all processor
 * can the full array be obtained.
 *
 * @author Jim Edwards, Ed Hartnett
 */

#include <config.h>
#include <pio.h>
#include <pio_internal.h>

/* 10MB default limit. */
PIO_Offset PIO_BUFFER_SIZE_LIMIT = 10485760;

/* Initial size of compute buffer. */
bufsize PIO_CNBUFFER_LIMIT = 33554432;

/* Global buffer pool pointer. */
static void *CN_bpool = NULL;

/* Maximum buffer usage. */
static PIO_Offset maxusage = 0;

/** Set the pio buffer size limit. This is the size of the data buffer
 * on the IO nodes.
 *
 * The pio_buffer_size_limit will only apply to files opened after
 * the setting is changed.
 *
 * @param limit the size of the buffer on the IO nodes
 *
 * @return The previous limit setting.
 */
PIO_Offset PIOc_set_buffer_size_limit(const PIO_Offset limit)
{
    PIO_Offset oldsize;
    oldsize = PIO_BUFFER_SIZE_LIMIT;
    if (limit > 0)
	PIO_BUFFER_SIZE_LIMIT = limit;
    return oldsize;
}

/** Initialize the compute buffer to size PIO_CNBUFFER_LIMIT.
 *
 * This routine initializes the compute buffer pool if the bget memory
 * management is used. If malloc is used (that is, PIO_USE_MALLOC is
 * non zero), this function does nothing.
 *
 * @param ios the iosystem descriptor which will use the new buffer
 */
void compute_buffer_init(iosystem_desc_t ios)
{
#if !PIO_USE_MALLOC

    if (!CN_bpool)
    {
	if (!(CN_bpool = malloc(PIO_CNBUFFER_LIMIT)))
	{
	    char errmsg[180];
	    sprintf(errmsg,"Unable to allocate a buffer pool of size %d on task %d:"
		    " try reducing PIO_CNBUFFER_LIMIT\n", PIO_CNBUFFER_LIMIT, ios.comp_rank);
	    piodie(errmsg, __FILE__, __LINE__);
	}
	
	bpool(CN_bpool, PIO_CNBUFFER_LIMIT);
	if (!CN_bpool)
	{
	    char errmsg[180];
	    sprintf(errmsg,"Unable to allocate a buffer pool of size %d on task %d:"
		    " try reducing PIO_CNBUFFER_LIMIT\n", PIO_CNBUFFER_LIMIT, ios.comp_rank);
	    piodie(errmsg, __FILE__, __LINE__);
	}
	
	bectl(NULL, malloc, free, PIO_CNBUFFER_LIMIT);
    }
#endif
}

/** Write a single distributed field to output. This routine is only
 * used if aggregation is off.
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param iodesc a pointer to the defined iodescriptor for the buffer
 * @param vid the variable id to be written
 * @param IOBUF the buffer to be written from this mpi task
 * @param fillvalue the optional fillvalue to be used for missing
 * data in this buffer
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int pio_write_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid,
			void *IOBUF, void *fillvalue)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;
    int ndims;             /* Number of dimensions according to iodesc. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int i;                 /* Loop counter. */
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int dsize;             /* Size of the type. */
    MPI_Status status;     /* Status from MPI_Recv calls. */
    PIO_Offset usage;      /* Size of current buffer. */
    int fndims;            /* Number of dims for variable according to netCDF. */
    PIO_Offset tdsize = 0; /* Total size. */

    LOG((1, "pio_write_array_nc vid = %d", vid));

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:write_darray_nc");
#endif

    /* Get the IO system info. */
    if (!(ios = file->iosystem))
	return PIO_EBADID;

    /* Get pointer to variable information. */
    if (!(vdesc = file->varlist + vid))
	return PIO_EBADID;

    ndims = iodesc->ndims;

    /* Get the number of dims for this var from netcdf. */
    ierr = PIOc_inq_varndims(file->fh, vid, &fndims);

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = 0;

	    if (ios->compmaster)
		mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

	    if (!mpierr)
		mpierr = MPI_Bcast(&file->fh, 1, MPI_INT, ios->compmaster, ios->intercomm);
	}
    }

    /* If this is an IO task, write the data. */
    if (ios->ioproc)
    {
	io_region *region;
	int regioncnt;
	int rrcnt;
	void *bufptr;
	void *tmp_buf = NULL;
	int tsize;            /* Type size. */
	size_t start[fndims]; /* Local start array for this task. */
	size_t count[fndims]; /* Local count array for this task. */
	int buflen;
	int j;                /* Loop counter. */

	PIO_Offset *startlist[iodesc->maxregions];
	PIO_Offset *countlist[iodesc->maxregions];

	/* Get the type size (again?) */
	MPI_Type_size(iodesc->basetype, &tsize);

	region = iodesc->firstregion;

	/* If this is a var with an unlimited dimension, and the
	 * iodesc ndims doesn't contain it, then add it to ndims. */
	if (vdesc->record >= 0 && ndims < fndims)
	    ndims++;

#ifdef _PNETCDF
	/* Make sure we have room in the buffer. */
	if (file->iotype == PIO_IOTYPE_PNETCDF)
	    flush_output_buffer(file, false, tsize * (iodesc->maxiobuflen));
#endif

	rrcnt = 0;
	/* For each region, figure out start/count arrays. */
	for (regioncnt = 0; regioncnt < iodesc->maxregions; regioncnt++)
	{
	    /* Init arrays to zeros. */
	    for (i = 0; i < ndims; i++)
	    {
		start[i] = 0;
		count[i] = 0;
	    }
	    
	    if (region)
	    {
		bufptr = (void *)((char *)IOBUF + tsize * region->loffset);
		if (vdesc->record >= 0)
		{
		    /* This is a record based multidimensional array. */

		    /* This does not look correct, but will work if
		     * unlimited dim is dim 0. */
		    start[0] = vdesc->record;

		    /* Set the local start and count arrays. */
		    for (i = 1; i < ndims; i++)
		    {
			start[i] = region->start[i - 1];
			count[i] = region->count[i - 1];
		    }

		    /* If there is data to be written, write one timestep. */
		    if (count[1] > 0)
			count[0] = 1;
		}
		else
		{
		    /* Array without unlimited dimension. */
		    for (i = 0; i < ndims; i++)
		    {
			start[i] = region->start[i];
			count[i] = region->count[i];
		    }
		}
	    }

	    switch(file->iotype)
	    {
#ifdef _NETCDF
#ifdef _NETCDF4
	    case PIO_IOTYPE_NETCDF4P:

		/* Use collective writes with this variable. */
		ierr = nc_var_par_access(file->fh, vid, NC_COLLECTIVE);

		/* Write the data. */
		if (iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8)
		    ierr = nc_put_vara_double(file->fh, vid, (size_t *)start, (size_t *)count,
					      (const double *)bufptr);
		else if (iodesc->basetype == MPI_INTEGER)
		    ierr = nc_put_vara_int(file->fh, vid, (size_t *)start, (size_t *)count,
					   (const int *)bufptr);
		else if (iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4)
		    ierr = nc_put_vara_float(file->fh, vid, (size_t *)start, (size_t *)count,
					     (const float *)bufptr);
		else
		    fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",
			    (int)iodesc->basetype);
		break;
	    case PIO_IOTYPE_NETCDF4C:
#endif /* _NETCDF4 */
	    case PIO_IOTYPE_NETCDF:
	    {
		/* Find the type size (again?) */
		mpierr = MPI_Type_size(iodesc->basetype, &dsize);
		
		size_t tstart[ndims], tcount[ndims];

		/* The IO master task does all the data writes, but
		 * sends the data to the other IO tasks (why?). */
		if (ios->io_rank == 0)
		{
		    for (i = 0; i < iodesc->num_aiotasks; i++)
		    {
			if (i == 0)
			{
			    buflen = 1;
			    for (j = 0; j < ndims; j++)
			    {
				tstart[j] =  start[j];
				tcount[j] =  count[j];
				buflen *= tcount[j];
				tmp_buf = bufptr;
			    }
			}
			else
			{
			    /* Handshake - tell the sending task I'm ready. */
			    mpierr = MPI_Send(&ierr, 1, MPI_INT, i, 0, ios->io_comm);  
			    mpierr = MPI_Recv(&buflen, 1, MPI_INT, i, 1, ios->io_comm, &status);
			    if (buflen > 0)
			    {
				mpierr = MPI_Recv(tstart, ndims, MPI_OFFSET, i, ios->num_iotasks+i,
						  ios->io_comm, &status);
				mpierr = MPI_Recv(tcount, ndims, MPI_OFFSET, i, 2 * ios->num_iotasks + i,
						  ios->io_comm, &status);
				tmp_buf = malloc(buflen * dsize);
				mpierr = MPI_Recv(tmp_buf, buflen, iodesc->basetype, i, i, ios->io_comm, &status);
			    }
			}

			if (buflen > 0)
			{
			    /* Write the data. */
			    if (iodesc->basetype == MPI_INTEGER)
				ierr = nc_put_vara_int(file->fh, vid, tstart, tcount, (const int *)tmp_buf);
			    else if (iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8)
				ierr = nc_put_vara_double(file->fh, vid, tstart, tcount, (const double *)tmp_buf);
			    else if (iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4)
				ierr = nc_put_vara_float(file->fh, vid, tstart, tcount, (const float *)tmp_buf);
			    else
				fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",
					(int)iodesc->basetype);

			    /* Was there an error from netCDF? */
			    if (ierr == PIO_EEDGE)
				for (i = 0; i < ndims; i++)
				    fprintf(stderr,"dim %d start %ld count %ld\n", i, tstart[i], tcount[i]);

			    /* Free the temporary buffer, if we don't need it any more. */
			    if (tmp_buf != bufptr)
				free(tmp_buf);
			}
		    }
		}
		else if (ios->io_rank < iodesc->num_aiotasks)
		{
		    buflen = 1;
		    for (i = 0; i < ndims; i++)
		    {
			tstart[i] = (size_t) start[i];
			tcount[i] = (size_t) count[i];
			buflen *= tcount[i];
			//               printf("%s %d %d %d %d\n",__FILE__,__LINE__,i,tstart[i],tcount[i]);
		    }
		    /*	     printf("%s %d %d %d %d %d %d %d %d %d\n",__FILE__,__LINE__,ios->io_rank,tstart[0],
			     tstart[1],tcount[0],tcount[1],buflen,ndims,fndims);*/
		    mpierr = MPI_Recv(&ierr, 1, MPI_INT, 0, 0, ios->io_comm, &status);  // task0 is ready to recieve
		    mpierr = MPI_Rsend(&buflen, 1, MPI_INT, 0, 1, ios->io_comm);
		    if (buflen > 0)
		    {
			mpierr = MPI_Rsend(tstart, ndims, MPI_OFFSET, 0, ios->num_iotasks+ios->io_rank,
					    ios->io_comm);
			mpierr = MPI_Rsend(tcount, ndims, MPI_OFFSET, 0,2*ios->num_iotasks+ios->io_rank,
					    ios->io_comm);
			mpierr = MPI_Rsend(bufptr, buflen, iodesc->basetype, 0, ios->io_rank, ios->io_comm);
		    }
		}
		break;
	    }
	    break;
#endif /* _NETCDF */
#ifdef _PNETCDF
	    case PIO_IOTYPE_PNETCDF:
		for (i = 0, dsize = 1; i < ndims; i++)
		    dsize *= count[i];

		tdsize += dsize;
		//	 if (dsize==1 && ndims==2)
		//	 printf("%s %d %d\n",__FILE__,__LINE__,iodesc->basetype);

		if (dsize > 0)
		{
		    //	   printf("%s %d %d %d\n",__FILE__,__LINE__,ios->io_rank,dsize);
		    startlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
		    countlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
		    for (i = 0; i < fndims; i++)
		    {
			startlist[rrcnt][i] = start[i];
			countlist[rrcnt][i] = count[i];
		    }
		    rrcnt++;
		}
		if (regioncnt == iodesc->maxregions - 1)
		{
		    // printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,ios->io_rank,iodesc->llen, tdsize);
		    //	   ierr = ncmpi_put_varn_all(file->fh, vid, iodesc->maxregions, startlist, countlist,
		    //			     IOBUF, iodesc->llen, iodesc->basetype);
		    int reqn = 0;

		    if (vdesc->nreqs % PIO_REQUEST_ALLOC_CHUNK == 0 )
		    {
			vdesc->request = realloc(vdesc->request,
						 sizeof(int) * (vdesc->nreqs + PIO_REQUEST_ALLOC_CHUNK));

			for (int i = vdesc->nreqs; i < vdesc->nreqs + PIO_REQUEST_ALLOC_CHUNK; i++)
			    vdesc->request[i] = NC_REQ_NULL;
			reqn = vdesc->nreqs;
		    }
		    else
			while(vdesc->request[reqn] != NC_REQ_NULL)
			    reqn++;

		    ierr = ncmpi_bput_varn(file->fh, vid, rrcnt, startlist, countlist,
					   IOBUF, iodesc->llen, iodesc->basetype, vdesc->request+reqn);
		    if (vdesc->request[reqn] == NC_REQ_NULL)
			vdesc->request[reqn] = PIO_REQ_NULL;  //keeps wait calls in sync
		    vdesc->nreqs = reqn;

		    //	   printf("%s %d %X %d\n",__FILE__,__LINE__,IOBUF,request);
		    for (i=0;i<rrcnt;i++)
		    {
			free(startlist[i]);
			free(countlist[i]);
		    }
		}
		break;
#endif /* _PNETCDF */
	    default:
		ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	    }

	    /* Move to the next region. */
	    if (region)
		region = region->next;
	} //    for (regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
    } // if (ios->ioproc)

    /* Check the error code returned by netCDF. */
    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);
    
#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:write_darray_nc");
#endif

    return ierr;
}

/** Write a set of one or more aggregated arrays to output file.
 *
 * This routine is used if aggregation is enabled, data is already on
 * the io-tasks
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param nvars the number of variables to be written with this
 * decomposition
 * @param vid: an array of the variable ids to be written
 * @param iodesc_ndims: the number of dimensions explicitly in the
 * iodesc
 * @param basetype the basic type of the minimal data unit
 * @param gsize array of the global dimensions of the field to
 * be written
 * @param maxregions max number of blocks to be written from
 * this iotask
 * @param firstregion pointer to the first element of a linked
 * list of region descriptions.
 * @param llen length of the iobuffer on this task for a single
 * field
 * @param maxiobuflen maximum llen participating
 * @param num_aiotasks actual number of iotasks participating
 * @param IOBUF the buffer to be written from this mpi task
 * @param frame the frame or record dimension for each of the nvars
 * variables in IOBUF
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int pio_write_darray_multi_nc(file_desc_t *file, const int nvars, const int *vid,
			      const int iodesc_ndims, MPI_Datatype basetype, const PIO_Offset *gsize,
			      const int maxregions, io_region *firstregion, const PIO_Offset llen,
			      const int maxiobuflen, const int num_aiotasks,
			      void *IOBUF, const int *frame)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;
    int ierr;
    int i;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int dsize;
    MPI_Status status;
    PIO_Offset usage;
    int fndims;
    PIO_Offset tdsize;
    int tsize;
    int ncid;
    tdsize=0;
    ierr = PIO_NOERR;

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:write_darray_multi_nc");
#endif

    ios = file->iosystem;
    if (ios == NULL)
    {
	fprintf(stderr,"Failed to find iosystem handle \n");
	return PIO_EBADID;
    }
    vdesc = (file->varlist)+vid[0];
    ncid = file->fh;

    if (vdesc == NULL)
    {
	fprintf(stderr,"Failed to find variable handle %d\n",vid[0]);
	return PIO_EBADID;
    }

    /* If async is in use, send message to IO master task. */
    if (ios->async_interface)
    {
	if (!ios->ioproc)
	{
	    int msg = 0;
	    if (ios->compmaster)
		mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

	    if (!mpierr)
		mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
	}
    }

    ierr = PIOc_inq_varndims(file->fh, vid[0], &fndims);
    MPI_Type_size(basetype, &tsize);

    if (ios->ioproc)
    {
	io_region *region;
	int regioncnt;
	int rrcnt;
	void *bufptr;
	int buflen, j;
	size_t start[fndims];
	size_t count[fndims];
	int ndims = iodesc_ndims;

	PIO_Offset *startlist[maxregions];
	PIO_Offset *countlist[maxregions];

	ncid = file->fh;
	region = firstregion;

	rrcnt = 0;
	for (regioncnt = 0; regioncnt < maxregions; regioncnt++)
	{
	    // printf("%s %d %d %d %d %d %d\n",__FILE__,__LINE__,region->start[0],region->count[0],ndims,fndims,vdesc->record);
	    for (i = 0; i < fndims; i++)
	    {
		start[i] = 0;
		count[i] = 0;
	    }
	    if (region)
	    {
		// this is a record based multidimensional array
		if (vdesc->record >= 0)
		{
		    for (i = fndims - ndims; i < fndims; i++)
		    {
			start[i] = region->start[i-(fndims-ndims)];
			count[i] = region->count[i-(fndims-ndims)];
		    }

		    if (fndims>1 && ndims<fndims && count[1]>0)
		    {
			count[0] = 1;
			start[0] = frame[0];
		    }
		    else if (fndims==ndims)
		    {
			start[0] += vdesc->record;
		    }
		    // Non-time dependent array
		}
		else
		{
		    for (i = 0; i < ndims; i++)
		    {
			start[i] = region->start[i];
			count[i] = region->count[i];
		    }
		}
	    }

	    switch(file->iotype)
	    {
#ifdef _NETCDF4
	    case PIO_IOTYPE_NETCDF4P:
		for (int nv = 0; nv < nvars; nv++)
		{
		    if (vdesc->record >= 0 && ndims < fndims)
		    {
			start[0] = frame[nv];
		    }
		    if (region)
		    {
			bufptr = (void *)((char *) IOBUF + tsize*(nv*llen + region->loffset));
		    }
		    ierr = nc_var_par_access(ncid, vid[nv], NC_COLLECTIVE);

		    if (basetype == MPI_DOUBLE ||basetype == MPI_REAL8)
		    {
			ierr = nc_put_vara_double (ncid, vid[nv],(size_t *) start,(size_t *) count,
						   (const double *)bufptr);
		    }
		    else if (basetype == MPI_INTEGER)
		    {
			ierr = nc_put_vara_int (ncid, vid[nv], (size_t *) start, (size_t *) count,
						(const int *)bufptr);
		    }
		    else if (basetype == MPI_FLOAT || basetype == MPI_REAL4)
		    {
			ierr = nc_put_vara_float (ncid, vid[nv], (size_t *) start, (size_t *) count,
						  (const float *)bufptr);
		    }
		    else
		    {
			fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",
				(int)basetype);
		    }
		}
		break;
#endif
#ifdef _PNETCDF
	    case PIO_IOTYPE_PNETCDF:
		for (i = 0, dsize = 1; i < fndims; i++)
		{
		    dsize *= count[i];
		}
		tdsize += dsize;

		if (dsize>0)
		{
		    //	   printf("%s %d %d %d\n",__FILE__,__LINE__,ios->io_rank,dsize);
		    startlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
		    countlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
		    for (i = 0; i < fndims; i++)
		    {
			startlist[rrcnt][i]=start[i];
			countlist[rrcnt][i]=count[i];
		    }
		    rrcnt++;
		}
		if (regioncnt==maxregions-1)
		{
		    //printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,ios->io_rank,iodesc->llen, tdsize);
		    //	   ierr = ncmpi_put_varn_all(ncid, vid, iodesc->maxregions, startlist, countlist,
		    //			     IOBUF, iodesc->llen, iodesc->basetype);

		    //printf("%s %d %ld \n",__FILE__,__LINE__,IOBUF);
		    for (int nv=0; nv<nvars; nv++)
		    {
			vdesc = (file->varlist)+vid[nv];
			if (vdesc->record >= 0 && ndims<fndims)
			{
			    for (int rc=0;rc<rrcnt;rc++)
			    {
				startlist[rc][0] = frame[nv];
			    }
			}
			bufptr = (void *)((char *) IOBUF + nv*tsize*llen);

			int reqn=0;
			if (vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 )
			{
			    vdesc->request = realloc(vdesc->request,
						     sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));

			    for (int i=vdesc->nreqs;i<vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK;i++)
			    {
				vdesc->request[i]=NC_REQ_NULL;
			    }
			    reqn = vdesc->nreqs;
			}
			else
			{
			    while(vdesc->request[reqn] != NC_REQ_NULL)
			    {
				reqn++;
			    }
			}
			ierr = ncmpi_iput_varn(ncid, vid[nv], rrcnt, startlist, countlist,
					       bufptr, llen, basetype, vdesc->request+reqn);
			/*
			  ierr = ncmpi_bput_varn(ncid, vid[nv], rrcnt, startlist, countlist,
			  bufptr, llen, basetype, &(vdesc->request));
			*/
			if (vdesc->request[reqn] == NC_REQ_NULL)
			{
			    vdesc->request[reqn] = PIO_REQ_NULL;  //keeps wait calls in sync
			}
			vdesc->nreqs += reqn+1;

			//	     printf("%s %d %d %d\n",__FILE__,__LINE__,vdesc->nreqs,vdesc->request[reqn]);
		    }
		    for (i=0;i<rrcnt;i++)
		    {
			if (ierr != PIO_NOERR)
			{
			    for (j=0;j<fndims;j++)
			    {
				printf("pio_darray: %d %d %d %ld %ld \n",__LINE__,i,j,startlist[i][j],countlist[i][j]);
			    }
			}
			free(startlist[i]);
			free(countlist[i]);
		    }
		}
		break;
#endif
	    default:
		ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	    }
	    if (region)
		region = region->next;
	} //    for (regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
    } // if (ios->ioproc)

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:write_darray_multi_nc");
#endif

    return ierr;
}

/** Write a set of one or more aggregated arrays to output file in
 * serial mode.
 *
 * This routine is used if aggregation is enabled, data is already on the
 * io-tasks
 *
 * @param file: a pointer to the open file descriptor for the file
 * that will be written to
 * @param nvars: the number of variables to be written with this
 * decomposition
 * @param vid: an array of the variable ids to be written
 * @param iodesc_ndims: the number of dimensions explicitly in the
 * iodesc
 * @param basetype : the basic type of the minimal data unit
 * @param gsize : array of the global dimensions of the field to be
 * written
 * @param maxregions : max number of blocks to be written from this
 * iotask
 * @param firstregion : pointer to the first element of a linked
 * list of region descriptions.
 * @param llen : length of the iobuffer on this task for a single
 * field
 * @param maxiobuflen : maximum llen participating
 * @param num_aiotasks : actual number of iotasks participating
 * @param IOBUF: the buffer to be written from this mpi task
 * @param frame : the frame or record dimension for each of the
 * nvars variables in IOBUF
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int pio_write_darray_multi_nc_serial(file_desc_t *file, const int nvars, const int *vid,
				     const int iodesc_ndims, MPI_Datatype basetype, const PIO_Offset *gsize,
				     const int maxregions, io_region *firstregion, const PIO_Offset llen,
				     const int maxiobuflen, const int num_aiotasks,
				     void *IOBUF, const int *frame)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;
    int ierr;
    int i;
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
    int dsize;
    MPI_Status status;
    PIO_Offset usage;
    int fndims;
    PIO_Offset tdsize;
    int tsize;
    int ncid;
    tdsize=0;
    ierr = PIO_NOERR;
#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:write_darray_multi_nc_serial");
#endif

    if (!(ios = file->iosystem))
    {
	fprintf(stderr,"Failed to find iosystem handle \n");
	return PIO_EBADID;
    }

    ncid = file->fh;

    if (!(vdesc = (file->varlist) + vid[0]))
    {
	fprintf(stderr,"Failed to find variable handle %d\n",vid[0]);
	return PIO_EBADID;
    }

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
	if (! ios->ioproc)
	{
	    int msg = 0;
	    
	    if (ios->comp_rank==0)
		mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);

	    if (!mpierr)
		mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
	}
    }

    ierr = PIOc_inq_varndims(file->fh, vid[0], &fndims);
    MPI_Type_size(basetype, &tsize);

    if (ios->ioproc)
    {
	io_region *region;
	int regioncnt;
	int rrcnt;
	void *bufptr;
	int buflen, j;
	size_t tmp_start[fndims*maxregions];
	size_t tmp_count[fndims*maxregions];

	int ndims = iodesc_ndims;

	ncid = file->fh;
	region = firstregion;


	rrcnt = 0;
	for (regioncnt = 0; regioncnt < maxregions; regioncnt++)
	{
	    for (i = 0; i < fndims; i++)
	    {
		tmp_start[i + regioncnt * fndims] = 0;
		tmp_count[i + regioncnt * fndims] = 0;
	    }
	    if (region)
	    {
		// this is a record based multidimensional array
		if (vdesc->record >= 0)
		{
		    for (i = fndims - ndims; i < fndims; i++)
		    {
			tmp_start[i + regioncnt * fndims] = region->start[i - (fndims - ndims)];
			tmp_count[i + regioncnt * fndims] = region->count[i - (fndims - ndims)];
		    }
		    // Non-time dependent array
		}
		else
		{
		    for (i = 0; i < ndims; i++)
		    {
			tmp_start[i + regioncnt * fndims] = region->start[i];
			tmp_count[i + regioncnt * fndims] = region->count[i];
		    }
		}
		region = region->next;
	    }
	}
	if (ios->io_rank > 0)
	{
	    mpierr = MPI_Recv(&ierr, 1, MPI_INT, 0, 0, ios->io_comm, &status);  // task0 is ready to recieve
	    MPI_Send(&llen, 1, MPI_OFFSET, 0, ios->io_rank, ios->io_comm);
	    if (llen>0)
	    {
		MPI_Send(&maxregions, 1, MPI_INT, 0, ios->io_rank+ios->num_iotasks, ios->io_comm);
		MPI_Send(tmp_start, maxregions*fndims, MPI_OFFSET, 0, ios->io_rank+2*ios->num_iotasks, ios->io_comm);
		MPI_Send(tmp_count, maxregions*fndims, MPI_OFFSET, 0, ios->io_rank+3*ios->num_iotasks, ios->io_comm);
		//	 printf("%s %d %ld\n",__FILE__,__LINE__,nvars*llen);
		MPI_Send(IOBUF, nvars*llen, basetype, 0, ios->io_rank+4*ios->num_iotasks, ios->io_comm);
	    }
	}
	else
	{
	    size_t rlen;
	    int rregions;
	    size_t start[fndims], count[fndims];
	    size_t loffset;
	    mpierr = MPI_Type_size(basetype, &dsize);

	    for (int rtask=0; rtask<ios->num_iotasks; rtask++)
	    {
		if (rtask>0)
		{
		    mpierr = MPI_Send(&ierr, 1, MPI_INT, rtask, 0, ios->io_comm);  // handshake - tell the sending task I'm ready
		    MPI_Recv(&rlen, 1, MPI_OFFSET, rtask, rtask, ios->io_comm, &status);
		    if (rlen>0){
			MPI_Recv(&rregions, 1, MPI_INT, rtask, rtask+ios->num_iotasks, ios->io_comm, &status);
			MPI_Recv(tmp_start, rregions*fndims, MPI_OFFSET, rtask, rtask+2*ios->num_iotasks, ios->io_comm, &status);
			MPI_Recv(tmp_count, rregions*fndims, MPI_OFFSET, rtask, rtask+3*ios->num_iotasks, ios->io_comm, &status);
			//	     printf("%s %d %d %ld\n",__FILE__,__LINE__,rtask,nvars*rlen);
			MPI_Recv(IOBUF, nvars*rlen, basetype, rtask, rtask+4*ios->num_iotasks, ios->io_comm, &status);
		    }
		}
		else
		{
		    rlen = llen;
		    rregions = maxregions;
		}
		if (rlen>0)
		{
		    loffset = 0;
		    for (regioncnt=0;regioncnt<rregions;regioncnt++)
		    {
			for (int i=0;i<fndims;i++)
			{
			    start[i] = tmp_start[i+regioncnt*fndims];
			    count[i] = tmp_count[i+regioncnt*fndims];
			}

			for (int nv=0; nv<nvars; nv++)
			{
			    bufptr = (void *)((char *) IOBUF+ tsize*(nv*rlen + loffset));

			    if (vdesc->record>=0)
			    {
				if (fndims>1 && ndims<fndims && count[1]>0)
				{
				    count[0] = 1;
				    start[0] = frame[nv];
				}
				else if (fndims==ndims)
				{
				    start[0]+=vdesc->record;
				}
			    }

			    if (basetype == MPI_INTEGER)
			    {
				ierr = nc_put_vara_int (ncid, vid[nv], start, count, (const int *) bufptr);
			    }
			    else if (basetype == MPI_DOUBLE || basetype == MPI_REAL8)
			    {
				ierr = nc_put_vara_double (ncid, vid[nv], start, count, (const double *) bufptr);
			    }
			    else if (basetype == MPI_FLOAT || basetype == MPI_REAL4)
			    {
				ierr = nc_put_vara_float (ncid,vid[nv], start, count, (const float *) bufptr);
			    }
			    else
			    {
				fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) basetype);
			    }

			    if (ierr != PIO_NOERR){
				for (i=0;i<fndims;i++)
				    fprintf(stderr,"vid %d dim %d start %ld count %ld \n",vid[nv],i,start[i],count[i]);
			    }
			}
			size_t tsize;
			tsize = 1;
			for (int i=0;i<fndims;i++)
			{
			    tsize*=count[i];
			}
			loffset += tsize;
		    }//    for (regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
		} // if (rlen>0)
	    } //        for (int rtask=0; rtask<ios->num_iotasks; rtask++){

	}
    } // if (ios->ioproc)

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:write_darray_multi_nc_serial");
#endif

    return ierr;
}

/** Write one or more arrays with the same IO decomposition to the file.
 *
 * @param ncid identifies the netCDF file
 * @param vid: an array of the variable ids to be written
 * @param ioid: the I/O description ID as passed back by
 * PIOc_InitDecomp().
 * @param nvars the number of variables to be written with this
 * decomposition
 * @param arraylen: the length of the array to be written. This
 * is the length of the distrubited array. That is, the length of
 * the portion of the data that is on the processor.
 * @param array: pointer to the data to be written. This is a
 * pointer to the distributed portion of the array that is on this
 * processor.
 * @param frame the frame or record dimension for each of the nvars
 * variables in IOBUF
 * @param fillvalue: pointer to the fill value to be used for
 * missing data.
 * @param flushtodisk
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_write_darray
 */
int PIOc_write_darray_multi(const int ncid, const int *vid, const int ioid,
			    const int nvars, const PIO_Offset arraylen,
			    void *array, const int *frame, void **fillvalue,
			    bool flushtodisk)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;
    io_desc_t *iodesc;

    int vsize, rlen;
    int ierr;
    var_desc_t *vdesc0;

    ierr = PIO_NOERR;

    file = pio_get_file_from_id(ncid);
    if (file == NULL)
    {
	fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
	return PIO_EBADID;
    }
    if (! (file->mode & PIO_WRITE))
    {
	fprintf(stderr,"ERROR:  Attempt to write to read-only file\n");
	return PIO_EPERM;
    }

    iodesc = pio_get_iodesc_from_id(ioid);
    if (iodesc == NULL)
    {
	//     print_trace(NULL);
	//fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
	return PIO_EBADID;
    }

    vdesc0 = file->varlist+vid[0];

    pioassert(nvars>0,"nvars <= 0",__FILE__,__LINE__);

    ios = file->iosystem;
    //   rlen = iodesc->llen*nvars;
    rlen=0;
    if (iodesc->llen>0)
    {
	rlen = iodesc->maxiobuflen*nvars;
    }
    if (vdesc0->iobuf)
    {
	piodie("Attempt to overwrite existing io buffer",__FILE__,__LINE__);
    }
    if (iodesc->rearranger>0)
    {
	if (rlen>0)
	{
	    MPI_Type_size(iodesc->basetype, &vsize);
	    //printf("rlen*vsize = %ld\n",rlen*vsize);

	    vdesc0->iobuf = bget((size_t) vsize* (size_t) rlen);
	    if (vdesc0->iobuf==NULL)
	    {
		printf("%s %d %d %ld\n",__FILE__,__LINE__,nvars,vsize*rlen);
		piomemerror(*ios,(size_t) rlen*(size_t) vsize, __FILE__,__LINE__);
	    }
	    if (iodesc->needsfill && iodesc->rearranger==PIO_REARR_BOX)
	    {
		if (vsize==4)
		{
		    for (int nv=0;nv < nvars; nv++)
		    {
			for (int i=0;i<iodesc->maxiobuflen;i++)
			{
			    ((float *) vdesc0->iobuf)[i+nv*(iodesc->maxiobuflen)] = ((float *)fillvalue)[nv];
			}
		    }
		}
		else if (vsize==8)
		{
		    for (int nv=0;nv < nvars; nv++)
		    {
			for (int i=0;i<iodesc->maxiobuflen;i++)
			{
			    ((double *)vdesc0->iobuf)[i+nv*(iodesc->maxiobuflen)] = ((double *)fillvalue)[nv];
			}
		    }
		}
	    }
	}
	
	ierr = rearrange_comp2io(*ios, iodesc, array, vdesc0->iobuf, nvars);
    }/*  this is wrong, need to think about it
	 else{
	 vdesc0->iobuf = array;
	 } */
    switch(file->iotype)
    {
    case PIO_IOTYPE_NETCDF4P:
    case PIO_IOTYPE_PNETCDF:
	ierr = pio_write_darray_multi_nc(file, nvars, vid,
					 iodesc->ndims, iodesc->basetype, iodesc->gsize,
					 iodesc->maxregions, iodesc->firstregion, iodesc->llen,
					 iodesc->maxiobuflen, iodesc->num_aiotasks,
					 vdesc0->iobuf, frame);
	break;
    case PIO_IOTYPE_NETCDF4C:
    case PIO_IOTYPE_NETCDF:
	ierr = pio_write_darray_multi_nc_serial(file, nvars, vid,
						iodesc->ndims, iodesc->basetype, iodesc->gsize,
						iodesc->maxregions, iodesc->firstregion, iodesc->llen,
						iodesc->maxiobuflen, iodesc->num_aiotasks,
						vdesc0->iobuf, frame);
	if (vdesc0->iobuf)
	{
	    brel(vdesc0->iobuf);
	    vdesc0->iobuf = NULL;
	}
	break;

    }

    if (iodesc->rearranger == PIO_REARR_SUBSET && iodesc->needsfill &&
       iodesc->holegridsize>0)
    {
	if (vdesc0->fillbuf)
	{
	    piodie("Attempt to overwrite existing buffer",__FILE__,__LINE__);
	}

	vdesc0->fillbuf = bget(iodesc->holegridsize*vsize*nvars);
	//printf("%s %d %x\n",__FILE__,__LINE__,vdesc0->fillbuf);
	if (vsize==4)
	{
	    for (int nv=0;nv<nvars;nv++)
	    {
		for (int i=0;i<iodesc->holegridsize;i++)
		{
		    ((float *) vdesc0->fillbuf)[i+nv*iodesc->holegridsize] = ((float *) fillvalue)[nv];
		}
	    }
	}
	else if (vsize==8)
	{
	    for (int nv=0;nv<nvars;nv++)
	    {
		for (int i=0;i<iodesc->holegridsize;i++)
		{
		    ((double *) vdesc0->fillbuf)[i+nv*iodesc->holegridsize] = ((double *) fillvalue)[nv];
		}
	    }
	}
	switch(file->iotype)
	{
	case PIO_IOTYPE_PNETCDF:
	    ierr = pio_write_darray_multi_nc(file, nvars, vid,
					     iodesc->ndims, iodesc->basetype, iodesc->gsize,
					     iodesc->maxfillregions, iodesc->fillregion, iodesc->holegridsize,
					     iodesc->holegridsize, iodesc->num_aiotasks,
					     vdesc0->fillbuf, frame);
	    break;
	case PIO_IOTYPE_NETCDF4P:
	case PIO_IOTYPE_NETCDF4C:
	case PIO_IOTYPE_NETCDF:
	    /*       ierr = pio_write_darray_multi_nc_serial(file, nvars, vid,
		     iodesc->ndims, iodesc->basetype, iodesc->gsize,
		     iodesc->maxfillregions, iodesc->fillregion, iodesc->holegridsize,
		     iodesc->holegridsize, iodesc->num_aiotasks,
		     vdesc0->fillbuf, frame);
	    */
	    /*       if (vdesc0->fillbuf != NULL){
		     printf("%s %d %x\n",__FILE__,__LINE__,vdesc0->fillbuf);
		     brel(vdesc0->fillbuf);
		     vdesc0->fillbuf = NULL;
		     }
	    */
	    break;
	}
    }

    flush_output_buffer(file, flushtodisk, 0);

    return ierr;
}

/** Write a distributed array to the output file.
 *
 * This routine aggregates output on the compute nodes and only sends
 * it to the IO nodes when the compute buffer is full or when a flush
 * is triggered.
 *
 * @param ncid: the ncid of the open netCDF file.
 * @param vid: the variable ID returned by PIOc_def_var().
 * @param ioid: the I/O description ID as passed back by
 * PIOc_InitDecomp().
 * @param arraylen: the length of the array to be written. This
 * is the length of the distrubited array. That is, the length of
 * the portion of the data that is on the processor.
 * @param array: pointer to the data to be written. This is a
 * pointer to the distributed portion of the array that is on this
 * processor.
 * @param fillvalue: pointer to the fill value to be used for
 * missing data.
 *
 * @returns 0 for success, non-zero error code for failure.
 * @ingroup PIO_write_darray
 */
int PIOc_write_darray(const int ncid, const int vid, const int ioid,
		      const PIO_Offset arraylen, void *array, void *fillvalue)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;
    io_desc_t *iodesc;
    var_desc_t *vdesc;
    void *bufptr;
    size_t rlen;
    int ierr;
    MPI_Datatype vtype;
    wmulti_buffer *wmb;
    int tsize;
    int *tptr;
    void *bptr;
    void *fptr;
    bool recordvar;
    int needsflush;
    bufsize totfree, maxfree;

    ierr = PIO_NOERR;
    needsflush = 0; // false
    file = pio_get_file_from_id(ncid);
    if (file == NULL)
    {
	fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
	return PIO_EBADID;
    }
    if (! (file->mode & PIO_WRITE))
    {
	fprintf(stderr,"ERROR:  Attempt to write to read-only file\n");
	return PIO_EPERM;
    }

    iodesc = pio_get_iodesc_from_id(ioid);
    if (iodesc == NULL)
    {
	fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
	return PIO_EBADID;
    }
    ios = file->iosystem;

    vdesc = (file->varlist)+vid;
    if (vdesc == NULL)
	return PIO_EBADID;

    /* Is this a record variable? */
    recordvar = vdesc->record < 0 ? true : false;
    
    if (iodesc->ndof != arraylen)
    {
	fprintf(stderr,"ndof=%ld, arraylen=%ld\n",iodesc->ndof,arraylen);
	piodie("ndof != arraylen",__FILE__,__LINE__);
    }
    wmb = &(file->buffer);
    if (wmb->ioid == -1)
    {
	if (recordvar)
	    wmb->ioid = ioid;
	else
	    wmb->ioid = -(ioid);
    }
    else
    {
	// separate record and non-record variables
	if (recordvar)
	{
	    while(wmb->next && wmb->ioid!=ioid)
		if (wmb->next!=NULL)
		    wmb = wmb->next;
#ifdef _PNETCDF
	    /* flush the previous record before starting a new one. this is collective */
	    //       if (vdesc->request != NULL && (vdesc->request[0] != NC_REQ_NULL) ||
	    //	  (wmb->frame != NULL && vdesc->record != wmb->frame[0])){
	    //	 needsflush = 2;  // flush to disk
	    //  }
#endif
	}
	else
	{
	    while(wmb->next && wmb->ioid!= -(ioid))
	    {
		if (wmb->next!=NULL)
		    wmb = wmb->next;
	    }
	}
    }
    if ((recordvar && wmb->ioid != ioid) || (!recordvar && wmb->ioid != -(ioid)))
    {
	wmb->next = (wmulti_buffer *) bget((bufsize) sizeof(wmulti_buffer));
	if (wmb->next == NULL)
	    piomemerror(*ios,sizeof(wmulti_buffer), __FILE__,__LINE__);
	wmb=wmb->next;
	wmb->next=NULL;
	if (recordvar)
	    wmb->ioid = ioid;
	else
	    wmb->ioid = -(ioid);
	wmb->validvars=0;
	wmb->arraylen=arraylen;
	wmb->vid=NULL;
	wmb->data=NULL;
	wmb->frame=NULL;
	wmb->fillvalue=NULL;
    }

    MPI_Type_size(iodesc->basetype, &tsize);
    // At this point wmb should be pointing to a new or existing buffer
    // so we can add the data
    //     printf("%s %d %X %d %d %d\n",__FILE__,__LINE__,wmb->data,wmb->validvars,arraylen,tsize);
    //    cn_buffer_report(*ios, true);
    bfreespace(&totfree, &maxfree);
    if (needsflush == 0)
	needsflush = (maxfree <= 1.1*(1+wmb->validvars)*arraylen*tsize );
    MPI_Allreduce(MPI_IN_PLACE, &needsflush, 1,  MPI_INT,  MPI_MAX, ios->comp_comm);

    if (needsflush > 0 )
    {
	// need to flush first
	//      printf("%s %d %ld %d %ld %ld\n",__FILE__,__LINE__,maxfree, wmb->validvars, (1+wmb->validvars)*arraylen*tsize,totfree);
	cn_buffer_report(*ios, true);

	flush_buffer(ncid, wmb, needsflush == 2);  // if needsflush == 2 flush to disk otherwise just flush to io node
    }
    
    if (arraylen > 0)
	if (!(wmb->data = bgetr(wmb->data, (1+wmb->validvars)*arraylen*tsize)))
	    piomemerror(*ios, (1+wmb->validvars)*arraylen*tsize, __FILE__, __LINE__);
    
    if (!(wmb->vid = (int *) bgetr(wmb->vid,sizeof(int)*(1+wmb->validvars))))
	piomemerror(*ios, (1+wmb->validvars)*sizeof(int), __FILE__, __LINE__);
    
    if (vdesc->record >= 0)
	if (!(wmb->frame = (int *)bgetr(wmb->frame, sizeof(int) * (1 + wmb->validvars))))
	    piomemerror(*ios, (1+wmb->validvars)*sizeof(int), __FILE__, __LINE__);

    if (iodesc->needsfill)
	if (!(wmb->fillvalue = bgetr(wmb->fillvalue,tsize*(1+wmb->validvars))))
	    piomemerror(*ios, (1+wmb->validvars)*tsize  , __FILE__,__LINE__);

    if (iodesc->needsfill)
    {
	if (fillvalue)
	{
	    memcpy((char *) wmb->fillvalue+tsize*wmb->validvars,fillvalue, tsize);
	}
	else
	{
	    vtype = (MPI_Datatype) iodesc->basetype;
	    if (vtype == MPI_INTEGER)
	    {
		int fill = PIO_FILL_INT;
		memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize);
	    }
	    else if (vtype == MPI_FLOAT || vtype == MPI_REAL4)
	    {
		float fill = PIO_FILL_FLOAT;
		memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize);
	    }
	    else if (vtype == MPI_DOUBLE || vtype == MPI_REAL8)
	    {
		double fill = PIO_FILL_DOUBLE;
		memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize);
	    }
	    else if (vtype == MPI_CHARACTER)
	    {
		char fill = PIO_FILL_CHAR;
		memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize);
	    }
	    else
	    {
		fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",vtype);
	    }
	}

    }

    wmb->arraylen = arraylen;
    wmb->vid[wmb->validvars]=vid;
    bufptr = (void *)((char *) wmb->data + arraylen*tsize*wmb->validvars);
    if (arraylen>0)
	memcpy(bufptr, array, arraylen*tsize);
    /*
      if (tsize==8){
      double asum=0.0;
      printf("%s %d %d %d %d\n",__FILE__,__LINE__,vid,arraylen,iodesc->ndof);
      for (int k=0;k<arraylen;k++){
      asum += ((double *) array)[k];
      }
      printf("%s %d %d %g\n",__FILE__,__LINE__,vid,asum);
      }
    */

    //   printf("%s %d %d %d %d %X\n",__FILE__,__LINE__,wmb->validvars,wmb->ioid,vid,bufptr);

    if (wmb->frame!=NULL)
	wmb->frame[wmb->validvars]=vdesc->record;
    wmb->validvars++;

    //   printf("%s %d %d %d %d %d\n",__FILE__,__LINE__,wmb->validvars,iodesc->maxbytes/tsize, iodesc->ndof, iodesc->llen);
    if (wmb->validvars >= iodesc->maxbytes/tsize)
	PIOc_sync(ncid);

    return ierr;
}

/** Read an array of data from a file to the (parallel) IO library.
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param iodesc a pointer to the defined iodescriptor for the buffer
 * @param vid the variable id to be read
 * @param IOBUF the buffer to be read into from this mpi task
 *
 * @return 0 on success, error code otherwise.
 * @ingroup PIO_read_darray
 */
int pio_read_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid,
		       void *IOBUF)
{
    int ierr=PIO_NOERR;
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;
    int ndims, fndims;
    MPI_Status status;
    int i;

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:read_darray_nc");
#endif
    ios = file->iosystem;
    if (ios == NULL)
	return PIO_EBADID;

    vdesc = (file->varlist)+vid;

    if (vdesc == NULL)
	return PIO_EBADID;

    ndims = iodesc->ndims;
    ierr = PIOc_inq_varndims(file->fh, vid, &fndims);

    if (fndims==ndims)
	vdesc->record=-1;

    if (ios->ioproc)
    {
	io_region *region;
	size_t start[fndims];
	size_t count[fndims];
	size_t tmp_start[fndims];
	size_t tmp_count[fndims];
	size_t tmp_bufsize=1;
	int regioncnt;
	void *bufptr;
	int tsize;

	int rrlen=0;
	PIO_Offset *startlist[iodesc->maxregions];
	PIO_Offset *countlist[iodesc->maxregions];

	// buffer is incremented by byte and loffset is in terms of the iodessc->basetype
	// so we need to multiply by the size of the basetype
	// We can potentially allow for one iodesc to have multiple datatypes by allowing the
	// calling program to change the basetype.
	region = iodesc->firstregion;
	MPI_Type_size(iodesc->basetype, &tsize);
	if (fndims>ndims)
	{
	    ndims++;
	    if (vdesc->record<0)
		vdesc->record=0;
	}
	for (regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++)
	{
	    //         printf("%s %d %d %ld %d %d\n",__FILE__,__LINE__,regioncnt,region,fndims,ndims);
	    tmp_bufsize=1;
	    if (region==NULL || iodesc->llen==0)
	    {
		for (i=0;i<fndims;i++)
		{
		    start[i] = 0;
		    count[i] = 0;
		}
		bufptr=NULL;
	    }
	    else
	    {
		if (regioncnt==0 || region==NULL)
		    bufptr = IOBUF;
		else
		    bufptr=(void *)((char *) IOBUF + tsize*region->loffset);

		//		printf("%s %d %d %d %d\n",__FILE__,__LINE__,iodesc->llen - region->loffset, iodesc->llen, region->loffset);

		if (vdesc->record >= 0 && fndims>1)
		{
		    start[0] = vdesc->record;
		    for (i=1;i<ndims;i++)
		    {
			start[i] = region->start[i-1];
			count[i] = region->count[i-1];
			//	    printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,i,start[i],count[i]);
		    }
		    if (count[1] > 0)
			count[0] = 1;
		}
		else
		{
		    // Non-time dependent array
		    for (i=0;i<ndims;i++)
		    {
			start[i] = region->start[i];
			count[i] = region->count[i];
			// printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,i,start[i],count[i]);
		    }
		}
	    }

	    switch(file->iotype)
	    {
#ifdef _NETCDF4
	    case PIO_IOTYPE_NETCDF4P:
		if (iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8)
		{
		    ierr = nc_get_vara_double (file->fh, vid,start,count, bufptr);
		}
		else if (iodesc->basetype == MPI_INTEGER)
		{
		    ierr = nc_get_vara_int (file->fh, vid, start, count,  bufptr);
		}
		else if (iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4)
		{
		    ierr = nc_get_vara_float (file->fh, vid, start,  count,  bufptr);
		}
		else
		{
		    fprintf(stderr,"Type not recognized %d in pioc_read_darray\n",(int) iodesc->basetype);
		}
		break;
#endif
#ifdef _PNETCDF
	    case PIO_IOTYPE_PNETCDF:
	    {
		tmp_bufsize=1;
		for (int j = 0; j < fndims; j++)
		    tmp_bufsize *= count[j];

		if (tmp_bufsize > 0)
		{
		    startlist[rrlen] = (PIO_Offset *) bget(fndims * sizeof(PIO_Offset));
		    countlist[rrlen] = (PIO_Offset *) bget(fndims * sizeof(PIO_Offset));

		    for (int j = 0; j < fndims; j++)
		    {
			startlist[rrlen][j] = start[j];
			countlist[rrlen][j] = count[j];
			/*	      printf("%s %d %d %d %d %ld %ld %ld\n",__FILE__,__LINE__,realregioncnt,
				      iodesc->maxregions, j,start[j],count[j],tmp_bufsize);*/
		    }
		    rrlen++;
		}
		if (regioncnt==iodesc->maxregions-1)
		{
		    ierr = ncmpi_get_varn_all(file->fh, vid, rrlen, startlist,
					      countlist, IOBUF, iodesc->llen, iodesc->basetype);
		    for (i=0;i<rrlen;i++)
		    {
			brel(startlist[i]);
			brel(countlist[i]);
		    }
		}
	    }
	    break;
#endif
	    default:
		ierr = iotype_error(file->iotype,__FILE__,__LINE__);

	    }
	    if (region)
		region = region->next;
	} // for (regioncnt=0;...)
    }

    ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:read_darray_nc");
#endif

    return ierr;
}

/** Read an array of data from a file to the (serial) IO library.
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param iodesc a pointer to the defined iodescriptor for the buffer
 * @param vid the variable id to be read.
 * @param IOBUF the buffer to be read into from this mpi task
 *
 * @returns
 * @ingroup PIO_read_darray
 */
int pio_read_darray_nc_serial(file_desc_t *file, io_desc_t *iodesc,
			      const int vid, void *IOBUF)
{
    int ierr=PIO_NOERR;
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    var_desc_t *vdesc;
    int ndims, fndims;
    MPI_Status status;
    int i;

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:read_darray_nc_serial");
#endif
    ios = file->iosystem;
    if (ios == NULL)
	return PIO_EBADID;

    vdesc = (file->varlist)+vid;

    if (vdesc == NULL)
	return PIO_EBADID;

    ndims = iodesc->ndims;
    ierr = PIOc_inq_varndims(file->fh, vid, &fndims);

    if (fndims==ndims)
	vdesc->record=-1;

    if (ios->ioproc)
    {
	io_region *region;
	size_t start[fndims];
	size_t count[fndims];
	size_t tmp_start[fndims * iodesc->maxregions];
	size_t tmp_count[fndims * iodesc->maxregions];
	size_t tmp_bufsize;
	int regioncnt;
	void *bufptr;
	int tsize;

	int rrlen = 0;

	// buffer is incremented by byte and loffset is in terms of the iodessc->basetype
	// so we need to multiply by the size of the basetype
	// We can potentially allow for one iodesc to have multiple datatypes by allowing the
	// calling program to change the basetype.
	region = iodesc->firstregion;
	MPI_Type_size(iodesc->basetype, &tsize);
	if (fndims>ndims)
	{
	    if (vdesc->record < 0)
		vdesc->record = 0;
	}
	for (regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++)
	{
	    if (region==NULL || iodesc->llen==0)
	    {
		for (i = 0; i < fndims; i++)
		{
		    tmp_start[i + regioncnt * fndims] = 0;
		    tmp_count[i + regioncnt * fndims] = 0;
		}
		bufptr=NULL;
	    }
	    else
	    {
		if (vdesc->record >= 0 && fndims>1)
		{
		    tmp_start[regioncnt*fndims] = vdesc->record;
		    for (i=1;i<fndims;i++)
		    {
			tmp_start[i+regioncnt*fndims] = region->start[i-1];
			tmp_count[i+regioncnt*fndims] = region->count[i-1];
		    }
		    if (tmp_count[1 + regioncnt * fndims] > 0)
			tmp_count[regioncnt * fndims] = 1;
		}
		else
		{
		    // Non-time dependent array
		    for (i = 0; i < fndims; i++)
		    {
			tmp_start[i + regioncnt * fndims] = region->start[i];
			tmp_count[i + regioncnt * fndims] = region->count[i];
		    }
		}
		/*	for (i=0;i<fndims;i++){
			printf("%d %d %d %ld %ld\n",__LINE__,regioncnt,i,tmp_start[i+regioncnt*fndims],
			tmp_count[i+regioncnt*fndims]);
			}*/

	    }
	    if (region)
		region = region->next;
	} // for (regioncnt=0;...)

	if (ios->io_rank>0)
	{
	    MPI_Send(&(iodesc->llen), 1, MPI_OFFSET, 0, ios->io_rank, ios->io_comm);
	    if (iodesc->llen > 0)
	    {
		MPI_Send(&(iodesc->maxregions), 1, MPI_INT, 0,
			  ios->num_iotasks + ios->io_rank, ios->io_comm);
		MPI_Send(tmp_count, iodesc->maxregions*fndims, MPI_OFFSET, 0,
			  2 * ios->num_iotasks + ios->io_rank, ios->io_comm);
		MPI_Send(tmp_start, iodesc->maxregions*fndims, MPI_OFFSET, 0,
			  3 * ios->num_iotasks + ios->io_rank, ios->io_comm);
		MPI_Recv(IOBUF, iodesc->llen, iodesc->basetype, 0,
			 4 * ios->num_iotasks+ios->io_rank, ios->io_comm, &status);
	    }
	}
	else if (ios->io_rank == 0)
	{
	    int maxregions=0;
	    size_t loffset, regionsize;
	    size_t this_start[fndims*iodesc->maxregions];
	    size_t this_count[fndims*iodesc->maxregions];
	    //      for (i=ios->num_iotasks-1; i>=0; i--){
	    for (int rtask = 1; rtask <= ios->num_iotasks; rtask++)
	    {
		if (rtask<ios->num_iotasks)
		{
		    MPI_Recv(&tmp_bufsize, 1, MPI_OFFSET, rtask, rtask, ios->io_comm, &status);
		    if (tmp_bufsize>0)
		    {
			MPI_Recv(&maxregions, 1, MPI_INT, rtask, ios->num_iotasks+rtask,
				 ios->io_comm, &status);
			MPI_Recv(this_count, maxregions*fndims, MPI_OFFSET, rtask,
				 2 * ios->num_iotasks + rtask, ios->io_comm, &status);
			MPI_Recv(this_start, maxregions*fndims, MPI_OFFSET, rtask,
				 3 * ios->num_iotasks + rtask, ios->io_comm, &status);
		    }
		}
		else
		{
		    maxregions=iodesc->maxregions;
		    tmp_bufsize=iodesc->llen;
		}
		loffset = 0;
		for (regioncnt=0;regioncnt<maxregions;regioncnt++)
		{
		    bufptr=(void *)((char *) IOBUF + tsize*loffset);
		    regionsize=1;
		    if (rtask<ios->num_iotasks)
		    {
			for (int m=0; m<fndims; m++)
			{
			    start[m] = this_start[m+regioncnt*fndims];
			    count[m] = this_count[m+regioncnt*fndims];
			    regionsize*=count[m];
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
		    loffset+=regionsize;

#ifdef _NETCDF
		    // Cant use switch here because MPI_DATATYPE may not be simple (openmpi)
		    if (iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8)
		    {
			ierr = nc_get_vara_double (file->fh, vid,start, count, bufptr);
		    }
		    else if (iodesc->basetype == MPI_INTEGER)
		    {
			ierr = nc_get_vara_int (file->fh, vid, start, count,  bufptr);
		    }
		    else if (iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4)
		    {
			ierr = nc_get_vara_float (file->fh, vid, start, count,  bufptr);
		    }
		    else
		    {
			fprintf(stderr,"Type not recognized %d in pioc_write_darray_nc_serial\n",
				(int)iodesc->basetype);
		    }

		    if (ierr != PIO_NOERR)
		    {
			for (int i = 0; i < fndims; i++)
			    fprintf(stderr,"vid %d dim %d start %ld count %ld err %d\n",
				    vid, i, start[i], count[i], ierr);

		    }

#endif
		}
		if (rtask < ios->num_iotasks)
		    MPI_Send(IOBUF, tmp_bufsize, iodesc->basetype, rtask,
			     4 * ios->num_iotasks + rtask, ios->io_comm);
	    }
	}
    }

    ierr = check_netcdf(file, ierr, __FILE__, __LINE__);

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:read_darray_nc_serial");
#endif

    return ierr;
}

/** Read a field from a file to the IO library.
 * @ingroup PIO_read_darray
 *
 * @param ncid identifies the netCDF file
 * @param vid the variable ID to be read
 * @param ioid: the I/O description ID as passed back by
 * PIOc_InitDecomp().
 * @param arraylen: the length of the array to be read. This
 * is the length of the distrubited array. That is, the length of
 * the portion of the data that is on the processor.
 * @param array: pointer to the data to be read. This is a
 * pointer to the distributed portion of the array that is on this
 * processor.
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_read_darray
 */
int PIOc_read_darray(const int ncid, const int vid, const int ioid,
		     const PIO_Offset arraylen, void *array)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;
    io_desc_t *iodesc;
    void *iobuf=NULL;
    size_t rlen=0;
    int ierr, tsize;
    MPI_Datatype vtype;

    file = pio_get_file_from_id(ncid);

    if (file == NULL)
    {
	fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
	return PIO_EBADID;
    }
    iodesc = pio_get_iodesc_from_id(ioid);
    if (iodesc == NULL)
    {
	fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
	return PIO_EBADID;
    }
    ios = file->iosystem;
    if (ios->iomaster)
    {
	rlen = iodesc->maxiobuflen;
    }
    else
    {
	rlen = iodesc->llen;
    }

    if (iodesc->rearranger > 0)
    {
	if (ios->ioproc && rlen>0)
	{
	    MPI_Type_size(iodesc->basetype, &tsize);
	    iobuf = bget(((size_t) tsize)*rlen);
	    if (iobuf==NULL)
	    {
		piomemerror(*ios,rlen*((size_t) tsize), __FILE__,__LINE__);
	    }
	}
    }
    else
    {
	iobuf = array;
    }

    switch(file->iotype)
    {
    case PIO_IOTYPE_NETCDF:
    case PIO_IOTYPE_NETCDF4C:
	ierr = pio_read_darray_nc_serial(file, iodesc, vid, iobuf);
	break;
    case PIO_IOTYPE_PNETCDF:
    case PIO_IOTYPE_NETCDF4P:
	ierr = pio_read_darray_nc(file, iodesc, vid, iobuf);
	break;
    default:
	ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
    if (iodesc->rearranger > 0)
    {
	ierr = rearrange_io2comp(*ios, iodesc, iobuf, array);

	if (rlen>0)
	    brel(iobuf);
    }

    return ierr;

}

/** Flush the output buffer. This is only relevant for files opened
 * with pnetcdf.
 *
 * @param file a pointer to the open file descriptor for the file
 * that will be written to
 * @param force true to force the flushing of the buffer
 * @param addsize additional size to add to buffer (in bytes)
 *
 * @return 0 for success, error code otherwise.
 * @private
 * @ingroup PIO_write_darray
 */
int flush_output_buffer(file_desc_t *file, bool force, PIO_Offset addsize)
{
    int ierr = PIO_NOERR;

#ifdef _PNETCDF
    var_desc_t *vdesc;
    int *status;
    PIO_Offset usage = 0;

#ifdef TIMING
    /* Start timing this function. */
    GPTLstart("PIO:flush_output_buffer");
#endif

    pioassert(file != NULL, "file pointer not defined", __FILE__,
	      __LINE__);

    /* Find out the buffer usage. */
    ierr = ncmpi_inq_buffer_usage(file->fh, &usage);

    /* If we are not forcing a flush, spread the usage to all IO
     * tasks. */
    if (!force && file->iosystem->io_comm != MPI_COMM_NULL)
    {
	usage += addsize;
	MPI_Allreduce(MPI_IN_PLACE, &usage, 1,  MPI_OFFSET,  MPI_MAX,
		      file->iosystem->io_comm);
    }

    /* Keep track of the maximum usage. */
    if (usage > maxusage)
	maxusage = usage;

    /* If the user forces it, or the buffer has exceeded the size
     * limit, then flush to disk. */
    if (force || usage >= PIO_BUFFER_SIZE_LIMIT)
    {
	int rcnt;
	bool prev_dist=false;
	int prev_record=-1;
	int prev_type=0;
	int  maxreq;
	int reqcnt;
	maxreq = 0;
	reqcnt=0;
	rcnt=0;
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
	    if (rcnt>0 && (prev_record != vdesc->record || vdesc->nreqs==0))
	    {
		ierr = ncmpi_wait_all(file->fh, rcnt,  request,status);
		rcnt=0;
	    }
	    prev_record = vdesc->record;
#endif
	    //      printf("%s %d %d %d %d \n",__FILE__,__LINE__,i,vdesc->nreqs,vdesc->request);
	    for (reqcnt=0;reqcnt<vdesc->nreqs;reqcnt++)
	    {
		request[rcnt++] = max(vdesc->request[reqcnt],NC_REQ_NULL);
	    }
	    free(vdesc->request);
	    vdesc->request = NULL;
	    vdesc->nreqs = 0;
	    //      if (file->iosystem->io_rank < 2) printf("%s %d varid=%d\n",__FILE__,__LINE__,i);
#ifdef FLUSH_EVERY_VAR
	    ierr = ncmpi_wait_all(file->fh, rcnt, request, status);
	    rcnt = 0;
#endif
	}
	//    if (file->iosystem->io_rank==0){
	//   printf("%s %d %d\n",__FILE__,__LINE__,rcnt);
	// }
	if (rcnt > 0)
	{
	    /*
	      if (file->iosystem->io_rank==0){
	      printf("%s %d %d ",__FILE__,__LINE__,rcnt);
	      for (int i=0; i<rcnt; i++){
	      printf("%d ",request[i]);
	      }
	      printf("\n");
	      }*/
	    ierr = ncmpi_wait_all(file->fh, rcnt, request, status);
	}
	for (int i = 0; i < PIO_MAX_VARS; i++)
	{
	    vdesc = file->varlist + i;
	    if (vdesc->iobuf)
	    {
		brel(vdesc->iobuf);
		vdesc->iobuf=NULL;
	    }
	    if (vdesc->fillbuf)
	    {
		brel(vdesc->fillbuf);
		vdesc->fillbuf=NULL;
	    }
	}

    }

#ifdef TIMING
    /* Stop timing this function. */
    GPTLstop("PIO:flush_output_buffer");
#endif

#endif /* _PNETCDF */
    return ierr;
}

/** Print out info about the buffer for debug purposes.
 *
 * @param ios the IO system structure
 * @param collective true if collective report is desired
 *
 * @private
 * @ingroup PIO_write_darray
 */
void cn_buffer_report(iosystem_desc_t ios, bool collective)
{

    if (CN_bpool)
    {
	long bget_stats[5];
	long bget_mins[5];
	long bget_maxs[5];

	bstats(bget_stats, bget_stats+1,bget_stats+2,bget_stats+3,bget_stats+4);
	if (collective)
	{
	    MPI_Reduce(bget_stats, bget_maxs, 5, MPI_LONG, MPI_MAX, 0, ios.comp_comm);
	    MPI_Reduce(bget_stats, bget_mins, 5, MPI_LONG, MPI_MIN, 0, ios.comp_comm);
	    if (ios.compmaster)
	    {
		printf("PIO: Currently allocated buffer space %ld %ld\n",
		       bget_mins[0], bget_maxs[0]);
		printf("PIO: Currently available buffer space %ld %ld\n",
		       bget_mins[1], bget_maxs[1]);
		printf("PIO: Current largest free block %ld %ld\n",
		       bget_mins[2], bget_maxs[2]);
		printf("PIO: Number of successful bget calls %ld %ld\n",
		       bget_mins[3], bget_maxs[3]);
		printf("PIO: Number of successful brel calls  %ld %ld\n",
		       bget_mins[4], bget_maxs[4]);
		//	print_trace(stdout);
	    }
	}
	else
	{
	    printf("%d: PIO: Currently allocated buffer space %ld \n",
		   ios.union_rank, bget_stats[0]) ;
	    printf("%d: PIO: Currently available buffer space %ld \n",
		   ios.union_rank, bget_stats[1]);
	    printf("%d: PIO: Current largest free block %ld \n",
		   ios.union_rank, bget_stats[2]);
	    printf("%d: PIO: Number of successful bget calls %ld \n",
		   ios.union_rank, bget_stats[3]);
	    printf("%d: PIO: Number of successful brel calls  %ld \n",
		   ios.union_rank, bget_stats[4]);
	}
    }
}

/** Free the buffer pool. If malloc is used (that is, PIO_USE_MALLOC is
 * non zero), this function does nothing.
 *
 * @param ios the IO system structure
 *
 * @private
 * @ingroup PIO_write_darray
 */
void free_cn_buffer_pool(iosystem_desc_t ios)
{
#if !PIO_USE_MALLOC
    if (CN_bpool)
    {
	cn_buffer_report(ios, true);
	bpoolrelease(CN_bpool);
	//    free(CN_bpool);
	CN_bpool = NULL;
    }
#endif /* !PIO_USE_MALLOC */
}

/** Flush the buffer. 
 *
 * @param ncid identifies the netCDF file
 * @param wmb
 * @param flushtodisk
 *
 * @private
 * @ingroup PIO_write_darray
 */
void flush_buffer(int ncid, wmulti_buffer *wmb, bool flushtodisk)
{
    if (wmb->validvars > 0)
    {
	PIOc_write_darray_multi(ncid, wmb->vid,  wmb->ioid, wmb->validvars,
				wmb->arraylen, wmb->data, wmb->frame,
				wmb->fillvalue, flushtodisk);
	wmb->validvars = 0;
	brel(wmb->vid);
	wmb->vid = NULL;
	brel(wmb->data);
	wmb->data = NULL;
	if (wmb->fillvalue)
	    brel(wmb->fillvalue);
	if (wmb->frame)
	    brel(wmb->frame);
	wmb->fillvalue = NULL;
	wmb->frame = NULL;
    }
}

/** Compute the maximum aggregate number of bytes. 
 *
 * @param ios the IO system structure
 * @param iodesc a pointer to the defined iodescriptor for the buffer
 *
 * @private
 * @ingroup PIO_write_darray
 */
void compute_maxaggregate_bytes(const iosystem_desc_t ios, io_desc_t *iodesc)
{
    int maxbytesoniotask = INT_MAX;
    int maxbytesoncomputetask = INT_MAX;
    int maxbytes;

    // printf("%s %d %d %d\n",__FILE__,__LINE__,iodesc->maxiobuflen, iodesc->ndof);

    if (ios.ioproc && iodesc->maxiobuflen > 0)
	maxbytesoniotask = PIO_BUFFER_SIZE_LIMIT / iodesc->maxiobuflen;

    if (ios.comp_rank >= 0 && iodesc->ndof > 0)
	maxbytesoncomputetask = PIO_CNBUFFER_LIMIT / iodesc->ndof;

    maxbytes = min(maxbytesoniotask, maxbytesoncomputetask);

    //  printf("%s %d %d %d\n",__FILE__,__LINE__,maxbytesoniotask, maxbytesoncomputetask);

    MPI_Allreduce(MPI_IN_PLACE, &maxbytes, 1, MPI_INT, MPI_MIN, ios.union_comm);
    iodesc->maxbytes = maxbytes;
    //  printf("%s %d %d %d\n",__FILE__,__LINE__,iodesc->maxbytes,iodesc->maxiobuflen);

}
