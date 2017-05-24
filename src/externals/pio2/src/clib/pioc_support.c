/** @file
 * Support functions for the PIO library.
 */
#include <config.h>
#if PIO_ENABLE_LOGGING
#include <stdarg.h>
#include <unistd.h>
#endif /* PIO_ENABLE_LOGGING */
#include <pio.h>
#include <pio_internal.h>

#include <execinfo.h>

#define VERSNO 2001

/* Some logging constants. */
#if PIO_ENABLE_LOGGING
#define MAX_LOG_MSG 1024
#define MAX_RANK_STR 12
#define ERROR_PREFIX "ERROR: "
#define NC_LEVEL_DIFF 3
int pio_log_level = 0;
int pio_log_ref_cnt = 0;
int my_rank;
FILE *LOG_FILE = NULL;
#endif /* PIO_ENABLE_LOGGING */

/**
 * The PIO library maintains its own set of ncids. This is the next
 * ncid number that will be assigned.
 */
extern int pio_next_ncid;

/** The default error handler used when iosystem cannot be located. */
extern int default_error_handler;

/**
 * Return a string description of an error code. If zero is passed,
 * the errmsg will be "No error".
 *
 * @param pioerr the error code returned by a PIO function call.
 * @param errmsg Pointer that will get the error message. The message
 * will be PIO_MAX_NAME chars or less.
 * @return 0 on success.
 */
int PIOc_strerror(int pioerr, char *errmsg)
{
    LOG((1, "PIOc_strerror pioerr = %d", pioerr));

    /* Caller must provide this. */
    pioassert(errmsg, "pointer to errmsg string must be provided", __FILE__, __LINE__);

    /* System error? NetCDF and pNetCDF errors are always negative. */
    if (pioerr > 0)
    {
        const char *cp = (const char *)strerror(pioerr);
        if (cp)
            strncpy(errmsg, cp, PIO_MAX_NAME);
        else
            strcpy(errmsg, "Unknown Error");
    }
    else if (pioerr == PIO_NOERR)
        strcpy(errmsg, "No error");
#if defined(_NETCDF)
    else if (pioerr <= NC2_ERR && pioerr >= NC4_LAST_ERROR)     /* NetCDF error? */
        strncpy(errmsg, nc_strerror(pioerr), NC_MAX_NAME);
#endif /* endif defined(_NETCDF) */
#if defined(_PNETCDF)
    else if (pioerr > PIO_FIRST_ERROR_CODE)     /* pNetCDF error? */
        strncpy(errmsg, ncmpi_strerror(pioerr), NC_MAX_NAME);
#endif /* defined( _PNETCDF) */
    else
        /* Handle PIO errors. */
        switch(pioerr)
        {
        case PIO_EBADIOTYPE:
            strcpy(errmsg, "Bad IO type");
            break;
        default:
            strcpy(errmsg, "Unknown Error: Unrecognized error code");
        }

    return PIO_NOERR;
}

/**
 * Set the logging level if PIO was built with
 * PIO_ENABLE_LOGGING. Set to -1 for nothing, 0 for errors only, 1 for
 * important logging, and so on. Log levels below 1 are only printed
 * on the io/component root.
 *
 * A log file is also produced for each task. The file is called
 * pio_log_X.txt, where X is the (0-based) task number.
 *
 * If the library is not built with logging, this function does
 * nothing.
 *
 * @param level the logging level, 0 for errors only, 5 for max
 * verbosity.
 * @returns 0 on success, error code otherwise.
 */
int PIOc_set_log_level(int level)
{
    int ret;

#if PIO_ENABLE_LOGGING
    /* Set the log level. */
    pio_log_level = level;

#if NETCDF_C_LOGGING_ENABLED
    /* If netcdf logging is available turn it on starting at level = 4. */
    if (level > NC_LEVEL_DIFF)
        if ((ret = nc_set_log_level(level - NC_LEVEL_DIFF)))
            return pio_err(NULL, NULL, ret, __FILE__, __LINE__);
#endif /* NETCDF_C_LOGGING_ENABLED */
#endif /* PIO_ENABLE_LOGGING */

    return PIO_NOERR;
}

/**
 * Initialize logging.  Open log file, if not opened yet, or increment
 * ref count if already open.
 */
void pio_init_logging(void)
{
#if PIO_ENABLE_LOGGING
    char log_filename[NC_MAX_NAME];

    if (!LOG_FILE)
    {
        /* Create a filename with the rank in it. */
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        sprintf(log_filename, "pio_log_%d.txt", my_rank);

        /* Open a file for this rank to log messages. */
        LOG_FILE = fopen(log_filename, "w");

        pio_log_ref_cnt = 1;
    }
    else
    {
        pio_log_ref_cnt++;
    }
#endif /* PIO_ENABLE_LOGGING */
}

/**
 * Finalize logging - close log files, if open.
 */
void pio_finalize_logging(void )
{
#if PIO_ENABLE_LOGGING
    pio_log_ref_cnt -= 1;
    if (LOG_FILE)
    {
        if (pio_log_ref_cnt == 0)
        {
            fclose(LOG_FILE);
            LOG_FILE = NULL;
        }
        else
            LOG((2, "pio_finalize_logging, postpone close, ref_cnt = %d",
                 pio_log_ref_cnt));
    }
#endif /* PIO_ENABLE_LOGGING */
}

#if PIO_ENABLE_LOGGING
/**
 * This function prints out a message, if the severity of the message
 * is lower than the global pio_log_level. To use it, do something
 * like this:
 *
 * pio_log(0, "this computer will explode in %d seconds", i);
 *
 * After the first arg (the severity), use the rest like a normal
 * printf statement. Output will appear on stdout.
 * This function is heavily based on the function in section 15.5 of
 * the C FAQ.
 *
 * In code this functions should be wrapped in the LOG(()) macro.
 *
 * @param severity the severity of the message, 0 for error messages,
 * then increasing levels of verbosity.
 * @param fmt the format string.
 * @param ... the arguments used in format string.
 */
void pio_log(int severity, const char *fmt, ...)
{
    va_list argp;
    int t;
    int rem_len = MAX_LOG_MSG;
    char msg[MAX_LOG_MSG];
    char *ptr = msg;
    char rank_str[MAX_RANK_STR];

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
    {
        strncpy(ptr, ERROR_PREFIX, (rem_len > 0) ? rem_len : 0);
        ptr += strlen(ERROR_PREFIX);
        rem_len -= strlen(ERROR_PREFIX);
    }
    for (t = 0; t < severity; t++)
    {
        strncpy(ptr++, "\t", (rem_len > 0) ? rem_len : 0);
        rem_len--;
    }

    /* Show the rank. */
    snprintf(rank_str, MAX_RANK_STR, "%d ", my_rank);
    strncpy(ptr, rank_str, (rem_len > 0) ? rem_len : 0);
    ptr += strlen(rank_str);
    rem_len -= strlen(rank_str);

    /* Print out the variable list of args with vprintf. */
    va_start(argp, fmt);
    vsnprintf(ptr, ((rem_len > 0) ? rem_len : 0), fmt, argp);
    va_end(argp);

    /* Put on a final linefeed. */
    ptr = msg + strlen(msg);
    rem_len = MAX_LOG_MSG - strlen(msg);
    strncpy(ptr, "\n\0", (rem_len > 0) ? rem_len : 0);

    /* Send message to stdout. */
    fprintf(stdout, "%s", msg);

    /* Send message to log file. */
    if (LOG_FILE)
        fprintf(LOG_FILE, "%s", msg);

    /* Ensure an immediate flush of stdout. */
    fflush(stdout);
    if (LOG_FILE)
        fflush(LOG_FILE);
}
#endif /* PIO_ENABLE_LOGGING */

/**
 * Obtain a backtrace and print it to stderr. This is appended to the
 * text decomposition file.
 *
 * Note from Jim:
 *
 * The stack trace can be used to identify the usage in
 * the model code of the particular decomposition in question and so
 * if using the pio performance tool leads to tuning that could be
 * applied in the model you know more or less where to do it.
 *
 * It's also useful if you have a model bug - then you have 20 or so
 * decomp files and you need to identify the one that was problematic.
 * So it's used as an add to the developer and not used at all by any
 * automated process or tools.
 *
 * @param fp file pointer to send output to
 */
void print_trace(FILE *fp)
{
    void *array[10];
    size_t size;
    char **strings;
    size_t i;

    /* Note that this won't actually work. */
    if (fp == NULL)
        fp = stderr;

    size = backtrace(array, 10);
    strings = backtrace_symbols(array, size);

    fprintf(fp,"Obtained %zd stack frames.\n", size);

    for (i = 0; i < size; i++)
        fprintf(fp,"%s\n", strings[i]);

    free(strings);
}

/**
 * Exit due to lack of memory.
 *
 * @param ios the iosystem description struct
 * @param req amount of memory that was being requested
 * @param fname name of code file where error occured
 * @param line the line of code where the error occurred.
 */
void piomemerror(iosystem_desc_t *ios, size_t req, char *fname, int line)
{
    char msg[80];
    sprintf(msg, "out of memory requesting: %ld", req);
    cn_buffer_report(ios, false);
    piodie(msg, fname, line);
}

/**
 * Abort program and call MPI_Abort().
 *
 * @param msg an error message
 * @param fname name of code file where error occured
 * @param line the line of code where the error occurred.
 */
void piodie(const char *msg, const char *fname, int line)
{
    fprintf(stderr,"Abort with message %s in file %s at line %d\n",
            msg ? msg : "_", fname ? fname : "_", line);

    print_trace(stderr);
#ifdef MPI_SERIAL
    abort();
#else
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
}

/**
 * Perform an assert. Note that this function does nothing if NDEBUG
 * is defined.
 *
 * @param expression the expression to be evaluated
 * @param msg an error message
 * @param fname name of code file where error occured
 * @param line the line of code where the error occurred.
 */
void pioassert(_Bool expression, const char *msg, const char *fname, int line)
{
#ifndef NDEBUG
    if (!expression)
        piodie(msg, fname, line);
#endif
}

/**
 * Handle MPI errors. An error message is sent to stderr, then the
 * check_netcdf() function is called with PIO_EIO.
 *
 * @param file pointer to the file_desc_t info. Ignored if NULL.
 * @param mpierr the MPI return code to handle
 * @param filename the name of the code file where error occured.
 * @param line the line of code where error occured.
 * @return PIO_NOERR for no error, otherwise PIO_EIO.
 */
int check_mpi(file_desc_t *file, int mpierr, const char *filename,
              int line)
{
    return check_mpi2(NULL, file, mpierr, filename, line);
}

/**
 * Handle MPI errors. An error message is sent to stderr, then the
 * check_netcdf() function is called with PIO_EIO. This version of the
 * function accepts an ios parameter, for the (rare) occasions where
 * we have an ios but not a file.
 *
 * @param ios pointer to the iosystem_info_t. May be NULL.
 * @param file pointer to the file_desc_t info. May be NULL.
 * @param mpierr the MPI return code to handle
 * @param filename the name of the code file where error occured.
 * @param line the line of code where error occured.
 * @return PIO_NOERR for no error, otherwise PIO_EIO.
 */
int check_mpi2(iosystem_desc_t *ios, file_desc_t *file, int mpierr,
               const char *filename, int line)
{
    if (mpierr)
    {
        char errstring[MPI_MAX_ERROR_STRING];
        int errstrlen;

        /* If we can get an error string from MPI, print it to stderr. */
        if (!MPI_Error_string(mpierr, errstring, &errstrlen))
            fprintf(stderr, "MPI ERROR: %s in file %s at line %d\n",
                    errstring, filename ? filename : "_", line);

        /* Handle all MPI errors as PIO_EIO. */
        return pio_err(ios, file, PIO_EIO, filename, line);
    }
    return PIO_NOERR;
}

/**
 * Check the result of a netCDF API call.
 *
 * @param file pointer to the PIO structure describing this
 * file. Ignored if NULL.
 * @param status the return value from the netCDF call.
 * @param fname the name of the code file.
 * @param line the line number of the netCDF call in the code.
 * @return the error code
 */
int check_netcdf(file_desc_t *file, int status, const char *fname, int line)
{
    return check_netcdf2(NULL, file, status, fname, line);
}

/**
 * Check the result of a netCDF API call. This is the same as
 * check_netcdf() except for the extra iosystem_desc_t pointer, which
 * is used to determine error handling when there is no file_desc_t
 * pointer.
 *
 * @param ios pointer to the iosystem description struct. Ignored if NULL.
 * @param file pointer to the PIO structure describing this file. Ignored if NULL.
 * @param status the return value from the netCDF call.
 * @param fname the name of the code file.
 * @param line the line number of the netCDF call in the code.
 * @return the error code
 */
int check_netcdf2(iosystem_desc_t *ios, file_desc_t *file, int status,
                  const char *fname, int line)
{
    int eh = default_error_handler; /* Error handler that will be used. */
    char errmsg[PIO_MAX_NAME + 1];  /* Error message. */

    /* User must provide this. */
    pioassert(fname, "code file name must be provided", __FILE__, __LINE__);

    /* No harm, no foul. */
    if (status == PIO_NOERR)
        return PIO_NOERR;

    LOG((1, "check_netcdf2 status = %d fname = %s line = %d", status, fname, line));

    /* Pick an error handler. */
    if (ios)
        eh = ios->error_handler;
    if (file)
        eh = file->iosystem->error_handler;
    pioassert(eh == PIO_INTERNAL_ERROR || eh == PIO_BCAST_ERROR || eh == PIO_RETURN_ERROR,
              "invalid error handler", __FILE__, __LINE__);
    LOG((2, "check_netcdf2 chose error handler = %d", eh));

    /* Get an error message. */
    if (eh != PIO_BCAST_ERROR && !PIOc_strerror(status, errmsg))
    {
        fprintf(stderr, "%s\n", errmsg);
        LOG((1, "check_netcdf2 errmsg = %s", errmsg));
    }

    /* Decide what to do based on the error handler. */
    if (eh == PIO_INTERNAL_ERROR)
        piodie(errmsg, fname, line);        /* Die! */
    else if (eh == PIO_BCAST_ERROR)
    {
	if (ios)
	    MPI_Bcast(&status, 1, MPI_INT, ios->ioroot, ios->my_comm);
	else if (file)
	    MPI_Bcast(&status, 1, MPI_INT, file->iosystem->ioroot, file->iosystem->my_comm);
    }

    /* For PIO_RETURN_ERROR, just return the error. */
    return status;
}

/**
 * Handle an error in PIO. This will consult the error handler
 * settings and either call MPI_Abort() or return an error code.
 *
 * The error hanlder has three settings:
 *
 * Errors cause abort: PIO_INTERNAL_ERROR.
 *
 * Error codes are broadcast to all tasks: PIO_BCAST_ERROR.
 *
 * Errors are returned to caller with no internal action:
 * PIO_RETURN_ERROR.
 *
 * @param ios pointer to the IO system info. Ignored if NULL.
 * @param file pointer to the file description data. Ignored if
 * NULL. If provided file->iosystem is used as ios pointer.
 * @param err_num the error code
 * @param fname name of code file where error occured.
 * @param line the line of code where the error occurred.
 * @returns err_num if abort is not called.
 */
int pio_err(iosystem_desc_t *ios, file_desc_t *file, int err_num, const char *fname,
            int line)
{
    char err_msg[PIO_MAX_NAME + 1];
    int err_handler = default_error_handler; /* Default error handler. */
    int ret;

    /* User must provide this. */
    pioassert(fname, "file name must be provided", __FILE__, __LINE__);

    /* No harm, no foul. */
    if (err_num == PIO_NOERR)
        return PIO_NOERR;

    /* Get the error message. */
    if ((ret = PIOc_strerror(err_num, err_msg)))
        return ret;

    /* If logging is in use, log an error message. */
    LOG((0, "%s err_num = %d fname = %s line = %d", err_msg, err_num, fname ? fname : '\0', line));

    /* What error handler should we use? */
    if (file)
        err_handler = file->iosystem->error_handler;
    else if (ios)
        err_handler = ios->error_handler;

    LOG((2, "pio_err chose error handler = %d", err_handler));

    /* Should we abort? */
    if (err_handler == PIO_INTERNAL_ERROR)
    {
        /* For debugging only, this will print a traceback of the call tree.  */
        print_trace(stderr);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /* What should we do here??? */
    if (err_handler == PIO_BCAST_ERROR)
    {
        /* ??? */
    }

    /* If abort was not called, we'll get here. */
    return err_num;
}

/**
 * Allocate an region.
 *
 * ndims the number of dimensions for the data in this region.
 * @returns a pointer to the newly allocated io_region struct.
 */
io_region *alloc_region(int ndims)
{
    io_region *region;

    /* Allocate memory for the io_region struct. */
    if (!(region = bget(sizeof(io_region))))
        return NULL;

    /* Allocate memory for the array of start indicies. */
    if (!(region->start = bget(ndims * sizeof(PIO_Offset))))
    {
        brel(region);
        return NULL;
    }

    /* Allocate memory for the array of counts. */
    if (!(region->count = bget(ndims * sizeof(PIO_Offset))))
    {
        brel(region);
        brel(region->start);
        return NULL;
    }

    region->loffset = 0;
    region->next = NULL;

    /* Initialize start and count arrays to zero. */
    for (int i = 0; i < ndims; i++)
    {
        region->start[i] = 0;
        region->count[i] = 0;
    }

    return region;
}

/**
 * Find the MPI type for a PIO type.
 *
 * @param pio_type a PIO type, PIO_INT, PIO_FLOAT, etc.
 * @param mpi_type a pointer to MPI_Datatype that will get the MPI
 * type that coresponds to the PIO type.
 * @returns 0 for success, error code otherwise.
 */
int find_mpi_type(int pio_type, MPI_Datatype *mpi_type)
{
    /* Check input. */
    pioassert(mpi_type, "invalid input", __FILE__, __LINE__);

    /* Decide on the base type. */
    switch(pio_type)
    {
    case PIO_BYTE:
        *mpi_type = MPI_BYTE;
        break;
    case PIO_CHAR:
        *mpi_type = MPI_CHAR;
        break;
    case PIO_SHORT:
        *mpi_type = MPI_SHORT;
        break;
    case PIO_INT:
        *mpi_type = MPI_INT;
        break;
    case PIO_FLOAT:
        *mpi_type = MPI_FLOAT;
        break;
    case PIO_DOUBLE:
        *mpi_type = MPI_DOUBLE;
        break;
#ifdef _NETCDF4
    case PIO_UBYTE:
        *mpi_type = MPI_UNSIGNED_CHAR;
        break;
    case PIO_USHORT:
        *mpi_type = MPI_UNSIGNED_SHORT;
        break;
    case PIO_UINT:
        *mpi_type = MPI_UNSIGNED;
        break;
    case PIO_INT64:
        *mpi_type = MPI_LONG_LONG;
        break;
    case PIO_UINT64:
        *mpi_type = MPI_UNSIGNED_LONG_LONG;
        break;
    case PIO_STRING:
        *mpi_type = MPI_CHAR;
        break;
#endif /* _NETCDF4 */
    default:
        return PIO_EBADTYPE;
    }

    return PIO_NOERR;
}

/**
 * Allocate space for an IO description struct.
 *
 * @param ios pointer to the IO system info.
 * @param piotype the PIO data type (ex. PIO_FLOAT, PIO_INT, etc.).
 * @param ndims the number of dimensions.
 * @iodesc pointer that gets a pointer to the newly allocated
 * io_desc_t or NULL if allocation failed.
 * @returns 0 for success, error code otherwise.
 */
int malloc_iodesc(iosystem_desc_t *ios, int piotype, int ndims,
                  io_desc_t **iodesc)
{
    MPI_Datatype mpi_type;
    int ret;

    /* Check input. */
    pioassert(ios && iodesc, "invalid input", __FILE__, __LINE__);

    /* Allocate space for the io_desc_t struct. */
    if (!(*iodesc = calloc(1, sizeof(io_desc_t))))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Get the MPI type corresponding with the PIO type. */
    if ((ret = find_mpi_type(piotype, &mpi_type)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Decide on the base type. */
    (*iodesc)->basetype = mpi_type;

    /* Initialize some values in the struct. */
    (*iodesc)->maxregions = 1;
    (*iodesc)->ioid = -1;
    (*iodesc)->ndims = ndims;
    (*iodesc)->firstregion = alloc_region(ndims);

    /* Set the swap memory settings to defaults. */
    (*iodesc)->rearr_opts = ios->rearr_opts;

    return PIO_NOERR;
}

/**
 * Free a region list.
 *
 * top a pointer to the start of the list to free.
 */
void free_region_list(io_region *top)
{
    io_region *ptr, *tptr;

    ptr = top;
    while (ptr)
    {
        if (ptr->start)
            brel(ptr->start);
        if (ptr->count)
            brel(ptr->count);
        tptr = ptr;
        ptr = ptr->next;
        brel(tptr);
    }
}

/**
 * Free a decomposition map.
 *
 * @param iosysid the IO system ID.
 * @param ioid the ID of the decomposition map to free.
 * @returns 0 for success, error code otherwise.
 */
int PIOc_freedecomp(int iosysid, int ioid)
{
    iosystem_desc_t *ios;
    io_desc_t *iodesc;
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function calls. */

    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    if (!(iodesc = pio_get_iodesc_from_id(ioid)))
        return pio_err(ios, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_FREEDECOMP; /* Message for async notification. */

            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&iosysid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&ioid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            LOG((2, "PIOc_freedecomp iosysid = %d ioid = %d", iosysid, ioid));
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(NULL, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    }

    /* Free the map. */
    free(iodesc->map);

    /* Free the dimlens. */
    free(iodesc->dimlen);

    if (iodesc->rfrom)
        free(iodesc->rfrom);

    if (iodesc->rtype)
    {
        for (int i = 0; i < iodesc->nrecvs; i++)
            if (iodesc->rtype[i] != PIO_DATATYPE_NULL)
                if ((mpierr = MPI_Type_free(iodesc->rtype + i)))
                    return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

        free(iodesc->rtype);
    }

    if (iodesc->stype)
    {
        for (int i = 0; i < iodesc->num_stypes; i++)
            if (iodesc->stype[i] != PIO_DATATYPE_NULL)
                if ((mpierr = MPI_Type_free(iodesc->stype + i)))
                    return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

        iodesc->num_stypes = 0;
        free(iodesc->stype);
    }

    if (iodesc->scount)
        free(iodesc->scount);

    if (iodesc->rcount)
        free(iodesc->rcount);

    if (iodesc->sindex)
        free(iodesc->sindex);

    if (iodesc->rindex)
        free(iodesc->rindex);

    if (iodesc->firstregion)
        free_region_list(iodesc->firstregion);

    if (iodesc->fillregion)
        free_region_list(iodesc->fillregion);

    if (iodesc->rearranger == PIO_REARR_SUBSET)
        if ((mpierr = MPI_Comm_free(&iodesc->subset_comm)))
            return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    return pio_delete_iodesc_from_list(ioid);
}

/**
 * Read a decomposition map from a file. The decomp file is only read
 * by task 0 in the communicator.
 *
 * @param file the filename
 * @param ndims pointer to an int with the number of dims.
 * @param gdims pointer to an array of dimension ids.
 * @param fmaplen
 * @param map
 * @param comm
 * @returns 0 for success, error code otherwise.
 */
int PIOc_readmap(const char *file, int *ndims, int **gdims, PIO_Offset *fmaplen,
                 PIO_Offset **map, MPI_Comm comm)
{
    int npes, myrank;
    int rnpes, rversno;
    int j;
    int *tdims;
    PIO_Offset *tmap;
    MPI_Status status;
    PIO_Offset maplen;
    int mpierr; /* Return code for MPI calls. */

    /* Check inputs. */
    if (!file || !ndims || !gdims || !fmaplen || !map)
        return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

    if ((mpierr = MPI_Comm_size(comm, &npes)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &myrank)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);

    if (myrank == 0)
    {
        FILE *fp = fopen(file, "r");
        if (!fp)
            pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

        fscanf(fp,"version %d npes %d ndims %d\n",&rversno, &rnpes, ndims);

        if (rversno != VERSNO)
            return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

        if (rnpes < 1 || rnpes > npes)
            return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

        if ((mpierr = MPI_Bcast(&rnpes, 1, MPI_INT, 0, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(ndims, 1, MPI_INT, 0, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        if (!(tdims = calloc(*ndims, sizeof(int))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        for (int i = 0; i < *ndims; i++)
            fscanf(fp,"%d ", tdims + i);

        if ((mpierr = MPI_Bcast(tdims, *ndims, MPI_INT, 0, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);

        for (int i = 0; i < rnpes; i++)
        {
            fscanf(fp, "%d %lld", &j, &maplen);
            if (j != i)  // Not sure how this could be possible
                return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);
            if (!(tmap = malloc(maplen * sizeof(PIO_Offset))))
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
            for (j = 0; j < maplen; j++)
                fscanf(fp, "%lld ", tmap+j);

            if (i > 0)
            {
                if ((mpierr = MPI_Send(&maplen, 1, PIO_OFFSET, i, i + npes, comm)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(tmap, maplen, PIO_OFFSET, i, i, comm)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                free(tmap);
            }
            else
            {
                *map = tmap;
                *fmaplen = maplen;
            }
        }
        fclose(fp);
    }
    else
    {
        if ((mpierr = MPI_Bcast(&rnpes, 1, MPI_INT, 0, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(ndims, 1, MPI_INT, 0, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        if (!(tdims = calloc(*ndims, sizeof(int))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(tdims, *ndims, MPI_INT, 0, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);

        if (myrank < rnpes)
        {
            if ((mpierr = MPI_Recv(&maplen, 1, PIO_OFFSET, 0, myrank + npes, comm, &status)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            if (!(tmap = malloc(maplen * sizeof(PIO_Offset))))
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
            if ((mpierr = MPI_Recv(tmap, maplen, PIO_OFFSET, 0, myrank, comm, &status)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            *map = tmap;
        }
        else
        {
            tmap = NULL;
            maplen = 0;
        }
        *fmaplen = maplen;
    }
    *gdims = tdims;
    return PIO_NOERR;
}

/**
 * Read a decomposition map from file.
 *
 * @param file the filename
 * @param ndims pointer to the number of dimensions
 * @param gdims pointer to an array of dimension ids
 * @param maplen pointer to the length of the map
 * @param map pointer to the map array
 * @param f90_comm
 * @returns 0 for success, error code otherwise.
 */
int PIOc_readmap_from_f90(const char *file, int *ndims, int **gdims, PIO_Offset *maplen,
                          PIO_Offset **map, int f90_comm)
{
    return PIOc_readmap(file, ndims, gdims, maplen, map, MPI_Comm_f2c(f90_comm));
}

/**
 * Write the decomposition map to a file using netCDF, everyones
 * favorite data format.
 *
 * @param iosysid the IO system ID.
 * @param filename the filename to be used.
 * @param cmode for PIOc_create(). Will be bitwise or'd with NC_WRITE.
 * @param ioid the ID of the IO description.
 * @param comm an MPI communicator.
 * @param title optial title attribute for the file. Must be less than
 * NC_MAX_NAME + 1 if provided. Ignored if NULL.
 * @param history optial history attribute for the file. Must be less
 * than NC_MAX_NAME + 1 if provided. Ignored if NULL.
 * @param fortran_order set to non-zero if fortran array ordering is
 * used, or to zero if C array ordering is used.
 * @returns 0 for success, error code otherwise.
 */
int PIOc_write_nc_decomp(int iosysid, const char *filename, int cmode, int ioid,
                         MPI_Comm comm, char *title, char *history, int fortran_order)
{
    iosystem_desc_t *ios;
    io_desc_t *iodesc;
    int npes;   /* Size of this communicator. */
    int myrank; /* Rank of this task. */
    int mpierr;
    int ret;

    /* Get the IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Check inputs. */
    if (!filename)
        return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((1, "PIOc_write_nc_decomp filename = %s iosysid = %d ioid = %d", filename,
         iosysid, ioid));

    /* Get the IO desc, which describes the decomposition. */
    if (!(iodesc = pio_get_iodesc_from_id(ioid)))
        return pio_err(ios, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Get the communicator size and task rank. */
    if ((mpierr = MPI_Comm_size(comm, &npes)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &myrank)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "npes = %d myrank = %d", npes, myrank));

    /* Allocate memory for the nmaplen. On task 0, this will contain
     * the length of the map on each task, for all tasks. */
    int task_maplen[npes];

    /* Gather maplens from all tasks and fill the task_maplen array on
     * all tasks. */
    if ((mpierr = MPI_Allgather(&iodesc->maplen, 1, MPI_INT, task_maplen, 1, MPI_INT, comm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    /* We will need to know the maximum maplen used for any task. */
    int max_maplen;

    /* Find the max maxplen. */
    if ((mpierr = MPI_Allreduce(&iodesc->maplen, &max_maplen, 1, MPI_INT, MPI_MAX, comm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((3, "max_maplen = %d", max_maplen));

    /* 2D array that, on task 0, will hold all the map information for
     * all tasks. */
    int full_map[npes][max_maplen];

    /* Fill local array with my map. Use the fill value for unused */
    /* elements at the end if max_maplen is longer than maplen. Also
     * subtract 1 because the iodesc->map is 1-based. */
    int my_map[max_maplen];
    for (int e = 0; e < max_maplen; e++)
    {
        my_map[e] = e < iodesc->maplen ? iodesc->map[e] - 1 : NC_FILL_INT;
        LOG((3, "my_map[%d] = %d", e, my_map[e]));
    }

    /* Gather my_map from all tasks and fill the full_map array. */
    if ((mpierr = MPI_Allgather(&my_map, max_maplen, MPI_INT, full_map, max_maplen,
                             MPI_INT, comm)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);

    for (int p = 0; p < npes; p++)
        for (int e = 0; e < max_maplen; e++)
            LOG((3, "full_map[%d][%d] = %d", p, e, full_map[p][e]));

    /* Write the netCDF decomp file. */
    if ((ret = pioc_write_nc_decomp_int(iosysid, filename, cmode, iodesc->ndims, iodesc->dimlen, npes,
                                        task_maplen, (int *)full_map, title, history, fortran_order)))
        return ret;


    return PIO_NOERR;
}

/**
 * Read the decomposition map from a netCDF decomp file produced by
 * PIOc_write_nc_decomp().
 *
 * @param iosysid the IO system ID.
 * @param filename the name of the decomp file.
 * @param ioid pointer that will get the newly-assigned ID of the IO
 * description. The ioid is needed to later free the decomposition.
 * @param comm an MPI communicator.
 * @param pio_type the PIO type to be used as the type for the data.
 * @param title pointer that will get optial title attribute for the
 * file. Will be less than NC_MAX_NAME + 1 if provided. Ignored if
 * NULL.
 * @param history pointer that will get optial history attribute for
 * the file. Will be less than NC_MAX_NAME + 1 if provided. Ignored if
 * NULL.
 * @param fortran_order pointer that gets set to 1 if fortran array
 * ordering is used, or to zero if C array ordering is used.
 * @returns 0 for success, error code otherwise.
 */
int PIOc_read_nc_decomp(int iosysid, const char *filename, int *ioidp, MPI_Comm comm,
                        int pio_type, char *title, char *history, int *fortran_order)
{
    iosystem_desc_t *ios; /* Pointer to the IO system info. */
    int ndims;            /* The number of data dims (except unlim). */
    int max_maplen;       /* The max maplen of any task. */
    int *global_dimlen;   /* An array with sizes of global dimensions. */
    int *task_maplen;     /* A map of one tasks mapping to global data. */
    int *full_map;        /* A map with the task maps of every task. */
    int num_tasks_decomp; /* The number of tasks for this decomp. */
    int size;             /* Size of comm. */
    int my_rank;          /* Task rank in comm. */
    MPI_Datatype mpi_type;             /* Will be used as the basetype in iodesc. */
    int mpi_type_int;                  /* int version of mpi_type. */
    char source_in[PIO_MAX_NAME + 1];  /* Text metadata in decomp file. */
    char version_in[PIO_MAX_NAME + 1]; /* Text metadata in decomp file. */
    int mpierr;
    int ret;

    /* Get the IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Check inputs. */
    if (!filename || !ioidp)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((1, "PIOc_read_nc_decomp filename = %s iosysid = %d pio_type = %d",
         filename, iosysid, pio_type));

    /* Get the MPI type. We need it as an int. */
    if ((ret = find_mpi_type(pio_type, &mpi_type)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    mpi_type_int = mpi_type;
    LOG((2, "mpi_type = %d mpi_type_int = %d", mpi_type, mpi_type_int));

    /* Get the communicator size and task rank. */
    if ((mpierr = MPI_Comm_size(comm, &size)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &my_rank)))
        return check_mpi2(ios, NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "size = %d my_rank = %d", size, my_rank));

    /* Read the file. This allocates three arrays that we have to
     * free. */
    if ((ret = pioc_read_nc_decomp_int(iosysid, filename, &ndims, &global_dimlen, &num_tasks_decomp,
                                       &task_maplen, &max_maplen, &full_map, title, history,
                                       source_in, version_in, fortran_order)))
        return ret;
    LOG((2, "ndims = %d num_tasks_decomp = %d max_maplen = %d", ndims, num_tasks_decomp,
         max_maplen));

    /* If the size does not match the number of tasks in the decomp,
     * that's an error. */
    if (size != num_tasks_decomp)
        ret = PIO_EINVAL;

    /* Now initialize the iodesc on each task for this decomposition. */
    if (!ret)
    {
        PIO_Offset compmap[task_maplen[my_rank]];

        /* Copy array into PIO_Offset array. Make it 1 based. */
        for (int e = 0; e < task_maplen[my_rank]; e++)
            compmap[e] = full_map[my_rank * max_maplen + e] + 1;

        /* Initialize the decomposition. */
        ret = PIOc_InitDecomp(iosysid, pio_type, ndims, global_dimlen, task_maplen[my_rank],
                              compmap, ioidp, NULL, NULL, NULL);
    }

    /* Free resources. */
    free(global_dimlen);
    free(task_maplen);
    free(full_map);

    return ret;
}

/* Write the decomp information in netCDF. This is an internal
 * function.
 *
 * @param iosysid the IO system ID.
 * @param filename the name the decomp file will have.
 * @param cmode for PIOc_create(). Will be bitwise or'd with NC_WRITE.
 * @param ndims number of dims in the data being described.
 * @param global_dimlen an array, of size ndims, with the size of the
 * global array in each dimension.
 * @param num_tasks the number of tasks the data are decomposed over.
 * @param task_maplen array of size num_tasks with the length of the
 * map for each task.
 * @param map pointer to a 2D array of size [num_tasks][max_maplen]
 * with the 0-based mapping from local to global array elements.
 * @param title null-terminated string that will be written as an
 * attribute. If provided, length must be < PIO_MAX_NAME + 1. Ignored
 * if NULL.
 * @param history null-terminated string that will be written as an
 * attribute. If provided, length must be < PIO_MAX_NAME + 1. Ignored
 * if NULL.
 * @param fortran_order set to non-zero if using fortran array
 * ordering, 0 for C array ordering.
 * @returns 0 for success, error code otherwise.
 */
int pioc_write_nc_decomp_int(int iosysid, const char *filename, int cmode, int ndims,
                             int *global_dimlen, int num_tasks, int *task_maplen, int *map,
                             const char *title, const char *history, int fortran_order)
{
    iosystem_desc_t *ios;
    int max_maplen = 0;
    int ncid;
    int ret;

    /* Get the IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Check inputs. */
    if (!filename || !global_dimlen || !task_maplen)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);
    if (title)
        if (strlen(title) > PIO_MAX_NAME)
            return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);
    if (history)
        if (strlen(history) > PIO_MAX_NAME)
            return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((2, "pioc_write_nc_decomp_int iosysid = %d filename = %s ndims = %d num_tasks = %d",
         iosysid, filename, ndims, num_tasks));

    /* Find the maximum maplen. */
    for (int t = 0; t < num_tasks; t++)
        if (task_maplen[t] > max_maplen)
            max_maplen = task_maplen[t];
    LOG((3, "max_maplen = %d", max_maplen));

    /* Create the netCDF decomp file. */
    if ((ret = PIOc_create(iosysid, filename, cmode | NC_WRITE, &ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Write an attribute with the version of this file. */
    char version[PIO_MAX_NAME + 1];
    sprintf(version, "%d.%d.%d", PIO_VERSION_MAJOR, PIO_VERSION_MINOR, PIO_VERSION_PATCH);
    if ((ret = PIOc_put_att_text(ncid, NC_GLOBAL, DECOMP_VERSION_ATT_NAME,
                                 strlen(version) + 1, version)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Write an attribute with the max map len. */
    if ((ret = PIOc_put_att_int(ncid, NC_GLOBAL, DECOMP_MAX_MAPLEN_ATT_NAME,
                                PIO_INT, 1, &max_maplen)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Write title attribute, if the user provided one. */
    if (title)
        if ((ret = PIOc_put_att_text(ncid, NC_GLOBAL, DECOMP_TITLE_ATT_NAME,
                                     strlen(title) + 1, title)))
            return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Write history attribute, if the user provided one. */
    if (history)
        if ((ret = PIOc_put_att_text(ncid, NC_GLOBAL, DECOMP_HISTORY_ATT_NAME,
                                     strlen(history) + 1, history)))
            return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Write a source attribute. */
    char source[] = "Decomposition file produced by PIO library.";
    if ((ret = PIOc_put_att_text(ncid, NC_GLOBAL, DECOMP_SOURCE_ATT_NAME,
                                 strlen(source) + 1, source)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Write an attribute with array ordering (C or Fortran). */
    char c_order_str[] = DECOMP_C_ORDER_STR;
    char fortran_order_str[] = DECOMP_FORTRAN_ORDER_STR;
    char *my_order_str = fortran_order ? fortran_order_str : c_order_str;
    if ((ret = PIOc_put_att_text(ncid, NC_GLOBAL, DECOMP_ORDER_ATT_NAME,
                                 strlen(my_order_str) + 1, my_order_str)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Write an attribute with the stack trace. This can be helpful
     * for debugging. */
    #define MAX_BACKTRACE 10
    void *bt[MAX_BACKTRACE];
    size_t bt_size;
    char **bt_strings;
    bt_size = backtrace(bt, MAX_BACKTRACE);
    bt_strings = backtrace_symbols(bt, bt_size);

    /* Find the max size. */
    int max_bt_size = 0;
    for (int b = 0; b < bt_size; b++)
        if (strlen(bt_strings[b]) > max_bt_size)
            max_bt_size = strlen(bt_strings[b]);
    if (max_bt_size > NC_MAX_NAME)
        max_bt_size = NC_MAX_NAME;

    /* Copy the backtrace into one long string. */
    char full_bt[max_bt_size * bt_size + bt_size + 1];
    full_bt[0] = '\0';
    for (int b = 0; b < bt_size; b++)
    {
        strncat(full_bt, bt_strings[b], max_bt_size);
        strcat(full_bt, "\n");
    }
    free(bt_strings);
    printf("full_bt = %s", full_bt);

    /* Write the stack trace as an attribute. */
    if ((ret = PIOc_put_att_text(ncid, NC_GLOBAL, DECOMP_BACKTRACE_ATT_NAME,
                                 strlen(full_bt) + 1, full_bt)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* We need a dimension for the dimensions in the data. (Example:
     * for 4D data we will need to store 4 dimension IDs.) */
    int dim_dimid;
    if ((ret = PIOc_def_dim(ncid, DECOMP_DIM_DIM, ndims, &dim_dimid)))
        return ret;

    /* We need a dimension for tasks. If we have 4 tasks, we need to
     * store an array of length 4 with the size of the local array on
     * each task. */
    int task_dimid;
    if ((ret = PIOc_def_dim(ncid, DECOMP_TASK_DIM_NAME, num_tasks, &task_dimid)))
        return ret;

    /* We need a dimension for the map. It's length may vary, we will
     * use the max_maplen for the dimension size. */
    int mapelem_dimid;
    if ((ret = PIOc_def_dim(ncid, DECOMP_MAPELEM_DIM_NAME, max_maplen,
                            &mapelem_dimid)))
        return ret;

    /* Define a var to hold the global size of the array for each
     * dimension. */
    int gsize_varid;
    if ((ret = PIOc_def_var(ncid, DECOMP_GLOBAL_SIZE_VAR_NAME, NC_INT, 1,
                            &dim_dimid, &gsize_varid)))
        return ret;

    /* Define a var to hold the length of the local array on each task. */
    int maplen_varid;
    if ((ret = PIOc_def_var(ncid, DECOMP_MAPLEN_VAR_NAME, NC_INT, 1, &task_dimid,
                            &maplen_varid)))
        return ret;

    /* Define a 2D var to hold the length of the local array on each task. */
    int map_varid;
    int map_dimids[2] = {task_dimid, mapelem_dimid};
    if ((ret = PIOc_def_var(ncid, DECOMP_MAP_VAR_NAME, NC_INT, 2, map_dimids,
                            &map_varid)))
        return ret;

    /* End define mode, to write data. */
    if ((ret = PIOc_enddef(ncid)))
        return ret;

    /* Write the global dimension sizes. */
    if ((PIOc_put_var_int(ncid, gsize_varid, global_dimlen)))
        return ret;

    /* Write the size of the local array on each task. */
    if ((PIOc_put_var_int(ncid, maplen_varid, task_maplen)))
        return ret;

    /* Write the map. */
    if ((PIOc_put_var_int(ncid, map_varid, map)))
        return ret;

    /* Close the netCDF decomp file. */
    if ((ret = PIOc_closefile(ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/* Read the decomp information from a netCDF decomp file. This is an
 * internal function.
 *
 * @param iosysid the IO system ID.
 * @param filename the name the decomp file will have.
 * @param ndims pointer to int that will get number of dims in the
 * data being described. Ignored if NULL.
 * @param global_dimlen a pointer that gets an array, of size ndims,
 * that will have the size of the global array in each
 * dimension. Ignored if NULL, otherwise must be freed by caller.
 * @param num_tasks pointer to int that gets the number of tasks the
 * data are decomposed over. Ignored if NULL.
 * @param task_maplen pointer that gets array of size num_tasks that
 * gets the length of the map for each task. Ignored if NULL,
 * otherwise must be freed by caller.
 * @param max_maplen pointer to int that gets the maximum maplen for
 * any task. Ignored if NULL.
 * @param map pointer that gets a 2D array of size [num_tasks][max_maplen]
 * that will have the 0-based mapping from local to global array
 * elements. Ignored if NULL, otherwise must be freed by caller.
 * @param title pointer that will get the contents of title attribute,
 * if present. If present, title will be < PIO_MAX_NAME + 1 in
 * length. Ignored if NULL.
 * @param history pointer that will get the contents of history attribute,
 * if present. If present, history will be < PIO_MAX_NAME + 1 in
 * length. Ignored if NULL.
 * @param source pointer that will get the contents of source
 * attribute. Source will be < PIO_MAX_NAME + 1 in length. Ignored if
 * NULL.
 * @param version pointer that will get the contents of version
 * attribute. It will be < PIO_MAX_NAME + 1 in length. Ignored if
 * NULL.
 * @param fortran_order int pointer that will get a 0 if this
 * decomposition file uses C array ordering, 1 if it uses Fortran
 * array ordering.
 * @returns 0 for success, error code otherwise.
 */
int pioc_read_nc_decomp_int(int iosysid, const char *filename, int *ndims, int **global_dimlen,
                            int *num_tasks, int **task_maplen, int *max_maplen, int **map, char *title,
                            char *history, char *source, char *version, int *fortran_order)
{
    iosystem_desc_t *ios;
    int ncid;
    int ret;

    /* Get the IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Check inputs. */
    if (!filename)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((1, "pioc_read_nc_decomp_int iosysid = %d filename = %s", iosysid, filename));

    /* Open the netCDF decomp file. */
    if ((ret = PIOc_open(iosysid, filename, NC_WRITE, &ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Read version attribute. */
    char version_in[PIO_MAX_NAME + 1];
    if ((ret = PIOc_get_att_text(ncid, NC_GLOBAL, DECOMP_VERSION_ATT_NAME, version_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    LOG((3, "version_in = %s", version_in));
    if (version)
        strncpy(version, version_in, PIO_MAX_NAME + 1);

    /* Read order attribute. */
    char order_in[PIO_MAX_NAME + 1];
    if ((ret = PIOc_get_att_text(ncid, NC_GLOBAL, DECOMP_ORDER_ATT_NAME, order_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    LOG((3, "order_in = %s", order_in));
    if (fortran_order)
    {
        if (!strncmp(order_in, DECOMP_C_ORDER_STR, PIO_MAX_NAME + 1))
            *fortran_order = 0;
        else if (!strncmp(order_in, DECOMP_FORTRAN_ORDER_STR, PIO_MAX_NAME + 1))
            *fortran_order = 1;
        else
            return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);
    }

    /* Read attribute with the max map len. */
    int max_maplen_in;
    if ((ret = PIOc_get_att_int(ncid, NC_GLOBAL, DECOMP_MAX_MAPLEN_ATT_NAME, &max_maplen_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    LOG((3, "max_maplen_in = %d", version_in));
    if (max_maplen)
        *max_maplen = max_maplen_in;

    /* Read title attribute, if it is in the file. */
    char title_in[NC_MAX_NAME + 1];
    ret = PIOc_get_att_text(ncid, NC_GLOBAL, DECOMP_TITLE_ATT_NAME, title_in);
    if (ret == PIO_NOERR)
    {
        /* If the caller wants it, copy the title for them. */
        if (title)
            strncpy(title, title_in, PIO_MAX_NAME + 1);
    }
    else if (ret == PIO_ENOTATT)
    {
        /* No title attribute. */
        if (title)
            title[0] = '\0';
    }
    else
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Read history attribute, if it is in the file. */
    char history_in[NC_MAX_NAME + 1];
    ret = PIOc_get_att_text(ncid, NC_GLOBAL, DECOMP_HISTORY_ATT_NAME, history_in);
    if (ret == PIO_NOERR)
    {
        /* If the caller wants it, copy the history for them. */
        if (history)
            strncpy(history, history_in, PIO_MAX_NAME + 1);
    }
    else if (ret == PIO_ENOTATT)
    {
        /* No history attribute. */
        if (history)
            history[0] = '\0';
    }
    else
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Read source attribute. */
    char source_in[NC_MAX_NAME + 1];
    if ((ret = PIOc_get_att_text(ncid, NC_GLOBAL, DECOMP_SOURCE_ATT_NAME, source_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if (source)
        strncpy(source, source_in, PIO_MAX_NAME + 1);

    /* Read dimension for the dimensions in the data. (Example: for 4D
     * data we will need to store 4 dimension IDs.) */
    int dim_dimid;
    PIO_Offset ndims_in;
    if ((ret = PIOc_inq_dimid(ncid, DECOMP_DIM_DIM, &dim_dimid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if ((ret = PIOc_inq_dim(ncid, dim_dimid, NULL, &ndims_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if (ndims)
        *ndims = ndims_in;

    /* Read the global sizes of the array. */
    int gsize_varid;
    int global_dimlen_in[ndims_in];
    if ((ret = PIOc_inq_varid(ncid, DECOMP_GLOBAL_SIZE_VAR_NAME, &gsize_varid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if ((ret = PIOc_get_var_int(ncid, gsize_varid, global_dimlen_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if (global_dimlen)
    {
        if (!(*global_dimlen = malloc(ndims_in * sizeof(int))))
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        for (int d = 0; d < ndims_in; d++)
            (*global_dimlen)[d] = global_dimlen_in[d];
    }

    /* Read dimension for tasks. If we have 4 tasks, we need to store
     * an array of length 4 with the size of the local array on each
     * task. */
    int task_dimid;
    PIO_Offset num_tasks_in;
    if ((ret = PIOc_inq_dimid(ncid, DECOMP_TASK_DIM_NAME, &task_dimid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if ((ret = PIOc_inq_dim(ncid, task_dimid, NULL, &num_tasks_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if (num_tasks)
        *num_tasks = num_tasks_in;

    /* Read the length if the local array on each task. */
    int maplen_varid;
    int task_maplen_in[num_tasks_in];
    if ((ret = PIOc_inq_varid(ncid, DECOMP_MAPLEN_VAR_NAME, &maplen_varid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if ((ret = PIOc_get_var_int(ncid, maplen_varid, task_maplen_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if (task_maplen)
    {
        if (!(*task_maplen = malloc(num_tasks_in * sizeof(int))))
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        for (int t = 0; t < num_tasks_in; t++)
            (*task_maplen)[t] = task_maplen_in[t];
    }

    /* Read the map. */
    int map_varid;
    int map_in[num_tasks_in][max_maplen_in];
    if ((ret = PIOc_inq_varid(ncid, DECOMP_MAP_VAR_NAME, &map_varid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if ((ret = PIOc_get_var_int(ncid, map_varid, (int *)map_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if (map)
    {
        if (!(*map = malloc(num_tasks_in * max_maplen_in * sizeof(int))))
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        for (int t = 0; t < num_tasks_in; t++)
            for (int l = 0; l < max_maplen_in; l++)
                (*map)[t * max_maplen_in + l] = map_in[t][l];
    }

    /* Close the netCDF decomp file. */
    if ((ret = PIOc_closefile(ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    return PIO_NOERR;
}

/**
 * Write the decomposition map to a file.
 *
 * @param file the filename to be used.
 * @param iosysid the IO system ID.
 * @param ioid the ID of the IO description.
 * @param comm an MPI communicator.
 * @returns 0 for success, error code otherwise.
 */
int PIOc_write_decomp(const char *file, int iosysid, int ioid, MPI_Comm comm)
{
    iosystem_desc_t *ios;
    io_desc_t *iodesc;

    LOG((1, "PIOc_write_decomp file = %s iosysid = %d ioid = %d", file, iosysid, ioid));

    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    if (!(iodesc = pio_get_iodesc_from_id(ioid)))
        return pio_err(ios, NULL, PIO_EBADID, __FILE__, __LINE__);

    return PIOc_writemap(file, iodesc->ndims, iodesc->dimlen, iodesc->maplen, iodesc->map,
                         comm);
}

/**
 * Write the decomposition map to a file.
 *
 * @param file the filename
 * @param ndims the number of dimensions
 * @param gdims an array of dimension ids
 * @param maplen the length of the map
 * @param map the map array
 * @param comm an MPI communicator.
 * @returns 0 for success, error code otherwise.
 */
int PIOc_writemap(const char *file, int ndims, const int *gdims, PIO_Offset maplen,
                  PIO_Offset *map, MPI_Comm comm)
{
    int npes, myrank;
    PIO_Offset *nmaplen = NULL;
    MPI_Status status;
    int i;
    PIO_Offset *nmap;
    int mpierr; /* Return code for MPI calls. */

    LOG((1, "PIOc_writemap file = %s ndims = %d maplen = %d", file, ndims, maplen));

    if ((mpierr = MPI_Comm_size(comm, &npes)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &myrank)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    LOG((2, "npes = %d myrank = %d", npes, myrank));

    /* Allocate memory for the nmaplen. */
    if (myrank == 0)
        if (!(nmaplen = malloc(npes * sizeof(PIO_Offset))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    if ((mpierr = MPI_Gather(&maplen, 1, PIO_OFFSET, nmaplen, 1, PIO_OFFSET, 0, comm)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);

    /* Only rank 0 writes the file. */
    if (myrank == 0)
    {
        FILE *fp;

        /* Open the file to write. */
        if (!(fp = fopen(file, "w")))
            return pio_err(NULL, NULL, PIO_EIO, __FILE__, __LINE__);

        /* Write the version and dimension info. */
        fprintf(fp,"version %d npes %d ndims %d \n", VERSNO, npes, ndims);
        for (i = 0; i < ndims; i++)
            fprintf(fp, "%d ", gdims[i]);
        fprintf(fp, "\n");

        /* Write the map. */
        fprintf(fp, "0 %lld\n", nmaplen[0]);
        for (i = 0; i < nmaplen[0]; i++)
            fprintf(fp, "%lld ", map[i]);
        fprintf(fp,"\n");

        for (i = 1; i < npes; i++)
        {
            LOG((2, "creating nmap for i = %d", i));
            nmap = (PIO_Offset *)malloc(nmaplen[i] * sizeof(PIO_Offset));

            if ((mpierr = MPI_Send(&i, 1, MPI_INT, i, npes + i, comm)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            if ((mpierr = MPI_Recv(nmap, nmaplen[i], PIO_OFFSET, i, i, comm, &status)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            LOG((2,"MPI_Recv map complete"));

            fprintf(fp, "%d %lld\n", i, nmaplen[i]);
            for (int j = 0; j < nmaplen[i]; j++)
                fprintf(fp, "%lld ", nmap[j]);
            fprintf(fp, "\n");

            free(nmap);
        }
        /* Free memory for the nmaplen. */
        free(nmaplen);
        fprintf(fp, "\n");
        print_trace(fp);

        /* Close the file. */
        fclose(fp);
        LOG((2,"decomp file closed."));
    }
    else
    {
        LOG((2,"ready to MPI_Recv..."));
        if ((mpierr = MPI_Recv(&i, 1, MPI_INT, 0, npes+myrank, comm, &status)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        LOG((2,"MPI_Recv got %d", i));
        if ((mpierr = MPI_Send(map, maplen, PIO_OFFSET, 0, myrank, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        LOG((2,"MPI_Send map complete"));
    }

    return PIO_NOERR;
}

/**
 * Write the decomposition map to a file for F90.
 *
 * @param file the filename
 * @param ndims the number of dimensions
 * @param gdims an array of dimension ids
 * @param maplen the length of the map
 * @param map the map array
 * @param comm an MPI communicator.
 * @returns 0 for success, error code otherwise.
 */
int PIOc_writemap_from_f90(const char *file, int ndims, const int *gdims,
                           PIO_Offset maplen, const PIO_Offset *map, int f90_comm)
{
    return PIOc_writemap(file, ndims, gdims, maplen, (PIO_Offset *)map,
                         MPI_Comm_f2c(f90_comm));
}

/**
 * Create a new file using pio. This is an internal function that is
 * called by both PIOc_create() and PIOc_createfile(). Input
 * parameters are read on comp task 0 and ignored elsewhere.
 *
 * @param iosysid A defined pio system ID, obtained from
 * PIOc_InitIntercomm() or PIOc_InitAsync().
 * @param ncidp A pointer that gets the ncid of the newly created
 * file.
 * @param iotype A pointer to a pio output format. Must be one of
 * PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF, PIO_IOTYPE_NETCDF4C, or
 * PIO_IOTYPE_NETCDF4P.
 * @param filename The filename to create.
 * @param mode The netcdf mode for the create operation.
 * @returns 0 for success, error code otherwise.
 * @ingroup PIO_createfile
 */
int PIOc_createfile_int(int iosysid, int *ncidp, int *iotype, const char *filename,
                        int mode)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
    int ierr;              /* Return code from function calls. */

    /* Get the IO system info from the iosysid. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* User must provide valid input for these parameters. */
    if (!ncidp || !iotype || !filename || strlen(filename) > NC_MAX_NAME)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* A valid iotype must be specified. */
    if (!iotype_is_valid(*iotype))
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((1, "PIOc_createfile iosysid = %d iotype = %d filename = %s mode = %d",
         iosysid, *iotype, filename, mode));

    /* Allocate space for the file info. */
    if (!(file = calloc(sizeof(file_desc_t), 1)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Fill in some file values. */
    file->fh = -1;
    file->iosystem = ios;
    file->iotype = *iotype;
    file->buffer.ioid = -1;
    for (int i = 0; i < PIO_MAX_VARS; i++)
    {
        file->varlist[i].record = -1;
        file->varlist[i].ndims = -1;
	file->varlist[i].iobuf = NULL;
	file->varlist[i].fillbuf = NULL;
    }
    file->mode = mode;

    /* Set to true if this task should participate in IO (only true for
     * one task with netcdf serial files. */
    if (file->iotype == PIO_IOTYPE_NETCDF4P || file->iotype == PIO_IOTYPE_PNETCDF ||
        ios->io_rank == 0)
        file->do_io = 1;

    LOG((2, "file->do_io = %d ios->async_interface = %d", file->do_io, ios->async_interface));

    /* If async is in use, and this is not an IO task, bcast the
     * parameters. */
    if (ios->async_interface)
    {
        int msg = PIO_MSG_CREATE_FILE;
        size_t len = strlen(filename);

        if (!ios->ioproc)
        {
            /* Send the message to the message handler. */
            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            /* Send the parameters of the function call. */
            if (!mpierr)
                mpierr = MPI_Bcast(&len, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&file->iotype, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&file->mode, 1, MPI_INT, ios->compmaster, ios->intercomm);
            LOG((2, "len = %d filename = %s iotype = %d mode = %d", len, filename,
                 file->iotype, file->mode));
        }

        /* Handle MPI errors. */
        LOG((2, "handling mpi errors mpierr = %d", mpierr));
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(file, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this task is in the IO component, do the IO. */
    if (ios->ioproc)
    {
        switch (file->iotype)
        {
#ifdef _NETCDF4
        case PIO_IOTYPE_NETCDF4P:
            file->mode = file->mode |  NC_MPIIO | NC_NETCDF4;
            LOG((2, "Calling nc_create_par io_comm = %d mode = %d fh = %d",
                 ios->io_comm, file->mode, file->fh));
            ierr = nc_create_par(filename, file->mode, ios->io_comm, ios->info, &file->fh);
            LOG((2, "nc_create_par returned %d file->fh = %d", ierr, file->fh));
            break;
        case PIO_IOTYPE_NETCDF4C:
            file->mode = file->mode | NC_NETCDF4;
#endif
        case PIO_IOTYPE_NETCDF:
            if (!ios->io_rank)
            {
                LOG((2, "Calling nc_create mode = %d", file->mode));
                ierr = nc_create(filename, file->mode, &file->fh);
            }
            break;
#ifdef _PNETCDF
        case PIO_IOTYPE_PNETCDF:
            LOG((2, "Calling ncmpi_create mode = %d", file->mode));
            ierr = ncmpi_create(ios->io_comm, filename, file->mode, ios->info, &file->fh);
            if (!ierr)
                ierr = ncmpi_buffer_attach(file->fh, pio_buffer_size_limit);
            break;
#endif
        }
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);

    /* If there was an error, free the memory we allocated and handle error. */
    if (ierr)
    {
        free(file);
        return check_netcdf2(ios, NULL, ierr, __FILE__, __LINE__);
    }

    /* Broadcast mode to all tasks. */
    if ((mpierr = MPI_Bcast(&file->mode, 1, MPI_INT, ios->ioroot, ios->union_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);

    /* This flag is implied by netcdf create functions but we need
       to know if its set. */
    file->mode = file->mode | PIO_WRITE;

    /* Assign the PIO ncid, necessary because files may be opened
     * on mutilple iosystems, causing the underlying library to
     * reuse ncids. Hilarious confusion ensues. */
    file->pio_ncid = pio_next_ncid++;
    LOG((2, "file->fh = %d file->pio_ncid = %d", file->fh, file->pio_ncid));

    /* Return the ncid to the caller. */
    *ncidp = file->pio_ncid;

    /* Add the struct with this files info to the global list of
     * open files. */
    pio_add_to_file_list(file);

    LOG((2, "Created file %s file->fh = %d file->pio_ncid = %d", filename,
         file->fh, file->pio_ncid));

    return ierr;
}

/**
 * Open an existing file using PIO library. This is an internal
 * function. Depending on the value of the retry parameter, a failed
 * open operation will be handled differently. If retry is non-zero,
 * then a failed attempt to open a file with netCDF-4 (serial or
 * parallel), or parallel-netcdf will be followed by an attempt to
 * open the file as a serial classic netCDF file. This is an important
 * feature to some NCAR users. The functionality is exposed to the
 * user as PIOc_openfile() (which does the retry), and PIOc_open()
 * (which does not do the retry).
 *
 * Input parameters are read on comp task 0 and ignored elsewhere.
 *
 * @param iosysid: A defined pio system descriptor (input)
 * @param ncidp: A pio file descriptor (output)
 * @param iotype: A pio output format (input)
 * @param filename: The filename to open
 * @param mode: The netcdf mode for the open operation
 * @param retry: non-zero to automatically retry with netCDF serial
 * classic.
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_openfile
 */
int PIOc_openfile_retry(int iosysid, int *ncidp, int *iotype, const char *filename,
                        int mode, int retry)
{
    iosystem_desc_t *ios;  /** Pointer to io system information. */
    file_desc_t *file;     /** Pointer to file information. */
    int imode;  /** internal mode val for netcdf4 file open */
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    int ierr = PIO_NOERR;  /** Return code from function calls. */

    /* Get the IO system info from the iosysid. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* User must provide valid input for these parameters. */
    if (!ncidp || !iotype || !filename)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);
    if (*iotype < PIO_IOTYPE_PNETCDF || *iotype > PIO_IOTYPE_NETCDF4P)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    LOG((2, "PIOc_openfile_retry iosysid = %d iotype = %d filename = %s mode = %d retry = %d",
         iosysid, *iotype, filename, mode, retry));

    /* Allocate space for the file info. */
    if (!(file = calloc(sizeof(*file), 1)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Fill in some file values. */
    file->fh = -1;
    file->iotype = *iotype;
    file->iosystem = ios;
    file->mode = mode;

    for (int i = 0; i < PIO_MAX_VARS; i++)
    {
        file->varlist[i].record = -1;
        file->varlist[i].ndims = -1;
    }

    /* Set to true if this task should participate in IO (only true
     * for one task with netcdf serial files. */
    if (file->iotype == PIO_IOTYPE_NETCDF4P || file->iotype == PIO_IOTYPE_PNETCDF ||
        ios->io_rank == 0)
        file->do_io = 1;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        int msg = PIO_MSG_OPEN_FILE;
        size_t len = strlen(filename);

        if (!ios->ioproc)
        {
            /* Send the message to the message handler. */
            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            /* Send the parameters of the function call. */
            if (!mpierr)
                mpierr = MPI_Bcast(&len, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast((void *)filename, len + 1, MPI_CHAR, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&file->iotype, 1, MPI_INT, ios->compmaster, ios->intercomm);
            if (!mpierr)
                mpierr = MPI_Bcast(&file->mode, 1, MPI_INT, ios->compmaster, ios->intercomm);
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
        switch (file->iotype)
        {
#ifdef _NETCDF4

        case PIO_IOTYPE_NETCDF4P:
#ifdef _MPISERIAL
            ierr = nc_open(filename, file->mode, &file->fh);
#else
            imode = file->mode |  NC_MPIIO;
            ierr = nc_open_par(filename, imode, ios->io_comm, ios->info, &file->fh);
            if (ierr == PIO_NOERR)
                file->mode = imode;
            LOG((2, "PIOc_openfile_retry:nc_open_par filename = %s mode = %d imode = %d ierr = %d",
                 filename, file->mode, imode, ierr));
#endif
            break;

        case PIO_IOTYPE_NETCDF4C:
            if (ios->io_rank == 0)
            {
                imode = file->mode | NC_NETCDF4;
                ierr = nc_open(filename, imode, &file->fh);
                if (ierr == PIO_NOERR)
                    file->mode = imode;
            }
            break;

#endif

        case PIO_IOTYPE_NETCDF:
            if (ios->io_rank == 0)
                ierr = nc_open(filename, file->mode, &file->fh);
            break;

#ifdef _PNETCDF
        case PIO_IOTYPE_PNETCDF:
            ierr = ncmpi_open(ios->io_comm, filename, file->mode, ios->info, &file->fh);

            // This should only be done with a file opened to append
            if (ierr == PIO_NOERR && (file->mode & PIO_WRITE))
            {
                if (ios->iomaster == MPI_ROOT)
                    LOG((2, "%d Setting IO buffer %ld", __LINE__, pio_buffer_size_limit));
                ierr = ncmpi_buffer_attach(file->fh, pio_buffer_size_limit);
            }
            LOG((2, "ncmpi_open(%s) : fd = %d", filename, file->fh));
            break;
#endif

        default:
            return pio_err(ios, file, PIO_EBADIOTYPE, __FILE__, __LINE__);
        }

        /* If the caller requested a retry, and we failed to open a
           file due to an incompatible type of NetCDF, try it once
           with just plain old basic NetCDF. */
        if (retry)
        {
            LOG((2, "retry error code ierr = %d io_rank %d", ierr, ios->io_rank));
            if ((ierr == NC_ENOTNC || ierr == NC_EINVAL) && (file->iotype != PIO_IOTYPE_NETCDF))
            {
                if (ios->iomaster == MPI_ROOT)
                    printf("PIO2 pio_file.c retry NETCDF\n");

                /* reset ierr on all tasks */
                ierr = PIO_NOERR;

                /* reset file markers for NETCDF on all tasks */
                file->iotype = PIO_IOTYPE_NETCDF;

                /* open netcdf file serially on main task */
                if (ios->io_rank == 0)
                    ierr = nc_open(filename, file->mode, &file->fh);
                else
                    file->do_io = 0;
            }
            LOG((2, "retry nc_open(%s) : fd = %d, iotype = %d, do_io = %d, ierr = %d", filename, file->fh, file->iotype, file->do_io, ierr));
        }
    }

    /* Broadcast and check the return code. */
    LOG((2, "Bcasting error code ierr = %d ios->ioroot = %d ios->my_comm = %d", ierr, ios->ioroot,
         ios->my_comm));
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);
    LOG((2, "Bcast error code ierr = %d", ierr));

    /* If there was an error, free allocated memory and deal with the error. */
    if (ierr)
    {
        free(file);
        return check_netcdf2(ios, NULL, ierr, __FILE__, __LINE__);
    }
    LOG((2, "error code Bcast complete ierr = %d ios->my_comm = %d", ierr, ios->my_comm));

    /* Broadcast results to all tasks. Ignore NULL parameters. */
    if ((mpierr = MPI_Bcast(&file->mode, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);

    /* Create the ncid that the user will see. This is necessary
     * because otherwise ncids will be reused if files are opened
     * on multiple iosystems. */
    file->pio_ncid = pio_next_ncid++;

    /* Return the PIO ncid to the user. */
    *ncidp = file->pio_ncid;

    /* Add this file to the list of currently open files. */
    pio_add_to_file_list(file);

    LOG((2, "Opened file %s file->pio_ncid = %d file->fh = %d ierr = %d",
         filename, file->pio_ncid, file->fh, ierr));

    return ierr;
}

/**
 * Internal function to provide inq_type function for pnetcdf.
 *
 * @param ncid ignored because pnetcdf does not have user-defined
 * types.
 * @param xtype type to learn about.
 * @param name pointer that gets name of type. Ignored if NULL.
 * @param sizep pointer that gets size of type. Ignored if NULL.
 * @returns 0 on success, error code otherwise.
 */
int pioc_pnetcdf_inq_type(int ncid, nc_type xtype, char *name,
                          PIO_Offset *sizep)
{
    int typelen;

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
    default:
        return PIO_EBADTYPE;
    }

    /* If pointers were supplied, copy results. */
    if (sizep)
        *sizep = typelen;
    if (name)
        strcpy(name, "some type");

    return PIO_NOERR;
}

/**
 * This is an internal function that handles both PIOc_enddef and
 * PIOc_redef.
 *
 * @param ncid the ncid of the file to enddef or redef
 * @param is_enddef set to non-zero for enddef, 0 for redef.
 * @returns PIO_NOERR on success, error code on failure. */
int pioc_change_def(int ncid, int is_enddef)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI functions. */

    LOG((2, "pioc_change_def ncid = %d is_enddef = %d", ncid, is_enddef));

    /* Find the info about this file. When I check the return code
     * here, some tests fail. ???*/
    if ((ierr = pio_get_file(ncid, &file)))
        return pio_err(NULL, NULL, ierr, __FILE__, __LINE__);
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async_interface)
    {
        if (!ios->ioproc)
        {
            int msg = is_enddef ? PIO_MSG_ENDDEF : PIO_MSG_REDEF;
            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&ncid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            LOG((3, "pioc_change_def ncid = %d mpierr = %d", ncid, mpierr));
        }

        /* Handle MPI errors. */
        LOG((3, "pioc_change_def handling MPI errors"));
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(file, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    LOG((3, "pioc_change_def ios->ioproc = %d", ios->ioproc));
    if (ios->ioproc)
    {
        LOG((3, "pioc_change_def calling netcdf function file->fh = %d file->do_io = %d iotype = %d",
             file->fh, file->do_io, file->iotype));
#ifdef _PNETCDF
        if (file->iotype == PIO_IOTYPE_PNETCDF)
        {
            if (is_enddef)
                ierr = ncmpi_enddef(file->fh);
            else
                ierr = ncmpi_redef(file->fh);
        }
#endif /* _PNETCDF */
        if (file->iotype != PIO_IOTYPE_PNETCDF && file->do_io)
        {
            if (is_enddef)
            {
                LOG((3, "pioc_change_def calling nc_enddef file->fh = %d", file->fh));
                ierr = nc_enddef(file->fh);
            }
            else
                ierr = nc_redef(file->fh);
        }
    }

    /* Broadcast and check the return code. */
    LOG((3, "pioc_change_def bcasting return code ierr = %d", ierr));
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(file, mpierr, __FILE__, __LINE__);
    if (ierr)
        return check_netcdf(file, ierr, __FILE__, __LINE__);
    LOG((3, "pioc_change_def succeeded"));

    return ierr;
}

/**
 * Check whether an IO type is valid for the build.
 *
 * @param iotype the IO type to check
 * @returns 0 if valid, non-zero otherwise.
 */
int iotype_is_valid(int iotype)
{
    /* Assume it's not valid. */
    int ret = 0;

    /* All builds include netCDF. */
    if (iotype == PIO_IOTYPE_NETCDF)
        ret++;

    /* Some builds include netCDF-4. */
#ifdef _NETCDF4
    if (iotype == PIO_IOTYPE_NETCDF4C || iotype == PIO_IOTYPE_NETCDF4P)
        ret++;
#endif /* _NETCDF4 */

    /* Some builds include pnetcdf. */
    if (iotype == PIO_IOTYPE_PNETCDF)
        ret++;
#ifdef _PNETCDF
#endif /* _PNETCDF */

    return ret;
}

/**
 * Internal function to compare rearranger flow control options.
 *
 * @param opt pointer to rearranger flow control options to compare.
 * @param exp_opt pointer to rearranger flow control options with
 * expected values.
 * @return true if values in opt == values in exp_opt, false
 * otherwise.
 */
bool cmp_rearr_comm_fc_opts(const rearr_comm_fc_opt_t *opt,
                            const rearr_comm_fc_opt_t *exp_opt)
{
    bool is_same = true;

    assert(opt && exp_opt);

    if (opt->enable_hs != exp_opt->enable_hs)
    {
        LOG((1, "Warning rearranger enable_hs = %s, expected = %s",
             opt->enable_hs ? "TRUE" : "FALSE", exp_opt->enable_hs ? "TRUE" : "FALSE"));
        is_same = false;
    }

    if (opt->enable_isend != exp_opt->enable_isend)
    {
        LOG((1, "Warning rearranger enable_isend = %s, expected = %s",
             opt->enable_isend ? "TRUE" : "FALSE", exp_opt->enable_isend ? "TRUE" : "FALSE"));
        is_same = false;
    }

    if (opt->max_pend_req != exp_opt->max_pend_req)
    {
        LOG((1, "Warning rearranger max_pend_req = %d, expected = %d",
             opt->max_pend_req, exp_opt->max_pend_req));
        is_same = false;
    }

    return is_same;
}

/**
 * Internal function to compare rearranger options.
 *
 * @param rearr_opts pointer to rearranger options to compare
 * @param exp_rearr_opts pointer to rearranger options with the
 * expected value
 * @return true if values in rearr_opts == values in exp_rearr_opts
 * false otherwise
 */
bool cmp_rearr_opts(const rearr_opt_t *rearr_opts, const rearr_opt_t *exp_rearr_opts)
{
    bool is_same = true;

    assert(rearr_opts && exp_rearr_opts);

    if (rearr_opts->comm_type != exp_rearr_opts->comm_type)
    {
        LOG((1, "Warning rearranger comm_type = %d, expected = %d. ", rearr_opts->comm_type,
             exp_rearr_opts->comm_type));
        is_same = false;
    }

    if (rearr_opts->fcd != exp_rearr_opts->fcd)
    {
        LOG((1, "Warning rearranger fcd = %d, expected = %d. ", rearr_opts->fcd,
             exp_rearr_opts->fcd));
        is_same = false;
    }

    is_same = is_same && cmp_rearr_comm_fc_opts(&(rearr_opts->comm_fc_opts_comp2io),
                                                &(exp_rearr_opts->comm_fc_opts_comp2io));
    is_same = is_same && cmp_rearr_comm_fc_opts(&(rearr_opts->comm_fc_opts_io2comp),
                                                &(exp_rearr_opts->comm_fc_opts_io2comp));

    return is_same;
}

/**
 * Internal function to reset rearranger opts in iosystem to valid values.
 * The old default for max pending requests was DEF_P2P_MAXREQ = 64.
 *
 * @param ios pointer to iosystem descriptor
 */
void check_and_reset_rearr_opts(iosystem_desc_t *ios)
{
    /* Disable handshake/isend and set max_pend_req to unlimited */
    const rearr_comm_fc_opt_t def_comm_nofc_opts =
        { false, false, PIO_REARR_COMM_UNLIMITED_PEND_REQ };
    /* Disable handshake /isend and set max_pend_req = 0 to turn off throttling */
    const rearr_comm_fc_opt_t def_coll_comm_fc_opts = { false, false, 0 };
    const rearr_opt_t def_coll_rearr_opts = {
        PIO_REARR_COMM_COLL,
        PIO_REARR_COMM_FC_2D_DISABLE,
        def_coll_comm_fc_opts,
        def_coll_comm_fc_opts
    };

    assert(ios);

    /* Reset to defaults, if needed (user did not set it correctly) */
    if (ios->rearr_opts.comm_type == PIO_REARR_COMM_COLL)
    {
        /* Compare and log the user and default rearr opts for coll. */
        cmp_rearr_opts(&(ios->rearr_opts), &def_coll_rearr_opts);
        /* Hard reset flow control options. */
        ios->rearr_opts = def_coll_rearr_opts;
    }
    else if (ios->rearr_opts.comm_type == PIO_REARR_COMM_P2P)
    {
        if (ios->rearr_opts.fcd == PIO_REARR_COMM_FC_2D_DISABLE)
        {
            /* Compare and log user and default opts. */
            cmp_rearr_comm_fc_opts(&(ios->rearr_opts.comm_fc_opts_comp2io),
                                   &def_comm_nofc_opts);
            cmp_rearr_comm_fc_opts(&(ios->rearr_opts.comm_fc_opts_io2comp),
                                   &def_comm_nofc_opts);
            /* Hard reset flow control opts to defaults. */
            ios->rearr_opts.comm_fc_opts_comp2io = def_comm_nofc_opts;
            ios->rearr_opts.comm_fc_opts_io2comp = def_comm_nofc_opts;
        }
        else if (ios->rearr_opts.fcd == PIO_REARR_COMM_FC_1D_COMP2IO)
        {
            /* Compare and log user and default opts. */
            cmp_rearr_comm_fc_opts(&(ios->rearr_opts.comm_fc_opts_io2comp),
                                   &def_comm_nofc_opts);
            /* Hard reset io2comp dir to defaults. */
            ios->rearr_opts.comm_fc_opts_io2comp = def_comm_nofc_opts;
        }
        else if (ios->rearr_opts.fcd == PIO_REARR_COMM_FC_1D_IO2COMP)
        {
            /* Compare and log user and default opts. */
            cmp_rearr_comm_fc_opts(&(ios->rearr_opts.comm_fc_opts_comp2io),
                                   &def_comm_nofc_opts);
            /* Hard reset comp2io dir to defaults. */
            ios->rearr_opts.comm_fc_opts_comp2io = def_comm_nofc_opts;
        }
        /* Don't reset if flow control is enabled in both directions
         * by user. */
    }
}

/**
 * Set the rearranger options associated with an iosystem
 *
 * @param comm_type Type of communication (pt2pt/coll) used
 * by the rearranger. See PIO_REARR_COMM_TYPE for more detail.
 * Possible values are :
 * PIO_REARR_COMM_P2P (Point to point communication)
 * PIO_REARR_COMM_COLL (Collective communication)
 * @param fcd Flow control direction for the rearranger.
 * See PIO_REARR_COMM_FC_DIR for more detail.
 * Possible values are :
 * PIO_REARR_COMM_FC_2D_ENABLE : Enable flow control from
 * compute processes to io processes and vice versa
 * PIO_REARR_COMM_FC_1D_COMP2IO : Enable flow control from
 * compute processes to io processes (only)
 * PIO_REARR_COMM_FC_1D_IO2COMP : Enable flow control from
 * io processes to compute processes (only)
 * PIO_REARR_COMM_FC_2D_DISABLE : Disable flow control from
 * compute processes to io processes and vice versa.
 * @param enable_hs_c2i Enable handshake while rearranging
 * data, from compute to io processes
 * @param enable_isend_c2i Enable isends while rearranging
 * data, from compute to io processes
 * @param max_pend_req_c2i Maximum pending requests during
 * data rearragment from compute processes to io processes
 * @param enable_hs_i2c Enable handshake while rearranging
 * data, from io to compute processes
 * @param enable_isend_i2c Enable isends while rearranging
 * data, from io to compute processes
 * @param max_pend_req_i2c Maximum pending requests during
 * data rearragment from io processes to compute processes
 * @param iosysidp index of the defined system descriptor
 * @return 0 on success, otherwise a PIO error code.
 */
int PIOc_set_rearr_opts(int iosysid, int comm_type, int fcd, bool enable_hs_c2i,
                        bool enable_isend_c2i, int max_pend_req_c2i,
                        bool enable_hs_i2c, bool enable_isend_i2c,
                        int max_pend_req_i2c)
{
    iosystem_desc_t *ios;
    int ret = PIO_NOERR;
    rearr_opt_t user_rearr_opts = {
        comm_type, fcd,
        {enable_hs_c2i,enable_isend_c2i, max_pend_req_c2i},
        {enable_hs_i2c, enable_isend_i2c, max_pend_req_i2c}
    };

    /* Get the IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    ios->rearr_opts = user_rearr_opts;

    /* Perform sanity checks on the user supplied values */
    check_and_reset_rearr_opts(ios);

    return ret;
}
