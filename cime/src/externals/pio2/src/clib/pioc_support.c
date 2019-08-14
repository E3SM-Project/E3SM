/** @file
 * Support functions for the PIO library.
 */
#include "config.h"
#if PIO_ENABLE_LOGGING
#include <stdarg.h>
#include <unistd.h>
#endif /* PIO_ENABLE_LOGGING */
#include <pio.h>
#include <pio_internal.h>

#include <execinfo.h>

/** This is used with text decomposition files. */
#define VERSNO 2001

/** In decomposition files, backtraces are included. This is the max
 * number of trace levels that will be used. */
#define MAX_BACKTRACE 10

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
 * Start the PIO timer.
 *
 * @param name name of the timer.
 * @return 0 for success, error code otherwise.
 * @author Ed Hartnett
 */
int
pio_start_timer(const char *name)
{
#ifdef TIMING
    GPTLstart(name);
#endif /* TIMING */
    return PIO_NOERR;
}

/**
 * Stop the PIO timer.
 *
 * @param name name of the timer.
 * @return 0 for success, error code otherwise.
 * @author Ed Hartnett
 */
int
pio_stop_timer(const char *name)
{
#ifdef TIMING
    GPTLstop(name);
#endif /* TIMING */
    return PIO_NOERR;
}

/**
 * Return a string description of an error code. If zero is passed,
 * the errmsg will be "No error".
 *
 * @param pioerr the error code returned by a PIO function call.
 * @param errmsg Pointer that will get the error message. The message
 * will be PIO_MAX_NAME chars or less.
 * @return 0 on success.
 * @author Jim Edwards
 */
int
PIOc_strerror(int pioerr, char *errmsg)
{
    PLOG((1, "PIOc_strerror pioerr = %d", pioerr));

    /* Caller must provide this. */
    pioassert(errmsg, "pointer to errmsg string must be provided", __FILE__,
              __LINE__);

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
        strncpy(errmsg, nc_strerror(pioerr), PIO_MAX_NAME);
#endif /* endif defined(_NETCDF) */
#if defined(_PNETCDF)
    else if (pioerr > PIO_FIRST_ERROR_CODE)     /* pNetCDF error? */
        strncpy(errmsg, ncmpi_strerror(pioerr), PIO_MAX_NAME);
#endif /* defined( _PNETCDF) */
    else
        /* Handle PIO errors. */
        switch(pioerr)
        {
        case PIO_EBADIOTYPE:
            strcpy(errmsg, "Bad IO type");
            break;
        case PIO_EVARDIMMISMATCH:
            strcpy(errmsg, "Variable dim mismatch in multivar call");
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
 * @author Ed Hartnett
 */
int
PIOc_set_log_level(int level)
{

#if PIO_ENABLE_LOGGING
    /* Set the log level. */
    pio_log_level = level;

#if NETCDF_C_LOGGING_ENABLED
    int ret;

    /* If netcdf logging is available turn it on starting at level = 4. */
    if (level > NC_LEVEL_DIFF)
        if ((ret = nc_set_log_level(level - NC_LEVEL_DIFF)))
            return pio_err(NULL, NULL, ret, __FILE__, __LINE__);
#endif /* NETCDF_C_LOGGING_ENABLED */
#endif /* PIO_ENABLE_LOGGING */

    return PIO_NOERR;
}

#ifdef USE_MPE

/* This array holds even numbers for MPE. */
int event_num[2][NUM_EVENTS];

/* This keeps track of whether MPE has been initialized. */
int mpe_logging_initialized = 0;

/** This will set up the MPE logging event numbers. The calling
 * program does not need to call MPE_Init_log(), that is done by the
 * mpe library in MPI_Init(). MPE must be installed, get it from
 * https://www.mcs.anl.gov/research/projects/perfvis/software/MPE/. PIO
 * and the whole I/O stack must be built with MPE.
 *
 * @param my_rank rank of processor in MPI_COMM_WORLD.
 * @author Ed Hartnett
 */
int
init_mpe(int my_rank)
{
    /* If we've already initialized MPE states, just return. */
    if (mpe_logging_initialized++)
        return 0;

    /* Get a bunch of event numbers. */
    event_num[START][INIT] = MPE_Log_get_event_number();
    event_num[END][INIT] = MPE_Log_get_event_number();
    event_num[START][DECOMP] = MPE_Log_get_event_number();
    event_num[END][DECOMP] = MPE_Log_get_event_number();
    event_num[START][CREATE] = MPE_Log_get_event_number();
    event_num[END][CREATE] = MPE_Log_get_event_number();
    event_num[START][OPEN] = MPE_Log_get_event_number();
    event_num[END][OPEN] = MPE_Log_get_event_number();
    event_num[START][DARRAY_WRITE] = MPE_Log_get_event_number();
    event_num[END][DARRAY_WRITE] = MPE_Log_get_event_number();
    event_num[START][CLOSE] = MPE_Log_get_event_number();
    event_num[END][CLOSE] = MPE_Log_get_event_number();
    event_num[START][DARRAY_READ] = MPE_Log_get_event_number();
    event_num[END][DARRAY_READ] = MPE_Log_get_event_number();

    /* On rank 0, set up the info states. */
    if (!my_rank)
    {
        /* Available colors: "white", "black", "red", "yellow", "green",
           "cyan", "blue", "magenta", "aquamarine", "forestgreen",
           "orange", "marroon", "brown", "pink", "coral", "gray" */
        MPE_Describe_info_state(event_num[START][INIT], event_num[END][INIT],
                                "PIO init", "green", "%s");
        MPE_Describe_info_state(event_num[START][DECOMP],
                                event_num[END][DECOMP], "PIO decomposition",
                                "cyan", "%s");
        MPE_Describe_info_state(event_num[START][CREATE], event_num[END][CREATE],
                                "PIO create file", "red", "%s");
        MPE_Describe_info_state(event_num[START][OPEN], event_num[END][OPEN],
                                "PIO open file", "orange", "%s");
        MPE_Describe_info_state(event_num[START][DARRAY_WRITE],
                                event_num[END][DARRAY_WRITE], "PIO darray write",
                                "pink", "%s");
        MPE_Describe_info_state(event_num[START][DARRAY_READ],
                                event_num[END][DARRAY_READ], "PIO darray read",
                                "magenta", "%s");
        MPE_Describe_info_state(event_num[START][CLOSE],
                                event_num[END][CLOSE], "PIO close",
                                "white", "%s");
    }
    return 0;
}

/**
 * Start MPE logging.
 *
 * @param state_num the MPE event state number to START (ex. INIT).
 * @author Ed Hartnett
 */
void
pio_start_mpe_log(int state)
{
    if (MPE_Log_event(event_num[START][state], 0, NULL))
        pio_err(NULL, NULL, PIO_EIO, __FILE__, __LINE__);
}

/**
 * End MPE logging.
 *
 * @param state one of the MPE states defined in pio_internal.h.
 * @param msg a text message to describe the state. Will be truncated
 * to MPE_MAX_MSG_LEN.
 * @author Ed Hartnett
 */
void
pio_stop_mpe_log(int state, const char *msg)
{
    MPE_LOG_BYTES bytebuf;
    int pos = 0;
    int msglen;
    int ret;

    /* Truncate messages longer than MPE_MAX_MSG_LEN. */
    msglen = strlen(msg) > MPE_MAX_MSG_LEN ? MPE_MAX_MSG_LEN : strlen(msg);

    /* Tell MPE to stop the state, with a message. */
    MPE_Log_pack(bytebuf, &pos, 's', msglen, msg);
    if ((ret = MPE_Log_event(event_num[END][state], 0, bytebuf)))
        pio_err(NULL, NULL, PIO_EIO, __FILE__, __LINE__);
}

#endif /* USE_MPE */

/**
 * Initialize logging.  Open log file, if not opened yet, or increment
 * ref count if already open.
 *
 * @author Jayesh Krishna, Ed Hartnett
 */
int
pio_init_logging(void)
{
    int ret = PIO_NOERR;

#ifdef USE_MPE
    {
        int mpe_rank;
        int mpierr;

        if ((mpierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpe_rank)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);

        if ((ret = init_mpe(mpe_rank)))
            return pio_err(NULL, NULL, ret, __FILE__, __LINE__);
    }
#endif /* USE_MPE */

#if PIO_ENABLE_LOGGING
    if (!LOG_FILE)
    {
        char log_filename[PIO_MAX_NAME];
        int mpierr;

        /* Create a filename with the rank in it. */
        if ((mpierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
        sprintf(log_filename, "pio_log_%d.log", my_rank);

        /* Open a file for this rank to log messages. */
        if (!(LOG_FILE = fopen(log_filename, "w")))
            return pio_err(NULL, NULL, PIO_EIO, __FILE__, __LINE__);

        pio_log_ref_cnt = 1;
    }
    else
    {
        pio_log_ref_cnt++;
    }
#endif /* PIO_ENABLE_LOGGING */

    return ret;
}

/**
 * Finalize logging - close log files, if open.
 */
void
pio_finalize_logging(void)
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
            PLOG((2, "pio_finalize_logging, postpone close, ref_cnt = %d",
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
 * In code this functions should be wrapped in the PLOG(()) macro.
 *
 * @param severity the severity of the message, 0 for error messages,
 * then increasing levels of verbosity.
 * @param fmt the format string.
 * @param ... the arguments used in format string.
 * @author Ed Hartnett
 */
void
pio_log(int severity, const char *fmt, ...)
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

    /* /\* Show the severity. *\/ */
    /* snprintf(rank_str, MAX_RANK_STR, ":%d ", severity); */
    /* strncpy(ptr, rank_str, (rem_len > 0) ? rem_len : 0); */
    /* ptr += strlen(rank_str); */
    /* rem_len -= strlen(rank_str); */

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
 * @author Jim Edwards
 */
void
print_trace(FILE *fp)
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
 * Abort program and call MPI_Abort().
 *
 * @param msg an error message
 * @param fname name of code file where error occured
 * @param line the line of code where the error occurred.
 * @author Jim Edwards
 */
void
piodie(const char *msg, const char *fname, int line)
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
 * @author Jim Edwards
 */
void
pioassert(_Bool expression, const char *msg, const char *fname, int line)
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
 * @param ios pointer to the iosystem_info_t. May be NULL.
 * @param file pointer to the file_desc_t info. Ignored if NULL.
 * @param mpierr the MPI return code to handle
 * @param filename the name of the code file where error occured.
 * @param line the line of code where error occured.
 * @return PIO_NOERR for no error, otherwise PIO_EIO.
 * @author Ed Hartnett
 */
int
check_mpi(iosystem_desc_t *ios, file_desc_t *file, int mpierr,
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
 * @author Ed Hartnett
 */
int
check_netcdf(file_desc_t *file, int status, const char *fname, int line)
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
 * @author Ed Hartnett
 */
int
check_netcdf2(iosystem_desc_t *ios, file_desc_t *file, int status,
              const char *fname, int line)
{
    int eh = default_error_handler; /* Error handler that will be used. */
    int rbuf;
    /* User must provide this. */
    pioassert(fname, "code file name must be provided", __FILE__, __LINE__);

    if (file && file->iosystem->ioproc &&
        (file->iotype == PIO_IOTYPE_PNETCDF || file->iotype == PIO_IOTYPE_NETCDF4P))
    {
        if (file->iosystem->io_rank == 0)
            MPI_Reduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_MIN, 0, file->iosystem->io_comm);
        else
            MPI_Reduce(&status, &rbuf, 1, MPI_INT, MPI_MIN, 0, file->iosystem->io_comm);
    }

    PLOG((1, "check_netcdf2 status = %d fname = %s line = %d", status, fname, line));

    /* Pick an error handler. */
    if (ios)
        eh = ios->error_handler;
    if (file)
        eh = file->iosystem->error_handler;
    pioassert(eh == PIO_INTERNAL_ERROR || eh == PIO_BCAST_ERROR || eh == PIO_RETURN_ERROR,
              "invalid error handler", __FILE__, __LINE__);
    PLOG((2, "check_netcdf2 chose error handler = %d", eh));

    /* Decide what to do based on the error handler. */
    if (eh == PIO_INTERNAL_ERROR && status != PIO_NOERR)
    {
        char errmsg[PIO_MAX_NAME + 1];  /* Error message. */
        PIOc_strerror(status, errmsg);
        piodie(errmsg, fname, line);        /* Die! */
    }
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
 * @author Jim Edwards
 */
int
pio_err(iosystem_desc_t *ios, file_desc_t *file, int err_num, const char *fname,
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
    PLOG((0, "%s err_num = %d fname = %s line = %d", err_msg, err_num, fname ? fname : '\0', line));

    /* What error handler should we use? */
    if (file)
        err_handler = file->iosystem->error_handler;
    else if (ios)
        err_handler = ios->error_handler;

    PLOG((2, "pio_err chose error handler = %d", err_handler));

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
 * Allocate a region struct, and initialize it.
 *
 * @param ios pointer to the IO system info, used for error
 * handling. Ignored if NULL.
 * @param ndims the number of dimensions for the data in this region.
 * @param regionp a pointer that gets a pointer to the newly allocated
 * io_region struct.
 * @returns 0 for success, error code otherwise.
 * @author Jim Edwards
 */
int
alloc_region2(iosystem_desc_t *ios, int ndims, io_region **regionp)
{
    io_region *region;

    /* Check inputs. */
    pioassert(ndims >= 0 && regionp, "invalid input", __FILE__, __LINE__);
    PLOG((1, "alloc_region2 ndims = %d sizeof(io_region) = %d", ndims,
          sizeof(io_region)));

    /* Allocate memory for the io_region struct. */
    if (!(region = calloc(1, sizeof(io_region))))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Allocate memory for the array of start indicies. */
    if (!(region->start = calloc(ndims, sizeof(PIO_Offset))))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Allocate memory for the array of counts. */
    if (!(region->count = calloc(ndims, sizeof(PIO_Offset))))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Return pointer to new region to caller. */
    *regionp = region;

    return PIO_NOERR;
}

/**
 * Given a PIO type, find the MPI type and the type size.
 *
 * @param pio_type a PIO type, PIO_INT, PIO_FLOAT, etc.
 * @param mpi_type a pointer to MPI_Datatype that will get the MPI
 * type that coresponds to the PIO type. Ignored if NULL.
 * @param type_size a pointer to int that will get the size of the
 * type, in bytes. (For example, 4 for PIO_INT). Ignored if NULL.
 * @returns 0 for success, error code otherwise.
 * @author Jim Edwards
 */
int
find_mpi_type(int pio_type, MPI_Datatype *mpi_type, int *type_size)
{
    MPI_Datatype my_mpi_type;
    int my_type_size;

    /* Decide on the base type. */
    switch(pio_type)
    {
    case PIO_BYTE:
        my_mpi_type = MPI_BYTE;
        my_type_size = NETCDF_CHAR_SIZE;
        break;
    case PIO_CHAR:
        my_mpi_type = MPI_CHAR;
        my_type_size = NETCDF_CHAR_SIZE;
        break;
    case PIO_SHORT:
        my_mpi_type = MPI_SHORT;
        my_type_size = NETCDF_SHORT_SIZE;
        break;
    case PIO_INT:
        my_mpi_type = MPI_INT;
        my_type_size = NETCDF_INT_FLOAT_SIZE;
        break;
    case PIO_FLOAT:
        my_mpi_type = MPI_FLOAT;
        my_type_size = NETCDF_INT_FLOAT_SIZE;
        break;
    case PIO_DOUBLE:
        my_mpi_type = MPI_DOUBLE;
        my_type_size = NETCDF_DOUBLE_INT64_SIZE;
        break;
#ifdef _NETCDF4
    case PIO_UBYTE:
        my_mpi_type = MPI_UNSIGNED_CHAR;
        my_type_size = NETCDF_CHAR_SIZE;
        break;
    case PIO_USHORT:
        my_mpi_type = MPI_UNSIGNED_SHORT;
        my_type_size = NETCDF_SHORT_SIZE;
        break;
    case PIO_UINT:
        my_mpi_type = MPI_UNSIGNED;
        my_type_size = NETCDF_INT_FLOAT_SIZE;
        break;
    case PIO_INT64:
        my_mpi_type = MPI_LONG_LONG;
        my_type_size = NETCDF_DOUBLE_INT64_SIZE;
        break;
    case PIO_UINT64:
        my_mpi_type = MPI_UNSIGNED_LONG_LONG;
        my_type_size = NETCDF_DOUBLE_INT64_SIZE;
        break;
    case PIO_STRING:
        my_mpi_type = MPI_CHAR;
        my_type_size = NETCDF_CHAR_SIZE;
        break;
#endif /* _NETCDF4 */
    default:
        return PIO_EBADTYPE;
    }

    /* If caller wants MPI type, set it. */
    if (mpi_type)
        *mpi_type = my_mpi_type;

    /* If caller wants type size, set it. */
    if (type_size)
        *type_size = my_type_size;

    return PIO_NOERR;
}

/**
 * Allocate space for an IO description struct, and initialize it.
 *
 * @param ios pointer to the IO system info, used for error
 * handling.
 * @param piotype the PIO data type (ex. PIO_FLOAT, PIO_INT, etc.).
 * @param ndims the number of dimensions.
 * @param iodesc pointer that gets the newly allocated io_desc_t.
 * @returns 0 for success, error code otherwise.
 * @author Jim Edwards
 */
int
malloc_iodesc(iosystem_desc_t *ios, int piotype, int ndims,
              io_desc_t **iodesc)
{
    MPI_Datatype mpi_type;
    PIO_Offset type_size;
    int mpierr;
    int ret;

    /* Check input. */
    pioassert(ios && piotype > 0 && ndims >= 0 && iodesc,
              "invalid input", __FILE__, __LINE__);

    PLOG((1, "malloc_iodesc piotype = %d ndims = %d", piotype, ndims));

    /* Get the MPI type corresponding with the PIO type. */
    if ((ret = find_mpi_type(piotype, &mpi_type, NULL)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* What is the size of the pio type? */
    if ((ret = pioc_pnetcdf_inq_type(0, piotype, NULL, &type_size)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Allocate space for the io_desc_t struct. */
    if (!(*iodesc = calloc(1, sizeof(io_desc_t))))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Remember the pio type and its size. */
    (*iodesc)->piotype = piotype;
    (*iodesc)->piotype_size = type_size;

    /* Remember the MPI type. */
    (*iodesc)->mpitype = mpi_type;

    /* Get the size of the type. */
    if (mpi_type == MPI_DATATYPE_NULL)
        (*iodesc)->mpitype_size = 0;
    else
        if ((mpierr = MPI_Type_size((*iodesc)->mpitype, &(*iodesc)->mpitype_size)))
            return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Initialize some values in the struct. */
    (*iodesc)->maxregions = 1;
    (*iodesc)->ioid = -1;
    (*iodesc)->ndims = ndims;

    /* Allocate space for, and initialize, the first region. */
    if ((ret = alloc_region2(ios, ndims, &((*iodesc)->firstregion))))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Set the swap memory settings to defaults for this IO system. */
    (*iodesc)->rearr_opts = ios->rearr_opts;

    return PIO_NOERR;
}

/**
 * Free a region list.
 *
 * top a pointer to the start of the list to free.
 * @author Jim Edwards
 */
void
free_region_list(io_region *top)
{
    io_region *ptr, *tptr;

    ptr = top;
    while (ptr)
    {
        if (ptr->start)
            free(ptr->start);
        if (ptr->count)
            free(ptr->count);
        tptr = ptr;
        ptr = ptr->next;
        free(tptr);
    }
}

/**
 * Free a decomposition map.
 *
 * @param iosysid the IO system ID.
 * @param ioid the ID of the decomposition map to free.
 * @returns 0 for success, error code otherwise.
 * @ingroup PIO_freedecomp_c
 * @author Jim Edwards
 */
int
PIOc_freedecomp(int iosysid, int ioid)
{
    iosystem_desc_t *ios;
    io_desc_t *iodesc;
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function calls. */

    PLOG((1, "PIOc_freedecomp iosysid = %d ioid = %d", iosysid, ioid));

    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    if (!(iodesc = pio_get_iodesc_from_id(ioid)))
        return pio_err(ios, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async)
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
            PLOG((2, "PIOc_freedecomp iosysid = %d ioid = %d", iosysid, ioid));
        }

        /* Handle MPI errors. */
        PLOG((3, "handline error mpierr %d ios->comproot %d", mpierr, ios->comproot));
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(NULL, NULL, mpierr2, __FILE__, __LINE__);
        PLOG((3, "handline error mpierr2 %d", mpierr2));
        if (mpierr)
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
    }

    PLOG((3, "freeing map, dimlen"));
    /* Free the map. */
    free(iodesc->map);

    /* Free the dimlens. */
    free(iodesc->dimlen);

    if (iodesc->remap)
        free(iodesc->remap);

    PLOG((3, "freeing rfrom, rtype"));
    if (iodesc->rfrom)
        free(iodesc->rfrom);

    if (iodesc->rtype)
    {
        for (int i = 0; i < iodesc->nrecvs; i++)
            if (iodesc->rtype[i] != PIO_DATATYPE_NULL)
                if ((mpierr = MPI_Type_free(&iodesc->rtype[i])))
                    return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);

        free(iodesc->rtype);
    }

    PLOG((3, "freeing stype, scount"));
    if (iodesc->stype)
    {
        for (int i = 0; i < iodesc->num_stypes; i++)
            if (iodesc->stype[i] != PIO_DATATYPE_NULL)
                if ((mpierr = MPI_Type_free(iodesc->stype + i)))
                    return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);

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

    PLOG((3, "freeing regions"));
    if (iodesc->firstregion)
        free_region_list(iodesc->firstregion);

    if (iodesc->fillregion)
        free_region_list(iodesc->fillregion);

    if (iodesc->rearranger == PIO_REARR_SUBSET)
        if ((mpierr = MPI_Comm_free(&iodesc->subset_comm)))
            return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);

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
 * @author Jim Edwards
 */
int
PIOc_readmap(const char *file, int *ndims, int **gdims, PIO_Offset *fmaplen,
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
        return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &myrank)))
        return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);

    if (myrank == 0)
    {
        FILE *fp = fopen(file, "r");
        if (!fp)
            pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

        if (fscanf(fp,"version %d npes %d ndims %d\n", &rversno, &rnpes, ndims) != 3)
            pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);
        if (rversno != VERSNO)
            return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

        if (rnpes < 1 || rnpes > npes)
            return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

        if ((mpierr = MPI_Bcast(&rnpes, 1, MPI_INT, 0, comm)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(ndims, 1, MPI_INT, 0, comm)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
        if (!(tdims = calloc(*ndims, sizeof(int))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        for (int i = 0; i < *ndims; i++)
            if (fscanf(fp,"%d ", tdims + i) != 1)
                pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

        if ((mpierr = MPI_Bcast(tdims, *ndims, MPI_INT, 0, comm)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);

        for (int i = 0; i < rnpes; i++)
        {
            if (fscanf(fp, "%d %lld", &j, &maplen) != 2)
                pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);
            if (j != i)  // Not sure how this could be possible
                return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);
            if (!(tmap = malloc(maplen * sizeof(PIO_Offset))))
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
            for (j = 0; j < maplen; j++)
                if (fscanf(fp, "%lld ", tmap + j) != 1)
                    pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

            if (i > 0)
            {
                if ((mpierr = MPI_Send(&maplen, 1, PIO_OFFSET, i, i + npes, comm)))
                    return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(tmap, maplen, PIO_OFFSET, i, i, comm)))
                    return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
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
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(ndims, 1, MPI_INT, 0, comm)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
        if (!(tdims = calloc(*ndims, sizeof(int))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(tdims, *ndims, MPI_INT, 0, comm)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);

        if (myrank < rnpes)
        {
            if ((mpierr = MPI_Recv(&maplen, 1, PIO_OFFSET, 0, myrank + npes, comm, &status)))
                return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
            if (!(tmap = malloc(maplen * sizeof(PIO_Offset))))
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
            if ((mpierr = MPI_Recv(tmap, maplen, PIO_OFFSET, 0, myrank, comm, &status)))
                return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
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
 * @author Jim Edwards
 */
int
PIOc_readmap_from_f90(const char *file, int *ndims, int **gdims, PIO_Offset *maplen,
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
 * @param title optial title attribute for the file. Must be less than
 * PIO_MAX_NAME + 1 if provided. Ignored if NULL.
 * @param history optial history attribute for the file. Must be less
 * than PIO_MAX_NAME + 1 if provided. Ignored if NULL.
 * @param fortran_order set to non-zero if fortran array ordering is
 * used, or to zero if C array ordering is used.
 * @returns 0 for success, error code otherwise.
 * @author Ed Hartnett
 */
int
PIOc_write_nc_decomp(int iosysid, const char *filename, int cmode, int ioid,
                     char *title, char *history, int fortran_order)
{
    iosystem_desc_t *ios; /* IO system info. */
    io_desc_t *iodesc;    /* Decomposition info. */
    int max_maplen;       /* The maximum maplen used for any task. */
    int *full_map;        /* 2D array holds all map info for all tasks. */
    int *my_map;          /* 1D array holds all map info for this task. */
    int mpierr;
    int ret;

    /* Get the IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Check inputs. */
    if (!filename)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);
    if (title)
        if (strlen(title) > PIO_MAX_NAME)
            return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);
    if (history)
        if (strlen(history) > PIO_MAX_NAME)
            return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    PLOG((1, "PIOc_write_nc_decomp filename = %s iosysid = %d ioid = %d "
          "ios->num_comptasks = %d", filename, iosysid, ioid, ios->num_comptasks));

    /* Get the IO desc, which describes the decomposition. */
    if (!(iodesc = pio_get_iodesc_from_id(ioid)))
        return pio_err(ios, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Allocate memory for array which will contain the length of the
     * map on each task, for all computation tasks. */
    int task_maplen[ios->num_comptasks];
    PLOG((3, "ios->num_comptasks = %d", ios->num_comptasks));

    /* Gather maplens from all computation tasks and fill the
     * task_maplen array on all tasks. */
    if ((mpierr = MPI_Allgather(&iodesc->maplen, 1, MPI_INT, task_maplen, 1, MPI_INT,
                                ios->comp_comm)))
        return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);

    /* Find the max maplen. */
    if ((mpierr = MPI_Allreduce(&iodesc->maplen, &max_maplen, 1, MPI_INT, MPI_MAX,
                                ios->comp_comm)))
        return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);
    PLOG((3, "max_maplen = %d", max_maplen));

    if (!(full_map = malloc(sizeof(int) * ios->num_comptasks * max_maplen)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    if (!(my_map = malloc(sizeof(int) * max_maplen)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Fill local array with my map. Use the fill value for unused */
    /* elements at the end if max_maplen is longer than maplen. Also
     * subtract 1 because the iodesc->map is 1-based. */
    for (int e = 0; e < max_maplen; e++)
    {
        my_map[e] = e < iodesc->maplen ? iodesc->map[e] - 1 : NC_FILL_INT;
        PLOG((3, "my_map[%d] = %d", e, my_map[e]));
    }

    /* Gather my_map from all computation tasks and fill the full_map array. */
    if ((mpierr = MPI_Allgather(my_map, max_maplen, MPI_INT, full_map, max_maplen,
                                MPI_INT, ios->comp_comm)))
        return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);

    free(my_map);

    for (int p = 0; p < ios->num_comptasks; p++)
        for (int e = 0; e < max_maplen; e++)
            PLOG((3, "full_map[%d][%d] = %d", p, e, full_map[p * max_maplen + e]));

    /* Write the netCDF decomp file. */
    if ((ret = pioc_write_nc_decomp_int(ios, filename, cmode, iodesc->ndims, iodesc->dimlen,
                                        ios->num_comptasks, task_maplen, full_map, title,
                                        history, fortran_order)))
        return ret;

    free(full_map);
    return PIO_NOERR;
}

/**
 * Read the decomposition map from a netCDF decomp file produced by
 * PIOc_write_nc_decomp().
 *
 * @param iosysid the IO system ID.
 * @param filename the name of the decomp file.
 * @param ioidp pointer that will get the newly-assigned ID of the IO
 * description. The ioid is needed to later free the decomposition.
 * @param comm an MPI communicator.
 * @param pio_type the PIO type to be used as the type for the data.
 * @param title pointer that will get optial title attribute for the
 * file. Will be less than PIO_MAX_NAME + 1 if provided. Ignored if
 * NULL.
 * @param history pointer that will get optial history attribute for
 * the file. Will be less than PIO_MAX_NAME + 1 if provided. Ignored if
 * NULL.
 * @param fortran_order pointer that gets set to 1 if fortran array
 * ordering is used, or to zero if C array ordering is used.
 * @returns 0 for success, error code otherwise.
 * @author Ed Hartnett
 */
int
PIOc_read_nc_decomp(int iosysid, const char *filename, int *ioidp, MPI_Comm comm,
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

    PLOG((1, "PIOc_read_nc_decomp filename = %s iosysid = %d pio_type = %d",
          filename, iosysid, pio_type));

    /* Get the communicator size and task rank. */
    if ((mpierr = MPI_Comm_size(comm, &size)))
        return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &my_rank)))
        return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);
    PLOG((2, "size = %d my_rank = %d", size, my_rank));

    /* Read the file. This allocates three arrays that we have to
     * free. */
    if ((ret = pioc_read_nc_decomp_int(iosysid, filename, &ndims, &global_dimlen, &num_tasks_decomp,
                                       &task_maplen, &max_maplen, &full_map, title, history,
                                       source_in, version_in, fortran_order)))
        return ret;
    PLOG((2, "ndims = %d num_tasks_decomp = %d max_maplen = %d", ndims, num_tasks_decomp,
          max_maplen));

    /* If the size does not match the number of tasks in the decomp,
     * that's an error. */
    if (size != num_tasks_decomp)
        ret = PIO_EINVAL;

    /* Now initialize the iodesc on each task for this decomposition. */
    if (!ret)
    {
        PIO_Offset *compmap;

        if (!(compmap = malloc(sizeof(PIO_Offset) * task_maplen[my_rank])))
            return PIO_ENOMEM;

        /* Copy array into PIO_Offset array. Make it 1 based. */
        for (int e = 0; e < task_maplen[my_rank]; e++)
            compmap[e] = full_map[my_rank * max_maplen + e] + 1;

        /* Initialize the decomposition. */
        ret = PIOc_InitDecomp(iosysid, pio_type, ndims, global_dimlen, task_maplen[my_rank],
                              compmap, ioidp, NULL, NULL, NULL);

        free(compmap);
    }

    /* Free resources. */
    free(global_dimlen);
    free(task_maplen);
    free(full_map);

    return ret;
}

/**
 * Write the decomp information in netCDF. This is an internal
 * function.
 *
 * @param ios pointer to io system info.
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
 * @author Ed Hartnett
 */
int
pioc_write_nc_decomp_int(iosystem_desc_t *ios, const char *filename, int cmode, int ndims,
                         int *global_dimlen, int num_tasks, int *task_maplen, int *map,
                         const char *title, const char *history, int fortran_order)
{
    int max_maplen = 0;
    int ncid;
    int ret;

    /* Check inputs. */
    pioassert(ios && filename && global_dimlen && task_maplen &&
              (!title || strlen(title) <= PIO_MAX_NAME) &&
              (!history || strlen(history) <= PIO_MAX_NAME), "invalid input",
              __FILE__, __LINE__);

    PLOG((2, "pioc_write_nc_decomp_int filename = %s ndims = %d num_tasks = %d", filename,
          ndims, num_tasks));

    /* Find the maximum maplen. */
    for (int t = 0; t < num_tasks; t++)
        if (task_maplen[t] > max_maplen)
            max_maplen = task_maplen[t];
    PLOG((3, "max_maplen = %d", max_maplen));

    /* Create the netCDF decomp file. */
    if ((ret = PIOc_create(ios->iosysid, filename, cmode | NC_WRITE, &ncid)))
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
    if (max_bt_size > PIO_MAX_NAME)
        max_bt_size = PIO_MAX_NAME;

    /* Copy the backtrace into one long string. */
    char full_bt[max_bt_size * bt_size + bt_size + 1];
    full_bt[0] = '\0';
    for (int b = 0; b < bt_size; b++)
    {
        strncat(full_bt, bt_strings[b], max_bt_size);
        strcat(full_bt, "\n");
    }
    free(bt_strings);

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

/**
 * Read the decomp information from a netCDF decomp file. This is an
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
 * @author Ed Hartnett
 */
int
pioc_read_nc_decomp_int(int iosysid, const char *filename, int *ndims, int **global_dimlen,
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

    PLOG((1, "pioc_read_nc_decomp_int iosysid = %d filename = %s", iosysid, filename));

    /* Open the netCDF decomp file. */
    if ((ret = PIOc_open(iosysid, filename, NC_WRITE, &ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);

    /* Read version attribute. */
    char version_in[PIO_MAX_NAME + 1];
    if ((ret = PIOc_get_att_text(ncid, NC_GLOBAL, DECOMP_VERSION_ATT_NAME, version_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    PLOG((3, "version_in = %s", version_in));
    if (version)
        strncpy(version, version_in, PIO_MAX_NAME + 1);

    /* Read order attribute. */
    char order_in[PIO_MAX_NAME + 1];
    if ((ret = PIOc_get_att_text(ncid, NC_GLOBAL, DECOMP_ORDER_ATT_NAME, order_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    PLOG((3, "order_in = %s", order_in));
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
    PLOG((3, "max_maplen_in = %d", max_maplen_in));
    if (max_maplen)
        *max_maplen = max_maplen_in;

    /* Read title attribute, if it is in the file. */
    char title_in[PIO_MAX_NAME + 1];
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
    char history_in[PIO_MAX_NAME + 1];
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
    char source_in[PIO_MAX_NAME + 1];
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
    int *map_in;

    if (!(map_in = malloc(sizeof(int) * num_tasks_in * max_maplen_in)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    if ((ret = PIOc_inq_varid(ncid, DECOMP_MAP_VAR_NAME, &map_varid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if ((ret = PIOc_get_var_int(ncid, map_varid, map_in)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    if (map)
    {
        if (!(*map = malloc(num_tasks_in * max_maplen_in * sizeof(int))))
            return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);
        for (int t = 0; t < num_tasks_in; t++)
            for (int l = 0; l < max_maplen_in; l++)
                (*map)[t * max_maplen_in + l] = map_in[t * max_maplen_in + l];
    }
    free(map_in);

    /* Close the netCDF decomp file. */
    PLOG((2, "pioc_read_nc_decomp_int about to close file ncid = %d", ncid));
    if ((ret = PIOc_closefile(ncid)))
        return pio_err(ios, NULL, ret, __FILE__, __LINE__);
    PLOG((2, "pioc_read_nc_decomp_int closed file"));

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
 * @author Jim Edwards
 */
int
PIOc_write_decomp(const char *file, int iosysid, int ioid, MPI_Comm comm)
{
    iosystem_desc_t *ios;
    io_desc_t *iodesc;

    PLOG((1, "PIOc_write_decomp file = %s iosysid = %d ioid = %d", file, iosysid, ioid));

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
 * @author Jim Edwards
 */
int
PIOc_writemap(const char *file, int ndims, const int *gdims, PIO_Offset maplen,
              PIO_Offset *map, MPI_Comm comm)
{
    int npes, myrank;
    PIO_Offset *nmaplen = NULL;
    MPI_Status status;
    int i;
    PIO_Offset *nmap;
    int mpierr; /* Return code for MPI calls. */

    PLOG((1, "PIOc_writemap file = %s ndims = %d maplen = %d", file, ndims, maplen));

    if ((mpierr = MPI_Comm_size(comm, &npes)))
        return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &myrank)))
        return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
    PLOG((2, "npes = %d myrank = %d", npes, myrank));

    /* Allocate memory for the nmaplen. */
    if (myrank == 0)
        if (!(nmaplen = malloc(npes * sizeof(PIO_Offset))))
            return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    if ((mpierr = MPI_Gather(&maplen, 1, PIO_OFFSET, nmaplen, 1, PIO_OFFSET, 0, comm)))
        return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);

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
            PLOG((2, "creating nmap for i = %d", i));
            nmap = (PIO_Offset *)malloc(nmaplen[i] * sizeof(PIO_Offset));

            if ((mpierr = MPI_Send(&i, 1, MPI_INT, i, npes + i, comm)))
                return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
            if ((mpierr = MPI_Recv(nmap, nmaplen[i], PIO_OFFSET, i, i, comm, &status)))
                return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
            PLOG((2,"MPI_Recv map complete"));

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
        PLOG((2,"decomp file closed."));
    }
    else
    {
        PLOG((2,"ready to MPI_Recv..."));
        if ((mpierr = MPI_Recv(&i, 1, MPI_INT, 0, npes+myrank, comm, &status)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
        PLOG((2,"MPI_Recv got %d", i));
        if ((mpierr = MPI_Send(map, maplen, PIO_OFFSET, 0, myrank, comm)))
            return check_mpi(NULL, NULL, mpierr, __FILE__, __LINE__);
        PLOG((2,"MPI_Send map complete"));
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
 * @param f90_comm an MPI communicator.
 * @returns 0 for success, error code otherwise.
 * @author Jim Edwards
 */
int
PIOc_writemap_from_f90(const char *file, int ndims, const int *gdims,
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
 * PIOc_Init_Intracomm() or PIOc_InitAsync().
 * @param ncidp A pointer that gets the ncid of the newly created
 * file. For NetCDF integration, this contains the ncid assigned by
 * the netCDF layer, which is used instead of a PIO-generated ncid.
 * @param iotype A pointer to a pio output format. Must be one of
 * PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF, PIO_IOTYPE_NETCDF4C, or
 * PIO_IOTYPE_NETCDF4P.
 * @param filename The filename to create.
 * @param mode The netcdf mode for the create operation.
 * @param use_ext_ncid non-zero to use an externally assigned ncid
 * (used in the netcdf integration layer).
 *
 * @returns 0 for success, error code otherwise.
 * @ingroup PIO_createfile_c
 * @author Ed Hartnett
 */
int
PIOc_createfile_int(int iosysid, int *ncidp, int *iotype, const char *filename,
                    int mode, int use_ext_ncid)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
    int ierr;              /* Return code from function calls. */

#ifdef USE_MPE
    pio_start_mpe_log(CREATE);
#endif /* USE_MPE */

    /* Get the IO system info from the iosysid. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* User must provide valid input for these parameters. */
    if (!ncidp || !iotype || !filename || strlen(filename) > PIO_MAX_NAME)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* A valid iotype must be specified. */
    if (!iotype_is_valid(*iotype))
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    PLOG((1, "PIOc_createfile_int iosysid = %d iotype = %d filename = %s mode = %d",
          iosysid, *iotype, filename, mode));

    /* Allocate space for the file info. */
    if (!(file = calloc(sizeof(file_desc_t), 1)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Fill in some file values. */
    file->fh = -1;
    file->iosystem = ios;
    file->iotype = *iotype;
    file->buffer = NULL;
    file->writable = 1;

    /* Set to true if this task should participate in IO (only true for
     * one task with netcdf serial files. */
    if (file->iotype == PIO_IOTYPE_NETCDF4P || file->iotype == PIO_IOTYPE_PNETCDF ||
        ios->io_rank == 0)
        file->do_io = 1;

    PLOG((2, "file->do_io = %d ios->async = %d", file->do_io, ios->async));

    /* If async is in use, and this is not an IO task, bcast the
     * parameters. */
    if (ios->async)
    {
        if (!ios->ioproc)
        {
            int msg = PIO_MSG_CREATE_FILE;
            size_t len = strlen(filename);

            /* Send the message to the message handler. */
            PLOG((3, "msg %d ios->union_comm %d MPI_COMM_NULL %d", msg, ios->union_comm, MPI_COMM_NULL));
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
                mpierr = MPI_Bcast(&mode, 1, MPI_INT, ios->compmaster, ios->intercomm);
            PLOG((2, "len = %d filename = %s iotype = %d mode = %d", len, filename,
                  file->iotype, mode));
        }

        /* Handle MPI errors. */
        PLOG((2, "handling mpi errors mpierr = %d", mpierr));
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
    }

    /* If this task is in the IO component, do the IO. */
    if (ios->ioproc)
    {
        switch (file->iotype)
        {
#ifdef _NETCDF4
        case PIO_IOTYPE_NETCDF4P:
            mode = mode |  NC_MPIIO | NC_NETCDF4;
            PLOG((2, "Calling nc_create_par io_comm = %d mode = %d fh = %d",
                  ios->io_comm, mode, file->fh));
            ierr = nc_create_par(filename, mode, ios->io_comm, ios->info, &file->fh);
            PLOG((2, "nc_create_par returned %d file->fh = %d", ierr, file->fh));
            break;
        case PIO_IOTYPE_NETCDF4C:
            mode = mode | NC_NETCDF4;
#endif
        case PIO_IOTYPE_NETCDF:
            if (!ios->io_rank)
            {
                PLOG((2, "Calling nc_create mode = %d", mode));
                ierr = nc_create(filename, mode, &file->fh);
            }
            break;
#ifdef _PNETCDF
        case PIO_IOTYPE_PNETCDF:
            PLOG((2, "Calling ncmpi_create mode = %d", mode));
            ierr = ncmpi_create(ios->io_comm, filename, mode, ios->info, &file->fh);
            if (!ierr)
                ierr = ncmpi_buffer_attach(file->fh, pio_buffer_size_limit);
            break;
#endif
        }
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);

    /* If there was an error, free the memory we allocated and handle error. */
    if (ierr)
    {
        free(file);
        return check_netcdf2(ios, NULL, ierr, __FILE__, __LINE__);
    }

    /* Broadcast writablility to all tasks. */
    if ((mpierr = MPI_Bcast(&file->writable, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);

    /* Broadcast next ncid to all tasks from io root, necessary
     * because files may be opened on mutilple iosystems, causing the
     * underlying library to reuse ncids. Hilarious confusion
     * ensues. */
    if (ios->async)
    {
        PLOG((3, "createfile bcasting pio_next_ncid %d", pio_next_ncid));
        if ((mpierr = MPI_Bcast(&pio_next_ncid, 1, MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
        PLOG((3, "createfile bcast pio_next_ncid %d", pio_next_ncid));
    }

    /* With the netCDF integration layer, the ncid is assigned for PIO
     * by the netCDF dispatch layer code. So it is passed in. In
     * normal PIO operation, the ncid is generated here. */
    if (use_ext_ncid)
    {
        /* Use the ncid passed in from the netCDF dispatch code. */
        file->pio_ncid = *ncidp;

        /* To prevent PIO from reusing the same ncid, if someone
         * starting mingling netcdf integration PIO and regular PIO
         * code. */
        pio_next_ncid = file->pio_ncid + 1;
    }
    else
    {
        /* Assign the PIO ncid. */
        file->pio_ncid = pio_next_ncid++;
        PLOG((2, "file->fh = %d file->pio_ncid = %d", file->fh, file->pio_ncid));

        /* Return the ncid to the caller. */
        *ncidp = file->pio_ncid;
    }

    /* Add the struct with this files info to the global list of
     * open files. */
    pio_add_to_file_list(file);

#ifdef USE_MPE
    pio_stop_mpe_log(CREATE, __func__);
#endif /* USE_MPE */
    PLOG((2, "Created file %s file->fh = %d file->pio_ncid = %d", filename,
          file->fh, file->pio_ncid));

    return ierr;
}

/**
 * Check that a file meets PIO requirements for use of unlimited
 * dimensions. This function is only called on netCDF-4 files. If the
 * file is found to violate PIO requirements it is closed.
 *
 * @param ncid the file->fh for this file (the real netCDF ncid, not
 * the pio_ncid).
 * @returns 0 if file is OK, error code otherwise.
 * @author Ed Hartnett
 */
int
check_unlim_use(int ncid)
{
#ifdef _NETCDF4
    int nunlimdims; /* Number of unlimited dims in file. */
    int nvars;       /* Number of vars in file. */
    int ierr;        /* Return code. */

    /* Are there 2 or more unlimited dims in this file? */
    if ((ierr = nc_inq_unlimdims(ncid, &nunlimdims, NULL)))
        return ierr;
    if (nunlimdims < 2)
        return PIO_NOERR;

    /* How many vars in file? */
    if ((ierr = nc_inq_nvars(ncid, &nvars)))
        return ierr;

    /* Check each var. */
    for (int v = 0; v < nvars && !ierr; v++)
    {
        int nvardims;
        if ((ierr = nc_inq_varndims(ncid, v, &nvardims)))
            return ierr;
        int vardimid[nvardims];
        if ((ierr = nc_inq_vardimid(ncid, v, vardimid)))
            return ierr;

        /* Check all var dimensions, except the first. If we find
         * unlimited, that's a problem. */
        for (int vd = 1; vd < nvardims; vd++)
        {
            size_t dimlen;
            if ((ierr = nc_inq_dimlen(ncid, vardimid[vd], &dimlen)))
                return ierr;
            if (dimlen == NC_UNLIMITED)
            {
                nc_close(ncid);
                return PIO_EINVAL;
            }
        }
    }
#endif /* _NETCDF4 */

    return PIO_NOERR;
}

/**
 * Internal function used when opening an existing file. This function
 * is called by PIOc_openfile_retry(). It learns some things about the
 * metadata in that file. The results end up in the file_desc_t.
 *
 * @param file pointer to the file_desc_t for this file.
 * @param ncid the ncid assigned to the file when opened.
 * @param iotype the iotype used to open the file.
 * @param nvars a pointer that gets the number of vars in the file.

 * @param rec_var gets an array (length nvars) of rec_var values for
 * each var in the file. This array must be freed by caller.
 * @param pio_type gets an array (length nvars) of pio_type values for
 * each var in the file. This array must be freed by caller.
 * @param pio_type_size gets an array (length nvars) of the size of
 * the PIO type for each var in the file. This array must be freed by
 * caller.
 * @param mpi_type gets an array (length nvars) of MPI type values for
 * each var in the file. This array must be freed by caller.
 * @param mpi_type_size gets an array (length nvars) of the size of
 * the MPI type for each var in the file. This array must be freed by
 * caller.
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_openfile_c
 * @author Ed Hartnett
 */
int
inq_file_metadata(file_desc_t *file, int ncid, int iotype, int *nvars, int **rec_var,
                  int **pio_type, int **pio_type_size, MPI_Datatype **mpi_type, int **mpi_type_size)
{
    int nunlimdims = 0;        /* The number of unlimited dimensions. */
    int unlimdimid;
    int *unlimdimids;
    int mpierr;
    int ret;

    if (iotype == PIO_IOTYPE_PNETCDF)
    {
#ifdef _PNETCDF
        if ((ret = ncmpi_inq_nvars(ncid, nvars)))
            return pio_err(NULL, file, PIO_ENOMEM, __FILE__, __LINE__);
#endif /* _PNETCDF */
    }
    else
    {
        if ((ret = nc_inq_nvars(ncid, nvars)))
            return pio_err(NULL, file, PIO_ENOMEM, __FILE__, __LINE__);
    }

    if (*nvars)
    {
        if (!(*rec_var = malloc(*nvars * sizeof(int))))
            return PIO_ENOMEM;
        if (!(*pio_type = malloc(*nvars * sizeof(int))))
            return PIO_ENOMEM;
        if (!(*pio_type_size = malloc(*nvars * sizeof(int))))
            return PIO_ENOMEM;
        if (!(*mpi_type = malloc(*nvars * sizeof(MPI_Datatype))))
            return PIO_ENOMEM;
        if (!(*mpi_type_size = malloc(*nvars * sizeof(int))))
            return PIO_ENOMEM;
    }

    /* How many unlimited dims for this file? */
    if (iotype == PIO_IOTYPE_PNETCDF)
    {
#ifdef _PNETCDF
        if ((ret = ncmpi_inq_unlimdim(ncid, &unlimdimid)))
            return pio_err(NULL, file, ret, __FILE__, __LINE__);
        nunlimdims = unlimdimid == -1 ? 0 : 1;
#endif /* _PNETCDF */
    }
    else if (iotype == PIO_IOTYPE_NETCDF)
    {
        if ((ret = nc_inq_unlimdim(ncid, &unlimdimid)))
            return pio_err(NULL, file, ret, __FILE__, __LINE__);
        nunlimdims = unlimdimid == -1 ? 0 : 1;
    }
    else
    {
#ifdef _NETCDF4
        if ((ret = nc_inq_unlimdims(ncid, &nunlimdims, NULL)))
            return pio_err(NULL, file, ret, __FILE__, __LINE__);
#endif /* _NETCDF4 */
    }

    /* Learn the unlimited dimension ID(s), if there are any. */
    if (nunlimdims)
    {
        if (!(unlimdimids = malloc(nunlimdims * sizeof(int))))
            return pio_err(NULL, file, PIO_ENOMEM, __FILE__, __LINE__);
        if (iotype == PIO_IOTYPE_PNETCDF || iotype == PIO_IOTYPE_NETCDF)
        {
            unlimdimids[0] = unlimdimid;
        }
        else
        {
#ifdef _NETCDF4
            if ((ret = nc_inq_unlimdims(ncid, NULL, unlimdimids)))
                return pio_err(NULL, file, ret, __FILE__, __LINE__);
#endif /* _NETCDF4 */
        }
    }

    /* Learn about each variable in the file. */
    for (int v = 0; v < *nvars; v++)
    {
        int var_ndims;   /* Number of dims for this var. */
        nc_type my_type;

        /* Find type of the var and number of dims in this var. Also
         * learn about type. */
        if (iotype == PIO_IOTYPE_PNETCDF)
        {
#ifdef _PNETCDF
            PIO_Offset type_size;

            if ((ret = ncmpi_inq_var(ncid, v, NULL, &my_type, &var_ndims, NULL, NULL)))
                return pio_err(NULL, file, ret, __FILE__, __LINE__);
            (*pio_type)[v] = (int)my_type;
            if ((ret = pioc_pnetcdf_inq_type(ncid, (*pio_type)[v], NULL, &type_size)))
                return check_netcdf(file, ret, __FILE__, __LINE__);
            (*pio_type_size)[v] = type_size;
#endif /* _PNETCDF */
        }
        else
        {
            size_t type_size;

            if ((ret = nc_inq_var(ncid, v, NULL, &my_type, &var_ndims, NULL, NULL)))
                return pio_err(NULL, file, ret, __FILE__, __LINE__);
            (*pio_type)[v] = (int)my_type;
            if ((ret = nc_inq_type(ncid, (*pio_type)[v], NULL, &type_size)))
                return check_netcdf(file, ret, __FILE__, __LINE__);
            (*pio_type_size)[v] = type_size;
        }

        /* Get the MPI type corresponding with the PIO type. */
        if ((ret = find_mpi_type((*pio_type)[v], &(*mpi_type)[v], NULL)))
            return pio_err(NULL, file, ret, __FILE__, __LINE__);

        /* Get the size of the MPI type. */
        if ((*mpi_type)[v] == MPI_DATATYPE_NULL)
            (*mpi_type_size)[v] = 0;
        else
            if ((mpierr = MPI_Type_size((*mpi_type)[v], &(*mpi_type_size)[v])))
                return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);

        /* What are the dimids associated with this var? */
        if (var_ndims)
        {
            int var_dimids[var_ndims];
            if (iotype == PIO_IOTYPE_PNETCDF)
            {
#ifdef _PNETCDF
                if ((ret = ncmpi_inq_vardimid(ncid, v, var_dimids)))
                    return pio_err(NULL, file, ret, __FILE__, __LINE__);
#endif /* _PNETCDF */
            }
            else
            {
                if ((ret = nc_inq_vardimid(ncid, v, var_dimids)))
                    return pio_err(NULL, file, ret, __FILE__, __LINE__);
            }

            /* Check against each variable dimid agains each unlimited
             * dimid. */
            for (int d = 0; d < var_ndims; d++)
            {
                int unlim_found = 0;

                /* Check against each unlimited dimid. */
                for (int ud = 0; ud < nunlimdims; ud++)
                {
                    if (var_dimids[d] == unlimdimids[ud])
                    {
                        unlim_found++;
                        break;
                    }
                }

                /* Only first dim may be unlimited, for PIO. */
                if (unlim_found)
                {
                    if (d == 0)
                        (*rec_var)[v] = 1;
                    else
                        return pio_err(NULL, file, PIO_EINVAL, __FILE__, __LINE__);

                }
                else
                    (*rec_var)[v] = 0;

            }
        }
    } /* next var */

    /* Free resources. */
    if (nunlimdims)
        free(unlimdimids);

    return PIO_NOERR;
}

/**
 * Find the appropriate IOTYPE from mode flags to nc_open().
 *
 * @param mode the mode flag from nc_open().
 * @param iotype pointer that gets the IOTYPE.
 *
 * @return 0 on success, error code otherwise.
 * @author Ed Hartnett
 */
int
find_iotype_from_omode(int mode, int *iotype)
{
    /* Check inputs. */
    pioassert(iotype, "pointer to iotype must be provided", __FILE__, __LINE__);

    /* Figure out the iotype. */
    if (mode & NC_NETCDF4)
    {
        if (mode & NC_MPIIO || mode & NC_MPIPOSIX)
            *iotype = PIO_IOTYPE_NETCDF4P;
        else
            *iotype = PIO_IOTYPE_NETCDF4C;
    }
    else
    {
        if (mode & NC_PNETCDF || mode & NC_MPIIO)
            *iotype = PIO_IOTYPE_PNETCDF;
        else
            *iotype = PIO_IOTYPE_NETCDF;
    }

    return PIO_NOERR;
}


/**
 * Find the appropriate IOTYPE from mode flags to nc_create().
 *
 * @param cmode the mode flag from nc_create().
 * @param iotype pointer that gets the IOTYPE.
 *
 * @return 0 on success, error code otherwise.
 * @author Ed Hartnett
 */
int
find_iotype_from_cmode(int cmode, int *iotype)
{
    /* Check inputs. */
    pioassert(iotype, "pointer to iotype must be provided", __FILE__, __LINE__);

    /* Figure out the iotype. */
    if (cmode & NC_NETCDF4)
    {
        if (cmode & NC_MPIIO || cmode & NC_MPIPOSIX)
            *iotype = PIO_IOTYPE_NETCDF4P;
        else
            *iotype = PIO_IOTYPE_NETCDF4C;
    }
    else
    {
        if (cmode & NC_PNETCDF || cmode & NC_MPIIO)
            *iotype = PIO_IOTYPE_PNETCDF;
        else
            *iotype = PIO_IOTYPE_NETCDF;
    }

    return PIO_NOERR;
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
 * @param iosysid a defined pio system descriptor.
 * @param ncidp a pio file descriptor.
 * @param iotype a pio output format.
 * @param filename the filename to open
 * @param mode the netcdf mode for the open operation
 * @param retry non-zero to automatically retry with netCDF serial
 * classic.
 * @param use_ext_ncid non-zero to use an externally assigned ncid
 * (used in the netcdf integration layer).
 *
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_openfile_c
 * @author Jim Edwards, Ed Hartnett
 */
int
PIOc_openfile_retry(int iosysid, int *ncidp, int *iotype, const char *filename,
                    int mode, int retry, int use_ext_ncid)
{
    iosystem_desc_t *ios;      /* Pointer to io system information. */
    file_desc_t *file;         /* Pointer to file information. */
    int imode;                 /* Internal mode val for netcdf4 file open. */
    int nvars = 0;
    int *rec_var = NULL;
    int *pio_type = NULL;
    int *pio_type_size = NULL;
    MPI_Datatype *mpi_type = NULL;
    int *mpi_type_size = NULL;
    int mpierr = MPI_SUCCESS, mpierr2;  /** Return code from MPI function codes. */
    int ierr = PIO_NOERR;      /* Return code from function calls. */

#ifdef USE_MPE
    pio_start_mpe_log(OPEN);
#endif /* USE_MPE */

    /* Get the IO system info from the iosysid. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* User must provide valid input for these parameters. */
    if (!ncidp || !iotype || !filename)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);
    if (*iotype < PIO_IOTYPE_PNETCDF || *iotype > PIO_IOTYPE_NETCDF4P)
        return pio_err(ios, NULL, PIO_EINVAL, __FILE__, __LINE__);

    PLOG((2, "PIOc_openfile_retry iosysid = %d iotype = %d filename = %s mode = %d retry = %d",
          iosysid, *iotype, filename, mode, retry));

    /* Allocate space for the file info. */
    if (!(file = calloc(sizeof(*file), 1)))
        return pio_err(ios, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    /* Fill in some file values. */
    file->fh = -1;
    file->iotype = *iotype;
    file->iosystem = ios;
    file->writable = (mode & PIO_WRITE) ? 1 : 0;

    /* Set to true if this task should participate in IO (only true
     * for one task with netcdf serial files. */
    if (file->iotype == PIO_IOTYPE_NETCDF4P || file->iotype == PIO_IOTYPE_PNETCDF ||
        ios->io_rank == 0)
        file->do_io = 1;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async)
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
                mpierr = MPI_Bcast(&mode, 1, MPI_INT, ios->compmaster, ios->intercomm);
        }

        /* Handle MPI errors. */
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
        switch (file->iotype)
        {
#ifdef _NETCDF4

        case PIO_IOTYPE_NETCDF4P:
#ifdef _MPISERIAL
            ierr = nc_open(filename, mode, &file->fh);
#else
            imode = mode |  NC_MPIIO;
            if ((ierr = nc_open_par(filename, imode, ios->io_comm, ios->info, &file->fh)))
                break;

            /* Check the vars for valid use of unlim dims. */
            if ((ierr = check_unlim_use(file->fh)))
                break;

            if ((ierr = inq_file_metadata(file, file->fh, PIO_IOTYPE_NETCDF4P, &nvars, &rec_var, &pio_type,
                                          &pio_type_size, &mpi_type, &mpi_type_size)))
                break;
            PLOG((2, "PIOc_openfile_retry:nc_open_par filename = %s mode = %d imode = %d ierr = %d",
                  filename, mode, imode, ierr));
#endif
            break;

        case PIO_IOTYPE_NETCDF4C:
            if (ios->io_rank == 0)
            {
                if ((ierr = nc_open(filename, mode, &file->fh)))
                    break;
                /* Check the vars for valid use of unlim dims. */
                if ((ierr = check_unlim_use(file->fh)))
                    break;
                ierr = inq_file_metadata(file, file->fh, PIO_IOTYPE_NETCDF4C, &nvars, &rec_var, &pio_type,
                                         &pio_type_size, &mpi_type, &mpi_type_size);
            }
            break;
#endif /* _NETCDF4 */

        case PIO_IOTYPE_NETCDF:
            if (ios->io_rank == 0)
            {
                if ((ierr = nc_open(filename, mode, &file->fh)))
                    break;
                ierr = inq_file_metadata(file, file->fh, PIO_IOTYPE_NETCDF, &nvars, &rec_var, &pio_type,
                                         &pio_type_size, &mpi_type, &mpi_type_size);
            }
            break;

#ifdef _PNETCDF
        case PIO_IOTYPE_PNETCDF:
            ierr = ncmpi_open(ios->io_comm, filename, mode, ios->info, &file->fh);

            // This should only be done with a file opened to append
            if (ierr == PIO_NOERR && (mode & PIO_WRITE))
            {
                if (ios->iomaster == MPI_ROOT)
                    PLOG((2, "%d Setting IO buffer %ld", __LINE__, pio_buffer_size_limit));
                ierr = ncmpi_buffer_attach(file->fh, pio_buffer_size_limit);
            }
            PLOG((2, "ncmpi_open(%s) : fd = %d", filename, file->fh));

            if (!ierr)
                ierr = inq_file_metadata(file, file->fh, PIO_IOTYPE_PNETCDF, &nvars, &rec_var, &pio_type,
                                         &pio_type_size, &mpi_type, &mpi_type_size);
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
            PLOG((2, "retry error code ierr = %d io_rank %d", ierr, ios->io_rank));
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
                {
                    ierr = nc_open(filename, mode, &file->fh);
                    if (ierr == PIO_NOERR)
                        ierr = inq_file_metadata(file, file->fh, PIO_IOTYPE_NETCDF, &nvars, &rec_var, &pio_type,
                                                 &pio_type_size, &mpi_type, &mpi_type_size);
                }
                else
                    file->do_io = 0;
            }
            PLOG((2, "retry nc_open(%s) : fd = %d, iotype = %d, do_io = %d, ierr = %d",
                  filename, file->fh, file->iotype, file->do_io, ierr));
        }
    }

    /* Broadcast and check the return code. */
    if (ios->ioroot == ios->union_rank)
        PLOG((2, "Bcasting error code ierr %d ios->ioroot %d ios->my_comm %d",
              ierr, ios->ioroot, ios->my_comm));
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
    PLOG((2, "Bcast openfile_retry error code ierr = %d", ierr));

    /* If there was an error, free allocated memory and deal with the error. */
    if (ierr)
    {
        free(file);
        return check_netcdf2(ios, NULL, ierr, __FILE__, __LINE__);
    }

    /* Broadcast writability to all tasks. */
    if ((mpierr = MPI_Bcast(&file->writable, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);

    /* Broadcast some values to all tasks from io root. */
    if (ios->async)
    {
        PLOG((3, "open bcasting pio_next_ncid %d ios->ioroot %d", pio_next_ncid, ios->ioroot));
        if ((mpierr = MPI_Bcast(&pio_next_ncid, 1, MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
    }

    if ((mpierr = MPI_Bcast(&nvars, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);

    /* Non io tasks need to allocate to store info about variables. */
    if (nvars && !rec_var)
    {
        if (!(rec_var = malloc(nvars * sizeof(int))))
            return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);
        if (!(pio_type = malloc(nvars * sizeof(int))))
            return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);
        if (!(pio_type_size = malloc(nvars * sizeof(int))))
            return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);
        if (!(mpi_type = malloc(nvars * sizeof(MPI_Datatype))))
            return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);
        if (!(mpi_type_size = malloc(nvars * sizeof(int))))
            return pio_err(ios, file, PIO_ENOMEM, __FILE__, __LINE__);
    }
    if (nvars)
    {
        if ((mpierr = MPI_Bcast(rec_var, nvars, MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(pio_type, nvars, MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(pio_type_size, nvars, MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(mpi_type, nvars*(int)(sizeof(MPI_Datatype)/sizeof(int)), MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Bcast(mpi_type_size, nvars, MPI_INT, ios->ioroot, ios->my_comm)))
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
    }

    /* With the netCDF integration layer, the ncid is assigned for PIO
     * by the netCDF dispatch layer code. So it is passed in. In
     * normal PIO operation, the ncid is generated here. */
    if (!use_ext_ncid)
    {
        /* Create the ncid that the user will see. This is necessary
         * because otherwise ncids will be reused if files are opened
         * on multiple iosystems. */
        file->pio_ncid = pio_next_ncid++;

        /* Return the PIO ncid to the user. */
        *ncidp = file->pio_ncid;
    }
    else
    {
        /* Use the ncid passed in from the netCDF dispatch code. */
        file->pio_ncid = *ncidp;

        /* To prevent PIO from reusing the same ncid, if someone
         * starting mingling netcdf integration PIO and regular PIO
         * code. */
        pio_next_ncid = file->pio_ncid + 1;
    }

    /* Add this file to the list of currently open files. */
    pio_add_to_file_list(file);

    /* Add info about the variables to the file_desc_t struct. */
    for (int v = 0; v < nvars; v++)
        if ((ierr = add_to_varlist(v, rec_var[v], pio_type[v], pio_type_size[v], mpi_type[v],
                                   mpi_type_size[v], &file->varlist)))
            return pio_err(ios, NULL, ierr, __FILE__, __LINE__);
    file->nvars = nvars;

    /* Free resources. */
    if (nvars)
    {
        if (rec_var)
            free(rec_var);
        if (pio_type)
            free(pio_type);
        if (pio_type_size)
            free(pio_type_size);
        if (mpi_type)
            free(mpi_type);
        if (mpi_type_size)
            free(mpi_type_size);
    }

#ifdef USE_MPE
    pio_stop_mpe_log(OPEN, __func__);
#endif /* USE_MPE */
    PLOG((2, "Opened file %s file->pio_ncid = %d file->fh = %d ierr = %d",
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
 * @author Ed Hartnett
 */
int
pioc_pnetcdf_inq_type(int ncid, nc_type xtype, char *name,
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
 * @returns PIO_NOERR on success, error code on failure.
 * @author Ed Hartnett
 */
int
pioc_change_def(int ncid, int is_enddef)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI functions. */

    PLOG((2, "pioc_change_def ncid = %d is_enddef = %d", ncid, is_enddef));

    /* Find the info about this file. When I check the return code
     * here, some tests fail. ???*/
    if ((ierr = pio_get_file(ncid, &file)))
        return pio_err(NULL, NULL, ierr, __FILE__, __LINE__);
    ios = file->iosystem;

    /* If async is in use, and this is not an IO task, bcast the parameters. */
    if (ios->async)
    {
        if (!ios->ioproc)
        {
            int msg = is_enddef ? PIO_MSG_ENDDEF : PIO_MSG_REDEF;
            if (ios->compmaster == MPI_ROOT)
                mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);

            if (!mpierr)
                mpierr = MPI_Bcast(&ncid, 1, MPI_INT, ios->compmaster, ios->intercomm);
            PLOG((3, "pioc_change_def ncid = %d mpierr = %d", ncid, mpierr));
        }

        /* Handle MPI errors. */
        PLOG((3, "pioc_change_def handling MPI errors"));
        if ((mpierr2 = MPI_Bcast(&mpierr, 1, MPI_INT, ios->comproot, ios->my_comm)))
            check_mpi(NULL, file, mpierr2, __FILE__, __LINE__);
        if (mpierr)
            return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
    }

    /* If this is an IO task, then call the netCDF function. */
    PLOG((3, "pioc_change_def ios->ioproc = %d", ios->ioproc));
    if (ios->ioproc)
    {
        PLOG((3, "pioc_change_def calling netcdf function file->fh = %d file->do_io = %d iotype = %d",
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
                PLOG((3, "pioc_change_def calling nc_enddef file->fh = %d", file->fh));
                ierr = nc_enddef(file->fh);
            }
            else
                ierr = nc_redef(file->fh);
        }
    }

    /* Broadcast and check the return code. */
    PLOG((3, "pioc_change_def bcasting return code ierr = %d", ierr));
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(NULL, file, mpierr, __FILE__, __LINE__);
    if (ierr)
        return check_netcdf(file, ierr, __FILE__, __LINE__);
    PLOG((3, "pioc_change_def succeeded"));

    return ierr;
}

/**
 * Check whether an IO type is valid for the build.
 *
 * @param iotype the IO type to check
 * @returns 0 if valid, non-zero otherwise.
 * @author Jim Edwards
 */
int
iotype_is_valid(int iotype)
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
 * Set the rearranger options associated with an iosystem
 *
 * @param iosysid a defined pio system descriptor.
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
 * @return 0 on success, otherwise a PIO error code.
 * @author Jayesh Krishna
 */
int
PIOc_set_rearr_opts(int iosysid, int comm_type, int fcd, bool enable_hs_c2i,
                    bool enable_isend_c2i, int max_pend_req_c2i,
                    bool enable_hs_i2c, bool enable_isend_i2c,
                    int max_pend_req_i2c)
{
    iosystem_desc_t *ios;
    rearr_opt_t user_rearr_opts = {
        comm_type, fcd,
        {enable_hs_c2i,enable_isend_c2i, max_pend_req_c2i},
        {enable_hs_i2c, enable_isend_i2c, max_pend_req_i2c}
    };

    /* Check inputs. */
    if ((comm_type != PIO_REARR_COMM_P2P && comm_type != PIO_REARR_COMM_FC_1D_COMP2IO) ||
        (fcd < 0 || fcd > PIO_REARR_COMM_FC_2D_DISABLE) ||
        (max_pend_req_c2i != PIO_REARR_COMM_UNLIMITED_PEND_REQ && max_pend_req_c2i < 0) ||
        (max_pend_req_i2c != PIO_REARR_COMM_UNLIMITED_PEND_REQ && max_pend_req_i2c < 0))
        return pio_err(NULL, NULL, PIO_EINVAL, __FILE__, __LINE__);

    /* Get the IO system info. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

    /* Set the options. */
    ios->rearr_opts = user_rearr_opts;

    return PIO_NOERR;
}

/**
 * This function determines which processes are assigned to the
 * different computation components. This function is called by
 * PIOc_init_async().
 *
 * The user may have passed a specification of tasks as array
 * proc_list, or it may be calculated by assigning processors starting
 * at the first one after the IO component, and assigning them in
 * order to each computation component.
 *
 * Note that memory is allocated for my_proc_list. This must be freed
 * by the caller.
 *
 * @param num_io_procs the number of IO processes.
 * @param component_count the number of computational components.
 * @param num_procs_per_comp array (length component_count) which
 * contains the number of processes to assign to each computation
 * component.
 * @param proc_list array (length component count) of arrays (length
 * num_procs_per_comp_array[cmp]) which contain the list of processes
 * for each computation component. May be NULL.
 * @param my_proc_list array (length component count) of arrays (length
 * num_procs_per_comp_array[cmp]) which will get the list of processes
 * for each computation component.
 * @returns 0 for success, error code otherwise
 * @author Ed Hartnett
 */
int
determine_procs(int num_io_procs, int component_count, int *num_procs_per_comp,
                int **proc_list, int **my_proc_list)
{
    /* If the user did not provide a list of processes for each
     * component, create one. */
    if (!proc_list)
    {
        int last_proc = num_io_procs;

        /* Fill the array of arrays. */
        for (int cmp = 0; cmp < component_count; cmp++)
        {
            PLOG((3, "calculating processors for component %d num_procs_per_comp[cmp] = %d",
                  cmp, num_procs_per_comp[cmp]));

            /* Allocate space for each array. */
            if (!(my_proc_list[cmp] = malloc(num_procs_per_comp[cmp] * sizeof(int))))
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

            int proc;
            for (proc = last_proc; proc < num_procs_per_comp[cmp] + last_proc; proc++)
            {
                my_proc_list[cmp][proc - last_proc] = proc;
                PLOG((3, "my_proc_list[%d][%d] = %d", cmp, proc - last_proc, proc));
            }
            last_proc = proc;
        }
    }
    else
    {
        for (int cmp = 0; cmp < component_count; cmp++)
        {
            /* Allocate space for each array. */
            if (!(my_proc_list[cmp] = malloc(num_procs_per_comp[cmp] * sizeof(int))))
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
            memcpy(my_proc_list[cmp], proc_list[cmp], num_procs_per_comp[cmp] * sizeof(int));
        }
    }
    return PIO_NOERR;
}
