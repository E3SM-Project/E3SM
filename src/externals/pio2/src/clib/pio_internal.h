/**
 * @file
 * Private headers and defines for the PIO C interface.
 * @author Jim Edwards
 * @date  2014
 *
 * @see http://code.google.com/p/parallelio/
 */

#ifndef __PIO_INTERNAL__
#define __PIO_INTERNAL__

#include <pio.h>

/* It seems that some versions of openmpi fail to define
 * MPI_OFFSET. */
#ifdef OMPI_OFFSET_DATATYPE
#ifndef MPI_OFFSET
#define MPI_OFFSET OMPI_OFFSET_DATATYPE
#endif
#endif
#ifndef MPI_Offset
#define MPI_Offset long long
#endif

#if defined(MPT_VERSION) || defined(OPEN_MPI)
/* Some MPI implementations do not allow passing MPI_DATATYPE_NULL to comm functions
 * even though the send or recv length is 0, in these cases we use MPI_CHAR */
#define PIO_DATATYPE_NULL MPI_CHAR
#else
#define PIO_DATATYPE_NULL MPI_DATATYPE_NULL
#endif

#include <bget.h>
#include <limits.h>
#include <math.h>
#ifdef TIMING
#include <gptl.h>
#endif
#include <assert.h>

#if PIO_ENABLE_LOGGING
void pio_log(int severity, const char *fmt, ...);
#define LOG(e) pio_log e
#else
#define LOG(e)
#endif /* PIO_ENABLE_LOGGING */

#define max(a,b)                                \
    ({ __typeof__ (a) _a = (a);                 \
        __typeof__ (b) _b = (b);                \
        _a > _b ? _a : _b; })

#define min(a,b)                                \
    ({ __typeof__ (a) _a = (a);                 \
        __typeof__ (b) _b = (b);                \
        _a < _b ? _a : _b; })

#define MAX_GATHER_BLOCK_SIZE 0
#define PIO_REQUEST_ALLOC_CHUNK 16

/** This is needed to handle _long() functions. It may not be used as
 * a data type when creating attributes or varaibles, it is only used
 * internally. */
#define PIO_LONG_INTERNAL 13

#if defined(__cplusplus)
extern "C" {
#endif

    extern PIO_Offset pio_buffer_size_limit;

    /** Used to sort map points in the subset rearranger. */
    typedef struct mapsort
    {
        int rfrom;
        PIO_Offset soffset;
        PIO_Offset iomap;
    } mapsort;

    /** swapm defaults. */
    typedef struct pio_swapm_defaults
    {
        int nreqs;
        bool handshake;
        bool isend;
    } pio_swapm_defaults;

    /* Handle an error in the PIO library. */
    int pio_err(iosystem_desc_t *ios, file_desc_t *file, int err_num, const char *fname,
                int line);

    /* For async cases, this runs on IO tasks and listens for messages. */
    int pio_msg_handler2(int io_rank, int component_count, iosystem_desc_t **iosys,
                         MPI_Comm io_comm);

    void pio_get_env(void);
    int  pio_add_to_iodesc_list(io_desc_t *iodesc);
    io_desc_t *pio_get_iodesc_from_id(int ioid);
    int pio_delete_iodesc_from_list(int ioid);
    int pio_num_iosystem(int *niosysid);

    int pio_get_file(int ncid, file_desc_t **filep);
    int pio_delete_file_from_list(int ncid);
    void pio_add_to_file_list(file_desc_t *file);
    void pio_push_request(file_desc_t *file, int request);

    /* Create a file (internal function). */
    int PIOc_createfile_int(int iosysid, int *ncidp, int *iotype, const char *filename, int mode);

    /* Open a file with optional retry as netCDF-classic if first
     * iotype does not work. */
    int PIOc_openfile_retry(int iosysid, int *ncidp, int *iotype,
                            const char *filename, int mode, int retry);

    iosystem_desc_t *pio_get_iosystem_from_id(int iosysid);
    int pio_add_to_iosystem_list(iosystem_desc_t *ios);

    /* Check the return code from a netCDF call. */
    int check_netcdf(file_desc_t *file, int status, const char *fname, int line);

    /* Check the return code from a netCDF call, with ios pointer. */
    int check_netcdf2(iosystem_desc_t *ios, file_desc_t *file, int status,
                      const char *fname, int line);

    /* Find the MPI type that matches a PIO type. */
    int find_mpi_type(int pio_type, MPI_Datatype *mpi_type);

    /* Check whether an IO type is valid for this build. */
    int iotype_is_valid(int iotype);

    /* Print error message and abort. */
    void piodie(const char *msg, const char *fname, int line);

    /* Assert that an expression is true. */
    void pioassert(bool exp, const char *msg, const char *fname, int line);

    int CalcStartandCount(int basetype, int ndims, const int *gdims, int num_io_procs,
                          int myiorank, PIO_Offset *start, PIO_Offset *kount);

    /* Check return from MPI function and print error message. */
    void CheckMPIReturn(int ierr, const char *file, int line);

    /* Like MPI_Alltoallw(), but with flow control. */
    int pio_swapm(void *sndbuf, int *sndlths, int *sdispls, MPI_Datatype *stypes,
                  void *rcvbuf, int *rcvlths, int *rdispls, MPI_Datatype *rtypes,
                  MPI_Comm comm, bool handshake, bool isend, int max_requests);

    long long lgcd_array(int nain, long long* ain);

    void PIO_Offset_size(MPI_Datatype *dtype, int *tsize);
    PIO_Offset GCDblocksize(int arrlen, const PIO_Offset *arr_in);

    /* Initialize the rearranger options. */
    void init_rearr_opts(iosystem_desc_t *iosys);

    /* Compare sets of rearranger options. */
    bool cmp_rearr_opts(const rearr_opt_t *rearr_opts, const rearr_opt_t *exp_rearr_opts);

    /* Reset rearranger opts in iosystem to valid values. */
    void check_and_reset_rearr_opts(iosystem_desc_t *ios);

    /* Compare rearranger flow control options. */
    bool cmp_rearr_comm_fc_opts(const rearr_comm_fc_opt_t *opt,
                                const rearr_comm_fc_opt_t *exp_opt);

    /* Create a subset rearranger. */
    int subset_rearrange_create(iosystem_desc_t *ios, int maplen, PIO_Offset *compmap, const int *gsize,
                                int ndim, io_desc_t *iodesc);


    /* Create a box rearranger. */
    int box_rearrange_create(iosystem_desc_t *ios, int maplen, const PIO_Offset *compmap, const int *gsize,
                             int ndim, io_desc_t *iodesc);


    /* Move data from IO tasks to compute tasks. */
    int rearrange_io2comp(iosystem_desc_t *ios, io_desc_t *iodesc, void *sbuf, void *rbuf);

    /* Move data from compute tasks to IO tasks. */
    int rearrange_comp2io(iosystem_desc_t *ios, io_desc_t *iodesc, void *sbuf, void *rbuf, int nvars);

    /* Allocate and initialize storage for decomposition information. */
    int malloc_iodesc(iosystem_desc_t *ios, int piotype, int ndims, io_desc_t **iodesc);
    void performance_tune_rearranger(iosystem_desc_t *ios, io_desc_t *iodesc);

    int flush_output_buffer(file_desc_t *file, bool force, PIO_Offset addsize);
    int compute_maxIObuffersize(MPI_Comm io_comm, io_desc_t *iodesc);
    io_region *alloc_region(int ndims);
    int pio_delete_iosystem_from_list(int piosysid);

    /* Find greatest commond divisor. */
    int gcd(int a, int b);

    /* Find greatest commond divisor for long long. */
    long long lgcd (long long a, long long b );

    /* Find greatest commond divisor in an array. */
    int gcd_array(int nain, int *ain);

    void free_region_list(io_region *top);

    /* Convert a global coordinate value into a local array index. */
    PIO_Offset coord_to_lindex(int ndims, const PIO_Offset *lcoord, const PIO_Offset *count);

    int ceil2(int i);
    int pair(int np, int p, int k);

    /* Create MPI datatypes used for comp2io and io2comp data transfers. */
    int define_iodesc_datatypes(iosystem_desc_t *ios, io_desc_t *iodesc);

    /* Create the derived MPI datatypes used for comp2io and io2comp
     * transfers. */
    int create_mpi_datatypes(MPI_Datatype basetype, int msgcnt, PIO_Offset dlen, const PIO_Offset *mindex,
                             const int *mcount, int *mfrom, MPI_Datatype *mtype);
    int compare_offsets(const void *a, const void *b) ;

    /* Print a trace statement, for debugging. */
    void print_trace (FILE *fp);

    void cn_buffer_report(iosystem_desc_t *ios, bool collective);

    /* Initialize the compute buffer. */
    int compute_buffer_init(iosystem_desc_t *ios);

    void free_cn_buffer_pool(iosystem_desc_t *ios);

    /* Flush PIO's data buffer. */
    int flush_buffer(int ncid, wmulti_buffer *wmb, bool flushtodisk);

    int compute_maxaggregate_bytes(iosystem_desc_t *ios, io_desc_t *iodesc);

    /* Announce a memory error with bget memory, and die. */
    void piomemerror(iosystem_desc_t *ios, size_t req, char *fname, int line);

    /* Check the return code from an MPI function call. */
    int check_mpi(file_desc_t *file, int mpierr, const char *filename, int line);

    /* Check the return code from an MPI function call. */
    int check_mpi2(iosystem_desc_t *ios, file_desc_t *file, int mpierr, const char *filename,
                   int line);

    /* Darray support functions. */

    /* Write aggregated arrays to file using parallel I/O (netCDF-4 parallel/pnetcdf) */
    int pio_write_darray_multi_nc(file_desc_t *file, int nvars, const int *vid, int iodesc_ndims,
                                  MPI_Datatype basetype, int maxregions, io_region *firstregion,
                                  PIO_Offset llen, int num_aiotasks, void *iobuf,
                                  const int *frame);

    /* Write aggregated arrays to file using serial I/O (netCDF-3/netCDF-4 serial) */
    int pio_write_darray_multi_nc_serial(file_desc_t *file, int nvars, const int *vid, int iodesc_ndims,
                                         MPI_Datatype basetype, int maxregions, io_region *firstregion,
                                         PIO_Offset llen, int num_aiotasks, void *iobuf,
                                         const int *frame);

    int pio_read_darray_nc(file_desc_t *file, io_desc_t *iodesc, int vid, void *iobuf);
    int pio_read_darray_nc_serial(file_desc_t *file, io_desc_t *iodesc, int vid, void *iobuf);

    /* Read atts with type conversion. */
    int PIOc_get_att_tc(int ncid, int varid, const char *name, nc_type memtype, void *ip);

    /* Write atts with type conversion. */
    int PIOc_put_att_tc(int ncid, int varid, const char *name, nc_type atttype,
                        PIO_Offset len, nc_type memtype, const void *op);

    /* Generalized get functions. */
    int PIOc_get_vars_tc(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
                         const PIO_Offset *stride, nc_type xtype, void *buf);
    int PIOc_get_var1_tc(int ncid, int varid, const PIO_Offset *index, nc_type xtype,
                         void *buf);

    /* Generalized put functions. */
    int PIOc_put_vars_tc(int ncid, int varid, const PIO_Offset *start, const PIO_Offset *count,
                         const PIO_Offset *stride, nc_type xtype, const void *buf);
    int PIOc_put_var1_tc(int ncid, int varid, const PIO_Offset *index, nc_type xtype,
                         const void *op);

    /* An internal replacement for a function pnetcdf does not
     * have. */
    int pioc_pnetcdf_inq_type(int ncid, nc_type xtype, char *name,
                              PIO_Offset *sizep);

    /* Handle end and re-defs. */
    int pioc_change_def(int ncid, int is_enddef);

    /* Initialize and finalize logging. */
    void pio_init_logging(void);
    void pio_finalize_logging(void );

    /* Write a netCDF decomp file. */
    int pioc_write_nc_decomp_int(int iosysid, const char *filename, int cmode, int ndims, int *global_dimlen,
                                 int num_tasks, int *task_maplen, int *map, const char *title,
                                 const char *history, int fortran_order);

    /* Read a netCDF decomp file. */
    int pioc_read_nc_decomp_int(int iosysid, const char *filename, int *ndims, int **global_dimlen,
                                int *num_tasks, int **task_maplen, int *max_maplen, int **map, char *title,
                                char *history, char *source, char *version, int *fortran_order);

#if defined(__cplusplus)
}
#endif

/** These are the messages that can be sent over the intercomm when
 * async is being used. */
enum PIO_MSG
{
    PIO_MSG_OPEN_FILE,
    PIO_MSG_CREATE_FILE,
    PIO_MSG_INQ_ATT,
    PIO_MSG_INQ_FORMAT,
    PIO_MSG_INQ_VARID,
    PIO_MSG_DEF_VAR,
    PIO_MSG_INQ_VAR,
    PIO_MSG_PUT_ATT_DOUBLE,
    PIO_MSG_PUT_ATT_INT,
    PIO_MSG_RENAME_ATT,
    PIO_MSG_DEL_ATT,
    PIO_MSG_INQ,
    PIO_MSG_GET_ATT_TEXT,
    PIO_MSG_GET_ATT_SHORT,
    PIO_MSG_PUT_ATT_LONG,
    PIO_MSG_REDEF,
    PIO_MSG_SET_FILL,
    PIO_MSG_ENDDEF,
    PIO_MSG_RENAME_VAR,
    PIO_MSG_PUT_ATT_SHORT,
    PIO_MSG_PUT_ATT_TEXT,
    PIO_MSG_INQ_ATTNAME,
    PIO_MSG_GET_ATT_ULONGLONG,
    PIO_MSG_GET_ATT_USHORT,
    PIO_MSG_PUT_ATT_ULONGLONG,
    PIO_MSG_GET_ATT_UINT,
    PIO_MSG_GET_ATT_LONGLONG,
    PIO_MSG_PUT_ATT_SCHAR,
    PIO_MSG_PUT_ATT_FLOAT,
    PIO_MSG_RENAME_DIM,
    PIO_MSG_GET_ATT_LONG,
    PIO_MSG_INQ_DIM,
    PIO_MSG_INQ_DIMID,
    PIO_MSG_PUT_ATT_USHORT,
    PIO_MSG_GET_ATT_FLOAT,
    PIO_MSG_SYNC,
    PIO_MSG_PUT_ATT_LONGLONG,
    PIO_MSG_PUT_ATT_UINT,
    PIO_MSG_GET_ATT_SCHAR,
    PIO_MSG_INQ_ATTID,
    PIO_MSG_DEF_DIM,
    PIO_MSG_GET_ATT_INT,
    PIO_MSG_GET_ATT_DOUBLE,
    PIO_MSG_PUT_ATT_UCHAR,
    PIO_MSG_GET_ATT_UCHAR,
    PIO_MSG_PUT_VARS_UCHAR,
    PIO_MSG_GET_VAR1_SCHAR,
    PIO_MSG_GET_VARS_ULONGLONG,
    PIO_MSG_GET_VARM_UCHAR,
    PIO_MSG_GET_VARM_SCHAR,
    PIO_MSG_GET_VARS_SHORT,
    PIO_MSG_GET_VAR_DOUBLE,
    PIO_MSG_GET_VARA_DOUBLE,
    PIO_MSG_GET_VAR_INT,
    PIO_MSG_GET_VAR_USHORT,
    PIO_MSG_PUT_VARS_USHORT,
    PIO_MSG_GET_VARA_TEXT,
    PIO_MSG_PUT_VARS_ULONGLONG,
    PIO_MSG_GET_VARA_INT,
    PIO_MSG_PUT_VARM,
    PIO_MSG_GET_VAR1_FLOAT,
    PIO_MSG_GET_VAR1_SHORT,
    PIO_MSG_GET_VARS_INT,
    PIO_MSG_PUT_VARS_UINT,
    PIO_MSG_GET_VAR_TEXT,
    PIO_MSG_GET_VARM_DOUBLE,
    PIO_MSG_PUT_VARM_UCHAR,
    PIO_MSG_PUT_VAR_USHORT,
    PIO_MSG_GET_VARS_SCHAR,
    PIO_MSG_GET_VARA_USHORT,
    PIO_MSG_PUT_VAR1_LONGLONG,
    PIO_MSG_PUT_VARA_UCHAR,
    PIO_MSG_PUT_VARM_SHORT,
    PIO_MSG_PUT_VAR1_LONG,
    PIO_MSG_PUT_VARS_LONG,
    PIO_MSG_GET_VAR1_USHORT,
    PIO_MSG_PUT_VAR_SHORT,
    PIO_MSG_PUT_VARA_INT,
    PIO_MSG_GET_VAR_FLOAT,
    PIO_MSG_PUT_VAR1_USHORT,
    PIO_MSG_PUT_VARA_TEXT,
    PIO_MSG_PUT_VARM_TEXT,
    PIO_MSG_GET_VARS_UCHAR,
    PIO_MSG_GET_VAR,
    PIO_MSG_PUT_VARM_USHORT,
    PIO_MSG_GET_VAR1_LONGLONG,
    PIO_MSG_GET_VARS_USHORT,
    PIO_MSG_GET_VAR_LONG,
    PIO_MSG_GET_VAR1_DOUBLE,
    PIO_MSG_PUT_VAR_ULONGLONG,
    PIO_MSG_PUT_VAR_INT,
    PIO_MSG_GET_VARA_UINT,
    PIO_MSG_PUT_VAR_LONGLONG,
    PIO_MSG_GET_VARS_LONGLONG,
    PIO_MSG_PUT_VAR_SCHAR,
    PIO_MSG_PUT_VAR_UINT,
    PIO_MSG_PUT_VAR,
    PIO_MSG_PUT_VARA_USHORT,
    PIO_MSG_GET_VAR_LONGLONG,
    PIO_MSG_GET_VARA_SHORT,
    PIO_MSG_PUT_VARS_SHORT,
    PIO_MSG_PUT_VARA_UINT,
    PIO_MSG_PUT_VARA_SCHAR,
    PIO_MSG_PUT_VARM_ULONGLONG,
    PIO_MSG_PUT_VAR1_UCHAR,
    PIO_MSG_PUT_VARM_INT,
    PIO_MSG_PUT_VARS_SCHAR,
    PIO_MSG_GET_VARA_LONG,
    PIO_MSG_PUT_VAR1,
    PIO_MSG_GET_VAR1_INT,
    PIO_MSG_GET_VAR1_ULONGLONG,
    PIO_MSG_GET_VAR_UCHAR,
    PIO_MSG_PUT_VARA_FLOAT,
    PIO_MSG_GET_VARA_UCHAR,
    PIO_MSG_GET_VARS_FLOAT,
    PIO_MSG_PUT_VAR1_FLOAT,
    PIO_MSG_PUT_VARM_FLOAT,
    PIO_MSG_PUT_VAR1_TEXT,
    PIO_MSG_PUT_VARS_TEXT,
    PIO_MSG_PUT_VARM_LONG,
    PIO_MSG_GET_VARS_LONG,
    PIO_MSG_PUT_VARS_DOUBLE,
    PIO_MSG_GET_VAR1,
    PIO_MSG_GET_VAR_UINT,
    PIO_MSG_PUT_VARA_LONGLONG,
    PIO_MSG_GET_VARA,
    PIO_MSG_PUT_VAR_DOUBLE,
    PIO_MSG_GET_VARA_SCHAR,
    PIO_MSG_PUT_VAR_FLOAT,
    PIO_MSG_GET_VAR1_UINT,
    PIO_MSG_GET_VARS_UINT,
    PIO_MSG_PUT_VAR1_ULONGLONG,
    PIO_MSG_PUT_VARM_UINT,
    PIO_MSG_PUT_VAR1_UINT,
    PIO_MSG_PUT_VAR1_INT,
    PIO_MSG_GET_VARA_FLOAT,
    PIO_MSG_GET_VARM_TEXT,
    PIO_MSG_PUT_VARS_FLOAT,
    PIO_MSG_GET_VAR1_TEXT,
    PIO_MSG_PUT_VARA_SHORT,
    PIO_MSG_PUT_VAR1_SCHAR,
    PIO_MSG_PUT_VARA_ULONGLONG,
    PIO_MSG_PUT_VARM_DOUBLE,
    PIO_MSG_GET_VARM_INT,
    PIO_MSG_PUT_VARA,
    PIO_MSG_PUT_VARA_LONG,
    PIO_MSG_GET_VARM_UINT,
    PIO_MSG_GET_VARM,
    PIO_MSG_PUT_VAR1_DOUBLE,
    PIO_MSG_GET_VARS_DOUBLE,
    PIO_MSG_GET_VARA_LONGLONG,
    PIO_MSG_GET_VAR_ULONGLONG,
    PIO_MSG_PUT_VARM_SCHAR,
    PIO_MSG_GET_VARA_ULONGLONG,
    PIO_MSG_GET_VAR_SHORT,
    PIO_MSG_GET_VARM_FLOAT,
    PIO_MSG_PUT_VAR_TEXT,
    PIO_MSG_PUT_VARS_INT,
    PIO_MSG_GET_VAR1_LONG,
    PIO_MSG_GET_VARM_LONG,
    PIO_MSG_GET_VARM_USHORT,
    PIO_MSG_PUT_VAR1_SHORT,
    PIO_MSG_PUT_VARS_LONGLONG,
    PIO_MSG_GET_VARM_LONGLONG,
    PIO_MSG_GET_VARS_TEXT,
    PIO_MSG_PUT_VARA_DOUBLE,
    PIO_MSG_PUT_VARS,
    PIO_MSG_PUT_VAR_UCHAR,
    PIO_MSG_GET_VAR1_UCHAR,
    PIO_MSG_PUT_VAR_LONG,
    PIO_MSG_GET_VARS,
    PIO_MSG_GET_VARM_SHORT,
    PIO_MSG_GET_VARM_ULONGLONG,
    PIO_MSG_PUT_VARM_LONGLONG,
    PIO_MSG_GET_VAR_SCHAR,
    PIO_MSG_GET_ATT_UBYTE,
    PIO_MSG_PUT_ATT_STRING,
    PIO_MSG_GET_ATT_STRING,
    PIO_MSG_PUT_ATT_UBYTE,
    PIO_MSG_INQ_VAR_FILL,
    PIO_MSG_DEF_VAR_FILL,
    PIO_MSG_DEF_VAR_DEFLATE,
    PIO_MSG_INQ_VAR_DEFLATE,
    PIO_MSG_INQ_VAR_SZIP,
    PIO_MSG_DEF_VAR_FLETCHER32,
    PIO_MSG_INQ_VAR_FLETCHER32,
    PIO_MSG_DEF_VAR_CHUNKING,
    PIO_MSG_INQ_VAR_CHUNKING,
    PIO_MSG_DEF_VAR_ENDIAN,
    PIO_MSG_INQ_VAR_ENDIAN,
    PIO_MSG_SET_CHUNK_CACHE,
    PIO_MSG_GET_CHUNK_CACHE,
    PIO_MSG_SET_VAR_CHUNK_CACHE,
    PIO_MSG_GET_VAR_CHUNK_CACHE,
    PIO_MSG_INITDECOMP_DOF,
    PIO_MSG_WRITEDARRAY,
    PIO_MSG_READDARRAY,
    PIO_MSG_SETERRORHANDLING,
    PIO_MSG_FREEDECOMP,
    PIO_MSG_CLOSE_FILE,
    PIO_MSG_DELETE_FILE,
    PIO_MSG_EXIT,
    PIO_MSG_GET_ATT,
    PIO_MSG_PUT_ATT,
    PIO_MSG_INQ_TYPE
};

#endif /* __PIO_INTERNAL__ */
