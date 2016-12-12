/**
 * @file 
 * Public headers for the PIO C interface.
 * @author Jim Edwards
 * @date  2014
 * 
 * @see http://code.google.com/p/parallelio/
 */

#ifndef _PIO_H_
#define _PIO_H_
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h> // memcpy
#include <mpi.h>

#ifdef _NETCDF
#include <netcdf.h>
#ifdef _NETCDF4
#include <netcdf_par.h>
#endif
#endif
#ifdef _PNETCDF
#include <pnetcdf.h>
#endif

#ifndef MPI_OFFSET
/** MPI_OFFSET is an integer type of size sufficient to represent the
 * size (in bytes) of the largest file supported by MPI. In some MPI
 * implementations MPI_OFFSET is not properly defined.  */
#define MPI_OFFSET  MPI_LONG_LONG
#endif

/** PIO_OFFSET is an integer type of size sufficient to represent the
 * size (in bytes) of the largest file supported by MPI. */
#define PIO_OFFSET MPI_OFFSET
#define PIO_Offset MPI_Offset

/** The maximum number of variables allowed in a netCDF file. */
#define PIO_MAX_VARS NC_MAX_VARS

/**
 * Variable description structure.
 */
typedef struct var_desc_t
{
    /** The unlimited dimension in the netCDF file (typically the time
     * dimension). -1 if there is no unlimited dimension. */
    int record;

    /** Number of dimensions for this variable. */
    int ndims;

    /** ID of each outstanding pnetcdf request for this variable. */
    int *request;

    /** Number of requests bending with pnetcdf. */
    int nreqs;

    /** Buffer that contains the fill value for this variable. */
    void *fillbuf;

    /** ??? */
    void *iobuf;
} var_desc_t;

/**
 * IO region structure.
 *
 * Each IO region is a unit of data which can be described using start and count
 * arrays. Each IO task may in general have multiple io regions per variable.  The 
 * box rearranger will have at most one io region per variable.
 * 
 */
typedef struct io_region
{
    int loffset;
    PIO_Offset *start;
    PIO_Offset *count;
    struct io_region *next;
} io_region;
 
/**
 * IO descriptor structure.
 *
 * This structure defines the mapping for a given variable between 
 * compute and IO decomposition.   
 */
typedef struct io_desc_t
{
    /** The ID of this io_desc_t. */
    int ioid;
    int async_id;
    int nrecvs;
    int ndof;
    int ndims;
    int num_aiotasks;
    int rearranger;
    int maxregions;
    bool needsfill;       // Does this decomp leave holes in the field (true) or write everywhere (false)

    /** The maximum number of bytes of this iodesc before flushing. */
    int maxbytes;        
    MPI_Datatype basetype;
    PIO_Offset llen;
    int maxiobuflen;
    PIO_Offset *gsize;

    int *rfrom;
    int *rcount;
    int *scount;
    PIO_Offset *sindex;
    PIO_Offset *rindex;

    MPI_Datatype *rtype;
    MPI_Datatype *stype;
    int num_stypes;
    int holegridsize;
    int maxfillregions;
    io_region *firstregion;
    io_region *fillregion;

    bool handshake;
    bool isend;
    int max_requests;
  
    MPI_Comm subset_comm;

    /** Pointer to the next io_desc_t in the list. */
    struct io_desc_t *next;
} io_desc_t;

/**
 * IO system descriptor structure.
 *
 * This structure contains the general IO subsystem data and MPI
 * structure
 */
typedef struct iosystem_desc_t 
{
    /** The ID of this iosystem_desc_t. This will be obtained by
     * calling PIOc_Init_Intercomm() or PIOc_Init_Intracomm(). */
    int iosysid;

    /** This is an MPI intra communicator that includes all the tasks in
     * both the IO and the computation communicators. */
    MPI_Comm union_comm;

    /** This is an MPI intra communicator that includes all the tasks
     * involved in IO. */
    MPI_Comm io_comm;

    /** This is an MPI intra communicator that includes all the tasks
     * involved in computation. */
    MPI_Comm comp_comm;

    /** This is an MPI inter communicator between IO communicator and
     * computation communicator. */
    MPI_Comm intercomm;

    /** This is a copy (but not an MPI copy) of either the comp (for
     * non-async) or the union (for async) communicator. */
    MPI_Comm my_comm;

    /** This MPI group contains the processors involved in
     * computation. */
    MPI_Group compgroup;
    
    /** This MPI group contains the processors involved in I/O. */
    MPI_Group iogroup;

    /** The number of tasks in the IO communicator. */
    int num_iotasks;

    /** The number of tasks in the computation communicator. */
    int num_comptasks;

    /** Rank of this task in the union communicator. */
    int union_rank;

    /** The rank of this process in the computation communicator, or -1
     * if this process is not part of the computation communicator. */
    int comp_rank;

    /** The rank of this process in the IO communicator, or -1 if this
     * process is not part of the IO communicator. */
    int io_rank;

    /** Set to MPI_ROOT if this task is the master of IO communicator, 0
     * otherwise. */
    int iomaster;

    /** Set to MPI_ROOT if this task is the master of comp communicator, 0
     * otherwise. */
    int compmaster;

    /** Rank of IO root task (which is rank 0 in io_comm) in the union
     * communicator. */
    int ioroot;

    /** Rank of computation root task (which is rank 0 in
     * comm_comms[cmp]) in the union communicator. */
    int comproot;

    /** An array of the ranks of all IO tasks within the union
     * communicator. */
    int *ioranks;

    /** Controls handling errors. */
    int error_handler;

    /** The rearranger decides which parts of a distributed array are
     * handled by which IO tasks. */
    int default_rearranger;

    /** True if asynchronous interface is in use. */
    bool async_interface;

    /** True if this task is a member of the IO communicator. */
    bool ioproc;

    /** MPI Info object. */
    MPI_Info info;

    /** Pointer to the next iosystem_desc_t in the list. */
    struct iosystem_desc_t *next;
} iosystem_desc_t;

/**
 * multi buffer.
 */
typedef struct wmulti_buffer
{
    int ioid;
    int validvars;
    int arraylen;
    int *vid;
    int *frame;
    void *fillvalue;
    void *data;
    struct wmulti_buffer *next;
} wmulti_buffer;

/**
 * File descriptor structure.
 *
 * This structure holds information associated with each open file
 */
typedef struct file_desc_t
{
    /** The IO system ID used to open this file. */
    iosystem_desc_t *iosystem;

    /** The buffersize does not seem to be used anywhere. */
    /* PIO_Offset buffsize;*/

    /** The ncid returned for this file by the underlying library
     * (netcdf or pnetcdf). */
    int fh;

    /** The PIO_TYPE value that was used to open this file. */
    int iotype;

    /** List of variables in this file. */
    struct var_desc_t varlist[PIO_MAX_VARS];

    /** ??? */
    int mode;

    /** ??? */
    struct wmulti_buffer buffer;

    /** Pointer to the next file_desc_t in the list of open files. */
    struct file_desc_t *next;

    /** True if this task should participate in IO (only true for one
     * task with netcdf serial files. */
    int do_io;
} file_desc_t;

/**
 * These are the supported methods of reading/writing netCDF
 * files. (Not all methods can be used with all netCDF files.)
 */
enum PIO_IOTYPE
{
    /** Parallel Netcdf  (parallel) */
    PIO_IOTYPE_PNETCDF = 1,

    /** Netcdf3 Classic format (serial) */
    PIO_IOTYPE_NETCDF = 2,

    /**  NetCDF4 (HDF5) compressed format (serial) */
    PIO_IOTYPE_NETCDF4C = 3,

    /** NetCDF4 (HDF5) parallel */
    PIO_IOTYPE_NETCDF4P = 4 
};

/**
 * These are the supported output data rearrangement methods.
 */
enum PIO_REARRANGERS
{
    /** Box rearranger. */
    PIO_REARR_BOX = 1,

    /** Subset rearranger. */
    PIO_REARR_SUBSET = 2
};

/**
 * These are the supported error handlers.
 */
enum PIO_ERROR_HANDLERS
{
    /** Errors cause abort. */
    PIO_INTERNAL_ERROR = (-51),

    /** Error codes are broadcast to all tasks. */
    PIO_BCAST_ERROR = (-52),

    /** Errors are returned to caller with no internal action. */
    PIO_RETURN_ERROR = (-53)  
};

/** Define the netCDF-based error codes. */
#if defined( _PNETCDF) || defined(_NETCDF)
#define PIO_GLOBAL NC_GLOBAL
#define PIO_UNLIMITED NC_UNLIMITED
#define PIO_DOUBLE NC_DOUBLE
#define PIO_REAL   NC_FLOAT
#define PIO_FLOAT  NC_FLOAT
#define PIO_INT    NC_INT
#define PIO_CHAR   NC_CHAR
#define PIO_NOERR  NC_NOERR
#define PIO_WRITE  NC_WRITE
#define PIO_NOWRITE  NC_NOWRITE
#define PIO_CLOBBER NC_CLOBBER	
#define PIO_NOCLOBBER NC_NOCLOBBER	
#define PIO_NOFILL NC_NOFILL
#define PIO_MAX_NAME NC_MAX_NAME
#define PIO_MAX_VAR_DIMS NC_MAX_VAR_DIMS
#define PIO_64BIT_OFFSET NC_64BIT_OFFSET
// NC_64BIT_DATA  This is a problem - need to define directly instead of using include file
#define PIO_64BIT_DATA 0x0010  
#define PIO_EBADID NC_EBADID
#define PIO_ENFILE NC_ENFILE
#define PIO_EEXIST NC_EEXIST
#define PIO_EINVAL NC_EINVAL
#define PIO_EPERM NC_EPERM
#define PIO_ENOTINDEFINE NC_ENOTINDEFINE
#define PIO_EINDEFINE NC_EINDEFINE
#define PIO_EINVALCOORDS NC_EINVALCOORDS
#define PIO_EMAXDIMS NC_EMAXDIMS
#define PIO_ENAMEINUSE NC_ENAMEINUSE
#define PIO_ENOTATT NC_ENOTATT
#define PIO_EMAXATTS NC_EMAXATTS
#define PIO_EBADTYPE NC_EBADTYPE
#define PIO_EBADDIM NC_EBADDIM
#define PIO_EUNLIMPOS NC_EUNLIMPOS
#define PIO_EMAXVARS NC_EMAXVARS
#define PIO_ENOTVAR NC_ENOTVAR
#define PIO_EGLOBAL NC_EGLOBAL
#define PIO_ENOTNC NC_ENOTNC
#define PIO_ESTS NC_ESTS
#define PIO_EMAXNAME NC_EMAXNAME
#define PIO_EUNLIMIT NC_EUNLIMIT
#define PIO_ENORECVARS NC_ENORECVARS
#define PIO_ECHAR NC_ECHAR
#define PIO_EEDGE NC_EEDGE
#define PIO_ESTRIDE NC_ESTRIDE
#define PIO_EBADNAME NC_EBADNAME
#define PIO_ERANGE NC_ERANGE
#define PIO_ENOMEM NC_ENOMEM
#define PIO_EVARSIZE NC_EVARSIZE
#define PIO_EDIMSIZE NC_EDIMSIZE
#define PIO_ETRUNC NC_ETRUNC
#define PIO_EAXISTYPE NC_EAXISTYPE
#define PIO_EDAP NC_EDAP
#define PIO_ECURL NC_ECURL
#define PIO_EIO NC_EIO
#define PIO_ENODATA NC_ENODATA
#define PIO_EDAPSVC NC_EDAPSVC
#define PIO_EDAS NC_EDAS
#define PIO_EDDS NC_EDDS
#define PIO_EDATADDS NC_EDATADDS
#define PIO_EDAPURL NC_EDAPURL
#define PIO_EDAPCONSTRAINT NC_EDAPCONSTRAINT
#define PIO_ETRANSLATION NC_ETRANSLATION
#define PIO_EHDFERR NC_EHDFERR
#define PIO_ECANTREAD NC_ECANTREAD
#define PIO_ECANTWRITE NC_ECANTWRITE
#define PIO_ECANTCREATE NC_ECANTCREATE
#define PIO_EFILEMETA NC_EFILEMETA
#define PIO_EDIMMETA NC_EDIMMETA
#define PIO_EATTMETA NC_EATTMETA
#define PIO_EVARMETA NC_EVARMETA
#define PIO_ENOCOMPOUND NC_ENOCOMPOUND
#define PIO_EATTEXISTS NC_EATTEXISTS
#define PIO_ENOTNC4 NC_ENOTNC4
#define PIO_ESTRICTNC3 NC_ESTRICTNC3
#define PIO_ENOTNC3 NC_ENOTNC3
#define PIO_ENOPAR NC_ENOPAR
#define PIO_EPARINIT NC_EPARINIT
#define PIO_EBADGRPID NC_EBADGRPID
#define PIO_EBADTYPID NC_EBADTYPID
#define PIO_ETYPDEFINED NC_ETYPDEFINED
#define PIO_EBADFIELD NC_EBADFIELD
#define PIO_EBADCLASS NC_EBADCLASS
#define PIO_EMAPTYPE NC_EMAPTYPE
#define PIO_ELATEFILL NC_ELATEFILL
#define PIO_ELATEDEF NC_ELATEDEF
#define PIO_EDIMSCALE NC_EDIMSCALE
#define PIO_ENOGRP NC_ENOGRP
#define PIO_ESTORAGE NC_ESTORAGE
#define PIO_EBADCHUNK NC_EBADCHUNK
#define PIO_ENOTBUILT NC_ENOTBUILT
#define PIO_EDISKLESS NC_EDISKLESS
#define PIO_FILL_DOUBLE NC_FILL_DOUBLE
#define PIO_FILL_FLOAT NC_FILL_FLOAT
#define PIO_FILL_INT NC_FILL_INT
#define PIO_FILL_CHAR NC_FILL_CHAR
#endif /*  defined( _PNETCDF) || defined(_NETCDF) */

/** Define the extra error codes for the parallel-netcdf library. */
#ifdef _PNETCDF
#define PIO_EINDEP  NC_EINDEP
#else  /* _PNETCDF */
#define PIO_EINDEP  (-203) 
#endif /* _PNETCDF */

/** Define error codes for PIO. */
#define PIO_EBADIOTYPE  -255

/** ??? */
#define PIO_REQ_NULL (NC_REQ_NULL-1)

#if defined(__cplusplus)
extern "C" {
#endif
    int PIOc_strerror(int pioerr, char *errstr);
    int PIOc_freedecomp(int iosysid, int ioid);
    int PIOc_inq_att (int ncid, int varid, const char *name, nc_type *xtypep, PIO_Offset *lenp); 
    int PIOc_inq_format (int ncid, int *formatp); 
    int PIOc_inq_varid (int ncid, const char *name, int *varidp); 
    int PIOc_inq_varnatts (int ncid, int varid, int *nattsp); 
    int PIOc_def_var (int ncid, const char *name, nc_type xtype,  int ndims, const int *dimidsp, int *varidp);
    int PIOc_def_var_deflate(int ncid, int varid, int shuffle, int deflate,
			     int deflate_level);
    int PIOc_inq_var_deflate(int ncid, int varid, int *shufflep, int *deflatep,
			     int *deflate_levelp);
    int PIOc_inq_var_szip(int ncid, int varid, int *options_maskp, int *pixels_per_blockp);
    int PIOc_def_var_chunking(int ncid, int varid, int storage, const PIO_Offset *chunksizesp);
    int PIOc_inq_var_chunking(int ncid, int varid, int *storagep, PIO_Offset *chunksizesp);
    int PIOc_def_var_fill(int ncid, int varid, int no_fill, const void *fill_value);
    int PIOc_inq_var_fill(int ncid, int varid, int *no_fill, void *fill_valuep);
    int PIOc_def_var_endian(int ncid, int varid, int endian);
    int PIOc_inq_var_endian(int ncid, int varid, int *endianp);
    int PIOc_set_chunk_cache(int iosysid, int iotype, PIO_Offset size, PIO_Offset nelems, float preemption);
    int PIOc_get_chunk_cache(int iosysid, int iotype, PIO_Offset *sizep, PIO_Offset *nelemsp, float *preemptionp);
    int PIOc_set_var_chunk_cache(int ncid, int varid, PIO_Offset size, PIO_Offset nelems,
				 float preemption);
    int PIOc_get_var_chunk_cache(int ncid, int varid, PIO_Offset *sizep, PIO_Offset *nelemsp,
				 float *preemptionp);
    int PIOc_inq_var (int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp, int *dimidsp, int *nattsp); 
    int PIOc_inq_varname (int ncid, int varid, char *name); 
    int PIOc_put_att_double (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const double *op); 
    int PIOc_put_att_int (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const int *op); 
    int PIOc_rename_att (int ncid, int varid, const char *name, const char *newname); 
    int PIOc_del_att (int ncid, int varid, const char *name); 
    int PIOc_inq_natts (int ncid, int *ngattsp); 
    int PIOc_inq (int ncid, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp);  
    int PIOc_get_att_text (int ncid, int varid, const char *name, char *ip); 
    int PIOc_get_att_short (int ncid, int varid, const char *name, short *ip); 
    int PIOc_put_att_long (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const long *op); 
    int PIOc_redef (int ncid); 
    int PIOc_set_fill (int ncid, int fillmode, int *old_modep); 
    int PIOc_enddef (int ncid); 
    int PIOc_rename_var (int ncid, int varid, const char *name); 
    int PIOc_put_att_short (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const short *op); 
    int PIOc_put_att_text (int ncid, int varid, const char *name, PIO_Offset len, const char *op); 
    int PIOc_inq_attname (int ncid, int varid, int attnum, char *name); 
    int PIOc_get_att_ulonglong (int ncid, int varid, const char *name, unsigned long long *ip); 
    int PIOc_get_att_ushort (int ncid, int varid, const char *name, unsigned short *ip); 
    int PIOc_put_att_ulonglong (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned long long *op); 
    int PIOc_inq_dimlen (int ncid, int dimid, PIO_Offset *lenp); 
    int PIOc_get_att_uint (int ncid, int varid, const char *name, unsigned int *ip); 
    int PIOc_get_att_longlong (int ncid, int varid, const char *name, long long *ip); 
    int PIOc_put_att_schar (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const signed char *op); 
    int PIOc_put_att_float (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const float *op); 
    int PIOc_inq_nvars (int ncid, int *nvarsp); 
    int PIOc_rename_dim (int ncid, int dimid, const char *name); 
    int PIOc_inq_varndims (int ncid, int varid, int *ndimsp); 
    int PIOc_get_att_long (int ncid, int varid, const char *name, long *ip); 
    int PIOc_inq_dim (int ncid, int dimid, char *name, PIO_Offset *lenp); 
    int PIOc_inq_dimid (int ncid, const char *name, int *idp); 
    int PIOc_inq_unlimdim (int ncid, int *unlimdimidp); 
    int PIOc_inq_vardimid (int ncid, int varid, int *dimidsp); 
    int PIOc_inq_attlen (int ncid, int varid, const char *name, PIO_Offset *lenp); 
    int PIOc_inq_dimname (int ncid, int dimid, char *name); 
    int PIOc_put_att_ushort (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned short *op); 
    int PIOc_get_att_float (int ncid, int varid, const char *name, float *ip); 
    int PIOc_sync (int ncid); 
    int PIOc_put_att_longlong (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const long long *op); 
    int PIOc_put_att_uint (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned int *op); 
    int PIOc_get_att_schar (int ncid, int varid, const char *name, signed char *ip); 
    int PIOc_inq_attid (int ncid, int varid, const char *name, int *idp); 
    int PIOc_def_dim (int ncid, const char *name, PIO_Offset len, int *idp); 
    int PIOc_inq_ndims (int ncid, int *ndimsp); 
    int PIOc_inq_vartype (int ncid, int varid, nc_type *xtypep); 
    int PIOc_get_att_int (int ncid, int varid, const char *name, int *ip); 
    int PIOc_get_att_double (int ncid, int varid, const char *name, double *ip); 
    int PIOc_inq_atttype (int ncid, int varid, const char *name, nc_type *xtypep); 
    int PIOc_put_att_uchar (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned char *op); 
    int PIOc_get_att_uchar (int ncid, int varid, const char *name, unsigned char *ip);              
    int PIOc_InitDecomp(const int iosysid, const int basetype,const int ndims, const int dims[], 
			const int maplen, const PIO_Offset *compmap, int *ioidp, const int *rearr, 
			const PIO_Offset *iostart,const PIO_Offset *iocount);
    int PIOc_Init_Intracomm(const MPI_Comm comp_comm, 
			    const int num_iotasks, const int stride, 
			    const int base, const int rearr, int *iosysidp);
    int PIOc_Init_Intercomm(int component_count, MPI_Comm peer_comm, MPI_Comm *comp_comms,
			    MPI_Comm io_comm, int *iosysidp);
    int PIOc_closefile(int ncid);
    int PIOc_createfile(const int iosysid, int *ncidp,  int *iotype,
			const char *fname, const int mode);
    int PIOc_openfile(const int iosysid, int *ncidp, int *iotype,
		      const char *fname, const int mode);
    int PIOc_write_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen,
			  void *array, void *fillvalue);
    int PIOc_write_darray_multi(const int ncid, const int vid[], const int ioid, const int nvars, const PIO_Offset arraylen,
				void *array, const int frame[], void *fillvalue[], bool flushtodisk);

    int PIOc_get_att_ubyte (int ncid, int varid, const char *name, unsigned char *ip);
    int PIOc_put_att_ubyte (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned char *op) ;
    int PIOc_set_blocksize(const int newblocksize);
    int PIOc_readmap(const char file[], int *ndims, int *gdims[], PIO_Offset *fmaplen, PIO_Offset *map[], const MPI_Comm comm);
    int PIOc_readmap_from_f90(const char file[],int *ndims, int *gdims[], PIO_Offset *maplen, PIO_Offset *map[], const int f90_comm);
    int PIOc_writemap(const char file[], const int ndims, const int gdims[], PIO_Offset maplen, PIO_Offset map[], const MPI_Comm comm);
    int PIOc_writemap_from_f90(const char file[], const int ndims, const int gdims[], const PIO_Offset maplen, const PIO_Offset map[], const int f90_comm);
    int PIOc_deletefile(const int iosysid, const char filename[]); 
    int PIOc_File_is_Open(int ncid);
    int PIOc_Set_File_Error_Handling(int ncid, int method);
    int PIOc_advanceframe(int ncid, int varid);
    int PIOc_setframe(const int ncid, const int varid,const int frame);
    int PIOc_get_numiotasks(int iosysid, int *numiotasks);
    int PIOc_get_iorank(int iosysid, int *iorank);
    int PIOc_get_local_array_size(int ioid);
    int PIOc_Set_IOSystem_Error_Handling(int iosysid, int method);
    int PIOc_set_hint(const int iosysid, char hint[], const char hintval[]);
    int PIOc_Init_Intracomm(const MPI_Comm comp_comm, 
			    const int num_iotasks, const int stride, 
			    const int base,const int rearr, int *iosysidp);
    int PIOc_finalize(const int iosysid);
    int PIOc_iam_iotask(const int iosysid, bool *ioproc);
    int PIOc_iotask_rank(const int iosysid, int *iorank);
    int PIOc_iosystem_is_active(const int iosysid, bool *active);
    int PIOc_put_vars_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned char *op) ;
    int PIOc_get_var1_schar (int ncid, int varid, const PIO_Offset index[], signed char *buf) ;
    int PIOc_put_vars_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned short *op) ;
    int pio_read_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid, void *IOBUF);
    int PIOc_put_vars_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned long long *op) ;
    int PIOc_get_vars_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned long long *buf) ;
    int PIOc_put_varm (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype)  ;
    int PIOc_read_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen, void *array);
    int PIOc_put_vars_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const unsigned int *op) ;
    int PIOc_get_varm_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], signed char *buf)  ;
    int PIOc_put_varm_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned char *op) ;
    int PIOc_put_var_ushort (int ncid, int varid, const unsigned short *op) ;
    int PIOc_get_vars_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], short *buf)  ;
    int PIOc_put_var1_longlong (int ncid, int varid, const PIO_Offset index[], const long long *op)  ;
    int PIOc_get_var_double (int ncid, int varid, double *buf) ;
    int PIOc_put_vara_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned char *op) ;
    int PIOc_put_varm_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const short *op)  ;
    int PIOc_get_vara_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], double *buf) ;
    int PIOc_put_var1_long (int ncid, int varid, const PIO_Offset index[], const long *ip) ;
    int PIOc_get_var_int (int ncid, int varid, int *buf) ;
    int PIOc_put_vars_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const long *op) ;
    int PIOc_put_var_short (int ncid, int varid, const short *op) ;
    int PIOc_get_vara_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], char *buf) ;
    int PIOc_put_vara_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const int *op)  ;

    int PIOc_put_var1_ushort (int ncid, int varid, const PIO_Offset index[], const unsigned short *op); 
    int PIOc_put_vara_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const char *op);  
    int PIOc_put_varm_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const char *op);  
    int PIOc_put_varm_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned short *op); 
    int PIOc_put_var_ulonglong (int ncid, int varid, const unsigned long long *op); 
    int PIOc_put_var_int (int ncid, int varid, const int *op); 
    int PIOc_put_var_longlong (int ncid, int varid, const long long *op); 
    int PIOc_put_var_schar (int ncid, int varid, const signed char *op); 
    int PIOc_put_var_uint (int ncid, int varid, const unsigned int *op); 
    int PIOc_put_var (int ncid, int varid, const void *buf, PIO_Offset bufcount, MPI_Datatype buftype); 
    int PIOc_put_vara_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned short *op); 
    int PIOc_put_vars_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const short *op);  
    int PIOc_put_vara_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned int *op); 
    int PIOc_put_vara_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const signed char *op);  
    int PIOc_put_varm_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned long long *op); 
    int PIOc_put_var1_uchar (int ncid, int varid, const PIO_Offset index[], const unsigned char *op); 
    int PIOc_put_varm_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const int *op);  
    int PIOc_put_vars_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const signed char *op);  
    int PIOc_put_var1 (int ncid, int varid, const PIO_Offset index[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype);  
    int PIOc_put_vara_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const float *op);  
    int PIOc_put_var1_float (int ncid, int varid, const PIO_Offset index[], const float *op);  
    int PIOc_put_varm_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const float *op);  
    int PIOc_put_var1_text (int ncid, int varid, const PIO_Offset index[], const char *op);  
    int PIOc_put_vars_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const char *op);  
    int PIOc_put_varm_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long *op); 
    int PIOc_put_vars_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const double *op);  
    int PIOc_put_vara_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const long long *op);  
    int PIOc_put_var_double (int ncid, int varid, const double *op); 
    int PIOc_put_var_float (int ncid, int varid, const float *op); 
    int PIOc_put_var1_ulonglong (int ncid, int varid, const PIO_Offset index[], const unsigned long long *op); 
    int PIOc_put_varm_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned int *op); 
    int PIOc_put_var1_uint (int ncid, int varid, const PIO_Offset index[], const unsigned int *op); 
    int PIOc_put_var1_int (int ncid, int varid, const PIO_Offset index[], const int *op);  
    int PIOc_put_vars_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const float *op);  
    int PIOc_put_vara_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const short *op);  
    int PIOc_put_var1_schar (int ncid, int varid, const PIO_Offset index[], const signed char *op);  
    int PIOc_put_vara_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const unsigned long long *op); 
    int PIOc_put_varm_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const double *op);  
    int PIOc_put_vara (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype);  
    int PIOc_put_vara_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const long *op); 
    int PIOc_put_var1_double (int ncid, int varid, const PIO_Offset index[], const double *op);  
    int PIOc_put_varm_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const signed char *op);  
    int PIOc_put_var_text (int ncid, int varid, const char *op); 
    int PIOc_put_vars_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const int *op);  
    int PIOc_put_var1_short (int ncid, int varid, const PIO_Offset index[], const short *op);  
    int PIOc_put_vars_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const long long *op);  
    int PIOc_put_vara_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const double *op);  
    int PIOc_put_vars (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype);  
    int PIOc_put_vars_tc(int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], nc_type xtype, const void *buf);  
    int PIOc_put_var_uchar (int ncid, int varid, const unsigned char *op); 
    int PIOc_put_var_long (int ncid, int varid, const long *op); 
    int PIOc_put_varm_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long long *op);
    int PIOc_get_vara_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], int *buf);
    int PIOc_get_var1_float (int ncid, int varid, const PIO_Offset index[], float *buf); 
    int PIOc_get_var1_short (int ncid, int varid, const PIO_Offset index[], short *buf); 
    int PIOc_get_vars_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], int *buf);                      
    int PIOc_get_var_text (int ncid, int varid, char *buf); 
    int PIOc_get_varm_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], double *buf); 
    int PIOc_get_vars_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], signed char *buf);  
    int PIOc_get_vara_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned short *buf); 
    int PIOc_get_var1_ushort (int ncid, int varid, const PIO_Offset index[], unsigned short *buf); 
    int PIOc_get_var_float (int ncid, int varid, float *buf); 
    int PIOc_get_vars_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned char *buf); 
    int PIOc_get_var (int ncid, int varid, void *buf, PIO_Offset bufcount, MPI_Datatype buftype); 
    int PIOc_get_var1_longlong (int ncid, int varid, const PIO_Offset index[], long long *buf); 
    int PIOc_get_vars_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned short *buf); 
    int PIOc_get_var_long (int ncid, int varid, long *buf); 
    int PIOc_get_var1_double (int ncid, int varid, const PIO_Offset index[], double *buf); 
    int PIOc_get_vara_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned int *buf); 
    int PIOc_get_vars_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], long long *buf); 
    int PIOc_get_var_longlong (int ncid, int varid, long long *buf); 
    int PIOc_get_vara_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], short *buf);  
    int PIOc_get_vara_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], long *buf); 
    int PIOc_get_var1_int (int ncid, int varid, const PIO_Offset index[], int *buf); 
    int PIOc_get_var1_ulonglong (int ncid, int varid, const PIO_Offset index[], unsigned long long *buf); 
    int PIOc_get_var_uchar (int ncid, int varid, unsigned char *buf); 
    int PIOc_get_vara_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned char *buf); 
    int PIOc_get_vars_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], float *buf);  
    int PIOc_get_vars_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], long *buf); 
    int PIOc_get_var1 (int ncid, int varid, const PIO_Offset index[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype); 
    int PIOc_get_var_uint (int ncid, int varid, unsigned int *buf); 
    int PIOc_get_vara (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype); 
    int PIOc_get_vara_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], signed char *buf); 
    int PIOc_get_var1_uint (int ncid, int varid, const PIO_Offset index[], unsigned int *buf); 
    int PIOc_get_vars_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned int *buf); 
    int PIOc_get_vara_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], float *buf);  
    int PIOc_get_varm_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], char *buf);  
    int PIOc_get_var1_text (int ncid, int varid, const PIO_Offset index[], char *buf); 
    int PIOc_get_varm_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], int *buf);  
    int PIOc_get_varm_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned int *buf);  
    int PIOc_get_varm (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype);  
    int PIOc_get_vars_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], double *buf); 
    int PIOc_get_vara_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], long long *buf); 
    int PIOc_get_var_ulonglong (int ncid, int varid, unsigned long long *buf); 
    int PIOc_get_vara_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned long long *buf); 
    int PIOc_get_var_short (int ncid, int varid, short *buf); 
    int PIOc_get_varm_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], float *buf);  
    int PIOc_get_var1_long (int ncid, int varid, const PIO_Offset index[], long *buf); 
    int PIOc_get_varm_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], long *buf); 
    int PIOc_get_varm_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned short *buf);  
    int PIOc_get_varm_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], long long *buf); 
    int PIOc_get_vars_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], char *buf);  
    int PIOc_get_var1_uchar (int ncid, int varid, const PIO_Offset index[], unsigned char *buf); 
    int PIOc_get_vars (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype);  
    int PIOc_get_vars_tc(int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], nc_type xtype, void *buf);  
    int PIOc_get_varm_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], short *buf);  
    int PIOc_get_varm_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned long long *buf);  
    int PIOc_get_var_schar (int ncid, int varid, signed char *buf); 
    int PIOc_iotype_available(const int iotype);
    int PIOc_set_log_level(int level);
    int PIOc_put_att(int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const void *op);
    int PIOc_get_att(int ncid, int varid, const char *name, void *ip);
    int PIOc_inq_type(int ncid, nc_type xtype, char *name, PIO_Offset *sizep);    
#if defined(__cplusplus)
}
#endif

#endif  // _PIO_H_

