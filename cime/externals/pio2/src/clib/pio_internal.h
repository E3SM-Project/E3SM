#ifndef __PIO_INTERNAL__
#define __PIO_INTERNAL__
#include <pio.h>
// It seems that some versions of openmpi fail to define MPI_OFFSET
#ifdef OMPI_OFFSET_DATATYPE
#ifndef MPI_OFFSET
#define MPI_OFFSET OMPI_OFFSET_DATATYPE
#endif
#endif
#ifndef MPI_Offset
#define MPI_Offset long long
#endif
#include <bget.h>
#include <limits.h>
#include <math.h>
#ifdef TIMING
#include <gptl.h>
#endif


#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


#define MAX_GATHER_BLOCK_SIZE 0
#define PIO_REQUEST_ALLOC_CHUNK 16

#if defined(__cplusplus)
extern "C" {
#endif

extern PIO_Offset PIO_BUFFER_SIZE_LIMIT;
extern bool PIO_Save_Decomps;


/**
 ** @brief Used to sort map points in the subset rearranger
*/
typedef struct mapsort
{
  int rfrom;
  PIO_Offset soffset;
  PIO_Offset iomap;
} mapsort;

/**
 * @brief swapm defaults.
 *
*/
typedef struct pio_swapm_defaults
{
  int nreqs;
  bool handshake;
  bool isend;
} pio_swapm_defaults;


  void pio_get_env(void);
  int  pio_add_to_iodesc_list(io_desc_t *iodesc);
  io_desc_t *pio_get_iodesc_from_id(int ioid);
  int pio_delete_iodesc_from_list(int ioid);

  file_desc_t *pio_get_file_from_id(int ncid);
  int pio_delete_file_from_list(int ncid);
  void pio_add_to_file_list(file_desc_t *file);
  void pio_push_request(file_desc_t *file, int request);

  iosystem_desc_t *pio_get_iosystem_from_id(int iosysid);
  int pio_add_to_iosystem_list(iosystem_desc_t *ios);
  
  int check_netcdf(file_desc_t *file,const int status, const char *fname, const int line);
  int iotype_error(const int iotype, const char *fname, const int line);
  void piodie(const char *msg,const char *fname, const int line);
  void pioassert(bool exp, const char *msg,const char *fname, const int line);
  int CalcStartandCount(const int basetype, const int ndims, const int *gdims, const int num_io_procs,
			const int myiorank, PIO_Offset *start, PIO_Offset *kount);
  void CheckMPIReturn(const int ierr,const char file[],const int line);
  int pio_fc_gather( void *sendbuf, const int sendcnt, const MPI_Datatype sendtype,
		     void *recvbuf, const int recvcnt, const MPI_Datatype recvtype, const int root, 
		     MPI_Comm comm, const int flow_cntl);
  int pio_fc_gatherv( void *sendbuf, const int sendcnt, const MPI_Datatype sendtype,
		      void *recvbuf, const int recvcnts[], const int recvdispl[], const MPI_Datatype recvtype, const int root, 
		      MPI_Comm comm, const int flow_cntl);
  
  int pio_fc_gatherv( void *sendbuf, const int sendcnt, const MPI_Datatype sendtype,
		      void *recvbuf, const int recvcnts[], const int rdispls[], const MPI_Datatype recvtype, const int root, 
		     MPI_Comm comm, const int flow_cntl);
  
  int pio_swapm(void *sndbuf, int sndlths[], int sdispls[], MPI_Datatype stypes[], 
		void *rcvbuf, int rcvlths[], int rdispls[], MPI_Datatype rtypes[], 
		MPI_Comm comm, const bool handshake, bool isend, const int max_requests);
  
  long long lgcd_array(int nain, long long*ain);
  
  void PIO_Offset_size(MPI_Datatype *dtype, int *tsize);
  PIO_Offset GCDblocksize(const int arrlen, const PIO_Offset arr_in[]);
  
  int subset_rearrange_create(const iosystem_desc_t ios,const int maplen, PIO_Offset compmap[], const int gsize[],
			      const int ndim, io_desc_t *iodesc);
  

  int box_rearrange_create(const iosystem_desc_t ios,const int maplen, const PIO_Offset compmap[], const int gsize[],
			   const int ndim, io_desc_t *iodesc);
  

  int rearrange_io2comp(const iosystem_desc_t ios, io_desc_t *iodesc, void *sbuf,
			void *rbuf);
  int rearrange_comp2io(const iosystem_desc_t ios, io_desc_t *iodesc, void *sbuf,
			void *rbuf, const int nvars);
  int calcdisplace(const int bsize, const int numblocks,const PIO_Offset map[],int displace[]);
  io_desc_t *malloc_iodesc(const int piotype, const int ndims);
  void performance_tune_rearranger(iosystem_desc_t ios, io_desc_t *iodesc);
  
  int flush_output_buffer(file_desc_t *file, bool force, PIO_Offset addsize);
  void compute_maxIObuffersize(MPI_Comm io_comm, io_desc_t *iodesc);
  io_region *alloc_region(const int ndims);
  int pio_delete_iosystem_from_list(int piosysid);
  int gcd(int a, int b); 
  long long lgcd (long long a,long long b );
  int gcd_array(int nain, int *ain);
  void free_region_list(io_region *top);
  void gindex_to_coord(const int ndims, const PIO_Offset gindex, const PIO_Offset gstride[], PIO_Offset *gcoord);
  PIO_Offset coord_to_lindex(const int ndims, const PIO_Offset lcoord[], const PIO_Offset count[]);

  int ceil2(const int i);
  int pair(const int np, const int p, const int k);
  int define_iodesc_datatypes(const iosystem_desc_t ios, io_desc_t *iodesc);

  int create_mpi_datatypes(const MPI_Datatype basetype,const int msgcnt,const PIO_Offset dlen, 
			   const PIO_Offset mindex[],const int mcount[],int *mfrom, MPI_Datatype mtype[]);
  int compare_offsets(const void *a,const void *b) ;

  int subset_rearrange_create(const iosystem_desc_t ios, int maplen, PIO_Offset compmap[], 
			      const int gsize[], const int ndims, io_desc_t *iodesc);
  void print_trace (FILE *fp);
  void cn_buffer_report(iosystem_desc_t ios, bool collective);
  void compute_buffer_init(iosystem_desc_t ios);
  void free_cn_buffer_pool(iosystem_desc_t ios);
void flush_buffer(int ncid, wmulti_buffer *wmb, bool flushtodisk);
  void piomemerror(iosystem_desc_t ios, size_t req, char *fname, const int line);
  void compute_maxaggregate_bytes(const iosystem_desc_t ios, io_desc_t *iodesc);

#ifdef BGQ
  void identity(MPI_Comm comm, int *iotask);
  void determineiotasks(const MPI_Comm comm, int *numiotasks,int *base, int *stride, int *rearr, 
			bool *iamIOtask);

#endif

#if defined(__cplusplus)
}
#endif

enum PIO_MSG{
  PIO_MSG_OPEN_FILE,
  PIO_MSG_CREATE_FILE,
  PIO_MSG_INQ_ATT,
  PIO_MSG_INQ_FORMAT,
  PIO_MSG_INQ_VARID,
  PIO_MSG_INQ_VARNATTS,
  PIO_MSG_DEF_VAR,
  PIO_MSG_INQ_VAR,
  PIO_MSG_INQ_VARNAME,
  PIO_MSG_PUT_ATT_DOUBLE,
  PIO_MSG_PUT_ATT_INT,
  PIO_MSG_RENAME_ATT,
  PIO_MSG_DEL_ATT,
  PIO_MSG_INQ_NATTS,
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
  PIO_MSG_INQ_DIMLEN,
  PIO_MSG_GET_ATT_UINT,
  PIO_MSG_GET_ATT_LONGLONG,
  PIO_MSG_PUT_ATT_SCHAR,
  PIO_MSG_PUT_ATT_FLOAT,
  PIO_MSG_INQ_NVARS,
  PIO_MSG_RENAME_DIM,
  PIO_MSG_INQ_VARNDIMS,
  PIO_MSG_GET_ATT_LONG,
  PIO_MSG_INQ_DIM,
  PIO_MSG_INQ_DIMID,
  PIO_MSG_INQ_UNLIMDIM,
  PIO_MSG_INQ_VARDIMID,
  PIO_MSG_INQ_ATTLEN,
  PIO_MSG_INQ_DIMNAME,
  PIO_MSG_PUT_ATT_USHORT,
  PIO_MSG_GET_ATT_FLOAT,
  PIO_MSG_SYNC,
  PIO_MSG_PUT_ATT_LONGLONG,
  PIO_MSG_PUT_ATT_UINT,
  PIO_MSG_GET_ATT_SCHAR,
  PIO_MSG_INQ_ATTID,
  PIO_MSG_DEF_DIM,
  PIO_MSG_INQ_NDIMS,
  PIO_MSG_INQ_VARTYPE,
  PIO_MSG_GET_ATT_INT,
  PIO_MSG_GET_ATT_DOUBLE,
  PIO_MSG_INQ_ATTTYPE,
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
  PIO_MSG_GET_VAR_CHUNK_CACHE
};

#endif
