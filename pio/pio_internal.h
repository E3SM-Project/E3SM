#ifndef __PIO_INTERNAL__
#define __PIO_INTERNAL__
#include <pio.h>

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#define MAX_GATHER_BLOCK_SIZE 32

int  pio_add_to_iodesc_list(io_desc_t *iodesc);
io_desc_t *pio_get_iodesc_from_id(int ioid);
int pio_delete_iodesc_from_list(int ioid);

file_desc_t *pio_get_file_from_id(int ncid);
int pio_delete_file_from_list(int ncid);
void pio_add_to_file_list(file_desc_t *file);

iosystem_desc_t *pio_get_iosystem_from_id(int iosysid);
int pio_add_to_iosystem_list(iosystem_desc_t *ios);

int check_netcdf(file_desc_t *file,const int status, const char *fname, const int line);
int iotype_error(const int iotype, const char *fname, const int line);
void piodie(const char *msg,const char *fname, const int line);
int CalcStartandCount(const int basetype, const int ndims, const int *gdims, const int num_io_procs,
		      const int myiorank, PIO_Offset *start, PIO_Offset *kount);
void CheckMPIReturn(const int ierr,const char file[],const int line);
int pio_fc_gather( void *sendbuf, const int sendcnt, const MPI_Datatype sendtype,
		   void *recvbuf, const int recvcnt, const MPI_Datatype recvtype, const int root, 
		   MPI_Comm comm, const int flow_cntl);

int pio_swapm(const int nprocs, const int mytask, void *sndbuf, const int sndlths[],
	      const int sdispls[], const MPI_Datatype stypes[], void *rcvbuf, const int rcvlths[], 
	      const int rdispls[], const MPI_Datatype rtypes[], MPI_Comm comm, const bool handshake, const bool isend, 
	      const int max_requests);
long long lgcd_array(int nain, long long*ain);
PIO_Offset GCDblocksize(const int arrlen, const PIO_Offset arr_in[]);

int box_rearrange_create(iosystem_desc_t *ios,const int maplen, const PIO_Offset compmap[], const int gsize[],
			 const int ndim, const int nioproc, io_desc_t *iodesc);
int box_rearrange_io2comp(iosystem_desc_t *ios, io_desc_t *iodesc, const int slen, void *sbuf,const int rlen,
			  void *rbuf, const int comm_option, const int fc_options);
int box_rearrange_comp2io(iosystem_desc_t *ios, io_desc_t *iodesc, const int slen, void *sbuf,const int rlen,
			  void *rbuf, const int comm_option, const int fc_options);

io_desc_t *malloc_iodesc(const int piotype, const int ndims);

#endif
