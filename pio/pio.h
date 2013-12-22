#include <mpi.h>

struct iosystem_desc_t
{
  MPI_Comm union_comm;
  MPI_Comm io_comm;
  MPI_Comm comp_comm;
  MPI_Comm intercomm;
  MPI_Comm my_comm;
  
  int num_tasks;
  int num_iotasks;
  int num_aiotasks;
  int num_comptasks;

  int union_rank;
  int comp_rank;
  int io_rank;

  int iomaster;
  int compmaster;

  int ioroot;
  int comproot;

  int error_handler;

  int async_interface;
  int ioproc;
  
  MPI_Info info;
  
};

struct file_desc_t
{
  struct iosystem_desc_t *iosystem;
  long int buffsize;
  int request_cnt;
  int fh;
  int iotype;
  int file_open;
};

enum PIO_IOTYPE{
  PIO_IOTYPE_PBINARY,
  PIO_IOYTPE_DIRECT_PBINARY,
  PIO_IOTYPE_PNETCDF,
  PIO_IOTYPE_NETCDF,
  PIO_IOTYPE_NETCDF4C,
  PIO_IOTYPE_NETCDF4P
};

#ifdef _PNETCDF
 #define PIO_GLOBAL NC_GLOBAL
 #define PIO_UNLIMITED NC_UNLIMITED
 #define PIO_DOUBLE NC_DOUBLE
 #define PIO_REAL   NC_REAL
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
 #define PIO_64BIT_DATA NC_64BIT_DATA
#endif
#define PIO_Offset MPI_Offset
