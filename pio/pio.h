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
  


  MPI_Info info;
  
}

struct file_desc_t
{
  struct iosystem_desc_t *iosystem;
  long int buffsize;
  int request_cnt;
  int fh;
  int iotype;
  int file_open;
}
