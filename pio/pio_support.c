#include <mpi.h>
#include <pio.h>
#include <string.h>
#include <stdio.h>

#ifdef _NETCDF
#include <netcdf.h>
#endif
#ifdef _PNETCDF
#include <pnetcdf.h>
#endif

int check_netcdf(file_desc_t *file,const int status, const char *fname, const int line){
  iosystem_desc_t *ios;
  int ierr;

  ios = file->iosystem;
  ierr = PIO_NOERR;

  switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
  case PIO_IOTYPE_NETCDF4P:
  case PIO_IOTYPE_NETCDF4C:
#endif
  case PIO_IOTYPE_NETCDF:
    if(status != NC_NOERR){
      if(ios->iomaster){
	fprintf(stderr,"NETCDF ERROR: %s\n",nc_strerror(status));
      }
      if(ios->error_handler == PIO_INTERNAL_ERROR){
	// abort
      }else if(ios->error_handler==PIO_BCAST_ERROR){
	ierr = mpi_bcast(status, 1, MPI_INT, ios->iomaster, ios->my_comm);
      }
    }
    break;
#endif
#ifdef _PNETCDF
  case PIO_IOTYPE_PNETCDF:
    if(status != NC_NOERR){
      if(ios->iomaster){
	fprintf(stderr,"PNETCDF ERROR: %s\n",ncmpi_strerror(status));
      }
      if(ios->error_handler == PIO_INTERNAL_ERROR){
	// abort
      }else if(ios->error_handler==PIO_BCAST_ERROR){
	ierr = mpi_bcast(status, 1, MPI_INT, ios->iomaster, ios->my_comm);
	  }
    }
    break;
#endif
  default:
    ierr = iotype_error(file->iotype,__FILE__,__LINE__);
  } 
  return ierr;
}

int iotype_error(const int iotype, const char *fname, const int line){
    fprintf(stderr, "ERROR: iotype %d not defined in build %s %d\n", iotype, fname,line);
    return PIO_NOERR;
}
