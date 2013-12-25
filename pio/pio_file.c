#include <pio.h>
#include <pio_internal.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

int PIOc_OpenFile(const int iosysid, int *ncidp, const int iotype,
		  char *fname, const int mode, const _Bool *checkmpi)
{
  int ierr;
  int msg;
  int mpierr;
  int amode;
  size_t len;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  msg = PIO_MSG_OPEN_FILE;
  amode = mode;
  ios = pio_get_iosystem_from_id(iosysid);
  file = (file_desc_t *) malloc(sizeof(file_desc_t));
  file->next = NULL;
  file->iotype = iotype;
  file->iosystem = ios;

#ifdef _NETCDF
  if(ios->num_iotasks==1 && file->iotype==PIO_IOTYPE_PNETCDF) {
    fprintf(stderr,"WARNING: only 1 iotask - changing iotype to netcdf\n");
    file->iotype = PIO_IOTYPE_NETCDF;
  }
#ifndef _NETCDF4
  if(iotype==PIO_IOTYPE_NETCDF4P || iotype==PIO_IOTYPE_NETCDF4C){
    if(ios->io_rank==0)
      printf(stderr,"WARNING: PIO was not built with NETCDF 4 support, changing iotype to netcdf\n");
    file->iotype = PIO_IOTYPE_NETCDF;
  }
#endif
#endif  

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    len = strlen(fname);
    mpierr = MPI_Bcast((void *) fname,len, MPI_CHAR, ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&(file->iotype), 1, MPI_INT,  ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&amode, 1, MPI_INT,  ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&checkmpi, 1, MPI_INT,  ios->compmaster, ios->intercomm);
  }
  

  if(ios->ioproc){
    switch(file->iotype){
    case PIO_IOYTPE_DIRECT_PBINARY:
    case PIO_IOTYPE_PBINARY:
      //      ierr = pio_open_mpiio(file, fname, checkmpi);
      break;
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
#ifdef _MPISERIAL      
      ierr = nc_open(fname, amode, &(file->fh));
#else
      ierr = nc_open_par(fname, amode, ios->io_comm,ios->info, &(file->fh));
#endif
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_open(fname, amode, &(file->fh));
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_open(ios->io_comm, fname, amode, ios->info, &(file->fh));
      break;
#endif
    default:
      ierr = iotype_error(iotype,__FILE__,__LINE__);
    }
  }

  pio_add_to_file_list(file);
  

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);


  return ierr;
}



int PIOc_CreateFile(const int iosysid, int *ncidp, const int iotype,
		 const char *fname, const int mode, const _Bool checkmpi)
{
  int ierr;
  int msg;
  int mpierr;
  int amode;
  size_t len;
  iosystem_desc_t *ios;
  file_desc_t *file;


  ierr = PIO_NOERR;

  ios = pio_get_iosystem_from_id(iosysid);
  file = (file_desc_t *) malloc(sizeof(file_desc_t));
  file->next = NULL;
  file->iotype = iotype;
  file->iosystem = ios;

  msg = PIO_MSG_CREATE_FILE;
  amode = mode;
#ifdef _NETCDF  
  if(ios->num_iotasks==1 && iotype==PIO_IOTYPE_PNETCDF) {
    fprintf(stderr,"WARNING: only 1 iotask - changing iotype to netcdf\n");
    file->iotype = PIO_IOTYPE_NETCDF;
  }
#ifndef _NETCDF4
  if(iotype==PIO_IOTYPE_NETCDF4P || iotype==PIO_IOTYPE_NETCDF4C){
    if(ios->io_rank==0)
      fprintf(stderr,"WARNING: PIO was not built with NETCDF 4 support, changing iotype to netcdf\n");
    file->iotype = PIO_IOTYPE_NETCDF;
  }
#endif
#endif

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = MPI_Send( &msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    len = strlen(fname);
    mpierr = MPI_Bcast((void *) fname,len, MPI_CHAR, ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&(file->iotype), 1, MPI_INT,  ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&amode, 1, MPI_INT,  ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&checkmpi, 1, MPI_INT,  ios->compmaster, ios->intercomm);
  }
  

  if(ios->ioproc){
    switch(iotype){
    case PIO_IOYTPE_DIRECT_PBINARY:
    case PIO_IOTYPE_PBINARY:
      //      ierr = pio_create_mpiio(file, fname, checkmpi);
      break;
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
#ifdef _MPISERIAL      
      ierr = nc_create(fname, amode, &(file->fh));
#else
      ierr = nc_create_par(fname, amode, ios->io_comm,ios->info, &(file->fh));
#endif
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_create(fname, amode, &(file->fh));
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_create(ios->io_comm, fname, amode, ios->info, &(file->fh));
      break;
#endif
    default:
      ierr = iotype_error(iotype,__FILE__,__LINE__);
    }
  }
  pio_add_to_file_list(file);

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
