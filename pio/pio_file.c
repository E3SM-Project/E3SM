#include <pio.h>
#include <pio_internal.h>
#include <string.h>
#include <stdio.h>

int PIOc_openfile(const int iosysid, int *ncidp, int *iotype,
		  char *fname, const int mode, _Bool checkmpi)
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
  if(ios==NULL){
    printf("bad iosysid %d\n",iosysid);
    return PIO_EBADID;
  }

  file = (file_desc_t *) malloc(sizeof(file_desc_t));
  file->next = NULL;
  file->iosystem = ios;
  
#ifdef _NETCDF
  if(ios->num_iotasks==1 && *iotype==PIO_IOTYPE_PNETCDF) {
    if(ios->io_rank==0)
      fprintf(stderr,"WARNING: only 1 iotask - changing iotype to netcdf\n");
    *iotype = PIO_IOTYPE_NETCDF;
  }
#ifndef _NETCDF4
  if(*iotype==PIO_IOTYPE_NETCDF4P || *iotype==PIO_IOTYPE_NETCDF4C){
    if(ios->io_rank==0)
      fprintf(stderr,"WARNING: PIO was not built with NETCDF 4 support, changing iotype to netcdf\n");
    *iotype = PIO_IOTYPE_NETCDF;
  }
#endif
#endif  
#ifndef _PNETCDF
  if(*iotype==PIO_IOTYPE_PNETCDF) {
    if(ios->io_rank==0)
      fprintf(stderr,"WARNING: pnetcdf not supported in build - changing iotype to netcdf\n");
    *iotype = PIO_IOTYPE_NETCDF;
  }
#endif    
  file->iotype = *iotype;

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
      amode = amode & PIO_64BIT_DATA;
      amode = amode  & PIO_64BIT_OFFSET;
      //printf("%d %d  \n",__LINE__,amode);
      amode = amode |  NC_MPIIO;

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
      // This should only be done with a file opened to append
      if(ierr == PIO_NOERR && (amode & PIO_WRITE)){
	if(ios->iomaster) printf("%d Setting IO buffer %d\n",__LINE__,PIO_BUFFER_SIZE_LIMIT);
	ierr = ncmpi_buffer_attach(file->fh, PIO_BUFFER_SIZE_LIMIT );
	for(int i=0;i<PIO_MAX_VARS;i++)
	  file->request[i]=MPI_REQUEST_NULL;
	file->nreq=0;
      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }


  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);
  if(ierr==PIO_NOERR){
    *ncidp = file->fh;
    mpierr = MPI_Bcast(ncidp, 1, MPI_INT, ios->ioroot, ios->my_comm);
    if((*ncidp) != file->fh) {
      *ncidp = -(*ncidp);
      file->fh = (*ncidp);
    }
    for(int i=0; i<PIO_MAX_VARS;i++){
      file->varlist[i].record = -1;
      file->varlist[i].ndims = -1;
      file->varlist[i].buffer = NULL;
    }

    pio_add_to_file_list(file);
  }

  return ierr;
}



int PIOc_createfile(const int iosysid, int *ncidp,  int *iotype,
		 const char *fname, const int mode)
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
  file->iosystem = ios;

  msg = PIO_MSG_CREATE_FILE;
  amode = mode;
#ifdef _NETCDF  
  if(ios->num_iotasks==1 && *iotype==PIO_IOTYPE_PNETCDF) {
    if(ios->io_rank==0)
      fprintf(stderr,"WARNING: only 1 iotask - changing iotype to netcdf\n");
    *iotype = PIO_IOTYPE_NETCDF;
  }
#ifndef _NETCDF4
  if(*iotype==PIO_IOTYPE_NETCDF4P || *iotype==PIO_IOTYPE_NETCDF4C){
    if(ios->io_rank==0)
      fprintf(stderr,"WARNING: PIO was not built with NETCDF 4 support, changing iotype to netcdf\n");
    *iotype = PIO_IOTYPE_NETCDF;
  }
#endif
#endif
#ifndef _PNETCDF
  if(*iotype==PIO_IOTYPE_PNETCDF) {
    if(ios->io_rank==0)
      fprintf(stderr,"WARNING: pnetcdf not supported in build - changing iotype to netcdf\n");
    *iotype = PIO_IOTYPE_NETCDF;
  }
#endif    
  file->iotype = *iotype;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = MPI_Send( &msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    len = strlen(fname);
    mpierr = MPI_Bcast((void *) fname,len, MPI_CHAR, ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&(file->iotype), 1, MPI_INT,  ios->compmaster, ios->intercomm);
    mpierr = MPI_Bcast(&amode, 1, MPI_INT,  ios->compmaster, ios->intercomm);
  }
  

  if(ios->ioproc){
    switch(file->iotype){
    case PIO_IOYTPE_DIRECT_PBINARY:
    case PIO_IOTYPE_PBINARY:
      //      ierr = pio_create_mpiio(file, fname);
      break;
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      //         The 64 bit options are not compatable with hdf5 format files
      //      printf("%d %d %d %d %d \n",__LINE__,amode,PIO_64BIT_DATA, PIO_64BIT_OFFSET, NC_MPIIO);
      amode = amode & PIO_64BIT_DATA;
      amode = amode  & PIO_64BIT_OFFSET;
      //printf("%d %d  \n",__LINE__,amode);
      amode = amode |  NC_MPIIO;
      printf("%d %d \n",__LINE__,amode);

      ierr = nc_create_par(fname, amode, ios->io_comm,ios->info  , &(file->fh));
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
      if(ierr == PIO_NOERR){
	if(ios->iomaster) printf("%d Setting IO buffer %d\n",__LINE__,PIO_BUFFER_SIZE_LIMIT);
	ierr = ncmpi_buffer_attach(file->fh, PIO_BUFFER_SIZE_LIMIT );
	for(int i=0;i<PIO_MAX_VARS;i++)
	  file->request[i]=MPI_REQUEST_NULL;
	file->nreq=0;

      }
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }
  
  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ierr == PIO_NOERR){
    *ncidp = file->fh;
    mpierr = MPI_Bcast(ncidp, 1, MPI_INT, ios->ioroot, ios->my_comm);
    if((*ncidp) != file->fh) {
      *ncidp = -(*ncidp);
      file->fh = (*ncidp);
    }
    for(int i=0; i<PIO_MAX_VARS;i++){
      file->varlist[i].record = -1;
      file->varlist[i].ndims = -1;
      file->varlist[i].buffer = NULL;
    }

    pio_add_to_file_list(file);
  }
  return ierr;
}

int PIOc_closefile(int ncid)
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_close(file->fh);
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_close(file->fh);
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      flush_output_buffer(file);
      ierr = ncmpi_buffer_detach(file->fh);
      ierr = ncmpi_close(file->fh);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  int iret =  pio_delete_file_from_list(ncid);


  return ierr;
}

int PIOc_deletefile(const int iosysid, const char fname[])
{
  int ierr;
  int msg;
  int mpierr;
  int chkerr;
  iosystem_desc_t *ios;

  ierr = PIO_NOERR;
  ios = pio_get_iosystem_from_id(iosysid);

  if(ios == NULL)
    return PIO_EBADID;

  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    //    mpierr = MPI_Bcast(iosysid,1, MPI_INT, ios->compmaster, ios->intercomm);
  }
  // The barriers are needed to assure that no task is trying to operate on the file while it is being deleted.
  if(ios->ioproc){
    MPI_Barrier(ios->io_comm);
#ifdef _NETCDF
    if(ios->io_rank==0)
      ierr = nc_delete(fname);
#else
#ifdef _PNETCDF
    ierr = ncmpi_delete(fname, ios->info);
#endif
#endif
    MPI_Barrier(ios->io_comm);
  }
  //   Special case - always broadcast the return from the  
  MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm);
  


  return ierr;
}
