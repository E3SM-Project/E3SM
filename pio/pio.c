#include <pio.h>
#include <pio_internal.h>

int PIOc_File_is_Open(int ncid)
{
  file_desc_t *file;
  file = pio_get_file_from_id(ncid);
  if(file==NULL)
    return 0;
  else
    return 1;
}

int PIOc_Set_File_Error_Handling(int ncid, int method)
{
  file_desc_t *file;
  int oldmethod;
  file = pio_get_file_from_id(ncid);
  oldmethod = file->iosystem->error_handler;
  file->iosystem->error_handler = method;
  return(oldmethod);
} 

 int PIOc_Set_IOSystem_Error_Handling(int iosysid, int method)
{
  iosystem_desc_t *ios;
  int oldmethod;
  ios = pio_get_iosystem_from_id(iosysid);
  oldmethod = ios->error_handler;
  ios->error_handler = method;
  return(oldmethod);
}  

int PIOc_InitDecomp(const int iosysid, const int basetype,const int *dims, 
		    const int lenblocks, const PIO_Offset *compdof, const PIO_Offset *iodofr, 
		    const PIO_Offset *iodofw,int *iodescp)
{
  return PIO_NOERR;
}

int PIOc_Init_Intracomm(const MPI_Comm comp_comm, 
			const int num_iotasks, const int stride, 
			const int base, int *iosysidp)
{
  iosystem_desc_t *iosys;
  int ierr;
  MPI_Group compgroup, iogroup;
  iosys = (iosystem_desc_t *) malloc(sizeof(iosystem_desc_t));

  iosys->union_comm = comp_comm;
  iosys->comp_comm = comp_comm;
  iosys->my_comm = comp_comm;
  iosys->io_comm = MPI_COMM_NULL;
  iosys->intercomm = MPI_COMM_NULL;
  iosys->error_handler = PIO_INTERNAL_ERROR;
  iosys->async_interface= false;
  iosys->iomaster = false;
  iosys->compmaster = false;
  iosys->ioproc = false;
#ifndef _MPISERIAL
  iosys->info = MPI_INFO_NULL;
#endif

  ierr = MPI_Comm_rank(comp_comm, &(iosys->comp_rank));
  if(iosys->comp_rank==0)
    iosys->compmaster = true;
  ierr = MPI_Comm_size(comp_comm, &(iosys->num_comptasks));
  
  if((num_iotasks < 1) || ((num_iotasks*stride) > iosys->num_comptasks)){
    return PIO_EBADID;
  }
#ifndef _MPISERIAL
  ierr = MPI_Info_create(&(iosys->info));
#endif

  iosys->ioranks = (int *) calloc(sizeof(int), iosys->num_iotasks);
  for(int i=0;i< num_iotasks; i++){
    iosys->ioranks[i] = (base + i*stride) % iosys->num_comptasks;
    if(iosys->ioranks[i] == iosys->comp_rank)
      iosys->ioproc = true;
  }
  iosys->ioroot = iosys->ioranks[0];

  ierr = MPI_Comm_group(comp_comm, &compgroup);
			
  ierr = MPI_Group_incl(compgroup, num_iotasks, iosys->ioranks, &iogroup);

  ierr = MPI_Comm_create(comp_comm, iogroup, &(iosys->io_comm));
  if(iosys->ioproc)
    ierr = MPI_Comm_rank(iosys->io_comm, &(iosys->io_rank));
  else
    iosys->io_rank = -1;

  iosys->num_aiotasks = iosys->num_iotasks;
  
  pio_add_to_iosystem_list(iosys);

  return PIO_NOERR;
}
  
  
  

    
