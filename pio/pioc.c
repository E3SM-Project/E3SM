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


int PIOc_InitDecomp(const int iosysid, const int basetype,const int ndims, const int dims[], 
		    const int maplen, const PIO_Offset *compmap, int *ioidp, PIO_Offset *iostart,PIO_Offset *iocount)
{
  iosystem_desc_t *ios;
  io_desc_t *iodesc;
  int mpierr;
  int ierr;
  int iosize;
  int lenblocks;
  int ndisp;

  for(int i=0;i<ndims;i++){
    if(dims[i]<=0){
      piodie("Invalid dims argument",__FILE__,__LINE__);
    }
  }
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  iodesc = (io_desc_t *) malloc(sizeof(io_desc_t));

  switch(basetype){
  case PIO_REAL:			
    iodesc->basetype=MPI_FLOAT;
    break;
  case PIO_DOUBLE:
    iodesc->basetype=MPI_DOUBLE;
    break;
  case PIO_CHAR:
    iodesc->basetype=MPI_CHAR;
    break;
  case PIO_INT:   
  defaut:
    iodesc->basetype = MPI_INTEGER;
    break;
  }    

  iodesc->ndims = ndims;
  iodesc->start = (PIO_Offset *) malloc(sizeof(PIO_Offset)*ndims);
  iodesc->count = (PIO_Offset *) malloc(sizeof(PIO_Offset)*ndims);
  for(int i=0;i<ndims;i++){
    iodesc->start[i] = 0;
    iodesc->count[i] = 0;
  }
  
  if(ios->ioproc){
    
    if((iostart != NULL) && (iocount != NULL)){ 
      for(int i=0;i<ndims;i++){
	iodesc->start[i] = iostart[i];
	iodesc->count[i] = iocount[i];
      }
      ios->num_aiotasks = ios->num_iotasks;
    }else{
      ios->num_aiotasks = CalcStartandCount(basetype, ndims, dims, ios->num_iotasks, ios->io_rank,
					    iodesc->start, iodesc->count);

    }
    iosize=1;
    for(int i=0;i<ndims;i++)
      iosize*=iodesc->count[i];

    iodesc->llen = iosize;
    CheckMPIReturn(MPI_Allreduce(&iosize, &(iodesc->maxiobuflen), 1, MPI_INT, MPI_MAX, ios->io_comm),__FILE__,__LINE__);
    
  }

  CheckMPIReturn(MPI_Bcast(&(ios->num_aiotasks), 1, MPI_INT, ios->iomaster,ios->my_comm),__FILE__,__LINE__);


  ierr = box_rearrange_create( ios, maplen, compmap, dims, ndims, ios->num_aiotasks, iodesc);


  lenblocks=1;
  for(int i=0;i<ndims;i++){
    if(iodesc->count[i] == dims[i]){
      lenblocks*=iodesc->count[i];
    }else{
      break;
    }
  }

  *ioidp = pio_add_to_iodesc_list(iodesc);

  
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
  iosys->compmaster = 0;
  iosys->ioproc = false;
#ifndef _MPISERIAL
  iosys->info = MPI_INFO_NULL;
#endif
  iosys->num_iotasks = num_iotasks;
  iosys->num_aiotasks = num_iotasks;

  CheckMPIReturn(MPI_Comm_rank(comp_comm, &(iosys->comp_rank)),__FILE__,__LINE__);
  CheckMPIReturn(MPI_Comm_size(comp_comm, &(iosys->num_comptasks)),__FILE__,__LINE__);
  
  if((num_iotasks < 1) || ((num_iotasks*stride) > iosys->num_comptasks)){
    return PIO_EBADID;
  }
#ifndef _MPISERIAL
  CheckMPIReturn(MPI_Info_create(&(iosys->info)),__FILE__,__LINE__);
#endif
  iosys->ioranks = (int *) calloc(sizeof(int), iosys->num_iotasks);
  for(int i=0;i< num_iotasks; i++){
    iosys->ioranks[i] = (base + i*stride) % iosys->num_comptasks;
    if(iosys->ioranks[i] == iosys->comp_rank)
      iosys->ioproc = true;
  }
  iosys->ioroot = iosys->ioranks[0];
  iosys->iomaster = iosys->ioranks[0];
  iosys->compmaster = 0;

  CheckMPIReturn(MPI_Comm_group(comp_comm, &compgroup),__FILE__,__LINE__);
			
  CheckMPIReturn(MPI_Group_incl(compgroup, num_iotasks, iosys->ioranks, &iogroup),__FILE__,__LINE__);

  CheckMPIReturn(MPI_Comm_create(comp_comm, iogroup, &(iosys->io_comm)),__FILE__,__LINE__);
  if(iosys->ioproc)
    CheckMPIReturn(MPI_Comm_rank(iosys->io_comm, &(iosys->io_rank)),__FILE__,__LINE__);
  else
    iosys->io_rank = -1;

  iosys->union_rank = iosys->comp_rank;


  
  *iosysidp = pio_add_to_iosystem_list(iosys);

  return PIO_NOERR;
}
  
int PIOc_Init_Intracomm_from_F90(int f90_comp_comm, 
			const int num_iotasks, const int stride, 
				 const int base, int *iosysidp){
  return PIOc_Init_Intracomm(MPI_Comm_f2c(f90_comp_comm), num_iotasks, stride,base,iosysidp);
}
  
