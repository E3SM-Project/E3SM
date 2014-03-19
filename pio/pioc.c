#include <pio.h>
#include <pio_internal.h>

int PIOc_iosystem_is_active(const int iosysid, bool *active)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  
  if(ios->comp_comm == MPI_COMM_NULL && ios->io_comm == MPI_COMM_NULL){
    *active = false;
  }else{
    *active = true;
  }
  return PIO_NOERR;
}

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

int PIOc_advanceframe(int ncid, int varid)
{
  file_desc_t *file;
  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;

  file->varlist[varid].record++;

  return(PIO_NOERR);
} 

int PIOc_setframe(const int ncid, const int varid,const int frame)
{
  file_desc_t *file;
  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;

  file->varlist[varid].record = frame;

  return(PIO_NOERR);
} 


int PIOc_get_numiotasks(int iosysid, int *numiotasks)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  *numiotasks = ios->num_iotasks;

  return PIO_NOERR;

}


int PIOc_get_iorank(int iosysid, int *iorank)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  *iorank = ios->io_rank;

  return PIO_NOERR;

}


int PIOc_get_local_array_size(int ioid)
{
  io_desc_t *iodesc;
  iodesc = pio_get_iodesc_from_id(ioid);
  return(iodesc->ndof);
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
		    const int maplen, const PIO_Offset *compmap, int *ioidp,const int *rearranger,  
		    const PIO_Offset *iostart,const PIO_Offset *iocount)
{
  iosystem_desc_t *ios;
  io_desc_t *iodesc;
  int mpierr;
  int ierr;
  int iosize;
  int ndisp;

  for(int i=0;i<ndims;i++){
    if(dims[i]<=0){
      piodie("Invalid dims argument",__FILE__,__LINE__);
    }
  }
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  
  iodesc = malloc_iodesc(basetype, ndims);

  if(rearranger == NULL)
    iodesc->rearranger = ios->default_rearranger;
  else
    iodesc->rearranger = *rearranger;

    
  if(iodesc->rearranger==PIO_REARR_SUBSET){
    if((iostart != NULL) && (iocount != NULL)){ 
      fprintf(stderr,"%s %s\n","Iostart and iocount arguments to PIOc_InitDecomp",
	      "are incompatable with subset rearrange method and will be ignored");
    }
    iodesc->num_aiotasks = ios->num_iotasks;
    ierr = subset_rearrange_create( *ios, maplen, compmap, dims, ndims, iodesc);
  }else{
      if(ios->ioproc){
      //  Unless the user specifies the start and count for each IO task compute it.    
	if((iostart != NULL) && (iocount != NULL)){ 
	  printf("iocount[0] = %ld %ld\n",iocount[0], iocount);
	  for(int i=0;i<ndims;i++){
	    iodesc->start[i] = iostart[i];
	    iodesc->count[i] = iocount[i];
	  }
	  iodesc->num_aiotasks = ios->num_iotasks;
	}else{
	
	  iodesc->num_aiotasks = CalcStartandCount(basetype, ndims, dims, 
						   ios->num_iotasks, ios->io_rank,
						   iodesc->start, iodesc->count);

	//	for(int j=0;j<ndims;j++)
	//	  printf("%d start[%d] %ld count %ld\n",ios->io_rank,j,iodesc->start[j],iodesc->count[j]);





      }
      compute_maxIObuffersize(ios->io_comm, iodesc);

    }
    // Depending on array size and io-blocksize the actual number of io tasks used may vary
    CheckMPIReturn(MPI_Bcast(&(iodesc->num_aiotasks), 1, MPI_INT, ios->ioroot,
			     ios->my_comm),__FILE__,__LINE__);
    // Compute the communications pattern for this decomposition
    if(iodesc->rearranger==PIO_REARR_BOX){   
      ierr = box_rearrange_create( *ios, maplen, compmap, dims, ndims, iodesc);
    }

    //  if(ios->ioproc)
    //   for(int i=0;i<ndims;i++)
    //     printf("%s %d dim %d start %ld count %ld\n",__FILE__,__LINE__,dims[i],iodesc->start[i],iodesc->count[i]);
  }

  *ioidp = pio_add_to_iodesc_list(iodesc);
  
  return PIO_NOERR;
}

int PIOc_Init_Intracomm(const MPI_Comm comp_comm, 
			const int num_iotasks, const int stride, 
			const int base,const int rearr, int *iosysidp)
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
  iosys->compmaster = false;
  iosys->iomaster = false;
  iosys->ioproc = false;
  iosys->default_rearranger = rearr;
  iosys->num_iotasks = num_iotasks;

  CheckMPIReturn(MPI_Comm_rank(comp_comm, &(iosys->comp_rank)),__FILE__,__LINE__);
  CheckMPIReturn(MPI_Comm_size(comp_comm, &(iosys->num_comptasks)),__FILE__,__LINE__);
  if((num_iotasks < 1) || ((num_iotasks*stride) > iosys->num_comptasks)){
    return PIO_EBADID;
  }
  if(iosys->comp_rank==0)
    iosys->compmaster = true;

#ifndef _MPISERIAL
  CheckMPIReturn(MPI_Info_create(&(iosys->info)),__FILE__,__LINE__);
  iosys->info = MPI_INFO_NULL;
#endif
  iosys->ioranks = (int *) calloc(sizeof(int), iosys->num_iotasks);
  for(int i=0;i< num_iotasks; i++){
    iosys->ioranks[i] = (base + i*stride) % iosys->num_comptasks;
    if(iosys->ioranks[i] == iosys->comp_rank)
      iosys->ioproc = true;
  }
  iosys->ioroot = iosys->ioranks[0];

  if(iosys->comp_rank == iosys->ioranks[0])
    iosys->iomaster = true;

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
				 const int base, const int rearr, int *iosysidp){
  return PIOc_Init_Intracomm(MPI_Comm_f2c(f90_comp_comm), num_iotasks, stride,base,rearr, iosysidp);
}
  

int PIOc_set_hint(const int iosysid, const char hint[], const char hintval[])
{
  iosystem_desc_t *ios;

  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  
  if(ios->ioproc)
    CheckMPIReturn( MPI_Info_set(ios->info, hint, hintval), __FILE__,__LINE__);

  return PIO_NOERR;

}

int PIOc_finalize(const int iosysid)
{
  iosystem_desc_t *ios, *nios;

  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  nios = ios;
  while(nios != NULL){
    if(nios->ioranks != NULL){
      free(nios->ioranks);
      nios->ioranks=NULL;
    }
    ios = nios;
    nios = nios->next;
    free(ios);
  }
  return PIO_NOERR;
}


int PIOc_iam_iotask(const int iosysid, bool *ioproc)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  
  *ioproc = ios->ioproc;
  return PIO_NOERR;
}

int PIOc_iotask_rank(const int iosysid, int *iorank)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  *iorank = ios->io_rank;

  return PIO_NOERR;
  
}

