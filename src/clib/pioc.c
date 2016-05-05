/**
 * @file 
 * @author Jim Edwards
 * @date  2014
 * @brief PIO C interface 
 *
 * @see http://code.google.com/p/parallelio/
 */


#include <pio.h>
#include <pio_internal.h>


static int counter=0;

/**
 ** @brief Check to see if PIO has been initialized.
 */
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
/**
 ** @brief Check to see if PIO file is open.
 */

int PIOc_File_is_Open(int ncid)
{
  file_desc_t *file;
  file = pio_get_file_from_id(ncid);
  if(file==NULL)
    return 0;
  else
    return 1;
}

/**
 ** @brief Set the error handling method to be used for subsequent 
 ** pio library calls, returns the previous method setting
 */
int PIOc_Set_File_Error_Handling(int ncid, int method)
{
  file_desc_t *file;
  int oldmethod;
  file = pio_get_file_from_id(ncid);
  oldmethod = file->iosystem->error_handler;
  file->iosystem->error_handler = method;
  return(oldmethod);
} 

/**
 ** @brief Increment the unlimited dimension of the given variable
 */
int PIOc_advanceframe(int ncid, int varid)
{
  file_desc_t *file;
  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;

  file->varlist[varid].record++;

  return(PIO_NOERR);
} 

/**
 * @ingroup PIO_setframe 
 * @brief Set the unlimited dimension of the given variable
 * 
 * @param ncid the ncid of the file.
 * @param varid the varid of the variable
 * @param frame the value of the unlimited dimension.  In c 0 for the
 * first record, 1 for the second
 *
 * @return PIO_NOERR for no error, or error code.
 */
int PIOc_setframe(const int ncid, const int varid, const int frame)
{
  file_desc_t *file;
  file = pio_get_file_from_id(ncid);
  if(file == NULL || varid<0 || varid>=PIO_MAX_VARS){
    return PIO_EBADID;
  }

  file->varlist[varid].record = frame;

  return(PIO_NOERR);
} 

/**
 ** @brief Get the number of IO tasks set.
 */
int PIOc_get_numiotasks(int iosysid, int *numiotasks)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  *numiotasks = ios->num_iotasks;

  return PIO_NOERR;

}


/**
 ** @brief Get the IO rank on the current task
 */
int PIOc_get_iorank(int iosysid, int *iorank)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  *iorank = ios->io_rank;

  return PIO_NOERR;

}

/**
 ** @brief Get the local size of the variable
 */

int PIOc_get_local_array_size(int ioid)
{
  io_desc_t *iodesc;
  iodesc = pio_get_iodesc_from_id(ioid);
  return(iodesc->ndof);
}

/**
 ** @ingroup PIO_error_method
 ** @brief Set the error handling method used for subsequent calls
 */

 int PIOc_Set_IOSystem_Error_Handling(int iosysid, int method)
{
  iosystem_desc_t *ios;
  int oldmethod;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios==NULL){
    fprintf(stderr,"%s %d Error setting eh method\n",__FILE__,__LINE__);
    print_trace(stderr);
    return PIO_EBADID;
  }	
  oldmethod = ios->error_handler;
  ios->error_handler = method;
  return(oldmethod);
}  

/**
 ** @ingroup PIO_initdecomp
 ** @brief C interface to the initdecomp
 ** @param  iosysid @copydoc iosystem_desc_t (input)
 ** @param  basetype the basic PIO data type used (input)
 ** @param  ndims the number of dimensions in the variable (input)
 ** @param  dims[] the global size of each dimension (input)
 ** @param  maplen the local length of the compmap array (input)
 ** @param  compmap[] a 1 based array of offsets into the array record on file.  A 0 in this array indicates a value which should not be transfered. (input)
 ** @param ioidp  the io description pointer (output)
 ** @param rearranger the rearranger to be used for this decomp or NULL to use the default (optional input)
 ** @param iostart An optional array of start values for block cyclic decompositions  (optional input)
 ** @param iocount An optional array of count values for block cyclic decompositions  (optional input)
 */


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
  

  if(PIO_Save_Decomps){
    char filename[30];
    if(ios->num_comptasks < 100) {
      sprintf(filename, "piodecomp%2.2dtasks%2.2ddims%2.2d.dat",ios->num_comptasks,ndims,counter);
    }else if(ios->num_comptasks < 10000) {
      sprintf(filename, "piodecomp%4.4dtasks%2.2ddims%2.2d.dat",ios->num_comptasks,ndims,counter);
    }else{
      sprintf(filename, "piodecomp%6.6dtasks%2.2ddims%2.2d.dat",ios->num_comptasks,ndims,counter);
    }
    PIOc_writemap(filename,ndims,dims,maplen,compmap,ios->comp_comm);
    counter++;
  }


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
	  //	  printf("iocount[0] = %ld %ld\n",iocount[0], iocount);
	  iodesc->maxiobuflen=1;
	  for(int i=0;i<ndims;i++){
	    iodesc->firstregion->start[i] = iostart[i];
	    iodesc->firstregion->count[i] = iocount[i];
	    compute_maxIObuffersize(ios->io_comm, iodesc);
	    
	  }
	  iodesc->num_aiotasks = ios->num_iotasks;
	}else{
	  iodesc->num_aiotasks = CalcStartandCount(basetype, ndims, dims, 
						   ios->num_iotasks, ios->io_rank,
						   iodesc->firstregion->start, iodesc->firstregion->count);
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
    /*
    if(ios->ioproc){
      io_region *ioregion = iodesc->firstregion;
      while(ioregion != NULL){
	for(int i=0;i<ndims;i++)
	  printf("%s %d i %d dim %d start %ld count %ld\n",__FILE__,__LINE__,i,dims[i],ioregion->start[i],ioregion->count[i]);
	ioregion = ioregion->next;
      }
    }
    */
  }

  *ioidp = pio_add_to_iodesc_list(iodesc);

  performance_tune_rearranger(*ios, iodesc);
  
  return PIO_NOERR;
}

/**
 ** @ingroup PIO_initdecomp
 ** This is a simplified initdecomp which can be used if the memory order of the data can be 
 ** expressed in terms of start and count on the file.
 ** in this case we compute the compdof and use the subset rearranger
 */


int PIOc_InitDecomp_bc(const int iosysid, const int basetype,const int ndims, const int dims[], 
		       const long int start[], const long int count[], int *ioidp)
		    
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
    if(start[i]<0 || count[i]< 0 || (start[i]+count[i])>dims[i]){
      piodie("Invalid start or count argument ",__FILE__,__LINE__);
    }
  }
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  int n, i, maplen=1;
    
  for( i=0;i<ndims;i++){
    maplen*=count[i];
  }
  PIO_Offset compmap[maplen], prod[ndims], loc[ndims];
    
  prod[ndims-1]=1;
  loc[ndims-1]=0;
  for(n=ndims-2;n>=0;n--){
    prod[n]=prod[n+1]*dims[n+1];  
    loc[n]=0;
  }
  for(i=0;i<maplen;i++){
    compmap[i]=0;
    for(n=ndims-1;n>=0;n--){
      compmap[i]+=(start[n]+loc[n])*prod[n];
    }
    n=ndims-1;
    loc[n]=(loc[n]+1)%count[n];
    while(loc[n]==0 && n>0){
      n--;
      loc[n]=(loc[n]+1)%count[n];
    }
  }
  int rearr = PIO_REARR_SUBSET;
  PIOc_InitDecomp( iosysid, basetype,ndims, dims, 
		   maplen,  compmap, ioidp, &rearr, NULL, NULL);  


  return PIO_NOERR;
}


/** @ingroup PIO_init 
 * Library initialization used when IO tasks are distinct from compute
 * tasks.
 *
 * This is a collective call.  Input parameters are read on
 * comp_rank=0 values on other tasks are ignored.  This variation of
 * PIO_init sets up a distinct set of tasks to handle IO, these tasks
 * do not return from this call.  Instead they go to an internal loop
 * and wait to receive further instructions from the computational
 * tasks.
 *
 * For 4 tasks, to have 2 of them be computational, and 2 of them
 * be IO, I would provide the following:
 *
 * component_count = 1
 *
 * peer_comm = MPI_COMM_WORLD
 *
 * comp_comms = an array with one element, an MPI (intra) communicator
 * that contains the two tasks designated to do computation
 * (processors 0, 1).

 * io_comm = an MPI (intra) communicator with the other two tasks (2,
 * 3).
 *
 * iosysidp = pointer that gets the IO system ID.
 *
 * Fortran function (from PIO1, in piolib_mod.F90) is:
 *
 * subroutine init_intercom(component_count, peer_comm, comp_comms,
 * io_comm, iosystem, rearr_opts)
 *
 * Some notes from Jim:
 *
 * Components and Component Count
 * ------------------------------
 *
 * It's a cesm thing - the cesm model is composed of several component
 * models (atm, ocn, ice, lnd, etc) that may or may not be collocated
 * on mpi tasks.  Since for intercomm the IOCOMM tasks are a subset of
 * the compute tasks for a given component we have a separate iocomm
 * for each model component.  and we call init_inracomm independently
 * for each component.
 *
 * When the IO tasks are independent of any model component then we
 * can have all of the components share one set of iotasks and we call
 * init_intercomm once with the information for all components.
 *
 * Inter vs Intra Communicators
 * ----------------------------
 *
 * ​For an intra you just need to provide the compute comm, pio creates
 * an io comm as a subset of that compute comm.
 *
 * For an inter you need to provide multiple comms - peer comm is the
 * communicator that is going to encompass all of the tasks - usually
 * this will be mpi_comm_world.  Then you need to provide a comm for
 * each component model that will share the io server, then an
 * io_comm.
 *
 * Example of Communicators
 * ------------------------
 *
 * Starting from MPI_COMM_WORLD the calling program will create an
 * IO_COMM and one or more COMP_COMMs, I think an example might be best:
 *
 * Suppose we have 10 tasks and 2 of them will be IO tasks.  Then 0:7
 * are in COMP_COMM and 8:9 are in IO_COMM In this case on tasks 0:7
 * COMP_COMM is defined and IO_COMM is MPI_COMM_NULL and on tasks 8:9
 * IO_COMM is defined and COMP_COMM is MPI_COMM_NULL The communicators
 * to handle communications between COMP_COMM and IO_COMM are defined
 * in init_intercomm and held in a pio internal data structure.
 *
 * Return or Not
 * -------------
 *
 * The io_comm tasks do not return from the init_intercomm routine.   
 *
 * Sequence of Events to do Asynch I/O
 * -----------------------------------
 *
 * Here is the sequence of events that needs to occur when an IO
 * operation is called from the collection of compute tasks.  I'm
 * going to use pio_put_var because write_darray has some special
 * characteristics that make it a bit more complicated...
 *
 * Compute tasks call pio_put_var with an integer argument
 *
 * The MPI_Send sends a message from comp_rank=0 to io_rank=0 on
 * union_comm (a comm defined as the union of io and compute tasks)
 * msg is an integer which indicates the function being called, in
 * this case the msg is PIO_MSG_PUT_VAR_INT
 *
 * The iotasks now know what additional arguments they should expect
 * to receive from the compute tasks, in this case a file handle, a
 * variable id, the length of the array and the array itself.
 *
 * The iotasks now have the information they need to complete the
 * operation and they call the pio_put_var routine.  (In pio1 this bit
 * of code is in pio_get_put_callbacks.F90.in)
 *
 * After the netcdf operation is completed (in the case of an inq or
 * get operation) the result is communicated back to the compute
 * tasks.
 *
 *
 * @param component_count The number of computational (ex. model)
 * components to associate with this IO component
 *
 * @param peer_comm The communicator from which all other communicator
 * arguments are derived
 *
 * @param comp_comms An array containing the computational
 * communicator for each of the computational components. The I/O
 * tasks pass MPI_COMM_NULL for this parameter.
 *
`* @param io_comm The io communicator. Processing tasks pass
 * MPI_COMM_NULL for this parameter.
 *
 * @param iosysidp An array of length component_count. It will get the
 * iosysid for each component.
 *
 * @return PIO_NOERR on success, error code otherwise.
 */
int PIOc_Init_Intercomm(int component_count, MPI_Comm peer_comm,
			MPI_Comm *comp_comms, MPI_Comm io_comm, int *iosysidp)
{
    iosystem_desc_t *iosys;
    iosystem_desc_t *my_iosys;
    int ierr = PIO_NOERR;
    int mpierr;
    int my_rank;
    int iam;
    int io_leader, comp_leader;
    int root;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Allocate struct to hold io system info for each component. */
    if (!(iosys = (iosystem_desc_t *) calloc(1, sizeof(iosystem_desc_t) * component_count)))
	ierr = PIO_ENOMEM;

    if (!ierr)
	for (int cmp = 0; cmp < component_count; cmp++)
	{
	    /* These are used when using the intercomm. */
	    int comp_master = MPI_PROC_NULL, io_master = MPI_PROC_NULL;

	    /* Get a pointer to the iosys struct */
	    my_iosys = &iosys[cmp];
    	    
	    /* Create an MPI info object. */
	    CheckMPIReturn(MPI_Info_create(&(my_iosys->info)),__FILE__,__LINE__);

	    /* This task is part of the computation communicator. */
	    if (comp_comms[cmp] != MPI_COMM_NULL)
	    {
		/* Copy the computation communicator. */
		mpierr = MPI_Comm_dup(comp_comms[cmp], &my_iosys->comp_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Create an MPI group with the computation tasks. */
		mpierr = MPI_Comm_group(my_iosys->comp_comm, &my_iosys->compgroup);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find out how many tasks are in this communicator. */
		mpierr = MPI_Comm_size(iosys->comp_comm, &my_iosys->num_comptasks);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Set the rank within the comp_comm. */
		mpierr = MPI_Comm_rank(my_iosys->comp_comm, &my_iosys->comp_rank);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the io leader in peer_comm. */
		iam = -1;
		mpierr = MPI_Allreduce(&iam, &io_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the comp leader in peer_comm. */
		if (!my_iosys->comp_rank)
		{
		    mpierr = MPI_Comm_rank(peer_comm, &iam);
		    CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		    if (mpierr)
			ierr = PIO_EIO;
		}
		else
		    iam = -1;

		/* Find the lucky comp_leader task. */
		mpierr = MPI_Allreduce(&iam, &comp_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Is this the compmaster? Only if the comp_rank is zero. */
		if (!my_iosys->comp_rank)
		{
		    my_iosys->compmaster = MPI_ROOT;
		    comp_master = MPI_ROOT;
		}
		else
		    my_iosys->compmaster = MPI_PROC_NULL;
		
		/* Set up the intercomm from the computation side. */
		mpierr = MPI_Intercomm_create(my_iosys->comp_comm, 0, peer_comm,
					      io_leader, cmp, &my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Create the union communicator. */
		mpierr = MPI_Intercomm_merge(my_iosys->intercomm, 0, &my_iosys->union_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;
	    }
	    else
	    {
		my_iosys->comp_comm = MPI_COMM_NULL;
		my_iosys->compgroup = MPI_GROUP_NULL;
		my_iosys->comp_rank = -1;
	    }

	    /* This task is part of the IO communicator, so set up the
	     * IO stuff. */
	    if (io_comm != MPI_COMM_NULL)
	    {
		/* Copy the IO communicator. */
		mpierr = MPI_Comm_dup(io_comm, &my_iosys->io_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Get an MPI group that includes the io tasks. */
		mpierr = MPI_Comm_group(my_iosys->io_comm, &my_iosys->iogroup);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find out how many tasks are in this communicator. */
		mpierr = MPI_Comm_size(iosys->io_comm, &my_iosys->num_iotasks);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Set the rank within the io_comm. */
		mpierr = MPI_Comm_rank(my_iosys->io_comm, &my_iosys->io_rank);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the io leader in peer_comm. */
		if (!my_iosys->io_rank)
		{
		    mpierr = MPI_Comm_rank(peer_comm, &iam);
		    CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		    if (mpierr)
			ierr = PIO_EIO;
		}
		else
		    iam = -1;

		/* Find the lucky io_leader task. */
		mpierr = MPI_Allreduce(&iam, &io_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Find the rank of the comp leader in peer_comm. */
		iam = -1;
		mpierr = MPI_Allreduce(&iam, &comp_leader, 1, MPI_INT, MPI_MAX, peer_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;
		
		/* This is an io task. */
		my_iosys->ioproc = true;

		/* Is this the iomaster? Only if the io_rank is zero. */
		if (!my_iosys->io_rank)
		{
		    my_iosys->iomaster = MPI_ROOT;
		    io_master = MPI_ROOT;
		}
		else
		    my_iosys->iomaster = 0;		

		/* Set up the intercomm from the I/O side. */
		mpierr = MPI_Intercomm_create(my_iosys->io_comm, 0, peer_comm,
					      comp_leader, cmp, &my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

		/* Create the union communicator. */
		mpierr = MPI_Intercomm_merge(my_iosys->intercomm, 0, &my_iosys->union_comm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);		
		if (mpierr)
		    ierr = PIO_EIO;

	    }
	    else
	    {
		my_iosys->io_comm = MPI_COMM_NULL;
		my_iosys->iogroup = MPI_GROUP_NULL;
		my_iosys->io_rank = -1;
		my_iosys->ioproc = false;
		my_iosys->iomaster = false;
	    }

	    /* my_comm points to the union communicator for async, and
	     * the comp_comm for non-async. It should not be freed
	     * since it is not a proper copy of the commuicator, just
	     * a copy of the reference to it. */
	    my_iosys->my_comm = my_iosys->union_comm;

	    /* Find rank in union communicator. */
	    mpierr = MPI_Comm_rank(my_iosys->union_comm, &my_iosys->union_rank);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);		
	    if (mpierr)
		ierr = PIO_EIO;

	    /* Find the rank of the io leader in the union communicator. */
	    if (!my_iosys->io_rank)
		my_iosys->ioroot = my_iosys->union_rank;
	    else
		my_iosys->ioroot = -1;

	    /* Distribute the answer to all tasks. */
	    mpierr = MPI_Allreduce(&my_iosys->ioroot, &root, 1, MPI_INT, MPI_MAX,
				   my_iosys->union_comm);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);		
	    if (mpierr)
		ierr = PIO_EIO;
	    my_iosys->ioroot = root;
	    
	    /* Find the rank of the computation leader in the union
	     * communicator. */
	    if (!my_iosys->comp_rank)
		my_iosys->comproot = my_iosys->union_rank;
	    else
		my_iosys->comproot = -1;

	    /* Distribute the answer to all tasks. */
	    mpierr = MPI_Allreduce(&my_iosys->comproot, &root, 1, MPI_INT, MPI_MAX,
				   my_iosys->union_comm);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);		
	    if (mpierr)
		ierr = PIO_EIO;
	    my_iosys->comproot = root;

	    /* Send the number of tasks in the IO and computation
	       communicators to each other over the intercomm. This is
	       a one-to-all bcast from the local task that passes
	       MPI_ROOT as the root (all other local tasks should pass
	       MPI_PROC_NULL as the root). The bcast is recieved by
	       all the members of the leaf group which each pass the
	       rank of the root relative to the root group. */
	    if (io_comm != MPI_COMM_NULL)
	    {
		comp_master = 0;
		mpierr = MPI_Bcast(&my_iosys->num_comptasks, 1, MPI_INT, comp_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);				
		mpierr = MPI_Bcast(&my_iosys->num_iotasks, 1, MPI_INT, io_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    }
	    else
	    {
		io_master = 0;
		mpierr = MPI_Bcast(&my_iosys->num_comptasks, 1, MPI_INT, comp_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);				
		mpierr = MPI_Bcast(&my_iosys->num_iotasks, 1, MPI_INT, io_master,
				   my_iosys->intercomm);
		CheckMPIReturn(mpierr, __FILE__, __LINE__);
	    }

	    /* Allocate an array to hold the ranks of the IO tasks
	     * within the union communicator. */
	    printf("%d allocating for %d iotasks\n", my_rank, my_iosys->num_iotasks);
	    if (!(my_iosys->ioranks = malloc(my_iosys->num_iotasks * sizeof(int))))
		return PIO_ENOMEM;
	    printf("%d allocated\n", my_rank);
	    
	    /* Allocate a temp array to help get the IO ranks. */
	    int *tmp_ioranks;
	    if (!(tmp_ioranks = malloc(my_iosys->num_iotasks * sizeof(int))))
		return PIO_ENOMEM;

	    /* Init array, then have IO tasks set their values, then
	     * use allreduce to distribute results to all tasks. */
	    for (int cnt = 0 ; cnt < my_iosys->num_iotasks; cnt++)
		tmp_ioranks[cnt] = -1;		
	    if (io_comm != MPI_COMM_NULL)
		tmp_ioranks[my_iosys->io_rank] = my_iosys->union_rank;
	    mpierr = MPI_Allreduce(tmp_ioranks, my_iosys->ioranks, my_iosys->num_iotasks, MPI_INT, MPI_MAX,
				   my_iosys->union_comm);
	    CheckMPIReturn(mpierr, __FILE__, __LINE__);

	    /* Free temp array. */
	    free(tmp_ioranks);

	    /* Set the default error handling. */
	    my_iosys->error_handler = PIO_INTERNAL_ERROR;

	    /* We do support asynch interface. */
	    my_iosys->async_interface = true;

	    /* For debug purposes, print the contents of the struct. */
	    for (int t = 0; t < my_iosys->num_iotasks + my_iosys->num_comptasks; t++)
	    {
		MPI_Barrier(my_iosys->union_comm);
		if (my_rank == t)
		    pio_iosys_print(my_rank, my_iosys);
	    }

	    /* Add this id to the list of PIO iosystem ids. */
	    iosysidp[cmp] = pio_add_to_iosystem_list(my_iosys);
	    printf("%d added to iosystem_list iosysid = %d\n", my_rank, iosysidp[cmp]);	    

	    /* Now call the function from which the IO tasks will not
	     * return until the PIO_MSG_EXIT message is sent. */
	    if (io_comm != MPI_COMM_NULL)
	    {
		printf("%d about to call pio_msg_handler\n", my_rank);
		if ((ierr = pio_msg_handler(my_iosys->io_rank, component_count, iosys)))
		    return ierr;
	    }
    
	}

    /* If there was an error, make sure all tasks see it. */
    if (ierr)
    {
	/*mpierr = MPI_Bcast(&ierr, 1, MPI_INT, iosys->my_rank, iosys->intercomm);*/
	mpierr = MPI_Bcast(&ierr, 1, MPI_INT, 0, iosys->intercomm);
	CheckMPIReturn(mpierr, __FILE__, __LINE__);
	if (mpierr)
	    ierr = PIO_EIO;
    }
    
    return ierr;
}

/** 
 ** @ingroup PIO_init
 ** @brief library initialization used when IO tasks are a subset of compute tasks
 ** @param comp_comm the MPI_Comm of the compute tasks
 ** @param num_iotasks the number of io tasks to use
 ** @param stride the offset between io tasks in the comp_comm
 ** @param base the comp_comm index of the first io task
 ** @param rearr the rearranger to use by default, this may be overriden in the @ref PIO_initdecomp
 ** @param iosysidp index of the defined system descriptor
 */

int PIOc_Init_Intracomm(const MPI_Comm comp_comm, 
			const int num_iotasks, const int stride, 
			const int base,const int rearr, int *iosysidp)
{
  iosystem_desc_t *iosys;
  int ierr = PIO_NOERR;
  int ustride;
  int lbase;
  int mpierr;

  iosys = (iosystem_desc_t *) malloc(sizeof(iosystem_desc_t));

  /* Copy the computation communicator into union_comm. */
  mpierr = MPI_Comm_dup(comp_comm, &iosys->union_comm);
  CheckMPIReturn(mpierr, __FILE__, __LINE__);		
  if (mpierr)
      ierr = PIO_EIO;

  /* Copy the computation communicator into comp_comm. */
  if (!ierr)
  {
      mpierr = MPI_Comm_dup(comp_comm, &iosys->comp_comm);
      CheckMPIReturn(mpierr, __FILE__, __LINE__);		
      if (mpierr)
	  ierr = PIO_EIO;
  }

  if (!ierr)
  {
      iosys->my_comm = iosys->comp_comm;
      iosys->io_comm = MPI_COMM_NULL;
      iosys->intercomm = MPI_COMM_NULL;
      iosys->error_handler = PIO_INTERNAL_ERROR;
      iosys->async_interface= false;
      iosys->compmaster = 0;
      iosys->iomaster = 0;
      iosys->ioproc = false;
      iosys->default_rearranger = rearr;
      iosys->num_iotasks = num_iotasks;

      ustride = stride;

      /* Find MPI rank and number of tasks in comp_comm communicator. */
      CheckMPIReturn(MPI_Comm_rank(iosys->comp_comm, &(iosys->comp_rank)),__FILE__,__LINE__);
      CheckMPIReturn(MPI_Comm_size(iosys->comp_comm, &(iosys->num_comptasks)),__FILE__,__LINE__);
      if(iosys->comp_rank==0)
	  iosys->compmaster = MPI_ROOT;  

      /* Ensure that settings for number of computation tasks, number
       * of IO tasks, and the stride are reasonable. */
      if((iosys->num_comptasks == 1) && (num_iotasks*ustride > 1)) {
	  // This is a serial run with a bad configuration. Set up a single task.
	  fprintf(stderr, "PIO_TP PIOc_Init_Intracomm reset stride and tasks.\n");
	  iosys->num_iotasks = 1;
	  ustride = 1;
      }
      if((iosys->num_iotasks < 1) || ((iosys->num_iotasks*ustride) > iosys->num_comptasks)){
	  fprintf(stderr, "PIO_TP PIOc_Init_Intracomm error\n");
	  fprintf(stderr, "num_iotasks=%d, ustride=%d, num_comptasks=%d\n", num_iotasks, ustride, iosys->num_comptasks);
	  return PIO_EBADID;
      }

      /* Create an array that holds the ranks of the tasks to be used for IO. */
      iosys->ioranks = (int *) calloc(sizeof(int), iosys->num_iotasks);
      for(int i=0;i< iosys->num_iotasks; i++){
	  iosys->ioranks[i] = (base + i*ustride) % iosys->num_comptasks;
	  if(iosys->ioranks[i] == iosys->comp_rank)
	      iosys->ioproc = true;
      }
      iosys->ioroot = iosys->ioranks[0];

      /* Create an MPI info object. */
      CheckMPIReturn(MPI_Info_create(&(iosys->info)),__FILE__,__LINE__);
      iosys->info = MPI_INFO_NULL;

      if(iosys->comp_rank == iosys->ioranks[0])
	  iosys->iomaster = MPI_ROOT;

      /* Create a group for the computation tasks. */
      CheckMPIReturn(MPI_Comm_group(iosys->comp_comm, &(iosys->compgroup)),__FILE__,__LINE__);
			
      /* Create a group for the IO tasks. */
      CheckMPIReturn(MPI_Group_incl(iosys->compgroup, iosys->num_iotasks, iosys->ioranks,
				    &(iosys->iogroup)),__FILE__,__LINE__);

      /* Create an MPI communicator for the IO tasks. */
      CheckMPIReturn(MPI_Comm_create(iosys->comp_comm, iosys->iogroup, &(iosys->io_comm)),__FILE__,__LINE__);

      /* For the tasks that are doing IO, get their rank. */
      if(iosys->ioproc)
	  CheckMPIReturn(MPI_Comm_rank(iosys->io_comm, &(iosys->io_rank)),__FILE__,__LINE__);
      else
	  iosys->io_rank = -1;

      iosys->union_rank = iosys->comp_rank;

      /* Add this iosys struct to the list in the PIO library. */
      *iosysidp = pio_add_to_iosystem_list(iosys);

      pio_get_env();

      /* allocate buffer space for compute nodes */
      compute_buffer_init(*iosys);
  }
  
  return ierr;
}

/**
 ** @internal 
 ** interface to call from pio_init from fortran
 ** @endinternal
 */
int PIOc_Init_Intracomm_from_F90(int f90_comp_comm, 
			const int num_iotasks, const int stride, 
				 const int base, const int rearr, int *iosysidp){
  return PIOc_Init_Intracomm(MPI_Comm_f2c(f90_comp_comm), num_iotasks, stride,base,rearr, iosysidp);
}
  
/**
 ** @brief Send a hint to the MPI-IO library 
 **
 */
int PIOc_set_hint(const int iosysid, char hint[], const char hintval[])
{
  iosystem_desc_t *ios;

  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  if(ios->ioproc)
    CheckMPIReturn( MPI_Info_set(ios->info, hint, hintval), __FILE__,__LINE__);

  return PIO_NOERR;

}

/** @ingroup PIO_finalize 
 * Clean up internal data structures, free MPI resources, and exit the
 * pio library.
 *
 * @param iosysid: the io system ID provided by PIOc_Init_Intracomm().
 *
 * @returns 0 for success or non-zero for error.
 */

int PIOc_finalize(const int iosysid)
{
  iosystem_desc_t *ios, *nios;
  int msg;
  int mpierr;

  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  
  /* If asynch IO is in use, send the PIO_MSG_EXIT message from the
   * comp master to the IO processes. */
  if (ios->async_interface && !ios->comp_rank)
  {
    msg = PIO_MSG_EXIT;
    mpierr = MPI_Send(&msg, 1, MPI_INT, ios->ioroot, 1, ios->union_comm);
    CheckMPIReturn(mpierr, __FILE__, __LINE__);		      
  }

  /* Free this memory that was allocated in init_intracomm. */
  if (ios->ioranks)
      free(ios->ioranks);

  /* Free the buffer pool. */
  free_cn_buffer_pool(*ios);

  /* Free the MPI groups. */
  if (ios->compgroup != MPI_GROUP_NULL)
    MPI_Group_free(&ios->compgroup);

  if (ios->iogroup != MPI_GROUP_NULL)
    MPI_Group_free(&(ios->iogroup));

  /* Free the MPI communicators. my_comm is just a copy (but not an
   * MPI copy), so does not have to have an MPI_Comm_free() call. */
  if(ios->intercomm != MPI_COMM_NULL){
    MPI_Comm_free(&(ios->intercomm));
  }
  if(ios->io_comm != MPI_COMM_NULL){
    MPI_Comm_free(&(ios->io_comm));
  }
  if(ios->comp_comm != MPI_COMM_NULL){
    MPI_Comm_free(&(ios->comp_comm));
  }
  if(ios->union_comm != MPI_COMM_NULL){
    MPI_Comm_free(&(ios->union_comm));
  }

  /* Delete the iosystem_desc_t data associated with this id. */
  return pio_delete_iosystem_from_list(iosysid);
}

/**
 ** @brief return a logical indicating whether this task is an iotask 
 */
int PIOc_iam_iotask(const int iosysid, bool *ioproc)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;
  
  *ioproc = ios->ioproc;
  return PIO_NOERR;
}

/**
 ** @brief return the rank of this task in the io comm or
 ** -1 if this task is not in the comm
 */
int PIOc_iotask_rank(const int iosysid, int *iorank)
{
  iosystem_desc_t *ios;
  ios = pio_get_iosystem_from_id(iosysid);
  if(ios == NULL)
    return PIO_EBADID;

  *iorank = ios->io_rank;

  return PIO_NOERR;
  
}

/**
 ** @brief return true if this iotype is supported in the build, 0 otherwise
 */

int PIOc_iotype_available(const int iotype)
{

  switch(iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
    case PIO_IOTYPE_NETCDF4C:
      return(1);
#endif
    case PIO_IOTYPE_NETCDF:
      return(1);
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      return(1);
      break;
#endif
    default:
      return(0);
    }

}
