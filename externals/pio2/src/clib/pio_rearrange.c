////
/// @file pio_rearrange.c
/// @author Jim Edwards
/// @date 2014
/// @brief Code to map IO to model decomposition
///
/// 
/// 
/// 
/// @see  http://code.google.com/p/parallelio/
////
#include <pio.h>
#include <pio_internal.h>

/** internal variable used for debugging */
int tmpioproc=-1;  


/** @internal
 ** Convert an index into a list of dimensions. E.g., for index 4 into a
 ** array defined as a[3][2], will return 1 1.
 ** @endinternal
*/
void idx_to_dim_list(const int ndims, const int gdims[], const PIO_Offset idx,
                     PIO_Offset dim_list[])
{
  int i, curr_idx, next_idx;
  curr_idx = idx;
  // Easiest to start from the right and move left.
  for (i = ndims-1; i >= 0; --i) {
    // This way of doing div/mod is slightly faster than using "/" and "%".
    next_idx = curr_idx / gdims[i];
    dim_list[i] = curr_idx - (next_idx*gdims[i]);
    curr_idx = next_idx;
  }
}
/**
 ** @internal
 ** Expand a region along dimension dim, by incrementing count[i] as much as
 ** possible, consistent with the map.
 **
 ** Once max_size is reached, the map is exhausted, or the next entries fail
 ** to match, expand_region updates the count and calls itself with the next
 ** outermost dimension, until the region has been expanded as much as
 ** possible along all dimensions.
 ** @endinternal
 */
void expand_region(const int dim, const int gdims[], const int maplen,
                   const PIO_Offset map[], const int region_size,
                   const int region_stride, const int max_size[],
                   PIO_Offset count[])
{
  int i, j, test_idx, expansion_done;
  // Precondition: maplen >= region_size (thus loop runs at least once).

  // Flag used to signal that we can no longer expand the region along
  // dimension dim.
  expansion_done = 0;

  // Expand no greater than max_size along this dimension.
  for (i = 1; i <= max_size[dim]; ++i) {
    // Count so far is at least i.
    count[dim] = i;

    // Now see if we can expand to i+1 by checking that the next
    // region_size elements are ahead by exactly region_stride.
    // Assuming monotonicity in the map, we could skip this for the
    // innermost dimension, but it's necessary past that because the
    // region does not necessarily comprise contiguous values.
    for (j = 0; j < region_size; ++j) {
      test_idx = j + i*region_size;
      // If we have exhausted the map, or the map no longer matches,
      // we are done, break out of both loops.
      if (test_idx >= maplen || map[test_idx] != map[j] + i*region_stride) {
        expansion_done = 1;
        break;
      }
    }
    if (expansion_done) break;
  }

  // Move on to next outermost dimension if there are more left, else return.
  if (dim > 0) {
    expand_region(dim-1, gdims, maplen, map, region_size*count[dim],
                  region_stride*gdims[dim], max_size, count);
  }
}
/**
 ** @internal
 ** Set start and count so that they describe the first region in map.
 ** @endinternal
 */
PIO_Offset find_region(const int ndims, const int gdims[],
                       const int maplen, const PIO_Offset map[],
                       PIO_Offset start[], PIO_Offset count[])
{
  int dim;
  int max_size[ndims];
  PIO_Offset regionlen=1;
  // Preconditions (which might be useful to check/assert):
  //   ndims is > 0
  //   maplen is > 0
  //   all elements of map are inside the bounds specified by gdims
  // The map array is 1 based, but calculations are 0 based
  idx_to_dim_list(ndims, gdims, map[0]-1, start);

  for (dim = 0; dim < ndims; ++dim) {
    // Can't expand beyond the array edge.
    max_size[dim] = gdims[dim] - start[dim];
  }

  // For each dimension, figure out how far we can expand in that dimension
  // while staying contiguous in the input array.
  //
  // Start with the innermost dimension (ndims-1), and it will recurse
  // through to the outermost dimensions.
  expand_region(ndims-1, gdims, maplen, map, 1, 1, max_size, count);

  for(dim=0;dim<ndims;dim++)
    regionlen*=count[dim];
  return(regionlen);

}

/**
 ** @internal
** Convert a global coordinate value into a local array index
** @endinternal
*/
PIO_Offset coord_to_lindex(const int ndims, const PIO_Offset lcoord[], const PIO_Offset count[])
{
  PIO_Offset lindex=0;
  PIO_Offset stride=1;

  for(int i=ndims-1; i>=0; i--){
    lindex += lcoord[i]*stride;
    stride = stride*count[i];
  }
  return lindex;

}

/**
 ** @internal
** Compute the max io buffersize needed for a given variable
** @endinternal
*/
void compute_maxIObuffersize(MPI_Comm io_comm, io_desc_t *iodesc)
{
  PIO_Offset iosize, totiosize;
  int i;
  io_region *region;

  //  compute the max io buffer size, for conveneance it is the combined size of all regions
  totiosize=0;
  region = iodesc->firstregion;
  while(region != NULL){
    if(region->count[0]>0){
      iosize=1;
      for(i=0;i<iodesc->ndims;i++)
	iosize*=region->count[i];
      totiosize+=iosize;
    }
    region = region->next;
  }
  // Share the max io buffer size with all io tasks
  CheckMPIReturn(MPI_Allreduce(MPI_IN_PLACE, &totiosize, 1, MPI_OFFSET, MPI_MAX, io_comm),__FILE__,__LINE__);
  iodesc->maxiobuflen = totiosize;
  if(iodesc->maxiobuflen<=0  ){
    fprintf(stderr,"%s %d %ld %ld %d %d %d\n",__FILE__,__LINE__,iodesc->maxiobuflen,totiosize,MPI_OFFSET,MPI_MAX,io_comm); 
    piodie("ERROR: maxiobuflen<=0",__FILE__,__LINE__);
  }

}
/**
 ** @internal
** Create the derived MPI datatypes used for comp2io and io2comp transfers
** @param basetype The type of data (int,real,double). 
** @param msgcnt The number of MPI messages/tasks to use.
** @param dlen The length of the data array.
** @param mindex An array of indexes into the data array from the comp map
** @param mcount The number of indexes to be put on each mpi message/task
** @param *mfrom A pointer to the previous structure in the read/write list
** @param mtype The final data structure sent through MPI to the read/write
** @endinternal
*/
int create_mpi_datatypes(const MPI_Datatype basetype,const int msgcnt,const PIO_Offset dlen, const PIO_Offset mindex[],const int mcount[],
			 int *mfrom, MPI_Datatype mtype[])
{
  PIO_Offset bsizeT[msgcnt];
  int blocksize;
  int numinds;
  PIO_Offset *lindex = NULL;

  numinds=0;
  for(int j=0;j<msgcnt;j++){
    numinds+=mcount[j];
  }

  pioassert(dlen>=0,"dlen < 0",__FILE__,__LINE__);
  pioassert(numinds>=0, "num inds < 0",__FILE__,__LINE__);


  if(mindex != NULL){
    // memcpy(lindex, mindex, (size_t) (dlen*sizeof(PIO_Offset)));
    lindex = malloc(numinds*sizeof(PIO_Offset));
    memcpy(lindex, mindex, (size_t) (numinds*sizeof(PIO_Offset)));
  }

  bsizeT[0]=0;
  mtype[0] = MPI_DATATYPE_NULL;
  int pos = 0;
  int ii = 0;
  if(msgcnt>0){

     if(mfrom == NULL){
      for(int i=0;i<msgcnt;i++){
	if(mcount[i]>0){
	  bsizeT[ii] = GCDblocksize(mcount[i], lindex+pos);
	  ii++;
	  pos+=mcount[i];

	}
      }
      blocksize = (int) lgcd_array(ii ,bsizeT);
    }else{
      blocksize=1;
    }

    pos = 0;
    for(int i=0;i< msgcnt; i++){
      if(mcount[i]>0){
	int len = mcount[i]/blocksize;
	int displace[len];
	if(blocksize==1){
	  if(mfrom == NULL){
	    for(int j=0;j<len;j++) {
		displace[j] = (int) (lindex[pos+j]);
	    }
	  }else{
	    int k=0;
	    for(int j=0;j<numinds;j++)
	      if(mfrom[j]==i) {
		displace[k++] = (int) (lindex[j]);
	      }
	  }
	    
	}else{
	  for(int j=0;j<mcount[i];j++)
	    (lindex+pos)[j]++;
	  for(int j=0;j<len;j++){
	    displace[j]= ((lindex+pos)[j*blocksize]-1);
	  }
	}

	CheckMPIReturn(MPI_Type_create_indexed_block(len, blocksize, displace, basetype, mtype+i),__FILE__,__LINE__);
	if(mtype[i] == MPI_DATATYPE_NULL){
	  piodie("Unexpected NULL MPI DATATYPE",__FILE__,__LINE__);
	}
	CheckMPIReturn(MPI_Type_commit(mtype+i), __FILE__,__LINE__);
	pos+=mcount[i];
      }
    }
    if (lindex)
	free(lindex);

  }

  return PIO_NOERR;

}

/**
 ** @internal
** Create the derived MPI datatypes used for comp2io and io2comp transfers
** @endinternal
*/
int define_iodesc_datatypes(const iosystem_desc_t ios, io_desc_t *iodesc)
{
  int i;
  if(ios.ioproc){
    if(iodesc->rtype==NULL){
      if(iodesc->nrecvs>0){
	iodesc->rtype = (MPI_Datatype *) bget(iodesc->nrecvs * sizeof(MPI_Datatype));
	if(iodesc->rtype == NULL){
	  piomemerror(ios,iodesc->nrecvs * sizeof(MPI_Datatype), __FILE__,__LINE__);
	}
	for(i=0; i<iodesc->nrecvs; i++){
	  iodesc->rtype[i] = MPI_DATATYPE_NULL;
	}
	if(iodesc->rearranger==PIO_REARR_SUBSET){
	  create_mpi_datatypes(iodesc->basetype, iodesc->nrecvs, iodesc->llen, iodesc->rindex, iodesc->rcount, iodesc->rfrom, iodesc->rtype);
	}else{
	  create_mpi_datatypes(iodesc->basetype, iodesc->nrecvs, iodesc->llen, iodesc->rindex, iodesc->rcount, NULL, iodesc->rtype);
	}
      }
            /*      if(tmpioproc==95)     {
	MPI_Aint lb;
	MPI_Aint extent;
	for(i=0;i<ntypes;i++){
	  MPI_Type_get_extent(iodesc->rtype[i], &lb, &extent);
	  printf("%s %d %d %d %d \n",__FILE__,__LINE__,i,lb,extent);
	  
	}
		}
	      */
    }
  }


  if(iodesc->stype==NULL){
    int ntypes;
    int ncnt;
    if(iodesc->rearranger==PIO_REARR_SUBSET){
      ntypes = 1;
      ncnt = iodesc->scount[0];
    }else{
      ntypes = ios.num_iotasks;
      ncnt = iodesc->ndof;
    }

    //  printf("COMP: %d\n",ntypes);


    iodesc->stype = (MPI_Datatype *) bget(ntypes * sizeof(MPI_Datatype));
    if(iodesc->stype == NULL){
      piomemerror(ios,ntypes * sizeof(MPI_Datatype), __FILE__,__LINE__);
    }
    for(i=0; i<ntypes; i++){
      iodesc->stype[i] = MPI_DATATYPE_NULL;
    }
    iodesc->num_stypes = ntypes;

    create_mpi_datatypes(iodesc->basetype, ntypes, ncnt, iodesc->sindex, iodesc->scount, NULL, iodesc->stype);
       /*    if(tmpioproc==95)   {
      MPI_Aint lb;
      MPI_Aint extent;
      for(i=0;i<ntypes;i++){
	MPI_Type_get_extent(iodesc->stype[i], &lb, &extent);
	printf("%s %d %d %d %d \n",__FILE__,__LINE__,i,lb,extent);
      }
            }
      */
  }

  return PIO_NOERR;

}


/**
 ** @internal
**  Completes the mapping for the box rearranger
** @endinternal
*/

int compute_counts(const iosystem_desc_t ios, io_desc_t *iodesc, const int maplen, 
		   const int dest_ioproc[], const PIO_Offset dest_ioindex[], MPI_Comm mycomm)
{

  int i;
  int iorank;

  int rank;
  int ntasks;

  MPI_Comm_rank(mycomm, &rank);
  MPI_Comm_size(mycomm, &ntasks);


  MPI_Datatype sr_types[ntasks];
  int send_counts[ntasks];
  int send_displs[ntasks];
  int recv_counts[ntasks];
  int recv_displs[ntasks];
  int *recv_buf=NULL;
  int nrecvs;
  int maxreq = MAX_GATHER_BLOCK_SIZE;
  int ierr;
  int io_comprank;
  int ioindex;
  int tsize;
  int numiotasks;
  PIO_Offset s2rindex[iodesc->ndof];


  
  if(iodesc->rearranger==PIO_REARR_BOX)
    numiotasks = ios.num_iotasks;
  else
    numiotasks=1;

  iodesc->scount = (int *) bget(numiotasks*sizeof(int));
  if(iodesc->scount == NULL){
    piomemerror(ios,numiotasks * sizeof(int), __FILE__,__LINE__);
  }
  for(i=0;i<numiotasks;i++)
    iodesc->scount[i]=0;

  // iodesc->scount is the amount of data sent to each task from the current task
  for(i=0;i<maplen; i++){
    if(dest_ioindex[i] >= 0){
      (iodesc->scount[dest_ioproc[i]])++;
    }
  }

  //  for(i=0;i<ios.num_iotasks;i++)
  //   printf("iodesc->scount = %d\n",iodesc->scount[i]);

  for(i=0;i<ntasks;i++){
    send_counts[i] = 0;
    send_displs[i] = 0;
    recv_counts[i] = 0;
    recv_displs[i] = 0;
    sr_types[i] = MPI_INT;
  }
  for(i=0;i<numiotasks;i++){
    int io_comprank;
    if(iodesc->rearranger==PIO_REARR_SUBSET)
      io_comprank=0;
    else
      io_comprank = ios.ioranks[i];
    send_counts[io_comprank] = 1;
    send_displs[io_comprank] = i*sizeof(int);
  }

  if(ios.ioproc){
    recv_buf = (int *) bget(ntasks * sizeof(int));
    if(recv_buf == NULL){
      piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
    }
    for(i=0;i<ntasks;i++){
      recv_buf[i] = 0;
      recv_counts[i] = 1;
      recv_displs[i] = i*sizeof(int);
    }
  }
  //  for(i=0;i<numiotasks;i++)
  //  printf("%s %d %d\n",__FILE__,__LINE__,iodesc->scount[i]);

  // Share the iodesc->scount from each compute task to all io tasks
  ierr = pio_swapm( iodesc->scount, send_counts, send_displs, sr_types, 
                    recv_buf,  recv_counts, recv_displs, sr_types,
		    mycomm, false, false, maxreq);

  nrecvs = 0;
  if(ios.ioproc){
    //       printf("recv_buf = ");
    for(i=0;i<ntasks; i++){
      //     printf(" %d ",recv_buf[i]);
      if(recv_buf[i] != 0)
	nrecvs++;
    }
    // printf("\n");

    iodesc->rcount = (int *) bget(max(1,nrecvs)*sizeof(int));
    if(iodesc->rcount == NULL){
      piomemerror(ios,max(1,nrecvs) * sizeof(int), __FILE__,__LINE__);
    }
    iodesc->rfrom = (int *)  bget(max(1,nrecvs)*sizeof(int));
    if(iodesc->rfrom == NULL){
      piomemerror(ios,max(1,nrecvs) * sizeof(int), __FILE__,__LINE__);
    }
    for(i=0;i<max(1,nrecvs);i++){
      iodesc->rcount[i]=0;
      iodesc->rfrom[i]=0;
    }

    nrecvs = 0;
    for(i=0;i<ntasks; i++){
      if(recv_buf[i] != 0){
	iodesc->rcount[nrecvs] = recv_buf[i];
	iodesc->rfrom[nrecvs] = i;
	nrecvs++;
      }

    }
    brel(recv_buf);
  }

  iodesc->nrecvs = nrecvs;
  if(iodesc->sindex == NULL && iodesc->ndof>0){
    iodesc->sindex = (PIO_Offset *) bget(iodesc->ndof * sizeof(PIO_Offset));
    if(iodesc->sindex == NULL){
      piomemerror(ios,iodesc->ndof * sizeof(PIO_Offset), __FILE__,__LINE__);
    }
    for(i=0;i<iodesc->ndof;i++)
      iodesc->sindex[i]=0;
  }
  int tempcount[numiotasks];
  int spos[numiotasks];

  spos[0]=0;
  tempcount[0]=0;
  for(i=1;i<numiotasks;i++){
    spos[i] = spos[i-1] + iodesc->scount[i-1];
    tempcount[i]=0;
  }

  for(i=0;i<maplen;i++){
    iorank =dest_ioproc[i]; 
    ioindex = dest_ioindex[i];
    if(iorank > -1){
      // this should be moved to create_box
      if(iodesc->rearranger==PIO_REARR_BOX)
	iodesc->sindex[spos[iorank]+tempcount[iorank]] = i;

      s2rindex[spos[iorank]+tempcount[iorank]] = ioindex;
      (tempcount[iorank])++;
    }
  }
    //    printf("%s %d %d %d %d %d\n",__FILE__,__LINE__,iodesc->llen,iodesc->ndof, maplen,spos[0]+tempcount[0]);

  for(i=0;i<ntasks;i++){
    send_counts[i] = 0;
    send_displs[i]  = 0;
    recv_counts[i] = 0;
    recv_displs[i]   =0;
  }
  MPI_Type_size(MPI_OFFSET, &tsize);

  for(i=0; i<ntasks; i++){
    sr_types[i] = MPI_OFFSET;
  }

  for(i=0;i<numiotasks;i++){
    if(iodesc->rearranger==PIO_REARR_BOX){
      io_comprank = ios.ioranks[i];
    }else{
      io_comprank=0;
    }
    send_counts[io_comprank] = iodesc->scount[i];
    if(send_counts[io_comprank]>0)
      send_displs[io_comprank]  = spos[i]*tsize ;
  }

  if(ios.ioproc){
    int totalrecv=0;
    for(i=0;i<nrecvs;i++){
      recv_counts[iodesc->rfrom[i]] = iodesc->rcount[i];
      totalrecv+=iodesc->rcount[i];
    }
    recv_displs[0] = 0;
    for(i=1;i<nrecvs;i++)
      recv_displs[iodesc->rfrom[i]] = recv_displs[iodesc->rfrom[i-1]]+iodesc->rcount[i-1]*tsize;


    if(totalrecv>0){
      //      printf("%s %d %d %d\n",__FILE__,__LINE__,totalrecv,iodesc->llen);
      totalrecv = iodesc->llen;  // can reduce memory usage here
      iodesc->rindex = (PIO_Offset *) bget(totalrecv*sizeof(PIO_Offset));
      if(iodesc->rindex == NULL){
	piomemerror(ios,totalrecv * sizeof(PIO_Offset), __FILE__,__LINE__);
      }
      for(i=0;i<totalrecv;i++)
	iodesc->rindex[i]=0;
    }
  }
  //   printf("%d rbuf_size %d\n",ios.comp_rank,rbuf_size);


  // s2rindex is the list of indeces on each compute task
  /*        
  printf("%d s2rindex: ", ios.comp_rank);
  for(i=0;i<iodesc->ndof;i++)
    printf("%ld ",s2rindex[i]);
  printf("\n");
  */
  //  printf("%s %d %d %d %d %d %d %d\n",__FILE__,__LINE__,send_counts[0],recv_counts[0],send_displs[0],recv_displs[0],sr_types[0],iodesc->llen);
  ierr = pio_swapm( s2rindex, send_counts, send_displs, sr_types, 
		    iodesc->rindex, recv_counts, recv_displs, sr_types,
  		    mycomm, false, false, 0);
  // printf("%s %d\n",__FILE__,__LINE__);

  //  rindex is an array of the indices of the data to be sent from
  //  this io task to each compute task. 
  /*   
  if(ios.ioproc){
    printf("%d rindex: ",ios.io_rank);
    for(int j=0;j<iodesc->llen;j++)
      printf(" %ld ",iodesc->rindex[j]);
    printf("\n");
  }
  */
  return ierr;

}

/** 
 ** @internal
 ** Moves data from compute tasks to IO tasks.
 ** @endinternal
 **
 */


int rearrange_comp2io(const iosystem_desc_t ios, io_desc_t *iodesc, void *sbuf,
			  void *rbuf, const int nvars)
{
  int ntasks;
  int niotasks;
  int *scount = iodesc->scount;
  int ivar;
  int i, tsize;
  int *sendcounts;
  int *recvcounts;
  int *sdispls;
  int *rdispls;
  MPI_Datatype *sendtypes;
  MPI_Datatype *recvtypes;
  MPI_Comm mycomm;

  //  maxreq = 0;
  //  printf("%s %d %d %d\n",__FILE__,__LINE__,nvars,iodesc->ndof);
#ifdef TIMING
  GPTLstart("PIO:rearrange_comp2io");
#endif
  
  if(iodesc->rearranger == PIO_REARR_BOX){
    mycomm = ios.union_comm;
    niotasks = ios.num_iotasks;
  }else{
    mycomm = iodesc->subset_comm;
    niotasks = 1;
  }  
  MPI_Comm_size(mycomm, &ntasks);

  pioassert(nvars>0,"nvars must be > 0",__FILE__,__LINE__);
  MPI_Type_size(iodesc->basetype, &tsize);

  define_iodesc_datatypes(ios, iodesc);

  sendcounts = (int *) bget(ntasks*sizeof(int));
  if(sendcounts == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  recvcounts = (int *) bget(ntasks*sizeof(int));
  if(recvcounts == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  sdispls = (int *) bget(ntasks*sizeof(int));
  if(sdispls == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  rdispls = (int *) bget(ntasks*sizeof(int));
  if(rdispls == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  sendtypes = (MPI_Datatype *) bget(ntasks*sizeof(MPI_Datatype));
  if(sendtypes == NULL){
    piomemerror(ios,ntasks * sizeof(MPI_Datatype), __FILE__,__LINE__);
  }
  recvtypes = (MPI_Datatype *) bget(ntasks*sizeof(MPI_Datatype));
  if(recvtypes == NULL){
    piomemerror(ios,ntasks * sizeof(MPI_Datatype), __FILE__,__LINE__);
  }
  //printf("%s %d %ld %ld\n",__FILE__,__LINE__,sendtypes,recvtypes);
  for(i=0;i<ntasks;i++){
    sendcounts[i] = 0;
    recvcounts[i] = 0; 
    sdispls[i] = 0; 
    rdispls[i] = 0;
    recvtypes[ i ] = MPI_DATATYPE_NULL;
    sendtypes[ i ] =  MPI_DATATYPE_NULL;
  }


  if(ios.ioproc && iodesc->nrecvs>0){
    //    printf("%d: rindex[%d] %d\n",ios.comp_rank,0,iodesc->rindex[0]);
    for( i=0;i<iodesc->nrecvs;i++){
      if(iodesc->rtype[i] != MPI_DATATYPE_NULL){
	if(iodesc->rearranger==PIO_REARR_SUBSET){
	  recvcounts[ i ] = 1;
	  // The stride here is the length of the collected array (llen)
          MPI_Type_hvector(nvars, 1, (MPI_Aint) iodesc->llen*tsize, iodesc->rtype[i], recvtypes+i);
	  if(recvtypes[i] == MPI_DATATYPE_NULL){
	    piodie("Unexpected NULL MPI DATATYPE",__FILE__,__LINE__);
	  }
          MPI_Type_commit(recvtypes+i);
	  //recvtypes[ i ] = iodesc->rtype[i];
	}else{
	  recvcounts[ iodesc->rfrom[i] ] = 1;
          MPI_Type_hvector(nvars, 1,(MPI_Aint) iodesc->llen*tsize,iodesc->rtype[i], recvtypes+iodesc->rfrom[i]);
	  if(recvtypes[iodesc->rfrom[i]] == MPI_DATATYPE_NULL){
	    piodie("Unexpected NULL MPI DATATYPE",__FILE__,__LINE__);
	  }
          MPI_Type_commit(recvtypes+iodesc->rfrom[i]);
	  // recvtypes[ iodesc->rfrom[i] ] = iodesc->rtype[i];
	  rdispls[ iodesc->rfrom[i] ] = 0;
	  //    printf("%d: rindex[%d] %d\n",ios.comp_rank,0,iodesc->rindex[0]);
	}
      }
      //   printf("%d: rindex[%d] %d\n",ios.comp_rank,i,iodesc->rindex[i]);
      
    }
  }

  for( i=0;i<niotasks; i++){
    int io_comprank = ios.ioranks[i];
    if(iodesc->rearranger==PIO_REARR_SUBSET)
      io_comprank=0;
    
    //    printf("scount[%d]=%d\n",i,scount[i]);
    if(scount[i] > 0 && sbuf != NULL) {
      sendcounts[io_comprank]=1;
      MPI_Type_hvector(nvars, 1, (MPI_Aint) iodesc->ndof*tsize, iodesc->stype[i], sendtypes+io_comprank);
      if(sendtypes[io_comprank] == MPI_DATATYPE_NULL){
	piodie("Unexpected NULL MPI DATATYPE",__FILE__,__LINE__);
      }
      MPI_Type_commit(sendtypes+io_comprank);
    }else{
      sendcounts[io_comprank]=0;
    }
  }      
  //  printf("%s %d %d\n",__FILE__,__LINE__,((int *)sbuf)[5]);
//printf("%s %d %ld %ld %ld %ld\n",__FILE__,__LINE__,sendtypes[0],recvtypes,sbuf,rbuf);
  // Data in sbuf on the compute nodes is sent to rbuf on the ionodes

  //  for(i=0;i<nvars*iodesc->ndof;i++)
  //   printf("%s %d %d %d \n",__FILE__,__LINE__,i,((int *)sbuf)[i]);

  pio_swapm( sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, 
	     mycomm, iodesc->handshake, iodesc->isend, iodesc->max_requests);
  //    if(ios.ioproc)
  //     for(i=0;i<nvars*iodesc->llen;i++)
  //	printf("%s %d %d %d\n",__FILE__,__LINE__,i,((int *)rbuf)[i]);
  //printf("%s %d %ld %ld\n",__FILE__,__LINE__,sendtypes[0],recvtypes);  

  brel(sendcounts);
  brel(recvcounts); 
  brel(sdispls);
  brel(rdispls);
//printf("%s %d %ld %ld\n",__FILE__,__LINE__,sendtypes[0],recvtypes);
  for( i=0;i<ntasks; i++){
    if(sendtypes[i] != MPI_DATATYPE_NULL){
      //printf("%s %d %d %ld\n",__FILE__,__LINE__,i,sendtypes[i]);
      MPI_Type_free(sendtypes+i);
    }
    if(recvtypes[i] != MPI_DATATYPE_NULL){     
      //printf("%s %d %d\n",__FILE__,__LINE__,i);    
      MPI_Type_free(recvtypes+i);
    }
  }
  
  brel(sendtypes);
  brel(recvtypes);

#ifdef TIMING
  GPTLstop("PIO:rearrange_comp2io");
#endif

  return PIO_NOERR;
}
/** 
 ** @internal
 ** Moves data from compute tasks to IO tasks.
 ** @endinternal
 **
 */

int rearrange_io2comp(const iosystem_desc_t ios, io_desc_t *iodesc, void *sbuf, void *rbuf)
{
  

  //  int maxreq = 0;
  MPI_Comm mycomm;

  int ntasks ;
  int niotasks;
  int *scount = iodesc->scount;

  int *sendcounts;
  int *recvcounts;
  int *sdispls;
  int *rdispls;
  MPI_Datatype *sendtypes;
  MPI_Datatype *recvtypes;

  int i, tsize;
#ifdef TIMING
  GPTLstart("PIO:rearrange_io2comp");
#endif
  if(iodesc->rearranger==PIO_REARR_BOX){
    mycomm = ios.union_comm;
    niotasks = ios.num_iotasks;
  }else{
    mycomm = iodesc->subset_comm;
    niotasks=1;
  }
  MPI_Comm_size(mycomm, &ntasks);

  define_iodesc_datatypes(ios, iodesc);

  sendcounts = (int *) bget(ntasks*sizeof(int));
  if(sendcounts == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  recvcounts = (int *) bget(ntasks*sizeof(int));
  if(recvcounts == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  sdispls = (int *) bget(ntasks*sizeof(int));
  if(sdispls == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  rdispls = (int *) bget(ntasks*sizeof(int));
  if(rdispls == NULL){
    piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
  }
  sendtypes = (MPI_Datatype *) bget(ntasks*sizeof(MPI_Datatype));
  if(sendtypes == NULL){
    piomemerror(ios,ntasks * sizeof(MPI_Datatype), __FILE__,__LINE__);
  }
  recvtypes = (MPI_Datatype *) bget(ntasks*sizeof(MPI_Datatype));
  if(recvtypes == NULL){
    piomemerror(ios,ntasks * sizeof(MPI_Datatype), __FILE__,__LINE__);
  }

  for( i=0;i< ntasks;i++){
    sendcounts[i] = 0;
    recvcounts[i] = 0;
    sdispls[i] = 0;
    rdispls[i] = 0;
    sendtypes[ i ] = MPI_DATATYPE_NULL;
    recvtypes[ i ] = MPI_DATATYPE_NULL;
  }
  if(ios.ioproc){
    for( i=0;i< iodesc->nrecvs;i++){
      if(iodesc->rtype[i] != MPI_DATATYPE_NULL){
	if(iodesc->rearranger==PIO_REARR_SUBSET){
	  if(sbuf != NULL){
	    sendcounts[ i ] = 1;
	    sendtypes[ i ] = iodesc->rtype[i];
	  }
	}else{
	  sendcounts[ iodesc->rfrom[i] ] = 1;
	  sendtypes[ iodesc->rfrom[i] ] = iodesc->rtype[i];
	}
      }
    }
  }

  for( i=0;i<niotasks; i++){
    int io_comprank = ios.ioranks[i];
    if(iodesc->rearranger==PIO_REARR_SUBSET)
      io_comprank=0;
    if(scount[i] > 0 && iodesc->stype[i] != MPI_DATATYPE_NULL) {
      recvcounts[io_comprank]=1;
      recvtypes[io_comprank]=iodesc->stype[i];
    }
  } 
  
  //
  // Data in sbuf on the ionodes is sent to rbuf on the compute nodes
  //
  #if DEBUG
  printf("%s %d \n",__FILE__,__LINE__);
  #endif
  pio_swapm( sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, 
	     mycomm, iodesc->handshake, iodesc->isend, iodesc->max_requests);
  #if DEBUG
  printf("%s %d \n",__FILE__,__LINE__);
  #endif
  brel(sendcounts);
  brel(recvcounts); 
  brel(sdispls);
  brel(rdispls);
  brel(sendtypes);
  brel(recvtypes);

#ifdef TIMING
  GPTLstop("PIO:rearrange_io2comp");
#endif
 
  return PIO_NOERR;

}

void determine_fill(iosystem_desc_t ios, io_desc_t *iodesc, const int gsize[], const PIO_Offset compmap[])
{
  PIO_Offset totalllen=0;
  PIO_Offset totalgridsize=1;
  int i;
  for( i=0;i<iodesc->ndims;i++){
    totalgridsize *= gsize[i];
  }
  if(iodesc->rearranger==PIO_REARR_SUBSET){
    totalllen=iodesc->llen;
  }else{
    for(i=0;i<iodesc->ndof;i++){
      if(compmap[i]>0){
	totalllen++;
      }
    }
  }
  //  printf("%s %d %ld %ld\n",__FILE__,__LINE__,totalllen,totalgridsize);

  MPI_Allreduce(MPI_IN_PLACE, &totalllen, 1, PIO_OFFSET, MPI_SUM, ios.union_comm);

  
  if(totalllen < totalgridsize){
    //if(ios.iomaster) printf("%s %d %ld %ld\n",__FILE__,__LINE__,totalllen,totalgridsize);
    iodesc->needsfill = true;
  }else{
    iodesc->needsfill = false;
    iodesc->gsize = NULL;
  }
  //  TURN OFF FILL for timing test
  //  iodesc->needsfill=false;

  if(iodesc->needsfill){
    iodesc->gsize = (PIO_Offset *) bget(iodesc->ndims * sizeof(PIO_Offset));
    if(iodesc->gsize==NULL){
      piomemerror(ios,iodesc->ndims * sizeof(MPI_Datatype), __FILE__,__LINE__);
    }
    for (i=0;i<iodesc->ndims;i++){
      iodesc->gsize[i] = gsize[i];
    }
  }


}


void iodesc_dump(io_desc_t *iodesc)
{
  printf("ioid= %d\n",iodesc->ioid);
  printf("async_id= %d\n",iodesc->async_id);
  printf("nrecvs= %d\n",iodesc->nrecvs);
  printf("ndof= %d\n",iodesc->ndof);
  printf("ndims= %d\n",iodesc->ndims);
  printf("num_aiotasks= %d\n",iodesc->num_aiotasks);
  printf("rearranger= %d\n",iodesc->rearranger);
  printf("maxregions= %d\n",iodesc->maxregions);
  printf("needsfill= %d\n",(int) iodesc->needsfill);

  printf("llen= %ld\n",iodesc->llen);
  printf("maxiobuflen= %d\n",iodesc->maxiobuflen);
  printf("gsize=");
  for(int i=0; i< iodesc->ndims; i++){
      printf(" %d",iodesc->ioid);
  }
  printf("\n");

  printf("rindex= ");
  for(int j=0;j<iodesc->llen;j++)
    printf(" %ld ",iodesc->rindex[j]);
  printf("\n");
  


}

/** 
 ** @internal
 ** The box rearranger computes a mapping between IO tasks and compute tasks such that the data
 ** on io tasks can be written with a single call to the underlying netcdf library.   This 
 ** may involve an all to all rearrangement in the mapping, but should minimize data movement in 
 ** lower level libraries
 ** @endinternal
 **
 */

int box_rearrange_create(const iosystem_desc_t ios,const int maplen, const PIO_Offset compmap[], const int gsize[],
			 const int ndims, io_desc_t *iodesc)
{
  int ierr=PIO_NOERR;
  int nprocs = ios.num_comptasks;
  int nioprocs = ios.num_iotasks;
  PIO_Offset gstride[ndims];
  PIO_Offset iomap;
  PIO_Offset start[ndims], count[ndims];
  int  tsize, i, j, k, llen;
  MPI_Datatype dtype;
  int dest_ioproc[maplen];
  PIO_Offset dest_ioindex[maplen];
  int sndlths[nprocs]; 
  int sdispls[nprocs];
  int recvlths[nprocs];
  int rdispls[nprocs];
  MPI_Datatype dtypes[nprocs];
  PIO_Offset iomaplen[nioprocs];
  int maxreq = MAX_GATHER_BLOCK_SIZE;

  iodesc->rearranger = PIO_REARR_BOX;

  iodesc->ndof = maplen;
  gstride[ndims-1]=1;
  for(int i=ndims-2;i>=0; i--)
    gstride[i]=gstride[i+1]*gsize[i+1];

  MPI_Type_size(MPI_OFFSET, &tsize);

  for(i=0; i< maplen; i++){
    dest_ioproc[i] = -1;
    dest_ioindex[i] = -1;
  }
  for(i=0;i<nprocs;i++){
    sndlths[i] = 0;
    sdispls[i] = 0;
    recvlths[i] = 0;
    rdispls[i] = 0;
    dtypes[i] = MPI_OFFSET;
  }
  iodesc->llen=0;
  if(ios.ioproc){
    for( i=0;i<nprocs;i++){
      sndlths[ i ] = 1;
    }
    // llen here is the number that will be read on each io task
    iodesc->llen=1;
    for(i=0;i<ndims;i++)
      iodesc->llen *= iodesc->firstregion->count[i];
  }
  determine_fill(ios, iodesc, gsize, compmap);

  /* 
  if(ios.ioproc){
    for(i=0; i<ndims; i++){
      printf("%d %d %d ",i,iodesc->firstregion->start[i],iodesc->firstregion->count[i]);
    }
    printf("\n%s %d\n",__FILE__,__LINE__);
  }
  */

  for( i=0;i<nioprocs; i++){
    int io_comprank = ios.ioranks[i];
    recvlths[ io_comprank ] = 1;
    rdispls[ io_comprank ] = i*tsize;
  }      

  //  The length of each iomap
  //  iomaplen = calloc(nioprocs, sizeof(PIO_Offset));
  pio_swapm(&(iodesc->llen), sndlths, sdispls, dtypes,
	    iomaplen, recvlths, rdispls, dtypes, 	
	    ios.union_comm, false, false, maxreq);




  /*
  printf("%s %d %d\n",__FILE__,__LINE__,nioprocs);
  for(i=0; i<nioprocs; i++){
    printf("%ld ",iomaplen[i]);
  }
  printf("\n");
  */
  for(i=0; i<nioprocs; i++){

    if(iomaplen[i]>0){
      int io_comprank = ios.ioranks[i];
      for( j=0; j<nprocs ; j++){
	sndlths[ j ] = 0;
	sdispls[ j ] = 0;
	rdispls[ j ] = 0;
	recvlths[ j ] = 0;
	if(ios.union_rank == io_comprank)
	  sndlths[ j ] = ndims;
      }
      recvlths[ io_comprank ] = ndims;
      
      // The count from iotask i is sent to all compute tasks
      
      pio_swapm(iodesc->firstregion->count,  sndlths, sdispls, dtypes,
		count, recvlths, rdispls, dtypes, 
		ios.union_comm, false, false, maxreq);
      
      // The start from iotask i is sent to all compute tasks
      pio_swapm(iodesc->firstregion->start,  sndlths, sdispls, dtypes,
		start, recvlths, rdispls, dtypes, 
		ios.union_comm, false, false, maxreq);

      for(k=0;k<maplen;k++){
	PIO_Offset gcoord[ndims], lcoord[ndims];
	bool found=true;
	// The compmap array is 1 based but calculations are 0 based
	idx_to_dim_list(ndims, gsize, compmap[k]-1, gcoord);


	for(j=0;j<ndims;j++){
	  if(gcoord[j] >= start[j] && gcoord[j] < start[j]+count[j]){
	    lcoord[j] = gcoord[j] - start[j];
	  }else{
	    found = false;
	    break;
	  }
	}
	if(found){
	  dest_ioindex[k] = coord_to_lindex(ndims, lcoord, count);
	  dest_ioproc[k] = i;
	}
      }
    }
  }

  for(k=0; k<maplen; k++){
    //       printf("%s %d %d %d\n",__FILE__,__LINE__,dest_ioproc[k],dest_ioindex[k]);

    if(dest_ioproc[k] < 0 && compmap[k]>0){
      fprintf(stderr,"No destination found for compmap[%d] = %ld\n",k,compmap[k]);
      piodie(" ",__FILE__,__LINE__);
    }
  }

  compute_counts(ios, iodesc, maplen, dest_ioproc, dest_ioindex, ios.union_comm);

  if(ios.ioproc){
    compute_maxIObuffersize(ios.io_comm, iodesc);
  }
  compute_maxaggregate_bytes(ios, iodesc);


#ifdef DEBUG  
  iodesc_dump(iodesc);
#endif
  return PIO_NOERR;
}
/** 
 ** @internal
 ** compare offsets is used by the sort in the subset rearrange
 ** @endinternal
 */
int compare_offsets(const void *a,const void *b) 
{
  mapsort *x = (mapsort *) a;
  mapsort *y = (mapsort *) b;
  return (int) (x->iomap - y->iomap);
}    

/**
 ** @internal
 ** Each region is a block of output which can be represented in a single call to the underlying 
 ** netcdf library.  This can be as small as a single data point, but we hope we've aggragated better than that. 
 ** @endinternal
 */
void get_start_and_count_regions(const int ndims, const int gdims[],const int maplen, 
				 const PIO_Offset map[], int *maxregions, io_region *firstregion)
{
  int i;
  int nmaplen;
  int regionlen;
  io_region *region, *prevregion;

  nmaplen = 0;
  region = firstregion;
  if(map != NULL){
    while(map[nmaplen++]<=0);
    nmaplen--;
  }
  region->loffset=nmaplen;

  *maxregions = 1;
  prevregion=NULL;

  while(nmaplen < maplen){
    // Here we find the largest region from the current offset into the iomap
    // regionlen is the size of that region and we step to that point in the map array
    // until we reach the end 
    for(i=0;i<ndims;i++){
      region->count[i]=1;
    }

    regionlen = find_region(ndims, gdims, maplen-nmaplen, 
				  map+nmaplen, region->start, region->count);
    //    printf("%s %d %d %d\n",__FILE__,__LINE__,region->start[0],region->count[0]);
    pioassert(region->start[0]>=0,"failed to find region",__FILE__,__LINE__);
    
    nmaplen = nmaplen+regionlen;
    if(region->next==NULL && nmaplen< maplen){
      region->next = alloc_region(ndims);
      // The offset into the local array buffer is the sum of the sizes of all of the previous regions (loffset) 
      region=region->next;
      region->loffset = nmaplen;
      // The calls to the io library are collective and so we must have the same number of regions on each
      // io task maxregions will be the total number of regions on this task 
      (*maxregions)++;
    }
  }


}

/** 
 ** @internal
 ** The subset rearranger needs a mapping from compute tasks to IO task, the only requirement is 
 ** that each compute task map to one and only one IO task.   This mapping groups by mpi task id
 ** others are possible and may be better for certain decompositions
 ** @endinternal
 **
 */
void default_subset_partition(const iosystem_desc_t ios, io_desc_t *iodesc)
{
  int taskratio = ios.num_comptasks/ios.num_iotasks;
  int color;
  int key;

  /* Create a new comm for each subset group with the io task in rank 0 and
     only 1 io task per group */

  if(ios.ioproc){
    key = 0;
    color= ios.io_rank;
  }else{
    key=max(1,ios.comp_rank%taskratio+1);
    color = min(ios.num_iotasks-1,ios.comp_rank/taskratio);
  }

  MPI_Comm_split(ios.comp_comm, color, key, &(iodesc->subset_comm));

}

/** 
 ** @internal
 ** The subset rearranger computes a mapping between IO tasks and compute tasks such that each compute
 ** task communicates with one and only one IO task.  
 ** @endinternal
 **
 */

int subset_rearrange_create(const iosystem_desc_t ios,const int maplen, PIO_Offset compmap[], 
			    const int gsize[], const int ndims, io_desc_t *iodesc)
{

  int taskratio;

  int i, j, jlast;
  bool hs=false;
  bool isend=false;
  PIO_Offset *iomap=NULL;
  int ierr = PIO_NOERR;
  mapsort *map=NULL;
  PIO_Offset totalgridsize;
  PIO_Offset *srcindex=NULL;
  PIO_Offset *myfillgrid = NULL;
  int maxregions;  
  int maxreq = MAX_GATHER_BLOCK_SIZE;
  int rank, ntasks, rcnt;
  size_t pio_offset_size=sizeof(PIO_Offset);
  
  
  tmpioproc = ios.io_rank;


  /* subset partitions each have exactly 1 io task which is task 0 of that subset_comm */ 
  /* TODO: introduce a mechanism for users to define partitions */
  default_subset_partition(ios, iodesc);
  iodesc->rearranger = PIO_REARR_SUBSET;
  
  MPI_Comm_rank(iodesc->subset_comm, &rank);
  MPI_Comm_size(iodesc->subset_comm, &ntasks);

  if(ios.ioproc){
    pioassert(rank==0,"Bad io rank in subset create",__FILE__,__LINE__);
  }else{
    pioassert(rank>0 && rank<ntasks,"Bad comp rank in subset create",__FILE__,__LINE__);
  }
  rcnt = 0;
  iodesc->ndof = maplen;
  if(ios.ioproc){
    iodesc->rcount = (int *) bget(ntasks *sizeof(int));
    if(iodesc->rcount == NULL){
      piomemerror(ios,ntasks * sizeof(int), __FILE__,__LINE__);
    }
    rcnt = 1;
  } 
  iodesc->scount = (int *) bget(sizeof(int));
  if(iodesc->scount == NULL){
    piomemerror(ios,sizeof(int), __FILE__,__LINE__);
  }
  iodesc->scount[0]=0;
  totalgridsize=1;
  for( i=0;i<ndims;i++){
    totalgridsize *= gsize[i];
  }

  for(i=0;i<maplen;i++){
    //  turns out this can be allowed in some cases 
    //    pioassert(compmap[i]>=0 && compmap[i]<=totalgridsize, "Compmap value out of bounds",__FILE__,__LINE__);
    if(compmap[i]>0){
      (iodesc->scount[0])++;
    }
  }
  if(iodesc->scount[0]>0){
    iodesc->sindex = (PIO_Offset *) bget(iodesc->scount[0]*pio_offset_size); 
    if(iodesc->sindex == NULL){
      piomemerror(ios,iodesc->scount[0]*pio_offset_size, __FILE__,__LINE__);
    }
    for(i=0;i<iodesc->scount[0];i++)
      iodesc->sindex[i]=0;
  }

  j=0;
  for(i=0;i<maplen;i++){
    if(compmap[i]>0){
      iodesc->sindex[j++]=i;

    }
  }

  // Pass the reduced maplen (without holes) from each compute task to its associated IO task
  //  printf("%s %d %ld\n",__FILE__,__LINE__,iodesc->scount);
  
  pio_fc_gather( (void *) iodesc->scount, 1, MPI_INT,
		 (void *) iodesc->rcount, rcnt, MPI_INT, 
		 0, iodesc->subset_comm, maxreq);

  iodesc->llen = 0;

  int rdispls[ntasks];
  int recvlths[ntasks];

  if(ios.ioproc){
    for(i=0;i<ntasks;i++){
      iodesc->llen+=iodesc->rcount[i];
      rdispls[i]=0;
      recvlths[i]= iodesc->rcount[ i ];
      if(i>0)
	rdispls[i] = rdispls[i-1]+ iodesc->rcount[ i-1 ];
    }
    //    printf("%s %d %ld %d %d\n",__FILE__,__LINE__,iodesc,iodesc->llen,maplen);
      
    if(iodesc->llen>0){
      srcindex = (PIO_Offset *) bget(iodesc->llen*pio_offset_size);
      if(srcindex == NULL){
	piomemerror(ios,iodesc->llen * pio_offset_size, __FILE__,__LINE__);
      }
      for(i=0;i < iodesc->llen; i++)
	srcindex[i]=0;
    }
  }else{
    for(i=0;i<ntasks;i++){
      recvlths[i]=0;
      rdispls[i]=0;
    }
  }
  determine_fill(ios, iodesc, gsize, compmap);

  // Pass the sindex from each compute task to its associated IO task

  pio_fc_gatherv((void *) iodesc->sindex, iodesc->scount[0], PIO_OFFSET,
		 (void *) srcindex, recvlths, rdispls, PIO_OFFSET,  
		 0, iodesc->subset_comm, maxreq);

  if(ios.ioproc && iodesc->llen>0){
    map = (mapsort *) bget(iodesc->llen * sizeof(mapsort));    
    if(map == NULL){
      piomemerror(ios,iodesc->llen * sizeof(mapsort), __FILE__,__LINE__);
    }
    iomap = (PIO_Offset *) bget(iodesc->llen*pio_offset_size);
    if(iomap == NULL){
      piomemerror(ios,iodesc->llen * pio_offset_size, __FILE__,__LINE__);
    }
    for(i=0;i<iodesc->llen;i++)
      iomap[i]=0;
  }

  // Now pass the compmap, skipping the holes

  PIO_Offset *shrtmap;
  if(maplen>iodesc->scount[0] && iodesc->scount[0]>0){
    shrtmap = (PIO_Offset *) bget(iodesc->scount[0]*pio_offset_size);
    if(shrtmap == NULL){
      piomemerror(ios,iodesc->scount[0] * pio_offset_size, __FILE__,__LINE__);
    }
    for(i=0;i<iodesc->scount[0];i++)
      shrtmap[i]=0;

    j=0;
    for(i=0;i<maplen;i++)
      if(compmap[i]>0)
	shrtmap[j++]=compmap[i];
  }else{
    shrtmap = compmap;
  }

  pio_fc_gatherv((void *) shrtmap, iodesc->scount[0], PIO_OFFSET,
		 (void *) iomap, recvlths, rdispls, PIO_OFFSET,  
		 0, iodesc->subset_comm, maxreq);


  if(shrtmap != compmap)
    brel(shrtmap);

  if(ios.ioproc && iodesc->llen>0){
    int pos=0;
    int k=0;
    mapsort *mptr;
    for(i=0;i<ntasks;i++){
      for(j=0;j<iodesc->rcount[i];j++){
	mptr = map+k;
	mptr->rfrom = i;
	mptr->soffset = srcindex[pos+j];
	mptr->iomap = iomap[pos+j];
	k++;
      }
      pos += iodesc->rcount[i];
    }
    // sort the mapping, this will transpose the data into IO order        
    qsort(map, iodesc->llen, sizeof(mapsort), compare_offsets); 

    iodesc->rindex = (PIO_Offset *) bget(iodesc->llen*pio_offset_size);
    if(iodesc->rindex == NULL){
      piomemerror(ios,iodesc->llen * pio_offset_size, __FILE__,__LINE__);
    }
    iodesc->rfrom = (int *) bget(iodesc->llen*sizeof(int));
    if(iodesc->rfrom == NULL){
      piomemerror(ios,iodesc->llen * sizeof(int), __FILE__,__LINE__);
    }
    for(i=0;i<iodesc->llen;i++){
      iodesc->rindex[i]=0;
      iodesc->rfrom[i]=0;
    }
  }
  int cnt[ntasks];
  int sndlths[ntasks];
  int sdispls[ntasks];
  MPI_Datatype dtypes[ntasks];
  for(i=0;i<ntasks;i++){
    cnt[i]=rdispls[i];

    /* offsets to swapm are in bytes */
    //    rdispls[i]*=pio_offset_size;
    sdispls[i]=0;
    sndlths[i]=0;
    dtypes[i]=PIO_OFFSET;
  }
  sndlths[0]=iodesc->scount[0];
  mapsort *mptr;
  for(i=0;i<iodesc->llen;i++){
    mptr = map+i;
    iodesc->rfrom[i] = mptr->rfrom;
    iodesc->rindex[i]=i;
    iomap[i] = mptr->iomap;
    srcindex[ (cnt[iodesc->rfrom[i]])++   ]=mptr->soffset;
  }

  if(ios.ioproc && iodesc->needsfill){
    /* we need the list of offsets which are not in the union of iomap */
    PIO_Offset thisgridsize[ios.num_iotasks];
    PIO_Offset thisgridmin[ios.num_iotasks], thisgridmax[ios.num_iotasks];
    int nio;
    PIO_Offset *myusegrid = NULL;
    int gcnt[ios.num_iotasks];
    int displs[ios.num_iotasks];

    thisgridmin[0] = 1;
    thisgridsize[0] =  totalgridsize/ios.num_iotasks;
    thisgridmax[0] = thisgridsize[0];
    int xtra = totalgridsize - thisgridsize[0]*ios.num_iotasks;
  
    for(nio=0;nio<ios.num_iotasks;nio++){
      int cnt=0;
      int imin=0;
      if(nio>0){
	thisgridsize[nio] =  totalgridsize/ios.num_iotasks;
	if(nio>=ios.num_iotasks-xtra) thisgridsize[nio]++;
	thisgridmin[nio]=thisgridmax[nio-1]+1;
	thisgridmax[nio]= thisgridmin[nio] + thisgridsize[nio]-1;
      }
      for(int i=0;i<iodesc->llen;i++){
	//	printf("%s %d %d %d %ld %d %d\n",__FILE__,__LINE__,i,nio,iomap[i],thisgridmin[nio],thisgridmax[nio]);
	if(iomap[i]>=thisgridmin[nio] && iomap[i]<=thisgridmax[nio]){
	  cnt++;
	  if(cnt==1) imin=i;
	}
      }
      //      printf("%s %d %d %d %d\n",__FILE__,__LINE__,nio,cnt,imin);
      MPI_Gather(&cnt, 1, MPI_INT, gcnt, 1, MPI_INT, nio, ios.io_comm);
      if(nio==ios.io_rank){
	displs[0]=0;
	//	printf("%s %d %d %d %d %d\n",__FILE__,__LINE__,nio,0,gcnt[0],displs[0]);
	for(i=1;i<ios.num_iotasks;i++){
	  displs[i] = displs[i-1]+gcnt[i-1];
	  //	  printf("%s %d %d %d %d %d\n",__FILE__,__LINE__,nio,i,gcnt[i],displs[i]);
	} 
	myusegrid = (PIO_Offset *) bget(thisgridsize[nio]*sizeof(PIO_Offset));
	for(i=0;i<thisgridsize[nio];i++){
	  myusegrid[i]=-1;
	}
      }
      pio_fc_gatherv((void *) (iomap+imin), cnt, PIO_OFFSET,
		 (void *) myusegrid, gcnt, displs, PIO_OFFSET,  
		 nio, ios.io_comm, maxreq);
    }
    PIO_Offset grid[thisgridsize[ios.io_rank]];
    for(i=0;i<thisgridsize[ios.io_rank];i++){
      grid[i]=0;
    }
    int cnt=0;
    for(i=0;i<thisgridsize[ios.io_rank];i++){
      int j = myusegrid[i] - thisgridmin[ios.io_rank];
      //      printf("%s %d %d %d %d %d\n",__FILE__,__LINE__,i,j,thisgridmin[ios.io_rank],myusegrid[i]);

      pioassert(j<thisgridsize[ios.io_rank] ,"out of bounds array index",__FILE__,__LINE__);
      if(j>=0){
	grid[j]=1;
	cnt++;
      }
    }
    if(myusegrid!=NULL){
      brel(myusegrid);
    }
    iodesc->holegridsize=thisgridsize[ios.io_rank]-cnt;
    if(iodesc->holegridsize>0){
      myfillgrid = (PIO_Offset *) bget(iodesc->holegridsize * sizeof(PIO_Offset));
    }
    for(i=0;i<iodesc->holegridsize;i++){
      myfillgrid[i]=-1;
    }
	
    j=0;
    for(i=0;i<thisgridsize[ios.io_rank];i++){
      if(grid[i]==0 ){
	if(myfillgrid[j] == -1){
	  myfillgrid[j++]=thisgridmin[ios.io_rank]+i;
	}else{
	  piodie("something wrong",__FILE__,__LINE__);
	  // printf("%s %d %d %d %d\n",__FILE__,__LINE__,i,j,myfillgrid[j]);
	}
      }
    } 
    maxregions = 0;	
    iodesc->maxfillregions=0;
    if(myfillgrid!=NULL){
      iodesc->fillregion = alloc_region(iodesc->ndims);
      get_start_and_count_regions(iodesc->ndims,gsize,iodesc->holegridsize, myfillgrid,&(iodesc->maxfillregions),
				  iodesc->fillregion);
      brel(myfillgrid);
      maxregions = iodesc->maxfillregions;
    }
    MPI_Allreduce(MPI_IN_PLACE,&maxregions,1, MPI_INT, MPI_MAX, ios.io_comm);
    iodesc->maxfillregions = maxregions;
  }

  CheckMPIReturn(MPI_Scatterv((void *) srcindex, recvlths, rdispls, PIO_OFFSET, 
			      (void *) iodesc->sindex, iodesc->scount[0],  PIO_OFFSET,
			      0, iodesc->subset_comm),__FILE__,__LINE__);



  if(ios.ioproc){

    /*
    printf("%s %d %d\n",__FILE__,__LINE__,ios.io_rank);
    for(i=0;i<iodesc->llen;i++)
      printf("%d ",iomap[i]);
    printf("\n");
    */
    iodesc->maxregions = 0;
    get_start_and_count_regions(iodesc->ndims,gsize,iodesc->llen, iomap,&(iodesc->maxregions),
				iodesc->firstregion);
    maxregions = iodesc->maxregions;
    MPI_Allreduce(MPI_IN_PLACE,&maxregions,1, MPI_INT, MPI_MAX, ios.io_comm);
    iodesc->maxregions = maxregions; 
    if(iomap != NULL)
      brel(iomap);
    
    
    if(map != NULL)
      brel(map);

    if(srcindex != NULL)
      brel(srcindex);
    
    compute_maxIObuffersize(ios.io_comm, iodesc);
    
    iodesc->nrecvs=ntasks;
#ifdef DEBUG
    iodesc_dump(iodesc);
#endif  

  }
  /* using maxiobuflen compute the maximum number of vars of this type that the io
     task buffer can handle */

  compute_maxaggregate_bytes(ios, iodesc);


  return ierr;


}
  
    

void performance_tune_rearranger(iosystem_desc_t ios, io_desc_t *iodesc)
{
  bool handshake;
  bool isend;
  int nreqs;
  int maxreqs;
  int nprocs;
  MPI_Comm mycomm;
#ifdef TIMING
#ifdef PERFTUNE  
  double *wall, usr[2], sys[2];
  void *cbuf, *ibuf;
  int tsize;
  int myrank;

  MPI_Type_size(iodesc->basetype, &tsize);
  cbuf = NULL;
  ibuf = NULL;
  if(iodesc->ndof>0){
    cbuf = bget( iodesc->ndof * tsize );
  }
  if(iodesc->llen>0){
    ibuf = bget( iodesc->llen * tsize );
  }

  if(iodesc->rearranger == PIO_REARR_BOX){
    mycomm = ios.union_comm;
  }else{
    mycomm = iodesc->subset_comm;
  }  

  MPI_Comm_size(mycomm, &nprocs);
  MPI_Comm_rank(mycomm, &myrank);

  int log2 = log(nprocs) / log(2) + 1;
  wall = bget(2*4*log2*sizeof(double));
  double mintime;
  int k=0;

  MPI_Barrier(mycomm);
  GPTLstamp( &wall[0], &usr[0], &sys[0]);
  rearrange_comp2io(ios, iodesc, cbuf, ibuf, 1);
  rearrange_io2comp(ios, iodesc, ibuf, cbuf);
  GPTLstamp( &wall[1], &usr[1], &sys[1]);
  mintime = wall[1]-wall[0];
  MPI_Allreduce(MPI_IN_PLACE, &mintime, 1, MPI_DOUBLE, MPI_MAX, mycomm);

  handshake = iodesc->handshake;
  isend = iodesc->isend;
  maxreqs = iodesc->max_requests;

  for(int i=0; i<2; i++){
    if(i==0){
      iodesc->handshake = false;
    }else{
      iodesc->handshake = true;
    }
    for(int j=0; j<2; j++){
      if(j==0){
	iodesc->isend = false;
      }else{
	iodesc->isend = true;
      }
      iodesc->max_requests = 0;

      for(nreqs=nprocs;nreqs>=2;nreqs/=2){
	iodesc->max_requests = nreqs;

	//	if(myrank==0){
	//	  printf("%s %d %d %d %d %f\n",__FILE__,__LINE__,nreqs,handshake,isend,mintime);
	//	}
	MPI_Barrier(mycomm);
	GPTLstamp( wall, usr, sys);
	rearrange_comp2io(ios, iodesc, cbuf, ibuf, 1);
	rearrange_io2comp(ios, iodesc, ibuf, cbuf);
	GPTLstamp( wall+1, usr, sys);
	wall[1]-=wall[0];
	MPI_Allreduce(MPI_IN_PLACE, wall+1,1,MPI_DOUBLE, MPI_MAX, mycomm);

	if(wall[1] < mintime*0.95){
	  handshake = iodesc->handshake;
	  isend = iodesc->isend;
	  maxreqs = nreqs;
	  mintime = wall[1];
	}else if(wall[1]> (mintime*1.05)){
	  exit;
	}
      }
    }
  }
  

  iodesc->handshake = handshake;
  iodesc->isend = isend;
  iodesc->max_requests = maxreqs;
  if(myrank==0){
    printf("spmd optimization: maxreqs: %d handshake:%d isend:%d mintime=%f\n",maxreqs,handshake,isend,mintime);
  }
  brel(wall);
  brel(cbuf);
  brel(ibuf);
#endif
#endif

}  
