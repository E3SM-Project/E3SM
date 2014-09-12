///
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
#include <limits.h>
#define DEF_P2P_MAXREQ 64

int tmpioproc=-1;

/// Convert a global array index to a global coordinate value
void gindex_to_coord(const int ndims, const PIO_Offset gindex, const PIO_Offset gstride[], PIO_Offset *gcoord)
{
  PIO_Offset tempindex;
  int i;

  tempindex = gindex;
  for(i=0;i<ndims-1;i++){
    gcoord[i] = tempindex/gstride[i];
    tempindex -= gcoord[i]*gstride[i];
  }
  gcoord[ndims-1] = tempindex;

}

/// Convert a global coordinate value into a local array index
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

void compute_maxIObuffersize(MPI_Comm io_comm, io_desc_t *iodesc)
{
  PIO_Offset iosize, totiosize;
  int i;
  io_region *region;

  //  compute the max io buffer size, for conveneance it is the combined size of all regions
  totiosize=0;
  region = iodesc->firstregion;
  while(region != NULL){
    iosize = 0;
    if(region->count[0]>0)
      iosize=1;
      for(i=0;i<iodesc->ndims;i++)
	iosize*=region->count[i];
      totiosize+=iosize;
    region = region->next;
  }
  // Share the max io buffer size with all io tasks
#ifndef _MPISERIAL
  CheckMPIReturn(MPI_Allreduce(MPI_IN_PLACE, &totiosize, 1, MPI_OFFSET, MPI_MAX, io_comm),__FILE__,__LINE__);
#endif
  
  iodesc->maxiobuflen = totiosize;
}


/// Expand a region with the given stride. Given an initial region size,
/// this simply checks to see whether the next entries in the map are ahead
/// by the given stride, and then the ones after that. Once max_size is
/// reached or the next entries fail to match, it returns how far it
/// expanded.
int expand_region(const int maplen, const PIO_Offset map[], const int region_size,
                  const int stride, const int max_size)
{
  int i, j, k, expansion;
  int can_expand;
  // Precondition: maplen >= region_size (thus loop runs at least once).

  can_expand = 1;
  
  for (i = 1; i <= max_size; ++i) {
    expansion = i;
    k=i*region_size;
    for (j = 0; j < region_size; ++j) {
      k++;
      pioassert(k<maplen,"out of region bounds",__FILE__,__LINE__);
      if (map[k] != map[j] + i*stride) {
        can_expand = 0;
        break;
      }
    }
    if (!can_expand) break;
  }
  //  printf("expansion %d\n",expansion);
  return expansion;
}

/// Set start and count so that they describe the first region in map.
int find_first_region(const int ndims, const int gdims[],
		      const int maplen, const PIO_Offset map[],
		      PIO_Offset start[], PIO_Offset count[])
{
  int i, region_size, max_size;
  PIO_Offset stride[ndims];
  // Preconditions (which might be useful to check/assert):
  //   ndims is > 0
  //   maplen is > 0
  //   all elements of map are inside the bounds specified by gdims

  stride[ndims-1] = 1;
  for(i=ndims-2;i>=0; --i)
    stride[i] = stride[i+1]*gdims[i+1];

  gindex_to_coord(ndims, map[0], stride, start);

  region_size = 1;
  // For each dimension, figure out how far we can expand in that dimension
  // while staying contiguous in the input array.
  //
  // To avoid something complicated involving recursion, calculate
  // the stride necessary to advance in a given dimension, and feed it into
  // the 1D expand_region function.
  for (i = ndims-1; i >= 0; --i) {
    // Can't expand beyond the array edge.
    max_size = gdims[i] - start[i];
    count[i] = expand_region(maplen, map, region_size, stride[i], max_size);
    region_size = region_size * count[i];
  }
  return region_size;
}




int create_mpi_datatypes(const MPI_Datatype basetype,const int msgcnt,const PIO_Offset dlen, const PIO_Offset mindex[],const int mcount[],
			 int *mfrom, MPI_Datatype mtype[])
{
  PIO_Offset bsizeT[msgcnt];
  int pos;
  int ii;
  PIO_Offset i8blocksize;
  int blocksize;
  PIO_Offset *lindex = NULL;
#ifdef _MPISERIAL
  mtype[0] = basetype * blocksize;
#else
  if(mindex != NULL){
    lindex = (PIO_Offset *) malloc(dlen * sizeof(PIO_Offset));
    memcpy(lindex, mindex, (size_t) (dlen*sizeof(PIO_Offset)));
  }
  bsizeT[0]=0;
  mtype[0] = MPI_DATATYPE_NULL;
  pos = 0;
  ii = 0;
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
	    for(int j=0;j<len;j++)
	      displace[j] = (int) (lindex[pos+j]);
	  }else{
	    int k=0;
	    for(int j=0;j<dlen;j++)
	      if(mfrom[j]==i)
		displace[k++] = (int) (lindex[j]);
	  }
	    
	}else{
	  for(int j=0;j<mcount[i];j++)
	    (lindex+pos)[j]++;
	  for(int j=0;j<len;j++){
	    displace[j]= ((lindex+pos)[j*blocksize]-1);
	  }
	}
	CheckMPIReturn(MPI_Type_create_indexed_block(len, blocksize, displace, basetype, mtype+i),__FILE__,__LINE__);
	CheckMPIReturn(MPI_Type_commit(mtype+i), __FILE__,__LINE__);
	pos+=mcount[i];

      }
    }

  }
  
  free(lindex);
#endif
  return PIO_NOERR;

}


int define_iodesc_datatypes(const iosystem_desc_t ios, io_desc_t *iodesc)
{
  int i;
  if(ios.ioproc){
    if(iodesc->rtype==NULL){
      int ntypes = iodesc->nrecvs;
      iodesc->rtype = (MPI_Datatype *) malloc(ntypes * sizeof(MPI_Datatype));
      for(i=0; i<ntypes; i++){
        iodesc->rtype[i] = MPI_DATATYPE_NULL;
      }
      iodesc->num_rtypes = ntypes;

      if(iodesc->rearranger==PIO_REARR_SUBSET){
	create_mpi_datatypes(iodesc->basetype, iodesc->nrecvs, iodesc->llen, iodesc->rindex, iodesc->rcount, iodesc->rfrom, iodesc->rtype);
      }else{
	create_mpi_datatypes(iodesc->basetype, iodesc->nrecvs, iodesc->llen, iodesc->rindex, iodesc->rcount, NULL, iodesc->rtype);
      }
#ifndef _MPISERIAL
      /*      {
	MPI_Aint lb;
	MPI_Aint extent;
	for(i=0;i<ntypes;i++){
	  MPI_Type_get_extent(iodesc->rtype[i], &lb, &extent);
	  printf("%s %d %d %d %d \n",__FILE__,__LINE__,i,lb,extent);
	  
	}
	}*/
#endif
    }
  }


  if(iodesc->stype==NULL){
    int ntypes;
    if(iodesc->rearranger==PIO_REARR_SUBSET)
      ntypes = 1;
    else
      ntypes = ios.num_iotasks;


    //  printf("COMP: %d\n",ntypes);


    iodesc->stype = (MPI_Datatype *) malloc(ntypes * sizeof(MPI_Datatype));
    for(i=0; i<ntypes; i++){
      iodesc->stype[i] = MPI_DATATYPE_NULL;
    }
    iodesc->num_stypes = ntypes;

    create_mpi_datatypes(iodesc->basetype, ntypes, iodesc->ndof, iodesc->sindex, iodesc->scount, NULL, iodesc->stype);
#ifndef _MPISERIAL
    /*    {
      MPI_Aint lb;
      MPI_Aint extent;
      for(i=0;i<ntypes;i++){
	MPI_Type_get_extent(iodesc->stype[i], &lb, &extent);
	printf("%s %d %d %d %d \n",__FILE__,__LINE__,i,lb,extent);
      }
      }*/
#endif
  }

  return PIO_NOERR;

}




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
  int maxreq = DEF_P2P_MAXREQ;
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

  iodesc->scount = (int *) calloc(numiotasks,sizeof(int));

  // iodesc->scount is the amount of data sent to each task from the current task
  for(i=0;i<maplen; i++){
    if(dest_ioindex[i] != -1){
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
    recv_buf = (int *) malloc(ntasks * sizeof(int));
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
  //  printf("%s %d\n",__FILE__,__LINE__);

  nrecvs = 0;
  if(ios.ioproc){
    //       printf("recv_buf = ");
    for(i=0;i<ntasks; i++){
      //     printf(" %d ",recv_buf[i]);
      if(recv_buf[i] != 0)
	nrecvs++;
    }
    // printf("\n");

    iodesc->rcount = (int *) calloc(max(1,nrecvs),sizeof(int));
    iodesc->rfrom = (int *) calloc(max(1,nrecvs),sizeof(int));
    

    nrecvs = 0;
    for(i=0;i<ntasks; i++){
      if(recv_buf[i] != 0){
	iodesc->rcount[nrecvs] = recv_buf[i];
	iodesc->rfrom[nrecvs] = i;
	nrecvs++;
      }

    }
    free(recv_buf);
  }

  iodesc->nrecvs = nrecvs;
  if(iodesc->sindex == NULL)
    iodesc->sindex = (PIO_Offset *) calloc(iodesc->ndof,sizeof(PIO_Offset));


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
#ifndef _MPISERIAL
  MPI_Type_size(MPI_OFFSET, &tsize);
#else
  tsize = sizeof(long long);
#endif
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
    for(i=0;i<nrecvs;i++)
      recv_counts[iodesc->rfrom[i]] = iodesc->rcount[i];
    recv_displs[0] = 0;
    for(i=1;i<nrecvs;i++)
      recv_displs[iodesc->rfrom[i]] = recv_displs[iodesc->rfrom[i-1]]+iodesc->rcount[i-1]*tsize;
    if(iodesc->llen>0)
      iodesc->rindex = (PIO_Offset *) calloc(iodesc->llen,sizeof(PIO_Offset));
  }

  //   printf("%d rbuf_size %d\n",ios.comp_rank,rbuf_size);


  // s2rindex is the list of indeces on each compute task
  /*        
  printf("%d s2rindex: ", ios.comp_rank);
  for(i=0;i<iodesc->ndof;i++)
    printf("%ld ",s2rindex[i]);
  printf("\n");
  */
  //  printf("%s %d %ld\n",__FILE__,__LINE__,iodesc->llen);
  //  printf("%s %d %d %d %d %d %d %d\n",__FILE__,__LINE__,send_counts[0],recv_counts[0],send_displs[0],recv_displs[0],sr_types[0],iodesc->llen);
  ierr = pio_swapm( s2rindex, send_counts, send_displs, sr_types, 
		    iodesc->rindex, recv_counts, recv_displs, sr_types,
  		    mycomm, true, false, maxreq);
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
    for(int j=0;j<nrecvs;j++){
      printf("%d rfrom %d ",ios.io_rank,iodesc->rfrom[j]);
      if(j==0)
	for(i=0;i<iodesc->rcount[j];i++)
	  printf("%ld ",iodesc->rindex[i]);
      else  
	for(i=0;i<iodesc->rcount[j];i++)
	  printf("%ld ",iodesc->rindex[rcount[j-1]+i]);
      printf("\n");
      }*/
  
  return ierr;

}

int rearrange_comp2io(const iosystem_desc_t ios, io_desc_t *iodesc, void *sbuf,
			  void *rbuf, const int comm_option, const int fc_options)
{

  bool handshake=true;
  bool isend = false;
  int maxreq = MAX_GATHER_BLOCK_SIZE;
  int ntasks;
  int niotasks;
  int *scount = iodesc->scount;

  int i, tsize;
  int *sendcounts;
  int *recvcounts;
  int *sdispls;
  int *rdispls;
  MPI_Datatype *sendtypes;
  MPI_Datatype *recvtypes;
  MPI_Comm mycomm;
  
  if(iodesc->rearranger == PIO_REARR_BOX){
    mycomm = ios.union_comm;
    niotasks = ios.num_iotasks;
  }else{
    mycomm = iodesc->subset_comm;
    niotasks = 1;
  }  
  MPI_Comm_size(mycomm, &ntasks);

#ifdef _MPISERIAL
  if(iodesc->basetype == 4){
    for(i=0;i<iodesc->llen;i++)
      ((int *) rbuf)[ iodesc->rindex[i] ] = ((int *)sbuf)[ iodesc->sindex[i]];
  }else{
    for(i=0;i<iodesc->llen;i++){
      ((double *) rbuf)[ iodesc->rindex[i] ] = ((double *)sbuf)[ iodesc->sindex[i]];   
    }
  }
#else
  define_iodesc_datatypes(ios, iodesc);
  
  sendcounts = (int *) malloc(ntasks*sizeof(int));
  recvcounts = (int *) malloc(ntasks*sizeof(int));
  sdispls = (int *) malloc(ntasks*sizeof(int));
  rdispls = (int *) malloc(ntasks*sizeof(int));
  sendtypes = (MPI_Datatype *) malloc(ntasks*sizeof(MPI_Datatype));
  recvtypes = (MPI_Datatype *) malloc(ntasks*sizeof(MPI_Datatype));

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
      if(iodesc->rearranger==PIO_REARR_SUBSET){
	recvcounts[ i ] = 1;
	recvtypes[ i ] = iodesc->rtype[i];
      }else{
	recvcounts[ iodesc->rfrom[0] ] = 1;
	recvtypes[ iodesc->rfrom[0] ] = iodesc->rtype[0];
	rdispls[ iodesc->rfrom[0] ] = 0;
	//    printf("%d: rindex[%d] %d\n",ios.comp_rank,0,iodesc->rindex[0]);
	for( i=1;i<iodesc->nrecvs;i++){
	  recvcounts[ iodesc->rfrom[i] ] = 1;
	  recvtypes[ iodesc->rfrom[i] ] = iodesc->rtype[i];

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
    if(scount[i] > 0) {
      sendcounts[io_comprank]=1;
      sendtypes[io_comprank]=iodesc->stype[i];
    }else{
      sendcounts[io_comprank]=0;
    }
  }      

  // Data in sbuf on the compute nodes is sent to rbuf on the ionodes

  pio_swapm( sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, 
	     mycomm, handshake, isend, maxreq);

  free(sendcounts);
  free(recvcounts); 
  free(sdispls);
  free(rdispls);
  free(sendtypes);
  free(recvtypes);
#endif
  return PIO_NOERR;
}

int rearrange_io2comp(const iosystem_desc_t ios, io_desc_t *iodesc, void *sbuf,
			  void *rbuf, const int comm_option, const int fc_options)
{
  

  bool handshake=true;
  bool isend = false;
  int maxreq = MAX_GATHER_BLOCK_SIZE;
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
  if(iodesc->rearranger==PIO_REARR_BOX){
    mycomm = ios.union_comm;
    niotasks = ios.num_iotasks;
  }else{
    mycomm = iodesc->subset_comm;
    niotasks=1;
  }
  MPI_Comm_size(mycomm, &ntasks);

#ifdef _MPISERIAL
  if(iodesc->basetype == 4){
    for(i=0;i<iodesc->llen;i++)
      ((int *) rbuf)[ iodesc->sindex[i] ] = ((int *)sbuf)[ iodesc->rindex[i]];
  }else{
    for(i=0;i<iodesc->llen;i++)
      ((double *) rbuf)[ iodesc->sindex[i] ] = ((double *)sbuf)[ iodesc->rindex[i]];
  }
#else  
  define_iodesc_datatypes(ios, iodesc);

  sendcounts = (int *) calloc(ntasks,sizeof(int));
  recvcounts = (int *) calloc(ntasks,sizeof(int));
  sdispls = (int *) calloc(ntasks,sizeof(int));
  rdispls = (int *) calloc(ntasks,sizeof(int));
  sendtypes = (MPI_Datatype *) malloc(ntasks*sizeof(MPI_Datatype));
  recvtypes = (MPI_Datatype *) malloc(ntasks*sizeof(MPI_Datatype));


  for( i=0;i< ntasks;i++){
    sendtypes[ i ] = MPI_DATATYPE_NULL;
    recvtypes[ i ] = MPI_DATATYPE_NULL;
  }
  if(ios.ioproc){
    for( i=0;i< iodesc->nrecvs;i++){
      if(iodesc->rearranger==PIO_REARR_SUBSET){
	sendcounts[ i ] = 1;
	sendtypes[ i ] = iodesc->rtype[i];
      }else{
	sendcounts[ iodesc->rfrom[i] ] = 1;
	sendtypes[ iodesc->rfrom[i] ] = iodesc->rtype[i];
      }
    }
  }

  for( i=0;i<niotasks; i++){
    int io_comprank = ios.ioranks[i];
    if(iodesc->rearranger==PIO_REARR_SUBSET)
      io_comprank=0;
    if(scount[i] > 0) {
      recvcounts[io_comprank]=1;
      recvtypes[io_comprank]=iodesc->stype[i];
    }
  } 
  
  //
  // Data in sbuf on the ionodes is sent to rbuf on the compute nodes
  //

  pio_swapm( sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, 
	     mycomm, handshake,isend, maxreq);

  free(sendcounts);
  free(recvcounts); 
  free(sdispls);
  free(rdispls);
  free(sendtypes);
  free(recvtypes);
#endif
 
  return PIO_NOERR;

}

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

  iodesc->rearranger = PIO_REARR_BOX;

  iodesc->ndof = maplen;
  gstride[ndims-1]=1;
  for(int i=ndims-2;i>=0; i--)
    gstride[i]=gstride[i+1]*gsize[i+1];

#ifndef _MPISERIAL
  MPI_Type_size(MPI_OFFSET, &tsize);
#endif

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
    iodesc->llen=1;
    for(i=0;i<ndims;i++)
      iodesc->llen *= iodesc->firstregion->count[i];
  }

  for( i=0;i<nioprocs; i++){
    int io_comprank = ios.ioranks[i];
    recvlths[ io_comprank ] = 1;
    rdispls[ io_comprank ] = i*tsize;
  }      

  //  The length of each iomap
  //  iomaplen = calloc(nioprocs, sizeof(PIO_Offset));
  pio_swapm(&(iodesc->llen), sndlths, sdispls, dtypes,
	    iomaplen, recvlths, rdispls, dtypes, 	
	    ios.union_comm, false, false, MAX_GATHER_BLOCK_SIZE);

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
		ios.union_comm, false, false, MAX_GATHER_BLOCK_SIZE);
      
      // The start from iotask i is sent to all compute tasks
      pio_swapm(iodesc->firstregion->start,  sndlths, sdispls, dtypes,
		start, recvlths, rdispls, dtypes, 
		ios.union_comm, false, false, MAX_GATHER_BLOCK_SIZE);

      for(k=0;k<maplen;k++){
	PIO_Offset gcoord[ndims], lcoord[ndims];
	bool found=true;
	gindex_to_coord(ndims, compmap[k], gstride, gcoord);
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

  {

  for(k=0; k<maplen; k++){
    if(dest_ioproc[k] == -1 && compmap[k]>=0){
      fprintf(stderr,"No destination found for compmap[%d] = %ld\n",k,compmap[k]);
      MPI_Abort(MPI_COMM_WORLD,0);
    }
  }
  }
  compute_counts(ios, iodesc, maplen, dest_ioproc, dest_ioindex, ios.union_comm);
  if(ios.ioproc){
    compute_maxIObuffersize(ios.io_comm, iodesc);
  }
  return PIO_NOERR;
}

int compare_offsets(const void *a,const void *b) 
{
  mapsort *x = (mapsort *) a;
  mapsort *y = (mapsort *) b;
  return (int) (x->iomap - y->iomap);
}    

//
// Find all regions. 
//
void get_start_and_count_regions(const MPI_Comm io_comm, io_desc_t *iodesc, const int gdims[],const PIO_Offset map[])
{
  int i;
  int nmaplen;
  int regionlen;
  io_region *region;

  nmaplen = 0;
  region = iodesc->firstregion;
  while(map[nmaplen++]<0);
  nmaplen--;
  region->loffset=nmaplen;

  iodesc->maxregions = 1;
  while(nmaplen < iodesc->llen){
    // Here we find the largest region from the current offset into the iomap
    // regionlen is the size of that region and we step to that point in the map array
    // until we reach the end 

    regionlen = find_first_region(iodesc->ndims, gdims, iodesc->llen-nmaplen, 
					map+nmaplen, region->start, region->count);

    pioassert(region->start[0]>=0,"failed to find region",__FILE__,__LINE__);

    nmaplen = nmaplen+regionlen;
    if(region->next==NULL && nmaplen<iodesc->llen){
      region->next = alloc_region(iodesc->ndims);
      // The offset into the local array buffer is the sum of the sizes of all of the previous regions (loffset) 
      region=region->next;
      region->loffset = nmaplen;
      // The calls to the io library are collective and so we must have the same number of regions on each
      // io task iodesc->maxregions will be the total number of regions on this task 
      iodesc->maxregions++;
    }
  }
  // pad maxregions on all tasks to the maximum and use this to assure that collective io calls are made.
#ifndef _MPISERIAL  
  MPI_Allreduce(MPI_IN_PLACE,&(iodesc->maxregions), 1, MPI_INTEGER, MPI_MAX, io_comm);
#endif
}

void default_subset_partition(const iosystem_desc_t ios, io_desc_t *iodesc)
{
  int taskratio = ios.num_comptasks/ios.num_iotasks;
  int color;
  int key;

  /* Create a new comm for each subset group with the io task in rank 0 and
     only 1 io task per group */

  if(ios.ioproc)
    key=0;
  else
    key=ios.comp_rank%taskratio+1;
  color = ios.comp_rank/taskratio;

  MPI_Comm_split(ios.comp_comm, color, key, &(iodesc->subset_comm));

}



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


  int maxreq = MAX_GATHER_BLOCK_SIZE;

  int rank, ntasks;
  size_t pio_offset_size=sizeof(PIO_Offset);

  /* subset partitions each have exactly 1 io task which is task 0 of that subset_comm */ 
  /* TODO: introduce a mechanism for users to define partitions */
  default_subset_partition(ios, iodesc);

  MPI_Comm_rank(iodesc->subset_comm, &rank);
  MPI_Comm_size(iodesc->subset_comm, &ntasks);

  if(ios.ioproc)
    pioassert(rank==0,"Bad io rank in subset create",__FILE__,__LINE__);
  else
    pioassert(rank>=0 && rank<ntasks,"Bad comp rank in subset create",__FILE__,__LINE__);
  totalgridsize=1;
  for(i=0;i<ndims;i++)
    totalgridsize*=gsize[i];

  
  iodesc->ndof = maplen;
  iodesc->scount = (int *) calloc(1,sizeof(int));

  for(i=0;i<maplen;i++){
    pioassert(compmap[i]>=-1 && compmap[i]<totalgridsize, "Compmap value out of bounds",__FILE__,__LINE__);
    if(compmap[i]>=0){
      (iodesc->scount[0])++;
    }
  }

  iodesc->rearranger = PIO_REARR_SUBSET;
  iodesc->sindex = (PIO_Offset *) calloc(iodesc->scount[0],pio_offset_size); 
  j=0;
  for(i=0;i<maplen;i++){
    if(compmap[i]>=0){
      iodesc->sindex[j++]=i;
    }
  }

  if(ios.ioproc){
    iodesc->rcount = (int *) malloc(ntasks *sizeof(int));
  }
  
  // Pass the reduced maplen (without holes) from each compute task to its associated IO task

  pio_fc_gather( (void *) iodesc->scount, 1, MPI_INT,
		 (void *) iodesc->rcount, 1, MPI_INT, 
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
    if(iodesc->llen>0){
      srcindex = (PIO_Offset *) calloc(iodesc->llen,pio_offset_size);
    }
  }else{
    for(i=0;i<ntasks;i++){
      recvlths[i]=0;
      rdispls[i]=0;
    }
  }
  // Pass the sindex from each compute task to its associated IO task

  pio_fc_gatherv((void *) iodesc->sindex, iodesc->scount[0], PIO_OFFSET,
		 (void *) srcindex, recvlths, rdispls, PIO_OFFSET,  
		 0, iodesc->subset_comm, maxreq);

  if(ios.ioproc){
    map = (mapsort *) malloc(iodesc->llen * sizeof(mapsort));    
    iomap = (PIO_Offset *) calloc(iodesc->llen,pio_offset_size);
  }

  // Now pass the compmap, skipping the holes

  PIO_Offset *shrtmap;
  if(maplen>iodesc->scount[0]){
    shrtmap = (PIO_Offset *) calloc(iodesc->scount[0],pio_offset_size);
    j=0;
    for(i=0;i<maplen;i++)
      if(compmap[i]>=0)
	shrtmap[j++]=compmap[i];
  }else{
    shrtmap = compmap;
  }

  pio_fc_gatherv((void *) shrtmap, iodesc->scount[0], PIO_OFFSET,
		 (void *) iomap, recvlths, rdispls, PIO_OFFSET,  
		 0, iodesc->subset_comm, maxreq);

  if(shrtmap != compmap)
    free(shrtmap);

  if(ios.ioproc){
    int pos=0;
    int k=0;
    for(i=0;i<ntasks;i++){
      for(j=0;j<iodesc->rcount[i];j++){
	(map+k)->rfrom = i;
	(map+k)->soffset = srcindex[pos+j];
	(map+k++)->iomap = iomap[pos+j];
      }
      pos += iodesc->rcount[i];
    }
    // sort the mapping, this will transpose the data into IO order        
    qsort(map, iodesc->llen, sizeof(mapsort), compare_offsets); 

    iodesc->rindex = (PIO_Offset *) calloc(iodesc->llen,pio_offset_size);
    iodesc->rfrom = (int *) calloc(iodesc->llen,sizeof(int));
  }
  int cnt[ntasks];
  int sndlths[ntasks];
  int sdispls[ntasks];
  MPI_Datatype dtypes[ntasks];
  for(i=0;i<ntasks;i++){
    cnt[i]=rdispls[i];

    /* offsets to swapm are in bytes */
    rdispls[i]*=pio_offset_size;
    sdispls[i]=0;
    sndlths[i]=0;
    dtypes[i]=PIO_OFFSET;
  }
  sndlths[0]=iodesc->scount[0];
  for(i=0;i<iodesc->llen;i++){
    iodesc->rfrom[i] = (map+i)->rfrom;
    iodesc->rindex[i]=i;
    iomap[i] = (map+i)->iomap;
    srcindex[ (cnt[iodesc->rfrom[i]])++   ]=(map+i)->soffset;
  }

  pio_swapm((void *) srcindex, recvlths, rdispls, dtypes, 
	    (void *) iodesc->sindex, sndlths, sdispls, dtypes, 
	    iodesc->subset_comm, hs, isend, maxreq);

  if(ios.ioproc){
    get_start_and_count_regions(ios.io_comm,iodesc,gsize,iomap);
  
    if(iomap != NULL)
      free(iomap);
    
    if(map != NULL)
      free(map);

    if(srcindex != NULL)
      free(srcindex);
    
    compute_maxIObuffersize(ios.io_comm, iodesc);

    iodesc->nrecvs=ntasks;
  }

  return ierr;


}
  
    
  
