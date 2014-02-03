#include <pio.h>
#include <pio_internal.h>
#define DEF_P2P_MAXREQ 64

int gindex_to_coord(const PIO_Offset gindex, const PIO_Offset gstride[], const int ndim, PIO_Offset *gcoord)
{
  PIO_Offset tempindex;

  tempindex = gindex;
  for(int i=ndim-1;i>0;i--){
    gcoord[i] = tempindex/gstride[i-1];
    tempindex -= gcoord[i]*gstride[i-1];
  }
  gcoord[0] = tempindex;
  return PIO_NOERR;
}


bool find_ioproc(const PIO_Offset gcoord[], const PIO_Offset lb[], const PIO_Offset ub[], const PIO_Offset lstride[],
		 const int ndim, const int nioproc, int *ioproc, PIO_Offset *ioindex)
{
  bool found[ndim];
  int ii, i, j, k;
  int decompstep[ndim];
  PIO_Offset lcoord[ndim];
  PIO_Offset lindex;

  *ioindex = -1;
  *ioproc = -1;


  for(j=0;j<ndim;j++){
    decompstep[j]=1;
    found[j]=false;
  }

  for(i=0,ii=0; i<nioproc; i++){
    for(j=0;j<ndim;j++,ii++){
      if(ub[ii] != ub[j]){
	decompstep[j]=i;
      }
    }
  }
  for(i=0,ii=0; i<nioproc; i++){
      for(j=0;j<ndim;j++,ii++){
	if(!found[j]){
	  if(lb[ii] <= gcoord[j] && gcoord[j] < ub[ii]){
	    found[j] = true;
	    *ioproc = i;
	    break;
	  }
	}
      }
  }
  *ioindex = -1;
  if(found[0]){
    *ioindex = gcoord[0] - lb[(*ioproc)*ndim];
    for(j=1;j<ndim;j++)
      if(found[j])
	*ioindex += gcoord[j] - lb[j+(*ioproc)*ndim]*lstride[j-1 + (*ioproc)*ndim] ;
  
  }


  return true;
}


int compute_dest(const int nmap, const PIO_Offset compmap[],const PIO_Offset *start, const PIO_Offset *count, const int gsize[],
		 const int ndim, const int nioproc, int *dest_ioproc, PIO_Offset *dest_ioindex)
{
  PIO_Offset lb[nioproc*ndim];
  PIO_Offset ub[nioproc*ndim];
  PIO_Offset lstride[nioproc*ndim];
  PIO_Offset gstride[ndim];
  PIO_Offset gcoord[ndim];
  PIO_Offset gindex;
  int ierr, i, j, ii;

  for( i=0, ii=0;i<nioproc;i++){
    for( j=0;j<ndim;j++, ii++){
      lb[ii] = start[ii];
      ub[ii] = lb[ii] + count[ii];
      lstride[ii]=1;
    }
  }
  
  gstride[0] = gsize[0];
  for( i=1;i<ndim;i++){
    gcoord[i] = -1;
    gstride[i] = gsize[i]*gstride[i-1];
  }
  
  for(i=0, ii=0;i<nioproc;i++){
    lstride[ii] = count[ii];
    for( j=1;j<ndim;j++,ii++){
      lstride[ii] = count[ii]*lstride[ii-1];
    }
  }
  for( i=0;i< nmap; i++){
    if(compmap[i]==0){  // sender hole
      dest_ioproc[i]=-1;
      dest_ioindex[i]=-1;
    }else{
      gindex = compmap[i]-1;
      ierr = gindex_to_coord(gindex, gstride, ndim, gcoord);


      if(! find_ioproc(gcoord, lb,  ub, lstride, 
		       ndim, nioproc, dest_ioproc+i, dest_ioindex+i)){
	fprintf(stderr,"No destination for compmap point = %ld\n",compmap[i]);
	piodie("Could not find ioproc",__FILE__,__LINE__);
      }
    }
  }
 
  return PIO_NOERR;
}

int create_mpi_datatypes(const MPI_Datatype basetype,const int msgcnt,PIO_Offset mindex[],const int mcount[],MPI_Datatype mtype[])
{
  PIO_Offset bsizeT[msgcnt];
  int pos;
  int ii;
  PIO_Offset *indptr;
  PIO_Offset i8blocksize;
  MPI_Datatype newtype;
  int blocksize;
  bool freenewtype=true;


  mtype[0] = MPI_DATATYPE_NULL;
  pos = 0;
  ii = 0;
  if(msgcnt>0){
    for(int i=0;i<msgcnt;i++){
      if(mcount[i]>0){
	indptr = mindex+pos;
	bsizeT[ii] = GCDblocksize(mcount[i], indptr);
	ii++;
	pos+=mcount[i];
      }
    }
    blocksize = (int) lgcd_array(ii ,bsizeT);

    CheckMPIReturn(MPI_Type_contiguous(blocksize, basetype, &newtype),__FILE__,__LINE__);
    CheckMPIReturn(MPI_Type_commit(&newtype), __FILE__,__LINE__);

    pos = 0;
    for(int i=0;i< msgcnt; i++){
      if(mcount[i]>0){
	int len = mcount[i]/blocksize;
	if(len>1){
	  int *displace = (int *) malloc(len * sizeof(int));
	  if(blocksize==1)
	    for(int j=0;j<len;j++)
	      displace[j] = (int) (mindex+pos)[j];
	  else
	    for(int j=0;j<len;j++)
	      (mindex+pos)[j]++;

	  CheckMPIReturn(MPI_Type_create_indexed_block(len, 1, displace, newtype, mtype+i),__FILE__,__LINE__);
	  CheckMPIReturn(MPI_Type_commit(mtype+i), __FILE__,__LINE__);
	  free(displace);
	}else{
	  mtype[i] = newtype;
	  freenewtype= false;
	}
	pos+=mcount[i];

      }
    }
    if(freenewtype) CheckMPIReturn(MPI_Type_free(&newtype),__FILE__,__LINE__);
  }
  return PIO_NOERR;

}

int compute_counts(iosystem_desc_t *ios, io_desc_t *iodesc,const int niomap)
{
  PIO_Offset *sindex;
  PIO_Offset *rindex;
  PIO_Offset *s2rindex;
  int i;
  int niotasks;
  int ncomptasks;
  int iorank;
  int *scount;
  MPI_Datatype *sr_types;
  int *send_counts;
  int *send_displs;
  int *recv_counts;
  int *recv_displs;
  int *recv_buf;
  int *rcount;
  int *rfrom;
  int nrecvs;
  int rbuf_size;
  int pio_maxreq;  
  int ierr;
  int io_comprank;
  MPI_Datatype dtype;
  int ioindex;
  int tsize, tsizei, tsizel;
  int ndof= iodesc->ndof;
  bool pio_hs;
  bool pio_isend;
  

  pio_hs = false;
  pio_isend = false;
  pio_maxreq = DEF_P2P_MAXREQ;
  niotasks = ios->num_iotasks;
  ncomptasks = ios->num_comptasks;
  scount = (int *) calloc(niotasks,sizeof(int));
  // scount is the amount of data sent to each task from the current task
  for(i=0;i<ndof; i++){
    iorank = iodesc->dest_ioproc[i];
    if(iorank != -1){
      (scount[iorank])++;
    }
  }

  iodesc->scount = scount;

  send_counts = (int *) calloc(ncomptasks,sizeof(int));
  send_displs = (int *) calloc(ncomptasks, sizeof(int));
  recv_counts = (int *) calloc(ncomptasks , sizeof(int));
  recv_displs = (int *) calloc(ncomptasks , sizeof(int));
  sr_types = (MPI_Datatype *) malloc(sizeof(MPI_Datatype)*ncomptasks);

  for(i=0;i<ncomptasks;i++){
    sr_types[i] = MPI_INT;
  }

  for(i=0;i<niotasks;i++){
    int io_comprank = ios->ioranks[i];
    send_counts[io_comprank] = 1;
    send_displs[io_comprank] = i*sizeof(int);
  }


  if(ios->ioproc){
    recv_buf = (int *) malloc(ncomptasks * sizeof(int));
    for(i=0;i<ncomptasks;i++){
      recv_buf[i] = 0;
      recv_counts[i] = 1;
      recv_displs[i] = i*sizeof(int);
    }
    rbuf_size = ncomptasks;
  }else{
    recv_buf = (int *) malloc(sizeof(int));
     rbuf_size = 1;
    recv_buf[0]=0;
  }
  // Share the scount from each compute task to all io tasks

  
  ierr = pio_swapm( ncomptasks, ios->comp_rank, scount, 
		    send_counts, send_displs, sr_types, recv_buf,  recv_counts, recv_displs,
		    sr_types,ios->union_comm, pio_hs, pio_isend, pio_maxreq);


  nrecvs = 0;
  if(ios->ioproc){
    for(i=0;i<ncomptasks; i++){
      if(recv_buf[i] != 0)
	nrecvs++;
    }
    rcount = (int *) calloc(max(1,nrecvs),sizeof(int));
    iodesc->rfrom = (int *) calloc(max(1,nrecvs),sizeof(int));
    
    rfrom = iodesc->rfrom;
    nrecvs = 0;
    for(i=0;i<ncomptasks; i++){
      if(recv_buf[i] != 0){
	rcount[nrecvs] = recv_buf[i];
	rfrom[nrecvs] = i;
	nrecvs++;
      }

    }
  }else{
    rcount = (int *) malloc(sizeof(int));
    rcount[0]=0;
  }
    
  iodesc->nrecvs = nrecvs;
  free(recv_buf);

  int *tempcount = (int *) calloc(niotasks, sizeof(int));
  int *spos = (int *) calloc(niotasks, sizeof(int));
  sindex = (PIO_Offset *) calloc(ndof,sizeof(PIO_Offset));
  s2rindex = (PIO_Offset *) calloc(ndof,sizeof(PIO_Offset));

  spos[0]=0;
  for(i=1;i<niotasks;i++)
    spos[i] = spos[i-1] + scount[i-1];
    

  for(i=0;i<ndof;i++){
    iorank =iodesc->dest_ioproc[i]; 
    ioindex = iodesc->dest_ioindex[i];
    if(iorank > -1){
      sindex[spos[iorank]+tempcount[iorank]] = i;
      s2rindex[spos[iorank]+tempcount[iorank]] = ioindex;
      (tempcount[iorank])++;
    }
  }
  free(tempcount);

  for(i=0;i<ncomptasks;i++){
    send_counts[i] = 0;
    send_displs[i]  = 0;
    recv_counts[i] = 0;
    recv_displs[i]   =0;
  }
  for(i=0;i<niotasks;i++){
    io_comprank = ios->ioranks[i];
    send_counts[io_comprank] = scount[i];
    send_displs[io_comprank]  = spos[i] -1;
  }
  free(spos);

  MPI_Type_size(MPI_LONG_LONG, &tsizel);
  MPI_Type_size(MPI_INT, &tsizei);
  if(sizeof(PIO_Offset) == tsizei){
    dtype = MPI_INT;
    tsize = tsizei;
  }else if(sizeof(PIO_Offset) == tsizel){
    dtype = MPI_LONG;
    tsize = tsizel;
  }

  for(i=0; i<ncomptasks; i++){
    sr_types[i] = dtype;
  }



  if(ios->ioproc){
    rbuf_size=max(1,niomap);
    
    for(i=0;i<nrecvs;i++)
      recv_counts[rfrom[i]] = rcount[i];
    recv_displs[0] = 0;
    for(i=1;i<nrecvs;i++)
      recv_displs[rfrom[i]] = recv_displs[rfrom[i-1]]+rcount[i-1]*tsize;
  }else{
    rbuf_size=1;
  }
  rindex = (PIO_Offset *) calloc(max(1,niomap),sizeof(PIO_Offset));


  ierr = pio_swapm( ncomptasks, ios->comp_rank, s2rindex, 
		    send_counts, send_displs, sr_types, rindex, recv_counts, recv_displs,
		    sr_types,ios->union_comm, pio_hs, pio_isend, pio_maxreq);


  free(s2rindex);
  free(sr_types);
  free(send_counts);
  free(send_displs);
  free(recv_counts);
  free(recv_displs);

  iodesc->rtype = NULL;
  iodesc->stype = NULL;

  if(ios->ioproc){
    iodesc->rtype = (MPI_Datatype *) malloc(max(1,nrecvs)*sizeof(MPI_Datatype));
    create_mpi_datatypes(iodesc->basetype, nrecvs, rindex, rcount, iodesc->rtype);
  }
    // for(i=0;i<nrecvs;i++)
    // iodesc->rtype[i] = iodesc->basetype;


  iodesc->stype = (MPI_Datatype *) malloc(niotasks*sizeof(MPI_Datatype));
  // for(i=0;i<niotasks;i++)
  //   iodesc->stype[i]=iodesc->basetype;
  create_mpi_datatypes(iodesc->basetype, niotasks, sindex, scount, iodesc->stype);



  free(rindex);
  free(rcount);
  free(sindex);
  return ierr;

}

int box_rearrange_comp2io(iosystem_desc_t *ios, io_desc_t *iodesc, void *sbuf,
			  void *rbuf, const int comm_option, const int fc_options)
{

  bool handshake=true;
  bool isend = false;
  int maxreq = MAX_GATHER_BLOCK_SIZE;
  int nprocs = ios->num_comptasks;
  int *scount = iodesc->scount;
  int *sendcounts = (int *) calloc(nprocs , sizeof(int));
  int *recvcounts = (int *) calloc(nprocs , sizeof(int));
  int *sdispls =  (int *) calloc(nprocs , sizeof(int));
  int *rdispls =  (int *) calloc(nprocs , sizeof(int));
  MPI_Datatype *sendtypes = (MPI_Datatype *) malloc(nprocs* sizeof(MPI_Datatype));
  MPI_Datatype *recvtypes = (MPI_Datatype *) malloc(nprocs* sizeof(MPI_Datatype));
  int i, tsize;

  for(i=0;i<nprocs;i++){
      recvtypes[ i ] = MPI_DATATYPE_NULL;
      sendtypes[ i] =  MPI_DATATYPE_NULL;
  }


  if(ios->ioproc){
    for( i=0;i<iodesc->nrecvs;i++){
      recvcounts[ iodesc->rfrom[i] ] = iodesc->llen;
      //      recvtypes[ iodesc->rfrom[i] ] = iodesc->rtype[i];
      recvtypes[ iodesc->rfrom[i] ] = iodesc->basetype;
      
      MPI_Type_size(iodesc->rtype[i],&tsize);
      rdispls[ iodesc->rfrom[i] ] = i*tsize;
      //recvtypes[ iodesc->rfrom[i] ] = iodesc->basetype;
      
    }
  }else{
    for( i=0;i<iodesc->nrecvs;i++){
      recvcounts[ iodesc->rfrom[i] ] = 0;
    }
  }  

  for( i=0;i<ios->num_iotasks; i++){
    int io_comprank = ios->ioranks[i];
    if(scount[i] > 0) {
      sendcounts[io_comprank]=1;
      sendtypes[io_comprank]=iodesc->stype[i];
      MPI_Type_size(iodesc->stype[i],&tsize);
      sdispls[io_comprank] = i*tsize;
    }else{
      sendcounts[io_comprank]=0;
    }
  }      

  pio_swapm( nprocs, ios->union_rank, sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, ios->union_comm, handshake, isend, maxreq);

  free(sendtypes);
  free(recvtypes);
  free(sendcounts);
  free(recvcounts);

  return PIO_NOERR;
}

int box_rearrange_io2comp(iosystem_desc_t *ios, io_desc_t *iodesc, void *sbuf,
			  void *rbuf, const int comm_option, const int fc_options)
{
  

  bool handshake=true;
  bool isend = false;
  int maxreq = MAX_GATHER_BLOCK_SIZE;
  int nprocs = ios->num_comptasks;
  int *scount = iodesc->scount;
  int *sendcounts = (int *) calloc(nprocs , sizeof(int));
  int *recvcounts = (int *) calloc(nprocs , sizeof(int));
  int *sdispls =  (int *) calloc(nprocs , sizeof(int));
  int *rdispls =  (int *) calloc(nprocs , sizeof(int));
  MPI_Datatype *sendtypes = (MPI_Datatype *) malloc(nprocs * sizeof(MPI_Datatype));
  MPI_Datatype *recvtypes = (MPI_Datatype *) malloc(nprocs * sizeof(MPI_Datatype));
  int i, tsize;
  

  for( i=0;i< nprocs;i++){
    sendtypes[ i ] = MPI_DATATYPE_NULL;
    recvtypes[ i ] = MPI_DATATYPE_NULL;
  }
  if(ios->ioproc){
    for( i=0;i< iodesc->nrecvs;i++){
      sendcounts[ iodesc->rfrom[i] ] = 1;
      sendtypes[ iodesc->rfrom[i] ] = iodesc->rtype[i];
      MPI_Type_size(iodesc->rtype[i], &tsize);
      sdispls[ iodesc->rfrom[i] ] = i*tsize;
    }
  }
    

  for( i=0;i<ios->num_iotasks; i++){
    int io_comprank = ios->ioranks[i];
    if(scount[i] > 0) {
      recvcounts[io_comprank]=1;
      recvtypes[io_comprank]=iodesc->stype[i];

      MPI_Type_size(iodesc->stype[i],&tsize);
      rdispls[ io_comprank ] = i*tsize;

    }else{
      recvcounts[io_comprank]=0;
      recvtypes[io_comprank]=MPI_DATATYPE_NULL;
    }
  }      

  pio_swapm( nprocs, ios->union_rank, sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, ios->union_comm, handshake, isend, maxreq);

  free(sendtypes);
  free(recvtypes);
  free(sendcounts);
  free(recvcounts);
  free(sdispls);
  free(rdispls);
  return PIO_NOERR;

}



int box_rearrange_create(iosystem_desc_t *ios,const int maplen, const PIO_Offset compmap[], const int gsize[],
			 const int ndim, const int nioproc, io_desc_t *iodesc)
{
  PIO_Offset *start;
  PIO_Offset *count;
  MPI_Datatype dtype;
  int ierr=PIO_NOERR;
  int niomap;
  int nprocs = ios->num_comptasks;
  int nioprocs = ios->num_iotasks;
  int i, tsizel, tsizei, tsize;
  iodesc->ndof = maplen;

  iodesc->dest_ioproc = (int *) malloc(maplen * sizeof(int));
  iodesc->dest_ioindex = (PIO_Offset *) malloc(maplen * sizeof(PIO_Offset));

  start = (PIO_Offset *) calloc(ndim*nioprocs, sizeof(PIO_Offset));
  count = (PIO_Offset *) calloc(ndim*nioprocs, sizeof(PIO_Offset));

  MPI_Type_size(MPI_LONG_LONG, &tsizel);
  MPI_Type_size(MPI_INT, &tsizei);
  if(sizeof(PIO_Offset) == tsizei){
    dtype = MPI_INT;
    tsize = tsizei;
  }else if(sizeof(PIO_Offset) == tsizel){
    dtype = MPI_LONG;
    tsize = tsizel;
  }

  int * sndlths = (int *) malloc( nprocs*sizeof(int));
  int * sdispls = (int *) malloc( nprocs*sizeof(int));
  int * recvlths = (int *) malloc( nprocs*sizeof(int));
  int * rdispls = (int *) malloc( nprocs*sizeof(int));
  MPI_Datatype * dtypes = (MPI_Datatype *) malloc(nprocs*sizeof(MPI_Datatype));

  for(i=0;i<nprocs;i++){
    sdispls[i] = 0;
    rdispls[i] = 0;
    recvlths[i] = 0;
    dtypes[ i ] = dtype;
    sndlths[ i ] = 0;
  }

  if(ios->ioproc){
    for( i=0;i<nprocs;i++){
      sndlths[ i ] = ndim;
    }
  }
  for( i=0;i<ios->num_iotasks; i++){
    int io_comprank = ios->ioranks[i];
    recvlths[ io_comprank ] = ndim;
    rdispls[ io_comprank ] = ndim*i*tsize;
  }      

  pio_swapm(ios->num_comptasks, ios->union_rank, iodesc->count, sndlths, sdispls, dtypes,
	    count, recvlths, rdispls, dtypes, ios->union_comm, false, false, MAX_GATHER_BLOCK_SIZE);


  pio_swapm(ios->num_comptasks, ios->union_rank, iodesc->start,  sndlths, sdispls, dtypes,
	    start, recvlths, rdispls, dtypes, ios->union_comm, false, false, MAX_GATHER_BLOCK_SIZE);


  ierr = compute_dest(maplen, compmap,  start, count, gsize, ndim, 
		      ios->num_iotasks, iodesc->dest_ioproc, iodesc->dest_ioindex);


  niomap = 1;
  for( i=0;i<ndim;i++)
    niomap *= iodesc->count[i];

  compute_counts(ios, iodesc, niomap);


  free(iodesc->dest_ioindex);
  free(iodesc->dest_ioproc);


  return PIO_NOERR;
}

