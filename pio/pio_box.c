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




/* Given the compmap and the start and count for each IOtask compute the destination ioproc and index */


int compute_dest(const int nmap, const PIO_Offset compmap[],const PIO_Offset *start, const PIO_Offset *count, const int gsize[],
		 const int ndims, const int nioproc, int *dest_ioproc, PIO_Offset *dest_ioindex)
{
  PIO_Offset iorank, ioindex;
  int ioproc;
  int gstride, lstride;
  int i, ii;

  for(i=0;i<nmap;i++){
    dest_ioproc[i] = -1;
    dest_ioindex[i] = -1;
  }
  for(ioproc=0, ii=0;ioproc<nioproc;ioproc++){
    gstride=1;
    lstride = 1;
    ioindex = 0;
    iorank=0;
    for(i=0;i<ndims;i++, ii++){
      for(iorank=start[ii]*lstride;iorank<(start[ii]+count[ii])*lstride;iorank+=lstride){
	for(int k=0;k<nmap;k++){
	  if(compmap[k] == iorank+1){
	    dest_ioproc[k] = ioproc;
	    dest_ioindex[k] = iorank;
	  }
	}
      }
      lstride *= count[ii-1];

    }

  }
  
  for(i=0;i<nmap;i++){
  //   printf("ioproc %d ioindex %d\n",dest_ioproc[i],dest_ioindex[i]);
    if(dest_ioproc[i]==-1 && compmap[i] >= 0){
      fprintf(stderr,"No destination found for %ld\n",compmap[i]);
    }
  }

 
  return PIO_NOERR;
}

int create_mpi_datatypes(const MPI_Datatype basetype,const int msgcnt,PIO_Offset mindex[],const int mcount[],MPI_Datatype mtype[])
{
  PIO_Offset bsizeT[msgcnt];
  int pos;
  int ii;
  PIO_Offset i8blocksize;
  MPI_Datatype newtype;
  int blocksize;
  bool freenewtype=true;

  bsizeT[0]=0;
  mtype[0] = MPI_DATATYPE_NULL;
  pos = 0;
  ii = 0;
  if(msgcnt>0){
    for(int i=0;i<msgcnt;i++){
      printf("mcount %d mindex %ld\n",mcount[i],mindex[i]);
      if(mcount[i]>0){
	bsizeT[ii] = GCDblocksize(mcount[i], mindex+pos);
		for(int j=0;j<mcount[i];j++)
		  printf(" %d ",mindex[pos+j]);
		printf("\n bsizet[%d] %ld\n",ii,bsizeT[ii]);
	ii++;
	pos+=mcount[i];
      }
    }
    blocksize = (int) lgcd_array(ii ,bsizeT);

    CheckMPIReturn(MPI_Type_contiguous(blocksize, basetype, &newtype),__FILE__,__LINE__);
    CheckMPIReturn(MPI_Type_commit(&newtype), __FILE__,__LINE__);
    
    printf("%d blocksize %d basetype %d newtype %d %d %d \n",msgcnt,blocksize,basetype,newtype,mindex[0],mcount[0]);


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

  //  for(i=0;i<niotasks;i++)
  //   printf("scount = %d\n",scount[i]);

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
  ierr = pio_swapm( scount, send_counts, send_displs, sr_types, 
                    recv_buf,  recv_counts, recv_displs, sr_types,
		    ios->union_comm, pio_hs, pio_isend, pio_maxreq);


  nrecvs = 0;
  if(ios->ioproc){
    //    printf("recv_buf = ");
    for(i=0;i<ncomptasks; i++){
      //   printf(" %d ",recv_buf[i]);
      if(recv_buf[i] != 0)
	nrecvs++;
    }
    // printf("\n");

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

    printf(" ior %d ioi %d\n",iorank,ioindex);

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
    send_displs[io_comprank]  = spos[i] ;
  }
  free(spos);

      printf("%s %d\n",__FILE__,__LINE__);
  MPI_Type_size(MPI_LONG_LONG, &tsizel);
      printf("%s %d\n",__FILE__,__LINE__);
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
  rindex = (PIO_Offset *) calloc(rbuf_size,sizeof(PIO_Offset));

  printf("%ld %ld %ld\n",s2rindex[0],s2rindex[1],s2rindex[2]);

  // s2rindex is the list of indeces on each compute task

  ierr = pio_swapm( s2rindex, send_counts, send_displs, sr_types, 
		    rindex, recv_counts, recv_displs, sr_types,
		    ios->union_comm, pio_hs, pio_isend, 0);


  //  rindex is an array of the indices of the data to be sent from
  //  this io task to each compute task. 

  for(i=0;i<rbuf_size;i++)
    printf("rindex[%d] %ld\n",i,rindex[i]);
  
  //  MPI_Abort(MPI_COMM_WORLD,0);

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
  
  iodesc->stype = (MPI_Datatype *) malloc(niotasks*sizeof(MPI_Datatype));

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
      // recvcounts[ iodesc->rfrom[i] ] = iodesc->llen;
      recvcounts[ iodesc->rfrom[i] ] = 1;
      recvtypes[ iodesc->rfrom[i] ] = iodesc->rtype[i];
      //     recvtypes[ iodesc->rfrom[i] ] = iodesc->basetype;
            printf("%s %d\n",__FILE__,__LINE__);

      MPI_Type_size(iodesc->rtype[i],&tsize);
      rdispls[ iodesc->rfrom[i] ] = i*tsize;
      
      
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
      printf("%s %d\n",__FILE__,__LINE__);
      MPI_Type_size(iodesc->stype[i],&tsize);
      printf("%s %d\n",__FILE__,__LINE__);
      sdispls[io_comprank] = i*tsize;
    }else{
      sendcounts[io_comprank]=0;
    }
  }      

  
  /*  printf("sbuf %d %d %d\n",((int *) sbuf)[0],((int *) sbuf)[1],((int *) sbuf)[2]);
  for(i=0;i<nprocs;i++){
    printf("%d sendcounts %d sdispls %d %d\n",i,sendcounts[i],sdispls[i],tsize);
  }
  for(i=0;i<nprocs;i++){
    printf("%d recvcounts %d sdispls %d %d\n",i,recvcounts[i],rdispls[i],tsize);
    }*/

  // Data in sbuf on the compute nodes is sent to rbuf on the ionodes


  pio_swapm( sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, 
	     ios->union_comm, handshake, isend, 0);

  if(ios->iomaster)
    printf("iobuf %d %d %d\n",((int *) rbuf)[0],((int *) rbuf)[1],((int *) rbuf)[2]);


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
      printf("%s %d\n",__FILE__,__LINE__);
      MPI_Type_size(iodesc->rtype[i], &tsize);
      sdispls[ iodesc->rfrom[i] ] = i*tsize;
    }
  }
    

  for( i=0;i<ios->num_iotasks; i++){
    int io_comprank = ios->ioranks[i];
    if(scount[i] > 0) {
      recvcounts[io_comprank]=1;
      recvtypes[io_comprank]=iodesc->stype[i];

      printf("%s %d\n",__FILE__,__LINE__);
      MPI_Type_size(iodesc->stype[i],&tsize);
      rdispls[ io_comprank ] = i*tsize;

    }else{
      recvcounts[io_comprank]=0;
      recvtypes[io_comprank]=MPI_DATATYPE_NULL;
    }
  } 
  //
  // Data in sbuf on the ionodes is sent to rbuf on the compute nodes
  //
  pio_swapm( sbuf,  sendcounts, sdispls, sendtypes,
	     rbuf, recvcounts, rdispls, recvtypes, 
	     ios->union_comm, handshake, isend, maxreq);

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
  int ierr=PIO_NOERR;
  int nprocs = ios->num_comptasks;
  int nioprocs = ios->num_iotasks;
  PIO_Offset gstride[ndim];
  PIO_Offset *iomap;
  int offset=0;
  int tsizei, tsizel, tsize, i, j, k;
  MPI_Datatype dtype;

  iodesc->ndof = maplen;
  iodesc->dest_ioproc = (int *) malloc(maplen * sizeof(int));
  iodesc->dest_ioindex = (PIO_Offset *) malloc(maplen * sizeof(PIO_Offset));


  if(ios->ioproc){
    iomap = (PIO_Offset *) calloc(max(1,iodesc->llen), sizeof(PIO_Offset));
    int i=0;

    gstride[0]=1;
    for(int i=1;i<ndim; i++)
    gstride[i]=gstride[i-1]*gsize[i-1];

    for(int idim=0; idim<ndim; idim++){
      iomap[0] += gstride[idim]*iodesc->start[idim];
    }
    for(i=1;i<iodesc->llen;i++){
      iomap[i]=iomap[i-1]+1;
      if((i % iodesc->count[0]) == 0){
	for(int idim=0; idim<ndim; idim++){
	  iomap[i] += gstride[idim]*iodesc->start[idim];
	}
      }
    }
      
  }  
      
  // share the iomap from each io task to all compute tasks

  
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
      sndlths[ i ] = 1;
    }
  }
  for( i=0;i<ios->num_iotasks; i++){
    int io_comprank = ios->ioranks[i];
    recvlths[ io_comprank ] = 1;
    rdispls[ io_comprank ] = i*tsize;
  }      
  PIO_Offset iomaplen[nioprocs];
  //  The length of each iomap
  pio_swapm(iodesc->llen, sndlths, sdispls, dtypes,
	    iomaplen, recvlths, rdispls, dtypes, 
	    ios->union_comm, false, false, MAX_GATHER_BLOCK_SIZE);


  for(i=0; i<nioprocs; i++){
    int io_comprank = ios->ioranks[i];
    for( j=0; j<nprocs ; j++){
      sndlths[ j ] = 0;
      recvlths[ j ] = iomaplen[ i ];
      rdispls[ j ] = 0;
    }
    sndlths[ io_comprank ] = iomaplen[ i ];
    PIO_Offset *liomap = (PIO_Offset *) calloc(iomaplen[i], sizeof(PIO_Offset));
    // The iomap from iotask i is sent to all compute tasks
    pio_swapm(iomap,  sndlths, sdispls, dtypes,
	      liomap, recvlths, rdispls, dtypes, 
	      ios->union_comm, false, false, MAX_GATHER_BLOCK_SIZE);
  

    for(k=0;k<maplen;k++){
      if(iodesc->dest_ioproc[k]=-1){
	for(j=0;j<iomaplen[i];j++){
	  if(compmap[k]==liomap[j]){
	    iodesc->dest_ioproc[k]=i;
	    iodesc->dest_ioindex[k]=j;
	  }
	}
      }
    }
    free(liomap);




  }  

  //  ierr = compute_dest(maplen, compmap,  start, count, gsize, ndim, 
  //		      ios->num_iotasks, iodesc->dest_ioproc, iodesc->dest_ioindex);


  int niomap = 1;
  for( i=0;i<ndim;i++)
    niomap *= iodesc->count[i];

  compute_counts(ios, iodesc, niomap);


  free(iodesc->dest_ioindex);
  free(iodesc->dest_ioproc);


  return PIO_NOERR;
}

