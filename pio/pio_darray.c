#include <pio.h>
#include <pio_internal.h>

PIO_Offset PIO_BUFFER_SIZE_LIMIT= 100000000; // 100MB default limit

#define MALLOC_FILL_ARRAY(type, n, fill, arr) \
  arr = malloc(n * sizeof (type));	      \
  if(fill != NULL)                                       \
    for(int _i=0; _i<n; _i++)			\
      ((type *) arr)[_i] = *((type *) fill)

// Changes to PIO_BUFFER_SIZE_LIMIT only apply to files opened after the change
PIO_Offset PIOc_set_buffer_size_limit(const PIO_Offset limit)
{
  PIO_Offset oldsize; 
  oldsize = PIO_BUFFER_SIZE_LIMIT;
  if(limit>0)
    PIO_BUFFER_SIZE_LIMIT=limit;
  return(oldsize);
}



int pio_write_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid, void *IOBUF, void *fillvalue)
{
  iosystem_desc_t *ios;
  var_desc_t *vdesc;
  int ndims;
  int ierr;
  int msg;
  int mpierr;
  int i;
  int dsize;
  MPI_Status status;
  PIO_Offset usage;
  MPI_Request request;
  int fndims;

  ierr = PIO_NOERR;

  ios = file->iosystem;
  if(ios == NULL){
    fprintf(stderr,"Failed to find iosystem handle \n");
    return PIO_EBADID;
  }
  vdesc = (file->varlist)+vid;

  if(vdesc == NULL){
    fprintf(stderr,"Failed to find variable handle %d\n",vid);
    return PIO_EBADID;
  }
  ndims = iodesc->ndims;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }

  ierr = PIOc_inq_varndims(file->fh, vid, &fndims);
  
  if(ios->ioproc){
    io_region *region;
    int ncid = file->fh;
    int regioncnt;
    void *bufptr;
    void *tmp_buf=NULL;
    int tsize;
    size_t start[ndims+1];
    size_t count[ndims+1];
    int buflen, j, i;
#ifndef _MPISERIAL
    MPI_Type_size(iodesc->basetype, &tsize);
#endif
    region = iodesc->firstregion;

    if(vdesc->record >= 0 && ndims<fndims)
      ndims++;
#ifdef _PNETCDF
    if(file->iotype == PIO_IOTYPE_PNETCDF){
      // make sure we have room in the buffer ;
	ierr = ncmpi_inq_buffer_usage(ncid, &usage);
	usage += tsize*(iodesc->maxiobuflen);
	MPI_Allreduce(MPI_IN_PLACE, &usage, 1,  MPI_LONG_LONG,  MPI_MAX, ios->io_comm);
	if(usage >= PIO_BUFFER_SIZE_LIMIT){
	  flush_output_buffer(file);
	}
	
    }
#endif


    for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
      if(region == NULL){
	  for(i=0;i<ndims;i++){
	    start[i] = 0;
	    count[i] = 0;
	  }
      }else{
	bufptr = IOBUF+tsize*region->loffset;
	// this is a time dependent multidimensional array
	if(vdesc->record >= 0){
	  start[0] = vdesc->record;
	  count[0] = 1;
	  for(int i=1;i<ndims;i++){
	    start[i] = region->start[i-1];
	    count[i] = region->count[i-1];
	  }
	// Non-time dependent array
	}else{
	  for( i=0;i<ndims;i++){
	    start[i] = region->start[i];
	    count[i] = region->count[i];
	  }
	}
      }

      switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
      case PIO_IOTYPE_NETCDF4P:
	ierr = nc_var_par_access(ncid, vid, NC_COLLECTIVE);
	switch(iodesc->basetype){
	case MPI_DOUBLE:
	case MPI_REAL8:
	  ierr = nc_put_vara_double (ncid, vid,(size_t *) start,(size_t *) count, (const double *) bufptr); 
	  break;
	case MPI_INTEGER:
	  ierr = nc_put_vara_int (ncid, vid, (size_t *) start, (size_t *) count, (const int *) bufptr); 
	  break;
	case MPI_FLOAT:
	case MPI_REAL4:
	  ierr = nc_put_vara_float (ncid, vid, (size_t *) start, (size_t *) count, (const float *) bufptr); 
	  break;
	default:
	  fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	}
      case PIO_IOTYPE_NETCDF4C:
#endif
      case PIO_IOTYPE_NETCDF:
	{
#ifndef _MPISERIAL
	  mpierr = MPI_Type_size(iodesc->basetype, &dsize);
#endif
	  size_t tstart[ndims], tcount[ndims];
	  if(ios->io_rank==0){
	    for(i=0;i<iodesc->num_aiotasks;i++){
	      if(i==0){	    
		for(j=0;j<ndims;j++){
		  tstart[j] =  start[j];
		  tcount[j] =  count[j];
		  tmp_buf = bufptr;
		}
	      }else{
		mpierr = MPI_Send( &ierr, 1, MPI_INT, i, 0, ios->io_comm);  // handshake - tell the sending task I'm ready
		mpierr = MPI_Recv( buflen, 1, MPI_INT, i, 1, ios->io_comm, &status);
		if(buflen>0){
		  mpierr = MPI_Recv( tstart, ndims, MPI_OFFSET, i, ios->num_iotasks+i, ios->io_comm, &status);
		  mpierr = MPI_Recv( tcount, ndims, MPI_OFFSET, i,2*ios->num_iotasks+i, ios->io_comm, &status);
		  tmp_buf = malloc(buflen * dsize);	
		  mpierr = MPI_Recv( tmp_buf, buflen, iodesc->basetype, i, i, ios->io_comm, &status);
		}
	      }
	      
	      if(buflen>0){
		if(iodesc->basetype == MPI_INTEGER){
		  ierr = nc_put_vara_int (ncid, vid, tstart, tcount, (const int *) tmp_buf); 
		}else if(iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8){
		  ierr = nc_put_vara_double (ncid, vid, tstart, tcount, (const double *) tmp_buf); 
		}else if(iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4){
		  ierr = nc_put_vara_float (ncid,vid, tstart, tcount, (const float *) tmp_buf); 
		}else{
		  fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
		}
		if(ierr == PIO_EEDGE){
		  for(i=0;i<ndims;i++)
		    fprintf(stderr,"dim %d start %ld count %ld\n",i,tstart[i],tcount[i]);
		}
		if(tmp_buf != bufptr)
		  free(tmp_buf);
	      }
	    }
	  }else if(ios->io_rank < iodesc->num_aiotasks ){
	    buflen=1;
	    for(i=0;i<ndims;i++){
	      tstart[i] = (size_t) start[i];
	      tcount[i] = (size_t) count[i];
	      buflen*=tcount[i];
	    }
	    mpierr = MPI_Recv( &ierr, 1, MPI_INT, 0, 0, ios->io_comm, &status);  // task0 is ready to recieve
	    mpierr = MPI_Rsend( buflen, 1, MPI_INT, 0, 1, ios->io_comm);
	    if(buflen>0) {
	      mpierr = MPI_Rsend( tstart, ndims, MPI_OFFSET, 0, ios->num_iotasks+ios->io_rank, ios->io_comm);
	      mpierr = MPI_Rsend( tcount, ndims, MPI_OFFSET, 0,2*ios->num_iotasks+ios->io_rank, ios->io_comm);
	      mpierr = MPI_Rsend( bufptr, buflen, iodesc->basetype, 0, ios->io_rank, ios->io_comm);
	    }
	  }
	  break;
	}
#endif
#ifdef _PNETCDF
      case PIO_IOTYPE_PNETCDF:
	for( i=0,dsize=1;i<ndims;i++)
	  dsize*=count[i];
	ierr = ncmpi_bput_vara(ncid, vid,  (PIO_Offset *) start,(PIO_Offset *) count, bufptr,
			       dsize, iodesc->basetype, &request);
	pio_push_request(file,request);
	if(ierr != PIO_NOERR){
	  printf("%s %d usage %ld %d %ld",__FILE__,__LINE__,usage,tsize,iodesc->maxiobuflen);
	for( i=0;i<ndims;i++)
	  printf(" %ld %ld ",start[i],count[i]);
	 printf("dsize %ld\n",dsize);
	}
	break;
#endif
      default:
	ierr = iotype_error(file->iotype,__FILE__,__LINE__);
      }
    
      if(region->next != NULL)
	region = region->next;
    } //    for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
  } // if(ios->ioproc)
  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);
  return ierr;
}

int PIOc_write_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen, void *array, void *fillvalue)
{
  iosystem_desc_t *ios;
  file_desc_t *file;
  io_desc_t *iodesc;
  void *iobuf;
  size_t vsize, rlen;
  int ierr;
  MPI_Datatype vtype;

  ierr = PIO_NOERR;


  file = pio_get_file_from_id(ncid);
  if(file == NULL){
    fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
    return PIO_EBADID;
  }
  iodesc = pio_get_iodesc_from_id(ioid);
  if(iodesc == NULL){
    fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
    return PIO_EBADID;
  }
  iobuf = NULL;

  ios = file->iosystem;

  rlen = iodesc->llen;
  if(iodesc->rearranger>0){
    if(rlen>0){
      vtype = (MPI_Datatype) iodesc->basetype;
      //    printf("rlen = %ld\n",rlen);
      if(vtype == MPI_INTEGER){
	MALLOC_FILL_ARRAY(int, rlen, fillvalue, iobuf);
      }else if(vtype == MPI_FLOAT || vtype == MPI_REAL4){
	MALLOC_FILL_ARRAY(float, rlen, fillvalue, iobuf);
      }else if(vtype == MPI_DOUBLE || vtype == MPI_REAL8){
	MALLOC_FILL_ARRAY(double, rlen, fillvalue, iobuf);
      }else if(vtype == MPI_CHARACTER){
	MALLOC_FILL_ARRAY(char, rlen, fillvalue, iobuf);
      }else{
	fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",vtype);
      }
    }
    //    printf(" rlen = %d %ld\n",rlen,iobuf); 

    //  }


    //    printf("%s %d %ld %d %d %d %ld\n",__FILE__,__LINE__,array,((int *)array)[0],((int *)array)[1],((int *)array)[2], fillvalue);

    ierr = box_rearrange_comp2io(*ios, iodesc, array, iobuf, 0, 0);
  }else{
    iobuf = array;
  }
  switch(file->iotype){
  case PIO_IOTYPE_PBINARY:
    break;
  case PIO_IOTYPE_PNETCDF:
  case PIO_IOTYPE_NETCDF:
  case PIO_IOTYPE_NETCDF4P:
  case PIO_IOTYPE_NETCDF4C:
    ierr = pio_write_darray_nc(file, iodesc, vid, iobuf, fillvalue);
  }

  if(iodesc->rearranger>0 && rlen>0)
    free(iobuf);

  return ierr;

}

int pio_read_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid, void *IOBUF)
{
  int ierr=PIO_NOERR;
  iosystem_desc_t *ios;
  var_desc_t *vdesc;
  int ndims, fndims;
  MPI_Status status;
  int i;
  ios = file->iosystem;
  if(ios == NULL)
    return PIO_EBADID;

  vdesc = (file->varlist)+vid;

  if(vdesc == NULL)
    return PIO_EBADID;

  ndims = iodesc->ndims;
  ierr = PIOc_inq_varndims(file->fh, vid, &fndims);
  
  if(ios->ioproc){
    io_region *region;
    size_t start[fndims];
    size_t count[fndims];
    size_t tmp_start[fndims];
    size_t tmp_count[fndims];
    size_t tmp_bufsize=1;
    int regioncnt;
    void *bufptr;
    int tsize;
    // buffer is incremented by byte and loffset is in terms of the iodessc->basetype
    // so we need to multiply by the size of the basetype
    // We can potentially allow for one iodesc to have multiple datatypes by allowing the
    // calling program to change the basetype.   
    region = iodesc->firstregion;
#ifndef _MPISERIAL
    MPI_Type_size(iodesc->basetype, &tsize);
#endif
    for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
      for(i=0;i<fndims;i++){
	start[i] = 0;
	count[i] = 0;
      }       
      if(region!=NULL){
	bufptr=IOBUF+tsize*region->loffset;

	if(fndims>ndims && vdesc->record >= 0){
	  ndims++;
	  start[0] = vdesc->record;
	  count[0] = 1;
	  for(i=1;i<ndims;i++){
	    start[i] = region->start[i-1];
	    count[i] = region->count[i-1];
	  } 
	}else{
	  // Non-time dependent array
	  for(i=0;i<ndims;i++){
	    start[i] = region->start[i];
	    count[i] = region->count[i];
	  }
	}
      }
    

      switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
      case PIO_IOTYPE_NETCDF4P:
	switch(iodesc->basetype){
	case MPI_DOUBLE:
	case MPI_REAL8:
	  ierr = nc_get_vara_double (file->fh, vid,start,count, bufptr); 
	  break;
	case MPI_INTEGER:
	  ierr = nc_get_vara_int (file->fh, vid, start, count,  bufptr); 
	  break;
	case MPI_FLOAT:
	case MPI_REAL4:
	  ierr = nc_get_vara_float (file->fh, vid, start,  count,  bufptr); 
	  break;
	default:
	  fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	}	
	break;
      case PIO_IOTYPE_NETCDF4C:
#endif
      case PIO_IOTYPE_NETCDF:
	if(ios->io_rank>0){
	  tmp_bufsize=1;
	  for( i=0;i<ndims; i++){
	    tmp_start[i] = start[i];
	    tmp_count[i] = count[i];
	    tmp_bufsize *= count[i];
	  }
	  MPI_Send( tmp_count, ndims, MPI_OFFSET, 0, ios->io_rank, ios->io_comm);
	  if(tmp_bufsize > 0){
	    MPI_Send( tmp_start, ndims, MPI_OFFSET, 0, ios->io_rank, ios->io_comm);
	    MPI_Recv( bufptr, tmp_bufsize, iodesc->basetype, 0, ios->io_rank, ios->io_comm, &status);
	  }
	}else if(ios->io_rank==0){
	  for( i=ios->num_iotasks-1; i>=0; i--){
	    if(i==0){
	      for(int k=0;k<ndims;k++)
		tmp_count[k] = count[k];
	    }else{
	      MPI_Recv(tmp_count, ndims, MPI_OFFSET, i, i, ios->io_comm, &status);
	    }
	    tmp_bufsize=1;
	    for(int j=0;j<ndims; j++){
	      tmp_bufsize *= tmp_count[j];
	    }
	    if(tmp_bufsize>0){
	      if(i==0){
		for(int k=0;k<ndims;k++)
		  tmp_start[k] = start[k];
	      }else{
		MPI_Recv(tmp_start, ndims, MPI_OFFSET, i, i, ios->io_comm, &status);
	      }
	      //	      for(int k=0;k<ndims;k++)
	      //	printf("%d %d %ld %ld %d\n",vid,k,tmp_start[k],tmp_count[k], ierr);
	      
	      if(iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8){
		ierr = nc_get_vara_double (file->fh, vid, tmp_start, tmp_count, bufptr); 
	      }else if(iodesc->basetype == MPI_INTEGER){
		ierr = nc_get_vara_int (file->fh, vid, tmp_start, tmp_count,  bufptr); 	     
	      }else if(iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4){
		ierr = nc_get_vara_float (file->fh, vid, tmp_start, tmp_count,  bufptr); 
	      }else{
		fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	      }	
	      /*
	      if(ierr != PIO_NOERR){
		printf("%s %d ",__FILE__,__LINE__);
		for(int j=0;j<ndims;j++)
		  printf(" %ld %ld",start[j],count[j]);
		printf("\n");
	      }
	      */
	      if(i>0){
		MPI_Rsend(bufptr, tmp_bufsize, iodesc->basetype, i, i, ios->io_comm);
	      }
	    }
	  }
    }
	  break;
#endif
#ifdef _PNETCDF
	case PIO_IOTYPE_PNETCDF:
	  {
	    PIO_Offset tmp_bufsize=1;
	    for(int j=0;j<ndims; j++){
	      tmp_bufsize *= count[j];
	      //	      printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,j,start[j],count[j]); 
	    }
	    // printf("%s %d %ld \n",__FILE__,__LINE__,tmp_bufsize); 
	    
	    ierr = ncmpi_get_vara_all(file->fh, vid,(PIO_Offset *) start,(PIO_Offset *) count, bufptr, tmp_bufsize, iodesc->basetype);
	  }
	  break;
#endif
	default:
	  ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	  
	}
      if(region->next != NULL)
	region = region->next;
    } // for(regioncnt=0;...)
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);
  return ierr;
}

int PIOc_read_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen, void *array)
{
  iosystem_desc_t *ios;
  file_desc_t *file;
  io_desc_t *iodesc;
  void *iobuf=NULL;
  size_t vsize=0, rlen=0;
  int ierr;
  MPI_Datatype vtype;
 

  file = pio_get_file_from_id(ncid);

  if(file == NULL){
    fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
    return PIO_EBADID;
  }
  iodesc = pio_get_iodesc_from_id(ioid);
  if(iodesc == NULL){
    fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
    return PIO_EBADID;
  }
  ios = file->iosystem;
  rlen = iodesc->llen;
  
  if(iodesc->rearranger > 0){
    if(ios->ioproc && rlen>0){
      vtype = (MPI_Datatype) iodesc->basetype;
      if(vtype == MPI_INTEGER){
	iobuf = malloc( rlen*sizeof(int));
      }else if(vtype == MPI_FLOAT || vtype == MPI_REAL4){
      iobuf = malloc( rlen*sizeof(float));
      }else if(vtype == MPI_DOUBLE || vtype == MPI_REAL8){
	iobuf = malloc( rlen*sizeof(double));
      }else if(vtype == MPI_CHARACTER){
	iobuf = malloc( rlen*sizeof(char));
      }else{
	fprintf(stderr,"Type not recognized %d in pioc_read_darray\n",vtype);
      }
      if(iobuf == NULL){
	fprintf(stderr,"malloc failed in pioc_read_darray\n");
	return PIO_ENOMEM;
      } 
    }
  }else{
    iobuf = array;
  }

  switch(file->iotype){
  case PIO_IOTYPE_PBINARY:
    break;
  case PIO_IOTYPE_PNETCDF:
  case PIO_IOTYPE_NETCDF:
  case PIO_IOTYPE_NETCDF4P:
  case PIO_IOTYPE_NETCDF4C:
    ierr = pio_read_darray_nc(file, iodesc, vid, iobuf);
  }
  if(iodesc->rearranger > 0){
    ierr = box_rearrange_io2comp(*ios, iodesc, iobuf, array, 0, 0);
    if(rlen>0)
      free(iobuf);
  }

  return ierr;

}

int flush_output_buffer(file_desc_t *file)
{
  var_desc_t *vardesc;
  int ierr=PIO_NOERR;
#ifdef _PNETCDF
  if(file->nreq==0)
    return ierr;

  if(file->nreq>PIO_MAX_VARS){
    fprintf(stderr,"Need to increase PIO_MAX_VARS %d\n",file->nreq);
  }
  // Since we don't know if the number of requests match across tasks we
  // need to use PIO_MAX_VARS here rather than file->nreq
  int status[file->nreq];
  ierr = ncmpi_wait_all(file->fh,file->nreq,  file->request,status);
  for(int i=0;i<file->nreq;i++){
    file->request[i]=MPI_REQUEST_NULL;
  }
  file->nreq = 0;

#endif
  return ierr;
}

