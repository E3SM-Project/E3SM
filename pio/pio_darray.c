#include <pio.h>
#include <pio_internal.h>


#define MALLOC_FILL_ARRAY(type, n, fill, arr) \
  arr = malloc(n * sizeof (type));	      \
  if(fill != NULL) \
    for(int _i=0; _i<n; _i++)			\
      ((type *) arr)[_i] = *((type *) fill)


int c_write_nc(const file_desc_t *file,const io_desc_t *iodesc, void *IOBUF, const int varid, const PIO_Offset start[], const PIO_Offset count[], MPI_Request *request)
{
  int ierr;
  int msg;
  int mpierr;
  int chkerr;
  iosystem_desc_t *ios;
  int dlen;
  int i;
  void *tmp_buf;
  int dsize;
  MPI_Status status;


  ios = file->iosystem;
  if(ios == NULL)
    return PIO_EBADID;
  msg = 0;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->comp_rank==0) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
  }


  if(ios->ioproc){
    int ndims = iodesc->ndims;
    int ncid = file->fh;
    size_t tstart[ndims], tcount[ndims];
    
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(ncid, varid, NC_COLLECTIVE);
      switch(iodesc->basetype){
      case MPI_DOUBLE:
      case MPI_REAL8:
	ierr = nc_put_vara_double (ncid, varid,(size_t *) start,(size_t *) count, (const double *) IOBUF); 
	break;
      case MPI_INTEGER:
	ierr = nc_put_vara_int (ncid, varid, (size_t *) start, (size_t *) count, (const int *) IOBUF); 
	break;
      case MPI_FLOAT:
      case MPI_REAL4:
	ierr = nc_put_vara_float (ncid, varid, (size_t *) start, (size_t *) count, (const float *) IOBUF); 
	break;
    default:
      fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
    }
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      mpierr = MPI_Type_size(iodesc->basetype, &dsize);
      if(ios->io_rank==0){
	for(i=0;i<ios->num_iotasks;i++){
	  if(i==0){	    
	    for(int i=0;i<ndims;i++){
	      tstart[i] = (int) start[i];
	      tcount[i] = (int) count[i];
	      tmp_buf = IOBUF;
	    }
	  }else{
	    if(i==1){
	      tmp_buf = malloc(iodesc->maxiobuflen * dsize);	
	    }
	    
	    mpierr = MPI_Send( &ierr, 1, MPI_INT, i, 0, ios->io_comm);  // handshake - tell the sending task I'm ready
	    mpierr = MPI_Recv( tstart, ndims, MPI_INT, i, ios->num_iotasks+i, ios->io_comm, &status);
	    mpierr = MPI_Recv( tcount, ndims, MPI_INT, i,2*ios->num_iotasks+i, ios->io_comm, &status);
	    mpierr = MPI_Recv( tmp_buf, iodesc->maxiobuflen, iodesc->basetype, i, i, ios->io_comm, &status);
	  }
	  if(iodesc->basetype == MPI_INTEGER){
	    ierr = nc_put_vara_int (ncid, varid, tstart, tcount, (const int *) tmp_buf); 
	  }else if(iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8){
	    ierr = nc_put_vara_double (ncid, varid, tstart, tcount, (const double *) tmp_buf); 
	  }else if(iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4){
	    ierr = nc_put_vara_float (ncid,varid, tstart, tcount, (const float *) tmp_buf); 
	  }else{
	    fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	  }
	}       
	free(tmp_buf);
      }else{
	for(i=0;i<ndims;i++){
	  tstart[i] = (size_t) start[i];
	  tcount[i] = (size_t) count[i];
	}
	mpierr = MPI_Recv( &ierr, 1, MPI_INT, 0, 0, ios->io_comm, &status);  // task0 is ready to recieve
	mpierr = MPI_Send( tstart, ndims, MPI_INT, 0, ios->num_iotasks+ios->io_rank, ios->io_comm);
	mpierr = MPI_Send( tcount, ndims, MPI_INT, 0,2*ios->num_iotasks+ios->io_rank, ios->io_comm);
	mpierr = MPI_Send( IOBUF, iodesc->maxiobuflen, iodesc->basetype, 0, ios->io_rank, ios->io_comm);

      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      for( i=0,dsize=1;i<ndims;i++)
	dsize*=count[i];
      // assert(dsize > 0)
      ierr = ncmpi_iput_vara(ncid, varid, (size_t *) start,(size_t *) count, IOBUF,
			     dsize, iodesc->basetype, request);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  chkerr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}







int pio_write_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid, const int vtype, 
			const PIO_Offset arraylen, void *array, void *fillvalue)
{
  iosystem_desc_t *ios;
  PIO_Offset *start, *count;
  var_desc_t *vdesc;
  int ndims;
  int ierr;

  ios = file->iosystem;
  vdesc = (file->varlist)+vid;
  ndims = iodesc->ndims;
  if(ios->ioproc){
    // this is a time dependent multidimensional array
    if(vdesc->record >= 0 && iodesc->start[1]>=0){
      start = (PIO_Offset *) calloc(ndims+1, sizeof(PIO_Offset));
      count = (PIO_Offset *) calloc(ndims+1, sizeof(PIO_Offset));
      for(int i=0;i<ndims;i++){
	start[i] = iodesc->start[i];
	count[i] = iodesc->count[i];
	start[ndims] = vdesc->record;
	count[ndims] = 1;
      }
      // this is a timedependent single value
    }else if(vdesc->record >= 0){
      start = (PIO_Offset *) malloc( sizeof(PIO_Offset));
      count = (PIO_Offset *) malloc( sizeof(PIO_Offset));
      start[0] = vdesc->record;
      count[0] = 1;
      // Non-time dependent array
    }else{
      start = (PIO_Offset *) calloc(ndims, sizeof(PIO_Offset));
      count = (PIO_Offset *) calloc(ndims, sizeof(PIO_Offset));
    }      

  }

  ierr = c_write_nc(file, iodesc, array, vid, start, count, &(vdesc->request));


  return PIO_NOERR;
}

int PIOc_write_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen, void *array, void *fillvalue)
{
  iosystem_desc_t *ios;
  file_desc_t *file;
  io_desc_t *iodesc;
  void *iobuf=NULL;
  size_t vsize, rlen;
  int ierr;
  MPI_Datatype vtype;

  file = pio_get_file_from_id(ncid);

  if(file == NULL)
    return PIO_EBADID;

  iodesc = pio_get_iodesc_from_id(ioid);
  if(iodesc == NULL)
    return PIO_EBADID;

  ios = file->iosystem;
  if(ios->ioproc){
    rlen = iodesc->llen;

    vtype = (MPI_Datatype) iodesc->basetype;
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
  ierr = box_rearrange_comp2io(ios, iodesc, arraylen, array, rlen, iobuf, 0, 0);


  
  switch(file->iotype){
  case PIO_IOTYPE_PBINARY:
    break;
  case PIO_IOTYPE_PNETCDF:
  case PIO_IOTYPE_NETCDF:
  case PIO_IOTYPE_NETCDF4P:
  case PIO_IOTYPE_NETCDF4C:
    ierr = pio_write_darray_nc(file, iodesc, vid, vtype, arraylen, iobuf, fillvalue);
  }
  return PIO_NOERR;

}

