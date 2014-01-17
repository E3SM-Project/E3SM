#include <pio.h>
#include <pio_internal.h>


#define MALLOC_FILL_ARRAY(type, n, fill, arr) \
  arr = malloc(n * sizeof (type));	      \
  if(fill != NULL) \
    for(int i=0; i<n; i++)			\
      ((type *) arr)[i] = *((type *) fill)

int pio_write_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid, const int vtype, 
			const PIO_Offset arraylen, void *array, void *fillvalue)
{
  iosystem_desc_t *ios;
  PIO_Offset *start, *count;
  var_desc_t *vdesc;
  int ndims;

  ios = file->iosystem;
  vdesc = (file->vardesc)+vid;
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
    switch(vtype){
    case MPI_INTEGER:
      MALLOC_FILL_ARRAY(int, rlen, fillvalue, iobuf);
      break;
    case MPI_REAL4:
      MALLOC_FILL_ARRAY(float, rlen, fillvalue, iobuf);
      break;
    case MPI_REAL8:
      MALLOC_FILL_ARRAY(double, rlen, fillvalue, iobuf);
      break;
    case MPI_CHARACTER:
      MALLOC_FILL_ARRAY(char, rlen, fillvalue, iobuf);
      break;
    default:
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

