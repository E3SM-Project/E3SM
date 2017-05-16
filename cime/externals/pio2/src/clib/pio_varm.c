#include <config.h>
#include <pio.h>
#include <pio_internal.h>

///
/// PIO interface to nc_put_varm
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;

  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm(file->fh, varid, start, count, stride, imap, buf, bufcount, buftype, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_uchar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_UCHAR;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;

  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_uchar(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_short
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_SHORT;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;

  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_short(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
///
/// PIO interface to nc_put_varm_text
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_TEXT;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_text(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_ushort
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_USHORT;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_ushort(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_ulonglong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_ULONGLONG;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_ulonglong(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
///
/// PIO interface to nc_put_varm_int
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_INT;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_int(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_float
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const float *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_FLOAT;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_float(file->fh, varid,(size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_float(file->fh, varid,(size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_float(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
///
/// PIO interface to nc_put_varm_long
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_LONG;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_long(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_uint
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_UINT;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_uint(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_double
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const double *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_DOUBLE;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_double(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}
///
/// PIO interface to nc_put_varm_schar
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const signed char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_SCHAR;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;



  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_schar(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

///
/// PIO interface to nc_put_varm_longlong
///
/// This routine is called collectively by all tasks in the communicator ios.union_comm.  
/// 
/// Refer to the <A HREF="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html"> netcdf documentation. </A>
///
int PIOc_put_varm_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], const long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  var_desc_t *vdesc;
  PIO_Offset usage;
  int *request;

  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_VARM_LONGLONG;

  /* Sorry, but varm functions are not supported by the async interface. */
  if(ios->async_interface)
      return PIO_EINVAL;

  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_var_par_access(file->fh, varid, NC_COLLECTIVE);
      ierr = nc_put_varm_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_varm_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      vdesc = file->varlist + varid;

      if(vdesc->nreqs%PIO_REQUEST_ALLOC_CHUNK == 0 ){
	vdesc->request = realloc(vdesc->request, 
				 sizeof(int)*(vdesc->nreqs+PIO_REQUEST_ALLOC_CHUNK));
      }
      request = vdesc->request+vdesc->nreqs;

      if(ios->io_rank==0){
	ierr = ncmpi_bput_varm_longlong(file->fh, varid, start, count, stride, imap, op, request);;
      }else{
	*request = PIO_REQ_NULL;
      }
      vdesc->nreqs++;
      flush_output_buffer(file, false, 0);
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_get_varm_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned char *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_UCHAR;
  ibuftype = MPI_UNSIGNED_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_uchar(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_uchar_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], signed char *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_SCHAR;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_schar(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_schar_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], double *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_DOUBLE;
  ibuftype = MPI_DOUBLE;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_double(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_double(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_double_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], char *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_TEXT;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_text(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_text_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], int *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_INT;
  ibuftype = MPI_INT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_int(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_int_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned int *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_UINT;
  ibuftype = MPI_UNSIGNED;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_uint(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_uint(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_uint_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM;
  ibufcnt = bufcount;
  ibuftype = buftype;
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm(file->fh, varid, start, count, stride, imap, buf, bufcount, buftype);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_all(file->fh, varid, start, count, stride, imap, buf, bufcount, buftype);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], float *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_FLOAT;
  ibuftype = MPI_FLOAT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_float(file->fh, varid,(size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_float(file->fh, varid,(size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_float(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_float_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], long *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_LONG;
  ibuftype = MPI_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_long(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_long_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned short *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_USHORT;
  ibuftype = MPI_UNSIGNED_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_ushort(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_ushort(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_ushort_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], long long *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_LONGLONG;
  ibuftype = MPI_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_longlong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_longlong(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_longlong_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], short *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_SHORT;
  ibuftype = MPI_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_short(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_short(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_short_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

int PIOc_get_varm_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], const PIO_Offset imap[], unsigned long long *buf) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  MPI_Datatype ibuftype;
  int ndims;
  int ibufcnt;
  bool bcast = false;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_VARM_ULONGLONG;
  ibuftype = MPI_UNSIGNED_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i]/stride[i];
  }
  ierr = PIO_NOERR;

  if(ios->async_interface && ! ios->ioproc){
    if(ios->compmaster) 
      mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
    mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, 0, ios->intercomm);
  }


  if(ios->ioproc){
    switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
    case PIO_IOTYPE_NETCDF4P:
      ierr = nc_get_varm_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->iomaster){
	ierr = nc_get_varm_ulonglong(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
#ifdef PNET_READ_AND_BCAST
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_varm_ulonglong(file->fh, varid, start, count, stride, imap,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
#else
      ierr = ncmpi_get_varm_ulonglong_all(file->fh, varid, start, count, stride, imap,  buf);;
#endif
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  if(ios->async_interface || bcast ||  
     (ios->num_iotasks < ios->num_comptasks)){
    MPI_Bcast(buf, ibufcnt, ibuftype, ios->ioroot, ios->my_comm);
  }

  return ierr;
}

