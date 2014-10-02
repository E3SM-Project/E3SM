#include <pio.h>
#include <pio_internal.h>

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_schar (int ncid, int varid, const PIO_Offset index[], signed char *buf) 
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
  msg = PIO_MSG_GET_VAR1_SCHAR;
  ibuftype = MPI_CHAR;
  ibufcnt = 1;
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
      ierr = nc_get_var1_schar(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_schar(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_schar(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned long long *buf) 
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
  msg = PIO_MSG_GET_VARS_ULONGLONG;
  ibuftype = MPI_UNSIGNED_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_ulonglong_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_uchar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_uchar_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_schar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_schar(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_schar_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], short *buf)  
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
  msg = PIO_MSG_GET_VARS_SHORT;
  ibuftype = MPI_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_short_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_double (int ncid, int varid, double *buf) 
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
  msg = PIO_MSG_GET_VAR_DOUBLE;
  ibuftype = MPI_DOUBLE;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_double(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_double(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_double_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], double *buf) 
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
  msg = PIO_MSG_GET_VARA_DOUBLE;
  ibuftype = MPI_DOUBLE;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_double(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_double(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_double_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_int (int ncid, int varid, int *buf) 
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
  msg = PIO_MSG_GET_VAR_INT;
  ibuftype = MPI_INT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_int(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_int(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_int_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_ushort (int ncid, int varid, unsigned short *buf) 
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
  msg = PIO_MSG_GET_VAR_USHORT;
  ibuftype = MPI_UNSIGNED_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_ushort(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_ushort(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_ushort_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], char *buf)  
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
  msg = PIO_MSG_GET_VARA_TEXT;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_text(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_text(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_text_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], int *buf)  
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
  msg = PIO_MSG_GET_VARA_INT;
  ibuftype = MPI_INT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_int(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_int(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_int_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_float (int ncid, int varid, const PIO_Offset index[], float *buf) 
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
  msg = PIO_MSG_GET_VAR1_FLOAT;
  ibuftype = MPI_FLOAT;
  ibufcnt = 1;
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
      ierr = nc_get_var1_float(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_float(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_float(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_short (int ncid, int varid, const PIO_Offset index[], short *buf) 
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
  msg = PIO_MSG_GET_VAR1_SHORT;
  ibuftype = MPI_SHORT;
  ibufcnt = 1;
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
      ierr = nc_get_var1_short(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_short(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_short(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_int (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], int *buf)  
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
  msg = PIO_MSG_GET_VARS_INT;
  ibuftype = MPI_INT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_int(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_int_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_text (int ncid, int varid, char *buf) 
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
  msg = PIO_MSG_GET_VAR_TEXT;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_text(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_text(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_text_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_double_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], signed char *buf)  
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
  msg = PIO_MSG_GET_VARS_SCHAR;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_schar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_schar_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned short *buf) 
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
  msg = PIO_MSG_GET_VARA_USHORT;
  ibuftype = MPI_UNSIGNED_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_ushort(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_ushort(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_ushort_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_ushort (int ncid, int varid, const PIO_Offset index[], unsigned short *buf) 
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
  msg = PIO_MSG_GET_VAR1_USHORT;
  ibuftype = MPI_UNSIGNED_SHORT;
  ibufcnt = 1;
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
      ierr = nc_get_var1_ushort(file->fh, varid, (size_t *) index,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_ushort(file->fh, varid, (size_t *) index,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_ushort(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_float (int ncid, int varid, float *buf) 
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
  msg = PIO_MSG_GET_VAR_FLOAT;
  ibuftype = MPI_FLOAT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_float(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_float(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_float_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned char *buf) 
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
  msg = PIO_MSG_GET_VARS_UCHAR;
  ibuftype = MPI_UNSIGNED_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_uchar(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_uchar_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var (int ncid, int varid, void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
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
  msg = PIO_MSG_GET_VAR;
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
      ierr = nc_get_var(file->fh, varid,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var(file->fh, varid,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_all(file->fh, varid, buf, bufcount, buftype);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_longlong (int ncid, int varid, const PIO_Offset index[], long long *buf) 
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
  msg = PIO_MSG_GET_VAR1_LONGLONG;
  ibuftype = MPI_LONG_LONG;
  ibufcnt = 1;
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
      ierr = nc_get_var1_longlong(file->fh, varid, (size_t *) index,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_longlong(file->fh, varid, (size_t *) index,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_longlong(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_ushort (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned short *buf) 
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
  msg = PIO_MSG_GET_VARS_USHORT;
  ibuftype = MPI_UNSIGNED_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_ushort_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_long (int ncid, int varid, long *buf) 
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
  msg = PIO_MSG_GET_VAR_LONG;
  ibuftype = MPI_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_long(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_long(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_long_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_double (int ncid, int varid, const PIO_Offset index[], double *buf) 
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
  msg = PIO_MSG_GET_VAR1_DOUBLE;
  ibuftype = MPI_DOUBLE;
  ibufcnt = 1;
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
      ierr = nc_get_var1_double(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_double(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_double(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned int *buf) 
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
  msg = PIO_MSG_GET_VARA_UINT;
  ibuftype = MPI_UNSIGNED;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_uint(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_uint(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_uint_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], long long *buf) 
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
  msg = PIO_MSG_GET_VARS_LONGLONG;
  ibuftype = MPI_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_longlong_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_longlong (int ncid, int varid, long long *buf) 
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
  msg = PIO_MSG_GET_VAR_LONGLONG;
  ibuftype = MPI_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_longlong(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_longlong(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_longlong_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_short (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], short *buf)  
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
  msg = PIO_MSG_GET_VARA_SHORT;
  ibuftype = MPI_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_short(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_short(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_short_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], long *buf) 
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
  msg = PIO_MSG_GET_VARA_LONG;
  ibuftype = MPI_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_long(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_long(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_long_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_int (int ncid, int varid, const PIO_Offset index[], int *buf) 
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
  msg = PIO_MSG_GET_VAR1_INT;
  ibuftype = MPI_INT;
  ibufcnt = 1;
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
      ierr = nc_get_var1_int(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_int(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_int(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_ulonglong (int ncid, int varid, const PIO_Offset index[], unsigned long long *buf) 
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
  msg = PIO_MSG_GET_VAR1_ULONGLONG;
  ibuftype = MPI_UNSIGNED_LONG_LONG;
  ibufcnt = 1;
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
      ierr = nc_get_var1_ulonglong(file->fh, varid, (size_t *) index,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_ulonglong(file->fh, varid, (size_t *) index,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_ulonglong(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_uchar (int ncid, int varid, unsigned char *buf) 
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
  msg = PIO_MSG_GET_VAR_UCHAR;
  ibuftype = MPI_UNSIGNED_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_uchar(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_uchar(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_uchar_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_uchar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned char *buf) 
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
  msg = PIO_MSG_GET_VARA_UCHAR;
  ibuftype = MPI_UNSIGNED_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_uchar(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_uchar(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_uchar_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], float *buf)  
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
  msg = PIO_MSG_GET_VARS_FLOAT;
  ibuftype = MPI_FLOAT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_float(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_float_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_long (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], long *buf) 
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
  msg = PIO_MSG_GET_VARS_LONG;
  ibuftype = MPI_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_long(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_long_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1 (int ncid, int varid, const PIO_Offset index[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
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
  msg = PIO_MSG_GET_VAR1;
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
      ierr = nc_get_var1(file->fh, varid, (size_t *) index,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1(file->fh, varid, (size_t *) index,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1(file->fh, varid, index, buf, bufcount, buftype);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_uint (int ncid, int varid, unsigned int *buf) 
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
  msg = PIO_MSG_GET_VAR_UINT;
  ibuftype = MPI_UNSIGNED;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_uint(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_uint(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_uint_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype) 
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
  msg = PIO_MSG_GET_VARA;
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
      ierr = nc_get_vara(file->fh, varid, (size_t *) start,  (size_t *) count,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara(file->fh, varid, (size_t *) start,  (size_t *) count,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_all(file->fh, varid, start, count, buf, bufcount, buftype);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_schar (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], signed char *buf) 
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
  msg = PIO_MSG_GET_VARA_SCHAR;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_schar(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_schar(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_schar_all(file->fh, varid, start, count,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_uint (int ncid, int varid, const PIO_Offset index[], unsigned int *buf) 
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
  msg = PIO_MSG_GET_VAR1_UINT;
  ibuftype = MPI_UNSIGNED;
  ibufcnt = 1;
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
      ierr = nc_get_var1_uint(file->fh, varid, (size_t *) index,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_uint(file->fh, varid, (size_t *) index,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_uint(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_uint (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], unsigned int *buf) 
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
  msg = PIO_MSG_GET_VARS_UINT;
  ibuftype = MPI_UNSIGNED;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_uint_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_float (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], float *buf)  
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
  msg = PIO_MSG_GET_VARA_FLOAT;
  ibuftype = MPI_FLOAT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_float(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_float(file->fh, varid, (size_t *) start, (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_float_all(file->fh, varid, start, count,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_text(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_text(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_text_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_text (int ncid, int varid, const PIO_Offset index[], char *buf) 
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
  msg = PIO_MSG_GET_VAR1_TEXT;
  ibuftype = MPI_CHAR;
  ibufcnt = 1;
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
      ierr = nc_get_var1_text(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_text(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_text(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_int(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_int(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_int_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_uint(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_uint_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
      ierr = nc_get_varm(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_all(file->fh, varid, start, count, stride, imap, buf, bufcount, buftype); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_double (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], double *buf) 
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
  msg = PIO_MSG_GET_VARS_DOUBLE;
  ibuftype = MPI_DOUBLE;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_double(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_double_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_longlong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], long long *buf) 
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
  msg = PIO_MSG_GET_VARA_LONGLONG;
  ibuftype = MPI_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_longlong(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_longlong(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_longlong_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_ulonglong (int ncid, int varid, unsigned long long *buf) 
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
  msg = PIO_MSG_GET_VAR_ULONGLONG;
  ibuftype = MPI_UNSIGNED_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_ulonglong(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_ulonglong(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_ulonglong_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vara_ulonglong (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], unsigned long long *buf) 
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
  msg = PIO_MSG_GET_VARA_ULONGLONG;
  ibuftype = MPI_UNSIGNED_LONG_LONG;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
  ibufcnt = 1;
  for(int i=0;i<ndims;i++){
    ibufcnt *= count[i];
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
      ierr = nc_get_vara_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vara_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vara_ulonglong_all(file->fh, varid, start, count,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_short (int ncid, int varid, short *buf) 
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
  msg = PIO_MSG_GET_VAR_SHORT;
  ibuftype = MPI_SHORT;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_short(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_short(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_short_all(file->fh, varid,  buf);;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_float(file->fh, varid,(size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_float(file->fh, varid,(size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_float_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_long (int ncid, int varid, const PIO_Offset index[], long *buf) 
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
  msg = PIO_MSG_GET_VAR1_LONG;
  ibuftype = MPI_LONG;
  ibufcnt = 1;
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
      ierr = nc_get_var1_long(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_long(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_long(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_long(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_long(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_long_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_ushort(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_ushort_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_longlong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_longlong_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars_text (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], char *buf)  
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
  msg = PIO_MSG_GET_VARS_TEXT;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars_text(file->fh, varid, (size_t *) start, (size_t *) count, (ptrdiff_t *) stride,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_text_all(file->fh, varid, start, count, stride,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var1_uchar (int ncid, int varid, const PIO_Offset index[], unsigned char *buf) 
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
  msg = PIO_MSG_GET_VAR1_UCHAR;
  ibuftype = MPI_UNSIGNED_CHAR;
  ibufcnt = 1;
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
      ierr = nc_get_var1_uchar(file->fh, varid, (size_t *) index,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var1_uchar(file->fh, varid, (size_t *) index,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
        ncmpi_begin_indep_data(file->fh);
        if(ios->iomaster){
      ierr = ncmpi_get_var1_uchar(file->fh, varid, index,  buf);;
        };
         ncmpi_end_indep_data(file->fh);
        bcast=true;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_vars (int ncid, int varid, const PIO_Offset start[], const PIO_Offset count[], const PIO_Offset stride[], void *buf, PIO_Offset bufcount, MPI_Datatype buftype)  
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
  msg = PIO_MSG_GET_VARS;
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
      ierr = nc_get_vars(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_vars(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,   buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_vars_all(file->fh, varid, start, count, stride, buf, bufcount, buftype); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_short(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride, (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_short_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
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
#ifdef READ_AND_BCAST
      if(ios->io_rank>0){
	   for(int idim=0;idim<ndims;idim++)
     	  count[idim]=0;
      }
      bcast=true;
#endif
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
      ierr = nc_get_varm_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_varm_ulonglong(file->fh, varid, (size_t *) start,  (size_t *) count, (ptrdiff_t *) stride,  (ptrdiff_t *) imap,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_varm_ulonglong_all(file->fh, varid, start, count, stride, imap,  buf); ;
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

#ifdef BGQ
#define READ_AND_BCAST 1
#endif
int PIOc_get_var_schar (int ncid, int varid, signed char *buf) 
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
  msg = PIO_MSG_GET_VAR_SCHAR;
  ibuftype = MPI_CHAR;
  ierr = PIOc_inq_varndims(file->fh, varid, &ndims);
  int dimid[ndims];
  PIO_Offset dimsize;
  ibufcnt = 1;
  PIOc_inq_vardimid(file->fh, varid, dimid);
  for(int i=0;i<ndims;i++){
    PIOc_inq_dimlen(file->fh, dimid[i], &dimsize);
    ibufcnt *= dimsize;
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
      ierr = nc_get_var_schar(file->fh, varid,  buf);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      bcast = true;
      if(ios->io_rank==0){
	ierr = nc_get_var_schar(file->fh, varid,  buf);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_var_schar_all(file->fh, varid,  buf);;
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

