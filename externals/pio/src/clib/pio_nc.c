/**
* @file   pio_nc.c
* @author Jim Edwards (jedwards@ucar.edu)
* @date     Feburary 2014 
* @brief    PIO interfaces to [NetCDF](http://www.unidata.ucar.edu/software/netcdf/docs/modules.html) support functions
* @details
*  This file provides an interface to the [NetCDF](http://www.unidata.ucar.edu/software/netcdf/docs/modules.html) support functions.
*  Each subroutine calls the underlying netcdf or pnetcdf or netcdf4 functions from 
*  the appropriate subset of mpi tasks (io_comm). Each routine must be called 
*  collectively from union_comm.
*  
*/
#include <pio.h>
#include <pio_internal.h>

/** 
 * @name    PIOc_inq_att
 * @brief   The PIO-C interface for the NetCDF function nc_inq_att.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_inq_att (int ncid, int varid, const char *name, nc_type *xtypep, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_ATT;

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
      ierr = nc_inq_att(file->fh, varid, name, xtypep, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_att(file->fh, varid, name, xtypep, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_att(file->fh, varid, name, xtypep, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(xtypep != NULL) mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(lenp != NULL) mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_format
 * @brief   The PIO-C interface for the NetCDF function nc_inq_format.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 */
int PIOc_inq_format (int ncid, int *formatp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_FORMAT;

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
      ierr = nc_inq_format(file->fh, formatp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_format(file->fh, formatp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_format(file->fh, formatp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(formatp , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_varid
 * @brief   The PIO-C interface for the NetCDF function nc_inq_varid.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_varid (int ncid, const char *name, int *varidp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VARID;

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
      ierr = nc_inq_varid(file->fh, name, varidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varid(file->fh, name, varidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varid(file->fh, name, varidp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(varidp,1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_varnatts
 * @brief   The PIO-C interface for the NetCDF function nc_inq_varnatts.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_varnatts (int ncid, int varid, int *nattsp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VARNATTS;

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
      ierr = nc_inq_varnatts(file->fh, varid, nattsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varnatts(file->fh, varid, nattsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varnatts(file->fh, varid, nattsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(nattsp,1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_def_var
 * @brief   The PIO-C interface for the NetCDF function nc_def_var.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_def_var (int ncid, const char *name, nc_type xtype, int ndims, const int *dimidsp, int *varidp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_DEF_VAR;

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
      ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
      if(ios->io_rank==0){
        ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);
        if(ierr == PIO_NOERR){
          ierr = nc_def_var_deflate(file->fh, *varidp, 0,1,1);
        }
      }
      break;
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_def_var(file->fh, name, xtype, ndims, dimidsp, varidp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(varidp , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_var
 * @brief   The PIO-C interface for the NetCDF function nc_inq_var.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_var (int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp, int *dimidsp, int *nattsp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VAR;

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
      ierr = nc_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_var(file->fh, varid, name, xtypep, ndimsp, dimidsp, nattsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(xtypep != NULL) mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);
    if(ndimsp != NULL){ mpierr = MPI_Bcast(ndimsp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);
      file->varlist[varid].ndims = (*ndimsp);}
      if(nattsp != NULL) mpierr = MPI_Bcast(nattsp,1, MPI_INT, ios->ioroot, ios->my_comm);
    if(name != NULL){ 
      int slen;
      if(ios->iomaster)
        slen = (int) strlen(name);
      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);
      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);
    }
    if(dimidsp != NULL) {int ndims;
      PIOc_inq_varndims(file->fh, varid, &ndims);
      mpierr = MPI_Bcast(dimidsp , ndims, MPI_INT, ios->ioroot, ios->my_comm);
     }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_varname
 * @brief   The PIO-C interface for the NetCDF function nc_inq_varname.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_varname (int ncid, int varid, char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VARNAME;

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
      ierr = nc_inq_varname(file->fh, varid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varname(file->fh, varid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varname(file->fh, varid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(name != NULL){
      int slen;
      if(ios->iomaster)
        slen = (int) strlen(name);
      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);
      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);
    }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_double
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_double.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_double (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const double *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_DOUBLE;

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
      ierr = nc_put_att_double(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_double(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_double(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_int
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_int.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_int (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_INT;

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
      ierr = nc_put_att_int(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_int(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_int(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_rename_att
 * @brief   The PIO-C interface for the NetCDF function nc_rename_att.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_rename_att (int ncid, int varid, const char *name, const char *newname) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_RENAME_ATT;

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
      ierr = nc_rename_att(file->fh, varid, name, newname);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_att(file->fh, varid, name, newname);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_att(file->fh, varid, name, newname);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_ubyte
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_ubyte.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_ubyte (int ncid, int varid, const char *name, unsigned char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_UBYTE;

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
      ierr = nc_get_att_ubyte(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_ubyte(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_ubyte(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_BYTE, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_natts
 * @brief   The PIO-C interface for the NetCDF function nc_inq_natts.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_inq_natts (int ncid, int *ngattsp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_NATTS;

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
      ierr = nc_inq_natts(file->fh, ngattsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_natts(file->fh, ngattsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_natts(file->fh, ngattsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(ngattsp,1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_del_att
 * @brief   The PIO-C interface for the NetCDF function nc_del_att.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_del_att (int ncid, int varid, const char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_DEL_ATT;

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
      ierr = nc_del_att(file->fh, varid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_del_att(file->fh, varid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_del_att(file->fh, varid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq
 * @brief   The PIO-C interface for the NetCDF function nc_inq.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 */
int PIOc_inq (int ncid, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ;

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
      ierr = nc_inq(file->fh, ndimsp, nvarsp, ngattsp, unlimdimidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq(file->fh, ndimsp, nvarsp, ngattsp, unlimdimidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq(file->fh, ndimsp, nvarsp, ngattsp, unlimdimidp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ndimsp != NULL)
        mpierr = MPI_Bcast(ndimsp, 1, MPI_INT, ios->ioroot, ios->my_comm);
       if(nvarsp != NULL)
        mpierr = MPI_Bcast(nvarsp, 1, MPI_INT, ios->ioroot, ios->my_comm);
       if(ngattsp != NULL)
        mpierr = MPI_Bcast(ngattsp, 1, MPI_INT, ios->ioroot, ios->my_comm);
       if(unlimdimidp != NULL)
        mpierr = MPI_Bcast(unlimdimidp, 1, MPI_INT, ios->ioroot, ios->my_comm);
   if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_text
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_text.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_text (int ncid, int varid, const char *name, char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_TEXT;

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
      ierr = nc_get_att_text(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_text(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_text(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_CHAR, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_short
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_short.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_short (int ncid, int varid, const char *name, short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_SHORT;

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
      ierr = nc_get_att_short(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_short(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_short(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_SHORT, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_long
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_long.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_long (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_LONG;

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
      ierr = nc_put_att_long(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_long(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_long(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_redef
 * @brief   The PIO-C interface for the NetCDF function nc_redef.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 */
int PIOc_redef (int ncid) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_REDEF;

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
      ierr = nc_redef(file->fh);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_redef(file->fh);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_redef(file->fh);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_set_fill
 * @brief   The PIO-C interface for the NetCDF function nc_set_fill.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 */
int PIOc_set_fill (int ncid, int fillmode, int *old_modep) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_SET_FILL;

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
      ierr = nc_set_fill(file->fh, fillmode, old_modep);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_set_fill(file->fh, fillmode, old_modep);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_set_fill(file->fh, fillmode, old_modep);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_enddef
 * @brief   The PIO-C interface for the NetCDF function nc_enddef.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html
 */
int PIOc_enddef (int ncid) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_ENDDEF;

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
      ierr = nc_enddef(file->fh);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_enddef(file->fh);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_enddef(file->fh);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_rename_var
 * @brief   The PIO-C interface for the NetCDF function nc_rename_var.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_rename_var (int ncid, int varid, const char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_RENAME_VAR;

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
      ierr = nc_rename_var(file->fh, varid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_var(file->fh, varid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_var(file->fh, varid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_short
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_short.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_short (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_SHORT;

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
      ierr = nc_put_att_short(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_short(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_short(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_text
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_text.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_text (int ncid, int varid, const char *name, PIO_Offset len, const char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_TEXT;

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
      ierr = nc_put_att_text(file->fh, varid, name, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_text(file->fh, varid, name, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_text(file->fh, varid, name, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_attname
 * @brief   The PIO-C interface for the NetCDF function nc_inq_attname.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_inq_attname (int ncid, int varid, int attnum, char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_ATTNAME;

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
      ierr = nc_inq_attname(file->fh, varid, attnum, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attname(file->fh, varid, attnum, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attname(file->fh, varid, attnum, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(name != NULL){
      int slen;
      if(ios->iomaster)
        slen = (int) strlen(name);
      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);
      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);
    }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_ulonglong
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_ulonglong.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_ulonglong (int ncid, int varid, const char *name, unsigned long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_ULONGLONG;

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
      ierr = nc_get_att_ulonglong(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_ulonglong(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_ulonglong(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_UNSIGNED_LONG_LONG, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_ushort
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_ushort.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_ushort (int ncid, int varid, const char *name, unsigned short *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_USHORT;

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
      ierr = nc_get_att_ushort(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_ushort(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_ushort(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_UNSIGNED_SHORT, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_ulonglong
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_ulonglong.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_ulonglong (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_ULONGLONG;

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
      ierr = nc_put_att_ulonglong(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_ulonglong(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_ulonglong(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_dimlen
 * @brief   The PIO-C interface for the NetCDF function nc_inq_dimlen.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_inq_dimlen (int ncid, int dimid, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_DIMLEN;

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
      ierr = nc_inq_dimlen(file->fh, dimid, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimlen(file->fh, dimid, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimlen(file->fh, dimid, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_var_fill
 * @brief   The PIO-C interface for the NetCDF function nc_inq_var_fill.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_var_fill (int ncid, int varid, int *no_fill, void *fill_value) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VAR_FILL;

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
      ierr = nc_inq_var_fill(file->fh, varid, no_fill, fill_value);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_var_fill(file->fh, varid, no_fill, fill_value);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_var_fill(file->fh, varid, no_fill, fill_value);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(fill_value, 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_uint
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_uint.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_uint (int ncid, int varid, const char *name, unsigned int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_UINT;

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
      ierr = nc_get_att_uint(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_uint(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_uint(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_UNSIGNED, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_longlong
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_longlong.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_longlong (int ncid, int varid, const char *name, long long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_LONGLONG;

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
      ierr = nc_get_att_longlong(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_longlong(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_longlong(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_LONG_LONG, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_schar
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_schar.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_schar (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const signed char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_SCHAR;

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
      ierr = nc_put_att_schar(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_schar(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_schar(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_float
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_float.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_float (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const float *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_FLOAT;

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
      ierr = nc_put_att_float(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_float(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_float(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_nvars
 * @brief   The PIO-C interface for the NetCDF function nc_inq_nvars.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_nvars (int ncid, int *nvarsp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_NVARS;

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
      ierr = nc_inq_nvars(file->fh, nvarsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_nvars(file->fh, nvarsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_nvars(file->fh, nvarsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(nvarsp,1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_rename_dim
 * @brief   The PIO-C interface for the NetCDF function nc_rename_dim.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_rename_dim (int ncid, int dimid, const char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_RENAME_DIM;

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
      ierr = nc_rename_dim(file->fh, dimid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_rename_dim(file->fh, dimid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_rename_dim(file->fh, dimid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_varndims
 * @brief   The PIO-C interface for the NetCDF function nc_inq_varndims.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_varndims (int ncid, int varid, int *ndimsp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VARNDIMS;
  if(file->varlist[varid].ndims > 0){
    (*ndimsp) = file->varlist[varid].ndims;
    return PIO_NOERR;
  }

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
      ierr = nc_inq_varndims(file->fh, varid, ndimsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_varndims(file->fh, varid, ndimsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_varndims(file->fh, varid, ndimsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(ndimsp,1, MPI_INT, ios->ioroot, ios->my_comm);
    file->varlist[varid].ndims = (*ndimsp);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_long
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_long.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_long (int ncid, int varid, const char *name, long *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_LONG;

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
      ierr = nc_get_att_long(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_long(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_long(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_LONG, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_dim
 * @brief   The PIO-C interface for the NetCDF function nc_inq_dim.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_inq_dim (int ncid, int dimid, char *name, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_DIM;

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
      ierr = nc_inq_dim(file->fh, dimid, name, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dim(file->fh, dimid, name, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dim(file->fh, dimid, name, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(name != NULL){ 
      int slen;
      if(ios->iomaster)
        slen = (int) strlen(name);
      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);
      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);
    }
      if(lenp != NULL) mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_dimid
 * @brief   The PIO-C interface for the NetCDF function nc_inq_dimid.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_inq_dimid (int ncid, const char *name, int *idp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_DIMID;

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
      ierr = nc_inq_dimid(file->fh, name, idp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimid(file->fh, name, idp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimid(file->fh, name, idp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_unlimdim
 * @brief   The PIO-C interface for the NetCDF function nc_inq_unlimdim.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_inq_unlimdim (int ncid, int *unlimdimidp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_UNLIMDIM;

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
      ierr = nc_inq_unlimdim(file->fh, unlimdimidp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_unlimdim(file->fh, unlimdimidp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_unlimdim(file->fh, unlimdimidp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(unlimdimidp,1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_vardimid
 * @brief   The PIO-C interface for the NetCDF function nc_inq_vardimid.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_vardimid (int ncid, int varid, int *dimidsp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VARDIMID;

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
      ierr = nc_inq_vardimid(file->fh, varid, dimidsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_vardimid(file->fh, varid, dimidsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_vardimid(file->fh, varid, dimidsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(ierr==PIO_NOERR){
      int ndims;
      PIOc_inq_varndims(file->fh, varid, &ndims);
      mpierr = MPI_Bcast(dimidsp , ndims, MPI_INT, ios->ioroot, ios->my_comm);
     }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_attlen
 * @brief   The PIO-C interface for the NetCDF function nc_inq_attlen.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_inq_attlen (int ncid, int varid, const char *name, PIO_Offset *lenp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_ATTLEN;

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
      ierr = nc_inq_attlen(file->fh, varid, name, (size_t *)lenp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attlen(file->fh, varid, name, (size_t *)lenp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attlen(file->fh, varid, name, lenp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(lenp , 1, MPI_OFFSET, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_dimname
 * @brief   The PIO-C interface for the NetCDF function nc_inq_dimname.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_inq_dimname (int ncid, int dimid, char *name) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_DIMNAME;

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
      ierr = nc_inq_dimname(file->fh, dimid, name);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_dimname(file->fh, dimid, name);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_dimname(file->fh, dimid, name);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    if(name != NULL){
      int slen;
      if(ios->iomaster)
        slen = (int) strlen(name);
      mpierr = MPI_Bcast(&slen, 1, MPI_INT, ios->ioroot, ios->my_comm);
      mpierr = MPI_Bcast(name, slen, MPI_CHAR, ios->ioroot, ios->my_comm);
    }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_ushort
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_ushort.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_ushort (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned short *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_USHORT;

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
      ierr = nc_put_att_ushort(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_ushort(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_ushort(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_float
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_float.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_float (int ncid, int varid, const char *name, float *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_FLOAT;

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
      ierr = nc_get_att_float(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_float(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_float(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_FLOAT, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_longlong
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_longlong.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_longlong (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const long long *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_LONGLONG;

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
      ierr = nc_put_att_longlong(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_longlong(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_longlong(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_uint
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_uint.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_uint (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned int *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_UINT;

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
      ierr = nc_put_att_uint(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_uint(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_uint(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_schar
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_schar.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_schar (int ncid, int varid, const char *name, signed char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_SCHAR;

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
      ierr = nc_get_att_schar(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_schar(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_schar(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_CHAR, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_attid
 * @brief   The PIO-C interface for the NetCDF function nc_inq_attid.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_inq_attid (int ncid, int varid, const char *name, int *idp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_ATTID;

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
      ierr = nc_inq_attid(file->fh, varid, name, idp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_attid(file->fh, varid, name, idp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_attid(file->fh, varid, name, idp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_def_dim
 * @brief   The PIO-C interface for the NetCDF function nc_def_dim.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_def_dim (int ncid, const char *name, PIO_Offset len, int *idp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_DEF_DIM;

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
      ierr = nc_def_dim(file->fh, name, (size_t)len, idp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_def_dim(file->fh, name, (size_t)len, idp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_def_dim(file->fh, name, len, idp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(idp , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_ndims
 * @brief   The PIO-C interface for the NetCDF function nc_inq_ndims.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__dimensions.html
 */
int PIOc_inq_ndims (int ncid, int *ndimsp) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_NDIMS;

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
      ierr = nc_inq_ndims(file->fh, ndimsp);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_ndims(file->fh, ndimsp);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_ndims(file->fh, ndimsp);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(ndimsp , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_vartype
 * @brief   The PIO-C interface for the NetCDF function nc_inq_vartype.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__variables.html
 */
int PIOc_inq_vartype (int ncid, int varid, nc_type *xtypep) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_VARTYPE;

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
      ierr = nc_inq_vartype(file->fh, varid, xtypep);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_vartype(file->fh, varid, xtypep);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_vartype(file->fh, varid, xtypep);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_int
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_int.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_int (int ncid, int varid, const char *name, int *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_INT;

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
      ierr = nc_get_att_int(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_int(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_int(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_INT, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_double
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_double.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_double (int ncid, int varid, const char *name, double *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_DOUBLE;

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
      ierr = nc_get_att_double(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_double(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_double(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_DOUBLE, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_ubyte
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_ubyte.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_ubyte (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_UBYTE;

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
      ierr = nc_put_att_ubyte(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_ubyte(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_ubyte(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_inq_atttype
 * @brief   The PIO-C interface for the NetCDF function nc_inq_atttype.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_inq_atttype (int ncid, int varid, const char *name, nc_type *xtypep) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_INQ_ATTTYPE;

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
      ierr = nc_inq_atttype(file->fh, varid, name, xtypep);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_inq_atttype(file->fh, varid, name, xtypep);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_inq_atttype(file->fh, varid, name, xtypep);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
    mpierr = MPI_Bcast(xtypep , 1, MPI_INT, ios->ioroot, ios->my_comm);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_put_att_uchar
 * @brief   The PIO-C interface for the NetCDF function nc_put_att_uchar.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_put_att_uchar (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned char *op) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_PUT_ATT_UCHAR;

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
      ierr = nc_put_att_uchar(file->fh, varid, name, xtype, (size_t)len, op);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_put_att_uchar(file->fh, varid, name, xtype, (size_t)len, op);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_put_att_uchar(file->fh, varid, name, xtype, len, op);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(__FILE__) + 20)* sizeof(char));
    sprintf(errstr,"in file %s",__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
  if(errstr != NULL) free(errstr);
  return ierr;
}

/** 
 * @name    PIOc_get_att_uchar
 * @brief   The PIO-C interface for the NetCDF function nc_get_att_uchar.
 * @details This routine is called collectively by all tasks in the communicator 
 *           ios.union_comm. For more information on the underlying NetCDF commmand
 *           please read about this function in the NetCDF documentation at: 
 *           http://www.unidata.ucar.edu/software/netcdf/docs/group__attributes.html
 */
int PIOc_get_att_uchar (int ncid, int varid, const char *name, unsigned char *ip) 
{
  int ierr;
  int msg;
  int mpierr;
  iosystem_desc_t *ios;
  file_desc_t *file;
  char *errstr;

  errstr = NULL;
  ierr = PIO_NOERR;

  file = pio_get_file_from_id(ncid);
  if(file == NULL)
    return PIO_EBADID;
  ios = file->iosystem;
  msg = PIO_MSG_GET_ATT_UCHAR;

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
      ierr = nc_get_att_uchar(file->fh, varid, name, ip);;
      break;
    case PIO_IOTYPE_NETCDF4C:
#endif
    case PIO_IOTYPE_NETCDF:
      if(ios->io_rank==0){
	ierr = nc_get_att_uchar(file->fh, varid, name, ip);;
      }
      break;
#endif
#ifdef _PNETCDF
    case PIO_IOTYPE_PNETCDF:
      ierr = ncmpi_get_att_uchar(file->fh, varid, name, ip);;
      break;
#endif
    default:
      ierr = iotype_error(file->iotype,__FILE__,__LINE__);
    }
  }

   if(ierr != PIO_NOERR){
    errstr = (char *) malloc((strlen(name)+strlen(__FILE__) + 40)* sizeof(char));
    sprintf(errstr,"name %s in file %s",name,__FILE__);
  }
  ierr = check_netcdf(file, ierr, errstr,__LINE__);
      if(ierr == PIO_NOERR){
        PIO_Offset attlen;
        PIOc_inq_attlen(file->fh, varid, name, &attlen);
        mpierr = MPI_Bcast(ip , (int) attlen, MPI_UNSIGNED_CHAR, ios->ioroot, ios->my_comm);
      }
  if(errstr != NULL) free(errstr);
  return ierr;
}

