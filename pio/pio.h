#ifndef _PIO_H_
#define _PIO_H_
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h> // memcpy

#ifdef _NETCDF
#include <netcdf.h>
#ifdef _NETCDF4
#include <netcdf_par.h>
#endif
#endif
#ifdef _PNETCDF
#include <pnetcdf.h>
#endif

#define PIO_Offset MPI_Offset
#define PIO_MAX_VARS 8192

int PIOc_freedecomp(int iosysid, int ioid);


typedef struct var_desc_t
{
  int record;
  void *buffer;
  MPI_Request request;
} var_desc_t;


typedef struct io_desc_t
{
  int ioid;
  int *dest_ioproc;
  int *rfrom;
  int *scount;
  int async_id;
  int nrecvs;
  int compsize;
  int maxiobuflen;
  int ndof;
  int ndims;
  MPI_Datatype basetype;
  MPI_Datatype *rtype;
  MPI_Datatype *stype;
  PIO_Offset *start;
  PIO_Offset *count;
  PIO_Offset llen;
  PIO_Offset glen;
  PIO_Offset *dest_ioindex;
  struct io_desc_t *next;
} io_desc_t;


typedef struct iosystem_desc_t 
{
  int iosysid;
  MPI_Comm union_comm;
  MPI_Comm io_comm;
  MPI_Comm comp_comm;
  MPI_Comm intercomm;
  MPI_Comm my_comm;
  
  int num_iotasks;
  int num_aiotasks;
  int num_comptasks;

  int union_rank;
  int comp_rank;
  int io_rank;

  bool iomaster;
  bool compmaster;

  int ioroot;
  int comproot;
  int *ioranks;

  int error_handler;

  bool async_interface;
  bool ioproc;
  
  MPI_Info info;
  struct iosystem_desc_t *next;
} iosystem_desc_t;

typedef struct file_desc_t
{
  iosystem_desc_t *iosystem;
  long int buffsize;
  int request_cnt;
  int fh;
  int iotype;
  struct var_desc_t varlist[PIO_MAX_VARS];
  struct file_desc_t *next;
} file_desc_t;

enum PIO_IOTYPE{
  PIO_IOTYPE_PBINARY=1,
  PIO_IOYTPE_DIRECT_PBINARY=2,
  PIO_IOTYPE_PNETCDF=5,
  PIO_IOTYPE_NETCDF=6,
  PIO_IOTYPE_NETCDF4C=7,
  PIO_IOTYPE_NETCDF4P=8
};

enum PIO_MSG{
  PIO_MSG_OPEN_FILE,
  PIO_MSG_CREATE_FILE,
  PIO_MSG_INQ_ATT,
  PIO_MSG_INQ_FORMAT,
  PIO_MSG_INQ_VARID,
  PIO_MSG_INQ_VARNATTS,
  PIO_MSG_DEF_VAR,
  PIO_MSG_INQ_VAR,
  PIO_MSG_INQ_VARNAME,
  PIO_MSG_PUT_ATT_DOUBLE,
  PIO_MSG_PUT_ATT_INT,
  PIO_MSG_RENAME_ATT,
  PIO_MSG_DEL_ATT,
  PIO_MSG_INQ_NATTS,
  PIO_MSG_INQ,
  PIO_MSG_GET_ATT_TEXT,
  PIO_MSG_GET_ATT_SHORT,
  PIO_MSG_PUT_ATT_LONG,
  PIO_MSG_REDEF,
  PIO_MSG_SET_FILL,
  PIO_MSG_ENDDEF,
  PIO_MSG_RENAME_VAR,
  PIO_MSG_PUT_ATT_SHORT,
  PIO_MSG_PUT_ATT_TEXT,
  PIO_MSG_INQ_ATTNAME,
  PIO_MSG_GET_ATT_ULONGLONG,
  PIO_MSG_GET_ATT_USHORT,
  PIO_MSG_PUT_ATT_ULONGLONG,
  PIO_MSG_INQ_DIMLEN,
  PIO_MSG_GET_ATT_UINT,
  PIO_MSG_GET_ATT_LONGLONG,
  PIO_MSG_PUT_ATT_SCHAR,
  PIO_MSG_PUT_ATT_FLOAT,
  PIO_MSG_INQ_NVARS,
  PIO_MSG_RENAME_DIM,
  PIO_MSG_INQ_VARNDIMS,
  PIO_MSG_GET_ATT_LONG,
  PIO_MSG_INQ_DIM,
  PIO_MSG_INQ_DIMID,
  PIO_MSG_INQ_UNLIMDIM,
  PIO_MSG_INQ_VARDIMID,
  PIO_MSG_INQ_ATTLEN,
  PIO_MSG_INQ_DIMNAME,
  PIO_MSG_PUT_ATT_USHORT,
  PIO_MSG_GET_ATT_FLOAT,
  PIO_MSG_SYNC,
  PIO_MSG_PUT_ATT_LONGLONG,
  PIO_MSG_PUT_ATT_UINT,
  PIO_MSG_GET_ATT_SCHAR,
  PIO_MSG_INQ_ATTID,
  PIO_MSG_DEF_DIM,
  PIO_MSG_INQ_NDIMS,
  PIO_MSG_INQ_VARTYPE,
  PIO_MSG_GET_ATT_INT,
  PIO_MSG_GET_ATT_DOUBLE,
  PIO_MSG_INQ_ATTTYPE,
  PIO_MSG_PUT_ATT_UCHAR,
  PIO_MSG_GET_ATT_UCHAR
};

enum PIO_ERROR_HANDLERS{
  PIO_INTERNAL_ERROR=(-51),
  PIO_BCAST_ERROR=(-52),
  PIO_RETURN_ERROR=(-53)
};

#if defined( _PNETCDF) || defined(_NETCDF)
#define PIO_GLOBAL NC_GLOBAL
#define PIO_UNLIMITED NC_UNLIMITED
#define PIO_DOUBLE NC_DOUBLE
#define PIO_REAL   NC_FLOAT
#define PIO_INT    NC_INT
#define PIO_CHAR   NC_CHAR
#define PIO_NOERR  NC_NOERR
#define PIO_WRITE  NC_WRITE
#define PIO_NOWRITE  NC_NOWRITE
#define PIO_CLOBBER NC_CLOBBER	
#define PIO_NOCLOBBER NC_NOCLOBBER	
#define PIO_NOFILL NC_NOFILL
#define PIO_MAX_NAME NC_MAX_NAME
#define PIO_MAX_VAR_DIMS NC_MAX_VAR_DIMS
#define PIO_64BIT_OFFSET NC_64BIT_OFFSET
// NC_64BIT_DATA  This is a problem - need to define directly instead of using include file
#define PIO_64BIT_DATA 0x0010  
#define PIO_EBADID NC_EBADID
#define PIO_ENFILE NC_ENFILE
#define PIO_EEXIST NC_EEXIST
#define PIO_EINVAL NC_EINVAL
#define PIO_EPERM NC_EPERM
#define PIO_ENOTINDEFINE NC_ENOTINDEFINE
#define PIO_EINDEFINE NC_EINDEFINE
#define PIO_EINVALCOORDS NC_EINVALCOORDS
#define PIO_EMAXDIMS NC_EMAXDIMS
#define PIO_ENAMEINUSE NC_ENAMEINUSE
#define PIO_ENOTATT NC_ENOTATT
#define PIO_EMAXATTS NC_EMAXATTS
#define PIO_EBADTYPE NC_EBADTYPE
#define PIO_EBADDIM NC_EBADDIM
#define PIO_EUNLIMPOS NC_EUNLIMPOS
#define PIO_EMAXVARS NC_EMAXVARS
#define PIO_ENOTVAR NC_ENOTVAR
#define PIO_EGLOBAL NC_EGLOBAL
#define PIO_ENOTNC NC_ENOTNC
#define PIO_ESTS NC_ESTS
#define PIO_EMAXNAME NC_EMAXNAME
#define PIO_EUNLIMIT NC_EUNLIMIT
#define PIO_ENORECVARS NC_ENORECVARS
#define PIO_ECHAR NC_ECHAR
#define PIO_EEDGE NC_EEDGE
#define PIO_ESTRIDE NC_ESTRIDE
#define PIO_EBADNAME NC_EBADNAME
#define PIO_ERANGE NC_ERANGE
#define PIO_ENOMEM NC_ENOMEM
#define PIO_EVARSIZE NC_EVARSIZE
#define PIO_EDIMSIZE NC_EDIMSIZE
#define PIO_ETRUNC NC_ETRUNC
#define PIO_EAXISTYPE NC_EAXISTYPE
#define PIO_EDAP NC_EDAP
#define PIO_ECURL NC_ECURL
#define PIO_EIO NC_EIO
#define PIO_ENODATA NC_ENODATA
#define PIO_EDAPSVC NC_EDAPSVC
#define PIO_EDAS NC_EDAS
#define PIO_EDDS NC_EDDS
#define PIO_EDATADDS NC_EDATADDS
#define PIO_EDAPURL NC_EDAPURL
#define PIO_EDAPCONSTRAINT NC_EDAPCONSTRAINT
#define PIO_ETRANSLATION NC_ETRANSLATION
#define PIO_EHDFERR NC_EHDFERR
#define PIO_ECANTREAD NC_ECANTREAD
#define PIO_ECANTWRITE NC_ECANTWRITE
#define PIO_ECANTCREATE NC_ECANTCREATE
#define PIO_EFILEMETA NC_EFILEMETA
#define PIO_EDIMMETA NC_EDIMMETA
#define PIO_EATTMETA NC_EATTMETA
#define PIO_EVARMETA NC_EVARMETA
#define PIO_ENOCOMPOUND NC_ENOCOMPOUND
#define PIO_EATTEXISTS NC_EATTEXISTS
#define PIO_ENOTNC4 NC_ENOTNC4
#define PIO_ESTRICTNC3 NC_ESTRICTNC3
#define PIO_ENOTNC3 NC_ENOTNC3
#define PIO_ENOPAR NC_ENOPAR
#define PIO_EPARINIT NC_EPARINIT
#define PIO_EBADGRPID NC_EBADGRPID
#define PIO_EBADTYPID NC_EBADTYPID
#define PIO_ETYPDEFINED NC_ETYPDEFINED
#define PIO_EBADFIELD NC_EBADFIELD
#define PIO_EBADCLASS NC_EBADCLASS
#define PIO_EMAPTYPE NC_EMAPTYPE
#define PIO_ELATEFILL NC_ELATEFILL
#define PIO_ELATEDEF NC_ELATEDEF
#define PIO_EDIMSCALE NC_EDIMSCALE
#define PIO_ENOGRP NC_ENOGRP
#define PIO_ESTORAGE NC_ESTORAGE
#define PIO_EBADCHUNK NC_EBADCHUNK
#define PIO_ENOTBUILT NC_ENOTBUILT
#define PIO_EDISKLESS NC_EDISKLESS



#endif

#endif  _PIO_H_
