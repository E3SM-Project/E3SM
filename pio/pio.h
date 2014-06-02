#ifndef _PIO_H_
#define _PIO_H_
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h> // memcpy
#include <mpi.h>

#ifdef _NETCDF
#include <netcdf.h>
#ifdef _NETCDF4
#include <netcdf_par.h>
#endif
#endif
#ifdef _PNETCDF
#include <pnetcdf.h>
#endif

#ifdef _MPISERIAL
typedef int MPI_Info;
typedef long long MPI_Offset;
#define MPI_OFFSET  ((MPI_Datatype)0x4c000844)
#define MPI_LONG_LONG ((MPI_Datatype)0x4c000809)
#define  MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_CHARACTER ((MPI_Datatype)1275068698)
#endif
#ifndef MPI_OFFSET
#define MPI_OFFSET  ((MPI_Datatype)0x4c000844)
#endif


#define PIO_OFFSET MPI_OFFSET
#define PIO_Offset MPI_Offset
#define PIO_MAX_VARS 2048

int PIOc_freedecomp(int iosysid, int ioid);


typedef struct var_desc_t
{
  int record;                      // the unlimited dimension index
  int ndims;                      // number of variable dimensions
  void *buffer;                  // data holder for buffered output
} var_desc_t;


typedef struct io_region
{
  int loffset;
  PIO_Offset *start;
  PIO_Offset *count;
  struct io_region *next;
} io_region;

typedef struct io_desc_t
{
  int ioid;
  int async_id;
  int nrecvs;
  int maxiobuflen;
  int ndof;
  int ndims;
  int num_aiotasks;
  int rearranger;
  int maxregions;

  MPI_Datatype basetype;
  PIO_Offset llen;


  int *rfrom;
  int *rcount;
  int *scount;
  PIO_Offset *sindex;
  PIO_Offset *rindex;

  MPI_Datatype *rtype;
  MPI_Datatype *stype;
  
  io_region *firstregion;

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
  int default_rearranger;

  bool async_interface;
  bool ioproc;
  
  MPI_Info info;
  struct iosystem_desc_t *next;
} iosystem_desc_t;

typedef struct file_desc_t
{
  iosystem_desc_t *iosystem;
  PIO_Offset buffsize;
  int fh;
  int iotype;
  struct var_desc_t varlist[PIO_MAX_VARS];
  MPI_Request request[PIO_MAX_VARS];   // request associated with buffered data for pnetcdf
  int nreq;   // next empty request slot to fill.
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

enum PIO_REARRANGERS{
  PIO_REARR_NONE = 0,
  PIO_REARR_BOX = 1,
  PIO_REARR_SUBSET = 2
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
#ifdef _PNETCDF
#define PIO_EINDEP  NC_EINDEP
#else
#define PIO_EINDEP  (-203) 
#endif
#if defined(__cplusplus)
extern "C" {
#endif

int PIOc_inq_att (int ncid, int varid, const char *name, nc_type *xtypep, PIO_Offset *lenp); 
int PIOc_inq_format (int ncid, int *formatp); 
int PIOc_inq_varid (int ncid, const char *name, int *varidp); 
int PIOc_inq_varnatts (int ncid, int varid, int *nattsp); 
int PIOc_def_var (int ncid, const char *name, nc_type xtype,  int ndims, const int *dimidsp, int *varidp); 
int PIOc_inq_var (int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp, int *dimidsp, int *nattsp); 
int PIOc_inq_varname (int ncid, int varid, char *name); 
int PIOc_put_att_double (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const double *op); 
int PIOc_put_att_int (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const int *op); 
int PIOc_rename_att (int ncid, int varid, const char *name, const char *newname); 
int PIOc_del_att (int ncid, int varid, const char *name); 
int PIOc_inq_natts (int ncid, int *ngattsp); 
int PIOc_inq (int ncid, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp);  
int PIOc_get_att_text (int ncid, int varid, const char *name, char *ip); 
int PIOc_get_att_short (int ncid, int varid, const char *name, short *ip); 
int PIOc_put_att_long (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const long *op); 
int PIOc_redef (int ncid); 
int PIOc_set_fill (int ncid, int fillmode, int *old_modep); 
int PIOc_enddef (int ncid); 
int PIOc_rename_var (int ncid, int varid, const char *name); 
int PIOc_put_att_short (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const short *op); 
int PIOc_put_att_text (int ncid, int varid, const char *name, PIO_Offset len, const char *op); 
int PIOc_inq_attname (int ncid, int varid, int attnum, char *name); 
int PIOc_get_att_ulonglong (int ncid, int varid, const char *name, unsigned long long *ip); 
int PIOc_get_att_ushort (int ncid, int varid, const char *name, unsigned short *ip); 
int PIOc_put_att_ulonglong (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned long long *op); 
int PIOc_inq_dimlen (int ncid, int dimid, PIO_Offset *lenp); 
int PIOc_get_att_uint (int ncid, int varid, const char *name, unsigned int *ip); 
int PIOc_get_att_longlong (int ncid, int varid, const char *name, long long *ip); 
int PIOc_put_att_schar (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const signed char *op); 
int PIOc_put_att_float (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const float *op); 
int PIOc_inq_nvars (int ncid, int *nvarsp); 
int PIOc_rename_dim (int ncid, int dimid, const char *name); 
int PIOc_inq_varndims (int ncid, int varid, int *ndimsp); 
int PIOc_get_att_long (int ncid, int varid, const char *name, long *ip); 
int PIOc_inq_dim (int ncid, int dimid, char *name, PIO_Offset *lenp); 
int PIOc_inq_dimid (int ncid, const char *name, int *idp); 
int PIOc_inq_unlimdim (int ncid, int *unlimdimidp); 
int PIOc_inq_vardimid (int ncid, int varid, int *dimidsp); 
int PIOc_inq_attlen (int ncid, int varid, const char *name, PIO_Offset *lenp); 
int PIOc_inq_dimname (int ncid, int dimid, char *name); 
int PIOc_put_att_ushort (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned short *op); 
int PIOc_get_att_float (int ncid, int varid, const char *name, float *ip); 
int PIOc_sync (int ncid); 
int PIOc_put_att_longlong (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const long long *op); 
int PIOc_put_att_uint (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned int *op); 
int PIOc_get_att_schar (int ncid, int varid, const char *name, signed char *ip); 
int PIOc_inq_attid (int ncid, int varid, const char *name, int *idp); 
int PIOc_def_dim (int ncid, const char *name, PIO_Offset len, int *idp); 
int PIOc_inq_ndims (int ncid, int *ndimsp); 
int PIOc_inq_vartype (int ncid, int varid, nc_type *xtypep); 
int PIOc_get_att_int (int ncid, int varid, const char *name, int *ip); 
int PIOc_get_att_double (int ncid, int varid, const char *name, double *ip); 
int PIOc_inq_atttype (int ncid, int varid, const char *name, nc_type *xtypep); 
int PIOc_put_att_uchar (int ncid, int varid, const char *name, nc_type xtype, PIO_Offset len, const unsigned char *op); 
int PIOc_get_att_uchar (int ncid, int varid, const char *name, unsigned char *ip);              
  int PIOc_InitDecomp(const int iosysid, const int basetype,const int ndims, const int dims[], 
		      const int maplen, const PIO_Offset *compmap, int *ioidp, const int *rearr, 
		      const PIO_Offset *iostart,const PIO_Offset *iocount);
  int PIOc_Init_Intracomm(const MPI_Comm comp_comm, 
			  const int num_iotasks, const int stride, 
			  const int base, const int rearr, int *iosysidp);
  int PIOc_closefile(int ncid);
int PIOc_createfile(const int iosysid, int *ncidp,  int *iotype,
		    const char *fname, const int mode);
int PIOc_openfile(const int iosysid, int *ncidp, int *iotype,
		  char *fname, const int mode, _Bool checkmpi);
  int PIOc_write_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen, void *array, void *fillvalue);








#if defined(__cplusplus)
}
#endif

#endif  // _PIO_H_

