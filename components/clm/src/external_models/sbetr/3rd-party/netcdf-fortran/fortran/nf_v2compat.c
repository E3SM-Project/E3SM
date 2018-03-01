/*
 *	Copyright 1996, University Corporation for Atmospheric Research
 *      See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */

/* $Id: fort-v2compat.c,v 1.33 2009/01/27 19:48:34 ed Exp $ */

/*
 *  Source for netCDF2 FORTRAN jacket library.
 */

/* Modified version of fort-v2compat.c used to provide required
 * C functions used in v2 compatability interface. This clone
 * was created to keep existing C code fort-v2compat.c pristine
 * and to make compiling easier. Note all cfortran.h stuff has  
 * been removed to make compiling easier and the functions
 * have been made external instead of static so that FORTRAN can
 * see them
 */

/* April, 2009
 * Modified by:  Richard Weed, Ph.D
 *               Center for Advanced Vehicular Systems
 *               Mississippi State University
 *               rweed@cavs.msstate.edu
 *
 *  C routines required for Fortran V2 compatability 
 */

/*
 * OVERVIEW
 *
 * This file contains jacket routines written in C for interfacing
 * Fortran netCDF-2 function calls to the actual C-binding netCDF
 * function call -- using either the netCDF-2 or netCDF-3 C API.
 * In general, these functions handle character-string parameter
 * conventions, convert between column-major-order arrays and
 * row-major-order arrays, and map between array indices beginning
 * at one and array indices beginning at zero.  They also adapt the
 * differing error handling mechanisms between version 2 and version 3.
 */

#include <config.h>

#ifndef NO_NETCDF_2

/* LINTLIBRARY */

#include	<ctype.h>
#include        <string.h>
#include	<stdlib.h>
#include	<stdio.h>
#include	"netcdf.h"
#include	"nfconfig.inc"

#ifndef USE_NETCDF4
#define NC_CLASSIC_MODEL 0
#else
/* There is a dependency error here;
NC_CLASSIC_MODEL will not be defined
if ../libsrc4/netcdf.h does not exist yet
(which it won't after a maintainer-clean).
So, define it here if not already defined.
*/
#ifndef NC_CLASSIC_MODEL
#define NC_CLASSIC_MODEL 0x0100
#endif
#endif


/*
 New function added by RW to support FORTRAN 2003 interfaces.
 Function to return C data type sizes to FORTRAN 2003 code for
 v2 imap conversion. Duplicates some code in f2c_vimap below
*/
extern size_t
v2data_size(nc_type datatype)
{
 size_t size;

 size = 0;
 switch (datatype)
  {

    case NC_CHAR:
      size = sizeof(char);
    break;
    case NC_BYTE:
#if NF_INT1_IS_C_SIGNED_CHAR
      size = sizeof(signed char);
#elif NF_INT1_IS_C_SHORT
      size = sizeof(short);
#elif NF_INT1_IS_C_INT
      size = sizeof(int);
#elif NF_INT1_IS_C_LONG
      size = sizeof(long);
#endif
    break;
    case NC_SHORT:
#if NF_INT2_IS_C_SHORT
      size = sizeof(short);
#elif NF_INT2_IS_C_INT
      size = sizeof(int);
#elif NF_INT2_IS_C_LONG
      size = sizeof(long);
#endif
    break;
    case NC_INT:
#if NF_INT_IS_C_INT
      size = sizeof(int);
#elif NF_INT_IS_C_LONG
      size = sizeof(long);
#endif
    break;
    case NC_FLOAT:
#if NF_REAL_IS_C_FLOAT
      size = sizeof(float);
#elif NF_REAL_IS_C_DOUBLE
      size = sizeof(double);
#endif
    break;
    case NC_DOUBLE:
#if NF_DOUBLEPRECISION_IS_C_FLOAT
      size = sizeof(float);
#elif NF_DOUBLEPRECISION_IS_C_DOUBLE
      size = sizeof(double);
#endif
    break;
    default:
      size = -1;
  }
 return size;
}

/**
 * Convert a Version 2 Fortran IMAP vector into a Version 3 C imap vector.
 */
extern ptrdiff_t*
f2c_v2imap(int ncid, int varid, const int* fimap, ptrdiff_t* cimap)
{
    int		rank;
    nc_type	datatype;

    if (nc_inq_vartype(ncid, varid, &datatype) ||
	nc_inq_varndims(ncid, varid, &rank) || rank <= 0)
    {
	return NULL;
    }

    /* else */
    if (fimap[0] == 0)
    {
	/*
	 * Special Fortran version 2 semantics: use external netCDF variable 
	 * structure.
	 */
	int		dimids[NC_MAX_VAR_DIMS];
	int		idim;
	size_t	total;

	if (nc_inq_vardimid(ncid, varid, dimids) != NC_NOERR)
	    return NULL;

	for (total = 1, idim = rank - 1; idim >= 0; --idim)
	{
	    size_t	length;

	    cimap[idim] = total;

	    if (nc_inq_dimlen(ncid, dimids[idim], &length) != NC_NOERR)
		return NULL;

	    total *= length;
	}
    }
    else
    {
	/*
	 * Regular Fortran version 2 semantics: convert byte counts to
	 * element counts.
	 */
	int	idim;
	size_t	size;

	switch (datatype)
	{

	    case NC_CHAR:
		size = sizeof(char);
		break;
	    case NC_BYTE:
#		if NF_INT1_IS_C_SIGNED_CHAR
		    size = sizeof(signed char);
#		elif NF_INT1_IS_C_SHORT
		    size = sizeof(short);
#		elif NF_INT1_IS_C_INT
		    size = sizeof(int);
#		elif NF_INT1_IS_C_LONG
		    size = sizeof(long);
#		endif
		break;
	    case NC_SHORT:
#		if NF_INT2_IS_C_SHORT
		    size = sizeof(short);
#		elif NF_INT2_IS_C_INT
		    size = sizeof(int);
#		elif NF_INT2_IS_C_LONG
		    size = sizeof(long);
#		endif
		break;
	    case NC_INT:
#		if NF_INT_IS_C_INT
		    size = sizeof(int);
#		elif NF_INT_IS_C_LONG
		    size = sizeof(long);
#		endif
		break;
	    case NC_FLOAT:
#		if NF_REAL_IS_C_FLOAT
		    size = sizeof(float);
#		elif NF_REAL_IS_C_DOUBLE
		    size = sizeof(double);
#		endif
		break;
	    case NC_DOUBLE:
#		if NF_DOUBLEPRECISION_IS_C_FLOAT
		    size = sizeof(float);
#		elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		    size = sizeof(double);
#		endif
		break;
	    default:
		return NULL;
	}

	for (idim = 0; idim < rank; ++idim)
	    cimap[idim] = fimap[rank - 1 - idim] / size;
    }

    return cimap;
}


/*
 * Compute the product of dimensional counts.
 */
static size_t
dimprod(const size_t* count, int rank)
{
    int		i;
    size_t	prod = 1;

    for (i = 0; i < rank; ++i)
	prod *= count[i];

    return prod;
}


/*
 * Set the C global variable ncopts.
 */
extern void
c_ncpopt(
    int val     /* NC_FATAL, NC_VERBOSE, or NC_FATAL|NC_VERBOSE */
)
{
    ncopts = val;
}

/*
 * Get the C global variable ncopts from FORTRAN.
 */
extern void
c_ncgopt(
    int	*val	/* NC_FATAL, NC_VERBOSE, or NC_FATAL|NC_VERBOSE */
)
{
    *val = ncopts;
}



/*
 * Create a new netCDF file, returning a netCDF ID.  New netCDF
 * file is placed in define mode.
 */
extern int
c_nccre(
    const char *pathname,	/* file name of new netCDF file */
    int clobmode,	/* either NCCLOB or NCNOCLOB */
    int *rcode		/* returned error code */
)
{
    int ncid = -1;

    if (pathname == NULL)
       *rcode = NC_EINVAL;
    else
    {
       *rcode = ((ncid = nccreate (pathname, clobmode)) == -1)
	  ? ncerr
	  : 0;
    }
    
    if (*rcode != 0)
    {
       nc_advise("NCCRE", *rcode, "");
       *rcode = ncerr;
    }

    return ncid;
}



/*
 * Open an existing netCDF file for access.
 */
extern int
c_ncopn(
    const char *pathname,	/* file name for netCDF to be opened */
    int rwmode,			/* either NCWRITE or NCNOWRIT */
    int *rcode			/* returned error code */
)
{
    int ncid = -1;

    /* Include NC_LOCK in check, in case NC_LOCK is ever implemented */
    if (rwmode < 0 ||
	rwmode > NC_WRITE + NC_SHARE + NC_CLASSIC_MODEL + NC_LOCK)
    {
        *rcode = NC_EINVAL;
        nc_advise("NCOPN", *rcode,
		"bad flag, did you forget to include netcdf.inc?");
    }
    else
    {
	if (pathname == NULL) {
	    *rcode = NC_EINVAL;
	}
	else
	{
	    *rcode = ((ncid = ncopen (pathname, rwmode)) == -1)
			? ncerr
			: 0;
	}

	if (*rcode != 0)
	{
	    nc_advise("NCOPN", *rcode, "");
	    *rcode = ncerr;
	}
    }

    return ncid;
}


/*
 * Add a new dimension to an open netCDF file in define mode.
 */
extern int
c_ncddef (
    int ncid,		/* netCDF ID */
    const char *dimname,/* dimension name */
    int dimlen,		/* size of dimension */
    int *rcode		/* returned error code */
)
{
    int dimid;

    if ((dimid = ncdimdef (ncid, dimname, (long)dimlen)) == -1)
	*rcode = ncerr;
    else
    {
	dimid++;
	*rcode = 0;
    }

    return dimid;
}


/*
 * Return the ID of a netCDF dimension, given the name of the dimension.
 */
extern int
c_ncdid (
    int ncid,		/* netCDF ID */
    const char *dimname,/* dimension name */
    int *rcode		/* returned error code */
)
{
    int dimid;

    if ((dimid = ncdimid (ncid, dimname)) == -1)
	*rcode = ncerr;
    else
    {
	dimid++;
	*rcode = 0;
    }

    return dimid;
}

/*
 * Add a new variable to an open netCDF file in define mode.
 */
extern int
c_ncvdef (
    int ncid,           /* netCDF ID */
    const char *varname,/* name of variable */
    nc_type datatype,   /* netCDF datatype of variable */
    int ndims,          /* number of dimensions of variable */
    int *dimids,        /* array of ndims dimensions IDs */
    int *rcode          /* returned error code */
)
{
    int varid, status;

    if ((status = nc_def_var(ncid, varname, datatype, ndims, dimids, &varid)))
    {
        nc_advise("NCVDEF", status, "");
        *rcode = ncerr;
        varid = -1;
    }
    else
    {
        varid++;
        *rcode = 0;
    }

    return varid;
}



/*
 * Return the ID of a netCDF variable given its name.
 */
extern int
c_ncvid (
    int ncid,		/* netCDF ID */
    const char *varname,/* variable name */
    int *rcode		/* returned error code */
)
{
    int varid;

    if ((varid = ncvarid (ncid, varname)) == -1)
	*rcode = ncerr;
    else
    {
	varid++;
	*rcode = 0;
    }

    return varid;
}


/*
 * Return number of bytes per netCDF data type.
 */
extern int
c_nctlen (
    nc_type datatype,	/* netCDF datatype */
    int* rcode		/* returned error code */
)
{
    int itype;

    *rcode = ((itype = (int) nctypelen (datatype)) == -1)
		?  ncerr
		: 0;

    return itype;
}

/*
 * Close an open netCDF file.
 */
extern void
c_ncclos (
    int ncid,		/* netCDF ID */
    int* rcode		/* returned error code */
)
{
    *rcode = ncclose(ncid) == -1
		? ncerr
		: 0;
}

/*
 * Put an open netCDF into define mode.
 */
extern void
c_ncredf (
    int ncid,		/* netCDF ID */
    int *rcode		/* returned error code */
)
{
    *rcode = ncredef(ncid) == -1
		? ncerr
		: 0;
}

/*
 * Take an open netCDF out of define mode.
 */
extern void
c_ncendf (
    int ncid,		/* netCDF ID */
    int *rcode		/* returned error code */
)
{
    *rcode = ncendef (ncid) == -1
		? ncerr
		: 0;
}

/*
 * Return information about an open netCDF file given its netCDF ID.
 */
extern void
c_ncinq (
    int ncid,		/* netCDF ID */
    int* indims,	/* returned number of dimensions */
    int* invars,	/* returned number of variables */
    int* inatts,	/* returned number of attributes */
    int* irecdim,	/* returned ID of the unlimited dimension */
    int* rcode		/* returned error code */
)
{
    *rcode = ncinquire(ncid, indims, invars, inatts, irecdim) == -1
		? ncerr
		: 0;
}

/*
 * Make sure that the disk copy of a netCDF file open for writing
 * is current.
 */
extern void
c_ncsnc(
    int ncid,		/* netCDF ID */
    int* rcode		/* returned error code */
)
{
    *rcode = ncsync (ncid) == -1
		? ncerr
		: 0;
}

/*
 * Restore the netCDF to a known consistent state in case anything
 * goes wrong during the definition of new dimensions, variables
 * or attributes.
 */
extern void
c_ncabor (
    int ncid,		/* netCDF ID */
    int* rcode		/* returned error code */
)
{
    *rcode = ncabort(ncid) == -1
		? ncerr
		: 0;
}


/*
 * Return the name and size of a dimension, given its ID.
 */
extern void
c_ncdinq (
    int ncid,			/* netCDF ID */
    int dimid,			/* dimension ID */
    char* dimname,		/* returned dimension name */
    int* size,			/* returned dimension size */
    int* rcode			/* returned error code */
)
{
    long siz;

    if (ncdiminq (ncid, dimid, dimname, &siz) == -1)
	*rcode = ncerr;
    else
    {
	*size = siz;
	*rcode = 0;
    }
}

/*
 * Rename an existing dimension in a netCDF open for writing.
 */
extern void
c_ncdren (
    int ncid,			/* netCDF ID */
    int dimid,			/* dimension ID */
    const char* dimname,	/* new name of dimension */
    int* rcode			/* returned error code */
)
{
    *rcode = ncdimrename(ncid, dimid, dimname) == -1
		? ncerr
		: 0;
}


/*
 * Return information about a netCDF variable, given its ID.
 */
extern void
c_ncvinq (
    int ncid,		/* netCDF ID */
    int varid,		/* variable ID */
    char* varname,	/* returned variable name */
    nc_type* datatype,	/* returned variable type */
    int* indims,	/* returned number of dimensions */
    int* dimarray,	/* returned array of ndims dimension IDs */
    int* inatts,	/* returned number of attributes */
    int* rcode		/* returned error code */
)
{
    *rcode = ncvarinq(ncid, varid, varname, datatype, indims,
		      dimarray, inatts) == -1
		? ncerr
		: 0;
}

/*
 * Put a single numeric data value into a variable of an open netCDF.
 */
extern void
c_ncvpt1 (
    int			ncid,	/* netCDF ID */
    int	 		varid,	/* variable ID */
    const size_t*	indices,/* multidim index of data to be written */
    const void*		value,	/* pointer to data value to be written */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = NC_ECHAR;
	    break;
	case NC_BYTE:
#	    if NF_INT1_IS_C_SIGNED_CHAR
		status = nc_put_var1_schar(ncid, varid, indices,
					   (const signed char*)value);
#	    elif NF_INT1_IS_C_SHORT
		status = nc_put_var1_short(ncid, varid, indices,
					   (const short*)value);
#	    elif NF_INT1_IS_C_INT
		status = nc_put_var1_int(ncid, varid, indices,
					   (const int*)value);
#	    elif NF_INT1_IS_C_LONG
		status = nc_put_var1_long(ncid, varid, indices,
					   (const long*)value);
#	    endif
	    break;
	case NC_SHORT:
#	    if NF_INT2_IS_C_SHORT
		status = nc_put_var1_short(ncid, varid, indices,
					   (const short*)value);
#	    elif NF_INT2_IS_C_INT
		status = nc_put_var1_int(ncid, varid, indices,
					   (const int*)value);
#	    elif NF_INT2_IS_C_LONG
		status = nc_put_var1_long(ncid, varid, indices,
					   (const long*)value);
#	    endif
	    break;
	case NC_INT:
#	    if NF_INT_IS_C_INT
		status = nc_put_var1_int(ncid, varid, indices,
					   (const int*)value);
#	    elif NF_INT_IS_C_LONG
		status = nc_put_var1_long(ncid, varid, indices,
					   (const long*)value);
#	    endif
	    break;
	case NC_FLOAT:
#	    if NF_REAL_IS_C_FLOAT
		status = nc_put_var1_float(ncid, varid, indices,
					   (const float*)value);
#	    elif NF_REAL_IS_C_DOUBLE
		status = nc_put_var1_double(ncid, varid, indices,
					   (const double*)value);
#	    endif
	    break;
	case NC_DOUBLE:
#	    if NF_DOUBLEPRECISION_IS_C_FLOAT
		status = nc_put_var1_float(ncid, varid, indices,
					   (const float*)value);
#	    elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		status = nc_put_var1_double(ncid, varid, indices,
					   (const double*)value);
#	    endif
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVPT1", status, "");
	*rcode = ncerr;
    }
}

/* 
 * Put a single character into an open netCDF file.
 */
extern void
c_ncvp1c(
    int			ncid,	/* netCDF ID */
    int	 		varid,	/* variable ID */
    const size_t*	indices,/* multidim index of data to be written */
    const char*		value,	/* pointer to data value to be written */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	status = datatype != NC_CHAR
		    ? NC_ECHAR
		    : nc_put_var1_text(ncid, varid, indices, value);
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVP1C", status, "");
	*rcode = ncerr;
    }
}

/*
 * Write a hypercube of numeric values into a netCDF variable of an open
 * netCDF file.
 */
extern void
c_ncvpt (
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    const void*		value,	/* block of data values to be written */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = NC_ECHAR;
	    break;
	case NC_BYTE:
#	    if NF_INT1_IS_C_SIGNED_CHAR
		status = nc_put_vara_schar(ncid, varid, start, count,
					   (const signed char*)value);
#	    elif NF_INT1_IS_C_SHORT
		status = nc_put_vara_short(ncid, varid, start, count,
					   (const short*)value);
#	    elif NF_INT1_IS_C_INT
		status = nc_put_vara_int(ncid, varid, start, count,
					   (const int*)value);
#	    elif NF_INT1_IS_C_LONG
		status = nc_put_vara_long(ncid, varid, start, count,
					   (const long*)value);
#	    endif
	    break;
	case NC_SHORT:
#	    if NF_INT2_IS_C_SHORT
		status = nc_put_vara_short(ncid, varid, start, count,
					   (const short*)value);
#	    elif NF_INT2_IS_C_INT
		status = nc_put_vara_int(ncid, varid, start, count,
					   (const int*)value);
#	    elif NF_INT2_IS_C_LONG
		status = nc_put_vara_long(ncid, varid, start, count,
					   (const long*)value);
#	    endif
	    break;
	case NC_INT:
#	    if NF_INT_IS_C_INT
		status = nc_put_vara_int(ncid, varid, start, count,
					   (const int*)value);
#	    elif NF_INT_IS_C_LONG
		status = nc_put_vara_long(ncid, varid, start, count,
					   (const long*)value);
#	    endif
	    break;
	case NC_FLOAT:
#	    if NF_REAL_IS_C_FLOAT
		status = nc_put_vara_float(ncid, varid, start, count,
					   (const float*)value);
#	    elif NF_REAL_IS_C_DOUBLE
		status = nc_put_vara_double(ncid, varid, start, count,
					   (const double*)value);
#	    endif
	    break;
	case NC_DOUBLE:
#	    if NF_DOUBLEPRECISION_IS_C_FLOAT
		status = nc_put_vara_float(ncid, varid, start, count,
					   (const float*)value);
#	    elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		status = nc_put_vara_double(ncid, varid, start, count,
					   (const double*)value);
#	    endif
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVPT", status, "");
	*rcode = ncerr;
    }
}


/*
 * Write a hypercube of character values into an open netCDF file.
 */
extern void
c_ncvptc(
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    const char*		value,	/* block of data values to be written */
    int			lenstr,	/* declared length of the data argument */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	if (datatype != NC_CHAR)
	    status = NC_ECHAR;
	else
	{
	    int	rank;

	    status = nc_inq_varndims(ncid, varid, &rank);
	    if (status == 0)
	    {
		if (dimprod(count, rank) > (size_t)lenstr)
		    status = NC_ESTS;
		else
		    status = nc_put_vara_text(ncid, varid, start, count, value);
	    }
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVPTC", status, "");
	*rcode = ncerr;
    }
}


/*
 * Write a generalized hypercube of numeric values into a netCDF variable of 
 * an open netCDF file.
 */
extern void
c_ncvptg (
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    const ptrdiff_t*	strides,/* netCDF variable access strides */
    const ptrdiff_t*	imap,	/* memory values access mapping vector */
    const void*		value,	/* block of data values to be written */
    int*		rcode	/* returned error code */
)
{
    int		status;
    int		rank;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0 &&
	(status = nc_inq_varndims(ncid, varid, &rank)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = NC_ECHAR;
	    break;
	case NC_BYTE:
#	    if NF_INT1_IS_C_SIGNED_CHAR
		status = nc_put_varm_schar(ncid, varid, start, count,
					   strides, imap,
					   (const signed char*)value);
#	    elif NF_INT1_IS_C_SHORT
		status = nc_put_varm_short(ncid, varid, start, count,
					   strides, imap,
					   (const short*)value);
#	    elif NF_INT1_IS_C_INT
		status = nc_put_varm_int(ncid, varid, start, count,
					   strides, imap,
					   (const int*)value);
#	    elif NF_INT1_IS_C_LONG
		status = nc_put_varm_long(ncid, varid, start, count,
					   strides, imap,
					   (const long*)value);
#	    endif
	    break;
	case NC_SHORT:
#	    if NF_INT2_IS_C_SHORT
		status = nc_put_varm_short(ncid, varid, start, count,
					   strides, imap,
					   (const short*)value);
#	    elif NF_INT2_IS_C_INT
		status = nc_put_varm_int(ncid, varid, start, count,
					   strides, imap,
					   (const int*)value);
#	    elif NF_INT2_IS_C_LONG
		status = nc_put_varm_long(ncid, varid, start, count,
					   strides, imap,
					   (const long*)value);
#	    endif
	    break;
	case NC_INT:
#	    if NF_INT_IS_C_INT
		status = nc_put_varm_int(ncid, varid, start, count,
					   strides, imap,
					   (const int*)value);
#	    elif NF_INT_IS_C_LONG
		status = nc_put_varm_long(ncid, varid, start, count,
					   strides, imap,
					   (const long*)value);
#	    endif
	    break;
	case NC_FLOAT:
#	    if NF_REAL_IS_C_FLOAT
		status = nc_put_varm_float(ncid, varid, start, count,
					   strides, imap,
					   (const float*)value);
#	    elif NF_REAL_IS_C_DOUBLE
		status = nc_put_varm_double(ncid, varid, start, count,
					   strides, imap,
					   (const double*)value);
#	    endif
	    break;
	case NC_DOUBLE:
#	    if NF_DOUBLEPRECISION_IS_C_FLOAT
		status = nc_put_varm_float(ncid, varid, start, count,
					   strides, imap,
					   (const float*)value);
#	    elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		status = nc_put_varm_double(ncid, varid, start, count,
					   strides, imap,
					   (const double*)value);
#	    endif
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVPTG", status, "");
	*rcode = ncerr;
    }
}


/*
 * Write a generalized hypercube of character values into a netCDF variable of 
 * an open netCDF file.
 */
extern void
c_ncvpgc(
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    const ptrdiff_t*	strides,/* netCDF variable access strides */
    const ptrdiff_t*	imap,	/* memory values access mapping vector */
    const char*		value,	/* block of data values to be written */
    int*		rcode	/* returned error code */
)
{
    int		status;
    int		rank;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0 &&
	(status = nc_inq_varndims(ncid, varid, &rank)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = nc_put_varm_text(ncid, varid, start, count,
				       strides, imap,
				       value);
	    break;
	default:
	    status = NC_ECHAR;
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVPGC", status, "");
	*rcode = ncerr;
    }
}


/*
 * Get a single numeric value from a variable of an open netCDF file.
 */
extern void
c_ncvgt1 (
    int			ncid,	/* netCDF ID */
    int	 		varid,	/* variable ID */
    const size_t*	indices,/* multidim index of data to be read */
    void*		value,	/* pointer to data value to be read */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = NC_ECHAR;
	    break;
	case NC_BYTE:
#	    if NF_INT1_IS_C_SIGNED_CHAR
		status = nc_get_var1_schar(ncid, varid, indices,
					   (signed char*)value);
#	    elif NF_INT1_IS_C_SHORT
		status = nc_get_var1_short(ncid, varid, indices,
					   (short*)value);
#	    elif NF_INT1_IS_C_INT
		status = nc_get_var1_int(ncid, varid, indices,
					   (int*)value);
#	    elif NF_INT1_IS_C_LONG
		status = nc_get_var1_long(ncid, varid, indices,
					   (long*)value);
#	    endif
	    break;
	case NC_SHORT:
#	    if NF_INT2_IS_C_SHORT
		status = nc_get_var1_short(ncid, varid, indices,
					   (short*)value);
#	    elif NF_INT2_IS_C_INT
		status = nc_get_var1_int(ncid, varid, indices,
					   (int*)value);
#	    elif NF_INT2_IS_C_LONG
		status = nc_get_var1_long(ncid, varid, indices,
					   (long*)value);
#	    endif
	    break;
	case NC_INT:
#	    if NF_INT_IS_C_INT
		status = nc_get_var1_int(ncid, varid, indices,
					   (int*)value);
#	    elif NF_INT_IS_C_LONG
		status = nc_get_var1_long(ncid, varid, indices,
					   (long*)value);
#	    endif
	    break;
	case NC_FLOAT:
#	    if NF_REAL_IS_C_FLOAT
		status = nc_get_var1_float(ncid, varid, indices,
					   (float*)value);
#	    elif NF_REAL_IS_C_DOUBLE
		status = nc_get_var1_double(ncid, varid, indices,
					   (double*)value);
#	    endif
	    break;
	case NC_DOUBLE:
#	    if NF_DOUBLEPRECISION_IS_C_FLOAT
		status = nc_get_var1_float(ncid, varid, indices,
					   (float*)value);
#	    elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		status = nc_get_var1_double(ncid, varid, indices,
					   (double*)value);
#	    endif
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVGT1", status, "");
	*rcode = ncerr;
    }
}


/*
 * Get a single character data value from a variable of an open
 * netCDF file.
 */
extern void
c_ncvg1c(
    int			ncid,	/* netCDF ID */
    int	 		varid,	/* variable ID */
    const size_t*	indices,/* multidim index of data to be read */
    char*		value,	/* pointer to data value to be read */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = nc_get_var1_text(ncid, varid, indices, value);
	    break;
	default:
	    status = NC_ECHAR;
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVG1C", status, "");
	*rcode = ncerr;
    }
}


/*
 * Read a hypercube of numeric values from a netCDF variable of an open
 * netCDF file.
 */
extern void
c_ncvgt(
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    void*		value,	/* block of data values to be read */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = NC_ECHAR;
	    break;
	case NC_BYTE:
#	    if NF_INT1_IS_C_SIGNED_CHAR
		status = nc_get_vara_schar(ncid, varid, start, count,
					   (signed char*)value);
#	    elif NF_INT1_IS_C_SHORT
		status = nc_get_vara_short(ncid, varid, start, count,
					   (short*)value);
#	    elif NF_INT1_IS_C_INT
		status = nc_get_vara_int(ncid, varid, start, count,
					   (int*)value);
#	    elif NF_INT1_IS_C_LONG
		status = nc_get_vara_long(ncid, varid, start, count,
					   (long*)value);
#	    endif
	    break;
	case NC_SHORT:
#	    if NF_INT2_IS_C_SHORT
		status = nc_get_vara_short(ncid, varid, start, count,
					   (short*)value);
#	    elif NF_INT2_IS_C_INT
		status = nc_get_vara_int(ncid, varid, start, count,
					   (int*)value);
#	    elif NF_INT2_IS_C_LONG
		status = nc_get_vara_long(ncid, varid, start, count,
					   (long*)value);
#	    endif
	    break;
	case NC_INT:
#	    if NF_INT_IS_C_INT
		status = nc_get_vara_int(ncid, varid, start, count,
					   (int*)value);
#	    elif NF_INT_IS_C_LONG
		status = nc_get_vara_long(ncid, varid, start, count,
					   (long*)value);
#	    endif
	    break;
	case NC_FLOAT:
#	    if NF_REAL_IS_C_FLOAT
		status = nc_get_vara_float(ncid, varid, start, count,
					   (float*)value);
#	    elif NF_REAL_IS_C_DOUBLE
		status = nc_get_vara_double(ncid, varid, start, count,
					   (double*)value);
#	    endif
	    break;
	case NC_DOUBLE:
#	    if NF_DOUBLEPRECISION_IS_C_FLOAT
		status = nc_get_vara_float(ncid, varid, start, count,
					   (float*)value);
#	    elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		status = nc_get_vara_double(ncid, varid, start, count,
					   (double*)value);
#	    endif
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVGT", status, "");
	*rcode = ncerr;
    }
}


/*
 * Read a hypercube of character values from a netCDF variable.
 */
extern void
c_ncvgtc(
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    char*		value,	/* block of data values to be read */
    int			lenstr,	/* declared length of the data argument */
    int*		rcode	/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0)
    {
	if (datatype != NC_CHAR)
	    status = NC_ECHAR;
	else if ((status = nc_get_vara_text(ncid, varid, start, count, value))
		 == 0)
	{
	    int	rank;

	    if ((status = nc_inq_varndims(ncid, varid, &rank)) == 0)
	    {
		size_t	total = dimprod(count, rank);

		(void) memset(value+total, ' ', lenstr - total);
	    }
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVGTC", status, "");
	*rcode = ncerr;
    }
}

/*
 * Read a generalized hypercube of numeric values from a netCDF variable of an 
 * open netCDF file.
 */
extern void
c_ncvgtg (
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    const ptrdiff_t*	strides,/* netCDF variable access strides */
    const ptrdiff_t*	imap,	/* memory values access basis vector */
    void*		value,	/* block of data values to be read */
    int*		rcode	/* returned error code */
)
{
    int		status;
    int		rank;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0 &&
	(status = nc_inq_varndims(ncid, varid, &rank)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = NC_ECHAR;
	    break;
	case NC_BYTE:
#	    if NF_INT1_IS_C_SIGNED_CHAR
		status = nc_get_varm_schar(ncid, varid, start, count,
					   strides, imap,
					   (signed char*)value);
#	    elif NF_INT1_IS_C_SHORT
		status = nc_get_varm_short(ncid, varid, start, count,
					   strides, imap,
					   (short*)value);
#	    elif NF_INT1_IS_C_INT
		status = nc_get_varm_int(ncid, varid, start, count,
					   strides, imap,
					   (int*)value);
#	    elif NF_INT1_IS_C_LONG
		status = nc_get_varm_long(ncid, varid, start, count,
					   strides, imap,
					   (long*)value);
#	    endif
	    break;
	case NC_SHORT:
#	    if NF_INT2_IS_C_SHORT
		status = nc_get_varm_short(ncid, varid, start, count,
					   strides, imap,
					   (short*)value);
#	    elif NF_INT2_IS_C_INT
		status = nc_get_varm_int(ncid, varid, start, count,
					   strides, imap,
					   (int*)value);
#	    elif NF_INT2_IS_C_LONG
		status = nc_get_varm_long(ncid, varid, start, count,
					   strides, imap,
					   (long*)value);
#	    endif
	    break;
	case NC_INT:
#	    if NF_INT_IS_C_INT
		status = nc_get_varm_int(ncid, varid, start, count,
					   strides, imap,
					   (int*)value);
#	    elif NF_INT_IS_C_LONG
		status = nc_get_varm_long(ncid, varid, start, count,
					   strides, imap,
					   (long*)value);
#	    endif
	    break;
	case NC_FLOAT:
#	    if NF_REAL_IS_C_FLOAT
		status = nc_get_varm_float(ncid, varid, start, count,
					   strides, imap,
					   (float*)value);
#	    elif NF_REAL_IS_C_DOUBLE
		status = nc_get_varm_double(ncid, varid, start, count,
					   strides, imap,
					   (double*)value);
#	    endif
	    break;
	case NC_DOUBLE:
#	    if NF_DOUBLEPRECISION_IS_C_FLOAT
		status = nc_get_varm_float(ncid, varid, start, count,
					   strides, imap,
					   (float*)value);
#	    elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		status = nc_get_varm_double(ncid, varid, start, count,
					   strides, imap,
					   (double*)value);
#	    endif
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVGTG", status, "");
	*rcode = ncerr;
    }
}

/*
 * Read a generalized hypercube of character values from a netCDF variable 
 * of an open netCDF file.
 */
extern void
c_ncvggc(
    int			ncid,	/* netCDF ID */
    int			varid,	/* variable ID */
    const size_t*	start,	/* multidimensional index of hypercube corner */
    const size_t*	count,	/* multidimensional hypercube edge lengths */
    const ptrdiff_t*	strides,/* netCDF variable access strides */
    const ptrdiff_t*	imap,	/* memory values access basis vector */
    char*		value,	/* block of data values to be written */
    int*		rcode	/* returned error code */
)
{
    int		status;
    int		rank;
    nc_type	datatype;

    if ((status = nc_inq_vartype(ncid, varid, &datatype)) == 0 &&
	(status = nc_inq_varndims(ncid, varid, &rank)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = nc_get_varm_text(ncid, varid, start, count,
				       strides, imap,
				       value);
	    break;
	default:
	    status = NC_ECHAR;
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCVGGC", status, "");
	*rcode = ncerr;
    }
}

/*
 * Change the name of a netCDF variable in an open netCDF file.
 */
extern void
c_ncvren (
    int ncid,		/* netCDF ID */
    int varid,		/* variable ID */
    const char* varname,/* new name for variable */
    int* rcode		/* returned error code */
)
{
    *rcode = ncvarrename (ncid, varid, varname) == -1
		? ncerr
		: 0;
}

/*
 * Add or changes a numeric variable or global attribute of an open
 * netCDF file.
 */
extern void
c_ncapt (
    int		ncid,		/* netCDF ID */
    int		varid,		/* variable ID */
    const char*	attname,	/* attribute name */
    nc_type	datatype,	/* attribute datatype */
    size_t	attlen,		/* attribute length */
    const void*	value,		/* pointer to data values */
    int*	rcode		/* returned error code */
)
{
    int		status;

    switch (datatype)
    {
    case NC_CHAR:
	status = NC_ECHAR;
	break;
    case NC_BYTE:
#	if NF_INT1_IS_C_SIGNED_CHAR
	    status = nc_put_att_schar(ncid, varid, attname, datatype,
				       attlen, (const signed char*)value);
#	elif NF_INT1_IS_C_SHORT
	    status = nc_put_att_short(ncid, varid, attname, datatype,
				       attlen, (const short*)value);
#	elif NF_INT1_IS_C_INT
	    status = nc_put_att_int(ncid, varid, attname, datatype,
				       attlen, (const int*)value);
#	elif NF_INT1_IS_C_LONG
	    status = nc_put_att_long(ncid, varid, attname, datatype,
				       attlen, (const long*)value);
#	endif
	break;
    case NC_SHORT:
#	if NF_INT2_IS_C_SHORT
	    status = nc_put_att_short(ncid, varid, attname, datatype,
				       attlen, (const short*)value);
#	elif NF_INT2_IS_C_INT
	    status = nc_put_att_int(ncid, varid, attname, datatype,
				       attlen, (const int*)value);
#	elif NF_INT2_IS_C_LONG
	    status = nc_put_att_long(ncid, varid, attname, datatype,
				       attlen, (const long*)value);
#	endif
	break;
    case NC_INT:
#	if NF_INT_IS_C_INT
	    status = nc_put_att_int(ncid, varid, attname, datatype,
				       attlen, (const int*)value);
#	elif NF_INT_IS_C_LONG
	    status = nc_put_att_long(ncid, varid, attname, datatype,
				       attlen, (const long*)value);
#	endif
	break;
    case NC_FLOAT:
#	if NF_REAL_IS_C_FLOAT
	    status = nc_put_att_float(ncid, varid, attname, datatype,
				       attlen, (const float*)value);
#	elif NF_REAL_IS_C_DOUBLE
	    status = nc_put_att_double(ncid, varid, attname, datatype,
				       attlen, (const double*)value);
#	endif
	break;
    case NC_DOUBLE:
#	if NF_DOUBLEPRECISION_IS_C_FLOAT
	    status = nc_put_att_float(ncid, varid, attname, datatype,
				       attlen, (const float*)value);
#	elif NF_DOUBLEPRECISION_IS_C_DOUBLE
	    status = nc_put_att_double(ncid, varid, attname, datatype,
				       attlen, (const double*)value);
#	endif
	break;
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCAPT", status, "");
	*rcode = ncerr;
    }
}

/*
 * Add or change a character attribute of an open netCDF file.
 */
extern void
c_ncaptc(
    int		ncid,		/* netCDF ID */
    int		varid,		/* variable ID */
    const char*	attname,	/* attribute name */
    nc_type	datatype,	/* attribute datatype */
    size_t	attlen,		/* attribute length */
    const char*	value,		/* pointer to data values */
    int*	rcode		/* returned error code */
)
{
    int		status;

    if (datatype != NC_CHAR)
	status = NC_ECHAR;
    else
	status = nc_put_att_text(ncid, varid, attname, attlen, value);

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCAPTC", status, "");
	*rcode = ncerr;
    }
}

/*
 * Return information about a netCDF attribute given its variable
 * ID and name.
 */
extern void
c_ncainq (
    int ncid,			/* netCDF ID */
    int varid,			/* variable ID */
    const char* attname,	/* attribute name */
    nc_type* datatype,		/* returned attribute datatype */
    int* attlen,		/* returned attribute length */
    int* rcode			/* returned error code */
)
{
    *rcode = ncattinq(ncid, varid, attname, datatype, attlen)
	     == -1
		? ncerr
		: 0;
}

/*
 * Get the value of a netCDF attribute given its variable ID and name.
 */
extern void
c_ncagt(
    int		ncid,		/* netCDF ID */
    int		varid,		/* variable ID */
    const char*	attname,	/* attribute name */
    void*	value,		/* pointer to data values */
    int*	rcode		/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_atttype(ncid, varid, attname, &datatype)) == 0)
    {
	switch (datatype)
	{
	case NC_CHAR:
	    status = NC_ECHAR;
	    break;
	case NC_BYTE:
#    	if NF_INT1_IS_C_SIGNED_CHAR
		status = nc_get_att_schar(ncid, varid, attname, 
					   (signed char*)value);
#    	elif NF_INT1_IS_C_SHORT
		status = nc_get_att_short(ncid, varid, attname, 
					   (short*)value);
#    	elif NF_INT1_IS_C_INT
		status = nc_get_att_int(ncid, varid, attname, 
					   (int*)value);
#    	elif NF_INT1_IS_C_LONG
		status = nc_get_att_long(ncid, varid, attname, 
					   (long*)value);
#    	endif
	    break;
	case NC_SHORT:
#    	if NF_INT2_IS_C_SHORT
		status = nc_get_att_short(ncid, varid, attname, 
					   (short*)value);
#    	elif NF_INT2_IS_C_INT
		status = nc_get_att_int(ncid, varid, attname, 
					   (int*)value);
#    	elif NF_INT2_IS_C_LONG
		status = nc_get_att_long(ncid, varid, attname, 
					   (long*)value);
#    	endif
	    break;
	case NC_INT:
#	    if NF_INT_IS_C_INT
		status = nc_get_att_int(ncid, varid, attname, 
					   (int*)value);
#	    elif NF_INT_IS_C_LONG
		status = nc_get_att_long(ncid, varid, attname, 
					   (long*)value);
#	    endif
	    break;
	case NC_FLOAT:
#    	if NF_REAL_IS_C_FLOAT
		status = nc_get_att_float(ncid, varid, attname, 
					   (float*)value);
#    	elif NF_REAL_IS_C_DOUBLE
		status = nc_get_att_double(ncid, varid, attname, 
					   (double*)value);
#    	endif
	    break;
	case NC_DOUBLE:
#    	if NF_DOUBLEPRECISION_IS_C_FLOAT
		status = nc_get_att_float(ncid, varid, attname, 
					   (float*)value);
#    	elif NF_DOUBLEPRECISION_IS_C_DOUBLE
		status = nc_get_att_double(ncid, varid, attname, 
					   (double*)value);
#    	endif
	    break;
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCAGT", status, "");
	*rcode = ncerr;
    }
}

/*
 * Get the value of a netCDF character attribute given its variable
 * ID and name.
 */
extern void
c_ncagtc(
    int		ncid,		/* netCDF ID */
    int		varid,		/* variable ID */
    const char*	attname,	/* attribute name */
    char*	value,		/* pointer to data values */
    int		attlen,		/* length of string argument */
    int*	rcode		/* returned error code */
)
{
    int		status;
    nc_type	datatype;

    if ((status = nc_inq_atttype(ncid, varid, attname, &datatype)) == 0)
    {
	if (datatype != NC_CHAR)
	    status = NC_ECHAR;
	else
	{
	    size_t	len;

	    status = nc_inq_attlen(ncid, varid, attname, &len);
	    if (status == 0)
	    {
		if (attlen < len)
		    status = NC_ESTS;
		else
		{
		    status = nc_get_att_text(ncid, varid, attname, 
					       value);
		    if (status == 0)
			(void) memset(value+len, ' ', attlen - len);
		}
	    }
	}
    }

    if (status == 0)
	*rcode = 0;
    else
    {
	nc_advise("NCAGTC", status, "");
	*rcode = ncerr;
    }
}

/*
 * Copy an attribute from one open netCDF file to another.
 */
extern void
c_ncacpy (
    int inncid,		/* input netCDF ID */
    int invarid,	/* variable ID of input netCDF or NC_GLOBAL */
    const char* attname,/* name of attribute in input netCDF to be copied */
    int outncid,	/* ID of output netCDF file for attribute */
    int outvarid,	/* ID of associated netCDF variable or NC_GLOBAL */
    int* rcode		/* returned error code */
)
{
    *rcode = ncattcopy(inncid, invarid, attname, outncid, outvarid)
	     == -1
		? ncerr
		: 0;
}

/*
 * Get the name of an attribute given its variable ID and number
 * as an attribute of that variable.
 */
extern void
c_ncanam (
    int ncid,		/* netCDF ID */
    int varid,		/* variable ID */
    int attnum,		/* attribute number */
    char* attname,	/* returned attribute name */
    int* rcode		/* returned error code */
)
{
    *rcode = ncattname(ncid, varid, attnum, attname) == -1
		? ncerr
		: 0;
}

/*
 * Rename an attribute in an open netCDF file.
 */
extern void
c_ncaren (
    int ncid,		/* netCDF ID */
    int varid,		/* variable ID */
    const char* attname,/* attribute name */
    const char* newname,/* new name */
    int* rcode		/* returned error code */
)
{
    *rcode = ncattrename(ncid, varid, attname, newname) == -1
		? ncerr
		: 0;
}

/*
 * Delete an attribute from an open netCDF file given the attribute name.
 */
extern void
c_ncadel (
    int ncid,		/* netCDF ID */
    int varid,		/* variable ID */
    const char* attname,/* attribute name */
    int* rcode		/* returned error code */
)
{
    *rcode = ncattdel(ncid, varid, attname) == -1
		? ncerr
		: 0;
}

/*
 * Set the fill mode of a netCDF file open for writing.
 */
extern int
c_ncsfil (
    int ncid,		/* netCDF ID */
    int fillmode,	/* fill mode, NCNOFILL or NCFILL */
    int* rcode		/* returned error code */
)
{
    int retval;

    *rcode = ((retval = ncsetfill(ncid, fillmode)) == -1)
		? ncerr
		: 0;

    return retval;
}

#endif /*!NO_NETCDF_2*/
