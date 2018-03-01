/*

Copyright 2006, University Corporation for Atmospheric Research. See
the COPYRIGHT file for copying and redistribution conditions.

$Id: fort-lib.c,v 1.15 2009/02/13 15:58:00 ed Exp $
*/

/*
 Extracted from of fort-lib.c. Used to supply four required netcdf4
 functions used in Fortran2003 interfaces. These need to be external
 and I don't want to mangle fort-lib.c just to define these four functions

 Version 1.0, April, 2009 - based on netcdf-4.0.1 source
 Version 2.0, April, 2010 - based on netcdf-4.1.1 source

 modified by: Richard Weed Ph.D.
              Center for Advanced Vehicular Systems
              Mississippi State University
              rweed@cavs.msstate.edu
*/
    
#include <config.h>
#include <stddef.h>	/* for NULL */
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "netcdf.h"

#ifdef USE_NETCDF4
/* These appear to only be defined in netcdf-4*/

/* Get the varids for a fortran function (i.e. add 1 to each
 * varid.) */
extern int
nc_inq_varids_f(int ncid, int *nvars, int *fvarids)
{
   int *varids, nvars1;
   int i, ret = NC_NOERR;

   /* Get the information from the C library. */
   if ((ret = nc_inq_varids(ncid, &nvars1, NULL)))
      return ret;
   if (!(varids = malloc(nvars1 * sizeof(int))))
      return NC_ENOMEM;
   if ((ret = nc_inq_varids(ncid, NULL, varids)))
      goto exit;

   /* Add one to each, for fortran. */
   for (i = 0; i < nvars1; i++)
      fvarids[i] = varids[i] + 1;

   /* Tell the user how many there are. */
   if (nvars)
      *nvars = nvars1;

  exit:
   free(varids);
   return ret;
}
/* Get the dimids for a fortran function (i.e. add 1 to each
 * dimid.) */
extern int
nc_inq_dimids_f(int ncid, int *ndims, int *fdimids, int parent)
{
   int *dimids, ndims1;
   int i, ret = NC_NOERR;

   /* Get the information from the C library. */
   if ((ret = nc_inq_dimids(ncid, &ndims1, NULL, parent)))
      return ret;
   if (!(dimids = malloc(ndims1 * sizeof(int))))
      return NC_ENOMEM;
   if ((ret = nc_inq_dimids(ncid, NULL, dimids, parent)))
      goto exit;

   /* Add one to each, for fortran. */
   for (i = 0; i < ndims1; i++)
      fdimids[i] = dimids[i] + 1;

   /* Tell the user how many there are. */
   if (ndims)
      *ndims = ndims1;

  exit:
   free(dimids);
   return ret;
}

/* Swap the dim sizes for fortran. */
extern int
nc_insert_array_compound_f(int ncid, int typeid, char *name, 
			 size_t offset, nc_type field_typeid,
			 int ndims, int *dim_sizesp)
{
   int *dim_sizes_f;
   int i, ret;

   if (ndims <= 0)
      return NC_EINVAL;

   /* Allocate some storage to hold ids. */
   if (!(dim_sizes_f = malloc(ndims * sizeof(int))))
      return NC_ENOMEM;

   /* Create a backwards list of dimension sizes. */
   for (i = 0; i < ndims; i++)
      dim_sizes_f[i] = dim_sizesp[ndims - i - 1];

   /* Call with backwards list. */
   ret = nc_insert_array_compound(ncid, typeid, name, offset, field_typeid, 
				  ndims, dim_sizes_f);

   /* Clean up. */
   free(dim_sizes_f);
   return ret;
}

extern int
nc_inq_compound_field_f(int ncid, nc_type xtype, int fieldid, char *name, 
			size_t *offsetp, nc_type *field_typeidp, int *ndimsp, 
			int *dim_sizesp)
{
   int ndims;
   int ret;

   /* Find out how many dims. */
   if ((ret = nc_inq_compound_field(ncid, xtype, fieldid, NULL, NULL, 
				    NULL, &ndims, NULL)))
      return ret;

   /* Call the function. */
   if ((ret = nc_inq_compound_field(ncid, xtype, fieldid, name, offsetp, 
				    field_typeidp, ndimsp, dim_sizesp)))
      return ret;

   /* Swap the order of the dimsizes. */
   if (ndims)
   {
      int *f, *b, temp;
      for (f = dim_sizesp, b = &dim_sizesp[ndims - 1]; f < b; f++, b--)
      {
	 temp = *f;
	 *f = *b;
	 *b = temp;
      }
   }  

   return NC_NOERR;
}

#endif /*USE_NETCDF4*/

/*
 add a dummy nc_rename_grp function if it is not supported. This is include
 here so we can build/test with netCDF < version 4.3.1 without 
*/

#ifndef NC_HAVE_RENAME_GRP
extern int
nc_rename_grp(int ncid, const char *name)
{
 printf("\n*** Warning - nc_rename_grp not supported in this netCDF version\n");
 printf("*** Update your netCDF C libraries to version 4.3.1 or higher\n");

 return NC_ENOGRP;

}
#endif

