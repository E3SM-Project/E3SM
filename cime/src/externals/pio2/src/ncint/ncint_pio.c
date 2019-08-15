/**
 * @file
 * @internal Additional nc_* functions to support netCDF integration.
 *
 * @author Ed Hartnett
 */

#include "config.h"
#include <stdlib.h>
#include <pio_internal.h>
#include "ncintdispatch.h"

/** This is te default io system id. */
extern int diosysid;

/** Have we initialized the netCDF integration layer? This is where we
 * register our dispatch layer with netcdf-c. */
extern int ncint_initialized;

/**
 * Same as PIOc_Init_Intracomm().
 *
 * @author Ed Hartnett
 */
int
nc_def_iosystemm(MPI_Comm comp_comm, int num_iotasks, int stride, int base,
                 int rearr, int *iosysidp)
{
    int ret;

    /* Call the PIOc_ function to initialize the intracomm. */
    if ((ret = PIOc_Init_Intracomm(comp_comm, num_iotasks, stride, base, rearr,
                                   iosysidp)))
        return ret;

    /* Remember the io system id. */
    diosysid = *iosysidp;

    return PIO_NOERR;
}

/**
 * Set the default iosystemID.
 *
 * @param iosysid The IO system ID to set.
 *
 * @return PIO_NOERR for success.
 * @author Ed Hartnett
 */
int
nc_set_iosystem(int iosysid)
{
    /* Remember the io system id. */
    diosysid = iosysid;

    return PIO_NOERR;
}

/**
 * Get the default iosystemID.
 *
 * @param iosysid Pointer that gets The IO system ID.
 *
 * @return PIO_NOERR for success.
 * @author Ed Hartnett
 */
int
nc_get_iosystem(int *iosysid)
{
    pioassert(iosysid, "pointer to iosysid must be provided", __FILE__,
              __LINE__);

    /* Remember the io system id. */
    *iosysid = diosysid;

    return PIO_NOERR;
}

/**
 * Same as PIOc_free_iosystem().
 *
 * @author Ed Hartnett
 */
int
nc_free_iosystem(int iosysid)
{
    return PIOc_free_iosystem(iosysid);
}

/**
 * Same as PIOc_init_decomp().
 *
 * @author Ed Hartnett
 */
int
nc_def_decomp(int iosysid, int pio_type, int ndims, const int *gdimlen,
               int maplen, const size_t *compmap, int *ioidp,
              int rearranger, const size_t *iostart,
              const size_t *iocount)
{
    return PIOc_init_decomp(iosysid, pio_type, ndims, gdimlen, maplen,
                            (const PIO_Offset *)compmap, ioidp, rearranger,
                            (const PIO_Offset *)iostart,
                            (const PIO_Offset *)iocount);
}

/**
 * Same as PIOc_freedecomp().
 *
 * @author Ed Hartnett
 */
int
nc_free_decomp(int ioid)
{
    PLOG((1, "nc_free_decomp ioid %d", ioid));
    return PIOc_freedecomp(diosysid, ioid);
}
