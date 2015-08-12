/*
 * NetCDF C Test for Parallel Support
 */
#include "netcdf_meta.h"

int main()
{
#ifdef NC_HAS_PARALLEL
	return 0;
#else
	XXX;
#endif
}
