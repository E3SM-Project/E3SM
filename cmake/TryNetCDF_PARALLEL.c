/*
 * NetCDF C Test for Parallel Support
 */
#include "netcdf_meta.h"

int main()
{
#if NC_HAS_PARALLEL==1
	return 0;
#else
	XXX;
#endif
}
