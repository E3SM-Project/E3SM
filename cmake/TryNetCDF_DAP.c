/*
 * NetCDF C Test for DAP Support
 */
#include "netcdf_meta.h"

int main()
{
#ifdef NC_HAS_DAP
	return 0;
#else
	XXX;
#endif
}
