/*
 * NetCDF C Test for DAP Support
 */
#include <stdio.h>
#include "netcdf_meta.h"

int main()
{
#ifdef NC_VERSION
	char version[20] = NC_VERSION;
	printf("%s",version);
	return 0;
#else
	XXX;
#endif
}
