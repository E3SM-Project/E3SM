/*
 * PnetCDF C Test for version string
 */
#include <stdio.h>
#include "pnetcdf.h"

int main()
{
#ifdef PNETCDF_VERSION_MAJOR
	int major = PNETCDF_VERSION_MAJOR;
	int minor = PNETCDF_VERSION_MINOR;
	int sub = PNETCDF_VERSION_SUB;

	printf("%u.%u.%u",major,minor,sub);
	return 0;
#else
	XXX;
#endif
}

