#ifdef _PNETCDF
#include <pnetcdf.h>

#ifdef FORTRANUNDERSCORE
#define pnetcdf_version_check pnetcdf_version_check_
#endif
#ifdef FORTRANCAPS
#define pnetcdf_version_check PNETCDF_VERSION_CHECK
#endif
#ifdef FORTRAN_GNUF2C
#define pnetcdf_version_check pnetcdf_version_check__
#endif



int pnetcdf_version_check()
{
  if(PNETCDF_VERSION_MAJOR >= 1 && PNETCDF_VERSION_MINOR >= 1){
    return 1;
  }else{
    return 0;
  }
}
#endif
