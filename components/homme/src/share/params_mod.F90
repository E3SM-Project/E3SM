#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module params_mod
   integer, public, parameter :: INTERNAL_EDGE  = 0
   integer, public, parameter :: EXTERNAL_EDGE  = 1

   integer, public, parameter :: RECURSIVE  = 0, &   ! Type of partitioning methods 
                                 KWAY       = 1, &
                                 VOLUME     = 2, &
                                 WRECURSIVE = 3, &
                                 SFCURVE    = 4

end module params_mod
