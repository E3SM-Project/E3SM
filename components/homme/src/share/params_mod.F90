#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module params_mod
   integer, public, parameter :: INTERNAL_EDGE  = 0
   integer, public, parameter :: EXTERNAL_EDGE  = 1

   ! Mesh partitioning methods. The space-filling curve is the only one supported:
   ! the METIS and Zoltan(2) methods have been removed. The value 4 is kept for
   ! backward compatibility with existing namelists (partmethod/se_partmethod=4),
   ! and the old METIS values (0-3) are still mapped to SFCURVE in namelist_mod.
   integer, public, parameter :: SFCURVE = 4


end module params_mod
