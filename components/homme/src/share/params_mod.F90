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
                                 SFCURVE    = 4, &
                                 ZOLTAN2RCB    = 5, &
                                 ZOLTAN2MJ    = 6, &
                                 ZOLTAN2RIB    = 7, &
                                 ZOLTAN2HSFC    = 8, &
                                 ZOLTAN2PATOH    = 9, &
                                 ZOLTAN2PHG   = 10, &
                                 ZOLTAN2METIS    = 11, &
                                 ZOLTAN2PARMETIS    = 12, &
                                 ZOLTAN2PARMA    = 13, &
                                 ZOLTAN2SCOTCH    = 14, &
                                 ZOLTAN2PTSCOTCH    = 15, &
                                 ZOLTAN2BLOCK    = 16, &
                                 ZOLTAN2CYCLIC    = 17, &
                                 ZOLTAN2RANDOM    = 18, &
                                 ZOLTAN2ZOLTAN    = 19, &
                                 ZOLTAN2ND    = 20

end module params_mod
