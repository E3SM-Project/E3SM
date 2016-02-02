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
                                 ZOLTAN2RCB    = 5, & !Geometric Recursive Bisection
                                 ZOLTAN2MJ    = 6, &  !Geometric Multijagged
                                 ZOLTAN2RIB    = 7, & !Recursive Inertial Bisection
                                 ZOLTAN2HSFC    = 8, & !Hilbert Space Filling Curve
                                 ZOLTAN2PATOH    = 9, & !Patoh Hypergraph partitioner->requires PATOH TPL
                                 ZOLTAN2PHG   = 10, & !PHG Hypergraph partitioner
                                 ZOLTAN2METIS    = 11, & !Metis Graph Partitioner -> requires metis TPL
                                 ZOLTAN2PARMETIS    = 12, & !Parmetis -> requires parmetis TPL
                                 ZOLTAN2PARMA    = 13, & !Parma Mesh Partitioner-> requires parma tpl
                                 ZOLTAN2SCOTCH    = 14, & !Scotch graph partitioner -> requires scotch tpl
                                 ZOLTAN2PTSCOTCH    = 15, & !PTscotch graph partitioner -> requires scotch tpl
                                 ZOLTAN2BLOCK    = 16, & 
                                 ZOLTAN2CYCLIC    = 17, &
                                 ZOLTAN2RANDOM    = 18, &
                                 ZOLTAN2ZOLTAN    = 19, &
                                 ZOLTAN2ND    = 20

end module params_mod
