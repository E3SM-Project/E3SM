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
                                 ZOLTAN2_1PHASEMAP       = 5, & !Zoltan2 Single Phase Mapper&Partitioner.
                                 ZOLTAN2RCB    = 6, &           !Geometric Recursive Bisection
                                 ZOLTAN2MJ    = 7, &            !Geometric Multijagged
                                 ZOLTAN2MJRCB    = 8, &         !Geometric Multijagged runned as recursive bisection.
                                 ZOLTAN2RIB    = 9, &           !Recursive Inertial Bisection
                                 ZOLTAN2HSFC   = 10, &          !Hilbert Space Filling Curve
                                 ZOLTAN2PATOH    = 11, &        !Patoh Hypergraph partitioner->requires PATOH TPL
                                 ZOLTAN2PHG   = 12, &           !PHG Hypergraph partitioner
                                 ZOLTAN2METIS    = 13, &        !Metis Graph Partitioner -> requires metis TPL
                                 ZOLTAN2PARMETIS    = 14, &     !Parmetis -> requires parmetis TPL
                                 ZOLTAN2PARMA    = 15, &        !Parma Mesh Partitioner-> requires parma tpl
                                 ZOLTAN2SCOTCH    = 16, &       !Scotch graph partitioner -> requires scotch tpl
                                 ZOLTAN2PTSCOTCH    = 17, &     !PTscotch graph partitioner -> requires scotch tpl
                                 ZOLTAN2BLOCK    = 18, &
                                 ZOLTAN2CYCLIC    = 19, &
                                 ZOLTAN2RANDOM    = 20, &
                                 ZOLTAN2ZOLTAN    = 21, &
                                 ZOLTAN2ND    = 22


   integer, public, parameter :: SPHERE_COORDS = 1, &
                                 CUBE_COORDS = 2, &
                                 FACE_2D_LB_COORDS = 3


   integer, public, parameter :: Z2_NO_TASK_MAPPING = 1, &
                                 Z2_TASK_MAPPING = 2, &
                                 Z2_OPTIMIZED_TASK_MAPPING = 3

end module params_mod
