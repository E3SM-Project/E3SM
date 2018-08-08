#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod
  use edge_mod_base, only: initLongEdgeBuffer, FreeLongEdgeBuffer, LongEdgeVpack, LongEdgeVunpackMIN, initEdgeBuffer, initEdgeSBuffer, FreeEdgeBuffer, &
                           edgeVpack, edgeVunpack, edgeVpack_nlyr, edgeVunpack_nlyr,       &
                           edgeVunpackMIN, edgeVunpackMAX, edgeDGVpack, edgeDGVunpack, edgeVunpackVert, edgeDefaultVal, initGhostBuffer3D, FreeGhostBuffer3D, &
                           ghostVpackfull, ghostVunpackfull, ghostVpack_unoriented, ghostVunpack_unoriented, ghostVpack3d, ghostVunpack3d, &
                           edgeSpack, edgeSunpackMin, edgeSunpackMax, edge_g
  implicit none
end module edge_mod
