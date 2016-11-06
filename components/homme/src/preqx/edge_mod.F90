#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod
  use edge_mod_base, only: initLongEdgeBuffer, FreeLongEdgeBuffer, LongEdgeVpack, LongEdgeVunpackMIN, initEdgeBuffer, initEdgeSBuffer, FreeEdgeBuffer, edgeVpack, edgeVunpack,       &
                           edgeVunpackMIN, edgeVunpackMAX, edgeDGVpack, edgeDGVunpack, edgeVunpackVert, edgeDefaultVal, initGhostBuffer3D, FreeGhostBuffer3D, &
                           ghostVpackfull, ghostVunpackfull, ghostVpack_unoriented, ghostVunpack_unoriented, ghostVpack3d, ghostVunpack3d, initGhostBufferTR, FreeGhostBufferTR,     &
                           ghostVpack, ghostVunpack, ghostVpackR, ghostVunpackR, ghostVpack2d, ghostVunpack2d, ghostVpack2d_single, ghostVunpack2d_single, ghostVpack2d_level,       &
                           ghostVunpack2d_level, edgeSpack, edgeSunpackMin, edgeSunpackMax
  implicit none
end module edge_mod
