#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod
  use edge_mod_base, only:&
    FreeEdgeBuffer,&
    FreeGhostBuffer3D,&
    FreeLongEdgeBuffer,&
    LongEdgeVpack,&
    LongEdgeVunpackMIN,&
    edgeDGVpack,&
    edgeDGVunpack,&
    edgeDefaultVal,&
    edgeSpack,&
    edgeSunpackMax,&
    edgeSunpackMin,&
    edgeVpack,&
    edgeVunpack,&
    edgeVunpackMAX,&
    edgeVunpackMIN,&
    edgeVunpackVert,&
    ghostVpack3d,&
    ghostVpack_unoriented,&
    ghostVpackfull,&
    ghostVunpack3d,&
    ghostVunpack_unoriented,&
    ghostVunpackfull,&
    initEdgeBuffer,&
    initEdgeSBuffer,&
    initGhostBuffer3D,&
    initLongEdgeBuffer

  implicit none
end module edge_mod
