#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advance_mod

  use control_mod,           only: qsplit,rsplit, use_moisture
  use derivative_mod,        only: derivative_t
  use dimensions_mod,        only: np, nlev, nlevp, nelemd, qsize, max_corner_elem
  use edgetype_mod,          only: EdgeDescriptor_t, EdgeBuffer_t
  use edge_mod,              only: edge_g, edgevpack_nlyr, edgevunpack_nlyr
  use element_mod,           only: element_t
  use hybrid_mod,            only: hybrid_t
  use hybvcoord_mod,         only: hvcoord_t
  use kinds,                 only: real_kind, iulog
  use perf_mod,              only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,          only: abortmp, parallel_t, iam
  use time_mod,              only: timelevel_t
  use prim_advance_mod_base, only: prim_advance_exp, prim_advance_init1, &
                                   applyCAMforcing_dynamics, applyCAMforcing, vertical_mesh_init2, set_prescribed_wind,&
                                   advance_hypervis_dp, compute_and_apply_rhs
  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init1, &
            applyCAMforcing_dynamics, applyCAMforcing, vertical_mesh_init2, set_prescribed_wind,&
            advance_hypervis_dp, compute_and_apply_rhs

contains



end module prim_advance_mod
