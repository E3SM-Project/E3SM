
module prim_advance_mod
  use prim_advance_mod_base, only: prim_advance_exp, prim_advance_init1, applyCAMforcing_dynamics, applyCAMforcing, &
                                   vertical_mesh_init2, edge3p1, ur_weights, set_prescribed_wind
  implicit none
  public :: prim_advance_exp, prim_advance_init1, applyCAMforcing_dynamics, applyCAMforcing, &
            vertical_mesh_init2, edge3p1, ur_weights, set_prescribed_wind
end module prim_advance_mod
