module scream_homme_interface_mod

  ! Import homme types to be used for module-level variables
  ! TODO: some of these are only needed cause homme's f90 interfaces expect them,
  !       but they are discarded and never used in scream (e.g., domain1d).
  !       We need to find a way to eliminate these useless dependencies
  ! NOTE: perhaps we should call the init_homme_f90_structures subroutine
  !       directly from atm_mct_init
  use domain_mod,   only: domain1d_t
  use element_mod,  only: element_t
  use parallel_mod, only: parallel_t
  use time_mod,     only: timelevel_t

  implicit none
  private 

  type (element_t), allocatable :: elem(:)
  type (domain1d_t)  :: dom_mt
  type (parallel_t)  :: par
  type (TimeLevel_t) :: tl

contains:

  subroutine init_homme_f90_structures () bind(c)
    use prim_driver_mod, only: prim_init1

    call prim_init1(elem,par,dom_mt,tl)
    call prim_init2(elem, hybrid, nets, nete, tl, hvcoord)
  end subroutine init_homme_f90_structures

end module scream_homme_interface_mod
