! Include Homme's config settings
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module homme_driver_mod
  use iso_c_binding,     only: c_ptr, c_int, c_double, c_bool, C_NULL_CHAR

  use parallel_mod,  only: abortmp

  implicit none
  private

  public :: prim_init_data_structures_f90
  public :: prim_complete_init1_phase_f90
  public :: prim_set_hvcoords_f90
  public :: prim_init_model_f90
  public :: prim_run_f90
  public :: prim_finalize_f90

contains

  subroutine prim_init_data_structures_f90 () bind(c)
    use prim_driver_mod,      only: prim_create_c_data_structures
    use prim_driver_base,     only: prim_init1_elem_arrays
    use prim_cxx_driver_base, only: setup_element_pointers
    use derivative_mod_base,  only: derivinit
    use time_mod,             only: TimeLevel_init
    use homme_context_mod,    only: is_geometry_inited, is_data_structures_inited, &
                                    par, elem, tl, deriv, hvcoord
    use kinds,                only: iulog

    if (is_data_structures_inited) then
      call abortmp ("Error! prim_init_data_structures_f90 was already called.\n")
    elseif (.not. is_geometry_inited) then
      call abortmp ("Error! 'homme_init_grids_f90 was not called yet.\n")
    endif

    if (par%masterproc) write(iulog, *) "Initing prim data structures..."

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv)

    ! ==================================
    ! Initialize element pointers
    ! ==================================
    call setup_element_pointers(elem)

    ! ==================================
    ! Initialize element arrays (fluxes and state)
    ! ==================================
    call prim_init1_elem_arrays(elem,par)

    ! Initialize the time levels
    call TimeLevel_init(tl)

    ! Initialize Kokkos data structures
    call prim_create_c_data_structures (tl, hvcoord, elem(1)%mp)

    is_data_structures_inited = .true.
  end subroutine prim_init_data_structures_f90

  subroutine prim_complete_init1_phase_f90 () bind(c)
    use prim_driver_base,  only: prim_init1_buffers, prim_init1_compose, prim_init1_cleanup
    use homme_context_mod, only: par, elem
    use compose_mod,       only: compose_control_kokkos_init_and_fin
    use prim_driver_mod,   only: prim_init_grid_views

    ! Compose is not in charge of init/finalize kokkos
    call compose_control_kokkos_init_and_fin(.false.)

    ! Init compose
    call prim_init1_compose(par,elem)

    ! ==================================
    ! Initialize the buffers for exchanges
    ! ==================================
    call prim_init1_buffers(elem,par)

    ! Cleanup the tmp stuff used in prim_init1_geometry
    call prim_init1_cleanup()

    call prim_init_grid_views (elem)

  end subroutine prim_complete_init1_phase_f90


  subroutine prim_set_hvcoords_f90 (ps0, hyai_ptr, hybi_ptr, hyam_ptr, hybm_ptr) bind(c)
    use iso_c_binding,     only: c_f_pointer
    use homme_context_mod, only: hvcoord, masterproc
    use dimensions_mod,    only: nlev, nlevp
    use hybvcoord_mod,     only: set_layer_locations
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: hyai_ptr, hybi_ptr, hyam_ptr, hybm_ptr
    real(kind=c_double), intent(in) :: ps0
    !
    ! Locals
    !
    real(kind=c_double), pointer :: hyai(:), hybi(:), hyam(:), hybm(:)

    call c_f_pointer (hyai_ptr, hyai, [nlevp])
    call c_f_pointer (hybi_ptr, hybi, [nlevp])
    call c_f_pointer (hyam_ptr, hyam, [nlev])
    call c_f_pointer (hybm_ptr, hybm, [nlev])

    hvcoord%ps0 = ps0
    hvcoord%hyai = hyai
    hvcoord%hybi = hybi
    hvcoord%hyam = hyam
    hvcoord%hybm = hybm

    call set_layer_locations(hvcoord,.true.,masterproc)
  end subroutine prim_set_hvcoords_f90

  subroutine prim_copy_cxx_to_f90 (copy_phis)
    use iso_c_binding,       only: c_ptr, c_loc
    use homme_context_mod,   only: tl, elem, deriv
    use dimensions_mod,      only: nlevp, nelemd, np
    use kinds,               only: real_kind
    use derivative_mod_base, only: gradient_sphere
    use theta_f2c_mod,       only: cxx_push_results_to_f90, init_geopotential_c
    use element_state,       only: elem_state_v, elem_state_w_i, elem_state_vtheta_dp,   &
                                   elem_state_phinh_i, elem_state_dp3d, elem_state_ps_v, &
                                   elem_state_Qdp, elem_state_Q, elem_derived_omega_p,   &
                                   elem_state_phis

    !
    ! Inputs
    !
    logical, intent(in) :: copy_phis
    !
    ! Locals
    !
    type (c_ptr) :: elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr, elem_state_phinh_i_ptr
    type (c_ptr) :: elem_state_dp3d_ptr, elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_derived_omega_p_ptr

    integer :: ie
    type (c_ptr) :: elem_state_phis_local_ptr, elem_derived_gradphis_local_ptr
    real (kind=real_kind), target, dimension(np,np)   :: elem_state_phis_local
    real (kind=real_kind), target, dimension(np,np,2) :: elem_gradphis_local

    ! Set pointers to states
    elem_state_v_ptr         = c_loc(elem_state_v)
    elem_state_w_i_ptr       = c_loc(elem_state_w_i)
    elem_state_vtheta_dp_ptr = c_loc(elem_state_vtheta_dp)
    elem_state_phinh_i_ptr   = c_loc(elem_state_phinh_i)
    elem_state_dp3d_ptr      = c_loc(elem_state_dp3d)
    elem_state_Qdp_ptr       = c_loc(elem_state_Qdp)
    elem_state_Q_ptr         = c_loc(elem_state_Q)
    elem_state_ps_v_ptr      = c_loc(elem_state_ps_v)
    elem_derived_omega_p_ptr = c_loc(elem_derived_omega_p)

    ! Copy data
    call cxx_push_results_to_f90(elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr,   &
                                 elem_state_phinh_i_ptr, elem_state_dp3d_ptr, elem_state_ps_v_ptr, &
                                 elem_state_Qdp_ptr, elem_state_Q_ptr, elem_derived_omega_p_ptr)
    if (copy_phis) then
      ! Set phis=phi(bottom)
      elem_state_phis(:,:,:) = elem_state_phinh_i(:,:,nlevp,tl%n0,:)

      ! Set geopotential views
      elem_state_phis_local_ptr       = c_loc(elem_state_phis_local)
      elem_derived_gradphis_local_ptr = c_loc(elem_gradphis_local)
      do ie=1,nelemd
        elem_state_phis_local = elem(ie)%state%phis
        elem_gradphis_local   = elem(ie)%derived%gradphis

        elem_gradphis_local = gradient_sphere(elem_state_phis_local, deriv, elem(ie)%Dinv)
        call init_geopotential_c (ie-1, elem_state_phis_local_ptr, elem_derived_gradphis_local_ptr)
      enddo
    endif
  end subroutine prim_copy_cxx_to_f90

  subroutine prim_init_model_f90 () bind(c)
    use prim_driver_mod,   only: prim_init_ref_states_views, &
                                 prim_init_diags_views, prim_init_kokkos_functors, &
                                 prim_init_state_views
    use prim_state_mod,    only: prim_printstate
    use model_init_mod,    only: model_init2
    use control_mod,       only: disable_diagnostics
    use dimensions_mod,    only: nelemd
    use homme_context_mod, only: is_model_inited, is_data_structures_inited, &
                                 elem, hybrid, hvcoord, deriv, tl

    ! Local variable
    logical(kind=c_bool), parameter :: allocate_buffer = .false.

    if (.not. is_data_structures_inited) then
      call abortmp ("Error! 'prim_init_data_structures_f90' has not been called yet.\n")
    elseif (is_model_inited) then
      call abortmp ("Error! 'prim_init_model_f90' has already been called.\n")
    endif

    ! Needed in model init, to have correct state values (or else EOS craps out)
    call prim_copy_cxx_to_f90 (.true.)

    ! Notably, this inits the ref states
    call model_init2(elem,hybrid,deriv,hvcoord,tl,1,nelemd)

    ! Initialize the C++ functors in the C++ context
    ! Here we set allocate_buffer=false since the AD
    ! allocates local memory for each process from a
    ! single buffer.
    call prim_init_kokkos_functors (allocate_buffer)

    ! Init ref_states views, and diags views
    call prim_init_ref_states_views (elem)
    call prim_init_diags_views (elem)

    ! In order to print up to date stuff in F90
    call prim_copy_cxx_to_f90 (.true.)
    if (.not. disable_diagnostics) then
      call prim_printstate(elem, tl, hybrid,hvcoord,1, nelemd)
    endif

    is_model_inited = .true.

  end subroutine prim_init_model_f90

  subroutine prim_run_f90 (nsplit_iteration) bind(c)
    use dimensions_mod,    only: nelemd
    use prim_driver_mod,   only: prim_run_subcycle
    use time_mod,          only: tstep
    use homme_context_mod, only: is_model_inited, elem, hybrid, tl, hvcoord

    integer(kind=c_int), value, intent(in) :: nsplit_iteration

    if (.not. is_model_inited) then
      call abortmp ("Error! prim_init_model_f90 was not called yet (or prim_finalize_f90 was already called).\n")
    endif

    if (tstep .le. 0) then
      call abortmp ("Error! No time step was set in Homme yet.\n")
    endif

    call prim_run_subcycle(elem,hybrid,1,nelemd,tstep,.false.,tl,hvcoord,nsplit_iteration)
  end subroutine prim_run_f90

  subroutine prim_finalize_f90 () bind(c)
    use homme_context_mod,    only: is_model_inited, elem, dom_mt, close_homme_log
    use prim_cxx_driver_base, only: prim_finalize

    if (.not. is_model_inited) then
      call abortmp ("Error! prim_init_model_f90 was not called yet (or prim_finalize_f90 was already called).\n")
    endif

    ! Cleanup some f90 stuff in Homme, and cleanup all cxx structures
    call prim_finalize()

    ! Deallocate the pointers
    deallocate (elem)
    deallocate (dom_mt)

    is_model_inited = .false.

    call close_homme_log()

  end subroutine prim_finalize_f90

end module homme_driver_mod
