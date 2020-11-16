#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module homme_driver_mod
  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool, C_NULL_CHAR

  use parallel_mod,  only: abortmp
  use perf_mod,      only: t_initf, t_finalizef, t_prf, t_startf, t_stopf

  implicit none
  private 

  public :: prim_init_data_structures_f90
  public :: prim_copy_cxx_to_f90
  public :: prim_init_model_f90
  public :: prim_run_f90
  public :: prim_finalize_f90

contains

  subroutine prim_init_data_structures_f90 () bind(c)
    use prim_driver_mod,      only: prim_create_c_data_structures, prim_init_kokkos_functors
    use prim_driver_base,     only: prim_init1_geometry, prim_init1_elem_arrays, &
                                    prim_init1_cleanup, prim_init1_buffers,      &
                                    MetaVertex, GridEdge
    use prim_cxx_driver_base, only: setup_element_pointers
    use derivative_mod_base,  only: derivinit
    use time_mod,             only: TimeLevel_init
    use hybvcoord_mod,        only: hvcoord_init
    use hybrid_mod,           only: hybrid_create
    use control_mod,          only: vfile_mid, vfile_int
    use homme_context_mod,    only: is_parallel_inited, is_geometry_inited, is_data_structures_inited, &
                                    par, dom_mt, elem, tl, deriv, hybrid, hvcoord, init_parallel_f90
    use dimensions_mod,       only: nelemd
    !
    ! Local(s)
    !
    integer :: ierr

    if (is_data_structures_inited) then
      call abortmp ("Error! prim_init_data_structures_f90 was already called.\n")
    elseif (.not. is_geometry_inited) then
      call abortmp ("Error! 'homme_init_geometry_f90 was not called yet.\n")
    endif

    print *, "Initing prim data structures..."

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv)

    ! ==================================
    ! Initialize the buffers for exchanges
    ! ==================================
    call prim_init1_buffers(elem,par)

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

    ! Init hvcoord
    hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%masterthread, ierr)

    ! Initialize Kokkos data structures
    call prim_create_c_data_structures (tl, hvcoord, elem(1)%mp)

    is_data_structures_inited = .true.
  end subroutine prim_init_data_structures_f90

  subroutine prim_copy_cxx_to_f90 () bind (c)
    use iso_c_binding, only: c_ptr, c_loc
    use theta_f2c_mod, only: cxx_push_results_to_f90
    use element_state, only: elem_state_v, elem_state_w_i, elem_state_vtheta_dp,     &
                             elem_state_phinh_i, elem_state_dp3d, elem_state_ps_v,   &
                             elem_state_Qdp, elem_state_Q, elem_derived_omega_p
    !
    ! Locals
    !
    type (c_ptr) :: elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr, elem_state_phinh_i_ptr
    type (c_ptr) :: elem_state_dp3d_ptr, elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_derived_omega_p_ptr
    
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
  end subroutine prim_copy_cxx_to_f90

  subroutine prim_init_model_f90 (standalone) bind(c)
    use prim_driver_mod,   only: prim_init_grid_views, prim_init_ref_states_views, &
                                 prim_init_diags_views, prim_init_kokkos_functors, &
                                 prim_init_state_views
  
    use prim_state_mod,    only: prim_printstate
    use model_init_mod,    only: model_init2
    use dimensions_mod,    only: nelemd
    use homme_context_mod, only: is_model_inited, is_data_structures_inited, &
                                 elem, hybrid, hvcoord, deriv, tl
    !
    ! Input(s)
    !
    logical(kind=c_bool), intent(in) :: standalone

    if (.not. is_data_structures_inited) then
      call abortmp ("Error! 'prim_init_data_structures_f90' has not been called yet.\n")
    elseif (is_model_inited) then
      call abortmp ("Error! 'prim_init_model_f90' has already been called.\n")
    endif

    ! Notably, this inits the ref states
    call model_init2(elem,hybrid,deriv,hvcoord,tl,1,nelemd)

    ! Initialize the C++ functors in the C++ context
    call prim_init_kokkos_functors ()

    ! Init grid views, ref_states views, and diags views
    call prim_init_grid_views (elem)
    call prim_init_ref_states_views (elem)
    call prim_init_diags_views (elem)

    ! If this is a standalone run, we need to copy initial f90 states to C++
    if (standalone) then
      call prim_init_state_views(elem)
    endif

    call prim_printstate(elem, tl, hybrid,hvcoord,1, nelemd)

    is_model_inited = .true.

  end subroutine prim_init_model_f90 

  subroutine prim_run_f90 (dt) bind(c)
    use dimensions_mod,    only: nelemd
    use prim_driver_mod,   only: prim_run_subcycle
    use time_mod,          only: tstep
    use homme_context_mod, only: is_model_inited, elem, hybrid, tl, hvcoord, par
    !
    ! Input(s)
    !
    real (kind=c_double), intent(in) :: dt

    if (.not. is_model_inited) then
      call abortmp ("Error! prim_init_model_f90 was not called yet (or prim_finalize_f90 was already called).\n")
    endif

    ! Set dt in the time mod
    tstep = dt

    if (par%masterproc) print *, "HOMME step: ", tl%nstep

    call prim_run_subcycle(elem,hybrid,1,nelemd,dt,.false.,tl,hvcoord,1)

  end subroutine prim_run_f90

  subroutine prim_finalize_f90 () bind(c)
    use homme_context_mod,    only: is_model_inited, elem, dom_mt, par
    use prim_cxx_driver_base, only: prim_finalize

    if (.not. is_model_inited) then
      call abortmp ("Error! prim_init_model_f90 was not called yet (or prim_finalize_f90 was already called).\n")
    endif

    ! Cleanup some f90 stuff in Homme, and cleanup all cxx structures
    call prim_finalize()

    ! Deallocate the pointers
    deallocate (elem)
    deallocate (dom_mt)

    call t_stopf('Total')

    ! Finalize and write the timings
    if(par%masterproc) print *,"writing timing data"
    call t_prf('HommeTime', par%comm)

    if(par%masterproc) print *,"calling t_finalizef"
    call t_finalizef()

    is_model_inited = .false.

  end subroutine prim_finalize_f90

end module homme_driver_mod
