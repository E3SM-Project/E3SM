#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
module prim_driver_mod

  use prim_driver_base, only:&
    deriv1,&
    smooth_topo_datasets

  implicit none


  private

  public :: prim_init1
  public :: prim_init2
  public :: prim_run_subcycle
  public :: prim_finalize

  private :: generate_global_to_local
  private :: init_cxx_connectivity_internal

#include <mpif.h>

  contains

  subroutine prim_init1(elem, par, dom_mt, tl)
    use iso_c_binding,    only : c_int, c_loc
    use derivative_mod,   only : derivinit
    use dimensions_mod,   only : nelemd, np
    use domain_mod,       only : domain1d_t
    use element_mod,      only : element_t
    use kinds,            only : iulog, real_kind
    use parallel_mod,     only : parallel_t
    use time_mod,         only : TimeLevel_t, TimeLevel_init
    use prim_driver_base, only : prim_init1_geometry, prim_init1_elem_arrays, prim_init1_cleanup, &
                                 MetaVertex, GridEdge, deriv1
#ifndef CAM
    use prim_driver_base, only : prim_init1_no_cam
#endif

    interface
      subroutine reset_cxx_comm (f_comm) bind(c)
        use iso_c_binding, only: c_int
        !
        ! Inputs
        !
        integer(kind=c_int), intent(in) :: f_comm
      end subroutine reset_cxx_comm
      subroutine initialize_hommexx_session() bind(c)
      end subroutine initialize_hommexx_session

    end interface
    !
    ! Inputs
    !
    type (element_t),   pointer     :: elem(:)
    type (parallel_t),  intent(in)  :: par
    type (domain1d_t),  pointer     :: dom_mt(:)
    type (timelevel_t), intent(out) :: tl

    ! Initialize MPI comm in C++
    call reset_cxx_comm (INT(par%comm,c_int))

    ! Initialize kokkos before any environment changes from the Fortran
    call initialize_hommexx_session()

#ifndef CAM
    call prim_init1_no_cam(par)
#endif

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv1)

    ! ==================================
    ! Initialize and partition the geometry
    ! ==================================
    call prim_init1_geometry(elem,par,dom_mt)

    ! ==================================
    ! Initialize element pointers
    ! ==================================
    call setup_element_pointers(elem)

    ! ==================================
    ! Initialize element arrays (fluxes and state)
    ! ==================================
    call prim_init1_elem_arrays(elem,par)

    ! ==================================
    ! Initialize C++ mpi communication structures
    ! ==================================
    call init_cxx_connectivity(nelemd,GridEdge,MetaVertex,par)

    ! Initialize the time levels
    call TimeLevel_init(tl)

    ! Cleanup the tmp stuff used in prim_init1_geometry
    call prim_init1_cleanup()

    if(par%masterproc) write(iulog,*) 'end of prim_init1'
  end subroutine prim_init1

  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)
    use iso_c_binding,    only : c_loc, c_ptr, c_bool
    use control_mod,      only : limiter_option, rsplit, qsplit, tstep_type, statefreq,  &
                                 nu, nu_p, nu_q, nu_s, nu_div, nu_top, vert_remap_q_alg, &
                                 hypervis_order, hypervis_subcycle, hypervis_scaling,    &
                                 ftype, prescribed_wind, moisture, disable_diagnostics,  &
                                 use_cpstar, use_semi_lagrange_transport
    use dimensions_mod,   only : qsize, nelemd, np, qsize
    use element_mod,      only : element_t
    use element_state,    only : elem_state_v, elem_state_temp, elem_state_dp3d,       &
                                 elem_state_Q, elem_state_Qdp, elem_state_ps_v,        &
                                 elem_accum_qvar, elem_accum_qmass, elem_accum_q1mass, &
                                 elem_accum_kener, elem_accum_pener, elem_accum_iener, &
                                 elem_accum_iener_wet
    use hybrid_mod,       only : hybrid_t
    use hybvcoord_mod,    only : hvcoord_t
    use time_mod,         only : timelevel_t
    use kinds,            only : real_kind
    use prim_driver_base, only : deriv1, prim_init2_base => prim_init2
    use prim_state_mod,   only : prim_printstate

    interface
      subroutine init_derivative_c (dvv_ptr) bind(c)
        use iso_c_binding, only : c_ptr
        !
        ! Inputs
        !
        type (c_ptr), intent(in) :: dvv_ptr
      end subroutine init_derivative_c
      subroutine init_simulation_params_c (remap_alg, limiter_option, rsplit, qsplit, time_step_type,    &
                                           qsize, state_frequency, nu, nu_p, nu_q, nu_s, nu_div, nu_top, &
                                           hypervis_order, hypervis_subcycle, hypervis_scaling,          &
                                           ftype, prescribed_wind, moisture, disable_diagnostics,        &
                                           use_cpstar, use_semi_lagrange_transport) bind(c)
        use iso_c_binding, only: c_int, c_bool, c_double
        !
        ! Inputs
        !
        integer(kind=c_int),  intent(in) :: remap_alg, limiter_option, rsplit, qsplit, time_step_type
        integer(kind=c_int),  intent(in) :: state_frequency, qsize
        real(kind=c_double),  intent(in) :: nu, nu_p, nu_q, nu_s, nu_div, nu_top, hypervis_scaling
        integer(kind=c_int),  intent(in) :: hypervis_order, hypervis_subcycle
        integer(kind=c_int),  intent(in) :: ftype
        logical(kind=c_bool), intent(in) :: prescribed_wind, moisture, disable_diagnostics, use_cpstar, use_semi_lagrange_transport
      end subroutine init_simulation_params_c
      subroutine init_elements_c (nelemd) bind(c)
        use iso_c_binding, only : c_int
        !
        ! Inputs
        !
        integer (kind=c_int), intent(in) :: nelemd
      end subroutine init_elements_c
      subroutine init_elements_2d_c (ie, D_ptr, Dinv_ptr, elem_fcor_ptr,                  &
                                     elem_mp_ptr, elem_spheremp_ptr, elem_rspheremp_ptr,      &
                                     elem_metdet_ptr, elem_metinv_ptr, phis_ptr,              &
                                     tensorvisc_ptr, vec_sph2cart_ptr) bind(c)
        use iso_c_binding, only : c_ptr, c_int
        !
        ! Inputs
        !
        integer (kind=c_int), intent(in) :: ie
        type (c_ptr) , intent(in) :: D_ptr, Dinv_ptr, elem_fcor_ptr
        type (c_ptr) , intent(in) :: elem_mp_ptr, elem_spheremp_ptr, elem_rspheremp_ptr
        type (c_ptr) , intent(in) :: elem_metdet_ptr, elem_metinv_ptr, phis_ptr
        type (c_ptr) , intent(in) :: tensorvisc_ptr, vec_sph2cart_ptr
      end subroutine init_elements_2d_c
      subroutine init_diagnostics_c (elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr,           &
                                     elem_accum_q1mass_ptr, elem_accum_iener_ptr, elem_accum_iener_wet_ptr, &
                                     elem_accum_kener_ptr, elem_accum_pener_ptr) bind(c)
        use iso_c_binding, only : c_ptr
        !
        ! Inputs
        !
        type (c_ptr), intent(in) :: elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr
        type (c_ptr), intent(in) :: elem_accum_q1mass_ptr, elem_accum_iener_ptr, elem_accum_iener_wet_ptr
        type (c_ptr), intent(in) :: elem_accum_kener_ptr, elem_accum_pener_ptr
      end subroutine init_diagnostics_c
      subroutine init_elements_states_c (elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr,   &
                                         elem_state_Qdp_ptr, elem_state_ps_v_ptr) bind(c)
        use iso_c_binding, only : c_ptr
        !
        ! Inputs
        !
        type (c_ptr) :: elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr
        type (c_ptr) :: elem_state_Qdp_ptr, elem_state_ps_v_ptr
      end subroutine init_elements_states_c
      subroutine init_boundary_exchanges_c () bind(c)
      end subroutine init_boundary_exchanges_c
      subroutine init_hvcoord_c (ps0,hybrid_am_ptr,hybrid_ai_ptr,hybrid_bm_ptr,hybrid_bi_ptr) bind(c)
        use iso_c_binding , only : c_ptr, c_double
        !
        ! Inputs
        !
        real (kind=c_double),  intent(in) :: ps0
        type (c_ptr),          intent(in) :: hybrid_am_ptr, hybrid_ai_ptr
        type (c_ptr),          intent(in) :: hybrid_bm_ptr, hybrid_bi_ptr
      end subroutine init_hvcoord_c
      subroutine init_time_level_c(nm1,n0,np1,nstep,nstep0) bind(c)
        use iso_c_binding, only: c_int
        !
        ! Inputs
        !
        integer(kind=c_int), intent(in) :: nm1, n0, np1, nstep, nstep0
      end subroutine init_time_level_c
    end interface

    !
    ! Inputs
    !
    type (element_t),   intent(inout)         :: elem(:)
    type (hybrid_t),    intent(in)            :: hybrid
    type (TimeLevel_t), intent(inout)         :: tl       ! time level struct
    type (hvcoord_t),   intent(inout), target :: hvcoord  ! hybrid vertical coordinate struct
    integer,            intent(in)            :: nets     ! starting thread element number (private)
    integer,            intent(in)            :: nete     ! ending thread element number   (private)
    !
    ! Locals
    !
    integer :: ie
    real (kind=real_kind), target :: dvv (np,np)

    real (kind=real_kind), target, dimension(np,np,2,2)     :: elem_D, elem_Dinv, elem_metinv, elem_tensorvisc
    real (kind=real_kind), target, dimension(np,np)         :: elem_mp, elem_fcor, elem_spheremp
    real (kind=real_kind), target, dimension(np,np)         :: elem_rspheremp, elem_metdet, elem_state_phis
    real (kind=real_kind), target, dimension(np,np,3,2)     :: elem_vec_sph2cart

    type (c_ptr) :: hybrid_am_ptr, hybrid_ai_ptr, hybrid_bm_ptr, hybrid_bi_ptr
    type (c_ptr) :: elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr
    type (c_ptr) :: elem_mp_ptr, elem_spheremp_ptr, elem_rspheremp_ptr
    type (c_ptr) :: elem_metdet_ptr, elem_metinv_ptr, elem_state_phis_ptr
    type (c_ptr) :: elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr
    type (c_ptr) :: elem_state_q_ptr, elem_state_Qdp_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_tensorvisc_ptr, elem_vec_sph2cart_ptr
    type (c_ptr) :: elem_accum_qvar_ptr, elem_accum_qmass_ptr, elem_accum_q1mass_ptr
    type (c_ptr) :: elem_accum_iener_ptr, elem_accum_iener_wet_ptr
    type (c_ptr) :: elem_accum_kener_ptr, elem_accum_pener_ptr

    ! Call the base version of prim_init2
    call prim_init2_base(elem,hybrid,nets,nete,tl,hvcoord)

    ! Initialize the C++ derivative structure
    dvv = deriv1%dvv
    call init_derivative_c(c_loc(dvv))

    ! Fill the simulation params structures in C++
    call init_simulation_params_c (vert_remap_q_alg, limiter_option, rsplit, qsplit, tstep_type,  &
                                   qsize, statefreq, nu, nu_p, nu_q, nu_s, nu_div, nu_top,        &
                                   hypervis_order, hypervis_subcycle, hypervis_scaling,           &
                                   ftype, LOGICAL(prescribed_wind==1,c_bool),                     &
                                   LOGICAL(moisture/="dry",c_bool),                               &
                                   LOGICAL(disable_diagnostics,c_bool),                           &
                                   LOGICAL(use_cpstar==1,c_bool),                           &
                                   LOGICAL(use_semi_lagrange_transport,c_bool))

    ! Initialize the C++ elements structure
    call init_elements_c (nelemd)

    ! Initialize the 2d element arrays in C++
    elem_D_ptr            = c_loc(elem_D)
    elem_Dinv_ptr         = c_loc(elem_Dinv)
    elem_fcor_ptr         = c_loc(elem_fcor)
    elem_mp_ptr           = c_loc(elem_mp)
    elem_spheremp_ptr     = c_loc(elem_spheremp)
    elem_rspheremp_ptr    = c_loc(elem_rspheremp)
    elem_metdet_ptr       = c_loc(elem_metdet)
    elem_metinv_ptr       = c_loc(elem_metinv)
    elem_tensorvisc_ptr   = c_loc(elem_tensorvisc)
    elem_vec_sph2cart_ptr = c_loc(elem_vec_sph2cart)
    elem_state_phis_ptr   = c_loc(elem_state_phis)

    do ie=1,nelemd
      elem_D            = elem(ie)%D
      elem_Dinv         = elem(ie)%Dinv
      elem_fcor         = elem(ie)%fcor
      elem_mp           = elem(ie)%mp
      elem_spheremp     = elem(ie)%spheremp
      elem_rspheremp    = elem(ie)%rspheremp
      elem_metdet       = elem(ie)%metdet
      elem_metinv       = elem(ie)%metinv
      elem_state_phis   = elem(ie)%state%phis
      elem_tensorvisc   = elem(ie)%tensorVisc
      elem_vec_sph2cart = elem(ie)%vec_sphere2cart
      call init_elements_2d_c (ie-1, elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr,        &
                               elem_mp_ptr, elem_spheremp_ptr, elem_rspheremp_ptr,    &
                               elem_metdet_ptr, elem_metinv_ptr, elem_state_phis_ptr, &
                               elem_tensorvisc_ptr, elem_vec_sph2cart_ptr)
    enddo

    ! Initialize the 3d element arrays in C++
    elem_state_v_ptr    = c_loc(elem_state_v)
    elem_state_temp_ptr = c_loc(elem_state_temp)
    elem_state_dp3d_ptr = c_loc(elem_state_dp3d)
    elem_state_Qdp_ptr  = c_loc(elem_state_Qdp)
    elem_state_ps_v_ptr = c_loc(elem_state_ps_v)
    call init_elements_states_c (elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr,   &
                                 elem_state_Qdp_ptr, elem_state_ps_v_ptr)

    ! Initialize the diagnostics arrays in C++
    elem_state_q_ptr         = c_loc(elem_state_q)
    elem_accum_qvar_ptr      = c_loc(elem_accum_qvar)
    elem_accum_qmass_ptr     = c_loc(elem_accum_qmass)
    elem_accum_q1mass_ptr    = c_loc(elem_accum_q1mass)
    elem_accum_iener_ptr     = c_loc(elem_accum_iener)
    elem_accum_iener_wet_ptr = c_loc(elem_accum_iener_wet)
    elem_accum_kener_ptr     = c_loc(elem_accum_kener)
    elem_accum_pener_ptr     = c_loc(elem_accum_pener)
    call init_diagnostics_c (elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr,           &
                             elem_accum_q1mass_ptr, elem_accum_iener_ptr, elem_accum_iener_wet_ptr, &
                             elem_accum_kener_ptr, elem_accum_pener_ptr)

    ! Initialize the hybrid vertical coordinate in C++
    hybrid_am_ptr = c_loc(hvcoord%hyam)
    hybrid_ai_ptr = c_loc(hvcoord%hyai)
    hybrid_bm_ptr = c_loc(hvcoord%hybm)
    hybrid_bi_ptr = c_loc(hvcoord%hybi)
    call init_hvcoord_c (hvcoord%ps0,hybrid_am_ptr,hybrid_ai_ptr,hybrid_bm_ptr,hybrid_bi_ptr)

    ! Initialize boundary exchange structure in C++
    call init_boundary_exchanges_c ()

    ! Initialize time level structure in C++
    call init_time_level_c(tl%nm1, tl%n0, tl%np1, tl%nstep, tl%nstep0)

  end subroutine prim_init2

  subroutine prim_run_subcycle(elem, hybrid, nets, nete, dt, single_column, tl, hvcoord,nsubstep)
    use iso_c_binding,  only : c_int, c_ptr, c_loc
    use control_mod,    only : qsplit, rsplit, statefreq
    use dimensions_mod, only : nelemd
    use element_mod,    only : element_t
    use element_state,  only : elem_state_v, elem_state_temp, elem_state_dp3d, &
                               elem_state_Qdp, elem_state_Q, elem_state_ps_v, elem_derived_omega_p
    use hybrid_mod,     only : hybrid_t
    use hybvcoord_mod,  only : hvcoord_t
    use kinds,          only : real_kind
    use time_mod,       only : timelevel_t, nextOutputStep, nsplit
    use control_mod,    only : statefreq
    use parallel_mod,   only : abortmp
    use prim_state_mod, only: prim_printstate
    interface
      subroutine prim_run_subcycle_c(tstep,nstep,nm1,n0,np1,next_output_step) bind(c)
        use iso_c_binding, only: c_int, c_double
        !
        ! Inputs
        !
        integer(kind=c_int),  intent(in) :: nstep, nm1, n0, np1, next_output_step
        real (kind=c_double), intent(in) :: tstep
      end subroutine prim_run_subcycle_c


      subroutine cxx_push_results_to_f90(elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr, &
                                         elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr,  &
                                         elem_derived_omega_p_ptr) bind(c)
        use iso_c_binding , only : c_ptr
        !
        ! Inputs
        !
        type (c_ptr),          intent(in) :: elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr
        type (c_ptr),          intent(in) :: elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
        type (c_ptr),          intent(in) :: elem_derived_omega_p_ptr
      end subroutine cxx_push_results_to_f90
    end interface
    !
    ! Inputs
    !
    type (element_t) ,    intent(inout) :: elem(:)
    type (hybrid_t),      intent(in)    :: hybrid                       ! distributed parallel structure (shared)
    type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets                         ! starting thread element number (private)
    integer,              intent(in)    :: nete                         ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt                           ! "timestep dependent" timestep
    logical,              intent(in)    :: single_column
    type (TimeLevel_t),   intent(inout) :: tl
    integer,              intent(in)    :: nsubstep                     ! nsubstep = 1 .. nsplit
    !
    ! Locals
    !
    logical :: compute_diagnostics
    integer (kind=c_int) :: nstep_end, nstep_c, nm1_c, n0_c, np1_c
    type (c_ptr) :: elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr
    type (c_ptr) :: elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_derived_omega_p_ptr

    if (nets/=1 .or. nete/=nelemd) then
      call abortmp ('We don''t allow to call C routines from a horizontally threaded region')
    endif

    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif
    compute_diagnostics   = .false.
    if (MODULO(nstep_end,statefreq)==0 .or. (tl%nstep <= tl%nstep0+(nstep_end-tl%nstep) )) then
      compute_diagnostics = .true.
    endif

    call prim_run_subcycle_c(dt,nstep_c,nm1_c,n0_c,np1_c,nextOutputStep)

    ! Set final timelevels from C into Fortran structure
    tl%nstep = nstep_c
    tl%nm1   = nm1_c + 1
    tl%n0    = n0_c  + 1
    tl%np1   = np1_c + 1

    if (MODULO(tl%nstep,statefreq)==0 .or. tl%nstep >= nextOutputStep) then
      ! Set pointers to states
      elem_state_v_ptr         = c_loc(elem_state_v)
      elem_state_temp_ptr      = c_loc(elem_state_temp)
      elem_state_dp3d_ptr      = c_loc(elem_state_dp3d)
      elem_state_Qdp_ptr       = c_loc(elem_state_Qdp)
      elem_state_Q_ptr         = c_loc(elem_state_Q)
      elem_state_ps_v_ptr      = c_loc(elem_state_ps_v)
      elem_derived_omega_p_ptr = c_loc(elem_derived_omega_p)

      ! Copy cxx arrays back to f90 structures
      call cxx_push_results_to_f90(elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr, &
                                   elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr, &
                                   elem_derived_omega_p_ptr)
    endif

    ! Print some diagnostic information

    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if

  end subroutine prim_run_subcycle

  subroutine prim_finalize ()
    interface
      subroutine finalize_hommexx_session() bind(c)
      end subroutine finalize_hommexx_session
    end interface

    ! This call lets C++ destroy the singleton containing all the views,
    ! and finalize the Kokkos execution space.
    call finalize_hommexx_session()
  end subroutine prim_finalize

!!!!!!!!!!!!!!!!!!!!!!! PRIVATE SUBROUTINES BELOW THIS LINE !!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_cxx_connectivity (nelemd, GridEdge, MetaVertex, par)
    use dimensions_mod, only : nelem
    use gridgraph_mod,  only : GridEdge_t
    use metagraph_mod,  only : MetaVertex_t
    use parallel_mod,   only : parallel_t
    !
    ! Inputs
    !
    integer, intent(in) :: nelemd
    type(GridEdge_t),   intent(in) :: GridEdge(:)
    type (parallel_t),  intent(in) :: par
    type(MetaVertex_t), intent(in) :: MetaVertex
    !
    ! Locals
    !
    integer :: Global2Local(nelem)

    call generate_global_to_local(MetaVertex,Global2Local,par)

    call init_cxx_connectivity_internal (nelemd, Global2Local, Gridedge)

  end subroutine init_cxx_connectivity

  subroutine generate_global_to_local (MetaVertex, Global2Local, par)
    use dimensions_mod, only : nelem
    use metagraph_mod,  only : MetaVertex_t
    use parallel_mod,   only : parallel_t, MPI_MAX, MPIinteger_t
    !
    ! Inputs
    !
    type (parallel_t),  intent(in) :: par
    type(MetaVertex_t), intent(in) :: MetaVertex
    integer, intent(out) :: Global2Local(nelem)
    !
    ! Locals
    !
    integer :: ie, ierr

    ! Defaults all local ids to 0 (meaning not on this process)
    Global2Local = 0
    do ie=1,SIZE(MetaVertex%members)
      Global2Local(MetaVertex%members(ie)%number) = ie
    enddo

    call MPI_Allreduce(MPI_IN_PLACE,Global2Local,nelem,MPIinteger_t,MPI_MAX,par%comm,ierr)

  end subroutine generate_global_to_local

  subroutine init_cxx_connectivity_internal (nelemd,Global2Local,GridEdge)
    use gridgraph_mod,  only : GridEdge_t
    use dimensions_mod, only : nelem
    !
    ! Interfaces
    !
    interface
      subroutine init_connectivity (num_local_elems) bind (c)
        use iso_c_binding, only : c_int
        !
        ! Inputs
        !
        integer (kind=c_int), intent(in) :: num_local_elems
      end subroutine init_connectivity

      subroutine finalize_connectivity () bind(c)
      end subroutine finalize_connectivity

      subroutine add_connection (first_lid,  first_gid,  first_pos,  first_pid, &
                                 second_lid, second_gid, second_pos, second_pid) bind(c)
        use iso_c_binding, only : c_int
        !
        ! Inputs
        !
        integer (kind=c_int), intent(in) :: first_lid,  first_gid,  first_pos,  first_pid
        integer (kind=c_int), intent(in) :: second_lid, second_gid, second_pos, second_pid
      end subroutine add_connection
    end interface
    !
    ! Inputs
    !
    integer, dimension(nelem), intent(in) :: Global2Local
    integer, intent(in) :: nelemd
    type(GridEdge_t), intent(in) :: GridEdge(:)
    !
    ! Locals
    !
    integer :: ie, num_edges
    type(GridEdge_t) :: e

    call init_connectivity(nelemd)

    num_edges = SIZE(GridEdge)
    do ie=1,num_edges
      e = GridEdge(ie)
      call add_connection(Global2Local(e%head%number),e%head%number,e%head_dir,e%head%processor_number, &
                          Global2Local(e%tail%number),e%tail%number,e%tail_dir,e%tail%processor_number)
    enddo

    call finalize_connectivity()
  end subroutine init_cxx_connectivity_internal

  subroutine setup_element_pointers (elem)
    use element_mod,    only : element_t
    use element_state,  only : allocate_element_arrays, elem_state_v, elem_state_temp, elem_state_dp3d,    &
                               elem_state_Qdp, elem_state_Q, elem_state_ps_v, elem_derived_omega_p,        &
                               elem_accum_pener, elem_accum_kener, elem_accum_iener, elem_accum_iener_wet, &
                               elem_accum_qvar, elem_accum_qmass, elem_accum_q1mass, elem_state_phis
    use dimensions_mod, only : nelemd
    !
    ! Inputs
    !
    type (element_t), intent(inout) :: elem(:)
    !
    ! Locals
    !
    integer :: ie

    call allocate_element_arrays(nelemd)

    do ie=1,nelemd
      elem(ie)%state%v         => elem_state_v(:,:,:,:,:,ie)
      elem(ie)%state%T         => elem_state_temp(:,:,:,:,ie)
      elem(ie)%state%dp3d      => elem_state_dp3d(:,:,:,:,ie)
      elem(ie)%state%ps_v      => elem_state_ps_v(:,:,:,ie)
      elem(ie)%state%Q         => elem_state_Q(:,:,:,:,ie)
      elem(ie)%state%Qdp       => elem_state_Qdp(:,:,:,:,:,ie)
      elem(ie)%state%phis      => elem_state_phis(:,:,ie)
      elem(ie)%derived%omega_p => elem_derived_omega_p(:,:,:,ie)

      elem(ie)%accum%KEner     => elem_accum_KEner    (:,:,:,ie)
      elem(ie)%accum%PEner     => elem_accum_PEner    (:,:,:,ie)
      elem(ie)%accum%IEner     => elem_accum_IEner    (:,:,:,ie)
      elem(ie)%accum%IEner_wet => elem_accum_IEner_wet(:,:,:,ie)
      elem(ie)%accum%Qvar      => elem_accum_Qvar     (:,:,:,:,ie)
      elem(ie)%accum%Qmass     => elem_accum_Qmass    (:,:,:,:,ie)
      elem(ie)%accum%Q1mass    => elem_accum_Q1mass   (:,:,:,ie)

    enddo

  end subroutine setup_element_pointers
end module

