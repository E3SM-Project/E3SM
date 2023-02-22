#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! This module contains constants and namelist variables used through out the model
! to avoid circular dependancies please do not 'use' any further modules here.
!
module control_mod
  use kinds, only : real_kind
  use physical_constants, only: dd_pi

  implicit none

  integer, public, parameter :: MAX_STRING_LEN=240
  integer, public, parameter :: MAX_FILE_LEN=240
  character(len=MAX_STRING_LEN)    , public :: integration    ! time integration (explicit, or full imp)

  ! Tracer transport algorithm type:
  !     0  spectral-element Eulerian
  !    12 interpolation semi-Lagrangian
  integer, public  :: transport_alg = 0
  ! Constrained density reconstructor for SL property preservation; not used if
  ! transport_alg = 0:
  !     0  none
  !     2  QLT
  !     3  CAAS
  !    20  QLT  with superlevels
  !    30  CAAS with superlevels
  integer, public  :: semi_lagrange_cdr_alg = 3
  ! If true, check mass conservation and shape preservation. The second
  ! implicitly checks tracer consistency.
  logical, public  :: semi_lagrange_cdr_check = .false.
  ! If > 0 and nu_q > 0, apply hyperviscosity to tracers 1 through this value,
  ! rather than just those that couple to the dynamics at the dynamical time
  ! step. These latter are 'active' tracers, in contrast to 'passive' tracers
  ! that directly couple only to the physics.
  integer, public  :: semi_lagrange_hv_q = 1
  ! If >= 1, then the SL algorithm may choose a nearby point inside the element
  ! halo available to it if the actual point is outside the halo. This is done
  ! in levels <= this parameter.
  integer, public :: semi_lagrange_nearest_point_lev = 256

! flag used by preqx, theta-l and theta-c models
! should be renamed to "hydrostatic_mode"
  logical, public :: theta_hydrostatic_mode


  integer, public  :: tstep_type= 5                           ! preqx timestepping options
  integer, public  :: rk_stage_user  = 0                      ! number of RK stages (shallow water model) 
  integer, public  :: ftype = 0                                ! Forcing Type
                                                               ! ftype = 0  HOMME ApplyColumn() type forcing process split
                                                               ! ftype = -1   ignore forcing  (used for testing energy balance)
                                              
  integer, public :: qsplit = 1           ! ratio of dynamics tsteps to tracer tsteps
  integer, public :: rsplit = 0           ! for vertically lagrangian dynamics, apply remap
                                          ! every rsplit tracer timesteps

  ! These factors replace rsplit and qsplit.
  !   If dt_remap_factor = 0, use vertically Eulerian dynamics.
  !   If dt_remap_factor > 0, the vertical remap time step is
  ! dt_remap_factor*tstep.
  !   The tracer transport time step is dt_tracer_factor*tstep.
  !   The smaller of dt_remap_factor and dt_tracer_factor must divide
  ! the larger.
  !   If dt_remap_factor >= dt_tracer_factor, then
  !     new dt_tracer_factor == old qsplit
  !     new dt_remap_factor == old dt_tracer_factor/dt_remap_factor
  ! Default values make qsplit and rsplit control the time steps.
  integer, public :: dt_remap_factor = -1, dt_tracer_factor = -1

  integer, public :: prim_step_type = -1 ! 1 = old code for EUL, 2 = prim_run_flexible for SL
                                          ! -1 means it wasn't set, error

  integer, public :: LFTfreq=0            ! leapfrog-trapazoidal frequency (shallow water only)
                                          ! interspace a lf-trapazoidal step every LFTfreq leapfrogs    
                                          ! 0 = disabled

! vert_remap_q_alg:   -1  PPM remap without monotone filter, used for some test cases
!                      0  Zerroukat monotonic splines
!                      1  PPM vertical remap with constant extension at the boundaries
!                     10  PPM with linear extrapolation at boundaries, with column limiter
!                     11  PPM with unlimited linear extrapolation at boundaries
 integer, public :: vert_remap_q_alg = 0    ! tracers
 integer, public :: vert_remap_u_alg = -2   ! remap for dynamics. default -2 means inherit vert_remap_q_alg

! advect theta 0: conservation form 
!              1: expanded divergence form (less noisy, non-conservative)
 integer, public :: theta_advect_form = 0
 real (kind=real_kind), public :: vtheta_thresh = 100.d0  ! threshold for virtual potential temperature minimum limiter
 real (kind=real_kind), public :: dp3d_thresh   = 0.125d0 ! threshold for dp3d minimum limiter

 integer, public :: pgrad_correction  = 0   ! 1=turn on theta model pressure gradient correction
 integer, public :: hv_ref_profiles   = 0   ! 1=turn on theta model HV reference profiles
 integer, public :: hv_theta_correction=0   ! 1=use HV on p-surface approximation for theta
 real (kind=real_kind), public :: hv_theta_thresh=.025d0  ! d(theta)/dp max threshold for HV correction term

 integer, public :: cubed_sphere_map = -1  ! -1 = chosen at run time
                                           !  0 = equi-angle Gnomonic (default)
                                           !  1 = equi-spaced Gnomonic (not yet coded)
                                           !  2 = element-local projection  (for var-res)
                                           !  3 = parametric (not yet coded)

!tolerance to define smth small, was introduced for lim 8 in 2d and 3d
  real (kind=real_kind), public, parameter :: tol_limiter=1e-13

  integer              , public :: limiter_option = 0
  character(len=MAX_STRING_LEN)    , public :: precon_method  ! if semi_implicit, type of preconditioner:
                                                  ! choices block_jacobi or identity

  integer              , public :: coord_transform_method   ! If zoltan2 is used, various ways of representing the coordinates methods
                                                            ! Instead of using the sphere coordinates, it is better to use cube or projected 2D coordinates for quality.
                                                            ! check params_mod for options.

  integer              , public :: z2_map_method            ! If zoltan2 is used,
                                                            ! Task mapping method to be used by zoltan2.
                                                            ! Z2_NO_TASK_MAPPING        (1) - is no task mapping
                                                            ! Z2_TASK_MAPPING           (2) - performs default task mapping of zoltan2.
                                                            ! Z2_OPTIMIZED_TASK_MAPPING (3) - includes network aware optimizations.
                                                            ! Use (3) if zoltan2 is enabled.

  integer              , public :: partmethod     ! partition methods
  character(len=MAX_STRING_LEN)    , public :: topology = "cube"       ! options: "cube", "plane"
  character(len=MAX_STRING_LEN)    , public :: geometry = "sphere"      ! options: "sphere", "plane"
  character(len=MAX_STRING_LEN)    , public :: test_case
  !most tests don't have forcing
  logical                          , public :: test_with_forcing = .false. 
  integer              , public :: statefreq      ! output frequency of synopsis of system state (steps)
  integer              , public :: restartfreq
  integer              , public :: runtype 
  integer              , public :: timerdetail 
  integer              , public :: numnodes 
  character(len=MAX_STRING_LEN)    , public :: restartfile 
  character(len=MAX_STRING_LEN)    , public :: restartdir

  ! flag used for "slice" planar tests (no variation in y-dir)
  logical, public :: planar_slice
  
! namelist variable set to dry,notdry,moist
! internally the code should use logical "use_moisture"
  character(len=MAX_STRING_LEN)    , public :: moisture  

  integer, public  :: use_cpstar=0          ! use cp or cp* in thermodynamics
  logical, public  :: use_moisture=.false.  ! use Q(:,:,:,1) to compute T_v

  
  integer              , public :: maxits         ! max iterations of solver
  real (kind=real_kind), public :: tol            ! solver tolerance (convergence criteria)
  integer              , public :: debug_level    ! debug level of CG solver


  character(len=MAX_STRING_LEN)    ,public  :: vfile_int=""   ! vertical formulation (ecmwf,ccm1)
  character(len=MAX_STRING_LEN)    ,public  :: vfile_mid=""   ! vertical grid spacing (equal,unequal)
  integer,                          public  :: vanalytic = 0  ! if 1, test initializes vertical coords
  real (kind=real_kind),            public  :: vtop = 0.1     ! top coordinate level for analytic vcoords

  real (kind=real_kind), public :: nu      = 7.0D5            ! viscosity (momentum equ)
  real (kind=real_kind), public :: nu_div  = -1               ! viscsoity (momentum equ, div component)
  real (kind=real_kind), public :: nu_s    = -1               ! default = nu   T equ. viscosity
  real (kind=real_kind), public :: nu_q    = -1               ! default = nu   tracer viscosity
  real (kind=real_kind), public :: nu_p    = -1               ! default = nu   ps equ. viscosity
  real (kind=real_kind), public :: nu_top  = 0.0D5            ! top-of-the-model viscosity
  real (kind=real_kind), public :: tom_sponge_start=0         ! start of sponge layer, in hPa

  integer, public :: hypervis_subcycle=1                      ! number of subcycles for hyper viscsosity timestep
  integer, public :: hypervis_subcycle_tom=0                  ! number of subcycles for TOM diffusion
                                                              !   0   apply together with hyperviscosity
                                                              !   >1  apply timesplit from hyperviscosity
  integer, public :: hypervis_subcycle_q=1                    ! number of subcycles for hyper viscsosity timestep on TRACERS
  integer, public :: hypervis_order=0                         ! laplace**hypervis_order.  0=not used  1=regular viscosity, 2=grad**4

  real (kind=real_kind), public :: hypervis_scaling=0         ! use tensor hyperviscosity

  !three types of hyper viscosity are supported right now:
  ! (1) const hv:    nu * del^2 del^2
  ! (2) tensor hv,   nu * ( \div * tensor * \grad ) * del^2
  !
  ! (1) hypervis_scaling=0
  ! (2) tensor HV var-res grids  
  !            tensor within each element:
  !            set hypervis_scaling > 0 (typical values would be 3.0)
  !            (\div * tensor * \grad) operator uses cartesian laplace
  !

  ! hyperviscosity parameters used for smoothing topography
  integer, public :: smooth_phis_numcycle = -1   ! -1 = disable
  integer, public :: smooth_phis_p2filt = -1     ! -1 = disable
  real (kind=real_kind), public :: smooth_phis_nudt = 0

  integer, public :: prescribed_wind=0    ! fix the velocities?


  integer, public, parameter :: west  = 1
  integer, public, parameter :: east  = 2
  integer, public, parameter :: south = 3
  integer, public, parameter :: north = 4
  integer, public, parameter :: swest = 5
  integer, public, parameter :: seast = 6
  integer, public, parameter :: nwest = 7
  integer, public, parameter :: neast = 8
  
  logical, public :: disable_diagnostics  = .FALSE.

  ! Physgrid parameters
  integer, public :: se_fv_phys_remap_alg = 1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test case global parameters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! generic test case parameter - can be used by any test case to define options
  integer, public :: sub_case = 1                  

  real (kind=real_kind), public :: initial_total_mass = 0    ! initial perturbation in JW test case
  real (kind=real_kind), public :: u_perturb   = 0         ! initial perturbation in JW test case
#ifndef CAM
  real (kind=real_kind), public :: pertlim = 0          !pertibation to temperature [like CESM]
#endif

  ! shallow water advection test paramters
  ! kmass = level index with density.  other levels contain test tracers
  integer, public  :: kmass  = -1
  integer, public  :: toy_chemistry = 0            !  1 = toy chemestry is turned on in 2D advection code
  real (kind=real_kind), public :: g_sw_output            	   = 9.80616D0          ! m s^-2

  ! parameters for dcmip12 test 2-0: steady state atmosphere with orography
  real(real_kind), public :: dcmip2_0_h0      = 2000.d0        ! height of mountain range        (meters)
  real(real_kind), public :: dcmip2_0_Rm      = 3.d0*dd_pi/4.d0   ! radius of mountain range        (radians)
  real(real_kind), public :: dcmip2_0_zetam   = dd_pi/16.d0       ! mountain oscillation half-width (radians)

  ! parameters for dcmip12 test 2-x: mountain waves
  real(real_kind), public :: dcmip2_x_ueq     = 20.d0          ! wind speed at equator (m/s)
  real(real_kind), public :: dcmip2_x_h0      = 250.0d0        ! mountain height       (m)
  real(real_kind), public :: dcmip2_x_d       = 5000.0d0       ! mountain half width   (m)
  real(real_kind), public :: dcmip2_x_xi      = 4000.0d0       ! mountain wavelength   (m)

  ! for dcmip 2014 test 4:
  integer,         public :: dcmip4_moist     = 1
  real(real_kind), public :: dcmip4_X         = 1.0d0 

  ! for dcmip 2016 test 2
  integer, public :: dcmip16_prec_type = 0;
  integer, public :: dcmip16_pbl_type  = 0;

  ! for dcmip 2016 test 3
  real (kind=real_kind), public :: dcmip16_mu      = 0        ! additional uniform viscosity (momentum)
  real (kind=real_kind), public :: dcmip16_mu_s    = 0        ! additional uniform viscosity (scalar dynamical variables)
  real (kind=real_kind), public :: dcmip16_mu_q    = -1       ! additional uniform viscosity (scalar tracers); -1 implies it defaults to dcmip16_mu_s value
  real (kind=real_kind), public :: interp_lon0     = 0.0d0

!PLANAR
  real (kind=real_kind), private, parameter :: tol_zero=1e-10 !tolerance to determine if lx,ly,sx,sy are set

  real (kind=real_kind), public :: bubble_T0 = 270.0       !bubble ref state
  real (kind=real_kind), public :: bubble_dT = 0.5         !bubble dTheta
  real (kind=real_kind), public :: bubble_xycenter = 0.0   !bubble xy position
  real (kind=real_kind), public :: bubble_zcenter = 3000.0 !bubble z position
  real (kind=real_kind), public :: bubble_ztop = 10000.0   !bubble z top
  real (kind=real_kind), public :: bubble_xyradius = 2000.0!bubble radius along x or y axis
  real (kind=real_kind), public :: bubble_zradius = 1500.0 !bubble radius along z axis
  logical,               public :: bubble_cosine  = .TRUE. !bubble uniform or cosine
  logical,               public :: bubble_moist  = .FALSE.    ! 
  real (kind=real_kind), public :: bubble_moist_drh = 0.0     !bubble dRH parameter
  real (kind=real_kind), public :: bubble_rh_background = 0.0 !bubble RH parameter
  integer,               public :: bubble_prec_type = 0       !0 kessler, 1 rj
  logical,               protected :: case_planar_bubble = .FALSE.

  public :: set_planar_defaults

contains

  function timestep_make_parameters_consistent(par, rsplit, qsplit, &
       dt_remap_factor, dt_tracer_factor, tstep, dtime, nsplit, nstep_factor, &
       abrtf, silent) result(status)

    ! Current and future development require a more flexibility in
    ! specifying time steps. This routine analyzes the settings and
    ! either sets unset ones consistently or provides an error message
    ! and aborts.
    !   A return value of 0 means success; <0 means there was an
    ! error. In the case of error, a message is written to iulog.
    !   If you want a value to be computed, set it to <0 on input.

    use parallel_mod, only: abortmp, parallel_t
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par
    integer, intent(inout) :: &
         ! Old method of specifying subcycles, in which the vertical
         ! remap time step is restricted to be at least as large as
         ! the tracer time step.
         rsplit, qsplit, &
         ! New method, which permits either time step to be the
         ! larger, subject to that one must divide the other.
         dt_remap_factor, dt_tracer_factor, &
         nsplit
    integer, intent(out) :: &
         ! On output, dtime/tstep.
         nstep_factor
    real(kind=real_kind), intent(inout) :: &
         ! Dynamics time step.
         tstep
    integer, intent(inout) :: &
         ! Physics-dynamics coupling time step.
         dtime
    logical, intent(in), optional :: abrtf, silent
    integer :: status

    real(kind=real_kind), parameter :: &
         zero = 0.0_real_kind, &
         eps = epsilon(1.0_real_kind), &
         divisible_tol = 1e3_real_kind*eps

    logical :: abort_in, silent_in

    abort_in = .true.
    if (present(abrtf)) abort_in = abrtf
    silent_in = .false.
    if (present(silent)) silent_in = silent

    status = timestep_make_subcycle_parameters_consistent( &
         par, rsplit, qsplit, dt_remap_factor, dt_tracer_factor, &
         abort_in, silent_in)
    if (status /= 0) return
    status = timestep_make_eam_parameters_consistent( &
         par, dt_remap_factor, dt_tracer_factor, nsplit, nstep_factor, tstep, dtime, &
         abort_in, silent_in)
  end function timestep_make_parameters_consistent

  function timestep_make_subcycle_parameters_consistent(par, rsplit, qsplit, &
       dt_remap_factor, dt_tracer_factor, abrtf, silent) result(status)

    use parallel_mod, only: abortmp, parallel_t
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par
    integer, intent(inout) :: rsplit, qsplit, dt_remap_factor, dt_tracer_factor
    logical, intent(in), optional :: abrtf, silent
    integer :: status

    integer :: qsplit_prev, rsplit_prev
    logical :: split_specified, factor_specified, split_is_master, abort_in, silent_in

    abort_in = .true.
    if (present(abrtf)) abort_in = abrtf
    silent_in = .false.
    if (present(silent)) silent_in = silent

    status = -1 ! error value for early returns on error

    split_specified = rsplit >= 0 .and. qsplit >= 1
    factor_specified = dt_remap_factor >= 0 .and. dt_tracer_factor >= 1

    if (.not. split_specified .and. .not. factor_specified) then
       if (par%masterproc .and. .not. silent_in) then
          write(iulog,*) 'Neither rsplit,qsplit nor dt_remap_factor,dt_tracer_factor &
               &are specified; one set must be.'
       end if
       if (abort_in) call abortmp('timestep_make_parameters_consistent: input error')
       return
    end if

    !! Process rsplit, qsplit, dt_remap_factor, dt_tracer_factor.

    ! To support namelists with defaulted qsplit, rsplit values, we
    ! permit (split_specified .and. factor_specified). In this case,
    ! factor_specified means factor values are used.

    split_is_master = .not. factor_specified

    if (split_is_master) then
       dt_remap_factor = rsplit*qsplit
       dt_tracer_factor = qsplit
    else
       if (dt_remap_factor > 0) then
          if (.not. (modulo(dt_remap_factor, dt_tracer_factor) == 0 .or. &
                     modulo(dt_tracer_factor, dt_remap_factor) == 0)) then
             if (par%masterproc .and. .not. silent_in) then
                write(iulog,*) 'dt_remap_factor and dt_tracer_factor were specified, &
                     &but neither divides the other.'
             end if
             if (abort_in) call abortmp('timestep_make_parameters_consistent: divisibility error')
             return
          end if
       end if
       qsplit_prev = qsplit
       rsplit_prev = rsplit
       qsplit = dt_tracer_factor
       ! This is the only inconsistent setting. But I want to keep
       ! this here b/c almost all uses of rsplit is simply for whether
       ! it's == 0, and I don't want to touch all those lines in this
       ! PR.
       if (dt_tracer_factor <= dt_remap_factor) then
          rsplit = dt_remap_factor/dt_tracer_factor
       else
          ! If rsplit cannot be set consistently (because
          ! dt_tracer_factor < dt_remap_factor), then just preserve
          ! the sign to distinguish between vertically Eulerian and
          ! Lagrangian methods.
          if (dt_remap_factor > 0) then
             rsplit = 1
          else
             rsplit = 0
          end if
       end if
       if (split_specified .and. (qsplit /= qsplit_prev .or. rsplit /= rsplit_prev) .and. &
            par%masterproc .and. .not. silent_in) then
          write(iulog,'(a,i2,a,i2,a,i2,a,i2,a)') &
               'dt_remap_factor and dt_tracer_factor were specified; changing qsplit from ', &
               qsplit_prev, ' to ', qsplit, ' and rsplit from ', rsplit_prev, ' to ', rsplit, '.'
       end if
    end if

    status = 0 ! success value
  end function timestep_make_subcycle_parameters_consistent

  function timestep_make_eam_parameters_consistent(par, dt_remap_factor, dt_tracer_factor, &
       nsplit, nstep_factor, tstep, dtime, abrtf, silent) result(status)

    use parallel_mod, only: abortmp, parallel_t
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par
    integer, intent(in) :: dt_remap_factor, dt_tracer_factor
    integer, intent(inout) :: nsplit
    integer, intent(out) :: nstep_factor
    real(kind=real_kind), intent(inout) :: tstep
    integer, intent(inout) :: dtime
    logical, intent(in), optional :: abrtf, silent
    integer :: status

    real(kind=real_kind), parameter :: &
         zero = 0.0_real_kind, &
         eps = epsilon(1.0_real_kind), &
         divisible_tol = 1e3_real_kind*eps

    real(kind=real_kind) :: nsplit_real, tmp
    integer :: dt_max_factor
    logical :: abort_in, silent_in

    abort_in = .true.
    if (present(abrtf)) abort_in = abrtf
    silent_in = .false.
    if (present(silent)) silent_in = silent

    status = -1 ! error value for early returns on error

    !! Process dtime, tstep, nsplit.
    dt_max_factor = max(dt_remap_factor, dt_tracer_factor)

    ! Every 'if' has an 'else', so every case is covered.
    if (nsplit > 0) nstep_factor = dt_max_factor*nsplit
    if (dtime > 0) then
       if (nsplit > zero) then
          tmp = real(dtime, real_kind)/real(nstep_factor, real_kind)
          if (tstep > zero) then
             if (abs(tstep - tmp) > divisible_tol*tmp) then
                if (par%masterproc .and. .not. silent_in) then
                   write(iulog,'(a,a,i6,a,i2,a,es11.4,a,i2)') &
                        'dtime, nsplit, tstep were all >0 on input, but they disagree: ', &
                        'dtime ', dtime, ' nsplit ', nsplit, ' tstep ', tstep, ' nstep_factor ', nstep_factor
                end if
                if (abort_in) call abortmp('timestep_make_parameters_consistent: divisibility error')
                return
             end if
          end if
          tstep = tmp
       elseif (tstep > zero) then
          nsplit_real = real(dtime, real_kind)/(dt_max_factor*tstep)
          nsplit = idnint(nsplit_real)
          nstep_factor = dt_max_factor*nsplit
          if (abs(nsplit_real - nsplit) > divisible_tol*nsplit_real) then
             if (par%masterproc .and. .not. silent_in) then
                write(iulog,'(a,es11.4,a,i7,a,es11.4,a,i2,a,i2,a,i2,a)') &
                     'nsplit was computed as ', nsplit_real, ' based on dtime ', dtime, &
                     ', tstep ', tstep, &
                     ', and dt_max_factor = max(dt_remap_factor, dt_tracer_factor) = max(', &
                     dt_remap_factor, ',', dt_tracer_factor, ') =', dt_max_factor, &
                     ', which is outside the divisibility tolerance. &
                     &Set tstep, dt_remap_factor, and dt_tracer_factor so that &
                     &tstep and dt_max_factor*tstep divide dtime.'
             end if
             if (abort_in) call abortmp('timestep_make_parameters_consistent: divisibility error')
             return
          end if
       else
          if (par%masterproc .and. .not. silent_in) then
             write(iulog,*) 'If dtime is set to >0, then either nsplit or tstep must be >0.'
          end if
          if (abort_in) call abortmp('timestep_make_parameters_consistent: input error')
          return
       end if
    else
       if (tstep > zero) then
          if (nsplit > 0) then
             dtime = tstep*nstep_factor
          else
#ifdef CAM
             if (par%masterproc .and. .not. silent_in) then
                write(iulog,*) 'If dtime is set to <=0 and tstep >0, then nsplit must be >0.'
             end if
             if (abort_in) call abortmp('timestep_make_parameters_consistent: input error')
             return
#endif
          end if
       else
          if (par%masterproc .and. .not. silent_in) then
             write(iulog,*) 'If dtime is set to <=0, then tstep must be >0.'
          end if
          if (abort_in) call abortmp('timestep_make_parameters_consistent: input error')
          return          
       end if
    end if

    status = 0 ! success value
  end function timestep_make_eam_parameters_consistent

  subroutine test_timestep_make_parameters_consistent(par, nerr)
    ! Test timestep_make_parameters_consistent.

    use parallel_mod, only: parallel_t
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par
    integer, intent(out) :: nerr

    real(real_kind), parameter :: eps = epsilon(1.0_real_kind), tol = 1e3_real_kind*eps

    real(real_kind) :: tstep
    integer :: i, rs, qs, drf, dtf, ns, nstep_fac, dtime
    logical :: a, s

    a = .false.
    nerr = 0

    !! Test backwards compatibility.
    dtime = 1800

    qs = 3; rs = 0; drf = -1; dtf = -1
    tstep = -1; ns = 2
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a)
    if (i /= 0 .or. drf /= 0 .or. dtf /= qs .or. nstep_fac /= qs*ns .or. &
         abs(tstep - dtime/(ns*qs)) > tol) &
         nerr = nerr + 1

    qs = 3; rs = 2; drf = -1; dtf = -1
    tstep = -1; ns = 3
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a)
    if (i /= 0 .or. drf /= qs*rs .or. dtf /= qs .or. nstep_fac /= qs*rs*ns .or. &
         abs(tstep - dtime/(qs*rs*ns)) > tol) &
         nerr = nerr + 1

    !! Test new interface.
    tstep = 300_real_kind

    qs = -1; rs = -1; drf = 0; dtf = 6
    dtime = -1; ns = 2
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a)
    if (i /= 0 .or. rs /= drf .or. qs /= dtf .or. nstep_fac /= dtf*ns .or. &
         abs(tstep - dtime/(ns*qs)) > tol) &
         nerr = nerr + 1

    qs = -1; rs = -1; drf = 12; dtf = 6
    dtime = -1; ns = 3
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a)
    if (i /= 0 .or. rs /= 2 .or. qs /= 6 .or. nstep_fac /= qs*rs*ns .or. &
         abs(tstep - dtime/(qs*rs*ns)) > tol) &
         nerr = nerr + 1

    qs = -1; rs = -1; drf = 12; dtf = 6
    tstep = 300_real_kind; dtime = 7200; ns = -1
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a)
    if (i /= 0 .or. rs /= 2 .or. qs /= 6 .or. nstep_fac /= qs*rs*ns .or. ns /= 2 .or. &
         abs(tstep - dtime/(qs*rs*ns)) > tol) &
         nerr = nerr + 1

    !! Test new interface with new time step flexibility.
    qs = -1; rs = -1; drf = 2; dtf = 6
    dtime = -1; ns = 3
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a)
    if (i /= 0 .or. rs == 0 .or. qs /= 6 .or. nstep_fac /= dtf*ns .or. &
         abs(tstep - dtime/(dtf*ns)) > tol) &
         nerr = nerr + 1

    !! Test error and warning conditions.
    ! Silence messages in this unit test since we're forcing them.
    s = .true.

    qs = -1; rs = 0; drf = -1; dtf = -1
    tstep = -1; ns = 2
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i == 0) nerr = nerr + 1

    qs = -1; rs = 0; drf = 3; dtf = 4
    tstep = -1; ns = 2
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i == 0) nerr = nerr + 1

    qs = -1; rs = 0; drf = 1; dtf = 4
    tstep = 300; dtime = 700; ns = 1
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i == 0) nerr = nerr + 1

    qs = -1; rs = 0; drf = 1; dtf = 4
    tstep = 300; dtime = 700; ns = -1
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i == 0) nerr = nerr + 1

    qs = -1; rs = 0; drf = 1; dtf = 4
    tstep = -1; dtime = 700; ns = -1
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i == 0) nerr = nerr + 1

    qs = -1; rs = 0; drf = 1; dtf = 4
    tstep = -1; dtime = -1; ns = 2
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i == 0) nerr = nerr + 1

    !! Test warning conditions.
    qs = 4; rs = 0; drf = 3; dtf = 6
    tstep = -1; dtime = 1800; ns = 2
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i /= 0 .or. qs /= dtf .or. rs /= 1) nerr = nerr + 1

    qs = 4; rs = 0; drf = 12; dtf = 6
    tstep = -1; dtime = 1800; ns = 2
    i = timestep_make_parameters_consistent(par,rs,qs,drf,dtf,tstep,dtime,ns,nstep_fac,a,s)
    if (i /= 0 .or. qs /= dtf .or. rs /= 2) nerr = nerr + 1

    if (par%masterproc .and. nerr > 0) &
         write(iulog,'(a,i2)') 'test_timestep_make_parameters_consistent nerr', nerr
  end subroutine test_timestep_make_parameters_consistent


subroutine set_planar_defaults()

use physical_constants, only: Lx, Ly, Sx, Sy
 
!since defaults here depend on test, they cannot be set before ctl_nl is read, unlike some other parameters, bubble_*, etc.        
!if true, most likely lx,ly,sx,sy weren't set in ctl_nl
    if (      abs(lx).le.tol_zero .and. abs(ly).le.tol_zero &
        .and. abs(sx).le.tol_zero .and. abs(sy).le.tol_zero )then
    if (test_case == "planar_dbl_vrtx") then
      Lx = 5000.0D0 * 1000.0D0
      Ly = 5000.0D0 * 1000.0D0
      Sx = 0.0D0
      Sy = 0.0D0
    else if (test_case == "planar_hydro_gravity_wave") then
       Lx = 6000.0D0 * 1000.0D0
       Ly = 6000.0D0 * 1000.0D0
       Sx = -3000.0D0 * 1000.0D0
       Sy = -3000.0D0 * 1000.0D0
    else if (test_case == "planar_nonhydro_gravity_wave") then
       Lx = 300.0D0 * 1000.0D0
       Ly = 300.0D0 * 1000.0D0
       Sx = -150.0D0 * 1000.0D0
       Sy = -150.0D0 * 1000.0D0
    else if (test_case == "planar_hydro_mtn_wave") then
       Lx = 240.0D0 * 1000.0D0
       Ly = 240.0D0 * 1000.0D0
       Sx = 0.0D0 * 1000.0D0
       Sy = 0.0D0 * 1000.0D0
    else if (test_case == "planar_nonhydro_mtn_wave") then
       Lx = 144.0D0 * 1000.0D0
       Ly = 144.0D0 * 1000.0D0
       Sx = 0.0D0 * 1000.0D0
       Sy = 0.0D0 * 1000.0D0
    else if (test_case == "planar_schar_mtn_wave") then
       Lx = 100.0D0 * 1000.0D0
       Ly = 100.0D0 * 1000.0D0
       Sx = 0.0D0 * 1000.0D0
       Sy = 0.0D0 * 1000.0D0
    else if (test_case == "planar_density_current" .OR. test_case == "planar_moist_density_current") then
       Lx = 51.2D0 * 1000.0D0
       Ly = 51.2D0 * 1000.0D0
       Sx = -25.6D0 * 1000.0D0
       Sy = -25.6D0 * 1000.0D0
    else if (test_case == "planar_rising_bubble" ) then
       Lx = 2.0D0 * 10000.0D0
       Ly = 2.0D0 * 10000.0D0
       Sx = -10000.0D0
       Sy = -10000.0D0
! THESE ARE WRONG AND NEED TO BE FIXED WHEN THESE CASES ARE ACTUALLY IMPLEMENTED....
!else if (test_case == "planar_baroclinic_instab" .OR. test_case == "planar_moist_baroclinic_instab") then
!       Lx = 5000.0D0 * 1000.0D0
!       Ly = 5000.0D0 * 1000.0D0
!       Sx = 0.0D0
!       Sy = 0.0D0
!    else if (test_case == "planar_tropical_cyclone") then
!       Lx = 5000.0D0 * 1000.0D0
!       Ly = 5000.0D0 * 1000.0D0
!       Sx = 0.0D0
!       Sy = 0.0D0
!    else if (test_case == "planar_supercell") then
!       Lx = 5000.0D0 * 1000.0D0
!       Ly = 5000.0D0 * 1000.0D0
!       Sx = 0.0D0
!       Sy = 0.0D0

    endif
    endif !if lx,ly,sx,sy are not set in nl

    if (test_case == "planar_rising_bubble" ) then
       case_planar_bubble = .TRUE.
    end if

end subroutine set_planar_defaults

end module control_mod
