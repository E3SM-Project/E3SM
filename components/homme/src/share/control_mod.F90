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

  integer, public :: LFTfreq=0            ! leapfrog-trapazoidal frequency (shallow water only)
                                          ! interspace a lf-trapazoidal step every LFTfreq leapfrogs
                                          ! 0 = disabled

! vert_remap_q_alg:   -1  remap without monotone filter, used for some test cases
!                      0  default value, Zerroukat monotonic splines
!                      1  PPM vertical remap with mirroring at the boundaries
!                         (solid wall bc's, high-order throughout)
!                      2  PPM vertical remap without mirroring at the boundaries
!                         (no bc's enforced, first-order at two cells bordering top and bottom boundaries)
 integer, public :: vert_remap_q_alg = 0

! advect theta 0: conservation form
!              1: expanded divergence form (less noisy, non-conservative)
 integer, public :: theta_advect_form = 0

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
  character(len=MAX_STRING_LEN)    , public :: test_case      ! options: if cube: "swtc1","swtc2",or "swtc6"
  integer              , public :: tasknum
  integer              , public :: statefreq      ! output frequency of synopsis of system state (steps)
  integer              , public :: restartfreq
  integer              , public :: runtype
  integer              , public :: timerdetail
  integer              , public :: numnodes
  logical              , public :: uselapi
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
  character(len=MAX_STRING_LEN)    ,public  :: vform = ""     ! vertical coordinate system (sigma,hybrid)
  integer,                          public  :: vanalytic = 0  ! if 1, test initializes vertical coords
  real (kind=real_kind),            public  :: vtop = 0.1     ! top coordinate level for analytic vcoords

  integer              , public :: fine_ne = -1               ! set for refined exodus meshes (variable viscosity)
  real (kind=real_kind), public :: max_hypervis_courant = 1d99! upper bound for Courant number
                                                              ! (only used for variable viscosity, recommend 1.9 in namelist)
  real (kind=real_kind), public :: nu      = 7.0D5            ! viscosity (momentum equ)
  real (kind=real_kind), public :: nu_div  = -1               ! viscsoity (momentum equ, div component)
  real (kind=real_kind), public :: nu_s    = -1               ! default = nu   T equ. viscosity
  real (kind=real_kind), public :: nu_q    = -1               ! default = nu   tracer viscosity
  real (kind=real_kind), public :: nu_p    = -1               ! default = nu   ps equ. viscosity
  real (kind=real_kind), public :: nu_top  = 0.0D5            ! top-of-the-model viscosity

  integer, public :: hypervis_subcycle=1                      ! number of subcycles for hyper viscsosity timestep
  integer, public :: hypervis_subcycle_tom=0                  ! number of subcycles for TOM diffusion
                                                              !   0   apply together with hyperviscosity
                                                              !   >1  apply timesplit from hyperviscosity
  integer, public :: hypervis_subcycle_q=1                    ! number of subcycles for hyper viscsosity timestep on TRACERS
  integer, public :: hypervis_order=0                         ! laplace**hypervis_order.  0=not used  1=regular viscosity, 2=grad**4
  integer, public :: psurf_vis = 0                            ! 0 = use laplace on eta surfaces
                                                              ! 1 = use (approx.) laplace on p surfaces

  real (kind=real_kind), public :: hypervis_power=0           ! if not 0, use variable hyperviscosity based on element area
  real (kind=real_kind), public :: hypervis_scaling=0         ! use tensor hyperviscosity

  !three types of hyper viscosity are supported right now:
  ! (1) const hv:    nu * del^2 del^2
  ! (2) scalar hv:   nu(lat,lon) * del^2 del^2
  ! (3) tensor hv,   nu * ( \div * tensor * \grad ) * del^2
  !
  ! (1) default:  hypervis_power=0, hypervis_scaling=0
  ! (2) Original version for var-res grids. (M. Levy)
  !            scalar coefficient within each element
  !            hypervisc_scaling=0
  !            set hypervis_power>0 and set fine_ne, max_hypervis_courant
  ! (3) tensor HV var-res grids
  !            tensor within each element:
  !            set hypervis_scaling > 0 (typical values would be 3.2 or 4.0)
  !            hypervis_power=0
  !            (\div * tensor * \grad) operator uses cartesian laplace
  !

  ! hyperviscosity parameters used for smoothing topography
  integer, public :: smooth_phis_numcycle = 0   ! 0 = disable
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

contains

  function timestep_make_parameters_consistent(par, rsplit, qsplit, &
       dt_remap_factor, dt_tracer_factor, tstep, dtime, nsplit, nstep_factor, &
       abort, silent) result(status)

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
    logical, intent(in), optional :: abort, silent
    integer :: status

    real(kind=real_kind), parameter :: &
         zero = 0.0_real_kind, &
         eps = epsilon(1.0_real_kind), &
         divisible_tol = 1e3_real_kind*eps

    logical :: abort_in, silent_in

    abort_in = .true.
    if (present(abort)) abort_in = abort
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
       dt_remap_factor, dt_tracer_factor, abort, silent) result(status)

    use parallel_mod, only: abortmp, parallel_t
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par
    integer, intent(inout) :: rsplit, qsplit, dt_remap_factor, dt_tracer_factor
    logical, intent(in), optional :: abort, silent
    integer :: status

    integer :: qsplit_prev, rsplit_prev
    logical :: split_specified, factor_specified, split_is_master, abort_in, silent_in

    abort_in = .true.
    if (present(abort)) abort_in = abort
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
       nsplit, nstep_factor, tstep, dtime, abort, silent) result(status)

    use parallel_mod, only: abortmp, parallel_t
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par
    integer, intent(in) :: dt_remap_factor, dt_tracer_factor
    integer, intent(inout) :: nsplit
    integer, intent(out) :: nstep_factor
    real(kind=real_kind), intent(inout) :: tstep
    integer, intent(inout) :: dtime
    logical, intent(in), optional :: abort, silent
    integer :: status

    real(kind=real_kind), parameter :: &
         zero = 0.0_real_kind, &
         eps = epsilon(1.0_real_kind), &
         divisible_tol = 1e3_real_kind*eps

    real(kind=real_kind) :: nsplit_real, tmp
    integer :: dt_max_factor
    logical :: abort_in, silent_in

    abort_in = .true.
    if (present(abort)) abort_in = abort
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
                write(iulog,'(a,es11.4,a,i7,a,es11.4,a)') &
                     'nsplit was computed as ', nsplit_real, ' based on dtime ', dtime, &
                     ' and tstep ', tstep, ', which is outside the divisibility tolerance. Set &
                     &tstep so that it divides dtime.'
             end if
             if (abort_in) call abortmp('timestep_make_parameters_consistent: divisibility error')
             return
          end if
       else
          print *,'um>',par%rank,par%masterproc,silent_in,nsplit,dtime,tstep
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

end module control_mod
