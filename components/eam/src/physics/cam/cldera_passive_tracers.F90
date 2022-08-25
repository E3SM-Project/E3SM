!===============================================================================
! CLDERA passive tracers
! provides dissipation rates and sources for diagnostic constituents
!
! AOA tracer implemented as defined in:
! Gupta, A., Gerber, E.P. and Lauritzen, P.H. (2020) 
! Numerical impacts on tracer transport: A proposed intercomparison test of Atmospheric 
! General Circulation Models, Quarterly Journal of the Royal Meteorological Society,
! doi:10.1002/qj.3881.
!
! E90 tracer implemented as defined in:
! Abalos, M. et al. (2017)
! Using the Artificial Tracer e90 to Examine Present and Future UTLS Tracer Transport in WACCM, 
! Journal of the Atmospheric Sciences,
! doi:10.1175/JAS-D-17-0135.1.
! 
! ST80_25 tracer implemented as defined in
! Eyring, V. et al. (2013)
! Overview of IGAC/SPARC Chemistry-Climate Model Initiative (CCMI) Community Simulations in 
! Support of Upcoming Ozone and Climate Assessments
!
!===============================================================================

module cldera_passive_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog
  use ref_pres,     only: pref_mid_norm

  implicit none
  private
  save

  ! Public interfaces
  public :: cldera_passive_tracers_register         ! register constituents
  public :: cldera_passive_tracers_implements_cnst  ! true if constituent is implemented by this package
  public :: cldera_passive_tracers_init_cnst        ! initialize constituent field
  public :: cldera_passive_tracers_init             ! initialize history fields, datasets
  public :: cldera_passive_tracers_timestep_init    ! place to perform per timestep initialization
  public :: cldera_passive_tracers_timestep_tend    ! calculate tendencies
  public :: cldera_passive_tracers_readnl           ! read namelist options

  ! Private module data

  integer, parameter :: ncnst=3  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'AOA', 'E90', 'ST80_25'/)

  integer :: ifirst ! global index of first constituent
  integer :: ixaoa  ! global index for AOA tracer
  integer :: ixe90  ! global index for E90 tracer
  integer :: ixst80 ! global index for ST80_25 tracer

  ! Data from namelist variables
  logical :: cldera_passive_tracers_flag  = .false.    ! true => turn on test tracer code, namelist variable
  logical :: cldera_passive_read_from_ic_file = .true. ! true => tracers initialized from IC file

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine cldera_passive_tracers_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cldera_passive_tracers_readnl'


    namelist /cldera_passive_tracers_nl/ cldera_passive_tracers_flag, cldera_passive_read_from_ic_file

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'cldera_passive_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, cldera_passive_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(cldera_passive_tracers_flag, 1, mpilog,  0, mpicom)
    call mpibcast(cldera_passive_read_from_ic_file, 1, mpilog,  0, mpicom)
#endif

  endsubroutine cldera_passive_tracers_readnl

!================================================================================

  subroutine cldera_passive_tracers_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents
    ! 
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixaoa,  readiv=cldera_passive_read_from_ic_file, &
                  longname='Age-of-air tracer')
    ifirst = ixaoa
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixe90,  readiv=cldera_passive_read_from_ic_file, &
                  longname='E90 tracer')
    call cnst_add(c_names(3), mwdry, cpair, 0._r8, ixst80, readiv=cldera_passive_read_from_ic_file, &
                  longname='ST80_25 tracer')

  end subroutine cldera_passive_tracers_register

!===============================================================================

  function cldera_passive_tracers_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: cldera_passive_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    cldera_passive_tracers_implements_cnst = .false.

    if (.not. cldera_passive_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          cldera_passive_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function cldera_passive_tracers_implements_cnst

!===============================================================================

  subroutine cldera_passive_tracers_init_cnst(name, q, gcid)

    !----------------------------------------------------------------------- 
    !
    ! Purpose: initialize test tracers mixing ratio fields 
    !  This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(out) :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, q, gcid)
       endif
    end do

  end subroutine cldera_passive_tracers_init_cnst

!===============================================================================

  subroutine cldera_passive_tracers_init

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default

    integer :: m, mm
    !-----------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) return

    ! Set names of tendencies and declare them as history variables

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld (cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))
       call add_default (cnst_name(mm), 1, ' ')
    end do

  end subroutine cldera_passive_tracers_init

!===============================================================================

  subroutine cldera_passive_tracers_timestep_init( phys_state )
    !-----------------------------------------------------------------------
    ! Provides a place to reinitialize diagnostic constituents
    ! Currently does nothing
    !-----------------------------------------------------------------------

    use time_manager,   only: get_curr_date
    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state

    type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state    


    integer c, i, k, ncol
    integer yr, mon, day, tod
    !--------------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) return

    call get_curr_date (yr,mon,day,tod)

    if ( day == 1 .and. tod == 0) then
       if (masterproc) then
         write(iulog,*) 'CLDERA_CLOCK_CONSTITUENTS: RE-INITIALIZING CONSTITUENTS'
       endif

       do c = begchunk, endchunk
          ncol = phys_state(c)%ncol
          do k = 1, pver
             do i = 1, ncol
                ! does nothing currently
                continue
             end do
          end do
       end do

    end if

  end subroutine cldera_passive_tracers_timestep_init

!===============================================================================

  subroutine cldera_passive_tracers_timestep_tend(state, ptend, dt)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,     only: get_rlat_all_p , get_lat_all_p
    use physconst,     only: gravit
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep
    use ref_pres,      only: pref_mid_norm
    use time_manager,  only: get_curr_time

    ! Arguments
    type(physics_state), intent(inout)    :: state           ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(in)    :: dt                 ! timestep
    !real(r8),           intent(inout) :: cflx(pcols,pcnst)  ! Surface constituent flux (kg/m^2/s)

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk             ! chunk identifier
    integer :: ncol              ! no. of column in chunk
    integer :: nstep             ! current timestep number

    logical  :: lq(pcnst)

    integer  :: day,sec          ! date variables
    real(r8) :: t                ! tracer boundary condition
    real(r8) :: aoa_scaling      ! scale AOA1 from nstep to time

    real(r8) :: efold_st80       ! e-folding timescale for e90 in s
    real(r8) :: efold_e90        ! e-folding timescale for e90 in s
    real(r8) :: mweight_e90      ! molecular weight of E90 in g/mol
    real(r8) :: sflx_e90         ! surface flux of E90 in mol cm^-2 s^-1

    !------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) then
       !Initialize an empty ptend for use with physics_update
       call physics_ptend_init(ptend,state%psetcols,'cldera_passive_trc_ts')
       return
    end if

    lq(:)      = .FALSE.
    lq(ixaoa)  = .TRUE.
    lq(ixe90)  = .TRUE.
    lq(ixst80) = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'cldera_passive_tracers', lq=lq)

    nstep = get_nstep()
    lchnk = state%lchnk
    ncol  = state%ncol

    ! ---- compute AOA time scaling (1 s in days)
    aoa_scaling = 1._r8/86400._r8

    ! ---- e-folding times for e90 (90 days in s), st80 (25 days in s)
    ! (reciprocals match WACCM 'E90_tau' and 'ST80_25_tau', and table A2 in 
    ! Tilmes+ (2016) Representation of CESM1 CAM4-chem within CCMI)
    efold_e90  = 90._r8 * 86400._r8
    efold_st80 = 25._r8 * 86400._r8

    ! ---- molecular weight and surface flux for e90
    ! (surface flux matches Abalos+ (2017), molecular weight is for CO, matches
    !  WACCM surface emissions specification in 
    !  emissions_E90global_surface_1750-2100_0.9x1.25_c20170322.nc)
    mweight_e90 = 28._r8                                   ! molecular weight of CO [g/mole]
    sflx_e90    = 2.7736e11_r8                             ! [molecules cm^-2 s^-1]
    sflx_e90    = sflx_e90 / 6.022141e23_r8                ! [moles cm^-2 s^-1]
    ! convert mol to g, g to kg, cm^-2 to m^-2, ==>  flux in [kg m^-2 s^-1]
    sflx_e90    = sflx_e90 * mweight_e90 * (1._r8/1000._r8) * 10000._r8


    ! -------------------- TRACER TENDENCIES --------------------
    do k = 1, pver
       do i = 1, ncol
       ! the explicit horizontal looping is not currently needed; maintained 
       ! in case spatial dependence in e.g. AOA is ever desired 

          ! ============ AOA ============
          ! clock tracer with a source of 1 day/day everywhere above ~700hPa
          if (pref_mid_norm(k) <= 0.7) then
              ptend%q(i,k,ixaoa) = 1.0_r8 * aoa_scaling
          else
              ptend%q(i,k,ixaoa) = 0.0_r8
              state%q(i,k,ixaoa) = 0.0_r8
          end if

          ! ============ E90 ============
          ! dissipates with e-folding time of 90 days
          ptend%q(i,k,ixe90) = -(1._r8 / efold_e90) * state%q(i, k, ixe90)
          ! apply surface flux to lowest model level
          if(k == pver) then
              state%q(i, k, ixe90) = state%q(i, k, ixe90) + &
                                     dt * gravit * 1._r8 / state%pdel(i, k) * sflx_e90
          end if

          ! ============ ST80_25 ============
          ! dissipates with e-folding time of 25 days below ~80 hPa, 
          ! constant concentration of 200 ppbv above
          if (pref_mid_norm(k) >= 0.08) then
              ptend%q(i,k,ixst80) = -(1._r8 / efold_st80) * state%q(i, k, ixst80)
          else
              ptend%q(i, k,ixst80) = 0._r8
              state%q(i, k, ixst80) = 200e-9_r8  ! 200 ppbv to vmr
          end if
       end do
    end do

    ! -------------------- TRACER FLUXES --------------------
    do i = 1, ncol

       ! ====== AOA ======
       ! no surface flux
       !cflx(i,ixaoa) = 0._r8

       ! ====== E90 ======
       !cflx(i,ixe90) = sflx_e90

       ! ====== ST80_25 ======
       ! no surface flux
       !cflx(i,ixst80) = 0._r8
       continue

    end do

  end subroutine cldera_passive_tracers_timestep_tend

!===========================================================================

  subroutine init_cnst_3d(m, q, gcid)

    use dyn_grid, only : get_horiz_grid_d, get_horiz_grid_dim_d
    use dycore,   only : dycore_is

    integer,  intent(in)  :: m       ! global constituent index
    real(r8), intent(out) :: q(:,:)  ! kg tracer/kg dry air (gcol,plev)
    integer,  intent(in)  :: gcid(:) ! global column id

    real(r8), allocatable :: lat(:)
    integer :: plon, plat, ngcols
    integer :: j, k, gsize
    !-----------------------------------------------------------------------

    if (masterproc) write(iulog,*) 'CLDERA_CLOCK CONSTITUENTS: INITIALIZING ',cnst_name(m),m

    ! ====== AOA ======
    if (m == ixaoa) then

       q(:,:) = 0.0_r8

    ! ====== E90 ======
    else if (m == ixe90) then

       q(:,:) = 0.0_r8

    ! ====== ST80_25 ======
    else if (m == ixst80) then

       q(:,:) = 0.0_r8

    end if

  end subroutine init_cnst_3d

!=====================================================================


end module cldera_passive_tracers
