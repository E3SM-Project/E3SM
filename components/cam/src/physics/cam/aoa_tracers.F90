!===============================================================================
! Age of air test tracers
! provides dissipation rate and surface fluxes for diagnostic constituents
!===============================================================================

module aoa_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog

  implicit none
  private
  save

  ! Public interfaces
  public :: aoa_tracers_register         ! register constituents
  public :: aoa_tracers_implements_cnst  ! true if named constituent is implemented by this package
  public :: aoa_tracers_init_cnst        ! initialize constituent field
  public :: aoa_tracers_init             ! initialize history fields, datasets
  public :: aoa_tracers_timestep_init    ! place to perform per timestep initialization
  public :: aoa_tracers_timestep_tend    ! calculate tendencies
  public :: aoa_tracers_readnl           ! read namelist options

  ! Private module data

  integer, parameter :: ncnst=4  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'AOA1', 'AOA2', 'HORZ', 'VERT'/)

  ! constituent source/sink names
  character(len=8), parameter :: src_names(ncnst) = (/'AOA1SRC', 'AOA2SRC', 'HORZSRC', 'VERTSRC'/)

  integer :: ifirst ! global index of first constituent
  integer :: ixaoa1 ! global index for AOA1 tracer
  integer :: ixaoa2 ! global index for AOA2 tracer
  integer :: ixht   ! global index for HORZ tracer
  integer :: ixvt   ! global index for VERT tracer

  ! Data from namelist variables
  logical :: aoa_tracers_flag  = .false.    ! true => turn on test tracer code, namelist variable
  logical :: aoa_read_from_ic_file = .true. ! true => tracers initialized from IC file

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine aoa_tracers_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aoa_tracers_readnl'


    namelist /aoa_tracers_nl/ aoa_tracers_flag, aoa_read_from_ic_file

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aoa_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aoa_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(aoa_tracers_flag, 1, mpilog,  0, mpicom)
    call mpibcast(aoa_read_from_ic_file, 1, mpilog,  0, mpicom)
#endif

  endsubroutine aoa_tracers_readnl

!================================================================================

  subroutine aoa_tracers_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents
    ! 
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixaoa1, readiv=aoa_read_from_ic_file, &
                  longname='Age-of_air tracer 1')
    ifirst = ixaoa1
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixaoa2, readiv=aoa_read_from_ic_file, &
                  longname='Age-of_air tracer 2')
    call cnst_add(c_names(3), mwdry, cpair, 1._r8, ixht,   readiv=aoa_read_from_ic_file, &
                  longname='horizontal tracer')
    call cnst_add(c_names(4), mwdry, cpair, 0._r8, ixvt,   readiv=aoa_read_from_ic_file, &
                  longname='vertical tracer')

  end subroutine aoa_tracers_register

!===============================================================================

  function aoa_tracers_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: aoa_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    aoa_tracers_implements_cnst = .false.

    if (.not. aoa_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          aoa_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function aoa_tracers_implements_cnst

!===============================================================================

  subroutine aoa_tracers_init_cnst(name, q, gcid)

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

    if (.not. aoa_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, q, gcid)
       endif
    end do

  end subroutine aoa_tracers_init_cnst

!===============================================================================

  subroutine aoa_tracers_init

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default

    integer :: m, mm
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    ! Set names of tendencies and declare them as history variables

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld (cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))
       call addfld (src_names(m), (/ 'lev' /), 'A',  'kg/kg/s', trim(cnst_name(mm))//' source/sink')

       call add_default (cnst_name(mm), 1, ' ')
       call add_default (src_names(m),  1, ' ')
    end do

  end subroutine aoa_tracers_init

!===============================================================================

  subroutine aoa_tracers_timestep_init( phys_state )
    !-----------------------------------------------------------------------
    ! Provides a place to reinitialize diagnostic constituents HORZ and VERT
    !-----------------------------------------------------------------------

    use time_manager,   only: get_curr_date
    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state

    type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state    


    integer c, i, k, ncol
    integer yr, mon, day, tod
    !--------------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    call get_curr_date (yr,mon,day,tod)

    if ( day == 1 .and. tod == 0) then
       if (masterproc) then
         write(iulog,*) 'AGE_OF_AIR_CONSTITUENTS: RE-INITIALIZING HORZ/VERT CONSTITUENTS'
       endif

       do c = begchunk, endchunk
          ncol = phys_state(c)%ncol
          do k = 1, pver
             do i = 1, ncol
                phys_state(c)%q(i,k,ixht) = 2._r8 + sin(phys_state(c)%lat(i))
                phys_state(c)%q(i,k,ixvt) = real(pver-k+1,r8)
             end do
          end do
       end do

    end if

  end subroutine aoa_tracers_timestep_init

!===============================================================================

  subroutine aoa_tracers_timestep_tend(state, ptend, cflx, landfrac, dt)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,     only: get_rlat_all_p , get_lat_all_p
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep

    ! Arguments
    type(physics_state), intent(in)    :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(inout) :: cflx(pcols,pcnst)  ! Surface constituent flux (kg/m^2/s)
    real(r8),            intent(in)    :: landfrac(pcols)    ! Land fraction
    real(r8),            intent(in)    :: dt                 ! timestep

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk                          ! chunk identifier
    integer :: ncol                           ! no. of column in chunk
    integer :: nstep                          ! current timestep number

    real(r8), parameter ::  teul=3.9e-7_r8       ! rate
    real(r8) :: qrel                          ! value to be relaxed to
    real(r8) :: xhorz                         ! updated value of HORZ
    real(r8) :: xvert                         ! updated value of VERT

    logical  :: lq(pcnst)
    !------------------------------------------------------------------



    if (.not. aoa_tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    lq(:)      = .FALSE.
    lq(ixaoa1) = .TRUE.
    lq(ixaoa2) = .TRUE.
    lq(ixht)   = .TRUE.
    lq(ixvt)   = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'aoa_tracers', lq=lq)

    nstep = get_nstep()
    lchnk = state%lchnk
    ncol  = state%ncol

    do k = 1, pver
       do i = 1, ncol

          ! AOA1
          ptend%q(i,k,ixaoa1) = 0.0_r8

          ! AOA2
          ptend%q(i,k,ixaoa2) = 0.0_r8

          ! HORZ
          qrel              = 2._r8 + sin(state%lat(i))
          xhorz             = (state%q(i,k,ixht) + dt*teul*qrel)/ (1._r8 + teul * dt)
          ptend%q(i,k,ixht) = (xhorz - state%q(i,k,ixht)) / dt

          ! VERT
          qrel              = real(pver-k+1,r8)
          xvert             = (state%q(i,k,ixvt) + dt*teul*qrel)/ (1._r8 + teul * dt)
          ptend%q(i,k,ixvt) = (xvert - state%q(i,k,ixvt)) / dt

       end do
    end do

    ! record tendencies on history files
    call outfld (src_names(1), ptend%q(:,:,ixaoa1), pcols, lchnk)
    call outfld (src_names(2), ptend%q(:,:,ixaoa2), pcols, lchnk)
    call outfld (src_names(3), ptend%q(:,:,ixht),   pcols, lchnk)
    call outfld (src_names(4), ptend%q(:,:,ixvt),   pcols, lchnk)

    ! Set tracer fluxes
    do i = 1, ncol

       ! AOA1
       cflx(i,ixaoa1) = 1.e-6_r8

       ! AOA2
       if (landfrac(i) .eq. 1._r8  .and.  state%lat(i) .gt. 0.35_r8) then
          cflx(i,ixaoa2) = 1.e-6_r8 + 1e-6_r8*0.0434_r8*real(nstep,r8)*dt/(86400._r8*365._r8)
       else
          cflx(i,ixaoa2) = 0._r8
       endif

       ! HORZ
       cflx(i,ixht) = 0._r8

       ! VERT
       cflx(i,ixvt) = 0._r8

    end do

  end subroutine aoa_tracers_timestep_tend

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

    if (masterproc) write(iulog,*) 'AGE-OF-AIR CONSTITUENTS: INITIALIZING ',cnst_name(m),m

    if (m == ixaoa1) then

       q(:,:) = 0.0_r8

    else if (m == ixaoa2) then

       q(:,:) = 0.0_r8

    else if (m == ixht) then

       call get_horiz_grid_dim_d( plon, plat )
       ngcols = plon*plat
       gsize = size(gcid)
       allocate(lat(ngcols))
       call get_horiz_grid_d(ngcols,clat_d_out=lat)
       do j = 1, gsize
          q(j,:) = 2._r8 + sin(lat(gcid(j)))
       end do
       deallocate(lat)

    else if (m == ixvt) then

       do k = 1, pver
          do j = 1, size(q,1)
             q(j,k) = real(pver-k+1,r8)
          end do
       end do

    end if

  end subroutine init_cnst_3d

!=====================================================================


end module aoa_tracers
