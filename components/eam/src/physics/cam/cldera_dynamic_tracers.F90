!===============================================================================
! CLDERA dynamic tracers
! provides interfaces for registering and adding histories for dynamic tracer constituents
!
! Potential vorticity tracer described in:
! Whitehead, J.P. and Jablonowski, C. (2015) 
! Potential vorticity: Measuring consistency between GCM dynamical cores and tracer advection 
! schemes: Potential Vorticity: Measuring Consistency, Quarterly Journal of the Royal 
! Meteorological Society, doi:10.1002/qj.2389
!
!===============================================================================

module cldera_dynamic_tracers

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
  public :: cldera_dynamic_tracers_register        ! register constituents
  public :: cldera_dynamic_tracers_implements_cnst ! true if constituent is implemented by this package
  public :: cldera_dynamic_tracers_is_pt
  public :: cldera_dynamic_tracers_is_pv
  public :: cldera_dynamic_tracers_is_enabled
  public :: cldera_dynamic_tracers_init            ! initialize history fields, datasets
  public :: cldera_dynamic_tracers_init_cnst       ! initialize constituent field
  public :: cldera_dynamic_tracers_timestep_tend
  public :: cldera_dynamic_tracers_readnl          ! read namelist options
  public :: cldera_dynamic_tracers_pv_idx
  public :: cldera_dynamic_tracers_pt_idx
  public :: cldera_dynamic_tracers_pv_offset

  real(r8), parameter :: cldera_dynamic_tracers_pv_offset = 1_r8

  ! ----- Private module data
  integer, parameter :: ncnst=2  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'PV_TRCR_   ', 'PT_TRCR    '/)
  character(len=8), parameter :: pv_out_name = 'PV_TRCR    '

  integer :: ifirst ! global index of first constituent
  integer :: ixpv  ! global index for PV tracer
  integer :: ixpt  ! global index for PT tracer

  ! Data from namelist variables
  logical :: cldera_dynamic_tracers_flag  = .false. ! true => activate this module, namelist variable

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine cldera_dynamic_tracers_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use mpishorthand
    use cam_abortutils, only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cldera_dynamic_tracers_readnl'


    namelist /cldera_dynamic_tracers_nl/ cldera_dynamic_tracers_flag

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'cldera_dynamic_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, cldera_dynamic_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(cldera_dynamic_tracers_flag, 1, mpilog,  0, mpicom)
#endif

  endsubroutine cldera_dynamic_tracers_readnl

!================================================================================

  subroutine cldera_dynamic_tracers_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents
    ! 
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    if (.not. cldera_dynamic_tracers_flag) return

    ! arguments, left to right:
    ! - constituent name, moleculat weight, specific heat at const p, min. value allowed
    ! As far as we can tell, the molecular weight and specific heat are used for nothing
    ! by default for a newly added constituent
    ! The default 'mixtype' defaults to 'wet' when not passed to cnst_add; this also seems
    ! unused by default for a newly added constituent (this might not be true for WACCM?)
    call cnst_add(c_names(1), 0._r8, 0._r8, 0._r8, ixpv,  readiv=.false., &
                  longname='Potential vorticity tracer', cam_outfld=.false.)
    ifirst = ixpv



    call cnst_add(c_names(2), 0._r8, 0._r8, 0._r8, ixpt,  readiv=.false., &
                  longname='Potential temperature tracer')


  end subroutine cldera_dynamic_tracers_register

!===============================================================================

  function cldera_dynamic_tracers_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: cldera_dynamic_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    cldera_dynamic_tracers_implements_cnst = .false.

    if (.not. cldera_dynamic_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          cldera_dynamic_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function cldera_dynamic_tracers_implements_cnst

  subroutine cldera_dynamic_tracers_pv_idx(idx)
    integer, intent(out) :: idx

    idx = 0
    if (.not. cldera_dynamic_tracers_flag) return

    idx = ixpv


  end subroutine cldera_dynamic_tracers_pv_idx
  subroutine cldera_dynamic_tracers_pt_idx(idx)
    integer, intent(out) :: idx

    idx = 0
    if (.not. cldera_dynamic_tracers_flag) return

    idx = ixpt


  end subroutine cldera_dynamic_tracers_pt_idx


  subroutine cldera_dynamic_tracers_timestep_tend(state, ptend, dt, ncol)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,     only: get_rlat_all_p , get_lat_all_p
    use physconst,     only: gravit, avogad
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep
    use ref_pres,      only: pref_mid_norm
    use time_manager,  only: get_curr_time

    ! Arguments
    type(physics_state), intent(inout) :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(in)    :: dt                 ! timestep
    integer,             intent(in)    :: ncol

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk             ! chunk identifier

    real(r8) :: pv_out(ncol,pver)
    if (.not. cldera_dynamic_tracers_flag) then
      return
    end if
    lchnk = state%lchnk
    do i=1,ncol
      do k=1,pver
        pv_out(i,k) = state%q(i,k,ixpv)
      end do
    end do
    call outfld(pv_out_name, pv_out(:,:)-cldera_dynamic_tracers_pv_offset, ncol, lchnk)           

  end subroutine cldera_dynamic_tracers_timestep_tend

!


! ===========================================================================

  function cldera_dynamic_tracers_is_pv(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is treated as pv
    ! 
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: cldera_dynamic_tracers_is_pv        ! return value

    !-----------------------------------------------------------------------

    cldera_dynamic_tracers_is_pv = (cldera_dynamic_tracers_flag) .and. (name == c_names(1))


  end function cldera_dynamic_tracers_is_pv

! ==========================================================================

 function cldera_dynamic_tracers_is_enabled()
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is treated as pv
    ! 
    !-----------------------------------------------------------------------

    logical :: cldera_dynamic_tracers_is_enabled        ! return value

    !-----------------------------------------------------------------------

    cldera_dynamic_tracers_is_enabled = cldera_dynamic_tracers_flag


  end function cldera_dynamic_tracers_is_enabled

! ==========================================================================

 
  function cldera_dynamic_tracers_is_pt(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is treated as pt
    ! 
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: cldera_dynamic_tracers_is_pt        ! return value


    cldera_dynamic_tracers_is_pt = (cldera_dynamic_tracers_flag) .and. (name == c_names(2))


  end function cldera_dynamic_tracers_is_pt


!===============================================================================
 

!===============================================================================
  
  subroutine cldera_dynamic_tracers_init

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default, horiz_only

    integer :: m, mm
    !-----------------------------------------------------------------------

    if (.not. cldera_dynamic_tracers_flag) return

    ! decalre PV, PT tracers as history variables
    call addfld (pv_out_name, (/ 'lev' /), 'I', 'm2 K/kg/s', cnst_longname(ixpv))
    call add_default (pv_out_name, 1, ' ')
    call addfld (cnst_name(ixpt), (/ 'lev' /), 'I', 'K', cnst_longname(ixpt))
    call add_default (cnst_name(ixpt), 1, ' ')
    
  end subroutine cldera_dynamic_tracers_init

!===============================================================================

  subroutine cldera_dynamic_tracers_init_cnst(name, q, gcid)

    !----------------------------------------------------------------------- 
    !
    ! Purpose: initialize tracers mixing ratio fields 
    !          This subroutine provided for testing only; PV, PT tracers
    !          need to be initialized by the dycore. This routine can be used
    !          instead to set tracers fields to zero. If so, this subroutine 
    !          is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(out) :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. cldera_dynamic_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, q, gcid)
       endif
    end do

  end subroutine cldera_dynamic_tracers_init_cnst

!===============================================================================

  subroutine init_cnst_3d(m, q, gcid)

    use dyn_grid,    only : get_horiz_grid_d, get_horiz_grid_dim_d
    use dycore,      only : dycore_is
    use ref_pres,    only : pref_mid_norm

    integer,  intent(in)  :: m       ! global constituent index
    real(r8), intent(out) :: q(:,:)  ! kg tracer/kg dry air (gcol,plev)
    integer,  intent(in)  :: gcid(:) ! global column id

    real(r8), allocatable :: lat(:)
    integer :: plon, plat, ngcols
    integer :: j, k, gsize
    !-----------------------------------------------------------------------
    ! initialization below is on the DYNAMICS grid; length of gcol will 
    ! be number of GLL points, not number of physics columns

    if (masterproc) write(iulog,*) 'CLDERA DYNAMIC CONSTITUENTS: INITIALIZING ',cnst_name(m),m

    ! ====== PV ======
    if (m == ixpv) then
       q(:,:) = 0.0_r8

    ! ====== PT ======
    else if (m == ixpt) then
       q(:,:) = 0.0_r8
    end if

  end subroutine init_cnst_3d

!=====================================================================


end module cldera_dynamic_tracers
