!
module gw_drag_xjb
!---------------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an 
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use spmd_utils,    only: masterproc
!=======Jinbo Xie=======
  use ppgrid,        only: pcols,pver,pverp,&
                           nvar_dirOA,nvar_dirOL,indexb,begchunk,endchunk
!=======Jinbo Xie=======
  !Jinbo Xie add another, pverp=pver+1
  use physics_types, only: physics_state, physics_ptend
  use cam_history,   only: outfld
  use scamMod,       only: single_column
  use cam_logfile,   only: iulog
  use abortutils,    only: endrun
  use gw_drag,       only: fcrit2       !czy20181120

!===========================
!   Jinbo Xie1 modification
!===========================
        use physconst,      only: rh2o,zvir,pi,rearth,r_universal  !zvir is the ep1 in wrf,rearth is the radius of earth (m),r_universal is the gas constant
!czy20181120        use module_bl_gwdo_gsd, only: gwdo_gsd,gwdo2d
!===========
! Jinbo Xie
!===========
!kpbl2d(!Layer index containing PBL top within or at the base interface)
        use hycoef,             only: hyai, hybi, hyam, hybm, etamid !get the znu,znw,p_top set to 0
        use hb_diff,            only: pblintd!rino_pub
!===========================
!   Jinbo Xie1 modification
!===========================


  implicit none
  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
!czy20181120  public gw_drag_readnl           ! Read namelist
  public gw_inti                  ! Initialization
  public gw_intr                  ! interface to actual parameterization

!+czy20181120 xjb
  public gwdo_gsd
  public gwdo2d
!-czy20181120 xjb
!
! PRIVATE: Rest of the data and interfaces are private to this module
!
!czy20181120  real(r8), parameter :: unset_r8 = huge(1._r8)

  ! fcrit2 has been made a namelist variable to facilitate backwards compatibility 
  ! with the CAM3 version of this parameterization.  In CAM3 fcrit2=0.5
!czy20181120  real(r8) :: fcrit2 = unset_r8   ! critical froude number squared
  integer, parameter :: pgwv = 0  ! number of waves allowed

  integer :: kbotbg, kbotoro      ! interface of gwd source
  integer :: ktopbg, ktoporo      ! top interface of gwd region

  real(r8) :: alpha(0:pver)       ! newtonian cooling coefficients
  real(r8) :: c(-pgwv:pgwv)       ! list of wave phase speeds
  real(r8) :: cpair               ! specific heat of dry air (constant p)
  real(r8) :: dback               ! background diffusivity
  real(r8) :: effkwv              ! effective wavenumber (fcrit2*kwv)
  real(r8) :: effgw_oro           ! tendency efficiency for orographic gw
  real(r8) :: effgw_spec=.125_r8  ! tendency efficiency for internal gw
  real(r8) :: fracldv             ! fraction of stress deposited in low level region
  real(r8) :: g                   ! acceleration of gravity
!==============Jinbo Xie===========
!   Jinbo Xie2 modification
!==================================
        !simply add par
        !for z,dz,from other files
        real(r8) :: ztop(pcols,pver)             ! top interface height asl (m)
        real(r8) :: zbot(pcols,pver)             ! bottom interface height asl (m)
        real(r8) :: zmid(pcols,pver)             ! middle interface height asl (m)
        real(r8) :: dz(pcols,pver)       ! model layer height

        !bulk richardson number from hb_diff
        !bulk at the surface
        !real(r8),parameter :: rino(pcols,nver)
        real(r8) :: rlat(pcols)

!==============Jinbo Xie===========
!   Jinbo Xie2 modification
!==================================


  real(r8) :: kwv                 ! effective horizontal wave number
  real(r8) :: lzmax               ! maximum vertical wavelength at source

  real(r8) :: mxasym              ! max asymmetry between tau(c) and tau(-c)
  real(r8) :: mxrange             ! max range of tau for all c
  real(r8) :: n2min               ! min value of bouyancy frequency
  real(r8) :: oroko2              ! 1/2 * horizontal wavenumber
  real(r8) :: orohmin             ! min surface displacment height for orographic waves
  real(r8) :: orovmin             ! min wind speed for orographic waves
  real(r8) :: r                   ! gas constant for dry air
  real(r8) :: rog                 ! r / g
  real(r8) :: taubgnd             ! background source strength (/tauscal)
  real(r8) :: taumin              ! minimum (nonzero) stress
  real(r8) :: tauscal             ! scale factor for background stress source
  real(r8) :: tndmax              ! maximum wind tendency
  real(r8) :: umcfac              ! factor to limit tendency to prevent reversing u-c
  real(r8) :: ubmc2mn             ! min (u-c)**2
  real(r8) :: zldvcon             ! constant for determining zldv from tau0

!===============================================================================
contains
!===============================================================================

!+czy20181120============================================================
#if 0
subroutine gw_drag_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'gw_drag_readnl'

   namelist /gw_drag_nl/ fcrit2
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'gw_drag_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, gw_drag_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(fcrit2, 1, mpir8, 0, mpicom)
#endif

   ! Error checking:

   ! check if fcrit2 was set.
   if (fcrit2 == unset_r8) then
      call endrun('gw_dray_readnl: fcrit2 must be set via the namelist')
   end if

end subroutine gw_drag_readnl
#endif
!-czy20181120============================================================

!================================================================================

  subroutine gw_inti (cpairx, cpwv, gx, rx, hypi)
!-----------------------------------------------------------------------
! Time independent initialization for multiple gravity wave parameterization.
!-----------------------------------------------------------------------
    use cam_history, only: addfld, add_default, phys_decomp
    use dycore,      only: dycore_is
    use phys_control,only: phys_getopts
!=======Jinbo Xie=======
    use startup_initialconds, only:topoGWD_file_get_id,setup_initialGWD,close_initial_fileGWD
#ifndef continuous
    use comsrf,              only:var,oc,oadir,ol,initialize_comsrf2!,dxydir
#else
    use comsrf,              only: var,oc,oadir,terrout,initialize_comsrf2
#endif
    !use ppgrid,              only: begchunk, endchunk,! pcols,nvar_dirOA,nvar_dirOL
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld
!=======Jinbo Xie=======
!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: cpairx                ! specific heat of dry air (constant p)
    real(r8), intent(in) :: cpwv                  ! specific heat of water vapor (constant p)
    real(r8), intent(in) :: gx                    ! acceleration of gravity
    real(r8), intent(in) :: rx                    ! gas constant for dry air
    real(r8), intent(in) :: hypi(pver+1)          ! reference interface pressures

!---------------------------Local storage-------------------------------
    integer :: k
    logical :: history_budget                     ! output tendencies and state variables for CAM4
                                                  ! temperature, water vapor, cloud ice and cloud
                                                  ! liquid budgets.
    integer :: history_budget_histfile_num        ! output history file number for budget fields

!-----------------------------------------------------------------------
!=================Jinbo Xie-3-18-2019============
  type(file_desc_t), pointer :: ncid_topoGWD  !Jinbo Xie add
  character(len=8) :: terroutchar(4)
  logical :: found=.false., found2=.false.
  integer :: i
!================================================
!=================Jinbo Xie-3-18-2019============
!#if 0
                call initialize_comsrf2()
                call setup_initialGWD()
                ncid_topoGWD=>topoGWD_file_get_id()
                call infld('SGH' ,ncid_topoGWD,'lon','lat',1,pcols,begchunk,&
                                endchunk,  var, found, grid_map='PHYS')
                call infld('OC', ncid_topoGWD,'lon','lat',1,pcols,begchunk, &
                                endchunk,  oc,  found, grid_map='PHYS')
                !keep the same interval of OA,OL
                call infld('OA', ncid_topoGWD,'lon','nvar_dirOA','lat',1,pcols,1,nvar_dirOA,begchunk, &
                                endchunk,  oadir(:,:,:),  found, grid_map='PHYS')
                call infld('OL', ncid_topoGWD,'lon','nvar_dirOL','lat',1,pcols,1,nvar_dirOL,begchunk, &
                                endchunk,  ol, found, grid_map='PHYS')
                if(.not. found) call endrun('ERROR: GWD topo file readerr')
                call close_initial_fileGWD()
!#endif
!=================Jinbo Xie readtopo-3-18-2019============
! Copy model constants
    cpair  = cpairx
    g      = gx
    r      = rx

! Set MGWD constants
    kwv    = 6.28e-5_r8          ! 100 km wave length
    dback  = 0.05_r8             ! background diffusivity
    tauscal= 0.001_r8            ! scale factor for background stress
    taubgnd= 6.4_r8              ! background stress amplitude

    zldvcon= 10._r8              ! constant for determining zldv
    lzmax  = 7.E3_r8             ! maximum vertical wavelength at source (m)
    if ( dycore_is ('LR') ) then
       effgw_oro = 0.125_r8
       fracldv= 0.0_r8
! Restore these to work with turulent mountain stress
!       effgw_oro = 1.0
!       fracldv= 0.7           ! fraction of tau0 diverged in low level region
    else
       effgw_oro = 0.125_r8
       fracldv= 0.0_r8
    endif
! Set phase speeds 
    do k = -pgwv, pgwv
       c(k)   = 10._r8 * k       ! 0, +/- 10, +/- 20, ... m/s
    end do

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'GW_INTI: pgwv = ', pgwv
       write(iulog,*) 'GW_INTI: c(l) = ', c
       write(iulog,*) 'GW_INTI: effgw_oro = ', effgw_oro
       write(iulog,*) 'GW_INTI: kwv = ', kwv
       write(iulog,*) 'GW_INTI: fcrit2 = ', fcrit2
       write(iulog,*) ' '
    end if

! Set radiative damping times
    do k = 0, pver
       alpha(k) = 1.e-6_r8       ! about 10 days.
    end do

! Min and max values to keep things reasonable
    mxasym = 0.1_r8              ! max factor of 10 from |tau(c)| to |tau(-c)|
    mxrange= 0.001_r8            ! factor of 100 from max to min |tau(c)|
    n2min  = 1.e-8_r8            ! min value of Brunt-Vaisalla freq squared
    orohmin= 10._r8              ! min surface displacement for orographic wave drag
    orovmin=  2._r8              ! min wind speed for orographic wave drag
    taumin = 1.e-10_r8           ! min stress considered > 0
    tndmax = 500._r8 / 86400._r8    ! max permitted tendency (500 m/s/day)
    umcfac = 0.5_r8              ! max permitted reduction in u-c
    ubmc2mn= 0.01_r8             ! min value of (u-c)^2

! Determine other derived constants
    oroko2 = 0.5_r8 * kwv
    effkwv = fcrit2 * kwv
    rog    = r/g

! Determine the bounds of the background and orographic stress regions
    ktopbg  = 0
    kbotoro = pver
    do k = 0, pver
       if (hypi(k+1) .lt. 10000._r8) kbotbg  = k    ! spectrum source at 100 mb
!!$       if (hypi(k+1) .lt.  3000.) ktoporo = k
    end do
    ktoporo = 0

    if (masterproc) then
       write(iulog,*) 'KTOPBG  =',ktopbg
       write(iulog,*) 'KBOTBG  =',kbotbg
       write(iulog,*) 'KTOPORO =',ktoporo
       write(iulog,*) 'KBOTORO =',kbotoro
    end if

! Declare history variables for orgraphic term
    call addfld ('TTGWORO ','K/s     ',pver, 'A','T tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('UTGWORO ','m/s2    ',pver, 'A','U tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('VTGWORO ','m/s2    ',pver, 'A','V tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('TAUGWX  ','N/m2    ',1,    'A','Zonal gravity wave surface stress',        phys_decomp)
    call addfld ('TAUGWY  ','N/m2    ',1,    'A','Meridional gravity wave surface stress',   phys_decomp)

!!========Jinbo Xie=========
    !add new field
    call addfld ('DTAUX3_LS','m/s2',pver,'A','U tendency - ls orographic drag',phys_decomp)
    call addfld ('DTAUY3_LS','m/s2',pver,'A','V tendency - ls orographic drag',phys_decomp)
    call addfld ('DTAUX3_BL','m/s2',pver,'A','U tendency - bl orographic drag',phys_decomp)
    call addfld ('DTAUY3_BL','m/s2',pver,'A','V tendency - bl orographic drag',phys_decomp)
    call addfld ('DTAUX3_SS','m/s2',pver,'A','U tendency - ss orographic drag',phys_decomp)
    call addfld ('DTAUY3_SS','m/s2',pver,'A','V tendency - ss orographic drag',phys_decomp)
    !call addfld ('DTAUX3_FD','m/s2',pver,'A','U tendency - fd orographic drag',phys_decomp)
    !call addfld ('DTAUY3_FD','m/s2',pver,'A','V tendency - fd orographic drag',phys_decomp)
        !stress
    call addfld ('DUSFC_LS',  'N/m2', 1 , 'A', 'ls zonal oro surface stress',phys_decomp)
    call addfld ('DVSFC_LS',  'N/m2', 1 , 'A', 'ls merio oro surface stress',phys_decomp)
    call addfld ('DUSFC_BL',  'N/m2', 1 , 'A', 'bl zonal oro surface stress',phys_decomp)
    call addfld ('DVSFC_BL',  'N/m2', 1 , 'A', 'bl merio oro surface stress',phys_decomp)
    call addfld ('DUSFC_SS',  'N/m2', 1 , 'A', 'ss zonal oro surface stress',phys_decomp)
    call addfld ('DVSFC_SS',  'N/m2', 1 , 'A', 'ss merio oro surface stress',phys_decomp)
    !call addfld ('DUSFC_FD',  'N/m2', 1 , 'A', 'fd zonal oro surface stress',phys_decomp)
    !call addfld ('DVSFC_FD',  'N/m2', 1 , 'A', 'fd merio oro surface stress',phys_decomp)
!!========Jinbo Xie=========
!!add default
    call add_default('DTAUX3_LS       ',   1,' ' )
    call add_default('DTAUY3_LS       ',   1,' ' )
    call add_default('DTAUX3_BL       ',   1,' ' )
    call add_default('DTAUY3_BL       ',   1,' ' )
    call add_default('DTAUX3_SS       ',   1,' ' )
    call add_default('DTAUY3_SS       ',   1,' ' )
    !call add_default('DTAUX3_FD       ',   1,' ' )
    !call add_default('DTAUY3_fD       ',   1,' ' )
    !!==============
    call add_default ('DUSFC_LS      ',   1 , ' ')
    call add_default ('DVSFC_LS      ',   1 , ' ')
    call add_default ('DUSFC_BL      ',   1 , ' ')
    call add_default ('DVSFC_BL      ',   1 , ' ')
    call add_default ('DUSFC_SS      ',   1 , ' ')
    call add_default ('DVSFC_SS      ',   1 , ' ')
    !call add_default ('DUSFC_FD      ',   1 , ' ')
    !call add_default ('DVSFC_FD      ',   1 , ' ')
!!========Jinbo Xie=========

    call phys_getopts(history_budget_out = history_budget, history_budget_histfile_num_out = history_budget_histfile_num)
    if ( history_budget ) then
       call add_default ('TTGWORO', history_budget_histfile_num, ' ')
    end if


! Declare history variables for spectrum
    if (pgwv > 0) then
       call addfld ('TTGWSPEC','K/s     ',pver, 'A','T tendency - gravity wave spectrum',       phys_decomp)
       call addfld ('UTGWSPEC','m/s2    ',pver, 'A','U tendency - gravity wave spectrum',       phys_decomp)
       call addfld ('VTGWSPEC','m/s2    ',pver, 'A','V tendency - gravity wave spectrum',       phys_decomp)
       call add_default ('UTGWSPEC', 1, ' ')
       call add_default ('VTGWSPEC', 1, ' ')
    end if
    return
  end  subroutine gw_inti

!===============================================================================

  subroutine gw_intr (state,  sgh,  pblh,   dt, ptend, landfrac, kvh)
!=====Jinbo Xie=====
        use phys_grid, only: get_rlat_all_p
!====Jinbo Xie=====

!-----------------------------------------------------------------------
! Interface for multiple gravity wave drag parameterization.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: sgh(pcols)            ! standard deviation of orography
    real(r8), intent(in) :: pblh(pcols)           ! planetary boundary layer height
    real(r8), intent(in) :: dt                    ! time step
    real(r8), intent(in) :: landfrac(pcols)        ! Land fraction
    real(r8), intent(in) :: kvh(pcols,0:pver)     ! molecular thermal diffusivity

    type(physics_state), intent(in) :: state      ! physics state structure
    type(physics_ptend), intent(inout):: ptend    ! parameterization tendency structure

!================
!  Jinbo Xie3
!================
!input par
integer :: kpbl2d_in(pcols)
!================
!  Jinbo Xie3
!================

!============================================
!Jinbo Xie
!locally added gw and bl drag
    real(r8) :: dtaux3_ls(pcols,pver)
    real(r8) :: dtauy3_ls(pcols,pver)
    real(r8) :: dtaux3_bl(pcols,pver)
    real(r8) :: dtauy3_bl(pcols,pver)
!
    real(r8) :: dtaux3_ss(pcols,pver)
    real(r8) :: dtauy3_ss(pcols,pver)
    !real(r8) :: dtaux3_fd(pcols,pver)
    !real(r8) :: dtauy3_fd(pcols,pver)

    real(r8) :: dusfc_ls(pcols)
    real(r8) :: dvsfc_ls(pcols)
    real(r8) :: dusfc_bl(pcols)
    real(r8) :: dvsfc_bl(pcols)
!
    real(r8) :: dusfc_ss(pcols)
    real(r8) :: dvsfc_ss(pcols)
    !real(r8) :: dusfc_fd(pcols)
    !real(r8) :: dvsfc_fd(pcols)
!Jinbo Xie
!============================================


!---------------------------Local storage-------------------------------
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns

    integer :: i,k                                ! loop indexes
    integer :: kldv(pcols)                        ! top interface of low level stress divergence region
    integer :: kldvmn                             ! min value of kldv
    integer :: ksrc(pcols)                        ! index of top interface of source region
    integer :: ksrcmn                             ! min value of ksrc

    real(r8) :: ttgw(pcols,pver)                  ! temperature tendency
    real(r8) :: utgw(pcols,pver)                  ! zonal wind tendency
    real(r8) :: vtgw(pcols,pver)                  ! meridional wind tendency

    real(r8) :: ni(pcols,0:pver)                  ! interface Brunt-Vaisalla frequency
    real(r8) :: nm(pcols,pver)                    ! midpoint Brunt-Vaisalla frequency
    real(r8) :: rdpldv(pcols)                     ! 1/dp across low level divergence region
    real(r8) :: rhoi(pcols,0:pver)                ! interface density
    real(r8) :: tau(pcols,-pgwv:pgwv,0:pver)      ! wave Reynolds stress
    real(r8) :: tau0x(pcols)                      ! c=0 sfc. stress (zonal)
    real(r8) :: tau0y(pcols)                      ! c=0 sfc. stress (meridional)
    real(r8) :: ti(pcols,0:pver)                  ! interface temperature
    real(r8) :: ubi(pcols,0:pver)                 ! projection of wind at interfaces
    real(r8) :: ubm(pcols,pver)                   ! projection of wind at midpoints
    real(r8) :: xv(pcols)                         ! unit vectors of source wind (x)
    real(r8) :: yv(pcols)                         ! unit vectors of source wind (y)

!============================================

!-----------------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol
! Profiles of background state variables
    call gw_prof(lchnk, ncol, &
         state%u   , state%v   , state%t   , state%pmid   , state%pint, &
         rhoi      , ni        , ti        , nm)

!-----------------------------------------------------------------------------
! Non-orographic backgound gravity wave spectrum
!-----------------------------------------------------------------------------
    if (pgwv >0) then

! Determine the wave source for a background spectrum at ~100 mb

       call gw_bgnd (lchnk          , ncol       ,                           &
            state%u    , state%v    , state%t    , state%pmid , state%pint , &
            state%pdel , state%rpdel, state%lnpint,kldv       , kldvmn     , &
            ksrc       , ksrcmn     , rdpldv     , tau        , ubi        , &
            ubm        , xv         , yv         , PGWV       , kbotbg     )

! Solve for the drag profile
       call gw_drag_prof (lchnk     , ncol       ,                           &
            PGWV       , kbotbg     , ktopbg     , state%u    , state%v    , &
            state%t    , state%pint , state%pdel , state%rpdel, state%lnpint,&
            rhoi       , ni         , ti         , nm         , dt         , &
            kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv     , &
            tau        , ubi        , ubm        , xv         , yv         , &
            effgw_spec , utgw       , vtgw       , tau0x      , tau0y      )

! Add the momentum tendencies to the output tendency arrays
       do k = 1, pver
          do i = 1, ncol
             ptend%u(i,k) = utgw(i,k)
             ptend%v(i,k) = vtgw(i,k)
          end do
       end do

! Write output fields to history file
       call outfld ('UTGWSPEC', utgw, pcols, lchnk)
       call outfld ('VTGWSPEC', vtgw, pcols, lchnk)

! zero net tendencies if no spectrum computed
    else
       ptend%u = 0._r8
       ptend%v = 0._r8
    end if
!-----------------------------------------------------------------------------
! Orographic stationary gravity wave
!-----------------------------------------------------------------------------



!==========================================================================
!  Jinbo Xie4  Modification
!  Present moment, use the profile adjustment in the code rather than cam's
!==========================================================================
!1. Replaced the basic units with cam's states
!===========================================
        !this is for z,dz,dx,dy
        ! add surface height (surface geopotential/gravity) to convert CAM
        ! heights based on geopotential above surface into height above sea
        ! level
        !taken from %%module cospsimulator_intr
        !CAM is top to surface, which may be opposite in WRF
        !fv is same dlat,dlon, so we do it directly
        !%%needs to decide which to reverse!!!!!!!
        !ztop and zbot are already reversed, start from bottom to top
        !dz needs no reverse also
        !zmid is different calculation process, 
        !so it needs reverse if to use
        ztop(1:ncol,1:pver)=0._r8
        zbot(1:ncol,1:pver)=0._r8
        zmid(1:ncol,1:pver)=0._r8
        do k=1,pverp-1
        ! assign values from top
        ztop(1:ncol,k)=state%zi(1:ncol,pverp-k)
        ! assign values from bottom           
        zbot(1:ncol,k)=state%zi(1:ncol,pverp-k+1)
        end do
        !transform adding the pressure
        !transfer from surface to sea level
        do k=1,pver
                do i=1,ncol
                ztop(i,k)=ztop(i,k)+state%phis(i)/g
                zbot(i,k)=zbot(i,k)+state%phis(i)/g
                zmid(i,k)=state%zm(i,k)+state%phis(i)/g
                !dz is from bottom to top already for gw_drag
                dz(i,k)=ztop(i,k)-zbot(i,k)
                end do
        end do
        !=======Jinbo Xie=========================
        !get the layer index of pblh in layer
        kpbl2d_in=0_r8
        do i=1,pcols
        kpbl2d_in(i)=pblh_get_level_idx(zbot(i,:)-state%phis(i)/g,pblh(i))
        end do
        !=========================================
!================
!p3d as state%pmid
!p3di as state%pint
!Take care
!Jinbo Xie
!===========
        call get_rlat_all_p(lchnk, ncol, rlat)
        !=====================
        !      Jinbo Xie            
        !=====================
        !Initialize
        utgw=0._r8
        vtgw=0._r8
        ttgw=0._r8

        call gwdo_gsd(&
        u3d=state%u(:,pver:1:-1),v3d=state%v(:,pver:1:-1),t3d=state%t(:,pver:1:-1),&
        qv3d=state%q(:,pver:1:-1,1),p3d=state%pmid(:,pver:1:-1),p3di=state%pint(:,pver:1:-1),&
        pi3d=state%exner(:,pver:1:-1),z=zbot,&
        rublten=utgw(:,pver:1:-1),rvblten=vtgw(:,pver:1:-1),rthblten=ttgw(:,pver:1:-1),&
        dtaux3d_ls=dtaux3_ls(:,pver:1:-1),dtauy3d_ls=dtauy3_ls(:,pver:1:-1),&
        dtaux3d_bl=dtaux3_bl(:,pver:1:-1),dtauy3d_bl=dtauy3_bl(:,pver:1:-1),&
        dtaux3d_ss=dtaux3_ss(:,pver:1:-1),dtauy3d_ss=dtauy3_ss(:,pver:1:-1),&
        dusfcg_ls=dusfc_ls,dvsfcg_ls=dvsfc_ls,&
        dusfcg_bl=dusfc_bl,dvsfcg_bl=dvsfc_bl,&
        dusfcg_ss=dusfc_ss,dvsfcg_ss=dvsfc_ss,&
        xland=landfrac,br=rino_pub,&
        var2d=state%var,oc12d=state%oc,&
        oa2d=state%oadir,&
        ol2d=state%ol,&!dxy2d=state%dxydir,&
        znu=etamid(pver:1:-1),dz=dz,pblh=pblh,&
        cp=cpair,g=g,rd=r,rv=rh2o,ep1=zvir,pi=pi,bnvbg=nm(:,pver:1:-1),&
        dt=dt,dx=rearth*cos(rlat)*(2._r8*pi/256._r8),dy=rearth*(pi/(128._r8-1._r8)),&
        kpbl2d=kpbl2d_in,itimestep=0,gwd_opt=0,&
        ids=1,ide=pcols,jds=0,jde=0,kds=1,kde=pver, &
        ims=1,ime=pcols,jms=0,jme=0,kms=1,kme=pver, &
        its=1,ite=pcols,jts=0,jte=0,kts=1,kte=pver, &
        gwd_ls=1,gwd_bl=1,gwd_ss=0,gwd_fd=0 )


! z and dz all above surface and sea level, no need to add a new layer
! (just need an empty),gwd_opt(no need in my, take out 33 option))
!(itimestep just needs an empty, number of timestep,0)
!p_top       pressure top of the model (pa), set to 0
!gwd_opt is a no need
!znu         eta values on half (mass) levels, this is needed, currently set to
!midpoint eta value (hybrid),either is ok
!znw         eta values on full (w) levels , no need set to 0

!we also turn the index around, since wrf is bot-top, and cam is top-bot
!xland is only needed for small scale GWD, so not set in the moment

!Jinbo Xie
    !output the tendency profile and drag
    call outfld ('DTAUX3_LS', dtaux3_ls,  pcols, lchnk)
    call outfld ('DTAUY3_LS', dtauy3_ls,  pcols, lchnk)
    call outfld ('DTAUX3_BL', dtaux3_bl,  pcols, lchnk)
    call outfld ('DTAUY3_BL', dtauy3_bl,  pcols, lchnk)
    call outfld ('DTAUX3_SS', dtaux3_ss,  pcols, lchnk)
    call outfld ('DTAUY3_SS', dtauy3_ss,  pcols, lchnk)
    !call outfld ('DTAUX3_FD', dtaux3_fd,  pcols, lchnk)
    !call outfld ('DTAUY3_FD', dtauy3_fd,  pcols, lchnk)
    call outfld ('DUSFC_LS', dusfc_ls,  pcols, lchnk)
    call outfld ('DVSFC_LS', dvsfc_ls,  pcols, lchnk)
    call outfld ('DUSFC_BL', dusfc_bl,  pcols, lchnk)
    call outfld ('DVSFC_BL', dvsfc_bl,  pcols, lchnk)
    call outfld ('DUSFC_SS', dusfc_ss,  pcols, lchnk)
    call outfld ('DVSFC_SS', dvsfc_ss,  pcols, lchnk)
    !call outfld ('DUSFC_FD', dusfc_fd,  pcols, lchnk)
    !call outfld ('DVSFC_FD', dvsfc_fd,  pcols, lchnk)
!==========================================================================
!  Jinbo Xie4 Modification
!  Present moment, use the profile adjustment in the code rather than cam's
!==========================================================================
!print*,"Jinbo Xie utgw max min",maxval(utgw),minval(utgw)
!print*,"Jinbo Xie vtgw max min",maxval(vtgw),minval(vtgw)
!#if 0
    do k = 1, pver
       do i = 1, ncol
          ptend%u(i,k) = ptend%u(i,k) + utgw(i,k) * landfrac(i)
          ptend%v(i,k) = ptend%v(i,k) + vtgw(i,k) * landfrac(i)
          ptend%s(i,k) = -(ptend%u(i,k) * (state%u(i,k) + ptend%u(i,k)*0.5_r8*dt) &
                          +ptend%v(i,k) * (state%v(i,k) + ptend%v(i,k)*0.5_r8*dt))
          ttgw(i,k) = ptend%s(i,k) / cpair
          utgw(i,k) = ptend%u(i,k)
          vtgw(i,k) = ptend%v(i,k)
       end do
    end do
!#endif


! Set flags for nonzero tendencies, q not yet affected by gwd
    ptend%name  = "vertical diffusion"
    ptend%lq(:) = .FALSE.
    ptend%ls    = .TRUE.
    ptend%lu    = .TRUE.
    ptend%lv    = .TRUE.

        !Jinbo Xie
        tau0x=0.0_r8
        tau0y=0.0_r8
                                !Jinbo Xie for base flux
                                !do not add small-scale gwd at the moment
        tau0x=dusfc_ls+dusfc_bl!+dusfc_ss 
        tau0y=dvsfc_ls+dvsfc_bl!+dvsfc_ss
        !Jinbo Xie
!#endif



! Write output fields to history file
    call outfld ('UTGWORO', utgw,  pcols, lchnk)
    call outfld ('VTGWORO', vtgw,  pcols, lchnk)
    call outfld ('TTGWORO', ttgw,  pcols, lchnk)
    call outfld ('TAUGWX',  tau0x, pcols, lchnk)
    call outfld ('TAUGWY',  tau0y, pcols, lchnk)
    call outfld ('SGH    ', sgh,   pcols, lchnk)

    return
  end  subroutine gw_intr

!===============================================================================
  subroutine gw_prof (lchnk, ncol, u, v, t, pm, pi, rhoi, ni, ti, nm)
!-----------------------------------------------------------------------
! Compute profiles of background state quantities for the multiple
! gravity wave drag parameterization.
! 
! The parameterization is assumed to operate only where water vapor 
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: pm(pcols,pver)        ! midpoint pressures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures

    real(r8), intent(out) :: rhoi(pcols,0:pver)   ! interface density
    real(r8), intent(out) :: ni(pcols,0:pver)     ! interface Brunt-Vaisalla frequency
    real(r8), intent(out) :: ti(pcols,0:pver)     ! interface temperature
    real(r8), intent(out) :: nm(pcols,pver)       ! midpoint Brunt-Vaisalla frequency

!---------------------------Local storage-------------------------------
    integer :: i,k                                ! loop indexes

    real(r8) :: dtdp
    real(r8) :: n2                                ! Brunt-Vaisalla frequency squared

!-----------------------------------------------------------------------------
! Determine the interface densities and Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------

! The top interface values are calculated assuming an isothermal atmosphere 
! above the top level.
    k = 0
    do i = 1, ncol
       ti(i,k) = t(i,k+1)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k) = sqrt (g*g / (cpair*ti(i,k)))
    end do

! Interior points use centered differences
    do k = 1, pver-1
       do i = 1, ncol
          ti(i,k) = 0.5_r8 * (t(i,k) + t(i,k+1))
          rhoi(i,k) = pi(i,k) / (r*ti(i,k))
          dtdp = (t(i,k+1)-t(i,k)) / (pm(i,k+1)-pm(i,k))
          n2 = g*g/ti(i,k) * (1._r8/cpair - rhoi(i,k)*dtdp)
          ni(i,k) = sqrt (max (n2min, n2))
       end do
    end do

! Bottom interface uses bottom level temperature, density; next interface
! B-V frequency.
    k = pver
    do i = 1, ncol
       ti(i,k) = t(i,k)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k) = ni(i,k-1)
    end do

!-----------------------------------------------------------------------------
! Determine the midpoint Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------
    do k=1,pver
       do i=1,ncol
          nm(i,k) = 0.5_r8 * (ni(i,k-1) + ni(i,k))
       end do
    end do

    return
  end subroutine gw_prof

!===============================================================================
  subroutine gw_oro (lchnk, ncol,                &
       u, v, t, sgh, pm, pi, dpm, zm, nm, pblh,  &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv)
!-----------------------------------------------------------------------
! Orographic source for multiple gravity wave drag parameterization.
! 
! The stress is returned for a single wave with c=0, over orography.
! For points where the orographic variance is small (including ocean),
! the returned stress is zero. 
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: sgh(pcols)            ! standard deviation of orography
    real(r8), intent(in) :: pm(pcols,pver)        ! midpoint pressures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures
    real(r8), intent(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real(r8), intent(in) :: zm(pcols,pver)        ! midpoint heights
    real(r8), intent(in) :: nm(pcols,pver)        ! midpoint Brunt-Vaisalla frequency
    real(r8), intent(in) :: pblh(pcols)           ! planetary boundary layer height

    integer, intent(out) :: kldv(pcols)           ! top interface of low level stress div region
    integer, intent(out) :: kldvmn                ! min value of kldv
    integer, intent(out) :: ksrc(pcols)           ! index of top interface of source region
    integer, intent(out) :: ksrcmn                ! min value of ksrc

    real(r8), intent(out) :: rdpldv(pcols)        ! 1/dp across low level divergence region
    real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
    real(r8), intent(out) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
    real(r8), intent(out) :: ubm(pcols,pver)      ! projection of wind at midpoints
    real(r8), intent(out) :: xv(pcols)            ! unit vectors of source wind (x)
    real(r8), intent(out) :: yv(pcols)            ! unit vectors of source wind (y)

!---------------------------Local storage-------------------------------
    integer :: i,k                                ! loop indexes

    real(r8) :: pil                               ! don't we have pi somewhere?
    real(r8) :: lzsrc                             ! vertical wavelength at source
    real(r8) :: hdsp(pcols)                       ! surface streamline displacment height (2*sgh)
    real(r8) :: sghmax                            ! max orographic sdv to use
    real(r8) :: tauoro(pcols)                     ! c=0 stress from orography
    real(r8) :: zldv(pcols)                       ! top of the low level stress divergence region
    real(r8) :: nsrc(pcols)                       ! b-f frequency averaged over source region
    real(r8) :: psrc(pcols)                       ! interface pressure at top of source region
    real(r8) :: rsrc(pcols)                       ! density averaged over source region
    real(r8) :: usrc(pcols)                       ! u wind averaged over source region
    real(r8) :: vsrc(pcols)                       ! v wind averaged over source region

!---------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the apropiate
! values of wind, stability, etc. for determining the wave source are 
! averages over the depth of the atmosphere pentrated by the typical mountain.
! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
!---------------------------------------------------------------------------
    k = pver
    do i = 1, ncol
       ksrc(i) = k-1
       psrc(i) = pi(i,k-1)
       rsrc(i) = pm(i,k)/(r*t(i,k)) * dpm(i,k)
       usrc(i) = u(i,k) * dpm(i,k)
       vsrc(i) = v(i,k) * dpm(i,k)
       nsrc(i) = nm(i,k)* dpm(i,k)
       hdsp(i) = 2.0_r8 * sgh(i)
    end do
    do k = pver-1, pver/2, -1
       do i = 1, ncol
          if (hdsp(i) > sqrt(zm(i,k)*zm(i,k+1))) then
             ksrc(i) = k-1
             psrc(i) = pi(i,k-1)
             rsrc(i) = rsrc(i) + pm(i,k) / (r*t(i,k))* dpm(i,k)
             usrc(i) = usrc(i) + u(i,k) * dpm(i,k)
             vsrc(i) = vsrc(i) + v(i,k) * dpm(i,k)
             nsrc(i) = nsrc(i) + nm(i,k)* dpm(i,k)
          end if
       end do
    end do
    do i = 1, ncol
       rsrc(i) = rsrc(i) / (pi(i,pver) - psrc(i))
       usrc(i) = usrc(i) / (pi(i,pver) - psrc(i))
       vsrc(i) = vsrc(i) / (pi(i,pver) - psrc(i))
       nsrc(i) = nsrc(i) / (pi(i,pver) - psrc(i))

       if (single_column) then
! needed the following fix when winds are identically 0
! orig ->  ubi(i,pver) = sqrt (usrc(i)**2 + vsrc(i)**2)

          ubi(i,pver) = max(sqrt (usrc(i)**2 + vsrc(i)**2),orovmin)

       else
          ubi(i,pver) = sqrt (usrc(i)**2 + vsrc(i)**2)
#if (defined BFB_CAM_SCAM_IOP )
          ubi(i,pver) = max(sqrt (usrc(i)**2 + vsrc(i)**2),orovmin)
#endif
       endif
!=============================== zhh ============================================
! important
       if ( abs(ubi(i,pver)) < 1E-12) then
          xv(i) = 0.0_r8
          yv(i) = 0.0_r8
       else
          xv(i) = usrc(i) / ubi(i,pver)
          yv(i) = vsrc(i) / ubi(i,pver)
       end if
!============================ 2012.01.14 =========================================
    end do

! Project the local wind at midpoints onto the source wind.
    do k = 1, pver
       do i = 1, ncol
          ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
       end do
    end do

! Compute the interface wind projection by averaging the midpoint winds.
! Use the top level wind at the top interface.
    do i = 1, ncol
       ubi(i,0) = ubm(i,1)
    end do
    do k = 1, pver-1
       do i = 1, ncol
          ubi(i,k) = 0.5_r8 * (ubm(i,k) + ubm(i,k+1))
       end do
    end do

!---------------------------------------------------------------------------
! Determine the depth of the low level stress divergence region, as
! the max of the source region depth and 1/2 the vertical wavelength at the
! source. 
!---------------------------------------------------------------------------
    pil = acos(-1._r8)
    do i = 1, ncol
       lzsrc   = min(2._r8 * pil * usrc(i) / nsrc(i), lzmax)
       zldv(i) = max(hdsp(i), 0.5_r8 * lzsrc)
    end do

! find the index of the top of the low level divergence region

    kldv(:) = pver-1
    do k = pver-1, pver/2, -1
       do i = 1, ncol
          if (zldv(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
             kldv(i)  = k-1
          end if
       end do
    end do

! Determine the orographic c=0 source term following McFarlane (1987).
! Set the source top interface index to pver, if the orographic term is zero.
    do i = 1, ncol
       if ((ubi(i,pver) .gt. orovmin) .and. (hdsp(i) .gt. orohmin)) then
          sghmax = fcrit2 * (ubi(i,pver) / nsrc(i))**2
          tauoro(i) = oroko2 * min(hdsp(i)**2, sghmax) * rsrc(i) * nsrc(i) * ubi(i,pver)
       else
          tauoro(i) = 0._r8
          ksrc(i) = pver
          kldv(i) = pver
       end if
    end do

! Set the phase speeds and wave numbers in the direction of the source wind.
! Set the source stress magnitude (positive only, note that the sign of the 
! stress is the same as (c-u).
    do i = 1, ncol
       tau(i,0,pver) = tauoro(i)
    end do

! Determine the min value of kldv and ksrc for limiting later loops
! and the pressure at the top interface of the low level stress divergence
! region.

    ksrcmn = pver
    kldvmn = pver
    do i = 1, ncol
       ksrcmn = min(ksrcmn, ksrc(i))
       kldvmn = min(kldvmn, kldv(i))
       if (kldv(i) .ne. pver) then
          rdpldv(i) = 1._r8 / (pi(i,kldv(i)) - pi(i,pver))
       end if
    end do
    if (fracldv .le. 0._r8) kldvmn = pver

    return
  end subroutine gw_oro

!===============================================================================
  subroutine gw_bgnd (lchnk, ncol, &
       u, v, t, pm, pi, dpm, rdpm, piln,     &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv, &
       ngwv, kbot)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
! 
! The parameterization is assumed to operate only where water vapor 
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------
    use phys_grid, only: get_rlat_all_p,get_rlon_all_p
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns

    integer, intent(in) :: kbot                   ! index of bottom (source) interface
    integer, intent(in) :: ngwv                   ! number of gravity waves to use

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: pm(pcols,pver)        ! midpoint pressures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures
    real(r8), intent(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real(r8), intent(in) :: rdpm(pcols,pver)      ! 1. / (pi(k)-pi(k-1))
    real(r8), intent(in) :: piln(pcols,0:pver)    ! ln(interface pressures)

    integer, intent(out) :: kldv(pcols)           ! top interface of low level stress divergence region
    integer, intent(out) :: kldvmn                ! min value of kldv
    integer, intent(out) :: ksrc(pcols)           ! index of top interface of source region
    integer, intent(out) :: ksrcmn                ! min value of ksrc

    real(r8), intent(in) :: rdpldv(pcols)        ! 1/dp across low level divergence region
    real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
    real(r8), intent(out) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
    real(r8), intent(out) :: ubm(pcols,pver)      ! projection of wind at midpoints
    real(r8), intent(out) :: xv(pcols)            ! unit vectors of source wind (x)
    real(r8), intent(out) :: yv(pcols)            ! unit vectors of source wind (y)

!---------------------------Local storage-------------------------------
    integer :: i,k,l                              ! loop indexes

    real(r8) :: rlat(pcols)                       ! latitude in radians for columns
    real(r8) :: tauback(pcols)                    ! background stress at c=0
    real(r8) :: usrc(pcols)                       ! u wind averaged over source region
    real(r8) :: vsrc(pcols)                       ! v wind averaged over source region
    real(r8) :: al0                               ! Used in lat dependence of GW spec. 
    real(r8) :: dlat0                             ! Used in lat dependence of GW spec.
    real(r8) :: flat_gw                           ! The actual lat dependence of GW spec.
    real(r8) :: pi_g                              ! 3.14........

!---------------------------------------------------------------------------
! Determine the source layer wind and unit vectors, then project winds.
!---------------------------------------------------------------------------

! Just use the source level interface values for the source
! wind speed and direction (unit vector).
    k = kbot
    do i = 1, ncol
       ksrc(i) = k
       kldv(i) = k
       usrc(i) = 0.5_r8*(u(i,k+1)+u(i,k))
       vsrc(i) = 0.5_r8*(v(i,k+1)+v(i,k))
       ubi(i,kbot) = sqrt (usrc(i)**2 + vsrc(i)**2)
!=============================== zhh ============================================
! important, revised by Zhang He, 2013-03-08
       if ( abs(ubi(i,k)) < 1E-12) then
          xv(i) = 0.0
          yv(i) = 0.0
       else
          xv(i) = usrc(i) / ubi(i,k)
          yv(i) = vsrc(i) / ubi(i,k)
       end if
!=============================== zhh ============================================
       xv(i) = usrc(i) / ubi(i,k)
       yv(i) = vsrc(i) / ubi(i,k)
    end do

! Project the local wind at midpoints onto the source wind.
    do k = 1, kbot
       do i = 1, ncol
          ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
       end do
    end do

! Compute the interface wind projection by averaging the midpoint winds.
! Use the top level wind at the top interface.
    do i = 1, ncol
       ubi(i,0) = ubm(i,1)
    end do
    do k = 1, kbot-1
       do i = 1, ncol
          ubi(i,k) = 0.5_r8 * (ubm(i,k) + ubm(i,k+1))
       end do
    end do

!-----------------------------------------------------------------------
! Gravity wave sources
!-----------------------------------------------------------------------

! Determine the background stress at c=0
    do i=1,ncol
      tauback(i) = taubgnd * tauscal
    enddo

!	Include dependence on latitude:
! 	The lat function was obtained by RR Garcia as 
!	currently used in his 2D model
!       [Added by F. Sassi on May 30, 1997]

    pi_g=4._r8*atan(1._r8)
    al0=40._r8*pi_g/180._r8
    dlat0=10._r8*pi_g/180._r8
!
    call get_rlat_all_p(lchnk, ncol, rlat)
    do i=1,ncol
      flat_gw= 0.5_r8*(1._r8+tanh((rlat(i)-al0)/dlat0)) + 0.5_r8*(1._r8+tanh(-(rlat(i)+al0)/dlat0)) 
      flat_gw=max(flat_gw,0.2_r8)
      tauback(i)=tauback(i)*flat_gw
    enddo

! Set the phase speeds and wave numbers in the direction of the source wind.
! Set the source stress magnitude (positive only, note that the sign of the 
! stress is the same as (c-u).

    do l = 1, ngwv
       do i = 1, ncol
          tau(i, l,kbot) = tauback(i) * exp(-(c(l)/30._r8)**2)
          tau(i,-l,kbot) = tau(i, l,kbot)
       end do
    end do
    do i = 1, ncol
       tau(i,0,kbot) = tauback(i)
    end do

! Determine the min value of kldv and ksrc for limiting later loops
! and the pressure at the top interface of the low level stress divergence
! region.

    ksrcmn = pver
    kldvmn = pver

    return
  end  subroutine gw_bgnd

!===============================================================================
  subroutine gw_drag_prof (lchnk, ncol,                           &
       ngwv, kbot, ktop, u, v, t,                                 &
       pi, dpm, rdpm, piln, rhoi, ni, ti, nm, dt,                 &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv, &
       effgw, ut, vt, tau0x, tau0y)
!-----------------------------------------------------------------------
! Solve for the drag profile from the multiple gravity wave drag
! parameterization.
! 1. scan up from the wave source to determine the stress profile
! 2. scan down the stress profile to determine the tendencies
!     => apply bounds to the tendency
!          a. from wkb solution
!          b. from computational stability constraints
!     => adjust stress on interface below to reflect actual bounded tendency
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: kbot                   ! index of bottom (source) interface
    integer, intent(in) :: ktop                   ! index of top interface of gwd region
    integer, intent(in) :: ngwv                   ! number of gravity waves to use
    integer, intent(in) :: kldv(pcols)            ! top interface of low level stress  divergence region
    integer, intent(in) :: kldvmn                 ! min value of kldv
    integer, intent(in) :: ksrc(pcols)            ! index of top interface of source region
    integer, intent(in) :: ksrcmn                 ! min value of ksrc

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures
    real(r8), intent(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real(r8), intent(in) :: rdpm(pcols,pver)      ! 1. / (pi(k)-pi(k-1))
    real(r8), intent(in) :: piln(pcols,0:pver)    ! ln(interface pressures)
    real(r8), intent(in) :: rhoi(pcols,0:pver)    ! interface density
    real(r8), intent(in) :: ni(pcols,0:pver)      ! interface Brunt-Vaisalla frequency
    real(r8), intent(in) :: ti(pcols,0:pver)      ! interface temperature
    real(r8), intent(in) :: nm(pcols,pver)        ! midpoint Brunt-Vaisalla frequency
    real(r8), intent(in) :: dt                    ! time step
    real(r8), intent(in) :: rdpldv(pcols)         ! 1/dp across low level divergence region
    real(r8), intent(in) :: ubi(pcols,0:pver)     ! projection of wind at interfaces
    real(r8), intent(in) :: ubm(pcols,pver)       ! projection of wind at midpoints
    real(r8), intent(in) :: xv(pcols)             ! unit vectors of source wind (x)
    real(r8), intent(in) :: yv(pcols)             ! unit vectors of source wind (y)
    real(r8), intent(in) :: effgw                 ! tendency efficiency

    real(r8), intent(inout) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress

    real(r8), intent(out) :: ut(pcols,pver)       ! zonal wind tendency
    real(r8), intent(out) :: vt(pcols,pver)       ! meridional wind tendency
    real(r8), intent(out) :: tau0x(pcols)         ! c=0 sfc. stress (zonal)
    real(r8), intent(out) :: tau0y(pcols)         ! c=0 sfc. stress (meridional)

!---------------------------Local storage-------------------------------
    integer :: i,k,l                              ! loop indexes

    real(r8) :: d(pcols)                          ! "total" diffusivity 
    real(r8) :: dsat(pcols,-pgwv:pgwv)            ! saturation diffusivity
    real(r8) :: dscal                             ! fraction of dsat to use
    real(r8) :: mi                                ! imaginary part of vertical wavenumber
    real(r8) :: taudmp                            ! stress after damping
    real(r8) :: taumax(pcols)                     ! max(tau) for any l
    real(r8) :: tausat(pcols,-pgwv:pgwv)          ! saturation stress
    real(r8) :: ubmc(pcols,-pgwv:pgwv)            ! (ub-c)
    real(r8) :: ubmc2                             ! (ub-c)**2
    real(r8) :: ubt(pcols,pver)                   ! ubar tendency
    real(r8) :: ubtl                              ! ubar tendency from wave l
    real(r8) :: ubtlsat                           ! saturation tendency

! Initialize gravity wave drag tendencies to zero

    do k=1,pver
       do i=1,pcols
          ut(i,k) = 0._r8
          vt(i,k) = 0._r8
       end do
    end do

!---------------------------------------------------------------------------
! Compute the stress profiles and diffusivities
!---------------------------------------------------------------------------

! Loop from bottom to top to get stress profiles      

    do k = kbot-1, ktop, -1

! Determine the absolute value of the saturation stress and the diffusivity
! for each wave.
! Define critical levels where the sign of (u-c) changes between interfaces.

       d(:ncol) = dback
       do l = -ngwv, ngwv
          do i = 1, ncol
             ubmc(i,l) = ubi(i,k) - c(l)
             tausat(i,l) = abs (effkwv * rhoi(i,k) * ubmc(i,l)**3 / (2._r8*ni(i,k)) )
             if (tausat(i,l) .le. taumin) tausat(i,l) = 0.0_r8
             if (single_column) then
! needed the following fix when winds are identically 0
! ie. thermal equilibrium case
                if ((ubi(i,k+1) .ne. c(l))) then
                   if (ubmc(i,l) / (ubi(i,k+1) - c(l)) .le. 0.0_r8) tausat(i,l) = 0.0_r8
                else
                   tausat(i,l) = 0.0_r8		
                end if
             else
!--------------------------------- zhh 2013-03-12 --------------------------------
!zhh                if (ubmc(i,l) / (ubi(i,k+1) - c(l)) .le. 0.0_r8) tausat(i,l) = 0.0_r8
                if (abs(ubi(i,k+1) - c(l)) < 1E-15)  then
                   tausat(i,l) = 0.0
                else if (ubmc(i,l) / (ubi(i,k+1) - c(l)) .le. 0.0) then
                   tausat(i,l) = 0.0
                end if
!--------------------------------- zhh 2013-03-12 --------------------------------
             end if
             dsat(i,l) = (ubmc(i,l) / ni(i,k))**2 * &
                  (effkwv * ubmc(i,l)**2 / (rog * ti(i,k) * ni(i,k)) - alpha(k))
             dscal = min (1.0_r8, tau(i,l,k+1) / (tausat(i,l)+taumin))
             d(i) = max( d(i), dscal * dsat(i,l))
          end do
       end do

! Compute stress for each wave. The stress at this level is the min of 
! the saturation stress and the stress at the level below reduced by damping.
! The sign of the stress must be the same as at the level below.

       do l = -ngwv, ngwv
          do i = 1, ncol
             ubmc2 = max(ubmc(i,l)**2, ubmc2mn)
             mi = ni(i,k) / (2._r8 * kwv * ubmc2) * (alpha(k) + ni(i,k)**2/ubmc2 * d(i))
             taudmp = tau(i,l,k+1) * exp(-2._r8*mi*rog*t(i,k+1)*(piln(i,k+1)-piln(i,k)))
             if (taudmp .le. taumin) taudmp = 0._r8
             tau(i,l,k) = min (taudmp, tausat(i,l))
          end do
       end do

! The orographic stress term does not change across the source region
! Note that k ge ksrcmn cannot occur without an orographic source term

       if (k .ge. ksrcmn) then
          do i = 1, ncol
             if (k .ge. ksrc(i)) then
                tau(i,0,k) = tau(i,0,pver) 
             end if
          end do
       end if

! Require that the orographic term decrease linearly (with pressure) 
! within the low level stress divergence region. This supersedes the above
! requirment of constant stress within the source region.
! Note that k ge kldvmn cannot occur without an orographic source term, since
! kldvmn=pver then and k<=pver-1

       if (k .ge. kldvmn) then
          do i = 1, ncol
             if (k .ge. kldv(i)) then
                tau(i,0,k) = min (tau(i,0,k), tau(i,0,pver)  * &
                     (1._r8 - fracldv * (pi(i,k)-pi(i,pver)) * rdpldv(i)))
             end if
          end do
       end if

! Apply lower bounds to the stress if ngwv > 0.

       if (ngwv .ge. 1) then

! Determine the max value of tau for any l

          do i = 1, ncol
             taumax(i) = tau(i,-ngwv,k)
          end do
          do l = -ngwv+1, ngwv
             do i = 1, ncol
                taumax(i) = max(taumax(i), tau(i,l,k))
             end do
          end do
          do i = 1, ncol
             taumax(i) = mxrange * taumax(i)
          end do

! Set the min value of tau for each wave to the max of mxrange*taumax or
! mxasym*tau(-c)

          do l = 1, ngwv
             do i = 1, ncol
                tau(i, l,k) = max(tau(i, l,k), taumax(i))
                tau(i, l,k) = max(tau(i, l,k), mxasym*tau(i,-l,k))
                tau(i,-l,k) = max(tau(i,-l,k), taumax(i))
                tau(i,-l,k) = max(tau(i,-l,k), mxasym*tau(i, l,k))
             end do
          end do
          l = 0
          do i = 1, ncol
             tau(i,l,k) = max(tau(i,l,k), mxasym * 0.5_r8 * (tau(i,l-1,k) + tau(i,l+1,k)))
          end do

       end if

    end do

! Put an upper bound on the stress at the top interface to force some stress
! divergence in the top layer. This prevents the easterlies from running
! away in the summer mesosphere, since most of the gravity wave activity
! will pass through a top interface at 75--80 km under strong easterly
! conditions. 
!++BAB fix to match ccm3.10
!!$    do l = -ngwv, ngwv
!!$       do i = 1, ncol
!!$          tau(i,l,0) = min(tau(i,l,0), 0.5*tau(i,l,1))
!!$       end do
!!$    end do
!--BAB fix to match ccm3.10

!---------------------------------------------------------------------------
! Compute the tendencies from the stress divergence.
!---------------------------------------------------------------------------

! Loop over levels from top to bottom
    do k = ktop+1, kbot

! Accumulate the mean wind tendency over wavenumber.
       do i = 1, ncol
          ubt (i,k) = 0.0_r8
       end do
       do l = -ngwv, ngwv
          do i = 1, ncol

! Determine the wind tendency including excess stress carried down from above.
             ubtl = g * (tau(i,l,k)-tau(i,l,k-1)) * rdpm(i,k)

! Require that the tendency be no larger than the analytic solution for
! a saturated region [proportional to (u-c)^3].
             ubtlsat = effkwv * abs((c(l)-ubm(i,k))**3.0_r8) /(2._r8*rog*t(i,k)*nm(i,k))
             ubtl = min(ubtl, ubtlsat)

! Apply tendency limits to maintain numerical stability.
! 1. du/dt < |c-u|/dt  so u-c cannot change sign (u^n+1 = u^n + du/dt * dt)
! 2. du/dt < tndmax    so that ridicuously large tendency are not permitted
             ubtl = min(ubtl, umcfac * abs(c(l)-ubm(i,k)) / dt)
             ubtl = min(ubtl, tndmax)

! Accumulate the mean wind tendency over wavenumber.
             ubt (i,k) = ubt (i,k) + sign(ubtl, c(l)-ubm(i,k))

! Redetermine the effective stress on the interface below from the wind 
! tendency. If the wind tendency was limited above, then the new stress
! will be small than the old stress and will cause stress divergence in
! the next layer down. This has the effect of smoothing large stress 
! divergences downward while conserving total stress.
             tau(i,l,k) = tau(i,l,k-1) + ubtl * dpm(i,k) / g
          end do
       end do

! Project the mean wind tendency onto the components and scale by "efficiency".
       do i = 1, ncol
          ut(i,k) = ubt(i,k) * xv(i) * effgw
          vt(i,k) = ubt(i,k) * yv(i) * effgw
       end do

! End of level loop
    end do

!-----------------------------------------------------------------------
! Project the c=0 stress (scaled) in the direction of the source wind for recording
! on the output file.
!-----------------------------------------------------------------------
    do i = 1, ncol
       tau0x(i) = tau(i,0,kbot) * xv(i) * effgw
       tau0y(i) = tau(i,0,kbot) * yv(i) * effgw
    end do

    return
  end subroutine gw_drag_prof


!===========Jinbo Xie===============================
function pblh_get_level_idx(height_array ,pblheight)
implicit none
real(8),intent(in),dimension(30) :: height_array
real(8),intent(in) :: pblheight
integer :: pblh_get_level_idx

!local
integer :: i
logical :: found

pblh_get_level_idx = -1
found=.False.

do i = 1, pver
        if((pblheight >= height_array(i).and.pblheight <height_array(i+1)))then
                pblh_get_level_idx =  i
                found=.True.
                return
        endif
enddo
end function
!================================Jinbo Xie====================
!+czy20181120
!================================Jinbo Xie====================
!-------------------------------------------------------------------------------
   subroutine gwdo_gsd(u3d,v3d,t3d,qv3d,p3d,p3di,pi3d,z,                       &
                  rublten,rvblten,rthblten,                                    &
                  dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,                 &
                  dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd,                 &
                  dusfcg_ls,dvsfcg_ls,dusfcg_bl,dvsfcg_bl,dusfcg_ss,dvsfcg_ss, &
                  dusfcg_fd,dvsfcg_fd,xland,br,                                &
                  var2d,oc12d,oa2d,ol2d,znu,znw,p_top,dz,pblh,          &
                  cp,g,rd,rv,ep1,pi,bnvbg,                                           &                        
                  dt,dx,dy,kpbl2d,itimestep,gwd_opt,                           &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte,                                 &
                    gwd_ls,gwd_bl,gwd_ss,gwd_fd)!Jinbo Xie added
!Jinbo Xie add dy, since global model is not dx=dy
!===========================
! Jinbo Xie0
!===========================
        !dt,dx,kpbl2d,itimestep,gwd_opt,
        !&
    !                ids,ide, jds,jde, kds,kde,                                 &
    !                ims,ime, jms,jme, kms,kme,                                 &
    !                its,ite, jts,jte, kts,kte))
!=================================
!   Jinbo Xie0 Modification
!changed the index used in WRF
!to chunk in CAM
!==============================
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!                                                                       
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- p3d         3d pressure (pa)
!-- p3di        3d pressure (pa) at interface level
!-- pi3d        3d exner function (dimensionless)
!-- rublten     u tendency due to pbl parameterization (m/s/s) 
!-- rvblten     v tendency due to pbl parameterization (m/s/s)
!-- rthblten    theta tendency due to pbl parameterization (K/s)
!-- znu         eta values (sigma values)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rd          gas constant for dry air (j/kg/k)
!-- z           height above sea level (m)
!-- rv          gas constant for water vapor (j/kg/k)
!-- dt          time step (s)
!-- dx          model grid interval (m)
!-- dz          height of model layers (m)
!-- xland       land mask (1 for land, 2 for water)
!-- br          bulk richardson number in surface layer
!-- pblh        planetary boundary layer height (m)
!-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!
!-------------------------------------------------------------------------------
  integer,  intent(in   )   ::      ids,ide, jds,jde, kds,kde,                 &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte
  integer,  intent(in   )   ::      itimestep,gwd_opt
!

  !real(r8),     intent(in   )   ::      dt,cp,g,rd,rv,ep1,pi!dt,dx,cp,g,rd,rv,ep1,pi
        real(r8),     intent(in   )   ::      cp,g,rd,rv,ep1,pi
!======Jinbo Xie=========
        real(r8),     intent(in), optional   ::  dt
        real(r8),     intent(in), dimension( ims:ime, kms:kme ),optional   ::  bnvbg
!======Jinbo Xie=========
!==========================
!Jinbo Xie add dy
        real(r8),     intent(in)   ::      dy
!variable
        !real(r8),   !dimension(:),  intent(in   )   ::      dx
        real(r8),    intent(in)   ::      dx(:)
!==========================
!
  real(r8),     dimension( ims:ime, kms:kme )             ,  &
        
            intent(in   )   ::                      qv3d, &
                                                              p3d, &
                                                             pi3d, &
                                                              t3d, &
                                                                z, &
                                                               dz
  real(r8),     dimension( ims:ime, kms:kme )                    , &
     intent(in   )   ::                                 p3di
  real(r8),     dimension( ims:ime, kms:kme )                    , &
           intent(inout)   ::                   rublten, &
                                                          rvblten, &
                                                          rthblten
  real(r8),     dimension( ims:ime, kms:kme ), optional                 , &
            intent(inout)   ::  dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,   &
                                dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd
!
  real(r8),     dimension( ims:ime, kms:kme)   ::                                    &
                                  dtaux2d_ls,dtauy2d_ls,dtaux2d_bl,dtauy2d_bl, &
                                  dtaux2d_ss,dtauy2d_ss,dtaux2d_fd,dtauy2d_fd

  real(r8),      dimension( ims:ime, kms:kme )                          , &
        
             intent(in   )   ::                        u3d, &
                                                                v3d
!
  integer,   dimension( ims:ime )                                   , &
             intent(in  )   ::             kpbl2d

  real(r8),   dimension( ims:ime )                                      , &
        intent(in  )   ::                                   pblh, &
                                                                 br, &
                                                                 xland

  real(r8),   dimension( ims:ime ), optional                            , &
             intent(inout  )   ::  dusfcg_ls,dvsfcg_ls,dusfcg_bl,dvsfcg_bl,    &
                                   dusfcg_ss,dvsfcg_ss,dusfcg_fd,dvsfcg_fd

  real(r8),   dimension( ims:ime ) ::  dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,        &
                                       dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd
!
!Jinbo Xie
!optional because other schemes may not need these par
!Jinbo Xie
  !real(r8),   dimension( ims:ime ),                                     , &
        real(r8),   dimension( ims:ime ),     optional                   , &
             intent(in  )   ::                                  var2d, &
                                                                oc12d

       !real(r8),dimension(ims:ime,nvar_dirOL),optional, intent(in) :: dxy2d
#ifndef continuous
       real(r8),dimension(ims:ime,nvar_dirOL),optional, intent(in) :: ol2d
#else
        real(r8),dimension(4,ims:ime,indexb),optional, intent(in) :: ol2d
#endif
        real(r8),dimension(ims:ime,nvar_dirOA),optional, intent(in) :: oa2d
!Jinbo Xie



!
  !real(r8),     dimension( kms:kme )                                             , &
            real(r8),   optional                                                         , &
            intent(in  )   ::                                             znu(:), &
                                                                          znw(:)
!
  real(r8),     optional, intent(in  )   ::                           p_top
!
!local
!
  real(r8),   dimension( its:ite, kts:kte )  ::                           delprsi, &
                                                                          pdh
  real(r8),   dimension( its:ite, kts:kte+1 )   ::     pdhi
  real(r8),   dimension( its:ite, nvar_dirOA )  ::     oa4
  !real(r8),   dimension( its:ite, nvar_dirOL )  ::     dxy42d
#ifndef continuous
  real(r8),   dimension( its:ite, nvar_dirOL )  ::     ol4
#else
  real(r8),   dimension( 4,its:ite, indexb )    ::     ol4
#endif
  integer ::  i,j,k,kdt,kpblmax
!
        !===========Jinbo Xie=========
        integer , intent(in) :: gwd_ls,gwd_bl,gwd_ss,gwd_fd!Jinbo Xie 
!real(r8),dimension(ims:ime, kms:kme) :: temp
        !===========Jinbo Xie=========


   do k = kts,kte
     if(znu(k).gt.0.6_r8) kpblmax = k + 1
   enddo
!
      do k = kts,kte+1
         do i = its,ite
            if(k.le.kte)pdh(i,k) = p3d(i,k)
             pdhi(i,k) = p3di(i,k)
         enddo
      enddo
!
      do k = kts,kte
        do i = its,ite
          delprsi(i,k) = pdhi(i,k)-pdhi(i,k+1)
        enddo
      enddo


!Jinbo Xie
!temp=t3d
!Jinbo Xie
!=======Jinbo Xie=================
!no need when there is no large drag
IF ( (gwd_ls .EQ. 1).and.(gwd_bl .EQ. 1)) then

        do i = its,ite
            oa4(i,:) = oa2d(i,:)
#ifndef continuous
            ol4(i,:) = ol2d(i,:)
            !dxy42d(i,:)=dxy2d(i,:)!rearth*(2*pi/360.)*dxy2d(i,:)
#else
            ol4(:,i,:) = ol2d(:,i,:)
            !turn terrx and terry into (m) unit
            !to be aligned with dxmter and dymeter (m) 
            ol4(2:3,i,:)=rearth*ol4(2:3,i,:)
#endif
        enddo
ENDIF
!========Jinbo Xie================

	!=================================================================
	! Jinbo Xie1  changed all ,j out, turn 3d into 2d, for cam formation
	!=================================================================
      call gwdo2d(dudt=rublten(ims,kms),dvdt=rvblten(ims,kms)                  &
             ,dthdt=rthblten(ims,kms)                                          &
              ,dtaux2d_ls=dtaux2d_ls,dtauy2d_ls=dtauy2d_ls                     &
              ,dtaux2d_bl=dtaux2d_bl,dtauy2d_bl=dtauy2d_bl                     &
              ,dtaux2d_ss=dtaux2d_ss,dtauy2d_ss=dtauy2d_ss                     &
              ,dtaux2d_fd=dtaux2d_fd,dtauy2d_fd=dtauy2d_fd                     &
              ,u1=u3d(ims,kms),v1=v3d(ims,kms)                                 &
              ,t1=t3d(ims,kms)                                                &
              ,q1=qv3d(ims,kms)                                                 &
              ,del=delprsi(its,kts)                                            &
              ,prsi=pdhi(its,kts)                                              &
              ,prsl=pdh(its,kts),prslk=pi3d(ims,kms)                           &
              ,zl=z(ims,kms),rcl=1.0_r8                                        &
              ,xland1=xland(ims),br1=br(ims),hpbl=pblh(ims)                    &
              ,bnv_in=bnvbg(ims,kms)                                         &
              ,dz2=dz(ims,kms)                                                 &
              ,kpblmax=kpblmax                                                 &
              ,dusfc_ls=dusfc_ls,dvsfc_ls=dvsfc_ls                             &
              ,dusfc_bl=dusfc_bl,dvsfc_bl=dvsfc_bl                             &
              ,dusfc_ss=dusfc_ss,dvsfc_ss=dvsfc_ss                             &
              ,dusfc_fd=dusfc_fd,dvsfc_fd=dvsfc_fd                             &
              ,var=var2d(ims),oc1=oc12d(ims)                                   &
              ,oa4=oa4,ol4=ol4                                                 &
              ,g=g,cp=cp,rd=rd,rv=rv,fv=ep1,pi=pi                              &
              ,dxmeter=dx,dymeter=dy,deltim=dt                                 &
              ,kpbl=kpbl2d(ims),kdt=itimestep,lat=j                            &
              ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde               &
              ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme               &
              ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte               &
              ,gsd_gwd_ls=gwd_ls,gsd_gwd_bl=gwd_bl,gsd_gwd_ss=gwd_ss,gsd_gwd_fd=gwd_fd)
                !=============Jinbo Xie 



!#if 0
!Jinbo Xie
!added dymeter in here
        !IF (gwd_opt == 33) then !research mode
!#if 1
!1
IF ( (gwd_ls .EQ. 1).and.(gwd_bl .EQ. 1)) then
                do i = its,ite
                dusfcg_ls(i)=dusfc_ls(i)
                dvsfcg_ls(i)=dvsfc_ls(i)
                dusfcg_bl(i)=dusfc_bl(i)
                dvsfcg_bl(i)=dvsfc_bl(i)
                enddo
        
             dtaux3d_ls=dtaux2d_ls
             dtaux3d_bl=dtaux2d_bl
             dtauy3d_ls=dtauy2d_ls
             dtauy3d_bl=dtauy2d_bl

!ENDIF
!Jinbo Xie temporarily for base flux
!IF (gwd_ss .EQ. 1) then
                do i = its,ite
                dusfcg_ss(i)=dusfc_ss(i)
                dvsfcg_ss(i)=dvsfc_ss(i)
                end do

                dtaux3d_ss=dtaux2d_ss
                dtauy3d_ss=dtauy2d_ss
!ENDIF
ENDIF
       
IF (gwd_fd .EQ. 1) then

                do i = its,ite
                dusfcg_fd(i)=dusfc_fd(i)
                dvsfcg_fd(i)=dvsfc_fd(i)
                enddo
                dtaux3d_fd=dtaux2d_fd
                dtauy3d_fd=dtauy2d_fd
ENDIF
!#endif
!
   end subroutine gwdo_gsd
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine gwdo2d(dudt,dvdt,dthdt,dtaux2d_ls,dtauy2d_ls,                    &
                    dtaux2d_bl,dtauy2d_bl,dtaux2d_ss,dtauy2d_ss,               &
                    dtaux2d_fd,dtauy2d_fd,u1,v1,t1,q1,                         &
                    del,                                                       &
                    prsi,prsl,prslk,zl,rcl,                                    &
                    xland1,br1,hpbl,bnv_in,dz2,                               &
                    kpblmax,dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,               &
                    dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,var,oc1,oa4,ol4,&
                    g,cp,rd,rv,fv,pi,dxmeter,dymeter,deltim,kpbl,kdt,lat,      &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte,                                 &
                    gsd_gwd_ls,gsd_gwd_bl,gsd_gwd_ss,gsd_gwd_fd)!Jinbo Xie 
!===============================
! Jinbo Xie add another dymeter
!===============================
!=====Jinbo Xie=====
use sub_xjb,only:OLgrid,dxygrid
!=====Jinbo Xie=====

!-------------------------------------------------------------------------------
!  
!  this code handles the time tendencies of u v due to the effect of mountain 
!  induced gravity wave drag from sub-grid scale orography. this routine 
!  not only treats the traditional upper-level wave breaking due to mountain 
!  variance (alpert 1988), but also the enhanced lower-tropospheric wave 
!  breaking due to mountain convexity and asymmetry (kim and arakawa 1995). 
!  thus, in addition to the terrain height data in a model grid gox, 
!  additional 10-2d topographic statistics files are needed, including 
!  orographic standard  deviation (var), convexity (oc1), asymmetry (oa4) 
!  and ol (ol4). these data sets are prepared based on the 30 sec usgs orography
!  hong (1999). the current scheme was implmented as in hong et al.(2008)
!
!  coded by song-you hong and young-joon kim and implemented by song-you hong
!
!  The Code is further added with flow-blocking (Kim et al. 2005), turbulent orographic form drag
!  (TOFD) (Beljaars et al. 2004), and small-scale drag (Tsiringakis et al.
!  2017).
!  program history log:
!    2018-11-15  Jinbo Xie implmented into CAS-ESM
!
!           Activation of each component is done by specifying the integer-parameters
!           (defined below) to 0: inactive or 1: active
!                    gsd_gwd_ls = 0 or 1: large-scale
!                    gsd_gwd_bl = 0 or 1: blocking drag 
!                    gsd_gwd_ss = 0 or 1: small-scale gravity wave drag
!                    gsd_gwd_fd = 0 or 1: topographic form drag
!
!  references:
!        hong et al. (2008), wea. and forecasting
!        kim and doyle (2005), Q. J. R. Meteor. Soc.
!        kim and arakawa (1995), j. atmos. sci.
!        alpet et al. (1988), NWP conference.
!        hong (1999), NCEP office note 424.
!        steeneveld et al (2008), JAMC
!        Tsiringakis et al. (2017), Q. J. R. Meteor. Soc.
!
!  notice : comparible or lower resolution orography files than model resolution
!           are desirable in preprocess (wps) to prevent weakening of the drag
!-------------------------------------------------------------------------------
!
!  input                                                                
!        dudt (ims:ime,kms:kme)  non-lin tendency for u wind component
!        dvdt (ims:ime,kms:kme)  non-lin tendency for v wind component
!        u1(ims:ime,kms:kme) zonal wind / sqrt(rcl)  m/sec  at t0-dt
!        v1(ims:ime,kms:kme) meridional wind / sqrt(rcl) m/sec at t0-dt
!        t1(ims:ime,kms:kme) temperature deg k at t0-dt
!        q1(ims:ime,kms:kme) specific humidity at t0-dt
!
!        rcl     a scaling factor = reciprocal of square of cos(lat)
!                for gmp.  rcl=1 if u1 and v1 are wind components.
!        deltim  time step    secs                                       
!        del(kts:kte)  positive increment of pressure across layer (pa)
!                                                                       
!  output
!        dudt, dvdt    wind tendency due to gwdo
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  kdt,lat,latd,lond,kpblmax,                         &
                            ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
!
   real(r8)                 ::  g,rd,rv,fv,cp,pi,deltim,rcl!dxmeter,deltim,rcl
!=======================
!!!!Jinbo Xie, add dymeter
        real(r8)                 :: dymeter
        real(r8),dimension(:)    :: dxmeter
!======================
   real(r8)                ::  dudt(ims:ime,kms:kme),dvdt(ims:ime,kms:kme),          &
                                dthdt(ims:ime,kms:kme),&
                            dtaux2d_ls(ims:ime,kms:kme),dtauy2d_ls(ims:ime,kms:kme), &
                            dtaux2d_bl(ims:ime,kms:kme),dtauy2d_bl(ims:ime,kms:kme), &
                            dtaux2d_ss(ims:ime,kms:kme),dtauy2d_ss(ims:ime,kms:kme), &
                            dtaux2d_fd(ims:ime,kms:kme),dtauy2d_fd(ims:ime,kms:kme), &
                            u1(ims:ime,kms:kme),v1(ims:ime,kms:kme),           & 
                            t1(ims:ime,kms:kme),q1(ims:ime,kms:kme),           &
                            zl(ims:ime,kms:kme),prsl(its:ite,kts:kte),         &
                            prslk(ims:ime,kms:kme)
   real(r8),intent(in)                ::  prsi(its:ite,kts:kte+1),del(its:ite,kts:kte)
!=========Jinbo Xie=======
   real(r8),intent(in),optional    ::  oa4(its:ite,nvar_dirOA)
   !real(r8),intent(in),optional    ::  dxy_in(its:ite,nvar_dirOL)
#ifndef continuous
   real(r8),intent(in),optional    ::  ol4(its:ite,nvar_dirOL)
#else
   real(r8),intent(in),optional    ::  ol4(4,its:ite,indexb)
   integer :: nind_non!number of nonfillvalue variables
#endif
!=========Jinbo Xie=======



!
! GSD surface drag options to regulate specific components
! Each component is tapered off automatically as a function of dx, so best to 
! keep them activated (=1).
!integer, parameter ::                                                          &
!   gsd_gwd_ls      = 1,       & ! large-scale gravity wave drag
!   gsd_gwd_bl      = 1,       & ! blocking drag 
!   gsd_gwd_ss      = 1,       & ! small-scale gravity wave drag (Steeneveld et al. 2008)
!   gsd_gwd_fd      = 0,       & ! form drag (Beljaars et al. 2004, QJRMS)

!
! added for small-scale orographic wave drag
   real(r8), dimension(its:ite,kts:kte)     :: utendwave,vtendwave,thx,thvx,za
   real(r8), dimension(ims:ime), intent(in) :: br1,hpbl,xland1
   real(r8), dimension(its:ite)             :: govrth
   real(r8), dimension(ims:ime,kms:kme), intent(in) :: dz2
   real(r8), dimension(its:ite,kts:kte+1)   :: zq
   real(r8)                 :: tauwavex0,tauwavey0,XNBV,density,tvcon,hpbl2
   integer              :: kpbl2,kvar
   real(r8), parameter      :: varmax = 200._r8
!
   integer              ::  kpbl(ims:ime)
   real(r8)                 ::  var(ims:ime),oc1(ims:ime),                         &
                            dusfc_ls(ims:ime),dvsfc_ls(ims:ime),               &
                            dusfc_bl(ims:ime),dvsfc_bl(ims:ime),               &
                            dusfc_ss(ims:ime),dvsfc_ss(ims:ime),               &
                            dusfc_fd(ims:ime),dvsfc_fd(ims:ime)
! Variables for scale-awareness:
! Small-scale GWD + turbulent form drag
   real(r8), parameter   :: dxmin_ss = 1000._r8, dxmax_ss = 12000._r8  ! min,max range of tapering (m)
! Large-scale GWD
   real(r8), parameter   :: dxmin_ls = 3000._r8, dxmax_ls = 13000._r8  ! min,max range of tapering (m)
!===========================
!!!!!Jinbo Xie
!!!!!Add y axis for taper consider
   real(r8), parameter   :: dymin_ls = 3000._r8, dymax_ls = 13000._r8  ! min,maxrange of tapering (m)
   real(r8), parameter   :: dymin_ss = 3000._r8, dymax_ss = 13000._r8  ! min,maxrange of tapering (m)
!==========================
   real(r8)              :: ss_taper, ls_taper  ! small- and large-scale tapering factors (-)


!
! added Beljaars orographic form drag
   real(r8), dimension(its:ite,kts:kte)     :: utendform,vtendform
   real(r8)                 :: a1,a2,wsp
! critical richardson number for wave breaking : ! larger drag with larger value
!


!==========Jinbo Xie=========
   !real(r8),parameter       ::  ric     = 0.25_r8 
real(r8),parameter       ::  ric     = 1._r8
real(r8),parameter       ::  ric_rig  = 0.25_r8
!==========Jinbo Xie==========

   real(r8),parameter       ::  dw2min  = 1._r8
   real(r8),parameter       ::  rimin   = -100._r8
   real(r8),parameter       ::  bnv2min = 1.0e-5_r8
   real(r8),parameter       ::  efmin   = 0.0_r8
   real(r8),parameter       ::  efmax   = 10.0_r8
   real(r8),parameter       ::  xl      = 4.0e4_r8
   real(r8),parameter       ::  critac  = 1.0e-5_r8
   real(r8),parameter       ::  gmax    = 1._r8
   real(r8),parameter       ::  veleps  = 1.0_r8                                              
   real(r8),parameter       ::  factop  = 0.5_r8                                               
   real(r8),parameter       ::  frc     = 1.0_r8   
   real(r8),parameter       ::  ce      = 0.8_r8  
   real(r8),parameter       ::  cg      = 0.5_r8
   integer,parameter    ::  kpblmin = 2
!
!  local variables
!
   integer              ::  j,i,k,lcap,lcapp1,nwd,idir,                          &
                            klcap,kp1,ikount,kk,nwd1!added nwd1 Jinbo Xie
!
   real(r8)                 ::  rcs,rclcs,csg,fdir,cleff,cs,rcsks,                 &
                            wdir,ti,rdz,temp,tem2,dw2,shr2,bvf2,rdelks,        &
                            wtkbj,tem,gfobnv,hd,fro,rim,temc,tem1,efact,       &
                            temv,dtaux,dtauy,eng0,eng1,theta,rad,wdir1!Jinbo Xie added theta,rad,wdir1
!=====Jinbo Xie=====
real(r8),dimension(its:ite,kts:kte),intent(in), optional :: bnv_in
!=====Jinbo Xie=====

!
   logical              ::  ldrag(its:ite),icrilv(its:ite),                    &
                            flag(its:ite),kloop1(its:ite)
!                                                                       
   real(r8)                 ::  taub(its:ite),taup(its:ite,kts:kte+1),             &
                            xn(its:ite),yn(its:ite),                           &
                            ubar(its:ite),vbar(its:ite),                       &
                            fr(its:ite),ulow(its:ite),                         &
                            rulow(its:ite),bnv(its:ite),                       &
                            oa1(its:ite),ol(its:ite),                          &
                            roll(its:ite),dtfac(its:ite),                      &
                            brvf(its:ite),xlinv(its:ite),                      &
                            delks(its:ite),delks1(its:ite),                    &
                            bnv2(its:ite,kts:kte),usqj(its:ite,kts:kte),       &
                            taud_ls(its:ite,kts:kte),taud_bl(its:ite,kts:kte), &
                            ro(its:ite,kts:kte),                               &
                            vtk(its:ite,kts:kte),vtj(its:ite,kts:kte),         &
                            zlowtop(its:ite),velco(its:ite,kts:kte-1),         &
                            coefm(its:ite)
!
   integer              ::  kbl(its:ite),klowtop(its:ite)
!
   logical :: iope
   integer,parameter    ::  mdir=2*nvar_dirOL!the number of directions

!integer              ::  nwdir(mdir)
!data nwdir/6,7,5,8,2,3,1,4/
!
!  variables for flow-blocking drag
!
   real(r8),parameter       :: frmax  = 10._r8
   real(r8),parameter       :: olmin  = 1.0e-5_r8
   real(r8),parameter       :: odmin  = 0.1_r8
   real(r8),parameter       :: odmax  = 10._r8
   real(r8),parameter       :: erad   = 6371.315e+3_r8
   integer              :: komax(its:ite)
   integer              :: kblk
   real(r8)                 :: cd
   real(r8)                 :: zblk,tautem

!Jinbo Xie
real(r8) :: zblk_col(its:ite)
real(r8) :: taub_xjb(its:ite)
real(r8) :: taufb_xjb(its:ite)
real(r8) :: wdir1_xjb(its:ite)
!Jinbo Xie


   real(r8)                 :: pe,ke 
!================================
   real(r8)                 :: dely,dxy4(its:ite,nvar_dirOL),&
                     delx(its:ite),dxy4p(its:ite,nvar_dirOL)
!=====Jinbo Xie==================
   real(r8)                 :: dxy(its:ite),dxyp(its:ite)
   real(r8)                 ::  olp(its:ite),&
                                 od(its:ite)
   real(r8)                 :: taufb(its:ite,kts:kte+1)
        !===========Jinbo Xie=========
        integer , intent(in) :: gsd_gwd_ls,gsd_gwd_bl,gsd_gwd_ss,gsd_gwd_fd        !Jinbo Xie 
        !===========Jinbo Xie=========


!=====Jinbo Xie=====
integer :: wdir_add_xjb(mdir,its:ite)
!=====Jinbo Xie=====


 !===================================
 !Jinbo Xie readdata
 !real(r8),allocatable :: need3(:,:)
 !real(r8) :: oc11(ims:ime)
 !real(r8) :: wind_xjb(kts:kte)
 !real(r8) :: shr2_xjb(its:ite,kts:kte)
 real(r8) :: l1,l2,S!,shrrok1,shrrok0,gamma1
 !
 logical  :: iint
 real(r8) :: zl_hint(its:ite)
 !===================================





!
!---- constants                                                         
!                                                                       
   rcs    = sqrt(rcl) 
   cs     = 1._r8 / sqrt(rcl)                                                     
   csg    = cs * g                                                      
   lcap   = kte                                                         
   lcapp1 = lcap + 1                                                 
   fdir   = mdir / (2.0_r8*pi)
!
!--- calculate scale-aware tapering factors
!
!=========================================
!!Jinbo Xie  add criteria for dymeter
!Taper for small GWD only, currently assumes equal length in both direction
!Taper matters not much
#if 0
if ( dxmeter .ge. dxmax_ls .and. dymeter .ge. dymax_ls) then
!=========================================
   ls_taper = 1.
else
   if ( dxmeter .le. dxmin_ls) then
      ls_taper = 0.
   else
      ls_taper = 0.5 * ( SIN(pi*(dxmeter-0.5*(dxmax_ls+dxmin_ls))/    &
                                (dxmax_ls-dxmin_ls)) + 1. )
   end if
end if
if ( dxmeter .ge. dxmax_ss ) then
   ss_taper = 1.
else
   if ( dxmeter .le. dxmin_ss) then
      ss_taper = 0.
   else
      ss_taper = dxmax_ss * (1. - dxmin_ss/dxmeter)/(dxmax_ss-dxmin_ss)
   end if
end if
#endif
!Jinbo Xie currently use smoothed topo, taper sets to none-taper
!Jinbo Xie maybe only used when using directly derived data from 30s
ls_taper=1._r8
ss_taper=1._r8

!
!--- calculate length of grid for flow-blocking drag
!
   delx   = dxmeter
!============
! Jinbo Xie2
!============
   dely=dymeter !Jinbo Xie, add dy, since global model dx/=dy
!============
! Jinbo Xie2
!============
!Jinbo Xie 
!varied delx,so everything needs add another dim
!
!
!-----initialize arrays                                                 
!                                                                       
   dtaux = 0.0_r8
   dtauy = 0.0_r8
   do i = its,ite                                                       
     klowtop(i)    = 0
     kbl(i)        = 0
   enddo                                                             
!
   do i = its,ite                                                       
     xn(i)         = 0.0_r8
     yn(i)         = 0.0_r8
     ubar (i)      = 0.0_r8
     vbar (i)      = 0.0_r8
     roll (i)      = 0.0_r8
     taub (i)      = 0.0_r8
     oa1(i)        = 0.0_r8
     ol(i)         = 0.0_r8
     ulow (i)      = 0.0_r8
     dtfac(i)      = 1.0_r8
     ldrag(i)      = .false.
     icrilv(i)     = .false. 
     flag(i)       = .true.
   enddo                                                             
!

   do k = kts,kte
     do i = its,ite
       usqj(i,k) = 0.0_r8
       bnv2(i,k) = 0.0_r8
       vtj(i,k)  = 0.0_r8
       vtk(i,k)  = 0.0_r8
       taup(i,k) = 0.0_r8
       taud_ls(i,k) = 0.0_r8
       taud_bl(i,k) = 0.0_r8
       dtaux2d_ls(i,k)= 0.0_r8
       dtauy2d_ls(i,k)= 0.0_r8
       dtaux2d_bl(i,k)= 0.0_r8
       dtauy2d_bl(i,k)= 0.0_r8
       dtaux2d_ss(i,k)= 0.0_r8
       dtauy2d_ss(i,k)= 0.0_r8
       dtaux2d_fd(i,k)= 0.0_r8
       dtauy2d_fd(i,k)= 0.0_r8
     enddo
   enddo
!

   do i = its,ite
     dusfc_ls(i) = 0.0_r8
     dvsfc_ls(i) = 0.0_r8
     dusfc_bl(i) = 0.0_r8
     dvsfc_bl(i) = 0.0_r8
     dusfc_ss(i) = 0.0_r8
     dvsfc_ss(i) = 0.0_r8
     dusfc_fd(i) = 0.0_r8
     dvsfc_fd(i) = 0.0_r8

!temp for base flux xjb
!taub_xjb(i)=0.0_r8
!taufb_xjb(i)=0.0_r8
!zblk_col(i)=0.0_r8

wdir1_xjb(i)=0.0_r8
   enddo
!
   do i = its,ite
     taup(i,kte+1) = 0.0_r8
     xlinv(i)     = 1.0_r8/xl                                                   
   enddo
!
!  initialize array for flow-blocking drag
!
   taufb(its:ite,kts:kte+1) = 0.0_r8
   komax(its:ite) = 0
!
   do k = kts,kte
     do i = its,ite
       vtj(i,k)  = t1(i,k)  * (1._r8+fv*q1(i,k))
       vtk(i,k)  = vtj(i,k) / prslk(i,k)
       ro(i,k)   = 1._r8/rd * prsl(i,k) / vtj(i,k) ! density kg/m**3
     enddo
   enddo

!
!  determine reference level: maximum of 2*var and pbl heights
!
   do i = its,ite
     zlowtop(i) = 2._r8 * var(i)
   enddo
!
   do i = its,ite
     kloop1(i) = .true.
   enddo
!
   do k = kts+1,kte
     do i = its,ite
         if(kloop1(i).and.zl(i,k)-zl(i,1).ge.zlowtop(i)) then
         klowtop(i) = k+1
         kloop1(i)  = .false.
         endif
     enddo
   enddo
!
   do i = its,ite
     kbl(i)   = max(kpbl(i), klowtop(i))
     kbl(i)   = max(min(kbl(i),kpblmax),kpblmin)
   enddo
!
!  determine the level of maximum orographic height
!
   komax(:) = kbl(:)
!
   do i = its,ite
     delks(i)  = 1.0_r8 / (prsi(i,1) - prsi(i,kbl(i)))
     delks1(i) = 1.0_r8 / (prsl(i,1) - prsl(i,kbl(i)))
   enddo


!
!  compute low level averages within pbl
!
   do k = kts,kpblmax
     do i = its,ite
       if (k.lt.kbl(i)) then
         rcsks   = rcs     * del(i,k) * delks(i)
         rdelks  = del(i,k)  * delks(i)
         ubar(i) = ubar(i) + rcsks  * u1(i,k)      ! pbl u  mean
         vbar(i) = vbar(i) + rcsks  * v1(i,k)      ! pbl v  mean
         roll(i) = roll(i) + rdelks * ro(i,k)      ! ro mean
       endif
     enddo
   enddo
!

!=======Jinbo Xie=======
!For ls and bl only
IF  ((gsd_gwd_ls .EQ. 1).or.(gsd_gwd_bl .EQ. 1)) then
!     figure out low-level horizontal wind direction 
!====Jinbo Xie order into a counterclockwise index instead====
!no more 1-8 index
   do i = its,ite                                                       
     wdir   = atan2(vbar(i),ubar(i)) + pi!changed into y/x Jinbo Xie
     !idir   = MOD(nint(fdir*wdir),mdir) + 1!starts from pi already
      !nwd   = idir
#ifndef continuous
        wdir1=wdir-pi
        if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
        nwd  = MOD(nint(fdir*wdir1),mdir) + 1
        else!(-pi,0)
        nwd  = MOD(nint(fdir*(wdir1+2._r8*pi)),mdir) + 1
        endif
        !turn backwords because start is pi
        !need turning
        rad=4.0_r8*atan(1.0_r8)/180.0_r8
        theta=(real(nwd,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))

wdir1_xjb(i)=wdir1/rad

        oa1(i)= oa4(i,1)*cos(theta*rad)+oa4(i,2)*sin(theta*rad)
        !select OL
        ol(i)  = ol4(i,MOD(nwd-1,int(mdir/2))+1)
        !calculate dxygrid, not so slow
        call dxygrid(dxmeter(i),dymeter,theta,dxy(i))



!----- compute orographic width along (ol) and perpendicular (olp)
!----- the direction of wind
!
!====Jinbo Xie====
!put wdir inside the (0,2*pi) section
!changing pi/2 either way is perpendicular
                !wdir1=wdir-pi
                if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
                nwd1  = MOD(nint(fdir*(wdir1+pi/2._r8)),mdir) + 1
                olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
                else!(-pi,0)
                nwd1  = MOD(nint(fdir*(wdir1-pi/2._r8+2._r8*pi)),mdir) + 1
                olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
                endif
                theta=(real(nwd1,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))
                call dxygrid(dxmeter(i),dymeter,theta,dxyp(i))
!====Jinbo Xie====
#else
            !determine oa by wind direction angle
            oa1(i)  = oa4(i,1)*cos(wdir-pi)+oa4(i,2)*sin(wdir-pi)!oalon,oalat
            !number of nonfillvalue variables
            !assume fillvalue 1.d36 in data
            !rather loose threshold to 10 below 1.d36
            nind_non=count(abs(ol4(1,i,:)).lt.1.d36-10)
            rad=4.0_r8*atan(1.0_r8)/180.0_r8
            wdir1=wdir-pi
                if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
                theta  = wdir1/rad
                else !(-pi,0) to (pi,2*pi)
                theta  = (wdir1+2._r8*pi)/rad
                endif
            call OLgrid(ol4(1,i,:nind_non),ol4(2,i,:nind_non),&
                        ol4(3,i,:nind_non),ol4(4,i,:nind_non),&
                        dxmeter(i),dymeter,&
                        nind_non,theta,1116.2_r8-0.878_r8*var(i),ol(i))
            call dxygrid(dxmeter(i),dymeter,theta,dxy(i))
                if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
                theta  = (wdir1+pi/2._r8)/rad
                else !(-pi,0)
                theta  = (wdir1-pi/2._r8+2._r8*pi)/rad
                endif
            call OLgrid(ol4(1,i,:nind_non),ol4(2,i,:nind_non),&
                        ol4(3,i,:nind_non),ol4(4,i,:nind_non),&
                        dxmeter(i),dymeter,&
                        nind_non,theta,1116.2_r8-0.878_r8*var(i),olp(i))
            call dxygrid(dxmeter(i),dymeter,theta,dxyp(i))
!Jinbo Xie
#endif
!
!----- compute orographic direction (horizontal orographic aspect ratio)
!
     od(i) = olp(i)/max(ol(i),olmin)
     od(i) = min(od(i),odmax)
     od(i) = max(od(i),odmin)
!
!----- compute length of grid in the along(dxy) and cross(dxyp) wind directions
!
!==========================================
!Jinbo Xie
   enddo
!===Jinbo Xie===
ENDIF
!===Jinbo Xie===

!Jinbo Xie Since variable grid,change dxy4 to larger
!==========================================
!
! END INITIALIZATION; BEGIN GWD CALCULATIONS:
!
IF ( ((gsd_gwd_ls .EQ. 1).or.(gsd_gwd_bl .EQ. 1)).and.   &
               (ls_taper .GT. 1.E-02) ) THEN   !====

!                                                                       
!---  saving richardson number in usqj for migwdi                       
!
   do k = kts,kte-1                                                     
     do i = its,ite                                                     
       ti        = 2.0_r8 / (t1(i,k)+t1(i,k+1)) 
       rdz       = 1._r8/(zl(i,k+1) - zl(i,k))
       tem1      = u1(i,k) - u1(i,k+1)
       tem2      = v1(i,k) - v1(i,k+1)   
       dw2       = rcl*(tem1*tem1 + tem2*tem2)
       shr2      = max(dw2,dw2min) * rdz * rdz
       bvf2      = g*(g/cp+rdz*(vtj(i,k+1)-vtj(i,k))) * ti                
       usqj(i,k) = max(bvf2/shr2,rimin)                            
       !bnv2(i,k) = 2.0_r8*g*rdz*(vtk(i,k+1)-vtk(i,k))/(vtk(i,k+1)+vtk(i,k))
       !bnv2(i,k) = max( bnv2(i,k), bnv2min )
       bnv2(i,k) = max(bnv_in(i,k)**2,bnv2min )
     enddo                                                          
   enddo                                                             

!
!----compute the "low level" or 1/3 wind magnitude (m/s)                
!                                                                       
   do i = its,ite                                                       
     ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0_r8)
     rulow(i) = 1._r8/ulow(i)
   enddo                                                             
!
   do k = kts,kte-1                                                    
     do i = its,ite                                                   
       velco(i,k)  = (0.5_r8*rcs) * ((u1(i,k)+u1(i,k+1)) * ubar(i)                &
                                   + (v1(i,k)+v1(i,k+1)) * vbar(i))                 
       velco(i,k)  = velco(i,k) * rulow(i)                               
       if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0._r8)) then
         velco(i,k) = veleps                                      
       endif
     enddo                                                          
   enddo                                                             
!                                                                       
!  no drag when critical level in the base layer                        
!                                                                       
   do i = its,ite                                                       
     ldrag(i) = velco(i,1).le.0._r8
   enddo                                                             
!
!  no drag when velco.lt.0                                               
!                             
                                          
   do k = kpblmin,kpblmax
     do i = its,ite                                                    
       if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. velco(i,k).le.0._r8
     enddo                                                          
   enddo                                                             
!                                                                       
!  no drag when bnv2.lt.0                                               
!                                                                       
   do k = kts,kpblmax
     do i = its,ite                                                    
       if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. bnv2(i,k).lt.0._r8
     enddo                                                          
   enddo                                                             

!                                                                       
!-----the low level weighted average ri is stored in usqj(1,1; im)      
!-----the low level weighted average n**2 is stored in bnv2(1,1; im)    
!---- this is called bnvl2 in phys_gwd_alpert_sub not bnv2                           
!---- rdelks (del(k)/delks) vert ave factor so we can * instead of /    
!                                                                       
   do i = its,ite                                                       
     wtkbj     = (prsl(i,1)-prsl(i,2)) * delks1(i)
     bnv2(i,1) = wtkbj * bnv2(i,1)                                
     usqj(i,1) = wtkbj * usqj(i,1)                                
   enddo                                                             
!
   do k = kpblmin,kpblmax                                                
     do i = its,ite                                                    
       if (k .lt. kbl(i)) then
         rdelks    = (prsl(i,k)-prsl(i,k+1)) * delks1(i)
         bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
         usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
       endif
     enddo                                                          
   enddo                                                             
!                                                                       
   do i = its,ite                                                       
     ldrag(i) = ldrag(i) .or. bnv2(i,1).le.0.0_r8                        
     ldrag(i) = ldrag(i) .or. ulow(i).eq.1.0_r8                       
     ldrag(i) = ldrag(i) .or. var(i) .le. 0.0_r8
   enddo                                                             

!                                                                       
!  set all ri low level values to the low level value          
!                                                                       
   do k = kpblmin,kpblmax
     do i = its,ite                                                    
       if (k .lt. kbl(i)) usqj(i,k) = usqj(i,1)
     enddo                                                          
   enddo                                                             
!
   do i = its,ite 
     if (.not.ldrag(i))   then   
       bnv(i) = sqrt( bnv2(i,1) )                                  
       fr(i) = bnv(i)  * rulow(i) * 2._r8 * var(i) * od(i)
       fr(i) = min(fr(i),frmax)
       xn(i)  = ubar(i) * rulow(i)
       yn(i)  = vbar(i) * rulow(i)
     endif
   enddo
!
!  compute the base level stress and store it in taub
!  calculate enhancement factor, number of mountains & aspect        
!  ratio const. use simplified relationship between standard            
!  deviation & critical hgt                                          
!

   do i = its,ite                                                       
     if (.not. ldrag(i))   then   
       efact    = (oa1(i) + 2._r8) ** (ce*fr(i)/frc)                         
       efact    = min( max(efact,efmin), efmax )                            
!!!!!!! cleff (effective grid length) is highly tunable parameter
!!!!!!! the bigger (smaller) value produce weaker (stronger) wave drag
       cleff    = sqrt(dxy(i)**2._r8 + dxyp(i)**2._r8)
!==============Jinbo Xie=============================================
       !cleff    = 3._r8 * max(dxmeter(i),cleff)!turned dxmeter to array
        cleff    = 3._r8 * max(dxmax_ls,cleff)
        !cleff    = 2._r8 * max(dxmax_ls,cleff)
!==============Jinbo Xie=============================================
       coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8) 
       xlinv(i) = coefm(i) / cleff                                             
       tem      = fr(i) * fr(i) * oc1(i)
       gfobnv   = gmax * tem / ((tem + cg)*bnv(i))   

       if ( gsd_gwd_ls .NE. 0 ) then
          taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)                       &
                   * ulow(i) * gfobnv * efact          
       else     ! We've gotten what we need for the blocking scheme
          taub(i) = 0.0_r8
       end if
!Jinbo Xie for base flux
taub_xjb(i)=taub(i)
!Jinbo Xie for base flux
     else                                                          
       taub(i) = 0.0_r8                         
       xn(i)   = 0.0_r8                         
       yn(i)   = 0.0_r8                         
     endif                                                         
   enddo

ENDIF   ! (gsd_gwd_ls .EQ. 1).or.(gsd_gwd_bl .EQ. 1)
!=========================================================
! add small-scale wavedrag for stable boundary layer
!=========================================================
  XNBV=0._r8
  tauwavex0=0._r8
  tauwavey0=0._r8
  density=1.2_r8
  utendwave=0._r8
  vtendwave=0._r8
  zq=0._r8
!
  IF ( (gsd_gwd_ss .EQ. 1).and.(ss_taper.GT.1.E-02) ) THEN
!
! declaring potential temperature
!
    do k = kts,kte
      do i = its,ite
        thx(i,k) = t1(i,k)/prslk(i,k)
      enddo
    enddo
!
    do k = kts,kte
      do i = its,ite
        tvcon = (1._r8+fv*q1(i,k))
        thvx(i,k) = thx(i,k)*tvcon
      enddo
    enddo
!
! Defining layer height
!
    do k = kts,kte
      do i = its,ite
        zq(i,k+1) = dz2(i,k)+zq(i,k)
      enddo
    enddo
!
    do k = kts,kte
      do i = its,ite
        za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
      enddo
    enddo

    do i=its,ite
       hpbl2 = hpbl(i)+10._r8
       kpbl2 = kpbl(i)
       kvar = 1
       do k=kts+1,MAX(kpbl(i),kts+1)
          IF (za(i,k)>300._r8) then
             kpbl2 = k
             IF (k == kpbl(i)) then
                hpbl2 = hpbl(i)+10._r8
             ELSE
                hpbl2 = za(i,k)+10._r8
             ENDIF
             exit
          ENDIF
       enddo

       if((xland1(i)-1.5_r8).le.0._r8 .and. 2._r8*var(i).le.hpbl(i))then
          if(br1(i).gt.0._r8 .and. thvx(i,kpbl2)-thvx(i,kts) > 0._r8)then
            cleff    = sqrt(dxy(i)**2_r8 + dxyp(i)**2_r8)
            cleff    = 2.0_r8 * max(dxmax_ss,cleff)
            coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8)
            xlinv(i) = coefm(i) / cleff
            govrth(i)=g/(0.5_r8*(thvx(i,kpbl2)+thvx(i,kts)))
            XNBV=sqrt(govrth(i)*(thvx(i,kpbl2)-thvx(i,kts))/hpbl2)
!
            if(abs(XNBV/u1(i,kpbl2)).gt.xlinv(i))then
              tauwavex0=0.5_r8*XNBV*xlinv(i)*(2._r8*MIN(var(i),varmax))**2_r8*ro(i,kvar)*u1(i,kvar)
              tauwavex0=tauwavex0*ss_taper   ! "Scale-awareness"
            else
              tauwavex0=0._r8
            endif
!
            if(abs(XNBV/v1(i,kpbl2)).gt.xlinv(i))then
              tauwavey0=0.5_r8*XNBV*xlinv(i)*(2._r8*MIN(var(i),varmax))**2._r8*ro(i,kvar)*v1(i,kvar)
              tauwavey0=tauwavey0*ss_taper   ! "Scale-awareness"
            else
              tauwavey0=0._r8
            endif
!
            do k=kts,kpbl(i) !MIN(kpbl2+1,kte-1)
              utendwave(i,k)=-1._r8*tauwavex0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
              vtendwave(i,k)=-1._r8*tauwavey0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
            enddo
          endif
       endif
    enddo ! end i loop

    do k = kts,kte
       do i = its,ite
         dudt(i,k)  = dudt(i,k) + utendwave(i,k)
         dvdt(i,k)  = dvdt(i,k) + vtendwave(i,k)
         dtaux2d_ss(i,k) = utendwave(i,k)
         dtauy2d_ss(i,k) = vtendwave(i,k)
         dusfc_ss(i) = dusfc_ss(i) + utendwave(i,k) * del(i,k)
         dvsfc_ss(i) = dvsfc_ss(i) + vtendwave(i,k) * del(i,k)
       enddo
    enddo

ENDIF  ! end if gsd_gwd_ss == 1
!================================================================
!add Beljaars et al. (2004, QJRMS, equ. 16) form drag:
!================================================================
IF ( (gsd_gwd_fd .EQ. 1).and.(ss_taper.GT.1.E-02) ) THEN

   utendform=0._r8
   vtendform=0._r8
   zq=0._r8

   IF ( (gsd_gwd_ss .NE. 1).and.(ss_taper.GT.1.E-02) ) THEN
      ! Defining layer height. This is already done above is small-scale GWD is used
      do k = kts,kte
        do i = its,ite
          zq(i,k+1) = dz2(i,k)+zq(i,k)
        enddo
      enddo

      do k = kts,kte
        do i = its,ite
          za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
        enddo
      enddo
   ENDIF

   DO i=its,ite
      IF ((xland1(i)-1.5) .le. 0.) then
          a1=0.00026615161_r8*var(i)**2_r8
          a2=a1*0.005363_r8
         DO k=kts,kte
            wsp=SQRT(u1(i,k)**2_r8 + v1(i,k)**2_r8)
            ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759 
            utendform(i,k)=-0.0759_r8*wsp*u1(i,k)* &
                           EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
            vtendform(i,k)=-0.0759_r8*wsp*v1(i,k)* &
                           EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
            !IF(za(i,k) > 4000.) exit
         ENDDO
      ENDIF
   ENDDO

   do k = kts,kte
      do i = its,ite
         dudt(i,k)  = dudt(i,k) + utendform(i,k)
         dvdt(i,k)  = dvdt(i,k) + vtendform(i,k)
         dtaux2d_fd(i,k) = utendform(i,k)
         dtauy2d_fd(i,k) = vtendform(i,k)
         dusfc_fd(i) = dusfc_fd(i) + utendform(i,k) * del(i,k)
         dvsfc_fd(i) = dvsfc_fd(i) + vtendform(i,k) * del(i,k)
      enddo
   enddo

ENDIF  ! end if gsd_gwd_fd == 1
!=======================================================
! More for the large-scale gwd component
IF ( (gsd_gwd_ls .EQ. 1).and.(ls_taper.GT.1.E-02) ) THEN
!                                                                       
!   now compute vertical structure of the stress.
!
   do k = kts,kpblmax
      do i = its,ite
         if (k .le. kbl(i)) taup(i,k) = taub(i)
      enddo
   enddo
!

!================================
!Jinbo Xie
!determination of the interface height
do i=its,ite
iint=.false.
        do k=kpblmin,kte-1
        if (k.gt.kbl(i).and.usqj(1,k)-usqj(1,k-1).lt.0.and.(.not.iint)) then
        iint=.true.
        zl_hint(i)=zl(i,k+1)
        endif
        enddo
enddo
!print*,"zl_hint",zl_hint
!!stop
!Jinbo Xie
!================================



   do k = kpblmin, kte-1                   ! vertical level k loop!
      kp1 = k + 1
      do i = its,ite
!
!   unstablelayer if ri < ric
!   unstable layer if upper air vel comp along surf vel <=0 (crit lay)
!   at (u-c)=0. crit layer exists and bit vector should be set (.le.)
!
         if (k .ge. kbl(i)) then
           !icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric)                  &
           !                      .or. (velco(i,k) .le. 0.0_r8)
!============================
!Jinbo Xie
!we modify the criteria for unstable layer
!that the lv is critical under 0.25
!while we keep wave breaking ric for
!other larger lv
           icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric_rig)&
                                 .or. (velco(i,k) .le. 0.0_r8)
!Jinbo Xie
!============================
           brvf(i)  = max(bnv2(i,k),bnv2min) ! brunt-vaisala frequency squared
           brvf(i)  = sqrt(brvf(i))          ! brunt-vaisala frequency
         endif
      enddo
!
      do i = its,ite
        if (k .ge. kbl(i) .and. (.not. ldrag(i)))   then   
          if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0_r8 ) then
            temv = 1.0_r8 / velco(i,k)
            tem1 = coefm(i)/dxy(i)*(ro(i,kp1)+ro(i,k))*brvf(i)*velco(i,k)*0.5_r8
            hd   = sqrt(taup(i,k) / tem1)
            fro  = brvf(i) * hd * temv

!
!  rim is the minimum-richardson number by shutts (1985)
!
            tem2   = sqrt(usqj(i,k))
            tem    = 1._r8 + tem2 * fro
            rim    = usqj(i,k) * (1._r8-fro) / (tem * tem)

!
!  check stability to employ the 'saturation hypothesis'
!  of lindzen (1981) except at tropospheric downstream regions
!
            if (rim .le. ric) then  ! saturation hypothesis!
              if ((oa1(i) .le. 0._r8).or.(kp1 .ge. kpblmin )) then
                temc = 2.0_r8 + 1.0_r8 / tem2
                hd   = velco(i,k) * (2.0_r8*sqrt(temc)-temc) / brvf(i)
                taup(i,kp1) = tem1 * hd * hd

!==============================================
!taup is restricted to monotoncally decrease
!to avoid unexpected high taup with taup cal
taup(i,kp1)=min(tem1*hd*hd,taup(i,k))
!add vertical decrease at low level below hint (Kim and Doyle 2005)
!where Ri first decreases
!#if 0
if (k.gt.klowtop(i).and.zl(i,k).le.zl_hint(i)) then
l1=(9.81_r8*bnv2(i,kp1)/velco(i,kp1)**2)!-(shr2_xjb(i,kp1)/velco(i,kp1))
l2=(9.81_r8*bnv2(i,k)/velco(i,k)**2)!-(shr2_xjb(i,k)/velco(i,k))
!print*,"l1,l2,l1/l2",l1,l2,l1/l2
taup(i,kp1)=min(taup(i,k),taup(i,k)*(l1/l2),tem1*hd*hd)
!taup(i,kp1)=max(0.2*taup(i,k),min(taup(i,k),taup(i,k)*(l1/l2),tem1*hd*hd))
!taup(i,k)*(l1/l2)
!print*,"taup(i,kp1)",taup(i,kp1)
!print*,"k",k
endif
!#endif
!==============================================
              endif
            else                    ! no wavebreaking!
              taup(i,kp1) = taup(i,k)
            endif
          endif
        endif
      enddo      
   enddo
!


   if(lcap.lt.kte) then                                               
      do klcap = lcapp1,kte                                          
         do i = its,ite                                                 
           taup(i,klcap) = prsi(i,klcap) / prsi(i,lcap) * taup(i,lcap)      
         enddo                                                       
      enddo                                                          
   endif      

ENDIF !END LARGE-SCALE TAU CALCULATION

!===============================================================
!COMPUTE BLOCKING COMPONENT 
!===============================================================
IF ( (gsd_gwd_bl .EQ. 1) .and. (ls_taper .GT. 1.E-02) ) THEN
                                                       
   do i = its,ite
      if(.not.ldrag(i)) then
!
!------- determine the height of flow-blocking layer
!
        kblk = 0
        pe = 0.0_r8

        do k = kte, kpblmin, -1
          if(kblk.eq.0 .and. k.le.komax(i)) then
!Jinbo Xie
!flow block appears within the reference level
!compare potential energy and kinetic energy
!divided by g*ro is to turn del(pa) into height
            pe = pe + bnv2(i,k)*(zl(i,komax(i))-zl(i,k))*del(i,k)/g/ro(i,k)
            ke = 0.5_r8*((rcs*u1(i,k))**2._r8+(rcs*v1(i,k))**2._r8)
!
!---------- apply flow-blocking drag when pe >= ke 
!
            if(pe.ge.ke) then
              kblk = k
              kblk = min(kblk,kbl(i))
              zblk = zl(i,kblk)-zl(i,kts)
!zblk_col(i)=zblk
            endif
          endif
        enddo


        if(kblk.ne.0) then
!
!--------- compute flow-blocking stress
!

!Jinbo Xie the max(dxmax_ls,dxy(i))**2
!Jinbo Xie here is a crude estimate since the cam is uneven 0.9*1.25deg
!dxmax_ls is different than the usual one
!because the taper is very different
!Jinbo Xie dxy is a length scale mostly in the direction of the flow to the ridge
!so it is good and not needed for an uneven grid area
!ref Lott and Miller (1997) original scheme
          cd = max(2.0_r8-1.0_r8/od(i),0.0_r8)
          taufb(i,kts) = 0.5_r8 * roll(i) * coefm(i) / max(dxmax_ls,dxy(i))**2 * cd * dxyp(i)   &
                         * olp(i) * zblk * ulow(i)**2
!Jinbo Xie for base flux
taufb_xjb(i)=taufb(i,kts)
!Jinbo Xie for base flux

        !changed grid box area into dy*dy
          tautem = taufb(i,kts)/float(kblk-kts)
          do k = kts+1, kblk
            taufb(i,k) = taufb(i,k-1) - tautem
          enddo
!
!----------sum orographic GW stress and flow-blocking stress
!
           !taup(i,:) = taup(i,:) + taufb(i,:)   ! Keep taup and taufb separate for now
        endif
      endif
   enddo 

ENDIF   ! end blocking drag
!===========================================================
IF ( (gsd_gwd_ls .EQ. 1 .OR. gsd_gwd_bl .EQ. 1) .and. (ls_taper .GT. 1.E-02) ) THEN

!                                                                       
!  calculate - (g)*d(tau)/d(pressure) and deceleration terms dtaux, dtauy
!
   do k = kts,kte                                                       
     do i = its,ite                                                       
       taud_ls(i,k) = 1._r8 * (taup(i,k+1) - taup(i,k)) * csg / del(i,k)
       taud_bl(i,k) = 1._r8 * (taufb(i,k+1) - taufb(i,k)) * csg / del(i,k)
     enddo                                                             
   enddo                                                             
!                                                                       
!  limit de-acceleration (momentum deposition ) at top to 1/2 value 
!  the idea is some stuff must go out the 'top'                     
!                                                                       

   do klcap = lcap,kte                                               
     do i = its,ite                                                    
       taud_ls(i,klcap) = taud_ls(i,klcap) * factop
       taud_bl(i,klcap) = taud_bl(i,klcap) * factop
     enddo                                                          
   enddo                                                             

!                                                                       
!  if the gravity wave drag would force a critical line             
!  in the lower ksmm1 layers during the next deltim timestep,     
!  then only apply drag until that critical line is reached.        
!                                                                       
   do k = kts,kpblmax-1
      do i = its,ite                                                    
         if (k .le. kbl(i)) then
           if((taud_ls(i,k)+taud_bl(i,k)).ne.0._r8)                      &
              dtfac(i) = min(dtfac(i),abs(velco(i,k)                     &
                   /(deltim*rcs*(taud_ls(i,k)+taud_bl(i,k)))))
         endif
      enddo
   enddo
!



   do k = kts,kte                                                       
      do i = its,ite 
         taud_ls(i,k)  = taud_ls(i,k) * dtfac(i) * ls_taper
         taud_bl(i,k)  = taud_bl(i,k) * dtfac(i) * ls_taper
         dtaux2d_ls(i,k) = taud_ls(i,k) * xn(i)
         dtauy2d_ls(i,k) = taud_ls(i,k) * yn(i)
         dtaux2d_bl(i,k) = taud_bl(i,k) * xn(i)
         dtauy2d_bl(i,k) = taud_bl(i,k) * yn(i)
         dudt(i,k)  = dtaux2d_ls(i,k) + dtaux2d_bl(i,k) + dudt(i,k)
         dvdt(i,k)  = dtauy2d_ls(i,k) + dtauy2d_bl(i,k) + dvdt(i,k)
         dusfc_ls(i)  = dusfc_ls(i) + dtaux2d_ls(i,k) * del(i,k)
         dvsfc_ls(i)  = dvsfc_ls(i) + dtauy2d_ls(i,k) * del(i,k)
         dusfc_bl(i)  = dusfc_bl(i) + dtaux2d_bl(i,k) * del(i,k)
         dvsfc_bl(i)  = dvsfc_bl(i) + dtauy2d_bl(i,k) * del(i,k)
      enddo                                                          
   enddo

ENDIF                                                             


!  Finalize dusfc and dvsfc diagnoses
do i = its,ite
   dusfc_ls(i) = (-1._r8/g*rcs) * dusfc_ls(i)
   dvsfc_ls(i) = (-1._r8/g*rcs) * dvsfc_ls(i)
   dusfc_bl(i) = (-1._r8/g*rcs) * dusfc_bl(i)
   dvsfc_bl(i) = (-1._r8/g*rcs) * dvsfc_bl(i)
   dusfc_ss(i) = (-1._r8/g*rcs) * dusfc_ss(i)
   dvsfc_ss(i) = (-1._r8/g*rcs) * dvsfc_ss(i)
   dusfc_fd(i) = (-1._r8/g*rcs) * dusfc_fd(i)
   dvsfc_fd(i) = (-1._r8/g*rcs) * dvsfc_fd(i)
enddo

!#Jinbo get base flux
do i = its,ite
dusfc_ss(i)=wdir1_xjb(i)
dvsfc_ss(i)=ol(i)
dtaux2d_ss(i,1)=oa1(i)
dtaux2d_ss(i,2)=ol(i)
dtaux2d_ss(i,3)=olp(i)
dtaux2d_ss(i,4)=dxy(i)
dtaux2d_ss(i,5)=dxyp(i)
dtaux2d_ss(i,6)=var(i)!hpbl(i)!zblk_xjb(i)!vtj(i,2)-vtj(i,1)
dtaux2d_ss(i,7)=klowtop(i)!pe_xjb(i)!bnv2(i,2)!vtk(i,k+1)-vtk(i,k)
dtaux2d_ss(i,8)=kpbl(i)
dtaux2d_ss(i,9)=zblk_col(i)
dtaux2d_ss(i,10)=oa4(i,1)
dtaux2d_ss(i,11)=oa4(i,2)
dtaux2d_ss(i,12)=oa4(i,3)

!print*,"wdir1_xjb(i)",wdir1_xjb(i)
!u1*q1
!v1*q1
enddo
!stop
#if 0
do i = its,ite
dusfc_ss(i)=taub_xjb(i)
dvsfc_ss(i)=taufb_xjb(i)
enddo
#endif

#if 0
do i = its,ite
dtaux2d_ls(i,1)=oa1(i)
dtaux2d_ls(i,2)=ol(i)
dtaux2d_ls(i,3)=olp(i)
dtaux2d_ls(i,4)=dxy(i)
dtaux2d_ls(i,5)=dxyp(i)
dtaux2d_ls(i,6)=var(i)!hpbl(i)!zblk_xjb(i)!vtj(i,2)-vtj(i,1)
dtaux2d_ls(i,7)=klowtop(i)!pe_xjb(i)!bnv2(i,2)!vtk(i,k+1)-vtk(i,k)
dtaux2d_ls(i,8)=kpbl(i)
dtaux2d_ls(i,9)=zblk_col(i)
dtaux2d_ls(i,10)=oa4(i,1)
dtaux2d_ls(i,11)=oa4(i,2)
dtaux2d_ls(i,12)=oa4(i,3)
enddo
#endif
   return                                                            
   end subroutine gwdo2d
!-------------------------------------------------------------------
!-czy20181120
end module gw_drag_xjb

