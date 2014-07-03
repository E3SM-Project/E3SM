!---------------------------------------------------------------------------------
! this file came from wa17 and modified by Fabrizio 07-02-2004
! standard gw_drag with modification (6) of latitude profile of gw spectrum
!---------------------------------------------------------------------------------

module gw_drag

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an 
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,   only : r8 => shr_kind_r8
  use ppgrid,         only : pcols, pver, pverp
  use constituents,   only : pcnst
  use physics_types,  only : physics_state, physics_ptend, physics_ptend_init
  use spmd_utils,     only : masterproc
  use cam_history,    only : outfld
  use physconst,      only : gravit, pi
  use cam_logfile,    only : iulog
  use abortutils,     only : endrun
  use molec_diff,     only : nbot_molec
 
  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public :: gw_inti                  ! Initialization
  public :: gw_intr                  ! interface to actual parameterization
  public :: gw_drag_readnl           ! Read namelist
  public :: gw_drag_register
  public :: idx_zmdt

!
! PRIVATE: Rest of the data and interfaces are private to this module
!
  integer,  parameter :: pgwv = 32           ! number of waves allowed (one-signed only)
  real(r8), parameter :: dc   = 2.5_r8       ! bin width for spectrum
  real(r8), parameter :: dc2  = 5.0_r8       ! bin width for output spectrum (for outputting tauxs)
  real(r8), parameter :: dca  = 0.1_r8       ! integration interval to get bin average
  real(r8), parameter :: spec_fac=10._r8/dc  ! number of bins compared to dc=10	!bab
  real(r8), parameter :: c0   = 30._r8       ! width of gaussian in phase speed
  real(r8) :: fav(-pgwv:pgwv)                ! average value of gaussian over bin
  integer :: kbotbg, kbotoro                 ! interface of gwd source
  integer :: ktopbg, ktoporo                 ! top interface of gwd region


!+Mods for Beres04
  integer :: bkbotbg                         ! interface of gwd beres source
  integer :: k700                            ! 700 mb index (used in the gw_bgnd_beres)
  real(r8), parameter :: effgw_beres=0.1_r8  ! beres source tendency efficiency
!-Mods for Beres04

!+Mods for CM
  real(r8) :: cref(-pgwv:pgwv)               ! Reference phase speed spectrum
  real(r8) :: cref2(-pgwv:pgwv)               ! Reference phase speed spectrum 2 (just for outputting tauxs)
  real(r8) :: frontgfc                       ! frontogenesis function critical threshold
!-Mods for CM


  real(r8) :: alpha(0:pver)                  ! newtonian cooling coefficients
  real(r8) :: cpair                          ! specific heat of dry air (constant p)
  real(r8) :: annphase                       ! phase of annual cycle for taubgnd (day of max)
  real(r8) :: annavg_n, annamp_n             ! NH annual average and amplitude for taubgnd
  real(r8) :: annavg_s, annamp_s             ! SH annual average and amplitude for taubgnd
  real(r8) :: dback                          ! background diffusivity
  real(r8) :: effkwv                         ! effective wavenumber (fcrit2*kwv)
  real(r8) :: effgw_spec                     ! efficiency factor for spectrum !bab
  real(r8) :: effgw_oro                      ! efficiency factor for orographic wave
  real(r8), parameter :: fracldv    = 0._r8  ! fraction of stress deposited in low level region
  real(r8) :: g                   ! acceleration of gravity
  real(r8) :: kwv                 ! effective horizontal wave number
  real(r8) :: lzmax               ! maximum vertical wavelength at source

  real(r8) :: mxasym              ! max asymmetry between tau(c) and tau(-c)
  real(r8) :: mxrange             ! max range of tau for all c

  real(r8) :: n2min               ! min value of bouyancy frequency
  real(r8) :: fcrit2              ! critical froude number
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

  integer, parameter :: nalph=66  ! no. of levels of pre-calculated newtonian cooling
  real(r8) :: alpha0(nalph)       ! newtonian cooling (1/day)
  real(r8) :: palph(nalph)        ! pressure levs of pre-calculated newtonian cooling (hpa)

  real(r8), parameter :: d2r  = pi/180._r8

  data alpha0 /  1.896007_r8   , 1.196965_r8   , 0.7251356_r8  ,  0.6397463_r8   , &
                 0.5777858_r8  , 0.5712274_r8  , 0.6836302_r8  ,  0.6678557_r8   , &
                 0.5683219_r8  , 0.4754283_r8  , 0.3960519_r8  ,  0.332022_r8    , &
                 0.2497581_r8  , 0.168667_r8   , 0.1323903_r8  ,  0.1257139_r8   , &
                 0.1069889_r8  , 0.09873954_r8 , 0.09215571_r8 ,  0.09398635_r8  , &
                 0.1061087_r8  , 0.1294598_r8  , 0.1544743_r8  ,  0.1648226_r8   , &
                 0.1687332_r8  , 0.1691513_r8  , 0.1664987_r8  ,  0.159048_r8    , &
                 0.149292_r8   , 0.1351563_r8  , 0.1174998_r8  ,  0.09913579_r8  , &
                 0.08300615_r8 , 0.0707_r8     , 0.0615588_r8  ,  0.0542623_r8   , &
                 0.0478562_r8  , 0.04132157_r8 , 0.03454087_r8 ,  0.02296682_r8  , &
                 0.006723819_r8, 0.02164464_r8 , 0.05756261_r8 ,  0.003844868_r8 , &
                 0.02929285_r8 , 0.006627098_r8, 0.04558291_r8 ,  0.02042176_r8  , &
                 0.00000000_r8 , 0.005880283_r8, 0.00689498_r8 ,  0.01343466_r8  , &
                 0.00000000_r8 , 0.03415992_r8 , 0.02855049_r8 ,  0.01688839_r8  , &
                 0.0272628_r8  , 0.02772121_r8 , 0.02135626_r8 ,  0.04863235_r8  , &
                 0.04568304_r8 , 0.00000000_r8 , 0.009604108_r8,  0.00000000_r8  , &
                 0.00000000_r8 , 0.00000000_r8 /

  data palph /   5.11075e-06_r8 , 9.8269e-06_r8   , 1.620185e-05_r8 , 2.671225e-05_r8 , & 
                 4.4041e-05_r8  , 7.261275e-05_r8 , 0.000119719_r8  , 0.00019738_r8   , &
                 0.0003254225_r8, 0.0005365325_r8 , 0.0008846025_r8 , 0.001458458_r8  , &
                 0.002404575_r8 , 0.00397825_r8   , 0.006556825_r8  , 0.01081382_r8   , &
                 0.017898_r8    , 0.02955775_r8   , 0.04873075_r8   , 0.07991075_r8   , &
                 0.1282732_r8   , 0.19812_r8      , 0.292025_r8     , 0.4101675_r8    , & 
                 0.55347_r8     , 0.73048_r8      , 0.9559475_r8    , 1.244795_r8     , &  
                 1.61285_r8     , 2.079325_r8     , 2.667425_r8     , 3.404875_r8     , &
                 4.324575_r8    , 5.4654_r8       , 6.87285_r8      , 8.599725_r8     , & 
                 10.70705_r8    , 13.26475_r8     , 16.35175_r8     , 20.05675_r8     , &
                 24.479_r8      , 29.728_r8       , 35.92325_r8     , 43.19375_r8     , &
                 51.6775_r8     , 61.5205_r8      , 72.8745_r8      , 85.65715_r8     , &
                 100.5147_r8    , 118.2503_r8     , 139.1154_r8     , 163.6621_r8     , &
                 192.5399_r8    , 226.5132_r8     , 266.4812_r8     , 313.5013_r8     , &
                 368.818_r8     , 433.8952_r8     , 510.4553_r8     , 600.5242_r8     , &
                 696.7963_r8    , 787.7021_r8     , 867.1607_r8     , 929.6489_r8     , & 
                 970.5548_r8    , 992.5561_r8 /

!+ Mods for Beres04
! BERES parameterization arrays
  integer, parameter :: MAXH = 20                 ! max values of heating depth
  integer, parameter :: MAXUH = 40                ! max value of mean wind in heating
  real(r8) :: mfcc(MAXH,-MAXUH:MAXUH,-PGWV:PGWV)  ! stores source spectra
  character(len=256) :: gw_drag_file = 'Beres04_file'
!-Mods for Beres04

  integer :: idx_zmdt
  integer :: pblh_idx=0
  integer :: kvt_idx=0

contains

!===============================================================================
  
  subroutine gw_drag_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'gw_drag_readnl'

    namelist /gw_drag_nl/ gw_drag_file
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
    call mpibcast (gw_drag_file, len(gw_drag_file), mpichar, 0, mpicom)
#endif

  end subroutine gw_drag_readnl

  subroutine gw_drag_register
    use ppgrid,       only : pver
    
    use physics_buffer, only : pbuf_add_field, dtype_r8

    ! add fields to physics buffer
    call pbuf_add_field('ZMDT','physpkg',dtype_r8,(/pcols,pver/),idx_zmdt)
  end subroutine gw_drag_register

  subroutine gw_inti (cpairx, cpwv, gx, rx, pref_edge)
!-----------------------------------------------------------------------
! Time independent initialization for multiple gravity wave parameterization.
!-----------------------------------------------------------------------
    use cam_history,      only : addfld, add_default, phys_decomp
    use dycore,           only : get_resolution
    use interpolate_data, only : lininterp
  
    use physics_buffer, only : pbuf_get_index
    
!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: cpairx                ! specific heat of dry air (constant p)
    real(r8), intent(in) :: cpwv                  ! specific heat of water vapor (constant p)
    real(r8), intent(in) :: gx                    ! acceleration of gravity
    real(r8), intent(in) :: rx                    ! gas constant for dry air
    real(r8), intent(in) :: pref_edge(pverp)      ! reference interface pressures

!---------------------------Local storage-------------------------------
    integer :: k
    integer :: m, n, i, l
    character(len=80) fnum, fname, dumc1, dumc2, dumc1x, dumc1y   ! dummies
    real(r8) :: ccn
    real(r8) :: cmn,cmx                           ! bin edges (max and min values of c)

!-----------------------------------------------------------------------
! Copy model constants
    cpair  = cpairx
    g      = gx
    r      = rx

! Set MGWD constants
    kwv    = 6.28e-5_r8       ! 100 km wave length
    dback  = 0.05_r8          ! background diffusivity
    fcrit2 = 1.0_r8           ! critical froude number squared		
    tauscal= 0.001_r8         ! scale factor for background stress

    pblh_idx = pbuf_get_index('pblh')
    kvt_idx  = pbuf_get_index('kvt')


    select case (get_resolution())
    case ('4x5')
       effgw_oro  = 0.125_r8
       effgw_spec = 1.0_r8	! rrg
       taubgnd    = 1.0_r8	! jhr
       frontgfc   = 7.5e-16_r8
    case ('1.9x2.5')
       effgw_oro  = 0.125_r8
!      effgw_spec = 1.0_r8 / spec_fac 
       effgw_spec = 1.0_r8	! fs (04/27/2008)
       taubgnd    = 1.5_r8
       frontgfc   = 1.25E-15_r8  !jhr
    case default
       effgw_oro  = 0.125_r8
       effgw_spec = 1.0_r8
       taubgnd    = 1.5_r8
       frontgfc   = 1.25E-15_r8  !jhr
    end select

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'gw_inti: gw spectrum taubgnd, effgw_spec = ', taubgnd, effgw_spec
       write(iulog,*) ' '
    end if

    zldvcon = 10._r8           ! constant for determining zldv
    lzmax   = 7.e3_r8          ! maximum vertical wavelength at source (m)

!+Mods for CM
    cref(:) = 0._r8
    cref2(:) = 0._r8
!-Mods for CM

! the following block by BAB
!+++ bab

! If a spectrum is being computed, integrate gaussian over bins

    if (pgwv > 0) then
       if (masterproc) write(iulog,*) 'GW_DRAG: pgwv = ', pgwv				! rrg
       do l = -pgwv, pgwv
          cref(l) = dc * l                                              ! fs
          cref2(l) = dc2 * l                                            ! jhr	
! Lower bound of bin
          cmn = cref(l) - 0.5_r8*dc
          cmx = cref(l) + 0.5_r8*dc
! Loop over integration intervals in bin
          fav(l) = 0.5_r8 * dca * (exp(-(cmn/c0)**2) + exp(-(cmx/c0)**2))
          do n = 1,nint(dc/dca)-1
             fav(l) = fav(l) + dca * exp(-((cmn+n*dca)/c0)**2)
          enddo
          fav(l) = fav(l) / dc
          if (masterproc) then
             write (iulog,*) 'GW_DRAG: c(l), fav(l)= ', cref(l), fav(l)	! rrg
          end if
       enddo
    else
       fav(0) = 1._r8
    end if

!--- bab

! pre-calculated newtonian damping: 
!     * convert to 1/s
!     * ensure it is not smaller than 1e-6
!     * convert palph from hpa to pa

  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400._r8
     if (alpha0(k) .lt. 1.e-6_r8) alpha0(k) = 1.e-6_r8
     palph(k) = palph(k) *1.e2_r8
  enddo

! interpolate to current vertical grid and obtain alpha

  call lininterp (alpha0  ,palph, nalph , alpha  , pref_edge , pverp)
  if (masterproc) then
     write (iulog,*) 'gw_inti: newtonian damping (1/day):'
     write (iulog,fmt='(a4,a12,a10)') ' k  ','  pref_edge      ','  alpha   '
     do k=0,pver
        write (iulog,fmt='(i4,1e12.5,1f10.2)') k,pref_edge(k+1),alpha(k)*86400._r8
     enddo
  endif

! set radiative damping times
!!$    do k = 0, pver
!!$       alpha(k) = 1.e-6       ! about 10 days.
!!$    end do

! Min and max values to keep things reasonable
    mxasym = 0._r8		     ! max factor of 10 from |tau(c)| to |tau(-c)|
    mxrange= 0._r8                   ! factor of 100 from max to min |tau(c)|
    n2min   = 1.e-8_r8               ! min value of Brunt-Vaisalla freq squared
    orohmin = 10._r8                 ! min surface displacement for orographic wave drag
    orovmin = 2._r8                  ! min wind speed for orographic wave drag
    taumin  = 1.e-10_r8              ! min stress considered > 0
    tndmax  = 400._r8 / 86400._r8    ! rrg: max permitted tendency (500 m/s/day)
    umcfac  = 0.5_r8                 ! max permitted reduction in u-c
    ubmc2mn = 0.01_r8                ! min value of (u-c)^2

! Determine other derived constants
    oroko2 = 0.5_r8 * kwv
    effkwv = fcrit2 * kwv
    rog    = r/g

! Determine the bounds of the background and orographic stress regions
    ktopbg  = 0
    kbotoro = pver
    do k = 0, pver
       if (pref_edge(k+1) .lt. 50000._r8) kbotbg  = k    ! spectrum source at 500 mb
       if (pref_edge(k+1) .lt. 50000._r8) kbotbg  = k    ! spectrum source at 500 mb
!+ Mods for Beres04
       if (pref_edge(k+1) .lt. 10000._r8) bkbotbg  = k+1    ! beres spectrum source at 100 mb
       if (pref_edge(k+1) .lt. 70000._r8) k700    = k+1     ! 700 mb index
!- Mods for Beres04
    end do
    ktoporo = 0

    if (masterproc) then
       write (iulog,*) 'KTOPBG  =',ktopbg
       write (iulog,*) 'KBOTBG  =',kbotbg
       write (iulog,*) 'KTOPORO =',ktoporo
       write (iulog,*) 'KBOTORO =',kbotoro
!+Mods for Beres04
       write (iulog,*) 'BKBOTBG =',bkbotbg
       write (iulog,*) 'K700    =',k700
!+Mods for Beres04
    end if

! Declare history variables for orgraphic term
    call addfld ('TTGW','K/s     ',pver, 'A','T tendency - gravity wave drag',phys_decomp)
    call addfld ('UTGWORO ','m/s2    ',pver, 'A','U tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('VTGWORO ','m/s2    ',pver, 'A','V tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('TAUGWX  ','N/m2    ',1,    'A','Zonal gravity wave surface stress',        phys_decomp)
    call addfld ('TAUGWY  ','N/m2    ',1,    'A','Meridional gravity wave surface stress',   phys_decomp)

!   call add_default ('TTGW', 1, ' ')
    call add_default ('UTGWORO ', 1, ' ')
    call add_default ('VTGWORO ', 1, ' ')
    call add_default ('TAUGWX  ', 1, ' ')
    call add_default ('TAUGWY  ', 1, ' ')

! Declare history variables for spectrum
    if (pgwv > 0) then
       call addfld ('UTGWSPEC','m/s2    ',pver, 'A','U tendency - gravity wave spectrum',       phys_decomp)
       call addfld ('VTGWSPEC','m/s2    ',pver, 'A','V tendency - gravity wave spectrum',       phys_decomp)
       call add_default ('UTGWSPEC', 1, ' ')
       call add_default ('VTGWSPEC', 1, ' ')

       call addfld ('BUTGWSPEC ','m/s2      ',pver, 'A','Beres U tendency',phys_decomp)
       call addfld ('BVTGWSPEC ','m/s2      ',pver, 'A','Beres V tendency',phys_decomp)
       call add_default ('BUTGWSPEC', 1, ' ')
       call add_default ('BVTGWSPEC', 1, ' ')

! Gravity Wave Tendencies Separated by wave phase speed

       call addfld ('UTEND1  ','m/s2    ',pver, 'A','U tendency   c < -40',                phys_decomp)
       call addfld ('UTEND2  ','m/s2    ',pver, 'A','U tendency  -40 < c < -15',                phys_decomp)
       call addfld ('UTEND3  ','m/s2    ',pver, 'A','U tendency  -15 < c <  15',                phys_decomp)
       call addfld ('UTEND4  ','m/s2    ',pver, 'A','U tendency   15 < c <  40',                phys_decomp)
       call addfld ('UTEND5  ','m/s2    ',pver, 'A','U tendency   40 < c ',                phys_decomp)


       call addfld ('BUTEND1  ','m/s2    ',pver, 'A','U tendency   c < -40',                phys_decomp)
       call addfld ('BUTEND2  ','m/s2    ',pver, 'A','U tendency  -40 < c < -15',                phys_decomp)
       call addfld ('BUTEND3  ','m/s2    ',pver, 'A','U tendency  -15 < c <  15',                phys_decomp)
       call addfld ('BUTEND4  ','m/s2    ',pver, 'A','U tendency   15 < c <  40',                phys_decomp)
       call addfld ('BUTEND5  ','m/s2    ',pver, 'A','U tendency   40 < c',                phys_decomp)


       call addfld ('NETDT  ','K/s   ',pver, 'A','Net heating rate',                             phys_decomp)
       call addfld ('MAXQ0  ','K/day   ',1  ,  'A','Max column heating rate',                      phys_decomp)   
       call addfld ('HDEPTH  ','km    ',1,    'A','Heating Depth',                                phys_decomp)
       call addfld ('BEMF  ','Pa    ',pver,    'A','Beres Eastward MF',                                phys_decomp)
       call addfld ('BWMF  ','Pa    ',pver,    'A','Beres Westward MF',                                phys_decomp)
       call addfld ('BNMF  ','Pa    ',pver,    'A','Beres Northward MF',                                phys_decomp)
       call addfld ('BSMF  ','Pa    ',pver,    'A','Beres Southward MF',                                phys_decomp)
       call addfld ('EMF  ','Pa    ',pver,    'A','Eastward MF',                                phys_decomp)
       call addfld ('WMF  ','Pa    ',pver,    'A','Westward MF',                                phys_decomp)
       call addfld ('NMF  ','Pa    ',pver,    'A','Northward MF',                                phys_decomp)
       call addfld ('SMF  ','Pa    ',pver,    'A','Southward MF',                                phys_decomp)

       call add_default ('UTEND1   ', 1, ' ')
       call add_default ('UTEND2   ', 1, ' ')
       call add_default ('UTEND3   ', 1, ' ')
       call add_default ('UTEND4   ', 1, ' ')
       call add_default ('UTEND5   ', 1, ' ')


       call add_default ('BUTEND1   ', 1, ' ')
       call add_default ('BUTEND2   ', 1, ' ')
       call add_default ('BUTEND3   ', 1, ' ')
       call add_default ('BUTEND4   ', 1, ' ')
       call add_default ('BUTEND5   ', 1, ' ')

       call add_default ('NETDT    ', 1, ' ')
       call add_default ('HDEPTH   ', 1, ' ')
       call add_default ('MAXQ0    ', 1, ' ')

    end if

    call addfld ('TTGWSDF' ,'K/s   ',pver, 'A','t tendency - gw spec: diffusion term',phys_decomp)
    call addfld ('TTGWSKE' ,'K/s   ',pver, 'A','t tendency - gw spec: kinetic energy conversion term',phys_decomp)
!   call add_default ('TTGWSDF', 1, ' ')
!   call add_default ('TTGWSKE', 1, ' ')

    call addfld ('EKGWSPEC' ,'M2/S   ',pverp, 'A','effective Kzz due to gw spectrum',phys_decomp)
    call add_default ('EKGWSPEC', 1, ' ')

    if (pgwv > 0) then
       call addfld ('TAUE' ,'Pa   ',pverp, 'A','CM Eastward Reynolds stress',phys_decomp)
       call addfld ('TAUW' ,'Pa   ',pverp, 'A','CM Westward Reynolds stress',phys_decomp)
       call addfld ('TAUNET' ,'Pa   ',pverp, 'A','CM E+W Reynolds stress',phys_decomp)
       call addfld ('TAUN' ,'Pa   ',pverp, 'A','CM Northward Reynolds stress',phys_decomp)
       call addfld ('TAUS' ,'Pa   ',pverp, 'A','CM Southward Reynolds stress',phys_decomp)
       call add_default ('TAUE', 1, ' ')
       call add_default ('TAUW', 1, ' ')
       call add_default ('TAUNET', 1, ' ')
       call add_default ('TAUN', 1, ' ')
       call add_default ('TAUS', 1, ' ')

       call addfld ('BTAUE' ,'Pa   ',pverp, 'A','Beres Eastward Reynolds stress',phys_decomp)
       call addfld ('BTAUW' ,'Pa   ',pverp, 'A','Beres Westward Reynolds stress',phys_decomp)
       call addfld ('BTAUN' ,'Pa   ',pverp, 'A','Beres Northward Reynolds stress',phys_decomp)
       call addfld ('BTAUS' ,'Pa   ',pverp, 'A','Beres Southward Reynolds stress',phys_decomp)
       call add_default ('BTAUE', 1, ' ')
       call add_default ('BTAUW', 1, ' ')
       call add_default ('BTAUN', 1, ' ')
       call add_default ('BTAUS', 1, ' ')

!+Mods for CM
       call addfld ('FRONTGF', 'K^2/M^2/S', pver, 'A', 'Frontogenesis function at gws src level', phys_decomp)
       call add_default ('FRONTGF', 1, ' ')
       call addfld ('FRONTGFA', 'K^2/M^2/S', pver, 'A', 'Frontogenesis function at gws src level', phys_decomp)
       call add_default ('FRONTGFA', 1, ' ')
!=Mods for CM



! Gravity wave source spectra by wave number

! Positive phase velocities
      l=-1
      do m=0,PGWV
         ccn=float(m)*dc2
         if (ccn <= cref2(pgwv)) then
            l=l+1
            write (fnum,fmt='(i3)') 100+l
            dumc1x='TAUXSp'//fnum(2:3)
            dumc1y='TAUYSp'//fnum(2:3)
            write (fnum,fmt='(f6.2,a5,f6.2,a4)') ccn,' m/s'
            dumc2='tau at c= '//fnum
	    call addfld (dumc1x(1:8),'Pa   ',pver, 'A',dumc2,phys_decomp)
            call addfld (dumc1y(1:8),'Pa   ',pver, 'A',dumc2,phys_decomp)

         end if
      end do


! Positive phase velocities
      l=-1
      do m=0,PGWV
         ccn=float(m)*dc
         if (ccn <= cref(pgwv)) then
            l=l+1
            write (fnum,fmt='(i3)') 100+l
            dumc1x='BTAUXSp'//fnum(2:3)
            dumc1y='BTAUYSp'//fnum(2:3)
            write (fnum,fmt='(f6.2,a5,f6.2,a4)') ccn,' m/s'
            dumc2='Beres tau at c= '//fnum
            call addfld (dumc1x(1:9),'Pa   ',pver, 'A',dumc2,phys_decomp)
            call addfld (dumc1y(1:9),'Pa   ',pver, 'A',dumc2,phys_decomp)

         end if
      end do


! Negative phase velocities
      l=0
      do m=-1,-PGWV,-1
         ccn=float(m)*dc2
         if (ccn >= cref2(-pgwv)) then
            l=l+1
            write (fnum,fmt='(i3)') 100+l
            dumc1x='TAUXSn'//fnum(2:3)
            dumc1y='TAUYSn'//fnum(2:3)
            write (fnum,fmt='(f6.2,a5,f6.2,a4)') ccn,' m/s'
            dumc2='tau at c= '//fnum
            call addfld (dumc1x(1:8),'Pa   ',pver, 'A',dumc2,phys_decomp)
            call addfld (dumc1y(1:8),'Pa   ',pver, 'A',dumc2,phys_decomp)

         end if
      end do



! Negative phase velocities
      l=0
      do m=-1,-PGWV,-1
         ccn=float(m)*dc
         if (ccn >= cref(-pgwv)) then
            l=l+1
            write (fnum,fmt='(i3)') 100+l
            dumc1x='BTAUXSn'//fnum(2:3)
            dumc1y='BTAUYSn'//fnum(2:3)
            write (fnum,fmt='(f6.2,a5,f6.2,a4)') ccn,' m/s'
            dumc2='Beres tau at c= '//fnum
            call addfld (dumc1x(1:9),'Pa   ',pver, 'A',dumc2,phys_decomp)
            call addfld (dumc1y(1:9),'Pa   ',pver, 'A',dumc2,phys_decomp)

         end if
      end do



    endif

!+Mods for Beres04
! Initialization of BERES' parameterization parameters
    call gw_inti_beres
!-Mods for Beres04

  end  subroutine gw_inti

!+Mods for Beres04
!==============================================================================
  subroutine gw_inti_beres

   use ioFileMod,        only: getfil
#if ( defined SPMD )
   use mpishorthand
#endif

   use netcdf

! VARIABLES NEEDED TO READ IN TABLE OF SPECTRA FROM FILE

   integer ih, iw,ip                ! height,wind, wave look-up table indices
   integer ncid, varid, ncstat       ! netcdf file id, variable id, result of netcdf call
   real(r8) :: tmpmfcc(20,-40:40,-40:40)  ! stores source spectra
   character(len=256) :: gw_drag_file_loc ! local filepath of gw_drag_file

   !----------------------------------------------------------------------
   ! read in look-up table for source spectra
   !-----------------------------------------------------------------------

   if (masterproc) then
      
     call getfil(gw_drag_file, gw_drag_file_loc)
     ncstat = NF90_OPEN (gw_drag_file_loc,0,ncid)

      if (ncstat .ne. 0) then
         write(iulog,*) 'Error reading in netcdf file ',gw_drag_file,'.  ', NF90_STRERROR(ncstat)
         write(iulog,*) 'Be sure that the variable gw_drag_file is set properly in gw_drag.F90.'
         call endrun
      endif

      ncstat = NF90_INQ_VARID (ncid,'mfcc',varid)

      if (ncstat .ne. 0) then
         write(iulog,*) 'Error reading data from ',gw_drag_file,'.  ', NF90_STRERROR(ncstat)
         call endrun
      endif

      ncstat = NF90_GET_VAR(ncid,varid,tmpmfcc)
      ncstat = NF90_CLOSE (ncid)

      mfcc(:,:,:) = tmpmfcc(:,:,-PGWV:PGWV)   ! reduce the amount of wave numbers needed

      write(iulog,*) 'Read-in source spectra from file'
      write (iulog,*) 'MFCC=',maxval(mfcc),minval(mfcc)

   endif
 !     Broadcast results
#ifdef SPMD
   call mpibcast (mfcc, MAXH*(2*PGWV+1)*(2*MAXUH+1),mpir8,0,mpicom)
#endif

 end subroutine gw_inti_beres

!-Mods for Beres04

!===============================================================================

  subroutine set_phase_speeds (lchnk, ncol, ngwv, ubi , c, bflag ) 

    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    integer, intent(in) :: ngwv

    logical, intent(in) :: bflag

    real(r8), intent(in) :: ubi(pcols,0:pver)
    
    real(r8), intent(out) :: c(pcols,-pgwv:pgwv)

!------------------------------------------------------------------------------
    integer l, i, n
    
    c(:,:) = 0._r8
    ! Check if dealing with the spectrum or the orographic component
    if ( ngwv /= 0) then

       ! Now check if dealing with the Beres04 component of the spectrum
       if (bflag) then

          do l=-ngwv,ngwv
             do i =1, ncol
                c(i,l) = cref(l)
             end do
          end do

       else

          do l=-ngwv,ngwv
             do i =1, ncol
                c(i,l) = abs(ubi(i,kbotbg)) + cref(l)
             end do
          end do

       end if

    else

       do i=1, ncol
          c(i,0) = 0._r8
       end do

    end if


  end subroutine set_phase_speeds

  

!===============================================================================

  subroutine gw_intr(state, sgh, pbuf, dt, ptend, landfrac)

!-----------------------------------------------------------------------
! Interface for multiple gravity wave drag parameterization.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field
    use time_manager,   only: get_nstep,get_step_size
    use physconst,      only: cpairv            ! Constituent dependent physconst variables (WACCM-X)
    use shr_assert_mod, only: shr_assert

    real(r8), intent(in) :: sgh(pcols)            ! standard deviation of orography
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(in) :: dt                    ! time step
    real(r8), intent(in) :: landfrac(pcols)       ! Land fraction

    type(physics_state), intent(in) :: state      ! physics state structure
    type(physics_ptend), intent(out):: ptend      ! parameterization tendency structure

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

    character(len=8) :: fnum,fname                ! dummy characters
    integer  ic,l,m                               ! dummy integers
    real(r8) qtgw(pcols,pver,pcnst)               ! constituents tendencies

    real(r8) :: taue(pcols,0:pver)                ! Reynolds stress for eastward propagating waves
    real(r8) :: tauw(pcols,0:pver)                ! Reynolds stress for westward propagating waves
    real(r8) :: taun(pcols,0:pver)                ! Reynolds stress for northward propagating waves
    real(r8) :: taus(pcols,0:pver)                ! Reynolds stress for southward propagating waves

    real(r8) :: btaue(pcols,0:pver)               ! Reynolds stress for eastward propagating waves Beres scheme
    real(r8) :: btauw(pcols,0:pver)               ! Reynolds stress for westward propagating waves Beres scheme
    real(r8) :: btaun(pcols,0:pver)               ! Reynolds stress for northward propagating waves Beres scheme
    real(r8) :: btaus(pcols,0:pver)               ! Reynolds stress for southward propagating waves Beres scheme

    real(r8) :: c(pcols,-pgwv:pgwv)                       ! spectrum phase speeds for each column

    real(r8), pointer :: pblh(:)
    real(r8), pointer :: kvt_in(:,:)
    real(r8) :: kvtt(pcols,0:pver)                ! molecular diffusivity 

! C.-C. Chen
    integer  :: src_level(pcols)                  ! index of gravity wave source

    real(r8) :: egwdffi(pcols,0:pver)             ! effective gw diffusivity at interfaces needed for output
    real(r8) :: egwdffi_tot(pcols,0:pver)         ! sum from the two types of spectral GW

    logical  :: lq(pcnst)

!-----------------------------------------------------------------------------
    ! Make sure frontogenesis has actually been set.
    call shr_assert(allocated(state%frontgf) .and. allocated(state%frontga), &
         msg="Dycore has not allocated frontogenesis fields before they are &
         &needed in gw_drag!")

    lchnk = state%lchnk
    ncol  = state%ncol

    lq(:) = .true.
    call physics_ptend_init(ptend, state%psetcols, "Gravity waves drag", ls=.true., lu=.true., lv=.true., lq=lq)

    taue(:,:) = 0._r8
    tauw(:,:) = 0._r8
    taun(:,:) = 0._r8
    taus(:,:) = 0._r8

    btaue(:,:) = 0._r8
    btauw(:,:) = 0._r8
    btaun(:,:) = 0._r8
    btaus(:,:) = 0._r8

! Profiles of background state variables
    call gw_prof(lchnk, ncol, &
         state%u   , state%v   , state%t   , state%pmid   , state%pint, &
         rhoi      , ni        , ti        , nm           )

!--------------------------------------------------------
! Initialize and calculate local molecular diffusivity
!--------------------------------------------------------

    call pbuf_get_field(pbuf, kvt_idx, kvt_in)  ! kvt_in(1:pcols,1:pverp)

    kvtt(:ncol,0:pver) = kvt_in(:ncol,:)

    do k = 2, nbot_molec
      kvtt(:ncol,k-1) = kvt_in(:ncol,k) / (cpairv(:ncol,k,lchnk)+cpairv(:ncol,k-1,lchnk)) * 2._r8
    enddo
    kvtt(:ncol,0) = kvt_in(:ncol,1) / (1.5_r8*cpairv(:ncol,1,lchnk)-0.5_r8*cpairv(:ncol,2,lchnk)) 		   

!-----------------------------------------------------------------------------
! Non-orographic background gravity wave spectra
!-----------------------------------------------------------------------------
    if (pgwv >0) then

!+Mods for Beres04
! Determine wave sources for Beres04 scheme
! set phase speeds for this chunk
       call set_phase_speeds (lchnk, ncol, PGWV,  ubi , c , .true.)

       call gw_bgnd_beres (lchnk          , ncol       ,                     &
            state%u    , state%v    , state%t    , state%pmid , state%pint , &
            state%pdel , state%rpdel, state%lnpint,kldv       , kldvmn     , &
            ksrc       , ksrcmn     , rdpldv     , tau        , ubi        , &
            ubm        , xv         , yv         , PGWV       , &
            state%zm, src_level, pbuf)

! C.-C. Chen
!       src_level(:) = bkbotbg

! Solve for the drag profile with Beres source spectrum
       call gw_drag_prof (lchnk     , ncol       ,                     &
            PGWV       , bkbotbg     , ktopbg     , state%u    , state%v     , &
            state%t    , state%pint , state%pdel , state%rpdel, state%lnpint ,&
            rhoi       , ni         , ti         , nm         , dt          , &
            kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv      , &
            tau        , ubi        , ubm        , xv         , yv          , &
            effgw_beres , utgw       , vtgw       , tau0x      , tau0y      , &
            kvtt       , state%pmid , state%q    , state%s    , ttgw       , qtgw, .true.     , &
            taue       , tauw       , taun       , taus       , c , &
            egwdffi    , src_level)
 
!  add the diffusion coefficients
       do k = 0, pver
          do i = 1, ncol
             egwdffi_tot(i,k) = egwdffi(i,k) 
          end do
       end do

! Store constituents tendencies
!caf       do m=1, pcnst+pnats
       do m=1, pcnst
          do k = 1, pver
             do i = 1 , ncol
                ptend%q(i,k,m) = qtgw(i,k,m)
             end do
          end do

       end do                          ! another constituent

! add the momentum tendencies to the output tendency arrays
       do k = 1, pver
          do i = 1, ncol
             ptend%u(i,k) = utgw(i,k)
             ptend%v(i,k) = vtgw(i,k)
             ptend%s(i,k) = ttgw(i,k)
          end do
       end do


! C.-C. Chen, momentum & energy conservation
       call momentum_energy_conservation(ncol, taue, tauw, taun, taus, &
            state%pint, state%pdel, state%t, state%u, state%v, ptend%u, ptend%v, ptend%s, &
            utgw, vtgw, ttgw, src_level, dt)

! Write output fields to history file
       call outfld ('BUTGWSPEC', utgw , pcols, lchnk)
       call outfld ('BVTGWSPEC', vtgw , pcols, lchnk)
       call outfld ('BTAUE', taue , pcols, lchnk)
       call outfld ('BTAUW', tauw , pcols, lchnk)
       call outfld ('BTAUN', taun , pcols, lchnk)
       call outfld ('BTAUS', taus , pcols, lchnk)


!-Mods for Beres04

! Determine the wave source for C&M background spectrum 

       call gw_bgnd (lchnk          , ncol       ,                           &
            state%u    , state%v    , state%t    , state%pmid , state%pint , &
            state%pdel , state%rpdel, state%lnpint,kldv       , kldvmn     , &
            ksrc       , ksrcmn     , rdpldv     , tau        , ubi        , &
            ubm        , xv         , yv         , PGWV       , kbotbg     , &
            state%frontgf      ,state%frontga     )

! set phase speeds for this chunk
       call set_phase_speeds (lchnk, ncol, PGWV,  ubi , c, .false. )

! C.-C. Chen
       src_level(:) = kbotbg

! Solve for the C&M drag profile
       call gw_drag_prof (lchnk     , ncol       ,                           &
            PGWV       , kbotbg     , ktopbg     , state%u    , state%v    , &
            state%t    , state%pint , state%pdel , state%rpdel, state%lnpint,&
            rhoi       , ni         , ti         , nm         , dt         , &
            kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv     , &
            tau        , ubi        , ubm        , xv         , yv         , &
            effgw_spec , utgw       , vtgw       , tau0x      , tau0y      , &
            kvtt       ,state%pmid  , state%q    , state%s    , ttgw       , qtgw , .false.  ,&
            taue       , tauw       , taun       , taus       ,c          , &
            egwdffi    , src_level)

!  add the diffusion coefficients
       do k = 0, pver
          do i = 1, ncol
             egwdffi_tot(i,k) = egwdffi_tot(i,k) + egwdffi(i,k) 
          end do
       end do


!Add the constituent tendencies
       do m=1, pcnst

          do k = 1, pver
             do i = 1 , ncol
                ptend%q(i,k,m) = ptend%q(i,k,m) + qtgw(i,k,m) 
             end do
          end do

       end do                          ! another constituent


! add the momentum tendencies to the output tendency arrays
       do k = 1, pver
          do i = 1, ncol
             ptend%u(i,k) = ptend%u(i,k) + utgw(i,k) 
             ptend%v(i,k) = ptend%v(i,k) + vtgw(i,k) 
             ptend%s(i,k) = ptend%s(i,k) + ttgw(i,k)
          end do
       end do

! C.-C. Chen, momentum & energy conservation
       call momentum_energy_conservation(ncol, taue, tauw, taun, taus, &
            state%pint, state%pdel, state%t, state%u, state%v, ptend%u, ptend%v, ptend%s, &
            utgw, vtgw, ttgw, src_level, dt)

! Write output fields to history file
       call outfld ('UTGWSPEC', utgw , pcols, lchnk)
       call outfld ('VTGWSPEC', vtgw , pcols, lchnk)
       call outfld ('TAUE', taue , pcols, lchnk)
       call outfld ('TAUW', tauw , pcols, lchnk)
       call outfld ('TAUNET', taue+tauw , pcols, lchnk)
       call outfld ('TAUN', taun , pcols, lchnk)
       call outfld ('TAUS', taus , pcols, lchnk)
       call outfld ('EKGWSPEC', egwdffi_tot , pcols, lchnk)  

    else

! zero net tendencies if no spectrum computed
       ptend%u = 0._r8
       ptend%v = 0._r8
       ptend%s = 0._r8
       ptend%q = 0._r8

    end if
!-----------------------------------------------------------------------------
! Orographic stationary gravity wave
!-----------------------------------------------------------------------------
    call pbuf_get_field(pbuf, pblh_idx, pblh)
! Determine the orographic wave source

    call gw_oro (lchnk, ncol,                                             &
         state%u    , state%v    , state%t    , sgh        , state%pmid , &
         state%pint , state%pdel , state%zm   , nm         , pblh       , &
         kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv     , &
         tau        , ubi        , ubm        , xv         , yv         )


! set phase speeds for this chunk
    call set_phase_speeds (lchnk, ncol, 0,  ubi , c, .false.)

! C.-C. Chen
       src_level(:) = kbotoro

! Solve for the orographic drag profile
       call gw_drag_prof (lchnk     , ncol       ,                           &
            0          , kbotoro    , ktoporo    , state%u    , state%v    , &
            state%t    , state%pint , state%pdel , state%rpdel, state%lnpint,&
            rhoi       , ni         , ti         , nm         , dt         , &
            kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv     , &
            tau        , ubi        , ubm        , xv         , yv         , &
            effgw_oro  , utgw       , vtgw       , tau0x      , tau0y      , &
            kvtt       ,state%pmid  , state%q    , state%s    , ttgw       , qtgw , .false.  ,&
            taue       , tauw       , taun       , taus       , c , &
            egwdffi    , src_level)



! Add the orographic tendencies to the spectrum tendencies
! Compute the temperature tendency from energy conservation (includes spectrum).
    do k = 1, pver
       do i = 1, ncol
          ptend%u(i,k) = ptend%u(i,k) + utgw(i,k) * landfrac(i)
          ptend%v(i,k) = ptend%v(i,k) + vtgw(i,k) * landfrac(i)
          ptend%s(i,k) = ptend%s(i,k) + ttgw(i,k)                                 &
                       -(ptend%u(i,k) * (state%u(i,k) + ptend%u(i,k)*0.5_r8*dt)   &
                       + ptend%v(i,k) * (state%v(i,k) + ptend%v(i,k)*0.5_r8*dt))
       end do
    end do

! Load TTGW with the total heating for output
    do k=1, pver
       do i=1, ncol
          ttgw(i,k) = ptend%s(i,k) 
       end do
    end do

    do m=1, pcnst

       do k = 1, pver
          do i = 1 , ncol
             ptend%q(i,k,m) = ptend%q(i,k,m) + qtgw(i,k,m) 
          end do
       end do

    end do                          ! another constituent

! Write output fields to history file
    call outfld ('UTGWORO', utgw,  pcols, lchnk)
    call outfld ('VTGWORO', vtgw,  pcols, lchnk)
    call outfld ('TTGW', ttgw / cpair,  pcols, lchnk)
    call outfld ('TAUGWX',  tau0x, pcols, lchnk)
    call outfld ('TAUGWY',  tau0y, pcols, lchnk)
    call outfld ('SGH    ', sgh,   pcols, lchnk)

    return
  end  subroutine gw_intr

!===============================================================================
       subroutine momentum_energy_conservation(ncol,taue, tauw, taun, taus, &
            pi, pdel, t, u, v, dudt, dvdt, dsdt, &
            utgw, vtgw, ttgw, src_level, dt)

! C.-C. Chen, momentum & energy conservation
       implicit none

       integer, intent(in) :: ncol
       real(r8), intent(in) :: taue(pcols,0:pver)
       real(r8), intent(in) :: tauw(pcols,0:pver)
       real(r8), intent(in) :: taun(pcols,0:pver)
       real(r8), intent(in) :: taus(pcols,0:pver)
       real(r8), intent(in) :: pi(pcols,1:pver)
       real(r8), intent(in) :: pdel(pcols,pver)
       real(r8), intent(in) :: t(pcols,pver)
       real(r8), intent(in) :: u(pcols,pver)
       real(r8), intent(in) :: v(pcols,pver)
       real(r8), intent(inout) :: dudt(pcols,pver)
       real(r8), intent(inout) :: dvdt(pcols,pver)
       real(r8), intent(inout) :: dsdt(pcols,pver)
       real(r8), intent(inout) :: utgw(pcols,pver)
       real(r8), intent(inout) :: vtgw(pcols,pver)
       real(r8), intent(inout) :: ttgw(pcols,pver)
       integer, intent(in) :: src_level(pcols)
       real(r8), intent(in) :: dt
! local variables
       integer :: i,l,k
       real(r8) :: dz, dE, ut_dz, vt_dz

       do i = 1, ncol

! total mass from ground to source level, rho*dz = dp/g
        dz = 0.0_r8
        do k=pver,src_level(i)+1,-1
         dz=dz+pdel(i,k)/g
        enddo 

! tendency for U & V below source level       
        do k=src_level(i)+1,pver
         ut_dz = -(taue(i,src_level(i))+tauw(i,src_level(i)))/dz
         vt_dz = -(taun(i,src_level(i))+taus(i,src_level(i)))/dz
         dudt(i,k)=dudt(i,k)+ut_dz
         dvdt(i,k)=dvdt(i,k)+vt_dz

         utgw(i,k) = utgw(i,k)+ut_dz
         vtgw(i,k) = vtgw(i,k)+vt_dz
        enddo

! net gain/loss of total energy in the column
        dE = 0.0_r8
        do k=1,pver
         dE=dE+pdel(i,k)*(dsdt(i,k)+dudt(i,k)*(u(i,k)+dudt(i,k)*0.5_r8*dt)+  &
                                    dvdt(i,k)*(v(i,k)+dvdt(i,k)*0.5_r8*dt) )
        enddo

        dE = dE/(pi(i,pver)-pi(i,src_level(i)))

! subtract net gain/loss of total energy below source level
        do k=src_level(i)+1,pver
         dsdt(i,k) = dsdt(i,k)-dE
         ttgw(i,k) = ttgw(i,k)-dE
        enddo
       enddo
             
       return 
       end subroutine momentum_energy_conservation

!===============================================================================
  subroutine gw_prof (lchnk, ncol, u, v, t, pm, pi, rhoi, ni, ti, nm )

    use phys_grid, only: get_lat_p, get_lon_p, get_rlat_p,get_rlon_p
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
    real(r8), intent(out) :: nm(pcols,pver)       ! midpoint brunt-vaisalla frequency

!---------------------------Local storage-------------------------------
    integer :: i,k                                ! loop indexes

    real(r8) :: dtdp
    real(r8) :: n2                                ! Brunt-Vaisalla frequency squared
    real(r8) :: ti1(pcols)                        ! Temperature at top interface
    real(r8) :: tm0(pcols)                        ! Temperature outside the model domain
    real(r8) :: pm0(pcols)                        ! pressure at which tm0 is provided (in hPa)
    real(r8) :: lat(pcols)                        ! lat indexes for current column
    real(r8) :: lon(pcols)                        ! lon indexes for current column

!-----------------------------------------------------------------------------
! Determine the interface densities and Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------

! The top interface values are calculated assuming an isothermal atmosphere 
! above the top level.
    ti1(:ncol) = t(:ncol,1)
    k = 0
    do i = 1, ncol
       ti(i,k)   = ti1(i)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k)   = sqrt (g*g / (cpair*ti(i,k)))
    end do

! Interior points use centered differences
    do k = 1, pver-1
       do i = 1, ncol
          ti(i,k)   = 0.5_r8 * (t(i,k) + t(i,k+1))
          rhoi(i,k) = pi(i,k) / (r*ti(i,k))
          dtdp      = (t(i,k+1) - t(i,k)) / (pm(i,k+1) - pm(i,k))
          n2        = g*g/ti(i,k) * (1._r8/cpair - rhoi(i,k)*dtdp)
          ni(i,k)   = sqrt (max (n2min, n2))
       end do
    end do

! Bottom interface uses bottom level temperature, density; next interface
! B-V frequency.
    k = pver
    do i = 1, ncol
       ti(i,k)   = t(i,k)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k)   = ni(i,k-1)
    end do

!-----------------------------------------------------------------------------
! Determine the midpoint Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------
    do k = 1,pver
       do i = 1,ncol
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

       ubi(i,pver) = sqrt (usrc(i)**2 + vsrc(i)**2)
       if (ubi(i,pver) /= 0._r8) then
          xv(i) = usrc(i) / ubi(i,pver)
          yv(i) = vsrc(i) / ubi(i,pver)
       else ! Magnitude 0 vector.
          xv(i) = 0._r8
          yv(i) = 0._r8
       end if
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
    pil = acos( -1._r8 )
    do i = 1, ncol
       lzsrc   = max(2._r8 * pil * usrc(i) / nsrc(i), lzmax)
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
       ngwv, kbot, frontgf, frontga)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
! 
! The parameterization is assumed to operate only where water vapor 
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------

    use phys_grid,     only: get_rlat_all_p
    use physconst,     only: pi_g => pi
!   use time_manager, only: get_nstep,get_step_size,get_curr_calday	! rrg 1/7/2005

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

    real(r8), intent(in) :: frontgf  (pcols,pver)        ! frontogenesis function
    real(r8), intent(in) :: frontga  (pcols,pver)        ! frontogenesis angle

    integer, intent(out) :: kldv(pcols)           ! top interface of low level stress divergence region
    integer, intent(out) :: kldvmn                ! min value of kldv
    integer, intent(out) :: ksrc(pcols)           ! index of top interface of source region
    integer, intent(out) :: ksrcmn                ! min value of ksrc

    real(r8), intent(in) :: rdpldv(pcols)        ! 1/dp across low level divergence region
    real(r8), intent(out) :: tau(pcols,-ngwv:ngwv,0:pver)! wave Reynolds stress
    real(r8), intent(out) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
    real(r8), intent(out) :: ubm(pcols,pver)      ! projection of wind at midpoints
    real(r8), intent(out) :: xv(pcols)            ! unit vectors of source wind (x)
    real(r8), intent(out) :: yv(pcols)            ! unit vectors of source wind (y)

!---------------------------Local storage-------------------------------
    integer :: i,k,l,m                            ! loop indexes

    real(r8) :: rlat(pcols)                       ! latitude in radians for columns
    real(r8) :: tauback(pcols)                    ! background stress at c=0
    real(r8) :: flat_gw                           ! The actual lat dependence of GW spec.

    real(r8) :: usrc(pcols)                       ! u wind averaged over source region
    real(r8) :: vsrc(pcols)                       ! v wind averaged over source region

    real(r8) :: xfrontg(pcols,pver) 

!---------------------------------------------------------------------------
! Determine the source layer wind and unit vectors, then project winds.
!---------------------------------------------------------------------------

! Just use the source level interface values for the source
! wind speed and direction (unit vector).
    do i = 1, ncol
       ksrc(i) = kbot
       kldv(i) = kbot
       usrc(i) = 0.5_r8*(u(i,kbot+1)+u(i,kbot))
       vsrc(i) = 0.5_r8*(v(i,kbot+1)+v(i,kbot))
       ubi(i,kbot) = sqrt (usrc(i)**2 + vsrc(i)**2)
       if (ubi(i,kbot) /= 0._r8) then
          xv(i) = usrc(i) / ubi(i,kbot)
          yv(i) = vsrc(i) / ubi(i,kbot)
       else ! Magnitude 0 vector.
          xv(i) = 0._r8
          yv(i) = 0._r8
       end if
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

   xfrontg(:,:) = 0._r8

   do i=1, ncol
      do k=1, pver
         xfrontg(i,k) = frontgf(i,k) 
      end do
   end do

   call outfld ('FRONTGF', xfrontg, pcols, lchnk)
   call outfld ('FRONTGFA', frontga, pcols, lchnk)

!-----------------------------------------------------------------------
! Gravity wave sources
!-----------------------------------------------------------------------

    tau(:,:,:) = 0._r8

! Determine the background stress at c=0
    do i=1,ncol
      tauback(i) = taubgnd * tauscal
    enddo

    do l = 1, ngwv
       do i = 1, ncol

          !
          ! GW generation depends  on frontogenesis at src level
          ! rrg: use frontgf 2 levels below launch (~600 mb) (note: hardwiring needs to be removed)
          !
          if ( xfrontg(i,kbot+2) > frontgfc) then	
             flat_gw = 1.0_r8
          else
             flat_gw = 0.0_r8
          end if
          tau(i, l,kbot) = tauback(i) * fav(l) * flat_gw
          tau(i,-l,kbot) = tau(i,l,kbot)

       end do
    end do

    do i = 1, ncol
       tau(i,0,kbot) = 0._r8
    end do

    
! Determine the min value of kldv and ksrc for limiting later loops
! and the pressure at the top interface of the low level stress divergence
! region.

    ksrcmn = pver
    kldvmn = pver

    return

  end  subroutine gw_bgnd

!===============================================================================
  subroutine gw_drag_prof (lchnk, ncol,                                  &
       ngwv,   kbot,   ktop,   u,      v,       t,                       &
       pi,     dpm,    rdpm,   piln,   rhoi,    ni,   ti,    nm,   dt,   &
       kldv,   kldvmn, ksrc,   ksrcmn, rdpldv,  tau,  ubi,   ubm,  xv,   &
       yv,     effgw,  ut,     vt,     tau0x,   tau0y,                   &
       kvtt,   pm,     q,      dse,    ttgw,    qtgw, bflag ,            &
       taue,   tauw,   taun,   taus,   c,       egwdffi,       src_level)
    use phys_grid,     only: get_rlat_all_p

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

    real(r8), intent(in) ::  kvtt(pcols,0:pver)                ! molecular thermal diffusivity
    real(r8), intent(in) ::  pm(pcols,pver)                    ! pressure at midpoints
    real(r8), intent(in) ::  q(pcols,pver,pcnst)               ! constituent array
    real(r8), intent(in) ::  dse(pcols,pver)                   ! dry static energy
    real(r8), intent(in)  :: c(pcols,-pgwv:pgwv)                       ! wave phase speeds for each column
    real(r8), intent(out) :: ttgw(pcols,pver)                  ! GW's heating
    real(r8), intent(out) :: qtgw(pcols,pver,pcnst)            ! GW's constituent tendency
    real(r8), intent(out) :: taue(pcols,0:pver)                ! Reynolds stress for eastward propagating waves
    real(r8), intent(out) :: tauw(pcols,0:pver)                ! Reynolds stress for westward propagating waves
    real(r8), intent(out) :: taun(pcols,0:pver)                ! Reynolds stress for northward propagating waves
    real(r8), intent(out) :: taus(pcols,0:pver)                ! Reynolds stress for southward propagating waves
    real(r8), intent(out) :: egwdffi(pcols,0:pver)             ! effective gw diffusivity at interfaces
                                                               ! egwdffi is now returned 
                                                               ! to the calling progam for output
! C.-C. Chen
    integer, intent(in) :: src_level(pcols)

!+Mods for Beres04
    logical, intent(in) :: bflag
!-Mods for Beres04

!---------------------------local storage-------------------------------
    integer :: i,k,l,m                            ! loop indexes

    real(r8) :: d(pcols)                          ! "total" diffusivity 
    real(r8) :: dsat(pcols,-ngwv:ngwv)            ! saturation diffusivity
    real(r8) :: dscal                             ! fraction of dsat to use
    real(r8) :: mi                                ! imaginary part of vertical wavenumber
    real(r8) :: taudmp                            ! stress after damping
    real(r8) :: taumax(pcols)                     ! max(tau) for any l
    real(r8) :: tausat(pcols,-ngwv:ngwv)          ! saturation stress
    real(r8) :: ubmc(pcols,-ngwv:ngwv)            ! (ub-c)
    real(r8) :: ubmc2                             ! (ub-c)**2
    real(r8) :: ubt(pcols,pver)                   ! ubar tendency
    real(r8) :: ubtl                              ! ubar tendency from wave l
    real(r8) :: ubtlsat                           ! saturation tendency
    real(r8) :: wrk
    real(r8) :: gwut(pcols,pver,-ngwv:ngwv)       ! gravity wave wind tendency for each wave

! C.-C. Chen
    real(r8) :: tausg(pcols,-ngwv:ngwv,0:pver)            ! signed wave Reynolds stress 

    real(r8) :: taulu(pcols,0:pver)                      ! Reynolds stress for waves propagating to the left
                                                  ! of the source zonal winds
    real(r8) :: tauru(pcols,0:pver)                      ! Reynolds stress for waves propagating to the right
                                                  ! of the source zonal winds
    real(r8) :: taulv(pcols,0:pver)                      ! Reynolds stress for waves propagating to the left
                                                  ! of the source meridional winds
    real(r8) :: taurv(pcols,0:pver)                      ! Reynolds stress for waves propagating to the right
                                                  ! of the source meridional winds
    real(r8) :: ut1(pcols,pver)                   ! zonal wind tendency for  c < -40
    real(r8) :: ut2(pcols,pver)                   ! zonal wind tendency for  -40 < c < -15
    real(r8) :: ut3(pcols,pver)                   ! zonal wind tendency for  -15 < c <  15
    real(r8) :: ut4(pcols,pver)                   ! zonal wind tendency for   15 < c <  40
    real(r8) :: ut5(pcols,pver)                   ! zonal wind tendency for   40 < c


    real(r8) :: EMF(pcols,pver)
    real(r8) :: WMF(pcols,pver)
    real(r8) :: NMF(pcols,pver)
    real(r8) :: SMF(pcols,pver)

    real(r8) :: cx,cy                             ! project wave phase speeds

    real(r8) :: rlat(pcols)                       ! latitude of the columns in the chunck (radians)

    real(r8) :: ptaper(pcols)                     ! rrg: polar taper
    real(r8) :: al0				  ! rrg: used in lat dependence of the taper
    real(r8) :: dlat0                             ! rrg: used in lat dependence of the taper


    real(r8) :: dummy(pcols,pver)
    real(r8) :: dummyx(pcols,pver)
    real(r8) :: dummyy(pcols,pver)
    real(r8) :: ccn
    character(len=8) :: fnum
    character(len=80) dumc1, dumc1x, dumc1y          ! dummies

    real(r8):: taux(pcols,-pgwv:pgwv,0:pver)      ! wave stress in zonal direction
    real(r8):: tauy(pcols,-pgwv:pgwv,0:pver)      ! wave stress in meridional direction

    integer  :: ix, iy                            ! dummies
    integer  :: usi,vsi                           ! dummies
    character(len=80) :: dum1,dum2,dum3,dum4      ! dummies
 
! Initialize gravity wave drag tendencies to zero

!+Mods for Beres04
    ut1(:,:) = 0._r8
    ut2(:,:) = 0._r8
    ut3(:,:) = 0._r8
    ut4(:,:) = 0._r8
    ut5(:,:) = 0._r8
 

!-Mods for Beres04

    call get_rlat_all_p(lchnk, ncol, rlat)

    do k = 1,pver
       do i = 1,pcols
          ut(i,k) = 0._r8
          vt(i,k) = 0._r8
       end do
    end do

    taue(:,:) = 0._r8
    tauw(:,:) = 0._r8
    taun(:,:) = 0._r8
    taus(:,:) = 0._r8

    EMF(:,:)=0._r8
    WMF(:,:)=0._r8
    SMF(:,:)=0._r8
    NMF(:,:)=0._r8


    gwut(:,:,:) = 0._r8

    taux(:,:,:)=0._r8
    tauy(:,:,:)=0._r8

    if ( ngwv > 0) then

       ! Load bottom values of Reynolds stresses
!       k=kbot

       do l=-ngwv, ngwv
          do i=1, ncol
             k = src_level(i)
             tausg(i,l,k) = sign( tau(i,l,k), c(i,l)-ubi(i,k) )
             end do
          end do

          taulu(:,:) = 0._r8
          tauru(:,:) = 0._r8
          do i=1, ncol
             do l=-ngwv, ngwv

                k = src_level(i)

                if ( c(i,l) > ubi(i,k) ) then

                   tauru(i,k) = tauru(i,k) + tausg(i,l,k)

                else if ( c(i,l) < ubi(i,k) ) then

                   taulu(i,k) = taulu(i,k) + tausg(i,l,k)

                end if

             end do
          end do

          do i=1, ncol

             k = src_level(i)

             if (xv(i) > 0._r8) then

                taue(i,k) = tauru(i,k) * xv(i)
                tauw(i,k) = taulu(i,k) * xv(i)

             else if ( xv(i) < 0._r8) then

                taue(i,k) = taulu(i,k) * xv(i)
                tauw(i,k) = tauru(i,k) * xv(i)

             end if



             if ( yv(i) > 0._r8) then

                taun(i,k) = tauru(i,k) * yv(i)
                taus(i,k) = taulu(i,k) * yv(i)

             else if ( yv(i) < 0._r8) then

                taun(i,k) = taulu(i,k) * yv(i)
                taus(i,k) = tauru(i,k) * yv(i)

             end if

          end do

       end if

!---------------------------------------------------------------------------
! Compute the stress profiles and diffusivities
!---------------------------------------------------------------------------

! Loop from bottom to top to get stress profiles      

! C.-C. Chen, move k loop into i loop
!    do k = kbot-1, ktop, -1

! Determine the absolute value of the saturation stress and the diffusivity
! for each wave.
! Define critical levels where the sign of (u-c) changes between interfaces.

       do l = -ngwv, ngwv
          do i = 1, ncol

          do k = src_level(i)-1, ktop, -1

             ubmc(i,l) = ubi(i,k) - c(i,l)
             tausat(i,l) = abs (effkwv * rhoi(i,k) * ubmc(i,l)**3 / (2._r8*ni(i,k)) )
             if (tausat(i,l) .le. taumin) tausat(i,l) = 0.0_r8

             if ( ubmc(i,l) / (ubi(i,k+1) - c(i,l) + 1e-8_r8) .le. 0.0_r8) tausat(i,l) = 0.0_r8

! Compute stress for each wave. The stress at this level is the min of 
! the saturation stress and the stress at the level below reduced by damping.
! The sign of the stress must be the same as at the level below.

       if (k <= nbot_molec) then

! Total diffusivity
!          d(:ncol) = dback
!          do i=1, ncol
             d(i) = dback + kvtt(i,k) 
!          end do

!          do l = -ngwv, ngwv
!             do i = 1, ncol
                ubmc2 = max(ubmc(i,l)**2, ubmc2mn)
                mi = ni(i,k) / (2._r8 * kwv * ubmc2) * (alpha(k) + ni(i,k)**2/ubmc2 * d(i))
                wrk = -2._r8*mi*rog*t(i,k+1)*(piln(i,k+1) - piln(i,k))
                if( wrk >= -150._r8 ) then
                   taudmp = tau(i,l,k+1) * exp( wrk )
                else
                   taudmp = 0._r8
                end if
                if( taudmp <= taumin ) then
		      taudmp = 0._r8
                end if
                tau(i,l,k) = min (taudmp, tausat(i,l))			
!             end do
!          end do
       else
!          do l= -ngwv, ngwv
!             do i=1, ncol
                tau(i,l,k) = min (tau(i,l,k+1), tausat(i,l))
!             end do
!          end do

       end if

! The orographic stress term does not change across the source region
! Note that k ge ksrcmn cannot occur without an orographic source term

       if((l.eq.0).and.(k.ge.ksrcmn).and.(k.ge.ksrc(i))) then
           tau(i,0,k) = tau(i,0,pver)
       endif

 !      if (k .ge. ksrcmn) then
 !         do i = 1, ncol
 !            if (k .ge. ksrc(i)) then
 !               tau(i,0,k) = tau(i,0,pver) 
 !            end if
 !         end do
 !      end if

! Require that the orographic term decrease linearly (with pressure) 
! within the low level stress divergence region. This supersedes the above
! requirment of constant stress within the source region.
! Note that k ge kldvmn cannot occur without an orographic source term, since
! kldvmn=pver then and k<=pver-1

       if((l.eq.0).and.(k.ge.kldvmn).and.(k.ge.kldv(i))) then
                tau(i,0,k) = min (tau(i,0,k), tau(i,0,pver)  * &
                     (1._r8 - fracldv * (pi(i,k) - pi(i,pver)) * rdpldv(i)))
       endif

!       if (k .ge. kldvmn) then
!          do i = 1, ncol
!             if (k .ge. kldv(i)) then
!                tau(i,0,k) = min (tau(i,0,k), tau(i,0,pver)  * &
!                     (1._r8 - fracldv * (pi(i,k) - pi(i,pver)) * rdpldv(i)))
!             end if
!          end do
!       end if

          enddo
          end do
       end do




! Diagnostic output of TAU
       if ( ngwv > 0) then
          do l=-ngwv, ngwv
             do i=1, ncol
             do k=src_level(i)-1,ktop,-1
                tausg(i,l,k) = sign( tau(i,l,k), c(i,l)-ubi(i,k) )

                taulu(i,k) = 0._r8
                tauru(i,k) = 0._r8
                taulv(i,k) = 0._r8
                taurv(i,k) = 0._r8

             enddo
             enddo
           enddo

          do l = -ngwv,ngwv
             do i = 1,ncol
             do k = src_level(i)-1,ktop,-1

                if ( c(i,l) < ubi(i,src_level(i)) ) then
                   taulu(i,k) = taulu(i,k) + tausg(i,l,k)
                   taulv(i,k) = taulv(i,k) + tausg(i,l,k)
                else if (  c(i,l) > ubi(i,src_level(i))  ) then
                   tauru(i,k) = tauru(i,k) + tausg(i,l,k)
                   taurv(i,k) = taurv(i,k) + tausg(i,l,k)
                end if
!             end do
!          end do
             
             enddo
             end do
          end do


          do i = 1, ncol
          do k = src_level(i)-1,ktop,-1
             if (xv(i) > 0._r8) then
                taue(i,k) = tauru(i,k) * xv(i)
                tauw(i,k) = taulu(i,k) * xv(i)
             else if ( xv(i) < 0._r8) then
                taue(i,k) = taulu(i,k) * xv(i)
                tauw(i,k) = tauru(i,k) * xv(i)
             end if

             if ( yv(i) > 0._r8) then
                taun(i,k) = taurv(i,k) * yv(i)
                taus(i,k) = taulv(i,k) * yv(i)
             else if ( yv(i) < 0._r8) then
                taun(i,k) = taulv(i,k) * yv(i)
                taus(i,k) = taurv(i,k) * yv(i)
             end if
          enddo
          end do

       end if

!    end do


!---------------------------------------------------------------------------
! Compute the tendencies from the stress divergence.
!---------------------------------------------------------------------------

!++rrg: define polar taper 

  al0   = 82.5_r8 * d2r ! 3.14159_r8/180._r8
  dlat0 =  5.0_r8 * d2r ! 3.14159_r8/180._r8

  if(ngwv > 0 .and. .not.bflag) then 	! taper CM only
    do i = 1, ncol
!     ptaper(i) = 0.25_r8 * (1.+tanh((rlat(i)+al0)/dlat0)) * (1.-tanh((rlat(i)-al0)/dlat0)) 
      ptaper(i) = cos(rlat(i))
    end do
  else					! do not taper other cases
      ptaper(:) = 1._r8
  endif

!--rrg

!----------------------------------------------------------------------------
!  do k = ktop+1, kbot			! Loop over levels from top to bottom
!----------------------------------------------------------------------------

! Accumulate the mean wind tendency over wavenumber.
       do i = 1, ncol
       do k = ktop+1,src_level(i)
          ubt(i,k) = 0.0_r8
       enddo
       end do

       do l = -ngwv, ngwv	! loop over wave
          do i = 1, ncol	! loop over column
          do k = ktop+1,src_level(i)

! Determine the wind tendency including excess stress carried down from above.
             ubtl = g * (tau(i,l,k)-tau(i,l,k-1)) * rdpm(i,k)

! Apply tendency limits to maintain numerical stability.
! 1. du/dt < |c-u|/dt  so u-c cannot change sign (u^n+1 = u^n + du/dt * dt)
! 2. du/dt < tndmax    so that ridicuously large tendencies are not permitted
             ubtl = min(ubtl, umcfac * abs(c(i,l)-ubm(i,k)) / dt)
             ubtl = min(ubtl, tndmax) 

! accumulate the mean wind tendency over wavenumber, aplying efficiency and taper:
             ubt (i,k) = ubt (i,k) + sign(ubtl, c(i,l)-ubm(i,k)) * effgw * ptaper(i)

! save tendency for each wave (for later computation of kzz), applying efficiency and taper:
             gwut(i,k,l) = sign(ubtl, c(i,l)-ubm(i,k)) * effgw * ptaper(i)

! Accumulate wind tendencies based on the phase speed of spectral components
             if (ngwv > 0 ) then

                IF (c(i,l).lt.-40._r8) THEN

                   ut1(i,k) = ut1(i,k) + sign(ubtl, c(i,l)-ubm(i,k)) * effgw * ptaper(i)

                ELSE IF ((c(i,l).ge.-40._r8).and.(c(i,l).lt.-15._r8)) THEN

                   ut2(i,k) = ut2(i,k) + sign(ubtl, c(i,l)-ubm(i,k)) * effgw * ptaper(i)

                ELSE IF ((c(i,l).ge.-15._r8).and.(c(i,l).le.15._r8)) THEN

                   ut3(i,k) = ut3(i,k) + sign(ubtl, c(i,l)-ubm(i,k)) * effgw * ptaper(i)

                ELSE IF ((c(i,l).gt.15._r8).and.(c(i,l).le.40._r8)) THEN

                   ut4(i,k) = ut4(i,k) + sign(ubtl, c(i,l)-ubm(i,k)) * effgw * ptaper(i)

                ELSE IF (c(i,l).gt.40._r8) THEN

                   ut5(i,k) = ut5(i,k) + sign(ubtl, c(i,l)-ubm(i,k)) * effgw * ptaper(i)

                END IF

             end if

! Redetermine the effective stress on the interface below from the wind 
! tendency. If the wind tendency was limited above, then the new stress
! will be smaller than the old stress and will cause stress divergence in
! the next layer down. This has the effect of smoothing large stress 
! divergences downward while conserving total stress.

             tau(i,l,k) = tau(i,l,k-1) + ubtl * dpm(i,k) / g

          enddo
          end do		! end loop over column
       end do			! end loop over wave


! rrg: project the mean wind tendency (modified by effgw) onto x,y components. 
       do i = 1, ncol
       do k = ktop+1, src_level(i)

           ut(i,k) = ubt(i,k) * xv(i) 
           vt(i,k) = ubt(i,k) * yv(i)

	if (ngwv > 0) then ! zonal tendencies for spectral bands

           ut1(i,k) = ut1(i,k) * xv(i) 
           ut2(i,k) = ut2(i,k) * xv(i) 
           ut3(i,k) = ut3(i,k) * xv(i)
           ut4(i,k) = ut4(i,k) * xv(i)
           ut5(i,k) = ut5(i,k) * xv(i)


        endif

       enddo
       end do

!---------------------------------------------------------------------------
!  end do 			! end loop over level 
!---------------------------------------------------------------------------


    if (ngwv > 0) then

       if (bflag) then

          ! Output wave tendencies binned by phase speed

          call outfld ('BUTEND1', ut1, pcols, lchnk)
          call outfld ('BUTEND2', ut2, pcols, lchnk)
          call outfld ('BUTEND3', ut3, pcols, lchnk)
          call outfld ('BUTEND4', ut4, pcols, lchnk)
          call outfld ('BUTEND5', ut5, pcols, lchnk)


       else

          call outfld ('UTEND1', ut1, pcols, lchnk)
          call outfld ('UTEND2', ut2, pcols, lchnk)
          call outfld ('UTEND3', ut3, pcols, lchnk)
          call outfld ('UTEND4', ut4, pcols, lchnk)
          call outfld ('UTEND5', ut5, pcols, lchnk)


       end if

       call gw_ediff (lchnk   ,ncol                               &
            ,ngwv     ,ktopbg   ,kbotbg  ,gwut      ,ubm          &
            ,nm       ,t       ,ksrc     ,egwdffi  ,c) 
       
       do m = 1, pcnst

          call gw_const_tend  (lchnk   ,ncol                      &
               ,ngwv         ,kbotbg      ,ktopbg     ,q(1,1,m)   &
               ,pi          ,pm           ,t          ,dt         &
               ,rdpm        ,qtgw(1,1,m)  ,egwdffi     )

       enddo

       call gw_temp_tend  (lchnk           ,ncol                 &
            ,ngwv     ,kbotbg     ,ktopbg  ,dse                  &
            ,pi       ,pm         ,t       ,gwut                 &
            ,dt       ,rdpm       ,egwdffi ,ttgw     ,c)

    else 

       qtgw(:,:,:) = 0._r8
       ttgw(:,:) = 0._r8
       
    end if

!-----------------------------------------------------------------------
! Project the c=0 stress (scaled) in the direction of the source wind for recording
! on the output file.
!-----------------------------------------------------------------------
    do i = 1, ncol
       tau0x(i) = tau(i,0,kbot) * xv(i) * effgw
       tau0y(i) = tau(i,0,kbot) * yv(i) * effgw
    end do


! Output Tau with height - output taux and tauy separately

    if (ngwv > 0) then

       do i=1,ncol
          do l=-ngwv,ngwv
             cx = c(i,l)*xv(i)
             cy = c(i,l)*yv(i)
             ix = NINT (cx/dc)
             iy = NINT (cy/dc)
             if (ix < -ngwv) then
                ix = -ngwv
             else if (ix > ngwv) then
                ix = ngwv
             endif
             if (iy < -ngwv) then
                iy = -ngwv
             else if (iy > ngwv) then
                iy = ngwv
             endif
             taux(i,ix,1:pver) = taux(i,ix,1:pver) + ABS(tau(i,l,1:pver)*xv(i))
             tauy(i,iy,1:pver) = tauy(i,iy,1:pver) + ABS(tau(i,l,1:pver)*yv(i))
          enddo
       enddo

       do i=1,ncol
          do k=1,pver
             usi = NINT (u(i,k)/dc)
             if (usi < -ngwv) then
                usi = -ngwv
             else if (usi > ngwv) then
                usi = ngwv
             endif

             vsi = NINT (v(i,k)/dc)
             if (vsi < -ngwv) then
                vsi = -ngwv
             else if (vsi > ngwv) then
                vsi = ngwv
             end if

             do l=-ngwv,usi-1
                WMF(i,k) = WMF(i,k) + (taux(i,l,k)*effgw)
             enddo

             do l=usi,ngwv
                EMF(i,k) = EMF(i,k) + (taux(i,l,k)*effgw)
             enddo

             do l=-ngwv,vsi-1
                SMF(i,k) = SMF(i,k) + (tauy(i,l,k)*effgw)
             enddo

             do l=vsi,ngwv
                NMF(i,k) = NMF(i,k) + (tauy(i,l,k)*effgw)
             enddo
          enddo
       enddo
       
       if (bflag) then
          call outfld ('BEMF',EMF(1:pcols,1:pver),pcols,lchnk)
          call outfld ('BWMF',WMF(1:pcols,1:pver),pcols,lchnk)
          call outfld ('BNMF',NMF(1:pcols,1:pver),pcols,lchnk)
          call outfld ('BSMF',SMF(1:pcols,1:pver),pcols,lchnk)
       else
          call outfld ('EMF',EMF(1:pcols,1:pver),pcols,lchnk)
          call outfld ('WMF',WMF(1:pcols,1:pver),pcols,lchnk)
          call outfld ('NMF',NMF(1:pcols,1:pver),pcols,lchnk)
          call outfld ('SMF',SMF(1:pcols,1:pver),pcols,lchnk)
       endif

! Positive phase velocities

       if (bflag) then
          do m=0,ngwv

             dummyx(:pcols,:pver)=0._r8
             dummyy(:pcols,:pver)=0._r8

             do i=1,ncol
                do l=-ngwv,ngwv
                   if(abs(cref(m)-c(i,l)) <= dc/2._r8) then

                      dummyx(i,1:pver)=taux(i,l,1:pver)
                      dummyy(i,1:pver)=tauy(i,l,1:pver)

                   end if

                enddo
             enddo
             write (fnum,fmt='(i3)') 100+m

             dumc1x='BTAUXSp'//fnum(2:3)
             dumc1y='BTAUYSp'//fnum(2:3)

             call outfld (dumc1x,dummyx(1:pcols,1:pver),pcols,lchnk)
             call outfld (dumc1y,dummyy(1:pcols,1:pver),pcols,lchnk)

          enddo
       else
          do m=0,ngwv

             dummyx(:pcols,:pver)=0._r8
             dummyy(:pcols,:pver)=0._r8

             do i=1,ncol
                do l=-ngwv,ngwv
                   if(abs(cref2(m)-c(i,l)) <= dc/2._r8) then

                      dummyx(i,1:pver)=taux(i,l,1:pver)
                      dummyy(i,1:pver)=tauy(i,l,1:pver)

                   end if

                enddo
             enddo
             write (fnum,fmt='(i3)') 100+m

             dumc1x='TAUXSp'//fnum(2:3)
             dumc1y='TAUYSp'//fnum(2:3)

             call outfld (dumc1x,dummyx(1:pcols,1:pver),pcols,lchnk)
             call outfld (dumc1y,dummyy(1:pcols,1:pver),pcols,lchnk)

          enddo
       endif

! Negative phase velocities

       if (bflag) then
          do m=-1,-ngwv, -1
             dummyx(:pcols,:pver)=0._r8
             dummyy(:pcols,:pver)=0._r8

             do i=1,ncol
                do l=-ngwv,ngwv
                   if(abs(cref(m)-c(i,l)) <= dc/2._r8) then

                      dummyx(i,1:pver)=taux(i,l,1:pver)
                      dummyy(i,1:pver)=tauy(i,l,1:pver)


                   end if

                enddo
             enddo
             write (fnum,fmt='(i3)') 100-m

             dumc1x='BTAUXSn'//fnum(2:3)
             dumc1y='BTAUYSn'//fnum(2:3)

             call outfld (dumc1x,dummyx(1:pcols,1:pver),pcols,lchnk)
             call outfld (dumc1y,dummyy(1:pcols,1:pver),pcols,lchnk)
          enddo
       else
          do m=-1,-ngwv, -1
             dummyx(:pcols,:pver)=0._r8
             dummyy(:pcols,:pver)=0._r8

             do i=1,ncol
                do l=-ngwv,ngwv
                   if(abs(cref2(m)-c(i,l)) <= dc/2._r8) then

                      dummyx(i,1:pver)=taux(i,l,1:pver)
                      dummyy(i,1:pver)=tauy(i,l,1:pver)

                   end if

                enddo
             enddo
             write (fnum,fmt='(i3)') 100-m

             dumc1x='TAUXSn'//fnum(2:3)
             dumc1y='TAUYSn'//fnum(2:3)

             call outfld (dumc1x,dummyx(1:pcols,1:pver),pcols,lchnk)
             call outfld (dumc1y,dummyy(1:pcols,1:pver),pcols,lchnk)
          enddo
       endif
    end if  !ngwv > 0


    return
  end subroutine gw_drag_prof

!=======================================================================================
subroutine gw_const_tend (lchnk   ,ncol                                            &
                           ,ngwv    ,kbot      ,ktop    ,q       ,pi      ,pm      &
                           ,t       ,dt        ,rdpm    ,dq      ,egwdffi )

!
!  calculates constituent tendencies due to gw breaking
!
!  method:
!  a constituent flux on interfaces is given by
!
!                                   
!              rho * (w'q') = rho * Deff qz
!
!
! where (all evaluated on interfaces):
!
!        rho   = density
!        qz    = constituent vertical gradient
!        Deff  = effective diffusivity
!
!
! an effective diffusivity is calculated by adding up the diffusivities from all waves
! (see gw_ediff)
! tendency is calculated by invoking  lu decomposition and solving as for a regular
! diffusion equation.
!
! author : sassi - jan 2001
!---------------------------------------------------------------------------------

 implicit none

!---------------------------input arguments---------------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: kbot                   ! index of bottom (source) midpoint
    integer, intent(in) :: ktop                   ! index of top midpoint of gwd region
    integer, intent(in) :: ngwv                   ! number of gravity waves to use

    real(r8), intent(in) :: q(pcols,pver)                 ! constituent
    real(r8), intent(in) :: t(pcols,pver)                 ! temperature at midpoints
    real(r8), intent(in) :: dt                            ! time step
    real(r8), intent(in) :: rdpm(pcols,pver)              ! 1./pdel  (thickness bet interfaces)
    real(r8), intent(in) :: pi(pcols,pverp)               ! interface pressure
    real(r8), intent(in) :: pm(pcols,pver)                ! midpoint  pressure
    real(r8), intent(in) :: egwdffi(pcols,pverp)          ! effective gw diffusivity at interfaces

!--------------------------output arguments-------------------------------------
 
    real(r8), intent(out) :: dq(pcols,pver)              ! constituent tendencies at midpoints

!--------------------------local worspace---------------------------------------

    integer i,k,l
    real(r8) tmpi2(pcols,pverp)                  ! factor
    real(r8) ca(pcols,pver)                       ! -upper diagonal
    real(r8) cc(pcols,pver)                       ! -lower diagonal
    real(r8) term(pcols,pver)                     ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
    real(r8) ze(pcols,pver)                       ! term in tri-diag. matrix system
    real(r8) qconst(pcols,pver)                   ! temporary storage for constituent
    real(r8) rho(pcols,pverp)
    real(r8) zero(pcols)

!-------------------------------------------------------------------------------
    zero       = 0.0_r8
    dq(:,:)    = 0.0_r8

! Compute factor for diffusion matrices at interfaces (except top), 
! Compute rho at interfaces p/RT,  Ti = (Tm_k + Tm_k-1)/2,  interface k-
    do k = 2, pver
       do i = 1, ncol
          rho(i,k) = pi(i,k) * 2._r8 / (r * (t(i,k) + t(i,k-1)))
       end do
    end do
    do i = 1, ncol
       rho(i,pverp) = pi (i,pverp) / (r * t(i,pver))
    end do
    do i = 1, ncol
       rho(i,1) = pi (i,1) / (r * t(i,1))
    end do

! Calculate dt * (g*rho)^2/dp at interior interfaces,  interface k-
    tmpi2(:,:) = 0.0_r8
    do k = ktop+2 , kbot+1
       do i = 1, ncol
          tmpi2(i,k) = dt * (g * rho(i,k))**2 / (pm(i,k) - pm(i,k-1))
       end do
    end do

   qconst = 0.0_r8
   do k = 1, pver
      do i = 1, ncol
         qconst(i,k) = q(i,k)
      enddo
   enddo

! Decompose and solve the diffusion matrix
     call vd_lu_decomp(                                                &
          lchnk      ,ncol       ,                                     &
          zero       ,egwdffi    ,tmpi2      ,rdpm       ,dt         ,zero , &
          ca         ,cc         ,term       ,ze         ,ktop+1     , &
          kbot+1 )
       call vd_lu_solve(                                                 &
            lchnk      ,ncol       ,                                     &
            qconst(1,1),ca         ,ze         ,term       ,             &
            ktop+1      ,kbot+1    ,zero   )

! Evaluate tendency to be reported back
      do k = ktop+1, kbot
         do i = 1, ncol
            dq(i,k) = (qconst(i,k) - q(i,k)) / dt
         enddo
     enddo
 
  return
end subroutine gw_const_tend


!=====================================================================================

subroutine gw_ediff (lchnk   ,ncol                             &
                    ,ngwv    ,ktop    ,kbot   ,gwut    ,ubm    &
                    ,nm      ,t                       &
                    ,ksrc    ,egwdffi ,c        )

! 
! Calculate effective diffusivity associated with GW forcing.
!
! Author: F. Sassi, Jan 31, 2001
!

implicit none

!-------------------------------Input Arguments------------------------------------------

    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: kbot                   ! index of bottom (source) midpoint
    integer, intent(in) :: ktop                   ! index of top midpoint of gwd region
    integer, intent(in) :: ngwv                   ! number of gravity waves to use
    integer, intent(in) :: ksrc(pcols)            ! index of top interface of source region

    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperature
    real(r8), intent(in) :: gwut(pcols,pver,-ngwv:ngwv)   ! GW zonal wind tendencies at midpoint
    real(r8), intent(in) :: ubm(pcols,pver)               ! projection of wind at midpoints
    real(r8), intent(in) :: nm(pcols,pver)                ! Brunt-Vaisalla frequency    
    real(r8), intent(in) :: c(pcols,-ngwv:ngwv)           ! wave phase speeds for each column

!-----------------------------Output Argument-------------------------------------------
    real(r8), intent(out) :: egwdffi(pcols,0:pver)  ! effective gw diffusivity at interfaces

!-----------------------------Local workspace-------------------------------------------

    real(r8) egwdffm(pcols,pver)                   ! effective gw diffusivity at midpoints
    integer i,k,l,m
!   real(r8), parameter :: prndl=0.50              ! Inverse Prandtl no.       !fs
    real(r8), parameter :: prndl=0.25_r8           ! Inverse Prandtl no.       !rrg 1/25/2005
    real(r8) :: coef
    real(r8) :: nms2(pver)                         ! zonal average nm^2
    real(r8) :: rho(pcols,pverp)                   ! density at interfaces
    real(r8), parameter :: dscale=7000._r8         ! density scale height

    real(r8) :: dummy(pcols,pver)
    real(r8) :: ccn,cmx,cmn
    character(len=8) :: fnum,fname,dumc1           ! dummy characters

!---------------------------------------------------------------------------------------

! zero
  
  egwdffi(:,:) = 0._r8
  egwdffm(:,:) = 0._r8

  if (ngwv == 0) return
  

! Calculate effective diffusivity at midpoints
    do k = ktop+1, kbot
       do i = 1, ncol
          do l = -ngwv, ngwv

             egwdffm(i,k) = egwdffm(i,k) + &
                    prndl * 0.5_r8 * gwut(i,k,l) * (c(i,l)-ubm(i,k)) / nm(i,k)**2  !rrg/bab

          enddo                 ! another wave
      enddo                     ! another long
    enddo                       ! another level


! Interpolate effective diffusivity to interfaces
   do i = 1, ncol
! Assume zero  at top and bottom interfaces
      do k=ktop+1, kbot-1
         egwdffi(i,k) = .5_r8*(egwdffm(i,k) + egwdffm(i,k+1))
      enddo
   enddo

! Limit diffusivity to some reasonable value
    do i=1, ncol
       do k=0, pver
          egwdffi(i,k) = min( 150._r8,egwdffi(i,k) )
       enddo
    enddo


! Do not calculate diffusivities below source level
   do k= ktop,kbot
      do i = 1, ncol
         if (k .ge. ksrc(i)) egwdffi(i,k) = 0.0_r8
      enddo
   enddo

   return
end subroutine gw_ediff

!========================================================================================

  subroutine gw_temp_tend (lchnk   ,ncol                                              &
                          ,ngwv    ,kbot   ,ktop    ,dse     ,pi       ,pm            &
                          ,t       ,gwut   ,dt      ,rdpm                    &
                          ,egwdffi ,ttgw   ,c)

!
!  Calculates temperature tendency due to GW breaking
!
!  Method:
!  GW mixing is calculated using the effective diffusivity obtained from gw_const_tend: 
!  it is used to diffuse dry static energy.
!
! Author : Sassi - Jan 2001
!---------------------------------------------------------------------------------

 implicit none

!---------------------------Input Arguments---------------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: kbot                   ! index of bottom (source) interface
    integer, intent(in) :: ktop                   ! index of top interface
    integer, intent(in) :: ngwv                   ! number of gravity waves to use

    real(r8), intent(in) :: gwut(pcols,pver,-ngwv:ngwv)   ! GW zonal wind tendencies at midpoint
    real(r8), intent(in) :: t(pcols,pver)                 ! temperature at midpoints
    real(r8), intent(in) :: dt                            ! time step
    real(r8), intent(in) :: rdpm(pcols,pver)              ! 1./pdel  (thickness bet interfaces)
    real(r8), intent(in) :: pi(pcols,pverp)               ! interface pressure
    real(r8), intent(in) :: pm(pcols,pver)                ! midpoint pressure
    real(r8), intent(in) :: egwdffi(pcols,pverp)          ! effective gw diffusivity at interfaces
    real(r8), intent(in) :: dse(pcols,pver)               ! midpoint dry static energy
    real(r8), intent(in) :: c(pcols,-ngwv:ngwv)                       ! wave phase speeds for each column

!--------------------------Output Arguments-------------------------------------

    real(r8), intent(out) :: ttgw(pcols,pver) 
 
!--------------------------Local Worspace---------------------------------------

    integer i,k,m,l
    real(r8) :: delta                        ! dissipation rate
    real(r8) :: tmpi2(pcols,pverp)           ! pi(k)/(.5*(tm1(k)+tm1(k-1))
    real(r8) :: ca(pcols,pver)               ! -upper diagonal
    real(r8) :: cc(pcols,pver)               ! -lower diagonal
    real(r8) :: term(pcols,pver)             ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
    real(r8) :: ze(pcols,pver)               ! term in tri-diag. matrix system
    real(r8) :: temp(pcols,pver)             ! temporary storage for constituents
    real(r8) :: rho(pcols,pverp)
    real(r8) :: zero(pcols)
    real(r8) :: dttdf(pcols,pver)            ! dse tendencies at midpoints: diffusion term
    real(r8) :: dttke(pcols,pver)            ! dse tendencies at midpoints: kin. energy conversion term

!-------------------------------------------------------------------------------
    zero(:)    = 0.0_r8
    dttdf(:,:) = 0.0_r8
    dttke(:,:) = 0.0_r8

    if (ngwv .eq. 0) return

!
! Compute factor for diffusion matrices at interior interfaces 
! Note that [ktop,kbot] are model interfaces (beginning at zero), 
! whereas in vd_lu_decomp they are expected as midpoints

! Compute rho at interfaces p/RT,  Ti = (Tm_k + Tm_k-1)/2,  interface k-
    do k = 2, pver
       do i = 1, ncol
          rho(i,k) = pi(i,k) * 2._r8 / (r * (t(i,k) + t(i,k-1)))
       end do
    end do
    do i = 1, ncol
       rho(i,pverp) = pi (i,pverp) / (r * t(i,pver))
    end do
    do i = 1, ncol
       rho(i,1) = pi (i,1) / (r * t(i,1))
    end do

! Calculate dt * (g*rho)^2/dp at interior interfaces,  interface k-
    tmpi2(:,:) = 0.0_r8
    do k = ktop+2 , kbot+1
       do i = 1, ncol
          tmpi2(i,k) = dt * (g * rho(i,k))**2 / (pm(i,k) - pm(i,k-1))
       end do
    end do

    do k = 1, pver
       do i = 1, ncol
          temp(i,k) =  dse(i,k)
       enddo
    enddo

! Decompose and solve the diffusion matrix
    call vd_lu_decomp(                                                &
         lchnk      ,ncol       ,                                     &
         zero       ,egwdffi    ,tmpi2       ,rdpm       ,dt         ,zero ,&
         ca         ,cc         ,term       ,ze         ,ktop+1      ,&
         kbot+1 )
    call vd_lu_solve(                                                 &
         lchnk      ,ncol       ,                                     &
         temp       ,ca         ,ze         ,term       ,             &
         ktop+1      ,kbot+1    ,zero   )

! Evaluated first part of tendency : diffusive term
   do k = 1, pver
      do i =1, ncol
         dttdf(i,k) = (temp(i,k) - dse(i,k)) / dt
      enddo
  enddo

! Evaluate second term: Conversion of kinetic energy into thermal
  do k = ktop+1, kbot
     do i = 1, ncol
        do l = -ngwv, ngwv
           dttke(i,k) = dttke(i,k) +  c(i,l) * gwut(i,k,l) 
        enddo
     enddo
  enddo

  do k=1, pver
     do i=1, ncol

        ttgw(i,k) = dttke(i,k) + dttdf(i,k)

     end do
  end do

  call outfld ('TTGWSDF', dttdf  / cpair, pcols, lchnk)
  call outfld ('TTGWSKE', dttke / cpair, pcols, lchnk)

  return
end subroutine gw_temp_tend


!==============================================================================


!+Mods for Beres04
subroutine gw_bgnd_beres (lchnk, ncol, &
     u, v, t, pm, pi, dpm, rdpm, piln,     &
     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv, &
     ngwv, zm, src_level, pbuf)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!
! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! specifying the gravity wave spectrum above convection based on latent
! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! pp. 324-337.
!
!-----------------------------------------------------------------------
  use phys_grid,    only : get_rlat_all_p, get_rlon_all_p
  use phys_grid,    only : get_lat_all_p, get_lon_all_p
  
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field

!------------------------------Arguments--------------------------------
  integer, intent(in) :: lchnk                  ! chunk identifier
  integer, intent(in) :: ncol                   ! number of atmospheric columns

  integer,  intent(in) :: ngwv                  ! number of gravity waves to use
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
  real(r8), intent(in) ::  zm(pcols,pver)        ! midpoint heights
! C.-C. Chen
  integer, intent(out) :: src_level(pcols)      ! index of source (top of heating)
  type(physics_buffer_desc), pointer :: pbuf(:)

!---------------------------Local storage-------------------------------
  integer :: i,k,l                              ! loop indexes

  real(r8) :: rlon(pcols)                       ! longitude in radians for columns
  real(r8) :: rlat(pcols)                       ! latitude in radians for columns
  real(r8) :: tauback(pcols)                    ! background stress at c=0
  real(r8) :: u700(pcols)                       ! u at 700mb
  real(r8) :: v700(pcols)                       ! v at 700mb
  real(r8) :: al0                               ! Used in lat dependence of GW spec.
  real(r8) :: dlat0                             ! Used in lat dependence of GW spec.
  real(r8) :: flat_gw                           ! The actual lat dependence of GW spec.
  real(r8) :: pi_g                              ! 3.14........

!-------------------------- More local storage required for the Beres scheme -----------
  real(r8) :: netdt(pcols,pver)                 ! heating rate
  integer  :: lons(pcols)
  integer  :: lats(pcols)
  real(r8) :: curq(pver)           ! heating rate in given column
  real(r8) :: q0                                ! maximum heating rate
  real(r8) :: curu(pver)                        ! mean wind in a given column
                                                  ! projected in direction of 700 mb wind
  real(r8) :: depth                             ! heating depth
  real(r8) :: hdepth(pcols)                     ! heating depth array
  real(r8) :: maxq0(pcols)                      ! maximum heating rate array
!  real(r8) :: BSEMF(pcols)                      ! source total Eastward MF
!  real(r8) :: BSWMF(pcols)                      ! source total Westward MF

  real(r8) :: htop                              ! heating top
  integer  :: mini,maxi                         ! min and max heating range index
  integer  :: Umini,Umaxi                       ! Critical level filteringrange
  integer  :: maxsi                               ! maximum heating depth height index
  real(r8) :: uh                                ! mean wind in heating region
  real(r8) :: Umin, Umax                        ! Min and Max values of U in column
  real(r8) :: tau0(-PGWV:PGWV)                  ! source level tau for one column
  integer  CS                                   ! speed of convective cells relative to storm
  integer  shift                                ! how much to shift spectra to ground relative


  real(r8) :: ustorm                             ! storm speed
  real(r8), parameter :: CF = 20._r8                ! Heating rate conversion factor

  real(r8), parameter :: AL = 1.0e5_r8               ! Averaging length

  real(r8) :: dummy(pcols)
  real(r8) :: ccn
  character(len=80) fnum,fmt ,dumc1          ! dummies
  integer :: m
  real(r8), pointer, dimension(:,:) :: zmdt

  call get_lon_all_p(lchnk, ncol, lons)
  call get_lat_all_p(lchnk, ncol, lats)

  call get_rlat_all_p(lchnk, ncol, rlat)
  call get_rlon_all_p(lchnk, ncol, rlon)

!       do i = 1, ncol

!          write(iulog,*) 'ncol: ',ncol
!          write(iulog,*) 'Current, lat & lon: ',lats(i), lons (i)
!          write(iulog,*) 'Current, lat & lon: ',rlat(i)*180./3.14, rlon (i)*180./3.14

!       end do

  pi_g=4._r8*atan(1._r8)
  al0=40._r8*pi_g/180._r8
  dlat0=10._r8*pi_g/180._r8

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------

!  BSEMF(:) = 0.0
!  BSWMF(:) = 0.0
  tau(:,:,:)=0.0_r8
  netdt(:,:) = 0.0_r8
  hdepth(:) = 0.0_r8
  maxq0(:) = 0.0_r8
  dummy(:) = 0.0_r8
  tau0(:)  = 0.0_r8

  !---------------------------------------------------------------------------
  ! Determine 700 mb layer wind and unit vectors, then project winds.
  !---------------------------------------------------------------------------

  ! Just use the 700 mb interface values for the source
  ! wind speed and direction (unit vector).

  do i = 1, ncol
     u700(i) = u(i,k700)
     v700(i) = v(i,k700)
     ubi(i,k700) = sqrt (u700(i)**2 + v700(i)**2)
     if (ubi(i,k700) /= 0._r8) then
        xv(i) = u700(i) / ubi(i,k700)
        yv(i) = v700(i) / ubi(i,k700)
     else ! Magnitude 0 vector.
        xv(i) = 0._r8
        yv(i) = 0._r8
     end if
  end do

  ! Project the local wind at midpoints onto the direction of 700mb wind.
  do k = 1, pver
     do i = 1, ncol
        ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
     end do
  end do


  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ! Use the bottom level wind at surface

  do i = 1, ncol
     ubi(i,0) = ubm(i,1)
  end do
  do k = 1, pver-1
     do i = 1, ncol
        ubi(i,k) = 0.5_r8 * (ubm(i,k) + ubm(i,k+1))
     end do
  end do
  do i = 1, ncol
     ubi(i,pver) = ubm(i,pver)
  end do

  !---------------------------------------------------------------------------
  ! Find the index of z = 20 km (maximum possible heating range)
  ! using just the first column value of zm
  !---------------------------------------------------------------------------
  k=pver
  do
     if (zm(1,k).ge. 20000._r8) EXIT
     k=k-1
  end do
  maxsi = k


  !---------------------------------------------------------------------------
  ! Set up heating
  !---------------------------------------------------------------------------
  call pbuf_get_field(pbuf, idx_zmdt, zmdt)

  netdt(:pcols,:pver) = zmdt(:pcols,:pver)

  call outfld ('NETDT', netdt, pcols, lchnk)


  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !!-----------------------------------------------------------------------
  !! Start LOOP over ALL COLUMNS
  !!-----------------------------------------------------------------------
  do i=1,ncol

     curq(1:pver) = netdt(i,1:pver)
     curu(1:pver) = ubm(i,1:pver)

     !-----------------------------------------------------------------------
     ! calculate HEATING DEPTH
     !
     ! heating depth is defined as the first  height range from the bottom
     ! in which heating rate is continuously positive
     !-----------------------------------------------------------------------
     k=pver


     ! Find first level of positive heating rate
     do
        if (k.le. maxsi) EXIT
        if (curq(k) .gt. 0.0_r8) EXIT
        k=k-1
     end do

     mini=k

     ! Find first level at which heating rate stops being positive
     do
        if (k.le. maxsi) EXIT
        if (curq(k) .le. 0.0_r8) EXIT
        k=k-1
     end do

     maxi=k
     depth = (zm(i,maxi)-zm(i,mini))/1000._r8
     htop  = zm(i,maxi)
     ! Depth is now in km

     hdepth(i) = depth
     ! Maximum heating rate
     q0 = MAXVAL(curq(maxi:mini))

     !output max heating rate in K/day
     maxq0(i) = q0*24._r8*3600._r8

     !-----------------------------------------------------------------------
     !Look up spectrum only if depth ge 2.5 km, else set tau0 = 0
     !-----------------------------------------------------------------------

     if ((depth .ge. 2.5_r8).AND.(ABS(rlat(i)).lt.(pi_g/2._r8))) then

        !-----------------------------------------------------------------------
        ! calculate MEAN WIND IN HEATING REGION relative to STORM CELLS
        !-----------------------------------------------------------------------

        ustorm = ubm(i,k700)               ! storm speed


        ! Calculate cell speed (CS) for ustrom > 10 m/s

        if (ABS(ustorm) .gt. 10._r8) then
           CS = SIGN(ABS(ustorm)-10._r8,ustorm)
        else
           CS = 0
        end if

        uh = 0.0_r8

        do k=maxi,mini
           uh = uh + ubm(i,k)/(mini-maxi+1)
        end do
        uh = uh - CS

        !-----------------------------------------------------------------------
        ! LOOK UP the spectrum using depth and uh
        ! Round both to the NEAREST WHOLE NUMBER
        ! SPECTRUM IS RELATIVE TO HEATING CELLS
        !-----------------------------------------------------------------------


        if(NINT(depth).gt.MAXH) then
           depth=MAXH*1.0_r8
           write(iulog,*) 'DEPTH exceeding table value'
        end if

        if(NINT(uh) > MAXUH .or. NINT(uh) < -MAXUH) then
           uh=MAXUH*1.0_r8
           write(iulog,*) 'UH exceeding table value'
        end if

            tau0(:)=mfcc(NINT(depth),NINT(uh),:)

        ! NOW SHIFT SPECTRUM SO IT's RELATIVE TO THE GROUND

        shift = -NINT(CS/dc)

        tau0=CSHIFT(tau0,shift)

        ! Adjust Magnitude

        q0 = q0 * CF                                  ! multiply by conversion factor

        tau0(:)=tau0(:)*q0*q0/AL

        ! Adjust Magnitude for latitudinal dependence


!!$        flat_gw= 0.5*(1.+tanh((rlat(i)-al0)/dlat0)) + 0.5*(1.+tanh(-(rlat(i)+al0)/dlat0))
!!$        tau0(:)=(1-flat_gw)*tau0(:)


        ! Adjust for critical level filtering

!        Umin = Minval(ubm(i,bkbotbg:mini))            ! bkbotbg should be maxi
!        Umax = Maxval(ubm(i,bkbotbg:mini))            ! if launching waves from
                                                      ! top of heating

        Umin = Minval(ubm(i,maxi:mini))         
        Umax = Maxval(ubm(i,maxi:mini))

        Umini = max(NINT(Umin/dc),-PGWV)
        Umaxi = min(NINT(Umax/dc),PGWV)

        if (Umaxi.gt.Umini) then
           do l=Umini,Umaxi
              tau0(l) = 0.0_r8
           end do
        end if


     else

        tau0(:)=0.0_r8

     end if                                      !END IF depth ge 2.5


     ! Set GW source to the top of the heating
     !   write(iulog,*) 'maxi: ',maxi

!     ksrc(i) = bkbotbg
!     kldv(i) = bkbotbg
! C.-C. Chen
     ksrc(i) = maxi
     kldv(i) = maxi
     src_level(i) = maxi

     !launch waves at the same level for now

!     tau(i,-ngwv:ngwv,bkbotbg) = tau0(-ngwv:ngwv)
     tau(i,-ngwv:ngwv,maxi) = tau0(-ngwv:ngwv)

! Diagnostic variables
!     do l=1,ngwv
!        if (xv(i) .gt.0) then
!           BSEMF(i) = BSEMF(i)+ABS(tau0(l)*xv(i))
!           BSWMF(i) = BSWMF(i)+ABS(tau0(-l)*xv(i))
!        else
!           BSWMF(i) = BSWMF(i)+ABS(tau0(l)*xv(i))
!           BSEMF(i) = BSEMF(i)+ABS(tau0(-l)*xv(i))
!        end if

!     end do



  enddo

  !!-----------------------------------------------------------------------
  !! END LOOP over ALL COLUMNS
  !!-----------------------------------------------------------------------

  call outfld ('HDEPTH', hdepth, pcols, lchnk)
  call outfld ('MAXQ0', maxq0, pcols, lchnk)
!    call outfld ('BSEMF', BSEMF, pcols, lchnk)
!    call outfld ('BSWMF', BSWMF, pcols, lchnk)


  ! Determine the min value of kldv and ksrc for limiting later loops
  ! and the pressure at the top interface of the low level stress divergence
  ! region.

  ksrcmn = pver
  kldvmn = pver

  return
end  subroutine gw_bgnd_beres
!-Mods for Beres04


!==============================================================================
  subroutine vd_lu_decomp(                                          &
       lchnk      ,ncol       ,                                     &
       ksrf       ,kv         ,tmpi       ,rpdel      ,ztodt      , cc_top, &
       ca         ,cc         ,dnom       ,ze         ,ntop       , &
       nbot       )
!-----------------------------------------------------------------------
! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
! tridiagonal diffusion matrix. 
! The diagonal elements (1+ca(k)+cc(k)) are not required by the solver.
! Also determine ze factor and denominator for ze and zf (see solver).
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns
    integer, intent(in) :: ntop                    ! top level to operate on
    integer, intent(in) :: nbot                    ! bottom level to operate on

    real(r8), intent(in) :: ksrf(pcols)            ! surface "drag" coefficient
    real(r8), intent(in) :: kv(pcols,pverp)        ! vertical diffusion coefficients
    real(r8), intent(in) :: tmpi(pcols,pverp)      ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel  (thickness bet interfaces)
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    real(r8), intent(in) :: cc_top(pcols)          ! lower diagonal on top interface (for fixed ubc only)

    real(r8), intent(out) :: ca(pcols,pver)        ! -upper diagonal
    real(r8), intent(out) :: cc(pcols,pver)        ! -lower diagonal
    real(r8), intent(out) :: dnom(pcols,pver)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
                                                   ! 1./(b(k) - c(k)*e(k-1))
    real(r8), intent(out) :: ze(pcols,pver)        ! term in tri-diag. matrix system

!---------------------------Local workspace-----------------------------
    integer :: i                                   ! longitude index
    integer :: k                                   ! vertical index

!-----------------------------------------------------------------------
! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
! a combination of ca and cc; they are not required by the solver.

    do k = nbot-1, ntop, -1
       do i = 1, ncol
          ca(i,k  ) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k  )
          cc(i,k+1) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k+1)
       end do
    end do

! The bottom element of the upper diagonal (ca) is zero (element not used).
! The subdiagonal (cc) is not needed in the solver.

    do i=1,ncol
       ca(i,nbot) = 0._r8
    end do

! Calculate e(k).  This term is 
! required in solution of tridiagonal matrix defined by implicit diffusion eqn.

    do i = 1, ncol
       dnom(i,nbot) = 1._r8/(1._r8 + cc(i,nbot) + ksrf(i)*ztodt*gravit*rpdel(i,nbot))
       ze(i,nbot)   = cc(i,nbot)*dnom(i,nbot)
    end do
    do k = nbot-1, ntop+1, -1
       do i=1,ncol
          dnom(i,k) = 1._r8/ (1._r8 + ca(i,k) + cc(i,k) - ca(i,k)*ze(i,k+1))
          ze(i,k) = cc(i,k)*dnom(i,k)
       end do
    end do
    do i=1,ncol
       dnom(i,ntop) = 1._r8/ (1._r8 + ca(i,ntop) + cc_top(i) - ca(i,ntop)*ze(i,ntop+1))
    end do

    return
  end subroutine vd_lu_decomp

!==============================================================================
  subroutine vd_lu_solve(                                           &
       lchnk      ,ncol       ,                                     &
       q          ,ca         ,ze         ,dnom       ,             &
       ntop       ,nbot       ,cd_top )
!-----------------------------------------------------------------------
! Solve the implicit vertical diffusion equation with zero flux boundary conditions.
! Actual surface fluxes are explicit (rather than implicit) and are applied separately. 
! Procedure for solution of the implicit equation follows Richtmyer and 
! Morton (1967,pp 198-200).

! The equation solved is
!
! -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k), 
!
! where d(k) is the input profile and q(k) is the output profile
!
! The solution has the form
!
! q(k)  = ze(k)*q(k-1) + zf(k)
!
! ze(k) = cc(k) * dnom(k)
!
! zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)
!
! dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)] =  1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]

! Note that the same routine is used for temperature, momentum and tracers,
! and that input variables are replaced.
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns
    integer, intent(in) :: ntop                    ! top level to operate on
    integer, intent(in) :: nbot                    ! bottom level to operate on

    real(r8), intent(in) :: ca(pcols,pver)         ! -upper diag coeff.of tri-diag matrix
    real(r8), intent(in) :: ze(pcols,pver)         ! term in tri-diag solution
    real(r8), intent(in) :: dnom(pcols,pver)       ! 1./(1. + ca(k) + cc(k) - ca(k)*ze(k+1))
    real(r8), intent(in) :: cd_top(pcols)          ! cc_top * ubc value

    real(r8), intent(inout) :: q(pcols,pver)       ! constituent field

!---------------------------Local workspace-----------------------------
    real(r8) :: zf(pcols,pver)                     ! term in tri-diag solution

    integer i,k                                    ! longitude, vertical indices

!-----------------------------------------------------------------------

! Calculate zf(k).  Terms zf(k) and ze(k) are required in solution of 
! tridiagonal matrix defined by implicit diffusion eqn.
! Note that only levels ntop through nbot need be solved for.

    do i = 1, ncol
       zf(i,nbot) = q(i,nbot)*dnom(i,nbot)
    end do
    do k = nbot-1, ntop+1, -1
       do i=1,ncol
          zf(i,k) = (q(i,k) + ca(i,k)*zf(i,k+1))*dnom(i,k)
       end do
    end do

! Include boundary condition on top element
    k = ntop
    do i=1, ncol
       zf(i,k) = (q(i,k) + cd_top(i) + ca(i,k)*zf(i,k+1))*dnom(i,k)
    end do

! Perform back substitution

    do i=1,ncol
       q(i,ntop) = zf(i,ntop)
    end do
    do k=ntop+1, nbot, +1
       do i=1,ncol
          q(i,k) = zf(i,k) + ze(i,k)*q(i,k-1)
       end do
    end do

    return
  end subroutine vd_lu_solve

!===============================================================================
end module gw_drag
