
module gw_drag

!--------------------------------------------------------------------------
! CAM and WACCM gravity wave parameterizations were merged by Sean Patrick
! Santos in Summer 2013, and at the same time, gw_drag was split into
! various modules. This is the CAM interface and driver module. The below
! notes are for the old CAM and WACCM versions of gw_drag.
!--------------------------------------------------------------------------
! This file came from wa17 and was modified by Fabrizio: 07-02-2004
! Standard gw_drag with modification (6) of latitude profile of gw spectrum
!--------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!--------------------------------------------------------------------------
!
! Nov 2020  O. Guba Option for energy fix in GWD
!
!--------------------------------------------------------------------------

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, pver
  use hycoef,        only: hyai, hybi, hyam, hybm, etamid
  use constituents,  only: pcnst
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use spmd_utils,    only: masterproc
  use cam_history,   only: outfld, hist_fld_active
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun

  use ref_pres,      only: do_molec_diff, ntop_molec, nbot_molec
  use physconst,     only: cpair, rh2o, zvir, pi, rearth, r_universal!zvir is the ep1 in wrf,rearth is the radius of earth(m),r_universal is the gas constant

  ! These are the actual switches for different gravity wave sources.
  ! The orographic control switches are also here
  use phys_control,  only: use_gw_oro, use_gw_front, use_gw_convect, use_gw_energy_fix, use_od_ls, use_od_bl, use_od_ss
  use od_common,     only: od_ls_ncleff, od_bl_ncd, od_ss_sncleff

! Typical module header
  implicit none
  private
  save

!
! PUBLIC: interfaces
!
  public :: gw_drag_readnl           ! Read namelist
  public :: gw_register              ! Register pbuf variables
  public :: gw_init                  ! Initialization
  public :: gw_tend                  ! interface to actual parameterization

!
! PRIVATE: Rest of the data and interfaces are private to this module
!
  ! Whether spectral waves are being used.
  logical :: do_spectral_waves

  real(r8), parameter :: unset_r8 = huge(1._r8)

  ! fcrit2 has been made a namelist variable to facilitate backwards
  ! compatibility with the CAM3 version of this parameterization.  In CAM3,
  ! fcrit2=0.5
  real(r8) :: fcrit2 = unset_r8   ! critical froude number squared

  ! Maximum wave number and width of spectrum bins.
  integer :: pgwv = -1
  real(r8) :: dc = unset_r8

  integer :: kbotbg      ! interface of gwd source
  ! Top level for gravity waves.
  integer, parameter :: ktop = 0

  ! Frontogenesis function critical threshold.
  real(r8) :: frontgfc = unset_r8

  ! Tendency efficiencies.
  real(r8) :: effgw_oro = unset_r8      ! Orographic waves.
  real(r8) :: effgw_cm = unset_r8       ! Waves from C&M scheme.
  real(r8) :: effgw_beres = unset_r8    ! Waves from Beres scheme.

  ! Effective horizontal wave number (100 km wavelength).
  real(r8), parameter :: kwv = 6.28e-5_r8

  ! Background stress source strength.
  real(r8) :: taubgnd = unset_r8

  ! Convective heating rate conversion factor, default is 20.0_r8
  real(r8) :: gw_convect_hcf = 20.0_r8

  ! Scaling factor for the heating depth
  real(r8) :: hdepth_scaling_factor = 1

  ! Whether or not to enforce an upper boundary condition of tau = 0.
  ! (Like many variables, this is only here to hold the value between
  ! the readnl phase and the init phase of the CAM physics; only gw_common
  ! should actually use it.)
  logical :: tau_0_ubc = .false.

  ! Beres parameterization array.
  ! Max heating depth value.
  integer, parameter :: maxh = 20
  ! Max value for mean wind in heating.
  integer, parameter :: maxuh = 40
  ! File to read source spectra from.
  character(len=256) :: gw_drag_file = 'Beres04_file'

  ! Indices into pbuf
  integer :: kvt_idx      = -1
  integer :: ttend_dp_idx = -1
  integer :: frontgf_idx  = -1
  integer :: frontga_idx  = -1

  ! Prefixes for history field names
  character(len=1), parameter :: cm_pf = " "
  character(len=1), parameter :: beres_pf = "B"

  ! namelist 
  logical          :: history_amwg                   ! output the variables used by the AMWG diag package

  logical  :: use_gw_convect_old          ! switch to enable legacy behavior
  real(r8) :: gw_convect_plev_src_wind    ! reference pressure level for source wind for convective GWD [Pa]
  real(r8) :: gw_convect_hdepth_min       ! minimum hdepth for for convective GWD spectrum lookup table [km]
  real(r8) :: gw_convect_storm_speed_min  ! minimum convective storm speed for convective GWD           [m/s]
  
!==========================================================================
contains
!==========================================================================

subroutine gw_drag_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'gw_drag_readnl'

  ! More specific name for dc to prevent a name clash or confusion in the
  ! namelist.
  real(r8) :: gw_dc = unset_r8

  namelist /gw_drag_nl/ pgwv, gw_dc, tau_0_ubc, effgw_beres, effgw_cm, &
      effgw_oro, fcrit2, frontgfc, gw_drag_file, taubgnd, gw_convect_hcf, &
      hdepth_scaling_factor, gw_convect_hdepth_min, &
      gw_convect_storm_speed_min, gw_convect_plev_src_wind, &
      use_gw_convect_old
  !----------------------------------------------------------------------

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
  call mpibcast(pgwv,        1, mpiint, 0, mpicom)
  call mpibcast(gw_dc,       1, mpir8,  0, mpicom)
  call mpibcast(tau_0_ubc,   1, mpilog, 0, mpicom)
  call mpibcast(effgw_beres, 1, mpir8,  0, mpicom)
  call mpibcast(effgw_cm,    1, mpir8,  0, mpicom)
  call mpibcast(effgw_oro,   1, mpir8,  0, mpicom)
  call mpibcast(fcrit2,      1, mpir8,  0, mpicom)
  call mpibcast(frontgfc,    1, mpir8,  0, mpicom)
  call mpibcast(taubgnd,     1, mpir8,  0, mpicom)
  call mpibcast(gw_drag_file, len(gw_drag_file), mpichar, 0, mpicom)
  call mpibcast(gw_convect_hcf,             1, mpir8,  0, mpicom)
  call mpibcast(hdepth_scaling_factor,      1, mpir8,  0, mpicom)
  call mpibcast(gw_convect_hdepth_min,      1, mpir8,  0, mpicom)
  call mpibcast(gw_convect_storm_speed_min, 1, mpir8,  0, mpicom)
  call mpibcast(gw_convect_plev_src_wind,   1, mpir8,  0, mpicom)
  call mpibcast(use_gw_convect_old,         1, mpilog, 0, mpicom)
#endif

  dc = gw_dc

  ! Check if pgwv was set.
  if (pgwv < 0) then
     call endrun('gw_drag_readnl: pgwv must be set via the namelist and &
          &non-negative.')
  end if

  ! Check if dc was set.
  if (dc == unset_r8) then
     call endrun('gw_drag_readnl: dc must be set via the namelist')
  end if

  ! Check if fcrit2 was set.
  if (fcrit2 == unset_r8) then
     call endrun('gw_drag_readnl: fcrit2 must be set via the namelist')
  end if

end subroutine gw_drag_readnl

!==========================================================================

subroutine gw_register()
  use od_common,  only: oro_drag_register

  call oro_drag_register()

end subroutine gw_register

!==========================================================================

subroutine gw_init(pbuf2d)
  !-----------------------------------------------------------------------
  ! Time independent initialization for multiple gravity wave
  ! parameterization.
  !-----------------------------------------------------------------------

  use cam_history,      only: addfld, horiz_only, add_default
  use interpolate_data, only: lininterp
  use phys_control,     only: phys_getopts
  use physics_buffer,   only: pbuf_get_index, physics_buffer_desc

  use ref_pres,   only: pref_edge
  use physconst,  only: gravit, rair

  use gw_common,  only: gw_common_init, orographic_only
  use gw_oro,     only: gw_oro_init
  use gw_front,   only: gw_front_init
  use gw_convect, only: gw_convect_init

  use od_common,  only: oro_drag_init
  !------------------------------Arguments--------------------------------
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  !---------------------------Local storage-------------------------------

  integer :: l, k

  ! Reference phase speed spectrum
  real(r8) :: cref(-pgwv:pgwv)

  ! Index for levels at specific pressures.
  integer :: kfront

  ! output tendencies and state variables for CAM4 temperature,
  ! water vapor, cloud ice and cloud liquid budgets.
  logical :: history_budget
  ! output history file number for budget fields
  integer :: history_budget_histfile_num
  ! output variables of interest in WACCM runs
  logical :: history_waccm

  ! Interpolated Newtonian cooling coefficients.
  real(r8) :: alpha(0:pver)

  ! Levels of pre-calculated Newtonian cooling (1/day).
  integer, parameter :: nalph=66
  real(r8) :: alpha0(nalph) = [ &
       1.896007_r8   , 1.196965_r8   , 0.7251356_r8  ,  0.6397463_r8   , &
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
       0.00000000_r8 , 0.00000000_r8 &
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real(r8) :: palph(nalph) = [ &
       5.11075e-6_r8 , 9.8269e-6_r8  , 1.620185e-5_r8, 2.671225e-5_r8  , &
       4.4041e-5_r8  , 7.261275e-5_r8, 1.19719e-4_r8 , 1.9738e-4_r8    , &
       3.254225e-4_r8, 5.365325e-4_r8, 8.846025e-4_r8, 0.001458458_r8  , &
       0.002404575_r8, 0.00397825_r8 , 0.006556825_r8, 0.01081382_r8   , &
       0.017898_r8   , 0.02955775_r8 , 0.04873075_r8 , 0.07991075_r8   , &
       0.1282732_r8  , 0.19812_r8    , 0.292025_r8   , 0.4101675_r8    , &
       0.55347_r8    , 0.73048_r8    , 0.9559475_r8  , 1.244795_r8     , &
       1.61285_r8    , 2.079325_r8   , 2.667425_r8   , 3.404875_r8     , &
       4.324575_r8   , 5.4654_r8     , 6.87285_r8    , 8.599725_r8     , &
       10.70705_r8   , 13.26475_r8   , 16.35175_r8   , 20.05675_r8     , &
       24.479_r8     , 29.728_r8     , 35.92325_r8   , 43.19375_r8     , &
       51.6775_r8    , 61.5205_r8    , 72.8745_r8    , 85.65715_r8     , &
       100.5147_r8   , 118.2503_r8   , 139.1154_r8   , 163.6621_r8     , &
       192.5399_r8   , 226.5132_r8   , 266.4812_r8   , 313.5013_r8     , &
       368.818_r8    , 433.8952_r8   , 510.4553_r8   , 600.5242_r8     , &
       696.7963_r8   , 787.7021_r8   , 867.1607_r8   , 929.6489_r8     , &
       970.5548_r8   , 992.5561_r8 &
       ]

  ! Spectra for gravity waves from convective sources.
  real(r8) :: mfcc(maxh,-maxuh:maxuh,-pgwv:pgwv)

  ! Allow reporting of error messages.
  character(len=128) :: errstring

  !-----------------------------------------------------------------------

  call oro_drag_init(pbuf2d)

  ! Set model flags.
  do_spectral_waves = (pgwv > 0 .and. (use_gw_front .or. use_gw_convect))
  orographic_only = (use_gw_oro .and. .not. do_spectral_waves)

  if (do_molec_diff) then
     kvt_idx     = pbuf_get_index('kvt')
  end if

  ! Set phase speeds
  cref = (/ (dc * l, l = -pgwv, pgwv) /)

  if (masterproc) then
     write(iulog,*) ' '
     write(iulog,*) "GW_DRAG: pgwv = ", pgwv
     do l = -pgwv, pgwv
        write (iulog,'(A,I4,A,F7.2)') "GW_DRAG: cref(",l,") = ",cref(l)
     enddo
     write(iulog,*) 'GW_DRAG: kwv = ', kwv
     write(iulog,*) 'GW_DRAG: fcrit2 = ', fcrit2
     write(iulog,*) ' '
  end if

  if (.not. orographic_only) then

     ! pre-calculated newtonian damping:
     !     * convert to 1/s
     !     * ensure it is not smaller than 1e-6
     !     * convert palph from hpa to pa

     do k=1,nalph
        alpha0(k) = alpha0(k) / 86400._r8
        alpha0(k) = max(alpha0(k), 1.e-6_r8)
        palph(k) = palph(k)*1.e2_r8
     end do

     ! interpolate to current vertical grid and obtain alpha

     call lininterp (alpha0  ,palph, nalph , alpha  , pref_edge , pver+1)
     if (masterproc) then
        write (iulog,*) 'gw_init: newtonian damping (1/day):'
        write (iulog,fmt='(a4,a12,a10)') ' k  ','  pref_edge      ', &
             '  alpha   '
        do k=0,pver
           write (iulog,fmt='(i4,1e12.5,1f10.2)') k,pref_edge(k+1), &
                alpha(k)*86400._r8
        end do
     end if

  else
     ! about 10 days
     alpha = 1.e-6_r8
  end if

  ! Determine the bounds of the background and orographic stress regions
  ! spectrum source at 500 mb
  kbotbg = maxloc(pref_edge, 1, (pref_edge < 50000._r8)) - 1

  if (masterproc) then
     write(iulog,*) 'KTOP    =',ktop
     write(iulog,*) 'KBOTBG  =',kbotbg
  end if

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_budget_out = history_budget, &
       history_budget_histfile_num_out = history_budget_histfile_num, &
       history_waccm_out = history_waccm, &
       history_amwg_out   = history_amwg  )

  ! Initialize subordinate modules.
  call gw_common_init(pver, pgwv, dc, cref, do_molec_diff, tau_0_ubc, &
       nbot_molec, ktop, kbotbg, fcrit2, kwv, gravit, rair, alpha, &
       errstring)
  if (trim(errstring) /= "") call endrun("gw_common_init: "//errstring)

  if (use_gw_oro.or.&
      use_od_ls.or.&
      use_od_bl.or.&
      use_od_ss) then
     !
     if (use_gw_oro.and.effgw_oro == unset_r8) then
        call endrun("gw_drag_init: Orographic gravity waves enabled, &
             &but effgw_oro was not set.")
     end if
     !
     if (use_gw_oro.and.use_od_ls) then
        call endrun("gw_drag_init: Both orographic gravity waves schemes are turned on, &
                             &please turn one off by setting use_gw_oro or use_od_ls as .false.")
     end if
     !
     call gw_oro_init(errstring)
     if (trim(errstring) /= "") call endrun("gw_oro_init: "//errstring)

     ! Declare history variables for orographic term
     call addfld ('TTGWORO',(/ 'lev' /), 'A','K/s', &
          'T tendency - orographic gravity wave drag')
     call addfld ('UTGWORO',(/ 'lev' /), 'A','m/s2', &
          'U tendency - orographic gravity wave drag')
     call addfld ('VTGWORO',(/ 'lev' /), 'A','m/s2', &
          'V tendency - orographic gravity wave drag')
     call addfld ('TAUGWX',horiz_only,    'A','N/m2', &
          'Zonal gravity wave surface stress')
     call addfld ('TAUGWY',horiz_only,    'A','N/m2', &
          'Meridional gravity wave surface stress')
    if (use_od_ls.or.&
        use_od_bl.or.&
        use_od_ss) then
    !added for orographic drag
    call addfld ('DTAUX3_LS',(/'lev'/),'A','m/s2','U tendency - ls orographic drag')
    call addfld ('DTAUY3_LS',(/'lev'/),'A','m/s2','V tendency - ls orographic drag')
    call addfld ('DTAUX3_BL',(/'lev'/),'A','m/s2','U tendency - bl orographic drag')
    call addfld ('DTAUY3_BL',(/'lev'/),'A','m/s2','V tendency - bl orographic drag')
    call addfld ('DTAUX3_SS',(/'lev'/),'A','m/s2','U tendency - ss orographic drag')
    call addfld ('DTAUY3_SS',(/'lev'/),'A','m/s2','V tendency - ss orographic drag')
    call addfld ('DUSFC_LS',horiz_only,'A', 'N/m2', 'ls zonal oro surface stress')
    call addfld ('DVSFC_LS',horiz_only,'A', 'N/m2', 'ls merio oro surface stress')
    call addfld ('DUSFC_BL',horiz_only,'A', 'N/m2', 'bl zonal oro surface stress')
    call addfld ('DVSFC_BL',horiz_only,'A', 'N/m2', 'bl merio oro surface stress')
    call addfld ('DUSFC_SS',horiz_only,'A', 'N/m2', 'ss zonal oro surface stress')
    call addfld ('DVSFC_SS',horiz_only,'A', 'N/m2', 'ss merio oro surface stress')
    call add_default('DTAUX3_LS      ',    1,' ')
    call add_default('DTAUY3_LS      ',    1,' ')
    call add_default('DTAUX3_BL      ',    1,' ')
    call add_default('DTAUY3_BL      ',    1,' ')
    call add_default('DTAUX3_SS      ',    1,' ')
    call add_default('DTAUY3_SS      ',    1,' ')
    call add_default ('DUSFC_LS      ',    1,' ')
    call add_default ('DVSFC_LS      ',    1,' ')
    call add_default ('DUSFC_BL      ',    1,' ')
    call add_default ('DVSFC_BL      ',    1,' ')
    call add_default ('DUSFC_SS      ',    1,' ')
    call add_default ('DVSFC_SS      ',    1,' ')
    !added for orographic drag output
    endif

     if (history_amwg) then
        call add_default('TAUGWX  ', 1, ' ')
        call add_default('TAUGWY  ', 1, ' ')
     end if

     if (history_budget ) then
        call add_default('TTGWORO', history_budget_histfile_num, ' ')
     end if

     if (history_waccm) then
        call add_default('UTGWORO ', 1, ' ')
        call add_default('VTGWORO ', 1, ' ')
        call add_default('TAUGWX  ', 1, ' ')
        call add_default('TAUGWY  ', 1, ' ')
     end if

  end if

  if (do_spectral_waves) then
     if (use_gw_front) then

        frontgf_idx = pbuf_get_index('FRONTGF')
        frontga_idx = pbuf_get_index('FRONTGA')

        if (any(unset_r8 == &
             (/ effgw_cm, frontgfc, taubgnd /))) then
           call endrun("gw_drag_init: Frontogenesis enabled, but not &
                &all required namelist variables were set!")
        end if

        do k = 0, pver
           ! Check frontogenesis at 600 hPa.
           if (pref_edge(k+1) < 60000._r8) kfront = k+1
        end do

        if (masterproc) then
           write (iulog,*) 'KFRONT  =',kfront
           write(iulog,*) ' '
           write(iulog,*) 'gw_init: gw spectrum taubgnd, ', &
                'effgw_cm = ',taubgnd, effgw_cm
           write(iulog,*) ' '
        end if

        call gw_front_init(taubgnd, frontgfc, kfront, errstring)
        if (trim(errstring) /= "") &
             call endrun("gw_front_init: "//errstring)

        ! Output for gravity waves from frontogenesis.
        call gw_spec_addflds(prefix=cm_pf, scheme="C&M", &
             history_waccm=history_waccm)

        call addfld ('FRONTGF', (/ 'lev' /), 'A', 'K^2/M^2/S', &
             'Frontogenesis function at gws src level')
        call addfld ('FRONTGFA', (/ 'lev' /), 'A', 'K^2/M^2/S', &
             'Frontogenesis function at gws src level')

        if (history_waccm) then
           call add_default('FRONTGF', 1, ' ')
           call add_default('FRONTGFA', 1, ' ')
        end if

     end if

     if (use_gw_convect) then

        ttend_dp_idx    = pbuf_get_index('TTEND_DP')

        ! Initialization of Beres' parameterization parameters
        call gw_init_beres(mfcc)
        call gw_convect_init(gw_convect_plev_src_wind, mfcc, errstring)
        if (trim(errstring) /= "") call endrun("gw_convect_init: "//errstring)

        ! Output for gravity waves from the Beres scheme.
        call gw_spec_addflds(prefix=beres_pf, scheme="Beres", &
             history_waccm=history_waccm)

        call addfld ('NETDT',(/ 'lev' /), 'A','K/s', &
             'Net heating rate')
        call addfld ('MAXQ0',horiz_only  ,  'A','K/day', &
             'Max column heating rate')
        call addfld ('HDEPTH',horiz_only,    'A','km', &
             'Heating Depth')

        if (history_waccm) then
           call add_default('NETDT    ', 1, ' ')
           call add_default('HDEPTH   ', 1, ' ')
           call add_default('MAXQ0    ', 1, ' ')
        end if

     end if

     call addfld ('EKGWSPEC' ,(/ 'ilev' /), 'A','M2/S', &
          'effective Kzz due to gw spectrum')

     if (history_waccm) then
        call add_default('EKGWSPEC', 1, ' ')
     end if

  end if

  ! Total temperature tendency output.
  call addfld ('TTGW',(/ 'lev' /), 'A','K/s', &
       'T tendency - gravity wave drag')

  if (masterproc) then
     write (iulog,*) 'using GW energy fix   =',use_gw_energy_fix
  end if

  if ( history_budget ) then
     call add_default ('TTGW', history_budget_histfile_num, ' ')
  end if

end subroutine gw_init

!==========================================================================

subroutine gw_init_beres(mfcc)

  use ioFileMod,        only: getfil
#if ( defined SPMD )
  use mpishorthand
#endif

  use netcdf

  real(r8), intent(out) :: mfcc(maxh,-maxuh:maxuh,-pgwv:pgwv)

  ! VARIABLES NEEDED TO READ IN TABLE OF SPECTRA FROM FILE

  ! netCDF file id, variable id, and error code.
  integer :: ncid, varid, ncstat
  ! The number of gravity waves in the file (i.e. the file has a gravity
  ! wave spectrum from -ngwv_file to ngwv_file). This should be determined
  ! from the file itself, rather than hard-coded.
  integer, parameter :: ngwv_file = 40
  character(len=256) :: gw_drag_file_loc ! local filepath of gw_drag_file

  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  if (masterproc) then

     call getfil(gw_drag_file, gw_drag_file_loc)
     ncstat = NF90_OPEN (gw_drag_file_loc,0,ncid)

     if (ncstat .ne. 0) then
        write(iulog,*) 'Error reading in netcdf file ',gw_drag_file,'.  ',&
             NF90_STRERROR(ncstat)
        write(iulog,*) 'Check that the namelist variable gw_drag_file is &
             &correct.'
        call endrun
     endif

     ncstat = NF90_INQ_VARID (ncid,'mfcc',varid)

     if (ncstat .ne. 0) then
        write(iulog,*) 'Error reading data from ',gw_drag_file,'.  ',&
             NF90_STRERROR(ncstat)
        call endrun
     endif

     ncstat = NF90_GET_VAR(ncid,varid,mfcc,start=[1,1,ngwv_file-pgwv+1])
     ncstat = NF90_CLOSE(ncid)

     write(iulog,*) 'Read-in source spectra from file'
     write(iulog,*) 'MFCC=',maxval(mfcc),minval(mfcc)

  endif
  !     Broadcast results
#ifdef SPMD
  call mpibcast (mfcc, maxh*(2*pgwv+1)*(2*maxuh+1),mpir8,0,mpicom)
#endif

end subroutine gw_init_beres

!==========================================================================

subroutine gw_tend(state, sgh, pbuf, dt, ptend, cam_in)
  !-----------------------------------------------------------------------
  ! Interface for multiple gravity wave drag parameterization.
  !-----------------------------------------------------------------------
  use physics_types,  only: physics_state_copy, set_dry_to_wet
  use constituents,   only: cnst_type
  use co2_cycle,      only: co2_cycle_set_cnst_type
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field
  use camsrfexch, only: cam_in_t
  ! Location-dependent cpair
  use physconst,  only: cpairv
  use od_common,  only: oro_drag_interface
  use gw_common,  only: gw_prof, momentum_energy_conservation, &
       gw_drag_prof
  use gw_oro,     only: gw_oro_src
  use gw_front,   only: gw_cm_src
  use gw_convect, only: gw_beres_src
  use dycore,     only: dycore_is
  use phys_grid,  only: get_rlat_all_p
  use physconst,  only: gravit,rair
  !------------------------------Arguments--------------------------------
  type(physics_state), intent(in) :: state      ! physics state structure
  ! Standard deviation of orography.
  real(r8), intent(in) :: sgh(pcols)
  type(physics_buffer_desc), pointer :: pbuf(:) ! Physics buffer
  real(r8), intent(in) :: dt                    ! time step
  ! Parameterization net tendencies.
  type(physics_ptend), intent(out):: ptend
  type(cam_in_t), intent(in) :: cam_in
  !locally added gw and bl drag
  real(r8) :: dtaux3_ls(pcols,pver)
  real(r8) :: dtauy3_ls(pcols,pver)
  real(r8) :: dtaux3_bl(pcols,pver)
  real(r8) :: dtauy3_bl(pcols,pver)
  real(r8) :: dtaux3_ss(pcols,pver)
  real(r8) :: dtauy3_ss(pcols,pver)
  real(r8) :: dummx3_fd(pcols,pver)
  real(r8) :: dummy3_fd(pcols,pver)
  !
  real(r8) :: dusfc_ls(pcols)
  real(r8) :: dvsfc_ls(pcols)
  real(r8) :: dusfc_bl(pcols)
  real(r8) :: dvsfc_bl(pcols)
  real(r8) :: dusfc_ss(pcols)
  real(r8) :: dvsfc_ss(pcols)
  real(r8) :: dummx_fd(pcols)
  real(r8) :: dummy_fd(pcols)
  !
  real(r8) :: dx(pcols),dy(pcols)

  !---------------------------Local storage-------------------------------

  type(physics_state) :: state1     ! Local copy of state variable

  integer :: lchnk                  ! chunk identifier
  integer :: ncol                   ! number of atmospheric columns

  integer :: i, k                   ! loop indices

  real(r8) :: ttgw(state%ncol,pver) ! temperature tendency
  real(r8) :: utgw(state%ncol,pver) ! zonal wind tendency
  real(r8) :: vtgw(state%ncol,pver) ! meridional wind tendency

  real(r8) :: ni(state%ncol,0:pver) ! interface Brunt-Vaisalla frequency
  real(r8) :: nm(state%ncol,pver)   ! midpoint Brunt-Vaisalla frequency
  real(r8) :: rhoi(state%ncol,0:pver)     ! interface density
  real(r8) :: tau(state%ncol,-pgwv:pgwv,0:pver)  ! wave Reynolds stress
  real(r8) :: tau0x(state%ncol)     ! c=0 sfc. stress (zonal)
  real(r8) :: tau0y(state%ncol)     ! c=0 sfc. stress (meridional)
  real(r8) :: ti(state%ncol,0:pver) ! interface temperature
  real(r8) :: ubi(state%ncol,0:pver)! projection of wind at interfaces
  real(r8) :: ubm(state%ncol,pver)  ! projection of wind at midpoints
  real(r8) :: xv(state%ncol)        ! unit vector of source wind (x)
  real(r8) :: yv(state%ncol)        ! unit vector of source wind (y)

  integer :: m                      ! dummy integers
  real(r8) :: qtgw(state%ncol,pver,pcnst) ! constituents tendencies

  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8) :: taucd(state%ncol,0:pver,4)

  ! gravity wave wind tendency for each wave
  real(r8) :: gwut(state%ncol,pver,-pgwv:pgwv)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8) :: dttdf(state%ncol,pver)
  real(r8) :: dttke(state%ncol,pver)

  ! spectrum phase speeds for each column
  real(r8) :: c(state%ncol,-pgwv:pgwv)

  ! pbuf fields
  ! Molecular diffusivity
  real(r8), pointer :: kvt_in(:,:)
  real(r8) :: kvtt(state%ncol,0:pver)

  ! Frontogenesis
  real(r8), pointer :: frontgf(:,:)
  real(r8), pointer :: frontga(:,:)

  ! Temperature change due to deep convection.
  real(r8), pointer, dimension(:,:) :: ttend_dp

  ! Indices of gravity wave source and lowest level where wind tendencies
  ! are allowed.
  integer :: src_level(state%ncol)
  integer :: tend_level(state%ncol)

  real(r8) :: hdepth(state%ncol)         ! heating depth array
  real(r8) :: maxq0(state%ncol)          ! maximum heating rate array

  ! effective gw diffusivity at interfaces needed for output
  real(r8) :: egwdffi(state%ncol,0:pver)
  ! sum from the two types of spectral GW
  real(r8) :: egwdffi_tot(state%ncol,0:pver)

  ! Which constituents are being affected by diffusion.
  logical  :: lq(pcnst)

  ! Contiguous copies of state arrays.
  real(r8) :: dse(state%ncol,pver)
  real(r8) :: t(state%ncol,pver)
  real(r8) :: u(state%ncol,pver)
  real(r8) :: v(state%ncol,pver)
  real(r8) :: q(state%ncol,pver,pcnst)
  real(r8) :: pmid(state%ncol,pver)
  real(r8) :: pint(state%ncol,pver+1)
  real(r8) :: piln(state%ncol,pver+1)
  real(r8) :: dpm(state%ncol,pver)
  real(r8) :: rdpm(state%ncol,pver)
  real(r8) :: zm(state%ncol,pver)

  real(r8) :: dE(state%ncol)

  ! local override option for constituents cnst_type
  character(len=3), dimension(pcnst) :: cnst_type_loc

  logical :: do_latitude_taper

  !------------------------------------------------------------------------

  ! Make local copy of input state.
  call physics_state_copy(state, state1)

  ! constituents are all treated as wet mmr
  ! don't convert co2 tracers to wet mixing ratios
  cnst_type_loc(:) = cnst_type(:)
  call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
  call set_dry_to_wet(state1, cnst_type_loc)

  lchnk = state1%lchnk
  ncol  = state1%ncol

  dse = state1%s(:ncol,:)
  t = state1%t(:ncol,:)
  u = state1%u(:ncol,:)
  v = state1%v(:ncol,:)
  q = state1%q(:ncol,:,:)
  pmid = state1%pmid(:ncol,:)
  pint = state1%pint(:ncol,:)
  piln = state1%lnpint(:ncol,:)
  dpm = state1%pdel(:ncol,:)
  rdpm = state1%rpdel(:ncol,:)
  zm = state1%zm(:ncol,:)

  lq = .true.
  call physics_ptend_init(ptend, state1%psetcols, "Grav_wave_drag", &
       ls=.true., lu=.true., lv=.true., lq=lq)

  ! Profiles of background state variables
  call gw_prof(ncol, cpair, &
       t         , pmid         , pint      , rhoi      , ti        , &
       nm        , ni)

  if (do_molec_diff) then
     !--------------------------------------------------------
     ! Initialize and calculate local molecular diffusivity
     !--------------------------------------------------------

     call pbuf_get_field(pbuf, kvt_idx, kvt_in)  ! kvt_in(1:pcols,1:pver+1)

     ! Set kvtt from pbuf field; kvtt still needs a factor of 1/cpairv.
     kvtt(:,0:pver) = kvt_in(:ncol,:)

     ! Use linear extrapolation of cpairv to top interface.
     kvtt(:,ntop_molec-1) = kvtt(:,ntop_molec-1) / &
          (1.5_r8*cpairv(:ncol,ntop_molec,lchnk) - &
          0.5_r8*cpairv(:ncol,ntop_molec+1,lchnk))

     ! Interpolate cpairv to other interfaces.
     do k = ntop_molec, nbot_molec-1
        kvtt(:,k) = kvtt(:,k) / &
             (cpairv(:ncol,k+1,lchnk)+cpairv(:ncol,k,lchnk)) * 2._r8
     enddo

  end if

  !------------------------------------------------------------------------
  ! Non-orographic background gravity wave spectra
  !------------------------------------------------------------------------

  if (do_spectral_waves) then

     egwdffi_tot = 0._r8

     if (use_gw_convect) then
        !------------------------------------------------------------------
        ! Convective gravity waves (Beres scheme)
        !------------------------------------------------------------------

        ! Set up heating
        call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)

        ! Determine wave sources for Beres04 scheme
        call gw_beres_src(ncol, pgwv, state1%lat(:ncol), u, v, ttend_dp, &
             zm, src_level, tend_level, tau, ubm, ubi, xv, yv, c, &
             hdepth, maxq0, gw_convect_hcf, hdepth_scaling_factor, &
             gw_convect_hdepth_min, gw_convect_storm_speed_min, &
             use_gw_convect_old)

        do_latitude_taper = .false.

        ! Solve for the drag profile with Beres source spectrum.
        call gw_drag_prof(ncol, pgwv, src_level, tend_level, do_latitude_taper, dt, &
             state1%lat(:ncol), t,    ti, pmid, pint, dpm,   rdpm, &
             piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
             effgw_beres, c,   kvtt, q,  dse,  tau,  utgw,  vtgw, &
             ttgw, qtgw,  taucd,     egwdffi,  gwut, dttdf, dttke)

        !  add the diffusion coefficients
        do k = 0, pver
           egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
        end do

        ! Store constituents tendencies
        do m=1, pcnst
           do k = 1, pver
              ptend%q(:ncol,k,m) = qtgw(:,k,m)
           end do
        end do

        ! add the momentum tendencies to the output tendency arrays
        do k = 1, pver
           ptend%u(:ncol,k) = utgw(:,k)
           ptend%v(:ncol,k) = vtgw(:,k)
           ptend%s(:ncol,k) = ttgw(:,k)
        end do


        ! C.-C. Chen, momentum & energy conservation
        call momentum_energy_conservation(ncol, tend_level, dt, taucd, &
             pint, dpm, u, v, ptend%u, ptend%v, ptend%s, utgw, vtgw, ttgw)

        call gw_spec_outflds(beres_pf, lchnk, ncol, pgwv, c, u, v, &
             xv, yv, gwut, dttdf, dttke, tau(:,:,1:), utgw, vtgw, taucd, state)

        ! Note: This is probably redundant, because ZMDT is already being
        ! output...
        call outfld ('NETDT', ttend_dp, pcols, lchnk)

        ! Diagnostic outputs.
        call outfld ('HDEPTH', hdepth, ncol, lchnk)
        call outfld ('MAXQ0', maxq0, ncol, lchnk)

     end if

     if (use_gw_front) then
        !------------------------------------------------------------------
        ! Frontally generated gravity waves
        !------------------------------------------------------------------

        ! Get frontogenesis physics buffer fields set by dynamics.
        call pbuf_get_field(pbuf, frontgf_idx, frontgf)
        call pbuf_get_field(pbuf, frontga_idx, frontga)

        ! Determine the wave source for C&M background spectrum
        call gw_cm_src(ncol, pgwv, kbotbg, u, v, frontgf, &
             src_level, tend_level, tau, ubm, ubi, xv, yv, c)

        if (dycore_is('UNSTRUCTURED')) then
          do_latitude_taper = .false.
        else
          do_latitude_taper = .true.
        end if

        ! Solve for the drag profile with C&M source spectrum.
        call gw_drag_prof(ncol, pgwv, src_level, tend_level, do_latitude_taper, dt, &
             state1%lat(:ncol), t,    ti, pmid, pint, dpm,   rdpm, &
             piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
             effgw_cm,    c,   kvtt, q,  dse,  tau,  utgw,  vtgw, &
             ttgw, qtgw,  taucd,     egwdffi,  gwut, dttdf, dttke)

        !  add the diffusion coefficients
        do k = 0, pver
           egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
        end do

        !Add the constituent tendencies
        do m=1, pcnst
           do k = 1, pver
              ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
           end do
        end do


        ! add the momentum tendencies to the output tendency arrays
        do k = 1, pver
           ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
           ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
           ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
        end do

        ! C.-C. Chen, momentum & energy conservation
        call momentum_energy_conservation(ncol, tend_level, dt, taucd, &
             pint, dpm, u, v, ptend%u, ptend%v, ptend%s, utgw, vtgw, ttgw)

        call gw_spec_outflds(cm_pf, lchnk, ncol, pgwv, c, u, v, &
             xv, yv, gwut, dttdf, dttke, tau(:,:,1:), utgw, vtgw, taucd, state)

        call outfld ('FRONTGF', frontgf, pcols, lchnk)
        call outfld ('FRONTGFA', frontga, pcols, lchnk)

     end if

     call outfld ('EKGWSPEC', egwdffi_tot , ncol, lchnk)

  end if

  if (use_gw_oro) then
     !---------------------------------------------------------------------
     ! Orographic stationary gravity waves
     !---------------------------------------------------------------------
     ! Determine the orographic wave source
     call gw_oro_src(ncol, &
          u, v, t, sgh(:ncol), pmid, pint, dpm, zm, nm, &
          src_level, tend_level, tau, ubm, ubi, xv, yv, c)

     do_latitude_taper = .false.

     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, 0, src_level, tend_level, do_latitude_taper, dt, &
          state1%lat(:ncol), t,    ti, pmid, pint, dpm,   rdpm, &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw_oro,   c,   kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw,  taucd,     egwdffi,  gwut(:,:,0:0), dttdf, dttke)
  endif
  !
  if ( use_od_ls .or. use_od_bl .or. use_od_ss) then
     utgw=0.0_r8
     vtgw=0.0_r8
     ttgw=0.0_r8
     call oro_drag_interface(state,cam_in,sgh,pbuf,dt,nm,&
                             use_od_ls,use_od_bl,use_od_ss,.false.,&
                             od_ls_ncleff,od_bl_ncd,od_ss_sncleff,&
                             utgw,vtgw,ttgw,&
                             dtaux3_ls=dtaux3_ls,dtauy3_ls=dtauy3_ls,&
                             dtaux3_bl=dtaux3_bl,dtauy3_bl=dtauy3_bl,&
                             dtaux3_ss=dtaux3_ss,dtauy3_ss=dtauy3_ss,&
                             dtaux3_fd=dummx3_fd,dtauy3_fd=dummy3_fd,&
                             dusfc_ls=dusfc_ls,dvsfc_ls=dvsfc_ls,&
                             dusfc_bl=dusfc_bl,dvsfc_bl=dvsfc_bl,&
                             dusfc_ss=dusfc_ss,dvsfc_ss=dvsfc_ss,&
                             dusfc_fd=dummx_fd,dvsfc_fd=dummy_fd)
  endif
  !
  ! Add the orographic tendencies to the spectrum tendencies
  ! Compute the temperature tendency from energy conservation
  ! (includes spectrum).
  ! both old and new gwd scheme will add the tendency to circulation
  !
  if (use_gw_oro.or.&
      use_od_ls .or.&
      use_od_bl .or.&
      use_od_ss) then
     if(.not. use_gw_energy_fix) then
        !original
        do k = 1, pver
           utgw(:,k) = utgw(:,k) * cam_in%landfrac(:ncol)
           ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
           vtgw(:,k) = vtgw(:,k) * cam_in%landfrac(:ncol)
           ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
           ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k) &
                            -(ptend%u(:ncol,k) * (u(:,k) + ptend%u(:ncol,k)*0.5_r8*dt) &
                             +ptend%v(:ncol,k) * (v(:,k) + ptend%v(:ncol,k)*0.5_r8*dt))
           ttgw(:,k) = ttgw(:,k) &
                            -(ptend%u(:ncol,k) * (u(:,k) + ptend%u(:ncol,k)*0.5_r8*dt) &
                             +ptend%v(:ncol,k) * (v(:,k) + ptend%v(:ncol,k)*0.5_r8*dt))
           ttgw(:,k) = ttgw(:,k) / cpairv(:ncol, k, lchnk)
        end do
    else
        do k = 1, pver
           utgw(:,k) = utgw(:,k) * cam_in%landfrac(:ncol)
           ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
           vtgw(:,k) = vtgw(:,k) * cam_in%landfrac(:ncol)
           ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
           ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
        enddo

        dE = 0.0
        do k = 1, pver
           dE(:ncol) = dE(:ncol) &
                     - dpm(:ncol,k)*(ptend%u(:ncol,k) * (u(:ncol,k)+ptend%u(:ncol,k)*0.5_r8*dt) &
                                    +ptend%v(:ncol,k) * (v(:ncol,k)+ptend%v(:ncol,k)*0.5_r8*dt) &
                                    +ptend%s(:ncol,k) )
        enddo
        dE(:ncol)=dE(:ncol) / (pint(:ncol,pver+1) - pint(:ncol,1))

        do k = 1, pver
           ptend%s(:ncol,k) = ptend%s(:ncol,k) + dE(:ncol)
           ttgw(:ncol,k) = ( ttgw(:ncol,k) + dE(:ncol) ) / cpairv(:ncol, k, lchnk)
        enddo
     endif ! gw energy fix 2

     do m = 1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Write output fields to history file
     call outfld('UTGWORO', utgw,  ncol, lchnk)
     call outfld('VTGWORO', vtgw,  ncol, lchnk)
     call outfld('TTGWORO', ttgw,  ncol, lchnk)
     !
     if (use_gw_oro) then
     !old gwd scheme
     tau0x = tau(:,0,pver) * xv * effgw_oro
     tau0y = tau(:,0,pver) * yv * effgw_oro
     call outfld('TAUGWX', tau0x, ncol, lchnk)
     call outfld('TAUGWY', tau0y, ncol, lchnk)
     endif
     !
     call outfld('SGH   ',   sgh,pcols, lchnk)
     !
     if (use_od_ls.or.&
         use_od_bl.or.&
         use_od_ss) then
     call outfld ('DTAUX3_LS', dtaux3_ls,  pcols, lchnk)
     call outfld ('DTAUY3_LS', dtauy3_ls,  pcols, lchnk)
     call outfld ('DTAUX3_BL', dtaux3_bl,  pcols, lchnk)
     call outfld ('DTAUY3_BL', dtauy3_bl,  pcols, lchnk)
     call outfld ('DTAUX3_SS', dtaux3_ss,  pcols, lchnk)
     call outfld ('DTAUY3_SS', dtauy3_ss,  pcols, lchnk)
     call outfld ('DUSFC_LS', dusfc_ls,  pcols, lchnk)
     call outfld ('DVSFC_LS', dvsfc_ls,  pcols, lchnk)
     call outfld ('DUSFC_BL', dusfc_bl,  pcols, lchnk)
     call outfld ('DVSFC_BL', dvsfc_bl,  pcols, lchnk)
     call outfld ('DUSFC_SS', dusfc_ss,  pcols, lchnk)
     call outfld ('DVSFC_SS', dvsfc_ss,  pcols, lchnk)
     endif
     !
  end if

  ! Convert the tendencies for the dry constituents to dry air basis.
  do m = 1, pcnst
     if (cnst_type(m).eq.'dry') then
        do k = 1, pver
           do i = 1, ncol
              ptend%q(i,k,m) = ptend%q(i,k,m)*state1%pdel(i,k)/state1%pdeldry(i,k)
           end do
        end do
     end if
  end do

  ! Write total temperature tendency to history file
  call outfld ('TTGW', ptend%s/cpairv(:,:,lchnk),  pcols, lchnk)

end subroutine gw_tend

!==========================================================================

! Add all history fields for a gravity wave spectrum source.
subroutine gw_spec_addflds(prefix, scheme, history_waccm)
  use cam_history, only: addfld, horiz_only, add_default
  use gw_common,   only: cref

  !------------------------------Arguments--------------------------------

  ! One character prefix prepended to output fields.
  character(len=1), intent(in) :: prefix
  ! Gravity wave scheme name prepended to output field descriptions.
  character(len=*), intent(in) :: scheme
  ! Whether or not to call add_default for fields output by WACCM.
  logical, intent(in) :: history_waccm

  !---------------------------Local storage-------------------------------

  integer :: l
  ! 7 chars is enough for "-100.00"
  character(len=7)  :: fnum
  ! 10 chars is enough for "BTAUXSn32"
  character(len=10) :: dumc1x, dumc1y
  ! Allow 80 chars for description
  character(len=80) dumc2

  !-----------------------------------------------------------------------

  ! Overall wind tendencies.
  call addfld (trim(prefix)//'UTGWSPEC',(/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency - gravity wave spectrum')
  call addfld (trim(prefix)//'VTGWSPEC',(/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' V tendency - gravity wave spectrum')

  ! Wind tendencies broken across five spectral bins.
  call addfld (trim(prefix)//'UTEND1',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency   c < -40')
  call addfld (trim(prefix)//'UTEND2',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency  -40 < c < -15')
  call addfld (trim(prefix)//'UTEND3',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency  -15 < c <  15')
  call addfld (trim(prefix)//'UTEND4',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency   15 < c <  40')
  call addfld (trim(prefix)//'UTEND5',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency   40 < c ')

  ! Reynold's stress toward each cardinal direction, and net zonal stress.
  call addfld (trim(prefix)//'TAUE' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Eastward Reynolds stress')
  call addfld (trim(prefix)//'TAUW' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Westward Reynolds stress')
  call addfld (trim(prefix)//'TAUNET' , (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' E+W Reynolds stress')
  call addfld (trim(prefix)//'TAUN' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Northward Reynolds stress')
  call addfld (trim(prefix)//'TAUS' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Southward Reynolds stress')

  ! Momentum flux in each direction.
  call addfld (trim(prefix)//'EMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Eastward MF')
  call addfld (trim(prefix)//'WMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Westward MF')
  call addfld (trim(prefix)//'NMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Northward MF')
  call addfld (trim(prefix)//'SMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Southward MF')

  ! Temperature tendency terms.
  call addfld (trim(prefix)//'TTGWSDF' , (/ 'lev' /), 'A','K/s', &
       trim(scheme)//' t tendency - diffusion term')
  call addfld (trim(prefix)//'TTGWSKE' , (/ 'lev' /), 'A','K/s', &
       trim(scheme)//' t tendency - kinetic energy conversion term')

  ! Gravity wave source spectra by wave number.
  do l=-pgwv,pgwv
     ! String containing reference speed.
     write (fnum,fmt='(f7.2)') cref(l)

     dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
     dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)
     dumc2 = trim(scheme)//" tau at c= "//trim(fnum)//" m/s"
     call addfld (trim(dumc1x),(/ 'lev' /), 'A','Pa',dumc2)
     call addfld (trim(dumc1y),(/ 'lev' /), 'A','Pa',dumc2)

     ! add zonal source spectra on specific pressure levels
     call addfld (trim(dumc1x)//'_100mb',horiz_only, 'A','Pa',dumc2//" at 100 mbar pressure surface")
     call addfld (trim(dumc1x)//'_50mb' ,horiz_only, 'A','Pa',dumc2//" at 50 mbar pressure surface")
     call addfld (trim(dumc1x)//'_30mb' ,horiz_only, 'A','Pa',dumc2//" at 30 mbar pressure surface")
     call addfld (trim(dumc1x)//'_10mb' ,horiz_only, 'A','Pa',dumc2//" at 10 mbar pressure surface")

  end do

  if (history_waccm) then
     call add_default(trim(prefix)//'UTGWSPEC', 1, ' ')
     call add_default(trim(prefix)//'VTGWSPEC', 1, ' ')

     call add_default(trim(prefix)//'UTEND1', 1, ' ')
     call add_default(trim(prefix)//'UTEND2', 1, ' ')
     call add_default(trim(prefix)//'UTEND3', 1, ' ')
     call add_default(trim(prefix)//'UTEND4', 1, ' ')
     call add_default(trim(prefix)//'UTEND5', 1, ' ')

     call add_default(trim(prefix)//'TAUE', 1, ' ')
     call add_default(trim(prefix)//'TAUW', 1, ' ')
     call add_default(trim(prefix)//'TAUNET', 1, ' ')
     call add_default(trim(prefix)//'TAUN', 1, ' ')
     call add_default(trim(prefix)//'TAUS', 1, ' ')
  end if

end subroutine gw_spec_addflds

!==========================================================================

! Outputs for spectral waves.
subroutine gw_spec_outflds(prefix, lchnk, ncol, ngwv, c, u, v, xv, yv, &
     gwut, dttdf, dttke, tau, utgw, vtgw, taucd, state)

  use gw_common, only: west, east, south, north
  use interpolate_data, only: vertinterp

  ! One-character prefix prepended to output fields.
  character(len=1), intent(in) :: prefix
  ! Chunk and number of columns in the chunk.
  integer, intent(in) :: lchnk
  integer, intent(in) :: ncol
  ! Max (+/-) wavenumber in gravity wave spectrum.
  integer, intent(in) :: ngwv
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-pgwv:pgwv)
  ! Winds at cell midpoints.
  real(r8), intent(in) :: u(ncol,pver)
  real(r8), intent(in) :: v(ncol,pver)
  ! Unit vector in the direction of wind at source level.
  real(r8), intent(in) :: xv(ncol)
  real(r8), intent(in) :: yv(ncol)
  ! Wind tendency for each wave.
  real(r8), intent(in) :: gwut(ncol,pver,-ngwv:ngwv)
  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8) :: dttdf(ncol,pver)
  real(r8) :: dttke(ncol,pver)
  ! Wave Reynolds stress.
  real(r8), intent(in) :: tau(ncol,-pgwv:pgwv,pver)
  ! Zonal and meridional total wind tendencies.
  real(r8), intent(in) :: utgw(ncol,pver)
  real(r8), intent(in) :: vtgw(ncol,pver)
  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8), intent(in) :: taucd(ncol,0:pver,4)
  type(physics_state), intent(in) :: state      ! physics state structure

  ! Indices
  integer :: i, k, l
  integer :: ix(ncol, -ngwv:ngwv), iy(ncol, -ngwv:ngwv)
  integer :: iu(ncol), iv(ncol)

  ! Zonal wind tendency, broken up into five bins.
  real(r8) :: utb(ncol, pver, 5)
  ! Definition of the bin boundaries.
  real(r8), parameter :: bounds(4) = (/ -40._r8, -15._r8, &
       15._r8, 40._r8 /)

  ! Momentum flux in the four cardinal directions.
  real(r8) :: mf(ncol, pver, 4)

  ! Wave stress in zonal/meridional direction
  real(r8) :: taux(ncol,-ngwv:ngwv,pver)
  real(r8) :: tauy(ncol,-ngwv:ngwv,pver)

  ! Temporaries for output
  real(r8) :: dummyx(pcols,pver)
  real(r8) :: dummyy(pcols,pver)
  real(r8) :: dummmy_p_surf(pcols) ! data interpolated to a pressure surface
  ! Variable names
  character(len=10) :: dumc1x, dumc1y


  ! Accumulate wind tendencies binned according to phase speed.

  utb = 0._r8

  ! Find which output bin the phase speed corresponds to.
  ix = find_bin(c(:, -ngwv:ngwv))

  ! Put the wind tendency in that bin.
  do l = -ngwv, ngwv
     do k = 1, pver
        do i = 1, ncol
           utb(i,k,ix(i,l)) = utb(i,k,ix(i,l)) + gwut(i,k,l)
        end do
     end do
  end do

  ! Find just the zonal part.
  do l = 1, 5
     do k = 1, pver
        utb(:, k, l) = utb(:, k, l) * xv
     end do
  end do

  call outfld(trim(prefix)//'UTEND1', utb(:,:,1), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND2', utb(:,:,2), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND3', utb(:,:,3), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND4', utb(:,:,4), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND5', utb(:,:,5), ncol, lchnk)

  ! Output temperature tendencies due to diffusion and from kinetic energy.
  call outfld(trim(prefix)//'TTGWSDF', dttdf / cpair, ncol, lchnk)
  call outfld(trim(prefix)//'TTGWSKE', dttke / cpair, ncol, lchnk)


  ! Output tau broken down into zonal and meridional components.

  taux = 0._r8
  tauy = 0._r8

  ! Project c, and convert each component to a wavenumber index.
  ! These are mappings from the wavenumber index of tau to those of taux
  ! and tauy, respectively.
  do l=-ngwv,ngwv
     ix(:,l) = c_to_l(c(:,l)*xv)
     iy(:,l) = c_to_l(c(:,l)*yv)
  end do

  ! Find projection of tau.
  do k = 1, pver
     do l = -ngwv,ngwv
        do i = 1, ncol
           taux(i,ix(i,l),k) = taux(i,ix(i,l),k) &
                + abs(tau(i,l,k)*xv(i))
           tauy(i,iy(i,l),k) = tauy(i,iy(i,l),k) &
                + abs(tau(i,l,k)*yv(i))
        end do
     end do
  end do

  do l=-ngwv,ngwv

     dummyx(1:ncol,1:pver) = taux(1:ncol,l,1:pver)
     dummyy(1:ncol,1:pver) = tauy(1:ncol,l,1:pver)

     dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
     dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)

     call outfld(dumc1x,dummyx,ncol,lchnk)
     call outfld(dumc1y,dummyy,ncol,lchnk)

     if (hist_fld_active(trim(dumc1x)//'_100mb')) then
        call vertinterp(ncol, pcols, pver, state%pmid, 10000._r8, dummyx, dummmy_p_surf)
        call outfld(trim(dumc1x)//'_100mb', dummmy_p_surf, pcols, lchnk)
     end if
     if (hist_fld_active(trim(dumc1x)//'_50mb')) then
        call vertinterp(ncol, pcols, pver, state%pmid, 5000._r8, dummyx, dummmy_p_surf)
        call outfld(trim(dumc1x)//'_50mb', dummmy_p_surf, pcols, lchnk)
     end if
     if (hist_fld_active(trim(dumc1x)//'_30mb')) then
        call vertinterp(ncol, pcols, pver, state%pmid, 3000._r8, dummyx, dummmy_p_surf)
        call outfld(trim(dumc1x)//'_30mb', dummmy_p_surf, pcols, lchnk)
     end if
     if (hist_fld_active(trim(dumc1x)//'_10mb')) then
        call vertinterp(ncol, pcols, pver, state%pmid, 1000._r8, dummyx, dummmy_p_surf)
        call outfld(trim(dumc1x)//'_10mb', dummmy_p_surf, pcols, lchnk)
     end if

  enddo


  ! Output momentum flux in each cardinal direction.
  mf = 0._r8

  do k = 1, pver

     ! Convert wind speed components to wavenumber indices.
     iu = c_to_l(u(:,k))
     iv = c_to_l(v(:,k))

     ! Sum tau components in each cardinal direction.
     ! Split west/east and north/south based on whether wave speed exceeds
     ! wind speed.
     do l = -ngwv, ngwv

        where (iu > l)
           mf(:,k,west) = mf(:,k,west) + taux(:,l,k)
        elsewhere
           mf(:,k,east) = mf(:,k,east) + taux(:,l,k)
        end where

        where (iv > l)
           mf(:,k,south) = mf(:,k,south) + tauy(:,l,k)
        elsewhere
           mf(:,k,north) = mf(:,k,north) + tauy(:,l,k)
        end where

     end do

  end do

  call outfld(trim(prefix)//'WMF',mf(:,:,west),ncol,lchnk)
  call outfld(trim(prefix)//'EMF',mf(:,:,east),ncol,lchnk)
  call outfld(trim(prefix)//'SMF',mf(:,:,south),ncol,lchnk)
  call outfld(trim(prefix)//'NMF',mf(:,:,north),ncol,lchnk)

  ! Simple output fields written to history file.
  ! Total wind tendencies.
  call outfld (trim(prefix)//'UTGWSPEC', utgw , ncol, lchnk)
  call outfld (trim(prefix)//'VTGWSPEC', vtgw , ncol, lchnk)

  ! Tau in each direction.
  call outfld (trim(prefix)//'TAUE', taucd(:,:,east), ncol, lchnk)
  call outfld (trim(prefix)//'TAUW', taucd(:,:,west), ncol, lchnk)
  call outfld (trim(prefix)//'TAUN', taucd(:,:,north), ncol, lchnk)
  call outfld (trim(prefix)//'TAUS', taucd(:,:,south), ncol, lchnk)

  call outfld (trim(prefix)//'TAUNET', taucd(:,:,east)+taucd(:,:,west), &
       ncol, lchnk)

contains

  ! Given a value, finds which bin marked by "bounds" the value falls
  ! into.
  elemental function find_bin(val) result(idx)
    real(r8), intent(in) :: val

    integer :: idx

    ! We just have to count how many bounds are exceeded.
    if (val >= 0._r8) then
       idx = count(val > bounds) + 1
    else
       idx = count(val >= bounds) + 1
    end if

  end function find_bin

  ! Convert a speed to a wavenumber between -ngwv and ngwv.
  elemental function c_to_l(c) result(l)
    real(r8), intent(in) :: c

    integer :: l

    l = min( max(int(c/dc),-ngwv), ngwv )

  end function c_to_l

end subroutine gw_spec_outflds

!==========================================================================

! Generates names for tau output across the wave spectrum (e.g.
! BTAUXSn01 or TAUYSp05).
! Probably this should use a wavenumber dimension on one field rather
! than creating a ton of numbered fields.
character(len=9) pure function tau_fld_name(l, prefix, x_not_y)
  ! Wavenumber
  integer, intent(in) :: l
  ! Single-character prefix for output
  character(len=1), intent(in) :: prefix
  ! X or Y?
  logical, intent(in) :: x_not_y

  character(len=2) :: num_str

  tau_fld_name = trim(prefix)

  tau_fld_name = trim(tau_fld_name)//"TAU"

  if (x_not_y) then
     tau_fld_name = trim(tau_fld_name)//"XS"
  else
     tau_fld_name = trim(tau_fld_name)//"YS"
  end if

  if (l < 0) then
     tau_fld_name = trim(tau_fld_name)//"n"
  else
     tau_fld_name = trim(tau_fld_name)//"p"
  end if

  write(num_str,'(I2.2)') abs(l)

  tau_fld_name = trim(tau_fld_name)//num_str

end function tau_fld_name

!==========================================================================

end module gw_drag
