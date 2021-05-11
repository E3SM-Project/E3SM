!===============================================================================
!
! !MODULE: seq_diag_mod -- computes spatial \& time averages of fluxed quatities
!
! !DESCRIPTION:
!    The coupler is required to do certain diagnostics, those calculations are
!    located in this module.
!
! !REMARKS:
!    CESM sign convention for fluxes is positive downward with hierarchy being
!       atm/glc/lnd/rof/ice/ocn
!    Sign convention:
!       positive value <=> the model is gaining water, heat, momentum, etc.
!    Unit convention:
!       heat flux     ~ W/m^2
!       momentum flux ~ N/m^2
!       water flux    ~ (kg/s)/m^2
!       salt  flux    ~ (kg/s)/m^2
!
! !REVISION HISTORY:
!    2012-aug-20 - T. Craig    - add rof component
!    2008-jul-10 - T. Craig    - updated budget implementation
!    2007-may-07 - B. Kauffman - initial port to cpl7.
!    2002-nov-21 - R. Jacob    - initial port to cpl6.
!    199x-mmm-dd - B. Kauffman - original version in cpl4.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_diag_mct
  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8, in=>shr_kind_in
  use shr_kind_mod, only: i8 => shr_kind_i8,  cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod, only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod, only : shr_mpi_max, shr_mpi_sum
  use shr_const_mod, only: shr_const_rearth, shr_const_pi, shr_const_latice, &
       shr_const_ice_ref_sal, shr_const_ocn_ref_sal, shr_const_isspval
  use mct_mod, only: mct_ggrid, mct_avect, mct_avect_lsize, mct_string, &
       mct_string_tochar, mct_gsmap, mct_aVect_indexRA, MCT_AVECT_NRATTR, &
       mct_string_clean, mct_avect_getrlist
  use esmf, only : esmf_clock
  use shr_log_mod, only: s_logunit=>shr_log_unit
  use seq_comm_mct, only: logunit, cplid, seq_comm_setptrs, seq_comm_clean
  use seq_timemgr_mod, only : seq_timemgr_EClockGetData
  use component_type_mod, only : COMPONENT_GET_DOM_CX, COMPONENT_GET_C2X_CX, &
       COMPONENT_GET_X2C_CX, COMPONENT_TYPE
  use seq_infodata_mod, only : seq_infodata_type, seq_infodata_getdata
  use shr_reprosum_mod, only: shr_reprosum_calc

  implicit none
  save
  private

  ! !PUBLIC TYPES:

  ! none

  !PUBLIC MEMBER FUNCTIONS:

  public seq_diag_zero_mct
  public seq_diag_atm_mct
  public seq_diag_lnd_mct
  public seq_diag_rof_mct
  public seq_diag_glc_mct
  public seq_diag_ocn_mct
  public seq_diag_ice_mct
  public seq_diag_accum_mct
  public seq_diag_sum0_mct
  public seq_diag_print_mct
  public seq_diag_avect_mct
  public seq_diag_avloc_mct
  public seq_diag_avdiff_mct

  !EOP

  !----------------------------------------------------------------------------
  ! Local data
  !----------------------------------------------------------------------------

  !----- local constants -----
  real(r8),parameter :: HFLXtoWFLX = & ! water flux implied by latent heat of fusion
       &  - (shr_const_ocn_ref_sal-shr_const_ice_ref_sal) / &
       &    (shr_const_ocn_ref_sal*shr_const_latice)

  real(r8),parameter :: SFLXtoWFLX = & ! water flux implied by salt flux
                                ! WFLX (kg/m^2s) = -SFLX (kg/m^2s)
                                !                  / ocn_ref_sal (psu) (34.7g/kg)
                                !		   / 1.e-3 kg/g
       -1._r8/(shr_const_ocn_ref_sal*1.e-3_r8)


  !--- C for component ---
  !--- "r" is recieve in the coupler, "s" is send from the coupler

  integer(in),parameter :: c_size = 22

  integer(in),parameter :: c_atm_as   = 1 ! model index: atm
  integer(in),parameter :: c_atm_ar   = 2 ! model index: atm
  integer(in),parameter :: c_inh_is   = 3 ! model index: ice, northern
  integer(in),parameter :: c_inh_ir   = 4 ! model index: ice, northern
  integer(in),parameter :: c_ish_is   = 5 ! model index: ice, southern
  integer(in),parameter :: c_ish_ir   = 6 ! model index: ice, southern
  integer(in),parameter :: c_lnd_ls   = 7 ! model index: lnd
  integer(in),parameter :: c_lnd_lr   = 8 ! model index: lnd
  integer(in),parameter :: c_ocn_os   = 9 ! model index: ocn
  integer(in),parameter :: c_ocn_or   =10 ! model index: ocn
  integer(in),parameter :: c_rof_rs   =11 ! model index: rof
  integer(in),parameter :: c_rof_rr   =12 ! model index: rof
  integer(in),parameter :: c_glc_gs   =13 ! model index: glc
  integer(in),parameter :: c_glc_gr   =14 ! model index: glc
  ! --- on atm grid ---
  integer(in),parameter :: c_inh_as   =15 ! model index: ice, northern
  integer(in),parameter :: c_inh_ar   =16 ! model index: ice, northern
  integer(in),parameter :: c_ish_as   =17 ! model index: ice, southern
  integer(in),parameter :: c_ish_ar   =18 ! model index: ice, southern
  integer(in),parameter :: c_lnd_as   =19 ! model index: lnd
  integer(in),parameter :: c_lnd_ar   =20 ! model index: lnd
  integer(in),parameter :: c_ocn_as   =21 ! model index: ocn
  integer(in),parameter :: c_ocn_ar   =22 ! model index: ocn

  character(len=8),parameter :: cname(c_size) = &
       (/' c2a_atm',' a2c_atm',' c2i_inh',' i2c_inh',' c2i_ish',' i2c_ish', &
       ' c2l_lnd',' l2c_lnd',' c2o_ocn',' o2c_ocn',' c2r_rof',' r2c_rof', &
       ' c2g_glc',' g2c_glc', &
       ' c2a_inh',' a2c_inh',' c2a_ish',' a2c_ish', &
       ' c2a_lnd',' a2c_lnd',' c2a_ocn',' a2c_ocn' /)

  !--- F for field ---

  integer(in),parameter :: f_area      = 1     ! area (wrt to unit sphere)
  integer(in),parameter :: f_hfrz      = 2     ! heat : latent, freezing
  integer(in),parameter :: f_hmelt     = 3     ! heat : latent, melting
  integer(in),parameter :: f_hswnet    = 4     ! heat : short wave, net
  integer(in),parameter :: f_hlwdn     = 5     ! heat : longwave down
  integer(in),parameter :: f_hlwup     = 6     ! heat : longwave up
  integer(in),parameter :: f_hlatv     = 7     ! heat : latent, vaporization
  integer(in),parameter :: f_hlatf     = 8     ! heat : latent, fusion, snow
  integer(in),parameter :: f_hioff     = 9     ! heat : latent, fusion, frozen runoff
  integer(in),parameter :: f_hsen      =10     ! heat : sensible
  integer(in),parameter :: f_wfrz      =11     ! water: freezing
  integer(in),parameter :: f_wmelt     =12     ! water: melting
  integer(in),parameter :: f_wrain     =13     ! water: precip, liquid
  integer(in),parameter :: f_wsnow     =14     ! water: precip, frozen
  integer(in),parameter :: f_wevap     =15     ! water: evaporation
  integer(in),parameter :: f_wsalt     =16     ! water: water equivalent of salt flux
  integer(in),parameter :: f_wroff     =17     ! water: runoff/flood
  integer(in),parameter :: f_wioff     =18     ! water: frozen runoff
  integer(in),parameter :: f_wfrz_16O  =19     ! water: freezing
  integer(in),parameter :: f_wmelt_16O =20     ! water: melting
  integer(in),parameter :: f_wrain_16O =21     ! water: precip, liquid
  integer(in),parameter :: f_wsnow_16O =22     ! water: precip, frozen
  integer(in),parameter :: f_wevap_16O =23     ! water: evaporation
  integer(in),parameter :: f_wroff_16O =24     ! water: runoff/flood
  integer(in),parameter :: f_wioff_16O =25     ! water: frozen runoff
  integer(in),parameter :: f_wfrz_18O  =26     ! water: freezing
  integer(in),parameter :: f_wmelt_18O =27     ! water: melting
  integer(in),parameter :: f_wrain_18O =28     ! water: precip, liquid
  integer(in),parameter :: f_wsnow_18O =29     ! water: precip, frozen
  integer(in),parameter :: f_wevap_18O =30     ! water: evaporation
  integer(in),parameter :: f_wroff_18O =31     ! water: runoff/flood
  integer(in),parameter :: f_wioff_18O =32     ! water: frozen runoff
  integer(in),parameter :: f_wfrz_HDO  =33     ! water: freezing
  integer(in),parameter :: f_wmelt_HDO =34     ! water: melting
  integer(in),parameter :: f_wrain_HDO =35     ! water: precip, liquid
  integer(in),parameter :: f_wsnow_HDO =36     ! water: precip, frozen
  integer(in),parameter :: f_wevap_HDO =37     ! water: evaporation
  integer(in),parameter :: f_wroff_HDO =38     ! water: runoff/flood
  integer(in),parameter :: f_wioff_HDO =39     ! water: frozen runoff

  integer(in),parameter :: f_size     = f_wioff_HDO   ! Total array size of all elements
  integer(in),parameter :: f_a        = f_area        ! 1st index for area
  integer(in),parameter :: f_a_end    = f_area        ! last index for area
  integer(in),parameter :: f_h        = f_hfrz        ! 1st index for heat
  integer(in),parameter :: f_h_end    = f_hsen        ! Last index for heat
  integer(in),parameter :: f_w        = f_wfrz        ! 1st index for water
  integer(in),parameter :: f_w_end    = f_wioff       ! Last index for water
  integer(in),parameter :: f_16O      = f_wfrz_16O    ! 1st index for 16O water isotope
  integer(in),parameter :: f_18O      = f_wfrz_18O    ! 1st index for 18O water isotope
  integer(in),parameter :: f_HDO      = f_wfrz_HDO    ! 1st index for HDO water isotope
  integer(in),parameter :: f_16O_end  = f_wioff_16O   ! Last index for 16O water isotope
  integer(in),parameter :: f_18O_end  = f_wioff_18O   ! Last index for 18O water isotope
  integer(in),parameter :: f_HDO_end  = f_wioff_HDO   ! Last index for HDO water isotope

  character(len=12),parameter :: fname(f_size) = &

       (/'        area','     hfreeze','       hmelt','      hnetsw','       hlwdn', &
       '       hlwup','     hlatvap','     hlatfus','      hiroff','        hsen', &
       '     wfreeze','       wmelt','       wrain','       wsnow',                &
       '       wevap','    weqsaltf','     wrunoff','     wfrzrof',                &
       ' wfreeze_16O','   wmelt_16O','   wrain_16O','   wsnow_16O',                &
       '   wevap_16O',' wrunoff_16O',' wfrzrof_16O',                               &
       ' wfreeze_18O','   wmelt_18O','   wrain_18O','   wsnow_18O',                &
       '   wevap_18O',' wrunoff_18O',' wfrzrof_18O',                               &
       ' wfreeze_HDO','   wmelt_HDO','   wrain_HDO','   wsnow_HDO',                &
       '   wevap_HDO',' wrunoff_HDO',' wfrzrof_HDO'/)

  !--- P for period ---

  integer(in),parameter :: p_size = 5

  integer(in),parameter :: p_inst = 1
  integer(in),parameter :: p_day  = 2
  integer(in),parameter :: p_mon  = 3
  integer(in),parameter :: p_ann  = 4
  integer(in),parameter :: p_inf  = 5

  character(len=8),parameter :: pname(p_size) = &
       (/'    inst','   daily',' monthly','  annual','all_time' /)

  logical          :: flds_wiso             ! If water isotope fields are active

  ! !PUBLIC DATA MEMBERS

  !--- time-averaged (annual?) global budge diagnostics ---
  !--- note: call sum0 then save budg_dataG and budg_ns on restart from/to root pe ---
  real(r8),public :: budg_dataL(f_size,c_size,p_size) ! local sum, valid on all pes
  real(r8),public :: budg_dataG(f_size,c_size,p_size) ! global sum, valid only on root pe
  real(r8),public :: budg_ns   (f_size,c_size,p_size) ! counter, valid only on root pe

  character(len=*),parameter :: afldname  = 'aream'
  character(len=*),parameter :: latname   = 'lat'
  character(len=*),parameter :: afracname = 'afrac'
  character(len=*),parameter :: lfracname = 'lfrac'
  character(len=*),parameter :: ofracname = 'ofrac'
  character(len=*),parameter :: ifracname = 'ifrac'

  character(*),parameter :: modName = "(seq_diag_mct) "

  integer(in),parameter :: debug = 0 ! internal debug level

  ! !PRIVATE DATA MEMBERS

  integer :: index_a2x_Faxa_swnet
  integer :: index_a2x_Faxa_lwdn
  integer :: index_a2x_Faxa_rainc
  integer :: index_a2x_Faxa_rainl
  integer :: index_a2x_Faxa_snowc
  integer :: index_a2x_Faxa_snowl

  integer :: index_x2a_Faxx_lwup
  integer :: index_x2a_Faxx_lat
  integer :: index_x2a_Faxx_sen
  integer :: index_x2a_Faxx_evap

  integer :: index_l2x_Fall_swnet
  integer :: index_l2x_Fall_lwup
  integer :: index_l2x_Fall_lat
  integer :: index_l2x_Fall_sen
  integer :: index_l2x_Fall_evap
  integer :: index_l2x_Flrl_rofsur
  integer :: index_l2x_Flrl_rofgwl
  integer :: index_l2x_Flrl_rofsub
  integer :: index_l2x_Flrl_rofdto
  integer :: index_l2x_Flrl_rofi
  integer :: index_l2x_Flrl_irrig

  integer :: index_x2l_Faxa_lwdn
  integer :: index_x2l_Faxa_rainc
  integer :: index_x2l_Faxa_rainl
  integer :: index_x2l_Faxa_snowc
  integer :: index_x2l_Faxa_snowl
  integer :: index_x2l_Flrr_flood

  integer :: index_r2x_Forr_rofl
  integer :: index_r2x_Forr_rofi
  integer :: index_r2x_Firr_rofi
  integer :: index_r2x_Flrr_flood

  integer :: index_x2r_Flrl_rofsur
  integer :: index_x2r_Flrl_rofgwl
  integer :: index_x2r_Flrl_rofsub
  integer :: index_x2r_Flrl_rofdto
  integer :: index_x2r_Flrl_rofi
  integer :: index_x2r_Flrl_irrig

  integer :: index_o2x_Fioo_frazil ! currently used by e3sm
  integer :: index_o2x_Fioo_q      ! currently used by cesm

  integer :: index_xao_Faox_lwup
  integer :: index_xao_Faox_lat
  integer :: index_xao_Faox_sen
  integer :: index_xao_Faox_evap

  integer :: index_x2o_Foxx_lwup
  integer :: index_x2o_Foxx_lat
  integer :: index_x2o_Foxx_sen
  integer :: index_x2o_Foxx_evap
  integer :: index_x2o_Foxx_swnet
  integer :: index_x2o_Foxx_rofl
  integer :: index_x2o_Foxx_rofi
  integer :: index_x2o_Faxa_lwdn
  integer :: index_x2o_Faxa_rain
  integer :: index_x2o_Faxa_snow
  integer :: index_x2o_Fioi_melth
  integer :: index_x2o_Fioi_meltw
  integer :: index_x2o_Fioi_bergh
  integer :: index_x2o_Fioi_bergw
  integer :: index_x2o_Fioi_salt

  integer :: index_i2x_Fioi_melth
  integer :: index_i2x_Fioi_meltw
  integer :: index_i2x_Fioi_salt
  integer :: index_i2x_Faii_swnet
  integer :: index_i2x_Fioi_swpen
  integer :: index_i2x_Faii_lwup
  integer :: index_i2x_Faii_lat
  integer :: index_i2x_Faii_sen
  integer :: index_i2x_Faii_evap

  integer :: index_x2i_Faxa_lwdn
  integer :: index_x2i_Faxa_rain
  integer :: index_x2i_Faxa_snow
  integer :: index_x2i_Fioo_frazil !currently used by e3sm
  integer :: index_x2i_Fioo_q      !currently used by cesm
  integer :: index_x2i_Fixx_rofi

  integer :: index_g2x_Fogg_rofl
  integer :: index_g2x_Fogg_rofi
  integer :: index_g2x_Figg_rofi

  integer :: index_x2o_Foxx_rofl_16O
  integer :: index_x2o_Foxx_rofi_16O
  integer :: index_x2o_Foxx_rofl_18O
  integer :: index_x2o_Foxx_rofi_18O
  integer :: index_x2o_Foxx_rofl_HDO
  integer :: index_x2o_Foxx_rofi_HDO

  integer :: index_a2x_Faxa_rainc_16O
  integer :: index_a2x_Faxa_rainc_18O
  integer :: index_a2x_Faxa_rainc_HDO
  integer :: index_a2x_Faxa_rainl_16O
  integer :: index_a2x_Faxa_rainl_18O
  integer :: index_a2x_Faxa_rainl_HDO
  integer :: index_a2x_Faxa_snowc_16O
  integer :: index_a2x_Faxa_snowc_18O
  integer :: index_a2x_Faxa_snowc_HDO
  integer :: index_a2x_Faxa_snowl_16O
  integer :: index_a2x_Faxa_snowl_18O
  integer :: index_a2x_Faxa_snowl_HDO

  integer :: index_x2a_Faxx_evap_16O
  integer :: index_x2a_Faxx_evap_18O
  integer :: index_x2a_Faxx_evap_HDO

  integer :: index_l2x_Fall_evap_16O
  integer :: index_l2x_Fall_evap_18O
  integer :: index_l2x_Fall_evap_HDO

  integer :: index_l2x_Flrl_rofl_16O
  integer :: index_l2x_Flrl_rofl_18O
  integer :: index_l2x_Flrl_rofl_HDO
  integer :: index_l2x_Flrl_rofi_16O
  integer :: index_l2x_Flrl_rofi_18O
  integer :: index_l2x_Flrl_rofi_HDO

  integer :: index_x2l_Faxa_rainc_16O
  integer :: index_x2l_Faxa_rainc_18O
  integer :: index_x2l_Faxa_rainc_HDO
  integer :: index_x2l_Faxa_rainl_16O
  integer :: index_x2l_Faxa_rainl_18O
  integer :: index_x2l_Faxa_rainl_HDO
  integer :: index_x2l_Faxa_snowc_16O
  integer :: index_x2l_Faxa_snowc_18O
  integer :: index_x2l_Faxa_snowc_HDO
  integer :: index_x2l_Faxa_snowl_16O
  integer :: index_x2l_Faxa_snowl_18O
  integer :: index_x2l_Faxa_snowl_HDO
  integer :: index_x2l_Flrr_flood_16O
  integer :: index_x2l_Flrr_flood_18O
  integer :: index_x2l_Flrr_flood_HDO

  integer :: index_r2x_Forr_rofl_16O
  integer :: index_r2x_Forr_rofl_18O
  integer :: index_r2x_Forr_rofl_HDO
  integer :: index_r2x_Forr_rofi_16O
  integer :: index_r2x_Forr_rofi_18O
  integer :: index_r2x_Forr_rofi_HDO
  integer :: index_r2x_Flrr_flood_16O
  integer :: index_r2x_Flrr_flood_18O
  integer :: index_r2x_Flrr_flood_HDO

  integer :: index_x2r_Flrl_rofl_16O
  integer :: index_x2r_Flrl_rofl_18O
  integer :: index_x2r_Flrl_rofl_HDO
  integer :: index_x2r_Flrl_rofi_16O
  integer :: index_x2r_Flrl_rofi_18O
  integer :: index_x2r_Flrl_rofi_HDO

  integer :: index_xao_Faox_evap_16O
  integer :: index_xao_Faox_evap_18O
  integer :: index_xao_Faox_evap_HDO

  integer :: index_x2o_Fioi_meltw_16O
  integer :: index_x2o_Fioi_meltw_18O
  integer :: index_x2o_Fioi_meltw_HDO
  integer :: index_x2o_Faxa_rain_16O
  integer :: index_x2o_Faxa_rain_18O
  integer :: index_x2o_Faxa_rain_HDO
  integer :: index_x2o_Faxa_snow_16O
  integer :: index_x2o_Faxa_snow_18O
  integer :: index_x2o_Faxa_snow_HDO

  integer :: index_i2x_Fioi_meltw_16O
  integer :: index_i2x_Fioi_meltw_18O
  integer :: index_i2x_Fioi_meltw_HDO
  integer :: index_i2x_Faii_evap_16O
  integer :: index_i2x_Faii_evap_18O
  integer :: index_i2x_Faii_evap_HDO

  integer :: index_x2i_Faxa_rain_16O
  integer :: index_x2i_Faxa_rain_18O
  integer :: index_x2i_Faxa_rain_HDO
  integer :: index_x2i_Faxa_snow_16O
  integer :: index_x2i_Faxa_snow_18O
  integer :: index_x2i_Faxa_snow_HDO

  !===============================================================================
contains
  !===============================================================================

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_zero_mct - zero out global budget diagnostic data.
  !
  ! !DESCRIPTION:
  !    Zero out global budget diagnostic data.
  !
  ! !REVISION HISTORY:
  !    2008-jul-11 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_zero_mct(EClock,mode)

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(in),optional :: EClock
    character(len=*), intent(in),optional :: mode

    !EOP

    integer(IN) :: ip,yr,mon,day,sec
    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_zero_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (.not. present(EClock) .and. .not. present(mode)) then
       call shr_sys_abort(subName//' ERROR EClock or mode should be present')
    endif

    if (present(EClock)) then
       call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
            curr_mon=mon,curr_day=day,curr_tod=sec)

       do ip = 1,p_size
          if (ip == p_inst) then
             budg_dataL(:,:,ip) = 0.0_r8
             budg_dataG(:,:,ip) = 0.0_r8
             budg_ns(:,:,ip) = 0.0_r8
          endif
          if (ip==p_day .and. sec==0) then
             budg_dataL(:,:,ip) = 0.0_r8
             budg_dataG(:,:,ip) = 0.0_r8
             budg_ns(:,:,ip) = 0.0_r8
          endif
          if (ip==p_mon .and. day==1 .and. sec==0) then
             budg_dataL(:,:,ip) = 0.0_r8
             budg_dataG(:,:,ip) = 0.0_r8
             budg_ns(:,:,ip) = 0.0_r8
          endif
          if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
             budg_dataL(:,:,ip) = 0.0_r8
             budg_dataG(:,:,ip) = 0.0_r8
             budg_ns(:,:,ip) = 0.0_r8
          endif
       enddo
    endif

    if (present(mode)) then
       if (trim(mode) == 'inst') then
          budg_dataL(:,:,p_inst) = 0.0_r8
          budg_dataG(:,:,p_inst) = 0.0_r8
          budg_ns(:,:,p_inst) = 0.0_r8
       elseif (trim(mode) == 'day') then
          budg_dataL(:,:,p_day) = 0.0_r8
          budg_dataG(:,:,p_day) = 0.0_r8
          budg_ns(:,:,p_day) = 0.0_r8
       elseif (trim(mode) == 'mon') then
          budg_dataL(:,:,p_mon) = 0.0_r8
          budg_dataG(:,:,p_mon) = 0.0_r8
          budg_ns(:,:,p_mon) = 0.0_r8
       elseif (trim(mode) == 'ann') then
          budg_dataL(:,:,p_ann) = 0.0_r8
          budg_dataG(:,:,p_ann) = 0.0_r8
          budg_ns(:,:,p_ann) = 0.0_r8
       elseif (trim(mode) == 'inf') then
          budg_dataL(:,:,p_inf) = 0.0_r8
          budg_dataG(:,:,p_inf) = 0.0_r8
          budg_ns(:,:,p_inf) = 0.0_r8
       elseif (trim(mode) == 'all') then
          budg_dataL(:,:,:) = 0.0_r8
          budg_dataG(:,:,:) = 0.0_r8
          budg_ns(:,:,:) = 0.0_r8
       else
          call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
       endif
    endif

  end subroutine seq_diag_zero_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_accum_mct - accum out global budget diagnostic data.
  !
  ! !DESCRIPTION:
  !    Accum out global budget diagnostic data.
  !
  ! !REVISION HISTORY:
  !    2008-jul-11 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_accum_mct()

    ! !INPUT/OUTPUT PARAMETERS:

    !EOP

    integer(in) :: ip

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_accum_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    do ip = p_inst+1,p_size
       budg_dataL(:,:,ip) = budg_dataL(:,:,ip) + budg_dataL(:,:,p_inst)
    enddo
    budg_ns(:,:,:) = budg_ns(:,:,:) + 1.0_r8

  end subroutine seq_diag_accum_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_sum0_mct - sum local to global on root
  !
  ! !DESCRIPTION:
  !    Sum local values to global on root
  !
  ! !REVISION HISTORY:
  !    2008-jul-19 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_sum0_mct()

    ! !INPUT/OUTPUT PARAMETERS:

    !EOP

    real(r8) :: budg_dataGtmp(f_size,c_size,p_size) ! temporary sum
    integer(in)      :: mpicom      ! mpi comm
    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_sum0_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,mpicom=mpicom)
    budg_dataGtmp = 0.0_r8
    call shr_mpi_sum(budg_dataL,budg_dataGtmp,mpicom,subName)
    budg_dataG = budg_dataG + budg_dataGtmp
    budg_dataL = 0.0_r8

  end subroutine seq_diag_sum0_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_atm_mct - compute global atm input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global atm input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2008-jul-10 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_atm_mct( atm, frac_a, infodata, do_a2x, do_x2a)

    ! !INPUT/OUTPUT PARAMETERS:

    type(component_type)    , intent(in)           :: atm    ! component type for instance1
    type(mct_aVect)         , intent(in)           :: frac_a ! frac bundle
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_a2x
    logical                 , intent(in), optional :: do_x2a

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: a2x_a        ! model to drv bundle
    type(mct_aVect), pointer :: x2a_a        ! drv to model bundle
    type(mct_ggrid), pointer :: dom_a
    integer(in)              :: k,n,ic,nf,ip      ! generic index
    integer(in)              :: kArea             ! index of area field in aVect
    integer(in)              :: kLat              ! index of lat field in aVect
    integer(in)              :: kl,ka,ko,ki       ! fraction indices
    integer(in)              :: lSize             ! size of aVect
    real(r8)                 :: ca_a       ! area of a grid cell
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_atm = .false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_atm_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    dom_a => component_get_dom_cx(atm)
    a2x_a => component_get_c2x_cx(atm)
    x2a_a => component_get_x2c_cx(atm)

    kArea = mct_aVect_indexRA(dom_a%data,afldname)
    kLat  = mct_aVect_indexRA(dom_a%data,latname)
    ka    = mct_aVect_indexRA(frac_a,afracname)
    kl    = mct_aVect_indexRA(frac_a,lfracname)
    ko    = mct_aVect_indexRA(frac_a,ofracname)
    ki    = mct_aVect_indexRA(frac_a,ifracname)

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    ip = p_inst

    if (present(do_a2x)) then
       if (first_time) then
          index_a2x_Faxa_swnet  = mct_aVect_indexRA(a2x_a,'Faxa_swnet')
          index_a2x_Faxa_lwdn   = mct_aVect_indexRA(a2x_a,'Faxa_lwdn')
          index_a2x_Faxa_rainc  = mct_aVect_indexRA(a2x_a,'Faxa_rainc')
          index_a2x_Faxa_rainl  = mct_aVect_indexRA(a2x_a,'Faxa_rainl')
          index_a2x_Faxa_snowc  = mct_aVect_indexRA(a2x_a,'Faxa_snowc')
          index_a2x_Faxa_snowl  = mct_aVect_indexRA(a2x_a,'Faxa_snowl')

          index_a2x_Faxa_rainc_16O   = mct_aVect_indexRA(a2x_a,'Faxa_rainc_16O',perrWith='quiet')
          if ( index_a2x_Faxa_rainc_16O /= 0 ) flds_wiso_atm = .true.
          if ( flds_wiso_atm )then
             flds_wiso = .true.
             index_a2x_Faxa_rainc_18O   = mct_aVect_indexRA(a2x_a,'Faxa_rainc_18O')
             index_a2x_Faxa_rainc_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_rainc_HDO')
             index_a2x_Faxa_rainl_16O   = mct_aVect_indexRA(a2x_a,'Faxa_rainl_16O')
             index_a2x_Faxa_rainl_18O   = mct_aVect_indexRA(a2x_a,'Faxa_rainl_18O')
             index_a2x_Faxa_rainl_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_rainl_HDO')
             index_a2x_Faxa_snowc_16O   = mct_aVect_indexRA(a2x_a,'Faxa_snowc_16O')
             index_a2x_Faxa_snowc_18O   = mct_aVect_indexRA(a2x_a,'Faxa_snowc_18O')
             index_a2x_Faxa_snowc_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_snowc_HDO')
             index_a2x_Faxa_snowl_16O   = mct_aVect_indexRA(a2x_a,'Faxa_snowl_16O')
             index_a2x_Faxa_snowl_18O   = mct_aVect_indexRA(a2x_a,'Faxa_snowl_18O')
             index_a2x_Faxa_snowl_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_snowl_HDO')
          end if

       end if

       lSize = mct_avect_lSize(a2x_a)
       do n=1,lSize
          do k=1,4

             if (k == 1) then
                ic = c_atm_ar
                ca_a = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
             elseif (k == 2) then
                ic = c_lnd_ar
                ca_a =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
             elseif (k == 3) then
                ic = c_ocn_ar
                ca_a =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
             elseif (k == 4) then
                if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
                   ic = c_inh_ar
                else
                   ic = c_ish_ar
                endif
                ca_a = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
             endif

             nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a
             nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*a2x_a%rAttr(index_a2x_Faxa_swnet,n)
             nf = f_hlwdn ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*a2x_a%rAttr(index_a2x_Faxa_lwdn,n)
             nf = f_wrain ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*a2x_a%rAttr(index_a2x_Faxa_rainc,n) &
                  + ca_a*a2x_a%rAttr(index_a2x_Faxa_rainl,n)
             nf = f_wsnow ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*a2x_a%rAttr(index_a2x_Faxa_snowc,n) &
                  + ca_a*a2x_a%rAttr(index_a2x_Faxa_snowl,n)
             if ( flds_wiso_atm )then
                nf = f_wrain_16O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_rainc_16O,n) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_rainl_16O,n)
                nf = f_wrain_18O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_rainc_18O,n) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_rainl_18O,n)
                nf = f_wrain_HDO;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_rainc_HDO,n) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_rainl_HDO,n)
                nf = f_wsnow_16O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_snowc_16O,n) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_snowl_16O,n)
                nf = f_wsnow_18O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_snowc_18O,n) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_snowl_18O,n)
                nf = f_wsnow_HDO;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_snowc_HDO,n) + &
                     ca_a*a2x_a%rAttr(index_a2x_Faxa_snowl_HDO,n)
             end if
          enddo
       enddo
       ! --- heat implied by snow flux ---
       ic = c_atm_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_lnd_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_ocn_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_inh_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_ish_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
    end if

    if (present(do_x2a)) then
       if (first_time) then
          index_x2a_Faxx_lwup   = mct_aVect_indexRA(x2a_a,'Faxx_lwup')
          index_x2a_Faxx_lat    = mct_aVect_indexRA(x2a_a,'Faxx_lat')
          index_x2a_Faxx_sen    = mct_aVect_indexRA(x2a_a,'Faxx_sen')
          index_x2a_Faxx_evap   = mct_aVect_indexRA(x2a_a,'Faxx_evap')

          if ( flds_wiso_atm )then
             index_x2a_Faxx_evap_16O = mct_aVect_indexRA(x2a_a,'Faxx_evap_16O')
             index_x2a_Faxx_evap_18O = mct_aVect_indexRA(x2a_a,'Faxx_evap_18O')
             index_x2a_Faxx_evap_HDO = mct_aVect_indexRA(x2a_a,'Faxx_evap_HDO')
          end if
       end if

       lSize = mct_avect_lSize(x2a_a)
       do n=1,lSize
          do k=1,4

             if (k == 1) then
                ic = c_atm_as
                ca_a = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
             elseif (k == 2) then
                ic = c_lnd_as
                ca_a =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
             elseif (k == 3) then
                ic = c_ocn_as
                ca_a =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
             elseif (k == 4) then
                if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
                   ic = c_inh_as
                else
                   ic = c_ish_as
                endif
                ca_a = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
             endif

             nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a
             nf = f_hlwup; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*x2a_a%rAttr(index_x2a_Faxx_lwup,n)
             nf = f_hlatv; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*x2a_a%rAttr(index_x2a_Faxx_lat,n)
             nf = f_hsen ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*x2a_a%rAttr(index_x2a_Faxx_sen,n)
             nf = f_wevap; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*x2a_a%rAttr(index_x2a_Faxx_evap,n)

             if ( flds_wiso_atm )then
                nf = f_wevap_16O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*x2a_a%rAttr(index_x2a_Faxx_evap_16O,n)
                nf = f_wevap_18O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*x2a_a%rAttr(index_x2a_Faxx_evap_18O,n)
                nf = f_wevap_HDO;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*x2a_a%rAttr(index_x2a_Faxx_evap_HDO,n)
             end if

          enddo
       enddo
    end if

    first_time = .false.

  end subroutine seq_diag_atm_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_lnd_mct - compute global lnd input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global lnd input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2008-jul-10 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_lnd_mct( lnd, frac_l, infodata, do_l2x, do_x2l)

    type(component_type)    , intent(in)           :: lnd    ! component type for instance1
    type(mct_aVect)         , intent(in)           :: frac_l ! frac bundle
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_l2x
    logical                 , intent(in), optional :: do_x2l

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: l2x_l        ! model to drv bundle
    type(mct_aVect), pointer :: x2l_l        ! drv to model bundle
    type(mct_ggrid), pointer :: dom_l
    integer(in)              :: n,ic,nf,ip ! generic index
    integer(in)              :: kArea        ! index of area field in aVect
    integer(in)              :: kl           ! fraction indices
    integer(in)              :: lSize        ! size of aVect
    real(r8)                 :: ca_l         ! area of a grid cell
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_lnd = .false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_lnd_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    dom_l => component_get_dom_cx(lnd)
    l2x_l => component_get_c2x_cx(lnd)
    x2l_l => component_get_x2c_cx(lnd)

    ip = p_inst

    kArea = mct_aVect_indexRA(dom_l%data,afldname)
    kl    = mct_aVect_indexRA(frac_l,lfracname)

    if (present(do_l2x)) then
       if (first_time) then
          index_l2x_Fall_swnet  = mct_aVect_indexRA(l2x_l,'Fall_swnet')
          index_l2x_Fall_lwup   = mct_aVect_indexRA(l2x_l,'Fall_lwup')
          index_l2x_Fall_lat    = mct_aVect_indexRA(l2x_l,'Fall_lat')
          index_l2x_Fall_sen    = mct_aVect_indexRA(l2x_l,'Fall_sen')
          index_l2x_Fall_evap   = mct_aVect_indexRA(l2x_l,'Fall_evap')
          index_l2x_Flrl_rofsur = mct_aVect_indexRA(l2x_l,'Flrl_rofsur')
          index_l2x_Flrl_rofgwl = mct_aVect_indexRA(l2x_l,'Flrl_rofgwl')
          index_l2x_Flrl_rofsub = mct_aVect_indexRA(l2x_l,'Flrl_rofsub')
          index_l2x_Flrl_rofdto = mct_aVect_indexRA(l2x_l,'Flrl_rofdto')
          index_l2x_Flrl_rofi   = mct_aVect_indexRA(l2x_l,'Flrl_rofi')
          index_l2x_Flrl_irrig  = mct_aVect_indexRA(l2x_l,'Flrl_irrig', perrWith='quiet')

          index_l2x_Fall_evap_16O    = mct_aVect_indexRA(l2x_l,'Fall_evap_16O',perrWith='quiet')
          if ( index_l2x_Fall_evap_16O /= 0 ) flds_wiso_lnd = .true.
          if ( flds_wiso_lnd )then
             flds_wiso = .true.
             index_l2x_Fall_evap_18O    = mct_aVect_indexRA(l2x_l,'Fall_evap_18O')
             index_l2x_Fall_evap_HDO    = mct_aVect_indexRA(l2x_l,'Fall_evap_HDO')
             index_l2x_Flrl_rofl_16O  = mct_aVect_indexRA(l2x_l,'Flrl_rofl_16O')
             index_l2x_Flrl_rofl_18O  = mct_aVect_indexRA(l2x_l,'Flrl_rofl_18O')
             index_l2x_Flrl_rofl_HDO  = mct_aVect_indexRA(l2x_l,'Flrl_rofl_HDO')
             index_l2x_Flrl_rofi_16O  = mct_aVect_indexRA(l2x_l,'Flrl_rofi_16O')
             index_l2x_Flrl_rofi_18O  = mct_aVect_indexRA(l2x_l,'Flrl_rofi_18O')
             index_l2x_Flrl_rofi_HDO  = mct_aVect_indexRA(l2x_l,'Flrl_rofi_HDO')
          end if
       end if

       lSize = mct_avect_lSize(l2x_l)
       ic = c_lnd_lr
       do n=1,lSize
          ca_l =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*l2x_l%rAttr(index_l2x_Fall_swnet,n)
          nf = f_hlwup ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*l2x_l%rAttr(index_l2x_Fall_lwup,n)
          nf = f_hlatv ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*l2x_l%rAttr(index_l2x_Fall_lat,n)
          nf = f_hsen  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*l2x_l%rAttr(index_l2x_Fall_sen,n)
          nf = f_wevap ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*l2x_l%rAttr(index_l2x_Fall_evap,n)
          nf = f_wroff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*l2x_l%rAttr(index_l2x_Flrl_rofsur,n) &
               - ca_l*l2x_l%rAttr(index_l2x_Flrl_rofgwl,n) &
               - ca_l*l2x_l%rAttr(index_l2x_Flrl_rofsub,n) &
               - ca_l*l2x_l%rAttr(index_l2x_Flrl_rofdto,n)
          if (index_l2x_Flrl_irrig /= 0) then
             nf = f_wroff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*l2x_l%rAttr(index_l2x_Flrl_irrig,n)
          end if
          nf = f_wioff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*l2x_l%rAttr(index_l2x_Flrl_rofi,n)

          if ( flds_wiso_lnd )then
             nf = f_wevap_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*l2x_l%rAttr(index_l2x_Fall_evap_16O,n)
             nf = f_wevap_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*l2x_l%rAttr(index_l2x_Fall_evap_18O,n)
             nf = f_wevap_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*l2x_l%rAttr(index_l2x_Fall_evap_HDO,n)

             nf = f_wroff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*l2x_l%rAttr(index_l2x_Flrl_rofl_16O,n)
             nf = f_wroff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*l2x_l%rAttr(index_l2x_Flrl_rofl_18O,n)
             nf = f_wroff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*l2x_l%rAttr(index_l2x_Flrl_rofl_HDO,n)

             nf = f_wioff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*l2x_l%rAttr(index_l2x_Flrl_rofi_16O,n)
             nf = f_wioff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*l2x_l%rAttr(index_l2x_Flrl_rofi_18O,n)
             nf = f_wioff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*l2x_l%rAttr(index_l2x_Flrl_rofi_HDO,n)
          end if
       end do
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
    end if

    if (present(do_x2l)) then
       if (first_time) then
          index_x2l_Faxa_lwdn   = mct_aVect_indexRA(x2l_l,'Faxa_lwdn')
          index_x2l_Faxa_rainc  = mct_aVect_indexRA(x2l_l,'Faxa_rainc')
          index_x2l_Faxa_rainl  = mct_aVect_indexRA(x2l_l,'Faxa_rainl')
          index_x2l_Faxa_snowc  = mct_aVect_indexRA(x2l_l,'Faxa_snowc')
          index_x2l_Faxa_snowl  = mct_aVect_indexRA(x2l_l,'Faxa_snowl')
          index_x2l_Flrr_flood  = mct_aVect_indexRA(x2l_l,'Flrr_flood')

          if ( flds_wiso_lnd )then
             index_x2l_Faxa_rainc_16O = mct_aVect_indexRA(x2l_l,'Faxa_rainc_16O')
             index_x2l_Faxa_rainc_18O = mct_aVect_indexRA(x2l_l,'Faxa_rainc_18O')
             index_x2l_Faxa_rainc_HDO = mct_aVect_indexRA(x2l_l,'Faxa_rainc_HDO')
             index_x2l_Faxa_rainl_16O = mct_aVect_indexRA(x2l_l,'Faxa_rainl_16O')
             index_x2l_Faxa_rainl_18O = mct_aVect_indexRA(x2l_l,'Faxa_rainl_18O')
             index_x2l_Faxa_rainl_HDO = mct_aVect_indexRA(x2l_l,'Faxa_rainl_HDO')
             index_x2l_Faxa_snowc_16O = mct_aVect_indexRA(x2l_l,'Faxa_snowc_16O')
             index_x2l_Faxa_snowc_18O = mct_aVect_indexRA(x2l_l,'Faxa_snowc_18O')
             index_x2l_Faxa_snowc_HDO = mct_aVect_indexRA(x2l_l,'Faxa_snowc_HDO')
             index_x2l_Faxa_snowl_16O = mct_aVect_indexRA(x2l_l,'Faxa_snowl_16O')
             index_x2l_Faxa_snowl_18O = mct_aVect_indexRA(x2l_l,'Faxa_snowl_18O')
             index_x2l_Faxa_snowl_HDO = mct_aVect_indexRA(x2l_l,'Faxa_snowl_HDO')
             index_x2l_Flrr_flood_16O = mct_aVect_indexRA(x2l_l,'Flrr_flood_16O')
             index_x2l_Flrr_flood_18O = mct_aVect_indexRA(x2l_l,'Flrr_flood_18O')
             index_x2l_Flrr_flood_HDO = mct_aVect_indexRA(x2l_l,'Flrr_flood_HDO')
          end if
       end if

       lSize = mct_avect_lSize(x2l_l)
       ic = c_lnd_ls
       do n=1,lSize
          ca_l =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
          nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          nf = f_hlwdn; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*x2l_l%rAttr(index_x2l_Faxa_lwdn,n)
          nf = f_wrain; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*x2l_l%rAttr(index_x2l_Faxa_rainc,n) &
               + ca_l*x2l_l%rAttr(index_x2l_Faxa_rainl,n)
          nf = f_wsnow; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*x2l_l%rAttr(index_x2l_Faxa_snowc,n) &
               + ca_l*x2l_l%rAttr(index_x2l_Faxa_snowl,n)
          nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*x2l_l%rAttr(index_x2l_Flrr_flood,n)

          if ( flds_wiso_lnd )then
             nf = f_wrain_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_rainc_16O,n) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_rainl_16O,n)
             nf = f_wrain_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_rainc_18O,n) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_rainl_18O,n)
             nf = f_wrain_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_rainc_HDO,n) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_rainl_HDO,n)

             nf = f_wsnow_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_snowc_16O,n) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_snowl_16O,n)
             nf = f_wsnow_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_snowc_18O,n) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_snowl_18O,n)
             nf = f_wsnow_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_snowc_HDO,n) + &
                  ca_l*x2l_l%rAttr(index_x2l_Faxa_snowl_HDO,n)

             nf = f_wroff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*x2l_l%rAttr(index_x2l_Flrr_flood_16O,n)
             nf = f_wroff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*x2l_l%rAttr(index_x2l_Flrr_flood_18O,n)
             nf = f_wroff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_l*x2l_l%rAttr(index_x2l_Flrr_flood_HDO,n)
          end if
       end do
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
    end if

    first_time = .false.

  end subroutine seq_diag_lnd_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_rof_mct - compute global rof input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global rof input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2008-jul-10 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_rof_mct( rof, frac_r, infodata)

    type(component_type)    , intent(in) :: rof    ! component type for instance1
    type(mct_aVect)         , intent(in) :: frac_r ! frac bundle
    type(seq_infodata_type) , intent(in) :: infodata

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: r2x_r
    type(mct_aVect), pointer :: x2r_r
    type(mct_ggrid), pointer :: dom_r
    integer(in)              :: n,ic,nf,ip      ! generic index
    integer(in)              :: kArea             ! index of area field in aVect
    integer(in)              :: lSize             ! size of aVect
    real(r8)                 :: ca_r    ! area of a grid cell
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_rof = .false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_rof_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    dom_r => component_get_dom_cx(rof)
    r2x_r => component_get_c2x_cx(rof)
    x2r_r => component_get_x2c_cx(rof)

    if (first_time) then
       index_x2r_Flrl_rofsur = mct_aVect_indexRA(x2r_r,'Flrl_rofsur')
       index_x2r_Flrl_rofgwl = mct_aVect_indexRA(x2r_r,'Flrl_rofgwl')
       index_x2r_Flrl_rofsub = mct_aVect_indexRA(x2r_r,'Flrl_rofsub')
       index_x2r_Flrl_rofdto = mct_aVect_indexRA(x2r_r,'Flrl_rofdto')
       index_x2r_Flrl_irrig  = mct_aVect_indexRA(x2r_r,'Flrl_irrig', perrWith='quiet')
       index_x2r_Flrl_rofi   = mct_aVect_indexRA(x2r_r,'Flrl_rofi')

       index_x2r_Flrl_rofl_16O = mct_aVect_indexRA(x2r_r,'Flrl_rofl_16O', perrWith='quiet')
       if ( index_x2r_Flrl_rofl_16O /= 0 ) flds_wiso_rof = .true.
       if ( flds_wiso_rof )then
          flds_wiso = .true.
          index_x2r_Flrl_rofl_18O = mct_aVect_indexRA(x2r_r,'Flrl_rofl_18O')
          index_x2r_Flrl_rofl_HDO = mct_aVect_indexRA(x2r_r,'Flrl_rofl_HDO')
          index_x2r_Flrl_rofi_16O = mct_aVect_indexRA(x2r_r,'Flrl_rofi_16O')
          index_x2r_Flrl_rofi_18O = mct_aVect_indexRA(x2r_r,'Flrl_rofi_18O')
          index_x2r_Flrl_rofi_HDO = mct_aVect_indexRA(x2r_r,'Flrl_rofi_HDO')
       end if
    end if

    ip = p_inst
    ic = c_rof_rr
    kArea = mct_aVect_indexRA(dom_r%data,afldname)
    lSize = mct_avect_lSize(x2r_r)
    do n=1,lSize
       ca_r =  dom_r%data%rAttr(kArea,n)
       nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*x2r_r%rAttr(index_x2r_Flrl_rofsur,n) &
            + ca_r*x2r_r%rAttr(index_x2r_Flrl_rofgwl,n) &
            + ca_r*x2r_r%rAttr(index_x2r_Flrl_rofsub,n) &
            + ca_r*x2r_r%rAttr(index_x2r_Flrl_rofdto,n)
       if (index_x2r_Flrl_irrig /= 0) then
          nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*x2r_r%rAttr(index_x2r_Flrl_irrig,n)
       end if

       nf = f_wioff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*x2r_r%rAttr(index_x2r_Flrl_rofi,n)

       if ( flds_wiso_rof )then
          nf = f_wroff_16O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*x2r_r%rAttr(index_x2r_Flrl_rofl_16O,n)
          nf = f_wroff_18O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*x2r_r%rAttr(index_x2r_Flrl_rofl_18O,n)
          nf = f_wroff_HDO;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*x2r_r%rAttr(index_x2r_Flrl_rofl_HDO,n)

          nf = f_wioff_16O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*x2r_r%rAttr(index_x2r_Flrl_rofi_16O,n)
          nf = f_wioff_18O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*x2r_r%rAttr(index_x2r_Flrl_rofi_18O,n)
          nf = f_wioff_HDO;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*x2r_r%rAttr(index_x2r_Flrl_rofi_HDO,n)
       end if
    end do
    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

    if (first_time) then
       index_r2x_Forr_rofl   = mct_aVect_indexRA(r2x_r,'Forr_rofl')
       index_r2x_Forr_rofi   = mct_aVect_indexRA(r2x_r,'Forr_rofi')
       index_r2x_Firr_rofi   = mct_aVect_indexRA(r2x_r,'Firr_rofi')
       index_r2x_Flrr_flood  = mct_aVect_indexRA(r2x_r,'Flrr_flood')

       if ( flds_wiso_rof )then
          index_r2x_Forr_rofl_16O   = mct_aVect_indexRA(r2x_r,'Forr_rofl_16O')
          index_r2x_Forr_rofl_18O   = mct_aVect_indexRA(r2x_r,'Forr_rofl_18O')
          index_r2x_Forr_rofl_HDO   = mct_aVect_indexRA(r2x_r,'Forr_rofl_HDO')
          index_r2x_Forr_rofi_16O   = mct_aVect_indexRA(r2x_r,'Forr_rofi_16O')
          index_r2x_Forr_rofi_18O   = mct_aVect_indexRA(r2x_r,'Forr_rofi_18O')
          index_r2x_Forr_rofi_HDO   = mct_aVect_indexRA(r2x_r,'Forr_rofi_HDO')
          index_r2x_Flrr_flood_16O  = mct_aVect_indexRA(r2x_r,'Flrr_flood_16O')
          index_r2x_Flrr_flood_18O  = mct_aVect_indexRA(r2x_r,'Flrr_flood_18O')
          index_r2x_Flrr_flood_HDO  = mct_aVect_indexRA(r2x_r,'Flrr_flood_HDO')
       end if
    end if

    ip = p_inst
    ic = c_rof_rs
    kArea = mct_aVect_indexRA(dom_r%data,afldname)
    lSize = mct_avect_lSize(r2x_r)
    do n=1,lSize
       ca_r =  dom_r%data%rAttr(kArea,n)
       nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*r2x_r%rAttr(index_r2x_Forr_rofl,n) &
            + ca_r*r2x_r%rAttr(index_r2x_Flrr_flood,n)
       nf = f_wioff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*r2x_r%rAttr(index_r2x_Forr_rofi,n) &
            - ca_r*r2x_r%rAttr(index_r2x_Firr_rofi,n)

       if ( flds_wiso_rof )then
          nf = f_wroff_16O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
               ca_r*r2x_r%rAttr(index_r2x_Forr_rofl_16O,n)
          nf = f_wroff_18O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
               ca_r*r2x_r%rAttr(index_r2x_Forr_rofl_18O,n)
          nf = f_wroff_HDO;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
               ca_r*r2x_r%rAttr(index_r2x_Forr_rofl_HDO,n)

          nf = f_wioff_16O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
               ca_r*r2x_r%rAttr(index_r2x_Forr_rofi_16O,n)
          nf = f_wioff_18O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
               ca_r*r2x_r%rAttr(index_r2x_Forr_rofi_18O,n)
          nf = f_wioff_HDO;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
               ca_r*r2x_r%rAttr(index_r2x_Forr_rofi_HDO,n)

          nf = f_wroff_16O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*r2x_r%rAttr(index_r2x_Flrr_flood_16O,n)
          nf = f_wroff_18O;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*r2x_r%rAttr(index_r2x_Flrr_flood_18O,n)
          nf = f_wroff_HDO;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               ca_r*r2x_r%rAttr(index_r2x_Flrr_flood_HDO,n)
       end if
    end do
    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

    first_time = .false.

  end subroutine seq_diag_rof_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_glc_mct - compute global glc input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global glc input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2008-jul-10 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_glc_mct( glc, frac_g, infodata)

    type(component_type)    , intent(in) :: glc    ! component type for instance1
    type(mct_aVect)         , intent(in) :: frac_g ! frac bundle
    type(seq_infodata_type) , intent(in) :: infodata

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: g2x_g
    type(mct_aVect), pointer :: x2g_g
    type(mct_ggrid), pointer :: dom_g
    integer(in)              :: n,ic,nf,ip      ! generic index
    integer(in)              :: kArea             ! index of area field in aVect
    integer(in)              :: lSize             ! size of aVect
    real(r8)                 :: ca_g ! area of a grid cell
    logical,save             :: first_time = .true.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_glc_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    dom_g => component_get_dom_cx(glc)
    g2x_g => component_get_c2x_cx(glc)
    x2g_g => component_get_x2c_cx(glc)

    if (first_time) then
       index_g2x_Fogg_rofl   = mct_aVect_indexRA(g2x_g,'Fogg_rofl')
       index_g2x_Fogg_rofi   = mct_aVect_indexRA(g2x_g,'Fogg_rofi')
       index_g2x_Figg_rofi   = mct_aVect_indexRA(g2x_g,'Figg_rofi')
    end if

    ip = p_inst
    ic = c_glc_gs
    kArea = mct_aVect_indexRA(dom_g%data,afldname)
    lSize = mct_avect_lSize(g2x_g)
    do n=1,lSize
       ca_g =  dom_g%data%rAttr(kArea,n)
       nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_g*g2x_g%rAttr(index_g2x_Fogg_rofl,n)
       nf = f_wioff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_g*g2x_g%rAttr(index_g2x_Fogg_rofi,n) &
            - ca_g*g2x_g%rAttr(index_g2x_Figg_rofi,n)
    end do
    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

    first_time = .false.

  end subroutine seq_diag_glc_mct

  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_ocn_mct - compute global ocn input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global ocn input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2008-jul-10 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_ocn_mct( ocn, xao_o, frac_o, infodata, do_o2x, do_x2o, do_xao)

    type(component_type)    , intent(in)          :: ocn    ! component type for instance1
    type(mct_aVect)         , intent(in)          :: frac_o ! frac bundle
    type(mct_aVect)         , intent(in)          :: xao_o
    type(seq_infodata_type) , intent(in)          :: infodata
    logical                 , intent(in),optional :: do_o2x
    logical                 , intent(in),optional :: do_x2o
    logical                 , intent(in),optional :: do_xao

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: o2x_o        ! model to drv bundle
    type(mct_aVect), pointer :: x2o_o        ! drv to model bundle
    type(mct_ggrid), pointer :: dom_o
    integer(in)              :: n,nf,ic,ip ! generic index
    integer(in)              :: kArea        ! index of area field in aVect
    integer(in)              :: ko,ki  ! fraction indices
    integer(in)              :: lSize        ! size of aVect
    real(r8)                 :: ca_i,ca_o  ! area of a grid cell
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_ocn = .false.
    character(len=cs)        :: cime_model

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_ocn_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (.not. present(do_o2x) .and. &
         .not. present(do_x2o) .and. &
         .not. present(do_xao)) then
       call shr_sys_abort(subName//"ERROR: must input a bundle")
    end if

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    dom_o => component_get_dom_cx(ocn)
    o2x_o => component_get_c2x_cx(ocn)
    x2o_o => component_get_x2c_cx(ocn)

    ip = p_inst

    kArea = mct_aVect_indexRA(dom_o%data,afldname)
    ko    = mct_aVect_indexRA(frac_o,ofracname)
    ki    = mct_aVect_indexRA(frac_o,ifracname)

    call seq_infodata_GetData(infodata, cime_model=cime_model)

    if (present(do_o2x)) then
       if (first_time) then
          if (trim(cime_model) == 'e3sm') then
             index_o2x_Fioo_frazil = mct_aVect_indexRA(o2x_o,'Fioo_frazil')
          else if (trim(cime_model) == 'cesm') then
             index_o2x_Fioo_q = mct_aVect_indexRA(o2x_o,'Fioo_q')
          end if
       end if

       lSize = mct_avect_lSize(o2x_o)
       ic = c_ocn_or
       do n=1,lSize
          ca_o =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
          ca_i =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
          nf = f_area; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o
          if (trim(cime_model) == 'e3sm') then
             nf = f_wfrz; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - (ca_o+ca_i)*max(0.0_r8,o2x_o%rAttr(index_o2x_Fioo_frazil,n))
          else if (trim(cime_model) == 'cesm') then
             nf = f_hfrz; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*max(0.0_r8,o2x_o%rAttr(index_o2x_Fioo_q,n))
          end if
       end do
       if (trim(cime_model) == 'e3sm') then
          budg_dataL(f_hfrz,ic,ip) = -budg_dataL(f_wfrz,ic,ip) * shr_const_latice
       else if (trim(cime_model) == 'cesm') then
          budg_dataL(f_wfrz,ic,ip) = budg_dataL(f_hfrz,ic,ip) * HFLXtoWFLX
       end if
    end if

    if (present(do_xao)) then
       if (first_time) then
          index_xao_Faox_lwup   = mct_aVect_indexRA(xao_o,'Faox_lwup')
          index_xao_Faox_lat    = mct_aVect_indexRA(xao_o,'Faox_lat')
          index_xao_Faox_sen    = mct_aVect_indexRA(xao_o,'Faox_sen')
          index_xao_Faox_evap   = mct_aVect_indexRA(xao_o,'Faox_evap')

          index_xao_Faox_evap_16O = mct_aVect_indexRA(xao_o,'Faox_evap_16O',perrWith='quiet')
          if ( index_xao_Faox_evap_16O /= 0 ) flds_wiso_ocn = .true.
          if ( flds_wiso_ocn )then
             flds_wiso = .true.
             index_xao_Faox_evap_18O = mct_aVect_indexRA(xao_o,'Faox_evap_18O')
             index_xao_Faox_evap_HDO = mct_aVect_indexRA(xao_o,'Faox_evap_HDO')
          end if
       end if

       lSize = mct_avect_lSize(xao_o)
       ic = c_ocn_or
       do n=1,lSize
          ca_o =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
          nf = f_hlwup; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*xao_o%rAttr(index_xao_Faox_lwup,n)
          nf = f_hlatv; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*xao_o%rAttr(index_xao_Faox_lat,n)
          nf = f_hsen ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*xao_o%rAttr(index_xao_Faox_sen,n)
          nf = f_wevap; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*xao_o%rAttr(index_xao_Faox_evap,n)

          if ( flds_wiso_ocn )then
             nf = f_wevap_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_o*xao_o%rAttr(index_xao_Faox_evap_16O,n)
             nf = f_wevap_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_o*xao_o%rAttr(index_xao_Faox_evap_18O,n)
             nf = f_wevap_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_o*xao_o%rAttr(index_xao_Faox_evap_HDO,n)
          end if

       end do
    end if

    if (present(do_x2o)) then
       if (first_time) then
          index_x2o_Fioi_melth  = mct_aVect_indexRA(x2o_o,'Fioi_melth')
          index_x2o_Fioi_meltw  = mct_aVect_indexRA(x2o_o,'Fioi_meltw')
          index_x2o_Fioi_bergh  = mct_aVect_indexRA(x2o_o,'PFioi_bergh', perrWith='quiet')
          index_x2o_Fioi_bergw  = mct_aVect_indexRA(x2o_o,'PFioi_bergw', perrWith='quiet')
          index_x2o_Fioi_salt   = mct_aVect_indexRA(x2o_o,'Fioi_salt')
          index_x2o_Foxx_swnet  = mct_aVect_indexRA(x2o_o,'Foxx_swnet')
          index_x2o_Faxa_lwdn   = mct_aVect_indexRA(x2o_o,'Faxa_lwdn')
          index_x2o_Faxa_rain   = mct_aVect_indexRA(x2o_o,'Faxa_rain')
          index_x2o_Faxa_snow   = mct_aVect_indexRA(x2o_o,'Faxa_snow')
          index_x2o_Foxx_lwup   = mct_aVect_indexRA(x2o_o,'Foxx_lwup')
          index_x2o_Foxx_lat    = mct_aVect_indexRA(x2o_o,'Foxx_lat')
          index_x2o_Foxx_sen    = mct_aVect_indexRA(x2o_o,'Foxx_sen')
          index_x2o_Foxx_evap   = mct_aVect_indexRA(x2o_o,'Foxx_evap')
          index_x2o_Foxx_rofl   = mct_aVect_indexRA(x2o_o,'Foxx_rofl')
          index_x2o_Foxx_rofi   = mct_aVect_indexRA(x2o_o,'Foxx_rofi')

          if ( flds_wiso_ocn )then
             index_x2o_Fioi_meltw_16O = mct_aVect_indexRA(x2o_o,'Fioi_meltw_16O')
             index_x2o_Fioi_meltw_18O = mct_aVect_indexRA(x2o_o,'Fioi_meltw_18O')
             index_x2o_Fioi_meltw_HDO = mct_aVect_indexRA(x2o_o,'Fioi_meltw_HDO')
             index_x2o_Faxa_rain_16O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_16O')
             index_x2o_Faxa_rain_18O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_18O')
             index_x2o_Faxa_rain_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_rain_HDO')
             index_x2o_Faxa_snow_16O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_16O')
             index_x2o_Faxa_snow_18O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_18O')
             index_x2o_Faxa_snow_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_snow_HDO')

             index_x2o_Foxx_rofl_16O  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_16O')
             index_x2o_Foxx_rofi_16O  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_16O')
             index_x2o_Foxx_rofl_18O  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_18O')
             index_x2o_Foxx_rofi_18O  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_18O')
             index_x2o_Foxx_rofl_HDO  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_HDO')
             index_x2o_Foxx_rofi_HDO  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_HDO')
          end if
       end if

       if (.not. present(do_xao)) then
          ! these are in x2o but they really are the atm/ocean flux
          ! computed in the coupler and are "like" an o2x
          lSize = mct_avect_lSize(x2o_o)
          ic = c_ocn_or
          do n=1,lSize
             ca_o =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
             ca_i =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
             nf = f_hlwup; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_lwup,n)
             nf = f_hlatv; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_lat,n)
             nf = f_hsen ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_sen,n)
             nf = f_wevap; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_evap,n)
          end do
       endif

       lSize = mct_avect_lSize(x2o_o)
       ic = c_ocn_os
       do n=1,lSize
          ca_o =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
          ca_i =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o

          if (index_x2o_Fioi_bergw == 0) then
             nf = f_wmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Fioi_meltw,n)
          else
             nf = f_wmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*(x2o_o%rAttr(index_x2o_Fioi_meltw,n)+x2o_o%rAttr(index_x2o_Fioi_bergw,n))
          endif

          if (index_x2o_Fioi_bergh == 0) then
             nf = f_hmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Fioi_melth,n)
          else
             nf = f_hmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*(x2o_o%rAttr(index_x2o_Fioi_melth,n)+x2o_o%rAttr(index_x2o_Fioi_bergh,n))
          endif

          if (trim(cime_model) == 'cesm') then
             nf = f_wsalt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Fioi_salt,n) * SFLXtoWFLX
          endif
          nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_swnet,n)
          nf = f_hlwdn ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_lwdn,n)
          nf = f_wrain ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_rain,n)
          nf = f_wsnow ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_snow,n)
          nf = f_wroff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofl,n)
          nf = f_wioff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofi,n)

          if ( flds_wiso_ocn )then
             nf = f_wmelt_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Fioi_meltw_16O,n)
             nf = f_wmelt_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Fioi_meltw_18O,n)
             nf = f_wmelt_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Fioi_meltw_HDO,n)

             nf = f_wrain_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_rain_16O,n)
             nf = f_wrain_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_rain_18O,n)
             nf = f_wrain_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_rain_HDO,n)

             nf = f_wsnow_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_snow_16O,n)
             nf = f_wsnow_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_snow_18O,n)
             nf = f_wsnow_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_snow_HDO,n)
             nf = f_wroff_16O ;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofl_16O,n)
             nf = f_wioff_16O ;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofi_16O,n)
             nf = f_wroff_18O ;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofl_18O,n)
             nf = f_wioff_18O ;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofi_18O,n)
             nf = f_wroff_HDO ;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofl_HDO,n)
             nf = f_wioff_HDO ;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Foxx_rofi_HDO,n)
          end if
       end do
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
    end if

    ! EBK -- isotope r2x_Forr_rofl/i?

    first_time = .false.

  end subroutine seq_diag_ocn_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_ice_mct - compute global ice input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global ice input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2008-jul-10 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_ice_mct( ice, frac_i, infodata, do_i2x, do_x2i)

    type(component_type)    , intent(in)           :: ice    ! component type for instance1
    type(mct_aVect)         , intent(in)           :: frac_i ! frac bundle
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_i2x
    logical                 , intent(in), optional :: do_x2i

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: i2x_i        ! model to drv bundle
    type(mct_aVect), pointer :: x2i_i        ! drv to model bundle
    type(mct_ggrid), pointer :: dom_i
    integer(in)              :: n,ic,nf,ip ! generic index
    integer(in)              :: kArea        ! index of area field in aVect
    integer(in)              :: kLat         ! index of lat field in aVect
    integer(in)              :: ko,ki  ! fraction indices
    integer(in)              :: lSize        ! size of aVect
    real(r8)                 :: ca_i,ca_o ! area of a grid cell
    logical,save             :: first_time        = .true.
    logical,save             :: flds_wiso_ice     = .false.
    logical,save             :: flds_wiso_ice_x2i = .false.
    character(len=cs)        :: cime_model

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_ice_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_infodata_GetData(infodata, cime_model=cime_model)

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    dom_i => component_get_dom_cx(ice)
    i2x_i => component_get_c2x_cx(ice)
    x2i_i => component_get_x2c_cx(ice)

    ip = p_inst

    kArea = mct_aVect_indexRA(dom_i%data,afldname)
    kLat  = mct_aVect_indexRA(dom_i%data,latname)
    ki    = mct_aVect_indexRA(frac_i,ifracname)
    ko    = mct_aVect_indexRA(frac_i,ofracname)

    if (present(do_i2x)) then
       index_i2x_Fioi_melth  = mct_aVect_indexRA(i2x_i,'Fioi_melth')
       index_i2x_Fioi_meltw  = mct_aVect_indexRA(i2x_i,'Fioi_meltw')
       index_i2x_Fioi_swpen  = mct_aVect_indexRA(i2x_i,'Fioi_swpen')
       index_i2x_Faii_swnet  = mct_aVect_indexRA(i2x_i,'Faii_swnet')
       index_i2x_Faii_lwup   = mct_aVect_indexRA(i2x_i,'Faii_lwup')
       index_i2x_Faii_lat    = mct_aVect_indexRA(i2x_i,'Faii_lat')
       index_i2x_Faii_sen    = mct_aVect_indexRA(i2x_i,'Faii_sen')
       index_i2x_Faii_evap   = mct_aVect_indexRA(i2x_i,'Faii_evap')
       index_i2x_Fioi_salt   = mct_aVect_indexRA(i2x_i,'Fioi_salt')

       index_i2x_Fioi_meltw_16O   = mct_aVect_indexRA(i2x_i,'Fioi_meltw_16O',perrWith='quiet')
       if ( index_i2x_Fioi_meltw_16O /= 0 ) flds_wiso_ice = .true.
       if ( flds_wiso_ice )then
          flds_wiso = .true.
          index_i2x_Fioi_meltw_18O   = mct_aVect_indexRA(i2x_i,'Fioi_meltw_18O')
          index_i2x_Fioi_meltw_HDO   = mct_aVect_indexRA(i2x_i,'Fioi_meltw_HDO')
          index_i2x_Faii_evap_16O    = mct_aVect_indexRA(i2x_i,'Faii_evap_16O')
          index_i2x_Faii_evap_18O    = mct_aVect_indexRA(i2x_i,'Faii_evap_18O')
          index_i2x_Faii_evap_HDO    = mct_aVect_indexRA(i2x_i,'Faii_evap_HDO')
       end if

       lSize = mct_avect_lSize(i2x_i)
       do n=1,lSize
          if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
             ic = c_inh_ir
          else
             ic = c_ish_ir
          endif
          ca_o =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ko,n)
          ca_i =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i
          nf = f_hmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*i2x_i%rAttr(index_i2x_Fioi_melth,n)
          nf = f_wmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*i2x_i%rAttr(index_i2x_Fioi_meltw,n)
          if (trim(cime_model) == 'cesm') then
             nf = f_wsalt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*i2x_i%rAttr(index_i2x_Fioi_salt,n) * SFLXtoWFLX
          endif
          nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*i2x_i%rAttr(index_i2x_Faii_swnet,n) &
               - ca_i*i2x_i%rAttr(index_i2x_Fioi_swpen,n)
          nf = f_hlwup ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*i2x_i%rAttr(index_i2x_Faii_lwup,n)
          nf = f_hlatv ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*i2x_i%rAttr(index_i2x_Faii_lat,n)
          nf = f_hsen  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*i2x_i%rAttr(index_i2x_Faii_sen,n)
          nf = f_wevap ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*i2x_i%rAttr(index_i2x_Faii_evap,n)

          if ( flds_wiso_ice )then
             nf = f_wmelt_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_i*i2x_i%rAttr(index_i2x_Fioi_meltw_16O,n)
             nf = f_wmelt_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_i*i2x_i%rAttr(index_i2x_Fioi_meltw_18O,n)
             nf = f_wmelt_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  ca_i*i2x_i%rAttr(index_i2x_Fioi_meltw_HDO,n)

             nf = f_wevap_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*i2x_i%rAttr(index_i2x_Faii_evap_16O,n)
             nf = f_wevap_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*i2x_i%rAttr(index_i2x_Faii_evap_18O,n)
             nf = f_wevap_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*i2x_i%rAttr(index_i2x_Faii_evap_HDO,n)
          end if
       end do
    end if

    if (present(do_x2i)) then
       if (first_time) then
          index_x2i_Faxa_lwdn   = mct_aVect_indexRA(x2i_i,'Faxa_lwdn')
          index_x2i_Faxa_rain   = mct_aVect_indexRA(x2i_i,'Faxa_rain')
          index_x2i_Faxa_snow   = mct_aVect_indexRA(x2i_i,'Faxa_snow')
          if (trim(cime_model) == 'e3sm') then
             index_x2i_Fioo_frazil = mct_aVect_indexRA(x2i_i,'Fioo_frazil')
          else if (trim(cime_model) == 'cesm') then
             index_x2i_Fioo_q      = mct_aVect_indexRA(x2i_i,'Fioo_q')
          end if
          index_x2i_Fixx_rofi   = mct_aVect_indexRA(x2i_i,'Fixx_rofi')

          index_x2i_Faxa_rain_16O   = mct_aVect_indexRA(x2i_i,'Faxa_rain_16O', perrWith='quiet')
          if ( index_x2i_Faxa_rain_16O /= 0 ) flds_wiso_ice_x2i = .true.
          if ( flds_wiso_ice_x2i )then
             flds_wiso = .true.
             index_x2i_Faxa_rain_18O   = mct_aVect_indexRA(x2i_i,'Faxa_rain_18O')
             index_x2i_Faxa_rain_HDO   = mct_aVect_indexRA(x2i_i,'Faxa_rain_HDO')
             index_x2i_Faxa_snow_16O   = mct_aVect_indexRA(x2i_i,'Faxa_snow_16O')
             index_x2i_Faxa_snow_18O   = mct_aVect_indexRA(x2i_i,'Faxa_snow_18O')
             index_x2i_Faxa_snow_HDO   = mct_aVect_indexRA(x2i_i,'Faxa_snow_HDO')
          end if
       end if

       lSize = mct_avect_lSize(x2i_i)
       do n=1,lSize
          if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
             ic = c_inh_is
          else
             ic = c_ish_is
          endif
          ca_o =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ko,n)
          ca_i =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
          nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i
          nf = f_hlwdn; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*x2i_i%rAttr(index_x2i_Faxa_lwdn,n)
          nf = f_wrain; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*x2i_i%rAttr(index_x2i_Faxa_rain,n)
          nf = f_wsnow; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*x2i_i%rAttr(index_x2i_Faxa_snow,n)
          nf = f_wioff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*x2i_i%rAttr(index_x2i_Fixx_rofi,n)

          if (trim(cime_model) == 'e3sm') then
             nf = f_wfrz ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  (ca_o+ca_i)*max(0.0_r8,x2i_i%rAttr(index_x2i_Fioo_frazil,n))
          else if (trim(cime_model) == 'cesm') then
             nf = f_hfrz ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
                  (ca_o+ca_i)*max(0.0_r8,x2i_i%rAttr(index_x2i_Fioo_q,n))
          end if
          if ( flds_wiso_ice_x2i )then
             nf  = f_wrain_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*x2i_i%rAttr(index_x2i_Faxa_rain_16O,n)
             nf  = f_wrain_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*x2i_i%rAttr(index_x2i_Faxa_rain_18O,n)
             nf  = f_wrain_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*x2i_i%rAttr(index_x2i_Faxa_rain_HDO,n)

             nf  = f_wsnow_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*x2i_i%rAttr(index_x2i_Faxa_snow_16O,n)
             nf  = f_wsnow_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*x2i_i%rAttr(index_x2i_Faxa_snow_18O,n)
             nf  = f_wsnow_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_i*x2i_i%rAttr(index_x2i_Faxa_snow_HDO,n)
          end if
       end do
       ic = c_inh_is
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
       if (trim(cime_model) == 'e3sm') then
          budg_dataL(f_hfrz ,ic,ip) = -budg_dataL(f_wfrz ,ic,ip)*shr_const_latice
       else if (trim(cime_model) == 'cesm') then
          budg_dataL(f_wfrz ,ic,ip) =  budg_dataL(f_hfrz ,ic,ip)*HFLXtoWFLX
       end if

       ic = c_ish_is
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
       if (trim(cime_model) == 'e3sm') then
          budg_dataL(f_hfrz ,ic,ip) = -budg_dataL(f_wfrz ,ic,ip)*shr_const_latice
       else if (trim(cime_model) == 'cesm') then
          budg_dataL(f_wfrz ,ic,ip) =  budg_dataL(f_hfrz ,ic,ip)*HFLXtoWFLX
       end if
    end if

    first_time = .false.

  end subroutine seq_diag_ice_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_print_mct - print global budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global budget diagnostics.
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diag_print_mct(EClock, stop_alarm, &
       budg_print_inst,  budg_print_daily,  budg_print_month,  &
       budg_print_ann,  budg_print_ltann,  budg_print_ltend, infodata)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) , intent(in) :: EClock
    logical          , intent(in) :: stop_alarm
    integer          , intent(in) :: budg_print_inst
    integer          , intent(in) :: budg_print_daily
    integer          , intent(in) :: budg_print_month
    integer          , intent(in) :: budg_print_ann
    integer          , intent(in) :: budg_print_ltann
    integer          , intent(in) :: budg_print_ltend
    type(seq_infodata_type) , intent(in) :: infodata

    !EOP

    !--- local ---
    integer(in)      :: ic,nf,ip,is ! data array indicies
    integer(in)      :: ica,icl,icn,ics,ico
    integer(in)      :: icar,icxs,icxr,icas
    integer(in)      :: cdate,sec   ! coded date, seconds
    integer(in)      :: yr,mon,day  ! date
    integer(in)      :: iam         ! pe number
    integer(in)      :: plev        ! print level
    logical          :: sumdone     ! has a sum been computed yet
    character(len=40):: str         ! string
    character(len=cs):: cime_model
    real(r8) :: dataGpr (f_size,c_size,p_size) ! values to print, scaled and such
    integer, parameter :: nisotopes = 3
    character(len=5), parameter :: isoname(nisotopes) = (/ 'H216O',   'H218O',   '  HDO'   /)
    integer, parameter          :: iso0(nisotopes)    = (/ f_16O,     f_18O,     f_hdO     /)
    integer, parameter          :: isof(nisotopes)    = (/ f_16O_end, f_18O_end, f_hdO_end /)

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_print_mct) '
    character(*),parameter :: F00   = "('(seq_diag_print_mct) ',4a)"

    !----- formats -----
    character(*),parameter :: FAH="(4a,i9,i6)"
    character(*),parameter :: FA0= "('    ',12x,6(6x,a8,1x))"
    character(*),parameter :: FA1= "('    ',a12,6f15.8)"
    character(*),parameter :: FA0r="('    ',12x,8(6x,a8,1x))"
    character(*),parameter :: FA1r="('    ',a12,8f15.8)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_infodata_GetData(infodata, cime_model=cime_model)

    !-------------------------------------------------------------------------------
    ! print instantaneous budget data
    !-------------------------------------------------------------------------------

    sumdone = .false.
    call seq_comm_setptrs(CPLID,iam=iam)
    call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
         curr_mon=mon,curr_day=day,curr_tod=sec)
    cdate = yr*10000+mon*100+day

    do ip = 1,p_size
       plev = 0
       if (ip == p_inst) then
          plev = max(plev,budg_print_inst)
       endif
       if (ip==p_day .and. sec==0) then
          plev = max(plev,budg_print_daily)
       endif
       if (ip==p_mon .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_month)
       endif
       if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_ann)
       endif
       if (ip==p_inf .and. mon==1 .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_ltann)
       endif
       if (ip==p_inf .and. stop_alarm) then
          plev = max(plev,budg_print_ltend)
       endif

       if (plev > 0) then
          ! ---- doprint ---- doprint ---- doprint ----

          if (.not.sumdone) then
             call seq_diag_sum0_mct()
             dataGpr = budg_dataG
             sumdone = .true.

             !  old budget normalizations (global area and 1e6 for water)
             dataGpr = dataGpr/(4.0_r8*shr_const_pi)
             dataGpr(f_w:f_w_end,:,:) = dataGpr(f_w:f_w_end,:,:) * 1.0e6_r8
             if ( flds_wiso )then
                dataGpr(iso0(1):isof(nisotopes),:,:) = dataGpr(iso0(1):isof(nisotopes),:,:) * 1.0e6_r8
             end if
             dataGpr = dataGpr/budg_ns

             if (iam /= 0) return
          endif

          ! ---------------------------------------------------------
          ! ---- detail atm budgets and breakdown into components ---
          ! ---------------------------------------------------------

          if (plev >= 3) then
             do ic = 1,2
                if (ic == 1) then
                   ica = c_atm_ar
                   icl = c_lnd_ar
                   icn = c_inh_ar
                   ics = c_ish_ar
                   ico = c_ocn_ar
                   str = "ATM_to_CPL"
                elseif (ic == 2) then
                   ica = c_atm_as
                   icl = c_lnd_as
                   icn = c_inh_as
                   ics = c_ish_as
                   ico = c_ocn_as
                   str = "CPL_TO_ATM"
                else
                   call shr_sys_abort(subname//' ERROR in ic index code 411')
                endif

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' AREA BUDGET (m2/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
                write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                do nf = f_a, f_a_end
                   write(logunit,FA1)    fname(nf),dataGpr(nf,ica,ip),dataGpr(nf,icl,ip), &
                        dataGpr(nf,icn,ip),dataGpr(nf,ics,ip),dataGpr(nf,ico,ip), &
                        dataGpr(nf,ica,ip)+dataGpr(nf,icl,ip)+ &
                        dataGpr(nf,icn,ip)+dataGpr(nf,ics,ip)+dataGpr(nf,ico,ip)
                enddo

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
                write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                do nf = f_h, f_h_end
                   write(logunit,FA1)    fname(nf),dataGpr(nf,ica,ip),dataGpr(nf,icl,ip), &
                        dataGpr(nf,icn,ip),dataGpr(nf,ics,ip),dataGpr(nf,ico,ip), &
                        dataGpr(nf,ica,ip)+dataGpr(nf,icl,ip)+ &
                        dataGpr(nf,icn,ip)+dataGpr(nf,ics,ip)+dataGpr(nf,ico,ip)
                enddo
                write(logunit,FA1)    '   *SUM*'   ,sum(dataGpr(f_h:f_h_end,ica,ip)),sum(dataGpr(f_h:f_h_end,icl,ip)), &
                     sum(dataGpr(f_h:f_h_end,icn,ip)),sum(dataGpr(f_h:f_h_end,ics,ip)),sum(dataGpr(f_h:f_h_end,ico,ip)), &
                     sum(dataGpr(f_h:f_h_end,ica,ip))+sum(dataGpr(f_h:f_h_end,icl,ip))+ &
                     sum(dataGpr(f_h:f_h_end,icn,ip))+sum(dataGpr(f_h:f_h_end,ics,ip))+sum(dataGpr(f_h:f_h_end,ico,ip))

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
                write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                do nf = f_w, f_w_end
                   if (nf == f_wsalt .and. trim(cime_model) == 'e3sm') cycle
                   write(logunit,FA1)    fname(nf),dataGpr(nf,ica,ip),dataGpr(nf,icl,ip), &
                        dataGpr(nf,icn,ip),dataGpr(nf,ics,ip),dataGpr(nf,ico,ip), &
                        dataGpr(nf,ica,ip)+dataGpr(nf,icl,ip)+ &
                        dataGpr(nf,icn,ip)+dataGpr(nf,ics,ip)+dataGpr(nf,ico,ip)
                enddo
                write(logunit,FA1)    '   *SUM*'   ,sum(dataGpr(f_w:f_w_end,ica,ip)),sum(dataGpr(f_w:f_w_end,icl,ip)), &
                     sum(dataGpr(f_w:f_w_end,icn,ip)),sum(dataGpr(f_w:f_w_end,ics,ip)),sum(dataGpr(f_w:f_w_end,ico,ip)), &
                     sum(dataGpr(f_w:f_w_end,ica,ip))+sum(dataGpr(f_w:f_w_end,icl,ip))+ &
                     sum(dataGpr(f_w:f_w_end,icn,ip))+sum(dataGpr(f_w:f_w_end,ics,ip))+sum(dataGpr(f_w:f_w_end,ico,ip))

                if ( flds_wiso )then
                   do is = 1, nisotopes
                      write(logunit,*) ' '
                      write(logunit,FAH) subname,trim(str)//' '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
                           trim(pname(ip)),': date = ',cdate,sec
                      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                      do nf = iso0(is), isof(is)
                         write(logunit,FA1)    fname(nf),dataGpr(nf,ica,ip),dataGpr(nf,icl,ip), &
                              dataGpr(nf,icn,ip),dataGpr(nf,ics,ip),dataGpr(nf,ico,ip), &
                              dataGpr(nf,ica,ip)+dataGpr(nf,icl,ip)+ &
                              dataGpr(nf,icn,ip)+dataGpr(nf,ics,ip)+dataGpr(nf,ico,ip)
                      enddo
                      write(logunit,FA1)    '   *SUM*', sum(dataGpr(iso0(is):isof(is),ica,ip)),&
                           sum(dataGpr(iso0(is):isof(is),icl,ip)), &
                           sum(dataGpr(iso0(is):isof(is),icn,ip)),&
                           sum(dataGpr(iso0(is):isof(is),ics,ip)), &
                           sum(dataGpr(iso0(is):isof(is),ico,ip)), &
                           sum(dataGpr(iso0(is):isof(is),ica,ip))+&
                           sum(dataGpr(iso0(is):isof(is),icl,ip))+ &
                           sum(dataGpr(iso0(is):isof(is),icn,ip))+&
                           sum(dataGpr(iso0(is):isof(is),ics,ip))+ &
                           sum(dataGpr(iso0(is):isof(is),ico,ip))
                   end do
                end if

             enddo
          endif   ! plev

          ! ---------------------------------------------------------
          ! ---- detail lnd/ocn/ice component budgets ----
          ! ---------------------------------------------------------

          if (plev >= 2) then
             do ic = 1,4
                if (ic == 1) then
                   icar = c_lnd_ar
                   icxs = c_lnd_ls
                   icxr = c_lnd_lr
                   icas = c_lnd_as
                   str = "LND"
                elseif (ic == 2) then
                   icar = c_ocn_ar
                   icxs = c_ocn_os
                   icxr = c_ocn_or
                   icas = c_ocn_as
                   str = "OCN"
                elseif (ic == 3) then
                   icar = c_inh_ar
                   icxs = c_inh_is
                   icxr = c_inh_ir
                   icas = c_inh_as
                   str = "ICE_NH"
                elseif (ic == 4) then
                   icar = c_ish_ar
                   icxs = c_ish_is
                   icxr = c_ish_ir
                   icas = c_ish_as
                   str = "ICE_SH"
                else
                   call shr_sys_abort(subname//' ERROR in ic index code 412')
                endif

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
                write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                do nf = f_h, f_h_end
                   write(logunit,FA1)    fname(nf),-dataGpr(nf,icar,ip),dataGpr(nf,icxs,ip), &
                        dataGpr(nf,icxr,ip),-dataGpr(nf,icas,ip), &
                        -dataGpr(nf,icar,ip)+dataGpr(nf,icxs,ip)+ &
                        dataGpr(nf,icxr,ip)-dataGpr(nf,icas,ip)
                enddo
                write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_h:f_h_end,icar,ip)),sum(dataGpr(f_h:f_h_end,icxs,ip)), &
                     sum(dataGpr(f_h:f_h_end,icxr,ip)),-sum(dataGpr(f_h:f_h_end,icas,ip)), &
                     -sum(dataGpr(f_h:f_h_end,icar,ip))+sum(dataGpr(f_h:f_h_end,icxs,ip))+ &
                     sum(dataGpr(f_h:f_h_end,icxr,ip))-sum(dataGpr(f_h:f_h_end,icas,ip))

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
                write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                do nf = f_w, f_w_end
                   if (nf == f_wsalt .and. trim(cime_model) == 'e3sm') cycle
                   write(logunit,FA1)    fname(nf),-dataGpr(nf,icar,ip),dataGpr(nf,icxs,ip), &
                        dataGpr(nf,icxr,ip),-dataGpr(nf,icas,ip), &
                        -dataGpr(nf,icar,ip)+dataGpr(nf,icxs,ip)+ &
                        dataGpr(nf,icxr,ip)-dataGpr(nf,icas,ip)
                enddo
                write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_w:f_w_end,icar,ip)),sum(dataGpr(f_w:f_w_end,icxs,ip)), &
                     sum(dataGpr(f_w:f_w_end,icxr,ip)),-sum(dataGpr(f_w:f_w_end,icas,ip)), &
                     -sum(dataGpr(f_w:f_w_end,icar,ip))+sum(dataGpr(f_w:f_w_end,icxs,ip))+ &
                     sum(dataGpr(f_w:f_w_end,icxr,ip))-sum(dataGpr(f_w:f_w_end,icas,ip))

                if ( flds_wiso ) then
                   do is = 1, nisotopes
                      write(logunit,*) ' '
                      write(logunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)), &
                           ': date = ',cdate,sec
                      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                      do nf = iso0(is), isof(is)
                         write(logunit,FA1)    fname(nf),-dataGpr(nf,icar,ip),dataGpr(nf,icxs,ip), &
                              dataGpr(nf,icxr,ip),-dataGpr(nf,icas,ip), &
                              -dataGpr(nf,icar,ip)+dataGpr(nf,icxs,ip)+ &
                              dataGpr(nf,icxr,ip)-dataGpr(nf,icas,ip)
                      enddo
                      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(iso0(is):isof(is),icar,ip)),&
                           sum(dataGpr(iso0(is):isof(is),icxs,ip)), &
                           sum(dataGpr(iso0(is):isof(is),icxr,ip)), &
                           -sum(dataGpr(iso0(is):isof(is),icas,ip)), &
                           -sum(dataGpr(iso0(is):isof(is),icar,ip)) &
                           +sum(dataGpr(iso0(is):isof(is),icxs,ip))+ &
                           sum(dataGpr(iso0(is):isof(is),icxr,ip)) &
                           -sum(dataGpr(iso0(is):isof(is),icas,ip))
                      write(logunit,*) ' '
                      write(logunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),&
                           ': date = ',cdate,sec
                      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                      do nf = iso0(is), isof(is)
                         write(logunit,FA1)    fname(nf),-dataGpr(nf,icar,ip),dataGpr(nf,icxs,ip), &
                              dataGpr(nf,icxr,ip),-dataGpr(nf,icas,ip), &
                              -dataGpr(nf,icar,ip)+dataGpr(nf,icxs,ip)+ &
                              dataGpr(nf,icxr,ip)-dataGpr(nf,icas,ip)
                      enddo
                      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(iso0(is):isof(is),icar,ip)),&
                           sum(dataGpr(iso0(is):isof(is),icxs,ip)), &
                           sum(dataGpr(iso0(is):isof(is),icxr,ip)), &
                           -sum(dataGpr(iso0(is):isof(is),icas,ip)), &
                           -sum(dataGpr(iso0(is):isof(is),icar,ip)) &
                           +sum(dataGpr(iso0(is):isof(is),icxs,ip))+ &
                           sum(dataGpr(iso0(is):isof(is),icxr,ip)) &
                           -sum(dataGpr(iso0(is):isof(is),icas,ip))
                   end do
                end if
             enddo
          endif   ! plev

          ! ---------------------------------------------------------
          ! ---- net summary budgets ----
          ! ---------------------------------------------------------

          if (plev >= 1) then

             write(logunit,*) ' '
             write(logunit,FAH) subname,'NET AREA BUDGET (m2/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
             write(logunit,FA0) '     atm','     lnd','     ocn','  ice nh','  ice sh',' *SUM*  '
             do nf = f_a,f_a_end
                write(logunit,FA1)    fname(nf),dataGpr(nf,c_atm_ar,ip), &
                     dataGpr(nf,c_lnd_lr,ip), &
                     dataGpr(nf,c_ocn_or,ip), &
                     dataGpr(nf,c_inh_ir,ip), &
                     dataGpr(nf,c_ish_ir,ip), &
                     dataGpr(nf,c_atm_ar,ip)+ &
                     dataGpr(nf,c_lnd_lr,ip)+ &
                     dataGpr(nf,c_ocn_or,ip)+ &
                     dataGpr(nf,c_inh_ir,ip)+ &
                     dataGpr(nf,c_ish_ir,ip)
             enddo

             write(logunit,*) ' '
             write(logunit,FAH) subname,'NET HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
             write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
             do nf = f_h, f_h_end
                write(logunit,FA1r)   fname(nf),dataGpr(nf,c_atm_ar,ip)+dataGpr(nf,c_atm_as,ip), &
                     dataGpr(nf,c_lnd_lr,ip)+dataGpr(nf,c_lnd_ls,ip), &
                     dataGpr(nf,c_rof_rr,ip)+dataGpr(nf,c_rof_rs,ip), &
                     dataGpr(nf,c_ocn_or,ip)+dataGpr(nf,c_ocn_os,ip), &
                     dataGpr(nf,c_inh_ir,ip)+dataGpr(nf,c_inh_is,ip), &
                     dataGpr(nf,c_ish_ir,ip)+dataGpr(nf,c_ish_is,ip), &
                     dataGpr(nf,c_glc_gr,ip)+dataGpr(nf,c_glc_gs,ip), &
                     dataGpr(nf,c_atm_ar,ip)+dataGpr(nf,c_atm_as,ip)+ &
                     dataGpr(nf,c_lnd_lr,ip)+dataGpr(nf,c_lnd_ls,ip)+ &
                     dataGpr(nf,c_rof_rr,ip)+dataGpr(nf,c_rof_rs,ip)+ &
                     dataGpr(nf,c_ocn_or,ip)+dataGpr(nf,c_ocn_os,ip)+ &
                     dataGpr(nf,c_inh_ir,ip)+dataGpr(nf,c_inh_is,ip)+ &
                     dataGpr(nf,c_ish_ir,ip)+dataGpr(nf,c_ish_is,ip)+ &
                     dataGpr(nf,c_glc_gr,ip)+dataGpr(nf,c_glc_gs,ip)
             enddo
             write(logunit,FA1r)'   *SUM*',sum(dataGpr(f_h:f_h_end,c_atm_ar,ip))+sum(dataGpr(f_h:f_h_end,c_atm_as,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_lnd_lr,ip))+sum(dataGpr(f_h:f_h_end,c_lnd_ls,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_rof_rr,ip))+sum(dataGpr(f_h:f_h_end,c_rof_rs,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_ocn_or,ip))+sum(dataGpr(f_h:f_h_end,c_ocn_os,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_inh_ir,ip))+sum(dataGpr(f_h:f_h_end,c_inh_is,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_ish_ir,ip))+sum(dataGpr(f_h:f_h_end,c_ish_is,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_glc_gr,ip))+sum(dataGpr(f_h:f_h_end,c_glc_gs,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_atm_ar,ip))+sum(dataGpr(f_h:f_h_end,c_atm_as,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_lnd_lr,ip))+sum(dataGpr(f_h:f_h_end,c_lnd_ls,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_rof_rr,ip))+sum(dataGpr(f_h:f_h_end,c_rof_rs,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_ocn_or,ip))+sum(dataGpr(f_h:f_h_end,c_ocn_os,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_inh_ir,ip))+sum(dataGpr(f_h:f_h_end,c_inh_is,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_ish_ir,ip))+sum(dataGpr(f_h:f_h_end,c_ish_is,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_glc_gr,ip))+sum(dataGpr(f_h:f_h_end,c_glc_gs,ip))

             write(logunit,*) ' '
             write(logunit,FAH) subname,'NET WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
             write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
             do nf = f_w, f_w_end
                if (nf == f_wsalt .and. trim(cime_model) == 'e3sm') cycle
                write(logunit,FA1r)   fname(nf),dataGpr(nf,c_atm_ar,ip)+dataGpr(nf,c_atm_as,ip), &
                     dataGpr(nf,c_lnd_lr,ip)+dataGpr(nf,c_lnd_ls,ip), &
                     dataGpr(nf,c_rof_rr,ip)+dataGpr(nf,c_rof_rs,ip), &
                     dataGpr(nf,c_ocn_or,ip)+dataGpr(nf,c_ocn_os,ip), &
                     dataGpr(nf,c_inh_ir,ip)+dataGpr(nf,c_inh_is,ip), &
                     dataGpr(nf,c_ish_ir,ip)+dataGpr(nf,c_ish_is,ip), &
                     dataGpr(nf,c_glc_gr,ip)+dataGpr(nf,c_glc_gs,ip), &
                     dataGpr(nf,c_atm_ar,ip)+dataGpr(nf,c_atm_as,ip)+ &
                     dataGpr(nf,c_lnd_lr,ip)+dataGpr(nf,c_lnd_ls,ip)+ &
                     dataGpr(nf,c_rof_rr,ip)+dataGpr(nf,c_rof_rs,ip)+ &
                     dataGpr(nf,c_ocn_or,ip)+dataGpr(nf,c_ocn_os,ip)+ &
                     dataGpr(nf,c_inh_ir,ip)+dataGpr(nf,c_inh_is,ip)+ &
                     dataGpr(nf,c_ish_ir,ip)+dataGpr(nf,c_ish_is,ip)+ &
                     dataGpr(nf,c_glc_gr,ip)+dataGpr(nf,c_glc_gs,ip)
             enddo
             write(logunit,FA1r)'   *SUM*',sum(dataGpr(f_w:f_w_end,c_atm_ar,ip))+sum(dataGpr(f_w:f_w_end,c_atm_as,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_lnd_lr,ip))+sum(dataGpr(f_w:f_w_end,c_lnd_ls,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_rof_rr,ip))+sum(dataGpr(f_w:f_w_end,c_rof_rs,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_ocn_or,ip))+sum(dataGpr(f_w:f_w_end,c_ocn_os,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_inh_ir,ip))+sum(dataGpr(f_w:f_w_end,c_inh_is,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_ish_ir,ip))+sum(dataGpr(f_w:f_w_end,c_ish_is,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_glc_gr,ip))+sum(dataGpr(f_w:f_w_end,c_glc_gs,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_atm_ar,ip))+sum(dataGpr(f_w:f_w_end,c_atm_as,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_lnd_lr,ip))+sum(dataGpr(f_w:f_w_end,c_lnd_ls,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_rof_rr,ip))+sum(dataGpr(f_w:f_w_end,c_rof_rs,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_ocn_or,ip))+sum(dataGpr(f_w:f_w_end,c_ocn_os,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_inh_ir,ip))+sum(dataGpr(f_w:f_w_end,c_inh_is,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_ish_ir,ip))+sum(dataGpr(f_w:f_w_end,c_ish_is,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_glc_gr,ip))+sum(dataGpr(f_w:f_w_end,c_glc_gs,ip))

             if ( flds_wiso ) then

                do is = 1, nisotopes
                   write(logunit,*) ' '
                   write(logunit,FAH) subname,'NET '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
                        trim(pname(ip)),': date = ',cdate,sec
                   write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
                   do nf = iso0(is), isof(is)
                      write(logunit,FA1r)   fname(nf),dataGpr(nf,c_atm_ar,ip)+dataGpr(nf,c_atm_as,ip), &
                           dataGpr(nf,c_lnd_lr,ip)+dataGpr(nf,c_lnd_ls,ip), &
                           dataGpr(nf,c_rof_rr,ip)+dataGpr(nf,c_rof_rs,ip), &
                           dataGpr(nf,c_ocn_or,ip)+dataGpr(nf,c_ocn_os,ip), &
                           dataGpr(nf,c_inh_ir,ip)+dataGpr(nf,c_inh_is,ip), &
                           dataGpr(nf,c_ish_ir,ip)+dataGpr(nf,c_ish_is,ip), &
                           dataGpr(nf,c_glc_gr,ip)+dataGpr(nf,c_glc_gs,ip), &
                           dataGpr(nf,c_atm_ar,ip)+dataGpr(nf,c_atm_as,ip)+ &
                           dataGpr(nf,c_lnd_lr,ip)+dataGpr(nf,c_lnd_ls,ip)+ &
                           dataGpr(nf,c_rof_rr,ip)+dataGpr(nf,c_rof_rs,ip)+ &
                           dataGpr(nf,c_ocn_or,ip)+dataGpr(nf,c_ocn_os,ip)+ &
                           dataGpr(nf,c_inh_ir,ip)+dataGpr(nf,c_inh_is,ip)+ &
                           dataGpr(nf,c_ish_ir,ip)+dataGpr(nf,c_ish_is,ip)+ &
                           dataGpr(nf,c_glc_gr,ip)+dataGpr(nf,c_glc_gs,ip)
                   enddo
                   write(logunit,FA1r)'   *SUM*',&
                        sum(dataGpr(iso0(is):isof(is),c_atm_ar,ip))+ &
                        sum(dataGpr(iso0(is):isof(is),c_atm_as,ip)),&
                        sum(dataGpr(iso0(is):isof(is),c_lnd_lr,ip))+ &
                        sum(dataGpr(iso0(is):isof(is),c_lnd_ls,ip)),&
                        sum(dataGpr(iso0(is):isof(is),c_rof_rr,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_rof_rs,ip)),&
                        sum(dataGpr(iso0(is):isof(is),c_ocn_or,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_ocn_os,ip)),&
                        sum(dataGpr(iso0(is):isof(is),c_inh_ir,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_inh_is,ip)),&
                        sum(dataGpr(iso0(is):isof(is),c_ish_ir,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_ish_is,ip)),&
                        sum(dataGpr(iso0(is):isof(is),c_glc_gr,ip))+ &
                        sum(dataGpr(iso0(is):isof(is),c_glc_gs,ip)),&
                        sum(dataGpr(iso0(is):isof(is),c_atm_ar,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_atm_as,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_lnd_lr,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_lnd_ls,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_rof_rr,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_rof_rs,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_ocn_or,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_ocn_os,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_inh_ir,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_inh_is,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_ish_ir,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_ish_is,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_glc_gr,ip))+&
                        sum(dataGpr(iso0(is):isof(is),c_glc_gs,ip))
                end do
             end if

          endif

          write(logunit,*) ' '
          ! ---- doprint ---- doprint ---- doprint ----
       endif  ! plev > 0
    enddo  ! ip = 1,p_size

  end subroutine seq_diag_print_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_avect_mct - print global budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global diagnostics for AV/ID.
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diag_avect_mct(infodata, id, av, dom, gsmap, comment)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(seq_infodata_type) , intent(in)           :: infodata
    integer(in)             , intent(in)           :: ID
    type(mct_aVect)         , intent(in)           :: av
    type(mct_gGrid)         , pointer              :: dom
    type(mct_gsMap)         , pointer              :: gsmap
    character(len=*)        , intent(in), optional :: comment

    !EOP

    !--- local ---
    logical                          :: bfbflag
    integer(in)                      :: n,k         ! counters
    integer(in)                      :: npts,nptsg  ! number of local/global pts in AV
    integer(in)                      :: kflds       ! number of fields in AV
    real(r8),                pointer :: sumbuf (:)  ! sum buffer
    real(r8),                pointer :: maxbuf (:)  ! max buffer
    real(r8),                pointer :: sumbufg(:)  ! sum buffer reduced
    real(r8),                pointer :: maxbufg(:)  ! max buffer reduced
    integer(i8),             pointer :: isumbuf (:) ! integer local sum
    integer(i8),             pointer :: isumbufg(:) ! integer global sum
    integer(i8)                      :: ihuge       ! huge
    integer(in)                      :: mpicom      ! mpi comm
    integer(in)                      :: iam         ! pe number
    integer(in)                      :: km,ka       ! field indices
    integer(in)                      :: ns          ! size of local AV
    integer(in)                      :: rcode       ! allocate return code
    real(r8),                pointer :: weight(:)   ! weight
    real(r8), allocatable            :: weighted_data(:,:) ! weighted data
    type(mct_string)                 :: mstring     ! mct char type
    character(CL)                    :: lcomment    ! should be long enough
    character(CL)                    :: itemc       ! string converted to char

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_avect_mct) '
    character(*),parameter :: F00   = "('(seq_diag_avect_mct) ',4a)"

    !-------------------------------------------------------------------------------
    ! print instantaneous budget data
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(ID, mpicom=mpicom, iam=iam)
    call seq_infodata_GetData(infodata, bfbflag=bfbflag)

    lcomment = ''
    if (present(comment)) then
       lcomment=trim(comment)
    endif

    ns = mct_aVect_lsize(AV)
    npts = mct_aVect_lsize(dom%data)
    if (ns /= npts) call shr_sys_abort(trim(subname)//' ERROR: size of AV,dom')
    km = mct_aVect_indexRA(dom%data,'mask')
    ka = mct_aVect_indexRA(dom%data,afldname)
    kflds = mct_aVect_nRattr(AV)
    allocate(sumbufg(kflds),stat=rcode)
    if (rcode /= 0) call shr_sys_abort(trim(subname)//' allocate sumbufg')

    npts = mct_aVect_lsize(AV)
    allocate(weight(npts),stat=rcode)
    if (rcode /= 0) call shr_sys_abort(trim(subname)//' allocate weight')

    weight(:) = 1.0_r8
    do n = 1,npts
       if (dom%data%rAttr(km,n) <= 1.0e-06_R8) then
          weight(n) = 0.0_r8
       else
          weight(n) = dom%data%rAttr(ka,n)*shr_const_rearth*shr_const_rearth
       endif
    enddo

    if (bfbflag) then
       allocate(weighted_data(npts,kflds),stat=rcode)
       if (rcode /= 0) call shr_sys_abort(trim(subname)//' allocate weighted_data')

       weighted_data = 0.0_r8
       do n = 1,npts
          do k = 1,kflds
             if (.not. shr_const_isspval(AV%rAttr(k,n))) then
                weighted_data(n,k) = AV%rAttr(k,n)*weight(n)
             endif
          enddo
       enddo

       call shr_reprosum_calc (weighted_data, sumbufg, npts, npts, kflds, &
                               commid=mpicom)

       deallocate(weighted_data)

    else
       allocate(sumbuf(kflds),stat=rcode)
       if (rcode /= 0) call shr_sys_abort(trim(subname)//' allocate sumbuf')
       sumbuf = 0.0_r8

       do n = 1,npts
          do k = 1,kflds
             if (.not. shr_const_isspval(AV%rAttr(k,n))) then
                sumbuf(k) = sumbuf(k) + AV%rAttr(k,n)*weight(n)
             endif
          enddo
       enddo

       !--- global reduction ---
       call shr_mpi_sum(sumbuf,sumbufg,mpicom,subname)

       deallocate(sumbuf)

    endif
    deallocate(weight)

    if (iam == 0) then
       !      write(logunit,*) 'sdAV: *** writing ',trim(lcomment),': k fld min/max/sum ***'
       do k = 1,kflds
          call mct_aVect_getRList(mstring,k,AV)
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
          if (len_trim(lcomment) > 0) then
             write(logunit,100) 'xxx','sorr',k,sumbufg(k),trim(lcomment),trim(itemc)
          else
             write(logunit,101) 'xxx','sorr',k,sumbufg(k),trim(itemc)
          endif
       enddo
       call shr_sys_flush(logunit)
    endif

    deallocate(sumbufg)

100 format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a,1x,a)
101 format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a)

  end subroutine seq_diag_avect_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_avloc_mct - print local budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print local diagnostics for AV/ID.
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diag_avloc_mct(av, comment)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect) , intent(in)           :: av
    character(len=*), intent(in), optional :: comment

    !EOP

    !--- local ---
    integer(in)                      :: n,k         ! counters
    integer(in)                      :: npts        ! number of local/global pts in AV
    integer(in)                      :: kflds       ! number of fields in AV
    real(r8),                pointer :: sumbuf (:)  ! sum buffer
    type(mct_string)                 :: mstring     ! mct char type
    character(CL)                    :: lcomment    ! should be long enough
    character(CL)                    :: itemc       ! string converted to char

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_avloc_mct) '
    character(*),parameter :: F00   = "('(seq_diag_avloc_mct) ',4a)"

    !-------------------------------------------------------------------------------
    ! print instantaneous budget data
    !-------------------------------------------------------------------------------

    lcomment = ''
    if (present(comment)) then
       lcomment=trim(comment)
    endif

    npts = mct_aVect_lsize(AV)
    kflds = mct_aVect_nRattr(AV)
    allocate(sumbuf(kflds))

    sumbuf = 0.0_r8
    do n = 1,npts
       do k = 1,kflds
          !      if (.not. shr_const_isspval(AV%rAttr(k,n))) then
          sumbuf(k) = sumbuf(k) + AV%rAttr(k,n)
          !      endif
       enddo
    enddo

    do k = 1,kflds
       call mct_aVect_getRList(mstring,k,AV)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       if (len_trim(lcomment) > 0) then
          write(logunit,100) 'xxx','sorr',k,sumbuf(k),trim(lcomment),trim(itemc)
       else
          write(logunit,101) 'xxx','sorr',k,sumbuf(k),trim(itemc)
       endif
    enddo
    call shr_sys_flush(logunit)

    deallocate(sumbuf)

100 format('avloc_diag ',a3,1x,a4,1x,i3,es26.19,1x,a,1x,a)
101 format('avloc_diag ',a3,1x,a4,1x,i3,es26.19,1x,a)

  end subroutine seq_diag_avloc_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_avdiff_mct - print global budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global diagnostics for AV/ID.
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diag_avdiff_mct(AV1,AV2,ID,comment)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect) , intent(in) :: AV1
    type(mct_aVect) , intent(in) :: AV2
    integer         , intent(in) :: ID
    character(len=*), intent(in), optional :: comment

    !EOP

    !--- local ---
    integer(in)      :: n,k,n1,k1,n2,k2         ! counters
    integer(in)      :: iam         ! pe number
    integer(in)      :: cnt         ! counter
    real(r8)         :: adiff,rdiff ! diff values
    type(mct_string) :: mstring     ! mct char type
    character(len=64):: lcomment    ! should be long enough

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_avdiff_mct) '
    character(*),parameter :: F00   = "('(seq_diag_avdiff_mct) ',4a)"

    !-------------------------------------------------------------------------------
    ! print instantaneous budget data
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(ID,iam=iam)

    lcomment = ''
    if (present(comment)) then
       lcomment=trim(comment)
    endif

    n1 = mct_aVect_lsize(AV1)
    k1 = mct_aVect_nRattr(AV1)
    n2 = mct_aVect_lsize(AV2)
    k2 = mct_aVect_nRattr(AV2)

    if (n1 /= n2 .or. k1 /= k2) then
       write(s_logunit,*) subname,trim(lcomment),' AV sizes different ',n1,n2,k1,k2
       return
    endif

    do k = 1,k1
       cnt = 0
       adiff = 0.
       rdiff = 0.
       do n = 1,n1
          if (AV1%rAttr(k,n) /= AV2%rAttr(k,n)) then
             cnt = cnt + 1
             adiff = max(adiff, abs(AV1%rAttr(k,n)-AV2%rAttr(k,n)))
             rdiff = max(rdiff, abs(AV1%rAttr(k,n)-AV2%rAttr(k,n))/(abs(AV1%rAttr(k,n))+abs(AV2%rAttr(k,n))))
          endif
       enddo
       if (cnt > 0) then
          call mct_aVect_getRList(mstring,k,AV1)
          write(s_logunit,*) subname,trim(lcomment),' AVs fld k diff ', &
               iam,mct_string_toChar(mstring),cnt,adiff,rdiff, &
               minval(AV1%rAttr(k,:)),minval(AV1%rAttr(k,:)), &
               maxval(AV1%rAttr(k,:)),maxval(AV2%rAttr(k,:))
          call mct_string_clean(mstring)
       endif
    enddo

  end subroutine seq_diag_avdiff_mct

end module seq_diag_mct
