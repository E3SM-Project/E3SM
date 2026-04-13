!===============================================================================
!
! !MODULE: seq_diag_mod -- computes spatial \& time averages of fluxed quatities
!
! !DESCRIPTION:
!    The coupler is required to do certain diagnostics, those calculations are
!    located in this module.
!
! !REMARKS:
!    E3SM sign convention for fluxes is positive downward with hierarchy being
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
!    2026-agr 12 - R. Jacob    - moab port
!    2012-aug-20 - T. Craig    - add rof component
!    2008-jul-10 - T. Craig    - updated budget implementation
!    2007-may-07 - B. Kauffman - initial port to cpl7.
!    2002-nov-21 - R. Jacob    - initial port to cpl6.
!    199x-mmm-dd - B. Kauffman - original version in cpl4.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_diag_moab
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
  use seq_comm_mct, only: logunit, cplid, seq_comm_setptrs, seq_comm_clean, &
       mbaxid, mblxid, mboxid, mbixid, mbrxid, mbofxid
  use shr_moab_mod, only: mbGetnCells, mbGetCellTagVals
  use seq_timemgr_mod, only : seq_timemgr_EClockGetData
  use component_type_mod, only : COMPONENT_GET_C2X_CX, COMPONENT_GET_X2C_CX, COMPONENT_TYPE
  use seq_infodata_mod, only : seq_infodata_type, seq_infodata_getdata
  use shr_reprosum_mod, only : shr_reprosum_calc
  use seq_diagBGC_moab,  only : seq_diagBGC_preprint_moab, seq_diagBGC_print_moab

  implicit none
  save
  private

  ! !PUBLIC TYPES:

  ! none

  !PUBLIC MEMBER FUNCTIONS:

  public seq_diag_zero_moab
  public seq_diag_atm_moab
  public seq_diag_lnd_moab
  public seq_diag_rof_moab
  public seq_diag_glc_moab
  public seq_diag_ocn_moab
  public seq_diag_ice_moab
  public seq_diag_accum_moab
  public seq_diag_sum0_moab
  public seq_diag_print_moab
  ! seq_diag_avect removed (not needed for MOAB)
  
  

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
  integer(in),parameter :: f_hpolar    =11     ! heat : AIS imbalance
  integer(in),parameter :: f_hh2ot     =12     ! heat : water temperature
  integer(in),parameter :: f_wfrz      =13     ! water: freezing
  integer(in),parameter :: f_wmelt     =14     ! water: melting
  integer(in),parameter :: f_wrain     =15     ! water: precip, liquid
  integer(in),parameter :: f_wsnow     =16     ! water: precip, frozen
  integer(in),parameter :: f_wpolar    =17     ! water: AIS imbalance
  integer(in),parameter :: f_wevap     =18     ! water: evaporation
  integer(in),parameter :: f_wroff     =19     ! water: runoff/flood
  integer(in),parameter :: f_wioff     =20     ! water: frozen runoff
  integer(in),parameter :: f_wirrig    =21     ! water: irrigation
  integer(in),parameter :: f_wfrz_16O  =22     ! water: freezing
  integer(in),parameter :: f_wmelt_16O =23     ! water: melting
  integer(in),parameter :: f_wrain_16O =24     ! water: precip, liquid
  integer(in),parameter :: f_wsnow_16O =25     ! water: precip, frozen
  integer(in),parameter :: f_wevap_16O =26     ! water: evaporation
  integer(in),parameter :: f_wroff_16O =27     ! water: runoff/flood
  integer(in),parameter :: f_wioff_16O =28     ! water: frozen runoff
  integer(in),parameter :: f_wfrz_18O  =29     ! water: freezing
  integer(in),parameter :: f_wmelt_18O =30     ! water: melting
  integer(in),parameter :: f_wrain_18O =31     ! water: precip, liquid
  integer(in),parameter :: f_wsnow_18O =32     ! water: precip, frozen
  integer(in),parameter :: f_wevap_18O =33     ! water: evaporation
  integer(in),parameter :: f_wroff_18O =34     ! water: runoff/flood
  integer(in),parameter :: f_wioff_18O =35     ! water: frozen runoff
  integer(in),parameter :: f_wfrz_HDO  =36     ! water: freezing
  integer(in),parameter :: f_wmelt_HDO =37     ! water: melting
  integer(in),parameter :: f_wrain_HDO =38     ! water: precip, liquid
  integer(in),parameter :: f_wsnow_HDO =39     ! water: precip, frozen
  integer(in),parameter :: f_wevap_HDO =40     ! water: evaporation
  integer(in),parameter :: f_wroff_HDO =41     ! water: runoff/flood
  integer(in),parameter :: f_wioff_HDO =42     ! water: frozen runoff

  integer(in),parameter :: f_size     = f_wioff_HDO   ! Total array size of all elements
  integer(in),parameter :: f_a        = f_area        ! 1st index for area
  integer(in),parameter :: f_a_end    = f_area        ! last index for area
  integer(in),parameter :: f_h        = f_hfrz        ! 1st index for heat
  integer(in),parameter :: f_h_end    = f_hh2ot       ! Last index for heat
  integer(in),parameter :: f_w        = f_wfrz        ! 1st index for water
  integer(in),parameter :: f_w_end    = f_wirrig      ! Last index for water
  integer(in),parameter :: f_16O      = f_wfrz_16O    ! 1st index for 16O water isotope
  integer(in),parameter :: f_18O      = f_wfrz_18O    ! 1st index for 18O water isotope
  integer(in),parameter :: f_HDO      = f_wfrz_HDO    ! 1st index for HDO water isotope
  integer(in),parameter :: f_16O_end  = f_wioff_16O   ! Last index for 16O water isotope
  integer(in),parameter :: f_18O_end  = f_wioff_18O   ! Last index for 18O water isotope
  integer(in),parameter :: f_HDO_end  = f_wioff_HDO   ! Last index for HDO water isotope

  character(len=12),parameter :: fname(f_size) = &

       (/'        area','     hfreeze','       hmelt','      hnetsw','       hlwdn', &
       '       hlwup','     hlatvap','     hlatfus','      hiroff','        hsen', &
       '      hpolar','    hh2otemp','     wfreeze','       wmelt','       wrain', &
       '       wsnow','      wpolar','       wevap','     wrunoff','     wfrzrof', &
       '      wirrig',                                                             &
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
  character(len=*),parameter :: lfrinname = 'lfrin'
  character(len=*),parameter :: ofracname = 'ofrac'
  character(len=*),parameter :: ifracname = 'ifrac'

  character(*),parameter :: modName = "(seq_diag_moab) "

  integer(in),parameter :: debug = 0 ! internal debug level

  ! !PRIVATE DATA MEMBERS

  integer :: index_i2x_Fioi_meltw_16O
  integer :: index_x2i_Faxa_rain_16O
  integer :: index_a2x_Faxa_rainc_16O
  integer :: index_x2o_Faxa_rain_16O
  integer :: index_l2x_Fall_evap_16O
  integer :: index_x2r_Flrl_rofl_16O
  integer :: index_l2x_Flrl_irrig
  integer :: index_x2r_Flrl_irrig
  integer :: index_o2x_Foxo_ismw

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

  subroutine seq_diag_zero_moab(EClock,mode)

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(in),optional :: EClock
    character(len=*), intent(in),optional :: mode

    !EOP

    integer(IN) :: ip,yr,mon,day,sec
    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_zero_moab) '

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

  end subroutine seq_diag_zero_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_accum_moab - accum out global budget diagnostic data.
  !
  ! !DESCRIPTION:
  !    Accum out global budget diagnostic data.
  !
  ! !REVISION HISTORY:
  !    2008-jul-11 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_accum_moab()

    ! !INPUT/OUTPUT PARAMETERS:

    !EOP

    integer(in) :: ip

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_accum_moab) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    do ip = p_inst+1,p_size
       budg_dataL(:,:,ip) = budg_dataL(:,:,ip) + budg_dataL(:,:,p_inst)
    enddo
    budg_ns(:,:,:) = budg_ns(:,:,:) + 1.0_r8

  end subroutine seq_diag_accum_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_sum0_moab - sum local to global on root
  !
  ! !DESCRIPTION:
  !    Sum local values to global on root
  !
  ! !REVISION HISTORY:
  !    2008-jul-19 - T. Craig - update
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_sum0_moab()

    ! !INPUT/OUTPUT PARAMETERS:

    !EOP

    real(r8) :: budg_dataGtmp(f_size,c_size,p_size) ! temporary sum
    integer(in)      :: mpicom      ! mpi comm
    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_sum0_moab) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,mpicom=mpicom)
    budg_dataGtmp = 0.0_r8
    call shr_mpi_sum(budg_dataL,budg_dataGtmp,mpicom,subName)
    budg_dataG = budg_dataG + budg_dataGtmp
    budg_dataL = 0.0_r8

  end subroutine seq_diag_sum0_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diag_atm_moab - compute global atm input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global atm input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2008-jul-10 - T. Craig - update
  !    2026-Apr-10 - R. Jacob - convert to moab
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diag_atm_moab( atm, infodata, do_a2x, do_x2a)

    ! !INPUT/OUTPUT PARAMETERS:

    type(component_type)    , intent(in)           :: atm    ! component type for instance1
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_a2x
    logical                 , intent(in), optional :: do_x2a

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp            ! for optional field detection
    character(CL)            :: atm_gnam          ! atm grid
    character(CL)            :: lnd_gnam          ! lnd grid
    integer(in)              :: k,n,ic,nf,ip      ! generic index
    integer(in)              :: lSize             ! size of mesh
    real(r8)                 :: ca_a              ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    real(r8), allocatable    :: lat_data(:)
    real(r8), allocatable    :: afrac_data(:)
    real(r8), allocatable    :: lfrac_data(:)
    real(r8), allocatable    :: ofrac_data(:)
    real(r8), allocatable    :: ifrac_data(:)
    real(r8), allocatable    :: fld_swnet(:)
    real(r8), allocatable    :: fld_lwdn(:)
    real(r8), allocatable    :: fld_rainc(:)
    real(r8), allocatable    :: fld_rainl(:)
    real(r8), allocatable    :: fld_snowc(:)
    real(r8), allocatable    :: fld_snowl(:)
    real(r8), allocatable    :: fld_rainc_16O(:)
    real(r8), allocatable    :: fld_rainc_18O(:)
    real(r8), allocatable    :: fld_rainc_HDO(:)
    real(r8), allocatable    :: fld_rainl_16O(:)
    real(r8), allocatable    :: fld_rainl_18O(:)
    real(r8), allocatable    :: fld_rainl_HDO(:)
    real(r8), allocatable    :: fld_snowc_16O(:)
    real(r8), allocatable    :: fld_snowc_18O(:)
    real(r8), allocatable    :: fld_snowc_HDO(:)
    real(r8), allocatable    :: fld_snowl_16O(:)
    real(r8), allocatable    :: fld_snowl_18O(:)
    real(r8), allocatable    :: fld_snowl_HDO(:)
    real(r8), allocatable    :: fld_lwup(:)
    real(r8), allocatable    :: fld_lat(:)
    real(r8), allocatable    :: fld_sen(:)
    real(r8), allocatable    :: fld_evap(:)
    real(r8), allocatable    :: fld_h2otemp(:)
    real(r8), allocatable    :: fld_evap_16O(:)
    real(r8), allocatable    :: fld_evap_18O(:)
    real(r8), allocatable    :: fld_evap_HDO(:)
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_atm = .false.
    logical,save             :: samegrid_al       ! samegrid atm and lnd

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_atm_moab) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lSize = mbGetnCells(mbaxid)
    allocate(area_data(lSize), lat_data(lSize))
    allocate(afrac_data(lSize), lfrac_data(lSize), ofrac_data(lSize), ifrac_data(lSize))

    call mbGetCellTagVals(mbaxid, afldname,  area_data,  lSize)
    call mbGetCellTagVals(mbaxid, latname,   lat_data,   lSize)
    call mbGetCellTagVals(mbaxid, afracname, afrac_data, lSize)
    call mbGetCellTagVals(mbaxid, ofracname, ofrac_data, lSize)
    call mbGetCellTagVals(mbaxid, ifracname, ifrac_data, lSize)

    if (first_time) then
       call seq_infodata_getData(infodata , &
            lnd_gnam=lnd_gnam             , &
            atm_gnam=atm_gnam             )
       samegrid_al = .true.
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.

       av_tmp => component_get_c2x_cx(atm)
       index_a2x_Faxa_rainc_16O = mct_aVect_indexRA(av_tmp,'Faxa_rainc_16O',perrWith='quiet')
       if (index_a2x_Faxa_rainc_16O /= 0) then
          flds_wiso_atm = .true.
          flds_wiso     = .true.
       end if
    end if

    if (samegrid_al) then
       call mbGetCellTagVals(mbaxid, lfracname, lfrac_data, lSize)
    else
       call mbGetCellTagVals(mbaxid, lfrinname, lfrac_data, lSize)
    endif

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    ip = p_inst

    if (present(do_a2x)) then
       allocate(fld_swnet(lSize), fld_lwdn(lSize))
       allocate(fld_rainc(lSize), fld_rainl(lSize), fld_snowc(lSize), fld_snowl(lSize))

       call mbGetCellTagVals(mbaxid, 'Faxa_swnet', fld_swnet, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_lwdn',  fld_lwdn,  lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_rainc', fld_rainc, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_rainl', fld_rainl, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_snowc', fld_snowc, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_snowl', fld_snowl, lSize)

       if (flds_wiso_atm) then
          allocate(fld_rainc_16O(lSize), fld_rainc_18O(lSize), fld_rainc_HDO(lSize))
          allocate(fld_rainl_16O(lSize), fld_rainl_18O(lSize), fld_rainl_HDO(lSize))
          allocate(fld_snowc_16O(lSize), fld_snowc_18O(lSize), fld_snowc_HDO(lSize))
          allocate(fld_snowl_16O(lSize), fld_snowl_18O(lSize), fld_snowl_HDO(lSize))

          call mbGetCellTagVals(mbaxid, 'Faxa_rainc_16O', fld_rainc_16O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_rainc_18O', fld_rainc_18O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_rainc_HDO', fld_rainc_HDO, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_rainl_16O', fld_rainl_16O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_rainl_18O', fld_rainl_18O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_rainl_HDO', fld_rainl_HDO, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_snowc_16O', fld_snowc_16O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_snowc_18O', fld_snowc_18O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_snowc_HDO', fld_snowc_HDO, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_snowl_16O', fld_snowl_16O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_snowl_18O', fld_snowl_18O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxa_snowl_HDO', fld_snowl_HDO, lSize)
       end if

       do n=1,lSize
          do k=1,4

             if (k == 1) then
                ic = c_atm_ar
                ca_a = -area_data(n) * afrac_data(n)
             elseif (k == 2) then
                ic = c_lnd_ar
                ca_a =  area_data(n) * lfrac_data(n)
             elseif (k == 3) then
                ic = c_ocn_ar
                ca_a =  area_data(n) * ofrac_data(n)
             elseif (k == 4) then
                if (lat_data(n) > 0.0_r8) then
                   ic = c_inh_ar
                else
                   ic = c_ish_ar
                endif
                ca_a = area_data(n) * ifrac_data(n)
             endif

             nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a
             nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_swnet(n)
             nf = f_hlwdn ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_lwdn(n)
             nf = f_wrain ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_rainc(n) &
                  + ca_a*fld_rainl(n)
             nf = f_wsnow ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_snowc(n) &
                  + ca_a*fld_snowl(n)
             if (flds_wiso_atm) then
                nf = f_wrain_16O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*fld_rainc_16O(n) + ca_a*fld_rainl_16O(n)
                nf = f_wrain_18O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*fld_rainc_18O(n) + ca_a*fld_rainl_18O(n)
                nf = f_wrain_HDO;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*fld_rainc_HDO(n) + ca_a*fld_rainl_HDO(n)
                nf = f_wsnow_16O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*fld_snowc_16O(n) + ca_a*fld_snowl_16O(n)
                nf = f_wsnow_18O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*fld_snowc_18O(n) + ca_a*fld_snowl_18O(n)
                nf = f_wsnow_HDO;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                     ca_a*fld_snowc_HDO(n) + ca_a*fld_snowl_HDO(n)
             end if
          enddo
       enddo
       ! --- heat implied by snow flux ---
       ic = c_atm_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_lnd_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_ocn_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_inh_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       ic = c_ish_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice

       deallocate(fld_swnet, fld_lwdn)
       deallocate(fld_rainc, fld_rainl, fld_snowc, fld_snowl)
       if (flds_wiso_atm) then
          deallocate(fld_rainc_16O, fld_rainc_18O, fld_rainc_HDO)
          deallocate(fld_rainl_16O, fld_rainl_18O, fld_rainl_HDO)
          deallocate(fld_snowc_16O, fld_snowc_18O, fld_snowc_HDO)
          deallocate(fld_snowl_16O, fld_snowl_18O, fld_snowl_HDO)
       end if
    end if

    if (present(do_x2a)) then
       allocate(fld_lwup(lSize), fld_lat(lSize), fld_sen(lSize), fld_evap(lSize), fld_h2otemp(lSize))

       call mbGetCellTagVals(mbaxid, 'Faxx_lwup',    fld_lwup,    lSize)
       call mbGetCellTagVals(mbaxid, 'Faxx_lat',     fld_lat,     lSize)
       call mbGetCellTagVals(mbaxid, 'Faxx_sen',     fld_sen,     lSize)
       call mbGetCellTagVals(mbaxid, 'Faxx_evap',    fld_evap,    lSize)
       call mbGetCellTagVals(mbaxid, 'Faoo_h2otemp', fld_h2otemp, lSize)

       if (flds_wiso_atm) then
          allocate(fld_evap_16O(lSize), fld_evap_18O(lSize), fld_evap_HDO(lSize))
          call mbGetCellTagVals(mbaxid, 'Faxx_evap_16O', fld_evap_16O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxx_evap_18O', fld_evap_18O, lSize)
          call mbGetCellTagVals(mbaxid, 'Faxx_evap_HDO', fld_evap_HDO, lSize)
       end if

       do n=1,lSize
          do k=1,4

             if (k == 1) then
                ic = c_atm_as
                ca_a = -area_data(n) * afrac_data(n)
             elseif (k == 2) then
                ic = c_lnd_as
                ca_a =  area_data(n) * lfrac_data(n)
             elseif (k == 3) then
                ic = c_ocn_as
                ca_a =  area_data(n) * ofrac_data(n)
             elseif (k == 4) then
                if (lat_data(n) > 0.0_r8) then
                   ic = c_inh_as
                else
                   ic = c_ish_as
                endif
                ca_a = area_data(n) * ifrac_data(n)
             endif

             nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a
             nf = f_hlwup; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_lwup(n)
             nf = f_hlatv; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_lat(n)
             nf = f_hsen ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_sen(n)
             nf = f_hh2ot; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_h2otemp(n)
             nf = f_wevap; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_evap(n)

             if (flds_wiso_atm) then
                nf = f_wevap_16O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_evap_16O(n)
                nf = f_wevap_18O;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_evap_18O(n)
                nf = f_wevap_HDO;
                budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_evap_HDO(n)
             end if

          enddo
       enddo

       deallocate(fld_lwup, fld_lat, fld_sen, fld_evap, fld_h2otemp)
       if (flds_wiso_atm) then
          deallocate(fld_evap_16O, fld_evap_18O, fld_evap_HDO)
       end if
    end if

    deallocate(area_data, lat_data)
    deallocate(afrac_data, lfrac_data, ofrac_data, ifrac_data)

    first_time = .false.

  end subroutine seq_diag_atm_moab

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

  subroutine seq_diag_lnd_moab( lnd, infodata, do_l2x, do_x2l)

    type(component_type)    , intent(in)           :: lnd    ! component type for instance1
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_l2x
    logical                 , intent(in), optional :: do_x2l

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp        ! for optional field detection
    integer(in)              :: n,ic,nf,ip    ! generic index
    integer(in)              :: lSize         ! size of mesh
    real(r8)                 :: ca_l          ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    real(r8), allocatable    :: lfrin_data(:)
    real(r8), allocatable    :: fld_swnet(:)
    real(r8), allocatable    :: fld_lwup(:)
    real(r8), allocatable    :: fld_lat(:)
    real(r8), allocatable    :: fld_sen(:)
    real(r8), allocatable    :: fld_evap(:)
    real(r8), allocatable    :: fld_rofsur(:)
    real(r8), allocatable    :: fld_rofgwl(:)
    real(r8), allocatable    :: fld_rofsub(:)
    real(r8), allocatable    :: fld_rofdto(:)
    real(r8), allocatable    :: fld_wslake(:)
    real(r8), allocatable    :: fld_irrig(:)
    real(r8), allocatable    :: fld_rofi(:)
    real(r8), allocatable    :: fld_evap_16O(:)
    real(r8), allocatable    :: fld_evap_18O(:)
    real(r8), allocatable    :: fld_evap_HDO(:)
    real(r8), allocatable    :: fld_rofl_16O(:)
    real(r8), allocatable    :: fld_rofl_18O(:)
    real(r8), allocatable    :: fld_rofl_HDO(:)
    real(r8), allocatable    :: fld_rofi_16O(:)
    real(r8), allocatable    :: fld_rofi_18O(:)
    real(r8), allocatable    :: fld_rofi_HDO(:)
    real(r8), allocatable    :: fld_lwdn(:)
    real(r8), allocatable    :: fld_rainc(:)
    real(r8), allocatable    :: fld_rainl(:)
    real(r8), allocatable    :: fld_snowc(:)
    real(r8), allocatable    :: fld_snowl(:)
    real(r8), allocatable    :: fld_flood(:)
    real(r8), allocatable    :: fld_supply(:)
    real(r8), allocatable    :: fld_rainc_16O(:)
    real(r8), allocatable    :: fld_rainc_18O(:)
    real(r8), allocatable    :: fld_rainc_HDO(:)
    real(r8), allocatable    :: fld_rainl_16O(:)
    real(r8), allocatable    :: fld_rainl_18O(:)
    real(r8), allocatable    :: fld_rainl_HDO(:)
    real(r8), allocatable    :: fld_snowc_16O(:)
    real(r8), allocatable    :: fld_snowc_18O(:)
    real(r8), allocatable    :: fld_snowc_HDO(:)
    real(r8), allocatable    :: fld_snowl_16O(:)
    real(r8), allocatable    :: fld_snowl_18O(:)
    real(r8), allocatable    :: fld_snowl_HDO(:)
    real(r8), allocatable    :: fld_flood_16O(:)
    real(r8), allocatable    :: fld_flood_18O(:)
    real(r8), allocatable    :: fld_flood_HDO(:)
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_lnd = .false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_lnd_moab) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    lSize = mbGetnCells(mblxid)
    allocate(area_data(lSize), lfrin_data(lSize))

    call mbGetCellTagVals(mblxid, afldname,  area_data,  lSize)
    call mbGetCellTagVals(mblxid, lfrinname, lfrin_data, lSize)

    if (first_time) then
       av_tmp => component_get_c2x_cx(lnd)
       index_l2x_Fall_evap_16O = mct_aVect_indexRA(av_tmp,'Fall_evap_16O',perrWith='quiet')
       if (index_l2x_Fall_evap_16O /= 0) then
          flds_wiso_lnd = .true.
          flds_wiso     = .true.
       end if
       index_l2x_Flrl_irrig = mct_aVect_indexRA(av_tmp,'Flrl_irrig',perrWith='quiet')
    end if

    ip = p_inst

    if (present(do_l2x)) then
       allocate(fld_swnet(lSize), fld_lwup(lSize), fld_lat(lSize), fld_sen(lSize), fld_evap(lSize))
       allocate(fld_rofsur(lSize), fld_rofgwl(lSize), fld_rofsub(lSize), fld_rofdto(lSize), fld_wslake(lSize))
       allocate(fld_rofi(lSize))

       call mbGetCellTagVals(mblxid, 'Fall_swnet',  fld_swnet,  lSize)
       call mbGetCellTagVals(mblxid, 'Fall_lwup',   fld_lwup,   lSize)
       call mbGetCellTagVals(mblxid, 'Fall_lat',    fld_lat,    lSize)
       call mbGetCellTagVals(mblxid, 'Fall_sen',    fld_sen,    lSize)
       call mbGetCellTagVals(mblxid, 'Fall_evap',   fld_evap,   lSize)
       call mbGetCellTagVals(mblxid, 'Flrl_rofsur', fld_rofsur, lSize)
       call mbGetCellTagVals(mblxid, 'Flrl_rofgwl', fld_rofgwl, lSize)
       call mbGetCellTagVals(mblxid, 'Flrl_rofsub', fld_rofsub, lSize)
       call mbGetCellTagVals(mblxid, 'Flrl_rofdto', fld_rofdto, lSize)
       call mbGetCellTagVals(mblxid, 'Flrl_wslake', fld_wslake, lSize)
       call mbGetCellTagVals(mblxid, 'Flrl_rofi',   fld_rofi,   lSize)

       if (index_l2x_Flrl_irrig /= 0) then
          allocate(fld_irrig(lSize))
          call mbGetCellTagVals(mblxid, 'Flrl_irrig', fld_irrig, lSize)
       end if

       if (flds_wiso_lnd) then
          allocate(fld_evap_16O(lSize), fld_evap_18O(lSize), fld_evap_HDO(lSize))
          allocate(fld_rofl_16O(lSize), fld_rofl_18O(lSize), fld_rofl_HDO(lSize))
          allocate(fld_rofi_16O(lSize), fld_rofi_18O(lSize), fld_rofi_HDO(lSize))
          call mbGetCellTagVals(mblxid, 'Fall_evap_16O', fld_evap_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Fall_evap_18O', fld_evap_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Fall_evap_HDO', fld_evap_HDO, lSize)
          call mbGetCellTagVals(mblxid, 'Flrl_rofl_16O', fld_rofl_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Flrl_rofl_18O', fld_rofl_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Flrl_rofl_HDO', fld_rofl_HDO, lSize)
          call mbGetCellTagVals(mblxid, 'Flrl_rofi_16O', fld_rofi_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Flrl_rofi_18O', fld_rofi_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Flrl_rofi_HDO', fld_rofi_HDO, lSize)
       end if

       ic = c_lnd_lr
       do n=1,lSize
          ca_l = area_data(n) * lfrin_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_swnet(n)
          nf = f_hlwup ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_lwup(n)
          nf = f_hlatv ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_lat(n)
          nf = f_hsen  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_sen(n)
          nf = f_wevap ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_evap(n)
          nf = f_wroff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofsur(n) &
               - ca_l*fld_rofgwl(n) &
               - ca_l*fld_rofsub(n) &
               - ca_l*fld_rofdto(n) &
               - ca_l*fld_wslake(n)
          if (index_l2x_Flrl_irrig /= 0) then
             nf = f_wroff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_irrig(n)
          end if
          nf = f_wioff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofi(n)

          if (flds_wiso_lnd) then
             nf = f_wevap_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_evap_16O(n)
             nf = f_wevap_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_evap_18O(n)
             nf = f_wevap_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_evap_HDO(n)

             nf = f_wroff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofl_16O(n)
             nf = f_wroff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofl_18O(n)
             nf = f_wroff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofl_HDO(n)

             nf = f_wioff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofi_16O(n)
             nf = f_wioff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofi_18O(n)
             nf = f_wioff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_rofi_HDO(n)
          end if
       end do
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

       deallocate(fld_swnet, fld_lwup, fld_lat, fld_sen, fld_evap)
       deallocate(fld_rofsur, fld_rofgwl, fld_rofsub, fld_rofdto, fld_wslake)
       deallocate(fld_rofi)
       if (index_l2x_Flrl_irrig /= 0) deallocate(fld_irrig)
       if (flds_wiso_lnd) then
          deallocate(fld_evap_16O, fld_evap_18O, fld_evap_HDO)
          deallocate(fld_rofl_16O, fld_rofl_18O, fld_rofl_HDO)
          deallocate(fld_rofi_16O, fld_rofi_18O, fld_rofi_HDO)
       end if
    end if

    if (present(do_x2l)) then
       allocate(fld_lwdn(lSize), fld_rainc(lSize), fld_rainl(lSize), fld_snowc(lSize), fld_snowl(lSize))
       allocate(fld_flood(lSize), fld_supply(lSize))

       call mbGetCellTagVals(mblxid, 'Faxa_lwdn',   fld_lwdn,   lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_rainc',  fld_rainc,  lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_rainl',  fld_rainl,  lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_snowc',  fld_snowc,  lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_snowl',  fld_snowl,  lSize)
       call mbGetCellTagVals(mblxid, 'Flrr_flood',  fld_flood,  lSize)
       call mbGetCellTagVals(mblxid, 'Flrr_supply', fld_supply, lSize)

       if (flds_wiso_lnd) then
          allocate(fld_rainc_16O(lSize), fld_rainc_18O(lSize), fld_rainc_HDO(lSize))
          allocate(fld_rainl_16O(lSize), fld_rainl_18O(lSize), fld_rainl_HDO(lSize))
          allocate(fld_snowc_16O(lSize), fld_snowc_18O(lSize), fld_snowc_HDO(lSize))
          allocate(fld_snowl_16O(lSize), fld_snowl_18O(lSize), fld_snowl_HDO(lSize))
          allocate(fld_flood_16O(lSize), fld_flood_18O(lSize), fld_flood_HDO(lSize))
          call mbGetCellTagVals(mblxid, 'Faxa_rainc_16O', fld_rainc_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_rainc_18O', fld_rainc_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_rainc_HDO', fld_rainc_HDO, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_rainl_16O', fld_rainl_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_rainl_18O', fld_rainl_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_rainl_HDO', fld_rainl_HDO, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_snowc_16O', fld_snowc_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_snowc_18O', fld_snowc_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_snowc_HDO', fld_snowc_HDO, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_snowl_16O', fld_snowl_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_snowl_18O', fld_snowl_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Faxa_snowl_HDO', fld_snowl_HDO, lSize)
          call mbGetCellTagVals(mblxid, 'Flrr_flood_16O', fld_flood_16O, lSize)
          call mbGetCellTagVals(mblxid, 'Flrr_flood_18O', fld_flood_18O, lSize)
          call mbGetCellTagVals(mblxid, 'Flrr_flood_HDO', fld_flood_HDO, lSize)
       end if

       ic = c_lnd_ls
       do n=1,lSize
          ca_l = area_data(n) * lfrin_data(n)
          nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          nf = f_hlwdn; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_lwdn(n)
          nf = f_wrain; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_rainc(n) &
               + ca_l*fld_rainl(n)
          nf = f_wsnow; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_snowc(n) &
               + ca_l*fld_snowl(n)
          nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_flood(n)
          nf = f_wirrig ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_supply(n)

          if (flds_wiso_lnd) then
             nf = f_wrain_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*fld_rainc_16O(n) + ca_l*fld_rainl_16O(n)
             nf = f_wrain_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*fld_rainc_18O(n) + ca_l*fld_rainl_18O(n)
             nf = f_wrain_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*fld_rainc_HDO(n) + ca_l*fld_rainl_HDO(n)

             nf = f_wsnow_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*fld_snowc_16O(n) + ca_l*fld_snowl_16O(n)
             nf = f_wsnow_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*fld_snowc_18O(n) + ca_l*fld_snowl_18O(n)
             nf = f_wsnow_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
                  ca_l*fld_snowc_HDO(n) + ca_l*fld_snowl_HDO(n)

             nf = f_wroff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_flood_16O(n)
             nf = f_wroff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_flood_18O(n)
             nf = f_wroff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_l*fld_flood_HDO(n)
          end if
       end do
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice

       deallocate(fld_lwdn, fld_rainc, fld_rainl, fld_snowc, fld_snowl)
       deallocate(fld_flood, fld_supply)
       if (flds_wiso_lnd) then
          deallocate(fld_rainc_16O, fld_rainc_18O, fld_rainc_HDO)
          deallocate(fld_rainl_16O, fld_rainl_18O, fld_rainl_HDO)
          deallocate(fld_snowc_16O, fld_snowc_18O, fld_snowc_HDO)
          deallocate(fld_snowl_16O, fld_snowl_18O, fld_snowl_HDO)
          deallocate(fld_flood_16O, fld_flood_18O, fld_flood_HDO)
       end if
    end if

    deallocate(area_data, lfrin_data)

    first_time = .false.

  end subroutine seq_diag_lnd_moab

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

  subroutine seq_diag_rof_moab( rof, infodata)

    type(component_type)    , intent(in) :: rof    ! component type for instance1
    type(seq_infodata_type) , intent(in) :: infodata

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp        ! for optional field detection
    integer(in)              :: n,ic,nf,ip    ! generic index
    integer(in)              :: lSize         ! size of mesh
    real(r8)                 :: ca_r          ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    real(r8), allocatable    :: fld_rofsur(:), fld_rofgwl(:), fld_rofsub(:), fld_rofdto(:), fld_rofi(:)
    real(r8), allocatable    :: fld_irrig(:)
    real(r8), allocatable    :: fld_rofl_16O(:), fld_rofl_18O(:), fld_rofl_HDO(:)
    real(r8), allocatable    :: fld_x2r_rofi_16O(:), fld_x2r_rofi_18O(:), fld_x2r_rofi_HDO(:)
    real(r8), allocatable    :: fld_r2x_rofl(:), fld_r2x_rofi(:), fld_firr_rofi(:)
    real(r8), allocatable    :: fld_flood(:), fld_supply(:)
    real(r8), allocatable    :: fld_r2x_rofl_16O(:), fld_r2x_rofl_18O(:), fld_r2x_rofl_HDO(:)
    real(r8), allocatable    :: fld_r2x_rofi_16O(:), fld_r2x_rofi_18O(:), fld_r2x_rofi_HDO(:)
    real(r8), allocatable    :: fld_flood_16O(:),    fld_flood_18O(:),    fld_flood_HDO(:)
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_rof = .false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_rof_moab) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    ip = p_inst

    lSize = mbGetnCells(mbrxid)
    allocate(area_data(lSize))
    call mbGetCellTagVals(mbrxid, afldname, area_data, lSize)

    if (first_time) then
       av_tmp => component_get_x2c_cx(rof)
       index_x2r_Flrl_irrig    = mct_aVect_indexRA(av_tmp,'Flrl_irrig',   perrWith='quiet')
       index_x2r_Flrl_rofl_16O = mct_aVect_indexRA(av_tmp,'Flrl_rofl_16O',perrWith='quiet')
       if ( index_x2r_Flrl_rofl_16O /= 0 ) then
          flds_wiso_rof = .true.
          flds_wiso     = .true.
       end if
    end if

    ! x2r block (c_rof_rr)
    ic = c_rof_rr
    allocate(fld_rofsur(lSize), fld_rofgwl(lSize), fld_rofsub(lSize), fld_rofdto(lSize), fld_rofi(lSize))
    call mbGetCellTagVals(mbrxid, 'Flrl_rofsur', fld_rofsur, lSize)
    call mbGetCellTagVals(mbrxid, 'Flrl_rofgwl', fld_rofgwl, lSize)
    call mbGetCellTagVals(mbrxid, 'Flrl_rofsub', fld_rofsub, lSize)
    call mbGetCellTagVals(mbrxid, 'Flrl_rofdto', fld_rofdto, lSize)
    call mbGetCellTagVals(mbrxid, 'Flrl_rofi',   fld_rofi,   lSize)
    if (index_x2r_Flrl_irrig /= 0) then
       allocate(fld_irrig(lSize))
       call mbGetCellTagVals(mbrxid, 'Flrl_irrig', fld_irrig, lSize)
    end if
    if ( flds_wiso_rof ) then
       allocate(fld_rofl_16O(lSize), fld_rofl_18O(lSize), fld_rofl_HDO(lSize))
       allocate(fld_x2r_rofi_16O(lSize), fld_x2r_rofi_18O(lSize), fld_x2r_rofi_HDO(lSize))
       call mbGetCellTagVals(mbrxid, 'Flrl_rofl_16O', fld_rofl_16O,      lSize)
       call mbGetCellTagVals(mbrxid, 'Flrl_rofl_18O', fld_rofl_18O,      lSize)
       call mbGetCellTagVals(mbrxid, 'Flrl_rofl_HDO', fld_rofl_HDO,      lSize)
       call mbGetCellTagVals(mbrxid, 'Flrl_rofi_16O', fld_x2r_rofi_16O,  lSize)
       call mbGetCellTagVals(mbrxid, 'Flrl_rofi_18O', fld_x2r_rofi_18O,  lSize)
       call mbGetCellTagVals(mbrxid, 'Flrl_rofi_HDO', fld_x2r_rofi_HDO,  lSize)
    end if
    do n=1,lSize
       ca_r = area_data(n)
       nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_rofsur(n) &
            + ca_r*fld_rofgwl(n) &
            + ca_r*fld_rofsub(n) &
            + ca_r*fld_rofdto(n)
       if (index_x2r_Flrl_irrig /= 0) then
          nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_irrig(n)
       end if
       nf = f_wioff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_rofi(n)
       if ( flds_wiso_rof ) then
          nf = f_wroff_16O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_rofl_16O(n)
          nf = f_wroff_18O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_rofl_18O(n)
          nf = f_wroff_HDO; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_rofl_HDO(n)
          nf = f_wioff_16O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_x2r_rofi_16O(n)
          nf = f_wioff_18O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_x2r_rofi_18O(n)
          nf = f_wioff_HDO; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_x2r_rofi_HDO(n)
       end if
    end do
    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
    deallocate(fld_rofsur, fld_rofgwl, fld_rofsub, fld_rofdto, fld_rofi)
    if (index_x2r_Flrl_irrig /= 0) deallocate(fld_irrig)
    if ( flds_wiso_rof ) then
       deallocate(fld_rofl_16O, fld_rofl_18O, fld_rofl_HDO)
       deallocate(fld_x2r_rofi_16O, fld_x2r_rofi_18O, fld_x2r_rofi_HDO)
    end if

    ! r2x block (c_rof_rs)
    ic = c_rof_rs
    allocate(fld_r2x_rofl(lSize), fld_r2x_rofi(lSize), fld_firr_rofi(lSize))
    allocate(fld_flood(lSize), fld_supply(lSize))
    call mbGetCellTagVals(mbrxid, 'Forr_rofl',   fld_r2x_rofl, lSize)
    call mbGetCellTagVals(mbrxid, 'Forr_rofi',   fld_r2x_rofi, lSize)
    call mbGetCellTagVals(mbrxid, 'Firr_rofi',   fld_firr_rofi, lSize)
    call mbGetCellTagVals(mbrxid, 'Flrr_flood',  fld_flood,    lSize)
    call mbGetCellTagVals(mbrxid, 'Flrr_supply', fld_supply,   lSize)
    if ( flds_wiso_rof ) then
       allocate(fld_r2x_rofl_16O(lSize), fld_r2x_rofl_18O(lSize), fld_r2x_rofl_HDO(lSize))
       allocate(fld_r2x_rofi_16O(lSize), fld_r2x_rofi_18O(lSize), fld_r2x_rofi_HDO(lSize))
       allocate(fld_flood_16O(lSize),    fld_flood_18O(lSize),    fld_flood_HDO(lSize))
       call mbGetCellTagVals(mbrxid, 'Forr_rofl_16O',  fld_r2x_rofl_16O, lSize)
       call mbGetCellTagVals(mbrxid, 'Forr_rofl_18O',  fld_r2x_rofl_18O, lSize)
       call mbGetCellTagVals(mbrxid, 'Forr_rofl_HDO',  fld_r2x_rofl_HDO, lSize)
       call mbGetCellTagVals(mbrxid, 'Forr_rofi_16O',  fld_r2x_rofi_16O, lSize)
       call mbGetCellTagVals(mbrxid, 'Forr_rofi_18O',  fld_r2x_rofi_18O, lSize)
       call mbGetCellTagVals(mbrxid, 'Forr_rofi_HDO',  fld_r2x_rofi_HDO, lSize)
       call mbGetCellTagVals(mbrxid, 'Flrr_flood_16O', fld_flood_16O,    lSize)
       call mbGetCellTagVals(mbrxid, 'Flrr_flood_18O', fld_flood_18O,    lSize)
       call mbGetCellTagVals(mbrxid, 'Flrr_flood_HDO', fld_flood_HDO,    lSize)
    end if
    do n=1,lSize
       ca_r = area_data(n)
       nf = f_wroff;  budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofl(n) &
            + ca_r*fld_flood(n)
       nf = f_wioff;  budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofi(n) &
            - ca_r*fld_firr_rofi(n)
       nf = f_wirrig; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_supply(n)
       if ( flds_wiso_rof ) then
          nf = f_wroff_16O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofl_16O(n)
          nf = f_wroff_18O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofl_18O(n)
          nf = f_wroff_HDO; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofl_HDO(n)
          nf = f_wioff_16O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofi_16O(n)
          nf = f_wioff_18O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofi_18O(n)
          nf = f_wioff_HDO; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_r*fld_r2x_rofi_HDO(n)
          nf = f_wroff_16O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_flood_16O(n)
          nf = f_wroff_18O; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_flood_18O(n)
          nf = f_wroff_HDO; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_r*fld_flood_HDO(n)
       end if
    end do
    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
    deallocate(fld_r2x_rofl, fld_r2x_rofi, fld_firr_rofi, fld_flood, fld_supply)
    if ( flds_wiso_rof ) then
       deallocate(fld_r2x_rofl_16O, fld_r2x_rofl_18O, fld_r2x_rofl_HDO)
       deallocate(fld_r2x_rofi_16O, fld_r2x_rofi_18O, fld_r2x_rofi_HDO)
       deallocate(fld_flood_16O,    fld_flood_18O,    fld_flood_HDO)
    end if

    deallocate(area_data)
    first_time = .false.

  end subroutine seq_diag_rof_moab

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

  subroutine seq_diag_glc_moab( glc, infodata)

    type(component_type)    , intent(in) :: glc    ! component type for instance1
    type(seq_infodata_type) , intent(in) :: infodata

    !EOP

    !----- local -----
    !TODO  change this to seq_comm_mct when GLC is ported
    integer(in)             :: mbgxid
    integer(in)              :: n,ic,nf,ip    ! generic index
    integer(in)              :: lSize         ! size of mesh
    real(r8)                 :: ca_g          ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    real(r8), allocatable    :: fld_rofl(:), fld_rofi(:), fld_irrofi(:)

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_glc_moab) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    ip = p_inst
    ic = c_glc_gs

    lSize = mbGetnCells(mbgxid)
    allocate(area_data(lSize), fld_rofl(lSize), fld_rofi(lSize), fld_irrofi(lSize))
    call mbGetCellTagVals(mbgxid, afldname,    area_data, lSize)
    call mbGetCellTagVals(mbgxid, 'Fogg_rofl', fld_rofl,  lSize)
    call mbGetCellTagVals(mbgxid, 'Fogg_rofi', fld_rofi,  lSize)
    call mbGetCellTagVals(mbgxid, 'Figg_rofi', fld_irrofi,lSize)
    do n=1,lSize
       ca_g = area_data(n)
       nf = f_wroff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_g*fld_rofl(n)
       nf = f_wioff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_g*fld_rofi(n) &
            - ca_g*fld_irrofi(n)
    end do
    budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
    deallocate(area_data, fld_rofl, fld_rofi, fld_irrofi)

  end subroutine seq_diag_glc_moab

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

  subroutine seq_diag_ocn_moab( ocn, infodata, do_o2x, do_x2o, do_xao)

    type(component_type)    , intent(in)          :: ocn    ! component type for instance1
    type(seq_infodata_type) , intent(in)          :: infodata
    logical                 , intent(in),optional :: do_o2x
    logical                 , intent(in),optional :: do_x2o
    logical                 , intent(in),optional :: do_xao

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp        ! for optional field detection
    integer(in)              :: n,nf,ic,ip    ! generic index
    integer(in)              :: lSize         ! size of mesh
    real(r8)                 :: ca_i,ca_o,ca_c  ! area of a grid cell
    real(r8), allocatable    :: area_data(:), ofrac_data(:), ifrac_data(:)
    real(r8), allocatable    :: fld_frazil(:), fld_q(:), fld_h2otemp(:)
    real(r8), allocatable    :: fld_frazil_li(:), fld_q_li(:)
    real(r8), allocatable    :: fld_ismw(:), fld_rrofl(:), fld_rrofi(:), fld_ismh(:), fld_rrofih(:)
    real(r8), allocatable    :: fld_lwup(:), fld_lat(:), fld_sen(:), fld_evap(:)
    real(r8), allocatable    :: fld_evap_16O(:), fld_evap_18O(:), fld_evap_HDO(:)
    real(r8), allocatable    :: fld_melth(:), fld_meltw(:), fld_bergh(:), fld_bergw(:)
    real(r8), allocatable    :: fld_swnet(:), fld_lwdn(:), fld_rain(:), fld_snow(:)
    real(r8), allocatable    :: fld_rofl(:), fld_rofi(:)
    real(r8), allocatable    :: fld_meltw_16O(:), fld_meltw_18O(:), fld_meltw_HDO(:)
    real(r8), allocatable    :: fld_rain_16O(:),  fld_rain_18O(:),  fld_rain_HDO(:)
    real(r8), allocatable    :: fld_snow_16O(:),  fld_snow_18O(:),  fld_snow_HDO(:)
    real(r8), allocatable    :: fld_rofl_16O(:),  fld_rofi_16O(:)
    real(r8), allocatable    :: fld_rofl_18O(:),  fld_rofi_18O(:)
    real(r8), allocatable    :: fld_rofl_HDO(:),  fld_rofi_HDO(:)
    logical,save             :: first_time    = .true.
    logical,save             :: flds_wiso_ocn = .false.
    logical,save             :: flds_polar    = .false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_ocn_moab) '

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

    ip = p_inst

    lSize = mbGetnCells(mboxid)
    allocate(area_data(lSize), ofrac_data(lSize), ifrac_data(lSize))
    call mbGetCellTagVals(mboxid, afldname,  area_data,  lSize)
    call mbGetCellTagVals(mboxid, ofracname, ofrac_data, lSize)
    call mbGetCellTagVals(mboxid, ifracname, ifrac_data, lSize)

    if (first_time) then
       av_tmp => component_get_c2x_cx(ocn)
       index_o2x_Foxo_ismw = mct_aVect_indexRA(av_tmp,'Foxo_ismw',perrWith='quiet')
       if ( index_o2x_Foxo_ismw /= 0 ) flds_polar = .true.
       av_tmp => component_get_x2c_cx(ocn)
       index_x2o_Faxa_rain_16O = mct_aVect_indexRA(av_tmp,'Faxa_rain_16O',perrWith='quiet')
       if ( index_x2o_Faxa_rain_16O /= 0 ) then
          flds_wiso_ocn = .true.
          flds_wiso     = .true.
       end if
    end if

    if (present(do_o2x)) then
       allocate(fld_frazil(lSize), fld_q(lSize), fld_h2otemp(lSize))
       call mbGetCellTagVals(mboxid, 'Fioo_frazil',  fld_frazil,  lSize)
       call mbGetCellTagVals(mboxid, 'Fioo_q',       fld_q,       lSize)
       call mbGetCellTagVals(mboxid, 'Faoo_h2otemp', fld_h2otemp, lSize)
       if (flds_polar) then
          allocate(fld_frazil_li(lSize), fld_q_li(lSize))
          allocate(fld_ismw(lSize), fld_rrofl(lSize), fld_rrofi(lSize), fld_ismh(lSize), fld_rrofih(lSize))
          call mbGetCellTagVals(mboxid, 'Foxo_frazil_li', fld_frazil_li, lSize)
          call mbGetCellTagVals(mboxid, 'Foxo_q_li',      fld_q_li,      lSize)
          call mbGetCellTagVals(mboxid, 'Foxo_ismw',      fld_ismw,      lSize)
          call mbGetCellTagVals(mboxid, 'Foxo_rrofl',     fld_rrofl,     lSize)
          call mbGetCellTagVals(mboxid, 'Foxo_rrofi',     fld_rrofi,     lSize)
          call mbGetCellTagVals(mboxid, 'Foxo_ismh',      fld_ismh,      lSize)
          call mbGetCellTagVals(mboxid, 'Foxo_rrofih',    fld_rrofih,    lSize)
       end if
       ic = c_ocn_or
       do n=1,lSize
          ca_o =  area_data(n) * ofrac_data(n)
          ca_i =  area_data(n) * ifrac_data(n)
          ca_c =  area_data(n)
          nf = f_area;  budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o
          nf = f_wfrz;  budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - (ca_o+ca_i)*max(0.0_r8,fld_frazil(n))
          nf = f_hfrz;  budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*max(0.0_r8,fld_q(n))
          nf = f_hh2ot; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_h2otemp(n)
          if (flds_polar) then
             nf = f_wpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_c*fld_frazil_li(n)
             nf = f_wpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_c*fld_ismw(n)
             nf = f_wpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - (ca_o+ca_i)*fld_rrofl(n)
             nf = f_wpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - (ca_o+ca_i)*fld_rrofi(n)
             nf = f_hpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_c*fld_q_li(n)
             nf = f_hpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_c*fld_ismh(n)
             nf = f_hpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_c*fld_rrofih(n)
          end if
       end do
       deallocate(fld_frazil, fld_q, fld_h2otemp)
       if (flds_polar) deallocate(fld_frazil_li, fld_q_li, fld_ismw, fld_rrofl, fld_rrofi, fld_ismh, fld_rrofih)
    end if

    if (present(do_xao)) then
       allocate(fld_lwup(lSize), fld_lat(lSize), fld_sen(lSize), fld_evap(lSize))
       call mbGetCellTagVals(mbofxid, 'Faox_lwup', fld_lwup, lSize)
       call mbGetCellTagVals(mbofxid, 'Faox_lat',  fld_lat,  lSize)
       call mbGetCellTagVals(mbofxid, 'Faox_sen',  fld_sen,  lSize)
       call mbGetCellTagVals(mbofxid, 'Faox_evap', fld_evap, lSize)
       if ( flds_wiso_ocn ) then
          allocate(fld_evap_16O(lSize), fld_evap_18O(lSize), fld_evap_HDO(lSize))
          call mbGetCellTagVals(mbofxid, 'Faox_evap_16O', fld_evap_16O, lSize)
          call mbGetCellTagVals(mbofxid, 'Faox_evap_18O', fld_evap_18O, lSize)
          call mbGetCellTagVals(mbofxid, 'Faox_evap_HDO', fld_evap_HDO, lSize)
       end if
       ic = c_ocn_or
       do n=1,lSize
          ca_o =  area_data(n) * ofrac_data(n)
          nf = f_hlwup; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*fld_lwup(n)
          nf = f_hlatv; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*fld_lat(n)
          nf = f_hsen ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*fld_sen(n)
          nf = f_wevap; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*fld_evap(n)
          if ( flds_wiso_ocn ) then
             nf = f_wevap_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*fld_evap_16O(n)
             nf = f_wevap_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*fld_evap_18O(n)
             nf = f_wevap_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o*fld_evap_HDO(n)
          end if
       end do
       deallocate(fld_lwup, fld_lat, fld_sen, fld_evap)
       if (flds_wiso_ocn) deallocate(fld_evap_16O, fld_evap_18O, fld_evap_HDO)
    end if

    if (present(do_x2o)) then
       if (.not. present(do_xao)) then
          ! these are in x2o but they really are the atm/ocean flux
          ! computed in the coupler and are "like" an o2x
          allocate(fld_lwup(lSize), fld_lat(lSize), fld_sen(lSize), fld_evap(lSize))
          call mbGetCellTagVals(mboxid, 'Foxx_lwup', fld_lwup, lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_lat',  fld_lat,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_sen',  fld_sen,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_evap', fld_evap, lSize)
          ic = c_ocn_or
          do n=1,lSize
             ca_o =  area_data(n) * ofrac_data(n)
             ca_i =  area_data(n) * ifrac_data(n)
             nf = f_hlwup; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_lwup(n)
             nf = f_hlatv; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_lat(n)
             nf = f_hsen ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_sen(n)
             nf = f_wevap; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_evap(n)
          end do
          deallocate(fld_lwup, fld_lat, fld_sen, fld_evap)
       endif

       allocate(fld_melth(lSize), fld_meltw(lSize), fld_bergh(lSize), fld_bergw(lSize))
       allocate(fld_swnet(lSize), fld_lwdn(lSize), fld_rain(lSize), fld_snow(lSize))
       allocate(fld_rofl(lSize), fld_rofi(lSize))
       call mbGetCellTagVals(mboxid, 'Fioi_melth',  fld_melth,  lSize)
       call mbGetCellTagVals(mboxid, 'Fioi_meltw',  fld_meltw,  lSize)
       call mbGetCellTagVals(mboxid, 'PFioi_bergh', fld_bergh,  lSize)
       call mbGetCellTagVals(mboxid, 'PFioi_bergw', fld_bergw,  lSize)
       call mbGetCellTagVals(mboxid, 'Foxx_swnet',  fld_swnet,  lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_lwdn',   fld_lwdn,   lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_rain',   fld_rain,   lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_snow',   fld_snow,   lSize)
       call mbGetCellTagVals(mboxid, 'Foxx_rofl',   fld_rofl,   lSize)
       call mbGetCellTagVals(mboxid, 'Foxx_rofi',   fld_rofi,   lSize)
       if ( flds_wiso_ocn ) then
          allocate(fld_meltw_16O(lSize), fld_meltw_18O(lSize), fld_meltw_HDO(lSize))
          allocate(fld_rain_16O(lSize),  fld_rain_18O(lSize),  fld_rain_HDO(lSize))
          allocate(fld_snow_16O(lSize),  fld_snow_18O(lSize),  fld_snow_HDO(lSize))
          allocate(fld_rofl_16O(lSize),  fld_rofi_16O(lSize))
          allocate(fld_rofl_18O(lSize),  fld_rofi_18O(lSize))
          allocate(fld_rofl_HDO(lSize),  fld_rofi_HDO(lSize))
          call mbGetCellTagVals(mboxid, 'Fioi_meltw_16O', fld_meltw_16O, lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_meltw_18O', fld_meltw_18O, lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_meltw_HDO', fld_meltw_HDO, lSize)
          call mbGetCellTagVals(mboxid, 'Faxa_rain_16O',  fld_rain_16O,  lSize)
          call mbGetCellTagVals(mboxid, 'Faxa_rain_18O',  fld_rain_18O,  lSize)
          call mbGetCellTagVals(mboxid, 'Faxa_rain_HDO',  fld_rain_HDO,  lSize)
          call mbGetCellTagVals(mboxid, 'Faxa_snow_16O',  fld_snow_16O,  lSize)
          call mbGetCellTagVals(mboxid, 'Faxa_snow_18O',  fld_snow_18O,  lSize)
          call mbGetCellTagVals(mboxid, 'Faxa_snow_HDO',  fld_snow_HDO,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_rofl_16O',  fld_rofl_16O,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_rofi_16O',  fld_rofi_16O,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_rofl_18O',  fld_rofl_18O,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_rofi_18O',  fld_rofi_18O,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_rofl_HDO',  fld_rofl_HDO,  lSize)
          call mbGetCellTagVals(mboxid, 'Foxx_rofi_HDO',  fld_rofi_HDO,  lSize)
       end if
       ic = c_ocn_os
       do n=1,lSize
          ca_o =  area_data(n) * ofrac_data(n)
          ca_i =  area_data(n) * ifrac_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o
          nf = f_hmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_melth(n)
          nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_swnet(n)
          nf = f_hlwdn ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_lwdn(n)
          nf = f_hpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_bergh(n)
          nf = f_wmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_meltw(n)
          nf = f_wrain ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rain(n)
          nf = f_wsnow ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_snow(n)
          nf = f_wpolar; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_bergw(n)
          nf = f_wroff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofl(n)
          nf = f_wioff ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofi(n)
          if ( flds_wiso_ocn ) then
             nf = f_wmelt_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_meltw_16O(n)
             nf = f_wmelt_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_meltw_18O(n)
             nf = f_wmelt_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_meltw_HDO(n)
             nf = f_wrain_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rain_16O(n)
             nf = f_wrain_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rain_18O(n)
             nf = f_wrain_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rain_HDO(n)
             nf = f_wsnow_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_snow_16O(n)
             nf = f_wsnow_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_snow_18O(n)
             nf = f_wsnow_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_snow_HDO(n)
             nf = f_wroff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofl_16O(n)
             nf = f_wioff_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofi_16O(n)
             nf = f_wroff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofl_18O(n)
             nf = f_wioff_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofi_18O(n)
             nf = f_wroff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofl_HDO(n)
             nf = f_wioff_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_rofi_HDO(n)
          end if
       end do
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
       deallocate(fld_melth, fld_meltw, fld_bergh, fld_bergw)
       deallocate(fld_swnet, fld_lwdn, fld_rain, fld_snow)
       deallocate(fld_rofl, fld_rofi)
       if (flds_wiso_ocn) then
          deallocate(fld_meltw_16O, fld_meltw_18O, fld_meltw_HDO)
          deallocate(fld_rain_16O,  fld_rain_18O,  fld_rain_HDO)
          deallocate(fld_snow_16O,  fld_snow_18O,  fld_snow_HDO)
          deallocate(fld_rofl_16O,  fld_rofi_16O)
          deallocate(fld_rofl_18O,  fld_rofi_18O)
          deallocate(fld_rofl_HDO,  fld_rofi_HDO)
       end if
    end if

    ! EBK -- isotope r2x_Forr_rofl/i?

    deallocate(area_data, ofrac_data, ifrac_data)
    first_time = .false.

  end subroutine seq_diag_ocn_moab

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

  subroutine seq_diag_ice_moab( ice, infodata, do_i2x, do_x2i)

    type(component_type)    , intent(in)           :: ice    ! component type for instance1
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_i2x
    logical                 , intent(in), optional :: do_x2i

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp        ! for optional field detection
    integer(in)              :: n,ic,nf,ip    ! generic index
    integer(in)              :: lSize         ! size of mesh
    real(r8)                 :: ca_i,ca_o     ! area of a grid cell
    real(r8), allocatable    :: area_data(:), lat_data(:)
    real(r8), allocatable    :: ofrac_data(:), ifrac_data(:)
    real(r8), allocatable    :: fld_melth(:), fld_meltw(:), fld_swpen(:), fld_swnet(:)
    real(r8), allocatable    :: fld_lwup(:), fld_lat(:), fld_sen(:), fld_evap(:)
    real(r8), allocatable    :: fld_meltw_16O(:), fld_meltw_18O(:), fld_meltw_HDO(:)
    real(r8), allocatable    :: fld_evap_16O(:),  fld_evap_18O(:),  fld_evap_HDO(:)
    real(r8), allocatable    :: fld_lwdn(:), fld_rain(:), fld_snow(:), fld_frazil(:), fld_q(:), fld_rofi(:)
    real(r8), allocatable    :: fld_rain_16O(:), fld_rain_18O(:), fld_rain_HDO(:)
    real(r8), allocatable    :: fld_snow_16O(:), fld_snow_18O(:), fld_snow_HDO(:)
    logical,save             :: first_time        = .true.
    logical,save             :: flds_wiso_ice     = .false.
    logical,save             :: flds_wiso_ice_x2i = .false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_ice_moab) '

    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    ip = p_inst

    lSize = mbGetnCells(mbixid)
    allocate(area_data(lSize), lat_data(lSize), ofrac_data(lSize), ifrac_data(lSize))
    call mbGetCellTagVals(mbixid, afldname,  area_data,  lSize)
    call mbGetCellTagVals(mbixid, latname,   lat_data,   lSize)
    call mbGetCellTagVals(mbixid, ofracname, ofrac_data, lSize)
    call mbGetCellTagVals(mbixid, ifracname, ifrac_data, lSize)

    if (first_time) then
       av_tmp => component_get_c2x_cx(ice)
       index_i2x_Fioi_meltw_16O = mct_aVect_indexRA(av_tmp,'Fioi_meltw_16O',perrWith='quiet')
       if ( index_i2x_Fioi_meltw_16O /= 0 ) then
          flds_wiso_ice = .true.
          flds_wiso     = .true.
       end if
       av_tmp => component_get_x2c_cx(ice)
       index_x2i_Faxa_rain_16O = mct_aVect_indexRA(av_tmp,'Faxa_rain_16O',perrWith='quiet')
       if ( index_x2i_Faxa_rain_16O /= 0 ) then
          flds_wiso_ice_x2i = .true.
          flds_wiso         = .true.
       end if
    end if

    if (present(do_i2x)) then
       allocate(fld_melth(lSize), fld_meltw(lSize), fld_swpen(lSize), fld_swnet(lSize))
       allocate(fld_lwup(lSize), fld_lat(lSize), fld_sen(lSize), fld_evap(lSize))
       call mbGetCellTagVals(mbixid, 'Fioi_melth', fld_melth, lSize)
       call mbGetCellTagVals(mbixid, 'Fioi_meltw', fld_meltw, lSize)
       call mbGetCellTagVals(mbixid, 'Fioi_swpen', fld_swpen, lSize)
       call mbGetCellTagVals(mbixid, 'Faii_swnet', fld_swnet, lSize)
       call mbGetCellTagVals(mbixid, 'Faii_lwup',  fld_lwup,  lSize)
       call mbGetCellTagVals(mbixid, 'Faii_lat',   fld_lat,   lSize)
       call mbGetCellTagVals(mbixid, 'Faii_sen',   fld_sen,   lSize)
       call mbGetCellTagVals(mbixid, 'Faii_evap',  fld_evap,  lSize)
       if ( flds_wiso_ice ) then
          allocate(fld_meltw_16O(lSize), fld_meltw_18O(lSize), fld_meltw_HDO(lSize))
          allocate(fld_evap_16O(lSize),  fld_evap_18O(lSize),  fld_evap_HDO(lSize))
          call mbGetCellTagVals(mbixid, 'Fioi_meltw_16O', fld_meltw_16O, lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_meltw_18O', fld_meltw_18O, lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_meltw_HDO', fld_meltw_HDO, lSize)
          call mbGetCellTagVals(mbixid, 'Faii_evap_16O',  fld_evap_16O,  lSize)
          call mbGetCellTagVals(mbixid, 'Faii_evap_18O',  fld_evap_18O,  lSize)
          call mbGetCellTagVals(mbixid, 'Faii_evap_HDO',  fld_evap_HDO,  lSize)
       end if
       do n=1,lSize
          if (lat_data(n) > 0.0_r8) then
             ic = c_inh_ir
          else
             ic = c_ish_ir
          endif
          ca_i =  area_data(n) * ifrac_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i
          nf = f_hmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*fld_melth(n)
          nf = f_hswnet; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_swnet(n) &
               - ca_i*fld_swpen(n)
          nf = f_hlwup ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_lwup(n)
          nf = f_hlatv ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_lat(n)
          nf = f_hsen  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_sen(n)
          nf = f_wmelt ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*fld_meltw(n)
          nf = f_wevap ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_evap(n)
          if ( flds_wiso_ice ) then
             nf = f_wmelt_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*fld_meltw_16O(n)
             nf = f_wmelt_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*fld_meltw_18O(n)
             nf = f_wmelt_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*fld_meltw_HDO(n)
             nf = f_wevap_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_evap_16O(n)
             nf = f_wevap_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_evap_18O(n)
             nf = f_wevap_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_evap_HDO(n)
          end if
       end do
       deallocate(fld_melth, fld_meltw, fld_swpen, fld_swnet)
       deallocate(fld_lwup, fld_lat, fld_sen, fld_evap)
       if (flds_wiso_ice) then
          deallocate(fld_meltw_16O, fld_meltw_18O, fld_meltw_HDO)
          deallocate(fld_evap_16O,  fld_evap_18O,  fld_evap_HDO)
       end if
    end if

    if (present(do_x2i)) then
       allocate(fld_lwdn(lSize), fld_rain(lSize), fld_snow(lSize))
       allocate(fld_frazil(lSize), fld_q(lSize), fld_rofi(lSize))
       call mbGetCellTagVals(mbixid, 'Faxa_lwdn',   fld_lwdn,   lSize)
       call mbGetCellTagVals(mbixid, 'Faxa_rain',   fld_rain,   lSize)
       call mbGetCellTagVals(mbixid, 'Faxa_snow',   fld_snow,   lSize)
       call mbGetCellTagVals(mbixid, 'Fioo_frazil', fld_frazil, lSize)
       call mbGetCellTagVals(mbixid, 'Fioo_q',      fld_q,      lSize)
       call mbGetCellTagVals(mbixid, 'Fixx_rofi',   fld_rofi,   lSize)
       if ( flds_wiso_ice_x2i ) then
          allocate(fld_rain_16O(lSize), fld_rain_18O(lSize), fld_rain_HDO(lSize))
          allocate(fld_snow_16O(lSize), fld_snow_18O(lSize), fld_snow_HDO(lSize))
          call mbGetCellTagVals(mbixid, 'Faxa_rain_16O', fld_rain_16O, lSize)
          call mbGetCellTagVals(mbixid, 'Faxa_rain_18O', fld_rain_18O, lSize)
          call mbGetCellTagVals(mbixid, 'Faxa_rain_HDO', fld_rain_HDO, lSize)
          call mbGetCellTagVals(mbixid, 'Faxa_snow_16O', fld_snow_16O, lSize)
          call mbGetCellTagVals(mbixid, 'Faxa_snow_18O', fld_snow_18O, lSize)
          call mbGetCellTagVals(mbixid, 'Faxa_snow_HDO', fld_snow_HDO, lSize)
       end if
       do n=1,lSize
          if (lat_data(n) > 0.0_r8) then
             ic = c_inh_is
          else
             ic = c_ish_is
          endif
          ca_o =  area_data(n) * ofrac_data(n)
          ca_i =  area_data(n) * ifrac_data(n)
          nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i
          nf = f_hlwdn; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_lwdn(n)
          nf = f_wrain; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_rain(n)
          nf = f_wsnow; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_snow(n)
          nf = f_wioff; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_rofi(n)
          nf = f_wfrz ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + &
               (ca_o+ca_i)*max(0.0_r8,fld_frazil(n))
          nf = f_hfrz ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - &
               (ca_o+ca_i)*max(0.0_r8,fld_q(n))
          if ( flds_wiso_ice_x2i ) then
             nf = f_wrain_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_rain_16O(n)
             nf = f_wrain_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_rain_18O(n)
             nf = f_wrain_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_rain_HDO(n)
             nf = f_wsnow_16O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_snow_16O(n)
             nf = f_wsnow_18O;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_snow_18O(n)
             nf = f_wsnow_HDO;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_snow_HDO(n)
          end if
       end do
       ic = c_inh_is
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

       ic = c_ish_is
       budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
       budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

       deallocate(fld_lwdn, fld_rain, fld_snow, fld_frazil, fld_q, fld_rofi)
       if (flds_wiso_ice_x2i) then
          deallocate(fld_rain_16O, fld_rain_18O, fld_rain_HDO)
          deallocate(fld_snow_16O, fld_snow_18O, fld_snow_HDO)
       end if
    end if

    deallocate(area_data, lat_data, ofrac_data, ifrac_data)
    first_time = .false.

  end subroutine seq_diag_ice_moab

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

  SUBROUTINE seq_diag_print_moab(EClock, stop_alarm, do_bgc_budg, &
       budg_print_inst,  budg_print_daily,  budg_print_month,  &
       budg_print_ann,  budg_print_ltann,  budg_print_ltend, infodata)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) , intent(in) :: EClock
    logical          , intent(in) :: stop_alarm
    logical          , intent(in) :: do_bgc_budg
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
    real(r8) :: dataGpr (f_size,c_size,p_size) ! values to print, scaled and such
    integer, parameter :: nisotopes = 3
    character(len=5), parameter :: isoname(nisotopes) = (/ 'H216O',   'H218O',   '  HDO'   /)
    integer, parameter          :: iso0(nisotopes)    = (/ f_16O,     f_18O,     f_hdO     /)
    integer, parameter          :: isof(nisotopes)    = (/ f_16O_end, f_18O_end, f_hdO_end /)

    !----- formats -----
    character(*),parameter :: subName = '(seq_diag_print_moab) '
    character(*),parameter :: F00   = "('(seq_diag_print_moab) ',4a)"

    !----- formats -----
    character(*),parameter :: FAH="(4a,i9,i6)"
    character(*),parameter :: FA0= "('    ',12x,6(6x,a8,1x))"
    character(*),parameter :: FA1= "('    ',a12,6f15.8)"
    character(*),parameter :: FA0r="('    ',12x,8(6x,a8,1x))"
    character(*),parameter :: FA1r="('    ',a12,8f15.8)"

    !-------------------------------------------------------------------------------

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

             if (do_bgc_budg) then
                call seq_diagBGC_preprint_moab()
             endif

             call seq_diag_sum0_moab()
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

          ! ---- doprint ---- doprint ---- doprint ----

          if (do_bgc_budg) then
             call seq_diagBGC_print_moab(EClock, ip, plev) 
          endif

       endif  ! plev > 0
    enddo  ! ip = 1,p_size

  end subroutine seq_diag_print_moab

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

end module seq_diag_moab
