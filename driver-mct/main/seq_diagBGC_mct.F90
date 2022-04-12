!===============================================================================
!
! !MODULE: seq_diagBGC_mod -- computes spatial \& time averages of fluxed BGC
!                              quatities
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
!       carbon flux   ~ (kg-C/s)/m^2
!
! !REVISION HISTORY:
!    2022-apr-12 - J. Wolfe    - initial budget implementation
!
! !INTERFACE: ------------------------------------------------------------------

module seq_diagBGC_mct
  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8, in=>shr_kind_in
  use shr_kind_mod, only: i8 => shr_kind_i8,  cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod, only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod, only : shr_mpi_max, shr_mpi_sum
  use shr_const_mod, only: shr_const_rearth, shr_const_pi, shr_const_isspval
  use shr_const_mod, only: shr_const_mwc, shr_const_mwco2
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

  public seq_diagBGC_zero_mct
  public seq_diagBGC_atm_mct
  public seq_diagBGC_lnd_mct
  public seq_diagBGC_rof_mct
  public seq_diagBGC_glc_mct
  public seq_diagBGC_ocn_mct
  public seq_diagBGC_ice_mct
  public seq_diagBGC_accum_mct
  public seq_diagBGC_sum0_mct
  public seq_diagBGC_preprint_mct
  public seq_diagBGC_print_mct
  public seq_diagBGC_avect_mct
  public seq_diagBGC_avloc_mct
  public seq_diagBGC_avdiff_mct

  !EOP

  !----------------------------------------------------------------------------
  ! Local data
  !----------------------------------------------------------------------------

  !----- local constants -----

  !--- C for component ---
  !--- "r" is receive in the coupler, "s" is send from the coupler

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
         ' c2g_glc',' g2c_glc',' c2a_inh',' a2c_inh',' c2a_ish',' a2c_ish', &
         ' c2a_lnd',' a2c_lnd',' c2a_ocn',' a2c_ocn' /)

  !--- F for field ---

  integer(in),parameter :: f_area      = 1     ! area (wrt to unit sphere)
  integer(in),parameter :: f_csurf     = 2     ! carbon : surface
  integer(in),parameter :: f_cblack    = 3     ! carbon : black carbon
  integer(in),parameter :: f_corgnc    = 4     ! carbon : organic

  integer(in),parameter :: f_size     = f_corgnc      ! Total array size of all elements
  integer(in),parameter :: f_a        = f_area        ! 1st index for area
  integer(in),parameter :: f_a_end    = f_area        ! last index for area
  integer(in),parameter :: f_c        = f_csurf       ! 1st index for carbon
  integer(in),parameter :: f_c_end    = f_corgnc      ! Last index for carbon

  character(len=12),parameter :: fname(f_size) = &
       (/'        area',' surface co2','black carbon','orgnc carbon' /)

  !--- P for period ---

  integer(in),parameter :: p_size = 5

  integer(in),parameter :: p_inst = 1
  integer(in),parameter :: p_day  = 2
  integer(in),parameter :: p_mon  = 3
  integer(in),parameter :: p_ann  = 4
  integer(in),parameter :: p_inf  = 5

  character(len=8),parameter :: pname(p_size) = &
       (/'    inst','   daily',' monthly','  annual','all_time' /)

  ! !PUBLIC DATA MEMBERS

  !--- time-averaged (annual?) global budget diagnostics ---
  !--- note: call sum0 then save budg_dataGBGC and budg_nsBGC on restart from/to root pe ---
  real(r8),public :: budg_dataL   (f_size,c_size,p_size) ! local sum, valid on all pes
  real(r8),public :: budg_dataGBGC(f_size,c_size,p_size) ! global sum, valid only on root pe
  real(r8),public :: budg_nsBGC   (f_size,c_size,p_size) ! counter, valid only on root pe
  real(r8),public :: dataGpr      (f_size,c_size,p_size) ! values to print, scaled and such

  character(len=*),parameter :: afldname  = 'aream'
  character(len=*),parameter :: latname   = 'lat'
  character(len=*),parameter :: afracname = 'afrac'
  character(len=*),parameter :: lfracname = 'lfrac'
  character(len=*),parameter :: lfrinname = 'lfrin'
  character(len=*),parameter :: ofracname = 'ofrac'
  character(len=*),parameter :: ifracname = 'ifrac'

  character(*),parameter :: modName = "(seq_diagBGC_mct) "

  integer(in),parameter :: debug = 0 ! internal debug level

  real(r8),parameter :: CO2toC = shr_const_mwc/shr_const_mwco2

  ! !PRIVATE DATA MEMBERS

  integer :: index_a2x_Faxa_bcphidry
  integer :: index_a2x_Faxa_bcphodry
  integer :: index_a2x_Faxa_bcphiwet
  integer :: index_a2x_Faxa_ocphidry
  integer :: index_a2x_Faxa_ocphodry
  integer :: index_a2x_Faxa_ocphiwet

  integer :: index_x2a_Fall_fco2_lnd
  integer :: index_x2a_Faoo_fco2_ocn

  integer :: index_l2x_Fall_fco2_lnd

  integer :: index_x2l_Faxa_bcphidry
  integer :: index_x2l_Faxa_bcphodry
  integer :: index_x2l_Faxa_bcphiwet
  integer :: index_x2l_Faxa_ocphidry
  integer :: index_x2l_Faxa_ocphodry
  integer :: index_x2l_Faxa_ocphiwet

  integer :: index_o2x_Faoo_fco2_ocn

  integer :: index_x2o_Faxa_bcphidry
  integer :: index_x2o_Faxa_bcphodry
  integer :: index_x2o_Faxa_bcphiwet
  integer :: index_x2o_Faxa_ocphidry
  integer :: index_x2o_Faxa_ocphodry
  integer :: index_x2o_Faxa_ocphiwet

  integer :: index_x2i_Faxa_bcphidry
  integer :: index_x2i_Faxa_bcphodry
  integer :: index_x2i_Faxa_bcphiwet
  integer :: index_x2i_Faxa_ocphidry
  integer :: index_x2i_Faxa_ocphodry
  integer :: index_x2i_Faxa_ocphiwet

  integer :: index_g2x_Fogg_rofl
  integer :: index_g2x_Fogg_rofi
  integer :: index_g2x_Figg_rofi

  !===============================================================================
contains
  !===============================================================================

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_zero_mct - zero out global budget diagnostic data.
  !
  ! !DESCRIPTION:
  !    Zero out global budget diagnostic data.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_zero_mct(EClock,mode)

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(in),optional :: EClock
    character(len=*), intent(in),optional :: mode

    !EOP

    integer(IN) :: ip,yr,mon,day,sec
    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_zero_mct) '

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
             budg_dataGBGC(:,:,ip) = 0.0_r8
             budg_nsBGC(:,:,ip) = 0.0_r8
          endif
          if (ip==p_day .and. sec==0) then
             budg_dataL(:,:,ip) = 0.0_r8
             budg_dataGBGC(:,:,ip) = 0.0_r8
             budg_nsBGC(:,:,ip) = 0.0_r8
          endif
          if (ip==p_mon .and. day==1 .and. sec==0) then
             budg_dataL(:,:,ip) = 0.0_r8
             budg_dataGBGC(:,:,ip) = 0.0_r8
             budg_nsBGC(:,:,ip) = 0.0_r8
          endif
          if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
             budg_dataL(:,:,ip) = 0.0_r8
             budg_dataGBGC(:,:,ip) = 0.0_r8
             budg_nsBGC(:,:,ip) = 0.0_r8
          endif
       enddo
    endif

    if (present(mode)) then
       if (trim(mode) == 'inst') then
          budg_dataL(:,:,p_inst) = 0.0_r8
          budg_dataGBGC(:,:,p_inst) = 0.0_r8
          budg_nsBGC(:,:,p_inst) = 0.0_r8
       elseif (trim(mode) == 'day') then
          budg_dataL(:,:,p_day) = 0.0_r8
          budg_dataGBGC(:,:,p_day) = 0.0_r8
          budg_nsBGC(:,:,p_day) = 0.0_r8
       elseif (trim(mode) == 'mon') then
          budg_dataL(:,:,p_mon) = 0.0_r8
          budg_dataGBGC(:,:,p_mon) = 0.0_r8
          budg_nsBGC(:,:,p_mon) = 0.0_r8
       elseif (trim(mode) == 'ann') then
          budg_dataL(:,:,p_ann) = 0.0_r8
          budg_dataGBGC(:,:,p_ann) = 0.0_r8
          budg_nsBGC(:,:,p_ann) = 0.0_r8
       elseif (trim(mode) == 'inf') then
          budg_dataL(:,:,p_inf) = 0.0_r8
          budg_dataGBGC(:,:,p_inf) = 0.0_r8
          budg_nsBGC(:,:,p_inf) = 0.0_r8
       elseif (trim(mode) == 'all') then
          budg_dataL(:,:,:) = 0.0_r8
          budg_dataGBGC(:,:,:) = 0.0_r8
          budg_nsBGC(:,:,:) = 0.0_r8
       else
          call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
       endif
    endif

  end subroutine seq_diagBGC_zero_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_accum_mct - accum out global budget diagnostic data.
  !
  ! !DESCRIPTION:
  !    Accum out global budget diagnostic data.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_accum_mct()

    ! !INPUT/OUTPUT PARAMETERS:

    !EOP

    integer(in) :: ip

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_accum_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    do ip = p_inst+1,p_size
       budg_dataL(:,:,ip) = budg_dataL(:,:,ip) + budg_dataL(:,:,p_inst)
    enddo
    budg_nsBGC(:,:,:) = budg_nsBGC(:,:,:) + 1.0_r8

  end subroutine seq_diagBGC_accum_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_sum0_mct - sum local to global on root
  !
  ! !DESCRIPTION:
  !    Sum local values to global on root
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_sum0_mct()

    ! !INPUT/OUTPUT PARAMETERS:

    !EOP

    real(r8) :: budg_dataGtmp(f_size,c_size,p_size) ! temporary sum
    integer(in)      :: mpicom      ! mpi comm
    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_sum0_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,mpicom=mpicom)
    budg_dataGtmp = 0.0_r8
    call shr_mpi_sum(budg_dataL,budg_dataGtmp,mpicom,subName)
    budg_dataGBGC = budg_dataGBGC + budg_dataGtmp
    budg_dataL = 0.0_r8

  end subroutine seq_diagBGC_sum0_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_atm_mct - compute global atm input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global atm input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_atm_mct( atm, frac_a, infodata, do_a2x, do_x2a)

    ! !INPUT/OUTPUT PARAMETERS:

    type(component_type)    , intent(in)           :: atm    ! component type for instance1
    type(mct_aVect)         , intent(in)           :: frac_a ! frac bundle
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_a2x
    logical                 , intent(in), optional :: do_x2a

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: a2x_a         ! model to drv bundle
    type(mct_aVect), pointer :: x2a_a         ! drv to model bundle
    type(mct_ggrid), pointer :: dom_a
    character(CL)            :: atm_gnam      ! atm grid
    character(CL)            :: lnd_gnam      ! lnd grid
    integer(in)              :: k,n,ic,nf,ip  ! generic index
    integer(in)              :: kArea         ! index of area field in aVect
    integer(in)              :: kLat          ! index of lat field in aVect
    integer(in)              :: kl,ka,ko,ki   ! fraction indices
    integer(in)              :: lSize         ! size of aVect
    real(r8)                 :: ca_a          ! area of a grid cell
    logical,save             :: first_time    = .true.
    logical,save             :: samegrid_al   ! samegrid atm and lnd

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_atm_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    dom_a => component_get_dom_cx(atm)
    a2x_a => component_get_c2x_cx(atm)
    x2a_a => component_get_x2c_cx(atm)

    kArea = mct_aVect_indexRA(dom_a%data,afldname)
    kLat  = mct_aVect_indexRA(dom_a%data,latname)
    ka    = mct_aVect_indexRA(frac_a,afracname)
    ko    = mct_aVect_indexRA(frac_a,ofracname)
    ki    = mct_aVect_indexRA(frac_a,ifracname)
    if (first_time) then
       call seq_infodata_getData(infodata , &
            lnd_gnam=lnd_gnam             , &
            atm_gnam=atm_gnam             )
       samegrid_al = .true.
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
    end if

    if (samegrid_al) then
       kl = mct_aVect_indexRA(frac_a,lfracname)
    else
       kl = mct_aVect_indexRA(frac_a,lfrinname)
    endif

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    ip = p_inst

    if (present(do_a2x)) then
       if (first_time) then
          index_a2x_Faxa_bcphidry = mct_aVect_indexRA(a2x_a,'Faxa_bcphidry')
          index_a2x_Faxa_bcphodry = mct_aVect_indexRA(a2x_a,'Faxa_bcphodry')
          index_a2x_Faxa_bcphiwet = mct_aVect_indexRA(a2x_a,'Faxa_bcphiwet')
          index_a2x_Faxa_ocphidry = mct_aVect_indexRA(a2x_a,'Faxa_ocphidry')
          index_a2x_Faxa_ocphodry = mct_aVect_indexRA(a2x_a,'Faxa_ocphodry')
          index_a2x_Faxa_ocphiwet = mct_aVect_indexRA(a2x_a,'Faxa_ocphiwet')
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
             nf = f_cblack; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*a2x_a%rAttr(index_a2x_Faxa_bcphidry,n) &
                                                                        + ca_a*a2x_a%rAttr(index_a2x_Faxa_bcphodry,n) &
                                                                        + ca_a*a2x_a%rAttr(index_a2x_Faxa_bcphiwet,n)
             nf = f_corgnc; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*a2x_a%rAttr(index_a2x_Faxa_ocphidry,n) &
                                                                        + ca_a*a2x_a%rAttr(index_a2x_Faxa_ocphodry,n) &
                                                                        + ca_a*a2x_a%rAttr(index_a2x_Faxa_ocphiwet,n)
          enddo
       enddo
    end if

    if (present(do_x2a)) then
       if (first_time) then
          index_x2a_Fall_fco2_lnd = mct_aVect_indexRA(x2a_a,'Fall_fco2_lnd',perrWith='quiet')
          index_x2a_Faoo_fco2_ocn = mct_aVect_indexRA(x2a_a,'Faoo_fco2_ocn',perrWith='quiet')
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

             nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a
             if (index_x2a_Fall_fco2_lnd /= 0) then
                nf = f_csurf ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                    + ca_a*x2a_a%rAttr(index_x2a_Fall_fco2_lnd,n)*CO2toC
             end if
             if (index_x2a_Faoo_fco2_ocn /= 0) then
                nf = f_csurf ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                    + ca_a*x2a_a%rAttr(index_x2a_Faoo_fco2_ocn,n)*CO2toC
             end if

          enddo
       enddo
    end if

    first_time = .false.

  end subroutine seq_diagBGC_atm_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_lnd_mct - compute global lnd input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global lnd input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_lnd_mct( lnd, frac_l, infodata, do_l2x, do_x2l)

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

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_lnd_mct) '

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
    kl    = mct_aVect_indexRA(frac_l,lfrinname)

    if (present(do_l2x)) then
       if (first_time) then
          index_l2x_Fall_fco2_lnd = mct_aVect_indexRA(l2x_l,'Fall_fco2_lnd',perrWith='quiet')
       end if

       lSize = mct_avect_lSize(l2x_l)
       ic = c_lnd_lr
       do n=1,lSize
          ca_l =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          if (index_x2a_Fall_fco2_lnd /= 0) then
             nf = f_csurf ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                 + ca_l*l2x_l%rAttr(index_l2x_Fall_fco2_lnd,n)*CO2toC
          end if
       end do
    end if

    if (present(do_x2l)) then
       if (first_time) then
          index_x2l_Faxa_bcphidry = mct_aVect_indexRA(x2l_l,'Faxa_bcphidry')
          index_x2l_Faxa_bcphodry = mct_aVect_indexRA(x2l_l,'Faxa_bcphodry')
          index_x2l_Faxa_bcphiwet = mct_aVect_indexRA(x2l_l,'Faxa_bcphiwet')
          index_x2l_Faxa_ocphidry = mct_aVect_indexRA(x2l_l,'Faxa_ocphidry')
          index_x2l_Faxa_ocphodry = mct_aVect_indexRA(x2l_l,'Faxa_ocphodry')
          index_x2l_Faxa_ocphiwet = mct_aVect_indexRA(x2l_l,'Faxa_ocphiwet')
       end if

       lSize = mct_avect_lSize(x2l_l)
       ic = c_lnd_ls
       do n=1,lSize
          ca_l =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          nf = f_cblack; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*x2l_l%rAttr(index_x2l_Faxa_bcphidry,n) &
                                                                     + ca_l*x2l_l%rAttr(index_x2l_Faxa_bcphodry,n) &
                                                                     + ca_l*x2l_l%rAttr(index_x2l_Faxa_bcphiwet,n)
          nf = f_corgnc; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*x2l_l%rAttr(index_x2l_Faxa_ocphidry,n) &
                                                                     + ca_l*x2l_l%rAttr(index_x2l_Faxa_ocphodry,n) &
                                                                     + ca_l*x2l_l%rAttr(index_x2l_Faxa_ocphiwet,n)
       end do
    end if

    first_time = .false.

  end subroutine seq_diagBGC_lnd_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_rof_mct - compute global rof input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global rof input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_rof_mct( rof, frac_r, infodata)

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

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_rof_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    dom_r => component_get_dom_cx(rof)
    r2x_r => component_get_c2x_cx(rof)
    x2r_r => component_get_x2c_cx(rof)

!    if (first_time) then
!    end if

    ip = p_inst
    ic = c_rof_rr
    kArea = mct_aVect_indexRA(dom_r%data,afldname)
    lSize = mct_avect_lSize(x2r_r)
    do n=1,lSize
       ca_r =  dom_r%data%rAttr(kArea,n)
    end do

!    if (first_time) then
!    end if

    ip = p_inst
    ic = c_rof_rs
    kArea = mct_aVect_indexRA(dom_r%data,afldname)
    lSize = mct_avect_lSize(r2x_r)
    do n=1,lSize
       ca_r =  dom_r%data%rAttr(kArea,n)
    end do

    first_time = .false.

  end subroutine seq_diagBGC_rof_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_glc_mct - compute global glc input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global glc input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_glc_mct( glc, frac_g, infodata)

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
    character(*),parameter :: subName = '(seq_diagBGC_glc_mct) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! add values found in this bundle to the budget table
    !---------------------------------------------------------------------------

    dom_g => component_get_dom_cx(glc)
    g2x_g => component_get_c2x_cx(glc)
    x2g_g => component_get_x2c_cx(glc)

!    if (first_time) then
!    end if

    ip = p_inst
    ic = c_glc_gs
    kArea = mct_aVect_indexRA(dom_g%data,afldname)
    lSize = mct_avect_lSize(g2x_g)
    do n=1,lSize
       ca_g =  dom_g%data%rAttr(kArea,n)
    end do

    first_time = .false.

  end subroutine seq_diagBGC_glc_mct

  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_ocn_mct - compute global ocn input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global ocn input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_ocn_mct( ocn, xao_o, frac_o, infodata, do_o2x, do_x2o, do_xao)

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

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_ocn_mct) '

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

    if (present(do_o2x)) then
       if (first_time) then
          index_o2x_Faoo_fco2_ocn = mct_aVect_indexRA(o2x_o,'Faoo_fco2_ocn',perrWith='quiet')
       end if

       lSize = mct_avect_lSize(o2x_o)
       ic = c_ocn_or
       do n=1,lSize
          ca_o =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
          ca_i =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
          nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o
          if (index_x2a_Faoo_fco2_ocn /= 0) then
             nf = f_csurf; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                + (ca_o+ca_i)*o2x_o%rAttr(index_o2x_Faoo_fco2_ocn,n)*CO2toC
          end if
       end do
    end if

    if (present(do_x2o)) then
       if (first_time) then
          index_x2o_Faxa_bcphidry = mct_aVect_indexRA(x2o_o,'Faxa_bcphidry')
          index_x2o_Faxa_bcphodry = mct_aVect_indexRA(x2o_o,'Faxa_bcphodry')
          index_x2o_Faxa_bcphiwet = mct_aVect_indexRA(x2o_o,'Faxa_bcphiwet')
          index_x2o_Faxa_ocphidry = mct_aVect_indexRA(x2o_o,'Faxa_ocphidry')
          index_x2o_Faxa_ocphodry = mct_aVect_indexRA(x2o_o,'Faxa_ocphodry')
          index_x2o_Faxa_ocphiwet = mct_aVect_indexRA(x2o_o,'Faxa_ocphiwet')
       end if

       lSize = mct_avect_lSize(x2o_o)
       ic = c_ocn_os
       do n=1,lSize
          ca_o =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
          ca_i =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o

          nf = f_cblack; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_bcphidry,n) &
                                                                     + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_bcphodry,n) &
                                                                     + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_bcphiwet,n)
          nf = f_corgnc; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_ocphidry,n) &
                                                                     + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_ocphodry,n) &
                                                                     + (ca_o+ca_i)*x2o_o%rAttr(index_x2o_Faxa_ocphiwet,n)
       end do
    end if

    first_time = .false.

  end subroutine seq_diagBGC_ocn_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_ice_mct - compute global ice input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global ice input/output flux diagnostics
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_ice_mct( ice, frac_i, infodata, do_i2x, do_x2i)

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

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_ice_mct) '

    !-------------------------------------------------------------------------------

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
!       if (first_time) then
!       endif

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
       end do
    end if

    if (present(do_x2i)) then
       if (first_time) then
          index_x2i_Faxa_bcphidry = mct_aVect_indexRA(x2i_i,'Faxa_bcphidry')
          index_x2i_Faxa_bcphodry = mct_aVect_indexRA(x2i_i,'Faxa_bcphodry')
          index_x2i_Faxa_bcphiwet = mct_aVect_indexRA(x2i_i,'Faxa_bcphiwet')
          index_x2i_Faxa_ocphidry = mct_aVect_indexRA(x2i_i,'Faxa_ocphidry')
          index_x2i_Faxa_ocphodry = mct_aVect_indexRA(x2i_i,'Faxa_ocphodry')
          index_x2i_Faxa_ocphiwet = mct_aVect_indexRA(x2i_i,'Faxa_ocphiwet')
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
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i
          nf = f_cblack; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*x2i_i%rAttr(index_x2i_Faxa_bcphidry,n) &
                                                                     + ca_i*x2i_i%rAttr(index_x2i_Faxa_bcphodry,n) &
                                                                     + ca_i*x2i_i%rAttr(index_x2i_Faxa_bcphiwet,n)
          nf = f_corgnc; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*x2i_i%rAttr(index_x2i_Faxa_ocphidry,n) &
                                                                     + ca_i*x2i_i%rAttr(index_x2i_Faxa_ocphodry,n) &
                                                                     + ca_i*x2i_i%rAttr(index_x2i_Faxa_ocphiwet,n)
       end do
    end if

    first_time = .false.

  end subroutine seq_diagBGC_ice_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_preprint_mct - print global BGC budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global BGC budget diagnostics.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diagBGC_preprint_mct

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    !EOP

    !--- local ---

    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! prepare instantaneous budget data for printing
    !-------------------------------------------------------------------------------

    call seq_diagBGC_sum0_mct()
    dataGpr = budg_dataGBGC

    !  old budget normalizations (global area and 1e10 for carbon)
    dataGpr = dataGpr/(4.0_r8*shr_const_pi)
    dataGpr(f_c:f_c_end,:,:) = dataGpr(f_c:f_c_end,:,:) * 1.0e10_r8
    dataGpr = dataGpr/budg_nsBGC

  end subroutine seq_diagBGC_preprint_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_print_mct - print global BGC budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global BGC budget diagnostics.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diagBGC_print_mct(EClock, ip, plev)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) , intent(in) :: EClock
    integer(in)      , intent(in) :: ip
    integer(in)      , intent(in) :: plev        ! print level

    !EOP

    !--- local ---
    integer(in)      :: ic,nf,is    ! data array indicies
    integer(in)      :: ica,icl,icn,ics,ico
    integer(in)      :: icar,icxs,icxr,icas
    integer(in)      :: cdate,sec   ! coded date, seconds
    integer(in)      :: yr,mon,day  ! date
    integer(in)      :: iam         ! pe number
    character(len=40):: str         ! string

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_print_mct) '
    character(*),parameter :: F00   = "('(seq_diagBGC_print_mct) ',4a)"

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

    call seq_comm_setptrs(CPLID,iam=iam)
    call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
         curr_mon=mon,curr_day=day,curr_tod=sec)
    cdate = yr*10000+mon*100+day

   if (plev > 0) then
      ! ---- doprint ---- doprint ---- doprint ----

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
            write(logunit,FAH) subname,trim(str)//' CARBON BUDGET (kg-C/m2s*1e10): period = ',trim(pname(ip)),': date = ',cdate,sec
            write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
            do nf = f_c, f_c_end
               write(logunit,FA1)    fname(nf),dataGpr(nf,ica,ip),dataGpr(nf,icl,ip), &
                    dataGpr(nf,icn,ip),dataGpr(nf,ics,ip),dataGpr(nf,ico,ip), &
                    dataGpr(nf,ica,ip)+dataGpr(nf,icl,ip)+ &
                    dataGpr(nf,icn,ip)+dataGpr(nf,ics,ip)+dataGpr(nf,ico,ip)
            enddo
            write(logunit,FA1)    '   *SUM*'   ,sum(dataGpr(f_c:f_c_end,ica,ip)),sum(dataGpr(f_c:f_c_end,icl,ip)), &
                 sum(dataGpr(f_c:f_c_end,icn,ip)),sum(dataGpr(f_c:f_c_end,ics,ip)),sum(dataGpr(f_c:f_c_end,ico,ip)), &
                 sum(dataGpr(f_c:f_c_end,ica,ip))+sum(dataGpr(f_c:f_c_end,icl,ip))+ &
                 sum(dataGpr(f_c:f_c_end,icn,ip))+sum(dataGpr(f_c:f_c_end,ics,ip))+sum(dataGpr(f_c:f_c_end,ico,ip))
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
            write(logunit,FAH) subname,trim(str)//' CARBON BUDGET (kg-C/m2s*1e10): period = ',trim(pname(ip)),': date = ',cdate,sec
            write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
            do nf = f_c, f_c_end
               write(logunit,FA1)    fname(nf),-dataGpr(nf,icar,ip),dataGpr(nf,icxs,ip), &
                    dataGpr(nf,icxr,ip),-dataGpr(nf,icas,ip), &
                    -dataGpr(nf,icar,ip)+dataGpr(nf,icxs,ip)+ &
                    dataGpr(nf,icxr,ip)-dataGpr(nf,icas,ip)
            enddo
            write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_c:f_c_end,icar,ip)),sum(dataGpr(f_c:f_c_end,icxs,ip)), &
                 sum(dataGpr(f_c:f_c_end,icxr,ip)),-sum(dataGpr(f_c:f_c_end,icas,ip)), &
                 -sum(dataGpr(f_c:f_c_end,icar,ip))+sum(dataGpr(f_c:f_c_end,icxs,ip))+ &
                 sum(dataGpr(f_c:f_c_end,icxr,ip))-sum(dataGpr(f_c:f_c_end,icas,ip))

         enddo
      endif   ! plev

      ! ---------------------------------------------------------
      ! ---- net summary budgets ----
      ! ---------------------------------------------------------

      if (plev >= 1) then

         write(logunit,*) ' '
         write(logunit,FAH) subname,'NET CARBON BUDGET (kg-C/m2s*1e10): period = ',trim(pname(ip)),': date = ',cdate,sec
         write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
         do nf = f_c, f_c_end
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
         write(logunit,FA1r)'   *SUM*',sum(dataGpr(f_c:f_c_end,c_atm_ar,ip))+sum(dataGpr(f_c:f_c_end,c_atm_as,ip)), &
              sum(dataGpr(f_c:f_c_end,c_lnd_lr,ip))+sum(dataGpr(f_c:f_c_end,c_lnd_ls,ip)), &
              sum(dataGpr(f_c:f_c_end,c_rof_rr,ip))+sum(dataGpr(f_c:f_c_end,c_rof_rs,ip)), &
              sum(dataGpr(f_c:f_c_end,c_ocn_or,ip))+sum(dataGpr(f_c:f_c_end,c_ocn_os,ip)), &
              sum(dataGpr(f_c:f_c_end,c_inh_ir,ip))+sum(dataGpr(f_c:f_c_end,c_inh_is,ip)), &
              sum(dataGpr(f_c:f_c_end,c_ish_ir,ip))+sum(dataGpr(f_c:f_c_end,c_ish_is,ip)), &
              sum(dataGpr(f_c:f_c_end,c_glc_gr,ip))+sum(dataGpr(f_c:f_c_end,c_glc_gs,ip)), &
              sum(dataGpr(f_c:f_c_end,c_atm_ar,ip))+sum(dataGpr(f_c:f_c_end,c_atm_as,ip))+ &
              sum(dataGpr(f_c:f_c_end,c_lnd_lr,ip))+sum(dataGpr(f_c:f_c_end,c_lnd_ls,ip))+ &
              sum(dataGpr(f_c:f_c_end,c_rof_rr,ip))+sum(dataGpr(f_c:f_c_end,c_rof_rs,ip))+ &
              sum(dataGpr(f_c:f_c_end,c_ocn_or,ip))+sum(dataGpr(f_c:f_c_end,c_ocn_os,ip))+ &
              sum(dataGpr(f_c:f_c_end,c_inh_ir,ip))+sum(dataGpr(f_c:f_c_end,c_inh_is,ip))+ &
              sum(dataGpr(f_c:f_c_end,c_ish_ir,ip))+sum(dataGpr(f_c:f_c_end,c_ish_is,ip))+ &
              sum(dataGpr(f_c:f_c_end,c_glc_gr,ip))+sum(dataGpr(f_c:f_c_end,c_glc_gs,ip))
      endif

      write(logunit,*) ' '
      ! ---- doprint ---- doprint ---- doprint ----
   endif  ! plev > 0

  end subroutine seq_diagBGC_print_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_avect_mct - print global budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global diagnostics for AV/ID.
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diagBGC_avect_mct(infodata, id, av, dom, gsmap, comment)

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
    character(*),parameter :: subName = '(seq_diagBGC_avect_mct) '
    character(*),parameter :: F00   = "('(seq_diagBGC_avect_mct) ',4a)"

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

  end subroutine seq_diagBGC_avect_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_avloc_mct - print local budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print local diagnostics for AV/ID.
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diagBGC_avloc_mct(av, comment)

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
    character(*),parameter :: subName = '(seq_diagBGC_avloc_mct) '
    character(*),parameter :: F00   = "('(seq_diagBGC_avloc_mct) ',4a)"

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

  end subroutine seq_diagBGC_avloc_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_avdiff_mct - print global budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global diagnostics for AV/ID.
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diagBGC_avdiff_mct(AV1,AV2,ID,comment)

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
    character(*),parameter :: subName = '(seq_diagBGC_avdiff_mct) '
    character(*),parameter :: F00   = "('(seq_diagBGC_avdiff_mct) ',4a)"

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

  end subroutine seq_diagBGC_avdiff_mct

end module seq_diagBGC_mct
