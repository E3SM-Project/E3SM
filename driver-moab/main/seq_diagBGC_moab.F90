!===============================================================================
!
! !MODULE: seq_diagBGC_moab -- computes spatial & time averages of fluxed BGC
!                              quantities (MOAB-native version)
!
! !DESCRIPTION:
!    The coupler is required to do certain diagnostics, those calculations are
!    located in this module.  This is the MOAB-native version of seq_diagBGC_mct.
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
!    2022-apr-12 - J. Wolfe    - initial budget implementation (MCT version)
!    2026-apr-10 - R. Jacob    - converted to MOAB-native version
!
! !INTERFACE: ------------------------------------------------------------------

module seq_diagBGC_moab
  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8, in=>shr_kind_in
  use shr_kind_mod, only: cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod, only : shr_sys_abort
  use shr_mpi_mod, only : shr_mpi_sum
  use shr_const_mod, only: shr_const_pi
  use shr_const_mod, only: shr_const_mwc, shr_const_mwco2
  use mct_mod, only: mct_avect, mct_aVect_indexRA
  use seq_comm_mct, only: logunit, cplid, seq_comm_setptrs, &
       mbaxid, mblxid, mboxid, mbixid, mbrxid
  use shr_moab_mod, only: mbGetnCells, mbGetCellTagVals
  use esmf, only : esmf_clock
  use seq_timemgr_mod, only : seq_timemgr_EClockGetData
  use component_type_mod, only : COMPONENT_GET_C2X_CX, &
       COMPONENT_GET_X2C_CX, COMPONENT_TYPE
  use seq_infodata_mod, only : seq_infodata_type, seq_infodata_getdata

  implicit none
  save
  private

  ! !PUBLIC TYPES:

  ! none

  !PUBLIC MEMBER FUNCTIONS:

  public seq_diagBGC_zero_moab
  public seq_diagBGC_atm_moab
  public seq_diagBGC_lnd_moab
  public seq_diagBGC_rof_moab
  public seq_diagBGC_glc_moab
  public seq_diagBGC_ocn_moab
  public seq_diagBGC_ice_moab
  public seq_diagBGC_accum_moab
  public seq_diagBGC_sum0_moab
  public seq_diagBGC_preprint_moab
  public seq_diagBGC_print_moab

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
  integer(in),parameter :: f_ioinor    = 5     ! carbon : ice-ocn inorganic
  integer(in),parameter :: f_ioorgn    = 6     ! carbon : ice-ocn organic

  integer(in),parameter :: f_size     = f_ioorgn      ! Total array size of all elements
  integer(in),parameter :: f_a        = f_area        ! 1st index for area
  integer(in),parameter :: f_a_end    = f_area        ! last index for area
  integer(in),parameter :: f_c        = f_csurf       ! 1st index for carbon
  integer(in),parameter :: f_c_end    = f_ioorgn      ! Last index for carbon

  character(len=12),parameter :: fname(f_size) = &
       (/'        area',' surface co2','black carbon','orgnc carbon' ,&
         ' i-o inorgnc',' i-o   orgnc'/)

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

  character(*),parameter :: modName = "(seq_diagBGC_moab) "

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
  integer :: index_x2o_Fioi_algae1
  integer :: index_x2o_Fioi_algae2
  integer :: index_x2o_Fioi_algae3
  integer :: index_x2o_Fioi_dic1
  integer :: index_x2o_Fioi_docr
  integer :: index_x2o_Fioi_doc1
  integer :: index_x2o_Fioi_doc2
  integer :: index_x2o_Fioi_doc3

  integer :: index_i2x_Fioi_algae1
  integer :: index_i2x_Fioi_algae2
  integer :: index_i2x_Fioi_algae3
  integer :: index_i2x_Fioi_dic1
  integer :: index_i2x_Fioi_docr
  integer :: index_i2x_Fioi_doc1
  integer :: index_i2x_Fioi_doc2
  integer :: index_i2x_Fioi_doc3

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
  ! !IROUTINE: seq_diagBGC_zero_moab - zero out global budget diagnostic data.
  !
  ! !DESCRIPTION:
  !    Zero out global budget diagnostic data.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_zero_moab(EClock,mode)

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(in),optional :: EClock
    character(len=*), intent(in),optional :: mode

    !EOP

    integer(IN) :: ip,yr,mon,day,sec
    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_zero_moab) '

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

  end subroutine seq_diagBGC_zero_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_accum_moab - accumulate global budget diagnostic data.
  !
  ! !DESCRIPTION:
  !    Accumulate global budget diagnostic data.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_accum_moab()

    !EOP

    integer(in) :: ip

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_accum_moab) '

    !-------------------------------------------------------------------------------

    do ip = p_inst+1,p_size
       budg_dataL(:,:,ip) = budg_dataL(:,:,ip) + budg_dataL(:,:,p_inst)
    enddo
    budg_nsBGC(:,:,:) = budg_nsBGC(:,:,:) + 1.0_r8

  end subroutine seq_diagBGC_accum_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_sum0_moab - sum local to global on root
  !
  ! !DESCRIPTION:
  !    Sum local values to global on root
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_sum0_moab()

    !EOP

    real(r8) :: budg_dataGtmp(f_size,c_size,p_size)
    integer(in)      :: mpicom
    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_sum0_moab) '

    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,mpicom=mpicom)
    budg_dataGtmp = 0.0_r8
    call shr_mpi_sum(budg_dataL,budg_dataGtmp,mpicom,subName)
    budg_dataGBGC = budg_dataGBGC + budg_dataGtmp
    budg_dataL = 0.0_r8

  end subroutine seq_diagBGC_sum0_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_atm_moab - compute global atm input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global atm input/output flux diagnostics using MOAB data structures.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_atm_moab( atm, infodata, do_a2x, do_x2a)

    ! !INPUT/OUTPUT PARAMETERS:

    type(component_type)    , intent(in)           :: atm      ! component type
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
    real(r8), allocatable    :: fld_bcphidry(:)
    real(r8), allocatable    :: fld_bcphodry(:)
    real(r8), allocatable    :: fld_bcphiwet(:)
    real(r8), allocatable    :: fld_ocphidry(:)
    real(r8), allocatable    :: fld_ocphodry(:)
    real(r8), allocatable    :: fld_ocphiwet(:)
    real(r8), allocatable    :: fld_fco2_lnd(:)
    real(r8), allocatable    :: fld_fco2_ocn(:)
    logical,save             :: first_time=.true. ! initialization flag
    logical,save             :: samegrid_al       ! samegrid atm and lnd

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_atm_moab) '

    !-------------------------------------------------------------------------------

    lSize = mbGetnCells(mbaxid)
    allocate(area_data(lSize), lat_data(lSize))
    allocate(afrac_data(lSize), lfrac_data(lSize), ofrac_data(lSize), ifrac_data(lSize))

    call mbGetCellTagVals(mbaxid, afldname,  area_data, lSize)
    call mbGetCellTagVals(mbaxid, latname,   lat_data,  lSize)
    call mbGetCellTagVals(mbaxid, afracname, afrac_data, lSize)
    call mbGetCellTagVals(mbaxid, ofracname, ofrac_data, lSize)
    call mbGetCellTagVals(mbaxid, ifracname, ifrac_data, lSize)

    if (first_time) then
       call seq_infodata_getData(infodata , &
            lnd_gnam=lnd_gnam             , &
            atm_gnam=atm_gnam             )
       samegrid_al = .true.
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
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
       allocate(fld_bcphidry(lSize), fld_bcphodry(lSize), fld_bcphiwet(lSize))
       allocate(fld_ocphidry(lSize), fld_ocphodry(lSize), fld_ocphiwet(lSize))

       call mbGetCellTagVals(mbaxid, 'Faxa_bcphidry', fld_bcphidry, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_bcphodry', fld_bcphodry, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_bcphiwet', fld_bcphiwet, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_ocphidry', fld_ocphidry, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_ocphodry', fld_ocphodry, lSize)
       call mbGetCellTagVals(mbaxid, 'Faxa_ocphiwet', fld_ocphiwet, lSize)

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
             nf = f_cblack; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_bcphidry(n) &
                                                                        + ca_a*fld_bcphodry(n) &
                                                                        + ca_a*fld_bcphiwet(n)
             nf = f_corgnc; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a*fld_ocphidry(n) &
                                                                        + ca_a*fld_ocphodry(n) &
                                                                        + ca_a*fld_ocphiwet(n)
          enddo
       enddo

       deallocate(fld_bcphidry, fld_bcphodry, fld_bcphiwet)
       deallocate(fld_ocphidry, fld_ocphodry, fld_ocphiwet)
    end if

    if (present(do_x2a)) then
       if (first_time) then
          av_tmp => component_get_x2c_cx(atm)
          index_x2a_Fall_fco2_lnd = mct_aVect_indexRA(av_tmp,'Fall_fco2_lnd',perrWith='quiet')
          index_x2a_Faoo_fco2_ocn = mct_aVect_indexRA(av_tmp,'Faoo_fco2_ocn',perrWith='quiet')
       end if

       if (index_x2a_Fall_fco2_lnd /= 0) then
          allocate(fld_fco2_lnd(lSize))
          call mbGetCellTagVals(mbaxid, 'Fall_fco2_lnd', fld_fco2_lnd, lSize)
       endif
       if (index_x2a_Faoo_fco2_ocn /= 0) then
          allocate(fld_fco2_ocn(lSize))
          call mbGetCellTagVals(mbaxid, 'Faoo_fco2_ocn', fld_fco2_ocn, lSize)
       endif
       ! area and fco2 accumulation loop
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

             nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_a

             if (index_x2a_Fall_fco2_lnd /= 0) then
                nf = f_csurf ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                  + ca_a*fld_fco2_lnd(n)*CO2toC
             endif
             if (index_x2a_Faoo_fco2_ocn /= 0) then
                nf = f_csurf ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                     + ca_a*fld_fco2_ocn(n)*CO2toC
             endif

          enddo
       enddo

       if (index_x2a_Fall_fco2_lnd /= 0) then
          deallocate(fld_fco2_lnd)
       end if

       if (index_x2a_Faoo_fco2_ocn /= 0) then
          deallocate(fld_fco2_ocn)
       end if

    end if

    deallocate(area_data, lat_data)
    deallocate(afrac_data, lfrac_data, ofrac_data, ifrac_data)

    first_time = .false.

  end subroutine seq_diagBGC_atm_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_lnd_moab - compute global lnd input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global lnd input/output flux diagnostics using MOAB data structures.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_lnd_moab( lnd, infodata, do_l2x, do_x2l)

    type(component_type)    , intent(in)           :: lnd      ! component type
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_l2x
    logical                 , intent(in), optional :: do_x2l

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp            ! for optional field detection
    integer(in)              :: n,ic,nf,ip        ! generic index
    integer(in)              :: lSize             ! size of mesh
    real(r8)                 :: ca_l              ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    real(r8), allocatable    :: lfrac_data(:)
    real(r8), allocatable    :: fld_fco2_lnd(:)
    real(r8), allocatable    :: fld_bcphidry(:)
    real(r8), allocatable    :: fld_bcphodry(:)
    real(r8), allocatable    :: fld_bcphiwet(:)
    real(r8), allocatable    :: fld_ocphidry(:)
    real(r8), allocatable    :: fld_ocphodry(:)
    real(r8), allocatable    :: fld_ocphiwet(:)
    logical,save             :: first_time=.true. ! initialization flag

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_lnd_moab) '

    !-------------------------------------------------------------------------------

    lSize = mbGetnCells(mblxid)
    allocate(area_data(lSize), lfrac_data(lSize))

    call mbGetCellTagVals(mblxid, afldname,  area_data,  lSize)
    call mbGetCellTagVals(mblxid, lfrinname, lfrac_data, lSize)

    ip = p_inst

    if (present(do_l2x)) then
       if (first_time) then
          av_tmp => component_get_c2x_cx(lnd)
          index_l2x_Fall_fco2_lnd = mct_aVect_indexRA(av_tmp,'Fall_fco2_lnd',perrWith='quiet')
       end if

       if (index_l2x_Fall_fco2_lnd /= 0) then
          allocate(fld_fco2_lnd(lSize))
          call mbGetCellTagVals(mblxid, 'Fall_fco2_lnd', fld_fco2_lnd, lSize)
       endif

       ic = c_lnd_lr
       do n=1,lSize
          ca_l = area_data(n) * lfrac_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          if (index_l2x_Fall_fco2_lnd /= 0) then
             nf = f_csurf ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                  + ca_l*fld_fco2_lnd(n)*CO2toC
          endif
       end do

       if (index_l2x_Fall_fco2_lnd /= 0) then
          deallocate(fld_fco2_lnd)
       end if
    end if

    if (present(do_x2l)) then
       allocate(fld_bcphidry(lSize), fld_bcphodry(lSize), fld_bcphiwet(lSize))
       allocate(fld_ocphidry(lSize), fld_ocphodry(lSize), fld_ocphiwet(lSize))

       call mbGetCellTagVals(mblxid, 'Faxa_bcphidry', fld_bcphidry, lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_bcphodry', fld_bcphodry, lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_bcphiwet', fld_bcphiwet, lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_ocphidry', fld_ocphidry, lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_ocphodry', fld_ocphodry, lSize)
       call mbGetCellTagVals(mblxid, 'Faxa_ocphiwet', fld_ocphiwet, lSize)

       ic = c_lnd_ls
       do n=1,lSize
          ca_l = area_data(n) * lfrac_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l
          nf = f_cblack; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_bcphidry(n) &
                                                                     + ca_l*fld_bcphodry(n) &
                                                                     + ca_l*fld_bcphiwet(n)
          nf = f_corgnc; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_l*fld_ocphidry(n) &
                                                                     + ca_l*fld_ocphodry(n) &
                                                                     + ca_l*fld_ocphiwet(n)
       end do

       deallocate(fld_bcphidry, fld_bcphodry, fld_bcphiwet)
       deallocate(fld_ocphidry, fld_ocphodry, fld_ocphiwet)
    end if

    deallocate(area_data, lfrac_data)

    first_time = .false.

  end subroutine seq_diagBGC_lnd_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_rof_moab - compute global rof input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global rof input/output flux diagnostics using MOAB data structures.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_rof_moab( rof, infodata)

    type(component_type)    , intent(in) :: rof      ! component type
    type(seq_infodata_type) , intent(in) :: infodata

    !EOP

    !----- local -----
    integer(in)              :: n,ic,nf,ip        ! generic index
    integer(in)              :: lSize             ! size of mesh
    real(r8)                 :: ca_r              ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    logical,save             :: first_time=.true. ! initialization flag

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_rof_moab) '

    !-------------------------------------------------------------------------------

    lSize = mbGetnCells(mbrxid)
    allocate(area_data(lSize))
    call mbGetCellTagVals(mbrxid, afldname, area_data, lSize)

    ip = p_inst
    ic = c_rof_rr
    do n=1,lSize
       ca_r = area_data(n)
    end do

    ip = p_inst
    ic = c_rof_rs
    do n=1,lSize
       ca_r = area_data(n)
    end do

    deallocate(area_data)

    first_time = .false.

  end subroutine seq_diagBGC_rof_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_glc_moab - compute global glc input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global glc input/output flux diagnostics.
  !     Delegates to the MCT version since there is no MOAB ID for glc.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_glc_moab( glc, infodata)

    type(component_type)    , intent(in) :: glc      ! component type
    type(seq_infodata_type) , intent(in) :: infodata

    !EOP

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_glc_moab) '

    !-------------------------------------------------------------------------------
    ! no BGC fields on glc mesh

  end subroutine seq_diagBGC_glc_moab

  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_ocn_moab - compute global ocn input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global ocn input/output flux diagnostics using MOAB data structures.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_ocn_moab( ocn, infodata, do_o2x, do_x2o)

    type(component_type)    , intent(in)          :: ocn      ! component type
    type(seq_infodata_type) , intent(in)          :: infodata
    logical                 , intent(in),optional :: do_o2x
    logical                 , intent(in),optional :: do_x2o

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp            ! for optional field detection
    integer(in)              :: n,nf,ic,ip        ! generic index
    integer(in)              :: lSize             ! size of mesh
    real(r8)                 :: ca_i,ca_o         ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    real(r8), allocatable    :: ofrac_data(:)
    real(r8), allocatable    :: ifrac_data(:)
    real(r8), allocatable    :: fld_fco2_ocn(:)
    real(r8), allocatable    :: fld_bcphidry(:)
    real(r8), allocatable    :: fld_bcphodry(:)
    real(r8), allocatable    :: fld_bcphiwet(:)
    real(r8), allocatable    :: fld_ocphidry(:)
    real(r8), allocatable    :: fld_ocphodry(:)
    real(r8), allocatable    :: fld_ocphiwet(:)
    real(r8), allocatable    :: fld_algae1(:)
    real(r8), allocatable    :: fld_algae2(:)
    real(r8), allocatable    :: fld_algae3(:)
    real(r8), allocatable    :: fld_dic1(:)
    real(r8), allocatable    :: fld_docr(:)
    real(r8), allocatable    :: fld_doc1(:)
    real(r8), allocatable    :: fld_doc2(:)
    real(r8), allocatable    :: fld_doc3(:)
    logical,save             :: first_time=.true. ! initialization flag
    logical,save             :: flds_c_oi=.false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_ocn_moab) '

    !-------------------------------------------------------------------------------

    if (.not. present(do_o2x) .and. &
        .not. present(do_x2o)) then
       call shr_sys_abort(subName//"ERROR: must input a bundle")
    end if

    lSize = mbGetnCells(mboxid)
    allocate(area_data(lSize), ofrac_data(lSize), ifrac_data(lSize))

    call mbGetCellTagVals(mboxid, afldname,  area_data,  lSize)
    call mbGetCellTagVals(mboxid, ofracname, ofrac_data, lSize)
    call mbGetCellTagVals(mboxid, ifracname, ifrac_data, lSize)

    ip = p_inst

    if (present(do_o2x)) then
       if (first_time) then
          av_tmp => component_get_c2x_cx(ocn)
          index_o2x_Faoo_fco2_ocn = mct_aVect_indexRA(av_tmp,'Faoo_fco2_ocn',perrWith='quiet')
       end if

       if (index_o2x_Faoo_fco2_ocn /= 0) then
          allocate(fld_fco2_ocn(lSize))
          call mbGetCellTagVals(mboxid, 'Faoo_fco2_ocn', fld_fco2_ocn, lSize)
       endif

       ic = c_ocn_or
       do n=1,lSize
          ca_o = area_data(n) * ofrac_data(n)
          ca_i = area_data(n) * ifrac_data(n)
          nf = f_area ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o
          if (index_o2x_Faoo_fco2_ocn /= 0) then
             nf = f_csurf; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) &
                                                + (ca_o+ca_i)*fld_fco2_ocn(n)*CO2toC
          endif
       end do

       if (index_o2x_Faoo_fco2_ocn /= 0) then
          deallocate(fld_fco2_ocn)
       end if
    end if

    if (present(do_x2o)) then
       if (first_time) then
          av_tmp => component_get_x2c_cx(ocn)
          index_x2o_Fioi_algae1 = mct_aVect_indexRA(av_tmp,'Fioi_algae1',perrWith='quiet')
          if ( index_x2o_Fioi_algae1 /= 0 ) flds_c_oi = .true.
          if ( flds_c_oi ) then
             index_x2o_Fioi_algae2 = mct_aVect_indexRA(av_tmp,'Fioi_algae2')
             index_x2o_Fioi_algae3 = mct_aVect_indexRA(av_tmp,'Fioi_algae3')
             index_x2o_Fioi_dic1   = mct_aVect_indexRA(av_tmp,'Fioi_dic1')
             index_x2o_Fioi_docr   = mct_aVect_indexRA(av_tmp,'Fioi_docr')
             index_x2o_Fioi_doc1   = mct_aVect_indexRA(av_tmp,'Fioi_doc1')
             index_x2o_Fioi_doc2   = mct_aVect_indexRA(av_tmp,'Fioi_doc2')
             index_x2o_Fioi_doc3   = mct_aVect_indexRA(av_tmp,'Fioi_doc3')
          end if
       end if

       allocate(fld_bcphidry(lSize), fld_bcphodry(lSize), fld_bcphiwet(lSize))
       allocate(fld_ocphidry(lSize), fld_ocphodry(lSize), fld_ocphiwet(lSize))

       call mbGetCellTagVals(mboxid, 'Faxa_bcphidry', fld_bcphidry, lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_bcphodry', fld_bcphodry, lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_bcphiwet', fld_bcphiwet, lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_ocphidry', fld_ocphidry, lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_ocphodry', fld_ocphodry, lSize)
       call mbGetCellTagVals(mboxid, 'Faxa_ocphiwet', fld_ocphiwet, lSize)

       ic = c_ocn_os
       do n=1,lSize
          ca_o = area_data(n) * ofrac_data(n)
          ca_i = area_data(n) * ifrac_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_o

          nf = f_cblack;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_bcphidry(n) &
                                                       + (ca_o+ca_i)*fld_bcphodry(n) &
                                                       + (ca_o+ca_i)*fld_bcphiwet(n)
          nf = f_corgnc;
          budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_ocphidry(n) &
                                                       + (ca_o+ca_i)*fld_ocphodry(n) &
                                                       + (ca_o+ca_i)*fld_ocphiwet(n)
       end do

       deallocate(fld_bcphidry, fld_bcphodry, fld_bcphiwet)
       deallocate(fld_ocphidry, fld_ocphodry, fld_ocphiwet)

       if ( flds_c_oi ) then
          allocate(fld_dic1(lSize))
          allocate(fld_algae1(lSize), fld_algae2(lSize), fld_algae3(lSize))
          allocate(fld_docr(lSize), fld_doc1(lSize), fld_doc2(lSize), fld_doc3(lSize))

          call mbGetCellTagVals(mboxid, 'Fioi_dic1',   fld_dic1,   lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_algae1', fld_algae1, lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_algae2', fld_algae2, lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_algae3', fld_algae3, lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_docr',   fld_docr,   lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_doc1',   fld_doc1,   lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_doc2',   fld_doc2,   lSize)
          call mbGetCellTagVals(mboxid, 'Fioi_doc3',   fld_doc3,   lSize)

          do n=1,lSize
             ca_o = area_data(n) * ofrac_data(n)
             ca_i = area_data(n) * ifrac_data(n)
             nf = f_ioinor;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + (ca_o+ca_i)*fld_dic1(n) &
                                                         *  shr_const_mwc * 1.0e-06_r8
             nf = f_ioorgn;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ((ca_o+ca_i)*fld_algae1(n) &
                                                          +  (ca_o+ca_i)*fld_algae2(n) &
                                                          +  (ca_o+ca_i)*fld_algae3(n) &
                                                          +  (ca_o+ca_i)*fld_docr(n)   &
                                                          +  (ca_o+ca_i)*fld_doc1(n)   &
                                                          +  (ca_o+ca_i)*fld_doc2(n)   &
                                                          +  (ca_o+ca_i)*fld_doc3(n))  &
                                                          *   shr_const_mwc * 1.0e-06_r8
          end do

          deallocate(fld_dic1)
          deallocate(fld_algae1, fld_algae2, fld_algae3)
          deallocate(fld_docr, fld_doc1, fld_doc2, fld_doc3)
       end if
    end if

    deallocate(area_data, ofrac_data, ifrac_data)

    first_time = .false.

  end subroutine seq_diagBGC_ocn_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_ice_moab - compute global ice input/output flux diagnostics
  !
  ! !DESCRIPTION:
  !     Compute global ice input/output flux diagnostics using MOAB data structures.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_diagBGC_ice_moab( ice, infodata, do_i2x, do_x2i)

    type(component_type)    , intent(in)           :: ice      ! component type
    type(seq_infodata_type) , intent(in)           :: infodata
    logical                 , intent(in), optional :: do_i2x
    logical                 , intent(in), optional :: do_x2i

    !EOP

    !----- local -----
    type(mct_aVect), pointer :: av_tmp            ! for optional field detection
    integer(in)              :: n,ic,nf,ip        ! generic index
    integer(in)              :: lSize             ! size of mesh
    real(r8)                 :: ca_i,ca_o         ! area of a grid cell
    real(r8), allocatable    :: area_data(:)
    real(r8), allocatable    :: lat_data(:)
    real(r8), allocatable    :: ofrac_data(:)
    real(r8), allocatable    :: ifrac_data(:)
    real(r8), allocatable    :: fld_bcphidry(:)
    real(r8), allocatable    :: fld_bcphodry(:)
    real(r8), allocatable    :: fld_bcphiwet(:)
    real(r8), allocatable    :: fld_ocphidry(:)
    real(r8), allocatable    :: fld_ocphodry(:)
    real(r8), allocatable    :: fld_ocphiwet(:)
    real(r8), allocatable    :: fld_dic1(:)
    real(r8), allocatable    :: fld_algae1(:)
    real(r8), allocatable    :: fld_algae2(:)
    real(r8), allocatable    :: fld_algae3(:)
    real(r8), allocatable    :: fld_docr(:)
    real(r8), allocatable    :: fld_doc1(:)
    real(r8), allocatable    :: fld_doc2(:)
    real(r8), allocatable    :: fld_doc3(:)
    logical,save             :: first_time=.true. ! initialization flag
    logical,save             :: flds_c_oi=.false.

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_ice_moab) '

    !-------------------------------------------------------------------------------

    lSize = mbGetnCells(mbixid)
    allocate(area_data(lSize), lat_data(lSize))
    allocate(ofrac_data(lSize), ifrac_data(lSize))

    call mbGetCellTagVals(mbixid, afldname,  area_data,  lSize)
    call mbGetCellTagVals(mbixid, latname,   lat_data,   lSize)
    call mbGetCellTagVals(mbixid, ofracname, ofrac_data, lSize)
    call mbGetCellTagVals(mbixid, ifracname, ifrac_data, lSize)

    ip = p_inst

    if (present(do_i2x)) then
       if (first_time) then
          av_tmp => component_get_c2x_cx(ice)
          index_i2x_Fioi_algae1 = mct_aVect_indexRA(av_tmp,'Fioi_algae1',perrWith='quiet')
          if ( index_i2x_Fioi_algae1 /= 0 ) flds_c_oi = .true.
       endif

       if ( flds_c_oi ) then
          allocate(fld_dic1(lSize))
          allocate(fld_algae1(lSize), fld_algae2(lSize), fld_algae3(lSize))
          allocate(fld_docr(lSize), fld_doc1(lSize), fld_doc2(lSize), fld_doc3(lSize))

          call mbGetCellTagVals(mbixid, 'Fioi_dic1',   fld_dic1,   lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_algae1', fld_algae1, lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_algae2', fld_algae2, lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_algae3', fld_algae3, lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_docr',   fld_docr,   lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_doc1',   fld_doc1,   lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_doc2',   fld_doc2,   lSize)
          call mbGetCellTagVals(mbixid, 'Fioi_doc3',   fld_doc3,   lSize)
       endif

       ! accumulation loop
       do n=1,lSize
          if (lat_data(n) > 0.0_r8) then
             ic = c_inh_ir
          else
             ic = c_ish_ir
          endif
          ca_i = area_data(n) * ifrac_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i
          if ( flds_c_oi ) then
             nf = f_ioinor;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - ca_i*fld_dic1(n) &
                                                          * shr_const_mwc * 1.0e-06_r8
             nf = f_ioorgn;
             budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) - (ca_i*fld_algae1(n) &
                                                          +  ca_i*fld_algae2(n) &
                                                          +  ca_i*fld_algae3(n) &
                                                          +  ca_i*fld_docr(n)   &
                                                          +  ca_i*fld_doc1(n)   &
                                                          +  ca_i*fld_doc2(n)   &
                                                          +  ca_i*fld_doc3(n))  &
                                                          *  shr_const_mwc * 1.0e-06_r8
          endif
       end do

       if ( flds_c_oi ) then
          deallocate(fld_dic1)
          deallocate(fld_algae1, fld_algae2, fld_algae3)
          deallocate(fld_docr, fld_doc1, fld_doc2, fld_doc3)
       end if
    end if

    if (present(do_x2i)) then
       if (first_time) then
          av_tmp => component_get_x2c_cx(ice)
          index_x2i_Faxa_bcphidry = mct_aVect_indexRA(av_tmp,'Faxa_bcphidry')
          index_x2i_Faxa_bcphodry = mct_aVect_indexRA(av_tmp,'Faxa_bcphodry')
          index_x2i_Faxa_bcphiwet = mct_aVect_indexRA(av_tmp,'Faxa_bcphiwet')
          index_x2i_Faxa_ocphidry = mct_aVect_indexRA(av_tmp,'Faxa_ocphidry')
          index_x2i_Faxa_ocphodry = mct_aVect_indexRA(av_tmp,'Faxa_ocphodry')
          index_x2i_Faxa_ocphiwet = mct_aVect_indexRA(av_tmp,'Faxa_ocphiwet')
       end if

       allocate(fld_bcphidry(lSize), fld_bcphodry(lSize), fld_bcphiwet(lSize))
       allocate(fld_ocphidry(lSize), fld_ocphodry(lSize), fld_ocphiwet(lSize))

       call mbGetCellTagVals(mbixid, 'Faxa_bcphidry', fld_bcphidry, lSize)
       call mbGetCellTagVals(mbixid, 'Faxa_bcphodry', fld_bcphodry, lSize)
       call mbGetCellTagVals(mbixid, 'Faxa_bcphiwet', fld_bcphiwet, lSize)
       call mbGetCellTagVals(mbixid, 'Faxa_ocphidry', fld_ocphidry, lSize)
       call mbGetCellTagVals(mbixid, 'Faxa_ocphodry', fld_ocphodry, lSize)
       call mbGetCellTagVals(mbixid, 'Faxa_ocphiwet', fld_ocphiwet, lSize)

       do n=1,lSize
          if (lat_data(n) > 0.0_r8) then
             ic = c_inh_is
          else
             ic = c_ish_is
          endif
          ca_o = area_data(n) * ofrac_data(n)
          ca_i = area_data(n) * ifrac_data(n)
          nf = f_area  ; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i
          nf = f_cblack; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_bcphidry(n) &
                                                                     + ca_i*fld_bcphodry(n) &
                                                                     + ca_i*fld_bcphiwet(n)
          nf = f_corgnc; budg_dataL(nf,ic,ip) = budg_dataL(nf,ic,ip) + ca_i*fld_ocphidry(n) &
                                                                     + ca_i*fld_ocphodry(n) &
                                                                     + ca_i*fld_ocphiwet(n)
       end do

       deallocate(fld_bcphidry, fld_bcphodry, fld_bcphiwet)
       deallocate(fld_ocphidry, fld_ocphodry, fld_ocphiwet)
    end if

    deallocate(area_data, lat_data, ofrac_data, ifrac_data)

    first_time = .false.

  end subroutine seq_diagBGC_ice_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_preprint_moab - prepare global BGC budget diagnostics for printing
  !
  ! !DESCRIPTION:
  !   Prepare global BGC budget diagnostics for printing.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diagBGC_preprint_moab

    implicit none

    !EOP

    !-------------------------------------------------------------------------------
    ! prepare instantaneous budget data for printing
    !-------------------------------------------------------------------------------

    call seq_diagBGC_sum0_moab()
    dataGpr = budg_dataGBGC

    !  old budget normalizations (global area and 1e10 for carbon)
    dataGpr = dataGpr/(4.0_r8*shr_const_pi)
    dataGpr(f_c:f_c_end,:,:) = dataGpr(f_c:f_c_end,:,:) * 1.0e10_r8
    dataGpr = dataGpr/budg_nsBGC

  end subroutine seq_diagBGC_preprint_moab

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_diagBGC_print_moab - print global BGC budget diagnostics
  !
  ! !DESCRIPTION:
  !   Print global BGC budget diagnostics.
  !
  ! !REVISION HISTORY:
  !    2022-apr-12 - J. Wolfe - initial budget implementation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE seq_diagBGC_print_moab(EClock, ip, plev)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) , intent(in) :: EClock
    integer(in)      , intent(in) :: ip
    integer(in)      , intent(in) :: plev        ! print level

    !EOP

    !--- local ---
    integer(in)      :: ic,nf       ! data array indices
    integer(in)      :: ica,icl,icn,ics,ico
    integer(in)      :: icar,icxs,icxr,icas
    integer(in)      :: cdate,sec   ! coded date, seconds
    integer(in)      :: yr,mon,day  ! date
    integer(in)      :: iam         ! pe number
    character(len=40):: str         ! string

    !----- formats -----
    character(*),parameter :: subName = '(seq_diagBGC_print_moab) '
    character(*),parameter :: F00   = "('(seq_diagBGC_print_moab) ',4a)"

    !----- formats -----
    character(*),parameter :: FAH="(4a,i9,i6)"
    character(*),parameter :: FA0= "('    ',12x,6(6x,a8,1x))"
    character(*),parameter :: FA1= "('    ',a12,6f15.8)"
    character(*),parameter :: FA0r="('    ',12x,8(6x,a8,1x))"
    character(*),parameter :: FA1r="('    ',a12,8f15.8)"

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

  end subroutine seq_diagBGC_print_moab

end module seq_diagBGC_moab
