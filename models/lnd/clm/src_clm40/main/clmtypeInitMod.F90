module clmtypeInitMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clmtypeInitMod
!
! !DESCRIPTION:
! Allocate clmtype components and initialize them to signaling NaN.
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use clmtype
  use clm_varpar    , only : maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
                             numpft, ndst, nlevurb, nlevsoi
  use clm_varctl  , only : use_c13, use_cn, use_cndv, use_crop
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initClmtype
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
! Modified by Colette L. Heald (05/06) for VOC emission factors
! 3/17/08 David Lawrence, changed nlevsoi to nlevgrnd where appropriate
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_pft_ecophys_constants
  private :: init_pft_DGVMecophys_constants
  private :: init_pft_pstate_type
  private :: init_pft_epv_type
  private :: init_pft_pdgvstate_type
  private :: init_pft_vstate_type
  private :: init_pft_estate_type
  private :: init_pft_wstate_type
  private :: init_pft_cstate_type
  private :: init_pft_nstate_type
  private :: init_pft_eflux_type
  private :: init_pft_mflux_type
  private :: init_pft_wflux_type
  private :: init_pft_cflux_type
  private :: init_pft_nflux_type
  private :: init_pft_vflux_type
  private :: init_pft_dflux_type
  private :: init_pft_depvd_type
  private :: init_column_pstate_type
  private :: init_column_estate_type
  private :: init_column_wstate_type
  private :: init_column_cstate_type
  private :: init_column_nstate_type
  private :: init_column_eflux_type
  private :: init_column_wflux_type
  private :: init_column_cflux_type
  private :: init_column_nflux_type
  private :: init_landunit_pstate_type
  private :: init_landunit_eflux_type
  private :: init_gridcell_pstate_type
  private :: init_gridcell_efstate_type
  private :: init_gridcell_wflux_type
!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initClmtype
!
! !INTERFACE:
  subroutine initClmtype()
!
! !DESCRIPTION:
! Initialize clmtype components to signaling nan
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
!    *%area, *%wt, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
!    *%ifspecial, *%ityplun, *%itype
!    *%pfti, *%pftf, *%pftn
!    *%coli, *%colf, *%coln
!    *%luni, *%lunf, *%lunn
!
! !USES:
    use abortutils, only : endrun
    use decompMod , only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    character(len=32), parameter :: subname = "initClmtype"
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    call init_pft_type     (begp, endp, pft)
    call init_column_type  (begc, endc, col)
    call init_landunit_type(begl, endl, lun)
    call init_gridcell_type(begg, endg, grc)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    ! pft DGVM-specific ecophysiological constants

    if (use_cndv) then
       call init_pft_DGVMecophys_constants()
    end if

    ! energy balance structures (all levels)

    call init_energy_balance_type(begp, endp, pebal)
    call init_energy_balance_type(begc, endc, cebal)

    ! water balance structures (all levels)

    call init_water_balance_type(begp, endp, pwbal)
    call init_water_balance_type(begc, endc, cwbal)

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(begp, endp, pcbal)
    call init_carbon_balance_type(begc, endc, ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(begp, endp, pnbal)
    call init_nitrogen_balance_type(begc, endc, cnbal)

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(begp, endp, pps)
    call init_pft_pstate_type(begc, endc, pps_a)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, pepv)

    ! pft DGVM state variables at pft level 

    if (use_cndv) then
       call init_pft_pdgvstate_type(begp, endp, pdgvs)
    end if
    call init_pft_vstate_type(begp, endp, pvs)

    ! pft energy state variables at the pft level

    call init_pft_estate_type(begp, endp, pes)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(begp, endp, pws)
    call init_pft_wstate_type(begc, endc, pws_a)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, pcs)
    call init_pft_cstate_type(begc, endc, pcs_a)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_pft_cstate_type(begp, endp, pc13s)
       call init_pft_cstate_type(begc, endc, pc13s_a)
       if (use_crop) then
          call endrun( trim(subname)//" ERROR:: CROP and C13 can NOT be on at the same time" )
       end if
    endif

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, pns)
    call init_pft_nstate_type(begc, endc, pns_a)

    ! pft energy flux variables at pft level 

    call init_pft_eflux_type(begp, endp, pef)

    ! pft momentum flux variables at pft level 

    call init_pft_mflux_type(begp, endp, pmf)

    ! pft water flux variables

    call init_pft_wflux_type(begp, endp, pwf)
    call init_pft_wflux_type(begc, endc, pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, pcf)
    call init_pft_cflux_type(begc, endc, pcf_a)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_pft_cflux_type(begp, endp, pc13f)
       call init_pft_cflux_type(begc, endc, pc13f_a)
    endif

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, pnf)
    call init_pft_nflux_type(begc, endc, pnf_a)

    ! pft VOC flux variables at pft level

    call init_pft_vflux_type(begp, endp, pvf)

    ! gridcell VOC emission factors (heald, 05/06)

    call init_gridcell_efstate_type(begg, endg, gve)

    ! pft dust flux variables at pft level 

    call init_pft_dflux_type(begp, endp, pdf)

    ! pft dry dep velocity variables at pft level 

    call init_pft_depvd_type(begp, endp, pdd)

    ! column physical state variables at column level 

    call init_column_pstate_type(begc, endc, cps)

    ! column energy state variables at column level 


    call init_column_estate_type(begc, endc, ces)

    ! column water state variables at column level 

    call init_column_wstate_type(begc, endc, cws)

    ! column carbon state variables at column level

    call init_column_cstate_type(begc, endc, ccs)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_column_cstate_type(begc, endc, cc13s)
    endif

    ! column nitrogen state variables at column level

    call init_column_nstate_type(begc, endc, cns)

    ! column energy flux variables at column level 

    call init_column_eflux_type(begc, endc, cef)

    ! column water flux variables at column level 

    call init_column_wflux_type(begc, endc, cwf)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, ccf)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_column_cflux_type(begc, endc, cc13f)
    endif

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, cnf)

    ! land unit physical state variables

    call init_landunit_pstate_type(begl, endl, lps)

    ! land unit energy flux variables 

    call init_landunit_eflux_type(begl, endl, lef)

    ! gridcell DGVM variables

    if (use_cndv) then
       call init_gridcell_dgvstate_type(begg, endg, gdgvs)
    end if

    ! gridcell physical state variables


    ! gridcell: water flux variables

    call init_gridcell_wflux_type(begg, endg, gwf)

    ! gridcell: energy flux variables

    call init_gridcell_eflux_type(begg, endg, gef)

    ! gridcell: water state variables

    call init_gridcell_wstate_type(begg, endg, gws)

    ! gridcell: energy state variables

    call init_gridcell_estate_type(begg, endg, ges)

  end subroutine initClmtype

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_type
!
! !INTERFACE:
  subroutine init_pft_type (beg, end, p)
!
! !DESCRIPTION:
! Initialize components of pft_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: p
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pft%gridcell(beg:end),&
             pft%wtgcell(beg:end), &
             pft%landunit(beg:end),&
             pft%wtlunit(beg:end), &
             pft%column(beg:end),  &
             pft%wtcol(beg:end),   &
             pft%itype(beg:end),   &
             pft%mxy(beg:end))

  end subroutine init_pft_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_type
!
! !INTERFACE:
  subroutine init_column_type (beg, end, c)
!
! !DESCRIPTION:
! Initialize components of column_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(column_type), intent(inout):: c
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(col%gridcell(beg:end),&
            col%wtgcell(beg:end), &
            col%landunit(beg:end),&
            col%wtlunit(beg:end), &
            col%pfti(beg:end),    &
            col%pftf(beg:end),    &
            col%npfts(beg:end),   &
            col%itype(beg:end))

  end subroutine init_column_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_type
!
! !INTERFACE:
  subroutine init_landunit_type (beg, end,l)
!
! !DESCRIPTION:
! Initialize components of landunit_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(landunit_type), intent(inout):: l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(lun%gridcell(beg:end), &
            lun%wtgcell(beg:end),  &
            lun%coli(beg:end),     &
            lun%colf(beg:end),     &
            lun%ncolumns(beg:end), &
            lun%pfti(beg:end),     &
            lun%pftf(beg:end),     &
            lun%npfts(beg:end),    & 
            lun%itype(beg:end),    &
            lun%ifspecial(beg:end),& 
            lun%lakpoi(beg:end),   &
            lun%urbpoi(beg:end),   &
            lun%glcmecpoi(beg:end))

   ! These should be moved to landunit physical state
   allocate(lun%canyon_hwr(beg:end),   &
            lun%wtroad_perv(beg:end),  &
            lun%ht_roof(beg:end),      &
            lun%wtlunit_roof(beg:end), &
            lun%wind_hgt_canyon(beg:end),&
            lun%z_0_town(beg:end),     &
            lun%z_d_town(beg:end))

   lun%canyon_hwr(beg:end)     = nan
   lun%wtroad_perv(beg:end)    = nan
   lun%ht_roof(beg:end)        = nan
   lun%wtlunit_roof(beg:end)   = nan
   lun%wind_hgt_canyon(beg:end)= nan
   lun%z_0_town(beg:end)       = nan
   lun%z_d_town(beg:end)       = nan

   lun%glcmecpoi(beg:end) = .false.

  end subroutine init_landunit_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_type
!
! !INTERFACE:
  subroutine init_gridcell_type (beg, end,g)
!
! !DESCRIPTION:
! Initialize components of gridcell_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: g
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(grc%luni(beg:end),      &
            grc%lunf(beg:end),      &
            grc%nlandunits(beg:end),&
            grc%coli(beg:end),      &
            grc%colf(beg:end),      &
            grc%ncolumns(beg:end),  &
            grc%pfti(beg:end),      &
            grc%pftf(beg:end),      &
            grc%npfts(beg:end),     &
            grc%gindex(beg:end),    &
            grc%area(beg:end),      &
            grc%lat(beg:end),       &
            grc%lon(beg:end),       &
            grc%latdeg(beg:end),    &
            grc%londeg(beg:end),    &
            grc%gris_mask(beg:end), &
            grc%gris_area(beg:end), &
            grc%aais_mask(beg:end), &
            grc%aais_area(beg:end))

   grc%gris_mask(beg:end) = nan
   grc%gris_area(beg:end) = nan
   grc%aais_mask(beg:end) = nan
   grc%aais_area(beg:end) = nan

  end subroutine init_gridcell_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_energy_balance_type
!
! !INTERFACE:
  subroutine init_energy_balance_type(beg, end, ebal)
!
! !DESCRIPTION:
! Initialize energy balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(energy_balance_type), intent(inout):: ebal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ebal%errsoi(beg:end))
    allocate(ebal%errseb(beg:end))
    allocate(ebal%errsol(beg:end))
    allocate(ebal%errlon(beg:end))

    ebal%errsoi(beg:end) = nan
    ebal%errseb(beg:end) = nan
    ebal%errsol(beg:end) = nan
    ebal%errlon(beg:end) = nan

  end subroutine init_energy_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_water_balance_type
!
! !INTERFACE:
  subroutine init_water_balance_type(beg, end, wbal)
!
! !DESCRIPTION:
! Initialize water balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(water_balance_type), intent(inout):: wbal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(wbal%begwb(beg:end))
    allocate(wbal%endwb(beg:end))
    allocate(wbal%errh2o(beg:end))

    wbal%begwb(beg:end) = nan
    wbal%endwb(beg:end) = nan
    wbal%errh2o(beg:end) = nan

  end subroutine init_water_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carbon_balance_type
!
! !INTERFACE:
  subroutine init_carbon_balance_type(beg, end, cbal)
!
! !DESCRIPTION:
! Initialize carbon balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(carbon_balance_type), intent(inout):: cbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(cbal%begcb(beg:end))
    allocate(cbal%endcb(beg:end))
    allocate(cbal%errcb(beg:end))

    cbal%begcb(beg:end) = nan
    cbal%endcb(beg:end) = nan
    cbal%errcb(beg:end) = nan

  end subroutine init_carbon_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nitrogen_balance_type
!
! !INTERFACE:
  subroutine init_nitrogen_balance_type(beg, end, nbal)
!
! !DESCRIPTION:
! Initialize nitrogen balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(nitrogen_balance_type), intent(inout):: nbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(nbal%begnb(beg:end))
    allocate(nbal%endnb(beg:end))
    allocate(nbal%errnb(beg:end))

    nbal%begnb(beg:end) = nan
    nbal%endnb(beg:end) = nan
    nbal%errnb(beg:end) = nan

  end subroutine init_nitrogen_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_ecophys_constants
!
! !INTERFACE:
  subroutine init_pft_ecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pftcon%noveg(0:numpft))
    allocate(pftcon%tree(0:numpft))
    allocate(pftcon%smpso(0:numpft)) 
    allocate(pftcon%smpsc(0:numpft)) 
    allocate(pftcon%fnitr(0:numpft))
    allocate(pftcon%foln(0:numpft))
    allocate(pftcon%dleaf(0:numpft))
    allocate(pftcon%c3psn(0:numpft))
    allocate(pftcon%mp(0:numpft))
    allocate(pftcon%qe25(0:numpft))
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft,numrad))
    allocate(pftcon%rhos(0:numpft,numrad))
    allocate(pftcon%taul(0:numpft,numrad))
    allocate(pftcon%taus(0:numpft,numrad))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
    allocate(pftcon%sla(0:numpft))
    allocate(pftcon%slatop(0:numpft))
    allocate(pftcon%dsladlai(0:numpft))
    allocate(pftcon%leafcn(0:numpft))
    allocate(pftcon%flnr(0:numpft))
    allocate(pftcon%woody(0:numpft))
    allocate(pftcon%lflitcn(0:numpft))
    allocate(pftcon%frootcn(0:numpft))
    allocate(pftcon%livewdcn(0:numpft))
    allocate(pftcon%deadwdcn(0:numpft))
    allocate(pftcon%graincn(0:numpft))
    allocate(pftcon%froot_leaf(0:numpft))
    allocate(pftcon%stem_leaf(0:numpft))
    allocate(pftcon%croot_stem(0:numpft))
    allocate(pftcon%flivewd(0:numpft))
    allocate(pftcon%fcur(0:numpft))
    allocate(pftcon%lf_flab(0:numpft))
    allocate(pftcon%lf_fcel(0:numpft))
    allocate(pftcon%lf_flig(0:numpft))
    allocate(pftcon%fr_flab(0:numpft))
    allocate(pftcon%fr_fcel(0:numpft))
    allocate(pftcon%fr_flig(0:numpft))
    allocate(pftcon%leaf_long(0:numpft))
    allocate(pftcon%evergreen(0:numpft))
    allocate(pftcon%stress_decid(0:numpft))
    allocate(pftcon%season_decid(0:numpft))
    allocate(pftcon%resist(0:numpft))
    allocate(pftcon%dwood(0:numpft))

    pftcon%noveg(:) = huge(1)
    pftcon%tree(:) = huge(1)
    pftcon%smpso(:) = nan
    pftcon%smpsc(:) = nan
    pftcon%fnitr(:) = nan
    pftcon%foln(:) = nan
    pftcon%dleaf(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%mp(:) = nan
    pftcon%qe25(:) = nan
    pftcon%xl(:) = nan
    pftcon%rhol(:,:numrad) = nan
    pftcon%rhos(:,:numrad) = nan
    pftcon%taul(:,:numrad) = nan
    pftcon%taus(:,:numrad) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%roota_par(:) = nan
    pftcon%rootb_par(:) = nan
    pftcon%sla(:) = nan
    pftcon%slatop(:) = nan
    pftcon%dsladlai(:) = nan
    pftcon%leafcn(:) = nan
    pftcon%flnr(:) = nan
    pftcon%woody(:) = nan
    pftcon%lflitcn(:) = nan
    pftcon%frootcn(:) = nan
    pftcon%livewdcn(:) = nan
    pftcon%deadwdcn(:) = nan
    pftcon%graincn(:) = nan
    pftcon%froot_leaf(:) = nan
    pftcon%stem_leaf(:) = nan
    pftcon%croot_stem(:) = nan
    pftcon%flivewd(:) = nan
    pftcon%fcur(:) = nan
    pftcon%lf_flab(:) = nan
    pftcon%lf_fcel(:) = nan
    pftcon%lf_flig(:) = nan
    pftcon%fr_flab(:) = nan
    pftcon%fr_fcel(:) = nan
    pftcon%fr_flig(:) = nan
    pftcon%leaf_long(:) = nan
    pftcon%evergreen(:) = nan
    pftcon%stress_decid(:) = nan
    pftcon%season_decid(:) = nan
    pftcon%resist(:) = nan
    pftcon%dwood(:) = nan

  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_DGVMecophys_constants
!
! !INTERFACE:
  subroutine init_pft_DGVMecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(dgv_pftcon%crownarea_max(0:numpft))
    allocate(dgv_pftcon%tcmin(0:numpft))
    allocate(dgv_pftcon%tcmax(0:numpft))
    allocate(dgv_pftcon%gddmin(0:numpft))
    allocate(dgv_pftcon%twmax(0:numpft))
    allocate(dgv_pftcon%reinickerp(0:numpft))
    allocate(dgv_pftcon%allom1(0:numpft))
    allocate(dgv_pftcon%allom2(0:numpft))
    allocate(dgv_pftcon%allom3(0:numpft))

    dgv_pftcon%crownarea_max(:) = nan
    dgv_pftcon%tcmin(:) = nan
    dgv_pftcon%tcmax(:) = nan
    dgv_pftcon%gddmin(:) = nan
    dgv_pftcon%twmax(:) = nan
    dgv_pftcon%reinickerp(:) = nan
    dgv_pftcon%allom1(:) = nan
    dgv_pftcon%allom2(:) = nan
    dgv_pftcon%allom3(:) = nan

  end subroutine init_pft_DGVMecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pstate_type
!
! !INTERFACE:
  subroutine init_pft_pstate_type(beg, end, pps)
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_pstate_type), intent(inout):: pps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pps%frac_veg_nosno(beg:end))
    allocate(pps%frac_veg_nosno_alb(beg:end))
    allocate(pps%emv(beg:end))
    allocate(pps%z0mv(beg:end))
    allocate(pps%z0hv(beg:end))
    allocate(pps%z0qv(beg:end))
    allocate(pps%rootfr(beg:end,1:nlevgrnd))
    allocate(pps%rootr(beg:end,1:nlevgrnd))
    allocate(pps%rresis(beg:end,1:nlevgrnd))
    allocate(pps%dewmx(beg:end))
    allocate(pps%rssun(beg:end))
    allocate(pps%rssha(beg:end))
    allocate(pps%laisun(beg:end))
    allocate(pps%laisha(beg:end))
    allocate(pps%btran(beg:end))
    allocate(pps%fsun(beg:end))
    allocate(pps%tlai(beg:end))
    allocate(pps%tsai(beg:end))
    allocate(pps%elai(beg:end))
    allocate(pps%esai(beg:end))
    allocate(pps%fwet(beg:end))
    allocate(pps%fdry(beg:end))
    allocate(pps%dt_veg(beg:end))
    allocate(pps%htop(beg:end))
    allocate(pps%hbot(beg:end))
    allocate(pps%z0m(beg:end))
    allocate(pps%displa(beg:end))
    allocate(pps%albd(beg:end,1:numrad))
    allocate(pps%albi(beg:end,1:numrad))
    allocate(pps%fabd(beg:end,1:numrad))
    allocate(pps%fabi(beg:end,1:numrad))
    allocate(pps%ftdd(beg:end,1:numrad))
    allocate(pps%ftid(beg:end,1:numrad))
    allocate(pps%ftii(beg:end,1:numrad))
    allocate(pps%u10(beg:end))
    allocate(pps%u10_clm(beg:end))
    allocate(pps%va(beg:end))
    allocate(pps%fv(beg:end))
    allocate(pps%ram1(beg:end))
    if ( crop_prog )then
       allocate(pps%hdidx(beg:end))
       allocate(pps%cumvd(beg:end))
       allocate(pps%htmx(beg:end))
       allocate(pps%vf(beg:end))
       allocate(pps%gddmaturity(beg:end))
       allocate(pps%gdd0(beg:end))
       allocate(pps%gdd8(beg:end))
       allocate(pps%gdd10(beg:end))
       allocate(pps%gdd020(beg:end))
       allocate(pps%gdd820(beg:end))
       allocate(pps%gdd1020(beg:end))
       allocate(pps%gddplant(beg:end))
       allocate(pps%gddtsoi(beg:end))
       allocate(pps%huileaf(beg:end))
       allocate(pps%huigrain(beg:end))
       allocate(pps%aleafi(beg:end))
       allocate(pps%astemi(beg:end))
       allocate(pps%aleaf(beg:end))
       allocate(pps%astem(beg:end))
       allocate(pps%croplive(beg:end))
       allocate(pps%cropplant(beg:end)) !,numpft)) ! make 2-D if using
       allocate(pps%harvdate(beg:end))  !,numpft)) ! crop rotation
       allocate(pps%idop(beg:end))
       allocate(pps%peaklai(beg:end))
    end if
    allocate(pps%vds(beg:end))
    allocate(pps%slasun(beg:end))
    allocate(pps%slasha(beg:end))
    allocate(pps%lncsun(beg:end))
    allocate(pps%lncsha(beg:end))
    allocate(pps%vcmxsun(beg:end))
    allocate(pps%vcmxsha(beg:end))
    allocate(pps%gdir(beg:end))
    allocate(pps%omega(beg:end,1:numrad))
    allocate(pps%eff_kid(beg:end,1:numrad))
    allocate(pps%eff_kii(beg:end,1:numrad))
    allocate(pps%sun_faid(beg:end,1:numrad))
    allocate(pps%sun_faii(beg:end,1:numrad))
    allocate(pps%sha_faid(beg:end,1:numrad))
    allocate(pps%sha_faii(beg:end,1:numrad))
    allocate(pps%forc_hgt_u_pft(beg:end))
    allocate(pps%forc_hgt_t_pft(beg:end))
    allocate(pps%forc_hgt_q_pft(beg:end))
    ! 4/14/05: PET
    ! Adding isotope code
    allocate(pps%cisun(beg:end))
    allocate(pps%cisha(beg:end))
    allocate(pps%alphapsnsun(beg:end))
    allocate(pps%alphapsnsha(beg:end))

    allocate(pps%sandfrac(beg:end))
    allocate(pps%clayfrac(beg:end))
    pps%sandfrac(beg:end) = nan
    pps%clayfrac(beg:end) = nan
    allocate(pps%mlaidiff(beg:end))
    allocate(pps%rb1(beg:end))
    allocate(pps%annlai(12,beg:end))
    pps%mlaidiff(beg:end) = nan
    pps%rb1(beg:end) = nan
    pps%annlai(:,:) = nan
    
    pps%frac_veg_nosno(beg:end) = huge(1)
    pps%frac_veg_nosno_alb(beg:end) = 0
    pps%emv(beg:end) = nan
    pps%z0mv(beg:end) = nan
    pps%z0hv(beg:end) = nan
    pps%z0qv(beg:end) = nan
    pps%rootfr(beg:end,:nlevgrnd) = spval
    pps%rootr (beg:end,:nlevgrnd) = spval
    pps%rresis(beg:end,:nlevgrnd) = spval
    pps%dewmx(beg:end) = nan
    pps%rssun(beg:end) = nan
    pps%rssha(beg:end) = nan
    pps%laisun(beg:end) = nan
    pps%laisha(beg:end) = nan
    pps%btran(beg:end) = spval
    pps%fsun(beg:end) = spval
    pps%tlai(beg:end) = 0._r8
    pps%tsai(beg:end) = 0._r8
    pps%elai(beg:end) = 0._r8
    pps%esai(beg:end) = 0._r8
    pps%fwet(beg:end) = nan
    pps%fdry(beg:end) = nan
    pps%dt_veg(beg:end) = nan
    pps%htop(beg:end) = 0._r8
    pps%hbot(beg:end) = 0._r8
    pps%z0m(beg:end) = nan
    pps%displa(beg:end) = nan
    pps%albd(beg:end,:numrad) = nan
    pps%albi(beg:end,:numrad) = nan
    pps%fabd(beg:end,:numrad) = nan
    pps%fabi(beg:end,:numrad) = nan
    pps%ftdd(beg:end,:numrad) = nan
    pps%ftid(beg:end,:numrad) = nan
    pps%ftii(beg:end,:numrad) = nan
    pps%u10(beg:end) = nan
    pps%u10_clm(beg:end) = nan
    pps%va(beg:end) = nan
    pps%fv(beg:end) = nan
    pps%ram1(beg:end) = nan
    if ( crop_prog )then
       pps%hdidx(beg:end)       = nan
       pps%cumvd(beg:end)       = nan
       pps%htmx(beg:end)        = 0.0_r8
       pps%vf(beg:end)          = 0.0_r8
       pps%gddmaturity(beg:end) = spval
       pps%gdd0(beg:end)        = spval
       pps%gdd8(beg:end)        = spval
       pps%gdd10(beg:end)       = spval
       pps%gdd020(beg:end)      = spval
       pps%gdd820(beg:end)      = spval
       pps%gdd1020(beg:end)     = spval
       pps%gddplant(beg:end)    = spval
       pps%gddtsoi(beg:end)     = spval
       pps%huileaf(beg:end)     = nan
       pps%huigrain(beg:end)    = nan
       pps%aleafi(beg:end)      = nan
       pps%astemi(beg:end)      = nan
       pps%aleaf(beg:end)       = nan
       pps%astem(beg:end)       = nan
       pps%croplive(beg:end)    = .false.
       pps%cropplant(beg:end)   = .false.
       pps%harvdate(beg:end)    = huge(1)
       pps%idop(beg:end)        = huge(1)
       pps%peaklai(beg:end)     = 0
    end if
    pps%vds(beg:end) = nan
    pps%slasun(beg:end) = nan
    pps%slasha(beg:end) = nan
    pps%lncsun(beg:end) = nan
    pps%lncsha(beg:end) = nan
    pps%vcmxsun(beg:end) = nan
    pps%vcmxsha(beg:end) = nan
    pps%gdir(beg:end) = nan
    pps%omega(beg:end,1:numrad) = nan
    pps%eff_kid(beg:end,1:numrad) = nan
    pps%eff_kii(beg:end,1:numrad) = nan
    pps%sun_faid(beg:end,1:numrad) = nan
    pps%sun_faii(beg:end,1:numrad) = nan
    pps%sha_faid(beg:end,1:numrad) = nan
    pps%sha_faii(beg:end,1:numrad) = nan
    pps%forc_hgt_u_pft(beg:end) = nan
    pps%forc_hgt_t_pft(beg:end) = nan
    pps%forc_hgt_q_pft(beg:end) = nan
    ! 4/14/05: PET
    ! Adding isotope code
    pps%cisun(beg:end) = nan
    pps%cisha(beg:end) = nan
    if (use_c13) then
       pps%alphapsnsun(beg:end) = nan
       pps%alphapsnsha(beg:end) = nan
    endif

  end subroutine init_pft_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_epv_type
!
! !INTERFACE:
  subroutine init_pft_epv_type(beg, end, pepv)
!
! !DESCRIPTION:
! Initialize pft ecophysiological variables
!
! !USES:
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_epv_type), intent(inout):: pepv
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pepv%dormant_flag(beg:end))
    allocate(pepv%days_active(beg:end))
    allocate(pepv%onset_flag(beg:end))
    allocate(pepv%onset_counter(beg:end))
    allocate(pepv%onset_gddflag(beg:end))
    allocate(pepv%onset_fdd(beg:end))
    allocate(pepv%onset_gdd(beg:end))
    allocate(pepv%onset_swi(beg:end))
    allocate(pepv%offset_flag(beg:end))
    allocate(pepv%offset_counter(beg:end))
    allocate(pepv%offset_fdd(beg:end))
    allocate(pepv%offset_swi(beg:end))
    allocate(pepv%lgsf(beg:end))
    allocate(pepv%bglfr(beg:end))
    allocate(pepv%bgtr(beg:end))
    allocate(pepv%dayl(beg:end))
    allocate(pepv%prev_dayl(beg:end))
    allocate(pepv%annavg_t2m(beg:end))
    allocate(pepv%tempavg_t2m(beg:end))
    allocate(pepv%gpp(beg:end))
    allocate(pepv%availc(beg:end))
    allocate(pepv%xsmrpool_recover(beg:end))
    allocate(pepv%xsmrpool_c13ratio(beg:end))
    allocate(pepv%alloc_pnow(beg:end))
    allocate(pepv%c_allometry(beg:end))
    allocate(pepv%n_allometry(beg:end))
    allocate(pepv%plant_ndemand(beg:end))
    allocate(pepv%tempsum_potential_gpp(beg:end))
    allocate(pepv%annsum_potential_gpp(beg:end))
    allocate(pepv%tempmax_retransn(beg:end))
    allocate(pepv%annmax_retransn(beg:end))
    allocate(pepv%avail_retransn(beg:end))
    allocate(pepv%plant_nalloc(beg:end))
    allocate(pepv%plant_calloc(beg:end))
    allocate(pepv%excess_cflux(beg:end))
    allocate(pepv%downreg(beg:end))
    allocate(pepv%prev_leafc_to_litter(beg:end))
    allocate(pepv%prev_frootc_to_litter(beg:end))
    allocate(pepv%tempsum_npp(beg:end))
    allocate(pepv%annsum_npp(beg:end))
    allocate(pepv%tempsum_litfall(beg:end))
    allocate(pepv%annsum_litfall(beg:end))
    ! 4/21/05, PET
    ! Adding isotope code
    allocate(pepv%rc13_canair(beg:end))
    allocate(pepv%rc13_psnsun(beg:end))
    allocate(pepv%rc13_psnsha(beg:end))

    pepv%dormant_flag(beg:end) = nan
    pepv%days_active(beg:end) = nan
    pepv%onset_flag(beg:end) = nan
    pepv%onset_counter(beg:end) = nan
    pepv%onset_gddflag(beg:end) = nan
    pepv%onset_fdd(beg:end) = nan
    pepv%onset_gdd(beg:end) = nan
    pepv%onset_swi(beg:end) = nan
    pepv%offset_flag(beg:end) = nan
    pepv%offset_counter(beg:end) = nan
    pepv%offset_fdd(beg:end) = nan
    pepv%offset_swi(beg:end) = nan
    pepv%lgsf(beg:end) = nan
    pepv%bglfr(beg:end) = nan
    pepv%bgtr(beg:end) = nan
    pepv%dayl(beg:end) = nan
    pepv%prev_dayl(beg:end) = nan
    pepv%annavg_t2m(beg:end) = nan
    pepv%tempavg_t2m(beg:end) = nan
    pepv%gpp(beg:end) = nan
    pepv%availc(beg:end) = nan
    pepv%xsmrpool_recover(beg:end) = nan
    if (use_c13) then
       pepv%xsmrpool_c13ratio(beg:end) = nan
    endif
    pepv%alloc_pnow(beg:end) = nan
    pepv%c_allometry(beg:end) = nan
    pepv%n_allometry(beg:end) = nan
    pepv%plant_ndemand(beg:end) = nan
    pepv%tempsum_potential_gpp(beg:end) = nan
    pepv%annsum_potential_gpp(beg:end) = nan
    pepv%tempmax_retransn(beg:end) = nan
    pepv%annmax_retransn(beg:end) = nan
    pepv%avail_retransn(beg:end) = nan
    pepv%plant_nalloc(beg:end) = nan
    pepv%plant_calloc(beg:end) = nan
    pepv%excess_cflux(beg:end) = nan
    pepv%downreg(beg:end) = nan
    pepv%prev_leafc_to_litter(beg:end) = nan
    pepv%prev_frootc_to_litter(beg:end) = nan
    pepv%tempsum_npp(beg:end) = nan
    pepv%annsum_npp(beg:end) = nan
    pepv%tempsum_litfall(beg:end) = nan
    pepv%annsum_litfall(beg:end) = nan
    if (use_c13) then
       ! 4/21/05, PET
       ! Adding isotope code
       pepv%rc13_canair(beg:end) = nan
       pepv%rc13_psnsun(beg:end) = nan
       pepv%rc13_psnsha(beg:end) = nan
    endif
    
  end subroutine init_pft_epv_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pdgvstate_type
!
! !INTERFACE:
  subroutine init_pft_pdgvstate_type(beg, end, pdgvs)
!
! !DESCRIPTION:
! Initialize pft DGVM state variables
!
! !USES:
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dgvstate_type), intent(inout):: pdgvs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdgvs%agddtw(beg:end))
    allocate(pdgvs%agdd(beg:end))
    allocate(pdgvs%t_mo(beg:end))
    allocate(pdgvs%t_mo_min(beg:end))
    allocate(pdgvs%prec365(beg:end))
    allocate(pdgvs%present(beg:end))
    allocate(pdgvs%pftmayexist(beg:end))
    allocate(pdgvs%nind(beg:end))
    allocate(pdgvs%lm_ind(beg:end))
    allocate(pdgvs%lai_ind(beg:end))
    allocate(pdgvs%fpcinc(beg:end))
    allocate(pdgvs%fpcgrid(beg:end))
    allocate(pdgvs%fpcgridold(beg:end))
    allocate(pdgvs%crownarea(beg:end))
    allocate(pdgvs%greffic(beg:end))
    allocate(pdgvs%heatstress(beg:end))

    pdgvs%agddtw(beg:end)           = nan
    pdgvs%agdd(beg:end)             = nan
    pdgvs%t_mo(beg:end)             = nan
    pdgvs%t_mo_min(beg:end)         = nan
    pdgvs%prec365(beg:end)          = nan
    pdgvs%present(beg:end)          = .false.
    pdgvs%pftmayexist(beg:end)      = .true.
    pdgvs%nind(beg:end)             = nan
    pdgvs%lm_ind(beg:end)           = nan
    pdgvs%lai_ind(beg:end)          = nan
    pdgvs%fpcinc(beg:end)           = nan
    pdgvs%fpcgrid(beg:end)          = nan
    pdgvs%fpcgridold(beg:end)       = nan
    pdgvs%crownarea(beg:end)        = nan
    pdgvs%greffic(beg:end)          = nan
    pdgvs%heatstress(beg:end)       = nan

  end subroutine init_pft_pdgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vstate_type
!
! !INTERFACE:
  subroutine init_pft_vstate_type(beg, end, pvs)
!
! !DESCRIPTION:
! Initialize pft VOC variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vstate_type), intent(inout):: pvs
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!------------------------------------------------------------------------

    allocate(pvs%t_veg24 (beg:end))
    allocate(pvs%t_veg240(beg:end))
    allocate(pvs%fsd24   (beg:end))
    allocate(pvs%fsd240  (beg:end))
    allocate(pvs%fsi24   (beg:end))
    allocate(pvs%fsi240  (beg:end))
    allocate(pvs%fsun24  (beg:end))
    allocate(pvs%fsun240 (beg:end))
    allocate(pvs%elai_p  (beg:end))

    pvs%t_veg24 (beg:end)   = spval
    pvs%t_veg240(beg:end)   = spval
    pvs%fsd24   (beg:end)   = spval
    pvs%fsd240  (beg:end)   = spval
    pvs%fsi24   (beg:end)   = spval
    pvs%fsi240  (beg:end)   = spval
    pvs%fsun24  (beg:end)   = spval
    pvs%fsun240 (beg:end)   = spval
    pvs%elai_p  (beg:end)   = spval
  end subroutine init_pft_vstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_estate_type
!
! !INTERFACE:
  subroutine init_pft_estate_type(beg, end, pes)
!
! !DESCRIPTION:
! Initialize pft energy state
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod, only : crop_prog
! !AGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    allocate(pes%t_ref2m(beg:end))
    allocate(pes%t_ref2m_min(beg:end))
    allocate(pes%t_ref2m_max(beg:end))
    allocate(pes%t_ref2m_min_inst(beg:end))
    allocate(pes%t_ref2m_max_inst(beg:end))
    allocate(pes%q_ref2m(beg:end))
    allocate(pes%t_ref2m_u(beg:end))
    allocate(pes%t_ref2m_r(beg:end))
    allocate(pes%t_ref2m_min_u(beg:end))
    allocate(pes%t_ref2m_min_r(beg:end))
    allocate(pes%t_ref2m_max_u(beg:end))
    allocate(pes%t_ref2m_max_r(beg:end))
    allocate(pes%t_ref2m_min_inst_u(beg:end))
    allocate(pes%t_ref2m_min_inst_r(beg:end))
    allocate(pes%t_ref2m_max_inst_u(beg:end))
    allocate(pes%t_ref2m_max_inst_r(beg:end))
    allocate(pes%t10(beg:end))
    if ( crop_prog )then
       allocate(pes%a10tmin(beg:end))
       allocate(pes%a5tmin(beg:end))
    end if
    allocate(pes%rh_ref2m(beg:end))
    allocate(pes%rh_ref2m_u(beg:end))
    allocate(pes%rh_ref2m_r(beg:end))
    allocate(pes%t_veg(beg:end))
    allocate(pes%thm(beg:end))

    pes%t_ref2m(beg:end) = nan
    pes%t_ref2m_min(beg:end) = nan
    pes%t_ref2m_max(beg:end) = nan
    pes%t_ref2m_min_inst(beg:end) = nan
    pes%t_ref2m_max_inst(beg:end) = nan
    pes%q_ref2m(beg:end) = nan
    pes%t_ref2m_u(beg:end) = nan
    pes%t_ref2m_r(beg:end) = nan
    pes%t_ref2m_min_u(beg:end) = nan
    pes%t_ref2m_min_r(beg:end) = nan
    pes%t_ref2m_max_u(beg:end) = nan
    pes%t_ref2m_max_r(beg:end) = nan
    pes%t_ref2m_min_inst_u(beg:end) = nan
    pes%t_ref2m_min_inst_r(beg:end) = nan
    pes%t_ref2m_max_inst_u(beg:end) = nan
    pes%t_ref2m_max_inst_r(beg:end) = nan
    pes%t10(beg:end)                = spval
    if ( crop_prog )then
       pes%a10tmin(beg:end)     = spval
       pes%a5tmin(beg:end)      = spval
    end if
    pes%rh_ref2m(beg:end) = nan
    pes%rh_ref2m_u(beg:end) = nan
    pes%rh_ref2m_r(beg:end) = nan
    pes%t_veg(beg:end) = nan
    pes%thm(beg:end) = nan

  end subroutine init_pft_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wstate_type
!
! !INTERFACE:
  subroutine init_pft_wstate_type(beg, end, pws)
!
! !DESCRIPTION:
! Initialize pft water state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wstate_type), intent(inout):: pws !pft water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pws%h2ocan(beg:end))
    pws%h2ocan(beg:end) = nan

  end subroutine init_pft_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cstate_type
!
! !INTERFACE:
  subroutine init_pft_cstate_type(beg, end, pcs)
!
! !DESCRIPTION:
! Initialize pft carbon state
!
! !USES:
    use surfrdMod, only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cstate_type), intent(inout):: pcs !pft carbon state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pcs%leafc(beg:end))
    allocate(pcs%leafc_storage(beg:end))
    allocate(pcs%leafc_xfer(beg:end))
    allocate(pcs%frootc(beg:end))
    allocate(pcs%frootc_storage(beg:end))
    allocate(pcs%frootc_xfer(beg:end))
    allocate(pcs%livestemc(beg:end))
    allocate(pcs%livestemc_storage(beg:end))
    allocate(pcs%livestemc_xfer(beg:end))
    allocate(pcs%deadstemc(beg:end))
    allocate(pcs%deadstemc_storage(beg:end))
    allocate(pcs%deadstemc_xfer(beg:end))
    allocate(pcs%livecrootc(beg:end))
    allocate(pcs%livecrootc_storage(beg:end))
    allocate(pcs%livecrootc_xfer(beg:end))
    allocate(pcs%deadcrootc(beg:end))
    allocate(pcs%deadcrootc_storage(beg:end))
    allocate(pcs%deadcrootc_xfer(beg:end))
    allocate(pcs%gresp_storage(beg:end))
    allocate(pcs%gresp_xfer(beg:end))
    allocate(pcs%cpool(beg:end))
    allocate(pcs%xsmrpool(beg:end))
    allocate(pcs%pft_ctrunc(beg:end))
    allocate(pcs%dispvegc(beg:end))
    allocate(pcs%storvegc(beg:end))
    allocate(pcs%totvegc(beg:end))
    allocate(pcs%totpftc(beg:end))
    allocate(pcs%leafcmax(beg:end))
    if ( crop_prog )then
       allocate(pcs%grainc(beg:end))
       allocate(pcs%grainc_storage(beg:end))
       allocate(pcs%grainc_xfer(beg:end))
    end if
    allocate(pcs%woodc(beg:end))

    pcs%leafc(beg:end) = nan
    pcs%leafc_storage(beg:end) = nan
    pcs%leafc_xfer(beg:end) = nan
    pcs%frootc(beg:end) = nan
    pcs%frootc_storage(beg:end) = nan
    pcs%frootc_xfer(beg:end) = nan
    pcs%livestemc(beg:end) = nan
    pcs%livestemc_storage(beg:end) = nan
    pcs%livestemc_xfer(beg:end) = nan
    pcs%deadstemc(beg:end) = nan
    pcs%deadstemc_storage(beg:end) = nan
    pcs%deadstemc_xfer(beg:end) = nan
    pcs%livecrootc(beg:end) = nan
    pcs%livecrootc_storage(beg:end) = nan
    pcs%livecrootc_xfer(beg:end) = nan
    pcs%deadcrootc(beg:end) = nan
    pcs%deadcrootc_storage(beg:end) = nan
    pcs%deadcrootc_xfer(beg:end) = nan
    pcs%gresp_storage(beg:end) = nan
    pcs%gresp_xfer(beg:end) = nan
    pcs%cpool(beg:end) = nan
    pcs%xsmrpool(beg:end) = nan
    pcs%pft_ctrunc(beg:end) = nan
    pcs%dispvegc(beg:end) = nan
    pcs%storvegc(beg:end) = nan
    pcs%totvegc(beg:end) = nan
    pcs%totpftc(beg:end) = nan
    pcs%leafcmax(beg:end) = nan
    if ( crop_prog )then
       pcs%grainc(beg:end)         = nan
       pcs%grainc_storage(beg:end) = nan
       pcs%grainc_xfer(beg:end)    = nan
    end if
    pcs%woodc(beg:end) = nan

  end subroutine init_pft_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nstate_type
!
! !INTERFACE:
  subroutine init_pft_nstate_type(beg, end, pns)
!
! !DESCRIPTION:
! Initialize pft nitrogen state
!
! !USES:
    use surfrdMod, only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nstate_type), intent(inout):: pns !pft nitrogen state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    if ( crop_prog )then
       allocate(pns%grainn(beg:end))
       allocate(pns%grainn_storage(beg:end))
       allocate(pns%grainn_xfer(beg:end))
    end if
    allocate(pns%leafn(beg:end))
    allocate(pns%leafn_storage(beg:end))
    allocate(pns%leafn_xfer(beg:end))
    allocate(pns%frootn(beg:end))
    allocate(pns%frootn_storage(beg:end))
    allocate(pns%frootn_xfer(beg:end))
    allocate(pns%livestemn(beg:end))
    allocate(pns%livestemn_storage(beg:end))
    allocate(pns%livestemn_xfer(beg:end))
    allocate(pns%deadstemn(beg:end))
    allocate(pns%deadstemn_storage(beg:end))
    allocate(pns%deadstemn_xfer(beg:end))
    allocate(pns%livecrootn(beg:end))
    allocate(pns%livecrootn_storage(beg:end))
    allocate(pns%livecrootn_xfer(beg:end))
    allocate(pns%deadcrootn(beg:end))
    allocate(pns%deadcrootn_storage(beg:end))
    allocate(pns%deadcrootn_xfer(beg:end))
    allocate(pns%retransn(beg:end))
    allocate(pns%npool(beg:end))
    allocate(pns%pft_ntrunc(beg:end))
    allocate(pns%dispvegn(beg:end))
    allocate(pns%storvegn(beg:end))
    allocate(pns%totvegn(beg:end))
    allocate(pns%totpftn(beg:end))

    if ( crop_prog )then
       pns%grainn(beg:end)         = nan
       pns%grainn_storage(beg:end) = nan
       pns%grainn_xfer(beg:end)    = nan
    end if
    pns%leafn(beg:end) = nan
    pns%leafn_storage(beg:end) = nan
    pns%leafn_xfer(beg:end) = nan
    pns%frootn(beg:end) = nan
    pns%frootn_storage(beg:end) = nan
    pns%frootn_xfer(beg:end) = nan
    pns%livestemn(beg:end) = nan
    pns%livestemn_storage(beg:end) = nan
    pns%livestemn_xfer(beg:end) = nan
    pns%deadstemn(beg:end) = nan
    pns%deadstemn_storage(beg:end) = nan
    pns%deadstemn_xfer(beg:end) = nan
    pns%livecrootn(beg:end) = nan
    pns%livecrootn_storage(beg:end) = nan
    pns%livecrootn_xfer(beg:end) = nan
    pns%deadcrootn(beg:end) = nan
    pns%deadcrootn_storage(beg:end) = nan
    pns%deadcrootn_xfer(beg:end) = nan
    pns%retransn(beg:end) = nan
    pns%npool(beg:end) = nan
    pns%pft_ntrunc(beg:end) = nan
    pns%dispvegn(beg:end) = nan
    pns%storvegn(beg:end) = nan
    pns%totvegn(beg:end) = nan
    pns%totpftn(beg:end) = nan

  end subroutine init_pft_nstate_type
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_eflux_type
!
! !INTERFACE:
  subroutine init_pft_eflux_type(beg, end, pef)
!
! !DESCRIPTION:
! Initialize pft energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_eflux_type), intent(inout):: pef
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pef%sabg(beg:end))
    allocate(pef%sabv(beg:end))
    allocate(pef%fsa(beg:end))
    allocate(pef%fsa_u(beg:end))
    allocate(pef%fsa_r(beg:end))
    allocate(pef%fsr(beg:end))
    allocate(pef%parsun(beg:end))
    allocate(pef%parsha(beg:end))
    allocate(pef%dlrad(beg:end))
    allocate(pef%ulrad(beg:end))
    allocate(pef%eflx_lh_tot(beg:end))
    allocate(pef%eflx_lh_tot_u(beg:end))
    allocate(pef%eflx_lh_tot_r(beg:end))
    allocate(pef%eflx_lh_grnd(beg:end))
    allocate(pef%eflx_soil_grnd(beg:end))
    allocate(pef%eflx_soil_grnd_u(beg:end))
    allocate(pef%eflx_soil_grnd_r(beg:end))
    allocate(pef%eflx_sh_tot(beg:end))
    allocate(pef%eflx_sh_tot_u(beg:end))
    allocate(pef%eflx_sh_tot_r(beg:end))
    allocate(pef%eflx_sh_grnd(beg:end))
    allocate(pef%eflx_sh_veg(beg:end))
    allocate(pef%eflx_lh_vege(beg:end))
    allocate(pef%eflx_lh_vegt(beg:end))
    allocate(pef%eflx_wasteheat_pft(beg:end))
    allocate(pef%eflx_heat_from_ac_pft(beg:end))
    allocate(pef%eflx_traffic_pft(beg:end))
    allocate(pef%eflx_anthro(beg:end))
    allocate(pef%cgrnd(beg:end))
    allocate(pef%cgrndl(beg:end))
    allocate(pef%cgrnds(beg:end))
    allocate(pef%eflx_gnet(beg:end))
    allocate(pef%dgnetdT(beg:end))
    allocate(pef%eflx_lwrad_out(beg:end))
    allocate(pef%eflx_lwrad_net(beg:end))
    allocate(pef%eflx_lwrad_net_u(beg:end))
    allocate(pef%eflx_lwrad_net_r(beg:end))
    allocate(pef%netrad(beg:end))
    allocate(pef%fsds_vis_d(beg:end))
    allocate(pef%fsds_nir_d(beg:end))
    allocate(pef%fsds_vis_i(beg:end))
    allocate(pef%fsds_nir_i(beg:end))
    allocate(pef%fsr_vis_d(beg:end))
    allocate(pef%fsr_nir_d(beg:end))
    allocate(pef%fsr_vis_i(beg:end))
    allocate(pef%fsr_nir_i(beg:end))
    allocate(pef%fsds_vis_d_ln(beg:end))
    allocate(pef%fsds_nir_d_ln(beg:end))
    allocate(pef%fsr_vis_d_ln(beg:end))
    allocate(pef%fsr_nir_d_ln(beg:end))
    allocate(pef%sun_add(beg:end,1:numrad))
    allocate(pef%tot_aid(beg:end,1:numrad))
    allocate(pef%sun_aid(beg:end,1:numrad))
    allocate(pef%sun_aii(beg:end,1:numrad))
    allocate(pef%sha_aid(beg:end,1:numrad))
    allocate(pef%sha_aii(beg:end,1:numrad))
    allocate(pef%sun_atot(beg:end,1:numrad))
    allocate(pef%sha_atot(beg:end,1:numrad))
    allocate(pef%sun_alf(beg:end,1:numrad))
    allocate(pef%sha_alf(beg:end,1:numrad))
    allocate(pef%sun_aperlai(beg:end,1:numrad))
    allocate(pef%sha_aperlai(beg:end,1:numrad))
    allocate(pef%sabg_lyr(beg:end,-nlevsno+1:1))
    allocate(pef%sfc_frc_aer(beg:end))
    allocate(pef%sfc_frc_bc(beg:end))
    allocate(pef%sfc_frc_oc(beg:end))
    allocate(pef%sfc_frc_dst(beg:end))
    allocate(pef%sfc_frc_aer_sno(beg:end))
    allocate(pef%sfc_frc_bc_sno(beg:end))
    allocate(pef%sfc_frc_oc_sno(beg:end))
    allocate(pef%sfc_frc_dst_sno(beg:end))
    allocate(pef%fsr_sno_vd(beg:end))
    allocate(pef%fsr_sno_nd(beg:end))
    allocate(pef%fsr_sno_vi(beg:end))
    allocate(pef%fsr_sno_ni(beg:end))
    allocate(pef%fsds_sno_vd(beg:end))
    allocate(pef%fsds_sno_nd(beg:end))
    allocate(pef%fsds_sno_vi(beg:end))
    allocate(pef%fsds_sno_ni(beg:end))

    pef%sabg(beg:end) = nan
    pef%sabv(beg:end) = nan
    pef%fsa(beg:end) = nan
    pef%fsa_u(beg:end) = nan
    pef%fsa_r(beg:end) = nan
    pef%fsr(beg:end) = nan
    pef%parsun(beg:end) = nan
    pef%parsha(beg:end) = nan
    pef%dlrad(beg:end) = nan
    pef%ulrad(beg:end) = nan
    pef%eflx_lh_tot(beg:end) = nan
    pef%eflx_lh_tot_u(beg:end) = nan
    pef%eflx_lh_tot_r(beg:end) = nan
    pef%eflx_lh_grnd(beg:end) = nan
    pef%eflx_soil_grnd(beg:end) = nan
    pef%eflx_soil_grnd_u(beg:end) = nan
    pef%eflx_soil_grnd_r(beg:end) = nan
    pef%eflx_sh_tot(beg:end) = nan
    pef%eflx_sh_tot_u(beg:end) = nan
    pef%eflx_sh_tot_r(beg:end) = nan
    pef%eflx_sh_grnd(beg:end) = nan
    pef%eflx_sh_veg(beg:end) = nan
    pef%eflx_lh_vege(beg:end) = nan
    pef%eflx_lh_vegt(beg:end) = nan
    pef%eflx_wasteheat_pft(beg:end) = nan
    pef%eflx_heat_from_ac_pft(beg:end) = nan
    pef%eflx_traffic_pft(beg:end) = nan
    pef%eflx_anthro(beg:end) = nan
    pef%cgrnd(beg:end) = nan
    pef%cgrndl(beg:end) = nan
    pef%cgrnds(beg:end) = nan
    pef%eflx_gnet(beg:end) = nan
    pef%dgnetdT(beg:end) = nan
    pef%eflx_lwrad_out(beg:end) = nan
    pef%eflx_lwrad_net(beg:end) = nan
    pef%eflx_lwrad_net_u(beg:end) = nan
    pef%eflx_lwrad_net_r(beg:end) = nan
    pef%netrad(beg:end) = nan
    pef%fsds_vis_d(beg:end) = nan
    pef%fsds_nir_d(beg:end) = nan
    pef%fsds_vis_i(beg:end) = nan
    pef%fsds_nir_i(beg:end) = nan
    pef%fsr_vis_d(beg:end) = nan
    pef%fsr_nir_d(beg:end) = nan
    pef%fsr_vis_i(beg:end) = nan
    pef%fsr_nir_i(beg:end) = nan
    pef%fsds_vis_d_ln(beg:end) = nan
    pef%fsds_nir_d_ln(beg:end) = nan
    pef%fsr_vis_d_ln(beg:end) = nan
    pef%fsr_nir_d_ln(beg:end) = nan
    pef%sun_add(beg:end,1:numrad) = nan
    pef%tot_aid(beg:end,1:numrad) = nan
    pef%sun_aid(beg:end,1:numrad) = nan
    pef%sun_aii(beg:end,1:numrad) = nan
    pef%sha_aid(beg:end,1:numrad) = nan
    pef%sha_aii(beg:end,1:numrad) = nan
    pef%sun_atot(beg:end,1:numrad) = nan
    pef%sha_atot(beg:end,1:numrad) = nan
    pef%sun_alf(beg:end,1:numrad) = nan
    pef%sha_alf(beg:end,1:numrad) = nan
    pef%sun_aperlai(beg:end,1:numrad) = nan
    pef%sha_aperlai(beg:end,1:numrad) = nan
    pef%sabg_lyr(beg:end,-nlevsno+1:1) = nan
    pef%sfc_frc_aer(beg:end) = nan
    pef%sfc_frc_bc(beg:end) = nan
    pef%sfc_frc_oc(beg:end) = nan
    pef%sfc_frc_dst(beg:end) = nan
    pef%sfc_frc_aer_sno(beg:end) = nan
    pef%sfc_frc_bc_sno(beg:end) = nan
    pef%sfc_frc_oc_sno(beg:end) = nan
    pef%sfc_frc_dst_sno(beg:end) = nan
    pef%fsr_sno_vd(beg:end) = nan
    pef%fsr_sno_nd(beg:end) = nan
    pef%fsr_sno_vi(beg:end) = nan
    pef%fsr_sno_ni(beg:end) = nan
    pef%fsds_sno_vd(beg:end) = nan
    pef%fsds_sno_nd(beg:end) = nan
    pef%fsds_sno_vi(beg:end) = nan
    pef%fsds_sno_ni(beg:end) = nan
  end subroutine init_pft_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_mflux_type
!
! !INTERFACE:
  subroutine init_pft_mflux_type(beg, end, pmf)
!
! !DESCRIPTION:
! Initialize pft momentum flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_mflux_type), intent(inout) :: pmf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pmf%taux(beg:end))
    allocate(pmf%tauy(beg:end))

    pmf%taux(beg:end) = nan
    pmf%tauy(beg:end) = nan

  end subroutine init_pft_mflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wflux_type
!
! !INTERFACE:
  subroutine init_pft_wflux_type(beg, end, pwf)
!
! !DESCRIPTION:
! Initialize pft water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wflux_type), intent(inout) :: pwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pwf%qflx_prec_intr(beg:end))
    allocate(pwf%qflx_prec_grnd(beg:end))
    allocate(pwf%qflx_rain_grnd(beg:end))
    allocate(pwf%qflx_snow_grnd(beg:end))
    allocate(pwf%qflx_snwcp_liq(beg:end))
    allocate(pwf%qflx_snwcp_ice(beg:end))
    allocate(pwf%qflx_evap_veg(beg:end))
    allocate(pwf%qflx_tran_veg(beg:end))
    allocate(pwf%qflx_evap_can(beg:end))
    allocate(pwf%qflx_evap_soi(beg:end))
    allocate(pwf%qflx_evap_tot(beg:end))
    allocate(pwf%qflx_evap_grnd(beg:end))
    allocate(pwf%qflx_dew_grnd(beg:end))
    allocate(pwf%qflx_sub_snow(beg:end))
    allocate(pwf%qflx_dew_snow(beg:end))

    pwf%qflx_prec_intr(beg:end) = nan
    pwf%qflx_prec_grnd(beg:end) = nan
    pwf%qflx_rain_grnd(beg:end) = nan
    pwf%qflx_snow_grnd(beg:end) = nan
    pwf%qflx_snwcp_liq(beg:end) = nan
    pwf%qflx_snwcp_ice(beg:end) = nan
    pwf%qflx_evap_veg(beg:end) = nan
    pwf%qflx_tran_veg(beg:end) = nan
    pwf%qflx_evap_can(beg:end) = nan
    pwf%qflx_evap_soi(beg:end) = nan
    pwf%qflx_evap_tot(beg:end) = nan
    pwf%qflx_evap_grnd(beg:end) = nan
    pwf%qflx_dew_grnd(beg:end) = nan
    pwf%qflx_sub_snow(beg:end) = nan
    pwf%qflx_dew_snow(beg:end) = nan

  end subroutine init_pft_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cflux_type
!
! !INTERFACE:
  subroutine init_pft_cflux_type(beg, end, pcf)
!
! !DESCRIPTION:
! Initialize pft carbon flux variables
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cflux_type), intent(inout) :: pcf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pcf%psnsun(beg:end))
    allocate(pcf%psnsha(beg:end))
    allocate(pcf%fpsn(beg:end))
    allocate(pcf%fco2(beg:end))

    allocate(pcf%m_leafc_to_litter(beg:end))
    allocate(pcf%m_frootc_to_litter(beg:end))
    allocate(pcf%m_leafc_storage_to_litter(beg:end))
    allocate(pcf%m_frootc_storage_to_litter(beg:end))
    allocate(pcf%m_livestemc_storage_to_litter(beg:end))
    allocate(pcf%m_deadstemc_storage_to_litter(beg:end))
    allocate(pcf%m_livecrootc_storage_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_litter(beg:end))
    allocate(pcf%m_leafc_xfer_to_litter(beg:end))
    allocate(pcf%m_frootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_litter(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_to_litter(beg:end))
    allocate(pcf%m_deadstemc_to_litter(beg:end))
    allocate(pcf%m_livecrootc_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_to_litter(beg:end))
    allocate(pcf%m_gresp_storage_to_litter(beg:end))
    allocate(pcf%m_gresp_xfer_to_litter(beg:end))
    allocate(pcf%hrv_leafc_to_litter(beg:end))             
    allocate(pcf%hrv_leafc_storage_to_litter(beg:end))     
    allocate(pcf%hrv_leafc_xfer_to_litter(beg:end))        
    allocate(pcf%hrv_frootc_to_litter(beg:end))            
    allocate(pcf%hrv_frootc_storage_to_litter(beg:end))    
    allocate(pcf%hrv_frootc_xfer_to_litter(beg:end))       
    allocate(pcf%hrv_livestemc_to_litter(beg:end))         
    allocate(pcf%hrv_livestemc_storage_to_litter(beg:end)) 
    allocate(pcf%hrv_livestemc_xfer_to_litter(beg:end))    
    allocate(pcf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(pcf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(pcf%hrv_deadstemc_storage_to_litter(beg:end)) 
    allocate(pcf%hrv_deadstemc_xfer_to_litter(beg:end))    
    allocate(pcf%hrv_livecrootc_to_litter(beg:end))        
    allocate(pcf%hrv_livecrootc_storage_to_litter(beg:end))
    allocate(pcf%hrv_livecrootc_xfer_to_litter(beg:end))   
    allocate(pcf%hrv_deadcrootc_to_litter(beg:end))        
    allocate(pcf%hrv_deadcrootc_storage_to_litter(beg:end))
    allocate(pcf%hrv_deadcrootc_xfer_to_litter(beg:end))   
    allocate(pcf%hrv_gresp_storage_to_litter(beg:end))     
    allocate(pcf%hrv_gresp_xfer_to_litter(beg:end))        
    allocate(pcf%hrv_xsmrpool_to_atm(beg:end))                 
    allocate(pcf%m_leafc_to_fire(beg:end))
    allocate(pcf%m_frootc_to_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_fire(beg:end))
    allocate(pcf%leafc_xfer_to_leafc(beg:end))
    allocate(pcf%frootc_xfer_to_frootc(beg:end))
    allocate(pcf%livestemc_xfer_to_livestemc(beg:end))
    allocate(pcf%deadstemc_xfer_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_xfer_to_livecrootc(beg:end))
    allocate(pcf%deadcrootc_xfer_to_deadcrootc(beg:end))
    allocate(pcf%leafc_to_litter(beg:end))
    allocate(pcf%frootc_to_litter(beg:end))
    allocate(pcf%leaf_mr(beg:end))
    allocate(pcf%froot_mr(beg:end))
    allocate(pcf%livestem_mr(beg:end))
    allocate(pcf%livecroot_mr(beg:end))
    allocate(pcf%leaf_curmr(beg:end))
    allocate(pcf%froot_curmr(beg:end))
    allocate(pcf%livestem_curmr(beg:end))
    allocate(pcf%livecroot_curmr(beg:end))
    allocate(pcf%leaf_xsmr(beg:end))
    allocate(pcf%froot_xsmr(beg:end))
    allocate(pcf%livestem_xsmr(beg:end))
    allocate(pcf%livecroot_xsmr(beg:end))
    allocate(pcf%psnsun_to_cpool(beg:end))
    allocate(pcf%psnshade_to_cpool(beg:end))
    allocate(pcf%cpool_to_xsmrpool(beg:end))
    allocate(pcf%cpool_to_leafc(beg:end))
    allocate(pcf%cpool_to_leafc_storage(beg:end))
    allocate(pcf%cpool_to_frootc(beg:end))
    allocate(pcf%cpool_to_frootc_storage(beg:end))
    allocate(pcf%cpool_to_livestemc(beg:end))
    allocate(pcf%cpool_to_livestemc_storage(beg:end))
    allocate(pcf%cpool_to_deadstemc(beg:end))
    allocate(pcf%cpool_to_deadstemc_storage(beg:end))
    allocate(pcf%cpool_to_livecrootc(beg:end))
    allocate(pcf%cpool_to_livecrootc_storage(beg:end))
    allocate(pcf%cpool_to_deadcrootc(beg:end))
    allocate(pcf%cpool_to_deadcrootc_storage(beg:end))
    allocate(pcf%cpool_to_gresp_storage(beg:end))
    allocate(pcf%cpool_leaf_gr(beg:end))
    allocate(pcf%cpool_leaf_storage_gr(beg:end))
    allocate(pcf%transfer_leaf_gr(beg:end))
    allocate(pcf%cpool_froot_gr(beg:end))
    allocate(pcf%cpool_froot_storage_gr(beg:end))
    allocate(pcf%transfer_froot_gr(beg:end))
    allocate(pcf%cpool_livestem_gr(beg:end))
    allocate(pcf%cpool_livestem_storage_gr(beg:end))
    allocate(pcf%transfer_livestem_gr(beg:end))
    allocate(pcf%cpool_deadstem_gr(beg:end))
    allocate(pcf%cpool_deadstem_storage_gr(beg:end))
    allocate(pcf%transfer_deadstem_gr(beg:end))
    allocate(pcf%cpool_livecroot_gr(beg:end))
    allocate(pcf%cpool_livecroot_storage_gr(beg:end))
    allocate(pcf%transfer_livecroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_storage_gr(beg:end))
    allocate(pcf%transfer_deadcroot_gr(beg:end))
    allocate(pcf%leafc_storage_to_xfer(beg:end))
    allocate(pcf%frootc_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_storage_to_xfer(beg:end))
    allocate(pcf%deadstemc_storage_to_xfer(beg:end))
    allocate(pcf%livecrootc_storage_to_xfer(beg:end))
    allocate(pcf%deadcrootc_storage_to_xfer(beg:end))
    allocate(pcf%gresp_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_to_deadcrootc(beg:end))
    allocate(pcf%gpp(beg:end))
    allocate(pcf%mr(beg:end))
    allocate(pcf%current_gr(beg:end))
    allocate(pcf%transfer_gr(beg:end))
    allocate(pcf%storage_gr(beg:end))
    allocate(pcf%gr(beg:end))
    allocate(pcf%ar(beg:end))
    allocate(pcf%rr(beg:end))
    allocate(pcf%npp(beg:end))
    allocate(pcf%agnpp(beg:end))
    allocate(pcf%bgnpp(beg:end))
    allocate(pcf%litfall(beg:end))
    allocate(pcf%vegfire(beg:end))
    allocate(pcf%wood_harvestc(beg:end))
    allocate(pcf%pft_cinputs(beg:end))
    allocate(pcf%pft_coutputs(beg:end))
    allocate(pcf%pft_fire_closs(beg:end))
    if ( crop_prog )then
       allocate(pcf%xsmrpool_to_atm(beg:end))
       allocate(pcf%grainc_xfer_to_grainc(beg:end))
       allocate(pcf%livestemc_to_litter(beg:end))
       allocate(pcf%grainc_to_food(beg:end))
       allocate(pcf%cpool_to_grainc(beg:end))
       allocate(pcf%cpool_to_grainc_storage(beg:end))
       allocate(pcf%cpool_grain_gr(beg:end))
       allocate(pcf%cpool_grain_storage_gr(beg:end))
       allocate(pcf%transfer_grain_gr(beg:end))
       allocate(pcf%grainc_storage_to_xfer(beg:end))
    end if
    if (use_cn) then
       allocate(pcf%frootc_alloc(beg:end))
       allocate(pcf%frootc_loss(beg:end))
       allocate(pcf%leafc_alloc(beg:end))
       allocate(pcf%leafc_loss(beg:end))
       allocate(pcf%woodc_alloc(beg:end))
       allocate(pcf%woodc_loss(beg:end))
    end if

    pcf%psnsun(beg:end) = nan
    pcf%psnsha(beg:end) = nan
    pcf%fpsn(beg:end) = spval
    pcf%fco2(beg:end) = 0._r8

    pcf%m_leafc_to_litter(beg:end) = nan
    pcf%m_frootc_to_litter(beg:end) = nan
    pcf%m_leafc_storage_to_litter(beg:end) = nan
    pcf%m_frootc_storage_to_litter(beg:end) = nan
    pcf%m_livestemc_storage_to_litter(beg:end) = nan
    pcf%m_deadstemc_storage_to_litter(beg:end) = nan
    pcf%m_livecrootc_storage_to_litter(beg:end) = nan
    pcf%m_deadcrootc_storage_to_litter(beg:end) = nan
    pcf%m_leafc_xfer_to_litter(beg:end) = nan
    pcf%m_frootc_xfer_to_litter(beg:end) = nan
    pcf%m_livestemc_xfer_to_litter(beg:end) = nan
    pcf%m_deadstemc_xfer_to_litter(beg:end) = nan
    pcf%m_livecrootc_xfer_to_litter(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_litter(beg:end) = nan
    pcf%m_livestemc_to_litter(beg:end) = nan
    pcf%m_deadstemc_to_litter(beg:end) = nan
    pcf%m_livecrootc_to_litter(beg:end) = nan
    pcf%m_deadcrootc_to_litter(beg:end) = nan
    pcf%m_gresp_storage_to_litter(beg:end) = nan
    pcf%m_gresp_xfer_to_litter(beg:end) = nan
    pcf%hrv_leafc_to_litter(beg:end) = nan             
    pcf%hrv_leafc_storage_to_litter(beg:end) = nan     
    pcf%hrv_leafc_xfer_to_litter(beg:end) = nan        
    pcf%hrv_frootc_to_litter(beg:end) = nan            
    pcf%hrv_frootc_storage_to_litter(beg:end) = nan    
    pcf%hrv_frootc_xfer_to_litter(beg:end) = nan       
    pcf%hrv_livestemc_to_litter(beg:end) = nan         
    pcf%hrv_livestemc_storage_to_litter(beg:end) = nan 
    pcf%hrv_livestemc_xfer_to_litter(beg:end) = nan    
    pcf%hrv_deadstemc_to_prod10c(beg:end) = nan        
    pcf%hrv_deadstemc_to_prod100c(beg:end) = nan       
    pcf%hrv_deadstemc_storage_to_litter(beg:end) = nan 
    pcf%hrv_deadstemc_xfer_to_litter(beg:end) = nan    
    pcf%hrv_livecrootc_to_litter(beg:end) = nan        
    pcf%hrv_livecrootc_storage_to_litter(beg:end) = nan
    pcf%hrv_livecrootc_xfer_to_litter(beg:end) = nan   
    pcf%hrv_deadcrootc_to_litter(beg:end) = nan        
    pcf%hrv_deadcrootc_storage_to_litter(beg:end) = nan
    pcf%hrv_deadcrootc_xfer_to_litter(beg:end) = nan   
    pcf%hrv_gresp_storage_to_litter(beg:end) = nan     
    pcf%hrv_gresp_xfer_to_litter(beg:end) = nan        
    pcf%hrv_xsmrpool_to_atm(beg:end) = nan                 
    pcf%m_leafc_to_fire(beg:end) = nan
    pcf%m_frootc_to_fire(beg:end) = nan
    pcf%m_leafc_storage_to_fire(beg:end) = nan
    pcf%m_frootc_storage_to_fire(beg:end) = nan
    pcf%m_livestemc_storage_to_fire(beg:end) = nan
    pcf%m_deadstemc_storage_to_fire(beg:end) = nan
    pcf%m_livecrootc_storage_to_fire(beg:end) = nan
    pcf%m_deadcrootc_storage_to_fire(beg:end) = nan
    pcf%m_leafc_xfer_to_fire(beg:end) = nan
    pcf%m_frootc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_xfer_to_fire(beg:end) = nan
    pcf%m_deadstemc_xfer_to_fire(beg:end) = nan
    pcf%m_livecrootc_xfer_to_fire(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_litter_fire(beg:end) = nan
    pcf%m_gresp_storage_to_fire(beg:end) = nan
    pcf%m_gresp_xfer_to_fire(beg:end) = nan
    pcf%leafc_xfer_to_leafc(beg:end) = nan
    pcf%frootc_xfer_to_frootc(beg:end) = nan
    pcf%livestemc_xfer_to_livestemc(beg:end) = nan
    pcf%deadstemc_xfer_to_deadstemc(beg:end) = nan
    pcf%livecrootc_xfer_to_livecrootc(beg:end) = nan
    pcf%deadcrootc_xfer_to_deadcrootc(beg:end) = nan
    pcf%leafc_to_litter(beg:end) = nan
    pcf%frootc_to_litter(beg:end) = nan
    pcf%leaf_mr(beg:end) = nan
    pcf%froot_mr(beg:end) = nan
    pcf%livestem_mr(beg:end) = nan
    pcf%livecroot_mr(beg:end) = nan
    pcf%leaf_curmr(beg:end) = nan
    pcf%froot_curmr(beg:end) = nan
    pcf%livestem_curmr(beg:end) = nan
    pcf%livecroot_curmr(beg:end) = nan
    pcf%leaf_xsmr(beg:end) = nan
    pcf%froot_xsmr(beg:end) = nan
    pcf%livestem_xsmr(beg:end) = nan
    pcf%livecroot_xsmr(beg:end) = nan
    pcf%psnsun_to_cpool(beg:end) = nan
    pcf%psnshade_to_cpool(beg:end) = nan
    pcf%cpool_to_xsmrpool(beg:end) = nan
    pcf%cpool_to_leafc(beg:end) = nan
    pcf%cpool_to_leafc_storage(beg:end) = nan
    pcf%cpool_to_frootc(beg:end) = nan
    pcf%cpool_to_frootc_storage(beg:end) = nan
    pcf%cpool_to_livestemc(beg:end) = nan
    pcf%cpool_to_livestemc_storage(beg:end) = nan
    pcf%cpool_to_deadstemc(beg:end) = nan
    pcf%cpool_to_deadstemc_storage(beg:end) = nan
    pcf%cpool_to_livecrootc(beg:end) = nan
    pcf%cpool_to_livecrootc_storage(beg:end) = nan
    pcf%cpool_to_deadcrootc(beg:end) = nan
    pcf%cpool_to_deadcrootc_storage(beg:end) = nan
    pcf%cpool_to_gresp_storage(beg:end) = nan
    pcf%cpool_leaf_gr(beg:end) = nan
    pcf%cpool_leaf_storage_gr(beg:end) = nan
    pcf%transfer_leaf_gr(beg:end) = nan
    pcf%cpool_froot_gr(beg:end) = nan
    pcf%cpool_froot_storage_gr(beg:end) = nan
    pcf%transfer_froot_gr(beg:end) = nan
    pcf%cpool_livestem_gr(beg:end) = nan
    pcf%cpool_livestem_storage_gr(beg:end) = nan
    pcf%transfer_livestem_gr(beg:end) = nan
    pcf%cpool_deadstem_gr(beg:end) = nan
    pcf%cpool_deadstem_storage_gr(beg:end) = nan
    pcf%transfer_deadstem_gr(beg:end) = nan
    pcf%cpool_livecroot_gr(beg:end) = nan
    pcf%cpool_livecroot_storage_gr(beg:end) = nan
    pcf%transfer_livecroot_gr(beg:end) = nan
    pcf%cpool_deadcroot_gr(beg:end) = nan
    pcf%cpool_deadcroot_storage_gr(beg:end) = nan
    pcf%transfer_deadcroot_gr(beg:end) = nan
    pcf%leafc_storage_to_xfer(beg:end) = nan
    pcf%frootc_storage_to_xfer(beg:end) = nan
    pcf%livestemc_storage_to_xfer(beg:end) = nan
    pcf%deadstemc_storage_to_xfer(beg:end) = nan
    pcf%livecrootc_storage_to_xfer(beg:end) = nan
    pcf%deadcrootc_storage_to_xfer(beg:end) = nan
    pcf%gresp_storage_to_xfer(beg:end) = nan
    pcf%livestemc_to_deadstemc(beg:end) = nan
    pcf%livecrootc_to_deadcrootc(beg:end) = nan
    pcf%gpp(beg:end) = nan
    pcf%mr(beg:end) = nan
    pcf%current_gr(beg:end) = nan
    pcf%transfer_gr(beg:end) = nan
    pcf%storage_gr(beg:end) = nan
    pcf%gr(beg:end) = nan
    pcf%ar(beg:end) = nan
    pcf%rr(beg:end) = nan
    pcf%npp(beg:end) = nan
    pcf%agnpp(beg:end) = nan
    pcf%bgnpp(beg:end) = nan
    pcf%litfall(beg:end) = nan
    pcf%vegfire(beg:end) = nan
    pcf%wood_harvestc(beg:end) = nan
    pcf%pft_cinputs(beg:end) = nan
    pcf%pft_coutputs(beg:end) = nan
    pcf%pft_fire_closs(beg:end) = nan
    if ( crop_prog )then
       pcf%xsmrpool_to_atm(beg:end)         = nan
       pcf%grainc_xfer_to_grainc(beg:end)   = nan
       pcf%livestemc_to_litter(beg:end)     = nan
       pcf%grainc_to_food(beg:end)          = nan
       pcf%cpool_to_grainc(beg:end)         = nan
       pcf%cpool_to_grainc_storage(beg:end) = nan
       pcf%cpool_grain_gr(beg:end)          = nan
       pcf%cpool_grain_storage_gr(beg:end)  = nan
       pcf%transfer_grain_gr(beg:end)       = nan
       pcf%grainc_storage_to_xfer(beg:end)  = nan
    end if
    if (use_cn) then
       pcf%frootc_alloc(beg:end) = nan
       pcf%frootc_loss(beg:end) = nan
       pcf%leafc_alloc(beg:end) = nan
       pcf%leafc_loss(beg:end) = nan
       pcf%woodc_alloc(beg:end) = nan
       pcf%woodc_loss(beg:end) = nan
    end if

  end subroutine init_pft_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nflux_type
!
! !INTERFACE:
  subroutine init_pft_nflux_type(beg, end, pnf)
!
! !DESCRIPTION:
! Initialize pft nitrogen flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nflux_type), intent(inout) :: pnf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pnf%m_leafn_to_litter(beg:end))
    allocate(pnf%m_frootn_to_litter(beg:end))
    allocate(pnf%m_leafn_storage_to_litter(beg:end))
    allocate(pnf%m_frootn_storage_to_litter(beg:end))
    allocate(pnf%m_livestemn_storage_to_litter(beg:end))
    allocate(pnf%m_deadstemn_storage_to_litter(beg:end))
    allocate(pnf%m_livecrootn_storage_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_litter(beg:end))
    allocate(pnf%m_leafn_xfer_to_litter(beg:end))
    allocate(pnf%m_frootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_litter(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_to_litter(beg:end))
    allocate(pnf%m_deadstemn_to_litter(beg:end))
    allocate(pnf%m_livecrootn_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_to_litter(beg:end))
    allocate(pnf%m_retransn_to_litter(beg:end))
    allocate(pnf%hrv_leafn_to_litter(beg:end))             
    allocate(pnf%hrv_frootn_to_litter(beg:end))            
    allocate(pnf%hrv_leafn_storage_to_litter(beg:end))     
    allocate(pnf%hrv_frootn_storage_to_litter(beg:end))    
    allocate(pnf%hrv_livestemn_storage_to_litter(beg:end)) 
    allocate(pnf%hrv_deadstemn_storage_to_litter(beg:end)) 
    allocate(pnf%hrv_livecrootn_storage_to_litter(beg:end))
    allocate(pnf%hrv_deadcrootn_storage_to_litter(beg:end))
    allocate(pnf%hrv_leafn_xfer_to_litter(beg:end))        
    allocate(pnf%hrv_frootn_xfer_to_litter(beg:end))       
    allocate(pnf%hrv_livestemn_xfer_to_litter(beg:end))    
    allocate(pnf%hrv_deadstemn_xfer_to_litter(beg:end))    
    allocate(pnf%hrv_livecrootn_xfer_to_litter(beg:end))   
    allocate(pnf%hrv_deadcrootn_xfer_to_litter(beg:end))   
    allocate(pnf%hrv_livestemn_to_litter(beg:end))         
    allocate(pnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(pnf%hrv_deadstemn_to_prod100n(beg:end))       
    allocate(pnf%hrv_livecrootn_to_litter(beg:end))        
    allocate(pnf%hrv_deadcrootn_to_litter(beg:end))        
    allocate(pnf%hrv_retransn_to_litter(beg:end))              
    allocate(pnf%m_leafn_to_fire(beg:end))
    allocate(pnf%m_frootn_to_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_litter_fire(beg:end))
    allocate(pnf%m_retransn_to_fire(beg:end))
    allocate(pnf%leafn_xfer_to_leafn(beg:end))
    allocate(pnf%frootn_xfer_to_frootn(beg:end))
    allocate(pnf%livestemn_xfer_to_livestemn(beg:end))
    allocate(pnf%deadstemn_xfer_to_deadstemn(beg:end))
    allocate(pnf%livecrootn_xfer_to_livecrootn(beg:end))
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(beg:end))
    allocate(pnf%leafn_to_litter(beg:end))
    allocate(pnf%leafn_to_retransn(beg:end))
    allocate(pnf%frootn_to_litter(beg:end))
    allocate(pnf%retransn_to_npool(beg:end))
    allocate(pnf%sminn_to_npool(beg:end))
    allocate(pnf%npool_to_leafn(beg:end))
    allocate(pnf%npool_to_leafn_storage(beg:end))
    allocate(pnf%npool_to_frootn(beg:end))
    allocate(pnf%npool_to_frootn_storage(beg:end))
    allocate(pnf%npool_to_livestemn(beg:end))
    allocate(pnf%npool_to_livestemn_storage(beg:end))
    allocate(pnf%npool_to_deadstemn(beg:end))
    allocate(pnf%npool_to_deadstemn_storage(beg:end))
    allocate(pnf%npool_to_livecrootn(beg:end))
    allocate(pnf%npool_to_livecrootn_storage(beg:end))
    allocate(pnf%npool_to_deadcrootn(beg:end))
    allocate(pnf%npool_to_deadcrootn_storage(beg:end))
    allocate(pnf%leafn_storage_to_xfer(beg:end))
    allocate(pnf%frootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_storage_to_xfer(beg:end))
    allocate(pnf%deadstemn_storage_to_xfer(beg:end))
    allocate(pnf%livecrootn_storage_to_xfer(beg:end))
    allocate(pnf%deadcrootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_to_deadstemn(beg:end))
    allocate(pnf%livestemn_to_retransn(beg:end))
    allocate(pnf%livecrootn_to_deadcrootn(beg:end))
    allocate(pnf%livecrootn_to_retransn(beg:end))
    allocate(pnf%ndeploy(beg:end))
    allocate(pnf%pft_ninputs(beg:end))
    allocate(pnf%pft_noutputs(beg:end))
    allocate(pnf%wood_harvestn(beg:end))
    allocate(pnf%pft_fire_nloss(beg:end))
    if ( crop_prog )then
       allocate(pnf%grainn_xfer_to_grainn(beg:end))
       allocate(pnf%livestemn_to_litter(beg:end))
       allocate(pnf%grainn_to_food(beg:end))
       allocate(pnf%npool_to_grainn(beg:end))
       allocate(pnf%npool_to_grainn_storage(beg:end))
       allocate(pnf%grainn_storage_to_xfer(beg:end))
    end if

    pnf%m_leafn_to_litter(beg:end) = nan
    pnf%m_frootn_to_litter(beg:end) = nan
    pnf%m_leafn_storage_to_litter(beg:end) = nan
    pnf%m_frootn_storage_to_litter(beg:end) = nan
    pnf%m_livestemn_storage_to_litter(beg:end) = nan
    pnf%m_deadstemn_storage_to_litter(beg:end) = nan
    pnf%m_livecrootn_storage_to_litter(beg:end) = nan
    pnf%m_deadcrootn_storage_to_litter(beg:end) = nan
    pnf%m_leafn_xfer_to_litter(beg:end) = nan
    pnf%m_frootn_xfer_to_litter(beg:end) = nan
    pnf%m_livestemn_xfer_to_litter(beg:end) = nan
    pnf%m_deadstemn_xfer_to_litter(beg:end) = nan
    pnf%m_livecrootn_xfer_to_litter(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_litter(beg:end) = nan
    pnf%m_livestemn_to_litter(beg:end) = nan
    pnf%m_deadstemn_to_litter(beg:end) = nan
    pnf%m_livecrootn_to_litter(beg:end) = nan
    pnf%m_deadcrootn_to_litter(beg:end) = nan
    pnf%m_retransn_to_litter(beg:end) = nan
    pnf%hrv_leafn_to_litter(beg:end) = nan             
    pnf%hrv_frootn_to_litter(beg:end) = nan            
    pnf%hrv_leafn_storage_to_litter(beg:end) = nan     
    pnf%hrv_frootn_storage_to_litter(beg:end) = nan    
    pnf%hrv_livestemn_storage_to_litter(beg:end) = nan 
    pnf%hrv_deadstemn_storage_to_litter(beg:end) = nan 
    pnf%hrv_livecrootn_storage_to_litter(beg:end) = nan
    pnf%hrv_deadcrootn_storage_to_litter(beg:end) = nan
    pnf%hrv_leafn_xfer_to_litter(beg:end) = nan        
    pnf%hrv_frootn_xfer_to_litter(beg:end) = nan       
    pnf%hrv_livestemn_xfer_to_litter(beg:end) = nan    
    pnf%hrv_deadstemn_xfer_to_litter(beg:end) = nan    
    pnf%hrv_livecrootn_xfer_to_litter(beg:end) = nan   
    pnf%hrv_deadcrootn_xfer_to_litter(beg:end) = nan   
    pnf%hrv_livestemn_to_litter(beg:end) = nan         
    pnf%hrv_deadstemn_to_prod10n(beg:end) = nan        
    pnf%hrv_deadstemn_to_prod100n(beg:end) = nan       
    pnf%hrv_livecrootn_to_litter(beg:end) = nan        
    pnf%hrv_deadcrootn_to_litter(beg:end) = nan        
    pnf%hrv_retransn_to_litter(beg:end) = nan           
    pnf%m_leafn_to_fire(beg:end) = nan
    pnf%m_frootn_to_fire(beg:end) = nan
    pnf%m_leafn_storage_to_fire(beg:end) = nan
    pnf%m_frootn_storage_to_fire(beg:end) = nan
    pnf%m_livestemn_storage_to_fire(beg:end) = nan
    pnf%m_deadstemn_storage_to_fire(beg:end) = nan
    pnf%m_livecrootn_storage_to_fire(beg:end) = nan
    pnf%m_deadcrootn_storage_to_fire(beg:end) = nan
    pnf%m_leafn_xfer_to_fire(beg:end) = nan
    pnf%m_frootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_xfer_to_fire(beg:end) = nan
    pnf%m_deadstemn_xfer_to_fire(beg:end) = nan
    pnf%m_livecrootn_xfer_to_fire(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_litter_fire(beg:end) = nan
    pnf%m_retransn_to_fire(beg:end) = nan
    pnf%leafn_xfer_to_leafn(beg:end) = nan
    pnf%frootn_xfer_to_frootn(beg:end) = nan
    pnf%livestemn_xfer_to_livestemn(beg:end) = nan
    pnf%deadstemn_xfer_to_deadstemn(beg:end) = nan
    pnf%livecrootn_xfer_to_livecrootn(beg:end) = nan
    pnf%deadcrootn_xfer_to_deadcrootn(beg:end) = nan
    pnf%leafn_to_litter(beg:end) = nan
    pnf%leafn_to_retransn(beg:end) = nan
    pnf%frootn_to_litter(beg:end) = nan
    pnf%retransn_to_npool(beg:end) = nan
    pnf%sminn_to_npool(beg:end) = nan
    pnf%npool_to_leafn(beg:end) = nan
    pnf%npool_to_leafn_storage(beg:end) = nan
    pnf%npool_to_frootn(beg:end) = nan
    pnf%npool_to_frootn_storage(beg:end) = nan
    pnf%npool_to_livestemn(beg:end) = nan
    pnf%npool_to_livestemn_storage(beg:end) = nan
    pnf%npool_to_deadstemn(beg:end) = nan
    pnf%npool_to_deadstemn_storage(beg:end) = nan
    pnf%npool_to_livecrootn(beg:end) = nan
    pnf%npool_to_livecrootn_storage(beg:end) = nan
    pnf%npool_to_deadcrootn(beg:end) = nan
    pnf%npool_to_deadcrootn_storage(beg:end) = nan
    pnf%leafn_storage_to_xfer(beg:end) = nan
    pnf%frootn_storage_to_xfer(beg:end) = nan
    pnf%livestemn_storage_to_xfer(beg:end) = nan
    pnf%deadstemn_storage_to_xfer(beg:end) = nan
    pnf%livecrootn_storage_to_xfer(beg:end) = nan
    pnf%deadcrootn_storage_to_xfer(beg:end) = nan
    pnf%livestemn_to_deadstemn(beg:end) = nan
    pnf%livestemn_to_retransn(beg:end) = nan
    pnf%livecrootn_to_deadcrootn(beg:end) = nan
    pnf%livecrootn_to_retransn(beg:end) = nan
    pnf%ndeploy(beg:end) = nan
    pnf%pft_ninputs(beg:end) = nan
    pnf%pft_noutputs(beg:end) = nan
    pnf%wood_harvestn(beg:end) = nan
    pnf%pft_fire_nloss(beg:end) = nan
    if ( crop_prog )then
       pnf%grainn_xfer_to_grainn(beg:end)   = nan
       pnf%livestemn_to_litter(beg:end)     = nan
       pnf%grainn_to_food(beg:end)          = nan
       pnf%npool_to_grainn(beg:end)         = nan
       pnf%npool_to_grainn_storage(beg:end) = nan
       pnf%grainn_storage_to_xfer(beg:end)  = nan
    end if

  end subroutine init_pft_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vflux_type
!
! !INTERFACE:
  subroutine init_pft_vflux_type(beg, end, pvf)
!
! !DESCRIPTION:
! Initialize pft VOC flux variables
!
    use clm_varcon, only : spval
    use shr_megan_mod, only: shr_megan_megcomps_n, shr_megan_mechcomps_n
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vflux_type), intent(inout) :: pvf

    integer :: i
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! (heald, 08/06)
!
!EOP
!------------------------------------------------------------------------

    if (shr_megan_mechcomps_n<1) return

    allocate(pvf%vocflx_tot(beg:end))
    allocate(pvf%vocflx(beg:end,1:shr_megan_mechcomps_n))
    allocate(pvf%Eopt_out(beg:end))
    allocate(pvf%topt_out(beg:end))
    allocate(pvf%alpha_out(beg:end))
    allocate(pvf%cp_out(beg:end))
    allocate(pvf%para_out(beg:end))
    allocate(pvf%par24a_out(beg:end))
    allocate(pvf%par240a_out(beg:end))
    allocate(pvf%paru_out(beg:end))
    allocate(pvf%par24u_out(beg:end))
    allocate(pvf%par240u_out(beg:end))
    allocate(pvf%gamma_out(beg:end))
    allocate(pvf%gammaL_out(beg:end))
    allocate(pvf%gammaT_out(beg:end))
    allocate(pvf%gammaP_out(beg:end))
    allocate(pvf%gammaA_out(beg:end))
    allocate(pvf%gammaS_out(beg:end))
    allocate(pvf%gammaC_out(beg:end))

    pvf%vocflx_tot(beg:end) = nan
    pvf%vocflx(beg:end,1:shr_megan_mechcomps_n) = nan
    pvf%Eopt_out(beg:end) = nan
    pvf%topt_out(beg:end) = nan
    pvf%alpha_out(beg:end) = nan
    pvf%cp_out(beg:end) = nan
    pvf%para_out(beg:end) = nan
    pvf%par24a_out(beg:end) = nan
    pvf%par240a_out(beg:end) = nan
    pvf%paru_out(beg:end) = nan
    pvf%par24u_out(beg:end) = nan
    pvf%par240u_out(beg:end) = nan
    pvf%gamma_out(beg:end) = nan
    pvf%gammaL_out(beg:end) = nan
    pvf%gammaT_out(beg:end) = nan
    pvf%gammaP_out(beg:end) = nan
    pvf%gammaA_out(beg:end) = nan
    pvf%gammaS_out(beg:end) = nan
    pvf%gammaC_out(beg:end) = nan

    allocate(pvf%meg(shr_megan_megcomps_n))

    do i=1,shr_megan_megcomps_n
       allocate(pvf%meg(i)%flux_out(beg:end))
       pvf%meg(i)%flux_out(beg:end) = nan
    enddo

  end subroutine init_pft_vflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_dflux_type
!
! !INTERFACE:
  subroutine init_pft_dflux_type(beg, end, pdf)
!
! !DESCRIPTION:
! Initialize pft dust flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dflux_type), intent(inout):: pdf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdf%flx_mss_vrt_dst(beg:end,1:ndst))
    allocate(pdf%flx_mss_vrt_dst_tot(beg:end))
    allocate(pdf%vlc_trb(beg:end,1:ndst))
    allocate(pdf%vlc_trb_1(beg:end))
    allocate(pdf%vlc_trb_2(beg:end))
    allocate(pdf%vlc_trb_3(beg:end))
    allocate(pdf%vlc_trb_4(beg:end))

    pdf%flx_mss_vrt_dst(beg:end,1:ndst) = nan
    pdf%flx_mss_vrt_dst_tot(beg:end) = nan
    pdf%vlc_trb(beg:end,1:ndst) = nan
    pdf%vlc_trb_1(beg:end) = nan
    pdf%vlc_trb_2(beg:end) = nan
    pdf%vlc_trb_3(beg:end) = nan
    pdf%vlc_trb_4(beg:end) = nan

  end subroutine init_pft_dflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_depvd_type
!
! !INTERFACE:
  subroutine init_pft_depvd_type(beg, end, pdd)

    use seq_drydep_mod, only:  n_drydep, drydep_method, DD_XLND
!
! !DESCRIPTION:
! Initialize pft dep velocity variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_depvd_type), intent(inout):: pdd
    integer :: i
!
! !REVISION HISTORY:
! Created by James Sulzman 541-929-6183
!
!EOP
!------------------------------------------------------------------------

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(pdd%drydepvel(beg:end,n_drydep))
       pdd%drydepvel = nan
    end if

  end subroutine init_pft_depvd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_pstate_type
!
! !INTERFACE:
  subroutine init_column_pstate_type(beg, end, cps)
!
! !DESCRIPTION:
! Initialize column physical state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_pstate_type), intent(inout):: cps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cps%snl(beg:end))      !* cannot be averaged up
    allocate(cps%isoicol(beg:end))  !* cannot be averaged up
    allocate(cps%bsw(beg:end,nlevgrnd))
    allocate(cps%watsat(beg:end,nlevgrnd))
    allocate(cps%watfc(beg:end,nlevgrnd))
    allocate(cps%watdry(beg:end,nlevgrnd)) 
    allocate(cps%watopt(beg:end,nlevgrnd)) 
    allocate(cps%hksat(beg:end,nlevgrnd))
    allocate(cps%sucsat(beg:end,nlevgrnd))
    allocate(cps%csol(beg:end,nlevgrnd))
    allocate(cps%tkmg(beg:end,nlevgrnd))
    allocate(cps%tkdry(beg:end,nlevgrnd))
    allocate(cps%tksatu(beg:end,nlevgrnd))
    allocate(cps%smpmin(beg:end))
    allocate(cps%hkdepth(beg:end))
    allocate(cps%wtfact(beg:end))
    allocate(cps%fracice(beg:end,nlevgrnd))
    allocate(cps%gwc_thr(beg:end))
    allocate(cps%mss_frc_cly_vld(beg:end))
    allocate(cps%mbl_bsn_fct(beg:end))
    allocate(cps%do_capsnow(beg:end))
    allocate(cps%snowdp(beg:end))
    allocate(cps%frac_sno (beg:end))
    allocate(cps%zi(beg:end,-nlevsno+0:nlevgrnd))
    allocate(cps%dz(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%z (beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%frac_iceold(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%imelt(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%eff_porosity(beg:end,nlevgrnd))
    allocate(cps%emg(beg:end))
    allocate(cps%z0mg(beg:end))
    allocate(cps%z0hg(beg:end))
    allocate(cps%z0qg(beg:end))
    allocate(cps%htvp(beg:end))
    allocate(cps%beta(beg:end))
    allocate(cps%zii(beg:end))
    allocate(cps%albgrd(beg:end,numrad))
    allocate(cps%albgri(beg:end,numrad))
    allocate(cps%rootr_column(beg:end,nlevgrnd))
    allocate(cps%rootfr_road_perv(beg:end,nlevgrnd))
    allocate(cps%rootr_road_perv(beg:end,nlevgrnd))
    allocate(cps%wf(beg:end))
!   allocate(cps%xirrig(beg:end))
    allocate(cps%max_dayl(beg:end))
    allocate(cps%bsw2(beg:end,nlevgrnd))
    allocate(cps%psisat(beg:end,nlevgrnd))
    allocate(cps%vwcsat(beg:end,nlevgrnd))
    allocate(cps%soilpsi(beg:end,nlevgrnd))
    allocate(cps%decl(beg:end))
    allocate(cps%coszen(beg:end))
    allocate(cps%fpi(beg:end))
    allocate(cps%fpg(beg:end))
    allocate(cps%annsum_counter(beg:end))
    allocate(cps%cannsum_npp(beg:end))
    allocate(cps%cannavg_t2m(beg:end))
    allocate(cps%me(beg:end))
    allocate(cps%fire_prob(beg:end))
    allocate(cps%mean_fire_prob(beg:end))
    allocate(cps%fireseasonl(beg:end))
    allocate(cps%farea_burned(beg:end))
    allocate(cps%ann_farea_burned(beg:end))
    allocate(cps%albsnd_hst(beg:end,numrad))
    allocate(cps%albsni_hst(beg:end,numrad))
    allocate(cps%albsod(beg:end,numrad))
    allocate(cps%albsoi(beg:end,numrad))
    allocate(cps%flx_absdv(beg:end,-nlevsno+1:1))
    allocate(cps%flx_absdn(beg:end,-nlevsno+1:1))
    allocate(cps%flx_absiv(beg:end,-nlevsno+1:1))
    allocate(cps%flx_absin(beg:end,-nlevsno+1:1))
    allocate(cps%snw_rds(beg:end,-nlevsno+1:0))
    allocate(cps%snw_rds_top(beg:end))
    allocate(cps%sno_liq_top(beg:end))
    allocate(cps%mss_bcpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_bcphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_bctot(beg:end,-nlevsno+1:0))
    allocate(cps%mss_bc_col(beg:end))
    allocate(cps%mss_bc_top(beg:end))
    allocate(cps%mss_ocpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_ocphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_octot(beg:end,-nlevsno+1:0))
    allocate(cps%mss_oc_col(beg:end))
    allocate(cps%mss_oc_top(beg:end))
    allocate(cps%mss_dst1(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst2(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst3(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst4(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dsttot(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst_col(beg:end))
    allocate(cps%mss_dst_top(beg:end))
    allocate(cps%h2osno_top(beg:end))
    allocate(cps%mss_cnc_bcphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_bcpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_ocphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_ocpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst1(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst2(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst3(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst4(beg:end,-nlevsno+1:0))
    allocate(cps%albgrd_pur(beg:end,numrad))
    allocate(cps%albgri_pur(beg:end,numrad))
    allocate(cps%albgrd_bc(beg:end,numrad))
    allocate(cps%albgri_bc(beg:end,numrad))
    allocate(cps%albgrd_oc(beg:end,numrad))
    allocate(cps%albgri_oc(beg:end,numrad))
    allocate(cps%albgrd_dst(beg:end,numrad))
    allocate(cps%albgri_dst(beg:end,numrad))
    allocate(cps%dTdz_top(beg:end))
    allocate(cps%snot_top(beg:end))
    allocate(cps%irrig_rate(beg:end))
    allocate(cps%n_irrig_steps_left(beg:end))
    allocate(cps%forc_pbot(beg:end))
    allocate(cps%forc_rho(beg:end))
    allocate(cps%glc_topo(beg:end))
    cps%isoicol(beg:end) = huge(1)
    cps%bsw(beg:end,1:nlevgrnd) = nan
    cps%watsat(beg:end,1:nlevgrnd) = nan
    cps%watfc(beg:end,1:nlevgrnd) = nan
    cps%watdry(beg:end,1:nlevgrnd) = nan
    cps%watopt(beg:end,1:nlevgrnd) = nan
    cps%hksat(beg:end,1:nlevgrnd) = nan
    cps%sucsat(beg:end,1:nlevgrnd) = nan
    cps%csol(beg:end,1:nlevgrnd) = nan
    cps%tkmg(beg:end,1:nlevgrnd) = nan
    cps%tkdry(beg:end,1:nlevgrnd) = nan
    cps%tksatu(beg:end,1:nlevgrnd) = nan
    cps%smpmin(beg:end) = nan
    cps%hkdepth(beg:end) = nan
    cps%wtfact(beg:end) = nan
    cps%fracice(beg:end,1:nlevgrnd) = nan
    cps%gwc_thr(beg:end) = nan
    cps%mss_frc_cly_vld(beg:end) = nan
    cps%mbl_bsn_fct(beg:end) = nan
    cps%do_capsnow (beg:end)= .false.
    cps%snowdp(beg:end) = nan
    cps%frac_sno(beg:end) = nan
    cps%zi(beg:end,-nlevsno+0:nlevgrnd) = nan
    cps%dz(beg:end,-nlevsno+1:nlevgrnd) = nan
    cps%z (beg:end,-nlevsno+1:nlevgrnd) = nan
    cps%frac_iceold(beg:end,-nlevsno+1:nlevgrnd) = spval
    cps%imelt(beg:end,-nlevsno+1:nlevgrnd) = huge(1)
    cps%eff_porosity(beg:end,1:nlevgrnd) = spval
    cps%emg(beg:end) = nan
    cps%z0mg(beg:end) = nan
    cps%z0hg(beg:end) = nan
    cps%z0qg(beg:end) = nan
    cps%htvp(beg:end) = nan
    cps%beta(beg:end) = nan
    cps%zii(beg:end) = nan
    cps%albgrd(beg:end,:numrad) = nan
    cps%albgri(beg:end,:numrad) = nan
    cps%rootr_column(beg:end,1:nlevgrnd) = spval
    cps%rootfr_road_perv(beg:end,1:nlevurb) = nan
    cps%rootr_road_perv(beg:end,1:nlevurb) = nan
    cps%wf(beg:end) = nan
!   cps%xirrig(beg:end) = 0._r8
    cps%bsw2(beg:end,1:nlevgrnd) = nan
    cps%psisat(beg:end,1:nlevgrnd) = nan
    cps%vwcsat(beg:end,1:nlevgrnd) = nan
    cps%soilpsi(beg:end,1:nlevgrnd) = spval
    cps%decl(beg:end) = nan
    cps%coszen(beg:end) = nan
    cps%fpi(beg:end) = nan
    cps%fpg(beg:end) = nan
    cps%annsum_counter(beg:end) = nan
    cps%cannsum_npp(beg:end) = nan
    cps%cannavg_t2m(beg:end) = nan
    cps%me(beg:end) = nan
    cps%fire_prob(beg:end) = nan
    cps%mean_fire_prob(beg:end) = nan
    cps%fireseasonl(beg:end) = nan
    cps%farea_burned(beg:end) = nan
    cps%ann_farea_burned(beg:end) = nan
    cps%albsnd_hst(beg:end,:numrad) = spval
    cps%albsni_hst(beg:end,:numrad) = spval
    cps%albsod(beg:end,:numrad) = nan
    cps%albsoi(beg:end,:numrad) = nan
    cps%flx_absdv(beg:end,-nlevsno+1:1) = spval
    cps%flx_absdn(beg:end,-nlevsno+1:1) = spval
    cps%flx_absiv(beg:end,-nlevsno+1:1) = spval
    cps%flx_absin(beg:end,-nlevsno+1:1) = spval
    cps%snw_rds(beg:end,-nlevsno+1:0) = nan
    cps%snw_rds_top(beg:end) = nan
    cps%sno_liq_top(beg:end) = nan
    cps%mss_bcpho(beg:end,-nlevsno+1:0) = nan
    cps%mss_bcphi(beg:end,-nlevsno+1:0) = nan
    cps%mss_bctot(beg:end,-nlevsno+1:0) = nan
    cps%mss_bc_col(beg:end) = nan
    cps%mss_bc_top(beg:end) = nan
    cps%mss_ocpho(beg:end,-nlevsno+1:0) = nan
    cps%mss_ocphi(beg:end,-nlevsno+1:0) = nan
    cps%mss_octot(beg:end,-nlevsno+1:0) = nan
    cps%mss_oc_col(beg:end) = nan
    cps%mss_oc_top(beg:end) = nan
    cps%mss_dst1(beg:end,-nlevsno+1:0) = nan
    cps%mss_dst2(beg:end,-nlevsno+1:0) = nan
    cps%mss_dst3(beg:end,-nlevsno+1:0) = nan
    cps%mss_dst4(beg:end,-nlevsno+1:0) = nan
    cps%mss_dsttot(beg:end,-nlevsno+1:0) = nan
    cps%mss_dst_col(beg:end) = nan
    cps%mss_dst_top(beg:end) = nan
    cps%h2osno_top(beg:end) = nan
    cps%mss_cnc_bcphi(beg:end,-nlevsno+1:0) = nan
    cps%mss_cnc_bcpho(beg:end,-nlevsno+1:0) = nan
    cps%mss_cnc_ocphi(beg:end,-nlevsno+1:0) = nan
    cps%mss_cnc_ocpho(beg:end,-nlevsno+1:0) = nan
    cps%mss_cnc_dst1(beg:end,-nlevsno+1:0) = nan
    cps%mss_cnc_dst2(beg:end,-nlevsno+1:0) = nan
    cps%mss_cnc_dst3(beg:end,-nlevsno+1:0) = nan
    cps%mss_cnc_dst4(beg:end,-nlevsno+1:0) = nan
    cps%albgrd_pur(beg:end,:numrad) = nan
    cps%albgri_pur(beg:end,:numrad) = nan
    cps%albgrd_bc(beg:end,:numrad) = nan
    cps%albgri_bc(beg:end,:numrad) = nan
    cps%albgrd_oc(beg:end,:numrad) = nan
    cps%albgri_oc(beg:end,:numrad) = nan 
    cps%albgrd_dst(beg:end,:numrad) = nan
    cps%albgri_dst(beg:end,:numrad) = nan
    cps%dTdz_top(beg:end) = nan
    cps%snot_top(beg:end) = nan
    cps%irrig_rate(beg:end) = nan
    cps%n_irrig_steps_left(beg:end) = 0
    cps%forc_pbot(beg:end) = nan
    cps%forc_rho(beg:end) = nan
    cps%glc_topo(beg:end) = nan

  end subroutine init_column_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_estate_type
!
! !INTERFACE:
  subroutine init_column_estate_type(beg, end, ces)
!
! !DESCRIPTION:
! Initialize column energy state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_estate_type), intent(inout):: ces
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(ces%t_grnd(beg:end))
    allocate(ces%t_grnd_u(beg:end))
    allocate(ces%t_grnd_r(beg:end))
    allocate(ces%dt_grnd(beg:end))
    allocate(ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd))
    allocate(ces%t_soi_10cm(beg:end))
    allocate(ces%t_lake(beg:end,1:nlevlak))
    allocate(ces%tssbef(beg:end,-nlevsno+1:nlevgrnd))
    allocate(ces%thv(beg:end))
    allocate(ces%hc_soi(beg:end))
    allocate(ces%hc_soisno(beg:end))
    allocate(ces%forc_t(beg:end))
    allocate(ces%forc_th(beg:end))

    ces%t_grnd(beg:end)    = nan
    ces%t_grnd_u(beg:end)  = nan
    ces%t_grnd_r(beg:end)  = nan
    ces%dt_grnd(beg:end)   = nan
    ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd) = spval
    ces%t_soi_10cm(beg:end) = spval
    ces%t_lake(beg:end,1:nlevlak)            = nan
    ces%tssbef(beg:end,-nlevsno+1:nlevgrnd)   = nan
    ces%thv(beg:end)       = nan
    ces%hc_soi(beg:end)    = nan
    ces%hc_soisno(beg:end) = nan
    ces%forc_t(beg:end) = nan
    ces%forc_th(beg:end) = nan

  end subroutine init_column_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wstate_type
!
! !INTERFACE:
  subroutine init_column_wstate_type(beg, end, cws)
!
! !DESCRIPTION:
! Initialize column water state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wstate_type), intent(inout):: cws !column water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cws%h2osno(beg:end))
    allocate(cws%errh2osno(beg:end))
    allocate(cws%snow_sources(beg:end))
    allocate(cws%snow_sinks(beg:end))
    allocate(cws%h2osoi_liq(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cws%h2osoi_ice(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cws%h2osoi_liqice_10cm(beg:end))
    allocate(cws%h2osoi_vol(beg:end,1:nlevgrnd))
    allocate(cws%h2osno_old(beg:end))
    allocate(cws%qg(beg:end))
    allocate(cws%dqgdT(beg:end))
    allocate(cws%snowice(beg:end))
    allocate(cws%snowliq(beg:end))
    allocate(cws%soilalpha(beg:end))
    allocate(cws%soilbeta(beg:end))
    allocate(cws%soilalpha_u(beg:end))
    allocate(cws%zwt(beg:end))
    allocate(cws%fcov(beg:end))
    allocate(cws%fsat(beg:end))
    allocate(cws%wa(beg:end))
    allocate(cws%wt(beg:end))
    allocate(cws%qcharge(beg:end))
    allocate(cws%smp_l(beg:end,1:nlevgrnd))
    allocate(cws%hk_l(beg:end,1:nlevgrnd))
    allocate(cws%forc_q(beg:end))

    cws%h2osno(beg:end) = nan
    cws%errh2osno(beg:end) = nan
    cws%snow_sources(beg:end) = nan
    cws%snow_sinks(beg:end) = nan
    cws%h2osoi_liq(beg:end,-nlevsno+1:nlevgrnd)= spval
    cws%h2osoi_ice(beg:end,-nlevsno+1:nlevgrnd) = spval
    cws%h2osoi_liqice_10cm(beg:end) = spval
    cws%h2osoi_vol(beg:end,1:nlevgrnd) = spval
    cws%h2osno_old(beg:end) = nan
    cws%qg(beg:end) = nan
    cws%dqgdT(beg:end) = nan
    cws%snowice(beg:end) = nan
    cws%snowliq(beg:end) = nan
    cws%soilalpha(beg:end) = nan
    cws%soilbeta(beg:end) = nan
    cws%soilalpha_u(beg:end) = nan
    cws%zwt(beg:end) = nan
    cws%fcov(beg:end) = nan
    cws%fsat(beg:end) = nan
    cws%wa(beg:end) = nan
    cws%wt(beg:end) = nan
    cws%qcharge(beg:end) = nan
    cws%smp_l(beg:end,1:nlevgrnd) = spval
    cws%hk_l(beg:end,1:nlevgrnd) = spval
    cws%forc_q(beg:end) = nan

  end subroutine init_column_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cstate_type
!
! !INTERFACE:
  subroutine init_column_cstate_type(beg, end, ccs)
!
! !DESCRIPTION:
! Initialize column carbon state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cstate_type), intent(inout):: ccs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ccs%soilc(beg:end))
    allocate(ccs%cwdc(beg:end))
    allocate(ccs%litr1c(beg:end))
    allocate(ccs%litr2c(beg:end))
    allocate(ccs%litr3c(beg:end))
    allocate(ccs%soil1c(beg:end))
    allocate(ccs%soil2c(beg:end))
    allocate(ccs%soil3c(beg:end))
    allocate(ccs%soil4c(beg:end))
    allocate(ccs%seedc(beg:end))
    allocate(ccs%col_ctrunc(beg:end))
    allocate(ccs%prod10c(beg:end))
    allocate(ccs%prod100c(beg:end))
    allocate(ccs%totprodc(beg:end))
    allocate(ccs%totlitc(beg:end))
    allocate(ccs%totsomc(beg:end))
    allocate(ccs%totecosysc(beg:end))
    allocate(ccs%totcolc(beg:end))

    ccs%soilc(beg:end) = nan
    ccs%cwdc(beg:end) = nan
    ccs%litr1c(beg:end) = nan
    ccs%litr2c(beg:end) = nan
    ccs%litr3c(beg:end) = nan
    ccs%soil1c(beg:end) = nan
    ccs%soil2c(beg:end) = nan
    ccs%soil3c(beg:end) = nan
    ccs%soil4c(beg:end) = nan
    ccs%seedc(beg:end) = nan
    ccs%col_ctrunc(beg:end) = nan
    ccs%prod10c(beg:end) = nan
    ccs%prod100c(beg:end) = nan
    ccs%totprodc(beg:end) = nan
    ccs%totlitc(beg:end) = nan
    ccs%totsomc(beg:end) = nan
    ccs%totecosysc(beg:end) = nan
    ccs%totcolc(beg:end) = nan

  end subroutine init_column_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nstate_type
!
! !INTERFACE:
  subroutine init_column_nstate_type(beg, end, cns)
!
! !DESCRIPTION:
! Initialize column nitrogen state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nstate_type), intent(inout):: cns
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cns%cwdn(beg:end))
    allocate(cns%litr1n(beg:end))
    allocate(cns%litr2n(beg:end))
    allocate(cns%litr3n(beg:end))
    allocate(cns%soil1n(beg:end))
    allocate(cns%soil2n(beg:end))
    allocate(cns%soil3n(beg:end))
    allocate(cns%soil4n(beg:end))
    allocate(cns%sminn(beg:end))
    allocate(cns%col_ntrunc(beg:end))
    allocate(cns%seedn(beg:end))
    allocate(cns%prod10n(beg:end))
    allocate(cns%prod100n(beg:end))
    allocate(cns%totprodn(beg:end))
    allocate(cns%totlitn(beg:end))
    allocate(cns%totsomn(beg:end))
    allocate(cns%totecosysn(beg:end))
    allocate(cns%totcoln(beg:end))

    cns%cwdn(beg:end) = nan
    cns%litr1n(beg:end) = nan
    cns%litr2n(beg:end) = nan
    cns%litr3n(beg:end) = nan
    cns%soil1n(beg:end) = nan
    cns%soil2n(beg:end) = nan
    cns%soil3n(beg:end) = nan
    cns%soil4n(beg:end) = nan
    cns%sminn(beg:end) = nan
    cns%col_ntrunc(beg:end) = nan
    cns%seedn(beg:end) = nan
    cns%prod10n(beg:end) = nan
    cns%prod100n(beg:end) = nan
    cns%totprodn(beg:end) = nan
    cns%totlitn(beg:end) = nan
    cns%totsomn(beg:end) = nan
    cns%totecosysn(beg:end) = nan
    cns%totcoln(beg:end) = nan

  end subroutine init_column_nstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_eflux_type
!
! !INTERFACE:
  subroutine init_column_eflux_type(beg, end, cef)
!
! !DESCRIPTION:
! Initialize column energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_eflux_type), intent(inout):: cef
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cef%eflx_snomelt(beg:end))
    allocate(cef%eflx_snomelt_u(beg:end))
    allocate(cef%eflx_snomelt_r(beg:end))
    allocate(cef%eflx_impsoil(beg:end))
    allocate(cef%eflx_fgr12(beg:end))
    allocate(cef%eflx_building_heat(beg:end))
    allocate(cef%eflx_urban_ac(beg:end))
    allocate(cef%eflx_urban_heat(beg:end))
    allocate(cef%eflx_bot(beg:end))

    cef%eflx_snomelt(beg:end)       = nan
    cef%eflx_snomelt_u(beg:end)       = nan
    cef%eflx_snomelt_r(beg:end)       = nan
    cef%eflx_impsoil(beg:end)       = nan
    cef%eflx_fgr12(beg:end)         = nan
    cef%eflx_building_heat(beg:end) = nan
    cef%eflx_urban_ac(beg:end) = nan
    cef%eflx_urban_heat(beg:end) = nan
    cef%eflx_bot(beg:end) = nan

  end subroutine init_column_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wflux_type
!
! !INTERFACE:
  subroutine init_column_wflux_type(beg, end, cwf)
!
! !DESCRIPTION:
! Initialize column water flux variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wflux_type), intent(inout):: cwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cwf%qflx_infl(beg:end))
    allocate(cwf%qflx_surf(beg:end))
    allocate(cwf%qflx_drain(beg:end))
    allocate(cwf%qflx_top_soil(beg:end))
    allocate(cwf%qflx_sl_top_soil(beg:end))
    allocate(cwf%qflx_snomelt(beg:end))
    allocate(cwf%qflx_qrgwl(beg:end))
    allocate(cwf%qflx_runoff(beg:end))
    allocate(cwf%qflx_runoff_u(beg:end))
    allocate(cwf%qflx_runoff_r(beg:end))
    allocate(cwf%qmelt(beg:end))
    allocate(cwf%h2ocan_loss(beg:end))
    allocate(cwf%qflx_rsub_sat(beg:end))
    allocate(cwf%flx_bc_dep_dry(beg:end))
    allocate(cwf%flx_bc_dep_wet(beg:end))
    allocate(cwf%flx_bc_dep_pho(beg:end))
    allocate(cwf%flx_bc_dep_phi(beg:end))
    allocate(cwf%flx_bc_dep(beg:end))
    allocate(cwf%flx_oc_dep_dry(beg:end))
    allocate(cwf%flx_oc_dep_wet(beg:end))
    allocate(cwf%flx_oc_dep_pho(beg:end))
    allocate(cwf%flx_oc_dep_phi(beg:end))
    allocate(cwf%flx_oc_dep(beg:end))
    allocate(cwf%flx_dst_dep_dry1(beg:end))
    allocate(cwf%flx_dst_dep_wet1(beg:end))
    allocate(cwf%flx_dst_dep_dry2(beg:end))
    allocate(cwf%flx_dst_dep_wet2(beg:end))
    allocate(cwf%flx_dst_dep_dry3(beg:end))
    allocate(cwf%flx_dst_dep_wet3(beg:end))
    allocate(cwf%flx_dst_dep_dry4(beg:end))
    allocate(cwf%flx_dst_dep_wet4(beg:end))
    allocate(cwf%flx_dst_dep(beg:end))
    allocate(cwf%qflx_snofrz_lyr(beg:end,-nlevsno+1:0))
    allocate(cwf%qflx_snofrz_col(beg:end))
    allocate(cwf%qflx_irrig(beg:end))
    allocate(cwf%qflx_glcice(beg:end))
    allocate(cwf%qflx_glcice_frz(beg:end))
    allocate(cwf%qflx_glcice_melt(beg:end))
    allocate(cwf%glc_rofi(beg:end))
    allocate(cwf%glc_rofl(beg:end))
    allocate(cwf%qflx_floodc(beg:end))
    allocate(cwf%qflx_snow_melt(beg:end))

    cwf%qflx_infl(beg:end) = nan
    cwf%qflx_surf(beg:end) = nan
    cwf%qflx_drain(beg:end) = nan
    cwf%qflx_top_soil(beg:end) = spval
    cwf%qflx_sl_top_soil(beg:end) = nan
    cwf%qflx_snomelt(beg:end) = nan
    cwf%qflx_qrgwl(beg:end) = nan
    cwf%qflx_runoff(beg:end) = nan
    cwf%qflx_runoff_u(beg:end) = nan
    cwf%qflx_runoff_r(beg:end) = nan
    cwf%qmelt(beg:end) = nan
    cwf%h2ocan_loss(beg:end) = nan
    cwf%qflx_rsub_sat(beg:end) = nan
    cwf%flx_bc_dep_dry(beg:end) = nan
    cwf%flx_bc_dep_wet(beg:end) = nan
    cwf%flx_bc_dep_pho(beg:end) = nan
    cwf%flx_bc_dep_phi(beg:end) = nan
    cwf%flx_bc_dep(beg:end) = nan
    cwf%flx_oc_dep_dry(beg:end) = nan
    cwf%flx_oc_dep_wet(beg:end) = nan
    cwf%flx_oc_dep_pho(beg:end) = nan
    cwf%flx_oc_dep_phi(beg:end) = nan
    cwf%flx_oc_dep(beg:end) = nan
    cwf%flx_dst_dep_dry1(beg:end) = nan
    cwf%flx_dst_dep_wet1(beg:end) = nan
    cwf%flx_dst_dep_dry2(beg:end) = nan
    cwf%flx_dst_dep_wet2(beg:end) = nan
    cwf%flx_dst_dep_dry3(beg:end) = nan
    cwf%flx_dst_dep_wet3(beg:end) = nan
    cwf%flx_dst_dep_dry4(beg:end) = nan
    cwf%flx_dst_dep_wet4(beg:end) = nan
    cwf%flx_dst_dep(beg:end) = nan
    cwf%qflx_snofrz_lyr(beg:end,-nlevsno+1:0) = spval
    cwf%qflx_snofrz_col(beg:end) = nan
    cwf%qflx_irrig(beg:end)  = nan
    cwf%qflx_glcice(beg:end) = nan
    cwf%qflx_glcice_frz(beg:end) = nan
    cwf%qflx_glcice_melt(beg:end) = nan
    cwf%glc_rofi(beg:end)    = nan
    cwf%glc_rofl(beg:end)    = nan
    cwf%qflx_floodc(beg:end) = spval
    cwf%qflx_snow_melt(beg:end) = spval

  end subroutine init_column_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cflux_type
!
! !INTERFACE:
  subroutine init_column_cflux_type(beg, end, ccf)
!
! !DESCRIPTION:
! Initialize column carbon flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cflux_type), intent(inout):: ccf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(ccf%m_leafc_to_litr1c(beg:end))
    allocate(ccf%m_leafc_to_litr2c(beg:end))
    allocate(ccf%m_leafc_to_litr3c(beg:end))
    allocate(ccf%m_frootc_to_litr1c(beg:end))
    allocate(ccf%m_frootc_to_litr2c(beg:end))
    allocate(ccf%m_frootc_to_litr3c(beg:end))
    allocate(ccf%m_leafc_storage_to_litr1c(beg:end))
    allocate(ccf%m_frootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_leafc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_frootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_to_cwdc(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc(beg:end))
    allocate(ccf%m_livecrootc_to_cwdc(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc(beg:end))
    allocate(ccf%m_gresp_storage_to_litr1c(beg:end))
    allocate(ccf%m_gresp_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc_fire(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc_fire(beg:end))
    allocate(ccf%hrv_leafc_to_litr1c(beg:end))             
    allocate(ccf%hrv_leafc_to_litr2c(beg:end))             
    allocate(ccf%hrv_leafc_to_litr3c(beg:end))             
    allocate(ccf%hrv_frootc_to_litr1c(beg:end))            
    allocate(ccf%hrv_frootc_to_litr2c(beg:end))            
    allocate(ccf%hrv_frootc_to_litr3c(beg:end))            
    allocate(ccf%hrv_livestemc_to_cwdc(beg:end))           
    allocate(ccf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(ccf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(ccf%hrv_livecrootc_to_cwdc(beg:end))          
    allocate(ccf%hrv_deadcrootc_to_cwdc(beg:end))          
    allocate(ccf%hrv_leafc_storage_to_litr1c(beg:end))     
    allocate(ccf%hrv_frootc_storage_to_litr1c(beg:end))    
    allocate(ccf%hrv_livestemc_storage_to_litr1c(beg:end)) 
    allocate(ccf%hrv_deadstemc_storage_to_litr1c(beg:end)) 
    allocate(ccf%hrv_livecrootc_storage_to_litr1c(beg:end))
    allocate(ccf%hrv_deadcrootc_storage_to_litr1c(beg:end))
    allocate(ccf%hrv_gresp_storage_to_litr1c(beg:end))     
    allocate(ccf%hrv_leafc_xfer_to_litr1c(beg:end))        
    allocate(ccf%hrv_frootc_xfer_to_litr1c(beg:end))       
    allocate(ccf%hrv_livestemc_xfer_to_litr1c(beg:end))    
    allocate(ccf%hrv_deadstemc_xfer_to_litr1c(beg:end))    
    allocate(ccf%hrv_livecrootc_xfer_to_litr1c(beg:end))   
    allocate(ccf%hrv_deadcrootc_xfer_to_litr1c(beg:end))   
    allocate(ccf%hrv_gresp_xfer_to_litr1c(beg:end))        
    allocate(ccf%m_litr1c_to_fire(beg:end))
    allocate(ccf%m_litr2c_to_fire(beg:end))
    allocate(ccf%m_litr3c_to_fire(beg:end))
    allocate(ccf%m_cwdc_to_fire(beg:end))
    if ( crop_prog )then
       allocate(ccf%grainc_to_litr1c(beg:end))
       allocate(ccf%grainc_to_litr2c(beg:end))
       allocate(ccf%grainc_to_litr3c(beg:end))
       allocate(ccf%livestemc_to_litr1c(beg:end))
       allocate(ccf%livestemc_to_litr2c(beg:end))
       allocate(ccf%livestemc_to_litr3c(beg:end))
    end if
    allocate(ccf%leafc_to_litr1c(beg:end))
    allocate(ccf%leafc_to_litr2c(beg:end))
    allocate(ccf%leafc_to_litr3c(beg:end))
    allocate(ccf%frootc_to_litr1c(beg:end))
    allocate(ccf%frootc_to_litr2c(beg:end))
    allocate(ccf%frootc_to_litr3c(beg:end))
    allocate(ccf%cwdc_to_litr2c(beg:end))
    allocate(ccf%cwdc_to_litr3c(beg:end))
    allocate(ccf%litr1_hr(beg:end))
    allocate(ccf%litr1c_to_soil1c(beg:end))
    allocate(ccf%litr2_hr(beg:end))
    allocate(ccf%litr2c_to_soil2c(beg:end))
    allocate(ccf%litr3_hr(beg:end))
    allocate(ccf%litr3c_to_soil3c(beg:end))
    allocate(ccf%soil1_hr(beg:end))
    allocate(ccf%soil1c_to_soil2c(beg:end))
    allocate(ccf%soil2_hr(beg:end))
    allocate(ccf%soil2c_to_soil3c(beg:end))
    allocate(ccf%soil3_hr(beg:end))
    allocate(ccf%soil3c_to_soil4c(beg:end))
    allocate(ccf%soil4_hr(beg:end))
    if (use_cn) then
       allocate(ccf%dwt_seedc_to_leaf(beg:end))
       allocate(ccf%dwt_seedc_to_deadstem(beg:end))
       allocate(ccf%dwt_conv_cflux(beg:end))
       allocate(ccf%dwt_prod10c_gain(beg:end))
       allocate(ccf%dwt_prod100c_gain(beg:end))
       allocate(ccf%dwt_frootc_to_litr1c(beg:end))
       allocate(ccf%dwt_frootc_to_litr2c(beg:end))
       allocate(ccf%dwt_frootc_to_litr3c(beg:end))
       allocate(ccf%dwt_livecrootc_to_cwdc(beg:end))
       allocate(ccf%dwt_deadcrootc_to_cwdc(beg:end))
       allocate(ccf%dwt_closs(beg:end))
       allocate(ccf%landuseflux(beg:end))
       allocate(ccf%landuptake(beg:end))
       allocate(ccf%prod10c_loss(beg:end))
       allocate(ccf%prod100c_loss(beg:end))
       allocate(ccf%product_closs(beg:end))
    end if
    allocate(ccf%lithr(beg:end))
    allocate(ccf%somhr(beg:end))
    allocate(ccf%hr(beg:end))
    allocate(ccf%sr(beg:end))
    allocate(ccf%er(beg:end))
    allocate(ccf%litfire(beg:end))
    allocate(ccf%somfire(beg:end))
    allocate(ccf%totfire(beg:end))
    allocate(ccf%nep(beg:end))
    allocate(ccf%nbp(beg:end))
    allocate(ccf%nee(beg:end))
    allocate(ccf%col_cinputs(beg:end))
    allocate(ccf%col_coutputs(beg:end))
    allocate(ccf%col_fire_closs(beg:end))

    if (use_cn) then
       allocate(ccf%cwdc_hr(beg:end))
       allocate(ccf%cwdc_loss(beg:end))
       allocate(ccf%litterc_loss(beg:end))
    end if

    ccf%m_leafc_to_litr1c(beg:end)                = nan
    ccf%m_leafc_to_litr2c(beg:end)                = nan
    ccf%m_leafc_to_litr3c(beg:end)                = nan
    ccf%m_frootc_to_litr1c(beg:end)               = nan
    ccf%m_frootc_to_litr2c(beg:end)               = nan
    ccf%m_frootc_to_litr3c(beg:end)               = nan
    ccf%m_leafc_storage_to_litr1c(beg:end)        = nan
    ccf%m_frootc_storage_to_litr1c(beg:end)       = nan
    ccf%m_livestemc_storage_to_litr1c(beg:end)    = nan
    ccf%m_deadstemc_storage_to_litr1c(beg:end)    = nan
    ccf%m_livecrootc_storage_to_litr1c(beg:end)   = nan
    ccf%m_deadcrootc_storage_to_litr1c(beg:end)   = nan
    ccf%m_leafc_xfer_to_litr1c(beg:end)           = nan
    ccf%m_frootc_xfer_to_litr1c(beg:end)          = nan
    ccf%m_livestemc_xfer_to_litr1c(beg:end)       = nan
    ccf%m_deadstemc_xfer_to_litr1c(beg:end)       = nan
    ccf%m_livecrootc_xfer_to_litr1c(beg:end)      = nan
    ccf%m_deadcrootc_xfer_to_litr1c(beg:end)      = nan
    ccf%m_livestemc_to_cwdc(beg:end)              = nan
    ccf%m_deadstemc_to_cwdc(beg:end)              = nan
    ccf%m_livecrootc_to_cwdc(beg:end)             = nan
    ccf%m_deadcrootc_to_cwdc(beg:end)             = nan
    ccf%m_gresp_storage_to_litr1c(beg:end)        = nan
    ccf%m_gresp_xfer_to_litr1c(beg:end)           = nan
    ccf%m_deadstemc_to_cwdc_fire(beg:end)         = nan
    ccf%m_deadcrootc_to_cwdc_fire(beg:end)        = nan
    ccf%hrv_leafc_to_litr1c(beg:end)              = nan             
    ccf%hrv_leafc_to_litr2c(beg:end)              = nan             
    ccf%hrv_leafc_to_litr3c(beg:end)              = nan             
    ccf%hrv_frootc_to_litr1c(beg:end)             = nan            
    ccf%hrv_frootc_to_litr2c(beg:end)             = nan            
    ccf%hrv_frootc_to_litr3c(beg:end)             = nan            
    ccf%hrv_livestemc_to_cwdc(beg:end)            = nan           
    ccf%hrv_deadstemc_to_prod10c(beg:end)         = nan        
    ccf%hrv_deadstemc_to_prod100c(beg:end)        = nan       
    ccf%hrv_livecrootc_to_cwdc(beg:end)           = nan          
    ccf%hrv_deadcrootc_to_cwdc(beg:end)           = nan          
    ccf%hrv_leafc_storage_to_litr1c(beg:end)      = nan     
    ccf%hrv_frootc_storage_to_litr1c(beg:end)     = nan    
    ccf%hrv_livestemc_storage_to_litr1c(beg:end)  = nan 
    ccf%hrv_deadstemc_storage_to_litr1c(beg:end)  = nan 
    ccf%hrv_livecrootc_storage_to_litr1c(beg:end) = nan
    ccf%hrv_deadcrootc_storage_to_litr1c(beg:end) = nan
    if ( crop_prog )then
       ccf%grainc_to_litr1c(beg:end)    = nan
       ccf%grainc_to_litr2c(beg:end)    = nan
       ccf%grainc_to_litr3c(beg:end)    = nan
       ccf%livestemc_to_litr1c(beg:end) = nan
       ccf%livestemc_to_litr2c(beg:end) = nan
       ccf%livestemc_to_litr3c(beg:end) = nan
    end if
    ccf%hrv_gresp_storage_to_litr1c(beg:end)      = nan     
    ccf%hrv_leafc_xfer_to_litr1c(beg:end)         = nan        
    ccf%hrv_frootc_xfer_to_litr1c(beg:end)        = nan       
    ccf%hrv_livestemc_xfer_to_litr1c(beg:end)     = nan    
    ccf%hrv_deadstemc_xfer_to_litr1c(beg:end)     = nan    
    ccf%hrv_livecrootc_xfer_to_litr1c(beg:end)    = nan   
    ccf%hrv_deadcrootc_xfer_to_litr1c(beg:end)    = nan   
    ccf%hrv_gresp_xfer_to_litr1c(beg:end)         = nan        
    ccf%m_litr1c_to_fire(beg:end)                 = nan
    ccf%m_litr2c_to_fire(beg:end)                 = nan
    ccf%m_litr3c_to_fire(beg:end)                 = nan
    ccf%m_cwdc_to_fire(beg:end)                   = nan
    ccf%leafc_to_litr1c(beg:end)                  = nan
    ccf%leafc_to_litr2c(beg:end)                  = nan
    ccf%leafc_to_litr3c(beg:end)                  = nan
    ccf%frootc_to_litr1c(beg:end)                 = nan
    ccf%frootc_to_litr2c(beg:end)                 = nan
    ccf%frootc_to_litr3c(beg:end)                 = nan
    ccf%cwdc_to_litr2c(beg:end)                   = nan
    ccf%cwdc_to_litr3c(beg:end)                   = nan
    ccf%litr1_hr(beg:end)                         = nan
    ccf%litr1c_to_soil1c(beg:end)                 = nan
    ccf%litr2_hr(beg:end)                         = nan
    ccf%litr2c_to_soil2c(beg:end)                 = nan
    ccf%litr3_hr(beg:end)                         = nan
    ccf%litr3c_to_soil3c(beg:end)                 = nan
    ccf%soil1_hr(beg:end)                         = nan
    ccf%soil1c_to_soil2c(beg:end)                 = nan
    ccf%soil2_hr(beg:end)                         = nan
    ccf%soil2c_to_soil3c(beg:end)                 = nan
    ccf%soil3_hr(beg:end)                         = nan
    ccf%soil3c_to_soil4c(beg:end)                 = nan
    ccf%soil4_hr(beg:end)                         = nan
    if (use_cn) then
       ccf%dwt_seedc_to_leaf(beg:end)                = nan
       ccf%dwt_seedc_to_deadstem(beg:end)            = nan
       ccf%dwt_conv_cflux(beg:end)                   = nan
       ccf%dwt_prod10c_gain(beg:end)                 = nan
       ccf%dwt_prod100c_gain(beg:end)                = nan
       ccf%dwt_frootc_to_litr1c(beg:end)             = nan
       ccf%dwt_frootc_to_litr2c(beg:end)             = nan
       ccf%dwt_frootc_to_litr3c(beg:end)             = nan
       ccf%dwt_livecrootc_to_cwdc(beg:end)           = nan
       ccf%dwt_deadcrootc_to_cwdc(beg:end)           = nan
       ccf%dwt_closs(beg:end)                        = nan
       ccf%landuseflux(beg:end)                      = nan
       ccf%landuptake(beg:end)                       = nan
       ccf%prod10c_loss(beg:end)                     = nan
       ccf%prod100c_loss(beg:end)                    = nan
       ccf%product_closs(beg:end)                    = nan
    end if
    ccf%lithr(beg:end)                            = nan
    ccf%somhr(beg:end)                            = nan
    ccf%hr(beg:end)                               = nan
    ccf%sr(beg:end)                               = nan
    ccf%er(beg:end)                               = nan
    ccf%litfire(beg:end)                          = nan
    ccf%somfire(beg:end)                          = nan
    ccf%totfire(beg:end)                          = nan
    ccf%nep(beg:end)                              = nan
    ccf%nbp(beg:end)                              = nan
    ccf%nee(beg:end)                              = nan
    ccf%col_cinputs(beg:end)                      = nan
    ccf%col_coutputs(beg:end)                     = nan
    ccf%col_fire_closs(beg:end)                   = nan

    if (use_cn) then
       ccf%cwdc_hr(beg:end)                          = nan
       ccf%cwdc_loss(beg:end)                        = nan
       ccf%litterc_loss(beg:end)                     = nan
    end if

  end subroutine init_column_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nflux_type
!
! !INTERFACE:
  subroutine init_column_nflux_type(beg, end, cnf)
!
! !DESCRIPTION:
! Initialize column nitrogen flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nflux_type), intent(inout):: cnf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cnf%ndep_to_sminn(beg:end))
    allocate(cnf%nfix_to_sminn(beg:end))
    allocate(cnf%m_leafn_to_litr1n(beg:end))
    allocate(cnf%m_leafn_to_litr2n(beg:end))
    allocate(cnf%m_leafn_to_litr3n(beg:end))
    allocate(cnf%m_frootn_to_litr1n(beg:end))
    allocate(cnf%m_frootn_to_litr2n(beg:end))
    allocate(cnf%m_frootn_to_litr3n(beg:end))
    allocate(cnf%m_leafn_storage_to_litr1n(beg:end))
    allocate(cnf%m_frootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_leafn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_frootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_to_cwdn(beg:end))
    allocate(cnf%m_deadstemn_to_cwdn(beg:end))
    allocate(cnf%m_livecrootn_to_cwdn(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn(beg:end))
    allocate(cnf%m_retransn_to_litr1n(beg:end))
    allocate(cnf%hrv_leafn_to_litr1n(beg:end))             
    allocate(cnf%hrv_leafn_to_litr2n(beg:end))             
    allocate(cnf%hrv_leafn_to_litr3n(beg:end))             
    allocate(cnf%hrv_frootn_to_litr1n(beg:end))            
    allocate(cnf%hrv_frootn_to_litr2n(beg:end))            
    allocate(cnf%hrv_frootn_to_litr3n(beg:end))            
    allocate(cnf%hrv_livestemn_to_cwdn(beg:end))           
    allocate(cnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(cnf%hrv_deadstemn_to_prod100n(beg:end))       
    allocate(cnf%hrv_livecrootn_to_cwdn(beg:end))          
    allocate(cnf%hrv_deadcrootn_to_cwdn(beg:end))          
    allocate(cnf%hrv_retransn_to_litr1n(beg:end))          
    allocate(cnf%hrv_leafn_storage_to_litr1n(beg:end))     
    allocate(cnf%hrv_frootn_storage_to_litr1n(beg:end))    
    allocate(cnf%hrv_livestemn_storage_to_litr1n(beg:end)) 
    allocate(cnf%hrv_deadstemn_storage_to_litr1n(beg:end)) 
    allocate(cnf%hrv_livecrootn_storage_to_litr1n(beg:end))
    allocate(cnf%hrv_deadcrootn_storage_to_litr1n(beg:end))
    allocate(cnf%hrv_leafn_xfer_to_litr1n(beg:end))        
    allocate(cnf%hrv_frootn_xfer_to_litr1n(beg:end))       
    allocate(cnf%hrv_livestemn_xfer_to_litr1n(beg:end))    
    allocate(cnf%hrv_deadstemn_xfer_to_litr1n(beg:end))    
    allocate(cnf%hrv_livecrootn_xfer_to_litr1n(beg:end))   
    allocate(cnf%hrv_deadcrootn_xfer_to_litr1n(beg:end))   
    allocate(cnf%m_deadstemn_to_cwdn_fire(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn_fire(beg:end))
    allocate(cnf%m_litr1n_to_fire(beg:end))
    allocate(cnf%m_litr2n_to_fire(beg:end))
    allocate(cnf%m_litr3n_to_fire(beg:end))
    allocate(cnf%m_cwdn_to_fire(beg:end))
    if ( crop_prog )then
       allocate(cnf%grainn_to_litr1n(beg:end))
       allocate(cnf%grainn_to_litr2n(beg:end))
       allocate(cnf%grainn_to_litr3n(beg:end))
       allocate(cnf%livestemn_to_litr1n(beg:end))
       allocate(cnf%livestemn_to_litr2n(beg:end))
       allocate(cnf%livestemn_to_litr3n(beg:end))
    end if
    allocate(cnf%leafn_to_litr1n(beg:end))
    allocate(cnf%leafn_to_litr2n(beg:end))
    allocate(cnf%leafn_to_litr3n(beg:end))
    allocate(cnf%frootn_to_litr1n(beg:end))
    allocate(cnf%frootn_to_litr2n(beg:end))
    allocate(cnf%frootn_to_litr3n(beg:end))
    allocate(cnf%cwdn_to_litr2n(beg:end))
    allocate(cnf%cwdn_to_litr3n(beg:end))
    allocate(cnf%litr1n_to_soil1n(beg:end))
    allocate(cnf%sminn_to_soil1n_l1(beg:end))
    allocate(cnf%litr2n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_l2(beg:end))
    allocate(cnf%litr3n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_l3(beg:end))
    allocate(cnf%soil1n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_s1(beg:end))
    allocate(cnf%soil2n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_s2(beg:end))
    allocate(cnf%soil3n_to_soil4n(beg:end))
    allocate(cnf%sminn_to_soil4n_s3(beg:end))
    allocate(cnf%soil4n_to_sminn(beg:end))
    allocate(cnf%sminn_to_denit_l1s1(beg:end))
    allocate(cnf%sminn_to_denit_l2s2(beg:end))
    allocate(cnf%sminn_to_denit_l3s3(beg:end))
    allocate(cnf%sminn_to_denit_s1s2(beg:end))
    allocate(cnf%sminn_to_denit_s2s3(beg:end))
    allocate(cnf%sminn_to_denit_s3s4(beg:end))
    allocate(cnf%sminn_to_denit_s4(beg:end))
    allocate(cnf%sminn_to_denit_excess(beg:end))
    allocate(cnf%sminn_leached(beg:end))
    allocate(cnf%dwt_seedn_to_leaf(beg:end))
    allocate(cnf%dwt_seedn_to_deadstem(beg:end))
    allocate(cnf%dwt_conv_nflux(beg:end))
    allocate(cnf%dwt_prod10n_gain(beg:end))
    allocate(cnf%dwt_prod100n_gain(beg:end))
    allocate(cnf%dwt_frootn_to_litr1n(beg:end))
    allocate(cnf%dwt_frootn_to_litr2n(beg:end))
    allocate(cnf%dwt_frootn_to_litr3n(beg:end))
    allocate(cnf%dwt_livecrootn_to_cwdn(beg:end))
    allocate(cnf%dwt_deadcrootn_to_cwdn(beg:end))
    allocate(cnf%dwt_nloss(beg:end))
    allocate(cnf%prod10n_loss(beg:end))
    allocate(cnf%prod100n_loss(beg:end))
    allocate(cnf%product_nloss(beg:end))
    allocate(cnf%potential_immob(beg:end))
    allocate(cnf%actual_immob(beg:end))
    allocate(cnf%sminn_to_plant(beg:end))
    allocate(cnf%supplement_to_sminn(beg:end))
    allocate(cnf%gross_nmin(beg:end))
    allocate(cnf%net_nmin(beg:end))
    allocate(cnf%denit(beg:end))
    allocate(cnf%col_ninputs(beg:end))
    allocate(cnf%col_noutputs(beg:end))
    allocate(cnf%col_fire_nloss(beg:end))

    cnf%ndep_to_sminn(beg:end) = nan
    cnf%nfix_to_sminn(beg:end) = nan
    cnf%m_leafn_to_litr1n(beg:end) = nan
    cnf%m_leafn_to_litr2n(beg:end) = nan
    cnf%m_leafn_to_litr3n(beg:end) = nan
    cnf%m_frootn_to_litr1n(beg:end) = nan
    cnf%m_frootn_to_litr2n(beg:end) = nan
    cnf%m_frootn_to_litr3n(beg:end) = nan
    cnf%m_leafn_storage_to_litr1n(beg:end) = nan
    cnf%m_frootn_storage_to_litr1n(beg:end) = nan
    cnf%m_livestemn_storage_to_litr1n(beg:end) = nan
    cnf%m_deadstemn_storage_to_litr1n(beg:end) = nan
    cnf%m_livecrootn_storage_to_litr1n(beg:end) = nan
    cnf%m_deadcrootn_storage_to_litr1n(beg:end) = nan
    cnf%m_leafn_xfer_to_litr1n(beg:end) = nan
    cnf%m_frootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livestemn_xfer_to_litr1n(beg:end) = nan
    cnf%m_deadstemn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livecrootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_deadcrootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livestemn_to_cwdn(beg:end) = nan
    cnf%m_deadstemn_to_cwdn(beg:end) = nan
    cnf%m_livecrootn_to_cwdn(beg:end) = nan
    cnf%m_deadcrootn_to_cwdn(beg:end) = nan
    cnf%m_retransn_to_litr1n(beg:end) = nan
    cnf%hrv_leafn_to_litr1n(beg:end) = nan             
    cnf%hrv_leafn_to_litr2n(beg:end) = nan             
    cnf%hrv_leafn_to_litr3n(beg:end) = nan             
    cnf%hrv_frootn_to_litr1n(beg:end) = nan            
    cnf%hrv_frootn_to_litr2n(beg:end) = nan            
    cnf%hrv_frootn_to_litr3n(beg:end) = nan            
    cnf%hrv_livestemn_to_cwdn(beg:end) = nan           
    cnf%hrv_deadstemn_to_prod10n(beg:end) = nan        
    cnf%hrv_deadstemn_to_prod100n(beg:end) = nan       
    cnf%hrv_livecrootn_to_cwdn(beg:end) = nan          
    cnf%hrv_deadcrootn_to_cwdn(beg:end) = nan          
    cnf%hrv_retransn_to_litr1n(beg:end) = nan          
    cnf%hrv_leafn_storage_to_litr1n(beg:end) = nan     
    cnf%hrv_frootn_storage_to_litr1n(beg:end) = nan    
    cnf%hrv_livestemn_storage_to_litr1n(beg:end) = nan 
    cnf%hrv_deadstemn_storage_to_litr1n(beg:end) = nan 
    cnf%hrv_livecrootn_storage_to_litr1n(beg:end) = nan
    cnf%hrv_deadcrootn_storage_to_litr1n(beg:end) = nan
    cnf%hrv_leafn_xfer_to_litr1n(beg:end) = nan        
    cnf%hrv_frootn_xfer_to_litr1n(beg:end) = nan       
    cnf%hrv_livestemn_xfer_to_litr1n(beg:end) = nan    
    cnf%hrv_deadstemn_xfer_to_litr1n(beg:end) = nan    
    cnf%hrv_livecrootn_xfer_to_litr1n(beg:end) = nan   
    cnf%hrv_deadcrootn_xfer_to_litr1n(beg:end) = nan   
    cnf%m_deadstemn_to_cwdn_fire(beg:end) = nan
    cnf%m_deadcrootn_to_cwdn_fire(beg:end) = nan
    cnf%m_litr1n_to_fire(beg:end) = nan
    cnf%m_litr2n_to_fire(beg:end) = nan
    cnf%m_litr3n_to_fire(beg:end) = nan
    cnf%m_cwdn_to_fire(beg:end) = nan
    if ( crop_prog )then
       cnf%grainn_to_litr1n(beg:end)    = nan
       cnf%grainn_to_litr2n(beg:end)    = nan
       cnf%grainn_to_litr3n(beg:end)    = nan
       cnf%livestemn_to_litr1n(beg:end) = nan
       cnf%livestemn_to_litr2n(beg:end) = nan
       cnf%livestemn_to_litr3n(beg:end) = nan
    end if
    cnf%leafn_to_litr1n(beg:end) = nan
    cnf%leafn_to_litr2n(beg:end) = nan
    cnf%leafn_to_litr3n(beg:end) = nan
    cnf%frootn_to_litr1n(beg:end) = nan
    cnf%frootn_to_litr2n(beg:end) = nan
    cnf%frootn_to_litr3n(beg:end) = nan
    cnf%cwdn_to_litr2n(beg:end) = nan
    cnf%cwdn_to_litr3n(beg:end) = nan
    cnf%litr1n_to_soil1n(beg:end) = nan
    cnf%sminn_to_soil1n_l1(beg:end) = nan
    cnf%litr2n_to_soil2n(beg:end) = nan
    cnf%sminn_to_soil2n_l2(beg:end) = nan
    cnf%litr3n_to_soil3n(beg:end) = nan
    cnf%sminn_to_soil3n_l3(beg:end) = nan
    cnf%soil1n_to_soil2n(beg:end) = nan
    cnf%sminn_to_soil2n_s1(beg:end) = nan
    cnf%soil2n_to_soil3n(beg:end) = nan
    cnf%sminn_to_soil3n_s2(beg:end) = nan
    cnf%soil3n_to_soil4n(beg:end) = nan
    cnf%sminn_to_soil4n_s3(beg:end) = nan
    cnf%soil4n_to_sminn(beg:end) = nan
    cnf%sminn_to_denit_l1s1(beg:end) = nan
    cnf%sminn_to_denit_l2s2(beg:end) = nan
    cnf%sminn_to_denit_l3s3(beg:end) = nan
    cnf%sminn_to_denit_s1s2(beg:end) = nan
    cnf%sminn_to_denit_s2s3(beg:end) = nan
    cnf%sminn_to_denit_s3s4(beg:end) = nan
    cnf%sminn_to_denit_s4(beg:end) = nan
    cnf%sminn_to_denit_excess(beg:end) = nan
    cnf%sminn_leached(beg:end) = nan
    cnf%dwt_seedn_to_leaf(beg:end) = nan
    cnf%dwt_seedn_to_deadstem(beg:end) = nan
    cnf%dwt_conv_nflux(beg:end) = nan
    cnf%dwt_prod10n_gain(beg:end) = nan
    cnf%dwt_prod100n_gain(beg:end) = nan
    cnf%dwt_frootn_to_litr1n(beg:end) = nan
    cnf%dwt_frootn_to_litr2n(beg:end) = nan
    cnf%dwt_frootn_to_litr3n(beg:end) = nan
    cnf%dwt_livecrootn_to_cwdn(beg:end) = nan
    cnf%dwt_deadcrootn_to_cwdn(beg:end) = nan
    cnf%dwt_nloss(beg:end) = nan
    cnf%prod10n_loss(beg:end) = nan
    cnf%prod100n_loss(beg:end) = nan
    cnf%product_nloss(beg:end) = nan
    cnf%potential_immob(beg:end) = nan
    cnf%actual_immob(beg:end) = nan
    cnf%sminn_to_plant(beg:end) = nan
    cnf%supplement_to_sminn(beg:end) = nan
    cnf%gross_nmin(beg:end) = nan
    cnf%net_nmin(beg:end) = nan
    cnf%denit(beg:end) = nan
    cnf%col_ninputs(beg:end) = nan
    cnf%col_noutputs(beg:end) = nan
    cnf%col_fire_nloss(beg:end) = nan

  end subroutine init_column_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_pstate_type
!
! !INTERFACE:
  subroutine init_landunit_pstate_type(beg, end, lps)
!
! !DESCRIPTION:
! Initialize landunit physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (landunit_pstate_type), intent(inout):: lps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(lps%t_building(beg:end))
    allocate(lps%t_building_max(beg:end))
    allocate(lps%t_building_min(beg:end))
    allocate(lps%tk_wall(beg:end,nlevurb))
    allocate(lps%tk_roof(beg:end,nlevurb))
    allocate(lps%tk_improad(beg:end,nlevgrnd))
    allocate(lps%cv_wall(beg:end,nlevurb))
    allocate(lps%cv_roof(beg:end,nlevurb))
    allocate(lps%cv_improad(beg:end,nlevgrnd))
    allocate(lps%thick_wall(beg:end))
    allocate(lps%thick_roof(beg:end))
    allocate(lps%nlev_improad(beg:end))
    allocate(lps%vf_sr(beg:end))
    allocate(lps%vf_wr(beg:end))
    allocate(lps%vf_sw(beg:end))
    allocate(lps%vf_rw(beg:end))
    allocate(lps%vf_ww(beg:end))
    allocate(lps%taf(beg:end))
    allocate(lps%qaf(beg:end))
    allocate(lps%sabs_roof_dir(beg:end,1:numrad))
    allocate(lps%sabs_roof_dif(beg:end,1:numrad))
    allocate(lps%sabs_sunwall_dir(beg:end,1:numrad))
    allocate(lps%sabs_sunwall_dif(beg:end,1:numrad))
    allocate(lps%sabs_shadewall_dir(beg:end,1:numrad))
    allocate(lps%sabs_shadewall_dif(beg:end,1:numrad))
    allocate(lps%sabs_improad_dir(beg:end,1:numrad))
    allocate(lps%sabs_improad_dif(beg:end,1:numrad))
    allocate(lps%sabs_perroad_dir(beg:end,1:numrad))
    allocate(lps%sabs_perroad_dif(beg:end,1:numrad))

    lps%t_building(beg:end) = nan
    lps%t_building_max(beg:end) = nan
    lps%t_building_min(beg:end) = nan
    lps%tk_wall(beg:end,1:nlevurb) = nan
    lps%tk_roof(beg:end,1:nlevurb) = nan
    lps%tk_improad(beg:end,1:nlevgrnd) = nan
    lps%cv_wall(beg:end,1:nlevurb) = nan
    lps%cv_roof(beg:end,1:nlevurb) = nan
    lps%cv_improad(beg:end,1:nlevgrnd) = nan
    lps%cv_improad(beg:end,1:5) = nan
    lps%thick_wall(beg:end) = nan
    lps%thick_roof(beg:end) = nan
    lps%nlev_improad(beg:end) = huge(1)
    lps%vf_sr(beg:end) = nan
    lps%vf_wr(beg:end) = nan
    lps%vf_sw(beg:end) = nan
    lps%vf_rw(beg:end) = nan
    lps%vf_ww(beg:end) = nan
    lps%taf(beg:end) = nan
    lps%qaf(beg:end) = nan
    lps%sabs_roof_dir(beg:end,1:numrad) = nan
    lps%sabs_roof_dif(beg:end,1:numrad) = nan
    lps%sabs_sunwall_dir(beg:end,1:numrad) = nan
    lps%sabs_sunwall_dif(beg:end,1:numrad) = nan
    lps%sabs_shadewall_dir(beg:end,1:numrad) = nan
    lps%sabs_shadewall_dif(beg:end,1:numrad) = nan
    lps%sabs_improad_dir(beg:end,1:numrad) = nan
    lps%sabs_improad_dif(beg:end,1:numrad) = nan
    lps%sabs_perroad_dir(beg:end,1:numrad) = nan
    lps%sabs_perroad_dif(beg:end,1:numrad) = nan

  end subroutine init_landunit_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_eflux_type
!
! !INTERFACE:
  subroutine init_landunit_eflux_type(beg, end, lef)
!
! !DESCRIPTION: 
! Initialize landunit energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end 
    type (landunit_eflux_type), intent(inout):: lef 
!
! !REVISION HISTORY:
! Created by Keith Oleson
!
!EOP
!------------------------------------------------------------------------

    allocate(lef%eflx_traffic(beg:end))
    allocate(lef%eflx_traffic_factor(beg:end))
    allocate(lef%eflx_wasteheat(beg:end))
    allocate(lef%eflx_heat_from_ac(beg:end))

    lef%eflx_traffic(beg:end) = nan
    lef%eflx_traffic_factor(beg:end) = nan
    lef%eflx_wasteheat(beg:end) = nan
    lef%eflx_heat_from_ac(beg:end) = nan

  end subroutine init_landunit_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_dgvstate_type
!
! !INTERFACE:
  subroutine init_gridcell_dgvstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell DGVM variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_dgvstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gps%agdd20(beg:end))
    allocate(gps%tmomin20(beg:end))
    allocate(gps%t10min(beg:end))

    gps%agdd20(beg:end) = nan
    gps%tmomin20(beg:end) = nan
    gps%t10min(beg:end) = nan

  end subroutine init_gridcell_dgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_pstate_type
!
! !INTERFACE:
  subroutine init_gridcell_pstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_pstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    
    
    !allocate(gps%bcphiwet2t(beg:end,1:2))
    !allocate(gps%bcphidry2t(beg:end,1:2))
    !allocate(gps%bcphodry2t(beg:end,1:2))
    !allocate(gps%ocphiwet2t(beg:end,1:2))
    !allocate(gps%ocphidry2t(beg:end,1:2))
    !allocate(gps%ocphodry2t(beg:end,1:2))
    !allocate(gps%dstx01wd2t(beg:end,1:2))
    !allocate(gps%dstx01dd2t(beg:end,1:2))
    !allocate(gps%dstx02wd2t(beg:end,1:2))
    !allocate(gps%dstx02dd2t(beg:end,1:2))
    !allocate(gps%dstx03wd2t(beg:end,1:2))
    !allocate(gps%dstx03dd2t(beg:end,1:2))
    !allocate(gps%dstx04wd2t(beg:end,1:2))
    !allocate(gps%dstx04dd2t(beg:end,1:2))
    
    !gps%bcphiwet2t(beg:end,1:2) = nan
    !gps%bcphidry2t(beg:end,1:2) = nan
    !gps%bcphodry2t(beg:end,1:2) = nan
    !gps%ocphiwet2t(beg:end,1:2) = nan
    !gps%ocphidry2t(beg:end,1:2) = nan
    !gps%ocphodry2t(beg:end,1:2) = nan
    !gps%dstx01wd2t(beg:end,1:2) = nan
    !gps%dstx01dd2t(beg:end,1:2) = nan
    !gps%dstx02wd2t(beg:end,1:2) = nan
    !gps%dstx02dd2t(beg:end,1:2) = nan
    !gps%dstx03wd2t(beg:end,1:2) = nan
    !gps%dstx03dd2t(beg:end,1:2) = nan
    !gps%dstx04wd2t(beg:end,1:2) = nan
    !gps%dstx04dd2t(beg:end,1:2) = nan

  end subroutine init_gridcell_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_efstate_type
!
! !INTERFACE:
  subroutine init_gridcell_efstate_type(beg, end, gve)
!
! !DESCRIPTION:
! Initialize gridcell isoprene emission factor variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_efstate_type), intent(inout) :: gve
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein (heald)
!
!EOP
!------------------------------------------------------------------------

    allocate(gve%efisop(6,beg:end))
    gve%efisop(:,beg:end) = nan

  end subroutine init_gridcell_efstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wflux_type
!
! !INTERFACE:
  subroutine init_gridcell_wflux_type(beg, end, gwf)
!
! !DESCRIPTION:
! Initialize gridcell water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wflux_type), intent(inout):: gwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gwf%qflx_runoffg(beg:end))
    allocate(gwf%qflx_snwcp_iceg(beg:end))
    allocate(gwf%qflx_liq_dynbal(beg:end))
    allocate(gwf%qflx_ice_dynbal(beg:end))
    allocate(gwf%qflx_floodg(beg:end))

    gwf%qflx_runoffg(beg:end) = 0._r8
    gwf%qflx_snwcp_iceg(beg:end) = 0._r8
    gwf%qflx_liq_dynbal(beg:end) = nan
    gwf%qflx_ice_dynbal(beg:end) = nan
    gwf%qflx_floodg(beg:end) = 0._r8 !rtm_flood: initialize to zero for 1st time step instead of nan

  end subroutine init_gridcell_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_eflux_type
!
! !INTERFACE:
  subroutine init_gridcell_eflux_type(beg, end, gef)
!
! !DESCRIPTION:
! Initialize gridcell energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_eflux_type), intent(inout):: gef
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(gef%eflx_sh_totg(beg:end))
    allocate(gef%eflx_dynbal(beg:end))

    gef%eflx_sh_totg(beg:end) = nan
    gef%eflx_dynbal(beg:end) = nan

  end subroutine init_gridcell_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wstate_type
!
! !INTERFACE:
  subroutine init_gridcell_wstate_type(beg, end, gws)
!
! !DESCRIPTION:
! Initialize gridcell water state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wstate_type), intent(inout):: gws
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(gws%gc_liq1(beg:end))
    allocate(gws%gc_liq2(beg:end))
    allocate(gws%gc_ice1(beg:end))     
    allocate(gws%gc_ice2(beg:end))    

    gws%gc_liq1(beg:end) = nan
    gws%gc_liq2(beg:end) = nan
    gws%gc_ice1(beg:end) = nan     
    gws%gc_ice2(beg:end) = nan    

  end subroutine init_gridcell_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_estate_type
!
! !INTERFACE:
  subroutine init_gridcell_estate_type(beg, end, ges)
!
! !DESCRIPTION:
! Initialize gridcell energy state variables     
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_estate_type), intent(inout):: ges
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(ges%gc_heat1(beg:end))     
    allocate(ges%gc_heat2(beg:end))    

    ges%gc_heat1(beg:end) = nan     
    ges%gc_heat2(beg:end) = nan    

  end subroutine init_gridcell_estate_type

end module clmtypeInitMod
