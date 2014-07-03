module clmtypeInitMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Allocate clmtype components and initialize them to signaling NaN.
  !
  ! !USES:
  use clmtype
  use shr_kind_mod, only: &
       r8 => shr_kind_r8
  use shr_infnan_mod, only: &
       nan => shr_infnan_nan, assignment(=)
  use clm_varcon, only : &
       spval, ispval, max_lunit
  use clm_varpar  , only : &
       maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
       numpft, ndst, nlevurb, nlevsoi, nlevdecomp, nlevdecomp_full, &
       ndecomp_cascade_transitions, ndecomp_pools, nlevcan, &
       nlayer, nlayert, crop_prog, ngases
  use shr_megan_mod, &
       only: shr_megan_megcomps_n, shr_megan_mechcomps_n
  use seq_drydep_mod, &
       only:  n_drydep, drydep_method, DD_XLND
  use decompMod, &
       only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initClmtype
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type

  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_decomp_cascade_constants

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
  private :: init_column_ch4_type
  private :: init_column_nflux_type

  private :: init_landunit_pstate_type
  private :: init_landunit_eflux_type

  private :: init_gridcell_efstate_type
  private :: init_gridcell_wflux_type
  private :: init_gridcell_ch4_type
  
  real(r8) :: nanr 
  !----------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initClmtype(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize clmtype components to signaling nan
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! LOCAL VARAIBLES:
    character(len=32), parameter :: subname = "initClmtype"
    !------------------------------------------------------------------------

    nanr = nan

    ! Determine necessary indices

    call init_pft_type     (bounds%begp, bounds%endp, pft)
    call init_column_type  (bounds%begc, bounds%endc, col)
    call init_landunit_type(bounds%begl, bounds%endl, lun)
    call init_gridcell_type(bounds%begg, bounds%endg, grc)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    call init_decomp_cascade_constants()

    ! pft DGVM-specific ecophysiological constants

    call init_pft_DGVMecophys_constants()

    ! energy balance structures (all levels)

    call init_energy_balance_type(bounds%begp, bounds%endp, pebal)
    call init_energy_balance_type(bounds%begc, bounds%endc, cebal)

    ! water balance structures (all levels)

    call init_water_balance_type(bounds%begp, bounds%endp, pwbal)
    call init_water_balance_type(bounds%begc, bounds%endc, cwbal)

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(bounds%begp, bounds%endp, pcbal)
    call init_carbon_balance_type(bounds%begc, bounds%endc, ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(bounds%begp, bounds%endp, pnbal)
    call init_nitrogen_balance_type(bounds%begc, bounds%endc, cnbal)

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(bounds%begp, bounds%endp, pps)
    call init_pft_pstate_type(bounds%begc, bounds%endc, pps_a)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(bounds%begp, bounds%endp, pepv)

    ! pft photosynthesis relevant variables
    call init_pft_psynstate_type(bounds%begp, bounds%endp, ppsyns)

    ! pft DGVM state variables at pft level 

    call init_pft_pdgvstate_type(bounds%begp, bounds%endp, pdgvs)
    call init_pft_vstate_type(bounds%begp, bounds%endp, pvs)

    ! pft energy state variables at the pft level 

    call init_pft_estate_type(bounds%begp, bounds%endp, pes)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(bounds%begp, bounds%endp, pws)
    call init_pft_wstate_type(bounds%begc, bounds%endc, pws_a)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(bounds%begp, bounds%endp, pcs)
    call init_pft_cstate_type(bounds%begc, bounds%endc, pcs_a)

    call init_pft_cstate_type(bounds%begp, bounds%endp, pc13s)
    call init_pft_cstate_type(bounds%begc, bounds%endc, pc13s_a)

    call init_pft_cstate_type(bounds%begp, bounds%endp, pc14s)
    call init_pft_cstate_type(bounds%begc, bounds%endc, pc14s_a)

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(bounds%begp, bounds%endp, pns)
    call init_pft_nstate_type(bounds%begc, bounds%endc, pns_a)

    ! pft energy flux variables at pft level

    call init_pft_eflux_type(bounds%begp, bounds%endp, pef)

    ! pft momentum flux variables at pft level 

    call init_pft_mflux_type(bounds%begp, bounds%endp, pmf)

    ! pft water flux variables

    call init_pft_wflux_type(bounds%begp, bounds%endp, pwf)
    call init_pft_wflux_type(bounds%begc, bounds%endc, pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(bounds%begp, bounds%endp, pcf)
    call init_pft_cflux_type(bounds%begc, bounds%endc, pcf_a)

    call init_pft_cflux_type(bounds%begp, bounds%endp, pc13f)
    call init_pft_cflux_type(bounds%begc, bounds%endc, pc13f_a)

    call init_pft_cflux_type(bounds%begp, bounds%endp, pc14f)
    call init_pft_cflux_type(bounds%begc, bounds%endc, pc14f_a)

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(bounds%begp, bounds%endp, pnf)
    call init_pft_nflux_type(bounds%begc, bounds%endc, pnf_a)

    ! pft VOC flux variables at pft level

    call init_pft_vflux_type(bounds%begp, bounds%endp, pvf)

    ! gridcell VOC emission factors (heald, 05/06)

    call init_gridcell_efstate_type(bounds%begg, bounds%endg, gve)

    ! pft dust flux variables at pft level 

    call init_pft_dflux_type(bounds%begp, bounds%endp, pdf)

    ! pft dry dep velocity variables at pft level 

    call init_pft_depvd_type(bounds%begp, bounds%endp, pdd)

    ! column physical state variables at column level 

    call init_column_pstate_type(bounds%begc, bounds%endc, cps)

    ! column energy state variables at column level

    call init_column_estate_type(bounds%begc, bounds%endc, ces)

    ! column water state variables at column level 

    call init_column_wstate_type(bounds%begc, bounds%endc, cws)

    ! column carbon state variables at column level

    call init_column_cstate_type(bounds%begc, bounds%endc, ccs)

    call init_column_cstate_type(bounds%begc, bounds%endc, cc13s)

    call init_column_cstate_type(bounds%begc, bounds%endc, cc14s)

    ! column nitrogen state variables at column level

    call init_column_nstate_type(bounds%begc, bounds%endc, cns)

    ! column energy flux variables at column level 

    call init_column_eflux_type(bounds%begc, bounds%endc, cef)

    ! column water flux variables at column level 

    call init_column_wflux_type(bounds%begc, bounds%endc, cwf)

    ! column carbon flux variables at column level

    call init_column_cflux_type(bounds%begc, bounds%endc, ccf)

    call init_column_cflux_type(bounds%begc, bounds%endc, cc13f)

    call init_column_cflux_type(bounds%begc, bounds%endc, cc14f)

    ! column CH4 flux variables at column level

    call init_column_ch4_type(bounds%begc, bounds%endc, cch4)

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(bounds%begc, bounds%endc, cnf)

    ! land unit physical state variables

    call init_landunit_pstate_type(bounds%begl, bounds%endl, lps)

    ! land unit energy flux variables 

    call init_landunit_eflux_type(bounds%begl, bounds%endl, lef)

    ! gridcell: water flux variables

    call init_gridcell_wflux_type(bounds%begg, bounds%endg, gwf)

    ! gridcell: energy flux variables

    call init_gridcell_eflux_type(bounds%begg, bounds%endg, gef)

    ! gridcell: water state variables

    call init_gridcell_wstate_type(bounds%begg, bounds%endg, gws)

    ! gridcell: energy state variables

    call init_gridcell_estate_type(bounds%begg, bounds%endg, ges)

    ! gridcell: physical state variables
    
    call init_gridcell_pstate_type(bounds%begg, bounds%endg, gps)

    ! gridcell: ch4 variables

    call init_gridcell_ch4_type(bounds%begg, bounds%endg, gch4)

  end subroutine initClmtype

  !------------------------------------------------------------------------
  subroutine init_pft_type (beg, end, pft)
    !
    ! !DESCRIPTION:
    ! Initialize components of pft_type structure
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: pft
    !------------------------------------------------------------------------

    allocate(pft%gridcell(beg:end),pft%wtgcell(beg:end))
    allocate(pft%landunit(beg:end),pft%wtlunit(beg:end))
    allocate(pft%column  (beg:end),pft%wtcol  (beg:end))

    allocate(pft%itype(beg:end))
    allocate(pft%mxy(beg:end))
    allocate(pft%active(beg:end))

  end subroutine init_pft_type

  !------------------------------------------------------------------------
  subroutine init_column_type (beg, end, c)
    !
    ! !DESCRIPTION:
    ! Initialize components of column_type structure
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(column_type), intent(inout):: c
    !------------------------------------------------------------------------

    allocate(col%gridcell(beg:end),col%wtgcell(beg:end))
    allocate(col%landunit(beg:end),col%wtlunit(beg:end))

    allocate(col%pfti(beg:end),col%pftf(beg:end),col%npfts(beg:end))

    allocate(col%itype(beg:end))
    allocate(col%active(beg:end))

  end subroutine init_column_type

  !------------------------------------------------------------------------
  subroutine init_landunit_type (beg, end,l)
    !
    ! !DESCRIPTION:
    ! Initialize components of landunit_type structure
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(landunit_type), intent(inout):: l
    !------------------------------------------------------------------------

    allocate(lun%gridcell(beg:end),lun%wtgcell(beg:end))

    allocate(lun%coli(beg:end),lun%colf(beg:end),lun%ncolumns(beg:end))
    allocate(lun%pfti(beg:end),lun%pftf(beg:end),lun%npfts   (beg:end))

    allocate(lun%itype(beg:end))
    allocate(lun%ifspecial(beg:end))
    allocate(lun%lakpoi(beg:end))
    allocate(lun%urbpoi(beg:end))
    allocate(lun%glcmecpoi(beg:end))
    lun%glcmecpoi(:)=.false.
    allocate(lun%active(beg:end))

    allocate(lun%canyon_hwr(beg:end))
    lun%canyon_hwr(:)=nanr
    allocate(lun%wtroad_perv(beg:end))
    lun%wtroad_perv(:)=nanr
    allocate(lun%ht_roof(beg:end))
    lun%ht_roof(:)=nanr
    allocate(lun%wtlunit_roof(beg:end))
    lun%wtlunit_roof(:)=nanr
    allocate(lun%z_0_town(beg:end))
    lun%z_0_town(:)=nanr
    allocate(lun%z_d_town(beg:end))
    lun%z_d_town(:)=nanr

  end subroutine init_landunit_type

  !------------------------------------------------------------------------
  subroutine init_gridcell_type (beg, end,grc)
    !
    ! !DESCRIPTION:
    ! Initialize components of gridcell_type structure
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: grc
    !------------------------------------------------------------------------

    allocate(grc%gindex(beg:end))
    allocate(grc%area(beg:end))
    allocate(grc%lat(beg:end))
    allocate(grc%lon(beg:end))
    allocate(grc%latdeg(beg:end))
    allocate(grc%londeg(beg:end))

    allocate(grc%landunit_indices(1:max_lunit, beg:end))
    grc%landunit_indices(:,:) = ispval

    allocate(grc%gris_mask(beg:end))
    grc%gris_mask(:)=nanr
    allocate(grc%gris_area(beg:end))
    grc%gris_area(:)=nanr
    allocate(grc%aais_mask(beg:end))
    grc%aais_mask(:)=nanr
    allocate(grc%aais_area(beg:end))
    grc%aais_area(:)=nanr
    allocate(grc%icemask(beg:end))
    grc%icemask(:)=nanr      
    allocate(grc%tws(beg:end))
    grc%tws(:)=nanr

  end subroutine init_gridcell_type

!------------------------------------------------------------------------
  subroutine init_energy_balance_type(beg, end, ebal)
    !
    ! !DESCRIPTION:
    ! Initialize energy balance variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(energy_balance_type), intent(inout):: ebal
    !------------------------------------------------------------------------

    allocate(ebal%errsoi(beg:end))
    ebal%errsoi(:)=nanr
    allocate(ebal%errseb(beg:end))
    ebal%errseb(:)=nanr
    allocate(ebal%errsol(beg:end))
    ebal%errsol(:)=nanr
    allocate(ebal%errlon(beg:end))
    ebal%errlon(:)=nanr

  end subroutine init_energy_balance_type

  !------------------------------------------------------------------------
  subroutine init_water_balance_type(beg, end, wbal)
    !
    ! !DESCRIPTION:
    ! Initialize water balance variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(water_balance_type), intent(inout):: wbal
    !------------------------------------------------------------------------

    allocate(wbal%begwb(beg:end))
    wbal%begwb(:)=nanr
    allocate(wbal%endwb(beg:end))
    wbal%endwb(:)=nanr
    allocate(wbal%errh2o(beg:end))
    wbal%errh2o(:)=nanr

  end subroutine init_water_balance_type

!------------------------------------------------------------------------
  subroutine init_carbon_balance_type(beg, end, cbal)
    !
    ! !DESCRIPTION:
    ! Initialize carbon balance variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(carbon_balance_type), intent(inout):: cbal
    !------------------------------------------------------------------------

    allocate(cbal%begcb(beg:end))
    cbal%begcb(:)=nanr
    allocate(cbal%endcb(beg:end))
    cbal%endcb(:)=nanr
    allocate(cbal%errcb(beg:end))
    cbal%errcb(:)=nanr

  end subroutine init_carbon_balance_type

  !------------------------------------------------------------------------
  subroutine init_nitrogen_balance_type(beg, end, nbal)
    !
    ! !DESCRIPTION:
    ! Initialize nitrogen balance variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(nitrogen_balance_type), intent(inout):: nbal
   !------------------------------------------------------------------------

    allocate(nbal%begnb(beg:end))
    nbal%begnb(:)=nanr
    allocate(nbal%endnb(beg:end))
    nbal%endnb(:)=nanr
    allocate(nbal%errnb(beg:end))
    nbal%errnb(:)=nanr

  end subroutine init_nitrogen_balance_type

  !------------------------------------------------------------------------
  subroutine init_pft_ecophys_constants()
    !
    ! !DESCRIPTION:
    ! Initialize pft physical state
    !
    ! !ARGUMENTS:
    implicit none
    !------------------------------------------------------------------------

    allocate(pftcon%noveg(0:numpft))
    pftcon%noveg(0:numpft)=huge(1)
    allocate(pftcon%tree(0:numpft))
    pftcon%tree(0:numpft)=huge(1)
    allocate(pftcon%smpso(0:numpft))
    pftcon%smpso(0:numpft)=nanr
    allocate(pftcon%smpsc(0:numpft))
    pftcon%smpsc(0:numpft)=nanr
    allocate(pftcon%fnitr(0:numpft))
    pftcon%fnitr(0:numpft)=nanr
    allocate(pftcon%foln(0:numpft))
    pftcon%foln(0:numpft)=nanr
    allocate(pftcon%dleaf(0:numpft))
    pftcon%dleaf(0:numpft)=nanr
    allocate(pftcon%c3psn(0:numpft))
    pftcon%c3psn(0:numpft)=nanr
    allocate(pftcon%xl(0:numpft))
    pftcon%xl(0:numpft)=nanr
    allocate(pftcon%rhol(0:numpft,numrad))
    pftcon%rhol(0:numpft,numrad)=nanr
    allocate(pftcon%rhos(0:numpft,numrad))
    pftcon%rhos(0:numpft,numrad)=nanr
    allocate(pftcon%taul(0:numpft,numrad))
    pftcon%taul(0:numpft,numrad)=nanr
    allocate(pftcon%taus(0:numpft,numrad))
    pftcon%taus(0:numpft,numrad)=nanr
    allocate(pftcon%z0mr(0:numpft))
    pftcon%z0mr(0:numpft)=nanr
    allocate(pftcon%displar(0:numpft))
    pftcon%displar(0:numpft)=nanr
    allocate(pftcon%roota_par(0:numpft))
    pftcon%roota_par(0:numpft)=nanr
    allocate(pftcon%rootb_par(0:numpft))
    pftcon%rootb_par(0:numpft)=nanr
    allocate(pftcon%slatop(0:numpft))
    pftcon%slatop(0:numpft)=nanr
    allocate(pftcon%dsladlai(0:numpft))
    pftcon%dsladlai(0:numpft)=nanr
    allocate(pftcon%leafcn(0:numpft))
    pftcon%leafcn(0:numpft)=nanr
    allocate(pftcon%flnr(0:numpft))
    pftcon%flnr(0:numpft)=nanr
    allocate(pftcon%woody(0:numpft))
    pftcon%woody(0:numpft)=nanr
    allocate(pftcon%lflitcn(0:numpft))
    pftcon%lflitcn(0:numpft)=nanr
    allocate(pftcon%frootcn(0:numpft))
    pftcon%frootcn(0:numpft)=nanr
    allocate(pftcon%livewdcn(0:numpft))
    pftcon%livewdcn(0:numpft)=nanr
    allocate(pftcon%deadwdcn(0:numpft))
    pftcon%deadwdcn(0:numpft)=nanr
    allocate(pftcon%graincn(0:numpft))
    pftcon%graincn(0:numpft)=nanr
    allocate(pftcon%froot_leaf(0:numpft))
    pftcon%froot_leaf(0:numpft)=nanr
    allocate(pftcon%stem_leaf(0:numpft))
    pftcon%stem_leaf(0:numpft)=nanr
    allocate(pftcon%croot_stem(0:numpft))
    pftcon%croot_stem(0:numpft)=nanr
    allocate(pftcon%flivewd(0:numpft))
    pftcon%flivewd(0:numpft)=nanr
    allocate(pftcon%fcur(0:numpft))
    pftcon%fcur(0:numpft)=nanr
    allocate(pftcon%lf_flab(0:numpft))
    pftcon%lf_flab(0:numpft)=nanr
    allocate(pftcon%lf_fcel(0:numpft))
    pftcon%lf_fcel(0:numpft)=nanr
    allocate(pftcon%lf_flig(0:numpft))
    pftcon%lf_flig(0:numpft)=nanr
    allocate(pftcon%fr_flab(0:numpft))
    pftcon%fr_flab(0:numpft)=nanr
    allocate(pftcon%fr_fcel(0:numpft))
    pftcon%fr_fcel(0:numpft)=nanr
    allocate(pftcon%fr_flig(0:numpft))
    pftcon%fr_flig(0:numpft)=nanr
    allocate(pftcon%leaf_long(0:numpft))
    pftcon%leaf_long(0:numpft)=nanr
    allocate(pftcon%evergreen(0:numpft))
    pftcon%evergreen(0:numpft)=nanr
    allocate(pftcon%stress_decid(0:numpft))
    pftcon%stress_decid(0:numpft)=nanr
    allocate(pftcon%season_decid(0:numpft))
    pftcon%season_decid(0:numpft)=nanr
    allocate(pftcon%dwood(0:numpft))
    pftcon%dwood(0:numpft)=nanr
    allocate(pftcon%rootprof_beta(0:numpft))
    pftcon%rootprof_beta(0:numpft)=nanr
    allocate(pftcon%fertnitro(0:numpft))
    pftcon%fertnitro(0:numpft)=nanr
    allocate(pftcon%fleafcn(0:numpft))
    pftcon%fleafcn(0:numpft)=nanr
    allocate(pftcon%ffrootcn(0:numpft))
    pftcon%ffrootcn(0:numpft)=nanr
    allocate(pftcon%fstemcn(0:numpft))
    pftcon%fstemcn(0:numpft)=nanr

  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
  subroutine init_decomp_cascade_constants()
    !
    ! !DESCRIPTION:
    ! Initialize decomposition cascade state
    !
    ! !ARGUMENTS:
    implicit none
    !------------------------------------------------------------------------

    !-- properties of each pathway along decomposition cascade 
    allocate(decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions))

    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_litter(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_soil(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cwd(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_stock(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_metabolic(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cellulose(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_lignin(0:ndecomp_pools))
    allocate(decomp_cascade_con%spinup_factor(0:ndecomp_pools))

    !-- properties of each pathway along decomposition cascade 
    decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions) = ''
    decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions) = 0
    decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions) = 0

    !-- properties of each decomposing pool
    decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools) = ''
    decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools) = ''
    decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools) = ''
    decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools) = ''
    decomp_cascade_con%is_litter(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_soil(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_cwd(0:ndecomp_pools) = .false.
    decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools) = nan
    decomp_cascade_con%initial_stock(0:ndecomp_pools) = nan
    decomp_cascade_con%is_metabolic(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_cellulose(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_lignin(0:ndecomp_pools) = .false.
    decomp_cascade_con%spinup_factor(0:ndecomp_pools) = nan

  end subroutine init_decomp_cascade_constants

  !------------------------------------------------------------------------
  subroutine init_pft_DGVMecophys_constants()
    !
    ! !DESCRIPTION:
    ! Initialize pft physical state
    !
    ! !ARGUMENTS:
    implicit none
    !------------------------------------------------------------------------

    allocate(dgv_pftcon%crownarea_max(0:numpft))
    dgv_pftcon%crownarea_max(0:numpft)=nanr
    allocate(dgv_pftcon%tcmin(0:numpft))
    dgv_pftcon%tcmin(0:numpft)=nanr
    allocate(dgv_pftcon%tcmax(0:numpft))
    dgv_pftcon%tcmax(0:numpft)=nanr
    allocate(dgv_pftcon%gddmin(0:numpft))
    dgv_pftcon%gddmin(0:numpft)=nanr
    allocate(dgv_pftcon%twmax(0:numpft))
    dgv_pftcon%twmax(0:numpft)=nanr
    allocate(dgv_pftcon%reinickerp(0:numpft))
    dgv_pftcon%reinickerp(0:numpft)=nanr
    allocate(dgv_pftcon%allom1(0:numpft))
    dgv_pftcon%allom1(0:numpft)=nanr
    allocate(dgv_pftcon%allom2(0:numpft))
    dgv_pftcon%allom2(0:numpft)=nanr
    allocate(dgv_pftcon%allom3(0:numpft))
    dgv_pftcon%allom3(0:numpft)=nanr       

  end subroutine init_pft_DGVMecophys_constants

  !------------------------------------------------------------------------
  subroutine init_pft_pstate_type(beg, end, pps)
    !
    ! !DESCRIPTION:
    ! Initialize pft physical state
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_pstate_type), intent(inout):: pps
    !------------------------------------------------------------------------

    allocate(pps%prec10(beg:end))
    pps%prec10(:)= nanr
    allocate(pps%prec60(beg:end))
    pps%prec60(:)= nanr
    allocate(pps%frac_veg_nosno(beg:end))
    pps%frac_veg_nosno(:)= huge(1)
    allocate(pps%frac_veg_nosno_alb(beg:end))
    pps%frac_veg_nosno_alb(:)= 0
    allocate(pps%emv(beg:end))
    pps%emv(:)= nanr
    allocate(pps%z0mv(beg:end))
    pps%z0mv(:)= nanr
    allocate(pps%z0hv(beg:end))
    pps%z0hv(:)= nanr
    allocate(pps%z0qv(beg:end))
    pps%z0qv(:)= nanr
    allocate(pps%rootfr(beg:end,1:nlevgrnd))
    pps%rootfr(:,:)= spval
    allocate(pps%rootr(beg:end,1:nlevgrnd))
    pps%rootr(:,:)= spval
    allocate(pps%rresis(beg:end,1:nlevgrnd))
    pps%rresis(:,:)= spval
    allocate(pps%dewmx(beg:end))
    pps%dewmx(:)= nanr
    allocate(pps%rssun(beg:end))
    pps%rssun(:)= nanr
    allocate(pps%rssha(beg:end))
    pps%rssha(:)= nanr
    allocate(pps%rhal(beg:end))
    pps%rhal(:)= nanr
    allocate(pps%vpdal(beg:end))
    pps%vpdal(:)= nanr
    allocate(pps%rssun_z(beg:end,1:nlevcan))
    pps%rssun_z(:,:)= nanr
    allocate(pps%rssha_z(beg:end,1:nlevcan))
    pps%rssha_z(:,:)= nanr
    allocate(pps%laisun(beg:end))
    pps%laisun(:)= nanr
    allocate(pps%laisha(beg:end))
    pps%laisha(:)= nanr
    allocate(pps%laisun_z(beg:end,1:nlevcan))
    pps%laisun_z(:,:)= nanr
    allocate(pps%laisha_z(beg:end,1:nlevcan))
    pps%laisha_z(:,:)= nanr
    allocate(pps%btran(beg:end))
    pps%btran(:)= spval
    allocate(pps%btran2(beg:end))
    pps%btran2(:)= spval
    allocate(pps%fsun(beg:end))
    pps%fsun(:)= spval
    allocate(pps%tlai(beg:end))
    pps%tlai(:)= 0._r8
    allocate(pps%tsai(beg:end))
    pps%tsai(:)= 0._r8
    allocate(pps%elai(beg:end))
    pps%elai(:)= 0._r8
    allocate(pps%esai(beg:end))
    pps%esai(:)= 0._r8
    allocate(pps%fwet(beg:end))
    pps%fwet(:)= nanr
    allocate(pps%fdry(beg:end))
    pps%fdry(:)= nanr
    allocate(pps%dt_veg(beg:end))
    pps%dt_veg(:)= nanr
    allocate(pps%htop(beg:end))
    pps%htop(:)= 0._r8
    allocate(pps%hbot(beg:end))
    pps%hbot(:)= 0._r8
    allocate(pps%z0m(beg:end))
    pps%z0m(:)= nanr
    allocate(pps%displa(beg:end))
    pps%displa(:)= nanr
    allocate(pps%albd(beg:end,1:numrad))
    pps%albd(:,:)= nanr
    allocate(pps%albi(beg:end,1:numrad))
    pps%albi(:,:)= nanr
    allocate(pps%fabd(beg:end,1:numrad))
    pps%fabd(:,:)= nanr
    allocate(pps%fabd_sun(beg:end,1:numrad))
    pps%fabd_sun(:,:)= nanr
    allocate(pps%fabd_sha(beg:end,1:numrad))
    pps%fabd_sha(:,:)= nanr
    allocate(pps%fabi(beg:end,1:numrad))
    pps%fabi(:,:)= nanr
    allocate(pps%fabi_sun(beg:end,1:numrad))
    pps%fabi_sun(:,:)= nanr
    allocate(pps%fabi_sha(beg:end,1:numrad))
    pps%fabi_sha(:,:)= nanr
    allocate(pps%ftdd(beg:end,1:numrad))
    pps%ftdd(:,:)= nanr
    allocate(pps%ftid(beg:end,1:numrad))
    pps%ftid(:,:)= nanr
    allocate(pps%ftii(beg:end,1:numrad))
    pps%ftii(:,:)= nanr
    allocate(pps%vcmaxcintsun(beg:end))
    pps%vcmaxcintsun(:)= nanr
    allocate(pps%vcmaxcintsha(beg:end))
    pps%vcmaxcintsha(:)= nanr
    allocate(pps%ncan(beg:end))
    pps%ncan(:)= 0
    allocate(pps%nrad(beg:end))
    pps%nrad(:)= 0
    allocate(pps%fabd_sun_z(beg:end,1:nlevcan))
    pps%fabd_sun_z(:,:)= 0._r8
    allocate(pps%fabd_sha_z(beg:end,1:nlevcan))
    pps%fabd_sha_z(:,:)= 0._r8
    allocate(pps%fabi_sun_z(beg:end,1:nlevcan))
    pps%fabi_sun_z(:,:)= 0._r8
    allocate(pps%fabi_sha_z(beg:end,1:nlevcan))
    pps%fabi_sha_z(:,:)= 0._r8
    allocate(pps%fsun_z(beg:end,1:nlevcan))
    pps%fsun_z(:,:)= 0._r8
    allocate(pps%tlai_z(beg:end,1:nlevcan))
    pps%tlai_z(:,:)= 0._r8
    allocate(pps%tsai_z(beg:end,1:nlevcan))
    pps%tsai_z(:,:)= 0._r8
    allocate(pps%u10(beg:end))
    pps%u10(:)= nanr
    allocate(pps%u10_clm(beg:end))
    pps%u10_clm(:)= nanr
    allocate(pps%va(beg:end))
    pps%va(:)= nanr
    allocate(pps%fv(beg:end))
    pps%fv(:)= nanr
    allocate(pps%ram1(beg:end))
    pps%ram1(:)= nanr
    allocate(pps%burndate(beg:end))
    pps%burndate(:)= ispval
    allocate(pps%ram1_lake(beg:end))
    pps%ram1_lake(:)= nanr
    allocate(pps%rh_leaf(beg:end))
    pps%rh_leaf(:)= spval
    allocate(pps%rhaf(beg:end))
    pps%rhaf(:)= spval
    allocate(pps%gddmaturity(beg:end))
    pps%gddmaturity(:)= spval
    allocate(pps%huileaf(beg:end))
    pps%huileaf(:)= nanr
    allocate(pps%huigrain(beg:end))
    pps%huigrain(:)= nanr
    allocate(pps%gddplant(beg:end))
    pps%gddplant(:)= spval
    allocate(pps%gddtsoi(beg:end))
    pps%gddtsoi(:)= spval
    allocate(pps%croplive(beg:end))
    pps%croplive(:)= .false.
    allocate(pps%peaklai(beg:end))
    pps%peaklai(:)= 0
    allocate(pps%aleafi(beg:end))
    pps%aleafi(:)= nanr
    allocate(pps%astemi(beg:end))
    pps%astemi(:)= nanr
    allocate(pps%aleaf(beg:end))
    pps%aleaf(:)= nanr
    allocate(pps%astem(beg:end))
    pps%astem(:)= nanr
    allocate(pps%gdd0(beg:end))
    pps%gdd0(:)= spval
    allocate(pps%gdd8(beg:end))
    pps%gdd8(:)= spval
    allocate(pps%gdd10(beg:end))
    pps%gdd10(:)= spval
    allocate(pps%gdd020(beg:end))
    pps%gdd020(:)= spval
    allocate(pps%gdd820(beg:end))
    pps%gdd820(:)= spval
    allocate(pps%gdd1020(beg:end))
    pps%gdd1020(:)= spval
    allocate(pps%harvdate(beg:end))
    pps%harvdate(:)= huge(1) 
    allocate(pps%htmx(beg:end))
    pps%htmx(:)= 0.0_r8
    allocate(pps%hdidx(beg:end))
    pps%hdidx(:)= nanr
    allocate(pps%cumvd(beg:end))
    pps%cumvd(:)= nanr
    allocate(pps%vf(beg:end))
    pps%vf(:)= 0.0_r8
    allocate(pps%cropplant(beg:end))
    pps%cropplant(:)= .false.
    allocate(pps%idop(beg:end))
    pps%idop(:)= huge(1)
    allocate(pps%vds(beg:end))
    pps%vds(:)= nanr
    allocate(pps%forc_hgt_u_pft(beg:end))
    pps%forc_hgt_u_pft(:)= nanr
    allocate(pps%forc_hgt_t_pft(beg:end))
    pps%forc_hgt_t_pft(:)= nanr
    allocate(pps%forc_hgt_q_pft(beg:end))
    pps%forc_hgt_q_pft(:)= nanr
    allocate(pps%lfpftd(beg:end))     

    allocate(pps%alphapsnsun(beg:end))
    pps%alphapsnsun(:)= spval
    allocate(pps%alphapsnsha(beg:end))
    pps%alphapsnsha(:)= spval

    allocate(pps%sandfrac(beg:end))
    pps%sandfrac(:)= nanr
    allocate(pps%clayfrac(beg:end))
    pps%clayfrac(:)= nanr
    allocate(pps%mlaidiff(beg:end))
    pps%mlaidiff(:)= nanr
    allocate(pps%rb1(beg:end))
    pps%rb1(:)= nanr
    allocate(pps%annlai(12,beg:end))
    pps%annlai(:,:)= nanr

    allocate(pps%irrig_rate(beg:end))
    pps%irrig_rate(:)= nanr
    allocate(pps%n_irrig_steps_left(beg:end))
    pps%n_irrig_steps_left(:)= 0

    allocate(pps%grnd_ch4_cond(beg:end))
    pps%grnd_ch4_cond(:)= nanr
    allocate(pps%canopy_cond(beg:end))
    pps%canopy_cond(:)= nanr

    allocate(pps%leaf_prof(beg:end,1:nlevdecomp_full))
    pps%leaf_prof(:,:)= spval
    allocate(pps%froot_prof(beg:end,1:nlevdecomp_full))
    pps%froot_prof(:,:)= spval
    allocate(pps%croot_prof(beg:end,1:nlevdecomp_full))
    pps%croot_prof(:,:)= spval
    allocate(pps%stem_prof(beg:end,1:nlevdecomp_full))
    pps%stem_prof(:,:)= spval

  end subroutine init_pft_pstate_type

  !------------------------------------------------------------------------
  subroutine init_pft_epv_type(beg, end, pepv)
    !
    ! !DESCRIPTION:
    ! Initialize pft ecophysiological variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_epv_type), intent(inout):: pepv
    !------------------------------------------------------------------------

    allocate(pepv%dormant_flag(beg:end))
    pepv%dormant_flag(:)= nanr
    allocate(pepv%days_active(beg:end))
    pepv%days_active(:)= nanr
    allocate(pepv%onset_flag(beg:end))
    pepv%onset_flag(:)= nanr
    allocate(pepv%onset_counter(beg:end))
    pepv%onset_counter(:)= nanr
    allocate(pepv%onset_gddflag(beg:end))
    pepv%onset_gddflag(:)= nanr
    allocate(pepv%onset_fdd(beg:end))
    pepv%onset_fdd(:)= nanr
    allocate(pepv%onset_gdd(beg:end))
    pepv%onset_gdd(:)= nanr
    allocate(pepv%onset_swi(beg:end))
    pepv%onset_swi(:)= nanr
    allocate(pepv%offset_flag(beg:end))
    pepv%offset_flag(:)= nanr
    allocate(pepv%offset_counter(beg:end))
    pepv%offset_counter(:)= nanr
    allocate(pepv%offset_fdd(beg:end))
    pepv%offset_fdd(:)= nanr
    allocate(pepv%offset_swi(beg:end))
    pepv%offset_swi(:)= nanr
    allocate(pepv%fert_counter(beg:end))
    pepv%fert_counter(:)= nanr
    allocate(pepv%grain_flag(beg:end))
    pepv%grain_flag(:)= nanr
    allocate(pepv%lgsf(beg:end))
    pepv%lgsf(:)= nanr
    allocate(pepv%bglfr(beg:end))
    pepv%bglfr(:)= nanr
    allocate(pepv%bgtr(beg:end))
    pepv%bgtr(:)= nanr
    allocate(pepv%annavg_t2m(beg:end))
    pepv%annavg_t2m(:)= nanr
    allocate(pepv%tempavg_t2m(beg:end))
    pepv%tempavg_t2m(:)= nanr
    allocate(pepv%gpp(beg:end))
    pepv%gpp(:)= nanr
    allocate(pepv%availc(beg:end))
    pepv%availc(:)= nanr
    allocate(pepv%xsmrpool_recover(beg:end))
    pepv%xsmrpool_recover(:)= nanr
    allocate(pepv%xsmrpool_c13ratio(beg:end))
    pepv%xsmrpool_c13ratio(:)= nanr
    allocate(pepv%alloc_pnow(beg:end))
    pepv%alloc_pnow(:)= nanr
    allocate(pepv%c_allometry(beg:end))
    pepv%c_allometry(:)= nanr
    allocate(pepv%n_allometry(beg:end))
    pepv%n_allometry(:)= nanr
    allocate(pepv%plant_ndemand(beg:end))
    pepv%plant_ndemand(:)= nanr
    allocate(pepv%tempsum_potential_gpp(beg:end))
    pepv%tempsum_potential_gpp(:)= nanr
    allocate(pepv%annsum_potential_gpp(beg:end))
    pepv%annsum_potential_gpp(:)= nanr
    allocate(pepv%tempmax_retransn(beg:end))
    pepv%tempmax_retransn(:)= nanr
    allocate(pepv%annmax_retransn(beg:end))
    pepv%annmax_retransn(:)= nanr
    allocate(pepv%avail_retransn(beg:end))
    pepv%avail_retransn(:)= nanr
    allocate(pepv%plant_nalloc(beg:end))
    pepv%plant_nalloc(:)= nanr
    allocate(pepv%plant_calloc(beg:end))
    pepv%plant_calloc(:)= nanr
    allocate(pepv%excess_cflux(beg:end))
    pepv%excess_cflux(:)= nanr
    allocate(pepv%downreg(beg:end))
    pepv%downreg(:)= nanr
    allocate(pepv%prev_leafc_to_litter(beg:end))
    pepv%prev_leafc_to_litter(:)= nanr
    allocate(pepv%prev_frootc_to_litter(beg:end))
    pepv%prev_frootc_to_litter(:)= nanr
    allocate(pepv%tempsum_npp(beg:end))
    pepv%tempsum_npp(:)= nanr
    allocate(pepv%annsum_npp(beg:end))
    pepv%annsum_npp(:)= nanr
    allocate(pepv%tempsum_litfall(beg:end))
    pepv%tempsum_litfall(:)= nanr
    allocate(pepv%annsum_litfall(beg:end))
    pepv%annsum_litfall(:)= nanr
    allocate(pepv%rc13_canair(beg:end))
    pepv%rc13_canair(:)= spval
    allocate(pepv%rc13_psnsun(beg:end))
    pepv%rc13_psnsun(:)= spval
    allocate(pepv%rc13_psnsha(beg:end))
    pepv%rc13_psnsha(:)= spval
    allocate(pepv%rc14_atm(beg:end))
    pepv%rc14_atm(:)= nanr

  end subroutine init_pft_epv_type

  !------------------------------------------------------------------------
  subroutine init_pft_pdgvstate_type(beg, end, pdgvs)
    !
    ! !DESCRIPTION:
    ! Initialize pft DGVM state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dgvstate_type), intent(inout):: pdgvs
    !------------------------------------------------------------------------

    allocate(pdgvs%agddtw(beg:end))
    pdgvs%agddtw(:)= nanr
    allocate(pdgvs%agdd(beg:end))
    pdgvs%agdd(:)= nanr
    allocate(pdgvs%t_mo(beg:end))
    pdgvs%t_mo(:)= nanr
    allocate(pdgvs%t_mo_min(beg:end))
    pdgvs%t_mo_min(:)= nanr
    allocate(pdgvs%agdd20(beg:end))
    pdgvs%agdd20(:)= nanr
    allocate(pdgvs%tmomin20(beg:end))
    pdgvs%tmomin20(:)= nanr
    allocate(pdgvs%prec365(beg:end))
    pdgvs%prec365(:)= nanr
    allocate(pdgvs%present(beg:end))
    pdgvs%present(:)= .false.
    allocate(pdgvs%pftmayexist(beg:end))
    pdgvs%pftmayexist(:)= .true.
    allocate(pdgvs%nind(beg:end))
    pdgvs%nind(:)= nanr
    allocate(pdgvs%lm_ind(beg:end))
    pdgvs%lm_ind(:)= nanr
    allocate(pdgvs%lai_ind(beg:end))
    pdgvs%lai_ind(:)= nanr
    allocate(pdgvs%fpcinc(beg:end))
    pdgvs%fpcinc(:)= nanr
    allocate(pdgvs%fpcgrid(beg:end))
    pdgvs%fpcgrid(:)= nanr
    allocate(pdgvs%fpcgridold(beg:end))
    pdgvs%fpcgridold(:)= nanr
    allocate(pdgvs%crownarea(beg:end))
    pdgvs%crownarea(:)= nanr
    allocate(pdgvs%greffic(beg:end))
    pdgvs%greffic(:)= nanr
    allocate(pdgvs%heatstress(beg:end))
    pdgvs%heatstress(:)= nanr

  end subroutine init_pft_pdgvstate_type

 !------------------------------------------------------------------------
  subroutine init_pft_vstate_type(beg, end, pvs)
    !
    ! !DESCRIPTION:
    ! Initialize pft VOC variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vstate_type), intent(inout):: pvs
    !------------------------------------------------------------------------

    allocate(pvs%t_veg24 (beg:end))
    pvs%t_veg24 (:)=spval
    allocate(pvs%t_veg240(beg:end))
    pvs%t_veg240(:)=spval
    allocate(pvs%fsd24   (beg:end))
    pvs%fsd24   (:)=spval
    allocate(pvs%fsd240  (beg:end))
    pvs%fsd240  (:)=spval
    allocate(pvs%fsi24   (beg:end))
    pvs%fsi24   (:)=spval
    allocate(pvs%fsi240  (beg:end))
    pvs%fsi240  (:)=spval
    allocate(pvs%fsun24  (beg:end))
    pvs%fsun24  (:)=spval
    allocate(pvs%fsun240 (beg:end))
    pvs%fsun240 (:)=spval

    allocate(pvs%elai_p  (beg:end))
    pvs%elai_p  (:)=spval

  end subroutine init_pft_vstate_type

  !------------------------------------------------------------------------
  subroutine init_pft_psynstate_type(beg, end, ppsyns)
    !
    ! !DESCRIPTION:
    ! Initialize pft energy state
    !
    ! !AGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_psynstate_type), intent(inout):: ppsyns
    !-----------------------------------------------------------------------

    allocate(ppsyns%c3flag(beg:end))
    ppsyns%c3flag(:)=.false.
    allocate(ppsyns%ac(beg:end,1:nlevcan))
    ppsyns%ac(:,:)= nanr
    allocate(ppsyns%aj(beg:end,1:nlevcan))
    ppsyns%aj(:,:)= nanr
    allocate(ppsyns%ap(beg:end,1:nlevcan))
    ppsyns%ap(:,:)= nanr
    allocate(ppsyns%ag(beg:end,1:nlevcan))
    ppsyns%ag(:,:)= nanr
    allocate(ppsyns%an(beg:end,1:nlevcan))
    ppsyns%an(:,:)= nanr
    allocate(ppsyns%vcmax_z(beg:end,1:nlevcan))
    ppsyns%vcmax_z(:,:)= nanr
    allocate(ppsyns%cp(beg:end))
    ppsyns%cp(:)= nanr
    allocate(ppsyns%kc(beg:end))
    ppsyns%kc(:)= nanr
    allocate(ppsyns%ko(beg:end))
    ppsyns%ko(:)= nanr
    allocate(ppsyns%qe(beg:end))
    ppsyns%qe(:)= nanr
    allocate(ppsyns%tpu_z(beg:end,1:nlevcan))
    ppsyns%tpu_z(:,:)= nanr
    allocate(ppsyns%kp_z(beg:end,1:nlevcan))
    ppsyns%kp_z(:,:)= nanr
    allocate(ppsyns%theta_cj(beg:end))
    ppsyns%theta_cj(:)= nanr
    allocate(ppsyns%bbb(beg:end))
    ppsyns%bbb(:)= nanr
    allocate(ppsyns%mbb(beg:end))
    ppsyns%mbb(:)= nanr
    allocate(ppsyns%gb_mol(beg:end))
    ppsyns%gb_mol(:)= nanr
    allocate(ppsyns%gs_mol(beg:end,1:nlevcan))
    ppsyns%gs_mol(:,:)= nanr

  end subroutine init_pft_psynstate_type

!------------------------------------------------------------------------
  subroutine init_pft_estate_type(beg, end, pes)
    !
    ! !DESCRIPTION:
    ! Initialize pft energy state
    !
    ! !AGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
    !-----------------------------------------------------------------------

    allocate(pes%t_ref2m(beg:end))
    pes%t_ref2m(:)= nanr
    allocate(pes%t_ref2m_min(beg:end))
    pes%t_ref2m_min(:)= nanr
    allocate(pes%t_ref2m_max(beg:end))
    pes%t_ref2m_max(:)= nanr
    allocate(pes%t_ref2m_min_inst(beg:end))
    pes%t_ref2m_min_inst(:)= nanr
    allocate(pes%t_ref2m_max_inst(beg:end))
    pes%t_ref2m_max_inst(:)= nanr
    allocate(pes%q_ref2m(beg:end))
    pes%q_ref2m(:)= nanr
    allocate(pes%t_ref2m_u(beg:end))
    pes%t_ref2m_u(:)= nanr
    allocate(pes%t_ref2m_r(beg:end))
    pes%t_ref2m_r(:)= nanr
    allocate(pes%t_ref2m_min_u(beg:end))
    pes%t_ref2m_min_u(:)= nanr
    allocate(pes%t_ref2m_min_r(beg:end))
    pes%t_ref2m_min_r(:)= nanr
    allocate(pes%t_ref2m_max_u(beg:end))
    pes%t_ref2m_max_u(:)= nanr
    allocate(pes%t_ref2m_max_r(beg:end))
    pes%t_ref2m_max_r(:)= nanr
    allocate(pes%t_ref2m_min_inst_u(beg:end))
    pes%t_ref2m_min_inst_u(:)= nanr
    allocate(pes%t_ref2m_min_inst_r(beg:end))
    pes%t_ref2m_min_inst_r(:)= nanr
    allocate(pes%t_ref2m_max_inst_u(beg:end))
    pes%t_ref2m_max_inst_u(:)= nanr
    allocate(pes%t_ref2m_max_inst_r(beg:end))
    pes%t_ref2m_max_inst_r(:)= nanr
    allocate(pes%t10(beg:end))
    pes%t10(:)= spval
    if ( crop_prog )then
       allocate(pes%a10tmin(beg:end))
       pes%a10tmin(:)= spval
       allocate(pes%a5tmin(beg:end))
       pes%a5tmin(:)= spval
    end if
    allocate(pes%rh_ref2m(beg:end))
    pes%rh_ref2m(:)= nanr
    allocate(pes%rh_ref2m_u(beg:end))
    pes%rh_ref2m_u(:)= nanr
    allocate(pes%rh_ref2m_r(beg:end))
    pes%rh_ref2m_r(:)= nanr
    allocate(pes%t_veg(beg:end))
    pes%t_veg(:)= nanr
    allocate(pes%thm(beg:end))
    pes%thm(:)= nanr

  end subroutine init_pft_estate_type

  !------------------------------------------------------------------------
  subroutine init_pft_wstate_type(beg, end, pws)
    !
    ! !DESCRIPTION:
    ! Initialize pft water state
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wstate_type), intent(inout):: pws !pft water state
    !------------------------------------------------------------------------

    allocate(pws%h2ocan(beg:end))
    pws%h2ocan(:)= spval

  end subroutine init_pft_wstate_type

  !------------------------------------------------------------------------
  subroutine init_pft_cstate_type(beg, end, pcs)
    !
    ! !DESCRIPTION:
    ! Initialize pft carbon state
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cstate_type), intent(inout):: pcs !pft carbon state
    !------------------------------------------------------------------------

    allocate(pcs%leafc(beg:end))
    pcs%leafc(:)= nanr
    allocate(pcs%leafc_storage(beg:end))
    pcs%leafc_storage(:)= nanr
    allocate(pcs%leafc_xfer(beg:end))
    pcs%leafc_xfer(:)= nanr
    allocate(pcs%frootc(beg:end))
    pcs%frootc(:)= nanr
    allocate(pcs%frootc_storage(beg:end))
    pcs%frootc_storage(:)= nanr
    allocate(pcs%frootc_xfer(beg:end))
    pcs%frootc_xfer(:)= nanr
    allocate(pcs%livestemc(beg:end))
    pcs%livestemc(:)= nanr
    allocate(pcs%livestemc_storage(beg:end))
    pcs%livestemc_storage(:)= nanr
    allocate(pcs%livestemc_xfer(beg:end))
    pcs%livestemc_xfer(:)= nanr
    allocate(pcs%deadstemc(beg:end))
    pcs%deadstemc(:)= nanr
    allocate(pcs%deadstemc_storage(beg:end))
    pcs%deadstemc_storage(:)= nanr
    allocate(pcs%deadstemc_xfer(beg:end))
    pcs%deadstemc_xfer(:)= nanr
    allocate(pcs%livecrootc(beg:end))
    pcs%livecrootc(:)= nanr
    allocate(pcs%livecrootc_storage(beg:end))
    pcs%livecrootc_storage(:)= nanr
    allocate(pcs%livecrootc_xfer(beg:end))
    pcs%livecrootc_xfer(:)= nanr
    allocate(pcs%deadcrootc(beg:end))
    pcs%deadcrootc(:)= nanr
    allocate(pcs%deadcrootc_storage(beg:end))
    pcs%deadcrootc_storage(:)= nanr
    allocate(pcs%deadcrootc_xfer(beg:end))
    pcs%deadcrootc_xfer(:)= nanr
    allocate(pcs%gresp_storage(beg:end))
    pcs%gresp_storage(:)= nanr
    allocate(pcs%gresp_xfer(beg:end))
    pcs%gresp_xfer(:)= nanr
    allocate(pcs%cpool(beg:end))
    pcs%cpool(:)= nanr
    allocate(pcs%xsmrpool(beg:end))
    pcs%xsmrpool(:)= nanr
    allocate(pcs%pft_ctrunc(beg:end))
    pcs%pft_ctrunc(:)= nanr
    allocate(pcs%dispvegc(beg:end))
    pcs%dispvegc(:)= nanr
    allocate(pcs%storvegc(beg:end))
    pcs%storvegc(:)= nanr
    allocate(pcs%totvegc(beg:end))
    pcs%totvegc(:)= nanr
    allocate(pcs%totpftc(beg:end))
    pcs%totpftc(:)= nanr
    allocate(pcs%leafcmax(beg:end))
    pcs%leafcmax(:)= nanr
    allocate(pcs%grainc(beg:end))
    pcs%grainc(:)= nanr
    allocate(pcs%grainc_storage(beg:end))
    pcs%grainc_storage(:)= nanr
    allocate(pcs%grainc_xfer(beg:end))
    pcs%grainc_xfer(:)= nanr
    allocate(pcs%woodc(beg:end))
    pcs%woodc(:)= nanr     

  end subroutine init_pft_cstate_type

  !------------------------------------------------------------------------
  subroutine init_pft_nstate_type(beg, end, pns)
    !
    ! !DESCRIPTION:
    ! Initialize pft nitrogen state
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nstate_type), intent(inout):: pns !pft nitrogen state
    !------------------------------------------------------------------------

    allocate(pns%grainn(beg:end))
    pns%grainn(:)= nanr
    allocate(pns%grainn_storage(beg:end))
    pns%grainn_storage(:)= nanr
    allocate(pns%grainn_xfer(beg:end))
    pns%grainn_xfer(:)= nanr
    allocate(pns%leafn(beg:end))
    pns%leafn(:)= nanr
    allocate(pns%leafn_storage(beg:end))
    pns%leafn_storage(:)= nanr
    allocate(pns%leafn_xfer(beg:end))
    pns%leafn_xfer(:)= nanr
    allocate(pns%frootn(beg:end))
    pns%frootn(:)= nanr
    allocate(pns%frootn_storage(beg:end))
    pns%frootn_storage(:)= nanr
    allocate(pns%frootn_xfer(beg:end))
    pns%frootn_xfer(:)= nanr
    allocate(pns%livestemn(beg:end))
    pns%livestemn(:)= nanr
    allocate(pns%livestemn_storage(beg:end))
    pns%livestemn_storage(:)= nanr
    allocate(pns%livestemn_xfer(beg:end))
    pns%livestemn_xfer(:)= nanr
    allocate(pns%deadstemn(beg:end))
    pns%deadstemn(:)= nanr
    allocate(pns%deadstemn_storage(beg:end))
    pns%deadstemn_storage(:)= nanr
    allocate(pns%deadstemn_xfer(beg:end))
    pns%deadstemn_xfer(:)= nanr
    allocate(pns%livecrootn(beg:end))
    pns%livecrootn(:)= nanr
    allocate(pns%livecrootn_storage(beg:end))
    pns%livecrootn_storage(:)= nanr
    allocate(pns%livecrootn_xfer(beg:end))
    pns%livecrootn_xfer(:)= nanr
    allocate(pns%deadcrootn(beg:end))
    pns%deadcrootn(:)= nanr
    allocate(pns%deadcrootn_storage(beg:end))
    pns%deadcrootn_storage(:)= nanr
    allocate(pns%deadcrootn_xfer(beg:end))
    pns%deadcrootn_xfer(:)= nanr
    allocate(pns%retransn(beg:end))
    pns%retransn(:)= nanr
    allocate(pns%npool(beg:end))
    pns%npool(:)= nanr
    allocate(pns%pft_ntrunc(beg:end))
    pns%pft_ntrunc(:)= nanr
    allocate(pns%dispvegn(beg:end))
    pns%dispvegn(:)= nanr
    allocate(pns%storvegn(beg:end))
    pns%storvegn(:)= nanr
    allocate(pns%totvegn(beg:end))
    pns%totvegn(:)= nanr
    allocate(pns%totpftn(beg:end))
    pns%totpftn(:)= nanr

  end subroutine init_pft_nstate_type

  !------------------------------------------------------------------------
  subroutine init_pft_eflux_type(beg, end, pef)
    !
    ! !DESCRIPTION:
    ! Initialize pft energy flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_eflux_type), intent(inout):: pef
    !------------------------------------------------------------------------

    allocate(pef%sabg(beg:end))
    pef%sabg(:)= nanr
    allocate(pef%sabv(beg:end))
    pef%sabv(:)= nanr
    allocate(pef%fsa(beg:end))
    pef%fsa(:)= nanr
    allocate(pef%fsa_u(beg:end))
    pef%fsa_u(:)= nanr
    allocate(pef%fsa_r(beg:end))
    pef%fsa_r(:)= nanr
    allocate(pef%fsr(beg:end))
    pef%fsr(:)= nanr
    allocate(pef%parsun_z(beg:end,1:nlevcan))
    pef%parsun_z(:,:)= nanr
    allocate(pef%parsha_z(beg:end,1:nlevcan))
    pef%parsha_z(:,:)= nanr
    allocate(pef%dlrad(beg:end))
    pef%dlrad(:)= nanr
    allocate(pef%ulrad(beg:end))
    pef%ulrad(:)= nanr
    allocate(pef%eflx_lh_tot(beg:end))
    pef%eflx_lh_tot(:)= nanr
    allocate(pef%eflx_lh_tot_u(beg:end))
    pef%eflx_lh_tot_u(:)= nanr
    allocate(pef%eflx_lh_tot_r(beg:end))
    pef%eflx_lh_tot_r(:)= nanr
    allocate(pef%eflx_lh_grnd(beg:end))
    pef%eflx_lh_grnd(:)= nanr
    allocate(pef%eflx_soil_grnd(beg:end))
    pef%eflx_soil_grnd(:)= nanr
    allocate(pef%eflx_soil_grnd_u(beg:end))
    pef%eflx_soil_grnd_u(:)= nanr
    allocate(pef%eflx_soil_grnd_r(beg:end))
    pef%eflx_soil_grnd_r(:)= nanr
    allocate(pef%eflx_sh_tot(beg:end))
    pef%eflx_sh_tot(:)= nanr
    allocate(pef%eflx_sh_tot_u(beg:end))
    pef%eflx_sh_tot_u(:)= nanr
    allocate(pef%eflx_sh_tot_r(beg:end))
    pef%eflx_sh_tot_r(:)= nanr
    allocate(pef%eflx_sh_grnd(beg:end))
    pef%eflx_sh_grnd(:)= nanr
    allocate(pef%eflx_sh_veg(beg:end))
    pef%eflx_sh_veg(:)= nanr
    allocate(pef%eflx_lh_vege(beg:end))
    pef%eflx_lh_vege(:)= nanr
    allocate(pef%eflx_lh_vegt(beg:end))
    pef%eflx_lh_vegt(:)= nanr
    allocate(pef%eflx_wasteheat_pft(beg:end))
    pef%eflx_wasteheat_pft(:)= nanr
    allocate(pef%eflx_heat_from_ac_pft(beg:end))
    pef%eflx_heat_from_ac_pft(:)= nanr
    allocate(pef%eflx_traffic_pft(beg:end))
    pef%eflx_traffic_pft(:)= nanr
    allocate(pef%eflx_anthro(beg:end))
    pef%eflx_anthro(:)= nanr
    allocate(pef%cgrnd(beg:end))
    pef%cgrnd(:)= nanr
    allocate(pef%cgrndl(beg:end))
    pef%cgrndl(:)= nanr
    allocate(pef%cgrnds(beg:end))
    pef%cgrnds(:)= nanr
    allocate(pef%eflx_gnet(beg:end))
    pef%eflx_gnet(:)= nanr
    allocate(pef%eflx_grnd_lake(beg:end))
    pef%eflx_grnd_lake(:)= nanr
    allocate(pef%dgnetdT(beg:end))
    pef%dgnetdT(:)= nanr
    allocate(pef%eflx_lwrad_out(beg:end))
    pef%eflx_lwrad_out(:)= nanr
    allocate(pef%eflx_lwrad_net(beg:end))
    pef%eflx_lwrad_net(:)= nanr
    allocate(pef%eflx_lwrad_net_u(beg:end))
    pef%eflx_lwrad_net_u(:)= nanr
    allocate(pef%eflx_lwrad_net_r(beg:end))
    pef%eflx_lwrad_net_r(:)= nanr
    allocate(pef%eflx_lwrad_out_u(beg:end))
    pef%eflx_lwrad_out_u(beg:end) = nanr
    allocate(pef%eflx_lwrad_out_r(beg:end))
    pef%eflx_lwrad_out_r(beg:end) = nanr
    allocate(pef%netrad(beg:end))
    pef%netrad(:)= nanr
    allocate(pef%fsds_vis_d(beg:end))
    pef%fsds_vis_d(:)= nanr
    allocate(pef%fsds_nir_d(beg:end))
    pef%fsds_nir_d(:)= nanr
    allocate(pef%fsds_vis_i(beg:end))
    pef%fsds_vis_i(:)= nanr
    allocate(pef%fsds_nir_i(beg:end))
    pef%fsds_nir_i(:)= nanr
    allocate(pef%fsr_vis_d(beg:end))
    pef%fsr_vis_d(:)= nanr
    allocate(pef%fsr_nir_d(beg:end))
    pef%fsr_nir_d(:)= nanr
    allocate(pef%fsr_vis_i(beg:end))
    pef%fsr_vis_i(:)= nanr
    allocate(pef%fsr_nir_i(beg:end))
    pef%fsr_nir_i(:)= nanr
    allocate(pef%fsds_vis_d_ln(beg:end))
    pef%fsds_vis_d_ln(:)= nanr
    allocate(pef%fsds_vis_i_ln(beg:end))
    pef%fsds_vis_i_ln(:)= nanr
    allocate(pef%parveg_ln(beg:end))
    pef%parveg_ln(:)= nanr
    allocate(pef%fsds_nir_d_ln(beg:end))
    pef%fsds_nir_d_ln(:)= nanr
    allocate(pef%fsr_vis_d_ln(beg:end))
    pef%fsr_vis_d_ln(:)= nanr
    allocate(pef%fsr_nir_d_ln(beg:end))
    pef%fsr_nir_d_ln(:)= nanr
    allocate(pef%sabg_lyr(beg:end,-nlevsno+1:1))
    pef%sabg_lyr(:,:)= spval
    allocate(pef%sabg_pen(beg:end))
    pef%sabg_pen(:)= nanr
    allocate(pef%sfc_frc_aer(beg:end))
    pef%sfc_frc_aer(:)= nanr
    allocate(pef%sfc_frc_bc(beg:end))
    pef%sfc_frc_bc(:)= nanr
    allocate(pef%sfc_frc_oc(beg:end))
    pef%sfc_frc_oc(:)= nanr
    allocate(pef%sfc_frc_dst(beg:end))
    pef%sfc_frc_dst(:)= nanr
    allocate(pef%sfc_frc_aer_sno(beg:end))
    pef%sfc_frc_aer_sno(:)= nanr
    allocate(pef%sfc_frc_bc_sno(beg:end))
    pef%sfc_frc_bc_sno(:)= nanr
    allocate(pef%sfc_frc_oc_sno(beg:end))
    pef%sfc_frc_oc_sno(:)= nanr
    allocate(pef%sfc_frc_dst_sno(beg:end))
    pef%sfc_frc_dst_sno(:)= nanr
    allocate(pef%fsr_sno_vd(beg:end))
    pef%fsr_sno_vd(:)= nanr
    allocate(pef%fsr_sno_nd(beg:end))
    pef%fsr_sno_nd(:)= nanr
    allocate(pef%fsr_sno_vi(beg:end))
    pef%fsr_sno_vi(:)= nanr
    allocate(pef%fsr_sno_ni(beg:end))
    pef%fsr_sno_ni(:)= nanr
    allocate(pef%fsds_sno_vd(beg:end))
    pef%fsds_sno_vd(:)= nanr
    allocate(pef%fsds_sno_nd(beg:end))
    pef%fsds_sno_nd(:)= nanr
    allocate(pef%fsds_sno_vi(beg:end))
    pef%fsds_sno_vi(:)= nanr
    allocate(pef%fsds_sno_ni(beg:end))
    pef%fsds_sno_ni(:)= nanr

    allocate(pef%eflx_sh_snow(beg:end))
    pef%eflx_sh_snow(:)= nanr
    allocate(pef%eflx_sh_soil(beg:end))
    pef%eflx_sh_soil(:)= nanr
    allocate(pef%eflx_sh_h2osfc(beg:end))
    pef%eflx_sh_h2osfc(:)= nanr

    allocate(pef%sabg_soil(beg:end))
    pef%sabg_soil(:)= nanr
    allocate(pef%sabg_snow(beg:end))
    pef%sabg_snow(:)= nanr
    allocate(pef%sabg_chk(beg:end))
    pef%sabg_chk(:)= nanr
  end subroutine init_pft_eflux_type

  !------------------------------------------------------------------------
  subroutine init_pft_mflux_type(beg, end, pmf)
    !
    ! !DESCRIPTION:
    ! Initialize pft momentum flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_mflux_type), intent(inout) :: pmf
    !------------------------------------------------------------------------

    allocate(pmf%taux(beg:end))
    pmf%taux(:)= nanr
    allocate(pmf%tauy(beg:end))
    pmf%tauy(:)= nanr

  end subroutine init_pft_mflux_type

  !------------------------------------------------------------------------
  subroutine init_pft_wflux_type(beg, end, pwf)
    !
    ! !DESCRIPTION:
    ! Initialize pft water flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wflux_type), intent(inout) :: pwf
    !------------------------------------------------------------------------

    allocate(pwf%qflx_prec_intr(beg:end))
    pwf%qflx_prec_intr(:)= nanr
    allocate(pwf%qflx_prec_grnd(beg:end))
    pwf%qflx_prec_grnd(:)= nanr
    allocate(pwf%qflx_rain_grnd(beg:end))
    pwf%qflx_rain_grnd(:)= nanr
    allocate(pwf%qflx_snow_grnd(beg:end))
    pwf%qflx_snow_grnd(:)= nanr
    allocate(pwf%qflx_snwcp_liq(beg:end))
    pwf%qflx_snwcp_liq(:)= nanr
    allocate(pwf%qflx_snwcp_ice(beg:end))
    pwf%qflx_snwcp_ice(:)= nanr
    allocate(pwf%qflx_evap_veg(beg:end))
    pwf%qflx_evap_veg(:)= nanr
    allocate(pwf%qflx_tran_veg(beg:end))
    pwf%qflx_tran_veg(:)= nanr
    allocate(pwf%qflx_evap_can(beg:end))
    pwf%qflx_evap_can(:)= nanr
    allocate(pwf%qflx_evap_soi(beg:end))
    pwf%qflx_evap_soi(:)= nanr
    allocate(pwf%qflx_evap_tot(beg:end))
    pwf%qflx_evap_tot(:)= nanr
    allocate(pwf%qflx_evap_grnd(beg:end))
    pwf%qflx_evap_grnd(:)= 0.0_r8
    allocate(pwf%qflx_dew_grnd(beg:end))
    pwf%qflx_dew_grnd(:)= 0.0_r8
    allocate(pwf%qflx_sub_snow(beg:end))
    pwf%qflx_sub_snow(:)= 0.0_r8
    allocate(pwf%qflx_dew_snow(beg:end))
    pwf%qflx_dew_snow(:)= 0.0_r8
    allocate(pwf%qflx_ev_snow(beg:end))
    pwf%qflx_ev_snow(:)= nanr
    allocate(pwf%qflx_ev_soil(beg:end))
    pwf%qflx_ev_soil(:)= nanr
    allocate(pwf%qflx_ev_h2osfc(beg:end))
    pwf%qflx_ev_h2osfc(:)= nanr
    allocate(pwf%qflx_irrig(beg:end))
    pwf%qflx_irrig(:)= nanr

  end subroutine init_pft_wflux_type

  !------------------------------------------------------------------------
  subroutine init_pft_cflux_type(beg, end, pcf)
    !
    ! !DESCRIPTION:
    ! Initialize pft carbon flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cflux_type), intent(inout) :: pcf
    !------------------------------------------------------------------------

    allocate(pcf%psnsun(beg:end))
    pcf%psnsun(:)= nanr
    allocate(pcf%psnsha(beg:end))
    pcf%psnsha(:)= nanr
    allocate(pcf%psnsun_z(beg:end,1:nlevcan))
    pcf%psnsun_z(:,:)= nanr
    allocate(pcf%psnsha_z(beg:end,1:nlevcan))
    pcf%psnsha_z(:,:)= nanr
    allocate(pcf%cisun_z(beg:end,1:nlevcan))
    pcf%cisun_z(:,:)= nanr
    allocate(pcf%cisha_z(beg:end,1:nlevcan))
    pcf%cisha_z(:,:)= nanr
    allocate(pcf%lmrsun(beg:end))
    pcf%lmrsun(:)= nanr
    allocate(pcf%lmrsha(beg:end))
    pcf%lmrsha(:)= nanr
    allocate(pcf%lmrsun_z(beg:end,1:nlevcan))
    pcf%lmrsun_z(:,:)= nanr
    allocate(pcf%lmrsha_z(beg:end,1:nlevcan))
    pcf%lmrsha_z(:,:)= nanr
    allocate(pcf%fpsn(beg:end))
    pcf%fpsn(:)= spval
    allocate(pcf%fco2(beg:end))
    pcf%fco2(:)= 0._r8
    allocate(pcf%psnsun_wc(beg:end))
    pcf%psnsun_wc(:)= nanr
    allocate(pcf%psnsha_wc(beg:end))
    pcf%psnsha_wc(:)= nanr
    allocate(pcf%fpsn_wc(beg:end))
    pcf%fpsn_wc(:)= nanr
    allocate(pcf%psnsun_wj(beg:end))
    pcf%psnsun_wj(:)= nanr
    allocate(pcf%psnsha_wj(beg:end))
    pcf%psnsha_wj(:)= nanr
    allocate(pcf%fpsn_wj(beg:end))
    pcf%fpsn_wj(:)= nanr
    allocate(pcf%psnsun_wp(beg:end))
    pcf%psnsun_wp(:)= nanr
    allocate(pcf%psnsha_wp(beg:end))
    pcf%psnsha_wp(:)= nanr
    allocate(pcf%fpsn_wp(beg:end))
    pcf%fpsn_wp(:)= nanr

    allocate(pcf%m_leafc_to_litter(beg:end))
    pcf%m_leafc_to_litter(:)= nanr
    allocate(pcf%m_frootc_to_litter(beg:end))
    pcf%m_frootc_to_litter(:)= nanr
    allocate(pcf%m_leafc_storage_to_litter(beg:end))
    pcf%m_leafc_storage_to_litter(:)= nanr
    allocate(pcf%m_frootc_storage_to_litter(beg:end))
    pcf%m_frootc_storage_to_litter(:)= nanr
    allocate(pcf%m_livestemc_storage_to_litter(beg:end))
    pcf%m_livestemc_storage_to_litter(:)= nanr
    allocate(pcf%m_deadstemc_storage_to_litter(beg:end))
    pcf%m_deadstemc_storage_to_litter(:)= nanr
    allocate(pcf%m_livecrootc_storage_to_litter(beg:end))
    pcf%m_livecrootc_storage_to_litter(:)= nanr
    allocate(pcf%m_deadcrootc_storage_to_litter(beg:end))
    pcf%m_deadcrootc_storage_to_litter(:)= nanr
    allocate(pcf%m_leafc_xfer_to_litter(beg:end))
    pcf%m_leafc_xfer_to_litter(:)= nanr
    allocate(pcf%m_frootc_xfer_to_litter(beg:end))
    pcf%m_frootc_xfer_to_litter(:)= nanr
    allocate(pcf%m_livestemc_xfer_to_litter(beg:end))
    pcf%m_livestemc_xfer_to_litter(:)= nanr
    allocate(pcf%m_deadstemc_xfer_to_litter(beg:end))
    pcf%m_deadstemc_xfer_to_litter(:)= nanr
    allocate(pcf%m_livecrootc_xfer_to_litter(beg:end))
    pcf%m_livecrootc_xfer_to_litter(:)= nanr
    allocate(pcf%m_deadcrootc_xfer_to_litter(beg:end))
    pcf%m_deadcrootc_xfer_to_litter(:)= nanr
    allocate(pcf%m_livestemc_to_litter(beg:end))
    pcf%m_livestemc_to_litter(:)= nanr
    allocate(pcf%m_deadstemc_to_litter(beg:end))
    pcf%m_deadstemc_to_litter(:)= nanr
    allocate(pcf%m_livecrootc_to_litter(beg:end))
    pcf%m_livecrootc_to_litter(:)= nanr
    allocate(pcf%m_deadcrootc_to_litter(beg:end))
    pcf%m_deadcrootc_to_litter(:)= nanr
    allocate(pcf%m_gresp_storage_to_litter(beg:end))
    pcf%m_gresp_storage_to_litter(:)= nanr
    allocate(pcf%m_gresp_xfer_to_litter(beg:end))
    pcf%m_gresp_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_leafc_to_litter(beg:end))
    pcf%hrv_leafc_to_litter(:)= nanr
    allocate(pcf%hrv_leafc_storage_to_litter(beg:end))
    pcf%hrv_leafc_storage_to_litter(:)= nanr
    allocate(pcf%hrv_leafc_xfer_to_litter(beg:end))
    pcf%hrv_leafc_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_frootc_to_litter(beg:end))
    pcf%hrv_frootc_to_litter(:)= nanr
    allocate(pcf%hrv_frootc_storage_to_litter(beg:end))
    pcf%hrv_frootc_storage_to_litter(:)= nanr
    allocate(pcf%hrv_frootc_xfer_to_litter(beg:end))
    pcf%hrv_frootc_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_livestemc_to_litter(beg:end))
    pcf%hrv_livestemc_to_litter(:)= nanr
    allocate(pcf%hrv_livestemc_storage_to_litter(beg:end))
    pcf%hrv_livestemc_storage_to_litter(:)= nanr
    allocate(pcf%hrv_livestemc_xfer_to_litter(beg:end))
    pcf%hrv_livestemc_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_deadstemc_to_prod10c(beg:end))
    pcf%hrv_deadstemc_to_prod10c(:)= nanr
    allocate(pcf%hrv_deadstemc_to_prod100c(beg:end))
    pcf%hrv_deadstemc_to_prod100c(:)= nanr
    allocate(pcf%hrv_deadstemc_storage_to_litter(beg:end))
    pcf%hrv_deadstemc_storage_to_litter(:)= nanr
    allocate(pcf%hrv_deadstemc_xfer_to_litter(beg:end))
    pcf%hrv_deadstemc_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_livecrootc_to_litter(beg:end))
    pcf%hrv_livecrootc_to_litter(:)= nanr
    allocate(pcf%hrv_livecrootc_storage_to_litter(beg:end))
    pcf%hrv_livecrootc_storage_to_litter(:)= nanr
    allocate(pcf%hrv_livecrootc_xfer_to_litter(beg:end))
    pcf%hrv_livecrootc_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_deadcrootc_to_litter(beg:end))
    pcf%hrv_deadcrootc_to_litter(:)= nanr
    allocate(pcf%hrv_deadcrootc_storage_to_litter(beg:end))
    pcf%hrv_deadcrootc_storage_to_litter(:)= nanr
    allocate(pcf%hrv_deadcrootc_xfer_to_litter(beg:end))
    pcf%hrv_deadcrootc_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_gresp_storage_to_litter(beg:end))
    pcf%hrv_gresp_storage_to_litter(:)= nanr
    allocate(pcf%hrv_gresp_xfer_to_litter(beg:end))
    pcf%hrv_gresp_xfer_to_litter(:)= nanr
    allocate(pcf%hrv_xsmrpool_to_atm(beg:end))
    pcf%hrv_xsmrpool_to_atm(:)= nanr

    ! fire related variables changed by F. Li and S. Levis           
    allocate(pcf%m_leafc_to_fire(beg:end))
    pcf%m_leafc_to_fire(:)= nanr
    allocate(pcf%m_leafc_storage_to_fire(beg:end))
    pcf%m_leafc_storage_to_fire(:)= nanr
    allocate(pcf%m_leafc_xfer_to_fire(beg:end))
    pcf%m_leafc_xfer_to_fire(:)= nanr
    allocate(pcf%m_livestemc_to_fire(beg:end))
    pcf%m_livestemc_to_fire(:)= nanr
    allocate(pcf%m_livestemc_storage_to_fire(beg:end))
    pcf%m_livestemc_storage_to_fire(:)= nanr
    allocate(pcf%m_livestemc_xfer_to_fire(beg:end))
    pcf%m_livestemc_xfer_to_fire(:)= nanr
    allocate(pcf%m_deadstemc_to_fire(beg:end))
    pcf%m_deadstemc_to_fire(:)= nanr
    allocate(pcf%m_deadstemc_storage_to_fire(beg:end))
    pcf%m_deadstemc_storage_to_fire(:)= nanr
    allocate(pcf%m_deadstemc_xfer_to_fire(beg:end))
    pcf%m_deadstemc_xfer_to_fire(:)= nanr
    allocate(pcf%m_frootc_to_fire(beg:end))
    pcf%m_frootc_to_fire(:)= nanr
    allocate(pcf%m_frootc_storage_to_fire(beg:end))
    pcf%m_frootc_storage_to_fire(:)= nanr
    allocate(pcf%m_frootc_xfer_to_fire(beg:end))
    pcf%m_frootc_xfer_to_fire(:)= nanr
    allocate(pcf%m_livecrootc_to_fire(beg:end))
    pcf%m_livecrootc_to_fire(:)= nanr
    allocate(pcf%m_livecrootc_storage_to_fire(beg:end))
    pcf%m_livecrootc_storage_to_fire(:)= nanr
    allocate(pcf%m_livecrootc_xfer_to_fire(beg:end))
    pcf%m_livecrootc_xfer_to_fire(:)= nanr
    allocate(pcf%m_deadcrootc_to_fire(beg:end))
    pcf%m_deadcrootc_to_fire(:)= nanr
    allocate(pcf%m_deadcrootc_storage_to_fire(beg:end))
    pcf%m_deadcrootc_storage_to_fire(:)= nanr
    allocate(pcf%m_deadcrootc_xfer_to_fire(beg:end))
    pcf%m_deadcrootc_xfer_to_fire(:)= nanr
    allocate(pcf%m_gresp_storage_to_fire(beg:end))
    pcf%m_gresp_storage_to_fire(:)= nanr
    allocate(pcf%m_gresp_xfer_to_fire(beg:end))
    pcf%m_gresp_xfer_to_fire(:)= nanr
    allocate(pcf%m_leafc_to_litter_fire(beg:end))
    pcf%m_leafc_to_litter_fire(:)= nanr
    allocate(pcf%m_leafc_storage_to_litter_fire(beg:end))
    pcf%m_leafc_storage_to_litter_fire(:)= nanr
    allocate(pcf%m_leafc_xfer_to_litter_fire(beg:end))
    pcf%m_leafc_xfer_to_litter_fire(:)= nanr
    allocate(pcf%m_livestemc_to_litter_fire(beg:end))
    pcf%m_livestemc_to_litter_fire(:)= nanr
    allocate(pcf%m_livestemc_storage_to_litter_fire(beg:end))
    pcf%m_livestemc_storage_to_litter_fire(:)= nanr
    allocate(pcf%m_livestemc_xfer_to_litter_fire(beg:end))
    pcf%m_livestemc_xfer_to_litter_fire(:)= nanr
    allocate(pcf%m_livestemc_to_deadstemc_fire(beg:end))
    pcf%m_livestemc_to_deadstemc_fire(:)= nanr
    allocate(pcf%m_deadstemc_to_litter_fire(beg:end))
    pcf%m_deadstemc_to_litter_fire(:)= nanr
    allocate(pcf%m_deadstemc_storage_to_litter_fire(beg:end))
    pcf%m_deadstemc_storage_to_litter_fire(:)= nanr
    allocate(pcf%m_deadstemc_xfer_to_litter_fire(beg:end))
    pcf%m_deadstemc_xfer_to_litter_fire(:)= nanr
    allocate(pcf%m_frootc_to_litter_fire(beg:end))
    pcf%m_frootc_to_litter_fire(:)= nanr
    allocate(pcf%m_frootc_storage_to_litter_fire(beg:end))
    pcf%m_frootc_storage_to_litter_fire(:)= nanr
    allocate(pcf%m_frootc_xfer_to_litter_fire(beg:end))
    pcf%m_frootc_xfer_to_litter_fire(:)= nanr
    allocate(pcf%m_livecrootc_to_litter_fire(beg:end))
    pcf%m_livecrootc_to_litter_fire(:)= nanr
    allocate(pcf%m_livecrootc_storage_to_litter_fire(beg:end))
    pcf%m_livecrootc_storage_to_litter_fire(:)= nanr
    allocate(pcf%m_livecrootc_xfer_to_litter_fire(beg:end))
    pcf%m_livecrootc_xfer_to_litter_fire(:)= nanr
    allocate(pcf%m_livecrootc_to_deadcrootc_fire(beg:end))
    pcf%m_livecrootc_to_deadcrootc_fire(:)= nanr
    allocate(pcf%m_deadcrootc_to_litter_fire(beg:end))
    pcf%m_deadcrootc_to_litter_fire(:)= nanr
    allocate(pcf%m_deadcrootc_storage_to_litter_fire(beg:end))
    pcf%m_deadcrootc_storage_to_litter_fire(:)= nanr
    allocate(pcf%m_deadcrootc_xfer_to_litter_fire(beg:end))
    pcf%m_deadcrootc_xfer_to_litter_fire(:)= nanr
    allocate(pcf%m_gresp_storage_to_litter_fire(beg:end))
    pcf%m_gresp_storage_to_litter_fire(:)= nanr
    allocate(pcf%m_gresp_xfer_to_litter_fire(beg:end))
    pcf%m_gresp_xfer_to_litter_fire(:)= nanr


    allocate(pcf%leafc_xfer_to_leafc(beg:end))
    pcf%leafc_xfer_to_leafc(:)= nanr
    allocate(pcf%frootc_xfer_to_frootc(beg:end))
    pcf%frootc_xfer_to_frootc(:)= nanr
    allocate(pcf%livestemc_xfer_to_livestemc(beg:end))
    pcf%livestemc_xfer_to_livestemc(:)= nanr
    allocate(pcf%deadstemc_xfer_to_deadstemc(beg:end))
    pcf%deadstemc_xfer_to_deadstemc(:)= nanr
    allocate(pcf%livecrootc_xfer_to_livecrootc(beg:end))
    pcf%livecrootc_xfer_to_livecrootc(:)= nanr
    allocate(pcf%deadcrootc_xfer_to_deadcrootc(beg:end))
    pcf%deadcrootc_xfer_to_deadcrootc(:)= nanr
    allocate(pcf%leafc_to_litter(beg:end))
    pcf%leafc_to_litter(:)= nanr
    allocate(pcf%frootc_to_litter(beg:end))
    pcf%frootc_to_litter(:)= nanr
    allocate(pcf%leaf_mr(beg:end))
    pcf%leaf_mr(:)= nanr
    allocate(pcf%froot_mr(beg:end))
    pcf%froot_mr(:)= nanr
    allocate(pcf%livestem_mr(beg:end))
    pcf%livestem_mr(:)= nanr
    allocate(pcf%livecroot_mr(beg:end))
    pcf%livecroot_mr(:)= nanr
    allocate(pcf%grain_mr(beg:end))
    pcf%grain_mr(:)= nanr
    allocate(pcf%leaf_curmr(beg:end))
    pcf%leaf_curmr(:)= nanr
    allocate(pcf%froot_curmr(beg:end))
    pcf%froot_curmr(:)= nanr
    allocate(pcf%livestem_curmr(beg:end))
    pcf%livestem_curmr(:)= nanr
    allocate(pcf%livecroot_curmr(beg:end))
    pcf%livecroot_curmr(:)= nanr
    allocate(pcf%grain_curmr(beg:end))
    pcf%grain_curmr(:)= nanr
    allocate(pcf%leaf_xsmr(beg:end))
    pcf%leaf_xsmr(:)= nanr
    allocate(pcf%froot_xsmr(beg:end))
    pcf%froot_xsmr(:)= nanr
    allocate(pcf%livestem_xsmr(beg:end))
    pcf%livestem_xsmr(:)= nanr
    allocate(pcf%livecroot_xsmr(beg:end))
    pcf%livecroot_xsmr(:)= nanr
    allocate(pcf%grain_xsmr(beg:end))
    pcf%grain_xsmr(:)= nanr
    allocate(pcf%psnsun_to_cpool(beg:end))
    pcf%psnsun_to_cpool(:)= nanr
    allocate(pcf%psnshade_to_cpool(beg:end))
    pcf%psnshade_to_cpool(:)= nanr
    allocate(pcf%cpool_to_xsmrpool(beg:end))
    pcf%cpool_to_xsmrpool(:)= nanr
    allocate(pcf%cpool_to_leafc(beg:end))
    pcf%cpool_to_leafc(:)= nanr
    allocate(pcf%cpool_to_leafc_storage(beg:end))
    pcf%cpool_to_leafc_storage(:)= nanr
    allocate(pcf%cpool_to_frootc(beg:end))
    pcf%cpool_to_frootc(:)= nanr
    allocate(pcf%cpool_to_frootc_storage(beg:end))
    pcf%cpool_to_frootc_storage(:)= nanr
    allocate(pcf%cpool_to_livestemc(beg:end))
    pcf%cpool_to_livestemc(:)= nanr
    allocate(pcf%cpool_to_livestemc_storage(beg:end))
    pcf%cpool_to_livestemc_storage(:)= nanr
    allocate(pcf%cpool_to_deadstemc(beg:end))
    pcf%cpool_to_deadstemc(:)= nanr
    allocate(pcf%cpool_to_deadstemc_storage(beg:end))
    pcf%cpool_to_deadstemc_storage(:)= nanr
    allocate(pcf%cpool_to_livecrootc(beg:end))
    pcf%cpool_to_livecrootc(:)= nanr
    allocate(pcf%cpool_to_livecrootc_storage(beg:end))
    pcf%cpool_to_livecrootc_storage(:)= nanr
    allocate(pcf%cpool_to_deadcrootc(beg:end))
    pcf%cpool_to_deadcrootc(:)= nanr
    allocate(pcf%cpool_to_deadcrootc_storage(beg:end))
    pcf%cpool_to_deadcrootc_storage(:)= nanr
    allocate(pcf%cpool_to_gresp_storage(beg:end))
    pcf%cpool_to_gresp_storage(:)= nanr
    allocate(pcf%cpool_leaf_gr(beg:end))
    pcf%cpool_leaf_gr(:)= nanr
    allocate(pcf%cpool_leaf_storage_gr(beg:end))
    pcf%cpool_leaf_storage_gr(:)= nanr
    allocate(pcf%transfer_leaf_gr(beg:end))
    pcf%transfer_leaf_gr(:)= nanr
    allocate(pcf%cpool_froot_gr(beg:end))
    pcf%cpool_froot_gr(:)= nanr
    allocate(pcf%cpool_froot_storage_gr(beg:end))
    pcf%cpool_froot_storage_gr(:)= nanr
    allocate(pcf%transfer_froot_gr(beg:end))
    pcf%transfer_froot_gr(:)= nanr
    allocate(pcf%cpool_livestem_gr(beg:end))
    pcf%cpool_livestem_gr(:)= nanr
    allocate(pcf%cpool_livestem_storage_gr(beg:end))
    pcf%cpool_livestem_storage_gr(:)= nanr
    allocate(pcf%transfer_livestem_gr(beg:end))
    pcf%transfer_livestem_gr(:)= nanr
    allocate(pcf%cpool_deadstem_gr(beg:end))
    pcf%cpool_deadstem_gr(:)= nanr
    allocate(pcf%cpool_deadstem_storage_gr(beg:end))
    pcf%cpool_deadstem_storage_gr(:)= nanr
    allocate(pcf%transfer_deadstem_gr(beg:end))
    pcf%transfer_deadstem_gr(:)= nanr
    allocate(pcf%cpool_livecroot_gr(beg:end))
    pcf%cpool_livecroot_gr(:)= nanr
    allocate(pcf%cpool_livecroot_storage_gr(beg:end))
    pcf%cpool_livecroot_storage_gr(:)= nanr
    allocate(pcf%transfer_livecroot_gr(beg:end))
    pcf%transfer_livecroot_gr(:)= nanr
    allocate(pcf%cpool_deadcroot_gr(beg:end))
    pcf%cpool_deadcroot_gr(:)= nanr
    allocate(pcf%cpool_deadcroot_storage_gr(beg:end))
    pcf%cpool_deadcroot_storage_gr(:)= nanr
    allocate(pcf%transfer_deadcroot_gr(beg:end))
    pcf%transfer_deadcroot_gr(:)= nanr
    allocate(pcf%leafc_storage_to_xfer(beg:end))
    pcf%leafc_storage_to_xfer(:)= nanr
    allocate(pcf%frootc_storage_to_xfer(beg:end))
    pcf%frootc_storage_to_xfer(:)= nanr
    allocate(pcf%livestemc_storage_to_xfer(beg:end))
    pcf%livestemc_storage_to_xfer(:)= nanr
    allocate(pcf%deadstemc_storage_to_xfer(beg:end))
    pcf%deadstemc_storage_to_xfer(:)= nanr
    allocate(pcf%livecrootc_storage_to_xfer(beg:end))
    pcf%livecrootc_storage_to_xfer(:)= nanr
    allocate(pcf%deadcrootc_storage_to_xfer(beg:end))
    pcf%deadcrootc_storage_to_xfer(:)= nanr
    allocate(pcf%gresp_storage_to_xfer(beg:end))
    pcf%gresp_storage_to_xfer(:)= nanr
    allocate(pcf%livestemc_to_deadstemc(beg:end))
    pcf%livestemc_to_deadstemc(:)= nanr
    allocate(pcf%livecrootc_to_deadcrootc(beg:end))
    pcf%livecrootc_to_deadcrootc(:)= nanr
    allocate(pcf%gpp(beg:end))
    pcf%gpp(:)= nanr
    allocate(pcf%mr(beg:end))
    pcf%mr(:)= nanr
    allocate(pcf%current_gr(beg:end))
    pcf%current_gr(:)= nanr
    allocate(pcf%transfer_gr(beg:end))
    pcf%transfer_gr(:)= nanr
    allocate(pcf%storage_gr(beg:end))
    pcf%storage_gr(:)= nanr
    allocate(pcf%gr(beg:end))
    pcf%gr(:)= nanr
    allocate(pcf%ar(beg:end))
    pcf%ar(:)= nanr
    allocate(pcf%rr(beg:end))
    pcf%rr(:)= nanr
    allocate(pcf%npp(beg:end))
    pcf%npp(:)= nanr
    allocate(pcf%agnpp(beg:end))
    pcf%agnpp(:)= nanr
    allocate(pcf%bgnpp(beg:end))
    pcf%bgnpp(:)= nanr
    allocate(pcf%litfall(beg:end))
    pcf%litfall(:)= nanr
    allocate(pcf%vegfire(beg:end))
    pcf%vegfire(:)= nanr
    allocate(pcf%wood_harvestc(beg:end))
    pcf%wood_harvestc(:)= nanr
    allocate(pcf%pft_cinputs(beg:end))
    pcf%pft_cinputs(:)= nanr
    allocate(pcf%pft_coutputs(beg:end))
    pcf%pft_coutputs(:)= nanr
    allocate(pcf%pft_fire_closs(beg:end))
    pcf%pft_fire_closs(:)= nanr
    allocate(pcf%cpool_to_grainc(beg:end))
    pcf%cpool_to_grainc(:)= nanr
    allocate(pcf%cpool_to_grainc_storage(beg:end))
    pcf%cpool_to_grainc_storage(:)= nanr
    allocate(pcf%livestemc_to_litter(beg:end))
    pcf%livestemc_to_litter(:)= nanr
    allocate(pcf%grainc_to_food(beg:end))
    pcf%grainc_to_food(:)= nanr
    allocate(pcf%grainc_xfer_to_grainc(beg:end))
    pcf%grainc_xfer_to_grainc(:)= nanr
    allocate(pcf%cpool_grain_gr(beg:end))
    pcf%cpool_grain_gr(:)= nanr
    allocate(pcf%cpool_grain_storage_gr(beg:end))
    pcf%cpool_grain_storage_gr(:)= nanr
    allocate(pcf%transfer_grain_gr(beg:end))
    pcf%transfer_grain_gr(:)= nanr
    allocate(pcf%xsmrpool_to_atm(beg:end))
    pcf%xsmrpool_to_atm(:)= nanr
    allocate(pcf%grainc_storage_to_xfer(beg:end))
    pcf%grainc_storage_to_xfer(:)= nanr
    allocate(pcf%frootc_alloc(beg:end))
    pcf%frootc_alloc(:)= nanr
    allocate(pcf%frootc_loss(beg:end))
    pcf%frootc_loss(:)= nanr
    allocate(pcf%leafc_alloc(beg:end))
    pcf%leafc_alloc(:)= nanr
    allocate(pcf%leafc_loss(beg:end))
    pcf%leafc_loss(:)= nanr
    allocate(pcf%woodc_alloc(beg:end))
    pcf%woodc_alloc(:)= nanr
    allocate(pcf%woodc_loss(beg:end))
    pcf%woodc_loss(:)= nanr       
    allocate(pcf%tempavg_agnpp(beg:end))
    pcf%tempavg_agnpp(:)= spval
    allocate(pcf%tempavg_bgnpp(beg:end))
    pcf%tempavg_bgnpp(:)= spval
    allocate(pcf%annavg_agnpp(beg:end))
    pcf%annavg_agnpp(:)= spval ! To detect first year
    allocate(pcf%annavg_bgnpp(beg:end))
    pcf%annavg_bgnpp(:)= spval ! To detect first year

  end subroutine init_pft_cflux_type

  !------------------------------------------------------------------------
  subroutine init_pft_nflux_type(beg, end, pnf)
    !
    ! !DESCRIPTION:
    ! Initialize pft nitrogen flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nflux_type), intent(inout) :: pnf
    !------------------------------------------------------------------------

    allocate(pnf%m_leafn_to_litter(beg:end))
    pnf%m_leafn_to_litter(:)= nanr
    allocate(pnf%m_frootn_to_litter(beg:end))
    pnf%m_frootn_to_litter(:)= nanr
    allocate(pnf%m_leafn_storage_to_litter(beg:end))
    pnf%m_leafn_storage_to_litter(:)= nanr
    allocate(pnf%m_frootn_storage_to_litter(beg:end))
    pnf%m_frootn_storage_to_litter(:)= nanr
    allocate(pnf%m_livestemn_storage_to_litter(beg:end))
    pnf%m_livestemn_storage_to_litter(:)= nanr
    allocate(pnf%m_deadstemn_storage_to_litter(beg:end))
    pnf%m_deadstemn_storage_to_litter(:)= nanr
    allocate(pnf%m_livecrootn_storage_to_litter(beg:end))
    pnf%m_livecrootn_storage_to_litter(:)= nanr
    allocate(pnf%m_deadcrootn_storage_to_litter(beg:end))
    pnf%m_deadcrootn_storage_to_litter(:)= nanr
    allocate(pnf%m_leafn_xfer_to_litter(beg:end))
    pnf%m_leafn_xfer_to_litter(:)= nanr
    allocate(pnf%m_frootn_xfer_to_litter(beg:end))
    pnf%m_frootn_xfer_to_litter(:)= nanr
    allocate(pnf%m_livestemn_xfer_to_litter(beg:end))
    pnf%m_livestemn_xfer_to_litter(:)= nanr
    allocate(pnf%m_deadstemn_xfer_to_litter(beg:end))
    pnf%m_deadstemn_xfer_to_litter(:)= nanr
    allocate(pnf%m_livecrootn_xfer_to_litter(beg:end))
    pnf%m_livecrootn_xfer_to_litter(:)= nanr
    allocate(pnf%m_deadcrootn_xfer_to_litter(beg:end))
    pnf%m_deadcrootn_xfer_to_litter(:)= nanr
    allocate(pnf%m_livestemn_to_litter(beg:end))
    pnf%m_livestemn_to_litter(:)= nanr
    allocate(pnf%m_deadstemn_to_litter(beg:end))
    pnf%m_deadstemn_to_litter(:)= nanr
    allocate(pnf%m_livecrootn_to_litter(beg:end))
    pnf%m_livecrootn_to_litter(:)= nanr
    allocate(pnf%m_deadcrootn_to_litter(beg:end))
    pnf%m_deadcrootn_to_litter(:)= nanr
    allocate(pnf%m_retransn_to_litter(beg:end))
    pnf%m_retransn_to_litter(:)= nanr
    allocate(pnf%hrv_leafn_to_litter(beg:end))
    pnf%hrv_leafn_to_litter(:)= nanr
    allocate(pnf%hrv_frootn_to_litter(beg:end))
    pnf%hrv_frootn_to_litter(:)= nanr
    allocate(pnf%hrv_leafn_storage_to_litter(beg:end))
    pnf%hrv_leafn_storage_to_litter(:)= nanr
    allocate(pnf%hrv_frootn_storage_to_litter(beg:end))
    pnf%hrv_frootn_storage_to_litter(:)= nanr
    allocate(pnf%hrv_livestemn_storage_to_litter(beg:end))
    pnf%hrv_livestemn_storage_to_litter(:)= nanr
    allocate(pnf%hrv_deadstemn_storage_to_litter(beg:end))
    pnf%hrv_deadstemn_storage_to_litter(:)= nanr
    allocate(pnf%hrv_livecrootn_storage_to_litter(beg:end))
    pnf%hrv_livecrootn_storage_to_litter(:)= nanr
    allocate(pnf%hrv_deadcrootn_storage_to_litter(beg:end))
    pnf%hrv_deadcrootn_storage_to_litter(:)= nanr
    allocate(pnf%hrv_leafn_xfer_to_litter(beg:end))
    pnf%hrv_leafn_xfer_to_litter(:)= nanr
    allocate(pnf%hrv_frootn_xfer_to_litter(beg:end))
    pnf%hrv_frootn_xfer_to_litter(:)= nanr
    allocate(pnf%hrv_livestemn_xfer_to_litter(beg:end))
    pnf%hrv_livestemn_xfer_to_litter(:)= nanr
    allocate(pnf%hrv_deadstemn_xfer_to_litter(beg:end))
    pnf%hrv_deadstemn_xfer_to_litter(:)= nanr
    allocate(pnf%hrv_livecrootn_xfer_to_litter(beg:end))
    pnf%hrv_livecrootn_xfer_to_litter(:)= nanr
    allocate(pnf%hrv_deadcrootn_xfer_to_litter(beg:end))
    pnf%hrv_deadcrootn_xfer_to_litter(:)= nanr
    allocate(pnf%hrv_livestemn_to_litter(beg:end))
    pnf%hrv_livestemn_to_litter(:)= nanr
    allocate(pnf%hrv_deadstemn_to_prod10n(beg:end))
    pnf%hrv_deadstemn_to_prod10n(:)= nanr
    allocate(pnf%hrv_deadstemn_to_prod100n(beg:end))
    pnf%hrv_deadstemn_to_prod100n(:)= nanr
    allocate(pnf%hrv_livecrootn_to_litter(beg:end))
    pnf%hrv_livecrootn_to_litter(:)= nanr
    allocate(pnf%hrv_deadcrootn_to_litter(beg:end))
    pnf%hrv_deadcrootn_to_litter(:)= nanr
    allocate(pnf%hrv_retransn_to_litter(beg:end))
    pnf%hrv_retransn_to_litter(:)= nanr

    ! fire variables changed by F. Li and S. Levis            
    allocate(pnf%m_leafn_to_fire(beg:end))
    pnf%m_leafn_to_fire(:)= nanr
    allocate(pnf%m_leafn_storage_to_fire(beg:end))
    pnf%m_leafn_storage_to_fire(:)= nanr
    allocate(pnf%m_leafn_xfer_to_fire(beg:end))
    pnf%m_leafn_xfer_to_fire(:)= nanr
    allocate(pnf%m_livestemn_to_fire(beg:end))
    pnf%m_livestemn_to_fire(:)= nanr
    allocate(pnf%m_livestemn_storage_to_fire(beg:end))
    pnf%m_livestemn_storage_to_fire(:)= nanr
    allocate(pnf%m_livestemn_xfer_to_fire(beg:end))
    pnf%m_livestemn_xfer_to_fire(:)= nanr
    allocate(pnf%m_deadstemn_to_fire(beg:end))
    pnf%m_deadstemn_to_fire(:)= nanr
    allocate(pnf%m_deadstemn_storage_to_fire(beg:end))
    pnf%m_deadstemn_storage_to_fire(:)= nanr
    allocate(pnf%m_deadstemn_xfer_to_fire(beg:end))
    pnf%m_deadstemn_xfer_to_fire(:)= nanr
    allocate(pnf%m_frootn_to_fire(beg:end))
    pnf%m_frootn_to_fire(:)= nanr
    allocate(pnf%m_frootn_storage_to_fire(beg:end))
    pnf%m_frootn_storage_to_fire(:)= nanr
    allocate(pnf%m_frootn_xfer_to_fire(beg:end))
    pnf%m_frootn_xfer_to_fire(:)= nanr
    allocate(pnf%m_livecrootn_to_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_fire(beg:end))
    pnf%m_livecrootn_storage_to_fire(:)= nanr
    allocate(pnf%m_livecrootn_xfer_to_fire(beg:end))
    pnf%m_livecrootn_xfer_to_fire(:)= nanr
    allocate(pnf%m_deadcrootn_to_fire(beg:end))
    pnf%m_deadcrootn_to_fire(:)= nanr
    allocate(pnf%m_deadcrootn_storage_to_fire(beg:end))
    pnf%m_deadcrootn_storage_to_fire(:)= nanr
    allocate(pnf%m_deadcrootn_xfer_to_fire(beg:end))
    pnf%m_deadcrootn_xfer_to_fire(:)= nanr
    allocate(pnf%m_retransn_to_fire(beg:end))
    pnf%m_retransn_to_fire(:)= nanr

    allocate(pnf%m_leafn_to_litter_fire(beg:end))
    pnf%m_leafn_to_litter_fire(:)= nanr
    allocate(pnf%m_leafn_storage_to_litter_fire(beg:end))
    pnf%m_leafn_storage_to_litter_fire(:)= nanr
    allocate(pnf%m_leafn_xfer_to_litter_fire(beg:end))
    pnf%m_leafn_xfer_to_litter_fire(:)= nanr
    allocate(pnf%m_livestemn_to_litter_fire(beg:end))
    pnf%m_livestemn_to_litter_fire(:)= nanr
    allocate(pnf%m_livestemn_storage_to_litter_fire(beg:end))
    pnf%m_livestemn_storage_to_litter_fire(:)= nanr
    allocate(pnf%m_livestemn_xfer_to_litter_fire(beg:end))
    pnf%m_livestemn_xfer_to_litter_fire(:)= nanr
    allocate(pnf%m_livestemn_to_deadstemn_fire(beg:end))
    pnf%m_livestemn_to_deadstemn_fire(:)= nanr
    allocate(pnf%m_deadstemn_to_litter_fire(beg:end))
    pnf%m_deadstemn_to_litter_fire(:)= nanr
    allocate(pnf%m_deadstemn_storage_to_litter_fire(beg:end))
    pnf%m_deadstemn_storage_to_litter_fire(:)= nanr
    allocate(pnf%m_deadstemn_xfer_to_litter_fire(beg:end))
    pnf%m_deadstemn_xfer_to_litter_fire(:)= nanr
    allocate(pnf%m_frootn_to_litter_fire(beg:end))
    pnf%m_frootn_to_litter_fire(:)= nanr
    allocate(pnf%m_frootn_storage_to_litter_fire(beg:end))
    pnf%m_frootn_storage_to_litter_fire(:)= nanr
    allocate(pnf%m_frootn_xfer_to_litter_fire(beg:end))
    pnf%m_frootn_xfer_to_litter_fire(:)= nanr
    allocate(pnf%m_livecrootn_to_litter_fire(beg:end))
    pnf%m_livecrootn_to_litter_fire(:)= nanr
    allocate(pnf%m_livecrootn_storage_to_litter_fire(beg:end))
    pnf%m_livecrootn_storage_to_litter_fire(:)= nanr
    allocate(pnf%m_livecrootn_xfer_to_litter_fire(beg:end))
    pnf%m_livecrootn_xfer_to_litter_fire(:)= nanr
    allocate(pnf%m_livecrootn_to_deadcrootn_fire(beg:end))
    pnf%m_livecrootn_to_deadcrootn_fire(:)= nanr
    allocate(pnf%m_deadcrootn_to_litter_fire(beg:end))
    pnf%m_deadcrootn_to_litter_fire(:)= nanr
    allocate(pnf%m_deadcrootn_storage_to_litter_fire(beg:end))
    pnf%m_deadcrootn_storage_to_litter_fire(:)= nanr
    allocate(pnf%m_deadcrootn_xfer_to_litter_fire(beg:end))
    pnf%m_deadcrootn_xfer_to_litter_fire(:)= nanr
    allocate(pnf%m_retransn_to_litter_fire(beg:end))
    pnf%m_retransn_to_litter_fire(:)= nanr



    allocate(pnf%leafn_xfer_to_leafn(beg:end))
    pnf%leafn_xfer_to_leafn(:)= nanr
    allocate(pnf%frootn_xfer_to_frootn(beg:end))
    pnf%frootn_xfer_to_frootn(:)= nanr
    allocate(pnf%livestemn_xfer_to_livestemn(beg:end))
    pnf%livestemn_xfer_to_livestemn(:)= nanr
    allocate(pnf%deadstemn_xfer_to_deadstemn(beg:end))
    pnf%deadstemn_xfer_to_deadstemn(:)= nanr
    allocate(pnf%livecrootn_xfer_to_livecrootn(beg:end))
    pnf%livecrootn_xfer_to_livecrootn(:)= nanr
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(beg:end))
    pnf%deadcrootn_xfer_to_deadcrootn(:)= nanr
    allocate(pnf%leafn_to_litter(beg:end))
    pnf%leafn_to_litter(:)= nanr
    allocate(pnf%leafn_to_retransn(beg:end))
    pnf%leafn_to_retransn(:)= nanr
    allocate(pnf%frootn_to_retransn(beg:end))
    pnf%frootn_to_retransn(:)= nanr
    allocate(pnf%frootn_to_litter(beg:end))
    pnf%frootn_to_litter(:)= nanr
    allocate(pnf%retransn_to_npool(beg:end))
    pnf%retransn_to_npool(:)= nanr
    allocate(pnf%sminn_to_npool(beg:end))
    pnf%sminn_to_npool(:)= nanr
    allocate(pnf%npool_to_leafn(beg:end))
    pnf%npool_to_leafn(:)= nanr
    allocate(pnf%npool_to_leafn_storage(beg:end))
    pnf%npool_to_leafn_storage(:)= nanr
    allocate(pnf%npool_to_frootn(beg:end))
    pnf%npool_to_frootn(:)= nanr
    allocate(pnf%npool_to_frootn_storage(beg:end))
    pnf%npool_to_frootn_storage(:)= nanr
    allocate(pnf%npool_to_livestemn(beg:end))
    pnf%npool_to_livestemn(:)= nanr
    allocate(pnf%npool_to_livestemn_storage(beg:end))
    pnf%npool_to_livestemn_storage(:)= nanr
    allocate(pnf%npool_to_deadstemn(beg:end))
    pnf%npool_to_deadstemn(:)= nanr
    allocate(pnf%npool_to_deadstemn_storage(beg:end))
    pnf%npool_to_deadstemn_storage(:)= nanr
    allocate(pnf%npool_to_livecrootn(beg:end))
    pnf%npool_to_livecrootn(:)= nanr
    allocate(pnf%npool_to_livecrootn_storage(beg:end))
    pnf%npool_to_livecrootn_storage(:)= nanr
    allocate(pnf%npool_to_deadcrootn(beg:end))
    pnf%npool_to_deadcrootn(:)= nanr
    allocate(pnf%npool_to_deadcrootn_storage(beg:end))
    pnf%npool_to_deadcrootn_storage(:)= nanr
    allocate(pnf%leafn_storage_to_xfer(beg:end))
    pnf%leafn_storage_to_xfer(:)= nanr
    allocate(pnf%frootn_storage_to_xfer(beg:end))
    pnf%frootn_storage_to_xfer(:)= nanr
    allocate(pnf%livestemn_storage_to_xfer(beg:end))
    pnf%livestemn_storage_to_xfer(:)= nanr
    allocate(pnf%deadstemn_storage_to_xfer(beg:end))
    pnf%deadstemn_storage_to_xfer(:)= nanr
    allocate(pnf%livecrootn_storage_to_xfer(beg:end))
    pnf%livecrootn_storage_to_xfer(:)= nanr
    allocate(pnf%deadcrootn_storage_to_xfer(beg:end))
    pnf%deadcrootn_storage_to_xfer(:)= nanr
    allocate(pnf%livestemn_to_deadstemn(beg:end))
    pnf%livestemn_to_deadstemn(:)= nanr
    allocate(pnf%livestemn_to_retransn(beg:end))
    pnf%livestemn_to_retransn(:)= nanr
    allocate(pnf%livecrootn_to_deadcrootn(beg:end))
    pnf%livecrootn_to_deadcrootn(:)= nanr
    allocate(pnf%livecrootn_to_retransn(beg:end))
    pnf%livecrootn_to_retransn(:)= nanr
    allocate(pnf%ndeploy(beg:end))
    pnf%ndeploy(:)= nanr
    allocate(pnf%pft_ninputs(beg:end))
    pnf%pft_ninputs(:)= nanr
    allocate(pnf%pft_noutputs(beg:end))
    pnf%pft_noutputs(:)= nanr
    allocate(pnf%wood_harvestn(beg:end))
    pnf%wood_harvestn(:)= nanr
    allocate(pnf%pft_fire_nloss(beg:end))
    pnf%pft_fire_nloss(:)= nanr
    allocate(pnf%npool_to_grainn(beg:end))
    pnf%npool_to_grainn(:)= nanr
    allocate(pnf%npool_to_grainn_storage(beg:end))
    pnf%npool_to_grainn_storage(:)= nanr
    allocate(pnf%livestemn_to_litter(beg:end))
    pnf%livestemn_to_litter(:)= nanr
    allocate(pnf%grainn_to_food(beg:end))
    pnf%grainn_to_food(:)= nanr
    allocate(pnf%grainn_xfer_to_grainn(beg:end))
    pnf%grainn_xfer_to_grainn(:)= nanr
    allocate(pnf%grainn_storage_to_xfer(beg:end))
    pnf%grainn_storage_to_xfer(:)= nanr
    allocate(pnf%fert(beg:end))
    pnf%fert(:)= nanr
    allocate(pnf%soyfixn(beg:end))
    pnf%soyfixn(:)= nanr     

  end subroutine init_pft_nflux_type

  !------------------------------------------------------------------------
  subroutine init_pft_vflux_type(beg, end, pvf)
    !
    ! !DESCRIPTION:
    ! Initialize pft VOC flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vflux_type), intent(inout) :: pvf

    integer :: i
    !------------------------------------------------------------------------

    if (shr_megan_mechcomps_n<1) return

    allocate(pvf%vocflx_tot(beg:end))
    pvf%vocflx_tot(:)= nanr
    allocate(pvf%vocflx(beg:end,1:shr_megan_mechcomps_n))
    pvf%vocflx(:,:)= nanr
    allocate(pvf%Eopt_out(beg:end))
    pvf%Eopt_out(:)= nanr
    allocate(pvf%topt_out(beg:end))
    pvf%topt_out(:)= nanr
    allocate(pvf%alpha_out(beg:end))
    pvf%alpha_out(:)= nanr
    allocate(pvf%cp_out(beg:end))
    pvf%cp_out(:)= nanr
    allocate(pvf%para_out(beg:end))
    pvf%para_out(:)= nanr
    allocate(pvf%par24a_out(beg:end))
    pvf%par24a_out(:)= nanr
    allocate(pvf%par240a_out(beg:end))
    pvf%par240a_out(:)= nanr
    allocate(pvf%paru_out(beg:end))
    pvf%paru_out(:)= nanr
    allocate(pvf%par24u_out(beg:end))
    pvf%par24u_out(:)= nanr
    allocate(pvf%par240u_out(beg:end))
    pvf%par240u_out(:)= nanr
    allocate(pvf%gamma_out(beg:end))
    pvf%gamma_out(:)= nanr
    allocate(pvf%gammaL_out(beg:end))
    pvf%gammaL_out(:)= nanr
    allocate(pvf%gammaT_out(beg:end))
    pvf%gammaT_out(:)= nanr
    allocate(pvf%gammaP_out(beg:end))
    pvf%gammaP_out(:)= nanr
    allocate(pvf%gammaA_out(beg:end))
    pvf%gammaA_out(:)= nanr
    allocate(pvf%gammaS_out(beg:end))
    pvf%gammaS_out(:)= nanr
    allocate(pvf%gammaC_out(beg:end))
    pvf%gammaC_out(:)= nanr
    allocate(pvf%meg(shr_megan_megcomps_n))

    do i=1,shr_megan_megcomps_n
       allocate(pvf%meg(i)%flux_out(beg:end))
       pvf%meg(i)%flux_out(:)= nanr    
    enddo

  end subroutine init_pft_vflux_type

  !------------------------------------------------------------------------
  subroutine init_pft_dflux_type(beg, end, pdf)
    !
    ! !DESCRIPTION:
    ! Initialize pft dust flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dflux_type), intent(inout):: pdf
    !------------------------------------------------------------------------

    allocate(pdf%flx_mss_vrt_dst(beg:end,1:ndst))
    pdf%flx_mss_vrt_dst(:,:)= nanr
    allocate(pdf%flx_mss_vrt_dst_tot(beg:end))
    pdf%flx_mss_vrt_dst_tot(:)= nanr
    allocate(pdf%vlc_trb(beg:end,1:ndst))
    pdf%vlc_trb(:,:)= nanr
    allocate(pdf%vlc_trb_1(beg:end))
    pdf%vlc_trb_1(:)= nanr
    allocate(pdf%vlc_trb_2(beg:end))
    pdf%vlc_trb_2(:)= nanr
    allocate(pdf%vlc_trb_3(beg:end))
    pdf%vlc_trb_3(:)= nanr
    allocate(pdf%vlc_trb_4(beg:end))
    pdf%vlc_trb_4(:)= nanr

  end subroutine init_pft_dflux_type

  !------------------------------------------------------------------------
  subroutine init_pft_depvd_type(beg, end, pdd)
    !
    ! !DESCRIPTION:
    ! Initialize pft dep velocity variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_depvd_type), intent(inout):: pdd
    integer :: i
    !------------------------------------------------------------------------

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(pdd%drydepvel(beg:end,n_drydep))
       pdd%drydepvel(:,:)=nanr
    end if

  end subroutine init_pft_depvd_type

!------------------------------------------------------------------------
  subroutine init_column_pstate_type(beg, end, cps)
    !
    ! !DESCRIPTION:
    ! Initialize column physical state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_pstate_type), intent(inout):: cps
    !------------------------------------------------------------------------

    allocate(cps%snl(beg:end))      !* cannot be averaged up
    allocate(cps%snow_layer_unity(beg:end,-nlevsno+1:0))
    cps%snow_layer_unity(:,:) = 1._r8

    allocate(cps%isoicol(beg:end))  !* cannot be averaged up

    !F. Li and S. Levis
    allocate(cps%gdp_lf(beg:end))  
    allocate(cps%peatf_lf(beg:end))  
    allocate(cps%abm_lf(beg:end))  
    allocate(cps%lgdp_col(beg:end))  
    allocate(cps%lgdp1_col(beg:end))
    allocate(cps%lpop_col(beg:end))  

    allocate(cps%bsw(beg:end,nlevgrnd))
    cps%bsw(:,:)= nanr
    allocate(cps%watsat(beg:end,nlevgrnd))
    cps%watsat(:,:)= nanr
    allocate(cps%watfc(beg:end,nlevgrnd))
    cps%watfc(:,:)= nanr
    allocate(cps%watdry(beg:end,nlevgrnd))
    cps%watdry(:,:)= nanr
    allocate(cps%watopt(beg:end,nlevgrnd))
    cps%watopt(:,:)= nanr
    allocate(cps%hksat(beg:end,nlevgrnd))
    cps%hksat(:,:)= nanr
    allocate(cps%sucsat(beg:end,nlevgrnd))
    cps%sucsat(:,:)= nanr
    allocate(cps%thk(beg:end, -nlevsno+1:nlevgrnd))
    cps%thk(:,:)= spval
    allocate(cps%csol(beg:end,nlevgrnd))
    cps%csol(:,:)= nanr
    allocate(cps%tkmg(beg:end,nlevgrnd))
    cps%tkmg(:,:)= nanr
    allocate(cps%tkdry(beg:end,nlevgrnd))
    cps%tkdry(:,:)= nanr
    allocate(cps%tksatu(beg:end,nlevgrnd))
    cps%tksatu(:,:)= nanr
    allocate(cps%smpmin(beg:end))
    cps%smpmin(:)= nanr
    allocate(cps%hkdepth(beg:end))
    cps%hkdepth(:)= nanr
    allocate(cps%wtfact(beg:end))
    cps%wtfact(:)= nanr
    allocate(cps%fracice(beg:end,nlevgrnd))
    cps%fracice(:,:)= nanr
    allocate(cps%icefrac(beg:end,nlevgrnd))
    cps%icefrac(:,:)= nanr
    allocate(cps%gwc_thr(beg:end))
    cps%gwc_thr(:)= nanr
    allocate(cps%mss_frc_cly_vld(beg:end))
    cps%mss_frc_cly_vld(:)= nanr
    allocate(cps%mbl_bsn_fct(beg:end))
    cps%mbl_bsn_fct(:)= nanr
    allocate(cps%do_capsnow(beg:end))
    allocate(cps%snow_depth(beg:end))
    cps%snow_depth(:)= nanr
    allocate(cps%snow_persistence(beg:end))
    cps%snow_persistence(:)= nanr   
    allocate(cps%snowdp(beg:end))
    cps%snowdp(:)= nanr
    allocate(cps%frac_sno(beg:end))
    cps%frac_sno(:)=nanr
    allocate(cps%zi(beg:end,-nlevsno+0:nlevgrnd))
    cps%zi(:,:)= nanr
    allocate(cps%dz(beg:end,-nlevsno+1:nlevgrnd))
    cps%dz(:,:)= nanr
    allocate(cps%z (beg:end,-nlevsno+1:nlevgrnd))
    cps%z (:,:)= nanr
    allocate(cps%frac_iceold(beg:end,-nlevsno+1:nlevgrnd))
    cps%frac_iceold(:,:)= spval
    allocate(cps%imelt(beg:end,-nlevsno+1:nlevgrnd))
    cps%imelt(:,:)= huge(1)
    allocate(cps%eff_porosity(beg:end,nlevgrnd))
    cps%eff_porosity(:,:)= spval
    allocate(cps%emg(beg:end))
    cps%emg(:)= nanr
    allocate(cps%z0mg(beg:end))
    cps%z0mg(:)= nanr
    allocate(cps%z0hg(beg:end))
    cps%z0hg(:)= nanr
    allocate(cps%z0qg(beg:end))
    cps%z0qg(:)= nanr
    allocate(cps%htvp(beg:end))
    cps%htvp(:)= nanr
    allocate(cps%beta(beg:end))
    cps%beta(:)= nanr
    allocate(cps%zii(beg:end))
    cps%zii(:)= nanr
    allocate(cps%albgrd(beg:end,numrad))
    cps%albgrd(:,:)= nanr
    allocate(cps%albgri(beg:end,numrad))
    cps%albgri(:,:)= nanr
    allocate(cps%rootr_column(beg:end,nlevgrnd))
    cps%rootr_column(:,:)= spval
    allocate(cps%rootfr_road_perv(beg:end,nlevgrnd))
    cps%rootfr_road_perv(:,:)= nanr
    allocate(cps%rootr_road_perv(beg:end,nlevgrnd))
    cps%rootr_road_perv(:,:)= nanr
    allocate(cps%wf(beg:end))
    cps%wf(:)= nanr
    allocate(cps%wf2(beg:end))
    allocate(cps%soilpsi(beg:end,nlevgrnd))
    cps%soilpsi(:,:)= spval
    allocate(cps%coszen(beg:end))
    cps%coszen(:)= nanr
    allocate(cps%bd(beg:end,nlevgrnd))
    cps%bd(:,:)= spval
    allocate(cps%fpi(beg:end))
    cps%fpi(:)= nanr
    allocate(cps%fpi_vr(beg:end,1:nlevdecomp_full))
    cps%fpi_vr(:,:)= nanr
    allocate(cps%fpg(beg:end))
    cps%fpg(:)= nanr
    allocate(cps%annsum_counter(beg:end))
    cps%annsum_counter(:)= nanr
    allocate(cps%cannsum_npp(beg:end))
    cps%cannsum_npp(:)= nanr
    allocate(cps%col_lag_npp(beg:end))
    cps%col_lag_npp(:)= spval
    allocate(cps%cannavg_t2m(beg:end))
    cps%cannavg_t2m(:)= nanr

    allocate(cps%nfire(beg:end))
    cps%nfire(:)= spval
    allocate(cps%farea_burned(beg:end))
    cps%farea_burned(:)= nanr
    allocate(cps%fsr_col(beg:end))
    cps%fsr_col(:)= nanr
    allocate(cps%fd_col(beg:end))
    cps%fd_col(:)= nanr
    allocate(cps%cropf_col(beg:end))
    cps%cropf_col(:)= nanr
    allocate(cps%prec10_col(beg:end))
    cps%prec10_col(:)= nanr
    allocate(cps%prec60_col(beg:end))
    cps%prec60_col(:)= nanr
    allocate(cps%lfc(beg:end))
    cps%lfc(:)= spval
    allocate(cps%lfc2(beg:end))
    cps%lfc2(:)= 0._r8
    allocate(cps%trotr1_col(beg:end))
    cps%trotr1_col(:)= 0._r8
    allocate(cps%trotr2_col(beg:end))
    cps%trotr2_col(:)= 0._r8
    allocate(cps%dtrotr_col(beg:end))
    cps%dtrotr_col(:)= 0._r8
    allocate(cps%baf_crop(beg:end))
    cps%baf_crop(:)= nanr
    allocate(cps%baf_peatf(beg:end))
    cps%baf_peatf(:)= nanr
    allocate(cps%fbac(beg:end))
    cps%fbac(:)= nanr
    allocate(cps%fbac1(beg:end))
    cps%fbac1(:)= nanr
    allocate(cps%btran_col(beg:end))
    cps%btran_col(:)= nanr
    allocate(cps%wtlf(beg:end))
    cps%wtlf(:)= nanr
    allocate(cps%lfwt(beg:end))
    cps%lfwt(:)= nanr

    allocate(cps%albsnd_hst(beg:end,numrad))
    cps%albsnd_hst(:,:)= spval
    allocate(cps%albsni_hst(beg:end,numrad))
    cps%albsni_hst(:,:)= spval
    allocate(cps%albsod(beg:end,numrad))
    cps%albsod(:,:)= spval
    allocate(cps%albsoi(beg:end,numrad))
    cps%albsoi(:,:)= spval
    allocate(cps%flx_absdv(beg:end,-nlevsno+1:1))
    cps%flx_absdv(:,:)= spval
    allocate(cps%flx_absdn(beg:end,-nlevsno+1:1))
    cps%flx_absdn(:,:)= spval
    allocate(cps%flx_absiv(beg:end,-nlevsno+1:1))
    cps%flx_absiv(:,:)= spval
    allocate(cps%flx_absin(beg:end,-nlevsno+1:1))
    cps%flx_absin(:,:)= spval
    allocate(cps%snw_rds(beg:end,-nlevsno+1:0))
    cps%snw_rds(:,:)= nanr
    allocate(cps%snw_rds_top(beg:end))
    cps%snw_rds_top(:)= nanr
    allocate(cps%sno_liq_top(beg:end))
    cps%sno_liq_top(:)= nanr
    allocate(cps%mss_bcpho(beg:end,-nlevsno+1:0))
    cps%mss_bcpho(:,:)= nanr
    allocate(cps%mss_bcphi(beg:end,-nlevsno+1:0))
    cps%mss_bcphi(:,:)= nanr
    allocate(cps%mss_bctot(beg:end,-nlevsno+1:0))
    cps%mss_bctot(:,:)= nanr
    allocate(cps%mss_bc_col(beg:end))
    cps%mss_bc_col(:)= nanr
    allocate(cps%mss_bc_top(beg:end))
    cps%mss_bc_top(:)= nanr
    allocate(cps%mss_ocpho(beg:end,-nlevsno+1:0))
    cps%mss_ocpho(:,:)= nanr
    allocate(cps%mss_ocphi(beg:end,-nlevsno+1:0))
    cps%mss_ocphi(:,:)= nanr
    allocate(cps%mss_octot(beg:end,-nlevsno+1:0))
    cps%mss_octot(:,:)= nanr
    allocate(cps%mss_oc_col(beg:end))
    cps%mss_oc_col(:)= nanr
    allocate(cps%mss_oc_top(beg:end))
    cps%mss_oc_top(:)= nanr
    allocate(cps%mss_dst1(beg:end,-nlevsno+1:0))
    cps%mss_dst1(:,:)= nanr
    allocate(cps%mss_dst2(beg:end,-nlevsno+1:0))
    cps%mss_dst2(:,:)= nanr
    allocate(cps%mss_dst3(beg:end,-nlevsno+1:0))
    cps%mss_dst3(:,:)= nanr
    allocate(cps%mss_dst4(beg:end,-nlevsno+1:0))
    cps%mss_dst4(:,:)= nanr
    allocate(cps%mss_dsttot(beg:end,-nlevsno+1:0))
    cps%mss_dsttot(:,:)= nanr
    allocate(cps%mss_dst_col(beg:end))
    cps%mss_dst_col(:)= nanr
    allocate(cps%mss_dst_top(beg:end))
    cps%mss_dst_top(:)= nanr
    allocate(cps%h2osno_top(beg:end))
    cps%h2osno_top(:)= nanr
    allocate(cps%mss_cnc_bcphi(beg:end,-nlevsno+1:0))
    cps%mss_cnc_bcphi(:,:)= nanr
    allocate(cps%mss_cnc_bcpho(beg:end,-nlevsno+1:0))
    cps%mss_cnc_bcpho(:,:)= nanr
    allocate(cps%mss_cnc_ocphi(beg:end,-nlevsno+1:0))
    cps%mss_cnc_ocphi(:,:)= nanr
    allocate(cps%mss_cnc_ocpho(beg:end,-nlevsno+1:0))
    cps%mss_cnc_ocpho(:,:)= nanr
    allocate(cps%mss_cnc_dst1(beg:end,-nlevsno+1:0))
    cps%mss_cnc_dst1(:,:)= nanr
    allocate(cps%mss_cnc_dst2(beg:end,-nlevsno+1:0))
    cps%mss_cnc_dst2(:,:)= nanr
    allocate(cps%mss_cnc_dst3(beg:end,-nlevsno+1:0))
    cps%mss_cnc_dst3(:,:)= nanr
    allocate(cps%mss_cnc_dst4(beg:end,-nlevsno+1:0))
    cps%mss_cnc_dst4(:,:)= nanr
    allocate(cps%albgrd_pur(beg:end,numrad))
    cps%albgrd_pur(:,:)= nanr
    allocate(cps%albgri_pur(beg:end,numrad))
    cps%albgri_pur(:,:)= nanr
    allocate(cps%albgrd_bc(beg:end,numrad))
    cps%albgrd_bc(:,:)= nanr
    allocate(cps%albgri_bc(beg:end,numrad))
    cps%albgri_bc(:,:)= nanr
    allocate(cps%albgrd_oc(beg:end,numrad))
    cps%albgrd_oc(:,:)= nanr
    allocate(cps%albgri_oc(beg:end,numrad))
    cps%albgri_oc(:,:)= nanr
    allocate(cps%albgrd_dst(beg:end,numrad))
    cps%albgrd_dst(:,:)= nanr
    allocate(cps%albgri_dst(beg:end,numrad))
    cps%albgri_dst(:,:)= nanr
    allocate(cps%dTdz_top(beg:end))
    cps%dTdz_top(:)= nanr
    allocate(cps%snot_top(beg:end))
    cps%snot_top(:)= nanr

    ! New variables for "S" Lakes
    allocate(cps%ws(beg:end))
    cps%ws(:)= nanr
    allocate(cps%ks(beg:end))
    cps%ks(:)= nanr
    allocate(cps%dz_lake(beg:end,nlevlak))
    cps%dz_lake(:,:)= nanr
    allocate(cps%z_lake(beg:end,nlevlak))
    cps%z_lake(:,:)= nanr
    allocate(cps%savedtke1(beg:end))
    cps%savedtke1(:)= spval  ! Initialize to spval so that c->g averaging will be done properl
    allocate(cps%cellsand(beg:end,nlevsoi))
    cps%cellsand(:,:)= nanr
    allocate(cps%cellclay(beg:end,nlevsoi))
    cps%cellclay(:,:)= nanr
    allocate(cps%cellorg(beg:end,nlevsoi))
    cps%cellorg(:,:)= nanr
    allocate(cps%lakedepth(beg:end))
    cps%lakedepth(:)= spval  ! Initialize to spval so that it can be a placeholde
    allocate(cps%etal(beg:end))
    cps%etal(:)= nanr
    allocate(cps%lakefetch(beg:end))
    cps%lakefetch(:)= nanr
    allocate(cps%ust_lake(beg:end))
    cps%ust_lake(:)= spval   ! Initial to spval to detect input from restart file if not arbini

    ! New variables for VIC hydrology 
    allocate(cps%b_infil(beg:end))
    cps%b_infil(:)= nanr
    allocate(cps%dsmax(beg:end))
    cps%dsmax(:)= nanr
    allocate(cps%ds(beg:end))
    cps%ds(:)= nanr
    allocate(cps%Wsvic(beg:end))
    cps%Wsvic(:)= nanr
    allocate(cps%c_param(beg:end))
    cps%c_param(:)= nanr
    allocate(cps%expt(beg:end, nlayer))
    cps%expt(:,:)= nanr
    allocate(cps%ksat(beg:end, nlayer))
    cps%ksat(:,:)= nanr
    allocate(cps%phi_s(beg:end, nlayer))
    cps%phi_s(:,:)= nanr
    allocate(cps%depth(beg:end, nlayert))
    cps%depth(:,:)= nanr
    allocate(cps%porosity(beg:end, nlayer))
    cps%porosity(:,:)= spval
    allocate(cps%max_moist(beg:end, nlayer))
    cps%max_moist(:,:)= nanr
    allocate(cps%vic_clm_fract(beg:end, nlayer, nlevsoi))
    cps%vic_clm_fract(:,:,:) = nanr

    ! New variable for finundated parameterization
    allocate(cps%zwt0(beg:end))
    cps%zwt0(:)= nanr
    allocate(cps%f0(beg:end))
    cps%f0(:)= nanr
    allocate(cps%p3(beg:end))
    cps%p3(:)= nanr

    ! New variable for methane
    allocate(cps%pH(beg:end))
    cps%pH(:)= nanr

    allocate(cps%glc_topo(beg:end))
    cps%glc_topo(:)= nanr
    allocate(cps%sub_surf_abs_SW(beg:end))
    cps%sub_surf_abs_SW(:)= nanr

    allocate(cps%rf_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    cps%rf_decomp_cascade(:,:,:)= nanr
    allocate(cps%pathfrac_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    cps%pathfrac_decomp_cascade(:,:,:)= nanr
    allocate(cps%nfixation_prof(beg:end,1:nlevdecomp_full))
    cps%nfixation_prof(:,:)= spval
    allocate(cps%ndep_prof(beg:end,1:nlevdecomp_full))
    cps%ndep_prof(:,:)= spval
    allocate(cps%alt(beg:end))
    cps%alt(:)= spval
    allocate(cps%altmax(beg:end))
    cps%altmax(:)= spval
    allocate(cps%altmax_lastyear(beg:end))
    cps%altmax_lastyear(:)= spval
    allocate(cps%alt_indx(beg:end))
    cps%alt_indx(:)= huge(1)
    allocate(cps%altmax_indx(beg:end))
    cps%altmax_indx(:)= huge(1)
    allocate(cps%altmax_lastyear_indx(beg:end))
    cps%altmax_lastyear_indx(:)= huge(1)
    allocate(cps%som_adv_coef(beg:end,1:nlevdecomp_full))
    cps%som_adv_coef(:,:)= spval
    allocate(cps%som_diffus_coef(beg:end,1:nlevdecomp_full))
    cps%som_diffus_coef(:,:)= spval

    allocate(cps%frac_sno_eff(beg:end))
    cps%frac_sno_eff(:)= spval
    allocate(cps%topo_std(beg:end))
    cps%topo_std(:)= nanr
    allocate(cps%topo_ndx(beg:end))
    cps%topo_ndx(:)= nanr
    allocate(cps%topo_slope(beg:end))
    cps%topo_slope(:)= nanr
    allocate(cps%hksat_min(beg:end,nlevgrnd))
    cps%hksat_min(:,:)= nanr
    allocate(cps%frac_h2osfc(beg:end))
    cps%frac_h2osfc(:)= spval
    allocate(cps%micro_sigma(beg:end))
    cps%micro_sigma(:)= nanr
    allocate(cps%h2osfc_thresh(beg:end))
    cps%h2osfc_thresh(:)= nanr
    allocate(cps%n_melt(beg:end))
    cps%n_melt(:)= nanr 

  end subroutine init_column_pstate_type

  !------------------------------------------------------------------------
  subroutine init_column_estate_type(beg, end, ces)
    !
    ! !DESCRIPTION:
    ! Initialize column energy state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_estate_type), intent(inout):: ces
    !------------------------------------------------------------------------

    allocate(ces%t_grnd(beg:end))
    ces%t_grnd(:)= nanr
    allocate(ces%t_grnd_u(beg:end))
    ces%t_grnd_u(:)= nanr
    allocate(ces%t_grnd_r(beg:end))
    ces%t_grnd_r(:)= nanr
    allocate(ces%dt_grnd(beg:end))
    ces%dt_grnd(:)= nanr
    allocate(ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd))
    ces%t_soisno(:,:)= spval
    allocate(ces%t_soi_10cm(beg:end))
    ces%t_soi_10cm(:)= spval
    allocate(ces%tsoi17(beg:end))
    ces%tsoi17(:)= spval
    allocate(ces%t_lake(beg:end,1:nlevlak))
    ces%t_lake(:,:)= nanr
    allocate(ces%tssbef(beg:end,-nlevsno+1:nlevgrnd))
    ces%tssbef(:,:)= nanr
    allocate(ces%thv(beg:end))
    ces%thv(:)= nanr
    allocate(ces%hc_soi(beg:end))
    ces%hc_soi(:)= nanr
    allocate(ces%hc_soisno(beg:end))
    ces%hc_soisno(:)= nanr
    allocate(ces%t_h2osfc(beg:end))
    ces%t_h2osfc(:)= spval
    allocate(ces%t_h2osfc_bef(beg:end))
    ces%t_h2osfc_bef(:)= nanr

  end subroutine init_column_estate_type

!------------------------------------------------------------------------
  subroutine init_column_wstate_type(beg, end, cws)
    !
    ! !DESCRIPTION:
    ! Initialize column water state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wstate_type), intent(inout):: cws !column water state
    !------------------------------------------------------------------------

    allocate(cws%h2osno(beg:end))
    cws%h2osno(:)= nanr
    allocate(cws%errh2osno(beg:end))
    cws%errh2osno(:)= nanr
    allocate(cws%snow_sources(beg:end))
    cws%snow_sources(:)= nanr
    allocate(cws%snow_sinks(beg:end))
    cws%snow_sinks(:)= nanr
    allocate(cws%h2osoi_liq(beg:end,-nlevsno+1:nlevgrnd))
    cws%h2osoi_liq(:,:)= spval
    allocate(cws%h2osoi_ice(beg:end,-nlevsno+1:nlevgrnd))
    cws%h2osoi_ice(:,:)= spval
    allocate(cws%h2osoi_liqice_10cm(beg:end))
    cws%h2osoi_liqice_10cm(:)= spval
    allocate(cws%h2osoi_vol(beg:end,1:nlevgrnd))
    cws%h2osoi_vol(:,:)= spval
    allocate(cws%bw(beg:end,-nlevsno+1:0))
    cws%bw(:,:)= spval
    allocate(cws%h2osno_old(beg:end))
    cws%h2osno_old(:)= nanr
    allocate(cws%qg(beg:end))
    cws%qg(:)= nanr
    allocate(cws%dqgdT(beg:end))
    cws%dqgdT(:)= nanr
    allocate(cws%snowice(beg:end))
    cws%snowice(:)= nanr
    allocate(cws%snowliq(beg:end))
    cws%snowliq(:)= nanr
    allocate(cws%soilalpha(beg:end))
    cws%soilalpha(:)= nanr
    allocate(cws%soilbeta(beg:end))
    cws%soilbeta(:)= nanr
    allocate(cws%soilalpha_u(beg:end))
    cws%soilalpha_u(:)= nanr
    allocate(cws%zwt(beg:end))
    cws%zwt(:)= nanr
    allocate(cws%fcov(beg:end))
    cws%fcov(:)= nanr
    allocate(cws%fsat(beg:end))
    cws%fsat(:)= nanr

    !New variable for methane code
    allocate(cws%finundated(beg:end))
    cws%finundated(:)= nanr
    allocate(cws%wa(beg:end))
    cws%wa(:)= spval
    allocate(cws%qcharge(beg:end))
    cws%qcharge(:)= nanr
    allocate(cws%smp_l(beg:end,1:nlevgrnd))
    cws%smp_l(:,:)= spval
    allocate(cws%hk_l(beg:end,1:nlevgrnd))
    cws%hk_l(:,:)= spval

    ! New variables for "S" lakes
    allocate(cws%lake_icefrac(beg:end,1:nlevlak))
    cws%lake_icefrac(:,:)= spval
    allocate(cws%lake_icethick(beg:end))
    cws%lake_icethick(:)= nanr

    allocate(cws%moist(beg:end,1:nlayert))
    cws%moist(:,:)= spval
    allocate(cws%ice(beg:end,1:nlayert))
    cws%ice(:,:)= spval
    allocate(cws%moist_vol(beg:end,1:nlayert))
    cws%moist_vol(:,:)= spval
    allocate(cws%max_infil(beg:end))
    cws%max_infil(:)= spval
    allocate(cws%i_0(beg:end))
    cws%i_0(:)= spval

    allocate(cws%h2osfc(beg:end))
    cws%h2osfc(:)= spval
    allocate(cws%qg_snow(beg:end))
    cws%qg_snow(:)= nanr
    allocate(cws%qg_soil(beg:end))
    cws%qg_soil(:)= nanr
    allocate(cws%qg_h2osfc(beg:end))
    cws%qg_h2osfc(:)= nanr
    allocate(cws%frost_table(beg:end))
    cws%frost_table(:)= spval
    allocate(cws%zwt_perched(beg:end))
    cws%zwt_perched(:)= spval
    allocate(cws%int_snow(beg:end))
    cws%int_snow(:)= spval
    allocate(cws%swe_old(beg:end,-nlevsno+1:0))
    cws%swe_old(:,:)= nanr      

  end subroutine init_column_wstate_type

  !------------------------------------------------------------------------
  subroutine init_column_cstate_type(beg, end, ccs)
    !
    ! !DESCRIPTION:
    ! Initialize column carbon state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cstate_type), intent(inout):: ccs
    !------------------------------------------------------------------------

    allocate(ccs%soilc(beg:end))
    ccs%soilc(:)= nanr
    allocate(ccs%cwdc(beg:end))
    ccs%cwdc(:)= nanr
    allocate(ccs%col_ctrunc(beg:end))
    ccs%col_ctrunc(:)= nanr
    allocate(ccs%decomp_cpools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    ccs%decomp_cpools_vr(:,:,:)= nanr
    allocate(ccs%decomp_cpools(beg:end,1:ndecomp_pools))
    ccs%decomp_cpools(:,:)= nanr
    allocate(ccs%decomp_cpools_1m(beg:end,1:ndecomp_pools))
    ccs%decomp_cpools_1m(:,:)= nanr
    allocate(ccs%col_ctrunc_vr(beg:end,1:nlevdecomp_full))
    ccs%col_ctrunc_vr(:,:)= nanr
    allocate(ccs%seedc(beg:end))
    ccs%seedc(:)= nanr
    allocate(ccs%prod10c(beg:end))
    ccs%prod10c(:)= nanr
    allocate(ccs%prod100c(beg:end))
    ccs%prod100c(:)= nanr
    allocate(ccs%totprodc(beg:end))
    ccs%totprodc(:)= nanr
    allocate(ccs%totlitc(beg:end))
    ccs%totlitc(:)= nanr
    allocate(ccs%totsomc(beg:end))
    ccs%totsomc(:)= nanr
    allocate(ccs%totlitc_1m(beg:end))
    ccs%totlitc_1m(:)= nanr
    allocate(ccs%totsomc_1m(beg:end))
    ccs%totsomc_1m(:)= nanr
    allocate(ccs%totecosysc(beg:end))
    ccs%totecosysc(:)= nanr
    allocate(ccs%totcolc(beg:end))
    ccs%totcolc(:)= nanr

    !F. Li and S. Levis
    allocate(ccs%rootc_col(beg:end))
    ccs%rootc_col(:)= nanr
    allocate(ccs%totvegc_col(beg:end))
    ccs%totvegc_col(:)= nanr
    allocate(ccs%leafc_col(beg:end))
    ccs%leafc_col(:)= nanr
    allocate(ccs%fuelc(beg:end))
    ccs%fuelc(:)= spval
    allocate(ccs%fuelc_crop(beg:end))
    ccs%fuelc_crop(:)= nanr

  end subroutine init_column_cstate_type

  !------------------------------------------------------------------------
  subroutine init_column_nstate_type(beg, end, cns)
    !
    ! !DESCRIPTION:
    ! Initialize column nitrogen state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nstate_type), intent(inout):: cns
    !------------------------------------------------------------------------

    allocate(cns%decomp_npools(beg:end,1:ndecomp_pools))
    cns%decomp_npools(:,:)= nanr
    allocate(cns%decomp_npools_1m(beg:end,1:ndecomp_pools))
    cns%decomp_npools_1m(:,:)= nanr
    allocate(cns%decomp_npools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    cns%decomp_npools_vr(:,:,:)= nanr
    allocate(cns%sminn_vr(beg:end,1:nlevdecomp_full))
    cns%sminn_vr(:,:)= nanr
    allocate(cns%col_ntrunc_vr(beg:end,1:nlevdecomp_full))
    cns%col_ntrunc_vr(:,:)= nanr
    allocate(cns%smin_no3_vr(beg:end,1:nlevdecomp_full))
    cns%smin_no3_vr(:,:)= nanr
    allocate(cns%smin_nh4_vr(beg:end,1:nlevdecomp_full))
    cns%smin_nh4_vr(:,:)= nanr
    allocate(cns%smin_no3(beg:end))
    cns%smin_no3(:)= nanr
    allocate(cns%smin_nh4(beg:end))
    cns%smin_nh4(:)= nanr
    allocate(cns%cwdn(beg:end))
    cns%cwdn(:)= nanr
    allocate(cns%sminn(beg:end))
    cns%sminn(:)= nanr
    allocate(cns%col_ntrunc(beg:end))
    cns%col_ntrunc(:)= nanr
    allocate(cns%seedn(beg:end))
    cns%seedn(:)= nanr
    allocate(cns%prod10n(beg:end))
    cns%prod10n(:)= nanr
    allocate(cns%prod100n(beg:end))
    cns%prod100n(:)= nanr
    allocate(cns%totprodn(beg:end))
    cns%totprodn(:)= nanr
    allocate(cns%totlitn(beg:end))
    cns%totlitn(:)= nanr
    allocate(cns%totsomn(beg:end))
    cns%totsomn(:)= nanr
    allocate(cns%totlitn_1m(beg:end))
    cns%totlitn_1m(:)= nanr
    allocate(cns%totsomn_1m(beg:end))
    cns%totsomn_1m(:)= nanr
    allocate(cns%totecosysn(beg:end))
    cns%totecosysn(:)= nanr
    allocate(cns%totcoln(beg:end))
    cns%totcoln(:)= nanr

  end subroutine init_column_nstate_type

  !------------------------------------------------------------------------
  subroutine init_column_eflux_type(beg, end, cef)
    !
    ! !DESCRIPTION:
    ! Initialize column energy flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_eflux_type), intent(inout):: cef
    !------------------------------------------------------------------------

    allocate(cef%eflx_snomelt(beg:end))
    cef%eflx_snomelt(:)= spval
    allocate(cef%eflx_snomelt_u(beg:end))
    cef%eflx_snomelt_u(:)= spval
    allocate(cef%eflx_snomelt_r(beg:end))
    cef%eflx_snomelt_r(:)= spval
    allocate(cef%eflx_impsoil(beg:end))
    cef%eflx_impsoil(:)= nanr
    allocate(cef%eflx_fgr12(beg:end))
    cef%eflx_fgr12(:)= nanr
    allocate(cef%eflx_fgr(beg:end, 1:nlevgrnd))
    cef%eflx_fgr(:,:)= nanr
    allocate(cef%eflx_building_heat(beg:end))
    cef%eflx_building_heat(:)= nanr
    allocate(cef%eflx_urban_ac(beg:end))
    cef%eflx_urban_ac(:)= nanr
    allocate(cef%eflx_urban_heat(beg:end))
    cef%eflx_urban_heat(:)= nanr
    allocate(cef%eflx_bot(beg:end))
    cef%eflx_bot(:)= nanr

  end subroutine init_column_eflux_type

  !------------------------------------------------------------------------
  subroutine init_column_wflux_type(beg, end, cwf)
    !
    ! !DESCRIPTION:
    ! Initialize column water flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wflux_type), intent(inout):: cwf
    !------------------------------------------------------------------------

    allocate(cwf%qflx_infl(beg:end))
    cwf%qflx_infl(:)= nanr
    allocate(cwf%qflx_surf(beg:end))
    cwf%qflx_surf(:)= nanr
    allocate(cwf%qflx_drain(beg:end))
    cwf%qflx_drain(:)= nanr
    allocate(cwf%qflx_top_soil(beg:end))
    cwf%qflx_top_soil(:)= spval
    allocate(cwf%qflx_sl_top_soil(beg:end))
    cwf%qflx_sl_top_soil(:)= nanr
    allocate(cwf%qflx_snomelt(beg:end))
    cwf%qflx_snomelt(:)= nanr
    allocate(cwf%qflx_qrgwl(beg:end))
    cwf%qflx_qrgwl(:)= nanr
    allocate(cwf%qflx_runoff(beg:end))
    cwf%qflx_runoff(:)= nanr
    allocate(cwf%qflx_runoff_u(beg:end))
    cwf%qflx_runoff_u(:)= nanr
    allocate(cwf%qflx_runoff_r(beg:end))
    cwf%qflx_runoff_r(:)= nanr
    allocate(cwf%qmelt(beg:end))
    cwf%qmelt(:)= nanr
    allocate(cwf%qflx_rsub_sat(beg:end))
    cwf%qflx_rsub_sat(:)= spval
    allocate(cwf%flx_bc_dep_dry(beg:end))
    cwf%flx_bc_dep_dry(:)= nanr
    allocate(cwf%flx_bc_dep_wet(beg:end))
    cwf%flx_bc_dep_wet(:)= nanr
    allocate(cwf%flx_bc_dep_pho(beg:end))
    cwf%flx_bc_dep_pho(:)= nanr
    allocate(cwf%flx_bc_dep_phi(beg:end))
    cwf%flx_bc_dep_phi(:)= nanr
    allocate(cwf%flx_bc_dep(beg:end))
    cwf%flx_bc_dep(:)= nanr
    allocate(cwf%flx_oc_dep_dry(beg:end))
    cwf%flx_oc_dep_dry(:)= nanr
    allocate(cwf%flx_oc_dep_wet(beg:end))
    cwf%flx_oc_dep_wet(:)= nanr
    allocate(cwf%flx_oc_dep_pho(beg:end))
    cwf%flx_oc_dep_pho(:)= nanr
    allocate(cwf%flx_oc_dep_phi(beg:end))
    cwf%flx_oc_dep_phi(:)= nanr
    allocate(cwf%flx_oc_dep(beg:end))
    cwf%flx_oc_dep(:)= nanr
    allocate(cwf%flx_dst_dep_dry1(beg:end))
    cwf%flx_dst_dep_dry1(:)= nanr
    allocate(cwf%flx_dst_dep_wet1(beg:end))
    cwf%flx_dst_dep_wet1(:)= nanr
    allocate(cwf%flx_dst_dep_dry2(beg:end))
    cwf%flx_dst_dep_dry2(:)= nanr
    allocate(cwf%flx_dst_dep_wet2(beg:end))
    cwf%flx_dst_dep_wet2(:)= nanr
    allocate(cwf%flx_dst_dep_dry3(beg:end))
    cwf%flx_dst_dep_dry3(:)= nanr
    allocate(cwf%flx_dst_dep_wet3(beg:end))
    cwf%flx_dst_dep_wet3(:)= nanr
    allocate(cwf%flx_dst_dep_dry4(beg:end))
    cwf%flx_dst_dep_dry4(:)= nanr
    allocate(cwf%flx_dst_dep_wet4(beg:end))
    cwf%flx_dst_dep_wet4(:)= nanr
    allocate(cwf%flx_dst_dep(beg:end))
    cwf%flx_dst_dep(:)= nanr
    allocate(cwf%qflx_snofrz_lyr(beg:end,-nlevsno+1:0))
    cwf%qflx_snofrz_lyr(:,:)=  spval
    allocate(cwf%qflx_snofrz_col(beg:end))
    cwf%qflx_snofrz_col(:)= nanr
    allocate(cwf%qflx_glcice(beg:end))
    cwf%qflx_glcice(:)= nanr
    allocate(cwf%qflx_glcice_frz(beg:end))
    cwf%qflx_glcice_frz(:)= nanr
    allocate(cwf%qflx_glcice_melt(beg:end))
    cwf%qflx_glcice_melt(:)= spval

    allocate(cwf%qflx_h2osfc_to_ice(beg:end))
    cwf%qflx_h2osfc_to_ice(:)= spval
    allocate(cwf%qflx_h2osfc_surf(beg:end))
    cwf%qflx_h2osfc_surf(:)= spval
    allocate(cwf%qflx_snow_h2osfc(beg:end))
    cwf%qflx_snow_h2osfc(:)= nanr      
    allocate(cwf%qflx_drain_perched(beg:end))
    cwf%qflx_drain_perched(:)= spval
    allocate(cwf%qflx_deficit(beg:end))
    cwf%qflx_deficit(:)= spval
    allocate(cwf%qflx_floodc(beg:end))
    cwf%qflx_floodc(:)= spval
    allocate(cwf%qflx_snow_melt(beg:end))
    cwf%qflx_snow_melt(:)= spval

  end subroutine init_column_wflux_type

  !------------------------------------------------------------------------
  subroutine init_column_cflux_type(beg, end, ccf)
    !
    ! !DESCRIPTION:
    ! Initialize column carbon flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cflux_type), intent(inout):: ccf
    !------------------------------------------------------------------------

    allocate(ccf%hrv_deadstemc_to_prod10c(beg:end))
    ccf%hrv_deadstemc_to_prod10c(:)= nanr
    allocate(ccf%hrv_deadstemc_to_prod100c(beg:end))
    ccf%hrv_deadstemc_to_prod100c(:)= nanr
    allocate(ccf%m_decomp_cpools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    ccf%m_decomp_cpools_to_fire_vr(:,:,:)= nanr
    allocate(ccf%m_decomp_cpools_to_fire(beg:end,1:ndecomp_pools))
    ccf%m_decomp_cpools_to_fire(:,:)= nanr
    allocate(ccf%decomp_cascade_hr_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    ccf%decomp_cascade_hr_vr(:,:,:)= spval
    allocate(ccf%decomp_cascade_hr(beg:end,1:ndecomp_cascade_transitions))
    ccf%decomp_cascade_hr(:,:)= nanr
    allocate(ccf%decomp_cascade_ctransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    ccf%decomp_cascade_ctransfer_vr(:,:,:)= nanr
    allocate(ccf%decomp_cascade_ctransfer(beg:end,1:ndecomp_cascade_transitions))
    ccf%decomp_cascade_ctransfer(:,:)= nanr
    allocate(ccf%decomp_cpools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    ccf%decomp_cpools_sourcesink(:,:,:)= nanr
    allocate(ccf%decomp_k(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    ccf%decomp_k(:,:,:)= spval
    allocate(ccf%t_scalar(beg:end,1:nlevdecomp_full))
    ccf%t_scalar(:,:)= spval
    allocate(ccf%w_scalar(beg:end,1:nlevdecomp_full))
    ccf%w_scalar(:,:)= spval
    allocate(ccf%hr_vr(beg:end,1:nlevdecomp_full))
    ccf%hr_vr(:,:)= nanr
    allocate(ccf%o_scalar(beg:end,1:nlevdecomp_full))
    ccf%o_scalar(:,:)= spval
    allocate(ccf%som_c_leached(beg:end))
    ccf%som_c_leached(:)= nanr
    allocate(ccf%decomp_cpools_leached(beg:end,1:ndecomp_pools))
    ccf%decomp_cpools_leached(:,:)= nanr
    allocate(ccf%decomp_cpools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    ccf%decomp_cpools_transport_tendency(:,:,:)= nanr
    allocate(ccf%phenology_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    ccf%phenology_c_to_litr_met_c(:,:)= nanr
    allocate(ccf%phenology_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    ccf%phenology_c_to_litr_cel_c(:,:)= nanr
    allocate(ccf%phenology_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    ccf%phenology_c_to_litr_lig_c(:,:)= nanr
    allocate(ccf%gap_mortality_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    ccf%gap_mortality_c_to_litr_met_c(:,:)= nanr
    allocate(ccf%gap_mortality_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    ccf%gap_mortality_c_to_litr_cel_c(:,:)= nanr
    allocate(ccf%gap_mortality_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    ccf%gap_mortality_c_to_litr_lig_c(:,:)= nanr
    allocate(ccf%gap_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full))
    ccf%gap_mortality_c_to_cwdc(:,:)= nanr
    allocate(ccf%fire_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full))
    ccf%fire_mortality_c_to_cwdc(:,:)= nanr
    allocate(ccf%m_c_to_litr_met_fire(beg:end,1:nlevdecomp_full))
    ccf%m_c_to_litr_met_fire(:,:)= nanr
    allocate(ccf%m_c_to_litr_cel_fire(beg:end,1:nlevdecomp_full))
    ccf%m_c_to_litr_cel_fire(:,:)= nanr
    allocate(ccf%m_c_to_litr_lig_fire(beg:end,1:nlevdecomp_full))
    ccf%m_c_to_litr_lig_fire(:,:)= nanr
    allocate(ccf%harvest_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    ccf%harvest_c_to_litr_met_c(:,:)= nanr
    allocate(ccf%harvest_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    ccf%harvest_c_to_litr_cel_c(:,:)= nanr
    allocate(ccf%harvest_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    ccf%harvest_c_to_litr_lig_c(:,:)= nanr
    allocate(ccf%harvest_c_to_cwdc(beg:end, 1:nlevdecomp_full))
    ccf%harvest_c_to_cwdc(:,:)= nanr
    allocate(ccf%phr_vr(beg:end,1:nlevdecomp_full))
    ccf%phr_vr(:,:)= nanr 
    allocate(ccf%somc_fire(beg:end))
    ccf%somc_fire(:)= nanr
    allocate(ccf%lf_conv_cflux(beg:end))
    ccf%lf_conv_cflux(:)= nanr
    allocate(ccf%dwt_seedc_to_leaf(beg:end))
    ccf%dwt_seedc_to_leaf(:)= nanr
    allocate(ccf%dwt_seedc_to_deadstem(beg:end))
    ccf%dwt_seedc_to_deadstem(:)= nanr
    allocate(ccf%dwt_conv_cflux(beg:end))
    ccf%dwt_conv_cflux(:)= nanr
    allocate(ccf%dwt_prod10c_gain(beg:end))
    ccf%dwt_prod10c_gain(:)= nanr
    allocate(ccf%dwt_prod100c_gain(beg:end))
    ccf%dwt_prod100c_gain(:)= nanr
    allocate(ccf%dwt_frootc_to_litr_met_c(beg:end,1:nlevdecomp_full))
    ccf%dwt_frootc_to_litr_met_c(:,:)= nanr
    allocate(ccf%dwt_frootc_to_litr_cel_c(beg:end,1:nlevdecomp_full))
    ccf%dwt_frootc_to_litr_cel_c(:,:)= nanr
    allocate(ccf%dwt_frootc_to_litr_lig_c(beg:end,1:nlevdecomp_full))
    ccf%dwt_frootc_to_litr_lig_c(:,:)= nanr
    allocate(ccf%dwt_livecrootc_to_cwdc(beg:end,1:nlevdecomp_full))
    ccf%dwt_livecrootc_to_cwdc(:,:)= nanr
    allocate(ccf%dwt_deadcrootc_to_cwdc(beg:end,1:nlevdecomp_full))
    ccf%dwt_deadcrootc_to_cwdc(:,:)= nanr
    allocate(ccf%dwt_closs(beg:end))
    ccf%dwt_closs(:)= nanr
    allocate(ccf%landuseflux(beg:end))
    ccf%landuseflux(:)= nanr
    allocate(ccf%landuptake(beg:end))
    ccf%landuptake(:)= nanr
    allocate(ccf%prod10c_loss(beg:end))
    ccf%prod10c_loss(:)= nanr
    allocate(ccf%prod100c_loss(beg:end))
    ccf%prod100c_loss(:)= nanr
    allocate(ccf%product_closs(beg:end))
    ccf%product_closs(:)= nanr

    allocate(ccf%lithr(beg:end))
    ccf%lithr(:)= nanr
    allocate(ccf%somhr(beg:end))
    ccf%somhr(:)= nanr
    allocate(ccf%hr(beg:end))
    ccf%hr(:)= nanr
    allocate(ccf%sr(beg:end))
    ccf%sr(:)= nanr
    allocate(ccf%er(beg:end))
    ccf%er(:)= nanr
    allocate(ccf%litfire(beg:end))
    ccf%litfire(:)= nanr
    allocate(ccf%somfire(beg:end))
    ccf%somfire(:)= nanr
    allocate(ccf%totfire(beg:end))
    ccf%totfire(:)= nanr
    allocate(ccf%nep(beg:end))
    ccf%nep(:)= nanr
    allocate(ccf%nbp(beg:end))
    ccf%nbp(:)= nanr
    allocate(ccf%nee(beg:end))
    ccf%nee(:)= nanr
    allocate(ccf%col_cinputs(beg:end))
    ccf%col_cinputs(:)= nanr
    allocate(ccf%col_coutputs(beg:end))
    ccf%col_coutputs(:)= nanr
    allocate(ccf%col_fire_closs(beg:end))
    ccf%col_fire_closs(:)= nanr

    allocate(ccf%cwdc_hr(beg:end))
    ccf%cwdc_hr(:)= nanr
    allocate(ccf%cwdc_loss(beg:end))
    ccf%cwdc_loss(:)= nanr
    allocate(ccf%litterc_loss(beg:end))
    ccf%litterc_loss(:)= nanr

  end subroutine init_column_cflux_type

  !------------------------------------------------------------------------
  subroutine init_column_ch4_type(beg, end, cch4)
    !
    ! !DESCRIPTION:
    ! Initialize column methane flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_ch4_type), intent(inout):: cch4
    !------------------------------------------------------------------------

    allocate(cch4%ch4_prod_depth_sat(beg:end,1:nlevgrnd))
    cch4%ch4_prod_depth_sat(:,:)= nanr
    allocate(cch4%ch4_prod_depth_unsat(beg:end,1:nlevgrnd))
    cch4%ch4_prod_depth_unsat(:,:)= nanr
    allocate(cch4%ch4_prod_depth_lake(beg:end,1:nlevgrnd))
    cch4%ch4_prod_depth_lake(:,:)= nanr
    allocate(cch4%ch4_oxid_depth_sat(beg:end,1:nlevgrnd))
    cch4%ch4_oxid_depth_sat(:,:)= nanr
    allocate(cch4%ch4_oxid_depth_unsat(beg:end,1:nlevgrnd))
    cch4%ch4_oxid_depth_unsat(:,:)= nanr
    allocate(cch4%ch4_oxid_depth_lake(beg:end,1:nlevgrnd))
    cch4%ch4_oxid_depth_lake(:,:)= nanr
    allocate(cch4%o2_oxid_depth_sat(beg:end,1:nlevgrnd))
    cch4%o2_oxid_depth_sat(:,:)= nanr
    allocate(cch4%o2_oxid_depth_unsat(beg:end,1:nlevgrnd))
    cch4%o2_oxid_depth_unsat(:,:)= nanr
    allocate(cch4%o2_decomp_depth_sat(beg:end,1:nlevgrnd))
    cch4%o2_decomp_depth_sat(:,:)= nanr
    allocate(cch4%o2_decomp_depth_unsat(beg:end,1:nlevgrnd))
    cch4%o2_decomp_depth_unsat(:,:)= spval ! To detect first time-step for denitrification cod
    allocate(cch4%o2_aere_depth_sat(beg:end,1:nlevgrnd))
    cch4%o2_aere_depth_sat(:,:)= nanr
    allocate(cch4%o2_aere_depth_unsat(beg:end,1:nlevgrnd))
    cch4%o2_aere_depth_unsat(:,:)= nanr
    allocate(cch4%co2_decomp_depth_sat(beg:end,1:nlevgrnd))
    cch4%co2_decomp_depth_sat(:,:)= nanr
    allocate(cch4%co2_decomp_depth_unsat(beg:end,1:nlevgrnd))
    cch4%co2_decomp_depth_unsat(:,:)= nanr
    allocate(cch4%co2_oxid_depth_sat(beg:end,1:nlevgrnd))
    cch4%co2_oxid_depth_sat(:,:)= nanr
    allocate(cch4%co2_oxid_depth_unsat(beg:end,1:nlevgrnd))
    cch4%co2_oxid_depth_unsat(:,:)= nanr
    allocate(cch4%ch4_aere_depth_sat(beg:end,1:nlevgrnd))
    cch4%ch4_aere_depth_sat(:,:)= nanr
    allocate(cch4%ch4_aere_depth_unsat(beg:end,1:nlevgrnd))
    cch4%ch4_aere_depth_unsat(:,:)= nanr
    allocate(cch4%ch4_tran_depth_sat(beg:end,1:nlevgrnd))
    cch4%ch4_tran_depth_sat(:,:)= nanr
    allocate(cch4%ch4_tran_depth_unsat(beg:end,1:nlevgrnd))
    cch4%ch4_tran_depth_unsat(:,:)= nanr
    allocate(cch4%co2_aere_depth_sat(beg:end,1:nlevgrnd))
    cch4%co2_aere_depth_sat(:,:)= nanr
    allocate(cch4%co2_aere_depth_unsat(beg:end,1:nlevgrnd))
    cch4%co2_aere_depth_unsat(:,:)= nanr
    allocate(cch4%ch4_surf_aere_sat(beg:end))
    cch4%ch4_surf_aere_sat(:)= nanr
    allocate(cch4%ch4_surf_aere_unsat(beg:end))
    cch4%ch4_surf_aere_unsat(:)= nanr
    allocate(cch4%ch4_ebul_depth_sat(beg:end,1:nlevgrnd))
    cch4%ch4_ebul_depth_sat(:,:)= nanr
    allocate(cch4%ch4_ebul_depth_unsat(beg:end,1:nlevgrnd))
    cch4%ch4_ebul_depth_unsat(:,:)= nanr
    allocate(cch4%ch4_ebul_total_sat(beg:end))
    cch4%ch4_ebul_total_sat(:)= nanr
    allocate(cch4%ch4_ebul_total_unsat(beg:end))
    cch4%ch4_ebul_total_unsat(:)= nanr
    allocate(cch4%ch4_surf_ebul_sat(beg:end))
    cch4%ch4_surf_ebul_sat(:)= nanr
    allocate(cch4%ch4_surf_ebul_unsat(beg:end))
    cch4%ch4_surf_ebul_unsat(:)= nanr
    allocate(cch4%ch4_surf_ebul_lake(beg:end))
    cch4%ch4_surf_ebul_lake(:)= nanr
    allocate(cch4%conc_ch4_sat(beg:end,1:nlevgrnd))
    cch4%conc_ch4_sat(:,:)= spval ! To detect file inpu
    allocate(cch4%conc_ch4_unsat(beg:end,1:nlevgrnd))
    cch4%conc_ch4_unsat(:,:)= spval ! To detect file inpu
    allocate(cch4%conc_ch4_lake(beg:end,1:nlevgrnd))
    cch4%conc_ch4_lake(:,:)= nanr ! Just a diagnostic, so nan is fin
    allocate(cch4%ch4_surf_diff_sat(beg:end))
    cch4%ch4_surf_diff_sat(:)= nanr
    allocate(cch4%ch4_surf_diff_unsat(beg:end))
    cch4%ch4_surf_diff_unsat(:)= nanr
    allocate(cch4%ch4_surf_diff_lake(beg:end))
    cch4%ch4_surf_diff_lake(:)= nanr
    allocate(cch4%conc_o2_sat(beg:end,1:nlevgrnd))
    cch4%conc_o2_sat(:,:)= spval ! To detect file inpu
    allocate(cch4%conc_o2_unsat(beg:end,1:nlevgrnd))
    cch4%conc_o2_unsat(:,:)= spval ! To detect file input and detect first time-step for denitrification cod
    allocate(cch4%conc_o2_lake(beg:end,1:nlevgrnd))
    cch4%conc_o2_lake(:,:)= nanr ! Just a diagnostic, so nan is fin
    allocate(cch4%ch4_dfsat_flux(beg:end))
    cch4%ch4_dfsat_flux(:)= nanr
    allocate(cch4%zwt_ch4_unsat(beg:end))
    cch4%zwt_ch4_unsat(:)= nanr
    allocate(cch4%fsat_bef(beg:end))
    cch4%fsat_bef(:)= spval ! To detect first time-ste
    allocate(cch4%lake_soilc(beg:end,1:nlevgrnd))
    cch4%lake_soilc(:,:)= spval ! To detect file inpu
    allocate(cch4%lake_raw(beg:end))
    cch4%lake_raw(:)= nanr
    allocate(cch4%totcolch4(beg:end))
    cch4%totcolch4(:)= spval ! To detect first time-ste
    allocate(cch4%fphr(beg:end,1:nlevgrnd))
    cch4%fphr(:,:)= nanr
    allocate(cch4%annsum_counter(beg:end))
    cch4%annsum_counter(:)= spval ! To detect first time-ste
    allocate(cch4%tempavg_somhr(beg:end))
    cch4%tempavg_somhr(:)= nanr
    allocate(cch4%annavg_somhr(beg:end))
    cch4%annavg_somhr(:)= spval ! To detect first yea
    allocate(cch4%tempavg_finrw(beg:end))
    cch4%tempavg_finrw(:)= nanr
    allocate(cch4%annavg_finrw(beg:end))
    cch4%annavg_finrw(:)= spval ! To detect first yea
    allocate(cch4%sif(beg:end))
    cch4%sif(:)= nanr
    allocate(cch4%o2stress_unsat(beg:end,1:nlevgrnd))
    cch4%o2stress_unsat(:,:)= spval ! To detect file inpu
    allocate(cch4%o2stress_sat(beg:end,1:nlevgrnd))
    cch4%o2stress_sat(:,:)= spval ! To detect file inpu
    allocate(cch4%ch4stress_unsat(beg:end,1:nlevgrnd))
    cch4%ch4stress_unsat(:,:)= nanr
    allocate(cch4%ch4stress_sat(beg:end,1:nlevgrnd))
    cch4%ch4stress_sat(:,:)= nanr    
    allocate(cch4%qflx_surf_lag(beg:end))
    cch4%qflx_surf_lag(:)= spval ! To detect file inpu
    allocate(cch4%finundated_lag(beg:end))
    cch4%finundated_lag(:)= spval ! To detect file inpu
    allocate(cch4%layer_sat_lag(beg:end,1:nlevgrnd))
    cch4%layer_sat_lag(:,:)= spval ! To detect file inpu

  end subroutine init_column_ch4_type

  !------------------------------------------------------------------------
  subroutine init_column_nflux_type(beg, end, cnf)
    !
    ! !DESCRIPTION:
    ! Initialize column nitrogen flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nflux_type), intent(inout):: cnf
    !------------------------------------------------------------------------

    allocate(cnf%ndep_to_sminn(beg:end))
    cnf%ndep_to_sminn(:)= nanr
    allocate(cnf%nfix_to_sminn(beg:end))
    cnf%nfix_to_sminn(:)= nanr
    allocate(cnf%fert_to_sminn(beg:end))
    cnf%fert_to_sminn(:)= nanr
    allocate(cnf%soyfixn_to_sminn(beg:end))
    cnf%soyfixn_to_sminn(:)= nanr
    allocate(cnf%hrv_deadstemn_to_prod10n(beg:end))
    cnf%hrv_deadstemn_to_prod10n(:)= nanr
    allocate(cnf%hrv_deadstemn_to_prod100n(beg:end))
    cnf%hrv_deadstemn_to_prod100n(:)= nanr

    allocate(cnf%m_n_to_litr_met_fire(beg:end,1:nlevdecomp_full))
    cnf%m_n_to_litr_met_fire(:,:)= nanr
    allocate(cnf%m_n_to_litr_cel_fire(beg:end,1:nlevdecomp_full))
    cnf%m_n_to_litr_cel_fire(:,:)= nanr
    allocate(cnf%m_n_to_litr_lig_fire(beg:end,1:nlevdecomp_full))
    cnf%m_n_to_litr_lig_fire(:,:)= nanr
    allocate(cnf%sminn_to_plant(beg:end))
    cnf%sminn_to_plant(:)= nanr
    allocate(cnf%potential_immob(beg:end))
    cnf%potential_immob(:)= nanr
    allocate(cnf%actual_immob(beg:end))
    cnf%actual_immob(:)= nanr
    allocate(cnf%gross_nmin(beg:end))
    cnf%gross_nmin(:)= nanr
    allocate(cnf%net_nmin(beg:end))
    cnf%net_nmin(:)= nanr
    allocate(cnf%denit(beg:end))
    cnf%denit(:)= nanr
    allocate(cnf%supplement_to_sminn(beg:end))
    cnf%supplement_to_sminn(:)= nanr
    allocate(cnf%m_decomp_npools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    cnf%m_decomp_npools_to_fire_vr(:,:,:)= nanr
    allocate(cnf%m_decomp_npools_to_fire(beg:end,1:ndecomp_pools))
    cnf%m_decomp_npools_to_fire(:,:)= nanr
    allocate(cnf%decomp_cascade_ntransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    cnf%decomp_cascade_ntransfer_vr(:,:,:)= nanr
    allocate(cnf%decomp_cascade_ntransfer(beg:end,1:ndecomp_cascade_transitions))
    cnf%decomp_cascade_ntransfer(:,:)= nanr
    allocate(cnf%decomp_cascade_sminn_flux_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    cnf%decomp_cascade_sminn_flux_vr(:,:,:)= nanr
    allocate(cnf%decomp_cascade_sminn_flux(beg:end,1:ndecomp_cascade_transitions))
    cnf%decomp_cascade_sminn_flux(:,:)= nanr
    allocate(cnf%decomp_npools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    cnf%decomp_npools_sourcesink(:,:,:)= nanr

    allocate(cnf%phenology_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    cnf%phenology_n_to_litr_met_n(:,:)= nanr
    allocate(cnf%phenology_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    cnf%phenology_n_to_litr_cel_n(:,:)= nanr
    allocate(cnf%phenology_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    cnf%phenology_n_to_litr_lig_n(:,:)= nanr
    allocate(cnf%gap_mortality_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    cnf%gap_mortality_n_to_litr_met_n(:,:)= nanr
    allocate(cnf%gap_mortality_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    cnf%gap_mortality_n_to_litr_cel_n(:,:)= nanr
    allocate(cnf%gap_mortality_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    cnf%gap_mortality_n_to_litr_lig_n(:,:)= nanr
    allocate(cnf%gap_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full))
    cnf%gap_mortality_n_to_cwdn(:,:)= nanr
    allocate(cnf%fire_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full))
    cnf%fire_mortality_n_to_cwdn(:,:)= nanr
    allocate(cnf%harvest_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    cnf%harvest_n_to_litr_met_n(:,:)= nanr
    allocate(cnf%harvest_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    cnf%harvest_n_to_litr_cel_n(:,:)= nanr
    allocate(cnf%harvest_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    cnf%harvest_n_to_litr_lig_n(:,:)= nanr
    allocate(cnf%harvest_n_to_cwdn(beg:end, 1:nlevdecomp_full))
    cnf%harvest_n_to_cwdn(:,:)= nanr

    allocate(cnf%sminn_to_denit_decomp_cascade_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    cnf%sminn_to_denit_decomp_cascade_vr(:,:,:)= nanr
    allocate(cnf%sminn_to_denit_decomp_cascade(beg:end,1:ndecomp_cascade_transitions))
    cnf%sminn_to_denit_decomp_cascade(:,:)= nanr
    allocate(cnf%sminn_to_denit_excess_vr(beg:end,1:nlevdecomp_full))
    cnf%sminn_to_denit_excess_vr(:,:)= nanr
    allocate(cnf%sminn_to_denit_excess(beg:end))
    cnf%sminn_to_denit_excess(:)= nanr
    allocate(cnf%sminn_leached_vr(beg:end,1:nlevdecomp_full))
    cnf%sminn_leached_vr(:,:)= nanr
    allocate(cnf%sminn_leached(beg:end))
    cnf%sminn_leached(:)= nanr

    allocate(cnf%f_nit_vr(beg:end,1:nlevdecomp_full))
    cnf%f_nit_vr(:,:)= nanr
    allocate(cnf%f_denit_vr(beg:end,1:nlevdecomp_full))
    cnf%f_denit_vr(:,:)= nanr
    allocate(cnf%smin_no3_leached_vr(beg:end,1:nlevdecomp_full))
    cnf%smin_no3_leached_vr(:,:)= nanr
    allocate(cnf%smin_no3_leached(beg:end))
    cnf%smin_no3_leached(:)= nanr
    allocate(cnf%smin_no3_runoff_vr(beg:end,1:nlevdecomp_full))
    cnf%smin_no3_runoff_vr(:,:)= nanr
    allocate(cnf%smin_no3_runoff(beg:end))
    cnf%smin_no3_runoff(:)= nanr
    allocate(cnf%pot_f_nit_vr(beg:end,1:nlevdecomp_full))
    cnf%pot_f_nit_vr(:,:)= nanr
    allocate(cnf%pot_f_nit(beg:end))
    cnf%pot_f_nit(:)= nanr
    allocate(cnf%pot_f_denit_vr(beg:end,1:nlevdecomp_full))
    cnf%pot_f_denit_vr(:,:)= nanr
    allocate(cnf%pot_f_denit(beg:end))
    cnf%pot_f_denit(:)= nanr
    allocate(cnf%actual_immob_no3_vr(beg:end,1:nlevdecomp_full))
    cnf%actual_immob_no3_vr(:,:)= nanr
    allocate(cnf%actual_immob_nh4_vr(beg:end,1:nlevdecomp_full))
    cnf%actual_immob_nh4_vr(:,:)= nanr
    allocate(cnf%smin_no3_to_plant_vr(beg:end,1:nlevdecomp_full))
    cnf%smin_no3_to_plant_vr(:,:)= nanr
    allocate(cnf%smin_nh4_to_plant_vr(beg:end,1:nlevdecomp_full))
    cnf%smin_nh4_to_plant_vr(:,:)= nanr
    allocate(cnf%f_nit(beg:end))
    cnf%f_nit(:)= nanr
    allocate(cnf%f_denit(beg:end))
    cnf%f_denit(:)= nanr
    allocate(cnf%n2_n2o_ratio_denit_vr(beg:end,1:nlevdecomp_full))
    cnf%n2_n2o_ratio_denit_vr(:,:)= nanr
    allocate(cnf%f_n2o_denit(beg:end))
    cnf%f_n2o_denit(:)= nanr
    allocate(cnf%f_n2o_denit_vr(beg:end,1:nlevdecomp_full))
    cnf%f_n2o_denit_vr(:,:)= nanr
    allocate(cnf%f_n2o_nit(beg:end))
    cnf%f_n2o_nit(:)= nanr
    allocate(cnf%f_n2o_nit_vr(beg:end,1:nlevdecomp_full))
    cnf%f_n2o_nit_vr(:,:)= nanr

    allocate(cnf%smin_no3_massdens_vr(beg:end,1:nlevdecomp_full))
    cnf%smin_no3_massdens_vr(:,:)= nanr
    allocate(cnf%soil_bulkdensity(beg:end,1:nlevdecomp_full))
    cnf%soil_bulkdensity(:,:)= nanr
    allocate(cnf%k_nitr_t_vr(beg:end,1:nlevdecomp_full))
    cnf%k_nitr_t_vr(:,:)= nanr
    allocate(cnf%k_nitr_ph_vr(beg:end,1:nlevdecomp_full))
    cnf%k_nitr_ph_vr(:,:)= nanr
    allocate(cnf%k_nitr_h2o_vr(beg:end,1:nlevdecomp_full))
    cnf%k_nitr_h2o_vr(:,:)= nanr
    allocate(cnf%k_nitr_vr(beg:end,1:nlevdecomp_full))
    cnf%k_nitr_vr(:,:)= nanr
    allocate(cnf%wfps_vr(beg:end,1:nlevdecomp_full))
    cnf%wfps_vr(:,:)= nanr
    allocate(cnf%fmax_denit_carbonsubstrate_vr(beg:end,1:nlevdecomp_full))
    cnf%fmax_denit_carbonsubstrate_vr(:,:)= nanr
    allocate(cnf%fmax_denit_nitrate_vr(beg:end,1:nlevdecomp_full))
    cnf%fmax_denit_nitrate_vr(:,:)= nanr
    allocate(cnf%f_denit_base_vr(beg:end,1:nlevdecomp_full))
    cnf%f_denit_base_vr(:,:)= nanr
    allocate(cnf%diffus(beg:end,1:nlevdecomp_full))
    cnf%diffus(:,:)= spval
    allocate(cnf%ratio_k1(beg:end,1:nlevdecomp_full))
    cnf%ratio_k1(:,:)= nanr
    allocate(cnf%ratio_no3_co2(beg:end,1:nlevdecomp_full))
    cnf%ratio_no3_co2(:,:)= spval
    allocate(cnf%soil_co2_prod(beg:end,1:nlevdecomp_full))
    cnf%soil_co2_prod(:,:)= nanr
    allocate(cnf%fr_WFPS(beg:end,1:nlevdecomp_full))
    cnf%fr_WFPS(:,:)= spval

    allocate(cnf%r_psi(beg:end,1:nlevdecomp_full))
    cnf%r_psi(:,:)= spval
    allocate(cnf%anaerobic_frac(beg:end,1:nlevdecomp_full))
    cnf%anaerobic_frac(:,:)= spval

    allocate(cnf%potential_immob_vr(beg:end,1:nlevdecomp_full))
    cnf%potential_immob_vr(:,:)= nanr
    allocate(cnf%actual_immob_vr(beg:end,1:nlevdecomp_full))
    cnf%actual_immob_vr(:,:)= nanr
    allocate(cnf%sminn_to_plant_vr(beg:end,1:nlevdecomp_full))
    cnf%sminn_to_plant_vr(:,:)= nanr
    allocate(cnf%supplement_to_sminn_vr(beg:end,1:nlevdecomp_full))
    cnf%supplement_to_sminn_vr(:,:)= nanr
    allocate(cnf%gross_nmin_vr(beg:end,1:nlevdecomp_full))
    cnf%gross_nmin_vr(:,:)= nanr
    allocate(cnf%net_nmin_vr(beg:end,1:nlevdecomp_full))
    cnf%net_nmin_vr(:,:)= nanr
    allocate(cnf%dwt_seedn_to_leaf(beg:end))
    cnf%dwt_seedn_to_leaf(:)= nanr
    allocate(cnf%dwt_seedn_to_deadstem(beg:end))
    cnf%dwt_seedn_to_deadstem(:)= nanr
    allocate(cnf%dwt_conv_nflux(beg:end))
    cnf%dwt_conv_nflux(:)= nanr
    allocate(cnf%dwt_prod10n_gain(beg:end))
    cnf%dwt_prod10n_gain(:)= nanr
    allocate(cnf%dwt_prod100n_gain(beg:end))
    cnf%dwt_prod100n_gain(:)= nanr
    allocate(cnf%dwt_frootn_to_litr_met_n(beg:end,1:nlevdecomp_full))
    cnf%dwt_frootn_to_litr_met_n(:,:)= nanr
    allocate(cnf%dwt_frootn_to_litr_cel_n(beg:end,1:nlevdecomp_full))
    cnf%dwt_frootn_to_litr_cel_n(:,:)= nanr
    allocate(cnf%dwt_frootn_to_litr_lig_n(beg:end,1:nlevdecomp_full))
    cnf%dwt_frootn_to_litr_lig_n(:,:)= nanr
    allocate(cnf%dwt_livecrootn_to_cwdn(beg:end,1:nlevdecomp_full))
    cnf%dwt_livecrootn_to_cwdn(:,:)= nanr
    allocate(cnf%dwt_deadcrootn_to_cwdn(beg:end,1:nlevdecomp_full))
    cnf%dwt_deadcrootn_to_cwdn(:,:)= nanr
    allocate(cnf%dwt_nloss(beg:end))
    cnf%dwt_nloss(:)= nanr
    allocate(cnf%prod10n_loss(beg:end))
    cnf%prod10n_loss(:)= nanr
    allocate(cnf%prod100n_loss(beg:end))
    cnf%prod100n_loss(:)= nanr
    allocate(cnf%product_nloss(beg:end))
    cnf%product_nloss(:)= nanr
    allocate(cnf%col_ninputs(beg:end))
    cnf%col_ninputs(:)= nanr
    allocate(cnf%col_noutputs(beg:end))
    cnf%col_noutputs(:)= nanr
    allocate(cnf%col_fire_nloss(beg:end))
    cnf%col_fire_nloss(:)= nanr
    allocate(cnf%som_n_leached(beg:end))
    cnf%som_n_leached(:)= nanr
    allocate(cnf%decomp_npools_leached(beg:end,1:ndecomp_pools))
    cnf%decomp_npools_leached(:,:)= nanr
    allocate(cnf%decomp_npools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    cnf%decomp_npools_transport_tendency(:,:,:)= nanr  

  end subroutine init_column_nflux_type

  !------------------------------------------------------------------------
  subroutine init_landunit_pstate_type(beg, end, lps)
    !
    ! !DESCRIPTION:
    ! Initialize landunit physical state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (landunit_pstate_type), intent(inout):: lps
    !------------------------------------------------------------------------

    allocate(lps%t_building(beg:end))
    lps%t_building(:)= nanr
    allocate(lps%t_building_max(beg:end))
    lps%t_building_max(:)= nanr
    allocate(lps%t_building_min(beg:end))
    lps%t_building_min(:)= nanr
    if ( nlevurb > 0 )then
       allocate(lps%tk_wall(beg:end,nlevurb))
       lps%tk_wall(:,:)= nanr
       allocate(lps%tk_roof(beg:end,nlevurb))
       lps%tk_roof(:,:)= nanr
       allocate(lps%cv_wall(beg:end,nlevurb))
       lps%cv_wall(:,:)= nanr
       allocate(lps%cv_roof(beg:end,nlevurb))
       lps%cv_roof(:,:)= nanr
    end if
    allocate(lps%tk_improad(beg:end,nlevurb))
    lps%tk_improad(:,:)= nanr
    allocate(lps%cv_improad(beg:end,nlevurb))
    lps%cv_improad(:,:)= nanr
    allocate(lps%thick_wall(beg:end))
    lps%thick_wall(:)= nanr
    allocate(lps%thick_roof(beg:end))
    lps%thick_roof(:)= nanr
    allocate(lps%nlev_improad(beg:end))
    lps%nlev_improad(:)= huge(1)
    allocate(lps%vf_sr(beg:end))
    lps%vf_sr(:)= nanr
    allocate(lps%vf_wr(beg:end))
    lps%vf_wr(:)= nanr
    allocate(lps%vf_sw(beg:end))
    lps%vf_sw(:)= nanr
    allocate(lps%vf_rw(beg:end))
    lps%vf_rw(:)= nanr
    allocate(lps%vf_ww(beg:end))
    lps%vf_ww(:)= nanr
    allocate(lps%taf(beg:end))
    lps%taf(:)= nanr
    allocate(lps%qaf(beg:end))
    lps%qaf(:)= nanr
    allocate(lps%sabs_roof_dir(beg:end,1:numrad))
    lps%sabs_roof_dir(:,:)= nanr
    allocate(lps%sabs_roof_dif(beg:end,1:numrad))
    lps%sabs_roof_dif(:,:)= nanr
    allocate(lps%sabs_sunwall_dir(beg:end,1:numrad))
    lps%sabs_sunwall_dir(:,:)= nanr
    allocate(lps%sabs_sunwall_dif(beg:end,1:numrad))
    lps%sabs_sunwall_dif(:,:)= nanr
    allocate(lps%sabs_shadewall_dir(beg:end,1:numrad))
    lps%sabs_shadewall_dir(:,:)= nanr
    allocate(lps%sabs_shadewall_dif(beg:end,1:numrad))
    lps%sabs_shadewall_dif(:,:)= nanr
    allocate(lps%sabs_improad_dir(beg:end,1:numrad))
    lps%sabs_improad_dir(:,:)= nanr
    allocate(lps%sabs_improad_dif(beg:end,1:numrad))
    lps%sabs_improad_dif(:,:)= nanr
    allocate(lps%sabs_perroad_dir(beg:end,1:numrad))
    lps%sabs_perroad_dir(:,:)= nanr
    allocate(lps%sabs_perroad_dif(beg:end,1:numrad))
    lps%sabs_perroad_dif(:,:)= nanr 

  end subroutine init_landunit_pstate_type

  !------------------------------------------------------------------------
  subroutine init_landunit_eflux_type(beg, end, lef)
    !
    ! !DESCRIPTION: 
    ! Initialize landunit energy flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end 
    type (landunit_eflux_type), intent(inout):: lef 
    !------------------------------------------------------------------------

    allocate(lef%eflx_traffic(beg:end))
    lef%eflx_traffic(:)= nanr
    allocate(lef%eflx_traffic_factor(beg:end))
    lef%eflx_traffic_factor(:)= nanr
    allocate(lef%eflx_wasteheat(beg:end))
    lef%eflx_wasteheat(:)= nanr
    allocate(lef%eflx_heat_from_ac(beg:end))
    lef%eflx_heat_from_ac(:)= nanr

  end subroutine init_landunit_eflux_type

  !------------------------------------------------------------------------
  subroutine init_gridcell_efstate_type(beg, end, gve)
    !
    ! !DESCRIPTION:
    ! Initialize gridcell isoprene emission factor variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_efstate_type), intent(inout) :: gve
    !------------------------------------------------------------------------

    allocate(gve%efisop(6,beg:end))
    gve%efisop(:,:)= nanr

  end subroutine init_gridcell_efstate_type

  !------------------------------------------------------------------------
  subroutine init_gridcell_wflux_type(beg, end, gwf)
    !
    ! !DESCRIPTION:
    ! Initialize gridcell water flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wflux_type), intent(inout):: gwf
    !------------------------------------------------------------------------

    allocate(gwf%qflx_runoffg(beg:end))
    gwf%qflx_runoffg(:)= 0._r8
    allocate(gwf%qflx_snwcp_iceg(beg:end))
    gwf%qflx_snwcp_iceg(:)= 0._r8
    allocate(gwf%qflx_liq_dynbal(beg:end))
    gwf%qflx_liq_dynbal(:)= nanr
    allocate(gwf%qflx_ice_dynbal(beg:end))
    gwf%qflx_ice_dynbal(:)= nanr

  end subroutine init_gridcell_wflux_type

  !------------------------------------------------------------------------
  subroutine init_gridcell_eflux_type(beg, end, gef)
    !
    ! !DESCRIPTION:
    ! Initialize gridcell energy flux variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_eflux_type), intent(inout):: gef
    !------------------------------------------------------------------------

    allocate(gef%eflx_sh_totg(beg:end))
    gef%eflx_sh_totg(:)= nanr
    allocate(gef%eflx_dynbal(beg:end))
    gef%eflx_dynbal(:)= nanr

  end subroutine init_gridcell_eflux_type

  !------------------------------------------------------------------------
  subroutine init_gridcell_wstate_type(beg, end, gws)
    !
    ! !DESCRIPTION:
    ! Initialize gridcell water state variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wstate_type), intent(inout):: gws
    !------------------------------------------------------------------------

    allocate(gws%gc_liq1(beg:end))
    gws%gc_liq1(:)= nanr
    allocate(gws%gc_liq2(beg:end))
    gws%gc_liq2(:)= nanr
    allocate(gws%gc_ice1(beg:end))
    gws%gc_ice1(:)= nanr
    allocate(gws%gc_ice2(beg:end))
    gws%gc_ice2(:)= nanr

  end subroutine init_gridcell_wstate_type

  !------------------------------------------------------------------------
  subroutine init_gridcell_estate_type(beg, end, ges)
    !
    ! !DESCRIPTION:
    ! Initialize gridcell energy state variables     
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_estate_type), intent(inout):: ges
    !------------------------------------------------------------------------

    allocate(ges%gc_heat1(beg:end))
    ges%gc_heat1(:)= nanr
    allocate(ges%gc_heat2(beg:end))
    ges%gc_heat2(:)= nanr

  end subroutine init_gridcell_estate_type

  !------------------------------------------------------------------------
  subroutine init_gridcell_pstate_type(beg, end, gps)
    !
    ! !DESCRIPTION:
    ! Initialize gridcell energy state variables     
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_pstate_type), intent(inout):: gps
    !------------------------------------------------------------------------

    allocate(gps%max_dayl(beg:end))
    gps%max_dayl(:)= nanr
    allocate(gps%dayl(beg:end))
    gps%dayl(:)= nanr
    allocate(gps%prev_dayl(beg:end))
    gps%prev_dayl(:)= nanr

  end subroutine init_gridcell_pstate_type


  !-------------------------------------------------------------------------
  subroutine init_gridcell_ch4_type(beg, end, gch4)
    !
    ! !DESCRIPTION:
    ! Initialize gridcell ch4 variables
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_ch4_type), intent(inout):: gch4
    !------------------------------------------------------------------------

    allocate(gch4%c_atm(beg:end,1:ngases))
    gch4%c_atm(:,:)= nanr
    allocate(gch4%ch4co2f(beg:end))
    gch4%ch4co2f(:)= nanr
    allocate(gch4%ch4prodg(beg:end))
    gch4%ch4prodg(:)= nanr
    allocate(gch4%nem(beg:end))
    gch4%nem(:)= nanr

  end subroutine init_gridcell_ch4_type

end module clmtypeInitMod
