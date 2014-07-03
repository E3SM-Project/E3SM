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
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan, bigint
  use clmtype
  use clm_varpar  , only : maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
                           numpft, ndst, nlevurb, nlevsoi, nlevdecomp, nlevdecomp_full, &
                           ndecomp_cascade_transitions, ndecomp_pools, nlevcan
#if (defined VICHYDRO)
  use clm_varpar  , only : nlayer, nlayert
#endif
  use clm_varctl  , only : use_c13, use_c14

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
!!F. Li and S. Levis (11/06/12)
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_pft_ecophys_constants
  private :: init_decomp_cascade_constants
#if (defined CNDV)
  private :: init_pft_DGVMecophys_constants
#endif
  private :: init_pft_pstate_type
  private :: init_pft_epv_type
#if (defined CNDV)
  private :: init_pft_pdgvstate_type
#endif
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
#ifdef LCH4
  private :: init_column_ch4_type
#endif
  private :: init_column_nflux_type
  private :: init_landunit_pstate_type
  private :: init_landunit_eflux_type
  private :: init_gridcell_efstate_type
  private :: init_gridcell_wflux_type
#ifdef LCH4
  private :: init_gridcell_ch4_type
#endif
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
!    *%area, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
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

    call init_pft_type     (begp, endp, clm3%g%l%c%p)
    call init_column_type  (begc, endc, clm3%g%l%c)
    call init_landunit_type(begl, endl, clm3%g%l)
    call init_gridcell_type(begg, endg, clm3%g)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    call init_decomp_cascade_constants()

#if (defined CNDV)
    ! pft DGVM-specific ecophysiological constants

    call init_pft_DGVMecophys_constants()
#endif

    ! energy balance structures (all levels)

    call init_energy_balance_type(begp, endp, clm3%g%l%c%p%pebal)
    call init_energy_balance_type(begc, endc, clm3%g%l%c%cebal)

    ! water balance structures (all levels)

    call init_water_balance_type(begp, endp, clm3%g%l%c%p%pwbal)
    call init_water_balance_type(begc, endc, clm3%g%l%c%cwbal)

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(begp, endp, clm3%g%l%c%p%pcbal)
    call init_carbon_balance_type(begc, endc, clm3%g%l%c%ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(begp, endp, clm3%g%l%c%p%pnbal)
    call init_nitrogen_balance_type(begc, endc, clm3%g%l%c%cnbal)

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(begp, endp, clm3%g%l%c%p%pps)
    call init_pft_pstate_type(begc, endc, clm3%g%l%c%cps%pps_a)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, clm3%g%l%c%p%pepv)

    !pft photosynthesis relevant variables
    call init_pft_psynstate_type(begp, endp, clm3%g%l%c%p%ppsyns)
#if (defined CNDV)
    ! pft DGVM state variables at pft level and averaged to column

    call init_pft_pdgvstate_type(begp, endp, clm3%g%l%c%p%pdgvs)
#endif
#if (defined CNDV)
    call init_pft_pdgvstate_type(begc, endc, clm3%g%l%c%cdgvs%pdgvs_a)
#endif
    call init_pft_vstate_type(begp, endp, clm3%g%l%c%p%pvs)

    ! pft energy state variables at the pft level and averaged to the column

    call init_pft_estate_type(begp, endp, clm3%g%l%c%p%pes)
    call init_pft_estate_type(begc, endc, clm3%g%l%c%ces%pes_a)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(begp, endp, clm3%g%l%c%p%pws)
    call init_pft_wstate_type(begc, endc, clm3%g%l%c%cws%pws_a)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pcs)
    call init_pft_cstate_type(begc, endc, clm3%g%l%c%ccs%pcs_a)
    
    if ( use_c13 ) then       
       call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pc13s)
       call init_pft_cstate_type(begc, endc, clm3%g%l%c%cc13s%pcs_a)
#ifdef CROP
       call endrun( trim(subname)//" ERROR:: CROP and C13 can NOT be on at the same time" )
#endif
    endif

    if ( use_c14 ) then
       call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pc14s)
       call init_pft_cstate_type(begc, endc, clm3%g%l%c%cc14s%pcs_a)
#ifdef CROP
       call endrun( trim(subname)//" ERROR:: CROP and C14 can NOT be on at the same time" )
#endif
    endif

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, clm3%g%l%c%p%pns)
    call init_pft_nstate_type(begc, endc, clm3%g%l%c%cns%pns_a)

    ! pft energy flux variables at pft level and averaged to column

    call init_pft_eflux_type(begp, endp, clm3%g%l%c%p%pef)
    call init_pft_eflux_type(begc, endc, clm3%g%l%c%cef%pef_a)

    ! pft momentum flux variables at pft level and averaged to the column

    call init_pft_mflux_type(begp, endp, clm3%g%l%c%p%pmf)
    call init_pft_mflux_type(begc, endc, clm3%g%l%c%cmf%pmf_a)

    ! pft water flux variables

    call init_pft_wflux_type(begp, endp, clm3%g%l%c%p%pwf)
    call init_pft_wflux_type(begc, endc, clm3%g%l%c%cwf%pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pcf)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%ccf%pcf_a)
    
    if ( use_c13 ) then       
       call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pc13f)
       call init_pft_cflux_type(begc, endc, clm3%g%l%c%cc13f%pcf_a)
    endif
    
    if ( use_c14 ) then
       call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pc14f)
       call init_pft_cflux_type(begc, endc, clm3%g%l%c%cc14f%pcf_a)
    endif

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, clm3%g%l%c%p%pnf)
    call init_pft_nflux_type(begc, endc, clm3%g%l%c%cnf%pnf_a)

    ! pft VOC flux variables at pft level

    call init_pft_vflux_type(begp, endp, clm3%g%l%c%p%pvf)

    ! gridcell VOC emission factors (heald, 05/06)

    call init_gridcell_efstate_type(begg, endg, clm3%g%gve)

    ! pft dust flux variables at pft level and averaged to column

    call init_pft_dflux_type(begp, endp, clm3%g%l%c%p%pdf)
    call init_pft_dflux_type(begc, endc, clm3%g%l%c%cdf%pdf_a)

    ! pft dry dep velocity variables at pft level and averaged to column

    call init_pft_depvd_type(begp, endp, clm3%g%l%c%p%pdd)

    ! column physical state variables at column level and averaged to
    ! the landunit

    call init_column_pstate_type(begc, endc, clm3%g%l%c%cps)
    call init_column_pstate_type(begl, endl, clm3%g%l%lps%cps_a)

    ! column energy state variables at column level and averaged to
    ! the gridcell

    call init_column_estate_type(begc, endc, clm3%g%l%c%ces)
    call init_column_estate_type(begg, endg, clm3%g%ges%ces_a)

    ! column water state variables at column level and averaged to
    ! the gridcell

    call init_column_wstate_type(begc, endc, clm3%g%l%c%cws)
    call init_column_wstate_type(begg, endg, clm3%g%gws%cws_a)

    ! column carbon state variables at column level

    call init_column_cstate_type(begc, endc, clm3%g%l%c%ccs)
    
    if ( use_c13 ) then       
       call init_column_cstate_type(begc, endc, clm3%g%l%c%cc13s)
    endif

    if ( use_c14 ) then       
       call init_column_cstate_type(begc, endc, clm3%g%l%c%cc14s)
    endif
    ! column nitrogen state variables at column level

    ! column nitrogen state variables at column level

    call init_column_nstate_type(begc, endc, clm3%g%l%c%cns)

    ! column energy flux variables at column level and averaged to
    ! the landunit and gridcell

    call init_column_eflux_type(begc, endc, clm3%g%l%c%cef)
    call init_column_eflux_type(begl, endl, clm3%g%l%lef%cef_a)
    call init_column_eflux_type(begg, endg, clm3%g%gef%cef_a)

    ! column water flux variables at column level and averaged to
    ! gridcell

    call init_column_wflux_type(begc, endc, clm3%g%l%c%cwf)
    call init_column_wflux_type(begg, endg, clm3%g%gwf%cwf_a)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, clm3%g%l%c%ccf)
    
    if ( use_c13 ) then       
       call init_column_cflux_type(begc, endc, clm3%g%l%c%cc13f)
    endif
    
    if ( use_c14 ) then       
       call init_column_cflux_type(begc, endc, clm3%g%l%c%cc14f)
    endif

#if (defined LCH4)
    ! column CH4 flux variables at column level
    call init_column_ch4_type(begc, endc, clm3%g%l%c%cch4)
#endif

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, clm3%g%l%c%cnf)

    ! land unit physical state variables

    call init_landunit_pstate_type(begl, endl, clm3%g%l%lps)

    ! land unit energy flux variables 

    call init_landunit_eflux_type(begl, endl, clm3%g%l%lef)

#if (defined CNDV)
    ! gridcell DGVM variables

    call init_gridcell_dgvstate_type(begg, endg, clm3%g%gdgvs)
#endif

    ! gridcell physical state variables


    ! gridcell: water flux variables

    call init_gridcell_wflux_type(begg, endg, clm3%g%gwf)

    ! gridcell: energy flux variables

    call init_gridcell_eflux_type(begg, endg, clm3%g%gef)

    ! gridcell: water state variables

    call init_gridcell_wstate_type(begg, endg, clm3%g%gws)

    ! gridcell: energy state variables

    call init_gridcell_estate_type(begg, endg, clm3%g%ges)

#if (defined LCH4)
    ! gridcell: ch4 variables

    call init_gridcell_ch4_type(begg, endg, clm3%g%gch4)
#endif

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

    allocate(p%gridcell(beg:end),p%wtgcell(beg:end))
    allocate(p%landunit(beg:end),p%wtlunit(beg:end))
    allocate(p%column  (beg:end),p%wtcol  (beg:end))

    allocate(p%itype(beg:end))
    allocate(p%mxy(beg:end))
    allocate(p%active(beg:end))

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

   allocate(c%gridcell(beg:end),c%wtgcell(beg:end))
   allocate(c%landunit(beg:end),c%wtlunit(beg:end))

   allocate(c%pfti(beg:end),c%pftf(beg:end),c%npfts(beg:end))

   allocate(c%itype(beg:end))
   allocate(c%active(beg:end))

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

   allocate(l%gridcell(beg:end),l%wtgcell(beg:end))

   allocate(l%coli(beg:end),l%colf(beg:end),l%ncolumns(beg:end))
   allocate(l%pfti(beg:end),l%pftf(beg:end),l%npfts   (beg:end))

   allocate(l%itype(beg:end))
   allocate(l%ifspecial(beg:end))
   allocate(l%lakpoi(beg:end))
   allocate(l%urbpoi(beg:end))
   allocate(l%glcmecpoi(beg:end))
   allocate(l%udenstype(beg:end))
   allocate(l%active(beg:end))

   ! MV - these should be moved to landunit physical state -MV
   allocate(l%canyon_hwr(beg:end))
   allocate(l%wtroad_perv(beg:end))
   allocate(l%ht_roof(beg:end))
   allocate(l%wtlunit_roof(beg:end))
   allocate(l%wind_hgt_canyon(beg:end))
   allocate(l%z_0_town(beg:end))
   allocate(l%z_d_town(beg:end))

   l%canyon_hwr(beg:end)  = nan
   l%wtroad_perv(beg:end) = nan
   l%ht_roof(beg:end) = nan
   l%wtlunit_roof(beg:end) = nan
   l%wind_hgt_canyon(beg:end) = nan
   l%z_0_town(beg:end) = nan
   l%z_d_town(beg:end) = nan

   l%glcmecpoi(beg:end) = .false.

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

   allocate(g%luni(beg:end),g%lunf(beg:end),g%nlandunits(beg:end))
   allocate(g%coli(beg:end),g%colf(beg:end),g%ncolumns  (beg:end))
   allocate(g%pfti(beg:end),g%pftf(beg:end),g%npfts     (beg:end))

   allocate(g%gindex(beg:end))
   allocate(g%area(beg:end))
   allocate(g%lat(beg:end))
   allocate(g%lon(beg:end))
   allocate(g%latdeg(beg:end))
   allocate(g%londeg(beg:end))
   allocate(g%gindex_a(beg:end))
   allocate(g%lat_a(beg:end))
   allocate(g%lon_a(beg:end))
   allocate(g%latdeg_a(beg:end))
   allocate(g%londeg_a(beg:end))

   allocate(g%gris_mask(beg:end))
   allocate(g%gris_area(beg:end))
   allocate(g%aais_mask(beg:end))
   allocate(g%aais_area(beg:end))
   allocate(g%tws(beg:end))
   g%gris_mask(beg:end) = nan
   g%gris_area(beg:end) = nan
   g%aais_mask(beg:end) = nan
   g%aais_area(beg:end) = nan
   g%tws(beg:end) = nan

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
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft,numrad))
    allocate(pftcon%rhos(0:numpft,numrad))
    allocate(pftcon%taul(0:numpft,numrad))
    allocate(pftcon%taus(0:numpft,numrad))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
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
    allocate(pftcon%dwood(0:numpft))
    allocate(pftcon%rootprof_beta(0:numpft))
    allocate(pftcon%fertnitro(0:numpft))
    allocate(pftcon%fleafcn(0:numpft))
    allocate(pftcon%ffrootcn(0:numpft))
    allocate(pftcon%fstemcn(0:numpft))

    pftcon%noveg(:) = bigint
    pftcon%tree(:) = bigint
    pftcon%smpso(:) = nan
    pftcon%smpsc(:) = nan
    pftcon%fnitr(:) = nan
    pftcon%foln(:) = nan
    pftcon%dleaf(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%xl(:) = nan
    pftcon%rhol(:,:numrad) = nan
    pftcon%rhos(:,:numrad) = nan
    pftcon%taul(:,:numrad) = nan
    pftcon%taus(:,:numrad) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%roota_par(:) = nan
    pftcon%rootb_par(:) = nan
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
    pftcon%dwood(:) = nan
    pftcon%rootprof_beta(:) = nan
    pftcon%fertnitro(:) = nan
    pftcon%fleafcn(:)   = nan
    pftcon%ffrootcn(:)  = nan
    pftcon%fstemcn(:)   = nan
  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_decomp_cascade_constants
!
! !INTERFACE:
  subroutine init_decomp_cascade_constants()
!
! !DESCRIPTION:
! Initialize decomposition cascade state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Charlie Koven
!
!EOP
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


#if (defined CNDV)
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
#endif

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
    use clm_varcon, only : spval,ispval
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
    allocate(pps%prec10(beg:end)) !F. Li and S. Levis
    allocate(pps%prec60(beg:end)) !F. Li and S. Levis
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
    allocate(pps%rhal(beg:end))
    allocate(pps%vpdal(beg:end))
    allocate(pps%rssun_z(beg:end,1:nlevcan))
    allocate(pps%rssha_z(beg:end,1:nlevcan))
    allocate(pps%laisun(beg:end))
    allocate(pps%laisha(beg:end))
    allocate(pps%laisun_z(beg:end,1:nlevcan))
    allocate(pps%laisha_z(beg:end,1:nlevcan))
    allocate(pps%btran(beg:end))
    allocate(pps%btran2(beg:end))   ! F. Li and S. Levis
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
    allocate(pps%fabd_sun(beg:end,1:numrad))
    allocate(pps%fabd_sha(beg:end,1:numrad))
    allocate(pps%fabi(beg:end,1:numrad))
    allocate(pps%fabi_sun(beg:end,1:numrad))
    allocate(pps%fabi_sha(beg:end,1:numrad))
    allocate(pps%ftdd(beg:end,1:numrad))
    allocate(pps%ftid(beg:end,1:numrad))
    allocate(pps%ftii(beg:end,1:numrad))
    allocate(pps%vcmaxcintsun(beg:end))
    allocate(pps%vcmaxcintsha(beg:end))
    allocate(pps%ncan(beg:end))
    allocate(pps%nrad(beg:end))
    allocate(pps%fabd_sun_z(beg:end,1:nlevcan))
    allocate(pps%fabd_sha_z(beg:end,1:nlevcan))
    allocate(pps%fabi_sun_z(beg:end,1:nlevcan))
    allocate(pps%fabi_sha_z(beg:end,1:nlevcan))
    allocate(pps%fsun_z(beg:end,1:nlevcan))
    allocate(pps%tlai_z(beg:end,1:nlevcan))
    allocate(pps%tsai_z(beg:end,1:nlevcan))
    allocate(pps%u10(beg:end))
    allocate(pps%u10_clm(beg:end))
    allocate(pps%va(beg:end))
    allocate(pps%fv(beg:end))
    allocate(pps%ram1(beg:end))
    allocate(pps%burndate(beg:end))   ! F. Li and S. Levis
    allocate(pps%ram1_lake(beg:end))
    allocate(pps%rh_leaf(beg:end))
    allocate(pps%rhaf(beg:end))
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
    allocate(pps%forc_hgt_u_pft(beg:end))
    allocate(pps%forc_hgt_t_pft(beg:end))
    allocate(pps%forc_hgt_q_pft(beg:end))
    allocate(pps%lfpftd(beg:end))      !F. Li and S. Levis

    ! 4/14/05: PET
    ! Adding isotope code
    
    if ( use_c13 ) then       
       allocate(pps%alphapsnsun(beg:end))
       allocate(pps%alphapsnsha(beg:end))
    endif

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
    
#if (defined LCH4)
    ! CH4 code
    allocate(pps%grnd_ch4_cond(beg:end))
    allocate(pps%canopy_cond(beg:end))
#endif
   ! and vertical profiles for calculating fluxes
    allocate(pps%leaf_prof(beg:end,1:nlevdecomp_full))
    allocate(pps%froot_prof(beg:end,1:nlevdecomp_full))
    allocate(pps%croot_prof(beg:end,1:nlevdecomp_full))
    allocate(pps%stem_prof(beg:end,1:nlevdecomp_full))
    pps%prec10(beg:end) = nan   ! F. Li and S. Levis
    pps%prec60(beg:end) = nan   ! F. Li and S. Levis
    pps%frac_veg_nosno(beg:end) = bigint
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
    pps%rhal(beg:end) = nan
    pps%vpdal(beg:end) = nan    
    pps%rssha(beg:end) = nan
    pps%rssun_z(beg:end,:nlevcan) = nan
    pps%rssha_z(beg:end,:nlevcan) = nan
    pps%laisun(beg:end) = nan
    pps%laisha(beg:end) = nan
    pps%laisun_z(beg:end,:nlevcan) = nan
    pps%laisha_z(beg:end,:nlevcan) = nan
    pps%btran(beg:end) = spval
    pps%btran2(beg:end) = spval       !F. Li and S. Levis
    pps%fsun(beg:end) = spval
    pps%fsun_z(beg:end,:nlevcan) = 0._r8
    pps%tlai(beg:end) = 0._r8
    pps%tsai(beg:end) = 0._r8
    pps%elai(beg:end) = 0._r8
    pps%tlai_z(beg:end,:nlevcan) = 0._r8
    pps%tsai_z(beg:end,:nlevcan) = 0._r8
    pps%esai(beg:end) = 0._r8
    pps%ncan(beg:end) = 0
    pps%nrad(beg:end) = 0
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
    pps%fabd_sun(beg:end,:numrad) = nan
    pps%fabd_sha(beg:end,:numrad) = nan
    pps%fabi(beg:end,:numrad) = nan
    pps%fabi_sun(beg:end,:numrad) = nan
    pps%fabi_sha(beg:end,:numrad) = nan
    pps%ftdd(beg:end,:numrad) = nan
    pps%ftid(beg:end,:numrad) = nan
    pps%ftii(beg:end,:numrad) = nan
    pps%vcmaxcintsun(beg:end) = nan
    pps%vcmaxcintsha(beg:end) = nan
    pps%fabd_sun_z(beg:end,:nlevcan) = 0._r8
    pps%fabd_sha_z(beg:end,:nlevcan) = 0._r8
    pps%fabi_sun_z(beg:end,:nlevcan) = 0._r8
    pps%fabi_sha_z(beg:end,:nlevcan) = 0._r8
    pps%u10(beg:end) = nan
    pps%u10_clm(beg:end) = nan
    pps%va(beg:end) = nan
    pps%fv(beg:end) = nan
    pps%ram1(beg:end) = nan
    pps%burndate(beg:end)    = ispval   ! F. Li and S. Levis
    pps%ram1_lake(beg:end) = nan
    pps%rh_leaf(beg:end) = spval
    pps%rhaf(beg:end)    = spval
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
       pps%harvdate(beg:end)    = bigint
       pps%idop(beg:end)        = bigint
       pps%peaklai(beg:end)     = 0
    end if
    pps%vds(beg:end) = nan
    pps%forc_hgt_u_pft(beg:end) = nan
    pps%forc_hgt_t_pft(beg:end) = nan
    pps%forc_hgt_q_pft(beg:end) = nan
    ! 4/14/05: PET
    ! Adding isotope code    ! EBK Check this!
    !!!pps%cisun(beg:end) = spval
    !!!pps%cisha(beg:end) = spval
    
    if ( use_c13 ) then       
       pps%alphapsnsun(beg:end) = spval
       pps%alphapsnsha(beg:end) = spval
    endif

#if defined (LCH4)
    ! CH4 code
    pps%grnd_ch4_cond(beg:end) = nan
    pps%canopy_cond(beg:end) = nan
#endif
   ! and vertical profiles for calculating fluxes
    pps%leaf_prof(beg:end,1:nlevdecomp_full) = spval
    pps%froot_prof(beg:end,1:nlevdecomp_full) = spval
    pps%croot_prof(beg:end,1:nlevdecomp_full) = spval
    pps%stem_prof(beg:end,1:nlevdecomp_full) = spval

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
    use clm_varcon, only : spval
!
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
    allocate(pepv%fert_counter(beg:end))
    allocate(pepv%grain_flag(beg:end))
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
    if ( use_c13 ) then
       allocate(pepv%xsmrpool_c13ratio(beg:end))
    endif
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
#if (defined CNDV)
    allocate(pepv%tempsum_litfall(beg:end))
    allocate(pepv%annsum_litfall(beg:end))
#endif
    if ( use_c13 ) then
       allocate(pepv%rc13_canair(beg:end))
       allocate(pepv%rc13_psnsun(beg:end))
       allocate(pepv%rc13_psnsha(beg:end))
    endif
    
    if ( use_c14 ) then
       allocate(pepv%rc14_atm(beg:end))
       ! allocate(pepv%rc14_canair(beg:end))
       ! allocate(pepv%rc14_psnsun(beg:end))
       ! allocate(pepv%rc14_psnsha(beg:end))
    endif

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
    pepv%fert_counter(beg:end) = nan
    pepv%grain_flag(beg:end) = nan
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
    if ( use_c13 ) then
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
#if (defined CNDV)
    pepv%tempsum_litfall(beg:end) = nan
    pepv%annsum_litfall(beg:end) = nan
#endif
    if ( use_c13 ) then
       pepv%rc13_canair(beg:end) = spval
       pepv%rc13_psnsun(beg:end) = spval
       pepv%rc13_psnsha(beg:end) = spval
    endif

    if ( use_c14 ) then
       pepv%rc14_atm(beg:end) = nan
       ! pepv%rc14_canair(beg:end) = nan
       ! pepv%rc14_psnsun(beg:end) = nan
       ! pepv%rc14_psnsha(beg:end) = nan
    endif
    
  end subroutine init_pft_epv_type

#if (defined CNDV)
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
#endif

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
! !IROUTINE: init_pft_psynstate_type
!
! !INTERFACE:
  subroutine init_pft_psynstate_type(beg, end, ppsyns)
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
    type (pft_psynstate_type), intent(inout):: ppsyns
!
! !REVISION HISTORY:
! Created by Jinyun Tang
!
!EOP
!-----------------------------------------------------------------------

   allocate(ppsyns%c3flag(beg:end))
   allocate(ppsyns%ac(beg:end,1:nlevcan))
   allocate(ppsyns%aj(beg:end,1:nlevcan))
   allocate(ppsyns%ap(beg:end,1:nlevcan))
   allocate(ppsyns%ag(beg:end,1:nlevcan))
   allocate(ppsyns%an(beg:end,1:nlevcan))
   allocate(ppsyns%vcmax_z(beg:end,1:nlevcan))
   allocate(ppsyns%cp(beg:end))
   allocate(ppsyns%kc(beg:end))
   allocate(ppsyns%ko(beg:end))
   allocate(ppsyns%qe(beg:end))
   allocate(ppsyns%tpu_z(beg:end,1:nlevcan))
   allocate(ppsyns%kp_z(beg:end,1:nlevcan))   
   allocate(ppsyns%theta_cj(beg:end))
   allocate(ppsyns%bbb(beg:end))
   allocate(ppsyns%mbb(beg:end))
   allocate(ppsyns%gb_mol(beg:end))
   allocate(ppsyns%gs_mol(beg:end,1:nlevcan))

   ppsyns%c3flag = .false.
   ppsyns%ac(beg:end,1:nlevcan) = nan
   ppsyns%aj(beg:end,1:nlevcan) = nan
   ppsyns%ap(beg:end,1:nlevcan) = nan
   ppsyns%ag(beg:end,1:nlevcan) = nan
   ppsyns%an(beg:end,1:nlevcan) = nan
   ppsyns%vcmax_z(beg:end,1:nlevcan) = nan
   ppsyns%cp(beg:end) = nan
   ppsyns%kc(beg:end) = nan
   ppsyns%ko(beg:end) = nan
   ppsyns%qe(beg:end) = nan
   ppsyns%tpu_z(beg:end,1:nlevcan) = nan
   ppsyns%kp_z(beg:end,1:nlevcan) = nan
   ppsyns%theta_cj(beg:end) = nan
   ppsyns%bbb(beg:end) = nan
   ppsyns%mbb(beg:end) = nan
   ppsyns%gb_mol(beg:end) = nan
   ppsyns%gs_mol(beg:end,1:nlevcan) = nan

  end subroutine init_pft_psynstate_type


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
! !USES:
    use clm_varcon, only : spval
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
    pws%h2ocan(beg:end) = spval

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
#ifdef CN
    allocate(pcs%woodc(beg:end))
#endif

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
#ifdef CN
    pcs%woodc(beg:end) = nan
#endif

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
    allocate(pef%parsun_z(beg:end,1:nlevcan))
    allocate(pef%parsha_z(beg:end,1:nlevcan))
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
    allocate(pef%eflx_grnd_lake(beg:end))
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
    allocate(pef%fsds_vis_i_ln(beg:end))
    allocate(pef%parveg_ln(beg:end))
    allocate(pef%fsds_nir_d_ln(beg:end))
    allocate(pef%fsr_vis_d_ln(beg:end))
    allocate(pef%fsr_nir_d_ln(beg:end))
    allocate(pef%sabg_lyr(beg:end,-nlevsno+1:1))
    allocate(pef%sabg_pen(beg:end))
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
    pef%parsun_z(beg:end,:nlevcan) = nan
    pef%parsha_z(beg:end,:nlevcan) = nan
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
    pef%eflx_grnd_lake(beg:end) = nan
    pef%dgnetdT(beg:end) = nan
    pef%eflx_lwrad_out(beg:end) = nan
    pef%eflx_lwrad_net(beg:end) = nan
    pef%eflx_lwrad_net_u(beg:end) = nan
    pef%eflx_lwrad_net_r(beg:end) = nan
    pef%netrad(beg:end) = nan
    pef%fsds_vis_d(beg:end) = nan
    pef%fsds_vis_i_ln(beg:end) = nan
    pef%parveg_ln(beg:end) = nan
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
    pef%sabg_lyr(beg:end,-nlevsno+1:1) = nan
    pef%sabg_pen(beg:end) = nan
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
    allocate(pef%eflx_sh_snow(beg:end))
    allocate(pef%eflx_sh_soil(beg:end))
    allocate(pef%eflx_sh_h2osfc(beg:end))
    pef%eflx_sh_snow(beg:end) = nan
    pef%eflx_sh_soil(beg:end) = nan
    pef%eflx_sh_h2osfc(beg:end) = nan

    allocate(pef%sabg_soil(beg:end))
    allocate(pef%sabg_snow(beg:end))
    allocate(pef%sabg_chk(beg:end))
    pef%sabg_soil(beg:end) = nan
    pef%sabg_snow(beg:end) = nan
    pef%sabg_chk(beg:end) = nan
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
    pwf%qflx_evap_grnd(beg:end) = 0.0_r8
    pwf%qflx_dew_grnd(beg:end) = 0.0_r8
    pwf%qflx_sub_snow(beg:end) = 0.0_r8
    pwf%qflx_dew_snow(beg:end) = 0.0_r8

    allocate(pwf%qflx_ev_snow(beg:end))
    allocate(pwf%qflx_ev_soil(beg:end))
    allocate(pwf%qflx_ev_h2osfc(beg:end))
    pwf%qflx_ev_snow(beg:end) = nan
    pwf%qflx_ev_soil(beg:end) = nan
    pwf%qflx_ev_h2osfc(beg:end) = nan
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
    allocate(pcf%psnsun_z(beg:end,1:nlevcan))
    allocate(pcf%psnsha_z(beg:end,1:nlevcan))
    allocate(pcf%cisun_z(beg:end,1:nlevcan))
    allocate(pcf%cisha_z(beg:end,1:nlevcan))
    allocate(pcf%lmrsun(beg:end))
    allocate(pcf%lmrsha(beg:end))
    allocate(pcf%lmrsun_z(beg:end,1:nlevcan))
    allocate(pcf%lmrsha_z(beg:end,1:nlevcan))
    allocate(pcf%fpsn(beg:end))
    allocate(pcf%fco2(beg:end))
    allocate(pcf%psnsun_wc(beg:end))
    allocate(pcf%psnsha_wc(beg:end))
    allocate(pcf%fpsn_wc(beg:end))
    allocate(pcf%psnsun_wj(beg:end))
    allocate(pcf%psnsha_wj(beg:end))
    allocate(pcf%fpsn_wj(beg:end))
    allocate(pcf%psnsun_wp(beg:end))
    allocate(pcf%psnsha_wp(beg:end))
    allocate(pcf%fpsn_wp(beg:end))

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
             
    ! fire related variables changed by F. Li and S. Levis           
    allocate(pcf%m_leafc_to_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_to_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_fire(beg:end))
    allocate(pcf%m_frootc_to_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livecrootc_to_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_fire(beg:end))
    allocate(pcf%m_leafc_to_litter_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_to_deadstemc_fire(beg:end))
    allocate(pcf%m_deadstemc_to_litter_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_frootc_to_litter_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_deadcrootc_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_litter_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_litter_fire(beg:end))


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
    allocate(pcf%grain_mr(beg:end))
    allocate(pcf%leaf_curmr(beg:end))
    allocate(pcf%froot_curmr(beg:end))
    allocate(pcf%livestem_curmr(beg:end))
    allocate(pcf%livecroot_curmr(beg:end))
    allocate(pcf%grain_curmr(beg:end))
    allocate(pcf%leaf_xsmr(beg:end))
    allocate(pcf%froot_xsmr(beg:end))
    allocate(pcf%livestem_xsmr(beg:end))
    allocate(pcf%livecroot_xsmr(beg:end))
    allocate(pcf%grain_xsmr(beg:end))
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
#ifdef CN
    allocate(pcf%frootc_alloc(beg:end))
    allocate(pcf%frootc_loss(beg:end))
    allocate(pcf%leafc_alloc(beg:end))
    allocate(pcf%leafc_loss(beg:end))
    allocate(pcf%woodc_alloc(beg:end))
    allocate(pcf%woodc_loss(beg:end))
#endif
#ifdef LCH4
    allocate(pcf%tempavg_agnpp(beg:end))
    allocate(pcf%tempavg_bgnpp(beg:end))
    allocate(pcf%annavg_agnpp(beg:end))
    allocate(pcf%annavg_bgnpp(beg:end))
#endif

    pcf%psnsun(beg:end) = nan
    pcf%psnsha(beg:end) = nan
    pcf%psnsun_z(beg:end,:nlevcan) = nan
    pcf%psnsha_z(beg:end,:nlevcan) = nan
    pcf%cisun_z(beg:end,:nlevcan) = nan
    pcf%cisha_z(beg:end,:nlevcan) = nan
    pcf%lmrsun(beg:end) = nan
    pcf%lmrsha(beg:end) = nan
    pcf%lmrsun_z(beg:end,:nlevcan) = nan
    pcf%lmrsha_z(beg:end,:nlevcan) = nan
    pcf%fpsn(beg:end) = spval
    pcf%fco2(beg:end) = 0._r8
    pcf%psnsun_wc(beg:end) = nan
    pcf%psnsha_wc(beg:end) = nan
    pcf%fpsn_wc(beg:end) = nan
    pcf%psnsun_wj(beg:end) = nan
    pcf%psnsha_wj(beg:end) = nan
    pcf%fpsn_wj(beg:end) = nan
    pcf%psnsun_wp(beg:end) = nan
    pcf%psnsha_wp(beg:end) = nan
    pcf%fpsn_wp(beg:end) = nan

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
     
    ! fire variable changed by F. Li and S. Levis
    pcf%m_leafc_to_fire(beg:end) = nan
    pcf%m_leafc_storage_to_fire(beg:end) = nan
    pcf%m_leafc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_to_fire(beg:end) = nan
    pcf%m_livestemc_storage_to_fire(beg:end) = nan
    pcf%m_livestemc_xfer_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_storage_to_fire(beg:end) = nan
    pcf%m_deadstemc_xfer_to_fire(beg:end) = nan
    pcf%m_frootc_to_fire(beg:end) = nan
    pcf%m_frootc_storage_to_fire(beg:end) = nan
    pcf%m_frootc_xfer_to_fire(beg:end) = nan
    pcf%m_livecrootc_to_fire(beg:end) = nan
    pcf%m_livecrootc_storage_to_fire(beg:end) = nan
    pcf%m_livecrootc_xfer_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_storage_to_fire(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_fire(beg:end) = nan
    pcf%m_gresp_storage_to_fire(beg:end) = nan
    pcf%m_gresp_xfer_to_fire(beg:end) = nan
    
    pcf%m_leafc_to_litter_fire(beg:end) = nan
    pcf%m_leafc_storage_to_litter_fire(beg:end) = nan
    pcf%m_leafc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_storage_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_to_deadstemc_fire(beg:end) = nan
    pcf%m_deadstemc_to_litter_fire(beg:end) = nan
    pcf%m_deadstemc_storage_to_litter_fire(beg:end) = nan
    pcf%m_deadstemc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_frootc_to_litter_fire(beg:end) = nan
    pcf%m_frootc_storage_to_litter_fire(beg:end) = nan
    pcf%m_frootc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_storage_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_to_deadcrootc_fire(beg:end) = nan
    pcf%m_deadcrootc_to_litter_fire(beg:end) = nan
    pcf%m_deadcrootc_storage_to_litter_fire(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_gresp_storage_to_litter_fire(beg:end) = nan
    pcf%m_gresp_xfer_to_litter_fire(beg:end) = nan


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
    pcf%grain_mr(beg:end) = nan
    pcf%leaf_curmr(beg:end) = nan
    pcf%froot_curmr(beg:end) = nan
    pcf%livestem_curmr(beg:end) = nan
    pcf%livecroot_curmr(beg:end) = nan
    pcf%grain_curmr(beg:end) = nan
    pcf%leaf_xsmr(beg:end) = nan
    pcf%froot_xsmr(beg:end) = nan
    pcf%livestem_xsmr(beg:end) = nan
    pcf%livecroot_xsmr(beg:end) = nan
    pcf%grain_xsmr(beg:end) = nan
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
#if (defined CN)
    pcf%frootc_alloc(beg:end) = nan
    pcf%frootc_loss(beg:end) = nan
    pcf%leafc_alloc(beg:end) = nan
    pcf%leafc_loss(beg:end) = nan
    pcf%woodc_alloc(beg:end) = nan
    pcf%woodc_loss(beg:end) = nan
#endif
#if (defined LCH4)
    pcf%tempavg_agnpp(beg:end) = spval ! For back-compatibility
    pcf%tempavg_bgnpp(beg:end) = spval ! For back-compatibility
    pcf%annavg_agnpp(beg:end) = spval ! To detect first year
    pcf%annavg_bgnpp(beg:end) = spval ! To detect first year
#endif


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
           
    ! fire variables changed by F. Li and S. Levis            
   allocate(pnf%m_leafn_to_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_to_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_fire(beg:end))
    allocate(pnf%m_frootn_to_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livecrootn_to_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_retransn_to_fire(beg:end))

    allocate(pnf%m_leafn_to_litter_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_to_deadstemn_fire(beg:end))
    allocate(pnf%m_deadstemn_to_litter_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_frootn_to_litter_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_deadcrootn_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_litter_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_retransn_to_litter_fire(beg:end))



    allocate(pnf%leafn_xfer_to_leafn(beg:end))
    allocate(pnf%frootn_xfer_to_frootn(beg:end))
    allocate(pnf%livestemn_xfer_to_livestemn(beg:end))
    allocate(pnf%deadstemn_xfer_to_deadstemn(beg:end))
    allocate(pnf%livecrootn_xfer_to_livecrootn(beg:end))
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(beg:end))
    allocate(pnf%leafn_to_litter(beg:end))
    allocate(pnf%leafn_to_retransn(beg:end))
    allocate(pnf%frootn_to_retransn(beg:end))
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
       allocate(pnf%fert(beg:end))
       allocate(pnf%soyfixn(beg:end))
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
        
    ! fire varibles changed by F. Li and S. Levis         
   pnf%m_leafn_to_fire(beg:end) = nan
    pnf%m_leafn_storage_to_fire(beg:end) = nan
    pnf%m_leafn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_to_fire(beg:end) = nan
    pnf%m_livestemn_storage_to_fire(beg:end) = nan
    pnf%m_livestemn_xfer_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_storage_to_fire(beg:end) = nan
    pnf%m_deadstemn_xfer_to_fire(beg:end) = nan
    pnf%m_frootn_to_fire(beg:end) = nan
    pnf%m_frootn_storage_to_fire(beg:end) = nan
    pnf%m_frootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_to_fire(beg:end) = nan
    pnf%m_livecrootn_storage_to_fire(beg:end) = nan
    pnf%m_livecrootn_xfer_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_storage_to_fire(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_fire(beg:end) = nan
    pnf%m_retransn_to_fire(beg:end) = nan

    pnf%m_leafn_to_litter_fire(beg:end) = nan
    pnf%m_leafn_storage_to_litter_fire(beg:end) = nan
    pnf%m_leafn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_storage_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_to_deadstemn_fire(beg:end) = nan
    pnf%m_deadstemn_to_litter_fire(beg:end) = nan
    pnf%m_deadstemn_storage_to_litter_fire(beg:end) = nan
    pnf%m_deadstemn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_frootn_to_litter_fire(beg:end) = nan
    pnf%m_frootn_storage_to_litter_fire(beg:end) = nan
    pnf%m_frootn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_storage_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_to_deadcrootn_fire(beg:end) = nan
    pnf%m_deadcrootn_to_litter_fire(beg:end) = nan
    pnf%m_deadcrootn_storage_to_litter_fire(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_retransn_to_litter_fire(beg:end) = nan



    pnf%leafn_xfer_to_leafn(beg:end) = nan
    pnf%frootn_xfer_to_frootn(beg:end) = nan
    pnf%livestemn_xfer_to_livestemn(beg:end) = nan
    pnf%deadstemn_xfer_to_deadstemn(beg:end) = nan
    pnf%livecrootn_xfer_to_livecrootn(beg:end) = nan
    pnf%deadcrootn_xfer_to_deadcrootn(beg:end) = nan
    pnf%leafn_to_litter(beg:end) = nan
    pnf%leafn_to_retransn(beg:end) = nan
    pnf%frootn_to_retransn(beg:end) = nan
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
       pnf%fert(beg:end)                    = nan
       pnf%soyfixn(beg:end)                 = nan
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
    use shr_megan_mod, only: shr_megan_megcomps_n
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

    if (shr_megan_megcomps_n<1) return

    allocate(pvf%vocflx_tot(beg:end))
    allocate(pvf%vocflx(beg:end,1:shr_megan_megcomps_n))
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
    pvf%vocflx(beg:end,1:shr_megan_megcomps_n) = nan
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
   
   !F. Li and S. Levis
    allocate(cps%gdp_lf(beg:end))  
    allocate(cps%peatf_lf(beg:end))  
    allocate(cps%abm_lf(beg:end))  
    allocate(cps%lgdp_col(beg:end))  
    allocate(cps%lgdp1_col(beg:end))
    allocate(cps%lpop_col(beg:end))  

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
    allocate(cps%snow_depth(beg:end))
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
    allocate(cps%wf2(beg:end))
!   allocate(cps%xirrig(beg:end))
    allocate(cps%max_dayl(beg:end))
    allocate(cps%soilpsi(beg:end,nlevgrnd))
    allocate(cps%decl(beg:end))
    allocate(cps%coszen(beg:end))
    allocate(cps%bd(beg:end,nlevgrnd))
    allocate(cps%fpi(beg:end))
    allocate(cps%fpi_vr(beg:end,1:nlevdecomp_full))
    allocate(cps%fpg(beg:end))
    allocate(cps%annsum_counter(beg:end))
    allocate(cps%cannsum_npp(beg:end))
    allocate(cps%col_lag_npp(beg:end))
    allocate(cps%cannavg_t2m(beg:end))
    

    ! fire-related variables changed by F. Li and S. Levis
    allocate(cps%nfire(beg:end))
    allocate(cps%farea_burned(beg:end))
    allocate(cps%fsr_col(beg:end))
    allocate(cps%fd_col(beg:end))
    allocate(cps%cropf_col(beg:end))
    allocate(cps%prec10_col(beg:end))
    allocate(cps%prec60_col(beg:end))
    allocate(cps%lfc(beg:end))
    allocate(cps%lfc2(beg:end))
    allocate(cps%trotr1_col(beg:end))
    allocate(cps%trotr2_col(beg:end))
    allocate(cps%dtrotr_col(beg:end))
    allocate(cps%baf_crop(beg:end))
    allocate(cps%baf_peatf(beg:end))
    allocate(cps%fbac(beg:end))
    allocate(cps%fbac1(beg:end))
    allocate(cps%btran_col(beg:end))
    allocate(cps%wtlf(beg:end))
    allocate(cps%lfwt(beg:end))

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
! New variables for "S" Lakes
    allocate(cps%ws(beg:end))
    allocate(cps%ks(beg:end))
    allocate(cps%dz_lake(beg:end,nlevlak))
    allocate(cps%z_lake(beg:end,nlevlak))
    allocate(cps%savedtke1(beg:end))
    allocate(cps%cellsand(beg:end,nlevsoi))
    allocate(cps%cellclay(beg:end,nlevsoi))
    allocate(cps%cellorg(beg:end,nlevsoi))
    allocate(cps%lakedepth(beg:end))
    allocate(cps%etal(beg:end))
    allocate(cps%lakefetch(beg:end))
    allocate(cps%ust_lake(beg:end))
! End new variables for S Lakes
#if (defined VICHYDRO)
  ! new variables for VIC hydrology 
    allocate(cps%b_infil(beg:end))
    allocate(cps%dsmax(beg:end))
    allocate(cps%ds(beg:end))
    allocate(cps%Wsvic(beg:end))
    allocate(cps%c_param(beg:end))
    allocate(cps%expt(beg:end, nlayer))
    allocate(cps%ksat(beg:end, nlayer))
    allocate(cps%phi_s(beg:end, nlayer))
    allocate(cps%depth(beg:end, nlayert))
    allocate(cps%porosity(beg:end, nlayer))
    allocate(cps%max_moist(beg:end, nlayer))
    allocate(cps%vic_clm_fract(beg:end, nlayer, nlevsoi))
#endif
#ifdef LCH4
! New variable for finundated parameterization
    allocate(cps%zwt0(beg:end))
    allocate(cps%f0(beg:end))
    allocate(cps%p3(beg:end))
! New variable for methane
    allocate(cps%pH(beg:end))
#endif
    allocate(cps%irrig_rate(beg:end))
    allocate(cps%n_irrig_steps_left(beg:end))
    allocate(cps%forc_pbot(beg:end))
    allocate(cps%forc_rho(beg:end))
    allocate(cps%glc_topo(beg:end))

    allocate(cps%rf_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cps%pathfrac_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cps%nfixation_prof(beg:end,1:nlevdecomp_full))
    allocate(cps%ndep_prof(beg:end,1:nlevdecomp_full))
    allocate(cps%alt(beg:end))
    allocate(cps%altmax(beg:end))
    allocate(cps%altmax_lastyear(beg:end))
    allocate(cps%alt_indx(beg:end))
    allocate(cps%altmax_indx(beg:end))
    allocate(cps%altmax_lastyear_indx(beg:end))
    allocate(cps%som_adv_coef(beg:end,1:nlevdecomp_full))
    allocate(cps%som_diffus_coef(beg:end,1:nlevdecomp_full))

    cps%isoicol(beg:end) = bigint

    !F. Li and S. Levis
    cps%gdp_lf(beg:end) = nan
    cps%peatf_lf(beg:end) = nan
    cps%abm_lf(beg:end) = 13 
    cps%lgdp_col(beg:end) = nan
    cps%lgdp1_col(beg:end) = nan
    cps%lpop_col(beg:end) = nan

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
    cps%snow_depth(beg:end) = nan
    cps%snowdp(beg:end) = nan
    cps%frac_sno(beg:end) = nan
    cps%zi(beg:end,-nlevsno+0:nlevgrnd) = nan
    cps%dz(beg:end,-nlevsno+1:nlevgrnd) = nan
    cps%z (beg:end,-nlevsno+1:nlevgrnd) = nan
    cps%frac_iceold(beg:end,-nlevsno+1:nlevgrnd) = spval
    cps%imelt(beg:end,-nlevsno+1:nlevgrnd) = bigint
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
    cps%rootfr_road_perv(beg:end,1:nlevgrnd) = nan
    cps%rootr_road_perv(beg:end,1:nlevgrnd) = nan
    cps%wf(beg:end) = nan
!   cps%xirrig(beg:end) = 0._r8
    cps%soilpsi(beg:end,1:nlevgrnd) = spval
    cps%decl(beg:end) = nan
    cps%coszen(beg:end) = nan
    cps%bd(beg:end,1:nlevgrnd) = spval
    cps%fpi(beg:end) = nan
    cps%fpi_vr(beg:end,1:nlevdecomp_full) = nan
    cps%fpg(beg:end) = nan
    cps%annsum_counter(beg:end) = nan
    cps%cannsum_npp(beg:end) = nan
    cps%col_lag_npp(beg:end) = spval
    cps%cannavg_t2m(beg:end) = nan


    ! fire-related varibles changed by F. Li and S. Levis
    cps%nfire(beg:end) = spval
    cps%farea_burned(beg:end) = nan
    cps%btran_col(beg:end) = nan
    cps%wtlf(beg:end) = nan
    cps%lfwt(beg:end) = nan
    cps%fsr_col(beg:end) = nan
    cps%fd_col(beg:end) = nan
    cps%cropf_col(beg:end) = nan
    cps%baf_crop(beg:end) = nan
    cps%baf_peatf(beg:end) = nan
    cps%fbac(beg:end) = nan
    cps%fbac1(beg:end) = nan
    cps%trotr1_col(beg:end) = 0._r8
    cps%trotr2_col(beg:end) = 0._r8
    cps%dtrotr_col(beg:end) = 0._r8
    cps%prec10_col(beg:end) = nan
    cps%prec60_col(beg:end) = nan
    cps%lfc(beg:end) = spval
    cps%lfc2(beg:end) = 0._r8

    cps%albsnd_hst(beg:end,:numrad) = spval
    cps%albsni_hst(beg:end,:numrad) = spval
    cps%albsod(beg:end,:numrad)     = spval
    cps%albsoi(beg:end,:numrad)     = spval
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
! New variables for "S" Lakes
    cps%ws(beg:end) = nan
    cps%ks(beg:end) = nan
    cps%dz_lake(beg:end,1:nlevlak) = nan
    cps%z_lake(beg:end,1:nlevlak) = nan
    cps%savedtke1(beg:end) = spval  ! Initialize to spval so that c->g averaging will be done properly
    cps%cellsand(beg:end,1:nlevsoi) = nan
    cps%cellclay(beg:end,1:nlevsoi) = nan
    cps%cellorg(beg:end,1:nlevsoi) = nan
    cps%lakedepth(beg:end) = spval  ! Initialize to spval so that it can be a placeholder
                                    ! for future file input
    cps%etal(beg:end) = nan
    cps%lakefetch(beg:end) = nan
    cps%ust_lake(beg:end) = spval   ! Initial to spval to detect input from restart file if not arbinit
! End new variables for S Lakes
#ifdef LCH4
    cps%zwt0(beg:end) = nan
    cps%f0(beg:end)   = nan
    cps%p3(beg:end)   = nan
! New variable for methane
    cps%pH(beg:end)   = nan
#endif

    cps%irrig_rate(beg:end) = nan
    cps%n_irrig_steps_left(beg:end) = 0
    cps%forc_pbot(beg:end) = nan
    cps%forc_rho(beg:end) = nan
    cps%glc_topo(beg:end) = nan

    cps%rf_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cps%pathfrac_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cps%nfixation_prof(beg:end,1:nlevdecomp_full) = spval
    cps%ndep_prof(beg:end,1:nlevdecomp_full) = spval
    cps%alt(beg:end) = spval
    cps%altmax(beg:end) = spval
    cps%altmax_lastyear(beg:end) = spval
    cps%alt_indx(beg:end) = bigint
    cps%altmax_indx(beg:end) = bigint
    cps%altmax_lastyear_indx(beg:end) = bigint
    cps%som_adv_coef(beg:end,1:nlevdecomp_full) = spval
    cps%som_diffus_coef(beg:end,1:nlevdecomp_full) = spval

    allocate(cps%frac_sno_eff(beg:end))
    allocate(cps%topo_std(beg:end))
    allocate(cps%topo_ndx(beg:end))
    allocate(cps%topo_slope(beg:end))
    allocate(cps%hksat_min(beg:end,nlevgrnd))
    allocate(cps%frac_h2osfc(beg:end))
    allocate(cps%micro_sigma(beg:end))
    allocate(cps%h2osfc_thresh(beg:end))
    allocate(cps%frac_h2osfc_temp(beg:end))
    allocate(cps%n_melt(beg:end))

    cps%frac_sno_eff(beg:end) = spval
    cps%topo_std(beg:end) = nan
    cps%topo_ndx(beg:end) = nan
    cps%topo_slope(beg:end) = nan
    cps%hksat_min(beg:end,1:nlevgrnd) = nan
    cps%frac_h2osfc(beg:end) = spval
    cps%micro_sigma(beg:end) = nan
    cps%h2osfc_thresh(beg:end) = nan
    cps%frac_h2osfc_temp(beg:end) = 0.0_r8
    cps%n_melt(beg:end) = nan
#if (defined VICHYDRO)
   ! new variables for VIC hydrology
    cps%b_infil(beg:end)  = nan
    cps%dsmax(beg:end)    = nan
    cps%ds(beg:end)       = nan
    cps%Wsvic(beg:end)    = nan
    cps%c_param(beg:end)  = nan
    cps%expt(beg:end, 1:nlayer)     = nan
    cps%ksat(beg:end, 1:nlayer)     = nan
    cps%phi_s(beg:end, 1:nlayer)    = nan
    cps%depth(beg:end, 1:nlayert)    = nan
    cps%porosity(beg:end, 1:nlayer)   = nan
    cps%max_moist(beg:end, 1:nlayer)   = nan
    cps%vic_clm_fract(beg:end,1:nlayer,1:nlevsoi) = nan
#endif

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
    allocate(ces%tsoi17(beg:end))
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
    ces%tsoi17(beg:end) = spval
    ces%t_lake(beg:end,1:nlevlak)            = nan
    ces%tssbef(beg:end,-nlevsno+1:nlevgrnd)   = nan
    ces%thv(beg:end)       = nan
    ces%hc_soi(beg:end)    = nan
    ces%hc_soisno(beg:end) = nan
    ces%forc_t(beg:end) = nan
    ces%forc_th(beg:end) = nan

    allocate(ces%t_h2osfc(beg:end))
    allocate(ces%t_h2osfc_bef(beg:end))

    ces%t_h2osfc(beg:end) = spval
    ces%t_h2osfc_bef(beg:end)   = nan
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
!New variable for methane code
#ifdef LCH4
    allocate(cws%finundated(beg:end))
#endif
    allocate(cws%wa(beg:end))
    allocate(cws%qcharge(beg:end))
    allocate(cws%smp_l(beg:end,1:nlevgrnd))
    allocate(cws%hk_l(beg:end,1:nlevgrnd))
! New variables for "S" lakes
    allocate(cws%lake_icefrac(beg:end,1:nlevlak))
    allocate(cws%lake_icethick(beg:end))
! End new variables for S lakes
    allocate(cws%forc_q(beg:end))
#if (defined VICHYDRO)
    allocate(cws%moist(beg:end,1:nlayert))
    allocate(cws%ice(beg:end,1:nlayert))
    allocate(cws%moist_vol(beg:end,1:nlayert))
    allocate(cws%max_infil(beg:end))
    allocate(cws%i_0(beg:end))
#endif

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
#ifdef LCH4
    cws%finundated(beg:end) = nan
#endif
    cws%wa(beg:end) = spval
    cws%qcharge(beg:end) = nan
    cws%smp_l(beg:end,1:nlevgrnd) = spval
    cws%hk_l(beg:end,1:nlevgrnd) = spval
! New variables for "S" lakes
    cws%lake_icefrac(beg:end,1:nlevlak) = spval  ! Initialize to spval for detection of whether
                                                 ! it has not been initialized by file input
                                                 ! and so c->g averaging will be done properly
    cws%lake_icethick(beg:end) = nan
! End new variables for S lakes
    cws%forc_q(beg:end) = nan

    allocate(cws%h2osfc(beg:end))
    allocate(cws%qg_snow(beg:end))
    allocate(cws%qg_soil(beg:end))
    allocate(cws%qg_h2osfc(beg:end))
    allocate(cws%frost_table(beg:end))
    allocate(cws%zwt_perched(beg:end))
    allocate(cws%int_snow(beg:end))
    allocate(cws%swe_old(beg:end,-nlevsno+1:0))

    cws%h2osfc(beg:end)      = spval
    cws%qg_snow(beg:end)     = nan
    cws%qg_soil(beg:end)     = nan
    cws%qg_h2osfc(beg:end)   = nan
    cws%frost_table(beg:end) = spval
    cws%zwt_perched(beg:end) = spval
    cws%int_snow(beg:end)    = spval
    cws%swe_old(beg:end,-nlevsno+1:0)= nan
#if (defined VICHYDRO)
    cws%moist(beg:end,1:nlayert) = spval
    cws%ice(beg:end,1:nlayert) = spval
    cws%moist_vol(beg:end,1:nlayert) = spval
    cws%max_infil(beg:end) = spval
    cws%i_0(beg:end) = spval
#endif

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
! !USES:
    use clm_varcon, only : spval
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
    allocate(ccs%col_ctrunc(beg:end))
    allocate(ccs%decomp_cpools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccs%decomp_cpools(beg:end,1:ndecomp_pools))
    allocate(ccs%decomp_cpools_1m(beg:end,1:ndecomp_pools))
    allocate(ccs%col_ctrunc_vr(beg:end,1:nlevdecomp_full))
    allocate(ccs%seedc(beg:end))
    allocate(ccs%prod10c(beg:end))
    allocate(ccs%prod100c(beg:end))
    allocate(ccs%totprodc(beg:end))
    allocate(ccs%totlitc(beg:end))
    allocate(ccs%totsomc(beg:end))
    allocate(ccs%totlitc_1m(beg:end))
    allocate(ccs%totsomc_1m(beg:end))
    allocate(ccs%totecosysc(beg:end))
    allocate(ccs%totcolc(beg:end))

    !F. Li and S. Levis
    allocate(ccs%rootc_col(beg:end))
    allocate(ccs%totvegc_col(beg:end))
    allocate(ccs%leafc_col(beg:end))
    allocate(ccs%fuelc(beg:end))
    allocate(ccs%fuelc_crop(beg:end))

    ccs%soilc(beg:end) = nan
    ccs%cwdc(beg:end) = nan
    ccs%decomp_cpools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan
    ccs%decomp_cpools(beg:end,1:ndecomp_pools) = nan
    ccs%decomp_cpools_1m(beg:end,1:ndecomp_pools) = nan
    ccs%col_ctrunc(beg:end) = nan
    ccs%col_ctrunc_vr(beg:end,1:nlevdecomp_full) = nan
    ccs%seedc(beg:end) = nan
    ccs%prod10c(beg:end) = nan
    ccs%prod100c(beg:end) = nan
    ccs%totprodc(beg:end) = nan
    ccs%totlitc(beg:end) = nan
    ccs%totsomc(beg:end) = nan
    ccs%totlitc_1m(beg:end) = nan
    ccs%totsomc_1m(beg:end) = nan
    ccs%totecosysc(beg:end) = nan
    ccs%totcolc(beg:end) = nan

    ccs%rootc_col(beg:end) = nan
    ccs%totvegc_col(beg:end) = nan
    ccs%leafc_col(beg:end) = nan
    ccs%fuelc(beg:end) = spval
    ccs%fuelc_crop(beg:end) = nan

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
! !USES:
    use clm_varcon, only : spval
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


    allocate(cns%decomp_npools(beg:end,1:ndecomp_pools))
    allocate(cns%decomp_npools_1m(beg:end,1:ndecomp_pools))
    allocate(cns%decomp_npools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(cns%sminn_vr(beg:end,1:nlevdecomp_full))
    allocate(cns%col_ntrunc_vr(beg:end,1:nlevdecomp_full))
#ifdef NITRIF_DENITRIF
    allocate(cns%smin_no3_vr(beg:end,1:nlevdecomp_full))
    allocate(cns%smin_nh4_vr(beg:end,1:nlevdecomp_full))
    allocate(cns%smin_no3(beg:end))
    allocate(cns%smin_nh4(beg:end))
#endif
    allocate(cns%cwdn(beg:end))
    allocate(cns%sminn(beg:end))
    allocate(cns%col_ntrunc(beg:end))
    allocate(cns%seedn(beg:end))
    allocate(cns%prod10n(beg:end))
    allocate(cns%prod100n(beg:end))
    allocate(cns%totprodn(beg:end))
    allocate(cns%totlitn(beg:end))
    allocate(cns%totsomn(beg:end))
    allocate(cns%totlitn_1m(beg:end))
    allocate(cns%totsomn_1m(beg:end))
    allocate(cns%totecosysn(beg:end))
    allocate(cns%totcoln(beg:end))

    cns%decomp_npools(beg:end,1:ndecomp_pools) = nan
    cns%decomp_npools_1m(beg:end,1:ndecomp_pools) = nan
    cns%decomp_npools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan
    cns%sminn_vr(beg:end,1:nlevdecomp_full) = nan
    cns%col_ntrunc_vr(beg:end,1:nlevdecomp_full) = nan
#ifdef NITRIF_DENITRIF
    cns%smin_no3_vr(beg:end,1:nlevdecomp_full) = nan
    cns%smin_nh4_vr(beg:end,1:nlevdecomp_full) = nan
    cns%smin_no3(beg:end) = nan
    cns%smin_nh4(beg:end) = nan
#endif
    cns%cwdn(beg:end) = nan
    cns%sminn(beg:end) = nan
    cns%col_ntrunc(beg:end) = nan
    cns%seedn(beg:end) = nan
    cns%prod10n(beg:end) = nan
    cns%prod100n(beg:end) = nan
    cns%totprodn(beg:end) = nan
    cns%totlitn(beg:end) = nan
    cns%totsomn(beg:end) = nan
    cns%totlitn_1m(beg:end) = nan
    cns%totsomn_1m(beg:end) = nan
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
! !USES:
    use clm_varcon, only : spval
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
    allocate(cef%eflx_fgr(beg:end, 1:nlevgrnd))
    allocate(cef%eflx_building_heat(beg:end))
    allocate(cef%eflx_urban_ac(beg:end))
    allocate(cef%eflx_urban_heat(beg:end))
    allocate(cef%eflx_bot(beg:end))

    cef%eflx_snomelt(beg:end)       = spval
    cef%eflx_snomelt_u(beg:end)       = spval
    cef%eflx_snomelt_r(beg:end)       = spval
    cef%eflx_impsoil(beg:end)       = nan
    cef%eflx_fgr12(beg:end)         = nan
    cef%eflx_fgr(beg:end, 1:nlevgrnd) = nan
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
    cwf%qflx_rsub_sat(beg:end) = spval
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
    cwf%qflx_irrig(beg:end)  = spval
    cwf%qflx_glcice(beg:end) = nan
    cwf%qflx_glcice_frz(beg:end) = nan
    cwf%qflx_glcice_melt(beg:end) = spval
    cwf%glc_rofi(beg:end)    = nan
    cwf%glc_rofl(beg:end)    = nan
    cwf%qflx_floodc(beg:end) = spval

    allocate(cwf%qflx_h2osfc_to_ice(beg:end))
    allocate(cwf%qflx_h2osfc_surf(beg:end))
    allocate(cwf%qflx_snow_h2osfc(beg:end))
    allocate(cwf%qflx_drain_perched(beg:end))
    allocate(cwf%qflx_floodc(beg:end))
    allocate(cwf%qflx_snow_melt(beg:end))

    cwf%qflx_h2osfc_to_ice(beg:end) = spval
    cwf%qflx_h2osfc_surf(beg:end) = spval
    cwf%qflx_snow_h2osfc(beg:end) = nan
    cwf%qflx_drain_perched(beg:end) = spval
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
  use clm_varcon, only : spval
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

    allocate(ccf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(ccf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(ccf%m_decomp_cpools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccf%m_decomp_cpools_to_fire(beg:end,1:ndecomp_pools))
    allocate(ccf%decomp_cascade_hr_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_hr(beg:end,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_ctransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_ctransfer(beg:end,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cpools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccf%decomp_k(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%t_scalar(beg:end,1:nlevdecomp_full))
    allocate(ccf%w_scalar(beg:end,1:nlevdecomp_full))
    allocate(ccf%hr_vr(beg:end,1:nlevdecomp_full))
    allocate(ccf%o_scalar(beg:end,1:nlevdecomp_full))
    allocate(ccf%som_c_leached(beg:end))
    allocate(ccf%decomp_cpools_leached(beg:end,1:ndecomp_pools))
    allocate(ccf%decomp_cpools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools))

    allocate(ccf%phenology_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%phenology_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%phenology_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full))
    allocate(ccf%fire_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_met_fire(beg:end,1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_cel_fire(beg:end,1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_lig_fire(beg:end,1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_cwdc(beg:end, 1:nlevdecomp_full))

#ifdef NITRIF_DENITRIF
    allocate(ccf%phr_vr(beg:end,1:nlevdecomp_full))
#endif

#ifdef CN
   !F. Li and S. Levis
    allocate(ccf%somc_fire(beg:end))
    allocate(ccf%lf_conv_cflux(beg:end))
    allocate(ccf%dwt_seedc_to_leaf(beg:end))
    allocate(ccf%dwt_seedc_to_deadstem(beg:end))
    allocate(ccf%dwt_conv_cflux(beg:end))
    allocate(ccf%dwt_prod10c_gain(beg:end))
    allocate(ccf%dwt_prod100c_gain(beg:end))
    allocate(ccf%dwt_frootc_to_litr_met_c(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_frootc_to_litr_cel_c(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_frootc_to_litr_lig_c(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_livecrootc_to_cwdc(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_deadcrootc_to_cwdc(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_closs(beg:end))
    allocate(ccf%landuseflux(beg:end))
    allocate(ccf%landuptake(beg:end))
    allocate(ccf%prod10c_loss(beg:end))
    allocate(ccf%prod100c_loss(beg:end))
    allocate(ccf%product_closs(beg:end))
#endif
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

#if (defined CN)
    allocate(ccf%cwdc_hr(beg:end))
    allocate(ccf%cwdc_loss(beg:end))
    allocate(ccf%litterc_loss(beg:end))
#endif

    ccf%m_c_to_litr_met_fire(beg:end,1:nlevdecomp_full)             = nan
    ccf%m_c_to_litr_cel_fire(beg:end,1:nlevdecomp_full)             = nan
    ccf%m_c_to_litr_lig_fire(beg:end,1:nlevdecomp_full)             = nan
    ccf%hrv_deadstemc_to_prod10c(beg:end)         = nan        
    ccf%hrv_deadstemc_to_prod100c(beg:end)        = nan       
    ccf%m_decomp_cpools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                = nan
    ccf%m_decomp_cpools_to_fire(beg:end,1:ndecomp_pools)                                     = nan
    ccf%decomp_cascade_hr_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions)        = nan
    ccf%decomp_cascade_hr(beg:end,1:ndecomp_cascade_transitions)                             = nan
    ccf%decomp_cascade_ctransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    ccf%decomp_cascade_ctransfer(beg:end,1:ndecomp_cascade_transitions)                      = nan
    ccf%decomp_cpools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                  = nan
    ccf%decomp_k(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions)                    = spval
! Initialize these four below to spval to allow history to not average over inactive points.
    ccf%t_scalar(beg:end,1:nlevdecomp_full)                         = spval
    ccf%w_scalar(beg:end,1:nlevdecomp_full)                         = spval
    ccf%hr_vr(beg:end, 1:nlevdecomp_full)                           = spval
    ccf%o_scalar(beg:end, 1:nlevdecomp_full)                        = spval
    ccf%som_c_leached(beg:end)                                                      = nan 
    ccf%decomp_cpools_leached(beg:end,1:ndecomp_pools)                              = nan
    ccf%decomp_cpools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan

    ccf%phenology_c_to_litr_met_c(beg:end, 1:nlevdecomp_full)                     = nan
    ccf%phenology_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full)                     = nan
    ccf%phenology_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full)                     = nan
    ccf%gap_mortality_c_to_litr_met_c(beg:end, 1:nlevdecomp_full)                 = nan
    ccf%gap_mortality_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full)                 = nan
    ccf%gap_mortality_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full)                 = nan
    ccf%gap_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%fire_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full)                      = nan
    ccf%harvest_c_to_litr_met_c(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%harvest_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%harvest_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%harvest_c_to_cwdc(beg:end, 1:nlevdecomp_full)                             = nan

#ifdef NITRIF_DENITRIF
    ccf%phr_vr(beg:end,1:nlevdecomp_full)                              = nan
#endif
#if (defined CN)
    !F. Li and S. Levis
    ccf%somc_fire(beg:end)                        = nan
    ccf%lf_conv_cflux(beg:end)                    = nan
    ccf%dwt_seedc_to_leaf(beg:end)                = nan
    ccf%dwt_seedc_to_deadstem(beg:end)            = nan
    ccf%dwt_conv_cflux(beg:end)                   = nan
    ccf%dwt_prod10c_gain(beg:end)                 = nan
    ccf%dwt_prod100c_gain(beg:end)                = nan
    ccf%dwt_frootc_to_litr_met_c(beg:end,1:nlevdecomp_full)             = nan
    ccf%dwt_frootc_to_litr_cel_c(beg:end,1:nlevdecomp_full)             = nan
    ccf%dwt_frootc_to_litr_lig_c(beg:end,1:nlevdecomp_full)             = nan
    ccf%dwt_livecrootc_to_cwdc(beg:end,1:nlevdecomp_full)           = nan
    ccf%dwt_deadcrootc_to_cwdc(beg:end,1:nlevdecomp_full)           = nan
    ccf%dwt_closs(beg:end)                        = nan
    ccf%landuseflux(beg:end)                      = nan
    ccf%landuptake(beg:end)                       = nan
    ccf%prod10c_loss(beg:end)                     = nan
    ccf%prod100c_loss(beg:end)                    = nan
    ccf%product_closs(beg:end)                    = nan
#endif
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

#if (defined CN)
    ccf%cwdc_hr(beg:end)                          = nan
    ccf%cwdc_loss(beg:end)                        = nan
    ccf%litterc_loss(beg:end)                     = nan
#endif

  end subroutine init_column_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_ch4_type
!
! !INTERFACE:
#ifdef LCH4
  subroutine init_column_ch4_type(beg, end, cch4)
!
! !DESCRIPTION:
! Initialize column methane flux variables
!
  use clm_varcon, only : spval
  use clm_varpar, only : ngases
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_ch4_type), intent(inout):: cch4
!
! !REVISION HISTORY:
! Created by William J. Riley
!
!EOP
!------------------------------------------------------------------------

    allocate(cch4%ch4_prod_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_prod_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_prod_depth_lake(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_lake(beg:end,1:nlevgrnd))
    allocate(cch4%o2_oxid_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_oxid_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_decomp_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_decomp_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_aere_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_aere_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_decomp_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_decomp_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_oxid_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_oxid_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_aere_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_aere_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_tran_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_tran_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_aere_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_aere_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_surf_aere_sat(beg:end))
    allocate(cch4%ch4_surf_aere_unsat(beg:end))
    allocate(cch4%ch4_ebul_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_ebul_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_ebul_total_sat(beg:end))
    allocate(cch4%ch4_ebul_total_unsat(beg:end))
    allocate(cch4%ch4_surf_ebul_sat(beg:end))
    allocate(cch4%ch4_surf_ebul_unsat(beg:end))
    allocate(cch4%ch4_surf_ebul_lake(beg:end))
    allocate(cch4%conc_ch4_sat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_ch4_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_ch4_lake(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_surf_diff_sat(beg:end))
    allocate(cch4%ch4_surf_diff_unsat(beg:end))
    allocate(cch4%ch4_surf_diff_lake(beg:end))
    allocate(cch4%conc_o2_sat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_o2_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_o2_lake(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_dfsat_flux(beg:end))
    allocate(cch4%zwt_ch4_unsat(beg:end))
    allocate(cch4%fsat_bef(beg:end))
    allocate(cch4%lake_soilc(beg:end,1:nlevgrnd))
    allocate(cch4%lake_raw(beg:end))
    allocate(cch4%totcolch4(beg:end))
    allocate(cch4%fphr(beg:end,1:nlevgrnd))
    allocate(cch4%annsum_counter(beg:end))
    allocate(cch4%tempavg_somhr(beg:end))
    allocate(cch4%annavg_somhr(beg:end))
    allocate(cch4%tempavg_finrw(beg:end))
    allocate(cch4%annavg_finrw(beg:end))
    allocate(cch4%sif(beg:end))
    allocate(cch4%o2stress_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%o2stress_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4stress_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4stress_sat(beg:end,1:nlevgrnd))
    allocate(cch4%qflx_surf_lag(beg:end))
    allocate(cch4%finundated_lag(beg:end))
    allocate(cch4%layer_sat_lag(beg:end,1:nlevgrnd))


    cch4%ch4_prod_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_prod_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_prod_depth_lake(beg:end,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_lake(beg:end,1:nlevgrnd) = nan
    cch4%o2_oxid_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%o2_oxid_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%o2_decomp_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%o2_decomp_depth_unsat(beg:end,1:nlevgrnd) = spval ! To detect first time-step for denitrification code
    cch4%o2_aere_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%o2_aere_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%co2_decomp_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%co2_decomp_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%co2_oxid_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%co2_oxid_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_aere_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_aere_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_tran_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_tran_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%co2_aere_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%co2_aere_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_surf_aere_sat(beg:end) = nan
    cch4%ch4_surf_aere_unsat(beg:end) = nan
    cch4%ch4_ebul_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_ebul_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_ebul_total_sat(beg:end) = nan
    cch4%ch4_ebul_total_unsat(beg:end) = nan
    cch4%ch4_surf_ebul_sat(beg:end) = nan
    cch4%ch4_surf_ebul_unsat(beg:end) = nan
    cch4%ch4_surf_ebul_lake(beg:end) = nan
    cch4%conc_ch4_sat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%conc_ch4_unsat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%conc_ch4_lake(beg:end,1:nlevgrnd) = nan ! Just a diagnostic, so nan is fine
    cch4%ch4_surf_diff_sat(beg:end) = nan
    cch4%ch4_surf_diff_unsat(beg:end) = nan
    cch4%ch4_surf_diff_lake(beg:end) = nan
    cch4%conc_o2_sat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%conc_o2_unsat(beg:end,1:nlevgrnd) = spval ! To detect file input and detect first time-step for denitrification code
    cch4%conc_o2_lake(beg:end,1:nlevgrnd) = nan ! Just a diagnostic, so nan is fine
    cch4%ch4_dfsat_flux(beg:end) = nan
    cch4%zwt_ch4_unsat(beg:end) = nan
    cch4%fsat_bef(beg:end) = spval ! To detect first time-step
    cch4%lake_soilc(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%lake_raw(beg:end) = nan
    cch4%totcolch4(beg:end) = spval ! To detect first time-step
    cch4%fphr(beg:end,1:nlevgrnd) = nan
    cch4%annsum_counter(beg:end) = spval ! To detect first time-step
    cch4%tempavg_somhr(beg:end) = nan
    cch4%annavg_somhr(beg:end) = spval ! To detect first year
    cch4%tempavg_finrw(beg:end) = nan
    cch4%annavg_finrw(beg:end) = spval ! To detect first year
    cch4%sif(beg:end) = nan
    cch4%o2stress_unsat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%o2stress_sat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%ch4stress_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4stress_sat(beg:end,1:nlevgrnd) = nan
    cch4%qflx_surf_lag(beg:end) = spval ! To detect file input
    cch4%finundated_lag(beg:end) = spval ! To detect file input
    cch4%layer_sat_lag(beg:end,1:nlevgrnd) = spval ! To detect file input


! def CH4

  end subroutine init_column_ch4_type
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nflux_type
!
! !INTERFACE:
  subroutine init_column_nflux_type(beg, end, cnf)
!
  use clm_varcon, only : spval
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
    allocate(cnf%fert_to_sminn(beg:end))
    allocate(cnf%soyfixn_to_sminn(beg:end))    
    allocate(cnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(cnf%hrv_deadstemn_to_prod100n(beg:end))       

    allocate(cnf%m_n_to_litr_met_fire(beg:end,1:nlevdecomp_full))
    allocate(cnf%m_n_to_litr_cel_fire(beg:end,1:nlevdecomp_full))
    allocate(cnf%m_n_to_litr_lig_fire(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_to_plant(beg:end))
    allocate(cnf%potential_immob(beg:end))
    allocate(cnf%actual_immob(beg:end))
    allocate(cnf%gross_nmin(beg:end))
    allocate(cnf%net_nmin(beg:end))
    allocate(cnf%denit(beg:end))
    allocate(cnf%supplement_to_sminn(beg:end))
    allocate(cnf%m_decomp_npools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(cnf%m_decomp_npools_to_fire(beg:end,1:ndecomp_pools))
    allocate(cnf%decomp_cascade_ntransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_ntransfer(beg:end,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_sminn_flux_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_sminn_flux(beg:end,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_npools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools))

    allocate(cnf%phenology_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%phenology_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%phenology_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full))
    allocate(cnf%fire_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_cwdn(beg:end, 1:nlevdecomp_full))

#ifndef NITRIF_DENITRIF
    allocate(cnf%sminn_to_denit_decomp_cascade_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%sminn_to_denit_decomp_cascade(beg:end,1:ndecomp_cascade_transitions))
    allocate(cnf%sminn_to_denit_excess_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_to_denit_excess(beg:end))
    allocate(cnf%sminn_leached_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_leached(beg:end))
#else
    allocate(cnf%f_nit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_leached_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_leached(beg:end))
    allocate(cnf%smin_no3_runoff_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_runoff(beg:end))
    allocate(cnf%pot_f_nit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%pot_f_nit(beg:end))
    allocate(cnf%pot_f_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%pot_f_denit(beg:end))
    allocate(cnf%actual_immob_no3_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%actual_immob_nh4_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_to_plant_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_nh4_to_plant_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_nit(beg:end))
    allocate(cnf%f_denit(beg:end))
    allocate(cnf%n2_n2o_ratio_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_n2o_denit(beg:end))
    allocate(cnf%f_n2o_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_n2o_nit(beg:end))
    allocate(cnf%f_n2o_nit_vr(beg:end,1:nlevdecomp_full))

    allocate(cnf%smin_no3_massdens_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%soil_bulkdensity(beg:end,1:nlevdecomp_full))
    allocate(cnf%k_nitr_t_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%k_nitr_ph_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%k_nitr_h2o_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%k_nitr_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%wfps_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%fmax_denit_carbonsubstrate_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%fmax_denit_nitrate_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%f_denit_base_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%diffus(beg:end,1:nlevdecomp_full))
    allocate(cnf%ratio_k1(beg:end,1:nlevdecomp_full))
    allocate(cnf%ratio_no3_co2(beg:end,1:nlevdecomp_full))
    allocate(cnf%soil_co2_prod(beg:end,1:nlevdecomp_full))
    allocate(cnf%fr_WFPS(beg:end,1:nlevdecomp_full))

    allocate(cnf%r_psi(beg:end,1:nlevdecomp_full))
    allocate(cnf%anaerobic_frac(beg:end,1:nlevdecomp_full))
#endif
    allocate(cnf%potential_immob_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%actual_immob_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_to_plant_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%supplement_to_sminn_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%gross_nmin_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%net_nmin_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_seedn_to_leaf(beg:end))
    allocate(cnf%dwt_seedn_to_deadstem(beg:end))
    allocate(cnf%dwt_conv_nflux(beg:end))
    allocate(cnf%dwt_prod10n_gain(beg:end))
    allocate(cnf%dwt_prod100n_gain(beg:end))
    allocate(cnf%dwt_frootn_to_litr_met_n(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_frootn_to_litr_cel_n(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_frootn_to_litr_lig_n(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_livecrootn_to_cwdn(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_deadcrootn_to_cwdn(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_nloss(beg:end))
    allocate(cnf%prod10n_loss(beg:end))
    allocate(cnf%prod100n_loss(beg:end))
    allocate(cnf%product_nloss(beg:end))
    allocate(cnf%col_ninputs(beg:end))
    allocate(cnf%col_noutputs(beg:end))
    allocate(cnf%col_fire_nloss(beg:end))
    allocate(cnf%som_n_leached(beg:end))
    allocate(cnf%decomp_npools_leached(beg:end,1:ndecomp_pools))
    allocate(cnf%decomp_npools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    
    cnf%ndep_to_sminn(beg:end) = nan
    cnf%nfix_to_sminn(beg:end) = nan
    cnf%fert_to_sminn(beg:end) = nan
    cnf%soyfixn_to_sminn(beg:end) = nan
    cnf%hrv_deadstemn_to_prod10n(beg:end) = nan        
    cnf%hrv_deadstemn_to_prod100n(beg:end) = nan       
    cnf%m_n_to_litr_met_fire(beg:end,1:nlevdecomp_full) = nan
    cnf%m_n_to_litr_cel_fire(beg:end,1:nlevdecomp_full) = nan
    cnf%m_n_to_litr_lig_fire(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_to_plant(beg:end) = nan
    cnf%potential_immob(beg:end) = nan
    cnf%actual_immob(beg:end) = nan
    cnf%gross_nmin(beg:end) = nan
    cnf%net_nmin(beg:end) = nan
    cnf%denit(beg:end) = nan
    cnf%supplement_to_sminn(beg:end) = nan
    cnf%m_decomp_npools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                 = nan
    cnf%m_decomp_npools_to_fire(beg:end,1:ndecomp_pools)                                      = nan
    cnf%decomp_cascade_ntransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions)  = nan
    cnf%decomp_cascade_ntransfer(beg:end,1:ndecomp_cascade_transitions)                       = nan
    cnf%decomp_cascade_sminn_flux_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cnf%decomp_cascade_sminn_flux(beg:end,1:ndecomp_cascade_transitions)                      = nan
    cnf%decomp_npools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                   = nan
    
    cnf%phenology_n_to_litr_met_n(beg:end, 1:nlevdecomp_full)                     = nan
    cnf%phenology_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full)                     = nan
    cnf%phenology_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full)                     = nan
    cnf%gap_mortality_n_to_litr_met_n(beg:end, 1:nlevdecomp_full)                 = nan
    cnf%gap_mortality_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full)                 = nan
    cnf%gap_mortality_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full)                 = nan
    cnf%gap_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%fire_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full)                      = nan
    cnf%harvest_n_to_litr_met_n(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%harvest_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%harvest_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%harvest_n_to_cwdn(beg:end, 1:nlevdecomp_full)                             = nan

#ifndef NITRIF_DENITRIF
    cnf%sminn_to_denit_decomp_cascade_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cnf%sminn_to_denit_decomp_cascade(beg:end,1:ndecomp_cascade_transitions) = nan
    cnf%sminn_to_denit_excess_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_to_denit_excess(beg:end) = nan
    cnf%sminn_leached_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_leached(beg:end) = nan
#else
    cnf%f_nit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_leached_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_leached(beg:end) = nan
    cnf%smin_no3_runoff_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_runoff(beg:end) = nan
    cnf%pot_f_nit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%pot_f_nit(beg:end) = nan
    cnf%pot_f_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%pot_f_denit(beg:end) = nan
    cnf%actual_immob_no3_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%actual_immob_nh4_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_to_plant_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_nh4_to_plant_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_nit(beg:end) = nan
    cnf%f_denit(beg:end) = nan
    cnf%n2_n2o_ratio_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_n2o_denit(beg:end) = nan
    cnf%f_n2o_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_n2o_nit(beg:end) = nan
    cnf%f_n2o_nit_vr(beg:end,1:nlevdecomp_full) = nan

    cnf%smin_no3_massdens_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%soil_bulkdensity(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_t_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_ph_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_h2o_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%wfps_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%fmax_denit_carbonsubstrate_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%fmax_denit_nitrate_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%f_denit_base_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%diffus(beg:end,1:nlevdecomp_full) = spval
    cnf%ratio_k1(beg:end,1:nlevdecomp_full) = nan
    cnf%ratio_no3_co2(beg:end,1:nlevdecomp_full) = spval
    cnf%soil_co2_prod(beg:end,1:nlevdecomp_full) = nan
    cnf%fr_WFPS(beg:end,1:nlevdecomp_full) = spval

    cnf%r_psi(beg:end,1:nlevdecomp_full) = spval
    cnf%anaerobic_frac(beg:end,1:nlevdecomp_full) = spval
#endif
    cnf%potential_immob_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%actual_immob_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_to_plant_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%supplement_to_sminn_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%gross_nmin_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%net_nmin_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_seedn_to_leaf(beg:end) = nan
    cnf%dwt_seedn_to_deadstem(beg:end) = nan
    cnf%dwt_conv_nflux(beg:end) = nan
    cnf%dwt_prod10n_gain(beg:end) = nan
    cnf%dwt_prod100n_gain(beg:end) = nan
    cnf%dwt_frootn_to_litr_met_n(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_frootn_to_litr_cel_n(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_frootn_to_litr_lig_n(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_livecrootn_to_cwdn(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_deadcrootn_to_cwdn(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_nloss(beg:end) = nan
    cnf%prod10n_loss(beg:end) = nan
    cnf%prod100n_loss(beg:end) = nan
    cnf%product_nloss(beg:end) = nan
    cnf%col_ninputs(beg:end) = nan
    cnf%col_noutputs(beg:end) = nan
    cnf%col_fire_nloss(beg:end) = nan
    cnf%som_n_leached(beg:end)                                                      = nan 
    cnf%decomp_npools_leached(beg:end,1:ndecomp_pools)                              = nan
    cnf%decomp_npools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan


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
    if ( nlevurb > 0 )then
      allocate(lps%tk_wall(beg:end,nlevurb))
      allocate(lps%tk_roof(beg:end,nlevurb))
      allocate(lps%cv_wall(beg:end,nlevurb))
      allocate(lps%cv_roof(beg:end,nlevurb))
    end if
    allocate(lps%tk_improad(beg:end,nlevurb))
    allocate(lps%cv_improad(beg:end,nlevurb))
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
    if ( nlevurb > 0 )then
      lps%tk_wall(beg:end,1:nlevurb) = nan
      lps%tk_roof(beg:end,1:nlevurb) = nan
      lps%cv_wall(beg:end,1:nlevurb) = nan
      lps%cv_roof(beg:end,1:nlevurb) = nan
    end if
    lps%tk_improad(beg:end,1:nlevurb) = nan
    lps%cv_improad(beg:end,1:nlevurb) = nan
    lps%thick_wall(beg:end) = nan
    lps%thick_roof(beg:end) = nan
    lps%nlev_improad(beg:end) = bigint
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

#if (defined CNDV)
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
#endif


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

    gwf%qflx_runoffg(beg:end) = 0._r8
    gwf%qflx_snwcp_iceg(beg:end) = 0._r8
    gwf%qflx_liq_dynbal(beg:end) = nan
    gwf%qflx_ice_dynbal(beg:end) = nan

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

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_ch4_type
!
! !INTERFACE:
#ifdef LCH4
  subroutine init_gridcell_ch4_type(beg, end, gch4)
!
! !DESCRIPTION:
! Initialize gridcell ch4 variables
!
  use clm_varpar, only: ngases
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_ch4_type), intent(inout):: gch4
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(gch4%c_atm(beg:end,1:ngases))
    allocate(gch4%ch4co2f(beg:end))
    allocate(gch4%ch4prodg(beg:end))
    allocate(gch4%nem(beg:end))

    gch4%c_atm(beg:end,1:ngases) = nan
    gch4%ch4co2f(beg:end) = nan
    gch4%ch4prodg(beg:end) = nan
    gch4%nem(beg:end) = nan

  end subroutine init_gridcell_ch4_type
#endif

end module clmtypeInitMod
