module seq_flux_mct
  
  use shr_kind_mod,      only: r8 => shr_kind_r8, in=>shr_kind_in
  use shr_sys_mod,       only: shr_sys_abort
  use shr_flux_mod,      only: shr_flux_atmocn, shr_flux_atmocn_diurnal
  use shr_orb_mod,       only: shr_orb_params, shr_orb_cosz, shr_orb_decl
  use shr_mct_mod,       only: shr_mct_queryConfigFile, shr_mct_sMatReaddnc

  use mct_mod
  use seq_flds_mod
  use seq_comm_mct
  use seq_infodata_mod

  use component_type_mod

  implicit none
  private 	
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public seq_flux_init_mct
  public seq_flux_initexch_mct

  public seq_flux_ocnalb_mct

  public seq_flux_atmocn_mct
  public seq_flux_atmocnexch_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  real(r8), pointer       :: lats(:)  ! latitudes  (degrees)
  real(r8), pointer       :: lons(:)  ! longitudes (degrees)
  integer(in),allocatable :: mask(:)  ! ocn domain mask: 0 <=> inactive cell
  integer(in),allocatable :: emask(:) ! ocn mask on exchange grid decomp

  real(r8), allocatable ::  uocn (:)  ! ocn velocity, zonal
  real(r8), allocatable ::  vocn (:)  ! ocn velocity, meridional
  real(r8), allocatable ::  tocn (:)  ! ocean temperature
  real(r8), allocatable ::  zbot (:)  ! atm level height
  real(r8), allocatable ::  ubot (:)  ! atm velocity, zonal     
  real(r8), allocatable ::  vbot (:)  ! atm velocity, meridional
  real(r8), allocatable ::  thbot(:)  ! atm potential T
  real(r8), allocatable ::  shum (:)  ! atm specific humidity
  real(r8), allocatable ::  shum_16O (:)  ! atm H2O tracer
  real(r8), allocatable ::  shum_HDO (:)  ! atm HDO tracer
  real(r8), allocatable ::  shum_18O (:)  ! atm H218O tracer
  real(r8), allocatable ::  roce_16O (:)  ! ocn H2O ratio 
  real(r8), allocatable ::  roce_HDO (:)  ! ocn HDO ratio 
  real(r8), allocatable ::  roce_18O (:)  ! ocn H218O ratio 
  real(r8), allocatable ::  dens (:)  ! atm density
  real(r8), allocatable ::  tbot (:)  ! atm bottom surface T
  real(r8), allocatable ::  sen  (:)  ! heat flux: sensible 
  real(r8), allocatable ::  lat  (:)  ! heat flux: latent   
  real(r8), allocatable ::  lwup (:)  ! lwup over ocean
  real(r8), allocatable ::  evap (:)  ! water flux: evaporation
  real(r8), allocatable ::  evap_16O (:) !H2O flux: evaporation
  real(r8), allocatable ::  evap_HDO (:) !HDO flux: evaporation
  real(r8), allocatable ::  evap_18O (:) !H218O flux: evaporation
  real(r8), allocatable ::  taux (:)  ! wind stress, zonal
  real(r8), allocatable ::  tauy (:)  ! wind stress, meridional
  real(r8), allocatable ::  tref (:)  ! diagnostic:  2m ref T
  real(r8), allocatable ::  qref (:)  ! diagnostic:  2m ref Q
  real(r8), allocatable :: duu10n(:)  ! diagnostic: 10m wind speed squared

  real(r8), allocatable :: fswpen (:) ! fraction of sw penetrating ocn surface layer
  real(r8), allocatable :: ocnsal (:) ! ocean salinity
  real(r8), allocatable :: uGust  (:) ! wind gust
  real(r8), allocatable :: lwdn   (:) ! long  wave, downward
  real(r8), allocatable :: swdn   (:) ! short wave, downward
  real(r8), allocatable :: swup   (:) ! short wave, upward
  real(r8), allocatable :: prec   (:) ! precip
  real(r8), allocatable :: prec_gust (:) ! atm precip for convective gustiness (kg/m^3)

  ! Diurnal cycle variables wrt flux
 
  real(r8), allocatable :: tbulk      (:) ! diagnostic: ocn bulk T  
  real(r8), allocatable :: tskin      (:) ! diagnostic: ocn skin T  
  real(r8), allocatable :: tskin_night(:) ! diagnostic: ocn skin T  
  real(r8), allocatable :: tskin_day  (:) ! diagnostic: ocn skin T  
  real(r8), allocatable :: cSkin      (:) ! diagnostic: ocn cool skin  
  real(r8), allocatable :: cSkin_night(:) ! diagnostic: ocn cool skin  
  real(r8), allocatable :: warm       (:) ! diagnostic: ocn warming  
  real(r8), allocatable :: salt       (:) ! diagnostic: ocn salting  
  real(r8), allocatable :: speed      (:) ! diagnostic: ocn speed    
  real(r8), allocatable :: regime     (:) ! diagnostic: ocn regime   
  real(r8), allocatable :: warmMax    (:) ! diagnostic: ocn warming, max daily value
  real(r8), allocatable :: windMax    (:) ! diagnostic: ocn wind   , max daily value
  real(r8), allocatable :: QsolAvg    (:) ! diagnostic: ocn Qsol   , daily avg
  real(r8), allocatable :: windAvg    (:) ! diagnostic: ocn wind   , daily avg
  real(r8), allocatable :: warmMaxInc (:) ! diagnostic: ocn warming, max daily value, increment
  real(r8), allocatable :: windMaxInc (:) ! diagnostic: ocn wind   , max daily value, increment
  real(r8), allocatable :: qSolInc    (:) ! diagnostic: ocn Qsol   , daily avg, increment
  real(r8), allocatable :: windInc    (:) ! diagnostic: ocn wind   , daily avg, increment
  real(r8), allocatable :: nInc       (:) ! diagnostic: a/o flux   , increment

  real(r8), allocatable ::  ustar(:)  ! saved ustar
  real(r8), allocatable ::  re   (:)  ! saved re
  real(r8), allocatable ::  ssq  (:)  ! saved sq

  ! Conversion from degrees to radians

  real(r8),parameter :: const_pi      = SHR_CONST_PI       ! pi
  real(r8),parameter :: const_deg2rad = const_pi/180.0_r8  ! deg to rads

  ! Coupler field indices

  integer :: index_a2x_Sa_z    
  integer :: index_a2x_Sa_u    
  integer :: index_a2x_Sa_v    
  integer :: index_a2x_Sa_tbot 
  integer :: index_a2x_Sa_ptem 
  integer :: index_a2x_Sa_shum 
  integer :: index_a2x_Sa_shum_16O
  integer :: index_a2x_Sa_shum_HDO
  integer :: index_a2x_Sa_shum_18O
  integer :: index_a2x_Sa_dens 
  integer :: index_a2x_Faxa_swndr
  integer :: index_a2x_Faxa_swndf
  integer :: index_a2x_Faxa_swvdr
  integer :: index_a2x_Faxa_swvdf
  integer :: index_a2x_Faxa_lwdn
  integer :: index_a2x_Faxa_rainc
  integer :: index_a2x_Faxa_rainl
  integer :: index_a2x_Faxa_snowc
  integer :: index_a2x_Faxa_snowl
  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_o2x_So_fswpen
  integer :: index_o2x_So_s
  integer :: index_o2x_So_roce_16O
  integer :: index_o2x_So_roce_HDO
  integer :: index_o2x_So_roce_18O
  integer :: index_xao_So_tref    
  integer :: index_xao_So_qref    
  integer :: index_xao_So_avsdr   
  integer :: index_xao_So_avsdf   
  integer :: index_xao_So_anidr   
  integer :: index_xao_So_anidf   
  integer :: index_xao_Faox_taux  
  integer :: index_xao_Faox_tauy   
  integer :: index_xao_Faox_lat   
  integer :: index_xao_Faox_sen   
  integer :: index_xao_Faox_evap 
  integer :: index_xao_Faox_evap_16O
  integer :: index_xao_Faox_evap_HDO
  integer :: index_xao_Faox_evap_18O
  integer :: index_xao_Faox_lwup
  integer :: index_xao_Faox_swdn
  integer :: index_xao_Faox_swup
  integer :: index_xao_So_ustar
  integer :: index_xao_So_re   
  integer :: index_xao_So_ssq  
  integer :: index_xao_So_duu10n 
  integer :: index_xao_So_u10
  integer :: index_xao_So_fswpen
  integer :: index_xao_So_warm_diurn
  integer :: index_xao_So_salt_diurn
  integer :: index_xao_So_speed_diurn
  integer :: index_xao_So_regime_diurn
  integer :: index_xao_So_tskin_diurn
  integer :: index_xao_So_tskin_day_diurn
  integer :: index_xao_So_tskin_night_diurn
  integer :: index_xao_So_cskin_diurn
  integer :: index_xao_So_cskin_night_diurn
  integer :: index_xao_So_tbulk_diurn
  integer :: index_xao_So_warmmax_diurn
  integer :: index_xao_So_windmax_diurn
  integer :: index_xao_So_qsolavg_diurn
  integer :: index_xao_So_windavg_diurn
  integer :: index_xao_So_warmmaxinc_diurn
  integer :: index_xao_So_windmaxinc_diurn
  integer :: index_xao_So_qsolinc_diurn
  integer :: index_xao_So_windinc_diurn
  integer :: index_xao_So_ninc_diurn

  character(len=16) :: fluxsetting = 'unknown'
  character(len=*),parameter  :: fluxsetting_atmocn = 'atmocn'
  character(len=*),parameter  :: fluxsetting_exchange = 'exchange'

  !--- for exchange grid ---
  type(mct_rearr) :: Re_a2e, Re_e2a, Re_o2e, Re_e2o  ! atm/ocn/exch rearrangers
  type(mct_sMat ) :: sMata2o, sMato2a                ! decomp sMat 
  type(mct_gsMap) :: gsmap_ae, gsmap_oe              ! gsmaps for atm/ocn on exch grid
  integer(in)     :: nloc_a2o,nloc_o2a,nloc_o,nloc_a,nloc_ae,nloc_oe 

!===============================================================================
contains
!===============================================================================

  subroutine seq_flux_init_mct(comp, fractions)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(component_type), intent(in) :: comp
    type(mct_aVect), intent(in)  :: fractions
    !
    ! Local variables
    !
    type(mct_gsMap), pointer :: gsMap
    type(mct_gGrid), pointer :: dom
    integer(in)              :: nloc
    integer                  :: ko,ki     ! fractions indices
    integer                  :: ier
    real(r8), pointer        :: rmask(:)  ! ocn domain mask
    character(*),parameter   :: subName =   '(seq_flux_init_mct) '
    !-----------------------------------------------------------------------

    gsmap => component_get_gsmap_cx(comp) 
    dom   => component_get_dom_cx(comp) 

    nloc = mct_avect_lsize(dom%data)

    ! Input fields atm
    allocate( zbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate zbot',ier)
    zbot = 0.0_r8
    allocate( ubot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ubot',ier)
    ubot = 0.0_r8
    allocate( vbot(nloc))
    if(ier/=0) call mct_die(subName,'allocate vbot',ier)
    vbot = 0.0_r8
    allocate(thbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate thbot',ier)
    thbot = 0.0_r8
    allocate(shum(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum',ier)
    shum = 0.0_r8
    allocate(shum_16O(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum_16O',ier)
    shum_16O = 0.0_r8
    allocate(shum_HDO(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum_HDO',ier)
    shum_HDO = 0.0_r8
    allocate(shum_18O(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum_18O',ier)
    shum_18O = 0.0_r8
    allocate(dens(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate dens',ier)
    dens = 0.0_r8
    allocate(tbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tbot',ier)
    tbot = 0.0_r8
    allocate(ustar(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ustar',ier)
    ustar = 0.0_r8
    allocate(re(nloc), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate re',ier)
    re = 0.0_r8
    allocate(ssq(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ssq',ier)
    ssq = 0.0_r8
    allocate( uocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate uocn',ier)
    uocn = 0.0_r8
    allocate( vocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate vocn',ier)
    vocn = 0.0_r8
    allocate( tocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tocn',ier)
    tocn = 0.0_r8
    allocate(roce_16O(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate roce_16O',ier)
    roce_16O = 0.0_r8
    allocate(roce_HDO(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate roce_HDO',ier)
    roce_HDO = 0.0_r8
    allocate(roce_18O(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate roce_18O',ier)
    roce_18O = 0.0_r8

    ! Output fields 
    allocate(sen (nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sen',ier)
    sen  = 0.0_r8
    allocate(lat (nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lat',ier)
    lat  = 0.0_r8
    allocate(evap(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap',ier)
    evap = 0.0_r8
    allocate(evap_16O(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap_16O',ier)
    evap_16O = 0.0_r8
    allocate(evap_HDO(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap_HDO',ier)
    evap_HDO = 0.0_r8
    allocate(evap_18O(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap_18O',ier)
    evap_18O = 0.0_r8
    allocate(lwup(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lwup',ier)
    lwup = 0.0_r8
    allocate(taux(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate taux',ier)
    taux = 0.0_r8
    allocate(tauy(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tauy',ier)
    tauy = 0.0_r8
    allocate(tref(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tref',ier)
    tref = 0.0_r8
    allocate(qref(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate qref',ier)
    qref = 0.0_r8
    allocate(duu10n(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate duu10n',ier)
    duu10n = 0.0_r8

    !--- flux_diurnal cycle flux fields ---
    allocate(uGust(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate uGust',ier)
    uGust = 0.0_r8
    allocate(lwdn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lwdn',ier)
    lwdn = 0.0_r8
    allocate(swdn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate swdn',ier)
    swdn = 0.0_r8
    allocate(swup(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate swup',ier)
    swup = 0.0_r8
    allocate(prec(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate prec',ier)
    prec = 0.0_r8
    allocate(prec_gust(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate prec_gust',ier)
    prec_gust = 0.0_r8
    allocate(fswpen(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate fswpen',ier)
    fswpen = 0.0_r8
    allocate(ocnsal(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ocnsal',ier)
    ocnsal = 0.0_r8

    allocate(tbulk(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tbulk',ier)
    tbulk = 0.0_r8
    allocate(tskin(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tskin',ier)
    tskin = 0.0_r8
    allocate(tskin_day(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tskin_day',ier)
    tskin_day = 0.0_r8
    allocate(tskin_night(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tskin_night',ier)
    tskin_night = 0.0_r8
    allocate(cskin(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate cskin',ier)
    cskin = 0.0_r8
    allocate(cskin_night(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate cskin_night',ier)
    cskin_night = 0.0_r8

    allocate(warm(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate warm',ier)
    warm = 0.0_r8
    allocate(salt(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate salt',ier)
    salt = 0.0_r8
    allocate(speed(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate speed',ier)
    speed = 0.0_r8
    allocate(regime(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate regime',ier)
    regime = 0.0_r8
    allocate(warmMax(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate warmMax',ier)
    warmMax = 0.0_r8
    allocate(windMax(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate windMax',ier)
    windMax = 0.0_r8
    allocate(qSolAvg(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate qSolAvg',ier)
    qSolAvg = 0.0_r8
    allocate(windAvg(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate windAvg',ier)
    windAvg = 0.0_r8

    allocate(warmMaxInc(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate warmMaxInc',ier)
    warmMaxInc = 0.0_r8
    allocate(windMaxInc(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate windMaxInc',ier)
    windMaxInc = 0.0_r8
    allocate(qSolInc(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate qSolInc',ier)
    qSolInc = 0.0_r8
    allocate(windInc(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate windInc',ier)
    windInc = 0.0_r8
    allocate(nInc     (nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate nInc',ier)
    nInc = 0.0_r8

    ! Grid fields
    allocate( lats(nloc),stat=ier )
    if(ier/=0) call mct_die(subName,'allocate lats',ier)
    lats = 0.0_r8
    allocate( lons(nloc),stat=ier )
    if(ier/=0) call mct_die(subName,'allocate lons',ier)
    lons = 0.0_r8
    allocate( emask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate emask',ier)
    emask = 0.0_r8
    allocate(mask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate mask',ier)
    mask = 0.0_r8
    
    ! Get lat, lon, mask, which is time-invariant
    allocate(rmask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate rmask',ier)
    call mct_gGrid_exportRAttr(dom, 'lat' , lats , nloc) 
    call mct_gGrid_exportRAttr(dom, 'lon' , lons , nloc) 

    ! setup the compute mask.
    ! prefer to compute just where ocean exists, so setup a mask here.
    ! this could be run with either the ocean or atm grid so need to be careful.
    ! really want the ocean mask on ocean grid or ocean mask mapped to atm grid,
    ! but do not have access to the ocean mask mapped to the atm grid.
    ! the dom mask is a good place to start, on ocean grid, it should be what we want,
    ! on the atm grid, it's just all 1's so not very useful.
    ! next look at ofrac+ifrac in fractions.  want to compute on all non-land points.
    ! using ofrac alone will exclude points that are currently all sea ice but that later
    ! could be less that 100% covered in ice.

    ! default compute everywhere, then "turn off" gridcells
    mask = 1
   
    ! use domain mask first
    call mct_gGrid_exportRAttr(dom, 'mask', rmask, nloc)
    where (rmask < 0.5_r8) mask = 0   ! like nint
    deallocate(rmask)

    ! then check ofrac + ifrac
    ko = mct_aVect_indexRA(fractions,"ofrac")
    ki = mct_aVect_indexRA(fractions,"ifrac")
    where (fractions%rAttr(ko,:)+fractions%rAttr(ki,:) <= 0.0_r8) mask(:) = 0

    emask = mask

    fluxsetting = trim(fluxsetting_atmocn)

  end subroutine seq_flux_init_mct

!===============================================================================

  subroutine seq_flux_initexch_mct(atm, ocn, mpicom_cplid, cplid)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(component_type), intent(in) :: atm 
    type(component_type), intent(in) :: ocn
    integer(in)         , intent(in) :: mpicom_cplid
    integer(in)         , intent(in) :: cplid
    !
    ! Local variables
    !
    type(mct_gsMap), pointer :: gsmap_a
    type(mct_gGrid), pointer :: dom_a
    type(mct_gsMap), pointer :: gsmap_o
    type(mct_gGrid), pointer :: dom_o
    integer(in)              :: kw,ka,ko,iw,ia,io,n
    character(len=128)       :: strat
    integer                  :: ier
    integer                  :: mytask
    integer(in)              :: kmsk            ! field indices
    character(len=128)       :: ConfigFileName  ! config file to read
    character(len=128)       :: MapLabel        ! map name
    character(len=128)       :: MapTypeLabel    ! map type
    character(len=256)       :: fileName
    character(len=1)         :: maptype
    character(len=3)         :: Smaptype
    type(mct_aVect)          :: avdom_oe
    type(mct_list)           :: sort_keys
    character(*),parameter :: subName =   '(seq_flux_initexch_mct) '
    !-----------------------------------------------------------------------

    gsmap_a => component_get_gsmap_cx(atm) ! gsmap_ax
    gsmap_o => component_get_gsmap_cx(ocn) ! gsmap_ox
    dom_a   => component_get_dom_cx(atm)   ! dom_ax
    dom_o   => component_get_dom_cx(ocn)   ! dom_ox

    call shr_mpi_commrank(mpicom_cplid, mytask)

    !--- Get mapping file info
    do n = 1,2
       ConfigFileName = "seq_maps.rc"
       if (n == 1) then
          MapLabel = "atm2ocn_fmapname:"
          MapTypeLabel = "atm2ocn_fmaptype:"
       elseif (n == 2) then
          MapLabel = "ocn2atm_fmapname:"
          MapTypeLabel = "ocn2atm_fmaptype:"
       else
          call shr_sys_abort(trim(subname)//' do error1')
       endif

       call shr_mct_queryConfigFile(mpicom_cplid, ConfigFilename, &
          trim(MapLabel),fileName,trim(MapTypeLabel),maptype)

       !--- hardwire decomposition to gsmap_o
       if (n == 1) then
          Smaptype = "src"
          call shr_mct_sMatReaddnc(sMata2o, gsmap_a, gsmap_o, Smaptype, &
             filename=fileName, mytask=mytask, mpicom=mpicom_cplid)
       elseif (n == 2) then
          Smaptype = "dst"
          call shr_mct_sMatReaddnc(sMato2a, gsmap_o, gsmap_a, Smaptype, &
             filename=fileName, mytask=mytask, mpicom=mpicom_cplid)
       else
          call shr_sys_abort(trim(subname)//' do error2')
       endif

    enddo

    !--- the two mapping files must have their local indices in identical order
    !--- sort the global indices as a starting point

    call mct_list_init(sort_keys,'grow:gcol')
    call mct_sMat_SortPermute(sMata2o,sort_keys)
    call mct_list_clean(sort_keys)
    call mct_list_init(sort_keys,'gcol:grow')
    call mct_sMat_SortPermute(sMato2a,sort_keys)
    call mct_list_clean(sort_keys)

    !--- now check that they are sorted properly

    nloc_a2o= mct_sMat_lsize(sMata2o)
    nloc_o2a= mct_sMat_lsize(sMato2a)

    if (nloc_a2o /= nloc_o2a) then
       write(logunit,*) trim(subname),' ERROR: sMat sizes',nloc_a2o,nloc_o2a
       call shr_sys_abort(trim(subname)//' ERROR in sMat sizes')
    endif
    ko = mct_sMat_indexIA(sMata2o,'grow')    ! local row (dst) index
    ka = mct_sMat_indexIA(sMato2a,'gcol')    ! local column (src) index
    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       ia = sMato2a%data%iAttr(ka,n)
       if (io /= ia) then
          write(logunit,*) trim(subname),' ERROR: sMat indices1 ',io,ia
          call shr_sys_abort(trim(subname)//' ERROR in sMat indices1')
       endif
    enddo
    ko = mct_sMat_indexIA(sMata2o,'gcol')    ! local column (src) index
    ka = mct_sMat_indexIA(sMato2a,'grow')    ! local row (dst) index
    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       ia = sMato2a%data%iAttr(ka,n)
       if (io /= ia) then
          write(logunit,*) trim(subname),' ERROR: sMat indices2 ',io,ia
          call shr_sys_abort(trim(subname)//' ERROR in sMat indices2')
       endif
    enddo

    !--- instantiate/create/compute various datatypes

    call mct_sMat_2XgsMap(sMata2o , gsmap_ae, 0, mpicom_cplid, cplid)
    call mct_sMat_2YgsMap(sMata2o , gsmap_oe, 0, mpicom_cplid, cplid)

    call mct_rearr_init(gsmap_a   , gsmap_ae, mpicom_cplid, Re_a2e)
    call mct_rearr_init(gsmap_ae  , gsmap_a,  mpicom_cplid, Re_e2a)
    call mct_rearr_init(gsmap_o   , gsmap_oe, mpicom_cplid, Re_o2e)
    call mct_rearr_init(gsmap_oe  , gsmap_o,  mpicom_cplid, Re_e2o)

    call mct_sMat_g2lMat(sMata2o  , gsmap_ae, 'column',mpicom_cplid)
    call mct_sMat_g2lMat(sMata2o  , gsmap_oe, 'row',   mpicom_cplid)
    call mct_sMat_g2lMat(sMato2a  , gsmap_ae, 'row',   mpicom_cplid)
    call mct_sMat_g2lMat(sMato2a  , gsmap_oe, 'column',mpicom_cplid)

    nloc_a  = mct_gsmap_lsize(gsmap_a  , mpicom_cplid)
    nloc_o  = mct_gsmap_lsize(gsmap_o  , mpicom_cplid)
    nloc_ae = mct_gsmap_lsize(gsmap_ae , mpicom_cplid)
    nloc_oe = mct_gsmap_lsize(gsmap_oe , mpicom_cplid)

    call mct_gsmap_clean(gsmap_ae)
    call mct_gsmap_clean(gsmap_oe)

    ! Input fields atm
    allocate( emask(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate emask',ier)
    allocate( zbot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate zbot',ier)
    allocate( ubot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ubot',ier)
    allocate( vbot(nloc_a2o))
    if(ier/=0) call mct_die(subName,'allocate vbot',ier)
    allocate(thbot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate thbot',ier)
    allocate(shum(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum',ier)
    allocate(shum_16O(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum_16O',ier)
    allocate(shum_HDO(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum_HDO',ier)
    allocate(shum_18O(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum_18O',ier)
    allocate(dens(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate dens',ier)
    allocate(tbot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tbot',ier)
    allocate(ustar(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ustar',ier)
    allocate(re(nloc_a2o), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate re',ier)
    allocate(ssq(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ssq',ier)
    allocate( uocn(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate uocn',ier)
    allocate( vocn(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate vocn',ier)
    allocate( tocn(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tocn',ier)

    ! Output fields 
    allocate(sen (nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sen',ier)
    allocate(lat (nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lat',ier)
    allocate(evap(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap',ier)
    allocate(evap_16O(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap_16O',ier)
    allocate(evap_HDO(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap_HDO',ier)
    allocate(evap_18O(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap_18O',ier)
    allocate(lwup(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lwup',ier)
    allocate(taux(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate taux',ier)
    allocate(tauy(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tauy',ier)
    allocate(tref(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tref',ier)
    allocate(qref(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate qref',ier)
    allocate(duu10n(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate duu10n',ier)

    ! set emask

    call mct_avect_init(avdom_oe,dom_o%data,lsize=nloc_oe)
    call mct_rearr_rearrange(dom_o%data, avdom_oe, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
    ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local dst index
    kmsk = mct_aVect_indexRA(avdom_oe,"mask",dieWith=subName)
    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       emask(n) = nint(avdom_oe%rAttr(kmsk,io))
       if (emask(n) == 0) then
          write(logunit,*) trim(subname),' ERROR: weights use masked ocean value'
          call shr_sys_abort(trim(subname)//' ERROR: weights use masked ocean value')
       endif
    enddo

    call mct_aVect_clean(avdom_oe)

    fluxsetting = trim(fluxsetting_exchange)

  end subroutine seq_flux_initexch_mct

!===============================================================================

  subroutine seq_flux_ocnalb_mct( infodata, ocn, a2x_o, fractions_o, xao_o )

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_infodata_type) , intent(in)    :: infodata
    type(component_type)    , intent(in)    :: ocn
    type(mct_aVect)         , intent(in)    :: a2x_o
    type(mct_aVect)         , intent(inout) :: fractions_o
    type(mct_aVect)         , intent(inout) :: xao_o
    !
    ! Local variables
    !
    type(mct_gGrid), pointer :: dom_o
    logical		:: flux_albav		! flux avg option
    integer(in)		:: n,i			! indices
    real(r8)		:: rlat			! gridcell latitude in radians
    real(r8)		:: rlon			! gridcell longitude in radians
    real(r8)		:: cosz			! Cosine of solar zenith angle
    real(r8)		:: eccen		! Earth orbit eccentricity
    real(r8)		:: mvelpp		! Earth orbit
    real(r8)		:: lambm0		! Earth orbit
    real(r8)		:: obliqr		! Earth orbit
    real(r8)		:: delta		! Solar declination angle  in radians
    real(r8)		:: eccf			! Earth orbit eccentricity factor
    real(r8)		:: calday		! calendar day including fraction, at 0e
    real(r8)		:: nextsw_cday		! calendar day of next atm shortwave
    real(r8)		:: anidr		! albedo: near infrared, direct
    real(r8)		:: avsdr		! albedo: visible      , direct
    real(r8)		:: anidf		! albedo: near infrared, diffuse
    real(r8)		:: avsdf		! albedo: visible      , diffuse
    real(r8)		:: swdnc		! temporary swdn
    real(r8)		:: swupc		! temporary swup
    integer(in)		:: ID			! comm ID
    integer(in)		:: ier			! error code
    integer(in)		:: kx,kr		! fractions indices
    integer(in)		:: klat,klon,kmsk	! field indices
    logical		:: update_alb		! was albedo updated
    logical,save	:: first_call = .true. 
    !
    real(r8),parameter :: albdif = 0.06_r8 ! 60 deg reference albedo, diffuse
    real(r8),parameter :: albdir = 0.07_r8 ! 60 deg reference albedo, direct 
    character(*),parameter :: subName =   '(seq_flux_ocnalb_mct) '
    !
    !-----------------------------------------------------------------------

    dom_o => component_get_dom_cx(ocn) ! dom_ox

    call seq_infodata_getData(infodata , &
         flux_albav=flux_albav)

    ! Determine indices

    update_alb = .false.

    if (first_call) then
       index_xao_So_anidr  = mct_aVect_indexRA(xao_o,'So_anidr')
       index_xao_So_anidf  = mct_aVect_indexRA(xao_o,'So_anidf')
       index_xao_So_avsdr  = mct_aVect_indexRA(xao_o,'So_avsdr')
       index_xao_So_avsdf  = mct_aVect_indexRA(xao_o,'So_avsdf')
       index_xao_Faox_swdn = mct_aVect_indexRA(xao_o,'Faox_swdn')
       index_xao_Faox_swup = mct_aVect_indexRA(xao_o,'Faox_swup')

       index_a2x_Faxa_swndr = mct_aVect_indexRA(a2x_o,'Faxa_swndr')
       index_a2x_Faxa_swndf = mct_aVect_indexRA(a2x_o,'Faxa_swndf')
       index_a2x_Faxa_swvdr = mct_aVect_indexRA(a2x_o,'Faxa_swvdr')
       index_a2x_Faxa_swvdf = mct_aVect_indexRA(a2x_o,'Faxa_swvdf')

       nloc_o  = mct_ggrid_lsize(dom_o)
       klat = mct_gGrid_indexRA(dom_o,"lat" ,dieWith=subName)
       klon = mct_gGrid_indexRA(dom_o,"lon" ,dieWith=subName)
       allocate( lats(nloc_o),stat=ier )
       if(ier/=0) call mct_die(subName,'allocate lats',ier)
       allocate( lons(nloc_o),stat=ier )
       if(ier/=0) call mct_die(subName,'allocate lons',ier)
       do n = 1,nloc_o
          lats(n) = dom_o%data%rAttr(klat,n)
          lons(n) = dom_o%data%rAttr(klon,n)
       enddo
       first_call = .false.
    endif

    if (flux_albav) then

       do n=1,nloc_o   
          anidr = albdir
          avsdr = albdir
          anidf = albdif
          avsdf = albdif

          ! Albedo is now function of latitude (will be new implementation)
          !rlat = const_deg2rad * lats(n)
          !anidr = 0.069_r8 - 0.011_r8 * cos(2._r8 * rlat)
          !avsdr = anidr
          !anidf = anidr
          !avsdf = anidr

          xao_o%rAttr(index_xao_So_avsdr,n) = avsdr
          xao_o%rAttr(index_xao_So_anidr,n) = anidr
          xao_o%rAttr(index_xao_So_avsdf,n) = avsdf
          xao_o%rAttr(index_xao_So_anidf,n) = anidf
       end do
       update_alb = .true.

    else

       !--- flux_atmocn needs swdn & swup = swdn*(-albedo)
       !--- swdn & albedos are time-aligned  BEFORE albedos get updated below ---
       do n=1,nloc_o
          avsdr = xao_o%rAttr(index_xao_So_avsdr,n)
          anidr = xao_o%rAttr(index_xao_So_anidr,n)
          avsdf = xao_o%rAttr(index_xao_So_avsdf,n)
          anidf = xao_o%rAttr(index_xao_So_anidf,n)
          swupc = a2x_o%rAttr(index_a2x_Faxa_swndr,n)*(-anidr) &
              & + a2x_o%rAttr(index_a2x_Faxa_swndf,n)*(-anidf) &
              & + a2x_o%rAttr(index_a2x_Faxa_swvdr,n)*(-avsdr) &
              & + a2x_o%rAttr(index_a2x_Faxa_swvdf,n)*(-avsdf)
          swdnc = a2x_o%rAttr(index_a2x_Faxa_swndr,n) &
              & + a2x_o%rAttr(index_a2x_Faxa_swndf,n) &
              & + a2x_o%rAttr(index_a2x_Faxa_swvdr,n) &
              & + a2x_o%rAttr(index_a2x_Faxa_swvdf,n)
          if ( anidr == 1.0_r8 ) then ! dark side of earth
             swupc = 0.0_r8
             swdnc = 0.0_r8
          end if
          xao_o%rAttr(index_xao_Faox_swdn,n) = swdnc
          xao_o%rAttr(index_xao_Faox_swup,n) = swupc
       end do

       ! Solar declination 
       ! Will only do albedo calculation if nextsw_cday is not -1.
       
       call seq_infodata_GetData(infodata,nextsw_cday=nextsw_cday,orb_eccen=eccen, &
          orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr)
       if (nextsw_cday >= -0.5_r8) then
          calday = nextsw_cday
          call shr_orb_decl(calday, eccen, mvelpp,lambm0, obliqr, delta, eccf)
          ! Compute albedos 
          do n=1,nloc_o
             rlat = const_deg2rad * lats(n)
             rlon = const_deg2rad * lons(n)
             cosz = shr_orb_cosz( calday, rlat, rlon, delta )
             if (cosz  >  0.0_r8) then !--- sun hit --
                anidr = (.026_r8/(cosz**1.7_r8 + 0.065_r8)) +   &
                        (.150_r8*(cosz         - 0.100_r8 ) *   &
                                 (cosz         - 0.500_r8 ) *   &
                                 (cosz         - 1.000_r8 )  )
                avsdr = anidr
                anidf = albdif
                avsdf = albdif
             else !--- dark side of earth ---
                anidr = 1.0_r8
                avsdr = 1.0_r8
                anidf = 1.0_r8
                avsdf = 1.0_r8
             end if

             xao_o%rAttr(index_xao_So_avsdr,n) = avsdr
             xao_o%rAttr(index_xao_So_anidr,n) = anidr
             xao_o%rAttr(index_xao_So_avsdf,n) = avsdf
             xao_o%rAttr(index_xao_So_anidf,n) = anidf

          end do   ! nloc_o
          update_alb = .true.
       endif    ! nextsw_cday
    end if   ! flux_albav

    !--- update current ifrad/ofrad values if albedo was updated

    if (update_alb) then
       kx = mct_aVect_indexRA(fractions_o,"ifrac")
       kr = mct_aVect_indexRA(fractions_o,"ifrad")
       fractions_o%rAttr(kr,:) = fractions_o%rAttr(kx,:)
       kx = mct_aVect_indexRA(fractions_o,"ofrac")
       kr = mct_aVect_indexRA(fractions_o,"ofrad")
       fractions_o%rAttr(kr,:) = fractions_o%rAttr(kx,:)
    endif
       
    end subroutine seq_flux_ocnalb_mct

!===============================================================================

  subroutine seq_flux_atmocnexch_mct( infodata, atm, ocn, fractions_a, fractions_o, &
       xao_a, xao_o)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_infodata_type) , intent(in)    :: infodata
    type(component_type)    , intent(in)    :: atm
    type(component_type)    , intent(in)    :: ocn
    type(mct_aVect)         , intent(in)    :: fractions_a
    type(mct_aVect)         , intent(in)    :: fractions_o
    type(mct_aVect)         , intent(inout) :: xao_a
    type(mct_aVect)         , intent(inout) :: xao_o
    !
    ! Local variables
    !
    type(mct_aVect) , pointer :: a2x
    type(mct_aVect) , pointer :: o2x
    type(mct_gsmap) , pointer :: gsmap_a
    type(mct_gsmap) , pointer :: gsmap_o

    type(mct_aVect) :: a2x_e
    type(mct_aVect) :: o2x_e
    type(mct_aVect) :: xaop_ae
    type(mct_aVect) :: xaop_oe
    type(mct_aVect) :: xaop_a
    type(mct_aVect) :: xaop_o
    type(mct_aVect) :: fractions_oe

    integer(in) :: kw,ka,ko,iw,ia,io,kf
    integer(in) :: n,i          ! indices
    logical     :: dead_comps   ! .true.  => dead components are used
    integer(in) :: index_tref  
    integer(in) :: index_qref  
    integer(in) :: index_duu10n
    integer(in) :: index_ustar 
    integer(in) :: index_ssq   
    integer(in) :: index_re    
    integer(in) :: index_u10   
    integer(in) :: index_taux  
    integer(in) :: index_tauy  
    integer(in) :: index_lat   
    integer(in) :: index_sen   
    integer(in) :: index_evap  
    integer(in) :: index_evap_16O
    integer(in) :: index_evap_HDO
    integer(in) :: index_evap_18O
    integer(in) :: index_lwup  
    integer(in) :: index_sumwt
    integer(in) :: atm_nx,atm_ny,ocn_nx,ocn_ny
    real(r8)    :: wt
    real(r8)    :: gust_fac = huge(1.0_r8) !wind gust factor
    integer(in) :: tod, dt
    logical,save:: first_call = .true. 
    logical     :: read_restart    ! .true. => model starting from restart
    logical     :: ocn_prognostic  ! .true. => ocn is prognostic 
    logical     :: flux_diurnal    ! .true. => turn on diurnal cycle in atm/ocn fluxes
    logical     :: cold_start      ! .true. to initialize internal fields in shr_flux diurnal
    character(len=256) :: fldlist  ! subset of xao fields
    !
    character(*),parameter :: subName =   '(seq_flux_atmocnexch_mct) '
    !
    !-----------------------------------------------------------------------

    gsmap_a => component_get_gsmap_cx(atm) 
    gsmap_o => component_get_gsmap_cx(ocn) 
    a2x     => component_get_c2x_cx(atm)  ! a2x_ax
    o2x     => component_get_c2x_cx(ocn)  ! o2x_ox 

    if (trim(fluxsetting) /= trim(fluxsetting_exchange)) then
       call shr_sys_abort(trim(subname)//' ERROR wrong fluxsetting')
    endif

    ! Update ocean surface fluxes 
    ! Must fabricate "reasonable" data (using dead components)

    call seq_infodata_GetData(infodata, &
         read_restart=read_restart, &
         dead_comps=dead_comps,         &
         atm_nx=atm_nx, atm_ny=atm_ny,  &
         ocn_nx=ocn_nx, ocn_ny=ocn_ny,  &
         ocn_prognostic=ocn_prognostic, &
         flux_diurnal=flux_diurnal,     &
         gust_fac = gust_fac            )

    cold_start = .false.   ! use restart data or data from last timestep

    if (first_call) then
       if (.not.read_restart) cold_start = .true.
       first_call = .false. 
    endif

    if (dead_comps) then
       do n = 1,nloc_a2o
          tocn(n) = 290.0_r8 ! ocn temperature            ~ Kelvin
          uocn(n) =   0.0_r8 ! ocn velocity, zonal        ~ m/s
          vocn(n) =   0.0_r8 ! ocn velocity, meridional   ~ m/s
          zbot(n) =  55.0_r8 ! atm height of bottom layer ~ m
          ubot(n) =   0.0_r8 ! atm velocity, zonal        ~ m/s
          vbot(n) =   2.0_r8 ! atm velocity, meridional   ~ m/s
          thbot(n)= 301.0_r8 ! atm potential temperature  ~ Kelvin
          shum(n) = 1.e-2_r8 ! atm specific humidity      ~ kg/kg
          shum_16O(n) = 1.e-2_r8 ! H216O specific humidity    ~ kg/kg
          shum_HDO(n) = 1.e-2_r8 ! HD16O specificy humidity   ~ kg/kg 
          shum_18O(n) = 1.e-2_r8 ! H218O specific humidity    ~ kg/kg
          roce_16O(n) = 1.0_r8 ! H216O ratio ~ mol/mol
          roce_HDO(n) = 1.0_r8 ! HD16O ratio ~ mol/mol 
          roce_18O(n) = 1.0_r8 ! H218O ratio ~ mol/mol
          dens(n) =   1.0_r8 ! atm density                ~ kg/m^3
          tbot(n) = 300.0_r8 ! atm temperature            ~ Kelvin
       enddo
    else        

       !--- instantiate exchange grid aVects
       call mct_AVect_init(a2x_e, a2x, nloc_ae)
       call mct_AVect_zero(a2x_e)
       call mct_AVect_init(o2x_e, o2x, nloc_oe)
       call mct_AVect_zero(o2x_e)

       !--- rearrange a2x and o2x into exchange grid

       call mct_rearr_rearrange(a2x, a2x_e, Re_a2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
       call mct_rearr_rearrange(o2x, o2x_e, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)

       !--- extract fields from a2x and o2x (_e) into local arrays on exchange grid

       ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local row index
       ka = mct_sMat_indexIA(sMata2o,'lcol')    ! local column index

       do n = 1,nloc_a2o
          io = sMata2o%data%iAttr(ko,n)
          ia = sMata2o%data%iAttr(ka,n)
          zbot(n) = a2x_e%rAttr(index_a2x_Sa_z   ,ia)
          ubot(n) = a2x_e%rAttr(index_a2x_Sa_u   ,ia)
          vbot(n) = a2x_e%rAttr(index_a2x_Sa_v   ,ia)
          thbot(n)= a2x_e%rAttr(index_a2x_Sa_ptem,ia)
          shum(n) = a2x_e%rAttr(index_a2x_Sa_shum,ia)
          shum_16O(n) = a2x_e%rAttr(index_a2x_Sa_shum_16O,ia)
          shum_HDO(n) = a2x_e%rAttr(index_a2x_Sa_shum_HDO,ia)
          shum_18O(n) = a2x_e%rAttr(index_a2x_Sa_shum_18O,ia)
          dens(n) = a2x_e%rAttr(index_a2x_Sa_dens,ia)
          tbot(n) = a2x_e%rAttr(index_a2x_Sa_tbot,ia)
          tocn(n) = o2x_e%rAttr(index_o2x_So_t   ,io)   
          uocn(n) = o2x_e%rAttr(index_o2x_So_u   ,io)
          vocn(n) = o2x_e%rAttr(index_o2x_So_v   ,io)
          roce_16O(n) = o2x_e%rAttr(index_o2x_So_roce_16O, io)
          roce_HDO(n) = o2x_e%rAttr(index_o2x_So_roce_HDO, io)
          roce_18O(n) = o2x_e%rAttr(index_o2x_So_roce_18O, io)
       enddo
       call mct_aVect_clean(a2x_e)
       call mct_aVect_clean(o2x_e)
    end if

    if (flux_diurnal) then
       call shr_flux_atmocn_diurnal (nloc_a2o , zbot , ubot, vbot, thbot, &
                          shum , shum_16O , shum_HDO, shum_18O, dens , tbot, uocn, vocn , &
                          tocn , emask, sen , lat , lwup , &
                          roce_16O, roce_HDO, roce_18O,    &
                          evap , evap_16O, evap_HDO, evap_18O, taux , tauy, tref, qref , &
                          uGust, lwdn , swdn , swup, prec, &
                          fswpen, ocnsal, ocn_prognostic, flux_diurnal,   &
                          lats , lons , warm , salt , speed, regime,      &
                          warmMax, windMax, qSolAvg, windAvg,             &
                          warmMaxInc, windMaxInc, qSolInc, windInc, nInc, &
                          tbulk, tskin, tskin_day, tskin_night, &
                          cskin, cskin_night, tod, dt,          &
                          duu10n,ustar, re  , ssq , missval = 0.0_r8, &
                          cold_start=cold_start)
    else
       call shr_flux_atmocn (nloc_a2o , zbot , ubot, vbot, thbot, prec_gust, gust_fac, &
                          shum , shum_16O , shum_HDO, shum_18O, dens , tbot, uocn, vocn , &
                          tocn , emask, sen , lat , lwup , &
                          roce_16O, roce_HDO, roce_18O,    & 
                          evap , evap_16O, evap_HDO, evap_18O, taux, tauy, tref, qref , &
                          duu10n,ustar, re  , ssq , missval = 0.0_r8 )
    endif

    !--- create temporary aVects on exchange, atm, or ocn decomp as needed

    fldlist = trim(seq_flds_xao_states)//":"//trim(seq_flds_xao_fluxes)//":sumwt"
    call mct_aVect_init(xaop_ae,rList=trim(fldlist),lsize=nloc_ae)
    call mct_aVect_zero(xaop_ae)
    call mct_aVect_init(xaop_oe,rList=trim(fldlist),lsize=nloc_oe)
    call mct_aVect_zero(xaop_oe)
    call mct_aVect_init(xaop_a, rList=trim(fldlist),lsize=nloc_a)
    call mct_aVect_zero(xaop_a)
    call mct_aVect_init(xaop_o, rList=trim(fldlist),lsize=nloc_o)
    call mct_aVect_zero(xaop_o)

    index_tref   = mct_aVect_indexRA(xaop_ae,"So_tref")
    index_qref   = mct_aVect_indexRA(xaop_ae,"So_qref")
    index_duu10n = mct_aVect_indexRA(xaop_ae,"So_duu10n")
    index_ustar  = mct_aVect_indexRA(xaop_ae,"So_ustar")
    index_ssq    = mct_aVect_indexRA(xaop_ae,"So_ssq")
    index_re     = mct_aVect_indexRA(xaop_ae,"So_re")
    index_u10    = mct_aVect_indexRA(xaop_ae,"So_u10")
    index_taux   = mct_aVect_indexRA(xaop_ae,"Faox_taux")
    index_tauy   = mct_aVect_indexRA(xaop_ae,"Faox_tauy")
    index_lat    = mct_aVect_indexRA(xaop_ae,"Faox_lat")
    index_sen    = mct_aVect_indexRA(xaop_ae,"Faox_sen")
    index_evap   = mct_aVect_indexRA(xaop_ae,"Faox_evap")
    index_evap_16O = mct_aVect_indexRA(xaop_ae,"Faox_evap_16O", perrWith='quiet')
    index_evap_HDO = mct_aVect_indexRA(xaop_ae,"Faox_evap_HDO", perrWith='quiet')
    index_evap_18O = mct_aVect_indexRA(xaop_ae,"Faox_evap_18O", perrWith='quiet')
    index_lwup   = mct_aVect_indexRA(xaop_ae,"Faox_lwup")
    index_sumwt  = mct_aVect_indexRA(xaop_ae,"sumwt")

    !--- aggregate ocean values locally based on exchange grid decomp

    ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local row index
    ka = mct_sMat_indexIA(sMata2o,'lcol')    ! local column index
    kw = mct_sMat_indexRA(sMata2o,'weight')  ! weight index

    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       ia = sMata2o%data%iAttr(ka,n)
       wt = sMata2o%data%rAttr(kw,n)
       xaop_oe%rAttr(index_sen   ,io) = xaop_oe%rAttr(index_sen   ,io) + sen(n) * wt
       xaop_oe%rAttr(index_lat   ,io) = xaop_oe%rAttr(index_lat   ,io) + lat(n) * wt
       xaop_oe%rAttr(index_taux  ,io) = xaop_oe%rAttr(index_taux  ,io) + taux(n)* wt
       xaop_oe%rAttr(index_tauy  ,io) = xaop_oe%rAttr(index_tauy  ,io) + tauy(n)* wt
       xaop_oe%rAttr(index_evap  ,io) = xaop_oe%rAttr(index_evap  ,io) + evap(n)* wt
       if ( index_evap_16O /= 0 ) xaop_oe%rAttr(index_evap_16O ,io) = xaop_oe%rAttr(index_evap_16O  ,io) + evap_16O(n)* wt
       if ( index_evap_HDO /= 0 ) xaop_oe%rAttr(index_evap_HDO ,io) = xaop_oe%rAttr(index_evap_HDO  ,io) + evap_HDO(n)* wt
       if ( index_evap_18O /= 0 ) xaop_oe%rAttr(index_evap_18O ,io) = xaop_oe%rAttr(index_evap_18O  ,io) + evap_18O(n)* wt
       xaop_oe%rAttr(index_tref  ,io) = xaop_oe%rAttr(index_tref  ,io) + tref(n)* wt
       xaop_oe%rAttr(index_qref  ,io) = xaop_oe%rAttr(index_qref  ,io) + qref(n)* wt
       xaop_oe%rAttr(index_ustar ,io) = xaop_oe%rAttr(index_ustar ,io) + ustar(n)*wt   ! friction velocity
       xaop_oe%rAttr(index_re    ,io) = xaop_oe%rAttr(index_re    ,io) + re(n)  * wt   ! reynolds number
       xaop_oe%rAttr(index_ssq   ,io) = xaop_oe%rAttr(index_ssq   ,io) + ssq(n) * wt   ! s.hum. saturation at Ts
       xaop_oe%rAttr(index_lwup  ,io) = xaop_oe%rAttr(index_lwup  ,io) + lwup(n)* wt   
       xaop_oe%rAttr(index_duu10n,io) = xaop_oe%rAttr(index_duu10n,io) + duu10n(n)*wt  
       xaop_oe%rAttr(index_u10   ,io) = xaop_oe%rAttr(index_u10   ,io) + sqrt(duu10n(n))*wt
       xaop_oe%rAttr(index_sumwt ,io) = xaop_oe%rAttr(index_sumwt ,io) + wt
    enddo

    !--- aggregate atm values locally based on exchange grid decomp

    ko = mct_sMat_indexIA(sMato2a,'lcol')    ! local column index
    ka = mct_sMat_indexIA(sMato2a,'lrow')    ! local row index
    kw = mct_sMat_indexRA(sMato2a,'weight')  ! weight index
    kf = mct_aVect_indexRA(fractions_o,"ofrac")

    !--- to apply fraction corrections, the indexing must be correct so rearrange
    call mct_avect_init(fractions_oe,fractions_o,lsize=nloc_oe)
    call mct_rearr_rearrange(fractions_o, fractions_oe, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
    do n = 1,nloc_o2a
       io = sMato2a%data%iAttr(ko,n)
       ia = sMato2a%data%iAttr(ka,n)
!tcx   wt = sMato2a%data%rAttr(kw,n)
       wt = sMato2a%data%rAttr(kw,n) * fractions_oe%rAttr(kf,io)
       xaop_ae%rAttr(index_sen   ,ia) = xaop_ae%rAttr(index_sen   ,ia) + sen(n) * wt
       xaop_ae%rAttr(index_lat   ,ia) = xaop_ae%rAttr(index_lat   ,ia) + lat(n) * wt
       xaop_ae%rAttr(index_taux  ,ia) = xaop_ae%rAttr(index_taux  ,ia) + taux(n)* wt
       xaop_ae%rAttr(index_tauy  ,ia) = xaop_ae%rAttr(index_tauy  ,ia) + tauy(n)* wt
       xaop_ae%rAttr(index_evap  ,ia) = xaop_ae%rAttr(index_evap  ,ia) + evap(n)* wt
       if ( index_evap_16O /= 0 ) xaop_ae%rAttr(index_evap_16O ,ia) = xaop_ae%rAttr(index_evap_16O  ,ia) + evap_16O(n)* wt
       if ( index_evap_HDO /= 0 ) xaop_ae%rAttr(index_evap_HDO ,ia) = xaop_ae%rAttr(index_evap_HDO  ,ia) + evap_HDO(n)* wt
       if ( index_evap_18O /= 0 ) xaop_ae%rAttr(index_evap_18O ,ia) = xaop_ae%rAttr(index_evap_18O  ,ia) + evap_18O(n)* wt
       xaop_ae%rAttr(index_tref  ,ia) = xaop_ae%rAttr(index_tref  ,ia) + tref(n)* wt
       xaop_ae%rAttr(index_qref  ,ia) = xaop_ae%rAttr(index_qref  ,ia) + qref(n)* wt
       xaop_ae%rAttr(index_ustar ,ia) = xaop_ae%rAttr(index_ustar ,ia) + ustar(n)*wt   ! friction velocity
       xaop_ae%rAttr(index_re    ,ia) = xaop_ae%rAttr(index_re    ,ia) + re(n)  * wt   ! reynolds number
       xaop_ae%rAttr(index_ssq   ,ia) = xaop_ae%rAttr(index_ssq   ,ia) + ssq(n) * wt   ! s.hum. saturation at Ts
       xaop_ae%rAttr(index_lwup  ,ia) = xaop_ae%rAttr(index_lwup  ,ia) + lwup(n)* wt   
       xaop_ae%rAttr(index_duu10n,ia) = xaop_ae%rAttr(index_duu10n,ia) + duu10n(n)*wt  
       xaop_ae%rAttr(index_u10   ,ia) = xaop_ae%rAttr(index_u10   ,ia) + sqrt(duu10n(n))*wt
       xaop_ae%rAttr(index_sumwt ,ia) = xaop_ae%rAttr(index_sumwt ,ia) + wt
    enddo

    call mct_aVect_clean(fractions_oe)

    !--- rearrange and sum from exchange grid to gsmap_a and gsmap_o decomps

    call mct_rearr_rearrange(xaop_ae, xaop_a, Re_e2a, sum=.true., &
         VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
    call mct_rearr_rearrange(xaop_oe, xaop_o, Re_e2o, sum=.true., &
         VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)

    !--- normalize by sum of wts associated with mapping

    do n = 1,nloc_a
       wt = xaop_a%rAttr(index_sumwt,n)
       if (wt /= 0.0_r8) then 
          wt = 1.0_r8/wt
       else
          wt = 1.0_r8
       endif
       xaop_a%rAttr(:,n) = xaop_a%rAttr(:,n) * wt
    enddo

    do n = 1,nloc_o
       wt = xaop_o%rAttr(index_sumwt,n)
       if (wt /= 0.0_r8) then 
          wt = 1.0_r8/wt
       else
          wt = 1.0_r8
       endif
       xaop_o%rAttr(:,n) = xaop_o%rAttr(:,n) * wt
    enddo

    !--- copy subset of fields to xao_a and xao_o and clean up

    call mct_avect_clean(xaop_ae)
    call mct_avect_clean(xaop_oe)

    call mct_avect_copy(xaop_a, xao_a)
    call mct_avect_copy(xaop_o, xao_o)

    call mct_avect_clean(xaop_a)
    call mct_avect_clean(xaop_o)

  end subroutine seq_flux_atmocnexch_mct

!===============================================================================

  subroutine seq_flux_atmocn_mct(infodata, tod, dt, a2x, o2x, xao)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_infodata_type) , intent(in)         :: infodata
    integer(in)             , intent(in)         :: tod,dt  ! NEW
    type(mct_aVect)         , intent(in)         :: a2x  ! a2x_ax or a2x_ox
    type(mct_aVect)         , intent(in)         :: o2x  ! o2x_ax or o2x_ox
    type(mct_aVect)         , intent(inout)      :: xao
    !
    ! Local variables
    !
    logical     :: flux_albav   ! flux avg option
    logical     :: dead_comps   ! .true.  => dead components are used
    integer(in) :: n,i          ! indices
    real(r8)    :: rlat         ! gridcell latitude in radians
    real(r8)    :: rlon         ! gridcell longitude in radians
    real(r8)    :: cosz         ! Cosine of solar zenith angle
    real(r8)    :: eccen        ! Earth orbit eccentricity
    real(r8)    :: mvelpp       ! Earth orbit
    real(r8)    :: lambm0       ! Earth orbit
    real(r8)    :: obliqr       ! Earth orbit
    real(r8)    :: delta        ! Solar declination angle  in radians
    real(r8)    :: eccf         ! Earth orbit eccentricity factor
    real(r8)    :: calday       ! calendar day including fraction, at 0e
    real(r8)    :: nextsw_cday  ! calendar day of next atm shortwave
    real(r8)    :: anidr        ! albedo: near infrared, direct
    real(r8)    :: avsdr        ! albedo: visible      , direct
    real(r8)    :: anidf        ! albedo: near infrared, diffuse
    real(r8)    :: avsdf        ! albedo: visible      , diffuse
    real(r8)    :: gust_fac = huge(1.0_r8) !wind gust factor
    integer(in) :: nloc, nloca, nloco    ! number of gridcells
    integer(in) :: ID           ! comm ID
    logical,save:: first_call = .true.
    logical     :: cold_start      ! .true. to initialize internal fields in shr_flux diurnal
    logical     :: read_restart    ! .true. => continue run
    logical     :: ocn_prognostic  ! .true. => ocn is prognostic 
    logical     :: flux_diurnal    ! .true. => turn on diurnal cycle in atm/ocn fluxes
    !
    real(r8),parameter :: albdif = 0.06_r8 ! 60 deg reference albedo, diffuse
    real(r8),parameter :: albdir = 0.07_r8 ! 60 deg reference albedo, direct 
    character(*),parameter :: subName =   '(seq_flux_atmocn_mct) '
    !
    !-----------------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         read_restart=read_restart, &
         flux_albav=flux_albav, &
         dead_comps=dead_comps, & 
         ocn_prognostic=ocn_prognostic, &
         flux_diurnal=flux_diurnal,     &
         gust_fac = gust_fac            )

    cold_start = .false.   ! use restart data or data from last timestep

    if (first_call) then
       if (.not.read_restart) cold_start = .true.
       index_xao_So_tref   = mct_aVect_indexRA(xao,'So_tref')
       index_xao_So_qref   = mct_aVect_indexRA(xao,'So_qref')
       index_xao_So_ustar  = mct_aVect_indexRA(xao,'So_ustar')  
       index_xao_So_re     = mct_aVect_indexRA(xao,'So_re')  
       index_xao_So_ssq    = mct_aVect_indexRA(xao,'So_ssq')
       index_xao_So_u10    = mct_aVect_indexRA(xao,'So_u10')
       index_xao_So_duu10n = mct_aVect_indexRA(xao,'So_duu10n')
       index_xao_Faox_taux = mct_aVect_indexRA(xao,'Faox_taux')
       index_xao_Faox_tauy = mct_aVect_indexRA(xao,'Faox_tauy')  
       index_xao_Faox_lat  = mct_aVect_indexRA(xao,'Faox_lat')   
       index_xao_Faox_sen  = mct_aVect_indexRA(xao,'Faox_sen')   
       index_xao_Faox_evap = mct_aVect_indexRA(xao,'Faox_evap')   
       index_xao_Faox_evap_16O = mct_aVect_indexRA(xao,'Faox_evap_16O', perrWith='quiet')
       index_xao_Faox_evap_HDO = mct_aVect_indexRA(xao,'Faox_evap_HDO', perrWith='quiet')
       index_xao_Faox_evap_18O = mct_aVect_indexRA(xao,'Faox_evap_18O', perrWith='quiet')
       index_xao_Faox_lwup = mct_aVect_indexRA(xao,'Faox_lwup')  
       index_xao_Faox_swdn = mct_aVect_indexRA(xao,'Faox_swdn')
       index_xao_Faox_swup = mct_aVect_indexRA(xao,'Faox_swup')
       index_xao_So_fswpen            = mct_aVect_indexRA(xao,'So_fswpen')
       index_xao_So_warm_diurn        = mct_aVect_indexRA(xao,'So_warm_diurn')
       index_xao_So_salt_diurn        = mct_aVect_indexRA(xao,'So_salt_diurn')
       index_xao_So_speed_diurn       = mct_aVect_indexRA(xao,'So_speed_diurn')
       index_xao_So_regime_diurn      = mct_aVect_indexRA(xao,'So_regime_diurn')
       index_xao_So_tskin_diurn       = mct_aVect_indexRA(xao,'So_tskin_diurn')
       index_xao_So_tskin_day_diurn   = mct_aVect_indexRA(xao,'So_tskin_day_diurn')
       index_xao_So_tskin_night_diurn = mct_aVect_indexRA(xao,'So_tskin_night_diurn')
       index_xao_So_cskin_diurn       = mct_aVect_indexRA(xao,'So_cskin_diurn')
       index_xao_So_cskin_night_diurn = mct_aVect_indexRA(xao,'So_cskin_night_diurn')
       index_xao_So_tbulk_diurn       = mct_aVect_indexRA(xao,'So_tbulk_diurn')
       index_xao_So_warmmax_diurn     = mct_aVect_indexRA(xao,'So_warmmax_diurn')
       index_xao_So_windmax_diurn     = mct_aVect_indexRA(xao,'So_windmax_diurn')
       index_xao_So_qsolavg_diurn     = mct_aVect_indexRA(xao,'So_qsolavg_diurn')
       index_xao_So_windavg_diurn     = mct_aVect_indexRA(xao,'So_windavg_diurn')
       index_xao_So_warmmaxinc_diurn  = mct_aVect_indexRA(xao,'So_warmmaxinc_diurn')
       index_xao_So_windmaxinc_diurn  = mct_aVect_indexRA(xao,'So_windmaxinc_diurn')
       index_xao_So_qsolinc_diurn     = mct_aVect_indexRA(xao,'So_qsolinc_diurn')
       index_xao_So_windinc_diurn     = mct_aVect_indexRA(xao,'So_windinc_diurn')
       index_xao_So_ninc_diurn        = mct_aVect_indexRA(xao,'So_ninc_diurn')
       
       index_a2x_Sa_z      = mct_aVect_indexRA(a2x,'Sa_z')
       index_a2x_Sa_u      = mct_aVect_indexRA(a2x,'Sa_u')
       index_a2x_Sa_v      = mct_aVect_indexRA(a2x,'Sa_v')
       index_a2x_Sa_tbot   = mct_aVect_indexRA(a2x,'Sa_tbot')
       index_a2x_Sa_ptem   = mct_aVect_indexRA(a2x,'Sa_ptem')
       index_a2x_Sa_shum   = mct_aVect_indexRA(a2x,'Sa_shum')
       index_a2x_Sa_shum_16O   = mct_aVect_indexRA(a2x,'Sa_shum_16O', perrWith='quiet')
       index_a2x_Sa_shum_HDO   = mct_aVect_indexRA(a2x,'Sa_shum_HDO', perrWith='quiet')
       index_a2x_Sa_shum_18O   = mct_aVect_indexRA(a2x,'Sa_shum_18O', perrWith='quiet')
       index_a2x_Sa_dens   = mct_aVect_indexRA(a2x,'Sa_dens')
       index_a2x_Faxa_lwdn = mct_aVect_indexRA(a2x,'Faxa_lwdn')
       index_a2x_Faxa_rainc= mct_aVect_indexRA(a2x,'Faxa_rainc')
       index_a2x_Faxa_rainl= mct_aVect_indexRA(a2x,'Faxa_rainl')
       index_a2x_Faxa_snowc= mct_aVect_indexRA(a2x,'Faxa_snowc')
       index_a2x_Faxa_snowl= mct_aVect_indexRA(a2x,'Faxa_snowl')
       
       index_o2x_So_t      = mct_aVect_indexRA(o2x,'So_t')
       index_o2x_So_u      = mct_aVect_indexRA(o2x,'So_u')
       index_o2x_So_v      = mct_aVect_indexRA(o2x,'So_v')
       index_o2x_So_fswpen = mct_aVect_indexRA(o2x,'So_fswpen')
       index_o2x_So_s      = mct_aVect_indexRA(o2x,'So_s')
       index_o2x_So_roce_16O = mct_aVect_indexRA(o2x,'So_roce_16O', perrWith='quiet')
       index_o2x_So_roce_HDO = mct_aVect_indexRA(o2x,'So_roce_HDO', perrWith='quiet')
       index_o2x_So_roce_18O = mct_aVect_indexRA(o2x,'So_roce_18O', perrWith='quiet')
       first_call = .false.
    end if
       
    if (trim(fluxsetting) /= trim(fluxsetting_atmocn)) then
       call shr_sys_abort(trim(subname)//' ERROR wrong fluxsetting')
    endif

    nloc = mct_aVect_lsize(xao)
    nloca = mct_aVect_lsize(a2x)
    nloco = mct_aVect_lsize(o2x)

    if (nloc /= nloca .or. nloc /= nloco) then
       call shr_sys_abort(trim(subname)//' ERROR nloc sizes do not match')
    endif

    ! Update ocean surface fluxes 
    ! Must fabricate "reasonable" data (when using dead components)

    emask = mask
    if (dead_comps) then
       do n = 1,nloc
          mask(n) =   1      ! ocn domain mask            ~ 0 <=> inactive cell
          tocn(n) = 290.0_r8 ! ocn temperature            ~ Kelvin
          uocn(n) =   0.0_r8 ! ocn velocity, zonal        ~ m/s
          vocn(n) =   0.0_r8 ! ocn velocity, meridional   ~ m/s
          zbot(n) =  55.0_r8 ! atm height of bottom layer ~ m
          ubot(n) =   0.0_r8 ! atm velocity, zonal        ~ m/s
          vbot(n) =   2.0_r8 ! atm velocity, meridional   ~ m/s
          thbot(n)= 301.0_r8 ! atm potential temperature  ~ Kelvin
          shum(n) = 1.e-2_r8 ! atm specific humidity      ~ kg/kg
!wiso note: shum_* should be multiplied by Rstd_* here?
          shum_16O(n) = 1.e-2_r8 ! H216O specific humidity ~ kg/kg
          shum_HDO(n) = 1.e-2_r8 ! HD16O specific humidity ~ kg/kg
          shum_18O(n) = 1.e-2_r8 ! H218O specific humidity ~ kg/kg
          roce_16O(n) = 1.0_r8   ! H216O surface ratio     ~ mol/mol
          roce_HDO(n) = 1.0_r8   ! HDO   surface ratio     ~ mol/mol 
          roce_18O(n) = 1.0_r8   ! H218O surface ratio     ~ mol/mol 
          dens(n) =   1.0_r8 ! atm density                ~ kg/m^3
          tbot(n) = 300.0_r8 ! atm temperature            ~ Kelvin
          uGust(n)=   0.0_r8
          lwdn(n) =   0.0_r8
          prec(n) =   0.0_r8
          prec_gust(n) =  0.0_r8
          fswpen(n)=  0.0_r8
          ocnsal(n)=  0.0_r8

          warm       (n) = 0.0_r8
          salt       (n) = 0.0_r8
          speed      (n) = 0.0_r8
          regime     (n) = 0.0_r8
          warmMax    (n) = 0.0_r8
          windMax    (n) = 0.0_r8
          qSolAvg    (n) = 0.0_r8
          windAvg    (n) = 0.0_r8
          warmMaxInc (n) = 0.0_r8
          windMaxInc (n) = 0.0_r8
          qSolInc    (n) = 0.0_r8
          windInc    (n) = 0.0_r8
          nInc       (n) = 0.0_r8
          tbulk      (n) = 0.0_r8
          tskin      (n) = 0.0_r8
          tskin_day  (n) = 0.0_r8
          tskin_night(n) = 0.0_r8
          cskin      (n) = 0.0_r8
          cskin_night(n) = 0.0_r8
          swdn       (n) = 0.0_r8
          swup       (n) = 0.0_r8
       enddo
    else
       do n = 1,nloc
          nInc(n) = 0._r8 ! needed for minval/maxval calculation  
          if (mask(n) /= 0) then	
             zbot(n) = a2x%rAttr(index_a2x_Sa_z   ,n)
             ubot(n) = a2x%rAttr(index_a2x_Sa_u   ,n)
             vbot(n) = a2x%rAttr(index_a2x_Sa_v   ,n)
             thbot(n)= a2x%rAttr(index_a2x_Sa_ptem,n)
             shum(n) = a2x%rAttr(index_a2x_Sa_shum,n)
             if ( index_a2x_Sa_shum_16O /= 0 ) shum_16O(n) = a2x%rAttr(index_a2x_Sa_shum_16O,n)
             if ( index_a2x_Sa_shum_HDO /= 0 ) shum_HDO(n) = a2x%rAttr(index_a2x_Sa_shum_HDO,n)
             if ( index_a2x_Sa_shum_18O /= 0 ) shum_18O(n) = a2x%rAttr(index_a2x_Sa_shum_18O,n)
             dens(n) = a2x%rAttr(index_a2x_Sa_dens,n)
             tbot(n) = a2x%rAttr(index_a2x_Sa_tbot,n)
             tocn(n) = o2x%rAttr(index_o2x_So_t   ,n)   
             uocn(n) = o2x%rAttr(index_o2x_So_u   ,n)
             vocn(n) = o2x%rAttr(index_o2x_So_v   ,n)
             if ( index_o2x_So_roce_16O /= 0 ) roce_16O(n) = o2x%rAttr(index_o2x_So_roce_16O, n)
             if ( index_o2x_So_roce_HDO /= 0 ) roce_HDO(n) = o2x%rAttr(index_o2x_So_roce_HDO, n)
             if ( index_o2x_So_roce_18O /= 0 ) roce_18O(n) = o2x%rAttr(index_o2x_So_roce_18O, n)
             !--- mask missing atm or ocn data if found
             if (dens(n) < 1.0e-12 .or. tocn(n) < 1.0) then
                emask(n) = 0
                !write(logunit,*) 'aoflux tcx1',n,dens(n),tocn(n)
             endif
!           !!uGust(n) = 1.5_r8*sqrt(uocn(n)**2 + vocn(n)**2) ! there is no wind gust data from ocn
             uGust(n) = 0.0_r8
             lwdn (n) = a2x%rAttr(index_a2x_Faxa_lwdn ,n)
             prec (n) = a2x%rAttr(index_a2x_Faxa_rainc,n) &
                    & + a2x%rAttr(index_a2x_Faxa_rainl,n) &
                    & + a2x%rAttr(index_a2x_Faxa_snowc,n) &
                    & + a2x%rAttr(index_a2x_Faxa_snowl,n)
             prec_gust (n) = a2x%rAttr(index_a2x_Faxa_rainc,n)
             fswpen(n)= o2x%rAttr(index_o2x_So_fswpen ,n)
             ocnsal(n)= o2x%rAttr(index_o2x_So_s      ,n)

             warm       (n) = xao%rAttr(index_xao_So_warm_diurn      ,n)
             salt       (n) = xao%rAttr(index_xao_So_salt_diurn      ,n)
             speed      (n) = xao%rAttr(index_xao_So_speed_diurn     ,n)
             regime     (n) = xao%rAttr(index_xao_So_regime_diurn    ,n)
             warmMax    (n) = xao%rAttr(index_xao_So_warmMax_diurn   ,n)
             windMax    (n) = xao%rAttr(index_xao_So_windMax_diurn   ,n)
             qSolAvg    (n) = xao%rAttr(index_xao_So_qsolavg_diurn   ,n)
             windAvg    (n) = xao%rAttr(index_xao_So_windavg_diurn   ,n)
             warmMaxInc (n) = xao%rAttr(index_xao_So_warmMaxInc_diurn,n)
             windMaxInc (n) = xao%rAttr(index_xao_So_windMaxInc_diurn,n)
             qSolInc    (n) = xao%rAttr(index_xao_So_qSolInc_diurn   ,n)
             windInc    (n) = xao%rAttr(index_xao_So_windInc_diurn   ,n)
             nInc       (n) = xao%rAttr(index_xao_So_nInc_diurn      ,n)
             tbulk      (n) = xao%rAttr(index_xao_So_tbulk_diurn     ,n)
             tskin      (n) = xao%rAttr(index_xao_So_tskin_diurn     ,n)
             tskin_day  (n) = xao%rAttr(index_xao_So_tskin_day_diurn ,n)
             tskin_night(n) = xao%rAttr(index_xao_So_tskin_night_diurn,n)
             cskin      (n) = xao%rAttr(index_xao_So_cskin_diurn     ,n)
             cskin_night(n) = xao%rAttr(index_xao_So_cskin_night_diurn,n)
             ! set in flux_ocnalb using data from previous timestep
             swdn       (n) = xao%rAttr(index_xao_Faox_swdn          ,n)
             swup       (n) = xao%rAttr(index_xao_Faox_swup          ,n)
          end if
       enddo
    end if

    if (flux_diurnal) then
       call shr_flux_atmocn_diurnal (nloc , zbot , ubot, vbot, thbot, &
                          shum , shum_16O , shum_HDO, shum_18O, dens , tbot, uocn, vocn , &
                          tocn , emask, sen , lat , lwup , &
                          roce_16O, roce_HDO, roce_18O,    &
                          evap , evap_16O, evap_HDO, evap_18O, taux , tauy, tref, qref , &
                          uGust, lwdn , swdn , swup, prec, &
                          fswpen, ocnsal, ocn_prognostic, flux_diurnal,    &
                          lats, lons , warm , salt , speed, regime,       &
                          warmMax, windMax, qSolAvg, windAvg,             &
                          warmMaxInc, windMaxInc, qSolInc, windInc, nInc, &
                          tbulk, tskin, tskin_day, tskin_night, &
                          cskin, cskin_night, tod, dt,          &
                          duu10n,ustar, re  , ssq, &
                          !missval should not be needed if flux calc 
                          !consistent with mrgx2a fraction
                          !duu10n,ustar, re  , ssq, missval = 0.0_r8 )
                          cold_start=cold_start)
    else
       call shr_flux_atmocn (nloc , zbot , ubot, vbot, thbot, prec_gust, gust_fac, &
                          shum , shum_16O , shum_HDO, shum_18O, dens , tbot, uocn, vocn , &
                          tocn , emask, sen , lat , lwup , &
                          roce_16O, roce_HDO, roce_18O,    &
                          evap , evap_16O, evap_HDO, evap_18O, taux , tauy, tref, qref , &
                          duu10n,ustar, re  , ssq)
                          !missval should not be needed if flux calc 
                          !consistent with mrgx2a fraction
                          !duu10n,ustar, re  , ssq, missval = 0.0_r8 )
    endif

    do n = 1,nloc
       if (mask(n) /= 0) then	
          xao%rAttr(index_xao_Faox_sen ,n) = sen(n)
          xao%rAttr(index_xao_Faox_lat ,n) = lat(n)
          xao%rAttr(index_xao_Faox_taux,n) = taux(n)
          xao%rAttr(index_xao_Faox_tauy,n) = tauy(n)
          xao%rAttr(index_xao_Faox_evap,n) = evap(n)
          if ( index_xao_Faox_evap_16O /= 0 ) xao%rAttr(index_xao_Faox_evap_16O,n) = evap_16O(n)
          if ( index_xao_Faox_evap_HDO /= 0 ) xao%rAttr(index_xao_Faox_evap_HDO,n) = evap_HDO(n)
          if ( index_xao_Faox_evap_18O /= 0 ) xao%rAttr(index_xao_Faox_evap_18O,n) = evap_18O(n)
          xao%rAttr(index_xao_So_tref  ,n) = tref(n)
	  xao%rAttr(index_xao_So_qref  ,n) = qref(n)
          xao%rAttr(index_xao_So_ustar ,n) = ustar(n)  ! friction velocity
          xao%rAttr(index_xao_So_re    ,n) = re(n)     ! reynolds number
          xao%rAttr(index_xao_So_ssq   ,n) = ssq(n)    ! s.hum. saturation at Ts
          xao%rAttr(index_xao_Faox_lwup,n) = lwup(n)   
          xao%rAttr(index_xao_So_duu10n,n) = duu10n(n)  
          xao%rAttr(index_xao_So_u10   ,n) = sqrt(duu10n(n))  
          xao%rAttr(index_xao_So_warm_diurn       ,n) = warm(n)
          xao%rAttr(index_xao_So_salt_diurn       ,n) = salt(n)
          xao%rAttr(index_xao_So_speed_diurn      ,n) = speed(n)
          xao%rAttr(index_xao_So_regime_diurn     ,n) = regime(n)
          xao%rAttr(index_xao_So_warmMax_diurn    ,n) = warmMax(n)
          xao%rAttr(index_xao_So_windMax_diurn    ,n) = windMax(n)
          xao%rAttr(index_xao_So_qSolAvg_diurn    ,n) = qSolAvg(n)
          xao%rAttr(index_xao_So_windAvg_diurn    ,n) = windAvg(n)
          xao%rAttr(index_xao_So_warmMaxInc_diurn ,n) = warmMaxInc(n)
          xao%rAttr(index_xao_So_windMaxInc_diurn ,n) = windMaxInc(n)
          xao%rAttr(index_xao_So_qSolInc_diurn    ,n) = qSolInc(n)
          xao%rAttr(index_xao_So_windInc_diurn    ,n) = windInc(n)
          xao%rAttr(index_xao_So_nInc_diurn       ,n) = nInc(n)
          xao%rAttr(index_xao_So_tbulk_diurn      ,n) = tbulk(n)
          xao%rAttr(index_xao_So_tskin_diurn      ,n) = tskin(n)
          xao%rAttr(index_xao_So_tskin_day_diurn  ,n) = tskin_day(n)
          xao%rAttr(index_xao_So_tskin_night_diurn,n) = tskin_night(n)
          xao%rAttr(index_xao_So_cskin_diurn      ,n) = cskin(n)
          xao%rAttr(index_xao_So_cskin_night_diurn,n) = cskin_night(n)
          xao%rAttr(index_xao_So_fswpen           ,n) = fswpen(n)
       end if
    enddo

  end subroutine seq_flux_atmocn_mct

!===============================================================================

end module seq_flux_mct
