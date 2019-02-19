module TopounitType
  
  ! -------------------------------------------------------- 
  ! ELM sub-grid hierarchy:
  ! Define topographic unit data types, with Init and Clean for each
  ! Init() calls allocate memory and set history fields
  ! -------------------------------------------------------- 
  ! 3 Aug 2015, PET
  
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use landunit_varcon, only : max_lunit
  use clm_varcon     , only : ispval, spval
  use clm_varpar     , only : numrad
  use clm_varctl     , only : iulog, use_cn, use_fates
  use decompMod      , only : bounds_type

  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! sub-grid topology and physical properties defined at the topographic unit level
  type, public :: topounit_physical_properties

    ! indices and weights for higher subgrid level (gridcell)
    integer , pointer :: gridcell   (:) => null() ! index into gridcell level quantities
    real(r8), pointer :: wtgcell    (:) => null() ! weight (relative to gridcell)

    ! Starting and ending indices for all subgrid types below the landunit level
    integer , pointer :: lndi       (:) => null() ! beginning landunit index for each topounit
    integer , pointer :: lndf       (:) => null() ! ending landunit index for each topounit
    integer , pointer :: nlandunits (:) => null() ! number of landunits for each topounit
    integer , pointer :: coli       (:) => null() ! beginning column index per landunit
    integer , pointer :: colf       (:) => null() ! ending column index for each landunit
    integer , pointer :: ncolumns   (:) => null() ! number of columns for each landunit
    integer , pointer :: pfti       (:) => null() ! beginning pft index for each landunit
    integer , pointer :: pftf       (:) => null() ! ending pft index for each landunit
    integer , pointer :: npfts      (:) => null() ! number of patches for each landunit

    ! indices into landunit-level arrays for landunits in this topounit (ispval implies
    ! this landunit doesn't exist on this topounit) [1:max_lunit, begt:endt]
    ! (note that the spatial dimension is last here, in contrast to most 2-d variables;
    ! this is for efficiency, since most loops will go over t in the outer loop, and
    ! landunit type in the inner loop)
    integer , pointer :: landunit_indices (:,:) => null() 
    
    ! physical properties
    real(r8), pointer :: area       (:) => null() ! land area (km^2)
    real(r8), pointer :: lat        (:) => null() ! mean latitude (radians)
    real(r8), pointer :: lon        (:) => null() ! mean longitude (radians)
    real(r8), pointer :: elevation  (:) => null() ! mean soil surface elevation, above mean sea level (m)
    real(r8), pointer :: slope      (:) => null() ! mean slope angle (radians)
    real(r8), pointer :: aspect     (:) => null() ! mean aspect angle, measured clockwise from north (radians)
    real(r8), pointer :: emissivity (:) => null() ! mean surface emissivity
    real(r8), pointer :: surfalb_dir(:,:) => null() ! (topunit,numrad) mean surface albedo (direct)
    real(r8), pointer :: surfalb_dif(:,:) => null() ! (topunit,numrad) mean surface albedo (diffuse)
  contains
    procedure, public :: Init  => init_top_pp
    procedure, public :: Clean => clean_top_pp  
  end type topounit_physical_properties
    
  !-----------------------------------------------------------------------
  ! Define the data structure where land model receives atmospheric state information.
  type, public :: topounit_atmospheric_state
    real(r8), pointer :: tbot       (:) => null() ! temperature of air at atmospheric forcing height (K)
    real(r8), pointer :: thbot      (:) => null() ! potential temperature of air at atmospheric forcing height (K)
    real(r8), pointer :: pbot       (:) => null() ! air pressure at atmospheric forcing height (Pa)
    real(r8), pointer :: rhobot     (:) => null() ! air density at atmospheric forcing height (kg/m**3)
    real(r8), pointer :: qbot       (:) => null() ! specific humidity at atmospheric forcing height (kg H2O/kg moist air)
    real(r8), pointer :: rhbot      (:) => null() ! relative humidity at atmospheric forcing height (%)
    real(r8), pointer :: ubot       (:) => null() ! wind speed in U (east) direction at atmospheric forcing height (m/s)
    real(r8), pointer :: vbot       (:) => null() ! wind speed in V (north) direction at atmospheric forcing height (m/s)
    real(r8), pointer :: windbot    (:) => null() ! horizontal component of wind at atmospheric forcing height (m/s)
    real(r8), pointer :: zbot       (:) => null() ! atmospheric forcing height (m)
    real(r8), pointer :: po2bot     (:) => null() ! partial pressure of O2 at atmospheric forcing height (Pa)
    real(r8), pointer :: pco2bot    (:) => null() ! partial pressure of CO2 at atmospheric forcing height (Pa) 
    real(r8), pointer :: pc13o2bot  (:) => null() ! partial pressure of C13O2 at atmospheric forcing height (Pa) 
    real(r8), pointer :: pch4bot    (:) => null() ! partial pressure of CH4 at atmospheric forcing height (Pa)
    ! Accumulated fields
    real(r8), pointer :: rh24h      (:) => null() ! 24-hour running mean of relative humidity at atmospheric forcing height (%)
    real(r8), pointer :: wind24h    (:) => null() ! 24-hour running mean of horizontal wind at atmospheric forcing height (m/s)
  contains
    procedure, public :: Init  => init_top_as
    procedure, public :: Clean => clean_top_as
    !procedure, public :: InitAccBuffer => init_acc_buffer_top_as
    !procedure, public :: InitAccVars   => init_acc_vars_top_as
    !procedure, public :: UpdateAccVars => update_acc_vars_top_as
  end type topounit_atmospheric_state

  !-----------------------------------------------------------------------
  ! Define the data structure that where land model receives atmospheric flux information.
  type, public :: topounit_atmospheric_flux
    real(r8), pointer :: rain      (:)   => null() ! rain rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: snow      (:)   => null() ! snow rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: solad     (:,:) => null() ! direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll) (W/m**2)
    real(r8), pointer :: solai     (:,:) => null() ! diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld) (W/m**2)
    real(r8), pointer :: solar     (:)   => null() ! incident solar radiation (W/m**2)
    real(r8), pointer :: lwrad     (:)   => null() ! atm downwrd IR longwave radiation (W/m**2) 
    ! Accumulated fields           
    real(r8), pointer :: prec24h   (:)   => null() ! 24-hour mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: prec10d   (:)   => null() ! 10-day mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: prec60d   (:)   => null() ! 60-day mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: fsd24h    (:)   => null() ! 24hr average of direct beam radiation (W/m**2)
    real(r8), pointer :: fsd240h   (:)   => null() ! 240hr average of direct beam radiation (W/m**2) 
    real(r8), pointer :: fsi24h    (:)   => null() ! 24hr average of diffuse beam radiation (W/m**2) 
    real(r8), pointer :: fsi240h   (:)   => null() ! 240hr average of diffuse beam radiation (W/m**2) 
    
  contains
    procedure, public :: Init  => init_top_af
    procedure, public :: Clean => clean_top_af
    !procedure, public :: InitAccBuffer => init_acc_buffer_top_af
    !procedure, public :: InitAccVars   => init_acc_vars_top_af
    !procedure, public :: UpdateAccVars => update_acc_vars_top_af
  end type topounit_atmospheric_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information for land at the level of topographic unit.
  type, public :: topounit_energy_state
    real(r8), pointer :: t_rad      (:) => null() ! mean radiative temperature of land surface (K)
  contains
    procedure, public :: Init  => init_top_es
    procedure, public :: Clean => clean_top_es
  end type topounit_energy_state

  !-----------------------------------------------------------------------
  ! declare the public instances of topounit types
  type(topounit_physical_properties),  public, target :: top_pp
  type(topounit_atmospheric_state),    public, target :: top_as
  type(topounit_atmospheric_flux),     public, target :: top_af
  type(topounit_energy_state),         public, target :: top_es
  
  contains
  
  !-----------------------------------------------------------------------
  subroutine init_top_pp(this, begt, endt)
    class(topounit_physical_properties) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index
    
    allocate(this%gridcell  (begt:endt)) ; this%gridcell  (:) = ispval
    allocate(this%wtgcell   (begt:endt)) ; this%wtgcell   (:) = nan
    allocate(this%lndi      (begt:endt)) ; this%lndi      (:) = ispval
    allocate(this%lndf      (begt:endt)) ; this%lndf      (:) = ispval
    allocate(this%nlandunits(begt:endt)) ; this%nlandunits(:) = ispval
    allocate(this%coli      (begt:endt)) ; this%coli      (:) = ispval
    allocate(this%colf      (begt:endt)) ; this%colf      (:) = ispval
    allocate(this%ncolumns  (begt:endt)) ; this%ncolumns  (:) = ispval
    allocate(this%pfti      (begt:endt)) ; this%pfti      (:) = ispval
    allocate(this%pftf      (begt:endt)) ; this%pftf      (:) = ispval
    allocate(this%npfts     (begt:endt)) ; this%npfts     (:) = ispval
    
    allocate(this%landunit_indices(1:max_lunit, begt:endt)); this%landunit_indices(:,:) = ispval

    allocate(this%area        (begt:endt)) ; this%area        (:) = nan
    allocate(this%lat         (begt:endt)) ; this%lat         (:) = nan
    allocate(this%lon         (begt:endt)) ; this%lon         (:) = nan
    allocate(this%elevation   (begt:endt)) ; this%elevation   (:) = nan
    allocate(this%slope       (begt:endt)) ; this%slope       (:) = nan
    allocate(this%aspect      (begt:endt)) ; this%aspect      (:) = nan
    allocate(this%emissivity  (begt:endt)) ; this%emissivity  (:) = nan
    allocate(this%surfalb_dir (begt:endt,1:numrad)) ; this%surfalb_dir(:,:) = nan
    allocate(this%surfalb_dif (begt:endt,1:numrad)) ; this%surfalb_dif(:,:) = nan
  end subroutine init_top_pp
  
  !-----------------------------------------------------------------------
  subroutine clean_top_pp(this)
    class(topounit_physical_properties) :: this
  
    deallocate(this%gridcell    )
    deallocate(this%wtgcell     )
    deallocate(this%lndi        )
    deallocate(this%lndf        )
    deallocate(this%nlandunits  )
    deallocate(this%coli        )
    deallocate(this%colf        )
    deallocate(this%ncolumns    )
    deallocate(this%pfti        )
    deallocate(this%pftf        )
    deallocate(this%npfts       )
    deallocate(this%landunit_indices )

    deallocate(this%area        )
    deallocate(this%lat         )
    deallocate(this%lon         )
    deallocate(this%elevation   )
    deallocate(this%slope       )
    deallocate(this%aspect      )
    deallocate(this%emissivity  )
    deallocate(this%surfalb_dir )
    deallocate(this%surfalb_dif )
  end subroutine clean_top_pp
  
  !-----------------------------------------------------------------------
  subroutine init_top_as(this, begt, endt)
    class(topounit_atmospheric_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index
   
    ! Allocate for atmospheric state forcing variables, initialize to special value
    allocate(this%tbot     (begt:endt)) ; this%tbot      (:) = spval
    allocate(this%thbot    (begt:endt)) ; this%thbot     (:) = spval
    allocate(this%pbot     (begt:endt)) ; this%pbot      (:) = spval
    allocate(this%rhobot   (begt:endt)) ; this%rhobot    (:) = spval
    allocate(this%qbot     (begt:endt)) ; this%qbot      (:) = spval
    allocate(this%rhbot    (begt:endt)) ; this%rhbot     (:) = spval
    allocate(this%ubot     (begt:endt)) ; this%ubot      (:) = spval
    allocate(this%vbot     (begt:endt)) ; this%vbot      (:) = spval
    allocate(this%windbot  (begt:endt)) ; this%windbot   (:) = spval
    allocate(this%zbot     (begt:endt)) ; this%zbot      (:) = spval
    allocate(this%po2bot   (begt:endt)) ; this%po2bot    (:) = spval
    allocate(this%pco2bot  (begt:endt)) ; this%pco2bot   (:) = spval
    allocate(this%pc13o2bot(begt:endt)) ; this%pc13o2bot (:) = spval
    allocate(this%pch4bot  (begt:endt)) ; this%pch4bot   (:) = spval
    if (use_fates) then
      allocate(this%rh24h  (begt:endt)) ; this%rh24h     (:) = spval
      allocate(this%wind24h(begt:endt)) ; this%wind24h   (:) = spval
    end if
    
    ! Set history fields for atmospheric state forcing variables
    !call hist_addfld1d (fname='TBOT', units='K',  &
    !     avgflag='A', long_name='temperature of air at atmospheric forcing height', &
    !     ptr_lnd=this%tbot)

    !call hist_addfld1d (fname='THBOT', units='K',  &
    !     avgflag='A', long_name='potential temperature of air at atmospheric forcing height', &
    !     ptr_lnd=this%thbot)
  end subroutine init_top_as

  !-----------------------------------------------------------------------
  subroutine clean_top_as(this, begt, endt)
    class(topounit_atmospheric_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    deallocate(this%tbot)
    deallocate(this%thbot)
    deallocate(this%pbot)
    deallocate(this%rhobot)
    deallocate(this%qbot)
    deallocate(this%rhbot)
    deallocate(this%ubot)
    deallocate(this%vbot)
    deallocate(this%windbot)
    deallocate(this%zbot)
    deallocate(this%po2bot)
    deallocate(this%pco2bot)
    deallocate(this%pc13o2bot)
    deallocate(this%pch4bot)
    if (use_fates) then
      deallocate(this%rh24h)
      deallocate(this%wind24h)
    end if
  end subroutine clean_top_as

  !-----------------------------------------------------------------------
  subroutine init_top_af(this, begt, endt)
    class(topounit_atmospheric_flux) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index
   
    ! Allocate for atmospheric flux forcing variables, initialize to special value
    allocate(this%rain     (begt:endt))          ; this%rain      (:) = spval
    allocate(this%snow     (begt:endt))          ; this%snow      (:) = spval
    allocate(this%solad    (begt:endt, numrad))  ; this%solad     (:,:) = spval
    allocate(this%solai    (begt:endt, numrad))  ; this%solai     (:,:) = spval
    allocate(this%solar    (begt:endt))          ; this%solar     (:) = spval
    allocate(this%lwrad    (begt:endt))          ; this%lwrad     (:) = spval
    if (use_fates) then
      allocate(this%prec24h  (begt:endt)) ; this%prec24h   (:) = spval
    end if
    if (use_cn) then
      allocate(this%prec10d  (begt:endt)) ; this%prec10d   (:) = spval
      allocate(this%prec60d  (begt:endt)) ; this%prec60d   (:) = spval
    end if
    allocate(this%fsd24h   (begt:endt))          ; this%fsd24h    (:) = spval
    allocate(this%fsd240h  (begt:endt))          ; this%fsd240h   (:) = spval
    allocate(this%fsi24h   (begt:endt))          ; this%fsi24h    (:) = spval
    allocate(this%fsi240h  (begt:endt))          ; this%fsi240h   (:) = spval
    
    ! Set history fields for atmospheric flux forcing variables
    !call hist_addfld1d (fname='TBOT', units='K',  &
    !     avgflag='A', long_name='temperature of air at atmospheric forcing height', &
    !     ptr_lnd=this%tbot)

    !call hist_addfld1d (fname='THBOT', units='K',  &
    !     avgflag='A', long_name='potential temperature of air at atmospheric forcing height', &
    !     ptr_lnd=this%thbot)
  end subroutine init_top_af

  !-----------------------------------------------------------------------
  subroutine clean_top_af(this, begt, endt)
    class(topounit_atmospheric_flux) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    deallocate(this%rain)
    deallocate(this%snow)
    deallocate(this%solad)
    deallocate(this%solai)
    deallocate(this%solar)
    deallocate(this%lwrad)
    if (use_fates) then
      deallocate(this%prec24h)
    end if
    if (use_cn) then
      deallocate(this%prec10d)
      deallocate(this%prec60d)
    end if
    deallocate(this%fsd24h)
    deallocate(this%fsd240h)
    deallocate(this%fsi24h)
    deallocate(this%fsi240h)
  end subroutine clean_top_af

  !-----------------------------------------------------------------------
  subroutine init_top_es(this, begt, endt)
    class(topounit_energy_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    allocate(this%t_rad   (begt:endt)) ; this%t_rad   (:) = nan
  end subroutine init_top_es
  
  !-----------------------------------------------------------------------
  subroutine clean_top_es(this, begt, endt)
    class(topounit_energy_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index
    
    deallocate(this%t_rad    )
  end subroutine clean_top_es
  
    
end module TopounitType
