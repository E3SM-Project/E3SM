module TopounitDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Topounit data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use clm_varcon     , only : spval, ispval
  use clm_varctl     , only : iulog, use_cn, use_fates, use_lch4
  use clm_varpar     , only : numrad
  use decompMod      , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

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
  ! declare the public instances of topounit data types
  type(topounit_atmospheric_state),    public, target :: top_as
  type(topounit_atmospheric_flux),     public, target :: top_af
  type(topounit_energy_state),         public, target :: top_es

  contains

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
  

end module TopounitDataType
