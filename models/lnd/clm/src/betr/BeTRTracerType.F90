module BeTRTracerType
  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module variables for calculating the tracer transport and bgc reactions.
  !
  ! !USES:
  use shr_kind_mod       , only: r8 => shr_kind_r8
  use decompMod          , only : bounds_type
  !
  implicit none
  private

  !----------------------------------------------------
  !betr tracer setup structure
  !----------------------------------------------------
  type, public :: BeTRtracer_type
   character(len=255) :: betr_simname                                   ! name of the simulation
   integer :: ntracers                                                  ! total number of tracers, gas/aqueous tracers + solid tracers that undergo active mineral protection
   integer :: ngwmobile_tracers                                         ! total number of tracers potentially undergoing gas/aqueous movement
   integer :: nvolatile_tracers                                         ! number of volatile_tracers
   integer :: nsolid_equil_tracers                                      ! number of tracers that undergo equilibrium adsorption in soil could include adsorbed doc, nh4(+)
   integer :: nsolid_passive_tracers                                    ! number of tracers that undergo active mineral protection   
   integer :: nco2_tags                                                 ! number of tagged co2 tracers
   integer :: nh2o_tracers                                              ! number of h2o tracers, this will be used to compute vapor gradient and thermal gradient driven isotopic flow
   logical :: is_oddstep = .true.                                       !this is not used now, originally was included to set up alternative numerical methods
   integer :: id_trc_n2         ! tag for n2
   integer :: id_trc_o2         ! tag for co2
   integer :: id_trc_ar         ! tag for ar
   integer :: id_trc_co2x       ! tag for co2 and its related species, co2x(CO2, H2CO3, HCO3(-), CO3(2-)),
   integer :: id_trc_ch4        ! tag for methane
   
   integer :: id_trc_no         ! tag for no
   integer :: id_trc_n2o        ! tag for n2o
   integer :: id_trc_air_co2x   ! tag for atmospheric co2
   integer :: id_trc_arrt_co2x  ! tag for autotrophic co2
   integer :: id_trc_hrsoi_co2x ! tag for heterotrophic co2   
   integer :: id_trc_nh3x       ! tag for nh3 and its related species, nh3x(NH3, NH4OH,NH4(+))
   integer :: id_trc_no3x       ! tag for no3 and its related species, no3x(HNO3,NO3(-))
   integer :: id_trc_no2x       ! tag for no2 and its related species, no2x(HNO2,NO2(-))
   integer :: id_trc_dom        ! tag for generic dissolved organic matter
   integer :: id_trc_doc        ! tag for generic dissolved organic carbon, used for testing single carbon pool model
   
   integer :: id_trc_o18_h2o    ! tag for H2O(18)
   integer :: id_trc_o17_h2o    ! tag for H2O(17)
   integer :: id_trc_o18_h2o_ice! tag for H2O(18) in ice
   integer :: id_trc_d_h2o      ! tag for DHO
   integer :: id_trc_d_h2o_ice  ! tag for DHO in ice
   integer :: id_trc_c13_co2x   ! tag for C(13)O2 and its related species
   integer :: id_trc_c14_co2x   ! tag for C(14)O2 and its related species
   integer :: id_trc_o18_co2x   ! tag for O(18)CO and its related species
   integer :: id_trc_o17_co2x   ! tag for O(17)CO and its related species
   
   integer :: id_trc_o18_o2     ! tag for O(18)O and its related species
   integer :: id_trc_o17_o2     ! tag for O(17)O and its related species  

   logical, pointer :: is_volatile(:)                                   !flag for volatile species,  true/false, (yes/no)
   logical, pointer :: is_adsorb(:)                                     !flag for adsorbable species, true/false (year/no), in equilibrium with aqueous phase and/or gaseous phase
   logical, pointer :: is_advective(:)                                  !flag for advective species, some species, like non-dissolved som does not undergo advection, rather bioturbation is the major mechanism for vertical transport
   logical, pointer :: is_mobile(:)                                     !flag indicating whether the tracer is mobile or inert, when it is innert, do not move it around
   logical, pointer :: is_h2o(:)                                        !flag for water isotope
   logical, pointer :: is_co2tag(:)                                     !tagged co2 tracer?
   logical, pointer :: is_dom(:)                                        !true if it is a dom tracer, place holder for rtm bgc
   logical, pointer :: is_isotope(:)
   integer, pointer :: refisoid(:)                                      !reference tracer for isotope calculation, this is setup only for non-h2o isotope now   
   integer, pointer :: adsorbid(:)                                      !which tracer is adsorbed
   integer, pointer :: volatileid(:)
   integer, pointer :: h2oid(:)   
   character(len=12),pointer:: tracernames(:)                           !array with tracer names
   real(r8),pointer :: gram_mole_wt(:)                                  !molecular weight of the master species, [g/mol]
   real(r8),pointer :: vtrans_scal(:)                                   !scaling factor for plant tracer uptake through transpiration, for non-water neutral aqueous tracers
  
  contains
     procedure, public  :: Init
     procedure, public  :: init_scalars
     procedure, private :: InitAllocate
  end type BeTRtracer_type
  
  
   
  contains
  

  subroutine Init(this)

    implicit none
    class(BeTRtracer_type) :: this
    
    call this%InitAllocate()
  end subroutine Init  
!--------------------------------------------------------------------------------       
  subroutine init_scalars(this)
  !
  ! DESCRIPTIONS
  !
  implicit none
  class(BeTRtracer_type) :: this
  
  this%ntracers               = 0                                    ! total number of tracers, gas/aqueous tracers + solid tracers that undergo active mineral protection
  this%ngwmobile_tracers      = 0                                    ! total number of tracers undergoing gas/aqueous movement
  this%nvolatile_tracers      = 0                                    ! number of volatile_tracers
  this%nsolid_equil_tracers   = 0                                    ! number of tracers that undergo equilibrium adsorption in soil could include adsorbed doc, nh4(+)
  this%nsolid_passive_tracers = 0                                    ! number of tracers that undergo active mineral protection   
  this%nco2_tags              = 0                                    ! number of tagged co2 tracers
  this%nh2o_tracers           = 0                                    ! number of h2o tracers, this will be used to compute vapor gradient and thermal gradient driven isotopic flow
  this%is_oddstep             = .true.                               !this is not used now, originally was included to set up alternative numerical methods  


  this%id_trc_ch4         = 0       ! tag for methane
  this%id_trc_o2          = 0       ! tag for co2
  this%id_trc_n2          = 0       ! tag for n2
  this%id_trc_no          = 0       ! tag for no
  this%id_trc_n2o         = 0       ! tag for n2o
  this%id_trc_ar          = 0       ! tag for ar
  this%id_trc_air_co2x    = 0       ! tag for atmospheric co2
  this%id_trc_arrt_co2x   = 0       ! tag for autotrophic co2
  this%id_trc_hrsoi_co2x  = 0       ! tag for heterotrophic co2
   
  this%id_trc_co2x        = 0       ! tag for co2 and its related species, co2x(CO2, H2CO3, HCO3(-), CO3(2-)),
  this%id_trc_nh3x        = 0       ! tag for nh3 and its related species, nh3x(NH3, NH4OH,NH4(+))
  this%id_trc_no3x        = 0       ! tag for no3 and its related species, no3x(HNO3,NO3(-))
  this%id_trc_no2x        = 0       ! tag for no2 and its related species, no2x(HNO2,NO2(-))
  this%id_trc_dom         = 0       ! tag for generic dissolved organic matter

   
  this%id_trc_o18_h2o     = 0    ! tag for H2O(18)
  this%id_trc_o17_h2o     = 0    ! tag for H2O(17)
  this%id_trc_d_h2o       = 0    ! tag for DHO
  this%id_trc_c13_co2x    = 0    ! tag for C(13)O2 and its related species
  this%id_trc_c14_co2x    = 0    ! tag for C(14)O2 and its related species
  this%id_trc_o18_co2x    = 0    ! tag for O(18)CO and its related species
  this%id_trc_o17_co2x    = 0    ! tag for O(17)CO and its related species
  this%id_trc_o18_h2o_ice = 0    ! tag for H2O(18) in ice
  this%id_trc_d_h2o_ice   = 0    ! tag for HDO in ice
  this%id_trc_o18_o2      = 0    ! tag for O(18)O and its related species
  this%id_trc_o17_o2      = 0    ! tag for O(17)O and its related species
  
  this%betr_simname = ''
  end subroutine init_scalars
  
  
!--------------------------------------------------------------------------------       
  subroutine InitAllocate(this)
  !
  ! DESCRIPTIONS
  !
  
  implicit none
  class(BeTRtracer_type) :: this  
  integer, parameter :: nanid=-1
  
  allocate(this%is_volatile        (this%ngwmobile_tracers));    this%is_volatile(:)      = .false.
  allocate(this%is_adsorb          (this%ngwmobile_tracers));    this%is_adsorb(:)        = .false.
  allocate(this%is_advective       (this%ngwmobile_tracers));    this%is_advective(:)     = .false.  
  allocate(this%is_mobile          (this%ntracers));             this%is_mobile(:)        = .false.
  allocate(this%is_h2o             (this%ngwmobile_tracers));    this%is_h2o(:)           = .false.
  allocate(this%is_co2tag          (this%ngwmobile_tracers));    this%is_co2tag(:)        = .false.
  allocate(this%is_dom             (this%ngwmobile_tracers));    this%is_dom(:)           = .false.
  allocate(this%is_isotope         (this%ngwmobile_tracers));    this%is_isotope(:)       = .false.
  
  allocate(this%adsorbid           (this%ngwmobile_tracers));    this%adsorbid(:)           = nanid
  allocate(this%volatileid         (this%ngwmobile_tracers));    this%volatileid(:)         = nanid
  allocate(this%h2oid              (this%nh2o_tracers));         this%h2oid(:)              = nanid
  allocate(this%tracernames        (this%ntracers));             this%tracernames(:)        = ''
  allocate(this%vtrans_scal        (this%ngwmobile_tracers));    this%vtrans_scal(:)        = 0._r8   !no transport through xylem transpiration
  
  end subroutine InitAllocate



  
end module BeTRTracerType
