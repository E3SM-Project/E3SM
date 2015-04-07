module HumanIndexMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: HumanIndexMod
!
! !DESCRIPTION:
! Calculates Wetbulb Temperature, Stull Wet Bulb Temperature,
!      Heat Index, Apparent Temperature, Simplified Wet Bulb 
!      Globe Temperature, Humidex, Discomfort Index, Stull 
!      Discomfort Index, Temperature Humidity Comfort Index, 
!      Temperature Humidity Physiology Index, Swamp Cooler 
!      Temperature, Kelvin to Celsius, Vapor Pressure, & QSat_2
!
! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use decompMod            , only : bounds_type
! !PUBLIC TYPES:
  implicit none
  save
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: HumanIndexReadNML    ! Read in the namelist for HumanIndex
  public :: Wet_Bulb             ! Wet bulb temperature
  public :: Wet_BulbS            ! Wet bulb temperature from relative humidity
  public :: HeatIndex            ! Heat index
  public :: AppTemp              ! Apparant temperature
  public :: swbgt                ! Simplified Wetbulb Globe temperature
  public :: hmdex                ! humidex, human discomfort based on heat and humidity
  public :: dis_coi              ! Discomfort index
  public :: dis_coiS             ! Discomfort index from relative humidity
  public :: THIndex              ! Temperature humidity index
  public :: SwampCoolEff         ! Swamp Cooling efficiency
  public :: KtoC                 ! Convert Kelvin to Celcius
  public :: VaporPres            ! Vapor pressure
  public :: QSat_2               ! Saturation mix. ratio and the change in sat. mix rat. with respect to Temp

!
! !PUBLIC MEMBER DATA:
  logical, public :: calc_human_stress_indices = .true.   ! If should calculate the set of human stress indices
  type,    public :: humanindex_type   
     real(r8), pointer :: tc_ref2m_patch              (:) ! Patch 2 m height surface air temperature (C)
     real(r8), pointer :: vap_ref2m_patch             (:) ! Patch 2 m height vapor pressure (Pa)
     real(r8), pointer :: appar_temp_ref2m_patch      (:) ! Patch 2 m apparent temperature (C)
     real(r8), pointer :: appar_temp_ref2m_r_patch    (:) ! Patch Rural 2 m apparent temperature (C)
     real(r8), pointer :: swbgt_ref2m_patch           (:) ! Patch 2 m Simplified Wetbulb Globe temperature (C)
     real(r8), pointer :: swbgt_ref2m_r_patch         (:) ! Patch Rural 2 m Simplified Wetbulb Globe temperature (C)
     real(r8), pointer :: humidex_ref2m_patch         (:) ! Patch 2 m Humidex (C)
     real(r8), pointer :: humidex_ref2m_r_patch       (:) ! Patch Rural 2 m Humidex (C)
     real(r8), pointer :: wbt_ref2m_patch             (:) ! Patch 2 m Stull Wet Bulb temperature (C)
     real(r8), pointer :: wbt_ref2m_r_patch           (:) ! Patch Rural 2 m Stull Wet Bulb temperature (C)
     real(r8), pointer :: wb_ref2m_patch              (:) ! Patch 2 m Wet Bulb temperature (C)
     real(r8), pointer :: wb_ref2m_r_patch            (:) ! Patch Rural 2 m Wet Bulb temperature (C)
     real(r8), pointer :: teq_ref2m_patch             (:) ! Patch 2 m height Equivalent temperature (K)
     real(r8), pointer :: teq_ref2m_r_patch           (:) ! Patch Rural 2 m Equivalent temperature (K)
     real(r8), pointer :: ept_ref2m_patch             (:) ! Patch 2 m height Equivalent Potential temperature (K)
     real(r8), pointer :: ept_ref2m_r_patch           (:) ! Patch Rural 2 m height Equivalent Potential temperature (K)
     real(r8), pointer :: discomf_index_ref2m_patch   (:) ! Patch 2 m Discomfort Index temperature (C)
     real(r8), pointer :: discomf_index_ref2m_r_patch (:) ! Patch Rural 2 m Discomfort Index temperature (C)
     real(r8), pointer :: discomf_index_ref2mS_patch  (:) ! Patch 2 m height Discomfort Index Stull temperature (C)
     real(r8), pointer :: discomf_index_ref2mS_r_patch(:) ! Patch Rural 2 m Discomfort Index Stull temperature (K)
     real(r8), pointer :: nws_hi_ref2m_patch          (:) ! Patch 2 m NWS Heat Index (C)
     real(r8), pointer :: nws_hi_ref2m_r_patch        (:) ! Patch Rural 2 m NWS Heat Index (C)
     real(r8), pointer :: thip_ref2m_patch            (:) ! Patch 2 m Temperature Humidity Index Physiology (C)
     real(r8), pointer :: thip_ref2m_r_patch          (:) ! Patch Rural 2 m Temperature Humidity Index Physiology (C)
     real(r8), pointer :: thic_ref2m_patch            (:) ! Patch 2 m Temperature Humidity Index Comfort (C)
     real(r8), pointer :: thic_ref2m_r_patch          (:) ! Patch Rural 2 m Temperature Humidity Index Comfort (C)
     real(r8), pointer :: swmp65_ref2m_patch          (:) ! Patch 2 m Swamp Cooler temperature 65% effi (C)
     real(r8), pointer :: swmp65_ref2m_r_patch        (:) ! Patch Rural 2 m Swamp Cooler temperature 65% effi (C)
     real(r8), pointer :: swmp80_ref2m_patch          (:) ! Patch 2 m Swamp Cooler temperature 80% effi (C)
     real(r8), pointer :: swmp80_ref2m_r_patch        (:) ! Patch Rural 2 m Swamp Cooler temperature 80% effi (C)
     real(r8), pointer :: appar_temp_ref2m_u_patch    (:) ! Patch Urban 2 m apparent temperature (C)
     real(r8), pointer :: swbgt_ref2m_u_patch         (:) ! Patch Urban 2 m Simplified Wetbulb Globe temperature (C)
     real(r8), pointer :: humidex_ref2m_u_patch       (:) ! Patch Urban 2 m Humidex (C)
     real(r8), pointer :: wbt_ref2m_u_patch           (:) ! Patch Urban 2 m Stull Wet Bulb temperature (C)
     real(r8), pointer :: wb_ref2m_u_patch            (:) ! Patch Urban 2 m Wet Bulb temperature (C)
     real(r8), pointer :: teq_ref2m_u_patch           (:) ! Patch Urban 2 m Equivalent
     real(r8), pointer :: ept_ref2m_u_patch           (:) ! Patch Urban 2 m height Equivalent Potential temperature (K)
     real(r8), pointer :: discomf_index_ref2m_u_patch (:) !Urban 2 m Discomfort Index temperature (C)
     real(r8), pointer :: discomf_index_ref2mS_u_patch(:) !Urban 2 m Discomfort Index Stull temperature (K)
     real(r8), pointer :: nws_hi_ref2m_u_patch        (:) !Urban 2 m NWS Heat Index (C)
     real(r8), pointer :: thip_ref2m_u_patch          (:) !Urban 2 m Temperature Humidity Index Physiology (C)
     real(r8), pointer :: thic_ref2m_u_patch          (:) !Urban 2 m Temperature Humidity Index Comfort (C)
     real(r8), pointer :: swmp65_ref2m_u_patch        (:) !Urban 2 m Swamp Cooler temperature 65% effi (C)
     real(r8), pointer :: swmp80_ref2m_u_patch        (:) !Urban 2 m Swamp Cooler temperature 80% effi (C) temperature (K)
  contains
     procedure, public  :: Init                ! Public initialization
     procedure, private :: InitAllocate        ! Private allocation method
     procedure, private :: InitHistory         ! Private history setup method
  end type humanindex_type
!
! !REVISION HISTORY:
! Created by Jonathan R Buzan 03-07-12
! Modified 03-14-12--- filter routines for WB
!
! Modified 08-12-12--- filter for below zero calculation. 
!    Added WB = T at 0 and below
! Modified 05-13-13--- Adding additional Metrics. 
!    Added Apparent Temperature (Australian BOM)
!    Added Simplified Wetbulb Globe Temperature
!    Added Humidex
!    Added Discomfort Index
!    The previous Metrics were from Keith Oleson
!    Added Temperature Humidity Index
!    Added Swamp Cooler Efficiency
!
! Modified 05-16-13--- Added Current Vapor Pressure and 
!    Kelvin to Celsius and converted all
!    equations that use these inputs
! Modified 08-30-13--- Finalized Comments.  Added a new 
!    qsat algorithm.  Changed wet bulb calculations
!    to calculate over the large range of atmospheric 
!    conditions.  
! Modified 03-21-14--- Changed Specific Humidity to Mixing
!    Ratio.
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init
!
! !INTERFACE:
subroutine Init(this, bounds )
!
! !DESCRIPTION: Initialize human index object
!       
! !USES:
! !ARGUMENTS:
    implicit none
    class(humanindex_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
! !LOCAL VARIABLES:
    type(bounds_type) :: bounds_tmp
!EOP
!-----------------------------------------------------------------------
    if ( calc_human_stress_indices ) then
       call this%InitAllocate ( bounds )
       call this%InitHistory (  bounds )
    else
       ! Associate statements need humanindex_inst to be allocated
       ! So allocate with size 1 when not being used
       bounds_tmp%begp = 1
       bounds_tmp%endp = 1
       call this%InitAllocate ( bounds )
    end if

end subroutine Init

!------------------------------------------------------------------------
subroutine InitAllocate(this, bounds)
!
! !DESCRIPTION:
! Initialize module data structure
!
! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
!
! !ARGUMENTS:
    class(humanindex_type) :: this
    type(bounds_type), intent(in) :: bounds  
!
! !LOCAL VARIABLES:
    integer :: begp, endp
!------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    allocate(this%vap_ref2m_patch          (begp:endp))                      ; this%vap_ref2m_patch          (:)   = nan
    allocate(this%humidex_ref2m_u_patch    (begp:endp))                      ; this%humidex_ref2m_u_patch    (:)   = nan
    allocate(this%humidex_ref2m_patch      (begp:endp))                      ; this%humidex_ref2m_patch      (:)   = nan
    allocate(this%humidex_ref2m_r_patch    (begp:endp))                      ; this%humidex_ref2m_r_patch    (:)   = nan
    allocate(this%nws_hi_ref2m_patch       (begp:endp))                      ; this%nws_hi_ref2m_patch       (:)   = nan
    allocate(this%nws_hi_ref2m_r_patch     (begp:endp))                      ; this%nws_hi_ref2m_r_patch     (:)   = nan
    allocate(this%thip_ref2m_patch         (begp:endp))                      ; this%thip_ref2m_patch         (:)   = nan
    allocate(this%thip_ref2m_r_patch       (begp:endp))                      ; this%thip_ref2m_r_patch       (:)   = nan
    allocate(this%thic_ref2m_patch         (begp:endp))                      ; this%thic_ref2m_patch         (:)   = nan
    allocate(this%thic_ref2m_r_patch       (begp:endp))                      ; this%thic_ref2m_r_patch       (:)   = nan
    allocate(this%nws_hi_ref2m_u_patch     (begp:endp))                      ; this%nws_hi_ref2m_u_patch     (:)   = nan
    allocate(this%thip_ref2m_u_patch       (begp:endp))                      ; this%thip_ref2m_u_patch       (:)   = nan
    allocate(this%thic_ref2m_u_patch       (begp:endp))                      ; this%thic_ref2m_u_patch       (:)   = nan
    allocate(this%tc_ref2m_patch            (begp:endp))                      ; this%tc_ref2m_patch             (:)   = nan
    allocate(this%appar_temp_ref2m_patch    (begp:endp))                      ; this%appar_temp_ref2m_patch     (:)   = nan
    allocate(this%appar_temp_ref2m_r_patch  (begp:endp))                      ; this%appar_temp_ref2m_r_patch   (:)   = nan
    allocate(this%swbgt_ref2m_patch         (begp:endp))                      ; this%swbgt_ref2m_patch          (:)   = nan
    allocate(this%swbgt_ref2m_r_patch       (begp:endp))                      ; this%swbgt_ref2m_r_patch        (:)   = nan
    allocate(this%wbt_ref2m_patch           (begp:endp))                      ; this%wbt_ref2m_patch            (:)   = nan
    allocate(this%wbt_ref2m_r_patch         (begp:endp))                      ; this%wbt_ref2m_r_patch          (:)   = nan
    allocate(this%wb_ref2m_patch            (begp:endp))                      ; this%wb_ref2m_patch             (:)   = nan
    allocate(this%wb_ref2m_r_patch          (begp:endp))                      ; this%wb_ref2m_r_patch           (:)   = nan
    allocate(this%teq_ref2m_patch           (begp:endp))                      ; this%teq_ref2m_patch            (:)   = nan
    allocate(this%teq_ref2m_r_patch         (begp:endp))                      ; this%teq_ref2m_r_patch          (:)   = nan
    allocate(this%ept_ref2m_patch           (begp:endp))                      ; this%ept_ref2m_patch            (:)   = nan
    allocate(this%ept_ref2m_r_patch         (begp:endp))                      ; this%ept_ref2m_r_patch          (:)   = nan
    allocate(this%discomf_index_ref2m_patch (begp:endp))                      ; this%discomf_index_ref2m_patch  (:)   = nan
    allocate(this%discomf_index_ref2m_r_patch(begp:endp))                     ; this%discomf_index_ref2m_r_patch(:)   = nan
    allocate(this%discomf_index_ref2mS_patch(begp:endp))                      ; this%discomf_index_ref2mS_patch (:)   = nan
    allocate(this%discomf_index_ref2mS_r_patch(begp:endp))                    ; this%discomf_index_ref2mS_r_patch(:)  = nan
    allocate(this%discomf_index_ref2mS_u_patch(begp:endp))                    ; this%discomf_index_ref2mS_u_patch(:)  = nan
    allocate(this%swmp65_ref2m_patch        (begp:endp))                      ; this%swmp65_ref2m_patch         (:)   = nan
    allocate(this%swmp65_ref2m_r_patch      (begp:endp))                      ; this%swmp65_ref2m_r_patch       (:)   = nan
    allocate(this%swmp80_ref2m_patch        (begp:endp))                      ; this%swmp80_ref2m_patch         (:)   = nan
    allocate(this%swmp80_ref2m_r_patch      (begp:endp))                      ; this%swmp80_ref2m_r_patch       (:)   = nan
    allocate(this%swmp80_ref2m_u_patch      (begp:endp))                      ; this%swmp80_ref2m_u_patch       (:)   = nan
    allocate(this%appar_temp_ref2m_u_patch  (begp:endp))                      ; this%appar_temp_ref2m_u_patch   (:)   = nan
    allocate(this%swbgt_ref2m_u_patch       (begp:endp))                      ; this%swbgt_ref2m_u_patch        (:)   = nan
    allocate(this%wbt_ref2m_u_patch         (begp:endp))                      ; this%wbt_ref2m_u_patch          (:)   = nan
    allocate(this%wbt_ref2m_u_patch         (begp:endp))                      ; this%wbt_ref2m_u_patch          (:)   = nan
    allocate(this%wb_ref2m_u_patch          (begp:endp))                      ; this%wb_ref2m_u_patch           (:)   = nan
    allocate(this%teq_ref2m_u_patch         (begp:endp))                      ; this%teq_ref2m_u_patch          (:)   = nan
    allocate(this%ept_ref2m_u_patch         (begp:endp))                      ; this%ept_ref2m_u_patch          (:)   = nan
    allocate(this%discomf_index_ref2m_u_patch(begp:endp))                     ; this%discomf_index_ref2m_u_patch(:)   = nan
    allocate(this%swmp65_ref2m_u_patch      (begp:endp))                      ; this%swmp65_ref2m_u_patch       (:)   = nan
end subroutine InitAllocate


!------------------------------------------------------------------------
subroutine InitHistory(this, bounds)
!
! !DESCRIPTION:
! Initialize history data
!
! !USES:
    use clm_varcon      , only : spval
    use histFileMod     , only : hist_addfld1d
!
! !ARGUMENTS:
    class(humanindex_type)        :: this
    type(bounds_type), intent(in) :: bounds  
!
! !LOCAL VARIABLES:
    integer :: begp, endp
!------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%appar_temp_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='APPAR_TEMP', units='C',  &
            avgflag='A', long_name='2 m apparent temperature', &
            ptr_patch=this%appar_temp_ref2m_patch)

    this%appar_temp_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='APPAR_TEMP_U', units='C',  &
            avgflag='A', long_name='Urban 2 m apparent temperature', &
            ptr_patch=this%appar_temp_ref2m_u_patch, set_nourb=spval)

    this%appar_temp_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='APPAR_TEMP_R', units='C',  &
            avgflag='A', long_name='Rural 2 m apparent temperature', &
            ptr_patch=this%appar_temp_ref2m_r_patch, set_spec=spval)

    this%swbgt_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWBGT', units='C',  &
            avgflag='A', long_name='2 m Simplified Wetbulb Globe Temp', &
            ptr_patch=this%swbgt_ref2m_patch)

    this%swbgt_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWBGT_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Simplified Wetbulb Globe Temp', &
            ptr_patch=this%swbgt_ref2m_u_patch, set_nourb=spval)

    this%swbgt_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWBGT_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Simplified Wetbulb Globe Temp', &
            ptr_patch=this%swbgt_ref2m_r_patch, set_spec=spval)

    this%humidex_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='HUMIDEX', units='C',  &
            avgflag='A', long_name='2 m Humidex', &
            ptr_patch=this%humidex_ref2m_patch)

    this%humidex_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='HUMIDEX_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Humidex', &
            ptr_patch=this%humidex_ref2m_u_patch, set_nourb=spval)

    this%humidex_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='HUMIDEX_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Humidex', &
            ptr_patch=this%humidex_ref2m_r_patch, set_spec=spval)

    this%wbt_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='WBT', units='C',  &
            avgflag='A', long_name='2 m Stull Wet Bulb', &
            ptr_patch=this%wbt_ref2m_u_patch, set_nourb=spval)

    this%wbt_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='WBT_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Stull Wet Bulb', &
            ptr_patch=this%wbt_ref2m_r_patch, set_spec=spval)

    this%wb_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='WBA', units='C',  &
            avgflag='A', long_name='2 m Wet Bulb', &
            ptr_patch=this%wb_ref2m_patch)

    this%wb_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='WBA_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Wet Bulb', &
            ptr_patch=this%wb_ref2m_u_patch, set_nourb=spval)

    this%wb_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='WBA_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Wet Bulb', &
            ptr_patch=this%wb_ref2m_r_patch, set_spec=spval)

    this%teq_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEQ', units='K',  &
            avgflag='A', long_name='2 m Equiv Temp', &
            ptr_patch=this%teq_ref2m_patch)

    this%teq_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEQ_U', units='K',  &
            avgflag='A', long_name='Urban 2 m Equiv Temp', &
            ptr_patch=this%teq_ref2m_u_patch, set_nourb=spval)

    this%teq_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEQ_R', units='K',  &
            avgflag='A', long_name='Rural 2 m Equiv Temp', &
            ptr_patch=this%teq_ref2m_r_patch, set_spec=spval)

    this%ept_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='EPT', units='K',  &
            avgflag='A', long_name='2 m Equiv Pot Temp', &
            ptr_patch=this%ept_ref2m_patch)

    this%ept_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='EPT_U', units='K',  &
            avgflag='A', long_name='Urban 2 m Equiv Pot Temp', &
            ptr_patch=this%ept_ref2m_u_patch, set_nourb=spval)

    this%ept_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='EPT_R', units='K',  &
            avgflag='A', long_name='Rural 2 m Equiv Pot Temp', &
            ptr_patch=this%ept_ref2m_r_patch, set_spec=spval)

    this%discomf_index_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISCOI', units='C',  &
            avgflag='A', long_name='2 m Discomfort Index', &
            ptr_patch=this%discomf_index_ref2m_patch)

    this%discomf_index_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISCOI_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Discomfort Index', &
            ptr_patch=this%discomf_index_ref2m_u_patch, set_nourb=spval)

    this%discomf_index_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISCOI_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Discomfort Index', &
            ptr_patch=this%discomf_index_ref2m_r_patch, set_spec=spval)

    this%discomf_index_ref2mS_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISCOIS', units='C',  &
            avgflag='A', long_name='2 m Stull Discomfort Index', &
            ptr_patch=this%discomf_index_ref2mS_patch)

    this%discomf_index_ref2mS_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISCOIS_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Stull Discomfort Index', &
            ptr_patch=this%discomf_index_ref2mS_u_patch, set_nourb=spval)

    this%discomf_index_ref2mS_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISCOIS_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Stull Discomfort Index', &
            ptr_patch=this%discomf_index_ref2mS_r_patch, set_spec=spval)

    this%nws_hi_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='HIA', units='C',  &
            avgflag='A', long_name='2 m NWS Heat Index', &
            ptr_patch=this%nws_hi_ref2m_patch)

    this%nws_hi_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='HIA_U', units='C',  &
            avgflag='A', long_name='Urban 2 m NWS Heat Index', &
            ptr_patch=this%nws_hi_ref2m_u_patch, set_nourb=spval)

    this%nws_hi_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='HIA_R', units='C',  &
            avgflag='A', long_name='Rural 2 m NWS Heat Index', &
            ptr_patch=this%nws_hi_ref2m_r_patch, set_spec=spval)

    this%thip_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='THIP', units='C',  &
            avgflag='A', long_name='2 m Temp Hum Index Physiology', &
            ptr_patch=this%thip_ref2m_patch)

    this%thip_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='THIP_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Temp Hum Index Physiology', &
            ptr_patch=this%thip_ref2m_u_patch, set_nourb=spval)

    this%thip_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='THIP_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Temp Hum Index Physiology', &
            ptr_patch=this%thip_ref2m_r_patch, set_spec=spval)

    this%thic_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='THIC', units='C',  &
            avgflag='A', long_name='2 m Temp Hum Index Comfort', &
            ptr_patch=this%thic_ref2m_patch)

    this%thic_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='THIC_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Temp Hum Index Comfort', &
            ptr_patch=this%thic_ref2m_u_patch, set_nourb=spval)

    this%thic_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='THIC_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Temp Hum Index Comfort', &
            ptr_patch=this%thic_ref2m_r_patch, set_spec=spval)

    this%swmp65_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWMP65', units='C',  &
            avgflag='A', long_name='2 m Swamp Cooler Temp 65% Eff', &
            ptr_patch=this%swmp65_ref2m_patch)

    this%swmp65_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWMP65_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Swamp Cooler Temp 65% Eff', &
            ptr_patch=this%swmp65_ref2m_u_patch, set_nourb=spval)

    this%swmp65_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWMP65_R', units='C',  &
            avgflag='A', long_name='Rural 2 m Swamp Cooler Temp 65% Eff', &
            ptr_patch=this%swmp65_ref2m_r_patch, set_spec=spval)

    this%swmp80_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWMP80', units='C',  &
            avgflag='A', long_name='2 m Swamp Cooler Temp 80% Eff', &
            ptr_patch=this%swmp80_ref2m_patch)

    this%swmp80_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWMP80_U', units='C',  &
            avgflag='A', long_name='Urban 2 m Swamp Cooler Temp 80% Eff', &
            ptr_patch=this%swmp80_ref2m_u_patch, set_nourb=spval)

    this%swmp80_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='SWMP80_R', units='C',  &
               avgflag='A', long_name='Rural 2 m Swamp Cooler Temp 80% Eff', &
               ptr_patch=this%swmp80_ref2m_r_patch, set_spec=spval)

end subroutine InitHistory

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HumanIndexReadNML
!
! !INTERFACE:
  subroutine HumanIndexReadNML( NLFilename )
!
! !DESCRIPTION:
!       
! !USES:
    use shr_mpi_mod   , only : shr_mpi_bcast
    use abortutils    , only : endrun
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use clm_varctl    , only : iulog
    use shr_log_mod   , only : errMsg => shr_log_errMsg
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'UrbanReadNML'  ! subroutine name
!EOP
!-----------------------------------------------------------------------
    namelist / clm_humanindex_inparm / calc_human_stress_indices

    ! ----------------------------------------------------------------------
    ! Read namelist from input namelist filename
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm_humanindex_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_humanindex_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_humanindex_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_humanindex_inparm namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(calc_human_stress_indices, mpicom)

  end subroutine HumanIndexReadNML

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: AppTemp
!
! !INTERFACE:
  subroutine AppTemp (Tc_1, vap_pres, u10_m, app_temp)
!
! !DESCRIPTION:
! Apparent Temperature (Australian BOM): Here we use equation 22 
!    where AT is a function of air temperature (C), water 
!    vapor pressure (kPa), and 10-m wind speed (m/s). vap_pres
!    from Erich Fischer (consistent with CLM equations)
!
! Reference:  Steadman, R.G., 1994: Norms of apparent temperature
!             in Australia, Aust. Met. Mag., 43, 1-16. 
!       
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_1     ! temperature (C)
    real(r8), intent(in)  :: vap_pres ! Vapor Pressure (pa)
    real(r8), intent(in)  :: u10_m    ! Winds at 10m (m/s)
    real(r8), intent(out) :: app_temp ! Apparent Temperature (C)
!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
    app_temp = Tc_1 + 3.30_r8*vap_pres/1000._r8 - 0.70_r8*u10_m - 4.0_r8

  end subroutine AppTemp
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: swbgt
!
! !INTERFACE:
  subroutine swbgt (Tc_2, vap_pres, s_wbgt)
!
! !DESCRIPTION:
! Simplified Wet Bulb Globe Temperature: 
!      Requires air temperature (C), water vapor pressure (hPa)
!
! Reference:  Willett, K.M., and S. Sherwood, 2010: Exceedance of heat
!       index thresholds for 15 regions under a warming 
!       climate using the wet-bulb globe temperature,
!       Int. J. Climatol., doi:10.1002/joc.2257
!       
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_2     ! temperature (C)
    real(r8), intent(in)  :: vap_pres ! Vapor Pressure (pa)
    real(r8), intent(out) :: s_wbgt   ! Simplified Wet Bulb Globe Temperature (C)

!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
    s_wbgt = 0.567_r8*(Tc_2)  + 0.393_r8*vap_pres/100._r8 + 3.94_r8

  end subroutine swbgt
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hmdex
!
! !INTERFACE:
  subroutine hmdex (Tc_3, vap_pres, humidex)
!
! !DESCRIPTION:
! Humidex:
!     Requires air temperature (C), water vapor pressure (hPa)
! Reference:  Masterson, J., and F. Richardson, 1979: Humidex, a 
!      method of quantifying human discomfort due to 
!      excessive heat and humidity, CLI 1-79, Environment 
!      Canada, Atmosheric Environment Service, Downsview, Ontario
!      
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_3     ! temperature (C)
    real(r8), intent(in)  :: vap_pres ! Vapor Pressure (Pa)
    real(r8), intent(out) :: humidex  ! Humidex (C)

!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
    humidex = Tc_3 + ((5._r8/9._r8) * (vap_pres/100._r8 - 10._r8))

  end subroutine hmdex
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: dis_coi
!
! !INTERFACE:
  subroutine dis_coi (Tc_4, wb_t, discoi)
!
! !DESCRIPTION:
! Discomfort Index
!      The wet bulb temperature is from Davies-Jones, 2008.
!      Requires air temperature (C), wet bulb temperature (C) 
! Reference:  Epstein, Y., and D.S. Moran, 2006: Thermal comfort and the heat stress indices,
!      Ind. Health, 44, 388-398.
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_4     ! temperature (C)
    real(r8), intent(in)  :: wb_t     ! Wet Bulb Temperature (C)
    real(r8), intent(out) :: discoi   ! Discomfort Index (C)

!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
    discoi = 0.5_r8*wb_t + 0.5_r8*Tc_4

  end subroutine dis_coi
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: dis_coiS
!
! !INTERFACE:
  subroutine dis_coiS (Tc_5, relhum, wbt_s, discois)
!
! !DESCRIPTION:
! Discomfort Index
! The wet bulb temperature is from Stull, 2011.
!       Requires air temperature (C), wet bulb temperature (C) 
! Reference:  Epstein, Y., and D.S. Moran, 2006: Thermal comfort and the heat stress indices,
!       Ind. Health, 44, 388-398.
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_5      ! temperature (C)
    real(r8), intent(in)  :: wbt_s     ! Wet Bulb Temperature (C)
    real(r8), intent(in)  :: relhum    ! Relative Humidity (%)
    real(r8), intent(out) :: discois   ! Discomfort Index (C)

!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) ::  Tc                    ! 2-m temperature with limit (C)
    real(r8) ::  rh                    ! 2-m relative humidity with limit (%)
    real(r8) ::  rh_min                ! Minimum 2-m relative humidity (%)

!
!-----------------------------------------------------------------------
    Tc = min(Tc_5,50._r8)
    rh = min(relhum,99._r8)
    rh = max(rh,5._r8)
    rh_min = Tc*(-2.27_r8)+27.7_r8
    if (Tc < -20._r8 .or. rh < rh_min) then
       ! wbt_s calculation invalid
       discois = Tc
    else
       ! wbt_s calculation valid
       discois = 0.5_r8*wbt_s + 0.5_r8*Tc
    end if

  end subroutine dis_coiS
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Wet_Bulb
!
! !INTERFACE:
  subroutine Wet_Bulb (Tin_1,vape,pin,relhum,qin,Teq,epott,wb_it)

!
! !DESCRIPTION:
! Calculates Wet Bulb Temperature, Theta_wb, Theta_e, Moist Pot Temp, 
!       Lifting Cond Temp, and Equiv Temp using Davies-Jones 2008 Method.
!       1st calculates the lifting cond temperature (Bolton 1980 eqn 22).  
!       Then calculates the moist pot temp (Bolton 1980 eqn 24). Then 
!       calculates Equivalent Potential Temperature (Bolton 1980 eqn 39).  
!       From equivalent pot temp, equiv temp and Theta_w (Davies-Jones 
!       2008 eqn 3.5-3.8).  An accurate 'first guess' of wet bulb temperature
!       is determined (Davies-Jones 2008 eqn 4.8-4.11). Newton-Raphson
!       is used for 2 iterations, determining final wet bulb temperature 
!       (Davies-Jones 2008 eqn 2.6).
! Requires Temperature,Vapor Pressure,Atmospheric Pressure,Relative Humidity,Mixing Ratio
! Reference:  Bolton: The computation of equivalent potential temperature. 
!       Monthly weather review (1980) vol. 108 (7) pp. 1046-1053
!       Davies-Jones: An efficient and accurate method for computing the 
!       wet-bulb temperature along pseudoadiabats. Monthly Weather Review 
!       (2008) vol. 136 (7) pp. 2764-2785
!       Flatau et al: Polynomial fits to saturation vapor pressure. 
!       Journal of Applied Meteorology (1992) vol. 31 pp. 1507-1513
! Note: Pressure needs to be in mb, mixing ratio needs to be in 
! kg/kg in some equations, and in g/kg in others.  
! Calculates Iteration via Newton-Raphson Method.  Only 2 iterations.
! Reference:  Davies-Jones: An efficient and accurate method for computing the 
!       wet-bulb temperature along pseudoadiabats. Monthly Weather Review 
!       (2008) vol. 136 (7) pp. 2764-2785
!       Flatau et al: Polynomial fits to saturation vapor pressure. 
!       Journal of Applied Meteorology (1992) vol. 31 pp. 1507-1513
! Note: Pressure needs to be in mb, mixing ratio needs to be in 
!       kg/kg in some equations. 
!
! !REVISION HISTORY:
!
! Created by Jonathan R Buzan 03-07-12
! Modified JRBuzan 06-29-13:  Major Revision.  Changes all Calculations to be based
!       upon Bolton eqn 39.  Uses Derivatives in Davies-Jones
!       2008 for calculation of vapor pressure.
! Modified JRBuzan 03-21-14:  Minor Revision.  Changed specific humidity to mixing
!       ratio.
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: Tin_1  ! 2-m air temperature (K)
    real(r8), intent(in) :: vape   ! Vapor Pressure (Pa)
    real(r8), intent(in) :: pin    ! Atmospheric Pressure (Pa)
    real(r8), intent(in) :: relhum ! Relative Humidity (%)
    real(r8), intent(in) :: qin    ! Specific Humidity (kg/kg)

    real(r8), intent(out) :: Teq   ! Equivalent Temperature (K)
    real(r8), intent(out) :: epott ! Equivalent Potential Temperature (K)
    real(r8), intent(out) :: wb_it ! Constant used for extreme cold temparatures (K)

!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: k1                  ! Quadratic Parameter (C)
    real(r8) :: k2                  ! Quadratic Parameter scaled by X (C) 
    real(r8) :: pmb                 ! Atmospheric Surface Pressure (mb)
    real(r8) :: D                   ! Linear Interpolation of X

    real(r8) :: constA = 2675._r8   ! Constant used for extreme cold temparatures (K)
    real(r8) :: grms = 1000._r8     ! Gram per Kilogram (g/kg)
    real(r8) :: p0 = 1000._r8       ! surface pressure (mb)
    real(r8) :: C                   ! Temperature of Freezing (K)

    real(r8) :: hot                 ! Dimensionless Quantity used for changing temperature regimes
    real(r8) :: cold                ! Dimensionless Quantity used for changing temperature regimes    

    real(r8) :: kappad = 0.2854_r8  ! Heat Capacity
    real(r8) :: T1                  ! Temperature (K)
    real(r8) :: vapemb              ! Vapor Pressure (mb)
    real(r8) :: mixr                ! Mixing Ratio (g/kg)

    real(r8) :: es_mb_teq           ! saturated vapour pressure for wrt TEQ (mb)
    real(r8) :: de_mbdTeq           ! Derivative of Saturated Vapor pressure wrt TEQ (mb/K)
    real(r8) :: dlnes_mbdTeq        ! Log derivative of the sat. vap pressure wrt TEQ (mb/K)
    real(r8) :: rs_teq              ! Mixing Ratio wrt TEQ (kg/kg)
    real(r8) :: rsdTeq              ! Derivative of Mixing Ratio wrt TEQ (kg/kg/K)
    real(r8) :: foftk_teq           ! Function of EPT wrt TEQ 
    real(r8) :: fdTeq               ! Derivative of Function of EPT wrt TEQ 

    real(r8) :: wb_temp             ! Wet Bulb Temperature First Guess (C)
    real(r8) :: es_mb_wb_temp       ! Vapor Pressure wrt Wet Bulb Temp (mb)
    real(r8) :: de_mbdwb_temp       ! Derivative of Sat. Vapor Pressure wrt WB Temp (mb/K)
    real(r8) :: dlnes_mbdwb_temp    ! Log Derivative of sat. vap. pressure wrt WB Temp (mb/K)
    real(r8) :: rs_wb_temp          ! Mixing Ratio wrt WB Temp (kg/kg)
    real(r8) :: rsdwb_temp          ! Derivative of Mixing Ratio wrt WB Temp (kg/kg/K)
    real(r8) :: foftk_wb_temp       ! Function of EPT wrt WB Temp
    real(r8) :: fdwb_temp           ! Derivative of function of EPT wrt WB Temp

    real(r8) :: tl                  ! Lifting Condensation Temperature (K)
    real(r8) :: theta_dl            ! Moist Potential Temperature (K)
    real(r8) :: pi                  ! Non dimensional Pressure
    real(r8) :: X                   ! Ratio of equivalent temperature to freezing scaled by Heat Capacity

    integer  :: j                   ! Iteration Step Number
!-----------------------------------------------------------------------

    C = SHR_CONST_TKFRZ             ! Freezing Temperature
    pmb = pin*0.01_r8               ! pa to mb
    vapemb = vape*0.01_r8           ! pa to mb
    T1 = Tin_1                      ! Use holder for T
    mixr = qin/(1._r8 - qin) * grms ! change specific humidity to mixing ratio (g/kg)
                  
    ! Calculate Equivalent Pot. Temp (pmb, T, mixing ratio (g/kg), pott, epott)
    ! Calculate Parameters for Wet Bulb Temp (epott, pmb)
    pi = (pmb/p0)**(kappad)
    D = (0.1859_r8*pmb/p0 + 0.6512)**(-1._r8)
    k1 = -38.5_r8*pi*pi +137.81_r8*pi -53.737_r8
    k2 = -4.392_r8*pi*pi +56.831_r8*pi -0.384_r8

    ! Calculate lifting condensation level.  first eqn 
    ! uses vapor pressure (mb)
    ! 2nd eqn uses relative humidity.  
    ! first equation: Bolton 1980 Eqn 21.
    !   tl = (2840._r8/(3.5_r8*log(T1) - log(vapemb) - 4.805_r8)) + 55._r8
    ! second equation: Bolton 1980 Eqn 22.  relhum = relative humidity
    tl = (1._r8/((1._r8/((T1 - 55._r8))) - (log(relhum/100._r8)/2840._r8))) + 55._r8

    ! Theta_DL: Bolton 1980 Eqn 24.
    theta_dl = T1*((p0/(pmb-vapemb))**kappad)*((T1/tl)**(mixr*0.00028_r8))

    ! EPT: Bolton 1980 Eqn 39.  
    epott = theta_dl*exp(((3.036_r8/tl)-0.00178_r8)*mixr*(1._r8 + 0.000448_r8*mixr))
    Teq = epott*pi    ! Equivalent Temperature at pressure
    X = (C/Teq)**3.504_r8
   
    ! Calculates the regime requirements of wet bulb equations.
    if (Teq > 355.15_r8) then
       hot = 1.0_r8
    else
       hot = 0.0_r8
    endif

    if ((X >= 1._r8) .AND. (X <= D)) then
       cold = 0._r8
    else
       cold = 1._r8
    endif

    ! Calculate Wet Bulb Temperature, initial guess
    ! Extremely cold regime if X.gt.D then need to 
    ! calculate dlnesTeqdTeq 
    if (X > D) then
       call QSat_2(Teq, pin, es_mb_teq, de_mbdTeq, dlnes_mbdTeq, rs_teq, rsdTeq, foftk_teq, fdTeq)
       wb_temp = Teq - C - ((constA*rs_teq)/(1._r8 + (constA*rs_teq*dlnes_mbdTeq)))
    else
       wb_temp = k1 - 1.21_r8 * cold - 1.45_r8 * hot - (k2 - 1.21_r8 * cold) * X + (0.58_r8 / X) * hot
    endif

    ! Newton-Raphson Method  2 iteration
    ! May need to put in a second iteration.  Probably best with a do loop.
    do j = 0, 1
       call QSat_2(wb_temp+C, pin, es_mb_wb_temp, de_mbdwb_temp, dlnes_mbdwb_temp, &
           rs_wb_temp, rsdwb_temp, foftk_wb_temp, fdwb_temp)
       wb_temp = wb_temp - ((foftk_wb_temp - X)/fdwb_temp)
       wb_it = wb_temp
    end do

  end subroutine Wet_Bulb
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Wet_BulbS
!
! !INTERFACE:
  subroutine Wet_BulbS (Tc_6,rh,wbt)

!
! !DESCRIPTION:
! Reference:  Stull, R., 2011: Wet-bulb temperature from relative humidity
!       and air temperature, J. Appl. Meteor. Climatol., doi:10.1175/JAMC-D-11-0143.1
!       Note: Requires air temperature (C) and relative humidity (%)
! Note: Pressure needs to be in mb, mixing ratio needs to be in 
!       kg/kg in some equations. 
! !REVISION HISTORY:
! Created by Jonathan R Buzan 03-07-12
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: Tc_6 ! Temperature (C)
    real(r8), intent(in) :: rh   ! Relative Humidity (%)
    real(r8), intent(out) :: wbt ! Wet Bulb Temperature (C)
!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!
    wbt = Tc_6 * atan(0.151977_r8*sqrt(rh + 8.313659_r8)) + &
          atan(Tc_6+rh) - atan(rh-1.676331_r8) + &
          0.00391838_r8*rh**(3._r8/2._r8)*atan(0.023101_r8*rh) - &
          4.686035_r8

  end subroutine Wet_BulbS
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HeatIndex
!
! !INTERFACE:
  subroutine HeatIndex (Tc_7, rh, hi)
!
! !DESCRIPTION:
! National Weather Service Heat Index 
! Requires air temperature (F), relative humidity (%)
! Valid for air temperatures above 20C. If below this set heatindex to air temperature.
! Reference: Steadman. The assessment of sultriness. Part I: 
!      A temperature-humidity index based on human physiology
!      and clothing science. J Appl Meteorol (1979) vol. 18 (7) pp. 861-873
!      Lans P. Rothfusz. "The heat index 'equation' (or
!      more than you ever wanted to know about heat index)", 
!      Scientific Services Division (NWS Southern Region Headquarters), 1 July 1990
! !REVISION HISTORY:
! Created by Jonathan R Buzan 03-07-12
! Modified JRBuzan 03-10-12
! Modified JRBuzan 05-14-13:  removed testing algorithm
!             Switched output to Celsius
!             Used Boundary Conditions from 
!             Keith Oleson
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_7     ! temperature (C)
    real(r8), intent(in)  :: rh       ! relative humidity (%)
    real(r8), intent(out) :: hi       ! Heat Index (C)
!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: tf
!
!-----------------------------------------------------------------------
    tf = (Tc_7) * 9._r8/5._r8 + 32._r8    ! fahrenheit

    if (tf < 68._r8) then
       hi = tf
    else
       hi = -42.379_r8 + 2.04901523_r8*tf                     &
                       + 10.14333127_r8*rh                    &
                       + (-0.22475541_r8*tf*rh)               &
                       + (-6.83783e-3_r8*tf**2._r8)           &
                       + (-5.481717e-2_r8*rh**2._r8)          &
                       + 1.22874e-3_r8*(tf**2._r8)*rh         &
                       + 8.5282e-4_r8*tf*rh**2._r8            &
                       + (-1.99e-6_r8*(tf**2._r8)*(rh**2._r8))
    endif
    hi = (hi - 32._r8) * 5._r8/9._r8     ! Celsius

  end subroutine HeatIndex
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: THIndex
!
! !INTERFACE:
  subroutine THIndex (Tc_8, wb_t, thic, thip)
!
! !DESCRIPTION:
! Temperature Humidity Index
! The wet bulb temperature is Davies-Jones 2008 (subroutine WetBulb)
! Requires air temperature (C), wet bulb temperature (C) 
! Calculates two forms of the index:  Comfort and Physiology
! Reference:  NWSCR (1976): Livestock hot weather stress. 
!      Regional operations manual letter C-31-76. 
!      National Weather Service Central Region, USA
!      Ingram: Evaporative cooling in the pig. Nature (1965)
! !REVISION HISTORY:
! Created by Jonathan R Buzan 03-15-13
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_8     ! temperature (C)
    real(r8), intent(in)  :: wb_t     ! Wet Bulb Temperature (C)
    real(r8), intent(out) :: thic     ! Temperature Humidity Index Comfort (C)
    real(r8), intent(out) :: thip     ! Temperature Humidity Index Physiology (C)

!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!
!    real(r8) :: 

!
!-----------------------------------------------------------------------
    thic = 0.72_r8*wb_t + 0.72_r8*(Tc_8) + 40.6_r8
    thip = 0.63_r8*wb_t + 1.17_r8*(Tc_8) + 32._r8

  end subroutine THIndex
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SwampCoolEff
!
! !INTERFACE:
  subroutine SwampCoolEff (Tc_9, wb_t, tswmp80, tswmp65)
!
! !DESCRIPTION:
! Swamp Cooler Efficiency
!       The wet bulb temperature is Davies-Jones 2008 (subroutine WetBulb)
!       Requires air temperature (C), wet bulb temperature (C) 
!       Assumes that the Swamp Cooler Efficiency 80% (properly maintained)
!       and 65% (improperly maintained).  
! Reference:  Koca et al: Evaporative cooling pads: test 
!       procedure and evaluation. Applied engineering
!       in agriculture (1991) vol. 7
! !REVISION HISTORY:
! Created by Jonathan R Buzan 03-15-13
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: Tc_9     ! temperature (C)
    real(r8), intent(in)  :: wb_t     ! Wet Bulb Temperature (C)
    real(r8), intent(out) :: tswmp80  ! Swamp Cooler Temp 80% Efficient (C)
    real(r8), intent(out) :: tswmp65  ! Swamp Cooler Temp 65% Efficient (C)

!
! !CALLED FROM:
! subroutine LakeFluxes in module LakeFluxesMod
! subroutine CanopyFluxes in module CanopyFluxesMod
! subroutine UrbanFluxes in module UrbanFluxesMod
! subroutine BareGroundFluxes in module BareGroundFluxesMod
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: neu80 = 0.80_r8  ! 80% Efficient
    real(r8) :: neu65 = 0.65_r8  ! 65% Efficient

!
!-----------------------------------------------------------------------
    tswmp80 = Tc_9 - neu80*(Tc_9 - wb_t)
    tswmp65 = Tc_9 - neu65*(Tc_9 - wb_t)

  end subroutine SwampCoolEff
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: KtoC
!
! !INTERFACE:
  subroutine KtoC (T_k, T_c)
!
! !DESCRIPTION:
! Converts Kelvins to Celsius
! !REVISION HISTORY:
! Created by Jonathan R Buzan 03-16-13
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T_k        ! temperature (K)
    real(r8), intent(out) :: T_c        ! temperature (C)

!
! !CALLED FROM:
! subroutines within this module
!
! !LOCAL VARIABLES:
!EOP
!
!    real(r8) :: 

!
!-----------------------------------------------------------------------
    T_c = T_k - SHR_CONST_TKFRZ

  end subroutine KtoC
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: VaporPres
!
! !INTERFACE:
  subroutine VaporPres (rh, e, erh)
!
! !DESCRIPTION:
! Calculates Vapor Pressure
!      Vapor Pressure from Erich Fischer (consistent with CLM 
!      equations, Keith Oleson)
! !REVISION HISTORY:
! Created by Jonathan R Buzan 03-16-13
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: rh        ! Relative Humidity (%)
    real(r8), intent(in)  :: e         ! Saturated Vapor Pressure (Pa)
    real(r8), intent(out) :: erh       ! Vapor Pressure (Pa)

!
! !CALLED FROM:
! subroutines within this module
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
    erh = (rh/100._r8) *e   ! Pa

  end subroutine VaporPres
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: QSat_2
!
! !INTERFACE:
  subroutine QSat_2 (T_k, p_t, es_mb, de_mbdT, dlnes_mbdT, rs, rsdT, foftk, fdT)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.  Uses Bolton eqn 10, 39.
! Davies-Jones eqns 2.3,A.1-A.10
! Reference:  Bolton: The computation of equivalent potential temperature. 
!      Monthly weather review (1980) vol. 108 (7) pp. 1046-1053
!      Davies-Jones: An efficient and accurate method for computing the 
!      wet-bulb temperature along pseudoadiabats. Monthly Weather Review 
!      (2008) vol. 136 (7) pp. 2764-2785
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T_k        ! temperature (K)
    real(r8), intent(in)  :: p_t        ! surface atmospheric pressure (pa)
    real(r8), intent(out) :: es_mb      ! vapor pressure (pa)
    real(r8), intent(out) :: de_mbdT    ! d(es)/d(T)
    real(r8), intent(out) :: dlnes_mbdT ! dln(es)/d(T)
    real(r8), intent(out) :: rs         ! humidity (kg/kg)
    real(r8), intent(out) :: rsdT       ! d(qs)/d(T)
    real(r8), intent(out) :: foftk      ! Davies-Jones eqn 2.3
    real(r8), intent(out) :: fdT        ! d(f)/d(T)

!
! !CALLED FROM:
! subroutines within this module
!
! !REVISION HISTORY:
! Created by: Jonathan R Buzan 08/08/13
!
! !LOCAL VARIABLES:
!EOP
!
!
    real(r8) :: lambd_a = 3.504_r8   ! Inverse of Heat Capacity
    real(r8) :: alpha = 17.67_r8     ! Constant to calculate vapour pressure
    real(r8) :: beta = 243.5_r8      ! Constant to calculate vapour pressure
    real(r8) :: epsilon = 0.6220_r8  ! Conversion between pressure/mixing ratio
    real(r8) :: es_C = 6.112_r8      ! Vapor Pressure at Freezing STD (mb)
    real(r8) :: vkp = 0.2854_r8      ! Heat Capacity
    real(r8) :: y0 = 3036._r8        ! constant
    real(r8) :: y1 = 1.78_r8         ! constant
    real(r8) :: y2 = 0.448_r8        ! constant
    real(r8) :: Cf = SHR_CONST_TKFRZ ! Freezing Temp (K)
    real(r8) :: refpres = 1000._r8   ! Reference Pressure (mb)
    real(r8) :: p_tmb                ! Pressure (mb)
    real(r8) :: ndimpress            ! Non-dimensional Pressure
    real(r8) :: prersdt              ! Place Holder for derivative humidity
    real(r8) :: pminuse              ! Vapor Pressure Difference (mb)
    real(r8) :: tcfbdiff             ! Temp diff ref (C)
    real(r8) :: p0ndplam             ! dimensionless pressure modified by ref pressure
    
    real(r8) :: rsy2rs2              ! Constant function of humidity
    real(r8) :: oty2rs               ! Constant function of humidity
    real(r8) :: y0tky1               ! Constant function of Temp

    real(r8) :: d2e_mbdT2            ! d2(es)/d(T)2
    real(r8) :: d2rsdT2              ! d2(r)/d(T)2
    real(r8) :: goftk                ! g(T) exponential in f(T)
    real(r8) :: gdT                  ! d(g)/d(T)
    real(r8) :: d2gdT2               ! d2(g)/d(T)2

    real(r8) :: d2fdT2               ! d2(f)/d(T)2  (K)
!
!-----------------------------------------------------------------------
    ! Constants used to calculate es(T)
    ! Clausius-Clapeyron
    p_tmb = p_t*0.01_r8
    tcfbdiff = T_k - Cf + beta
    es_mb = es_C*exp(alpha*(T_k - Cf)/(tcfbdiff))
    dlnes_mbdT = alpha*beta/((tcfbdiff)*(tcfbdiff))
    pminuse = p_tmb - es_mb
    de_mbdT = es_mb*dlnes_mbdT
    d2e_mbdT2 = dlnes_mbdT*(de_mbdT - 2*es_mb/(tcfbdiff))

    ! Constants used to calculate rs(T)
    ndimpress = (p_tmb/refpres)**vkp
    p0ndplam = refpres*ndimpress**lambd_a
    rs = epsilon*es_mb/(p0ndplam - es_mb)
    prersdt = epsilon*p_tmb/((pminuse)*(pminuse))
    rsdT = prersdt*de_mbdT
    d2rsdT2 = prersdt*(d2e_mbdT2 -de_mbdT*de_mbdT*(2._r8/(pminuse)))

    ! Constants used to calculate g(T)
    rsy2rs2 = rs + y2*rs*rs
    oty2rs = 1._r8 + 2._r8*y2*rs
    y0tky1 = y0/T_k - y1
    goftk = y0tky1*(rs + y2*rs*rs)
    gdT = - y0*(rsy2rs2)/(T_k*T_k) + (y0tky1)*(oty2rs)*rsdT
    d2gdT2 = 2._r8*y0*rsy2rs2/(T_k*T_k*T_k) - 2._r8*y0*rsy2rs2*(oty2rs)*rsdT + &
    y0tky1*2._r8*y2*rsdT*rsdT + y0tky1*oty2rs*d2rsdT2

    ! Calculations for used to calculate f(T,ndimpress)
    foftk = ((Cf/T_k)**lambd_a)*(1._r8 - es_mb/p0ndplam)**(vkp*lambd_a)* &
           exp(-lambd_a*goftk)
    fdT = -lambd_a*(1._r8/T_k + vkp*de_mbdT/pminuse + gdT)
    d2fdT2 = lambd_a*(1._r8/(T_k*T_k) - vkp*de_mbdT*de_mbdT/(pminuse*pminuse) - &
             vkp*d2e_mbdT2/pminuse - d2gdT2)

  end subroutine QSat_2
!EOP
!-----------------------------------------------------------------------

end module HumanIndexMod

