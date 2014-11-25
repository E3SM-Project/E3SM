module EDBioType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! ED/CLM public data type
  !

  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use clm_varpar     , only : nlevgrnd
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC DATA TYPES:
  type, public :: EDbio_type

     real(r8), pointer :: trimming_patch             (:) 
     real(r8), pointer :: area_plant_patch           (:) 
     real(r8), pointer :: area_trees_patch           (:) 
     real(r8), pointer :: canopy_spread_patch        (:) 
     real(r8), pointer :: GCcanopy_patch             (:)   ! mmol m-2 s-1
     real(r8), pointer :: PFTbiomass_patch           (:,:) ! total biomass of each pft
     real(r8), pointer :: PFTleafbiomass_patch       (:,:) ! total biomass of each pft   
     real(r8), pointer :: PFTstorebiomass_patch      (:,:) ! total biomass of each pft   
     real(r8), pointer :: PFTnindivs_patch           (:,:) ! total biomass of each pft 

     real(r8), pointer :: nesterov_fire_danger_patch (:)   ! total biomass of each pft 
     real(r8), pointer :: spitfire_ROS_patch         (:)   ! total biomass of each pft 
     real(r8), pointer :: effect_wspeed_patch        (:)   ! total biomass of each pft 
     real(r8), pointer :: TFC_ROS_patch              (:)   ! total biomass of each pft 
     real(r8), pointer :: fire_intensity_patch       (:)   ! total biomass of each pft 
     real(r8), pointer :: fire_area_patch            (:)   ! total biomass of each pft 
     real(r8), pointer :: scorch_height_patch        (:)   ! total biomass of each pft 
     real(r8), pointer :: fire_fuel_bulkd_patch      (:)   ! total biomass of each pft 
     real(r8), pointer :: fire_fuel_eff_moist_patch  (:)   ! total biomass of each pft 
     real(r8), pointer :: fire_fuel_sav_patch        (:)   ! total biomass of each pft       
     real(r8), pointer :: fire_fuel_mef_patch        (:)   ! total biomass of each pft 
     real(r8), pointer :: sum_fuel_patch             (:)   ! total biomass of each pft 

     real(r8), pointer :: litter_in_patch            (:)   ! total biomass of each pft 
     real(r8), pointer :: litter_out_patch           (:)   ! total biomass of each pft    
     real(r8), pointer :: efpot_patch                (:)   ! potential transpiration
     real(r8), pointer :: rb_patch                   (:)   ! boundary layer conductance

     real(r8), pointer :: daily_temp_patch           (:)   ! daily temperature for fire and phenology models
     real(r8), pointer :: daily_rh_patch             (:)   ! daily RH for fire model
     real(r8), pointer :: daily_prec_patch           (:)   ! daily rain for fire and phenology models. 

     !seed model. Aggregated to gridcell for now. 

     real(r8), pointer :: seed_bank_patch            (:)   ! Mass of seeds.                 kGC/m2
     real(r8), pointer :: seeds_in_patch             (:)   ! Production of seed mass.       kGC/m2/year
     real(r8), pointer :: seed_decay_patch           (:)   ! Decay of seed mass.            kGC/m2/year 
     real(r8), pointer :: seed_germination_patch     (:)   ! Germiantion rate of seed mass. kGC/m2/year

     real(r8), pointer :: ED_bstore_patch            (:)   ! Total stored biomass. kGC/m2
     real(r8), pointer :: ED_bdead_patch             (:)   ! Total dead biomass.   kGC/m2
     real(r8), pointer :: ED_balive_patch            (:)   ! Total alive biomass.  kGC/m2
     real(r8), pointer :: ED_bleaf_patch             (:)   ! Total leaf biomass.   kGC/m2
     real(r8), pointer :: ED_biomass_patch           (:)   ! Total biomass.        kGC/m2

     real(r8), pointer :: ED_GDD_patch               (:)   ! ED Phenology growing degree days. 
                                                           !      This could and should be site-level. RF 
     integer , pointer :: phen_cd_status_patch       (:)   ! ED Phenology cold deciduous status

     ! These variables are also in CNCarbonflux_vars and CNcarbonstate_vars and CNnitrogenstate_vars
     ! They have unique named definitions here - so this assumes currently that both use_ed and use_cn
     ! cannot both be true

   contains

     procedure , public  :: Init   
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory
     procedure , public  :: SetValues

  end type EDbio_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(EDbio_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (EDbio_type) :: this 
    type(bounds_type), intent(in)    :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%trimming_patch             (begp:endp))            ; this%trimming_patch             (:)   = 0.0_r8    
    allocate(this%canopy_spread_patch        (begp:endp))            ; this%canopy_spread_patch        (:)   = 0.0_r8    
    allocate(this%GCcanopy_patch             (begp:endp))            ; this%GCcanopy_patch             (:)   = 0.0_r8    
    allocate(this%area_plant_patch           (begp:endp))            ; this%area_plant_patch           (:)   = 0.0_r8    
    allocate(this%area_trees_patch           (begp:endp))            ; this%area_trees_patch           (:)   = 0.0_r8    
    allocate(this%PFTbiomass_patch           (begp:endp,1:nlevgrnd)) ; this%PFTbiomass_patch           (:,:) = 0.0_r8    
    allocate(this%PFTleafbiomass_patch       (begp:endp,1:nlevgrnd)) ; this%PFTleafbiomass_patch       (:,:) = 0.0_r8    
    allocate(this%PFTstorebiomass_patch      (begp:endp,1:nlevgrnd)) ; this%PFTstorebiomass_patch      (:,:) = 0.0_r8    
    allocate(this%PFTnindivs_patch           (begp:endp,1:nlevgrnd)) ; this%PFTnindivs_patch           (:,:) = 0.0_r8    
    allocate(this%nesterov_fire_danger_patch (begp:endp))            ; this%nesterov_fire_danger_patch (:)   = 0.0_r8    
    allocate(this%spitfire_ROS_patch         (begp:endp))            ; this%spitfire_ROS_patch         (:)   = 0.0_r8    
    allocate(this%effect_wspeed_patch        (begp:endp))            ; this%effect_wspeed_patch        (:)   = 0.0_r8    
    allocate(this%TFC_ROS_patch              (begp:endp))            ; this%TFC_ROS_patch              (:)   = 0.0_r8    
    allocate(this%fire_intensity_patch       (begp:endp))            ; this%fire_intensity_patch       (:)   = 0.0_r8    
    allocate(this%fire_area_patch            (begp:endp))            ; this%fire_area_patch            (:)   = 0.0_r8    
    allocate(this%scorch_height_patch        (begp:endp))            ; this%scorch_height_patch        (:)   = 0.0_r8    
    allocate(this%fire_fuel_bulkd_patch      (begp:endp))            ; this%fire_fuel_bulkd_patch      (:)   = 0.0_r8    
    allocate(this%fire_fuel_eff_moist_patch  (begp:endp))            ; this%fire_fuel_eff_moist_patch  (:)   = 0.0_r8    
    allocate(this%fire_fuel_sav_patch        (begp:endp))            ; this%fire_fuel_sav_patch        (:)   = 0.0_r8    
    allocate(this%fire_fuel_mef_patch        (begp:endp))            ; this%fire_fuel_mef_patch        (:)   = 0.0_r8    
    allocate(this%sum_fuel_patch             (begp:endp))            ; this%sum_fuel_patch             (:)   = 0.0_r8    
    allocate(this%litter_in_patch            (begp:endp))            ; this%litter_in_patch            (:)   = 0.0_r8    
    allocate(this%litter_out_patch           (begp:endp))            ; this%litter_out_patch           (:)   = 0.0_r8    
    allocate(this%efpot_patch                (begp:endp))            ; this%efpot_patch                (:)   = 0.0_r8    
    allocate(this%rb_patch                   (begp:endp))            ; this%rb_patch                   (:)   = 0.0_r8    
    allocate(this%seed_bank_patch            (begp:endp))            ; this%seed_bank_patch            (:)   = 0.0_r8    
    allocate(this%seed_decay_patch           (begp:endp))            ; this%seed_decay_patch           (:)   = 0.0_r8    
    allocate(this%seeds_in_patch             (begp:endp))            ; this%seeds_in_patch             (:)   = 0.0_r8    
    allocate(this%seed_germination_patch     (begp:endp))            ; this%seed_germination_patch     (:)   = 0.0_r8    
    allocate(this%ED_bstore_patch            (begp:endp))            ; this%ED_bstore_patch            (:)   = 0.0_r8    
    allocate(this%ED_bdead_patch             (begp:endp))            ; this%ED_bdead_patch             (:)   = 0.0_r8    
    allocate(this%ED_balive_patch            (begp:endp))            ; this%ED_balive_patch            (:)   = 0.0_r8    
    allocate(this%ED_bleaf_patch             (begp:endp))            ; this%ED_bleaf_patch             (:)   = 0.0_r8    
    allocate(this%ED_biomass_patch           (begp:endp))            ; this%ED_biomass_patch           (:)   = 0.0_r8    
    allocate(this%ED_GDD_patch               (begp:endp))            ; this%ED_GDD_patch               (:)   = 0.0_r8    
    allocate(this%phen_cd_status_patch       (begp:endp))            ; this%phen_cd_status_patch       (:)   = 0.0_r8    

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use clm_varpar , only : nlevdecomp, nlevdecomp_full, crop_prog, nlevgrnd
    use clm_varcon , only : spval
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class(EDbio_type) :: this    
    type(bounds_type)         , intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(8)      :: vr_suffix
    character(10)     :: active
    integer           :: begp,endp
    integer           :: begc,endc
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    call hist_addfld1d (fname='TRIMMING', units='none',  &
         avgflag='A', long_name='Degree to which canopy expansion is limited by leaf economics', &
         ptr_patch=this%trimming_patch, set_lake=0._r8, set_urb=0._r8)  

    call hist_addfld1d (fname='AREA_PLANT', units='m2',  &
         avgflag='A', long_name='area occupied by all plants', &
         ptr_patch=this%area_plant_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='AREA_TREES', units='m2',  &
         avgflag='A', long_name='area occupied by woody plants', &
         ptr_patch=this%area_trees_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='CANOPY_SPREAD', units='none',  &   
         avgflag='A', long_name='Scaling factor between tree basal area and canopy area', &
         ptr_patch=this%canopy_spread_patch, set_lake=0._r8, set_urb=0._r8)   

    call hist_addfld1d (fname='GCCANOPY', units='none',  &
         avgflag='A', long_name='Canopy Conductance: mmol m-2 s-1', &
         ptr_patch=this%GCcanopy_patch, set_lake=0._r8, set_urb=0._r8)  

    call hist_addfld2d (fname='PFTbiomass',  units='kgC/m2', type2d='levgrnd', &
         avgflag='A', long_name='total PFT level biomass', &
         ptr_patch=this%PFTbiomass_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld2d (fname='PFTleafbiomass',  units='kgC/m2', type2d='levgrnd', &
         avgflag='A', long_name='total PFT level biomass', &
         ptr_patch=this%PFTleafbiomass_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld2d (fname='PFTstorebiomass',  units='kgC/m2', type2d='levgrnd', &
         avgflag='A', long_name='total PFT level biomass', &
         ptr_patch=this%PFTstorebiomass_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld2d (fname='PFTnindivs',  units='kgC/m2', type2d='levgrnd', &
         avgflag='A', long_name='total PFT level biomass', &
         ptr_patch=this%PFTnindivs_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FIRE_NESTEROV_INDEX', units='none',  &
         avgflag='A', long_name='nesterov_fire_danger index', &
         ptr_patch=this%nesterov_fire_danger_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FIRE_ROS', units='m/min',  &
         avgflag='A', long_name='fire rate of spread m/min', &
         ptr_patch=this%spitfire_ROS_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='EFFECT_WSPEED', units='none',  &
         avgflag='A', long_name='effective windspeed for fire spread', &
         ptr_patch=this%effect_wspeed_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FIRE_TFC_ROS', units='none',  &
         avgflag='A', long_name='total fuel consumed', &
         ptr_patch=this%TFC_ROS_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FIRE_INTENSITY', units='kJ/m/s',  &
         avgflag='A', long_name='spitfire fire intensity: kJ/m/s', &
         ptr_patch=this%fire_intensity_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FIRE_AREA', units='fraction',  &
         avgflag='A', long_name='spitfire fire area:m2', &
         ptr_patch=this%fire_area_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='SCORCH_HEIGHT', units='m',  &
         avgflag='A', long_name='spitfire fire area:m2', &
         ptr_patch=this%scorch_height_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='fire_fuel_mef', units='m',  &
         avgflag='A', long_name='spitfire fuel moisture', &
         ptr_patch=this%fire_fuel_mef_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='fire_fuel_bulkd', units='m',  &
         avgflag='A', long_name='spitfire fuel bulk density', &
         ptr_patch=this%fire_fuel_bulkd_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='fire_fuel_eff_moist', units='m',  &
         avgflag='A', long_name='spitfire fuel moisture', &
         ptr_patch=this%fire_fuel_eff_moist_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='fire_fuel_sav', units='m',  &
         avgflag='A', long_name='spitfire fuel surface/volume ', &
         ptr_patch=this%fire_fuel_sav_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='TFC_ROS', units='m',  &
         avgflag='A', long_name='spitfire fuel surface/volume ', &
         ptr_patch=this%TFC_ROS_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='SUM_FUEL', units=' KgC m-2 y-1',  &
         avgflag='A', long_name='Litter flux in leaves', &
         ptr_patch=this%sum_fuel_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='LITTER_IN', units=' KgC m-2 y-1',  &
         avgflag='A', long_name='Litter flux in leaves', &
         ptr_patch=this%litter_in_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='LITTER_OUT', units=' KgC m-2 y-1',  &
         avgflag='A', long_name='Litter flux out leaves', &
         ptr_patch=this%litter_out_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='SEED_BANK', units=' KgC m-2',  &
         avgflag='A', long_name='Total Seed Mass of all PFTs', &
         ptr_patch=this%seed_bank_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='SEEDS_IN', units=' KgC m-2 y-1',  &
         avgflag='A', long_name='Seed Production Rate', &
         ptr_patch=this%seeds_in_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='SEED_GERMINATION', units=' KgC m-2 y-1',  &
         avgflag='A', long_name='Seed mass converted into new cohorts', &
         ptr_patch=this%seed_germination_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='SEED_DECAY', units=' KgC m-2 y-1',  &
         avgflag='A', long_name='Seed mass decay', &
         ptr_patch=this%seed_decay_patch, set_lake=0._r8, set_urb=0._r8)              

    call hist_addfld1d (fname='ED_bstore', units=' KgC m-2',  &
         avgflag='A', long_name='ED stored biomass', &
         ptr_patch=this%ED_bstore_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='ED_bdead', units=' KgC m-2',  &
         avgflag='A', long_name='ED dead biomass', &
         ptr_patch=this%ED_bdead_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='ED_balive', units=' KgC m-2',  &
         avgflag='A', long_name='ED live biomass', &
         ptr_patch=this%ED_balive_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='ED_bleaf', units=' KgC m-2',  &
         avgflag='A', long_name='ED leaf biomass', &
         ptr_patch=this%ED_bleaf_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='ED_biomass', units=' KgC m-2',  &
         avgflag='A', long_name='ED total biomass', &
         ptr_patch=this%ED_biomass_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='ED_GDD', units='deg C',  &
         avgflag='A', long_name='ED phenology growing degree days', &
         ptr_patch=this%ED_GDD_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='RB', units=' s m-1',  &
         avgflag='A', long_name='leaf boundary resistance', &
         ptr_patch=this%rb_patch, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='EFPOT', units='',  &
         avgflag='A', long_name='potential evap', &
         ptr_patch=this%efpot_patch, set_lake=0._r8, set_urb=0._r8)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine SetValues( this, val)
    !
    ! !ARGUMENTS:
    class(EDbio_type) :: this

    real(r8), intent(in) :: val

    !
    ! FIX(SPM,082714) - commenting these lines out while merging ED branch to CLM
    ! trunk.  Commented out by RF to work out science issues
    !
    !this%trimming_patch        (:)   = val
    !this%canopy_spread_patch   (:)   = val
    !this%PFTbiomass_patch      (:,:) = val
    !this%PFTleafbiomass_patch  (:,:) = val
    !this%PFTstorebiomass_patch (:,:) = val
    !this%PFTnindivs_patch      (:,:) = val
    this%efpot_patch           (:)   = val
    this%rb_patch              (:)   = val

  end subroutine SetValues

end module EDBioType
