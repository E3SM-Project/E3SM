module dynHarvestMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the harvest data, as well as the state updates that happen as a
  ! result of harvest.
  
  ! Currently, it is assumed that the harvest data are on the flanduse_timeseries file. However, this
  ! could theoretically be changed so that the harvest data were separated from the
  ! pftdyn data, allowing them to differ in the years over which they apply.
  
  ! !USES:
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use abortutils            , only : endrun
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use CNStateType           , only : cnstate_type
  use VegetationPropertiesType        , only : veg_vp
  use elm_varcon            , only : grlnd
  use ColumnType            , only : col_pp
  use ColumnDataType        , only : col_cf, col_nf, col_pf  
  use VegetationType        , only : veg_pp                
  use VegetationDataType    , only : veg_cs, veg_cf, veg_ns, veg_nf  
  use topounit_varcon       , only : max_topounits
  use VegetationDataType    , only : veg_ps, veg_pf
  use elm_varctl            , only : use_cn, use_fates, iulog
  use FatesConstantsMod      , only : hlm_harvest_area_fraction
  use FatesConstantsMod      , only : hlm_harvest_carbon

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynHarvest_init    ! initialize data structures for harvest information
  public :: dynHarvest_interp_harvest_types  ! get harvest data for each harvest type, if needed
  public :: CNHarvest          ! harvest mortality routine for CN code (non-FATES)
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNHarvestPftToColumn   ! gather pft-level harvest fluxes to the column level
  
  ! !PRIVATE TYPES:

  ! Note that, since we have our own dynHarvest_file object (distinct from dynpft_file),
  ! we could theoretically have a separate file providing harvest data from that providing
  ! the pftdyn data
  type(dyn_file_type), target :: dynHarvest_file ! information for the file containing harvest data

  ! Define the underlying input harvest variables
  ! these are hardcoded to match the LUH input data
  ! these are fraction of vegetated area harvested, split into five harvest type variables
  ! in next iteration (v3) change these to fraction of grid cell
  ! HARVEST_VH1 = harvest from primary forest
  ! HARVEST_VH2 = harvest from primary non-forest
  ! HARVEST_SH1 = harvest from secondary mature forest
  ! HARVEST_SH2 = harvest from secondary young forest
  ! HARVEST_SH3 = harvest from secondary non-forest
  ! note that these are summed for CNHarvest in harvest(:) and applied to total forest area only
  ! the individual category rates are passed to fates in harvest_rates(:,:)

  ! for FATES: capacity for passing harvest data in units of carbon harvested per year (per grid cell) has been added
  ! but these data are not yet included in the input file
  ! the code here can be changed to wood_harvest_units = harvest_carbon to pass carbon data to FATES if:
  ! the carbon data are in the same input variables as listed below
  ! and the carbon units in the input file match that expected by FATES

  integer, public, parameter :: num_harvest_vars = 5
  character(len=64), public, parameter :: harvest_varnames(num_harvest_vars) = &
       [character(len=64) :: 'HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 'HARVEST_SH2', 'HARVEST_SH3']

  type(dyn_var_time_uninterp_type) :: harvest_vars(num_harvest_vars)   ! value of each harvest variable

  ! the units flag must match the units of harvest_varnames
  ! set this here because dynHarvest_init is called after alm_fates%init
  ! this flag is accessed only if namelist do_harvest is TRUE

  integer, public, parameter    :: wood_harvest_units = 2    ! 1 = area fraction, 2 = carbon
  real(r8), allocatable, public :: harvest_rates(:,:) ! harvest rates
  logical, private              :: do_harvest ! whether we're in a period when we should do harvest
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynHarvest_init(bounds, harvest_filename)
    
    ! !DESCRIPTION:
    ! Initialize data structures for harvest information.
    ! This should be called once, during model initialization.
     
    ! !USES:
    use elm_varctl            , only : use_cn, use_fates
    use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
    use dynTimeInfoMod        , only : YEAR_POSITION_START_OF_TIMESTEP
    use dynTimeInfoMod        , only : YEAR_POSITION_END_OF_TIMESTEP
    
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds           ! proc-level bounds
    character(len=*) , intent(in) :: harvest_filename ! name of file containing harvest information
    
    ! !LOCAL VARIABLES:
    integer :: varnum     ! counter for harvest variables
    integer :: harvest_shape(1)  ! harvest shape 
    integer :: num_points ! number of spatial points
    integer :: ier        ! error code
    
    character(len=*), parameter :: subname = 'dynHarvest_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    allocate(harvest_rates(num_harvest_vars,bounds%begg:bounds%endg),stat=ier)
    harvest_rates(:,bounds%begg:bounds%endg) = 0._r8
    if (ier /= 0) then
       call endrun(msg=' allocation error for harvest_rates'//errMsg(__FILE__, __LINE__))
    end if

    !dynHarvest_file = dyn_file_type(harvest_filename, YEAR_POSITION_START_OF_TIMESTEP)
    dynHarvest_file = dyn_file_type(harvest_filename, YEAR_POSITION_END_OF_TIMESTEP)
    
    ! Get initial harvest data
    if (use_cn .or. use_fates) then
       num_points = (bounds%endg - bounds%begg + 1)
       harvest_shape(1) = num_points
       do varnum = 1, num_harvest_vars
          harvest_vars(varnum) = dyn_var_time_uninterp_type( &
               dyn_file=dynHarvest_file, varname=harvest_varnames(varnum), &
               dim1name=grlnd, conversion_factor=1.0_r8, &
               do_check_sums_equal_1=.false., data_shape=harvest_shape)
       end do
    end if
    
  end subroutine dynHarvest_init

!-----------------------------------------------------------------------
  subroutine dynHarvest_interp_harvest_types(bounds)
    !
    ! !DESCRIPTION:
    ! Get harvest data for model time, by land category, when needed.
    ! this function called only for CN (non-FATES) or FATES
    !
    ! Note that harvest data are stored as rates (not weights) and so time interpolation
    ! is not necessary - the harvest rate is held constant through the year.  This is
    ! consistent with the treatment of changing PFT weights, where interpolation of the
    ! annual endpoint weights leads to a constant rate of change in PFT weight through the
    ! year, with abrupt changes in the rate at annual boundaries.
    !
    ! Note the difference between this and the old dynHarvest_interp is that here, we keep the different
    ! forcing sets distinct (e.g., for passing to FATES which has distinct primary and secondary lands)
    ! and thus store it in harvest_rates
    !
    ! !USES:
    use dynTimeInfoMod , only : time_info_type
    use elm_varctl     , only : use_cn, use_fates
    
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    
    ! !LOCAL VARIABLES:
    integer               :: varnum       ! counter for harvest variables
    real(r8), allocatable :: this_data(:) ! data for a single harvest variable
    character(len=*), parameter :: subname = 'dynHarvest_interp_harvest_types'
    !-----------------------------------------------------------------------
    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! input harvest data for current year are stored in year+1 in the file
    call dynHarvest_file%time_info%set_current_year_get_year(1)
    if (use_cn .or. use_fates) then
       harvest_rates(1:num_harvest_vars,bounds%begg:bounds%endg) = 0._r8

       if (dynHarvest_file%time_info%is_before_time_series()) then
          ! Turn off harvest before the start of the harvest time series
          do_harvest = .false.
       else
          ! Note that do_harvest stays true even past the end of the time series. This
          ! means that harvest rates will be maintained at the rate given in the last
          ! year of the file for all years past the end of this specified time series.
          do_harvest = .true.
          ! Right now we don't account for the topounit in plant harvest
          allocate(this_data(bounds%begg:bounds%endg))
          do varnum = 1, num_harvest_vars
             call harvest_vars(varnum)%get_current_data(this_data)
             harvest_rates(varnum,bounds%begg:bounds%endg) = this_data(bounds%begg:bounds%endg)
          end do
          deallocate(this_data)
       end if
    end if
  end subroutine dynHarvest_interp_harvest_types


  !-----------------------------------------------------------------------
  subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp, cnstate_vars)
    !
    ! !DESCRIPTION:
    ! Harvest mortality routine for coupled carbon-nitrogen code (CN)
    
    ! !USES:
    use pftvarcon       , only : noveg, nbrdlf_evr_shrub, pprodharv10
    use elm_varcon      , only : secspday
    use clm_time_manager, only : get_days_per_year
    use GridcellType   , only : grc_pp
    
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    ! !LOCAL VARIABLES:
    integer :: p                         ! patch index
    integer :: t,ti,topi                 ! topounit indices TKT
    integer :: g                         ! gridcell index
    integer :: fp                        ! patch filter index
    real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
    real(r8):: m                         ! rate for fractional harvest mortality (1/s)
    real(r8):: days_per_year             ! days per year
    integer :: varnum                    ! counter for harvest variables
    !-----------------------------------------------------------------------

   associate(& 
   pgridcell                           =>    veg_pp%gridcell                                , & ! Input:  [integer (:)]  pft-level index into gridcell-level quantities     
   ivt                                 =>    veg_pp%itype                                   , & ! Input:  [integer (:)]  pft vegetation type                                

   leafc                               =>    veg_cs%leafc                                   , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   frootc                              =>    veg_cs%frootc                                  , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   livestemc                           =>    veg_cs%livestemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   deadstemc                           =>    veg_cs%deadstemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   livecrootc                          =>    veg_cs%livecrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   deadcrootc                          =>    veg_cs%deadcrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   xsmrpool                            =>    veg_cs%xsmrpool                                , & ! Input:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
   leafc_storage                       =>    veg_cs%leafc_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   frootc_storage                      =>    veg_cs%frootc_storage                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   livestemc_storage                   =>    veg_cs%livestemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   deadstemc_storage                   =>    veg_cs%deadstemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   livecrootc_storage                  =>    veg_cs%livecrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   deadcrootc_storage                  =>    veg_cs%deadcrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   gresp_storage                       =>    veg_cs%gresp_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   leafc_xfer                          =>    veg_cs%leafc_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   frootc_xfer                         =>    veg_cs%frootc_xfer                             , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livestemc_xfer                      =>    veg_cs%livestemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc_xfer                      =>    veg_cs%deadstemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   livecrootc_xfer                     =>    veg_cs%livecrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc_xfer                     =>    veg_cs%deadcrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   gresp_xfer                          =>    veg_cs%gresp_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   cpool                               =>    veg_cs%cpool                                   , & ! Input:  [real(r8) (:)]  (gC/m2) Plant C storage

   leafn                               =>    veg_ns%leafn                                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   frootn                              =>    veg_ns%frootn                                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   livestemn                           =>    veg_ns%livestemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   deadstemn                           =>    veg_ns%deadstemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
   livecrootn                          =>    veg_ns%livecrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   deadcrootn                          =>    veg_ns%deadcrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   retransn                            =>    veg_ns%retransn                                , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   npool                               =>    veg_ns%npool                                   , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of stored N 
   leafn_storage                       =>    veg_ns%leafn_storage                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   frootn_storage                      =>    veg_ns%frootn_storage                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   livestemn_storage                   =>    veg_ns%livestemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   deadstemn_storage                   =>    veg_ns%deadstemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   livecrootn_storage                  =>    veg_ns%livecrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   deadcrootn_storage                  =>    veg_ns%deadcrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   leafn_xfer                          =>    veg_ns%leafn_xfer                              , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   frootn_xfer                         =>    veg_ns%frootn_xfer                             , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   livestemn_xfer                      =>    veg_ns%livestemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   deadstemn_xfer                      =>    veg_ns%deadstemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   livecrootn_xfer                     =>    veg_ns%livecrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   deadcrootn_xfer                     =>    veg_ns%deadcrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   
   ! add phosphorus
   leafp                               =>    veg_ps%leafp                                   , & ! Input:  [real(r8) (:)]  (gP/m2) leaf P                                    
   frootp                              =>    veg_ps%frootp                                  , & ! Input:  [real(r8) (:)]  (gP/m2) fine root P                               
   livestemp                           =>    veg_ps%livestemp                               , & ! Input:  [real(r8) (:)]  (gP/m2) live stem P                               
   deadstemp                           =>    veg_ps%deadstemp                               , & ! Input:  [real(r8) (:)]  (gP/m2) dead stem P                               
   livecrootp                          =>    veg_ps%livecrootp                              , & ! Input:  [real(r8) (:)]  (gP/m2) live coarse root P                        
   deadcrootp                          =>    veg_ps%deadcrootp                              , & ! Input:  [real(r8) (:)]  (gP/m2) dead coarse root P                        
   retransp                            =>    veg_ps%retransp                                , & ! Input:  [real(r8) (:)]  (gP/m2) plant pool of retranslocated P            
   ppool                               =>    veg_ps%ppool                                   , & ! Input:  [real(r8) (:)]  (gP/m2) plant pool of storage P  
   leafp_storage                       =>    veg_ps%leafp_storage                           , & ! Input:  [real(r8) (:)]  (gP/m2) leaf P storage                            
   frootp_storage                      =>    veg_ps%frootp_storage                          , & ! Input:  [real(r8) (:)]  (gP/m2) fine root P storage                       
   livestemp_storage                   =>    veg_ps%livestemp_storage                       , & ! Input:  [real(r8) (:)]  (gP/m2) live stem P storage                       
   deadstemp_storage                   =>    veg_ps%deadstemp_storage                       , & ! Input:  [real(r8) (:)]  (gP/m2) dead stem P storage                       
   livecrootp_storage                  =>    veg_ps%livecrootp_storage                      , & ! Input:  [real(r8) (:)]  (gP/m2) live coarse root P storage                
   deadcrootp_storage                  =>    veg_ps%deadcrootp_storage                      , & ! Input:  [real(r8) (:)]  (gP/m2) dead coarse root P storage                
   leafp_xfer                          =>    veg_ps%leafp_xfer                              , & ! Input:  [real(r8) (:)]  (gP/m2) leaf P transfer                           
   frootp_xfer                         =>    veg_ps%frootp_xfer                             , & ! Input:  [real(r8) (:)]  (gP/m2) fine root P transfer                      
   livestemp_xfer                      =>    veg_ps%livestemp_xfer                          , & ! Input:  [real(r8) (:)]  (gP/m2) live stem P transfer                      
   deadstemp_xfer                      =>    veg_ps%deadstemp_xfer                          , & ! Input:  [real(r8) (:)]  (gP/m2) dead stem P transfer                      
   livecrootp_xfer                     =>    veg_ps%livecrootp_xfer                         , & ! Input:  [real(r8) (:)]  (gP/m2) live coarse root P transfer               
   deadcrootp_xfer                     =>    veg_ps%deadcrootp_xfer                         , & ! Input:  [real(r8) (:)]  (gP/m2) dead coarse root P transfer               

   hrv_leafc_to_litter                 =>    veg_cf%hrv_leafc_to_litter                     , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_to_litter                =>    veg_cf%hrv_frootc_to_litter                    , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_to_litter             =>    veg_cf%hrv_livestemc_to_litter                 , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod10c            =>    veg_cf%hrv_deadstemc_to_prod10c                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod100c           =>    veg_cf%hrv_deadstemc_to_prod100c               , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_to_litter            =>    veg_cf%hrv_livecrootc_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_to_litter            =>    veg_cf%hrv_deadcrootc_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_xsmrpool_to_atm                 =>    veg_cf%hrv_xsmrpool_to_atm                     , & ! Output: [real(r8) (:)]                                                    
   hrv_leafc_storage_to_litter         =>    veg_cf%hrv_leafc_storage_to_litter             , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_storage_to_litter        =>    veg_cf%hrv_frootc_storage_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_storage_to_litter     =>    veg_cf%hrv_livestemc_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_storage_to_litter     =>    veg_cf%hrv_deadstemc_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_storage_to_litter    =>    veg_cf%hrv_livecrootc_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_storage_to_litter    =>    veg_cf%hrv_deadcrootc_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_gresp_storage_to_litter         =>    veg_cf%hrv_gresp_storage_to_litter             , & ! Output: [real(r8) (:)]                                                    
   hrv_leafc_xfer_to_litter            =>    veg_cf%hrv_leafc_xfer_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_xfer_to_litter           =>    veg_cf%hrv_frootc_xfer_to_litter               , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_xfer_to_litter        =>    veg_cf%hrv_livestemc_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_xfer_to_litter        =>    veg_cf%hrv_deadstemc_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_xfer_to_litter       =>    veg_cf%hrv_livecrootc_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_xfer_to_litter       =>    veg_cf%hrv_deadcrootc_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    
   hrv_gresp_xfer_to_litter            =>    veg_cf%hrv_gresp_xfer_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_cpool_to_litter                 =>    veg_cf%hrv_cpool_to_litter                     , & ! Output: [real(r8) (:)]             

   hrv_leafn_to_litter                 =>    veg_nf%hrv_leafn_to_litter                     , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_to_litter                =>    veg_nf%hrv_frootn_to_litter                    , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_to_litter             =>    veg_nf%hrv_livestemn_to_litter                 , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod10n            =>    veg_nf%hrv_deadstemn_to_prod10n                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod100n           =>    veg_nf%hrv_deadstemn_to_prod100n               , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_to_litter            =>    veg_nf%hrv_livecrootn_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_to_litter            =>    veg_nf%hrv_deadcrootn_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_retransn_to_litter              =>    veg_nf%hrv_retransn_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   hrv_npool_to_litter                 =>    veg_nf%hrv_npool_to_litter                     , & ! Output: [real(r8) (:)]   
   hrv_leafn_storage_to_litter         =>    veg_nf%hrv_leafn_storage_to_litter             , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_storage_to_litter        =>    veg_nf%hrv_frootn_storage_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_storage_to_litter     =>    veg_nf%hrv_livestemn_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_storage_to_litter     =>    veg_nf%hrv_deadstemn_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_storage_to_litter    =>    veg_nf%hrv_livecrootn_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_storage_to_litter    =>    veg_nf%hrv_deadcrootn_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_leafn_xfer_to_litter            =>    veg_nf%hrv_leafn_xfer_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_xfer_to_litter           =>    veg_nf%hrv_frootn_xfer_to_litter               , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_xfer_to_litter        =>    veg_nf%hrv_livestemn_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_xfer_to_litter        =>    veg_nf%hrv_deadstemn_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_xfer_to_litter       =>    veg_nf%hrv_livecrootn_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_xfer_to_litter       =>    veg_nf%hrv_deadcrootn_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    

   hrv_leafp_to_litter                 =>    veg_pf%hrv_leafp_to_litter                     , & ! Output: [real(r8) (:)]                                                    
   hrv_frootp_to_litter                =>    veg_pf%hrv_frootp_to_litter                    , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemp_to_litter             =>    veg_pf%hrv_livestemp_to_litter                 , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemp_to_prod10p            =>    veg_pf%hrv_deadstemp_to_prod10p                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemp_to_prod100p           =>    veg_pf%hrv_deadstemp_to_prod100p               , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootp_to_litter            =>    veg_pf%hrv_livecrootp_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootp_to_litter            =>    veg_pf%hrv_deadcrootp_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_retransp_to_litter              =>    veg_pf%hrv_retransp_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   hrv_ppool_to_litter                 =>    veg_pf%hrv_ppool_to_litter                     , & ! Output: [real(r8) (:)] 
   hrv_leafp_storage_to_litter         =>    veg_pf%hrv_leafp_storage_to_litter             , & ! Output: [real(r8) (:)]                                                    
   hrv_frootp_storage_to_litter        =>    veg_pf%hrv_frootp_storage_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemp_storage_to_litter     =>    veg_pf%hrv_livestemp_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemp_storage_to_litter     =>    veg_pf%hrv_deadstemp_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootp_storage_to_litter    =>    veg_pf%hrv_livecrootp_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootp_storage_to_litter    =>    veg_pf%hrv_deadcrootp_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_leafp_xfer_to_litter            =>    veg_pf%hrv_leafp_xfer_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_frootp_xfer_to_litter           =>    veg_pf%hrv_frootp_xfer_to_litter               , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemp_xfer_to_litter        =>    veg_pf%hrv_livestemp_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemp_xfer_to_litter        =>    veg_pf%hrv_deadstemp_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootp_xfer_to_litter       =>    veg_pf%hrv_livecrootp_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootp_xfer_to_litter       =>    veg_pf%hrv_deadcrootp_xfer_to_litter             & ! Output: [real(r8) (:)]                                                    
   )


   days_per_year = get_days_per_year()

   ! patch loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      g = pgridcell(p)
      t = veg_pp%topounit(p)
      topi = grc_pp%topi(g)
      ti = t - topi + 1
      
      ! If this is a tree pft, then
      ! get the annual harvest "mortality" rate (am) from harvest array
      ! and convert to rate per second
      if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then

         if (do_harvest) then
            am = 0._r8
            do varnum = 1, num_harvest_vars
               am = am + harvest_rates(varnum,g)
            end do
            m  = am/(days_per_year * secspday)
         else
            m = 0._r8
         end if   

         ! pft-level harvest carbon fluxes
         ! displayed pools
         hrv_leafc_to_litter(p)               = leafc(p)               * m
         hrv_frootc_to_litter(p)              = frootc(p)              * m
         hrv_livestemc_to_litter(p)           = livestemc(p)           * m
         hrv_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
         hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
         hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

         ! storage pools
         hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
         hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
         hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
         hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
         hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
         hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
         hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m
         hrv_cpool_to_litter(p)               = cpool(p)               * m

         ! transfer pools
         hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
         hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
         hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
         hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
         hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
         hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
         hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

         ! pft-level harvest mortality nitrogen fluxes
         ! displayed pools
         hrv_leafn_to_litter(p)               = leafn(p)               * m
         hrv_frootn_to_litter(p)              = frootn(p)              * m
         hrv_livestemn_to_litter(p)           = livestemn(p)           * m
         hrv_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
         hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
         hrv_retransn_to_litter(p)            = retransn(p)            * m
         hrv_npool_to_litter(p)               = npool(p)               * m

         ! storage pools
         hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
         hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
         hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
         hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
         hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
         hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

         ! transfer pools
         hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
         hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
         hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
         hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
         hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
         hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m
         
         ! pft-level harvest mortality phosphorus fluxes
         ! displayed pools
         hrv_leafp_to_litter(p)               = leafp(p)               * m
         hrv_frootp_to_litter(p)              = frootp(p)              * m
         hrv_livestemp_to_litter(p)           = livestemp(p)           * m
         hrv_deadstemp_to_prod10p(p)          = deadstemp(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemp_to_prod100p(p)         = deadstemp(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootp_to_litter(p)          = livecrootp(p)          * m
         hrv_deadcrootp_to_litter(p)          = deadcrootp(p)          * m
         hrv_retransp_to_litter(p)            = retransp(p)            * m
         hrv_ppool_to_litter(p)               = ppool(p)               * m

         ! storage pools
         hrv_leafp_storage_to_litter(p)       = leafp_storage(p)       * m
         hrv_frootp_storage_to_litter(p)      = frootp_storage(p)      * m
         hrv_livestemp_storage_to_litter(p)   = livestemp_storage(p)   * m
         hrv_deadstemp_storage_to_litter(p)   = deadstemp_storage(p)   * m
         hrv_livecrootp_storage_to_litter(p)  = livecrootp_storage(p)  * m
         hrv_deadcrootp_storage_to_litter(p)  = deadcrootp_storage(p)  * m

         ! transfer pools
         hrv_leafp_xfer_to_litter(p)          = leafp_xfer(p)          * m
         hrv_frootp_xfer_to_litter(p)         = frootp_xfer(p)         * m
         hrv_livestemp_xfer_to_litter(p)      = livestemp_xfer(p)      * m
         hrv_deadstemp_xfer_to_litter(p)      = deadstemp_xfer(p)      * m
         hrv_livecrootp_xfer_to_litter(p)     = livecrootp_xfer(p)     * m
         hrv_deadcrootp_xfer_to_litter(p)     = deadcrootp_xfer(p)     * m

      end if  ! end tree block

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes from harvest to the column
   ! for litter C and N inputs
   ! and litter P inputs

   call CNHarvestPftToColumn(num_soilc, filter_soilc, &
      cnstate_vars)

    end associate 
 end subroutine CNHarvest

 !-----------------------------------------------------------------------
 subroutine CNHarvestPftToColumn (num_soilc, filter_soilc, &
      cnstate_vars)
   !
   ! !DESCRIPTION:
   ! called at the end of CNHarvest to gather all pft-level harvest litterfall fluxes
   ! to the column level and assign them to the three litter pools
   
   ! !USES:
   use elm_varpar , only : maxpatch_pft, nlevdecomp
   !
   ! !ARGUMENTS:
   integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                   , intent(in)    :: filter_soilc(:) ! soil column filter
   type(cnstate_type)        , intent(in)    :: cnstate_vars
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
   !-----------------------------------------------------------------------

   associate(                                                                                        & 
        ivt                              =>    veg_pp%itype                                                , & ! Input:  [integer  (:)   ]  pft vegetation type                                
        wtcol                            =>    veg_pp%wtcol                                                , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
        
        lf_flab                          =>    veg_vp%lf_flab                                       , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
        lf_fcel                          =>    veg_vp%lf_fcel                                       , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
        lf_flig                          =>    veg_vp%lf_flig                                       , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
        fr_flab                          =>    veg_vp%fr_flab                                       , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
        fr_fcel                          =>    veg_vp%fr_fcel                                       , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
        fr_flig                          =>    veg_vp%fr_flig                                       , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  
        
        leaf_prof                        =>    cnstate_vars%leaf_prof_patch                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
        froot_prof                       =>    cnstate_vars%froot_prof_patch                            , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
        croot_prof                       =>    cnstate_vars%croot_prof_patch                            , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
        stem_prof                        =>    cnstate_vars%stem_prof_patch                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
        
        hrv_leafc_to_litter              =>    veg_cf%hrv_leafc_to_litter                , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_to_litter             =>    veg_cf%hrv_frootc_to_litter               , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_to_litter          =>    veg_cf%hrv_livestemc_to_litter            , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemc_to_prod10c        =>    veg_cf%hrv_deadstemc_to_prod10c           , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemc_to_prod100c       =>    veg_cf%hrv_deadstemc_to_prod100c          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_to_litter         =>    veg_cf%hrv_livecrootc_to_litter           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_to_litter         =>    veg_cf%hrv_deadcrootc_to_litter           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_storage_to_litter      =>    veg_cf%hrv_leafc_storage_to_litter        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_storage_to_litter     =>    veg_cf%hrv_frootc_storage_to_litter       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_storage_to_litter  =>    veg_cf%hrv_livestemc_storage_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_storage_to_litter  =>    veg_cf%hrv_deadstemc_storage_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_storage_to_litter =>    veg_cf%hrv_livecrootc_storage_to_litter   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_storage_to_litter =>    veg_cf%hrv_deadcrootc_storage_to_litter   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_storage_to_litter      =>    veg_cf%hrv_gresp_storage_to_litter        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_xfer_to_litter         =>    veg_cf%hrv_leafc_xfer_to_litter           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_xfer_to_litter        =>    veg_cf%hrv_frootc_xfer_to_litter          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_xfer_to_litter     =>    veg_cf%hrv_livestemc_xfer_to_litter       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_xfer_to_litter     =>    veg_cf%hrv_deadstemc_xfer_to_litter       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_xfer_to_litter    =>    veg_cf%hrv_livecrootc_xfer_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_xfer_to_litter    =>    veg_cf%hrv_deadcrootc_xfer_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_xfer_to_litter         =>    veg_cf%hrv_gresp_xfer_to_litter           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_cpool_to_litter              =>    veg_cf%hrv_cpool_to_litter                , & ! Input:  [real(r8) (:)   ]       
        chrv_deadstemc_to_prod10c        =>    col_cf%hrv_deadstemc_to_prod10c             , & ! InOut:  [real(r8) (:)   ]                                                    
        chrv_deadstemc_to_prod100c       =>    col_cf%hrv_deadstemc_to_prod100c            , & ! InOut:  [real(r8) (:)   ]                                                    
        harvest_c_to_litr_met_c          =>    col_cf%harvest_c_to_litr_met_c              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
        harvest_c_to_litr_cel_c          =>    col_cf%harvest_c_to_litr_cel_c              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
        harvest_c_to_litr_lig_c          =>    col_cf%harvest_c_to_litr_lig_c              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
        harvest_c_to_cwdc                =>    col_cf%harvest_c_to_cwdc                    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
        
        hrv_leafn_to_litter              =>    veg_nf%hrv_leafn_to_litter              , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_to_litter             =>    veg_nf%hrv_frootn_to_litter             , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_to_litter          =>    veg_nf%hrv_livestemn_to_litter          , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemn_to_prod10n        =>    veg_nf%hrv_deadstemn_to_prod10n         , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemn_to_prod100n       =>    veg_nf%hrv_deadstemn_to_prod100n        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_to_litter         =>    veg_nf%hrv_livecrootn_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_to_litter         =>    veg_nf%hrv_deadcrootn_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_retransn_to_litter           =>    veg_nf%hrv_retransn_to_litter           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_npool_to_litter              =>    veg_nf%hrv_npool_to_litter              , & ! Input:  [real(r8) (:)   ]
        hrv_leafn_storage_to_litter      =>    veg_nf%hrv_leafn_storage_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_storage_to_litter     =>    veg_nf%hrv_frootn_storage_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_storage_to_litter  =>    veg_nf%hrv_livestemn_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_storage_to_litter  =>    veg_nf%hrv_deadstemn_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_storage_to_litter =>    veg_nf%hrv_livecrootn_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_storage_to_litter =>    veg_nf%hrv_deadcrootn_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_xfer_to_litter         =>    veg_nf%hrv_leafn_xfer_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_xfer_to_litter        =>    veg_nf%hrv_frootn_xfer_to_litter        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_xfer_to_litter     =>    veg_nf%hrv_livestemn_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_xfer_to_litter     =>    veg_nf%hrv_deadstemn_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_xfer_to_litter    =>    veg_nf%hrv_livecrootn_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_xfer_to_litter    =>    veg_nf%hrv_deadcrootn_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
        chrv_deadstemn_to_prod10n        =>    col_nf%hrv_deadstemn_to_prod10n           , & ! InOut:  [real(r8) (:)   ]                                                    
        chrv_deadstemn_to_prod100n       =>    col_nf%hrv_deadstemn_to_prod100n          , & ! InOut:  [real(r8) (:)   ]                                                    
        harvest_n_to_litr_met_n          =>    col_nf%harvest_n_to_litr_met_n            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
        harvest_n_to_litr_cel_n          =>    col_nf%harvest_n_to_litr_cel_n            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
        harvest_n_to_litr_lig_n          =>    col_nf%harvest_n_to_litr_lig_n            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter lignin pool (gN/m3/s)
        harvest_n_to_cwdn                =>    col_nf%harvest_n_to_cwdn                  ,  & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to CWD pool (gN/m3/s)
        
        ! add P harvest fluxes 
        hrv_leafp_to_litter              =>    veg_pf%hrv_leafp_to_litter              , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootp_to_litter             =>    veg_pf%hrv_frootp_to_litter             , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemp_to_litter          =>    veg_pf%hrv_livestemp_to_litter          , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemp_to_prod10p        =>    veg_pf%hrv_deadstemp_to_prod10p         , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemp_to_prod100p       =>    veg_pf%hrv_deadstemp_to_prod100p        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootp_to_litter         =>    veg_pf%hrv_livecrootp_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootp_to_litter         =>    veg_pf%hrv_deadcrootp_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_retransp_to_litter           =>    veg_pf%hrv_retransp_to_litter           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_ppool_to_litter              =>    veg_pf%hrv_ppool_to_litter              , & ! Input:  [real(r8) (:)   ]   
        hrv_leafp_storage_to_litter      =>    veg_pf%hrv_leafp_storage_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootp_storage_to_litter     =>    veg_pf%hrv_frootp_storage_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemp_storage_to_litter  =>    veg_pf%hrv_livestemp_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemp_storage_to_litter  =>    veg_pf%hrv_deadstemp_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootp_storage_to_litter =>    veg_pf%hrv_livecrootp_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootp_storage_to_litter =>    veg_pf%hrv_deadcrootp_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafp_xfer_to_litter         =>    veg_pf%hrv_leafp_xfer_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootp_xfer_to_litter        =>    veg_pf%hrv_frootp_xfer_to_litter        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemp_xfer_to_litter     =>    veg_pf%hrv_livestemp_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemp_xfer_to_litter     =>    veg_pf%hrv_deadstemp_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootp_xfer_to_litter    =>    veg_pf%hrv_livecrootp_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootp_xfer_to_litter    =>    veg_pf%hrv_deadcrootp_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
        chrv_deadstemp_to_prod10p        =>    col_pf%hrv_deadstemp_to_prod10p           , & ! InOut:  [real(r8) (:)   ]                                                    
        chrv_deadstemp_to_prod100p       =>    col_pf%hrv_deadstemp_to_prod100p          , & ! InOut:  [real(r8) (:)   ]                                                    
        harvest_p_to_litr_met_p          =>    col_pf%harvest_p_to_litr_met_p            , & ! InOut:  [real(r8) (:,:) ]  P fluxes associated with harvest to litter metabolic pool (gP/m3/s)
        harvest_p_to_litr_cel_p          =>    col_pf%harvest_p_to_litr_cel_p            , & ! InOut:  [real(r8) (:,:) ]  P fluxes associated with harvest to litter cellulose pool (gP/m3/s)
        harvest_p_to_litr_lig_p          =>    col_pf%harvest_p_to_litr_lig_p            , & ! InOut:  [real(r8) (:,:) ]  P fluxes associated with harvest to litter lignin pool (gP/m3/s)
        harvest_p_to_cwdp                =>    col_pf%harvest_p_to_cwdp                    & ! InOut:  [real(r8) (:,:) ]  P fluxes associated with harvest to CWD pool (gP/m3/s)

        )

     do j = 1, nlevdecomp
        do pi = 1,maxpatch_pft
           do fc = 1,num_soilc
              c = filter_soilc(fc)

              if (pi <=  col_pp%npfts(c)) then
                 p = col_pp%pfti(c) + pi - 1

                 if (veg_pp%active(p)) then

                    ! leaf harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                         hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                         hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                    ! fine root harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                         hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                         hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                    ! wood harvest mortality carbon fluxes
                    harvest_c_to_cwdc(c,j)  = harvest_c_to_cwdc(c,j)  + &
                         hrv_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j) 
                    harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                         hrv_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                         hrv_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 

                    ! storage harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_leafc_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                         hrv_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_livestemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_deadstemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_cpool_to_litter(p)      * wtcol(p) * leaf_prof(p,j)


                    ! transfer harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_leafc_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                         hrv_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_livestemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_deadstemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)

                    ! leaf harvest mortality nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                         hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                         hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                    ! fine root litter nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                         hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                         hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                    ! wood harvest mortality nitrogen fluxes
                    harvest_n_to_cwdn(c,j)  = harvest_n_to_cwdn(c,j)  + &
                         hrv_livestemn_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                         hrv_livecrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                         hrv_deadcrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! retranslocated N pool harvest mortality fluxes
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_npool_to_litter(p) * wtcol(p) * leaf_prof(p,j)

                    ! storage harvest mortality nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                         hrv_leafn_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                         hrv_frootn_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_livestemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_deadstemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_livecrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_deadcrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! transfer harvest mortality nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                         hrv_leafn_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                         hrv_frootn_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_livestemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_deadstemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_livecrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_deadcrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)


                    ! leaf harvest mortality phosphorus fluxes
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_leafp_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_p_to_litr_cel_p(c,j) = harvest_p_to_litr_cel_p(c,j) + &
                         hrv_leafp_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_p_to_litr_lig_p(c,j) = harvest_p_to_litr_lig_p(c,j) + &
                         hrv_leafp_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                    ! fine root litter phosphorus fluxes
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_frootp_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_p_to_litr_cel_p(c,j) = harvest_p_to_litr_cel_p(c,j) + &
                         hrv_frootp_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_p_to_litr_lig_p(c,j) = harvest_p_to_litr_lig_p(c,j) + &
                         hrv_frootp_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                    ! wood harvest mortality phosphorus fluxes
                    harvest_p_to_cwdp(c,j)  = harvest_p_to_cwdp(c,j)  + &
                         hrv_livestemp_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_p_to_cwdp(c,j) = harvest_p_to_cwdp(c,j) + &
                         hrv_livecrootp_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_p_to_cwdp(c,j) = harvest_p_to_cwdp(c,j) + &
                         hrv_deadcrootp_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! retranslocated phosphorus pool harvest mortality fluxes
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_retransp_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_ppool_to_litter(p) * wtcol(p) * leaf_prof(p,j)

                    ! storage harvest mortality phosphorus fluxes
                    harvest_p_to_litr_met_p(c,j)      = harvest_p_to_litr_met_p(c,j)      + &
                         hrv_leafp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_p_to_litr_met_p(c,j)     = harvest_p_to_litr_met_p(c,j)     + &
                         hrv_frootp_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_p_to_litr_met_p(c,j)  = harvest_p_to_litr_met_p(c,j)  + &
                         hrv_livestemp_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_p_to_litr_met_p(c,j)  = harvest_p_to_litr_met_p(c,j)  + &
                         hrv_deadstemp_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_livecrootp_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_deadcrootp_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! transfer harvest mortality phosphorus fluxes
                    harvest_p_to_litr_met_p(c,j)      = harvest_p_to_litr_met_p(c,j)      + &
                         hrv_leafp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_p_to_litr_met_p(c,j)     = harvest_p_to_litr_met_p(c,j)     + &
                         hrv_frootp_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_p_to_litr_met_p(c,j)  = harvest_p_to_litr_met_p(c,j)  + &
                         hrv_livestemp_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_p_to_litr_met_p(c,j)  = harvest_p_to_litr_met_p(c,j)  + &
                         hrv_deadstemp_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_livecrootp_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_p_to_litr_met_p(c,j) = harvest_p_to_litr_met_p(c,j) + &
                         hrv_deadcrootp_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)

                 end if
              end if

           end do

        end do
     end do
   
     do pi = 1,maxpatch_pft
        do fc = 1,num_soilc
           c = filter_soilc(fc)

           if (pi <=  col_pp%npfts(c)) then
              p = col_pp%pfti(c) + pi - 1

              if (veg_pp%active(p)) then


                 ! wood harvest mortality carbon fluxes to product pools
                 chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                      phrv_deadstemc_to_prod10c(p)  * wtcol(p)
                 chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                      phrv_deadstemc_to_prod100c(p)  * wtcol(p)


                 ! wood harvest mortality nitrogen fluxes to product pools
                 chrv_deadstemn_to_prod10n(c)  = chrv_deadstemn_to_prod10n(c)  + &
                      phrv_deadstemn_to_prod10n(p)  * wtcol(p)
                 chrv_deadstemn_to_prod100n(c)  = chrv_deadstemn_to_prod100n(c)  + &
                      phrv_deadstemn_to_prod100n(p)  * wtcol(p)


                 ! wood harvest mortality phosphorus fluxes to product pools
                 chrv_deadstemp_to_prod10p(c)  = chrv_deadstemp_to_prod10p(c)  + &
                      phrv_deadstemp_to_prod10p(p)  * wtcol(p)
                 chrv_deadstemp_to_prod100p(c)  = chrv_deadstemp_to_prod100p(c)  + &
                      phrv_deadstemp_to_prod100p(p)  * wtcol(p)

              end if
           end if

        end do

     end do

   end associate 

 end subroutine CNHarvestPftToColumn

end module dynHarvestMod
