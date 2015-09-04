module dynHarvestMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the harvest data, as well as the state updates that happen as a
  ! result of harvest.
  !
  ! Currently, it is assumed that the harvest data are on the flanduse_timeseries file. However, this
  ! could theoretically be changed so that the harvest data were separated from the
  ! pftdyn data, allowing them to differ in the years over which they apply.
  !
  ! !USES:
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use abortutils            , only : endrun
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use CNStateType           , only : cnstate_type
  use CNCarbonStateType     , only : carbonstate_type
  use CNNitrogenStateType   , only : nitrogenstate_type
  use CNCarbonFluxType      , only : carbonflux_type
  use CNNitrogenFluxType    , only : nitrogenflux_type
  use EcophysConType        , only : ecophyscon
  use clm_varcon            , only : grlnd
  use ColumnType            , only : col                
  use PatchType             , only : pft                
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynHarvest_init    ! initialize data structures for harvest information
  public :: dynHarvest_interp  ! get harvest data for current time step, if needed
  public :: CNHarvest          ! harvest mortality routine for CN code
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNHarvestPftToColumn   ! gather pft-level harvest fluxes to the column level
  !
  ! !PRIVATE TYPES:

  ! Note that, since we have our own dynHarvest_file object (distinct from dynpft_file),
  ! we could theoretically have a separate file providing harvest data from that providing
  ! the pftdyn data
  type(dyn_file_type), target :: dynHarvest_file ! information for the file containing harvest data

  ! Define the underlying harvest variables
  integer, parameter :: num_harvest_vars = 5
  character(len=64), parameter :: harvest_varnames(num_harvest_vars) = &
       [character(len=64) :: 'HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 'HARVEST_SH2', 'HARVEST_SH3']
  
  type(dyn_var_time_uninterp_type) :: harvest_vars(num_harvest_vars)   ! value of each harvest variable

  real(r8) , allocatable   :: harvest(:) ! harvest rates
  logical                  :: do_harvest ! whether we're in a period when we should do harvest
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynHarvest_init(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize data structures for harvest information.
    ! This should be called once, during model initialization.
    ! 
    ! This also calls dynHarvest_interp for the initial time
    !
    ! !USES:
    use clm_varctl, only : use_cn, flanduse_timeseries
    use dynVarTimeUninterpMod   , only : dyn_var_time_uninterp_type
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer :: varnum     ! counter for harvest variables
    integer :: num_points ! number of spatial points
    integer :: ier        ! error code
    
    character(len=*), parameter :: subname = 'dynHarvest_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    allocate(harvest(bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for harvest'//errMsg(__FILE__, __LINE__))
    end if

    dynHarvest_file = dyn_file_type(flanduse_timeseries)
    
    ! Get initial harvest data
    if (use_cn) then
       num_points = (bounds%endg - bounds%begg + 1)
       do varnum = 1, num_harvest_vars
          harvest_vars(varnum) = dyn_var_time_uninterp_type( &
               dyn_file=dynHarvest_file, varname=harvest_varnames(varnum), &
               dim1name=grlnd, conversion_factor=1.0_r8, &
               do_check_sums_equal_1=.false., data_shape=[num_points])
       end do
       call dynHarvest_interp(bounds)
    end if
    
  end subroutine dynHarvest_init


  !-----------------------------------------------------------------------
  subroutine dynHarvest_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get harvest data for model time, when needed.
    !
    ! Note that harvest data are stored as rates (not weights) and so time interpolation
    ! is not necessary - the harvest rate is held constant through the year.  This is
    ! consistent with the treatment of changing PFT weights, where interpolation of the
    ! annual endpoint weights leads to a constant rate of change in PFT weight through the
    ! year, with abrupt changes in the rate at annual boundaries.
    !
    ! !USES:
    use clm_varctl     , only : use_cn
    use dynTimeInfoMod , only : time_info_type
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: varnum       ! counter for harvest variables
    real(r8), allocatable :: this_data(:) ! data for a single harvest variable

    character(len=*), parameter :: subname = 'dynHarvest_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! As a workaround for an internal compiler error with ifort 13.1.2 on goldbach, call
    ! the specific name of this procedure rather than using its generic name
    call dynHarvest_file%time_info%set_current_year_get_year()

    ! Get total harvest for this time step
    if (use_cn) then
       harvest(bounds%begg:bounds%endg) = 0._r8

       if (dynHarvest_file%time_info%is_before_time_series()) then
          ! Turn off harvest before the start of the harvest time series
          do_harvest = .false.
       else
          ! Note that do_harvest stays true even past the end of the time series. This
          ! means that harvest rates will be maintained at the rate given in the last
          ! year of the file for all years past the end of this specified time series.
          do_harvest = .true.
          allocate(this_data(bounds%begg:bounds%endg))
          do varnum = 1, num_harvest_vars
             call harvest_vars(varnum)%get_current_data(this_data)
             harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + &
                                                this_data(bounds%begg:bounds%endg)
          end do
          deallocate(this_data)
       end if
    end if

  end subroutine dynHarvest_interp


  !-----------------------------------------------------------------------
  subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! Harvest mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use pftvarcon       , only : noveg, nbrdlf_evr_shrub, pprodharv10
    use clm_varcon      , only : secspday
    use clm_time_manager, only : get_days_per_year
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p                         ! patch index
    integer :: g                         ! gridcell index
    integer :: fp                        ! patch filter index
    real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
    real(r8):: m                         ! rate for fractional harvest mortality (1/s)
    real(r8):: days_per_year             ! days per year
    !-----------------------------------------------------------------------

   associate(& 
   pgridcell                           =>    pft%gridcell                                , & ! Input:  [integer (:)]  pft-level index into gridcell-level quantities     
   ivt                                 =>    pft%itype                                   , & ! Input:  [integer (:)]  pft vegetation type                                

   leafc                               =>    carbonstate_vars%leafc_patch                                   , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   frootc                              =>    carbonstate_vars%frootc_patch                                  , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   livestemc                           =>    carbonstate_vars%livestemc_patch                               , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   deadstemc                           =>    carbonstate_vars%deadstemc_patch                               , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   livecrootc                          =>    carbonstate_vars%livecrootc_patch                              , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   deadcrootc                          =>    carbonstate_vars%deadcrootc_patch                              , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   xsmrpool                            =>    carbonstate_vars%xsmrpool_patch                                , & ! Input:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
   leafc_storage                       =>    carbonstate_vars%leafc_storage_patch                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   frootc_storage                      =>    carbonstate_vars%frootc_storage_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   livestemc_storage                   =>    carbonstate_vars%livestemc_storage_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   deadstemc_storage                   =>    carbonstate_vars%deadstemc_storage_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   livecrootc_storage                  =>    carbonstate_vars%livecrootc_storage_patch                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   deadcrootc_storage                  =>    carbonstate_vars%deadcrootc_storage_patch                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   gresp_storage                       =>    carbonstate_vars%gresp_storage_patch                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   leafc_xfer                          =>    carbonstate_vars%leafc_xfer_patch                              , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   frootc_xfer                         =>    carbonstate_vars%frootc_xfer_patch                             , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livestemc_xfer                      =>    carbonstate_vars%livestemc_xfer_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc_xfer                      =>    carbonstate_vars%deadstemc_xfer_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   livecrootc_xfer                     =>    carbonstate_vars%livecrootc_xfer_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc_xfer                     =>    carbonstate_vars%deadcrootc_xfer_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   gresp_xfer                          =>    carbonstate_vars%gresp_xfer_patch                              , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               

   leafn                               =>    nitrogenstate_vars%leafn_patch                                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   frootn                              =>    nitrogenstate_vars%frootn_patch                                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   livestemn                           =>    nitrogenstate_vars%livestemn_patch                               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   deadstemn                           =>    nitrogenstate_vars%deadstemn_patch                               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
   livecrootn                          =>    nitrogenstate_vars%livecrootn_patch                              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   deadcrootn                          =>    nitrogenstate_vars%deadcrootn_patch                              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   retransn                            =>    nitrogenstate_vars%retransn_patch                                , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   leafn_storage                       =>    nitrogenstate_vars%leafn_storage_patch                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   frootn_storage                      =>    nitrogenstate_vars%frootn_storage_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   livestemn_storage                   =>    nitrogenstate_vars%livestemn_storage_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   deadstemn_storage                   =>    nitrogenstate_vars%deadstemn_storage_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   livecrootn_storage                  =>    nitrogenstate_vars%livecrootn_storage_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   deadcrootn_storage                  =>    nitrogenstate_vars%deadcrootn_storage_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   leafn_xfer                          =>    nitrogenstate_vars%leafn_xfer_patch                              , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   frootn_xfer                         =>    nitrogenstate_vars%frootn_xfer_patch                             , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   livestemn_xfer                      =>    nitrogenstate_vars%livestemn_xfer_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   deadstemn_xfer                      =>    nitrogenstate_vars%deadstemn_xfer_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   livecrootn_xfer                     =>    nitrogenstate_vars%livecrootn_xfer_patch                         , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   deadcrootn_xfer                     =>    nitrogenstate_vars%deadcrootn_xfer_patch                         , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               

   hrv_leafc_to_litter                 =>    carbonflux_vars%hrv_leafc_to_litter_patch                     , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_to_litter                =>    carbonflux_vars%hrv_frootc_to_litter_patch                    , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_to_litter             =>    carbonflux_vars%hrv_livestemc_to_litter_patch                 , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod10c            =>    carbonflux_vars%hrv_deadstemc_to_prod10c_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod100c           =>    carbonflux_vars%hrv_deadstemc_to_prod100c_patch               , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_to_litter            =>    carbonflux_vars%hrv_livecrootc_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_to_litter            =>    carbonflux_vars%hrv_deadcrootc_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_xsmrpool_to_atm                 =>    carbonflux_vars%hrv_xsmrpool_to_atm_patch                     , & ! Output: [real(r8) (:)]                                                    
   hrv_leafc_storage_to_litter         =>    carbonflux_vars%hrv_leafc_storage_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_storage_to_litter        =>    carbonflux_vars%hrv_frootc_storage_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_storage_to_litter     =>    carbonflux_vars%hrv_livestemc_storage_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_storage_to_litter     =>    carbonflux_vars%hrv_deadstemc_storage_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_storage_to_litter    =>    carbonflux_vars%hrv_livecrootc_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_storage_to_litter    =>    carbonflux_vars%hrv_deadcrootc_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
   hrv_gresp_storage_to_litter         =>    carbonflux_vars%hrv_gresp_storage_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
   hrv_leafc_xfer_to_litter            =>    carbonflux_vars%hrv_leafc_xfer_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_xfer_to_litter           =>    carbonflux_vars%hrv_frootc_xfer_to_litter_patch               , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_xfer_to_litter        =>    carbonflux_vars%hrv_livestemc_xfer_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_xfer_to_litter        =>    carbonflux_vars%hrv_deadstemc_xfer_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_xfer_to_litter       =>    carbonflux_vars%hrv_livecrootc_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_xfer_to_litter       =>    carbonflux_vars%hrv_deadcrootc_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
   hrv_gresp_xfer_to_litter            =>    carbonflux_vars%hrv_gresp_xfer_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    

   hrv_leafn_to_litter                 =>    nitrogenflux_vars%hrv_leafn_to_litter_patch                     , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_to_litter                =>    nitrogenflux_vars%hrv_frootn_to_litter_patch                    , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_to_litter             =>    nitrogenflux_vars%hrv_livestemn_to_litter_patch                 , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod10n            =>    nitrogenflux_vars%hrv_deadstemn_to_prod10n_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod100n           =>    nitrogenflux_vars%hrv_deadstemn_to_prod100n_patch               , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_to_litter            =>    nitrogenflux_vars%hrv_livecrootn_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_to_litter            =>    nitrogenflux_vars%hrv_deadcrootn_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_retransn_to_litter              =>    nitrogenflux_vars%hrv_retransn_to_litter_patch                  , & ! Output: [real(r8) (:)]                                                    
   hrv_leafn_storage_to_litter         =>    nitrogenflux_vars%hrv_leafn_storage_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_storage_to_litter        =>    nitrogenflux_vars%hrv_frootn_storage_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_storage_to_litter     =>    nitrogenflux_vars%hrv_livestemn_storage_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_storage_to_litter     =>    nitrogenflux_vars%hrv_deadstemn_storage_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_storage_to_litter    =>    nitrogenflux_vars%hrv_livecrootn_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_storage_to_litter    =>    nitrogenflux_vars%hrv_deadcrootn_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
   hrv_leafn_xfer_to_litter            =>    nitrogenflux_vars%hrv_leafn_xfer_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_xfer_to_litter           =>    nitrogenflux_vars%hrv_frootn_xfer_to_litter_patch               , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_xfer_to_litter        =>    nitrogenflux_vars%hrv_livestemn_xfer_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_xfer_to_litter        =>    nitrogenflux_vars%hrv_deadstemn_xfer_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_xfer_to_litter       =>    nitrogenflux_vars%hrv_livecrootn_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_xfer_to_litter       =>    nitrogenflux_vars%hrv_deadcrootn_xfer_to_litter_patch             & ! Output: [real(r8) (:)]                                                    
   )


   days_per_year = get_days_per_year()

   ! patch loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      g = pgridcell(p)
      
      ! If this is a tree pft, then
      ! get the annual harvest "mortality" rate (am) from harvest array
      ! and convert to rate per second
      if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then

         if (do_harvest) then
            am = harvest(g)
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
         
      end if  ! end tree block

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes from harvest to the column
   ! for litter C and N inputs

   call CNHarvestPftToColumn(num_soilc, filter_soilc, &
      cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars)

    end associate 
 end subroutine CNHarvest

 !-----------------------------------------------------------------------
 subroutine CNHarvestPftToColumn (num_soilc, filter_soilc, &
      cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars)
   !
   ! !DESCRIPTION:
   ! called at the end of CNHarvest to gather all pft-level harvest litterfall fluxes
   ! to the column level and assign them to the three litter pools
   !
   ! !USES:
   use clm_varpar , only : maxpatch_pft, nlevdecomp
   !
   ! !ARGUMENTS:
   integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                  , intent(in)    :: filter_soilc(:) ! soil column filter
   type(cnstate_type)       , intent(in)    :: cnstate_vars
   type(carbonstate_type)   , intent(in)    :: carbonstate_vars
   type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
   type(carbonflux_type)    , intent(inout) :: carbonflux_vars
   type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
   !
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
   !-----------------------------------------------------------------------

   associate(                                                                                        & 
        ivt                              =>    pft%itype                                                , & ! Input:  [integer  (:)   ]  pft vegetation type                                
        wtcol                            =>    pft%wtcol                                                , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
        
        lf_flab                          =>    ecophyscon%lf_flab                                       , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
        lf_fcel                          =>    ecophyscon%lf_fcel                                       , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
        lf_flig                          =>    ecophyscon%lf_flig                                       , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
        fr_flab                          =>    ecophyscon%fr_flab                                       , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
        fr_fcel                          =>    ecophyscon%fr_fcel                                       , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
        fr_flig                          =>    ecophyscon%fr_flig                                       , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  
        
        leaf_prof                        =>    cnstate_vars%leaf_prof_patch                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
        froot_prof                       =>    cnstate_vars%froot_prof_patch                            , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
        croot_prof                       =>    cnstate_vars%croot_prof_patch                            , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
        stem_prof                        =>    cnstate_vars%stem_prof_patch                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
        
        hrv_leafc_to_litter              =>    carbonflux_vars%hrv_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_to_litter             =>    carbonflux_vars%hrv_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_to_litter          =>    carbonflux_vars%hrv_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemc_to_prod10c        =>    carbonflux_vars%hrv_deadstemc_to_prod10c_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemc_to_prod100c       =>    carbonflux_vars%hrv_deadstemc_to_prod100c_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_to_litter         =>    carbonflux_vars%hrv_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_to_litter         =>    carbonflux_vars%hrv_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_storage_to_litter      =>    carbonflux_vars%hrv_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_storage_to_litter     =>    carbonflux_vars%hrv_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_storage_to_litter  =>    carbonflux_vars%hrv_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_storage_to_litter  =>    carbonflux_vars%hrv_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_storage_to_litter =>    carbonflux_vars%hrv_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_storage_to_litter =>    carbonflux_vars%hrv_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_storage_to_litter      =>    carbonflux_vars%hrv_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_xfer_to_litter         =>    carbonflux_vars%hrv_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_xfer_to_litter        =>    carbonflux_vars%hrv_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_xfer_to_litter     =>    carbonflux_vars%hrv_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_xfer_to_litter     =>    carbonflux_vars%hrv_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_xfer_to_litter    =>    carbonflux_vars%hrv_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_xfer_to_litter    =>    carbonflux_vars%hrv_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_xfer_to_litter         =>    carbonflux_vars%hrv_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        chrv_deadstemc_to_prod10c        =>    carbonflux_vars%hrv_deadstemc_to_prod10c_col             , & ! InOut:  [real(r8) (:)   ]                                                    
        chrv_deadstemc_to_prod100c       =>    carbonflux_vars%hrv_deadstemc_to_prod100c_col            , & ! InOut:  [real(r8) (:)   ]                                                    
        harvest_c_to_litr_met_c          =>    carbonflux_vars%harvest_c_to_litr_met_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
        harvest_c_to_litr_cel_c          =>    carbonflux_vars%harvest_c_to_litr_cel_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
        harvest_c_to_litr_lig_c          =>    carbonflux_vars%harvest_c_to_litr_lig_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
        harvest_c_to_cwdc                =>    carbonflux_vars%harvest_c_to_cwdc_col                    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
        
        hrv_leafn_to_litter              =>    nitrogenflux_vars%hrv_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_to_litter             =>    nitrogenflux_vars%hrv_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_to_litter          =>    nitrogenflux_vars%hrv_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemn_to_prod10n        =>    nitrogenflux_vars%hrv_deadstemn_to_prod10n_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemn_to_prod100n       =>    nitrogenflux_vars%hrv_deadstemn_to_prod100n_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_to_litter         =>    nitrogenflux_vars%hrv_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_to_litter         =>    nitrogenflux_vars%hrv_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_retransn_to_litter           =>    nitrogenflux_vars%hrv_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_storage_to_litter      =>    nitrogenflux_vars%hrv_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_storage_to_litter     =>    nitrogenflux_vars%hrv_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_storage_to_litter  =>    nitrogenflux_vars%hrv_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_storage_to_litter  =>    nitrogenflux_vars%hrv_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_storage_to_litter =>    nitrogenflux_vars%hrv_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_storage_to_litter =>    nitrogenflux_vars%hrv_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_xfer_to_litter         =>    nitrogenflux_vars%hrv_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_xfer_to_litter        =>    nitrogenflux_vars%hrv_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_xfer_to_litter     =>    nitrogenflux_vars%hrv_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_xfer_to_litter     =>    nitrogenflux_vars%hrv_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_xfer_to_litter    =>    nitrogenflux_vars%hrv_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_xfer_to_litter    =>    nitrogenflux_vars%hrv_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        chrv_deadstemn_to_prod10n        =>    nitrogenflux_vars%hrv_deadstemn_to_prod10n_col           , & ! InOut:  [real(r8) (:)   ]                                                    
        chrv_deadstemn_to_prod100n       =>    nitrogenflux_vars%hrv_deadstemn_to_prod100n_col          , & ! InOut:  [real(r8) (:)   ]                                                    
        harvest_n_to_litr_met_n          =>    nitrogenflux_vars%harvest_n_to_litr_met_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
        harvest_n_to_litr_cel_n          =>    nitrogenflux_vars%harvest_n_to_litr_cel_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
        harvest_n_to_litr_lig_n          =>    nitrogenflux_vars%harvest_n_to_litr_lig_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter lignin pool (gN/m3/s)
        harvest_n_to_cwdn                =>    nitrogenflux_vars%harvest_n_to_cwdn_col                    & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to CWD pool (gN/m3/s)
        
        )

     do j = 1, nlevdecomp
        do pi = 1,maxpatch_pft
           do fc = 1,num_soilc
              c = filter_soilc(fc)

              if (pi <=  col%npfts(c)) then
                 p = col%pfti(c) + pi - 1

                 if (pft%active(p)) then

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

                 end if
              end if

           end do

        end do
     end do
   
     do pi = 1,maxpatch_pft
        do fc = 1,num_soilc
           c = filter_soilc(fc)

           if (pi <=  col%npfts(c)) then
              p = col%pfti(c) + pi - 1

              if (pft%active(p)) then


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
              end if
           end if

        end do

     end do

   end associate 

 end subroutine CNHarvestPftToColumn

end module dynHarvestMod
