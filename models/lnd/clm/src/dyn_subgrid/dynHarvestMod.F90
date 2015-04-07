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
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use decompMod               , only : bounds_type, BOUNDS_LEVEL_PROC
  use abortutils              , only : endrun
  use dynFileMod              , only : dyn_file_type
  use dynVarTimeUninterpMod   , only : dyn_var_time_uninterp_type
  use CNVegCarbonStateType    , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType     , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType  , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType   , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType , only : soilbiogeochem_state_type
  use pftconMod               , only : pftcon
  use clm_varcon              , only : grlnd
  use ColumnType              , only : col                
  use PatchType               , only : patch                
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynHarvest_init    ! initialize data structures for harvest information
  public :: dynHarvest_interp  ! get harvest data for current time step, if needed
  public :: CNHarvest          ! harvest mortality routine for CN code
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNHarvestPftToColumn   ! gather patch-level harvest fluxes to the column level
  !
  ! !PRIVATE TYPES:

  ! Note that, since we have our own dynHarvest_file object (distinct from dynpft_file),
  ! we could theoretically have a separate file providing harvest data from that providing
  ! the pftdyn data
  type(dyn_file_type), target :: dynHarvest_file ! information for the file containing harvest data

  ! Define the underlying harvest variables
  integer, parameter :: num_harvest_inst = 5
  character(len=64), parameter :: harvest_varnames(num_harvest_inst) = &
       [character(len=64) :: 'HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 'HARVEST_SH2', 'HARVEST_SH3']
  
  type(dyn_var_time_uninterp_type) :: harvest_inst(num_harvest_inst)   ! value of each harvest variable

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
    use clm_varctl            , only : use_cn, flanduse_timeseries
    use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
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
       do varnum = 1, num_harvest_inst
          harvest_inst(varnum) = dyn_var_time_uninterp_type( &
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
          do varnum = 1, num_harvest_inst
             call harvest_inst(varnum)%get_current_data(this_data)
             harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + &
                                                this_data(bounds%begg:bounds%endg)
          end do
          deallocate(this_data)
       end if
    end if

  end subroutine dynHarvest_interp


  !-----------------------------------------------------------------------
  subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Harvest mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use pftconMod       , only : noveg, nbrdlf_evr_shrub
    use clm_varcon      , only : secspday
    use clm_time_manager, only : get_days_per_year
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
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
         ivt                                 =>    patch%itype                                                      , & ! Input:  [integer (:)]  pft vegetation type                                
         
         leafc                               =>    cnveg_carbonstate_inst%leafc_patch                             , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
         frootc                              =>    cnveg_carbonstate_inst%frootc_patch                            , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
         livestemc                           =>    cnveg_carbonstate_inst%livestemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
         deadstemc                           =>    cnveg_carbonstate_inst%deadstemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
         livecrootc                          =>    cnveg_carbonstate_inst%livecrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
         deadcrootc                          =>    cnveg_carbonstate_inst%deadcrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
         xsmrpool                            =>    cnveg_carbonstate_inst%xsmrpool_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
         gresp_xfer                          =>    cnveg_carbonstate_inst%gresp_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
         
         leafn                               =>    cnveg_nitrogenstate_inst%leafn_patch                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
         frootn                              =>    cnveg_nitrogenstate_inst%frootn_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
         livestemn                           =>    cnveg_nitrogenstate_inst%livestemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
         deadstemn                           =>    cnveg_nitrogenstate_inst%deadstemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
         livecrootn                          =>    cnveg_nitrogenstate_inst%livecrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
         deadcrootn                          =>    cnveg_nitrogenstate_inst%deadcrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
         retransn                            =>    cnveg_nitrogenstate_inst%retransn_patch                        , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
         
         hrv_leafc_to_litter                 =>    cnveg_carbonflux_inst%hrv_leafc_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_to_litter                =>    cnveg_carbonflux_inst%hrv_frootc_to_litter_patch               , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_to_litter             =>    cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_to_prod10c            =>    cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_to_prod100c           =>    cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_patch          , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_xsmrpool_to_atm                 =>    cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_patch                , & ! Output: [real(r8) (:)]                                                    
         hrv_leafc_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_storage_to_litter        =>    cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootc_storage_to_litter    =>    cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_storage_to_litter    =>    cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         hrv_gresp_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_leafc_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_xfer_to_litter           =>    cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_gresp_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         
         hrv_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%hrv_leafn_to_litter_patch              , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_to_litter                =>    cnveg_nitrogenflux_inst%hrv_frootn_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_to_litter             =>    cnveg_nitrogenflux_inst%hrv_livestemn_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_to_prod10n            =>    cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod10n_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_to_prod100n           =>    cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod100n_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_livecrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_retransn_to_litter              =>    cnveg_nitrogenflux_inst%hrv_retransn_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_leafn_storage_to_litter         =>    cnveg_nitrogenflux_inst%hrv_leafn_storage_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_storage_to_litter        =>    cnveg_nitrogenflux_inst%hrv_frootn_storage_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_livestemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_deadstemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%hrv_livecrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         hrv_leafn_xfer_to_litter            =>    cnveg_nitrogenflux_inst%hrv_leafn_xfer_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_xfer_to_litter           =>    cnveg_nitrogenflux_inst%hrv_frootn_xfer_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_livestemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_deadstemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%hrv_livecrootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_xfer_to_litter_patch      & ! Output: [real(r8) (:)]                                                    
         )

      
      days_per_year = get_days_per_year()

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         g = patch%gridcell(p)

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

            ! patch-level harvest carbon fluxes
            ! displayed pools
            hrv_leafc_to_litter(p)               = leafc(p)               * m
            hrv_frootc_to_litter(p)              = frootc(p)              * m
            hrv_livestemc_to_litter(p)           = livestemc(p)           * m
            hrv_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * pftcon%pprodharv10(ivt(p))
            hrv_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * (1.0_r8 - pftcon%pprodharv10(ivt(p)))
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

            ! patch-level harvest mortality nitrogen fluxes
            ! displayed pools
            hrv_leafn_to_litter(p)               = leafn(p)               * m
            hrv_frootn_to_litter(p)              = frootn(p)              * m
            hrv_livestemn_to_litter(p)           = livestemn(p)           * m
            hrv_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * pftcon%pprodharv10(ivt(p))
            hrv_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * (1.0_r8 - pftcon%pprodharv10(ivt(p)))
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

      ! gather all patch-level litterfall fluxes from harvest to the column
      ! for litter C and N inputs

      call CNHarvestPftToColumn(num_soilc, filter_soilc, &
           soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

    end associate 

  end subroutine CNHarvest

 !-----------------------------------------------------------------------
 subroutine CNHarvestPftToColumn (num_soilc, filter_soilc, &
      soilbiogeochem_state_inst, CNVeg_carbonflux_inst, cnveg_nitrogenflux_inst)
   !
   ! !DESCRIPTION:
   ! called at the end of CNHarvest to gather all patch-level harvest litterfall fluxes
   ! to the column level and assign them to the three litter pools
   !
   ! !USES:
   use clm_varpar , only : maxpatch_pft, nlevdecomp
   !
   ! !ARGUMENTS:
   integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                         , intent(in)    :: filter_soilc(:) ! soil column filter
   type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
   type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
   !
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
   !-----------------------------------------------------------------------

   associate(                                                                                                   & 
        ivt                              =>    patch%itype                                                      , & ! Input:  [integer  (:)   ]  pft vegetation type                                
        wtcol                            =>    patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
        
        lf_flab                          =>    pftcon%lf_flab                                                 , & ! Input:  leaf litter labile fraction                       
        lf_fcel                          =>    pftcon%lf_fcel                                                 , & ! Input:  leaf litter cellulose fraction                    
        lf_flig                          =>    pftcon%lf_flig                                                 , & ! Input:  leaf litter lignin fraction                       
        fr_flab                          =>    pftcon%fr_flab                                                 , & ! Input:  fine root litter labile fraction                  
        fr_fcel                          =>    pftcon%fr_fcel                                                 , & ! Input:  fine root litter cellulose fraction               
        fr_flig                          =>    pftcon%fr_flig                                                 , & ! Input:  fine root litter lignin fraction                  
        
        leaf_prof                        =>    soilbiogeochem_state_inst%leaf_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
        froot_prof                       =>    soilbiogeochem_state_inst%froot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
        croot_prof                       =>    soilbiogeochem_state_inst%croot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
        stem_prof                        =>    soilbiogeochem_state_inst%stem_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
        
        hrv_leafc_to_litter              =>    cnveg_carbonflux_inst%hrv_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_to_litter             =>    cnveg_carbonflux_inst%hrv_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_to_litter          =>    cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemc_to_prod10c        =>    cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemc_to_prod100c       =>    cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_to_litter         =>    cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_to_litter         =>    cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_storage_to_litter      =>    cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_storage_to_litter  =>    cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_storage_to_litter  =>    cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_storage_to_litter =>    cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_storage_to_litter =>    cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_storage_to_litter      =>    cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_xfer_to_litter         =>    cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_xfer_to_litter     =>    cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_xfer_to_litter     =>    cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_xfer_to_litter    =>    cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_xfer_to_litter    =>    cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_xfer_to_litter         =>    cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        chrv_deadstemc_to_prod10c        =>    cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_col             , & ! InOut:  [real(r8) (:)   ]                                                    
        chrv_deadstemc_to_prod100c       =>    cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_col            , & ! InOut:  [real(r8) (:)   ]                                                    
        harvest_c_to_litr_met_c          =>    cnveg_carbonflux_inst%harvest_c_to_litr_met_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
        harvest_c_to_litr_cel_c          =>    cnveg_carbonflux_inst%harvest_c_to_litr_cel_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
        harvest_c_to_litr_lig_c          =>    cnveg_carbonflux_inst%harvest_c_to_litr_lig_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
        harvest_c_to_cwdc                =>    cnveg_carbonflux_inst%harvest_c_to_cwdc_col                    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
        
        hrv_leafn_to_litter              =>    cnveg_nitrogenflux_inst%hrv_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_to_litter             =>    cnveg_nitrogenflux_inst%hrv_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_to_litter          =>    cnveg_nitrogenflux_inst%hrv_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemn_to_prod10n        =>    cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod10n_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        phrv_deadstemn_to_prod100n       =>    cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod100n_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_to_litter         =>    cnveg_nitrogenflux_inst%hrv_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_to_litter         =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_retransn_to_litter           =>    cnveg_nitrogenflux_inst%hrv_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_storage_to_litter      =>    cnveg_nitrogenflux_inst%hrv_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_storage_to_litter  =>    cnveg_nitrogenflux_inst%hrv_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_storage_to_litter  =>    cnveg_nitrogenflux_inst%hrv_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_storage_to_litter =>    cnveg_nitrogenflux_inst%hrv_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_storage_to_litter =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_xfer_to_litter         =>    cnveg_nitrogenflux_inst%hrv_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_xfer_to_litter     =>    cnveg_nitrogenflux_inst%hrv_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_xfer_to_litter     =>    cnveg_nitrogenflux_inst%hrv_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%hrv_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        chrv_deadstemn_to_prod10n        =>    cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod10n_col           , & ! InOut:  [real(r8) (:)   ]                                                    
        chrv_deadstemn_to_prod100n       =>    cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod100n_col          , & ! InOut:  [real(r8) (:)   ]                                                    
        harvest_n_to_litr_met_n          =>    cnveg_nitrogenflux_inst%harvest_n_to_litr_met_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
        harvest_n_to_litr_cel_n          =>    cnveg_nitrogenflux_inst%harvest_n_to_litr_cel_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
        harvest_n_to_litr_lig_n          =>    cnveg_nitrogenflux_inst%harvest_n_to_litr_lig_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter lignin pool (gN/m3/s)
        harvest_n_to_cwdn                =>    cnveg_nitrogenflux_inst%harvest_n_to_cwdn_col                    & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to CWD pool (gN/m3/s)
        )

     do j = 1, nlevdecomp
        do pi = 1,maxpatch_pft
           do fc = 1,num_soilc
              c = filter_soilc(fc)

              if (pi <=  col%npatches(c)) then
                 p = col%patchi(c) + pi - 1

                 if (patch%active(p)) then

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

           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1

              if (patch%active(p)) then


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
