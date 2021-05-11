module CarbonIsoFluxMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon isotopic flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varpar             , only : ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  use elm_varpar             , only : max_patch_per_col, maxpatch_pft
  use abortutils             , only : endrun
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use CNCarbonFluxType       , only : carbonflux_type
  use CNCarbonStateType      , only : carbonstate_type
  use CNStateType            , only : cnstate_type
  use ColumnType             , only : col_pp 
  use ColumnDataType         , only : column_carbon_state, column_carbon_flux
  use ColumnDataType         , only : col_cs, col_cf  
  use VegetationType         , only : veg_pp                
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_carbon_flux 
  use VegetationDataType     , only : veg_cs, veg_cf  
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: CarbonIsoFlux1
  public  :: CarbonIsoFlux2
  public  :: CarbonIsoFlux2h
  public  :: CarbonIsoFlux3
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNCIsoLitterToColumn
  private :: CNCIsoGapPftToColumn
  private :: CNCIsoHarvestPftToColumn
  private :: CarbonIsoFluxCalc
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       isotopeflux_vars, isotopestate_vars, isotope, isocol_cs, isoveg_cs, isocol_cf, isoveg_cf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic flux
    ! variables (except for gap-phase mortality and fire fluxes)
    use tracer_varcon, only : is_active_betr_bgc  
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(in)    :: carbonstate_vars
    type(carbonflux_type)  , intent(inout) :: isotopeflux_vars
    type(carbonstate_type) , intent(in)    :: isotopestate_vars
    character(len=*)       , intent(in)    :: isotope         ! 'c13' or 'c14'
    type(column_carbon_state),intent(in)   :: isocol_cs
    type(vegetation_carbon_state),intent(in)   :: isoveg_cs
    type(column_carbon_flux),intent(inout)    :: isocol_cf
    type(vegetation_carbon_flux),intent(inout):: isoveg_cf
    !
    ! !LOCAL VARIABLES:
    integer :: fp,pi,l,fc,cc,j
    !-----------------------------------------------------------------------

    associate(&
    cascade_donor_pool  =>  decomp_cascade_con%cascade_donor_pool  & !  [integer (:)]  which pool is C taken from for a given decomposition step 
    )

      ! patch-level non-mortality fluxes
   
      ! Note: if the variables which are arguments to CarbonIsoFluxCalc are ever changed to NOT be
      ! pointers, then the CarbonIsoFluxCalc routine will need to be changed to declare the bounds
      ! of each argument, these bounds will need to be passed in, and - importantly for
      ! threading to work properly - the subroutine calls will need to be changed so that
      ! instead of 'call CarbonIsoFluxCalc(foo, ...)' we have 'call CarbonIsoFluxCalc(foo(begp:endp), ...)'.
      
      call CarbonIsoFluxCalc(&
           isoveg_cf%leafc_xfer_to_leafc           , veg_cf%leafc_xfer_to_leafc, &
           isoveg_cs%leafc_xfer                   , veg_cs%leafc_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%frootc_xfer_to_frootc         , veg_cf%frootc_xfer_to_frootc, &
           isoveg_cs%frootc_xfer                  , veg_cs%frootc_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livestemc_xfer_to_livestemc   , veg_cf%livestemc_xfer_to_livestemc, &
           isoveg_cs%livestemc_xfer               , veg_cs%livestemc_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%deadstemc_xfer_to_deadstemc   , veg_cf%deadstemc_xfer_to_deadstemc, &
           isoveg_cs%deadstemc_xfer               , veg_cs%deadstemc_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livecrootc_xfer_to_livecrootc , veg_cf%livecrootc_xfer_to_livecrootc, &
           isoveg_cs%livecrootc_xfer              , veg_cs%livecrootc_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%deadcrootc_xfer_to_deadcrootc , veg_cf%deadcrootc_xfer_to_deadcrootc, &
           isoveg_cs%deadcrootc_xfer              , veg_cs%deadcrootc_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%leafc_to_litter               , veg_cf%leafc_to_litter, &
           isoveg_cs%leafc                        , veg_cs%leafc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%frootc_to_litter              , veg_cf%frootc_to_litter, &
           isoveg_cs%frootc                       , veg_cs%frootc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livestemc_to_deadstemc        , veg_cf%livestemc_to_deadstemc, &
           isoveg_cs%livestemc                    , veg_cs%livestemc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livecrootc_to_deadcrootc      , veg_cf%livecrootc_to_deadcrootc, &
           isoveg_cs%livecrootc                   , veg_cs%livecrootc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%leaf_curmr                    , veg_cf%leaf_curmr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%froot_curmr                   , veg_cf%froot_curmr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livestem_curmr                , veg_cf%livestem_curmr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livecroot_curmr               , veg_cf%livecroot_curmr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%xr                            , veg_cf%xr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)  

      call CarbonIsoFluxCalc(&
           isoveg_cf%leaf_xsmr                     , veg_cf%leaf_xsmr, &
           isoveg_cs%totvegc                      , veg_cs%totvegc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%froot_xsmr                    , veg_cf%froot_xsmr, &
           isoveg_cs%totvegc                      , veg_cs%totvegc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livestem_xsmr                 , veg_cf%livestem_xsmr, &
           isoveg_cs%totvegc                      , veg_cs%totvegc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livecroot_xsmr                , veg_cf%livecroot_xsmr, &
           isoveg_cs%totvegc                      , veg_cs%totvegc, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_xsmrpool             , veg_cf%cpool_to_xsmrpool, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_leafc                , veg_cf%cpool_to_leafc, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_leafc_storage        , veg_cf%cpool_to_leafc_storage, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_frootc               , veg_cf%cpool_to_frootc, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_frootc_storage       , veg_cf%cpool_to_frootc_storage, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_livestemc            , veg_cf%cpool_to_livestemc, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_livestemc_storage    , veg_cf%cpool_to_livestemc_storage, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_deadstemc            , veg_cf%cpool_to_deadstemc, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_deadstemc_storage    , veg_cf%cpool_to_deadstemc_storage, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_livecrootc           , veg_cf%cpool_to_livecrootc, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_livecrootc_storage   , veg_cf%cpool_to_livecrootc_storage, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_deadcrootc           , veg_cf%cpool_to_deadcrootc, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_deadcrootc_storage   , veg_cf%cpool_to_deadcrootc_storage, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_leaf_gr                 , veg_cf%cpool_leaf_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_froot_gr                , veg_cf%cpool_froot_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_livestem_gr             , veg_cf%cpool_livestem_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_deadstem_gr             , veg_cf%cpool_deadstem_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_livecroot_gr            , veg_cf%cpool_livecroot_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_deadcroot_gr            , veg_cf%cpool_deadcroot_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_leaf_storage_gr         , veg_cf%cpool_leaf_storage_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_froot_storage_gr        , veg_cf%cpool_froot_storage_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_livestem_storage_gr     , veg_cf%cpool_livestem_storage_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_deadstem_storage_gr     , veg_cf%cpool_deadstem_storage_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_livecroot_storage_gr    , veg_cf%cpool_livecroot_storage_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_deadcroot_storage_gr    , veg_cf%cpool_deadcroot_storage_gr, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%cpool_to_gresp_storage        , veg_cf%cpool_to_gresp_storage, &
           isoveg_cs%cpool                        , veg_cs%cpool, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%transfer_leaf_gr              , veg_cf%transfer_leaf_gr, &
           isoveg_cs%gresp_xfer                   , veg_cs%gresp_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%transfer_froot_gr             , veg_cf%transfer_froot_gr, &
           isoveg_cs%gresp_xfer                   , veg_cs%gresp_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%transfer_livestem_gr          , veg_cf%transfer_livestem_gr, &
           isoveg_cs%gresp_xfer                   , veg_cs%gresp_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%transfer_deadstem_gr          , veg_cf%transfer_deadstem_gr, &
           isoveg_cs%gresp_xfer                   , veg_cs%gresp_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%transfer_livecroot_gr         , veg_cf%transfer_livecroot_gr, &
           isoveg_cs%gresp_xfer                   , veg_cs%gresp_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%transfer_deadcroot_gr         , veg_cf%transfer_deadcroot_gr, &
           isoveg_cs%gresp_xfer                   , veg_cs%gresp_xfer, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%leafc_storage_to_xfer         , veg_cf%leafc_storage_to_xfer, &
           isoveg_cs%leafc_storage                , veg_cs%leafc_storage, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%frootc_storage_to_xfer        , veg_cf%frootc_storage_to_xfer, &
           isoveg_cs%frootc_storage               , veg_cs%frootc_storage, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livestemc_storage_to_xfer     , veg_cf%livestemc_storage_to_xfer, &
           isoveg_cs%livestemc_storage            , veg_cs%livestemc_storage, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%deadstemc_storage_to_xfer     , veg_cf%deadstemc_storage_to_xfer, &
           isoveg_cs%deadstemc_storage            , veg_cs%deadstemc_storage, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%livecrootc_storage_to_xfer    , veg_cf%livecrootc_storage_to_xfer, &
           isoveg_cs%livecrootc_storage           , veg_cs%livecrootc_storage, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%deadcrootc_storage_to_xfer    , veg_cf%deadcrootc_storage_to_xfer, &
           isoveg_cs%deadcrootc_storage           , veg_cs%deadcrootc_storage, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%gresp_storage_to_xfer         , veg_cf%gresp_storage_to_xfer, &
           isoveg_cs%gresp_storage                , veg_cs%gresp_storage, &
           num_soilp                                            , filter_soilp, 1._r8, 0, isotope)

      ! call routine to shift patch-level litterfall fluxes to column, for isotopes
      ! the non-isotope version of this routine is called in PhenologyMod.F90.F90
      ! For later clean-up, it would be possible to generalize this function to operate on a single 
      ! patch-to-column flux.

      call CNCIsoLitterToColumn(num_soilc, filter_soilc, cnstate_vars, isotopeflux_vars, isocol_cf, isoveg_cf)

      if (.not. is_active_betr_bgc) then

         ! column-level non-mortality fluxes

         do fc = 1,num_soilc
            cc = filter_soilc(fc)
            do j = 1, nlevdecomp
               do l = 1, ndecomp_cascade_transitions
                  if ( col_cs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) /= 0._r8) then
                     isocol_cf%decomp_cascade_hr_vr(cc,j,l)  =  &
                          col_cf%decomp_cascade_hr_vr(cc,j,l) * &
                          (isocol_cs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) &
                         / col_cs%decomp_cpools_vr(cc,j,cascade_donor_pool(l))) * 1._r8
                  else
                     isocol_cf%decomp_cascade_hr_vr(cc,j,l) = 0._r8
                  end if
               end do
            end do
         end do

         do fc = 1,num_soilc
            cc = filter_soilc(fc)
            do j = 1, nlevdecomp
               do l = 1, ndecomp_cascade_transitions
                  if ( col_cs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) /= 0._r8) then
                     isocol_cf%decomp_cascade_ctransfer_vr(cc,j,l)  =  &
                          col_cf%decomp_cascade_ctransfer_vr(cc,j,l) * &
                          (isocol_cs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) &
                          / col_cs%decomp_cpools_vr(cc,j,cascade_donor_pool(l))) * 1._r8
                  else
                     isocol_cf%decomp_cascade_ctransfer_vr(cc,j,l) = 0._r8
                  end if
               end do
            end do
         end do
    endif
    end associate

  end subroutine CarbonIsoFlux1

  !-----------------------------------------------------------------------
  subroutine CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       isotopeflux_vars, isotopestate_vars, isotope, isocol_cs, isoveg_cs, isocol_cf, isoveg_cf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic fluxes for gap mortality
    !
    use tracer_varcon, only : is_active_betr_bgc      
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(in)    :: carbonstate_vars
    type(carbonflux_type)  , intent(inout) :: isotopeflux_vars
    type(carbonstate_type) , intent(in)    :: isotopestate_vars
    character(len=*)       , intent(in)    :: isotope         ! 'c13' or 'c14'
    type(column_carbon_state),intent(in)   :: isocol_cs
    type(vegetation_carbon_state),intent(in)   :: isoveg_cs
    type(column_carbon_flux),intent(inout)    :: isocol_cf
    type(vegetation_carbon_flux),intent(inout):: isoveg_cf
    !
    ! !LOCAL VARIABLES:
    integer :: fp,pi
    !-----------------------------------------------------------------------

    ! patch-level gap mortality fluxes
   
    call CarbonIsoFluxCalc(&
         isoveg_cf%m_leafc_to_litter                    , veg_cf%m_leafc_to_litter, &
         isoveg_cs%leafc                               , veg_cs%leafc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_leafc_storage_to_litter            , veg_cf%m_leafc_storage_to_litter, &
         isoveg_cs%leafc_storage                       , veg_cs%leafc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_leafc_xfer_to_litter               , veg_cf%m_leafc_xfer_to_litter, &
         isoveg_cs%leafc_xfer                          , veg_cs%leafc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_frootc_to_litter                   , veg_cf%m_frootc_to_litter, &
         isoveg_cs%frootc                              , veg_cs%frootc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_frootc_storage_to_litter           , veg_cf%m_frootc_storage_to_litter, &
         isoveg_cs%frootc_storage                      , veg_cs%frootc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_frootc_xfer_to_litter              , veg_cf%m_frootc_xfer_to_litter, &
         isoveg_cs%frootc_xfer                         , veg_cs%frootc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_livestemc_to_litter                , veg_cf%m_livestemc_to_litter, &
         isoveg_cs%livestemc                           , veg_cs%livestemc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_livestemc_storage_to_litter        , veg_cf%m_livestemc_storage_to_litter, &
         isoveg_cs%livestemc_storage                   , veg_cs%livestemc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_livestemc_xfer_to_litter           , veg_cf%m_livestemc_xfer_to_litter, &
         isoveg_cs%livestemc_xfer                      , veg_cs%livestemc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_deadstemc_to_litter                , veg_cf%m_deadstemc_to_litter, &
         isoveg_cs%deadstemc                           , veg_cs%deadstemc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_deadstemc_storage_to_litter        , veg_cf%m_deadstemc_storage_to_litter, &
         isoveg_cs%deadstemc_storage                   , veg_cs%deadstemc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_deadstemc_xfer_to_litter           , veg_cf%m_deadstemc_xfer_to_litter, &
         isoveg_cs%deadstemc_xfer                      , veg_cs%deadstemc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_livecrootc_to_litter               , veg_cf%m_livecrootc_to_litter, &
         isoveg_cs%livecrootc                          , veg_cs%livecrootc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_livecrootc_storage_to_litter       , veg_cf%m_livecrootc_storage_to_litter, &
         isoveg_cs%livecrootc_storage                  , veg_cs%livecrootc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_livecrootc_xfer_to_litter          , veg_cf%m_livecrootc_xfer_to_litter, &
         isoveg_cs%livecrootc_xfer                     , veg_cs%livecrootc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_deadcrootc_to_litter               , veg_cf%m_deadcrootc_to_litter, &
         isoveg_cs%deadcrootc                          , veg_cs%deadcrootc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_deadcrootc_storage_to_litter       , veg_cf%m_deadcrootc_storage_to_litter, &
         isoveg_cs%deadcrootc_storage                  , veg_cs%deadcrootc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_deadcrootc_xfer_to_litter          , veg_cf%m_deadcrootc_xfer_to_litter, &
         isoveg_cs%deadcrootc_xfer                     , veg_cs%deadcrootc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_gresp_storage_to_litter            , veg_cf%m_gresp_storage_to_litter, &
         isoveg_cs%gresp_storage                       , veg_cs%gresp_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_gresp_xfer_to_litter               , veg_cf%m_gresp_xfer_to_litter, &
         isoveg_cs%gresp_xfer                          , veg_cs%gresp_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%m_cpool_to_litter               , veg_cf%m_cpool_to_litter, &
         isoveg_cs%cpool                          , veg_cs%cpool, &
         num_soilp                                              , filter_soilp, 1._r8, 0, isotope)            

    ! call routine to shift patch-level gap mortality fluxes to column , for isotopes
    ! the non-isotope version of this routine is in GapMortalityMod.F90.

    call CNCIsoGapPftToColumn(num_soilc, filter_soilc, cnstate_vars, isotopeflux_vars, isocol_cf, isoveg_cf)

  end subroutine CarbonIsoFlux2

  !-----------------------------------------------------------------------
  subroutine CarbonIsoFlux2h(num_soilc , filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       isotopeflux_vars, isotopestate_vars, isotope, isocol_cs, isoveg_cs, isocol_cf, isoveg_cf)
    !
    ! !DESCRIPTION:
    ! set the carbon isotopic fluxes for harvest mortality
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(in)    :: carbonstate_vars
    type(carbonflux_type)  , intent(inout) :: isotopeflux_vars
    type(carbonstate_type) , intent(in)    :: isotopestate_vars
    character(len=*)       , intent(in)    :: isotope         ! 'c13' or 'c14'
    type(column_carbon_state),intent(in)   :: isocol_cs
    type(vegetation_carbon_state),intent(in)   :: isoveg_cs
    type(column_carbon_flux),intent(inout)    :: isocol_cf
    type(vegetation_carbon_flux),intent(inout):: isoveg_cf
    !-----------------------------------------------------------------------

    ! patch-level gap mortality fluxes
   

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_leafc_to_litter                  , veg_cf%hrv_leafc_to_litter, &
         isoveg_cs%leafc                               , veg_cs%leafc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_leafc_storage_to_litter          , veg_cf%hrv_leafc_storage_to_litter, &
         isoveg_cs%leafc_storage                       , veg_cs%leafc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_leafc_xfer_to_litter             , veg_cf%hrv_leafc_xfer_to_litter, &
         isoveg_cs%leafc_xfer                          , veg_cs%leafc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_frootc_to_litter                 , veg_cf%hrv_frootc_to_litter, &
         isoveg_cs%frootc                              , veg_cs%frootc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_frootc_storage_to_litter         , veg_cf%hrv_frootc_storage_to_litter, &
         isoveg_cs%frootc_storage                      , veg_cs%frootc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_frootc_xfer_to_litter            , veg_cf%hrv_frootc_xfer_to_litter, &
         isoveg_cs%frootc_xfer                         , veg_cs%frootc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_livestemc_to_litter              , veg_cf%hrv_livestemc_to_litter, &
         isoveg_cs%livestemc                           , veg_cs%livestemc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_livestemc_storage_to_litter      , veg_cf%hrv_livestemc_storage_to_litter, &
         isoveg_cs%livestemc_storage                   , veg_cs%livestemc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_livestemc_xfer_to_litter         , veg_cf%hrv_livestemc_xfer_to_litter, &
         isoveg_cs%livestemc_xfer                      , veg_cs%livestemc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_deadstemc_to_prod10c             , veg_cf%hrv_deadstemc_to_prod10c, &
         isoveg_cs%deadstemc                           , veg_cs%deadstemc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_deadstemc_to_prod100c            , veg_cf%hrv_deadstemc_to_prod100c, &
         isoveg_cs%deadstemc                           , veg_cs%deadstemc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_deadstemc_storage_to_litter      , veg_cf%hrv_deadstemc_storage_to_litter, &
         isoveg_cs%deadstemc_storage                   , veg_cs%deadstemc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_deadstemc_xfer_to_litter         , veg_cf%hrv_deadstemc_xfer_to_litter, &
         isoveg_cs%deadstemc_xfer                      , veg_cs%deadstemc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_livecrootc_to_litter             , veg_cf%hrv_livecrootc_to_litter, &
         isoveg_cs%livecrootc                          , veg_cs%livecrootc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_livecrootc_storage_to_litter     , veg_cf%hrv_livecrootc_storage_to_litter, &
         isoveg_cs%livecrootc_storage                  , veg_cs%livecrootc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_livecrootc_xfer_to_litter        , veg_cf%hrv_livecrootc_xfer_to_litter, &
         isoveg_cs%livecrootc_xfer                     , veg_cs%livecrootc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_deadcrootc_to_litter             , veg_cf%hrv_deadcrootc_to_litter, &
         isoveg_cs%deadcrootc                          , veg_cs%deadcrootc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_deadcrootc_storage_to_litter     , veg_cf%hrv_deadcrootc_storage_to_litter, &
         isoveg_cs%deadcrootc_storage                  , veg_cs%deadcrootc_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_deadcrootc_xfer_to_litter        , veg_cf%hrv_deadcrootc_xfer_to_litter, &
         isoveg_cs%deadcrootc_xfer                     , veg_cs%deadcrootc_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_gresp_storage_to_litter          , veg_cf%hrv_gresp_storage_to_litter, &
         isoveg_cs%gresp_storage                       , veg_cs%gresp_storage, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_gresp_xfer_to_litter             , veg_cf%hrv_gresp_xfer_to_litter, &
         isoveg_cs%gresp_xfer                          , veg_cs%gresp_xfer, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_xsmrpool_to_atm                  , veg_cf%hrv_xsmrpool_to_atm, &
         isoveg_cs%totvegc                             , veg_cs%totvegc, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    call CarbonIsoFluxCalc(&
         isoveg_cf%hrv_cpool_to_litter                  , veg_cf%hrv_cpool_to_litter, &
         isoveg_cs%cpool                               , veg_cs%cpool, &
         num_soilp                                                   , filter_soilp, 1._r8, 0, isotope)

    ! call routine to shift patch-level gap mortality fluxes to column, for isotopes
    ! the non-isotope version of this routine is in GapMortalityMod.F90.

    call CNCIsoHarvestPftToColumn(num_soilc, filter_soilc, cnstate_vars, isotopeflux_vars, isocol_cf, isoveg_cf)

  end subroutine CarbonIsoFlux2h

  !-----------------------------------------------------------------------
  subroutine CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       isotopeflux_vars, isotopestate_vars, isotope, isocol_cs, isoveg_cs, isocol_cf, isoveg_cf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic fluxes for fire mortality
    use tracer_varcon, only : is_active_betr_bgc  
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(in)    :: carbonstate_vars
    type(carbonflux_type)  , intent(inout) :: isotopeflux_vars
    type(carbonstate_type) , intent(in)    :: isotopestate_vars
    character(len=*)       , intent(in)    :: isotope         ! 'c13' or 'c14'
    type(column_carbon_state),intent(in)   :: isocol_cs
    type(vegetation_carbon_state),intent(in)   :: isoveg_cs
    type(column_carbon_flux),intent(inout)    :: isocol_cf
    type(vegetation_carbon_flux),intent(inout):: isoveg_cf
    !
    ! !LOCAL VARIABLES:
    integer :: pi,pp,l,fc,cc,j
    !-----------------------------------------------------------------------

    associate(                                           &
         croot_prof =>   cnstate_vars%croot_prof_patch , & !  [real(r8) (:,:) ]  (1/m) profile of coarse roots                          
         stem_prof  =>   cnstate_vars%stem_prof_patch,   & !  [real(r8) (:,:) ]  (1/m) profile of stems                                 
         leaf_prof  =>   cnstate_vars%leaf_prof_patch    & !  [real(r8) (:,:) ]  (1/m) profile of leaves      
         )

      ! patch-level fire mortality fluxes

     call CarbonIsoFluxCalc(&
           isoveg_cf%m_leafc_to_fire              , veg_cf%m_leafc_to_fire, &
           isoveg_cs%leafc                       , veg_cs%leafc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_leafc_storage_to_fire      , veg_cf%m_leafc_storage_to_fire, &
           isoveg_cs%leafc_storage               , veg_cs%leafc_storage, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_leafc_xfer_to_fire         , veg_cf%m_leafc_xfer_to_fire, &
           isoveg_cs%leafc_xfer                  , veg_cs%leafc_xfer, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_frootc_to_fire             , veg_cf%m_frootc_to_fire, &
           isoveg_cs%frootc                      , veg_cs%frootc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_frootc_storage_to_fire     , veg_cf%m_frootc_storage_to_fire, &
           isoveg_cs%frootc_storage              , veg_cs%frootc_storage, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_frootc_xfer_to_fire        , veg_cf%m_frootc_xfer_to_fire, &
           isoveg_cs%frootc_xfer                 , veg_cs%frootc_xfer, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_livestemc_to_fire          , veg_cf%m_livestemc_to_fire, &
           isoveg_cs%livestemc                   , veg_cs%livestemc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_livestemc_storage_to_fire  , veg_cf%m_livestemc_storage_to_fire, &
           isoveg_cs%livestemc_storage           , veg_cs%livestemc_storage, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_livestemc_xfer_to_fire     , veg_cf%m_livestemc_xfer_to_fire, &
           isoveg_cs%livestemc_xfer              , veg_cs%livestemc_xfer, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadstemc_to_fire          , veg_cf%m_deadstemc_to_fire, &
           isoveg_cs%deadstemc                   , veg_cs%deadstemc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadstemc_to_litter_fire   , veg_cf%m_deadstemc_to_litter_fire, &
           isoveg_cs%deadstemc                   , veg_cs%deadstemc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadstemc_storage_to_fire  , veg_cf%m_deadstemc_storage_to_fire, &
           isoveg_cs%deadstemc_storage           , veg_cs%deadstemc_storage, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadstemc_xfer_to_fire     , veg_cf%m_deadstemc_xfer_to_fire, &
           isoveg_cs%deadstemc_xfer              , veg_cs%deadstemc_xfer, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_livecrootc_to_fire         , veg_cf%m_livecrootc_to_fire, &
           isoveg_cs%livecrootc                  , veg_cs%livecrootc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_livecrootc_storage_to_fire , veg_cf%m_livecrootc_storage_to_fire, &
           isoveg_cs%livecrootc_storage          , veg_cs%livecrootc_storage, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_livecrootc_xfer_to_fire    , veg_cf%m_livecrootc_xfer_to_fire, &
           isoveg_cs%livecrootc_xfer             , veg_cs%livecrootc_xfer, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadcrootc_to_fire         , veg_cf%m_deadcrootc_to_fire, &
           isoveg_cs%deadcrootc                  , veg_cs%deadcrootc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadcrootc_to_litter_fire  , veg_cf%m_deadcrootc_to_litter_fire, &
           isoveg_cs%deadcrootc                  , veg_cs%deadcrootc, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadcrootc_storage_to_fire , veg_cf%m_deadcrootc_storage_to_fire, &
           isoveg_cs%deadcrootc_storage          , veg_cs%deadcrootc_storage, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_deadcrootc_xfer_to_fire    , veg_cf%m_deadcrootc_xfer_to_fire, &
           isoveg_cs%deadcrootc_xfer             , veg_cs%deadcrootc_xfer, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_gresp_storage_to_fire      , veg_cf%m_gresp_storage_to_fire, &
           isoveg_cs%gresp_storage               , veg_cs%gresp_storage, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CarbonIsoFluxCalc(&
           isoveg_cf%m_gresp_xfer_to_fire         , veg_cf%m_gresp_xfer_to_fire, &
           isoveg_cs%gresp_xfer                  , veg_cs%gresp_xfer, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

     call CarbonIsoFluxCalc(&
           isoveg_cf%m_cpool_to_fire              , veg_cf%m_cpool_to_fire, &
           isoveg_cs%cpool                       , veg_cs%cpool, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)   

     call CarbonIsoFluxCalc(&
           isoveg_cf%m_cpool_to_litter_fire       , veg_cf%m_cpool_to_litter_fire, &
           isoveg_cs%cpool                       , veg_cs%cpool, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope) 


      if (.not. is_active_betr_bgc) then

         ! calculate the column-level flux of deadstem and deadcrootc to cwdc as the result of fire mortality.

         do pi = 1,max_patch_per_col
            do fc = 1,num_soilc
               cc = filter_soilc(fc)
               if ( pi <=  col_pp%npfts(cc) ) then
                  pp = col_pp%pfti(cc) + pi - 1
                  if (veg_pp%active(pp)) then
                     do j = 1, nlevdecomp
                        isocol_cf%fire_mortality_c_to_cwdc(cc,j) = &
                             isocol_cf%fire_mortality_c_to_cwdc(cc,j) + &
                             isoveg_cf%m_deadstemc_to_litter_fire(pp) * veg_pp%wtcol(pp) * stem_prof(pp,j)
                        isocol_cf%fire_mortality_c_to_cwdc(cc,j) = &
                             isocol_cf%fire_mortality_c_to_cwdc(cc,j) + &
                             isoveg_cf%m_deadcrootc_to_litter_fire(pp) * veg_pp%wtcol(pp) * croot_prof(pp,j)
                        isocol_cf%fire_mortality_c_to_cwdc(cc,j) = &
                             isocol_cf%fire_mortality_c_to_cwdc(cc,j) + &
                             isoveg_cf%m_cpool_to_litter_fire(pp) * veg_pp%wtcol(pp) * leaf_prof(pp,j)

                     end do
                  end if
               end if
            end do
         end do


         do fc = 1,num_soilc
            cc = filter_soilc(fc)
            do j = 1, nlevdecomp
               do l = 1, ndecomp_pools
                  if ( col_cs%decomp_cpools_vr(cc,j,l) /= 0._r8) then
                     isocol_cf%m_decomp_cpools_to_fire_vr(cc,j,l)  =  &
                          col_cf%m_decomp_cpools_to_fire_vr(cc,j,l) * &
                          (isocol_cs%decomp_cpools_vr(cc,j,l) / col_cs%decomp_cpools_vr(cc,j,l)) * 1._r8
                  else
                     isocol_cf%m_decomp_cpools_to_fire_vr(cc,j,l) = 0._r8
                  end if
               end do
            end do
         end do
    endif
    end associate
  end subroutine CarbonIsoFlux3

  !-----------------------------------------------------------------------
  subroutine CNCIsoLitterToColumn (num_soilc, filter_soilc, &
       cnstate_vars, carbonflux_vars, col_cf, veg_cf)
    !
    ! !DESCRIPTION:
    ! called at the end of cn_phenology to gather all patch-level litterfall fluxes
    ! to the column level and assign them to the three litter pools
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars
     type(column_carbon_flux), intent(inout) :: col_cf
     type(vegetation_carbon_flux), intent(inout) :: veg_cf
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j
    !-----------------------------------------------------------------------

    associate(                                                                           & 
         ivt                       =>    veg_pp%itype                                     , & ! Input:  [integer  (:)   ]  pft vegetation type                                
         wtcol                     =>    veg_pp%wtcol                                     , & ! Input:  [real(r8) (:)   ]  weight (relative to column) for this pft (0-1)    
         
         lf_flab                   =>    veg_vp%lf_flab                            , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
         lf_fcel                   =>    veg_vp%lf_fcel                            , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
         lf_flig                   =>    veg_vp%lf_flig                            , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
         fr_flab                   =>    veg_vp%fr_flab                            , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
         fr_fcel                   =>    veg_vp%fr_fcel                            , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
         fr_flig                   =>    veg_vp%fr_flig                            , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  

         leaf_prof                 =>    cnstate_vars%leaf_prof_patch                  , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                =>    cnstate_vars%froot_prof_patch                 , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         
         leafc_to_litter           =>    veg_cf%leafc_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
         frootc_to_litter          =>    veg_cf%frootc_to_litter        , & ! Input:  [real(r8) (:)   ]                                                    
         phenology_c_to_litr_met_c =>    col_cf%phenology_c_to_litr_met_c , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
         phenology_c_to_litr_cel_c =>    col_cf%phenology_c_to_litr_cel_c , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
         phenology_c_to_litr_lig_c =>    col_cf%phenology_c_to_litr_lig_c   & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)
         )

      do j = 1, nlevdecomp
         do pi = 1,max_patch_per_col
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if ( pi <=  col_pp%npfts(c) ) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then
                     ! leaf litter carbon fluxes
                     phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                          + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                          + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                          + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! fine root litter carbon fluxes
                     phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                          + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                          + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                          + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  end if
               end if

            end do
         end do

      end do

    end associate
   end subroutine CNCIsoLitterToColumn

   !-----------------------------------------------------------------------
   subroutine CNCIsoGapPftToColumn (num_soilc, filter_soilc, &
        cnstate_vars, carbonflux_vars, col_cf, veg_cf)
     !
     ! !DESCRIPTION:
     ! gather all patch-level gap mortality fluxes
     ! to the column level and assign them to the three litter pools (+ cwd pool)
     !
     ! !ARGUMENTS:
     integer               , intent(in)    :: num_soilc         ! number of soil columns in filter
     integer               , intent(in)    :: filter_soilc(:)   ! soil column filter
     type(cnstate_type)    , intent(in)    :: cnstate_vars
     type(carbonflux_type) , intent(inout) :: carbonflux_vars
     type(column_carbon_flux), intent(inout) :: col_cf
     type(vegetation_carbon_flux), intent(inout) :: veg_cf
     !
     ! !LOCAL VARIABLES:
     integer :: fc,c,pi,p,j               ! indices
     !-----------------------------------------------------------------------

     associate(                                                                                       & 
          ivt                            =>    veg_pp%itype                                            , & ! Input:  [integer  (:)   ]  pft vegetation type                                
          wtcol                          =>    veg_pp%wtcol                                            , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
          
          lf_flab                        =>    veg_vp%lf_flab                                   , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
          lf_fcel                        =>    veg_vp%lf_fcel                                   , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
          lf_flig                        =>    veg_vp%lf_flig                                   , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
          fr_flab                        =>    veg_vp%fr_flab                                   , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
          fr_fcel                        =>    veg_vp%fr_fcel                                   , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
          fr_flig                        =>    veg_vp%fr_flig                                   , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  

          leaf_prof                      =>    cnstate_vars%leaf_prof_patch                         , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
          froot_prof                     =>    cnstate_vars%froot_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
          croot_prof                     =>    cnstate_vars%croot_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
          stem_prof                      =>    cnstate_vars%stem_prof_patch                         , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
          
          m_leafc_to_litter              =>    veg_cf%m_leafc_to_litter              , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_to_litter             =>    veg_cf%m_frootc_to_litter             , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_to_litter          =>    veg_cf%m_livestemc_to_litter          , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_to_litter          =>    veg_cf%m_deadstemc_to_litter          , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_to_litter         =>    veg_cf%m_livecrootc_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_to_litter         =>    veg_cf%m_deadcrootc_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          m_leafc_storage_to_litter      =>    veg_cf%m_leafc_storage_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_storage_to_litter     =>    veg_cf%m_frootc_storage_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_storage_to_litter  =>    veg_cf%m_livestemc_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_storage_to_litter  =>    veg_cf%m_deadstemc_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_storage_to_litter =>    veg_cf%m_livecrootc_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_storage_to_litter =>    veg_cf%m_deadcrootc_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
          m_gresp_storage_to_litter      =>    veg_cf%m_gresp_storage_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
          m_leafc_xfer_to_litter         =>    veg_cf%m_leafc_xfer_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_xfer_to_litter        =>    veg_cf%m_frootc_xfer_to_litter        , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_xfer_to_litter     =>    veg_cf%m_livestemc_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_xfer_to_litter     =>    veg_cf%m_deadstemc_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_xfer_to_litter    =>    veg_cf%m_livecrootc_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_xfer_to_litter    =>    veg_cf%m_deadcrootc_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
          m_gresp_xfer_to_litter         =>    veg_cf%m_gresp_xfer_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          m_cpool_to_litter              =>    veg_cf%m_cpool_to_litter              , & ! Input:  [real(r8) (:)   ]  
          gap_mortality_c_to_litr_met_c  =>    col_cf%gap_mortality_c_to_litr_met_c    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
          gap_mortality_c_to_litr_cel_c  =>    col_cf%gap_mortality_c_to_litr_cel_c    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
          gap_mortality_c_to_litr_lig_c  =>    col_cf%gap_mortality_c_to_litr_lig_c    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
          gap_mortality_c_to_cwdc        =>    col_cf%gap_mortality_c_to_cwdc            & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
          )
          
       do j = 1, nlevdecomp
          do pi = 1,maxpatch_pft
             do fc = 1,num_soilc
                c = filter_soilc(fc)

                if (pi <=  col_pp%npfts(c)) then
                   p = col_pp%pfti(c) + pi - 1

                   if (veg_pp%active(p)) then

                      ! leaf gap mortality carbon fluxes
                      gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                           m_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                      gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                           m_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                      gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                           m_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                      ! fine root gap mortality carbon fluxes
                      gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                           m_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                      gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                           m_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                      gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                           m_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                      ! wood gap mortality carbon fluxes
                      gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                           m_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                      gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                           m_deadstemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                      gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                           m_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                      gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                           m_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)

                      ! storage gap mortality carbon fluxes
                      gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                           m_leafc_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                           m_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                           m_livestemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                           m_deadstemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                           m_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                           m_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                           m_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                           m_cpool_to_litter(p)      * wtcol(p) * leaf_prof(p,j)


                      ! transfer gap mortality carbon fluxes
                      gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                           m_leafc_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                           m_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                           m_livestemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                           m_deadstemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                           m_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                           m_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                      gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                           m_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)

                   end if
                end if

             end do

          end do
       end do

     end associate

   end subroutine CNCIsoGapPftToColumn

   !-----------------------------------------------------------------------
   subroutine CNCIsoHarvestPftToColumn (&
        num_soilc, filter_soilc, &
        cnstate_vars, carbonflux_vars, col_cf, veg_cf)
     !
     ! !DESCRIPTION:
     ! gather all patch-level harvest mortality fluxes
     ! to the column level and assign them to the litter, cwd, and wood product pools
     !
     ! !ARGUMENTS:
     integer               , intent(in)    :: num_soilc         ! number of soil columns in filter
     integer               , intent(in)    :: filter_soilc(:)   ! soil column filter
     type(cnstate_type)    , intent(in)    :: cnstate_vars
     type(carbonflux_type) , intent(inout) :: carbonflux_vars
     type(column_carbon_flux), intent(inout) :: col_cf
     type(vegetation_carbon_flux), intent(inout) :: veg_cf
     !
     ! !LOCAL VARIABLES:
     integer :: fc,c,pi,p,j               ! indices
     !-----------------------------------------------------------------------

     associate(                                                                                           & 
          ivt                              =>    veg_pp%itype                                              , & ! Input:  [integer  (:)   ]  pft vegetation type                                
          wtcol                            =>    veg_pp%wtcol                                              , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
          
          lf_flab                          =>    veg_vp%lf_flab                                     , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
          lf_fcel                          =>    veg_vp%lf_fcel                                     , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
          lf_flig                          =>    veg_vp%lf_flig                                     , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
          fr_flab                          =>    veg_vp%fr_flab                                     , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
          fr_fcel                          =>    veg_vp%fr_fcel                                     , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
          fr_flig                          =>    veg_vp%fr_flig                                     , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  
          
          leaf_prof                        =>    cnstate_vars%leaf_prof_patch                           , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
          froot_prof                       =>    cnstate_vars%froot_prof_patch                          , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
          croot_prof                       =>    cnstate_vars%croot_prof_patch                          , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
          stem_prof                        =>    cnstate_vars%stem_prof_patch                           , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
          
          hrv_leafc_to_litter              =>    veg_cf%hrv_leafc_to_litter              , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_to_litter             =>    veg_cf%hrv_frootc_to_litter             , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_to_litter          =>    veg_cf%hrv_livestemc_to_litter          , & ! Input:  [real(r8) (:)   ]                                                    
          phrv_deadstemc_to_prod10c        =>    veg_cf%hrv_deadstemc_to_prod10c         , & ! Input:  [real(r8) (:)   ]                                                    
          phrv_deadstemc_to_prod100c       =>    veg_cf%hrv_deadstemc_to_prod100c        , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_to_litter         =>    veg_cf%hrv_livecrootc_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_to_litter         =>    veg_cf%hrv_deadcrootc_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_leafc_storage_to_litter      =>    veg_cf%hrv_leafc_storage_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_storage_to_litter     =>    veg_cf%hrv_frootc_storage_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_storage_to_litter  =>    veg_cf%hrv_livestemc_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadstemc_storage_to_litter  =>    veg_cf%hrv_deadstemc_storage_to_litter  , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_storage_to_litter =>    veg_cf%hrv_livecrootc_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_storage_to_litter =>    veg_cf%hrv_deadcrootc_storage_to_litter , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_gresp_storage_to_litter      =>    veg_cf%hrv_gresp_storage_to_litter      , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_leafc_xfer_to_litter         =>    veg_cf%hrv_leafc_xfer_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_xfer_to_litter        =>    veg_cf%hrv_frootc_xfer_to_litter        , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_xfer_to_litter     =>    veg_cf%hrv_livestemc_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadstemc_xfer_to_litter     =>    veg_cf%hrv_deadstemc_xfer_to_litter     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_xfer_to_litter    =>    veg_cf%hrv_livecrootc_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_xfer_to_litter    =>    veg_cf%hrv_deadcrootc_xfer_to_litter    , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_gresp_xfer_to_litter         =>    veg_cf%hrv_gresp_xfer_to_litter         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_cpool_to_litter              =>    veg_cf%hrv_cpool_to_litter              , & ! Input:  [real(r8) (:)   ]      

          chrv_deadstemc_to_prod10c        =>    col_cf%hrv_deadstemc_to_prod10c           , & ! InOut:  [real(r8) (:)   ]                                                    
          chrv_deadstemc_to_prod100c       =>    col_cf%hrv_deadstemc_to_prod100c          , & ! InOut:  [real(r8) (:)   ]                                                    
          harvest_c_to_litr_met_c          =>    col_cf%harvest_c_to_litr_met_c            , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
          harvest_c_to_litr_cel_c          =>    col_cf%harvest_c_to_litr_cel_c            , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
          harvest_c_to_litr_lig_c          =>    col_cf%harvest_c_to_litr_lig_c            , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
          harvest_c_to_cwdc                =>    col_cf%harvest_c_to_cwdc                    & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
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
                   chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                        phrv_deadstemc_to_prod10c(p)  * wtcol(p)
                   chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                        phrv_deadstemc_to_prod100c(p)  * wtcol(p)
                end if
             end if
          end do
       end do

     end associate 

   end subroutine CNCIsoHarvestPftToColumn

   !-----------------------------------------------------------------------
   subroutine CarbonIsoFluxCalc(&
        ciso_flux, ctot_flux, &
        ciso_state, ctot_state, &
        num, filter, frax_c13, diag, isotope)
     !
     ! !DESCRIPTION:
     ! On the radiation time step, set the carbon isotopic flux
     ! variables (except for gap-phase mortality and fire fluxes)
     !
     ! !ARGUMENTS:
     real(r8)         , intent(inout), pointer :: ciso_flux(:)  ! isoC flux
     real(r8)         , intent(in)   , pointer :: ctot_flux(:)  ! totC flux
     real(r8)         , intent(in)   , pointer :: ciso_state(:) ! isoC state, upstream pool
     real(r8)         , intent(in)   , pointer :: ctot_state(:) ! totC state, upstream pool
     real(r8)         , intent(in)             :: frax_c13      ! fractionation factor (1 = no fractionation) for C13
     integer          , intent(in)             :: num           ! number of filter members
     integer          , intent(in)             :: filter(:)     ! filter indices
     integer          , intent(in)             :: diag          ! 0=no diagnostics, 1=print diagnostics
     character(len=*) , intent(in)             :: isotope       ! 'c13' or 'c14'
     !
     ! ! LOCAL VARIABLES:
     integer  :: i,f     ! indices
     real(r8) :: temp
     real(r8) :: frax
     !-----------------------------------------------------------------------

     ! if C14, double the fractionation
     select case (isotope)
     case ('c14')
        frax = 1._r8 + (1._r8 - frax_c13) * 2._r8
     case ('c13')
        frax = frax_c13
     case default
        call endrun(msg='CarbonIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
     end select

     ! loop over the supplied filter
     do f = 1,num
        i = filter(f)
        if (ctot_state(i) /= 0._r8) then
           ciso_flux(i) = ctot_flux(i) * (ciso_state(i)/ctot_state(i)) * frax
        else
           ciso_flux(i) = 0._r8
        end if

        if (diag == 1) then
           ! put diagnostic print statements here for isoC flux calculations
        end if
     end do

   end subroutine CarbonIsoFluxCalc

end module CarbonIsoFluxMod
