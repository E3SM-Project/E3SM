module CNCIsoFluxMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon isotopic flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  use clm_varpar                         , only : max_patch_per_col, maxpatch_pft
  use abortutils                         , only : endrun
  use pftconMod                          , only : pftcon
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use ColumnType                         , only : col                
  use PatchType                          , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: CIsoFlux1
  public  :: CIsoFlux2
  public  :: CIsoFlux2h
  public  :: CIsoFlux3
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNCIsoLitterToColumn
  private :: CNCIsoGapPftToColumn
  private :: CNCIsoHarvestPftToColumn
  private :: CIsoFluxCalc
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp,         &
       soilbiogeochem_state_inst,                                                &
       soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,         &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                            &
       iso_soilbiogeochem_carbonflux_inst,  iso_soilbiogeochem_carbonstate_inst, &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst,                    &
       isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic flux
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                               , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type)       , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonflux_type)  , intent(in)    :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)  , intent(inout) :: iso_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: iso_soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: iso_cnveg_carbonstate_inst
    character(len=*)                      , intent(in)    :: isotope         ! 'c13' or 'c14'
    !
    ! !LOCAL VARIABLES:
    integer :: fp,pi,l,fc,cc,j
    integer :: cdp 
    !-----------------------------------------------------------------------

    associate(                                                            &
         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool , & 
         soilbiogeochem_cs     => soilbiogeochem_carbonstate_inst       , &
         soilbiogeochem_cf     => soilbiogeochem_carbonflux_inst        , &
         cnveg_cf              => cnveg_carbonflux_inst                 , &
         cnveg_cs              => cnveg_carbonstate_inst                , &
         iso_cnveg_cf          => iso_cnveg_carbonflux_inst             , &
         iso_cnveg_cs          => iso_cnveg_carbonstate_inst            , &
         iso_soilbiogeochem_cs => iso_soilbiogeochem_carbonstate_inst   , &
         iso_soilbiogeochem_cf => iso_soilbiogeochem_carbonflux_inst      &
         )

      ! patch-level non-mortality fluxes
      
      ! Note: if the variables which are arguments to CIsoFluxCalc are ever changed to NOT be
      ! pointers, then the CIsoFluxCalc routine will need to be changed to declare the bounds
      ! of each argument, these bounds will need to be passed in, and - importantly for
      ! threading to work properly - the subroutine calls will need to be changed so that
      ! instead of 'call CIsoFluxCalc(foo, ...)' we have 'call CIsoFluxCalc(foo(begp:endp), ...)'.
      
      call CIsoFluxCalc(&
           iso_cnveg_cf%leafc_xfer_to_leafc_patch           , cnveg_cf%leafc_xfer_to_leafc_patch, &
           iso_cnveg_cs%leafc_xfer_patch                    , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%frootc_xfer_to_frootc_patch         , cnveg_cf%frootc_xfer_to_frootc_patch, &
           iso_cnveg_cs%frootc_xfer_patch                   , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestemc_xfer_to_livestemc_patch   , cnveg_cf%livestemc_xfer_to_livestemc_patch, &
           iso_cnveg_cs%livestemc_xfer_patch                , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadstemc_xfer_to_deadstemc_patch   , cnveg_cf%deadstemc_xfer_to_deadstemc_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch                , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecrootc_xfer_to_livecrootc_patch , cnveg_cf%livecrootc_xfer_to_livecrootc_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch               , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadcrootc_xfer_to_deadcrootc_patch , cnveg_cf%deadcrootc_xfer_to_deadcrootc_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch               , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leafc_to_litter_patch               , cnveg_cf%leafc_to_litter_patch, &
           iso_cnveg_cs%leafc_patch                         , cnveg_cs%leafc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%frootc_to_litter_patch              , cnveg_cf%frootc_to_litter_patch, &
           iso_cnveg_cs%frootc_patch                        , cnveg_cs%frootc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestemc_to_deadstemc_patch        , cnveg_cf%livestemc_to_deadstemc_patch, &
           iso_cnveg_cs%livestemc_patch                     , cnveg_cs%livestemc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecrootc_to_deadcrootc_patch      , cnveg_cf%livecrootc_to_deadcrootc_patch, &
           iso_cnveg_cs%livecrootc_patch                    , cnveg_cs%livecrootc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leaf_curmr_patch                    , cnveg_cf%leaf_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%froot_curmr_patch                   , cnveg_cf%froot_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestem_curmr_patch                , cnveg_cf%livestem_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecroot_curmr_patch               , cnveg_cf%livecroot_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leaf_xsmr_patch                     , cnveg_cf%leaf_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%froot_xsmr_patch                    , cnveg_cf%froot_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestem_xsmr_patch                 , cnveg_cf%livestem_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecroot_xsmr_patch                , cnveg_cf%livecroot_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_xsmrpool_patch             , cnveg_cf%cpool_to_xsmrpool_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_leafc_patch                , cnveg_cf%cpool_to_leafc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_leafc_storage_patch        , cnveg_cf%cpool_to_leafc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_frootc_patch               , cnveg_cf%cpool_to_frootc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_frootc_storage_patch       , cnveg_cf%cpool_to_frootc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livestemc_patch            , cnveg_cf%cpool_to_livestemc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livestemc_storage_patch    , cnveg_cf%cpool_to_livestemc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadstemc_patch            , cnveg_cf%cpool_to_deadstemc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadstemc_storage_patch    , cnveg_cf%cpool_to_deadstemc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livecrootc_patch           , cnveg_cf%cpool_to_livecrootc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livecrootc_storage_patch   , cnveg_cf%cpool_to_livecrootc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadcrootc_patch           , cnveg_cf%cpool_to_deadcrootc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadcrootc_storage_patch   , cnveg_cf%cpool_to_deadcrootc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_leaf_gr_patch                 , cnveg_cf%cpool_leaf_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_froot_gr_patch                , cnveg_cf%cpool_froot_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livestem_gr_patch             , cnveg_cf%cpool_livestem_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadstem_gr_patch             , cnveg_cf%cpool_deadstem_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livecroot_gr_patch            , cnveg_cf%cpool_livecroot_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadcroot_gr_patch            , cnveg_cf%cpool_deadcroot_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_leaf_storage_gr_patch         , cnveg_cf%cpool_leaf_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_froot_storage_gr_patch        , cnveg_cf%cpool_froot_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livestem_storage_gr_patch     , cnveg_cf%cpool_livestem_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadstem_storage_gr_patch     , cnveg_cf%cpool_deadstem_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livecroot_storage_gr_patch    , cnveg_cf%cpool_livecroot_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadcroot_storage_gr_patch    , cnveg_cf%cpool_deadcroot_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_gresp_storage_patch        , cnveg_cf%cpool_to_gresp_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_leaf_gr_patch              , cnveg_cf%transfer_leaf_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_froot_gr_patch             , cnveg_cf%transfer_froot_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_livestem_gr_patch          , cnveg_cf%transfer_livestem_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_deadstem_gr_patch          , cnveg_cf%transfer_deadstem_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_livecroot_gr_patch         , cnveg_cf%transfer_livecroot_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_deadcroot_gr_patch         , cnveg_cf%transfer_deadcroot_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leafc_storage_to_xfer_patch         , cnveg_cf%leafc_storage_to_xfer_patch, &
           iso_cnveg_cs%leafc_storage_patch                 , cnveg_cs%leafc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%frootc_storage_to_xfer_patch        , cnveg_cf%frootc_storage_to_xfer_patch, &
           iso_cnveg_cs%frootc_storage_patch                , cnveg_cs%frootc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestemc_storage_to_xfer_patch     , cnveg_cf%livestemc_storage_to_xfer_patch, &
           iso_cnveg_cs%livestemc_storage_patch             , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadstemc_storage_to_xfer_patch     , cnveg_cf%deadstemc_storage_to_xfer_patch, &
           iso_cnveg_cs%deadstemc_storage_patch             , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecrootc_storage_to_xfer_patch    , cnveg_cf%livecrootc_storage_to_xfer_patch, &
           iso_cnveg_cs%livecrootc_storage_patch            , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadcrootc_storage_to_xfer_patch    , cnveg_cf%deadcrootc_storage_to_xfer_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch            , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gresp_storage_to_xfer_patch         , cnveg_cf%gresp_storage_to_xfer_patch, &
           iso_cnveg_cs%gresp_storage_patch                 , cnveg_cs%gresp_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      ! call routine to shift patch-level litterfall fluxes to column, for isotopes
      ! the non-isotope version of this routine is called in CNPhenologyMod.F90
      ! For later clean-up, it would be possible to generalize this function to operate on a single 
      ! patch-to-column flux.

      call CNCIsoLitterToColumn(num_soilc, filter_soilc, soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)

      ! column-level non-mortality fluxes

      do fc = 1,num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) /= 0._r8) then
                  iso_soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l)  =  &
                      soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l) * &
                      (iso_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) &
                         / soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp)) * 1._r8
               else
                  iso_soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

      do fc = 1,num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) /= 0._r8) then
                  iso_soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l)  =  &
                      soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) * &
                      (iso_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) &
                      / soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp)) * 1._r8
               else
                  iso_soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

    end associate

  end subroutine CIsoFlux1

  !-----------------------------------------------------------------------
  subroutine CIsoFlux2(num_soilc, filter_soilc, num_soilp  , filter_soilp, &
       soilbiogeochem_state_inst, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst, isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic fluxes for gap mortality
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                         , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)     , intent(in)    :: cnveg_carbonflux_inst 
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: iso_cnveg_carbonstate_inst
    character(len=*)                , intent(in)    :: isotope         ! 'c13' or 'c14'
    !
    ! !LOCAL VARIABLES:
    integer :: fp,pi
    !-----------------------------------------------------------------------

    associate(                                               &
         cnveg_cf     => cnveg_carbonflux_inst     , &
         cnveg_cs     => cnveg_carbonstate_inst    , &
         iso_cnveg_cf => iso_cnveg_carbonflux_inst , &
         iso_cnveg_cs => iso_cnveg_carbonstate_inst        &
         )

      ! patch-level gap mortality fluxes
      
      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_to_litter_patch                          , cnveg_cf%m_leafc_to_litter_patch, &
           iso_cnveg_cs%leafc_patch                                      , cnveg_cs%leafc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_storage_to_litter_patch                  , cnveg_cf%m_leafc_storage_to_litter_patch, &
           iso_cnveg_cs%leafc_storage_patch                              , cnveg_cs%leafc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_xfer_to_litter_patch                     , cnveg_cf%m_leafc_xfer_to_litter_patch, &
           iso_cnveg_cs%leafc_xfer_patch                                 , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_to_litter_patch                         , cnveg_cf%m_frootc_to_litter_patch, &
           iso_cnveg_cs%frootc_patch                                     , cnveg_cs%frootc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_storage_to_litter_patch                 , cnveg_cf%m_frootc_storage_to_litter_patch, &
           iso_cnveg_cs%frootc_storage_patch                             , cnveg_cs%frootc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_xfer_to_litter_patch                    , cnveg_cf%m_frootc_xfer_to_litter_patch, &
           iso_cnveg_cs%frootc_xfer_patch                                , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_to_litter_patch                      , cnveg_cf%m_livestemc_to_litter_patch, &
           iso_cnveg_cs%livestemc_patch                                  , cnveg_cs%livestemc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_storage_to_litter_patch              , cnveg_cf%m_livestemc_storage_to_litter_patch, &
           iso_cnveg_cs%livestemc_storage_patch                          , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_xfer_to_litter_patch                 , cnveg_cf%m_livestemc_xfer_to_litter_patch, &
           iso_cnveg_cs%livestemc_xfer_patch                             , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_to_litter_patch                      , cnveg_cf%m_deadstemc_to_litter_patch, &
           iso_cnveg_cs%deadstemc_patch                                  , cnveg_cs%deadstemc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_storage_to_litter_patch              , cnveg_cf%m_deadstemc_storage_to_litter_patch, &
           iso_cnveg_cs%deadstemc_storage_patch                          , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_xfer_to_litter_patch                 , cnveg_cf%m_deadstemc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch                             , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_to_litter_patch                     , cnveg_cf%m_livecrootc_to_litter_patch, &
           iso_cnveg_cs%livecrootc_patch                                 , cnveg_cs%livecrootc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_storage_to_litter_patch             , cnveg_cf%m_livecrootc_storage_to_litter_patch, &
           iso_cnveg_cs%livecrootc_storage_patch                         , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_xfer_to_litter_patch                , cnveg_cf%m_livecrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch                            , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_to_litter_patch                     , cnveg_cf%m_deadcrootc_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_patch                                 , cnveg_cs%deadcrootc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_storage_to_litter_patch             , cnveg_cf%m_deadcrootc_storage_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch                         , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_xfer_to_litter_patch                , cnveg_cf%m_deadcrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch                            , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_storage_to_litter_patch                  , cnveg_cf%m_gresp_storage_to_litter_patch, &
           iso_cnveg_cs%gresp_storage_patch                              , cnveg_cs%gresp_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_xfer_to_litter_patch                     , cnveg_cf%m_gresp_xfer_to_litter_patch, &
           iso_cnveg_cs%gresp_xfer_patch                                 , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      ! call routine to shift patch-level gap mortality fluxes to column , for isotopes
      ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

      call CNCIsoGapPftToColumn(num_soilc, filter_soilc, soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)

    end associate

  end subroutine CIsoFlux2

  !-----------------------------------------------------------------------
  subroutine CIsoFlux2h(num_soilc , filter_soilc, num_soilp  , filter_soilp, &
       soilbiogeochem_state_inst,                                            &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                        &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst, isotope) 
    !
    ! !DESCRIPTION:
    ! set the carbon isotopic fluxes for harvest mortality
    !
    ! !ARGUMENTS:
    integer                           , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                           , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                           , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                           , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type)   , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)       , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)      , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)       , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)      , intent(in)    :: iso_cnveg_carbonstate_inst
    character(len=*)                  , intent(in)    :: isotope         ! 'c13' or 'c14'
    !-----------------------------------------------------------------------

    associate(                                               &
         cnveg_cf     => cnveg_carbonflux_inst           , &
         cnveg_cs     => cnveg_carbonstate_inst          , &
         iso_cnveg_cf => iso_cnveg_carbonflux_inst       , &
         iso_cnveg_cs => iso_cnveg_carbonstate_inst        &
         )

      ! patch-level gap mortality fluxes
   
      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_leafc_to_litter_patch              , cnveg_cf%hrv_leafc_to_litter_patch, &
           iso_cnveg_cs%leafc_patch                            , cnveg_cs%leafc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_leafc_storage_to_litter_patch      , cnveg_cf%hrv_leafc_storage_to_litter_patch, &
           iso_cnveg_cs%leafc_storage_patch                    , cnveg_cs%leafc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_leafc_xfer_to_litter_patch         , cnveg_cf%hrv_leafc_xfer_to_litter_patch, &
           iso_cnveg_cs%leafc_xfer_patch                       , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_frootc_to_litter_patch             , cnveg_cf%hrv_frootc_to_litter_patch, &
           iso_cnveg_cs%frootc_patch                           , cnveg_cs%frootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_frootc_storage_to_litter_patch     , cnveg_cf%hrv_frootc_storage_to_litter_patch, &
           iso_cnveg_cs%frootc_storage_patch                   , cnveg_cs%frootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_frootc_xfer_to_litter_patch        , cnveg_cf%hrv_frootc_xfer_to_litter_patch, &
           iso_cnveg_cs%frootc_xfer_patch                      , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livestemc_to_litter_patch          , cnveg_cf%hrv_livestemc_to_litter_patch, &
           iso_cnveg_cs%livestemc_patch                        , cnveg_cs%livestemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livestemc_storage_to_litter_patch  , cnveg_cf%hrv_livestemc_storage_to_litter_patch, &
           iso_cnveg_cs%livestemc_storage_patch                , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livestemc_xfer_to_litter_patch     , cnveg_cf%hrv_livestemc_xfer_to_litter_patch, &
           iso_cnveg_cs%livestemc_xfer_patch                   , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadstemc_to_prod10c_patch         , cnveg_cf%hrv_deadstemc_to_prod10c_patch, &
           iso_cnveg_cs%deadstemc_patch                        , cnveg_cs%deadstemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadstemc_to_prod100c_patch        , cnveg_cf%hrv_deadstemc_to_prod100c_patch, &
           iso_cnveg_cs%deadstemc_patch                        , cnveg_cs%deadstemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadstemc_storage_to_litter_patch  , cnveg_cf%hrv_deadstemc_storage_to_litter_patch, &
           iso_cnveg_cs%deadstemc_storage_patch                , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadstemc_xfer_to_litter_patch     , cnveg_cf%hrv_deadstemc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch                   , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livecrootc_to_litter_patch         , cnveg_cf%hrv_livecrootc_to_litter_patch, &
           iso_cnveg_cs%livecrootc_patch                       , cnveg_cs%livecrootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livecrootc_storage_to_litter_patch , cnveg_cf%hrv_livecrootc_storage_to_litter_patch, &
           iso_cnveg_cs%livecrootc_storage_patch               , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livecrootc_xfer_to_litter_patch    , cnveg_cf%hrv_livecrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch                  , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadcrootc_to_litter_patch         , cnveg_cf%hrv_deadcrootc_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_patch                       , cnveg_cs%deadcrootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadcrootc_storage_to_litter_patch , cnveg_cf%hrv_deadcrootc_storage_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch               , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadcrootc_xfer_to_litter_patch    , cnveg_cf%hrv_deadcrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch                  , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_gresp_storage_to_litter_patch      , cnveg_cf%hrv_gresp_storage_to_litter_patch, &
           iso_cnveg_cs%gresp_storage_patch                    , cnveg_cs%gresp_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_gresp_xfer_to_litter_patch         , cnveg_cf%hrv_gresp_xfer_to_litter_patch, &
           iso_cnveg_cs%gresp_xfer_patch                       , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_xsmrpool_to_atm_patch              , cnveg_cf%hrv_xsmrpool_to_atm_patch, &
           cnveg_cs%totvegc_patch                          , cnveg_cs%totvegc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      ! call routine to shift patch-level gap mortality fluxes to column, 
      ! for isotopes the non-isotope version of this routine is in CNGapMortalityMod.F90.

      call CNCIsoHarvestPftToColumn(num_soilc, filter_soilc, soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)

    end associate

  end subroutine CIsoFlux2h

  !-----------------------------------------------------------------------
  subroutine CIsoFlux3(num_soilc , filter_soilc, num_soilp  , filter_soilp, &
       soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst,               &
       iso_soilbiogeochem_carbonstate_inst, isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic fluxes for fire mortality
    !
    ! !ARGUMENTS:
    integer                               , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type)       , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: iso_cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: iso_soilbiogeochem_carbonstate_inst
    character(len=*)                      , intent(in)    :: isotope         ! 'c13' or 'c14'
    !
    ! !LOCAL VARIABLES:
    integer :: pi,pp,l,fc,cc,j
    !-----------------------------------------------------------------------

    associate(                                                                 &
         croot_prof            => soilbiogeochem_state_inst%croot_prof_patch , & ! Input: [real(r8) (:,:) ]  (1/m) profile of coarse roots                          
         stem_prof             => soilbiogeochem_state_inst%stem_prof_patch  , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                                 
         soilbiogeochem_cs     => soilbiogeochem_carbonstate_inst            , &
         cnveg_cf              => cnveg_carbonflux_inst                      , &
         cnveg_cs              => cnveg_carbonstate_inst                     , &
         iso_cnveg_cf          => iso_cnveg_carbonflux_inst                  , &
         iso_cnveg_cs          => iso_cnveg_carbonstate_inst                 , &
         iso_soilbiogeochem_cs => iso_soilbiogeochem_carbonstate_inst          &
         )

      ! patch-level fire mortality fluxes

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_to_fire_patch              , cnveg_cf%m_leafc_to_fire_patch, &
           iso_cnveg_cs%leafc_patch                        , cnveg_cs%leafc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_storage_to_fire_patch      , cnveg_cf%m_leafc_storage_to_fire_patch, &
           iso_cnveg_cs%leafc_storage_patch                , cnveg_cs%leafc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_xfer_to_fire_patch         , cnveg_cf%m_leafc_xfer_to_fire_patch, &
           iso_cnveg_cs%leafc_xfer_patch                   , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_to_fire_patch             , cnveg_cf%m_frootc_to_fire_patch, &
           iso_cnveg_cs%frootc_patch                       , cnveg_cs%frootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_storage_to_fire_patch     , cnveg_cf%m_frootc_storage_to_fire_patch, &
           iso_cnveg_cs%frootc_storage_patch               , cnveg_cs%frootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_xfer_to_fire_patch        , cnveg_cf%m_frootc_xfer_to_fire_patch, &
           iso_cnveg_cs%frootc_xfer_patch                  , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_to_fire_patch          , cnveg_cf%m_livestemc_to_fire_patch, &
           iso_cnveg_cs%livestemc_patch                    , cnveg_cs%livestemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_storage_to_fire_patch  , cnveg_cf%m_livestemc_storage_to_fire_patch, &
           iso_cnveg_cs%livestemc_storage_patch            , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_xfer_to_fire_patch     , cnveg_cf%m_livestemc_xfer_to_fire_patch, &
           iso_cnveg_cs%livestemc_xfer_patch               , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_to_fire_patch          , cnveg_cf%m_deadstemc_to_fire_patch, &
           iso_cnveg_cs%deadstemc_patch                    , cnveg_cs%deadstemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_to_litter_fire_patch   , cnveg_cf%m_deadstemc_to_litter_fire_patch, &
           iso_cnveg_cs%deadstemc_patch                    , cnveg_cs%deadstemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_storage_to_fire_patch  , cnveg_cf%m_deadstemc_storage_to_fire_patch, &
           iso_cnveg_cs%deadstemc_storage_patch            , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_xfer_to_fire_patch     , cnveg_cf%m_deadstemc_xfer_to_fire_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch               , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_to_fire_patch         , cnveg_cf%m_livecrootc_to_fire_patch, &
           iso_cnveg_cs%livecrootc_patch                   , cnveg_cs%livecrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_storage_to_fire_patch , cnveg_cf%m_livecrootc_storage_to_fire_patch, &
           iso_cnveg_cs%livecrootc_storage_patch           , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_xfer_to_fire_patch    , cnveg_cf%m_livecrootc_xfer_to_fire_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch              , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_to_fire_patch         , cnveg_cf%m_deadcrootc_to_fire_patch, &
           iso_cnveg_cs%deadcrootc_patch                   , cnveg_cs%deadcrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_to_litter_fire_patch  , cnveg_cf%m_deadcrootc_to_litter_fire_patch, &
           iso_cnveg_cs%deadcrootc_patch                   , cnveg_cs%deadcrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_storage_to_fire_patch , cnveg_cf%m_deadcrootc_storage_to_fire_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch           , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_xfer_to_fire_patch    , cnveg_cf%m_deadcrootc_xfer_to_fire_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch              , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_storage_to_fire_patch      , cnveg_cf%m_gresp_storage_to_fire_patch, &
           iso_cnveg_cs%gresp_storage_patch                , cnveg_cs%gresp_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_xfer_to_fire_patch         , cnveg_cf%m_gresp_xfer_to_fire_patch, &
           iso_cnveg_cs%gresp_xfer_patch                   , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      ! calculate the column-level flux of deadstem and deadcrootc to cwdc as the result of fire mortality.
      do pi = 1,max_patch_per_col
         do fc = 1,num_soilc
            cc = filter_soilc(fc)
            if ( pi <=  col%npatches(cc) ) then
               pp = col%patchi(cc) + pi - 1
               if (patch%active(pp)) then
                  do j = 1, nlevdecomp
                     iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) = &
                          iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) + &
                          iso_cnveg_cf%m_deadstemc_to_litter_fire_patch(pp) * patch%wtcol(pp) * stem_prof(pp,j)
                     iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) = &
                          iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) + &
                          iso_cnveg_cf%m_deadcrootc_to_litter_fire_patch(pp) * patch%wtcol(pp) * croot_prof(pp,j)
                  end do
               end if
            end if
         end do
      end do


      do fc = 1,num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_pools
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,l) /= 0._r8) then
                  iso_cnveg_cf%m_decomp_cpools_to_fire_vr_col(cc,j,l)  =  &
                      cnveg_cf%m_decomp_cpools_to_fire_vr_col(cc,j,l) * &
                      (iso_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,l) / &
                           soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,l)) * 1._r8
               else
                  iso_cnveg_cf%m_decomp_cpools_to_fire_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

    end associate

  end subroutine CIsoFlux3

  !-----------------------------------------------------------------------
  subroutine CNCIsoLitterToColumn (num_soilc, filter_soilc, &
       soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! called at the end of cn_phenology to gather all patch-level litterfall fluxes
    ! to the column level and assign them to the three litter pools
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j
    !-----------------------------------------------------------------------

    associate(                                                                                     & 
         ivt                       => patch%itype                                             , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                     => patch%wtcol                                             , & ! Input:  [real(r8) (:)   ]  weight (relative to column) for this patch (0-1)    
         
         lf_flab                   => pftcon%lf_flab                                          , & ! Input:  leaf litter labile fraction                       
         lf_fcel                   => pftcon%lf_fcel                                          , & ! Input:  leaf litter cellulose fraction                    
         lf_flig                   => pftcon%lf_flig                                          , & ! Input:  leaf litter lignin fraction                       
         fr_flab                   => pftcon%fr_flab                                          , & ! Input:  fine root litter labile fraction                  
         fr_fcel                   => pftcon%fr_fcel                                          , & ! Input:  fine root litter cellulose fraction               
         fr_flig                   => pftcon%fr_flig                                          , & ! Input:  fine root litter lignin fraction                  

         leaf_prof                 => soilbiogeochem_state_inst%leaf_prof_patch               , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                => soilbiogeochem_state_inst%froot_prof_patch              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         
         leafc_to_litter           => iso_cnveg_carbonflux_inst%leafc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         frootc_to_litter          => iso_cnveg_carbonflux_inst%frootc_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         phenology_c_to_litr_met_c => iso_cnveg_carbonflux_inst%phenology_c_to_litr_met_c_col , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
         phenology_c_to_litr_cel_c => iso_cnveg_carbonflux_inst%phenology_c_to_litr_cel_c_col , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
         phenology_c_to_litr_lig_c => iso_cnveg_carbonflux_inst%phenology_c_to_litr_lig_c_col   & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)
         )

      do j = 1, nlevdecomp
         do pi = 1,max_patch_per_col
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if ( pi <=  col%npatches(c) ) then
                  p = col%patchi(c) + pi - 1
                  if (patch%active(p)) then
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
        soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)
     !
     ! !DESCRIPTION:
     ! gather all patch-level gap mortality fluxes
     ! to the column level and assign them to the three litter pools (+ cwd pool)
     !
     ! !ARGUMENTS:
     integer                         , intent(in)    :: num_soilc         ! number of soil columns in filter
     integer                         , intent(in)    :: filter_soilc(:)   ! soil column filter
     type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
     type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fc,c,pi,p,j               ! indices
     !-----------------------------------------------------------------------

     associate(                                                                                             & 
          ivt                            => patch%itype                                                  , & ! Input:  [integer  (:)   ]  patch vegetation type                                
          wtcol                          => patch%wtcol                                                  , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
          
          lf_flab                        => pftcon%lf_flab                                             , & ! Input:  leaf litter labile fraction                       
          lf_fcel                        => pftcon%lf_fcel                                             , & ! Input:  leaf litter cellulose fraction                    
          lf_flig                        => pftcon%lf_flig                                             , & ! Input:  leaf litter lignin fraction                       
          fr_flab                        => pftcon%fr_flab                                             , & ! Input:  fine root litter labile fraction                  
          fr_fcel                        => pftcon%fr_fcel                                             , & ! Input:  fine root litter cellulose fraction               
          fr_flig                        => pftcon%fr_flig                                             , & ! Input:  fine root litter lignin fraction                  

          leaf_prof                      => soilbiogeochem_state_inst%leaf_prof_patch                  , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
          froot_prof                     => soilbiogeochem_state_inst%froot_prof_patch                 , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
          croot_prof                     => soilbiogeochem_state_inst%croot_prof_patch                 , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
          stem_prof                      => soilbiogeochem_state_inst%stem_prof_patch                  , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
          
          m_leafc_to_litter              => iso_cnveg_carbonflux_inst%m_leafc_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_to_litter             => iso_cnveg_carbonflux_inst%m_frootc_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_to_litter          => iso_cnveg_carbonflux_inst%m_livestemc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_to_litter          => iso_cnveg_carbonflux_inst%m_deadstemc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_to_litter         => iso_cnveg_carbonflux_inst%m_livecrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_to_litter         => iso_cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          m_leafc_storage_to_litter      => iso_cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_storage_to_litter     => iso_cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_storage_to_litter  => iso_cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_storage_to_litter  => iso_cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_storage_to_litter => iso_cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_storage_to_litter => iso_cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          m_gresp_storage_to_litter      => iso_cnveg_carbonflux_inst%m_gresp_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          m_leafc_xfer_to_litter         => iso_cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_xfer_to_litter        => iso_cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          m_gresp_xfer_to_litter         => iso_cnveg_carbonflux_inst%m_gresp_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          
          gap_mortality_c_to_litr_met_c  => iso_cnveg_carbonflux_inst%gap_mortality_c_to_litr_met_c_col    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
          gap_mortality_c_to_litr_cel_c  => iso_cnveg_carbonflux_inst%gap_mortality_c_to_litr_cel_c_col    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
          gap_mortality_c_to_litr_lig_c  => iso_cnveg_carbonflux_inst%gap_mortality_c_to_litr_lig_c_col    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
          gap_mortality_c_to_cwdc        => iso_cnveg_carbonflux_inst%gap_mortality_c_to_cwdc_col            & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
          )
          
       do j = 1, nlevdecomp
          do pi = 1,maxpatch_pft
             do fc = 1,num_soilc
                c = filter_soilc(fc)

                if (pi <=  col%npatches(c)) then
                   p = col%patchi(c) + pi - 1

                   if (patch%active(p)) then

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
   subroutine CNCIsoHarvestPftToColumn (num_soilc, filter_soilc, &
        soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)
     !
     ! !DESCRIPTION:
     ! gather all patch-level harvest mortality fluxes
     ! to the column level and assign them to the litter, cwd, and wood product pools
     !
     ! !ARGUMENTS:
     integer                         , intent(in)    :: num_soilc         ! number of soil columns in filter
     integer                         , intent(in)    :: filter_soilc(:)   ! soil column filter
     type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
     type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fc,c,pi,p,j               ! indices
     !-----------------------------------------------------------------------

     associate(                                                                                                  & 
          ivt                              => patch%itype                                                      , & ! Input:  [integer  (:)   ]  patch vegetation type                                
          wtcol                            => patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
          
          lf_flab                          => pftcon%lf_flab                                                   , & ! Input:  leaf litter labile fraction                       
          lf_fcel                          => pftcon%lf_fcel                                                   , & ! Input:  leaf litter cellulose fraction                    
          lf_flig                          => pftcon%lf_flig                                                   , & ! Input:  leaf litter lignin fraction                       
          fr_flab                          => pftcon%fr_flab                                                   , & ! Input:  fine root litter labile fraction                  
          fr_fcel                          => pftcon%fr_fcel                                                   , & ! Input:  fine root litter cellulose fraction               
          fr_flig                          => pftcon%fr_flig                                                   , & ! Input:  fine root litter lignin fraction                  
          
          leaf_prof                        => soilbiogeochem_state_inst%leaf_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
          froot_prof                       => soilbiogeochem_state_inst%froot_prof_patch                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
          croot_prof                       => soilbiogeochem_state_inst%croot_prof_patch                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
          stem_prof                        => soilbiogeochem_state_inst%stem_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
          
          hrv_leafc_to_litter              => iso_cnveg_carbonflux_inst%hrv_leafc_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_to_litter             => iso_cnveg_carbonflux_inst%hrv_frootc_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_to_litter          => iso_cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
          phrv_deadstemc_to_prod10c        => iso_cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          phrv_deadstemc_to_prod100c       => iso_cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_to_litter         => iso_cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_to_litter         => iso_cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_leafc_storage_to_litter      => iso_cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_storage_to_litter     => iso_cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_storage_to_litter  => iso_cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadstemc_storage_to_litter  => iso_cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_storage_to_litter => iso_cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_storage_to_litter => iso_cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_gresp_storage_to_litter      => iso_cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_leafc_xfer_to_litter         => iso_cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_xfer_to_litter        => iso_cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadstemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_gresp_xfer_to_litter         => iso_cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          chrv_deadstemc_to_prod10c        => iso_cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_col           , & ! Output: [real(r8) (:)   ]                                                    
          chrv_deadstemc_to_prod100c       => iso_cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_col          , & ! Output: [real(r8) (:)   ]                                                    
          harvest_c_to_litr_met_c          => iso_cnveg_carbonflux_inst%harvest_c_to_litr_met_c_col            , & ! Output: [real(r8) (:,:) ]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
          harvest_c_to_litr_cel_c          => iso_cnveg_carbonflux_inst%harvest_c_to_litr_cel_c_col            , & ! Output: [real(r8) (:,:) ]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
          harvest_c_to_litr_lig_c          => iso_cnveg_carbonflux_inst%harvest_c_to_litr_lig_c_col            , & ! Output: [real(r8) (:,:) ]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
          harvest_c_to_cwdc                => iso_cnveg_carbonflux_inst%harvest_c_to_cwdc_col                    & ! Output: [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
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
   subroutine CIsoFluxCalc(&
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
        call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
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

   end subroutine CIsoFluxCalc

end module CNCIsoFluxMod
