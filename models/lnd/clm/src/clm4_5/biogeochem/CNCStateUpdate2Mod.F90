module CNCStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use abortutils   , only: endrun
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate2
  public:: CStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by gap-phase mortality fluxes
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_step_size
    use clm_varpar   , only: nlevsoi, nlevdecomp
    use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
    !
    ! !LOCAL VARIABLES:
    type(pft_cflux_type), pointer :: pcisof
    type(pft_cstate_type), pointer :: pcisos
    type(column_cflux_type), pointer :: ccisof
    type(column_cstate_type), pointer :: ccisos
    integer :: c,p,j      ! indices
    integer :: fp,fc    ! lake filter indices
    real(r8):: dt       ! radiation time step (seconds)
    !-----------------------------------------------------------------------

   ! select which isotope
   select case (isotope)
   case ('bulk')
      pcisof =>  pcf
      pcisos =>  pcs
      ccisof =>  ccf
      ccisos =>  ccs
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
      ccisof =>  cc14f
      ccisos =>  cc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
      ccisof =>  cc13f
      ccisos =>  cc13s
   case default
      call endrun(msg='CNCIsoStateUpdate2Mod: iso must be bulk, c13 or c14'//&
           errMsg(__FILE__, __LINE__))
   end select

   associate(& 
   gap_mortality_c_to_litr_met_c       =>    ccisof%gap_mortality_c_to_litr_met_c        , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
   gap_mortality_c_to_litr_cel_c       =>    ccisof%gap_mortality_c_to_litr_cel_c        , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
   gap_mortality_c_to_litr_lig_c       =>    ccisof%gap_mortality_c_to_litr_lig_c        , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
   gap_mortality_c_to_cwdc             =>    ccisof%gap_mortality_c_to_cwdc              , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
   decomp_cpools_vr                    =>    ccisos%decomp_cpools_vr                     , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   m_deadcrootc_storage_to_litter      =>    pcisof%m_deadcrootc_storage_to_litter       , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_litter              =>    pcisof%m_deadcrootc_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter         =>    pcisof%m_deadcrootc_xfer_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter       =>    pcisof%m_deadstemc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_litter               =>    pcisof%m_deadstemc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter          =>    pcisof%m_deadstemc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_litter          =>    pcisof%m_frootc_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_litter                  =>    pcisof%m_frootc_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter             =>    pcisof%m_frootc_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_litter           =>    pcisof%m_gresp_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter              =>    pcisof%m_gresp_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_litter           =>    pcisof%m_leafc_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_to_litter                   =>    pcisof%m_leafc_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter              =>    pcisof%m_leafc_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter      =>    pcisof%m_livecrootc_storage_to_litter       , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_litter              =>    pcisof%m_livecrootc_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter         =>    pcisof%m_livecrootc_xfer_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter       =>    pcisof%m_livestemc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_litter               =>    pcisof%m_livestemc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter          =>    pcisof%m_livestemc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   deadcrootc                          =>    pcisos%deadcrootc                           , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage                  =>    pcisos%deadcrootc_storage                   , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer                     =>    pcisos%deadcrootc_xfer                      , & ! InOut:  [real(r8) (:)] (gC/m2) dead coarse root C transfer                
   deadstemc                           =>    pcisos%deadstemc                            , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc_storage                   =>    pcisos%deadstemc_storage                    , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc_xfer                      =>    pcisos%deadstemc_xfer                       , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   frootc                              =>    pcisos%frootc                               , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage                      =>    pcisos%frootc_storage                       , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                         =>    pcisos%frootc_xfer                          , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   gresp_storage                       =>    pcisos%gresp_storage                        , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   gresp_xfer                          =>    pcisos%gresp_xfer                           , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafc                               =>    pcisos%leafc                                , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage                       =>    pcisos%leafc_storage                        , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                          =>    pcisos%leafc_xfer                           , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   livecrootc                          =>    pcisos%livecrootc                           , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage                  =>    pcisos%livecrootc_storage                   , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer                     =>    pcisos%livecrootc_xfer                      , & ! InOut:  [real(r8) (:)] (gC/m2) live coarse root C transfer                
   livestemc                           =>    pcisos%livestemc                            , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc_storage                   =>    pcisos%livestemc_storage                    , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc_xfer                      =>    pcisos%livestemc_xfer                         & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   )

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! column level carbon fluxes from gap-phase mortality
    do j = 1,nlevdecomp
       ! column loop
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          
          ! column gap mortality fluxes
          decomp_cpools_vr(c,j,i_met_lit) = decomp_cpools_vr(c,j,i_met_lit) + gap_mortality_c_to_litr_met_c(c,j) * dt
          decomp_cpools_vr(c,j,i_cel_lit) = decomp_cpools_vr(c,j,i_cel_lit) + gap_mortality_c_to_litr_cel_c(c,j) * dt
          decomp_cpools_vr(c,j,i_lig_lit) = decomp_cpools_vr(c,j,i_lig_lit) + gap_mortality_c_to_litr_lig_c(c,j) * dt
          decomp_cpools_vr(c,j,i_cwd) = decomp_cpools_vr(c,j,i_cwd) + gap_mortality_c_to_cwdc(c,j) * dt
          
       end do
    end do 

    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! pft-level carbon fluxes from gap-phase mortality
       ! displayed pools
       leafc(p)               = leafc(p)              - m_leafc_to_litter(p)              * dt
       frootc(p)              = frootc(p)             - m_frootc_to_litter(p)             * dt
       livestemc(p)           = livestemc(p)          - m_livestemc_to_litter(p)          * dt
       deadstemc(p)           = deadstemc(p)          - m_deadstemc_to_litter(p)          * dt
       livecrootc(p)          = livecrootc(p)         - m_livecrootc_to_litter(p)         * dt
       deadcrootc(p)          = deadcrootc(p)         - m_deadcrootc_to_litter(p)         * dt

       ! storage pools
       leafc_storage(p)       = leafc_storage(p)      - m_leafc_storage_to_litter(p)      * dt
       frootc_storage(p)      = frootc_storage(p)     - m_frootc_storage_to_litter(p)     * dt
       livestemc_storage(p)   = livestemc_storage(p)  - m_livestemc_storage_to_litter(p)  * dt
       deadstemc_storage(p)   = deadstemc_storage(p)  - m_deadstemc_storage_to_litter(p)  * dt
       livecrootc_storage(p)  = livecrootc_storage(p) - m_livecrootc_storage_to_litter(p) * dt
       deadcrootc_storage(p)  = deadcrootc_storage(p) - m_deadcrootc_storage_to_litter(p) * dt
       gresp_storage(p)       = gresp_storage(p)      - m_gresp_storage_to_litter(p)      * dt

       ! transfer pools
       leafc_xfer(p)          = leafc_xfer(p)         - m_leafc_xfer_to_litter(p)         * dt
       frootc_xfer(p)         = frootc_xfer(p)        - m_frootc_xfer_to_litter(p)        * dt
       livestemc_xfer(p)      = livestemc_xfer(p)     - m_livestemc_xfer_to_litter(p)     * dt
       deadstemc_xfer(p)      = deadstemc_xfer(p)     - m_deadstemc_xfer_to_litter(p)     * dt
       livecrootc_xfer(p)     = livecrootc_xfer(p)    - m_livecrootc_xfer_to_litter(p)    * dt
       deadcrootc_xfer(p)     = deadcrootc_xfer(p)    - m_deadcrootc_xfer_to_litter(p)    * dt
       gresp_xfer(p)          = gresp_xfer(p)         - m_gresp_xfer_to_litter(p)         * dt
    end do ! end of pft loop

    end associate 
 end subroutine CStateUpdate2

 !-----------------------------------------------------------------------
 subroutine CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
   !
   ! !DESCRIPTION:
   ! Update all the prognostic carbon state
   ! variables affected by harvest mortality fluxes
   !
   ! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varpar   , only: nlevdecomp
   use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
   !
   ! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
   !
   ! !LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   type(column_cflux_type), pointer :: ccisof
   type(column_cstate_type), pointer :: ccisos
   integer :: c,p,j,k,l       ! indices
   integer :: fp,fc    ! lake filter indices
   real(r8):: dt       ! radiation time step (seconds)
   !-----------------------------------------------------------------------

   ! select which isotope
   select case (isotope)
   case ('bulk')
      pcisof =>  pcf
      pcisos =>  pcs
      ccisof =>  ccf
      ccisos =>  ccs
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
      ccisof =>  cc14f
      ccisos =>  cc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
      ccisof =>  cc13f
      ccisos =>  cc13s
   case default
      call endrun(msg='CNCIsoStateUpdate2Mod: iso must be bulk, c13 or c14'//&
           errMsg(__FILE__, __LINE__))
   end select

   associate(& 
   harvest_c_to_litr_met_c             =>    ccisof%harvest_c_to_litr_met_c              , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
   harvest_c_to_litr_cel_c             =>    ccisof%harvest_c_to_litr_cel_c              , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
   harvest_c_to_litr_lig_c             =>    ccisof%harvest_c_to_litr_lig_c              , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
   harvest_c_to_cwdc                   =>    ccisof%harvest_c_to_cwdc                    , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to CWD pool (gC/m3/s)
   decomp_cpools_vr                    =>    ccisos%decomp_cpools_vr                     , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

   hrv_deadcrootc_storage_to_litter    =>    pcisof%hrv_deadcrootc_storage_to_litter     , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_to_litter            =>    pcisof%hrv_deadcrootc_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_xfer_to_litter       =>    pcisof%hrv_deadcrootc_xfer_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_storage_to_litter     =>    pcisof%hrv_deadstemc_storage_to_litter      , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod10c            =>    pcisof%hrv_deadstemc_to_prod10c             , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod100c           =>    pcisof%hrv_deadstemc_to_prod100c            , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_xfer_to_litter        =>    pcisof%hrv_deadstemc_xfer_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_storage_to_litter        =>    pcisof%hrv_frootc_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_to_litter                =>    pcisof%hrv_frootc_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_xfer_to_litter           =>    pcisof%hrv_frootc_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_storage_to_litter         =>    pcisof%hrv_gresp_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_xfer_to_litter            =>    pcisof%hrv_gresp_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_storage_to_litter         =>    pcisof%hrv_leafc_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_to_litter                 =>    pcisof%hrv_leafc_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_xfer_to_litter            =>    pcisof%hrv_leafc_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_storage_to_litter    =>    pcisof%hrv_livecrootc_storage_to_litter     , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_to_litter            =>    pcisof%hrv_livecrootc_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_xfer_to_litter       =>    pcisof%hrv_livecrootc_xfer_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_storage_to_litter     =>    pcisof%hrv_livestemc_storage_to_litter      , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_to_litter             =>    pcisof%hrv_livestemc_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_xfer_to_litter        =>    pcisof%hrv_livestemc_xfer_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_xsmrpool_to_atm                 =>    pcisof%hrv_xsmrpool_to_atm                  , & ! Input:  [real(r8) (:)]                                                    
   deadcrootc                          =>    pcisos%deadcrootc                           , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage                  =>    pcisos%deadcrootc_storage                   , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer                     =>    pcisos%deadcrootc_xfer                      , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   deadstemc                           =>    pcisos%deadstemc                            , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc_storage                   =>    pcisos%deadstemc_storage                    , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc_xfer                      =>    pcisos%deadstemc_xfer                       , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   frootc                              =>    pcisos%frootc                               , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage                      =>    pcisos%frootc_storage                       , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                         =>    pcisos%frootc_xfer                          , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   gresp_storage                       =>    pcisos%gresp_storage                        , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   gresp_xfer                          =>    pcisos%gresp_xfer                           , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafc                               =>    pcisos%leafc                                , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage                       =>    pcisos%leafc_storage                        , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                          =>    pcisos%leafc_xfer                           , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   livecrootc                          =>    pcisos%livecrootc                           , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage                  =>    pcisos%livecrootc_storage                   , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer                     =>    pcisos%livecrootc_xfer                      , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   livestemc                           =>    pcisos%livestemc                            , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc_storage                   =>    pcisos%livestemc_storage                    , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc_xfer                      =>    pcisos%livestemc_xfer                       , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   xsmrpool                            =>    pcisos%xsmrpool                               & ! InOut:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
   )

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! column level carbon fluxes from harvest mortality
    do j = 1, nlevdecomp
       ! column loop
       do fc = 1,num_soilc
          c = filter_soilc(fc)

          ! column harvest fluxes
          decomp_cpools_vr(c,j,i_met_lit) = decomp_cpools_vr(c,j,i_met_lit) + harvest_c_to_litr_met_c(c,j) * dt
          decomp_cpools_vr(c,j,i_cel_lit) = decomp_cpools_vr(c,j,i_cel_lit) + harvest_c_to_litr_cel_c(c,j) * dt
          decomp_cpools_vr(c,j,i_lig_lit) = decomp_cpools_vr(c,j,i_lig_lit) + harvest_c_to_litr_lig_c(c,j) * dt
          decomp_cpools_vr(c,j,i_cwd) = decomp_cpools_vr(c,j,i_cwd) + harvest_c_to_cwdc(c,j)  * dt

          ! wood to product pools - states updated in CNWoodProducts()
       end do
    end do
    
    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! pft-level carbon fluxes from harvest mortality
       ! displayed pools
       leafc(p)               = leafc(p)              - hrv_leafc_to_litter(p)              * dt
       frootc(p)              = frootc(p)             - hrv_frootc_to_litter(p)             * dt
       livestemc(p)           = livestemc(p)          - hrv_livestemc_to_litter(p)          * dt
       deadstemc(p)           = deadstemc(p)          - hrv_deadstemc_to_prod10c(p)         * dt
       deadstemc(p)           = deadstemc(p)          - hrv_deadstemc_to_prod100c(p)        * dt
       livecrootc(p)          = livecrootc(p)         - hrv_livecrootc_to_litter(p)         * dt
       deadcrootc(p)          = deadcrootc(p)         - hrv_deadcrootc_to_litter(p)         * dt
       
       ! xsmrpool
       xsmrpool(p)            = xsmrpool(p)           - hrv_xsmrpool_to_atm(p)              * dt

       ! storage pools
       leafc_storage(p)       = leafc_storage(p)      - hrv_leafc_storage_to_litter(p)      * dt
       frootc_storage(p)      = frootc_storage(p)     - hrv_frootc_storage_to_litter(p)     * dt
       livestemc_storage(p)   = livestemc_storage(p)  - hrv_livestemc_storage_to_litter(p)  * dt
       deadstemc_storage(p)   = deadstemc_storage(p)  - hrv_deadstemc_storage_to_litter(p)  * dt
       livecrootc_storage(p)  = livecrootc_storage(p) - hrv_livecrootc_storage_to_litter(p) * dt
       deadcrootc_storage(p)  = deadcrootc_storage(p) - hrv_deadcrootc_storage_to_litter(p) * dt
       gresp_storage(p)       = gresp_storage(p)      - hrv_gresp_storage_to_litter(p)      * dt

       ! transfer pools
       leafc_xfer(p)          = leafc_xfer(p)         - hrv_leafc_xfer_to_litter(p)         * dt
       frootc_xfer(p)         = frootc_xfer(p)        - hrv_frootc_xfer_to_litter(p)        * dt
       livestemc_xfer(p)      = livestemc_xfer(p)     - hrv_livestemc_xfer_to_litter(p)     * dt
       deadstemc_xfer(p)      = deadstemc_xfer(p)     - hrv_deadstemc_xfer_to_litter(p)     * dt
       livecrootc_xfer(p)     = livecrootc_xfer(p)    - hrv_livecrootc_xfer_to_litter(p)    * dt
       deadcrootc_xfer(p)     = deadcrootc_xfer(p)    - hrv_deadcrootc_xfer_to_litter(p)    * dt
       gresp_xfer(p)          = gresp_xfer(p)         - hrv_gresp_xfer_to_litter(p)         * dt

    end do ! end of pft loop

    end associate 
 end subroutine CStateUpdate2h

end module CNCStateUpdate2Mod
