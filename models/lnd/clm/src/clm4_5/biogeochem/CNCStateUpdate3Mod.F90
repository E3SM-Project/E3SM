module CNCStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by fire fluxes
    !
    ! !USES:
    use clmtype
    use clm_time_manager , only: get_step_size
    use clm_varpar       , only: nlevdecomp, ndecomp_pools
    use clm_varpar       , only: i_cwd, i_met_lit, i_cel_lit, i_lig_lit
    use abortutils       , only: endrun
    use shr_log_mod      , only: errMsg => shr_log_errMsg
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
    integer :: c,p,j,l,k      ! indices
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
      call endrun(msg='CNCIsoStateUpdate3Mod: iso must be bulk, c13 or c14'//&
           errMsg(__FILE__, __LINE__))
   end select

   associate(& 
   fire_mortality_c_to_cwdc            =>    ccisof%fire_mortality_c_to_cwdc             , & ! Input:  [real(r8) (:,:)]  C fluxes associated with fire mortality to CWD pool (gC/m3/s)
   m_decomp_cpools_to_fire_vr          =>    ccisof%m_decomp_cpools_to_fire_vr           , & ! Input:  [real(r8) (:,:,:)]  vertically-resolved decomposing C fire loss (gC/m3/s)
   decomp_cpools_vr                    =>    ccisos%decomp_cpools_vr                     , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   m_c_to_litr_met_fire                =>    ccisof%m_c_to_litr_met_fire                 , & ! Input:  [real(r8) (:,:)]                                                  
   m_c_to_litr_cel_fire                =>    ccisof%m_c_to_litr_cel_fire                 , & ! Input:  [real(r8) (:,:)]                                                  
   m_c_to_litr_lig_fire                =>    ccisof%m_c_to_litr_lig_fire                 , & ! Input:  [real(r8) (:,:)]                                                  
   m_leafc_to_fire                     =>    pcisof%m_leafc_to_fire                      , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_fire             =>    pcisof%m_leafc_storage_to_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_fire                =>    pcisof%m_leafc_xfer_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_fire                 =>    pcisof%m_livestemc_to_fire                  , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_fire         =>    pcisof%m_livestemc_storage_to_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_fire            =>    pcisof%m_livestemc_xfer_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_fire                 =>    pcisof%m_deadstemc_to_fire                  , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_fire         =>    pcisof%m_deadstemc_storage_to_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_fire            =>    pcisof%m_deadstemc_xfer_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_fire                    =>    pcisof%m_frootc_to_fire                     , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_fire            =>    pcisof%m_frootc_storage_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_fire               =>    pcisof%m_frootc_xfer_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_fire                =>    pcisof%m_livecrootc_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_fire        =>    pcisof%m_livecrootc_storage_to_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_fire           =>    pcisof%m_livecrootc_xfer_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_fire                =>    pcisof%m_deadcrootc_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_fire        =>    pcisof%m_deadcrootc_storage_to_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_fire           =>    pcisof%m_deadcrootc_xfer_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_fire             =>    pcisof%m_gresp_storage_to_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_fire                =>    pcisof%m_gresp_xfer_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_to_litter_fire              =>    pcisof%m_leafc_to_litter_fire               , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_litter_fire      =>    pcisof%m_leafc_storage_to_litter_fire       , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter_fire         =>    pcisof%m_leafc_xfer_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_litter_fire          =>    pcisof%m_livestemc_to_litter_fire           , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter_fire  =>    pcisof%m_livestemc_storage_to_litter_fire   , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter_fire     =>    pcisof%m_livestemc_xfer_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_deadstemc_fire       =>    pcisof%m_livestemc_to_deadstemc_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_litter_fire          =>    pcisof%m_deadstemc_to_litter_fire           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter_fire  =>    pcisof%m_deadstemc_storage_to_litter_fire   , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter_fire     =>    pcisof%m_deadstemc_xfer_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_litter_fire             =>    pcisof%m_frootc_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_litter_fire     =>    pcisof%m_frootc_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter_fire        =>    pcisof%m_frootc_xfer_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_litter_fire         =>    pcisof%m_livecrootc_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter_fire  =>    pcisof%m_livecrootc_storage_to_litter_fire  , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter_fire    =>    pcisof%m_livecrootc_xfer_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_deadcrootc_fire     =>    pcisof%m_livecrootc_to_deadcrootc_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_litter_fire         =>    pcisof%m_deadcrootc_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_litter_fire  =>    pcisof%m_deadcrootc_storage_to_litter_fire  , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter_fire    =>    pcisof%m_deadcrootc_xfer_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_litter_fire      =>    pcisof%m_gresp_storage_to_litter_fire       , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter_fire         =>    pcisof%m_gresp_xfer_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
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
   livestemc_xfer                      =>    pcisos%livestemc_xfer                         & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   )

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! column level carbon fluxes from fire
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          ! pft-level wood to column-level CWD (uncombusted wood)
          decomp_cpools_vr(c,j,i_cwd) = decomp_cpools_vr(c,j,i_cwd) + fire_mortality_c_to_cwdc(c,j) * dt
          
          ! pft-level wood to column-level litter (uncombusted wood)
          decomp_cpools_vr(c,j,i_met_lit) = decomp_cpools_vr(c,j,i_met_lit) + m_c_to_litr_met_fire(c,j)* dt
          decomp_cpools_vr(c,j,i_cel_lit) = decomp_cpools_vr(c,j,i_cel_lit) + m_c_to_litr_cel_fire(c,j)* dt
          decomp_cpools_vr(c,j,i_lig_lit) = decomp_cpools_vr(c,j,i_lig_lit) + m_c_to_litr_lig_fire(c,j)* dt
       end do
    end do
    
    ! litter and CWD losses to fire
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             decomp_cpools_vr(c,j,l) = decomp_cpools_vr(c,j,l) - m_decomp_cpools_to_fire_vr(c,j,l) * dt
          end do
       end do
    end do
    
    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! pft-level carbon fluxes from fire
       ! displayed pools
       leafc(p)               = leafc(p)               - m_leafc_to_fire(p)            * dt
       leafc(p)               = leafc(p)               - m_leafc_to_litter_fire(p)     * dt
       frootc(p)              = frootc(p)              - m_frootc_to_fire(p)           * dt
       frootc(p)              = frootc(p)              - m_frootc_to_litter_fire(p)    * dt
       livestemc(p)           = livestemc(p)           - m_livestemc_to_fire(p)        * dt
       livestemc(p)           = livestemc(p)           - m_livestemc_to_litter_fire(p) * dt
       deadstemc(p)           = deadstemc(p)           - m_deadstemc_to_fire(p)        * dt
       deadstemc(p)           = deadstemc(p)           - m_deadstemc_to_litter_fire(p) * dt
       livecrootc(p)          = livecrootc(p)          - m_livecrootc_to_fire(p)       * dt
       livecrootc(p)          = livecrootc(p)          - m_livecrootc_to_litter_fire(p)* dt
       deadcrootc(p)          = deadcrootc(p)          - m_deadcrootc_to_fire(p)       * dt
       deadcrootc(p)          = deadcrootc(p)          - m_deadcrootc_to_litter_fire(p)* dt

       ! storage pools
       leafc_storage(p)       = leafc_storage(p)       - m_leafc_storage_to_fire(p)            * dt
       leafc_storage(p)       = leafc_storage(p)       - m_leafc_storage_to_litter_fire(p)     * dt
       frootc_storage(p)      = frootc_storage(p)      - m_frootc_storage_to_fire(p)           * dt
       frootc_storage(p)      = frootc_storage(p)      - m_frootc_storage_to_litter_fire(p)    * dt
       livestemc_storage(p)   = livestemc_storage(p)   - m_livestemc_storage_to_fire(p)        * dt
       livestemc_storage(p)   = livestemc_storage(p)   - m_livestemc_storage_to_litter_fire(p) * dt
       deadstemc_storage(p)   = deadstemc_storage(p)   - m_deadstemc_storage_to_fire(p)        * dt
       deadstemc_storage(p)   = deadstemc_storage(p)   - m_deadstemc_storage_to_litter_fire(p) * dt
       livecrootc_storage(p)  = livecrootc_storage(p)  - m_livecrootc_storage_to_fire(p)       * dt
       livecrootc_storage(p)  = livecrootc_storage(p)  - m_livecrootc_storage_to_litter_fire(p)* dt
       deadcrootc_storage(p)  = deadcrootc_storage(p)  - m_deadcrootc_storage_to_fire(p)       * dt
       deadcrootc_storage(p)  = deadcrootc_storage(p)  - m_deadcrootc_storage_to_litter_fire(p)* dt
       gresp_storage(p)       = gresp_storage(p)       - m_gresp_storage_to_fire(p)            * dt
       gresp_storage(p)       = gresp_storage(p)       - m_gresp_storage_to_litter_fire(p)     * dt

       ! transfer pools
       leafc_xfer(p)      = leafc_xfer(p)      - m_leafc_xfer_to_fire(p)            * dt
       leafc_xfer(p)      = leafc_xfer(p)      - m_leafc_xfer_to_litter_fire(p)     * dt
       frootc_xfer(p)     = frootc_xfer(p)     - m_frootc_xfer_to_fire(p)           * dt
       frootc_xfer(p)     = frootc_xfer(p)     - m_frootc_xfer_to_litter_fire(p)    * dt
       livestemc_xfer(p)  = livestemc_xfer(p)  - m_livestemc_xfer_to_fire(p)        * dt
       livestemc_xfer(p)  = livestemc_xfer(p)  - m_livestemc_xfer_to_litter_fire(p) * dt
       deadstemc_xfer(p)  = deadstemc_xfer(p)  - m_deadstemc_xfer_to_fire(p)        * dt
       deadstemc_xfer(p)  = deadstemc_xfer(p)  - m_deadstemc_xfer_to_litter_fire(p) * dt
       livecrootc_xfer(p) = livecrootc_xfer(p) - m_livecrootc_xfer_to_fire(p)       * dt
       livecrootc_xfer(p) = livecrootc_xfer(p) - m_livecrootc_xfer_to_litter_fire(p)* dt
       deadcrootc_xfer(p) = deadcrootc_xfer(p) - m_deadcrootc_xfer_to_fire(p)       * dt
       deadcrootc_xfer(p) = deadcrootc_xfer(p) - m_deadcrootc_xfer_to_litter_fire(p)* dt
       gresp_xfer(p)      = gresp_xfer(p)      - m_gresp_xfer_to_fire(p)            * dt
       gresp_xfer(p)      = gresp_xfer(p)      - m_gresp_xfer_to_litter_fire(p)     * dt

    end do ! end of pft loop

    end associate 
 end subroutine CStateUpdate3

end module CNCStateUpdate3Mod
