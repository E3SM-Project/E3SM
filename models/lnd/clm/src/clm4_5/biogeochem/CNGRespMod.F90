module CNGRespMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for growth respiration fluxes,
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNGResp
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNGResp(num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables
    !
    ! !USES:
    use clmtype
    use pftvarcon, only : npcropmin, grperc, grpnow
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    !
    ! !LOCAL VARIABLES:
    integer :: p                ! indices
    integer :: fp               ! lake filter pft index
    !-----------------------------------------------------------------------

   associate(& 
   ivt                                 =>   pft%itype                                    , & ! Input:  [integer (:)]  pft vegetation type                                
   cpool_to_leafc                      =>    pcf%cpool_to_leafc                          , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_leafc_storage              =>    pcf%cpool_to_leafc_storage                  , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_frootc                     =>    pcf%cpool_to_frootc                         , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_frootc_storage             =>    pcf%cpool_to_frootc_storage                 , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livestemc                  =>    pcf%cpool_to_livestemc                      , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livestemc_storage          =>    pcf%cpool_to_livestemc_storage              , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_deadstemc                  =>    pcf%cpool_to_deadstemc                      , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_deadstemc_storage          =>    pcf%cpool_to_deadstemc_storage              , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livecrootc                 =>    pcf%cpool_to_livecrootc                     , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livecrootc_storage         =>    pcf%cpool_to_livecrootc_storage             , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_deadcrootc                 =>    pcf%cpool_to_deadcrootc                     , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C (gC/m2/s)        
   cpool_to_deadcrootc_storage         =>    pcf%cpool_to_deadcrootc_storage             , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C storage (gC/m2/s)
   cpool_to_grainc                     =>    pcf%cpool_to_grainc                         , & ! Input:  [real(r8) (:)]  allocation to grain C (gC/m2/s)                   
   cpool_to_grainc_storage             =>    pcf%cpool_to_grainc_storage                 , & ! Input:  [real(r8) (:)]  allocation to grain C storage (gC/m2/s)           
   grainc_xfer_to_grainc               =>    pcf%grainc_xfer_to_grainc                   , & ! Input:  [real(r8) (:)]  grain C growth from storage (gC/m2/s)             
   leafc_xfer_to_leafc                 =>    pcf%leafc_xfer_to_leafc                     , & ! Input:  [real(r8) (:)]  leaf C growth from storage (gC/m2/s)              
   frootc_xfer_to_frootc               =>    pcf%frootc_xfer_to_frootc                   , & ! Input:  [real(r8) (:)]  fine root C growth from storage (gC/m2/s)         
   livestemc_xfer_to_livestemc         =>    pcf%livestemc_xfer_to_livestemc             , & ! Input:  [real(r8) (:)]  live stem C growth from storage (gC/m2/s)         
   deadstemc_xfer_to_deadstemc         =>    pcf%deadstemc_xfer_to_deadstemc             , & ! Input:  [real(r8) (:)]  dead stem C growth from storage (gC/m2/s)         
   livecrootc_xfer_to_livecrootc       =>    pcf%livecrootc_xfer_to_livecrootc           , & ! Input:  [real(r8) (:)]  live coarse root C growth from storage (gC/m2/s)  
   deadcrootc_xfer_to_deadcrootc       =>    pcf%deadcrootc_xfer_to_deadcrootc           , & ! Input:  [real(r8) (:)]  dead coarse root C growth from storage (gC/m2/s)  
   woody                               =>    pftcon%woody                                , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)
   cpool_grain_gr                      =>    pcf%cpool_grain_gr                          , & ! InOut:  [real(r8) (:)]                                                    
   cpool_grain_storage_gr              =>    pcf%cpool_grain_storage_gr                  , & ! InOut:  [real(r8) (:)]                                                    
   transfer_grain_gr                   =>    pcf%transfer_grain_gr                       , & ! InOut:  [real(r8) (:)]                                                    
   cpool_leaf_gr                       =>    pcf%cpool_leaf_gr                           , & ! InOut:  [real(r8) (:)]                                                    
   cpool_leaf_storage_gr               =>    pcf%cpool_leaf_storage_gr                   , & ! InOut:  [real(r8) (:)]                                                    
   transfer_leaf_gr                    =>    pcf%transfer_leaf_gr                        , & ! InOut:  [real(r8) (:)]                                                    
   cpool_froot_gr                      =>    pcf%cpool_froot_gr                          , & ! InOut:  [real(r8) (:)]                                                    
   cpool_froot_storage_gr              =>    pcf%cpool_froot_storage_gr                  , & ! InOut:  [real(r8) (:)]                                                    
   transfer_froot_gr                   =>    pcf%transfer_froot_gr                       , & ! InOut:  [real(r8) (:)]                                                    
   cpool_livestem_gr                   =>    pcf%cpool_livestem_gr                       , & ! InOut:  [real(r8) (:)]                                                    
   cpool_livestem_storage_gr           =>    pcf%cpool_livestem_storage_gr               , & ! InOut:  [real(r8) (:)]                                                    
   transfer_livestem_gr                =>    pcf%transfer_livestem_gr                    , & ! InOut:  [real(r8) (:)]                                                    
   cpool_deadstem_gr                   =>    pcf%cpool_deadstem_gr                       , & ! InOut:  [real(r8) (:)]                                                    
   cpool_deadstem_storage_gr           =>    pcf%cpool_deadstem_storage_gr               , & ! InOut:  [real(r8) (:)]                                                    
   transfer_deadstem_gr                =>    pcf%transfer_deadstem_gr                    , & ! InOut:  [real(r8) (:)]                                                    
   cpool_livecroot_gr                  =>    pcf%cpool_livecroot_gr                      , & ! InOut:  [real(r8) (:)]                                                    
   cpool_livecroot_storage_gr          =>    pcf%cpool_livecroot_storage_gr              , & ! InOut:  [real(r8) (:)]                                                    
   transfer_livecroot_gr               =>    pcf%transfer_livecroot_gr                   , & ! InOut:  [real(r8) (:)]                                                    
   cpool_deadcroot_gr                  =>    pcf%cpool_deadcroot_gr                      , & ! InOut:  [real(r8) (:)]                                                    
   cpool_deadcroot_storage_gr          =>    pcf%cpool_deadcroot_storage_gr              , & ! InOut:  [real(r8) (:)]                                                    
   transfer_deadcroot_gr               =>    pcf%transfer_deadcroot_gr                     & ! InOut:  [real(r8) (:)]                                                    
   )

   ! Loop through pfts
   ! start pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         cpool_livestem_gr(p)          = cpool_to_livestemc(p) * grperc(ivt(p))
         cpool_livestem_storage_gr(p)  = cpool_to_livestemc_storage(p) * &
                                         grperc(ivt(p)) * grpnow(ivt(p))
         transfer_livestem_gr(p)       = livestemc_xfer_to_livestemc(p) * &
                                         grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
         cpool_grain_gr(p)             = cpool_to_grainc(p) * grperc(ivt(p))
         cpool_grain_storage_gr(p)     = cpool_to_grainc_storage(p) * &
                                         grperc(ivt(p)) * grpnow(ivt(p))
         transfer_grain_gr(p)          = grainc_xfer_to_grainc(p) * grperc(ivt(p)) &
                                         * (1._r8 - grpnow(ivt(p)))
      end if

      ! leaf and fine root growth respiration
      cpool_leaf_gr(p)          = cpool_to_leafc(p) * grperc(ivt(p))
      cpool_leaf_storage_gr(p)  = cpool_to_leafc_storage(p) * grperc(ivt(p)) * &
                                  grpnow(ivt(p))
      transfer_leaf_gr(p)       = leafc_xfer_to_leafc(p) * grperc(ivt(p)) * &
                                  (1._r8 - grpnow(ivt(p)))
      cpool_froot_gr(p)         = cpool_to_frootc(p) * grperc(ivt(p))
      cpool_froot_storage_gr(p) = cpool_to_frootc_storage(p) * grperc(ivt(p)) * &
                                  grpnow(ivt(p))
      transfer_froot_gr(p)      = frootc_xfer_to_frootc(p) * grperc(ivt(p)) * &
                                  (1._r8 - grpnow(ivt(p)))

      if (woody(ivt(p)) == 1._r8) then
          cpool_livestem_gr(p)          = cpool_to_livestemc(p) * grperc(ivt(p))
          cpool_livestem_storage_gr(p)  = cpool_to_livestemc_storage(p) * &
                                          grperc(ivt(p)) * grpnow(ivt(p))
          transfer_livestem_gr(p)       = livestemc_xfer_to_livestemc(p) * &
                                          grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
          cpool_deadstem_gr(p)          = cpool_to_deadstemc(p) * grperc(ivt(p))
          cpool_deadstem_storage_gr(p)  = cpool_to_deadstemc_storage(p) * &
                                          grperc(ivt(p)) * grpnow(ivt(p))
          transfer_deadstem_gr(p)       = deadstemc_xfer_to_deadstemc(p) * &
                                          grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
          cpool_livecroot_gr(p)         = cpool_to_livecrootc(p) * grperc(ivt(p))
          cpool_livecroot_storage_gr(p) = cpool_to_livecrootc_storage(p) * &
                                          grperc(ivt(p)) * grpnow(ivt(p))
          transfer_livecroot_gr(p)      = livecrootc_xfer_to_livecrootc(p) * &
                                          grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
          cpool_deadcroot_gr(p)         = cpool_to_deadcrootc(p) * grperc(ivt(p))
          cpool_deadcroot_storage_gr(p) = cpool_to_deadcrootc_storage(p) * &
                                          grperc(ivt(p)) * grpnow(ivt(p))
          transfer_deadcroot_gr(p)      = deadcrootc_xfer_to_deadcrootc(p) * &
                                          grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
      end if

   end do

    end associate 
 end subroutine CNGResp

end module CNGRespMod
