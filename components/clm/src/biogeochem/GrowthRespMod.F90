module GrowthRespMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for growth respiration fluxes,
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use pftvarcon        , only : grperc, grpnow, npcropmin
  use VegetationPropertiesType   , only : veg_vp
  use CNCarbonFluxType , only : carbonflux_type
  use VegetationType        , only : veg_pp                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: GrowthResp
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine GrowthResp(num_soilp, filter_soilp, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_soilp       ! number of soil patches in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil patches
    type(carbonflux_type), intent(inout) :: carbonflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p                ! indices
    integer :: fp               ! lake filter pft index
    !-----------------------------------------------------------------------

    associate(                                                                                     & 
         ivt                           =>    veg_pp%itype                                           , & ! Input:  [integer (:)]  pft vegetation type                                

         woody                         =>    veg_vp%woody                                    , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)

         cpool_to_leafc                =>    carbonflux_vars%cpool_to_leafc_patch                , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_leafc_storage        =>    carbonflux_vars%cpool_to_leafc_storage_patch        , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_frootc               =>    carbonflux_vars%cpool_to_frootc_patch               , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_frootc_storage       =>    carbonflux_vars%cpool_to_frootc_storage_patch       , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livestemc            =>    carbonflux_vars%cpool_to_livestemc_patch            , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livestemc_storage    =>    carbonflux_vars%cpool_to_livestemc_storage_patch    , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_deadstemc            =>    carbonflux_vars%cpool_to_deadstemc_patch            , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_deadstemc_storage    =>    carbonflux_vars%cpool_to_deadstemc_storage_patch    , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livecrootc           =>    carbonflux_vars%cpool_to_livecrootc_patch           , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livecrootc_storage   =>    carbonflux_vars%cpool_to_livecrootc_storage_patch   , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_deadcrootc           =>    carbonflux_vars%cpool_to_deadcrootc_patch           , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C (gC/m2/s)        
         cpool_to_deadcrootc_storage   =>    carbonflux_vars%cpool_to_deadcrootc_storage_patch   , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C storage (gC/m2/s)
         cpool_to_grainc               =>    carbonflux_vars%cpool_to_grainc_patch               , & ! Input:  [real(r8) (:)]  allocation to grain C (gC/m2/s)                   
         cpool_to_grainc_storage       =>    carbonflux_vars%cpool_to_grainc_storage_patch       , & ! Input:  [real(r8) (:)]  allocation to grain C storage (gC/m2/s)           
         grainc_xfer_to_grainc         =>    carbonflux_vars%grainc_xfer_to_grainc_patch         , & ! Input:  [real(r8) (:)]  grain C growth from storage (gC/m2/s)             
         leafc_xfer_to_leafc           =>    carbonflux_vars%leafc_xfer_to_leafc_patch           , & ! Input:  [real(r8) (:)]  leaf C growth from storage (gC/m2/s)              
         frootc_xfer_to_frootc         =>    carbonflux_vars%frootc_xfer_to_frootc_patch         , & ! Input:  [real(r8) (:)]  fine root C growth from storage (gC/m2/s)         
         livestemc_xfer_to_livestemc   =>    carbonflux_vars%livestemc_xfer_to_livestemc_patch   , & ! Input:  [real(r8) (:)]  live stem C growth from storage (gC/m2/s)         
         deadstemc_xfer_to_deadstemc   =>    carbonflux_vars%deadstemc_xfer_to_deadstemc_patch   , & ! Input:  [real(r8) (:)]  dead stem C growth from storage (gC/m2/s)         
         livecrootc_xfer_to_livecrootc =>    carbonflux_vars%livecrootc_xfer_to_livecrootc_patch , & ! Input:  [real(r8) (:)]  live coarse root C growth from storage (gC/m2/s)  
         deadcrootc_xfer_to_deadcrootc =>    carbonflux_vars%deadcrootc_xfer_to_deadcrootc_patch , & ! Input:  [real(r8) (:)]  dead coarse root C growth from storage (gC/m2/s)  
         cpool_grain_gr                =>    carbonflux_vars%cpool_grain_gr_patch                , & ! InOut:  [real(r8) (:)]                                                    
         cpool_grain_storage_gr        =>    carbonflux_vars%cpool_grain_storage_gr_patch        , & ! InOut:  [real(r8) (:)]                                                    
         transfer_grain_gr             =>    carbonflux_vars%transfer_grain_gr_patch             , & ! InOut:  [real(r8) (:)]                                                    
         cpool_leaf_gr                 =>    carbonflux_vars%cpool_leaf_gr_patch                 , & ! InOut:  [real(r8) (:)]                                                    
         cpool_leaf_storage_gr         =>    carbonflux_vars%cpool_leaf_storage_gr_patch         , & ! InOut:  [real(r8) (:)]                                                    
         transfer_leaf_gr              =>    carbonflux_vars%transfer_leaf_gr_patch              , & ! InOut:  [real(r8) (:)]                                                    
         cpool_froot_gr                =>    carbonflux_vars%cpool_froot_gr_patch                , & ! InOut:  [real(r8) (:)]                                                    
         cpool_froot_storage_gr        =>    carbonflux_vars%cpool_froot_storage_gr_patch        , & ! InOut:  [real(r8) (:)]                                                    
         transfer_froot_gr             =>    carbonflux_vars%transfer_froot_gr_patch             , & ! InOut:  [real(r8) (:)]                                                    
         cpool_livestem_gr             =>    carbonflux_vars%cpool_livestem_gr_patch             , & ! InOut:  [real(r8) (:)]                                                    
         cpool_livestem_storage_gr     =>    carbonflux_vars%cpool_livestem_storage_gr_patch     , & ! InOut:  [real(r8) (:)]                                                    
         transfer_livestem_gr          =>    carbonflux_vars%transfer_livestem_gr_patch          , & ! InOut:  [real(r8) (:)]                                                    
         cpool_deadstem_gr             =>    carbonflux_vars%cpool_deadstem_gr_patch             , & ! InOut:  [real(r8) (:)]                                                    
         cpool_deadstem_storage_gr     =>    carbonflux_vars%cpool_deadstem_storage_gr_patch     , & ! InOut:  [real(r8) (:)]                                                    
         transfer_deadstem_gr          =>    carbonflux_vars%transfer_deadstem_gr_patch          , & ! InOut:  [real(r8) (:)]                                                    
         cpool_livecroot_gr            =>    carbonflux_vars%cpool_livecroot_gr_patch            , & ! InOut:  [real(r8) (:)]                                                    
         cpool_livecroot_storage_gr    =>    carbonflux_vars%cpool_livecroot_storage_gr_patch    , & ! InOut:  [real(r8) (:)]                                                    
         transfer_livecroot_gr         =>    carbonflux_vars%transfer_livecroot_gr_patch         , & ! InOut:  [real(r8) (:)]                                                    
         cpool_deadcroot_gr            =>    carbonflux_vars%cpool_deadcroot_gr_patch            , & ! InOut:  [real(r8) (:)]                                                    
         cpool_deadcroot_storage_gr    =>    carbonflux_vars%cpool_deadcroot_storage_gr_patch    , & ! InOut:  [real(r8) (:)]                                                    
         transfer_deadcroot_gr         =>    carbonflux_vars%transfer_deadcroot_gr_patch           & ! InOut:  [real(r8) (:)]                                                    
         )
      
      ! Loop through patches
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

  end subroutine GrowthResp

end module GrowthRespMod
