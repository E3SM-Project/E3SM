
module CNCStateUpdate3Mod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CStateUpdate3Mod
!
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
!
! !REVISION HISTORY:
! 7/27/2004: Created by Peter Thornton
! F. Li and S. Levis (11/06/12)
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CStateUpdate3
!
! !INTERFACE:
subroutine CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic carbon state
! variables affected by fire fluxes
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varpar   , only: nlevdecomp, ndecomp_pools
   use clm_varpar   , only: i_cwd, i_met_lit, i_cel_lit, i_lig_lit
   use abortutils  , only: endrun
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 3/29/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real(r8), pointer :: fire_mortality_c_to_cwdc(:,:)                   ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)
   real(r8), pointer :: m_decomp_cpools_to_fire_vr(:,:,:)               ! vertically-resolved decomposing C fire loss (gC/m3/s)
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:)
   real(r8), pointer :: m_deadcrootc_to_fire(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemc_storage_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_fire(:)
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:)
   real(r8), pointer :: m_frootc_storage_to_fire(:)
   real(r8), pointer :: m_frootc_to_fire(:)
   real(r8), pointer :: m_frootc_xfer_to_fire(:)
   real(r8), pointer :: m_gresp_storage_to_fire(:)
   real(r8), pointer :: m_gresp_xfer_to_fire(:)
   real(r8), pointer :: m_leafc_storage_to_fire(:)
   real(r8), pointer :: m_leafc_to_fire(:)
   real(r8), pointer :: m_leafc_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootc_storage_to_fire(:)
   real(r8), pointer :: m_livecrootc_to_fire(:)
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)
   real(r8), pointer :: m_livestemc_storage_to_fire(:)
   real(r8), pointer :: m_livestemc_to_fire(:)
   real(r8), pointer :: m_livestemc_xfer_to_fire(:)
   !
   real(r8), pointer :: m_leafc_to_litter_fire(:)   
   real(r8), pointer :: m_leafc_storage_to_litter_fire(:)                
   real(r8), pointer :: m_leafc_xfer_to_litter_fire(:)  
   real(r8), pointer :: m_livestemc_to_litter_fire(:)    
   real(r8), pointer :: m_livestemc_storage_to_litter_fire(:)        
   real(r8), pointer :: m_livestemc_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_livestemc_to_deadstemc_fire(:)    
   real(r8), pointer :: m_deadstemc_to_litter_fire(:) 
   real(r8), pointer :: m_deadstemc_storage_to_litter_fire(:)           
   real(r8), pointer :: m_deadstemc_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_frootc_to_litter_fire(:)        
   real(r8), pointer :: m_frootc_storage_to_litter_fire(:)  
   real(r8), pointer :: m_frootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootc_to_litter_fire(:)    
   real(r8), pointer :: m_livecrootc_storage_to_litter_fire(:)      
   real(r8), pointer :: m_livecrootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootc_to_deadcrootc_fire(:)    
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)        
   real(r8), pointer :: m_deadcrootc_storage_to_litter_fire(:)  
   real(r8), pointer :: m_deadcrootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_gresp_storage_to_litter_fire(:)      
   real(r8), pointer :: m_gresp_xfer_to_litter_fire(:) 
   !
   real(r8), pointer :: m_c_to_litr_met_fire(:,:)
   real(r8), pointer :: m_c_to_litr_cel_fire(:,:)
   real(r8), pointer :: m_c_to_litr_lig_fire(:,:)
! local pointers to implicit in/out arrays
   real(r8), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
!
! local pointers to implicit out arrays
!
! !OTHER LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   type(column_cflux_type), pointer :: ccisof
   type(column_cstate_type), pointer :: ccisos
   integer :: c,p,j,l,k      ! indices
   integer :: fp,fc    ! lake filter indices
   real(r8):: dt       ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------
   ! select which isotope
   select case (isotope)
   case ('bulk')
      pcisof => clm3%g%l%c%p%pcf
      pcisos => clm3%g%l%c%p%pcs
      ccisof => clm3%g%l%c%ccf
      ccisos => clm3%g%l%c%ccs
   case ('c14')
      pcisof => clm3%g%l%c%p%pc14f
      pcisos => clm3%g%l%c%p%pc14s
      ccisof => clm3%g%l%c%cc14f
      ccisos => clm3%g%l%c%cc14s
   case ('c13')
      pcisof => clm3%g%l%c%p%pc13f
      pcisos => clm3%g%l%c%p%pc13s
      ccisof => clm3%g%l%c%cc13f
      ccisos => clm3%g%l%c%cc13s
   case default
      call endrun('CNCIsoStateUpdate3Mod: iso must be bulk, c13 or c14')
   end select

    ! assign local pointers at the column level
    fire_mortality_c_to_cwdc       => ccisof%fire_mortality_c_to_cwdc
    m_decomp_cpools_to_fire_vr     => ccisof%m_decomp_cpools_to_fire_vr
    decomp_cpools_vr               => ccisos%decomp_cpools_vr
    m_c_to_litr_met_fire           => ccisof%m_c_to_litr_met_fire
    m_c_to_litr_cel_fire           => ccisof%m_c_to_litr_cel_fire
    m_c_to_litr_lig_fire           => ccisof%m_c_to_litr_lig_fire

    ! assign local pointers at the pft level
    m_leafc_to_fire                => pcisof%m_leafc_to_fire
    m_leafc_storage_to_fire        => pcisof%m_leafc_storage_to_fire
    m_leafc_xfer_to_fire           => pcisof%m_leafc_xfer_to_fire
    m_livestemc_to_fire            => pcisof%m_livestemc_to_fire
    m_livestemc_storage_to_fire    => pcisof%m_livestemc_storage_to_fire
    m_livestemc_xfer_to_fire       => pcisof%m_livestemc_xfer_to_fire
    m_deadstemc_to_fire            => pcisof%m_deadstemc_to_fire
    m_deadstemc_storage_to_fire    => pcisof%m_deadstemc_storage_to_fire
    m_deadstemc_xfer_to_fire       => pcisof%m_deadstemc_xfer_to_fire
    m_frootc_to_fire               => pcisof%m_frootc_to_fire
    m_frootc_storage_to_fire       => pcisof%m_frootc_storage_to_fire
    m_frootc_xfer_to_fire          => pcisof%m_frootc_xfer_to_fire
    m_livecrootc_to_fire           => pcisof%m_livecrootc_to_fire
    m_livecrootc_storage_to_fire   => pcisof%m_livecrootc_storage_to_fire
    m_livecrootc_xfer_to_fire      => pcisof%m_livecrootc_xfer_to_fire
    m_deadcrootc_to_fire           => pcisof%m_deadcrootc_to_fire
    m_deadcrootc_storage_to_fire   => pcisof%m_deadcrootc_storage_to_fire
    m_deadcrootc_xfer_to_fire      => pcisof%m_deadcrootc_xfer_to_fire
    m_gresp_storage_to_fire        => pcisof%m_gresp_storage_to_fire
    m_gresp_xfer_to_fire           => pcisof%m_gresp_xfer_to_fire

    m_leafc_to_litter_fire               => pcisof%m_leafc_to_litter_fire
    m_leafc_storage_to_litter_fire      => pcisof%m_leafc_storage_to_litter_fire
    m_leafc_xfer_to_litter_fire         => pcisof%m_leafc_xfer_to_litter_fire
    m_livestemc_to_litter_fire          => pcisof%m_livestemc_to_litter_fire
    m_livestemc_storage_to_litter_fire  => pcisof%m_livestemc_storage_to_litter_fire
    m_livestemc_xfer_to_litter_fire     => pcisof%m_livestemc_xfer_to_litter_fire
    m_livestemc_to_deadstemc_fire       => pcisof%m_livestemc_to_deadstemc_fire
    m_deadstemc_to_litter_fire          => pcisof%m_deadstemc_to_litter_fire
    m_deadstemc_storage_to_litter_fire  => pcisof%m_deadstemc_storage_to_litter_fire
    m_deadstemc_xfer_to_litter_fire     => pcisof%m_deadstemc_xfer_to_litter_fire
    m_frootc_to_litter_fire             => pcisof%m_frootc_to_litter_fire
    m_frootc_storage_to_litter_fire     => pcisof%m_frootc_storage_to_litter_fire
    m_frootc_xfer_to_litter_fire        => pcisof%m_frootc_xfer_to_litter_fire
    m_livecrootc_to_litter_fire         => pcisof%m_livecrootc_to_litter_fire
    m_livecrootc_storage_to_litter_fire => pcisof%m_livecrootc_storage_to_litter_fire
    m_livecrootc_xfer_to_litter_fire    => pcisof%m_livecrootc_xfer_to_litter_fire
    m_livecrootc_to_deadcrootc_fire     => pcisof%m_livecrootc_to_deadcrootc_fire
    m_deadcrootc_to_litter_fire         => pcisof%m_deadcrootc_to_litter_fire
    m_deadcrootc_storage_to_litter_fire => pcisof%m_deadcrootc_storage_to_litter_fire
    m_deadcrootc_xfer_to_litter_fire    => pcisof%m_deadcrootc_xfer_to_litter_fire
    m_gresp_storage_to_litter_fire      => pcisof%m_gresp_storage_to_litter_fire
    m_gresp_xfer_to_litter_fire         => pcisof%m_gresp_xfer_to_litter_fire
    

    deadcrootc                     => pcisos%deadcrootc
    deadcrootc_storage             => pcisos%deadcrootc_storage
    deadcrootc_xfer                => pcisos%deadcrootc_xfer
    deadstemc                      => pcisos%deadstemc
    deadstemc_storage              => pcisos%deadstemc_storage
    deadstemc_xfer                 => pcisos%deadstemc_xfer
    frootc                         => pcisos%frootc
    frootc_storage                 => pcisos%frootc_storage
    frootc_xfer                    => pcisos%frootc_xfer
    gresp_storage                  => pcisos%gresp_storage
    gresp_xfer                     => pcisos%gresp_xfer
    leafc                          => pcisos%leafc
    leafc_storage                  => pcisos%leafc_storage
    leafc_xfer                     => pcisos%leafc_xfer
    livecrootc                     => pcisos%livecrootc
    livecrootc_storage             => pcisos%livecrootc_storage
    livecrootc_xfer                => pcisos%livecrootc_xfer
    livestemc                      => pcisos%livestemc
    livestemc_storage              => pcisos%livestemc_storage
    livestemc_xfer                 => pcisos%livestemc_xfer

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

end subroutine CStateUpdate3
!-----------------------------------------------------------------------
#endif

end module CNCStateUpdate3Mod
