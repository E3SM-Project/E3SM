
module CNNStateUpdate3Mod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NStateUpdate3Mod
!
! !DESCRIPTION:
! Module for nitrogen state variable update, mortality fluxes.
! Also, sminn leaching flux.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar   , only: nlevdecomp, ndecomp_pools
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: NStateUpdate3
!
! !REVISION HISTORY:
! 7/27/2004: Created by Peter Thornton
! F. Li and S. Levis (11/06/12)
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NStateUpdate3
!
! !INTERFACE:
subroutine NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic nitrogen state
! variables affected by gap-phase mortality fluxes. Also the Sminn leaching flux.
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varctl  , only: iulog
   use clm_varpar   , only: i_cwd, i_met_lit, i_cel_lit, i_lig_lit
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
#ifndef NITRIF_DENITRIF
   real(r8), pointer :: sminn_leached_vr(:,:) 
#else
   real(r8), pointer :: smin_no3_leached_vr(:,:) 
   real(r8), pointer :: smin_no3_runoff_vr(:,:)     ! vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
   real(r8), pointer :: smin_no3_vr(:,:) 
   real(r8), pointer :: smin_nh4_vr(:,:) 
#endif
   real(r8), pointer :: m_leafn_to_fire(:)             
   real(r8), pointer :: m_leafn_storage_to_fire(:)     
   real(r8), pointer :: m_leafn_xfer_to_fire(:)       
   real(r8), pointer :: m_livestemn_to_fire(:)         
   real(r8), pointer :: m_livestemn_storage_to_fire(:) 
   real(r8), pointer :: m_livestemn_xfer_to_fire(:)    
   real(r8), pointer :: m_deadstemn_to_fire(:)         
   real(r8), pointer :: m_deadstemn_storage_to_fire(:) 
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:)    
   real(r8), pointer :: m_frootn_to_fire(:)            
   real(r8), pointer :: m_frootn_storage_to_fire(:)    
   real(r8), pointer :: m_frootn_xfer_to_fire(:)       
   real(r8), pointer :: m_livecrootn_to_fire(:)       
   real(r8), pointer :: m_livecrootn_storage_to_fire(:)
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)  
   real(r8), pointer :: m_deadcrootn_to_fire(:)  
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)   
   real(r8), pointer :: m_retransn_to_fire(:)          
   real(r8), pointer :: m_decomp_npools_to_fire_vr(:,:,:) 
 
   real(r8), pointer :: m_leafn_to_litter_fire(:)   
   real(r8), pointer :: m_leafn_storage_to_litter_fire(:)                
   real(r8), pointer :: m_leafn_xfer_to_litter_fire(:)  
   real(r8), pointer :: m_livestemn_to_litter_fire(:)    
   real(r8), pointer :: m_livestemn_storage_to_litter_fire(:)        
   real(r8), pointer :: m_livestemn_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_livestemn_to_deadstemn_fire(:)    
   real(r8), pointer :: m_deadstemn_to_litter_fire(:) 
   real(r8), pointer :: m_deadstemn_storage_to_litter_fire(:)           
   real(r8), pointer :: m_deadstemn_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_frootn_to_litter_fire(:)        
   real(r8), pointer :: m_frootn_storage_to_litter_fire(:)  
   real(r8), pointer :: m_frootn_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootn_to_litter_fire(:)    
   real(r8), pointer :: m_livecrootn_storage_to_litter_fire(:)      
   real(r8), pointer :: m_livecrootn_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootn_to_deadcrootn_fire(:)    
   real(r8), pointer :: m_deadcrootn_to_litter_fire(:)        
   real(r8), pointer :: m_deadcrootn_storage_to_litter_fire(:)  
   real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire(:)
   real(r8), pointer :: m_retransn_to_litter_fire(:)      
   real(r8), pointer :: m_n_to_litr_met_fire(:,:)
   real(r8), pointer :: m_n_to_litr_cel_fire(:,:)
   real(r8), pointer :: m_n_to_litr_lig_fire(:,:)
   real(r8), pointer :: fire_mortality_n_to_cwdn(:,:)              ! N fluxes associated with fire mortality to CWD pool (gN/m3/s)
!
! local pointers to implicit in/out scalars
   real(r8), pointer :: sminn_vr(:,:)              ! (gN/m3) soil mineral N
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,j,l,k        ! indices
   integer :: fp,fc      ! lake filter indices
   real(r8):: dt         ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------

    ! assign local pointers at the column level
    fire_mortality_n_to_cwdn       => clm3%g%l%c%cnf%fire_mortality_n_to_cwdn
#ifndef NITRIF_DENITRIF
    sminn_leached_vr                  => clm3%g%l%c%cnf%sminn_leached_vr
#else
    smin_no3_leached_vr               => clm3%g%l%c%cnf%smin_no3_leached_vr
    smin_no3_runoff_vr                => clm3%g%l%c%cnf%smin_no3_runoff_vr
    smin_no3_vr                       => clm3%g%l%c%cns%smin_no3_vr
    smin_nh4_vr                       => clm3%g%l%c%cns%smin_nh4_vr
#endif
    m_decomp_npools_to_fire_vr            => clm3%g%l%c%cnf%m_decomp_npools_to_fire_vr
    m_n_to_litr_met_fire                  => clm3%g%l%c%cnf%m_n_to_litr_met_fire
    m_n_to_litr_cel_fire                  => clm3%g%l%c%cnf%m_n_to_litr_cel_fire
    m_n_to_litr_lig_fire                  => clm3%g%l%c%cnf%m_n_to_litr_lig_fire
   
    decomp_npools_vr                           => clm3%g%l%c%cns%decomp_npools_vr
    sminn_vr                          => clm3%g%l%c%cns%sminn_vr

    ! assign local pointers at the pft level
   m_leafn_to_fire                => clm3%g%l%c%p%pnf%m_leafn_to_fire
    m_leafn_storage_to_fire        => clm3%g%l%c%p%pnf%m_leafn_storage_to_fire
    m_leafn_xfer_to_fire           => clm3%g%l%c%p%pnf%m_leafn_xfer_to_fire
    m_livestemn_to_fire            => clm3%g%l%c%p%pnf%m_livestemn_to_fire
    m_livestemn_storage_to_fire    => clm3%g%l%c%p%pnf%m_livestemn_storage_to_fire
    m_livestemn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_fire
    m_deadstemn_to_fire            => clm3%g%l%c%p%pnf%m_deadstemn_to_fire
    m_deadstemn_storage_to_fire    => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_fire
    m_deadstemn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_fire
    m_frootn_to_fire               => clm3%g%l%c%p%pnf%m_frootn_to_fire
    m_frootn_storage_to_fire       => clm3%g%l%c%p%pnf%m_frootn_storage_to_fire
    m_frootn_xfer_to_fire          => clm3%g%l%c%p%pnf%m_frootn_xfer_to_fire
    m_livecrootn_to_fire           => clm3%g%l%c%p%pnf%m_livecrootn_to_fire
    m_livecrootn_storage_to_fire   => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_fire
    m_livecrootn_xfer_to_fire      => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_fire
    m_deadcrootn_to_fire           => clm3%g%l%c%p%pnf%m_deadcrootn_to_fire
    m_deadcrootn_storage_to_fire   => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_fire
    m_deadcrootn_xfer_to_fire      => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_fire
    m_retransn_to_fire             => clm3%g%l%c%p%pnf%m_retransn_to_fire
    m_leafn_to_litter_fire                => clm3%g%l%c%p%pnf%m_leafn_to_litter_fire
    m_leafn_storage_to_litter_fire        => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter_fire
    m_leafn_xfer_to_litter_fire           => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter_fire
    m_livestemn_to_litter_fire            => clm3%g%l%c%p%pnf%m_livestemn_to_litter_fire
    m_livestemn_storage_to_litter_fire    => clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter_fire
    m_livestemn_xfer_to_litter_fire       => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter_fire
    m_livestemn_to_deadstemn_fire         => clm3%g%l%c%p%pnf%m_livestemn_to_deadstemn_fire
    m_deadstemn_to_litter_fire            => clm3%g%l%c%p%pnf%m_deadstemn_to_litter_fire
    m_deadstemn_storage_to_litter_fire    => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter_fire
    m_deadstemn_xfer_to_litter_fire       =>clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter_fire
    m_frootn_to_litter_fire               => clm3%g%l%c%p%pnf%m_frootn_to_litter_fire
    m_frootn_storage_to_litter_fire       => clm3%g%l%c%p%pnf%m_frootn_storage_to_litter_fire
    m_frootn_xfer_to_litter_fire          => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter_fire
    m_livecrootn_to_litter_fire           => clm3%g%l%c%p%pnf%m_livecrootn_to_litter_fire
    m_livecrootn_storage_to_litter_fire   => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter_fire
    m_livecrootn_xfer_to_litter_fire      => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter_fire
    m_livecrootn_to_deadcrootn_fire       => clm3%g%l%c%p%pnf%m_livecrootn_to_deadcrootn_fire
    m_deadcrootn_to_litter_fire           => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter_fire
    m_deadcrootn_storage_to_litter_fire   => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter_fire
    m_deadcrootn_xfer_to_litter_fire      => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter_fire
    m_retransn_to_litter_fire             => clm3%g%l%c%p%pnf%m_retransn_to_litter_fire

    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn

   ! set time steps
   dt = real( get_step_size(), r8 )

   do j = 1, nlevdecomp
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)

#ifndef NITRIF_DENITRIF         
         ! mineral N loss due to leaching
         sminn_vr(c,j) = sminn_vr(c,j) - sminn_leached_vr(c,j) * dt
#else
         ! mineral N loss due to leaching and runoff
         smin_no3_vr(c,j) = max(smin_no3_vr(c,j) - ( smin_no3_leached_vr(c,j) + smin_no3_runoff_vr(c,j) ) * dt, 0._r8)
         sminn_vr(c,j) = smin_no3_vr(c,j) + smin_nh4_vr(c,j)
#endif         
         
         ! column level nitrogen fluxes from fire
         ! pft-level wood to column-level CWD (uncombusted wood)
         decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + fire_mortality_n_to_cwdn(c,j) * dt
         
         ! pft-level wood to column-level litter (uncombusted wood)
         decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + m_n_to_litr_met_fire(c,j)* dt
         decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + m_n_to_litr_cel_fire(c,j)* dt
         decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + m_n_to_litr_lig_fire(c,j)* dt
      end do ! end of column loop
   end do
   
   ! litter and CWD losses to fire
   do l = 1, ndecomp_pools
      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_npools_vr(c,j,l) = decomp_npools_vr(c,j,l) - m_decomp_npools_to_fire_vr(c,j,l) * dt
         end do
      end do
   end do

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from fire
      ! displayed pools
      leafn(p)      = leafn(p)      - m_leafn_to_fire(p)             * dt
      frootn(p)     = frootn(p)     - m_frootn_to_fire(p)            * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_fire(p)         * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_fire(p)         * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_fire(p)        * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_fire(p)        * dt
      
      leafn(p)      = leafn(p)      - m_leafn_to_litter_fire(p)           * dt
      frootn(p)     = frootn(p)     - m_frootn_to_litter_fire(p)          * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_litter_fire(p)       * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_litter_fire(p)       * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_litter_fire(p)      * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_litter_fire(p)      * dt
     
      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - m_leafn_storage_to_fire(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - m_frootn_storage_to_fire(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - m_livestemn_storage_to_fire(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - m_deadstemn_storage_to_fire(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - m_livecrootn_storage_to_fire(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - m_deadcrootn_storage_to_fire(p) * dt

      leafn_storage(p)      = leafn_storage(p)      - m_leafn_storage_to_litter_fire(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - m_frootn_storage_to_litter_fire(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - m_livestemn_storage_to_litter_fire(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - m_deadstemn_storage_to_litter_fire(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - m_livecrootn_storage_to_litter_fire(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - m_deadcrootn_storage_to_litter_fire(p) * dt


      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - m_leafn_xfer_to_fire(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - m_frootn_xfer_to_fire(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - m_livestemn_xfer_to_fire(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - m_deadstemn_xfer_to_fire(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - m_livecrootn_xfer_to_fire(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - m_deadcrootn_xfer_to_fire(p) * dt
  
      leafn_xfer(p)      = leafn_xfer(p)      - m_leafn_xfer_to_litter_fire(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - m_frootn_xfer_to_litter_fire(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - m_livestemn_xfer_to_litter_fire(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - m_deadstemn_xfer_to_litter_fire(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - m_livecrootn_xfer_to_litter_fire(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - m_deadcrootn_xfer_to_litter_fire(p) * dt

      ! retranslocated N pool
      retransn(p) = retransn(p) - m_retransn_to_fire(p)        * dt
      retransn(p) = retransn(p) - m_retransn_to_litter_fire(p) * dt

   end do

end subroutine NStateUpdate3
!-----------------------------------------------------------------------
#endif

end module CNNStateUpdate3Mod
