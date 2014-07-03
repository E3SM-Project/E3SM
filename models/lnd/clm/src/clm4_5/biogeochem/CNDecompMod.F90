module CNDecompMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDecompMod
!
! !DESCRIPTION:
! Module holding routines used in litter and soil decomposition model
! for coupled carbon-nitrogen code.
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_TKFRZ
    use clm_varctl  , only: iulog
    use clm_varcon, only: dzsoi_decomp
#ifndef CENTURY_DECOMP
    use CNDecompCascadeMod_BGC, only : decomp_rate_constants
#else
    use CNDecompCascadeMod_CENTURY, only : decomp_rate_constants
#endif
#ifdef NITRIF_DENITRIF
    use CNNitrifDenitrifMod, only: nitrif_denitrif
#endif
    use CNVerticalProfileMod, only: decomp_vertprofiles

   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public:: CNDecompAlloc
   
!
! !REVISION HISTORY:
! 8/15/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated to vector data structures
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNDecompAlloc
!
! !INTERFACE:
subroutine CNDecompAlloc (lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
   num_soilp, filter_soilp)
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use CNAllocationMod , only: CNAllocation
   use clm_time_manager, only: get_step_size
   use clm_varpar   , only: nlevsoi,nlevgrnd,nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
   use pft2colMod      , only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbp, ubp        ! pft-index bounds
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 8/15/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! all c pools involved in decomposition
   real(r8), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   ! all n pools involved in decomposition
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools

   integer, pointer :: clandunit(:)      ! index into landunit level quantities
   integer , pointer :: itypelun(:)      ! landunit type
   ! pft level
   real(r8), pointer :: rootfr(:,:)      ! fraction of roots in each soil layer  (nlevgrnd)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: fpi_vr(:,:)                            ! fraction of potential immobilization (no units)
   real(r8), pointer :: decomp_cascade_hr_vr(:,:,:)            ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: decomp_cascade_ctransfer_vr(:,:,:)     ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: decomp_pools_hr(:,:)                   ! het. resp. from decomposing C pools (gC/m2/s)
   real(r8), pointer :: decomp_cascade_ctransfer(:,:)          ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
   real(r8), pointer :: decomp_cascade_ntransfer_vr(:,:,:)     ! vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
   real(r8), pointer :: decomp_cascade_sminn_flux_vr(:,:,:)    ! vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)

   real(r8), pointer :: potential_immob_vr(:,:)
#ifndef NITRIF_DENITRIF
   real(r8), pointer :: sminn_to_denit_decomp_cascade_vr(:,:,:)
   real(r8), pointer :: sminn_to_denit_excess_vr(:,:)
#endif
   real(r8), pointer :: gross_nmin_vr(:,:)
   real(r8), pointer :: net_nmin_vr(:,:)
   real(r8), pointer :: gross_nmin(:)            ! gross rate of N mineralization (gN/m2/s)
   real(r8), pointer :: net_nmin(:)              ! net rate of N mineralization (gN/m2/s)
   ! For methane code
#ifdef LCH4
   real(r8), pointer :: fphr(:,:)                ! fraction of potential SOM + LITTER heterotrophic respiration
   real(r8), pointer :: w_scalar(:,:)            ! fraction by which decomposition is limited by moisture availability
#endif
   real(r8), pointer :: decomp_k(:,:,:)                       ! rate constant for decomposition (1./sec)
   real(r8), pointer :: rf_decomp_cascade(:,:,:)              ! respired fraction in decomposition step (frac)
   integer,  pointer :: cascade_donor_pool(:)                 ! which pool is C taken from for a given decomposition step
   integer,  pointer :: cascade_receiver_pool(:)              ! which pool is C added to for a given decomposition step
   real(r8), pointer :: pathfrac_decomp_cascade(:,:,:)        ! what fraction of C leaving a given pool passes through a given transition (frac)
   logical,  pointer :: floating_cn_ratio_decomp_pools(:)     ! TRUE => pool has fixed C:N ratio

!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,j,k,l,m          !indices
   integer :: fc           !lake filter column index
   real(r8):: p_decomp_cpool_loss(lbc:ubc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
   real(r8):: pmnf_decomp_cascade(lbc:ubc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
   real(r8):: immob(lbc:ubc,1:nlevdecomp)        !potential N immobilization
   real(r8):: ratio        !temporary variable
   real(r8):: dnp          !denitrification proportion
   real(r8):: cn_decomp_pools(lbc:ubc,1:nlevdecomp,1:ndecomp_pools)
   real(r8), pointer :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
   integer, parameter :: i_atm = 0
   integer, pointer :: altmax_indx(:)                  ! maximum annual depth of thaw
   integer, pointer :: altmax_lastyear_indx(:)         ! prior year maximum annual depth of thaw


   ! For methane code
#ifndef NITRIF_DENITRIF
   real(r8):: phr_vr(lbc:ubc,1:nlevdecomp)       !potential HR (gC/m3/s)
#else
   real(r8), pointer :: phr_vr(:,:)              !potential HR (gC/m3/s)
#endif
   real(r8):: hrsum(lbc:ubc,1:nlevdecomp)        !sum of HR (gC/m2/s)
   
   !EOP
   !-----------------------------------------------------------------------
   
   decomp_cpools_vr              => clm3%g%l%c%ccs%decomp_cpools_vr
   decomp_cascade_hr_vr          => clm3%g%l%c%ccf%decomp_cascade_hr_vr
   decomp_cascade_ctransfer_vr   => clm3%g%l%c%ccf%decomp_cascade_ctransfer_vr
   decomp_npools_vr              => clm3%g%l%c%cns%decomp_npools_vr
   decomp_cascade_ntransfer_vr   => clm3%g%l%c%cnf%decomp_cascade_ntransfer_vr
   decomp_cascade_sminn_flux_vr  => clm3%g%l%c%cnf%decomp_cascade_sminn_flux_vr
   fpi_vr                => clm3%g%l%c%cps%fpi_vr
   potential_immob_vr    => clm3%g%l%c%cnf%potential_immob_vr
   
   decomp_k                        => clm3%g%l%c%ccf%decomp_k
   rf_decomp_cascade               => clm3%g%l%c%cps%rf_decomp_cascade
   cascade_donor_pool              => decomp_cascade_con%cascade_donor_pool
   cascade_receiver_pool           => decomp_cascade_con%cascade_receiver_pool
   pathfrac_decomp_cascade         => clm3%g%l%c%cps%pathfrac_decomp_cascade
   floating_cn_ratio_decomp_pools  => decomp_cascade_con%floating_cn_ratio_decomp_pools
   initial_cn_ratio                => decomp_cascade_con%initial_cn_ratio
   altmax_indx                     => clm3%g%l%c%cps%altmax_indx
   altmax_lastyear_indx            => clm3%g%l%c%cps%altmax_lastyear_indx
   
#ifndef NITRIF_DENITRIF
   sminn_to_denit_decomp_cascade_vr   => clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade_vr
   sminn_to_denit_excess_vr => clm3%g%l%c%cnf%sminn_to_denit_excess_vr
#else
   phr_vr                   => clm3%g%l%c%ccf%phr_vr
#endif
   gross_nmin_vr            => clm3%g%l%c%cnf%gross_nmin_vr
   net_nmin_vr              => clm3%g%l%c%cnf%net_nmin_vr
   gross_nmin               => clm3%g%l%c%cnf%gross_nmin
   net_nmin                 => clm3%g%l%c%cnf%net_nmin
   ! For methane code
#ifdef LCH4
   fphr               => clm3%g%l%c%cch4%fphr
   w_scalar           => clm3%g%l%c%ccf%w_scalar
#endif
   
   rootfr                => clm3%g%l%c%p%pps%rootfr
   clandunit             => clm3%g%l%c%landunit
   itypelun              => clm3%g%l%itype
   
   
   
   call decomp_rate_constants(lbc, ubc, num_soilc, filter_soilc)
   
   
   ! set initial values for potential C and N fluxes
   p_decomp_cpool_loss(:,:,:) = 0._r8
   pmnf_decomp_cascade(:,:,:) = 0._r8
   
   ! column loop to calculate potential decomp rates and total immobilization
   ! demand.
   
   
!!! calculate c:n ratios of applicable pools
   do l = 1, ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( decomp_npools_vr(c,j,l) .gt. 0._r8 ) then
                  cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
               end if
            end do
         end do
      else
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
            end do
         end do
      end if
   end do
   
   
   ! calculate the non-nitrogen-limited fluxes
   ! these fluxes include the  "/ dt" term to put them on a
   ! per second basis, since the rate constants have been
   ! calculated on a per timestep basis.
   
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. decomp_k(c,j,cascade_donor_pool(k)) .gt. 0._r8 ) then
               p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
               
               if ( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then  !! not transition of cwd to litter
                  
                  if (cascade_receiver_pool(k) .ne. i_atm ) then  ! not 100% respiration
                     ratio = 0._r8
                     
                     if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8) then
                        ratio = cn_decomp_pools(c,j,cascade_receiver_pool(k))/cn_decomp_pools(c,j,cascade_donor_pool(k))
                     endif
                     
                     pmnf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(c,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k) - ratio) &
                          / cn_decomp_pools(c,j,cascade_receiver_pool(k)) )
                     
                  else   ! 100% respiration
                     pmnf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
                  endif
                  
               else   ! CWD -> litter
                  pmnf_decomp_cascade(c,j,k) = 0._r8
               end if
            end if
         end do
         
      end do
   end do
   
   ! Sum up all the potential immobilization fluxes (positive pmnf flux)
   ! and all the mineralization fluxes (negative pmnf flux)
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         immob(c,j) = 0._r8
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (pmnf_decomp_cascade(c,j,k) > 0._r8) then
               immob(c,j) = immob(c,j) + pmnf_decomp_cascade(c,j,k)
            else
               gross_nmin_vr(c,j) = gross_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
            end if
         end do
      end do
   end do
   
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         potential_immob_vr(c,j) = immob(c,j)
      end do
   end do
   
   ! Add up potential hr for methane calculations
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         phr_vr(c,j) = 0._r8
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            phr_vr(c,j) = phr_vr(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
         end do
      end do
   end do
   
   call decomp_vertprofiles(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp,filter_soilp)
   
#ifdef NITRIF_DENITRIF
   ! calculate nitrification and denitrification rates
   call nitrif_denitrif(lbc, ubc, num_soilc, filter_soilc)
#endif
   
   
   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.
   
   call CNAllocation(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp, &
        filter_soilp)
   
   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N
   
   dnp = 0.01_r8
   
   ! calculate c:n ratios of applicable pools
   do l = 1, ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( decomp_npools_vr(c,j,l) .gt. 0._r8 ) then
                  cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
               end if
            end do
         end do
      else
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
            end do
         end do
      end if
   end do
   
   ! upon return from CNAllocation, the fraction of potential immobilization
   ! has been set (cps%fpi_vr). now finish the decomp calculations.
   ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
   ! Also calculate denitrification losses as a simple proportion
   ! of mineralization flux.
   
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) .gt. 0._r8) then
               if ( pmnf_decomp_cascade(c,j,k) .gt. 0._r8 ) then
                  p_decomp_cpool_loss(c,j,k) = p_decomp_cpool_loss(c,j,k) * fpi_vr(c,j)
                  pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
#ifndef NITRIF_DENITRIF
                  sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
               else
                  sminn_to_denit_decomp_cascade_vr(c,j,k) = -dnp * pmnf_decomp_cascade(c,j,k)
#endif
               end if
               decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
               decomp_cascade_ctransfer_vr(c,j,k) = (1._r8 - rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
               if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. cascade_receiver_pool(k) .ne. i_atm) then
                  decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
               else
                  decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
               endif
               if ( cascade_receiver_pool(k) .ne. 0 ) then
                  decomp_cascade_sminn_flux_vr(c,j,k) = pmnf_decomp_cascade(c,j,k)
               else  ! keep sign convention negative for terminal pools
                  decomp_cascade_sminn_flux_vr(c,j,k) = - pmnf_decomp_cascade(c,j,k)
               endif
               net_nmin_vr(c,j) = net_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
            else
               decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
#ifndef NITRIF_DENITRIF
               sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
#endif
               decomp_cascade_sminn_flux_vr(c,j,k) = 0._r8
            end if
            
         end do
      end do
   end do
   
#ifdef LCH4
   ! Calculate total fraction of potential HR, for methane code
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         hrsum(c,j) = 0._r8
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            hrsum(c,j) = hrsum(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
         end do
      end do
   end do
   
   ! Nitrogen limitation / (low)-moisture limitation
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         if (phr_vr(c,j) > 0._r8) then
            fphr(c,j) = hrsum(c,j) / phr_vr(c,j) * w_scalar(c,j)
            fphr(c,j) = max(fphr(c,j), 0.01_r8) ! Prevent overflow errors for 0 respiration
         else
            fphr(c,j) = 1._r8
         end if
      end do
   end do
#endif
   
   ! vertically integrate net and gross mineralization fluxes for diagnostic output
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
         gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)   
      end do
   end do
   
 end subroutine CNDecompAlloc
 
 
#endif
 
end module CNDecompMod
