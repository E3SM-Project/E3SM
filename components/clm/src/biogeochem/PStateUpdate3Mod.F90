module PStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for phosphorus state variable update, mortality fluxes.
  ! Also, sminn leaching flux.
  ! X.YANG
  ! !USES:
  use shr_kind_mod        , only: r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use clm_varpar          , only: nlevdecomp,ndecomp_pools,ndecomp_cascade_transitions
  use clm_time_manager    , only : get_step_size
  use clm_varctl          , only : iulog, use_nitrif_denitrif
  use clm_varpar          , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType         , only : cnstate_type
  use PhosphorusStateType , only : phosphorusstate_type
  use PhosphorusFLuxType  , only : phosphorusflux_type
  use soilorder_varcon    , only : smax,ks_sorption
  !! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode
  use clm_varctl          , only : nu_com
  use EcophysConType      , only : ecophyscon 
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: PStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine PStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars,phosphorusflux_vars, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic phosphorus state
    ! variables affected by gap-phase mortality fluxes. Also the Sminn leaching flux.
    ! NOTE - associate statements have been removed where there are
    ! no science equatiops. This increases readability and maintainability.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columps in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columps
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(phosphorusflux_type)  , intent(inout)    :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(cnstate_type)         , intent(in)    :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k        ! indices
    integer :: fp,fc      ! lake filter indices
    real(r8):: dt         ! radiation time step (seconds)

   real(r8):: smax_c       ! parameter(gP/m2), maximum amount of sorbed P in equilibrium with solution P
   real(r8):: ks_sorption_c ! parameter(gP/m2), empirical constant for sorbed P in equilibrium with solution P 
   real(r8):: flux_mineralization(bounds%begc:bounds%endc,1:nlevdecomp)   !! local temperary variable
   real(r8):: temp_solutionp(bounds%begc:bounds%endc,1:nlevdecomp)
   real(r8):: aa,bb,cc ! solve quadratic function

    !-----------------------------------------------------------------------

    associate(& 
         isoilorder     => cnstate_vars%isoilorder ,&
         pdep_prof      => cnstate_vars%pdep_prof_col ,&
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool ,&
         pf => phosphorusflux_vars  , &
         ps => phosphorusstate_vars , &
         vmax_minsurf_p_vr => ecophyscon%vmax_minsurf_p_vr , &
         km_minsurf_p_vr   => ecophyscon%km_minsurf_p_vr     &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      !! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes
      !! - X.YANG
      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            flux_mineralization(c,j) = 0._r8
         enddo
      enddo

      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                    flux_mineralization(c,j) = flux_mineralization(c,j) - &
                                               pf%decomp_cascade_sminp_flux_vr_col(c,j,k)*dt
               end do
            end do
         else
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                    flux_mineralization(c,j) = flux_mineralization(c,j) + &
                                               pf%decomp_cascade_sminp_flux_vr_col(c,j,k)*dt

               end do
            end do
         endif
      end do


      do j = 1, nlevdecomp
              ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            flux_mineralization(c,j) = flux_mineralization(c,j) + &
                                       pf%biochem_pmin_vr_col(c,j)

         end do
      end do

    if (nu_com .eq. 'RD') then
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
               ! assign read in parameter values
               smax_c = smax( isoilorder(c) )
               ks_sorption_c = ks_sorption( isoilorder(c) )
               temp_solutionp(c,j) = ps%solutionp_vr_col(c,j)
               ps%solutionp_vr_col(c,j)      = ps%solutionp_vr_col(c,j)  + ( flux_mineralization(c,j) + pf%primp_to_labilep_vr_col(c,j)*dt &
                    + pf%secondp_to_labilep_vr_col(c,j)*dt &
                    + pf%supplement_to_sminp_vr_col(c,j)*dt - pf%sminp_to_plant_vr_col(c,j)*dt&
                    - pf%labilep_to_secondp_vr_col(c,j)*dt - pf%sminp_leached_vr_col(c,j)*dt ) / &
                    (1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+ps%solutionp_vr_col(c,j))**2._r8)

               ps%labilep_vr_col(c,j) = ps%labilep_vr_col(c,j) + ((smax_c*ks_sorption_c)/(ks_sorption_c+temp_solutionp(c,j))**2._r8 ) * &
                    ( flux_mineralization(c,j) + pf%primp_to_labilep_vr_col(c,j)*dt + pf%secondp_to_labilep_vr_col(c,j)*dt &
                    + pf%supplement_to_sminp_vr_col(c,j)*dt - pf%sminp_to_plant_vr_col(c,j)*dt &
                    - pf%labilep_to_secondp_vr_col(c,j)*dt - pf%sminp_leached_vr_col(c,j)*dt ) / &
                    ( 1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+temp_solutionp(c,j))**2._r8 )
                             
               pf%desorb_to_solutionp_vr(c,j) = ( flux_mineralization(c,j)/dt + pf%primp_to_labilep_vr_col(c,j) &
                                + pf%secondp_to_labilep_vr_col(c,j) &
                                + pf%supplement_to_sminp_vr_col(c,j) - pf%sminp_to_plant_vr_col(c,j) &
                                - pf%labilep_to_secondp_vr_col(c,j) - pf%sminp_leached_vr_col(c,j) ) / &
                                (1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+ps%solutionp_vr_col(c,j))**2._r8)
            
                pf%adsorb_to_labilep_vr(c,j) = ((smax_c*ks_sorption_c)/(ks_sorption_c+temp_solutionp(c,j))**2._r8 ) * &
                             ( flux_mineralization(c,j)/dt + pf%primp_to_labilep_vr_col(c,j) + pf%secondp_to_labilep_vr_col(c,j) &
                             + pf%supplement_to_sminp_vr_col(c,j) - pf%sminp_to_plant_vr_col(c,j) &
                             - pf%labilep_to_secondp_vr_col(c,j) - pf%sminp_leached_vr_col(c,j) ) / &
                             ( 1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+temp_solutionp(c,j))**2._r8 )
            end do
         end do
      else ! ECA  
        do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               smax_c = vmax_minsurf_p_vr(isoilorder(c),j)
               ks_sorption_c = km_minsurf_p_vr(isoilorder(c),j)
               temp_solutionp(c,j) = ( ps%solutionp_vr_col(c,j) + ps%labilep_vr_col(c,j) + &
                            (flux_mineralization(c,j) + pf%primp_to_labilep_vr_col(c,j)*dt + &
                            pf%secondp_to_labilep_vr_col(c,j)*dt + pf%supplement_to_sminp_vr_col(c,j)*dt - &
                            pf%sminp_to_plant_vr_col(c,j)*dt - pf%labilep_to_secondp_vr_col(c,j)*dt - &
                            pf%sminp_leached_vr_col(c,j)*dt ))
                if (temp_solutionp(c,j) < 0.0_r8) then
                  pf%labilep_to_secondp_vr_col(c,j) = pf%labilep_to_secondp_vr_col(c,j)/(pf%labilep_to_secondp_vr_col(c,j)+pf%sminp_leached_vr_col(c,j))*(temp_solutionp(c,j) + pf%labilep_to_secondp_vr_col(c,j)*dt + pf%sminp_leached_vr_col(c,j)*dt) /dt
                  pf%sminp_leached_vr_col(c,j) = pf%sminp_leached_vr_col(c,j)/(pf%labilep_to_secondp_vr_col(c,j)+pf%sminp_leached_vr_col(c,j))*(temp_solutionp(c,j) + pf%labilep_to_secondp_vr_col(c,j)*dt + pf%sminp_leached_vr_col(c,j)*dt) /dt
                  temp_solutionp(c,j) = 0.0_r8
                  ps%solutionp_vr_col(c,j) = 0.0_r8
                  ps%labilep_vr_col(c,j) = 0.0_r8
                else
                  ! sorbp = smax*solutionp/(ks+solutionp)
                  ! sorbp + solutionp = smax*solutionp/(ks+solutionp) + solutionp = total p pool after competition
                  ! solve quadratic function to get equilibrium solutionp and adsorbp pools
                  aa = 1;
                  bb = smax_c + ks_sorption_c - temp_solutionp(c,j)
                  cc = -1.0 * ks_sorption_c *  temp_solutionp(c,j)
                  ps%solutionp_vr_col(c,j)  = (-bb+(bb*bb-4.0*aa*cc)**0.5)/(2.0*aa)
                  ps%labilep_vr_col(c,j) = temp_solutionp(c,j) - ps%solutionp_vr_col(c,j)
               end if
            enddo
         enddo
      end if
             
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            do l = 1, ndecomp_pools
                ps%decomp_ppools_vr_col(c,j,l) = ps%decomp_ppools_vr_col(c,j,l)- pf%biochem_pmin_ppools_vr_col(c,j,l)*dt
            end do
 
            ps%secondp_vr_col(c,j) = ps%secondp_vr_col(c,j) + ( pf%labilep_to_secondp_vr_col(c,j) - pf%secondp_to_labilep_vr_col(c,j) &
                                     - pf%secondp_to_occlp_vr_col(c,j) )*dt

            ps%occlp_vr_col(c,j)   = ps%occlp_vr_col(c,j) + ( pf%secondp_to_occlp_vr_col(c,j) ) * dt

            ps%primp_vr_col(c,j)   = ps%primp_vr_col(c,j) - ( pf%primp_to_labilep_vr_col(c,j) )*dt + pf%pdep_to_sminp_col(c)*dt * pdep_prof(c,j)
         end do
      enddo

      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! column level phosphorus fluxes from fire
            ! pft-level wood to column-level CWD (uncombusted wood)
            ps%decomp_ppools_vr_col(c,j,i_cwd) = ps%decomp_ppools_vr_col(c,j,i_cwd) + pf%fire_mortality_p_to_cwdp_col(c,j) * dt

            ! pft-level wood to column-level litter (uncombusted wood)
            ps%decomp_ppools_vr_col(c,j,i_met_lit) = ps%decomp_ppools_vr_col(c,j,i_met_lit) + pf%m_p_to_litr_met_fire_col(c,j)* dt
            ps%decomp_ppools_vr_col(c,j,i_cel_lit) = ps%decomp_ppools_vr_col(c,j,i_cel_lit) + pf%m_p_to_litr_cel_fire_col(c,j)* dt
            ps%decomp_ppools_vr_col(c,j,i_lig_lit) = ps%decomp_ppools_vr_col(c,j,i_lig_lit) + pf%m_p_to_litr_lig_fire_col(c,j)* dt
         end do ! end of column loop
      end do

      ! litter and CWD losses to fire
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               ps%decomp_ppools_vr_col(c,j,l) = ps%decomp_ppools_vr_col(c,j,l) - pf%m_decomp_ppools_to_fire_vr_col(c,j,l) * dt
            end do
         end do
      end do


      ! patch-level phosphorus fluxes 

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         !from fire displayed pools
         ps%leafp_patch(p)              =  ps%leafp_patch(p)      - pf%m_leafp_to_fire_patch(p)      * dt
         ps%frootp_patch(p)             =  ps%frootp_patch(p)     - pf%m_frootp_to_fire_patch(p)     * dt
         ps%livestemp_patch(p)          =  ps%livestemp_patch(p)  - pf%m_livestemp_to_fire_patch(p)  * dt
         ps%deadstemp_patch(p)          =  ps%deadstemp_patch(p)  - pf%m_deadstemp_to_fire_patch(p)  * dt
         ps%livecrootp_patch(p)         =  ps%livecrootp_patch(p) - pf%m_livecrootp_to_fire_patch(p) * dt
         ps%deadcrootp_patch(p)         =  ps%deadcrootp_patch(p) - pf%m_deadcrootp_to_fire_patch(p) * dt

         ps%leafp_patch(p)              =  ps%leafp_patch(p)      - pf%m_leafp_to_litter_fire_patch(p)           * dt
         ps%frootp_patch(p)             =  ps%frootp_patch(p)     - pf%m_frootp_to_litter_fire_patch(p)          * dt
         ps%livestemp_patch(p)          =  ps%livestemp_patch(p)  - pf%m_livestemp_to_litter_fire_patch(p)       * dt
         ps%deadstemp_patch(p)          =  ps%deadstemp_patch(p)  - pf%m_deadstemp_to_litter_fire_patch(p)       * dt
         ps%livecrootp_patch(p)         =  ps%livecrootp_patch(p) - pf%m_livecrootp_to_litter_fire_patch(p)      * dt
         ps%deadcrootp_patch(p)         =  ps%deadcrootp_patch(p) - pf%m_deadcrootp_to_litter_fire_patch(p)      * dt

         ! storage pools
         ps%leafp_storage_patch(p)      =  ps%leafp_storage_patch(p)      - pf%m_leafp_storage_to_fire_patch(p)      * dt
         ps%frootp_storage_patch(p)     =  ps%frootp_storage_patch(p)     - pf%m_frootp_storage_to_fire_patch(p)     * dt
         ps%livestemp_storage_patch(p)  =  ps%livestemp_storage_patch(p)  - pf%m_livestemp_storage_to_fire_patch(p)  * dt
         ps%deadstemp_storage_patch(p)  =  ps%deadstemp_storage_patch(p)  - pf%m_deadstemp_storage_to_fire_patch(p)  * dt
         ps%livecrootp_storage_patch(p) =  ps%livecrootp_storage_patch(p) - pf%m_livecrootp_storage_to_fire_patch(p) * dt
         ps%deadcrootp_storage_patch(p) =  ps%deadcrootp_storage_patch(p) - pf%m_deadcrootp_storage_to_fire_patch(p) * dt

         ps%leafp_storage_patch(p)      =  ps%leafp_storage_patch(p)      - pf%m_leafp_storage_to_litter_fire_patch(p)      * dt
         ps%frootp_storage_patch(p)     =  ps%frootp_storage_patch(p)     - pf%m_frootp_storage_to_litter_fire_patch(p)     * dt
         ps%livestemp_storage_patch(p)  =  ps%livestemp_storage_patch(p)  - pf%m_livestemp_storage_to_litter_fire_patch(p)  * dt
         ps%deadstemp_storage_patch(p)  =  ps%deadstemp_storage_patch(p)  - pf%m_deadstemp_storage_to_litter_fire_patch(p)  * dt
         ps%livecrootp_storage_patch(p) =  ps%livecrootp_storage_patch(p) - pf%m_livecrootp_storage_to_litter_fire_patch(p) * dt
         ps%deadcrootp_storage_patch(p) =  ps%deadcrootp_storage_patch(p) - pf%m_deadcrootp_storage_to_litter_fire_patch(p) * dt


         ! trapsfer pools
         ps%leafp_xfer_patch(p)         =  ps%leafp_xfer_patch(p)      - pf%m_leafp_xfer_to_fire_patch(p)      * dt
         ps%frootp_xfer_patch(p)        =  ps%frootp_xfer_patch(p)     - pf%m_frootp_xfer_to_fire_patch(p)     * dt
         ps%livestemp_xfer_patch(p)     =  ps%livestemp_xfer_patch(p)  - pf%m_livestemp_xfer_to_fire_patch(p)  * dt
         ps%deadstemp_xfer_patch(p)     =  ps%deadstemp_xfer_patch(p)  - pf%m_deadstemp_xfer_to_fire_patch(p)  * dt
         ps%livecrootp_xfer_patch(p)    =  ps%livecrootp_xfer_patch(p) - pf%m_livecrootp_xfer_to_fire_patch(p) * dt
         ps%deadcrootp_xfer_patch(p)    =  ps%deadcrootp_xfer_patch(p) - pf%m_deadcrootp_xfer_to_fire_patch(p) * dt

         ps%leafp_xfer_patch(p)         =  ps%leafp_xfer_patch(p)      - pf%m_leafp_xfer_to_litter_fire_patch(p)      * dt
         ps%frootp_xfer_patch(p)        =  ps%frootp_xfer_patch(p)     - pf%m_frootp_xfer_to_litter_fire_patch(p)     * dt
         ps%livestemp_xfer_patch(p)     =  ps%livestemp_xfer_patch(p)  - pf%m_livestemp_xfer_to_litter_fire_patch(p)  * dt
         ps%deadstemp_xfer_patch(p)     =  ps%deadstemp_xfer_patch(p)  - pf%m_deadstemp_xfer_to_litter_fire_patch(p)  * dt
         ps%livecrootp_xfer_patch(p)    =  ps%livecrootp_xfer_patch(p) - pf%m_livecrootp_xfer_to_litter_fire_patch(p) * dt
         ps%deadcrootp_xfer_patch(p)    =  ps%deadcrootp_xfer_patch(p) - pf%m_deadcrootp_xfer_to_litter_fire_patch(p) * dt

         ! retranslocated N pool
         ps%retransp_patch(p)           =  ps%retransp_patch(p) - pf%m_retransp_to_fire_patch(p)        * dt
         ps%retransp_patch(p)           =  ps%retransp_patch(p) - pf%m_retransp_to_litter_fire_patch(p) * dt
      end do

    end associate 

  end subroutine PStateUpdate3

end module PStateUpdate3Mod
