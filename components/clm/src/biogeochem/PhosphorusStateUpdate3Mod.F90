module PhosphorusStateUpdate3Mod

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
  use tracer_varcon       , only : is_active_betr_bgc
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode
  use clm_varctl          , only : nu_com
  use clm_varctl          , only : ECA_Pconst_RGspin
  use VegetationPropertiesType      , only : veg_vp 
  use ColumnDataType      , only : col_ps
  use VegetationDataType  , only : veg_ps
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: PhosphorusStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine PhosphorusStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
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
         vmax_minsurf_p_vr => veg_vp%vmax_minsurf_p_vr , &
         km_minsurf_p_vr   => veg_vp%km_minsurf_p_vr     &
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
      if(is_active_betr_bgc)then
        do j = 1, nlevdecomp
          do fc = 1,num_soilc
            c = filter_soilc(fc)
            col_ps%primp_vr(c,j)   = col_ps%primp_vr(c,j) - pf%primp_to_labilep_vr_col(c,j) *dt &
                 + pf%pdep_to_sminp_col(c)*dt * pdep_prof(c,j)
          end do
        enddo
      else
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
                                       pf%biochem_pmin_vr_col(c,j)*dt
           end do
        end do

      if (nu_com .eq. 'RD') then
        do j = 1, nlevdecomp
          do fc = 1,num_soilc
            c = filter_soilc(fc)
               ! assign read in parameter values
               smax_c = smax( isoilorder(c) )
               ks_sorption_c = ks_sorption( isoilorder(c) )
               temp_solutionp(c,j) = col_ps%solutionp_vr(c,j)
               col_ps%solutionp_vr(c,j)      = col_ps%solutionp_vr(c,j)  + ( flux_mineralization(c,j) &
                    + pf%primp_to_labilep_vr_col(c,j)*dt &
                    + pf%secondp_to_labilep_vr_col(c,j)*dt &
                    + pf%supplement_to_sminp_vr_col(c,j)*dt - pf%sminp_to_plant_vr_col(c,j)*dt&
                    - pf%labilep_to_secondp_vr_col(c,j)*dt - pf%sminp_leached_vr_col(c,j)*dt ) / &
                    (1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+col_ps%solutionp_vr(c,j))**2._r8)

               col_ps%labilep_vr(c,j) = col_ps%labilep_vr(c,j) + ((smax_c*ks_sorption_c)&
                    /(ks_sorption_c+temp_solutionp(c,j))**2._r8 ) * &
                    ( flux_mineralization(c,j) + pf%primp_to_labilep_vr_col(c,j)*dt + pf%secondp_to_labilep_vr_col(c,j)*dt &
                    + pf%supplement_to_sminp_vr_col(c,j)*dt - pf%sminp_to_plant_vr_col(c,j)*dt &
                    - pf%labilep_to_secondp_vr_col(c,j)*dt - pf%sminp_leached_vr_col(c,j)*dt ) / &
                    ( 1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+temp_solutionp(c,j))**2._r8 )
                             
               pf%desorb_to_solutionp_vr(c,j) = ( flux_mineralization(c,j)/dt + pf%primp_to_labilep_vr_col(c,j) &
                                + pf%secondp_to_labilep_vr_col(c,j) &
                                + pf%supplement_to_sminp_vr_col(c,j) - pf%sminp_to_plant_vr_col(c,j) &
                                - pf%labilep_to_secondp_vr_col(c,j) - pf%sminp_leached_vr_col(c,j) ) / &
                                (1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+col_ps%solutionp_vr(c,j))**2._r8)
            
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
                col_ps%solutionp_vr_prev(c,j) = col_ps%solutionp_vr(c,j)
                col_ps%labilep_vr_prev(c,j) = col_ps%labilep_vr(c,j)

                smax_c = vmax_minsurf_p_vr(isoilorder(c),j)
                ks_sorption_c = km_minsurf_p_vr(isoilorder(c),j)
                temp_solutionp(c,j) = ( col_ps%solutionp_vr(c,j) + col_ps%labilep_vr(c,j) + &
                            (flux_mineralization(c,j) + pf%primp_to_labilep_vr_col(c,j)*dt + &
                            pf%secondp_to_labilep_vr_col(c,j)*dt + pf%supplement_to_sminp_vr_col(c,j)*dt - &
                            pf%sminp_to_plant_vr_col(c,j)*dt - pf%labilep_to_secondp_vr_col(c,j)*dt - &
                            pf%sminp_leached_vr_col(c,j)*dt ))

                 if (temp_solutionp(c,j) < 0.0_r8) then
                    pf%labilep_to_secondp_vr_col(c,j) = pf%labilep_to_secondp_vr_col(c,j)/ &
                            (pf%labilep_to_secondp_vr_col(c,j)+pf%sminp_leached_vr_col(c,j))* &
                            (temp_solutionp(c,j) + pf%labilep_to_secondp_vr_col(c,j)*dt + &
                            pf%sminp_leached_vr_col(c,j)*dt) /dt
                    pf%sminp_leached_vr_col(c,j) = pf%sminp_leached_vr_col(c,j)/ &
                            (pf%labilep_to_secondp_vr_col(c,j)+pf%sminp_leached_vr_col(c,j))* &
                            (temp_solutionp(c,j) + pf%labilep_to_secondp_vr_col(c,j)*dt + &
                            pf%sminp_leached_vr_col(c,j)*dt) /dt
                       temp_solutionp(c,j) = 0.0_r8
                       col_ps%solutionp_vr(c,j) = 0.0_r8
                       col_ps%labilep_vr(c,j) = 0.0_r8
                 else
                       ! sorbp = smax*solutionp/(ks+solutionp)
                       ! sorbp + solutionp = smax*solutionp/(ks+solutionp) + solutionp = total p pool after competition
                       ! solve quadratic function to get equilibrium solutionp and adsorbp pools
                       aa = 1;
                       bb = smax_c + ks_sorption_c - temp_solutionp(c,j)
                       cc = -1.0_r8 * ks_sorption_c *  temp_solutionp(c,j)
                       col_ps%solutionp_vr(c,j)  = (-bb+(bb*bb-4.0_r8*aa*cc)**0.5_r8)/(2.0_r8*aa)
                       col_ps%labilep_vr(c,j) = temp_solutionp(c,j) - col_ps%solutionp_vr(c,j)
                 end if

                 col_ps%solutionp_vr_cur(c,j) = col_ps%solutionp_vr(c,j)
                 col_ps%labilep_vr_cur(c,j) = col_ps%labilep_vr(c,j)
              enddo
           enddo
         end if

         if (nu_com .eq. 'RD') then
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                do l = 1, ndecomp_pools

                   col_ps%decomp_ppools_vr(c,j,l) = col_ps%decomp_ppools_vr(c,j,l)- pf%biochem_pmin_ppools_vr_col(c,j,l)*dt

                end do
             end do
          end do
         end if

         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               col_ps%secondp_vr_prev(c,j) = col_ps%secondp_vr(c,j)
               col_ps%occlp_vr_prev(c,j)   = col_ps%occlp_vr(c,j)
               col_ps%primp_vr_prev(c,j)   = col_ps%primp_vr(c,j)

               col_ps%secondp_vr(c,j) = col_ps%secondp_vr(c,j) + ( pf%labilep_to_secondp_vr_col(c,j) &
                    - pf%secondp_to_labilep_vr_col(c,j) &
                                     - pf%secondp_to_occlp_vr_col(c,j) )*dt
               col_ps%occlp_vr(c,j)   = col_ps%occlp_vr(c,j) + ( pf%secondp_to_occlp_vr_col(c,j) ) * dt
               col_ps%primp_vr(c,j)   = col_ps%primp_vr(c,j) - ( pf%primp_to_labilep_vr_col(c,j) )*dt + pf%pdep_to_sminp_col(c)*dt &
                    * pdep_prof(c,j)

               col_ps%secondp_vr_cur(c,j) = col_ps%secondp_vr(c,j)
               col_ps%occlp_vr_cur(c,j)   = col_ps%occlp_vr(c,j)
               col_ps%primp_vr_cur(c,j)   = col_ps%primp_vr(c,j)
            end do
         enddo
         
         ! phosphorus pools do not change during RG spinup, but fluxes are still calculated to drive soil/plant P cycles
         ! rationale: observed P pools should be our best representation of present-day soil P conditions
         ! If we use observed P to initialize regular spinup, soil P pools will dramatically deplete during the spinup
         ! Then, the transient simulation will start with a much lower soil phosphorus availability that is inconsistent with obs 
         if ((nu_com .ne. 'RD') .and. ECA_Pconst_RGspin ) then
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  col_ps%solutionp_vr(c,j) = col_ps%solutionp_vr_prev(c,j)
                  col_ps%labilep_vr(c,j) = col_ps%labilep_vr_prev(c,j)
                  col_ps%secondp_vr(c,j) = col_ps%secondp_vr_prev(c,j)
                  col_ps%occlp_vr(c,j)   = col_ps%occlp_vr_prev(c,j)
                  col_ps%primp_vr(c,j)   = col_ps%primp_vr_prev(c,j)
               end do
            end do
         end if

      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! column level phosphorus fluxes from fire
            ! pft-level wood to column-level CWD (uncombusted wood)
            col_ps%decomp_ppools_vr(c,j,i_cwd) = col_ps%decomp_ppools_vr(c,j,i_cwd) + pf%fire_mortality_p_to_cwdp_col(c,j) * dt

            ! pft-level wood to column-level litter (uncombusted wood)
            col_ps%decomp_ppools_vr(c,j,i_met_lit) = col_ps%decomp_ppools_vr(c,j,i_met_lit) + pf%m_p_to_litr_met_fire_col(c,j)* dt
            col_ps%decomp_ppools_vr(c,j,i_cel_lit) = col_ps%decomp_ppools_vr(c,j,i_cel_lit) + pf%m_p_to_litr_cel_fire_col(c,j)* dt
            col_ps%decomp_ppools_vr(c,j,i_lig_lit) = col_ps%decomp_ppools_vr(c,j,i_lig_lit) + pf%m_p_to_litr_lig_fire_col(c,j)* dt
         end do ! end of column loop
      end do

      ! litter and CWD losses to fire
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               col_ps%decomp_ppools_vr(c,j,l) = col_ps%decomp_ppools_vr(c,j,l) - pf%m_decomp_ppools_to_fire_vr_col(c,j,l) * dt
            end do
         end do
      end do

    endif !is_active_betr_bgc
      ! patch-level phosphorus fluxes 

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         !from fire displayed pools
         veg_ps%leafp(p)              =  veg_ps%leafp(p)      - pf%m_leafp_to_fire_patch(p)      * dt
         veg_ps%frootp(p)             =  veg_ps%frootp(p)     - pf%m_frootp_to_fire_patch(p)     * dt
         veg_ps%livestemp(p)          =  veg_ps%livestemp(p)  - pf%m_livestemp_to_fire_patch(p)  * dt
         veg_ps%deadstemp(p)          =  veg_ps%deadstemp(p)  - pf%m_deadstemp_to_fire_patch(p)  * dt
         veg_ps%livecrootp(p)         =  veg_ps%livecrootp(p) - pf%m_livecrootp_to_fire_patch(p) * dt
         veg_ps%deadcrootp(p)         =  veg_ps%deadcrootp(p) - pf%m_deadcrootp_to_fire_patch(p) * dt

         veg_ps%leafp(p)              =  veg_ps%leafp(p)      - pf%m_leafp_to_litter_fire_patch(p)           * dt
         veg_ps%frootp(p)             =  veg_ps%frootp(p)     - pf%m_frootp_to_litter_fire_patch(p)          * dt
         veg_ps%livestemp(p)          =  veg_ps%livestemp(p)  - pf%m_livestemp_to_litter_fire_patch(p)       * dt
         veg_ps%deadstemp(p)          =  veg_ps%deadstemp(p)  - pf%m_deadstemp_to_litter_fire_patch(p)       * dt
         veg_ps%livecrootp(p)         =  veg_ps%livecrootp(p) - pf%m_livecrootp_to_litter_fire_patch(p)      * dt
         veg_ps%deadcrootp(p)         =  veg_ps%deadcrootp(p) - pf%m_deadcrootp_to_litter_fire_patch(p)      * dt

         ! storage pools
         veg_ps%leafp_storage(p)      =  veg_ps%leafp_storage(p)      - pf%m_leafp_storage_to_fire_patch(p)      * dt
         veg_ps%frootp_storage(p)     =  veg_ps%frootp_storage(p)     - pf%m_frootp_storage_to_fire_patch(p)     * dt
         veg_ps%livestemp_storage(p)  =  veg_ps%livestemp_storage(p)  - pf%m_livestemp_storage_to_fire_patch(p)  * dt
         veg_ps%deadstemp_storage(p)  =  veg_ps%deadstemp_storage(p)  - pf%m_deadstemp_storage_to_fire_patch(p)  * dt
         veg_ps%livecrootp_storage(p) =  veg_ps%livecrootp_storage(p) - pf%m_livecrootp_storage_to_fire_patch(p) * dt
         veg_ps%deadcrootp_storage(p) =  veg_ps%deadcrootp_storage(p) - pf%m_deadcrootp_storage_to_fire_patch(p) * dt

         veg_ps%leafp_storage(p)      =  veg_ps%leafp_storage(p)      - pf%m_leafp_storage_to_litter_fire_patch(p)      * dt
         veg_ps%frootp_storage(p)     =  veg_ps%frootp_storage(p)     - pf%m_frootp_storage_to_litter_fire_patch(p)     * dt
         veg_ps%livestemp_storage(p)  =  veg_ps%livestemp_storage(p)  - pf%m_livestemp_storage_to_litter_fire_patch(p)  * dt
         veg_ps%deadstemp_storage(p)  =  veg_ps%deadstemp_storage(p)  - pf%m_deadstemp_storage_to_litter_fire_patch(p)  * dt
         veg_ps%livecrootp_storage(p) =  veg_ps%livecrootp_storage(p) - pf%m_livecrootp_storage_to_litter_fire_patch(p) * dt
         veg_ps%deadcrootp_storage(p) =  veg_ps%deadcrootp_storage(p) - pf%m_deadcrootp_storage_to_litter_fire_patch(p) * dt


         ! trapsfer pools
         veg_ps%leafp_xfer(p)         =  veg_ps%leafp_xfer(p)      - pf%m_leafp_xfer_to_fire_patch(p)      * dt
         veg_ps%frootp_xfer(p)        =  veg_ps%frootp_xfer(p)     - pf%m_frootp_xfer_to_fire_patch(p)     * dt
         veg_ps%livestemp_xfer(p)     =  veg_ps%livestemp_xfer(p)  - pf%m_livestemp_xfer_to_fire_patch(p)  * dt
         veg_ps%deadstemp_xfer(p)     =  veg_ps%deadstemp_xfer(p)  - pf%m_deadstemp_xfer_to_fire_patch(p)  * dt
         veg_ps%livecrootp_xfer(p)    =  veg_ps%livecrootp_xfer(p) - pf%m_livecrootp_xfer_to_fire_patch(p) * dt
         veg_ps%deadcrootp_xfer(p)    =  veg_ps%deadcrootp_xfer(p) - pf%m_deadcrootp_xfer_to_fire_patch(p) * dt

         veg_ps%leafp_xfer(p)         =  veg_ps%leafp_xfer(p)      - pf%m_leafp_xfer_to_litter_fire_patch(p)      * dt
         veg_ps%frootp_xfer(p)        =  veg_ps%frootp_xfer(p)     - pf%m_frootp_xfer_to_litter_fire_patch(p)     * dt
         veg_ps%livestemp_xfer(p)     =  veg_ps%livestemp_xfer(p)  - pf%m_livestemp_xfer_to_litter_fire_patch(p)  * dt
         veg_ps%deadstemp_xfer(p)     =  veg_ps%deadstemp_xfer(p)  - pf%m_deadstemp_xfer_to_litter_fire_patch(p)  * dt
         veg_ps%livecrootp_xfer(p)    =  veg_ps%livecrootp_xfer(p) - pf%m_livecrootp_xfer_to_litter_fire_patch(p) * dt
         veg_ps%deadcrootp_xfer(p)    =  veg_ps%deadcrootp_xfer(p) - pf%m_deadcrootp_xfer_to_litter_fire_patch(p) * dt

         ! retranslocated N pool
         veg_ps%retransp(p)           =  veg_ps%retransp(p) - pf%m_retransp_to_fire_patch(p)        * dt
         veg_ps%retransp(p)           =  veg_ps%retransp(p) - pf%m_retransp_to_litter_fire_patch(p) * dt
         veg_ps%ppool(p)              =  veg_ps%ppool(p) - pf%m_ppool_to_fire_patch(p)              * dt
         veg_ps%ppool(p)              =  veg_ps%ppool(p) - pf%m_ppool_to_litter_fire_patch(p)       * dt
      end do

    end associate 

  end subroutine PhosphorusStateUpdate3

end module PhosphorusStateUpdate3Mod
