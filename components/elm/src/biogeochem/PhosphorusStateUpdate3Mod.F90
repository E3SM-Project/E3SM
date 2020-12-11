module PhosphorusStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for phosphorus state variable update, mortality fluxes.
  ! Also, sminn leaching flux.
  ! X.YANG
  ! !USES:
  use shr_kind_mod        , only: r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use elm_varpar          , only: nlevdecomp,ndecomp_pools,ndecomp_cascade_transitions
  use clm_time_manager    , only : get_step_size
  use elm_varctl          , only : iulog, use_nitrif_denitrif
  use elm_varpar          , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use elm_varctl          , only : use_erosion, ero_ccycle
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType         , only : cnstate_type
  use PhosphorusStateType , only : phosphorusstate_type
  use PhosphorusFLuxType  , only : phosphorusflux_type
  use soilorder_varcon    , only : smax,ks_sorption
  use tracer_varcon       , only : is_active_betr_bgc
  ! bgc interface & pflotran:
  use elm_varctl          , only : use_pflotran, pf_cmode
  use elm_varctl          , only : nu_com
  use elm_varctl          , only : ECA_Pconst_RGspin
  use VegetationPropertiesType      , only : veg_vp 
  use ColumnDataType      , only : col_ps, col_pf
  use VegetationDataType  , only : veg_ps, veg_pf
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
    ! Also the erosion flux.
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
            col_ps%primp_vr(c,j)   = col_ps%primp_vr(c,j) - col_pf%primp_to_labilep_vr(c,j) *dt &
                 + col_pf%pdep_to_sminp(c)*dt * pdep_prof(c,j)
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
                                               col_pf%decomp_cascade_sminp_flux_vr(c,j,k)*dt
               end do
             end do
           else
             do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                    flux_mineralization(c,j) = flux_mineralization(c,j) + &
                                               col_pf%decomp_cascade_sminp_flux_vr(c,j,k)*dt

               end do
             end do
           endif
        end do
   

        do j = 1, nlevdecomp
              ! column loop
           do fc = 1,num_soilc
             c = filter_soilc(fc)
             flux_mineralization(c,j) = flux_mineralization(c,j) + &
                                       col_pf%biochem_pmin_vr(c,j)*dt
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
                    + col_pf%primp_to_labilep_vr(c,j)*dt &
                    + col_pf%secondp_to_labilep_vr(c,j)*dt &
                    + col_pf%supplement_to_sminp_vr(c,j)*dt - col_pf%sminp_to_plant_vr(c,j)*dt&
                    - col_pf%labilep_to_secondp_vr(c,j)*dt - col_pf%sminp_leached_vr(c,j)*dt ) / &
                    (1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+col_ps%solutionp_vr(c,j))**2._r8)

               col_ps%labilep_vr(c,j) = col_ps%labilep_vr(c,j) + ((smax_c*ks_sorption_c)&
                    /(ks_sorption_c+temp_solutionp(c,j))**2._r8 ) * &
                    ( flux_mineralization(c,j) + col_pf%primp_to_labilep_vr(c,j)*dt + col_pf%secondp_to_labilep_vr(c,j)*dt &
                    + col_pf%supplement_to_sminp_vr(c,j)*dt - col_pf%sminp_to_plant_vr(c,j)*dt &
                    - col_pf%labilep_to_secondp_vr(c,j)*dt - col_pf%sminp_leached_vr(c,j)*dt ) / &
                    ( 1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+temp_solutionp(c,j))**2._r8 )
                             
               col_pf%desorb_to_solutionp_vr(c,j) = ( flux_mineralization(c,j)/dt + col_pf%primp_to_labilep_vr(c,j) &
                                + col_pf%secondp_to_labilep_vr(c,j) &
                                + col_pf%supplement_to_sminp_vr(c,j) - col_pf%sminp_to_plant_vr(c,j) &
                                - col_pf%labilep_to_secondp_vr(c,j) - col_pf%sminp_leached_vr(c,j) ) / &
                                (1._r8+(smax_c*ks_sorption_c)/(ks_sorption_c+col_ps%solutionp_vr(c,j))**2._r8)
            
               col_pf%adsorb_to_labilep_vr(c,j) = ((smax_c*ks_sorption_c)/(ks_sorption_c+temp_solutionp(c,j))**2._r8 ) * &
                             ( flux_mineralization(c,j)/dt + col_pf%primp_to_labilep_vr(c,j) + col_pf%secondp_to_labilep_vr(c,j) &
                             + col_pf%supplement_to_sminp_vr(c,j) - col_pf%sminp_to_plant_vr(c,j) &
                             - col_pf%labilep_to_secondp_vr(c,j) - col_pf%sminp_leached_vr(c,j) ) / &
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
                            (flux_mineralization(c,j) + col_pf%primp_to_labilep_vr(c,j)*dt + &
                            col_pf%secondp_to_labilep_vr(c,j)*dt + col_pf%supplement_to_sminp_vr(c,j)*dt - &
                            col_pf%sminp_to_plant_vr(c,j)*dt - col_pf%labilep_to_secondp_vr(c,j)*dt - &
                            col_pf%sminp_leached_vr(c,j)*dt ))

                 if (temp_solutionp(c,j) < 0.0_r8) then
                    col_pf%labilep_to_secondp_vr(c,j) = col_pf%labilep_to_secondp_vr(c,j)/ &
                            (col_pf%labilep_to_secondp_vr(c,j)+col_pf%sminp_leached_vr(c,j))* &
                            (temp_solutionp(c,j) + col_pf%labilep_to_secondp_vr(c,j)*dt + &
                            col_pf%sminp_leached_vr(c,j)*dt) /dt
                    col_pf%sminp_leached_vr(c,j) = col_pf%sminp_leached_vr(c,j)/ &
                            (col_pf%labilep_to_secondp_vr(c,j)+col_pf%sminp_leached_vr(c,j))* &
                            (temp_solutionp(c,j) + col_pf%labilep_to_secondp_vr(c,j)*dt + &
                            col_pf%sminp_leached_vr(c,j)*dt) /dt
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

                   col_ps%decomp_ppools_vr(c,j,l) = col_ps%decomp_ppools_vr(c,j,l)- col_pf%biochem_pmin_ppools_vr(c,j,l)*dt

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

               col_ps%secondp_vr(c,j) = col_ps%secondp_vr(c,j) + ( col_pf%labilep_to_secondp_vr(c,j) &
                    - col_pf%secondp_to_labilep_vr(c,j) &
                                     - col_pf%secondp_to_occlp_vr(c,j) )*dt
               col_ps%occlp_vr(c,j)   = col_ps%occlp_vr(c,j) + ( col_pf%secondp_to_occlp_vr(c,j) ) * dt
               col_ps%primp_vr(c,j)   = col_ps%primp_vr(c,j) - ( col_pf%primp_to_labilep_vr(c,j) )*dt + col_pf%pdep_to_sminp(c)*dt &
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
            col_ps%decomp_ppools_vr(c,j,i_cwd) = col_ps%decomp_ppools_vr(c,j,i_cwd) + col_pf%fire_mortality_p_to_cwdp(c,j) * dt

            ! pft-level wood to column-level litter (uncombusted wood)
            col_ps%decomp_ppools_vr(c,j,i_met_lit) = col_ps%decomp_ppools_vr(c,j,i_met_lit) + col_pf%m_p_to_litr_met_fire(c,j)* dt
            col_ps%decomp_ppools_vr(c,j,i_cel_lit) = col_ps%decomp_ppools_vr(c,j,i_cel_lit) + col_pf%m_p_to_litr_cel_fire(c,j)* dt
            col_ps%decomp_ppools_vr(c,j,i_lig_lit) = col_ps%decomp_ppools_vr(c,j,i_lig_lit) + col_pf%m_p_to_litr_lig_fire(c,j)* dt
         end do ! end of column loop
      end do

      ! litter and CWD losses to fire
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               col_ps%decomp_ppools_vr(c,j,l) = col_ps%decomp_ppools_vr(c,j,l) - col_pf%m_decomp_ppools_to_fire_vr(c,j,l) * dt
            end do
         end do
      end do

    endif !is_active_betr_bgc

    ! soil P loss due to soil erosion
    if ( ero_ccycle ) then
      ! column loop
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         do j = 1, nlevdecomp
            ! pool loop
            do l = 1, ndecomp_pools
               if ( decomp_cascade_con%is_soil(l) ) then
                  col_ps%decomp_ppools_vr(c,j,l) = col_ps%decomp_ppools_vr(c,j,l) - col_pf%decomp_ppools_yield_vr(c,j,l) * dt
               end if
            end do
            col_ps%labilep_vr(c,j) = col_ps%labilep_vr(c,j) - col_pf%labilep_yield_vr(c,j) * dt
            col_ps%secondp_vr(c,j) = col_ps%secondp_vr(c,j) - col_pf%secondp_yield_vr(c,j) * dt
            col_ps%occlp_vr(c,j)   = col_ps%occlp_vr(c,j)   - col_pf%occlp_yield_vr(c,j) * dt
            col_ps%primp_vr(c,j)   = col_ps%primp_vr(c,j)   - col_pf%primp_yield_vr(c,j) * dt
         end do
      end do
    end if

      ! patch-level phosphorus fluxes 

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         !from fire displayed pools
         veg_ps%leafp(p)              =  veg_ps%leafp(p)      - veg_pf%m_leafp_to_fire(p)      * dt
         veg_ps%frootp(p)             =  veg_ps%frootp(p)     - veg_pf%m_frootp_to_fire(p)     * dt
         veg_ps%livestemp(p)          =  veg_ps%livestemp(p)  - veg_pf%m_livestemp_to_fire(p)  * dt
         veg_ps%deadstemp(p)          =  veg_ps%deadstemp(p)  - veg_pf%m_deadstemp_to_fire(p)  * dt
         veg_ps%livecrootp(p)         =  veg_ps%livecrootp(p) - veg_pf%m_livecrootp_to_fire(p) * dt
         veg_ps%deadcrootp(p)         =  veg_ps%deadcrootp(p) - veg_pf%m_deadcrootp_to_fire(p) * dt

         veg_ps%leafp(p)              =  veg_ps%leafp(p)      - veg_pf%m_leafp_to_litter_fire(p)           * dt
         veg_ps%frootp(p)             =  veg_ps%frootp(p)     - veg_pf%m_frootp_to_litter_fire(p)          * dt
         veg_ps%livestemp(p)          =  veg_ps%livestemp(p)  - veg_pf%m_livestemp_to_litter_fire(p)       * dt
         veg_ps%deadstemp(p)          =  veg_ps%deadstemp(p)  - veg_pf%m_deadstemp_to_litter_fire(p)       * dt
         veg_ps%livecrootp(p)         =  veg_ps%livecrootp(p) - veg_pf%m_livecrootp_to_litter_fire(p)      * dt
         veg_ps%deadcrootp(p)         =  veg_ps%deadcrootp(p) - veg_pf%m_deadcrootp_to_litter_fire(p)      * dt

         ! storage pools
         veg_ps%leafp_storage(p)      =  veg_ps%leafp_storage(p)      - veg_pf%m_leafp_storage_to_fire(p)      * dt
         veg_ps%frootp_storage(p)     =  veg_ps%frootp_storage(p)     - veg_pf%m_frootp_storage_to_fire(p)     * dt
         veg_ps%livestemp_storage(p)  =  veg_ps%livestemp_storage(p)  - veg_pf%m_livestemp_storage_to_fire(p)  * dt
         veg_ps%deadstemp_storage(p)  =  veg_ps%deadstemp_storage(p)  - veg_pf%m_deadstemp_storage_to_fire(p)  * dt
         veg_ps%livecrootp_storage(p) =  veg_ps%livecrootp_storage(p) - veg_pf%m_livecrootp_storage_to_fire(p) * dt
         veg_ps%deadcrootp_storage(p) =  veg_ps%deadcrootp_storage(p) - veg_pf%m_deadcrootp_storage_to_fire(p) * dt

         veg_ps%leafp_storage(p)      =  veg_ps%leafp_storage(p)      - veg_pf%m_leafp_storage_to_litter_fire(p)      * dt
         veg_ps%frootp_storage(p)     =  veg_ps%frootp_storage(p)     - veg_pf%m_frootp_storage_to_litter_fire(p)     * dt
         veg_ps%livestemp_storage(p)  =  veg_ps%livestemp_storage(p)  - veg_pf%m_livestemp_storage_to_litter_fire(p)  * dt
         veg_ps%deadstemp_storage(p)  =  veg_ps%deadstemp_storage(p)  - veg_pf%m_deadstemp_storage_to_litter_fire(p)  * dt
         veg_ps%livecrootp_storage(p) =  veg_ps%livecrootp_storage(p) - veg_pf%m_livecrootp_storage_to_litter_fire(p) * dt
         veg_ps%deadcrootp_storage(p) =  veg_ps%deadcrootp_storage(p) - veg_pf%m_deadcrootp_storage_to_litter_fire(p) * dt


         ! trapsfer pools
         veg_ps%leafp_xfer(p)         =  veg_ps%leafp_xfer(p)      - veg_pf%m_leafp_xfer_to_fire(p)      * dt
         veg_ps%frootp_xfer(p)        =  veg_ps%frootp_xfer(p)     - veg_pf%m_frootp_xfer_to_fire(p)     * dt
         veg_ps%livestemp_xfer(p)     =  veg_ps%livestemp_xfer(p)  - veg_pf%m_livestemp_xfer_to_fire(p)  * dt
         veg_ps%deadstemp_xfer(p)     =  veg_ps%deadstemp_xfer(p)  - veg_pf%m_deadstemp_xfer_to_fire(p)  * dt
         veg_ps%livecrootp_xfer(p)    =  veg_ps%livecrootp_xfer(p) - veg_pf%m_livecrootp_xfer_to_fire(p) * dt
         veg_ps%deadcrootp_xfer(p)    =  veg_ps%deadcrootp_xfer(p) - veg_pf%m_deadcrootp_xfer_to_fire(p) * dt

         veg_ps%leafp_xfer(p)         =  veg_ps%leafp_xfer(p)      - veg_pf%m_leafp_xfer_to_litter_fire(p)      * dt
         veg_ps%frootp_xfer(p)        =  veg_ps%frootp_xfer(p)     - veg_pf%m_frootp_xfer_to_litter_fire(p)     * dt
         veg_ps%livestemp_xfer(p)     =  veg_ps%livestemp_xfer(p)  - veg_pf%m_livestemp_xfer_to_litter_fire(p)  * dt
         veg_ps%deadstemp_xfer(p)     =  veg_ps%deadstemp_xfer(p)  - veg_pf%m_deadstemp_xfer_to_litter_fire(p)  * dt
         veg_ps%livecrootp_xfer(p)    =  veg_ps%livecrootp_xfer(p) - veg_pf%m_livecrootp_xfer_to_litter_fire(p) * dt
         veg_ps%deadcrootp_xfer(p)    =  veg_ps%deadcrootp_xfer(p) - veg_pf%m_deadcrootp_xfer_to_litter_fire(p) * dt

         ! retranslocated N pool
         veg_ps%retransp(p)           =  veg_ps%retransp(p) - veg_pf%m_retransp_to_fire(p)        * dt
         veg_ps%retransp(p)           =  veg_ps%retransp(p) - veg_pf%m_retransp_to_litter_fire(p) * dt
         veg_ps%ppool(p)              =  veg_ps%ppool(p) - veg_pf%m_ppool_to_fire(p)              * dt
         veg_ps%ppool(p)              =  veg_ps%ppool(p) - veg_pf%m_ppool_to_litter_fire(p)       * dt
      end do

    end associate 

  end subroutine PhosphorusStateUpdate3

end module PhosphorusStateUpdate3Mod
