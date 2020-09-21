module NitrogenStateUpdate1Mod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod           , only: r8 => shr_kind_r8
  use clm_time_manager       , only : get_step_size
  use elm_varpar             , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use elm_varpar             , only : crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl             , only : iulog, use_nitrif_denitrif
  use elm_varcon             , only : nitrif_n2o_loss_frac
  use pftvarcon              , only : npcropmin, nc3crop
  use VegetationPropertiesType         , only : veg_vp
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use GridcellDataType       , only : grc_ns, grc_nf
  use ColumnDataType         , only : col_ns, col_nf
  use VegetationType         , only : veg_pp
  use VegetationDataType     , only : veg_ns, veg_nf
  use tracer_varcon          , only : is_active_betr_bgc
  ! bgc interface & pflotran:
  use clm_varctl             , only : use_pflotran, pf_cmode
  ! forest fertilization experiment
  use clm_time_manager       , only : get_curr_date
  use CNStateType            , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use clm_varctl             , only : forest_fert_exp
  use clm_varctl             , only : nu_com
  use clm_varctl             , only : NFIX_PTASE_plant
  use decompMod              , only : bounds_type
  use elm_varcon             , only : dzsoi_decomp
  use clm_varctl             , only : use_fates
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: NitrogenStateUpdateDynPatch
  public :: NitrogenStateUpdate1
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdateDynPatch(bounds, num_soilc_with_inactive, &
       filter_soilc_with_inactive)
    !
    ! !DESCRIPTION:
    ! Update nitrogen states based on fluxes from dyn_cnbal_patch
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc_with_inactive       ! number of columns in soil filter
    integer                  , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    !
    ! !LOCAL VARIABLES:
    integer                                  :: c                             ! column index
    integer                                  :: fc                            ! column filter index
    integer                                  :: g                             ! gridcell index
    integer                                  :: j                             ! level index
    real(r8)                                 :: dt                            ! time step (seconds)

    character(len=*)         , parameter     :: subname = 'NitrogenStateUpdateDynPatch'
    !-----------------------------------------------------------------------

      dt = real( get_step_size(), r8 )

      if (.not.use_fates) then

         do g = bounds%begg, bounds%endg
            grc_ns%seedn(g) = grc_ns%seedn(g) &
                 - grc_nf%dwt_seedn_to_leaf(g)     * dt &
                 - grc_nf%dwt_seedn_to_deadstem(g) * dt &
                 - grc_nf%dwt_seedn_to_npool(g)    * dt
         end do

         do fc = 1, num_soilc_with_inactive
            
            c = filter_soilc_with_inactive(fc)
            col_ns%prod10n(c) = col_ns%prod10n(c) + col_nf%dwt_prod10n_gain(c)*dt
            col_ns%prod100n(c) = col_ns%prod100n(c) + col_nf%dwt_prod100n_gain(c)*dt
            col_ns%prod1n(c) = col_ns%prod1n(c) + col_nf%dwt_crop_productn_gain(c)*dt

            do j = 1,nlevdecomp

               col_ns%decomp_npools_vr(c,j,i_met_lit) = col_ns%decomp_npools_vr(c,j,i_met_lit) + &
                    col_nf%dwt_frootn_to_litr_met_n(c,j) * dt
               col_ns%decomp_npools_vr(c,j,i_cel_lit) = col_ns%decomp_npools_vr(c,j,i_cel_lit) + &
                    col_nf%dwt_frootn_to_litr_cel_n(c,j) * dt
               col_ns%decomp_npools_vr(c,j,i_lig_lit) = col_ns%decomp_npools_vr(c,j,i_lig_lit) + &
                    col_nf%dwt_frootn_to_litr_lig_n(c,j) * dt
               col_ns%decomp_npools_vr(c,j,i_cwd) = col_ns%decomp_npools_vr(c,j,i_cwd) + &
                    ( col_nf%dwt_livecrootn_to_cwdn(c,j) + col_nf%dwt_deadcrootn_to_cwdn(c,j) ) * dt

            end do
         end do
      end if

  end subroutine NitrogenStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, nitrogenflux_vars, nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    use tracer_varcon, only : is_active_betr_bgc      
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    real(r8), parameter :: frootc_nfix_thc = 10._r8  !threshold fine root carbon for nitrogen fixation gC/m2

    integer:: kyr                     ! current year 
    integer:: kmo                     ! month of year  (1, ..., 12)
    integer:: kda                     ! day of month   (1, ..., 31) 
    integer:: mcsec                   ! seconds of day (0, ..., seconds/day) 
    !-----------------------------------------------------------------------

    associate(                                                                                           & 
         ivt                   => veg_pp%itype                                , & ! Input:  [integer  (:)     ]  pft vegetation type                                

         woody                 => veg_vp%woody                         , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)

         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool    , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         ndep_prof             => cnstate_vars%ndep_prof_col               , & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
         nfixation_prof        => cnstate_vars%nfixation_prof_col          , & ! Input:  [real(r8) (:,:)   ]  profile over which N fixation is distributed through column (1/m)
         
         nf                    => nitrogenflux_vars                        , &
         ns                    => nitrogenstate_vars &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column-level fluxes

      if (.not. is_active_betr_bgc .and. .not.(use_pflotran .and. pf_cmode)) then

         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               
               if (.not. use_nitrif_denitrif) then
                  
                  ! N deposition and fixation
                  col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) + col_nf%ndep_to_sminn(c)*dt * ndep_prof(c,j)
                  col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) + col_nf%nfix_to_sminn(c)*dt * nfixation_prof(c,j)
                  
               else

                  ! N deposition and fixation (put all into NH4 pool)
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%ndep_to_sminn(c)*dt * ndep_prof(c,j)
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%nfix_to_sminn(c)*dt * nfixation_prof(c,j)
                  
               end if

               ! plant to litter fluxes
               ! phenology and dynamic landcover fluxes
               col_nf%decomp_npools_sourcesink(c,j,i_met_lit) = &
                    col_nf%phenology_n_to_litr_met_n(c,j) * dt
               
               col_nf%decomp_npools_sourcesink(c,j,i_cel_lit) = &
                    col_nf%phenology_n_to_litr_cel_n(c,j) * dt
               
               col_nf%decomp_npools_sourcesink(c,j,i_lig_lit) = &
                    col_nf%phenology_n_to_litr_lig_n(c,j) * dt
            end do
         end do
         
         ! repeating N dep and fixation for crops
         if ( crop_prog )then
            do j = 1, nlevdecomp
               
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (.not. use_nitrif_denitrif) then
                     
                     ! N deposition and fixation
                     col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) + col_nf%fert_to_sminn(c)*dt * ndep_prof(c,j)
                     col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) + col_nf%soyfixn_to_sminn(c)*dt * nfixation_prof(c,j)
                  else
                     
                     ! N deposition and fixation (put all into NH4 pool)
                     col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%fert_to_sminn(c)*dt * ndep_prof(c,j)
                     col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%soyfixn_to_sminn(c)*dt * nfixation_prof(c,j)
                     
                  end if
               end do
            end do
         end if
         
         ! decomposition fluxes
         do k = 1, ndecomp_cascade_transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  col_nf%decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) = &
                       col_nf%decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) - &
                       col_nf%decomp_cascade_ntransfer_vr(c,j,k) * dt
               end do
            end do
         end do

         do k = 1, ndecomp_cascade_transitions
            if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     
                     col_nf%decomp_npools_sourcesink(c,j,cascade_receiver_pool(k)) = &
                          col_nf%decomp_npools_sourcesink(c,j,cascade_receiver_pool(k)) + &
                          (col_nf%decomp_cascade_ntransfer_vr(c,j,k) + col_nf%decomp_cascade_sminn_flux_vr(c,j,k)) * dt
                  end do
               end do
            else  ! terminal transitions
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     col_nf%decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) = &
                          col_nf%decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) - &
                          col_nf%decomp_cascade_sminn_flux_vr(c,j,k) * dt
                  end do
               end do
            end if
         end do
         
         if (.not. use_nitrif_denitrif) then
            
            !--------------------------------------------------------
            !-------------    NITRIF_DENITRIF OFF -------------------
            !--------------------------------------------------------
            
            ! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes and denitrification fluxes
            do k = 1, ndecomp_cascade_transitions
               if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
                  do j = 1, nlevdecomp
                     ! column loop
                     do fc = 1,num_soilc
                        c = filter_soilc(fc)
                        col_ns%sminn_vr(c,j)  = col_ns%sminn_vr(c,j) - &
                             (col_nf%sminn_to_denit_decomp_cascade_vr(c,j,k) + col_nf%decomp_cascade_sminn_flux_vr(c,j,k))* dt
                     end do
                  end do
               else
                  do j = 1, nlevdecomp
                     ! column loop
                     do fc = 1,num_soilc
                        c = filter_soilc(fc)
                        col_ns%sminn_vr(c,j)  = col_ns%sminn_vr(c,j) - col_nf%sminn_to_denit_decomp_cascade_vr(c,j,k)* dt
                        
                        col_ns%sminn_vr(c,j)  = col_ns%sminn_vr(c,j) + col_nf%decomp_cascade_sminn_flux_vr(c,j,k)* dt
                        
                     end do
                  end do
               endif
            end do
            
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  ! "bulk denitrification"
                  col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) - col_nf%sminn_to_denit_excess_vr(c,j) * dt
                  
                  ! total plant uptake from mineral N
                  col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) - col_nf%sminn_to_plant_vr(c,j)*dt
                  
                  ! flux that prevents N limitation (when Carbon_only is set)
                  col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) + col_nf%supplement_to_sminn_vr(c,j)*dt
               end do
            end do
            
         else   
            
            !--------------------------------------------------------
            !-------------    NITRIF_DENITRIF ON --------------------
            !--------------------------------------------------------
            
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  
                  ! mineralization fluxes (divert a fraction of this stream to nitrification flux, add the rest to NH4 pool)
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%gross_nmin_vr(c,j)*dt
                  
                  ! immobilization fluxes
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) - col_nf%actual_immob_nh4_vr(c,j)*dt
                  
                  col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) - col_nf%actual_immob_no3_vr(c,j)*dt
                  
                  ! plant uptake fluxes
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) - col_nf%smin_nh4_to_plant_vr(c,j)*dt
                  
                  col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) - col_nf%smin_no3_to_plant_vr(c,j)*dt
                  
                  ! Account for nitrification fluxes
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) - col_nf%f_nit_vr(c,j) * dt
                  
                  col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) + col_nf%f_nit_vr(c,j) * dt * (1._r8 - nitrif_n2o_loss_frac)
                  
                  ! Account for denitrification fluxes
                  col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) - col_nf%f_denit_vr(c,j) * dt
                  
                  ! flux that prevents N limitation (when Carbon_only is set; put all into NH4)
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%supplement_to_sminn_vr(c,j)*dt
                  
                  ! update diagnostic total
                  col_ns%sminn_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_ns%smin_no3_vr(c,j)
                  
               end do ! end of column loop
            end do
            
         end if
      endif  !end if is_active_betr_bgc 

      ! forest fertilization
      call get_curr_date(kyr, kmo, kda, mcsec)
      if (forest_fert_exp) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 1 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               if (.not. use_nitrif_denitrif) then
                  do j = 1, nlevdecomp
                     col_ns%sminn_vr(c,j) = col_ns%sminn_vr(c,j) + fert_dose(c,kmo)*ndep_prof(c,j)
                  end do
               else
                  do j = 1, nlevdecomp
                     col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + fert_dose(c,kmo)/2._r8*ndep_prof(c,j)
                     col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) + fert_dose(c,kmo)/2._r8*ndep_prof(c,j)
                     col_ns%sminn_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_ns%smin_no3_vr(c,j)
                  end do
               end if
            end if
         end do
      end if

      ! patch loop
      
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! phenology: transfer growth fluxes
         veg_ns%leafn(p)       = veg_ns%leafn(p)       + veg_nf%leafn_xfer_to_leafn(p)*dt
         veg_ns%leafn_xfer(p)  = veg_ns%leafn_xfer(p)  - veg_nf%leafn_xfer_to_leafn(p)*dt
         veg_ns%frootn(p)      = veg_ns%frootn(p)      + veg_nf%frootn_xfer_to_frootn(p)*dt
         veg_ns%frootn_xfer(p) = veg_ns%frootn_xfer(p) - veg_nf%frootn_xfer_to_frootn(p)*dt

         if (woody(ivt(p)) == 1.0_r8) then
            veg_ns%livestemn(p)       = veg_ns%livestemn(p)       + veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%livestemn_xfer(p)  = veg_ns%livestemn_xfer(p)  - veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%deadstemn(p)       = veg_ns%deadstemn(p)       + veg_nf%deadstemn_xfer_to_deadstemn(p)*dt
            veg_ns%deadstemn_xfer(p)  = veg_ns%deadstemn_xfer(p)  - veg_nf%deadstemn_xfer_to_deadstemn(p)*dt
            veg_ns%livecrootn(p)      = veg_ns%livecrootn(p)      + veg_nf%livecrootn_xfer_to_livecrootn(p)*dt
            veg_ns%livecrootn_xfer(p) = veg_ns%livecrootn_xfer(p) - veg_nf%livecrootn_xfer_to_livecrootn(p)*dt
            veg_ns%deadcrootn(p)      = veg_ns%deadcrootn(p)      + veg_nf%deadcrootn_xfer_to_deadcrootn(p)*dt
            veg_ns%deadcrootn_xfer(p) = veg_ns%deadcrootn_xfer(p) - veg_nf%deadcrootn_xfer_to_deadcrootn(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            veg_ns%livestemn(p)       = veg_ns%livestemn(p)      + veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%livestemn_xfer(p)  = veg_ns%livestemn_xfer(p) - veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%grainn(p)          = veg_ns%grainn(p)         + veg_nf%grainn_xfer_to_grainn(p)*dt
            veg_ns%grainn_xfer(p)     = veg_ns%grainn_xfer(p)    - veg_nf%grainn_xfer_to_grainn(p)*dt
         end if

         ! phenology: litterfall and retranslocation fluxes
         veg_ns%leafn(p)    = veg_ns%leafn(p)    - veg_nf%leafn_to_litter(p)*dt
         veg_ns%frootn(p)   = veg_ns%frootn(p)   - veg_nf%frootn_to_litter(p)*dt
         veg_ns%leafn(p)    = veg_ns%leafn(p)    - veg_nf%leafn_to_retransn(p)*dt
         veg_ns%retransn(p) = veg_ns%retransn(p) + veg_nf%leafn_to_retransn(p)*dt

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt(p)) == 1._r8) then
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_deadstemn(p)*dt
            veg_ns%deadstemn(p)  = veg_ns%deadstemn(p)  + veg_nf%livestemn_to_deadstemn(p)*dt
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%livecrootn(p) = veg_ns%livecrootn(p) - veg_nf%livecrootn_to_deadcrootn(p)*dt
            veg_ns%deadcrootn(p) = veg_ns%deadcrootn(p) + veg_nf%livecrootn_to_deadcrootn(p)*dt
            veg_ns%livecrootn(p) = veg_ns%livecrootn(p) - veg_nf%livecrootn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%livecrootn_to_retransn(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
            veg_ns%frootn(p)     = veg_ns%frootn(p)     - veg_nf%frootn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%frootn_to_retransn(p)*dt
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_litter(p)*dt
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%grainn(p)     = veg_ns%grainn(p)     - veg_nf%grainn_to_food(p)*dt

            veg_ns%cropseedn_deficit(p) = veg_ns%cropseedn_deficit(p) &
                 - veg_nf%crop_seedn_to_leaf(p) * dt
         end if

         ! uptake from soil mineral N pool
         veg_ns%npool(p) = &
              veg_ns%npool(p) + veg_nf%sminn_to_npool(p)*dt
         if (nu_com .ne. 'RD') veg_ns%npool(p) = veg_ns%npool(p) + veg_nf%supplement_to_plantn(p)*dt
         if (NFIX_PTASE_plant) veg_ns%npool(p) = veg_ns%npool(p) + veg_nf%nfix_to_plantn(p)*dt

         ! deployment from retranslocation pool
         veg_ns%npool(p)    = veg_ns%npool(p)    + veg_nf%retransn_to_npool(p)*dt
         veg_ns%retransn(p) = veg_ns%retransn(p) - veg_nf%retransn_to_npool(p)*dt

         ! allocation fluxes
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_leafn(p)*dt
         veg_ns%leafn(p)           = veg_ns%leafn(p)          + veg_nf%npool_to_leafn(p)*dt
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_leafn_storage(p)*dt
         veg_ns%leafn_storage(p)   = veg_ns%leafn_storage(p)  + veg_nf%npool_to_leafn_storage(p)*dt
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_frootn(p)*dt
         veg_ns%frootn(p)          = veg_ns%frootn(p)         + veg_nf%npool_to_frootn(p)*dt
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_frootn_storage(p)*dt
         veg_ns%frootn_storage(p)  = veg_ns%frootn_storage(p) + veg_nf%npool_to_frootn_storage(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn(p)*dt
            veg_ns%livestemn(p)          = veg_ns%livestemn(p)          + veg_nf%npool_to_livestemn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p)  + veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadstemn(p)*dt
            veg_ns%deadstemn(p)          = veg_ns%deadstemn(p)          + veg_nf%npool_to_deadstemn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadstemn_storage(p)*dt
            veg_ns%deadstemn_storage(p)  = veg_ns%deadstemn_storage(p)  + veg_nf%npool_to_deadstemn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livecrootn(p)*dt
            veg_ns%livecrootn(p)         = veg_ns%livecrootn(p)         + veg_nf%npool_to_livecrootn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livecrootn_storage(p)*dt
            veg_ns%livecrootn_storage(p) = veg_ns%livecrootn_storage(p) + veg_nf%npool_to_livecrootn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadcrootn(p)*dt
            veg_ns%deadcrootn(p)         = veg_ns%deadcrootn(p)         + veg_nf%npool_to_deadcrootn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadcrootn_storage(p)*dt
            veg_ns%deadcrootn_storage(p) = veg_ns%deadcrootn_storage(p) + veg_nf%npool_to_deadcrootn_storage(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn(p)*dt
            veg_ns%livestemn(p)          = veg_ns%livestemn(p)          + veg_nf%npool_to_livestemn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p)  + veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_grainn(p)*dt
            veg_ns%grainn(p)             = veg_ns%grainn(p)             + veg_nf%npool_to_grainn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_grainn_storage(p)*dt
            veg_ns%grainn_storage(p)     = veg_ns%grainn_storage(p)     + veg_nf%npool_to_grainn_storage(p)*dt
         end if

         ! move storage pools into transfer pools
         veg_ns%leafn_storage(p)  = veg_ns%leafn_storage(p)  - veg_nf%leafn_storage_to_xfer(p)*dt
         veg_ns%leafn_xfer(p)     = veg_ns%leafn_xfer(p)     + veg_nf%leafn_storage_to_xfer(p)*dt
         veg_ns%frootn_storage(p) = veg_ns%frootn_storage(p) - veg_nf%frootn_storage_to_xfer(p)*dt
         veg_ns%frootn_xfer(p)    = veg_ns%frootn_xfer(p)    + veg_nf%frootn_storage_to_xfer(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p)  - veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%livestemn_xfer(p)     = veg_ns%livestemn_xfer(p)     + veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%deadstemn_storage(p)  = veg_ns%deadstemn_storage(p)  - veg_nf%deadstemn_storage_to_xfer(p)*dt
            veg_ns%deadstemn_xfer(p)     = veg_ns%deadstemn_xfer(p)     + veg_nf%deadstemn_storage_to_xfer(p)*dt
            veg_ns%livecrootn_storage(p) = veg_ns%livecrootn_storage(p) - veg_nf%livecrootn_storage_to_xfer(p)*dt
            veg_ns%livecrootn_xfer(p)    = veg_ns%livecrootn_xfer(p)    + veg_nf%livecrootn_storage_to_xfer(p)*dt
            veg_ns%deadcrootn_storage(p) = veg_ns%deadcrootn_storage(p) - veg_nf%deadcrootn_storage_to_xfer(p)*dt
            veg_ns%deadcrootn_xfer(p)    = veg_ns%deadcrootn_xfer(p)    + veg_nf%deadcrootn_storage_to_xfer(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p) - veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%livestemn_xfer(p)     = veg_ns%livestemn_xfer(p)    + veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%grainn_storage(p)     = veg_ns%grainn_storage(p)    - veg_nf%grainn_storage_to_xfer(p)*dt
            veg_ns%grainn_xfer(p)        = veg_ns%grainn_xfer(p)       + veg_nf%grainn_storage_to_xfer(p)*dt
         end if

      end do

    end associate

  end subroutine NitrogenStateUpdate1

end module NitrogenStateUpdate1Mod
