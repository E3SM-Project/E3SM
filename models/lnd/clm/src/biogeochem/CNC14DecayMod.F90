module CNC14DecayMod

  !-----------------------------------------------------------------------
  ! Module for 14-carbon flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use clm_time_manager                   , only : get_step_size, get_days_per_year
  use clm_varpar                         , only : ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  use clm_varcon                         , only : secspday
  use clm_varctl                         , only : spinup_state
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: C14Decay
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine C14Decay( num_soilc, filter_soilc, num_soilp, filter_soilp, &
       c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the radioactive decay of C14
    !
    ! !ARGUMENTS:
    integer                               , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(CNVeg_carbonstate_type)          , intent(inout) :: c14_cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type) , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,j,l,p,fc,c,i
    real(r8) :: dt            ! radiation time step (seconds)
    real(r8) :: half_life
    real(r8) :: decay_const
    real(r8) :: days_per_year ! days per year
    real(r8) :: spinup_term   ! spinup accelerated decomposition factor, used to accelerate transport as well
    !-----------------------------------------------------------------------

    associate(                                                                               & 
         spinup_factor      =>    decomp_cascade_con%spinup_factor                         , & ! Input:   [real(r8) (:)     ]  factor for AD spinup associated with each pool

         decomp_cpools_vr   =>    c14_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col , & ! Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

         seedc              =>    c14_cnveg_carbonstate_inst%seedc_col                     , & ! Output:  [real(r8) (:)     ]                                          
         cpool              =>    c14_cnveg_carbonstate_inst%cpool_patch                   , & ! Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool  
         xsmrpool           =>    c14_cnveg_carbonstate_inst%xsmrpool_patch                , & ! Output:  [real(r8) (:)     ]  (gC/m2) execss maint resp C pool        
         deadcrootc         =>    c14_cnveg_carbonstate_inst%deadcrootc_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C              
         deadcrootc_storage =>    c14_cnveg_carbonstate_inst%deadcrootc_storage_patch      , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage      
         deadcrootc_xfer    =>    c14_cnveg_carbonstate_inst%deadcrootc_xfer_patch         , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer     
         deadstemc          =>    c14_cnveg_carbonstate_inst%deadstemc_patch               , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                     
         deadstemc_storage  =>    c14_cnveg_carbonstate_inst%deadstemc_storage_patch       , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage             
         deadstemc_xfer     =>    c14_cnveg_carbonstate_inst%deadstemc_xfer_patch          , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer            
         frootc             =>    c14_cnveg_carbonstate_inst%frootc_patch                  , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C                     
         frootc_storage     =>    c14_cnveg_carbonstate_inst%frootc_storage_patch          , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage             
         frootc_xfer        =>    c14_cnveg_carbonstate_inst%frootc_xfer_patch             , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer            
         gresp_storage      =>    c14_cnveg_carbonstate_inst%gresp_storage_patch           , & ! Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage      
         gresp_xfer         =>    c14_cnveg_carbonstate_inst%gresp_xfer_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer     
         leafc              =>    c14_cnveg_carbonstate_inst%leafc_patch                   , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C                          
         leafc_storage      =>    c14_cnveg_carbonstate_inst%leafc_storage_patch           , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                  
         leafc_xfer         =>    c14_cnveg_carbonstate_inst%leafc_xfer_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                 
         livecrootc         =>    c14_cnveg_carbonstate_inst%livecrootc_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C              
         livecrootc_storage =>    c14_cnveg_carbonstate_inst%livecrootc_storage_patch      , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage      
         livecrootc_xfer    =>    c14_cnveg_carbonstate_inst%livecrootc_xfer_patch         , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer     
         livestemc          =>    c14_cnveg_carbonstate_inst%livestemc_patch               , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C                     
         livestemc_storage  =>    c14_cnveg_carbonstate_inst%livestemc_storage_patch       , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage             
         livestemc_xfer     =>    c14_cnveg_carbonstate_inst%livestemc_xfer_patch          , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer            
         pft_ctrunc         =>    c14_cnveg_carbonstate_inst%ctrunc_patch                    & ! Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation 
         )

      ! set time steps
      dt = real( get_step_size(), r8 )
      days_per_year = get_days_per_year()

      half_life = 5568._r8 * secspday * days_per_year  !! libby half-life value, for comparison against ages calculated with this value
      ! half_life = 5730._r8 * secspday * days_per_year  !! recent half-life value
      decay_const = - log(0.5_r8) / half_life

      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         seedc(c) = seedc(c) *  (1._r8 - decay_const * dt)
      end do ! end of columns loop

      do l = 1, ndecomp_pools
         if ( spinup_state == 1) then
            ! speed up radioactive decay by the same factor as decomposition so tat SOM ages prematurely in all respects
            spinup_term = spinup_factor(l) 
         else
            spinup_term = 1.
         endif
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               decomp_cpools_vr(c,j,l) = decomp_cpools_vr(c,j,l) * (1._r8 - decay_const * spinup_term * dt)
            end do
         end do
      end do ! end of columns loop

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         cpool(p)              = cpool(p)               * (1._r8 - decay_const * dt)
         xsmrpool(p)           = xsmrpool(p)            * (1._r8 - decay_const * dt)
         deadcrootc(p)         = deadcrootc(p)          * (1._r8 - decay_const * dt)
         deadcrootc_storage(p) = deadcrootc_storage(p)  * (1._r8 - decay_const * dt)
         deadcrootc_xfer(p)    = deadcrootc_xfer(p)     * (1._r8 - decay_const * dt)
         deadstemc(p)          = deadstemc(p)           * (1._r8 - decay_const * dt)
         deadstemc_storage(p)  = deadstemc_storage(p)   * (1._r8 - decay_const * dt)
         deadstemc_xfer(p)     = deadstemc_xfer(p)      * (1._r8 - decay_const * dt)
         frootc(p)             = frootc(p)              * (1._r8 - decay_const * dt)
         frootc_storage(p)     = frootc_storage(p)      * (1._r8 - decay_const * dt)
         frootc_xfer(p)        = frootc_xfer(p)         * (1._r8 - decay_const * dt)
         gresp_storage(p)      = gresp_storage(p)       * (1._r8 - decay_const * dt)
         gresp_xfer(p)         = gresp_xfer(p)          * (1._r8 - decay_const * dt)
         leafc(p)              = leafc(p)               * (1._r8 - decay_const * dt)
         leafc_storage(p)      = leafc_storage(p)       * (1._r8 - decay_const * dt)
         leafc_xfer(p)         = leafc_xfer(p)          * (1._r8 - decay_const * dt)
         livecrootc(p)         = livecrootc(p)          * (1._r8 - decay_const * dt)
         livecrootc_storage(p) = livecrootc_storage(p)  * (1._r8 - decay_const * dt)
         livecrootc_xfer(p)    = livecrootc_xfer(p)     * (1._r8 - decay_const * dt)
         livestemc(p)          = livestemc(p)           * (1._r8 - decay_const * dt)
         livestemc_storage(p)  = livestemc_storage(p)   * (1._r8 - decay_const * dt)
         livestemc_xfer(p)     = livestemc_xfer(p)      * (1._r8 - decay_const * dt)
         pft_ctrunc(p)         = pft_ctrunc(p)          * (1._r8 - decay_const * dt)
      end do

    end associate

  end subroutine C14Decay

end module CNC14DecayMod
