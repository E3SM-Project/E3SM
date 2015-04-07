module EDAccumulateFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This routine accumulates NPP, GPP and respiration of each cohort over the course of each 24 hour period. 
  ! The fluxes are stored per cohort, and the npp_clm (etc) fluxes are calcualted in EDPhotosynthesis
  ! This routine cannot be in EDPhotosynthesis because EDPhotosynthesis is a loop and therefore would
  ! erroneously add these things up multiple times. 
  ! Rosie Fisher. March 2014. 
  !
  ! !USES:
  implicit none
  !
  public :: AccumulateFluxes_ED
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine AccumulateFluxes_ED(bounds, p, ed_allsites_inst, photosyns_inst)
    !
    ! !DESCRIPTION:
    ! see above
    !
    ! !USES:
    use shr_kind_mod      , only : r8 => shr_kind_r8
    use decompMod         , only : bounds_type
    use EDTypesMod        , only : ed_patch_type, ed_cohort_type, ed_site_type, map_clmpatch_to_edpatch
    use PatchType         , only : patch
    use PhotosynthesisMod , only : photosyns_type
    !
    ! !ARGUMENTS    
    type(bounds_type)    , intent(in)            :: bounds  
    integer              , intent(in)            :: p     !patch/'p'
    type(ed_site_type)   , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    type(photosyns_type) , intent(inout)         :: photosyns_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer  :: currentCohort ! current cohort
    type(ed_patch_type) , pointer  :: currentPatch ! current patch
    integer :: iv !leaf layer
    integer :: g  !gridcell
    !----------------------------------------------------------------------

    associate(& 
         fpsn      => photosyns_inst%fpsn_patch      , & ! Output: [real(r8) (:)]   photosynthesis (umol CO2 /m**2 /s)
         psncanopy => photosyns_inst%psncanopy_patch   & ! Output: [real(r8) (:,:)] canopy scale photosynthesis umol CO2 /m**2/ s
         )

      fpsn(p) = psncanopy(p)

      if (patch%is_veg(p)) then

         g = patch%gridcell(p)
         currentPatch => map_clmpatch_to_edpatch(ed_allsites_inst(g), p) 
         currentCohort => currentPatch%shortest

         do while(associated(currentCohort))

            ! Accumulate fluxes from hourly to daily values. 
            ! _clm fluxes are KgC/indiv/timestep _acc are KgC/indiv/day

            currentCohort%npp_acc  = currentCohort%npp_acc  + currentCohort%npp_clm 
            currentCohort%gpp_acc  = currentCohort%gpp_acc  + currentCohort%gpp_clm 
            currentCohort%resp_acc = currentCohort%resp_acc + currentCohort%resp_clm

            do iv=1,currentCohort%nv
               if(currentCohort%year_net_uptake(iv) == 999._r8)then ! note that there were leaves in this layer this year. 
                  currentCohort%year_net_uptake(iv) = 0._r8
               end if
               currentCohort%year_net_uptake(iv) = currentCohort%year_net_uptake(iv) + currentCohort%ts_net_uptake(iv)
            enddo

            currentCohort => currentCohort%taller
         enddo ! while(associated(currentCohort)

      end if !is_veg

    end associate

  end subroutine AccumulateFluxes_ED

end module EDAccumulateFluxesMod
