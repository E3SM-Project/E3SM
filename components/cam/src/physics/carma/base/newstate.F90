! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine manages the calculations that update state variables
!! of the model with new values at the current simulation time.
!!
!! @author Bardeen
!! @version Jan 2012
subroutine newstate(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  
  real(kind=f)                    :: pc_orig(NZ,NBIN,NELEM)
  real(kind=f)                    :: gc_orig(NZ,NGAS)
  real(kind=f)                    :: t_orig(NZ)
  real(kind=f)                    :: cldfrc_orig(NZ)
  real(kind=f)                    :: scale_cldfrc(NZ)
  real(kind=f)                    :: pc_cloudy(NZ,NBIN,NELEM)
  real(kind=f)                    :: gc_cloudy(NZ,NGAS)
  real(kind=f)                    :: t_cloudy(NZ)
  real(kind=f)                    :: rlheat_cloudy(NZ)
  real(kind=f)                    :: partheat_cloudy(NZ)
  real(kind=f)                    :: zsubsteps_cloudy(NZ)
  real(kind=f)                    :: pc_clear(NZ,NBIN,NELEM)
  real(kind=f)                    :: gc_clear(NZ,NGAS)
  real(kind=f)                    :: t_clear(NZ)
  real(kind=f)                    :: rlheat_clear(NZ)
  real(kind=f)                    :: partheat_clear(NZ)
  real(kind=f)                    :: zsubsteps_clear(NZ)
  real(kind=f)                    :: scale_threshold(NZ)
  integer                         :: igroup
  integer                         :: igas
  integer                         :: ielem
  integer                         :: ibin
  integer                         :: iz


  ! Calculate changes due to vertical transport
  if (do_vtran) then
  
    call vertical(carma, cstate, rc)
    if (rc < RC_OK) return
  endif

  
  ! There can be two phases to the microphysics: in-cloud and clear sky. Particles
  ! that are tagged as "In-cloud" will only be processed in the in-cloud loop, and their
  ! concentrations will be scaled by the cloud fraction since it is assumed to be all
  ! in-cloud. Other particle types will be process in-cloud and out of cloud; however,
  ! their mass is assumed to be a gridbox average.

  ! If doing doing in-cloud processing, then scale the parameters for in-cloud concentrations.
  ! 
  ! NOTE: Don't want to do this before sedimentation, since sedimentation doesn't take into
  ! account the varying cloud fractions, and thus a particle scaled at one level and cloud
  ! fraction would be scaled inappropriately at another level and cloud fraction.
  !
  ! NOTE: All detrainment also happens only in the in-cloud portion.
  if (do_incloud) then
  
    ! First do the in-cloud processing.
  
    ! Convert "cloud" particles to in-cloud values.
    !
    ! NOTE: If a particle is a "cloud" particle, it means that the entire mass of the
    ! particle is in the incloud portion of the grid box. Particle that are not "cloud
    ! particles" have their mass spread throughout the grid box.
    pc_orig(:,:,:) = pc(:,:,:)
    gc_orig(:,:)   = gc(:,:)
    t_orig(:)      = t(:)
    
    ! If the cloud fraction gets too small it causes the microphysics to require a
    ! lot of substeps. Enforce a minimum cloud fraction for the purposes of scaling
    ! to incloud values.
    scale_cldfrc(:) = max(CLDFRC_MIN, cldfrc(:))
    scale_cldfrc(:) = min(1._f - CLDFRC_MIN, scale_cldfrc(:))

    do ielem = 1, NELEM
      igroup = igelem(ielem)
      
      if (is_grp_cloud(igroup)) then
        do ibin = 1, NBIN
          pc(:, ibin, ielem)  = pc(:, ibin, ielem)  / scale_cldfrc(:)
          pcd(:, ibin, ielem) = pcd(:, ibin, ielem) / scale_cldfrc(:)
        end do
      end if
    end do

    call newstate_calc(carma, cstate, scale_cldfrc(:), rc)
    if (rc < RC_OK) return
    
    ! Save the new in-cloud values for the gas, particle and temperature fields.
    pc_cloudy(:,:,:)    = pc(:,:,:)
    gc_cloudy(:,:)      = gc(:,:)
    t_cloudy(:)         = t(:)
    rlheat_cloudy(:)    = rlheat(:)
    partheat_cloudy(:)  = partheat(:)
    
    if (do_substep) zsubsteps_cloudy(:) = zsubsteps(:)



    ! Now do the clear sky portion, using the original gridbox average concentrations.
    ! This is optional. If clear sky is not selected then all of the microphysics is
    ! done in-cloud.
    pc(:,:,:) = pc_orig(:,:,:)
    gc(:,:)   = gc_orig(:,:)
    t(:)      = t_orig(:)

    if (do_clearsky) then
      
      ! Convert "cloud" particles to clear sky values.
      !
      ! NOTE: If a particle is a "cloud" particle, it means that the entire mass of the
      ! particle is in the in-cloud portion of the grid box. They have no mass in the
      ! clear sky portion.
      do ielem = 1, NELEM
        igroup = igelem(ielem)
        
        if (is_grp_cloud(igroup)) then
          pc(:, :, ielem)  = 0._f
          pcd(:, :, ielem) = 0._f
        end if
      end do
      
      ! Don't let the supersaturation be scaled by setting the cloud fraction used
      ! by the saturation code to 1.0. Any clouds formed in-situ in the clear sky
      ! are assumed to fill the grid box.
      cldfrc_orig(:) = cldfrc(:)
      cldfrc(:)      = 1._f
  
      ! Recalculate supersaturation.
      do igas = 1, NGAS
        do iz = 1, NZ
          call supersat(carma, cstate, iz, igas, rc)
          if (rc < RC_OK) return
        end do
      end do

      call newstate_calc(carma, cstate, (1._f - scale_cldfrc(:)), rc)
      if (rc < RC_OK) return
  
      ! Restore the cloud fraction
      cldfrc(:) = cldfrc_orig(:)
      
      ! Save the new clear sky values for the gas, particle and temperature fields.
      pc_clear(:,:,:)     = pc(:,:,:)
      gc_clear(:,:)       = gc(:,:)
      t_clear(:)          = t(:)
      rlheat_clear(:)     = rlheat(:)
      partheat_clear(:)   = partheat(:)

      if (do_substep) zsubsteps_clear(:) = zsubsteps(:)

    ! If not doing a clear sky calculation, then the clear sky portion reamins
    ! the same except for any contribution from advection.
    else
    
      ! NOTE: If a particle is a "cloud" particle, it means that the entire mass of the
      ! particle is in the in-cloud portion of the grid box. They have no mass in the
      ! clear sky portion.
      do ielem = 1, NELEM
        igroup = igelem(ielem)
        
        if (is_grp_cloud(igroup)) then
          pc_clear(:, :, ielem)  = 0._f
        else
          pc_clear(:, :, ielem)  =  pc(:, :, ielem)
        end if
      end do

      do igas = 1, NGAS
        gc_clear(:,:)     = gc(:,:)
      end do
      
      t_clear(:)          = t(:)
      rlheat_clear(:)     = 0._f
      partheat_clear(:)   = 0._f

      ! If substepping, then add the advected part that is being doled out over
      ! the substeps.
      if (do_substep) then
        do igas = 1, NGAS
          gc_clear(:, igas)     = gc_clear(:, igas) + d_gc(:, igas)
        end do
        t_clear(:)          = t_clear(:) + d_t(:)

        zsubsteps_clear(:) = 0._f
      end if
    end if
    
    
    ! Add up the changes to the particle from the cloudy and clear sky components.
    do ielem = 1, NELEM
      igroup = igelem(ielem)
      
      do ibin = 1, NBIN
        pc(:, ibin, ielem)   = (1._f - scale_cldfrc(:)) * pc_clear(:, ibin, ielem) + scale_cldfrc(:) * pc_cloudy(:, ibin, ielem)
      end do
    end do
        
    t(:) = (1._f - scale_cldfrc(:)) * t_clear(:) + scale_cldfrc(:) * t_cloudy(:)
    
    if (do_grow) then
      rlheat(:)   = (1._f - scale_cldfrc(:)) * rlheat_clear(:)   + scale_cldfrc(:) * rlheat_cloudy(:)
      partheat(:) = (1._f - scale_cldfrc(:)) * partheat_clear(:) + scale_cldfrc(:) * partheat_cloudy(:)
    end if

    do igas = 1, NGAS
      gc(:, igas) = (1._f - scale_cldfrc(:)) * gc_clear(:, igas) + scale_cldfrc(:) * gc_cloudy(:, igas)

      ! Recalculate gridbox average supersaturation.
      do iz = 1, NZ
        call supersat(carma, cstate, iz, igas, rc)
        if (rc < RC_OK) return
      end do
    end do

    if (do_substep) zsubsteps(:) = zsubsteps_clear(:) + zsubsteps_cloudy(:)
  
  
  
  ! No special in-cloud/clear sky processing, everything is gridbox average.
  else
    scale_threshold(:) = 1._f
    call newstate_calc(carma, cstate, scale_threshold, rc)
    if (rc < RC_OK) return
  end if
    
  ! Return to caller with new state computed 
  return
end
