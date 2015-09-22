#include "carma_globaer.h"

!! This routine handles all preliminary setup at the beginning
!! of every timestep.  Things that would appropriately be done
!! here include:
!!  Input or otherwise define interface quantities from other submodels.
!!  Save any model state that is needed to compute tendencies.
!!  Save any model state that might be needed for comparison at end of step.
!!  Update timestep counter and simulation time.
!!
!! @author Bill McKie
!! @version Oct-1995
subroutine prestep(carma, cstate, rc)

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

  integer                        :: iz      ! z index
  integer                        :: igrp    ! group index
  integer                        :: igas    ! gas index
  integer                        :: ibin    ! bin index
  integer                        :: ielem   ! element index
  integer                        :: iep
  real(kind=f)                   :: tmp_gc(NZ, NGAS)
  real(kind=f)                   :: tmp_t(NZ)


  ! If substepping is enabled, then determine how much the
  ! gas concentration and temperature changed during this time step.
  if (do_substep) then
    if (NGAS > 0) then
      d_gc(:,:) = gc(:,:) - gcl(:,:)
      
      do igas = 1, NGAS
        do iz = 1, NZ
    
          ! NOTE: When d_gc is negative, you can get into problems with overshoot
          ! to negative gas concentrations. To prevent that, when gc is negative
          ! apply it all in the first step. Only substep gc when gc is increasing.
          !
          ! NOTE: Perhaps there should be a limit, so that small changes happen
          ! over the course of the timestep, but large changes get applied on the
          ! first step. For now, doing it all on the first step should be the most
          ! stable.
          !
          ! NOTE: The case that is problematic is when the particle is growing
          ! (i.e. supersaturated) and d_gc is negative. For better performance,
          ! substep the gas unless both of these are true. This might run into
          ! trouble if d_t is large and negative.
!          if (d_gc(iz, igas) < 0._f) then
          if ((d_gc(iz, igas) < 0._f) .and. ((supsatiold(iz, igas) > 0._f) &
               .or. (supsati(iz, igas) > 0._f))) then
          
            ! Start from the new state and don't step the gas.
            d_gc(iz, igas) = 0._f
            gcl(iz, igas) = gc(iz, igas)
          else
    
            ! Start the step from the old state and step the gas.
            gc(iz, igas) = gcl(iz, igas)
          end if
        end do
      end do
    end if

    ! Start the temperature from the old state.
    d_t(:) = t(:) - told(:)
    t(:) = told(:)
  endif

   
  ! Don't allow particle concentrations to get too small.
  do iz = 1, NZ
    do ibin = 1, NBIN
      do ielem = 1, NELEM
        call smallconc(carma, cstate, iz, ibin, ielem, rc)
      end do
    end do
  end do
  
  ! Set <pcl> to <pc> from previous time step. This is needed by coagulation
  ! as well as substepping.
  if (do_substep .or. do_coag) then
    pcl(:,:,:) = pc(:,:,:)
  endif

  ! Find maximum particle concentrations.
  do iz = 1, NZ
    call maxconc(carma, cstate, iz, rc)
  end do

  ! Return to caller with preliminary timestep things completed.
  return
end
