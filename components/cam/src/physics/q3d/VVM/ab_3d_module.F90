MODULE ab_3d_module
! Controls the integration of individual CRM (one timestep calculation)

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: channel_t

USE CONSTLD,        only: notherm,buoy,physics,turbulence
      
! Subroutines being called
USE vort_3d_module,       only: vort_3d,vort_3d_corec
USE wind_module,          only: wind_3d
USE turb_3d_module,       only: turb_3d,turb_3d_therm,turb_3d_vort
USE buoyf_module,         only: buoyf_3d
USE q_chk_module,         only: q_chk_3d
USE update_thermo_module, only: update_thermodynamics,revise_thermodynamics
USE rcalc_module,         only: rcalc_3d

IMPLICIT NONE
PRIVATE

PUBLIC :: ab_3d

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE ab_3d
!-----------------------------------------------------------------------
!=======================================================================
   SUBROUTINE AB_3D (N1, N2, ITT, channel)
!=======================================================================
      type(channel_t), intent(inout) :: channel   ! channel data
      
      INTEGER, INTENT(IN) :: &
         ITT,       & ! time step count
         N1,        & ! AB forcing time index for previous timestep
         N2           ! AB forcing time index for current timestep

!     Predict thermodynamic fields (advection, microphysics, radiation)
!     ITT is related to how often radiation is updated.  
      IF (.NOT.NOTHERM) CALL RCALC_3D (N1, N2, ITT, channel)

!     Calculate the nonlinear turbulence coefficients
!     ITT is related to how often surface flux is calculated.
      IF (TURBULENCE) CALL TURB_3D (channel)

      IF (.NOT.NOTHERM) THEN

!     Calculate the thermodynamic tendencies due to turbulence
      IF (TURBULENCE) CALL TURB_3D_THERM(channel)

!     Re-update the thermodynamic variables
      CALL REVISE_THERMODYNAMICS(channel)

      IF (PHYSICS) THEN
!       Fill negative values & Saturation adjustment
        CALL Q_CHK_3D (channel)
      ENDIF

      ENDIF

!     Calculate the buoyancy based on the predicted thermodynamic fields 
      CALL BUOYF_3D (BUOY, channel)

!     Predict the vorticity components
      CALL VORT_3D  (N1, N2, channel)

!     Calculate the vorticity tendencies due to turbulence
      IF (TURBULENCE) CALL TURB_3D_VORT(channel)

!     Re-update the vorticity components
      CALL VORT_3D_COREC(channel)

!     Diagnose the wind components
      CALL WIND_3D  (N1, N2, .TRUE., .TRUE., channel)  ! UCHANGE=.T., WCHANGE=.T. 

   END SUBROUTINE ab_3d

END MODULE ab_3d_module
