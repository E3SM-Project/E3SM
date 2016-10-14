
module RichardsMod

#ifdef USE_PETSC_LIB

  ! !USES:
  use mpp_varctl  , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: RichardsFlux
  public :: RichardsFluxDerivativeWrtTemperature

  !------------------------------------------------------------------------------

contains

  subroutine RichardsFlux(Pres_up, &
       kr_up,                      &
       dkr_dP_up,                  &
       den_up,                     &
       dden_dP_up,                 &
       vis_up,                     &
       dvis_dP_up,                 &
       perm_vec_up,                &
       Pres_dn,                    &
       kr_dn,                      &
       dkr_dP_dn,                  &
       den_dn,                     &
       dden_dP_dn,                 &
       vis_dn,                     &
       dvis_dP_dn,                 &
       perm_vec_dn,                &
       area,                       &
       dist_up,                    &
       dist_dn,                    &
       dist_unitvec,               &
       compute_deriv,              &
       internal_conn,              &
       cond_type,                  &
       flux,                       &
       dflux_dP_up,                &
       dflux_dP_dn                 &
       )

    !
    ! !DESCRIPTION:
    ! Based on the primary (P) and secondary (kr, den, etc) values of upwind and downwind
    ! control volumes, this subroutine computes:
    !  - Two-point flux for Richards equation, and
    !  - (optinal) Derivative of the flux w.r.t. to upwind and downwind pressure.
    !
    ! Negative flux implies flow occurs from upwind to downwind control volume.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GRAVITY_CONSTANT, PRESSURE_REF
    use MultiPhysicsProbConstants, only : COND_DIRICHLET, COND_MASS_FLUX, COND_MASS_RATE
    use MultiPhysicsProbConstants, only : COND_DIRICHLET_FRM_OTR_GOVEQ, COND_SEEPAGE_BC
    use MultiPhysicsProbConstants, only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: Pres_up
    PetscReal, intent(in)  :: kr_up
    PetscReal, intent(in)  :: dkr_dP_up
    PetscReal, intent(in)  :: den_up
    PetscReal, intent(in)  :: dden_dP_up
    PetscReal, intent(in)  :: vis_up
    PetscReal, intent(in)  :: dvis_dP_up
    PetscReal, intent(in)  :: perm_vec_up(3)
    PetscReal, intent(in)  :: Pres_dn
    PetscReal, intent(in)  :: kr_dn
    PetscReal, intent(in)  :: dkr_dP_dn
    PetscReal, intent(in)  :: den_dn
    PetscReal, intent(in)  :: dden_dP_dn
    PetscReal, intent(in)  :: vis_dn
    PetscReal, intent(in)  :: dvis_dP_dn
    PetscReal, intent(in)  :: perm_vec_dn(3)
    PetscReal, intent(in)  :: area
    PetscReal, intent(in)  :: dist_up
    PetscReal, intent(in)  :: dist_dn
    PetscReal, intent(in)  :: dist_unitvec(3)
    PetscBool, intent(in)  :: compute_deriv
    PetscBool, intent(in)  :: internal_conn
    PetscInt,  intent(in)  :: cond_type
    PetscReal, intent(out) :: flux
    PetscReal, intent(out) :: dflux_dP_up
    PetscReal, intent(out) :: dflux_dP_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: upweight
    PetscReal :: Dq
    PetscReal :: grav_vec(3)
    PetscReal :: udist_dot_ugrav
    PetscReal :: dist_gravity
    PetscReal :: gravityterm
    PetscReal :: den_ave
    PetscReal :: dphi
    PetscReal :: ukvr
    PetscReal :: v_darcy
    PetscReal :: q
    PetscReal :: perm_up
    PetscReal :: perm_dn

    PetscReal :: dden_ave_dP_up
    PetscReal :: dden_ave_dP_dn
    PetscReal :: dgravityterm_dden_up
    PetscReal :: dgravityterm_dden_dn
    PetscReal :: dphi_dP_up
    PetscReal :: dphi_dP_dn
    PetscReal :: dukvr_dP_up
    PetscReal :: dukvr_dP_dn
    PetscReal :: dq_dP_up
    PetscReal :: dq_dP_dn

    PetscBool :: seepage_bc_update

    perm_up = dabs(dist_unitvec(1))*perm_vec_up(1) + &
         dabs(dist_unitvec(2))*perm_vec_up(2) + &
         dabs(dist_unitvec(3))*perm_vec_up(3)

    perm_dn = dabs(dist_unitvec(1))*perm_vec_dn(1) + &
         dabs(dist_unitvec(2))*perm_vec_dn(2) + &
         dabs(dist_unitvec(3))*perm_vec_dn(3)

    if (internal_conn) then
       upweight = dist_up/(dist_up + dist_dn)
       Dq       = (perm_up * perm_dn)/(dist_up*perm_dn + dist_dn*perm_up)
    else

       select case(cond_type)
       case (COND_DIRICHLET, COND_MASS_FLUX, COND_SEEPAGE_BC)
          upweight = 0.d0
          Dq       = perm_dn/(dist_up + dist_dn)
       case (COND_DIRICHLET_FRM_OTR_GOVEQ)
          upweight = dist_up/(dist_up + dist_dn)
          Dq       = (perm_up * perm_dn)/(dist_up*perm_dn + dist_dn*perm_up)
       case default
          write(iulog,*)'RichardsFlux: Unknown cond_type Add additional code.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    endif

    ! Gravity effects.
    !   Note make dot-operation explicit, to emphasize gravity vector does not
    ! need to align with vertical coordinate.
    grav_vec(1) = 0.d0
    grav_vec(2) = 0.d0
    grav_vec(3) = -GRAVITY_CONSTANT
    udist_dot_ugrav = dist_unitvec(1)*grav_vec(1) + &
         dist_unitvec(2)*grav_vec(2) + &
         dist_unitvec(3)*grav_vec(3)

    dist_gravity = (dist_up + dist_dn)*udist_dot_ugrav

    den_ave      = upweight*den_up + (1.d0 - upweight)*den_dn

    ! gravityterm = (rho*g*z)
    gravityterm  = (upweight*den_up + (1.d0 - upweight)*den_dn)*FMWH2O*dist_gravity;

    dphi         = Pres_up - Pres_dn + gravityterm;

    seepage_bc_update = PETSC_FALSE
    if (.not.(internal_conn) .and. (cond_type == COND_SEEPAGE_BC)) then
       if (dphi > 0.d0 .and. Pres_up <= PRESSURE_REF) then
          seepage_bc_update = PETSC_TRUE
       endif
    endif
    if (seepage_bc_update) dphi = 0.d0

    if (dphi>=0.d0) then
       ukvr = kr_up/vis_up
    else
       ukvr = kr_dn/vis_dn
    endif

    if (.not.(internal_conn) .and. (cond_type == COND_MASS_FLUX)) then
       v_darcy = 0.d0
    else
       v_darcy = -Dq * ukvr * dphi
    endif


    q    = v_darcy * area
    flux = q * den_ave

    if (compute_deriv) then

       dden_ave_dP_up       = upweight*dden_dP_up
       dden_ave_dP_dn       = (1.d0 - upweight)*dden_dP_dn

       dgravityterm_dden_up = upweight*dist_gravity*FMWH2O
       dgravityterm_dden_dn = (1.d0-upweight)*dist_gravity*FMWH2O

       dphi_dP_up           =  1.d0 + dgravityterm_dden_up*dden_dP_up;
       dphi_dP_dn           = -1.d0 + dgravityterm_dden_dn*dden_dP_dn;

       if (seepage_bc_update) dphi_dP_dn = 0.d0

       if (dphi>=0) then
          dukvr_dP_up = dkr_dP_up/vis_up - kr_up/(vis_up*vis_up)*dvis_dP_up
          dukvr_dP_dn = 0.d0
       else
          dukvr_dP_up = 0.d0
          dukvr_dP_dn = dkr_dP_dn/vis_dn - kr_dn/(vis_dn*vis_dn)*dvis_dP_dn
       endif

       dq_dP_up = Dq*(dukvr_dP_up*dphi + ukvr*dphi_dP_up)*area;
       dq_dP_dn = Dq*(dukvr_dP_dn*dphi + ukvr*dphi_dP_dn)*area;

       if (.not.(internal_conn) .and. (cond_type == COND_MASS_FLUX)) then
          dflux_dP_up = 0.d0
          dflux_dP_dn = 0.d0
       else
          dflux_dP_up = (dq_dp_up*den_ave - q*dden_ave_dp_up);
          dflux_dP_dn = (dq_dp_dn*den_ave - q*dden_ave_dp_dn);

       endif

    endif

  end subroutine RichardsFlux

  !------------------------------------------------------------------------------
  subroutine RichardsFluxDerivativeWrtTemperature(Pres_up, &
       kr_up,                      &
       den_up,                     &
       dden_dT_up,                 &
       vis_up,                     &
       dvis_dT_up,                 &
       perm_vec_up,                &
       Pres_dn,                    &
       kr_dn,                      &
       den_dn,                     &
       dden_dT_dn,                 &
       vis_dn,                     &
       dvis_dT_dn,                 &
       perm_vec_dn,                &
       area,                       &
       dist_up,                    &
       dist_dn,                    &
       dist_unitvec,               &
       compute_deriv,              &
       internal_conn,              &
       cond_type,                  &
       flux,                       &
       dflux_dT_up,                &
       dflux_dT_dn                 &
       )
    !
    ! !DESCRIPTION:
    ! Based on the primary (P) and secondary (kr, den, etc) values of upwind and downwind
    ! control volumes, this subroutine computes:
    !  - Two-point flux for Richards equation, and
    !  - (optinal) Derivative of the flux w.r.t. to upwind and downwind pressure.
    !
    ! Negative flux implies flow occurs from upwind to downwind control volume.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GRAVITY_CONSTANT, PRESSURE_REF
    use MultiPhysicsProbConstants, only : COND_DIRICHLET, COND_MASS_FLUX, COND_MASS_RATE
    use MultiPhysicsProbConstants, only : COND_DIRICHLET_FRM_OTR_GOVEQ, COND_SEEPAGE_BC
    use MultiPhysicsProbConstants, only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: Pres_up
    PetscReal, intent(in)  :: kr_up
    PetscReal, intent(in)  :: den_up
    PetscReal, intent(in)  :: dden_dT_up
    PetscReal, intent(in)  :: vis_up
    PetscReal, intent(in)  :: dvis_dT_up
    PetscReal, intent(in)  :: perm_vec_up(3)
    PetscReal, intent(in)  :: Pres_dn
    PetscReal, intent(in)  :: kr_dn
    PetscReal, intent(in)  :: den_dn
    PetscReal, intent(in)  :: dden_dT_dn
    PetscReal, intent(in)  :: vis_dn
    PetscReal, intent(in)  :: dvis_dT_dn
    PetscReal, intent(in)  :: perm_vec_dn(3)
    PetscReal, intent(in)  :: area
    PetscReal, intent(in)  :: dist_up
    PetscReal, intent(in)  :: dist_dn
    PetscReal, intent(in)  :: dist_unitvec(3)
    PetscBool, intent(in)  :: compute_deriv
    PetscBool, intent(in)  :: internal_conn
    PetscInt,  intent(in)  :: cond_type
    PetscReal, intent(out) :: flux
    PetscReal, intent(out) :: dflux_dT_up
    PetscReal, intent(out) :: dflux_dT_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: upweight
    PetscReal :: Dq
    PetscReal :: grav_vec(3)
    PetscReal :: udist_dot_ugrav
    PetscReal :: dist_gravity
    PetscReal :: gravityterm
    PetscReal :: den_ave
    PetscReal :: dphi
    PetscReal :: ukvr
    PetscReal :: v_darcy
    PetscReal :: q
    PetscReal :: perm_up
    PetscReal :: perm_dn

    PetscReal :: dden_ave_dT_up
    PetscReal :: dden_ave_dT_dn
    PetscReal :: dgravityterm_dden_up
    PetscReal :: dgravityterm_dden_dn
    PetscReal :: dphi_dT_up
    PetscReal :: dphi_dT_dn
    PetscReal :: dukvr_dT_up
    PetscReal :: dukvr_dT_dn
    PetscReal :: dq_dT_up
    PetscReal :: dq_dT_dn

    PetscBool :: seepage_bc_update

    perm_up = dabs(dist_unitvec(1))*perm_vec_up(1) + &
         dabs(dist_unitvec(2))*perm_vec_up(2) + &
         dabs(dist_unitvec(3))*perm_vec_up(3)

    perm_dn = dabs(dist_unitvec(1))*perm_vec_dn(1) + &
         dabs(dist_unitvec(2))*perm_vec_dn(2) + &
         dabs(dist_unitvec(3))*perm_vec_dn(3)

    if (internal_conn) then
       upweight = dist_up/(dist_up + dist_dn)
       Dq       = (perm_up * perm_dn)/(dist_up*perm_dn + dist_dn*perm_up)
    else

       select case(cond_type)
       case (COND_DIRICHLET, COND_MASS_FLUX, COND_SEEPAGE_BC)
          upweight = 0.d0
          Dq       = perm_dn/(dist_up + dist_dn)
       case (COND_DIRICHLET_FRM_OTR_GOVEQ)
          upweight = dist_up/(dist_up + dist_dn)
          Dq       = (perm_up * perm_dn)/(dist_up*perm_dn + dist_dn*perm_up)
       case default
          write(iulog,*)'RichardsFlux: Unknown cond_type Add additional code.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    endif

    ! Gravity effects.
    !   Note make dot-operation explicit, to emphasize gravity vector does not
    ! need to align with vertical coordinate.
    grav_vec(1) = 0.d0
    grav_vec(2) = 0.d0
    grav_vec(3) = -GRAVITY_CONSTANT
    udist_dot_ugrav = dist_unitvec(1)*grav_vec(1) + &
         dist_unitvec(2)*grav_vec(2) + &
         dist_unitvec(3)*grav_vec(3)

    dist_gravity = (dist_up + dist_dn)*udist_dot_ugrav

    den_ave      = upweight*den_up + (1.d0 - upweight)*den_dn

    ! gravityterm = (rho*g*z)
    gravityterm  = (upweight*den_up + (1.d0 - upweight)*den_dn)*FMWH2O*dist_gravity;

    dphi         = Pres_up - Pres_dn + gravityterm;

    seepage_bc_update = PETSC_FALSE
    if (.not.(internal_conn) .and. (cond_type == COND_SEEPAGE_BC)) then
       if (dphi > 0.d0 .and. Pres_up <= PRESSURE_REF) then
          seepage_bc_update = PETSC_TRUE
       endif
    endif
    if (seepage_bc_update) dphi = 0.d0

    if (dphi>=0.d0) then
       ukvr = kr_up/vis_up
    else
       ukvr = kr_dn/vis_dn
    endif

    if (.not.(internal_conn) .and. (cond_type == COND_MASS_FLUX)) then
       v_darcy = 0.d0
    else
       v_darcy = -Dq * ukvr * dphi
    endif

    q    = v_darcy * area
    flux = q * den_ave

    if (compute_deriv) then

       dden_ave_dT_up       = upweight*dden_dT_up
       dden_ave_dT_dn       = (1.d0 - upweight)*dden_dT_dn

       dgravityterm_dden_up = upweight*dist_gravity*FMWH2O
       dgravityterm_dden_dn = (1.d0-upweight)*dist_gravity*FMWH2O

       dphi_dT_up           = dgravityterm_dden_up*dden_dT_up
       dphi_dT_dn           = dgravityterm_dden_dn*dden_dT_dn

       if (seepage_bc_update) dphi_dT_dn = 0.d0

       if (dphi>=0) then
          dukvr_dT_up = - kr_up/(vis_up*vis_up)*dvis_dT_up
          dukvr_dT_dn = 0.d0
       else
          dukvr_dT_up = 0.d0
          dukvr_dT_dn = - kr_dn/(vis_dn*vis_dn)*dvis_dT_dn
       endif

       dq_dT_up = Dq*(dukvr_dT_up*dphi + ukvr*dphi_dT_up)*area;
       dq_dT_dn = Dq*(dukvr_dT_dn*dphi + ukvr*dphi_dT_dn)*area;

       if (.not.(internal_conn) .and. (cond_type == COND_MASS_FLUX)) then
          dflux_dT_up = 0.d0
          dflux_dT_dn = 0.d0
       else
          dflux_dT_up = (dq_dT_up*den_ave - q*dden_ave_dT_up);
          dflux_dT_dn = (dq_dT_dn*den_ave - q*dden_ave_dT_dn);
       endif

    endif

    dflux_dT_up = -dflux_dT_up
    dflux_dT_dn = -dflux_dT_dn

  end subroutine RichardsFluxDerivativeWrtTemperature

#endif

end module RichardsMod
