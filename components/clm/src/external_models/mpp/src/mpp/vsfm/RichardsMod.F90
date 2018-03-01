
module RichardsMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>
  use petscsys

  ! !USES:
  use mpp_varctl  , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: RichardsFlux
  public :: RichardsFluxDerivativeWrtTemperature
  public :: RichardsFluxConductanceModel
  !
  !------------------------------------------------------------------------------

contains

  subroutine RichardsFlux( &
       aux_var_up,         &
       aux_var_dn,         &
       conn_up2dn,         &
       compute_deriv,      &
       internal_conn,      &
       cond_type,          &
       flux,               &
       dflux_dP_up,        &
       dflux_dP_dn         &
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
    use MultiPhysicsProbConstants  , only : GRAVITY_CONSTANT, PRESSURE_REF
    use MultiPhysicsProbConstants  , only : COND_DIRICHLET, COND_MASS_FLUX, COND_MASS_RATE
    use MultiPhysicsProbConstants  , only : COND_DIRICHLET_FRM_OTR_GOVEQ, COND_SEEPAGE_BC
    use MultiPhysicsProbConstants  , only : FMWH2O
    use ConnectionSetType          , only : connection_type
    use RichardsODEPressureAuxType , only : rich_ode_pres_auxvar_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) , intent(in)  :: aux_var_up
    class(rich_ode_pres_auxvar_type) , intent(in)  :: aux_var_dn
    type(connection_type)            , intent(in)  :: conn_up2dn
    PetscBool                        , intent(in)  :: compute_deriv
    PetscBool                        , intent(in)  :: internal_conn
    PetscInt                         , intent(in)  :: cond_type
    PetscReal                        , intent(out) :: flux
    PetscReal                        , intent(out) :: dflux_dP_up
    PetscReal                        , intent(out) :: dflux_dP_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: area
    PetscReal :: dist_up
    PetscReal :: dist_dn
    PetscReal :: dist_unitvec(3)
    PetscReal :: Pres_up
    PetscReal :: kr_up
    PetscReal :: dkr_dP_up
    PetscReal :: den_up
    PetscReal :: dden_dP_up
    PetscReal :: vis_up
    PetscReal :: dvis_dP_up
    PetscReal :: perm_vec_up(3)
    PetscReal :: Pres_dn
    PetscReal :: kr_dn
    PetscReal :: dkr_dP_dn
    PetscReal :: den_dn
    PetscReal :: dden_dP_dn
    PetscReal :: vis_dn
    PetscReal :: dvis_dP_dn
    PetscReal :: perm_vec_dn(3)
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

    ! Get geometric attributes about connection
    area            = conn_up2dn%GetArea()
    dist_up         = conn_up2dn%GetDistUp()
    dist_dn         = conn_up2dn%GetDistDn()
    dist_unitvec(1) = conn_up2dn%GetDistUnitVecX()
    dist_unitvec(2) = conn_up2dn%GetDistUnitVecY()
    dist_unitvec(3) = conn_up2dn%GetDistUnitVecZ()

    ! Get variables related to upwind cell
    Pres_up         = aux_var_up%pressure
    kr_up           = aux_var_up%kr
    dkr_dP_up       = aux_var_up%dkr_dP
    den_up          = aux_var_up%den
    dden_dP_up      = aux_var_up%dden_dP
    vis_up          = aux_var_up%vis
    dvis_dP_up      = aux_var_up%dvis_dP
    perm_vec_up     = aux_var_up%perm

    ! Get variables related to downwind cell
    Pres_dn         = aux_var_dn%pressure
    kr_dn           = aux_var_dn%kr
    dkr_dP_dn       = aux_var_dn%dkr_dP
    den_dn          = aux_var_dn%den
    dden_dP_dn      = aux_var_dn%dden_dP
    vis_dn          = aux_var_dn%vis
    dvis_dP_dn      = aux_var_dn%dvis_dP
    perm_vec_dn     = aux_var_dn%perm

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
  subroutine RichardsFluxDerivativeWrtTemperature( &
       aux_var_up,                                 &
       aux_var_dn,                                 &
       conn_up2dn,                                 &
       compute_deriv,                              &
       internal_conn,                              &
       cond_type,                                  &
       flux,                                       &
       dflux_dT_up,                                &
       dflux_dT_dn                                 &
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
    use MultiPhysicsProbConstants  , only : GRAVITY_CONSTANT, PRESSURE_REF
    use MultiPhysicsProbConstants  , only : COND_DIRICHLET, COND_MASS_FLUX, COND_MASS_RATE
    use MultiPhysicsProbConstants  , only : COND_DIRICHLET_FRM_OTR_GOVEQ, COND_SEEPAGE_BC
    use MultiPhysicsProbConstants  , only : FMWH2O
    use ConnectionSetType          , only : connection_type
    use RichardsODEPressureAuxType , only : rich_ode_pres_auxvar_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) , intent(in)  :: aux_var_up
    class(rich_ode_pres_auxvar_type) , intent(in)  :: aux_var_dn
    type(connection_type)            , intent(in)  :: conn_up2dn
    PetscBool                        , intent(in)  :: compute_deriv
    PetscBool                        , intent(in)  :: internal_conn
    PetscInt                         , intent(in)  :: cond_type
    PetscReal                        , intent(out) :: flux
    PetscReal                        , intent(out) :: dflux_dT_up
    PetscReal                        , intent(out) :: dflux_dT_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: area
    PetscReal :: dist_up
    PetscReal :: dist_dn
    PetscReal :: dist_unitvec(3)
    PetscReal :: Pres_up
    PetscReal :: kr_up
    PetscReal :: den_up
    PetscReal :: dden_dT_up
    PetscReal :: vis_up
    PetscReal :: dvis_dT_up
    PetscReal :: perm_vec_up(3)
    PetscReal :: Pres_dn
    PetscReal :: kr_dn
    PetscReal :: den_dn
    PetscReal :: dden_dT_dn
    PetscReal :: vis_dn
    PetscReal :: dvis_dT_dn
    PetscReal :: perm_vec_dn(3)
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

    area            = conn_up2dn%GetArea()
    dist_up         = conn_up2dn%GetDistUp()
    dist_dn         = conn_up2dn%GetDistDn()
    dist_unitvec(1) = conn_up2dn%GetDistUnitVecX()
    dist_unitvec(2) = conn_up2dn%GetDistUnitVecY()
    dist_unitvec(3) = conn_up2dn%GetDistUnitVecZ()

    ! Get variables related to upwind cell
    Pres_up         = aux_var_up%pressure
    kr_up           = aux_var_up%kr
    den_up          = aux_var_up%den
    dden_dT_up      = aux_var_up%dden_dT
    vis_up          = aux_var_up%vis
    dvis_dT_up      = aux_var_up%dvis_dT
    perm_vec_up     = aux_var_up%perm

    ! Get variables related to downwind cell
    Pres_dn         = aux_var_dn%pressure
    kr_dn           = aux_var_dn%kr
    den_dn          = aux_var_dn%den
    dden_dT_dn      = aux_var_dn%dden_dT
    vis_dn          = aux_var_dn%vis
    dvis_dT_dn      = aux_var_dn%dvis_dT
    perm_vec_dn     = aux_var_dn%perm
    
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

  !------------------------------------------------------------------------------
  subroutine RichardsFluxConductanceModel( &
       aux_var_up,                         &
       aux_var_dn,                         &
       aux_var_conn,                       &
       conn_up2dn,                         &
       compute_deriv,                      &
       internal_conn,                      &
       cond_type,                          &
       flux,                               &
       dflux_dP_up,                        &
       dflux_dP_dn                         &
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
    use MultiPhysicsProbConstants      , only : GRAVITY_CONSTANT, PRESSURE_REF
    use MultiPhysicsProbConstants      , only : COND_DIRICHLET, COND_MASS_FLUX, COND_MASS_RATE
    use MultiPhysicsProbConstants      , only : COND_DIRICHLET_FRM_OTR_GOVEQ, COND_SEEPAGE_BC
    use MultiPhysicsProbConstants      , only : FMWH2O
    use ConnectionSetType              , only : connection_type
    use RichardsODEPressureAuxType     , only : rich_ode_pres_auxvar_type
    use RichardsODEPressureConnAuxType , only : rich_ode_pres_conn_auxvar_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type)      , intent(in)  :: aux_var_up
    class(rich_ode_pres_auxvar_type)      , intent(in)  :: aux_var_dn
    class(rich_ode_pres_conn_auxvar_type) , intent(in)  :: aux_var_conn
    type(connection_type)                 , intent(in)  :: conn_up2dn
    PetscBool                             , intent(in)  :: compute_deriv
    PetscBool                             , intent(in)  :: internal_conn
    PetscInt                              , intent(in)  :: cond_type
    PetscReal                             , intent(out) :: flux
    PetscReal                             , intent(out) :: dflux_dP_up
    PetscReal                             , intent(out) :: dflux_dP_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: area
    PetscReal :: Pres_up
    PetscReal :: den_up
    PetscReal :: dden_dP_up
    PetscReal :: Pres_dn
    PetscReal :: den_dn
    PetscReal :: dden_dP_dn
    PetscReal :: krg
    PetscReal :: dkrg_dP_up
    PetscReal :: dkrg_dP_dn
    PetscReal :: den_ave
    PetscReal :: upweight
    PetscReal :: dphi
    PetscReal :: dden_ave_dP_up, dden_ave_dP_dn
    PetscReal :: dphi_dP_up, dphi_dP_dn

    upweight    = 0.5d0

    ! Get geometric attribute about connection
    area        = conn_up2dn%GetArea()

    ! Get variables about upwind cell
    Pres_up     = aux_var_up%pressure
    den_up      = aux_var_up%den
    dden_dP_up  = aux_var_up%dden_dP

    ! Get variables about downwind cell
    Pres_dn     = aux_var_dn%pressure
    den_dn      = aux_var_dn%den
    dden_dP_dn  = aux_var_dn%dden_dP

    ! Get variables about connection
    krg         = aux_var_conn%krg
    dkrg_dP_up  = aux_var_conn%dkrg_dP_up
    dkrg_dP_dn  = aux_var_conn%dkrg_dP_dn

    den_ave = upweight*den_up + (1.d0 - upweight)*den_dn
    dphi = (Pres_up - Pres_dn)

    flux    = -den_ave * krg * dphi * area

    if (compute_deriv) then

       dden_ave_dP_up       = upweight*dden_dP_up
       dden_ave_dP_dn       = (1.d0 - upweight)*dden_dP_dn

       dphi_dP_up           =  1.d0
       dphi_dP_dn           = -1.d0

       if (.not.(internal_conn) .and. (cond_type == COND_MASS_FLUX)) then
          dflux_dP_up = 0.d0
          dflux_dP_dn = 0.d0
       else
          dflux_dP_up = + dden_ave_dP_up * krg        * dphi       * area &
                        + den_ave        * dkrg_dP_up * dphi       * area &
                        + den_ave        * krg        * dphi_dP_up * area

          dflux_dP_dn = + dden_ave_dP_dn * krg        * dphi       * area &
                        + den_ave        * dkrg_dP_dn * dphi       * area &
                        + den_ave        * krg        * dphi_dP_dn * area
       endif

    endif

  end subroutine RichardsFluxConductanceModel
#endif

end module RichardsMod
