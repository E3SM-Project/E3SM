module InlineSurface_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use InlineSurface_Aux_module
  use Global_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use Option_module

  implicit none

  private


  public ::                           &
       InlineSurfaceAccumulation    , &
       InlineSurfaceFlux            , &
       InlineSurfaceBCFlux          , &
       InlineSurfaceAccumulationJac , &
       InlineSurfaceFluxJac         , &
       InlineSurfaceBCFluxJac

contains

  !************************************************************************** !

  subroutine InlineSurfaceAccumulation(auxvar,material_auxvar,option,Res)
    ! 
    ! Compute the inline surface accumulation term to be used in
    ! Richards mode
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    ! 
    implicit none
    type(inlinesurface_auxvar_type) :: auxvar
    type(material_auxvar_type)      :: material_auxvar
    type(option_type)               :: option
    PetscReal                       :: Res(option%nflowdof),area

    area   = 0.5d0 * material_auxvar%volume / auxvar%half_cell_height
    Res(1) = auxvar%density * auxvar%surface_water_depth / option%flow_dt * area

  end subroutine InlineSurfaceAccumulation

  !************************************************************************** !

  subroutine InlineSurfaceAccumulationJac(auxvar,material_auxvar,option,J)
    ! 
    ! Compute the Jacobian of the inline surface accumulation term to
    ! be used in Richards mode
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    !
    use PFLOTRAN_Constants_module, only : FMWH2O
    implicit none
    type(inlinesurface_auxvar_type) :: auxvar
    type(material_auxvar_type)      :: material_auxvar
    type(option_type)               :: option
    PetscReal                       :: J(option%nflowdof,option%nflowdof),area

    J(1,1) = 0.0d0
    if ( auxvar%surface_water_depth > 0.d0 ) then
      area   = 0.5d0 * material_auxvar%volume / auxvar%half_cell_height
      J(1,1) = area / ( option%flow_dt * FMWH2O * ABS(option%gravity(3)) )
    endif

  end subroutine InlineSurfaceAccumulationJac

  !************************************************************************** !

  subroutine InlineSurfaceFlux(auxvar_up,auxvar_dn,area,dist,Res)
    ! 
    ! Compute the inline surface flux term to be used in Richards mode
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    ! 
    implicit none
    type(inlinesurface_auxvar_type) :: auxvar_up,auxvar_dn
    PetscReal :: area,dist(-1:3),Res(1)

#if 0 
    ! full diffusive wave approximation with no upwinding
    PetscReal :: rho,dz,Cf,gradZ,gradH,epsilon,v,havg
    epsilon = 1.0d-4
    rho     = 0.5d0*(auxvar_up%density             + auxvar_dn%density) 
    dz      = 0.5d0*(auxvar_up%half_cell_height    + auxvar_dn%half_cell_height) 
    Cf      = 0.5d0*(auxvar_up%Mannings_coeff      + auxvar_dn%Mannings_coeff)
    havg    = 0.5d0*(auxvar_up%surface_water_depth + auxvar_dn%surface_water_depth)
    havg    = MAX(havg,0.0d0)
    gradH   = (auxvar_up%surface_water_depth - auxvar_dn%surface_water_depth &
         + dist(3)*dist(0))/dist(0)
    v       = -havg**(2./3.)/(Cf*SQRT(SQRT(gradH*gradH+epsilon*epsilon)))*gradH
    Res(1)  = rho*havg*v*(area/(2.0d0*dz))

#else
    ! hybrid approximation with Scott's upwinding
    PetscReal :: Cf,dphi,dz,epsilon,have,rho,rhohu,slope
    epsilon = 1.0d-4
    slope   = ABS(dist(3)) + epsilon
    dphi    = auxvar_up%surface_water_depth - auxvar_dn%surface_water_depth &
         - dist(3)*dist(0)
    if (dphi > 0) then
      have = auxvar_up%surface_water_depth - MAX(+dist(3)*dist(0),0.0d0)
    else
      have = auxvar_dn%surface_water_depth - MAX(-dist(3)*dist(0),0.0d0)
    endif
    rho = 0.5d0*(auxvar_up%density         +auxvar_dn%density) 
    dz  = 0.5d0*(auxvar_up%half_cell_height+auxvar_dn%half_cell_height) 
    Cf  = 0.5d0*(auxvar_up%Mannings_coeff  +auxvar_dn%Mannings_coeff)
    if (have > 0.0d0) then
      rhohu  = have**(5./3.)*dphi*rho/(Cf*dist(0)*SQRT(slope))
      Res(1) = 0.5d0*rhohu*area/dz
    else
      Res(1) = 0.0d0
    endif
#endif

  end subroutine InlineSurfaceFlux

  !************************************************************************** !

  subroutine InlineSurfaceBCFlux(ibndtype,auxvar_up,auxvar_dn,area,dist,Res)
    ! 
    ! Compute the inline surface boundary flux term to be used in Richards mode
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    !
    use PFLOTRAN_Constants_module
    implicit none
    PetscInt  :: ibndtype(:)
    type(inlinesurface_auxvar_type) :: auxvar_up,auxvar_dn
    PetscReal :: area,dist(-1:3),Res(1)

    select case(ibndtype(TH_PRESSURE_DOF))
      case(SURFACE_DIRICHLET)
        call InlineSurfaceFlux(auxvar_up,auxvar_dn,area,dist,Res)
      case(SURFACE_SPILLOVER)
        call InlineSurfaceFlux(auxvar_up,auxvar_dn,area,dist,Res)
        if (auxvar_dn%surface_water_depth < auxvar_up%surface_water_depth) then
          Res(TH_PRESSURE_DOF) = 0.0d0
        endif
      case(SURFACE_ZERO_GRADHEIGHT)
        auxvar_dn%surface_water_depth = auxvar_up%surface_water_depth
        call InlineSurfaceFlux(auxvar_up,auxvar_dn,area,dist,Res)
    end select

  end subroutine InlineSurfaceBCFlux

  !************************************************************************** !

  subroutine InlineSurfaceFluxJac(auxvar_up,auxvar_dn,area,dist,option,Jup,Jdn)
    ! 
    ! Compute the Jacobian of the inline surface flux term to be used
    ! in Richards mode
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    !
    use Option_module
    use PFLOTRAN_Constants_module, only : FMWH2O
    implicit none
    type(inlinesurface_auxvar_type) :: auxvar_up,auxvar_dn
    type(option_type)               :: option
    PetscReal :: area,dist(-1:3),Jup(1,1),Jdn(1,1)
    PetscReal :: Cf,dphi,dz,epsilon,have,rho,slope,const,deriv

    epsilon  = 1.0d-4
    slope    = ABS(dist(3)) + epsilon
    dphi     = auxvar_up%surface_water_depth - auxvar_dn%surface_water_depth &
         - dist(3)*dist(0)
    rho      = 0.5d0*(auxvar_up%density         +auxvar_dn%density) 
    dz       = 0.5d0*(auxvar_up%half_cell_height+auxvar_dn%half_cell_height) 
    Cf       = 0.5d0*(auxvar_up%Mannings_coeff  +auxvar_dn%Mannings_coeff)
    const    = 1.0d0/(SQRT(slope)*dist(0)*Cf)/ABS(option%gravity(3))/ &
         FMWH2O*area/(2.0d0*dz)
    Jup(1,1) = 0.0d0
    Jdn(1,1) = 0.0d0
    if (dphi > 0) then
      have = auxvar_up%surface_water_depth - MAX(+dist(3)*dist(0),0.0d0)
      if (have > 0.0d0) then
        deriv = have**(2./3.) * ( 5./3. *dphi + have )
        Jup(1,1) = deriv*const
        deriv = -have**(5./3.)
        Jdn(1,1) = deriv*const
      endif
    else
      have = auxvar_dn%surface_water_depth - MAX(-dist(3)*dist(0),0.0d0)
      if (have > 0.0d0) then
        deriv = have**(5./3.)
        Jup(1,1) = deriv*const
        deriv = have**(2./3.) * ( 5./3. *dphi - have )
        Jdn(1,1) = deriv*const
      endif
    endif

  end subroutine InlineSurfaceFluxJac

  !************************************************************************** !

  subroutine InlineSurfaceBCFluxJac(ibndtype,auxvar_up,auxvar_dn,area,dist,option,Jdn)
    ! 
    ! Compute the Jacobian of the inline surface boundary flux term to
    ! be used in Richards mode
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    !
    use Option_module
    use PFLOTRAN_Constants_module
    implicit none
    PetscInt :: ibndtype(:)
    type(inlinesurface_auxvar_type) :: auxvar_up,auxvar_dn
    type(option_type)               :: option
    PetscReal :: area,dist(-1:3),Jup(1,1),Jdn(1,1)

    select case(ibndtype(TH_PRESSURE_DOF))
      case(SURFACE_DIRICHLET)
        call InlineSurfaceFluxJac(auxvar_up,auxvar_dn,area,dist,option,Jup,Jdn)
      case(SURFACE_SPILLOVER)
        call InlineSurfaceFluxJac(auxvar_up,auxvar_dn,area,dist,option,Jup,Jdn)
        if (auxvar_dn%surface_water_depth < auxvar_up%surface_water_depth) then
          Jdn = 0.d0
        endif
      case(SURFACE_ZERO_GRADHEIGHT)
        auxvar_dn%surface_water_depth = auxvar_up%surface_water_depth
        call InlineSurfaceFluxJac(auxvar_up,auxvar_dn,area,dist,option,Jup,Jdn)
    end select

  end subroutine InlineSurfaceBCFluxJac

end module InlineSurface_module
