module Transport_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Matrix_Block_Aux_module  

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  
  public :: TDispersion, &
            TDispersionBC, &
            TFlux, &
            TFluxDerivative, &
            TFluxCoef, &
            TSrcSinkCoef, &
            TFlux_CD, &
            TFluxDerivative_CD, &
            TFluxCoef_CD, &
            TFluxTVD
  
  ! this interface is required for the pointer to procedure employed
  ! for flux limiters below
  interface
    function TFluxLimiterDummy(d)
      PetscReal :: d
      PetscReal :: TFluxLimiterDummy
    end function TFluxLimiterDummy
  end interface
  
  public :: TFluxLimiterDummy, &
            TFluxLimiter, &
            TFluxLimitUpwind, &
            TFluxLimitMinmod, &
            TFluxLimitMC, &
            TFluxLimitSuperBee, &
            TFluxLimitVanLeer
  
  PetscInt, parameter, public :: TVD_LIMITER_UPWIND = 1
  PetscInt, parameter, public :: TVD_LIMITER_MC = 2
  PetscInt, parameter, public :: TVD_LIMITER_MINMOD = 3
  PetscInt, parameter, public :: TVD_LIMITER_SUPERBEE = 4
  PetscInt, parameter, public :: TVD_LIMITER_VAN_LEER = 5
              
contains

! ************************************************************************** !

subroutine TDispersion(global_auxvar_up,material_auxvar_up, &
                      cell_centered_velocity_up,dispersivity_up, &
                      global_auxvar_dn,material_auxvar_dn, &
                      cell_centered_velocity_dn,dispersivity_dn,dist, &
                      rt_parameter,option,qdarcy, &
                      harmonic_tran_coefs_over_dist)
  ! 
  ! Computes a single coefficient representing:
  !   (saturation * porosity * 
  !   (mechanical_dispersion + tortuosity * molecular_diffusion) / 
  !    distance [between cell centers] through a harmonic average
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/31/16
  ! 

  use Option_module
  use Connection_module

  implicit none
  
  type(option_type) :: option
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn 
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  PetscReal :: dispersivity_up(3), dispersivity_dn(3)
  PetscReal :: cell_centered_velocity_up(3,2), &
               cell_centered_velocity_dn(3,2)
  PetscReal :: dist(-1:3)
  PetscReal :: qdarcy(*)
  type(reactive_transport_param_type) :: rt_parameter
  PetscReal :: harmonic_tran_coefs_over_dist(rt_parameter%naqcomp, &
                                             rt_parameter%nphase)
  
  PetscInt :: iphase, nphase
  PetscReal :: dist_up, dist_dn
  PetscReal :: sat_up, sat_dn
  PetscReal :: velocity_dn(3), velocity_up(3)
  PetscReal :: distance_gravity, upwind_weight ! both are dummy variables
  PetscReal :: q
  PetscReal :: Dxx_up, Dyy_up, Dzz_up
  PetscReal :: Dxx_dn, Dyy_dn, Dzz_dn
  PetscInt, parameter :: LONGITUDINAL = 1
  PetscInt, parameter :: TRANSVERSE_HORIZONTAL = 2
  PetscInt, parameter :: TRANSVERSE_VERTICAL = 3
  PetscReal :: mechanical_dispersion_up
  PetscReal :: mechanical_dispersion_dn
  PetscReal :: molecular_diffusion_up(rt_parameter%naqcomp)
  PetscReal :: molecular_diffusion_dn(rt_parameter%naqcomp)
  PetscReal :: hydrodynamic_dispersion_up(rt_parameter%naqcomp)
  PetscReal :: hydrodynamic_dispersion_dn(rt_parameter%naqcomp)
  PetscReal :: t_ref_inv

  nphase = rt_parameter%nphase
  
  harmonic_tran_coefs_over_dist(:,:) = 0.d0    

  call ConnectionCalculateDistances(dist,option%gravity,dist_up, &
                                    dist_dn,distance_gravity, &
                                    upwind_weight)
  do iphase = 1, nphase
    sat_up = global_auxvar_up%sat(iphase)
    sat_dn = global_auxvar_dn%sat(iphase)
    ! skip phase if it does not exist on either side of the connection
    if (sat_up < eps .or. sat_dn < eps) cycle 
    if (rt_parameter%temperature_dependent_diffusion) then
      select case(iphase)
        case(LIQUID_PHASE)
          t_ref_inv = 1.d0/298.15d0 ! 1. / (25. + 273.15)
          molecular_diffusion_up(:) = &
            rt_parameter%diffusion_coefficient(:,iphase) * &
            exp(rt_parameter%diffusion_activation_energy(:) / &
                IDEAL_GAS_CONSTANT* &
                (t_ref_inv - 1.d0/(global_auxvar_up%temp + 273.15d0)))
          molecular_diffusion_dn(:) = &
            rt_parameter%diffusion_coefficient(:,iphase) * &
            exp(rt_parameter%diffusion_activation_energy(:) / &
                IDEAL_GAS_CONSTANT* &
                (t_ref_inv - 1.d0/(global_auxvar_dn%temp + 273.15d0)))
        case(GAS_PHASE)
          ! if gas phase exists, gas pressure %pres(GAS_PHASE) should be total
          ! pressure
          molecular_diffusion_up(:) = &
            rt_parameter%diffusion_coefficient(:,iphase) * &
            ((((global_auxvar_up%temp+273.15d0)/273.15d0)**1.8d0) * &
             (101325.d0/global_auxvar_up%pres(GAS_PHASE)))
          molecular_diffusion_dn(:) = &
            rt_parameter%diffusion_coefficient(:,iphase) * &
            ((((global_auxvar_dn%temp+273.15d0)/273.15d0)**1.8d0) * &
             (101325.d0/global_auxvar_dn%pres(GAS_PHASE)))
      end select
    else
      molecular_diffusion_up(:) = &
            rt_parameter%diffusion_coefficient(:,iphase)
      molecular_diffusion_dn(:) = &
            rt_parameter%diffusion_coefficient(:,iphase)
    endif
    q = qdarcy(iphase)
    if (rt_parameter%calculate_transverse_dispersion) then
      velocity_dn = q*dist(1:3) + (1.d0-dist(1:3))* &
                    cell_centered_velocity_dn(:,iphase)
      velocity_up = q*dist(1:3) + (1.d0-dist(1:3))* &
                    cell_centered_velocity_up(:,iphase)
      Dxx_up = dispersivity_up(LONGITUDINAL)* &
                 dabs(velocity_up(X_DIRECTION)) + &
               dispersivity_up(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_up(Y_DIRECTION)) + &
               dispersivity_up(TRANSVERSE_VERTICAL)* &
                 dabs(velocity_up(Z_DIRECTION))
      Dxx_dn = dispersivity_dn(LONGITUDINAL)* &
                 dabs(velocity_dn(X_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(Y_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_VERTICAL)* &
                 dabs(velocity_dn(Z_DIRECTION))
      Dyy_up = dispersivity_up(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_up(X_DIRECTION)) + &
               dispersivity_up(LONGITUDINAL)* &
                 dabs(velocity_up(Y_DIRECTION)) + &
               dispersivity_up(TRANSVERSE_VERTICAL)* &
                 dabs(velocity_up(Z_DIRECTION))
      Dyy_dn = dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(X_DIRECTION)) + &
               dispersivity_dn(LONGITUDINAL)* &
                 dabs(velocity_dn(Y_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_VERTICAL)* &
                 dabs(velocity_dn(Z_DIRECTION))
      Dzz_up = dispersivity_up(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_up(X_DIRECTION)) + &
               dispersivity_up(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_up(Y_DIRECTION)) + &
               dispersivity_up(LONGITUDINAL)* &
                 dabs(velocity_up(Z_DIRECTION))
      Dzz_dn = dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(X_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(Y_DIRECTION)) + &
               dispersivity_dn(LONGITUDINAL)* &
                 dabs(velocity_dn(Z_DIRECTION))
      ! dot product on unit direction vector
      mechanical_dispersion_up = &
        max(dist(1)*Dxx_up+dist(2)*Dyy_up+dist(3)*Dzz_up,1.d-40)
      mechanical_dispersion_dn = &
        max(dist(1)*Dxx_dn+dist(2)*Dyy_dn+dist(3)*Dzz_dn,1.d-40)
    else
      mechanical_dispersion_up = dispersivity_up(LONGITUDINAL)*dabs(q) 
      mechanical_dispersion_dn = dispersivity_dn(LONGITUDINAL)*dabs(q) 
    endif
    ! hydrodynamic dispersion = mechanical disperson + &
    !   saturation * porosity * tortuosity * molecular diffusion
    hydrodynamic_dispersion_up(:) = &
      max(mechanical_dispersion_up + &
          sat_up * material_auxvar_up%porosity * &
          material_auxvar_up%tortuosity * molecular_diffusion_up(:), &
          1.d-40)
    hydrodynamic_dispersion_dn(:) = &
      max(mechanical_dispersion_dn + &
          sat_dn * material_auxvar_dn%porosity * &
          material_auxvar_dn%tortuosity * molecular_diffusion_dn(:), &
          1.d-40)
    ! harmonic average of hydrodynamic dispersion divided by distance
    harmonic_tran_coefs_over_dist(:,iphase) = &
      (hydrodynamic_dispersion_up(:)*hydrodynamic_dispersion_dn(:))/ &
      (hydrodynamic_dispersion_up(:)*dist_dn + &
       hydrodynamic_dispersion_dn(:)*dist_up)
  enddo

end subroutine TDispersion

! ************************************************************************** !

subroutine TDispersionBC(ibndtype, &
                          global_auxvar_up, &
                          global_auxvar_dn,material_auxvar_dn, &
                          cell_centered_velocity_dn,dispersivity_dn, &
                          dist_dn, &
                          rt_parameter,option,qdarcy, &
                          tran_coefs_over_dist)
  ! 
  ! Computes a single coefficient representing:
  !   (saturation * porosity * 
  !   (mechanical_dispersion + tortuosity * molecular_diffusion) / 
  !    distance [between cell centers] 
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/31/16
  ! 

  use Option_module
  use Connection_module

  implicit none
  
  PetscInt :: ibndtype
  type(option_type) :: option
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn 
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: dispersivity_dn(3)
  PetscReal :: cell_centered_velocity_dn(3,2)
  PetscReal :: dist_dn(-1:3)
  PetscReal :: qdarcy(*)
  type(reactive_transport_param_type) :: rt_parameter
  PetscReal :: tran_coefs_over_dist(rt_parameter%naqcomp, &
                                             rt_parameter%nphase)
  
  PetscInt :: iphase, nphase
  PetscReal :: spD
  PetscReal :: sat_up, sat_dn
  PetscReal :: velocity_dn(3)
  PetscReal :: distance_gravity, upwind_weight ! both are dummy variables
  PetscReal :: q
  PetscReal :: Dxx, Dyy, Dzz
  PetscInt, parameter :: LONGITUDINAL = 1
  PetscInt, parameter :: TRANSVERSE_HORIZONTAL = 2
  PetscInt, parameter :: TRANSVERSE_VERTICAL = 3
  PetscReal :: mechanical_dispersion
  PetscReal :: molecular_diffusion(rt_parameter%naqcomp)
  PetscReal :: hydrodynamic_dispersion(rt_parameter%naqcomp)
  PetscReal :: t_ref_inv

  nphase = rt_parameter%nphase
  
  tran_coefs_over_dist(:,:) = 0.d0    

  do iphase = 1, nphase
    ! we use upwind saturation as that is the saturation at the boundary face
    sat_up = global_auxvar_up%sat(iphase)
    sat_dn = global_auxvar_dn%sat(iphase)
    if (sat_up < eps .or. sat_dn < eps) cycle 
    if (rt_parameter%temperature_dependent_diffusion) then
      select case(iphase)
        case(LIQUID_PHASE)
          t_ref_inv = 1.d0/298.15d0 ! 1. / (25. + 273.15)
          molecular_diffusion(:) = &
            rt_parameter%diffusion_coefficient(:,iphase) * &
            exp(rt_parameter%diffusion_activation_energy(:) / &
                IDEAL_GAS_CONSTANT* &
                (t_ref_inv - 1.d0/(global_auxvar_up%temp + 273.15d0)))
        case(GAS_PHASE)
          ! if gas phase exists, gas pressure %pres(GAS_PHASE) should be total
          ! pressure
          molecular_diffusion(:) = &
            rt_parameter%diffusion_coefficient(:,iphase) * &
            ((((global_auxvar_up%temp+273.15d0)/273.15d0)**1.8d0) * &
             (101325.d0/global_auxvar_up%pres(GAS_PHASE)))
      end select
    else
      molecular_diffusion(:) = &
            rt_parameter%diffusion_coefficient(:,iphase)
    endif
    q = qdarcy(iphase)
    if (rt_parameter%calculate_transverse_dispersion) then
      velocity_dn = q*dist_dn(1:3) + (1.d0-dist_dn(1:3))* &
                    cell_centered_velocity_dn(:,iphase)
      Dxx = dispersivity_dn(LONGITUDINAL)* &
                 dabs(velocity_dn(X_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(Y_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_VERTICAL)* &
                 dabs(velocity_dn(Z_DIRECTION))
      Dyy = dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(X_DIRECTION)) + &
               dispersivity_dn(LONGITUDINAL)* &
                 dabs(velocity_dn(Y_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_VERTICAL)* &
                 dabs(velocity_dn(Z_DIRECTION))
      Dzz = dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(X_DIRECTION)) + &
               dispersivity_dn(TRANSVERSE_HORIZONTAL)* &
                 dabs(velocity_dn(Y_DIRECTION)) + &
               dispersivity_dn(LONGITUDINAL)* &
                 dabs(velocity_dn(Z_DIRECTION))
      ! dot product on unit direction vector
      mechanical_dispersion = &
        max(dist_dn(1)*Dxx+dist_dn(2)*Dyy+dist_dn(3)*Dzz,1.d-40)
    else
      mechanical_dispersion = dispersivity_dn(LONGITUDINAL)*dabs(q) 
    endif

    select case(ibndtype)
      case(DIRICHLET_BC,DIRICHLET_ZERO_GRADIENT_BC)
        ! if outflow, skip
        if (ibndtype == DIRICHLET_ZERO_GRADIENT_BC .and. q < 0.d0) cycle
        ! hydrodynamic dispersion = mechanical disperson + &
        !   saturation * porosity * tortuosity * molecular diffusion
        hydrodynamic_dispersion(:) = &
          max(mechanical_dispersion + &
              sat_up * material_auxvar_dn%porosity * &
              material_auxvar_dn%tortuosity * molecular_diffusion(:), &
              1.d-40)
        ! hydrodynamic dispersion divided by distance
        ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
        tran_coefs_over_dist(:,iphase) =  &
          hydrodynamic_dispersion(:)/dist_dn(0)
      case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
    end select
  enddo

end subroutine TDispersionBC

! ************************************************************************** !

subroutine TFlux(rt_parameter, &
                 rt_auxvar_up,global_auxvar_up, & 
                 rt_auxvar_dn,global_auxvar_dn, & 
                 coef_up,coef_dn,option,Res)
  ! 
  ! Computes flux term in residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn 
  type(option_type) :: option
  PetscReal :: coef_up(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_dn(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: Res(rt_parameter%ncomp)
  
  PetscInt :: iphase
  PetscInt :: idof
  PetscInt :: ndof
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: icollcomp
  PetscInt :: icoll
  PetscInt :: iaqcomp
  PetscInt :: neg999

  iphase = 1
  ndof = rt_parameter%naqcomp
  
  Res = 0.d0
  
  ! units = (L water/sec)*(mol/L) = mol/s
  ! total = mol/L water
  Res(1:ndof) = coef_up(1:ndof,iphase)*rt_auxvar_up%total(1:ndof,iphase) + &
                coef_dn(1:ndof,iphase)*rt_auxvar_dn%total(1:ndof,iphase)

  if (rt_parameter%ncoll > 0) then
    neg999 = -999
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_colloid + icoll
      Res(idof) = &
       ! conc_mob = mol/L water
        !geh: neg999 set to catch use of colloids.
        coef_up(neg999,iphase)*rt_auxvar_up%colloid%conc_mob(icoll)+ &
        coef_dn(neg999,iphase)*rt_auxvar_dn%colloid%conc_mob(icoll)
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    neg999 = -999
    do icollcomp = 1, rt_parameter%ncollcomp
      iaqcomp = rt_parameter%coll_spec_to_pri_spec(icollcomp)
      ! total_eq_mob = mol/L water
      Res(iaqcomp) = Res(iaqcomp) + &
        coef_up(neg999,iphase)*rt_auxvar_up%colloid%total_eq_mob(icollcomp) + &
        coef_dn(neg999,iphase)*rt_auxvar_dn%colloid%total_eq_mob(icollcomp)
    enddo
  endif
  if (rt_parameter%ngas > 0) then
    iphase = 2
    Res(1:ndof) = Res(1:ndof) + &
                  coef_up(1:ndof,iphase)*rt_auxvar_up%total(1:ndof,iphase) + &
                  coef_dn(1:ndof,iphase)*rt_auxvar_dn%total(1:ndof,iphase)
  endif

end subroutine TFlux

! ************************************************************************** !

subroutine TFlux_CD(rt_parameter, &
                 rt_auxvar_up,global_auxvar_up, & 
                 rt_auxvar_dn,global_auxvar_dn, & 
                 coef_11,coef_12,coef_21,coef_22,option,Res_1,Res_2)
  ! 
  ! TFlux: Computes flux term in residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn 
  type(option_type) :: option
  PetscReal :: coef_11(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_12(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_21(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_22(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: Res_1(rt_parameter%ncomp)
  PetscReal :: Res_2(rt_parameter%ncomp)
  
  PetscInt :: iphase
  PetscInt :: idof
  PetscInt :: ndof
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: icollcomp
  PetscInt :: icoll
  PetscInt :: iaqcomp
  PetscInt :: neg999

  iphase = 1
  ndof = rt_parameter%naqcomp
  
  Res_1 = 0.d0
  Res_2 = 0.d0
  
  ! units = (L water/sec)*(mol/L) = mol/s
  ! total(:,1) = mol/L water
  Res_1(1:ndof) = coef_11(1:ndof,iphase)*rt_auxvar_up%total(1:ndof,iphase) + &
                  coef_12(1:ndof,iphase)*rt_auxvar_dn%total(1:ndof,iphase)
  Res_2(1:ndof) = coef_21(1:ndof,iphase)*rt_auxvar_up%total(1:ndof,iphase) + &
                  coef_22(1:ndof,iphase)*rt_auxvar_dn%total(1:ndof,iphase)

  if (rt_parameter%ncoll > 0) then
    neg999 = -999
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_colloid + icoll
       ! conc_mob = mol/L water
      Res_1(idof) = coef_11(neg999,iphase)*rt_auxvar_up%colloid%conc_mob(icoll)+ &
                    coef_12(neg999,iphase)*rt_auxvar_dn%colloid%conc_mob(icoll)
      Res_2(idof) = coef_21(neg999,iphase)*rt_auxvar_up%colloid%conc_mob(icoll)+ &
                    coef_22(neg999,iphase)*rt_auxvar_dn%colloid%conc_mob(icoll)
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    neg999 = -999
    do icollcomp = 1, rt_parameter%ncollcomp
      iaqcomp = rt_parameter%coll_spec_to_pri_spec(icollcomp)
      ! total_eq_mob = mol/L water
      Res_1(iaqcomp) = Res_1(iaqcomp) + &
        coef_11(neg999,iphase)*rt_auxvar_up%colloid%total_eq_mob(icollcomp) + &
        coef_12(neg999,iphase)*rt_auxvar_dn%colloid%total_eq_mob(icollcomp)
      Res_2(iaqcomp) = Res_2(iaqcomp) + &
        coef_21(neg999,iphase)*rt_auxvar_up%colloid%total_eq_mob(icollcomp) + &
        coef_22(neg999,iphase)*rt_auxvar_dn%colloid%total_eq_mob(icollcomp)
    enddo
  endif
  if (rt_parameter%ngas > 0) then
    iphase = 2
  ! total(:,2) = mol/L gas
    Res_1(1:ndof) = coef_11(1:ndof,iphase)*rt_auxvar_up%total(1:ndof,iphase) + &
                    coef_12(1:ndof,iphase)*rt_auxvar_dn%total(1:ndof,iphase)
    Res_2(1:ndof) = coef_21(1:ndof,iphase)*rt_auxvar_up%total(1:ndof,iphase) + &
                    coef_22(1:ndof,iphase)*rt_auxvar_dn%total(1:ndof,iphase)
  endif
  
end subroutine TFlux_CD

! ************************************************************************** !

subroutine TFluxDerivative(rt_parameter, &
                           rt_auxvar_up,global_auxvar_up, & 
                           rt_auxvar_dn,global_auxvar_dn, & 
                           coef_up,coef_dn,option,J_up,J_dn)
  ! 
  ! Computes derivatives of flux term in residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(option_type) :: option
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn 
  PetscReal :: coef_up(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_dn(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: J_up(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_dn(rt_parameter%ncomp,rt_parameter%ncomp)
  
  PetscInt :: iphase
  PetscInt :: icomp
  PetscInt :: icoll
  PetscInt :: idof
  PetscInt :: istart
  PetscInt :: iendaq
  PetscInt :: nphase
  PetscInt :: irow
  PetscInt :: neg999
 
  nphase = rt_parameter%nphase
  
  ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water)
  !       = kg water/sec
  istart = 1
  iendaq = rt_parameter%naqcomp
  J_up = 0.d0
  J_dn = 0.d0
  do iphase = 1, nphase
    if (associated(rt_auxvar_dn%aqueous%dtotal)) then
      do irow = istart, iendaq
        J_up(irow,istart:iendaq) = &
          J_up(irow,istart:iendaq) + &
          rt_auxvar_up%aqueous%dtotal(irow,:,iphase)* &
          coef_up(irow,iphase)
        J_dn(irow,istart:iendaq) = &
          J_dn(irow,istart:iendaq) + &
          rt_auxvar_dn%aqueous%dtotal(irow,:,iphase)* &
          coef_dn(irow,iphase)
      enddo
    else  
      do icomp = istart, iendaq
        J_up(icomp,icomp) = J_up(icomp,icomp) + coef_up(icomp,iphase)* &
                            global_auxvar_up%den_kg(iphase)*1.d-3
        J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn(icomp,iphase)* &
                            global_auxvar_dn%den_kg(iphase)*1.d-3
      enddo
    endif
  enddo

  iphase = 1
  if (rt_parameter%ncoll > 0) then
    neg999 = -999
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_colloid + icoll
      J_up(idof,idof) = coef_up(neg999,iphase) * &
        global_auxvar_up%den_kg(iphase)*1.d-3
      J_dn(idof,idof) = coef_dn(neg999,iphase) * &
        global_auxvar_dn%den_kg(iphase)*1.d-3
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    ! dRj_dCj - mobile
    ! istart & iend same as above
    do irow = istart, iendaq
      J_up(irow,istart:iendaq) = J_up(irow,istart:iendaq) + &
        rt_auxvar_up%colloid%dRj_dCj%dtotal(irow,:,iphase)* &
        coef_up(irow,iphase)
      J_dn(irow,istart:iendaq) = J_dn(irow,istart:iendaq) + &
        rt_auxvar_dn%colloid%dRj_dCj%dtotal(irow,:,iphase)* &
        coef_dn(irow,iphase)
    enddo
    ! need the below
    ! dRj_dSic
    ! dRic_dSic
    ! dRic_dCj
  endif

end subroutine TFluxDerivative

! ************************************************************************** !

subroutine TFluxDerivative_CD(rt_parameter, &
                           rt_auxvar_up,global_auxvar_up, & 
                           rt_auxvar_dn,global_auxvar_dn, & 
                           coef_11,coef_12,coef_21,coef_22,option, &
                           J_11,J_12,J_21,J_22)
  ! 
  ! TFluxDerivative: Computes derivatives of flux term in residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(option_type) :: option
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn 
  PetscReal :: coef_11(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_12(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_21(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: coef_22(rt_parameter%naqcomp,rt_parameter%nphase) 
  PetscReal :: J_11(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_12(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_21(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_22(rt_parameter%ncomp,rt_parameter%ncomp)
  
  PetscInt :: iphase
  PetscInt :: icoll
  PetscInt :: idof
  PetscInt :: istart
  PetscInt :: iendaq
  PetscInt :: nphase
  PetscInt :: irow
  PetscInt :: icomp
  PetscInt :: neg999
  
  nphase = rt_parameter%nphase
  
  ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  istart = 1
  iendaq = rt_parameter%naqcomp
  J_11 = 0.d0
  J_12 = 0.d0
  J_21 = 0.d0
  J_22 = 0.d0
  do iphase = 1, nphase
    if (associated(rt_auxvar_dn%aqueous%dtotal)) then
      do irow = istart, iendaq
        J_11(irow,istart:iendaq) = &
          J_11(irow,istart:iendaq) + &
          rt_auxvar_up%aqueous%dtotal(irow,:,iphase)* &
          coef_11(irow,iphase)
        J_12(irow,istart:iendaq) = &
          J_12(irow,istart:iendaq) + &
          rt_auxvar_dn%aqueous%dtotal(irow,:,iphase)* &
          coef_12(irow,iphase)
        J_21(irow,istart:iendaq) = &
          J_21(irow,istart:iendaq) + &
          rt_auxvar_up%aqueous%dtotal(irow,:,iphase)* &
          coef_21(irow,iphase)
        J_22(irow,istart:iendaq) = &
          J_22(irow,istart:iendaq) + &
          rt_auxvar_dn%aqueous%dtotal(irow,:,iphase)* &
          coef_22(irow,iphase)
      enddo
    else  
      do icomp = istart, iendaq
        J_11(icomp,icomp) = &
          J_11(icomp,icomp) + &
          coef_11(icomp,iphase)*global_auxvar_up%den_kg(iphase)*1.d-3
        J_12(icomp,icomp) = &
          J_12(icomp,icomp) + &
          coef_12(icomp,iphase)*global_auxvar_dn%den_kg(iphase)*1.d-3
        J_21(icomp,icomp) = &
          J_21(icomp,icomp) + &
          coef_21(icomp,iphase)*global_auxvar_up%den_kg(iphase)*1.d-3
        J_22(icomp,icomp) = &
          J_22(icomp,icomp) + &
          coef_22(icomp,iphase)*global_auxvar_dn%den_kg(iphase)*1.d-3
      enddo
    endif
  enddo

  iphase = 1
  if (rt_parameter%ncoll > 0) then
    neg999 = -999
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_colloid + icoll
      J_11(idof,idof) = coef_11(neg999,iphase)* &
                        global_auxvar_up%den_kg(iphase)*1.d-3
      J_12(idof,idof) = coef_12(neg999,iphase)* &
                        global_auxvar_dn%den_kg(iphase)*1.d-3
      J_21(idof,idof) = coef_21(neg999,iphase)* &
                        global_auxvar_up%den_kg(iphase)*1.d-3
      J_22(idof,idof) = coef_22(neg999,iphase)* &
                        global_auxvar_dn%den_kg(iphase)*1.d-3
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    neg999 = -999
    ! dRj_dCj - mobile
    ! istart & iend same as above
    do irow = istart, iendaq
      J_11(irow,istart:iendaq) = J_11(irow,istart:iendaq) + &
        rt_auxvar_up%colloid%dRj_dCj%dtotal(irow,:,iphase)* &
        coef_11(neg999,iphase)
      J_12(irow,istart:iendaq) = J_12(irow,istart:iendaq) + &
        rt_auxvar_dn%colloid%dRj_dCj%dtotal(irow,:,iphase)* &
        coef_12(neg999,iphase)
      J_21(irow,istart:iendaq) = J_21(irow,istart:iendaq) + &
        rt_auxvar_up%colloid%dRj_dCj%dtotal(irow,:,iphase)* &
        coef_21(neg999,iphase)
      J_22(irow,istart:iendaq) = J_22(irow,istart:iendaq) + &
        rt_auxvar_dn%colloid%dRj_dCj%dtotal(irow,:,iphase)* &
        coef_22(neg999,iphase)
      ! need the below
      ! dRj_dSic
      ! dRic_dSic
      ! dRic_dCj
    enddo
  endif

end subroutine TFluxDerivative_CD

! ************************************************************************** !

subroutine TFluxCoef(rt_parameter,option,area,velocity, &
                     tran_coefs_over_dist, &
                     fraction_upwind,T_up,T_dn)
  ! 
  ! Computes flux coefficients for transport matrix
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/10
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter  
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: velocity(*)
  ! this is the harmonic mean of saturation * porosity * (mechanical 
  !   dispersion + tortuosity * molecular_diffusion) / distance
  PetscReal :: tran_coefs_over_dist(rt_parameter%naqcomp, &
                                    rt_parameter%nphase)
  PetscReal :: fraction_upwind
  PetscReal :: T_up(rt_parameter%naqcomp,rt_parameter%nphase)
  PetscReal :: T_dn(rt_parameter%naqcomp,rt_parameter%nphase)

  PetscInt :: iphase, nphase
  PetscReal :: coef_up(rt_parameter%naqcomp)
  PetscReal :: coef_dn(rt_parameter%naqcomp)
  PetscReal :: q
  
  nphase = rt_parameter%nphase
  
  do iphase = 1, nphase
    q = velocity(iphase)

    if (option%use_upwinding) then
      ! upstream weighting
      ! units = (m^3 water/m^2 bulk/sec)
      if (q > 0.d0) then
        coef_up(:) =  tran_coefs_over_dist(:,iphase)+q
        coef_dn(:) = -tran_coefs_over_dist(:,iphase)
      else
        coef_up(:) =  tran_coefs_over_dist(:,iphase)
        coef_dn(:) = -tran_coefs_over_dist(:,iphase)+q
      endif
    else
      ! central difference, currently assuming uniform grid spacing
      ! units = (m^3 water/m^2 bulk/sec)
      ! 
      coef_up(:) =  tran_coefs_over_dist(:,iphase) + &
                    (1.d0-fraction_upwind)*q
      coef_dn(:) = -tran_coefs_over_dist(:,iphase) + &
                    fraction_upwind*q
    endif  
  
    ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
    !       = L water/sec
    T_up(:,iphase) = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
    T_dn(:,iphase) = coef_dn*area*1000.d0
  enddo
    
end subroutine TFluxCoef

! ************************************************************************** !

subroutine TFluxCoef_CD(rt_parameter,option,area,velocity, &
                        tran_coefs_over_dist, &
                        fraction_upwind,T_11,T_12,T_21,T_22)
  ! 
  ! Computes flux coefficients for transport matrix
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/10
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter  
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: velocity(*)
  ! this is the harmonic mean of saturation * porosity * (mechanical 
  !   dispersion + tortuosity * molecular_diffusion) / distance
  PetscReal :: tran_coefs_over_dist(rt_parameter%naqcomp, &
                                    rt_parameter%nphase)
  PetscReal :: fraction_upwind
  PetscReal :: T_11(rt_parameter%naqcomp,rt_parameter%nphase)
  PetscReal :: T_12(rt_parameter%naqcomp,rt_parameter%nphase)
  PetscReal :: T_21(rt_parameter%naqcomp,rt_parameter%nphase)
  PetscReal :: T_22(rt_parameter%naqcomp,rt_parameter%nphase)

  PetscInt :: iphase
  PetscReal :: coef_up(rt_parameter%naqcomp), coef_dn(rt_parameter%naqcomp)
  PetscReal :: tempreal
  PetscReal :: advection_upwind(rt_parameter%nphase)
  PetscReal :: advection_downwind(rt_parameter%nphase)
  PetscReal :: q
  
  ! T_11 = diagonal term for upwind cell (row)
  ! T_12 = off diagonal term for upwind cell (row)
  ! T_21 = off diagonal term for downwind cell (row)
  ! T_22 = diagonal term for downwind cell (row)

  ! Advection
  if (option%use_upwinding) then
    ! upstream weighting
    ! units = (m^3 water/m^2 bulk/sec)
    do iphase = 1, rt_parameter%nphase
      if (velocity(iphase) > 0.d0) then
        advection_upwind(iphase) = velocity(iphase)
        advection_downwind(iphase) = 0.d0
      else
        advection_upwind(iphase) = 0.d0
        advection_downwind(iphase) = velocity(iphase)
      endif
    enddo
  else
    ! central difference
    do iphase = 1, rt_parameter%nphase
      advection_upwind(iphase) = (1.d0-fraction_upwind)*velocity(iphase)
      advection_downwind(iphase) = fraction_upwind*velocity(iphase)
    enddo
  endif
    
  tempreal = area*1000.d0
  do iphase = 1, rt_parameter%nphase
    T_11(:,iphase) = (tran_coefs_over_dist(:,iphase) + &
                      advection_upwind(iphase))*tempreal
    T_12(:,iphase) = (-tran_coefs_over_dist(:,iphase) + &
                      advection_downwind(iphase))*tempreal
!    T_21(:,iphase) = -(tran_coefs_over_dist(:,iphase) + &
!                       advection_upwind(iphase))*tempreal
!    T_22(:,iphase) = (tran_coefs_over_dist(:,iphase) - &
!                       advection_downwind(iphase))*tempreal
    T_21(:,iphase) = -T_11(:,iphase)
    T_22(:,iphase) = -T_12(:,iphase)
  enddo

end subroutine TFluxCoef_CD

! ************************************************************************** !

subroutine TSrcSinkCoef(rt_parameter,qsrc,tran_src_sink_type,T_in,T_out)
  ! 
  ! Computes src/sink coefficients for transport matrix
  ! Here qsrc [m^3/sec] provided by flow.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module

  implicit none

  type(reactive_transport_param_type) :: rt_parameter
  PetscReal :: qsrc(2)
  PetscInt :: tran_src_sink_type
  PetscReal :: T_in(2) ! coefficient that scales concentration at cell
  PetscReal :: T_out(2) ! concentration that scales external concentration
  
  PetscInt :: iphase
      
  T_in = 0.d0 
  T_out = 0.d0
     
  select case(tran_src_sink_type)
    case(EQUILIBRIUM_SS)
      ! units should be mol/sec
      ! 1.d-3 is a relatively large rate designed to equilibrate 
      ! the aqueous concentration with the concentrations specified at
      ! the src/sink
      T_in = 1.d-3 ! units L water/sec
      T_out = -1.d0*T_in
    case(MASS_RATE_SS)
      ! in this case, rt_auxvar_bc%total actually holds the mass rate
      T_in = 0.d0
      T_out = -1.d0
    case default
      ! qsrc always in m^3/sec
      do iphase = 1, rt_parameter%nphase
        if (qsrc(iphase) > 0.d0) then ! injection
          T_in(iphase) = 0.d0
          T_out(iphase) = -1.d0*qsrc(iphase)*1000.d0 ! m^3/sec * 1000 L/m^3 -> L/s
        else
          T_out(iphase) = 0.d0
          T_in(iphase) = -1.d0*qsrc(iphase)*1000.d0 ! m^3/sec * 1000 L/m^3 -> L/s
        endif
      enddo
  end select

  ! Units of Tin & Tout should be L/s.  When multiplied by Total (M) you get
  ! moles/sec, the units of the residual.  To get the units of the Jacobian
  ! kg/sec, one must either scale by dtotal or den/1000. (kg/L).

end subroutine TSrcSinkCoef

! ************************************************************************** !

subroutine TFluxTVD(rt_parameter,velocity,area,dist, &
                    total_up2,rt_auxvar_up, & 
                    rt_auxvar_dn,total_dn2, & 
                    TFluxLimitPtr, &
                    option,flux)
  ! 
  ! Computes TVD flux term
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  PetscReal :: velocity(:), area
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  PetscReal, pointer :: total_up2(:,:), total_dn2(:,:)
  type(option_type) :: option
  PetscReal :: flux(rt_parameter%ncomp)
  procedure (TFluxLimiterDummy), pointer :: TFluxLimitPtr
  PetscReal :: dist(-1:3)    ! list of distance vectors, size(-1:3,num_connections) where
                            !   -1 = fraction upwind
                            !   0 = magnitude of distance 
                            !   1-3 = components of unit vector 
                            
  PetscInt :: iphase
  PetscInt :: idof, ndof
  PetscReal :: dc, theta, correction, nu, velocity_area
  
  ndof = rt_parameter%naqcomp

  flux = 0.d0
  
  ! flux should be in mol/sec
  
  do iphase = 1, rt_parameter%nphase
    nu = velocity(iphase)*option%tran_dt/dist(0)
    ! L/sec = m/sec * m^2 * 1000 [L/m^3]
    velocity_area = velocity(iphase)*area*1000.d0
    if (velocity_area >= 0.d0) then
      ! mol/sec = L/sec * mol/L
      flux = velocity_area*rt_auxvar_up%total(1:rt_parameter%naqcomp,iphase)
      if (associated(total_up2)) then
        do idof = 1, ndof
          dc = rt_auxvar_dn%total(idof,iphase) - &
               rt_auxvar_up%total(idof,iphase)
          if (dabs(dc) < 1.d-20) then
            theta = 1.d0
          else
            theta = (rt_auxvar_up%total(idof,iphase) - &
                    total_up2(idof,iphase)) / &
                    dc
          endif
          ! mol/sec = L/sec * mol/L
          correction = 0.5d0*velocity_area*(1.d0-nu)* &
                       TFluxLimitPtr(theta)* &
                       dc
          flux(idof) = flux(idof) + correction
        enddo
      endif
    else
      flux = velocity_area*rt_auxvar_dn%total(1:rt_parameter%naqcomp,iphase)
      if (associated(total_dn2)) then
        do idof = 1, ndof
          dc = rt_auxvar_dn%total(idof,iphase) - &
               rt_auxvar_up%total(idof,iphase)
          if (dabs(dc) < 1.d-20) then
            theta = 1.d0
          else
            theta = (total_dn2(idof,iphase) - &
                     rt_auxvar_dn%total(idof,iphase)) / &
                    dc
          endif
          correction = 0.5d0*velocity_area*(1.d0+nu)* &
                       TFluxLimitPtr(theta)* &
                       dc
          flux(idof) = flux(idof) + correction
        enddo
      endif
    endif
  enddo

end subroutine TFluxTVD

! ************************************************************************** !

function TFluxLimiter(theta)
  ! 
  ! Applies flux limiter
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  implicit none
  
  PetscReal :: theta
  
  PetscReal :: TFluxLimiter
  
  ! Linear
  !---------
  ! upwind
  TFluxLimiter = 0.d0
  ! Lax-Wendroff
  !TFluxLimiter = 1.d0
  ! Beam-Warming
  !TFluxLimiter = theta
  ! Fromm
  !TFluxLimiter = 0.5d0*(1.d0+theta)

  ! Higher-order
  !---------
  ! minmod
  !TFluxLimiter = max(0.d0,min(1.d0,theta))
  ! superbee
  !TFluxLimiter = max(0.d0,min(1.d0,2.d0*theta),min(2.d0,theta))
  ! MC
  !TFluxLimiter = max(0.d0,min((1.d0+theta)/2.d0,2.d0,2.d0*theta))
  ! van Leer
  !TFluxLimiter = (theta+dabs(theta))/(1.d0+dabs(theta)

end function TFluxLimiter

! ************************************************************************** !

function TFluxLimitUpwind(theta)
  ! 
  ! Applies an upwind flux limiter
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  implicit none
  
  PetscReal :: theta
  
  PetscReal :: TFluxLimitUpwind
  
  ! upwind
  TFluxLimitUpwind = 0.d0

end function TFluxLimitUpwind

! ************************************************************************** !

function TFluxLimitMinmod(theta)
  ! 
  ! Applies a minmod flux limiter
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  implicit none
  
  PetscReal :: theta
  
  PetscReal :: TFluxLimitMinmod
  
  ! minmod
  TFluxLimitMinmod = max(0.d0,min(1.d0,theta))

end function TFluxLimitMinmod

! ************************************************************************** !

function TFluxLimitMC(theta)
  ! 
  ! Applies an MC flux limiter
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  implicit none
  
  PetscReal :: theta
  
  PetscReal :: TFluxLimitMC
  
 ! MC
  TFluxLimitMC = max(0.d0,min((1.d0+theta)/2.d0,2.d0,2.d0*theta))

end function TFluxLimitMC

! ************************************************************************** !

function TFluxLimitSuperBee(theta)
  ! 
  ! Applies an superbee flux limiter
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  implicit none
  
  PetscReal :: theta
  
  PetscReal :: TFluxLimitSuperBee
  
  ! superbee
  TFluxLimitSuperBee =  max(0.d0,min(1.d0,2.d0*theta),min(2.d0,theta))

end function TFluxLimitSuperBee

! ************************************************************************** !

function TFluxLimitVanLeer(theta)
  ! 
  ! Applies an van Leer flux limiter
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  implicit none
  
  PetscReal :: theta
  
  PetscReal :: TFluxLimitVanLeer
  
  ! superbee
  TFluxLimitVanLeer = (theta+dabs(theta))/(1.d0+dabs(theta))

end function TFluxLimitVanLeer

! ************************************************************************** !

end module Transport_module
