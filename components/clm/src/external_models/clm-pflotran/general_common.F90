module General_Common_module

  use General_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

#define CONVECTION
#define LIQUID_DARCY_FLUX
#define GAS_DARCY_FLUX
#define DIFFUSION
#define LIQUID_DIFFUSION
#define GAS_DIFFUSION
#define CONDUCTION

#define WATER_SRCSINK
#define AIR_SRCSINK
#define ENERGY_SRCSINK
  
!#define DEBUG_GENERAL_FILEOUTPUT
!#define DEBUG_FLUXES  

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

#ifdef DEBUG_GENERAL_FILEOUTPUT
  PetscInt, parameter :: debug_unit = 87
  PetscInt, parameter :: debug_info_unit = 86
  character(len=MAXWORDLENGTH) :: debug_filename
  PetscInt :: debug_flag = 0
  PetscInt :: debug_iteration_count
  PetscInt :: debug_timestep_cut_count
  PetscInt :: debug_timestep_count
#endif

  public :: GeneralAccumulation, &
            GeneralFlux, &
            GeneralBCFlux, &
            GeneralSrcSink, &
            GeneralAccumDerivative, &
            GeneralFluxDerivative, &
            GeneralBCFluxDerivative, &
            GeneralSrcSinkDerivative

contains

! ************************************************************************** !

subroutine GeneralAccumulation(gen_auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res,Jac, &
                               analytical_derivatives,debug_cell)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none

  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_cell
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: porosity
  PetscReal :: volume_over_dt
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use gen_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = gen_auxvar%effective_porosity
  
  ! accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] * 
    !                           xmol[kmol comp/kmol phase]
    do icomp = 1, option%nflowspec
#ifdef DEBUG_GENERAL
      ! for debug version, aux var entries are initialized to NaNs.  even if
      ! saturation is zero, density may be a NaN.  So the conditional prevents
      ! this calculation.  For non-debug, aux var entries are initialized to
      ! 0.d0
      if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
      Res(icomp) = Res(icomp) + gen_auxvar%sat(iphase) * &
                                gen_auxvar%den(iphase) * &
                                gen_auxvar%xmol(icomp,iphase)
#ifdef DEBUG_GENERAL
      endif
#endif
    enddo
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * volume_over_dt

  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
#ifdef DEBUG_GENERAL
    ! for debug version, aux var entries are initialized to NaNs.  even if
    ! saturation is zero, density may be a NaN.  So the conditional prevents
    ! this calculation.  For non-debug, aux var entries are initialized to
    ! 0.d0
    if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
    Res(energy_id) = Res(energy_id) + gen_auxvar%sat(iphase) * &
                                      gen_auxvar%den(iphase) * &
                                      gen_auxvar%U(iphase)
#ifdef DEBUG_GENERAL
    endif
#endif
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * gen_auxvar%temp) * volume_over_dt
  
  if (analytical_derivatives) then
    Jac = 0.d0
    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        ! satl = 1
        ! ----------
        ! Water Equation
        ! por * satl * denl * Xwl
        ! ---
        ! w/respect to liquid pressure
        ! dpor_dp * denl * Xwl + 
        ! por * ddenl_dpl * Xwl
        Jac(1,1) = &
          gen_auxvar%d%por_p * gen_auxvar%den(1) * gen_auxvar%xmol(1,1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%xmol(1,1)
        ! w/respect to air mole fraction
        ! liquid phase density is indepenent of air mole fraction
        ! por * denl * dXwl_dXal
        ! Xwl = 1. - Xal
        ! dXwl_dXal = -1.
        Jac(1,2) = porosity * gen_auxvar%den(1) * (-1.d0)
        ! w/repect to temperature
        ! por * ddenl_dT * Xwl
        Jac(1,3) = porosity * gen_auxvar%d%denl_T * gen_auxvar%xmol(1,1)
        ! ----------
        ! Air Equation
        ! por * satl * denl * Xal
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Xal + 
        ! por * ddenl_dpl * Xal
        Jac(2,1) = &
          gen_auxvar%d%por_p * gen_auxvar%den(1) * gen_auxvar%xmol(2,1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%xmol(2,1)
        ! w/respect to air mole fraction
        Jac(2,2) = porosity * gen_auxvar%den(1)
        ! w/repect to temperature
        ! por * ddenl_dT * Xwl
        Jac(2,3) = porosity * gen_auxvar%d%denl_T * gen_auxvar%xmol(2,1)
        ! ----------
        ! Energy Equation
        ! por * satl * denl * Ul + (1-por) * denr * Cp * T
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Ul + 
        ! por * ddenl_dpl * Ul + 
        ! por * denl * dUl_dpl + 
        ! -dpor_dpl * dens * Cp * T
        Jac(3,1) = &
          gen_auxvar%d%por_p * gen_auxvar%den(1) * gen_auxvar%U(1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%U(1) + &
          porosity * gen_auxvar%den(1) * gen_auxvar%d%Ul_pl + &
          (-1.d0) * gen_auxvar%d%por_p * &
            material_auxvar%soil_particle_density * &
            soil_heat_capacity * gen_auxvar%temp
        ! w/respect to air mole fraction
        Jac(3,2) = 0.d0
        ! w/respect to temperature
        ! por * ddenl_dT * Ul + 
        ! por * denl * dUl_dT + 
        ! (1-por) * dens * Cp
        Jac(3,3) = &
          porosity * gen_auxvar%d%denl_T * gen_auxvar%U(1) + &
          porosity * gen_auxvar%den(1) * gen_auxvar%d%Ul_T + &
          (1.d0 - porosity) * material_auxvar%soil_particle_density * &
            soil_heat_capacity
      case(GAS_STATE)
        ! satg = 1
        ! ----------
        ! Water Equation
        ! por * satg * deng * Xwg
        ! ---
        ! w/respect to gas pressure
        ! dpor_dp * deng * Xwg + 
        ! por * ddeng_dpg * Xwg
        Jac(1,1) = &
          gen_auxvar%d%por_p * gen_auxvar%den(2) * gen_auxvar%xmol(1,2) + &
          porosity * gen_auxvar%d%deng_pg * gen_auxvar%xmol(1,2) + &
          ! xmol_p(1,2) = dXwg_dpg
          porosity * gen_auxvar%den(2) * gen_auxvar%d%xmol_p(1,2)
        ! w/respect to air pressure
        ! por * deng_pa * Xwg + 
        ! por * deng * Xwg_pa
        Jac(1,2) = porosity * gen_auxvar%d%deng_pa * gen_auxvar%xmol(1,2) + &
          ! liquid phase index hijacked for air pressure derivative
          ! xmol_p(1,1) = dXwg_dpa
          porosity * gen_auxvar%den(2) * gen_auxvar%d%xmol_p(1,1)
        ! w/repect to temperature
        ! por * ddenl_dT * Xwl
        Jac(1,3) = porosity * gen_auxvar%d%deng_T * gen_auxvar%xmol(1,2)
        ! ----------
        ! Air Equation
        ! por * satg * deng * Xag
        ! w/respect to gas pressure
        ! dpor_dpl * denl * Xal + 
        ! por * ddenl_dpl * Xal
        Jac(2,1) = &
          gen_auxvar%d%por_p * gen_auxvar%den(2) * gen_auxvar%xmol(2,2) + &
          porosity * gen_auxvar%d%deng_pg * gen_auxvar%xmol(2,2) + &
          ! por * ddeng_dpg * Xag
          ! xmol_p(2,2) = dXag_dpg
          porosity *  gen_auxvar%den(2) * gen_auxvar%d%xmol_p(2,2) 
        ! w/respect to air pressure
        Jac(2,2) = porosity * gen_auxvar%d%deng_pa * gen_auxvar%xmol(2,2) + &
          ! liquid phase index hijacked for air pressure derivative
          ! por * ddeng_dpa * Xag
          ! xmol_p(2,1) = dXag_dpa
          porosity * gen_auxvar%den(2) * gen_auxvar%d%xmol_p(2,1)
        ! w/repect to temperature
        ! por * ddeng_dT * Xag
        Jac(2,3) = porosity * gen_auxvar%d%deng_T * gen_auxvar%xmol(2,2)
        ! ----------
        ! Energy Equation
        ! por * satg * deng * Ug + (1-por) * denr * Cp * T
        ! w/respect to gas pressure
        ! dpor_dpg * deng * Ug + 
        ! por * ddeng_dpg * Ug + 
        ! por * deng * dUg_dpg + 
        ! -dpor_dpg * denr * Cp * T
        Jac(3,1) = &
          gen_auxvar%d%por_p * gen_auxvar%den(2) * gen_auxvar%U(2) + &
          porosity * gen_auxvar%d%deng_pg * gen_auxvar%U(2) + &
          porosity * gen_auxvar%den(2) * gen_auxvar%d%Ug_pg + &
          (-1.d0) * gen_auxvar%d%por_p * &
            material_auxvar%soil_particle_density * &
            soil_heat_capacity * gen_auxvar%temp
        ! w/respect to air pressure
        ! por * ddeng_dpa * Ug + 
        ! por * deng * dUg_dpa
        Jac(3,2) = porosity * gen_auxvar%d%deng_pa * gen_auxvar%U(2) + &
          porosity * gen_auxvar%den(2) * gen_auxvar%d%Ug_pa
        ! w/respect to temperature
        ! por * ddeng_dT * Ug + 
        ! por * deng * dUg_dT + 
        ! (1-por) * dens * Cp
        Jac(3,3) = &
          porosity * gen_auxvar%d%denl_T * gen_auxvar%U(2) + &
          porosity * gen_auxvar%den(2) * gen_auxvar%d%Ug_T + &
          (1.d0 - porosity) * material_auxvar%soil_particle_density * &
            soil_heat_capacity
      case(TWO_PHASE_STATE)
        ! ----------
        ! Water Equation
        ! por * (satl * denl * Xwl + satg * deng * Xwg)
        ! ---
        ! w/respect to gas pressure
        ! dpor_dp * (satl * denl * Xwl + satg * deng * Xwg) +
        ! por * (satl * ddenl_dpg * Xwl + satg * ddeng_dpg * Xwg) + 
        ! por * (satl * denl * dXwl_dpg + satg * deng * dXwg_dpg)
        Jac(1,1) = &
          gen_auxvar%d%por_p * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%xmol(1,1) + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%xmol(1,2)) + &
          porosity * &
                                 ! denl_pl = denl_pg
            (gen_auxvar%sat(1) * gen_auxvar%d%denl_pl * gen_auxvar%xmol(1,1) + &
             gen_auxvar%sat(2) * gen_auxvar%d%deng_pg * gen_auxvar%xmol(1,2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%d%xmol_p(1,1) + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%d%xmol_p(1,2))
        ! w/respect to gas saturation
        ! porosity, density, mole fraction are independent of gas saturation
        ! por * (dsatl_dsatg * denl * Xwl + dsatg_dsatg * deng * Xwg)
        ! dsatl_dsatg = -1.
        ! dsatg_dsatg = 1.
        Jac(1,2) = porosity * &
          (-1.d0 * gen_auxvar%den(1) * gen_auxvar%xmol(1,1) + &
            1.d0 * gen_auxvar%den(2) * gen_auxvar%xmol(1,2))
        ! w/repect to temperature
        ! porosity, saturation, mole fraction are independent of temperature
        ! por * (satl * ddenl_T * Xwl + satg * ddeng_T * Xwg) +
        ! por * (satl * denl * dXwl_dpg + satg * deng * dXwg_dpg)
        Jac(1,3) = porosity * &
          (gen_auxvar%sat(1) * gen_auxvar%d%denl_T * gen_auxvar%xmol(1,1) + &
           gen_auxvar%sat(2) * gen_auxvar%d%deng_T * gen_auxvar%xmol(1,2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%d%xmol_T(1,1) + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%d%xmol_T(1,2))
        ! ----------
        ! Gas Equation
        ! por * (satl * denl * Xal + satg * deng * Xag)
        ! ---
        ! w/respect to gas pressure
        ! dpor_dp * (satl * denl * Xal + satg * deng * Xag) +
        ! por * (satl * ddenl_dpg * Xal + satg * ddeng_dpg * Xag) + 
        ! por * (satl * denl * dXal_dpg + satg * deng * dXag_dpg)
        Jac(2,1) = &
          gen_auxvar%d%por_p * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%xmol(2,1) + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%xmol(2,2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%d%denl_pl * gen_auxvar%xmol(2,1) + &
             gen_auxvar%sat(2) * gen_auxvar%d%deng_pg * gen_auxvar%xmol(2,2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%d%xmol_p(2,1) + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%d%xmol_p(2,2))
        ! w/respect to gas saturation
        ! porosity, density, mole fraction are independent of gas saturation
        ! por * (dsatl_dsatg * denl * Xal + dsatg_dsatg * deng * Xag)
        ! dsatl_dsatg = -1.
        ! dsatg_dsatg = 1.
        Jac(2,2) = porosity * &
          (-1.d0 * gen_auxvar%den(1) * gen_auxvar%xmol(2,1) + &
            1.d0 * gen_auxvar%den(2) * gen_auxvar%xmol(2,2))
        ! w/repect to temperature
        ! porosity, saturation, mole fraction are independent of temperature
        ! por * (satl * ddenl_T * Xal + satg * ddeng_T * Xag)
        Jac(2,3) = porosity * &
          (gen_auxvar%sat(1) * gen_auxvar%d%denl_T * gen_auxvar%xmol(2,1) + &
           gen_auxvar%sat(2) * gen_auxvar%d%deng_T * gen_auxvar%xmol(2,2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%d%xmol_T(2,1) + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%d%xmol_T(2,2))           
        ! ----------
        ! Energy Equation
        ! por * (satl * denl * Ul + satg * deng * Ug) + 
        ! (1-por) * denr * Cpr * T
        ! ---
        ! w/respect to gas pressure
        ! dpor_dp * (satl * denl * Ul + satg * deng * Ug) + 
        ! por * (satl * denl_dp * Ul + satg * deng_dp * Ug) + 
        ! por * (satl * denl * dUl_dp + satg * deng * dUg_dp) +
        ! -dpor_dp * denr * Cpr * T
        Jac(3,1) = &
          gen_auxvar%d%por_p * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%U(1) + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%U(2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%d%denl_pl * gen_auxvar%U(1) + &
             gen_auxvar%sat(2) * gen_auxvar%d%deng_pg * gen_auxvar%U(2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%d%Ul_pl + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%d%Ug_pg) + &
          (-1.d0) * gen_auxvar%d%por_p * &
            material_auxvar%soil_particle_density * &
            soil_heat_capacity * gen_auxvar%temp
        ! w/respect to gas saturation
        ! por * (dsatl_dsatg * denl * Ul + dsatg_dsatg * deng * Ug)
        Jac(3,2) = porosity * &
          (-1.d0 * gen_auxvar%den(1) * gen_auxvar%U(1) + &
            1.d0 * gen_auxvar%den(2) * gen_auxvar%U(2))
        ! w/respect to temperature
        Jac(3,3) = &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%d%denl_T * gen_auxvar%U(1) + &
             gen_auxvar%sat(2) * gen_auxvar%d%deng_T * gen_auxvar%U(2)) + &
          porosity * &
            (gen_auxvar%sat(1) * gen_auxvar%den(1) * gen_auxvar%d%Ul_T + &
             gen_auxvar%sat(2) * gen_auxvar%den(2) * gen_auxvar%d%Ug_T) + &
          (1.d0 - porosity) * &
            material_auxvar%soil_particle_density * &
            soil_heat_capacity             
    end select
    Jac = Jac * volume_over_dt
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,7es24.15)') 'accum:', Res
  endif
#endif                    

end subroutine GeneralAccumulation

! ************************************************************************** !

subroutine GeneralFlux(gen_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       thermal_conductivity_up, &
                       gen_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area, dist, general_parameter, &
                       option,v_darcy,Res,Jup,Jdn, &
                       analytical_derivatives, &
                       debug_connection)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  use Fracture_module
  use Klinkenberg_module
  
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: iphase
  
  PetscReal :: xmol(option%nflowspec)
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass
  PetscReal :: delta_X_whatever
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, q
  PetscReal :: tot_mole_flux, wat_mole_flux, air_mole_flux
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  
  PetscReal :: dummy_dperm_up, dummy_dperm_dn
  PetscReal :: temp_perm_up, temp_perm_dn

  ! Darcy flux
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dpaup, ddelta_pressure_dpadn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  
  PetscReal :: up_scale, dn_scale
  PetscReal :: tot_mole_flux_ddel_pressure
  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dtot_mole_flux_dp, dtot_mole_flux_dT, dtot_mole_flux_dsatg
  PetscReal :: dpl_dsatg
  PetscReal :: ddelta_pressure_pl

  ! Diffusion
  PetscReal :: dstpd_up_dporup, dstpd_dn_dpordn
  PetscReal :: dstpd_up_dsatup, dstpd_dn_dsatdn
  PetscReal :: dstpd_up_ddenup, dstpd_dn_ddendn
  PetscReal :: dsatup, dsatdn
  PetscReal :: delta_X_whatever_dxmolup, delta_X_whatever_dxmoldn
  PetscReal :: dxmass_air_up_dxmol_air_up, dxmass_air_dn_dxmol_air_dn
  PetscReal :: dtot_mole_flux_dstpd, dtot_mole_flux_ddeltaX
  PetscReal :: dtot_mole_flux_ddenave
  PetscReal :: diffusion_scale
  PetscReal :: ddiffusion_coef_dTup, ddiffusion_coef_dTdn
  PetscReal :: ddiffusion_coef_dpup, ddiffusion_coef_dpdn
  PetscReal :: dtot_mole_flux_ddiffusion_coef
  PetscReal :: dstpd_ave_over_dist_dstpd_up, dstpd_ave_over_dist_dstpd_dn
  
  ! Conduction
  PetscReal :: dkeff_up_dsatlup, dkeff_dn_dsatldn
  PetscReal :: dkeff_ave_dkeffup, dkeff_ave_dkeffdn
  PetscReal :: dheat_flux_ddelta_temp, dheat_flux_dkeff_ave
  
  ! DELETE
  
  PetscReal :: Jlup(3,3), Jldn(3,3)
  PetscReal :: Jgup(3,3), Jgdn(3,3)
  PetscReal :: Jcup(3,3), Jcdn(3,3)
   
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_up%fracture)) then
    call FracturePermEvaluate(material_auxvar_up,perm_up,temp_perm_up, &
                              dummy_dperm_up,dist)
    perm_up = temp_perm_up
  endif
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif
  
  if (associated(klinkenberg)) then
    perm_ave_over_dist(1) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
    temp_perm_up = klinkenberg%Evaluate(perm_up, &
                                         gen_auxvar_up%pres(option%gas_phase))
    temp_perm_dn = klinkenberg%Evaluate(perm_dn, &
                                         gen_auxvar_dn%pres(option%gas_phase))
    perm_ave_over_dist(2) = (temp_perm_up * temp_perm_dn) / &
                            (dist_up*temp_perm_dn + dist_dn*temp_perm_up)
  else
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
  endif
      
  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0  
  
  v_darcy = 0.d0

#ifdef CONVECTION
#ifdef LIQUID_DARCY_FLUX
  iphase = LIQUID_PHASE
  if (gen_auxvar_up%mobility(iphase) + &
      gen_auxvar_dn%mobility(iphase) > eps) then
    
    density_kg_ave = GeneralAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           gen_auxvar_up%den_kg, &
                                           gen_auxvar_dn%den_kg, &
                                           ddensity_kg_ave_dden_kg_up, &
                                           ddensity_kg_ave_dden_kg_dn)

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = gen_auxvar_up%pres(iphase) - &
                     gen_auxvar_dn%pres(iphase) + &
                     gravity_term
    if (analytical_derivatives) then
      ddelta_pressure_dpup = 1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_up * &
                             gen_auxvar_up%d%denl_pl * fmw_comp(iphase)
      ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_dn * &
                             gen_auxvar_dn%d%denl_pl * fmw_comp(iphase)
      ddelta_pressure_dTup = dist_gravity * ddensity_kg_ave_dden_kg_up * &
                             gen_auxvar_up%d%denl_T * fmw_comp(iphase)
      ddelta_pressure_dTdn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                             gen_auxvar_dn%d%denl_T * fmw_comp(iphase)
    endif
    up_scale = 0.d0
    dn_scale = 0.d0
    if (delta_pressure >= 0.d0) then
      up_scale = 1.d0
      mobility = gen_auxvar_up%mobility(iphase)
      xmol(:) = gen_auxvar_up%xmol(:,iphase)
      uH = gen_auxvar_up%H(iphase)
    else
      dn_scale = 1.d0
      mobility = gen_auxvar_dn%mobility(iphase)
      xmol(:) = gen_auxvar_dn%xmol(:,iphase)
      uH = gen_auxvar_dn%H(iphase)
    endif      

    if (mobility > floweps ) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = GeneralAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          gen_auxvar_up%den, &
                                          gen_auxvar_dn%den, &
                                          ddensity_ave_dden_up, &
                                          ddensity_ave_dden_dn)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      tot_mole_flux = q*density_ave
      tot_mole_flux_ddel_pressure = perm_ave_over_dist(iphase) * &
                                       mobility * area * density_ave
      ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
      !                                 xmol[kmol comp/kmol phase]
      wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
      air_mole_flux = tot_mole_flux * xmol(air_comp_id)
      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      Res(energy_id) = Res(energy_id) + tot_mole_flux * uH
      
      if (analytical_derivatives) then
        Jlup = 0.d0
        Jldn = 0.d0
        select case(global_auxvar_up%istate)
          case(LIQUID_STATE)
            ! derivative wrt liquid pressure
            ! derivative total mole flux wrt liquid pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%denl_pl + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityl_pl + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
            ! derivative water wrt liquid pressure
            Jlup(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp
            ! derivative air wrt liquid pressure
            Jlup(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp
            ! derivative energy wrt liquid pressure
            Jlup(3,1) = uH * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hl_pl
            
            ! derivative wrt air mole fraction
            ! derivative water wrt air mole fraction
            Jlup(1,2) = -1.d0 * up_scale * tot_mole_flux
            ! derivative air wrt air mole fraction
            Jlup(2,2) = 1.d0 * up_scale * tot_mole_flux
            ! derivative energy wrt air mole fraction
            ! Jlup(3,2) = 0.d0
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%denl_T + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityl_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTup
            ! derivative water wrt temperature
            Jlup(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jlup(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jlup(3,3) = uH * dtot_mole_flux_dT + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hl_T
                     
          case(GAS_STATE)
            ! derivative wrt gas pressure
            ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
            !   liquid pressure derivatives.
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%denl_pl + &
              ! mole fraction has to be added in below since it differs for air 
              ! and water            
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityl_pl + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
            ! derivative water wrt gas pressure
            Jlup(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jlup(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jlup(3,1) = uH * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hl_pl

            ! derivative wrt air pressure
            ! derivative water wrt air saturation
            ! Jlup(1,2) = 0.d0
            ! derivative air wrt air saturation
            ! Jlup(2,2) = 0.d0
            ! derivative energy wrt air saturation
            ! Jlup(3,2) = 0.d0
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%denl_T + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityl_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTup
              ! there is no derivative of mole fraction wrt temperature in
              ! gas state
            ! derivative water wrt temperature
            Jlup(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jlup(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jlup(3,3) = dtot_mole_flux_dT * uH + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hl_T 
                     
          case(TWO_PHASE_STATE)
            ! derivative wrt gas pressure
            ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
            !   liquid pressure derivatives.
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%denl_pl + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityl_pl + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
            ! derivative water wrt gas pressure
            Jlup(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jlup(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jlup(3,1) = uH * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hl_pl
            
            ! derivative wrt gas saturation
            ! pl = pg - pc(satg)
            dpl_dsatg = -1.d0 * gen_auxvar_up%d%pc_satg
            ! delta pressure = plup - pldn
            ddelta_pressure_pl = 1.d0          
            ! derivative total mole flux wrt gas saturation
            dtot_mole_flux_dsatg = &
              ! liquid viscosity
              ! since liquid viscosity in a two phase state is a function
              ! of total pressure (gas pressure), there is not derivative
              ! wrt gas saturation
              !up_scale * &
              !tot_mole_flux / mobility * &
              !gen_auxvar_up%d%mobilityl_pl * dpl_dsatg + &
              ! relative permeability
              up_scale * &
              tot_mole_flux / mobility * &
              gen_auxvar_up%d%mobilityl_satg + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_pl * dpl_dsatg         
            ! derivative water wrt gas saturation
            Jlup(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dsatg
            ! derivative air wrt gas saturation
            Jlup(2,2) = xmol(air_comp_id) * dtot_mole_flux_dsatg
            ! derivative energy wrt gas saturation
            Jlup(3,2) = dtot_mole_flux_dsatg * uH
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%denl_T + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityl_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTup
            ! derivative water wrt temperature
            Jlup(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_T(wat_comp_id,iphase)
            ! derivative air wrt temperature
            Jlup(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_T(air_comp_id,iphase)
            ! derivative energy wrt temperature
            Jlup(3,3) = dtot_mole_flux_dT * uH + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hl_T        
        end select
        select case(global_auxvar_dn%istate)
          case(LIQUID_STATE)
            ! derivative wrt liquid pressure
            ! derivative total mole flux wrt liquid pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn *gen_auxvar_dn%d%denl_pl + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityl_pl + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
            ! derivative water wrt liquid pressure
            Jldn(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp
            ! derivative air wrt liquid pressure
            Jldn(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp
            ! derivative energy wrt liquid pressure
            Jldn(3,1) = uH * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hl_pl
            
            ! derivative wrt air mole fraction
            ! derivative water wrt air mole fraction
            Jldn(1,2) = -1.d0 * dn_scale * tot_mole_flux
            ! derivative air wrt air mole fraction
            Jldn(2,2) = 1.d0 * dn_scale * tot_mole_flux
            ! derivative energy wrt air mole fraction
            ! Jldn(3,2) = 0.d0
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_T + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityl_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
            ! derivative water wrt temperature
            Jldn(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jldn(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jldn(3,3) = uH * dtot_mole_flux_dT + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hl_T
                     
          case(GAS_STATE)
            ! derivative wrt gas pressure
            ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
            !   liquid pressure derivatives.
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_pl + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityl_pl + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
            ! derivative water wrt gas pressure
            Jldn(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jldn(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jldn(3,1) = uH * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hl_pl

            ! derivative wrt air pressure
            ! derivative water wrt air saturation
            ! Jldn(1,2) = 0.d0
            ! derivative air wrt air saturation
            ! Jldn(2,2) = 0.d0
            ! derivative energy wrt air saturation
            ! Jldn(3,2) = 0.d0
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_T + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityl_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
              ! there is no derivative of mole fraction wrt temperature in
              ! gas state            
            ! derivative water wrt temperature
            Jldn(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jldn(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jldn(3,3) = dtot_mole_flux_dT * uH + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hl_T 
                     
          case(TWO_PHASE_STATE)
            ! derivative wrt gas pressure
            ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
            !   liquid pressure derivatives.
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn *gen_auxvar_dn%d%denl_pl + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityl_pl + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
            ! derivative water wrt gas pressure
            Jldn(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jldn(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jldn(3,1) = uH * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hl_pl
            
            ! derivative wrt gas saturation
            ! pl = pg - pc(satg)
            dpl_dsatg = -1.d0 * gen_auxvar_dn%d%pc_satg
            ! delta pressure = plup - pldn
            ddelta_pressure_pl = -1.d0
            ! derivative total mole flux wrt gas saturation
            dtot_mole_flux_dsatg = &
              ! liquid viscosity
              ! since liquid viscosity in a two phase state is a function
              ! of total pressure (gas pressure), there is not derivative
              ! wrt gas saturation
              !dn_scale * &
              !tot_mole_flux / mobility * &
              !gen_auxvar_dn%d%mobilityl_pl * dpl_dsatg + &
              ! relative permeability
              dn_scale * &
              tot_mole_flux / mobility * &
              gen_auxvar_dn%d%mobilityl_satg + &
              !pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_pl * dpl_dsatg
            ! derivative water wrt gas saturation
            Jldn(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dsatg
            ! derivative air wrt gas saturation
            Jldn(2,2) = xmol(air_comp_id) * dtot_mole_flux_dsatg
            ! derivative energy wrt gas saturation
            Jldn(3,2) = dtot_mole_flux_dsatg * uH
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_T + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityl_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
            ! derivative water wrt temperature
            Jldn(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_T(wat_comp_id,iphase)
            ! derivative air wrt temperature
            Jldn(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_T(air_comp_id,iphase)
            ! derivative energy wrt temperature
            Jldn(3,3) = dtot_mole_flux_dT * uH + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hl_T        
        end select
        Jup = Jup + Jlup
        Jdn = Jdn + Jldn
      endif
    endif                   
  endif
#endif
#ifdef GAS_DARCY_FLUX
  iphase = GAS_PHASE
  if (gen_auxvar_up%mobility(iphase) + &
      gen_auxvar_dn%mobility(iphase) > eps) then
    
    density_kg_ave = GeneralAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           gen_auxvar_up%den_kg, &
                                           gen_auxvar_dn%den_kg, &
                                           ddensity_kg_ave_dden_kg_up, &
                                           ddensity_kg_ave_dden_kg_dn)

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = gen_auxvar_up%pres(iphase) - &
                     gen_auxvar_dn%pres(iphase) + &
                     gravity_term
    ! if a gas phase does not exist on either side of the connection, the gas
    ! phase properties from the opposite side are used.
    if (analytical_derivatives) then
      ddelta_pressure_dpup = 1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_up * &
                             gen_auxvar_up%d%deng_pg * fmw_comp(iphase)
      ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_dn * &
                             gen_auxvar_dn%d%deng_pg * fmw_comp(iphase)
      ddelta_pressure_dpaup = dist_gravity * ddensity_kg_ave_dden_kg_up * &
                              gen_auxvar_up%d%deng_pa * fmw_comp(iphase)
      ddelta_pressure_dpadn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                              gen_auxvar_dn%d%deng_pa * fmw_comp(iphase)
      ddelta_pressure_dTup = dist_gravity * ddensity_kg_ave_dden_kg_up * &
                             gen_auxvar_up%d%deng_T * fmw_comp(iphase)
      ddelta_pressure_dTdn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                             gen_auxvar_dn%d%deng_T * fmw_comp(iphase)
    endif
    up_scale = 0.d0
    dn_scale = 0.d0
    if (delta_pressure >= 0.d0) then
      up_scale = 1.d0
      mobility = gen_auxvar_up%mobility(iphase)
      xmol(:) = gen_auxvar_up%xmol(:,iphase)
      uH = gen_auxvar_up%H(iphase)
    else
      dn_scale = 1.d0
      mobility = gen_auxvar_dn%mobility(iphase)
      xmol(:) = gen_auxvar_dn%xmol(:,iphase)
      uH = gen_auxvar_dn%H(iphase)
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = GeneralAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          gen_auxvar_up%den, &
                                          gen_auxvar_dn%den, &
                                          ddensity_ave_dden_up, &
                                          ddensity_ave_dden_dn)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      tot_mole_flux = q*density_ave
      tot_mole_flux_ddel_pressure = perm_ave_over_dist(iphase) * &
                                       mobility * area * density_ave      
      ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
      !                                 xmol[kmol comp/kmol phase]
      wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
      air_mole_flux = tot_mole_flux * xmol(air_comp_id)
      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      Res(energy_id) = Res(energy_id) + tot_mole_flux * uH

      if (analytical_derivatives) then
      
        Jgup = 0.d0
        Jgdn = 0.d0
        select case(global_auxvar_up%istate)
          case(LIQUID_STATE)
            ! derivative wrt liquid pressure
            ! derivative total mole flux wrt liquid pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%deng_pg + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_pg + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
            ! derivative water wrt liquid pressure
            Jgup(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp
            ! derivative air wrt liquid pressure
            Jgup(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp
            ! derivative energy wrt liquid pressure
            Jgup(3,1) = uH * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hg_pg
            
            ! derivative wrt air mole fraction
            ! derivative water wrt air mole fraction
            Jgup(1,2) = -1.d0 * up_scale * tot_mole_flux
            ! derivative air wrt air mole fraction
            Jgup(2,2) = 1.d0 * up_scale * tot_mole_flux
            ! derivative energy wrt air mole fraction
            ! Jgup(3,2) = 0.d0
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%deng_T + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTup
            ! derivative water wrt temperature
            Jgup(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jgup(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jgup(3,3) = uH * dtot_mole_flux_dT + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hg_T
                     
          case(GAS_STATE)
            ! derivative wrt gas pressure
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%deng_pg + &
              ! mole fraction has to be added in below since it differs for air 
              ! and water
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_pg + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
            ! derivative water wrt gas pressure
            Jgup(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jgup(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jgup(3,1) = uH * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hg_pg

            ! derivative wrt air pressure
            ! derivative water wrt air saturation
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%deng_pa + &
              ! mole fraction has to be added in below since it differs for air 
              ! and water
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_pa + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpaup
            Jgup(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       ! dXwg_pa for gas phase is stored in liquid phase of xmol_p
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(wat_comp_id,LIQUID_PHASE)
            ! derivative air wrt air saturation
            Jgup(2,2) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       ! dXag_pa for gas phase is stored in liquid phase of xmol_p
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(air_comp_id,LIQUID_PHASE)
            ! derivative energy wrt air saturation
            Jgup(3,2) = uH * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hg_pa
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%deng_T + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTup
              ! there is no derivative of mole fraction wrt temperature in
              ! gas state            
            ! derivative water wrt temperature
            Jgup(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jgup(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jgup(3,3) = dtot_mole_flux_dT * uH + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hg_T 
                     
          case(TWO_PHASE_STATE)
            ! derivative wrt gas pressure
            ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
            !   liquid pressure derivatives.
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%deng_pg + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_pg + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
            ! derivative water wrt gas pressure
            Jgup(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jgup(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jgup(3,1) = uH * dtot_mole_flux_dp + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hg_pg
            
            ! derivative wrt gas saturation
            ! derivative total mole flux wrt gas saturation
            dtot_mole_flux_dsatg = &
              ! relative permeability
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_satg
            ! derivative water wrt gas saturation
            Jgup(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dsatg
            ! derivative air wrt gas saturation
            Jgup(2,2) = xmol(air_comp_id) * dtot_mole_flux_dsatg
            ! derivative energy wrt gas saturation
            Jgup(3,2) = dtot_mole_flux_dsatg * uH
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_up * gen_auxvar_up%d%deng_T + &
              ! liquid mobility
              up_scale * &
              tot_mole_flux / mobility * gen_auxvar_up%d%mobilityg_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTup
            ! derivative water wrt temperature
            Jgup(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_T(wat_comp_id,iphase)
            ! derivative air wrt temperature
            Jgup(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%xmol_T(air_comp_id,iphase)
            ! derivative energy wrt temperature
            Jgup(3,3) = dtot_mole_flux_dT * uH + &
                       up_scale * &
                       tot_mole_flux * gen_auxvar_up%d%Hg_T        
        end select
        select case(global_auxvar_dn%istate)
          case(LIQUID_STATE)
            ! derivative wrt liquid pressure
            ! derivative total mole flux wrt liquid pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pg + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_pg + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
            ! derivative water wrt liquid pressure
            Jgdn(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp
            ! derivative air wrt liquid pressure
            Jgdn(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp
            ! derivative energy wrt liquid pressure
            Jgdn(3,1) = uH * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hg_pg
            
            ! derivative wrt air mole fraction
            ! derivative water wrt air mole fraction
            Jgdn(1,2) = -1.d0 * dn_scale * tot_mole_flux
            ! derivative air wrt air mole fraction
            Jgdn(2,2) = 1.d0 * dn_scale * tot_mole_flux
            ! derivative energy wrt air mole fraction
            ! Jgdn(3,2) = 0.d0
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_T + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
            ! derivative water wrt temperature
            Jgdn(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jgdn(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jgdn(3,3) = uH * dtot_mole_flux_dT + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hg_T
                     
          case(GAS_STATE)
            ! derivative wrt gas pressure
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pg + &
              ! mole fraction has to be added in below since it differs for air 
              ! and water
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_pg + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
            ! derivative water wrt gas pressure
            Jgdn(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jgdn(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jgdn(3,1) = uH * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hg_pg

            ! derivative wrt air pressure
            ! derivative water wrt air saturation
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pa + &
              ! mole fraction has to be added in below since it differs for air 
              ! and water
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_pa + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpadn
            Jgdn(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       ! dXwg_pa for gas phase is stored in liquid phase of xmol_p
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,LIQUID_PHASE)
            ! derivative air wrt air saturation
            Jgdn(2,2) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       ! dXag_pa for gas phase is stored in liquid phase of xmol_p
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,LIQUID_PHASE)
            ! derivative energy wrt air saturation
            Jgdn(3,2) = uH * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hg_pa
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_T + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
              ! there is no derivative of mole fraction wrt temperature in
              ! gas state            
            ! derivative water wrt temperature
            Jgdn(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
            ! derivative air wrt temperature
            Jgdn(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
            ! derivative energy wrt temperature
            Jgdn(3,3) = dtot_mole_flux_dT * uH + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hg_T 
                     
          case(TWO_PHASE_STATE)
            ! derivative wrt gas pressure
            ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
            !   liquid pressure derivatives.
            ! derivative total mole flux wrt gas pressure
            dtot_mole_flux_dp = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pg + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_pg + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
            ! derivative water wrt gas pressure
            Jgdn(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
            ! derivative air wrt gas pressure
            Jgdn(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
            ! derivative energy wrt gas pressure
            Jgdn(3,1) = uH * dtot_mole_flux_dp + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hg_pg
            
            ! derivative wrt gas saturation
            ! derivative total mole flux wrt gas saturation
            dtot_mole_flux_dsatg = &
              ! relative permeability
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_satg
            ! derivative water wrt gas saturation
            Jgdn(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dsatg
            ! derivative air wrt gas saturation
            Jgdn(2,2) = xmol(air_comp_id) * dtot_mole_flux_dsatg
            ! derivative energy wrt gas saturation
            Jgdn(3,2) = dtot_mole_flux_dsatg * uH
          
            ! derivative wrt temperature
            ! derivative total mole flux wrt temperature
            dtot_mole_flux_dT = &
              ! ave. liquid density
              q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_T + &
              ! liquid mobility
              dn_scale * &
              tot_mole_flux / mobility * gen_auxvar_dn%d%mobilityg_T + &
              ! pressure gradient
              tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
            ! derivative water wrt temperature
            Jgdn(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_T(wat_comp_id,iphase)
            ! derivative air wrt temperature
            Jgdn(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%xmol_T(air_comp_id,iphase)
            ! derivative energy wrt temperature
            Jgdn(3,3) = dtot_mole_flux_dT * uH + &
                       dn_scale * &
                       tot_mole_flux * gen_auxvar_dn%d%Hg_T        
        end select
        Jup = Jup + Jgup
        Jdn = Jdn + Jgdn
      endif
    endif               
  endif
#endif  
  ! CONVECTION
#endif

#ifdef DIFFUSION
  ! add in gas component diffusion in gas and liquid phases
!#if 0
#ifdef LIQUID_DIFFUSION  
  iphase = LIQUID_PHASE
  sat_up = gen_auxvar_up%sat(iphase)
  sat_dn = gen_auxvar_dn%sat(iphase)
  dsatup = 1.d0
  dsatdn = 1.d0
  ! by changing #if 1 -> 0, gas component is allowed to diffuse in liquid
  ! phase even if the phase does not exist.
#if 1
  if (sqrt(sat_up*sat_dn) > eps) then
#else
  if (sat_up > eps .or. sat_dn > eps) then
    ! for now, if liquid state neighboring gas, we allow for minute
    ! diffusion in liquid phase.
    if (iphase == option%liquid_phase) then
      if ((sat_up > eps .or. sat_dn > eps)) then
        ! sat_up = max(sat_up,eps)
        if (sat_up < eps) then
          sat_up = eps
          dsatup = 0.d0
        endif
        ! sat_dn = max(sat_dn,eps)
        if (sat_dn < eps) then
          sat_dn = eps
          dsatdn = 0.d0
        endif
      endif
    endif
#endif
    if (general_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_up = gen_auxvar_up%den(iphase)
      den_dn = gen_auxvar_dn%den(iphase)
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_up = 1.d0
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      density_ave = GeneralAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      gen_auxvar_up%den, &
                                      gen_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_up = sat_up*material_auxvar_up%tortuosity* &
              gen_auxvar_up%effective_porosity*den_up
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              gen_auxvar_dn%effective_porosity*den_dn
              
    dstpd_up_dporup = stpd_up / gen_auxvar_up%effective_porosity
    dstpd_dn_dpordn = stpd_dn / gen_auxvar_dn%effective_porosity
    dstpd_up_dsatup = stpd_up / sat_up
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_up_ddenup = tempreal * stpd_up / den_up
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    tempreal = stpd_up*dist_dn+stpd_dn*dist_up
    stpd_ave_over_dist = stpd_up*stpd_dn/tempreal
    dstpd_ave_over_dist_dstpd_up = (stpd_dn-stpd_ave_over_dist*dist_dn)/tempreal
    dstpd_ave_over_dist_dstpd_dn = (stpd_up-stpd_ave_over_dist*dist_up)/tempreal
    
    if (general_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                   gen_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmolup = 1.d0
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      dxmass_air_up_dxmol_air_up = (fmw_comp(2) - xmass_air_up * (fmw_comp(2) - fmw_comp(1))) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      dxmass_air_dn_dxmol_air_dn = (fmw_comp(2) - xmass_air_dn * (fmw_comp(2) - fmw_comp(1))) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmolup = 1.d0 * dxmass_air_up_dxmol_air_up
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             general_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
    
    if (analytical_derivatives) then
    
      Jlup = 0.d0
      Jldn = 0.d0
      select case(global_auxvar_up%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            (dstpd_up_dporup * gen_auxvar_up%d%por_p + &
            ! if density harmonic averaged
             dstpd_up_ddenup * gen_auxvar_up%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%denl_pl
          ! derivative water wrt liquid pressure
          Jlup(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jlup(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jlup(3,1) = 0.d0
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jlup(1,2) = -1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup
          ! derivative air wrt air mole fraction
          Jlup(2,2) = 1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup
          ! derivative energy wrt air mole fraction
          ! Jlup(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_ddenup * gen_auxvar_up%d%denl_T + &
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%denl_T
            ! diffusion coefficient derivative wrt temperature
          ! derivative water wrt temperature
          Jlup(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jlup(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jlup(3,3) = 0.d0
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            (dstpd_up_dporup * gen_auxvar_up%d%por_p + &
             dstpd_up_ddenup * gen_auxvar_up%d%denl_pl) + &
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%denl_pl
          ! derivative water wrt gas pressure
          Jlup(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jlup(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jlup(3,1) = 0.d0

          ! derivative wrt air pressure
          dtot_mole_flux_dp = &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup / &
            gen_auxvar_up%d%Hc
          ! derivative water wrt air saturation
           Jlup(1,2) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt air saturation
           Jlup(2,2) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt air saturation
           Jlup(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_ddenup * gen_auxvar_up%d%denl_T + &
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%denl_T + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup  * &
            (-1.d0) * gen_auxvar_up%xmol(air_comp_id,LIQUID_PHASE) / &
            gen_auxvar_up%d%Hc * gen_auxvar_up%d%Hc_T
          ! diffusion coefficient derivative wrt temperature          
          ! derivative water wrt temperature
          Jlup(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jlup(2,3) = dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jlup(3,3) = 0.d0
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            (dstpd_up_dporup * gen_auxvar_up%d%por_p + &
             dstpd_up_ddenup * gen_auxvar_up%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%denl_pl + &
            ! air mole fraction
            1.d0 * & ! xmolup - xmoldn, not -1 in docs
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup * &
            gen_auxvar_up%d%xmol_p(air_comp_id,LIQUID_PHASE)
          ! derivative water wrt gas pressure
          Jlup(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jlup(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jlup(3,1) = 0.d0        
            
          ! derivative wrt gas saturation
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_dsatup * dsatup * (-1.d0) ! satl -> satg
          ! derivative water wrt gas saturation
          Jlup(1,2) = -1.d0 * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jlup(2,2) = 1.d0 * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jlup(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_ddenup * gen_auxvar_up%d%denl_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%denl_T + &
            ! dispersion coefficient
            ! air mole fraction
            1.d0 * & ! xmolup - xmoldn, not -1 in docs
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup * &
            gen_auxvar_up%d%xmol_T(air_comp_id,LIQUID_PHASE)          
          ! derivative water wrt temperature
          Jlup(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jlup(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jlup(3,3) = 0.d0
      end select
      select case(global_auxvar_dn%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_pl
          ! derivative water wrt liquid pressure
          Jldn(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jldn(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jldn(3,1) = 0.d0
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jldn(1,2) = -1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative air wrt air mole fraction
          Jldn(2,2) = 1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative energy wrt air mole fraction
          ! Jldn(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%denl_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_T
            ! diffusion coefficient derivative wrt temperature
          ! derivative water wrt temperature
          Jldn(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jldn(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jldn(3,3) = 0.d0
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_pl
          ! derivative water wrt gas pressure
          Jldn(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jldn(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jldn(3,1) = 0.d0

          ! derivative wrt air pressure
          dtot_mole_flux_dp = &
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn / &
            gen_auxvar_dn%d%Hc
          ! derivative water wrt air saturation
           Jldn(1,2) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt air saturation
           Jldn(2,2) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt air saturation
           Jldn(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%denl_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_T + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn  * &
            (-1.d0) * gen_auxvar_dn%xmol(air_comp_id,LIQUID_PHASE) / &
            gen_auxvar_dn%d%Hc * gen_auxvar_dn%d%Hc_T          
          ! diffusion coefficient derivative wrt temperature          
          ! derivative water wrt temperature
          Jldn(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jldn(2,3) = dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jldn(3,3) = 0.d0
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_pl + &
            ! air mole fraction
            1.d0 * & ! xmoldn - xmoldn, not -1 in docs
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_p(air_comp_id,LIQUID_PHASE)
          ! derivative water wrt gas pressure
          Jldn(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jldn(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jldn(3,1) = 0.d0        
            
          ! derivative wrt gas saturation
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_dsatdn * dsatdn * (-1.d0) ! satl -> satg
          ! derivative water wrt gas saturation
          Jldn(1,2) = -1.d0 * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jldn(2,2) = 1.d0 * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jldn(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%denl_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_T + &
            ! dispersion coefficient
            ! air mole fraction
            1.d0 * & ! xmoldn - xmoldn, not -1 in docs
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_T(air_comp_id,LIQUID_PHASE)          
          ! derivative water wrt temperature
          Jldn(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jldn(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jldn(3,3) = 0.d0
      end select
      Jup = Jup + Jlup
      Jdn = Jdn + Jldn
    endif
  endif
#endif
!#if 0
#ifdef GAS_DIFFUSION
  iphase = GAS_PHASE
  sat_up = gen_auxvar_up%sat(iphase)
  sat_dn = gen_auxvar_dn%sat(iphase)
  !geh: i am not sure why both of these conditionals were included.  seems
  !     like the latter would never be false.
  if (sqrt(sat_up*sat_dn) > eps) then
    dsatup = 1.d0
    dsatdn = 1.d0
    if (general_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_up = gen_auxvar_up%den(iphase)
      den_dn = gen_auxvar_dn%den(iphase)
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_up = 1.d0
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      density_ave = GeneralAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      gen_auxvar_up%den, &
                                      gen_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_up = sat_up*material_auxvar_up%tortuosity* &
              gen_auxvar_up%effective_porosity*den_up
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              gen_auxvar_dn%effective_porosity*den_dn
              
    dstpd_up_dporup = stpd_up / gen_auxvar_up%effective_porosity
    dstpd_dn_dpordn = stpd_dn / gen_auxvar_dn%effective_porosity
    dstpd_up_dsatup = stpd_up / sat_up
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_up_ddenup = tempreal * stpd_up / den_up
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    tempreal = stpd_up*dist_dn+stpd_dn*dist_up
    stpd_ave_over_dist = stpd_up*stpd_dn/tempreal
    dstpd_ave_over_dist_dstpd_up = (stpd_dn-stpd_ave_over_dist*dist_dn)/tempreal
    dstpd_ave_over_dist_dstpd_dn = (stpd_up-stpd_ave_over_dist*dist_up)/tempreal
    
    if (general_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                   gen_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmolup = 1.d0
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      dxmass_air_up_dxmol_air_up = (fmw_comp(2) - xmass_air_up * (fmw_comp(2) - fmw_comp(1))) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      dxmass_air_dn_dxmol_air_dn = (fmw_comp(2) - xmass_air_dn * (fmw_comp(2) - fmw_comp(1))) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmolup = 1.d0 * dxmass_air_up_dxmol_air_up
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    ! need to account for multiple phases
    ! Eq. 1.9b.  The gas density is added below
    if (general_temp_dep_gas_air_diff) then
      temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
      pressure_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                            gen_auxvar_dn%pres(iphase))
      tempreal = (temp_ave+273.15d0)/273.15d0
      diffusion_scale = tempreal**1.8d0 * 101325.d0 / pressure_ave
                             ! 0.9d0 = 0.5 * 1.8
      ddiffusion_coef_dTup = 0.9d0 * diffusion_scale / (tempreal * 273.15d0)
      ddiffusion_coef_dTdn = ddiffusion_coef_dTup
      ddiffusion_coef_dpup = -1.d0 * diffusion_scale / pressure_ave * 0.5d0
      ddiffusion_coef_dpdn = ddiffusion_coef_dpup
    else
      diffusion_scale = 1.d0
      ddiffusion_coef_dTup = 0.d0
      ddiffusion_coef_dTdn = 0.d0
      ddiffusion_coef_dpup = 0.d0
      ddiffusion_coef_dpdn = 0.d0
    endif
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             diffusion_scale * &
                             general_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddiffusion_coef = tot_mole_flux / diffusion_scale
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave    
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
    
    if (analytical_derivatives) then

      Jgup = 0.d0
      Jgdn = 0.d0
      select case(global_auxvar_up%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            (dstpd_up_dporup * gen_auxvar_up%d%por_p + &
            ! if density harmonic averaged
             dstpd_up_ddenup * gen_auxvar_up%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpup
          ! derivative water wrt liquid pressure
          Jgup(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jgup(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jgup(3,1) = 0.d0
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jgup(1,2) = -1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup
          ! derivative air wrt air mole fraction
          Jgup(2,2) = 1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup
          ! derivative energy wrt air mole fraction
          ! Jgup(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_ddenup * gen_auxvar_up%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTup
          ! derivative water wrt temperature
          Jgup(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jgup(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jgup(3,3) = 0.d0
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            (dstpd_up_dporup * gen_auxvar_up%d%por_p + &
             dstpd_up_ddenup * gen_auxvar_up%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpup + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup * &
            gen_auxvar_up%d%xmol_p(air_comp_id,GAS_PHASE)
          ! derivative water wrt gas pressure
          Jgup(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jgup(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jgup(3,1) = 0.d0

          ! derivative wrt air pressure
          dtot_mole_flux_dp = &
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
             dstpd_up_ddenup * gen_auxvar_up%d%deng_pa + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%deng_pa + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup * &
                              ! liquid phase is hijacked to store \dpa
            gen_auxvar_up%d%xmol_p(air_comp_id,LIQUID_PHASE)
          ! derivative water wrt air saturation
           Jgup(1,2) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt air saturation
           Jgup(2,2) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt air saturation
           Jgup(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_ddenup * gen_auxvar_up%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTup + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup * &
            gen_auxvar_up%d%xmol_T(air_comp_id,GAS_PHASE)  
          ! derivative water wrt temperature
          Jgup(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jgup(2,3) = dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jgup(3,3) = 0.d0
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            (dstpd_up_dporup * gen_auxvar_up%d%por_p + &
             dstpd_up_ddenup * gen_auxvar_up%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpup + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup * &
            gen_auxvar_up%d%xmol_p(air_comp_id,GAS_PHASE)
          ! derivative water wrt gas pressure
          Jgup(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jgup(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jgup(3,1) = 0.d0        
            
          ! derivative wrt gas saturation
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_dsatup * dsatup ! satg
          ! derivative water wrt gas saturation
          Jgup(1,2) = -1.d0 * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jgup(2,2) = 1.d0 * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jgup(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_up * &
            dstpd_up_ddenup * gen_auxvar_up%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_up * &
            gen_auxvar_up%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTup  + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmolup * &
            gen_auxvar_up%d%xmol_T(air_comp_id,GAS_PHASE)          
          ! derivative water wrt temperature
          Jgup(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jgup(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jgup(3,3) = 0.d0
      end select
      select case(global_auxvar_dn%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpdn
          ! derivative water wrt liquid pressure
          Jgdn(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jgdn(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jgdn(3,1) = 0.d0
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jgdn(1,2) = -1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative air wrt air mole fraction
          Jgdn(2,2) = 1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative energy wrt air mole fraction
          ! Jgdn(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTdn
          ! derivative water wrt temperature
          Jgdn(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jgdn(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jgdn(3,3) = 0.d0
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpdn + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_p(air_comp_id,GAS_PHASE)
          ! derivative water wrt gas pressure
          Jgdn(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jgdn(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jgdn(3,1) = 0.d0

          ! derivative wrt air pressure
          dtot_mole_flux_dp = &
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pa + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pa + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
                              ! liquid phase is hijacked to store \dpa
            gen_auxvar_dn%d%xmol_p(air_comp_id,LIQUID_PHASE)
          ! derivative water wrt air saturation
           Jgdn(1,2) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt air saturation
           Jgdn(2,2) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt air saturation
           Jgdn(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTdn + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_T(air_comp_id,GAS_PHASE)  
          ! derivative water wrt temperature
          Jgdn(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jgdn(2,3) = dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jgdn(3,3) = 0.d0
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpdn + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_p(air_comp_id,GAS_PHASE)
          ! derivative water wrt gas pressure
          Jgdn(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jgdn(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jgdn(3,1) = 0.d0        
            
          ! derivative wrt gas saturation
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_dsatdn * dsatdn ! satg
          ! derivative water wrt gas saturation
          Jgdn(1,2) = -1.d0 * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jgdn(2,2) = 1.d0 * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jgdn(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTdn  + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_T(air_comp_id,GAS_PHASE)          
          ! derivative water wrt temperature
          Jgdn(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jgdn(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jgdn(3,3) = 0.d0
      end select
      Jup = Jup + Jgup
      Jdn = Jdn + Jgdn    
    endif
  endif
#endif
! DIFFUSION
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_wet-k_dry)
  ! 1 = dry
  ! 2 = wet
  sat_up = gen_auxvar_up%sat(option%liquid_phase)
  sat_dn = gen_auxvar_dn%sat(option%liquid_phase)
  
  tempreal = sqrt(sat_up) * &
             (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  k_eff_up = thermal_conductivity_up(1) + tempreal
  dkeff_up_dsatlup = 0.5d0 * tempreal / sat_up
  
  tempreal = sqrt(sat_dn) * &
             (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  k_eff_dn = thermal_conductivity_dn(1) + tempreal
  dkeff_dn_dsatldn = 0.5d0 * tempreal / sat_dn

  if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
    tempreal = k_eff_up*dist_dn+k_eff_dn*dist_up
    k_eff_ave = k_eff_up*k_eff_dn/tempreal
    dkeff_ave_dkeffup = (k_eff_dn-k_eff_ave*dist_dn)/tempreal
    dkeff_ave_dkeffdn = (k_eff_up-k_eff_ave*dist_up)/tempreal
  else
    k_eff_ave = 0.d0
    dkeff_ave_dkeffup = 0.d0
    dkeff_ave_dkeffdn = 0.d0
  endif
  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = gen_auxvar_up%temp - gen_auxvar_dn%temp
  dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s
  heat_flux = dheat_flux_ddelta_temp * delta_temp
  dheat_flux_dkeff_ave = area * 1.d-6 * delta_temp
  ! MJ/s or MW
  Res(energy_id) = Res(energy_id) + heat_flux
  
  if (analytical_derivatives) then
    Jcup = 0.d0
    Jcdn = 0.d0
    select case(global_auxvar_up%istate)
      case(LIQUID_STATE,GAS_STATE)
        ! only derivative is energy wrt temperature
        ! derivative energy wrt temperature
        ! positive for upwind
        Jcup(3,3) = 1.d0 * dheat_flux_ddelta_temp
                     
      case(TWO_PHASE_STATE)
        ! only derivatives are energy wrt saturation and temperature
        ! derivative energy wrt gas saturation
        Jcup(3,2) = dheat_flux_dkeff_ave * dkeff_ave_dkeffup * &
                    dkeff_up_dsatlup * (-1.d0) ! satl -> satg
        ! derivative energy wrt temperature
        ! positive for upwind
        Jcup(3,3) = 1.d0 * dheat_flux_ddelta_temp
    end select
    select case(global_auxvar_dn%istate)
      case(LIQUID_STATE,GAS_STATE)
        ! only derivative is energy wrt temperature
        ! derivative energy wrt temperature
        ! positive for upwind
        Jcdn(3,3) = -1.d0 * dheat_flux_ddelta_temp
                     
      case(TWO_PHASE_STATE)
        ! only derivatives are energy wrt saturation and temperature
        ! derivative energy wrt gas saturation
        Jcdn(3,2) = dheat_flux_dkeff_ave * dkeff_ave_dkeffdn * &
                    dkeff_dn_dsatldn * (-1.d0) ! satl -> satg
        ! derivative energy wrt temperature
        ! positive for upwind
        Jcdn(3,3) = -1.d0 * dheat_flux_ddelta_temp
    end select
    Jup = Jup + Jcup
    Jdn = Jdn + Jcdn  
  endif
! CONDUCTION
#endif

end subroutine GeneralFlux

! ************************************************************************** !

subroutine GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         gen_auxvar_up,global_auxvar_up, &
                         gen_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         thermal_conductivity_dn, &
                         area,dist,general_parameter, &
                         option,v_darcy,Res,J, &
                         analytical_derivatives, &
                         debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module                              
  use Material_Aux_class
  use Fracture_module
  use Klinkenberg_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: J(3,3)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_connection
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  PetscInt :: bc_type
  PetscReal :: xmol(option%nflowspec)  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: gravity_term
  PetscReal :: mobility, q 
  PetscReal :: tot_mole_flux
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: boundary_pressure
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass  
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: tempreal
  PetscReal :: delta_X_whatever
  PetscReal :: wat_mole_flux, air_mole_flux

  ! Darcy flux
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dpaup, ddelta_pressure_dpadn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  PetscReal :: dv_darcy_ddelta_pressure
  PetscReal :: dv_darcy_dmobility
  
  PetscReal :: up_scale, dn_scale
  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dtot_mole_flux_dp, dtot_mole_flux_dT, dtot_mole_flux_dsatg
  PetscReal :: dpl_dsatg
  PetscReal :: ddelta_pressure_pl
  PetscReal :: tot_mole_flux_ddel_pressure, tot_mole_flux_dmobility
  PetscReal :: xmol_bool

  ! Diffusion
  PetscReal :: stpd_dn
  PetscReal :: dstpd_dn_dpordn
  PetscReal :: dstpd_dn_dsatdn
  PetscReal :: dstpd_dn_ddendn
  PetscReal :: dsatdn
  PetscReal :: delta_X_whatever_dxmoldn
  PetscReal :: dxmass_air_dn_dxmol_air_dn
  PetscReal :: dtot_mole_flux_dstpd, dtot_mole_flux_ddeltaX
  PetscReal :: dtot_mole_flux_ddenave
  PetscReal :: diffusion_scale
  PetscReal :: ddiffusion_coef_dTdn
  PetscReal :: ddiffusion_coef_dpdn
  PetscReal :: dtot_mole_flux_ddiffusion_coef
  PetscReal :: dstpd_ave_over_dist_dstpd_dn
  PetscReal :: pressure_ave
  
  ! Conduction
  PetscReal :: dkeff_up_dsatlup, dkeff_dn_dsatldn
  PetscReal :: dkeff_ave_dkeffup, dkeff_ave_dkeffdn
  PetscReal :: dheat_flux_ddelta_temp, dheat_flux_dkeff_ave
  
  ! DELETE
  
  PetscReal :: Jl(3,3)
  PetscReal :: Jg(3,3)
  PetscReal :: Jc(3,3)
  
  PetscInt :: idof
  
  PetscReal :: temp_perm_dn
  PetscReal :: dummy_dperm_dn
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  Res = 0.d0
  J = 0.d0
  v_darcy = 0.d0  

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif  
  
  if (associated(klinkenberg)) then
    perm_dn_adj(1) = perm_dn
                                          
    perm_dn_adj(2) = klinkenberg%Evaluate(perm_dn, &
                                          gen_auxvar_dn%pres(option%gas_phase))
  else
    perm_dn_adj(:) = perm_dn
  endif
  
#ifdef CONVECTION  
#ifdef LIQUID_DARCY_FLUX
  iphase = LIQUID_PHASE
  mobility = 0.d0
  xmol_bool = 1.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
      if (gen_auxvar_up%mobility(iphase) + &
          gen_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(GENERAL_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(GENERAL_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        boundary_pressure = gen_auxvar_up%pres(iphase)
        if (iphase == LIQUID_PHASE .and. &
            global_auxvar_up%istate == GAS_STATE) then
          ! the idea here is to accommodate a free surface boundary
          ! face.  this will not work for an interior grid cell as
          ! there should be capillary pressure in force.
          boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
        endif
        density_kg_ave = GeneralAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                gen_auxvar_up%den_kg, &
                                                gen_auxvar_dn%den_kg, &
                                                ddensity_kg_ave_dden_kg_up, &
                                                ddensity_kg_ave_dden_kg_dn)
        ddensity_kg_ave_dden_kg_up = 0.d0 ! always
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          gen_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (analytical_derivatives) then
          ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                                 ddensity_kg_ave_dden_kg_dn * &
                                 gen_auxvar_dn%d%denl_pl * fmw_comp(iphase)
          ddelta_pressure_dTdn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                                 gen_auxvar_dn%d%denl_T * fmw_comp(iphase)
        endif
        if (bc_type == SEEPAGE_BC .or. &
            bc_type == CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              gen_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
            if (analytical_derivatives) then
              option%io_buffer = 'CONDUCTANCE_BC and SEEPAGE_BC need to be &
                &Verified in GeneralBCFlux().'
              call printErrMsg(option)
              ddelta_pressure_dpdn = 0.d0
              ddelta_pressure_dTdn = 0.d0
            endif
          endif
        endif
        dn_scale = 0.d0
        if (delta_pressure >= 0.d0) then
          mobility = gen_auxvar_up%mobility(iphase)
          xmol(:) = gen_auxvar_up%xmol(:,iphase)
          uH = gen_auxvar_up%H(iphase)
        else
          dn_scale = 1.d0        
          mobility = gen_auxvar_dn%mobility(iphase)
          xmol(:) = gen_auxvar_dn%xmol(:,iphase)
          uH = gen_auxvar_dn%H(iphase)
        endif      
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        dv_darcy_ddelta_pressure = perm_ave_over_dist * mobility
        v_darcy(iphase) = dv_darcy_ddelta_pressure * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = GeneralAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            gen_auxvar_up%den, &
                                            gen_auxvar_dn%den, &
                                            ddensity_ave_dden_up, &
                                            ddensity_ave_dden_dn)    
        ddensity_ave_dden_up = 0.d0 ! always
        dv_darcy_dmobility = perm_ave_over_dist * delta_pressure
      endif
    case(NEUMANN_BC)
      xmol_bool = 0.d0
      dv_darcy_ddelta_pressure = 0.d0
      dv_darcy_dmobility = 0.d0
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 0.d0
      ddelta_pressure_dpdn = 0.d0
      ddelta_pressure_dTdn = 0.d0
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(GENERAL_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(GENERAL_GAS_FLUX_INDEX)
      end select
      xmol = 0.d0
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      xmol(iphase) = 1.d0
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = gen_auxvar_up%den(iphase)
          uH = gen_auxvar_up%H(iphase)
        else 
          dn_scale = 1.d0
          density_ave = gen_auxvar_dn%den(iphase)
          uH = gen_auxvar_dn%H(iphase)
          ddensity_ave_dden_dn = 1.d0
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in GeneralBCFlux phase loop.'
      call printErrMsg(option)
  end select
  if (dabs(v_darcy(iphase)) > 0.d0 .or. mobility > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    tot_mole_flux_ddel_pressure = dv_darcy_ddelta_pressure * area * &
                                  density_ave
    tot_mole_flux_dmobility = dv_darcy_dmobility * area * density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
    air_mole_flux = tot_mole_flux * xmol(air_comp_id)
    Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
    Res(energy_id) = Res(energy_id) + tot_mole_flux * uH
      
    if (analytical_derivatives) then
      Jl = 0.d0
      select case(global_auxvar_dn%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn *gen_auxvar_dn%d%denl_pl + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityl_pl + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
          ! derivative water wrt liquid pressure
          Jl(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jl(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jl(3,1) = uH * dtot_mole_flux_dp + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hl_pl
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jl(1,2) = -1.d0 * dn_scale * tot_mole_flux * xmol_bool
          ! derivative air wrt air mole fraction
          Jl(2,2) = 1.d0 * dn_scale * tot_mole_flux * xmol_bool
          ! derivative energy wrt air mole fraction
          ! Jl(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_T + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityl_T + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
          ! derivative water wrt temperature
          Jl(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jl(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jl(3,3) = uH * dtot_mole_flux_dT + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hl_T
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
          !   liquid pressure derivatives.
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_pl + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityl_pl + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
          ! derivative water wrt gas pressure
          Jl(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
          ! derivative air wrt gas pressure
          Jl(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
          ! derivative energy wrt gas pressure
          Jl(3,1) = uH * dtot_mole_flux_dp + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hl_pl

          ! derivative wrt air pressure
          ! derivative water wrt air saturation
          ! Jl(1,2) = 0.d0
          ! derivative air wrt air saturation
          ! Jl(2,2) = 0.d0
          ! derivative energy wrt air saturation
          ! Jl(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_T + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityl_T + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
            ! there is no derivative of mole fraction wrt temperature in
            ! gas state            
          ! derivative water wrt temperature
          Jl(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jl(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jl(3,3) = dtot_mole_flux_dT * uH + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hl_T 
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
          !   liquid pressure derivatives.
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn *gen_auxvar_dn%d%denl_pl + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityl_pl + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
          ! derivative water wrt gas pressure
          Jl(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
          ! derivative air wrt gas pressure
          Jl(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
          ! derivative energy wrt gas pressure
          Jl(3,1) = uH * dtot_mole_flux_dp + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hl_pl
            
          ! derivative wrt gas saturation
          ! pl = pg - pc(satg)
          dpl_dsatg = -1.d0 * gen_auxvar_dn%d%pc_satg
          ! delta pressure = plup - pldn
          ddelta_pressure_pl = -1.d0
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            ! liquid viscosity
            ! since liquid viscosity in a two phase state is a function
            ! of total pressure (gas pressure), there is no derivative
            ! wrt gas saturation
            !dn_scale * &
            !tot_mole_flux_dmobility * &
            !gen_auxvar_dn%d%mobilityl_pl * dpl_dsatg + &
            ! relative permeability
            dn_scale * &
            tot_mole_flux_dmobility * &
            gen_auxvar_dn%d%mobilityl_satg + &
            !pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_pl * dpl_dsatg
          ! derivative water wrt gas saturation
          Jl(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jl(2,2) = xmol(air_comp_id) * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jl(3,2) = dtot_mole_flux_dsatg * uH
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%denl_T + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityl_T + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
          ! derivative water wrt temperature
          Jl(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_T(wat_comp_id,iphase)
          ! derivative air wrt temperature
          Jl(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_T(air_comp_id,iphase)
          ! derivative energy wrt temperature
          Jl(3,3) = dtot_mole_flux_dT * uH + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hl_T        
      end select
      J = J + Jl
    endif
  endif                   
#endif
#ifdef GAS_DARCY_FLUX
  iphase = GAS_PHASE
  mobility = 0.d0
  xmol_bool = 1.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
      if (gen_auxvar_up%mobility(iphase) + &
          gen_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(GENERAL_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(GENERAL_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        boundary_pressure = gen_auxvar_up%pres(iphase)
        if (iphase == LIQUID_PHASE .and. &
            global_auxvar_up%istate == GAS_STATE) then
          ! the idea here is to accommodate a free surface boundary
          ! face.  this will not work for an interior grid cell as
          ! there should be capillary pressure in force.
          boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
        endif
        density_kg_ave = GeneralAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                gen_auxvar_up%den_kg, &
                                                gen_auxvar_dn%den_kg, &
                                                ddensity_kg_ave_dden_kg_up, &
                                                ddensity_kg_ave_dden_kg_dn)
        ddensity_kg_ave_dden_kg_up = 0.d0 ! always
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          gen_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (analytical_derivatives) then
          ddelta_pressure_dpadn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                                  gen_auxvar_dn%d%deng_pa * fmw_comp(iphase)
          ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                                 ddensity_kg_ave_dden_kg_dn * &
                                 gen_auxvar_dn%d%deng_pg * fmw_comp(iphase)
          ddelta_pressure_dTdn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                                 gen_auxvar_dn%d%deng_T * fmw_comp(iphase)
        endif
        if (bc_type == SEEPAGE_BC .or. &
            bc_type == CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              gen_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
            if (analytical_derivatives) then
              ddelta_pressure_dpdn = 0.d0
              ddelta_pressure_dTdn = 0.d0
            endif
          endif
        endif
        dn_scale = 0.d0
        ! don't expect the derivative to match precisely at delta_pressure = 0
        ! due to potential switch in direction for numerically perturbed
        ! residual
        if (delta_pressure >= 0.d0) then
          mobility = gen_auxvar_up%mobility(iphase)
          xmol(:) = gen_auxvar_up%xmol(:,iphase)
          uH = gen_auxvar_up%H(iphase)
        else
          dn_scale = 1.d0        
          mobility = gen_auxvar_dn%mobility(iphase)
          xmol(:) = gen_auxvar_dn%xmol(:,iphase)
          uH = gen_auxvar_dn%H(iphase)
        endif  
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        dv_darcy_ddelta_pressure = perm_ave_over_dist * mobility
        v_darcy(iphase) = dv_darcy_ddelta_pressure * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = GeneralAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            gen_auxvar_up%den, &
                                            gen_auxvar_dn%den, &
                                            ddensity_ave_dden_up, &
                                            ddensity_ave_dden_dn)    
        ddensity_ave_dden_up = 0.d0 ! always
        dv_darcy_dmobility = perm_ave_over_dist * delta_pressure
      endif
    case(NEUMANN_BC)
      xmol_bool = 0.d0
      dv_darcy_ddelta_pressure = 0.d0
      dv_darcy_dmobility = 0.d0
      ddensity_ave_dden_up = 0.d0 ! always
      ddensity_ave_dden_dn = 0.d0
      ddelta_pressure_dpdn = 0.d0
      ddelta_pressure_dpadn = 0.d0
      ddelta_pressure_dTdn = 0.d0      
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(GENERAL_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(GENERAL_GAS_FLUX_INDEX)
      end select
      xmol = 0.d0
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      xmol(iphase) = 1.d0
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = gen_auxvar_up%den(iphase)
          uH = gen_auxvar_up%H(iphase)
        else 
          dn_scale = 1.d0
          density_ave = gen_auxvar_dn%den(iphase)
          uH = gen_auxvar_dn%H(iphase)
          ddensity_ave_dden_dn = 1.d0
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in GeneralBCFlux phase loop.'
      call printErrMsg(option)
  end select

  if (dabs(v_darcy(iphase)) > 0.d0 .or. mobility > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    tot_mole_flux_ddel_pressure = dv_darcy_ddelta_pressure * area * &
                                  density_ave
    tot_mole_flux_dmobility = dv_darcy_dmobility * area * density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
    air_mole_flux = tot_mole_flux * xmol(air_comp_id)
    Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
    Res(energy_id) = Res(energy_id) + tot_mole_flux * uH
      
    if (analytical_derivatives) then
      Jg = 0.d0
      select case(global_auxvar_dn%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pg + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_pg + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
          ! derivative water wrt liquid pressure
          Jg(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jg(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jg(3,1) = uH * dtot_mole_flux_dp + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hg_pg
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jg(1,2) = -1.d0 * dn_scale * tot_mole_flux * xmol_bool
          ! derivative air wrt air mole fraction
          Jg(2,2) = 1.d0 * dn_scale * tot_mole_flux * xmol_bool
          ! derivative energy wrt air mole fraction
          ! Jg(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_T + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_T + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
          ! derivative water wrt temperature
          Jg(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jg(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jg(3,3) = uH * dtot_mole_flux_dT + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hg_T
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pg + &
            ! mole fraction has to be added in below since it differs for air 
            ! and water
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_pg + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
          ! derivative water wrt gas pressure
          Jg(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
          ! derivative air wrt gas pressure
          Jg(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
          ! derivative energy wrt gas pressure
          Jg(3,1) = uH * dtot_mole_flux_dp + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hg_pg

          ! derivative wrt air pressure
          ! derivative water wrt air saturation
          dtot_mole_flux_dp = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pa + &
            ! mole fraction has to be added in below since it differs for air 
            ! and water
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_pa + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dpadn
          Jg(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      ! dXwg_pa for gas phase is stored in liquid phase of xmol_p
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,LIQUID_PHASE)
          ! derivative air wrt air saturation
          Jg(2,2) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      ! dXag_pa for gas phase is stored in liquid phase of xmol_p
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,LIQUID_PHASE)
          ! derivative energy wrt air saturation
          Jg(3,2) = uH * dtot_mole_flux_dp + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hg_pa
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_T + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_T + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
            ! there is no derivative of mole fraction wrt temperature in
            ! gas state            
          ! derivative water wrt temperature
          Jg(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jg(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jg(3,3) = dtot_mole_flux_dT * uH + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hg_T 
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
          !   liquid pressure derivatives.
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_pg + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_pg + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
          ! derivative water wrt gas pressure
          Jg(1,1) = xmol(wat_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(wat_comp_id,iphase)
          ! derivative air wrt gas pressure
          Jg(2,1) = xmol(air_comp_id) * dtot_mole_flux_dp + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_p(air_comp_id,iphase)
          ! derivative energy wrt gas pressure
          Jg(3,1) = uH * dtot_mole_flux_dp + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hg_pg
            
          ! derivative wrt gas saturation
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            ! relative permeability
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_satg
          ! derivative water wrt gas saturation
          Jg(1,2) = xmol(wat_comp_id) * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jg(2,2) = xmol(air_comp_id) * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jg(3,2) = dtot_mole_flux_dsatg * uH
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! ave. liquid density
            q * ddensity_ave_dden_dn * gen_auxvar_dn%d%deng_T + &
            ! liquid mobility
            dn_scale * &
            tot_mole_flux_dmobility * gen_auxvar_dn%d%mobilityg_T + &
            ! pressure gradient
            tot_mole_flux_ddel_pressure * ddelta_pressure_dTdn
          ! derivative water wrt temperature
          Jg(1,3) = xmol(wat_comp_id) * dtot_mole_flux_dT + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_T(wat_comp_id,iphase)
          ! derivative air wrt temperature
          Jg(2,3) = xmol(air_comp_id) * dtot_mole_flux_dT + &
                      dn_scale * xmol_bool * &
                      tot_mole_flux * gen_auxvar_dn%d%xmol_T(air_comp_id,iphase)
          ! derivative energy wrt temperature
          Jg(3,3) = dtot_mole_flux_dT * uH + &
                      dn_scale * &
                      tot_mole_flux * gen_auxvar_dn%d%Hg_T        
      end select
      J = J + Jg
    endif
  endif                   
#endif  
! CONVECTION
#endif
  
#ifdef DIFFUSION
#ifdef LIQUID_DIFFUSION  
  iphase = LIQUID_PHASE
  dsatdn = 1.d0
  ! diffusion all depends upon the downwind cell.  phase diffusion only
  ! occurs if a phase exists in both auxvars (boundary and internal) or
  ! a liquid phase exists in the internal cell. so, one could say that
  ! liquid diffusion always exists as the internal cell has a liquid phase,
  ! but gas phase diffusion only occurs if the internal cell has a gas
  ! phase.
  sat_dn = gen_auxvar_dn%sat(iphase)
  if (sat_dn > eps .and. ibndtype(iphase) /= NEUMANN_BC) then
    if (general_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_dn = gen_auxvar_dn%den(iphase)
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      density_ave = GeneralAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      gen_auxvar_up%den, &
                                      gen_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ddensity_ave_dden_up = 0.d0
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              gen_auxvar_dn%effective_porosity*den_dn
              
    dstpd_dn_dpordn = stpd_dn / gen_auxvar_dn%effective_porosity
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    dstpd_ave_over_dist_dstpd_dn = 1.d0 / dist(0)
    stpd_ave_over_dist = stpd_dn * dstpd_ave_over_dist_dstpd_dn
    
    if (general_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                   gen_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      dxmass_air_dn_dxmol_air_dn = (fmw_comp(2) - xmass_air_dn * (fmw_comp(2) - fmw_comp(1))) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             general_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
    
    if (analytical_derivatives) then
      Jl = 0.d0
      select case(global_auxvar_dn%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_pl
          ! derivative water wrt liquid pressure
          Jl(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jl(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jl(3,1) = 0.d0
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jl(1,2) = -1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative air wrt air mole fraction
          Jl(2,2) = 1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative energy wrt air mole fraction
          ! Jl(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%denl_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_T
            ! diffusion coefficient derivative wrt temperature
          ! derivative water wrt temperature
          Jl(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jl(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jl(3,3) = 0.d0
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_pl
          ! derivative water wrt gas pressure
          Jl(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jl(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jl(3,1) = 0.d0

          ! derivative wrt air pressure
          dtot_mole_flux_dp = &
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn / &
            gen_auxvar_dn%d%Hc
          ! derivative water wrt air saturation
           Jl(1,2) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt air saturation
           Jl(2,2) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt air saturation
           Jl(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%denl_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_T + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn  * &
            (-1.d0) * gen_auxvar_dn%xmol(air_comp_id,LIQUID_PHASE) / &
            gen_auxvar_dn%d%Hc * gen_auxvar_dn%d%Hc_T          
          ! diffusion coefficient derivative wrt temperature          
          ! derivative water wrt temperature
          Jl(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jl(2,3) = dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jl(3,3) = 0.d0
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%denl_pl) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_pl + &
            ! air mole fraction
            1.d0 * & ! xmoldn - xmoldn, not -1 in docs
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_p(air_comp_id,LIQUID_PHASE)
          ! derivative water wrt gas pressure
          Jl(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jl(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jl(3,1) = 0.d0        
            
          ! derivative wrt gas saturation
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_dsatdn * dsatdn * (-1.d0) ! satl -> satg
          ! derivative water wrt gas saturation
          Jl(1,2) = -1.d0 * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jl(2,2) = 1.d0 * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jl(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%denl_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%denl_T + &
            ! dispersion coefficient
            ! air mole fraction
            1.d0 * & ! xmoldn - xmoldn, not -1 in docs
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_T(air_comp_id,LIQUID_PHASE)          
          ! derivative water wrt temperature
          Jl(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jl(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jl(3,3) = 0.d0
      end select
      J = J + Jl
    endif
  endif
#endif
#ifdef GAS_DIFFUSION
  iphase = GAS_PHASE
  sat_dn = gen_auxvar_dn%sat(iphase)
  !geh: i am not sure why both of these conditionals were included.  seems
  !     like the latter would never be false.
  if (sat_dn > eps .and. ibndtype(iphase) /= NEUMANN_BC) then
    dsatdn = 1.d0
    if (general_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_dn = gen_auxvar_dn%den(iphase)
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      !TODO(geh): why are we averaging density here?
      density_ave = GeneralAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      gen_auxvar_up%den, &
                                      gen_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ddensity_ave_dden_up = 0.d0                                      
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              gen_auxvar_dn%effective_porosity*den_dn
              
    dstpd_dn_dpordn = stpd_dn / gen_auxvar_dn%effective_porosity
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    dstpd_ave_over_dist_dstpd_dn = 1.d0 / dist(0)
    stpd_ave_over_dist = stpd_dn * dstpd_ave_over_dist_dstpd_dn    
    
    if (general_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                   gen_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    ! need to account for multiple phases
    ! Eq. 1.9b.  The gas density is added below
    if (general_temp_dep_gas_air_diff) then
      temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
      pressure_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                            gen_auxvar_dn%pres(iphase))
      tempreal = (temp_ave+273.15d0)/273.15d0
      diffusion_scale = tempreal**1.8d0 * 101325.d0 / pressure_ave
                             ! 0.9d0 = 0.5 * 1.8
      ddiffusion_coef_dTdn = 0.9d0 * diffusion_scale / (tempreal * 273.15d0)
      ddiffusion_coef_dpdn = -1.d0 * diffusion_scale / pressure_ave * 0.5d0
    else
      diffusion_scale = 1.d0
      ddiffusion_coef_dTdn = 0.d0
      ddiffusion_coef_dpdn = 0.d0
    endif
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             diffusion_scale * &
                             general_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddiffusion_coef = tot_mole_flux / diffusion_scale
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave    
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
    
    if (analytical_derivatives) then

      Jg = 0.d0
      select case(global_auxvar_dn%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative total mole flux wrt liquid pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpdn
          ! derivative water wrt liquid pressure
          Jg(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt liquid pressure
          Jg(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt liquid pressure
          Jg(3,1) = 0.d0
            
          ! derivative wrt air mole fraction
          ! derivative water wrt air mole fraction
          Jg(1,2) = -1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative air wrt air mole fraction
          Jg(2,2) = 1.d0 * dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn
          ! derivative energy wrt air mole fraction
          ! Jg(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTdn
          ! derivative water wrt temperature
          Jg(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jg(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jg(3,3) = 0.d0
                     
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpdn + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_p(air_comp_id,GAS_PHASE)
          ! derivative water wrt gas pressure
          Jg(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jg(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jg(3,1) = 0.d0

          ! derivative wrt air pressure
          dtot_mole_flux_dp = &
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pa + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pa + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
                              ! liquid phase is hijacked to store \dpa
            gen_auxvar_dn%d%xmol_p(air_comp_id,LIQUID_PHASE)
          ! derivative water wrt air saturation
           Jg(1,2) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt air saturation
           Jg(2,2) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt air saturation
           Jg(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = & 
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTdn + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_T(air_comp_id,GAS_PHASE)  
          ! derivative water wrt temperature
          Jg(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jg(2,3) = dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jg(3,3) = 0.d0
                     
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          ! derivative total mole flux wrt gas pressure
          dtot_mole_flux_dp = & 
            ! liquid density and porosity
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            (dstpd_dn_dpordn * gen_auxvar_dn%d%por_p + &
             dstpd_dn_ddendn * gen_auxvar_dn%d%deng_pg) + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_pg + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dpdn + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_p(air_comp_id,GAS_PHASE)
          ! derivative water wrt gas pressure
          Jg(1,1) = -1.d0 * dtot_mole_flux_dp
          ! derivative air wrt gas pressure
          Jg(2,1) = 1.d0 * dtot_mole_flux_dp
          ! derivative energy wrt gas pressure
          Jg(3,1) = 0.d0        
            
          ! derivative wrt gas saturation
          ! derivative total mole flux wrt gas saturation
          dtot_mole_flux_dsatg = &
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_dsatdn * dsatdn ! satg
          ! derivative water wrt gas saturation
          Jg(1,2) = -1.d0 * dtot_mole_flux_dsatg
          ! derivative air wrt gas saturation
          Jg(2,2) = 1.d0 * dtot_mole_flux_dsatg
          ! derivative energy wrt gas saturation
          Jg(3,2) = 0.d0
          
          ! derivative wrt temperature
          ! derivative total mole flux wrt temperature
          dtot_mole_flux_dT = &
            ! liquid density
            dtot_mole_flux_dstpd * dstpd_ave_over_dist_dstpd_dn * &
            dstpd_dn_ddendn * gen_auxvar_dn%d%deng_T + &
            ! if density arithmetically averaged
            dtot_mole_flux_ddenave * ddensity_ave_dden_dn * &
            gen_auxvar_dn%d%deng_T + &
            ! diffusion coefficient
            dtot_mole_flux_ddiffusion_coef * ddiffusion_coef_dTdn  + &
            ! air mole fraction
            dtot_mole_flux_ddeltaX * delta_X_whatever_dxmoldn * &
            gen_auxvar_dn%d%xmol_T(air_comp_id,GAS_PHASE)          
          ! derivative water wrt temperature
          Jg(1,3) = -1.d0 * dtot_mole_flux_dT
          ! derivative air wrt temperature
          Jg(2,3) = 1.d0 * dtot_mole_flux_dT
          ! derivative energy wrt temperature
          Jg(3,3) = 0.d0
      end select
      J = J + Jg
    endif
  endif
#endif
! DIFFUSION
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_wet-k_dry)
  ! 1 = dry
  ! 2 = wet
  heat_flux = 0.d0
  select case (ibndtype(GENERAL_ENERGY_EQUATION_INDEX))
    case (DIRICHLET_BC)
      sat_dn = gen_auxvar_dn%sat(option%liquid_phase)
      tempreal = sqrt(sat_dn) * &
                 (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
      k_eff_dn = thermal_conductivity_dn(1) + tempreal
      dkeff_dn_dsatldn = 0.5d0 * tempreal / sat_dn

      dkeff_ave_dkeffdn = 1.d0 / dist(0)
      k_eff_ave = k_eff_dn * dkeff_ave_dkeffdn
      ! units:
      ! k_eff = W/K-m = J/s/K-m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = k_eff * delta_temp * area = J/s
      delta_temp = gen_auxvar_up%temp - gen_auxvar_dn%temp
      dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s
      heat_flux = dheat_flux_ddelta_temp * delta_temp
      dheat_flux_dkeff_ave = area * 1.d-6 * delta_temp
    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = auxvars(auxvar_mapping(GENERAL_ENERGY_FLUX_INDEX)) * area
      dheat_flux_ddelta_temp = 0.d0
      dkeff_dn_dsatldn = 0.d0
      dkeff_ave_dkeffdn = 0.d0
      dheat_flux_dkeff_ave = 0.d0
    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'GeneralBCFlux heat conduction loop.'
      call printErrMsg(option)
  end select
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux
  
  if (analytical_derivatives) then
    Jc = 0.d0
    select case(global_auxvar_dn%istate)
      case(LIQUID_STATE,GAS_STATE)
        ! only derivative is energy wrt temperature
        ! derivative energy wrt temperature
        ! positive for upwind
        Jc(3,3) = -1.d0 * dheat_flux_ddelta_temp
                     
      case(TWO_PHASE_STATE)
        ! only derivatives are energy wrt saturation and temperature
        ! derivative energy wrt gas saturation
        Jc(3,2) = dheat_flux_dkeff_ave * dkeff_ave_dkeffdn * &
                  dkeff_dn_dsatldn * (-1.d0) ! satl -> satg
        ! derivative energy wrt temperature
        ! positive for upwind
        Jc(3,3) = -1.d0 * dheat_flux_ddelta_temp
    end select
    J = J + Jc
  endif
! CONDUCTION
#endif

end subroutine GeneralBCFlux

! ************************************************************************** !

subroutine GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                          gen_auxvar,global_auxvar,ss_flow_vol_flux, &
                          scale,Res,J,analytical_derivatives,debug_cell)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
  PetscReal :: J(option%nflowdof,option%nflowdof)  
  PetscBool :: analytical_derivatives  
  PetscBool :: debug_cell
      
  PetscReal :: qsrc_mol
  PetscReal :: enthalpy, internal_energy
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscReal :: Jl(option%nflowdof,option%nflowdof)  
  PetscReal :: Jg(option%nflowdof,option%nflowdof)  
  PetscReal :: Je(option%nflowdof,option%nflowdof)  
  PetscReal :: dden_bool
  PetscReal :: hw_dp, hw_dT, ha_dp, ha_dT, dum1, dum2, dum3
  PetscErrorCode :: ierr

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  Res = 0.d0
  J = 0.d0
  
#ifdef WATER_SRCSINK
  qsrc_mol = 0.d0
  dden_bool = 0.d0
  select case(flow_src_sink_type)
    case(MASS_RATE_SS)
      qsrc_mol = qsrc(wat_comp_id)/fmw_comp(wat_comp_id) ! kg/sec -> kmol/sec
    case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
      qsrc_mol = qsrc(wat_comp_id)/fmw_comp(wat_comp_id)*scale 
    case(VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec
      qsrc_mol = qsrc(wat_comp_id)*gen_auxvar%den(wat_comp_id) ! den = kmol/m^3
      dden_bool = 1.d0
    case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
      qsrc_mol = qsrc(wat_comp_id)*gen_auxvar%den(wat_comp_id)*scale
      dden_bool = 1.d0
  end select
  ss_flow_vol_flux(wat_comp_id) = qsrc_mol/gen_auxvar%den(wat_comp_id)
  Res(wat_comp_id) = qsrc_mol
  if (analytical_derivatives) then
    Jl = 0.d0
    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        ! derivative wrt liquid pressure
        Jl(1,1) = dden_bool * qsrc(wat_comp_id) * gen_auxvar%d%denl_pl
        ! derivative wrt air mole fraction
        ! derivative wrt temperature
        Jl(1,3) = dden_bool * qsrc(wat_comp_id) * gen_auxvar%d%denl_T
      case(GAS_STATE)
        option%io_buffer = 'Water injection not set up for gas state in &
          &GeneralSrcSink.'
        call printErrMsg(option)
        ! derivative wrt gas pressure
        ! derivative wrt air pressure
        ! derivative wrt temperature
      case(TWO_PHASE_STATE)
        ! derivative wrt gas pressure
        Jl(1,1) = dden_bool * qsrc(wat_comp_id) * gen_auxvar%d%denl_pl
        ! derivative wrt gas saturation
        ! derivative wrt temperature
        Jl(1,1) = dden_bool * qsrc(wat_comp_id) * gen_auxvar%d%denl_T
    end select
    J = J + Jl
  endif
#endif
#ifdef AIR_SRCSINK
  qsrc_mol = 0.d0
  dden_bool = 0.d0
  select case(flow_src_sink_type)
    case(MASS_RATE_SS)
      qsrc_mol = qsrc(air_comp_id)/fmw_comp(air_comp_id) ! kg/sec -> kmol/sec
    case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
      qsrc_mol = qsrc(air_comp_id)/fmw_comp(air_comp_id)*scale 
    case(VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec
      qsrc_mol = qsrc(air_comp_id)*gen_auxvar%den(air_comp_id) ! den = kmol/m^3
      dden_bool = 1.d0
    case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
      qsrc_mol = qsrc(air_comp_id)*gen_auxvar%den(air_comp_id)*scale
      dden_bool = 1.d0
  end select
  ss_flow_vol_flux(air_comp_id) = qsrc_mol/gen_auxvar%den(air_comp_id)
  Res(air_comp_id) = qsrc_mol
  if (analytical_derivatives) then
    Jg = 0.d0
    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        option%io_buffer = 'Air injection not set up for liquid state in &
          &GeneralSrcSink.'
        call printErrMsg(option)
        ! derivative wrt liquid pressure
        ! derivative wrt air mole fraction
        ! derivative wrt temperature
      case(GAS_STATE)
        ! derivative wrt gas pressure
        Jg(2,1) = dden_bool * qsrc(air_comp_id) * gen_auxvar%d%deng_pg
        ! derivative wrt air pressure
        Jg(2,2) = dden_bool * qsrc(air_comp_id) * gen_auxvar%d%deng_pa
        ! derivative wrt temperature
        Jg(2,3) = dden_bool * qsrc(air_comp_id) * gen_auxvar%d%deng_T
      case(TWO_PHASE_STATE)
        ! derivative wrt gas pressure
        Jg(2,1) = dden_bool * qsrc(air_comp_id) * gen_auxvar%d%deng_pg
        ! derivative wrt gas saturation
        ! derivative wrt temperature
        Jg(2,3) = dden_bool * qsrc(air_comp_id) * gen_auxvar%d%deng_T
    end select
    J = J + Jg
  endif
#endif
  if (dabs(qsrc(TWO_INTEGER)) < 1.d-40 .and. &
      qsrc(ONE_INTEGER) < 0.d0) then ! extraction only
    ! Res(1) holds qsrc_mol for water.  If the src/sink value for air is zero,
    ! remove/add the equivalent mole fraction of air in the liquid phase.
    qsrc_mol = Res(wat_comp_id)*gen_auxvar%xmol(TWO_INTEGER,ONE_INTEGER)
    Res(TWO_INTEGER) = qsrc_mol
    ss_flow_vol_flux(air_comp_id) = qsrc_mol/gen_auxvar%den(TWO_INTEGER)
    if (analytical_derivatives) then
      !Jg = 0.d0
      select case(global_auxvar%istate)
        case(LIQUID_STATE)
          ! derivative wrt liquid pressure
          ! derivative wrt air mole fraction
          Jg(2,2) = Jg(2,2) + Res(ONE_INTEGER)
          ! derivative wrt temperature
        case(GAS_STATE)
          ! derivative wrt gas pressure
          ! derivative wrt air pressure
          ! derivative wrt temperature
        case(TWO_PHASE_STATE)
          ! derivative wrt gas pressure
          Jg(2,1) = Jg(2,1) + &
                    dden_bool * qsrc(wat_comp_id) * gen_auxvar%d%denl_pl * &
                    gen_auxvar%xmol(TWO_INTEGER,ONE_INTEGER) + &
                    Res(ONE_INTEGER) * gen_auxvar%d%xmol_p(2,1)
          ! derivative wrt gas saturation
          ! derivative wrt temperature
          Jg(2,3) = Jg(2,3) + &
                    dden_bool * qsrc(wat_comp_id) * gen_auxvar%d%denl_T * &
                    gen_auxvar%xmol(TWO_INTEGER,ONE_INTEGER) + &
                    Res(ONE_INTEGER) * gen_auxvar%d%xmol_T(2,1)
      end select
      J = J + Jg
    endif
  endif
  ! energy units: MJ/sec
  if (size(qsrc) == THREE_INTEGER) then
    if (dabs(qsrc(THREE_INTEGER)) < 1.d-40) then
      cell_pressure = &
        maxval(gen_auxvar%pres(option%liquid_phase:option%gas_phase))
      if (dabs(qsrc(ONE_INTEGER)) > 0.d0) then
        if (associated(gen_auxvar%d)) then
          call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,enthalpy, &
                                hw_dp,hw_dT,ierr)
          hw_dp = hw_dp * 1.d-6
          hw_dT = hw_dT * 1.d-6
        else
          call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,enthalpy,ierr)
        endif
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
        ! enthalpy units: MJ/kmol                       ! water component mass
        Res(option%energy_id) = Res(option%energy_id) + Res(ONE_INTEGER) * &
                                                        enthalpy
        if (analytical_derivatives) then
          Je = 0.d0
          Je(3,1) = Jl(1,1) * enthalpy + &
                    Res(ONE_INTEGER) * hw_dp
          Je(3,3) = Jl(1,3) * enthalpy + &
                    Res(ONE_INTEGER) * hw_dT
        endif
        J = J + Je
      endif
      if (dabs(qsrc(TWO_INTEGER)) > 0.d0) then
        ! this is pure air, we use the enthalpy of air, NOT the air/water
        ! mixture in gas
        ! air enthalpy is only a function of temperature and the 
        dummy_pressure = 0.d0
        if (associated(gen_auxvar%d)) then
          call EOSGasEnergy(gen_auxvar%temp,dummy_pressure,enthalpy, &
                            ha_dT,ha_dp,dum1,dum2,dum3,ierr)
          ha_dp = ha_dp * 1.d-6
          ha_dT = ha_dT * 1.d-6
        else
          call EOSGasEnergy(gen_auxvar%temp,dummy_pressure, &
                            enthalpy,internal_energy,ierr)
        endif
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> MJ/kmol                                  
        ! enthalpy units: MJ/kmol                       ! air component mass
        Res(option%energy_id) = Res(option%energy_id) + Res(TWO_INTEGER) * &
                                                        enthalpy
        if (analytical_derivatives) then
          Je = 0.d0
          Je(3,1) = Jl(2,1) * enthalpy + &
                    Res(TWO_INTEGER) * ha_dp
          Je(3,2) = Jl(2,2) * enthalpy
          Je(3,3) = Jl(2,3) * enthalpy + &
                    Res(TWO_INTEGER) * ha_dT
        endif
        J = J + Je
      endif
    else
      Res(option%energy_id) = qsrc(THREE_INTEGER)*scale ! MJ/s
      ! no derivative
    endif
  endif
  
end subroutine GeneralSrcSink

! ************************************************************************** !

subroutine GeneralAccumDerivative(gen_auxvar,global_auxvar,material_auxvar, &
                                  soil_heat_capacity,option,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none

  type(general_auxvar_type) :: gen_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: jac(option%nflowdof,option%nflowdof)
  PetscReal :: jac_pert(option%nflowdof,option%nflowdof)
  PetscInt :: idof, irow

!geh:print *, 'GeneralAccumDerivative'

  call GeneralAccumulation(gen_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar,soil_heat_capacity,option, &
                           res,jac,general_analytical_derivatives, &
                           PETSC_FALSE)
                           
  do idof = 1, option%nflowdof
    call GeneralAccumulation(gen_auxvar(idof), &
                             global_auxvar, &
                             material_auxvar,soil_heat_capacity, &
                             option,res_pert,jac_pert,PETSC_FALSE,PETSC_FALSE)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar(idof)%pert
!geh:print *, irow, idof, J(irow,idof), gen_auxvar(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_analytical_derivatives) then
    J = jac
  endif

  if (general_isothermal) then
    J(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    J(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    J(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    J(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,10es24.15)') 'accum deriv:', J
  endif
#endif

end subroutine GeneralAccumDerivative

! ************************************************************************** !

subroutine GeneralFluxDerivative(gen_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 thermal_conductivity_up, &
                                 gen_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 thermal_conductivity_dn, &
                                 area, dist, &
                                 general_parameter, &
                                 option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up(0:), gen_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_up(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_dn(option%nflowdof,option%nflowdof)
  PetscReal :: Jdummy(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0
  
!geh:print *, 'GeneralFluxDerivative'
  option%iflag = -2
  call GeneralFlux(gen_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up, &
                   thermal_conductivity_up, &
                   gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn, &
                   thermal_conductivity_dn, &
                   area,dist,general_parameter, &
                   option,v_darcy,res,Janal_up,Janal_dn,&
                   general_analytical_derivatives,PETSC_FALSE)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_auxvar_up(idof),global_auxvar_up, &
                     material_auxvar_up, &
                     thermal_conductivity_up, &
                     gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res_pert,Jdummy,Jdummy, &
                     PETSC_FALSE,PETSC_FALSE)
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_up(idof)%pert
!geh:print *, 'up: ', irow, idof, Jup(irow,idof), gen_auxvar_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                     material_auxvar_up, &
                     thermal_conductivity_up, &
                     gen_auxvar_dn(idof),global_auxvar_dn, &
                     material_auxvar_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res_pert,Jdummy,Jdummy, &
                     PETSC_FALSE,PETSC_FALSE)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_dn(idof)%pert
!geh:print *, 'dn: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_isothermal) then
    Jup(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jup(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
    Jdn(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jup(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jup(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
    Jdn(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,20es24.15)') 'flux deriv:', Jup, Jdn
  endif
#endif
  
end subroutine GeneralFluxDerivative

! ************************************************************************** !

subroutine GeneralBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   gen_auxvar_up, &
                                   global_auxvar_up, &
                                   gen_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   thermal_conductivity_dn, &
                                   area,dist,general_parameter, &
                                   option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module 
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow
  PetscReal :: Jdum(option%nflowdof,option%nflowdof)

  Jdn = 0.d0
!geh:print *, 'GeneralBCFluxDerivative'

  option%iflag = -2
  call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     gen_auxvar_up,global_auxvar_up, &
                     gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res,Jdum, &
                     PETSC_FALSE,PETSC_FALSE)                     
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                       gen_auxvar_up,global_auxvar_up, &
                       gen_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area,dist,general_parameter, &
                       option,v_darcy,res_pert,Jdum, &
                       PETSC_FALSE,PETSC_FALSE)   
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_dn(idof)%pert
!print *, 'bc: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_isothermal) then
    Jdn(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jdn(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,10es24.15)') 'bc flux deriv:', Jdn
  endif
#endif
  
end subroutine GeneralBCFluxDerivative

! ************************************************************************** !

subroutine GeneralSrcSinkDerivative(option,qsrc,flow_src_sink_type, &
                                    gen_auxvars,global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(general_auxvar_type) :: gen_auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow
  PetscReal :: Jdum(option%nflowdof,option%nflowdof)  

  option%iflag = -3
  call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                      gen_auxvars(ZERO_INTEGER),global_auxvar,dummy_real, &
                      scale,res,Jdum,PETSC_FALSE,PETSC_FALSE)
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                        gen_auxvars(idof),global_auxvar,dummy_real, &
                        scale,res_pert,Jdum,PETSC_FALSE,PETSC_FALSE)            
    do irow = 1, option%nflowdof
      Jac(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvars(idof)%pert
    enddo !irow
  enddo ! idof
  
  if (general_isothermal) then
    Jac(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jac(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jac(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jac(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,20es24.15)') 'src/sink deriv:', Jac
  endif
#endif

end subroutine GeneralSrcSinkDerivative

! ************************************************************************** !

function GeneralAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn,dden_up,dden_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/14
  ! 

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)
  PetscReal :: dden_up, dden_dn

  PetscReal :: GeneralAverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0
  if (iphase == LIQUID_PHASE) then
    if (istate_up == GAS_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0
    else if (istate_dn == GAS_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == LIQUID_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0      
    else if (istate_dn == LIQUID_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0      
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0      
    endif
  endif

end function GeneralAverageDensity

end module General_Common_module
