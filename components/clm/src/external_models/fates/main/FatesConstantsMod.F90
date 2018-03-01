module FatesConstantsMod
  ! This module is used to define global _immutable_ data. Everything in
  ! this module must have the parameter attribute.

  implicit none

  public

  ! kinds
  integer, parameter :: fates_r8 = selected_real_kind(12) ! 8 byte real

  ! string lengths
  integer, parameter :: fates_avg_flag_length = 3
  integer, parameter :: fates_short_string_length = 32
  integer, parameter :: fates_long_string_length = 199

  ! Unset and various other 'special' values
  integer, parameter :: fates_unset_int = -9999

  
  ! Integer equivalent of true  (in case some compilers dont auto convert)
  integer, parameter :: itrue = 1

  ! Integer equivalent of false (in case come compilers dont auto convert)
  integer, parameter :: ifalse = 0

  
  ! Unit conversion constants:

  ! Conversion factor umols of Carbon -> kg of Carbon (1 mol = 12g)
  ! We do not use umolC_per_kg because it is a non-terminating decimal
  real(fates_r8), parameter :: umolC_to_kgC = 12.0E-9_fates_r8
  
  ! Conversion factor: grams per kilograms
  real(fates_r8), parameter :: g_per_kg = 1000.0_fates_r8
  
  ! Conversion factor: miligrams per grams
  real(fates_r8), parameter :: mg_per_g = 1000.0_fates_r8

  ! Conversion factor: micromoles per milimole
  real(fates_r8), parameter :: umol_per_mmol = 1000.0_fates_r8

  ! Conversion factor: milimoles per mole
  real(fates_r8), parameter :: mmol_per_mol = 1000.0_fates_r8
  
  ! Conversion factor: micromoles per mole
  real(fates_r8), parameter :: umol_per_mol = 1.0E6_fates_r8

  ! Conversion factor: m2 per ha
  real(fates_r8), parameter :: m2_per_ha = 1.0e4_fates_r8

  ! Conversion factor :: ha per m2
  real(fates_r8), parameter :: ha_per_m2 = 1.0e-4_fates_r8

  ! Conversion: seconds per minute
  real(fates_r8), parameter :: sec_per_min = 60.0_fates_r8

  ! Conversion: seconds per day
  real(fates_r8), parameter :: sec_per_day = 86400.0_fates_r8

  ! Conversion: days per second
  real(fates_r8), parameter :: days_per_sec = 1.0_fates_r8/86400.0_fates_r8

  ! Conversion: days per year. assume HLM uses 365 day calendar.  
  ! If we need to link to 365.25-day-calendared HLM, rewire to pass through interface
  real(fates_r8), parameter :: days_per_year = 365.00_fates_r8
  
  ! Conversion: years per day. assume HLM uses 365 day calendar.  
  ! If we need to link to 365.25-day-calendared HLM, rewire to pass through interface
  real(fates_r8), parameter :: years_per_day = 1.0_fates_r8/365.00_fates_r8
  
  ! Physical constants

  ! universal gas constant [J/K/kmol]
  real(fates_r8), parameter :: rgas_J_K_kmol          = 8314.4598_fates_r8

  ! freezing point of water at 1 atm (K)
  real(fates_r8), parameter :: t_water_freeze_k_1atm   = 273.15_fates_r8     

  ! freezing point of water at triple point (K)
  real(fates_r8), parameter :: t_water_freeze_k_triple = 273.16_fates_r8      

  ! Density of fresh liquid water (kg/m3)
  real(fates_r8), parameter :: dens_fresh_liquid_water = 1.0E3_fates_r8

  ! Gravity constant on earth [m/s]
  real(fates_r8), parameter :: grav_earth = 9.8_fates_r8


  ! For numerical inquiry
  real(fates_r8), parameter :: fates_huge = huge(g_per_kg)

  real(fates_r8), parameter :: fates_tiny = tiny(g_per_kg)

  ! Geometric Constants
  
  ! PI
  real(fates_r8), parameter :: pi_const = 3.14159265359_fates_r8

end module FatesConstantsMod
