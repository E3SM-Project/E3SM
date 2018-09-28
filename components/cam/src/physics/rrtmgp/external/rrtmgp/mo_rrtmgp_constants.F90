! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Physical and mathematical constants used in gas optics calculation
!   If the host model in which RRTGMP is embedded has defined these constants elsewhere
!   the model definitions can be used instead by renaming. For example,
! use  mo_model_constants, only k_boltz => boltzman_k, ...
!   where the syntax is local_name => original_name
!   and all the local names need to be defined
!
! "Constants" specific to the earth's atmosphere should also be made consistent with the
!   host model but may be changed in a call to init_constants(), normally at initialization
! -------------------------------------------------------------------------------------------------
module mo_rrtmgp_constants
  use mo_rte_kind, only: wp

  ! -----------------------------------------
  ! Physical constants
  ! Boltzmann constant [J/K] = [(kg m^2)/(K s^2)]
  real(wp), parameter :: k_boltz = 1.380650e-23_wp

  !  molecular weight of water [kg/mol]
  real(wp), parameter :: m_h2o =  0.018016_wp

  ! Avogadro's number [molec/mol]
  real(wp), parameter :: avogad = 6.02214199e23_wp

  ! Universal gas constant [J/(mol K)]
  real(wp), parameter :: R_univ_gconst = 8.314472_wp

  ! -----------------------------------------
  !
  ! Constants specific to the earth's atmosphere -- changeable in init() because they
  !   might be different on e.g. other planets

  ! molecular weight of dry air [kg/mol]
  real(wp), protected :: m_dry = 0.028964_wp

  ! Gravity at Earth's surface [m/s2]
  real(wp), protected :: grav = 9.80665_wp

  ! Specific heat at constant pressure for dry air [J/(K kg)]
  real(wp), protected :: cp_dry = 1004.64_wp

contains
  ! -----------------------------------------
  subroutine init_constants(gravity, mol_weight_dry_air, heat_capacity_dry_air)
    real(wp), optional, intent(in) :: gravity, mol_weight_dry_air, heat_capacity_dry_air

    if(present(gravity))               grav   = gravity
    if(present(mol_weight_dry_air))    m_dry  = mol_weight_dry_air
    if(present(heat_capacity_dry_air)) cp_dry = heat_capacity_dry_air

  end subroutine init_constants
  ! -----------------------------------------
end module mo_rrtmgp_constants
