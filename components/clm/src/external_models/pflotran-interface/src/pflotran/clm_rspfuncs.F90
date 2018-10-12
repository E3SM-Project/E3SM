module CLM_RspFuncs_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

! constants  
  PetscReal, parameter, public :: rpi = 3.14159265358979323846d0

! temperature response function
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_CLMCN = 1
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_Q10  = 2
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_DLEM = 3

! moisture response function
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_CLMCN = 1
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_DLEM = 2

! pH response function
  PetscInt, parameter, public :: PH_RESPONSE_FUNCTION_CENTURY = 1 
  PetscInt, parameter, public :: PH_RESPONSE_FUNCTION_DLEM = 2

! molecular weight
  PetscReal, parameter, public :: N_molecular_weight = 14.0067d0
  PetscReal, parameter, public :: C_molecular_weight = 12.0110d0
  PetscReal, parameter, public :: CN_ratio_mass_to_mol = 1.166156023644992d0
  ! A NOTE here: when coupled with CLM-CN, make sure that the above constants ARE consistent with CLM
  ! otherwise may cause some tiny but detectable mass-balance errors due to unit conversion.

  public :: GetTemperatureResponse, &
            GetMoistureResponse, &
            GetpHResponse, &
            FuncMonod

contains
! ************************************************************************** !
! temperature response function 

function GetTemperatureResponse(tc, itype, Q10)

  implicit none
  
  PetscInt  :: itype
  PetscReal :: Ft, tc, Q10, tk

  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal :: GetTemperatureResponse
  PetscReal, parameter :: Frz_Q10 = 2.0d0

  select case(itype)
!     CLM4.5 temperature response function
      case(TEMPERATURE_RESPONSE_FUNCTION_Q10)
        if(tc>0.d0) then
          Ft = Q10 ** ((tc - 25.0d0) / 10.0d0)
        else
          Ft = (Q10**(-25.0d0/10.0d0))*(Frz_Q10**((tc/10.0d0)))
        endif

!     CLM-CN
!     Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))
      case(TEMPERATURE_RESPONSE_FUNCTION_CLMCN)
  
          tk = tc + 273.15d0

          if(tk > 227.15d0) then
              Ft = exp(308.56d0*(one_over_71_02 - 1.d0/(tk - 227.13d0)))
          else
              Ft = 0.d0
          endif
!     DLEM temperature response function for methane oxidation
!     Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 12
      case(TEMPERATURE_RESPONSE_FUNCTION_DLEM)
          if(tc < -5.0d0) then
            Ft = 0.0d0
          elseif(tc >= 30.0d0) then
            Ft = 1.0d0
          else
            Ft = Q10 ** ((tc - 30.0d0) / 10.0d0)
          endif
       case default
            Ft = 1.0d0
  end select
  GetTemperatureResponse = Ft 

end function GetTemperatureResponse
  
! ************************************************************************** !

Function GetMoistureResponse(theta, ghosted_id, itype)


#ifdef CLM_PFLOTRAN
#include "petsc/finclude/petscvec.h"
  use petscvec
  use clm_pflotran_interface_data
#endif
  
  implicit none
  
  PetscReal :: F_theta
  ! theta IS soil VWC: 0 - porosity (not adjusted)
  PetscReal :: theta
  PetscReal :: GetMoistureResponse
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0

  PetscInt :: ghosted_id, itype
  PetscReal :: maxpsi, psi, lsat
  PetscReal, parameter :: minpsi = -10.0d6    ! Pa

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: sucsat_pf_loc(:)    !
  PetscScalar, pointer :: watfc_pf_loc(:)     !
  PetscScalar, pointer :: porosity_pf_loc(:)  !
  PetscScalar, pointer :: bd_dry_pf_loc(:)    !
  PetscScalar, pointer :: bsw_pf_loc(:)    !
  PetscReal :: thetar, thetas, se
#endif

  PetscErrorCode :: ierr

#ifdef CLM_PFLOTRAN

  select case(itype)
!   CLM-CN
    case(MOISTURE_RESPONSE_FUNCTION_CLMCN)
      call VecGetArrayReadF90(clm_pf_idata%sucsat_pfs, sucsat_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(clm_pf_idata%bulkdensity_dry_pfs, bd_dry_pf_loc, ierr)   ! 'bd' (kg/m3)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(clm_pf_idata%bsw_pfs, bsw_pf_loc, ierr)
      CHKERRQ(ierr)
      ! sucsat [mm of H20] from CLM is the suction (positive) at water saturated (called air-entry pressure)
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      maxpsi = sucsat_pf_loc(ghosted_id) * (-EARTH_GRAVITY)                         ! mmH2O --> -Pa
      lsat = theta/min(1.d0, 1.d0-min(0.9999d0,bd_dry_pf_loc(ghosted_id)/2.70d3))     ! bd = (1._r8-dry_porosity)*2.7d3

      ! soil matric potential by Clapp-Hornburger method (this is the default used by CLM)
      psi = sucsat_pf_loc(ghosted_id) * (-EARTH_GRAVITY) * (lsat**(-bsw_pf_loc(ghosted_id)))  ! mmH2O --> -Pa
      psi = min(psi, maxpsi)
      if(psi > minpsi) then
        F_theta = log(minpsi/psi)/log(minpsi/maxpsi)
        ! very wet soil (close to saturated)
        if(psi>(maxpsi-1.d02)) then
          F_theta = F_theta*0.10d0   ! 0.10 is an arbitrary value here (but NOT totaly shut off decomp)
        endif
      else
        F_theta = 0.0d0
      endif

      call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pfs, sucsat_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(clm_pf_idata%bulkdensity_dry_pfs, bd_dry_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(clm_pf_idata%bsw_pfs, bsw_pf_loc, ierr)
      CHKERRQ(ierr)

!     DLEM
!     Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 13
    case(MOISTURE_RESPONSE_FUNCTION_DLEM) 
      call VecGetArrayReadF90(clm_pf_idata%effporosity_pfs, porosity_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(clm_pf_idata%watfc_pfs, watfc_pf_loc, ierr)
      CHKERRQ(ierr)
      thetas = porosity_pf_loc(ghosted_id)
      thetar = watfc_pf_loc(ghosted_id)
      if(theta >= thetas) then
        F_theta = 1.0d0
      elseif (theta <= thetar) then
        F_theta = 0.0d0
      else
        se = (theta - thetar)/(thetas - thetar)
        F_theta = 1.0 - se * se * 0.368 * exp(se)

        if(F_theta < 0.0d0) then
           F_theta = 0.0d0
        endif

        if(F_theta > 1.0d0) then
           F_theta = 1.0d0
        endif
      endif
      call VecRestoreArrayReadF90(clm_pf_idata%effporosity_pfs, porosity_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(clm_pf_idata%watfc_pfs, watfc_pf_loc, ierr)
      CHKERRQ(ierr)
    case default
        F_theta = 1.0d0
  end select
#else

  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = 1.0d0 !log(theta_min/max(theta_min,theta)) * one_over_log_theta_min 
  
#endif

  GetMoistureResponse = F_theta

end function GetMoistureResponse

! ************************************************************************** !

Function GetpHResponse(pH, itype)

  PetscReal :: f_ph, pH, GetpHResponse
  PetscInt :: itype

! ph function from Parton et al., (2001, 1996)
!  k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi
#ifdef CLM_PFLOTRAN
  select case(itype)
      case(PH_RESPONSE_FUNCTION_CENTURY)
          f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + pH))/rpi

!     DLEM temperature response function for methane oxidation
! Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 12
      case(PH_RESPONSE_FUNCTION_DLEM)
          if(pH <= 4.0d0 .or. pH >= 10.0d0) then
            f_ph = 0.0d0
          elseif(pH < 7.0d0) then
            f_ph = 1.02d0 /(1.0d0 + 1.0d6 * exp(-2.5d0 * pH))
          else
            f_ph = 1.02d0 /(1.0d0 + 1.0d6 * exp(-2.5d0 * (14.0d0 - pH)))
          endif
      case default
          f_ph = 1.0d0
  end select
#else
  f_ph = 1.0
#endif

  if(f_ph < 0.0d0) then
     f_ph = 0.0d0
  endif

  if(f_ph > 1.0d0) then
     f_ph = 1.0d0
  endif
  GetpHResponse = f_ph

end function GetpHResponse

! ************************************************************************** !
! Monod function
Function funcMonod(conc, conc_halfsat, compute_derivative)

  implicit none

  PetscBool :: compute_derivative
  PetscReal :: conc, conc_halfsat
  PetscReal :: funcMonod

  !----------------------------------------------------------
  if (.not.compute_derivative) then
    funcMonod = conc/(conc+conc_halfsat)
  else
    funcMonod = conc_halfsat/(conc+conc_halfsat)/(conc+conc_halfsat)
  endif

end function funcMonod

! ************************************************************************** !

end module CLM_RspFuncs_module
