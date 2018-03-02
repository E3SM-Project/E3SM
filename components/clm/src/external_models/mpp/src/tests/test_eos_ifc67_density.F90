program test_eos_ifc67

  use EOSWaterMod, only : Density
  use EOSWaterMod, only : DENSITY_IFC67

  implicit none

#include <petsc/finclude/petsc.h>
  
  PetscReal            :: den
  PetscReal            :: dden_dp
  PetscReal            :: dden_dT
  PetscReal            :: val, val_bmark, err_thres
  character(len=10)    :: val_name
  !
  PetscReal, parameter :: p                         = 120000.d0                ! [Pa]
  PetscReal, parameter :: t_K                       = 300.d0                   ! [K]
  !
  PetscReal, parameter :: den_bmark                 = 55.323696656461536d0     ! [kmol m^{-3}]
  PetscReal, parameter :: dden_dp_bmark             = 2.4854904480147891d-008  ! [kmol m^{-3} Pa^{-1}]
  PetscReal, parameter :: dden_dT_bmark             = -1.5298638598102345d-002 ! [kmol m^{-3} K^{-1}]
  !
  PetscReal, parameter :: den_abs_err_threshold     = 1.d-11
  PetscReal, parameter :: dden_dp_abs_err_threshold = 1.d-16
  PetscReal, parameter :: dden_dT_abs_err_threshold = 1.d-15

  call Density(p, t_K, DENSITY_IFC67, den, dden_dp, dden_dT)

  ! Check density
  val       = den
  val_bmark = den_bmark
  err_thres = den_abs_err_threshold
  val_name  = 'den'
  if ( abs(val - val_bmark) > err_thres) then
     write(*,*) val_name // ' error greater than ',err_thres
     stop 1
  endif
  
  ! Check derivative of density w.r.t. pressure
  val       = dden_dp
  val_bmark = dden_dp_bmark
  err_thres = dden_dp_abs_err_threshold
  val_name  = 'dden_dp'
  if ( abs(val - val_bmark) > err_thres) then
     write(*,*) val_name // ' error greater than ',err_thres
     stop 1
  endif
  
  ! Check derivative of density w.r.t. temperature
  val       = dden_dT
  val_bmark = dden_dT_bmark
  err_thres = dden_dT_abs_err_threshold
  val_name  = 'dden_dT'
  if ( abs(val - val_bmark) > err_thres) then
     write(*,*) val_name // ' error greater than ',err_thres
     stop 1
  endif
  
  
end program test_eos_ifc67
