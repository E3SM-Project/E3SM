module Gas_EOS_module
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  public
 
contains

! ************************************************************************** !

subroutine ideal_gaseos_noderiv(p,tc,d,h,u)
    
  PetscReal,intent(in) :: p  ! [Pa]
  PetscReal,intent(in) :: tc ! [C]
  PetscReal, intent(out) :: d ! [kmol/m^3]
  PetscReal, intent(out) :: h ! [J/kmol]
  PetscReal, intent(out) :: u ! [J/kmol]
   
  PetscReal, parameter :: Rg=8.31415 
  ! Cpg units: J/mol-K
  PetscReal, parameter :: Cv_air = 20.85 ! J/(K mol) heat capacity wiki
  PetscReal  tk
    
  tk = tc + 273.15
  d = p / tk / Rg * 1.d-3 ! mol/m^3 -> kmol/m^3
  h = Cv_air * tk * 1.d3  ! J/mol -> J/kmol
  u = (Cv_air - Rg) * tk * 1.d3 ! J/mol -> J/kmol
   
end subroutine ideal_gaseos_noderiv

! ************************************************************************** !

subroutine ideal_gaseos(p,tc,d,d_p,d_t,h,h_p,h_t,u,u_p,u_t)
    
  PetscReal,intent(in) :: p  ! [Pa]
  PetscReal,intent(in) :: tc ! [C]
  PetscReal, intent(out) :: d ! [kmol/m^3]
  PetscReal, intent(out) :: h ! [J/kmol]
  PetscReal, intent(out) :: u ! [J/kmol]
  PetscReal, intent(out) :: d_p,d_t,h_p,h_t,u_p,u_t

  PetscReal, parameter :: Rg=8.31415 
  ! Cpg units: J/mol-K
  PetscReal, parameter :: Cv_air = 20.85 ! J/(K mol) heat capacity wiki
  PetscReal  tk

  tk = tc + 273.15
  d = p / tk / Rg * 1.d-3 ! mol/m^3 -> kmol/m^3
  h = Cv_air * tk * 1.d3  ! J/mol -> J/kmol
  u = (Cv_air - Rg) * tk * 1.d3 ! J/mol -> J/kmol

  d_p =  d / p
  d_t = -d / tk
  h_p = 0.d0
  h_t = Cv_air * 1.d3
  u_p = 0.d0
  u_t = (Cv_air - Rg) * 1.d3

end subroutine ideal_gaseos

! ************************************************************************** !

subroutine visgas_noderiv(t,p_air,p_gas,d_air,visg)
  ! 
  ! REFERENCES
  ! THIS ROUTINE IS LARGELY ADAPTED FROM THE TOUGH CODE.
  ! this routine computes the viscosity of vapor-air mixtures.
  ! it uses a modified version of a formulation based on kinetic
  ! gas theory, as given by j.o. hirschfelder, c.f. curtiss, and
  ! r.b. bird, molecular theory of gases and liquids, john wiley
  ! & sons, 1954, pp. 528-530.
  ! the modification made to the hirschfelder et al. expressions is
  ! that for vapor viscosity accurate (empirical) values are used,
  ! rather than the first order expression of kinetic theory.
  ! the formulation matches experimental data on viscosities of
  ! vapor-air mixtures in the temperature range from 100 to 150
  ! deg. c, for all compositions, to better than 4%.
  ! 
  PetscReal, intent(in) :: t     ! [C]
  PetscReal, intent(in) :: p_air ! [Pa]
  PetscReal, intent(in) :: p_gas ! [Pa]
  PetscReal, intent(in) :: d_air ! [kmol/m^3]
  PetscReal, intent(out) :: visg ! [Pa-s]

  PetscReal ::  fair,fwat,cair,cwat

      data  fair,   fwat,    cair,  cwat &
           /97.d0, 363.d0, 3.617d0, 2.655d0/
 
      PetscReal fmix,cmix,d,xga,xg1,tk,trd1,trd3,ome1,ome3,ard,fmw3,vis1, &
             v1,vs,vis2,vis3,z1,g,h,e,z2,z3


!c======================================================================

      fmix = sqrt (fair*fwat)
      cmix = (cair+cwat)*0.5d0

!      do k = 1,nb
 !       if (iphas(k).eq.2 .or. iphas(k).eq.0) then

          d   = d_air *FMWAIR       
          xga = p_air / p_gas ! for debug, set x constant
          xg1 = 1.D0 - xga
          tk  = t +273.15d0

          trd1 = tk/fair
          trd3 = tk/fmix
          ome1 = (1.188d0-0.051d0*trd1)/trd1
          ome3 = (1.480d0-0.412d0*log(trd3))/trd3
          ard  = 1.095d0/trd3
          fmw3 = 2.d0*FMWAIR*FMWH2O/(FMWAIR+FMWH2O)
          vis1 = 266.93d-7*sqrt(FMWAIR*trd1*fair)/(cair*cair*ome1*trd1)
 
          v1 = .407d0*t +80.4d0
          if (t .le.350.d0) then
            vs = 1.d-7*(v1-d*(1858.d0-5.9d0*t )*1.d-3)
          else
!             if (t .gt.350.d0) 
!cpcl .      vs = 1.d-7*(v1 + 0.353d0*d + 676.5d-6*d**2 + 102.1d-9*d**3)
           vs = 1.d-7*(v1 + (0.353d0 + (676.5d-6 + 102.1d-9*d)*d)*d)
          endif

          vis2 = 10.d0*vs
          vis3 = 266.93d-7*sqrt(fmw3*trd3*fmix)/(cmix*cmix*ome3*trd3)
          z1   = xga*xga/vis1+2.d0*xg1*xga/vis3+xg1*xg1/vis2
          g    = xga*xga*FMWAIR/FMWH2O
          h    = xg1*xg1*FMWH2O/FMWAIR
          e    = (2.d0*xga*xg1*FMWAIR*FMWH2O/fmw3**2)*vis3/(vis1*vis2)
          z2   = 0.6d0*ard*(g/vis1+e+h/vis2)
          z3   = 0.6d0*ard*(g+e*(vis1+vis2)-2.d0*xga*xg1+h)
          visg  = (1.d0+z3)/(z1+z2)*.1d0
           
end subroutine visgas_noderiv

! ************************************************************************** !

subroutine Henry_air_noderiv(p,tc,ps,Henry)
! Calculate Henry Coefficient for N2
! t in K
! Henry have the same unit as p and ps, then make it dimensionless by
! devide it with p

    implicit none
    PetscReal,intent(in) ::  p,tc,ps
    PetscReal,intent(out) ::  Henry

    PetscReal  Tr,tao,tmp,t
    PetscReal, parameter :: a=-9.67578, b=4.72162, c=11.70585
    PetscReal, parameter :: Tcl=647.096 ! H2O critical temp(K) from IAPWS(1995b)

    t=tc+273.15D0
    Tr=t/Tcl
    tao=1.D0-Tr
    tmp= a/Tr + B * tao**0.355/Tr + c * (Tr**(-0.41)) * exp(tao)
    Henry=exp(tmp)*ps

   return 
end subroutine Henry_air_noderiv

! ************************************************************************** !

subroutine Henry_air(p,tc,ps,ps_p,ps_t,Henry,Henry_p,Henry_t)
   implicit none
    PetscReal,intent(in) ::  p,tc,ps,ps_p,ps_t
    PetscReal,intent(out) ::  Henry,Henry_p,Henry_t
! note t/K, p/Pa, Henry/Pa 

    PetscReal  Tr,tao,tmp,t
    PetscReal, parameter :: a=-9.67578, b=4.72162, c=11.70585
    PetscReal, parameter :: Tcl=647.096 ! H2O critical temp from IAPWS(1995b)

    t=tc+273.15D0
    Tr=t/Tcl
    tao=1.D0-Tr
    tmp= a/Tr + b * tao**0.355/Tr + c * (Tr**(-0.41)) * exp(tao)
    Henry=exp(tmp)*ps

    tmp =((-a/Tr+b*(-0.355*tao**(-0.645)-tao**0.355/Tr))/Tr - &
         c*exp(tao)*(tao**(-.41))*(0.41/Tr-1.))/Tcl
    Henry_t=Henry*(tmp +ps_t/ps)
    Henry_p=ps_p*Henry/ps

  
   return 
end subroutine Henry_air

end module Gas_EOS_module
