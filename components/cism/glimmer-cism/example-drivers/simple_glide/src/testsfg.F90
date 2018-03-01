module testsFG
!   Copyright (C) 2005-2007 Ed Bueler
!  
!   This file is part of PISM.
!  
!   PISM is free software; you can redistribute it and/or modify it under the
!   terms of the GNU General Public License as published by the Free Software
!   Foundation; either version 2 of the License, or (at your option) any later
!   version.
!  
!   PISM is distributed in the hope that it will be useful, but WITHOUT ANY
!   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
!   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
!   details.
!  
!   You should have received a copy of the GNU General Public License
!   along with PISM; if not, write to the Free Software
!   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TESTSFG is a Fortran 90 implementation of two exact solutions for a
! thermocoupled ice sheet.  Reference:
!
!    E. Bueler, J. Brown, and C. Lingle (2007).  "Exact solutions to the 
!    thermomechanically coupled shallow ice approximation: effective tools
!    for verification", J. Glaciol., J. Glaciol., vol. 53 no. 182, 499--516.
!
! ELB 3/29/05; 7/27/07; 7/29/08
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use glimmer_global, only : dp
   public :: testF, testG
   private :: bothexact, p3, p4

   ! DOUBLE PRECISION DESIRABLE:
   integer, parameter, public :: kind=dp

   real(kind), parameter, public :: SperA = 31556926.0  ! 365.2422 days
   real(kind), parameter, public :: g=9.81       ! m/s^2; accel of gravity
   real(kind), parameter, public :: Rgas=8.314   ! J/(mol K)

   ! ice properties; parameters which appear in constitutive relation
   real(kind), parameter, public :: rho=910.0    ! kg/m^3; density
   real(kind), parameter, public :: k=2.1        ! J/m K s; thermal conductivity
   real(kind), parameter, public :: cpheat=2009.0! J/kg K; specific heat capacity
   real(kind), parameter, public :: n=3          ! Glen exponent
   ! next two are EISMINT II values; Paterson-Budd for T<263
   real(kind), parameter, public :: A=3.615e-13  ! Pa^-3 s^-1
   real(kind), parameter, public :: Q=6.0e4      ! J/mol

   ! EISMINT II temperature boundary condition (Experiment F):
   real(kind), parameter, public :: Ggeo=.042    ! J/m^2 s; geo. heat flux
   real(kind), parameter, public :: ST=1.67e-5   ! K m^-1
   real(kind), parameter, public :: Tmin=223.15  ! K

   ! parameters describing extent of sheet
   real(kind), parameter, public :: H0=3000.0    ! m
   real(kind), parameter, public :: L=750000.0   ! m

   ! period and magnitude of perturbation; inactive in Test F:
   real(kind), parameter, public :: Tp=2000.0*SperA  ! s
   real(kind), parameter, public :: Cp=200.0     ! m

contains

   subroutine testF(r,z,H,T,U,w,Sig,M,Sigc)
      real(kind), intent(in) :: r, z
      real(kind), intent(out) :: H, T, U, w, Sig, M, Sigc
      call bothexact(0.0_kind,r,z,0.0_kind,H,T,U,w,Sig,M,Sigc)
   end subroutine testF

   subroutine testG(t,r,z,H,TT,U,w,Sig,M,Sigc)
      real(kind), intent(in) :: t, r, z
      real(kind), intent(out) :: H, TT, U, w, Sig, M, Sigc
      call bothexact(t,r,z,Cp,H,TT,U,w,Sig,M,Sigc)
   end subroutine testG

   subroutine bothexact(t,r,z,Cp,H,TT,U,w,Sig,M,Sigc)
      real(kind), intent(in) :: t, r, z, Cp
      real(kind), intent(out) :: H, TT, U, w, Sig, M, Sigc

      real(kind), parameter :: pi = 3.14159265358979
      real(kind), parameter :: Kcond=k/(rho*cpheat)  ! constant in temp eqn

      ! declare all temporary quantities real(kind); computed in blocks below
      real(kind) pow, Hconst, s, lamhat, f, goft, Ts, nusqrt, nu
      real(kind) lamhatr, fr, Hr, mu, I3, surfArr, Uconst, omega
      real(kind) Sigmu, lamhatrr, frr, Hrr, Tsr, nur, mur, phi, gamma, I4
      real(kind) I4H, divQ, Ht, nut, dTt, Tr, Tz, Tzz

      if (r<=0 .or. r>=L) then
         print *,'code and derivation assume 0<r<L'
         stop
      end if

      ! compute H from analytical steady state Hs (Test D) plus perturbation
      pow = n/(2*n+2)
      Hconst = H0/((1-1/n)**pow)
      s = r/L
      lamhat = (1+1/n)*s - (1/n) + (1-s)**(1+1/n) - s**(1+1/n);
      if (r>0.3*L .and. r<0.9*L) then
         f = ( cos(pi*(r-0.6*L)/(0.6*L)) )**2
      else
         f = 0.0
      end if
      goft = Cp*sin(2.0*pi*t/Tp)
      H = Hconst*(lamhat)**pow + goft*f

      ! compute TT = temperature
      Ts = Tmin+ST*r
      nusqrt = sqrt( 1 + (4.0*H*Ggeo)/(k*Ts) )
      nu = ( k*Ts/(2.0*Ggeo) )*( 1 + nusqrt )
      TT = Ts * (nu+H) / (nu+z)

      ! compute surface slope and horizontal velocity
      lamhatr = ((1+1/n)/L)*( 1 - (1-s)**(1/n) - s**(1/n) )
      if (r>0.3*L .and. r<0.9*L) then
         fr = -(pi/(0.6*L)) * sin(2.0*pi*(r-0.6*L)/(0.6*L))
      else
         fr = 0.0
      end if
      Hr = Hconst * pow * lamhat**(pow-1) * lamhatr + goft*fr   ! chain rule
      if (Hr>0) then
         print *,'code and derivation assume H_r negative for all 0<r<L'
         stop
      end if
      mu = Q/(Rgas*Ts*(nu+H))
      I3 = p3(mu*H) * exp(mu*H) - p3(mu*(H-z)) * exp(mu*(H-z))
      surfArr = exp(-Q/(Rgas*Ts))
      Uconst = 2.0 * (rho*g)**n * A
      omega = Uconst * (-Hr)**n * surfArr * mu**(-n-1)
      U = omega * I3

      ! compute strain heating
      Sigmu = -(Q*(nu+z)) / (Rgas*Ts*(nu+H))
      Sig = (Uconst*g/cpheat) * exp(Sigmu) * ( abs(Hr)*(H-z) )**(n+1)

      ! compute vertical velocity
      lamhatrr = ((1+1/n) / (n*L*L)) * ( (1-s)**((1/n)-1) - s**((1/n)-1) )
      if (r>0.3*L .and. r<0.9*L) then
         frr = -(2.0*pi*pi/(0.36*L*L)) * cos(2.0*pi*(r-0.6*L)/(0.6*L))
      else
         frr = 0.0
      end if
      Hrr = Hconst*pow*(pow-1)*(lamhat)**(pow-2) * lamhatr**2  + &
         Hconst*pow*(lamhat)**(pow-1)*lamhatrr + goft*frr
      Tsr = ST
      nur = (k*Tsr/(2.0*Ggeo)) * (1 + nusqrt) + &
          (1/Ts) * (Hr*Ts-H*Tsr) / nusqrt
      mur = (-Q/(Rgas*Ts*Ts*(nu+H)**2)) * (Tsr*(nu+H)+Ts*(nur+Hr))
      phi = 1/r + n*Hrr/Hr + Q*Tsr/(Rgas*Ts*Ts) - (n+1)*mur/mu   ! division by r
      gamma = mu**n * exp(mu*H) * (mur*H+mu*Hr) * H**n
      I4 = p4(mu*H) * exp(mu*H) - p4(mu*(H-z)) * exp(mu*(H-z))
      w = omega * ((mur/mu - phi)*I4/mu + (phi*(H-z)+Hr)*I3 - gamma*z)

      ! compute compensatory accumulation M
      I4H = p4(mu*H) * exp(mu*H) - 24
      divQ = - omega * (mur/mu - phi) * I4H / mu + omega * gamma * H
      Ht = (Cp*2.0*pi/Tp) * cos(2.0*pi*t/Tp) * f
      M = Ht + divQ

      ! compute compensatory heating
      nut = Ht/nusqrt
      dTt = Ts * ((nut+Ht)*(nu+z)-(nu+H)*nut) * (nu+z)**(-2)
      Tr = Tsr*(nu+H)/(nu+z) + Ts * ((nur+Hr)*(nu+z)-(nu+H)*nur) * (nu+z)**(-2)
      Tz = -Ts * (nu+H) * (nu+z)**(-2)
      Tzz = 2.0 * Ts * (nu+H) * (nu+z)**(-3)
      Sigc = dTt + U*Tr + w*Tz - Kcond*Tzz - Sig
   end subroutine bothexact

   function p3(x)
      real(kind), intent(in) :: x
      !real(kind)             :: p3
      ! p_3=x^3-3*x^2+6*x-6, using Horner's
      p3 = -6 + x*(6 + x*(-3 + x))
   end function p3

   function p4(x)
      real(kind), intent(in) :: x
      !real(kind)             :: p4
      ! p_4=x^4-4*x^3+12*x^2-24*x+24, using Horner's
      p4 = 24 + x*(-24 + x*(12 + x*(-4 + x)))
   end function p4
   
   subroutine model_exact(t,r,z,Hh,H0,TT,U,w,Sig,M,Sigc)
      real(kind), intent(in) :: t, r, z, Hh, H0 
      real(kind), intent(inout) :: TT
      real(kind), intent(out) :: U, w, Sig, M, Sigc

      real(kind), parameter :: pi = 3.14159265358979
      real(kind), parameter :: Kcond=k/(rho*cpheat)  ! constant in temp eqn

      ! declare all temporary quantities real(kind); computed in blocks below
      real(kind) pow, Hconst, s, lamhat, f, goft, Ts, nusqrt, nu
      real(kind) lamhatr, fr, Hr, mu, I3, surfArr, Uconst, omega
      real(kind) Sigmu, lamhatrr, frr, Hrr, Tsr, nur, mur, phi, gamma, I4
      real(kind) I4H, divQ, Ht, nut, dTt, Tr, Tz, Tzz
      real Cp, H

      if (r<=0 .or. r>=L) then
         print *,'code and derivation assume 0<r<L'
         stop
      end if
      H = Hh
      !Don't need the perturbation 
      Cp = 0.0
      ! compute H from analytical steady state Hs (Test D) plus perturbation
      pow = n/(2*n+2)
      Hconst = H0/((1-1/n)**pow)
      s = r/L
      lamhat = (1+1/n)*s - (1/n) + (1-s)**(1+1/n) - s**(1+1/n);
      if (r>0.3*L .and. r<0.9*L) then
         f = ( cos(pi*(r-0.6*L)/(0.6*L)) )**2
      else
         f = 0.0
      end if
      goft = Cp*sin(2.0*pi*t/Tp)
      !H = Hconst*(lamhat)**pow + goft*f
      if (H .gt. Hconst*(lamhat)**pow + goft*f) then
          H = Hconst*(lamhat)**pow + goft*f
      end if
      ! compute TT = temperature
      Ts = Tmin+ST*r
      nusqrt = sqrt( 1 + (4.0*H*Ggeo)/(k*Ts) )
      nu = ( k*Ts/(2.0*Ggeo) )*( 1 + nusqrt )
      if(TT .eq. 0.0) then
          TT = Ts * (nu+H) / (nu+z)
      end if

      ! compute surface slope and horizontal velocity
      lamhatr = ((1+1/n)/L)*( 1 - (1-s)**(1/n) - s**(1/n) )
      if (r>0.3*L .and. r<0.9*L) then
         fr = -(pi/(0.6*L)) * sin(2.0*pi*(r-0.6*L)/(0.6*L))
      else
         fr = 0.0
      end if
      Hr = Hconst * pow * lamhat**(pow-1) * lamhatr + goft*fr   ! chain rule
      if (Hr>0) then
         print *,'code and derivation assume H_r negative for all 0<r<L'
         stop
      end if
      mu = Q/(Rgas*Ts*(nu+H))
      I3 = p3(mu*H) * exp(mu*H) - p3(mu*(H-z)) * exp(mu*(H-z))
      surfArr = exp(-Q/(Rgas*Ts))
      Uconst = 2.0 * (rho*g)**n * A
      omega = Uconst * (-Hr)**n * surfArr * mu**(-n-1)
      U = omega * I3

      ! compute strain heating
      Sigmu = -(Q*(nu+z)) / (Rgas*Ts*(nu+H))
      Sig = (Uconst*g/cpheat) * exp(Sigmu) * ( abs(Hr)*(H-z) )**(n+1)

      ! compute vertical velocity
      lamhatrr = ((1+1/n) / (n*L*L)) * ( (1-s)**((1/n)-1) - s**((1/n)-1) )
      if (r>0.3*L .and. r<0.9*L) then
         frr = -(2.0*pi*pi/(0.36*L*L)) * cos(2.0*pi*(r-0.6*L)/(0.6*L))
      else
         frr = 0.0
      end if
      Hrr = Hconst*pow*(pow-1)*(lamhat)**(pow-2) * lamhatr**2  + &
         Hconst*pow*(lamhat)**(pow-1)*lamhatrr + goft*frr
      Tsr = ST
      nur = (k*Tsr/(2.0*Ggeo)) * (1 + nusqrt) + &
          (1/Ts) * (Hr*Ts-H*Tsr) / nusqrt
      mur = (-Q/(Rgas*Ts*Ts*(nu+H)**2)) * (Tsr*(nu+H)+Ts*(nur+Hr))
      phi = 1/r + n*Hrr/Hr + Q*Tsr/(Rgas*Ts*Ts) - (n+1)*mur/mu   ! division by r
      gamma = mu**n * exp(mu*H) * (mur*H+mu*Hr) * H**n
      I4 = p4(mu*H) * exp(mu*H) - p4(mu*(H-z)) * exp(mu*(H-z))
      w = omega * ((mur/mu - phi)*I4/mu + (phi*(H-z)+Hr)*I3 - gamma*z)

      ! compute compensatory accumulation M
      I4H = p4(mu*H) * exp(mu*H) - 24
      divQ = - omega * (mur/mu - phi) * I4H / mu + omega * gamma * H
      Ht = (Cp*2.0*pi/Tp) * cos(2.0*pi*t/Tp) * f
      M = Ht + divQ

      ! compute compensatory heating
      nut = Ht/nusqrt
      dTt = Ts * ((nut+Ht)*(nu+z)-(nu+H)*nut) * (nu+z)**(-2)
      Tr = Tsr*(nu+H)/(nu+z) + Ts * ((nur+Hr)*(nu+z)-(nu+H)*nur) * (nu+z)**(-2)
      Tz = -Ts * (nu+H) * (nu+z)**(-2)
      Tzz = 2.0 * Ts * (nu+H) * (nu+z)**(-3)
      Sigc = dTt + U*Tr + w*Tz - Kcond*Tzz - Sig
   end subroutine model_exact
end module testsFG
