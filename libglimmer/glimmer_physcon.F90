
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_physcon.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> Contains physical constants required by the ice model.
module glimmer_physcon

  use glimmer_global, only : dp

#ifdef CCSMCOUPLED 

  !Note: These CESM values are current as of September 2015.
  use shr_const_mod, only: pi=>     SHR_CONST_PI,&       ! 3.14159265358979323846
                           rhoi=>   SHR_CONST_RHOICE,&   ! 0.917e3
                           rhoo=>   SHR_CONST_RHOSW,&    ! 1.026e3
                           rhow=>   SHR_CONST_RHOFW,&    ! 1.000e3
                           rearth=> SHR_CONST_REARTH,&   ! 6.37122e6
                           grav=>   SHR_CONST_G,&        ! 9.80616
                           shci=>   SHR_CONST_CPICE,&    ! 2.11727e3
                           lhci=>   SHR_CONST_LATICE,&   ! 3.337e5
                           trpt=>   SHR_CONST_TKTRIP     ! 273.16
  implicit none  
  save
 
#else   

  implicit none  
  save
                                         
  real(dp),parameter :: pi = 3.14159265358979d0  !< Value of \f$\pi\f$.
  real(dp),parameter :: rhow = 1000.0d0          !< The density of fresh water (kg m<SUP>-3</SUP>)  
  real(dp),parameter :: rearth  = 6.37122d6      ! radius of earth (m)  

  !NOTE: The following are not parameters because the default values can be overridden in the config file.
  !      This may be desirable for test problems such as MISMIP that specify values different from CISM's default values. 
  !      It may also be useful to set these to CESM values, for comparing CESM-coupled to standalone runs.
  real(dp) :: rhoi = 910.d0                      !< The density of ice (kg m<SUP>-3</SUP>)  
  real(dp) :: rhoo = 1028.0d0                    !< The density of the ocean (kg m<SUP>-3</SUP>)
  real(dp) :: grav = 9.81d0                      !< The acceleration due to gravity (m s<SUP>-2</SUP>)
  real(dp) :: shci = 2009.0d0                    !< Specific heat capacity of ice (J kg<SUP>-1</SUP> K<SUP>-1</SUP>) 
  real(dp) :: lhci = 335.0d3                     !< Latent heat of melting of ice (J kg<SUP>-1</SUP>)  
  real(dp) :: trpt = 273.15d0                    !< Triple point of water (K)    
#endif

  real(dp),parameter :: scyr = 31536000.d0       !< Number of seconds in a year of exactly 365 days
  real(dp),parameter :: rhom = 3300.0d0          !< The density of magma(?) (kg m<SUP>-3</SUP>) 
  real(dp),parameter :: rhos = 2600.0d0          !< The density of solid till (kg m$^{-3}$) 
  integer, parameter :: gn = 3                   !< The power dependency of Glen's flow law.
  real(dp),parameter :: actenh = 139.0d3         !< Activation energy in Glen's flow law for \f$T^{*}\geq263\f$K. (J mol<SUP>-1</SUP>)
  real(dp),parameter :: actenl = 60.0d3          !< Activation energy in Glen's flow law for \f$T^{*}<263\f$K. (J mol<SUP>-1</SUP>)  
  real(dp),parameter :: arrmlh = 1.733d3         !< Constant of proportionality in Arrhenius relation
                                                 !< in \texttt{patebudd}, for \f$T^{*}\geq263\f$K.
                                                 !< (Pa<SUP>-3</SUP> s<SUP>-1</SUP>) 
  real(dp),parameter :: arrmll = 3.613d-13       !< Constant of proportionality in Arrhenius relation
                                                 !< in \texttt{patebudd}, for \f$T^{*}<263\f$K.
                                                 !< (Pa<SUP>-3</SUP> s<SUP>-1</SUP>) 
  real(dp),parameter :: gascon = 8.314d0         !< The gas ideal constant \f$R\f$ (J mol<SUP>-1</SUP> K<SUP>-1</SUP>)
  real(dp),parameter :: coni = 2.1d0             !< Thermal conductivity of ice (W m<SUP>-1</SUP> K<SUP>-1</SUP>)
  real(dp),parameter :: pmlt = 9.7456d-8         !< Factor for dependence of melting point on pressure (K Pa<SUP>-1</SUP>)
  real(dp),parameter :: tocnfrz_sfc = -1.92d0    !< Freezing temperature of seawater (deg C) at surface pressure, S = 35 PSU
  real(dp),parameter :: dtocnfrz_dh = -7.53d-4   !< Rate of change of freezing temperature of seawater with depth (deg/m), given S = 35 PSU
                                                 !< These values are from the Ocean Water Freezing Point Calculator,
                                                 !< http://www.csgnetwork.com/h2ofreezecalc.html (25 Nov. 2014)
end module glimmer_physcon
