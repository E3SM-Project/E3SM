! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2015 - D. Swales - Original version
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_CALIPSO_INTERFACE
  USE COSP_KINDS,              ONLY: wp
  USE MOD_LIDAR_SIMULATOR,     ONLY: alpha,beta,gamma
  IMPLICIT NONE
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE calipso_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  type calipso_IN
     integer,pointer ::       &
          Npoints,      & ! Number of gridpoints.
          Ncolumns,     & ! Number of columns.
          Nlevels         ! Number of levels.

     real(wp),dimension(:,:),pointer :: &
          beta_mol,     & ! Molecular backscatter coefficient
          tau_mol         ! Molecular optical depth
     real(wp),dimension(:,:,:),pointer :: &
          betatot,      & ! 
          tautot,       & ! Optical thickess integrated from top
          betatot_ice,  & ! Backscatter coefficient for ice particles
          betatot_liq,  & ! Backscatter coefficient for liquid particles
          tautot_ice,   & ! Total optical thickness of ice
          tautot_liq      ! Total optical thickness of liq
     real(wp),dimension(:,:,:,:),pointer :: &
          taupart
  end type calipso_IN

CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_calipso_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_calipso_init() 
    
    ! Polynomial coefficients (Alpha, Beta, Gamma) which allow to compute the 
    ! ATBperpendicular as a function of the ATB for ice or liquid cloud particles 
    ! derived from CALIPSO-GOCCP observations at 120m vertical grid 
    ! (Cesana and Chepfer, JGR, 2013).
    !
    ! Relationship between ATBice and ATBperp,ice for ice particles:
    !                ATBperp,ice = Alpha*ATBice 
    ! Relationship between ATBice and ATBperp,ice for liquid particles:
    !          ATBperp,ice = Beta*ATBice^2 + Gamma*ATBice
    Alpha = 0.2904_wp
    Beta  = 0.4099_wp
    Gamma = 0.009_wp    
    
  end subroutine cosp_calipso_init

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !	END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_CALIPSO_INTERFACE
