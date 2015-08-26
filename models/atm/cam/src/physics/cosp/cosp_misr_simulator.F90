! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Nov 2008 - A. Bodas-Salcedo - Initial version
!
!

MODULE MOD_COSP_MISR_SIMULATOR
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------- SUBROUTINE COSP_MISR_SIMULATOR -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_MISR_SIMULATOR(gbx,sgx,y)
  
  ! Arguments
  type(cosp_gridbox),intent(in) :: gbx  ! Gridbox info
  type(cosp_subgrid),intent(in) :: sgx  ! Subgridbox info
  type(cosp_misr),intent(inout) :: y    ! MISR simulator output
  
  ! Local variables 
  integer :: Nlevels,Npoints
  real :: dtau_s(gbx%Npoints, gbx%Nlevels)
  real :: dtau_s_snow(gbx%Npoints, gbx%Nlevels)
  real :: dtau_c(gbx%Npoints, gbx%Nlevels)
  real :: at(gbx%Npoints, gbx%Nlevels)
  real :: frac_out(gbx%Npoints, gbx%Ncolumns, gbx%Nlevels)
  real :: prec_frac(gbx%Npoints, gbx%Ncolumns, gbx%Nlevels)
  integer :: sunlit(gbx%Npoints)
  
  real :: zfull(gbx%Npoints, gbx%Nlevels) !  height (in meters) of full model levels (i.e. midpoints)
                                          !  zfull(npoints,1)    is    top level of model
                                          !  zfull(npoints,nlev) is bottom level of model
     
    
  Nlevels = gbx%Nlevels
  Npoints = gbx%Npoints
  ! Levels from TOA to surface
  zfull  = gbx%zlev(:,Nlevels:1:-1)
  at     = gbx%T(:,Nlevels:1:-1) 
  dtau_s = gbx%dtau_s(:,Nlevels:1:-1) 
  dtau_s_snow = gbx%dtau_s_snow(:,Nlevels:1:-1) 
  dtau_c = gbx%dtau_c(:,Nlevels:1:-1) 
  frac_out(1:Npoints,:,1:Nlevels) = sgx%frac_out(1:Npoints,:,Nlevels:1:-1)
  prec_frac(1:Npoints,:,1:Nlevels) = sgx%prec_frac(1:Npoints,:,Nlevels:1:-1)
  sunlit = int(gbx%sunlit)
 
  call MISR_simulator(gbx%npoints,gbx%nlevels,gbx%ncolumns,&
                     sunlit,zfull,at,dtau_s,dtau_c,dtau_s_snow,frac_out,prec_frac, R_UNDEF, &
                     y%fq_MISR,y%MISR_dist_model_layertops,y%MISR_meanztop,y%MISR_cldarea)
            
END SUBROUTINE COSP_MISR_SIMULATOR

END MODULE MOD_COSP_MISR_SIMULATOR
