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
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! May 2015 - Dustin Swales    - Modified for COSPv2.0
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_UTILS
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG
  IMPLICIT NONE

CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_PRECIP_MXRATIO --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_PRECIP_MXRATIO(Npoints,Nlevels,Ncolumns,p,T,prec_frac,prec_type, &
                          n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,gamma3,gamma4, &
                          flux,mxratio,reff)

    ! Input arguments, (IN)
    integer,intent(in) :: Npoints,Nlevels,Ncolumns
    real(wp),intent(in),dimension(Npoints,Nlevels) :: p,T,flux
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: prec_frac
    real(wp),intent(in) :: n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,gamma3,gamma4,prec_type
    ! Input arguments, (OUT)
    real(wp),intent(out),dimension(Npoints,Ncolumns,Nlevels) :: mxratio
    real(wp),intent(inout),dimension(Npoints,Ncolumns,Nlevels) :: reff
    ! Local variables
    integer :: i,j,k
    real(wp) :: sigma,one_over_xip1,xi,rho0,rho,lambda_x,gamma_4_3_2,delta
    
    mxratio = 0.0

    if (n_ax >= 0.0) then ! N_ax is used to control which hydrometeors need to be computed
        xi      = d_x/(alpha_x + b_x - n_bx + 1._wp)
        rho0    = 1.29_wp
        sigma   = (gamma2/(gamma1*c_x))*(n_ax*a_x*gamma2)**xi
        one_over_xip1 = 1._wp/(xi + 1._wp)
        gamma_4_3_2 = 0.5_wp*gamma4/gamma3
        delta = (alpha_x + b_x + d_x - n_bx + 1._wp)
        
        do k=1,Nlevels
            do j=1,Ncolumns
                do i=1,Npoints
                    if ((prec_frac(i,j,k)==prec_type).or.(prec_frac(i,j,k)==3.)) then
                        rho = p(i,k)/(287.05_wp*T(i,k))
                        mxratio(i,j,k)=(flux(i,k)*((rho/rho0)**g_x)*sigma)**one_over_xip1
                        mxratio(i,j,k)=mxratio(i,j,k)/rho
                        ! Compute effective radius
                        if ((reff(i,j,k) <= 0._wp).and.(flux(i,k) /= 0._wp)) then
                           lambda_x = (a_x*c_x*((rho0/rho)**g_x)*n_ax*gamma1/flux(i,k))**(1._wp/delta)
                           reff(i,j,k) = gamma_4_3_2/lambda_x
                        endif
                    endif
                enddo
            enddo
        enddo
    endif
END SUBROUTINE COSP_PRECIP_MXRATIO


END MODULE MOD_COSP_UTILS
