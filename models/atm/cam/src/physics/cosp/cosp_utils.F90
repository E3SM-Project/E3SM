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
! Jul 2007 - A. Bodas-Salcedo - Initial version
!

MODULE MOD_COSP_UTILS
  USE MOD_COSP_CONSTANTS
  IMPLICIT NONE

  INTERFACE Z_TO_DBZ
    MODULE PROCEDURE Z_TO_DBZ_2D,Z_TO_DBZ_3D,Z_TO_DBZ_4D
  END INTERFACE

  INTERFACE COSP_CHECK_INPUT
    MODULE PROCEDURE COSP_CHECK_INPUT_1D,COSP_CHECK_INPUT_2D,COSP_CHECK_INPUT_3D
  END INTERFACE
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_PRECIP_MXRATIO --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_PRECIP_MXRATIO(Npoints,Nlevels,Ncolumns,p,T,prec_frac,prec_type, &
                          n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2, &
                          flux,mxratio)

    ! Input arguments, (IN)
    integer,intent(in) :: Npoints,Nlevels,Ncolumns
    real,intent(in),dimension(Npoints,Nlevels) :: p,T,flux
    real,intent(in),dimension(Npoints,Ncolumns,Nlevels) :: prec_frac
    real,intent(in) :: n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,prec_type
    ! Input arguments, (OUT)
    real,intent(out),dimension(Npoints,Ncolumns,Nlevels) :: mxratio
    ! Local variables
    integer :: i,j,k
    real :: sigma,one_over_xip1,xi,rho0,rho
    
    mxratio = 0.0

    if (n_ax >= 0.0) then ! N_ax is used to control which hydrometeors need to be computed
        !gamma1  = gamma(alpha_x + b_x + d_x + 1.0)
        !gamma2  = gamma(alpha_x + b_x + 1.0)
        xi      = d_x/(alpha_x + b_x - n_bx + 1.0)
        rho0    = 1.29
        sigma   = (gamma2/(gamma1*c_x))*(n_ax*a_x*gamma2)**xi
        one_over_xip1 = 1.0/(xi + 1.0)
        
        do k=1,Nlevels
            do j=1,Ncolumns
                do i=1,Npoints
                    if ((prec_frac(i,j,k)==prec_type).or.(prec_frac(i,j,k)==3.)) then
                        rho = p(i,k)/(287.05*T(i,k))
                        mxratio(i,j,k)=(flux(i,k)*((rho/rho0)**g_x)*sigma)**one_over_xip1
                        mxratio(i,j,k)=mxratio(i,j,k)/rho
                    endif
                enddo
            enddo
        enddo
    endif
END SUBROUTINE COSP_PRECIP_MXRATIO


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE ZERO_INT -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELEMENTAL SUBROUTINE ZERO_INT(x,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                 y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                 y21,y22,y23,y24,y25,y26,y27,y28,y29,y30)

  integer,intent(inout) :: x
  integer,intent(inout),optional :: y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                    y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                    y21,y22,y23,y24,y25,y26,y27,y28,y29,y30
  x = 0
  if (present(y01)) y01 = 0
  if (present(y02)) y02 = 0
  if (present(y03)) y03 = 0
  if (present(y04)) y04 = 0
  if (present(y05)) y05 = 0
  if (present(y06)) y06 = 0
  if (present(y07)) y07 = 0
  if (present(y08)) y08 = 0
  if (present(y09)) y09 = 0
  if (present(y10)) y10 = 0
  if (present(y11)) y11 = 0
  if (present(y12)) y12 = 0
  if (present(y13)) y13 = 0
  if (present(y14)) y14 = 0
  if (present(y15)) y15 = 0
  if (present(y16)) y16 = 0
  if (present(y17)) y17 = 0
  if (present(y18)) y18 = 0
  if (present(y19)) y19 = 0
  if (present(y20)) y20 = 0
  if (present(y21)) y21 = 0
  if (present(y22)) y22 = 0
  if (present(y23)) y23 = 0
  if (present(y24)) y24 = 0
  if (present(y25)) y25 = 0
  if (present(y26)) y26 = 0
  if (present(y27)) y27 = 0
  if (present(y28)) y28 = 0
  if (present(y29)) y29 = 0
  if (present(y30)) y30 = 0
END SUBROUTINE  ZERO_INT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE ZERO_REAL ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELEMENTAL SUBROUTINE ZERO_REAL(x,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                 y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                 y21,y22,y23,y24,y25,y26,y27,y28,y29,y30)

  real,intent(inout) :: x
  real,intent(inout),optional :: y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                 y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                 y21,y22,y23,y24,y25,y26,y27,y28,y29,y30
  x = 0.0
  if (present(y01)) y01 = 0.0
  if (present(y02)) y02 = 0.0
  if (present(y03)) y03 = 0.0
  if (present(y04)) y04 = 0.0
  if (present(y05)) y05 = 0.0
  if (present(y06)) y06 = 0.0
  if (present(y07)) y07 = 0.0
  if (present(y08)) y08 = 0.0
  if (present(y09)) y09 = 0.0
  if (present(y10)) y10 = 0.0
  if (present(y11)) y11 = 0.0
  if (present(y12)) y12 = 0.0
  if (present(y13)) y13 = 0.0
  if (present(y14)) y14 = 0.0
  if (present(y15)) y15 = 0.0
  if (present(y16)) y16 = 0.0
  if (present(y17)) y17 = 0.0
  if (present(y18)) y18 = 0.0
  if (present(y19)) y19 = 0.0
  if (present(y20)) y20 = 0.0
  if (present(y21)) y21 = 0.0
  if (present(y22)) y22 = 0.0
  if (present(y23)) y23 = 0.0
  if (present(y24)) y24 = 0.0
  if (present(y25)) y25 = 0.0
  if (present(y26)) y26 = 0.0
  if (present(y27)) y27 = 0.0
  if (present(y28)) y28 = 0.0
  if (present(y29)) y29 = 0.0
  if (present(y30)) y30 = 0.0
END SUBROUTINE  ZERO_REAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_2D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_2D(mdi,z)
    real,intent(in) :: mdi
    real,dimension(:,:),intent(inout) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]
    
    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18*z
    where (z > 1.0e-6) ! Limit to -60 dBZ
      z = 10.0*log10(z)
    elsewhere
      z = mdi
    end where  
  END SUBROUTINE Z_TO_DBZ_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_3D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_3D(mdi,z)
    real,intent(in) :: mdi
    real,dimension(:,:,:),intent(inout) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]
    
    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18*z
    where (z > 1.0e-6) ! Limit to -60 dBZ
      z = 10.0*log10(z)
    elsewhere
      z = mdi
    end where  
  END SUBROUTINE Z_TO_DBZ_3D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_4D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_4D(mdi,z)
    real,intent(in) :: mdi
    real,dimension(:,:,:,:),intent(inout) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]
    
    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18*z
    where (z > 1.0e-6) ! Limit to -60 dBZ
      z = 10.0*log10(z)
    elsewhere
      z = mdi
    end where  
  END SUBROUTINE Z_TO_DBZ_4D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_1D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_1D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real,intent(inout) :: x(:)
    real,intent(in),optional :: min_val,max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_1D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_1D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_2D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_2D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real,intent(inout) :: x(:,:)
    real,intent(in),optional :: min_val,max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_2D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_3D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_3D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real,intent(inout) :: x(:,:,:)
    real,intent(in),optional :: min_val,max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_3D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_3D


END MODULE MOD_COSP_UTILS
