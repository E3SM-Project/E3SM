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
MODULE MOD_COSP_ISCCP_INTERFACE
  USE COSP_KINDS,      ONLY: wp
  USE mod_icarus,      ONLY: isccp_top_height,isccp_top_height_direction
  IMPLICIT NONE
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                  TYPE isccp_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Derived input type for ISCCP simulator
  type isccp_IN
     integer,pointer  ::       &
          Npoints,             & ! Number of gridpoints.
          Ncolumns,            & ! Number of columns.
          Nlevels,             & ! Number of levels.
          top_height,          & !
          top_height_direction   !
     integer,pointer ::        &
          sunlit(:)              ! Sunlit points (npoints)
     real(wp),pointer ::       &
          emsfc_lw
     real(wp),pointer ::       &
          skt(:)                 ! Surface temperature (npoints)
     real(wp),pointer ::       &
          at(:,:),             & ! Temperature (npoint,nlev)
          pfull(:,:),          & ! Pressure (npoints,nlev)
          qv(:,:)                ! Specific humidity (npoints,nlev)
     real(wp),pointer  ::       &          
          phalf(:,:)             ! Pressure at half levels (npoints,nlev+1)
     real(wp),pointer ::       &
          frac_out(:,:,:),     & ! Cloud fraction (npoints,ncolumns,nlevels)
          dtau(:,:,:),         & ! Optical depth (npoints,ncolumns,nlevels)
          dem(:,:,:)             ! Emissivity (npoints,ncolumns,nlevels)
  end type isccp_IN

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE cosp_isccp_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ISCCP_INIT(top_height,top_height_direction)
     integer,intent(in) :: &
         top_height, &
         top_height_direction

    ! Cloud-top height determination
    isccp_top_height           = top_height
    isccp_top_height_direction = top_height_direction

  END SUBROUTINE COSP_ISCCP_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_ISCCP_INTERFACE
