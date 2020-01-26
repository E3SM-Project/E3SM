! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, Centre National de la Recherche Scientifique
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
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - optimization for vectorization
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_parasol
  USE COSP_KINDS,          ONLY: wp
  USE COSP_MATH_CONSTANTS, ONLY: pi
  use mod_cosp_config,     ONLY: R_UNDEF,PARASOL_NREFL,PARASOL_NTAU,PARASOL_TAU,PARASOL_SZA,rlumA,rlumB
  implicit none

contains
  SUBROUTINE parasol_subcolumn(npoints,nrefl,tautot_S_liq,tautot_S_ice,refl)
    ! ##########################################################################
    ! Purpose: To compute Parasol reflectance signal from model-simulated profiles 
    !          of cloud water and cloud fraction in each sub-column of each model 
    !          gridbox.
    !
    !
    ! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
    ! - optimization for vectorization
    !
    ! Version 2.0 (October 2008)
    ! Version 2.1 (December 2008)
    ! ##########################################################################
    
    ! INPUTS
    INTEGER,intent(in) :: &
         npoints,              & ! Number of horizontal gridpoints
         nrefl                   ! Number of angles for which the reflectance is computed
    REAL(WP),intent(inout),dimension(npoints) :: &
         tautot_S_liq,         & ! Liquid water optical thickness, from TOA to SFC
         tautot_S_ice            ! Ice water optical thickness, from TOA to SFC
    ! OUTPUTS
    REAL(WP),intent(inout),dimension(npoints,nrefl) :: &
         refl                    ! Parasol reflectances
    
    ! LOCAL VARIABLES
    REAL(WP),dimension(npoints) :: &
         tautot_S,             & ! Cloud optical thickness, from TOA to surface
         frac_taucol_liq,      & !
         frac_taucol_ice         !
    
    ! Look up table variables:
    INTEGER                            :: ny,it 
    REAL(WP),dimension(PARASOL_NREFL)         :: r_norm
    REAL(WP),dimension(PARASOL_NREFL,PARASOL_NTAU-1) :: aa,ab,ba,bb
    REAL(WP),dimension(npoints,5)      :: rlumA_mod,rlumB_mod
    
    !--------------------------------------------------------------------------------
    ! Lum_norm=f(PARASOL_SZA,tau_cloud) derived from adding-doubling calculations
    !        valid ONLY ABOVE OCEAN (albedo_sfce=5%)
    !        valid only in one viewing direction (theta_v=30�, phi_s-phi_v=320�)
    !        based on adding-doubling radiative transfer computation
    !        for PARASOL_TAU values (0 to 100) and for PARASOL_SZA values (0 to 80)
    !        for 2 scattering phase functions: liquid spherical, ice non spherical
    
    ! Initialize
    rlumA_mod(1:npoints,1:5) = 0._wp
    rlumB_mod(1:npoints,1:5) = 0._wp

    r_norm(1:PARASOL_NREFL)=1._wp/ cos(pi/180._wp*PARASOL_SZA(1:PARASOL_NREFL))
    
    tautot_S_liq(1:npoints) = max(tautot_S_liq(1:npoints),PARASOL_TAU(1))
    tautot_S_ice(1:npoints) = max(tautot_S_ice(1:npoints),PARASOL_TAU(1))
    tautot_S(1:npoints)     = tautot_S_ice(1:npoints) + tautot_S_liq(1:npoints)

    ! Relative fraction of the opt. thick due to liquid or ice clouds
    WHERE (tautot_S(1:npoints) .gt. 0.)
       frac_taucol_liq(1:npoints) = tautot_S_liq(1:npoints) / tautot_S(1:npoints)
       frac_taucol_ice(1:npoints) = tautot_S_ice(1:npoints) / tautot_S(1:npoints)
    ELSEWHERE
       frac_taucol_liq(1:npoints) = 1._wp
       frac_taucol_ice(1:npoints) = 0._wp
    END WHERE
    tautot_S(1:npoints)=MIN(tautot_S(1:npoints),PARASOL_TAU(PARASOL_NTAU))
    
    ! Linear interpolation    
    DO ny=1,PARASOL_NTAU-1
       ! Microphysics A (liquid clouds) 
       aA(1:PARASOL_NREFL,ny) = (rlumA(1:PARASOL_NREFL,ny+1)-rlumA(1:PARASOL_NREFL,ny))/(PARASOL_TAU(ny+1)-PARASOL_TAU(ny))
       bA(1:PARASOL_NREFL,ny) = rlumA(1:PARASOL_NREFL,ny) - aA(1:PARASOL_NREFL,ny)*PARASOL_TAU(ny)
       ! Microphysics B (ice clouds)
       aB(1:PARASOL_NREFL,ny) = (rlumB(1:PARASOL_NREFL,ny+1)-rlumB(1:PARASOL_NREFL,ny))/(PARASOL_TAU(ny+1)-PARASOL_TAU(ny))
       bB(1:PARASOL_NREFL,ny) = rlumB(1:PARASOL_NREFL,ny) - aB(1:PARASOL_NREFL,ny)*PARASOL_TAU(ny)
    ENDDO
    
    DO it=1,PARASOL_NREFL
       DO ny=1,PARASOL_NTAU-1
          WHERE (tautot_S(1:npoints) .ge. PARASOL_TAU(ny).and. &
                 tautot_S(1:npoints) .le. PARASOL_TAU(ny+1))
             rlumA_mod(1:npoints,it) = aA(it,ny)*tautot_S(1:npoints) + bA(it,ny)
             rlumB_mod(1:npoints,it) = aB(it,ny)*tautot_S(1:npoints) + bB(it,ny)
          END WHERE
       END DO
    END DO
    
    DO it=1,PARASOL_NREFL
       refl(1:npoints,it) = frac_taucol_liq(1:npoints) * rlumA_mod(1:npoints,it) &
            + frac_taucol_ice(1:npoints) * rlumB_mod(1:npoints,it)
       ! Normalized radiance -> reflectance: 
       refl(1:npoints,it) = refl(1:npoints,it) * r_norm(it)
    ENDDO
    
    RETURN
  END SUBROUTINE parasol_subcolumn
  ! ######################################################################################
  ! SUBROUTINE parasol_gridbox
  ! ######################################################################################
  subroutine parasol_column(npoints,nrefl,ncol,land,refl,parasolrefl)

    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal grid points
         ncol,    & ! Number of subcolumns
         nrefl      ! Number of solar zenith angles for parasol reflectances
    real(wp),intent(in),dimension(npoints) :: &
         land       ! Landmask [0 - Ocean, 1 - Land]
    real(wp),intent(in),dimension(npoints,ncol,nrefl) :: &
         refl       ! Subgrid parasol reflectance ! parasol

    ! Outputs
    real(wp),intent(out),dimension(npoints,nrefl) :: &
         parasolrefl   ! Grid-averaged parasol reflectance

    ! Local variables
    integer :: k,ic

    ! Compute grid-box averaged Parasol reflectances
    parasolrefl(:,:) = 0._wp
    do k = 1, nrefl
       do ic = 1, ncol
          parasolrefl(:,k) = parasolrefl(:,k) + refl(:,ic,k)
       enddo
    enddo
    
    do k = 1, nrefl
       parasolrefl(:,k) = parasolrefl(:,k) / float(ncol)
       ! if land=1 -> parasolrefl=R_UNDEF
       ! if land=0 -> parasolrefl=parasolrefl
       parasolrefl(:,k) = parasolrefl(:,k) * MAX(1._wp-land(:),0.0) &
            + (1._wp - MAX(1._wp-land(:),0.0))*R_UNDEF
    enddo
  end subroutine parasol_column

end module mod_parasol
