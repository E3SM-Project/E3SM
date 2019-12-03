! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, Roger Marchand, version 1.2
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
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_MISR_SIMULATOR
  use cosp_kinds,      only: wp
  use MOD_COSP_STATS,  ONLY: hist2D
  use mod_cosp_config, ONLY: R_UNDEF,numMISRHgtBins,numMISRTauBins,misr_histHgt,  &
                             misr_histTau
  implicit none 

  ! Parameters
  real(wp),parameter :: &
       misr_taumin = 0.3_wp,            & ! Minimum optical depth for joint-histogram
       tauchk      = -1.*log(0.9999999)   ! Lower limit on optical depth

contains

  ! ######################################################################################
  ! SUBROUTINE misr_subcolumn
  ! ######################################################################################
  SUBROUTINE MISR_SUBCOLUMN(npoints,ncol,nlev,dtau,zfull,at,sunlit,tauOUT,               &
                        dist_model_layertops,box_MISR_ztop)
    ! INPUTS
    INTEGER, intent(in) :: &
         npoints,    & ! Number of horizontal gridpoints
         ncol,       & ! Number of subcolumns
         nlev          ! Number of vertical layers
    INTEGER, intent(in),dimension(npoints) :: &
         sunlit        ! 1 for day points, 0 for night time
    REAL(WP),intent(in),dimension(npoints,ncol,nlev) :: &
         dtau          ! Optical thickness
    REAL(WP),intent(in),dimension(npoints,nlev) :: &
         zfull,      & ! Height of full model levels (i.e. midpoints), [nlev] is bottom
         at            ! Temperature (K)

    ! OUTPUTS
    REAL(WP),intent(out),dimension(npoints,ncol) :: &
         box_MISR_ztop,     & ! Cloud-top height in each column
         tauOUT               ! Optical depth in each column
    REAL(WP),intent(out),dimension(npoints,numMISRHgtBins) :: &
         dist_model_layertops ! 

    ! INTERNAL VARIABLES
    INTEGER :: ilev,j,loop,ibox,thres_crossed_MISR
    INTEGER :: iMISR_ztop
    REAL(WP) :: cloud_dtau,MISR_penetration_height,ztest

    ! ############################################################################
    ! Initialize
    box_MISR_ztop(1:npoints,1:ncol) = 0._wp  

    do j=1,npoints

       ! Estimate distribution of Model layer tops
       dist_model_layertops(j,:)=0
       do ilev=1,nlev
          ! Define location of "layer top"
          if(ilev.eq.1 .or. ilev.eq.nlev) then
             ztest=zfull(j,ilev)
          else
             ztest=0.5_wp*(zfull(j,ilev)+zfull(j,ilev-1))
          endif

          ! Find MISR layer that contains this level
          ! *NOTE* the first MISR level is "no height" level
          iMISR_ztop=2
          do loop=2,numMISRHgtBins
             if ( ztest .gt. 1000*misr_histHgt(loop+1) ) then
                iMISR_ztop=loop+1
             endif
          enddo
          
          dist_model_layertops(j,iMISR_ztop) = dist_model_layertops(j,iMISR_ztop)+1
       enddo

       ! For each GCM cell or horizontal model grid point   
       do ibox=1,ncol
          ! Compute optical depth as a cummulative distribution in the vertical (nlev).
          tauOUT(j,ibox)=sum(dtau(j,ibox,1:nlev))

          thres_crossed_MISR=0
          do ilev=1,nlev
             ! If there a cloud, start the counter and store this height
             if(thres_crossed_MISR .eq. 0 .and. dtau(j,ibox,ilev) .gt. 0.) then
                ! First encountered a "cloud"
                thres_crossed_MISR = 1  
                cloud_dtau         = 0            
             endif

             if( thres_crossed_MISR .lt. 99 .and. thres_crossed_MISR .gt. 0 ) then
                if( dtau(j,ibox,ilev) .eq. 0.) then
                   ! We have come to the end of the current cloud layer without yet 
                   ! selecting a CTH boundary. Restart cloud tau counter 
                   cloud_dtau=0
                else
                   ! Add current optical depth to count for the current cloud layer
                   cloud_dtau=cloud_dtau+dtau(j,ibox,ilev)
                endif
                
                ! If the cloud is continuous but optically thin (< 1) from above the 
                ! current layer cloud top to the current level then MISR will like 
                ! see a top below the top of the current layer.
                if( dtau(j,ibox,ilev).gt.0 .and. (cloud_dtau-dtau(j,ibox,ilev)) .lt. 1) then
                   if(dtau(j,ibox,ilev) .lt. 1 .or. ilev.eq.1 .or. ilev.eq.nlev) then
                      ! MISR will likely penetrate to some point within this layer ... the middle
                      MISR_penetration_height=zfull(j,ilev)
                   else
                      ! Take the OD = 1.0 level into this layer
                      MISR_penetration_height=0.5_wp*(zfull(j,ilev)+zfull(j,ilev-1)) - &
                           0.5_wp*(zfull(j,ilev-1)-zfull(j,ilev+1))/dtau(j,ibox,ilev) 
                   endif
                   box_MISR_ztop(j,ibox)=MISR_penetration_height
                endif
                
                ! Check for a distinctive water layer
                if(dtau(j,ibox,ilev) .gt. 1 .and. at(j,ilev) .gt. 273 ) then
                   ! Must be a water cloud, take this as CTH level
                   thres_crossed_MISR=99
                endif
                
                ! If the total column optical depth is "large" than MISR can't see
                ! anything else. Set current point as CTH level
                if(sum(dtau(j,ibox,1:ilev)) .gt. 5) then
                   thres_crossed_MISR=99           
                endif
             endif
          enddo  
          
          ! Check to see if there was a cloud for which we didn't 
          ! set a MISR cloud top boundary
          if( thres_crossed_MISR .eq. 1) then
             ! If the cloud has a total optical depth of greater
             ! than ~ 0.5 MISR will still likely pick up this cloud
             ! with a height near the true cloud top
             ! otherwise there should be no CTH
             if(sum(dtau(j,ibox,1:nlev)) .gt. 0.5) then
                ! keep MISR detected CTH
             elseif(sum(dtau(j,ibox,1:nlev)) .gt. 0.2) then
                ! MISR may detect but wont likley have a good height
                box_MISR_ztop(j,ibox)=-1
             else
                ! MISR not likely to even detect.
                ! so set as not cloudy
                box_MISR_ztop(j,ibox)=0
             endif
          endif
       enddo  ! loop of subcolumns
       
    enddo    ! loop of gridpoints
    
    ! Modify MISR CTH for satellite spatial / pattern matcher effects
    ! Code in this region added by roj 5/2006 to account
    ! for spatial effect of the MISR pattern matcher.
    ! Basically, if a column is found between two neighbors
    ! at the same CTH, and that column has no hieght or
    ! a lower CTH, THEN misr will tend to but place the
    ! odd column at the same height as it neighbors.
    
    ! This setup assumes the columns represent a about a 1 to 4 km scale
    ! it will need to be modified significantly, otherwise
!	! DS2015: Add loop over gridpoints and index accordingly.
!    if(ncol.eq.1) then
!       ! Adjust based on neightboring points.
!       do j=2,npoints-1   
!          if(box_MISR_ztop(j-1,1) .gt. 0                             .and. &
!             box_MISR_ztop(j+1,1) .gt. 0                             .and. &
!             abs(box_MISR_ztop(j-1,1)-box_MISR_ztop(j+1,1)) .lt. 500 .and. &
!             box_MISR_ztop(j,1) .lt. box_MISR_ztop(j+1,1)) then
!             box_MISR_ztop(j,1) = box_MISR_ztop(j+1,1)    
!          endif
!       enddo
!    else
!       ! Adjust based on neighboring subcolumns.
!       do j=1,npoints
!          do ibox=2,ncol-1  
!                 if(box_MISR_ztop(j,ibox-1) .gt. 0                                .and. &
!                 box_MISR_ztop(j,ibox+1) .gt. 0                                .and. &
!                 abs(box_MISR_ztop(j,ibox-1)-box_MISR_ztop(j,ibox+1)) .lt. 500 .and. &
!                 box_MISR_ztop(j,ibox) .lt. box_MISR_ztop(j,ibox+1)) then
!                 box_MISR_ztop(j,ibox) = box_MISR_ztop(j,ibox+1)    
!               endif
!          enddo
!       enddo
!    endif
!    ! DS2015 END
     
    ! Fill dark scenes 
    do j=1,numMISRHgtBins
       where(sunlit .ne. 1) dist_model_layertops(1:npoints,j) = R_UNDEF
    enddo

  end SUBROUTINE MISR_SUBCOLUMN

  ! ######################################################################################
  ! SUBROUTINE misr_column
  ! ######################################################################################
  SUBROUTINE MISR_COLUMN(npoints,ncol,box_MISR_ztop,sunlit,tau,MISR_cldarea,MISR_mean_ztop,fq_MISR_TAU_v_CTH)

    ! INPUTS
    INTEGER, intent(in) :: &
         npoints,        & ! Number of horizontal gridpoints
         ncol              ! Number of subcolumns
    INTEGER, intent(in),dimension(npoints) :: &
         sunlit            ! 1 for day points, 0 for night time
    REAL(WP),intent(in),dimension(npoints,ncol) :: &
         box_MISR_ztop,  & ! Cloud-top height in each column
         tau               ! Column optical thickness

    ! OUTPUTS
    REAL(WP),intent(inout),dimension(npoints) :: &
         MISR_cldarea,   & ! Fraction area covered by clouds
         MISR_mean_ztop    ! Mean cloud top height MISR would observe
    REAL(WP),intent(inout),dimension(npoints,7,numMISRHgtBins) :: &
         fq_MISR_TAU_v_CTH ! Joint histogram of cloud-cover and tau

    ! INTERNAL VARIABLES
    INTEGER :: j
    LOGICAL,dimension(ncol) :: box_cloudy 
    real(wp),dimension(npoints,ncol) :: tauWRK,box_MISR_ztopWRK
    ! ############################################################################

    ! Compute column quantities and joint-histogram
    MISR_cldarea(1:npoints)                       = 0._wp
    MISR_mean_ztop(1:npoints)                     = 0._wp
    fq_MISR_TAU_v_CTH(1:npoints,1:7,1:numMISRHgtBins) = 0._wp
    tauWRK(1:npoints,1:ncol)                      = tau(1:npoints,1:ncol)
    box_MISR_ztopWRK(1:npoints,1:ncol)            = box_MISR_ztop(1:npoints,1:ncol)
    do j=1,npoints

       ! Subcolumns that are cloudy(true) and not(false)
       box_cloudy(1:ncol) = merge(.true.,.false.,tau(j,1:ncol) .gt. tauchk)

       ! Fill optically thin clouds with fill value
       where(.not. box_cloudy(1:ncol)) tauWRK(j,1:ncol)  = -999._wp
       where(box_MISR_ztopWRK(j,1:ncol) .eq. 0) box_MISR_ztopWRK(j,1:ncol)=-999._wp

       ! Compute joint histogram and column quantities for points that are sunlit and cloudy
       if (sunlit(j) .eq. 1) then 
          ! Joint histogram
          call hist2D(tauWRK(j,1:ncol),box_MISR_ztopWRK(j,1:ncol),ncol,misr_histTau,numMISRTauBins,&
               1000*misr_histHgt,numMISRHgtBins,fq_MISR_TAU_v_CTH(j,1:numMISRTauBins,1:numMISRHgtBins))
          fq_MISR_TAU_v_CTH(j,1:numMISRTauBins,1:numMISRHgtBins) =                       &
             100._wp*fq_MISR_TAU_v_CTH(j,1:numMISRTauBins,1:numMISRHgtBins)/ncol

          ! Column cloud area
          MISR_cldarea(j)=real(count(box_MISR_ztopWRK(j,1:ncol) .ne. -999.))/ncol

          ! Column cloud-top height
          if ( count(box_MISR_ztopWRK(j,1:ncol) .ne. -999.) .ne. 0 ) then
             MISR_mean_ztop(j) = sum(box_MISR_ztopWRK(j,1:ncol),box_MISR_ztopWRK(j,1:ncol) .ne. -999.)/ &
                  count(box_MISR_ztopWRK(j,1:ncol) .ne. -999.)
          else
             MISR_mean_ztop(j) = R_UNDEF
          endif

       else
          MISR_cldarea(j)         = R_UNDEF
          MISR_mean_ztop(npoints) = R_UNDEF
       endif
    enddo

  end SUBROUTINE MISR_COLUMN

end MODULE MOD_MISR_SIMULATOR
