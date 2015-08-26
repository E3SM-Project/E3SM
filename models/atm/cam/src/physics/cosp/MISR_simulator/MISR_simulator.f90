! 
! Copyright (c) 2009,  Roger Marchand, version 1.2
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list of 
!       conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list 
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the University of Washington nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
! BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
! SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

      SUBROUTINE MISR_simulator(
     &     npoints,
     &     nlev,
     &     ncol,
     &     sunlit,
     &     zfull,
     &     at,
     &     dtau_s,
     &     dtau_c,
     &     dtau_s_snow,
     &     frac_out,
     &     prec_frac,
     &     missing_value,
     &     fq_MISR_TAU_v_CTH,
     &     dist_model_layertops,
     &     MISR_mean_ztop,
     &     MISR_cldarea
     & )
    

      implicit none
      integer n_MISR_CTH
      parameter(n_MISR_CTH=16)
         
!     -----
!     Input 
!     -----

      INTEGER npoints                   !  if ncol ==1, the number of model points in the horizontal grid  
                            !   else    the number of GCM grid points
                            
      INTEGER nlev                      !  number of model vertical levels
      
      INTEGER ncol                      !  number of model sub columns 
                        !  (must already be generated in via scops and passed to this
                        !   routine via the variable frac_out )
  
      INTEGER sunlit(npoints)           !  1 for day points, 0 for night time

      REAL zfull(npoints,nlev)          !  height (in meters) of full model levels (i.e. midpoints)
                                        !  zfull(npoints,1)    is    top level of model
                                        !  zfull(npoints,nlev) is bottom level of model (closest point to surface)  

      REAL at(npoints,nlev)             !  temperature in each model level (K)
 
      REAL dtau_s(npoints,nlev)         !  visible wavelength cloud optical depth ... for "stratiform" condensate
                                        !  NOTE:  this the cloud optical depth of only the
                    			!     the model cell (i,j)

                    
      REAL dtau_c(npoints,nlev)         !  visible wavelength cloud optical depth ... for "convective" condensate
                                        !  NOTE:  this the cloud optical depth of only the
                    			!     the model cell (i,j)
      
      REAL dtau_s_snow(npoints,nlev)      !  visible wavelength SNOW optical depth ... 
                                     
      REAL frac_out(npoints,ncol,nlev)  !  NOTE: only need if columns>1 ... subgrid scheme in use.
      REAL prec_frac(npoints,ncol,nlev)  !  same as frac_out but for precipitation -- for Steve's snow scheme 
      
      REAL missing_value
                                 
!     ------
!     Outputs
!     ------
            
      REAL fq_MISR_TAU_v_CTH(npoints,7,n_MISR_CTH)      
      REAL dist_model_layertops(npoints,n_MISR_CTH)
      REAL MISR_cldarea(npoints)               ! fractional area coverged by clouds 
      REAL MISR_mean_ztop(npoints)             ! mean cloud top hieght(m) MISR would observe
                                   ! NOTE: == 0 if area ==0
                            

!     ------
!     Working variables 
!     ------

      REAL tau(npoints,ncol)        ! total column optical depth ... 

      INTEGER j,ilev,ilev2,ibox,k
      INTEGER itau
         
      LOGICAL box_cloudy(npoints,ncol)
      
      real isccp_taumin
      real boxarea
      real tauchk
      REAL box_MISR_ztop(npoints,ncol)  ! cloud top hieght(m) MISR would observe
      
      integer thres_crossed_MISR 
      integer loop,iMISR_ztop
      
      real dtau, cloud_dtau, MISR_penetration_height,ztest     
      
      real MISR_CTH_boundaries(n_MISR_CTH+1)
      
      DATA MISR_CTH_boundaries / -99, 0, 0.5, 1, 1.5, 2, 2.5, 3,
     c                    4, 5, 7, 9, 11, 13, 15, 17, 99 /
      
      DATA isccp_taumin / 0.3 /
    
      tauchk = -1.*log(0.9999999)
        
      !
      ! For each GCM cell or horizontal model grid point ...
      ! 
      do j=1,npoints    

         !
         !  estimate distribution of Model layer tops
         !  
         dist_model_layertops(j,:)=0

       do ilev=1,nlev 
            
        ! define location of "layer top"
        if(ilev.eq.1 .or. ilev.eq.nlev) then
            ztest=zfull(j,ilev)
        else
            ztest=0.5*(zfull(j,ilev)+zfull(j,ilev-1)) 
        endif   

        ! find MISR layer that contains this level
        ! note, the first MISR level is "no height" level
        iMISR_ztop=2
        do loop=2,n_MISR_CTH
        
            if ( ztest .gt.
     &                1000*MISR_CTH_boundaries(loop+1) ) then
        
                iMISR_ztop=loop+1
            endif
        enddo

        dist_model_layertops(j,iMISR_ztop)=
     &          dist_model_layertops(j,iMISR_ztop)+1
       enddo
    
    
         !
         ! compute total cloud optical depth for each column
         !       
       do ibox=1,ncol     
       
        ! Initialize tau to zero in each subcolum
            tau(j,ibox)=0. 
        box_cloudy(j,ibox)=.false.
        box_MISR_ztop(j,ibox)=0  
        
        ! initialize threshold detection for each sub column 
        thres_crossed_MISR=0;
       
        do ilev=1,nlev
     
             dtau=0
             
                 if (frac_out(j,ibox,ilev).eq.1) then
                        dtau = dtau_s(j,ilev)
                 endif
                 
                 if (frac_out(j,ibox,ilev).eq.2) then
                        dtau = dtau_c(j,ilev)
                 end if 
                 
		 if ((prec_frac(j,ibox,ilev).eq.1) .or.
     &               (prec_frac(j,ibox,ilev).eq.3)) then
                         dtau = dtau + dtau_s_snow(j,ilev)
                 end if    
 
             tau(j,ibox)=tau(j,ibox)+ dtau
              
                     
        ! NOW for MISR ..
        ! if there a cloud ... start the counter ... store this height
        if(thres_crossed_MISR .eq. 0 .and. dtau .gt. 0.) then
        
            ! first encountered a "cloud"
            thres_crossed_MISR=1  
            cloud_dtau=0            
        endif   
                
        if( thres_crossed_MISR .lt. 99 .and.
     &              thres_crossed_MISR .gt. 0 ) then
     
                if( dtau .eq. 0.) then
        
                    ! we have come to the end of the current cloud
                ! layer without yet selecting a CTH boundary.
                ! ... restart cloud tau counter 
                cloud_dtau=0
            else
                ! add current optical depth to count for 
                ! the current cloud layer
                cloud_dtau=cloud_dtau+dtau
            endif
                
            ! if the cloud is continuous but optically thin (< 1)
            ! from above the current layer cloud top to the current level
            ! then MISR will like see a top below the top of the current 
            ! layer
            if( dtau.gt.0 .and. (cloud_dtau-dtau) .lt. 1) then
            
                if(dtau .lt. 1 .or. ilev.eq.1 .or. ilev.eq.nlev) then

                    ! MISR will likely penetrate to some point
                    ! within this layer ... the middle
                    MISR_penetration_height=zfull(j,ilev)

                else
                    ! take the OD = 1.0 level into this layer
                    MISR_penetration_height=
     &                     0.5*(zfull(j,ilev)+zfull(j,ilev-1)) - 
     &                     0.5*(zfull(j,ilev-1)-zfull(j,ilev+1))
     &                  /dtau 
                endif   

                box_MISR_ztop(j,ibox)=MISR_penetration_height
                
            endif
        
            ! check for a distinctive water layer
            if(dtau .gt. 1 .and. at(j,ilev).gt.273 ) then
     
                    ! must be a water cloud ... 
                ! take this as CTH level
                thres_crossed_MISR=99
            endif
        
            ! if the total column optical depth is "large" than
            ! MISR can't seen anything else ... set current point as CTH level
            if(tau(j,ibox) .gt. 5) then 

                thres_crossed_MISR=99           
            endif

        endif ! MISR CTH booundary not set
        
        enddo  !ilev - loop over vertical levesl
    
        ! written by roj 5/2006
        ! check to see if there was a cloud for which we didn't 
        ! set a MISR cloud top boundary
        if( thres_crossed_MISR .eq. 1) then
    
        ! if the cloud has a total optical depth of greater
        ! than ~ 0.5 MISR will still likely pick up this cloud
        ! with a height near the true cloud top
        ! otherwise there should be no CTH
        if( tau(j,ibox) .gt. 0.5) then

            ! keep MISR detected CTH
            
        elseif(tau(j,ibox) .gt. 0.2) then

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
       

        !     
        !   Modify MISR CTH for satellite spatial / pattern matcher effects
    !
    !   Code in this region added by roj 5/2006 to account
    !   for spatial effect of the MISR pattern matcher.
    !   Basically, if a column is found between two neighbors
    !   at the same CTH, and that column has no hieght or
    !   a lower CTH, THEN misr will tend to but place the
    !   odd column at the same height as it neighbors.
    !
    !   This setup assumes the columns represent a about a 1 to 4 km scale
    !   it will need to be modified significantly, otherwise
        if(ncol.eq.1) then
    
       ! adjust based on neightboring points ... i.e. only 2D grid was input
           do j=2,npoints-1
            
            if(box_MISR_ztop(j-1,1).gt.0 .and. 
     &             box_MISR_ztop(j+1,1).gt.0       ) then

                if( abs( box_MISR_ztop(j-1,1) -  
     &                   box_MISR_ztop(j+1,1) ) .lt. 500 
     &              .and.
     &                   box_MISR_ztop(j,1) .lt. 
     &                   box_MISR_ztop(j+1,1)     ) then
            
                    box_MISR_ztop(j,1) =
     &                      box_MISR_ztop(j+1,1)    
                endif

            endif
         enddo
        else
         
         ! adjust based on neighboring subcolumns ....
         do ibox=2,ncol-1
            
            if(box_MISR_ztop(1,ibox-1).gt.0 .and. 
     &             box_MISR_ztop(1,ibox+1).gt.0        ) then

                if( abs( box_MISR_ztop(1,ibox-1) -  
     &                   box_MISR_ztop(1,ibox+1) ) .lt. 500 
     &              .and.
     &                   box_MISR_ztop(1,ibox) .lt. 
     &                   box_MISR_ztop(1,ibox+1)     ) then
            
                    box_MISR_ztop(1,ibox) =
     &                      box_MISR_ztop(1,ibox+1)    
                endif

            endif
         enddo
      
        endif

        !     
    !     DETERMINE CLOUD TYPE FREQUENCIES
    !
    !     Now that ztop and tau have been determined, 
    !     determine amount of each cloud type
        boxarea=1./real(ncol)  
        do j=1,npoints 

         ! reset frequencies -- modified loop structure, roj 5/2006 
         do ilev=1,7  ! "tau loop"  
            do  ilev2=1,n_MISR_CTH                      
            fq_MISR_TAU_v_CTH(j,ilev,ilev2)=0.     
            enddo
         enddo
           
         MISR_cldarea(j)=0.
         MISR_mean_ztop(j)=0.

         do ibox=1,ncol

            if (tau(j,ibox) .gt. (tauchk)) then
               box_cloudy(j,ibox)=.true.
            endif
  
            itau = 0
        
            if (box_cloudy(j,ibox)) then
    
          !determine optical depth category
              if (tau(j,ibox) .lt. isccp_taumin) then
                  itau=1
              else if (tau(j,ibox) .ge. isccp_taumin                                    
     &          .and. tau(j,ibox) .lt. 1.3) then
                  itau=2
              else if (tau(j,ibox) .ge. 1.3 
     &          .and. tau(j,ibox) .lt. 3.6) then
                  itau=3
              else if (tau(j,ibox) .ge. 3.6 
     &          .and. tau(j,ibox) .lt. 9.4) then
                  itau=4
              else if (tau(j,ibox) .ge. 9.4 
     &          .and. tau(j,ibox) .lt. 23.) then
                  itau=5
              else if (tau(j,ibox) .ge. 23. 
     &          .and. tau(j,ibox) .lt. 60.) then
                  itau=6
              else if (tau(j,ibox) .ge. 60.) then
                  itau=7
              endif
              
             endif  

       ! update MISR histograms and summary metrics - roj 5/2005
       if (sunlit(j).eq.1) then 
                     
              !if cloudy added by roj 5/2005
          if( box_MISR_ztop(j,ibox).eq.0) then
          
            ! no cloud detected
            iMISR_ztop=0

          elseif( box_MISR_ztop(j,ibox).eq.-1) then

            ! cloud can be detected but too thin to get CTH
            iMISR_ztop=1    

            fq_MISR_TAU_v_CTH(j,itau,iMISR_ztop)=
     &            fq_MISR_TAU_v_CTH(j,itau,iMISR_ztop) + boxarea

          else
            
            !
            ! determine index for MISR bin set
            !

            iMISR_ztop=2
            
            do loop=2,n_MISR_CTH
        
                if ( box_MISR_ztop(j,ibox) .gt.
     &                1000*MISR_CTH_boundaries(loop+1) ) then
        
                  iMISR_ztop=loop+1

                endif
            enddo
          
            if(box_cloudy(j,ibox)) then
            
               ! there is an isccp clouds so itau(j) is defined
               fq_MISR_TAU_v_CTH(j,itau,iMISR_ztop)=
     &            fq_MISR_TAU_v_CTH(j,itau,iMISR_ztop) + boxarea
     
            else
                ! MISR CTH resolution is trying to fill in a
                ! broken cloud scene where there is no condensate.
                ! The MISR CTH-1D-OD product will only put in a cloud
                ! if the MISR cloud mask indicates cloud.
                ! therefore we will not include this column in the histogram
                ! in reality aerosoal and 3D effects or bright surfaces
                ! could fool the MISR cloud mask

                ! the alternative is to count as very thin cloud ??
!               fq_MISR_TAU_v_CTH(1,iMISR_ztop)=
!     &                     fq_MISR_TAU_v_CTH(1,iMISR_ztop) + boxarea
            endif


            MISR_mean_ztop(j)=MISR_mean_ztop(j)+
     &                       box_MISR_ztop(j,ibox)*boxarea          

            MISR_cldarea(j)=MISR_cldarea(j) + boxarea 
 
          endif
       else
          ! Set to issing data. A. Bodas - 14/05/2010
          do loop=1,n_MISR_CTH
             do k=1,7
                fq_MISR_TAU_v_CTH(j,k,loop) = missing_value
             enddo
             dist_model_layertops(j,loop) = missing_value
          enddo
          MISR_cldarea(j) = missing_value
          MISR_mean_ztop(npoints) = missing_value

       endif ! is sunlight ?
       
       enddo ! ibox - loop over subcolumns          
      
       if( MISR_cldarea(j) .gt. 0.) then
        MISR_mean_ztop(j)= MISR_mean_ztop(j) / MISR_cldarea(j)   ! roj 5/2006
       endif

       enddo  ! loop over grid points

      return
      end 
