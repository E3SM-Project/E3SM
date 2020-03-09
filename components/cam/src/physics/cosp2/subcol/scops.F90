! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, British Crown Copyright, the Met Office
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
module mod_scops
  USE COSP_KINDS,     ONLY: wp
  USE MOD_RNG!,        ONLY: rng_state,get_rng
  use mod_cosp_error, ONLY: errorMessage

  implicit none

  integer,parameter :: default_overlap = 3 ! Used when invalid overlap assumption is provided.
  
contains
  subroutine scops(npoints,nlev,ncol,rngs,cc,conv,overlap,frac_out,ncolprint)
    INTEGER :: npoints,    &    ! Number of model points in the horizontal
               nlev,       &    ! Number of model levels in column
               ncol,       &    ! Number of subcolumns
               overlap          ! Overlap type (1=max, 2=rand, 3=max/rand)
    type(rng_state),dimension(npoints) :: rngs            
    INTEGER, parameter :: huge32 = 2147483647
    INTEGER, parameter :: i2_16  = 65536
    INTEGER :: i,j,ilev,ibox,ncolprint,ilev2

    REAL(WP), dimension(npoints,nlev) ::  &
         cc,         &    ! Input cloud cover in each model level (fraction)
                          ! NOTE:  This is the HORIZONTAL area of each
                          !        grid box covered by clouds
         conv,       &    ! Input convective cloud cover in each model level (fraction)
                          ! NOTE:  This is the HORIZONTAL area of each
                          !        grid box covered by convective clouds
         tca              ! Total cloud cover in each model level (fraction)
                          ! with extra layer of zeroes on top
                          ! in this version this just contains the values input
                          ! from cc but with an extra level
    REAL(WP),intent(inout), dimension(npoints,ncol,nlev) :: &
         frac_out         ! Boxes gridbox divided up into equivalent of BOX in 
                          ! original version, but indexed by column then row, rather than
                          ! by row then column
    REAL(WP), dimension(npoints,ncol) :: &
         threshold,   &   ! pointer to position in gridbox
         maxocc,      &   ! Flag for max overlapped conv cld
         maxosc,      &   ! Flag for max overlapped strat cld
         boxpos,      &   ! ordered pointer to position in gridbox
         threshold_min    ! minimum value to define range in with new threshold is chosen.
    REAL(WP), dimension(npoints) :: &
         ran              ! vector of random numbers

    ! Test for valid input overlap assumption
    if (overlap .ne. 1 .and. overlap .ne. 2 .and. overlap .ne. 3) then
       overlap=default_overlap
       call errorMessage('ERROR(scops): Invalid overlap assumption provided. Using default overlap assumption (max/ran)')
    endif

    boxpos = spread(([(i, i=1,ncol)]-0.5)/ncol,1,npoints)
    
    ! #######################################################################
    ! Initialize working variables
    ! #######################################################################
    
    ! Initialize frac_out to zero
    frac_out(1:npoints,1:ncol,1:nlev)=0.0     
    
    ! Assign 2d tca array using 1d input array cc
    tca(1:npoints,1:nlev) = cc(1:npoints,1:nlev)
    
    if (ncolprint.ne.0) then
       write (6,'(a)') 'frac_out_pp_rev:'
       do j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)
       enddo
       write (6,'(a)') 'ncol:'
       write (6,'(I3)') ncol
    endif
    if (ncolprint.ne.0) then
       write (6,'(a)') 'last_frac_pp:'
       do j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') (tca(j,1))
       enddo
    endif
    
    ! #######################################################################
    ! ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
    ! frac_out is the array that contains the information 
    ! where 0 is no cloud, 1 is a stratiform cloud and 2 is a
    ! convective cloud
    ! #######################################################################
    
    ! Loop over vertical levels
    DO ilev = 1,nlev
       
       ! Initialise threshold
       IF (ilev.eq.1) then
          ! If max overlap 
          IF (overlap.eq.1) then
             ! Select pixels spread evenly across the gridbox
             threshold(1:npoints,1:ncol)=boxpos(1:npoints,1:ncol)
          ELSE
             DO ibox=1,ncol
                !include 'congvec.f90'
                ran(1:npoints) = get_rng(RNGS)
                ! select random pixels from the non-convective
                ! part the gridbox ( some will be converted into
                ! convective pixels below )
                threshold(1:npoints,ibox) = conv(1:npoints,ilev)+(1-conv(1:npoints,ilev))*ran(npoints)
             enddo
          ENDIF
          IF (ncolprint.ne.0) then
             write (6,'(a)') 'threshold_nsf2:'
             do j=1,npoints,1000
                write(6,'(a10)') 'j='
                write(6,'(8I10)') j
                write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
             enddo
          ENDIF
       ENDIF
       
       IF (ncolprint.ne.0) then
          write (6,'(a)') 'ilev:'
          write (6,'(I2)') ilev
       ENDIF
       
       DO ibox=1,ncol
          ! All versions
          !maxocc(1:npoints,ibox) = merge(1,0,boxpos(1:npoints,ibox) .le. conv(1:npoints,ilev))
          !maxocc(1:npoints,ibox) = merge(1,0, conv(1:npoints,ilev) .gt. boxpos(1:npoints,ibox))
          do j=1,npoints
             if (boxpos(j,ibox).le.conv(j,ilev)) then
                maxocc(j,ibox) = 1
             else
                maxocc(j,ibox) = 0
             end if
          enddo
          
          ! Max overlap
          if (overlap.eq.1) then 
             threshold_min(1:npoints,ibox) = conv(1:npoints,ilev)
             maxosc(1:npoints,ibox)        = 1               
          endif
          
          ! Random overlap
          if (overlap.eq.2) then 
             threshold_min(1:npoints,ibox) = conv(1:npoints,ilev)
             maxosc(1:npoints,ibox)        = 0
          endif
          ! Max/Random overlap
          if (overlap.eq.3) then 
             ! DS2014 START: The bounds on tca are not valid when ilev=1.
             !threshold_min(1:npoints,ibox) = max(conv(1:npoints,ilev),min(tca(1:npoints,ilev-1),tca(1:npoints,ilev)))
             !maxosc(1:npoints,ibox) = merge(1,0,threshold(1:npoints,ibox) .lt. &
             !     min(tca(1:npoints,ilev-1),tca(1:npoints,ilev)) .and. &
             !     (threshold(1:npoints,ibox).gt.conv(1:npoints,ilev)))
             if (ilev .ne. 1) then
                threshold_min(1:npoints,ibox) = max(conv(1:npoints,ilev),min(tca(1:npoints,ilev-1),tca(1:npoints,ilev)))
                maxosc(1:npoints,ibox) = merge(1,0,threshold(1:npoints,ibox) .lt. &
                     min(tca(1:npoints,ilev-1),tca(1:npoints,ilev)) .and. &
                     (threshold(1:npoints,ibox).gt.conv(1:npoints,ilev)))
             else
                threshold_min(1:npoints,ibox) = max(conv(1:npoints,ilev),min(0._wp,tca(1:npoints,ilev)))
                maxosc(1:npoints,ibox) = merge(1,0,threshold(1:npoints,ibox) .lt. &
                     min(0._wp,tca(1:npoints,ilev)) .and. &
                     (threshold(1:npoints,ibox).gt.conv(1:npoints,ilev)))
             endif
          endif
          
          ! Reset threshold 
          !include 'congvec.f90'
          ran(1:npoints) = get_rng(RNGS)
          
          threshold(1:npoints,ibox)= maxocc(1:npoints,ibox)*(boxpos(1:npoints,ibox)) +            &
               (1-maxocc(1:npoints,ibox))*((maxosc(1:npoints,ibox))*(threshold(1:npoints,ibox)) + &
               (1-maxosc(1:npoints,ibox))*(threshold_min(1:npoints,ibox)+                         &
               (1-threshold_min(1:npoints,ibox))*ran(1:npoints)))
          
          ! Fill frac_out with 1's where tca is greater than the threshold
          frac_out(1:npoints,ibox,ilev) = merge(1,0,tca(1:npoints,ilev).gt.threshold(1:npoints,ibox))
          
          ! Code to partition boxes into startiform and convective parts goes here
          where(threshold(1:npoints,ibox).le.conv(1:npoints,ilev) .and. conv(1:npoints,ilev).gt.0.) frac_out(1:npoints,ibox,ilev)=2
       ENDDO ! ibox
       
       
       ! Set last_frac to tca at this level, so as to be tca from last level next time round
       if (ncolprint.ne.0) then
          do j=1,npoints ,1000
             write(6,'(a10)') 'j='
             write(6,'(8I10)') j
             write (6,'(a)') 'last_frac:'
             write (6,'(8f5.2)') (tca(j,ilev))
             write (6,'(a)') 'conv:'
             write (6,'(8f5.2)') (conv(j,ilev),ibox=1,ncolprint)
             write (6,'(a)') 'max_overlap_cc:'
             write (6,'(8f5.2)') (maxocc(j,ibox),ibox=1,ncolprint)
             write (6,'(a)') 'max_overlap_sc:'
             write (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)
             write (6,'(a)') 'threshold_min_nsf2:'
             write (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)
             write (6,'(a)') 'threshold_nsf2:'
             write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
             write (6,'(a)') 'frac_out_pp_rev:'
             write (6,'(8f5.2)') ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
          enddo
       endif
       
    enddo ! Loop over nlev
    
    ! END
  end subroutine scops
end module mod_scops
