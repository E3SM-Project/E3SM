      subroutine scops(npoints,nlev,ncol,seed,cc,conv,
     &                 overlap,frac_out,ncolprint)


! *****************************COPYRIGHT****************************
! (c) British Crown Copyright 2009, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the
! following conditions are met:
! 
!     * Redistributions of source code must retain the above 
!       copyright  notice, this list of conditions and the following 
!       disclaimer.
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its 
!       contributors may be used to endorse or promote products
!       derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
! 
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************

      implicit none

      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns


      INTEGER overlap         !  overlap type
                              !  1=max
                              !  2=rand
                              !  3=max/rand
      REAL cc(npoints,nlev)
                  !  input cloud cover in each model level (fraction)
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds

      REAL conv(npoints,nlev)
                  !  input convective cloud cover in each model
                  !   level (fraction)
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by convective clouds

      INTEGER i,j,ilev,ibox,ncolprint,ilev2

      REAL frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column


      INTEGER seed(npoints)
      !  seed values for marsaglia  random number generator
      !  It is recommended that the seed is set
      !  to a different value for each model
      !  gridbox it is called on, as it is
      !  possible that the choice of the same
      !  seed value every time may introduce some
      !  statistical bias in the results, particularly
      !  for low values of NCOL.

      REAL tca(npoints,0:nlev) ! total cloud cover in each model level (fraction)
                                        ! with extra layer of zeroes on top
                                        ! in this version this just contains the values input
                                        ! from cc but with an extra level

      REAL threshold(npoints,ncol)   ! pointer to position in gridbox
      REAL maxocc(npoints,ncol)      ! Flag for max overlapped conv cld
      REAL maxosc(npoints,ncol)      ! Flag for max overlapped strat cld

      REAL boxpos(npoints,ncol)      ! ordered pointer to position in gridbox

      REAL threshold_min(npoints,ncol) ! minimum value to define range in with new threshold
                                        ! is chosen

      REAL ran(npoints)                 ! vector of random numbers

      INTEGER i8,irand,i2_16,huge32,overflow_32  ! variables for RNG
      logical bigendian
      PARAMETER(huge32=2147483647)
      PARAMETER(i8 = selected_int_kind(13))

      bigendian = (transfer(1_i8, 1) == 0)
      i2_16=65536

      do ibox=1,ncol
        do j=1,npoints 
        boxpos(j,ibox)=(ibox-.5)/ncol
        enddo
      enddo

!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!     Initialised frac_out to zero

      do ilev=1,nlev
        do ibox=1,ncol
          do j=1,npoints
          frac_out(j,ibox,ilev)=0.0
          enddo
        enddo
      enddo

!     assign 2d tca array using 1d input array cc

      do j=1,npoints
        tca(j,0)=0
      enddo

      do ilev=1,nlev
        do j=1,npoints
          tca(j,ilev)=cc(j,ilev)
        enddo
      enddo

      if (ncolprint.ne.0) then
        write (6,'(a)') 'frac_out_pp_rev:'
          do j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') 
     &     ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)

          enddo
        write (6,'(a)') 'ncol:'
        write (6,'(I3)') ncol
      endif
      if (ncolprint.ne.0) then
        write (6,'(a)') 'last_frac_pp:'
          do j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') (tca(j,0))
          enddo
      endif

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!     frac_out is the array that contains the information 
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud
      
      !loop over vertical levels
      DO 200 ilev = 1,nlev
                                  
!     Initialise threshold

        IF (ilev.eq.1) then
          ! If max overlap 
          IF (overlap.eq.1) then
            ! select pixels spread evenly
            ! across the gridbox
              DO ibox=1,ncol
                do j=1,npoints
                  threshold(j,ibox)=boxpos(j,ibox)
                enddo
              enddo
          ELSE
              DO ibox=1,ncol
                include 'congvec.inc'
                ! select random pixels from the non-convective
                ! part the gridbox ( some will be converted into
                ! convective pixels below )
                do j=1,npoints
                  threshold(j,ibox)=
     &            conv(j,ilev)+(1-conv(j,ilev))*ran(j)
                enddo
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
          do j=1,npoints
            if (boxpos(j,ibox).le.conv(j,ilev)) then
              maxocc(j,ibox) = 1
            else
              maxocc(j,ibox) = 0
            end if
          enddo

          ! Max overlap
          if (overlap.eq.1) then 
            do j=1,npoints
              threshold_min(j,ibox)=conv(j,ilev)
              maxosc(j,ibox)=1
            enddo
          endif

          ! Random overlap
          if (overlap.eq.2) then 
            do j=1,npoints
              threshold_min(j,ibox)=conv(j,ilev)
              maxosc(j,ibox)=0
            enddo
          endif

          ! Max/Random overlap
          if (overlap.eq.3) then 
            do j=1,npoints
              threshold_min(j,ibox)=max(conv(j,ilev),
     &          min(tca(j,ilev-1),tca(j,ilev)))
              if (threshold(j,ibox)
     &          .lt.min(tca(j,ilev-1),tca(j,ilev))
     &          .and.(threshold(j,ibox).gt.conv(j,ilev))) then
                   maxosc(j,ibox)= 1
              else
                   maxosc(j,ibox)= 0
              end if
            enddo
          endif
    
          ! Reset threshold 

          include 'congvec.inc'

          do j=1,npoints
            threshold(j,ibox)=
              !if max overlapped conv cloud
     &        maxocc(j,ibox) * (                                       
     &            boxpos(j,ibox)                                               
     &        ) +                                                      
              !else
     &        (1-maxocc(j,ibox)) * (                                   
                  !if max overlapped strat cloud
     &            (maxosc(j,ibox)) * (                                 
                      !threshold=boxpos
     &                threshold(j,ibox)                                        
     &            ) +                                                  
                  !else
     &            (1-maxosc(j,ibox)) * (                               
                      !threshold_min=random[thrmin,1]
     &                threshold_min(j,ibox)+
     &                  (1-threshold_min(j,ibox))*ran(j)  
     &           ) 
     &        )
          enddo

        ENDDO ! ibox

!          Fill frac_out with 1's where tca is greater than the threshold

           DO ibox=1,ncol
             do j=1,npoints 
               if (tca(j,ilev).gt.threshold(j,ibox)) then
               frac_out(j,ibox,ilev)=1
               else
               frac_out(j,ibox,ilev)=0
               end if               
             enddo
           ENDDO

!         Code to partition boxes into startiform and convective parts
!         goes here

           DO ibox=1,ncol
             do j=1,npoints 
                if (threshold(j,ibox).le.conv(j,ilev)) then
                    ! = 2 IF threshold le conv(j)
                    frac_out(j,ibox,ilev) = 2 
                else
                    ! = the same IF NOT threshold le conv(j) 
                    frac_out(j,ibox,ilev) = frac_out(j,ibox,ilev)
                end if
             enddo
           ENDDO

!         Set last_frac to tca at this level, so as to be tca 
!         from last level next time round

          if (ncolprint.ne.0) then

            do j=1,npoints ,1000
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write (6,'(a)') 'last_frac:'
            write (6,'(8f5.2)') (tca(j,ilev-1))
    
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
            write (6,'(8f5.2)') 
     &       ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
          enddo
          endif

200   CONTINUE    !loop over nlev


      end

