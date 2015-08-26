! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
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
!     * Neither the name of the Lawrence Livermore National Security Limited Liability Corporation 
!       nor the names of its contributors may be used to endorse or promote products derived from 
!       this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      
      subroutine prec_scops(npoints,nlev,ncol,ls_p_rate,cv_p_rate,
     &                      frac_out,prec_frac)


      implicit none

      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns

      INTEGER i,j,ilev,ibox,cv_col
      
      REAL ls_p_rate(npoints,nlev),cv_p_rate(npoints,nlev)

      REAL frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column
                              !TOA to SURFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL prec_frac(npoints,ncol,nlev) ! 0 -> clear sky
                                        ! 1 -> LS precipitation
                                        ! 2 -> CONV precipitation
                    ! 3 -> both
                                        !TOA to SURFACE!!!!!!!!!!!!!!!!!!
                    
      INTEGER flag_ls, flag_cv
      INTEGER frac_out_ls(npoints,ncol),frac_out_cv(npoints,ncol) !flag variables for 
                       ! stratiform cloud and convective cloud in the vertical column

      cv_col = 0.05*ncol
      if (cv_col .eq. 0) cv_col=1
 
      do ilev=1,nlev
        do ibox=1,ncol
          do j=1,npoints 
            prec_frac(j,ibox,ilev) = 0
          enddo
        enddo
      enddo
      
      do j=1,npoints
       do ibox=1,ncol
        frac_out_ls(j,ibox)=0
        frac_out_cv(j,ibox)=0
        flag_ls=0
        flag_cv=0
        do ilev=1,nlev
          if (frac_out(j,ibox,ilev) .eq. 1) then 
            flag_ls=1
          endif
          if (frac_out(j,ibox,ilev) .eq. 2) then 
            flag_cv=1
          endif
        enddo !loop over nlev
        if (flag_ls .eq. 1) then
           frac_out_ls(j,ibox)=1
        endif
        if (flag_cv .eq. 1) then
           frac_out_cv(j,ibox)=1
        endif
       enddo  ! loop over ncol
      enddo ! loop over npoints

!      initialize the top layer      
       do j=1,npoints
        flag_ls=0
        flag_cv=0
    
        if (ls_p_rate(j,1) .gt. 0.) then 
            do ibox=1,ncol ! possibility ONE
                if (frac_out(j,ibox,1) .eq. 1) then 
                    prec_frac(j,ibox,1) = 1
                    flag_ls=1
                endif
            enddo ! loop over ncol
            if (flag_ls .eq. 0) then ! possibility THREE
                do ibox=1,ncol
                    if (frac_out(j,ibox,2) .eq. 1) then 
                        prec_frac(j,ibox,1) = 1
                        flag_ls=1
                    endif
                enddo ! loop over ncol
            endif
        if (flag_ls .eq. 0) then ! possibility Four
        do ibox=1,ncol
        if (frac_out_ls(j,ibox) .eq. 1) then 
            prec_frac(j,ibox,1) = 1
            flag_ls=1
        endif
        enddo ! loop over ncol
        endif
        if (flag_ls .eq. 0) then ! possibility Five
        do ibox=1,ncol
    !     prec_frac(j,1:ncol,1) = 1
        prec_frac(j,ibox,1) = 1
        enddo ! loop over ncol
            endif
        endif
       ! There is large scale precipitation
     
        if (cv_p_rate(j,1) .gt. 0.) then 
         do ibox=1,ncol ! possibility ONE
          if (frac_out(j,ibox,1) .eq. 2) then 
           if (prec_frac(j,ibox,1) .eq. 0) then
        prec_frac(j,ibox,1) = 2
       else
        prec_frac(j,ibox,1) = 3
       endif
       flag_cv=1
      endif
        enddo ! loop over ncol
        if (flag_cv .eq. 0) then ! possibility THREE
        do ibox=1,ncol
        if (frac_out(j,ibox,2) .eq. 2) then 
                if (prec_frac(j,ibox,1) .eq. 0) then
            prec_frac(j,ibox,1) = 2
            else
            prec_frac(j,ibox,1) = 3
            endif
            flag_cv=1
        endif
        enddo ! loop over ncol
        endif
        if (flag_cv .eq. 0) then ! possibility Four
        do ibox=1,ncol
        if (frac_out_cv(j,ibox) .eq. 1) then 
                if (prec_frac(j,ibox,1) .eq. 0) then
            prec_frac(j,ibox,1) = 2
            else
            prec_frac(j,ibox,1) = 3
            endif
            flag_cv=1
        endif
        enddo ! loop over ncol
        endif
        if (flag_cv .eq. 0) then  ! possibility Five
        do ibox=1,cv_col
                if (prec_frac(j,ibox,1) .eq. 0) then
            prec_frac(j,ibox,1) = 2
            else
            prec_frac(j,ibox,1) = 3
            endif 
        enddo !loop over cv_col
            endif 
        endif 
        ! There is convective precipitation
        
        enddo ! loop over npoints
!      end of initializing the top layer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     working on the levels from top to surface
      do ilev=2,nlev
       do j=1,npoints
        flag_ls=0
        flag_cv=0
    
        if (ls_p_rate(j,ilev) .gt. 0.) then 
         do ibox=1,ncol ! possibility ONE&TWO
          if ((frac_out(j,ibox,ilev) .eq. 1) .or. 
     &       ((prec_frac(j,ibox,ilev-1) .eq. 1) 
     &       .or. (prec_frac(j,ibox,ilev-1) .eq. 3))) then 
           prec_frac(j,ibox,ilev) = 1
           flag_ls=1
          endif
        enddo ! loop over ncol
        if ((flag_ls .eq. 0) .and. (ilev .lt. nlev)) then ! possibility THREE
        do ibox=1,ncol
        if (frac_out(j,ibox,ilev+1) .eq. 1) then 
            prec_frac(j,ibox,ilev) = 1
            flag_ls=1
        endif
        enddo ! loop over ncol
        endif
        if (flag_ls .eq. 0) then ! possibility Four
        do ibox=1,ncol
        if (frac_out_ls(j,ibox) .eq. 1) then 
            prec_frac(j,ibox,ilev) = 1
            flag_ls=1
        endif
        enddo ! loop over ncol
        endif
        if (flag_ls .eq. 0) then ! possibility Five
        do ibox=1,ncol
!     prec_frac(j,1:ncol,ilev) = 1
        prec_frac(j,ibox,ilev) = 1
        enddo ! loop over ncol
         endif
      endif ! There is large scale precipitation
    
        if (cv_p_rate(j,ilev) .gt. 0.) then 
         do ibox=1,ncol ! possibility ONE&TWO
          if ((frac_out(j,ibox,ilev) .eq. 2) .or. 
     &       ((prec_frac(j,ibox,ilev-1) .eq. 2) 
     &       .or. (prec_frac(j,ibox,ilev-1) .eq. 3))) then 
            if (prec_frac(j,ibox,ilev) .eq. 0) then
         prec_frac(j,ibox,ilev) = 2
        else
         prec_frac(j,ibox,ilev) = 3
        endif 
        flag_cv=1
        endif
       enddo ! loop over ncol
        if ((flag_cv .eq. 0) .and. (ilev .lt. nlev)) then ! possibility THREE
        do ibox=1,ncol
        if (frac_out(j,ibox,ilev+1) .eq. 2) then 
                if (prec_frac(j,ibox,ilev) .eq. 0) then
            prec_frac(j,ibox,ilev) = 2
            else
            prec_frac(j,ibox,ilev) = 3
            endif
            flag_cv=1
        endif
        enddo ! loop over ncol
        endif
        if (flag_cv .eq. 0) then ! possibility Four
        do ibox=1,ncol
        if (frac_out_cv(j,ibox) .eq. 1) then 
                if (prec_frac(j,ibox,ilev) .eq. 0) then
            prec_frac(j,ibox,ilev) = 2
            else
            prec_frac(j,ibox,ilev) = 3
            endif
            flag_cv=1
        endif
        enddo ! loop over ncol
        endif
        if (flag_cv .eq. 0) then  ! possibility Five 
        do ibox=1,cv_col
                if (prec_frac(j,ibox,ilev) .eq. 0) then
            prec_frac(j,ibox,ilev) = 2
            else
            prec_frac(j,ibox,ilev) = 3
            endif 
        enddo !loop over cv_col 
            endif 
        endif ! There is convective precipitation
    
        enddo ! loop over npoints
        enddo ! loop over nlev

      end

