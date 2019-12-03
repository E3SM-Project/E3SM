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
! Jul 2008 - A. Bodas-Salcedo - Added capability of producing outputs in standard grid
! Oct 2008 - J.-L. Dufresne   - Bug fixed. Assignment of Npoints,Nlevels,Nhydro,Ncolumns 
!                               in COSP_STATS
! Oct 2008 - H. Chepfer       - Added PARASOL reflectance arguments
! Jun 2010 - T. Yokohata, T. Nishimura and K. Ogochi - Added NEC SXs optimisations
! Jan 2013 - G. Cesana        - Added betaperp and temperature arguments 
!                             - Added phase 3D/3Dtemperature/Map output variables in diag_lidar 
! May 2015 - D. Swales        - Modified for cosp2.0 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_STATS
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF,R_GROUND
  IMPLICIT NONE
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_CHANGE_VERTICAL_GRID(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,Nglevels,newgrid_bot,newgrid_top,r,log_units)
   implicit none
   ! Input arguments
   integer,intent(in) :: Npoints  !# of grid points
   integer,intent(in) :: Nlevels  !# of levels
   integer,intent(in) :: Ncolumns !# of columns
   real(wp),dimension(Npoints,Nlevels),intent(in) :: zfull ! Height at model levels [m] (Bottom of model layer)
   real(wp),dimension(Npoints,Nlevels),intent(in) :: zhalf ! Height at half model levels [m] (Bottom of model layer)
   real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: y     ! Variable to be changed to a different grid
   integer,intent(in) :: Nglevels  !# levels in the new grid
   real(wp),dimension(Nglevels),intent(in) :: newgrid_bot ! Lower boundary of new levels  [m]
   real(wp),dimension(Nglevels),intent(in) :: newgrid_top ! Upper boundary of new levels  [m]
   logical,optional,intent(in) :: log_units ! log units, need to convert to linear units
   ! Output
   real(wp),dimension(Npoints,Ncolumns,Nglevels),intent(out) :: r ! Variable on new grid

   ! Local variables
   integer :: i,j,k
   logical :: lunits
   integer :: l
   real(wp) :: w ! Weight
   real(wp) :: dbb, dtb, dbt, dtt ! Distances between edges of both grids
   integer :: Nw  ! Number of weights
   real(wp) :: wt  ! Sum of weights
   real(wp),dimension(Nlevels) :: oldgrid_bot,oldgrid_top ! Lower and upper boundaries of model grid
   real(wp) :: yp ! Local copy of y at a particular point.
              ! This allows for change of units.

   lunits=.false.
   if (present(log_units)) lunits=log_units

   r = 0._wp

   do i=1,Npoints
     ! Calculate tops and bottoms of new and old grids
     oldgrid_bot = zhalf(i,:)
     oldgrid_top(1:Nlevels-1) = oldgrid_bot(2:Nlevels)
     oldgrid_top(Nlevels) = zfull(i,Nlevels) +  zfull(i,Nlevels) - zhalf(i,Nlevels) ! Top level symmetric
     l = 0 ! Index of level in the old grid
     ! Loop over levels in the new grid
     do k = 1,Nglevels
       Nw = 0 ! Number of weigths
       wt = 0._wp ! Sum of weights
       ! Loop over levels in the old grid and accumulate total for weighted average
       do
         l = l + 1
         w = 0.0 ! Initialise weight to 0
         ! Distances between edges of both grids
         dbb = oldgrid_bot(l) - newgrid_bot(k)
         dtb = oldgrid_top(l) - newgrid_bot(k)
         dbt = oldgrid_bot(l) - newgrid_top(k)
         dtt = oldgrid_top(l) - newgrid_top(k)
         if (dbt >= 0.0) exit ! Do next level in the new grid
         if (dtb > 0.0) then
           if (dbb <= 0.0) then
             if (dtt <= 0) then
               w = dtb
             else
               w = newgrid_top(k) - newgrid_bot(k)
             endif
           else
             if (dtt <= 0) then
               w = oldgrid_top(l) - oldgrid_bot(l)
             else
               w = -dbt
             endif
           endif
           ! If layers overlap (w/=0), then accumulate
           if (w /= 0.0) then
             Nw = Nw + 1
             wt = wt + w
             do j=1,Ncolumns
               if (lunits) then
                 if (y(i,j,l) /= R_UNDEF) then
                   yp = 10._wp**(y(i,j,l)/10._wp)
                 else
                   yp = 0._wp
                 endif
               else
                 yp = y(i,j,l)
               endif
               r(i,j,k) = r(i,j,k) + w*yp
             enddo
           endif
         endif
       enddo
       l = l - 2
       if (l < 1) l = 0
       ! Calculate average in new grid
       if (Nw > 0) then
         do j=1,Ncolumns
           r(i,j,k) = r(i,j,k)/wt
         enddo
       endif
     enddo
   enddo

   ! Set points under surface to R_UNDEF, and change to dBZ if necessary
   do k=1,Nglevels
     do j=1,Ncolumns
       do i=1,Npoints
         if (newgrid_top(k) > zhalf(i,1)) then ! Level above model bottom level
           if (lunits) then
             if (r(i,j,k) <= 0.0) then
               r(i,j,k) = R_UNDEF
             else
               r(i,j,k) = 10._wp*log10(r(i,j,k))
             endif
           endif
         else ! Level below surface
           r(i,j,k) = R_GROUND
         endif
       enddo
     enddo
   enddo

END SUBROUTINE COSP_CHANGE_VERTICAL_GRID

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------- SUBROUTINE COSP_LIDAR_ONLY_CLOUD -----------------
  ! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
  ! All rights reserved.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_LIDAR_ONLY_CLOUD(Npoints,Ncolumns,Nlevels,beta_tot, &
                                   beta_mol,Ze_tot,lidar_only_freq_cloud,tcc)
    ! Inputs
    integer,intent(in) :: &
         Npoints,       & ! Number of horizontal gridpoints
         Ncolumns,      & ! Number of subcolumns
         Nlevels          ! Number of vertical layers
    real(wp),dimension(Npoints,Nlevels),intent(in) :: &
         beta_mol         ! Molecular backscatter
    real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         beta_tot,      & ! Total backscattered signal
         Ze_tot           ! Radar reflectivity
    ! Outputs
    real(wp),dimension(Npoints,Nlevels),intent(out) :: &
         lidar_only_freq_cloud
    real(wp),dimension(Npoints),intent(out) ::&
         tcc
    
    ! local variables
    real(wp) :: sc_ratio
    real(wp),parameter :: &
         s_cld=5.0, &
         s_att=0.01
    integer :: flag_sat,flag_cld,pr,i,j
    
    lidar_only_freq_cloud = 0._wp
    tcc = 0._wp
    do pr=1,Npoints
       do i=1,Ncolumns
          flag_sat = 0
          flag_cld = 0
          do j=1,Nlevels
             sc_ratio = beta_tot(pr,i,j)/beta_mol(pr,j)
             if ((sc_ratio .le. s_att) .and. (flag_sat .eq. 0)) flag_sat = j
             if (Ze_tot(pr,i,j) .lt. -30.) then  !radar can't detect cloud
                if ( (sc_ratio .gt. s_cld) .or. (flag_sat .eq. j) ) then  !lidar sense cloud
                   lidar_only_freq_cloud(pr,j)=lidar_only_freq_cloud(pr,j)+1. !top->surf
                   flag_cld=1
                endif
             else  !radar sense cloud (z%Ze_tot(pr,i,j) .ge. -30.)
                flag_cld=1
             endif
          enddo !levels
          if (flag_cld .eq. 1) tcc(pr)=tcc(pr)+1._wp
       enddo !columns
    enddo !points
    lidar_only_freq_cloud=lidar_only_freq_cloud/Ncolumns
    tcc=tcc/Ncolumns
    
    ! Unit conversion
    where(lidar_only_freq_cloud /= R_UNDEF) &
            lidar_only_freq_cloud = lidar_only_freq_cloud*100._wp
    where(tcc /= R_UNDEF) tcc = tcc*100._wp
    
  END SUBROUTINE COSP_LIDAR_ONLY_CLOUD
  
  ! ######################################################################################
  ! FUNCTION hist1D
  ! ######################################################################################
  function hist1d(Npoints,var,nbins,bins)
    ! Inputs
    integer,intent(in) :: &
         Npoints, & ! Number of points in input array
         Nbins      ! Number of bins for sorting
    real(wp),intent(in),dimension(Npoints) :: &
         var        ! Input variable to be sorted
    real(wp),intent(in),dimension(Nbins+1) :: &
         bins       ! Histogram bins [lowest,binTops]  
    ! Outputs
    real(wp),dimension(Nbins) :: &
         hist1d     ! Output histogram      
    ! Local variables
    integer :: ij
    
    do ij=2,Nbins+1  
       hist1D(ij-1) = count(var .ge. bins(ij-1) .and. var .lt. bins(ij))
       if (count(var .eq. R_GROUND) .ge. 1) hist1D(ij-1)=R_UNDEF
    enddo
    
  end function hist1D
  
  ! ######################################################################################
  ! SUBROUTINE hist2D
  ! ######################################################################################
  subroutine hist2D(var1,var2,npts,bin1,nbin1,bin2,nbin2,jointHist)
    implicit none
    
    ! INPUTS
    integer, intent(in) :: &
         npts,  & ! Number of data points to be sorted
         nbin1, & ! Number of bins in histogram direction 1 
         nbin2    ! Number of bins in histogram direction 2
    real(wp),intent(in),dimension(npts) :: &
         var1,  & ! Variable 1 to be sorted into bins
         var2     ! variable 2 to be sorted into bins
    real(wp),intent(in),dimension(nbin1+1) :: &
         bin1     ! Histogram bin 1 boundaries
    real(wp),intent(in),dimension(nbin2+1) :: &
         bin2     ! Histogram bin 2 boundaries
    ! OUTPUTS
    real(wp),intent(out),dimension(nbin1,nbin2) :: &
         jointHist
    
    ! LOCAL VARIABLES
    integer :: ij,ik
    
    do ij=2,nbin1+1
       do ik=2,nbin2+1
          jointHist(ij-1,ik-1)=count(var1 .ge. bin1(ij-1) .and. var1 .lt. bin1(ij) .and. &
               var2 .ge. bin2(ik-1) .and. var2 .lt. bin2(ik))        
       enddo
    enddo
  end subroutine hist2D
END MODULE MOD_COSP_STATS
