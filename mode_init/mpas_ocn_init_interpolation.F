! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_vertical_grids
!
!> \brief MPAS ocean vertical grid generator
!> \author Xylar Asay-Davis
!> \date   10/30/2015
!> \details
!>  This module contains the routines for generating
!>  vertical grids.
!
!-----------------------------------------------------------------------
module ocn_init_interpolation

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer

   use ocn_constants

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_interpolation_linear_vert, &
             ocn_init_interpolation_nearest_horiz, &
             ocn_init_interpolation_bilinear_horiz

   interface ocn_init_interpolation_nearest_horiz
     module procedure ocn_init_interpolation_nearest_horiz_2D
     module procedure ocn_init_interpolation_nearest_horiz_3D
   end interface

   interface ocn_init_interpolation_bilinear_horiz
     module procedure ocn_init_interpolation_bilinear_horiz_2D
     module procedure ocn_init_interpolation_bilinear_horiz_3D
   end interface

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_interpolation_linear_vert
!
!> \brief   linearly interpolate a field in the vertical
!> \author  Xylar Asay-Davis
!> \date    10/12/2015
!> \details
!>  Perform vertical linear interpolation of a field from a reference field
!>  with data located at inZ to new locations outZ.  Bu default, out-of-range values of
!>  outZ are clamped to the nearest in-range value.  That is, wherever
!>  outZ > inZ(1), outField = inField(1); wherever outZ < inZ(inNVertLevels),
!>  outFiled = inField(inNVertLevels).  If the optional extrapolate argument
!>  is present and set to .true., linear extrapolation is perfored outside the
!>  bounds of inZ.

!-----------------------------------------------------------------------

   subroutine ocn_init_interpolation_linear_vert(inZ, inField, inNVertLevels, outZ, outField, outNVertLevels, extrapolate)!{{{

   !--------------------------------------------------------------------

      integer, intent(in) :: inNVertLevels, outNVertLevels
      real (kind=RKIND), dimension(inNVertLevels), intent(in) :: inZ, inField
      real (kind=RKIND), dimension(outNVertLevels), intent(in) :: outZ
      real (kind=RKIND), dimension(outNVertLevels), intent(out) :: outField
      logical, optional, intent(in) :: extrapolate

      ! Define variable pointers
      integer :: outK, inK
      real (kind=RKIND) :: z, frac

      logical :: doExtrapolate

      if(outNVertLevels <= 0) return

      ! handle single level of input data as a special case
      if(inNVertLevels == 1) then
        outField(1:outNVertLevels) = inField(1)
        return
      end if

      doExtrapolate = .false.
      if(present(extrapolate)) then
        doExtrapolate = extrapolate
      end if

      if(doExtrapolate) then
        do outK=1,outNVertLevels
          z = outZ(outK)
          if(z >= inZ(1)) then
            inK = 1
          else if(z <= inZ(inNVertLevels)) then
            inK = inNVertLevels-1
          else
            do inK = 1, inNVertLevels-1
              if(z >= inZ(inK+1)) exit
            end do
          end if
          ! frac can be outside [0,1) if we're extrapolating
          frac = (z - inZ(inK))/(inZ(inK+1) - inZ(inK))
          outField(outK) = (1.0_RKIND - frac)*inField(inK) + frac*inField(inK+1)
        end do
      else ! not extrapolating
        do outK=1,outNVertLevels
          z = outZ(outK)
          if(z >= inZ(1)) then
            inK = 1
            frac = 0.0_RKIND
          else if(z <= inZ(inNVertLevels)) then
            inK = inNVertLevels-1
            frac = 1.0_RKIND
          else
            do inK = 1, inNVertLevels-1
              if(z >= inZ(inK+1)) exit
            end do
            ! frac should always be inside [0,1)
            frac = (z - inZ(inK))/(inZ(inK+1) - inZ(inK))
          end if
          outField(outK) = (1.0_RKIND - frac)*inField(inK) + frac*inField(inK+1)
        end do
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_init_interpolation_linear_vert!}}}

!***********************************************************************
!
!  routine ocn_init_interpolation_nearest_horiz_2D
!
!> \brief   nearest-neighbor interpolation in horiz.
!> \author  Xylar Asay-Davis
!> \date    10/30/2015
!> \details
!>  Perform horizontal nearest-neighbor interpolation of a field from
!>  values on a logically rectangular grid.

!-----------------------------------------------------------------------

   subroutine ocn_init_interpolation_nearest_horiz_2D(inX, inY, inField, inNx, inNy, &
                                                   outX, outY, outField, outN, &
                                                   inXPeriod, inYPeriod)!{{{

   !--------------------------------------------------------------------

      integer, intent(in) :: inNx, inNy, outN
      real (kind=RKIND), dimension(inNx), intent(in) :: inX
      real (kind=RKIND), dimension(inNy), intent(in) :: inY
      real (kind=RKIND), dimension(inNx,inNy), intent(in) :: inField
      real (kind=RKIND), dimension(outN), intent(in) :: outX, outY
      real (kind=RKIND), dimension(outN), intent(out) :: outField
      real (kind=RKIND), intent(in), optional :: inXPeriod, inYPeriod

      ! Define variable pointers
      integer :: outIndex, xSearch, ySearch, searchIdx
      real (kind=RKIND) :: currentX, currentY, minDist, dist

      do outIndex = 1, outN
         currentX = outX(outIndex)
         if(present(inXPeriod)) then
            ! put currentX in the range of [inX(1),inX(1)+inXPeriod)
            currentX = mod(currentX-inX(1),inXPeriod) + inX(1)
         end if
         currentY = outY(outIndex)
         if(present(inYPeriod)) then
            ! put currentY in the range of [inY(1),inY(1)+inYPeriod)
            currentY = mod(currentY-inY(1),inYPeriod) + inY(1)
         end if

         xSearch = 1
         minDist = 1e34
         do searchIdx = 1, inNx
            dist = abs(currentX - inX(searchIdx))
            if (dist < minDist) then
               minDist = dist
               xSearch = searchIdx
            end if
         end do

         ySearch = 1
         minDist = 1e34
         do searchIdx = 1, inNy
            dist = abs(currentY - inY(searchIdx))
            if (dist < minDist) then
               minDist = dist
               ySearch = searchIdx
            end if
         end do

         outField(outIndex) = inField(xSearch, ySearch)

      end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_interpolation_nearest_horiz_2D!}}}

!***********************************************************************
!
!  routine ocn_init_interpolation_nearest_horiz_3D
!
!> \brief   nearest-neighbor interpolation in horiz.
!> \author  Xylar Asay-Davis
!> \date    10/30/2015
!> \details
!>  Perform horizontal nearest-neighbor interpolation of a field from
!>  values on a logically rectangular grid.

!-----------------------------------------------------------------------

   subroutine ocn_init_interpolation_nearest_horiz_3D(inX, inY, inField, inNx, inNy, &
                                                   outX, outY, outField, outN, &
                                                   inXPeriod, inYPeriod)!{{{

   !--------------------------------------------------------------------

      integer, intent(in) :: inNx, inNy, outN
      real (kind=RKIND), dimension(inNx), intent(in) :: inX
      real (kind=RKIND), dimension(inNy), intent(in) :: inY
      real (kind=RKIND), dimension(:,:,:), intent(in) :: inField
      real (kind=RKIND), dimension(outN), intent(in) :: outX, outY
      real (kind=RKIND), dimension(:,:), intent(out) :: outField
      real (kind=RKIND), intent(in), optional :: inXPeriod, inYPeriod

      ! Define variable pointers
      integer :: outIndex, xSearch, ySearch, searchIdx
      real (kind=RKIND) :: currentX, currentY, minDist, dist

      do outIndex = 1, outN
         currentX = outX(outIndex)
         if(present(inXPeriod)) then
            ! put currentX in the range of [inX(1),inX(1)+inXPeriod)
            currentX = mod(currentX-inX(1),inXPeriod) + inX(1)
         end if
         currentY = outY(outIndex)
         if(present(inYPeriod)) then
            ! put currentY in the range of [inY(1),inY(1)+inYPeriod)
            currentY = mod(currentY-inY(1),inYPeriod) + inY(1)
         end if

         xSearch = 1
         minDist = 1e34
         do searchIdx = 1, inNx
            dist = abs(currentX - inX(searchIdx))
            if (dist < minDist) then
               minDist = dist
               xSearch = searchIdx
            end if
         end do

         ySearch = 1
         minDist = 1e34
         do searchIdx = 1, inNy
            dist = abs(currentY - inY(searchIdx))
            if (dist < minDist) then
               minDist = dist
               ySearch = searchIdx
            end if
         end do

         outField(:,outIndex) = inField(xSearch, ySearch,:)

      end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_interpolation_nearest_horiz_3D!}}}

!***********************************************************************
!
!  routine ocn_init_interpolation_bilinear_horiz_2D
!
!> \brief   bilinear interpolation in horiz.
!> \author  Xylar Asay-Davis
!> \date    10/30/2015
!> \details
!>  Perform horizontal bilinear interpolation of a field from
!>  values on a logically rectangular grid. Optional parameters
!>  inXPeriod and inYPeriod are used to specify the period of the
!>  input grid.  If either or both are omitted, the grid is not
!>  treated as being periodic in either or both dimensions.

!-----------------------------------------------------------------------

   subroutine ocn_init_interpolation_bilinear_horiz_2D(inX, inY, inField, inNx, inNy, &
                                                   outX, outY, outField, outN, &
                                                   inXPeriod, inYPeriod, extrapX, extrapY)!{{{

   !--------------------------------------------------------------------

      integer, intent(in) :: inNx, inNy, outN
      real (kind=RKIND), dimension(inNx), intent(in) :: inX
      real (kind=RKIND), dimension(inNy), intent(in) :: inY
      real (kind=RKIND), dimension(inNx,inNy), intent(in) :: inField
      real (kind=RKIND), dimension(outN), intent(in) :: outX, outY
      real (kind=RKIND), dimension(outN), intent(out) :: outField
      real (kind=RKIND), intent(in), optional :: inXPeriod, inYPeriod
      logical, intent(in), optional :: extrapX, extrapY

      ! Define variable pointers
      integer :: outIndex, xInd1, xInd2, yInd1, yInd2, k
      real (kind=RKIND) :: x, y, xFrac, yFrac

      do outIndex = 1, outN
         x = outX(outIndex)
         y = outY(outIndex)

         if(present(inXPeriod)) then
            call getLinearCoeffs(x, inX, inNx, xInd1, xInd2, xFrac, inXPeriod)
         else
            call getLinearCoeffs(x, inX, inNx, xInd1, xInd2, xFrac)
         end if

         ! if we're not extrapolating, limit xFrac
         if(present(extrapX)) then
            if(.not. extrapX) then
               xFrac = min(1.0_RKIND,max(0.0_RKIND,xFrac))
            end if
         else
            ! by default, we don't extrapolate
            xFrac = min(1.0_RKIND,max(0.0_RKIND,xFrac))
         end if

         if(present(inYPeriod)) then
            call getLinearCoeffs(y, inY, inNy, yInd1, yInd2, yFrac, inYPeriod)
         else
            call getLinearCoeffs(y, inY, inNy, yInd1, yInd2, yFrac)
         end if

         ! if we're not extrapolating, limit yFrac
         if(present(extrapY)) then
            if(.not. extrapY) then
               yFrac = min(1.0_RKIND,max(0.0_RKIND,yFrac))
            end if
         else
            ! by default, we don't extrapolate
            yFrac = min(1.0_RKIND,max(0.0_RKIND,yFrac))
         end if

         outField(outIndex) =  &
            (1.0_RKIND-xFrac)*(1.0_RKIND-yFrac)*inField(xInd1,yInd1) &
            + xFrac*(1.0_RKIND-yFrac)*inField(xInd2,yInd1) &
            + (1.0_RKIND-xFrac)*yFrac*inField(xInd1,yInd2) &
            + xFrac*yFrac*inField(xInd2,yInd2)
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_interpolation_bilinear_horiz_2D!}}}

!***********************************************************************
!
!  routine ocn_init_interpolation_bilinear_horiz_3D
!
!> \brief   bilinear interpolation in horiz.
!> \author  Xylar Asay-Davis
!> \date    10/30/2015
!> \details
!>  Perform horizontal bilinear interpolation of a field from
!>  values on a logically rectangular grid. Optional parameters
!>  inXPeriod and inYPeriod are used to specify the period of the
!>  input grid.  If either or both are omitted, the grid is not
!>  treated as being periodic in either or both dimensions.

!-----------------------------------------------------------------------

   subroutine ocn_init_interpolation_bilinear_horiz_3D(inX, inY, inField, inNx, inNy, &
                                                   outX, outY, outField, outN, &
                                                   inXPeriod, inYPeriod, extrapX, extrapY)!{{{

   !--------------------------------------------------------------------

      integer, intent(in) :: inNx, inNy, outN
      real (kind=RKIND), dimension(inNx), intent(in) :: inX
      real (kind=RKIND), dimension(inNy), intent(in) :: inY
      real (kind=RKIND), dimension(:,:,:), intent(in) :: inField
      real (kind=RKIND), dimension(outN), intent(in) :: outX, outY
      real (kind=RKIND), dimension(:,:), intent(out) :: outField
      real (kind=RKIND), intent(in), optional :: inXPeriod, inYPeriod
      logical, intent(in), optional :: extrapX, extrapY

      ! Define variable pointers
      integer :: outIndex, xInd1, xInd2, yInd1, yInd2, k
      real (kind=RKIND) :: x, y, xFrac, yFrac

      do outIndex = 1, outN
         x = outX(outIndex)
         y = outY(outIndex)

         if(present(inXPeriod)) then
            call getLinearCoeffs(x, inX, inNx, xInd1, xInd2, xFrac, inXPeriod)
         else
            call getLinearCoeffs(x, inX, inNx, xInd1, xInd2, xFrac)
         end if

         ! if we're not extrapolating, limit xFrac
         if(present(extrapX)) then
            if(.not. extrapX) then
               xFrac = min(1.0_RKIND,max(0.0_RKIND,xFrac))
            end if
         else
            ! by default, we don't extrapolate
            xFrac = min(1.0_RKIND,max(0.0_RKIND,xFrac))
         end if

         if(present(inYPeriod)) then
            call getLinearCoeffs(y, inY, inNy, yInd1, yInd2, yFrac, inYPeriod)
         else
            call getLinearCoeffs(y, inY, inNy, yInd1, yInd2, yFrac)
         end if

         ! if we're not extrapolating, limit yFrac
         if(present(extrapY)) then
            if(.not. extrapY) then
               yFrac = min(1.0_RKIND,max(0.0_RKIND,yFrac))
            end if
         else
            ! by default, we don't extrapolate
            yFrac = min(1.0_RKIND,max(0.0_RKIND,yFrac))
         end if

         outField(:,outIndex) =  &
            (1.0_RKIND-xFrac)*(1.0_RKIND-yFrac)*inField(xInd1,yInd1,:) &
            + xFrac*(1.0_RKIND-yFrac)*inField(xInd2,yInd1,:) &
            + (1.0_RKIND-xFrac)*yFrac*inField(xInd1,yInd2,:) &
            + xFrac*yFrac*inField(xInd2,yInd2,:)
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_interpolation_bilinear_horiz_3D!}}}

!***********************************************************************
!
!  routine getLinearCoeffs
!
!> \brief   compute coefficients for linear interpolation
!> \author  Xylar Asay-Davis
!> \date    10/30/2015
!> \details
!>  Given a point and an array of locations, returns the indices
!>  the point is bounded by, as a fraction, where the point lies
!>  between the nearest two points in array of locations.
!>  Optional parameter period is used to indicate that the array
!>  of locations is periodic with the given period.

!-----------------------------------------------------------------------

   subroutine getLinearCoeffs(xValue, xArray, nx, index1, index2, frac, period)!{{{

      integer, intent(in) :: nx
      real (kind=RKIND), intent(in) :: xValue
      real (kind=RKIND), dimension(nx), intent(in) :: xArray
      integer, intent(out) :: index1, index2
      real (kind=RKIND), intent(out) :: frac
      real (kind=RKIND), intent(in), optional :: period

      integer :: xIndex
      real (kind=RKIND) :: x

      x = xValue

      if(present(period)) then
         ! Set up bilinear interpolation indices in x, watching for periodic boundary
         ! shift x to be within the range of [xArray(1),xArray(1)+period)
         x = modulo(x-xArray(1),period) + xArray(1)
         if (x >= xArray(nx)) then
            ! at the periodic boundary so treat as special case
            index1 = nx
            index2 = 1
            frac = (x-xArray(index1))/(xArray(index2)+period-xArray(index1))
         else
            do xIndex = 1, nx-1
               if (x .le. xArray(xIndex+1)) then
                  index1 = xIndex
                  index2 = xIndex+1
                  frac = (x-xArray(index1))/(xArray(index2)-xArray(index1))
                  exit
               end if
            end do
         end if
      else
         ! not periodic
         index1 = nx-1
         do xIndex = 1, nx-2
            if (x .le. xArray(xIndex+1)) then
               index1 = xIndex
               exit
            end if
         end do
         index2 = index1+1
         frac = (x-xArray(index1))/(xArray(index2)-xArray(index1))
      end if

   end subroutine getLinearCoeffs!}}}


!***********************************************************************

end module ocn_init_interpolation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
