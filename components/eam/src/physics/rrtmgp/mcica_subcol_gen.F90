module mcica_subcol_gen

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
! 
! Purpose: Create McICA stochastic arrays for cloud optical properties.
! Input cloud optical properties directly: cloud optical depth, single
! scattering albedo and asymmetry parameter.  Output will be stochastic
! arrays of these variables.  (longwave scattering is not yet available)
!
! Original code: From RRTMG based on Raisanen et al., QJRMS, 2004.
! 
! Uses the KISS random number generator.
!
! Overlap assumption: maximum-random.
! 
!----------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use shr_RandNum_mod,  only: ShrKissRandGen
use assertions

implicit none
private
save

public :: mcica_subcol_mask

!========================================================================================
contains
!========================================================================================

subroutine mcica_subcol_mask( &
   ngpt, ncol, pver, changeseed, &
   pmid, cldfrac, iscloudy)

   ! Arrays use CAM vertical index convention: index increases from top to bottom.
   ! This index ordering is assumed in the maximum-random overlap algorithm which starts
   ! at the top of a column and marches down, with each layer depending on the state
   ! of the layer above it.
   !
   ! For GCM mode, changeseed must be offset between LW and SW by at least the
   ! number of subcolumns

   ! arguments
   integer, intent(in) :: ngpt                      ! number of subcolumns (g-point intervals)
   integer, intent(in) :: ncol                      ! number of columns
   integer, intent(in) :: pver                      ! number of levels
   integer, intent(in) :: changeseed                ! if the subcolumn generator is called multiple times, 
                                                    ! permute the seed between each call.
   real(r8), intent(in) :: pmid(ncol,pver)         ! layer pressures (Pa)
   real(r8), intent(in) :: cldfrac(ncol,pver)      ! layer cloud fraction
   logical, intent(out) :: iscloudy(ngpt,ncol,pver) ! flag that says whether a gridbox is cloudy

   ! Local vars
   integer :: i, isubcol, k, n
   real(r8), parameter :: cldmin = 1.0e-80_r8  ! min cloud fraction
   real(r8) :: cldf(ncol,pver)      ! cloud fraction clipped to cldmin

   type(ShrKissRandGen) :: kiss_gen  ! KISS RNG object
   integer  :: kiss_seed(ncol,4)
   real(r8) :: rand_num_1d(ncol,1)   ! random number (kissvec)
   real(r8) :: rand_num(ncol,pver)   ! random number (kissvec)

   real(r8) :: cdf(ngpt,ncol,pver)   ! random numbers
   character(len=128) :: sub_name = 'mcica_subcol_gen_mask'
   !------------------------------------------------------------------------------------------ 

   ! Make sure inputs are TOA to surface
   call assert(pmid(1,1) < pmid(1,2), trim(sub_name) // ': inputs not TOA to SFC')

   ! clip cloud fraction
   cldf(:,:) = cldfrac(:ncol,:)
   where (cldf(:,:) < cldmin)
      cldf(:,:) = 0._r8
   end where

   ! Create a seed that depends on the state of the columns.
   ! Use pmid from bottom four layers. 
   do i = 1, ncol
      kiss_seed(i,1) = (pmid(i,pver)   - int(pmid(i,pver)))    * 1000000000
      kiss_seed(i,2) = (pmid(i,pver-1) - int(pmid(i,pver-1)))  * 1000000000
      kiss_seed(i,3) = (pmid(i,pver-2) - int(pmid(i,pver-2)))  * 1000000000
      kiss_seed(i,4) = (pmid(i,pver-3) - int(pmid(i,pver-3)))  * 1000000000
   end do

   ! create the RNG object
   kiss_gen = ShrKissRandGen(kiss_seed)

   ! Advance randum number generator by changeseed values
   do i = 1, changeseed
      call kiss_gen%random(rand_num_1d)
   end do

   ! Generate random numbers in each subcolumn at every level
   do isubcol = 1,ngpt
      call kiss_gen%random(rand_num)
      cdf(isubcol,:,:) = rand_num(:,:)
   enddo

   ! Maximum-Random overlap
   ! i) pick a random number for top layer.
   ! ii) walk down the column: 
   !    - if the layer above is cloudy, use the same random number as in the layer above
   !    - if the layer above is clear, use a new random number 
   do k = 2, pver
      do i = 1, ncol
         do isubcol = 1, ngpt
            if (cdf(isubcol,i,k-1) > 1._r8 - cldf(i,k-1) ) then
               cdf(isubcol,i,k) = cdf(isubcol,i,k-1) 
            else
               cdf(isubcol,i,k) = cdf(isubcol,i,k) * (1._r8 - cldf(i,k-1))
            end if
         end do
      end do
   end do
 
   ! Set binary cloudy flag
   do k = 1, pver
      iscloudy(:,:,k) = (cdf(:,:,k) >= 1._r8 - spread(cldf(:,k), dim=1, nCopies=ngpt) )
   end do

   call kiss_gen%finalize()

end subroutine mcica_subcol_mask

end module mcica_subcol_gen
