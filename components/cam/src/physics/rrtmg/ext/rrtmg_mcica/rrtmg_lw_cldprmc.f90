!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_cldprmc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.5 $
!     created:   $Date: 2007/09/07 21:40:50 $
!
      module rrtmg_lw_cldprmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb 
      use parrrtm, only : ngptlw
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3
      use rrlw_wvn, only: ngb
      use rrlw_vsn, only: hvrclc, hnamclc

      implicit none

      contains

! ------------------------------------------------------------------------------
      subroutine cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, &
                         ciwpmc, clwpmc, reicmc, dgesmc, relqmc, ncbands, taucmc)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

! ------- Input -------

      integer, intent(in) :: nlayers                    ! total number of layers
      integer, intent(in) :: inflag                     ! see definitions
      integer, intent(in) :: iceflag                    ! see definitions
      integer, intent(in) :: liqflag                    ! see definitions

      real(kind=r8), intent(in) :: cldfmc(:,:)          ! cloud fraction [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: ciwpmc(:,:)          ! cloud ice water path [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: clwpmc(:,:)          ! cloud liquid water path [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: relqmc(:)            ! liquid particle effective radius (microns)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: reicmc(:)            ! ice particle effective radius (microns)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: dgesmc(:)            ! ice particle generalized effective size (microns)
                                                        !    Dimensions: (nlayers)

! ------- Output -------

      integer, intent(out) :: ncbands                   ! number of cloud spectral bands
      real(kind=r8), intent(inout) :: taucmc(:,:)       ! cloud optical depth [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)

! ------- Local -------

      integer :: lay                            ! Layer index
      integer :: ib                             ! spectral band index
      integer :: ig                             ! g-point interval index
      integer :: index 

      real(kind=r8) :: abscoice(ngptlw)         ! ice absorption coefficients
      real(kind=r8) :: abscoliq(ngptlw)         ! liquid absorption coefficients
      real(kind=r8) :: cwp                      ! cloud water path
      real(kind=r8) :: radice                   ! cloud ice effective radius (microns)
      real(kind=r8) :: dgeice                   ! cloud ice generalized effective size
      real(kind=r8) :: factor                   ! 
      real(kind=r8) :: fint                     ! 
      real(kind=r8) :: radliq                   ! cloud liquid droplet radius (microns)
      real(kind=r8), parameter :: eps = 1.e-6_r8     ! epsilon
      real(kind=r8), parameter :: cldmin = 1.e-80_r8 ! minimum value for cloud quantities

! ------- Definitions -------

!     Explanation of the method for each value of INFLAG.  Values of
!     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
!     INFLAG = 2 does distinguish between liquid and ice clouds, and
!     requires further user input to specify the method to be used to 
!     compute the aborption due to each.
!     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
!                  optical depth are input.  
!     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
!                  water path (g/m2) are input.  The (gray) cloud optical 
!                  depth is computed as in CAM3.
!     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
!                  water path (g/m2), and cloud ice fraction are input.
!       ICEFLAG = 0:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in CAM3.
!       ICEFLAG = 1:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in 
!                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The 
!                     spectral regions in this work have been matched with
!                     the spectral bands in RRTM to as great an extent 
!                     as possible:  
!                     E&C 1      IB = 5      RRTM bands 9-16
!                     E&C 2      IB = 4      RRTM bands 6-8
!                     E&C 3      IB = 3      RRTM bands 3-5
!                     E&C 4      IB = 2      RRTM band 2
!                     E&C 5      IB = 1      RRTM band 1
!       ICEFLAG = 2:  The ice effective radius (microns) is input and the
!                     optical properties due to ice clouds are computed from
!                     the optical properties stored in the RT code,
!                     STREAMER v3.0 (Reference: Key. J., Streamer 
!                     User's Guide, Cooperative Institute for
!                     Meteorological Satellite Studies, 2001, 96 pp.).
!                     Valid range of values for re are between 5.0 and
!                     131.0 micron.
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!       LIQFLAG = 0:  The optical depths due to water clouds are computed as
!                     in CAM3.
!       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
!                     and the optical depths due to water clouds are computed 
!                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                     The values for absorption coefficients appropriate for
!                     the spectral bands in RRTM have been obtained for a 
!                     range of effective radii by an averaging procedure 
!                     based on the work of J. Pinto (private communication).
!                     Linear interpolation is used to get the absorption 
!                     coefficients for the input effective radius.

      hvrclc = '$Revision: 1.5 $'

      ncbands = 1

! This initialization is done in rrtmg_lw_subcol.F90.
!      do lay = 1, nlayers
!         do ig = 1, ngptlw
!            taucmc(ig,lay) = 0.0_r8
!         enddo
!      enddo

! Main layer loop
      do lay = 1, nlayers

        do ig = 1, ngptlw
          cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
          if (cldfmc(ig,lay) .ge. cldmin .and. &
             (cwp .ge. cldmin .or. taucmc(ig,lay) .ge. cldmin)) then

! Ice clouds and water clouds combined.
            if (inflag .eq. 0) then
! Cloud optical depth already defined in taucmc, return to main program
               return

            elseif(inflag .eq. 1) then 
                stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'
!               cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
!               taucmc(ig,lay) = abscld1 * cwp

! Separate treatement of ice clouds and water clouds.
            elseif(inflag .eq. 2) then
               radice = reicmc(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwpmc(ig,lay) .eq. 0.0_r8) then
                  abscoice(ig) = 0.0_r8

               elseif (iceflag .eq. 0) then
                  if (radice .lt. 10.0_r8) stop 'ICE RADIUS TOO SMALL'
                  abscoice(ig) = absice0(1) + absice0(2)/radice

               elseif (iceflag .eq. 1) then
! mji - turn off limits to mimic CAM3
!                  if (radice .lt. 13.0_r8 .or. radice .gt. 130._r8) stop &
!                      'ICE RADIUS OUT OF BOUNDS'
                  ncbands = 5
                  ib = ngb(ig)
                  abscoice(ig) = absice1(1,ib) + absice1(2,ib)/radice

! For iceflag=2 option, combine with iceflag=0 option to handle out of bounds 
! particle sizes.
! Use iceflag=2 option for ice particle effective radii from 5.0 and 131.0 microns
! and use iceflag=0 option for ice particles greater than 131.0 microns.
! *** NOTE: Transition between two methods has not been smoothed. 

               elseif (iceflag .eq. 2) then
                  if (radice .lt. 5.0_r8) stop 'ICE RADIUS OUT OF BOUNDS'
                  if (radice .ge. 5.0_r8 .and. radice .le. 131._r8) then
                     ncbands = 16
                     factor = (radice - 2._r8)/3._r8
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice(ig) = &
                         absice2(index,ib) + fint * &
                         (absice2(index+1,ib) - (absice2(index,ib))) 
                  elseif (radice .gt. 131._r8) then
                     abscoice(ig) = absice0(1) + absice0(2)/radice
                  endif
               
! For iceflag=3 option, combine with iceflag=0 option to handle large particle sizes.
! Use iceflag=3 option for ice particle effective radii from 3.2 and 91.0 microns
! (generalized effective size, dge, from 5 to 140 microns), and use iceflag=0 option
! for ice particle effective radii greater than 91.0 microns (dge = 140 microns).
! *** NOTE: Fu parameterization requires particle size in generalized effective size.
! *** NOTE: Transition between two methods has not been smoothed. 

               elseif (iceflag .eq. 3) then
                  dgeice = dgesmc(lay)
                  if (dgeice .lt. 5.0_r8) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                  if (dgeice .ge. 5.0_r8 .and. dgeice .le. 140._r8) then
                     ncbands = 16
                     factor = (dgeice - 2._r8)/3._r8
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice(ig) = &
                         absice3(index,ib) + fint * &
                         (absice3(index+1,ib) - (absice3(index,ib)))
                  elseif (dgeice .gt. 140._r8) then
                     abscoice(ig) = absice0(1) + absice0(2)/radice
                  endif
   
               endif
                  
! Calculation of absorption coefficients due to water clouds.
               if (clwpmc(ig,lay) .eq. 0.0_r8) then
                  abscoliq(ig) = 0.0_r8

               elseif (liqflag .eq. 0) then
                   abscoliq(ig) = absliq0

               elseif (liqflag .eq. 1) then
                  radliq = relqmc(lay)
                  if (radliq .lt. 1.5_r8 .or. radliq .gt. 60._r8) stop &
                       'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                  index = radliq - 1.5_r8
                  if (index .eq. 58) index = 57
                  if (index .eq. 0) index = 1
                  fint = radliq - 1.5_r8 - index
                  ib = ngb(ig)
                  abscoliq(ig) = &
                        absliq1(index,ib) + fint * &
                        (absliq1(index+1,ib) - (absliq1(index,ib)))
               endif

               taucmc(ig,lay) = ciwpmc(ig,lay) * abscoice(ig) + &
                                clwpmc(ig,lay) * abscoliq(ig)

            endif
         endif
         enddo
      enddo

      end subroutine cldprmc

      end module rrtmg_lw_cldprmc
