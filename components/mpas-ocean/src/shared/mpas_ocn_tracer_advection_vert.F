! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_vert
!
!> \brief MPAS standard tracer advection
!> \author Doug Jacobsen, Phil Jones
!> \date   03/09/12, separated and revised 7/2021
!> \details
!>  This module contains routines for computing vertical advection
!>  fluxes for use in overall tracer advection tendency routines.
!
!-----------------------------------------------------------------------

module ocn_tracer_advection_vert

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_log

   use ocn_config
   use ocn_mesh

   implicit none
   private
   save

   ! private module variables
   real (kind=RKIND) ::  &
      coef3rdOrder        !< high-order coefficient

   integer :: vertOrder   !< choice of order for vertical advection
   integer, parameter :: &! enumerator for supported vertical adv order
      vertOrder2=2,      &!< 2nd order
      vertOrder3=3,      &!< 3rd order
      vertOrder4=4        !< 4th order

   ! public method interfaces
   public :: ocn_tracer_advection_vert_flx, &
             ocn_tracer_advection_vert_init

!***********************************************************************

   contains

!***********************************************************************
!
!  routine ocn_tracer_advection_vert_flx
!
!> \brief Computes vertical tracer advection fluxes
!> \author Doug Jacobsen
!> \date   03/09/12, separated and revised 7/2021
!> \details
!>  This routine computes the high-order vertical tracer advection
!>  based on requested order.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_vert_flx(tracer, w, layerThick, &
                                            vertFlx)

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertFlx         !< [out] high-order vertical advection flux

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tracer          !< [in] current tracer values

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         w,             &!< [in] Vertical velocity
         layerThick      !< [in] Layer thickness to use for advection

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell,            &! cell address
         k, kmin, kmax      ! current, min, max vertical layer indices

      real (kind=RKIND) :: &
         verticalWeightK,  &! weights for vertical advection
         verticalWeightKm1  ! weights for vertical advection

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Compute the high order vertical flux.

      ! Initialize return flux to zero
      ! Also ensures that top and bottom fluxes are zero

#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) present(vertFlx)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCellsAll
      do k=1,nVertLevels+1
         vertFlx(k,iCell) = 0.0_RKIND
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! Compute flux for interior layers at high order based on
      ! user-requested order

      select case (vertOrder)
      case (vertOrder4)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(minLevelCell, maxLevelCell, w, &
         !$acc            vertFlx, tracer) &
         !$acc    private(k, kmin, kmax)
#else
         !$omp do schedule(runtime) &
         !$omp    private(k, kmin, kmax)
#endif
         do iCell = 1, nCellsAll
            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)
            do k=kmin+2,kmax-1
               vertFlx(k, iCell) = w(k,iCell)* &
                                   ( 7.0_RKIND*(tracer(k  ,iCell)+ &
                                                tracer(k-1,iCell)) &
                                             - (tracer(k+1,iCell)+ &
                                                tracer(k-2,iCell)))/&
                                                12.0_RKIND
            end do ! vertical loop
         end do ! iCell Loop
#ifndef MPAS_OPENACC
         !$omp end do
#endif

      case (vertOrder3)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(minLevelCell, maxLevelCell, w, &
         !$acc            vertFlx, tracer) &
         !$acc    private(k, kmin, kmax)
#else
         !$omp do schedule(runtime) &
         !$omp    private(k, kmin, kmax)
#endif
         do iCell = 1, nCellsAll
            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)
            do k=kmin+2,kmax-1
               vertFlx(k, iCell) = (w(k,iCell)* &
                                   (7.0_RKIND * (tracer(k  ,iCell)+ &
                                                 tracer(k-1,iCell)) - &
                                                (tracer(k+1,iCell)+ &
                                                 tracer(k-2,iCell))) - &
                                coef3rdOrder*abs(w(k,iCell))* &
                                               ((tracer(k+1,iCell)- &
                                                 tracer(k-2,iCell)) - &
                                      3.0_RKIND*(tracer(k  ,iCell)- &
                                                 tracer(k-1,iCell))))/ &
                                                 12.0_RKIND
            end do ! vertical loop
         end do ! iCell Loop
#ifndef MPAS_OPENACC
         !$omp end do
#endif

      case (vertOrder2)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(minLevelCell, maxLevelCell, w, &
         !$acc            layerThick, vertFlx, tracer) &
         !$acc    private(k, kmin, kmax, &
         !$acc            verticalWeightK, verticalWeightKm1)
#else
         !$omp do schedule(runtime) &
         !$omp    private(k, kmin, kmax, &
         !$omp            verticalWeightK, verticalWeightKm1)
#endif
         do iCell = 1, nCellsAll
            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)
            do k=kmin+2,kmax-1
               verticalWeightK   = layerThick(k-1,iCell) / &
                                  (layerThick(k  ,iCell) + &
                                   layerThick(k-1,iCell))
               verticalWeightKm1 = layerThick(k  ,iCell) / &
                                  (layerThick(k  ,iCell) + &
                                   layerThick(k-1,iCell))
               vertFlx(k,iCell) = w(k, iCell) * &
                           (verticalWeightK  *tracer(k  ,iCell) + &
                            verticalWeightKm1*tracer(k-1,iCell))
            end do ! vertical loop
         end do ! iCell Loop
#ifndef MPAS_OPENACC
         !$omp end do
#endif

      end select ! vertOrder

      ! Now take care of the edge cases, reducing order for
      ! layers near the top and bottom where high-order forms
      ! can not be computed.

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(minLevelCell, maxLevelCell, w, &
      !$acc            layerThick, vertFlx, tracer) &
      !$acc    private(k, kmin, kmax, &
      !$acc            verticalWeightK, verticalWeightKm1)
#else
      !$omp do schedule(runtime) &
      !$omp    private(k, kmin, kmax, &
      !$omp            verticalWeightK, verticalWeightKm1)
#endif
      do iCell = 1, nCellsAll
         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)
         ! at top, flux is zero (already initialized)
         ! at next-to-top (kmin+1), reduce to 2nd order
         !   but avoid case where 0 or 1 active layer (kmax <= kmin)
         if (kmax > kmin) then
            k = kmin+1
            verticalWeightK   = layerThick(k-1,iCell) / &
                               (layerThick(k  ,iCell) + &
                                layerThick(k-1,iCell))
            verticalWeightKm1 = layerThick(k  ,iCell) / &
                               (layerThick(k  ,iCell) + &
                                layerThick(k-1,iCell))
            vertFlx(k,iCell) = w(k,iCell)&
                          *(verticalWeightK  *tracer(k  ,iCell) + &
                            verticalWeightKm1*tracer(k-1,iCell))
            ! Deepest active level also at 2nd order
            k = kmax
            verticalWeightK   = layerThick(k-1,iCell) / &
                               (layerThick(k  ,iCell) + &
                                layerThick(k-1,iCell))
            verticalWeightKm1 = layerThick(k  ,iCell) / &
                               (layerThick(k  ,iCell) + &
                                layerThick(k-1,iCell))
            vertFlx(k,iCell) = w(k,iCell)* &
                          (verticalWeightK  *tracer(k  ,iCell) + &
                           verticalWeightKm1*tracer(k-1,iCell))
         end if ! kmax > kmin
      end do ! iCell Loop
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_advection_vert_flx !}}}

!***********************************************************************
!
!  routine ocn_tracer_advection_vert_init
!
!> \brief Initializes vertical tracer advection flux computation.
!> \author Doug Jacobsen, Phil Jones
!> \date   03/09/12, updated May 2019
!> \details
!>  This routine initializes constants and choices for the computation
!>  of vertical tracer advection fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_vert_init(err) !{{{

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      ! end of preamble
      !----------------
      ! begin code

      err = 0 ! set error code to success

      ! set choice of vertical advection order
      select case (config_vert_tracer_adv_flux_order)
      case (2)
         vertOrder = vertOrder2
         coef3rdOrder = 0.0_RKIND
      case (3)
         vertOrder = vertOrder3
         coef3rdOrder = config_coef_3rd_order
      case (4)
         vertOrder = vertOrder4
         coef3rdOrder = 0.0_RKIND
      case default
         vertOrder = vertOrder2
         coef3rdOrder = 0.0_RKIND
         call mpas_log_write( &
         'Invalid value for vertical advection order, defaulting to 2nd order',&
         MPAS_LOG_WARN)
      end select ! vertAdvOrder

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_advection_vert_init!}}}

!***********************************************************************

end module ocn_tracer_advection_vert

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
