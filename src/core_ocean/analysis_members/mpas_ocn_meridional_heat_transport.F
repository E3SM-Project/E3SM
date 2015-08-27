! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_meridional_heat_transport
!
!> \brief MPAS ocean analysis core member: meridional_heat_transport
!> \author Mark Petersen
!> \date   March 2014
!> \details
!>  MPAS ocean analysis core member: meridional_heat_transport
!>  Compute zonal means of selected variables
!
!-----------------------------------------------------------------------

module ocn_meridional_heat_transport

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

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

   public :: ocn_init_meridional_heat_transport, &
             ocn_compute_meridional_heat_transport, &
             ocn_restart_meridional_heat_transport, &
             ocn_finalize_meridional_heat_transport

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: nMerHeatTransBinsUsed

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_meridional_heat_transport
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_meridional_heat_transport(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meridionalHeatTransportAMPool
      type (mpas_pool_type), pointer :: meshPool

      integer ::  iBin
      integer, pointer ::  nMerHeatTransBins

      real (kind=RKIND) :: binWidth
      ! These are array size 1 because mpas_dmpar_min_real_array calls require arrays.
      real (kind=RKIND), dimension(1) :: minBin, maxBin, minBinDomain, maxBinDomain
      real (kind=RKIND), dimension(:), pointer ::  binBoundaryMerHeatTrans, binVariable

      integer, pointer :: config_AM_meridionalHeatTransport_num_bins
      real (kind=RKIND), pointer :: config_AM_meridionalHeatTransport_min_bin, config_AM_meridionalHeatTransport_max_bin

      logical, pointer :: on_a_sphere

      dminfo = domain % dminfo

      err = 0

      minBin =  1.0e34_RKIND
      maxBin = -1.0e34_RKIND

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nMerHeatTransBins', nMerHeatTransBins)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'meridionalHeatTransportAM', meridionalHeatTransportAMPool)

      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_num_bins', config_AM_meridionalHeatTransport_num_bins)
      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_min_bin', config_AM_meridionalHeatTransport_min_bin)
      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_max_bin', config_AM_meridionalHeatTransport_max_bin)

      nMerHeatTransBinsUsed = config_AM_meridionalHeatTransport_num_bins

      call mpas_pool_get_array(meridionalHeatTransportAMPool, 'binBoundaryMerHeatTrans', binBoundaryMerHeatTrans)

      ! Find min and max values of binning variable.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

         ! Bin by latitude on a sphere, by yCell otherwise.
         if (on_a_sphere) then
            call mpas_pool_get_array(meshPool, 'latCell', binVariable)
         else
            call mpas_pool_get_array(meshPool, 'yCell', binVariable)
         end if

         minBin = min(minBin, minval(binVariable) )
         maxBin = max(maxBin, maxval(binVariable) )

         block => block % next
      end do

      call mpas_dmpar_min_real_array(dminfo, 1, minBin, minBinDomain)
      call mpas_dmpar_max_real_array(dminfo, 1, maxBin, maxBinDomain)

      ! Set up bins. 
      binBoundaryMerHeatTrans = -1.0e34_RKIND

      ! Change min and max bin bounds to configuration settings, if applicable.
      if (config_AM_meridionalHeatTransport_min_bin > -1.0e33_RKIND) then
         minBinDomain(1) = config_AM_meridionalHeatTransport_min_bin
      else
         ! use measured min value, but decrease slightly to include least value.
         minBinDomain(1) = minBinDomain(1) - 1.0e-10_RKIND * abs(minBinDomain(1))
      end if

      if (config_AM_meridionalHeatTransport_max_bin > -1.0e33_RKIND) then
         maxBinDomain(1) = config_AM_meridionalHeatTransport_max_bin
      else
         ! use measured max value, but increase slightly to include max value.
         maxBinDomain(1) = maxBinDomain(1) + 1.0e-10_RKIND * abs(maxBinDomain(1))
      end if

      binBoundaryMerHeatTrans(1) = minBinDomain(1)
      binWidth = (maxBinDomain(1) - minBinDomain(1)) / nMerHeatTransBinsUsed

      do iBin = 2, nMerHeatTransBinsUsed
         binBoundaryMerHeatTrans(iBin) = binBoundaryMerHeatTrans(iBin-1) + binWidth
      end do
      binBoundaryMerHeatTrans(nMerHeatTransBinsUsed+1) = binBoundaryMerHeatTrans(nMerHeatTransBinsUsed) + binWidth

   end subroutine ocn_init_meridional_heat_transport!}}}

!***********************************************************************
!
!  routine ocn_compute_meridional_heat_transport
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_meridional_heat_transport(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meridionalHeatTransportAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      integer :: iTracer, k, iCell, kMax, i, iEdge
      integer :: iBin, iField, nMerHeatTransVariables
      integer, pointer :: nCellsSolve, nVertLevels, nMerHeatTransBins, indexTemperature
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgeSignOnCell, cellsOnEdge, edgesOnCell

      real (kind=RKIND) :: div_huT
      real (kind=RKIND), dimension(:), pointer ::  areaCell, binVariable, binBoundaryMerHeatTrans, dvEdge
      real (kind=RKIND), dimension(:), pointer ::  meridionalHeatTransportLat
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, normalTransportVelocity
      real (kind=RKIND), dimension(:,:), pointer :: meridionalHeatTransportLatZ
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:,:), allocatable :: mht_meridional_integral
      real (kind=RKIND), dimension(:,:,:), allocatable :: sumMerHeatTrans, totalSumMerHeatTrans

      logical, pointer :: on_a_sphere

      err = 0
      dminfo = domain % dminfo

      call mpas_pool_get_subpool(domain % blocklist % structs, 'meridionalHeatTransportAM', meridionalHeatTransportAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      nMerHeatTransVariables = 1

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nMerHeatTransBins', nMerHeatTransBins)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meridionalHeatTransportAMPool, 'binBoundaryMerHeatTrans', binBoundaryMerHeatTrans)

      allocate(sumMerHeatTrans(nMerHeatTransVariables,nVertLevels,nMerHeatTransBinsUsed))
      allocate(totalSumMerHeatTrans(nMerHeatTransVariables,nVertLevels,nMerHeatTransBinsUsed))
      allocate(mht_meridional_integral(nVertLevels,nMerHeatTransBinsUsed))

      sumMerHeatTrans = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)

         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

         ! Bin by latitude on a sphere, by yCell otherwise.
         if (on_a_sphere) then
            call mpas_pool_get_array(meshPool, 'latCell', binVariable)
         else
            call mpas_pool_get_array(meshPool, 'yCell', binVariable)
         end if

         do iCell = 1, nCellsSolve
            kMax = maxLevelCell(iCell)

            if (binVariable(iCell) .lt. binBoundaryMerHeatTrans(1)) cycle

            do iBin = 1, nMerHeatTransBinsUsed
               if (binVariable(iCell) .lt. binBoundaryMerHeatTrans(iBin+1) ) then

                  do k = 1, kMax

                     ! Compute divergence of huT, i.e. layerThicknessEdge * normalVelocity * temperature, at an edge
                     ! for meridional heat transport.  Here we use a centered difference to compute the temperature at
                     ! the edge, which is an approximation to the actual edge temperature used in the horizontal
                     ! advection scheme (for example, FCT).  We expect that the error in this approximation is small.
                     ! Here we do not divide by the area, as one normally does in a divergence calculation, so that
                     ! div_huT is weighted by area here.
                     iField = 1
                     div_huT = 0.0_RKIND
                     do i = 1, nEdgesOnCell(iCell)
                        iEdge = edgesOnCell(i, iCell)
                        div_huT = div_huT - layerThicknessEdge(k, iEdge) * normalTransportVelocity(k, iEdge) &
                             * 0.5_RKIND * (tracers(indexTemperature,k,cellsOnEdge(1,iEdge)) + tracers(indexTemperature,k,cellsOnEdge(2,iEdge))) &
                             * edgeSignOnCell(i, iCell) * dvEdge(iEdge)
                     end do
                     sumMerHeatTrans(iField,k,iBin) = sumMerHeatTrans(iField,k,iBin) + div_huT

                  end do
                  exit

               endif
            end do

         end do

         block => block % next
      end do

      ! mpi summation over all processors
      ! Note the input and output arrays are of the same dimension, so summation is
      ! over the domain decompositon (by processor), not over an array index.
      call mpas_dmpar_sum_real_array(dminfo, nVertLevels*nMerHeatTransBinsUsed*nMerHeatTransVariables, sumMerHeatTrans, totalSumMerHeatTrans)

      ! Even though these variables do not include an index that is decomposed amongst
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_dimension(block % dimensions, 'nMerHeatTransBins', nMerHeatTransBins)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         call mpas_pool_get_subpool(block % structs, 'meridionalHeatTransportAM', meridionalHeatTransportAMPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)

         call mpas_pool_get_array(meridionalHeatTransportAMPool, 'meridionalHeatTransportLat', meridionalHeatTransportLat)
         call mpas_pool_get_array(meridionalHeatTransportAMPool, 'meridionalHeatTransportLatZ', meridionalHeatTransportLatZ)

         do iBin = 1, nMerHeatTransBinsUsed
            do k = 1, nVertLevels

               ! MHT = sum ( div(huT) A ) * rho c_p, in PW
               ! where the sum is over each latitude bin
               ! Here we simply multiply by (rho c_p) and convert to PW:
               iField = 1
               mht_meridional_integral(k,iBin) = totalSumMerHeatTrans(iField,k,iBin)*rho_sw*cp_sw*1.0e-15_RKIND

            end do
         end do

         ! Compute integral of ( sum ( div(huT) A ) * rho c_p ) from southernmost latitude to bin boundary.
         ! Note that mht_meridional_integral is indexed by bin, spanning 1:nMerHeatTransBinsUsed, while
         ! meridionalHeatTransportLatZ (second dimension) is indexed by bin boundary, spanning 1:nMerHeatTransBinsUsed+1
         meridionalHeatTransportLatZ(:,1) = 0.0_RKIND
         do iBin = 2, nMerHeatTransBinsUsed+1
            meridionalHeatTransportLatZ(:,iBin) = meridionalHeatTransportLatZ(:,iBin-1) + mht_meridional_integral(:,iBin-1)
         end do

         ! meridionalHeatTransportLatZ is a function of depth.  Sum in vertical to get
         ! meridionalHeatTransportLat, a single value for each latitude bin boundary.
         ! meridionalHeatTransportLat is indexed by bin boundary, spanning 1:nMerHeatTransBinsUsed+1
         do iBin = 1, nMerHeatTransBinsUsed+1
            meridionalHeatTransportLat(iBin) = sum(meridionalHeatTransportLatZ(:,iBin))
         end do


         block => block % next
      end do

      deallocate(sumMerHeatTrans)
      deallocate(totalSumMerHeatTrans)
      deallocate(mht_meridional_integral)

   end subroutine ocn_compute_meridional_heat_transport!}}}

!***********************************************************************
!
!  routine ocn_restart_meridional_heat_transport
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_meridional_heat_transport(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_meridional_heat_transport!}}}

!***********************************************************************
!
!  routine ocn_finalize_meridional_heat_transport
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_meridional_heat_transport(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_meridional_heat_transport!}}}

end module ocn_meridional_heat_transport

! vim: foldmethod=marker
