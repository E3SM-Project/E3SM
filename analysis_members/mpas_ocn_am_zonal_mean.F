! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_am_zonal_mean
!
!> \brief MPAS ocean analysis core member: zonal_mean
!> \author Mark Petersen
!> \date   March 2014
!> \details
!>  MPAS ocean analysis core member: zonal_mean
!>  Compute zonal means of selected variables
!
!-----------------------------------------------------------------------

module ocn_am_zonal_mean

   use mpas_configure
   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar

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

   public :: ocn_am_setup_packages_zonal_mean, &
             ocn_am_init_zonal_mean, &
             ocn_am_compute_zonal_mean, &
             ocn_am_write_zonal_mean, &
             ocn_am_restart_zonal_mean, &
             ocn_am_finalize_zonal_mean

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: amZonalMeanTimer

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_am_setup_packages_zonal_mean
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_am_setup_packages_zonal_mean(err)!{{{

      use mpas_packages

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

      ! turn on package for this analysis member
      amZonalMeanActive = .true.

   end subroutine ocn_am_setup_packages_zonal_mean!}}}


!***********************************************************************
!
!  routine ocn_am_init_zonal_mean
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_init_zonal_mean(domain, err)!{{{

      use mpas_packages

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

      type (amZonalMean_type), pointer :: amZonalMean
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mesh_type), pointer :: mesh

      integer ::  nZonalMeanBins, iBin

      real (kind=RKIND) :: binWidth
      ! These are array size 1 because mpas_dmpar_min_real_array calls require arrays.
      real (kind=RKIND), dimension(1) :: minBin, maxBin, minBinDomain, maxBinDomain
      real (kind=RKIND), dimension(:), pointer ::  binCenterZonalMean, binBoundaryZonalMean, binVariable

      dminfo = domain % dminfo

      err = 0

      minBin =  1.0e34
      maxBin = -1.0e34

      nZonalMeanBins = domain % blocklist % mesh % nZonalMeanBins
      amZonalMean => domain % blocklist % amZonalMean

      binCenterZonalMean => amZonalMean % binCenterZonalMean % array
      binBoundaryZonalMean => amZonalMean % binBoundaryZonalMean % array

      block => domain % blocklist
      do while (associated(block))
         mesh => block % mesh

         binCenterZonalMean => amZonalMean % binCenterZonalMean % array
         binBoundaryZonalMean => amZonalMean % binBoundaryZonalMean % array

         if (config_ocn_am_zonal_mean_binning_variable=='latCell') then
            binVariable  => mesh % latCell % array
         elseif (config_ocn_am_zonal_mean_binning_variable=='yCell') then
            binVariable  => mesh % yCell % array
         else
            write (stderrUnit,'(a)') ' Error: incorrect choice of config_ocn_am_zonal_mean_binning_variable'
            call mpas_dmpar_abort(dminfo)
         end if 

         minBin = min(minBin, minval(binVariable) )
         maxBin = max(maxBin, maxval(binVariable) )

         block => block % next
      end do

      call mpas_dmpar_min_real_array(dminfo, 1, minBin, minBinDomain)
      call mpas_dmpar_max_real_array(dminfo, 1, maxBin, maxBindomain)

      ! set up bins. 
      binWidth = (maxBinDomain(1) - minBinDomain(1)) / nZonalMeanBins
      binBoundaryZonalMean(1) = minBinDomain(1)
      binCenterZonalMean(1) = minBinDomain(1) + binWidth/2.0
      do iBin = 2,nZonalMeanBins
         binBoundaryZonalMean(iBin) = binBoundaryZonalMean(iBin-1) + binWidth
         binCenterZonalMean(iBin) = binCenterZonalMean(iBin-1) + binWidth
      end do
      binBoundaryZonalMean(nZonalMeanBins+1) = binBoundaryZonalMean(nZonalMeanBins) + binWidth

   end subroutine ocn_am_init_zonal_mean!}}}

!***********************************************************************
!
!  routine ocn_am_compute_zonal_mean
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_compute_zonal_mean(domain, timeLevel, err)!{{{

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

      type (amZonalMean_type), pointer :: amZonalMean
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (state_type), pointer :: state
      type (mesh_type), pointer :: mesh
      type (scratch_type), pointer :: scratch
      type (diagnostics_type), pointer :: diagnostics

      integer :: num_tracers, nCellsSolve, nVertLevels
      integer :: iTracer, k, iCell, kMax
      integer ::  nZonalMeanBins, iBin, nZonalMeanVariables, iField
      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:), pointer ::  areaCell, binVariable, binBoundaryZonalMean
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityZonal, normalVelocityMeridional
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:,:,:), allocatable :: sumZonalMean, totalSumZonalMean

      err = 0
      dminfo = domain % dminfo

      call mpas_timer_start("compute_zonal_mean", .false., amZonalMeanTimer)

      amZonalMean => domain % blocklist % amZonalMean
      num_tracers = domain % blocklist % state % time_levs(timeLevel) % state % num_tracers
      nZonalMeanVariables = num_tracers + 3
      nZonalMeanBins = domain % blocklist % mesh % nZonalMeanBins
      nVertLevels = domain % blocklist % mesh % nVertLevels

      binBoundaryZonalMean => amZonalMean % binBoundaryZonalMean % array

      allocate(sumZonalMean(nZonalMeanVariables,nVertLevels,nZonalMeanBins),totalSumZonalMean(nZonalMeanVariables,nVertLevels,nZonalMeanBins))

      sumZonalMean = 0.0

      block => domain % blocklist
      do while (associated(block))
         state => block % state % time_levs(timeLevel) % state
         mesh => block % mesh
         scratch => block % scratch
         diagnostics => block % diagnostics

         num_tracers = state % num_tracers

         nCellsSolve = mesh % nCellsSolve

         areaCell          => mesh % areaCell % array
         maxLevelCell      => mesh % maxLevelCell % array
         tracers           => state % tracers % array
         normalVelocityZonal => diagnostics % normalVelocityZonal % array
         normalVelocityMeridional => diagnostics % normalVelocityMeridional % array

         if (config_ocn_am_zonal_mean_binning_variable=='latCell') then
            binVariable  => mesh % latCell % array
         elseif (config_ocn_am_zonal_mean_binning_variable=='yCell') then
            binVariable  => mesh % yCell % array
         else
            write (stderrUnit,'(a)') ' Error: incorrect choice of config_ocn_am_zonal_mean_binning_variable'
            call mpas_dmpar_abort(dminfo)
         end if 

         ! note that sum is for each vertical index, which is a little wrong for z-star and very wrong for PBCs.
         do iCell = 1,nCellsSolve
            kMax = maxLevelCell(iCell)

            do iBin=1,nZonalMeanBins
               if (binVariable(iCell).lt.binBoundaryZonalMean(iBin+1) ) then

                  do k=1,kMax

                     ! Field 1 is the total area in this bin, which can vary by level due to land.
                     sumZonalMean(1,k,iBin) = sumZonalMean(1,k,iBin) + areaCell(iCell)

                     do iField = 1,num_tracers
                        sumZonalMean(iField+1,k,iBin) = sumZonalMean(iField+1,k,iBin) + tracers(iField,k,iCell)*areaCell(iCell)
                     enddo

                     iField = num_tracers+2
                     sumZonalMean(iField,k,iBin) = sumZonalMean(iField,k,iBin) + normalVelocityZonal(k,iCell)*areaCell(iCell)
                     iField = iField+1
                     sumZonalMean(iField,k,iBin) = sumZonalMean(iField,k,iBin) + normalVelocityMeridional(k,iCell)*areaCell(iCell)

                  end do
                  exit

               endif
            end do

         end do

         block => block % next
      end do

      ! mpi summation over all processors
      call mpas_dmpar_sum_real_array(dminfo, nVertLevels*nZonalMeanBins*nZonalMeanVariables, sumZonalMean, totalSumZonalMean)

      ! normalize by area
      do iField=1,nZonalMeanVariables
         do iBin=1,nZonalMeanBins
            do k=1,nVertLevels
               totalSumZonalMean(:,k,iBin) = totalSumZonalMean(:,k,iBin) / (totalSumZonalMean(1,k,iBin)+1.0e-34)
            end do
         end do
      end do

      ! Even though these variables do not include an index that is decomposed amongst 
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         amZonalMean => block % amZonalMean

         do iBin=1,nZonalMeanBins
            do k=1,nVertLevels

               do iField = 1,num_tracers
                  amZonalMean % tracersZonalMean % array(iField,k,iBin) = totalSumZonalMean(iField+1,k,iBin)
               enddo

               iField = num_tracers+2
               amZonalMean % zonalVelocityZonalMean % array(k,iBin) = totalSumZonalMean(iField,k,iBin)
               iField = iField+1
               amZonalMean % meridionalVelocityZonalMean % array(k,iBin) = totalSumZonalMean(iField,k,iBin)

            end do
         end do

         block => block % next
      end do

      deallocate(sumZonalMean,totalSumZonalMean)

      call mpas_timer_stop("zonal_mean", amZonalMeanTimer)

   end subroutine ocn_am_compute_zonal_mean!}}}

!***********************************************************************
!
!  routine ocn_am_restart_zonal_mean
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_restart_zonal_mean(domain, err)!{{{

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

   end subroutine ocn_am_restart_zonal_mean!}}}

!***********************************************************************
!
!  routine ocn_am_write_zonal_mean
!
!> \brief   Driver for MPAS-Ocean analysis output
!> \author  Mark Petersen
!> \date    March 2014
!> \details 
!>  This routine writes all output for this MPAS-Ocean analysis member.
!>  At this time this is just a stub, and all analysis output is written
!>  to the output file specified by config_output_name.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_write_zonal_mean(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

   end subroutine ocn_am_write_zonal_mean!}}}

!***********************************************************************
!
!  routine ocn_am_finalize_zonal_mean
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_finalize_zonal_mean(domain, err)!{{{

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

   end subroutine ocn_am_finalize_zonal_mean!}}}

end module ocn_am_zonal_mean

! vim: foldmethod=marker
