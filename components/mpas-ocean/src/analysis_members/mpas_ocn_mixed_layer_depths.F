! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_mixed_layer_depths
!
!> \brief MPAS ocean analysis mode member: mixed_layer_depths
!> \author Luke Van Roekel
!> \date   August 2015
!> \details
!>  MPAS ocean analysis mode member: mixed_layer_depths
!>
!     Computes mixed layer depths via a gradient method and threshold method
!     may add more methods from Holte and Talley (2009) at a future time
!-----------------------------------------------------------------------

module ocn_mixed_layer_depths

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_eddy_parameterization_helpers

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

   public :: ocn_init_mixed_layer_depths, &
             ocn_compute_mixed_layer_depths, &
             ocn_restart_mixed_layer_depths, &
             ocn_finalize_mixed_layer_depths

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_mixed_layer_depths
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    August 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_mixed_layer_depths(domain, err)!{{{

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

   end subroutine ocn_init_mixed_layer_depths!}}}

!***********************************************************************
!
!  routine ocn_compute_mixed_layer_depths
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    August 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_mixed_layer_depths(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: mixedLayerDepthsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: mixedLayerDepthsAM
      type (mpas_pool_type), pointer :: forcingPool

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels
      integer :: k, iCell, i, refIndex, nCells, refLevel(1)
      integer, pointer :: index_temperature
      integer, dimension(:), pointer :: maxLevelCell, nCellsArray

      logical :: found_temp_mld, found_den_mld
!      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:,:), allocatable :: pressureAdjustedForLandIce 

      real (kind=RKIND), dimension(:), pointer :: tThreshMLD, tGradientMLD, landIceDraft
      real (kind=RKIND), dimension(:), pointer :: dGradientMLD, landIcePressure
      integer, dimension(:), pointer :: landIceMask
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer :: interp_local
      real (kind=RKIND), allocatable, dimension(:,:) :: densityGradient, temperatureGradient
      real (kind=RKIND) :: mldTemp,dTempThres, dTempGrad, dDenGrad
      real (kind=RKIND) :: dz,temp_ref_lev, dV, dVm1, dVp1, localVals(6)
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mixedLayerDepthsAM', mixedLayerDepthsAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

      if (config_AM_mixedLayerDepths_interp_method == 1) interp_local = 1
      if (config_AM_mixedLayerDepths_interp_method == 2) interp_local = 2
      if (config_AM_mixedLayerDepths_interp_method == 3) interp_local = 3

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mixedLayerDepthsAM', mixedLayerDepthsAMPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         !Use work array for storing of pressure, potentially modified by landIce
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
         call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

         !need to call the mixed layer depth computation for compute on startup 
         call ocn_eddy_compute_mixed_layer_depth(statePool, forcingPool)

         nCells = nCellsArray( size(nCellsArray) )
         allocate(pressureAdjustedForLandIce(nVertLevels, size(pressure,2)))

         pressureAdjustedForLandIce(:,:) = pressure(:,:)
         !If landice cavity remove land ice pressure to search for ML depth
         if ( associated(landIcePressure) ) then
           !$omp parallel
           !$omp do schedule(runtime) private(k)
           do iCell = 1,nCells
             do k = 1,maxLevelCell(iCell)
                pressureAdjustedForLandIce(k,iCell) = pressureAdjustedForLandIce(k,iCell)   &
                                                        - landIcePressure(iCell)
             end do
           end do
           !$omp end do
           !$omp end parallel
         end if

         if(config_AM_mixedLayerDepths_Tthreshold) then
            call mpas_pool_get_array(mixedLayerDepthsAMPool, 'tThreshMLD',tThreshMLD)
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(refIndex, found_temp_mld, k, localvals, temp_ref_lev, dVp1, dV, mldTemp)
            do iCell = 1,nCells

               !Initialize RefIndex for cases of very shallow columns
               refIndex = maxLevelCell(iCell)

               found_temp_mld = .false.

               do k=1, maxLevelCell(iCell)-1
                  if(pressureAdjustedForLandIce(k+1,iCell) > config_AM_mixedLayerDepths_reference_pressure) then
                     localvals(2:3)=tracers(index_temperature,k:k+1,iCell)
                     localvals(5:6)=pressureAdjustedForLandIce(k:k+1,iCell)

                     call interp_bw_levels(localVals(2),localVals(3),localVals(5),localVals(6), &
                                   config_AM_mixedLayerDepths_reference_pressure,interp_local,temp_ref_lev)
                     refIndex=k
                     exit
                   end if
               end do

               do k=refIndex,maxLevelCell(iCell)-1

                    if(.not. found_temp_mld .and. abs(tracers(index_temperature,k+1,iCell) - temp_ref_lev) &
                        .ge. config_AM_mixedLayerDepths_crit_temp_threshold) then
                        dVp1 = abs(tracers(index_temperature,k+1,iCell) - temp_ref_lev)
                        dV   = abs(tracers(index_temperature,k  ,iCell) - temp_ref_lev)

                        localVals(2:3)=zMid(k:k+1,iCell)
                        call interp_bw_levels(localVals(2),localVals(3), dV, dVp1, config_AM_mixedLayerDepths_crit_temp_threshold,     &
                                   interp_local, mldTemp)!,dVm1, localVals(1))
                        mldTemp=max(mldTemp,zMid(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                        tThreshMLD(iCell)=abs(min(mldTemp,zMid(k,iCell))) !MLD should be deeper than zMid(k)
                        found_temp_mld = .true.
                        exit
                     end if
                end do

! if the mixed layer depth is not found, it is set to the depth of the bottom most level
                   if(.not. found_temp_mld) tThreshMLD(iCell) = abs(zMid(maxLevelCell(iCell),iCell))
               end do !iCell
               !$omp end do
               !$omp end parallel

               if (associated(landIceMask) ) then
                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iCell = 1, nCells
                     tThreshMLD(iCell) = tThreshMLD(iCell) - abs(landIceDraft(iCell))*landIceMask(iCell)
                  end do
                  !$omp end do
                  !$omp end parallel
               end if

          end if !end tThresh MLD search

! Compute the mixed layer depth based on a gradient threshold in temperature and density
        if(config_AM_mixedLayerDepths_Tgradient) then
           call mpas_pool_get_array(mixedLayerDepthsAMPool, 'tGradMLD',tGradientMLD)

           allocate(temperatureGradient(nVertLevels,2))

           !$omp parallel
           !$omp do schedule(runtime) &
           !$omp private(temperatureGradient, found_temp_mld, k, dz, mldTemp, refLevel)
           do iCell = 1,nCells

            temperatureGradient(:,1) = 0.0_RKIND
            temperatureGradient(1,2) = 1

               found_temp_mld=.false.

               do k=2,maxLevelCell(iCell)-1
                     dz=abs(pressureAdjustedForLandIce(k-1,iCell)-pressureAdjustedForLandIce(k,iCell)) + 1.0E-15_RKIND
                     temperatureGradient(k,1) = abs(tracers(index_temperature,k-1,iCell) - tracers(index_temperature,k,iCell)) / dz
                     temperatureGradient(k,2) = k
               end do

! smooth the gradients to eliminate reduce single point maxima

               do k=2,maxLevelCell(iCell)-1
                  temperatureGradient(k,1) = (temperatureGradient(k-1,1) + temperatureGradient(k,1) &
                                           + temperatureGradient(k+1,1)) / float(3)
               end do


               do k=2, maxLevelCell(iCell)-1
                  if(.not. found_temp_mld .and. temperatureGradient(k+1,1) .ge. config_AM_mixedLayerDepths_temp_gradient_threshold) then
                    call interp_bw_levels(zTop(k,iCell),zTop(k+1,iCell),temperatureGradient(k,1),temperatureGradient(k+1,1), &
                           config_AM_mixedLayerDepths_temp_gradient_threshold, interp_local,mldTemp,temperatureGradient(k-1,1),zTop(k-1,iCell))

                        mldTemp=max(mldTemp,zTop(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                        tGradientMLD(iCell)=abs(min(mldTemp,zTop(k,iCell))) !MLD should be deeper than zMid(k)

                    found_temp_mld=.true.
                    exit
                  end if

              end do !maxLevelCell

                if(.not. found_temp_mld) then
                    refLevel=maxloc(temperatureGradient(:,1))
                    tGradientMLD(iCell) = abs(zTop(refLevel(1),iCell))
                end if

         end do !iCell
         !$omp end do
         !$omp end parallel

         !normalize MLD to top of ice cavity
         if (associated(landIceMask) ) then
           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
             tGradientMLD(iCell) = tGradientMLD(iCell) - abs(landIceDraft(iCell))*landIceMask(iCell)
           end do
           !$omp end do
           !$omp end parallel
         end if


         deallocate(temperatureGradient)

       end if !if(temperaturegradientflag)

           if(config_AM_mixedLayerDepths_Dgradient) then
                 call mpas_pool_get_array(mixedLayerDepthsAMPool, 'dGradMLD',dGradientMLD)

             allocate(densityGradient(nVertLevels,2))

           !$omp parallel
           !$omp do schedule(runtime) &
           !$omp private(densityGradient, found_den_mld, k, dz, mldTemp, refLevel)
           do iCell = 1,nCells

            densityGradient(:,:)=0.0_RKIND
            densityGradient(1,2) = 1

               found_den_mld=.false.

               do k=2,maxLevelCell(iCell)-1
                     dz=abs(pressureAdjustedForLandIce(k-1,iCell)-pressureAdjustedForLandIce(k,iCell)) + 1.0E-15_RKIND
                     densityGradient(k,1) = abs(potentialDensity(k-1,iCell)-potentialDensity(k,iCell)) / dz
                     densityGradient(k,2) = k
               end do

! smooth the gradients to eliminate reduce single point maxima

               do k=2,maxLevelCell(iCell)-1
                  densityGradient(k,1) = (densityGradient(k-1,1) + densityGradient(k,1) + densityGradient(k+1,1)) / float(3)
               end do


               do k=2, maxLevelCell(iCell)-1
                  if(.not. found_den_mld .and. densityGradient(k+1,1) .ge. config_AM_mixedLayerDepths_den_gradient_threshold) then
                    call interp_bw_levels(zTop(k,iCell),zTop(k+1,iCell),densityGradient(k,1),densityGradient(k+1,1),           &
                           config_AM_mixedLayerDepths_den_gradient_threshold, interp_local,mldTemp,densityGradient(k-1,1),zTop(k-1,iCell))
                        mldTemp=max(mldTemp,zTop(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                        dGradientMLD(iCell)=abs(min(mldTemp,zTop(k,iCell))) !MLD should be deeper than zMid(k)
                    found_den_mld=.true.
                    exit
                  end if

              end do !maxLevelCell


                if(.not. found_den_mld) then
                    refLevel=maxloc(densityGradient(:,2))
                    dGradientMLD(iCell) = abs(zTop(refLevel(1),iCell))
                end if

         end do !iCell
         !$omp end do
         !$omp end parallel

         if (associated(landIceMask) ) then
           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
             dGradientMLD(iCell) = dGradientMLD(iCell) - abs(landIceDraft(iCell))*landIceMask(iCell)
           end do
           !$omp end do
           !$omp end parallel
         end if

         deallocate(densityGradient)
       end if !if(densitygradientflag)

         deallocate(pressureAdjustedForLandIce)
         block => block % next
      end do

      ! Even though some variables do not include an index that is decomposed amongst
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
!      block => domain % blocklist
!      do while (associated(block))
!         call mpas_pool_get_subpool(block % structs, 'temPlateAM', temPlateAMPool)
!
!         ! assignment of final temPlateAM variables could occur here.
!
!         block => block % next
!      end do

   end subroutine ocn_compute_mixed_layer_depths!}}}

!***********************************************************************
!
!  routine interp_bw_levels
!
!> \brief   Interpolates between model layers
!> \author  Luke Van Roekel
!> \date    September 2015
!> \details
!>  This routine conducts computations to compute various field values
!>     between model levels (in pressure or depth) or could interpolate
!>      between temperature/salinity/density values.  Interpolations are
!>      of the form
!>       y = coeffs(1)*x^3 + coeffs(2)*x^2 + coeffs(3)*x + coeffs(4)
!
!-----------------------------------------------------------------------

   subroutine interp_bw_levels(y0,y1,x0,x1,xT,interp_f,yT,xm1,ym1)!{{{

   integer,intent(in) :: interp_f  ! linear, quadratic, or spline
   real(kind=RKIND),intent(in)  :: y0,y1,x0,x1,xT
   real(kind=RKIND),intent(inout) :: yT
   real(kind=RKIND),optional,intent(in) :: xm1,ym1
                ! these values are to match the slope at a given point

!------------------------------------------------------------------------
!
!  Local variables for the interpolations
!
!------------------------------------------------------------------------

   real(kind=RKIND) :: coeffs(4)   ! stores the coefficients for the interp
   real(kind=RKIND) :: Minv(4,4)   ! holds values for computing quad and spline
   real(kind=RKIND) :: det
   real(kind=RKIND) :: rhs(4)
   integer :: k,k2

   coeffs(:) = 0.0_RKIND
   Minv(:,:) = 0.0_RKIND
   rhs(:)    = 0.0_RKIND


   select case (interp_f)

       case (1)

          coeffs(2) = (y1-y0)/(x1-x0)
          coeffs(1) = y0 - coeffs(2)*x0
       case (2)

          det = -(x1-x0)**2
          rhs(1) = y0
          rhs(2) = y1

          if(present(xm1) .and. present(ym1)) then
             rhs(3) = (y0-ym1)/(x0-xm1)
          else
             rhs(3) = 0.0_RKIND
          end if

          Minv(1,1) = -1.0_RKIND/det
          Minv(1,2) = 1.0_RKIND/det
          Minv(1,3) = -1.0_RKIND/(x1-x0)
          Minv(2,1) = 2.0_RKIND*x0/det
          Minv(2,2) = -2.0_RKIND*x0/det
          Minv(2,3) = (x1+x0)/(x1-x0)
          Minv(3,1) = -(x0**2)/det
          Minv(3,2) = x1*(2.0_RKIND*x0-x1)/det
          Minv(3,3) = -x1*x0/(x1-x0)

          do k=1,3
             do k2=1,3
                coeffs(k2) = coeffs(k2) + Minv(4-k2,k)*rhs(k)
             end do
          end do

        case (3)
          det = -(x1-x0)**3
          rhs(1) = y1
          rhs(2) = y0
          if(present(xm1) .and. present(ym1)) then
             rhs(3) = (y0-ym1)/(x0-xm1)
          else
             rhs(3) = 0.0_RKIND
          end if

          rhs(4) = (y1-y0)/(x1-x0)

          Minv(1,1) = 2.0_RKIND/det
          Minv(1,2) = -2.0_RKIND/det
          Minv(1,3) = (x0-x1)/det
          Minv(1,4) = (x0-x1)/det
          Minv(2,1) = -3.0_RKIND * (x1+x0)/det
          Minv(2,2) = 3.0_RKIND*(x1+x0)/det
          Minv(2,3) = (x1-x0)*(2.0_RKIND*x1+x0)/det
          Minv(2,4) = (x1-x0)*(2.0_RKIND*x0+x1)/det
          Minv(3,1) = 6.0_RKIND*x1*x0/det
          Minv(3,2) = -6.0_RKIND*x1*x0/det
          Minv(3,3) = -x1*(x1-x0)*(2.0_RKIND*x0+x1)/det
          Minv(3,4) = -x0*(x1-x0)*(2.0_RKIND*x1+x0)/det
          Minv(4,1) = -(x0**2)*(3.0_RKIND*x1-x0)/det
          Minv(4,2) = -(x1**2)*(-3.0_RKIND*x0+x1)/det
          Minv(4,3) = x0*(x1**2)*(x1-x0)/det
          Minv(4,4) = x1*(x0**2)*(x1-x0)/det

          do k=1,4
             do k2=1,4
                coeffs(k2) = coeffs(k2) + Minv(5-k2,k)*rhs(k)
             end do
          end do

     end select

     yT = coeffs(4)*xT**3 + coeffs(3)*xT**2 + coeffs(2)*xT + coeffs(1)
   end subroutine interp_bw_levels!}}}

!***********************************************************************
!
!  routine ocn_restart_mixed_layer_depths
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    September 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_mixed_layer_depths(domain, err)!{{{

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

   end subroutine ocn_restart_mixed_layer_depths!}}}

!***********************************************************************
!
!  routine ocn_finalize_mixed_layer_depths
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    August 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_mixed_layer_depths(domain, err)!{{{

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

   end subroutine ocn_finalize_mixed_layer_depths!}}}

end module ocn_mixed_layer_depths

! vim: foldmethod=marker
