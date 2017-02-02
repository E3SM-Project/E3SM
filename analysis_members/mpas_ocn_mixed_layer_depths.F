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
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: mixedLayerDepthsAM

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve
      integer :: k, iCell, i, refIndex, refLevel(1)
      integer, pointer :: index_temperature
      integer, dimension(:), pointer :: maxLevelCell

      logical :: found_temp_mld, found_den_mld
      logical,pointer :: tThresholdFlag, dThresholdFlag
      logical,pointer :: tGradientFlag, dGradientFlag
!      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:), pointer :: tThresholdMLD, tGradientMLD
      real (kind=RKIND), dimension(:), pointer :: dThresholdMLD, dGradientMLD
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:,:), pointer :: zTop, zMid, pressure
      real (kind=RKIND), dimension(:,:), pointer :: potentialDensity
      real (kind=RKIND), pointer :: tempThresh
      real (kind=RKIND), pointer :: tempGrad
      real (kind=RKIND), pointer :: denThresh
      real (kind=RKIND), pointer :: denGrad
      integer, pointer :: interp_type
      integer :: interp_local
      real (kind=RKIND), pointer :: refPress
      real (kind=RKIND), allocatable, dimension(:,:) :: densityGradient, temperatureGradient
      real (kind=RKIND) :: mldTemp,dTempThres, dDenThres, dTempGrad, dDenGrad
      real (kind=RKIND) :: dz,temp_ref_lev, den_ref_lev, dV, dVm1, dVp1, localVals(6)
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mixedLayerDepthsAM', mixedLayerDepthsAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_Tthreshold', tThresholdFlag)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_Dthreshold', dThresholdFlag)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_Tgradient', tGradientFlag)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_Dgradient', dGradientFlag)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_crit_temp_threshold', tempThresh)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_crit_dens_threshold', denThresh)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_temp_gradient_threshold', tempGrad)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_den_gradient_threshold', denGrad)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_interp_method', interp_type)
      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_reference_pressure', refPress)

      if (interp_type == 1) interp_local = 1
      if (interp_type == 2) interp_local = 2
      if (interp_type == 3) interp_local = 3

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'mixedLayerDepthsAM', mixedLayerDepthsAMPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', &
                        potentialDensity)
         call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)


         if(tThresholdFlag) then
            call mpas_pool_get_array(mixedLayerDepthsAMPool, 'tThreshMLD',tThresholdMLD)
            do iCell = 1,nCellsSolve

               found_temp_mld = .false.

               do k=1, maxLevelCell(iCell)-1
                  if(pressure(k+1,iCell) > refPress) then
                     localvals(2:3)=tracers(index_temperature,k:k+1,iCell)
                     localvals(5:6)=pressure(k:k+1,iCell)

                     call interp_bw_levels(localVals(2),localVals(3), &
                                   localVals(5),localVals(6),refPress,interp_local,                &
                                   temp_ref_lev)
                     refIndex=k
                     exit
                   endif
               enddo

               do k=refIndex,maxLevelCell(iCell)-1

                    if(.not. found_temp_mld .and. abs(tracers(index_temperature,k+1,iCell) - temp_ref_lev) .ge. tempThresh) then
                        dVp1 = abs(tracers(index_temperature,k+1,iCell) - temp_ref_lev)
                        dV   = abs(tracers(index_temperature,k  ,iCell) - temp_ref_lev)

                        localVals(2:3)=zMid(k:k+1,iCell)
                        call interp_bw_levels(localVals(2),localVals(3), dV, dVp1, tempThresh,     &
                                   interp_local, mldTemp)!,dVm1, localVals(1))
                        mldTemp=max(mldTemp,zMid(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                        tThresholdMLD(iCell)=abs(min(mldTemp,zMid(k,iCell))) !MLD should be deeper than zMid(k)
                        found_temp_mld = .true.
                        exit
                     endif
                enddo

! if the mixed layer depth is not found, it is set to the depth of the bottom most level
                   if(.not. found_temp_mld) tThresholdMLD(iCell) = abs(zMid(maxLevelCell(iCell),iCell))
               enddo !iCell
          endif !end tThresh MLD search

         if(dThresholdFlag) then
            call mpas_pool_get_array(mixedLayerDepthsAMPool, 'dThreshMLD',dThresholdMLD)

            do iCell = 1,nCellsSolve

               found_den_mld = .false.

               do k=1, maxLevelCell(iCell)-1
                  if(pressure(k+1,iCell) > refPress) then
                     localvals(2:3)=potentialDensity(k:k+1,iCell)
                     localvals(5:6)=pressure(k:k+1,iCell)

                     call interp_bw_levels(localVals(2),localVals(3), &
                                   localVals(5),localVals(6),refPress,interp_local,                &
                                   den_ref_lev)
                     refIndex=k
                    exit
                   endif
               enddo

               do k=refIndex,maxLevelCell(iCell)-1

                    if(.not. found_den_mld .and. abs(potentialDensity(k+1,iCell) - den_ref_lev) .ge. denThresh) then
                        dVp1 = abs(potentialDensity(k+1,iCell) - den_ref_lev)
                        dV   = abs(potentialDensity(k  ,iCell) - den_ref_lev)
                        localVals(2:3)=zMid(k:k+1,iCell)
                        call interp_bw_levels(localVals(2),localVals(3), dV, dVp1, denThresh,     &
                                   interp_local, mldTemp)!, dVm1, localVals(1))
                        mldTemp=max(mldTemp,zMid(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                        dThresholdMLD(iCell)=abs(min(mldTemp,zMid(k,iCell))) !MLD should be deeper than zMid(k)
                        found_den_mld = .true.
                        exit
                     endif
               enddo

! if the mixed layer depth is not found, it is set to the depth of the bottom most level

                    if(.not. found_den_mld) dThresholdMLD(iCell) = abs(zMid(maxLevelCell(iCell),iCell))

             enddo !iCell

        endif !end dThresh MLD search


! Compute the mixed layer depth based on a gradient threshold in temperature and density
        if(tGradientFlag) then
           call mpas_pool_get_array(mixedLayerDepthsAMPool, 'tGradMLD',tGradientMLD)

           allocate(temperatureGradient(nVertLevels,2))

           do iCell = 1,nCellsSolve

            temperatureGradient(:,1) = 0.0_RKIND
            temperatureGradient(1,2) = 1

               found_temp_mld=.false.

               do k=2,maxLevelCell(iCell)-1
                     dz=abs(pressure(k-1,iCell)-pressure(k,iCell))
                     temperatureGradient(k,1) = abs(tracers(index_temperature,k-1,iCell) - tracers(index_temperature,k,iCell)) / dz
                     temperatureGradient(k,2) = k
               enddo

! smooth the gradients to eliminate reduce single point maxima

               do k=2,maxLevelCell(iCell)-1
                  temperatureGradient(k,1) = (temperatureGradient(k-1,1) + temperatureGradient(k,1) &
                                           + temperatureGradient(k+1,1)) / float(3)
               enddo


               do k=2, maxLevelCell(iCell)-1
                  if(.not. found_temp_mld .and. temperatureGradient(k+1,1) .ge. tempGrad) then
                    call interp_bw_levels(zTop(k,iCell),zTop(k+1,iCell),temperatureGradient(k,1),temperatureGradient(k+1,1), &
                           tempGrad, interp_local,mldTemp,temperatureGradient(k-1,1),zTop(k-1,iCell))

                        mldTemp=max(mldTemp,zTop(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                        tGradientMLD(iCell)=abs(min(mldTemp,zTop(k,iCell))) !MLD should be deeper than zMid(k)

                    found_temp_mld=.true.
                    exit
                  endif

              enddo !maxLevelCell

                if(.not. found_temp_mld) then
                    refLevel=maxloc(temperatureGradient(:,1))
                    tGradientMLD(iCell) = abs(zTop(refLevel(1),iCell))
                endif

         enddo !icell

         deallocate(temperatureGradient)

       endif !if(temperaturegradientflag)

           if(dGradientFlag) then
                 call mpas_pool_get_array(mixedLayerDepthsAMPool, 'dGradMLD',dGradientMLD)

             allocate(densityGradient(nVertLevels,2))

           do iCell = 1,nCellsSolve

            densityGradient(:,1)=0.0_RKIND
            densityGradient(1,2) = 1

               found_den_mld=.false.

               do k=2,maxLevelCell(iCell)-1
                     dz=abs(pressure(k-1,iCell)-pressure(k,iCell))
                     densityGradient(k,1) = abs(potentialDensity(k-1,iCell)-potentialDensity(k,iCell)) / dz
                     densityGradient(k,2) = k
               enddo

! smooth the gradients to eliminate reduce single point maxima

               do k=2,maxLevelCell(iCell)-1
                  densityGradient(k,1) = (densityGradient(k-1,1) + densityGradient(k,1) + densityGradient(k+1,1)) / float(3)
               enddo


               do k=2, maxLevelCell(iCell)-1
                  if(.not. found_den_mld .and. densityGradient(k+1,1) .ge. denGrad) then
                    call interp_bw_levels(zTop(k,iCell),zTop(k+1,iCell),densityGradient(k,1),densityGradient(k+1,1),           &
                           denGrad, interp_local,mldTemp,densityGradient(k-1,1),zTop(k-1,iCell))
                        mldTemp=max(mldTemp,zTop(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                        dGradientMLD(iCell)=abs(min(mldTemp,zTop(k,iCell))) !MLD should be deeper than zMid(k)
                    found_den_mld=.true.
                    exit
                  endif

              enddo !maxLevelCell


                if(.not. found_den_mld) then
                    refLevel=maxloc(densityGradient(:,2))
                    dGradientMLD(iCell) = abs(zTop(refLevel(1),iCell))
                endif

         enddo !icell

         deallocate(densityGradient)
       endif !if(densitygradientflag)

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
          endif

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
             enddo
          enddo

        case (3)
          det = -(x1-x0)**3
          rhs(1) = y1
          rhs(2) = y0
          if(present(xm1) .and. present(ym1)) then
             rhs(3) = (y0-ym1)/(x0-xm1)
          else
             rhs(3) = 0.0_RKIND
          endif

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
             enddo
          enddo

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
