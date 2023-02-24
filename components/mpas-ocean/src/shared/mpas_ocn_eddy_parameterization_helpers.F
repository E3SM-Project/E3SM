! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!

module ocn_eddy_parameterization_helpers

   use mpas_pool_routines
   use mpas_derived_types
   use mpas_constants
   use mpas_threading
   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_mesh

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

   public :: ocn_eddy_compute_buoyancy_gradient, &
             ocn_eddy_compute_mixed_layer_depth


  !***********************************************************************

  contains

  !***********************************************************************
  !
  !  routine ocn_eddy_compute_buoyancy_gradient
  !
  !> \brief Computes fixed depth horizontal buoyancy gradient 
  !> \details
  !>  Computes the horizontal buoyancy gradient on fixed depth levels. 
  !>  This quantity is extracted from mpas_ocn_gm as it is needed for both
  !>  the GM eddy parameterization and the submesoscale eddy parameterization
  !***********************************************************************

  subroutine ocn_eddy_compute_buoyancy_gradient()

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     integer :: k, nEdges, nCells, iCell, iEdge, cell1, cell2
     ! gradDensityEdge: Normal gradient of density
     !           units: none
     real(kind=RKIND), dimension(:, :), allocatable :: gradDensityEdge, gradDensityTopOfEdge, &
        dDensityDzTopOfCell, dDensityDzTopOfEdge, &
        gradZMidEdge, gradZMidTopOfEdge
     real(kind=RKIND) :: rtmp, h1, h2
     real(kind=RKIND), parameter :: epsGM = 1.0E-12_RKIND

     nCells = nCellsAll
     nEdges = nEdgesAll
     allocate(gradDensityEdge(nVertLevels, nEdges), &
              dDensityDzTopOfCell(nVertLevels, nCells), &
              gradDensityTopOfEdge(nVertLevels, nEdges), &
              dDensityDzTopOfEdge(nVertLevels, nEdges), &
              gradZMidEdge(nVertLevels, nEdges), &
              gradZMidTopOfEdge(nVertLevels, nEdges))

     !$omp parallel
     !$omp do schedule(runtime) private(k)
     do iEdge=1,nEdgesAll
        do k=1,nVertLevels
           gradDensityEdge(k,iEdge) = 0.0_RKIND
        end do
     end do
     !$omp end do
     !$omp end parallel

     nEdges = nEdgesHalo(2)
     !$omp parallel
     !$omp do schedule(runtime) private(k, rtmp)
     do iCell = 1, nCells
        do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
           rtmp = (displacedDensity(k-1,iCell) - density(k,iCell)) / (zMid(k-1,iCell) - zMid(k,iCell))
           dDensityDzTopOfCell(k,iCell) = min(rtmp, -epsGM)
        end do
        ! Approximation of dDensityDzTopOfCell on the top and bottom interfaces through the idea of having
        ! ghost cells above the top and below the bottom layers of the same depths and density.
        ! Essentially, this enforces the boundary condition (d density)/dz = 0 at the top and bottom.
        dDensityDzTopOfCell(1:minLevelCell(iCell),iCell) = 0.0_RKIND
     end do
     !$omp end do
     !$omp end parallel

     nEdges = nEdgesHalo( 2 )

     ! Interpolate dDensityDzTopOfCell to edge and layer interface
     !$omp parallel
     !$omp do schedule(runtime) private(k, cell1, cell2)
     do iEdge = 1, nEdges
        do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           dDensityDzTopOfEdge(k,iEdge) = 0.5_RKIND * (dDensityDzTopOfCell(k,cell1) + dDensityDzTopOfCell(k,cell2))
        end do
     end do
     !$omp end do
     !$omp end parallel

     ! Compute density gradient (gradDensityEdge) and gradient of zMid (gradZMidEdge)
     ! along the constant coordinate surface.
     ! The computed variables lives at edge and mid-layer depth
     !$omp parallel
     !$omp do schedule(runtime) private(cell1, cell2, k)
     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        do k=minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
           gradDensityEdge(k,iEdge) = (density(k,cell2) - density(k,cell1)) / dcEdge(iEdge)
           gradZMidEdge(k,iEdge) = (zMid(k,cell2) - zMid(k,cell1)) / dcEdge(iEdge)
        end do
     end do
     !$omp end do

     !$omp do schedule(runtime) private(k, h1, h2)
     do iEdge = 1, nEdges
        ! The interpolation can only be carried out on non-boundary edges
        if (maxLevelEdgeTop(iEdge) .GE. minLevelEdgeBot(iEdge)) then
           do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
              h1 = layerThickEdgeMean(k-1,iEdge)
              h2 = layerThickEdgeMean(k,iEdge)
              ! Using second-order interpolation below
              gradDensityTopOfEdge(k,iEdge) = (h2 * gradDensityEdge(k-1,iEdge) + h1 * &
                                               gradDensityEdge(k,iEdge)) / (h1 + h2)
              gradZMidTopOfEdge(k,iEdge) = (h2 * gradZMidEdge(k-1,iEdge) + h1 * &
                                            gradZMidEdge(k,iEdge)) / (h1 + h2)
           end do

           ! Approximation of values on the top and bottom interfaces through the idea of having ghost cells
           ! above the top and below the bottom layers of the same depths and density.
           gradDensityTopOfEdge(minLevelEdgeBot(iEdge),iEdge) = gradDensityEdge(minLevelEdgeBot(iEdge),iEdge)
           gradZMidTopOfEdge(minLevelEdgeBot(iEdge),iEdge) = gradZMidEdge(minLevelEdgeBot(iEdge),iEdge)
        end if
     end do
     !$omp end do

     !$omp do schedule(runtime) private(k)
     do iEdge = 1, nEdges
        if (maxLevelEdgeTop(iEdge) .GE. minLevelEdgeBot(iEdge)) then
           do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
              gradDensityEddy(k,iEdge) = (gradDensityTopOfEdge(k,iEdge) - dDensityDzTopOfEdge(k,iEdge) &
                                             * gradZMidTopOfEdge(k,iEdge))
           end do
        end if
     end do
     !$omp end do
     !$omp end parallel
     deallocate(gradDensityEdge, &
                dDensityDzTopOfCell, &
                gradDensityTopOfEdge, &
                dDensityDzTopOfEdge, &
                gradZMidEdge, &
                gradZMidTopOfEdge)

  end subroutine ocn_eddy_compute_buoyancy_gradient

  !***********************************************************************
  !
  !  routine ocn_eddy_compute_mixed_layer_depth
  !
  !> \brief   Computes the density threshold based mixed layer depth
  !> \details
  !>   Computes the density threshold based mixed layer depth.  This 
  !>   quantity is required by the submesoscale and redi parameterizations
  !>   The calculation is moved out of analysis members for this reason
  !***********************************************************************

  subroutine ocn_eddy_compute_mixed_layer_depth(statePool, forcingPool)

     type (mpas_pool_type), pointer, intent(in) :: statePool, forcingPool

     integer :: iEdge, nEdges, refIndex, cell1, cell2, nCells, k, iCell
     real(kind=RKIND) :: dDenThres, den_ref_lev
     real(kind=RKIND),dimension(:), allocatable :: depth
     integer, dimension(:), pointer :: landIceMask
     real (kind=RKIND), dimension(:,:), pointer :: layerThickness
     real (kind=RKIND), dimension(:), pointer :: landIceDraft, landIcePressure
     logical :: found_den_mld
     real (kind=RKIND) :: dV, dVp1, refDepth, coeffs(2), mldTemp, dz
     real (kind=RKIND), dimension(:,:), allocatable :: pressureAdjustedForLandIce

     nCells = nCellsAll

     if ( config_eddyMLD_use_old ) then
       call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
       call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
       call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

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

       !$omp parallel
       !$omp do schedule(runtime) &
       !$omp private(refIndex, found_den_mld, k, coeffs, den_ref_lev, dVp1, dV, mldTemp)
       do iCell = 1,nCells

          !Initialize RefIndex for cases of very shallow columns
          refIndex = maxLevelCell(iCell)
          found_den_mld = .false.

          do k=1, maxLevelCell(iCell)-1
             if(pressureAdjustedForLandIce(k+1,iCell) > config_eddyMLD_reference_pressure) then

                coeffs(2) = (potentialDensity(k+1,iCell) - potentialDensity(k,iCell)) /  &
                            (pressureAdjustedForLandIce(k+1,iCell) - pressureAdjustedForLandIce(k,iCell))
                coeffs(1) = potentialDensity(k,iCell) - coeffs(2)*pressureAdjustedForLandIce(k,iCell)

                den_ref_lev = coeffs(2)*config_eddyMLD_reference_pressure + coeffs(1)

                refIndex=k
                exit
             end if
          end do

          do k=refIndex,maxLevelCell(iCell)-1

            if(.not. found_den_mld .and. abs(potentialDensity(k+1,iCell) - den_ref_lev) .ge. &
                     config_eddyMLD_dens_threshold) then
                dVp1 = abs(potentialDensity(k+1,iCell) - den_ref_lev)
                dV   = abs(potentialDensity(k  ,iCell) - den_ref_lev)

                coeffs(2) = (zMid(k+1,iCell) - zMid(k,iCell)) / (dVp1 - dV)
                coeffs(1) = zMid(k,iCell) - coeffs(2)*dV

                mldTemp = coeffs(2)*config_eddyMLD_dens_threshold + coeffs(1)

                mldTemp=max(mldTemp,zMid(k+1,iCell)) !make sure MLD isn't deeper than zMid(k+1)
                dThreshMLD(iCell)=abs(min(mldTemp,zMid(k,iCell))) !MLD should be deeper than zMid(k)
                indMLD(iCell) = k+1
                found_den_mld = .true.
                exit
             end if
          end do

          ! if the mixed layer depth is not found, it is set to the depth of the bottom most level

          if(.not. found_den_mld) then
             dThreshMLD(iCell) = abs(zMid(maxLevelCell(iCell),iCell))
             indMLD(iCell) = maxLevelCell(iCell)
          end if
        end do !iCell
        !$omp end do
        !$omp end parallel

        if (associated(landIceMask) ) then
          !$omp parallel
          !$omp do schedule(runtime)
          do iCell = 1, nCells
             dThreshMLD(iCell) = dThreshMLD(iCell) - abs(landIceDraft(iCell))*landIceMask(iCell)
          end do
          !$omp end do
          !$omp end parallel
        end if

        deallocate(pressureAdjustedForLandIce)

    else ! fixed mixed layer depth code
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
       allocate(depth(nVertLevels))
       refDepth = config_eddyMLD_reference_depth

       !$omp parallel
       !$omp do schedule(runtime) &
       !$omp private(depth, refIndex, coeffs, found_den_mld, k, den_ref_lev, dVp1, dV, mldTemp)
       do iCell = 1,nCells

          found_den_mld = .false.

          depth(1) = 0.5_RKIND*layerThickness(1,iCell)
          do k=2,maxLevelCell(iCell)
             depth(k) = depth(k-1) + 0.5_RKIND*(layerThickness(k-1,iCell) + &
                                                layerThickness(k  ,iCell))
          end do

          do k=1, maxLevelCell(iCell)-1
             if(depth(k+1) >= refDepth) then
                coeffs(2) = (potentialDensity(k+1,iCell) - potentialDensity(k,iCell)) &
                             / (depth(k+1) - depth(k))
                coeffs(1) = potentialDensity(k,iCell) - coeffs(2)*depth(k)

                den_ref_lev = coeffs(2)*refDepth + coeffs(1)

                refIndex=k
                exit
             end if
          end do

          do k=refIndex,maxLevelCell(iCell)-1

             if(.not. found_den_mld .and. (potentialDensity(k+1,iCell) - den_ref_lev) .ge. &
                      config_eddyMLD_dens_threshold) then
                dVp1 = (potentialDensity(k+1,iCell) - den_ref_lev)
                dV   = (potentialDensity(k  ,iCell) - den_ref_lev)

                coeffs(2) = (depth(k+1) - depth(k)) / (dVp1 - dV)
                coeffs(1) = depth(k) - coeffs(2)*dV

                mldTemp = coeffs(2)*config_eddyMLD_dens_threshold + coeffs(1)
                mldTemp=min(mldTemp,depth(k+1)) !make sure MLD isn't deeper than zMid(k+1)
                dThreshMLD(iCell)=abs(max(mldTemp,depth(k))) !MLD should be deeper than zMid(k)
                indMLD(iCell) = k+1
                found_den_mld = .true.
                exit
             end if
          end do

          ! if the mixed layer depth is not found, it is set to the depth of the bottom most level
          if(.not. found_den_mld) then
             dThreshMLD(iCell) = depth(maxLevelCell(iCell))
             indMLD(iCell) = maxLevelCell(iCell)
          end if
       end do !iCell
       !$omp end do
       !$omp end parallel

       deallocate(depth)

    end if !end mld calculation choice

end subroutine ocn_eddy_compute_mixed_layer_depth

!***********************************************************************
!
!  routine interp_mld
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

   subroutine interp_mld(y0,y1,x0,x1,xT,interp_f,yT,xm1,ym1)

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

          coeffs(2) = (y1-y0)/(x1-x0)
          coeffs(1) = y0 - coeffs(2)*x0

     yT = coeffs(4)*xT**3 + coeffs(3)*xT**2 + coeffs(2)*xT + coeffs(1)
   end subroutine interp_mld

end module ocn_eddy_parameterization_helpers
