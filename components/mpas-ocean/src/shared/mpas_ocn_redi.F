! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
module ocn_redi

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_threading

   use ocn_mesh
   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables

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

   public :: ocn_redi_compute_isopycnal_slopes, &
             ocn_redi_init

   real(kind=RKIND) :: slopeTaperFactor, sfcTaperFactor
   real(kind=RKIND) :: cavityFix, oldLimiter

!***********************************************************************
!
!  routine ocn_redi_compute_isopycnal_slopes
!
!> \brief   Computes slope triads for redi mixing
!> \details
!>  This routine is the main driver for the slope triad scheme from
!
!   Isoneutral diffusion in a z-coordinate model - by Griffies et al (1998)
!   https://doi.org/10.1175/1520-0485(1998)028%3C0805:IDIAZC%3E2.0.CO;2
!
!   Important differences are
!     1. a correction term is added to accomadate steeply sloping model surfaces
!     2. a safe_slope function is added to limit slope when d rho/dz is small
!-----------------------------------------------------------------------

  contains

     function safe_slope(num, den) result(slope)

        real(kind=RKIND), intent(in) :: num, den
        real(kind=RKIND) :: slope

        ! compute isopycnal slopes safely wrt. finite precision
        ! small slope assumption, so don't want |s| > 1

        if (abs(den) > abs(num)) then
           slope = num / den
        else
           slope = sign(1.0_RKIND, den * num)
        end if

     end function

!***********************************************************************

   subroutine ocn_redi_compute_isopycnal_slopes(statePool, &
                                                timeLevelIn)

     !-----------------------------------------------------------------
     !
     ! input/output variables
     !
     !-----------------------------------------------------------------

     type(mpas_pool_type), intent(in) :: &
        statePool ! pool containing state variables

     integer, intent(in), optional :: &
                     timeLevelIn          ! time level for state variables

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     type(mpas_pool_type), pointer :: tracersPool
     real(kind=RKIND), dimension(:, :, :), pointer :: activeTracers
     integer, pointer :: indexTemperaturePtr, indexSalinityPtr
     integer :: indexTemperature, indexSalinity, i, iCellSelf
     real(kind=RKIND), dimension(:,:), allocatable :: dzTop, dTdzTop, dSdzTop
     real(kind=RKIND), dimension(:), allocatable :: k33Norm

     real(kind=RKIND), dimension(:,:), pointer :: layerThickness
     real(kind=RKIND), dimension(:), pointer :: ssh

     real(kind=RKIND) :: bldedge, mldedge, thickslope, scalefactor
     real(kind=RKIND) :: dcEdgeInv, areaEdge, drhoDT, drhoDS, dTdx
     real(kind=RKIND) :: dSdx, drhoDx, tempSlopeUp, tempSlopeDown
     real(kind=RKIND) :: slopeTaperUp, slopeTaperDown, sfcTaper
     real(kind=RKIND) :: zMLD, BruntVaisalaFreqTopEdge, maxN, sfcTaperTemp
     integer :: cell1, cell2, iCell, iEdge, k, timeLevel, nCells, nEdges

     real(kind=RKIND), parameter :: epsGM = 1.0e-12_RKIND
     if (present(timeLevelIn)) then
        timeLevel = timeLevelIn
     else
        timeLevel = 1
     end if

     call mpas_timer_start('redi compute slopes')

     call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
     call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)
     call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
     call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperaturePtr)
     call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinityPtr)
     indexTemperature = indexTemperaturePtr
     indexSalinity    = indexSalinityPtr

     !Initialize RediKappa to the constant coefficient value and modify if equalGM or N2_dependent
     !$omp parallel
     !$omp do schedule(runtime) private(k)
     do iEdge = 1, nEdgesAll
        do k = 1, nVertLevels
           RediKappaSfcTaper(k,iEdge) = 1.0_RKIND
        end do
     end do
     !$omp end do
     !$omp end parallel

     if ( config_redi_use_new_sfc_taper ) then
        !$omp parallel
        !$omp do schedule(runtime) private(k)
        do iEdge = 1, nEdgesAll
           do k = 1, nVertLevels
              RediKappa(k,iEdge) = config_Redi_constant_kappa*RediHorizontalTaper(iEdge)
           end do
        end do
        !$omp end do
        !$omp end parallel
     end if

     !$omp parallel
     !$omp do schedule(runtime) private(k)
     do iCell = 1, nCellsAll + 1
        do k = 1, nVertLevels + 1
           k33(k, iCell) = 0.0_RKIND
        end do
     end do
     !$omp end do
     !$omp end parallel

     allocate(dzTop(nVertLevels + 1,nCellsAll))
     allocate(dTdzTop(nVertLevels + 1,nCellsAll))
     allocate(dSdzTop(nVertLevels + 1,nCellsAll))
     allocate(k33Norm(nVertLevels + 1))

     nEdges = nEdgesAll
     if ( config_Redi_use_surface_taper ) then
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(k, cell1, cell2, mldEdge, bldEdge, thickSlope, scaleFactor)
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)

            bldEdge = 0.5_RKIND*(boundaryLayerDepth(cell1) + boundaryLayerDepth(cell2))
            mldEdge = config_redi_transition_layer_bld_fraction*bldEdge

            thickSlope = 1.0_RKIND/max(1.0E-15, 1.0_RKIND - min(config_redi_transition_layer_bld_fraction,1.0_RKIND))

            do k=minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
               scaleFactor = thickslope*(0.5_RKIND*abs(zMid(k,cell1)+zMid(k,cell2))/bldEdge - &
                             config_redi_transition_layer_bld_fraction)
               scaleFactor = min(1.0_RKIND,max(0.0_RKIND,scaleFactor))
               RediKappaSfcTaper(k,iEdge) = scaleFactor
            end do
         end do
         !$omp end do
         !$omp end parallel
     end if

     nCells = nCellsAll
     !$omp parallel
     !$omp do schedule(runtime)  &
     !$omp private(k)
     do iCell = 1, nCells
        ! prep dz, dTdz and dSdz for this column
        do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
           dzTop(k,iCell) = 0.5_RKIND*(layerThickness(k - 1, iCell) + layerThickness(k, iCell))
           dTdzTop(k,iCell) = (activeTracers(indexTemperature, k - 1, iCell) &
                              - activeTracers(indexTemperature, k, iCell)) &
                              /dzTop(k,iCell)
           dSdzTop(k,iCell) = (activeTracers(indexSalinity, k - 1, iCell) &
                              - activeTracers(indexSalinity, k, iCell)) &
                              /dzTop(k,iCell)
        end do
        dzTop(1:minLevelCell(iCell),iCell) = -1e-15_RKIND
        dTdzTop(1:minLevelCell(iCell),iCell) = -1e-15_RKIND
        dSdzTop(1:minLevelCell(iCell),iCell) = -1e-15_RKIND
        dzTop(maxLevelCell(iCell) + 1,iCell) = -1e-15_RKIND
        dTdzTop(maxLevelCell(iCell) + 1,iCell) = -1e-15_RKIND
        dSdzTop(maxLevelCell(iCell) + 1,iCell) = -1e-15_RKIND
     end do
     !$omp end do
     !$omp end parallel

     nCells = nCellsHalo(2)
     !$omp parallel
     !$omp do schedule(runtime)  &
     !$omp private(k33Norm, k, i, iEdge, cell1, cell2, &
     !$omp         iCellSelf, dcEdgeInv, areaEdge, drhoDT, drhoDS, dTdx, dSdx, drhoDx, &
     !$omp         slopeTaperUp, slopeTaperDown, tempSlopeUp, tempSlopeDown, sfcTaperTemp, &
     !$omp         sfcTaper)
     do iCell = 1, nCells
        k33(1:maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
        k33Norm(1:maxLevelCell(iCell) + 1) = epsGM

        do i = 1, nEdgesOnCell(iCell)
           iEdge = edgesOnCell(i, iCell)
           cell1 = cellsOnEdge(1, iEdge)
           cell2 = cellsOnEdge(2, iEdge)
           if (cell1 == iCell) then
              iCellSelf = 1
           else  ! cell2 == iCell
              iCellSelf = 2
           end if
           dcEdgeInv = 1.0_RKIND/dcEdge(iEdge)
           areaEdge = dcEdge(iEdge)*dvEdge(iEdge)

           do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
              drhoDT = -thermExpCoeff(k, iCell)
              drhoDS = salineContractCoeff(k, iCell)
              dTdx = (activeTracers(indexTemperature, k, cell2) &
                     - activeTracers(indexTemperature, k, cell1)) &
                     *dcEdgeInv + cavityFix*0.5_RKIND*(dTdzTop(k,cell2) + &
                      dTdzTop(k,cell1))*(zMid(k,cell2) - zMid(k,cell1))* &
                            dcEdgeInv
              dSdx = (activeTracers(indexSalinity, k, cell2) &
                      - activeTracers(indexSalinity, k, cell1)) &
                       *dcEdgeInv + cavityFix*0.5_RKIND*(dSdzTop(k,cell1) + &
                       dSdzTop(k,cell2))*(zMid(k,cell2) - zMid(k,cell1))* &
                       dcEdgeInv
              drhoDx = drhoDT*dTdx + drhoDS*dSdx

              tempSlopeUp = drhoDT*dTdzTop(k,iCell) + drhoDS*dSdzTop(k,iCell)
              tempSlopeDown = drhoDT*dTdzTop(k+1,iCell) + drhoDS*dSdzTop(k+1,iCell)

              ! Always compute *Up on the top cell and *Down on the bottom
              ! cell, even though they are never used. This avoids an if
              ! statement or separate computation for top and bottom.
              slopeTriadUp(k, iCellSelf, iEdge) = &
                       safe_slope( -drhoDx, tempSlopeUp )
              slopeTriadDown(k, iCellSelf, iEdge) = &
                       safe_slope( -drhoDx, tempSlopeDown )

              ! set taper of slope ('F' function from Danabasoglu and McWilliams 95)
              if (abs(slopeTriadDown(k, iCellSelf, iEdge)) > 0.6_RKIND*config_redi_maximum_slope) then
                 slopeTaperDown = 0.0_RKIND
              else if (abs(slopeTriadDown(k, iCellSelf, iEdge)) < 0.2_RKIND*config_redi_maximum_slope) then
                 slopeTaperDown = 1.0_RKIND
              else
                 slopeTaperDown = 0.5_RKIND*(1.0_RKIND - ((2.5_RKIND*abs(slopeTriadDown(k, iCellSelf, iEdge)))/ &
                                  config_Redi_maximum_slope - 1.0_RKIND)*(4.0_RKIND - abs(10.0_RKIND* &
                                  abs(slopeTriadDown(k, iCellSelf, iEdge))/config_Redi_maximum_slope - 4.0_RKIND)))
              end if
              if (abs(slopeTriadUp(k, iCellSelf, iEdge)) > 0.6_RKIND*config_redi_maximum_slope) then
                 slopeTaperUp = 0.0_RKIND
              else if (abs(slopeTriadUp(k, iCellSelf, iEdge)) < 0.2_RKIND*config_redi_maximum_slope) then
                 slopeTaperUp = 1.0_RKIND
              else
                 slopeTaperUp = 0.5_RKIND*(1.0_RKIND - ((2.5_RKIND*abs(slopeTriadUp(k, iCellSelf, iEdge)))/ &
                                  config_Redi_maximum_slope - 1.0_RKIND)*(4.0_RKIND - abs(10.0_RKIND* &
                                  abs(slopeTriadUp(k, iCellSelf, iEdge))/config_Redi_maximum_slope - 4.0_RKIND)))
              end if

              slopeTaperUp = 1.0_RKIND + slopeTaperFactor*(slopeTaperUp - 1.0_RKIND)
              slopeTaperDown = 1.0_RKIND + slopeTaperFactor*(slopeTaperDown - 1.0_RKIND)

              if (config_Redi_use_new_sfc_taper) then
                 sfcTaperTemp = RediKappaSfcTaper(k,iEdge)
              else
                 if (k < min( indMLD(cell1), indMLD(cell2))) then
                    sfcTaperTemp = 0.0_RKIND
                 else
                    sfcTaperTemp = 1.0_RKIND
                 end if
              end if

              sfcTaper = 1.0_RKIND + sfcTaperFactor*(sfcTaperTemp - 1.0_RKIND)

              !For Redi term 1 no limiting in the boundary layer and slope limiting below
              !in the boundary layer the first term is zero, so 1.0 is added to allow term 1 to stay active
              !the second term is 1 in the boundary layer and 0 below the boundary layer
              limiterUp(k,iCellSelf,iEdge) = (1.0_RKIND - oldlimiter)*(1.0_RKIND +  &
                                             sfcTaper*(slopeTaperUp - 1.0_RKIND)) + &
                                             oldlimiter*slopeTaperUp*sfcTaper
              limiterDown(k,iCellSelf,iEdge) = (1.0_RKIND - oldlimiter)*(1.0_RKIND + &
                                               sfcTaper*(slopeTaperDown - 1.0_RKIND)) + &
                                               oldlimiter*slopeTaperDown*sfcTaper

              !removing sfcTaperUp/Down as they are identical
              if (k > minLevelCell(iCell)) then
                 k33(k, iCell) = k33(k, iCell) +  &
                                 slopeTaperUp*sfcTaper*areaEdge*dzTop(k,iCell)*slopeTriadUp(k, iCellSelf, iEdge)**2
                 k33Norm(k) = k33Norm(k) + areaEdge*dzTop(k,iCell)
              end if

              !Taper Redi by tapering the slopes
              k33(k + 1, iCell) = k33(k + 1, iCell) +  &
                                  slopeTaperDown*sfcTaper*areaEdge*dzTop(k + 1,iCell)*slopeTriadDown(k, iCellSelf, iEdge)**2
              k33Norm(k + 1) = k33Norm(k + 1) + areaEdge*dzTop(k + 1,iCell)

              slopeTriadUp(k, iCellSelf, iEdge) = &
                       slopeTaperUp*sfcTaper*slopeTriadUp(k, iCellSelf, iEdge)
              slopeTriadDown(k, iCellSelf, iEdge) = &
                       slopeTaperDown*sfcTaper*slopeTriadDown(k, iCellSelf, iEdge)

           end do ! maxLevelEdgeTop(iEdge)
        end do ! nEdgesOnCell(iCell)

        ! Normalize k33
        do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
           k33(k, iCell) = k33(k, iCell)/k33Norm(k)
        end do
        k33(1:minLevelCell(iCell), iCell) = 0.0_RKIND
        k33(maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
     end do ! iCell
     !$omp end do
     !$omp end parallel

     deallocate (dzTop)
     deallocate (dTdzTop)
     deallocate (dSdzTop)
     deallocate (k33Norm)

     ! allow disabling of K33 for testing
     if (config_disable_redi_k33) then
        nCells = nCellsAll
        !$omp parallel
        !$omp do schedule(runtime)
        do iCell = 1, nCells
           k33(:, iCell) = 0.0_RKIND
        end do
        !$omp end do
        !$omp end parallel
     end if

     ! add a N2 dependent redi here

     if (config_Redi_closure == 'N2_dependent') then
        !$omp parallel
        !$omp do schedule(runtime) private(i, k, BruntVaisalaFreqTopEdge, maxN, cell1, cell2)
        do iEdge=1,nEdges
           k=minLevelEdgeBot(iEdge)
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           i = max(k,min(int(indexBoundaryLayerDepth(cell1)),int(indexBoundaryLayerDepth(cell2))))
           BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTop(i,cell1), &
                       BruntVaisalaFreqTop(i,cell2))
           maxN = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
           !Scan from the boundary layer depth only
           do k = i+1, maxLevelEdgeTop(iEdge)
              !According to Danabasoglu and Marshall (2007) take max BVF
              BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTop(k,cell1), &
                                            BruntVaisalaFreqTop(k,cell2))
              maxN = max(maxN,max(BruntVaisalaFreqTopEdge, 0.0_RKIND))
           enddo

           do k = i+1, maxLevelEdgeTop(iEdge)
              BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTop(k,cell1), &
                                            BruntVaisalaFreqTop(k,cell2))
              BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)

              RediKappa(k, iEdge) = min(max(config_redi_spatially_variable_min_kappa /   &
                                        config_redi_spatially_variable_max_kappa,        &
                                        BruntVaisalaFreqTopEdge/(maxN + 1.0E-10_RKIND)), &
                                                   1.0_RKIND)*config_Redi_constant_kappa
           end do
       enddo
       !$omp end do
       !$omp end parallel
     end if

     if (config_Redi_closure == 'equalGM') then
        ! Set Redi closure to be the same as the GM closure, i.e. if
        ! config_GM_closure is EdenGreatbatch or Visbeck, you will see those
        ! same closure fields in RediKappa. The horizontal tapering still
        ! comes from the config_Redi_horizontal_taper flag, not GM.
        !$omp parallel
        !$omp do schedule(runtime) private(k)
        do iEdge=1, nEdges
           do k = minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
              RediKappa(k,iEdge) = RediHorizontalTaper(iEdge) * gmBolusKappa(iEdge)
           end do
        end do
        !$omp end do
        !$omp end parallel
     else if ( config_Redi_closure == 'constant'     .and. config_Redi_horizontal_taper == 'RossbyRadius' &
               .or. &
               config_Redi_closure == 'N2_dependent' .and. config_Redi_horizontal_taper == 'RossbyRadius' &
               .or. &
               config_Redi_closure == 'N2_dependent' .and. config_Redi_horizontal_taper == 'ramp' &
             ) then
        !$omp parallel
        !$omp do schedule(runtime) private(k)
        do iEdge=1, nEdges
           do k = minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
              RediKappa(k,iEdge) = RediHorizontalTaper(iEdge) * RediKappa(k,iEdge)
           end do
        end do
        !$omp end do
        !$omp end parallel
     end if

     call mpas_timer_stop('redi compute slopes')

   end subroutine ocn_redi_compute_isopycnal_slopes

   subroutine ocn_redi_init(ierr)

     integer, intent(in) :: ierr

     oldlimiter = 1.0_RKIND
     if (config_Redi_use_new_sfc_taper) oldlimiter = 0.0_RKIND

     cavityFix = 0.0_RKIND
     if (config_Redi_use_cavity_fix) cavityFix = 1.0_RKIND

     if (config_Redi_use_slope_taper) then
        slopeTaperFactor = 1.0_RKIND
     else
        slopeTaperFactor = 0.0_RKIND
     end if

     if (config_Redi_use_surface_taper) then
        sfcTaperFactor = 1.0_RKIND
     else
        sfcTaperFactor = 0.0_RKIND
     end if

   end subroutine ocn_redi_init

end module ocn_redi
