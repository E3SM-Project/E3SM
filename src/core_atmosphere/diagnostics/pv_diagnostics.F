! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module pv_diagnostics

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
    use mpas_kind_types, only : RKIND

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: tend
    type (MPAS_pool_type), pointer :: tend_physics

    type (MPAS_clock_type), pointer :: clock

    public :: pv_diagnostics_setup, &
              pv_diagnostics_compute

    private

    logical :: need_ertel_pv, need_u_pv, need_v_pv, need_theta_pv, need_vort_pv, need_iLev_DT, &
               need_tend_lw, need_tend_sw, need_tend_bl, need_tend_cu, need_tend_mix, need_tend_mp, &
               need_tend_diab, need_tend_fric, need_tend_diab_pv, need_tend_fric_pv, need_dtheta_mp


    contains


    !-----------------------------------------------------------------------
    !  routine pv_diagnostics_setup
    !
    !> \brief Initialize the diagnostic
    !> \author 
    !> \date   
    !> \details
    !>  Initialize the diagnostic
    !
    !-----------------------------------------------------------------------
    subroutine pv_diagnostics_setup(all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type, MPAS_STREAM_OUTPUT, MPAS_STREAM_INPUT, &
                                       MPAS_STREAM_INPUT_OUTPUT
        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock


        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'tend', tend)
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)

        clock => simulation_clock
   
    end subroutine pv_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine pv_diagnostics_compute
    !
    !> \brief Compute diagnostic before model output is written
    !> \author 
    !> \date   
    !> \details
    !>  Compute diagnostic before model output is written
    !
    !-----------------------------------------------------------------------
    subroutine pv_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        logical :: need_any_diags, need_any_budget

        need_any_diags = .false.
        need_any_budget = .false.


        need_ertel_pv = MPAS_field_will_be_written('ertel_pv')
        need_any_diags = need_any_diags .or. need_ertel_pv
        need_u_pv = MPAS_field_will_be_written('u_pv')
        need_any_diags = need_any_diags .or. need_u_pv
        need_v_pv = MPAS_field_will_be_written('v_pv')
        need_any_diags = need_any_diags .or. need_v_pv
        need_theta_pv = MPAS_field_will_be_written('theta_pv')
        need_any_diags = need_any_diags .or. need_theta_pv
        need_vort_pv = MPAS_field_will_be_written('vort_pv')
        need_any_diags = need_any_diags .or. need_vort_pv
        need_iLev_DT = MPAS_field_will_be_written('iLev_DT')
        need_any_diags = need_any_diags .or. need_iLev_DT

        need_tend_lw = MPAS_field_will_be_written('depv_dt_lw')
        need_any_diags = need_any_diags .or. need_tend_lw
        need_any_budget = need_any_budget .or. need_tend_lw
        need_tend_sw = MPAS_field_will_be_written('depv_dt_sw')
        need_any_diags = need_any_diags .or. need_tend_sw
        need_any_budget = need_any_budget .or. need_tend_sw
        need_tend_bl = MPAS_field_will_be_written('depv_dt_bl')
        need_any_diags = need_any_diags .or. need_tend_bl
        need_any_budget = need_any_budget .or. need_tend_bl
        need_tend_cu = MPAS_field_will_be_written('depv_dt_cu')
        need_any_diags = need_any_diags .or. need_tend_cu
        need_any_budget = need_any_budget .or. need_tend_cu
        need_tend_mix = MPAS_field_will_be_written('depv_dt_mix')
        need_any_diags = need_any_diags .or. need_tend_mix
        need_any_budget = need_any_budget .or. need_tend_mix
        need_dtheta_mp = MPAS_field_will_be_written('dtheta_dt_mp')
        need_any_diags = need_any_diags .or. need_dtheta_mp
        need_any_budget = need_any_budget .or. need_dtheta_mp
        need_tend_mp = MPAS_field_will_be_written('depv_dt_mp')
        need_any_diags = need_any_diags .or. need_tend_mp
        need_any_budget = need_any_budget .or. need_tend_mp
        need_tend_diab = MPAS_field_will_be_written('depv_dt_diab')
        need_any_diags = need_any_diags .or. need_tend_diab
        need_any_budget = need_any_budget .or. need_tend_diab
        need_tend_fric = MPAS_field_will_be_written('depv_dt_fric')
        need_any_diags = need_any_diags .or. need_tend_fric
        need_any_budget = need_any_budget .or. need_tend_fric
        need_tend_diab_pv = MPAS_field_will_be_written('depv_dt_diab_pv')
        need_any_diags = need_any_diags .or. need_tend_diab_pv
        need_any_budget = need_any_budget .or. need_tend_diab_pv
        need_tend_fric_pv = MPAS_field_will_be_written('depv_dt_fric_pv')
        need_any_diags = need_any_diags .or. need_tend_fric_pv
        need_any_budget = need_any_budget .or. need_tend_fric_pv

        if (need_any_diags) then
            call atm_compute_pv_diagnostics(state, 1, diag, mesh)
        end if
        if (need_any_budget) then
            call atm_compute_pvBudget_diagnostics(state, 1, diag, mesh, tend, tend_physics)
        end if
   
    end subroutine pv_diagnostics_compute


   real(kind=RKIND) function dotProduct(a, b, sz)

      implicit none

      real(kind=RKIND), dimension(:), intent(in) :: a, b
      integer, intent(in) :: sz

      integer :: i
      real(kind=RKIND) :: rsum

      rsum = 0.0_RKIND

      do i=1,sz
         rsum = rsum + a(i)*b(i)
      end do

      dotProduct = rsum
   end function dotProduct

   integer function elementIndexInArray(val, array, sz)

      implicit none

      integer, intent(in) :: val
      integer, dimension(:), intent(in) :: array
      integer, intent(in) :: sz

      integer :: i, ind
      ind = -1
      do i=1,sz
         if (array(i)==val) then
            ind = i
            elementIndexInArray = ind !This returns, right?
            exit !just in case :)
         end if
      end do
      elementIndexInArray = ind
   end function elementIndexInArray
   
   real(kind=RKIND) function formErtelPV(gradxu, gradtheta, density, unitX, unitY, unitZ)

      use mpas_constants, only : omega_e => omega

      implicit none

      real(kind=RKIND), dimension(3), intent(inout) :: gradxu
      real(kind=RKIND), dimension(3), intent(in) :: gradtheta
      real(kind=RKIND), intent(in) :: density
      real(kind=RKIND), dimension(3), intent(in) :: unitX, unitY, unitZ

      real(kind=RKIND) :: epv, eVort
      real(kind=RKIND), dimension(3) :: eVortDir, eVortComponents

      !earth vorticity is in +z-direction in global Cartesian space
      eVort = 2.0 * omega_e
      eVortDir(1) = 0.0_RKIND
      eVortDir(2) = 0.0_RKIND
      eVortDir(3) = eVort

      eVortComponents(1) = dotProduct(eVortDir, unitX,3)
      eVortComponents(2) = dotProduct(eVortDir, unitY,3)
      eVortComponents(3) = dotProduct(eVortDir, unitZ,3)

      gradxu(:) = gradxu(:) + eVortComponents(:)

      epv = dotProduct(gradxu, gradtheta,3) / density

      epv = epv * 1.0e6 !SI to PVUs
    
      formErtelPV = epv
   end function formErtelPV
   
   subroutine local2FullVorticity(gradxu, unitX, unitY, unitZ)
      !given gradxu, return gradxu+earthVort
      
      use mpas_constants, only : omega_e => omega

      implicit none

      real(kind=RKIND), dimension(3), intent(inout) :: gradxu
      real(kind=RKIND), dimension(3), intent(in) :: unitX, unitY, unitZ
      
      real(kind=RKIND) :: eVort
      real(kind=RKIND), dimension(3) :: eVortDir, eVortComponents

      !earth vorticity is in z-direction in global Cartesian space
      eVort = 2.0 * omega_e
      eVortDir(1) = 0.0_RKIND
      eVortDir(2) = 0.0_RKIND
      eVortDir(3) = eVort

      eVortComponents(1) = dotProduct(eVortDir, unitX,3)
      eVortComponents(2) = dotProduct(eVortDir, unitY,3)
      eVortComponents(3) = dotProduct(eVortDir, unitZ,3)

      gradxu(:) = gradxu(:) + eVortComponents(:)
   end subroutine local2FullVorticity
   
   real(kind=RKIND) function calc_verticalVorticity_cell(c0, level, nVerticesOnCell, verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell, vVortVertex)
      !area weighted average of vorticity at vertices to cell center for the specified cell
      !
      implicit none

      real(kind=RKIND), intent(in) :: areaCell
      integer, intent(in) :: c0, level, nVerticesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnVertex
      real(kind=RKIND), dimension(:,:), intent(in) :: kiteAreasOnVertex, vVortVertex

      real(kind=RKIND) :: vVortCell
      integer :: i, iVertex, cellIndOnVertex

      vVortCell = 0.0_RKIND
      do i = 1,nVerticesOnCell
         iVertex = verticesOnCell(i,c0)
         cellIndOnVertex = elementIndexInArray(c0, cellsOnVertex(:,iVertex), 3)
         vVortCell = vVortCell + kiteAreasOnVertex(cellIndOnVertex, iVertex)*vVortVertex(level, iVertex)/areaCell
      end do

      calc_verticalVorticity_cell = vVortCell
   end function calc_verticalVorticity_cell

   subroutine coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, c0, xyz)

      implicit none

      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors
      integer, intent(in) :: c0
      real(kind=RKIND), dimension(3,3), intent(out) :: xyz

      integer :: i

      xyz(:,1) = cellTangentPlane(:,1,c0) !are these guaranteed unit vectors?
      xyz(:,2) = cellTangentPlane(:,2,c0)
      xyz(:,3) = localVerticalUnitVectors(:,c0)
      do i=1,2
         call normalizeVector(xyz(:,i), 3)
      end do
   end subroutine coordinateSystem_cell

   real(kind=RKIND) function fluxSign(c0, iEdge, cellsOnEdge)
      
      !For finite volume computations, we'll use a normal pointing out of the cell
      implicit none

      integer, intent(in) :: c0
      integer, intent(in) :: iEdge
      integer, dimension(:,:), intent(in) :: cellsOnEdge

      if (c0 == cellsOnEdge(1,iEdge)) then
         fluxSign = 1.0_RKIND
      else
         fluxSign = -1.0_RKIND
      end if
   end function fluxSign

   real(kind=RKIND) function calc_heightCellCenter(c0, level, zgrid)

      implicit none

      integer, intent(in) :: c0, level
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid

      calc_heightCellCenter = 0.5*(zgrid(level,c0)+zgrid(level+1,c0))
   end function calc_heightCellCenter

   real(kind=RKIND) function calc_heightVerticalEdge(c0, c1, level, zgrid)

      implicit none

      integer, intent(in) :: c0, c1, level
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid

      real(kind=RKIND) :: hTop, hBottom

      hTop = .5*(zgrid(level+1,c0)+zgrid(level+1,c1))
      hBottom = .5*(zgrid(level,c0)+zgrid(level,c1))

      calc_heightVerticalEdge = hTop-hBottom
   end function calc_heightVerticalEdge

   subroutine normalizeVector(vals, sz)
      !normalize a vector to unit magnitude
      implicit none

      real (kind=RKIND), dimension(:), intent(inout) :: vals
      integer, intent(in) :: sz

      integer :: i
      real (kind=RKIND) :: mag

      mag = 0.0_RKIND !sqrt(sum(squares))
      do i=1,sz
         mag = mag+vals(i)*vals(i)
      end do
      mag = sqrt(mag)
      vals(:) = vals(:)/mag
   end subroutine normalizeVector

   real(kind=RKIND) function calcVolumeCell(areaCell, nEdges, hEdge)

      implicit none

      integer, intent(in) :: nEdges
      real(kind=RKIND), intent(in) :: areaCell
      real(kind=RKIND), dimension(nEdges), intent(in) :: hEdge

      integer :: i
      real(kind=RKIND) :: avgHt, vol

      avgHt = 0.0_RKIND
      do i=1,nEdges
         avgHt = avgHt + hEdge(i)
      end do
      avgHt = avgHt/nEdges

      vol = areaCell*avgHt
      calcVolumeCell = vol
   end function calcVolumeCell

   real(kind=RKIND) function calc_horizDeriv_fv(valEdges, nNbrs, dvEdge, dhEdge, &
                                                normalEdge, unitDeriv, volumeCell)
      !normals to edges point out of cell
      implicit none

      integer, intent(in) :: nNbrs
      real(kind=RKIND), dimension(:), intent(in) :: valEdges, dvEdge, dhEdge
      real(kind=RKIND), dimension(3,nNbrs), intent(in) :: normalEdge
      real(kind=RKIND), dimension(3), intent(in) :: unitDeriv
      real(kind=RKIND), intent(in) :: volumeCell

      integer :: i
      real(kind=RKIND) :: vale, rsum, areaFace
      real(kind=RKIND), dimension(3) :: unitNormalEdge

      rsum = 0.0_RKIND
      do i=1,nNbrs
         vale = valEdges(i) !0.5 * (val0 + valNbrs(i))
         areaFace = dvEdge(i) * dhEdge(i)
         unitNormalEdge(:) = normalEdge(:,i)
         call normalizeVector(unitNormalEdge,3)
         areaFace = areaFace*dotProduct(unitNormalEdge, unitDeriv,3)  !* abs(dotProduct(unitNormalEdge, unitDeriv,3))
         rsum = rsum + vale * areaFace
      end do
      rsum = rsum / volumeCell

      calc_horizDeriv_fv = rsum
   end function calc_horizDeriv_fv

   !cell centers are halfway between w faces
   real(kind=RKIND) function calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)

      implicit none

      real(kind=RKIND), intent(in) :: val0, valp, valm, z0,zp,zm !center, plus, minus
      
      real(kind=RKIND) :: dval_dzp, dval_dzm

      !Average 1 sided differences to below and above since not equally spaced pts
      dval_dzp = calc_vertDeriv_one(valp, val0, zp-z0)
      dval_dzm = calc_vertDeriv_one(val0, valm, z0-zm)
      calc_vertDeriv_center = 0.5*(dval_dzp+dval_dzm)

   end function calc_vertDeriv_center

   real(kind=RKIND) function calc_vertDeriv_one(valp, valm, dz)
      !1 sided finite difference

      implicit none

      real(kind=RKIND), intent(in) :: valp, valm, dz

      calc_vertDeriv_one = (valp - valm) / dz

   end function calc_vertDeriv_one
   
   subroutine floodFill_strato(mesh, diag, pvuVal, stratoPV)
      !Searching down each column from TOA to find 2pvu surface is buggy with stratospheric wave breaking,
      !since will find 2 pvu at a higher level than "tropopause". This looks to be worse as mesh gets finer and vertical vorticity jumps.
      !Note that stratospheric blobs may persist for long times w/ slow mixing downstream of mountains or deep convection.
      !A few quicker fixes (make sure <2pvu for a number of layers; search down from 10PVU instead of TOA) are hacky and not robust.
      
      !To alleviate the (hopefully) pockets of wave breaking, we can flood fill from a known
      !stratosphere region (e.g., model top > 2pvu) and hopefully filter down around any trouble regions.
      !The problem w/ using only the flood fill is that strong surface PV anomalies can connect to 2pvu, 
      !and the resulting "flood-filled 2 pvu" can have sizeable areas that are just at the surface while there is clearly a tropopause above (e.g., in a cross-section).
      !To address large surface blobs, take the flood fill mask and try to go up from the surface to 10 pvu w/in column. If can, all stratosphere. Else, disconnect "surface blob".
      
      !The "output" is iLev_DT, which is the vertical index for the level >= pvuVal. If >nVertLevels, pvuVal above column. If <2, pvuVal below column.
      !Communication between blocks during the flood fill may be needed to treat some edge cases appropriately.

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
     
      implicit none
      
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      real(kind=RKIND), intent(in) :: pvuVal, stratoPV
      
      integer :: iCell, k, nChanged, iNbr, iCellNbr
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
      integer, dimension(:,:), pointer :: cellsOnCell
      
      real(kind=RKIND) :: sgnHemi, sgn
      real(kind=RKIND),dimension(:),pointer:: latCell
      real(kind=RKIND), dimension(:,:), pointer :: ertel_pv
      
      integer, dimension(:,:), allocatable :: candInStrato, inStrato
      
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'latCell', latCell)

      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
      
      allocate(candInStrato(nVertLevels, nCells+1))
      allocate(inStrato(nVertLevels, nCells+1))
      candInStrato(:,:) = 0
      inStrato(:,:) = 0
      !store whether each level above DT to avoid repeating logic. we'll use candInStrato as a isVisited marker further below.
      do iCell=1,nCells
         sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
         if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
         do k=1,nVertLevels
            sgn = ertel_pv(k,iCell)*sgnHemi-pvuVal
            if (sgn .GE. 0) candInStrato(k,iCell) = 1
         end do
      end do
      
      !seed flood fill with model top that's above DT.
      !can have model top below 2pvu (eg, tropics)
      nChanged = 0
      do iCell=1,nCells
         do k=nVertLevels-5,nVertLevels
            if (candInStrato(k,iCell) .GT. 0) then
               inStrato(k,iCell) = 1
               candInStrato(k,iCell) = 0
               nChanged = nChanged+1
            end if
         end do
      end do
      
      !flood fill from the given seeds. since I don't know enough fortran,
      !we'll just brute force a continuing loop rather than queue.
      do while(nChanged .GT. 0)
        nChanged = 0
        do iCell=1,nCells
          do k=nVertLevels,1,-1
             !update if candidate and neighbor in strato
             if (candInStrato(k,iCell) .GT. 0) then
                !nbr above
                if (k .LT. nVertLevels) then
                  if (inStrato(k+1,iCell) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end if
                
                !side nbrs
                do iNbr = 1, nEdgesOnCell(iCell)
                  iCellNbr = cellsOnCell(iNbr,iCell)
                  if (inStrato(k,iCellNbr) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end do
                
                !nbr below
                if (k .GT. 1) then
                  if (inStrato(k-1,iCell) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end if
             end if !candInStrato
          end do !levels
        end do !cells
      end do !while
      
      !Detach high surface PV blobs w/o vertical connection to "stratosphere"
      do iCell=1,nCells
        if (inStrato(1,iCell) .GT. 0) then
          !see how high up we can walk in the column
          do k=2,nVertLevels
            if (inStrato(k,iCell) .LT. 1) then
              exit
            end if !k is highest connected level to sfc
            sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
            if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
            sgn = ertel_pv(k,iCell)*sgnHemi-stratoPV
            if (sgn .LT. 0) then !not actually connected to "stratosphere"
              inStrato(1:k,iCell) = 0
            end if
          end do !k
        end if !inStrato at sfc
      end do !iCell
      
      !Fill iLev_DT with the lowest level above the tropopause (If DT above column, iLev>nVertLevels. If DT below column, iLev=0.
      nChanged = 0
      do iCell=1,nCells
        do k=1,nVertLevels
          if (inStrato(k,iCell) .GT. 0) then
            nChanged = 1
            exit
          end if
        end do !k
        if (nChanged .GT. 0) then !found lowest level
          if (k .EQ. 1) then 
            sgnHemi = sign(1.0_RKIND, latCell(iCell))
            sgn = ertel_pv(k,iCell)*sgnHemi-pvuVal
            if (sgn .GT. 0) then !whole column above DT
              iLev_DT(iCell) = 0
            end if
          else
            iLev_DT(iCell) = k
          end if
        else !whole column below DT
          iLev_DT(iCell) = nVertLevels+2
        end if
      end do !iCell
     
   end subroutine floodFill_strato
   
   subroutine interp_pv_diagnostics(mesh, diag, pvuVal, missingVal)
      !compute various fields on 2pvu surface using calculated PVU field
      !potential temperature, uZonal, uMeridional, vertical vorticity

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      
      implicit none
      
      type (mpas_pool_type), intent(in)  :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      real(kind=RKIND) ::  pvuVal, missingVal
      
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
                                          
      real(kind=RKIND),dimension(:),pointer:: areaCell, latCell, u_pv, v_pv, theta_pv, vort_pv
      real(kind=RKIND),dimension(:,:),pointer:: uZonal, uMeridional, vorticity, theta, ertel_pv, &
                                                kiteAreasOnVertex
      real(kind=RKIND), dimension(:,:), allocatable :: vVort
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCells)
      
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'uReconstructZonal', uZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uMeridional)
      call mpas_pool_get_array(diag, 'u_pv', u_pv)
      call mpas_pool_get_array(diag, 'v_pv', v_pv)
      call mpas_pool_get_array(diag, 'theta_pv', theta_pv)
      call mpas_pool_get_array(diag, 'vort_pv', vort_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
      
      !write(0,*) 'Interpolating u,v,theta,vort to pv '
      
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, uZonal, u_pv, missingVal, iLev_DT)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, uMeridional, v_pv, missingVal, iLev_DT)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, theta, theta_pv, missingVal, iLev_DT)
                     
      allocate(vVort(nVertLevels, nCells+1))
      do iCell=1,nCells
         do k=1,nVertLevels
            vVort(k,iCell) = calc_verticalVorticity_cell(iCell, k, nEdgesOnCell(iCell), verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell(iCell), vorticity)
         end do
      end do
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, vVort, vort_pv, missingVal, iLev_DT)
      deallocate(vVort)
      !write(0,*) 'Done interpolating '
   end subroutine interp_pv_diagnostics     
   
   subroutine interp_pvBudget_diagnostics(mesh, diag, pvuVal, missingVal)
      !compute various fields on 2pvu surface using calculated PVU field
      !tend_diab, tend_fric

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      
      implicit none
      
      type (mpas_pool_type), intent(in)  :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      real(kind=RKIND) ::  pvuVal, missingVal
      
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: iLev_DT
                                          
      real(kind=RKIND),dimension(:),pointer:: latCell, depv_dt_diab_pv, depv_dt_fric_pv
      real(kind=RKIND),dimension(:,:),pointer:: depv_dt_diab, depv_dt_fric, ertel_pv
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCells)
      
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
      call mpas_pool_get_array(diag, 'depv_dt_diab_pv', depv_dt_diab_pv)
      call mpas_pool_get_array(diag, 'depv_dt_fric_pv', depv_dt_fric_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
      
      !write(0,*) 'Interpolating u,v,theta,vort to pv '
      
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, depv_dt_diab, depv_dt_diab_pv, missingVal, iLev_DT)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, depv_dt_fric, depv_dt_fric_pv, missingVal, iLev_DT)
      !write(0,*) 'Done interpolating '
   end subroutine interp_pvBudget_diagnostics
   
   subroutine interp_pv( nCells, nLevels, interpVal, &
                         latCell, field0, field1,field, &
                         missingVal, iLev_DT)

      implicit none
      !linear-in-PV interpolate columns of field1 to where field0 is interpVal*sign(lat)
      !using level above tropopause already diagnosed
      
      ! input

      integer :: nCells, nLevels
      integer, intent(in) :: iLev_DT(nCells)
      real(kind=RKIND) ::  interpVal, missingVal
      real(kind=RKIND), intent(in) ::latCell(nCells)
      real(kind=RKIND), intent(in) :: field0(nLevels,nCells), field1(nLevels,nCells)
      real(kind=RKIND), intent(out) :: field(nCells)

      !  local
      
      integer :: iCell, iLev, levInd, indlNbr
      real(kind=RKIND) :: valh, vall, vallNbr, sgnh, sgnl, sgnlNbr
      real(kind=RKIND) :: dv_dl, levFrac, valInterpCell, sgnHemi

      do iCell = 1, nCells
        !starting from top, trap val if values on opposite side
        levInd = -1 !what should happen with missing values?
        levFrac = 0.0
        sgnHemi = sign(1.0_RKIND, latCell(iCell)) !problem at the equator...is sign(0)=0?
        if (sgnHemi .EQ. 0.0) sgnHemi = 1.0
        valInterpCell = interpVal*sgnHemi
        
        iLev = iLev_DT(iCell)
        if (iLev .GT. nLevels) then
          levInd = -1
          sgnl = -1.0
        else if (iLev .LT. 1) then
          levInd = -1
          sgnl = 1.0
        else
          valh = field0(iLev,iCell)
          vall = field0(iLev-1,iCell)
          !sandwiched value. equal in case val0 is a vals[l].
          !get linear interpolation: val0 = vals[l]+dvals/dl * dl
          !Avoid divide by 0 by just assuming value is 
          !halfway between...
   
          dv_dl = valh-vall;
          if (abs(dv_dl)<1.e-6) then
            levFrac = 0.5;
          else
            levFrac = (valInterpCell-vall)/dv_dl
          end if
          
          levInd = iLev-1
        end if !iLev in column

        !find value of field using index we just found
        if (levInd<0) then !didn't trap value
          if (sgnl>0.0) then !column above value, take surface
            field(iCell) = field1(1,iCell)
          else !column below value, take top
            !field(iCell) = missingVal
            field(iCell) = field1(nLevels,iCell)
          end if
        else
          valh = field1(levInd+1,iCell)
          vall = field1(levInd,iCell)
        
          dv_dl = valh-vall
          field(iCell) = vall+dv_dl*levFrac
        end if
      end do
      
   end subroutine interp_pv
   
   subroutine calc_gradxu_cell(gradxu, addEarthVort, &
                             iCell, level, nVertLevels, nEdgesCell0, verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell0, &
                             uReconstructX, uReconstructY, uReconstructZ, w,vorticity)
      implicit none
      
      real(kind=RKIND), dimension(3), intent(out) :: gradxu
      integer, intent(in) :: addEarthVort, iCell, level, nVertLevels, nEdgesCell0
      real(kind=RKIND), intent(in) :: areaCell0
      real(kind=RKIND), dimension(:), intent(in) :: dvEdge
      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors, edgeNormalVectors
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid,uReconstructX, uReconstructY, uReconstructZ, &
                                                      w, vorticity, kiteAreasOnVertex
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      
      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: val0, valNbr, volumeCell, areaFactor, z0, zp, zm, valp, valm, dw_dx, dw_dy, du_dz, dv_dz
      real(kind=RKIND), dimension(3) :: unitDeriv, velCell0, velCellp, velCellm
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(nEdgesCell0) :: valEdges, dvEdgeCell, dhEdge
      real(kind=RKIND), dimension(3,nEdgesCell0) :: normalEdgeCell
     
     !local coordinate system
      call coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, iCell, xyzLocal)
      !normal vectors at voronoi polygon edges pointing out of cell
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         !dhEdge(i) = calc_heightVerticalEdge(iCell, iNbr, level, zgrid) !vertical thickness of that face
         !if don't want to consider 3d cell since we haven't calculated the cell
         !volume well, set all thicknesses to be the same
         dhEdge(i) = 100.0_RKIND

         iEdge = edgesOnCell(i,iCell)
         dvEdgeCell(i) = dvEdge(iEdge)
         val0 = fluxSign(iCell, iEdge, cellsOnEdge)
         normalEdgeCell(:,i) = edgeNormalVectors(:,iEdge)
         call normalizeVector(normalEdgeCell(:,i),3)
         normalEdgeCell(:,i) = normalEdgeCell(:,i)*val0
      end do

      volumeCell = calcVolumeCell(areaCell0, nEdgesCell0, dhEdge)
      
      !w
      val0 = .5*(w(level+1, iCell)+w(level, iCell))
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = .5*(w(level+1, iNbr)+w(level, iNbr))
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      dw_dx = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      dw_dy = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !vertical derivatives
      !calc_heightCellCenter(c0, level, zgrid) calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      !du/dz and dv/dz
      velCell0(1) = uReconstructX(level,iCell)
      velCell0(2) = uReconstructY(level,iCell)
      velCell0(3) = uReconstructZ(level,iCell)
      z0 = calc_heightCellCenter(iCell, level, zgrid)
      if (level>1) then
         !have cell beneath
         velCellm(1) = uReconstructX(level-1,iCell)
         velCellm(2) = uReconstructY(level-1,iCell)
         velCellm(3) = uReconstructZ(level-1,iCell)
         zm = calc_heightCellCenter(iCell, level-1, zgrid)
      end if
      if (level<nVertLevels) then
         !have cell above
         velCellp(1) = uReconstructX(level+1,iCell)
         velCellp(2) = uReconstructY(level+1,iCell)
         velCellp(3) = uReconstructZ(level+1,iCell)
         zp = calc_heightCellCenter(iCell, level+1, zgrid)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valp = dotProduct(velCellp, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_one(valp, val0, zp-z0)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valp = dotProduct(velCellp, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valm = dotProduct(velCellm, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_one(val0, valm, z0-zm)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valm = dotProduct(velCellp, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_one(val0, valm, z0-zm)
      else
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valp = dotProduct(velCellp, xyzLocal(:,1),3)
         valm = dotProduct(velCellm, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valp = dotProduct(velCellp, xyzLocal(:,2),3)
         valm = dotProduct(velCellm, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if

      gradxu(3) = calc_verticalVorticity_cell(iCell, level, nEdgesCell0, verticesOnCell, cellsOnVertex, &
                                              kiteAreasOnVertex, areaCell0, vorticity)

      gradxu(1) = dw_dy-dv_dz
      gradxu(2) = du_dz-dw_dx
      
      if (addEarthVort>0) then
        call local2FullVorticity(gradxu, xyzLocal(:,1), xyzLocal(:,2), xyzLocal(:,3))
      end if
     
   end subroutine calc_gradxu_cell
   
   subroutine calc_grad_cell(gradtheta, &
                             iCell, level, nVertLevels, nEdgesCell0, verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell0, &
                             theta)
      !
      implicit none
      
      real(kind=RKIND), dimension(3), intent(out) :: gradtheta
      real(kind=RKIND), intent(in) :: areaCell0
      real(kind=RKIND), dimension(:), intent(in) :: dvEdge
      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors, edgeNormalVectors
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid, theta, kiteAreasOnVertex
      integer, intent(in) :: iCell, level, nVertLevels, nEdgesCell0
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      
      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: val0, valNbr, volumeCell, areaFactor, z0, zp, zm, valp, valm
      real(kind=RKIND), dimension(3) :: unitDeriv, velCell0, velCellp, velCellm
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(nEdgesCell0) :: valEdges, dvEdgeCell, dhEdge
      real(kind=RKIND), dimension(3,nEdgesCell0) :: normalEdgeCell

      !local coordinate system
      call coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, iCell, xyzLocal)
      !normal vectors at voronoi polygon edges pointing out of cell
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         !dhEdge(i) = calc_heightVerticalEdge(iCell, iNbr, level, zgrid) !vertical thickness of that face
         !if don't want to consider 3d cell since we haven't calculated the cell
         !volume well, set all thicknesses to be the same
         dhEdge(i) = 100.0_RKIND

         iEdge = edgesOnCell(i,iCell)
         dvEdgeCell(i) = dvEdge(iEdge)
         val0 = fluxSign(iCell, iEdge, cellsOnEdge)
         normalEdgeCell(:,i) = edgeNormalVectors(:,iEdge)
         call normalizeVector(normalEdgeCell(:,i),3)
         normalEdgeCell(:,i) = normalEdgeCell(:,i)*val0
      end do

      volumeCell = calcVolumeCell(areaCell0, nEdgesCell0, dhEdge)

      !Need to get 3d curl and grad theta
      !horizontal derivatives
      !calc_horizDeriv_fv(valEdges, nNbrs, dvEdge, dhEdge, &
      !                                         normalEdge, unitDeriv, volumeCell)
      !theta
      val0 = theta(level, iCell)
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = theta(level,iNbr)
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      gradtheta(1) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      gradtheta(2) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !vertical derivatives
      !calc_heightCellCenter(c0, level, zgrid) calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      !theta
      gradtheta(3) = 0.0_RKIND
      z0 = calc_heightCellCenter(iCell, level, zgrid)
      val0 = theta(level, iCell)
      if (level>1) then
         !have cell beneath
         valm = theta(level-1, iCell)
         zm = calc_heightCellCenter(iCell, level-1, zgrid)
      end if
      if (level<nVertLevels) then
         !have cell above
         valp = theta(level+1, iCell)
         zp = calc_heightCellCenter(iCell, level+1, zgrid)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         gradtheta(3) = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         gradtheta(3) = calc_vertDeriv_one(val0, valm, z0-zm)
      else
         gradtheta(3) = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if
   
   end subroutine calc_grad_cell
   
   subroutine calc_vertical_curl(vorticity, u, dcEdge, areaTriangle, verticesOnEdge, nEdges, nVertices)
      ! Adapted from computation of circulation and relative vorticity at each vertex in atm_compute_solve_diagnostics()
      !This takes scvt face values and computes finite volume curl at scvt vertices (triangle cell centers), but
      !only works on 1 horizontal level at a time
      
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: vorticity
      integer, intent(in) :: nEdges, nVertices
      integer, dimension(:,:), intent(in) :: verticesOnEdge
      real (kind=RKIND), dimension(:), intent(in) :: dcEdge, areaTriangle
      real (kind=RKIND), dimension(:), intent(in) :: u
      
      integer :: iEdge, iVertex
      
      vorticity(:) = 0.0
      do iEdge=1,nEdges
            vorticity(verticesOnEdge(1,iEdge)) = vorticity(verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(iEdge)
            vorticity(verticesOnEdge(2,iEdge)) = vorticity(verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(iEdge)
      end do
      do iVertex=1,nVertices
            vorticity(iVertex) = vorticity(iVertex) / areaTriangle(iVertex)
      end do

   end subroutine calc_vertical_curl
   
   subroutine calc_epv(mesh, time_lev, state, diag)
      
      !EPV= absoluteVorticity/density . grad(theta)

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      
      implicit none
      
      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh

      integer :: iCell, k
      integer, pointer :: nCellsSolve, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
      !real(kind=RKIND) :: rvord
      real(kind=RKIND), dimension(3) :: gradxu, gradtheta
      real(kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real(kind=RKIND), dimension(:,:), pointer :: w, rho, vorticity, zgrid, &
                                                   localVerticalUnitVectors, edgeNormalVectors, kiteAreasOnVertex, &
                                                   theta, uReconstructX, uReconstructY, uReconstructZ, &
                                                   ertel_pv
      real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellTangentPlane', cellTangentPlane)
      call mpas_pool_get_array(mesh, 'localVerticalUnitVectors', localVerticalUnitVectors)
      call mpas_pool_get_array(mesh, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      
      call mpas_pool_get_array(state, 'w', w, time_lev)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      
      !epv and diabatic component ----------------------
      do iCell=1,nCellsSolve
         do k=1,nVertLevels
            !vort1/rho1
            call calc_gradxu_cell(gradxu, 1, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             uReconstructX, uReconstructY, uReconstructZ, w,vorticity)
            gradxu(:) = gradxu(:)/rho(k,iCell)
            
            !epv1 -------------
            call calc_grad_cell(gradtheta, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             theta)

            ertel_pv(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6 !SI to PVUs
         end do
      end do
   end subroutine calc_epv
   
   subroutine atm_compute_pv_diagnostics(state, time_lev, diag, mesh)
   ! diagnose epv and some fields on horizontal surfaces
   
      use mpas_constants
      use mpas_derived_types, only : field2DReal
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
      use mpas_dmpar, only : mpas_dmpar_exch_halo_field
   
      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal, stratoPV
      real (kind=RKIND), dimension(:,:), pointer :: theta, rho, theta_m, rho_zz, zz, dtheta_dt_mix, tend_theta_euler
      type (field2DReal), pointer :: theta_f, uReconstructX_f, uReconstructY_f, uReconstructZ_f, w_f
      type (field2DReal), pointer :: rthratenlw_f, rthratensw_f, rthcuten_f, rthblten_f, dtheta_dt_mp_f, theta_euler_f, dtheta_dt_mix_f
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
      call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)

      call mpas_pool_get_array(mesh, 'zz', zz)

      do iCell=1,nCells
         do k=1,nVertLevels
            theta(k,iCell) = theta_m(k,iCell) / (1._RKIND + rvord * scalars(index_qv,k,iCell))
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
         end do
      end do
      
      !nick szapiro
      write(0,*) "Calculating epv"
      
      !need halo cells for everything w/ horizontal derivative
      call mpas_pool_get_field(state, 'w', w_f, time_lev)
      call mpas_pool_get_field(diag, 'uReconstructX', uReconstructX_f)
      call mpas_pool_get_field(diag, 'uReconstructY', uReconstructY_f)
      call mpas_pool_get_field(diag, 'uReconstructZ', uReconstructZ_f)
      call mpas_pool_get_field(diag, 'theta', theta_f)

      call mpas_dmpar_exch_halo_field(theta_f)
      call mpas_dmpar_exch_halo_field(uReconstructX_f)
      call mpas_dmpar_exch_halo_field(uReconstructY_f)
      call mpas_dmpar_exch_halo_field(uReconstructZ_f)
      call mpas_dmpar_exch_halo_field(w_f)
      
      call calc_epv(mesh, time_lev, state, diag)
      
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      stratoPV = 10.0_RKIND
      call floodFill_strato(mesh, diag, pvuVal, stratoPV)
      call interp_pv_diagnostics(mesh, diag, pvuVal, missingVal)
   
   end subroutine atm_compute_pv_diagnostics
   
   subroutine calc_pvBudget(state, time_lev, diag, mesh, tend, tend_physics)
      
      use mpas_vector_reconstruction
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      
      implicit none
      
      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: tend_physics
      type (mpas_pool_type), intent(inout) :: tend !modify tend_w_euler to uncouple with density

      integer :: iCell, k, iEdge
      integer, pointer :: nCellsSolve, nVertLevels, nVertices, nCells, nEdges
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex, verticesOnEdge
      !real(kind=RKIND) :: rvord
      real(kind=RKIND), dimension(3) :: gradxu, gradtheta, gradxf
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, areaTriangle, dcEdge
      real(kind=RKIND), dimension(:,:), pointer :: w, rho, vorticity, zgrid, &
                                                   localVerticalUnitVectors, edgeNormalVectors, kiteAreasOnVertex, &
                                                   theta, uReconstructX, uReconstructY, uReconstructZ
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw, depv_dt_sw, depv_dt_bl, depv_dt_cu, depv_dt_mp, depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_diab, depv_dt_fric
      real(kind=RKIND), dimension(:,:), pointer :: tend_u_phys, tend_u_euler, rho_edge, tend_w_euler
      real(kind=RKIND), dimension(:,:), pointer :: rthblten, rthcuten, rthratenlw, rthratensw, &
                                                   dtheta_dt_mp, dtheta_dt_mix
      real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
      
      real(kind=RKIND), dimension(:,:), allocatable :: varVerts, tenduX, tenduY, tenduZ, tenduZonal,tendUMerid
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellTangentPlane', cellTangentPlane)
      call mpas_pool_get_array(mesh, 'localVerticalUnitVectors', localVerticalUnitVectors)
      call mpas_pool_get_array(mesh, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      
      call mpas_pool_get_array(state, 'w', w, time_lev)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      
      call mpas_pool_get_array(tend_physics, 'rthblten', rthblten)
      call mpas_pool_get_array(tend_physics, 'rthcuten', rthcuten)
      call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
      call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)
      call mpas_pool_get_array(diag, 'dtheta_dt_mp', dtheta_dt_mp)
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)
      
      call mpas_pool_get_array(diag, 'depv_dt_lw', depv_dt_lw)
      call mpas_pool_get_array(diag, 'depv_dt_sw', depv_dt_sw)
      call mpas_pool_get_array(diag, 'depv_dt_bl', depv_dt_bl)
      call mpas_pool_get_array(diag, 'depv_dt_cu', depv_dt_cu)
      call mpas_pool_get_array(diag, 'depv_dt_mp', depv_dt_mp)
      call mpas_pool_get_array(diag, 'depv_dt_mix', depv_dt_mix)
      
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
      
      call mpas_pool_get_array(diag , 'tend_u_phys', tend_u_phys)
      call mpas_pool_get_array(diag , 'rho_edge', rho_edge)
      call mpas_pool_get_array(tend, 'u_euler', tend_u_euler)
      call mpas_pool_get_array(tend, 'w_euler', tend_w_euler)
      
      allocate(varVerts(nVertLevels,nVertices+1))
      allocate(tenduX(nVertLevels,nCells+1))
      allocate(tenduY(nVertLevels,nCells+1))
      allocate(tenduZ(nVertLevels,nCells+1))
      allocate(tenduZonal(nVertLevels,nCells+1))
      allocate(tenduMerid(nVertLevels,nCells+1))
      
      !diabatic component ----------------------
      do iCell=1,nCellsSolve
         do k=1,nVertLevels
            !vort1/rho1
            call calc_gradxu_cell(gradxu, 1, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             uReconstructX, uReconstructY, uReconstructZ, w,vorticity)
            gradxu(:) = gradxu(:)/rho(k,iCell)
            
            !depv_dt_lw/sw/mp/ -------------
            !absolute vorticity here should maybe be from before taking timestep (from field that generated that tendency...)
            if (associated(rthratenlw)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthratenlw)
               depv_dt_lw(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6 !SI to PVUs
            else
               depv_dt_lw(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(rthratensw)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthratensw)
               depv_dt_sw(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_sw(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(rthblten)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthblten)
               depv_dt_bl(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_bl(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(rthcuten)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthcuten)
               depv_dt_cu(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_cu(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(dtheta_dt_mp)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                dtheta_dt_mp)
               depv_dt_mp(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_mp(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(dtheta_dt_mp)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                dtheta_dt_mix)
               depv_dt_mix(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_mix(k,iCell) = 0.0_RKIND
            end if
         end do 
      end do
      depv_dt_diab = depv_dt_lw + depv_dt_sw + depv_dt_bl + depv_dt_cu + depv_dt_mp + depv_dt_mix
      
      !frictional component ----------------------
      !vertical curl at vertices ( like SAT analogies tend_u:varVerts :: u:vorticity)
      do iEdge=1,nEdges
         do k=1,nVertLevels
            tend_u_phys(k,iEdge) = tend_u_phys(k,iEdge)+tend_u_euler(k,iEdge)/rho_edge(k,iEdge)
         end do
      end do
      !tend_u_phys = tend_u_phys + tend_u_euler/rho_edge
      do k=1,nVertLevels
         call calc_vertical_curl(varVerts(k,:), tend_u_phys(k,:), dcEdge, areaTriangle, verticesOnEdge, nEdges, nVertices)
      end do
      
      !tend_u at cell centers
      call mpas_reconstruct(mesh, tend_u_phys,         &
                               tenduX,tenduY,tenduZ,   &
                               tenduZonal,tenduMerid)
      !uncouple tend_w_euler
      do iCell=1,nCells
         do k=2,nVertLevels
            !average density to vertical interfaces between cells
            !top of lowest cell is interface 2
            tend_w_euler(k,iCell) = tend_w_euler(k,iCell)/( .5*( rho(k-1,iCell)+rho(k,iCell) ) )
         end do
      end do
      !constant extrapolation for density above and below cell centers
      tend_w_euler(1,1:nCells) = tend_w_euler(1,1:nCells)/rho(1,1:nCells)
      tend_w_euler(nVertLevels+1,1:nCells) = tend_w_euler(nVertLevels+1,1:nCells)/rho(nVertLevels,1:nCells)
      
      do iCell=1,nCellsSolve
         do k=1,nVertLevels
            !calculating grad(theta)/rho . (grad x F/rho)
            
            !gradtheta term
            call calc_grad_cell(gradtheta, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             theta)
            !
            gradtheta(:) = gradtheta(:)/rho(k,iCell)
            
            !we can call calc_gradxu_cell where:
            !w: tend_w     uReconstruct{X,Y,Z}: tend_u to cell centers     vorticity: tend_u at vertices
            !
            call calc_gradxu_cell(gradxf, 0, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             tenduX, tenduY,tenduZ, tend_w_euler,varVerts)
            
            depv_dt_fric(k,iCell) = dotProduct(gradxf, gradtheta,3)* 1.0e6
         end do
      end do
      
      deallocate(varVerts)
      deallocate(tenduX)
      deallocate(tenduY)
      deallocate(tenduZ)
      deallocate(tenduZonal)
      deallocate(tenduMerid)
      
   end subroutine calc_pvBudget
   
   subroutine atm_compute_pvBudget_diagnostics(state, time_lev, diag, mesh, tend, tend_physics)
      ! after calculating epv,
      !pv budget in the "classic" formulation: (e.g., Pedlosky) 
      !Depv_Dt = Thermo+Friction = vort3d/rho . grad(Dtheta/Dt) + gradTheta/rho . grad x F/rho
      ! The thermo term gets calculated just like epv but theta replaced w/ Dtheta/Dt
      ! F/rho is tend_{u,v,w} and we'll calculate a cell's vertical and horizontal curl separately.
      
      use mpas_constants
      use mpas_derived_types, only : field2DReal
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
      use mpas_dmpar, only : mpas_dmpar_exch_halo_field
      
      implicit none
      
      type (mpas_pool_type), intent(inout) :: diag, tend
      type (mpas_pool_type), intent(in) :: state, mesh, tend_physics
      integer, intent(in) :: time_lev            ! which time level to use from state
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal
      real (kind=RKIND), dimension(:,:), pointer :: dtheta_dt_mix, tend_theta_euler
      type (field2DReal), pointer :: rthratenlw_f, rthratensw_f, rthcuten_f, rthblten_f, dtheta_dt_mp_f, theta_euler_f, dtheta_dt_mix_f
      type (field2DReal), pointer :: tend_u_phys_f, tend_u_euler_f, tend_w_euler_f
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
      
      !nick szapiro
      write(0,*) "Calculating pvBudget"
      
      !need halo cells for everything w/ horizontal derivative
      !Dtheta/Dt
      call mpas_pool_get_array(tend, 'theta_euler', tend_theta_euler)
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)
      do iCell=1,nCells
         do k=1,nVertLevels
            !with modified moist potential temperature being the model state variable being mixed,
            ! assume qv field is not mixed and so there's no tend_qv to consider
            dtheta_dt_mix(k,iCell) = tend_theta_euler(k,iCell)/( 1._RKIND + rvord*scalars(index_qv,k,iCell) )
         end do
      end do
      call mpas_pool_get_field(tend_physics, 'rthratenlw', rthratenlw_f)
      call mpas_pool_get_field(tend_physics, 'rthratensw', rthratensw_f)
      call mpas_pool_get_field(tend_physics, 'rthcuten', rthcuten_f)
      call mpas_pool_get_field(tend_physics, 'rthblten', rthblten_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_mp', dtheta_dt_mp_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_mix', dtheta_dt_mix_f)
      
      call mpas_dmpar_exch_halo_field(rthratenlw_f)
      call mpas_dmpar_exch_halo_field(rthratensw_f)
      call mpas_dmpar_exch_halo_field(rthcuten_f)
      call mpas_dmpar_exch_halo_field(rthblten_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mp_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mix_f)
      
      !friction
      call mpas_pool_get_field(diag , 'tend_u_phys', tend_u_phys_f)
      call mpas_pool_get_field(tend, 'u_euler', tend_u_euler_f)
      call mpas_pool_get_field(tend, 'w_euler', tend_w_euler_f)
      call mpas_dmpar_exch_halo_field(tend_u_phys_f)
      call mpas_dmpar_exch_halo_field(tend_u_euler_f)
      call mpas_dmpar_exch_halo_field(tend_w_euler_f)
      
      call calc_pvBudget(state, time_lev, diag, mesh, tend, tend_physics)
      
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      call interp_pvBudget_diagnostics(mesh, diag, pvuVal, missingVal)
   
   end subroutine atm_compute_pvBudget_diagnostics

end module pv_diagnostics
