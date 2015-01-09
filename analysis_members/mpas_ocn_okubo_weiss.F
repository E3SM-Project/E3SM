! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_okubo_weiss
!
!> \brief MPAS ocean analysis core member: okubo_weiss
!> \author Andre Schmeisser
!> \date   August 2014
!> \details
!>  MPAS ocean analysis core member: okubo_weiss
!
!-----------------------------------------------------------------------

module Sort
   use iso_c_binding
   use mpas_grid_types

   implicit none

#ifdef SINGLE_PRECISION
   integer, parameter :: C_REAL = C_FLOAT
   integer, parameter :: SIZE_REAL = 4
#else
   integer, parameter :: C_REAL = C_DOUBLE
   integer, parameter :: SIZE_REAL = 8
#endif

   interface
      subroutine qsort(array, elem_count, elem_size, compare) bind(C, name="qsort")
         import
         type(c_ptr), value       :: array
         integer(c_size_t), value :: elem_count
         integer(c_size_t), value :: elem_size
         type(c_funptr), value    :: compare ! int (*compare)(const void*, const void*)
      end subroutine qsort ! standard C library qsort
   end interface

   public :: compareIntDescending, compareRealDescending

   contains

   integer(c_int) function compareIntDescending(a, b) bind(C)
      integer(c_int), intent(in) :: a, b

      if (a < b) then
         compareIntDescending = 1
      else if (a == b) then
         compareIntDescending = 0
      else
         compareIntDescending = -1
      end if
   end function compareIntDescending

   integer(c_int) function compareRealDescending( a, b ) bind(C)
      real (C_REAL), intent(in) :: a, b

      if (a < b) then
         compareRealDescending = 1
      else if (a == b) then
         compareRealDescending = 0
      else
         compareRealDescending = -1
      end if
   end function compareRealDescending

end module Sort



module ocn_okubo_weiss

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping

   use ocn_constants
   use ocn_diagnostics_routines
   use mpas_tensor_operations
   use mpas_matrix_operations
   use iso_c_binding

   implicit none
   private
   save

   interface
      subroutine compute_ev_2(A, wr, wi)
         use iso_c_binding, only: c_double
         real (c_double), dimension(2,2) :: A
         real (c_double), dimension(2) :: wr
         real (c_double), dimension(2) :: wi
      end subroutine compute_ev_2
   end interface

   interface
      subroutine compute_ev_3(A, wr, wi)
         use iso_c_binding, only: c_double
         real (c_double), dimension(3,3) :: A
         real (c_double), dimension(3) :: wr
         real (c_double), dimension(3) :: wi
      end subroutine compute_ev_3
   end interface

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

   public :: ocn_setup_packages_okubo_weiss, &
             ocn_init_okubo_weiss, &
             ocn_compute_okubo_weiss, &
             ocn_write_okubo_weiss, &
             ocn_restart_okubo_weiss, &
             ocn_finalize_okubo_weiss

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: amOkuboWeissTimer
   type (timer_node), pointer :: amOWCCTimer, &
      CCLocalTimer, CCComLoopTimer, CCStatsTimer, CCStatsPerProcT, &
      CCAggregateStatsT, CCAggGatherT, CCAggSortT, CCOutputT
   integer :: nCellsGlobal
   integer :: nTotalCellsGlobal
   logical :: statsFileCreated


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_setup_packages_okubo_weiss
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_okubo_weiss(configPool, packagePool, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

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
      logical, pointer :: amOkuboWeissActive

      err = 0

      call mpas_pool_get_package(packagePool, 'amOkuboWeissActive', amOkuboWeissActive)

      ! turn on package for this analysis member
      amOkuboWeissActive = .true.

   end subroutine ocn_setup_packages_okubo_weiss!}}}


!***********************************************************************
!
!  routine ocn_init_okubo_weiss
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_okubo_weiss(domain, err)!{{{

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      integer, pointer :: nCellsSolve, nVertLevels
      integer, dimension(1) :: cellsArray, cellsSumArray

      !-----------------------------------------------------------------
      !
      ! Compute the global number of cells by doing an MPI sum over all
      ! local number of cells per domain.
      !
      ! TO DO: Replace this with a global constant to avoid this unnecessary
      !        computation
      !
      !-----------------------------------------------------------------

      cellsArray(1) = 0
      dminfo = domain % dminfo
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         cellsArray(1) = cellsArray(1)  + nCellsSolve

         block => block % next
      end do

      call mpas_dmpar_sum_int_array(dminfo, 1, cellsArray, cellsSumArray)
      nCellsGlobal = cellsSumArray(1)
      nTotalCellsGlobal = nCellsGlobal*nVertLevels

#ifdef SINGLE_PRECISION
      ! The MPI communication represents eddy IDs as reals instead
      ! of integers so the IDs don't have to be send separately from
      ! the real-valued variables.
      ! Check that the cell IDs are actually exactly representable in
      ! floating point format.
      if (nTotalCellsGlobal >= 16777216) then
         write(stderrUnit,*) "Error: Maximum number of cells is not exactly representable " // &
                             "as float. Compile in double precision or change OW eddy " // &
                             "tracking method to use integers"
         err = 1
      end if
#endif

      statsFileCreated = .false.
      err = 0

   end subroutine ocn_init_okubo_weiss!}}}


!***********************************************************************
!
!  routine ocn_compute_OW_values
!
!> \brief   Compute values of Okubo-Weiss and Lamba_2 field
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!> Computes the values for Okubo-Weiss field OW, as well as the second
!> eigenvalue of the strain rate tensor Lambda_2.
!> The OW and Lam2 values use the x/y components of the velocity gradient
!> only. On the sphere, this is only correct if the velocity gradient is
!> rotated to the local tangential plane first.
!> Lam2_R3 uses the full three dimensional tensor.
!>
!> Note that the lambda_2 values are multiplied by 4 to be on the same
!> scale as the OW value!
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_OW_values(dvel, nVertLevels, nCells, maxLevelCell, &
                                    OW, OW_norm, Lam2, Lam2_R3, Lam2_norm, &
                                    S, om, Lam1)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: dvel
      real (kind=RKIND), intent(in) :: OW_norm, Lam2_norm
      integer, intent(in) :: nVertLevels, nCells
      integer, dimension(:), intent(in) :: maxLevelCell

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: S, om, OW, Lam1, Lam2, Lam2_R3

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer k, iCell, kmax
      real (kind=RKIND) :: sn, ss
      real (kind=RKIND), dimension(2,2) :: sym2, asym2
      real (kind=RKIND), dimension(2,2) :: T2
      real (kind=RKIND), dimension(3,3) :: sym3, asym3
      real (kind=RKIND), dimension(3,3) :: T3
      real (kind=RKIND), dimension(2) :: lambda, lambdaI
      real (kind=RKIND), dimension(3) :: lambda3, lambdaI3

      ! Compute OW according to (12): OW = s_n^2 + s_s^2 - \omega^2
      ! This only considers the x/y component of the velocity gradient.
      ! TO DO:
      ! To get correct values on a sphere, the velocity gradient needs to
      ! be rotated such that the local tangential plane is the x/y plane.
      do iCell = 1, nCells
         kmax = maxLevelCell(iCell)
         do k = 1, kmax
            sn = dvel(1, 1, k, iCell) - dvel(2, 2, k, iCell)
            ss = dvel(1, 2, k, iCell) + dvel(2, 1, k, iCell)
            S(k, iCell) = sn*sn + ss*ss
            om(k, iCell) = dvel(1, 2, k, iCell) - dvel(2, 1, k, iCell)
            OW(k, iCell) = S(k, iCell) - om(k, iCell)*om(k, iCell)
         end do
         do k = kmax+1, nVertLevels
            S(k, iCell) = 0
            om(k, iCell) = 0
            OW(k, iCell) = 0
         end do
      end do
      OW = OW / OW_norm;

      ! Compute Lambda_2 parameter
      ! Lam2 only considers the x/y components of the velocity gradient,
      ! analogously to the OW computation above.
      ! Lam2_R3 considers the full 3-dimensional velocity gradient.
      do iCell = 1, nCells
         kmax = maxLevelCell(iCell)
         do k = 1, kmax
            !sym =  0.5 * (dvel(1:2, 1:2, k, iCell) + Transpose(dvel(1:2, 1:2, k, iCell)))
            !asym = 0.5 * (dvel(1:2, 1:2, k, iCell) - Transpose(dvel(1:2, 1:2, k, iCell)))
            sym3 =  0.5 * (dvel(:, :, k, iCell) + Transpose(dvel(:, :, k, iCell)))
            asym3 = 0.5 * (dvel(:, :, k, iCell) - Transpose(dvel(:, :, k, iCell)))
            sym2 =  sym3(1:2, 1:2)
            asym2 = asym3(1:2, 1:2)
            T2 = matmul(sym2, sym2) + matmul(asym2, asym2)
            T3 = matmul(sym3, sym3) + matmul(asym3, asym3)

            ! Compute eigen-values of 2x2 matrix
            call compute_ev_2(T2, lambda, lambdaI)
            ! Take the second eigen-value, multiply by 4 to be on the same scale as OW
            Lam1(k, iCell) = 4*lambda(1)
            Lam2(k, iCell) = 4*lambda(2)

            ! Compute eigen-values of 3x3 matrix
            call compute_ev_3(T3, lambda3, lambdaI3)
            ! Take the second eigen-value, multiply by 4 to be on the same scale as OW
            Lam2_R3(k, iCell) = 4*lambda3(2)
         end do
         do k = kmax+1, nVertLevels
            Lam1(k, iCell) = 0
            Lam2(k, iCell) = 0
            Lam2_R3(k, iCell) = 0
         end do
      end do

      Lam2 = Lam2 / Lam2_norm;
      Lam2_R3 = Lam2_R3 / Lam2_norm;

   end subroutine ocn_compute_OW_values!}}}


!***********************************************************************
!
!  routine ocn_threshold_OW
!
!> \brief   Threshold OW values
!> \author  Andre Schmeisser
!> \date    August 2014
!
!-----------------------------------------------------------------------

   subroutine ocn_threshold_OW(nVertLevels, nCells, threshold, OW, OW_thresh)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: OW
      real (kind=RKIND), intent(in) :: threshold
      integer, intent(in) :: nVertLevels, nCells

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:,:), intent(out) :: OW_thresh

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer k, iCell

      do k = 1, nVertLevels
         do iCell = 1, nCells
            if (OW(k, iCell) < threshold) then
               OW_thresh(k,iCell) = 1
            else
               OW_thresh(k,iCell) = 0
            end if
         end do
      end do

   end subroutine ocn_threshold_OW!}}}



!***********************************************************************
!
!  function find
!
!> \brief   The "find" part of a Union-Find data-structure
!> \author  Andre Schmeisser
!> \date    August 2014
!
!-----------------------------------------------------------------------

   recursive integer function find(unionMap, x) result(res)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: x

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), intent(inout) :: unionMap

      integer xEquiv

      xEquiv = unionMap(x)
      if (xEquiv < 0 .or. xEquiv == x) then
         res = x
         return
      end if

      ! follow chain to root and do path compression
      unionMap(x) = find(unionMap, xEquiv)
      res = unionMap(x)

   end function find!}}}

!***********************************************************************
!
!  routine union
!
!> \brief   The "union" part of a Union-Find data-structure
!> \author  Andre Schmeisser
!> \date    August 2014
!
!-----------------------------------------------------------------------

   subroutine union(unionMap, x, y)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: x, y

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), intent(inout) :: unionMap

      integer xRoot, yRoot

      ! find canonical elements for sets x and y
      xRoot = find(unionMap, x)
      yRoot = find(unionMap, y)

      if (xRoot == yRoot) then
         return
      end if

      ! join sets
      if (xRoot < yRoot) then ! this check enforces a monotony in the indices
         unionMap(yRoot) = xRoot
      else
         unionMap(xRoot) = yRoot
      end if

   end subroutine union!}}}

!***********************************************************************
!
!  routine ocn_compute_OW_component_IDs
!
!> \brief   Compute connected components IDs of the thresholded OW field
!> \author  Andre Schmeisser
!> \date    August 2014
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_OW_component_IDs(dminfo, block, meshPool, processorId, nVertLevels, &
                                           nCells, OW_thresh, OW_cc_id)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (dm_info), intent(in) :: dminfo
      type (block_type), intent(in) :: block
      integer, intent(in) :: nVertLevels, nCells, processorId
      type (mpas_pool_type), intent(in) :: meshPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:,:), intent(inout) :: OW_thresh

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:,:), intent(out) :: OW_cc_id

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: amOkuboWeissPool

      integer, dimension(:), pointer :: nEdgesOnCell, indexToCellID
      integer, dimension(:,:), pointer :: cellsOnCell

      type (field2DInteger), pointer :: OW_cc_id_field

      integer, pointer :: nCellsSolve
      integer k, iCell, i, iOtherCell, cellIdx, otherCellIdx
      integer origK, origCell, nLocalCCs
      integer, dimension(nVertLevels*nCells) :: unionMap
      integer communcationLoop
      integer, dimension(1) :: updateArray, globalUpdateArray

      call mpas_timer_start("CC local", .false., CCLocalTimer)

      call mpas_pool_get_subpool(block % structs, 'amOkuboWeiss', amOkuboWeissPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_field(amOkuboWeissPool, 'eddyID', OW_cc_id_field)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      unionMap(:) = -1

      ! Compute local connected components for this domain
      !
      ! TODO: To be completely correct, this should go to maxVertLevels only
      ! and be in the inner loop. But then, the check for adjacent cells needs
      ! to check if these are valid for that vert level, too.
      ! It is ok to do it like this because for land cells the OW is zero anyway,
      ! thus OW_thresh is not set.
      do k = 1, nVertLevels
         do iCell = 1, nCells
            cellIdx = iCell + (k-1)*nCells

            if (OW_thresh(k, iCell) > 0) then
               unionMap(cellIdx) = cellIdx
               ! Check all neighboring cells, if any of them is "set", mark
               ! them as belonging to the same component

               ! Check adjacent cells on same height level
               do i = 1, nEdgesOnCell(iCell)
                  iOtherCell = cellsOnCell(i, iCell)
                  if (iOtherCell > 0 .and. iOtherCell <= nCells) then
                     if (OW_thresh(k, iOtherCell) > 0) then
                        otherCellIdx = iOtherCell + (k-1)*nCells
                        call union(unionMap, cellIdx, otherCellIdx)
                     end if
                  end if
               end do

               ! Check cell with lower k-level (higher level has not yet been processed)
               if (k > 1) then
                  if (OW_thresh(k-1, iCell) > 0) then
                     otherCellIdx = cellIdx - nCells
                     call union(unionMap, cellIdx, otherCellIdx)
                  end if
               end if
            end if
         end do
      end do

      ! The number of local connected components is an upper bound for the
      ! number of globally connected components in this block.
      ! (Two separate local CCs might be globally connected through other blocks)
      nLocalCCs = ocn_count_local_connected_components(unionMap)

      do k = 1, nVertLevels
         do iCell = 1, nCells
            cellIdx = iCell + (k-1)*nCells
            if (unionMap(cellIdx) >= 0) then
               OW_cc_id(k, iCell) = find(unionMap, cellIdx)
               ! Convert local ID in the range [1, nVertLevels*nCells] to
               ! global ID in the range [1, nVertLevels*nCellsGlobal]
               origK = OW_cc_id(k, iCell) / nCells + 1
               origCell = mod(OW_cc_id(k, iCell)-1, nCells)+1
               OW_cc_id(k, iCell) = (nCellsGlobal*(origK-1))+indexToCellID(origCell)
            else
               OW_cc_id(k, iCell) = -1
            end if
         end do
      end do

      call mpas_timer_stop("CC local", CCLocalTimer)
      call mpas_timer_start("CC communcation loop", .false., CCComLoopTimer)

      ! Merge components that are connected over domain boundaries
      ! Communicate the CC IDs over the halo field.
      ! If a component is connected over domains, the ID assigned on this domain
      ! and the one from a neighboring domain differ. Use the lower ID as the
      ! canonical one and update the affected IDs appropriately.
      ! Iterate this process until no more changes occur to propagate the correct
      ! results.
      ! This is not the most efficient way to globally merge components, but
      ! uses the existing infrastructure to communicate with neighboring domains.
      ! By the nature of eddies being local, components should not stretch over
      ! many components and thus this quickly converges.
      do communcationLoop = 1, 20
         call mpas_dmpar_exch_halo_field(OW_cc_id_field)

         ! Loop over the ghost cells updated with data from adjacent blocks
         ! and merge
         updateArray(1) = 0
         do k = 1, nVertLevels
            do iCell = nCellsSolve+1, nCells
               ! Check if ghost cell is part of an eddy
               if (OW_cc_id(k, iCell) > 0) then
                  ! Loop over neighbors of ghost cell
                  do i = 1, nEdgesOnCell(iCell)
                     iOtherCell = cellsOnCell(i, iCell)
                     if (iOtherCell > 0 .and. iOtherCell <= nCellsSolve) then
                        ! Check if neighbor is part of an eddy
                        if (OW_cc_id(k, iOtherCell) > 0) then
                           if (OW_cc_id(k, iCell) < OW_cc_id(k, iOtherCell)) then
                              ! Ghost cell has lower ID, it takes precedence
                              ! over this block.
                              updateArray(1) = 1
                              call ocn_update_components(OW_cc_id, k, iCell, iOtherCell, nCells, nVertLevels)
                           else if (OW_cc_id(k, iCell) > OW_cc_id(k, iOtherCell)) then
                              ! Ghost cell has higher ID, this block takes precedende.
                              OW_cc_id(k, iCell) = OW_cc_id(k, iOtherCell)
                              updateArray(1) = 1
                           end if
                        end if
                     end if
                  end do
               end if
            end do
         end do

         write (*,*) "Something updated in communication process:", communcationLoop, updateArray(1)
         if (communcationLoop > 2) then
            call mpas_dmpar_max_int_array(dminfo, 1, updateArray, globalUpdateArray)
            if (globalUpdateArray(1) == 0) then
               write (*,*) "Nothing changed globally, exiting after round", communcationLoop
               exit
            end if
         end if
      end do
      call mpas_timer_stop("CC communcation loop", CCComLoopTimer)

      call mpas_timer_start("CC eddy stats", .false., CCStatsTimer)
      call ocn_compute_eddy_stats(dminfo, block, nVertLevels, nCells, nCellsSolve, nLocalCCs, &
                                  nEdgesOnCell, cellsOnCell, OW_cc_id, OW_thresh)
      call mpas_timer_stop("CC local", CCStatsTimer)

   end subroutine ocn_compute_OW_component_IDs!}}}

!***********************************************************************
!
!  function ocn_count_local_connected_components
!
!> \brief   Count the number of local connected components
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!> Local connected components are defined by being mapped to the same ID.
!> For each component, there is one entry in the Union-Find data structure
!> that is representative of the union, by mapping to itself. Thus this
!> function simply counts the entries "unionMap(i) == i" that map to
!> themselves.
!
!-----------------------------------------------------------------------

   integer function ocn_count_local_connected_components(unionMap) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), intent(in) :: unionMap

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer i, count

      count = 0
      do i = 1, size(unionMap)
         if (unionMap(i) == i) then
            count = count+1
         end if
      end do

      ocn_count_local_connected_components = count

   end function ocn_count_local_connected_components!}}}


!***********************************************************************
!
!  routine ocn_update_components
!
!> \brief   Update CC ID for all cells with a given ID to that of another ID
!> \author  Andre Schmeisser
!> \date    August 2014
!
!-----------------------------------------------------------------------

   subroutine ocn_update_components(OW_cc_id, srcK, iCell, iOtherCell, nCells, nVertLevels)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: srcK, iCell, iOtherCell, nVertLevels, nCells

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:,:), intent(inout) :: OW_cc_id

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer k, i
      integer fromId, replaceId

      fromId = OW_cc_id(srcK, iOtherCell)
      replaceId = OW_cc_id(srcK, iCell)

      do k = 1, nVertLevels
         do i = 1, nCells
            if (OW_cc_id(k, i) == fromId) then
               OW_cc_id(k, i) = replaceId
            end if
         end do
      end do

   end subroutine ocn_update_components!}}}

!***********************************************************************
!
!  routine ocn_compute_eddy_stats
!
!> \brief   Aggregate statistics of partial eddies
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!> Each (partial) eddy has an ID used for sorting and several statistics
!> associated with it. For this, statistics are stored in an "array of structs"
!> fashion, where the first value represents the ID.
!> QSort is used to sort these blocks of size elemSize.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_eddy_stats(dminfo, block, nVertLevels, nCells, nCellsSolve, &
                                     nLocalCCs, nEdgesOnCell, cellsOnCell, OW_cc_id, OW_thresh)!{{{


      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (dm_info), intent(in) :: dminfo
      type (block_type), intent(in) :: block
      integer, intent(in) :: nVertLevels, nCells, nCellsSolve, nLocalCCs
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: cellsOnCell
      integer, dimension(:,:), intent(in) :: OW_thresh, OW_cc_id

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:,:), pointer ::  layerThickness
      real (kind=RKIND), dimension(:), pointer :: posX, posY, posZ
      real (kind=RKIND), dimension(:,:), pointer :: velX, velY, velZ
      logical, pointer :: on_a_sphere, use_lat_lon_coords

      integer, dimension(nVertLevels*nCells) :: nextCellIdxK, nextCellIdxI
      real (kind=RKIND), dimension(nLocalCCs) :: sumVol, numCells, origCCId, wsVelX, wsVelY, wsVelZ, &
            wsPosX, wsPosY, wsPosZ
      real (kind=RKIND) :: vol
      integer k, iCell, i, curCC, iIdx, curI, curK, iOtherCell
      integer, dimension(1) :: numCCsArray, maxNumCCsArray

      integer, dimension(size(OW_thresh,1), size(OW_thresh,2)) :: OW_thresh_local

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(ocnConfigs, 'config_okubo_weiss_use_lat_lon_coords', use_lat_lon_coords)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      if (.not. on_a_sphere) then
         use_lat_lon_coords = .false.
      end if

      if (use_lat_lon_coords) then
         call mpas_pool_get_array(meshPool, 'lonCell', posX)
         call mpas_pool_get_array(meshPool, 'latCell',  posY)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velX)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velY)
      else
         call mpas_pool_get_array(meshPool, 'xCell', posX)
         call mpas_pool_get_array(meshPool, 'yCell', posY)
         call mpas_pool_get_array(meshPool, 'zCell', posZ)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityX', velX)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityY', velY)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZ', velZ)
      end if

      call mpas_timer_start("stats per proc", .false., CCStatsPerProcT)

      OW_thresh_local = OW_thresh

      ! Compute stats on each connected component of an eddy
      !
      ! Loop through cells, marking cells as processed on the go
      ! by resetting the OW_thresh_local variable
      ! If a cell belongs to a component (and has not been marked as
      ! processed), give it a new eddy ID and start doing statistics on
      ! all connected cells. This is done by doing a "flood fill" of
      ! the connected cells, marking each neighboring cell as processed
      ! first (so that it's processed only once) and then putting it into
      ! a list of cells to be processed (nextCellIdxI, nextCellIdxK).
      !
      ! Ghost cells are not processed in the statistics so that they are not
      ! counted multiple times. They are processed for finding all the neighbors
      ! in the flood fill process, though. Otherwise, more local eddies parts
      ! might be created than there is memory alloated for statistics arrays,
      ! as the counting of connected components further up also includes these
      ! cells.
      curCC = 0
      do k = 1, nVertLevels
         do iCell = 1, nCellsSolve
            if (OW_thresh_local(k, iCell) > 0) then
               ! Found a new eddy.
               curCC = curCC + 1
               if (curCC > nLocalCCs) then
                  write (stderrUnit,*) "ERROR: curCC > nLocalCCs", curCC, nLocalCCs
                  write (stderrUnit,*) "THIS IS A BUG!"
               end if
               ! Put this cell into the list of cells to be processed as
               ! seed for the flood fill of this eddy
               OW_thresh_local(k, iCell) = 0
               iIdx = 1
               nextCellIdxK(1) = k
               nextCellIdxI(1) = iCell
               sumVol(curCC) = 0
               numCells(curCC) = 0
               wsPosX(curCC) = 0
               wsPosY(curCC) = 0
               wsPosZ(curCC) = 0
               wsVelX(curCC) = 0
               wsVelY(curCC) = 0
               wsVelZ(curCC) = 0
               origCCId(curCC) = OW_cc_id(k, iCell)
               ! Flood fill: loop over list of cells to be processed
               ! iIdx is the position in the list, while curI/curK are the
               ! actual indices for the cell being processed
               do while (iIdx > 0)
                  curI = nextCellIdxI(iIdx)
                  curK = nextCellIdxK(iIdx)
                  iIdx = iIdx - 1
                  if (curI <= nCellsSolve) then
                     ! Compute stats on cell
                     vol = areaCell(curI)*layerThickness(curK, curI)
                     sumVol(curCC) = sumVol(curCC) + vol
                     numCells(curCC) = numCells(curCC) + 1
                     wsPosX(curCC) = wsPosX(curCC) + vol * posX(curI)
                     wsPosY(curCC) = wsPosY(curCC) + vol * posY(curI)
                     wsVelX(curCC) = wsVelX(curCC) + vol * velX(curK, curI)
                     wsVelY(curCC) = wsVelY(curCC) + vol * velY(curK, curI)
                     if (use_lat_lon_coords) then
                        wsPosZ(curCC) = wsPosZ(curCC) + vol * curK
                     else
                        wsPosZ(curCC) = wsPosZ(curCC) + vol * posZ(curI)
                        wsVelZ(curCC) = wsVelZ(curCC) + vol * velZ(curK, curI)
                     end if
                  end if

                  ! Check neighbors of cell and put them into list:
                  ! Adjacent cells on same height level
                  do i = 1, nEdgesOnCell(curI)
                     iOtherCell = cellsOnCell(i, curI)
                     if (iOtherCell > 0 .and. iOtherCell <= nCells) then
                        if (OW_thresh_local(curK, iOtherCell) > 0) then
                           iIdx = iIdx + 1
                           nextCellIdxK(iIdx) = curK
                           nextCellIdxI(iIdx) = iOtherCell
                           OW_thresh_local(curK, iOtherCell) = 0
                        end if
                     end if
                  end do
                  ! Cell on lower k-level
                  if (curK > 1) then
                     if (OW_thresh_local(curK-1, curI) > 0) then
                        iIdx = iIdx + 1
                        nextCellIdxK(iIdx) = curK-1
                        nextCellIdxI(iIdx) = curI
                        OW_thresh_local(curK-1, curI) = 0
                     end if
                  end if
                  ! Cell on higher k-level
                  if (curK < nVertLevels) then
                     if (OW_thresh_local(curK+1, curI) > 0) then
                        iIdx = iIdx + 1
                        nextCellIdxK(iIdx) = curK+1
                        nextCellIdxI(iIdx) = curI
                        OW_thresh_local(curK+1, curI) = 0
                     end if
                  end if
               end do
            end if
         end do
      end do

      call mpas_timer_stop("stats per proc", CCStatsPerProcT)

      numCCsArray(1) = curCC
      call mpas_dmpar_max_int_array(dminfo, 1, numCCsArray, maxNumCCsArray)
      write (*,*) "Number of (partial) eddies on this domain:", numCCsArray(1)
      write (*,*) "Global maximum of (partial) eddies on a domain:", maxNumCCsArray(1)

      call mpas_timer_start("aggregate stats", .false., CCAggregateStatsT)
      call ocn_aggregate_eddy_stats(dmInfo, use_lat_lon_coords, curCC, maxNumCCsArray(1), &
                                    origCCId, sumVol, numCells, wsPosX, wsPosY, wsPosZ, &
                                    wsVelX, wsVelY, wsVelZ)
      call mpas_timer_stop("aggregate stats", CCAggregateStatsT)

   end subroutine ocn_compute_eddy_stats!}}}


!***********************************************************************
!
!  routine ocn_aggregate_eddy_stats
!
!> \brief   Aggregate statistics of partial eddies
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!> Send all partial eddy statistics to IO node and merge them into
!> statistics of complete eddies.
!
!-----------------------------------------------------------------------

   subroutine ocn_aggregate_eddy_stats(dmInfo, use_lat_lon_coords,  numCCs, maxCCsPerDomain, &
                                       origCCId, sumVol, numCells, wsPosX, wsPosY, wsPosZ, &
                                       wsVelX, wsVelY, wsVelZ)!{{{


      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (dm_info), intent(in) :: dminfo
      logical, intent(in) :: use_lat_lon_coords
      integer, intent(in) :: numCCs, maxCCsPerDomain
      real (kind=RKIND), dimension(:), intent(in) :: origCCId, sumVol, numCells, &
                             wsPosX, wsPosY, wsPosZ, wsVelX, wsVelY, wsVelZ

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(9 * maxCCsPerDomain) :: stats
      real (kind=RKIND), dimension(9 * maxCCsPerDomain * dmInfo % nprocs) :: &
         globalStats, aggregated

      integer i, j, ierr, curId, listIdx, totalCCs

      call mpas_timer_start("gather", .false., CCAggGatherT)
      ! Multiplex all eddy statistics into one array so it can be send with
      ! one MPI communication call
      stats = 0
      do i = 1, numCCs
         stats(9*(i-1) + 1) = origCCId(i)
         stats(9*(i-1) + 2) = sumVol(i)
         stats(9*(i-1) + 3) = wsPosX(i)
         stats(9*(i-1) + 4) = wsPosY(i)
         stats(9*(i-1) + 5) = wsPosZ(i)
         stats(9*(i-1) + 6) = wsVelX(i)
         stats(9*(i-1) + 7) = wsVelY(i)
         stats(9*(i-1) + 8) = wsVelZ(i)
         stats(9*(i-1) + 9) = numCells(i)
      end do

      ! Send all local statistics to IO node
      call MPI_Gather(stats, 9*maxCCsPerDomain, MPI_REALKIND, &
                      globalStats, 9*maxCCsPerDomain, MPI_REALKIND, &
                      IO_NODE, dmInfo % comm, ierr)
      call mpas_timer_stop("gather", CCAggGatherT)

      if (dminfo % my_proc_id == IO_NODE) then
         call mpas_timer_start("sort/reduce", .false., CCAggSortT)

         ! Sort partial eddy stats according to eddy ID for efficient merging
         totalCCs = maxCCsPerDomain * dmInfo % nprocs
         if (totalCCs > 0) then
            call ocn_sort_eddy_stats(globalStats, totalCCs, 9)
         end if

         ! Aggregate the results
         listIdx = 0
         aggregated = 0
         curId = -1
         do i = 1, totalCCs
            if (globalStats(9*(i-1)+1) /= curId) then
               curId = globalStats(9*(i-1)+1)
               if (curId == 0) then
                  exit
               end if
               listIdx = listIdx + 1
            end if
            aggregated(9*(listIdx-1)+1) = curId
            ! Merge stats about weighted position, velocity, etc. by summing
            do j = 2, 9
               aggregated(9*(listIdx-1)+j) = aggregated(9*(listIdx-1)+j) &
                                           + globalStats(9*(i-1)+j)
            end do
         end do
         write (*,*) "Number of aggregated eddies:", listIdx
         call mpas_timer_stop("sort/reduce", CCAggSortT)

         ! Output aggregated
         call mpas_timer_start("output", .false., CCOutputT)
         call ocn_output_eddy_stats(listIdx, aggregated, use_lat_lon_coords)
         call mpas_timer_stop("output", CCOutputT)
      end if

   end subroutine ocn_aggregate_eddy_stats!}}}


!***********************************************************************
!
!  routine ocn_sort_eddy_stats
!
!> \brief   Sorts eddy statistics according to eddy ID.
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!> Each (partial) eddy has an ID used for sorting and several statistics
!> associated with it. For this, statistics are stored in an "array of structs"
!> fashion (blocks of size elemSize), where the first value represents the ID.
!> QSort is used to sort these blocks.
!> Eddies split over domains will be represented by several partial eddies
!> with the same ID, which after sorting will be adjacent in memory,
!> allowing for efficient aggregation.
!
!-----------------------------------------------------------------------

   subroutine ocn_sort_eddy_stats(stats, numRecords, recordLength)
      use iso_c_binding
      use Sort

      integer, intent(in) :: numRecords, recordLength
      real (kind=RKIND), dimension(numRecords*recordLength), intent(inout), target :: stats
      integer(c_size_t) elemCount, elemSize

      elemCount = numRecords
      elemSize = recordLength * SIZE_REAL

      call qsort(c_loc(stats(1)), elemCount, elemSize, c_funloc(compareRealDescending))

   end subroutine ocn_sort_eddy_stats


!***********************************************************************
!
!  routine ocn_output_eddy_stats
!
!> \brief   Outputs eddy statistics
!> \author  Andre Schmeisser
!> \date    August 2014
!
!-----------------------------------------------------------------------

   subroutine ocn_output_eddy_stats(numEddies, aggData, use_lat_lon_coords)!{{{


      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: numEddies
      real (kind=RKIND), dimension(:), intent(in) :: aggData
      logical, intent(in) :: use_lat_lon_coords

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND) :: v
      integer i, fileID

      fileID = mpas_get_free_unit()

      ! Create output file.
      ! If this is not the first output step, open for appending.
      if (.not. statsFileCreated) then
          open(fileID, file='okubo_weiss_eddy_stats.txt', STATUS='UNKNOWN', POSITION='rewind')
          statsFileCreated = .true.
      else
          open(fileID, file='okubo_weiss_eddy_stats.txt', STATUS='UNKNOWN', POSITION='append')
          write(fileID, *) ""
      endif 

      ! Output number of eddies and statistics for each eddy
      write (fileID, *) "Number of eddies:", numEddies
      do i = 1, numEddies
         write (fileID, *) "Eddy ID:", i, int(aggData(9*(i-1)+1))
         v = aggData(9*(i-1)+2)
         write (fileID, *) "Volume sum:", v
         write (fileID, *) "Number of cells:", aggData(9*(i-1)+9)
         write (fileID, *) "Weighted position:", aggData(9*(i-1)+3), aggData(9*(i-1)+4), aggData(9*(i-1)+5)
         write (fileID, *) "Average position:", aggData(9*(i-1)+3)/v, aggData(9*(i-1)+4)/v, aggData(9*(i-1)+5)/v
         if (use_lat_lon_coords) then
            write (fileID, *) "Average velocity:", aggData(9*(i-1)+6)/v, aggData(9*(i-1)+7)/v
         else
            write (fileID, *) "Average Velocity:", aggData(9*(i-1)+6)/v, aggData(9*(i-1)+7)/v, aggData(9*(i-1)+8)/v
         end if
      end do

      close(fileID)

   end subroutine ocn_output_eddy_stats!}}}

!***********************************************************************
!
!  function mpas_get_free_unit
!
!-----------------------------------------------------------------------

   integer function mpas_get_free_unit()!{{{
      implicit none

      integer :: index
      logical :: isOpened

      mpas_get_free_unit = 0
      do index = 1,99
         if((index /= 5) .and. (index /= 6)) then
            inquire(unit = index, opened = isOpened)
            if( .not. isOpened) then
               mpas_get_free_unit = index
               return
            end if
         end if
      end do
   end function mpas_get_free_unit!}}}

!***********************************************************************
!
!  routine mpas_velocity_gradient_R3Cell
!
!> \brief   Computes the velocity gradient at cell centers, in R3
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  This routine computes the velocity gradient at cell centers using the weak 
!>  derivative.  Output is an R3 velocity gradient tensor in 3x3 format.
!
!-----------------------------------------------------------------------

   subroutine mpas_velocity_gradient_R3Cell(normalVelocity, tangentialVelocity, &
      meshPool, edgeSignOnCell, edgeTangentVectors, includeHalo, &
      nVertLevels, nEdges, &
      velocityGradient)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         edgeTangentVectors,   &!< Input: unit vector tangent to an edge
         normalVelocity,      &!< Input: Horizontal velocity normal to edge
         tangentialVelocity    !< Input: Horizontal velocity tangent to edge

      integer, dimension(:,:), intent(in) :: &
         edgeSignOnCell        !< Input: Direction of vector connecting cells

      integer, intent(in) :: &
         nVertLevels, nEdges

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      logical, intent(in) :: & 
         includeHalo !< Input: If true, halo cells and edges are included in computation

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: &
         velocityGradient   !< Output: strain rate tensor at cell center, R3, in symmetric 6-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, iCell, i, j, k
      integer, pointer :: nCellsCompute

      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: edgesOnCell


      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND), dimension(3,3) :: outerProductEdge3x3
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      real (kind=RKIND), dimension(3,3, nVertLevels,nEdges) :: outerProductEdgeFull

      if (includeHalo) then
         call mpas_pool_get_dimension(meshPool, 'nCells', nCellsCompute)
      else 
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsCompute)
      end if

      !call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      do iEdge=1,nEdges
         do k=1,nVertLevels
           do i=1,3
             do j=1,3
               ! outer produce at each edge:
               ! u_e n_e n_e* + v_e n_e \tilde{n}_e* 
               outerProductEdge3x3(i,j) = edgeNormalVectors(i,iEdge) &
                       *(  normalVelocity(k,iEdge)    *edgeNormalVectors(j,iEdge) &
                         + tangentialVelocity(k,iEdge)*edgeTangentVectors(j,iEdge) &
                           )
             enddo
           enddo
           outerProductEdgeFull(:,:,k,iEdge) = outerProductEdge3x3(:,:)
         enddo
      enddo

      velocityGradient = 0.0
      do iCell = 1, nCellsCompute
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               ! edgeSignOnCell is to get outward unit normal on edgeNormalVectors
               ! minus sign in front is to match form on divergence operator
               velocityGradient(:,:,k,iCell) = velocityGradient(:,:,k,iCell) &
                 - edgeSignOnCell(i,iCell)*outerProductEdgeFull(:,:,k,iEdge)*invAreaCell*dvEdge(iEdge) 
            end do
         end do
      end do

   end subroutine mpas_velocity_gradient_R3Cell!}}}


!***********************************************************************
!
!  routine ocn_compute_okubo_weiss
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_okubo_weiss(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: amOkuboWeissPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: amOkuboWeiss

      integer, pointer :: nVertLevels, nCells, nEdges
      integer, dimension(:), pointer :: maxLevelCell    
      integer, dimension(:,:), pointer :: edgeSignOnCell

      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: tangentialVelocity
      real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors
      
      real (kind=RKIND), dimension(:,:), pointer :: om, OW 
      integer, dimension(:,:), pointer :: OW_cc_id

      type(field2DReal), pointer :: SField, Lam1Field, Lam2Field, Lam2_R3Field
      type(field4DReal), pointer :: velocityGradientField
      type(field2DInteger), pointer :: OW_threshField

      real (kind=RKIND), dimension(:,:), pointer :: S, Lam1, Lam2, Lam2_R3
      real (kind=RKIND), dimension(:,:,:,:), pointer :: velocityGradient
      integer, dimension(:,:), pointer :: OW_thresh

      logical, pointer :: config_compute_OW_CC
      real (kind=RKIND), pointer :: OW_normalization, Lam2_normalization, threshold

      call mpas_timer_start("compute_okubo_weiss", .false., amOkuboWeissTimer)
      err = 0
      dminfo = domain % dminfo

      call mpas_pool_get_config(ocnConfigs, 'config_okubo_weiss_normalization', OW_normalization)
      call mpas_pool_get_config(ocnConfigs, 'config_okubo_weiss_lambda2_normalization', Lam2_normalization)
      call mpas_pool_get_config(ocnConfigs, 'config_okubo_weiss_compute_connected_components', config_compute_OW_CC)
      call mpas_pool_get_config(ocnConfigs, 'config_okubo_weiss_threshold_value', threshold)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'okuboWeissScratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'amOkuboWeiss', amOkuboWeissPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
         call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)
         
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)

         call mpas_pool_get_array(amOkuboWeissPool, 'okuboWeiss', OW)
         call mpas_pool_get_array(amOkuboWeissPool, 'eddyID', OW_cc_id)
         call mpas_pool_get_array(amOkuboWeissPool, 'vorticity', om)

         call mpas_pool_get_field(scratchPool, 'velocityGradient', velocityGradientField)         
         call mpas_pool_get_field(scratchPool, 'thresholdedOkuboWeiss', OW_threshField)
         call mpas_pool_get_field(scratchPool, 'shearAndStrain', SField)
         call mpas_pool_get_field(scratchPool, 'lambda1', Lam1Field)
         call mpas_pool_get_field(scratchPool, 'lambda2', Lam2Field)
         call mpas_pool_get_field(scratchPool, 'lambda2R3', Lam2_R3Field)

         call mpas_allocate_scratch_field(velocityGradientField, .true.)
         call mpas_allocate_scratch_field(OW_threshField, .true.)
         call mpas_allocate_scratch_field(SField, .true.)
         call mpas_allocate_scratch_field(Lam1Field, .true.)
         call mpas_allocate_scratch_field(Lam2Field, .true.)
         call mpas_allocate_scratch_field(Lam2_R3Field, .true.)

         velocityGradient => velocityGradientField % array
         OW_thresh => OW_threshField % array
         S => SField % array
         Lam1 => Lam1Field % array
         Lam2 => Lam2Field % array
         Lam2_R3 => Lam2_R3Field % array

         ! Compute velocity gradient         
         call mpas_velocity_gradient_R3Cell(normalVelocity, tangentialVelocity, &
               meshPool, edgeSignOnCell, edgeTangentVectors, .true., &
               nVertLevels, nEdges, velocityGradient)

         ! Compute Okubo-Weiss and Lambda 2 values
         call ocn_compute_OW_values(velocityGradient, nVertLevels, nCells, maxLevelCell, OW, OW_normalization, &
               Lam2, Lam2_R3, Lam2_normalization, S, om, Lam1)

         ! Threshold field OW
         call ocn_threshold_OW(nVertLevels, nCells, threshold, OW, OW_thresh)

         ! Compute connected components of thresholded field
         if (config_compute_OW_CC) then
            call mpas_timer_start("OW connected components", .false., amOWCCTimer)
            call ocn_compute_OW_component_IDs(dminfo, block, meshPool, dminfo % my_proc_id, nVertLevels, &
                                              nCells, OW_thresh, OW_cc_id)
            call mpas_timer_stop("OW connected components", amOWCCTimer)
         end if

         call mpas_deallocate_scratch_field(velocityGradientField, .true.)
         call mpas_deallocate_scratch_field(OW_threshField, .true.)
         call mpas_deallocate_scratch_field(SField, .true.)
         call mpas_deallocate_scratch_field(Lam1Field, .true.)
         call mpas_deallocate_scratch_field(Lam2Field, .true.)
         call mpas_deallocate_scratch_field(Lam2_R3Field, .true.)

         block => block % next
      end do

      call mpas_timer_stop("compute_okubo_weiss", amOkuboWeissTimer)

   end subroutine ocn_compute_okubo_weiss!}}}

!***********************************************************************
!
!  routine ocn_restart_okubo_weiss
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_okubo_weiss(domain, err)!{{{

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

   end subroutine ocn_restart_okubo_weiss!}}}

!***********************************************************************
!
!  routine ocn_write_okubo_weiss
!
!> \brief   MPAS-Ocean analysis output
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  This routine writes all output for this MPAS-Ocean analysis member.
!>  At this time this is just a stub, and all analysis output is written
!>  to the output file specified by config_output_name.
!
!-----------------------------------------------------------------------

   subroutine ocn_write_okubo_weiss(domain, err)!{{{

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

   end subroutine ocn_write_okubo_weiss!}}}

!***********************************************************************
!
!  routine ocn_finalize_okubo_weiss
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_okubo_weiss(domain, err)!{{{

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

   end subroutine ocn_finalize_okubo_weiss!}}}

end module ocn_okubo_weiss

! vim: foldmethod=marker
