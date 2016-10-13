! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_threading2


    contains


    !-----------------------------------------------------------------------
    !  routine mpas_threading_init
    !
    !> \brief  Pre-computes thread loop bounds for cell, edge, and vertex elements
    !> \author Michael Duda
    !> \date   6 July 2015
    !> \details
    !>  This routine is responsible for computing thread loop bounds for cell,
    !>  edge, and vertex elements in each block of the input blocklist argument.
    !>  Starting and ending loop bounds are computed for these three element
    !>  types for all elements (e.g., nCells) as well as owned elements (e.g.,
    !>  nCellsSolve).
    !>
    !>  When MPAS is compiled without OpenMP support, this routine computes loop
    !>  bounds as though there is just a single thread; otherwise, it is assumed
    !>  that all threads (given by OMP_get_num_threads()) will be used to
    !>  decompose each of the element ranges.
    !>
    !>  At present, a return value of 0 is always returned in the optional
    !>  output argument, ierr.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_threading_init(blocklist, ierr) 

        use mpas_derived_types, only : block_type
        use mpas_pool_routines, only : mpas_pool_get_dimension
#ifdef _OPENMP
        use omp_lib
#endif

        implicit none

        type (block_type), pointer :: blocklist
        integer, intent(out), optional :: ierr

        type (block_type), pointer :: block
        integer :: threadid
        integer, pointer :: nCells, nCellsSolve, nEdges, nEdgesSolve, nVertices, nVerticesSolve

        
        if (present(ierr)) ierr = 0

        block => blocklist
        do while (associated(block))
#ifdef _OPENMP
!$OMP PARALLEL
!$OMP MASTER
            block % nThreads = OMP_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
#else 
            block % nThreads = 1
#endif 

            allocate(block % cellThreadStart(block % nThreads))
            allocate(block % cellThreadEnd(block % nThreads))
            allocate(block % cellSolveThreadStart(block % nThreads))
            allocate(block % cellSolveThreadEnd(block % nThreads))
            allocate(block % edgeThreadStart(block % nThreads))
            allocate(block % edgeThreadEnd(block % nThreads))
            allocate(block % edgeSolveThreadStart(block % nThreads))
            allocate(block % edgeSolveThreadEnd(block % nThreads))
            allocate(block % vertexThreadStart(block % nThreads))
            allocate(block % vertexThreadEnd(block % nThreads))
            allocate(block % vertexSolveThreadStart(block % nThreads))
            allocate(block % vertexSolveThreadEnd(block % nThreads))

            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
            call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
            call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

#ifdef _OPENMP
!$OMP PARALLEL PRIVATE(threadid)
            threadid = OMP_get_thread_num()

            block % cellThreadStart(threadid+1) = (threadid * nCells / block % nThreads) + 1
            block % cellThreadEnd(threadid+1)   = ((threadid+1) * nCells / block % nThreads)
            block % cellSolveThreadStart(threadid+1) = (threadid * nCellsSolve / block % nThreads) + 1
            block % cellSolveThreadEnd(threadid+1)   = ((threadid+1) * nCellsSolve / block % nThreads)
            block % edgeThreadStart(threadid+1) = (threadid * nEdges / block % nThreads) + 1
            block % edgeThreadEnd(threadid+1)   = ((threadid+1) * nEdges / block % nThreads)
            block % edgeSolveThreadStart(threadid+1) = (threadid * nEdgesSolve / block % nThreads) + 1
            block % edgeSolveThreadEnd(threadid+1)   = ((threadid+1) * nEdgesSolve / block % nThreads)
            block % vertexThreadStart(threadid+1) = (threadid * nVertices / block % nThreads) + 1
            block % vertexThreadEnd(threadid+1)   = ((threadid+1) * nVertices / block % nThreads)
            block % vertexSolveThreadStart(threadid+1) = (threadid * nVerticesSolve / block % nThreads) + 1
            block % vertexSolveThreadEnd(threadid+1)   = ((threadid+1) * nVerticesSolve / block % nThreads)
!$OMP END PARALLEL
#else 
            block % cellThreadStart(1) = 1
            block % cellThreadEnd(1)   = nCells
            block % cellSolveThreadStart(1) = 1
            block % cellSolveThreadEnd(1)   = nCellsSolve
            block % edgeThreadStart(1) = 1
            block % edgeThreadEnd(1)   = nEdges
            block % edgeSolveThreadStart(1) = 1
            block % edgeSolveThreadEnd(1)   = nEdgesSolve
            block % vertexThreadStart(1) = 1
            block % vertexThreadEnd(1)   = nVertices
            block % vertexSolveThreadStart(1) = 1
            block % vertexSolveThreadEnd(1)   = nVerticesSolve
#endif 

            block => block % next
        end do

    end subroutine mpas_threading_init


    !-----------------------------------------------------------------------
    !  routine mpas_threading_finalize
    !
    !> \brief  Deallocates memory associated with threading in MPAS
    !> \author Michael Duda
    !> \date   6 July 2015
    !> \details
    !>  This routine deallocates any memory that was allocated in the call to
    !>  mpas_threading_init() for each block in the input block list.
    !>
    !>  At present, a return value of 0 is always returned in the optional
    !>  output argument, ierr.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_threading_finalize(blocklist, ierr) 

        use mpas_derived_types, only : block_type

        implicit none

        type (block_type), pointer :: blocklist
        integer, intent(out), optional :: ierr

        type (block_type), pointer :: block


        if (present(ierr)) ierr = 0

        block => blocklist
        do while (associated(block))
            deallocate(block % cellThreadStart)
            deallocate(block % cellThreadEnd)
            deallocate(block % cellSolveThreadStart)
            deallocate(block % cellSolveThreadEnd)
            deallocate(block % edgeThreadStart)
            deallocate(block % edgeThreadEnd)
            deallocate(block % edgeSolveThreadStart)
            deallocate(block % edgeSolveThreadEnd)
            deallocate(block % vertexThreadStart)
            deallocate(block % vertexThreadEnd)
            deallocate(block % vertexSolveThreadStart)
            deallocate(block % vertexSolveThreadEnd)

            block => block % next
        end do

    end subroutine mpas_threading_finalize
 
end module mpas_threading2
