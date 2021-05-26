! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!#define HALO_EXCH_DEBUG

module test_core_halo_exch

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_field_routines
   use mpas_dmpar
   use mpas_threading
   use mpas_log
   use mpas_timer

   implicit none
   private

   public :: test_core_halo_exch_test

   contains

   !***********************************************************************
   !
   !  routine test_core_halo_exch_test
   !
   !> \brief   MPAS Test Core halo exchange test
   !> \author  Doug Jacobsen
   !> \date    10/21/2015
   !> \details 
   !>  This routine performs tests related to halo exchanges.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_halo_exch_test(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      integer :: threadNum
      integer :: iErr

      err = 0

      threadNum = mpas_threading_get_thread_num()

      call mpas_timer_start('halo exch tests')
      if ( threadNum == 0 ) then
         call mpas_log_write(' - Performing exchange group tests')
      end if
      call test_core_halo_exch_group_test(domain, threadErrs, iErr)
      call mpas_threading_barrier()
      if ( threadNum == 0 ) then
         call mpas_log_write('    -- Return code: $i', intArgs=(/iErr/))
         err = ior(err, iErr)
      end if

      if ( threadNum == 0 ) then
         call mpas_log_write(' - Performing single field halo exchange tests')
      end if
      call test_core_halo_exch_single_field_test(domain, threadErrs, iErr)
      call mpas_threading_barrier()
      if ( threadNum == 0 ) then
         call mpas_log_write('    -- Return code: $i', intArgs=(/iErr/))
         err = ior(err, iErr)
      end if

      if ( threadNum == 0 ) then
         call mpas_log_write(' - Performing old halo exchange tests')
      end if
      call test_core_halo_exch_full_test(domain, threadErrs, iErr)
      call mpas_threading_barrier()
      if ( threadNum == 0 ) then
         call mpas_log_write('    -- Return code: $i', intArgs=(/iErr/))
         err = ior(err, iErr)
      end if

      call mpas_timer_stop('halo exch tests')

   end subroutine test_core_halo_exch_test!}}}

   !***********************************************************************
   !
   !  routine test_core_halo_exch_full_test
   !
   !> \brief   MPAS Test Core halo exchange full test
   !> \author  Doug Jacobsen
   !> \date    10/21/2015
   !> \details 
   !>  This routine performs tests related to halo exchanges. It creates
   !>  fields, and fills their zero halos with testable values. Then performs halo
   !>  exchanges, and differences the full result with the expect value. If the
   !>  test passes, it returns a zero, if it fails it returns 1.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_halo_exch_full_test(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool, haloExchTestPool

      type (field5DReal), pointer :: real5DField
      type (field4DReal), pointer :: real4DField
      type (field3DReal), pointer :: real3DField
      type (field2DReal), pointer :: real2DField
      type (field1DReal), pointer :: real1DField
      type (field3DInteger), pointer :: int3DField
      type (field2DInteger), pointer :: int2DField
      type (field1DInteger), pointer :: int1DField

      real (kind=RKIND), dimension(:, :, :, :, :), pointer :: real5D
      real (kind=RKIND), dimension(:, :, :, :), pointer :: real4D
      real (kind=RKIND), dimension(:, :, :), pointer :: real3D
      real (kind=RKIND), dimension(:, :), pointer :: real2D
      real (kind=RKIND), dimension(:), pointer :: real1D

      real (kind=RKIND) :: realValue
      integer :: integerValue

      integer, dimension(:, :, :), pointer :: int3D
      integer, dimension(:, :), pointer :: int2D
      integer, dimension(:), pointer :: int1D

      integer :: i, j, k, l, m
      integer :: iDim1, iDim2, iDim3, iDim4, iDim5
      integer, pointer :: nCells, nEdges, nVertices
      integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve
      integer, dimension(:), pointer :: indexToCellID
      integer, dimension(:), pointer :: indexToEdgeID
      integer, dimension(:), pointer :: indexToVertexID

      integer :: threadNum

      threadNum = mpas_threading_get_thread_num() + 1

      threadErrs(threadNum) = 0

      call test_core_halo_exch_setup_fields(domain, threadErrs, err)

      call mpas_threading_barrier()

      call mpas_timer_start('old halo exchanges')
      ! Perform halo exchanges
      call mpas_pool_get_subpool(domain % blocklist % structs, 'haloExchTest', haloExchTestPool)

      ! Exchange persistent cell fields
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistReal5D', real5DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistReal4D', real4DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistReal3D', real3DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistReal2D', real2DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistReal1D', real1DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistInt3D', int3DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistInt2D', int2DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellPersistInt1D', int1DField)

      call mpas_dmpar_exch_halo_field(real5DField)
      call mpas_dmpar_exch_halo_field(real4DField)
      call mpas_dmpar_exch_halo_field(real3DField)
      call mpas_dmpar_exch_halo_field(real2DField)
      call mpas_dmpar_exch_halo_field(real1DField)
      call mpas_dmpar_exch_halo_field(int3DField)
      call mpas_dmpar_exch_halo_field(int2DField)
      call mpas_dmpar_exch_halo_field(int1DField)

      ! Exchange persistent edge fields
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistReal5D', real5DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistReal4D', real4DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistReal3D', real3DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistReal2D', real2DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistReal1D', real1DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistInt3D', int3DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistInt2D', int2DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgePersistInt1D', int1DField)

      call mpas_dmpar_exch_halo_field(real5DField)
      call mpas_dmpar_exch_halo_field(real4DField)
      call mpas_dmpar_exch_halo_field(real3DField)
      call mpas_dmpar_exch_halo_field(real2DField)
      call mpas_dmpar_exch_halo_field(real1DField)
      call mpas_dmpar_exch_halo_field(int3DField)
      call mpas_dmpar_exch_halo_field(int2DField)
      call mpas_dmpar_exch_halo_field(int1DField)

      ! Exchange persistent vertex fields
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistReal5D', real5DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistReal4D', real4DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistReal3D', real3DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistReal2D', real2DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistReal1D', real1DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistInt3D', int3DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistInt2D', int2DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexPersistInt1D', int1DField)

      call mpas_dmpar_exch_halo_field(real5DField)
      call mpas_dmpar_exch_halo_field(real4DField)
      call mpas_dmpar_exch_halo_field(real3DField)
      call mpas_dmpar_exch_halo_field(real2DField)
      call mpas_dmpar_exch_halo_field(real1DField)
      call mpas_dmpar_exch_halo_field(int3DField)
      call mpas_dmpar_exch_halo_field(int2DField)
      call mpas_dmpar_exch_halo_field(int1DField)

      ! Exchange scratch cell fields
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal5D', real5DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal4D', real4DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal3D', real3DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal2D', real2DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal1D', real1DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt3D', int3DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt2D', int2DField)
      call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt1D', int1DField)

      call mpas_dmpar_exch_halo_field(real5DField)
      call mpas_dmpar_exch_halo_field(real4DField)
      call mpas_dmpar_exch_halo_field(real3DField)
      call mpas_dmpar_exch_halo_field(real2DField)
      call mpas_dmpar_exch_halo_field(real1DField)
      call mpas_dmpar_exch_halo_field(int3DField)
      call mpas_dmpar_exch_halo_field(int2DField)
      call mpas_dmpar_exch_halo_field(int1DField)

      ! Exchange edge scratch fields
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal5D', real5DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal4D', real4DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal3D', real3DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal2D', real2DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal1D', real1DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt3D', int3DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt2D', int2DField)
      call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt1D', int1DField)

      call mpas_dmpar_exch_halo_field(real5DField)
      call mpas_dmpar_exch_halo_field(real4DField)
      call mpas_dmpar_exch_halo_field(real3DField)
      call mpas_dmpar_exch_halo_field(real2DField)
      call mpas_dmpar_exch_halo_field(real1DField)
      call mpas_dmpar_exch_halo_field(int3DField)
      call mpas_dmpar_exch_halo_field(int2DField)
      call mpas_dmpar_exch_halo_field(int1DField)

      ! Exchange scratch vertex fields
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal5D', real5DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal4D', real4DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal3D', real3DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal2D', real2DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal1D', real1DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt3D', int3DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt2D', int2DField)
      call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt1D', int1DField)

      call mpas_dmpar_exch_halo_field(real5DField)
      call mpas_dmpar_exch_halo_field(real4DField)
      call mpas_dmpar_exch_halo_field(real3DField)
      call mpas_dmpar_exch_halo_field(real2DField)
      call mpas_dmpar_exch_halo_field(real1DField)
      call mpas_dmpar_exch_halo_field(int3DField)
      call mpas_dmpar_exch_halo_field(int2DField)
      call mpas_dmpar_exch_halo_field(int1DField)

      call mpas_timer_stop('old halo exchanges')

      call mpas_threading_barrier()

      call test_core_halo_exch_validate_fields(domain, threadErrs, err)

   end subroutine test_core_halo_exch_full_test!}}}


   !***********************************************************************
   !
   !  routine test_core_halo_exch_single_field_test
   !
   !> \brief   MPAS Test Core halo exchange single field test
   !> \author  Doug Jacobsen
   !> \date    01/12/2016
   !> \details 
   !>  This routine performs tests of single field halo exchanges. It creates
   !>  fields, and fills their zero halos with testable values. Then performs halo
   !>  exchanges, and differences the full result with the expect value. If the
   !>  test passes, it returns a zero, if it fails it returns 1.
   !>  To perform the halo exchanges, it uses the single field halo exchange
   !>  routines rather than exchange groups.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_halo_exch_single_field_test(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      integer :: threadNum

      threadNum = mpas_threading_get_thread_num() + 1

      threadErrs(threadNum) = 0

      call test_core_halo_exch_setup_fields(domain, threadErrs, err)

      call mpas_threading_barrier()

      call mpas_timer_start('single halo exchanges')
      ! Perform halo exchanges

      ! Exchange persistent cell fields
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistReal5D')
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistReal4D')
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistReal3D')
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistReal2D')
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistReal1D')
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistInt3D')
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistInt2D')
      call mpas_dmpar_field_halo_exch(domain, 'cellPersistInt1D')

      ! Exchange persistent edge fields
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistReal5D')
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistReal4D')
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistReal3D')
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistReal2D')
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistReal1D')
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistInt3D')
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistInt2D')
      call mpas_dmpar_field_halo_exch(domain, 'edgePersistInt1D')

      ! Exchange persistent vertex fields
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistReal5D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistReal4D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistReal3D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistReal2D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistReal1D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistInt3D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistInt2D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexPersistInt1D')

      ! Exchange scratch cell fields
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchReal5D')
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchReal4D')
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchReal3D')
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchReal2D')
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchReal1D')
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchInt3D')
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchInt2D')
      call mpas_dmpar_field_halo_exch(domain, 'cellScratchInt1D')

      ! Exchange edge scratch fields
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchReal5D')
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchReal4D')
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchReal3D')
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchReal2D')
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchReal1D')
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchInt3D')
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchInt2D')
      call mpas_dmpar_field_halo_exch(domain, 'edgeScratchInt1D')

      ! Exchange scratch vertex fields
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchReal5D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchReal4D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchReal3D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchReal2D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchReal1D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchInt3D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchInt2D')
      call mpas_dmpar_field_halo_exch(domain, 'vertexScratchInt1D')

      call mpas_timer_stop('single halo exchanges')

      call mpas_threading_barrier()

      call test_core_halo_exch_validate_fields(domain, threadErrs, err)

   end subroutine test_core_halo_exch_single_field_test!}}}


   !***********************************************************************
   !
   !  routine test_core_halo_exch_group_test
   !
   !> \brief   MPAS Test Core halo exchange group test
   !> \author  Doug Jacobsen
   !> \date    01/05/2016
   !> \details 
   !>  This routine performs tests related to halo exchanges using exchange groups.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_halo_exch_group_test(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      character (len=StrKIND) :: groupName
      character (len=StrKIND) :: fieldName

      call mpas_timer_start('build exchange groups')

      groupName = 'cellPersistGroup'
      call mpas_dmpar_exch_group_create(domain, groupName)

      fieldName = 'cellPersistInt1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      fieldName = 'cellPersistInt2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      fieldName = 'cellPersistInt3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      fieldName = 'cellPersistReal1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      fieldName = 'cellPersistReal2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      fieldName = 'cellPersistReal3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      fieldName = 'cellPersistReal4D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      fieldName = 'cellPersistReal5D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)

      groupName = 'edgePersistGroup'
      call mpas_dmpar_exch_group_create(domain, groupName)
 
      fieldName = 'edgePersistInt1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgePersistInt2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgePersistInt3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgePersistReal1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgePersistReal2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgePersistReal3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgePersistReal4D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgePersistReal5D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      groupName = 'vertexPersistGroup'
      call mpas_dmpar_exch_group_create(domain, groupName)
 
      fieldName = 'vertexPersistInt1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexPersistInt2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexPersistInt3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexPersistReal1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexPersistReal2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexPersistReal3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexPersistReal4D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexPersistReal5D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      groupName = 'cellScratchGroup'
      call mpas_dmpar_exch_group_create(domain, groupName)
 
      fieldName = 'cellScratchInt1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'cellScratchInt2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'cellScratchInt3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'cellScratchReal1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'cellScratchReal2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'cellScratchReal3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'cellScratchReal4D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'cellScratchReal5D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      groupName = 'edgeScratchGroup'
      call mpas_dmpar_exch_group_create(domain, groupName)
 
      fieldName = 'edgeScratchInt1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgeScratchInt2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgeScratchInt3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgeScratchReal1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgeScratchReal2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgeScratchReal3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgeScratchReal4D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'edgeScratchReal5D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      groupName = 'vertexScratchGroup'
      call mpas_dmpar_exch_group_create(domain, groupName)
 
      fieldName = 'vertexScratchInt1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexScratchInt2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexScratchInt3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexScratchReal1D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexScratchReal2D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexScratchReal3D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexScratchReal4D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)
 
      fieldName = 'vertexScratchReal5D'
      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName)


      call mpas_timer_stop('build exchange groups')

      call test_core_halo_exch_setup_fields(domain, threadErrs, err)

      call mpas_threading_barrier()

      call mpas_timer_start('group halo exchanges')

      groupName = 'cellPersistGroup'
      call mpas_threading_barrier()
      call mpas_dmpar_exch_group_begin_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_local_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_end_halo_exch(domain, groupName)

      groupName = 'edgePersistGroup'
      call mpas_threading_barrier()
      call mpas_dmpar_exch_group_begin_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_local_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_end_halo_exch(domain, groupName)
 
      groupName = 'vertexPersistGroup'
      call mpas_threading_barrier()
      call mpas_dmpar_exch_group_begin_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_local_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_end_halo_exch(domain, groupName)
 
      groupName = 'cellScratchGroup'
      call mpas_threading_barrier()
      call mpas_dmpar_exch_group_begin_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_local_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_end_halo_exch(domain, groupName)
 
      groupName = 'edgeScratchGroup'
      call mpas_threading_barrier()
      call mpas_dmpar_exch_group_begin_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_local_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_end_halo_exch(domain, groupName)
 
      groupName = 'vertexScratchGroup'
      call mpas_threading_barrier()
      call mpas_dmpar_exch_group_begin_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_local_halo_exch(domain, groupName)
      call mpas_dmpar_exch_group_end_halo_exch(domain, groupName)

      call mpas_timer_stop('group halo exchanges')

      call test_core_halo_exch_validate_fields(domain, threadErrs, err)

   end subroutine test_core_halo_exch_group_test!}}}


   !***********************************************************************
   !
   !  routine test_core_halo_exch_setup_fields
   !
   !> \brief   MPAS Test Core halo exchange field setup routine
   !> \author  Doug Jacobsen
   !> \date    01/06/2016
   !> \details 
   !>  This routine sets up fields for a halo exchange test. It initializes
   !>  their values based on indexTo*ID, where * is cell, edge, or vertex.
   !>  Additionally, it allocates the scratch fields and fills them with the
   !>  correct values.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_halo_exch_setup_fields(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool, haloExchTestPool

      type (field5DReal), pointer :: real5DField
      type (field4DReal), pointer :: real4DField
      type (field3DReal), pointer :: real3DField
      type (field2DReal), pointer :: real2DField
      type (field1DReal), pointer :: real1DField
      type (field3DInteger), pointer :: int3DField
      type (field2DInteger), pointer :: int2DField
      type (field1DInteger), pointer :: int1DField

      real (kind=RKIND), dimension(:, :, :, :, :), pointer :: real5D
      real (kind=RKIND), dimension(:, :, :, :), pointer :: real4D
      real (kind=RKIND), dimension(:, :, :), pointer :: real3D
      real (kind=RKIND), dimension(:, :), pointer :: real2D
      real (kind=RKIND), dimension(:), pointer :: real1D

      real (kind=RKIND) :: realValue
      integer :: integerValue

      integer, dimension(:, :, :), pointer :: int3D
      integer, dimension(:, :), pointer :: int2D
      integer, dimension(:), pointer :: int1D

      integer :: i, j, k, l, m
      integer :: iDim1, iDim2, iDim3, iDim4, iDim5
      integer, pointer :: nCells, nEdges, nVertices
      integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve
      integer, dimension(:), pointer :: indexToCellID
      integer, dimension(:), pointer :: indexToEdgeID
      integer, dimension(:), pointer :: indexToVertexID

      integer :: threadNum

      threadNum = mpas_threading_get_thread_num() + 1

      threadErrs(threadNum) = 0

      ! Allocate all scratch fields
      ! Fill all fields with index values
      block => domain % blocklist
      do while ( associated(block) )
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'haloExchTest', haloExchTestPool)

         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
         call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
         call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

         ! Fill persistent cell fields
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistInt1D', int1D)

         iDim1 = nCellsSolve
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToCellID(i), kind=RKIND)
            integerValue = indexToCellID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = realValue
                     end do
                     real4D(l, k, j, i) = realValue
                  end do
                  real3D(k, j, i) = realValue
                  int3D(k, j, i) = integerValue
               end do
               real2D(j, i) = realValue
               int2D(j, i) = integerValue
            end do
            real1D(i) = realValue
            int1D(i) = integerValue
         end do
         !$omp end do

         ! Fill persistent edge fields
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistInt1D', int1D)

         iDim1 = nEdgesSolve
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToEdgeID(i), kind=RKIND)
            integerValue = indexToEdgeID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = realValue
                     end do
                     real4D(l, k, j, i) = realValue
                  end do
                  real3D(k, j, i) = realValue
                  int3D(k, j, i) = integerValue
               end do
               real2D(j, i) = realValue
               int2D(j, i) = integerValue
            end do
            real1D(i) = realValue
            int1D(i) = integerValue
         end do
         !$omp end do

         ! Fill persistent vertex fields
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistInt1D', int1D)

         iDim1 = nVerticesSolve
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToVertexID(i), kind=RKIND)
            integerValue = indexToVertexID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = realValue
                     end do
                     real4D(l, k, j, i) = realValue
                  end do
                  real3D(k, j, i) = realValue
                  int3D(k, j, i) = integerValue
               end do
               real2D(j, i) = realValue
               int2D(j, i) = integerValue
            end do
            real1D(i) = realValue
            int1D(i) = integerValue
         end do
         !$omp end do

         ! Allocate scratch cell fields
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal5D', real5DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal4D', real4DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal3D', real3DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal2D', real2DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal1D', real1DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt3D', int3DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt2D', int2DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt1D', int1DField)

         call mpas_allocate_scratch_field(real5DField, .true.)
         call mpas_allocate_scratch_field(real4DField, .true.)
         call mpas_allocate_scratch_field(real3DField, .true.)
         call mpas_allocate_scratch_field(real2DField, .true.)
         call mpas_allocate_scratch_field(real1DField, .true.)
         call mpas_allocate_scratch_field(int3DField, .true.)
         call mpas_allocate_scratch_field(int2DField, .true.)
         call mpas_allocate_scratch_field(int1DField, .true.)

         call mpas_threading_barrier()

         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchInt1D', int1D)

         iDim1 = nCellsSolve
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToCellID(i), kind=RKIND)
            integerValue = indexToCellID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = realValue
                     end do
                     real4D(l, k, j, i) = realValue
                  end do
                  real3D(k, j, i) = realValue
                  int3D(k, j, i) = integerValue
               end do
               real2D(j, i) = realValue
               int2D(j, i) = integerValue
            end do
            real1D(i) = realValue
            int1D(i) = integerValue
         end do
         !$omp end do

         ! Allocate scratch edge fields
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal5D', real5DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal4D', real4DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal3D', real3DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal2D', real2DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal1D', real1DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt3D', int3DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt2D', int2DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt1D', int1DField)

         call mpas_allocate_scratch_field(real5DField, .true.)
         call mpas_allocate_scratch_field(real4DField, .true.)
         call mpas_allocate_scratch_field(real3DField, .true.)
         call mpas_allocate_scratch_field(real2DField, .true.)
         call mpas_allocate_scratch_field(real1DField, .true.)
         call mpas_allocate_scratch_field(int3DField, .true.)
         call mpas_allocate_scratch_field(int2DField, .true.)
         call mpas_allocate_scratch_field(int1DField, .true.)

         call mpas_threading_barrier()

         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchInt1D', int1D)

         iDim1 = nEdgesSolve
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToEdgeID(i), kind=RKIND)
            integerValue = indexToEdgeID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = realValue
                     end do
                     real4D(l, k, j, i) = realValue
                  end do
                  real3D(k, j, i) = realValue
                  int3D(k, j, i) = integerValue
               end do
               real2D(j, i) = realValue
               int2D(j, i) = integerValue
            end do
            real1D(i) = realValue
            int1D(i) = integerValue
         end do
         !$omp end do

         ! Allocate scratch vertex fields
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal5D', real5DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal4D', real4DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal3D', real3DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal2D', real2DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal1D', real1DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt3D', int3DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt2D', int2DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt1D', int1DField)

         call mpas_allocate_scratch_field(real5DField, .true.)
         call mpas_allocate_scratch_field(real4DField, .true.)
         call mpas_allocate_scratch_field(real3DField, .true.)
         call mpas_allocate_scratch_field(real2DField, .true.)
         call mpas_allocate_scratch_field(real1DField, .true.)
         call mpas_allocate_scratch_field(int3DField, .true.)
         call mpas_allocate_scratch_field(int2DField, .true.)
         call mpas_allocate_scratch_field(int1DField, .true.)

         call mpas_threading_barrier()

         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchInt1D', int1D)

         iDim1 = nVerticesSolve
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToVertexID(i), kind=RKIND)
            integerValue = indexToVertexID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = realValue
                     end do
                     real4D(l, k, j, i) = realValue
                  end do
                  real3D(k, j, i) = realValue
                  int3D(k, j, i) = integerValue
               end do
               real2D(j, i) = realValue
               int2D(j, i) = integerValue
            end do
            real1D(i) = realValue
            int1D(i) = integerValue
         end do
         !$omp end do

         block => block % next
      end do

      call mpas_threading_barrier()

   end subroutine test_core_halo_exch_setup_fields!}}}

   !***********************************************************************
   !
   !  routine test_core_halo_exch_validate_fields
   !
   !> \brief   MPAS Test Core halo exchange validate fields routine
   !> \author  Doug Jacobsen
   !> \date    01/06/2016
   !> \details 
   !>  This routine validates the fields after a halo exchange. It checks that
   !>  every entry in the field has the expected value, when compared with
   !>  indexTo*ID, when indexTo*ID is cell, edge, or vertex.
   !>  If the validation fails, a 1 is return, and if the validation passes a 0 is return.
   !>  Additionally, this routine deallocates scratch fields.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_halo_exch_validate_fields(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool, haloExchTestPool

      type (field5DReal), pointer :: real5DField
      type (field4DReal), pointer :: real4DField
      type (field3DReal), pointer :: real3DField
      type (field2DReal), pointer :: real2DField
      type (field1DReal), pointer :: real1DField
      type (field3DInteger), pointer :: int3DField
      type (field2DInteger), pointer :: int2DField
      type (field1DInteger), pointer :: int1DField

      real (kind=RKIND), dimension(:, :, :, :, :), pointer :: real5D
      real (kind=RKIND), dimension(:, :, :, :), pointer :: real4D
      real (kind=RKIND), dimension(:, :, :), pointer :: real3D
      real (kind=RKIND), dimension(:, :), pointer :: real2D
      real (kind=RKIND), dimension(:), pointer :: real1D

      real (kind=RKIND) :: realValue
      integer :: integerValue

      integer, dimension(:, :, :), pointer :: int3D
      integer, dimension(:, :), pointer :: int2D
      integer, dimension(:), pointer :: int1D

      integer :: i, j, k, l, m
      integer :: iDim1, iDim2, iDim3, iDim4, iDim5
      integer, pointer :: nCells, nEdges, nVertices
      integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve
      integer, dimension(:), pointer :: indexToCellID
      integer, dimension(:), pointer :: indexToEdgeID
      integer, dimension(:), pointer :: indexToVertexID

      integer :: threadNum

      threadNum = mpas_threading_get_thread_num() + 1

      threadErrs(threadNum) = 0

      ! Validate results
      block => domain % blocklist
      do while ( associated(block) )
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'haloExchTest', haloExchTestPool)

         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
         call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
         call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

         indexToCellID( nCells + 1 ) = 0
         indexToEdgeID( nEdges + 1 ) = 0
         indexToVertexID( nVertices + 1 ) = 0

         ! Compare persistent cell fields
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellPersistInt1D', int1D)

         ! Subtract index from all peristent cell fields
         iDim1 = size(real5D, dim=5)
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToCellID(i), kind=RKIND)
            integerValue = indexToCellID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = real5D(m, l, k, j, i) - realValue
                     end do
                     real4D(l, k, j, i) = real4D(l, k, j, i) - realValue
                  end do
                  real3D(k, j, i) = real3D(k, j, i) - realValue
                  int3D(k, j, i) = int3D(k, j, i) - integerValue
               end do
               real2D(j, i) = real2D(j, i) - realValue
               int2D(j, i) = int2D(j, i) - integerValue
            end do
            real1D(i) = real1D(i) - realValue
            int1D(i) = int1D(i) - integerValue
         end do
         !$omp end do

         ! Validate that all differences are zero.
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Testing persistent cell fields')
#endif
         if ( sum(real5D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real4D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real3D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real2D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real1D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int3D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int2D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int1D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Test result: $i', intArgs=(/threadErrs(threadNum)/))
#endif

         ! Compare persistent edge fields
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgePersistInt1D', int1D)

         ! Subtract index from all peristent edge fields
         iDim1 = size(real5D, dim=5)
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToEdgeID(i), kind=RKIND)
            integerValue = indexToEdgeID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = real5D(m, l, k, j, i) - realValue
                     end do
                     real4D(l, k, j, i) = real4D(l, k, j, i) - realValue
                  end do
                  real3D(k, j, i) = real3D(k, j, i) - realValue
                  int3D(k, j, i) = int3D(k, j, i) - integerValue
               end do
               real2D(j, i) = real2D(j, i) - realValue
               int2D(j, i) = int2D(j, i) - integerValue
            end do
            real1D(i) = real1D(i) - realValue
            int1D(i) = int1D(i) - integerValue
         end do
         !$omp end do

         ! Validate that all differences are zero.
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Testing persistent Edge fields')
#endif
         if ( sum(real5D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real4D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real3D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real2D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real1D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int3D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int2D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int1D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Test result: $i', intArgs=(/threadErrs(threadNum)/))
#endif

         ! Compare persistent vertex fields
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexPersistInt1D', int1D)

         ! Subtract index from all peristent vertex fields
         iDim1 = size(real5D, dim=5)
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToVertexID(i), kind=RKIND)
            integerValue = indexToVertexID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = real5D(m, l, k, j, i) - realValue
                     end do
                     real4D(l, k, j, i) = real4D(l, k, j, i) - realValue
                  end do
                  real3D(k, j, i) = real3D(k, j, i) - realValue
                  int3D(k, j, i) = int3D(k, j, i) - integerValue
               end do
               real2D(j, i) = real2D(j, i) - realValue
               int2D(j, i) = int2D(j, i) - integerValue
            end do
            real1D(i) = real1D(i) - realValue
            int1D(i) = int1D(i) - integerValue
         end do
         !$omp end do

         ! Validate that all differences are zero.
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Testing persistent Vertex fields')
#endif
         if ( sum(real5D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real4D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real3D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real2D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real1D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int3D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Test result: $i', intArgs=(/threadErrs(threadNum)/))
#endif

         ! Compare scratch cell fields
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'cellScratchInt1D', int1D)

         ! Subtract index from all peristent cell fields
         iDim1 = size(real5D, dim=5)
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToCellID(i), kind=RKIND)
            integerValue = indexToCellID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = real5D(m, l, k, j, i) - realValue
                     end do
                     real4D(l, k, j, i) = real4D(l, k, j, i) - realValue
                  end do
                  real3D(k, j, i) = real3D(k, j, i) - realValue
                  int3D(k, j, i) = int3D(k, j, i) - integerValue
               end do
               real2D(j, i) = real2D(j, i) - realValue
               int2D(j, i) = int2D(j, i) - integerValue
            end do
            real1D(i) = real1D(i) - realValue
            int1D(i) = int1D(i) - integerValue
         end do
         !$omp end do

         ! Validate that all differences are zero.
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Testing scratch cell fields')
#endif
         if ( sum(real5D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real4D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real3D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real2D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real1D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int3D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int2D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int1D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Test result: $i', intArgs=(/threadErrs(threadNum)/))
#endif

         ! Compare scratch edge fields
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'edgeScratchInt1D', int1D)

         ! Subtract index from all peristent edge fields
         iDim1 = size(real5D, dim=5)
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m)
         do i = 1, iDim1
            realValue = real(indexToEdgeID(i), kind=RKIND)
            integerValue = indexToEdgeID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim5
                        real5D(m, l, k, j, i) = real5D(m, l, k, j, i) - realValue
                     end do
                     real4D(l, k, j, i) = real4D(l, k, j, i) - realValue
                  end do
                  real3D(k, j, i) = real3D(k, j, i) - realValue
                  int3D(k, j, i) = int3D(k, j, i) - integerValue
               end do
               real2D(j, i) = real2D(j, i) - realValue
               int2D(j, i) = int2D(j, i) - integerValue
            end do
            real1D(i) = real1D(i) - realValue
            int1D(i) = int1D(i) - integerValue
         end do
         !$omp end do

         ! Validate that all differences are zero.
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Testing scratch edge fields')
#endif
         if ( sum(real5D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real4D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real3D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real2D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real1D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int3D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int2D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int1D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Test result: $i', intArgs=(/threadErrs(threadNum)/)
#endif

         ! Compare scratch vertex fields
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal5D', real5D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal4D', real4D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal3D', real3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal2D', real2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchReal1D', real1D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchInt3D', int3D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchInt2D', int2D)
         call mpas_pool_get_array(haloExchTestPool, 'vertexScratchInt1D', int1D)

         ! Subtract index from all peristent vertex fields
         iDim1 = size(real5D, dim=5)
         iDim2 = size(real5D, dim=4)
         iDim3 = size(real5D, dim=3)
         iDim4 = size(real5D, dim=2)
         iDim5 = size(real5D, dim=1)

         !$omp do schedule(runtime) private(j, k, l, m, realValue, integerValue)
         do i = 1, iDim1
            realValue = real(indexToVertexID(i), kind=RKIND)
            integerValue = indexToVertexID(i)
            do j = 1, iDim2
               do k = 1, iDim3
                  do l = 1, iDim4
                     do m = 1, iDim4
                        real5D(m, l, k, j, i) = real5D(m, l, k, j, i) - realValue
                     end do
                     real4D(l, k, j, i) = real4D(l, k, j, i) - realValue
                  end do
                  real3D(k, j, i) = real3D(k, j, i) - realValue
                  int3D(k, j, i) = int3D(k, j, i) - integerValue
               end do
               real2D(j, i) = real2D(j, i) - realValue
               int2D(j, i) = int2D(j, i) - integerValue
            end do
            real1D(i) = real1D(i) - realValue
            int1D(i) = int1D(i) - integerValue
         end do
         !$omp end do

         ! Validate that all differences are zero.
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Testing scratch vertex fields')
#endif
         if ( sum(real5D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real4D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real3D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real2D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(real1D) /= 0.0_RKIND ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int3D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int2D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if

         if ( sum(int1D) /= 0 ) then
            threadErrs(threadNum) = 1
         end if
#ifdef HALO_EXCH_DEBUG
         call mpas_log_write('     -- Test result: $i', intArgs=(/threadErrs(threadNum)/))
#endif

         block => block % next
      end do

      call mpas_threading_barrier()

      ! Deallocatae all scratch fields
      block => domain % blocklist
      do while ( associated(block) )
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         ! Allocate scratch cell fields
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal5D', real5DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal4D', real4DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal3D', real3DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal2D', real2DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchReal1D', real1DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt3D', int3DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt2D', int2DField)
         call mpas_pool_get_field(haloExchTestPool, 'cellScratchInt1D', int1DField)

         call mpas_deallocate_scratch_field(real5DField, .true.)
         call mpas_deallocate_scratch_field(real4DField, .true.)
         call mpas_deallocate_scratch_field(real3DField, .true.)
         call mpas_deallocate_scratch_field(real2DField, .true.)
         call mpas_deallocate_scratch_field(real1DField, .true.)
         call mpas_deallocate_scratch_field(int3DField, .true.)
         call mpas_deallocate_scratch_field(int2DField, .true.)
         call mpas_deallocate_scratch_field(int1DField, .true.)

         ! Allocate scratch edge fields
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal5D', real5DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal4D', real4DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal3D', real3DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal2D', real2DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchReal1D', real1DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt3D', int3DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt2D', int2DField)
         call mpas_pool_get_field(haloExchTestPool, 'edgeScratchInt1D', int1DField)

         call mpas_deallocate_scratch_field(real5DField, .true.)
         call mpas_deallocate_scratch_field(real4DField, .true.)
         call mpas_deallocate_scratch_field(real3DField, .true.)
         call mpas_deallocate_scratch_field(real2DField, .true.)
         call mpas_deallocate_scratch_field(real1DField, .true.)
         call mpas_deallocate_scratch_field(int3DField, .true.)
         call mpas_deallocate_scratch_field(int2DField, .true.)
         call mpas_deallocate_scratch_field(int1DField, .true.)

         ! Allocate scratch vertex fields
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal5D', real5DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal4D', real4DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal3D', real3DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal2D', real2DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchReal1D', real1DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt3D', int3DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt2D', int2DField)
         call mpas_pool_get_field(haloExchTestPool, 'vertexScratchInt1D', int1DField)

         call mpas_deallocate_scratch_field(real5DField, .true.)
         call mpas_deallocate_scratch_field(real4DField, .true.)
         call mpas_deallocate_scratch_field(real3DField, .true.)
         call mpas_deallocate_scratch_field(real2DField, .true.)
         call mpas_deallocate_scratch_field(real1DField, .true.)
         call mpas_deallocate_scratch_field(int3DField, .true.)
         call mpas_deallocate_scratch_field(int2DField, .true.)
         call mpas_deallocate_scratch_field(int1DField, .true.)

         block => block % next
      end do

      err = sum(threadErrs)

   end subroutine test_core_halo_exch_validate_fields!}}}


end module test_core_halo_exch
