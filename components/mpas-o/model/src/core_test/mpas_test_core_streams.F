! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module test_core_streams

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_field_routines
   use mpas_dmpar
   use mpas_log


   contains


   !***********************************************************************
   !
   !  routine test_core_streams_test
   !
   !> \brief   tests reading/writing single- and double-precision streams
   !> \author  Michael Duda
   !> \date    2 February 2016
   !> \details 
   !>  This routine verifies that both single- and double-precision streams
   !>  can be written via the mpas_stream_manager module. Only real-valued
   !>  fields are read/written, under the assumption that any precision
   !>  changes would not affect logical, character, or integer data.
   !>
   !>  It is assumed that there is a var_struct containing fields 
   !>  cellPersistReal{0,1,2,3,4,5}D that can be added to the streams
   !>  that are created by this routine.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_streams_test(domain, threadErrs, ierr)

      use mpas_stream_manager

      implicit none

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: ierr

      integer :: i, j, k, l, m
      integer :: iDim1, iDim2, iDim3, iDim4, iDim5
      integer :: local_ierr
      integer, pointer :: nCellsSolve
      integer, dimension(:), pointer :: indexToCellID
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: haloExchTestPool

      real (kind=RKIND), dimension(:, :, :, :, :), pointer :: real5D
      real (kind=RKIND), dimension(:, :, :, :), pointer :: real4D
      real (kind=RKIND), dimension(:, :, :), pointer :: real3D
      real (kind=RKIND), dimension(:, :), pointer :: real2D
      real (kind=RKIND), dimension(:), pointer :: real1D
      real (kind=RKIND), pointer :: real0D
      real (kind=RKIND) :: realValue
      real (kind=RKIND), dimension(5) :: d


      ierr = 0

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'haloExchTest', haloExchTestPool)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal5D', real5D)
      call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal4D', real4D)
      call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal3D', real3D)
      call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal2D', real2D)
      call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal1D', real1D)
      call mpas_pool_get_array(haloExchTestPool, 'cellPersistReal0D', real0D)


      !
      ! Initialize fields before writing them out
      !
      iDim1 = nCellsSolve
      iDim2 = size(real5D, dim=4)
      iDim3 = size(real5D, dim=3)
      iDim4 = size(real5D, dim=2)
      iDim5 = size(real5D, dim=1)

      do i = 1, iDim1
         realValue = real(indexToCellID(i), kind=RKIND)/10.0_RKIND
         do j = 1, iDim2
            do k = 1, iDim3
               do l = 1, iDim4
                  do m = 1, iDim5
                     real5D(m, l, k, j, i) = realValue
                  end do
                  real4D(l, k, j, i) = realValue
               end do
               real3D(k, j, i) = realValue
            end do
            real2D(j, i) = realValue
         end do
         real1D(i) = realValue
      end do
      real0D = 2.0_RKIND * asin(1.0_RKIND)


      !
      ! Create output streams and write real-valued fields in both r4 and r8 precision
      !
      call MPAS_stream_mgr_create_stream(domain % streamManager, 'R4_stream', MPAS_STREAM_OUTPUT, 'r4_data.nc', &
                                         realPrecision=MPAS_IO_SINGLE_PRECISION, &
                                         clobberMode=MPAS_STREAM_CLOBBER_TRUNCATE, &
                                         ierr=local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error creating ''R4_stream''.')
         return
      end if

      call MPAS_stream_mgr_create_stream(domain % streamManager, 'R8_stream', MPAS_STREAM_OUTPUT, 'r8_data.nc', &
                                         realPrecision=MPAS_IO_DOUBLE_PRECISION, &
                                         clobberMode=MPAS_STREAM_CLOBBER_TRUNCATE, &
                                         ierr=local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error creating ''R8_stream''.')
         return
      end if

      call MPAS_stream_mgr_create_stream(domain % streamManager, 'R8_time_stream', MPAS_STREAM_OUTPUT, 'r8_data.$Y-$M-$D-$d_$h.$m.$s.nc', &
                                         filenameInterval="0001-00-00_00:00:00", &
                                         realPrecision=MPAS_IO_DOUBLE_PRECISION, &
                                         clobberMode=MPAS_STREAM_CLOBBER_TRUNCATE, &
                                         ierr=local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error creating ''R8_stream''.')
         return
      end if


      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'xtime', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal5D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal4D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal3D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal2D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal1D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal0D', ierr=local_ierr)

      call MPAS_stream_mgr_add_stream_fields(domain % streamManager, 'R8_stream', 'R4_stream', ierr=local_ierr)

      call MPAS_stream_mgr_add_stream_fields(domain % streamManager, 'R8_time_stream', 'R4_stream', ierr=local_ierr)

      call MPAS_stream_mgr_write(domain % streamManager, 'R4_stream', forceWriteNow=.true., ierr=local_ierr)
      call MPAS_stream_mgr_write(domain % streamManager, 'R8_stream', forceWriteNow=.true., ierr=local_ierr)
      call MPAS_stream_mgr_write(domain % streamManager, 'R8_time_stream', forceWriteNow=.true., ierr=local_ierr)
      call MPAS_stream_mgr_write(domain % streamManager, 'R8_time_stream', forceWriteNow=.true., writeTime="9999-01-01_00:00:00", ierr=local_ierr)


      !
      ! Destroy streams
      !
      call MPAS_stream_mgr_destroy_stream(domain % streamManager, 'R4_stream', local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error destroying ''R4_stream''.', MPAS_LOG_ERR)
         return
      end if

      call MPAS_stream_mgr_destroy_stream(domain % streamManager, 'R8_stream', local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error destroying ''R8_stream''.', MPAS_LOG_ERR)
         return
      end if




      !
      ! Create input streams and read real-valued fields in both r4 and r8 precision
      !
      call MPAS_stream_mgr_create_stream(domain % streamManager, 'R4_stream', MPAS_STREAM_INPUT, 'r4_data.nc', &
                                         realPrecision=MPAS_IO_SINGLE_PRECISION, &
                                         clobberMode=MPAS_STREAM_CLOBBER_TRUNCATE, &
                                         ierr=local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error creating ''R4_stream''.', MPAS_LOG_ERR)
         return
      end if

      call MPAS_stream_mgr_create_stream(domain % streamManager, 'R8_stream', MPAS_STREAM_INPUT, 'r8_data.nc', &
                                         realPrecision=MPAS_IO_DOUBLE_PRECISION, &
                                         clobberMode=MPAS_STREAM_CLOBBER_TRUNCATE, &
                                         ierr=local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error creating ''R8_stream''.', MPAS_LOG_ERR)
         return
      end if

      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'xtime', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal5D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal4D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal3D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal2D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal1D', ierr=local_ierr)
      call MPAS_stream_mgr_add_field(domain % streamManager, 'R4_stream', 'cellPersistReal0D', ierr=local_ierr)

      call MPAS_stream_mgr_add_stream_fields(domain % streamManager, 'R8_stream', 'R4_stream', ierr=local_ierr)

      !
      ! Reset contents of real arrays before reading back from streams
      !
      do i = 1, iDim1
         realValue = -1.0
         do j = 1, iDim2
            do k = 1, iDim3
               do l = 1, iDim4
                  do m = 1, iDim5
                     real5D(m, l, k, j, i) = realValue
                  end do
                  real4D(l, k, j, i) = realValue
               end do
               real3D(k, j, i) = realValue
            end do
            real2D(j, i) = realValue
         end do
         real1D(i) = realValue
      end do
      real0D = -1.0

      call MPAS_stream_mgr_read(domain % streamManager, 'R4_stream', rightNow=.true., ierr=local_ierr)

      !
      ! Verify contents of real arrays after reading back from streams
      !
      d(:) = 0.0
      do i = 1, iDim1
         realValue = real(real(indexToCellID(i), kind=RKIND)/10.0_RKIND,kind=R4KIND)
         do j = 1, iDim2
            do k = 1, iDim3
               do l = 1, iDim4
                  do m = 1, iDim5
                     d(5) = d(5) + (real5D(m, l, k, j, i) - realValue)
                  end do
                  d(4) = d(4) + (real4D(l, k, j, i) - realValue)
               end do
               d(3) = d(3) + (real3D(k, j, i) - realValue)
            end do
            d(2) = d(2) + (real2D(j, i) - realValue)
         end do
         d(1) = d(1) + (real1D(i) - realValue)
      end do
      if (real0D /= real(2.0_RKIND * asin(1.0_RKIND),kind=R4KIND)) then
         call mpas_log_write('   Difference detected when reading back 0-d field from single-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 0-d field from single-precision stream - SUCCESS')
      end if
      if (d(1) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 1-d field from single-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 1-d field from single-precision stream - SUCCESS')
      end if
      if (d(2) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 2-d field from single-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 2-d field from single-precision stream - SUCCESS')
      end if
      if (d(3) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 3-d field from single-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 3-d field from single-precision stream - SUCCESS')
      end if
      if (d(4) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 4-d field from single-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 4-d field from single-precision stream - SUCCESS')
      end if
      if (d(5) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 5-d field from single-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 5-d field from single-precision stream - SUCCESS')
      end if

      !
      ! Reset contents of real arrays before reading back from streams
      !
      do i = 1, iDim1
         realValue = -1.0
         do j = 1, iDim2
            do k = 1, iDim3
               do l = 1, iDim4
                  do m = 1, iDim5
                     real5D(m, l, k, j, i) = realValue
                  end do
                  real4D(l, k, j, i) = realValue
               end do
               real3D(k, j, i) = realValue
            end do
            real2D(j, i) = realValue
         end do
         real1D(i) = realValue
      end do
      real0D = -1.0

      call MPAS_stream_mgr_read(domain % streamManager, 'R8_stream', rightNow=.true., ierr=local_ierr)

      !
      ! Verify contents of real arrays after reading back from streams
      !
      d(:) = 0.0
      do i = 1, iDim1
         realValue = real(indexToCellID(i), kind=RKIND)/10.0_RKIND
         do j = 1, iDim2
            do k = 1, iDim3
               do l = 1, iDim4
                  do m = 1, iDim5
                     d(5) = d(5) + (real5D(m, l, k, j, i) - realValue)
                  end do
                  d(4) = d(4) + (real4D(l, k, j, i) - realValue)
               end do
               d(3) = d(3) + (real3D(k, j, i) - realValue)
            end do
            d(2) = d(2) + (real2D(j, i) - realValue)
         end do
         d(1) = d(1) + (real1D(i) - realValue)
      end do
      if (real0D /= (2.0_RKIND * asin(1.0_RKIND))) then
         call mpas_log_write('   Difference detected when reading back 0-d field from double-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 0-d field from double-precision stream - SUCCESS')
      end if
      if (d(1) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 1-d field from double-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 1-d field from double-precision stream - SUCCESS')
      end if
      if (d(2) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 2-d field from double-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 2-d field from double-precision stream - SUCCESS')
      end if
      if (d(3) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 3-d field from double-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 3-d field from double-precision stream - SUCCESS')
      end if
      if (d(4) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 4-d field from double-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 4-d field from double-precision stream - SUCCESS')
      end if
      if (d(5) /= 0.0) then
         call mpas_log_write('   Difference detected when reading back 5-d field from double-precision stream - FAILURE', MPAS_LOG_ERR)
         ierr = ierr + 1
      else
         call mpas_log_write('   Reading 5-d field from double-precision stream - SUCCESS')
      end if


      !
      ! Destroy streams
      !
      call MPAS_stream_mgr_destroy_stream(domain % streamManager, 'R4_stream', local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error destroying ''R4_stream''.', MPAS_LOG_ERR)
         return
      end if

      call MPAS_stream_mgr_destroy_stream(domain % streamManager, 'R8_stream', local_ierr)
      if (local_ierr /= MPAS_STREAM_MGR_NOERR) then
         ierr = 1
         call mpas_log_write('Error destroying ''R8_stream''.', MPAS_LOG_ERR)
         return
      end if

   end subroutine test_core_streams_test

end module test_core_streams
