module test_shr_mpi_mod
  use shr_mpi_mod,  only: shr_mpi_gathScatVInit, &
                          shr_mpi_gatherV,       &
                          shr_mpi_scatterv,      &
                          shr_mpi_commrank,      &
                          shr_mpi_chkerr,        &
                          shr_mpi_commsize,      &
                          shr_mpi_send,          &
                          shr_mpi_recv,          &
                          shr_mpi_barrier
  use shr_kind_mod, only: r8 => SHR_KIND_R8
  use shr_sys_mod,  only: shr_sys_abort
  implicit none
#include <mpif.h>

  private

  public :: test_gathScat
  public :: test_gathScatDiffPES

  contains

logical function test_gathScat( mpicom, rootid, locArr )
   use shr_kind_mod,  only: SHR_KIND_IN
   use shr_const_mod, only: SHR_CONST_SPVAL
   implicit none
   integer(SHR_KIND_IN), intent(IN) :: mpicom
   integer(SHR_KIND_IN), intent(IN) :: rootid
   real(r8), pointer    :: locArr(:)

   real(r8), pointer :: glob1DArr(:), glob1DArrBack(:)
   integer(SHR_KIND_IN),  pointer :: globSize(:), displs(:)
   integer(SHR_KIND_IN),  pointer :: globSizeBack(:), displsBack(:)
   real(r8), pointer :: locArrBack(:)
   integer :: rank, npes, ierr
   logical, pointer :: results(:)

   if ( .not. associated(locArr) )then
      test_gathScat = .false.
      return
   end if
   allocate( locArrBack(size(locArr)) )
   locArrBack(:) = SHR_CONST_SPVAL
   call shr_mpi_gathScatvInit( mpicom, rootid, locArr, glob1DArr, globSize, displs )
   call shr_mpi_gathScatvInit( mpicom, rootid, locArrBack, glob1DArrBack, &
                               globSizeBack, displsBack )
   call shr_mpi_gatherv( locarr, size(locArr), glob1DArr, globSize, displs, rootid, &
                         mpicom )
   call shr_mpi_commrank( mpicom, rank )
   call shr_mpi_commsize( mpicom, npes )
   if ( rank == rootid ) glob1DArrBack(:) = glob1DArr(:)
   call shr_mpi_scatterv( locarrBack, size(locArrBack), glob1DArrBack, globSizeBack,    &
                          displsBack, rootid, mpicom )
   ! Test that original local array and array from gather/scatter are same
   if ( all(locArr == locArrBack) .and. all(locArrBack /= SHR_CONST_SPVAL) )then
      test_gathScat = .true.
   else
      test_gathScat = .false.
   end if
   ! Now check that global arrays are the same after the gather
   if ( rank == rootid .and. test_gathScat ) glob1DArrBack(:) = SHR_CONST_SPVAL
   call shr_mpi_gatherv( locarr, size(locArr), glob1DArrBack, globSize, displs, rootid, &
                         mpicom )
   if ( rank == rootid .and. test_gathScat )then
      if ( all(glob1DArr(:) == glob1DArrBack(:)) .and. all(glob1DArrBack(:) /= SHR_CONST_SPVAL) )then
         test_gathScat = .true.
      else
         test_gathScat = .false.
      end if
   end if
   deallocate( glob1DArr,     globSize,     displs     )
   deallocate( glob1DArrBack, globSizeBack, displsBack )
   return
end function test_gathScat

logical function test_gathScatDiffPES( mpicom, mpicom2, rootid, locArr )
   use shr_kind_mod,  only: SHR_KIND_IN
   use shr_const_mod, only: SHR_CONST_SPVAL
   implicit none
   integer(SHR_KIND_IN), intent(IN) :: mpicom
   integer(SHR_KIND_IN), intent(IN) :: mpicom2
   integer(SHR_KIND_IN), intent(IN) :: rootid
   real(r8), pointer    :: locArr(:)

   real(r8), pointer :: glob1DArr(:)
   integer(SHR_KIND_IN),  pointer :: globSize(:), displs(:)
   integer :: rank, npes, ierr, rank2, npes2, nsize, i
   integer, pointer :: lsize(:)
   logical, pointer :: results(:)
   real(r8), pointer    :: locArr2(:)
   real(r8), pointer :: glob1DArr2(:)
   integer(SHR_KIND_IN),  pointer :: globSize2(:), displs2(:)

   if ( .not. associated(locArr) )then
      test_gathScatDiffPES = .false.
      return
   end if
   ! First gather the local array into a global array that you keep
   call shr_mpi_gathScatvInit( mpicom, rootid, locArr, glob1DArr, globSize, displs )
   call shr_mpi_gatherv( locarr, size(locArr), glob1DArr, globSize, displs, rootid, &
                         mpicom )
   ! Then scatter/gather using the other communicator -- make sure global array identical
   call shr_mpi_commrank( mpicom, rank )
   if ( mpicom2 /= MPI_COMM_NULL )then
      call shr_mpi_commsize( mpicom2, npes2 )
      ! Figure out size for each local array and send to each processor in group
      if ( rank == rootid )then
         nsize = size(glob1DArr) / npes2
         allocate( lsize(0:npes2-1) )
         lsize(0:npes2-2) = nsize
         lsize(npes2-1)   = size(glob1DArr) - sum(lsize(0:npes2-2))
         do i = 1, npes2-1
            write(6,*) "lsize, peid = ", lsize(i), i
            call shr_mpi_send( lsize(i), i, 1055, mpicom2 )
         end do
         deallocate( lsize )
      else
         call shr_mpi_recv( nsize, rootid, 1055, mpicom2 )
      end if
      allocate( locArr2(nsize) )
      call shr_mpi_gathScatvInit( mpicom2, rootid, locArr2, glob1DArr2, globSize2, &
                                  displs2 )
      call shr_mpi_scatterv( locarr2, size(locArr2), glob1DArr, globSize2,    &
                          displs2, rootid, mpicom2 )
      glob1DArr2(:) = SHR_CONST_SPVAL
      call shr_mpi_gatherv( locarr2, size(locArr2), glob1DArr2, globSize2, displs2, &
                            rootid, mpicom2 )
      call shr_mpi_commrank( mpicom, rank2 )
      if ( (rank == rootid) .and. (rank2 == rootid) )then
         if ( all(glob1DArr(:) == glob1DArr2(:)) .and. &
              all(glob1DArr2(:) /= SHR_CONST_SPVAL) )then
            test_gathScatDiffPES = .true.
         else
            test_gathScatDiffPES = .false.
         end if
      end if
      deallocate( glob1DArr2,   globSize2,    displs2    )
   end if
   deallocate( glob1DArr,    globSize,     displs     )
   return
end function test_gathScatDiffPES

end module test_shr_mpi_mod

program test_shr_mpi

  use test_shr_mpi_mod, only: test_gathScat, test_gathScatDiffPES
  use shr_mpi_mod,      only: shr_mpi_init,          &
                              shr_mpi_finalize,      &
                              shr_mpi_commrank,      &
                              shr_mpi_commsize,      &
                              shr_mpi_chkerr,        &
                              shr_mpi_barrier
  use shr_kind_mod,     only: r8 => SHR_KIND_R8
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  implicit none
#include <mpif.h>
  integer            :: mpicom = MPI_COMM_WORLD
  integer, parameter :: rootid = 0
  real(r8), pointer  :: locArr(:)
  integer            :: i, gsize, rank, npes, npe1, npe2
  integer,  pointer  :: seed(:)
  integer            :: seedSize
  character(len=80)  :: TestType
  real(r8) :: x
  logical :: masterproc
  integer :: mpicom1, mpicom2
  integer :: mpigrp, mpigrp1, mpigrp2, ierr

  call shr_mpi_init( )
  call shr_mpi_commrank( mpicom, rank )
  call shr_mpi_commsize( mpicom, npes )
  masterproc = rank == rootid
  if ( masterproc ) write(6,*) "shr_mpi_mod unit test"
  call random_seed( size=seedSize )
  allocate( seed(seedSize) )
  seed(:) = rank*1000 + 1444
  call random_seed( put=seed        )
  deallocate( seed )
  ! Get communicators for a subset of the processors
  if ( npes > 3 )then
      ! Create new groups of 1 and 2 processors
      ! Must include rank 0 in both...
      call mpi_comm_group( mpicom, mpigrp, ierr )
      call shr_mpi_chkerr( ierr, "Error getting mpi group" )
      call mpi_group_incl( mpigrp, 1, (/0/), mpigrp1, ierr )
      call shr_mpi_chkerr( ierr, "Error getting mpi group-1" )
      call mpi_comm_create( mpicom, mpigrp1, mpicom1, ierr )
      call shr_mpi_chkerr( ierr, "Error creating new comm group with 1 processor" )
      call mpi_group_incl( mpigrp, 2, (/0,2/), mpigrp2, ierr )
      call shr_mpi_chkerr( ierr, "Error getting mpi group-2" )
      call mpi_comm_create( mpicom, mpigrp2, mpicom2, ierr )
      call shr_mpi_chkerr( ierr, "Error creating new comm group with 2 processors" )
      ! Initialize gather/scatter for new communicator groups
      call shr_mpi_barrier( mpicom  )
      if ( mpicom1 /= MPI_COMM_NULL )then
         call shr_mpi_barrier( mpicom1 )
         call shr_mpi_commsize( mpicom1, npe1 )
         if ( npe1 /= 1 )  call shr_sys_abort( "mpicom1 wrong size" )
      end if
      if ( mpicom2 /= MPI_COMM_NULL )then
          call shr_mpi_barrier( mpicom2 )
          call shr_mpi_commsize( mpicom2, npe2 )
          if ( npe2 /= 2 )  call shr_sys_abort( "mpicom2 wrong size" )
      end if
  end if
  do i = 1, 4
     if (      i == 1 )then
        TestType = "same sizes, random values"
        gsize = 10
        call fillArrayRandom( gsize, locArr )
     else if ( i == 2 )then
        TestType = "same sizes, ordered values"
        gsize = 100
        call fillArrayOrdered( gsize, locArr, rank )
     else if ( i == 3 )then
        TestType = "random sizes, random values"
        call random_number( x )
        gsize = nint( x*100._r8 ) + 100
        call fillArrayRandom( gsize, locArr )
     else if ( i == 4 )then
        TestType = "random sizes, ordered values"
        call random_number( x )
        gsize = nint( x*200._r8 ) + 50
        call fillArrayOrdered( gsize, locArr, rank )
     else
        call shr_sys_abort( "Bad index number for test" )
     end if
     if ( masterproc ) write(6,*) "Gather/scatter test for: ", trim(TestType)
     write(6,*) 'rank, size, locarr = ', rank, gsize, locArr
     call shr_sys_flush(6)
     if ( .not. test_gathScat( mpicom, rootid, locArr ) )then
        call shr_sys_abort( "Error in doing scatter/gather" )
     end if
     call shr_mpi_barrier( mpicom )
     if ( masterproc ) write(6,*) "PASS"
     if ( npes > 3 )then
        if ( masterproc ) write(6,*) "Gather/scatter test on mpicom1 for: ", trim(TestType)
        call shr_sys_flush(6)
        if ( .not. test_gathScatDiffPES( mpicom, mpicom1, rootid, locArr ) )then
           call shr_sys_abort( "Error in reconstructing array with mpicom1" )
        end if
        call shr_mpi_barrier( mpicom )
        if ( masterproc ) write(6,*) "PASS"
        call shr_mpi_barrier( mpicom )
        if ( masterproc ) write(6,*) "PASS"
        if ( masterproc ) write(6,*) "Gather/scatter test on mpicom2 for: ", trim(TestType)
        call shr_sys_flush(6)
        if ( .not. test_gathScatDiffPES( mpicom, mpicom2, rootid, locArr ) )then
           call shr_sys_abort( "Error in reconstructing array with mpicom2" )
        end if
        call shr_mpi_barrier( mpicom )
        if ( masterproc ) write(6,*) "PASS"
     end if
     deallocate( locArr )
  end do
  call shr_mpi_finalize( )
  if ( masterproc ) write(6,*) "SUCCESS!"
  if ( masterproc ) write(6,*) "PASS"

contains

subroutine fillArrayRandom( gsize, locArr )
  integer, intent(in) :: gsize
  real(r8), pointer :: locArr(:)

  real(r8) :: x
  integer :: g

  allocate( locArr(gsize) )
  do g = 1, gsize
     call random_number( x )
     locArr(g) = x * 1000.0_r8
  end do
end subroutine fillArrayRandom

subroutine fillArrayOrdered( gsize, locArr, rank )
  integer, intent(in) :: gsize
  integer, intent(in) :: rank
  real(r8), pointer :: locArr(:)

  real(r8) :: x
  integer :: g

  allocate( locArr(gsize) )
  do g = 1, gsize
     locArr(g) = real( g, r8 ) + rank*1000.0_r8
  end do
end subroutine fillArrayOrdered

end program test_shr_mpi
