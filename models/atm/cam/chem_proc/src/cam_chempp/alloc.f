
      module REALLOC
!-----------------------------------------------------------------------
!        ... Reallocation module for integer, real, and character arrays
!-----------------------------------------------------------------------

      implicit none

      interface REALLOCATE
         module PROCEDURE R3_ALLOCATE, R2_ALLOCATE, R1_ALLOCATE
     $,                   I3_ALLOCATE, I2_ALLOCATE, I1_ALLOCATE
     $,                   C2_ALLOCATE, C1_ALLOCATE
      end interface

      CONTAINS

      subroutine R3_ALLOCATE( array )
!-----------------------------------------------------------------------
!        ... Real realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      real, pointer :: array(:,:,:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz(3)
      real, allocatable :: temp(:,:,:)

      siz(1) = SIZE( array,1 )
      siz(2) = SIZE( array,2 )
      siz(3) = SIZE( array,3 )
      ALLOCATE( temp(siz(1),siz(2),siz(3)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'R3_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '             Size = ',siz(1),' x ',siz(2),' x ',siz(3)
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz(1),siz(2),siz(3)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'R3_ALLOCATE: Failed to allocate new array'
	 write(*,*) '             Size = ',siz(1),' x ',siz(2),' x ',siz(3)
	 stop
      end if
      array = 0.
      array(:siz(1),:,:) = temp(:siz(1),:,:)
      DEALLOCATE( temp )
 
      end subroutine R3_ALLOCATE

      subroutine R2_ALLOCATE( array )
!-----------------------------------------------------------------------
!        ... Real realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      real, pointer :: array(:,:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz(2)
      real, allocatable :: temp(:,:)

      siz(1) = SIZE( array,1 )
      siz(2) = SIZE( array,2 )
      ALLOCATE( temp(siz(1),siz(2)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'R2_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '             Size = ',siz(1),' x ',siz(2)
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz(1),siz(2)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'R2_ALLOCATE: Failed to allocate new array'
	 write(*,*) '             Size = ',2*siz(1),' x ',siz(2)
	 stop
      end if
      array = 0.
      array(:siz(1),:) = temp(:siz(1),:)
      DEALLOCATE( temp )
 
      end subroutine R2_ALLOCATE

      subroutine R1_ALLOCATE( array )
!-----------------------------------------------------------------------
!        ... Real realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      real, pointer :: array(:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz
      real, allocatable :: temp(:)

      siz = SIZE( array )
      ALLOCATE( temp(siz),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'R1_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '             Size = ',siz
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'R1_ALLOCATE: Failed to allocate new array'
	 write(*,*) '             Size = ',2*siz
	 stop
      end if
      array = 0.
      array(:siz) = temp(:siz)
      DEALLOCATE( temp )
 
      end subroutine R1_ALLOCATE

      subroutine I3_ALLOCATE( array )
!-----------------------------------------------------------------------
!        ... Integer realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, pointer :: array(:,:,:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz(3)
      integer, allocatable :: temp(:,:,:)

      siz(1) = SIZE( array,1 )
      siz(2) = SIZE( array,2 )
      siz(3) = SIZE( array,3 )
      ALLOCATE( temp(siz(1),siz(2),siz(3)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'I3_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '             Size = ',siz(1),' x ',siz(2),' x ',siz(3)
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz(1),siz(2),siz(3)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'I3_ALLOCATE: Failed to allocate new array'
	 write(*,*) '             Size = ',siz(1),' x ',siz(2),' x ',siz(3)
	 stop
      end if
      array = 0
      array(:siz(1),:,:) = temp(:siz(1),:,:)
      DEALLOCATE( temp )
 
      end subroutine I3_ALLOCATE

      subroutine I2_ALLOCATE( array )
!-----------------------------------------------------------------------
!        ... Integer realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, pointer :: array(:,:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz(2)
      integer, allocatable :: temp(:,:)

      siz(1) = SIZE( array,1 )
      siz(2) = SIZE( array,2 )
      ALLOCATE( temp(siz(1),siz(2)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'I2_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '             Size = ',siz(1),' x ',siz(2)
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz(1),siz(2)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'I2_ALLOCATE: Failed to allocate new array'
	 write(*,*) '             Size = ',2*siz(1),' x ',siz(2)
	 stop
      end if
      array = 0
      array(:siz(1),:) = temp(:siz(1),:)
      DEALLOCATE( temp )
 
      end subroutine I2_ALLOCATE

      subroutine I1_ALLOCATE( array )
!-----------------------------------------------------------------------
!        ... Integer realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, pointer :: array(:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz
      integer, allocatable :: temp(:)

      siz = SIZE( array )
      ALLOCATE( temp(siz),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'I1_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '            Size = ',siz
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'I1_ALLOCATE: Failed to allocate new array'
	 write(*,*) '            Size = ',2*siz
	 stop
      end if
      array = 0
      array(:siz) = temp(:siz)
      DEALLOCATE( temp )
 
      end subroutine I1_ALLOCATE

      subroutine C1_ALLOCATE( array, clen )
!-----------------------------------------------------------------------
!        ... Character realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: clen                     ! length of array elements
      character(len=clen), pointer :: array(:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz
      character(len=clen), allocatable :: temp(:)

      siz = SIZE( array )
      ALLOCATE( temp(siz),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'C1_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '            Size = ',siz
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'C1_ALLOCATE: Failed to allocate new array'
	 write(*,*) '            Size = ',2*siz
	 stop
      end if
      array = ' '
      array(:siz) = temp(:siz)
      DEALLOCATE( temp )
 
      end subroutine C1_ALLOCATE

      subroutine C2_ALLOCATE( array, clen )
!-----------------------------------------------------------------------
!        ... Character realloc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: clen                     ! length of array elements
      character(len=clen), pointer :: array(:,:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat
      integer :: siz(2)
      character(len=clen), allocatable :: temp(:,:)

      siz(1) = SIZE( array,1 )
      siz(2) = SIZE( array,2 )
      ALLOCATE( temp(siz(1),siz(2)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'C2_ALLOCATE: Failed to allocate temp array'
	 write(*,*) '             Size = ',siz(1),' x ',siz(2)
	 stop
      end if
      temp = array
      DEALLOCATE( array )
      ALLOCATE( array(2*siz(1),siz(2)),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'C2_ALLOCATE: Failed to allocate new array'
	 write(*,*) '             Size = ',2*siz(1),' x ',siz(2)
	 stop
      end if
      array = ' '
      array(:siz(1),:) = temp(:siz(1),:)
      DEALLOCATE( temp )
 
      end subroutine C2_ALLOCATE

      end module REALLOC
