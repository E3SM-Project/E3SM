
      subroutine SUB_SCAN( filecnt, &
                           lib_files, &
                           usr_paths, &
                           usr_files, &
                           sub_cnt )
!---------------------------------------------------------------------
!	... This subroutine scans a "standard" library routine for
!           matching user routines which will override the library
!           routine.
!---------------------------------------------------------------------

      implicit none

!---------------------------------------------------------------------
!	... Dummy args
!---------------------------------------------------------------------
      integer, intent(in)    ::     filecnt          ! count of library subroutines
      integer, intent(inout) ::     sub_cnt          ! count of "user" subroutines
      character(len=*), intent(inout) ::  lib_files(*)  ! library filenames
      character(len=*), intent(inout) ::  usr_files(*)  ! user filenames
      character(len=*), intent(inout) ::  usr_paths(*)  ! user filepaths

!---------------------------------------------------------------------
!	... Local variables
!---------------------------------------------------------------------
      integer ::   i, j, alls, cnt
      integer, allocatable :: mark(:)
      character(len=128), allocatable :: wrk_files(:)  ! work space
      character(len=128), allocatable :: wrk_paths(:)  ! work space
      logical, allocatable :: keep(:)

      if( filecnt > 0 .and. sub_cnt > 0 ) then
	 ALLOCATE( wrk_files(sub_cnt), wrk_paths(sub_cnt), keep(sub_cnt), stat = alls )
	 if( alls /= 0 ) then
	    write(*,*) ' SUB_SCAN : Failed to allocated wrk array'
	    stop 'Alloc err'
	 end if
	 ALLOCATE( mark(filecnt), stat = alls )
	 if( alls /= 0 ) then
	    write(*,*) ' SUB_SCAN : Failed to allocated wrk array'
	    stop 'Alloc err'
	 end if
	 wrk_files(:sub_cnt) = usr_files(:sub_cnt)
	 wrk_paths(:sub_cnt) = usr_paths(:sub_cnt)
	 do i = 1,filecnt
	    mark(i) = INDEX( lib_files(i)(:LEN_TRIM(lib_files(i))), &
			     '/', back = .true. ) + 1
	 end do
	 keep(:sub_cnt) = .true.
         do i = 1,filecnt
	    do j = 1,sub_cnt
	       if( keep(j) .and. &
	           lib_files(i)(mark(i):LEN_TRIM(lib_files(i))) == &
		   wrk_files(j)(:LEN_TRIM(wrk_files(j))) ) then
	          lib_files(i) = wrk_paths(j)(:LEN_TRIM(usr_paths(j))) &
			         // wrk_files(j)(:LEN_TRIM(wrk_files(j)))
		  keep(j) = .false.
	          exit
	       end if
	    end do
         end do

         cnt = COUNT( keep(:sub_cnt) )
         if( cnt /= sub_cnt ) then
	    usr_files(:cnt) = PACK( wrk_files(:sub_cnt), mask = keep(:sub_cnt) )
	    usr_paths(:cnt) = PACK( wrk_paths(:sub_cnt), mask = keep(:sub_cnt) )
	    sub_cnt = cnt
	 end if
	 DEALLOCATE( wrk_files )
	 DEALLOCATE( wrk_paths )
	 DEALLOCATE( keep )
	 DEALLOCATE( mark )
      end if

      end subroutine SUB_SCAN
