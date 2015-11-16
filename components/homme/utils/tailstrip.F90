! =============================================================
! A routine useful in stripping tails from filenames, etc.
! =============================================================

      subroutine tailstrip(name,tail)
      implicit none
!------------------------------Arguments--------------------------------
!
! Input
!
      character(len=*) name       ! Input character string
!
! Output
!
      character(len=*) tail        ! Returned character string between first two dots found
!
!---------------------------Local variables-----------------------------
!
      integer i               ! Loop index
      integer istart, iend    ! start and end characters of 
!
!-----------------------------------------------------------------------
!
! Find first . character in string
!
      istart = 1
      do i = 1,LEN( name )
        if (name(i:i).eq.'.') then
           istart = i+1
           go to 10
       end if
      end do
 10   continue

!
! Find second . character in string
!

      iend = istart
      do i=istart,LEN(name)
        if (name(i:i).eq.'.') then
           iend = i-1
           go to 20
        endif
      end do
      if (iend .eq. istart) iend = len(name)
 20   continue

      tail(1:(iend-istart+1))=name(istart:iend)
      tail = tail(1:(iend-istart+1))//""

      return
      end
