#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module utils_mod

implicit none
private
      
       public :: InsertIntoArray
       public :: RemoveFromArray

contains 

   subroutine LinearFind(array,value,found,indx)
      integer, intent(in)  :: array(:)
      integer, intent(in)  :: value
      logical, intent(out) :: found
      integer, intent(out) :: indx

      integer :: n,nz,i
      n = SIZE(array)
      ! =============================================
      ! Array of all zeros... insert at the begining
      ! =============================================
      found = .FALSE.
      if((array(1) == 0) .or. (value < array(1)) ) then
        found = .FALSE.
        indx = 1
        return
      endif
      nz = COUNT(array .gt. 0)
      if(nz == n) then
          print *,'LinearFind... No free spaace in array'
          return
      endif

      do i=1,nz
         if(array(i) == value) then
            ! ===============
            ! Already in list
            ! ===============
            found = .TRUE.
            indx = i 
            return
         else if((array(i) < value) .and. (value < array(i+1)) ) then
            ! =====================================
            ! Insert it into the middle of the list
            ! =====================================
            found = .FALSE.
            indx = i+1
            return
         else if((array(i) < value) .and. (array(i+1) == 0) )then
            ! =====================================
            ! Insert it into the end of the List
            ! =====================================
            found = .FALSE.
            indx = i+1
            return
         endif
      enddo

   end subroutine LinearFind

   subroutine InsertIntoArray(array,value,ierr)
        integer, intent(inout) :: array(:)
        integer, intent(in)    :: value
        integer, intent(out) :: ierr
        integer, allocatable :: tmp(:)
        logical  :: found
        integer  :: indx
        integer  :: n,i,nz
        logical, parameter  :: Debug = .FALSE.

        call LinearFind(array,value,found,indx)

        if(Debug) print *,'InsertIntoArray: array is :',array
        if(Debug) print *,'InsertIntoArray: value is :',value
        if(Debug) print *,'InsertIntoArray: found is :',found
        if(Debug) print *,'InsertIntoArray: indx is :',indx

        if(.not. found) then

          n = SIZE(array)
          allocate(tmp(n))
          tmp = array

          if(indx== 1) then
             array(1) = value
             array(2:n) = tmp(1:n-1)
          else
             array(1:indx-1) = tmp(1:indx-1)
             array(indx) = value
             array(indx+1:n) = tmp(indx:n-1)
          endif

          deallocate(tmp)

        endif

   end subroutine InsertIntoArray

   subroutine RemoveFromArray(array,value,ierr)
        integer, intent(inout) :: array(:)
        integer, intent(in)    :: value
        integer, intent(out)   :: ierr
        integer, allocatable   :: tmp(:)

        integer :: it,ia,n

        n = SIZE(array)
        allocate(tmp(n))
        tmp = array

        array = 0
        it = 1
        ia = 1
        do while (it < n)
           if(tmp(it) .ne.  value ) then
                array(ia) = tmp(it)
                ia = ia+1
                it = it+1
           else
                it = it+1
           endif
        enddo

        deallocate(tmp)
               
   end subroutine RemoveFromArray


end module utils_mod
