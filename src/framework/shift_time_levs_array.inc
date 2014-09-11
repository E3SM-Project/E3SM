      !!!!!!!!!!!!!!!!!!!!!!!
      ! Implementation note:
      !
      ! In this subroutine, we use an array of fields as a ready-made array 
      !    of field pointers; these pointers exist in the field types as 'next' pointers
      !!!!!!!!!!!!!!!!!!!!!!!

  
      nlevs = size(fldarr)       
      allocate(fldarr_ptr(nlevs))

      !
      ! Initialize pointers to first block of all time levels
      !
      do i=1,nlevs
         fldarr_ptr(i) % next => fldarr(i)
      end do


      !
      ! Loop over all blocks
      !
      do while (associated(fldarr_ptr(1) % next))

         !
         ! Shift time levels for this block
         !
         arr_ptr => fldarr_ptr(1) % next % array
         do i=1,nlevs-1
            fldarr_ptr(i) % next % array => fldarr_ptr(i+1) % next % array
         end do
         fldarr_ptr(nlevs) % next % array => arr_ptr

         ! Advance pointers to next block
         do i=1,nlevs
            fldarr_ptr(i) % next => fldarr_ptr(i) % next % next
         end do
      end do

      deallocate(fldarr_ptr)
