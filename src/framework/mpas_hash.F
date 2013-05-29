! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_hash
!
!> \brief   MPAS Hash table module
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This module provides A dictionary/hashtable with insert, search, and remove routines.
!
!-----------------------------------------------------------------------
module mpas_hash

   ! Parameters
   integer, parameter :: TABLESIZE=27183     !< Number of spaces in the table (the number of linked lists)
 
   type hashnode
      integer :: key
      type (hashnode), pointer :: next
   end type hashnode
 
   type hashnode_ptr
      type (hashnode), pointer :: p        !< Pointer to a list of entries
   end type hashnode_ptr
 
   type hashtable
      integer :: size
      type (hashnode_ptr), dimension(TABLESIZE) :: table !< The hashtable array
   end type hashtable

 
   contains
 
!***********************************************************************
!
!  routine mpas_hash_init
!
!> \brief   MPAS Hash table init routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine initializes a hashtable.
!
!-----------------------------------------------------------------------
   subroutine mpas_hash_init(h)!{{{
   
     implicit none
 
     ! Arguments
     type (hashtable), intent(inout) :: h !< Input/Output: Hash table
 
     ! Local variables
     integer :: i

     h%size = 0
 
     do i=1,TABLESIZE
        nullify(h%table(i)%p)
     end do
 
   end subroutine mpas_hash_init!}}}
 
 
!***********************************************************************
!
!  routine mpas_hash_insert
!
!> \brief   MPAS Hash table insert routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine inserts a key into a hashtable. If the key already exists in the hash table,
!> a second copy of the key is added to the table.
!
!-----------------------------------------------------------------------
   subroutine mpas_hash_insert(h, key)!{{{
   
     implicit none
 
     ! Arguments
     integer, intent(in) :: key !< Input: Key
     type (hashtable), intent(inout) :: h !< Input/Output: Hashtable
 
     ! Local variables
     integer :: hashval, i
     type (hashnode), pointer :: hn 
 
     hashval = mod(key, TABLESIZE) + 1  
    
     allocate(hn) 
     hn%key = key
     hn%next => h%table(hashval)%p
     h%table(hashval)%p => hn 

     h%size = h%size + 1
 
   end subroutine mpas_hash_insert!}}}
 
!***********************************************************************
!
!  logical function mpas_hash_search
!
!> \brief   MPAS Hash table search routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This function searches for a key within a hashtable. If the key is found TRUE is returned, and FALSE is return otherwise.
!
!-----------------------------------------------------------------------
   logical function mpas_hash_search(h, key)!{{{
   
      implicit none
  
      ! Arguments
      integer, intent(in) :: key !< Input: Key
      type (hashtable), intent(inout) :: h !< Input/Output: Hashtable
  
      ! Local variables
      integer :: hashval, i
      type (hashnode), pointer :: cursor 
  
      mpas_hash_search = .false.
  
      hashval = mod(key, TABLESIZE) + 1  

      cursor => h%table(hashval)%p
      do while(associated(cursor))
         if (cursor%key == key) then
            mpas_hash_search = .true.
            return 
         else
            cursor => cursor%next 
         end if
      end do
    
      return
 
   end function mpas_hash_search


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: hash_size
   !
   ! Purpose: Returns the number of items in the hash table h.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer function mpas_hash_size(h)

      implicit none

      ! Arguments
      type (hashtable) :: h

      mpas_hash_size = h%size

      return

   end function mpas_hash_size
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: hash_destroy
   !
   ! Purpose: Frees all memory associated with hashtable h. This routine may be
   !   used to remove all entries from a hashtable.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine mpas_hash_destroy(h)
   
      implicit none
  
      ! Arguments
      type (hashtable), intent(inout) :: h
  
      ! Local variables
      integer :: i
      type (hashnode), pointer :: cursor, cursor_prev
     
      do i=1,TABLESIZE
         cursor => h%table(i)%p
         do while(associated(cursor))
            cursor_prev => cursor
            cursor => cursor%next
            deallocate(cursor_prev)
         end do 
         nullify(h%table(i)%p)
      end do 

      h%size = 0
 
   end subroutine mpas_hash_destroy!}}}
 
end module mpas_hash
