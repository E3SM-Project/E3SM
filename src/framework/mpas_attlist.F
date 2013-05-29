! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_attlist
!
!> \brief   MPAS Attribute list module
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This module provides type definitions and subroutines for working with attribute lists.
!
!-----------------------------------------------------------------------

module mpas_attlist

   use mpas_kind_types

   ! Derived type for holding field attributes
   type att_list_type
      character (len=StrKIND) :: attName
      integer :: attType
      integer :: attValueInt
      integer, dimension(:), pointer :: attValueIntA => null()
      real (kind=RKIND) :: attValueReal
      real (kind=RKIND), dimension(:), pointer :: attValueRealA => null()
      character (len=StrKIND) :: attValueText
      type (att_list_type), pointer :: next => null()
   end type att_list_type

   interface mpas_add_att
      module procedure mpas_add_att_int0d
      module procedure mpas_add_att_int1d
      module procedure mpas_add_att_real0d
      module procedure mpas_add_att_real1d
      module procedure mpas_add_att_text
   end interface mpas_add_att

   interface mpas_get_att
      module procedure mpas_get_att_int0d
      module procedure mpas_get_att_int1d
      module procedure mpas_get_att_real0d
      module procedure mpas_get_att_real1d
      module procedure mpas_get_att_text
   end interface mpas_get_att


   !!!!! PRIVATE? !!!!!

   integer, parameter :: ATT_INT   = 1
   integer, parameter :: ATT_INTA  = 2
   integer, parameter :: ATT_REAL  = 3
   integer, parameter :: ATT_REALA = 4
   integer, parameter :: ATT_TEXT  = 5


contains

!***********************************************************************
!
!  routine mpas_add_att_int0d
!
!> \brief   MPAS Add 0D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 0D integer attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine mpas_add_att_int0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_INT
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueInt = attValue

   end subroutine mpas_add_att_int0d!}}}

!***********************************************************************
!
!  routine mpas_add_att_int1d
!
!> \brief   MPAS Add 1D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 1D integer attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine mpas_add_att_int1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, dimension(:), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_INTA
      allocate(cursor % attValueIntA(size(attValue)))
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueIntA(:) = attValue(:)

   end subroutine mpas_add_att_int1d!}}}

!***********************************************************************
!
!  routine mpas_add_att_real0d
!
!> \brief   MPAS Add 0D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 0D real attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine mpas_add_att_real0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_REAL
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueReal = attValue

   end subroutine mpas_add_att_real0d!}}}

!***********************************************************************
!
!  routine mpas_add_att_real1d
!
!> \brief   MPAS Add 1D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 1D real attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine mpas_add_att_real1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), dimension(:), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_REALA
      allocate(cursor % attValueRealA(size(attValue)))
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueRealA(:) = attValue(:)

   end subroutine mpas_add_att_real1d!}}}

!***********************************************************************
!
!  routine mpas_add_att_text
!
!> \brief   MPAS Add text attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a text attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine mpas_add_att_text(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      character (len=*), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_TEXT
      write(cursor % attName,'(a)') trim(attName)
      write(cursor % attValueText,'(a)') trim(attValue)

   end subroutine mpas_add_att_text!}}}

!***********************************************************************
!
!  routine mpas_get_att_int0d
!
!> \brief   MPAS get 0D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 0D integer attribute.
!
!-----------------------------------------------------------------------
   subroutine mpas_get_att_int0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, intent(out) :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_INT) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               attValue = cursor % attValueInt
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_int0d!}}}

!***********************************************************************
!
!  routine mpas_get_att_int1d
!
!> \brief   MPAS get 1D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 1D integer attribute.
!
!-----------------------------------------------------------------------
   subroutine mpas_get_att_int1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, dimension(:), pointer :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_INTA) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               allocate(attValue(size(cursor % attValueIntA)))
               attValue(:) = cursor % attValueIntA(:)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_int1d!}}}

!***********************************************************************
!
!  routine mpas_get_att_real0d
!
!> \brief   MPAS get 0D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 0D real attribute.
!
!-----------------------------------------------------------------------
   subroutine mpas_get_att_real0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), intent(out) :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_REAL) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               attValue = cursor % attValueReal
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_real0d!}}}

!***********************************************************************
!
!  routine mpas_get_att_real1d
!
!> \brief   MPAS get 1D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 1D real attribute.
!
!-----------------------------------------------------------------------
   subroutine mpas_get_att_real1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), dimension(:), pointer :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_REALA) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               allocate(attValue(size(cursor % attValueRealA)))
               attValue(:) = cursor % attValueRealA(:)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_real1d!}}}

!***********************************************************************
!
!  routine mpas_get_att_text
!
!> \brief   MPAS get text attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a text attribute.
!
!-----------------------------------------------------------------------
   subroutine mpas_get_att_text(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      character (len=*), intent(out) :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_TEXT) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               write(attValue,'(a)') trim(cursor % attValueText)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_text!}}}

!***********************************************************************
!
!  routine mpas_remove_att
!
!> \brief   MPAS remove attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine removes an attribute from an attribute list.
!
!-----------------------------------------------------------------------
   subroutine mpas_remove_att(attList, attName, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor, cursor_prev

      if (present(ierr)) ierr = 0

      cursor => attList

      ! Item is at the head of the list
      if (trim(attName) == trim(cursor % attName)) then
         attList => cursor % next
         if (cursor % attType == ATT_REALA) then
            deallocate(cursor % attValueRealA)
         else if (cursor % attType == ATT_INTA) then
            deallocate(cursor % attValueIntA)
         end if
         deallocate(cursor)
         return
      end if

      cursor_prev => cursor
      cursor => cursor % next
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            cursor_prev % next => cursor % next

            if (cursor % attType == ATT_REALA) then
               deallocate(cursor % attValueRealA)
            else if (cursor % attType == ATT_INTA) then
               deallocate(cursor % attValueIntA)
            end if
            deallocate(cursor)
            
            return
         end if 

         cursor_prev => cursor
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_remove_att!}}}

!***********************************************************************
!
!  routine mpas_deallocate_attlist
!
!> \brief   MPAS attribute list deallocation routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine deallocates an attribute list.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_attlist(attList, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         attList => attList % next
         if (cursor % attType == ATT_REALA) then
            deallocate(cursor % attValueRealA)
         else if (cursor % attType == ATT_INTA) then
            deallocate(cursor % attValueIntA)
         end if
         deallocate(cursor)
         cursor => attList
      end do

   end subroutine mpas_deallocate_attlist!}}}

end module mpas_attlist
