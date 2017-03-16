! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_decomp
!
!> \brief MPAS Decomposition list
!> \author Doug Jacobsen
!> \date   03/10/2015
!> \details
!>  This module will contain the mpas_decomp_list type, along with routines to register and
!>  access a decomposition function.
!
!-----------------------------------------------------------------------

#define COMMA ,
#define DECOMP_DEBUG_WRITE(M) ! call mpas_log_write(M)
#define DECOMP_WARN_WRITE(M) call mpas_log_write( M , messageType=MPAS_LOG_WARN)
#define DECOMP_ERROR_WRITE(M) call mpas_log_write( M , messageType=MPAS_LOG_CRIT)

module mpas_decomp

   use mpas_kind_types
   use mpas_derived_types
   use mpas_stream_manager
   use mpas_log

   implicit none

   contains


!-----------------------------------------------------------------------
!  routine mpas_decomp_create_decomp_list
!
!> \brief MPAS decomp list create routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine creates a decomposition list.
!
!-----------------------------------------------------------------------
   subroutine mpas_decomp_create_decomp_list(decompList)!{{{
      type (mpas_decomp_list), pointer :: decompList

      integer :: errLocal

      DECOMP_DEBUG_WRITE('* Creating decomposition list.')

      if ( .not. associated(decompList) ) then
         allocate(decompList)
         decompList % decompName = ''
         decompList % nameLen = -1
         nullify(decompList % decompFunc)
         nullify(decompList % next)
      end if

   end subroutine!}}}


!-----------------------------------------------------------------------
!  routine mpas_decomp_destroy_decomp_list
!
!> \brief MPAS decomp list destruction routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine destroys a decomposition list.
!
!-----------------------------------------------------------------------
   subroutine mpas_decomp_destroy_decomp_list(decompList)!{{{
      type (mpas_decomp_list), pointer :: decompList

      type (mpas_decomp_list), pointer :: decompCursor

      DECOMP_DEBUG_WRITE('* Destroying decomposition list.')

      do while (associated(decompList % next))
         decompCursor => decompList % next

         if ( associated(decompCursor % next) ) then
            decompList % next => decompCursor % next
         else
            nullify(decompList % next)
         end if

         deallocate(decompCursor)
      end do

      deallocate(decompList)
   end subroutine mpas_decomp_destroy_decomp_list!}}}


!-----------------------------------------------------------------------
!  routine mpas_decomp_register_method
!
!> \brief MPAS decomp list register routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine registers a decomposition function in a decomposition list.
!>  If present, the optional argument iErr is set to MPAS_DECOMP_ERROR if 
!>  the function is already registered and to MPAS_DECOMP_NOERR otherwise.
!
!-----------------------------------------------------------------------
   subroutine mpas_decomp_register_method(decompList, decompName, decompFunc, iErr)!{{{
      type (mpas_decomp_list), pointer :: decompList
      character (len=*), intent(in) :: decompName
      procedure (mpas_decomp_function), pointer :: decompFunc
      integer, intent(out) :: iErr

      type (mpas_decomp_list), pointer :: decompCursor
      type (mpas_decomp_list), pointer :: decompCursorPrev
      integer :: nameLen

      iErr = MPAS_DECOMP_NOERR

      nameLen = len_trim(decompName)

      decompCursor => decompList
      nullify(decompCursorPrev)
      do while (associated(decompCursor))
         if ( nameLen == decompCursor % nameLen ) then
            if ( trim(decompName) == trim(decompCursor % decompName) ) then
               DECOMP_DEBUG_WRITE('Decomposition ' // trim(decompName) // ' is already registered.')
               iErr = MPAS_DECOMP_ERROR
               return
            end if
         end if
         decompCursorPrev => decompCursor
         decompCursor => decompCursor % next
      end do

      decompCursor => decompCursorPrev

      DECOMP_DEBUG_WRITE('* Registering decomposition '//trim(decompName))

      allocate(decompCursor % next)
      decompCursor => decompCursor % next
      decompCursor % decompName = decompName
      decompCursor % decompFunc => decompFunc
      decompCursor % nameLen = nameLen
      nullify(decompCursor % next)
   end subroutine mpas_decomp_register_method!}}}


!-----------------------------------------------------------------------
!  routine mpas_decomp_get_method
!
!> \brief MPAS decomp list query routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine querys a decomposition list for a specific decomposition function.
!>  If present, the optional argument iErr is set to MPAS_DECOMP_ERROR if 
!>  the requested function is not found and to MPAS_DECOMP_NOERR otherwise.
!
!-----------------------------------------------------------------------
   subroutine mpas_decomp_get_method(decompList, decompName, decompFunc, iErr)!{{{
      type (mpas_decomp_list), pointer :: decompList
      character (len=*), intent(in) :: decompName
      procedure (mpas_decomp_function), pointer :: decompFunc
      integer, intent(out) :: iErr

      type (mpas_decomp_list), pointer :: decompCursor
      integer :: nameLen

      iErr = MPAS_DECOMP_NOERR

      nameLen = len_trim(decompName)

      decompCursor => decompList
      do while (associated(decompCursor))
         if ( nameLen == decompCursor % nameLen ) then
            if ( trim(decompCursor % decompName) == trim(decompName) ) then
               DECOMP_DEBUG_WRITE('* Returning decomposition '//trim(decompName))
               decompFunc => decompCursor % decompFunc
               return
            end if
         end if
         decompCursor => decompCursor % next
      end do

      iErr = MPAS_DECOMP_ERROR
      DECOMP_DEBUG_WRITE('Decomposition '//trim(decompName)//' not found in decomposition list.')

   end subroutine mpas_decomp_get_method!}}}


!-----------------------------------------------------------------------
!  routine mpas_decomp_remove_method
!
!> \brief MPAS decomp list removal routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine removes a specific decomposition function from a  decomposition list.
!
!-----------------------------------------------------------------------
   subroutine mpas_decomp_remove_method(decompList, decompName, iErr)!{{{
      type (mpas_decomp_list), pointer :: decompList
      character (len=*), intent(in) :: decompName
      integer, intent(out) :: iErr

      integer :: nameLen
      type (mpas_decomp_list), pointer :: decompCursor, decompCursorPrev

      iErr = MPAS_DECOMP_NOERR

      nameLen = len_trim(decompName)

      decompCursor => decompList
      nullify(decompCursorPrev)
      do while (associated(decompCursor))
         if ( nameLen == decompCursor % nameLen ) then
            if ( trim(decompCursor % decompName) == trim(decompName) ) then
               DECOMP_DEBUG_WRITE('* Removing decomposition '//trim(decompName))

               if ( associated(decompCursor % next) ) then
                  decompCursorPrev % next => decompCursor % next
               else
                  nullify(decompCursorPrev % next)
               end if

               deallocate(decompCursor)
               return
            end if
         end if
         decompCursorPrev => decompCursor
         decompCursor => decompCursor % next
      end do
   end subroutine mpas_decomp_remove_method!}}}


!-----------------------------------------------------------------------
!  function mpas_uniform_decomp
!
!> \brief MPAS uniform decomposition function
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This function provides the owned indices for a uniform decomposition method.
!
!-----------------------------------------------------------------------
   function mpas_uniform_decomp(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)!{{{
      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, intent(in) :: globalDimSize
      integer, intent(in) :: numBlocks

      integer, dimension(:), pointer :: ownedIndices

      integer :: iErr

      integer :: blockDimSize, blockStart, i

      iErr = 0

      DECOMP_DEBUG_WRITE('* Providing a uniform decomposition')

      blockDimSize = globalDimSize / numBlocks
      if ( block % blockID < mod(globalDimSize, numBlocks) ) then
         blockDimSize = blockDimSize + 1
      end if

      blockStart = block % blockID * blockDimSize
      if ( block % blockID >= mod(globalDimSize, numBlocks) ) then
         blockStart = blockStart + mod(globalDimSize, numBlocks)
      end if

      allocate(ownedIndices(blockDimSize))
      do i = 1, blockDimSize
         ownedIndices(i) = blockStart + i
      end do

   end function mpas_uniform_decomp!}}}

end module mpas_decomp
