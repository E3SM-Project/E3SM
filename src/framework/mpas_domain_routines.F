! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_domain_routines
!
!> \brief   MPAS Domain routines
!> \author  Doug Jacobsen, Michael Duda
!> \date    03/10/2015
!> \details 
!> This module defines routines related to the domain derived data type, and
!> the block derived data type.
!
!-----------------------------------------------------------------------
module mpas_domain_routines

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar

   contains

!***********************************************************************
!
!  routine mpas_allocate_domain
!
!> \brief   MPAS Domain allocation routine
!> \author  Michael Duda
!> \date    04/02/13
!> \details 
!> This routine allocates a domain structure.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_domain(dom)!{{{

      implicit none

      type (domain_type), pointer :: dom !< Input/Output: Domain structure

      allocate(dom % dminfo)
      nullify(dom % blocklist)

      allocate(dom % configs)
      allocate(dom % packages)
      allocate(dom % clock)
      allocate(dom % streamManager)
      allocate(dom % ioContext)

      call mpas_pool_create_pool(dom % configs)
      call mpas_pool_create_pool(dom % packages)

   end subroutine mpas_allocate_domain!}}}


!***********************************************************************
!
!  routine mpas_allocate_block
!
!> \brief   MPAS Block allocation routine
!> \author  Michael Duda
!> \date    04/02/13
!> \details 
!> This routine allocates a block structure. It calls routines to allocate the variable structures
!> that are members of the block type.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_block(nHaloLayers, b, dom, blockID) !{{{

      implicit none

      integer, intent(in) :: nHaloLayers !< Input: Number of halo laters
      type (block_type), pointer :: b !< Input/Output: Block structure
      type (domain_type), pointer :: dom !< Input: Domain structure
      integer, intent(in) :: blockID !< Input: Global ID of block

      integer :: i

      b % blockID = blockID

      allocate(b % parinfo)

      b % domain => dom

      allocate(b % structs)
      allocate(b % dimensions)
      allocate(b % allFields)
      call mpas_pool_create_pool(b % structs)
      call mpas_pool_create_pool(b % dimensions)
      call mpas_pool_create_pool(b % allFields)
      call mpas_pool_create_pool(b % allStructs)

      b % configs => dom % configs
      b % packages => dom % packages

   end subroutine mpas_allocate_block!}}}


!***********************************************************************
!
!  routine mpas_deallocate_domain
!
!> \brief   MPAS Domain deallocation routine
!> \author  Michael Duda
!> \date    04/02/13
!> \details 
!> This routine deallocates a domain structure. 
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_domain(dom)!{{{

      implicit none

      type (domain_type), pointer :: dom !< Input/Output: Domain to deallocate

      type (block_type), pointer :: block_ptr
      type (mpas_exchange_group), pointer :: exchGroupPtr

      block_ptr => dom % blocklist
      do while (associated(block_ptr))
         call mpas_deallocate_block(block_ptr)
         block_ptr => block_ptr % next
      end do

      exchGroupPtr => dom % exchangeGroups
      do while (associated(exchGroupPtr))
         call mpas_dmpar_exch_group_destroy(dom, exchGroupPtr % groupName)
         exchGroupPtr => dom % exchangeGroups
      end do

      call mpas_pool_destroy_pool(dom % configs)
      call mpas_pool_destroy_pool(dom % packages)
      deallocate(dom % clock)
      deallocate(dom % ioContext)

   end subroutine mpas_deallocate_domain!}}}


!***********************************************************************
!
!  routine mpas_deallocate_block
!
!> \brief   MPAS Block deallocation routine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a block structure.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_block(b)!{{{
 
      implicit none

      type (block_type), intent(inout) :: b !< Input/Output: Block to be deallocated.

      integer :: i

      ! BUG: It seems like we should be deallocating the exchange lists before we 
      !      deallocate the array of head pointers and the parinfo type...
      !      It also seems like these deallocations should happen with mpas_dmpar_destroy_multihalo_exchange_list

      call mpas_pool_destroy_pool(b % structs)
      call mpas_pool_destroy_pool(b % dimensions)

      deallocate(b % parinfo % cellsToSend)
      deallocate(b % parinfo % cellsToRecv)
      deallocate(b % parinfo % cellsToCopy)

      deallocate(b % parinfo % edgesToSend)
      deallocate(b % parinfo % edgesToRecv)
      deallocate(b % parinfo % edgesToCopy)

      deallocate(b % parinfo % verticesToSend)
      deallocate(b % parinfo % verticesToRecv)
      deallocate(b % parinfo % verticesToCopy)

      deallocate(b % parinfo)

   end subroutine mpas_deallocate_block!}}}


!***********************************************************************
!
!  routine mpas_link_fields
!
!> \brief   MPAS Link Fields routine
!> \author  Doug Jacobsen
!> \date    04/24/2014
!> \details 
!> Sets up field links across blocks. This routine should be called prior to
!> any halo exchanges. Also links parinfo into each field.
!
!-----------------------------------------------------------------------
   subroutine mpas_link_fields(domain)!{{{
      type (domain_type), intent(in) :: domain

      type (block_type), pointer :: blockPtr
      type (block_type), pointer :: prevBlock
      type (block_type), pointer :: nextBlock

      blockPtr => domain % blocklist
      do while(associated(blockPtr))
         if (associated(blockPtr % prev)) then
            prevBlock => blockPtr % prev
         else
            nullify(prevBlock)
         end if

         if (associated(blockPtr % next)) then
            nextBlock => blockPtr % next
         else
            nullify(nextBlock)
         end if

         if (associated(prevBlock) .and. associated(nextBlock)) then
            call mpas_pool_link_pools(blockPtr % structs, prevBlock % structs, nextBlock % structs)
         else if (associated(prevBlock)) then
            call mpas_pool_link_pools(blockPtr % structs, prevBlock % structs)
         else if (associated(nextBlock)) then
            call mpas_pool_link_pools(blockPtr % structs, nextPool=nextBlock % structs)
         else
            call mpas_pool_link_pools(blockPtr % structs)
         end if

         call mpas_pool_link_parinfo(blockPtr, blockPtr % structs)

         blockPtr => blockPtr % next
      end do

   end subroutine mpas_link_fields!}}}


!***********************************************************************
!
!  funtion mpas_dimension_num_garbage_elements
!
!> \brief   MPAS Number of garbage elements for dimension
!> \author  Doug Jacobsen
!> \date    02/28/2015
!> \details 
!> This routine takes as input the name of a dimension and a dimension pool.
!> Using this information, it returns the number of garbage elements that should
!> be padded on each array allocated with this dimension size.
!
!-----------------------------------------------------------------------
   integer function mpas_dimension_num_garbage_elements(dimName, iErr) result(numElements)!{{{
      character (len=*), intent(in) :: dimName
      integer, intent(in) :: iErr

      numElements = 0

      if ( trim(dimName) == "nCells" ) then
         numElements = 1
      else if ( trim(dimName) == "nVertices" ) then
         numElements = 1
      else if ( trim(dimName) == "nEdges" ) then
         numElements = 1
      end if

   end function mpas_dimension_num_garbage_elements!}}}

end module mpas_domain_routines
