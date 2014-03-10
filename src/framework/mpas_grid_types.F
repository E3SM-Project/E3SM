! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_grid_types
!
!> \brief   MPAS Grid and field type defintion module
!> \author  Michael Duda, Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This module defines derived data types related to fields, and variable structures.
!> It also includes routines for allocating and deallocating these types.
!
!-----------------------------------------------------------------------
module mpas_grid_types

   use mpas_kind_types
   use mpas_dmpar_types
   use mpas_attlist
   use mpas_packages

   integer, parameter :: nTimeLevs = 2

   ! Derived type describing info for doing I/O specific to a field
   type io_info
      character (len=StrKIND) :: fieldName
      character (len=StrKIND) :: units
      character (len=StrKIND) :: description
      integer, dimension(4) :: start
      integer, dimension(4) :: count
      logical :: input
      logical :: sfc
      logical :: restart
      logical :: output
   end type io_info

   ! Derived type for storing fields
   type field5DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(5) :: dimNames
      integer, dimension(5) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field5DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field5DReal


   ! Derived type for storing fields
   type field4DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(4) :: dimNames
      integer, dimension(4) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field4DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field4DReal



   ! Derived type for storing fields
   type field3DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DReal


   ! Derived type for storing fields
   type field2DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DReal


   ! Derived type for storing fields
   type field1DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DReal


   ! Derived type for storing fields
   type field0DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND) :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DReal


   ! Derived type for storing fields
   type field3DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DInteger


   ! Derived type for storing fields
   type field2DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DInteger


   ! Derived type for storing fields
   type field1DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DInteger


   ! Derived type for storing fields
   type field0DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DInteger


   ! Derived type for storing fields
   type field1DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DChar


   ! Derived type for storing fields
   type field0DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND) :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DChar


   ! Derived type for storing grid meta-data
   type mesh_type

      type (block_type), pointer :: block

#include "field_dimensions.inc"

      logical :: on_a_sphere
      real (kind=RKIND) :: sphere_radius

#include "time_invariant_fields.inc"

   end type mesh_type


#include "variable_groups.inc"


   ! Type for storing (possibly architecture specific) information concerning to parallelism
   type parallel_info
      type (mpas_multihalo_exchange_list), pointer :: cellsToSend            ! List of types describing which cells to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToRecv            ! List of types describing which cells to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToCopy            ! List of types describing which cells to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: edgesToSend            ! List of types describing which edges to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToRecv            ! List of types describing which edges to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToCopy            ! List of types describing which edges to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: verticesToSend         ! List of types describing which vertices to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToRecv         ! List of types describing which vertices to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToCopy         ! List of types describing which vertices to copy from other blocks
   end type parallel_info


   ! Derived type for storing part of a domain; used as a basic unit of work for a process
   type block_type

#include "block_group_members.inc"

      integer :: blockID   ! Unique global ID number for this block
      integer :: localBlockID  ! Unique local ID number for this block

      type (domain_type), pointer :: domain

      type (parallel_info), pointer :: parinfo

      type (block_type), pointer :: prev, next
   end type block_type


   ! Derived type for storing list of blocks from a domain to be handled by a process
   type domain_type
      type (block_type), pointer :: blocklist
   
      ! Also store parallelization info here
      type (dm_info), pointer :: dminfo
#include "model_variables.inc"
      character (len=StrKIND*2) :: history !< History attribute, read in from input file.
      character (len=StrKIND) :: mesh_id !< mesh_id attribute, randomly generated
      character (len=StrKIND) :: Conventions !< Conventions attribute, read in from input file.
      character (len=StrKIND) :: source !< source attribute, read in from input file.
      character (len=StrKIND) :: mesh_spec !< mesh_spec attribute, read in from input file.
      character (len=StrKIND) :: parent_id !< parent_id attribute, read in from input file.
      character (len=StrKIND) :: git_version = 'MPAS_GIT_VERSION' !< version from git describe.
   end type domain_type

   interface mpas_allocate_scratch_field
     module procedure mpas_allocate_scratch_field1d_integer
     module procedure mpas_allocate_scratch_field2d_integer
     module procedure mpas_allocate_scratch_field3d_integer
     module procedure mpas_allocate_scratch_field1d_real
     module procedure mpas_allocate_scratch_field2d_real
     module procedure mpas_allocate_scratch_field3d_real
     module procedure mpas_allocate_scratch_field4d_real
     module procedure mpas_allocate_scratch_field5d_real
     module procedure mpas_allocate_scratch_field1d_char
   end interface

   interface mpas_deallocate_scratch_field
     module procedure mpas_deallocate_scratch_field1d_integer
     module procedure mpas_deallocate_scratch_field2d_integer
     module procedure mpas_deallocate_scratch_field3d_integer
     module procedure mpas_deallocate_scratch_field1d_real
     module procedure mpas_deallocate_scratch_field2d_real
     module procedure mpas_deallocate_scratch_field3d_real
     module procedure mpas_deallocate_scratch_field4d_real
     module procedure mpas_deallocate_scratch_field5d_real
     module procedure mpas_deallocate_scratch_field1d_char
   end interface

   interface mpas_deallocate_field
     module procedure mpas_deallocate_field0d_integer
     module procedure mpas_deallocate_field1d_integer
     module procedure mpas_deallocate_field2d_integer
     module procedure mpas_deallocate_field3d_integer
     module procedure mpas_deallocate_field0d_real
     module procedure mpas_deallocate_field1d_real
     module procedure mpas_deallocate_field2d_real
     module procedure mpas_deallocate_field3d_real
     module procedure mpas_deallocate_field4d_real
     module procedure mpas_deallocate_field5d_real
     module procedure mpas_deallocate_field0d_char
     module procedure mpas_deallocate_field1d_char
   end interface


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
   subroutine mpas_allocate_domain(dom, dminfo)!{{{

      implicit none

      type (domain_type), pointer :: dom !< Input/Output: Domain structure
      type (dm_info), pointer :: dminfo !< Input: Domain Information

      nullify(dom % blocklist)
      dom % dminfo => dminfo

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
   subroutine mpas_allocate_block(nHaloLayers, b, dom, blockID, &!{{{
#include "dim_dummy_args.inc"
                            )

      implicit none

      integer, intent(in) :: nHaloLayers !< Input: Number of halo laters
      type (block_type), pointer :: b !< Input/Output: Block structure
      type (domain_type), pointer :: dom !< Input: Domain structure
      integer, intent(in) :: blockID !< Input: Global ID of block
#include "dim_dummy_decls.inc"



      integer :: i

      b % blockID = blockID

      allocate(b % parinfo)

      b % domain => dom

#include "block_allocs.inc"

   end subroutine mpas_allocate_block!}}}


#include "group_alloc_routines.inc"

#include "provis_alloc_routines.inc"

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

      block_ptr => dom % blocklist
      do while (associated(block_ptr))
         call mpas_deallocate_block(block_ptr)
         block_ptr => block_ptr % next
      end do

   end subroutine mpas_deallocate_domain!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field1d_integer
!
!> \brief   MPAS 1D Scratch integer allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 1D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field1d_integer(f, single_block_in)!{{{
       type (field1dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated or all blocks.
       logical :: single_block
       type (field1dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_integer!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field2d_integer
!
!> \brief   MPAS 2D Scratch integer allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 2D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field2d_integer(f, single_block_in)!{{{
       type (field2dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field2dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field2d_integer!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field3d_integer
!
!> \brief   MPAS 3D Scratch integer allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 3D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field3d_integer(f, single_block_in)!{{{
       type (field3dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field3dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field3d_integer!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field1d_real
!
!> \brief   MPAS 1D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 1D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field1d_real(f, single_block_in)!{{{
       type (field1dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field1dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field2d_real
!
!> \brief   MPAS 2D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 2D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field2d_real(f, single_block_in)!{{{
       type (field2dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field2dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field2d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field3d_real
!
!> \brief   MPAS 3D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 3D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field3d_real(f, single_block_in)!{{{
       type (field3dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field3dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field3d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field4D_real
!
!> \brief   MPAS 4D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 4D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field4d_real(f, single_block_in)!{{{
       type (field4dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field4dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3), f_cursor % dimSizes(4)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3), f % dimSizes(4)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field4d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field5D_real
!
!> \brief   MPAS 5D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 5D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field5d_real(f, single_block_in)!{{{
       type (field5dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field5dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3), f_cursor % dimSizes(4), f_cursor % dimSizes(5)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3), f % dimSizes(4), f % dimSizes(5)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field5d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field1D_char
!
!> \brief   MPAS 1D Scratch character deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 1D scratch character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field1d_char(f, single_block_in)!{{{
       type (field1dChar), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field1dChar), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_char!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field1D_integer
!
!> \brief   MPAS 1D Scratch integer deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field1d_integer(f, single_block_in)!{{{
       type (field1dInteger), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field1dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field2D_integer
!
!> \brief   MPAS 2D Scratch integer deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field2d_integer(f, single_block_in)!{{{
       type (field2dInteger), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field2dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field2d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field3D_integer
!
!> \brief   MPAS 3D Scratch integer deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field3d_integer(f, single_block_in)!{{{
       type (field3dInteger), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field3dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field3d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field1D_real
!
!> \brief   MPAS 1D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field1d_real(f, single_block_in)!{{{
       type (field1dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field1dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field2D_real
!
!> \brief   MPAS 2D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field2d_real(f, single_block_in)!{{{
       type (field2dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field2dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field2d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field3D_real
!
!> \brief   MPAS 3D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field3d_real(f, single_block_in)!{{{
       type (field3dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field3dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field3d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field4D_real
!
!> \brief   MPAS 4D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 4D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field4d_real(f, single_block_in)!{{{
       type (field4dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field4dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field4d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field5D_real
!
!> \brief   MPAS 5D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 5D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field5d_real(f, single_block_in)!{{{
       type (field5dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field5dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field5d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field1D_char
!
!> \brief   MPAS 1D Scratch character deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D scratch character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field1d_char(f, single_block_in)!{{{
       type (field1dChar), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field1dChar), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_char!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field0d_integer
!
!> \brief   MPAS 0D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 0D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field0d_integer(f)!{{{
       type (field0dInteger), pointer :: f !< Input: Field to deallocate
       type (field0dInteger), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)
         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field1D_integer
!
!> \brief   MPAS 1D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field1d_integer(f)!{{{
       type (field1dInteger), pointer :: f !< Input: Field to deallocate
       type (field1dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field2D_integer
!
!> \brief   MPAS 2D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field2d_integer(f)!{{{
       type (field2dInteger), pointer :: f !< Input: Field to deallocate
       type (field2dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field2d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field3D_integer
!
!> \brief   MPAS 3D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field3d_integer(f)!{{{
       type (field3dInteger), pointer :: f !< Input: Field to deallocate
       type (field3dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field3d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field0d_real
!
!> \brief   MPAS 0D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 0D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field0d_real(f)!{{{
       type (field0dReal), pointer :: f !< Input: Field to deallocate
       type (field0dReal), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field1D_real
!
!> \brief   MPAS 1D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field1d_real(f)!{{{
       type (field1dReal), pointer :: f !< Input: Field to deallocate
       type (field1dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field2D_real
!
!> \brief   MPAS 2D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D real field.
!

   subroutine mpas_deallocate_field2d_real(f)!{{{
       type (field2dReal), pointer :: f !< Input: Field to deallocate
       type (field2dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field2d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field3D_real
!
!> \brief   MPAS 3D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field3d_real(f)!{{{
       type (field3dReal), pointer :: f !< Input: Field to deallocate
       type (field3dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field3d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field4D_real
!
!> \brief   MPAS 4D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 4D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field4d_real(f)!{{{
       type (field4dReal), pointer :: f !< Input: Field to deallocate
       type (field4dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field4d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field5D_real
!
!> \brief   MPAS 5D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 5D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field5d_real(f)!{{{
       type (field5dReal), pointer :: f !< Input: Field to deallocate
       type (field5dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field5d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field0D_char
!
!> \brief   MPAS 0D character deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 0D character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field0d_char(f)!{{{
       type (field0dChar), pointer :: f !< Input: Field to deallocate
       type (field0dChar), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)
         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_char!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field1D_char
!
!> \brief   MPAS 1D character deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field1d_char(f)!{{{
       type (field1dChar), pointer :: f !< Input: Field to deallocate
       type (field1dChar), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_char!}}}

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

#include "block_deallocs.inc"

   end subroutine mpas_deallocate_block!}}}


#include "group_dealloc_routines.inc"


#include "group_copy_routines.inc"


#include "group_shift_level_routines.inc"


#include "field_links.inc"


end module mpas_grid_types
