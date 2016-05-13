! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_field_routines
!
!> \brief   MPAS Field Routines
!> \author  Doug Jacobsen, Michael Duda
!> \date    03/10/2015
!> \details 
!> This module defines routines related to MPAS field types (defined in mpas_data_types).
!> All routines defined in this module are non-blocking for threads.
!
!-----------------------------------------------------------------------
module mpas_field_routines

   use mpas_kind_types
   use mpas_derived_types
   use mpas_threading
   use mpas_attlist

   interface mpas_allocate_mold
      module procedure mpas_allocate_mold_1dreal
      module procedure mpas_allocate_mold_2dreal
      module procedure mpas_allocate_mold_3dreal
      module procedure mpas_allocate_mold_4dreal
      module procedure mpas_allocate_mold_5dreal
      module procedure mpas_allocate_mold_1dinteger
      module procedure mpas_allocate_mold_2dinteger
      module procedure mpas_allocate_mold_3dinteger
      module procedure mpas_allocate_mold_1dchar
   end interface

   interface mpas_duplicate_field
      module procedure mpas_duplicate_field0d_real
      module procedure mpas_duplicate_field1d_real
      module procedure mpas_duplicate_field2d_real
      module procedure mpas_duplicate_field3d_real
      module procedure mpas_duplicate_field4d_real
      module procedure mpas_duplicate_field5d_real
      module procedure mpas_duplicate_field0d_integer
      module procedure mpas_duplicate_field1d_integer
      module procedure mpas_duplicate_field2d_integer
      module procedure mpas_duplicate_field3d_integer
      module procedure mpas_duplicate_field0d_char
      module procedure mpas_duplicate_field1d_char
      module procedure mpas_duplicate_field0d_logical
   end interface

   interface mpas_shift_time_levs
      module procedure mpas_shift_time_levs_0dreal
      module procedure mpas_shift_time_levs_1dreal
      module procedure mpas_shift_time_levs_2dreal
      module procedure mpas_shift_time_levs_3dreal
      module procedure mpas_shift_time_levs_4dreal
      module procedure mpas_shift_time_levs_5dreal
      module procedure mpas_shift_time_levs_0dinteger
      module procedure mpas_shift_time_levs_1dinteger
      module procedure mpas_shift_time_levs_2dinteger
      module procedure mpas_shift_time_levs_3dinteger
      module procedure mpas_shift_time_levs_0dchar
      module procedure mpas_shift_time_levs_1dchar
      module procedure mpas_shift_time_levs_0dlogical
   end interface

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
      module procedure mpas_deallocate_field0d_logical
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
!  routine mpas_allocate_scratch_field1d_integer
!
!> \brief   MPAS 1D Scratch integer allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 1D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field1d_integer(f, single_block_in, init_array_in)!{{{
       type (field1dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field1dInteger), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1)))
                   if ( init_array ) then
                      f_cursor % array(:) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1)))
               if ( init_array ) then
                  f % array(:) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field2d_integer(f, single_block_in, init_array_in)!{{{
       type (field2dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field2dInteger), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
                   if ( init_array ) then
                      f_cursor % array(:, :) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
               if ( init_array ) then
                  f % array(:, :) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field3d_integer(f, single_block_in, init_array_in)!{{{
       type (field3dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field3dInteger), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
                   if ( init_array ) then
                      f_cursor % array(:, :, :) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
               if ( init_array ) then
                  f % array(:, :, :) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field1d_real(f, single_block_in, init_array_in)!{{{
       type (field1dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field1dReal), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1)))
                   if ( init_array ) then
                      f_cursor % array(:) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1)))
               if ( init_array ) then
                  f % array(:) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field2d_real(f, single_block_in, init_array_in)!{{{
       type (field2dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field2dReal), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
                   if ( init_array ) then
                      f_cursor % array(:, :) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
               if ( init_array ) then
                  f % array(:, :) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field3d_real(f, single_block_in, init_array_in)!{{{
       type (field3dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field3dReal), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
                   if ( init_array ) then
                      f_cursor % array(:, :, :) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
               if ( init_array ) then
                  f % array(:, :, :) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field4d_real(f, single_block_in, init_array_in)!{{{
       type (field4dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field4dReal), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3), f_cursor % dimSizes(4)))
                   if ( init_array ) then
                      f_cursor % array(:, :, :, :) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3), f % dimSizes(4)))
               if ( init_array ) then
                  f % array(:, :, :, :) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field5d_real(f, single_block_in, init_array_in)!{{{
       type (field5dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field5dReal), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3), f_cursor % dimSizes(4), f_cursor % dimSizes(5)))
                   if ( init_array ) then
                      f_cursor % array(:, :, :, :, :) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3), f % dimSizes(4), f % dimSizes(5)))
               if ( init_array ) then
                  f % array(:, :, :, :, :) = f % defaultValue
               end if
             end if
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
   subroutine mpas_allocate_scratch_field1d_char(f, single_block_in, init_array_in)!{{{
       type (field1dChar), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical, intent(in), optional :: init_array_in !< Input: Logical flag that determines if allocated arrays are initialized
       logical :: single_block, init_array
       type (field1dChar), pointer :: f_cursor
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if (present(init_array_in)) then
          init_array = init_array_in
       else
          init_array = .true.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          if(.not. single_block) then
             f_cursor => f
             do while(associated(f_cursor))
                if(.not.associated(f_cursor % array)) then
                   allocate(f_cursor % array(f_cursor % dimSizes(1)))
                   if ( init_array ) then
                      f_cursor % array(:) = f_cursor % defaultValue
                   end if
                end if
                f_cursor => f_cursor % next
             end do
          else
             if(.not.associated(f % array)) then
               allocate(f % array(f % dimSizes(1)))
               if ( init_array ) then
                  f % array(:) = f % defaultValue
               end if
             end if
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       integer :: threadNum

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
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
       end if

   end subroutine mpas_deallocate_scratch_field1d_char!}}}


!***********************************************************************
!
!  routine mpas_deallocate_field0d_logical
!
!> \brief   MPAS 0D logical deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 0D logical field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field0d_logical(f)!{{{
       type (field0dLogical), pointer :: f !< Input: Field to deallocate
       type (field0dLogical), pointer :: f_cursor
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()
   
       if ( threadNum == 0 ) then
          f_cursor => f
   
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if
   
            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)
            f_cursor => f
          end do
       end if

   end subroutine mpas_deallocate_field0d_logical!}}}


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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()
   
       if ( threadNum == 0 ) then
          f_cursor => f
   
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if
   
            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)
            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if
   
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)
   
            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if
   
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)
   
            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if
   
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)
   
            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       f_cursor => f

       if ( threadNum == 0 ) then
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if
   
            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)
   
            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if

            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if

            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)

            f_cursor => f
          end do
       end if

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
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field2d_real(f)!{{{
       type (field2dReal), pointer :: f !< Input: Field to deallocate
       type (field2dReal), pointer :: f_cursor
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if

            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if

            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)

            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if

            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if

            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)

            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if

            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if

            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)

            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if

            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if

            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)

            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f

          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if

            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)
            f_cursor => f
          end do
       end if

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
       integer :: threadNum
       integer :: i, iErr

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f % next)) then
              f => f % next
            else
              nullify(f)
            end if

            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if

            if ( associated(f_cursor % attLists) ) then
               do i = 1, size(f_cursor % attLists, dim=1)
                  call mpas_deallocate_attlist(f_cursor % attLists(i) % attList, iErr)
               end do
               deallocate(f_cursor % attLists)
            end if

            deallocate(f_cursor)

            f_cursor => f
          end do
       end if

   end subroutine mpas_deallocate_field1d_char!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_1dreal
!
!> \brief   Allocates a 1-d real array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_1dreal(dst, src)!{{{

      implicit none

      real(kind=RKIND), dimension(:), pointer :: dst
      real(kind=RKIND), dimension(:), intent(in) :: src
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()
  
      if ( threadNum == 0 ) then
         allocate(dst(size(src)))
      end if

   end subroutine mpas_allocate_mold_1dreal!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_2dreal
!
!> \brief   Allocates a 2-d real array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_2dreal(dst, src)!{{{

      implicit none

      real(kind=RKIND), dimension(:,:), pointer :: dst
      real(kind=RKIND), dimension(:,:), intent(in) :: src
  
      integer, dimension(2) :: dims
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      dims = shape(src)

      if ( threadNum == 0 ) then
         allocate(dst(dims(1),dims(2)))
      end if

   end subroutine mpas_allocate_mold_2dreal!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_3dreal
!
!> \brief   Allocates a 3-d real array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_3dreal(dst, src)!{{{

      implicit none

      real(kind=RKIND), dimension(:,:,:), pointer :: dst
      real(kind=RKIND), dimension(:,:,:), intent(in) :: src
  
      integer, dimension(3) :: dims
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      dims = shape(src)

      if ( threadNum == 0 ) then
         allocate(dst(dims(1),dims(2),dims(3)))
      end if

   end subroutine mpas_allocate_mold_3dreal!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_4dreal
!
!> \brief   Allocates a 4-d real array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_4dreal(dst, src)!{{{

      implicit none

      real(kind=RKIND), dimension(:,:,:,:), pointer :: dst
      real(kind=RKIND), dimension(:,:,:,:), intent(in) :: src
  
      integer, dimension(4) :: dims
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      dims = shape(src)

      if ( threadNum == 0 ) then
         allocate(dst(dims(1),dims(2),dims(3),dims(4)))
      end if

   end subroutine mpas_allocate_mold_4dreal!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_5dreal
!
!> \brief   Allocates a 5-d real array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_5dreal(dst, src)!{{{

      implicit none

      real(kind=RKIND), dimension(:,:,:,:,:), pointer :: dst
      real(kind=RKIND), dimension(:,:,:,:,:), intent(in) :: src
  
      integer, dimension(5) :: dims
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      dims = shape(src)

      if ( threadNum == 0 ) then
         allocate(dst(dims(1),dims(2),dims(3),dims(4),dims(5)))
      end if

   end subroutine mpas_allocate_mold_5dreal!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_1dinteger
!
!> \brief   Allocates a 1-d integer array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_1dinteger(dst, src)!{{{

      implicit none

      integer, dimension(:), pointer :: dst
      integer, dimension(:), intent(in) :: src
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()
  
      if ( threadNum == 0 ) then
         allocate(dst(size(src)))
      end if

   end subroutine mpas_allocate_mold_1dinteger!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_2dinteger
!
!> \brief   Allocates a 2-d integer array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_2dinteger(dst, src)!{{{

      implicit none

      integer, dimension(:,:), pointer :: dst
      integer, dimension(:,:), intent(in) :: src
  
      integer, dimension(2) :: dims
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      dims = shape(src)

      if ( threadNum == 0 ) then
         allocate(dst(dims(1),dims(2)))
      end if

   end subroutine mpas_allocate_mold_2dinteger!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_3dinteger
!
!> \brief   Allocates a 3-d integer array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_3dinteger(dst, src)!{{{

      implicit none

      integer, dimension(:,:,:), pointer :: dst
      integer, dimension(:,:,:), intent(in) :: src
  
      integer, dimension(3) :: dims
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      dims = shape(src)

      if ( threadNum == 0 ) then
         allocate(dst(dims(1),dims(2),dims(3)))
      end if

   end subroutine mpas_allocate_mold_3dinteger!}}}


!***********************************************************************
!
!  routine mpas_allocate_mold_1dchar
!
!> \brief   Allocates a 1-d character array using the dimensions of another array
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Allocates the array dst to have the same dimensions as the array src.
!> This routine exists to provide the same functionality as F2008's 
!> ALLOCATE(A,MOLD=B) functionality, or a similar functionality to F2003's
!> ALLOCATE(A,SOURCE=B) but without actually copying the source values.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_mold_1dchar(dst, src)!{{{

      implicit none

      character(len=StrKIND), dimension(:), pointer :: dst
      character(len=StrKIND), dimension(:), intent(in) :: src
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()
  
      if ( threadNum == 0 ) then
         allocate(dst(size(src)))
      end if

   end subroutine mpas_allocate_mold_1dchar!}}}

!***********************************************************************
!
!  routine mpas_duplicate_field0d_real
!
!> \brief   MPAS 0D real field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field0d_real(src, dst, copy_array_only) !{{{

      implicit none

      type (field0DReal), intent(in), target :: src     !< Input: Field to be duplicated
      type (field0DReal), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field0DReal), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_scalar.inc"

   end subroutine mpas_duplicate_field0d_real !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field1d_real
!
!> \brief   MPAS 1D real field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field1d_real(src, dst, copy_array_only) !{{{

      implicit none

      type (field1DReal), intent(in), target :: src     !< Input: Field to be duplicated
      type (field1DReal), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field1DReal), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field1d_real !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field2d_real
!
!> \brief   MPAS 2D real field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field2d_real(src, dst, copy_array_only) !{{{

      implicit none

      type (field2DReal), intent(in), target :: src     !< Input: Field to be duplicated
      type (field2DReal), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field2DReal), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field2d_real !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field3d_real
!
!> \brief   MPAS 3D real field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field3d_real(src, dst, copy_array_only) !{{{

      implicit none

      type (field3DReal), intent(in), target :: src     !< Input: Field to be duplicated
      type (field3DReal), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field3DReal), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field3d_real !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field4d_real
!
!> \brief   MPAS 4D real field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field4d_real(src, dst, copy_array_only) !{{{

      implicit none

      type (field4DReal), intent(in), target :: src     !< Input: Field to be duplicated
      type (field4DReal), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field4DReal), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field4d_real !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field5d_real
!
!> \brief   MPAS 5D real field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field5d_real(src, dst, copy_array_only) !{{{

      implicit none

      type (field5DReal), intent(in), target :: src     !< Input: Field to be duplicated
      type (field5DReal), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field5DReal), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field5d_real !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field0d_integer
!
!> \brief   MPAS 0D integer field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field0d_integer(src, dst, copy_array_only) !{{{

      implicit none

      type (field0DInteger), intent(in), target :: src  !< Input: Field to be duplicated
      type (field0DInteger), pointer :: dst             !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field0DInteger), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_scalar.inc"

   end subroutine mpas_duplicate_field0d_integer !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field1d_integer
!
!> \brief   MPAS 1D integer field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field1d_integer(src, dst, copy_array_only) !{{{

      implicit none

      type (field1DInteger), intent(in), target :: src  !< Input: Field to be duplicated
      type (field1DInteger), pointer :: dst             !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field1DInteger), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field1d_integer !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field2d_integer
!
!> \brief   MPAS 2D integer field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field2d_integer(src, dst, copy_array_only) !{{{

      implicit none

      type (field2DInteger), intent(in), target :: src  !< Input: Field to be duplicated
      type (field2DInteger), pointer :: dst             !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field2DInteger), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field2d_integer !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field3d_integer
!
!> \brief   MPAS 3D integer field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field3d_integer(src, dst, copy_array_only) !{{{

      implicit none

      type (field3DInteger), intent(in), target :: src  !< Input: Field to be duplicated
      type (field3DInteger), pointer :: dst             !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field3DInteger), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field3d_integer !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field0d_char
!
!> \brief   MPAS 0D character field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field0d_char(src, dst, copy_array_only) !{{{

      implicit none

      type (field0DChar), intent(in), target :: src     !< Input: Field to be duplicated
      type (field0DChar), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field0DChar), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_scalar.inc"

   end subroutine mpas_duplicate_field0d_char !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field1d_char
!
!> \brief   MPAS 1D character field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field1d_char(src, dst, copy_array_only) !{{{

      implicit none

      type (field1DChar), intent(in), target :: src     !< Input: Field to be duplicated
      type (field1DChar), pointer :: dst                !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field1DChar), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_array.inc"

   end subroutine mpas_duplicate_field1d_char !}}}


!***********************************************************************
!
!  routine mpas_duplicate_field0d_logical
!
!> \brief   MPAS 0D logical field duplication routine.
!> \author  Michael Duda
!> \date    04/12/14
!> \details 
!> Creates a duplicate of the source field.
!
!-----------------------------------------------------------------------
   subroutine mpas_duplicate_field0d_logical(src, dst, copy_array_only) !{{{

      implicit none

      type (field0DLogical), intent(in), target :: src  !< Input: Field to be duplicated
      type (field0DLogical), pointer :: dst             !< Output: Field to contain the duplicate
      logical, intent(in), optional :: copy_array_only  !< Input: whether to assume that dst exists, and only copy array data

      type (field0DLogical), pointer :: src_cursor, dst_cursor
      logical :: local_copy_only
      integer :: threadNum, iConstituent

#include "duplicate_field_scalar.inc"

   end subroutine mpas_duplicate_field0d_logical !}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_0dreal
!
!> \brief   MPAS 0D real time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_0dreal(fldarr)!{{{

      implicit none

      type (field0DReal), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field0DReal), dimension(:), pointer :: fldarr_ptr
      real(kind=RKIND) :: scalar
      
#include "shift_time_levs_scalar.inc"

   end subroutine mpas_shift_time_levs_0dreal!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_1dreal
!
!> \brief   MPAS 1D real time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_1dreal(fldarr)!{{{

      implicit none

      type (field1DReal), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field1DReal), dimension(:), pointer :: fldarr_ptr
      real(kind=RKIND), dimension(:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_1dreal!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_2dreal
!
!> \brief   MPAS 2D real time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_2dreal(fldarr)!{{{

      implicit none

      type (field2DReal), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field2DReal), dimension(:), pointer :: fldarr_ptr
      real(kind=RKIND), dimension(:,:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_2dreal!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_3dreal
!
!> \brief   MPAS 3D real time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_3dreal(fldarr)!{{{

      implicit none

      type (field3DReal), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field3DReal), dimension(:), pointer :: fldarr_ptr
      real(kind=RKIND), dimension(:,:,:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_3dreal!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_4dreal
!
!> \brief   MPAS 4D real time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_4dreal(fldarr)!{{{

      implicit none

      type (field4DReal), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field4DReal), dimension(:), pointer :: fldarr_ptr
      real(kind=RKIND), dimension(:,:,:,:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_4dreal!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_5dreal
!
!> \brief   MPAS 5D real time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_5dreal(fldarr)!{{{

      implicit none

      type (field5DReal), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field5DReal), dimension(:), pointer :: fldarr_ptr
      real(kind=RKIND), dimension(:,:,:,:,:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_5dreal!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_0dinteger
!
!> \brief   MPAS 0D integer time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_0dinteger(fldarr)!{{{

      implicit none

      type (field0DInteger), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field0DInteger), dimension(:), pointer :: fldarr_ptr
      integer :: scalar
      
#include "shift_time_levs_scalar.inc"

   end subroutine mpas_shift_time_levs_0dinteger!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_1dinteger
!
!> \brief   MPAS 1D integer time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_1dinteger(fldarr)!{{{

      implicit none

      type (field1DInteger), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field1DInteger), dimension(:), pointer :: fldarr_ptr
      integer, dimension(:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_1dinteger!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_2dinteger
!
!> \brief   MPAS 2D integer time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_2dinteger(fldarr)!{{{

      implicit none

      type (field2DInteger), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field2DInteger), dimension(:), pointer :: fldarr_ptr
      integer, dimension(:,:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_2dinteger!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_3dinteger
!
!> \brief   MPAS 3D integer time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_3dinteger(fldarr)!{{{

      implicit none

      type (field3DInteger), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field3DInteger), dimension(:), pointer :: fldarr_ptr
      integer, dimension(:,:,:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_3dinteger!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_0dchar
!
!> \brief   MPAS 0D character time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_0dchar(fldarr)!{{{

      implicit none

      type (field0DChar), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field0DChar), dimension(:), pointer :: fldarr_ptr
      character (len=StrKIND) :: scalar
      
#include "shift_time_levs_scalar.inc"

   end subroutine mpas_shift_time_levs_0dchar!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_1dchar
!
!> \brief   MPAS 1D character time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_1dchar(fldarr)!{{{

      implicit none

      type (field1DChar), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field1DChar), dimension(:), pointer :: fldarr_ptr
      character (len=StrKIND), dimension(:), pointer :: arr_ptr
      
#include "shift_time_levs_array.inc"

   end subroutine mpas_shift_time_levs_1dchar!}}}


!***********************************************************************
!
!  routine mpas_shift_time_levs_0dlogical
!
!> \brief   MPAS 0D logical time-level shift routine
!> \author  Michael Duda
!> \date    04/14/14
!> \details 
!> Shifts the contents of the array of fields provided by the input argument.
!> After returning, the storage for fldarr(n) will point to what was the storage 
!> for fldarr(n+1) in a period fashion, so that, for N time levels, the storage
!> for fldarr(N) will point to what was the storage for fldarr(1).
!
!-----------------------------------------------------------------------
   subroutine mpas_shift_time_levs_0dlogical(fldarr)!{{{

      implicit none

      type (field0DLogical), dimension(:), pointer :: fldarr

      integer :: i, nlevs, threadNum
      type (field0DLogical), dimension(:), pointer :: fldarr_ptr
      logical :: scalar
      
#include "shift_time_levs_scalar.inc"

   end subroutine mpas_shift_time_levs_0dlogical!}}}


end module mpas_field_routines
