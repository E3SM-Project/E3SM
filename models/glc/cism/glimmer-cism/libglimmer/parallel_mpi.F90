!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   parallel_mpi.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module parallel

  use netcdf
  implicit none

!PW - Repeat from glimmer_horiz_bcs_parallel.F90; needs to be set some place
!     which can be 'used' by both the parallel and glimmer_horiz_bcs modules
  integer, parameter, private :: HORIZ_BCS_WALL_SLIP = 0
  integer, parameter, private :: HORIZ_BCS_CYCLIC = 1

  integer, parameter, private :: horiz_bcs_type_north = HORIZ_BCS_CYCLIC
  integer, parameter, private :: horiz_bcs_type_south = HORIZ_BCS_CYCLIC
  integer, parameter, private :: horiz_bcs_type_east  = HORIZ_BCS_CYCLIC
  integer, parameter, private :: horiz_bcs_type_west  = HORIZ_BCS_CYCLIC
!PW - End of repeat

  ! Debug and Verification Level
  integer,parameter :: DEBUG_LEVEL = 1 
	! If > 0, then debug code executed.  Added for parallel_halo_verify()

  !NOTE: The glam/glissade dycore currently requires nhalo = 2,
  !       whereas the glide dycore requires nhalo = 0.
  !      For glide simulations, we set nhalo = 0 by calling distributed_grid 
  !       with optional argument nhalo = 0.

  integer, save :: nhalo = 2

  !TODO - If we will always have lhalo = uhalo = nhalo, then we should 
  !       define lhalo and uhalo in terms of nhalo.

  integer, save :: lhalo = 2
  integer, save :: uhalo = 2

  ! halo widths for staggered grid
!  integer,parameter :: staggered_lhalo = lhalo
!  integer,parameter :: staggered_uhalo = uhalo-1
  integer, save :: staggered_lhalo = 2
  integer, save :: staggered_uhalo = 1

!TODO - Remove staggered_whalo/shalo/ehalo/nhalo here and in other parts of the code
!  integer,parameter :: staggered_whalo = lhalo
!  integer,parameter :: staggered_shalo = lhalo
!  integer,parameter :: staggered_ehalo = uhalo-1
!  integer,parameter :: staggered_nhalo = uhalo-1
  integer, save :: staggered_whalo = 2
  integer, save :: staggered_shalo = 2
  integer, save :: staggered_ehalo = 1
  integer, save :: staggered_nhalo = 1

  integer,save :: main_rank
  logical,save :: main_task
  integer,save :: comm, tasks, this_rank

  ! distributed grid
  integer,save :: global_ewn,global_nsn,local_ewn,local_nsn,own_ewn,own_nsn
  integer,save :: global_col_offset, global_row_offset

  integer,save :: ewlb,ewub,nslb,nsub
  integer,save :: east,north,south,west

  ! common work space
  integer,dimension(4),save :: d_gs_mybounds
  integer,dimension(:,:),allocatable,save :: d_gs_bounds

  ! distributed gather flow control parameter
  integer,parameter :: max_gather_block_size = 64 ! max and default

  ! global IDs
  integer,save :: ProcsEW

  !TODO - Remove these declarations

  ! JEFF Declarations for undistributed variables on main_task.
  ! Later move to separate module?  These are only temporary until code is completely distributed.
  real(8),dimension(:,:,:),allocatable :: gathered_efvs  ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:,:),allocatable :: gathered_efvs2  ! Variable for testing that scatter/gather are inverses
  real(8),dimension(:,:,:),allocatable :: gathered_uvel  ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:,:),allocatable :: gathered_vvel  ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:),allocatable :: gathered_uflx    ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:),allocatable :: gathered_vflx    ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:,:),allocatable :: gathered_velnorm  ! Variable calculated in run_ho_diagnostic(), is this used?
  real(8),dimension(:,:),allocatable :: gathered_thck    ! Used in horizontal_remap_in()
  real(8),dimension(:,:),allocatable :: gathered_stagthck ! Used in horizontal_remap_in()
  real(4),dimension(:,:),allocatable :: gathered_acab    ! Used in horizontal_remap_in()
  real(8),dimension(:,:,:),allocatable :: gathered_temp  ! Used in horizontal_remap_in()
  real(8),dimension(:,:),allocatable :: gathered_dusrfdew  ! Used in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_dusrfdns  ! Used in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_dthckdew  ! Used in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_dthckdns  ! Used in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauxx   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauyy   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauxy   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauscalar   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauxz   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauyz   ! Calculated in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_topg  ! Bedrock topology, Used in glide_set_mask()
  integer,dimension(:,:),allocatable :: gathered_thkmask  ! Calculated in glide_set_mask()
  real(8),dimension(:,:),allocatable :: gathered_marine_bc_normal  ! Calculated in glide_marine_margin_normal()
  real(8),dimension(:,:,:),allocatable :: gathered_surfvel   ! Used in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_gline_flux   ! Calculated in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_ubas   ! Used in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_vbas   ! Used in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_relx   ! Used in glide_marinlim()
  real(8),dimension(:,:,:),allocatable :: gathered_flwa   ! Used in glide_marinlim()
  real(4),dimension(:,:),allocatable :: gathered_calving   ! Used in glide_marinlim()
  real(4),dimension(:,:),allocatable :: gathered_backstress   ! Used in glide_marinlim()
  real(8),dimension(:,:),allocatable :: gathered_usrf   ! Used in glide_marinlim()
  logical,dimension(:,:),allocatable :: gathered_backstressmap ! Used in glide_marinlim()
  real(8),dimension(:,:),allocatable :: gathered_tau_x   ! Calculated in calc_basal_shear()
  real(8),dimension(:,:),allocatable :: gathered_tau_y   ! Calculated in calc_basal_shear()
  real(8),dimension(:,:),allocatable :: gathered_lsrf   ! Used in glide_marinlim()

  interface broadcast
     module procedure broadcast_character
     module procedure broadcast_integer
     module procedure broadcast_integer_1d
     module procedure broadcast_logical
     module procedure broadcast_real4
     module procedure broadcast_real4_1d
     module procedure broadcast_real8     
     module procedure broadcast_real8_1d
  end interface

  interface distributed_gather_var
     module procedure distributed_gather_var_integer_2d
     module procedure distributed_gather_var_logical_2d
     module procedure distributed_gather_var_real4_2d
     module procedure distributed_gather_var_real4_3d
     module procedure distributed_gather_var_real8_2d
     module procedure distributed_gather_var_real8_3d
  end interface

  interface distributed_get_var
     module procedure distributed_get_var_integer_2d
     module procedure distributed_get_var_real4_1d
     module procedure distributed_get_var_real4_2d
     module procedure distributed_get_var_real8_2d
     module procedure distributed_get_var_real8_3d
  end interface

  interface distributed_print
     ! Gathers a distributed variable and writes to file
     module procedure distributed_print_integer_2d
     module procedure distributed_print_real8_2d
     module procedure distributed_print_real8_3d
  end interface

  interface distributed_put_var
     module procedure distributed_put_var_integer_2d
     module procedure distributed_put_var_real4_1d
     module procedure distributed_put_var_real4_2d
     module procedure distributed_put_var_real8_2d
     module procedure distributed_put_var_real8_3d

     !TODO - Should the parallel_put_var routines be part of this interface?
     module procedure parallel_put_var_real4
     module procedure parallel_put_var_real8
  end interface

  interface distributed_scatter_var
     module procedure distributed_scatter_var_integer_2d
     module procedure distributed_scatter_var_logical_2d
     module procedure distributed_scatter_var_real4_2d
     module procedure distributed_scatter_var_real4_3d
     module procedure distributed_scatter_var_real8_2d
     module procedure distributed_scatter_var_real8_3d
  end interface

  interface parallel_def_var
     module procedure parallel_def_var_dimids
     module procedure parallel_def_var_nodimids
  end interface

  interface parallel_get_att
     module procedure parallel_get_att_character
     module procedure parallel_get_att_real4
     module procedure parallel_get_att_real4_1d
     module procedure parallel_get_att_real8
     module procedure parallel_get_att_real8_1d
  end interface

  interface parallel_get_var
     module procedure parallel_get_var_integer_1d
     module procedure parallel_get_var_real4_1d
     module procedure parallel_get_var_real8_1d
  end interface

  interface parallel_halo
     module procedure parallel_halo_integer_2d
     module procedure parallel_halo_logical_2d
     module procedure parallel_halo_real4_2d
     module procedure parallel_halo_real8_2d
     module procedure parallel_halo_real8_3d
  end interface

  interface parallel_halo_verify
     module procedure parallel_halo_verify_integer_2d
     module procedure parallel_halo_verify_real8_2d
     module procedure parallel_halo_verify_real8_3d
  end interface

  interface staggered_parallel_halo
     module procedure staggered_parallel_halo_integer_2d
     module procedure staggered_parallel_halo_integer_3d
     module procedure staggered_parallel_halo_real8_2d
     module procedure staggered_parallel_halo_real8_3d
     module procedure staggered_parallel_halo_real8_6d
  end interface

  interface parallel_print
     ! Writes a parallel (same on all processors) variable to file by just writing from main_task
     module procedure parallel_print_integer_2d
     module procedure parallel_print_real8_2d
     module procedure parallel_print_real8_3d
  end interface

  interface parallel_put_att
     module procedure parallel_put_att_character
     module procedure parallel_put_att_real4
     module procedure parallel_put_att_real4_1d
     module procedure parallel_put_att_real8
     module procedure parallel_put_att_real8_1d
  end interface

  interface parallel_put_var
     module procedure parallel_put_var_real4
     module procedure parallel_put_var_real8
     module procedure parallel_put_var_real8_1d
  end interface

  interface parallel_reduce_max
     module procedure parallel_reduce_max_integer
     module procedure parallel_reduce_max_real4
     module procedure parallel_reduce_max_real8
  end interface

contains

  subroutine broadcast_character(c)
    use mpi_mod
    implicit none
    character(len=*) :: c
    integer :: ierror,n
    ! begin
    n = len(c)
    call mpi_bcast(c,n,mpi_character,main_rank,comm,ierror)
  end subroutine broadcast_character

  subroutine broadcast_integer(i)
    use mpi_mod
    implicit none
    integer :: i,ierror
    ! begin
    call mpi_bcast(i,1,mpi_integer,main_rank,comm,ierror)
  end subroutine broadcast_integer

  subroutine broadcast_integer_1d(a)
    use mpi_mod
    implicit none
    integer,dimension(:) :: a
    integer :: ierror
    ! begin
    call mpi_bcast(a,size(a),mpi_integer,main_rank,comm,ierror)
  end subroutine broadcast_integer_1d

  subroutine broadcast_logical(l)
    use mpi_mod
    implicit none
    logical :: l
    integer :: ierror
    ! begin
    call mpi_bcast(l,1,mpi_logical,main_rank,comm,ierror)
  end subroutine broadcast_logical

  subroutine broadcast_real4(r)
    use mpi_mod
    implicit none
    integer :: ierror
    real(4) :: r
    ! begin
    call mpi_bcast(r,1,mpi_real4,main_rank,comm,ierror)
  end subroutine broadcast_real4

  subroutine broadcast_real4_1d(a)
    use mpi_mod
    implicit none
    real(4),dimension(:) :: a
    integer :: ierror
    ! begin
    call mpi_bcast(a,size(a),mpi_real4,main_rank,comm,ierror)
  end subroutine broadcast_real4_1d

  subroutine broadcast_real8(r)
    use mpi_mod
    implicit none
    integer :: ierror
    real(8) :: r
    ! begin
    call mpi_bcast(r,1,mpi_real8,main_rank,comm,ierror)
  end subroutine broadcast_real8

  subroutine broadcast_real8_1d(a)
    use mpi_mod
    implicit none
    real(8),dimension(:) :: a
    integer :: ierror
    ! begin
    call mpi_bcast(a,size(a),mpi_real8,main_rank,comm,ierror)
  end subroutine broadcast_real8_1d

  function distributed_execution()
     ! Returns if running distributed or not.
     logical distributed_execution

     distributed_execution = .true.
  end function distributed_execution

  subroutine distributed_gather_var_integer_2d(values, global_values)

    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    integer,dimension(:,:),intent(in) :: values
    integer,dimension(:,:),allocatable,intent(inout) :: global_values

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: sendbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(&
                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1) &
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_int(sendbuf,size(sendbuf),mpi_integer,&
       recvbuf,recvcounts,displs,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if
    ! automatic deallocation
  end subroutine distributed_gather_var_integer_2d

  subroutine distributed_gather_var_logical_2d(values, global_values)

    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    logical,dimension(:,:),intent(in) :: values
    logical,dimension(:,:),allocatable,intent(inout) :: global_values

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    logical,dimension(:),allocatable :: recvbuf
    logical,dimension(:,:),allocatable :: sendbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(&
                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       global_values(:,:) = .false.
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_log(sendbuf,size(sendbuf),mpi_logical,&
         recvbuf,recvcounts,displs,mpi_logical,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if
    ! automatic deallocation
  end subroutine distributed_gather_var_logical_2d

  subroutine distributed_gather_var_real4_2d(values, global_values)

    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(4),dimension(:,:),intent(in) :: values
    real(4),dimension(:,:),allocatable,intent(inout) :: global_values

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(4),dimension(:),allocatable :: recvbuf
    real(4),dimension(:,:),allocatable :: sendbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(&
                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1) &
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real4(sendbuf,size(sendbuf),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if
    ! automatic deallocation
  end subroutine distributed_gather_var_real4_2d

  subroutine distributed_gather_var_real4_3d(values, global_values, ld1, ud1)

    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(4),dimension(:,:,:),intent(in) :: values
    real(4),dimension(:,:,:),allocatable,intent(inout) :: global_values
    integer,optional,intent(in) :: ld1, ud1

    integer :: i,ierror,j,k,d1l,d1u
    integer,dimension(:),allocatable :: displs,recvcounts
    real(4),dimension(:),allocatable :: recvbuf
    real(4),dimension(:,:,:),allocatable :: sendbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       if (present(ld1)) then
         d1l = ld1
       else
         d1l = 1
       endif
       if (present(ud1)) then
         d1u = ud1
       else
         d1u = size(values,1)-(d1l-1)
       endif
       if (size(values,1) /= d1u-d1l+1) then
          write(*,*) "size(values,1) .ne. d1u-d1l+1 in gather call"
          call parallel_stop(__FILE__, __LINE__)
       endif
       allocate(global_values(d1l:d1u,&
                              minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
                              minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       global_values(:,:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)&
                      *size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:,:) = values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real4(sendbuf,size(sendbuf),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(:,&
                        d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/size(values,1),&
                       d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if
    ! automatic deallocation
  end subroutine distributed_gather_var_real4_3d

  subroutine distributed_gather_var_real8_2d(values, global_values)

    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(8),dimension(:,:),intent(in) :: values
    real(8),dimension(:,:),allocatable,intent(inout) :: global_values

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:),allocatable :: sendbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(&
                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if
    ! automatic deallocation
  end subroutine distributed_gather_var_real8_2d

  subroutine distributed_gather_var_real8_3d(values, global_values, ld1, ud1)

    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(8),dimension(:,:,:),intent(in) :: values
    real(8),dimension(:,:,:),allocatable,intent(inout) :: global_values
    integer,optional,intent(in) :: ld1, ud1

    integer :: i,ierror,j,k,d1l,d1u
    integer,dimension(:),allocatable :: displs,recvcounts
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:,:),allocatable :: sendbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       if (present(ld1)) then
         d1l = ld1
       else
         d1l = 1
       endif
       if (present(ud1)) then
         d1u = ud1
       else
         d1u = size(values,1)-(d1l-1)
       endif
       if (size(values,1) /= d1u-d1l+1) then
          write(*,*) "size(values,1) .ne. d1u-d1l+1 in gather call"
          call parallel_stop(__FILE__, __LINE__)
       endif
       allocate(global_values(d1l:d1u,&
                              minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
                              minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       global_values(:,:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)&
                      *size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:,:) = values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(:,&
                        d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/size(values,1),&
                       d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if
    ! automatic deallocation
  end subroutine distributed_gather_var_real8_3d

  function distributed_get_var_integer_2d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_get_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: sendbuf
    integer,dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       distributed_get_var_integer_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_integer_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_integer,&
         recvbuf,size(recvbuf),mpi_integer,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))
    !automatic deallocation
  end function distributed_get_var_integer_2d

  function distributed_get_var_real4_1d(ncid,varid,values,start)
    use mpi_mod
    use netcdf
    implicit none
    integer :: distributed_get_var_real4_1d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:) :: values

    integer :: i,ierror,myn,status,x1id,y1id
    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: global_values,sendbuf

    ! begin

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       global_values(:) = 0
       distributed_get_var_real4_1d = &
            nf90_get_var(ncid,varid,global_values(1:myn),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
               global_values(bounds(1,i):bounds(2,i))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real4_1d)
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         values,size(values),mpi_real4,main_rank,comm,ierror)
    !automatic deallocation
  end function distributed_get_var_real4_1d

  function distributed_get_var_real4_2d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_get_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: sendbuf
    real(4),dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       distributed_get_var_real4_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real4_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         recvbuf,size(recvbuf),mpi_real4,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))
    !automatic deallocation
  end function distributed_get_var_real4_2d

  function distributed_get_var_real8_2d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_get_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: sendbuf
    real(8),dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       distributed_get_var_real8_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real8_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))
    !automatic deallocation

  end function distributed_get_var_real8_2d

  function distributed_get_var_real8_3d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_get_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: sendbuf
    real(8),dimension(:,:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:)),size(values,3)))
       global_values(:,:,:) = 0
       distributed_get_var_real8_3d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns,:),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*&
            (bounds(4,:)-bounds(3,:)+1)*size(values,3)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(global_values(&
               bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i),:),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real8_3d)
    allocate(recvbuf(local_ewn,local_nsn,size(values,3)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,:,:) = recvbuf(:size(values,1),:size(values,2),:)
    !automatic deallocation
  end function distributed_get_var_real8_3d

  function distributed_isparallel()
     implicit none
     logical :: distributed_isparallel

     distributed_isparallel = .true.
  end function distributed_isparallel


  subroutine distributed_grid(ewn, nsn, nhalo_in)

    implicit none
    integer, intent(inout) :: ewn, nsn        ! global grid dimensions
    integer, intent(in), optional :: nhalo_in ! number of rows of halo cells
    integer :: best,i,j,metric
    integer :: ewrank,ewtasks,nsrank,nstasks
    real(8) :: rewtasks,rnstasks

    ! begin

    ! Optionally, change the halo values
    ! Note: The higher-order dycores (glam, glissade) currently require nhalo = 2.
    !       The Glide SIA dycore requires nhalo = 0.
    ! The default halo values at the top of the module are appropriate for
    !  the higher-order dycores.  Here they can be reset to zero for Glide.

    if (present(nhalo_in)) then
       if (main_task) then
          write(*,*) 'Setting halo values: nhalo =', nhalo_in
          if (nhalo_in < 0) then
             write(*,*) 'ERROR: nhalo must be >= 0'
             call parallel_stop(__FILE__, __LINE__)
          elseif (nhalo_in /= 2) then
             write(*,*) 'WARNING: parallel dycores tested only with nhalo = 2'
          endif
       endif 
       nhalo = nhalo_in
       lhalo = nhalo
       uhalo = nhalo
       staggered_lhalo = lhalo
       staggered_uhalo = max(uhalo-1, 0)
       !TODO - Remove the following variables
       staggered_whalo = lhalo
       staggered_shalo = lhalo
       staggered_ehalo = max(uhalo-1, 0)
       staggered_nhalo = max(uhalo-1, 0)
    endif

    global_ewn = ewn
    global_nsn = nsn

    ewtasks = 0
    nstasks = 0
    best = huge(best)
    do i = 1,min(tasks,global_ewn)
       j = tasks/i
       if (j<=global_nsn.and.i*j==tasks) then ! try to use all tasks
          metric = abs(i*global_nsn-j*global_ewn) ! zero if ewn/nsn == i/j
          if (metric<best) then
             best = metric
             ewtasks = i
             nstasks = j
          end if
       end if
    end do
    if (ewtasks*nstasks/=tasks) call parallel_stop(__FILE__,__LINE__)

    ! Store critical value for creating global IDs.  Defines grid distribution.
    ProcsEW = ewtasks

    ! For globalID calculations determine processor's global grid index offsets
    ! sum block sizes for row blocks preceding this_rank
    ! Do not include halo offsets in global calculations
    ! (There are ProcsEW processors per row.)
    global_col_offset = 0
    do ewrank=0,mod(this_rank, ProcsEW)-1
      rewtasks = 1/real(ewtasks,8)
      ewlb = nint(ewrank*global_ewn*rewtasks)+1
      ewub = nint((ewrank+1)*global_ewn*rewtasks)
      own_ewn = ewub-ewlb+1
      global_col_offset = global_col_offset + own_ewn
    enddo

    ! sum block sizes for column blocks preceding this_rank
    ! (Integer division required for this_rank/ProcsEW)
    global_row_offset = 0
    do nsrank=0,(this_rank/ProcsEW)-1
      rnstasks = 1/real(nstasks,8)
      nslb = nint(nsrank*global_nsn*rnstasks)+1
      nsub = nint((nsrank+1)*global_nsn*rnstasks)
      own_nsn = nsub-nslb+1
      global_row_offset = global_row_offset + own_nsn
    enddo

    ! Set local processor's grid indices, including halo offsets
    ewrank = mod(this_rank,ewtasks)
    rewtasks = 1/real(ewtasks,8)
    ewlb = nint(ewrank*global_ewn*rewtasks)+1-lhalo
    ewub = nint((ewrank+1)*global_ewn*rewtasks)+uhalo
    local_ewn = ewub-ewlb+1
    own_ewn = local_ewn-lhalo-uhalo
    ewn = local_ewn

    nsrank = this_rank/ewtasks
    rnstasks = 1/real(nstasks,8)
    nslb = nint(nsrank*global_nsn*rnstasks)+1-lhalo
    nsub = nint((nsrank+1)*global_nsn*rnstasks)+uhalo
    local_nsn = nsub-nslb+1
    own_nsn = local_nsn-lhalo-uhalo
    nsn = local_nsn

    west = this_rank-1
    if ((west/ewtasks<this_rank/ewtasks).or.(west<0)) west = west+ewtasks
    east = this_rank+1
    if (east/ewtasks>this_rank/ewtasks) east = east-ewtasks
    south = this_rank-ewtasks
    if (south<0) south = south+tasks
    north = this_rank+ewtasks
    if (north>=tasks) north = north-tasks

    ! Check that haven't split up the problem too much.  Idea is that do not want halos overlapping in either dimension.
    ! local_* - lhalo - uhalo is the actual number of non-halo cells on a processor.
    if ((local_nsn - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "NS halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if ((local_ewn  - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "EW halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    ! Print grid geometry
!    write(*,*) "Process ", this_rank, " Total = ", tasks, " ewtasks = ", ewtasks, " nstasks = ", nstasks
!    write(*,*) "Process ", this_rank, " ewrank = ", ewrank, " nsrank = ", nsrank
!    write(*,*) "Process ", this_rank, " l_ewn = ", local_ewn, " o_ewn = ", own_ewn
!    write(*,*) "Process ", this_rank, " l_nsn = ", local_nsn, " o_nsn = ", own_nsn
!    write(*,*) "Process ", this_rank, " ewlb = ", ewlb, " ewub = ", ewub
!    write(*,*) "Process ", this_rank, " nslb = ", nslb, " nsub = ", nsub
!    write(*,*) "Process ", this_rank, " east = ", east, " west = ", west
!    write(*,*) "Process ", this_rank, " north = ", north, " south = ", south
!    write(*,*) "Process ", this_rank, " ew_vars = ", own_ewn, " ns_vars = ", own_nsn
    call distributed_print_grid(own_ewn, own_nsn)

  end subroutine distributed_grid

  function distributed_owner(ew,ewn,ns,nsn)
    implicit none
    logical :: distributed_owner
    integer :: ew,ewn,ns,nsn
    ! begin
    distributed_owner = (ew>lhalo.and.ew<=local_ewn-uhalo.and.&
         ns>lhalo.and.ns<=local_nsn-uhalo)
  end function distributed_owner

  subroutine distributed_print_grid(l_ewn,l_nsn)
    ! Gathers and prints the overall grid layout by processor counts.
    use mpi_mod
    implicit none

    integer :: l_ewn, l_nsn
    integer :: i,j,curr_count
    integer,dimension(2) :: mybounds
    integer,dimension(:,:),allocatable :: bounds

    ! begin
    mybounds(1) = l_ewn
    mybounds(2) = l_nsn

    if (main_task) then
       allocate(bounds(2,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          if (bounds(1,i) .ne. -1) then
             ! total up number of processors with matching distribution
             curr_count = 1
             do j = i+1,tasks
                if ((bounds(1,i) .eq. bounds(1,j)) .and. (bounds(2,i) .eq. bounds(2,j))) then
                   ! if matching current distribution, increment counter
                   curr_count = curr_count + 1
                   bounds(1,j) = -1  ! mark so not counted later
                   bounds(2,j) = -1
                endif
             enddo
             write(*,*) "Layout(EW,NS) = ", bounds(1,i), bounds(2,i), " total procs = ", curr_count
          endif
       end do
    end if
    ! automatic deallocation
  end subroutine distributed_print_grid

  subroutine distributed_print_integer_2d(name,values)
    use mpi_mod
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: global_values,sendbuf

    ! begin
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_int(sendbuf,size(sendbuf),mpi_integer,&
       recvbuf,recvcounts,displs,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,1)<local_ewn) then
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if
    ! automatic deallocation
  end subroutine distributed_print_integer_2d

  subroutine distributed_print_real8_2d(name,values)
    use mpi_mod
    implicit none
    character(*) :: name
    real(8),dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,1)<local_ewn) then
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if
    ! automatic deallocation
  end subroutine distributed_print_real8_2d

  subroutine distributed_print_real8_3d(name,values)
    use mpi_mod
    implicit none
    character(*) :: name
    real(8),dimension(:,:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:,:),allocatable :: global_values,sendbuf

    ! begin
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(size(values,1),minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)*size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(size(values,1),mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:,:) = values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    sendbuf(:,mybounds(1):mybounds(2),mybounds(3):mybounds(4)) = sendbuf(:,mybounds(1):mybounds(2),mybounds(3):mybounds(4))
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(:,bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/size(values,1),bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,2)<local_ewn) then
          do j = lbound(global_values,3),ubound(global_values,3)
             do i = lbound(global_values,2),ubound(global_values,2)
                write(u,'(2i6,100g15.5e3)') j,i,global_values(:,i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,3),ubound(global_values,3)
             do i = lbound(global_values,2),ubound(global_values,2)
                write(u,'(2i6,100g15.5e3)') j,i,global_values(:,i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if
    ! automatic deallocation
  end subroutine distributed_print_real8_3d

  function distributed_put_var_integer_2d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_put_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_int(sendbuf,size(sendbuf),mpi_integer,&
       recvbuf,recvcounts,displs,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_integer_2d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
    end if
    call broadcast(distributed_put_var_integer_2d)
    !automatic deallocation
  end function distributed_put_var_integer_2d

  function distributed_put_var_real4_1d(ncid,varid,values)
    use mpi_mod
    use netcdf
    implicit none
    integer :: distributed_put_var_real4_1d,ncid,varid
    real(4),dimension(:) :: values

    integer :: i,ierror,myn,status,x0id,x1id,y0id,y1id

    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: global_values,recvbuf

    ! begin

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(x1id)
    call broadcast(y0id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       global_values(:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    call fc_gatherv_real4(values,size(values),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i)) = &
               recvbuf(displs(i)+1:displs(i+1))
       end do
       distributed_put_var_real4_1d = &
            nf90_put_var(ncid,varid,global_values(1:myn))
    end if
    call broadcast(distributed_put_var_real4_1d)
    !automatic deallocation
  end function distributed_put_var_real4_1d

  function distributed_put_var_real4_2d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_put_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: recvbuf
    real(4),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real4(sendbuf,size(sendbuf),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_real4_2d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
    end if
    call broadcast(distributed_put_var_real4_2d)
    !automatic deallocation
  end function distributed_put_var_real4_2d

  function distributed_put_var_real8_2d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_put_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_real8_2d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
    end if
    call broadcast(distributed_put_var_real8_2d)
    !automatic deallocation
  end function distributed_put_var_real8_2d

  function distributed_put_var_real8_3d(ncid,varid,values,start)
    use mpi_mod
    implicit none
    integer :: distributed_put_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:,:) :: values

    integer :: ew,i,ierror,ns,nz
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:,:),allocatable :: global_values,sendbuf

    ! begin

    nz = size(values,3)
    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:)),nz))
       global_values(:,:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)&
            *nz
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4),nz))
    sendbuf(:,:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo,:)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i),:) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1,nz/))
       end do
       distributed_put_var_real8_3d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns,:),start)
    end if
    call broadcast(distributed_put_var_real8_3d)
    !automatic deallocation
  end function distributed_put_var_real8_3d

  subroutine distributed_scatter_var_integer_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    use mpi_mod
    implicit none
    integer,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    integer,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:),allocatable :: sendbuf
    integer,dimension(:,:),allocatable :: recvbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_integer,&
         recvbuf,size(recvbuf),mpi_integer,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine distributed_scatter_var_integer_2d

  subroutine distributed_scatter_var_logical_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    use mpi_mod
    implicit none
    logical,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    logical,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    logical,dimension(:),allocatable :: sendbuf
    logical,dimension(:,:),allocatable :: recvbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_logical,&
         recvbuf,size(recvbuf),mpi_logical,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine distributed_scatter_var_logical_2d

  subroutine distributed_scatter_var_real4_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    use mpi_mod
    implicit none
    real(4),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(4),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(4),dimension(:),allocatable :: sendbuf
    real(4),dimension(:,:),allocatable :: recvbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                     d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                     (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         recvbuf,size(recvbuf),mpi_real4,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine distributed_scatter_var_real4_2d

  subroutine distributed_scatter_var_real4_3d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    use mpi_mod
    implicit none
    real(4),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(4),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(4),dimension(:),allocatable :: sendbuf
    real(4),dimension(:,:,:),allocatable :: recvbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)*size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(:,&
                                   d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         recvbuf,size(recvbuf),mpi_real4,main_rank,comm,ierror)
    values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine distributed_scatter_var_real4_3d

  subroutine distributed_scatter_var_real8_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    use mpi_mod
    implicit none
    real(8),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(8),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(8),dimension(:),allocatable :: sendbuf
    real(8),dimension(:,:),allocatable :: recvbuf

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine distributed_scatter_var_real8_2d

  subroutine distributed_scatter_var_real8_3d(values, global_values, deallocflag)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    use mpi_mod
    implicit none
    real(8),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(8),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    logical,optional :: deallocflag
    logical :: deallocmem

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(8),dimension(:),allocatable :: sendbuf
    real(8),dimension(:,:,:),allocatable :: recvbuf

    if (present(deallocflag)) then
       deallocmem = deallocflag
    else
       deallocmem = .true.
    endif

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif

       d_gs_mybounds(1) = ewlb+lhalo
       d_gs_mybounds(2) = ewub-uhalo
       d_gs_mybounds(3) = nslb+lhalo
       d_gs_mybounds(4) = nsub-uhalo
       call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
          mpi_integer,main_rank,comm)
    endif

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)&
                      *size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(:,&
                                   d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:,:)

    if (deallocmem) deallocate(global_values)
    ! automatic deallocation
  end subroutine distributed_scatter_var_real8_3d

  subroutine global_sum(x)
    use mpi_mod
    implicit none
    real(8),dimension(:) :: x
    
    integer :: ierror
    real(8),dimension(size(x)) :: sum
    ! begin
    call mpi_allreduce(x,sum,size(x),mpi_real8,mpi_sum,comm,ierror)
    x(:) = sum(:)
  end subroutine global_sum

  subroutine not_parallel(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    call parallel_stop(file,line)
  end subroutine not_parallel

  subroutine parallel_barrier
    use mpi_mod
    implicit none
    integer :: ierror
    ! begin
    call mpi_barrier(comm,ierror)
  end subroutine parallel_barrier

  function parallel_boundary(ew,ewn,ns,nsn)
    implicit none
    logical :: parallel_boundary
    integer :: ew,ewn,ns,nsn
    ! begin
    parallel_boundary = (ewlb<1.and.ew==1+lhalo).or.&
         (ewub>global_ewn.and.ew==ewn-uhalo).or.&
         (nslb<1.and.ns==1+lhalo).or.&
         (nsub>global_nsn.and.ns==nsn-uhalo)
  end function parallel_boundary

  function parallel_close(ncid)
    implicit none
    integer :: ncid,parallel_close
    ! begin
    if (main_task) parallel_close = nf90_close(ncid)
    call broadcast(parallel_close)
  end function parallel_close

  function parallel_create(path,cmode,ncid)
    implicit none
    integer :: cmode,ncid,parallel_create
    character(len=*) :: path
    ! begin
    if (main_task) parallel_create = nf90_create(path,cmode,ncid)
    call broadcast(parallel_create)
    call broadcast(ncid)
  end function parallel_create
    
  function parallel_def_dim(ncid,name,len,dimid)
    use netcdf
    implicit none
    integer :: dimid,len,ncid,parallel_def_dim
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_dim = nf90_def_dim(ncid,name,len,dimid)
    call broadcast(parallel_def_dim)
    call broadcast(dimid)
  end function parallel_def_dim

  function parallel_def_var_dimids(ncid,name,xtype,dimids,varid)
    implicit none
    integer :: ncid,parallel_def_var_dimids,varid,xtype
    integer,dimension(:) :: dimids
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_var_dimids = &
         nf90_def_var(ncid,name,xtype,dimids,varid)
    call broadcast(parallel_def_var_dimids)
    call broadcast(varid)
  end function parallel_def_var_dimids

  function parallel_def_var_nodimids(ncid,name,xtype,varid)
    implicit none
    integer :: ncid,parallel_def_var_nodimids,varid,xtype
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_var_nodimids = &
         nf90_def_var(ncid,name,xtype,varid)
    call broadcast(parallel_def_var_nodimids)
    call broadcast(varid)
  end function parallel_def_var_nodimids

  function parallel_enddef(ncid)
    implicit none
    integer :: ncid,parallel_enddef
    ! begin
    if (main_task) parallel_enddef = nf90_enddef(ncid)
    call broadcast(parallel_enddef)
  end function parallel_enddef

  subroutine parallel_finalise
    use mpi_mod
    implicit none
    integer :: ierror
    ! begin
    call mpi_finalize(ierror)
  end subroutine parallel_finalise

  function parallel_get_att_character(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_character,varid
    character(len=*) :: name,values
    ! begin
    if (main_task) parallel_get_att_character = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_character)
    call broadcast(values)
  end function parallel_get_att_character

  function parallel_get_att_real4(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real4,varid
    character(len=*) :: name
    real(4) :: values
    ! begin
    if (main_task) parallel_get_att_real4 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4)
    call broadcast(values)
  end function parallel_get_att_real4

  function parallel_get_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real4_1d,varid
    character(len=*) :: name
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real4_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4_1d)
    call broadcast(values)
  end function parallel_get_att_real4_1d

  function parallel_get_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8,varid
    character(len=*) :: name
    real(8) :: values
    ! begin
    if (main_task) parallel_get_att_real8 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8)
    call broadcast(values)
  end function parallel_get_att_real8

  function parallel_get_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8_1d,varid
    character(len=*) :: name
    real(8),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real8_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8_1d)
    call broadcast(values)
  end function parallel_get_att_real8_1d

  function parallel_get_var_integer_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_integer_1d,varid
    integer,dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_integer_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_integer_1d)
    call broadcast(values)
  end function parallel_get_var_integer_1d

  function parallel_get_var_real4_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_real4_1d,varid
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_real4_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real4_1d)
    call broadcast(values)
  end function parallel_get_var_real4_1d

  function parallel_get_var_real8_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_real8_1d,varid
    real(8),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_real8_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real8_1d)
    call broadcast(values)
  end function parallel_get_var_real8_1d

  !TODO - Pass locew in first position and locns in second position?
  !TODO - Remove his function if no longer needed?

  function parallel_globalID(locns, locew, upstride)
    ! Returns a unique ID for a given row and column reference that is identical across all processors.
    ! For instance if Proc 0: (17,16) is the same global cell as Proc 3: (17,1), then the globalID will be the same for both.
    ! These IDs are spaced upstride apart.  upstride = number of vertical layers.  Typically (upn) + number of ghost layers (2 = top and bottom)
    integer,intent(IN) :: locns, locew, upstride
    integer :: parallel_globalID
    ! locns is local NS (row) grid index
    ! locew is local EW (col) grid index
    integer :: global_row, global_col, global_ID
    character(len=40) :: local_coord

    ! including global domain halo adds lhalo to offsets
    global_row = (locns - lhalo) + (global_row_offset + lhalo)
    global_col = (locew - lhalo) + (global_col_offset + lhalo)

    ! if halo cell and if using periodic boundary conditions,
    ! define global ID to be associated non-halo cell
    if (global_row .le. lhalo) then
      if (horiz_bcs_type_south .eq. HORIZ_BCS_CYCLIC) then
        global_row = global_row + global_nsn
      endif
    endif

    if (global_row > (global_nsn+lhalo)) then
      if  (horiz_bcs_type_north .eq. HORIZ_BCS_CYCLIC) then
        global_row = global_row - global_nsn
      endif
    endif

    if (global_col .le. lhalo) then
      if (horiz_bcs_type_west .eq. HORIZ_BCS_CYCLIC) then
        global_col = global_col + global_ewn
      endif
    endif

    if (global_col > (global_ewn+lhalo)) then
      if (horiz_bcs_type_east .eq. HORIZ_BCS_CYCLIC) then
        global_col = global_col - global_ewn
      endif
    endif

    ! including global domain halo adds (lhalo + uhalo) to global_ewn
    global_ID = ((global_row - 1) * (global_ewn + lhalo + uhalo) + (global_col - 1)) * upstride + 1

    ! JEFF Testing Code
    ! write(local_coord, "A13,I10.1,A2,I10.1,A1") " (NS, EW) = (", locns, ", ", locew, ")"
    ! write(*,*) "Processor reference ", this_rank, local_coord, " globalID = ", global_ID

    !return value
    parallel_globalID = global_ID

  end function parallel_globalID


  function parallel_globalID_scalar(locew, locns, upstride)

    !WHL - This function is similar to parallel_globalID, but assigns 0's to cells outside the global domain

    ! Returns a unique ID for a given row and column reference that is identical across all processors.
    ! For instance if Proc 0: (17,16) is the same global cell as Proc 3: (17,1), then the globalID will be the same for both.
    ! These IDs are spaced upstride apart.  upstride = number of vertical layers.
    integer,intent(IN) :: locns, locew, upstride
    integer :: parallel_globalID_scalar
    ! locns is local NS (row) grid index
    ! locew is local EW (col) grid index
    integer :: global_row, global_col, global_ID
    character(len=40) :: local_coord

    ! including global domain halo adds lhalo to offsets
    global_row = (locns - lhalo) + global_row_offset
    global_col = (locew - lhalo) + global_col_offset

    ! including global domain halo adds (lhalo + uhalo) to global_ewn
    global_ID = ((global_row - 1)*(global_ewn) + (global_col - 1)) * upstride + 1

    ! JEFF Testing Code
    ! write(local_coord, "A13,I10.1,A2,I10.1,A1") " (NS, EW) = (", locns, ", ", locew, ")"
    ! write(*,*) "Processor reference ", this_rank, local_coord, " globalID = ", global_ID

    !return value
    parallel_globalID_scalar = global_ID

  end function parallel_globalID_scalar

  subroutine parallel_halo_integer_2d(a)
    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    
    integer :: erequest,ierror,nrequest,srequest,wrequest
    integer,dimension(lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(local_ewn,lhalo) :: nsend,srecv
    integer,dimension(local_ewn,uhalo) :: nrecv,ssend

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
      a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)
  end subroutine parallel_halo_integer_2d

  subroutine parallel_halo_logical_2d(a)
    use mpi_mod
    implicit none
    logical,dimension(:,:) :: a

    integer :: erequest,ierror,nrequest,srequest,wrequest
    logical,dimension(lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    logical,dimension(uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    logical,dimension(local_ewn,lhalo) :: nsend,srecv
    logical,dimension(local_ewn,uhalo) :: nrecv,ssend

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_logical,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_logical,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_logical,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_logical,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_logical,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_logical,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_logical,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_logical,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)
  end subroutine parallel_halo_logical_2d

  subroutine parallel_halo_real4_2d(a)
    use mpi_mod
    implicit none
    real(4),dimension(:,:) :: a

    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(4),dimension(lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(4),dimension(uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(4),dimension(local_ewn,lhalo) :: nsend,srecv
    real(4),dimension(local_ewn,uhalo) :: nrecv,ssend

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real4,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real4,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real4,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real4,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real4,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real4,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real4,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real4,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)
  end subroutine parallel_halo_real4_2d


  subroutine parallel_halo_real8_2d(a, periodic_offset_ew, periodic_offset_ns)

    !WHL - added optional arguments for periodic offsets, to support ismip-hom test cases

    use mpi_mod
    implicit none
    real(8),dimension(:,:) :: a
    real(8), intent(in), optional :: &
       periodic_offset_ew,  &! offset halo values by this amount
                             ! if positive, the offset is positive for W halo, negative for E halo
       periodic_offset_ns    ! offset halo values by this amount
                             ! if positive, the offset is positive for S halo, negative for N halo
    
    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(8),dimension(lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(local_ewn,lhalo) :: nsend,srecv
    real(8),dimension(local_ewn,uhalo) :: nrecv,ssend

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    if (present(periodic_offset_ew)) then
       if (periodic_offset_ew /= 0.d0) then
          if (this_rank <= west) then   ! this proc lies at the west edge of the global domain
!             print*, 'Offset at west edge: this_rank, west =', this_rank, west
             a(:lhalo,1+lhalo:local_nsn-uhalo) =   &
                a(:lhalo,1+lhalo:local_nsn-uhalo) + periodic_offset_ew
          endif
          if (this_rank >= east) then   ! this proc lies at the east edge of the global domain
!             print*, 'Offset at east edge: this_rank, east =', this_rank, east
             a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) =    &
                a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) - periodic_offset_ew
          endif
       endif
    endif

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (present(periodic_offset_ns)) then
       if (periodic_offset_ns /= 0.d0) then
          if (this_rank <= south) then  ! this proc lies at the south edge of the global domain
!             print*, 'Offset at south edge: this_rank, south =', this_rank, south
             a(:,:lhalo) = a(:,:lhalo) + periodic_offset_ns
          endif
          if (this_rank >= north) then  ! this proc lies at the north edge of the global domain
!             print*, 'Offset at north edge: this_rank, north =', this_rank, north
             a(:,local_nsn-uhalo+1:) = a(:,local_nsn-uhalo+1:) - periodic_offset_ns
          endif
       endif
    endif

  end subroutine parallel_halo_real8_2d

  subroutine parallel_halo_real8_3d(a)

    use mpi_mod
    implicit none
    real(8),dimension(:,:,:) :: a

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(8),dimension(size(a,1),lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(size(a,1),uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(size(a,1),local_ewn,lhalo) :: nsend,srecv
    real(8),dimension(size(a,1),local_ewn,uhalo) :: nrecv,ssend

    ! begin

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ") &
                 and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,local_nsn-uhalo+1:) = nrecv(:,:,:)

  end subroutine parallel_halo_real8_3d

  !TODO - Remove this subroutine?
  subroutine parallel_halo_temperature(a)
    !JEFF This routine is for updating the halo for the variable model%temper%temp.
    ! This variable is two larger in each dimension, because of the current advection code.
    ! Per Bill L, we will remove this difference when we update the remapping code.
    use mpi_mod
    implicit none
    real(8),dimension(:,:,:) :: a

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(8),dimension(size(a,1),lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(size(a,1),uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(size(a,1),local_ewn,lhalo) :: nsend,srecv
    real(8),dimension(size(a,1),local_ewn,uhalo) :: nrecv,ssend

    ! begin
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
      a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    !JEFF Change middle index to put halo into correct place for temperature.
    a(:,2:lhalo+1,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    !JEFF Put an upper bound on middle index to prevent overrunning index
    a(:,local_ewn-uhalo+1:local_ewn-uhalo+2,1+lhalo:local_nsn-uhalo) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    !JEFF Change middle to move one index further in
    ssend(:,:,:) = a(:,2:local_ewn+1,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    !JEFF Change middle to move one index further in
    a(:,2:local_ewn+1,:lhalo) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    !JEFF Limit last index to size of halo and middle to one index further in
    a(:,2:local_ewn+1,local_nsn-uhalo+1:local_nsn) = nrecv(:,:,:)
  end subroutine parallel_halo_temperature

  function parallel_halo_verify_integer_2d(a)
    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a

    integer :: erequest,ierror,nrequest,srequest,wrequest
    integer,dimension(lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(local_ewn,lhalo) :: nsend,srecv
    integer,dimension(local_ewn,uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_integer_2d

    ! begin

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    ! ANY True if any value is true (LOGICAL)
    notverify_flag = ANY(a(:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:))

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:lhalo) /= srecv(:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,local_nsn-uhalo+1:) /= nrecv(:,:))

    ! if notverify_flag is TRUE, then there was some difference detected
    if (notverify_flag) then
         write(*,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_integer_2d = .NOT. notverify_flag  ! return if verified (True) or not verified (False)
  end function parallel_halo_verify_integer_2d

  function parallel_halo_verify_real8_2d(a)
    use mpi_mod
    implicit none
    real(8),dimension(:,:) :: a
    
    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(8),dimension(lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(local_ewn,lhalo) :: nsend,srecv
    real(8),dimension(local_ewn,uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_real8_2d

    ! begin

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    notverify_flag = ANY(a(:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:))

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:lhalo) /= srecv(:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,local_nsn-uhalo+1:) /= nrecv(:,:))

    if (notverify_flag) then
         write(*,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_real8_2d = .NOT. notverify_flag
  end function parallel_halo_verify_real8_2d

  function parallel_halo_verify_real8_3d(a)
    use mpi_mod
    implicit none
    real(8),dimension(:,:,:) :: a
    
    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(8),dimension(size(a,1),lhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(size(a,1),uhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(size(a,1),local_ewn,lhalo) :: nsend,srecv
    real(8),dimension(size(a,1),local_ewn,uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_real8_3d

    ! begin

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    notverify_flag = ANY(a(:,:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:,:))

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:,:lhalo) /= srecv(:,:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:,local_nsn-uhalo+1:) /= nrecv(:,:,:))

    if (notverify_flag) then
         write(*,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_real8_3d = .NOT. notverify_flag
  end function parallel_halo_verify_real8_3d

  ! parallel_initialise should generally just be called by standalone cism drivers
  ! When cism is nested inside a climate model (so mpi_init has already been called) use parallel_set_info instead

  subroutine parallel_initialise
    use mpi_mod 
    implicit none
    integer :: ierror 
    integer, parameter :: my_main_rank = 0
    ! begin 
    call mpi_init(ierror)
    call parallel_set_info(mpi_comm_world, my_main_rank)
  end subroutine parallel_initialise

  ! parallel_set_info should be called directly when cism is nested inside a climate model
  ! (then, mpi_init has already been called, so do NOT use parallel_initialise)

  subroutine parallel_set_info(my_comm, my_main_rank)
    use mpi_mod
    implicit none
    integer, intent(in) :: my_comm       ! CISM's global communicator
    integer, intent(in) :: my_main_rank  ! rank of the master task
    integer :: ierror 
    ! begin
    comm = my_comm
    main_rank = my_main_rank
    call mpi_comm_size(comm,tasks,ierror)
    call mpi_comm_rank(comm,this_rank,ierror)
    main_task = (this_rank==main_rank)
  end subroutine parallel_set_info

  function parallel_inq_attname(ncid,varid,attnum,name)
    implicit none
    integer :: attnum,ncid,parallel_inq_attname,varid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_attname = &
         nf90_inq_attname(ncid,varid,attnum,name)
    call broadcast(parallel_inq_attname)
    call broadcast(name)
  end function parallel_inq_attname

  function parallel_inq_dimid(ncid,name,dimid)
    implicit none
    integer :: dimid,ncid,parallel_inq_dimid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_dimid = nf90_inq_dimid(ncid,name,dimid)
    call broadcast(parallel_inq_dimid)
    call broadcast(dimid)
  end function parallel_inq_dimid

  function parallel_inq_varid(ncid,name,varid)
    implicit none
    integer :: ncid,parallel_inq_varid,varid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_varid = nf90_inq_varid(ncid,name,varid)
    call broadcast(parallel_inq_varid)
    call broadcast(varid)
  end function parallel_inq_varid

  function parallel_inquire(ncid,nvariables)
    implicit none
    integer :: ncid,parallel_inquire,nvariables
    ! begin
    if (main_task) parallel_inquire = nf90_inquire(ncid,nvariables=nvariables)
    call broadcast(parallel_inquire)
    call broadcast(nvariables)
  end function parallel_inquire

  function parallel_inquire_dimension(ncid,dimid,name,len)
    implicit none
    integer :: dimid,ncid,parallel_inquire_dimension
    integer,optional :: len
    character(len=*),optional :: name
    
    integer :: l
    
    ! begin

    if (present(name)) then
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,name,len=l)
       call broadcast(name)
    else
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,len=l)
    end if
    call broadcast(parallel_inquire_dimension)
    if (present(len)) then
       call broadcast(l)
       len = l
    end if
  end function parallel_inquire_dimension

  function parallel_inquire_variable(ncid,varid,name,ndims,dimids,natts)
    implicit none
    integer :: ncid,parallel_inquire_variable,varid
    integer,optional :: ndims,natts
    character(len=*),optional :: name
    integer,dimension(:),optional :: dimids

    integer :: nd,na
    ! begin
    if (present(name)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,name=name)
       call broadcast(parallel_inquire_variable)
       call broadcast(name)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (present(dimids)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,dimids=dimids)
       call broadcast(parallel_inquire_variable)
       call broadcast(dimids)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (main_task) parallel_inquire_variable = &
         nf90_inquire_variable(ncid,varid,ndims=nd,natts=na)
    call broadcast(parallel_inquire_variable)
    if (present(ndims)) then
       call broadcast(nd)
       ndims = nd
    end if
    if (present(natts)) then
       call broadcast(na)
       natts = na
    end if
  end function parallel_inquire_variable

  function parallel_open(path,mode,ncid)
    implicit none
    integer :: mode,ncid,parallel_open
    character(len=*) :: path
    ! begin
    if (main_task) parallel_open = nf90_open(path,mode,ncid)
    call broadcast(parallel_open)
  end function parallel_open

  subroutine parallel_print_all(name,values)
    implicit none
    character(*) :: name
    real(8),dimension(:,:,:) :: values

    integer,parameter :: u = 33
    integer :: i,j,t
    ! begin
    if (main_task) then
       open(unit=u,file=name,form="formatted",status="replace")
       close(u)
    end if
    do t = 0,tasks-1
       call parallel_barrier
       if (t==this_rank) then
          open(unit=u,file=name,form="formatted",position="append")
          do j = 1,size(values,3)
             do i = 1,size(values,2)
                write(u,'(2i5,100g15.5e3)') nslb+j-1,ewlb+i-1,values(:,i,j)
             end do
             write(u,'()')
          end do
          write(u,'(//)')
          close(u)
       end if
    end do
  end subroutine parallel_print_all

  subroutine parallel_print_integer_2d(name,values)
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j

    ! begin
    if (main_task) then
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
       close(u)
    end if

    call parallel_barrier  ! Only the main_task writes the variable.  Rest wait here.
    ! automatic deallocation
  end subroutine parallel_print_integer_2d

  subroutine parallel_print_real8_2d(name,values)
    implicit none
    character(*) :: name
    real(8),dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j

    ! begin
    if (main_task) then
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
       close(u)
    end if

    call parallel_barrier  ! Only the main_task writes the variable.  Rest wait here.
  end subroutine parallel_print_real8_2d

  subroutine parallel_print_real8_3d(name,values)
    implicit none
    character(*) :: name
    real(8),dimension(:,:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j

    ! begin
    if (main_task) then
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       do j = lbound(values,3),ubound(values,3)
          do i = lbound(values,2),ubound(values,2)
             write(u,'(2i6,100g15.5e3)') j,i,values(:,i,j)
          end do
          write(u,'()')
       end do
       close(u)
    end if

    call parallel_barrier  ! Only the main_task writes the variable.  Rest wait here.
  end subroutine parallel_print_real8_3d

  function parallel_put_att_character(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_character,varid
    character(len=*) :: name,values
    ! begin
    if (main_task) parallel_put_att_character = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_character)
  end function parallel_put_att_character

  function parallel_put_att_real4(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real4,varid
    character(len=*) :: name
    real(4) :: values
    ! begin
    if (main_task) parallel_put_att_real4 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4)
  end function parallel_put_att_real4

  function parallel_put_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real4_1d,varid
    character(len=*) :: name
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real4_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4_1d)
  end function parallel_put_att_real4_1d

  function parallel_put_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8,varid
    character(len=*) :: name
    real(8) :: values
    ! begin
    if (main_task) parallel_put_att_real8 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8)
  end function parallel_put_att_real8

  function parallel_put_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8_1d,varid
    character(len=*) :: name
    real(8),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real8_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8_1d)
  end function parallel_put_att_real8_1d

  function parallel_put_var_real4(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real4,varid
    integer,dimension(:) :: start
    real(4) :: values
    ! begin
    if (main_task) parallel_put_var_real4 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real4)
  end function parallel_put_var_real4

  function parallel_put_var_real8(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8,varid
    integer,dimension(:) :: start
    real(8) :: values
    ! begin
    if (main_task) parallel_put_var_real8 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real8)
  end function parallel_put_var_real8

  function parallel_put_var_real8_1d(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8_1d,varid
    integer,dimension(:),optional :: start
    real(8),dimension(:) :: values
    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values)
       end if
    end if
    call broadcast(parallel_put_var_real8_1d)
  end function parallel_put_var_real8_1d

  function parallel_redef(ncid)
    implicit none
    integer :: ncid,parallel_redef
    ! begin
    if (main_task) parallel_redef = nf90_redef(ncid)
    call broadcast(parallel_redef)
  end function parallel_redef

  function parallel_reduce_sum(x)
    use mpi_mod
    implicit none
    real(8) :: x

    integer :: ierror
    real(8) :: recvbuf,sendbuf, parallel_reduce_sum
    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_sum,comm,ierror)
    parallel_reduce_sum = recvbuf
    return
  end function parallel_reduce_sum

  function parallel_reduce_max_integer(x)
    use mpi_mod
    implicit none
    integer :: x

    integer :: ierror
    integer :: recvbuf,sendbuf, parallel_reduce_max_integer
    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer,mpi_max,comm,ierror)
    parallel_reduce_max_integer = recvbuf
    return
  end function parallel_reduce_max_integer

  function parallel_reduce_max_real4(x)
    use mpi_mod
    implicit none
    real(4) :: x

    integer :: ierror
    real(4) :: recvbuf,sendbuf, parallel_reduce_max_real4
    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real4,mpi_max,comm,ierror)
    parallel_reduce_max_real4 = recvbuf
    return
  end function parallel_reduce_max_real4

  function parallel_reduce_max_real8(x)
    use mpi_mod
    implicit none
    real(8) :: x

    integer :: ierror
    real(8) :: recvbuf,sendbuf, parallel_reduce_max_real8
    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_max,comm,ierror)
    parallel_reduce_max_real8 = recvbuf
    return
  end function parallel_reduce_max_real8

  ! Andy removed support for returnownedvector in October 2011.
  ! subroutine parallel_set_trilinos_return_vect
    ! Trilinos can return the full solution to each node or just the owned portion
    ! For parallel_mpi mode only the owned portion is expected
  !  call returnownedvector()  ! in trilinosLinearSolver.cpp
  ! end subroutine parallel_set_trilinos_return_vect

  subroutine parallel_show_minmax(label,values)
    use mpi_mod
    implicit none
    character(*) :: label
    real(8),dimension(:,:,:) :: values
    
    integer :: ierror
    real(8) :: allmin,allmax,mymin,mymax
    ! begin
    mymin = minval(values(:,1+lhalo:size(values,2)-uhalo,&
         1+lhalo:size(values,3)-uhalo))
    mymax = maxval(values(:,1+lhalo:size(values,2)-uhalo,&
         1+lhalo:size(values,3)-uhalo))
    call mpi_reduce(mymin,allmin,1,mpi_real8,mpi_min,main_rank,comm,ierror)
    call mpi_reduce(mymax,allmax,1,mpi_real8,mpi_max,main_rank,comm,ierror)
    if (main_task) print *,label,allmin,allmax
  end subroutine parallel_show_minmax

  subroutine parallel_stop(file,line)
    use mpi_mod
    implicit none
    integer :: line
    character(len=*) :: file
    integer :: ierror
    ! begin
    if (main_task) write(0,*) "PARALLEL STOP in ",file," at line ",line
    call mpi_finalize(ierror)
    stop "PARALLEL STOP"
  end subroutine parallel_stop

  function parallel_sync(ncid)
    implicit none
    integer :: ncid,parallel_sync
    ! begin
    if (main_task) parallel_sync = nf90_sync(ncid)
    call broadcast(parallel_sync)
  end function parallel_sync

  !TODO - This subroutine is called only from periodic_boundaries subroutine, which is not
  !       currently used.  Remove it?

  subroutine parallel_velo_halo(a)
    use mpi_mod
    implicit none
    real(8),dimension(:,:) :: a

    integer :: ierror,nrequest,erequest
    real(8),dimension(size(a,2)-lhalo-uhalo+1) :: wsend,erecv
    real(8),dimension(size(a,1)-lhalo) :: ssend,nrecv

    ! begin
    if (size(a,1)/=local_ewn-1.or.size(a,2)/=local_nsn-1) &
         call parallel_stop(__FILE__,__LINE__)

    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    wsend(:) = a(1+lhalo,1+lhalo:size(a,2)-uhalo+1)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(size(a,1),1+lhalo:size(a,2)-uhalo+1) = erecv(:)

    ssend(:) = a(1+lhalo:,1+lhalo)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(1+lhalo:,size(a,2)) = nrecv(:)

  end subroutine parallel_velo_halo


  subroutine staggered_parallel_halo_integer_2d(a)
    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a

    ! Implements a staggered grid halo update.
    ! As the grid is staggered, the array 'a' is one smaller in both dimensions than an unstaggered array.

    ! The grid is laid out from the SW, and the lower left corner is assigned to this_rank = 0.
    ! It's eastern nbhr is task_id = 1, proceeding rowwise and starting from the western edge.
    ! The South-most processes own one additional row of stagggered variables on the southern edge
    ! and have one less 'southern' halo row than other processes. Likewise, the West-most processes own one 
    ! additional column of staggered variables on the western edge and have one less 'western' halo column. 
    ! This is implemented by a modification to the staggered_lhalo value on these processes. 

    !WHL - I don't think we need to say that the South-most processes "own" an additional row of
    !      staggered variables on the southern edge.  I think we can treat the southern edge as a halo row
    !      and still enforce the various global BC we want.

    ! Maintaining global boundary conditions are not addressed within this routine (yet).

    ! integer :: erequest,ierror,one,nrequest,srequest,wrequest
    integer :: ierror,nrequest,srequest,erequest,wrequest

    integer,dimension(staggered_lhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: esend,wrecv
    integer,dimension(staggered_uhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: erecv,wsend
    integer,dimension(size(a,1),staggered_lhalo) :: nsend,srecv
    integer,dimension(size(a,1),staggered_uhalo) :: nrecv,ssend

    !WHL - I defined a logical variable to determine whether or not to fill halo cells
    !      at the edge of the global domain.  I am setting it to true by default to support
    !      cyclic global BCs.

    !TODO - Set to true in all cases? (Or simply remove the 'if's?)  

    logical :: fill_global_halos = .true.

    ! begin

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east .or. fill_global_halos) then
      call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
!      wsend(:,1:size(a,2)-staggered_shalo-staggered_nhalo) = &
!        a(1+staggered_whalo:1+staggered_whalo+staggered_ehalo-1, &
!            1+staggered_shalo:size(a,2)-staggered_nhalo)
      wsend(:,1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
        a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
            1+staggered_lhalo:size(a,2)-staggered_uhalo)
      call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
!      esend(:,1:size(a,2)-staggered_shalo-staggered_nhalo) = &
!        a(size(a,1)-staggered_ehalo-staggered_whalo+1:size(a,1)-staggered_ehalo, &
!            1+staggered_shalo:size(a,2)-staggered_nhalo)
      esend(:,1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
        a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
            1+staggered_lhalo:size(a,2)-staggered_uhalo)
      call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
!      a(size(a,1)-staggered_ehalo+1:size(a,1), &
!          1+staggered_shalo:size(a,2)-staggered_nhalo) = &
!          erecv(:,1:size(a,2)-staggered_shalo-staggered_nhalo)
      a(size(a,1)-staggered_uhalo+1:size(a,1), &
          1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
          erecv(:,1:size(a,2)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
!      a(1:staggered_whalo, &
!          1+staggered_shalo:size(a,2)-staggered_nhalo) = &
!          wrecv(:,1:size(a,2)-staggered_shalo-staggered_nhalo)
      a(1:staggered_lhalo, &
          1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
          wrecv(:,1:size(a,2)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > south .or. fill_global_halos) then
      ssend(:,:) = &
!        a(:,1+staggered_shalo:1+staggered_shalo+staggered_nhalo-1)
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      nsend(:,:) = &
!        a(:,size(a,2)-staggered_nhalo-staggered_shalo+1:size(a,2)-staggered_nhalo)
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
!      a(:,size(a,2)-staggered_nhalo+1:size(a,2)) = nrecv(:,:)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2)) = nrecv(:,:)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
!      a(:,1:staggered_shalo) = srecv(:,:)
      a(:,1:staggered_lhalo) = srecv(:,:)
    endif

  end subroutine staggered_parallel_halo_integer_2d


  subroutine staggered_parallel_halo_integer_3d(a)
    use mpi_mod
    implicit none
    integer,dimension(:,:,:) :: a

    ! Implements a staggered grid halo update for a 3D field.
    ! As the grid is staggered, the array 'a' is one smaller in both dimensions than an unstaggered array.
    ! The vertical dimension is assumed to be the first index, i.e., a(k,i,j).

    ! The grid is laid out from the SW, and the lower left corner is assigned to this_rank = 0.
    ! It's eastern nbhr is task_id = 1, proceeding rowwise and starting from the western edge.
    ! The South-most processes own one additional row of stagggered variables on the southern edge
    ! and have one less 'southern' halo row than other processes. Likewise, the West-most processes own one 
    ! additional column of staggered variables on the western edge and have one less 'western' halo column. 
    ! This is implemented by a modification to the staggered_lhalo value on these processes. 

    ! Maintaining global boundary conditions are not addressed within this routine (yet).

    ! integer :: erequest,ierror,one,nrequest,srequest,wrequest
    integer :: ierror,nrequest,srequest,erequest,wrequest

    integer,dimension(size(a,1),staggered_lhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: esend,wrecv
    integer,dimension(size(a,1),staggered_uhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: erecv,wsend
    integer,dimension(size(a,1),size(a,2),staggered_lhalo) :: nsend,srecv
    integer,dimension(size(a,1),size(a,2),staggered_uhalo) :: nrecv,ssend

    !WHL - I defined a logical variable to determine whether or not to fill halo cells
    !      at the edge of the global domain.  I am setting it to true by default to support
    !      cyclic global BCs.

    !TODO - Set to true in all cases? (Or simply remove the 'if's?)  

    logical :: fill_global_halos = .true.

    ! begin

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1.or.size(a,3)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. fill_global_halos) then
      call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
      wsend(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
            1+staggered_lhalo:size(a,3)-staggered_uhalo)
      call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
      esend(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
            1+staggered_lhalo:size(a,3)-staggered_uhalo)
      call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2), &
          1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
          erecv(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:,1:staggered_lhalo, &
          1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
          wrecv(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > south .or. fill_global_halos) then
      ssend(:,:,:) = &
        a(:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      nsend(:,:,:) = &
        a(:,:,size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = nrecv(:,:,:)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,1:staggered_lhalo) = srecv(:,:,:)
    endif

  end subroutine staggered_parallel_halo_integer_3d


  !WHL - Edited the original subroutine so that values from N and E edges
  !      of global domain can be written to halo cells at the S and W edges,
  !      to allow cyclic BCs for staggered variables

  subroutine staggered_parallel_halo_real8_2d(a)

    use mpi_mod
    implicit none
    real(8),dimension(:,:) :: a

    ! Implements a staggered grid halo update.
    ! As the grid is staggered, the array 'a' is one smaller in both dimensions than an unstaggered array.

    ! The grid is laid out from the SW, and the lower left corner is assigned to this_rank = 0.
    ! It's eastern nbhr is task_id = 1, proceeding rowwise and starting from the western edge.
    ! The South-most processes own one additional row of stagggered variables on the southern edge
    ! and have one less 'southern' halo row than other processes. Likewise, the West-most processes own one 
    ! additional column of staggered variables on the western edge and have one less 'western' halo column. 
    ! This is implemented by a modification to the staggered_lhalo value on these processes. 

    ! Maintaining global boundary conditions are not addressed within this routine (yet).

    ! integer :: erequest,ierror,one,nrequest,srequest,wrequest
    integer :: ierror,nrequest,srequest,erequest,wrequest
!    real(8),dimension(staggered_whalo,size(a,2)-staggered_shalo-staggered_nhalo) :: esend,wrecv
!    real(8),dimension(staggered_ehalo,size(a,2)-staggered_shalo-staggered_nhalo) :: erecv,wsend
!    real(8),dimension(size(a,1),staggered_shalo) :: nsend,srecv
!    real(8),dimension(size(a,1),staggered_nhalo) :: nrecv,ssend
    real(8),dimension(staggered_lhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: esend,wrecv
    real(8),dimension(staggered_uhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: erecv,wsend
    real(8),dimension(size(a,1),staggered_lhalo) :: nsend,srecv
    real(8),dimension(size(a,1),staggered_uhalo) :: nrecv,ssend

    !WHL - I defined a logical variable to determine whether or not to fill halo cells
    !      at the edge of the global domain.  I am setting it to true by default to support
    !      cyclic global BCs.

    !TODO - Set to true in all cases? (Or simply remove the 'if's?)  

    logical :: fill_global_halos = .true.

    ! begin

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east .or. fill_global_halos) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
!      wsend(:,1:size(a,2)-staggered_shalo-staggered_nhalo) = &
!        a(1+staggered_whalo:1+staggered_whalo+staggered_ehalo-1, &
!            1+staggered_shalo:size(a,2)-staggered_nhalo)
      wsend(:,1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
        a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
            1+staggered_lhalo:size(a,2)-staggered_uhalo)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
!      esend(:,1:size(a,2)-staggered_shalo-staggered_nhalo) = &
!        a(size(a,1)-staggered_ehalo-staggered_whalo+1:size(a,1)-staggered_ehalo, &
!            1+staggered_shalo:size(a,2)-staggered_nhalo)
      esend(:,1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
        a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
            1+staggered_lhalo:size(a,2)-staggered_uhalo)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
!      a(size(a,1)-staggered_ehalo+1:size(a,1), &
!          1+staggered_shalo:size(a,2)-staggered_nhalo) = &
!          erecv(:,1:size(a,2)-staggered_shalo-staggered_nhalo)
      a(size(a,1)-staggered_uhalo+1:size(a,1), &
          1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
          erecv(:,1:size(a,2)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
!      a(1:staggered_whalo, &
!          1+staggered_shalo:size(a,2)-staggered_nhalo) = &
!          wrecv(:,1:size(a,2)-staggered_shalo-staggered_nhalo)
      a(1:staggered_lhalo, &
          1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
          wrecv(:,1:size(a,2)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > south .or. fill_global_halos) then
      ssend(:,:) = &
!        a(:,1+staggered_shalo:1+staggered_shalo+staggered_nhalo-1)
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      nsend(:,:) = &
!        a(:,size(a,2)-staggered_nhalo-staggered_shalo+1:size(a,2)-staggered_nhalo)
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
!      a(:,size(a,2)-staggered_nhalo+1:size(a,2)) = nrecv(:,:)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2)) = nrecv(:,:)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
!      a(:,1:staggered_shalo) = srecv(:,:)
      a(:,1:staggered_lhalo) = srecv(:,:)
    endif

  end subroutine staggered_parallel_halo_real8_2d

  !WHL - Edited the original subroutine so that values from N and E edges
  !      of global domain can be written to halo cells at the S and W edges,
  !      to allow cyclic BCs for staggered variables

  subroutine staggered_parallel_halo_real8_3d(a)

    use mpi_mod
    implicit none
    real(8),dimension(:,:,:) :: a

    ! Implements a staggered grid halo update for a 3D field.
    ! As the grid is staggered, the array 'a' is one smaller in both dimensions than an unstaggered array.
    ! The vertical dimension is assumed to be the first index, i.e., a(k,i,j).

    ! The grid is laid out from the SW, and the lower left corner is assigned to this_rank = 0.
    ! It's eastern nbhr is task_id = 1, proceeding rowwise and starting from the western edge.
    ! The South-most processes own one additional row of stagggered variables on the southern edge
    ! and have one less 'southern' halo row than other processes. Likewise, the West-most processes own one 
    ! additional column of staggered variables on the western edge and have one less 'western' halo column. 
    ! This is implemented by a modification to the staggered_lhalo value on these processes. 

    ! Maintaining global boundary conditions are not addressed within this routine (yet).

    ! integer :: erequest,ierror,one,nrequest,srequest,wrequest
    integer :: ierror,nrequest,srequest,erequest,wrequest

!    real(8),dimension(size(a,1),staggered_whalo,size(a,3)-staggered_shalo-staggered_nhalo) :: esend,wrecv
!    real(8),dimension(size(a,1),staggered_ehalo,size(a,3)-staggered_shalo-staggered_nhalo) :: erecv,wsend
!    real(8),dimension(size(a,1),size(a,2),staggered_shalo) :: nsend,srecv
!    real(8),dimension(size(a,1),size(a,2),staggered_nhalo) :: nrecv,ssend
    real(8),dimension(size(a,1),staggered_lhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: esend,wrecv
    real(8),dimension(size(a,1),staggered_uhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: erecv,wsend
    real(8),dimension(size(a,1),size(a,2),staggered_lhalo) :: nsend,srecv
    real(8),dimension(size(a,1),size(a,2),staggered_uhalo) :: nrecv,ssend

    !WHL - I defined a logical variable to determine whether or not to fill halo cells
    !      at the edge of the global domain.  I am setting it to true by default to support
    !      cyclic global BCs.

    !TODO - Set to true in all cases? (Or simply remove the 'if's?)  

    logical :: fill_global_halos = .true.

    ! begin

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1 .or. size(a,3)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. fill_global_halos) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
!      wsend(:,:,1:size(a,3)-staggered_shalo-staggered_nhalo) = &
!        a(:,1+staggered_whalo:1+staggered_whalo+staggered_ehalo-1, &
!            1+staggered_shalo:size(a,3)-staggered_nhalo)
      wsend(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
            1+staggered_lhalo:size(a,3)-staggered_uhalo)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
!      esend(:,:,1:size(a,3)-staggered_shalo-staggered_nhalo) = &
!        a(:,size(a,2)-staggered_ehalo-staggered_whalo+1:size(a,2)-staggered_ehalo, &
!            1+staggered_shalo:size(a,3)-staggered_nhalo)
      esend(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
            1+staggered_lhalo:size(a,3)-staggered_uhalo)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
!      a(:,size(a,2)-staggered_ehalo+1:size(a,2), &
!          1+staggered_shalo:size(a,3)-staggered_nhalo) = &
!          erecv(:,:,1:size(a,3)-staggered_shalo-staggered_nhalo)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2), &
          1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
          erecv(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
!      a(:,1:staggered_whalo, &
!          1+staggered_shalo:size(a,3)-staggered_nhalo) = &
!          wrecv(:,:,1:size(a,3)-staggered_shalo-staggered_nhalo)
      a(:,1:staggered_lhalo, &
          1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
          wrecv(:,:,1:size(a,3)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > south .or. fill_global_halos) then
!      ssend(:,:,:) = &
!        a(:,:,1+staggered_shalo:1+staggered_shalo+staggered_nhalo-1)
      ssend(:,:,:) = &
        a(:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      nsend(:,:,:) = &
!        a(:,:,size(a,3)-staggered_nhalo-staggered_shalo+1:size(a,3)-staggered_nhalo)
        a(:,:,size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
!      a(:,:,size(a,3)-staggered_nhalo+1:size(a,3)) = nrecv(:,:,:)
      a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = nrecv(:,:,:)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
!      a(:,:,1:staggered_shalo) = srecv(:,:,:)
      a(:,:,1:staggered_lhalo) = srecv(:,:,:)
    endif

  end subroutine staggered_parallel_halo_real8_3d

  subroutine staggered_parallel_halo_real8_6d(a)

    use mpi_mod
    implicit none
    real(8),dimension(-1:,-1:,-1:,:,:,:) :: a

    ! Implements a staggered grid halo update for a 6D field.
    ! This subroutine is custom-made for the 6D arrays that hold matrix entries.

    ! As the grid is staggered, the array 'a' is one smaller in both dimensions than an unstaggered array.
    ! The vertical dimension is assumed to precede the i and j indices, i.e., a(:,:,:,k,i,j).

    ! NOTE: The first three dimensions are -1:1.
    ! The subroutine is specifically designed for matrix arrays with this structure.

    ! The grid is laid out from the SW, and the lower left corner is assigned to this_rank = 0.
    ! It's eastern nbhr is task_id = 1, proceeding rowwise and starting from the western edge.
    ! The South-most processes own one additional row of stagggered variables on the southern edge
    ! and have one less 'southern' halo row than other processes. Likewise, the West-most processes own one 
    ! additional column of staggered variables on the western edge and have one less 'western' halo column. 
    ! This is implemented by a modification to the staggered_lhalo value on these processes. 

    ! Maintaining global boundary conditions are not addressed within this routine (yet).

    ! integer :: erequest,ierror,one,nrequest,srequest,wrequest
    integer :: ierror,nrequest,srequest,erequest,wrequest

    real(8),dimension(-1:1,-1:1,-1:1,size(a,4), &
                      staggered_lhalo,size(a,6)-staggered_lhalo-staggered_uhalo) :: esend,wrecv
    real(8),dimension(-1:1,-1:1,-1:1,size(a,4), &
                      staggered_uhalo,size(a,6)-staggered_lhalo-staggered_uhalo) :: erecv,wsend
    real(8),dimension(-1:1,-1:1,-1:1,size(a,4),size(a,5),staggered_lhalo) :: nsend,srecv
    real(8),dimension(-1:1,-1:1,-1:1,size(a,4),size(a,5),staggered_uhalo) :: nrecv,ssend

    !WHL - I defined a logical variable to determine whether or not to fill halo cells
    !      at the edge of the global domain.  I am setting it to true by default to support
    !      cyclic global BCs.

    !TODO - Set to true in all cases? (Or simply remove the 'if's?)  

    logical :: fill_global_halos = .true.

    ! begin

    ! Confirm staggered array
    if (size(a,5)/=local_ewn-1 .or. size(a,6)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. fill_global_halos) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. fill_global_halos) then
      wsend(:,:,:,:,:,1:size(a,6)-staggered_lhalo-staggered_uhalo) = &
        a(:,:,:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
                  1+staggered_lhalo:size(a,6)-staggered_uhalo)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
      esend(:,:,:,:,:,1:size(a,6)-staggered_lhalo-staggered_uhalo) = &
        a(:,:,:,:,size(a,5)-staggered_uhalo-staggered_lhalo+1:size(a,5)-staggered_uhalo, &
                  1+staggered_lhalo:size(a,6)-staggered_uhalo)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. fill_global_halos) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:,:,:,:,size(a,5)-staggered_uhalo+1:size(a,5), &
                1+staggered_lhalo:size(a,6)-staggered_uhalo) = &
          erecv(:,:,:,:,:,1:size(a,6)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > west .or. fill_global_halos) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:,:,:,:,1:staggered_lhalo, &
                1+staggered_lhalo:size(a,6)-staggered_uhalo) = &
          wrecv(:,:,:,:,:,1:size(a,6)-staggered_lhalo-staggered_uhalo)
    endif

    if (this_rank > south .or. fill_global_halos) then
      ssend(:,:,:,:,:,:) = &
        a(:,:,:,:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      nsend(:,:,:,:,:,:) = &
        a(:,:,:,:,:,size(a,6)-staggered_uhalo-staggered_lhalo+1:size(a,6)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. fill_global_halos) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,:,:,:,size(a,6)-staggered_uhalo+1:size(a,6)) = nrecv(:,:,:,:,:,:)
    endif

    if (this_rank > south .or. fill_global_halos) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,:,:,:,1:staggered_lhalo) = srecv(:,:,:,:,:,:)
    endif

  end subroutine staggered_parallel_halo_real8_6d

! Following routines imported from the Community Earth System Model
! (models/utils/mct/mpeu.m_FcComms.F90)
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gather_int - Gather an array of type integer
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em integer} 
! to the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gather
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
  subroutine fc_gather_int (sendbuf, sendcnt, sendtype, &
                            recvbuf, recvcnt, recvtype, &
                            root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      integer,               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer,               intent(in)  :: recvcnt
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      integer,               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   integer :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, i, count, displs
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               if (recvcnt > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  displs = p*recvcnt
                  call mpi_irecv ( recvbuf(displs+1), recvcnt, &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         displs = mytid*recvcnt
         do i=1,sendcnt
            recvbuf(displs+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

      endif

   else
 
      call mpi_gather (sendbuf, sendcnt, sendtype, &
                       recvbuf, recvcnt, recvtype, &
                       root, comm, ier)
   endif

   return
  end subroutine fc_gather_int

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_int - Gather an array of type integer
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em integer} 
! to the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_int (sendbuf, sendcnt, sendtype, &
                              recvbuf, recvcnts, displs, recvtype, &
                              root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      integer,               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      integer,               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   integer :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
        fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

     endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_int

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_real4 - Gather an array of type real*4
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em real*4} to
! the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_real4 (sendbuf, sendcnt, sendtype, &
                                recvbuf, recvcnts, displs, recvtype, &
                                root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      real(4),               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      real(4),               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   real(4) :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1.0
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

      endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_real4

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_real8 - Gather an array of type real*4
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em real*8} to
! the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_real8 (sendbuf, sendcnt, sendtype, &
                                recvbuf, recvcnts, displs, recvtype, &
                                root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      real(8),               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      real(8),               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   real(8) :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1.0
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

      endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_real8

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_log - Gather an array of type logical
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em logical} 
! to the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_log (sendbuf, sendcnt, sendtype, &
                              recvbuf, recvcnts, displs, recvtype, &
                              root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      logical,               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      logical,               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   logical :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = .true.
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
        fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

     endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_log

end module parallel
