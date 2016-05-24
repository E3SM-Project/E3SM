!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_horiz_bcs_parallel.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!WHL - Commenting out a bunch of references to 'comm' so the code will compile

!Stores routines to support the maintainance of the desired horizontal boundary conditions on the global domain.

!Currently, these routines only support closed BC's with a free-slip transverse velocity. Once these are tested,
!more will be added with namelist variables.

!Each routine has been tested with a field to make sure that at least, they aren't indexing out of bounds. Testing
!for correctness is to come soon.

module glimmer_horiz_bcs
  implicit none !If a variable isn't defined, throw a compiler error
  private       !Must declare all public routines in the header.

!WHL - Added a logical parameter to turn on and off the various subroutines in this module
!      If set to true, all subroutines execute as designed.
!      If false, subroutines do nothing and return.

  logical, parameter :: enabled_global_horiz_bcs = .false.
!  logical, parameter :: enabled_global_horiz_bcs = .true.

!WHL - This is a hack to get the code to compile, because it cannot find
!      these variables in parallel_mpi
  integer,parameter :: staggered_whalo = 2
  integer,parameter :: staggered_shalo = 2
  integer,parameter :: staggered_ehalo = 1
  integer,parameter :: staggered_nhalo = 1

  integer, parameter, public :: HORIZ_BCS_WALL_SLIP = 0
  integer, parameter, public :: HORIZ_BCS_CYCLIC = 1

  integer, parameter :: horiz_bcs_type_north = HORIZ_BCS_CYCLIC
  integer, parameter :: horiz_bcs_type_south = HORIZ_BCS_CYCLIC
  integer, parameter :: horiz_bcs_type_east  = HORIZ_BCS_CYCLIC
  integer, parameter :: horiz_bcs_type_west  = HORIZ_BCS_CYCLIC

  !Enforce boundary conditions for a variety of variable types
  public :: horiz_bcs_unstag_scalar   !Unstaggered scalar variables
  public :: horiz_bcs_stag_vector_ew  !Staggered ew-direction vector
  public :: horiz_bcs_stag_vector_ns  !Staggered ns-direction vector
  public :: horiz_bcs_stag_scalar   !Unstaggered scalar variables

  interface horiz_bcs_unstag_scalar
    module procedure horiz_bcs_unstag_scalar_logical_2d
    module procedure horiz_bcs_unstag_scalar_integer_2d
    module procedure horiz_bcs_unstag_scalar_real4_2d
    module procedure horiz_bcs_unstag_scalar_real8_2d
    module procedure horiz_bcs_unstag_scalar_real8_3d
  end interface

  interface horiz_bcs_stag_vector_ew
    module procedure horiz_bcs_stag_vector_ew_real8_2d
    module procedure horiz_bcs_stag_vector_ew_real8_3d
  end interface

  interface horiz_bcs_stag_vector_ns
    module procedure horiz_bcs_stag_vector_ns_real8_2d
    module procedure horiz_bcs_stag_vector_ns_real8_3d
  end interface

  interface horiz_bcs_stag_scalar
    module procedure horiz_bcs_stag_scalar_integer_2d
    module procedure horiz_bcs_stag_scalar_real8_2d
    module procedure horiz_bcs_stag_scalar_real8_3d
  end interface

  integer, parameter, public :: ghost_shift = 0


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, lhalo, uhalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !! , comm
    use mpi_mod 
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    if ( nsub > global_nsn ) then   !I am on the north boundary
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(:,lhalo+own_nsn+i) = a(:,lhalo+own_nsn+1-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),uhalo))
          allocate(sbuf_n(size(a,1),lhalo))
          rbuf_n = a(:,lhalo+own_nsn+1      -ghost_shift:lhalo+own_nsn+uhalo-ghost_shift)
          sbuf_n = a(:,lhalo+own_nsn+1-lhalo-ghost_shift:lhalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then    !I am on the east boundary
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(lhalo+own_ewn+i,:) = a(lhalo+own_ewn+1-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(uhalo,size(a,2)))
          allocate(sbuf_e(lhalo,size(a,2)))
          rbuf_e = a(lhalo+own_ewn+1      -ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:)
          sbuf_e = a(lhalo+own_ewn+1-lhalo-ghost_shift:lhalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then   !I am on the south boundary
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(:,lhalo+1-i) = a(:,lhalo+i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),lhalo))
          allocate(sbuf_s(size(a,1),uhalo))
          rbuf_s = a(:,1      +ghost_shift:lhalo      +ghost_shift)
          sbuf_s = a(:,lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then    !I am on the west boundary
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(lhalo+1-i,:) = a(lhalo+i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(lhalo,size(a,2)))
          allocate(sbuf_w(uhalo,size(a,2)))
          rbuf_w = a(1      +ghost_shift:lhalo      +ghost_shift,:)
          sbuf_w = a(lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,lhalo+own_nsn+1-ghost_shift:lhalo+own_nsn+uhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(lhalo+own_ewn+1-ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:lhalo+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:lhalo+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif
  end subroutine horiz_bcs_unstag_scalar_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_integer_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, lhalo, uhalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !!, comm
    use mpi_mod 
    implicit none
    integer,dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    integer,dimension(:,:), allocatable :: sbuf_n, rbuf_n
    integer,dimension(:,:), allocatable :: sbuf_e, rbuf_e
    integer,dimension(:,:), allocatable :: sbuf_s, rbuf_s
    integer,dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    if ( nsub > global_nsn ) then   !I am on the north boundary
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(:,lhalo+own_nsn+i) = a(:,lhalo+own_nsn+1-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),uhalo))
          allocate(sbuf_n(size(a,1),lhalo))
          rbuf_n = a(:,lhalo+own_nsn+1      -ghost_shift:lhalo+own_nsn+uhalo-ghost_shift)
          sbuf_n = a(:,lhalo+own_nsn+1-lhalo-ghost_shift:lhalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_integer , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_integer , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then    !I am on the east boundary
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(lhalo+own_ewn+i,:) = a(lhalo+own_ewn+1-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(uhalo,size(a,2)))
          allocate(sbuf_e(lhalo,size(a,2)))
          rbuf_e = a(lhalo+own_ewn+1      -ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:)
          sbuf_e = a(lhalo+own_ewn+1-lhalo-ghost_shift:lhalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_integer , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_integer , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then   !I am on the south boundary
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(:,lhalo+1-i) = a(:,lhalo+i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),lhalo))
          allocate(sbuf_s(size(a,1),uhalo))
          rbuf_s = a(:,1      +ghost_shift:lhalo      +ghost_shift)
          sbuf_s = a(:,lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_integer , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_integer , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then    !I am on the west boundary
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(lhalo+1-i,:) = a(lhalo+i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(lhalo,size(a,2)))
          allocate(sbuf_w(uhalo,size(a,2)))
          rbuf_w = a(1      +ghost_shift:lhalo      +ghost_shift,:)
          sbuf_w = a(lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_integer , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_integer , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,lhalo+own_nsn+1-ghost_shift:lhalo+own_nsn+uhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(lhalo+own_ewn+1-ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:lhalo+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:lhalo+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif
  end subroutine horiz_bcs_unstag_scalar_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_logical_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, lhalo, uhalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !!, comm
    use mpi_mod 
    implicit none
    logical,dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    logical,dimension(:,:), allocatable :: sbuf_n, rbuf_n
    logical,dimension(:,:), allocatable :: sbuf_e, rbuf_e
    logical,dimension(:,:), allocatable :: sbuf_s, rbuf_s
    logical,dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    if ( nsub > global_nsn ) then   !I am on the north boundary
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(:,lhalo+own_nsn+i) = a(:,lhalo+own_nsn+1-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),uhalo))
          allocate(sbuf_n(size(a,1),lhalo))
          rbuf_n = a(:,lhalo+own_nsn+1      -ghost_shift:lhalo+own_nsn+uhalo-ghost_shift)
          sbuf_n = a(:,lhalo+own_nsn+1-lhalo-ghost_shift:lhalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_logical , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_logical , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then    !I am on the east boundary
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(lhalo+own_ewn+i,:) = a(lhalo+own_ewn+1-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(uhalo,size(a,2)))
          allocate(sbuf_e(lhalo,size(a,2)))
          rbuf_e = a(lhalo+own_ewn+1      -ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:)
          sbuf_e = a(lhalo+own_ewn+1-lhalo-ghost_shift:lhalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_logical , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_logical , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then   !I am on the south boundary
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(:,lhalo+1-i) = a(:,lhalo+i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),lhalo))
          allocate(sbuf_s(size(a,1),uhalo))
          rbuf_s = a(:,1      +ghost_shift:lhalo      +ghost_shift)
          sbuf_s = a(:,lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_logical , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_logical , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then    !I am on the west boundary
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(lhalo+1-i,:) = a(lhalo+i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(lhalo,size(a,2)))
          allocate(sbuf_w(uhalo,size(a,2)))
          rbuf_w = a(1      +ghost_shift:lhalo      +ghost_shift,:)
          sbuf_w = a(lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_logical , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_logical , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,lhalo+own_nsn+1-ghost_shift:lhalo+own_nsn+uhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(lhalo+own_ewn+1-ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:lhalo+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:lhalo+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif
  end subroutine horiz_bcs_unstag_scalar_logical_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real4_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, lhalo, uhalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !!, comm
    use mpi_mod 
    implicit none
    real(4),dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(4),dimension(:,:), allocatable :: sbuf_n, rbuf_n
    real(4),dimension(:,:), allocatable :: sbuf_e, rbuf_e
    real(4),dimension(:,:), allocatable :: sbuf_s, rbuf_s
    real(4),dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    if ( nsub > global_nsn ) then   !I am on the north boundary
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(:,lhalo+own_nsn+i) = a(:,lhalo+own_nsn+1-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),uhalo))
          allocate(sbuf_n(size(a,1),lhalo))
          rbuf_n = a(:,lhalo+own_nsn+1      -ghost_shift:lhalo+own_nsn+uhalo-ghost_shift)
          sbuf_n = a(:,lhalo+own_nsn+1-lhalo-ghost_shift:lhalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then    !I am on the east boundary
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(lhalo+own_ewn+i,:) = a(lhalo+own_ewn+1-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(uhalo,size(a,2)))
          allocate(sbuf_e(lhalo,size(a,2)))
          rbuf_e = a(lhalo+own_ewn+1      -ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:)
          sbuf_e = a(lhalo+own_ewn+1-lhalo-ghost_shift:lhalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then   !I am on the south boundary
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(:,lhalo+1-i) = a(:,lhalo+i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),lhalo))
          allocate(sbuf_s(size(a,1),uhalo))
          rbuf_s = a(:,1      +ghost_shift:lhalo      +ghost_shift)
          sbuf_s = a(:,lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then    !I am on the west boundary
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(lhalo+1-i,:) = a(lhalo+i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(lhalo,size(a,2)))
          allocate(sbuf_w(uhalo,size(a,2)))
          rbuf_w = a(1      +ghost_shift:lhalo      +ghost_shift,:)
          sbuf_w = a(lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,lhalo+own_nsn+1-ghost_shift:lhalo+own_nsn+uhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(lhalo+own_ewn+1-ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:lhalo+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:lhalo+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif
  end subroutine horiz_bcs_unstag_scalar_real4_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_3d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, lhalo, uhalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !!, comm
    use mpi_mod 
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    if ( nsub > global_nsn ) then   !I am on the north boundary
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(:,:,lhalo+own_nsn+i) = a(:,:,lhalo+own_nsn+1-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),size(a,2),uhalo))
          allocate(sbuf_n(size(a,1),size(a,2),lhalo))
          rbuf_n = a(:,:,lhalo+own_nsn+1      -ghost_shift:lhalo+own_nsn+uhalo-ghost_shift)
          sbuf_n = a(:,:,lhalo+own_nsn+1-lhalo-ghost_shift:lhalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then    !I am on the east boundary
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , uhalo
            a(:,lhalo+own_ewn+i,:) = a(:,lhalo+own_ewn+1-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(size(a,1),uhalo,size(a,3)))
          allocate(sbuf_e(size(a,1),lhalo,size(a,3)))
          rbuf_e = a(:,lhalo+own_ewn+1      -ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:)
          sbuf_e = a(:,lhalo+own_ewn+1-lhalo-ghost_shift:lhalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then   !I am on the south boundary
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(:,:,lhalo+1-i) = a(:,:,lhalo+i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),size(a,2),lhalo))
          allocate(sbuf_s(size(a,1),size(a,2),uhalo))
          rbuf_s = a(:,:,1      +ghost_shift:lhalo      +ghost_shift)
          sbuf_s = a(:,:,lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then    !I am on the west boundary
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , lhalo
            a(:,lhalo+1-i,:) = a(:,lhalo+i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(size(a,1),lhalo,size(a,3)))
          allocate(sbuf_w(size(a,1),uhalo,size(a,3)))
          rbuf_w = a(:,1      +ghost_shift:lhalo      +ghost_shift,:)
          sbuf_w = a(:,lhalo+1+ghost_shift:lhalo+uhalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,:,lhalo+own_nsn+1-ghost_shift:lhalo+own_nsn+uhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(:,lhalo+own_ewn+1-ghost_shift:lhalo+own_ewn+uhalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,:,1+ghost_shift:lhalo+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:lhalo+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif
  end subroutine horiz_bcs_unstag_scalar_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ew_real8_3d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, &
!!                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !!, comm
    use mpi_mod
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( nsub > global_nsn ) then
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , nhalo
            a(:,:,shalo+ns_npts+i) = a(:,:,shalo+ns_npts-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),size(a,2),nhalo  ))
          allocate(sbuf_n(size(a,1),size(a,2),shalo+1))
          rbuf_n = a(:,:,shalo+own_nsn+1          -ghost_shift:shalo+own_nsn+nhalo-ghost_shift)
          sbuf_n = a(:,:,shalo+own_nsn+1-(shalo+1)-ghost_shift:shalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , ehalo
            a(:,whalo+ew_npts+i,:) = -a(:,whalo+ew_npts-i,:)
          enddo
          !Normal velocities are zero at boundaries
          a(:,whalo+ew_npts,:) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(size(a,1),ehalo  ,size(a,3)))
          allocate(sbuf_e(size(a,1),whalo+1,size(a,3)))
          rbuf_e = a(:,whalo+own_ewn+1          -ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:)
          sbuf_e = a(:,whalo+own_ewn+1-(whalo+1)-ghost_shift:whalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , shalo
            a(:,:,shalo+1-i) = a(:,:,shalo+1+i)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !southern boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: shalo >= 1 & nsn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using shalo, shalo+2, and shalo+3
          a(:,:,shalo+1) = a(:,:,shalo) / 3.D0 + a(:,:,shalo+2) - a(:,:,shalo+3) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),size(a,2),shalo+1))
          allocate(sbuf_s(size(a,1),size(a,2),nhalo  ))
          rbuf_s = a(:,:,1      +ghost_shift:shalo+1      +ghost_shift)
          sbuf_s = a(:,:,shalo+2+ghost_shift:shalo+1+nhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , whalo
            a(:,whalo+1-i,:) = -a(:,whalo+1+i,:)
          enddo
          !Normal velocities are zero at boundaries
          a(:,whalo+1      ,:) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(size(a,1),whalo+1,size(a,3)))
          allocate(sbuf_w(size(a,1),ehalo  ,size(a,3)))
          rbuf_w = a(:,1      +ghost_shift:whalo+1      +ghost_shift,:)
          sbuf_w = a(:,whalo+2+ghost_shift:whalo+1+ehalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,:,shalo+own_nsn+1-ghost_shift:shalo+own_nsn+nhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(:,whalo+own_ewn+1-ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,:,1+ghost_shift:shalo+1+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:whalo+1+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif

  end subroutine horiz_bcs_stag_vector_ew_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ew_real8_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, &
!!                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !!, comm
    use mpi_mod
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( nsub > global_nsn ) then
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , nhalo
            a(:,shalo+ns_npts+i) = a(:,shalo+ns_npts-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),nhalo  ))
          allocate(sbuf_n(size(a,1),shalo+1))
          rbuf_n = a(:,shalo+own_nsn+1          -ghost_shift:shalo+own_nsn+nhalo-ghost_shift)
          sbuf_n = a(:,shalo+own_nsn+1-(shalo+1)-ghost_shift:shalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , ehalo
            a(whalo+ew_npts+i,:) = -a(whalo+ew_npts-i,:)
          enddo
          !Normal velocities are zero at boundaries
          a(whalo+ew_npts,:) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(ehalo  ,size(a,2)))
          allocate(sbuf_e(whalo+1,size(a,2)))
          rbuf_e = a(whalo+own_ewn+1          -ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:)
          sbuf_e = a(whalo+own_ewn+1-(whalo+1)-ghost_shift:whalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , shalo
            a(:,shalo+1-i) = a(:,shalo+1+i)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !southern boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: shalo >= 1 & nsn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using shalo, shalo+2, and shalo+3
          a(:,shalo+1) = a(:,shalo) / 3.D0 + a(:,shalo+2) - a(:,shalo+3) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),shalo+1))
          allocate(sbuf_s(size(a,1),nhalo  ))
          rbuf_s = a(:,1      +ghost_shift:shalo+1      +ghost_shift)
          sbuf_s = a(:,shalo+2+ghost_shift:shalo+1+nhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , whalo
            a(whalo+1-i,:) = -a(whalo+1+i,:)
          enddo
          !Normal velocities are zero at boundaries
          a(whalo+1      ,:) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(whalo+1,size(a,2)))
          allocate(sbuf_w(ehalo  ,size(a,2)))
          rbuf_w = a(1      +ghost_shift:whalo+1      +ghost_shift,:)
          sbuf_w = a(whalo+2+ghost_shift:whalo+1+ehalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,shalo+own_nsn+1-ghost_shift:shalo+own_nsn+nhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(whalo+own_ewn+1-ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:shalo+1+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:whalo+1+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif

  end subroutine horiz_bcs_stag_vector_ew_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ns_real8_3d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, &
!!                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks   !!, comm
    use mpi_mod
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( nsub > global_nsn ) then
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , nhalo
            a(:,:,shalo+ns_npts+i) = -a(:,:,shalo+ns_npts-i)
          enddo
          !Normal velocities are zero at boundaries
          a(:,:,shalo+ns_npts) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),size(a,2),nhalo  ))
          allocate(sbuf_n(size(a,1),size(a,2),shalo+1))
          rbuf_n = a(:,:,shalo+own_nsn+1          -ghost_shift:shalo+own_nsn+nhalo-ghost_shift)
          sbuf_n = a(:,:,shalo+own_nsn+1-(shalo+1)-ghost_shift:shalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , ehalo
            a(:,whalo+ew_npts+i,:) = a(:,whalo+ew_npts-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(size(a,1),ehalo  ,size(a,3)))
          allocate(sbuf_e(size(a,1),whalo+1,size(a,3)))
          rbuf_e = a(:,whalo+own_ewn+1          -ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:)
          sbuf_e = a(:,whalo+own_ewn+1-(whalo+1)-ghost_shift:whalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , shalo
            a(:,:,shalo+1-i) = -a(:,:,shalo+1+i)
          enddo
          !Normal velocities are zero at boundaries
          a(:,:,shalo+1      ) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),size(a,2),shalo+1))
          allocate(sbuf_s(size(a,1),size(a,2),nhalo  ))
          rbuf_s = a(:,:,1      +ghost_shift:shalo+1      +ghost_shift)
          sbuf_s = a(:,:,shalo+2+ghost_shift:shalo+1+nhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , whalo
            a(:,whalo+1-i,:) = a(:,whalo+1+i,:)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(:,whalo+1,:) = a(:,whalo,:) / 3.D0 + a(:,whalo+2,:) - a(:,whalo+3,:) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(size(a,1),whalo+1,size(a,3)))
          allocate(sbuf_w(size(a,1),ehalo  ,size(a,3)))
          rbuf_w = a(:,1      +ghost_shift:whalo+1      +ghost_shift,:)
          sbuf_w = a(:,whalo+2+ghost_shift:whalo+1+ehalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,:,shalo+own_nsn+1-ghost_shift:shalo+own_nsn+nhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(:,whalo+own_ewn+1-ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,:,1+ghost_shift:shalo+1+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:whalo+1+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif

  end subroutine horiz_bcs_stag_vector_ns_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ns_real8_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, &
!!                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks   !!, comm
    use mpi_mod
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( nsub > global_nsn ) then
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , nhalo
            a(:,shalo+ns_npts+i) = -a(:,shalo+ns_npts-i)
          enddo
          !Normal velocities are zero at boundaries
          a(:,shalo+ns_npts) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),nhalo  ))
          allocate(sbuf_n(size(a,1),shalo+1))
          rbuf_n = a(:,shalo+own_nsn+1          -ghost_shift:shalo+own_nsn+nhalo-ghost_shift)
          sbuf_n = a(:,shalo+own_nsn+1-(shalo+1)-ghost_shift:shalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , ehalo
            a(whalo+ew_npts+i,:) = a(whalo+ew_npts-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(ehalo  ,size(a,2)))
          allocate(sbuf_e(whalo+1,size(a,2)))
          rbuf_e = a(whalo+own_ewn+1          -ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:)
          sbuf_e = a(whalo+own_ewn+1-(whalo+1)-ghost_shift:whalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , shalo
            a(:,shalo+1-i) = -a(:,shalo+1+i)
          enddo
          !Normal velocities are zero at boundaries
          a(:,shalo+1      ) = 0.D0
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),shalo+1))
          allocate(sbuf_s(size(a,1),nhalo  ))
          rbuf_s = a(:,1      +ghost_shift:shalo+1      +ghost_shift)
          sbuf_s = a(:,shalo+2+ghost_shift:shalo+1+nhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , whalo
            a(whalo+1-i,:) = a(whalo+1+i,:)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(whalo+1,:) = a(whalo,:) / 3.D0 + a(whalo+2,:) - a(whalo+3,:) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(whalo+1,size(a,2)))
          allocate(sbuf_w(ehalo  ,size(a,2)))
          rbuf_w = a(1      +ghost_shift:whalo+1      +ghost_shift,:)
          sbuf_w = a(whalo+2+ghost_shift:whalo+1+ehalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,shalo+own_nsn+1-ghost_shift:shalo+own_nsn+nhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(whalo+own_ewn+1-ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:shalo+1+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:whalo+1+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif

  end subroutine horiz_bcs_stag_vector_ns_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_integer_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, &
!!                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks  !!, comm
    use mpi_mod
    implicit none
    integer,dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    integer,dimension(:,:), allocatable :: sbuf_n, rbuf_n
    integer,dimension(:,:), allocatable :: sbuf_e, rbuf_e
    integer,dimension(:,:), allocatable :: sbuf_s, rbuf_s
    integer,dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( nsub > global_nsn ) then
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , nhalo
            a(:,shalo+ns_npts+i) = a(:,shalo+ns_npts-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),nhalo  ))
          allocate(sbuf_n(size(a,1),shalo+1))
          rbuf_n = a(:,shalo+own_nsn+1          -ghost_shift:shalo+own_nsn+nhalo-ghost_shift)
          sbuf_n = a(:,shalo+own_nsn+1-(shalo+1)-ghost_shift:shalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_integer , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_integer , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , ehalo
            a(whalo+ew_npts+i,:) = a(whalo+ew_npts-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(ehalo  ,size(a,2)))
          allocate(sbuf_e(whalo+1,size(a,2)))
          rbuf_e = a(whalo+own_ewn+1          -ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:)
          sbuf_e = a(whalo+own_ewn+1-(whalo+1)-ghost_shift:whalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_integer , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_integer , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , shalo
            a(:,shalo+1-i) = a(:,shalo+1+i)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(:,shalo+1) = a(:,shalo) / 3.D0 + a(:,shalo+2) - a(:,shalo+3) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),shalo+1))
          allocate(sbuf_s(size(a,1),nhalo  ))
          rbuf_s = a(:,1      +ghost_shift:shalo+1      +ghost_shift)
          sbuf_s = a(:,shalo+2+ghost_shift:shalo+1+nhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_integer , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_integer , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , whalo
            a(whalo+1-i,:) = a(whalo+1+i,:)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(whalo+1,:) = a(whalo,:) / 3.D0 + a(whalo+2,:) - a(whalo+3,:) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(whalo+1,size(a,2)))
          allocate(sbuf_w(ehalo  ,size(a,2)))
          rbuf_w = a(1      +ghost_shift:whalo+1      +ghost_shift,:)
          sbuf_w = a(whalo+2+ghost_shift:whalo+1+ehalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_integer , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_integer , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,shalo+own_nsn+1-ghost_shift:shalo+own_nsn+nhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(whalo+own_ewn+1-ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:shalo+1+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:whalo+1+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif

  end subroutine horiz_bcs_stag_scalar_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_real8_2d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, &
!!                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks   !!, comm
    use mpi_mod
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( nsub > global_nsn ) then
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , nhalo
            a(:,shalo+ns_npts+i) = a(:,shalo+ns_npts-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),nhalo  ))
          allocate(sbuf_n(size(a,1),shalo+1))
          rbuf_n = a(:,shalo+own_nsn+1          -ghost_shift:shalo+own_nsn+nhalo-ghost_shift)
          sbuf_n = a(:,shalo+own_nsn+1-(shalo+1)-ghost_shift:shalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , ehalo
            a(whalo+ew_npts+i,:) = a(whalo+ew_npts-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(ehalo  ,size(a,2)))
          allocate(sbuf_e(whalo+1,size(a,2)))
          rbuf_e = a(whalo+own_ewn+1          -ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:)
          sbuf_e = a(whalo+own_ewn+1-(whalo+1)-ghost_shift:whalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , shalo
            a(:,shalo+1-i) = a(:,shalo+1+i)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(:,shalo+1) = a(:,shalo) / 3.D0 + a(:,shalo+2) - a(:,shalo+3) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),shalo+1))
          allocate(sbuf_s(size(a,1),nhalo  ))
          rbuf_s = a(:,1      +ghost_shift:shalo+1      +ghost_shift)
          sbuf_s = a(:,shalo+2+ghost_shift:shalo+1+nhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , whalo
            a(whalo+1-i,:) = a(whalo+1+i,:)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(whalo+1,:) = a(whalo,:) / 3.D0 + a(whalo+2,:) - a(whalo+3,:) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(whalo+1,size(a,2)))
          allocate(sbuf_w(ehalo  ,size(a,2)))
          rbuf_w = a(1      +ghost_shift:whalo+1      +ghost_shift,:)
          sbuf_w = a(whalo+2+ghost_shift:whalo+1+ehalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,shalo+own_nsn+1-ghost_shift:shalo+own_nsn+nhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(whalo+own_ewn+1-ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:shalo+1+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(1+ghost_shift:whalo+1+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif

  end subroutine horiz_bcs_stag_scalar_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_real8_3d( a )
    use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, &
!!                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo, &
                        rank => this_rank, ewtasks => ProcsEW, tasks   !!, comm
    use mpi_mod
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    integer :: i, partner, nstasks, ierr
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    integer :: send_req_n, recv_req_n
    integer :: send_req_e, recv_req_e
    integer :: send_req_s, recv_req_s
    integer :: send_req_w, recv_req_w
    real(8),dimension(:,:,:), allocatable :: sbuf_n, rbuf_n
    real(8),dimension(:,:,:), allocatable :: sbuf_e, rbuf_e
    real(8),dimension(:,:,:), allocatable :: sbuf_s, rbuf_s
    real(8),dimension(:,:,:), allocatable :: sbuf_w, rbuf_w

!WHL - optional return
    if (.not. enabled_global_horiz_bcs) return

    nstasks = tasks / ewtasks
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( nsub > global_nsn ) then
      select case (horiz_bcs_type_north)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , nhalo
            a(:,:,shalo+ns_npts+i) = a(:,:,shalo+ns_npts-i)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = mod(rank,ewtasks)
          allocate(rbuf_n(size(a,1),size(a,2),nhalo  ))
          allocate(sbuf_n(size(a,1),size(a,2),shalo+1))
          rbuf_n = a(:,:,shalo+own_nsn+1          -ghost_shift:shalo+own_nsn+nhalo-ghost_shift)
          sbuf_n = a(:,:,shalo+own_nsn+1-(shalo+1)-ghost_shift:shalo+own_nsn      -ghost_shift)
!!          call mpi_irecv( rbuf_n , size( rbuf_n ) , mpi_real8 , partner , 1 , comm , recv_req_n , ierr )
!!          call mpi_isend( sbuf_n , size( sbuf_n ) , mpi_real8 , partner , 2 , comm , send_req_n , ierr )
      endselect
    endif
    if ( ewub > global_ewn ) then
      select case (horiz_bcs_type_east)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , ehalo
            a(:,whalo+ew_npts+i,:) = a(:,whalo+ew_npts-i,:)
          enddo
        case (HORIZ_BCS_CYCLIC)
          partner = rank - (ewtasks-1)
          allocate(rbuf_e(size(a,1),ehalo  ,size(a,3)))
          allocate(sbuf_e(size(a,1),whalo+1,size(a,3)))
          rbuf_e = a(:,whalo+own_ewn+1          -ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:)
          sbuf_e = a(:,whalo+own_ewn+1-(whalo+1)-ghost_shift:whalo+own_ewn      -ghost_shift,:)
!!          call mpi_irecv( rbuf_e , size( rbuf_e ) , mpi_real8 , partner , 3 , comm , recv_req_e , ierr )
!!          call mpi_isend( sbuf_e , size( sbuf_e ) , mpi_real8 , partner , 4 , comm , send_req_e , ierr )
      endselect
    endif
    if ( nslb < 1 ) then
      select case (horiz_bcs_type_south)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , shalo
            a(:,:,shalo+1-i) = a(:,:,shalo+1+i)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(:,:,shalo+1) = a(:,:,shalo) / 3.D0 + a(:,:,shalo+2) - a(:,:,shalo+3) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = (nstasks-1)*ewtasks+mod(rank,ewtasks)
          allocate(rbuf_s(size(a,1),size(a,2),shalo+1))
          allocate(sbuf_s(size(a,1),size(a,2),nhalo  ))
          rbuf_s = a(:,:,1      +ghost_shift:shalo+1      +ghost_shift)
          sbuf_s = a(:,:,shalo+2+ghost_shift:shalo+1+nhalo+ghost_shift)
!!          call mpi_irecv( rbuf_s , size( rbuf_s ) , mpi_real8 , partner , 2 , comm , recv_req_s , ierr )
!!          call mpi_isend( sbuf_s , size( sbuf_s ) , mpi_real8 , partner , 1 , comm , send_req_s , ierr )
      endselect
    endif
    if ( ewlb < 1 ) then
      select case (horiz_bcs_type_west)
        case (HORIZ_BCS_WALL_SLIP)
          do i = 1 , whalo
            a(:,whalo+1-i,:) = a(:,whalo+1+i,:)
          enddo
          !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
          !western boundary, however, and that must be interpolated.
          !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
          !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
          a(:,whalo+1,:) = a(:,whalo,:) / 3.D0 + a(:,whalo+2,:) - a(:,whalo+3,:) / 3.D0
        case (HORIZ_BCS_CYCLIC)
          partner = rank + (ewtasks-1)
          allocate(rbuf_w(size(a,1),whalo+1,size(a,3)))
          allocate(sbuf_w(size(a,1),ehalo  ,size(a,3)))
          rbuf_w = a(:,1      +ghost_shift:whalo+1      +ghost_shift,:)
          sbuf_w = a(:,whalo+2+ghost_shift:whalo+1+ehalo+ghost_shift,:)
!!          call mpi_irecv( rbuf_w , size( rbuf_w ) , mpi_real8 , partner , 4 , comm , recv_req_w , ierr )
!!          call mpi_isend( sbuf_w , size( sbuf_w ) , mpi_real8 , partner , 3 , comm , send_req_w , ierr )
      endselect
    endif

    if ( ( nsub > global_nsn ) .and. ( horiz_bcs_type_north == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_n , mpi_status_ignore , ierr )
      call mpi_wait( send_req_n , mpi_status_ignore , ierr )
      a(:,:,shalo+own_nsn+1-ghost_shift:shalo+own_nsn+nhalo-ghost_shift) = rbuf_n
      deallocate(rbuf_n)
      deallocate(sbuf_n)
    endif
    if ( ( ewub > global_ewn ) .and. ( horiz_bcs_type_east  == HORIZ_BCS_CYCLIC ) ) then   !I am on the north boundary
      call mpi_wait( recv_req_e , mpi_status_ignore , ierr )
      call mpi_wait( send_req_e , mpi_status_ignore , ierr )
      a(:,whalo+own_ewn+1-ghost_shift:whalo+own_ewn+ehalo-ghost_shift,:) = rbuf_e
      deallocate(rbuf_e)
      deallocate(sbuf_e)
    endif
    if ( ( nslb < 1          ) .and. ( horiz_bcs_type_south == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_s , mpi_status_ignore , ierr )
      call mpi_wait( send_req_s , mpi_status_ignore , ierr )
      a(:,:,1+ghost_shift:shalo+1+ghost_shift) = rbuf_s
      deallocate(rbuf_s)
      deallocate(sbuf_s)
    endif
    if ( ( ewlb < 1          ) .and. ( horiz_bcs_type_west  == HORIZ_BCS_CYCLIC ) ) then   !I am on the south boundary
      call mpi_wait( recv_req_w , mpi_status_ignore , ierr )
      call mpi_wait( send_req_w , mpi_status_ignore , ierr )
      a(:,1+ghost_shift:whalo+1+ghost_shift,:) = rbuf_w
      deallocate(rbuf_w)
      deallocate(sbuf_w)
    endif

  end subroutine horiz_bcs_stag_scalar_real8_3d

end module glimmer_horiz_bcs

