!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_test.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This module holds some test subroutines for the Glissade dynamical core
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glissade_test

  use glimmer_global, only: dp
  use glimmer_log
  use glide_types

  implicit none

  private
  public :: glissade_test_halo, glissade_test_transport

contains

!=======================================================================

  subroutine glissade_test_halo(model)

    use parallel

    ! various tests of parallel halo updates

    ! print statements are formatted for a 30x30 global array of scalars
    ! (34x34 with nhalo = 2), as for the standard dome problem

    type(glide_global_type), intent(inout) :: model      ! model instance

    integer, dimension (:,:), allocatable    ::  pgID       ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDr4     ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDr8     ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDr8_3d  ! unique global ID for parallel runs  

    integer, dimension (:,:), allocatable    ::  pgIDstagi  ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDstagr  ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDstagr3 ! unique global ID for parallel runs  

    logical,  dimension(:,:), allocatable   :: logvar
    integer,  dimension(:,:), allocatable   :: intvar
    real,     dimension(:,:), allocatable   :: r4var
    real(dp), dimension(:,:), allocatable   :: r8var
    real(dp), dimension(:,:,:), allocatable :: r8var_3d

    integer :: i, j, k
    integer :: nx, ny, nz

    integer, parameter :: rdiag = 0         ! rank for diagnostic prints 

    print*, ' '
    print*, 'In glissade_test_halo, this_rank =', this_rank

    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn

    allocate(logvar(nx,ny))
    allocate(intvar(nx,ny))
    allocate(r4var(nx,ny))
    allocate(r8var(nx,ny))
    allocate(r8var_3d(nz,nx,ny))

    allocate(pgID(nx,ny))
    allocate(pgIDr4(nx,ny))
    allocate(pgIDr8(nx,ny))
    allocate(pgIDr8_3d(nz,nx,ny))
    allocate(pgIDstagi(nx-1,ny-1))
    allocate(pgIDstagr(nx-1,ny-1))
    allocate(pgIDstagr3(nz,nx-1,ny-1))

    if (main_task) then
       print*, ' '
       print*, 'nx, ny, nz =', nx, ny, nz
       print*, 'uhalo, lhalo =', uhalo, lhalo
       print*, 'global_ewn, global_nsn =', global_ewn, global_nsn
       print*, ' '
    endif

    print*, 'this_rank, global_row/col offset =', this_rank, global_row_offset, global_col_offset

    ! Test the 5 standard parallel_halo routines for scalars: logical_2d, integer_2d, real4_2d, real8_2d, real8_3d

    ! logical 2D field

    logvar(:,:) = .false.

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       logvar(i,j) = .true.
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'logvar: this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34l3)') logvar(:,j)
       enddo
    endif

    call parallel_halo(logvar)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34l3)') logvar(:,j)
       enddo
    endif

    ! The next part of the code concerns parallel global IDs

    ! Compute parallel global ID for each grid cell

    pgID(:,:) = 0   ! integer

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgID(i,j) = parallel_globalID_scalar(i,j,nz)    ! function in parallel_mpi.F90
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (integer), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34i5)') pgID(:,j)
       enddo
    endif

    call parallel_halo(pgID)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34i5)') pgID(:,j)
       enddo
    endif

    ! real 2D
    
    pgIDr4(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr4(i,j) = real(parallel_globalID_scalar(i,j,nz))
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (r4 2D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr4(:,j)
       enddo
    endif

    call parallel_halo(pgIDr4)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr4(:,j)
       enddo
    endif

    ! double 2D
    
    pgIDr8(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr8(i,j) = real(parallel_globalID_scalar(i,j,nz), dp)
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (r8 2D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr8(:,j)
       enddo
    endif

    call parallel_halo(pgIDr8)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr8(:,j)
       enddo
    endif

    ! double 3D

    pgIDr8_3d(:,:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          pgIDr8_3d(k,i,j) = real(parallel_globalID_scalar(i,j,nz),dp) + real(k,dp)    ! function in parallel_mpi.F90
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (real 3D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr8_3d(k,:,j)
       enddo
    endif

    call parallel_halo(pgIDr8_3d)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr8_3d(k,:,j)
       enddo
    endif

    ! Repeat for staggered variables

    ! First for an integer 2D field

    pgIDstagi(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagi(i,j) = parallel_globalID_scalar(i,j,nz)    ! function in parallel_mpi.F90
    enddo
    enddo

    ! Print
    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (integer), this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33i5)') pgIDstagi(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagi)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33i5)') pgIDstagi(:,j)
       enddo
    endif

    ! Then for a real 2D field

    pgIDstagr(:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagr(i,j) = real(parallel_globalID_scalar(i,j,nz),dp)    ! function in parallel_mpi.F90
    enddo
    enddo

    ! Print
    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (real 2D), this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr(:,j)
       enddo
    endif

    ! Then for a real 3D field

    pgIDstagr3(:,:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          pgIDstagr3(k,i,j) = real(parallel_globalID_scalar(i,j,nz),dp) + real(k,dp)    ! function in parallel_mpi.F90
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (real 3D), k, this_rank =', k, this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr3)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    deallocate(logvar)
    deallocate(intvar)
    deallocate(r4var)
    deallocate(r8var)
    deallocate(r8var_3d)
    deallocate(pgID)
    deallocate(pgIDr4)
    deallocate(pgIDr8)
    deallocate(pgIDr8_3d)
    deallocate(pgIDstagi)
    deallocate(pgIDstagr)
    deallocate(pgIDstagr3)

  end subroutine glissade_test_halo

!=======================================================================

  subroutine glissade_test_transport(model)

    use parallel
    use glissade_transport, only: glissade_transport_driver, &
         glissade_transport_setup_tracers, glissade_transport_finish_tracers
    use glimmer_paramets, only: len0, thk0, tim0
    use glimmer_physcon, only: pi, scyr

    use glide_diagnostics

    !-------------------------------------------------------------------
    ! Test transport of a cylinder or block of ice once around the domain and
    ! back to the starting point, assuming uniform motion in a straight line.
    !
    ! Instructions for use:
    ! (1) At the top of this module, set test_transport = .true. and choose a
    !     value for thk_init.
    ! (2) Set the config file to run with the glissade dycore for one timestep, 
    !     with CF output written every timestep. 
    !     Note: Whatever the initial ice geometry in the config file
    !     (e.g., the dome test case), the ice extent will be preserved, but
    !     the thickness will be set everywhere to thk_init, giving a steep
    !     front at the ice margin that is challenging for transport schemes.
    ! (3) Comment out the call to glissade_diagnostic_variable_solve in
    !     cism_driver or simple_glide.  (It probable won't hurt to
    !     have it turned on, but will just use extra cycles.)
    !
    ! During the run, the following will happen:
    ! (1) During glissade_initialise, the ice thickness will be set to
    !     thk_init (everywhere there is ice).
    ! (2) During the first call to glissade_tstep, this subroutine 
    !     (glissade_test_transport) will be called.  The ice will be transported 
    !     around to its initial starting point (assuming periodic BCs).  
    ! (3) Then the model will return from glissade_tstep before doing any 
    !     other computations. Output will be written to netCDF and the code 
    !     will complete.
    !
    ! There should be two time slices in the output netCDF file.
    ! Compare the ice geometry at these two slices to see how diffusive 
    !  and/or dispersive the transport scheme is. A perfect scheme
    !  would leave the geometry unchanged.  
    !-------------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model      ! model instance

    integer :: ntracers   ! number of tracers to be transported

    real(dp), dimension(:,:,:), allocatable :: uvel, vvel   ! uniform velocity field (m/yr)

    integer :: i, j, k, n
    integer :: nx, ny, nz
    real(dp) :: dx, dy

    integer, parameter :: rdiag = 0         ! rank for diagnostic prints 

    real(dp), parameter :: umag = 100.     ! uniform speed (m/yr)

    ! Set angle of motion
    !WHL - Tested all of these angles (eastward, northward, and northeastward)
    real(dp), parameter :: theta = 0.d0      ! eastward
!    real(dp), parameter :: theta = pi/2.d0   ! northward
!    real(dp), parameter :: theta = pi/4.d0   ! northeastward

    real(dp), parameter :: thk = 500.d0

    real(dp) :: dt          ! time step in yr

    integer :: ntstep       ! run for this number of timesteps
    real(dp) :: theta_c    ! angle dividing paths through EW walls from paths thru NS walls
    real(dp) :: lenx       ! length of shortest path thru EW walls
    real(dp) :: leny       ! length of shortest path thru NS walls
    real(dp) :: len_path   ! length of path back to starting point
    real(dp) :: adv_cfl    ! advective CFL number

    logical :: do_upwind_transport  ! if true, do upwind transport

    ! Initialize

    dx = model%numerics%dew * len0
    dy = model%numerics%dns * len0

    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn

    allocate(uvel(nz,nx-1,ny-1))
    allocate(vvel(nz,nx-1,ny-1))

    ! Find the length of the path around the domain and back to the starting point

    lenx = global_ewn * dx
    leny = global_nsn * dy
    theta_c = atan(leny/lenx)   ! 0 <= theta_c <= pi/2

    if ( (theta >= -theta_c   .and. theta <= theta_c) .or.   &
         (theta >= pi-theta_c .and. theta <= pi+theta_c) ) then
       ! path will cross east/west wall
       len_path = lenx / abs(cos(theta))
    else
       ! path will cross north/south wall
       len_path = leny / abs(sin(theta))
    endif

    ! Choose the number of time steps such that the ice will travel
    ! less than one grid cell per time step

    ntstep = nint(len_path/dx) + 10   ! 10 is arbitrary

    ! Choose the time step such that the ice will go around the domain
    ! exactly once in the chosen number of time steps.

    dt = len_path / (umag*ntstep)

    ! CFL check, just to be sure

    adv_cfl = max (dt*umag*cos(theta)/dx, dt*umag*sin(theta)/dy)
    
    if (adv_cfl >= 1.d0) then
       print*, 'dt is too big for advective CFL; increase ntstep to', ntstep * adv_cfl
       stop
    endif

    ! Print some diagnostics

    if (main_task) then
       print*, ' '
       print*, 'In glissade_test_transport'
       print*, 'nx, ny, nz =', nx, ny, nz
       print*, 'len_path =', len_path
       print*, 'umag (m/yr) =', umag
       print*, 'dt (yr) =', dt
       print*, 'ntstep =', ntstep
       print*, 'theta (deg) =', theta * 180.d0/pi
    endif

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Initial thck'
       do j = ny, 1, -1
          write(6,'(19f7.2)') model%geometry%thck(1:19,j) * thk0
       enddo
    endif

    ! Set uniform ice speed everywhere

    do j = 1, ny-1
    do i = 1, nx-1
       do k = 1, nz
          uvel(k,i,j) = umag * cos(theta)
          vvel(k,i,j) = umag * sin(theta)
       enddo
    enddo
    enddo

    ! Determine which transport scheme

    if (model%options%whichevol == EVOL_UPWIND) then
       do_upwind_transport = .true.
    else
       do_upwind_transport = .false.
    endif

    ! Transport the ice around the domain

    do n = 1, ntstep

       ! call transport scheme
       ! Note: glissade_transport_driver expects dt in seconds, uvel/vvel in m/s

       call glissade_transport_setup_tracers(model)

       call glissade_transport_driver(dt*scyr,                                              &
                                      dx,                        dy,                        &
                                      nx,                        ny,                        &
                                      nz-1,                      model%numerics%sigma,      &
                                      uvel(:,:,:)/scyr,          vvel(:,:,:)/scyr,          &
                                      model%geometry%thck(:,:),                             &
                                      model%climate%acab(:,:),                              &
                                      model%temper%bmlt_ground(:,:),                        &
                                      model%geometry%ntracers,                              &
                                      model%geometry%tracers(:,:,:,:),                      &
                                      model%geometry%tracers_usrf(:,:,:),                   &
                                      model%geometry%tracers_lsrf(:,:,:),                   &
                                      model%options%which_ho_vertical_remap,                &
                                      upwind_transport_in = do_upwind_transport)

       call glissade_transport_finish_tracers(model)

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'New thck, n =', n
          do j = ny, 1, -1
             write(6,'(19f7.2)') model%geometry%thck(1:19,j) * thk0
          enddo
       endif

    enddo  ! ntstep

    if (main_task) print*, 'Done in glissade_test_parallel'

    deallocate(uvel)
    deallocate(vvel)

  end subroutine glissade_test_transport

!=======================================================================

  end module glissade_test

!=======================================================================
