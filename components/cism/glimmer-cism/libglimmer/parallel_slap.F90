!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   parallel_slap.F90 - part of the Community Ice Sheet Model (CISM)  
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

module parallel

  use netcdf
  use glimmer_global, only : dp, sp
  implicit none

  !NOTE: The glam/glissade dycore currently requires nhalo = 2,
  !       whereas the glide dycore requires nhalo = 0.
  !      For glide simulations, we set nhalo = 0 by calling distributed_grid 
  !       with optional argument nhalo = 0.

  integer, save :: nhalo = 2

  !TODO - Define lhalo and uhalo in terms of nhalo.

  integer, save :: lhalo = 2
  integer, save :: uhalo = 2

  integer, save :: staggered_lhalo = 2
  integer, save :: staggered_uhalo = 1

#ifdef _USE_MPI_WITH_SLAP
  logical,save :: main_task
  integer,save :: this_rank
  integer,save :: tasks
  integer,save :: comm
#else
  logical,parameter :: main_task = .true.
  integer,parameter :: this_rank = 0
  integer,parameter :: tasks = 1
#endif

  ! distributed grid
  integer,save :: global_ewn,global_nsn,local_ewn,local_nsn,own_ewn,own_nsn
  integer,save :: global_col_offset, global_row_offset

  integer,save :: ewlb,ewub,nslb,nsub
  integer,save :: east,north,south,west
  integer,save :: ewtasks,nstasks

  !Note: Currently, periodic_bc is always T, since periodic BCs are always applied.                                                             
  !      Outflow BCs require some additional operations on scalars after periodic BCs are applied.
  !      No-penetration BCs are enforced by setting velocity masks, without altering halo updates.

  ! global boundary conditions
  logical,save :: periodic_bc  ! doubly periodic
  logical,save :: outflow_bc   ! if true, set scalars in global halo to zero
                               ! does not apply to staggered variables (e.g., uvel, vvel)

  ! global IDs
  integer,parameter :: ProcsEW = 1

  !TODO - Remove these gathered_* declarations.  No longer used.

  ! JEFF Declarations for undistributed variables on main_task.
  ! Later move to separate module?  These are only temporary until code is completely distributed.
  real(dp),dimension(:,:,:),allocatable :: gathered_efvs  ! Output var from glam_velo_fordsiapstr(), used often
  real(dp),dimension(:,:,:),allocatable :: gathered_efvs2  ! Variable for testing that scatter/gather are inverses
  real(dp),dimension(:,:,:),allocatable :: gathered_uvel  ! Output var from glam_velo_fordsiapstr(), used often
  real(dp),dimension(:,:,:),allocatable :: gathered_vvel  ! Output var from glam_velo_fordsiapstr(), used often
  real(dp),dimension(:,:),allocatable :: gathered_uflx    ! Output var from glam_velo_fordsiapstr(), used often
  real(dp),dimension(:,:),allocatable :: gathered_vflx    ! Output var from glam_velo_fordsiapstr(), used often
  real(dp),dimension(:,:,:),allocatable :: gathered_velnorm  ! Variable calculated in run_ho_diagnostic(), is this used?
  real(dp),dimension(:,:),allocatable :: gathered_thck    ! Used in horizontal_remap_in()
  real(dp),dimension(:,:),allocatable :: gathered_stagthck ! Used in horizontal_remap_in()
  real(sp),dimension(:,:),allocatable :: gathered_acab    ! Used in horizontal_remap_in()
  real(dp),dimension(:,:,:),allocatable :: gathered_temp  ! Used in horizontal_remap_in()
  real(dp),dimension(:,:),allocatable :: gathered_dusrfdew  ! Used in glide_stress()
  real(dp),dimension(:,:),allocatable :: gathered_dusrfdns  ! Used in glide_stress()
  real(dp),dimension(:,:),allocatable :: gathered_dthckdew  ! Used in glide_stress()
  real(dp),dimension(:,:),allocatable :: gathered_dthckdns  ! Used in glide_stress()
  real(dp),dimension(:,:,:),allocatable :: gathered_tauxx   ! Calculated in glide_stress()
  real(dp),dimension(:,:,:),allocatable :: gathered_tauyy   ! Calculated in glide_stress()
  real(dp),dimension(:,:,:),allocatable :: gathered_tauxy   ! Calculated in glide_stress()
  real(dp),dimension(:,:,:),allocatable :: gathered_tauscalar   ! Calculated in glide_stress()
  real(dp),dimension(:,:,:),allocatable :: gathered_tauxz   ! Calculated in glide_stress()
  real(dp),dimension(:,:,:),allocatable :: gathered_tauyz   ! Calculated in glide_stress()
  real(dp),dimension(:,:),allocatable :: gathered_topg  ! Bedrock topology, Used in glide_set_mask()
  integer,dimension(:,:),allocatable :: gathered_thkmask  ! Calculated in glide_set_mask()
  real(dp),dimension(:,:),allocatable :: gathered_marine_bc_normal  ! Calculated in glide_marine_margin_normal()
  real(dp),dimension(:,:,:),allocatable :: gathered_surfvel   ! Used in calc_gline_flux()
  real(dp),dimension(:,:),allocatable :: gathered_gline_flux   ! Calculated in calc_gline_flux()
  real(dp),dimension(:,:),allocatable :: gathered_ubas   ! Used in calc_gline_flux()
  real(dp),dimension(:,:),allocatable :: gathered_vbas   ! Used in calc_gline_flux()
  real(dp),dimension(:,:),allocatable :: gathered_relx   ! Used in glide_marinlim()
  real(dp),dimension(:,:,:),allocatable :: gathered_flwa   ! Used in glide_marinlim()
  real(sp),dimension(:,:),allocatable :: gathered_calving   ! Used in glide_marinlim()
  real(sp),dimension(:,:),allocatable :: gathered_backstress   ! Used in glide_marinlim()
  real(dp),dimension(:,:),allocatable :: gathered_usrf   ! Used in glide_marinlim()
  logical,dimension(:,:),allocatable :: gathered_backstressmap ! Used in glide_marinlim()
  real(dp),dimension(:,:),allocatable :: gathered_tau_x   ! Calculated in calc_basal_shear()
  real(dp),dimension(:,:),allocatable :: gathered_tau_y   ! Calculated in calc_basal_shear()
  real(dp),dimension(:,:),allocatable :: gathered_lsrf   ! Used in glide_marinlim()

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
     module procedure distributed_get_var_real8_1d
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
     module procedure distributed_put_var_real8_1d
     module procedure distributed_put_var_real8_2d
     module procedure distributed_put_var_real8_3d

     !TODO - Should the parallel_put_var routines be part of this interface?
     module procedure parallel_put_var_real4
     module procedure parallel_put_var_real8
  end interface

  interface parallel_convert_haloed_to_nonhaloed
     module procedure parallel_convert_haloed_to_nonhaloed_real4_2d
     module procedure parallel_convert_haloed_to_nonhaloed_real8_2d
  end interface parallel_convert_haloed_to_nonhaloed

  interface parallel_convert_nonhaloed_to_haloed
     module procedure parallel_convert_nonhaloed_to_haloed_real4_2d
     module procedure parallel_convert_nonhaloed_to_haloed_real8_2d
  end interface parallel_convert_nonhaloed_to_haloed
  
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

  interface distributed_scatter_var
     module procedure distributed_scatter_var_integer_2d
     module procedure distributed_scatter_var_logical_2d
     module procedure distributed_scatter_var_real4_2d
     module procedure distributed_scatter_var_real4_3d
     module procedure distributed_scatter_var_real8_2d
     module procedure distributed_scatter_var_real8_3d
  end interface

  interface global_sum
     module procedure global_sum_real8_scalar
     module procedure global_sum_real8_1d
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

  interface parallel_print
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

  interface parallel_reduce_sum
     module procedure parallel_reduce_sum_integer
     module procedure parallel_reduce_sum_real4
     module procedure parallel_reduce_sum_real8
     module procedure parallel_reduce_sum_real8_nvar
  end interface

  interface parallel_reduce_max
     module procedure parallel_reduce_max_integer
     module procedure parallel_reduce_max_real4
     module procedure parallel_reduce_max_real8
  end interface

  interface parallel_reduce_min
     module procedure parallel_reduce_min_integer
     module procedure parallel_reduce_min_real4
     module procedure parallel_reduce_min_real8
  end interface

  ! This reduce interface determines the global min value and the processor on which it occurs
  interface parallel_reduce_maxloc
     module procedure parallel_reduce_maxloc_integer
     module procedure parallel_reduce_maxloc_real4
     module procedure parallel_reduce_maxloc_real8
  end interface

  ! This reduce interface determines the global min value and the processor on which it occurs
  interface parallel_reduce_minloc
     module procedure parallel_reduce_minloc_integer
     module procedure parallel_reduce_minloc_real4
     module procedure parallel_reduce_minloc_real8
  end interface

  interface staggered_parallel_halo
     module procedure staggered_parallel_halo_integer_2d
     module procedure staggered_parallel_halo_integer_3d
     module procedure staggered_parallel_halo_real8_2d
     module procedure staggered_parallel_halo_real8_3d
     module procedure staggered_parallel_halo_real8_4d
  end interface

  interface staggered_parallel_halo_extrapolate
     module procedure staggered_parallel_halo_extrapolate_integer_2d
     module procedure staggered_parallel_halo_extrapolate_real8_2d
  end interface

contains

  subroutine broadcast_character(c, proc)
    implicit none
    character(len=*) :: c
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_character

  subroutine broadcast_integer(i, proc)
    implicit none
    integer :: i
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_integer

  subroutine broadcast_integer_1d(a, proc)
    implicit none
    integer,dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_integer_1d

  subroutine broadcast_logical(l, proc)
    implicit none
    logical :: l
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_logical

  subroutine broadcast_real4(r, proc)
    implicit none
    real(sp) :: r
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_real4

  subroutine broadcast_real4_1d(a, proc)
    real(sp),dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_real4_1d

  subroutine broadcast_real8(r, proc)
    implicit none
    real(dp) :: r
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_real8

  subroutine broadcast_real8_1d(a, proc)
    implicit none
    real(dp),dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from - not relevant to serial version
  end subroutine broadcast_real8_1d

  function distributed_get_var_integer_2d(ncid,varid,values,start)

    implicit none
    integer :: distributed_get_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    integer :: ilo, ihi, jlo, jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + staggered_lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + staggered_lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_get_var_integer_2d =  &
          nf90_get_var(ncid,varid,values(ilo:ihi,jlo:jhi),start)

    endif

  end function distributed_get_var_integer_2d

  function distributed_get_var_real4_1d(ncid,varid,values,start)

    implicit none
    integer :: distributed_get_var_real4_1d,ncid,varid
    integer,dimension(:) :: start
    real(sp),dimension(:) :: values

    integer :: status, x1id, y1id
    integer :: ilo, ihi

    ! begin

    if (main_task) then

       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
       if (varid==x1id) then
          ilo = 1+lhalo
          ihi = local_ewn - uhalo
       else if (varid==y1id) then
          ilo = 1+lhalo
          ihi = local_nsn - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_get_var_real4_1d = &
              nf90_get_var(ncid,varid,values(ilo:ihi),start)

    endif

  end function distributed_get_var_real4_1d

  function distributed_get_var_real4_2d(ncid,varid,values,start)

    implicit none
    integer :: distributed_get_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(sp),dimension(:,:) :: values

    integer :: ilo, ihi, jlo, jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_get_var_real4_2d =  &
          nf90_get_var(ncid,varid,values(ilo:ihi,jlo:jhi),start)

    endif

  end function distributed_get_var_real4_2d

  !WHL - added this function

  function distributed_get_var_real8_1d(ncid,varid,values,start)

    implicit none
    integer :: distributed_get_var_real8_1d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:) :: values

    integer :: status, x1id, y1id
    integer :: ilo, ihi

    ! begin

    if (main_task) then

       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
       if (varid==x1id) then
          ilo = 1+lhalo
          ihi = local_ewn - uhalo
       else if (varid==y1id) then
          ilo = 1+lhalo
          ihi = local_nsn - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_get_var_real8_1d = &
              nf90_get_var(ncid,varid,values(ilo:ihi),start)

    endif

  end function distributed_get_var_real8_1d

  function distributed_get_var_real8_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_get_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:) :: values

    integer :: ilo, ihi, jlo, jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_get_var_real8_2d =  &
          nf90_get_var(ncid,varid,values(ilo:ihi,jlo:jhi),start)

    endif

  end function distributed_get_var_real8_2d

  function distributed_get_var_real8_3d(ncid,varid,values,start)

    implicit none
    integer :: distributed_get_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:,:) :: values

    integer :: ilo, ihi, jlo, jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_get_var_real8_3d =  &
            nf90_get_var(ncid,varid,values(ilo:ihi,jlo:jhi,:),start)

    endif

  end function distributed_get_var_real8_3d


  subroutine distributed_grid(ewn,      nsn,  &
                              nhalo_in, outflow_bc_in)

    implicit none

    integer, intent(inout) :: ewn, nsn          ! global grid dimensions
    integer, intent(in), optional :: nhalo_in   ! number of rows of halo cells
    logical, intent(in), optional :: outflow_bc_in  ! true for outflow global BCs
                                                    ! (scalars in global halo set to zero)
           
    integer :: ewrank,nsrank

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
          endif
       endif
       nhalo = nhalo_in
       lhalo = nhalo
       uhalo = nhalo
       staggered_lhalo = lhalo
       staggered_uhalo = max(uhalo-1, 0)
    endif

    ! initialize some grid quantities to be consistent with parallel_mpi

    global_ewn = ewn
    global_nsn = nsn

    global_row_offset = 0
    global_col_offset = 0
  
    ewrank = 0
    nsrank = 0
    ewtasks = 1
    nstasks = 1

    east = 0      ! all halo updates are local copies by the main task
    west = 0
    north = 0
    south = 0

! Trey's original code
!    ewlb = 1
!    ewub = global_ewn
!    local_ewn = ewub-ewlb+1
!    own_ewn = local_ewn-lhalo-uhalo
!    ewn = local_ewn

!    nslb = 1
!    nsub = global_nsn
!    local_nsn = nsub-nslb+1
!    own_nsn = local_nsn-lhalo-uhalo
!    nsn = local_nsn

!WHL - modified code for nonzero halo values
    ewlb = 1 - lhalo
    ewub = global_ewn + uhalo
    local_ewn = ewub - ewlb + 1
    own_ewn = local_ewn - lhalo - uhalo
    ewn = local_ewn

    nslb = 1 - lhalo
    nsub = global_nsn + uhalo
    local_nsn = nsub - nslb + 1
    own_nsn = local_nsn - lhalo - uhalo
    nsn = local_nsn

    ! set periodic BC as the default
    ! Note: Currently, it is not strictly necessary to set periodic_bc = T, because
    !       all halo updates carried out for periodic BC are also carried out for other BCs.
    !       (Outflow BCs simply add some additional operations on scalars.)                                                                                        !       But setting it here as the default in case that changes. 
    periodic_bc = .true.

    if (present(outflow_bc_in)) then
       outflow_bc = outflow_bc_in
    else
       outflow_bc = .false.
    endif
    
    ! Print grid geometry
    write(*,*) "Process ", this_rank, " Total = ", tasks, " ewtasks = ", ewtasks, " nstasks = ", nstasks
    write(*,*) "Process ", this_rank, " ewrank = ", ewrank, " nsrank = ", nsrank
    write(*,*) "Process ", this_rank, " l_ewn = ", local_ewn, " o_ewn = ", own_ewn
    write(*,*) "Process ", this_rank, " l_nsn = ", local_nsn, " o_nsn = ", own_nsn
    write(*,*) "Process ", this_rank, " ewlb = ", ewlb, " ewub = ", ewub
    write(*,*) "Process ", this_rank, " nslb = ", nslb, " nsub = ", nsub
    write(*,*) "Process ", this_rank, " east = ", east, " west = ", west
    write(*,*) "Process ", this_rank, " north = ", north, " south = ", south
    write(*,*) "Process ", this_rank, " ew_vars = ", own_ewn, " ns_vars = ", own_nsn

  end subroutine distributed_grid

  function distributed_execution()
     ! Returns if running distributed or not.
     logical distributed_execution

     distributed_execution = .false.
  end function distributed_execution

  subroutine distributed_gather_var_integer_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    integer,dimension(:,:),intent(in) :: values
    integer,dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

      !WHL - Commented code will not work if the local arrays include halo cells
!!    allocate(global_values(size(values,1), size(values,2)))
!!    global_values(:,:) = values(:,:)
    allocate(global_values(size(values,1)-uhalo-lhalo, size(values,2)-uhalo-lhalo))
    global_values(:,:) = values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine distributed_gather_var_integer_2d

  subroutine distributed_gather_var_logical_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    logical,dimension(:,:),intent(in) :: values
    logical,dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

      !WHL - Commented code will not work if the local arrays include halo cells
!!    allocate(global_values(size(values,1), size(values,2)))
!!    global_values(:,:) = values(:,:)
    allocate(global_values(size(values,1)-uhalo-lhalo, size(values,2)-uhalo-lhalo))
    global_values(:,:) = values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine distributed_gather_var_logical_2d

  subroutine distributed_gather_var_real4_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    real(sp),dimension(:,:),intent(in) :: values
    real(sp),dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

      !WHL - Commented code will not work if the local arrays include halo cells
!!    allocate(global_values(size(values,1), size(values,2)))
!!    global_values(:,:) = values(:,:)
    allocate(global_values(size(values,1)-uhalo-lhalo, size(values,2)-uhalo-lhalo))
    global_values(:,:) = values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine distributed_gather_var_real4_2d

  subroutine distributed_gather_var_real4_3d(values, global_values, ld1, ud1)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    real(sp),dimension(:,:,:),intent(in) :: values
    real(sp),dimension(:,:,:),allocatable,intent(inout) :: global_values
    integer,optional,intent(in) :: ld1, ud1

    integer :: d1l,d1u

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
       d1u = size(values,1)
    endif
    if (size(values,1) /= d1u-d1l+1) then
       write(*,*) "size(values,1) .ne. d1u-d1l+1 in gather call"
       call parallel_stop(__FILE__, __LINE__)
    endif

    !WHL - Commented code will not work if the local arrays include halo cells
!!    allocate(global_values(d1l:d1u, size(values,2), size(values,3)))
!!    global_values(d1l:d1u,:,:) = values(1:size(values,1),:,:)
    allocate(global_values(d1l:d1u, size(values,2)-uhalo-lhalo, size(values,3)-uhalo-lhalo))
    global_values(d1l:d1u,:,:) = values(1:size(values,1), 1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine distributed_gather_var_real4_3d

  subroutine distributed_gather_var_real8_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    real(dp),dimension(:,:),intent(in) :: values
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

      !WHL - Commented code will not work if the local arrays include halo cells
!!    allocate(global_values(size(values,1), size(values,2)))
!!    global_values(:,:) = values(:,:)
    allocate(global_values(size(values,1)-uhalo-lhalo, size(values,2)-uhalo-lhalo))
    global_values(:,:) = values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine distributed_gather_var_real8_2d

  subroutine distributed_gather_var_real8_3d(values, global_values, ld1, ud1)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    real(dp),dimension(:,:,:),intent(in) :: values
    real(dp),dimension(:,:,:),allocatable,intent(inout) :: global_values
    integer,optional,intent(in) :: ld1, ud1

    integer :: d1l,d1u

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
       d1u = size(values,1)
    endif
    if (size(values,1) /= d1u-d1l+1) then
       write(*,*) "size(values,1) .ne. d1u-d1l+1 in gather call"
       call parallel_stop(__FILE__, __LINE__)
    endif

    !WHL - Commented code will not work if the local arrays include halo cells
!!    allocate(global_values(d1l:d1u, size(values,2), size(values,3)))
!!    global_values(d1l:d1u,:,:) = values(1:size(values,1),:,:)
    allocate(global_values(d1l:d1u, size(values,2)-uhalo-lhalo, size(values,3)-uhalo-lhalo))
    global_values(d1l:d1u,:,:) = values(1:size(values,1), 1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine distributed_gather_var_real8_3d

  function distributed_isparallel()
     implicit none
     logical :: distributed_isparallel

     distributed_isparallel = .false.
  end function distributed_isparallel

  function distributed_owner(ew,ewn,ns,nsn)
    implicit none
    logical :: distributed_owner
    integer :: ew,ewn,ns,nsn
    ! begin
    distributed_owner = .true.
  end function distributed_owner

  subroutine distributed_print_integer_2d(name,values)
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j,k

    write(ts,'(i3.3)') tasks
    open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
    if (size(values,1)<local_ewn) then
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    else
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    end if
    close(u)
  end subroutine distributed_print_integer_2d

  subroutine distributed_print_real8_2d(name,values)
    implicit none
    character(*) :: name
    real(dp),dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j,k

    write(ts,'(i3.3)') tasks
    open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
    if (size(values,1)<local_ewn) then
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    else
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    end if
    close(u)
  end subroutine distributed_print_real8_2d

  subroutine distributed_print_real8_3d(name,values)
    implicit none
    character(*) :: name
    real(dp),dimension(:,:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j,k

    write(ts,'(i3.3)') tasks
    open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
    if (size(values,2)<local_ewn) then
       do j = lbound(values,3),ubound(values,3)
          do i = lbound(values,2),ubound(values,2)
             write(u,'(2i6,100g15.5e3)') j,i,values(:,i,j)
          end do
          write(u,'()')
       end do
    else
       do j = lbound(values,3),ubound(values,3)
          do i = lbound(values,2),ubound(values,2)
             write(u,'(2i6,100g15.5e3)') j,i,values(:,i,j)
          end do
          write(u,'()')
       end do
    end if
    close(u)
  end subroutine distributed_print_real8_3d


  function distributed_put_var_integer_2d(ncid,varid,values,start)

    implicit none
    integer :: distributed_put_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    integer :: ilo,ihi,jlo,jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_put_var_integer_2d = &
          nf90_put_var(ncid,varid,values(ilo:ihi,jlo:jhi),start)

     endif

  end function distributed_put_var_integer_2d

  function distributed_put_var_real4_1d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_real4_1d,ncid,varid
    integer,dimension(:),optional :: start
    real(sp),dimension(:) :: values

    integer :: status, x0id, x1id, y0id, y1id
    integer :: ilo, ihi

    ! begin

    if (main_task) then

       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"y1",y1id)

       if (varid==x0id) then         ! staggered grid
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
       else if (varid==x1id) then    ! unstaggered grid
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
       else if (varid==y0id) then    ! staggered grid
          ilo = 1 + lhalo
          ihi = local_nsn - 1 - uhalo
       else if (varid==y1id) then    ! unstaggered grid
          ilo = 1 + lhalo
          ihi = local_nsn - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       if (present(start)) then
         distributed_put_var_real4_1d = nf90_put_var(ncid,varid,values(ilo:ihi),start)
       else
         distributed_put_var_real4_1d = nf90_put_var(ncid,varid,values(ilo:ihi))
       endif

    endif

  end function distributed_put_var_real4_1d


  function distributed_put_var_real4_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(sp),dimension(:,:) :: values

    integer :: ilo,ihi,jlo,jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_put_var_real4_2d = &
          nf90_put_var(ncid,varid,values(ilo:ihi,jlo:jhi),start)

    endif

  end function distributed_put_var_real4_2d

  !WHL - added this function
  function distributed_put_var_real8_1d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_real8_1d,ncid,varid
    integer,dimension(:),optional :: start
    real(dp),dimension(:) :: values

    integer :: status, x0id, x1id, y0id, y1id
    integer :: ilo, ihi

    ! begin

    if (main_task) then

       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"y1",y1id)

       if (varid==x0id) then         ! staggered grid
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
       else if (varid==x1id) then    ! unstaggered grid
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
       else if (varid==y0id) then    ! staggered grid
          ilo = 1 + lhalo
          ihi = local_nsn - 1 - uhalo
       else if (varid==y1id) then    ! unstaggered grid
          ilo = 1 + lhalo
          ihi = local_nsn - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       if (present(start)) then
         distributed_put_var_real8_1d = nf90_put_var(ncid,varid,values(ilo:ihi),start)
       else
         distributed_put_var_real8_1d = nf90_put_var(ncid,varid,values(ilo:ihi))
       endif

    endif

  end function distributed_put_var_real8_1d


  function distributed_put_var_real8_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:) :: values

    integer :: ilo,ihi,jlo,jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_put_var_real8_2d = &
          nf90_put_var(ncid,varid,values(ilo:ihi,jlo:jhi),start)

    endif

  end function distributed_put_var_real8_2d

  !TODO - Should we assume that ewn is the first index?
  !       Typically it is the 2nd index.
  function distributed_put_var_real8_3d(ncid,varid,values,start)

    implicit none
    integer :: distributed_put_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:,:) :: values

    integer :: ilo,ihi,jlo,jhi

    ! begin

    if (main_task) then

       if (size(values,1)==local_ewn) then
          ilo = 1 + lhalo
          ihi = local_ewn - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - uhalo
       else if (size(values,1)==local_ewn-1) then
          ilo = 1 + lhalo
          ihi = local_ewn - 1 - uhalo
          jlo = 1 + lhalo
          jhi = local_nsn - 1 - uhalo
       else
          call parallel_stop(__FILE__,__LINE__)
       end if

       distributed_put_var_real8_3d =   &
          nf90_put_var(ncid,varid,values(ilo:ihi,jlo:jhi,:),start)

    endif

  end function distributed_put_var_real8_3d

  subroutine distributed_scatter_var_integer_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    integer,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    integer,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = global_values(:,:)

    deallocate(global_values)  ! automatic deallocation
  end subroutine distributed_scatter_var_integer_2d

  subroutine distributed_scatter_var_logical_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    logical,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    logical,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = global_values(:,:)

    deallocate(global_values)  ! automatic deallocation
  end subroutine distributed_scatter_var_logical_2d

  subroutine distributed_scatter_var_real4_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    implicit none
    real(sp),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(sp),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = global_values(:,:)

    deallocate(global_values)  ! automatic deallocation
  end subroutine distributed_scatter_var_real4_2d

  subroutine distributed_scatter_var_real4_3d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! NOTE: The horizontal dimensions are assumed to be the second and third dimensions.
    !       Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    real(sp),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(sp),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(:, 1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = global_values(:,:,:)

    deallocate(global_values)  ! automatic deallocation
  end subroutine distributed_scatter_var_real4_3d

  subroutine distributed_scatter_var_real8_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    real(dp),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = global_values(:,:)

    deallocate(global_values)  ! automatic deallocation
  end subroutine distributed_scatter_var_real8_2d

  subroutine distributed_scatter_var_real8_3d(values, global_values, deallocflag)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! NOTE: The horizontal dimensions are assumed to be the second and third dimensions.
    !       Variables are assumed to lie on the scalar grid (at cell centers).
    implicit none
    real(dp),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(dp),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    logical,optional :: deallocflag
    logical :: deallocmem

    !TODO - Why does this distributed_scatter subroutine have a deallocmem flag, but not the others?
    if (present(deallocflag)) then
       deallocmem = deallocflag
    else
       deallocmem = .true.
    endif

    ! begin
    values(:, 1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = global_values(:,:,:)

    if (deallocmem) deallocate(global_values) ! automatic deallocation
  end subroutine distributed_scatter_var_real8_3d

  subroutine global_sum_real8_scalar(x)
    implicit none
    real(dp) :: x
  end subroutine global_sum_real8_scalar

  subroutine global_sum_real8_1d(x)
    implicit none
    real(dp),dimension(:) :: x
  end subroutine global_sum_real8_1d

  subroutine not_parallel(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    write(0,*) "WARNING: not parallel in ",file," at line ",line
  end subroutine not_parallel

  subroutine parallel_barrier
    implicit none
  end subroutine parallel_barrier

  function parallel_boundary(ew,ewn,ns,nsn)
    implicit none
    logical :: parallel_boundary
    integer :: ew,ewn,ns,nsn
    ! begin
    parallel_boundary = (ew==1.or.ew==ewn.or.ns==1.or.ns==nsn)
  end function parallel_boundary

  function parallel_close(ncid)
    implicit none
    integer :: ncid,parallel_close
    ! begin
    if (main_task) parallel_close = nf90_close(ncid)
    call broadcast(parallel_close)
  end function parallel_close

  subroutine parallel_convert_haloed_to_nonhaloed_real4_2d(input_with_halo, output_no_halo)
    ! Given an input array that has halo cells, return an output array without halo cells
    real(sp),dimension(:,:), intent(in)  :: input_with_halo
    real(sp),dimension(:,:), intent(out) :: output_no_halo

    if (size(input_with_halo,1) /= local_ewn .or. size(input_with_halo,2) /= local_nsn) then
       write(*,*) "Unexpected size for input_with_halo: ", &
            size(input_with_halo,1), size(input_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_no_halo,1) /= own_ewn .or. size(output_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for output_no_halo: ", &
            size(output_no_halo,1), size(output_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_no_halo(1:own_ewn, 1:own_nsn) = &
         input_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine parallel_convert_haloed_to_nonhaloed_real4_2d

  subroutine parallel_convert_haloed_to_nonhaloed_real8_2d(input_with_halo, output_no_halo)
    ! Given an input array that has halo cells, return an output array without halo cells
    real(dp),dimension(:,:), intent(in)  :: input_with_halo
    real(dp),dimension(:,:), intent(out) :: output_no_halo

    if (size(input_with_halo,1) /= local_ewn .or. size(input_with_halo,2) /= local_nsn) then
       write(*,*) "Unexpected size for input_with_halo: ", &
            size(input_with_halo,1), size(input_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_no_halo,1) /= own_ewn .or. size(output_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for output_no_halo: ", &
            size(output_no_halo,1), size(output_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_no_halo(1:own_ewn, 1:own_nsn) = &
         input_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

  end subroutine parallel_convert_haloed_to_nonhaloed_real8_2d

  subroutine parallel_convert_nonhaloed_to_haloed_real4_2d(input_no_halo, output_with_halo)
    ! Given an input array without halo cells, return an output array with halo cells
    real(sp),dimension(:,:), intent(in)  :: input_no_halo
    real(sp),dimension(:,:), intent(out) :: output_with_halo
    
    if (size(input_no_halo,1) /= own_ewn .or. size(input_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for input_no_halo: ", &
            size(input_no_halo,1), size(input_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_with_halo,1) /= local_ewn .or. size(output_with_halo,2) /= local_nsn) then
       write(*,*) "Unexpected size for output_with_halo: ", &
            size(output_with_halo,1), size(output_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = &
         input_no_halo(1:own_ewn, 1:own_nsn)

    call parallel_halo(output_with_halo)
    
  end subroutine parallel_convert_nonhaloed_to_haloed_real4_2d

  subroutine parallel_convert_nonhaloed_to_haloed_real8_2d(input_no_halo, output_with_halo)
    ! Given an input array without halo cells, return an output array with halo cells
    real(dp),dimension(:,:), intent(in)  :: input_no_halo
    real(dp),dimension(:,:), intent(out) :: output_with_halo
    
    if (size(input_no_halo,1) /= own_ewn .or. size(input_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for input_no_halo: ", &
            size(input_no_halo,1), size(input_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_with_halo,1) /= local_ewn .or. size(output_with_halo,2) /= local_nsn) then
       write(*,*) "Unexpected size for output_with_halo: ", &
            size(output_with_halo,1), size(output_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = &
         input_no_halo(1:own_ewn, 1:own_nsn)

    call parallel_halo(output_with_halo)
    
  end subroutine parallel_convert_nonhaloed_to_haloed_real8_2d

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

#ifdef _USE_MPI_WITH_SLAP
  subroutine parallel_finalise
    use mpi_mod
    implicit none
    integer :: ierror 
    ! begin 
    call mpi_finalize(ierror)
  end subroutine
#else
  subroutine parallel_finalise
    implicit none
  end subroutine parallel_finalise
#endif

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
    real(sp) :: values
    ! begin
    if (main_task) parallel_get_att_real4 = &
         nf90_get_att(ncid,varid,name,values)
  end function parallel_get_att_real4

  function parallel_get_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real4_1d,varid
    character(len=*) :: name
    real(sp),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real4_1d = &
         nf90_get_att(ncid,varid,name,values)
  end function parallel_get_att_real4_1d

  function parallel_get_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8,varid
    character(len=*) :: name
    real(dp) :: values
    ! begin
    if (main_task) parallel_get_att_real8 = &
         nf90_get_att(ncid,varid,name,values)
  end function parallel_get_att_real8

  function parallel_get_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8_1d,varid
    character(len=*) :: name
    real(dp),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real8_1d = &
         nf90_get_att(ncid,varid,name,values)
  end function parallel_get_att_real8_1d

  function parallel_get_var_integer_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_integer_1d,varid
    integer,dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_integer_1d = &
         nf90_get_var(ncid,varid,values)
  end function parallel_get_var_integer_1d

  function parallel_get_var_real4_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_real4_1d,varid
    real(sp),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_real4_1d = &
         nf90_get_var(ncid,varid,values)
  end function parallel_get_var_real4_1d

  function parallel_get_var_real8_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_real8_1d,varid
    real(dp),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_real8_1d = &
         nf90_get_var(ncid,varid,values)
  end function parallel_get_var_real8_1d

  !TODO - Is function parallel_globalID still needed?

  function parallel_globalID(locns, locew, upstride)
    ! Returns a unique ID for a given row and column reference that is identical across all processors.
    ! For instance if Proc 2: (17,16) is the same global cell as Proc 3: (17,1), then the globalID will be the same for both.
    ! These IDs are spaced upstride apart.  upstride = number of vertical layers.  Typically (upn) + number of ghost layers (2 = top and bottom)
    integer,intent(in) :: locns, locew, upstride
    integer :: parallel_globalID
    ! locns is local NS (row) grid index
    ! locew is local EW (col) grid index
    integer :: global_row, global_col, global_ID

    global_row = (locns - uhalo) + this_rank/ProcsEW * own_nsn
    	! Integer division required for this_rank/ProcsEW
    global_col = (locew - lhalo) + mod(this_rank, ProcsEW) * own_ewn
        ! There are ProcsEW processors per row.

    global_ID = ((global_row - 1) * global_ewn + (global_col - 1)) * upstride + 1

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

    ! including global domain halo adds lhalo to offsets
    global_row = locns - lhalo
    global_col = locew - lhalo

    ! including global domain halo adds (lhalo + uhalo) to global_ewn
    global_ID = ((global_row - 1)*(global_ewn) + (global_col - 1)) * upstride + 1

    !return value
    parallel_globalID_scalar = global_ID

  end function parallel_globalID_scalar


  subroutine parallel_globalindex(ilocal, jlocal, iglobal, jglobal)
    ! Calculates the global i,j indices from the local i,j indices
    integer,intent(in)  :: ilocal,  jlocal  ! These include the halos
    integer,intent(out) :: iglobal, jglobal ! These do NOT include halos

   ! No check is currently made for being located in the global (periodic) halo
    iglobal = (ilocal - lhalo)
    jglobal = (jlocal - lhalo)
  end subroutine parallel_globalindex

  subroutine parallel_localindex(iglobal, jglobal, ilocal, jlocal, rlocal)
    ! Calculates the local i,j indices and rank from the global i,j indices                                                                                               
    integer,intent(in) :: iglobal, jglobal 
    integer,intent(out)  :: ilocal, jlocal, rlocal

    ilocal = iglobal + lhalo 
    jlocal = jglobal + lhalo 
    rlocal = this_rank

    ! Check whether these are valid values of ilocal and jlocal
    if ( (ilocal > lhalo .and. ilocal <= lhalo + own_ewn)  &
                         .and.                             &
         (jlocal > lhalo .and. jlocal <= lhalo + own_nsn) ) then
         ! global indices are valid
    else ! global indices are invalid 
       if (main_task) then
          write(*,*) 'Invalid global indices: iglobal, jglobal =', iglobal, jglobal
          call parallel_stop(__FILE__,__LINE__)
       endif
    endif
  end subroutine parallel_localindex


  subroutine parallel_halo_integer_2d(a)

    implicit none
    integer,dimension(:,:) :: a

    integer,dimension(lhalo,local_nsn-lhalo-uhalo) :: ecopy
    integer,dimension(uhalo,local_nsn-lhalo-uhalo) :: wcopy
    integer,dimension(local_ewn,lhalo) :: ncopy
    integer,dimension(local_ewn,uhalo) :: scopy

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    if (outflow_bc) then

       a(:lhalo,1+lhalo:local_nsn-uhalo) = 0
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = 0
       a(:,:lhalo) = 0
       a(:,local_nsn-uhalo+1:) = 0

    else    ! periodic BC

       ecopy(:,:) = a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
       wcopy(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
       a(:lhalo,1+lhalo:local_nsn-uhalo) = ecopy(:,:)
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wcopy(:,:)

       ncopy(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
       scopy(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
       a(:,:lhalo) = ncopy(:,:)
       a(:,local_nsn-uhalo+1:) = scopy(:,:)

    endif

  end subroutine parallel_halo_integer_2d


  subroutine parallel_halo_logical_2d(a)

    implicit none
    logical,dimension(:,:) :: a

    logical,dimension(lhalo,local_nsn-lhalo-uhalo) :: ecopy
    logical,dimension(uhalo,local_nsn-lhalo-uhalo) :: wcopy
    logical,dimension(local_ewn,lhalo) :: ncopy
    logical,dimension(local_ewn,uhalo) :: scopy

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    if (outflow_bc) then

       a(:lhalo,1+lhalo:local_nsn-uhalo) = .false.
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = .false.
       a(:,:lhalo) = .false.
       a(:,local_nsn-uhalo+1:) = .false.

    else    ! periodic BC

       ecopy(:,:) = a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
       wcopy(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
       a(:lhalo,1+lhalo:local_nsn-uhalo) = ecopy(:,:)
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wcopy(:,:)

       ncopy(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
       scopy(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
       a(:,:lhalo) = ncopy(:,:)
       a(:,local_nsn-uhalo+1:) = scopy(:,:)

    endif

  end subroutine parallel_halo_logical_2d


  subroutine parallel_halo_real4_2d(a)

    implicit none
    real(sp),dimension(:,:) :: a

    real(sp),dimension(lhalo,local_nsn-lhalo-uhalo) :: ecopy
    real(sp),dimension(uhalo,local_nsn-lhalo-uhalo) :: wcopy
    real(sp),dimension(local_ewn,lhalo) :: ncopy
    real(sp),dimension(local_ewn,uhalo) :: scopy

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    if (outflow_bc) then

       a(:lhalo,1+lhalo:local_nsn-uhalo) = 0.
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = 0.
       a(:,:lhalo) = 0.
       a(:,local_nsn-uhalo+1:) = 0.

    else    ! periodic BC

       ecopy(:,:) = a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
       wcopy(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
       a(:lhalo,1+lhalo:local_nsn-uhalo) = ecopy(:,:)
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wcopy(:,:)

       ncopy(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
       scopy(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
       a(:,:lhalo) = ncopy(:,:)
       a(:,local_nsn-uhalo+1:) = scopy(:,:)

    endif

  end subroutine parallel_halo_real4_2d


  subroutine parallel_halo_real8_2d(a, periodic_offset_ew, periodic_offset_ns)

    !WHL - added optional arguments for periodic offsets, to support ismip-hom test cases

    implicit none
    real(dp),dimension(:,:) :: a
    real(dp), intent(in), optional :: &
       periodic_offset_ew,  &! offset halo values by this amount
                             ! if positive, the offset is positive for W halo, negative for E halo 
       periodic_offset_ns    ! offset halo values by this amount
                             ! if positive, the offset is positive for S halo, negative for N halo

    real(dp),dimension(lhalo,local_nsn-lhalo-uhalo) :: ecopy
    real(dp),dimension(uhalo,local_nsn-lhalo-uhalo) :: wcopy
    real(dp),dimension(local_ewn,lhalo) :: ncopy
    real(dp),dimension(local_ewn,uhalo) :: scopy

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    if (outflow_bc) then

       a(:lhalo,1+lhalo:local_nsn-uhalo) = 0.d0
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = 0.d0
       a(:,:lhalo) = 0.d0
       a(:,local_nsn-uhalo+1:) = 0.d0

    else    ! periodic BC

       ecopy(:,:) = a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
       wcopy(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
       a(:lhalo,1+lhalo:local_nsn-uhalo) = ecopy(:,:)
       a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wcopy(:,:)

       if (present(periodic_offset_ew)) then
          if (periodic_offset_ew /= 0.d0) then
             a(:lhalo,1+lhalo:local_nsn-uhalo) =   &
                  a(:lhalo,1+lhalo:local_nsn-uhalo) + periodic_offset_ew
             a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) =    &
                  a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) - periodic_offset_ew
          endif
       endif

       ncopy(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
       scopy(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
       a(:,:lhalo) = ncopy(:,:)
       a(:,local_nsn-uhalo+1:) = scopy(:,:)

       if (present(periodic_offset_ns)) then
          if (periodic_offset_ns /= 0.d0) then
             a(:,:lhalo) = a(:,:lhalo) + periodic_offset_ns
             a(:,local_nsn-uhalo+1:) = a(:,local_nsn-uhalo+1:) - periodic_offset_ns
          endif
       endif

    endif  ! open or periodic BC

  end subroutine parallel_halo_real8_2d


  subroutine parallel_halo_real8_3d(a)

    implicit none
    real(dp),dimension(:,:,:) :: a

    real(dp),dimension(size(a,1),lhalo,local_nsn-lhalo-uhalo) :: ecopy
    real(dp),dimension(size(a,1),uhalo,local_nsn-lhalo-uhalo) :: wcopy
    real(dp),dimension(size(a,1),local_ewn,lhalo) :: ncopy
    real(dp),dimension(size(a,1),local_ewn,uhalo) :: scopy

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn .or. size(a,3)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,2), ",", size(a,3), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    if (outflow_bc) then

       a(:,:lhalo,1+lhalo:local_nsn-uhalo) = 0.d0
       a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = 0.d0
       a(:,:,:lhalo) = 0.d0
       a(:,:,local_nsn-uhalo+1:) = 0.d0

    else    ! periodic BC

       ecopy(:,:,:) = a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
       wcopy(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
       a(:,:lhalo,1+lhalo:local_nsn-uhalo) = ecopy(:,:,:)
       a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wcopy(:,:,:)

       ncopy(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
       scopy(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
       a(:,:,:lhalo) = ncopy(:,:,:)
       a(:,:,local_nsn-uhalo+1:) = scopy(:,:,:)

    endif

  end subroutine parallel_halo_real8_3d

  function parallel_halo_verify_integer_2d(a)
    implicit none
    integer,dimension(:,:) :: a
    logical :: parallel_halo_verify_integer_2d
    parallel_halo_verify_integer_2d = .true.
  end function parallel_halo_verify_integer_2d

  function parallel_halo_verify_real8_2d(a)
    implicit none
    real(dp),dimension(:,:) :: a
    logical :: parallel_halo_verify_real8_2d
    parallel_halo_verify_real8_2d = .true.
  end function parallel_halo_verify_real8_2d

  function parallel_halo_verify_real8_3d(a)
    implicit none
    real(dp),dimension(:,:,:) :: a
    logical :: parallel_halo_verify_real8_3d
    parallel_halo_verify_real8_3d = .true.
  end function parallel_halo_verify_real8_3d

#ifdef _USE_MPI_WITH_SLAP
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
    integer, intent(in) :: my_main_rank  ! rank of the master task (ignored for parallel_slap)
    integer :: ierror 
    ! begin
    comm = my_comm
    call mpi_comm_size(comm,tasks,ierror)
    call mpi_comm_rank(comm,this_rank,ierror)
    main_task = .true. !For parallel_slap, each node duplicates all of the calculations.
  end subroutine parallel_set_info

#else
  subroutine parallel_initialise
    implicit none
  end subroutine parallel_initialise

  subroutine parallel_set_info(my_comm, my_main_rank)
    implicit none
    integer, intent(in) :: my_comm       ! CISM's global communicator (IGNORED)
    integer, intent(in) :: my_main_rank  ! rank of the master task (IGNORED)
  end subroutine parallel_set_info

#endif

  subroutine parallel_print_integer_2d(name,values)
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(values,2)
       do i = 1,size(values,1)
          write(u,*) j,i,values(i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine parallel_print_integer_2d

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

  subroutine parallel_print_real8_2d(name,values)
    implicit none
    character(*) :: name
    real(dp),dimension(:,:) :: values
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(values,2)
       do i = 1,size(values,1)
          write(u,*) j,i,values(i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine parallel_print_real8_2d

  subroutine parallel_print_real8_3d(name,values)
    implicit none
    character(*) :: name
    real(dp),dimension(:,:,:) :: values
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(values,3)
       do i = 1,size(values,2)
          write(u,'(2i6,100g15.5e3)') j,i,values(:,i,j)
       end do
       write(u,'()')
    end do
    close(u)
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
    real(sp) :: values
    ! begin
    if (main_task) parallel_put_att_real4 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4)
  end function parallel_put_att_real4

  function parallel_put_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real4_1d,varid
    character(len=*) :: name
    real(sp),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real4_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4_1d)
  end function parallel_put_att_real4_1d

  function parallel_put_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8,varid
    character(len=*) :: name
    real(dp) :: values
    ! begin
    if (main_task) parallel_put_att_real8 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8)
  end function parallel_put_att_real8

  function parallel_put_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8_1d,varid
    character(len=*) :: name
    real(dp),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real8_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8_1d)
  end function parallel_put_att_real8_1d

  function parallel_put_var_real4(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real4,varid
    integer,dimension(:) :: start
    real(sp) :: values
    ! begin
    if (main_task) parallel_put_var_real4 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real4)
  end function parallel_put_var_real4

  function parallel_put_var_real8(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8,varid
    integer,dimension(:) :: start
    real(dp) :: values
    ! begin
    if (main_task) parallel_put_var_real8 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real8)
  end function parallel_put_var_real8

  function parallel_put_var_real8_1d(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8_1d,varid
    integer,dimension(:),optional :: start
    real(dp),dimension(:) :: values
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

! ------------------------------------------
! functions for parallel_reduce_sum interface
! ------------------------------------------
  function parallel_reduce_sum_integer(x)
    ! Sum x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    integer :: x, parallel_reduce_sum_integer

    parallel_reduce_sum_integer = x
    return
  end function parallel_reduce_sum_integer

  function parallel_reduce_sum_real4(x)
    ! Sum x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    real(sp) :: x, parallel_reduce_sum_real4

    parallel_reduce_sum_real4 = x
    return
  end function parallel_reduce_sum_real4

  function parallel_reduce_sum_real8(x)
    ! Sum x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    real(dp) :: x, parallel_reduce_sum_real8

    parallel_reduce_sum_real8 = x
    return
  end function parallel_reduce_sum_real8

  function parallel_reduce_sum_real8_nvar(x)
    ! Sum x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    real(dp) :: x(:)
    real(dp), dimension(size(x)) :: parallel_reduce_sum_real8_nvar

    parallel_reduce_sum_real8_nvar(:) = x(:)
    return
  end function parallel_reduce_sum_real8_nvar

! ------------------------------------------
! functions for parallel_reduce_max interface
! ------------------------------------------
  function parallel_reduce_max_integer(x)
    ! Max x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    integer :: x, parallel_reduce_max_integer

    parallel_reduce_max_integer = x
    return
  end function parallel_reduce_max_integer

  function parallel_reduce_max_real4(x)
    ! Max x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    real(sp) :: x, parallel_reduce_max_real4

    parallel_reduce_max_real4 = x
    return
  end function parallel_reduce_max_real4

  function parallel_reduce_max_real8(x)
    ! Max x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    real(dp) :: x, parallel_reduce_max_real8

    parallel_reduce_max_real8 = x
    return
  end function parallel_reduce_max_real8

! ------------------------------------------
! routines for parallel_reduce_maxloc interface
! ------------------------------------------
  subroutine parallel_reduce_maxloc_integer(xin, xout, xprocout)
    ! Max x across all of the nodes and its proc number
    ! In parallel_slap mode just return x.
    implicit none
    integer, intent(in) :: xin         ! variable to reduce
    integer, intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    xout = xin
    xprocout = this_rank
  end subroutine parallel_reduce_maxloc_integer

  subroutine parallel_reduce_maxloc_real4(xin, xout, xprocout)
    ! Max x across all of the nodes and its proc number
    ! In parallel_slap mode just return x.
    implicit none
    real(sp), intent(in) :: xin         ! variable to reduce
    real(sp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    xout = xin
    xprocout = this_rank
  end subroutine parallel_reduce_maxloc_real4

  subroutine parallel_reduce_maxloc_real8(xin, xout, xprocout)
    ! Max x across all of the nodes and its proc number
    ! In parallel_slap mode just return x.
    implicit none
    real(dp), intent(in) :: xin         ! variable to reduce
    real(dp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    xout = xin
    xprocout = this_rank
  end subroutine parallel_reduce_maxloc_real8

! ------------------------------------------
! functions for parallel_reduce_min interface
! ------------------------------------------
  function parallel_reduce_min_integer(x)
    ! Min x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    integer :: x, parallel_reduce_min_integer

    parallel_reduce_min_integer = x
    return
  end function parallel_reduce_min_integer

  function parallel_reduce_min_real4(x)
    ! Min x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    real(sp) :: x, parallel_reduce_min_real4

    parallel_reduce_min_real4 = x
    return
  end function parallel_reduce_min_real4

  function parallel_reduce_min_real8(x)
    ! Min x across all of the nodes.
    ! In parallel_slap mode just return x.
    implicit none
    real(dp) :: x, parallel_reduce_min_real8

    parallel_reduce_min_real8 = x
    return
  end function parallel_reduce_min_real8

! ------------------------------------------
! routines for parallel_reduce_minloc interface
! ------------------------------------------
  subroutine parallel_reduce_minloc_integer(xin, xout, xprocout)
    ! Min x across all of the nodes and its proc number
    ! In parallel_slap mode just return x.
    implicit none
    integer, intent(in) :: xin         ! variable to reduce
    integer, intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    xout = xin
    xprocout = this_rank
  end subroutine parallel_reduce_minloc_integer

  subroutine parallel_reduce_minloc_real4(xin, xout, xprocout)
    ! Min x across all of the nodes and its proc number
    ! In parallel_slap mode just return x.
    implicit none
    real(sp), intent(in) :: xin         ! variable to reduce
    real(sp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    xout = xin
    xprocout = this_rank
  end subroutine parallel_reduce_minloc_real4

  subroutine parallel_reduce_minloc_real8(xin, xout, xprocout)
    ! Min x across all of the nodes and its proc number
    ! In parallel_slap mode just return x.
    implicit none
    real(dp), intent(in) :: xin         ! variable to reduce
    real(dp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    xout = xin
    xprocout = this_rank
  end subroutine parallel_reduce_minloc_real8


  subroutine parallel_show_minmax(label,values)
    implicit none
    character(*) :: label
    real(dp),dimension(:,:,:) :: values
    ! begin
    print *,label,minval(values),maxval(values)
  end subroutine parallel_show_minmax

  subroutine parallel_stop(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    write(0,*) "STOP in ",file," at line ",line
    stop
  end subroutine parallel_stop

  function parallel_sync(ncid)
    implicit none
    integer :: ncid,parallel_sync
    ! begin
    if (main_task) parallel_sync = nf90_sync(ncid)
    call broadcast(parallel_sync)
  end function parallel_sync

  subroutine staggered_no_penetration_mask(umask, vmask)

    implicit none
    integer,dimension(:,:) :: umask, vmask  ! mask set to 1 wherever the outflow velocity should be zero

    ! initialize the no-penetration masks to 0
    umask(:,:) = 0
    vmask(:,:) = 0

    ! set u velocity mask = 1 at the east global boundary and vertices eastward
    umask(local_ewn-uhalo:,:) = 1

    ! set u velocity mask = 1 at the west global boundary and vertices westward
    umask(:lhalo,:) = 1
    
    ! set v velocity mask = 1 at the north global boundary and vertices northward
    vmask(:,local_nsn-uhalo:) = 1

    ! set v velocity mask = 1 at the south global boundary and vertices southward
    vmask(:,:lhalo) = 1

    call staggered_parallel_halo(umask)
    call staggered_parallel_halo(vmask)

  end subroutine staggered_no_penetration_mask

  subroutine staggered_parallel_halo_extrapolate_integer_2d(a)

    implicit none
    integer,dimension(:,:) :: a
    integer :: i, j

    ! begin

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Extrapolate the staggered field into halo cells along the global boundary.
    ! Currently this is used only for kinbcmask.
    ! Note: The extrapolation region includes locally owned cells along
    !       the north and east boundaries of the global domain.

    ! extrapolate westward
    do i = 1, staggered_lhalo
       a(i, staggered_lhalo+1:size(a,2)-staggered_uhalo-1) = &
            a(staggered_lhalo+1, staggered_lhalo+1:size(a,2)-staggered_uhalo-1)
    enddo

    ! extrapolate eastward
    do i = size(a,1)-staggered_uhalo, size(a,1)
       a(i, staggered_lhalo+1:size(a,2)-staggered_uhalo-1) = &
            a(size(a,1)-staggered_uhalo-1, staggered_lhalo+1:size(a,2)-staggered_uhalo-1)
    enddo

    ! extrapolate southward
    do j = 1, staggered_lhalo
       a(1:size(a,1), j) = a(1:size(a,1), staggered_lhalo+1)
    enddo

    ! extrapolate northward
    do j = size(a,2)-staggered_uhalo, size(a,2)
       a(1:size(a,1), j) = a(1:size(a,1), size(a,2)-staggered_uhalo-1)
    enddo

  end subroutine staggered_parallel_halo_extrapolate_integer_2d


  subroutine staggered_parallel_halo_extrapolate_real8_2d(a)

    implicit none
    real(dp),dimension(:,:) :: a
    integer :: i, j

    ! begin

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Extrapolate the staggered field into halo cells along the global boundary.
    ! Currently this is used only for kinbcmask.
    ! Note: The extrapolation region includes locally owned cells along
    !       the north and east boundaries of the global domain.

    ! extrapolate westward
    do i = 1, staggered_lhalo
       a(i, staggered_lhalo+1:size(a,2)-staggered_uhalo-1) = &
            a(staggered_lhalo+1, staggered_lhalo+1:size(a,2)-staggered_uhalo-1)
    enddo

    ! extrapolate eastward
    do i = size(a,1)-staggered_uhalo, size(a,1)
       a(i, staggered_lhalo+1:size(a,2)-staggered_uhalo-1) = &
            a(size(a,1)-staggered_uhalo-1, staggered_lhalo+1:size(a,2)-staggered_uhalo-1)
    enddo

    ! extrapolate southward
    do j = 1, staggered_lhalo
       a(1:size(a,1), j) = a(1:size(a,1), staggered_lhalo+1)
    enddo

    ! extrapolate northward
    do j = size(a,2)-staggered_uhalo, size(a,2)
       a(1:size(a,1), j) = a(1:size(a,1), size(a,2)-staggered_uhalo-1)
    enddo

  end subroutine staggered_parallel_halo_extrapolate_real8_2d


  subroutine staggered_parallel_halo_integer_2d(a)

    implicit none
    integer,dimension(:,:) :: a

    integer,dimension(staggered_lhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: ecopy
    integer,dimension(staggered_uhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: wcopy
    integer,dimension(size(a,1),staggered_lhalo) :: ncopy
    integer,dimension(size(a,1),staggered_uhalo) :: scopy

    ! begin

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    wcopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
       a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
         1+staggered_lhalo:size(a,2)-staggered_uhalo)

    ecopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
       a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
         1+staggered_lhalo:size(a,2)-staggered_uhalo)

    a(size(a,1)-staggered_uhalo+1:size(a,1), 1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
       wcopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo)

    a(1:staggered_lhalo, 1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
       ecopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo)

    scopy(:,:) = a(:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
    ncopy(:,:) = a(:, size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)

    a(:, size(a,2)-staggered_uhalo+1:size(a,2)) = scopy(:,:)
    a(:, 1:staggered_lhalo) = ncopy(:,:)

  end subroutine staggered_parallel_halo_integer_2d


  subroutine staggered_parallel_halo_integer_3d(a)

    implicit none
    integer,dimension(:,:,:) :: a

    integer,dimension(size(a,1),staggered_lhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: ecopy
    integer,dimension(size(a,1),staggered_uhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: wcopy
    integer,dimension(size(a,1),size(a,2),staggered_lhalo) :: ncopy
    integer,dimension(size(a,1),size(a,2),staggered_uhalo) :: scopy

    ! begin

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1 .or. size(a,3)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    wcopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
       a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
           1+staggered_lhalo:size(a,3)-staggered_uhalo)

    ecopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
       a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
           1+staggered_lhalo:size(a,3)-staggered_uhalo)

    a(:, size(a,2)-staggered_uhalo+1:size(a,2), 1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
       wcopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo)

    a(:, 1:staggered_lhalo, 1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
       ecopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo)

    scopy(:,:,:) = a(:,:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
    ncopy(:,:,:) = a(:,:, size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)

    a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = scopy(:,:,:)
    a(:,:,1:staggered_lhalo) = ncopy(:,:,:)

  end subroutine staggered_parallel_halo_integer_3d


  subroutine staggered_parallel_halo_real8_2d(a)

    implicit none
    real(dp),dimension(:,:) :: a

    real(dp),dimension(staggered_lhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: ecopy
    real(dp),dimension(staggered_uhalo,size(a,2)-staggered_lhalo-staggered_uhalo) :: wcopy
    real(dp),dimension(size(a,1),staggered_lhalo) :: ncopy
    real(dp),dimension(size(a,1),staggered_uhalo) :: scopy

    ! begin

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    wcopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
       a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
         1+staggered_lhalo:size(a,2)-staggered_uhalo)

    ecopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo) = &
       a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
         1+staggered_lhalo:size(a,2)-staggered_uhalo)

    a(size(a,1)-staggered_uhalo+1:size(a,1), 1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
       wcopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo)

    a(1:staggered_lhalo, 1+staggered_lhalo:size(a,2)-staggered_uhalo) = &
       ecopy(:, 1:size(a,2)-staggered_lhalo-staggered_uhalo)

    scopy(:,:) = a(:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
    ncopy(:,:) = a(:, size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)

    a(:, size(a,2)-staggered_uhalo+1:size(a,2)) = scopy(:,:)
    a(:, 1:staggered_lhalo) = ncopy(:,:)

  end subroutine staggered_parallel_halo_real8_2d


  subroutine staggered_parallel_halo_real8_3d(a)

    implicit none
    real(dp),dimension(:,:,:) :: a

    real(dp),dimension(size(a,1),staggered_lhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: ecopy
    real(dp),dimension(size(a,1),staggered_uhalo,size(a,3)-staggered_lhalo-staggered_uhalo) :: wcopy
    real(dp),dimension(size(a,1),size(a,2),staggered_lhalo) :: ncopy
    real(dp),dimension(size(a,1),size(a,2),staggered_uhalo) :: scopy

    ! begin

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1 .or. size(a,3)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    wcopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
       a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
           1+staggered_lhalo:size(a,3)-staggered_uhalo)

    ecopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo) = &
       a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
           1+staggered_lhalo:size(a,3)-staggered_uhalo)

    a(:, size(a,2)-staggered_uhalo+1:size(a,2), 1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
       wcopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo)

    a(:, 1:staggered_lhalo, 1+staggered_lhalo:size(a,3)-staggered_uhalo) = &
       ecopy(:,:, 1:size(a,3)-staggered_lhalo-staggered_uhalo)

    scopy(:,:,:) = a(:,:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
    ncopy(:,:,:) = a(:,:, size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)

    a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = scopy(:,:,:)
    a(:,:,1:staggered_lhalo) = ncopy(:,:,:)

  end subroutine staggered_parallel_halo_real8_3d

!WHL - New subroutine for 4D arrays
  subroutine staggered_parallel_halo_real8_4d(a)

    ! Implements a staggered grid halo update for a 4D field.
    ! This subroutine is used for the 4D arrays that hold matrix entries.

    ! As the grid is staggered, the array 'a' is one smaller in x and y dimensions than an unstaggered array.
    ! The vertical dimension is assumed to precede the i and j indices, i.e., a(:,k,i,j).
    ! The first dimension holds matrix elements for a single row.

    implicit none
    real(dp),dimension(:,:,:,:) :: a

    real(dp),dimension(size(a,1),size(a,2),staggered_lhalo,size(a,4)-staggered_lhalo-staggered_uhalo) :: ecopy
    real(dp),dimension(size(a,1),size(a,2),staggered_uhalo,size(a,4)-staggered_lhalo-staggered_uhalo) :: wcopy
    real(dp),dimension(size(a,1),size(a,2),size(a,3),staggered_lhalo) :: ncopy
    real(dp),dimension(size(a,1),size(a,2),size(a,3),staggered_uhalo) :: scopy

    ! begin

    ! Confirm staggered array
    if (size(a,3)/=local_ewn-1 .or. size(a,4)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    wcopy(:,:,:, 1:size(a,4)-staggered_lhalo-staggered_uhalo) = &
       a(:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             1+staggered_lhalo:size(a,4)-staggered_uhalo)

    ecopy(:,:,:, 1:size(a,4)-staggered_lhalo-staggered_uhalo) = &
       a(:,:,size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo, &
             1+staggered_lhalo:size(a,4)-staggered_uhalo)

    a(:,:, size(a,3)-staggered_uhalo+1:size(a,3), 1+staggered_lhalo:size(a,4)-staggered_uhalo) = &
       wcopy(:,:,:, 1:size(a,4)-staggered_lhalo-staggered_uhalo)

    a(:,:, 1:staggered_lhalo, 1+staggered_lhalo:size(a,4)-staggered_uhalo) = &
       ecopy(:,:,:, 1:size(a,4)-staggered_lhalo-staggered_uhalo)

    scopy(:,:,:,:) = a(:,:,:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
    ncopy(:,:,:,:) = a(:,:,:, size(a,4)-staggered_uhalo-staggered_lhalo+1:size(a,4)-staggered_uhalo)

    a(:,:,:,size(a,4)-staggered_uhalo+1:size(a,4)) = scopy(:,:,:,:)
    a(:,:,:,1:staggered_lhalo) = ncopy(:,:,:,:)

  end subroutine staggered_parallel_halo_real8_4d

end module parallel
