!
!-------------------------------------------------------------------------------
! dynamics - physics coupling for elevation classes
!-------------------------------------------------------------------------------
module ec_coupling

  use shr_kind_mod,   only: r8 => shr_kind_r8, SHR_KIND_CL
  use constituents,   only: pcnst
  use elev_classes,   only: max_elevation_classes
  use ppgrid,         only: begchunk, endchunk
  use phys_grid,      only: grid_chunk_s, grid_chunk_e
  use physics_types,  only: physics_state

  implicit none
  private
  save

  ! phys_column_t represents all elevation class info for one physics column
  !       A physics column corresponds to the grid cell scale
  ! elevation classes are packed in physics-column order
  ! ec_loc specifies the chunk and column for all the elevation classes
  !       which make up this physics column (elevation class set).
  !    An elevation class set is all the elevation classes making up one
  !          physics grid cell (typically the same as a GLL cell).
  !    The first dimension of ec_loc matches area and elevation
  type, public :: phys_column_t
    integer                          :: num_elevation_classes
    real(r8),            allocatable :: area(:)                ! fraction?
    real(r8),            allocatable :: elevation(:)           ! m
    ! The second dimension of ec_loc are for the chunk index and column
    integer,             allocatable :: ec_loc(:,:)            ! (#ecs, 2)
    integer                          :: grid_mean_loc(2)
  end type phys_column_t

  ! Public module variables
  ! ec_sets holds information about the elevation classes on a PE (task)
  !    A single task can hold one or more sets of elevation classes.
  !    A set of elevation classes is all the elevation classes making up one
  !          physics grid cell (typically the same as a GLL cell).
  !    
  type(phys_column_t), allocatable :: ec_sets(:) ! # EC sets on this PE (task)

  ! We need to keep track of the previous grid-mean dynamics state
  type(physics_state), allocatable :: prev_dyn_state(:)

  ! Public interface functions
  public avg_elevation_classes_to_phys_state
  public dyn_state_to_elevation_classes
  public elevation_classes_to_dyn_tend

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !---------------------------------------------------------------------------
  !
  !  avg_elevation_classes_to_phys_state
  !
  !  Average a set of elevation class state or tendencies to the grid cell mean
  !
  !---------------------------------------------------------------------------
  subroutine avg_elevation_classes_to_phys_state(phys_state)
    use ppgrid,         only: pcols, pver

    ! Dummy arguments
    type(physics_state), intent(inout) :: phys_state(begchunk:grid_chunk_e)

    ! Local variables
    integer                         :: set, k, ic
    integer                         :: gchnk, lchnk, gcol, lcol
    real(r8)                        :: accum(pver)
    real(r8)                        :: area, atemp

    ! Loop over our ec_sets, averaging the elevation class columns
    do set = 1, size(ec_sets, 1)
      accum(:) = 0._r8
      area = 0._r8
      if (ec_sets(set)%num_elevation_classes > 1) then
        ! If # classes == 1, then average and EC column are the same
        do ic = 1, ec_sets(set)%num_elevation_classes
          lchnk = ec_sets(set)%ec_loc(ic, 1)
          lcol = ec_sets(set)%ec_loc(ic, 2)
          atemp = ec_sets(set)%area(ic)
          ! Average physics state U (weighted by area?)
          accum(:) = accum(:) + (phys_state(lchnk)%u(lcol,:) * atemp)
          area = area + atemp
        end do ! End loop over elevation classes
        !$omp parallel do private(k)
        do k = 1, pver
          gchnk = ec_sets(set)%grid_mean_loc(1)
          gcol = ec_sets(set)%grid_mean_loc(2)
          phys_state(gchnk)%u(gcol, k) = accum(k) / area
        end do
      end if
    end do
  end subroutine avg_elevation_classes_to_phys_state

  !---------------------------------------------------------------------------
  !
  !  dyn_state_to_elevation_classes
  !
  !  Create a dynamics tendency and apply it to phys_state
  !
  !---------------------------------------------------------------------------
  subroutine dyn_state_to_elevation_classes(dt, curr_dyn_state, phys_state)
    use physics_types, only: physics_state_copy
    use physconst,     only: rair, gravit
    use hycoef,        only: hyam, hybm, hyai, hybi, ps0

    ! Dummy arguments
    real(r8),            intent(in)    :: dt
    type(physics_state), intent(in)    :: curr_dyn_state(grid_chunk_s:grid_chunk_e)
    type(physics_state), intent(inout) :: phys_state(begchunk:grid_chunk_e)

    ! Local variables
    integer                            :: set, m, k, ic, nlev
    integer                            :: gchnk, lchnk, gcol, lcol
    real(r8), allocatable              :: tend(:)
    real(r8)                           :: work1(max_elevation_classes), work2

    nlev = size(curr_dyn_state(grid_chunk_s)%u, 2)
    allocate(tend(nlev))

    ! Traverse the EC sets and apply dynamics update
    do set = 1, size(ec_sets, 1)
      gchnk = ec_sets(set)%grid_mean_loc(1)
      gcol = ec_sets(set)%grid_mean_loc(2)
      ! Calculate the T tendency
      if (allocated(prev_dyn_state)) then
        tend(:) = curr_dyn_state(gchnk)%T(gcol,:) - prev_dyn_state(gchnk)%T(gcol,:)
      else
        tend(:) = curr_dyn_state(gchnk)%T(gcol,:)
      end if
      ! Calculate the work
      work2 = 0.0_r8
      do ic = 1, ec_sets(set)%num_elevation_classes
        lchnk = ec_sets(set)%ec_loc(ic, 1)
        lcol = ec_sets(set)%ec_loc(ic, 2)
        work1(ic) = exp(-((gravit*ec_sets(set)%elevation(ic)) -              &
             curr_dyn_state(gchnk)%phis(gcol))                /              &
             (rair * curr_dyn_state(gchnk)%T(gcol,nlev)))
        work2 = work2 + ec_sets(set)%area(ic) * work1(ic)
      end do
      ! Update elevation classes
      do ic = 1, ec_sets(set)%num_elevation_classes
        lchnk = ec_sets(set)%ec_loc(ic, 1)
        lcol = ec_sets(set)%ec_loc(ic, 2)
        ! Apply U (equally to each elevation class)
        phys_state(lchnk)%u(lcol,:) = curr_dyn_state(gchnk)%u(gcol,:)
        ! Apply V (equally to each elevation class)
        phys_state(lchnk)%v(lcol,:) = curr_dyn_state(gchnk)%v(gcol,:)
        ! Apply T dynamics tendency
        phys_state(lchnk)%T(lcol,:) = phys_state(lchnk)%T(lcol,:) + tend(:)
        ! Apply omega to each elevation class
        phys_state(lchnk)%omega(lcol,:) = curr_dyn_state(gchnk)%omega(gcol,:)
        ! Apply PHIS 
        phys_state(lchnk)%phis(lcol) = gravit*ec_sets(set)%elevation(ic)
        ! Apply work to PS
        phys_state(lchnk)%ps(lcol) = curr_dyn_state(gchnk)%ps(gcol) * work1(ic) / work2
      end do
    end do

    ! Apply Q dynamics tendency
    do m = 1, pcnst
      ! Traverse the EC sets and apply dynamics update
      do set = 1, size(ec_sets, 1)
        gchnk = ec_sets(set)%grid_mean_loc(1)
        gcol = ec_sets(set)%grid_mean_loc(2)
        ! Find the tracer tendency
        if (allocated(prev_dyn_state)) then
          tend(:) = curr_dyn_state(gchnk)%q(gcol,:,m) - prev_dyn_state(gchnk)%q(gcol,:,m) ! change due to dynamics
        else
          tend(:) = curr_dyn_state(gchnk)%q(gcol,:,m) ! just use current value
        end if
        ! Update elevation classes
        do ic = 1, ec_sets(set)%num_elevation_classes
          lchnk = ec_sets(set)%ec_loc(ic, 1)
          lcol = ec_sets(set)%ec_loc(ic, 2)
          do k = 1, nlev
!!XXgoldyXX: No use of area?
            if(tend(k) <= 0._r8) then
              if (allocated(prev_dyn_state)) then
                if(prev_dyn_state(gchnk)%q(gcol,k,m) > 0._r8) then
                  phys_state(lchnk)%q(lcol,k,m) = phys_state(lchnk)%q(lcol,k,m) + &
                       tend(k) * phys_state(lchnk)%q(lcol,k,m) / prev_dyn_state(gchnk)%q(gcol,k,m)
                else
                  ! prev_dyn_state exists but the tendency was <= 0
!!XXgoldyXX: Really?
                  phys_state(lchnk)%q(lcol,k,m) = phys_state(lchnk)%q(lcol,k,m)
                end if
              else
                ! prev_dyn_state does not exist, just use the current
                phys_state(lchnk)%q(lcol,k,m) = curr_dyn_state(gchnk)%q(gcol,k,m)*ec_sets(set)%area(lcol)
              end if
            else
!!XXgoldyXX: Really?
              phys_state(lchnk)%q(lcol,k,m) = phys_state(lchnk)%q(lcol,k,m)+tend(k)
            end if
          end do
        end do
      end do
    end do

    ! Store the current dynamics state for the next iteration
    if (.not. allocated(prev_dyn_state)) then
      allocate(prev_dyn_state(grid_chunk_s:grid_chunk_e))
    end if
    do ic = grid_chunk_s, grid_chunk_e
      call physics_state_copy(curr_dyn_state(ic), prev_dyn_state(ic))
    end do
    ! Cleanup
    deallocate(tend)

  end subroutine dyn_state_to_elevation_classes

  !---------------------------------------------------------------------------
  !
  !  elevation_classes_to_dyn_tend
  !
  !  Average the physics tendencies in phys_state into tendency for dynamics
  !
  !---------------------------------------------------------------------------
  subroutine elevation_classes_to_dyn_tend(phys_state, phys_tend)
    use physics_types,  only: physics_state, physics_tend

    ! Dummy arguments
    type(physics_state),     intent(inout) :: phys_state(begchunk:grid_chunk_e)
    type(physics_tend),      intent(inout) :: phys_tend(begchunk:grid_chunk_e)

    ! Local variables
    integer                            :: set, ic, k, m, nlev
    integer                            :: gchnk, lchnk, gcol, lcol
    real(r8), allocatable              :: tend(:,:)
    real(r8)                           :: area, atemp
    integer, parameter                 :: uind = 1 ! u tend index
    integer, parameter                 :: vind = 2 ! v tend index
    integer, parameter                 :: tind = 3 ! t tend index
    integer, parameter                 :: lind = 3 ! last tend index before Q

    nlev = size(phys_tend(begchunk)%dudt, 2)
    allocate(tend(nlev,1:pcnst+lind))

    ! Traverse the EC sets and apply dynamics update
    do set = 1, size(ec_sets, 1)
      gchnk = ec_sets(set)%grid_mean_loc(1)
      gcol = ec_sets(set)%grid_mean_loc(2)
      tend(:,:) = 0._r8
      area = 0._r8
      do ic = 1, ec_sets(set)%num_elevation_classes
        lchnk = ec_sets(set)%ec_loc(ic, 1)
        lcol = ec_sets(set)%ec_loc(ic, 2)
        atemp = ec_sets(set)%area(ic)
        tend(:,uind) = tend(:,uind) + (phys_tend(lchnk)%dudt(lcol,:) * atemp)
        tend(:,vind) = tend(:,vind) + (phys_tend(lchnk)%dvdt(lcol,:) * atemp)
        tend(:,tind) = tend(:,tind) + (phys_tend(lchnk)%dtdt(lcol,:) * atemp)
        do m = 1, pcnst
          tend(:,m+lind) = tend(:,m+lind) + (phys_state(lchnk)%Q(lcol,:,m) * atemp)
        end do
        area = area + atemp
      end do ! End loop over elevation classes
      !$omp parallel do private(k,m)
      do k = 1, nlev
        ! Average U physics tendency (weighted by area?)
        phys_tend(gchnk)%dudt(gcol,k) = tend(k,uind) / area
        phys_tend(gchnk)%dvdt(gcol,k) = tend(k,vind) / area
        phys_tend(gchnk)%dtdt(gcol,k) = tend(k,tind) / area
        do m = 1, pcnst
          phys_state(gchnk)%Q(gcol,k,m) = tend(k,m+lind) / area
        end do
      end do
    end do ! End loop over sets

    ! Cleanup
    deallocate(tend)

  end subroutine elevation_classes_to_dyn_tend

end module ec_coupling
