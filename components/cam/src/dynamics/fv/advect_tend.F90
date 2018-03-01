!----------------------------------------------------------------------
! this module computes the total advection tendencies of advected
! constituents for the finite volume dycore
!----------------------------------------------------------------------
module advect_tend

  use shr_kind_mod, only : r8 => shr_kind_r8

  save
  private

  public :: compute_adv_tends_xyz

  real(r8), allocatable :: adv_tendxyz(:,:,:,:)

contains

  !----------------------------------------------------------------------
  ! computes the total advective tendencies
  ! called twice each time step:
  !   - first call sets the initial mixing ratios
  !   - second call computes and outputs the tendencies
  !----------------------------------------------------------------------
  subroutine compute_adv_tends_xyz( grid, tracer )
    use dynamics_vars, only : T_FVDYCORE_GRID
    use cam_history,   only : outfld
    use time_manager,  only : get_step_size
    use constituents,  only : tottnam

    implicit none

    type (T_FVDYCORE_GRID), intent(in) :: grid
    real (r8) ::  tracer(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km, grid%nq )

    real(r8) :: dt,idt
    integer  :: iq, idim, i, j,ic
    logical  :: init
    real(r8) :: tmpxy(grid%ifirstxy:grid%ilastxy,grid%km)


    init = .false.

    if ( .not. allocated( adv_tendxyz ) ) then
       init = .true.
       allocate( adv_tendxyz(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km, grid%nq ) )
       adv_tendxyz(:,:,:,:) = 0._r8
    endif

!!$    adv_tendxyz(:,:,:,:grid%nq) = q(:,:,:,:grid%nq) - adv_tendxyz(:,:,:,:grid%nq)

    do ic=1,grid%nq
       adv_tendxyz(:,:,:,ic) = tracer(:,:,:,ic) - adv_tendxyz(:,:,:,ic)
    enddo

    if ( .not. init ) then
       dt = get_step_size()
       idt = 1._r8/dt

       do i = 1,grid%nq
          ! call outfld
          do j = grid%jfirstxy, grid%jlastxy
             idim = grid%ilastxy - grid%ifirstxy + 1
             tmpxy(:,:) = adv_tendxyz(:,j,:,i)*idt

             call outfld( tottnam(i),  tmpxy, idim, j)
          enddo
       enddo

       deallocate(adv_tendxyz)
    endif

  end subroutine compute_adv_tends_xyz

end module advect_tend
