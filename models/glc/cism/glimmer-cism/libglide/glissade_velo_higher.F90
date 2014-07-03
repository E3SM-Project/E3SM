!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! THE FUNCTIONALITY OF THIS MODULE HAS BEEN REMOVED FROM THE RELEASE CODE. 
!
! ONLY STUBS REMAIN HERE. TRYING TO CALL ANY OF THESE ROUTINES WILL RESULT IN A FATAL ERROR.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_velo_higher

    use glimmer_global, only: sp, dp
    use glimmer_log

    implicit none

    private
    public :: glissade_velo_higher_init, glissade_velo_higher_solve
    contains

!****************************************************************************

  subroutine glissade_velo_higher_init

    call write_log("glissade_velo_higher_init has been removed from the release code", GM_FATAL)

  end subroutine glissade_velo_higher_init

!****************************************************************************

  subroutine glissade_velo_higher_solve(nx,         ny,           &
                                        nz,         sigma,        &
                                        nhalo,                    &
                                        dx,         dy,           &
                                        thck,       usrf,         &
                                        topg,       eus,          &
!!                                        stagthck, stagusrf,     &
                                        thklim,                   &
                                        flwa,                     &
                                        uvel,     vvel,           &
!                                        beta,                    &
!                                        whichbabc,               &
                                        whichefvs,                &
                                        whichresid,               &
                                        whichnonlinear,           &
                                        whichsparse,              &
                                        whichapprox)

    integer, intent(in) ::   &
       nx, ny,               &  ! number of grid cells in each direction
       nz,                   &  ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)
       nhalo                    ! number of rows/columns of halo cells
 
    real(dp), dimension(:), intent(in) :: &
       sigma

    real(dp), intent(in) ::  &
       dx,  dy                  ! grid cell length and width (m)
                                ! assumed to have the same value for each grid cell

    real(dp), dimension(:,:), intent(in) ::  &
       thck,                 &  ! ice thickness (m)
       usrf,                 &  ! upper surface elevation (m)
       topg                     ! elevation of topography (m)

    real(dp), intent(in) ::   &
       eus                      ! eustatic sea level (m)
                                ! = 0. by default

    real(dp), intent(in) ::   & 
       thklim                 ! minimum ice thickness for active cells (m)

    real (dp), dimension(:,:,:), intent(in) ::  &
       flwa                   ! flow factor in units of Pa^(-n) s^(-1)

    real (dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel             ! velocity components (m/s)

    !TODO - Uncomment and pass in
!!    real (dp), dimension(:,:), intent(in) ::  &
!!       beta                   ! basal traction parameter

    integer, intent(in) :: whichefvs      ! option for effective viscosity calculation 
                                          ! (calculate it or make it uniform)
    integer, intent(in) :: whichresid     ! option for method to use when calculating residual
    integer, intent(in) :: whichnonlinear ! options for which nonlinear method (Picard or JFNK)
    integer, intent(in) :: whichsparse    ! options for which method for doing elliptic solve
                                          ! (BiCG, GMRES, standalone Trilinos, etc.)
    integer, intent(in), optional :: whichapprox    ! options for which Stokes approximation to use
                                                    ! 0 = SIA, 1 = SSA, 2 = Blatter-Pattyn HO
                                                    ! default = 2

    call write_log("glissade_velo_higher_solve has been removed from the release code", GM_FATAL)

  end subroutine glissade_velo_higher_solve

!****************************************************************************

end module glissade_velo_higher

!****************************************************************************
