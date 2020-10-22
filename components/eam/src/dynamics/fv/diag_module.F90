!-----------------------------------------------------------------------
!
! BOP
!
! !MODULE: 
module diag_module
!
! !DESCRIPTION
!   Utilities which perform special calculations for diagnostics
!   Currently only compute\_vdot\_grad to calculate total derivative
!
! !PUBLIC MEMBER FUNCTIONS
  public :: compute_vdot_gradp

! !REVISION HISTORY:
!   05.09.10  Rasch   Creation of compute_vdot_gradp
!   05.10.18  Sawyer  Revisions for 2D decomp, placed in module
!   07.01.29  Chen    Removed pft2d calculation for OMGA (is in cd_core)
!
! EOP
!-----------------------------------------------------------------------
  private

CONTAINS

  subroutine compute_vdot_gradp(grid, dt, frac, cx, cy, pexy, omgaxy )

  use shr_kind_mod, only :  r8 => shr_kind_r8
  use dynamics_vars, only : T_FVDYCORE_GRID
#if defined( SPMD )
  use mod_comm, only: mp_send3d, mp_recv3d, &
                      mp_sendirr, mp_recvirr
#endif

  implicit none

! !INPUT PARAMETERS:
  type (T_FVDYCORE_GRID), intent(in) :: grid
  real(r8), intent(in):: dt
  real(r8), intent(in):: frac

  real(r8), intent(in):: cx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
  real(r8), intent(in):: cy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)
  real(r8), target, intent(in)::   &
    pexy(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy) ! P (pascal) at layer edges
  real(r8), target, intent(inout):: &
    omgaxy(grid%ifirstxy:grid%ilastxy,grid%km,grid%jfirstxy:grid%jlastxy) ! vert. press. velocity (pa/sec)

! Local 
  integer :: im       ! dimension in east-west
  integer :: jm       ! dimension in North-South
  integer :: km       ! number of Lagrangian layers
  integer :: jfirst   ! starting latitude index for MPI
  integer :: jlast    ! ending latitude index for MPI
  integer :: kfirst   ! starting level index for MPI
  integer :: klast    ! ending level index for MPI
  integer :: js2g0    ! j==1 not included
  integer :: jn2g0    ! j==jm not included

  real(r8) :: pm(grid%im, grid%jfirst-1:grid%jlast+1)
  real(r8) :: grad(grid%im, grid%jfirst:grid%jlast+1)
  real(r8) :: fac, sum1

  real(r8), pointer :: pe(:,:,:)      ! YZ version of edge pressures
  real(r8), pointer :: omga(:,:,:)    ! YZ version of vert. vel.

  real(r8), parameter :: half =  0.5_r8
  real(r8), parameter :: zero =  0.0_r8

  integer  :: i,j,k

#if defined( SPMD )
  integer  :: iam, dest, src, npr_y, npes_yz
  real(r8) :: penorth(grid%im, grid%kfirst:grid%klast+1)
  real(r8) :: pesouth(grid%im, grid%kfirst:grid%klast+1)
#endif

  im     = grid%im
  jm     = grid%jm
  km     = grid%km
  jfirst = grid%jfirst
  jlast  = grid%jlast
  kfirst = grid%kfirst
  klast  = grid%klast
  js2g0  = grid%js2g0
  jn2g0  = grid%jn2g0

  fac = half / (dt * frac)

#if defined( SPMD )
  if ( grid%twod_decomp == 1 ) then
    allocate(pe(im,kfirst:klast+1,jfirst:jlast))
    allocate(omga(im,kfirst:klast,jfirst:jlast))
    call mp_sendirr( grid%commxy, grid%ikj_xy_to_yz%SendDesc,                  &
                     grid%ikj_xy_to_yz%RecvDesc, omgaxy, omga,                 &
                     modc=grid%modc_dynrun )
    call mp_recvirr( grid%commxy, grid%ikj_xy_to_yz%SendDesc,                  &
                     grid%ikj_xy_to_yz%RecvDesc, omgaxy, omga,                 &
                     modc=grid%modc_dynrun )
    call mp_sendirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                     grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                     modc=grid%modc_dynrun )
    call mp_recvirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                     grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                     modc=grid%modc_dynrun )
  else
    pe => pexy
    omga => omgaxy
  endif
  iam   = grid%iam
  npes_yz   = grid%npes_yz
 if (iam .lt. npes_yz) then
  npr_y = grid%npr_y
  dest  = iam+1
  src   = iam-1
  if ( mod(iam+1,npr_y) == 0 ) dest = -1
  if ( mod(iam,npr_y) == 0 ) src = -1

!
! Have to give more thought to the source and destination for 2D
!
  call mp_send3d(grid%commyz, dest, src, im, km+1, jm,       &
                 1, im, kfirst, klast+1, jfirst, jlast,     &
                 1, im, kfirst, klast+1, jlast, jlast, pe)
  call mp_recv3d(grid%commyz, src, im, km+1, jm,              &
                 1, im, kfirst, klast+1, jfirst-1, jfirst-1, &
                 1, im, kfirst, klast+1, jfirst-1, jfirst-1, pesouth)
  call mp_send3d(grid%commyz, src, dest, im, km+1, jm,       &
                 1, im, kfirst, klast+1, jfirst, jlast,     &
                 1, im, kfirst, klast+1, jfirst, jfirst, pe)
  call mp_recv3d(grid%commyz, dest, im, km+1, jm,            &
                 1, im, kfirst, klast+1, jlast+1, jlast+1,  &
                 1, im, kfirst, klast+1, jlast+1, jlast+1, penorth)
 end if  !  (iam .lt. npes_yz)
#else
  pe => pexy
  omga => omgaxy
#endif

!$omp parallel do private(i,j,k,pm,grad, sum1)
  do k=kfirst,klast

! Compute layer mean p
     do j=jfirst,jlast
        do i=1,im
           pm(i,j) = half * ( pe(i,k,j) + pe(i,k+1,j) )
        enddo
     enddo

#if defined( SPMD )
     if ( jfirst/=1 ) then
        do i=1,im
           pm(i,jfirst-1) = half * ( pesouth(i,k) + pesouth(i,k+1))
        enddo
     endif

     if ( jlast/=jm ) then
        do i=1,im
           pm(i,jlast+1) = half * ( penorth(i,k) + penorth(i,k+1))
        enddo
     endif
#endif

     do j=js2g0,jn2g0
           i=1
           grad(i,j) = fac * cx(i,j,k) * (pm(i,j)-pm(im,j)) 
        do i=2,im
           grad(i,j) = fac * cx(i,j,k) * (pm(i,j)-pm(i-1,j)) 
        enddo
     enddo

     do j=js2g0,jn2g0
        do i=1,im-1
           omga(i,k,j) = omga(i,k,j) + grad(i,j) + grad(i+1,j)
        enddo
           i=im
           omga(i,k,j) = omga(i,k,j) + grad(i,j) + grad(1,j)
     enddo

     do j=js2g0,min(jm,jlast+1)
        do i=1,im
           grad(i,j) = fac * cy(i,j,k) * (pm(i,j)-pm(i,j-1)) 
        enddo
     enddo

     do j=js2g0,jn2g0
        do i=1,im
           omga(i,k,j) = omga(i,k,j) + grad(i,j) + grad(i,j+1)
        enddo
     enddo

! Note: Since V*grad(P) at poles are harder to compute accurately we use the average of sourding points
!       to be used as input to physics.

     if ( jfirst==1 ) then
        sum1 = zero
        do i=1,im
           sum1 = sum1 + omga(i,k,2)
        enddo
        sum1 = sum1 / real(im,r8)
        do i=1,im
           omga(i,k,1) = sum1
        enddo
     endif

     if ( jlast==jm ) then
        sum1 = zero
        do i=1,im
           sum1 = sum1 + omga(i,k,jm-1)
        enddo
        sum1 = sum1 / real(im,r8)
        do i=1,im
           omga(i,k,jm) = sum1
        enddo
     endif
  enddo

#if defined( SPMD)
  if ( grid%twod_decomp == 1 ) then
!
! Transpose back to XY  (if 1D, the changes to omgaxy were made in place)
!
    call mp_sendirr( grid%commxy, grid%ikj_yz_to_xy%SendDesc,                  &
                     grid%ikj_yz_to_xy%RecvDesc, omga, omgaxy,                 &
                     modc=grid%modc_dynrun )
    call mp_recvirr( grid%commxy, grid%ikj_yz_to_xy%SendDesc,                  &
                     grid%ikj_yz_to_xy%RecvDesc, omga, omgaxy,                 &
                     modc=grid%modc_dynrun )
    deallocate( pe )
    deallocate( omga )
  endif
#endif

  end subroutine compute_vdot_gradp

end module diag_module
