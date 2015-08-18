module vdiff_lu_solver

! Low level vertical diffusion module containing tridiagonal solver.

! This module was created solely to share vd_lu_decomp and vd_lu_solve
! between gw_drag and diffusion_solver.
use module_perturb

implicit none
private
save

! Public interface
public :: vd_lu_decomp
public :: vd_lu_solve

public :: lu_decomp

! 8-byte real.
integer, parameter :: r8 = selected_real_kind(12)

! Type to hold the sparse matrix decomposition from vd_lu_decomp.
type :: lu_decomp
   integer :: ncol = 0
   integer :: pver = 0
   ! Upper and lower diagonals.
   real(r8), allocatable :: ca(:,:)
   real(r8), allocatable :: cc(:,:)
   ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
   real(r8), allocatable :: dnom(:,:)
   ! Term in tri-diag. matrix system.
   real(r8), allocatable :: ze(:,:)
contains
  procedure, pass(decomp) :: finalize => lu_decomp_dealloc
end type lu_decomp

! LU decomposition constructor.
interface lu_decomp
   module procedure lu_decomp_alloc
end interface

contains

! ========================================================================!

subroutine vd_lu_decomp(pcols, pver, ncol,                          &
                        ksrf,  kv,   tmpi,   rpdel,  ztodt,  gravit,&
                        cc_top, ntop,  nbot, decomp, cpairv, flb,lchnk)

  !---------------------------------------------------------------------- !
  ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the !
  ! tridiagonal diffusion matrix.                                         !
  ! The diagonal elements (1+ca(k)+cc(k)) are not required by the solver. !
  ! Also determine ze factor and denominator for ze and zf (see solver).  !
  !---------------------------------------------------------------------- !

  ! ---------------------- !
  ! Input-Output Arguments !
  ! ---------------------- !
  integer,  intent(in),optional  :: flb,lchnk
  ! Allocated column and level dimensions.
  integer,  intent(in)  :: pcols, pver
  ! Columns actually used during computation.
  integer,  intent(in)  :: ncol
  ! Top and bottom levels to operate on.
  integer,  intent(in)  :: ntop, nbot

  ! Surface "drag" coefficient. [ kg/s/m2 ]
  real(r8), intent(in)  :: ksrf(pcols)
  ! Vertical diffusion coefficients. [ m2/s ]
  real(r8), intent(in)  :: kv(pcols,pver+1)
  ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
  real(r8), intent(in)  :: tmpi(pcols,pver+1)
  ! 1./pdel (pdel is thickness between interfaces). [ Pa^-1 ]
  real(r8), intent(in)  :: rpdel(pcols,pver)
  ! 2 delta-t [ s ]
  real(r8), intent(in)  :: ztodt
  ! Acceleration due to gravity. [ m/s2 ]
  real(r8), intent(in)  :: gravit
  ! Lower diagonal on top interface (for fixed ubc only).
  real(r8), intent(in)  :: cc_top(pcols)

  ! Output decomposition.
  type(lu_decomp), intent(out) :: decomp

  ! "Variable" specific heat at constant pressure.
  real(r8), intent(in), optional  :: cpairv(pcols,pver)

  ! --------------- !
  ! Local Variables !
  ! --------------- !

  ! Column and level indices.
  integer :: i, k

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  decomp = lu_decomp(ncol, pver)

  ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the
  ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
  ! a combination of ca and cc; they are not required by the solver.

  !
  ! If switch present and set true, then input kv is kvt for use in
  ! diagonal calculations
  !

  if ( present(cpairv) ) then
     do k = nbot-1, ntop, -1
        do i = 1, ncol
           decomp%ca(i,k  ) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k  ) &
                / cpairv(i,k)
           decomp%cc(i,k+1) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k+1) &
                / cpairv(i,k+1)
        end do
     end do
  else
     do k = nbot - 1, ntop, -1
        do i = 1, ncol
           decomp%ca(i,k  ) = kv(i,k+1) * tmpi(i,k+1) * rpdel(i,k  )
           decomp%cc(i,k+1) = kv(i,k+1) * tmpi(i,k+1) * rpdel(i,k+1)

           if(present(lchnk) .and. present(flb)) then
              if(flb==1 .and. icolprnt(lchnk) ==i .and. k==kprnt)then
                 !decomp%ca(i,k  ) = 0.0_r8 * tmpi(i,k+1) * rpdel(i,k  )!BALLI WRONG
                 !decomp%cc(i,k+1) = 0.0_r8 * tmpi(i,k+1) * rpdel(i,k+1)!BALLI WRONG
                 write(202,*)'--->luvdecomp_2:',decomp%ca(i,k),kv(i,k+1),tmpi(i,k+1),rpdel(i,k),rpdel(i,k +1),'WRONG'
              endif
           endif
           
        end do
     end do
  endif

  ! The bottom element of the upper diagonal (ca) is zero (not used).
  ! The subdiagonal (cc) is not needed in the solver.
  decomp%ca(:,nbot) = 0._r8

  ! Calculate e(k). This term is required in the solution of the
  ! tridiagonal matrix as defined by the implicit diffusion equation.

  do i = 1, ncol
     decomp%dnom(i,nbot) = 1._r8/ &
          (1._r8 + decomp%cc(i,nbot) + ksrf(i)*ztodt*gravit*rpdel(i,nbot))
     decomp%ze(i,nbot)   = decomp%cc(i,nbot) * decomp%dnom(i,nbot)
  end do

  do k = nbot - 1, ntop + 1, -1
     do i = 1, ncol
        decomp%dnom(i,k) = 1._r8/ &
             (1._r8 + decomp%ca(i,k) + decomp%cc(i,k) - &
             decomp%ca(i,k)*decomp%ze(i,k+1))

        if(present(lchnk) .and. present(flb)) then
           if(flb==1 .and. icolprnt(lchnk) ==i .and. k==kprnt)then
              write(202,*)'--->luvdecomp_1:',decomp%dnom(i,k),decomp%ca(i,k),decomp%cc(i,k),decomp%ze(i,k+1)
           endif
        endif
        decomp%ze(i,k)   = decomp%cc(i,k) * decomp%dnom(i,k)
     end do
  end do

  do i = 1, ncol
     decomp%dnom(i,ntop) = 1._r8/ &
          (1._r8 + decomp%ca(i,ntop) + cc_top(i) - &
          decomp%ca(i,ntop)*decomp%ze(i,ntop+1))
  end do

end subroutine vd_lu_decomp

! ========================================================================!

subroutine vd_lu_solve(pcols, pver,   ncol, &
                       q,     decomp, ntop, nbot, cd_top, flb,lchnk)
  !-----------------------------------------------------------------------!
  ! Solve the implicit vertical diffusion equation with zero flux         !
  ! boundary conditions. Actual surface fluxes are explicit (rather than  !
  ! implicit) and applied separately. Procedure for solution of the       !
  ! implicit equation follows Richtmyer and Morton (1967,pp 198-200).     !
  !                                                                       !
  ! The equation solved is                                                !
  !                                                                       !
  !     -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k),                 !
  !                                                                       !
  ! where d(k) is the input profile and q(k) is the output profile        !
  !                                                                       !
  ! The solution has the form                                             !
  !                                                                       !
  !     q(k) = ze(k)*q(k-1) + zf(k)                                       !
  !                                                                       !
  !     ze(k) = cc(k) * dnom(k)                                           !
  !                                                                       !
  !     zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)                          !
  !                                                                       !
  !     dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)]                               !
  !             = 1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]                   !
  !                                                                       !
  ! Note that the same routine is used for temperature, momentum and      !
  ! tracers, and that input variables are replaced.                       !
  ! ----------------------------------------------------------------------!

  ! ---------------------- !
  ! Input-Output Arguments !
  ! ---------------------- !
  integer, optional :: flb,lchnk
  ! Allocated column and level dimensions.
  integer,  intent(in)  :: pcols, pver
  ! Columns actually used during computation.
  integer,  intent(in)  :: ncol
  ! Top and bottom levels to operate on.
  integer,  intent(in)  :: ntop, nbot

  ! LU decomposition information.
  type(lu_decomp), intent(in) :: decomp

  ! cc_top * ubc value
  real(r8), intent(in)    :: cd_top(pcols)

  ! Constituent field.
  real(r8), intent(inout) :: q(pcols,pver)

  ! --------------- !
  ! Local Variables !
  ! --------------- !

  ! Term in tri-diag solution.
  real(r8) :: zf(ncol,pver)
  ! Column and level indices.
  integer :: i, k

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  ! Calculate zf(k). Terms zf(k) and ze(k) are required in solution of
  ! tridiagonal matrix defined by implicit diffusion equation.
  ! Note that only levels ntop through nbot need be solved for.

  do i = 1, ncol
     zf(i,nbot) = q(i,nbot)*decomp%dnom(i,nbot)
  end do

  do k = nbot - 1, ntop + 1, -1
     do i = 1, ncol
        zf(i,k) = (q(i,k) + decomp%ca(i,k)*zf(i,k+1)) * decomp%dnom(i,k)
        if(present(lchnk) .and. present(flb)) then
           if(flb==1 .and. icolprnt(lchnk) ==i .and. k==kprnt)then
              !zf(i,k) = (q(i,k) + 0.0_r8*zf(i,k+1)) * 1.0_r8 !BALLI WRONG
              write(202,*)'--->luvdiff_2:',zf(i,k), q(i,k), decomp%ca(i,k),zf(i,k+1),&
                decomp%dnom(i,k)!, 'WRONG!!!'
           endif
        endif
     end do
  end do

  ! Include boundary condition on top element

  k = ntop
  do i = 1, ncol
     zf(i,k) = (q(i,k) + cd_top(i) + decomp%ca(i,k)*zf(i,k+1)) * &
          decomp%dnom(i,k)
  end do

  ! Perform back substitution

  do i = 1, ncol
     q(i,ntop) = zf(i,ntop)
  end do

  do k = ntop + 1, nbot, +1
     do i = 1, ncol
        q(i,k) = zf(i,k) + decomp%ze(i,k)*q(i,k-1)
        if(present(lchnk) .and. present(flb)) then
           if(flb==1 .and.icolprnt(lchnk) ==i .and. k==kprnt)write(202,*)'--->luvdiff_1:',q(i,k),zf(i,k),decomp%ze(i,k), q(i,k-1)
        endif
     end do
  end do

end subroutine vd_lu_solve

! LU decomposition allocation.
pure function lu_decomp_alloc(ncol, pver) result(new_decomp)
  integer, intent(in) :: ncol
  integer, intent(in) :: pver
  type(lu_decomp) :: new_decomp

  new_decomp%ncol = ncol
  new_decomp%pver = pver

  ! Simple allocation with no error checking.
  allocate(new_decomp%ca(ncol,pver))
  allocate(new_decomp%cc(ncol,pver))
  allocate(new_decomp%dnom(ncol,pver))
  allocate(new_decomp%ze(ncol,pver))

end function lu_decomp_alloc

! LU decomposition deallocation.
pure subroutine lu_decomp_dealloc(decomp)
  class(lu_decomp), intent(inout) :: decomp

  decomp%ncol = 0
  decomp%pver = 0

  ! Simple deallocation with no error checking.
  deallocate(decomp%ca)
  deallocate(decomp%cc)
  deallocate(decomp%dnom)
  deallocate(decomp%ze)

end subroutine lu_decomp_dealloc

end module vdiff_lu_solver
