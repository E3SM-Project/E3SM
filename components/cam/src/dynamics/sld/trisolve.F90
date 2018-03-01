subroutine trisolve(n       ,atri    ,btrie   ,ctri    ,dtri    , &
                    dnmn    )
!-----------------------------------------------------------------------
!
! Purpose:
! Solve tri-diagonal system of semi-implicit diverence equations in
! Normal Mode space.
! NOTE:  Storage in the vectors assumed to be along columns ("N")
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: n          ! length of complex vector
  real(r8), intent(in)   :: atri (2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(in)   :: btrie(2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(in)   :: ctri (2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(in)   :: dtri (2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(out)  :: dnmn (2*n) ! divergence solution in Normal Mode space.
!
!---------------------------Local workspace-----------------------------
!
  integer nn             ! n-wavenumber index
  integer nnm2           ! nn-2
  integer nnp2           ! nn+2
  real(r8) tmp  (2)      ! tmp workspace (complex)
  real(r8) denom(2)      ! tmp workspace (complex)
  real(r8) numer(2)      ! tmp workspace (complex)
  real(r8) e    (2*pmax) ! tmp space in solving tri-diag matrix
  real(r8) f    (2*pmax) ! tmp space in solving tri-diag matrix
  real(r8) denom1        ! tmp space in solving tri-diag matrix
!
!-----------------------------------------------------------------------
!
  denom1 =  btrie(1)*btrie(1) + btrie(2)*btrie(2)
  if(n .gt. 1) then
     e(1) = (atri(1)*btrie(1) + atri(2)*btrie(2))/denom1
     e(2) = (atri(2)*btrie(1) - atri(1)*btrie(2))/denom1
  endif
  f(1) = (dtri(1)*btrie(1) + dtri(2)*btrie(2))/denom1
  f(2) = (dtri(2)*btrie(1) - dtri(1)*btrie(2))/denom1
!
! Begin solution by traveling down (by 2's) the sub-diagonal and
! cancelling every other element
!
  if (n .ge. 3) then
     do nn = 3,n,2
        nnm2 = nn-2
        tmp(1)   = ctri (2*nn-1)*e(2*nnm2-1) - ctri(2*nn  )*e(2*nnm2  )
        tmp(2)   = ctri (2*nn-1)*e(2*nnm2  ) + ctri(2*nn  )*e(2*nnm2-1)
        denom(1) = btrie(2*nn-1) - tmp(1)
        denom(2) = btrie(2*nn  ) - tmp(2)
        tmp(1)   = ctri (2*nn-1)*f(2*nnm2-1) - ctri(2*nn  )*f(2*nnm2  )
        tmp(2)   = ctri (2*nn-1)*f(2*nnm2  ) + ctri(2*nn  )*f(2*nnm2-1)
        numer(1) = dtri (2*nn-1) + tmp(1)
        numer(2) = dtri (2*nn  ) + tmp(2)
        denom1   = denom(1)*denom(1) + denom(2)*denom(2)
        if(nn .ne. n) then
           e(2*nn-1) = (atri(2*nn-1)*denom(1) + atri(2*nn  )*denom(2))/denom1
           e(2*nn  ) = (atri(2*nn  )*denom(1) - atri(2*nn-1)*denom(2))/denom1
        endif
        f(2*nn-1) = (numer(1)*denom(1) + numer(2)*denom(2))/denom1
        f(2*nn  ) = (numer(2)*denom(1) - numer(1)*denom(2))/denom1
     end do
  endif
!
! Solve for Nth (or Nth-1) divergence element
!
  dnmn(2*n-1) = f(2*n-1)
  dnmn(2*n  ) = f(2*n  )
!
! Perform back-substitution, getting the solution for every other
! element in the divergence vector
!
  if (n .ge. 3) then
     do nn = n-2,1,-2
        nnp2 = nn+2
        tmp(1) = e(2*nn-1)*dnmn(2*nnp2-1) - e(2*nn  )*dnmn(2*nnp2  )
        tmp(2) = e(2*nn-1)*dnmn(2*nnp2  ) + e(2*nn  )*dnmn(2*nnp2-1)
        dnmn(2*nn-1) = f(2*nn-1) + tmp(1)
        dnmn(2*nn  ) = f(2*nn  ) + tmp(2)
     end do
  endif
!
  return
end subroutine trisolve
