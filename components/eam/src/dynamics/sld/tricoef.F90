subroutine tricoef(n       ,a0nm    ,bpnm    ,bmnm    ,atri    , &
                   btri    ,ctri    )
!-----------------------------------------------------------------------
!
! Purpose:
! Compute coefficients associated with solving the tri-diagonal system
! of semi-implicit diverence equations in Normal Mode space.
! NOTE 1:  Storage in the vectors assumed to be along columns ("N")
! NOTE 2:  Eigenvalue part of "btri" not added here.  To be added later
! in "TSTEP"
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
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: n         ! length of complex vector
  real(r8), intent(in)   :: a0nm(2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(in)   :: bpnm(2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(in)   :: bmnm(2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(out)  :: atri(2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(out)  :: btri(2*n) ! wave # coefs (use in vert normal mode space)
  real(r8), intent(out)  :: ctri(2*n) ! wave # coefs (use in vert normal mode space)
!
!---------------------------Local workspace-----------------------------
!
  integer  nn       ! n-wavenumber index
  integer  nnm1     ! nn - 1
  integer  nnp1     ! nn + 1
  real(r8) tmp      ! real(r8) temp space
  real(r8) tmpb1(2) ! real/imaginary temp spaces
  real(r8) tmpb2(2) ! real/imaginary temp spaces
  real(r8) denom    ! denominator
!
!-----------------------------------------------------------------------
!
! Perform "complex" arithmetic
! NOTE:  Eigenvalue part of "btri" not added here.  To be added later
! in "VERTNM"
!
  do nn = 1,n
     nnm1 = nn - 1
     nnp1 = nn + 1
     atri(2*nn-1) = 0._r8
     atri(2*nn  ) = 0._r8
     ctri(2*nn-1) = 0._r8
     ctri(2*nn  ) = 0._r8
     tmpb1(1)     = 0._r8
     tmpb1(2)     = 0._r8
     tmpb2(1)     = 0._r8
     tmpb2(2)     = 0._r8
     if(nn .ne. 1) then
        tmp          =  bmnm(2*nn  -1)*bpnm(2*nnm1-1)
        denom        =  a0nm(2*nnm1-1)*a0nm(2*nnm1-1) + a0nm(2*nnm1  )*a0nm(2*nnm1  )
        tmpb2(1)     =  a0nm(2*nnm1-1)*tmp/denom
        tmpb2(2)     = -a0nm(2*nnm1  )*tmp/denom
        tmp          =  bmnm(2*nn  -1)*bmnm(2*nnm1-1)
        ctri(2*nn-1) = -a0nm(2*nnm1-1)*tmp/denom
        ctri(2*nn  ) =  a0nm(2*nnm1  )*tmp/denom
     endif
     if(nn .ne. n) then
        tmp          =  bpnm(2*nn  -1)*bmnm(2*nnp1-1)
        denom        =  a0nm(2*nnp1-1)*a0nm(2*nnp1-1) + a0nm(2*nnp1  )*a0nm(2*nnp1  )
        tmpb1(1)     =  a0nm(2*nnp1-1)*tmp/denom
        tmpb1(2)     = -a0nm(2*nnp1  )*tmp/denom
        tmp          =  bpnm(2*nn  -1)*bpnm(2*nnp1-1)
        atri(2*nn-1) = -a0nm(2*nnp1-1)*tmp/denom
        atri(2*nn  ) =  a0nm(2*nnp1  )*tmp/denom
     endif
!
     btri(2*nn-1) = a0nm(2*nn-1) + tmpb1(1) + tmpb2(1)
     btri(2*nn  ) = a0nm(2*nn  ) + tmpb1(2) + tmpb2(2)
  end do
!
  return
end subroutine tricoef
