subroutine nmmatrix(zb      ,zcr1    ,bm1     ,bmi     )
!-----------------------------------------------------------------------
!
! Purpose:
! Determine eigenvalues and left/right eigenvectors of the hydrostatic
! matrix.  Results used to solve for two-time-level SLD, semi-implicit
! Coriolis term
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
  use sgexx,        only: sgeev
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog

  implicit none

!------------------------------Parameters-------------------------------
!
  integer, parameter :: lwork = 100*plev    ! work array dimension
!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: zb  (plev,plev) ! semi-implicit matrix in d equation
  real(r8), intent(out)  :: zcr1(plev)      ! eigenvalues (real)
  real(r8), intent(out)  :: bm1 (plev,plev) ! transpose of right eigenvector (normalized)
  real(r8), intent(out)  :: bmi (plev,plev) ! transpose of inverse of bm1
!                                           ! ( = left eigenvectors --  normalized)
!
!---------------------------Local variables-----------------------------
!
  integer  k                    ! index
  integer  kk                   ! index
  integer  kkk                  ! index
  integer  istat                ! status of SGEEV call
  real(r8) eps                  ! used for epsilon tests
  real(r8) sum                  ! tmp
  real(r8) fmax                 ! tmp

  real(r8) x      (plev,plev)   ! local copy of "zb"
  real(r8) bmlft  (plev,plev)   ! left  eigenvector
  real(r8) bmrgt  (plev,plev)   ! right eigenvector
  real(r8) zci    (plev)        ! eigenvalues (complex)
  real(r8) work   (lwork)       ! work array
  real(r8) diag                 ! tmp

  logical  lcomp                ! error flag
  logical  lneg                 ! error flag
  logical  lsame                ! error flag
  logical  lnonzero             ! error flag
  logical  flag(plev)           ! flag used in sorting eigenvalues/vectors
!
!-----------------------------------------------------------------------
!
! Matrix "zb" is stored 1st index = column; 2nd index = row of matrix.
! Transpose "zb" and store in work array before calling library routine
!
  do k = 1,plev
     do kk = 1,plev
        x(k,kk) = zb(kk,k)
     end do
  end do
!
! Determine eigenvalues and left/right eigenvectors of "zb"
!
  call sgeev('v'     ,'v'     ,plev    ,x       ,plev    , &
             zcr1    ,zci     ,bmlft   ,plev    ,bmrgt   , &
             plev    ,work    ,lwork   ,istat   )
!
  if (istat .ne. 0) then
     write(iulog,*) 'ERROR in NMMATRIX:  SGEEV returned "istat" = ',istat
     call endrun ()
  endif
  if (int(work(1)) .gt. lwork) then
     write(iulog,*) 'ERROR in NMMATRIX: dimension of work array must be at least', int(work(1))
     call endrun ()
  endif
!
! Find "eps" based upon the real eigenvalue with the largest magnitude
!
  fmax = -1.e36_r8
  do kk = 1,plev
     if(zcr1(kk) .gt. fmax) fmax = zcr1(kk)
  end do
  eps = fmax*epsilon(fmax)*10._r8
!
! Test that all eigenvalues are positive, real, and distinct
!
  lcomp    = .false.
  lneg     = .false.
  lsame    = .false.
  do k = 1,plev
     if(    zcr1(k)  .lt.  0._r8) lneg  = .true.
     if(abs(zci(k  )) .gt. eps) lcomp = .true.
     do kk = k+1,plev
        if(abs( (zcr1(k) - zcr1(kk)) ) .lt. eps) lsame = .true.
     end do
  end do
  if (lneg) then
     call endrun ('ERROR in NMMATRIX: negative eigenvalue(s)')
  endif
  if (lcomp) then
     call endrun ('ERROR in NMMATRIX: complex eigenvalue(s)')
  endif
  if (lsame) then
     call endrun ('ERROR in NMMATRIX: non-distinct eigenvalue(s)')
  endif
!
! Re-order the eigenvalues (and associated eigenvectors) in descending
! order.  Use "zci", "bmi", and "bm1" as temporary arrays.
!
  do k = 1,plev
     flag(k) = .true.
  end do
  do k = 1,plev
     fmax = -1.e35_r8
     do kk = 1,plev
        if(flag(kk) .and. zcr1(kk) .gt. fmax) then
           kkk = kk
           fmax = zcr1(kk)
        endif
     end do
     flag(kkk) = .false.
     zci(k) = zcr1(kkk)
     do kk = 1,plev
        bmi(kk,k) = bmlft(kk,kkk)
        bm1(kk,k) = bmrgt(kk,kkk)
     end do
  end do
!
! Copy arrays from temporaries back to original arrays
!
  do k = 1,plev
     zcr1(k) = zci(k)
     do kk = 1,plev
        bmlft(kk,k) = bmi(kk,k)
        bmrgt(kk,k) = bm1(kk,k)
     end do
  end do
!
! Make sure the element of largest magnitude within each eigenvector
! is positive (in this case, we will choose the LEFT eigenvector).
! Adjust the signs of all other elements of the left AND associated
! RIGHT eigenvector, accordingly.
!
! NOTE:  This is done only to ensure that the signs of all eigenvectors
! will be consistent across all computing platforms.  This step is
! not necessary except that it makes debugging the "vertical normal
! mode" code MUCH easier.
!
  do k = 1,plev
     fmax = -1.e35_r8
     do kk = 1,plev
        if(abs(bmlft(kk,k)) .gt. fmax) then
           kkk = kk
           fmax = abs(bmlft(kk,k))
        endif
     end do
     if(bmlft(kkk,k) .lt. 0._r8) then
        do kk = 1,plev
           bmlft(kk,k) = -bmlft(kk,k)
           bmrgt(kk,k) = -bmrgt(kk,k)
        end do
     endif
  end do
!
! Normalize "bm1" and determine its inverse (and confirm).
!
! NOTE:  both "bm1" and "bmi" are stored in transpose form (1st index =
! column; 2nd index = row of matrix) for later use in coding vector
! inner products using these matrices.
!
  lnonzero = .false.
  do k = 1,plev
     sum = 0._r8
     do kk = 1,plev
        sum = sum + bmlft(kk,k)*bmrgt(kk,k)
     end do
!
     if(sum .lt. 0._r8) then
        sum = sqrt(-sum)
        do kk = 1,plev
           bmi(kk,k) = -bmlft(kk,k)/sum
           bm1(k,kk) =  bmrgt(kk,k)/sum
        end do
     else
        sum = sqrt(sum)
        do kk = 1,plev
           bmi(kk,k) =  bmlft(kk,k)/sum
           bm1(k,kk) =  bmrgt(kk,k)/sum
        end do
     endif
  end do
  do k = 1,plev
     do kk = 1,plev
        diag = 0._r8
        do kkk = 1,plev
           diag = diag + bmi(kkk,kk)*bm1(k,kkk)
        end do
        if(kk .ne. k) then
           if(diag .gt. eps) lnonzero = .true.
        endif
     end do
  end do
  if (lnonzero) then
     call endrun ('NMMATRIX: Y(T)X has non-zero off-diagonal elements')
  endif
!
  return
end subroutine nmmatrix
