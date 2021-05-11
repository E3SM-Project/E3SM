module mo_util

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private
  public :: rebin, rebin_fast

contains

  subroutine rebin( nsrc, ntrg, src_x, trg_x, src, trg )
    !---------------------------------------------------------------
    !	... rebin src to trg
    !---------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer, intent(in)   :: nsrc                  ! dimension source array
    integer, intent(in)   :: ntrg                  ! dimension target array
    real(r8), intent(in)      :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)      :: trg_x(ntrg+1)         ! target coordinates
    real(r8), intent(in)      :: src(nsrc)             ! source array
    real(r8), intent(out)     :: trg(ntrg)             ! target array

    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    integer  :: i, l
    integer  :: si, si1
    integer  :: sil, siu
    real(r8)     :: y
    real(r8)     :: sl, su
    real(r8)     :: tl, tu

    !---------------------------------------------------------------
    !	... check interval overlap
    !---------------------------------------------------------------
    !     if( trg_x(1) < src_x(1) .or. trg_x(ntrg+1) > src_x(nsrc+1) ) then
    !        write(iulog,*) 'rebin: target grid is outside source grid'
    !        write(iulog,*) '       target grid from ',trg_x(1),' to ',trg_x(ntrg+1)
    !        write(iulog,*) '       source grid from ',src_x(1),' to ',src_x(nsrc+1)
    !        call endrun
    !     end if

    do i = 1,ntrg
       tl = trg_x(i)
       if( tl < src_x(nsrc+1) ) then
          do sil = 1,nsrc+1
             if( tl <= src_x(sil) ) then
                exit
             end if
          end do
          tu = trg_x(i+1)
          do siu = 1,nsrc+1
             if( tu <= src_x(siu) ) then
                exit
             end if
          end do
          y   = 0._r8
          sil = max( sil,2 )
          siu = min( siu,nsrc+1 )
          do si = sil,siu
             si1 = si - 1
             sl  = max( tl,src_x(si1) )
             su  = min( tu,src_x(si) )
             y   = y + (su - sl)*src(si1)
          end do
          trg(i) = y/(trg_x(i+1) - trg_x(i))
       else
          trg(i) = 0._r8
       end if
    end do
  end subroutine rebin

  subroutine rebin_fast( nsrc, ntrg, ncol, src_x, trg_x, src, trg, status )
    !---------------------------------------------------------------
    !	... rebin src to trg
    !
    ! This version increments the pointers on the source and target
    ! grids so that the overall algorithm cost is
    !     O(max(nsrc,ntgt))
    ! instead of
    !     O(nsrc*ntgt)
    ! like original rebin. It also takes multiple src and trg
    ! columns so the remap bookkeeping can be amortized over
    ! multiple quantities. The implementation is intended to be BFB
    ! with rebin_orig if src_x and trg_x are each monotonically
    ! increasing.
    ! ---------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer , intent(in)  :: nsrc            ! spatial dimension source array
    integer , intent(in)  :: ntrg            ! spatial dimension target array
    integer , intent(in)  :: ncol            ! number of columns to rebin
    real(r8), intent(in)  :: src_x(nsrc+1)   ! source coordinates
    real(r8), intent(in)  :: trg_x(ntrg+1)   ! target coordinates
    real(r8), intent(in)  :: src(ncol,nsrc)  ! source array
    real(r8), intent(out) :: trg(ncol,ntrg)  ! target array
    integer , intent(out) :: status          !  0 if all is well;
                                             ! -1 if trg grid is not monotonically increasing
                                             !  1 if src grid is not monotonically increasing

    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    integer  :: i
    integer  :: si, si1
    integer  :: sil, siu
    real(r8) :: y(ncol)
    real(r8) :: sl, su
    real(r8) :: tl, tu

    sil = 1
    siu = 1
    do i = 1,ntrg
       tl = trg_x(i)
       if (tl < src_x(nsrc+1) ) then
          do while (sil < nsrc+1)
             if (tl <= src_x(sil)) exit
             sil = sil + 1
          end do
          tu = trg_x(i+1)
          if (tl >= tu) then
             status = -1
             return
          end if
          do while (siu < nsrc+1)
             if (src_x(siu) >= src_x(siu+1)) then
                status = 1
                return
             end if
             if (tu <= src_x(siu)) exit
             siu = siu + 1
          end do
          y   = 0._r8
          sil = max(sil, 2     )
          siu = min(siu, nsrc+1)
          do si = sil,siu
             si1 = si - 1
             sl = max(tl, src_x(si1))
             su = min(tu, src_x(si ))
             y  = y + (su - sl)*src(:,si1)
          end do
          trg(:,i) = y/(trg_x(i+1) - trg_x(i))
       else
          trg(:,i) = 0._r8
       end if
    end do
    status = 0
  end subroutine rebin_fast

  subroutine test_rebin()
    integer :: nzs, nzt, ncol, ne

    ne = 0
    ne = ne + test_rebin1(1, 1, 1)
    ne = ne + test_rebin1(1, 1, 11)
    do nzs = 7, 337, 41
       do nzt = 5, 331, 57
          do ncol = 1, 17, 16
             ne = ne + test_rebin1(nzs, nzt, ncol)
          end do
       end do
    end do
    ne = ne + test_rebin1(32, 32, 4, bad_src=.true.)
    ne = ne + test_rebin1(32, 32, 4, bad_tgt=.true.)
    if (ne > 0) print *, 'test_rebin FAILed: nerr', ne
  end subroutine test_rebin

  function test_rebin1(nzs, nzt, ncol, bad_src, bad_tgt) result(nerr)
    integer, intent(in) :: nzs, nzt, ncol
    logical, intent(in), optional :: bad_src, bad_tgt

    real(r8), allocatable :: zs(:), zt(:), vs(:,:), vt(:,:), vso(:,:), vto(:,:)
    real(r8) :: u, maxerr, zt1, i_src, i_tgt
    integer :: i, j, nerr, status
    logical :: bs, bt

    bs = .false.
    bt = .false.
    if (present(bad_src)) bs = bad_src
    if (present(bad_tgt)) bt = bad_tgt

    nerr = 0

    allocate(zs(nzs+1), zt(nzt+1), vso(nzs,ncol), vto(nzt,ncol), vs(ncol,nzs), vt(ncol,nzt))

    do i = 1, nzs+1
       call random_number(zs(i))
    end do
    if (bs) zs(2) = -zs(2)
    do i = 2, nzs+1
       zs(i) = zs(i) + zs(i-1)
    end do
    do i = 1, nzt+1
       call random_number(zt(i))
    end do
    if (bt) zt(2) = -zt(2)
    do i = 2, nzt+1
       zt(i) = zt(i) + zt(i-1)
    end do
    do i = 1, nzs
       do j = 1, ncol
          call random_number(vs(j,i))
          call random_number(u)
          if (u < 0.1) vs(j,i) = -vs(j,i)
          vso(i,j) = vs(j,i)
       end do
    end do

    do j = 1, ncol
       call rebin(nzs, nzt, zs, zt, vso(:,j), vto(:,j))
    end do
    call rebin_fast(nzs, nzt, ncol, zs, zt, vs, vt, status)
    ! Test that the output status from rebin_fast is the expected one.
    if (bs .and. status /=  1) nerr = nerr + 1
    if (bt .and. status /= -1) nerr = nerr + 1
    if (.not. (bs .or. bt) .and. status /= 0) nerr = nerr + 1
    ! Test that rebin_fast is BFB with rebin.
    if (.not. (bs .or. bt)) then
       maxerr = 0
       do i = 1, nzt
          do j = 1, ncol
             if (vt(j,i) /= vto(i,j)) then
                nerr = nerr + 1
                maxerr = max(maxerr, abs(vt(j,i) - vto(i,j)))
             end if
          end do
       end do
    end if
    if (.not. (bs .or. bt)) then
       ! Test the conservation of the integral.
       ! Make the grids cover same interval.
       u = (zs(nzs+1) - zs(1))/(zt(nzt+1) - zt(1))
       zt1 = zt(1)
       do i = 1, nzt+1
          zt(i) = zs(1) + u*(zt(i) - zt1)
       end do
       zt(nzt+1) = zs(nzs+1)
       if (zt(nzt) >= zt(nzt+1)) nerr = nerr + 1
       call rebin_fast(nzs, nzt, ncol, zs, zt, vs, vt, status)
       do j = 1, ncol
          i_src = 0
          do i = 1, nzs
             i_src = i_src + (zs(i+1) - zs(i))*vs(j,i)
          end do
          i_tgt = 0
          do i = 1, nzt
             i_tgt = i_tgt + (zt(i+1) - zt(i))*vt(j,i)
          end do
          if (abs(i_tgt - i_src) > 1.e3_r8*max(abs(i_tgt), abs(i_src))*epsilon(1._r8)) then
             nerr = nerr + 1
             print '(a,es10.3,es10.3,es10.3)', 'integral error:', i_src, i_tgt, i_tgt-i_src
          end if
       end do
    end if

    deallocate(zs, zt, vs, vt, vso, vto)

    if (nerr > 0) then
       print '(a,i4,i4,i4,l2,l2,i4,es9.2)', '> ', nzs, nzt, ncol, bs, bt, nerr, maxerr
    end if
  end function test_rebin1

end module mo_util
