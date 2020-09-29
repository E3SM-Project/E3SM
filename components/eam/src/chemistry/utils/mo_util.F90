module mo_util

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private
  public :: rebin

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


end module mo_util
