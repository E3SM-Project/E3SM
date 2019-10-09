
      module mo_wavelen

      use mo_params,    only : kw
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

      integer     :: nw
      real(r8)    :: wl(kw), wc(kw), wu(kw)
      real(r8)    :: deltaw(kw)
      real(r8)    :: delw_bin(kw)
      real(r8)    :: sflx(kw)

      end module mo_wavelen
