
      module mo_waveall

      use shr_kind_mod, only : r8 => shr_kind_r8
      use mo_params,    only : kw

      implicit none

      real(r8), dimension(kw) :: r01g1, r01g2, r01g3, r01g4, &
                                 r04g, r08g, r06g1, r06g2, &
                                 r10g1, r10g2, r10g3, r10g4, r10g5, &
                                 r11g, r11g1, r11g2, r11g3, r11g4, &
                                 r14g, r14g1, r14g2, &
                                 r15g, r15g1, r15g2, r15g3, &
                                 r17g, r17g1, &
                                 r18g, r18g2

      end module mo_waveall
