













      module mo_lu_factor

      private
      public :: lu_fac

      contains
                                                                        
      subroutine lu_fac01( lu )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) ::   lu(:)
                                                                        
         lu(1) = 1._r8 / lu(1)
                                                                        
         lu(2) = 1._r8 / lu(2)
                                                                        
         lu(3) = 1._r8 / lu(3)
                                                                        
         lu(4) = 1._r8 / lu(4)
                                                                        
         lu(5) = 1._r8 / lu(5)
                                                                        
         lu(6) = 1._r8 / lu(6)
                                                                        
         lu(7) = 1._r8 / lu(7)
                                                                        
         lu(8) = 1._r8 / lu(8)
                                                                        
         lu(9) = 1._r8 / lu(9)
                                                                        
         lu(10) = 1._r8 / lu(10)
                                                                        
         lu(11) = 1._r8 / lu(11)
                                                                        
         lu(12) = 1._r8 / lu(12)
                                                                        
         lu(13) = 1._r8 / lu(13)
                                                                        
         lu(14) = 1._r8 / lu(14)
                                                                        
         lu(15) = 1._r8 / lu(15)
                                                                        
         lu(16) = 1._r8 / lu(16)
                                                                        
         lu(17) = 1._r8 / lu(17)
                                                                        
         lu(18) = 1._r8 / lu(18)
                                                                        
         lu(19) = 1._r8 / lu(19)
                                                                        
         lu(20) = 1._r8 / lu(20)
                                                                        
         lu(21) = 1._r8 / lu(21)
                                                                        
         lu(22) = 1._r8 / lu(22)
                                                                        
         lu(23) = 1._r8 / lu(23)
                                                                        
         lu(24) = 1._r8 / lu(24)
                                                                        
         lu(25) = 1._r8 / lu(25)
                                                                        
         lu(26) = 1._r8 / lu(26)
                                                                        
         lu(27) = 1._r8 / lu(27)
                                                                        
         lu(28) = 1._r8 / lu(28)
                                                                        
         lu(30) = 1._r8 / lu(30)
         lu(31) = lu(31) * lu(30)
         lu(35) = lu(35) - lu(31) * lu(33)
         lu(85) = lu(85) - lu(31) * lu(76)
                                                                        
         lu(34) = 1._r8 / lu(34)
         lu(35) = lu(35) * lu(34)
         lu(36) = lu(36) * lu(34)
         lu(85) = lu(85) - lu(35) * lu(77)
         lu(86) = lu(86) - lu(36) * lu(77)
                                                                        
         lu(37) = 1._r8 / lu(37)
         lu(38) = lu(38) * lu(37)
         lu(39) = lu(39) * lu(37)
         lu(85) = lu(85) - lu(38) * lu(78)
         lu(86) = lu(86) - lu(39) * lu(78)
         lu(93) = lu(93) - lu(38) * lu(88)
         lu(94) = lu(94) - lu(39) * lu(88)
                                                                        
         lu(40) = 1._r8 / lu(40)
         lu(41) = lu(41) * lu(40)
         lu(42) = lu(42) * lu(40)
         lu(43) = lu(43) * lu(40)
         lu(44) = lu(44) * lu(40)
         lu(45) = lu(45) * lu(40)
         lu(80) = lu(80) - lu(41) * lu(79)
         lu(84) = lu(84) - lu(42) * lu(79)
         lu(85) = lu(85) - lu(43) * lu(79)
         lu(86) = lu(86) - lu(44) * lu(79)
         lu(87) = lu(87) - lu(45) * lu(79)
         lu(97) = lu(97) - lu(41) * lu(96)
         lu(100) = lu(100) - lu(42) * lu(96)
         lu(101) = lu(101) - lu(43) * lu(96)
         lu(102) = lu(102) - lu(44) * lu(96)
         lu(103) = lu(103) - lu(45) * lu(96)
                                                                        
         lu(46) = 1._r8 / lu(46)
         lu(47) = lu(47) * lu(46)
         lu(48) = lu(48) * lu(46)
         lu(52) = lu(52) - lu(47) * lu(49)
         lu(53) = lu(53) - lu(48) * lu(49)
         lu(63) = lu(63) - lu(47) * lu(59)
         lu(64) = lu(64) - lu(48) * lu(59)
         lu(71) = - lu(47) * lu(66)
         lu(72) = lu(72) - lu(48) * lu(66)
         lu(85) = lu(85) - lu(47) * lu(80)
         lu(86) = lu(86) - lu(48) * lu(80)
         lu(101) = lu(101) - lu(47) * lu(97)
         lu(102) = lu(102) - lu(48) * lu(97)
                                                                        
         lu(50) = 1._r8 / lu(50)
         lu(51) = lu(51) * lu(50)
         lu(52) = lu(52) * lu(50)
         lu(53) = lu(53) * lu(50)
         lu(70) = lu(70) - lu(51) * lu(67)
         lu(71) = lu(71) - lu(52) * lu(67)
         lu(72) = lu(72) - lu(53) * lu(67)
         lu(84) = lu(84) - lu(51) * lu(81)
         lu(85) = lu(85) - lu(52) * lu(81)
         lu(86) = lu(86) - lu(53) * lu(81)
         lu(92) = lu(92) - lu(51) * lu(89)
         lu(93) = lu(93) - lu(52) * lu(89)
         lu(94) = lu(94) - lu(53) * lu(89)
                                                                        
         lu(55) = 1._r8 / lu(55)
         lu(56) = lu(56) * lu(55)
         lu(57) = lu(57) * lu(55)
         lu(58) = lu(58) * lu(55)
         lu(61) = lu(61) - lu(56) * lu(60)
         lu(63) = lu(63) - lu(57) * lu(60)
         lu(65) = lu(65) - lu(58) * lu(60)
         lu(69) = lu(69) - lu(56) * lu(68)
         lu(71) = lu(71) - lu(57) * lu(68)
         lu(73) = - lu(58) * lu(68)
         lu(83) = - lu(56) * lu(82)
         lu(85) = lu(85) - lu(57) * lu(82)
         lu(87) = lu(87) - lu(58) * lu(82)
         lu(91) = lu(91) - lu(56) * lu(90)
         lu(93) = lu(93) - lu(57) * lu(90)
         lu(95) = lu(95) - lu(58) * lu(90)
         lu(99) = lu(99) - lu(56) * lu(98)
         lu(101) = lu(101) - lu(57) * lu(98)
         lu(103) = lu(103) - lu(58) * lu(98)
                                                                        
                                                                        
      end subroutine lu_fac01
                                                                        
      subroutine lu_fac02( lu )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) ::   lu(:)
                                                                        
         lu(61) = 1._r8 / lu(61)
         lu(62) = lu(62) * lu(61)
         lu(63) = lu(63) * lu(61)
         lu(64) = lu(64) * lu(61)
         lu(65) = lu(65) * lu(61)
         lu(70) = lu(70) - lu(62) * lu(69)
         lu(71) = lu(71) - lu(63) * lu(69)
         lu(72) = lu(72) - lu(64) * lu(69)
         lu(73) = lu(73) - lu(65) * lu(69)
         lu(84) = lu(84) - lu(62) * lu(83)
         lu(85) = lu(85) - lu(63) * lu(83)
         lu(86) = lu(86) - lu(64) * lu(83)
         lu(87) = lu(87) - lu(65) * lu(83)
         lu(92) = lu(92) - lu(62) * lu(91)
         lu(93) = lu(93) - lu(63) * lu(91)
         lu(94) = lu(94) - lu(64) * lu(91)
         lu(95) = lu(95) - lu(65) * lu(91)
         lu(100) = lu(100) - lu(62) * lu(99)
         lu(101) = lu(101) - lu(63) * lu(99)
         lu(102) = lu(102) - lu(64) * lu(99)
         lu(103) = lu(103) - lu(65) * lu(99)
                                                                        
         lu(70) = 1._r8 / lu(70)
         lu(71) = lu(71) * lu(70)
         lu(72) = lu(72) * lu(70)
         lu(73) = lu(73) * lu(70)
         lu(85) = lu(85) - lu(71) * lu(84)
         lu(86) = lu(86) - lu(72) * lu(84)
         lu(87) = lu(87) - lu(73) * lu(84)
         lu(93) = lu(93) - lu(71) * lu(92)
         lu(94) = lu(94) - lu(72) * lu(92)
         lu(95) = lu(95) - lu(73) * lu(92)
         lu(101) = lu(101) - lu(71) * lu(100)
         lu(102) = lu(102) - lu(72) * lu(100)
         lu(103) = lu(103) - lu(73) * lu(100)
                                                                        
         lu(85) = 1._r8 / lu(85)
         lu(86) = lu(86) * lu(85)
         lu(87) = lu(87) * lu(85)
         lu(94) = lu(94) - lu(86) * lu(93)
         lu(95) = lu(95) - lu(87) * lu(93)
         lu(102) = lu(102) - lu(86) * lu(101)
         lu(103) = lu(103) - lu(87) * lu(101)
                                                                        
         lu(94) = 1._r8 / lu(94)
         lu(95) = lu(95) * lu(94)
         lu(103) = lu(103) - lu(95) * lu(102)
                                                                        
         lu(103) = 1._r8 / lu(103)
                                                                        
                                                                        
      end subroutine lu_fac02
                                                                        
      subroutine lu_fac( lu )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) ::   lu(:)
                                                                        
      call lu_fac01( lu )
      call lu_fac02( lu )
                                                                        
      end subroutine lu_fac
                                                                        
      end module mo_lu_factor
