program endian_convert
!
!  convert a HOMME vertical coordinate file from big endian to little endia
!  or vice-versa
!
! for G95:  
!   setenv G95_UNIT_ENDIAN_11 big
!   setenv G95_UNIT_ENDIAN_12 little
!   g95 endian_convert.F90 ; ./a.out
!
      real*8 :: AK,BK
      integer k
      integer nlev
      character(len=80) :: fname

      nlev=26

      fname="cami-26.fbin"
      open(UNIT=11,FILE=fname,form='unformatted',status='unknown',iostat=ios)
      fname="cami-26.fbin.converted"
      open(UNIT=12,FILE=fname,form='unformatted',status='unknown',iostat=ios)

      do k=1,nlev+1
         read(11) AK
         read(11) BK
         print *,k,AK,BK
         write(12) AK
         write(12) BK
      enddo
      close(11)
      close(12)

      fname="camm-26.fbin"
      open(UNIT=11,FILE=fname,form='unformatted',status='unknown',iostat=ios)
      fname="camm-26.fbin.converted"
      open(UNIT=12,FILE=fname,form='unformatted',status='unknown',iostat=ios)

      do k=1,nlev
         read(11) AK
         read(11) BK
         print *,k,AK,BK
         write(12) AK
         write(12) BK
      enddo
      close(11)
      close(12)

end program

