module noise_mod
  real :: thn_a, qvn_a
  integer :: nz_noise, nfr_noise
contains                                                                        

  SUBROUTINE noise(ff, nx, nz, nzn) 
    integer, intent(in) :: nx, nz, nzn
    real, intent(out) :: ff(nx, nz) 
    ! random generater seed
    integer, save :: idum=-10

    integer :: i, k
    !       do i=1,103                                                      
    !       aa=ran2(idum)                                                   
    !       enddo                                                           


    DO i = 1, nx 
       DO k = 1, nz 
          ff(i, k) = 0. 
       enddo
    enddo

    DO i = 1, nx - 1 
       DO k = 2, nzn 
          ff(i, k) = 2. * (ran2(idum) - .5) 
       enddo
    enddo
    DO k = 1, nzn 
       ff(nx, k) = ff(1, k) 
    enddo

    RETURN 
  END SUBROUTINE noise

  FUNCTION ran2(idum)
    implicit none
    double precision :: ran2
    integer, intent(inout) :: idum
    INTEGER :: im1, im2, imm1, ia1, ia2, iq1, iq2, ir1, ir2, ntab, &
         ndiv
    real ::  am, eps, rnmx
    PARAMETER(im1 = 2147483563, im2 = 2147483399, am = 1. / im1,     &
         imm1 = im1 - 1, ia1 = 40014, ia2 = 40692, iq1 = 53668, iq2 =      &
         52774, ir1 = 12211, ir2 = 3791, ntab = 32, ndiv = 1 + imm1 / ntab,&
         eps = 1.2e-7, rnmx = 1. - eps)                                    

    INTEGER idum2, j, k, iv(ntab), iy 
    SAVE iv, iy, idum2 
    DATA idum2 / 123456789 /, iv / ntab * 0 /, iy / 0 / 
    IF(idum.le.0) then 
       idum = max( - idum, 1) 
       idum2 = idum 
       DO j = ntab + 8, 1, - 1 
          k = idum / iq1 
          idum = ia1 * (idum - k * iq1) - k * ir1 
          IF (idum.lt.0) idum = idum + im1 
          IF (j.le.ntab) iv(j) = idum 
       enddo
       iy = iv(1) 
    ENDIF
    k = idum / iq1 
    idum = ia1 * (idum - k * iq1) - k * ir1 
    IF(idum.lt.0) idum = idum + im1 
    k = idum2 / iq2 
    idum2 = ia2 * (idum2 - k * iq2) - k * ir2 
    IF(idum2.lt.0) idum2 = idum2 + im2 
    j = 1 + iy / ndiv 
    iy = iv(j) - idum2 
    iv(j) = idum 
    IF(iy.lt.1) iy = iy + imm1 
    ran2 = min(am * iy, rnmx) 
    RETURN 
  END FUNCTION ran2


end module noise_mod
