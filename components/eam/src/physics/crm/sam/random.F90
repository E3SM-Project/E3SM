module random_mod
  implicit none

contains

  ! Simple randaom number generator in the range [0,1]
  ! ranset_(iseed) initializes with iseed
  ! ranf_() returns next random numer

  real(crm_rknd) function ranf_()
    use params, only: crm_rknd
    implicit none
    ! ranf_ = rand_(0)
    call random_number(ranf_)
    return
  end

  subroutine ranset_(iseed)
    use params, only: crm_rknd
    implicit none
    integer iseed, i, m, nsteps
    ! i = rand_(1) ! reinitialize (reset)
    nsteps = iseed*10000
    do i = 1,nsteps
      m = ranf_()
      !	m = rand_(0)
    end do
    return
  end

  real(crm_rknd) function rand_(iseed)
    use params, only: crm_rknd
    implicit none
    integer iseed
    integer ia1, ia0, ia1ma0, ic, ix1, ix0, iy0, iy1
    save ia1, ia0, ia1ma0, ic, ix1, ix0
    data ix1, ix0, ia1, ia0, ia1ma0, ic/0,0,1536,1029,507,1731/
    if (iseed.ne.0) then
      ia1 = 1536
      ia0 = 1029
      ia1ma0 = 507
      ic = 1731
      ix1 = 0
      ix0 = 0
      rand_ = 0
    else
      iy0 = ia0*ix0
      iy1 = ia1*ix1 + ia1ma0*(ix0-ix1) + iy0
      iy0 = iy0 + ic
      ix0 = mod (iy0, 2048)
      iy1 = iy1 + (iy0-ix0)/2048
      ix1 = mod (iy1, 2048)
      rand_ = ix1*2048 + ix0
      rand_ = rand_ / 4194304.D0
    end if
    return
  end

end module random_mod
