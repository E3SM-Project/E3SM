program zm_average_gathered_test

  use shr_kind_mod, only: r8 => shr_kind_r8

  use zm_average_gathered, only: &
       zm_gath_avg, &
       zm_gath_avg_init, &
       zm_gath_avg_accum, &
       zm_gath_avg_output, &
       pcols, &
       pver, &
       pverp

  implicit none

  real(r8) :: pdel(pcols, pver)

  real(r8) :: mu(pcols, pver)
  real(r8) :: eu(pcols, pver)
  real(r8) :: du(pcols, pver)
  real(r8) :: md(pcols, pver)
  real(r8) :: ed(pcols, pver)
  real(r8) :: dp(pcols, pver)
  real(r8) :: dsubcld(pcols)
  integer :: jt(pcols)
  integer :: maxg(pcols)
  integer :: ideep(pcols)
  integer :: lengath

  type(zm_gath_avg) :: zga
  
  real(r8) :: con_count(pcols, pver)

  integer :: i, k

  do k = 1, pver
     do i = 1, pcols
        pdel(i,k) = k + i*0.5
     end do
  end do

  call zm_gath_avg_init(pdel, zga)

  lengath = 3
  ideep(1) = 2
  ideep(2) = 4
  ideep(3) = 5

  jt(1) = 3
  maxg(1) = 6

  mu(1, :) = 2.
  eu(1, :) = 2. - 1.
  du(1, :) = 2. - 2.
  md(1, :) = 2. - 3.
  ed(1, :) = 2. - 4.

  dp(1, :) = 0.01 * pdel(2, :)

  jt(2) = 3
  maxg(2) = 6

  mu(2, :) = 1.
  eu(2, :) = 1. - 1.
  du(2, :) = 1. - 2.
  md(2, :) = 1. - 3.
  ed(2, :) = 1. - 4.

  dp(2, :) = 0.01 * pdel(4, :)

  jt(3) = 3
  maxg(3) = 6

  mu(3, :) = 3.
  eu(3, :) = 3. - 1.
  du(3, :) = 3. - 2.
  md(3, :) = 3. - 3.
  ed(3, :) = 3. - 4.

  dp(3, :) = 0.01 * pdel(5, :)

  dsubcld = 0.
  do k = 1, pver
     do i = 1, lengath
        if (k >= maxg(i)) then
           dsubcld(i) = dsubcld(i) + dp(i,k)
        end if
     end do
  end do

  call zm_gath_avg_accum( &
       lengath, ideep, mu, eu, du, &
       md, ed, dp, dsubcld, jt, &
       maxg, zga)

  lengath = 3
  ideep(1) = 4
  ideep(2) = 5
  ideep(3) = 7

  jt(1) = 4
  maxg(1) = 5

  mu(1, :) = 2.
  eu(1, :) = 2. - 1.
  du(1, :) = 2. - 2.
  md(1, :) = 2. - 3.
  ed(1, :) = 2. - 4.

  dp(1, :) = 0.01 * pdel(4, :)

  jt(2) = 2
  maxg(2) = 7

  mu(2, :) = 1.
  eu(2, :) = 1. - 1.
  du(2, :) = 1. - 2.
  md(2, :) = 1. - 3.
  ed(2, :) = 1. - 4.

  dp(2, :) = 0.01 * pdel(5, :)

  jt(3) = 3
  maxg(3) = 6

  mu(3, :) = 4.
  eu(3, :) = 4. - 1.
  du(3, :) = 4. - 2.
  md(3, :) = 4. - 3.
  ed(3, :) = 4. - 4.

  dp(3, :) = 0.01 * pdel(7, :)

  dsubcld = 0.
  do k = 1, pver
     do i = 1, lengath
        if (k >= maxg(i)) then
           dsubcld(i) = dsubcld(i) + dp(i,k)
        end if
     end do
  end do

  call zm_gath_avg_accum( &
       lengath, ideep, mu, eu, du, &
       md, ed, dp, dsubcld, jt, &
       maxg, zga)

  dp = 0.

  call zm_gath_avg_output( &
       zga, lengath, ideep, mu, eu, &
       du, md, ed, dp, dsubcld, &
       jt, maxg)

  if (lengath /= 4) then
     print *, "Wrong lengath: ", lengath
  end if

  if (any(ideep(1:4) /= [2, 4, 5, 7])) then
     print *, "Wrong ideep: ", ideep
  end if

  if (any(jt(1:4) /= [3, 3, 2, 3])) then
     print *, "Wrong jt: ", jt
  end if

  if (any(maxg(1:4) /= [6, 6, 7, 6])) then
     print *, "Wrong maxg: ", maxg
  end if

  if (any(mu(1, :) /= 1.)) then
     print *, "Wrong mu1: ", mu(1, :)
  end if

  if (any(mu(2, :) /= 1.5)) then
     print *, "Wrong mu2: ", mu(2, :)
  end if

  if (any(mu(3, :) /= 2.)) then
     print *, "Wrong mu3: ", mu(3, :)
  end if

  if (any(mu(4, :) /= 2.)) then
     print *, "Wrong mu4: ", mu(4, :)
  end if

  con_count = 0.
  con_count(1,:) = 0.5
  con_count(2,:) = 1
  con_count(3,:) = 1.
  con_count(4,:) = 0.5

  if (any(mu(1:4,:) /= eu(1:4,:) + 1.*con_count(1:4,:))) then
     print *, "Wrong eu: ", mu(1:4,:), eu(1:4,:) + 1.*con_count(1:4,:)
  end if

  if (any(mu(1:4,:) /= du(1:4,:) + 2.*con_count(1:4,:))) then
     print *, "Wrong du: ", mu(1:4,:), du(1:4,:)
  end if

  if (any(mu(1:4,:) /= md(1:4,:) + 3.*con_count(1:4,:))) then
     print *, "Wrong md: ", mu(1:4,:), md(1:4,:)
  end if

  if (any(mu(1:4,:) /= ed(1:4,:) + 4.*con_count(1:4,:))) then
     print *, "Wrong ed: ", mu(1:4,:), ed(1:4,:)
  end if

  if (any(dp(1,:) /= 0.01 * pdel(2,:))) then
     print *, "Wrong dp1: ", dp(1,:)
  end if

  if (any(dp(2,:) /= 0.01 * pdel(4,:))) then
     print *, "Wrong dp2: ", dp(2,:)
  end if

  if (any(dp(3,:) /= 0.01 * pdel(5,:))) then
     print *, "Wrong dp3: ", dp(3,:)
  end if

  if (any(dp(4,:) /= 0.01 * pdel(7,:))) then
     print *, "Wrong dp4: ", dp(4,:)
  end if

  if (dsubcld(1) /= 0.01 * pdel(2, 6) + 0.01 * pdel(2, 7) + 0.01 * pdel(2, 8) ) then
     print *, "Wrong dsubcld1: ", dsubcld(1)
  end if

  if (dsubcld(2) /= 0.01 * pdel(4, 6) + 0.01 * pdel(4, 7) + 0.01 * pdel(4, 8) ) then
     print *, "Wrong dsubcld2: ", dsubcld(2)
  end if

  if (dsubcld(3) /= 0.01 * pdel(5, 7) + 0.01 * pdel(5, 8) ) then
     print *, "Wrong dsubcld3: ", dsubcld(3)
  end if

  if (dsubcld(4) /= 0.01 * pdel(7, 6) + 0.01 * pdel(7, 7) + 0.01 * pdel(7, 8) ) then
     print *, "Wrong dsubcld4: ", dsubcld(4)
  end if

end program zm_average_gathered_test
