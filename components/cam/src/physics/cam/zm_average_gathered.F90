module zm_average_gathered

  use shr_kind_mod, only: r8 => shr_kind_r8
#ifndef TEST_MODE
  use ppgrid, only: pcols, pver, pverp
  use shr_assert_mod, only: shr_assert, shr_assert_in_domain
#endif

  implicit none
  private

  public zm_gath_avg
  public zm_gath_avg_init
  public zm_gath_avg_accum
  public zm_gath_avg_output

#ifdef TEST_MODE
  ! For testing
  integer, parameter, public :: pcols = 8
  integer, parameter, public :: pver = 8
  integer, parameter, public :: pverp = pver + 1
#endif

  type zm_gath_avg
     ! How often has each column had convection?
     integer :: con_count(pcols)
     ! How many substeps have occurred?
     integer :: step_count
     ! Un-gathered fields to average and gather
     real(r8) :: mu(pcols, pver)
     real(r8) :: eu(pcols, pver)
     real(r8) :: du(pcols, pver)
     real(r8) :: md(pcols, pver)
     real(r8) :: ed(pcols, pver)
     ! layer thickness
     real(r8) :: dp(pcols, pver)
     ! top and bottom of convective region
     integer :: jt(pcols)
     integer :: maxg(pcols)
  end type zm_gath_avg

contains

  subroutine zm_gath_avg_init(pdel, zga)
    real(r8), intent(in) :: pdel(pcols, pver)
    type(zm_gath_avg), intent(out) :: zga

    zga%con_count = 0
    zga%step_count = 0

    zga%mu = 0.
    zga%eu = 0.
    zga%du = 0.
    zga%md = 0.
    zga%ed = 0.

    zga%dp = 0.01 * pdel

    zga%jt = pver
    zga%maxg = 1
  end subroutine zm_gath_avg_init

  subroutine zm_gath_avg_accum( &
       lengath, ideep, mu, eu, du, &
       md, ed, dp, dsubcld, jt, &
       maxg, zga)
    integer, intent(in) :: lengath
    integer, intent(in) :: ideep(pcols)
    real(r8), intent(in) :: mu(pcols, pver)
    real(r8), intent(in) :: eu(pcols, pver)
    real(r8), intent(in) :: du(pcols, pver)
    real(r8), intent(in) :: md(pcols, pver)
    real(r8), intent(in) :: ed(pcols, pver)
    real(r8), intent(in) :: dp(pcols, pver)
    real(r8), intent(in) :: dsubcld(pcols)
    integer, intent(in) :: jt(pcols)
    integer, intent(in) :: maxg(pcols)

    type(zm_gath_avg), intent(inout) :: zga

    integer :: i, k

    zga%step_count = zga%step_count + 1

    do i = 1, lengath
       zga%con_count(ideep(i)) = zga%con_count(ideep(i)) + 1

       zga%mu(ideep(i),:) = zga%mu(ideep(i),:) + mu(i,:)
       zga%eu(ideep(i),:) = zga%eu(ideep(i),:) + eu(i,:)
       zga%du(ideep(i),:) = zga%du(ideep(i),:) + du(i,:)
       zga%md(ideep(i),:) = zga%md(ideep(i),:) + md(i,:)
       zga%ed(ideep(i),:) = zga%ed(ideep(i),:) + ed(i,:)

#ifdef TEST_MODE
       if (any(zga%dp(ideep(i),:) /= dp(i,:))) then
          stop 1
       end if
       if (dsubcld(i) /= sum(zga%dp(ideep(i),maxg(i):pver))) then
          stop 2
       end if
#else
       ! Just to test that I have correctly traced this during development.
       ! Do not put this in production.
       do k = 1, pver
          call shr_assert_in_domain(abs(zga%dp(ideep(i), k) - dp(i,k)), &
               lt=1.e-2_r8, varname="dp error")
       end do
       call shr_assert_in_domain(abs(dsubcld(i) - sum(zga%dp(ideep(i),maxg(i):pver))), &
            lt=1.e-2_r8, varname="dsubcld error")
#endif

       zga%jt(ideep(i)) = min(zga%jt(ideep(i)), jt(i))
       zga%maxg(ideep(i)) = max(zga%maxg(ideep(i)), maxg(i))
    end do

  end subroutine zm_gath_avg_accum

  subroutine zm_gath_avg_output( &
       zga, lengath, ideep, mu, eu, &
       du, md, ed, dp, dsubcld, &
       jt, maxg)
    type(zm_gath_avg), intent(in) :: zga

    integer, intent(out) :: lengath
    integer, intent(out) :: ideep(pcols)
    real(r8), intent(out) :: mu(pcols, pver)
    real(r8), intent(out) :: eu(pcols, pver)
    real(r8), intent(out) :: du(pcols, pver)
    real(r8), intent(out) :: md(pcols, pver)
    real(r8), intent(out) :: ed(pcols, pver)
    real(r8), intent(out) :: dp(pcols, pver)
    real(r8), intent(out) :: dsubcld(pcols)
    integer, intent(out) :: jt(pcols)
    integer, intent(out) :: maxg(pcols)

    integer :: i

    lengath = 0
    do i = 1, pcols
       if (zga%con_count(i) > 0) then
          lengath = lengath + 1
          ideep(lengath) = i
       end if
    end do

    do i = 1, lengath
       mu(i,:) = zga%mu(ideep(i),:) / zga%step_count
       eu(i,:) = zga%eu(ideep(i),:) / zga%step_count
       du(i,:) = zga%du(ideep(i),:) / zga%step_count
       md(i,:) = zga%md(ideep(i),:) / zga%step_count
       ed(i,:) = zga%ed(ideep(i),:) / zga%step_count

       jt(i) = zga%jt(ideep(i))
       maxg(i) = zga%maxg(ideep(i))

       dp(i,:) = zga%dp(ideep(i),:)
       dsubcld(i) = sum(dp(i,maxg(i):pver))
    end do

  end subroutine zm_gath_avg_output

end module zm_average_gathered
