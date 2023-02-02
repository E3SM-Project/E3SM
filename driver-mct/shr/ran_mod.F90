module ran_mod
! module contains three functions
! ran1 returns a uniform random number between 0-1
! spread returns random number between min - max
! normal returns a normal distribution
use shr_kind_mod, only : r8=> shr_kind_r8
!use mod_kinds        , only : r8
public :: ran1, spread, normal,init_random_seed,randperm


contains
    function ran1()  !returns random number between 0 - 1
!        use numz
        implicit none
        real(r8):: ran1,x
!        call init_random_seed()
        call random_number(x) ! built in fortran 90 random number function
        ran1=x
    end function ran1

    function spread(min,max)  !returns random number between min - max
!        use numz
        implicit none
        real(r8):: spread
        real(r8):: min,max
        spread=(max - min) * ran1() + min
    end function spread

    function normal() !returns a normal distribution of (0,1.0)
 !       use numz
        implicit none
!        real(r8),intent(in):: mean,sigma
        real(r8):: normal,tmp
        integer :: flag
        real(r8):: fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0
            do while(rsq.ge.1.0.or.rsq.eq.0.0) ! new from for do
                r1=2.0*ran1()-1.0
                r2=2.0*ran1()-1.0
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp
        return
    end function normal


   SUBROUTINE init_random_seed()
     implicit none
     INTEGER :: i, n, clock
     INTEGER, DIMENSION(:), ALLOCATABLE :: seed
     CALL RANDOM_SEED(size = n)
     ALLOCATE(seed(n))
     CALL SYSTEM_CLOCK(COUNT=clock)
     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     CALL RANDOM_SEED(PUT = seed)
     DEALLOCATE(seed)
     END SUBROUTINE init_random_seed



function randperm(num)
!  use data_type, only : IB, RP
  implicit none
  integer(kind=4), intent(in) :: num
  integer(kind=4) :: number, i, j, k
  integer(kind=4), dimension(num) :: randperm
  real(r8), dimension(num) :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
               number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
               number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
end function randperm


end module ran_mod
