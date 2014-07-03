#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module domain_mod
  use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
  implicit none
  private

  type, public :: domain1D_t
     integer :: start                     ! domain start index
     integer :: end                       ! domain end index
  end type domain1D_t

  type, public :: domain2D_t
     type (domain1D_t) :: x
     type (domain1D_t) :: y
  end type domain2D_t

! ==========================================
! Public Interfaces
! ==========================================
  
  public :: decompose
contains           
  
  function decompose(start,end,ndomains,ipe) result(domain)

    integer, intent(in) :: start      ! starting index
    integer, intent(in) :: end        ! ending   index
    integer, intent(in) :: ndomains   ! number of domains
    integer, intent(in) :: ipe        ! domain number

    type (domain1D_t) :: domain

    ! Local stuff 

    integer :: beg(0:ndomains)
    integer :: length
    integer :: n

    call t_startf('decompose')
    length = end-start+1
    beg(0)=start
    do n=1,ndomains-1
       if (n.le.mod(length,ndomains)) then
          beg(n)=beg(n-1)+(length-1)/ndomains+1
       else
          beg(n)=beg(n-1)+length/ndomains
       end if
    end do
    beg(ndomains)=start+length

    domain%start=beg(ipe)
    domain%end  =beg(ipe+1)-1

    call t_stopf('decompose')
  end function decompose


end module domain_mod





