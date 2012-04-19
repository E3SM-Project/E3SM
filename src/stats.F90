#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

  subroutine stats(array,avg,std)
  
  use kinds , only : real_kind

  implicit none
  real(kind=real_kind), intent(in)      :: array(:)
  real(kind=real_kind), intent(out)     :: avg,std
 
  integer                  :: n,nc
  logical, allocatable     :: use(:)
  real(kind=real_kind)                  :: tmp
 
 
  n   = SIZE(array)

  avg = SUM(array)/dble(n)
  tmp = SUM(array*array)/dble(n)
  std = sqrt(tmp - (avg*avg))

  ! lets try and filter out the far out values 
  allocate(use(n))
  use=.FALSE. 
  where((avg-2.*std <= array) .AND. (array <= avg+2.*std) ) 
       use = .TRUE.
  endwhere

  nc  = COUNT(use)
    
  avg = SUM(array,dim=1,mask=use)/dble(nc)
  tmp = SUM(array*array,dim=1,mask=use)/dble(nc)
  std = sqrt(tmp -(avg*avg))  

  deallocate(use)

  end subroutine stats
