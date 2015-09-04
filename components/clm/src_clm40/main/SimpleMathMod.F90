module SimpleMathMod

#include "shr_assert.h"
  !------------------------------------------------------------------------------
  !
  ! DESCRIPTIONS:
  ! module contains simple mathematical functions for arrays
  ! Created by Jinyun Tang, Feb., 2014

implicit none
  
  interface array_normalization
    module procedure array_normalization_2d, array_normalization_2d_filter
  end interface array_normalization
  
  interface array_div_vector
    module procedure array_div_vector_filter, array_div_vector_nofilter
  end interface array_div_vector
contains
!--------------------------------------------------------------------------------
  subroutine array_normalization_2d(which_dim, arr2d_inout)
  !
  !DESCRIPTIONS
  !do normalization for the input array along dimension which_dim
  !
  !USES
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  
  integer,  intent(in) :: which_dim     !do normalization along which dimension?
  real(r8), intent(inout) :: arr2d_inout(:,:)   !input 2d array    

  
  !local variables
  integer  :: sz1, sz2     !array size
  integer  :: j1, j2       !indices
  real(r8) :: arr_sum
  
  sz1 = size(arr2d_inout,1)
  sz2 = size(arr2d_inout,2)
  
  if(which_dim==1)then
    !normalize along dimension 1, so loop along dimension 2     
    do j2 = 1, sz2
      !obtain the total
      arr_sum=0._r8
      do j1 = 1, sz1
        arr_sum=arr_sum+arr2d_inout(j1,j2)
      enddo
      !normalize with the total if arr_sum is non-zero
      if(arr_sum/=0._r8)then
        do j1 = 1, sz1
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr_sum
        enddo
      endif
    enddo
  elseif(which_dim==2)then
    !normalize along dimension 2, so loop along dimension 1
    do j1 = 1, sz1
      !obtain the total
      arr_sum=0._r8
      do j2 = 1, sz2
        arr_sum=arr_sum+arr2d_inout(j1,j2)
      enddo
      !normalize with the total if arr_sum is non-zero
      !I think there should be a safer mask for this to screen off spval values
      !Jinyun Tang, May 30, 2014
      if(arr_sum>0._r8 .or. arr_sum < 0._r8)then
        do j2 = 1, sz2
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr_sum
        enddo
      endif
    enddo
  endif
  return
  end subroutine array_normalization_2d
  
!--------------------------------------------------------------------------------
  subroutine array_normalization_2d_filter(lbj1, ubj1, lbj2, ubj2, numf, filter, arr2d_inout)
  !
  !DESCRIPTIONS
  !do normalization with filter for the input array along dimension 2
  
  !
  !USES
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg  
  implicit none
  integer,  intent(in) :: lbj1         !left bound of dim 1
  integer,  intent(in) :: lbj2         !left bound of dim 2
  integer,  intent(in) :: ubj1         !right bound of dim 1
  integer,  intent(in) :: ubj2         !right bound of dim 2
  integer,  intent(in) :: numf         !filter size
  integer,  intent(in) :: filter(:)    !filter
  real(r8), intent(inout) :: arr2d_inout(lbj1: , lbj2: )   !input 2d array
  
  
  !local variables
  integer  :: sz1, sz2     !array size
  integer  :: j2           !indices
  integer  :: f, p         !indices
  real(r8) :: arr_sum(lbj1:ubj1)

  ! Enforce expected array sizes
  SHR_ASSERT_ALL((ubound(arr2d_inout) == (/ubj1, ubj2/)),      errMsg(__FILE__, __LINE__))
  

  arr_sum(:) = 0._r8  
  do j2 = lbj2, ubj2  
    do f = 1, numf
      p = filter(f)
      !obtain the total      
      arr_sum(p)=arr_sum(p)+arr2d_inout(p,j2)
    enddo
  enddo
  
    !normalize with the total if arr_sum is non-zero
  do j2 = lbj2, ubj2
    do f = 1, numf
      p = filter(f)
      !I found I have to ensure >0._r8 because of some unknown reason, jyt May 23, 2014
      !I will test this later with arr_sum(p)/=0._r8
      if(arr_sum(p)>0._r8 .or. arr_sum(p)<0._r8)then
        arr2d_inout(p,j2) = arr2d_inout(p,j2)/arr_sum(p)
      endif
    enddo
  enddo 
  return
  end subroutine array_normalization_2d_filter
!--------------------------------------------------------------------------------  
  
  subroutine array_div_vector_filter(lbj1, ubj1, lbj2, ubj2, &
       arr1d_in, fn, filter,  arr2d_inout)
  !
  !DESCRIPTIONS
  !array divided by a vector, arr2d_in is divided by one
  !element in arr1d_in  
  !It always assumes the filter is along with dimenion 1
  !
  ! USES
  !
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg   
  implicit none
  integer,  intent(in) :: lbj1         !left bound of dim 1
  integer,  intent(in) :: lbj2         !left bound of dim 2
  integer,  intent(in) :: ubj1         !right bound of dim 1
  integer,  intent(in) :: ubj2         !right bound of dim 2 
  real(r8), intent(in) :: arr1d_in(lbj1: )   !1d scaling factor
  integer , intent(in) :: fn
  integer , intent(in) :: filter(:)     !filter
  real(r8), intent(inout) :: arr2d_inout(lbj1: ,lbj2: ) !2d array to be scaled  

  integer :: sz
  integer :: j, f, p
  
  ! Enforce expected array sizes
  SHR_ASSERT_ALL((ubound(arr2d_inout) == (/ubj1, ubj2/)),      errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(arr1d_in) == (/ubj1/)),            errMsg(__FILE__, __LINE__))


  do j = lbj2, ubj2
     do f = 1, fn
        p = filter(f)
        if (arr1d_in(p) > 0._r8 .or. arr1d_in(p) < 0._r8) then
           arr2d_inout(p,j) = arr2d_inout(p,j)/arr1d_in(p)
        else
           arr2d_inout(p,j) = 0._r8
        end if
     end do
  end do
  return
  end subroutine array_div_vector_filter
  
!--------------------------------------------------------------------------------  
  
  subroutine array_div_vector_nofilter(arr1d_in, which_dim, arr2d_inout)
  !
  !DESCRIPTIONS
  !array divided by a vector, each row in arr2d_in is divided by one
  !element in arr1d_in
  !
  !USES
  !
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_assert_mod , only : shr_assert
  use shr_log_mod    , only : errMsg => shr_log_errMsg  
  implicit none
  real(r8), intent(in) :: arr1d_in(:)     !scaling factor
  integer,  intent(in) :: which_dim        !which dimension is scaled
  real(r8), intent(inout) :: arr2d_inout(:,:)   !2d array to be scaled
  
  integer :: sz1, sz2
  integer :: j1, j2

  sz1=size(arr2d_inout,1)  
  sz2=size(arr2d_inout,2)
  
  if(which_dim==1)then
    ! Enforce expected array sizes   
    call shr_assert(sz1    == size(arr1d_in), errMsg(__FILE__, __LINE__))
    
    do j2 = 1, sz2
      do j1 = 1, sz1
        if(arr1d_in(j1)>0._r8)then
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr1d_in(j1)
        endif
      enddo
    enddo
  else
    ! Enforce expected array sizes   
    call shr_assert(sz2    == size(arr1d_in), errMsg(__FILE__, __LINE__))  

    do j2 = 1, sz2
      do j1 = 1, sz1
        if(arr1d_in(j2)>0._r8 .or. arr1d_in(j2)<0._r8)then
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr1d_in(j2)
        endif
      enddo
    enddo    

  endif
  return
  end subroutine array_div_vector_nofilter
  
end module SimpleMathMod
