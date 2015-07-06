!=============================================================================
! expression parser utility -- 
!   for parsing simple linear mathematical expressions of the form 
!   X = a*Y + b*Z + ... 
!
!=============================================================================
module shr_expr_parser_mod
  use shr_kind_mod,only : r8 => shr_kind_r8
  use shr_kind_mod,only : cx => shr_kind_cx

  implicit none
  private

  public :: shr_exp_parse     ! parses simple strings which contain expressions 
  public :: shr_exp_item_t    ! user defined type which contains an expression component
  public :: shr_exp_maxitems  ! number of objects contained in arry returned by shr_exp_parse

  integer, parameter :: shr_exp_maxitems = 256
  integer, parameter :: lrg_num = 256

  ! contains componets of expression
  type shr_exp_item_t
     character(len=64) :: name
     character(len=64) :: vars(lrg_num)
     real(r8)          :: coeffs(lrg_num)
     integer           :: n_terms = 0
  end type shr_exp_item_t

contains

  ! -----------------------------------------------------------------
  ! parses expressions provided in array of strings
  ! -----------------------------------------------------------------
  function shr_exp_parse( exp_array, nitems ) result(exp_items)

    character(len=*), intent(in)   :: exp_array(:) ! contains a expressions
    integer, optional, intent(out) :: nitems       ! number of expressions parsed
    type(shr_exp_item_t)           :: exp_items(shr_exp_maxitems) ! array of items returned
    
    integer :: i,j, nmax, n_exp_items
    character(len=cx) :: tmp_str

    n_exp_items = 0
    nmax = size( exp_array )
    
    do i = 1,nmax
       if (len_trim(exp_array(i))>0) then

          j = scan( exp_array(i), '=' )

          if ( j>0 ) then 

             n_exp_items = n_exp_items + 1
             exp_items(n_exp_items)%n_terms = 0

             exp_items(n_exp_items)%name = trim(adjustl(exp_array(i)(:j-1)))

             tmp_str = trim(adjustl(exp_array(i)(j+1:)))

             j = scan( tmp_str, '+' )

             if (j>0) then
                call set_coefvar( tmp_str(:j-1), exp_items(n_exp_items) )
                tmp_str = tmp_str(j-1:)
             else
                call set_coefvar( tmp_str, exp_items(n_exp_items) )
             endif

          else

             tmp_str = trim(adjustl(exp_array(i))) ! assumed to begin with '+'

          endif

          ! at this point tmp_str begins with '+'
          j = scan( tmp_str, '+' )

          if (j>0) then

             ! remove the leading + ...
             tmp_str = tmp_str(j+1:)
             j = scan( tmp_str, '+' )

             do while(j>0)

                call set_coefvar( tmp_str(:j-1), exp_items(n_exp_items) )

                tmp_str = tmp_str(j+1:)
                j = scan( tmp_str, '+' )

             enddo

             call set_coefvar( tmp_str, exp_items(n_exp_items) )

          endif

       endif
    enddo

    if ( present(nitems) ) then
       nitems = n_exp_items
    endif

  end function shr_exp_parse

 !==========================
 ! Private Methods 

  ! -----------------------------------------------------------------
  ! -----------------------------------------------------------------
  subroutine set_coefvar( term, item )
    character(len=*), intent(in)  :: term
    type(shr_exp_item_t) , intent(inout) :: item

    integer :: k, n

    item%n_terms = item%n_terms + 1
    n = item%n_terms

    k = scan( term, '*' )
    if (k>0) then
       item%vars(n) = trim(adjustl(term(k+1:)))
       read( term(:k-1), *) item%coeffs(n)
    else
       item%vars(n) = trim(adjustl(term))
       item%coeffs(n) = 1.0_r8
    endif

  end subroutine set_coefvar

end module shr_expr_parser_mod
