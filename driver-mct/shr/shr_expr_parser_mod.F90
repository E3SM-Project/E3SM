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

  ! contains componets of expression
  type shr_exp_item_t
     character(len=64) :: name
     character(len=64),pointer :: vars(:) => null()
     real(r8)         ,pointer :: coeffs(:) => null()
     integer           :: n_terms = 0
     type(shr_exp_item_t), pointer :: next_item
  end type shr_exp_item_t

contains

  ! -----------------------------------------------------------------
  ! parses expressions provided in array of strings
  ! -----------------------------------------------------------------
  function shr_exp_parse( exp_array, nitems ) result(exp_items_list)

    character(len=*), intent(in)   :: exp_array(:) ! contains a expressions
    integer, optional, intent(out) :: nitems       ! number of expressions parsed
    type(shr_exp_item_t), pointer  :: exp_items_list ! linked list of items returned
    
    integer :: i,j, jj, nmax, nterms, n_exp_items
    character(len=cx) :: tmp_str
    type(shr_exp_item_t), pointer :: exp_item, list_item

    nullify( exp_items_list )

    n_exp_items = 0
    nmax = size( exp_array )
    
    do i = 1,nmax
       if (len_trim(exp_array(i))>0) then

          j = scan( exp_array(i), '=' )

          if ( j>0 ) then 

             n_exp_items = n_exp_items + 1

             allocate( exp_item )
             exp_item%n_terms = 0
             exp_item%name = trim(adjustl(exp_array(i)(:j-1)))

             tmp_str = trim(adjustl(exp_array(i)(j+1:)))

             nterms = 1             
             jj = scan( tmp_str, '+' )
             do while(jj>0)
                nterms = nterms + 1
                tmp_str = tmp_str(jj+1:)
                jj = scan( tmp_str, '+' )
             enddo
             
             allocate( exp_item%vars(nterms) )
             allocate( exp_item%coeffs(nterms) )

             tmp_str = trim(adjustl(exp_array(i)(j+1:)))

             j = scan( tmp_str, '+' )

             if (j>0) then
                call set_coefvar( tmp_str(:j-1), exp_item )
                tmp_str = tmp_str(j-1:)
             else
                call set_coefvar( tmp_str, exp_item )
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

                call set_coefvar( tmp_str(:j-1), exp_item )

                tmp_str = tmp_str(j+1:)
                j = scan( tmp_str, '+' )

             enddo

             call set_coefvar( tmp_str, exp_item )

          endif

          
          if (associated(exp_item)) then
             if (associated(exp_items_list)) then
                list_item => exp_items_list
                do while(associated(list_item%next_item))
                   list_item => list_item%next_item
                enddo
                list_item%next_item => exp_item
             else
                exp_items_list => exp_item
             endif
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
