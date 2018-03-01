module FatesUtilsMod
  
  ! This module contains helper functions and subroutines which are general in nature.
  ! Think string parsing, timing, maybe numerics, etc.
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals, only      : fates_log

contains
  
  
  function check_hlm_list(hlms,hlm_name) result(astatus)
    
    ! ---------------------------------------------------------------------------------
    ! This simple function compares a string of HLM tags to see if any of the names
    ! match the name of the currently active HLM. If any do, return true, if any
    ! don't, if any don't its a big secret.
    ! ---------------------------------------------------------------------------------
    
    character(len=*),intent(in) :: hlms
    character(len=*),intent(in) :: hlm_name

    integer :: index
    logical :: astatus

    astatus = .false.
    index = scan(trim(hlms),trim(hlm_name))
 
    if(index>0)then
       astatus=.true.
    end if
    return

  end function check_hlm_list
  
  ! =====================================================================================

  subroutine check_var_real(r8_var, var_name, return_code)

     real(r8),intent(in)         :: r8_var
     character(len=*),intent(in) :: var_name
     integer,intent(out)         :: return_code

     real(r8), parameter :: r8_type  = 1.0
     real(r8), parameter :: overflow  = huge(r8_type)
     real(r8), parameter :: underflow = tiny(r8_type)

     return_code = 0

     ! NaN check
     if (r8_var /= r8_var) then
        write(fates_log(),*) 'NaN detected, ',trim(var_name),': ',r8_var
        return_code = 1
     end if
     
     ! Overflow check (within 100th of max precision)
     if (abs(r8_var) > 0.01*overflow) then
        write(fates_log(),*) 'Nigh overflow detected, ',trim(var_name),': ',r8_var
        return_code = return_code + 10
     end if

     ! Underflow check (within 100x of min precision)
     if (abs(r8_var) < 100.0_r8*underflow) then
        write(fates_log(),*) 'Nigh underflow detected, ',trim(var_name),': ',r8_var
        return_code = return_code + 100
     end if


  end subroutine check_var_real


end module FatesUtilsMod
