module BeTRHistVarType

implicit none


type, public :: betr_hist_var_type
  character(len=64) :: varname
  character(len=18) :: units
  character(len=128):: long_name
  character(len=4)  :: avg_flag
  character(len=12) :: use_default
end type betr_hist_var_type
  public :: hist_var_copy
contains


  subroutine hist_var_copy(hist_var2, hist_var1)

  implicit none
  type(betr_hist_var_type), intent(in) :: hist_var1
  type(betr_hist_var_type), intent(out) :: hist_var2

  hist_var2%varname= hist_var1%varname
  hist_var2%units= hist_var1%units
  hist_var2%long_name=hist_var1%long_name
  hist_var2%avg_flag = hist_var1%avg_flag
  hist_var2%use_default=hist_var1%use_default

  end subroutine hist_var_copy

end module BeTRHistVarType
