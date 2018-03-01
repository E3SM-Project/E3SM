module debuginfoMod
implicit none
  private
  public :: write_subname
  logical:: ldebug=.true.
contains

  subroutine write_subname(subname)
  !
  ! DESCRIPTIONS
  ! write subroutine name when the model is in debugging mode
  implicit none
  character(len=*), intent(in) :: subname
  
  
  
  if(ldebug)write(*,*)subname
  
  end subroutine write_subname  
  

end module debuginfoMod
