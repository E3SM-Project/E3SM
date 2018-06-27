module mkabortutils

  ! Abort the model for abnormal termination
   private
   save

   public  :: endrun

CONTAINS

  subroutine endrun(msg)
    implicit none
    character(len=*), intent(in), optional :: msg    ! string to be printed
    
    if (present (msg)) then
       write(6,*)'ENDRUN:', msg
    else
       write(6,*)'ENDRUN: called without a message string'
    end if
    stop
    
  end subroutine endrun
  
end module mkabortutils
