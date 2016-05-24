subroutine prune_trailing_numbers(instr,outstr)
implicit none

  character(len=*)instr
  character(len=*)outstr
  integer inlen
  integer i;

  inlen=LEN(TRIM(instr))

  i=inlen
  do while(instr(i:i) == "0" .or. &
           instr(i:i) == "1" .or. &
           instr(i:i) == "2" .or. &
           instr(i:i) == "3" .or. &
           instr(i:i) == "4" .or. &
           instr(i:i) == "5" .or. &
           instr(i:i) == "6" .or. &
           instr(i:i) == "7" .or. &
           instr(i:i) == "8" .or. &
           instr(i:i) == "9" .or. &
           instr(i:i) == "." )
      i=i-1
  end do

  outstr=instr(1:i)

end subroutine prune_trailing_numbers
