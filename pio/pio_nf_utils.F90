module pio_nf_utils
  use pio_types, only : file_desc_t, var_desc_t
  use nf_mod, only : pio_inq_vartype=>inq_vartype
  use pionfget_mod, only : pio_get_var=>get_var
  use pionfput_mod, only : pio_put_var=>put_var
  use pio_types, only : pio_int, pio_real, pio_double, pio_char
  use pio_kinds, only : i4, r4, r8
  use pio_support, only : piodie

  implicit none
  private
  public :: copy_pio_var
  interface copy_pio_var
     module procedure copy_pio_var01d
     module procedure copy_pio_var2d
  end interface
contains

subroutine copy_pio_var01d(ifh, ofh, ivid, ovid, length, strlength)
  type(File_Desc_t) :: Ifh, Ofh
  type(Var_Desc_t) ::  ivid, ovid
  integer, intent(in) :: length
  integer, intent(in), optional :: strlength
  integer(i4), allocatable :: ival(:)
  real(r4), allocatable :: rval(:)
  real(r8), allocatable :: dval(:)
  character(len=length), allocatable :: cval(:)

  integer :: ierr
  
!  ierr = pio_inq_vartype(ifh, ivid, itype)
!  ierr = pio_inq_vartype(ofh, ovid, otype)

  
  if( ivid%type.ne.ovid%type) then
     write(6,*) 'WARNING: copy_pio_var coercing type ', ivid%type, ' to ',ovid%type
  end if
  select case(ivid%type)
  case (PIO_int)
     allocate(ival(length))
     ierr = pio_get_var(ifh, ivid%varid, ival)
     ierr = pio_put_var(ofh, ovid%varid, ival)
     deallocate(ival)
  case (PIO_real)
     allocate(rval(length))
     ierr = pio_get_var(ifh, ivid%varid, rval)
     ierr = pio_put_var(ofh, ovid%varid, rval)
     deallocate(rval)
  case (PIO_double)
     allocate(dval(length))
     ierr = pio_get_var(ifh, ivid%varid, dval)
     ierr = pio_put_var(ofh, ovid%varid, dval)
     deallocate(dval)
  case (PIO_char)
     if(present(strlength)) then
        allocate(cval(strlength))
     else
        allocate(cval(1))
     end if
     ierr = pio_get_var(ifh, ivid%varid, cval)
     ierr = pio_put_var(ofh, ovid%varid, cval)

     deallocate(cval)
  end select
end subroutine copy_pio_var01d

subroutine copy_pio_var2d(ifh, ofh, ivid, ovid, length)
  type(File_Desc_t) :: Ifh, Ofh
  type(Var_Desc_t) ::  ivid, ovid
  integer, intent(in) :: length(:)
!  integer, intent(in), optional :: strlength
  integer(i4), allocatable :: ival(:,:)
  real(r4), allocatable :: rval(:,:)
  real(r8), allocatable :: dval(:,:)
!  character(len=length), allocatable :: cval(:,:)

  integer :: ierr
  
!  ierr = pio_inq_vartype(ifh, ivid, itype)
!  ierr = pio_inq_vartype(ofh, ovid, otype)

  
  if( ivid%type.ne.ovid%type) then
     write(6,*) 'WARNING: copy_pio_var coercing type ', ivid%type, ' to ',ovid%type
  end if
  select case(ivid%type)
  case (PIO_int)
     allocate(ival(length(1),length(2)))
     ierr = pio_get_var(ifh, ivid%varid, ival)
     ierr = pio_put_var(ofh, ovid%varid, ival)
     deallocate(ival)
  case (PIO_real)
     allocate(rval(length(1),length(2)))
     ierr = pio_get_var(ifh, ivid%varid, rval)
     ierr = pio_put_var(ofh, ovid%varid, rval)
     deallocate(rval)
  case (PIO_double)
     allocate(dval(length(1),length(2)))
     ierr = pio_get_var(ifh, ivid%varid, dval)
     ierr = pio_put_var(ofh, ovid%varid, dval)
     deallocate(dval)
  case (PIO_char)
     !        if(present(strlength)) then
     !           allocate(cval(strlength))
     !        else
     !           allocate(cval(1))
     !        end if
     !        ierr = pio_get_var(ifh, ivid%varid, cval)
     !        ierr = pio_put_var(ofh, ovid%varid, cval)
     !        deallocate(cval)
  end select
end subroutine copy_pio_var2d


end module pio_nf_utils
