#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module common_io_mod
  use control_mod, only : MAX_STRING_LEN         !HOMME Specific: MAX_STRING_LEN
#if defined( PIO ) || defined ( PIO_INTERP )
  use pio, only : var_desc_t, file_desc_t, io_desc_t, nfsizekind=>PIO_OffSet, iosystem_desc_t, & ! _EXTERNAL
       nf_double=>pio_double, nf_int=>pio_int, unlim_dim=>pio_unlimited, nf_noerr=>pio_noerr
#endif


  implicit none
  private

  public :: nfsizekind, nf_noerr

  integer, parameter, public :: max_output_streams=5
  integer, parameter, public :: max_output_variables=50
  integer, parameter, public :: varname_len=16
  ! output analysis namelists, 


  integer, parameter, public :: MAX_INFILES=120
!  integer, parameter, public :: MAX_VECVARS=10
  
  character(len=160), public, target :: infilenames(MAX_INFILES)
!  character(len=10), public :: vector_uvars(MAX_VECVARS), vector_vvars(MAX_VECVARS)



  integer, public :: output_timeunits(max_output_streams) 
  !  0=steps, 1=pday, 2=phour
  integer, public :: output_start_time(max_output_streams)
  integer, public :: output_end_time(max_output_streams)
  integer, public :: output_frequency(max_output_streams)
  character(len=MAX_STRING_LEN), public :: output_prefix
  character(len=MAX_STRING_LEN), public :: output_dir
  character(len=9), public :: output_type 

  character(len=varname_len), public, target,  dimension(max_output_variables) :: output_varnames1   
  character(len=varname_len), public, target,  dimension(max_output_variables) :: output_varnames2
  character(len=varname_len), public, target,  dimension(max_output_variables) :: output_varnames3
  character(len=varname_len), public, target,  dimension(max_output_variables) :: output_varnames4
  character(len=varname_len), public, target,  dimension(max_output_variables) :: output_varnames5

  ! these are used for PIO output method
  integer, public :: num_io_procs
  integer, public :: num_agg
  integer, public :: io_stride

  ! end of analysis_nl namelist variables
  public :: nf_selectedvar
  public :: get_varindex
  public :: get_current_varnames
  public :: nf_addrequiredvar
  public :: nf_dim, nf_variable, nf_handle, nf_int, nf_double
  public :: unlim_dim




  integer, parameter, public :: beginstate=1, dimsstate=2,varsstate=3,readystate=4 

  type nf_variable
#ifdef SEQUENCE
     sequence
#endif
     character(len=varname_len) :: varname
     integer :: ivarID
     integer :: ndims
     logical :: required
     logical :: timedependent
     integer :: dimsid
     type(Var_desc_t) :: Vardesc   
  end type nf_variable

  type nf_dim
#ifdef SEQUENCE
     sequence
#endif
     character(len=varname_len) :: dimname
     integer :: dimID
     integer :: dsize
  end type nf_dim
  !
  ! output file information is contained in the linked list nf_handle, that list in turn contains
  !  a linked list of output variable information for that file.


  public :: nf_decomp
  type nf_decomp
#ifdef SEQUENCE
     sequence
#endif
     type(IO_desc_t) :: Iodesc
     ! This is an identifier that should be unique to the combination of dims that describe this decomposition
     integer :: dimsid     
  end type nf_decomp

  type nf_handle
     !     private
#ifdef SEQUENCE
     sequence
#endif
     type(File_desc_t) :: FileID
     type(nf_variable), pointer :: varlist(:)
     type(nf_dim), pointer :: dimlist(:)
     type(nf_decomp), pointer :: decomplist(:) 
     integer :: ncFileID
     integer :: timeDimID
     integer :: timevarID
     integer :: iframe     ! frame count for this file
     integer :: state      ! used for sanity check

  end type nf_handle

  type(iosystem_desc_t), public :: PIOFS   ! just one of these needed 


contains




  !
  !  Returns a pointer to the list of variables to be written to file ios
  !

  function get_current_varnames(ios)
    integer, intent(in) :: ios
    character(len=varname_len), pointer :: get_current_varnames(:)
    select case(ios)
    case (1)
       get_current_varnames=>output_varnames1
    case (2)
       get_current_varnames=>output_varnames2
    case (3)
       get_current_varnames=>output_varnames3
    case (4)
       get_current_varnames=>output_varnames4
    case (5)
       get_current_varnames=>output_varnames5
    end select
  end function get_current_varnames

  !  
  !  Given a var name this gives an index into an nf_variable array varlist that matchs that 
  !  variable, returns -1 if no match is found 
  !
  function get_varindex(name, varlist)
    character(len=*), intent(in) :: name
    type(nf_variable), intent(in) :: varlist(:)
    integer :: get_varindex
    integer :: i, vmax

    get_varindex=-1
    vmax = size(varlist)
    do i=1,vmax
       if(name .eq. varlist(i)%varname) then
          get_varindex=i
          exit
       end if
    end do
  end function get_varindex

  function nf_selectedvar(var, output_list) result(selected)
    character*(*), intent(in) :: var, output_list(:)
    logical ::  selected
    integer :: i
    selected=.false.

    do i=1,max_output_variables
       if(var.eq.output_list(i)) then
          selected=.true.
          exit
       end if
    enddo

  end function nf_selectedvar


  subroutine nf_addrequiredvar(ncdf_list,name)
    use parallel_mod, only : abortmp
    character(len=*), intent(in) :: name
    type(nf_handle), intent(in), target :: ncdf_list(:)
    type(nf_handle), pointer :: ncdf
    integer :: ios, i
    character(len=varname_len), pointer :: output_varnames(:)

! add a variable to the list of required variables
    do ios=1,max_output_streams
       ncdf=>ncdf_list(ios)
       if(ncdf%state == dimsstate .or. ncdf%state == varsstate) then
          output_varnames => get_current_varnames(ios)
          
          do i=1,max_output_variables
             if(len(trim(output_varnames(i)))==0) then
                output_varnames(i)=name
                !                print *, 'variable ',name, ' added at ',i
                exit
             end if
          end do
       end if
    end do

    
  end subroutine nf_addrequiredvar


end module common_io_mod
