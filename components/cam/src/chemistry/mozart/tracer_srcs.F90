!-------------------------------------------------------------------
! manages reading and interpolation of offline tracer sources
! Created by: Francis Vitt -- 2 May 2006
!-------------------------------------------------------------------
module tracer_srcs

  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,   only : masterproc

  use tracer_data,  only : trfld,trfile,MAXTRCRS
  use cam_logfile,  only : iulog

  implicit none

  private  ! all unless made public
  save 

  public :: tracer_srcs_init
  public :: num_tracer_srcs
  public :: tracer_src_flds
  public :: tracer_srcs_adv
  public :: get_srcs_data
  public :: write_tracer_srcs_restart
  public :: read_tracer_srcs_restart
  public :: tracer_srcs_defaultopts
  public :: tracer_srcs_setopts
  public :: init_tracer_srcs_restart

  type(trfld), pointer :: fields(:) => null()
  type(trfile) :: file

  integer :: num_tracer_srcs
  character(len=16), allocatable :: tracer_src_flds(:)

  character(len=64)  :: specifier(MAXTRCRS) = ''
  character(len=256) :: filename = 'tracer_srcs_file'
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_init()

    use mo_chem_utls, only : get_extfrc_ndx
    use tracer_data,  only : trcdata_init
    use cam_history,  only : addfld
    use ppgrid,       only : pver
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: i ,ndx

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

    num_tracer_srcs = 0
    if (associated(fields)) num_tracer_srcs = size( fields )

    if( num_tracer_srcs < 1 ) then
       
       if (masterproc) then
          write(iulog,*) 'There are no offline tracer sources'
          write(iulog,*) ' '
       end if
       return
    end if

    allocate( tracer_src_flds(num_tracer_srcs))

    do i = 1, num_tracer_srcs

       ndx = get_extfrc_ndx( fields(i)%fldnam )

       if (ndx < 1) then
          write(iulog,*) fields(i)%fldnam//' is not configured to have an external source'
          call endrun('tracer_srcs_init')
       endif

       tracer_src_flds(i) = fields(i)%fldnam
 
       call addfld(trim(fields(i)%fldnam)//'_trsrc', (/ 'lev' /), 'I','/cm3/s', 'tracer source rate' )

    enddo 

  end subroutine tracer_srcs_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_setopts(       &
       tracer_srcs_file_in,      &
       tracer_srcs_filelist_in,  &
       tracer_srcs_datapath_in,  &
       tracer_srcs_type_in,      &
       tracer_srcs_specifier_in, &
       tracer_srcs_rmfile_in,    &
       tracer_srcs_cycle_yr_in,  &
       tracer_srcs_fixed_ymd_in, &
       tracer_srcs_fixed_tod_in  &
       )

    implicit none

    character(len=*), intent(in), optional :: tracer_srcs_file_in
    character(len=*), intent(in), optional :: tracer_srcs_filelist_in
    character(len=*), intent(in), optional :: tracer_srcs_datapath_in
    character(len=*), intent(in), optional :: tracer_srcs_type_in
    character(len=*), intent(in), optional :: tracer_srcs_specifier_in(:)
    logical,          intent(in), optional :: tracer_srcs_rmfile_in
    integer,          intent(in), optional :: tracer_srcs_cycle_yr_in
    integer,          intent(in), optional :: tracer_srcs_fixed_ymd_in
    integer,          intent(in), optional :: tracer_srcs_fixed_tod_in

    if ( present(tracer_srcs_file_in) ) then
       filename = tracer_srcs_file_in
    endif
    if ( present(tracer_srcs_filelist_in) ) then
       filelist = tracer_srcs_filelist_in
    endif
    if ( present(tracer_srcs_datapath_in) ) then
       datapath = tracer_srcs_datapath_in
    endif
    if ( present(tracer_srcs_type_in) ) then
       data_type = tracer_srcs_type_in
    endif
    if ( present(tracer_srcs_specifier_in) ) then
       specifier = tracer_srcs_specifier_in
    endif
    if ( present(tracer_srcs_rmfile_in) ) then
       rmv_file = tracer_srcs_rmfile_in
    endif
    if ( present(tracer_srcs_cycle_yr_in) ) then
       cycle_yr = tracer_srcs_cycle_yr_in
    endif
    if ( present(tracer_srcs_fixed_ymd_in) ) then
       fixed_ymd = tracer_srcs_fixed_ymd_in
    endif
    if ( present(tracer_srcs_fixed_tod_in) ) then
       fixed_tod = tracer_srcs_fixed_tod_in
    endif

  endsubroutine tracer_srcs_setopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_defaultopts(   &
       tracer_srcs_file_out,     &
       tracer_srcs_filelist_out, &
       tracer_srcs_datapath_out, &
       tracer_srcs_type_out,     &
       tracer_srcs_specifier_out,&
       tracer_srcs_rmfile_out,   &
       tracer_srcs_cycle_yr_out, &
       tracer_srcs_fixed_ymd_out,&
       tracer_srcs_fixed_tod_out &
       ) 

    implicit none

    character(len=*), intent(out), optional :: tracer_srcs_file_out
    character(len=*), intent(out), optional :: tracer_srcs_filelist_out
    character(len=*), intent(out), optional :: tracer_srcs_datapath_out
    character(len=*), intent(out), optional :: tracer_srcs_type_out
    character(len=*), intent(out), optional :: tracer_srcs_specifier_out(:)
    logical,          intent(out), optional :: tracer_srcs_rmfile_out
    integer,          intent(out), optional :: tracer_srcs_cycle_yr_out
    integer,          intent(out), optional :: tracer_srcs_fixed_ymd_out
    integer,          intent(out), optional :: tracer_srcs_fixed_tod_out

    if ( present(tracer_srcs_file_out) ) then
       tracer_srcs_file_out = filename
    endif
    if ( present(tracer_srcs_filelist_out) ) then
       tracer_srcs_filelist_out = filelist
    endif
    if ( present(tracer_srcs_datapath_out) ) then
       tracer_srcs_datapath_out = datapath
    endif
    if ( present(tracer_srcs_type_out) ) then
       tracer_srcs_type_out = data_type
    endif
    if ( present(tracer_srcs_specifier_out) ) then
       tracer_srcs_specifier_out = specifier
    endif
    if ( present(tracer_srcs_rmfile_out) ) then
       tracer_srcs_rmfile_out = rmv_file
    endif
    if ( present(tracer_srcs_cycle_yr_out) ) then
       tracer_srcs_cycle_yr_out = cycle_yr
    endif
    if ( present(tracer_srcs_fixed_ymd_out) ) then
       tracer_srcs_fixed_ymd_out = fixed_ymd
    endif
    if ( present(tracer_srcs_fixed_tod_out) ) then
       tracer_srcs_fixed_tod_out = fixed_tod
    endif

  endsubroutine tracer_srcs_defaultopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_adv( pbuf2d, state )

    use tracer_data, only : advance_trcdata
    use ppgrid,      only : begchunk, endchunk
    use physics_types,only : physics_state
    use cam_history, only : outfld
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: i,c,ncol

    if( num_tracer_srcs < 1 ) return

    call advance_trcdata( fields, file, state, pbuf2d )

    do c = begchunk,endchunk
       ncol = state(c)%ncol
       do i = 1,num_tracer_srcs
          call outfld( trim(fields(i)%fldnam)//'_trsrc', fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk  )
       enddo
    enddo

  end subroutine tracer_srcs_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_srcs_data( field_name, data, ncol, lchnk, pbuf  )

    use tracer_data, only : get_fld_data
    use physics_buffer, only : physics_buffer_desc

    implicit none

    character(len=*), intent(in) :: field_name
    real(r8), intent(out) :: data(:,:)
    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    if( num_tracer_srcs < 1 ) return

    call get_fld_data( fields, field_name, data, ncol, lchnk, pbuf )

  end subroutine get_srcs_data

!-------------------------------------------------------------------

  subroutine init_tracer_srcs_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'tracer_srcs', piofile, file )

  end subroutine init_tracer_srcs_restart
!-------------------------------------------------------------------
  subroutine write_tracer_srcs_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_tracer_srcs_restart

!-------------------------------------------------------------------

  subroutine read_tracer_srcs_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'tracer_srcs', piofile, file )

  end subroutine read_tracer_srcs_restart


end module tracer_srcs
