!-------------------------------------------------------------------
! manages reading and interpolation of offline tracer fields
! Created by: Francis Vitt -- 2 May 2006
!-------------------------------------------------------------------
module tracer_cnst

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld,trfile,MAXTRCRS
  use cam_logfile,  only : iulog

  implicit none

  private  ! all unless made public
  save 

  public :: tracer_cnst_init
  public :: num_tracer_cnst
  public :: tracer_cnst_flds
  public :: tracer_cnst_adv
  public :: get_cnst_data
  public :: get_cnst_data_ptr
  public :: write_tracer_cnst_restart
  public :: read_tracer_cnst_restart
  public :: tracer_cnst_defaultopts
  public :: tracer_cnst_setopts
  public :: init_tracer_cnst_restart

  type(trfld), pointer :: fields(:) => null()
  type(trfile) :: file

  integer :: num_tracer_cnst
  character(len=16), pointer :: tracer_cnst_flds(:) => null()
  real(r8), allocatable, target, dimension(:,:,:,:) :: data_q  ! constituent mass mixing ratios

  character(len=64)  :: specifier(MAXTRCRS) = ''
  character(len=256) :: filename = 'tracer_cnst_file'
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_cnst_init()

    use mo_chem_utls,only : get_inv_ndx
    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: i ,ndx, istat

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

    num_tracer_cnst = 0
    if (associated(fields)) num_tracer_cnst = size( fields )

    if( num_tracer_cnst < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'There are no offline invariant species'
          write(iulog,*) ' '
       endif
       return
    end if

    allocate( tracer_cnst_flds(num_tracer_cnst), stat=istat)
    call handle_err(istat, 'tracer_cnst_init: ERROR allocating tracer_cnst_flds')


    do i = 1, num_tracer_cnst

       ndx = get_inv_ndx( fields(i)%fldnam )

       if (ndx < 1) then
          write(iulog,*) fields(i)%fldnam//' is not an invariant'
          call endrun('tracer_cnst_init')
       endif

       tracer_cnst_flds(i) = fields(i)%fldnam

       call addfld(trim(fields(i)%fldnam),'mol/mol ', pver, &
                   'I', 'prescribed tracer constituent', phys_decomp )
    enddo 

    allocate(data_q(pcols,pver,num_tracer_cnst,begchunk:endchunk), stat=istat)
    call handle_err(istat, 'tracer_cnst_init: ERROR allocating data_q')

  end subroutine tracer_cnst_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_cnst_setopts(       &
       tracer_cnst_file_in,      &
       tracer_cnst_filelist_in,  &
       tracer_cnst_datapath_in,  &
       tracer_cnst_type_in,      &
       tracer_cnst_specifier_in, &
       tracer_cnst_rmfile_in,    &
       tracer_cnst_cycle_yr_in,  &
       tracer_cnst_fixed_ymd_in, &
       tracer_cnst_fixed_tod_in  &
       )

    implicit none

    character(len=*), intent(in), optional :: tracer_cnst_file_in
    character(len=*), intent(in), optional :: tracer_cnst_filelist_in
    character(len=*), intent(in), optional :: tracer_cnst_datapath_in
    character(len=*), intent(in), optional :: tracer_cnst_type_in
    character(len=*), intent(in), optional :: tracer_cnst_specifier_in(:)
    logical,          intent(in), optional :: tracer_cnst_rmfile_in
    integer,          intent(in), optional :: tracer_cnst_cycle_yr_in
    integer,          intent(in), optional :: tracer_cnst_fixed_ymd_in
    integer,          intent(in), optional :: tracer_cnst_fixed_tod_in

    if ( present(tracer_cnst_file_in) ) then
       filename = tracer_cnst_file_in
    endif
    if ( present(tracer_cnst_filelist_in) ) then
       filelist = tracer_cnst_filelist_in
    endif
    if ( present(tracer_cnst_datapath_in) ) then
       datapath = tracer_cnst_datapath_in
    endif
    if ( present(tracer_cnst_type_in) ) then
       data_type = tracer_cnst_type_in
    endif
    if ( present(tracer_cnst_specifier_in) ) then
       specifier = tracer_cnst_specifier_in
    endif
    if ( present(tracer_cnst_rmfile_in) ) then
       rmv_file = tracer_cnst_rmfile_in
    endif
    if ( present(tracer_cnst_cycle_yr_in) ) then
       cycle_yr = tracer_cnst_cycle_yr_in
    endif
    if ( present(tracer_cnst_fixed_ymd_in) ) then
       fixed_ymd = tracer_cnst_fixed_ymd_in
    endif
    if ( present(tracer_cnst_fixed_tod_in) ) then
       fixed_tod = tracer_cnst_fixed_tod_in
    endif

  endsubroutine tracer_cnst_setopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_cnst_defaultopts(   &
       tracer_cnst_file_out,     &
       tracer_cnst_filelist_out, &
       tracer_cnst_datapath_out, &
       tracer_cnst_type_out,     &
       tracer_cnst_specifier_out,&
       tracer_cnst_rmfile_out,   &
       tracer_cnst_cycle_yr_out, &
       tracer_cnst_fixed_ymd_out,&
       tracer_cnst_fixed_tod_out &
       ) 

    implicit none

    character(len=*), intent(out), optional :: tracer_cnst_file_out
    character(len=*), intent(out), optional :: tracer_cnst_filelist_out
    character(len=*), intent(out), optional :: tracer_cnst_datapath_out
    character(len=*), intent(out), optional :: tracer_cnst_type_out
    character(len=*), intent(out), optional :: tracer_cnst_specifier_out(:)
    logical,          intent(out), optional :: tracer_cnst_rmfile_out
    integer,          intent(out), optional :: tracer_cnst_cycle_yr_out
    integer,          intent(out), optional :: tracer_cnst_fixed_ymd_out
    integer,          intent(out), optional :: tracer_cnst_fixed_tod_out

    if ( present(tracer_cnst_file_out) ) then
       tracer_cnst_file_out = filename
    endif
    if ( present(tracer_cnst_filelist_out) ) then
       tracer_cnst_filelist_out = filelist
    endif
    if ( present(tracer_cnst_datapath_out) ) then
       tracer_cnst_datapath_out = datapath
    endif
    if ( present(tracer_cnst_type_out) ) then
       tracer_cnst_type_out = data_type
    endif
    if ( present(tracer_cnst_specifier_out) ) then
       tracer_cnst_specifier_out = specifier
    endif
    if ( present(tracer_cnst_rmfile_out) ) then
       tracer_cnst_rmfile_out = rmv_file
    endif
    if ( present(tracer_cnst_cycle_yr_out) ) then
       tracer_cnst_cycle_yr_out = cycle_yr
    endif
    if ( present(tracer_cnst_fixed_ymd_out) ) then
       tracer_cnst_fixed_ymd_out = fixed_ymd
    endif
    if ( present(tracer_cnst_fixed_tod_out) ) then
       tracer_cnst_fixed_tod_out = fixed_tod
    endif

  endsubroutine tracer_cnst_defaultopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_cnst_adv( pbuf2d, state )

    use physics_buffer, only : physics_buffer_desc
    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use chem_mods,    only : fix_mass
    use mo_chem_utls, only : get_inv_ndx
    use cam_history,  only : outfld
    use physconst,    only: mwdry   ! molecular weight dry air ~ kg/kmole
    use physconst,    only: boltz

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: i,ind,c,ncol
    real(r8) :: to_vmr(pcols,pver)

    if( num_tracer_cnst < 1 ) return

    call advance_trcdata( fields, file, state, pbuf2d )

    ! copy prescribed tracer fields into state variable with the correct units

    do i = 1,num_tracer_cnst
       ind = get_inv_ndx( tracer_cnst_flds(i) )
       do c = begchunk,endchunk
          ncol = state(c)%ncol

          select case ( to_lower(trim(fields(i)%units(:GLC(fields(i)%units)))) )
          case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
             to_vmr(:ncol,:) = (1.e6_r8*boltz*state(c)%t(:ncol,:))/(state(c)%pmiddry(:ncol,:))
          case ('kg/kg','mmr')
             to_vmr(:ncol,:) = mwdry/fix_mass(ind)
          case ('mol/mol','mole/mole','vmr')
             to_vmr(:ncol,:) = 1._r8
          case default
             write(iulog,*) 'tracer_cnst_adv: units = ',trim(fields(i)%units) ,' are not recognized'
             call endrun('tracer_cnst_adv: units are not recognized')
          end select

          fields(i)%data(:ncol,:,c) = to_vmr(:ncol,:) * fields(i)%data(:ncol,:,c)      ! vmr
          call outfld( trim(tracer_cnst_flds(i)), fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk )

       enddo
    enddo

  end subroutine tracer_cnst_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_cnst_data( field_name, data, ncol, lchnk, pbuf  )

    use tracer_data, only : get_fld_data
    use physics_buffer, only : physics_buffer_desc

    implicit none

    character(len=*), intent(in) :: field_name
    real(r8), intent(out) :: data(:,:)
    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    if( num_tracer_cnst < 1 ) return

    call get_fld_data( fields, field_name, data, ncol, lchnk, pbuf  )

  end subroutine get_cnst_data

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_cnst_data_ptr(name, state, q, pbuf)

    use tracer_data,   only : get_fld_data, get_fld_ndx
    use physconst,     only : mwdry            ! molecular weight dry air ~ kg/kmole
    use chem_mods,     only : fix_mass
    use mo_chem_utls,  only : get_inv_ndx
    use physics_types, only : physics_state
    use ppgrid,        only : pcols, pver
    use physics_buffer, only : physics_buffer_desc

    implicit none

    character(len=*),    intent(in) :: name
    type(physics_state), intent(in) :: state
    real(r8), pointer, dimension(:,:) :: q     ! constituent mass mixing ratio
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: lchnk
    integer :: ncol
    integer :: inv_id, idx

    lchnk = state%lchnk
    ncol  = state%ncol

    ! make sure the requested constituent can be provided
    inv_id = get_inv_ndx(name) 
    if (.not. inv_id > 0) then
       if (masterproc) then
          write(iulog,*) 'get_cnst_data_ptr: '//name//' is not a prescribed tracer constituent'
       endif
       return
    endif


    call get_fld_ndx( fields, name, idx  )
    call get_fld_data( fields, name, data_q(:,:,idx,lchnk), ncol, lchnk, pbuf  )

    data_q(:ncol,:,idx,lchnk) = data_q(:ncol,:,idx,lchnk)*fix_mass(inv_id)/mwdry   ! vmr->mmr
    q => data_q(:,:,idx,lchnk)

  end subroutine get_cnst_data_ptr

!-------------------------------------------------------------------

  subroutine init_tracer_cnst_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'tracer_cnst', piofile, file )

  end subroutine init_tracer_cnst_restart
!-------------------------------------------------------------------
  subroutine write_tracer_cnst_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_tracer_cnst_restart

!-------------------------------------------------------------------
  subroutine read_tracer_cnst_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'tracer_cnst', piofile, file )

  end subroutine read_tracer_cnst_restart

end module tracer_cnst
