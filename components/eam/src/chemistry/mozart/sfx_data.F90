!-------------------------------------------------------------------
! reading surface fluxes from data similar to pg2 trcdata (zhang73)
!-------------------------------------------------------------------
module sfx_data
  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use cam_abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld,trfile
  use cam_logfile,  only : iulog
  !use ppgrid       , only : begchunk, endchunk
  !use physics_types, only : physics_state

  implicit none

  private  ! all unless made public
  save 

  public :: fields
  public :: sfx_data_init
  public :: sfx_data_adv
  public :: has_sfx_data
  public :: sfx_data_defaultopts
  public :: sfx_data_setopts
  public :: sfx_data_readnl                    ! read chem namelist 

  type(trfld), pointer :: fields(:) => null()
  type(trfile) :: file
  logical :: has_sfx_data 

  integer, parameter, public :: N_FLDS = 3
  integer :: number_flds

  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'CYCLICAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

! The named variables in fields have to match those in the sfx nc file; if one of them is mis-spelled, the pointer will point to a wrong one without error messages.
  character(len=16), dimension(N_FLDS), parameter :: fld_names = & ! data field names
       (/'lhf','shf','lwup'/)

  character(len=16), dimension(N_FLDS), parameter :: fld_units = & ! data field names
       (/'W m-2','W m-2','W m-2'/)

  integer :: index_map(N_FLDS)

  integer, public, parameter ::  lhf_ndx    =      1
  integer, public, parameter ::  shf_ndx   =      2
  integer, public, parameter ::  lwup_ndx   =      2

  
contains

!================================================================================================

!--(zhang73) mic components/eam/src/chemistry/mozart/chemistry.F90 <chem_readnl> 

  subroutine sfx_data_readnl(nlfile)

    ! Read chem namelist group.

    use cam_abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    !use sfx_data,       only: sfx_data_defaultopts,  sfx_data_setopts

   ! args

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    ! sfx data
    character(len=shr_kind_cl) :: sfx_data_file               ! prescribed data file
    character(len=shr_kind_cl) :: sfx_data_filelist           ! list of prescribed data files (series of files)
    character(len=shr_kind_cl) :: sfx_data_path               ! absolute path of prescribed data files 
    character(len=24)  :: sfx_data_type               ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)
    logical            :: sfx_data_rmfile             ! remove data file from local disk (default .false.)
    integer            :: sfx_data_cycle_yr
    integer            :: sfx_data_fixed_ymd
    integer            :: sfx_data_fixed_tod

    ! sfx inputs

    namelist /sfxdata_inparm/ &
         sfx_data_file, sfx_data_filelist, sfx_data_path, &
         sfx_data_type, &
         sfx_data_rmfile, sfx_data_cycle_yr, sfx_data_fixed_ymd, sfx_data_fixed_tod

    ! get the default settings

    call sfx_data_defaultopts( &
         sfx_data_file_out      = sfx_data_file,      &
         sfx_data_filelist_out  = sfx_data_filelist,  &
         sfx_data_path_out      = sfx_data_path,      &
         sfx_data_type_out      = sfx_data_type,      &
         sfx_data_rmfile_out    = sfx_data_rmfile,    &
         sfx_data_cycle_yr_out  = sfx_data_cycle_yr,  &
         sfx_data_fixed_ymd_out = sfx_data_fixed_ymd, &
         sfx_data_fixed_tod_out = sfx_data_fixed_tod  ) 

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'sfxdata_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, sfxdata_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('sfx_data_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables

    ! sfx data

    call mpibcast (sfx_data_file,      len(sfx_data_file),       mpichar, 0, mpicom)
    call mpibcast (sfx_data_filelist,  len(sfx_data_filelist),   mpichar, 0, mpicom)
    call mpibcast (sfx_data_path,      len(sfx_data_path),       mpichar, 0, mpicom)
    call mpibcast (sfx_data_type,      len(sfx_data_type),       mpichar, 0, mpicom)
    call mpibcast (sfx_data_rmfile,    1,                        mpilog, 0,  mpicom)
    call mpibcast (sfx_data_cycle_yr,  1,                        mpiint, 0,  mpicom)
    call mpibcast (sfx_data_fixed_ymd, 1,                        mpiint, 0,  mpicom)
    call mpibcast (sfx_data_fixed_tod, 1,                        mpiint, 0,  mpicom)
#endif

    ! set the options

   call sfx_data_setopts( &
        sfx_data_file_in      = sfx_data_file,      &
        sfx_data_filelist_in  = sfx_data_filelist,  &
        sfx_data_path_in      = sfx_data_path,      &
        sfx_data_type_in      = sfx_data_type,      &
        sfx_data_rmfile_in    = sfx_data_rmfile,    &
        sfx_data_cycle_yr_in  = sfx_data_cycle_yr,  &
        sfx_data_fixed_ymd_in = sfx_data_fixed_ymd, &
        sfx_data_fixed_tod_in = sfx_data_fixed_tod )

  endsubroutine sfx_data_readnl

!================================================================================================

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sfx_data_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, horiz_only
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: ndx, istat, i
    
    if ( has_sfx_data ) then
       if ( masterproc ) then
          write(iulog,*) 'sfx_data_ini: sfx data :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(fld_names)))
    file%in_pbuf(:) = .false.
    call trcdata_init( fld_names, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
        
    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'sfx_data_init: There are no sfx data'
          write(iulog,*) ' '
       endif
       return
    end if
    !write(iulog,*)'number_flds=',number_flds
    do i = 1,number_flds
       ndx = get_ndx( fields(i)%fldnam )
       index_map(i) = ndx
       write(iulog,*)'i=',i,'ndx=',ndx,'index_map(i)=',index_map(i)
       if (ndx < 1) then
          call endrun('sfx_data_init: '//trim(fields(i)%fldnam)//' is not one of the named sfx data fields ')
       endif
       !write(iulog,*)'fld_names(i)=',fld_names(i)

       call addfld(fld_names(i), horiz_only, 'I', fld_units(i), 'sfx data' )
       call addfld(trim(fld_names(i))//'_tot', horiz_only, 'I', fld_units(i), 'sfx data' )
    enddo


  end subroutine sfx_data_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sfx_data_setopts(&
       sfx_data_file_in,      &
       sfx_data_filelist_in,  &
       sfx_data_path_in,      &
       sfx_data_type_in,      &
       sfx_data_rmfile_in,    &
       sfx_data_cycle_yr_in,  &
       sfx_data_fixed_ymd_in, &
       sfx_data_fixed_tod_in  &
       )

    implicit none

    character(len=*), intent(in), optional :: sfx_data_file_in
    character(len=*), intent(in), optional :: sfx_data_filelist_in
    character(len=*), intent(in), optional :: sfx_data_path_in
    character(len=*), intent(in), optional :: sfx_data_type_in
    logical,          intent(in), optional :: sfx_data_rmfile_in
    integer,          intent(in), optional :: sfx_data_cycle_yr_in
    integer,          intent(in), optional :: sfx_data_fixed_ymd_in
    integer,          intent(in), optional :: sfx_data_fixed_tod_in

    if ( present(sfx_data_file_in) ) then
       filename = sfx_data_file_in
     endif
    if ( present(sfx_data_filelist_in) ) then
       filelist = sfx_data_filelist_in
     endif
    if ( present(sfx_data_path_in) ) then
       datapath = sfx_data_path_in
     endif
    if ( present(sfx_data_type_in) ) then
       datatype = sfx_data_type_in
     endif
    if ( present(sfx_data_rmfile_in) ) then
       rmv_file = sfx_data_rmfile_in
     endif
    if ( present(sfx_data_cycle_yr_in) ) then
       cycle_yr = sfx_data_cycle_yr_in
     endif
    if ( present(sfx_data_fixed_ymd_in) ) then
       fixed_ymd = sfx_data_fixed_ymd_in
     endif
    if ( present(sfx_data_fixed_tod_in) ) then
       fixed_tod = sfx_data_fixed_tod_in
    endif
    
      has_sfx_data =  .false.
      if (len_trim(filename) > 0 .or. len_trim(filelist) >0 )then
         if ((filename(1:3) .eq. 'sfx') .or. (filelist(1:3) .eq. 'sfx'))then
            has_sfx_data=   .true.
         endif
      endif
      if ( masterproc ) then
       write(iulog,*)'has_sfx_data=', has_sfx_data
       write(iulog,*)'filename=', filename
      endif

  endsubroutine sfx_data_setopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
 

  subroutine sfx_data_defaultopts(   &
       sfx_data_file_out,     &
       sfx_data_filelist_out, &
       sfx_data_path_out,     &
       sfx_data_type_out,     &
       sfx_data_rmfile_out,   &
       sfx_data_cycle_yr_out, &
       sfx_data_fixed_ymd_out,&
       sfx_data_fixed_tod_out &
       ) 

    implicit none

    character(len=*), intent(out), optional :: sfx_data_file_out
    character(len=*), intent(out), optional :: sfx_data_filelist_out
    character(len=*), intent(out), optional :: sfx_data_path_out
    character(len=*), intent(out), optional :: sfx_data_type_out
    logical,          intent(out), optional :: sfx_data_rmfile_out
    integer,          intent(out), optional :: sfx_data_cycle_yr_out
    integer,          intent(out), optional :: sfx_data_fixed_ymd_out
    integer,          intent(out), optional :: sfx_data_fixed_tod_out

    if ( present(sfx_data_file_out) ) then
       sfx_data_file_out = filename
    endif
    if ( present(sfx_data_filelist_out) ) then
       sfx_data_filelist_out = filelist
    endif
    if ( present(sfx_data_path_out) ) then
       sfx_data_path_out = datapath
    endif
    if ( present(sfx_data_type_out) ) then
       sfx_data_type_out = datatype
    endif
    if ( present(sfx_data_rmfile_out) ) then
       sfx_data_rmfile_out = rmv_file
    endif
    if ( present(sfx_data_cycle_yr_out) ) then
       sfx_data_cycle_yr_out = cycle_yr
    endif
    if ( present(sfx_data_fixed_ymd_out) ) then
       sfx_data_fixed_ymd_out = fixed_ymd
    endif
    if ( present(sfx_data_fixed_tod_out) ) then
       sfx_data_fixed_tod_out = fixed_tod
    endif

  endsubroutine sfx_data_defaultopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sfx_data_adv( pbuf2d, state, cam_in)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use cam_history,  only : outfld
    use physconst,    only : boltz                ! J/K/molecule
    use physics_buffer, only : physics_buffer_desc
    use camsrfexch,   only: cam_in_t        !(zhang73)       
    use physconst,    only: stebol, latvap  !(zhang73)
                                                              
    implicit none

  ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)  ! import state

  ! local vars
    integer :: ind,c,ncol,i
    real(r8) :: to_mmr(pcols,pver)

    if( .not. has_sfx_data ) return

    call advance_trcdata( fields, file, state, pbuf2d  )
    
    ! set the tracer fields with the correct units
    do i = 1,number_flds
       ind = index_map(i)
       do c = begchunk,endchunk
          ncol = state(c)%ncol
          call outfld( fields(i)%fldnam, fields(i)%data(:ncol,1,c), ncol, state(c)%lchnk )

          !(zhang73) update cam_in%lhf,shf,cflx.. here
          if ( trim(fld_names(i)) == 'lhf' )then
              cam_in(c)%lhf(:ncol) = cam_in(c)%lhf(:ncol) + fields(i)%data(:ncol,1,c)
              cam_in(c)%cflx(:ncol,1) = cam_in(c)%cflx(:ncol,1) + fields(i)%data(:ncol,1,c)/latvap
              call outfld( 'lhf_tot', cam_in(c)%lhf(:ncol), ncol, state(c)%lchnk )
          end if
          if ( trim(fld_names(i)) == 'shf' )then
              cam_in(c)%shf(:ncol) = cam_in(c)%shf(:ncol) + fields(i)%data(:ncol,1,c)
              call outfld( 'shf_tot', cam_in(c)%shf(:ncol), ncol, state(c)%lchnk )
          end if
          if ( trim(fld_names(i)) == 'lwup' )then
              cam_in(c)%lwup(:ncol) = cam_in(c)%lwup(:ncol) + fields(i)%data(:ncol,1,c)
              call outfld( 'lwup_tot', cam_in(c)%lwup(:ncol), ncol, state(c)%lchnk )
          end if
       enddo
    enddo

  end subroutine sfx_data_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  integer function get_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_ndx = 0
    do i = 1,N_FLDS
      if ( trim(name) == trim(fld_names(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module sfx_data
