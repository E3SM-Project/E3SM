!-------------------------------------------------------------------
! manages reading and interpolation of prescribed ghg tracers
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_ghg

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_ghg_init
  public :: prescribed_ghg_adv
  public :: write_prescribed_ghg_restart
  public :: read_prescribed_ghg_restart
  public :: has_prescribed_ghg
  public :: prescribed_ghg_register
  public :: init_prescribed_ghg_restart
  public :: prescribed_ghg_readnl

  logical :: has_prescribed_ghg = .false.
  integer, parameter, public :: N_GHG = 5
  integer :: number_flds

  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=16)  :: specifier(N_GHG) = ''

  character(len=8)    :: ghg_names(N_GHG) = (/ 'prsd_co2',  'prsd_ch4',  'prsd_n2o',  'prsd_f11',  'prsd_f12'  /)
  real(r8), parameter :: molmass(N_GHG)   = (/ 44.00980_r8, 16.04060_r8, 44.01288_r8, 137.3675_r8, 120.9132_r8 /)

  integer :: index_map(N_GHG)

contains


!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_ghg_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_ghg_readnl'

   character(len=16)  :: prescribed_ghg_specifier(N_GHG)
   character(len=256) :: prescribed_ghg_file
   character(len=256) :: prescribed_ghg_filelist
   character(len=256) :: prescribed_ghg_datapath
   character(len=32)  :: prescribed_ghg_type
   logical            :: prescribed_ghg_rmfile
   integer            :: prescribed_ghg_cycle_yr
   integer            :: prescribed_ghg_fixed_ymd
   integer            :: prescribed_ghg_fixed_tod

   namelist /prescribed_ghg_nl/ &
      prescribed_ghg_specifier, &
      prescribed_ghg_file,      &
      prescribed_ghg_filelist,  &
      prescribed_ghg_datapath,  &
      prescribed_ghg_type,      &
      prescribed_ghg_rmfile,    &
      prescribed_ghg_cycle_yr,  &
      prescribed_ghg_fixed_ymd, &
      prescribed_ghg_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_ghg_specifier= specifier
   prescribed_ghg_file     = filename
   prescribed_ghg_filelist = filelist
   prescribed_ghg_datapath = datapath
   prescribed_ghg_type     = datatype
   prescribed_ghg_rmfile   = rmv_file
   prescribed_ghg_cycle_yr = cycle_yr
   prescribed_ghg_fixed_ymd= fixed_ymd
   prescribed_ghg_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_ghg_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_ghg_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_ghg_specifier,len(prescribed_ghg_specifier(1))*N_GHG,     mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_file,     len(prescribed_ghg_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_filelist, len(prescribed_ghg_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_datapath, len(prescribed_ghg_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_type,     len(prescribed_ghg_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_ghg_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_ghg_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_ghg_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = prescribed_ghg_specifier
   filename   = prescribed_ghg_file
   filelist   = prescribed_ghg_filelist
   datapath   = prescribed_ghg_datapath
   datatype   = prescribed_ghg_type
   rmv_file   = prescribed_ghg_rmfile
   cycle_yr   = prescribed_ghg_cycle_yr
   fixed_ymd  = prescribed_ghg_fixed_ymd
   fixed_tod  = prescribed_ghg_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_ghg = .true.

end subroutine prescribed_ghg_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ghg_register()
    use ppgrid,         only: pver, pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: i,idx

    if (has_prescribed_ghg) then
       do i = 1,N_GHG
          call pbuf_add_field(ghg_names(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
       enddo
    endif

  endsubroutine prescribed_ghg_register
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ghg_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: ndx, istat, i
    
    if ( has_prescribed_ghg ) then
       if ( masterproc ) then
          write(iulog,*) 'ghg is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
        
    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'There are no prescribed ghg tracers'
          write(iulog,*) ' '
       endif
       return
    end if

    do i = 1,number_flds
       ndx = get_ndx( fields(i)%fldnam )
       index_map(i) = ndx

       if (ndx < 1) then
          call endrun('prescribed_ghg_init: '//trim(fields(i)%fldnam)//' is not one of the named ghg fields in pbuf2d')
       endif
       call addfld( fields(i)%fldnam, (/ 'lev' /), 'I','kg/kg', 'prescribed ghg' )
    enddo

  end subroutine prescribed_ghg_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ghg_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz                ! J/K/molecule
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_set_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    integer :: ind,c,ncol,i
    real(r8) :: to_mmr(pcols,pver)
    real(r8) :: outdata(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)

    character(len=32) :: units_str

    if( .not. has_prescribed_ghg ) return

    call advance_trcdata( fields, file, state, pbuf2d )
    
    ! set the correct units and invoke history outfld
    do i = 1,number_flds
       ind = index_map(i)

       units_str = trim(to_lower(trim(fields(i)%units(:GLC(fields(i)%units)))))

!$OMP PARALLEL DO PRIVATE (C, NCOL, OUTDATA, TO_MMR, tmpptr, pbuf_chnk)
       do c = begchunk,endchunk
          ncol = state(c)%ncol

          select case ( units_str )
          case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
             to_mmr(:ncol,:) = (molmass(ind)*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
          case ('kg/kg','mmr')
             to_mmr(:ncol,:) = 1._r8
          case ('mol/mol','mole/mole','vmr','fraction')
             to_mmr(:ncol,:) = molmass(ind)/mwdry
          case default
             print*, 'prescribed_ghg_adv: units = ',trim(fields(i)%units) ,' are not recognized'
             call endrun('prescribed_ghg_adv: units are not recognized')
          end select

          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          call pbuf_get_field(pbuf_chnk, fields(i)%pbuf_ndx, tmpptr )

          tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)

          outdata(:ncol,:) = tmpptr(:ncol,:) 
          call outfld( fields(1)%fldnam, outdata(:ncol,:), ncol, state(c)%lchnk )

       enddo
    enddo

  end subroutine prescribed_ghg_adv

!-------------------------------------------------------------------

!-------------------------------------------------------------------
  subroutine init_prescribed_ghg_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_ghg', piofile, file )

  end subroutine init_prescribed_ghg_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_ghg_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_ghg_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_ghg_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_ghg', piofile, file )

  end subroutine read_prescribed_ghg_restart
!-------------------------------------------------------------------
  integer function get_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_ndx = 0
    do i = 1,N_GHG
      if ( trim(name) == trim(ghg_names(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module prescribed_ghg
