!-------------------------------------------------------------------
! manages reading and interpolation of prescribed volcanic aerosol
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_volcaero

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog
  use radconstants,   only: nswbands, nlwbands
  use volc_rad_data,  only: volc_rad_data_init, advance_volc_rad_data
  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_volcaero_readnl
  public :: prescribed_volcaero_register
  public :: prescribed_volcaero_init
  public :: prescribed_volcaero_adv
  public :: write_prescribed_volcaero_restart
  public :: read_prescribed_volcaero_restart
  public :: has_prescribed_volcaero
  public :: init_prescribed_volcaero_restart


  logical :: has_prescribed_volcaero = .false.
  character(len=8), parameter :: volcaero_name = 'VOLC_MMR'
  character(len=13), parameter :: volcrad_name = 'VOLC_RAD_GEOM'
  character(len=9), parameter :: volcmass_name = 'VOLC_MASS'
  character(len=11), parameter :: volcmass_column_name = 'VOLC_MASS_C'
  character(len=32), allocatable :: specifier(:)
  integer :: band_dim(6)
  ! Variables for CMIP6 style volcanic aerosols file
  character(len=12), parameter :: volc_ext_sun_name     = 'VOLC_EXT_SUN'    !extinction coefficient of solar bands in 1/km BALLI- change units appropriately 
  character(len=14), parameter :: volc_omega_sun_name   = 'VOLC_OMEGA_SUN'  !single scattering albedo of solar bands
  character(len=10), parameter :: volc_g_sun_name       = 'VOLC_G_SUN'      !asymmetry factor of solar bands

  character(len=14), parameter :: volc_ext_earth_name   = 'VOLC_EXT_EARTH'   !extinction coefficient of terrestrial bands in 1/km BALLI- change units appropriately
  character(len=16), parameter :: volc_omega_earth_name = 'VOLC_OMEGA_EARTH' !single scattering albedo of terrestrial bands
  character(len=12), parameter :: volc_g_earth_name     = 'VOLC_G_EARTH'     !asymmetry factor of terrestrial bands

  ! These variables are settable via the namelist (with longer names)
  character(len=16)  :: fld_name = 'MMRVOLC'
  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  character(len=32)  :: file_type = 'VOLC_MIXING_RATIO'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  integer            :: radius_ndx

  ! Variables settable via the namelist for CMIP6 style volcanic aerosols file
  character(len=16)  :: ext_sun_name     = 'ext_sun'
  character(len=16)  :: omega_sun_name   = 'omega_sun'
  character(len=16)  :: g_sun_name       = 'g_sun'

  character(len=16)  :: ext_earth_name   = 'ext_earth'
  character(len=16)  :: omega_earth_name = 'omega_earth'
  character(len=16)  :: g_earth_name     = 'g_earth'

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_volcaero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_volcaero_readnl'

   character(len=16)  :: prescribed_volcaero_name
   character(len=256) :: prescribed_volcaero_file
   character(len=256) :: prescribed_volcaero_filelist
   character(len=256) :: prescribed_volcaero_datapath
   character(len=32)  :: prescribed_volcaero_type
   character(len=32)  :: prescribed_volcaero_filetype
   logical            :: prescribed_volcaero_rmfile
   integer            :: prescribed_volcaero_cycle_yr
   integer            :: prescribed_volcaero_fixed_ymd
   integer            :: prescribed_volcaero_fixed_tod
   
   !local variables for CMIP6 style volcanic file

   character(len=16)  :: prescribed_volc_ext_sun_name   
   character(len=16)  :: prescribed_volc_omega_sun_name 
   character(len=16)  :: prescribed_volc_g_sun_name     
   
   character(len=16)  :: prescribed_volc_ext_earth_name  
   character(len=16)  :: prescribed_volc_omega_earth_name
   character(len=16)  :: prescribed_volc_g_earth_name    

   

   namelist /prescribed_volcaero_nl/    &
      prescribed_volcaero_name,         &
      prescribed_volcaero_file,         &
      prescribed_volcaero_filelist,     &
      prescribed_volcaero_datapath,     &
      prescribed_volcaero_type,         &
      prescribed_volcaero_rmfile,       &
      prescribed_volcaero_cycle_yr,     &
      prescribed_volcaero_fixed_ymd,    &
      prescribed_volcaero_fixed_tod,    &
      prescribed_volcaero_filetype,     &
      prescribed_volc_ext_sun_name,     &
      prescribed_volc_omega_sun_name,   &
      prescribed_volc_g_sun_name,       &
      prescribed_volc_ext_earth_name,   & 
      prescribed_volc_omega_earth_name, &
      prescribed_volc_g_earth_name

   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_volcaero_name     = fld_name
   prescribed_volcaero_file     = filename
   prescribed_volcaero_filelist = filelist
   prescribed_volcaero_datapath = datapath
   prescribed_volcaero_type     = data_type
   prescribed_volcaero_filetype = file_type
   prescribed_volcaero_rmfile   = rmv_file
   prescribed_volcaero_cycle_yr = cycle_yr
   prescribed_volcaero_fixed_ymd= fixed_ymd
   prescribed_volcaero_fixed_tod= fixed_tod

   prescribed_volc_ext_sun_name     = ext_sun_name
   prescribed_volc_omega_sun_name   = omega_sun_name
   prescribed_volc_g_sun_name       = g_sun_name 
   prescribed_volc_ext_earth_name   = ext_earth_name
   prescribed_volc_omega_earth_name = omega_earth_name
   prescribed_volc_g_earth_name     = g_earth_name


   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_volcaero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_volcaero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_volcaero_name,     len(prescribed_volcaero_name),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_file,     len(prescribed_volcaero_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_filelist, len(prescribed_volcaero_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_datapath, len(prescribed_volcaero_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_type,     len(prescribed_volcaero_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_filetype, len(prescribed_volcaero_filetype), mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_volcaero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volcaero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volcaero_fixed_tod,1, mpiint,  0, mpicom)

   call mpibcast(prescribed_volc_ext_sun_name,      len(prescribed_volc_ext_sun_name),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_omega_sun_name,    len(prescribed_volc_omega_sun_name),   mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_g_sun_name,        len(prescribed_volc_g_sun_name),       mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_ext_earth_name,    len(prescribed_volc_ext_earth_name),   mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_omega_earth_name,  len(prescribed_volc_omega_earth_name), mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_g_earth_name,      len(prescribed_volc_g_earth_name),     mpichar, 0, mpicom)
#endif

   ! Update module variables with user settings.
   fld_name   = prescribed_volcaero_name
   filename   = prescribed_volcaero_file
   filelist   = prescribed_volcaero_filelist
   datapath   = prescribed_volcaero_datapath
   data_type  = prescribed_volcaero_type
   file_type  = prescribed_volcaero_filetype
   rmv_file   = prescribed_volcaero_rmfile
   cycle_yr   = prescribed_volcaero_cycle_yr
   fixed_ymd  = prescribed_volcaero_fixed_ymd
   fixed_tod  = prescribed_volcaero_fixed_tod

   ext_sun_name     = prescribed_volc_ext_sun_name
   omega_sun_name   = prescribed_volc_omega_sun_name
   g_sun_name       = prescribed_volc_g_sun_name
   ext_earth_name   = prescribed_volc_ext_earth_name
   omega_earth_name = prescribed_volc_omega_earth_name
   g_earth_name     = prescribed_volc_g_earth_name



   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_volcaero = .true.

end subroutine prescribed_volcaero_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_register()
    use ppgrid,         only: pver,pcols
    use physics_buffer, only: pbuf_add_field, dtype_r8

    integer :: idx

    if (has_prescribed_volcaero) then
       if (trim(adjustl(file_type))== 'VOLC_CMIP6') then
          call pbuf_add_field(ext_sun_name,     'physpkg',dtype_r8,(/nswbands,pcols,pver/),idx) !BALLI- change hardwired 14 16 etc...
          call pbuf_add_field(omega_sun_name,   'physpkg',dtype_r8,(/nswbands,pcols,pver/),idx)
          call pbuf_add_field(g_sun_name,       'physpkg',dtype_r8,(/nswbands,pcols,pver/),idx)
          
          call pbuf_add_field(ext_earth_name,   'physpkg',dtype_r8,(/nlwbands,pcols,pver/),idx)
          call pbuf_add_field(omega_earth_name, 'physpkg',dtype_r8,(/nlwbands,pcols,pver/),idx)
          call pbuf_add_field(g_earth_name,     'physpkg',dtype_r8,(/nlwbands,pcols,pver/),idx)
          
       endif

       call pbuf_add_field(volcaero_name,'physpkg',dtype_r8,(/pcols,pver/),idx) !BALLI- we have to initialize it for radiation codes....but why????
       call pbuf_add_field(volcrad_name, 'physpkg',dtype_r8,(/pcols,pver/),idx)

    endif

  endsubroutine prescribed_volcaero_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_init()

    use tracer_data,    only: trcdata_init
    use cam_history,    only: addfld, horiz_only
    use ppgrid,         only: pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index

    implicit none

    integer :: ndx, istat
    integer :: errcode, ispf
    
    if ( has_prescribed_volcaero ) then
       if ( masterproc ) then
          write(iulog,*) 'volcanic aerosol is prescribed in :'//trim(filename)
       endif
    else
       return
    endif
    ispf = 1
    if (trim(adjustl(file_type))== 'VOLC_CMIP6') then
       allocate(specifier(6))
       
       specifier(ispf) = trim(adjustl(ext_sun_name))
       band_dim(ispf)  = nswbands; ispf = ispf + 1

       specifier(ispf) = trim(adjustl(omega_sun_name))
       band_dim(ispf)  = nswbands; ispf = ispf + 1

       specifier(ispf) = trim(adjustl(g_sun_name))
       band_dim(ispf)  = nswbands; ispf = ispf + 1

       specifier(ispf) = trim(adjustl(ext_earth_name))
       band_dim(ispf)  = nlwbands; ispf = ispf + 1

       specifier(ispf) = trim(adjustl(omega_earth_name))
       band_dim(ispf)  = nlwbands; ispf = ispf + 1

       specifier(ispf) = trim(adjustl(g_earth_name))
       band_dim(ispf)  = nlwbands

       
       !BALLI-add comments!!
       allocate(file%in_pbuf(size(specifier)))
       file%in_pbuf(:) = .true.
       call volc_rad_data_init(specifier, filename, datapath, band_dim, data_type)

    else if(trim(adjustl(file_type))== 'VOLC_MIXING_RATIO') then
       allocate(specifier(1))
       specifier(1) = trim(volcaero_name)//':'//trim(fld_name)
       
       allocate(file%in_pbuf(size(specifier)))
       file%in_pbuf(:) = .true.
       call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
            rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)
       
       call addfld(volcaero_name, (/ 'lev' /), 'I','kg/kg', 'prescribed volcanic aerosol dry mass mixing ratio' )
       call addfld(volcrad_name, (/ 'lev' /), 'I','m', 'volcanic aerosol geometric-mean radius' )
       call addfld(volcmass_name, (/ 'lev' /), 'I','kg/m^2', 'volcanic aerosol vertical mass path in layer' )
       call addfld(volcmass_column_name, horiz_only, 'I','kg/m^2', 'volcanic aerosol column mass' )
       
       radius_ndx = pbuf_get_index(volcrad_name, errcode)
    else
       call endrun('prescribed_volcaero_init: Invalid volcanic file type')
    endif



  end subroutine prescribed_volcaero_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz, gravit        ! J/K/molecule
    use tropopause,   only : tropopause_find, TROP_ALG_TWMO, TROP_ALG_CLIMATE
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    
    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer :: c,ncol,i,k
    real(r8) :: to_mmr(pcols,pver)
    real(r8), parameter :: molmass = 47.9981995_r8
    real(r8) :: ptrop
    real(r8) :: concvolc ! micrograms of wetted aerosol per cubic centimeter
    real(r8) :: volcmass(pcols,pver)
    real(r8) :: columnmass(pcols)
    real(r8) :: mmrvolc
    integer  :: tropLev(pcols)

    real(r8) :: outdata(pcols,pver)
    real(r8), pointer :: data(:,:)
    real(r8), pointer :: radius(:,:)

    !WACCM-derived relation between mass concentration and wet aerosol radius in meters
    real(r8),parameter :: radius_conversion = 1.9e-4_r8

    if( .not. has_prescribed_volcaero ) return

    if (trim(adjustl(file_type))== 'VOLC_CMIP6') then
       call advance_volc_rad_data (specifier, band_dim, state, pbuf2d)
    else if(trim(adjustl(file_type))== 'VOLC_MIXING_RATIO') then

       call advance_trcdata( fields, file, state, pbuf2d )
       
       ! copy prescribed tracer fields into state svariable with the correct units
       do c = begchunk,endchunk
          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          call pbuf_get_field(pbuf_chnk, radius_ndx, radius)
          radius(:,:) = 0._r8
          ncol = state(c)%ncol
          select case ( to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))) )
          case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
             to_mmr(:ncol,:) = (molmass*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
          case ('kg/kg','mmr','kg kg-1')
             to_mmr(:ncol,:) = 1._r8
          case ('mol/mol','mole/mole','vmr','fraction')
             to_mmr(:ncol,:) = molmass/mwdry
          case default
             write(iulog,*) 'prescribed_volcaero_adv: units = ',trim(fields(1)%units) ,' are not recognized'
             call endrun('prescribed_volcaero_adv: units are not recognized')
          end select
          
          call pbuf_get_field(pbuf_chnk, fields(1)%pbuf_ndx, data)
          data(:ncol,:) = to_mmr(:ncol,:) * data(:ncol,:) ! mmr
          
          call tropopause_find(state(c), tropLev, primary=TROP_ALG_TWMO, backup=TROP_ALG_CLIMATE)
          do i = 1,ncol
             do k = 1,pver
                ! set to zero below tropopause
                if ( k >= tropLev(i) ) then
                   data(i,k) = 0._r8
                endif
                mmrvolc = data(i,k)
                if (mmrvolc > 0._r8) then
                   concvolc = (mmrvolc * state(c)%pdel(i,k))/(gravit * state(c)%zm(i,k))
                   radius(i,k) = radius_conversion*(concvolc**(1._r8/3._r8))
                endif
             enddo
          enddo
          
          volcmass(:ncol,:) = data(:ncol,:)*state(c)%pdel(:ncol,:)/gravit
          columnmass(:ncol) = sum(volcmass(:ncol,:), 2)
          
          call outfld( volcaero_name,        data(:,:),     pcols, state(c)%lchnk)
          call outfld( volcrad_name,         radius(:,:),   pcols, state(c)%lchnk)
          call outfld( volcmass_name,        volcmass(:,:), pcols, state(c)%lchnk)
          call outfld( volcmass_column_name, columnmass(:), pcols, state(c)%lchnk)
          
       enddo
    endif

  end subroutine prescribed_volcaero_adv

!-------------------------------------------------------------------
  subroutine init_prescribed_volcaero_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_volcaero', piofile, file )

  end subroutine init_prescribed_volcaero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_volcaero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_volcaero_restart
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_volcaero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_volcaero', piofile, file )

  end subroutine read_prescribed_volcaero_restart

end module prescribed_volcaero
