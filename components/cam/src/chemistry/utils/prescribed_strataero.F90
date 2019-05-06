!-------------------------------------------------------------------
! manages reading and interpolation of prescribed stratospheric aerosols
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_strataero

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld, trfile
  use cam_logfile,      only : iulog

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_strataero_readnl
  public :: prescribed_strataero_register
  public :: prescribed_strataero_init
  public :: prescribed_strataero_adv
  public :: write_prescribed_strataero_restart
  public :: read_prescribed_strataero_restart
  public :: has_prescribed_strataero
  public :: init_prescribed_strataero_restart

  logical :: has_prescribed_strataero = .false.
  character(len=16), parameter :: mmr_name = 'VOLC_MMR'
  character(len=16), parameter :: rad_name = 'VOLC_RAD_GEOM'
  character(len=16), parameter :: sad_name = 'VOLC_SAD'
  character(len=16), parameter :: mass_name = 'VOLC_MASS'
  character(len=16), parameter :: mass_column_name = 'VOLC_MASS_C'
  character(len=16), parameter :: dens_name = 'VOLC_DENS'

  character(len=16), parameter :: mmr_name1 = 'VOLC_MMR1'
  character(len=16), parameter :: mmr_name2 = 'VOLC_MMR2'
  character(len=16), parameter :: mmr_name3 = 'VOLC_MMR3'
  character(len=16), parameter :: rad_name1 = 'VOLC_RAD_GEOM1'
  character(len=16), parameter :: rad_name2 = 'VOLC_RAD_GEOM2'
  character(len=16), parameter :: rad_name3 = 'VOLC_RAD_GEOM3'
  character(len=16), parameter :: mass_name1 = 'VOLC_MASS1'
  character(len=16), parameter :: mass_name2 = 'VOLC_MASS2'
  character(len=16), parameter :: mass_name3 = 'VOLC_MASS3'
  character(len=16), parameter :: mass_column_name1 = 'VOLC_MASS_C1'
  character(len=16), parameter :: mass_column_name2 = 'VOLC_MASS_C2'
  character(len=16), parameter :: mass_column_name3 = 'VOLC_MASS_C3'
  character(len=16), parameter :: dens_name1 = 'VOLC_DENS1'
  character(len=16), parameter :: dens_name2 = 'VOLC_DENS2'
  character(len=16), parameter :: dens_name3 = 'VOLC_DENS3'

  ! These variables are settable via the namelist (with longer names)
  character(len=32)  :: specifier(7) = ' '
  character(len=256) :: filename = 'NONE'
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  integer            :: rad_ndx1 = -1
  integer            :: rad_ndx2 = -1
  integer            :: rad_ndx3 = -1
  integer            :: sad_ndx = -1
  integer            :: mmr_ndx1 = -1
  integer            :: mmr_ndx2 = -1
  integer            :: mmr_ndx3 = -1

  logical            :: prescribed_strataero_use_chemtrop = .false.
  logical            :: three_mode = .true.
  integer :: rad_fld_no=-1, sad_fld_no=-1

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_strataero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_strataero_readnl'

   character(len=32)  :: prescribed_strataero_specifier(7)
   character(len=256) :: prescribed_strataero_file
   character(len=256) :: prescribed_strataero_filelist
   character(len=256) :: prescribed_strataero_datapath
   character(len=32)  :: prescribed_strataero_type
   logical            :: prescribed_strataero_rmfile
   integer            :: prescribed_strataero_cycle_yr
   integer            :: prescribed_strataero_fixed_ymd
   integer            :: prescribed_strataero_fixed_tod

   namelist /prescribed_strataero_nl/ &
      prescribed_strataero_specifier, &
      prescribed_strataero_file,      &
      prescribed_strataero_filelist,  &
      prescribed_strataero_datapath,  &
      prescribed_strataero_type,      &
      prescribed_strataero_rmfile,    &
      prescribed_strataero_cycle_yr,  &
      prescribed_strataero_fixed_ymd, &
      prescribed_strataero_fixed_tod, &
      prescribed_strataero_use_chemtrop
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_strataero_specifier= specifier
   prescribed_strataero_file     = filename
   prescribed_strataero_filelist = filelist
   prescribed_strataero_datapath = datapath
   prescribed_strataero_type     = data_type
   prescribed_strataero_rmfile   = rmv_file
   prescribed_strataero_cycle_yr = cycle_yr
   prescribed_strataero_fixed_ymd= fixed_ymd
   prescribed_strataero_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_strataero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_strataero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_strataero_specifier,len(prescribed_strataero_specifier)*7, mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_file,     len(prescribed_strataero_file),        mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_filelist, len(prescribed_strataero_filelist),    mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_datapath, len(prescribed_strataero_datapath),    mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_type,     len(prescribed_strataero_type),        mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_strataero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_strataero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_strataero_fixed_tod,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_strataero_use_chemtrop, 1, mpilog,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier(:) = prescribed_strataero_specifier(:)
   filename   = prescribed_strataero_file
   filelist   = prescribed_strataero_filelist
   datapath   = prescribed_strataero_datapath
   data_type  = prescribed_strataero_type
   rmv_file   = prescribed_strataero_rmfile
   cycle_yr   = prescribed_strataero_cycle_yr
   fixed_ymd  = prescribed_strataero_fixed_ymd
   fixed_tod  = prescribed_strataero_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 .and. filename.ne.'NONE') has_prescribed_strataero = .true.

end subroutine prescribed_strataero_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_strataero_register()
    use ppgrid,         only: pver,pcols
    use physics_buffer, only: pbuf_add_field, dtype_r8
    use pio,            only: var_desc_t, file_desc_t, pio_closefile, pio_inq_varid, pio_seterrorhandling, &
                              PIO_INTERNAL_ERROR, PIO_BCAST_ERROR, PIO_NOERR
    use cam_pio_utils,  only: cam_pio_openfile
    use ioFileMod, only : getfil

    type(var_desc_t)  :: varid
    type(file_desc_t) :: file_handle
    character(len=256) :: filepath, filen
    integer :: ierr

    if (has_prescribed_strataero) then

       filepath = trim(datapath)//'/'//trim(filename)

       call getfil( filepath, filen, 0 )
       call cam_pio_openfile( file_handle, filen, 0 )

       call pio_seterrorhandling(file_handle, PIO_BCAST_ERROR)

       ierr = pio_inq_varid( file_handle, 'so4mass_a1', varid )
       three_mode = three_mode .and. (ierr.eq.PIO_NOERR)
       ierr = pio_inq_varid( file_handle, 'so4mass_a2', varid )
       three_mode = three_mode .and. (ierr.eq.PIO_NOERR)
       ierr = pio_inq_varid( file_handle, 'so4mass_a3', varid )
       three_mode = three_mode .and. (ierr.eq.PIO_NOERR)
       ierr = pio_inq_varid( file_handle, 'diamwet_a1', varid )
       three_mode = three_mode .and. (ierr.eq.PIO_NOERR)
       ierr = pio_inq_varid( file_handle, 'diamwet_a2', varid )
       three_mode = three_mode .and. (ierr.eq.PIO_NOERR)
       ierr = pio_inq_varid( file_handle, 'diamwet_a3', varid )
       three_mode = three_mode .and. (ierr.eq.PIO_NOERR)

       call pio_seterrorhandling(file_handle, PIO_INTERNAL_ERROR)

       call pio_closefile( file_handle )

       if (three_mode) then
          call pbuf_add_field(mmr_name1, 'physpkg', dtype_r8,(/pcols,pver/), mmr_ndx1)
          call pbuf_add_field(mmr_name2, 'physpkg', dtype_r8,(/pcols,pver/), mmr_ndx2)
          call pbuf_add_field(mmr_name3, 'physpkg', dtype_r8,(/pcols,pver/), mmr_ndx3)
          call pbuf_add_field(rad_name1, 'physpkg', dtype_r8,(/pcols,pver/), rad_ndx1)
          call pbuf_add_field(rad_name2, 'physpkg', dtype_r8,(/pcols,pver/), rad_ndx2)
          call pbuf_add_field(rad_name3, 'physpkg', dtype_r8,(/pcols,pver/), rad_ndx3)
          call pbuf_add_field(sad_name, 'physpkg', dtype_r8,(/pcols,pver/), sad_ndx)
          specifier(1:7) = (/'VOLC_MMR1:so4mass_a1            ', &
                             'VOLC_MMR2:so4mass_a2            ', &
                             'VOLC_MMR3:so4mass_a3            ', &
                             'VOLC_RAD_GEOM1:diamwet_a1       ', &
                             'VOLC_RAD_GEOM2:diamwet_a2       ', &
                             'VOLC_RAD_GEOM3:diamwet_a3       ', &
                             'VOLC_SAD:SAD_AERO               ' /)
          rad_fld_no = 4
          sad_fld_no = 7
       else
          if (masterproc) then
             write(iulog, *) ' pbuf add mmr_name = '//trim(mmr_name)
          end if
          call pbuf_add_field(mmr_name, 'physpkg', dtype_r8,(/pcols,pver/), mmr_ndx1)
          call pbuf_add_field(rad_name, 'physpkg', dtype_r8,(/pcols,pver/), rad_ndx1)
          call pbuf_add_field(sad_name, 'physpkg', dtype_r8,(/pcols,pver/), sad_ndx)
          specifier(1:3) = (/'VOLC_MMR:H2SO4_mass             ', &
                             'VOLC_RAD_GEOM:rmode             ', &
                             'VOLC_SAD:sad                    ' /)
          rad_fld_no = 2
          sad_fld_no = 3
       endif
    endif

  endsubroutine prescribed_strataero_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_strataero_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, horiz_only
    use error_messages, only: handle_err
    
    if ( has_prescribed_strataero ) then
       if ( masterproc ) then
          write(iulog,*) 'stratospheric aerosol is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    file%geop_alt = .true.

    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

    if (three_mode) then
       call addfld(dens_name1, (/ 'lev' /), 'I','molecules/cm3', 'prescribed volcanic aerosol number density in Mode 1' )
       call addfld(dens_name2, (/ 'lev' /), 'I','molecules/cm3', 'prescribed volcanic aerosol number density in Mode 2' )
       call addfld(dens_name3, (/ 'lev' /), 'I','molecules/cm3', 'prescribed volcanic aerosol number density in Mode 3' )
       call addfld(mmr_name1, (/ 'lev' /), 'I','kg/kg', 'prescribed volcanic aerosol dry mass mixing ratio in Mode 1' )
       call addfld(mmr_name2, (/ 'lev' /), 'I','kg/kg', 'prescribed volcanic aerosol dry mass mixing ratio in Mode 2' )
       call addfld(mmr_name3, (/ 'lev' /), 'I','kg/kg', 'prescribed volcanic aerosol dry mass mixing ratio in Mode 3' )
       call addfld(rad_name1, (/ 'lev' /), 'I','m', 'volcanic aerosol geometric-mode radius in Mode 1' )
       call addfld(rad_name2, (/ 'lev' /), 'I','m', 'volcanic aerosol geometric-mode radius in Mode 2' )
       call addfld(rad_name3, (/ 'lev' /), 'I','m', 'volcanic aerosol geometric-mode radius in Mode 3' )
       call addfld(mass_name1, (/ 'lev' /), 'I','kg/m^2', 'volcanic aerosol vertical mass path in layer in Mode 1' )
       call addfld(mass_name2, (/ 'lev' /), 'I','kg/m^2', 'volcanic aerosol vertical mass path in layer in Mode 2' )
       call addfld(mass_name3, (/ 'lev' /), 'I','kg/m^2', 'volcanic aerosol vertical mass path in layer in Mode 3' )
       call addfld(mass_column_name1, horiz_only, 'I','kg/m^2', 'volcanic aerosol column mass in Mode 1' )
       call addfld(mass_column_name2, horiz_only, 'I','kg/m^2', 'volcanic aerosol column mass in Mode 2' )
       call addfld(mass_column_name3, horiz_only, 'I','kg/m^2', 'volcanic aerosol column mass IN Mode 3' )
    else
       call addfld(dens_name, (/ 'lev' /), 'I','molecules/cm3', 'prescribed volcanic aerosol number density' )
       call addfld(mmr_name,  (/ 'lev' /), 'I','kg/kg', 'prescribed volcanic aerosol dry mass mixing ratio' )
       call addfld(rad_name,  (/ 'lev' /), 'I','m', 'volcanic aerosol geometric-mode radius' )
       call addfld(mass_name, (/ 'lev' /), 'I','kg/m^2', 'volcanic aerosol vertical mass path in layer' )
       call addfld(mass_column_name, horiz_only, 'I','kg/m^2', 'volcanic aerosol column mass' )
    endif
    call addfld(sad_name, (/ 'lev' /), 'I','cm2/cm3', 'stratospheric aerosol surface area density' )

  end subroutine prescribed_strataero_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_strataero_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz, gravit        ! J/K/molecule
    use tropopause,   only : tropopause_findChemTrop

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use physconst,      only : pi

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer :: c,ncol,i,k
    real(r8) :: to_mmr(pcols,pver)
    real(r8), parameter :: molmass = 4._r8/3._r8*98.0_r8 !convert dry mass to wet mass of h2so4 
    real(r8) :: volcmass1(pcols,pver)
    real(r8) :: volcmass2(pcols,pver)
    real(r8) :: volcmass3(pcols,pver)
    real(r8) :: columnmass1(pcols)
    real(r8) :: columnmass2(pcols)
    real(r8) :: columnmass3(pcols)
    integer  :: tropLev(pcols)
    real(r8) :: area_fact, radius_fact

    real(r8), pointer :: mass1(:,:)
    real(r8), pointer :: mass2(:,:)
    real(r8), pointer :: mass3(:,:)
    real(r8), pointer :: area(:,:)
    real(r8), pointer :: radius1(:,:)
    real(r8), pointer :: radius2(:,:)
    real(r8), pointer :: radius3(:,:)

    !WACCM-derived relation between mass concentration and wet aerosol radius in meters
    real(r8),parameter :: radius_conversion = 1.9e-4_r8

    logical :: zero_aerosols
    real(r8), parameter :: rad2deg = 180._r8/pi                ! radians to degrees conversion factor

    if( .not. has_prescribed_strataero ) return

    call advance_trcdata( fields, file, state, pbuf2d )

    ! copy prescribed tracer fields into state svariable with the correct units
    do c = begchunk,endchunk

       pbuf_chnk => pbuf_get_chunk(pbuf2d, c)

       ncol = state(c)%ncol

       select case ( to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))) )
       case ("molecules/cm3air", "molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
          to_mmr(:ncol,:) = (molmass*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
       case ('kg/kg','mmr','kg kg-1')
          to_mmr(:ncol,:) = 1._r8 ! input file must have converted to wet sulfate mass (=4/3*dry mass)
       case ('mol/mol','mole/mole','vmr','fraction')
          to_mmr(:ncol,:) = molmass/mwdry
       case default
          write(iulog,*) 'prescribed_strataero_adv: mass units = ',trim(fields(1)%units) ,' are not recognized'
          call endrun('prescribed_strataero_adv: mass units are not recognized')
       end select

       if (mmr_ndx1>0) call pbuf_get_field(pbuf_chnk, mmr_ndx1, mass1)
       if (mmr_ndx2>0) call pbuf_get_field(pbuf_chnk, mmr_ndx2, mass2)
       if (mmr_ndx3>0) call pbuf_get_field(pbuf_chnk, mmr_ndx3, mass3)

       if (three_mode) then
          call outfld( dens_name1, mass1(:,:), pcols, state(c)%lchnk)
          call outfld( dens_name2, mass2(:,:), pcols, state(c)%lchnk)
          call outfld( dens_name3, mass3(:,:), pcols, state(c)%lchnk)
       else
          call outfld( dens_name, mass1(:,:), pcols, state(c)%lchnk)
       endif

       if (mmr_ndx1>0) mass1(:ncol,:) = to_mmr(:ncol,:) * mass1(:ncol,:) ! mmr
       if (mmr_ndx2>0) mass2(:ncol,:) = to_mmr(:ncol,:) * mass2(:ncol,:) ! mmr
       if (mmr_ndx3>0) mass3(:ncol,:) = to_mmr(:ncol,:) * mass3(:ncol,:) ! mmr

       if (rad_ndx1>0) call pbuf_get_field(pbuf_chnk, rad_ndx1, radius1)
       if (rad_ndx2>0) call pbuf_get_field(pbuf_chnk, rad_ndx2, radius2)
       if (rad_ndx3>0) call pbuf_get_field(pbuf_chnk, rad_ndx3, radius3)

       select case ( to_lower(trim(fields(rad_fld_no)%units(:GLC(fields(rad_fld_no)%units)))) )
       case ("m","meters")
          radius_fact = 1._r8
       case ("cm","centimeters")
          radius_fact = 1.e-2_r8
       case default
          write(iulog,*) 'prescribed_strataero_adv: radius units = ',trim(fields(rad_fld_no)%units) ,' are not recognized'
          call endrun('prescribed_strataero_adv: radius units are not recognized')
       end select

       !MAM output is diamter so we need to half the value
       if (three_mode) then
          radius1(:ncol,:) = radius_fact*radius1(:ncol,:)*0.5_r8
          radius2(:ncol,:) = radius_fact*radius2(:ncol,:)*0.5_r8
          radius3(:ncol,:) = radius_fact*radius3(:ncol,:)*0.5_r8
       else
          radius1(:ncol,:) = radius_fact*radius1(:ncol,:)
       endif

       call pbuf_get_field(pbuf_chnk, sad_ndx, area)

       select case ( to_lower(trim(fields(sad_fld_no)%units(:7))) )
       case ("um2/cm3")
          area_fact = 1.e-8_r8
       case ("cm2/cm3")
          area_fact = 1._r8
       case default
          write(iulog,*) 'prescribed_strataero_adv: surface area density units = ',&
                         trim(fields(rad_fld_no)%units) ,' are not recognized'
          call endrun('prescribed_strataero_adv: surface area density units are not recognized')
       end select
       area(:ncol,:) = area_fact*area(:ncol,:)

       ! this definition of tropopause is consistent with what is used in chemistry
       call tropopause_findChemTrop(state(c), tropLev)

       do i = 1,ncol
          do k = 1,pver
             zero_aerosols = k >= tropLev(i)
             if ( .not.prescribed_strataero_use_chemtrop .and. abs( state(c)%lat(i)*rad2deg ) > 50._r8 ) then
                zero_aerosols = state(c)%pmid(i,k) >= 30000._r8
             endif
             ! set to zero at and below tropopause
             if ( zero_aerosols ) then
                if (mmr_ndx1>0) mass1(i,k) = 0._r8
                if (mmr_ndx2>0) mass2(i,k) = 0._r8
                if (mmr_ndx3>0) mass3(i,k) = 0._r8
                if (rad_ndx1>0) radius1(i,k) = 0._r8
                if (rad_ndx2>0) radius2(i,k) = 0._r8
                if (rad_ndx3>0) radius3(i,k) = 0._r8
                area(i,k) = 0._r8
             endif
          enddo
       enddo

       volcmass1(:ncol,:) = mass1(:ncol,:)*state(c)%pdel(:ncol,:)/gravit
       columnmass1(:ncol) = sum(volcmass1(:ncol,:), 2)

       if (three_mode) then
          volcmass2(:ncol,:) = mass2(:ncol,:)*state(c)%pdel(:ncol,:)/gravit
          volcmass3(:ncol,:) = mass3(:ncol,:)*state(c)%pdel(:ncol,:)/gravit
          columnmass2(:ncol) = sum(volcmass2(:ncol,:), 2)
          columnmass3(:ncol) = sum(volcmass3(:ncol,:), 2)
          call outfld( mmr_name1,         mass1(:,:),     pcols, state(c)%lchnk)
          call outfld( mmr_name2,         mass2(:,:),     pcols, state(c)%lchnk)
          call outfld( mmr_name3,         mass3(:,:),     pcols, state(c)%lchnk)
          call outfld( mass_name1,        volcmass1(:,:), pcols, state(c)%lchnk)
          call outfld( mass_name2,        volcmass2(:,:), pcols, state(c)%lchnk)
          call outfld( mass_name3,        volcmass3(:,:), pcols, state(c)%lchnk)
          call outfld( mass_column_name1, columnmass1(:), pcols, state(c)%lchnk)
          call outfld( mass_column_name2, columnmass2(:), pcols, state(c)%lchnk)
          call outfld( mass_column_name3, columnmass3(:), pcols, state(c)%lchnk)
          call outfld( rad_name1,         radius1(:,:),   pcols, state(c)%lchnk)
          call outfld( rad_name2,         radius2(:,:),   pcols, state(c)%lchnk)
          call outfld( rad_name3,         radius3(:,:),   pcols, state(c)%lchnk)
       else
          call outfld( mmr_name,         mass1(:,:),     pcols, state(c)%lchnk)
          call outfld( mass_name,        volcmass1(:,:), pcols, state(c)%lchnk)
          call outfld( mass_column_name, columnmass1(:), pcols, state(c)%lchnk)
          call outfld( rad_name,         radius1(:,:),   pcols, state(c)%lchnk)
       endif

       call outfld( sad_name,         area(:,:),     pcols, state(c)%lchnk)

    enddo

  end subroutine prescribed_strataero_adv

!-------------------------------------------------------------------
  subroutine init_prescribed_strataero_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart

    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_strataero', piofile, file )

  end subroutine init_prescribed_strataero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_strataero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_strataero_restart
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_strataero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_strataero', piofile, file )

  end subroutine read_prescribed_strataero_restart

end module prescribed_strataero
