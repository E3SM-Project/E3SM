
module waccm_forcing
!================================================================================================
!
! Provides WACCM forcing data for use without interactive chemistry
! -- for GHG chemistry
!
! FVITT 21 Mar 2011  -- creation
!
!================================================================================================

  use shr_kind_mod,       only: r8 => shr_kind_r8
  use cam_abortutils,         only: endrun
  use cam_logfile,        only: iulog

  use tracer_data,  only : trfld, trfile
  use ppgrid,       only : pcols, pver
  use ppgrid,       only : begchunk, endchunk
  use spmd_utils,   only : masterproc

  implicit none
  private
  save

! Public interfaces
  public :: waccm_forcing_init
  public :: waccm_forcing_adv
  public :: get_cnst               ! return prescribed constituents for nlte
  public :: get_solar              ! return prescribed net solar heating rate
  public :: waccm_forcing_readnl

! Private module data

  type(trfld), pointer :: fields(:) => null()
  type(trfile)         :: file

  integer, parameter :: N_FLDS = 7
  integer, parameter :: N_MMRS = 6

  integer :: number_flds

  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'CYCLICAL'
  logical            :: rmv_file = .false.

  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  integer ::         o1_ndx=1, o2_ndx=2, o3_ndx=3, no_ndx=4, h_ndx=5, co2_ndx=6, qrs_ndx=7
  character(len=16)  :: specifier(N_FLDS) = (/ 'O      ','O2     ','O3     ', 'NO     ', 'H      ','CO2    ', 'QRS_TOT' /)
  real(r8), parameter :: molmass(N_MMRS) = (/ 15.99940_r8, 31.99880_r8, 47.99820_r8, 30.00614_r8, 1.007400_r8, 44.00980_r8 /)

!================================================================================================
contains
!================================================================================================

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine waccm_forcing_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'prescribed_aero_readnl'

    character(len=16)  :: waccm_forcing_specifier(N_FLDS)
    character(len=256) :: waccm_forcing_file
    character(len=256) :: waccm_forcing_filelist
    character(len=256) :: waccm_forcing_datapath
    character(len=32)  :: waccm_forcing_type
    logical            :: waccm_forcing_rmfile
    integer            :: waccm_forcing_cycle_yr
    integer            :: waccm_forcing_fixed_ymd
    integer            :: waccm_forcing_fixed_tod

    namelist /waccm_forcing_nl/ &
         waccm_forcing_specifier, &
         waccm_forcing_file,      &
         waccm_forcing_filelist,  &
         waccm_forcing_datapath,  &
         waccm_forcing_type,      &
         waccm_forcing_rmfile,    &
         waccm_forcing_cycle_yr,  &
         waccm_forcing_fixed_ymd, &
         waccm_forcing_fixed_tod      
    !-----------------------------------------------------------------------------

    ! Initialize namelist variables from local module variables.
    waccm_forcing_specifier= specifier
    waccm_forcing_file     = filename
    waccm_forcing_filelist = filelist
    waccm_forcing_datapath = datapath
    waccm_forcing_type     = datatype
    waccm_forcing_rmfile   = rmv_file
    waccm_forcing_cycle_yr = cycle_yr
    waccm_forcing_fixed_ymd= fixed_ymd
    waccm_forcing_fixed_tod= fixed_tod

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'waccm_forcing_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, waccm_forcing_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(waccm_forcing_specifier,len(waccm_forcing_specifier(1))*N_FLDS, mpichar, 0, mpicom)
    call mpibcast(waccm_forcing_file,     len(waccm_forcing_file),     mpichar, 0, mpicom)
    call mpibcast(waccm_forcing_filelist, len(waccm_forcing_filelist), mpichar, 0, mpicom)
    call mpibcast(waccm_forcing_datapath, len(waccm_forcing_datapath), mpichar, 0, mpicom)
    call mpibcast(waccm_forcing_type,     len(waccm_forcing_type),     mpichar, 0, mpicom)
    call mpibcast(waccm_forcing_rmfile,   1, mpilog,  0, mpicom)
    call mpibcast(waccm_forcing_cycle_yr, 1, mpiint,  0, mpicom)
    call mpibcast(waccm_forcing_fixed_ymd,1, mpiint,  0, mpicom)
    call mpibcast(waccm_forcing_fixed_tod,1, mpiint,  0, mpicom)
#endif

    ! Update module variables with user settings.
    specifier  = waccm_forcing_specifier
    filename   = waccm_forcing_file
    filelist   = waccm_forcing_filelist
    datapath   = waccm_forcing_datapath
    datatype   = waccm_forcing_type
    rmv_file   = waccm_forcing_rmfile
    cycle_yr   = waccm_forcing_cycle_yr
    fixed_ymd  = waccm_forcing_fixed_ymd
    fixed_tod  = waccm_forcing_fixed_tod

  end subroutine waccm_forcing_readnl

 !------------------------------------------------------------------------
 !------------------------------------------------------------------------
  subroutine waccm_forcing_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: i


    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

    do i = 1,N_FLDS
       call addfld( 'WFRC_'//trim(fields(i)%fldnam), (/ 'lev' /), 'I', fields(i)%units, 'for waccm forcing' )
    enddo

    return
  end subroutine waccm_forcing_init

!=======================================================================

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine waccm_forcing_adv (state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry           ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz           ! J/K/molecule
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    integer :: c,ncol,i
    real(r8) :: to_mmr(pcols,pver)

    call advance_trcdata( fields, file, state, pbuf2d  )
    
    ! set the tracer fields with the correct units
    do i = 1,N_FLDS

       do c = begchunk,endchunk
          ncol = state(c)%ncol

          if ( i<=N_MMRS ) then

             select case ( to_lower(trim(fields(i)%units(:GLC(fields(i)%units)))) )
             case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
                to_mmr(:ncol,:) = (molmass(i)*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
             case ('kg/kg','mmr')
                to_mmr(:ncol,:) = 1._r8
             case ('mol/mol','mole/mole','vmr','fraction')
                to_mmr(:ncol,:) = molmass(i)/mwdry
             case default
                print*, 'waccm_forcing_adv: units = ',trim(fields(i)%units) ,' are not recognized'
                call endrun('waccm_forcing_adv: units are not recognized')
             end select

             call outfld( 'WFRC_'//trim(fields(i)%fldnam), fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk )
             
             fields(i)%data(:ncol,:,c) = to_mmr(:ncol,:) * fields(i)%data(:ncol,:,c)
          else
             call outfld( 'WFRC_'//trim(fields(i)%fldnam), fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk )
          endif

       enddo
    enddo

    return
  end subroutine waccm_forcing_adv


!================================================================================================

  subroutine get_cnst (ncol, lchnk, co2, o1, o2, no, h, o3)
!
! Get mass mixing ratios specified from input dataset for used in Fomichev routines
!-------------------------------------------------------------------------

! Arguments
    integer,  intent(in)  :: ncol                   ! no. of columns in chunk
    integer,  intent(in)  :: lchnk                  ! chunk identifier

    real(r8), optional, pointer, dimension(:,:) :: co2
    real(r8), optional, pointer, dimension(:,:) :: o1
    real(r8), optional, pointer, dimension(:,:) :: o2
    real(r8), optional, pointer, dimension(:,:) :: no
    real(r8), optional, pointer, dimension(:,:) :: h
    real(r8), optional, pointer, dimension(:,:) :: o3

!------------------------------------------------------------------------

    if (present(co2)) then
       co2 => fields(co2_ndx)%data(:,:,lchnk)
    endif
    if (present(o1)) then
       o1  => fields(o1_ndx )%data(:,:,lchnk)
    endif
    if (present(o2)) then
       o2  => fields(o2_ndx )%data(:,:,lchnk)
    endif
    if (present(no)) then
       no  => fields(no_ndx )%data(:,:,lchnk)
    endif
    if (present(h)) then
       h  => fields(no_ndx )%data(:,:,lchnk)
    endif
    if (present(o3)) then
       o3 => fields(o3_ndx)%data(:,:,lchnk)
    endif

  end subroutine get_cnst

!================================================================================================

  subroutine get_solar (ncol, lchnk, qrs_mlt)
!
! Get M/LT solar heating rates specified from input dataset
!-------------------------------------------------------------------------

! Arguments
    integer,  intent(in)  :: ncol                   ! no. of columns in chunk
    integer,  intent(in)  :: lchnk                  ! chunk identifier
    real(r8), intent(out) :: qrs_mlt(:,:)    ! M/LT solar heating rates

! Local workspace
    integer :: k

!------------------------------------------------------------------------

   qrs_mlt(:ncol,:) = fields(qrs_ndx)%data(:ncol,:,lchnk)

  end subroutine get_solar

end module waccm_forcing
