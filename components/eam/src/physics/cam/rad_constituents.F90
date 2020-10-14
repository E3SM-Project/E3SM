
module rad_constituents

!------------------------------------------------------------------------------------------------
!
! Provide constituent distributions and properties to the radiation and 
! cloud microphysics routines.
! 
! The logic to control which constituents are used in the climate calculations
! and which are used in diagnostic radiation calculations is contained in this module.
!
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: rga
use physics_types,  only: physics_state
use constituents,   only: cnst_name, cnst_get_ind
use radconstants,   only: gasnamelength, nradgas, rad_gas_index, ot_length
use phys_prop,      only: physprop_accum_unique_files, physprop_init, &
                          physprop_get_id, physprop_get
use cam_history,    only: addfld, horiz_only, fieldname_len, add_default, outfld
use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index


use error_messages, only: alloc_err   
use cam_abortutils,     only: endrun
use cam_logfile,    only: iulog

implicit none
private
save

! Public interfaces

public :: &
   rad_cnst_readnl,             &! read namelist values and parse
   rad_cnst_init,               &! find optics files and all constituents
   rad_cnst_get_info,           &! return info about climate/diagnostic lists
   rad_cnst_get_mode_idx,       &! return mode index of specified mode type
   rad_cnst_get_spec_idx,       &! return specie index of specified specie type
   rad_cnst_get_gas,            &! return pointer to mmr for gasses
   rad_cnst_get_aer_mmr,        &! return pointer to mmr for aerosols
   rad_cnst_get_mam_mmr_idx,    &! get constituent index of mam specie mmr (climate list only)
   rad_cnst_get_aer_props,      &! return physical properties for aerosols
   rad_cnst_get_mode_props,     &! return physical properties for aerosol modes
   rad_cnst_get_mode_num,       &! return mode number mixing ratio
   rad_cnst_get_mode_num_idx,   &! get constituent index of mode number m.r. (climate list only)
   rad_cnst_out,                &! output constituent diagnostics (mass per layer and column burden)
   rad_cnst_get_call_list        ! return list of active climate/diagnostic calls to radiation

integer, parameter :: cs1 = 256
integer, public, parameter :: N_DIAG = 10
character(len=cs1), public :: iceopticsfile, liqopticsfile
character(len=32),  public :: icecldoptics,liqcldoptics
logical,            public :: oldcldoptics = .false.

! Private module data

! max number of strings in mode definitions
integer, parameter :: n_mode_str = 100    ! max number of strings in mode definitions

! max number of externally mixed entities in the climate/diag lists
integer, parameter :: n_rad_cnst = N_RAD_CNST

! Namelist variables
character(len=cs1), dimension(n_mode_str) :: mode_defs   = ' '
character(len=cs1) :: rad_climate(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_1(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_2(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_3(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_4(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_5(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_6(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_7(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_8(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_9(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_10(n_rad_cnst) = ' '

! type to provide access to the components of a mode
type :: mode_component_t
   integer :: nspec
   ! For "source" variables below, value is:
   ! 'N' if in pbuf (non-advected)
   ! 'A' if in state (advected)
   character(len=  1) :: source_num_a  ! source of interstitial number conc field
   character(len= 32) :: camname_num_a ! name registered in pbuf or constituents for number mixing ratio of interstitial species
   character(len=  1) :: source_num_c  ! source of cloud borne number conc field
   character(len= 32) :: camname_num_c ! name registered in pbuf or constituents for number mixing ratio of cloud borne species
   character(len=  1), pointer :: source_mmr_a(:)  ! source of interstitial specie mmr fields
   character(len= 32), pointer :: camname_mmr_a(:) ! name registered in pbuf or constituents for mmr of interstitial components
   character(len=  1), pointer :: source_mmr_c(:)  ! source of cloud borne specie mmr fields
   character(len= 32), pointer :: camname_mmr_c(:) ! name registered in pbuf or constituents for mmr of cloud borne components
   character(len= 32), pointer :: type(:)          ! specie type (as used in MAM code)
   character(len=cs1), pointer :: props(:)         ! file containing specie properties
   integer          :: idx_num_a    ! index in pbuf or constituents for number mixing ratio of interstitial species
   integer          :: idx_num_c    ! index in pbuf for number mixing ratio of interstitial species
   integer, pointer :: idx_mmr_a(:) ! index in pbuf or constituents for mmr of interstitial species
   integer, pointer :: idx_mmr_c(:) ! index in pbuf for mmr of interstitial species
   integer, pointer :: idx_props(:) ! ID used to access physical properties of mode species from phys_prop module
end type mode_component_t

! type to provide access to all modes
type :: modes_t
   integer :: nmodes
   character(len= 32),     pointer :: names(:) ! names used to identify a mode in the climate/diag lists
   character(len= 32),     pointer :: types(:) ! type of mode (as used in MAM code)
   type(mode_component_t), pointer :: comps(:) ! components which define the mode
end type modes_t

type(modes_t), target :: modes  ! mode definitions

! type to provide access to the data parsed from the rad_climate and rad_diag_* strings
type :: rad_cnst_namelist_t
   integer :: ncnst
   character(len=  1), pointer :: source(:)  ! 'A' for state (advected), 'N' for pbuf (non-advected), 
                                             ! 'M' for mode, 'Z' for zero
   character(len= 64), pointer :: camname(:) ! name registered in pbuf or constituents
   character(len=cs1), pointer :: radname(:) ! radname is the name as identfied in radiation,
                                             ! must be one of (rgaslist if a gas) or
                                             ! (/fullpath/filename.nc if an aerosol)
   character(len=  1), pointer :: type(:)    ! 'A' if aerosol, 'G' if gas, 'M' if mode
end type rad_cnst_namelist_t

type(rad_cnst_namelist_t) :: namelist(0:N_DIAG) ! gas, bulk aerosol, and modal components used in
                                                ! climate/diagnostic calculations

logical :: active_calls(0:N_DIAG)     ! active_calls(i) is true if the i-th call to radiation is 
                                      ! specified.  Note that the 0th call is for the climate
                                      ! calculation which is always made.

! Storage for gas components in the climate/diagnostic lists

type :: gas_t
   character(len=1)  :: source       ! A for state (advected), N for pbuf (non-advected), Z for zero
   character(len=64) :: camname      ! name of constituent in physics state or buffer
   character(len=32) :: mass_name    ! name for mass per layer field in history output
   integer           :: idx          ! index from constituents or from pbuf
end type gas_t

type :: gaslist_t
   integer                :: ngas
   character(len=2)       :: list_id  ! set to "  " for climate list, or two character integer
                                      ! (include leading zero) to identify diagnostic list
   type(gas_t), pointer   :: gas(:)   ! dimension(ngas) where ngas = nradgas is from radconstants
end type gaslist_t

type(gaslist_t), target :: gaslist(0:N_DIAG)  ! gasses used in climate/diagnostic calculations

! Storage for bulk aerosol components in the climate/diagnostic lists

type :: aerosol_t
   character(len=1)   :: source         ! A for state (advected), N for pbuf (non-advected), Z for zero
   character(len=64)  :: camname        ! name of constituent in physics state or buffer
   character(len=cs1) :: physprop_file  ! physprop filename
   character(len=32)  :: mass_name      ! name for mass per layer field in history output
   integer            :: idx            ! index of constituent in physics state or buffer
   integer            :: physprop_id    ! ID used to access physical properties from phys_prop module
end type aerosol_t

type :: aerlist_t
   integer                  :: numaerosols  ! number of aerosols
   character(len=2)         :: list_id      ! set to "  " for climate list, or two character integer
                                            ! (include leading zero) to identify diagnostic list
   type(aerosol_t), pointer :: aer(:)       ! dimension(numaerosols)
end type aerlist_t

type(aerlist_t), target :: aerosollist(0:N_DIAG) ! list of aerosols used in climate/diagnostic calcs

! storage for modal aerosol components in the climate/diagnostic lists

type :: modelist_t
   integer          :: nmodes              ! number of modes
   character(len=2) :: list_id             ! set to "  " for climate list, or two character integer
                                           ! (include leading zero) to identify diagnostic list
   integer,   pointer :: idx(:)            ! index of the mode in the mode definition object
   character(len=cs1), pointer :: physprop_files(:) ! physprop filename
   integer,   pointer :: idx_props(:)      ! index of the mode properties in the physprop object
end type modelist_t

type(modelist_t), target :: ma_list(0:N_DIAG) ! list of aerosol modes used in climate/diagnostic calcs


! values for constituents with requested value of zero
real(r8), allocatable, target :: zero_cols(:,:) 

! define generic interface routines
interface rad_cnst_get_info
   module procedure rad_cnst_get_info
   module procedure rad_cnst_get_info_by_mode
   module procedure rad_cnst_get_info_by_mode_spec
   module procedure rad_cnst_get_info_by_spectype
end interface

interface rad_cnst_get_aer_mmr
   module procedure rad_cnst_get_aer_mmr_by_idx
   module procedure rad_cnst_get_mam_mmr_by_idx
end interface

interface rad_cnst_get_aer_props
   module procedure rad_cnst_get_aer_props_by_idx
   module procedure rad_cnst_get_mam_props_by_idx
end interface

logical :: verbose = .true.
character(len=1), parameter :: nl = achar(10)

#if ( defined MODAL_AERO_9MODE )
integer, parameter :: num_mode_types = 10
integer, parameter :: num_spec_types = 11
character(len=14), parameter :: mode_type_names(num_mode_types) = (/ &
   'accum         ', 'aitken        ', 'primary_carbon', 'fine_seasalt  ', &
   'fine_dust     ', 'coarse        ', 'coarse_seasalt', 'coarse_dust   ', &
   'accum_marine  ', 'aitken_marine ' /)
character(len=9), parameter :: spec_type_names(num_spec_types) = (/ &
   'sulfate  ', 'ammonium ', 'nitrate  ', 'p-organic', &
   's-organic', 'black-c  ', 'seasalt  ', 'dust     ', &
   'm-poly   ', 'm-prot   ', 'm-lip    ' /)
#elif ( defined MODAL_AERO_4MODE_MOM )
integer, parameter :: num_mode_types = 8
integer, parameter :: num_spec_types = 9
character(len=14), parameter :: mode_type_names(num_mode_types) = (/ &
   'accum         ', 'aitken        ', 'primary_carbon', 'fine_seasalt  ', &
   'fine_dust     ', 'coarse        ', 'coarse_seasalt', 'coarse_dust   '  /)
character(len=9), parameter :: spec_type_names(num_spec_types) = (/ &
   'sulfate  ', 'ammonium ', 'nitrate  ', 'p-organic', &
   's-organic', 'black-c  ', 'seasalt  ', 'dust     ', &
   'm-organic' /)
#else
integer, parameter :: num_mode_types = 8
integer, parameter :: num_spec_types = 8
character(len=14), parameter :: mode_type_names(num_mode_types) = (/ &
   'accum         ', 'aitken        ', 'primary_carbon', 'fine_seasalt  ', &
   'fine_dust     ', 'coarse        ', 'coarse_seasalt', 'coarse_dust   '  /)
character(len=9), parameter :: spec_type_names(num_spec_types) = (/ &
   'sulfate  ', 'ammonium ', 'nitrate  ', 'p-organic', &
   's-organic', 'black-c  ', 'seasalt  ', 'dust     '/)
#endif


!==============================================================================
contains
!==============================================================================

subroutine rad_cnst_readnl(nlfile)

   ! Read rad_cnst_nl namelist group.  Parse input.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=2) :: suffix
   character(len=1), pointer   :: ctype(:)
   character(len=*), parameter :: subname = 'rad_cnst_readnl'

   namelist /rad_cnst_nl/ mode_defs,     &
                          rad_climate,   &
                          rad_diag_1,    &
                          rad_diag_2,    &
                          rad_diag_3,    &
                          rad_diag_4,    &
                          rad_diag_5,    &
                          rad_diag_6,    &
                          rad_diag_7,    &
                          rad_diag_8,    &
                          rad_diag_9,    &
                          rad_diag_10,   &
                          iceopticsfile, &
                          liqopticsfile, &
                          icecldoptics,  &
                          liqcldoptics,  &
                          oldcldoptics

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rad_cnst_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rad_cnst_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast (mode_defs,     len(mode_defs(1))*n_mode_str,     mpichar, 0, mpicom)
   call mpibcast (rad_climate,   len(rad_climate(1))*n_rad_cnst,   mpichar, 0, mpicom)
   call mpibcast (rad_diag_1,    len(rad_diag_1(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_2,    len(rad_diag_2(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_3,    len(rad_diag_3(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_4,    len(rad_diag_4(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_5,    len(rad_diag_5(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_6,    len(rad_diag_6(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_7,    len(rad_diag_7(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_8,    len(rad_diag_8(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_9,    len(rad_diag_9(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_10,   len(rad_diag_10(1))*n_rad_cnst,   mpichar, 0, mpicom)
   call mpibcast (iceopticsfile, len(iceopticsfile),               mpichar, 0, mpicom)
   call mpibcast (liqopticsfile, len(liqopticsfile),               mpichar, 0, mpicom)
   call mpibcast (liqcldoptics,  len(liqcldoptics),                mpichar, 0, mpicom)
   call mpibcast (icecldoptics,  len(icecldoptics),                mpichar, 0, mpicom)
   call mpibcast (oldcldoptics,  1,                                mpilog , 0, mpicom)
#endif

   ! Parse the namelist input strings

   ! Mode definition stings
   call parse_mode_defs(mode_defs, modes)
   
   ! Lists of externally mixed entities for climate and diagnostic calculations
   do i = 0,N_DIAG
      select case (i)
      case(0)
         call parse_rad_specifier(rad_climate, namelist(i))
      case (1)
         call parse_rad_specifier(rad_diag_1, namelist(i))
      case (2)
         call parse_rad_specifier(rad_diag_2, namelist(i))
      case (3)
         call parse_rad_specifier(rad_diag_3, namelist(i))
      case (4)
         call parse_rad_specifier(rad_diag_4, namelist(i))
      case (5)
         call parse_rad_specifier(rad_diag_5, namelist(i))
      case (6)
         call parse_rad_specifier(rad_diag_6, namelist(i))
      case (7)
         call parse_rad_specifier(rad_diag_7, namelist(i))
      case (8)
         call parse_rad_specifier(rad_diag_8, namelist(i))
      case (9)
         call parse_rad_specifier(rad_diag_9, namelist(i))
      case (10)
         call parse_rad_specifier(rad_diag_10, namelist(i))
      end select
   enddo

   ! were there any constituents specified for the nth diagnostic call?
   ! if so, radiation will make a call with those consituents
   active_calls(:) = (namelist(:)%ncnst > 0)
   	
   ! Initialize the gas and aerosol lists with the information from the
   ! namelist.  This is done here so that this information is available via
   ! the query functions at the time when the register methods are called.

   ! Set the list_id fields which distinquish the climate and diagnostic lists
   do i = 0, N_DIAG
      if (active_calls(i)) then
         if (i > 0) then
            write(suffix, fmt = '(i2.2)') i
         else
            suffix='  '
         end if
         aerosollist(i)%list_id = suffix
         gaslist(i)%list_id     = suffix
         ma_list(i)%list_id     = suffix
      end if
   end do

   ! Create a list of the unique set of filenames containing property data

   ! Start with the bulk aerosol species in the climate/diagnostic lists.
   ! The physprop_accum_unique_files routine has the side effect of returning the number
   ! of bulk aerosols in each list (they're identified by type='A').
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call physprop_accum_unique_files(namelist(i)%radname, namelist(i)%type)
      endif
   enddo

   ! Add physprop files for the species from the mode definitions.
   do i = 1, modes%nmodes
      allocate(ctype(modes%comps(i)%nspec))
      ctype = 'A'
      call physprop_accum_unique_files(modes%comps(i)%props, ctype)
      deallocate(ctype)
   end do

   ! Initialize the gas, bulk aerosol, and modal aerosol lists.  This step splits the
   ! input climate/diagnostic lists into the corresponding gas, bulk and modal aerosol
   ! lists.
   if (masterproc) write(iulog,*) nl//subname//': Radiation constituent lists:'
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call list_init1(namelist(i), gaslist(i), aerosollist(i), ma_list(i))

         if (masterproc .and. verbose) then
            call print_lists(gaslist(i), aerosollist(i), ma_list(i))
         end if

      end if
   end do

   if (masterproc .and. verbose) call print_modes(modes)

end subroutine rad_cnst_readnl

!================================================================================================

subroutine rad_cnst_init()

   ! The initialization of the gas and aerosol lists is finished by
   ! 1) read the physprop files
   ! 2) find the index of each constituent in the constituent or physics buffer arrays
   ! 3) find the index of the aerosol constituents used to access its properties from the
   !    physprop module.

   integer :: i
   integer :: num_aerosols
   logical, parameter :: stricttest = .true.
   character(len=*), parameter :: subname = 'rad_cnst_init'
   !-----------------------------------------------------------------------------

   ! memory to point to if zero value requested
   allocate(zero_cols(pcols,pver))
   zero_cols = 0._r8

   ! Allocate storage for the physical properties of each aerosol; read properties from
   ! the data files.
   call physprop_init()

   ! Start checking that specified radiative constituents are present in the constituent
   ! or physics buffer arrays.
   if (masterproc) write(iulog,*) nl//subname//': checking for radiative constituents'

   ! Finish initializing the mode definitions.
   call init_mode_comps(modes)

   ! Finish initializing the gas, bulk aerosol, and mode lists.
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call list_init2(gaslist(i), aerosollist(i), ma_list(i))
      end if
   end do

   ! Check that all gases supported by the radiative transfer code have been specified.
   if (stricttest) then
      do i = 1, nradgas
         if (gaslist(0)%gas(i)%source .eq. 'Z' ) then
            call endrun(subname//': list of radiative gasses must include all radiation gasses for the climate specication')
         endif
      enddo
   endif

   ! Initialize history output of climate diagnostic quantities
   call rad_gas_diag_init(gaslist(0))
   call rad_aer_diag_init(aerosollist(0))


end subroutine rad_cnst_init

!================================================================================================

subroutine rad_cnst_get_gas(list_idx, gasname, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the gas from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),            intent(in) :: gasname
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: lchnk
   integer :: igas
   integer :: idx
   character(len=1) :: source
   type(gaslist_t), pointer :: list
   character(len=*), parameter :: subname = 'rad_cnst_get_gas'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      list => gaslist(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif
      
   lchnk = state%lchnk

   ! Get index of gas in internal arrays.  rad_gas_index will abort if the 
   ! specified gasname is not recognized by the radiative transfer code.
   igas = rad_gas_index(trim(gasname))
 
   ! Get data source
   source = list%gas(igas)%source
   idx    = list%gas(igas)%idx
   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_gas

!================================================================================================

subroutine rad_cnst_get_info(list_idx, gasnames, aernames, &
                             use_data_o3, ngas, naero, nmodes)

   ! Return info about gas and aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=64), optional, intent(out) :: gasnames(:)
   character(len=64), optional, intent(out) :: aernames(:)
   logical,           optional, intent(out) :: use_data_o3
   integer,           optional, intent(out) :: naero
   integer,           optional, intent(out) :: ngas
   integer,           optional, intent(out) :: nmodes

   ! Local variables
   type(gaslist_t),  pointer :: g_list ! local pointer to gas list of interest
   type(aerlist_t),  pointer :: a_list ! local pointer to aerosol list of interest
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer          :: i
   integer          :: arrlen  ! length of assumed shape array
   integer          :: gaslen  ! length of assumed shape array
   integer          :: igas    ! index of a gas in the gas list
   character(len=1) :: source  ! A for state, N for pbuf, Z for zero

   character(len=*), parameter :: subname = 'rad_cnst_get_info'
   !-----------------------------------------------------------------------------

   g_list => gaslist(list_idx)
   a_list => aerosollist(list_idx)
   m_list => ma_list(list_idx)

   ! number of bulk aerosols in list
   if (present(naero)) then
      naero = a_list%numaerosols
   endif

   ! number of aerosol modes in list
   if (present(nmodes)) then
      nmodes = m_list%nmodes
   endif

   ! number of gases in list
   if (present(ngas)) then
      ngas = g_list%ngas
   endif

   ! names of aerosols in list
   if (present(aernames)) then

      ! check that output array is long enough
      arrlen = size(aernames)
      if (arrlen < a_list%numaerosols) then
         write(iulog,*) subname//': ERROR: naero=', a_list%numaerosols, '  arrlen=', arrlen
         call endrun(subname//': ERROR: aernames too short')
      end if

      do i = 1, a_list%numaerosols
         aernames(i) = a_list%aer(i)%camname
      end do

   end if

   ! names of gas in list
   if (present(gasnames)) then

      ! check that output array is long enough
      gaslen = size(gasnames)
      if (gaslen < g_list%ngas) then
         write(iulog,*) subname//': ERROR: ngas=', g_list%ngas, '  gaslen=', gaslen
         call endrun(subname//': ERROR: gasnames too short')
      end if

      do i = 1, g_list%ngas
         gasnames(i) = g_list%gas(i)%camname
      end do

   end if

   ! Does the climate calculation use data ozone?
   if (present(use_data_o3)) then

      ! get index of O3 in gas list
      igas = rad_gas_index('O3')
 
      ! Get data source
      source = g_list%gas(igas)%source

      use_data_o3 = .false.
      if (source == 'N') use_data_o3 = .true.
   endif

end subroutine rad_cnst_get_info

!================================================================================================

subroutine rad_cnst_get_info_by_mode(list_idx, m_idx, &
   mode_type, num_name, num_name_cw, nspec)

   ! Return info about modal aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of mode in the specified list
   character(len=32), optional, intent(out) :: mode_type   ! type of mode (as used in MAM code)
   character(len=32), optional, intent(out) :: num_name    ! name of interstitial number mixing ratio
   character(len=32), optional, intent(out) :: num_name_cw ! name of cloud borne number mixing ratio
   integer,           optional, intent(out) :: nspec       ! number of species in the mode

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer          :: nmodes
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_mode'
   !-----------------------------------------------------------------------------

   m_list => ma_list(list_idx)

   ! check for valid mode index
   nmodes = m_list%nmodes
   if (m_idx < 1 .or. m_idx > nmodes) then
      write(iulog,*) subname//': ERROR - invalid mode index: ', m_idx
      call endrun(subname//': ERROR - invalid mode index')
   end if

   ! get index into the mode definition object
   mm = m_list%idx(m_idx)

   ! mode type
   if (present(mode_type)) then
      mode_type = modes%types(mm)
   endif

   ! number of species in the mode
   if (present(nspec)) then
      nspec = modes%comps(mm)%nspec
   endif

   ! name of interstitial number mixing ratio
   if (present(num_name)) then
      num_name = modes%comps(mm)%camname_num_a
   endif

   ! name of cloud borne number mixing ratio
   if (present(num_name_cw)) then
      num_name_cw = modes%comps(mm)%camname_num_c
   endif

end subroutine rad_cnst_get_info_by_mode

!================================================================================================

subroutine rad_cnst_get_info_by_mode_spec(list_idx, m_idx, s_idx, &
   spec_type, spec_name, spec_name_cw)

   ! Return info about modal aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of mode in the specified list
   integer,                     intent(in)  :: s_idx       ! index of specie in the specified mode
   character(len=32), optional, intent(out) :: spec_type   ! type of specie
   character(len=32), optional, intent(out) :: spec_name   ! name of interstitial specie
   character(len=32), optional, intent(out) :: spec_name_cw ! name of cloud borne specie

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer          :: nmodes
   integer          :: nspec
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_mode_spec'
   !-----------------------------------------------------------------------------

   m_list => ma_list(list_idx)

   ! check for valid mode index
   nmodes = m_list%nmodes
   if (m_idx < 1 .or. m_idx > nmodes) then
      write(iulog,*) subname//': ERROR - invalid mode index: ', m_idx
      call endrun(subname//': ERROR - invalid mode index')
   end if

   ! get index into the mode definition object
   mm = m_list%idx(m_idx)

   ! check for valid specie index
   nspec = modes%comps(mm)%nspec
   if (s_idx < 1 .or. s_idx > nspec) then
      write(iulog,*) subname//': ERROR - invalid specie index: ', s_idx
      call endrun(subname//': ERROR - invalid specie index')
   end if

   ! specie type
   if (present(spec_type)) then
      spec_type = modes%comps(mm)%type(s_idx)
   endif

   ! interstitial specie name
   if (present(spec_name)) then
      spec_name = modes%comps(mm)%camname_mmr_a(s_idx)
   endif

   ! cloud borne specie name
   if (present(spec_name_cw)) then
      spec_name_cw = modes%comps(mm)%camname_mmr_c(s_idx)
   endif

end subroutine rad_cnst_get_info_by_mode_spec

!================================================================================================

subroutine rad_cnst_get_info_by_spectype(list_idx, spectype, mode_idx, spec_idx)

   ! Return info about modes in the specified climate/diagnostics list

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),            intent(in)  :: spectype    ! species type
   integer,           optional, intent(out) :: mode_idx    ! index of a mode that contains a specie of spectype
   integer,           optional, intent(out) :: spec_idx    ! index of the species of spectype

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer  :: i, nmodes, m_idx, nspec, ispec
   logical  :: found_spectype

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_spectype'
   !-----------------------------------------------------------------------------

   m_list => ma_list(list_idx)

   ! number of modes in specified list
   nmodes = m_list%nmodes

   ! loop through modes in specified climate/diagnostic list
   found_spectype = .false.
   do i = 1, nmodes

      ! get index of the mode in the definition object
      m_idx = m_list%idx(i)

      ! number of species in the mode
      nspec = modes%comps(m_idx)%nspec

      ! loop through species looking for spectype
      do ispec = 1, nspec

         if (trim(modes%comps(m_idx)%type(ispec)) == trim(spectype)) then
            if (present(mode_idx)) mode_idx = i
            if (present(spec_idx)) spec_idx = ispec
            found_spectype = .true.
            exit
         end if
      end do

      if (found_spectype) exit
   end do

   if (.not. found_spectype) then
      if (present(mode_idx)) mode_idx = -1
      if (present(spec_idx)) spec_idx = -1
   end if

end subroutine rad_cnst_get_info_by_spectype

!================================================================================================

function rad_cnst_get_mode_idx(list_idx, mode_type) result(mode_idx)

   ! Return mode index of the specified type in the specified climate/diagnostics list.
   ! Return -1 if not found.

   ! Arguments
   integer,           intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),  intent(in)  :: mode_type   ! mode type

   ! Return value
   integer                        :: mode_idx    ! mode index

   ! Local variables
   type(modelist_t), pointer :: m_list

   integer  :: i, nmodes, m_idx

   character(len=*), parameter :: subname = 'rad_cnst_get_mode_idx'
   !-----------------------------------------------------------------------------

   ! if mode type not found return -1
   mode_idx = -1

   ! specified mode list
   m_list => ma_list(list_idx)

   ! number of modes in specified list
   nmodes = m_list%nmodes

   ! loop through modes in specified climate/diagnostic list
   do i = 1, nmodes

      ! get index of the mode in the definition object
      m_idx = m_list%idx(i)

      ! look in mode definition object (modes) for the mode types
      if (trim(modes%types(m_idx)) == trim(mode_type)) then
         mode_idx = i
         exit
      end if
   end do

end function rad_cnst_get_mode_idx

!================================================================================================

function rad_cnst_get_spec_idx(list_idx, mode_idx, spec_type) result(spec_idx)

   ! Return specie index of the specified type in the specified mode of the specified
   ! climate/diagnostics list.  Return -1 if not found.

   ! Arguments
   integer,           intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,           intent(in)  :: mode_idx    ! mode index
   character(len=*),  intent(in)  :: spec_type   ! specie type

   ! Return value
   integer                        :: spec_idx    ! specie index

   ! Local variables
   type(modelist_t),       pointer :: m_list
   type(mode_component_t), pointer :: mode_comps

   integer  :: i, m_idx, nspec

   character(len=*), parameter :: subname = 'rad_cnst_get_spec_idx'
   !-----------------------------------------------------------------------------

   ! if specie type not found return -1
   spec_idx = -1

   ! modes in specified list
   m_list => ma_list(list_idx)

   ! get index of the specified mode in the definition object
   m_idx = m_list%idx(mode_idx)

   ! object containing the components of the mode
   mode_comps => modes%comps(m_idx)

   ! number of species in specified mode
   nspec = mode_comps%nspec

   ! loop through species in specified mode
   do i = 1, nspec

      ! look in mode definition object (modes) for the mode types
      if (trim(mode_comps%type(i)) == trim(spec_type)) then
         spec_idx = i
         exit
      end if
   end do

end function rad_cnst_get_spec_idx
!================================================================================================

subroutine rad_cnst_get_call_list(call_list)

   ! Return info about which climate/diagnostic calculations are requested

   ! Arguments
   logical, intent(out) :: call_list(0:N_DIAG)
   !-----------------------------------------------------------------------------

   call_list(:) = active_calls(:)

end subroutine rad_cnst_get_call_list

!================================================================================================

subroutine rad_cnst_out(list_idx, state, pbuf)

   ! Output the mass per layer, and total column burdens for gas and aerosol
   ! constituents in either the climate or diagnostic lists

   ! Arguments
   integer,                     intent(in) :: list_idx
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)


   ! Local variables
   integer :: i, naer, ngas, lchnk, ncol
   integer :: idx
   character(len=1)  :: source
   character(len=32) :: name, cbname
   real(r8)          :: mass(pcols,pver)
   real(r8)          :: cb(pcols)
   real(r8), pointer :: mmr(:,:)
   type(aerlist_t), pointer :: aerlist
   type(gaslist_t), pointer :: g_list
   character(len=*), parameter :: subname = 'rad_cnst_out'
   !-----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointer with requested aerosol list
   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   naer = aerlist%numaerosols
   do i = 1, naer

      source = aerlist%aer(i)%source
      idx    = aerlist%aer(i)%idx
      name   = aerlist%aer(i)%mass_name
      ! construct name for column burden field by replacing the 'm_' prefix by 'cb_'
      cbname = 'cb_' // name(3:len_trim(name))

      select case( source )
      case ('A')
         mmr => state%q(:,:,idx)
      case ('N')
         call pbuf_get_field(pbuf, idx, mmr)
      end select

      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga
      call outfld(trim(name), mass, pcols, lchnk)

      cb(:ncol) = sum(mass(:ncol,:),2)
      call outfld(trim(cbname), cb, pcols, lchnk)

   end do

   ! Associate pointer with requested gas list
   g_list => gaslist(list_idx)

   ngas = g_list%ngas
   do i = 1, ngas

      source = g_list%gas(i)%source
      idx    = g_list%gas(i)%idx
      name   = g_list%gas(i)%mass_name
      cbname = 'cb_' // name(3:len_trim(name))
      select case( source )
      case ('A')
         mmr => state%q(:,:,idx)
      case ('N')
         call pbuf_get_field(pbuf, idx, mmr)
      end select

      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga
      call outfld(trim(name), mass, pcols, lchnk)

      cb(:ncol) = sum(mass(:ncol,:),2)
      call outfld(trim(cbname), cb, pcols, lchnk)

   end do

end subroutine rad_cnst_out

!================================================================================================
! Private methods
!================================================================================================

subroutine init_mode_comps(modes)

   ! Initialize the mode definitions by looking up the relevent indices in the
   ! constituent and pbuf arrays, and getting the physprop IDs

   ! Arguments
   type(modes_t), intent(inout) :: modes

   ! Local variables
   integer :: m, ispec, nspec

   character(len=*), parameter :: routine = 'init_modes'
   !-----------------------------------------------------------------------------

   do m = 1, modes%nmodes

      ! indices for number mixing ratio components
      modes%comps(m)%idx_num_a = get_cam_idx(modes%comps(m)%source_num_a, modes%comps(m)%camname_num_a, routine)
      modes%comps(m)%idx_num_c = get_cam_idx(modes%comps(m)%source_num_c, modes%comps(m)%camname_num_c, routine)

      ! allocate memory for species
      nspec = modes%comps(m)%nspec
      allocate( &
         modes%comps(m)%idx_mmr_a(nspec), &
         modes%comps(m)%idx_mmr_c(nspec), &
         modes%comps(m)%idx_props(nspec)  )

      do ispec = 1, nspec

         ! indices for species mixing ratio components
         modes%comps(m)%idx_mmr_a(ispec) = get_cam_idx(modes%comps(m)%source_mmr_a(ispec), &
                                                   modes%comps(m)%camname_mmr_a(ispec), routine)
         modes%comps(m)%idx_mmr_c(ispec) = get_cam_idx(modes%comps(m)%source_mmr_c(ispec), &
                                                   modes%comps(m)%camname_mmr_c(ispec), routine)

         ! get physprop ID
         modes%comps(m)%idx_props(ispec) = physprop_get_id(modes%comps(m)%props(ispec)) 
         if (modes%comps(m)%idx_props(ispec) == -1) then
            call endrun(routine//' : ERROR idx not found for '//trim(modes%comps(m)%props(ispec)))
         end if

      end do

   end do

end subroutine init_mode_comps

!================================================================================================

integer function get_cam_idx(source, name, routine)

   ! get index of name in internal CAM array; either the constituent array
   ! or the physics buffer

   character(len=*), intent(in) :: source
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: routine  ! name of calling routine

   integer :: idx
   integer :: errcode
   !-----------------------------------------------------------------------------
   
   if (source(1:1) == 'N') then

      idx = pbuf_get_index(trim(name),errcode)
      if (errcode < 0) then
         call endrun(routine//' ERROR: cannot find physics buffer field '//trim(name))
      end if

   else if (source(1:1) == 'A') then

      call cnst_get_ind(trim(name), idx, abort=.false.)
      if (idx < 0) then
         call endrun(routine//' ERROR: cannot find constituent field '//trim(name))
      end if

   else if (source(1:1) == 'Z') then

      idx = -1

   else

      call endrun(routine//' ERROR: invalid source for specie '//trim(name))

   end if
  
   get_cam_idx = idx

end function get_cam_idx

!================================================================================================

subroutine list_init1(namelist, gaslist, aerlist, ma_list)

   ! Initialize the gas and bulk and modal aerosol lists with the 
   ! entities specified in the climate or diagnostic lists.

   ! This first phase initialization just sets the information that
   ! is available at the time the namelist is read.

   type(rad_cnst_namelist_t), intent(in) :: namelist ! parsed namelist input for climate or diagnostic lists

   type(gaslist_t),        intent(inout) :: gaslist
   type(aerlist_t),        intent(inout) :: aerlist
   type(modelist_t),       intent(inout) :: ma_list


   ! Local variables
   integer :: ii, idx, m, naero, nmodes
   integer :: igas, ifileindex, ba_idx, ma_idx
   integer :: istat
   character(len=*), parameter :: routine = 'list_init1'
   !-----------------------------------------------------------------------------

   ! nradgas is set by the radiative transfer code
   gaslist%ngas = nradgas

   ! Determine the number of bulk aerosols and aerosol modes in the list
   naero = 0
   nmodes = 0
   do ii = 1, namelist%ncnst
      if (trim(namelist%type(ii)) == 'A') naero  = naero + 1
      if (trim(namelist%type(ii)) == 'M') nmodes = nmodes + 1
   end do
   aerlist%numaerosols = naero
   ma_list%nmodes      = nmodes

   ! allocate storage for the aerosol, gas, and mode lists
   allocate( &
      aerlist%aer(aerlist%numaerosols),      &
      gaslist%gas(gaslist%ngas),             &
      ma_list%idx(ma_list%nmodes),           &
      ma_list%physprop_files(ma_list%nmodes), &
      ma_list%idx_props(ma_list%nmodes),     &
      stat=istat)
   if (istat /= 0) call endrun(routine//': allocate ERROR; aero and gas list components')

   if (masterproc .and. verbose) then
      if (len_trim(gaslist%list_id) == 0) then
         write(iulog,*) nl//' '//routine//': namelist input for climate list'
      else
         write(iulog,*) nl//' '//routine//': namelist input for diagnostic list:'//gaslist%list_id
      end if
   end if

   ! Loop over the radiatively active components specified in the namelist
   ba_idx = 0
   ma_idx = 0
   do ii = 1, namelist%ncnst

      if (masterproc .and. verbose) &
         write(iulog,*) "  rad namelist spec: "// trim(namelist%source(ii)) &
         //":"//trim(namelist%camname(ii))//":"//trim(namelist%radname(ii))

      ! Check that the source specifier is legal.
      if (namelist%source(ii) /= 'A' .and. namelist%source(ii) /= 'M' .and. &
          namelist%source(ii) /= 'N' .and. namelist%source(ii) /= 'Z' ) then
         call endrun(routine//": source must either be A, M, N or Z:"//&
                     " illegal specifier in namelist input: "//namelist%source(ii))
      end if

      ! Add component to appropriate list (gas, modal or bulk aerosol)
      if (namelist%type(ii) == 'A') then 

         ! Add to bulk aerosol list
         ba_idx = ba_idx + 1

         aerlist%aer(ba_idx)%source        = namelist%source(ii)
         aerlist%aer(ba_idx)%camname       = namelist%camname(ii)
         aerlist%aer(ba_idx)%physprop_file = namelist%radname(ii)

      else if (namelist%type(ii) == 'M') then 

         ! Add to modal aerosol list
         ma_idx = ma_idx + 1

         ! Look through the mode definitions for the name of the specified mode.  The
         ! index into the modes object all the information relevent to the mode definition.
         ma_list%idx(ma_idx) = -1
         do m = 1, modes%nmodes
            if (trim(namelist%camname(ii)) == trim(modes%names(m))) then
               ma_list%idx(ma_idx) = m
               exit
            end if
         end do
         if (ma_list%idx(ma_idx) == -1) &
            call endrun(routine//' ERROR cannot find mode name '//trim(namelist%camname(ii)))

         ! Also save the name of the physprop file
         ma_list%physprop_files(ma_idx) = namelist%radname(ii)

      else 

         ! Add to gas list

         ! The radiative transfer code requires the input of a specific set of gases
         ! which is hardwired into the code.  The CAM interface to the RT code uses
         ! the names in the radconstants module to refer to these gases.  The user
         ! interface (namelist) also uses these names to identify the gases treated
         ! by the RT code.  We use the index order set in radconstants for convenience
         ! only.

         ! First check that the gas name specified by the user is allowed.
         ! rad_gas_index will abort on illegal names.
         igas = rad_gas_index(namelist%radname(ii))

         ! Set values in the igas index
         gaslist%gas(igas)%source  = namelist%source(ii)
         gaslist%gas(igas)%camname = namelist%camname(ii)

      end if
   end do

end subroutine list_init1

!================================================================================================

subroutine list_init2(gaslist, aerlist, ma_list)

   ! Final initialization phase gets the component indices in the constituent array
   ! and the physics buffer, and indices into physprop module.

   type(gaslist_t),        intent(inout) :: gaslist
   type(aerlist_t),        intent(inout) :: aerlist
   type(modelist_t),       intent(inout) :: ma_list

   ! Local variables
   integer :: i
   character(len=*), parameter :: routine = 'list_init2'
   !-----------------------------------------------------------------------------

   ! Loop over gases
   do i = 1, gaslist%ngas

      ! locate the specie mixing ratio in the pbuf or state
      gaslist%gas(i)%idx = get_cam_idx(gaslist%gas(i)%source, gaslist%gas(i)%camname, routine)

   end do

   ! Loop over bulk aerosols
   do i = 1, aerlist%numaerosols

      ! locate the specie mixing ratio in the pbuf or state
      aerlist%aer(i)%idx = get_cam_idx(aerlist%aer(i)%source, aerlist%aer(i)%camname, routine)

      ! get the physprop_id from the phys_prop module
      aerlist%aer(i)%physprop_id = physprop_get_id(aerlist%aer(i)%physprop_file)

   end do

   ! Loop over modes
   do i = 1, ma_list%nmodes

      ! get the physprop_id from the phys_prop module
      ma_list%idx_props(i) = physprop_get_id(ma_list%physprop_files(i))

   end do

end subroutine list_init2

!================================================================================================

subroutine rad_gas_diag_init(glist)

! Add diagnostic fields to the master fieldlist.

   type(gaslist_t), intent(inout) :: glist

   integer :: i, ngas
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   ngas = glist%ngas
   if (ngas == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = glist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, ngas

      ! construct names for mass per layer diagnostics
      name = 'm_' // trim(glist%gas(i)%camname) // trim(suffix)
      glist%gas(i)%mass_name = name
      long_name = trim(glist%gas(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), (/ 'lev' /), 'A', 'kg/m^2', trim(long_name))

      ! construct names for column burden diagnostics
      name = 'cb_' // trim(glist%gas(i)%camname) // trim(suffix)
      long_name = trim(glist%gas(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), horiz_only, 'A', 'kg/m^2', trim(long_name))

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_gas_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_gas_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_gas_diag_init

!================================================================================================

subroutine rad_aer_diag_init(alist)

! Add diagnostic fields to the master fieldlist.

   type(aerlist_t), intent(inout) :: alist

   integer :: i, naer
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   naer = alist%numaerosols
   if (naer == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = alist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, naer

      ! construct names for mass per layer diagnostic fields
      name = 'm_' // trim(alist%aer(i)%camname) // trim(suffix)
      alist%aer(i)%mass_name = name
      long_name = trim(alist%aer(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), (/ 'lev' /), 'A', 'kg/m^2', trim(long_name))

      ! construct names for column burden diagnostic fields
      name = 'cb_' // trim(alist%aer(i)%camname) // trim(suffix)
      long_name = trim(alist%aer(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), horiz_only, 'A', 'kg/m^2', trim(long_name))

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_aer_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_aer_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_aer_diag_init


!================================================================================================

subroutine parse_mode_defs(nl_in, modes)

   ! Parse the mode definition specifiers.  The specifiers are of the form:
   ! 
   ! 'mode_name:mode_type:=',
   !  'source_num_a:camname_num_a:source_num_c:camname_num_c:num_mr:+',
   !  'source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file[:+]'[,]
   !  ['source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file][:+][']
   !
   ! where the ':' separated fields are:
   ! mode_name -- name of the mode.
   ! mode_type -- type of mode.  Valid values are from the MAM code.
   ! =         -- this line terminator identifies the initial string in a
   !              mode definition
   ! +         -- this line terminator indicates that the mode definition is
   !              continued in the next string
   ! source_num_a  -- Source of interstitial number mixing ratio,  'A', 'N', or 'Z'
   ! camname_num_a -- the name of the interstitial number component.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_num_c  -- Source of cloud borne number mixing ratio,  'A', 'N', or 'Z'
   ! camname_num_c -- the name of the cloud borne number component.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_mmr_a  -- Source of interstitial specie mass mixing ratio,  'A', 'N' or 'Z'
   ! camname_mmr_a -- the name of the interstitial specie.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_mmr_c  -- Source of cloud borne specie mass mixing ratio,  'A', 'N' or 'Z'
   ! camname_mmr_c -- the name of the cloud borne specie.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! spec_type -- species type.  Valid values far from the MAM code, except that
   !              the value 'num_mr' designates a number mixing ratio and has no
   !              associated field for the prop_file.  There can only be one entry
   !              with the num_mr type in a mode definition.
   ! prop_file -- For aerosol species this is a filename, which is
   !              identified by a ".nc" suffix.  The file contains optical and 
   !              other physical properties of the aerosol.
   !
   ! A mode definition must contain only 1 string for the number mixing ratio components
   ! and at least 1 string for the species.


   character(len=*), intent(inout) :: nl_in(:)    ! namelist input (blanks are removed on output)
   type(modes_t),    intent(inout) :: modes       ! structure containing parsed input

   ! Local variables
   integer :: i, m
   integer :: istat
   integer :: nmodes, nstr, istr
   integer :: mbeg, mcur
   integer :: nspec, ispec
   integer :: strlen, ibeg, iend, ipos
   logical :: num_mr_found
   character(len=*), parameter :: routine = 'parse_mode_defs'
   character(len=len(nl_in(1))) :: tmpstr
   character(len=1)  :: tmp_src_a
   character(len=32) :: tmp_name_a
   character(len=1)  :: tmp_src_c
   character(len=32) :: tmp_name_c
   character(len=32) :: tmp_type
   !-------------------------------------------------------------------------
  
   ! Determine number of modes defined by counting number of strings that are
   ! terminated by ':='
   ! (algorithm stops counting at first blank element).
   nmodes = 0
   nstr = 0
   do m = 1, n_mode_str

      if (len_trim(nl_in(m)) == 0) exit
      nstr = nstr + 1
      
      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(nl_in(m))
      nl_in(m) = tmpstr
      do
         strlen = len_trim(nl_in(m))
         ipos = index(nl_in(m), ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = nl_in(m)(:ipos-1) // nl_in(m)(ipos+1:strlen)
         nl_in(m) = tmpstr
      end do
      ! count strings with ':=' terminator
      if (nl_in(m)(strlen-1:strlen) == ':=') nmodes = nmodes + 1

   end do
   modes%nmodes = nmodes

   ! return if no modes defined
   if (nmodes == 0) return

   ! allocate components that depend on nmodes
   allocate( &
      modes%names(nmodes),  &
      modes%types(nmodes),  &
      modes%comps(nmodes),  &
      stat=istat )
   if (istat > 0) then
      write(iulog,*) routine//': ERROR: cannot allocate storage for modes.  nmodes=', nmodes
      call endrun(routine//': ERROR allocating storage for modes')
   end if
   

   mcur = 1              ! index of current string being processed

   ! loop over modes
   do m = 1, nmodes

      mbeg = mcur  ! remember the first string of a mode

      ! check that first string in mode definition is ':=' terminated
      iend = len_trim(nl_in(mcur))
      if (nl_in(mcur)(iend-1:iend) /= ':=') call parse_error('= not found', nl_in(mcur))

      ! count species in mode definition.  definition will contain 1 string with
      ! with a ':+' terminator for each specie
      nspec = 0
      mcur = mcur + 1
      do
         iend = len_trim(nl_in(mcur))
         if (nl_in(mcur)(iend-1:iend) /= ':+') exit
         nspec = nspec + 1
         mcur = mcur + 1
      end do
      
      ! a mode must have at least one specie
      if (nspec == 0) call parse_error('mode must have at least one specie', nl_in(mbeg))

      ! allocate components that depend on number of species
      allocate( &
         modes%comps(m)%source_mmr_a(nspec),  &
         modes%comps(m)%camname_mmr_a(nspec), &
         modes%comps(m)%source_mmr_c(nspec),  &
         modes%comps(m)%camname_mmr_c(nspec), &
         modes%comps(m)%type(nspec),          &
         modes%comps(m)%props(nspec),         &
         stat=istat)

      if (istat > 0) then
         write(iulog,*) routine//': ERROR: cannot allocate storage for species.  nspec=', nspec
         call endrun(routine//': ERROR allocating storage for species')
      end if

      ! initialize components
      modes%comps(m)%nspec         = nspec
      modes%comps(m)%source_num_a  = ' '
      modes%comps(m)%camname_num_a = ' '
      modes%comps(m)%source_num_c  = ' '
      modes%comps(m)%camname_num_c = ' '
      do ispec = 1, nspec
         modes%comps(m)%source_mmr_a(ispec)  = ' '
         modes%comps(m)%camname_mmr_a(ispec) = ' '
         modes%comps(m)%source_mmr_c(ispec)  = ' '
         modes%comps(m)%camname_mmr_c(ispec) = ' '
         modes%comps(m)%type(ispec)          = ' '
         modes%comps(m)%props(ispec)         = ' '
      end do

      ! return to first string in mode definition
      mcur = mbeg
      tmpstr = nl_in(mcur)
      
      ! mode name
      ipos = index(tmpstr, ':')
      if (ipos < 2) call parse_error('mode name not found', tmpstr)
      modes%names(m) = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! mode type
      ipos = index(tmpstr, ':')
      if (ipos == 0) call parse_error('mode type not found', tmpstr)
      ! check for valid mode type
      call check_mode_type(tmpstr, 1, ipos-1)
      modes%types(m) = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! mode type must be followed by '='
      if (tmpstr(1:1) /= '=') call parse_error('= not found', tmpstr)

      ! move to next string
      mcur = mcur + 1
      tmpstr = nl_in(mcur)

      ! process mode component strings
      num_mr_found = .false.   ! keep track of whether number mixing ratio component is found
      ispec = 0                ! keep track of the number of species found
      do

         ! source of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find source field first', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_a = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_a = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! source of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find a source field', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_c = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_c = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! component type
         ipos = scan(tmpstr, ': ')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)

         if (tmpstr(:ipos-1) == 'num_mr') then

            ! there can only be one number mixing ratio component
            if (num_mr_found) call parse_error('more than 1 number component', nl_in(mcur))

            num_mr_found = .true.
            modes%comps(m)%source_num_a  = tmp_src_a
            modes%comps(m)%camname_num_a = tmp_name_a
            modes%comps(m)%source_num_c  = tmp_src_c
            modes%comps(m)%camname_num_c = tmp_name_c
            tmpstr                       = tmpstr(ipos+1:)

         else

            ! check for valid specie type
            call check_specie_type(tmpstr, 1, ipos-1)
            tmp_type = tmpstr(:ipos-1)
            tmpstr   = tmpstr(ipos+1:)

            ! get the properties file
            ipos = scan(tmpstr, ': ')
            if (ipos == 0) call parse_error('next separator not found', tmpstr)
            ! check for valid filename -- must have .nc extension
            if (tmpstr(ipos-3:ipos-1) /= '.nc') &
               call parse_error('filename not valid', tmpstr)

            ispec = ispec + 1
            modes%comps(m)%source_mmr_a(ispec)  = tmp_src_a
            modes%comps(m)%camname_mmr_a(ispec) = tmp_name_a
            modes%comps(m)%source_mmr_c(ispec)  = tmp_src_c
            modes%comps(m)%camname_mmr_c(ispec) = tmp_name_c
            modes%comps(m)%type(ispec)          = tmp_type
            modes%comps(m)%props(ispec)         = tmpstr(:ipos-1)
            tmpstr                              = tmpstr(ipos+1:)
         end if

         ! check if there are more components.  either the current character is
         ! a ' ' which means this string is the final mode component, or the character
         ! is a '+' which means there are more components
         if (tmpstr(1:1) == ' ') exit

         if (tmpstr(1:1) /= '+') &
               call parse_error('+ field not found', tmpstr)

         ! continue to next component...
         mcur = mcur + 1
         tmpstr = nl_in(mcur)
      end do

      ! check that a number component was found
      if (.not. num_mr_found) call parse_error('number component not found', nl_in(mbeg))

      ! check that the right number of species were found
      if (ispec /= nspec) call parse_error('component parsing got wrong number of species', nl_in(mbeg))

      ! continue to next mode...
      mcur = mcur + 1
      tmpstr = nl_in(mcur)
   end do

   !------------------------------------------------------------------------------------------------
   contains
   !------------------------------------------------------------------------------------------------

   ! internal subroutines used for error checking and reporting

   subroutine parse_error(msg, str)

      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: str

      write(iulog,*) routine//': ERROR: '//msg
      write(iulog,*) ' input string: '//trim(str)
      call endrun(routine//': ERROR: '//msg)

   end subroutine parse_error

   !------------------------------------------------------------------------------------------------

   subroutine check_specie_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie
   
      integer :: i

      do i = 1, num_spec_types
         if (str(ib:ie) == trim(spec_type_names(i))) return
      end do

      call parse_error('specie type not valid', str(ib:ie))

   end subroutine check_specie_type

   !------------------------------------------------------------------------------------------------

   subroutine check_mode_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie  ! begin, end character of mode type substring
   
      integer :: i

      do i = 1, num_mode_types
         if (str(ib:ie) == trim(mode_type_names(i))) return
      end do

      call parse_error('mode type not valid', str(ib:ie))

   end subroutine check_mode_type

   !------------------------------------------------------------------------------------------------

end subroutine parse_mode_defs

!================================================================================================

subroutine parse_rad_specifier(specifier, namelist_data)

!-----------------------------------------------------------------------------
! Private method for parsing the radiation namelist specifiers.  The specifiers
! are of the form 'source_camname:radname' where:
! source  -- either 'N' for pbuf (non-advected) or 'A' for state (advected)
! camname -- the name of a constituent that must be found in the constituent
!            component of the state when source=A or in the physics buffer
!            when source=N
! radname -- For gases this is a name that identifies the constituent to the
!            radiative transfer codes.  These names are contained in the
!            radconstants module.  For aerosols this is a filename, which is
!            identified by a ".nc" suffix.  The file contains optical and 
!            other physical properties of the aerosol.
!
! This code also identifies whether the constituent is a gas or an aerosol
! and adds that info to a structure that stores the parsed data.
!-----------------------------------------------------------------------------

    character(len=*), dimension(:), intent(in) :: specifier
    type(rad_cnst_namelist_t),   intent(inout) :: namelist_data

    ! Local variables
    integer            :: number, i, j
    integer            :: ipos, strlen
    integer            :: astat
    character(len=cs1) :: tmpstr
    character(len=1)   :: source(n_rad_cnst)
    character(len=64)  :: camname(n_rad_cnst)
    character(len=cs1) :: radname(n_rad_cnst)
    character(len=1)   :: type(n_rad_cnst)
    !-------------------------------------------------------------------------
  
    number = 0

    parse_loop: do i = 1, n_rad_cnst
      if ( len_trim(specifier(i)) == 0 ) then 
         exit parse_loop
      endif

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(specifier(i))
      do
         strlen = len_trim(tmpstr)
         ipos = index(tmpstr, ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = tmpstr(:ipos-1) // tmpstr(ipos+1:strlen)
      end do

      ! Locate the ':' separating source from camname.
      j = index(tmpstr, ':')
      source(i) = tmpstr(:j-1)
      tmpstr = tmpstr(j+1:)

      ! locate the ':' separating camname from radname
      j = scan(tmpstr, ':')
 
      camname(i) = tmpstr(:j-1)
      radname(i) = tmpstr(j+1:)

      ! determine the type of constituent
      if (source(i) == 'M') then 
         type(i) = 'M'
      else if(index(radname(i),".nc") .gt. 0) then
         type(i) = 'A'
      else
         type(i) = 'G'
      end if

      number = number+1    
    end do parse_loop

    namelist_data%ncnst = number

    if (number == 0) return

    allocate(namelist_data%source (number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%source')
    allocate(namelist_data%camname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%camname')
    allocate(namelist_data%radname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%radname')
    allocate(namelist_data%type(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%type')

    namelist_data%source(:namelist_data%ncnst)  = source (:namelist_data%ncnst)
    namelist_data%camname(:namelist_data%ncnst) = camname(:namelist_data%ncnst)
    namelist_data%radname(:namelist_data%ncnst) = radname(:namelist_data%ncnst)
    namelist_data%type(:namelist_data%ncnst)    = type(:namelist_data%ncnst)

end subroutine parse_rad_specifier

!================================================================================================

subroutine rad_cnst_get_aer_mmr_by_idx(list_idx, aer_idx, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: aer_idx
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: lchnk
   integer :: idx
   character(len=1) :: source
   type(aerlist_t), pointer :: aerlist
   character(len=*), parameter :: subname = 'rad_cnst_get_aer_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   lchnk = state%lchnk

   ! Check for valid input aerosol index
   if (aer_idx < 1  .or.  aer_idx > aerlist%numaerosols) then
      write(iulog,*) subname//': aer_idx= ', aer_idx, '  numaerosols= ', aerlist%numaerosols
      call endrun(subname//': aerosol list index out of range')
   end if

   ! Get data source
   source = aerlist%aer(aer_idx)%source
   idx    = aerlist%aer(aer_idx)%idx
   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_aer_mmr_by_idx

!================================================================================================

subroutine rad_cnst_get_mam_mmr_by_idx(list_idx, mode_idx, spec_idx, phase, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the modal aerosol specie from the specified
   ! climate or diagnostic list.  

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: mode_idx    ! mode index
   integer,                     intent(in) :: spec_idx    ! index of specie in the mode
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > modes%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', modes%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Get data source
   if (phase == 'a') then
      source = modes%comps(m_idx)%source_mmr_a(spec_idx)
      idx    = modes%comps(m_idx)%idx_mmr_a(spec_idx)
   else if (phase == 'c') then
      source = modes%comps(m_idx)%source_mmr_c(spec_idx)
      idx    = modes%comps(m_idx)%idx_mmr_c(spec_idx)
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_mam_mmr_by_idx

!================================================================================================

subroutine rad_cnst_get_mam_mmr_idx(mode_idx, spec_idx, idx)

   ! Return constituent index of mam specie mass mixing ratio for aerosol modes in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the 
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   ! Arguments
   integer, intent(in)  :: mode_idx    ! mode index
   integer, intent(in)  :: spec_idx    ! index of specie in the mode
   integer, intent(out) :: idx         ! index of specie in the constituent array

   ! Local variables
   integer :: m_idx
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_mmr_idx'
   !-----------------------------------------------------------------------------

   ! assume climate list (i.e., species are in the constituent array)
   mlist => ma_list(0)

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > modes%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', modes%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Assume data source is interstitial since that's what's in the constituent array
   idx    = modes%comps(m_idx)%idx_mmr_a(spec_idx)

end subroutine rad_cnst_get_mam_mmr_idx

!================================================================================================

subroutine rad_cnst_get_mode_num(list_idx, mode_idx, phase, state, pbuf, num)

   ! Return pointer to number mixing ratio for the aerosol mode from the specified
   ! climate or diagnostic list.  

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: mode_idx    ! mode index
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: num(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_num'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Get data source
   if (phase == 'a') then
      source = modes%comps(m_idx)%source_num_a
      idx    = modes%comps(m_idx)%idx_num_a
   else if (phase == 'c') then
      source = modes%comps(m_idx)%source_num_c
      idx    = modes%comps(m_idx)%idx_num_c
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      num => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, num)
   case ('Z')
      num => zero_cols
   end select

end subroutine rad_cnst_get_mode_num

!================================================================================================

subroutine rad_cnst_get_mode_num_idx(mode_idx, cnst_idx)

   ! Return constituent index of mode number mixing ratio for the aerosol mode in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the 
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   ! Arguments
   integer,  intent(in)  :: mode_idx    ! mode index
   integer,  intent(out) :: cnst_idx    ! constituent index

   ! Local variables
   integer :: m_idx
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_num'
   !-----------------------------------------------------------------------------

   ! assume climate list
   mlist => ma_list(0)

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check that source is 'A' which means the index is for the constituent array
   source = modes%comps(m_idx)%source_num_a
   if (source /= 'A') then
      write(iulog,*) subname//': source= ', source
      call endrun(subname//': requested mode number index not in constituent array')
   end if

   ! Return index in constituent array
   cnst_idx = modes%comps(m_idx)%idx_num_a

end subroutine rad_cnst_get_mode_num_idx

!================================================================================================

integer function rad_cnst_get_aer_idx(list_idx, aer_name)

   ! Return the index of aerosol aer_name in the list specified by list_idx.

    ! Arguments
   integer,             intent(in) :: list_idx    ! 0 for climate list, 1-N_DIAG for diagnostic lists
   character(len=*),    intent(in) :: aer_name    ! aerosol name (in state or pbuf)

   ! Local variables
   integer :: i, aer_idx
   type(aerlist_t), pointer :: aerlist
   character(len=*), parameter :: subname = "rad_cnst_get_aer_idx"
   !-------------------------------------------------------------------------
   
   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Get index in aerosol list for requested name
   aer_idx = -1
   do i = 1, aerlist%numaerosols
      if (trim(aer_name) == trim(aerlist%aer(i)%camname)) then
         aer_idx = i
         exit
      end if
   end do

   if (aer_idx == -1) call endrun(subname//": ERROR - name not found")
  
   rad_cnst_get_aer_idx = aer_idx

end function rad_cnst_get_aer_idx

!================================================================================================

subroutine rad_cnst_get_aer_props_by_idx(list_idx, &
   aer_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, num_to_mass_aer)

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   use phys_prop, only: physprop_get


   ! Arguments
   integer,                     intent(in)  :: list_idx ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: aer_idx  ! index of the aerosol
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:) 
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:) 
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:) 
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)         
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)         
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)
   character(len=20), optional, intent(out) :: aername           
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer        
   real(r8),          optional, intent(out) :: dispersion_aer    
   real(r8),          optional, intent(out) :: num_to_mass_aer   

   real(r8),          optional, pointer     :: r_sw_ext(:,:)         
   real(r8),          optional, pointer     :: r_sw_scat(:,:)         
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)         
   real(r8),          optional, pointer     :: r_lw_abs(:,:)         
   real(r8),          optional, pointer     :: mu(:)         

   ! Local variables
   integer :: id
   character(len=*), parameter :: subname = 'rad_cnst_get_aer_props_by_idx'
   type(aerlist_t), pointer :: aerlist
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   if (aer_idx < 1 .or. aer_idx > aerlist%numaerosols) then
      write(iulog,*) subname//': aerosol list index out of range: ', aer_idx ,' list index: ',list_idx
      call endrun(subname//': aer_idx out of range')
   end if

   id = aerlist%aer(aer_idx)%physprop_id

   if (present(opticstype))        call physprop_get(id, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(id, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(id, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(id, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(id, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(id, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(id, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(id, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(id, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(id, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(id, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(id, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(id, refindex_aer_lw=refindex_aer_lw)

   if (present(aername))           call physprop_get(id, aername=aername)
   if (present(density_aer))       call physprop_get(id, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(id, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(id, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(id, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(id, num_to_mass_aer=num_to_mass_aer)

   if (present(r_lw_abs))          call physprop_get(id, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(id, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(id, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(id, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(id, mu=mu)

end subroutine rad_cnst_get_aer_props_by_idx

!================================================================================================

subroutine rad_cnst_get_mam_props_by_idx(list_idx, &
   mode_idx, spec_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, &
   num_to_mass_aer, spectype)

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   use phys_prop, only: physprop_get

   ! Arguments
   integer,                     intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: mode_idx  ! mode index
   integer,                     intent(in)  :: spec_idx  ! index of specie in the mode
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:) 
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:) 
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:) 
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)         
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)         
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)

   real(r8),          optional, pointer     :: r_sw_ext(:,:)         
   real(r8),          optional, pointer     :: r_sw_scat(:,:)         
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)         
   real(r8),          optional, pointer     :: r_lw_abs(:,:)         
   real(r8),          optional, pointer     :: mu(:)         

   character(len=20), optional, intent(out) :: aername           
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer        
   real(r8),          optional, intent(out) :: dispersion_aer    
   real(r8),          optional, intent(out) :: num_to_mass_aer   
   character(len=32), optional, intent(out) :: spectype

   ! Local variables
   integer :: m_idx, id
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_props_by_idx'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > modes%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', modes%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   id = modes%comps(m_idx)%idx_props(spec_idx)

   if (present(opticstype))        call physprop_get(id, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(id, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(id, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(id, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(id, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(id, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(id, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(id, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(id, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(id, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(id, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(id, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(id, refindex_aer_lw=refindex_aer_lw)

   if (present(r_lw_abs))          call physprop_get(id, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(id, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(id, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(id, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(id, mu=mu)

   if (present(aername))           call physprop_get(id, aername=aername)
   if (present(density_aer))       call physprop_get(id, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(id, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(id, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(id, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(id, num_to_mass_aer=num_to_mass_aer)

   if (present(spectype)) spectype = modes%comps(m_idx)%type(spec_idx)

end subroutine rad_cnst_get_mam_props_by_idx

!================================================================================================

subroutine rad_cnst_get_mode_props(list_idx, mode_idx, &
   extpsw, abspsw, asmpsw, absplw, refrtabsw, &
   refitabsw, refrtablw, refitablw, ncoef, prefr, &
   prefi, sigmag, dgnum, dgnumlo, dgnumhi, &
   rhcrystal, rhdeliques)

   ! Return requested properties for the mode from the specified
   ! climate or diagnostic list.

   use phys_prop, only: physprop_get

   ! Arguments
   integer,             intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,             intent(in)  :: mode_idx  ! mode index

   real(r8),  optional, pointer     :: extpsw(:,:,:,:)
   real(r8),  optional, pointer     :: abspsw(:,:,:,:)
   real(r8),  optional, pointer     :: asmpsw(:,:,:,:)
   real(r8),  optional, pointer     :: absplw(:,:,:,:)
   real(r8),  optional, pointer     :: refrtabsw(:,:)
   real(r8),  optional, pointer     :: refitabsw(:,:)
   real(r8),  optional, pointer     :: refrtablw(:,:)
   real(r8),  optional, pointer     :: refitablw(:,:)
   integer,   optional, intent(out) :: ncoef
   integer,   optional, intent(out) :: prefr
   integer,   optional, intent(out) :: prefi
   real(r8),  optional, intent(out) :: sigmag
   real(r8),  optional, intent(out) :: dgnum
   real(r8),  optional, intent(out) :: dgnumlo
   real(r8),  optional, intent(out) :: dgnumhi
   real(r8),  optional, intent(out) :: rhcrystal
   real(r8),  optional, intent(out) :: rhdeliques

   ! Local variables
   integer :: id
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_props'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the physprop index for the requested mode
   id = mlist%idx_props(mode_idx)

   if (present(extpsw))      call physprop_get(id, extpsw=extpsw)
   if (present(abspsw))      call physprop_get(id, abspsw=abspsw)
   if (present(asmpsw))      call physprop_get(id, asmpsw=asmpsw)
   if (present(absplw))      call physprop_get(id, absplw=absplw)

   if (present(refrtabsw))   call physprop_get(id, refrtabsw=refrtabsw)
   if (present(refitabsw))   call physprop_get(id, refitabsw=refitabsw)
   if (present(refrtablw))   call physprop_get(id, refrtablw=refrtablw)
   if (present(refitablw))   call physprop_get(id, refitablw=refitablw)

   if (present(ncoef))       call physprop_get(id, ncoef=ncoef)
   if (present(prefr))       call physprop_get(id, prefr=prefr)
   if (present(prefi))       call physprop_get(id, prefi=prefi)
   if (present(sigmag))      call physprop_get(id, sigmag=sigmag)
   if (present(dgnum))       call physprop_get(id, dgnum=dgnum)
   if (present(dgnumlo))     call physprop_get(id, dgnumlo=dgnumlo)
   if (present(dgnumhi))     call physprop_get(id, dgnumhi=dgnumhi)
   if (present(rhcrystal))   call physprop_get(id, rhcrystal=rhcrystal)
   if (present(rhdeliques))  call physprop_get(id, rhdeliques=rhdeliques)

end subroutine rad_cnst_get_mode_props

!================================================================================================

subroutine print_modes(modes)

   type(modes_t), intent(inout) :: modes

   integer :: i, m
   !---------------------------------------------------------------------------------------------

   write(iulog,*)' Mode Definitions'

   do m = 1, modes%nmodes

      write(iulog,*) nl//' name=',trim(modes%names(m)),'  type=',trim(modes%types(m))
      write(iulog,*) ' src_a=',trim(modes%comps(m)%source_num_a),'  num_a=',trim(modes%comps(m)%camname_num_a), &
                     ' src_c=',trim(modes%comps(m)%source_num_c),'  num_c=',trim(modes%comps(m)%camname_num_c)

      do i = 1, modes%comps(m)%nspec

         write(iulog,*) ' src_a=',trim(modes%comps(m)%source_mmr_a(i)), '  mmr_a=',trim(modes%comps(m)%camname_mmr_a(i)), &
                       '  src_c=',trim(modes%comps(m)%source_mmr_c(i)), '  mmr_c=',trim(modes%comps(m)%camname_mmr_c(i)), &
                       '  type=',trim(modes%comps(m)%type(i))
         write(iulog,*) '     prop file=', trim(modes%comps(m)%props(i))
      end do

   end do

end subroutine print_modes

!================================================================================================

subroutine print_lists(gas_list, aer_list, ma_list)

   ! Print summary of gas, bulk and modal aerosol lists.  This is just the information
   ! read from the namelist.

   use radconstants, only: gascnst=>gaslist

   type(aerlist_t),  intent(in) :: aer_list
   type(gaslist_t),  intent(in) :: gas_list
   type(modelist_t), intent(in) :: ma_list

   integer :: i, id

   if (len_trim(gas_list%list_id) == 0) then
      write(iulog,*) nl//' gas list for climate calculations'
   else
      write(iulog,*) nl//' gas list for diag'//gas_list%list_id//' calculations'
   end if

   do i = 1, nradgas
      if (gas_list%gas(i)%source .eq. 'N') then
         write(iulog,*) '  '//gas_list%gas(i)%source//':'//gascnst(i)//' has pbuf name:'//&
                        trim(gas_list%gas(i)%camname)
      else if (gas_list%gas(i)%source .eq. 'A') then
         write(iulog,*) '  '//gas_list%gas(i)%source//':'//gascnst(i)//' has constituents name:'//&
                        trim(gas_list%gas(i)%camname)
      endif
   enddo

   if (len_trim(aer_list%list_id) == 0) then
      write(iulog,*) nl//' bulk aerosol list for climate calculations'
   else
      write(iulog,*) nl//' bulk aerosol list for diag'//aer_list%list_id//' calculations'
   end if

   do i = 1, aer_list%numaerosols
      write(iulog,*) '  '//trim(aer_list%aer(i)%source)//':'//trim(aer_list%aer(i)%camname)//&
                     ' optics and phys props in :'//trim(aer_list%aer(i)%physprop_file)
   enddo

   if (len_trim(ma_list%list_id) == 0) then
      write(iulog,*) nl//' modal aerosol list for climate calculations'
   else
      write(iulog,*) nl//' modal aerosol list for diag'//ma_list%list_id//' calculations'
   end if

   do i = 1, ma_list%nmodes
      id = ma_list%idx(i)
      write(iulog,*) '  '//trim(modes%names(id))
   enddo

end subroutine print_lists

!================================================================================================

end module rad_constituents
