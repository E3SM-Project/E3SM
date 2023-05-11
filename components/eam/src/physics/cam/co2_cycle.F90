

module co2_cycle

!------------------------------------------------------------------------------------------------
 
! CO2 was used in radiation calculation.
!
! Purpose:
! Provides distributions of CO2_LND, CO2_OCN, CO2_FF, CO2
! Read co2 flux from ocn and fossil fuel.
! Get  co2 flux from lnd through coupler. 
!
! Authors: Jeff Lee, Keith Lindsay, Balwinder Singh 
! 
!                                              
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8, cxx =>SHR_KIND_CXX, cl =>SHR_KIND_CL
use co2_data_flux,  only: co2_data_flux_type

implicit none
private
save
 
! Public interfaces
public co2_cycle_readnl              ! read the namelist 
public co2_register                  ! register consituents
public co2_transport                 ! turn on co2 tracers transport
public co2_implements_cnst           ! returns true if consituent is implemented by this package
public co2_init_cnst                 ! initialize mixing ratios if not read from initial file
public co2_init                      ! initialize (history) variables
public co2_time_interp_ocn           ! time interpolate co2 flux
public co2_time_interp_fuel          ! time interpolate co2 flux
public co2_cycle_set_ptend           ! set tendency from aircraft emissions
public co2_cycle_set_cnst_type       ! set co2 tracers mixing type for local versions of cnst_type

! Public data
 
public data_flux_ocn                 ! data read in for co2 flux from ocn
public data_flux_fuel                ! data read in for co2 flux from fuel
 
TYPE(co2_data_flux_type) :: data_flux_ocn                  
TYPE(co2_data_flux_type) :: data_flux_fuel
                         
public c_i                           ! global index for new constituents
public co2_readFlux_ocn              ! read co2 flux from data file 
public co2_readFlux_fuel             ! read co2 flux from data file 
public co2_readFlux_aircraft         ! read co2 flux from data file
public co2_print_diags_timestep      ! print out co2 conservation diagnostics every timestep
public co2_print_diags_monthly       ! print out co2 conservation diagnostics monthly
public co2_print_diags_total         ! print out a co2 conservation check for the full run
public co2_conserv_error_tol_per_year! error tolerance for co2 conservation


! Namelist variables
logical :: co2_flag              = .false.       ! true => turn on co2 code, namelist variable
logical :: co2_readFlux_ocn      = .false.       ! true => read co2 flux from ocn,  namelist variable
logical :: co2_readFlux_fuel     = .false.       ! true => read co2 flux from fuel, namelist variable
logical :: co2_readFlux_aircraft = .false.       ! true => read aircraft co2 flux from date file, namelist variable
logical :: co2_print_diags_timestep = .false.    ! true => print out co2 conservation diagnostics every timestep
logical :: co2_print_diags_monthly  = .false.    ! true => print out co2 conservation diagnostics monthly
logical :: co2_print_diags_total    = .false.    ! true => print out a co2 conservation check for the full run
real(r8) :: co2_conserv_error_tol_per_year       ! error tolerance for co2 conservation

character(len=cl) :: co2flux_ocn_file  = 'unset' ! co2 flux from ocn
character(len=cl) :: co2flux_fuel_file = 'unset' ! co2 flux from fossil fuel                       

!-----------------------------------------------------------------------
! new constituents
integer, parameter :: ncnst = 4                      ! number of constituents implemented

character(len=7), dimension(ncnst), parameter :: & ! constituent names
     c_names = (/'CO2_OCN', 'CO2_FFF', 'CO2_LND', 'CO2    '/)

integer :: co2_ocn_glo_ind ! global index of 'CO2_OCN'
integer :: co2_fff_glo_ind ! global index of 'CO2_FFF'
integer :: co2_lnd_glo_ind ! global index of 'CO2_LND'
integer :: co2_glo_ind     ! global index of 'CO2'

integer, dimension(ncnst) :: c_i                   ! global index

!================================================================================================
contains
!================================================================================================

subroutine co2_cycle_readnl(nlfile)

   !-----------------------------------------------------------------------
   ! Read co2_cycle_nl namelist group.
   !
   ! Called by:
   !    runtime_opts.F90
   !-----------------------------------------------------------------------

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: masterproc, mpicom, mstrid=>masterprocid, mpi_logical, & 
                              mpi_character, mpi_real8
   use srf_field_check, only: active_Faoo_fco2_ocn
   use shr_log_mod ,    only: errMsg => shr_log_errMsg
   use cam_abortutils,  only: endrun

   !arguments
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer            :: unitn, ierr
   character(len=cxx) :: err_str
   character(len=*), parameter :: subname = 'co2_cycle_readnl'

   namelist /co2_cycle_nl/ co2_flag, co2_readFlux_ocn, co2_readFlux_fuel, co2_readFlux_aircraft, &
                           co2flux_ocn_file, co2flux_fuel_file, & 
                           co2_print_diags_timestep, co2_print_diags_monthly, co2_print_diags_total, &
                           co2_conserv_error_tol_per_year
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'co2_cycle_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, co2_cycle_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist'//errmsg(__FILE__,__LINE__))
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD

   !Broadcast namelist variables
   call mpi_bcast(co2_flag,                               1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_flag"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2_readFlux_ocn,                       1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_readFlux_ocn"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2_readFlux_fuel,                      1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_readFlux_fuel"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2_readFlux_aircraft,                  1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_readFlux_aircraft"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2flux_ocn_file,   len(co2flux_ocn_file),   mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_ocn_file"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2flux_fuel_file, len(co2flux_fuel_file),   mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_file"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2_print_diags_timestep,               1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_print_diags_timestep"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2_print_diags_monthly,                1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_print_diags_monthly"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2_print_diags_total,                  1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_print_diags_total"//errmsg(__FILE__,__LINE__))
   call mpi_bcast(co2_conserv_error_tol_per_year,         1,   mpi_real8,     mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_conserv_error_tol_per_year"//errmsg(__FILE__,__LINE__))
#endif


   ! Consistency check
   if (co2_readFlux_ocn .and. active_Faoo_fco2_ocn) then
      err_str = subname//': ERROR: reading ocn flux dataset is enabled, but coupler is setting'&
           //' the ocn co2 flux.  Cannot do both.'
      call endrun(trim(err_str)//errmsg(__FILE__,__LINE__))
   end if
   

end subroutine co2_cycle_readnl

!================================================================================================

subroutine co2_register
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: register advected constituents 
  ! 
  ! Called by:
  !    physpkg.F90
  !-----------------------------------------------------------------------



  use constituents,   only: cnst_add
  use physconst,      only: mwdry, mwco2, gravit, cpair
  
  !local variables
  integer  :: idx
  
  if (.not. co2_flag) return
  
  ! CO2 as dry tracer
  do idx = 1, ncnst
     !mwco2(molecular weight), cpair (heat capacities),1.e-20_r8 (minimum allowed mmr)
     call cnst_add(c_names(idx), mwco2, cpair, 1.e-20_r8, c_i(idx), longname=c_names(idx), mixtype='dry')
     
     select case (trim(c_names(idx)))
     case ('CO2_OCN')
        co2_ocn_glo_ind = c_i(idx)
     case ('CO2_FFF')
        co2_fff_glo_ind = c_i(idx)
     case ('CO2_LND')
        co2_lnd_glo_ind = c_i(idx)
     case ('CO2')
        co2_glo_ind     = c_i(idx)
     end select
     
  end do

end subroutine co2_register

!================================================================================================

function co2_transport()

!-----------------------------------------------------------------------
! 
! Purpose: return true if this package is active
!
! Called by:
!    camsrfexch.F90 atm_import_export.F90 cam_diagnostics.F90 
!    physpkg.F90 restart_physics.F90
!-----------------------------------------------------------------------
   logical :: co2_transport
!-----------------------------------------------------------------------

   co2_transport = co2_flag

end function co2_transport

!================================================================================================

function co2_implements_cnst(name)

!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
!
! Called by:
!    inidat.F90 
!-----------------------------------------------------------------------
    implicit none
!-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name  ! constituent name
    logical :: co2_implements_cnst        ! return value

    integer :: m_ind     
      
    co2_implements_cnst = .false.
 
    if (.not. co2_flag) return
 
    do m_ind = 1, ncnst
       if (name == c_names(m_ind)) then
          co2_implements_cnst = .true.
          return
       end if
    end do
  end function co2_implements_cnst

!===============================================================================  

subroutine co2_init

!----------------------------------------------------------------------- 
! 
! Purpose: initialize co2,
!          declare history variables,
!          read co2 flux form ocn,  as data_flux_ocn
!          read co2 flux form fule, as data_flux_fuel
!
! Called by:
!    physpkg.F90
!-----------------------------------------------------------------------

    use cam_history,    only: addfld, horiz_only, add_default
    use co2_data_flux,  only: co2_data_flux_init
    use constituents,   only: cnst_get_ind, cnst_name, cnst_longname, sflxnam
    !Local variables
    integer :: m_ind, mm


    if (.not. co2_flag) return
 
    ! Add constituents and fluxes to history file
    do m_ind = 1, ncnst
        mm = c_i(m_ind)

       call addfld(trim(cnst_name(mm))//'_BOT',     horiz_only, 'A', 'kg/kg', trim(cnst_longname(mm))//', Bottom Layer')
       call addfld(cnst_name(mm),  (/ 'lev' /), 'A',               'kg/kg', cnst_longname(mm))
       call addfld(sflxnam(mm),   horiz_only, 'A',                 'kg/m2/s', trim(cnst_name(mm))//' surface flux')

       call add_default(cnst_name(mm), 1, ' ')
       call add_default(sflxnam(mm),   1, ' ')

       ! The addfld call for the 'TM*' fields are made by default in the 
       ! constituent_burden module.
       call add_default('TM'//trim(cnst_name(mm)), 1, ' ')
    end do

    ! Options to output both 3D and 2D aircraft emissions (2D is column-integrated 3D field)
    ! Default is only to output the 2D field. -BEH
    call addfld('AF'//trim(cnst_name(c_i(4))), (/ 'lev' /), 'A', 'kg/m2/s', trim(cnst_longname(co2_glo_ind))//' column aircraft flux')
    call addfld('TAF'//trim(cnst_name(c_i(4))), horiz_only, 'A', 'kg/m2/s', trim(cnst_longname(co2_glo_ind))//' column-integrated aircraft flux')
    call add_default('TAF'//trim(cnst_name(c_i(4))), 1, ' ')

 
    ! Read flux data
    if (co2_readFlux_ocn) then
       call co2_data_flux_init ( co2flux_ocn_file,  'CO2_flux', data_flux_ocn )
    end if
 
    if (co2_readFlux_fuel) then
       call co2_data_flux_init ( co2flux_fuel_file, 'CO2_flux', data_flux_fuel )
    end if
 
  end subroutine co2_init

!==========================================================================================

  subroutine co2_time_interp_ocn              
 
!-----------------------------------------------------------------------
!
! Purpose: Time interpolate co2 flux to current time.
!          Read in new monthly data if necessary
!
! Called by:
!    atm_import_export.F90
!-----------------------------------------------------------------------

   use time_manager,   only: is_first_step
   use co2_data_flux,  only: co2_data_flux_advance

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return
 
   if (co2_readFlux_ocn)  then
      call co2_data_flux_advance ( data_flux_ocn )
   endif
 
  end subroutine co2_time_interp_ocn

!===========================================================================================
 
  subroutine co2_time_interp_fuel             
 
!-----------------------------------------------------------------------
!
! Purpose: Time interpolate co2 flux to current time.
!          Read in new monthly data if necessary
!
! Called by:
!    atm_import_export.F90
!-----------------------------------------------------------------------

   use time_manager,   only: is_first_step
   use co2_data_flux,  only: co2_data_flux_advance

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return
 
   if (co2_readFlux_fuel) then
      call co2_data_flux_advance ( data_flux_fuel )
   endif
 
  end subroutine co2_time_interp_fuel

!===========================================================================================

subroutine co2_init_cnst(name, q, gcid)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set initial values of CO2_OCN, CO2_FFF, CO2_LND, CO2
! Need to be called from process_inidat in inidat.F90
! (or, initialize co2 in co2_timestep_init)        
!
! Called by:
!    inidat.F90
!-----------------------------------------------------------------------

  use chem_surfvals,  only: chem_surfvals_get

  ! Arguments
   character(len=*), intent(in) :: name         ! constituent name
   real(r8), intent(out) :: q(:,:)   !  mass mixing ratio
   integer, intent(in) :: gcid(:)    ! global column id
   !-----------------------------------------------------------------------

   if (.not. co2_flag) return
 
   select case (name)
   case ('CO2_OCN')
      q = chem_surfvals_get('CO2MMR')
   case ('CO2_FFF')
      q = chem_surfvals_get('CO2MMR')
   case ('CO2_LND')
      q = chem_surfvals_get('CO2MMR')
   case ('CO2')
      q = chem_surfvals_get('CO2MMR')
   end select

end subroutine co2_init_cnst
!===============================================================================


!===============================================================================

subroutine co2_cycle_set_ptend(state, pbuf, ptend)

!-------------------------------------------------------------------------------
! Purpose:
! Set ptend, using aircraft CO2 emissions in ac_CO2 from pbuf
!
! Called by:
!    physpkg.F90
!-------------------------------------------------------------------------------

   use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
   use constituents,   only: pcnst
   use ppgrid,         only: pver
   use physconst,      only: gravit

   ! Arguments
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_ptend), intent(out)   :: ptend     ! indivdual parameterization tendencies

   ! Local variables
   logical :: lq(pcnst)
   integer :: ifld, ncol, k
   real(r8), pointer :: ac_CO2(:,:)

   !----------------------------------------------------------------------------
   if (.not. co2_flag .or. .not. co2_readFlux_aircraft) then
      call physics_ptend_init(ptend, state%psetcols, 'none')
      return
   end if

   ! aircraft fluxes are added to 'CO2_FFF' and 'CO2' tendencies
   lq(:)               = .false.
   lq(co2_fff_glo_ind) = .true.
   lq(co2_glo_ind)     = .true.

   call physics_ptend_init(ptend, state%psetcols, 'co2_cycle_ac', lq=lq)

   ifld = pbuf_get_index('ac_CO2')   
   call pbuf_get_field(pbuf, ifld, ac_CO2)

   ! [ac_CO2] = 'kg m-2 s-1'
   ! [ptend%q] = 'kg kg-1 s-1'
   ncol = state%ncol
   do k = 1, pver
      ptend%q(:ncol,k,co2_fff_glo_ind) = gravit * state%rpdeldry(:ncol,k) * ac_CO2(:ncol,k)
      ptend%q(:ncol,k,co2_glo_ind)     = gravit * state%rpdeldry(:ncol,k) * ac_CO2(:ncol,k)
   end do

end subroutine co2_cycle_set_ptend

!===============================================================================

subroutine co2_cycle_set_cnst_type(cnst_type_loc, cnst_type_val)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set a local copy of cnst_type to be 'wet' or 'dry'       
!
!-----------------------------------------------------------------------
   use constituents, only: pcnst

! Arguments
   character(len=3), intent(inout) :: cnst_type_loc(pcnst) ! a local copy of cnst_type
   character(len=3), intent(in)    :: cnst_type_val        ! set mmr type: 'wet' or 'dry'
   integer                         :: m                    ! loop index
!-----------------------------------------------------------------------

   if (.not. co2_flag) return

   ! set cnst_type_loc for each CO2 tracer
   do m = 1, ncnst
      cnst_type_loc(c_i(m)) = cnst_type_val
   end do

end subroutine co2_cycle_set_cnst_type
!===============================================================================
 
end module co2_cycle
