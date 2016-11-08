

module co2_cycle

!------------------------------------------------------------------------------------------------
 
! CO2 was used in radiation calculation.
!
! Purpose:
! Provides distributions of CO2_LND, CO2_OCN, CO2_FF, CO2
! Read co2 flux from ocn and fossil fuel.
! Get  co2 flux from lnd through coupler. 
!
! Author: Jeff Lee              
!                                              
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pver
use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
use physconst,      only: mwdry, mwco2, gravit, cpair
use constituents,   only: cnst_add, cnst_get_ind, cnst_name, cnst_longname, sflxnam
use chem_surfvals,  only: chem_surfvals_get
use co2_data_flux
use cam_abortutils,     only: endrun

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

! Public data
 
public data_flux_ocn                 ! data read in for co2 flux from ocn
public data_flux_fuel                ! data read in for co2 flux from fuel
 
TYPE(read_interp) :: data_flux_ocn                  
TYPE(read_interp) :: data_flux_fuel
                         
public c_i                           ! global index for new constituents
public co2_readFlux_ocn              ! read co2 flux from data file 
public co2_readFlux_fuel             ! read co2 flux from data file 


! Namelist variables
logical :: co2_flag            = .false.      ! true => turn on co2 code, namelist variable
logical :: co2_readFlux_ocn    = .false.      ! true => read co2 flux from ocn,  namelist variable
logical :: co2_readFlux_fuel   = .false.      ! true => read co2 flux from fuel, namelist variable
character(len=256) :: co2flux_ocn_file  = 'unset' ! co2 flux from ocn
character(len=256) :: co2flux_fuel_file = 'unset' ! co2 flux from fossil fuel                       

!-----------------------------------------------------------------------
! new constituents
integer, parameter :: ncnst=4                      ! number of constituents implemented

character(len=7), dimension(ncnst), parameter :: & ! constituent names
   c_names = (/'CO2_OCN', 'CO2_FFF', 'CO2_LND', 'CO2    '/)
real(r8), dimension(ncnst), parameter :: &         ! molecular weights
   c_mw = (/mwco2, mwco2, mwco2, mwco2/)
real(r8), dimension(ncnst), parameter :: &         ! heat capacities
   c_cp = (/cpair, cpair, cpair, cpair/)
real(r8), dimension(ncnst), parameter :: &         ! minimum mmr
   c_qmin = (/1.e-20_r8, 1.e-20_r8, 1.e-20_r8, 1.e-20_r8/)
integer, dimension(ncnst) :: c_i                   ! global index

!================================================================================================
contains
!================================================================================================

subroutine co2_cycle_readnl(nlfile)

   ! Read co2_cycle_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=*), parameter :: subname = 'co2_cycle_readnl'

   namelist /co2_cycle_nl/ co2_flag, co2_readFlux_ocn, co2_readFlux_fuel, &
                           co2flux_ocn_file, co2flux_fuel_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'co2_cycle_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, co2_cycle_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast (co2_flag,                               1,   mpilog,  0, mpicom)
   call mpibcast (co2_readFlux_ocn,                       1,   mpilog,  0, mpicom)
   call mpibcast (co2_readFlux_fuel,                      1,   mpilog,  0, mpicom)
   call mpibcast (co2flux_ocn_file,   len(co2flux_ocn_file),   mpichar, 0, mpicom)
   call mpibcast (co2flux_fuel_file, len(co2flux_fuel_file),   mpichar, 0, mpicom)
#endif

end subroutine co2_cycle_readnl

!================================================================================================

subroutine co2_register

!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents 
! 
!-----------------------------------------------------------------------
   integer  :: i

   if (.not. co2_flag) return
 
! CO2 as dry tracer
   do i = 1, ncnst
      call cnst_add(c_names(i), c_mw(i), c_cp(i), c_qmin(i), c_i(i), longname=c_names(i), mixtype='dry')
   end do

end subroutine co2_register

!================================================================================================

function co2_transport()

!-----------------------------------------------------------------------
 
! Purpose: return true if this package is active

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
!-----------------------------------------------------------------------
    implicit none
!-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name  ! constituent name
    logical :: co2_implements_cnst        ! return value

    integer :: m     
      
    co2_implements_cnst = .false.
 
    if (.not. co2_flag) return
 
    do m = 1, ncnst
       if (name == c_names(m)) then
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
!-----------------------------------------------------------------------

    use cam_history, only: addfld, horiz_only, add_default
 
    integer :: m, mm
      
    if (.not. co2_flag) return
 
    ! Add constituents and fluxes to history file
    do m = 1, ncnst
       call cnst_get_ind(c_names(m), mm)

       call addfld(trim(cnst_name(mm))//'_BOT',     horiz_only, 'A', 'kg/kg', trim(cnst_longname(mm))//', Bottom Layer')
       call addfld(cnst_name(mm),  (/ 'lev' /), 'A',               'kg/kg', cnst_longname(mm))
       call addfld(sflxnam(mm),   horiz_only, 'A',                 'kg/m2/s', trim(cnst_name(mm))//' surface flux')

       call add_default(cnst_name(mm), 1, ' ')
       call add_default(sflxnam(mm),   1, ' ')

       ! The addfld call for the 'TM*' fields are made by default in the 
       ! constituent_burden module.
       call add_default('TM'//trim(cnst_name(mm)), 1, ' ')
    end do
 
    ! Read flux data
    if (co2_readFlux_ocn) then
       call read_data_flux ( co2flux_ocn_file,  data_flux_ocn )
    end if
 
    if (co2_readFlux_fuel) then
       call read_data_flux ( co2flux_fuel_file, data_flux_fuel )
    end if
 
  end subroutine co2_init

!==========================================================================================

  subroutine co2_time_interp_ocn              
 
!-----------------------------------------------------------------------
!
! Purpose: Time interpolate co2 flux to current time.
!          Read in new monthly data if necessary
!
!-----------------------------------------------------------------------

   use time_manager,   only: is_first_step

   if (.not. co2_flag) return
 
   if (co2_readFlux_ocn)  then
      if (is_first_step()) then
         call interp_time_flux ( data_flux_ocn,  prev_timestep=.true.  )
      else
         call interp_time_flux ( data_flux_ocn,  prev_timestep=.false. )
      end if
   endif
 
  end subroutine co2_time_interp_ocn

!===========================================================================================
 
  subroutine co2_time_interp_fuel             
 
!-----------------------------------------------------------------------
!
! Purpose: Time interpolate co2 flux to current time.
!          Read in new monthly data if necessary
!
!-----------------------------------------------------------------------

   use time_manager,   only: is_first_step

   if (.not. co2_flag) return
 
   if (co2_readFlux_fuel) then
      if (is_first_step()) then
         call interp_time_flux ( data_flux_fuel, prev_timestep=.true.  )
      else
         call interp_time_flux ( data_flux_fuel, prev_timestep=.false. )
      endif
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
!-----------------------------------------------------------------------
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
 
end module co2_cycle
