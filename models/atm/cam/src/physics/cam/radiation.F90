module radiation

!---------------------------------------------------------------------------------
! Purpose:
!
! Provides the CAM interface to the radiation code
!
! Revision history:
! May  2004, D. B. Coleman,  Initial version of interface module.
! July 2004, B. Eaton,       Use interfaces from new shortwave, longwave, and ozone modules.
! Feb  2005, B. Eaton,       Add namelist variables and control of when calcs are done.
! Mar  2008, J. Kay, 	     Add FLDS (downwelling LW rad at surface) as an outfld.
! May  2008, J. Kay, 	     ADD FSUTOA (upwelling SW radiation at TOA) as an outfld.
! Sep  2009, R. Neale, 	     Add FLDSC (clear sky downwelling LW rad at surface) as an outfld.
! Jan  2010, J. Kay          Add COSP simulator calls
! Jan  2010, J. Kay          Add cloud optical depth diagnostics
!---------------------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8
use spmd_utils,      only: masterproc
use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
use physics_types,   only: physics_state, physics_ptend
use physconst,       only: cpair, cappa
use time_manager,    only: get_nstep, is_first_restart_step
use abortutils,      only: endrun
use error_messages,  only: handle_err
use cam_control_mod, only: lambm0, obliqr, mvelpp, eccen
use scamMod,         only: scm_crm_mode, single_column,have_cld,cldobs,&
                           have_clwp,clwpobs,have_tg,tground
use perf_mod,        only: t_startf, t_stopf
use cam_logfile,     only: iulog

implicit none
private
save

public :: &
   radiation_register,    &! registers radiation physics buffer fields
   radiation_defaultopts, &! set default values of namelist variables in runtime_opts
   radiation_setopts,     &! set namelist values from runtime_opts
   radiation_printopts,   &! print namelist values to log
   radiation_get,         &! provide read access to private module data
   radiation_nextsw_cday, &! calendar day of next radiation calculation
   radiation_do,          &! query which radiation calcs are done this timestep
   radiation_init,        &! calls radini
   radiation_tend          ! moved from radctl.F90

integer,public, allocatable :: cosp_cnt(:)       ! counter for cosp
integer,public              :: cosp_cnt_init = 0 !initial value for cosp counter

! Private module data

! Default values for namelist variables

integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).
integer :: iradae = -12    ! frequency of absorp/emis calc in time steps (positive)
                           ! or hours (negative).
integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                           ! or hours (negative) SW/LW radiation will be
                           ! run continuously from the start of an
                           ! initial or restart run
logical :: spectralflux  = .false. ! calculate fluxes (up and down) per band.

!Physics buffer indices
integer :: qrs_idx      = 0
integer :: qrl_idx      = 0 
integer :: cld_idx      = 0
integer :: concld_idx   = 0
integer :: rel_idx      = 0
integer :: rei_idx      = 0
integer :: cicewp_idx = -1
integer :: cliqwp_idx = -1
integer :: cldemis_idx = -1
integer :: cldtau_idx = -1
integer :: nmxrgn_idx = -1
integer :: pmxrgn_idx = -1

!===============================================================================
contains
!===============================================================================

  subroutine radiation_register
!-----------------------------------------------------------------------
! 
! Register radiation fields in the physics buffer
!
!-----------------------------------------------------------------------

    use physics_buffer,  only: pbuf_add_field, dtype_r8

    call pbuf_add_field('QRS' , 'global',dtype_r8,(/pcols,pver/), qrs_idx) ! shortwave radiative heating rate 
    call pbuf_add_field('QRL' , 'global',dtype_r8,(/pcols,pver/), qrl_idx) ! longwave  radiative heating rate 

  end subroutine radiation_register

!================================================================================================

subroutine radiation_defaultopts(iradsw_out, iradlw_out, iradae_out, irad_always_out, spectralflux_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   integer, intent(out), optional :: iradsw_out
   integer, intent(out), optional :: iradlw_out
   integer, intent(out), optional :: iradae_out
   integer, intent(out), optional :: irad_always_out
   logical, intent(out), optional :: spectralflux_out
   !-----------------------------------------------------------------------

   if ( present(iradsw_out) )      iradsw_out = iradsw
   if ( present(iradlw_out) )      iradlw_out = iradlw
   if ( present(iradae_out) )      iradae_out = iradae
   if ( present(irad_always_out) ) irad_always_out = irad_always
   if ( present(spectralflux_out) ) spectralflux_out = spectralflux

end subroutine radiation_defaultopts

!================================================================================================

subroutine radiation_setopts(dtime, nhtfrq, iradsw_in, iradlw_in, iradae_in, &
   irad_always_in, spectralflux_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! *** NOTE *** This routine needs information about dtime (init by dycore) 
!              and nhtfrq (init by history) to do its checking.  Being called
!              from runtime_opts provides these values possibly before they
!              have been set in the modules responsible for them.
!-----------------------------------------------------------------------

   integer, intent(in)           :: dtime           ! timestep size (s)
   integer, intent(in)           :: nhtfrq          ! output frequency of primary history file
   integer, intent(in), optional :: iradsw_in
   integer, intent(in), optional :: iradlw_in
   integer, intent(in), optional :: iradae_in
   integer, intent(in), optional :: irad_always_in
   logical, intent(in), optional :: spectralflux_in

   ! Local
   integer :: ntspdy   ! no. timesteps per day
   integer :: nhtfrq1  ! local copy of input arg nhtfrq
!-----------------------------------------------------------------------

   if ( present(iradsw_in) )      iradsw = iradsw_in
   if ( present(iradlw_in) )      iradlw = iradlw_in
   if ( present(iradae_in) )      iradae = iradae_in
   if ( present(irad_always_in) ) irad_always = irad_always_in
   if ( present(spectralflux_in) ) spectralflux = spectralflux_in

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   ! Convert iradae from hours to timesteps if necessary and check that
   ! iradae must be an even multiple of iradlw
   if (iradae < 0) iradae = nint((-iradae*3600._r8)/dtime)
   if (mod(iradae,iradlw)/=0) then
      write(iulog,*)'radiation_setopts: iradae must be an even multiple of iradlw.'
      write(iulog,*)'     iradae = ',iradae,', iradlw = ',iradlw
      call endrun('radiation_setopts: iradae must be an even multiple of iradlw.')
   end if

   ! Do absorptivities/emissivities have to go on a restart dataset?
   ! The value of nhtfrq from the namelist may need conversion.
   nhtfrq1 = nhtfrq
   if (nhtfrq1 < 0) nhtfrq1 = nint((-nhtfrq1*3600._r8)/dtime)
   ntspdy = nint(86400._r8/dtime)
   if (nhtfrq1 /= 0) then
      if (masterproc .and. mod(nhtfrq1,iradae)/=0) then
         write(iulog,*)'radiation_setopts: *** NOTE: Extra overhead invoked putting',  &
            ' a/e numbers on restart dataset. ***   ',         &
            ' To avoid, make mod(nhtfrq,iradae) = 0'
      end if
   else
      if (masterproc) then
         if (mod(ntspdy,iradae) /= 0 .or. iradae > ntspdy) then
            write(iulog,*)'radiation_setopts: *** NOTE: Extra overhead invoked',  &
               ' putting a/e numbers on restart dataset. ***'
            write(iulog,*)' To avoid, make mod(timesteps per day,iradae)= 0'
         end if
      end if
   end if

end subroutine radiation_setopts

!===============================================================================

subroutine radiation_get(iradsw_out, iradlw_out, iradae_out, irad_always_out, spectralflux_out)
!----------------------------------------------------------------------- 
! Purpose: Provide access to private module data.  (This should be eliminated.)
!-----------------------------------------------------------------------

   integer, intent(out), optional :: iradsw_out
   integer, intent(out), optional :: iradlw_out
   integer, intent(out), optional :: iradae_out
   integer, intent(out), optional :: irad_always_out
   logical, intent(out), optional :: spectralflux_out
   !-----------------------------------------------------------------------

   if ( present(iradsw_out) )      iradsw_out = iradsw
   if ( present(iradlw_out) )      iradlw_out = iradlw
   if ( present(iradae_out) )      iradae_out = iradae
   if ( present(irad_always_out) ) irad_always_out = irad_always
   if ( present(spectralflux_out) ) spectralflux_out = spectralflux_out

end subroutine radiation_get

!================================================================================================

subroutine radiation_printopts
!----------------------------------------------------------------------- 
! Purpose: Print runtime options to log.
!-----------------------------------------------------------------------


   if(irad_always /= 0) write(iulog,10) irad_always
   write(iulog,20) iradsw,iradlw,iradae
10 format(' Execute SW/LW radiation continuously for the first ',i5, ' timestep(s) of this run')
20 format(' Frequency of Shortwave Radiation calc. (IRADSW)     ',i5/, &
          ' Frequency of Longwave Radiation calc. (IRADLW)      ',i5/,  &
          ' Frequency of Absorptivity/Emissivity calc. (IRADAE) ',i5)

end subroutine radiation_printopts

!================================================================================================

function radiation_do(op, timestep)
!----------------------------------------------------------------------- 
! Purpose: Returns true if the specified operation is done this timestep.
!-----------------------------------------------------------------------

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value

   ! Local variables
   integer :: nstep             ! current timestep number
   !-----------------------------------------------------------------------

   if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

   select case (op)

   case ('sw') ! do a shortwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradsw == 1                     &
                    .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('lw') ! do a longwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradlw == 1                     &
                    .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('absems') ! do an absorptivity/emissivity calculation this timestep?
      radiation_do = nstep == 0  .or.  iradae == 1                     &
                    .or. (mod(nstep-1,iradae) == 0  .and.  nstep /= 1)

   case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
      radiation_do = mod(nstep,iradae) /= 0
         
   case default
      call endrun('radiation_do: unknown operation:'//op)

   end select
end function radiation_do

!================================================================================================

real(r8) function radiation_nextsw_cday()
  
!----------------------------------------------------------------------- 
! Purpose: Returns calendar day of next sw radiation calculation
!-----------------------------------------------------------------------

   use time_manager, only: get_curr_calday, get_nstep, get_step_size

   ! Local variables
   integer :: nstep      ! timestep counter
   logical :: dosw       ! true => do shosrtwave calc   
   integer :: offset     ! offset for calendar day calculation
   integer :: dTime      ! integer timestep size
   real(r8):: calday     ! calendar day of 
   !-----------------------------------------------------------------------

   radiation_nextsw_cday = -1._r8
   dosw   = .false.
   nstep  = get_nstep()
   dtime  = get_step_size()
   offset = 0
   do while (.not. dosw)
      nstep = nstep + 1
      offset = offset + dtime
      if (radiation_do('sw', nstep)) then
         radiation_nextsw_cday = get_curr_calday(offset=offset) 
         dosw = .true.
      end if
   end do
   if(radiation_nextsw_cday == -1._r8) then
      call endrun('error in radiation_nextsw_cday')
   end if
        
end function radiation_nextsw_cday

!================================================================================================

  subroutine radiation_init()
!-----------------------------------------------------------------------
!
! Initialize the radiation parameterization, add fields to the history buffer
! 
!-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default, phys_decomp
    use physconst,      only: gravit, cpair, epsilo, stebol, &
                             pstd, mwdry, mwco2, mwo3
    
    use physics_buffer, only: pbuf_get_index
    use radsw,          only: radsw_init
    use radlw,          only: radlw_init
    use radae,          only: radae_init
    use radconstants,   only: radconstants_init
    use radiation_data, only: init_rad_data
    use phys_control,   only: phys_getopts
    use cospsimulator_intr, only: docosp, cospsimulator_intr_init

    integer :: nstep                       ! current timestep number
    logical :: history_amwg                ! output the variables used by the AMWG diag package
    logical :: history_budget              ! output tendencies and state variables for CAM4
                                           ! temperature, water vapor, cloud ice and cloud
                                           ! liquid budgets.
    integer :: history_budget_histfile_num ! output history file number for budget fields
    !-----------------------------------------------------------------------

    call radconstants_init()
    call init_rad_data()

    call radsw_init(gravit)
    call radlw_init(gravit, stebol)
    call radae_init(gravit, epsilo, stebol, pstd, mwdry, mwco2, mwo3)

    ! Get physics buffer indices
    cld_idx    = pbuf_get_index('CLD')
    concld_idx = pbuf_get_index('CONCLD')
    rel_idx    = pbuf_get_index('REL')
    rei_idx    = pbuf_get_index('REI')

    ! "irad_always" is number of time steps to execute radiation continuously from start of
    ! initial OR restart run

    nstep = get_nstep()
    if ( irad_always > 0) then
       nstep       = get_nstep()
       irad_always = irad_always + nstep
    end if
    if (docosp) call cospsimulator_intr_init

    allocate(cosp_cnt(begchunk:endchunk))
    if (is_first_restart_step()) then
      cosp_cnt(begchunk:endchunk)=cosp_cnt_init
    else
      cosp_cnt(begchunk:endchunk)=0 	
    end if

    ! Shortwave radiation
    call addfld ('SOLIN   ','W/m2    ',1,    'A','Solar insolation',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLL    ','W/m2    ',1,    'A','Solar downward near infrared direct  to surface',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('SOLS    ','W/m2    ',1,    'A','Solar downward visible direct  to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLLD   ','W/m2    ',1,    'A','Solar downward near infrared diffuse to surface',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('SOLSD   ','W/m2    ',1,    'A','Solar downward visible diffuse to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRS     ','K/s     ',pver, 'A','Solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRSC    ','K/s     ',pver, 'A','Clearsky solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNS    ','W/m2    ',1,    'A','Net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNT    ','W/m2    ',1,    'A','Net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOA  ','W/m2    ',1,    'A','Net solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSUTOA  ','W/m2    ',1,    'A','Upwelling solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOAC ','W/m2    ',1,    'A','Clearsky net solar flux at top of atmosphere',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('FSDTOA  ','W/m2    ',1,    'A','Downwelling solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSN200  ','W/m2    ',1,    'A','Net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSN200C ','W/m2    ',1,    'A','Clearsky net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTC   ','W/m2    ',1,    'A','Clearsky net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNSC   ','W/m2    ',1,    'A','Clearsky net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSDSC   ','W/m2    ',1,    'A','Clearsky downwelling solar flux at surface',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('FSDS    ','W/m2    ',1,    'A','Downwelling solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUS     ','W/m2    ',pverp,'I','Shortwave upward flux',phys_decomp)
    call addfld ('FDS     ','W/m2    ',pverp,'I','Shortwave downward flux',phys_decomp)
    call addfld ('FUSC    ','W/m2    ',pverp,'I','Shortwave clear-sky upward flux',phys_decomp)
    call addfld ('FDSC    ','W/m2    ',pverp,'I','Shortwave clear-sky downward flux',phys_decomp)
    call addfld ('FSNIRTOA','W/m2    ',1,    'A','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAC','W/m2    ',1,    'A','Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAS','W/m2    ',1,    'A','Net near-infrared flux (>= 0.7 microns) at top of atmosphere',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('SWCF    ','W/m2    ',1,    'A','Shortwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')

    call addfld ('TOT_CLD_VISTAU    ','1',pver,    'A','Total gbx cloud visible sw optical depth',phys_decomp, &
         sampling_seq='rad_lwsw',flag_xyfill=.true.)
    call addfld ('TOT_ICLD_VISTAU   ','1',pver,    'A','Total in-cloud visible sw optical depth',phys_decomp, &
         sampling_seq='rad_lwsw',flag_xyfill=.true.)
    call addfld ('LIQ_ICLD_VISTAU   ','1',pver,    'A','Liquid in-cloud visible sw optical depth',phys_decomp, &
         sampling_seq='rad_lwsw',flag_xyfill=.true.)
    call addfld ('ICE_ICLD_VISTAU   ','1',pver,    'A','Ice in-cloud visible sw optical depth',phys_decomp, &
         sampling_seq='rad_lwsw',flag_xyfill=.true.)

    ! Longwave radiation
    call addfld ('QRL     ','K/s     ',pver, 'A','Longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRLC    ','K/s     ',pver, 'A','Clearsky longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNS    ','W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLDS    ','W/m2    ',1,    'A','Downwelling longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNT    ','W/m2    ',1,    'A','Net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUT    ','W/m2    ',1,    'A','Upwelling longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUTC   ','W/m2    ',1,    'A','Clearsky upwelling longwave flux at top of model',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('FLNTC   ','W/m2    ',1,    'A','Clearsky net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200  ','W/m2    ',1,    'A','Net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200C ','W/m2    ',1,    'A','Clearsky net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNSC   ','W/m2    ',1,    'A','Clearsky net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLDSC   ','W/m2    ',1,    'A','Clearsky downwelling longwave flux at surface',phys_decomp, &
                                                                                       sampling_seq='rad_lwsw')
    call addfld ('LWCF    ','W/m2    ',1,    'A','Longwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUL     ','W/m2    ',pverp,'I','Longwave upward flux',phys_decomp)
    call addfld ('FDL     ','W/m2    ',pverp,'I','Longwave downward flux',phys_decomp)
    call addfld ('FULC    ','W/m2    ',pverp,'I','Longwave clear-sky upward flux',phys_decomp)
    call addfld ('FDLC    ','W/m2    ',pverp,'I','Longwave clear-sky downward flux',phys_decomp)

    ! Heating rate needed for d(theta)/dt computation
    call addfld ('HR      ','K/s     ',pver, 'A','Heating rate needed for d(theta)/dt computation',phys_decomp)
   
   ! determine default variables
    call phys_getopts(history_amwg_out   = history_amwg  , &
                      history_budget_out = history_budget, &
                      history_budget_histfile_num_out = history_budget_histfile_num)

    if (history_amwg) then
       ! Shortwave variables
       call add_default ('SOLIN   ', 1, ' ')
       call add_default ('QRS     ', 1, ' ')
       call add_default ('FSNS    ', 1, ' ')
       call add_default ('FSNT    ', 1, ' ')
       call add_default ('FSDTOA  ', 1, ' ')
       call add_default ('FSNTOA  ', 1, ' ')
       call add_default ('FSUTOA  ', 1, ' ')
       call add_default ('FSNTOAC ', 1, ' ')
       call add_default ('FSNTC   ', 1, ' ')
       call add_default ('FSNSC   ', 1, ' ')
       call add_default ('FSDSC   ', 1, ' ')
       call add_default ('FSDS    ', 1, ' ')
       call add_default ('SWCF    ', 1, ' ')
       ! Longwave variables
       call add_default ('QRL     ', 1, ' ')
       call add_default ('FLNS    ', 1, ' ')
       call add_default ('FLDS    ', 1, ' ')
       call add_default ('FLNT    ', 1, ' ')
       call add_default ('FLUT    ', 1, ' ')
       call add_default ('FLUTC   ', 1, ' ')
       call add_default ('FLNTC   ', 1, ' ')
       call add_default ('FLNSC   ', 1, ' ')
       call add_default ('FLDSC   ', 1, ' ')
       call add_default ('LWCF    ', 1, ' ')
    endif
    if (single_column.and.scm_crm_mode) then
       ! Shortwave variables
       call add_default ('FUS     ', 1, ' ')
       call add_default ('FUSC    ', 1, ' ')
       call add_default ('FDS     ', 1, ' ')
       call add_default ('FDSC    ', 1, ' ')
       ! Longwave variables
       call add_default ('FUL     ', 1, ' ')
       call add_default ('FULC    ', 1, ' ')
       call add_default ('FDL     ', 1, ' ')
       call add_default ('FDLC    ', 1, ' ')
    endif

    if ( history_budget .and. history_budget_histfile_num > 1 ) then
       call add_default ('QRL     ', history_budget_histfile_num, ' ')
       call add_default ('QRS     ', history_budget_histfile_num, ' ')
    end if
    
    cicewp_idx = pbuf_get_index('CICEWP')
    cliqwp_idx = pbuf_get_index('CLIQWP')
    cldemis_idx= pbuf_get_index('CLDEMIS')
    cldtau_idx = pbuf_get_index('CLDTAU')
    nmxrgn_idx = pbuf_get_index('NMXRGN')
    pmxrgn_idx = pbuf_get_index('PMXRGN')

  end subroutine radiation_init

!===============================================================================
  
  subroutine radiation_tend(state,ptend, pbuf, &
       cam_out, cam_in, &
       landfrac,landm,icefrac,snowh, &
       fsns,    fsnt, flns,    flnt,  &
       fsds, net_flx)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Driver for radiation computation.
    ! 
    ! Method: 
    ! Radiation uses cgs units, so conversions must be done from
    ! model fields to radiation fields.
    !
    ! Revision history:
    ! May 2004    D.B. Coleman     Merge of code from radctl.F90 and parts of tphysbc.F90.
    ! 2004-08-09  B. Eaton         Add pointer variables for constituents.
    ! 2004-08-24  B. Eaton         Access O3 and GHG constituents from chem_get_cnst.
    ! 2004-08-30  B. Eaton         Replace chem_get_cnst by rad_constituent_get.
    !-----------------------------------------------------------------------

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    
    use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
    use physics_types,   only: physics_state, physics_ptend
    use cospsimulator_intr, only: docosp, cospsimulator_intr_run, cosp_nradsteps
    use time_manager,    only: get_curr_calday
    use camsrfexch,      only: cam_out_t, cam_in_t    
    use cam_history,     only: outfld
    use cam_history_support, only : fillvalue
    use radheat,         only: radheat_tend
    use ppgrid
    use pspect
    use physconst,        only: cpair, stebol
    use radconstants,     only: nlwbands, nswbands
    use radsw,            only: radcswmx
    use radlw,            only: radclwmx
    use rad_constituents, only: rad_cnst_get_gas, rad_cnst_out
    use aer_rad_props,    only: aer_rad_props_sw, aer_rad_props_lw
    use interpolate_data, only: vertinterp
    use radiation_data,   only: output_rad_data
    use cloud_cover_diags,only: cloud_cover_diags_out


    ! Arguments
    real(r8), intent(in)    :: landfrac(pcols)  ! land fraction
    real(r8), intent(in)    :: landm(pcols)     ! land fraction ramp
    real(r8), intent(in)    :: icefrac(pcols)   ! land fraction
    real(r8), intent(in)    :: snowh(pcols)     ! Snow depth (liquid water equivalent)
    real(r8), intent(inout) :: fsns(pcols)      ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)      ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)      ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)      ! Net outgoing lw flux at model top
    real(r8), intent(inout) :: fsds(pcols)      ! Surface solar down flux
    real(r8), intent(inout) :: net_flx(pcols)

    type(physics_state), intent(in), target :: state
    type(physics_ptend), intent(out)        :: ptend
    
    type(physics_buffer_desc), pointer      :: pbuf(:)
    type(cam_out_t),     intent(inout)      :: cam_out
    type(cam_in_t),      intent(in)         :: cam_in

    ! Local variables

    logical :: dosw, dolw, doabsems
    integer, pointer :: nmxrgn(:)              ! pbuf pointer to Number of maximally overlapped regions
    real(r8),pointer :: pmxrgn(:,:)            ! Maximum values of pressure for each
                                               !    maximally overlapped region.
                                               !    0->pmxrgn(i,1) is range of pressure for
                                               !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
                                               !    2nd region, etc
    real(r8),pointer :: emis(:,:)              ! Cloud longwave emissivity
    real(r8),pointer :: cldtau(:,:)            ! Cloud longwave optical depth
    real(r8),pointer :: cicewp(:,:)            ! in-cloud cloud ice water path
    real(r8),pointer :: cliqwp(:,:)            ! in-cloud cloud liquid water path

    real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
    real(r8) cllow(pcols)                      !       "     low  cloud cover
    real(r8) clmed(pcols)                      !       "     mid  cloud cover
    real(r8) clhgh(pcols)                      !       "     hgh  cloud cover

    real(r8) :: ftem(pcols,pver)               ! Temporary workspace for outfld variables

    integer :: itim
    real(r8), pointer, dimension(:,:) :: rel     ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei     ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction
    real(r8), pointer, dimension(:,:) :: concld     ! convective cloud fraction
    real(r8), pointer, dimension(:,:) :: qrs        ! shortwave radiative heating rate 
    real(r8), pointer, dimension(:,:) :: qrl        ! longwave  radiative heating rate 
    real(r8) :: qrsc(pcols,pver)                    ! clearsky shortwave radiative heating rate 
    real(r8) :: qrlc(pcols,pver)                    ! clearsky longwave  radiative heating rate 

    integer lchnk, ncol
    real(r8) :: calday                        ! current calendar day
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)
    real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
    logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps

    integer i, k                  ! index
    integer :: istat
    real(r8) solin(pcols)         ! Solar incident flux
    real(r8) fsntoa(pcols)        ! Net solar flux at TOA
    real(r8) fsutoa(pcols)        ! Upwelling solar flux at TOA
    real(r8) fsntoac(pcols)       ! Clear sky net solar flux at TOA
    real(r8) fsnirt(pcols)        ! Near-IR flux absorbed at toa
    real(r8) fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
    real(r8) fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8) fsntc(pcols)         ! Clear sky total column abs solar flux
    real(r8) fsdtoa(pcols)        ! Solar input = Flux Solar Downward Top of Atmosphere
    real(r8) fsnsc(pcols)         ! Clear sky surface abs solar flux
    real(r8) fsdsc(pcols)         ! Clear sky surface downwelling solar flux
    real(r8) flut(pcols)          ! Upward flux at top of model
    real(r8) lwcf(pcols)          ! longwave cloud forcing
    real(r8) swcf(pcols)          ! shortwave cloud forcing
    real(r8) tot_cld_vistau(pcols,pver)   ! gbx water+ice cloud optical depth (only during day, night = fillvalue)
    real(r8) tot_icld_vistau(pcols,pver)  ! in-cld water+ice cloud optical depth (only during day, night = fillvalue)
    real(r8) liq_icld_vistau(pcols,pver)  ! in-cld liq cloud optical depth (only during day, night = fillvalue)
    real(r8) ice_icld_vistau(pcols,pver)  ! in-cld ice cloud optical depth (only during day, night = fillvalue)
    real(r8) flutc(pcols)         ! Upward Clear Sky flux at top of model
    real(r8) flntc(pcols)         ! Clear sky lw flux at model top
    real(r8) fldsc(pcols)         ! Clear sky down lw flux at srf 
    real(r8) flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
    real(r8) fln200(pcols)        ! net longwave flux interpolated to 200 mb
    real(r8) fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
    real(r8) flds(pcols)          ! Srf downwelling longwave flux
    real(r8) fns(pcols,pverp)     ! net shortwave flux
    real(r8) fcns(pcols,pverp)    ! net clear-sky shortwave flux
    real(r8) fsn200(pcols)        ! fns interpolated to 200 mb
    real(r8) fsn200c(pcols)       ! fcns interpolated to 200 mb
    real(r8) fnl(pcols,pverp)     ! net longwave flux
    real(r8) fcnl(pcols,pverp)    ! net clear-sky longwave flux

    real(r8) pbr(pcols,pver)      ! Model mid-level pressures (dynes/cm2)
    real(r8) pnm(pcols,pverp)     ! Model interface pressures (dynes/cm2)
    real(r8) eccf                 ! Earth/sun distance factor
    real(r8) lwupcgs(pcols)       ! Upward longwave flux in cgs units

    real(r8), parameter :: cgs2mks = 1.e-3_r8
 
    real(r8), pointer, dimension(:,:) :: n2o    ! nitrous oxide mass mixing ratio
    real(r8), pointer, dimension(:,:) :: ch4    ! methane mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc11  ! cfc11 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc12  ! cfc12 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o3     ! Ozone mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o2     ! Oxygen mass mixing ratio
    real(r8), dimension(pcols) :: o2_col        ! column oxygen mmr
    real(r8), pointer, dimension(:,:) :: co2    ! co2   mass mixing ratio
    real(r8), dimension(pcols) :: co2_col_mean  ! co2 column mean mmr
    real(r8), pointer, dimension(:,:) :: sp_hum ! specific humidity

    ! Aerosol shortwave radiative properties
    real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
    real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
    real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
    real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau

    ! Aerosol longwave absorption optical depth
    real(r8) :: odap_aer(pcols,pver,nlwbands)

    ! Gathered indicies of day and night columns 
    !  chunk_column_index = IdxDay(daylight_column_index)
    integer :: Nday                      ! Number of daylight columns
    integer :: Nnite                     ! Number of night columns
    integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
    integer, dimension(pcols) :: IdxNite ! Indicies of night coumns

    character(*), parameter :: name = 'radiation_tend'
!----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol = state%ncol

    calday = get_curr_calday()

    itim = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, concld_idx, concld, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, qrs_idx,qrs)
    call pbuf_get_field(pbuf, qrl_idx,qrl)

    call pbuf_get_field(pbuf, rel_idx, rel    )
    call pbuf_get_field(pbuf, rei_idx, rei    )
   
    !  For CRM, make cloud equal to input observations:
    if (single_column.and.scm_crm_mode.and.have_cld) then
       do k = 1,pver
          cld(:ncol,k)= cldobs(k)
       enddo
    endif

    ! Cosine solar zenith angle for current time step
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call zenith (calday, clat, clon, coszrs, ncol)

    call output_rad_data(  pbuf, state, cam_in, landm, coszrs )

    ! Gather night/day column indices.
    Nday = 0
    Nnite = 0
    do i = 1, ncol
       if ( coszrs(i) > 0.0_r8 ) then
          Nday = Nday + 1
          IdxDay(Nday) = i
       else
          Nnite = Nnite + 1
          IdxNite(Nnite) = i
       end if
    end do

    dosw     = radiation_do('sw')      ! do shortwave heating calc this timestep?
    dolw     = radiation_do('lw')      ! do longwave heating calc this timestep?
    doabsems = radiation_do('absems')  ! do absorptivity/emissivity calc this timestep?

    if (dosw .or. dolw) then

       ! pbuf cloud properties set in cloud_diagnostics
       call pbuf_get_field(pbuf, cicewp_idx, cicewp)
       call pbuf_get_field(pbuf, cliqwp_idx, cliqwp)
       call pbuf_get_field(pbuf, cldemis_idx, emis)
       call pbuf_get_field(pbuf, cldtau_idx, cldtau)
       call pbuf_get_field(pbuf, pmxrgn_idx, pmxrgn)
       call pbuf_get_field(pbuf, nmxrgn_idx, nmxrgn)

       ! For CRM, make cloud liquid water path equal to input observations
       if(single_column.and.scm_crm_mode.and.have_clwp)then
          do k=1,pver
             cliqwp(:ncol,k) = clwpobs(k)
          end do
       endif

       ! Get specific humidity
       call rad_cnst_get_gas(0,'H2O', state, pbuf,  sp_hum)

       ! Get ozone mass mixing ratio.
       call rad_cnst_get_gas(0,'O3',  state, pbuf,  o3)

       ! Get CO2 mass mixing ratio and compute column mean values
       call rad_cnst_get_gas(0,'CO2', state, pbuf,  co2)
       call calc_col_mean(state, co2, co2_col_mean)

       ! construct cgs unit reps of pmid and pint and get "eccf" - earthsundistancefactor
       call radinp(ncol, state%pmid, state%pint, pbr, pnm, eccf)

       ! Solar radiation computation

       if (dosw) then

          call t_startf('rad_sw')

          ! Get Oxygen mass mixing ratio.
          call rad_cnst_get_gas(0,'O2', state, pbuf,  o2)
          call calc_col_mean(state, o2, o2_col)
   
          ! Get aerosol radiative properties.
          call t_startf('aero_optics_sw')
          call aer_rad_props_sw(0, state, pbuf,  nnite, idxnite, &
               aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)
          call t_stopf('aero_optics_sw')

          call radcswmx(lchnk, &
               ncol,       pnm,        pbr,        sp_hum,     o3,         &
               o2_col,     cld,        cicewp,     cliqwp,     rel,        &
               rei,        eccf,       coszrs,     solin,      &
               cam_in%asdir, cam_in%asdif, cam_in%aldir, cam_in%aldif, nmxrgn, &
               pmxrgn,     qrs,        qrsc,       fsnt,       fsntc,      fsdtoa, &
               fsntoa,     fsutoa,     fsntoac,    fsnirt,     fsnrtc,     fsnirtsq,   &
               fsns,       fsnsc,      fsdsc,      fsds,       cam_out%sols, &
               cam_out%soll, cam_out%solsd, cam_out%solld, fns, fcns,      &
               Nday,       Nnite,      IdxDay,     IdxNite,    co2_col_mean, &
               aer_tau,    aer_tau_w,  aer_tau_w_g, aer_tau_w_f , liq_icld_vistau, ice_icld_vistau  ) 

          call t_stopf('rad_sw')

          !  Output net fluxes at 200 mb
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, fsn200c)
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns, fsn200)

          !
          ! Convert units of shortwave fields needed by rest of model from CGS to MKS
          !
          do i=1,ncol
             solin(i) = solin(i)*cgs2mks
             fsds(i)  = fsds(i)*cgs2mks
             fsnirt(i)= fsnirt(i)*cgs2mks
             fsnrtc(i)= fsnrtc(i)*cgs2mks
             fsnirtsq(i)= fsnirtsq(i)*cgs2mks
             fsnt(i)  = fsnt(i) *cgs2mks
             fsdtoa(i)= fsdtoa(i) *cgs2mks
             fsns(i)  = fsns(i) *cgs2mks
             fsntc(i) = fsntc(i)*cgs2mks
             fsnsc(i) = fsnsc(i)*cgs2mks
             fsdsc(i) = fsdsc(i)*cgs2mks
             fsntoa(i)=fsntoa(i)*cgs2mks
             fsutoa(i)=fsutoa(i)*cgs2mks
             fsntoac(i)=fsntoac(i)*cgs2mks
             fsn200(i)  = fsn200(i)*cgs2mks
             fsn200c(i) = fsn200c(i)*cgs2mks
             swcf(i)=fsntoa(i) - fsntoac(i)
          end do
          ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair


          ! Dump shortwave radiation information to history buffer (diagnostics)
          call outfld('QRS     ',ftem  ,pcols,lchnk)
          ftem(:ncol,:pver) = qrsc(:ncol,:pver)/cpair
          call outfld('QRSC    ',ftem  ,pcols,lchnk)
          call outfld('SOLIN   ',solin ,pcols,lchnk)
          call outfld('FSDS    ',fsds  ,pcols,lchnk)
          call outfld('FSNIRTOA',fsnirt,pcols,lchnk)
          call outfld('FSNRTOAC',fsnrtc,pcols,lchnk)
          call outfld('FSNRTOAS',fsnirtsq,pcols,lchnk)
          call outfld('FSNT    ',fsnt  ,pcols,lchnk)
          call outfld('FSDTOA  ',fsdtoa,pcols,lchnk)
          call outfld('FSNS    ',fsns  ,pcols,lchnk)
          call outfld('FSNTC   ',fsntc ,pcols,lchnk)
          call outfld('FSNSC   ',fsnsc ,pcols,lchnk)
          call outfld('FSDSC   ',fsdsc ,pcols,lchnk)
          call outfld('FSNTOA  ',fsntoa,pcols,lchnk)
          call outfld('FSUTOA  ',fsutoa,pcols,lchnk)
          call outfld('FSNTOAC ',fsntoac,pcols,lchnk)
          call outfld('SOLS    ',cam_out%sols  ,pcols,lchnk)
          call outfld('SOLL    ',cam_out%soll  ,pcols,lchnk)
          call outfld('SOLSD   ',cam_out%solsd ,pcols,lchnk)
          call outfld('SOLLD   ',cam_out%solld ,pcols,lchnk)
          call outfld('FSN200  ',fsn200,pcols,lchnk)
          call outfld('FSN200C ',fsn200c,pcols,lchnk)
          call outfld('SWCF    ',swcf  ,pcols,lchnk)


	  !! initialize tau_cld_vistau and tau_icld_vistau as fillvalue, they will stay fillvalue for night columns
          tot_icld_vistau(1:pcols,1:pver)=fillvalue
          tot_cld_vistau(1:pcols,1:pver)=fillvalue

	  !! only do calcs for tot_cld_vistau and tot_icld_vistau on daytime columns
          do i=1,Nday
	     !! sum the water and ice optical depths to get total in-cloud cloud optical depth
             tot_icld_vistau(IdxDay(i),1:pver)=liq_icld_vistau(IdxDay(i),1:pver)+ice_icld_vistau(IdxDay(i),1:pver)

	     !! sum wat and ice, multiply by cloud fraction to get grid-box value
             tot_cld_vistau(IdxDay(i),1:pver)=(liq_icld_vistau(IdxDay(i),1:pver)+ &
                  ice_icld_vistau(IdxDay(i),1:pver))*cld(IdxDay(i),1:pver)
          end do

	  ! add fillvalue for night columns
          do i = 1, Nnite
              liq_icld_vistau(IdxNite(i),:)  = fillvalue
              ice_icld_vistau(IdxNite(i),:)  = fillvalue
          end do

	  !! use idx_sw_diag, idx_sw_diag is the index to the visible band of the shortwave
          call outfld ('TOT_CLD_VISTAU    ',tot_cld_vistau  ,pcols,lchnk)
          call outfld ('TOT_ICLD_VISTAU   ',tot_icld_vistau ,pcols,lchnk)
          call outfld ('LIQ_ICLD_VISTAU   ',liq_icld_vistau ,pcols,lchnk)
          call outfld ('ICE_ICLD_VISTAU   ',ice_icld_vistau ,pcols,lchnk)
       end if   ! dosw

       ! Longwave radiation computation

       if (dolw) then

          call t_startf("rad_lw")

          ! Convert upward longwave flux units to CGS
          do i=1,ncol
             lwupcgs(i) = cam_in%lwup(i)*1000._r8
             if(single_column.and.scm_crm_mode.and.have_tg) &
                  lwupcgs(i) = 1000*stebol*tground(1)**4
          end do

          ! Get gas phase constituents.
          call rad_cnst_get_gas(0,'N2O',   state, pbuf,  n2o)
          call rad_cnst_get_gas(0,'CH4',   state, pbuf,  ch4)
          call rad_cnst_get_gas(0,'CFC11', state, pbuf,  cfc11)
          call rad_cnst_get_gas(0,'CFC12', state, pbuf,  cfc12)

          ! absems requires lw absorption optical depth and transmission through aerosols
          call t_startf('aero_optics_lw')
          if (doabsems) call aer_rad_props_lw(0, state, pbuf,  odap_aer)
          call t_stopf('aero_optics_lw')

          call radclwmx(lchnk, ncol, doabsems, &
             lwupcgs, state%t, sp_hum, o3, pbr, &
             pnm, state%lnpmid, state%lnpint, n2o, ch4, &
             cfc11, cfc12, cld, emis, pmxrgn, &
             nmxrgn, qrl, qrlc, flns, flnt, flnsc, &
             flntc, cam_out%flwds, fldsc, flut, flutc, &
             fnl, fcnl, co2_col_mean, odap_aer)

          call t_stopf("rad_lw")

          !  Output fluxes at 200 mb
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl, fln200)
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, fln200c)
          !
          ! Convert units of longwave fields needed by rest of model from CGS to MKS
          !
          do i=1,ncol
             flnt(i)  = flnt(i)*cgs2mks
             flut(i)  = flut(i)*cgs2mks
             flutc(i) = flutc(i)*cgs2mks
             flns(i)  = flns(i)*cgs2mks
             flds(i)  = cam_out%flwds(i)*cgs2mks
             fldsc(i) = fldsc(i)*cgs2mks
             flntc(i) = flntc(i)*cgs2mks
             fln200(i)  = fln200(i)*cgs2mks
             fln200c(i) = fln200c(i)*cgs2mks
             flnsc(i) = flnsc(i)*cgs2mks
             cam_out%flwds(i) = cam_out%flwds(i)*cgs2mks
             lwcf(i)=flutc(i) - flut(i)
          end do

          ! Dump longwave radiation information to history tape buffer (diagnostics)
          call outfld('QRL     ',qrl (:ncol,:)/cpair,ncol,lchnk)
          call outfld('QRLC    ',qrlc(:ncol,:)/cpair,ncol,lchnk)
          call outfld('FLNT    ',flnt  ,pcols,lchnk)
          call outfld('FLUT    ',flut  ,pcols,lchnk)
          call outfld('FLUTC   ',flutc ,pcols,lchnk)
          call outfld('FLNTC   ',flntc ,pcols,lchnk)
          call outfld('FLNS    ',flns  ,pcols,lchnk)
          call outfld('FLDS    ',flds  ,pcols,lchnk)
          call outfld('FLNSC   ',flnsc ,pcols,lchnk)
          call outfld('FLDSC   ',fldsc ,pcols,lchnk)
          call outfld('LWCF    ',lwcf  ,pcols,lchnk)
          call outfld('FLN200  ',fln200,pcols,lchnk)
          call outfld('FLN200C ',fln200c,pcols,lchnk)

       end if  !dolw

       ! Output aerosol mmr
       call rad_cnst_out(0, state, pbuf)

       ! Cloud cover diagnostics
       ! radsw can change pmxrgn and nmxrgn so cldsav needs to follow radsw
       call cloud_cover_diags_out(lchnk, ncol, cld, state%pmid, nmxrgn, pmxrgn )

       if (docosp) then
	   !! cosp_cnt referenced for each chunk... cosp_cnt(lchnk)
	   !! advance counter for this timestep
	   cosp_cnt(lchnk) = cosp_cnt(lchnk) + 1

	   !! if counter is the same as cosp_nradsteps, run cosp and reset counter
           if (cosp_nradsteps .eq. cosp_cnt(lchnk)) then
              !call should be compatible with rrtm radiation.F90 interface too, should be with (in),optional
              call cospsimulator_intr_run(state,  pbuf, &
                   cam_in, emis, coszrs, cliqwp_in=cliqwp, cicewp_in=cicewp)
	      cosp_cnt(lchnk) = 0  !! reset counter
           end if
       end if

    else   !  if (dosw .or. dolw) then

       ! convert radiative heating rates from Q*dp to Q for energy conservation
       if (conserve_energy) then
!DIR$ CONCURRENT
          do k =1 , pver
!DIR$ CONCURRENT
             do i = 1, ncol
                qrs(i,k) = qrs(i,k)/state%pdel(i,k)
                qrl(i,k) = qrl(i,k)/state%pdel(i,k)
             end do
          end do
       end if

    end if   !  if (dosw .or. dolw) then

    ! Compute net radiative heating tendency
    call radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                      fsnt, flns, flnt, cam_in%asdir, net_flx)

    ! Compute heating rate for dtheta/dt 
    do k=1,pver
       do i=1,ncol
          ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
       end do
    end do
    call outfld('HR      ',ftem    ,pcols   ,lchnk   )

    ! convert radiative heating rates to Q*dp for energy conservation
    if (conserve_energy) then
!DIR$ CONCURRENT
       do k =1 , pver
!DIR$ CONCURRENT
          do i = 1, ncol
             qrs(i,k) = qrs(i,k)*state%pdel(i,k)
             qrl(i,k) = qrl(i,k)*state%pdel(i,k)
          end do
       end do
    end if
 
    cam_out%netsw(:ncol) = fsns(:ncol)

 end subroutine radiation_tend

!===============================================================================

subroutine radinp(ncol, pmid, pint, pmidrd, pintrd, eccf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
! Convert model pressures to cgs.
! 
! Author: CCM1, CMS Contact J. Kiehl
!-----------------------------------------------------------------------
   use shr_orb_mod
   use time_manager, only: get_curr_calday

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)    ! Pressure at model mid-levels (pascals)
   real(r8), intent(in) :: pint(pcols,pverp)   ! Pressure at model interfaces (pascals)
!
! Output arguments
!
   real(r8), intent(out) :: pmidrd(pcols,pver)  ! Pressure at mid-levels (dynes/cm*2)
   real(r8), intent(out) :: pintrd(pcols,pverp) ! Pressure at interfaces (dynes/cm*2)
   real(r8), intent(out) :: eccf                ! Earth-sun distance factor

!
!---------------------------Local variables-----------------------------
!
   integer i                ! Longitude loop index
   integer k                ! Vertical loop index

   real(r8) :: calday       ! current calendar day
   real(r8) :: delta        ! Solar declination angle
!-----------------------------------------------------------------------
!
   calday = get_curr_calday()
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf)

!
! Convert pressure from pascals to dynes/cm2
!
   do k=1,pver
      do i=1,ncol
         pmidrd(i,k) = pmid(i,k)*10.0_r8
         pintrd(i,k) = pint(i,k)*10.0_r8
      end do
   end do
   do i=1,ncol
      pintrd(i,pverp) = pint(i,pverp)*10.0_r8
   end do

end subroutine radinp

!===============================================================================

subroutine calc_col_mean(state, mmr_pointer, mean_value)
!----------------------------------------------------------------------- 
! 
! Radiation only knows how to work with the column-mean co2 value.
! Compute the column mean.  
!
!-----------------------------------------------------------------------

   use cam_logfile,  only: iulog

   type(physics_state),        intent(in)  :: state
   real(r8), dimension(:,:),   pointer     :: mmr_pointer  ! mass mixing ratio (lev)
   real(r8), dimension(pcols), intent(out) :: mean_value   ! column mean mmr

   integer  :: i, k, ncol
   real(r8) :: ptot(pcols)
   !-----------------------------------------------------------------------

   ncol         = state%ncol
   mean_value   = 0.0_r8
   ptot         = 0.0_r8

   do k=1,pver
      do i=1,ncol
         mean_value(i) = mean_value(i) + mmr_pointer(i,k)*state%pdeldry(i,k)
         ptot(i)         = ptot(i) + state%pdeldry(i,k)
      end do
   end do
   do i=1,ncol
      mean_value(i) = mean_value(i) / ptot(i)
   end do

end subroutine calc_col_mean

!===============================================================================

end module radiation

