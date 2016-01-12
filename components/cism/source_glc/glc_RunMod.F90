!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_RunMod

!BOP
! !MODULE: glc_RunMod

! !DESCRIPTION:
!  Contains the routine for stepping the model forward one timestep
!
! !REVISION HISTORY:
!  SVN:$Id: step_mod.F90 2019 2006-09-29 22:00:15Z njn01 $
!  Adapted by William Lipscomb from step_mod.F90 in POP 2.0 and from 
!   glint_example.F90 in GLIMMER
!
! !USES:

   use glc_kinds_mod
   use glc_time_management, only:  thour, time_manager, check_time_flag, init_time_flag
   use shr_sys_mod
   use glc_communicate, only: my_task, master_task
   use glc_constants, only: verbose, stdout, glc_smb
   use glc_exit_mod, only : exit_glc, sigAbort
   
   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glc_run

!----------------------------------------------------------------------
!
!   module variables
!
!----------------------------------------------------------------------

   integer (i4) ::   &
      cpl_stop_now      ,&! flag id for stop_now flag
      tavg_flag           ! flag to access tavg frequencies

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_run
! !INTERFACE:

 subroutine glc_run(EClock)

! !DESCRIPTION:
!  This routine advances the simulation one timestep.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use glad_main
   use glimmer_log
   use glc_fields 
   use glc_history, only : glc_history_write
   use esmf, only : ESMF_Clock

! !ARGUMENTS:
   type(ESMF_Clock),     intent(in)    :: EClock
   
   
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local or common variables:
!
!-----------------------------------------------------------------------

   logical, save :: first_call = .true.        ! flag for initializing timers

  character(fname_length) ::  & 
     paramfile     ! Name of the top-level configuration file
 
  ! Scalars which hold information about the global grid
 
  integer (i4) ::  &
     nx,ny          ! Size of global glc_grid 

  ! Scalar model outputs
 
  ! Other variables
 
  !TODO - Remove?  Currently not used
  logical ::  &
     ice_tstep    ,&! true if ice timestep was done
     outflag        ! output flag

  integer (i4) ::  & 
     i,j,n          ! indices 

!-----------------------------------------------------------------------
!  things to do on first call
!-----------------------------------------------------------------------

   if (first_call) then
      write(stdout,*) 'In glc_run, first_call =', first_call
      ! this line should set cpl_stop_now = 1 (flag id index)
      cpl_stop_now  = init_time_flag('stop_now',default=.false.)
      tavg_flag     = init_time_flag('tavg')      
      first_call = .false.
   endif

!-----------------------------------------------------------------------
!
!  Take one GLAD time step 
!  Note: For SMB scheme, tsfc = ground surface temperature (Celsius)
!                        qsmb = flux of new glacier ice (kg/m^2s)
!
!        For PDD scheme, tsfc = 2m reference temperature (Celsius)
!                        qsmb = precipitation (kg/m^2/s)
!-----------------------------------------------------------------------

     if (glc_smb) then

         if (verbose .and. my_task==master_task) then 
            write(stdout,*) ' '
            write(stdout,*) 'Call glad, thour =', thour
            write(stdout,*) ' '
         endif

         ! TODO(wjs, 2015-03-23) We will need a loop over instances, either here or
         ! around the call to glc_run
         
         call glad_gcm (params = ice_sheet, instance_index = 1,        &
                        time = nint(thour),                            &
                        qsmb = qsmb, tsfc = tsfc,                      &
                        ice_covered = ice_covered, topo = topo,        &
                        rofi = rofi, rofl = rofl, hflx = hflx,         &
                        ice_sheet_grid_mask=ice_sheet_grid_mask,       &
                        icemask_coupled_fluxes=icemask_coupled_fluxes, &
                        ice_tstep = ice_tstep)

     else    ! use PDD scheme

!TODO - Implement and test PDD option
        call exit_glc(sigAbort, 'ERROR: attempt to use PDD scheme, which has not been implemented')

     endif   ! glc_smb

!-----------------------------------------------------------------------
!
!  update timestep counter, set corresponding model time, set
!  time-dependent logical switches to determine program flow.
!
!-----------------------------------------------------------------------

   call time_manager

   if (verbose .and. my_task==master_task) then
      write(stdout,*) 'Called time manager: new hour =', thour 
   endif

   !-----------------------------------------------------------------------
   ! Write a history file if it's time to do so
   !-----------------------------------------------------------------------

   ! TODO loop over instances
   call glc_history_write(ice_sheet%instances(1), EClock)
   
!-----------------------------------------------------------------------
!EOC

   end subroutine glc_run

!***********************************************************************

 end module glc_RunMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
