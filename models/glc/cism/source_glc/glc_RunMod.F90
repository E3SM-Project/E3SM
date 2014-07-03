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
   use glc_constants, only: verbose, stdout, glc_nec, glc_smb

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

 subroutine glc_run

! !DESCRIPTION:
!  This routine advances the simulation one timestep.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use glint_main
   use glimmer_log
   use glint_global_interp
   use glint_example_clim
   use glc_global_fields 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local or common variables:
!
!-----------------------------------------------------------------------

   logical, save ::    &
      first_call = .true.,        &! flag for initializing timers
      first_global_budget = .true.

  character(fname_length) ::  & 
     paramfile     ! Name of the top-level configuration file
 
  ! Scalars which hold information about the global grid
 
  integer (i4) ::  &
     nx,ny          ! Size of global glc_grid 

  ! Scalar model outputs
 
  !TODO - These are needed only for PDD option (not yet implemented)
  real(r8) ::      & 
     twin         ,&! Timestep-integrated input water flux (kg) 
     twout        ,&! Timestep-integrated output water flux (kg) 
     ice_vol        ! Total ice volume (m^3) 
 
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
!  Take one GLINT time step 
!  Note: For SMB scheme, tsfc = ground surface temperature (Celsius)
!                        qsmb = flux of new glacier ice (kg/m^2s)
!
!        For PDD scheme, tsfc = 2m reference temperature (Celsius)
!                        qsmb = precipitation (kg/m^2/s)
!-----------------------------------------------------------------------

     if (glc_smb) then

         if (verbose .and. my_task==master_task) then 
            write(stdout,*) ' '
            write(stdout,*) 'Call glint, thour =', thour
            write(stdout,*) ' '
             !TODO - Make sure iglint_global and jglint_global are defined appropriately for the global grid
             !      (Currently hardwired in glint_type)
!            write(stdout,*) 'Global fields from CLM to Glint'
!            do n = 1, glc_nec
!               i = iglint_global
!               j = jglint_global   ! N to S global indexing as in Glint
!               write(stdout,*) ' '
!               write(stdout,*) 'i, j, n =', i, j, n
!               write(stdout,*) 'tsfc(n) =', tsfc(i,j,n)
!               write(stdout,*) 'topo(n) =', topo(i,j,n)
!               write(stdout,*) 'qsmb(n) =', qsmb(i,j,n)
!            enddo
         endif

         call glint_gcm (ice_sheet,        nint(thour),     &
                         qsmb,             tsfc,            &
                         topo,                              &
                         ice_tstep = ice_tstep,             & 
                         gfrac = gfrac,    gtopo = gtopo,   &
                         grofi = grofi,    grofl = grofl,   &
                         ghflx = ghflx)

         if (verbose .and. my_task==master_task) then
!            write(stdout,*) ' '
!            write(stdout,*) 'Global fields from GLINT to CLM:'
             !TODO - Make sure iglint_global and jglint_global are defined appropriately for the global grid
!            do n = 1, glc_nec
!               i = iglint_global
!               j = jglint_global   ! N to S global indexing as in GLINT
!               write(stdout,*) ' '
!               write(stdout,*) 'i, j, n =', i, j, n
!               write(stdout,*) 'gfrac(n) =', gfrac(i,j,n)
!               write(stdout,*) 'gtopo(n) =', gtopo(i,j,n)
!               write(stdout,*) 'grofi(n) =', grofi(i,j,n)
!               write(stdout,*) 'grofl(n) =', grofl(i,j,n)
!               write(stdout,*) 'ghflx(n) =', ghflx(i,j,n)
!            enddo
         endif

     else    ! use PDD scheme

!TODO - Implement and test PDD option
         write(stdout,*) 'Using positive-degree-day scheme'
         write(stdout,*) 'WARNING: This has not been tested!'

         call glint (ice_sheet,                  &
                     nint(thour),                &
                     tsfc(:,:,1),                &  ! 2-m air temp
                     qsmb(:,:,1),                &  ! precip
                     orog,                       &
                     output_flag     = outflag,  &
                     ice_frac        = ice_frac, &
                     water_out       = fw,       &
                     water_in        = fw_in,    &
                     total_water_in  = twin,     &
                     total_water_out = twout,    &
                     ice_volume      = ice_vol)

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
!EOC

   end subroutine glc_run

!***********************************************************************

 end module glc_RunMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
