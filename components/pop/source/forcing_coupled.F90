!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing_coupled

!BOP
!MODULE: forcing_coupled

! !DESCRIPTION:
!  This module contains all the routines necessary for coupling POP to
!  atmosphere and sea ice models using the NCAR CCSM flux coupler.  To
!  enable the routines in this module, the coupled ifdef option must
!  be specified during the make process.
!
! !REVISION HISTORY:
!  SVN:$Id: forcing_coupled.F90 26603 2011-01-28 23:09:02Z njn01 $
!
! !USES:
 
   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block, get_block
   use domain_size
   use domain
   use io_types, only: stdout, nml_in

   use communicate
   use global_reductions
   use constants
   use io
   use time_management
   use grid
   use prognostic
   use exit_mod
   use ice, only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, tlast_ice
   use forcing_shf
   use forcing_sfwf
   use forcing_ws, only: ws_data_type
   use forcing_fields
   use timers

   !*** ccsm
   use ms_balance
   use tavg
   use registry
   use named_field_mod, only: named_field_register, named_field_get_index, &
       named_field_set, named_field_get
   use forcing_fields
      
   implicit none
   save

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::   &
      coupled_freq_iopt,   &! coupler frequency option
      coupled_freq,        &! frequency of coupling
      ncouple_per_day       ! num of coupler comms per day


#if CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  ids for tavg diagnostics computed from forcing_coupled
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_EVAP_F,       &! tavg id for evaporation flux
      tavg_PREC_F,       &! tavg id for precipitation flux (rain + snow)
      tavg_SNOW_F,       &! tavg id for snow          flux
      tavg_MELT_F,       &! tavg id for melt          flux
      tavg_ROFF_F,       &! tavg id for river runoff  flux
      tavg_IOFF_F,       &! tavg id for ice   runoff  flux due to land-model snow capping
      tavg_SALT_F,       &! tavg id for salt          flux
      tavg_SENH_F,       &! tavg id for sensible heat flux
      tavg_LWUP_F,       &! tavg id for longwave heat flux up
      tavg_LWDN_F,       &! tavg id for longwave heat flux dn
      tavg_MELTH_F,      &! tavg id for melt     heat flux
      tavg_IFRAC          ! tavg id for ice fraction

#endif

!-----------------------------------------------------------------------
!
!  Options for distributing net shortwave heat flux over a coupling
!  interval. All options preserve time-integrated flux.
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      qsw_distrb_iopt_const   = 1, &! qsw constant over a coupling interval
      qsw_distrb_iopt_12hr    = 2, &! qsw smoothly spread over 12 hour window
                                    !    only works for daily coupling
      qsw_distrb_iopt_cosz    = 3   ! qsw proportional to cos of solar zenith angle

   integer (int_kind) :: qsw_distrb_iopt

   real (r8), dimension(:), allocatable ::  &
      qsw_12hr_factor

!-----------------------------------------------------------------------
!  variables for qsw cosz option
!-----------------------------------------------------------------------

   integer (int_kind) :: timer_compute_cosz

   real (r8) ::  &
      tday00_interval_beg,    & ! model time at beginning of coupling interval
      orb_eccen,              & ! Earth eccentricity
      orb_obliqr,             & ! Earth Obliquity
      orb_lambm0,             & ! longitude of perihelion at v-equinox
      orb_mvelpp                ! Earths Moving vernal equinox of orbit +pi

   real (r8), dimension(:,:,:), allocatable :: &
      QSW_COSZ_WGHT,      & ! weights
      QSW_COSZ_WGHT_NORM    ! normalization for QSW_COSZ_WGHT

 
   integer (int_kind), private ::   &
      cpl_ts                ! flag id for coupled_ts flag


!EOC
!***********************************************************************

 contains

!***********************************************************************

!BOP
! !IROUTINE: pop_init_coupled 
! !INTERFACE:

 subroutine pop_init_coupled

! !DESCRIPTION:
!  This routine sets up everything necessary for coupling with CCSM4. 
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) ::  &
      coupled_freq_opt, qsw_distrb_opt

   namelist /coupled_nml/ coupled_freq_opt, coupled_freq, qsw_distrb_opt

   integer (int_kind) ::   &
      k, iblock, nsend,    &
      nml_error             ! namelist i/o error flag

   type (block) ::       &
      this_block          ! block information for current block

!-----------------------------------------------------------------------
!
!  variables associated with qsw 12hr
!
!-----------------------------------------------------------------------

   real (r8) ::  &
      time_for_forcing,   &! time of day for surface forcing
      frac_day_forcing,   &! fraction of day based on time_for_forcing
      cycle_function,     &! intermediate result
      weight_forcing,     &! forcing weights
      sum_forcing          ! sum of forcing weights

   integer (int_kind) ::  &
      count_forcing        ! time step counter (== nsteps_this_interval+1)

   integer (int_kind) ::  &
      i,j,n


!-----------------------------------------------------------------------
!
!  read coupled_nml namelist to start coupling and determine
!  coupling frequency
!
!-----------------------------------------------------------------------
      
   coupled_freq_opt  = 'never'
   coupled_freq_iopt = freq_opt_never
   coupled_freq      = 100000
   qsw_distrb_opt    = 'const'
      
   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=coupled_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading coupled_nml')
   endif

      if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
          write(stdout,*) ' Coupling:'
          write(stdout,blank_fmt)
          write(stdout,*) ' coupled_nml namelist settings:'
          write(stdout,blank_fmt)
          write(stdout, coupled_nml)
          write(stdout,blank_fmt)
      endif

      if (my_task == master_task) then
        select case (coupled_freq_opt)

        case ('nyear')
          coupled_freq_iopt = -1000

        case ('nmonth')
          coupled_freq_iopt = -1000

        case ('nday')
          if (coupled_freq == 1) then
            coupled_freq_iopt = freq_opt_nday
            ncouple_per_day = 1
          else
            coupled_freq_iopt = -1000
          endif

        case ('nhour')
          if (coupled_freq <= 24) then
            coupled_freq_iopt = freq_opt_nhour
            ncouple_per_day = 24/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('nsecond')
          if (coupled_freq <= seconds_in_day) then
            coupled_freq_iopt = freq_opt_nsecond
            ncouple_per_day = seconds_in_day/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('nstep')
          if (coupled_freq <= nsteps_per_day) then
            coupled_freq_iopt = freq_opt_nstep
            ncouple_per_day = nsteps_per_day/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('never')
          coupled_freq_iopt = -9999

        case default
          coupled_freq_iopt = -2000
        end select

        select case (qsw_distrb_opt)
        case ('const')
          qsw_distrb_iopt = qsw_distrb_iopt_const
        case ('12hr')
          qsw_distrb_iopt = qsw_distrb_iopt_12hr
        case ('cosz')
          qsw_distrb_iopt = qsw_distrb_iopt_cosz
          call register_string('qsw_distrb_iopt_cosz')
        case default
          qsw_distrb_iopt = -1000
        end select

      endif
            
      call broadcast_scalar(coupled_freq_iopt, master_task)
      call broadcast_scalar(coupled_freq     , master_task)
      call broadcast_scalar(qsw_distrb_iopt  , master_task)
      call broadcast_scalar(ncouple_per_day  , master_task)

      if (coupled_freq_iopt == -1000) then
        call exit_POP(sigAbort,  &
                 'ERROR: Coupling frequency must be at least once per day')
      else if (coupled_freq_iopt == -2000) then
        call exit_POP(sigAbort,  &
                 'ERROR: Unknown option for coupling frequency')
      endif

      if (registry_match('lcoupled') .eqv. (coupled_freq_iopt == -9999)  ) then
       call exit_POP(sigAbort,  &
       'ERROR: inconsistency between lcoupled and coupled_freq_iopt settings')
      endif

      if (qsw_distrb_iopt == -1000) then
        call exit_POP(sigAbort,  &
                 'ERROR: Unknown option for qsw_distrb_opt')
      endif

!-----------------------------------------------------------------------
!
!  check consistency of the qsw_distrb_iopt option with various
!  time manager options
!
!-----------------------------------------------------------------------

     if ( (qsw_distrb_iopt == qsw_distrb_iopt_12hr) .or. &
           (qsw_distrb_iopt == qsw_distrb_iopt_cosz) ) then
        if ( tmix_iopt /= tmix_avgfit )  &
          call exit_POP(sigAbort,   &
               'ERROR: time_mix_opt must be set to avgfit for qsw_distrb_opt '/&
            &/ 'of 12hr or cosz')

        if ( dttxcel(1) /= c1  .or.  dtuxcel /= c1 )   &
          call exit_POP(sigAbort,   &
               'ERROR: using the specified accelerated integration '/&
            &/ 'technique may not be appropriate for qsw_distrb_opt '/&
            &/ 'of 12hr or cosz')
      endif

!-----------------------------------------------------------------------
!
!  allocate and compute the short wave heat flux multiplier for qsw 12hr
!
!-----------------------------------------------------------------------

      allocate ( qsw_12hr_factor(nsteps_per_interval))
      
      qsw_12hr_factor = c1
      if ( qsw_distrb_iopt == qsw_distrb_iopt_12hr ) then

!       mimic a day

        time_for_forcing = c0 
        count_forcing    =  1
        sum_forcing      = c0

        do n=1,nsteps_per_interval
          frac_day_forcing = time_for_forcing / seconds_in_day 
          cycle_function = cos( pi * ( c2 * frac_day_forcing - c1 ) )
          qsw_12hr_factor(n) = c2 * (      cycle_function      &
                                     + abs(cycle_function) )   &
                                     * cycle_function
          weight_forcing = c1
          if (  count_forcing == 2  .or.             &
            mod(count_forcing,time_mix_freq) == 0 )  &
            weight_forcing = p5
          time_for_forcing = time_for_forcing + weight_forcing * dt(1)
          sum_forcing = sum_forcing   &
                    + weight_forcing * dt(1) * qsw_12hr_factor(n)
          count_forcing = count_forcing + 1
        enddo

        qsw_12hr_factor = qsw_12hr_factor * seconds_in_day   &
                              / sum_forcing

!       check the final integral

        count_forcing =  1
        sum_forcing   = c0

        do n=1,nsteps_per_interval
          weight_forcing = c1
          if (  count_forcing == 2  .or.             &
            mod(count_forcing,time_mix_freq) == 0 )  &
            weight_forcing = p5
          sum_forcing = sum_forcing  &
                    + weight_forcing * dt(1) * qsw_12hr_factor(n)
          count_forcing = count_forcing + 1
        enddo

        if ( sum_forcing < (seconds_in_day - 1.0e-5_r8)  .or.  &
             sum_forcing > (seconds_in_day + 1.0e-5_r8) )      &
          call exit_POP (sigAbort, &
              'ERROR: qsw 12hr temporal integral is incorrect')

      endif

!-----------------------------------------------------------------------
!
!  allocate space for qsw cosz fields
!
!-----------------------------------------------------------------------

      if ( qsw_distrb_iopt == qsw_distrb_iopt_cosz ) then
        allocate( &
          QSW_COSZ_WGHT(nx_block,ny_block,nblocks_clinic), &
          QSW_COSZ_WGHT_NORM(nx_block,ny_block,nblocks_clinic))
      endif

#if CCSMCOUPLED

!-----------------------------------------------------------------------
!
!  define tavg fields computed from forcing_coupled routines
!
!-----------------------------------------------------------------------


   call define_tavg_field(tavg_EVAP_F,'EVAP_F',2,                              &
                          long_name='Evaporation Flux from Coupler',           &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_PREC_F,'PREC_F',2,                              &
                          long_name='Precipitation Flux from Cpl (rain+snow)', &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_SNOW_F,'SNOW_F',2,                              &
                          long_name='Snow Flux from Coupler',                  &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_MELT_F,'MELT_F',2,                              &
                          long_name='Melt Flux from Coupler',                  &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_ROFF_F,'ROFF_F',2,                              &
                          long_name='Runoff Flux from Coupler',                &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_IOFF_F,'IOFF_F',2,                              &
                          long_name='Ice Runoff Flux from Coupler due to Land-Model Snow Capping',            &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_SALT_F,'SALT_F',2,                              &
                          long_name='Salt Flux from Coupler (kg of salt/m^2/s)',&
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_SENH_F,'SENH_F',2,                              &
                          long_name='Sensible Heat Flux from Coupler',         &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_LWUP_F,'LWUP_F',2,                              &
                          long_name='Longwave Heat Flux (up) from Coupler',    &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_LWDN_F,'LWDN_F',2,                              &
                          long_name='Longwave Heat Flux (dn) from Coupler',    &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_MELTH_F,'MELTH_F',2,                            &
                          long_name='Melt Heat Flux from Coupler',             &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_IFRAC,'IFRAC',2,                                &
                          long_name='Ice Fraction from Coupler',               &
                          units='fraction', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')


!-----------------------------------------------------------------------
!
!  Initialize flags and shortwave absorption profile
!  Note that the cpl_write_xxx flags have _no_ default value;
!  therefore, they must be explicitly set .true. and .false.
!  at the appropriate times
!
!-----------------------------------------------------------------------

      call init_time_flag('coupled_ts', cpl_ts,           &
                           owner='pop_init_coupled',      &
                           freq_opt = coupled_freq_iopt,  &
                           freq     = coupled_freq)

!-----------------------------------------------------------------------
!
!  If this is a restart, then read_restart knows the last timestep was
!  a coupled timestep and has registered the string 'coupled_ts_last_true'
!  (read_restart was called prior to the initialization of coupled_ts)
!
!-----------------------------------------------------------------------

      if (registry_match('coupled_ts_last_true') )  &
          call override_time_flag (cpl_ts, old_value=.true.)

      lsmft_avail = .true.


!-----------------------------------------------------------------------
!
!  initialize timer for computing cosz
!
!-----------------------------------------------------------------------
      
   if ( qsw_distrb_iopt == qsw_distrb_iopt_cosz ) then
      call get_timer (timer_compute_cosz, 'COMPUTE_COSZ', nblocks_clinic, &
                                          distrb_clinic%nprocs)
   endif
!-----------------------------------------------------------------------
!
!  register this subroutine
!
!-----------------------------------------------------------------------

   call register_string('pop_init_coupled')


#endif
!-----------------------------------------------------------------------
!EOC

      call flushm (stdout)

 end subroutine pop_init_coupled

!***********************************************************************
!BOP
! !IROUTINE: pop_init_partially_coupled
! !INTERFACE:

 subroutine pop_init_partially_coupled

! !DESCRIPTION:
!  This routine initializes and allocates arrays for the partially-coupled
!  option
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

#if CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) ::   &
      lcoupled

   character (char_len) ::  &
      message

   integer (int_kind) ::   &
      number_of_fatal_errors

   lcoupled = registry_match('lcoupled')

   if ( lcoupled .and. shf_formulation /= 'partially-coupled' ) then
     shf_num_comps = 1
     shf_comp_qsw  = 1

     allocate(SHF_COMP(nx_block,ny_block,max_blocks_clinic,shf_num_comps))
     SHF_COMP = c0
    endif

!-----------------------------------------------------------------------
!
!  initialize and allocate some partially coupled variables
!
!-----------------------------------------------------------------------

   if ( lcoupled                                          &
        .and. sfwf_formulation /= 'partially-coupled'     &
        .and. sfc_layer_type == sfc_layer_varthick .and.  &
        .not. lfw_as_salt_flx .and. liceform ) then

     sfwf_num_comps = 1
     sfwf_comp_cpl  = 1
     tfw_num_comps  = 1
     tfw_comp_cpl   = 1

     allocate(SFWF_COMP(nx_block,ny_block,   max_blocks_clinic,sfwf_num_comps))
     allocate( TFW_COMP(nx_block,ny_block,nt,max_blocks_clinic, tfw_num_comps))

     SFWF_COMP = c0
     TFW_COMP  = c0
   endif

!-----------------------------------------------------------------------
!
!  check compatibility of partially-coupled option with other options
!
!-----------------------------------------------------------------------

   number_of_fatal_errors = 0

   if (.not. lcoupled .and. (shf_formulation  == 'partially-coupled' .or.  &
                             sfwf_formulation == 'partially-coupled' ) ) then
     message =   &
         'ERROR: partially-coupled option is allowed only when coupled'
     write(stdout,*) message
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

   if (lcoupled .and. (shf_formulation  == 'partially-coupled' .and.  &
                      sfwf_formulation /= 'partially-coupled') .or.  &
                      (shf_formulation  /= 'partially-coupled' .and.  &
                       sfwf_formulation == 'partially-coupled') ) then
     message =   &
        'partially-coupled must be used for both shf and sfwf'
     write(stdout,*) message
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

   if (lcoupled .and. shf_formulation /= 'partially-coupled' .and.  &
       shf_data_type /= 'none') then
     message =   &
         'shf_data_type must be set to none or '/&
      &/ 'shf_formulation must be partially_coupled when lcoupled is true'
     write(stdout,*) message
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

   if (lcoupled .and. sfwf_formulation /= 'partially-coupled' .and.  &
       sfwf_data_type /= 'none') then
     message =   &
        'sfwf_data_type must be set to none or '/&
     &/ 'sfwf_formulation must be partially_coupled when lcoupled is true'
     write(stdout,*) message
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!     check coupled compatibility with other forcing options
!
!-----------------------------------------------------------------------

   if (lcoupled .and. ws_data_type /= 'none') then
     message =   &
       'ws_data_type must be set to none in coupled mode'
     write(stdout,*) message
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif
   
   if (number_of_fatal_errors /= 0)  &
      call exit_POP(sigAbort,'subroutine pop_init_partially_coupled')

#endif
!-----------------------------------------------------------------------
!EOC

 call flushm (stdout)

 end subroutine pop_init_partially_coupled


!***********************************************************************

!BOP
! !IROUTINE: pop_set_coupled_forcing
! !INTERFACE:

 subroutine pop_set_coupled_forcing


! !DESCRIPTION:
!  This routine is called immediately following the receipt of fluxes
!  from the coupler. It combines fluxes received from the coupler into 
!  the STF array and converts from W/m**2 into model units. It also 
!  balances salt/freshwater in marginal seas and sets SHF_QSW_RAW 
!  and SHF_COMP. Compute QSW_COSZ_WGHT_NORM if needed.
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC
     
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#if CCSMCOUPLED
   integer (int_kind) :: n, nn, iblock

   real (r8) :: cosz_day ! time where cosz is computed

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
      WORK1, WORK2        ! local work space
 
!-----------------------------------------------------------------------
!
!  combine heat flux components into STF array and convert from W/m**2
!        (note: latent heat flux = evaporation*latent_heat_vapor_mks)
!        (note: snow melt heat flux = - snow_f*latent_heat_fusion_mks)
!
!-----------------------------------------------------------------------

   !*** need to zero out any padded cells
   WORK1 = c0
   WORK2 = c0


   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      STF(:,:,1,iblock) = (EVAP_F(:,:,iblock)*latent_heat_vapor_mks         &
                           + SENH_F(:,:,iblock) + LWUP_F(:,:,iblock)        &
                           + LWDN_F(:,:,iblock) + MELTH_F(:,:,iblock)       &
                           -(SNOW_F(:,:,iblock)+IOFF_F(:,:,iblock)) * latent_heat_fusion_mks)*  &
                             RCALCT(:,:,iblock)*hflux_factor 
   enddo
   !$OMP END PARALLEL DO
                                        
!-----------------------------------------------------------------------
!
!  combine freshwater flux components
!
!  for variable thickness surface layer, compute fresh water and
!  salt fluxes
!
!-----------------------------------------------------------------------

   if (sfc_layer_type == sfc_layer_varthick .and.   &
       .not. lfw_as_salt_flx) then

      !*** compute fresh water flux (cm/s)

      !$OMP PARALLEL DO PRIVATE(iblock,n)
      do iblock = 1, nblocks_clinic

         FW(:,:,iblock) = RCALCT(:,:,iblock) * &
            ( PREC_F(:,:,iblock)+EVAP_F(:,:,iblock)  &
             +ROFF_F(:,:,iblock)+IOFF_F(:,:,iblock))*fwmass_to_fwflux

         WORK1(:,:,iblock) = RCALCT(:,:,iblock) *   &
            MELT_F(:,:,iblock) * fwmass_to_fwflux

        !*** compute tracer concentration in fresh water
        !*** in principle, temperature of each water flux
        !*** could be different. e.g.
        !TFW(:,:,1,iblock) = RCALCT(:,:,iblock)*fwmass_to_fwflux  &
        !                    (PREC_F(:,:,iblock)*TEMP_PREC(:,:,iblock) +  &
        !                     EVAP_F(:,:,iblock)*TEMP_EVAP(:,:,iblock) +  &
        !                     MELT_F(:,:,iblock)*TEMP_MELT(:,:,iblock) +  &
        !                     ROFF_F(:,:,iblock)*TEMP_ROFF(:,:,iblock))
        !*** currently assume water comes in at sea surface temp

         call tmelt(WORK2(:,:,iblock),TRACER(:,:,1,2,curtime,iblock))
         TFW(:,:,1,iblock) = FW(:,:,iblock)*TRACER(:,:,1,1,curtime,iblock)  &
                    + WORK1(:,:,iblock) * WORK2(:,:,iblock)

         FW(:,:,iblock) = FW(:,:,iblock) + WORK1(:,:,iblock)

        !*** compute salt flux
        !*** again, salinity could be different for each
        !***   component of water flux
        !TFW(:,:,2,iblock) = RCALCT(:,:,iblock)*fwmass_to_fwflux  &
        !                    (PREC_F(:,:,iblock)*SALT_PREC(:,:,iblock) +  &
        !                     EVAP_F(:,:,iblock)*SALT_EVAP(:,:,iblock) +  &
        !                     MELT_F(:,:,iblock)*SALT_MELT(:,:,iblock) +  &
        !                     ROFF_F(:,:,iblock)*SALT_ROFF(:,:,iblock))
        !*** currently assume prec, evap and roff are fresh
        !*** and all salt come from ice melt

         where (MELT_F(:,:,iblock) /= c0)
            WORK1(:,:,iblock) =   &
               SALT_F(:,:,iblock)/MELT_F(:,:,iblock) ! salinity (msu) of melt water
         elsewhere
            WORK1(:,:,iblock) = c0
         end where

         TFW(:,:,2,iblock) =  RCALCT(:,:,iblock)*MELT_F(:,:,iblock)*  &
                         fwmass_to_fwflux*WORK1(:,:,iblock)
           ! + PREC_F(:,:,iblock)*c0 + EVAP_F(:,:,iblock)*c0 + ROFF_F(:,:,iblock)*c0 + IOFF_F(:,:,iblock)*c0

         do n=3,nt
            TFW(:,:,n,iblock) = c0  ! no additional tracers in fresh water
         end do

      enddo
      !$OMP END PARALLEL DO

   else  ! convert fresh water to virtual salinity flux

!-----------------------------------------------------------------------
!
!  if not a variable thickness surface layer or if fw_as_salt_flx
!  flag is on, convert fresh and salt inputs to a virtual salinity flux
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1, nblocks_clinic
        STF(:,:,2,iblock) = RCALCT(:,:,iblock)*(  &
                     (PREC_F(:,:,iblock)+EVAP_F(:,:,iblock)+  &
                      MELT_F(:,:,iblock)+ROFF_F(:,:,iblock)+IOFF_F(:,:,iblock))*salinity_factor   &
                    + SALT_F(:,:,iblock)*sflux_factor)  
      enddo
      !$OMP END PARALLEL DO
 
!-----------------------------------------------------------------------
!
!  balance salt/freshwater in marginal seas
!
!-----------------------------------------------------------------------
 
      if  (lms_balance .and. sfwf_formulation /= 'partially-coupled' ) then
       call ms_balancing (STF(:,:,2,:),EVAP_F, PREC_F, MELT_F,ROFF_F,IOFF_F,   &
                          SALT_F, QFLUX, 'salt')
      endif
 
   endif
 

   !$OMP PARALLEL DO PRIVATE(iblock,n)
   do iblock = 1, nblocks_clinic

      SHF_QSW_RAW(:,:,iblock) = SHF_QSW(:,:,iblock)

      if ( shf_formulation == 'partially-coupled' ) then
        SHF_COMP(:,:,iblock,shf_comp_cpl) = STF(:,:,1,iblock) 
        if ( .not. lms_balance ) then
          SHF_COMP(:,:,iblock,shf_comp_cpl) =   &
                  SHF_COMP(:,:,iblock,shf_comp_cpl) * MASK_SR(:,:,iblock)
          SHF_QSW(:,:,iblock) = SHF_QSW(:,:,iblock) * MASK_SR(:,:,iblock)
        endif
      endif
 
      SHF_COMP(:,:,iblock,shf_comp_qsw) = SHF_QSW(:,:,iblock)

      if ( sfwf_formulation == 'partially-coupled' ) then

        if (sfc_layer_type == sfc_layer_varthick .and.  &
            .not. lfw_as_salt_flx) then
          SFWF_COMP(:,:,iblock,sfwf_comp_cpl) =  &
             FW(:,:,iblock) * MASK_SR(:,:,iblock)
             do n=1,nt
              TFW_COMP(:,:,n,iblock,tfw_comp_cpl) = &
                   TFW(:,:,n,iblock) * MASK_SR(:,:,iblock)
             enddo
        else
          SFWF_COMP(:,:,iblock,sfwf_comp_cpl) = &
             STF(:,:,2,iblock) * MASK_SR(:,:,iblock)
        endif

      else

        if ( sfc_layer_type == sfc_layer_varthick .and.  &
             .not. lfw_as_salt_flx .and. liceform ) then
          SFWF_COMP(:,:,iblock,sfwf_comp_cpl)  = FW(:,:,iblock)
          TFW_COMP (:,:,:,iblock,tfw_comp_cpl) = TFW(:,:,:,iblock)
        endif

      endif
 
      if ( luse_cpl_ifrac ) then
        OCN_WGT(:,:,iblock) = (c1-IFRAC(:,:,iblock)) * RCALCT(:,:,iblock)
      endif


    enddo
    !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  Compute QSW_COSZ_WGHT_NORM.
!-----------------------------------------------------------------------

    if ( qsw_distrb_iopt == qsw_distrb_iopt_cosz ) then
       tday00_interval_beg = tday00

       !$OMP PARALLEL DO PRIVATE(iblock,nn,cosz_day)
       do iblock = 1, nblocks_clinic

          QSW_COSZ_WGHT_NORM(:,:,iblock) = c0

          do nn = 1, nsteps_per_interval
             cosz_day = tday00_interval_beg + interval_cum_dayfrac(nn-1) &
                - interval_cum_dayfrac(nsteps_per_interval)

             call compute_cosz(cosz_day, iblock, QSW_COSZ_WGHT(:,:,iblock))

             if (interval_avg_ts(nn)) then
                QSW_COSZ_WGHT_NORM(:,:,iblock) = &
                   QSW_COSZ_WGHT_NORM(:,:,iblock) &
                   + p5 * QSW_COSZ_WGHT(:,:,iblock)
             else
                QSW_COSZ_WGHT_NORM(:,:,iblock) = &
                   QSW_COSZ_WGHT_NORM(:,:,iblock) &
                   + QSW_COSZ_WGHT(:,:,iblock)
             endif

          enddo

          where (QSW_COSZ_WGHT_NORM(:,:,iblock) > c0) &
             QSW_COSZ_WGHT_NORM(:,:,iblock) = &
                (fullsteps_per_interval + p5 * halfsteps_per_interval) &
                / QSW_COSZ_WGHT_NORM(:,:,iblock)

       enddo
       !$OMP END PARALLEL DO
    endif
 
#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine pop_set_coupled_forcing

!***********************************************************************

!BOP
! !IROUTINE: set_combined_forcing
! !INTERFACE:

 subroutine set_combined_forcing (STF,FW,TFW)

! !DESCRIPTION:
!
! This routine combines heat flux components into the STF array and
! converts from W/m**2, then combines terms when the "partially-coupled"
! has been selected
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
 
   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      intent(inout) :: &
      STF,             &! surface tracer fluxes at current timestep
      TFW               !  tracer concentration in water flux

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),  &
      intent(inout) ::  &
      FW                !  fresh water flux


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  integer (int_kind) ::  &
     iblock,             &! local address of current block
     n                    ! index

#if CCSMCOUPLED
  real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
     WORK1, WORK2        ! local work arrays

   !*** need to zero out any padded cells
   WORK1 = c0
   WORK2 = c0

   if ( shf_formulation == 'partially-coupled' ) then
     !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock=1,nblocks_clinic
       STF(:,:,1,iblock) =  SHF_COMP(:,:,iblock,shf_comp_wrest)     &
                   + SHF_COMP(:,:,iblock,shf_comp_srest)     &
                   + SHF_COMP(:,:,iblock,shf_comp_cpl)
     enddo
     !$OMP END PARALLEL DO
   endif

   if ( sfwf_formulation == 'partially-coupled' ) then
     if (sfc_layer_type == sfc_layer_varthick .and.  &
         .not. lfw_as_salt_flx) then
       !$OMP PARALLEL DO PRIVATE(iblock)
       do iblock=1,nblocks_clinic
          STF(:,:,2,iblock) =  SFWF_COMP(:,:,  iblock,sfwf_comp_wrest) &
                             + SFWF_COMP(:,:,  iblock,sfwf_comp_srest)
          FW(:,:,iblock)    =  SFWF_COMP(:,:,  iblock,sfwf_comp_cpl)   &
                             + SFWF_COMP(:,:,  iblock,sfwf_comp_flxio)
          TFW(:,:,:,iblock) =   TFW_COMP(:,:,:,iblock, tfw_comp_cpl)    &
                             +  TFW_COMP(:,:,:,iblock, tfw_comp_flxio)
       enddo
       !$OMP END PARALLEL DO
     else
       if ( lms_balance ) then

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock=1,nblocks_clinic
           WORK1(:,:,iblock) = SFWF_COMP(:,:,iblock,sfwf_comp_flxio) /  &
                                salinity_factor
           WORK2(:,:,iblock) = SFWF_COMP(:,:,iblock,sfwf_comp_cpl)
         enddo
         !$OMP END PARALLEL DO

         call ms_balancing (WORK2, EVAP_F,PREC_F, MELT_F, ROFF_F, IOFF_F,   &
                            SALT_F, QFLUX, 'salt', ICEOCN_F=WORK1)

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock=1,nblocks_clinic
           STF(:,:,2,iblock) =  SFWF_COMP(:,:,iblock,sfwf_comp_wrest)  &
                              + SFWF_COMP(:,:,iblock,sfwf_comp_srest)  &
                              + WORK2(:,:,iblock)                      &
                              + SFWF_COMP(:,:,iblock,sfwf_comp_flxio)* &
                                MASK_SR(:,:,iblock)
         enddo
         !$OMP END PARALLEL DO

       else     

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock=1,nblocks_clinic
           STF(:,:,2,iblock) =  SFWF_COMP(:,:,iblock,sfwf_comp_wrest)  &
                              + SFWF_COMP(:,:,iblock,sfwf_comp_srest)  &
                              + SFWF_COMP(:,:,iblock,sfwf_comp_cpl)    &
                              + SFWF_COMP(:,:,iblock,sfwf_comp_flxio) 
         enddo
         !$OMP END PARALLEL DO
 
       endif
     endif
   endif


#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine set_combined_forcing
 
 
!***********************************************************************
!BOP
! !IROUTINE: tavg_coupled_forcing
! !INTERFACE:

 subroutine tavg_coupled_forcing

! !DESCRIPTION:
!  This routine accumulates tavg diagnostics related to forcing_coupled
!  forcing.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
#if CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

   type (block) ::       &
      this_block          ! block information for current block

   real (r8), dimension(nx_block,ny_block) :: &
      WORK                ! local temp space for tavg diagnostics

!-----------------------------------------------------------------------
!
!  compute and accumulate tavg forcing diagnostics
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

         call accumulate_tavg_field(EVAP_F(:,:,iblock), tavg_EVAP_F,iblock,1)
         call accumulate_tavg_field(PREC_F(:,:,iblock), tavg_PREC_F,iblock,1)
         call accumulate_tavg_field(SNOW_F(:,:,iblock), tavg_SNOW_F,iblock,1)
         call accumulate_tavg_field(MELT_F(:,:,iblock), tavg_MELT_F,iblock,1)
         call accumulate_tavg_field(ROFF_F(:,:,iblock), tavg_ROFF_F,iblock,1)
         call accumulate_tavg_field(IOFF_F(:,:,iblock), tavg_IOFF_F,iblock,1)
         call accumulate_tavg_field(SALT_F(:,:,iblock), tavg_SALT_F,iblock,1)
         call accumulate_tavg_field(SENH_F(:,:,iblock), tavg_SENH_F,iblock,1)
         call accumulate_tavg_field(LWUP_F(:,:,iblock), tavg_LWUP_F,iblock,1)
         call accumulate_tavg_field(LWDN_F(:,:,iblock), tavg_LWDN_F,iblock,1)
         call accumulate_tavg_field(MELTH_F(:,:,iblock),tavg_MELTH_F,iblock,1)
         call accumulate_tavg_field(IFRAC(:,:,iblock),  tavg_IFRAC,iblock,1)
   end do

   !$OMP END PARALLEL DO

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_coupled_forcing

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells_coupler_fluxes
! !INTERFACE:

   subroutine update_ghost_cells_coupler_fluxes(errorCode)

! !DESCRIPTION:
!  This routine accumulates tavg diagnostics related to forcing_coupled
!  forcing.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode  

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  update halos for all coupler fields
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

#if CCSMCOUPLED
   call POP_HaloUpdate(SNOW_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating SNOW_F')
      return
   endif

   call POP_HaloUpdate(PREC_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating PREC_F')
      return
   endif

   call POP_HaloUpdate(EVAP_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating EVAP_F')
      return
   endif

   call POP_HaloUpdate(MELT_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating MELT_F')
      return
   endif

   call POP_HaloUpdate(ROFF_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating ROFF_F')
      return
   endif

   call POP_HaloUpdate(IOFF_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating IOFF_F')
      return
   endif

   call POP_HaloUpdate(SALT_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating SALT_F')
      return
   endif

   call POP_HaloUpdate(SENH_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating SENH_F')
      return
   endif

   call POP_HaloUpdate(LWUP_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating LWUP_F')
      return
   endif

   call POP_HaloUpdate(LWDN_F,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating LWDN_F')
      return
   endif

   call POP_HaloUpdate(MELTH_F,POP_haloClinic,         &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating MELTH_F')
      return
   endif

   call POP_HaloUpdate(SHF_QSW,POP_haloClinic,         &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating SHF_QSW')
      return
   endif

   call POP_HaloUpdate(IFRAC,POP_haloClinic,           &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating IFRAC')
      return
   endif

   call POP_HaloUpdate(ATM_PRESS,POP_haloClinic,       &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating ATM_PRESS')
      return
   endif

   call POP_HaloUpdate(U10_SQR,POP_haloClinic,         &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'update_ghost_cells_coupler: error updating U10_SQR')
      return
   endif

#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine update_ghost_cells_coupler_fluxes

 
!***********************************************************************
!BOP
! !IROUTINE: rotate_wind_stress
! !INTERFACE:

   subroutine rotate_wind_stress (WORK1,WORK2)

! !DESCRIPTION:
!   This subroutine rotates true zonal/meridional wind stress into local
!   coordinates, converts to dyne/cm**2, and shifts SMFT to the U grid
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) ::   &
      WORK1, WORK2        ! contains taux and tauy from coupler


!EOP
!BOC
#if CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (kind=int_kind) :: iblock  
   integer (POP_i4)        :: errorCode  

!-----------------------------------------------------------------------
!
!  rotate and convert
!
!-----------------------------------------------------------------------

   SMFT(:,:,1,:) = (WORK1(:,:,:)*cos(ANGLET(:,:,:)) +  &
                    WORK2(:,:,:)*sin(ANGLET(:,:,:)))*  &
                    RCALCT(:,:,:)*momentum_factor
   SMFT(:,:,2,:) = (WORK2(:,:,:)*cos(ANGLET(:,:,:)) -  &
                    WORK1(:,:,:)*sin(ANGLET(:,:,:)))*  &
                    RCALCT(:,:,:)*momentum_factor

!-----------------------------------------------------------------------
!
!  perform halo updates following the vector rotation
!
!-----------------------------------------------------------------------

   call POP_HaloUpdate(SMFT(:,:,1,:),POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(SMFT(:,:,2,:),POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_POP_r8)

!-----------------------------------------------------------------------
!
!  shift SMFT to U grid
!
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
   call tgrid_to_ugrid(SMF(:,:,1,iblock),SMFT(:,:,1,iblock),iblock)
   call tgrid_to_ugrid(SMF(:,:,2,iblock),SMFT(:,:,2,iblock),iblock)
   enddo ! iblock


#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine rotate_wind_stress

!***********************************************************************
!BOP
! !IROUTINE: compute_cosz
! !INTERFACE:

 subroutine compute_cosz(tday, iblock, COSZ)

! !DESCRIPTION:
!  This subroutine computes cos of the solar zenith angle.
!  Negative values are set to zero.
!
! !REVISION HISTORY:
!  same as module
!
! !USES:

   use shr_orb_mod, only: shr_orb_decl, shr_orb_cosz

! !INPUT PARAMETERS:

   real (r8), intent(in) :: tday
   integer (int_kind), intent(in) :: iblock

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:), intent(out) :: COSZ

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::   &
      i, j            ! loop indices

   real (r8) :: &
      calday,       & ! Calendar day, including fraction
      delta,        & ! Solar declination angle in rad
      eccf            ! Earth-sun distance factor (ie. (1/r)**2)

!-----------------------------------------------------------------------

   call timer_start(timer_compute_cosz, block_id=iblock)

!  shr_orb code assumes Jan 1 = calday 1, unlike Jan 1 = tday 0
   calday = tday + c1

   call shr_orb_decl(calday, orb_eccen, orb_mvelpp, orb_lambm0, &
                     orb_obliqr, delta, eccf)

   do j = 1, ny_block
      do i = 1, nx_block
         COSZ(i,j) = shr_orb_cosz(calday, TLAT(i,j,iblock), &
                                  TLON(i,j,iblock), delta)
         COSZ(i,j) = max(c0, COSZ(i,j))
      enddo
   enddo

   call timer_stop(timer_compute_cosz, block_id=iblock)

!-----------------------------------------------------------------------
!EOC

 end subroutine compute_cosz

!***********************************************************************

 end module forcing_coupled

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

