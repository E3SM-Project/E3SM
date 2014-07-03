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
!  SVN:$Id: forcing_coupled.F90 901 2006-05-08 20:47:20Z njn01 $
!
! !USES:
 
   use kinds_mod
   use blocks, only: nx_block, ny_block, block, get_block
   use domain_size
   use domain
   use communicate
   use global_reductions
   use boundary
   use constants
   use io
   use time_management
   use grid
   use prognostic
   use exit_mod
   use ice, only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, tlast_ice
   use forcing_shf
   use forcing_sfwf
   use timers

   !*** ccsm
   use ms_balance
   use tavg
   use cpl_contract_mod
   use cpl_interface_mod
   use cpl_fields_mod
   use registry
   use shr_sys_mod
      
   implicit none
   save

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   logical (log_kind) ::   &
      lcoupled,            &! flag for coupled forcing
      ldiag_cpl = .false. ,&
      lccsm                 ! flag to denote ccsm-specific code


   integer (int_kind) ::   &
      coupled_freq_iopt,   &! coupler frequency option
      coupled_freq          ! frequency of coupling

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
      tavg_SALT_F,       &! tavg id for salt          flux
      tavg_SENH_F,       &! tavg id for sensible heat flux
      tavg_LWUP_F,       &! tavg id for longwave heat flux up
      tavg_LWDN_F,       &! tavg id for longwave heat flux dn
      tavg_MELTH_F        ! tavg id for melt     heat flux

!-----------------------------------------------------------------------
!
!  diurnal cycle switch for the net shortwave heat flux
!
!     qsw_diurnal_cycle = .T.  diurnal cycle is ON
!                       = .F.  diurnal cycle is OFF
!
!-----------------------------------------------------------------------

   logical (log_kind) :: qsw_diurnal_cycle

      real (r8), dimension(:), allocatable ::  &
        diurnal_cycle_factor


#if coupled

   integer (int_kind) ::   &
      timer_send_to_cpl,   &
      timer_recv_from_cpl, &
      timer_recv_to_send,  &
      timer_send_to_recv 
 
   integer (int_kind), private ::   &
      cpl_stop_now,        &! flag id for stop_now flag
      cpl_ts,              &! flag id for coupled_ts flag
      cpl_write_restart,   &! flag id for write restart
      cpl_write_history,   &! flag id for write history
      cpl_write_tavg,      &! flag id for write tavg      
      cpl_diag_global,     &! flag id for computing diagnostics
      cpl_diag_transp       ! flag id for computing diagnostics

   integer (int_kind), dimension(cpl_fields_ibuf_total) ::  &
      isbuf,               &! integer control buffer for sends
      irbuf                 ! integer control buffer for receives
 
   type(cpl_contract) ::  &
      contractS,          &! contract for sends to coupler
      contractR            ! contract for receives from coupler
 
   real (r8), dimension(:,:), allocatable ::  &
      sbuf                 ! temporary send/recv buffer

!-----------------------------------------------------------------------
!  The following variables are used in the exchange of 
!  information between cpl6 and the ocean code.
!
!  ocn --> cpl6
!  ============
!    cpl_fields_ibuf_total -- length of integer ocean "send buffer" vector (isbuf)
!    nsend --  total number of 2D fields sent to cpl6 from ocn
!
!    integer send buffer indices (isbuf in subroutine init_coupled):  
!
!     o  cpl_fields_ibuf_cdate   -- ocean's character date string (yyyymmdd)
!     o  cpl_fields_ibuf_sec     -- ocean's character time string (seconds)
!     o  cpl_fields_ibuf_precadj -- precipitation adjustment factor*1e6
!     o  cpl_fields_ibuf_lsize   -- (iphys_e-iphys_b+1)*(jphys_e-jphys_b+1)
!     o  cpl_fields_ibuf_lisize  -- (iphys_e-iphys_b+1)
!     o  cpl_fields_ibuf_ljsize  -- (jphys_e-jphys_b+1)
!     o  cpl_fields_ibuf_gsize   -- nx_global*ny_global
!     o  cpl_fields_ibuf_gisize  -- nx_global
!     o  cpl_fields_ibuf_gjsize  -- ny_global
!     o  cpl_fields_ibuf_ncpl    -- ncouple_per_day
!     o  cpl_fields_ibuf_nfields -- cpl_fields_grid_total
!     o  cpl_fields_ibuf_dead    --  0 ==>  not a "dead" model
!
!    real send buffer indices (sbuf in subroutine init_coupled):
!
!     o  cpl_fields_grid_lon   -- radian*TLON(i,j)
!     o  cpl_fields_grid_lat   -- radian*TLAT(i,j)
!     o  cpl_fields_grid_area  -- TAREA(i,j)/(radius*radius)
!     o  cpl_fields_grid_mask  -- float(REGION_MASK(i,j))
!     o  cpl_fields_grid_index -- (j_global(j)-1)*(nx_global)+i_global(i)
!
!    real send buffer indices (sbuf in subroutine send_to_coupler):
!
!      o  index_o2c_So_u     -- surface u velocity
!      o  index_o2c_So_v     -- surface v velocity
!      o  index_o2c_So_t     -- surface temperature
!      o  index_o2c_So_s     -- surface salinity
!      o  index_o2c_So_dhdx  -- e,w surface slope
!      o  index_o2c_So_dhdy  -- n,s surface slope
!      o  index_o2c_Fioo_q     -- qflux
!
!
!    cpl6 --> ocn  
!    ============
!
!    cpl_fields_ibuf_total -- length of integer ocean "receive buffer" vector (irbuf)
!
!    integer receive buffer indices (irbuf in subroutine recv_from_coupler):
!
!     o  cpl_fields_ibuf_stopnow  -- stop ocean integration now
!     o  cpl_fields_ibuf_infobug  -- write ocean/coupler diagnostics now  
!     o  cpl_fields_ibuf_resteod  -- write ocean restart files at end of day
!     o  cpl_fields_ibuf_histeod  -- write ocean history files at end of day
!     o  cpl_fields_ibuf_histtavg -- write ocean "tavg"  files at end of day
!     o  cpl_fields_ibuf_diageod  -- write ocean diagnostics   at end of day
!
!    real receive buffer indices (sbuf in subroutine recv_from_coupler):
!
!     o  index_c2o_Foxx_taux   -- zonal wind stress (taux)
!     o  index_c2o_Foxx_tauy   -- meridonal wind stress (tauy)
!     o  index_c2o_Foxx_snow   -- water flux due to snow
!     o  index_c2o_Foxx_rain   -- water flux due to rain
!     o  index_c2o_Foxx_evap   -- evaporation flux
!     o  index_c2o_Foxx_meltw  -- snow melt flux
!     o  index_c2o_Foxx_salt   -- salt
!     o  index_c2o_Foxx_swnet  -- net short-wave heat flux
!     o  index_c2o_Foxx_sen    -- sensible heat flux
!     o  index_c2o_Foxx_lwup   -- longwave radiation (up)
!     o  index_c2o_Foxx_lwdn   -- longwave radiation (down)
!     o  index_c2o_Foxx_melth  -- heat flux from snow&ice melt
!     o  index_c2o_Si_ifrac    -- ice fraction
!     o  index_c2o_Sa_pslv     -- sea-level pressure
!     o  index_c2o_Faoc_duu10n -- 10m wind speed squared
!     o  index_c2o_Forr_roff   -- river runoff flux
!
!
!-----------------------------------------------------------------------

   integer(kind=int_kind) :: index_o2c_So_t        ! temperature
   integer(kind=int_kind) :: index_o2c_So_u        ! velocity, zonal
   integer(kind=int_kind) :: index_o2c_So_v        ! velocity, meridional
   integer(kind=int_kind) :: index_o2c_So_s        ! salinity
   integer(kind=int_kind) :: index_o2c_So_dhdx     ! surface slope, zonal
   integer(kind=int_kind) :: index_o2c_So_dhdy     ! surface slope, meridional
   integer(kind=int_kind) :: index_o2c_Fioo_q      ! heat of fusion (q>0) melt pot (q<0)

   integer(kind=int_kind) :: index_c2o_Si_ifrac    ! state: ice fraction
   integer(kind=int_kind) :: index_c2o_Sa_pslv     ! state: sea level pressure
   integer(kind=int_kind) :: index_c2o_Faoc_duu10n ! state: 10m wind speed squared
   integer(kind=int_kind) :: index_c2o_Foxx_taux   ! wind stress: zonal
   integer(kind=int_kind) :: index_c2o_Foxx_tauy   ! wind stress: meridional
   integer(kind=int_kind) :: index_c2o_Foxx_swnet  ! heat flux: shortwave net
   integer(kind=int_kind) :: index_c2o_Foxx_lat    ! heat flux: latent
   integer(kind=int_kind) :: index_c2o_Foxx_sen    ! heat flux: sensible
   integer(kind=int_kind) :: index_c2o_Foxx_lwup   ! heat flux: long-wave up
   integer(kind=int_kind) :: index_c2o_Foxx_lwdn   ! heat flux: long-wave dow
   integer(kind=int_kind) :: index_c2o_Foxx_melth  ! heat flux: melt
   integer(kind=int_kind) :: index_c2o_Foxx_salt   ! salt flux
   integer(kind=int_kind) :: index_c2o_Foxx_prec   ! water flux: rain+snow
   integer(kind=int_kind) :: index_c2o_Foxx_snow   ! water flux: snow
   integer(kind=int_kind) :: index_c2o_Foxx_rain   ! water flux: rain
   integer(kind=int_kind) :: index_c2o_Foxx_evap   ! water flux: evap
   integer(kind=int_kind) :: index_c2o_Foxx_meltw  ! water flux: melt
   integer(kind=int_kind) :: index_c2o_Forr_roff   ! water flux: runoff

   real (r8) ::  &
      tlast_coupled

   real (r8),   &
      dimension(:,:,:,:), allocatable ::  &
      SBUFF_SUM           ! accumulated sum of send buffer quantities
                          ! for averaging before being sent

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      EVAP_F = c0,       &! evaporation   flux    from cpl (kg/m2/s)
      PREC_F = c0,       &! precipitation flux    from cpl (kg/m2/s)
                          ! (rain + snow)
      SNOW_F = c0,       &! snow          flux    from cpl (kg/m2/s)
      MELT_F = c0,       &! melt          flux    from cpl (kg/m2/s)
      ROFF_F = c0,       &! river runoff  flux    from cpl (kg/m2/s)
      SALT_F = c0,       &! salt          flux    from cpl (kg(salt)/m2/s)
      SENH_F = c0,       &! sensible heat flux    from cpl (W/m2   )
      LWUP_F = c0,       &! longwave heat flux up from cpl (W/m2   )
      LWDN_F = c0,       &! longwave heat flux dn from cpl (W/m2   )
      MELTH_F= c0         ! melt     heat flux    from cpl (W/m2   )
 

#endif
!EOC
!***********************************************************************

 contains

!***********************************************************************

!BOP
! !IROUTINE: init_coupled
! !INTERFACE:

 subroutine init_coupled(SMF, SMFT, STF, SHF_QSW, lsmft_avail)

! !DESCRIPTION:
!  This routine sets up everything necessary for coupling with
!  the CCSM2 flux coupler, version 6 (cpl6)
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic),  &
      intent(out) ::   &
      SMF,             &!  surface momentum fluxes (wind stress)
      SMFT              !  surface momentum fluxes at T points

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic),  &
      intent(out) ::   &
      STF               !  surface tracer fluxes

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),  &
      intent(out) ::  &
      SHF_QSW          !  penetrative solar heat flux

   logical (log_kind), intent(out) ::  &
      lsmft_avail        ! true if SMFT is an available field

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) ::  &
      coupled_freq_opt

   namelist /coupled_nml/ coupled_freq_opt, coupled_freq,  &
                          qsw_diurnal_cycle, lccsm

   integer (int_kind) ::   &
      k, iblock, nsend,    &
      ncouple_per_day,     &! num of coupler comms per day
      nml_error             ! namelist i/o error flag

   type (block) ::       &
      this_block          ! block information for current block

!-----------------------------------------------------------------------
!
!  variables associated with the solar diurnal cycle
!
!-----------------------------------------------------------------------

   real (r8) ::  &
      time_for_forcing,   &! time of day for surface forcing
      frac_day_forcing,   &! fraction of day based on time_for_forcing
      cycle_function,     &! intermediate result of the diurnal cycle function
      weight_forcing,     &! forcing weights
      sum_forcing          ! sum of forcing weights

   integer (int_kind) ::  &
      count_forcing        ! time step counter (== nsteps_this_interval+1)

   integer (int_kind) ::  &
      i,j,n

!-----------------------------------------------------------------------
!  for now:
!  ONLY ALLOW 1 BLOCK PER PROCEESOR
!-----------------------------------------------------------------------

   if (nblocks_clinic /= 1) then
      call exit_POP(sigAbort,'ERROR init_coupled requires nblocks_clinic = 1')
   endif

!-----------------------------------------------------------------------
!
!  read coupled_nml namelist to start coupling and determine
!  coupling frequency
!
!-----------------------------------------------------------------------
      
   lcoupled          = .false.
   lccsm             = .false.
   coupled_freq_opt  = 'never'
   coupled_freq_iopt = freq_opt_never
   coupled_freq      = 100000
   qsw_diurnal_cycle = .false.
      
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
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nday
            ncouple_per_day = 1
          else
            coupled_freq_iopt = -1000
          endif

        case ('nhour')
          if (coupled_freq <= 24) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nhour
            ncouple_per_day = 24/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('nsecond')
          if (coupled_freq <= seconds_in_day) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nsecond
            ncouple_per_day = seconds_in_day/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('nstep')
          if (coupled_freq <= nsteps_per_day) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nstep
            ncouple_per_day = nsteps_per_day/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('never')
          lcoupled = .false.

        case default
          coupled_freq_iopt = -2000
        end select

      endif
            
      call broadcast_scalar(lcoupled,          master_task)
      call broadcast_scalar(lccsm,             master_task)
      call broadcast_scalar(coupled_freq_iopt, master_task)
      call broadcast_scalar(coupled_freq     , master_task)
      call broadcast_scalar(qsw_diurnal_cycle, master_task)

      if (coupled_freq_iopt == -1000) then
        call exit_POP(sigAbort,  &
                 'ERROR: Coupling frequency must be at least once per day')
      else if (coupled_freq_iopt == -2000) then
        call exit_POP(sigAbort,  &
                 'ERROR: Unknown option for coupling frequency')
      endif

!-----------------------------------------------------------------------
!
!     register lcoupled if running with the flux coupler
!
!-----------------------------------------------------------------------

      if (lcoupled) call register_string('lcoupled')
      call register_string('init_coupled')

!-----------------------------------------------------------------------
!
!  check consistency of the qsw_diurnal_cycle option with various
!  time manager options
!
!-----------------------------------------------------------------------

      if ( qsw_diurnal_cycle ) then
        if ( tmix_iopt /= tmix_avgfit )  &
          call exit_POP(sigAbort,   &
               'ERROR: time_mix_opt must be set to avgfit for solar diurnal cycle') 

        if ( dttxcel(1) /= c1  .or.  dtuxcel /= c1 )   &
          call exit_POP(sigAbort,   &
               'ERROR: using the specified accelerated integration '/& 
            &/ ' technique may not be appropriate for solar diurnal cycle')
      endif

!-----------------------------------------------------------------------
!
!  allocate and compute the short wave heat flux multiplier for the 
!  diurnal cycle
!
!-----------------------------------------------------------------------

      allocate ( diurnal_cycle_factor(nsteps_per_interval))
      
      diurnal_cycle_factor = c1
      if ( qsw_diurnal_cycle ) then

!       mimic a day

        time_for_forcing = c0 
        count_forcing    =  1
        sum_forcing      = c0

        do n=1,nsteps_per_interval
          frac_day_forcing = time_for_forcing / seconds_in_day 
          cycle_function = cos( pi * ( c2 * frac_day_forcing - c1 ) )
          diurnal_cycle_factor(n) = c2 * ( cycle_function      &
                                     + abs(cycle_function) )   &
                                     * cycle_function
          weight_forcing = c1
          if (  count_forcing == 2  .or.             &
            mod(count_forcing,time_mix_freq) == 0 )  &
            weight_forcing = p5
          time_for_forcing = time_for_forcing + weight_forcing * dt(1)
          sum_forcing = sum_forcing   &
                    + weight_forcing * dt(1) * diurnal_cycle_factor(n)
          count_forcing = count_forcing + 1
        enddo

        diurnal_cycle_factor = diurnal_cycle_factor * seconds_in_day   &
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
                    + weight_forcing * dt(1) * diurnal_cycle_factor(n)
          count_forcing = count_forcing + 1
        enddo

        if ( sum_forcing < (seconds_in_day - 1.0e-5_r8)  .or.  &
             sum_forcing > (seconds_in_day + 1.0e-5_r8) )      &
          call exit_POP (sigAbort, &
              'ERROR: qsw diurnal cycle temporal integral is incorrect')

      endif
      
#if coupled
      if (.not. lcoupled) then
        call exit_POP(sigAbort,   &
             'ERROR: Coupled ifdef option enabled but lcoupled=false')
      endif


!-----------------------------------------------------------------------
!
!  define tavg fields computed from forcing_coupled routines
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_EVAP_F,'EVAP_F',2,                              &
                          long_name='Evaporation Flux from Coupler',           &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          missing_value=undefined_nf_r4,                       &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_PREC_F,'PREC_F',2,                              &
                          long_name='Precipitation Flux from Cpl (rain+snow)', &
                          missing_value=undefined_nf_r4,                       &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_SNOW_F,'SNOW_F',2,                              &
                          long_name='Snow Flux from Coupler',                  &
                          missing_value=undefined_nf_r4,                       &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_MELT_F,'MELT_F',2,                              &
                          long_name='Melt Flux from Coupler',                  &
                          missing_value=undefined_nf_r4,                       &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_ROFF_F,'ROFF_F',2,                              &
                          long_name='Runoff Flux from Coupler',                &
                          missing_value=undefined_nf_r4,                       &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_SALT_F,'SALT_F',2,                              &
                          long_name='Salt Flux from Coupler (kg of salt/m^2/s)',&
                          missing_value=undefined_nf_r4,                       &
                          units='kg/m^2/s', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_SENH_F,'SENH_F',2,                              &
                          long_name='Sensible Heat Flux from Coupler',         &
                          missing_value=undefined_nf_r4,                       &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_LWUP_F,'LWUP_F',2,                              &
                          long_name='Longwave Heat Flux (up) from Coupler',    &
                          missing_value=undefined_nf_r4,                       &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_LWDN_F,'LWDN_F',2,                              &
                          long_name='Longwave Heat Flux (dn) from Coupler',    &
                          missing_value=undefined_nf_r4,                       &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_MELTH_F,'MELTH_F',2,                            &
                          long_name='Melt Heat Flux from Coupler',             &
                          missing_value=undefined_nf_r4,                       &
                          units='watt/m^2', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')


!-----------------------------------------------------------------------
!
!  Initialize flags and shortwave absorption profile
!  Note that the cpl_write_xxx flags have _no_ default value;
!  therefore, they must be explicitly set .true. and .false.
!  at the appropriate times
!
!-----------------------------------------------------------------------

      cpl_stop_now      = init_time_flag('stop_now',default=.false.)
      cpl_ts            = init_time_flag('coupled_ts',                  &
                                         freq_opt = coupled_freq_iopt,  &
                                         freq     = coupled_freq)
      cpl_write_restart = init_time_flag('cpl_write_restart')
      cpl_write_history = init_time_flag('cpl_write_history')
      cpl_write_tavg    = init_time_flag('cpl_write_tavg'   )
      cpl_diag_global   = init_time_flag('cpl_diag_global')
      cpl_diag_transp   = init_time_flag('cpl_diag_transp')

      lsmft_avail = .true.
      tlast_coupled = c0


!-----------------------------------------------------------------------
!
!   initialize and send buffer
!
!-----------------------------------------------------------------------

      isbuf = 0

      isbuf(cpl_fields_ibuf_cdate  ) = iyear*10000 + imonth*100 + iday
      isbuf(cpl_fields_ibuf_sec    ) =   &
         ihour*seconds_in_hour + iminute*seconds_in_minute + isecond

!maltrud  ASSUME NBLOCKS_CLINIC = 1
      iblock = 1
      this_block = get_block(blocks_clinic(iblock),iblock)

      isbuf(cpl_fields_ibuf_lsize  ) = (this_block%ie-this_block%ib+1)*  &
                                       (this_block%je-this_block%jb+1)
      isbuf(cpl_fields_ibuf_lisize ) = (this_block%ie-this_block%ib+1)
      isbuf(cpl_fields_ibuf_ljsize ) = (this_block%je-this_block%jb+1)
      isbuf(cpl_fields_ibuf_gsize  ) = nx_global*ny_global
      isbuf(cpl_fields_ibuf_gisize ) = nx_global
      isbuf(cpl_fields_ibuf_gjsize ) = ny_global
      isbuf(cpl_fields_ibuf_ncpl   ) = ncouple_per_day
      isbuf(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
      isbuf(cpl_fields_ibuf_dead   ) = 0           ! not a dead model

      allocate(sbuf((this_block%ie-this_block%ib+1)*(this_block%je-this_block%jb+1)  &
      ,        cpl_fields_grid_total))
      sbuf = -888.0
      n=0

   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n=n+1
         sbuf(n,cpl_fields_grid_lon  ) = radian*TLON(i,j,iblock)
         sbuf(n,cpl_fields_grid_lat  ) = radian*TLAT(i,j,iblock)
         sbuf(n,cpl_fields_grid_area ) = TAREA(i,j,iblock)/(radius*radius)
         sbuf(n,cpl_fields_grid_mask ) = float(REGION_MASK(i,j,iblock))
         sbuf(n,cpl_fields_grid_index) =     &
            (this_block%j_glob(j)-1)*(nx_global) + this_block%i_glob(i)
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!  initialize the contracts
!-----------------------------------------------------------------------

      call cpl_interface_contractInit(contractS,cpl_fields_ocnname,  &
         cpl_fields_cplname,cpl_fields_o2c_fields,isbuf,sbuf)
      call cpl_interface_contractInit(contractR,cpl_fields_ocnname,  &
         cpl_fields_cplname,cpl_fields_c2o_fields,isbuf,sbuf)

      write(stdout,*) '(ocn) Initialized contracts with coupler'
      call shr_sys_flush(stdout)

      deallocate(sbuf)

      !--- allocate SBUFF_SUM
      nsend = cpl_interface_contractNumatt(contractS)
      allocate (SBUFF_SUM(nx_block,ny_block,max_blocks_clinic,nsend))

      !--- determine send indices
      index_o2c_So_u     = cpl_interface_contractIndex(contractS,'So_u')
      index_o2c_So_v     = cpl_interface_contractIndex(contractS,'So_v')
      index_o2c_So_t     = cpl_interface_contractIndex(contractS,'So_t')
      index_o2c_So_s     = cpl_interface_contractIndex(contractS,'So_s')
      index_o2c_So_dhdx  = cpl_interface_contractIndex(contractS,'So_dhdx')
      index_o2c_So_dhdy  = cpl_interface_contractIndex(contractS,'So_dhdy')
      index_o2c_Fioo_q   = cpl_interface_contractIndex(contractS,'Fioo_q')

      !--- determine receive indices
      index_c2o_Foxx_taux   = cpl_interface_contractIndex(contractR,'Foxx_taux')
      index_c2o_Foxx_tauy   = cpl_interface_contractIndex(contractR,'Foxx_tauy')
      index_c2o_Foxx_snow   = cpl_interface_contractIndex(contractR,'Foxx_snow')
      index_c2o_Foxx_rain   = cpl_interface_contractIndex(contractR,'Foxx_rain')
      index_c2o_Foxx_evap   = cpl_interface_contractIndex(contractR,'Foxx_evap')
      index_c2o_Foxx_meltw  = cpl_interface_contractIndex(contractR,'Foxx_meltw')
      index_c2o_Foxx_salt   = cpl_interface_contractIndex(contractR,'Foxx_salt')
      index_c2o_Foxx_swnet  = cpl_interface_contractIndex(contractR,'Foxx_swnet')
      index_c2o_Foxx_sen    = cpl_interface_contractIndex(contractR,'Foxx_sen')
      index_c2o_Foxx_lwup   = cpl_interface_contractIndex(contractR,'Foxx_lwup')
      index_c2o_Foxx_lwdn   = cpl_interface_contractIndex(contractR,'Foxx_lwdn')
      index_c2o_Foxx_melth  = cpl_interface_contractIndex(contractR,'Foxx_melth')
      index_c2o_Si_ifrac    = cpl_interface_contractIndex(contractR,'Si_ifrac')
      index_c2o_Sa_pslv     = cpl_interface_contractIndex(contractR,'Sa_pslv')
      index_c2o_Faoc_duu10n = cpl_interface_contractIndex(contractR,'Faoc_duu10n')
      index_c2o_Forr_roff   = cpl_interface_contractIndex(contractR,'Forr_roff')

      !--- receive initial message from coupler
      call cpl_interface_ibufRecv(cpl_fields_cplname,irbuf)

!-----------------------------------------------------------------------
!
!  send initial state info to coupler
!
!-----------------------------------------------------------------------

      call sum_buffer

      call send_to_coupler

!-----------------------------------------------------------------------
!
!  initialize timers for coupled model
!
!-----------------------------------------------------------------------
      
   call get_timer (timer_send_to_cpl  , 'SEND'        , 1, &
                                         distrb_clinic%nprocs)
   call get_timer (timer_recv_from_cpl, 'RECV'        , 1, &
                                         distrb_clinic%nprocs)
   call get_timer (timer_recv_to_send , 'RECV to SEND', 1, &
                                         distrb_clinic%nprocs)
   call get_timer (timer_send_to_recv , 'SEND to RECV', 1, &
                                         distrb_clinic%nprocs)

#endif
!-----------------------------------------------------------------------
!EOC

      call flushm (stdout)

 end subroutine init_coupled

!***********************************************************************

!BOP
! !IROUTINE: set_coupled_forcing
! !INTERFACE:

 subroutine set_coupled_forcing(SMF,SMFT,STF,SHF_QSW,FW,TFW,IFRAC,  &
        ATM_PRESS, U10_SQR)


! !DESCRIPTION:
!  This routine call coupler communication routines to set
!  surface forcing data
!  Note: We are using intent "inout" for SMF,SMFT, STF, SHF_QSW
!        and IFRAC in order to preserve their values inbetween
!        coupling timesteps.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
 
   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic),  &
      intent(inout) ::  &
      SMF,              &!  surface momentum fluxes (wind stress)
      SMFT               !  surface momentum fluxes at T points

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic),  &
         intent(inout) ::  &
      STF,             &!  surface tracer fluxes
      TFW               !  tracer concentration in water flux

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),  &
         intent(inout) ::  &
      SHF_QSW,          &!  penetrative solar heat flux
      FW,               &!  fresh water flux
      IFRAC,            &!  fractional ice coverage
      ATM_PRESS,        &!  atmospheric pressure forcing
      U10_SQR            !  10m wind speed squared

!EOP
!BOC
     
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n, iblock

           
 
#if coupled
!-----------------------------------------------------------------------
!
!     if it is time to couple, exchange data with flux coupler
!     be sure to trigger communication on very first time step
!
!-----------------------------------------------------------------------

   if (nsteps_run /= 0) call sum_buffer

   if (check_time_flag(cpl_ts) .or. nsteps_run == 0) then

     !*** send state variables at end of coupling interval 
 

     if (nsteps_run /= 0) then
        call timer_stop  (timer_recv_to_send)
 
        call timer_start (timer_send_to_cpl)
        call send_to_coupler
        call timer_stop  (timer_send_to_cpl)
     endif

     call timer_start (timer_send_to_recv)
     call timer_stop  (timer_send_to_recv) 
 
     !*** recv data to advance next time step

     call timer_start (timer_recv_from_cpl)
     call recv_from_coupler(SMF,SMFT,STF,SHF_QSW,FW,TFW,IFRAC,ATM_PRESS,  &
                            U10_SQR)
     call timer_stop  (timer_recv_from_cpl)

     call timer_start (timer_recv_to_send)

     !$OMP PARALLEL DO PRIVATE(iblock,n)
 
     do iblock = 1, nblocks_clinic
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
 
   endif

#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine set_coupled_forcing

!***********************************************************************

!BOP
! !IROUTINE: set_combined_forcing
! !INTERFACE:

 subroutine set_combined_forcing (STF,FW,TFW)

! !DESCRIPTION:
!
! This routine combines terms when the "partially-coupled"
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

  real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
     WORK1, WORK2        ! local work arrays

#if coupled

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
       !$OMP PARALLEL DO PRIVATE(iblock,n)
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

         !$OMP PARALLEL DO PRIVATE(iblock,WORK1,WORK2)
         do iblock=1,nblocks_clinic
           WORK1(:,:,iblock) = SFWF_COMP(:,:,iblock,sfwf_comp_flxio) /  &
                                salinity_factor
           WORK2(:,:,iblock) = SFWF_COMP(:,:,iblock,sfwf_comp_cpl)
         enddo
         !$OMP END PARALLEL DO

         call ms_balancing (WORK2, EVAP_F,PREC_F, MELT_F, ROFF_F,    &
                            SALT_F, QFLUX, 'salt', ICEOCN_F=WORK1)

         !$OMP PARALLEL DO PRIVATE(iblock,WORK2)
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
 
#if coupled
!***********************************************************************

!BOP
! !IROUTINE: recv_from_coupler
! !INTERFACE:

 subroutine recv_from_coupler(SMF,SMFT,STF,SHF_QSW,FW,TFW,   &
        IFRAC,ATM_PRESS,U10_SQR)

! !DESCRIPTION:
!  This routine receives message from coupler with surface flux data
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), intent(out) ::  &
      SMF,              &!  surface momentum fluxes (wind stress)
      SMFT               !  surface momentum fluxes at T points

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), intent(out) ::  &
      STF,             &!  surface tracer fluxes
      TFW               !  tracer concentration in water flux

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(out) ::  &
      SHF_QSW,          &!  penetrative solar heat flux
      FW,               &!  fresh water flux
      IFRAC,            &!  fractional ice coverage
      ATM_PRESS,        &!  atmospheric pressure forcing
      U10_SQR            !  10m wind speed squared

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len)    :: label
 
   integer (int_kind) ::  &
      nrecv,              &
      i,j,k,n,iblock

   real (r8), dimension(nx_block,ny_block) ::  &
      WORKB

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
      WORK1, WORK2,      &! local work space
      WORK3, WORK4,      &! local work space
      WORK5               ! local work space

   real (r8) ::  &
      m2percm2,  &
      gsum

   type (block) :: this_block ! local block info

!-----------------------------------------------------------------------
!
!  receive message from coupler and check for terminate signal
!
!-----------------------------------------------------------------------

   nrecv = cpl_interface_contractNumatt(contractR)

!maltrud  ASSUME NBLOCKS_CLINIC = 1
   iblock = 1
   this_block = get_block(blocks_clinic(iblock),iblock)

   allocate(sbuf((this_block%ie-this_block%ib+1)*(this_block%je-this_block%jb+1), &
      nrecv))
   call cpl_interface_contractRecv(cpl_fields_cplname,contractR, &
      irbuf,sbuf)

!-----------------------------------------------------------------------
!
!  check all coupler flags and respond appropriately
!
!-----------------------------------------------------------------------

   if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
     call set_time_flag(cpl_stop_now,.true.)
     if (my_task == master_task) then
       call int_to_char (4,iyear   , cyear  )
       call int_to_char (2,imonth  , cmonth )
       call int_to_char (2,iday    , cday   )
       call int_to_char (2,ihour   , chour  )
       call int_to_char (2,iminute , cminute)
       call int_to_char (2,isecond , csecond)
       write(stdout,*) '(recv_from_coupler) ',  &
        'cpl requests termination now: ', &
        cyear,'/',cmonth,'/',cday,' ', chour,':',cminute,':',csecond
     endif
     RETURN
   endif


   if (irbuf(cpl_fields_ibuf_infobug) >= 2) then
     ldiag_cpl = .true. 
   else
     ldiag_cpl = .false.
   endif

   if (irbuf(cpl_fields_ibuf_resteod) == 1) then
     call set_time_flag(cpl_write_restart,.true.)
     if (my_task == master_task) then
       write(stdout,*) '(recv_from_coupler) ', &
         'cpl requests restart file at eod  ',cyear,'/',cmonth,'/',cday
     endif
   endif
  
!   if (irbuf(cpl_fields_ibuf_histeod) == 1) then
!    ignore for now
!     call set_time_flag(cpl_write_history,.true.)
!     call int_to_char (4,iyear   , cyear )
!     call int_to_char (2,imonth  ,cmonth )
!     call int_to_char (2,iday    ,cday   )
!     call int_to_char (2,ihour   ,chour  )
!     call int_to_char (2,iminute ,cminute)
!     call int_to_char (2,isecond ,csecond)
!     if (my_task == master_task) then
!     write(stdout,*) ' cpl requests history file at eod '  &
!  ,                     ' ', cyear,'/',cmonth,'/',cday, '  '
!     endif
!   endif
  
  
!   if (irbuf(cpl_fields_ibuf_histtavg) == 1) then
!    ignore for now
!     call set_time_flag(cpl_write_tavg, .true.)
!     call int_to_char (4,iyear   , cyear )
!     call int_to_char (2,imonth  ,cmonth )
!     call int_to_char (2,iday    ,cday   )
!     call int_to_char (2,ihour   ,chour  )
!     call int_to_char (2,iminute ,cminute)
!     call int_to_char (2,isecond ,csecond)
!     if (my_task == master_task) then
!     write(stdout,*) ' cpl requests tavg file at eod '  &
!  ,                     ' ', cyear,'/',cmonth,'/',cday, '  '
!     endif
!   endif
  
    if (irbuf(cpl_fields_ibuf_diageod) == 1) then
      call set_time_flag(cpl_diag_global,.true.)
      call set_time_flag(cpl_diag_transp,.true.)

      call int_to_char (4,iyear   ,cyear  )
      call int_to_char (2,imonth  ,cmonth )
      call int_to_char (2,iday    ,cday   )
      call int_to_char (2,ihour   ,chour  )
      call int_to_char (2,iminute ,cminute)
      call int_to_char (2,isecond ,csecond)

      if (my_task == master_task) then
        write(stdout,*) ' cpl requests diagnostics at eod ' , &
                        ' ', cyear,'/',cmonth,'/',cday, '  '
      endif
    endif

!-----------------------------------------------------------------------
!
!  unpack and distribute wind stress, then convert to correct units
!  and rotate components to local coordinates
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         WORK1(i,j,iblock) = sbuf(n,index_c2o_Foxx_taux)
         WORK2(i,j,iblock) = sbuf(n,index_c2o_Foxx_tauy)
      enddo
      enddo
   enddo

   !***
   !*** do boundary updates now to ensure correct T->U grid
   !***

   call update_ghost_cells(WORK1, bndy_clinic, &
                           field_loc_center, field_type_vector)
   call update_ghost_cells(WORK2, bndy_clinic, &
                           field_loc_center, field_type_vector)

   n = 0
   do iblock = 1, nblocks_clinic

      !***
      !*** Rotate true zonal/meridional wind stress into local
      !*** coordinates and convert to dyne/cm**2
      !***

      SMFT(:,:,1,iblock) = (WORK1(:,:,iblock)*cos(ANGLET(:,:,iblock)) +   &
                            WORK2(:,:,iblock)*sin(ANGLET(:,:,iblock)))*  &
                           RCALCT(:,:,iblock)*momentum_factor
      SMFT(:,:,2,iblock) = (WORK2(:,:,iblock)*cos(ANGLET(:,:,iblock)) -   &
                            WORK1(:,:,iblock)*sin(ANGLET(:,:,iblock)))*  &
                           RCALCT(:,:,iblock)*momentum_factor
 
      !***
      !*** Shift SMFT to U grid
      !***

      call tgrid_to_ugrid(SMF(:,:,1,iblock),SMFT(:,:,1,iblock),iblock)
      call tgrid_to_ugrid(SMF(:,:,2,iblock),SMFT(:,:,2,iblock),iblock)

!-----------------------------------------------------------------------
!
!  unpack and distribute fresh water flux and salt flux
!
!-----------------------------------------------------------------------

      this_block = get_block(blocks_clinic(iblock),iblock)

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         SNOW_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_snow)
         WORKB (i,j       ) = sbuf(n,index_c2o_Foxx_rain)
         EVAP_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_evap)
         MELT_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_meltw)
         ROFF_F(i,j,iblock) = sbuf(n,index_c2o_Forr_roff)
         SALT_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_salt)

         PREC_F(i,j,iblock) = WORKB(i,j) + SNOW_F(i,j,iblock)    ! rain + snow

         WORKB(i,j        ) = sbuf(n,index_c2o_Foxx_swnet)
         SHF_QSW(i,j,iblock) = WORKB(i,j)*  &
            RCALCT(i,j,iblock)*hflux_factor  !  convert from W/m**2

         SENH_F(i,j,iblock)  = sbuf(n,index_c2o_Foxx_sen)
         LWUP_F(i,j,iblock)  = sbuf(n,index_c2o_Foxx_lwup)
         LWDN_F(i,j,iblock)  = sbuf(n,index_c2o_Foxx_lwdn)
         MELTH_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_melth)

         WORKB(i,j       ) = sbuf(n,index_c2o_Si_ifrac)
         IFRAC(i,j,iblock) = WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from Pa to dynes/cm**2
         WORKB(i,j       ) = sbuf(n,index_c2o_Sa_pslv)
         ATM_PRESS(i,j,iblock) = c10 * WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from m**2/s**2 to cm**2/s**2
         WORKB(i,j       ) = sbuf(n,index_c2o_Faoc_duu10n)
         U10_SQR(i,j,iblock) = cmperm * cmperm * WORKB(i,j) * RCALCT(i,j,iblock)

      enddo
      enddo

   enddo

   call update_ghost_cells(SNOW_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(PREC_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(EVAP_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(MELT_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(ROFF_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(SALT_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)

   call update_ghost_cells(SENH_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(LWUP_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(LWDN_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(MELTH_F, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(SHF_QSW, bndy_clinic, &
                           field_loc_center, field_type_scalar)

   call update_ghost_cells(IFRAC, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(ATM_PRESS, bndy_clinic, &
                           field_loc_center, field_type_scalar)
   call update_ghost_cells(U10_SQR, bndy_clinic, &
                           field_loc_center, field_type_scalar)

!-----------------------------------------------------------------------
!
!  combine heat flux components into STF array and convert from W/m**2
!        (note: latent heat flux = evaporation*latent_heat_vapor)
!        (note: snow melt heat flux = - snow_f*latent_heat_fusion_mks)
!
!-----------------------------------------------------------------------


   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      STF(:,:,1,iblock) = (EVAP_F(:,:,iblock)*latent_heat_vapor             &
                           + SENH_F(:,:,iblock) + LWUP_F(:,:,iblock)        &
                           + LWDN_F(:,:,iblock) + MELTH_F(:,:,iblock)       &
                           - SNOW_F(:,:,iblock) * latent_heat_fusion_mks)*  &
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
             +ROFF_F(:,:,iblock))*fwmass_to_fwflux

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
           ! + PREC_F(:,:,iblock)*c0 + EVAP_F(:,:,iblock)*c0 + ROFF_F(:,:,iblock)*c0

         do n=3,nt
            TFW(:,:,n,iblock) = c0  ! no additional tracers in fresh water
         end do

      enddo
      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  if not a variable thickness surface layer or if fw_as_salt_flx
!  flag is on, convert fresh and salt inputs to a virtual salinity flux
!
!-----------------------------------------------------------------------

   else  ! convert fresh water to virtual salinity flux

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1, nblocks_clinic
        STF(:,:,2,iblock) = RCALCT(:,:,iblock)*(  &
                     (PREC_F(:,:,iblock)+EVAP_F(:,:,iblock)+  &
                      MELT_F(:,:,iblock)+ROFF_F(:,:,iblock))*salinity_factor   &
                    + SALT_F(:,:,iblock)*sflux_factor)  
      enddo
      !$OMP END PARALLEL DO
 
!-----------------------------------------------------------------------
!
!  balance salt/freshwater in marginal seas
!
!-----------------------------------------------------------------------
 
      if  (lms_balance .and. sfwf_formulation /= 'partially-coupled' ) then
       call ms_balancing (STF(:,:,2,:),EVAP_F, PREC_F, MELT_F,ROFF_F,   &
                          SALT_F, QFLUX, 'salt')
      endif
 
   endif
 
!-----------------------------------------------------------------------
!
!  diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then

      if  (my_task == master_task) then
         call int_to_char (4,iyear   ,cyear  )
         call int_to_char (2,imonth  ,cmonth )
         call int_to_char (2,iday    ,cday   )
         call int_to_char (2,ihour   ,chour  )
         call int_to_char (2,iminute ,cminute)
         call int_to_char (2,isecond ,csecond)
         write(stdout,*)' Global averages of fluxes received from cpl',  &
                         ' at ', cyear,'/',cmonth ,'/',cday,             &
                            ' ', chour,':',cminute,':',csecond
         call shr_sys_flush(stdout)
      endif
 
      m2percm2  = mpercm*mpercm
      do k = 1,nrecv

         n = 0
         !$OMP PARALLEL DO PRIVATE(iblock,n)
         do iblock = 1, nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)

            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               n = n + 1
               WORK1(i,j,iblock) = sbuf(n,k)  ! mult. by TAREA in global_sum_prod
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

!maltrud do we need this update
!        call update_ghost_cells(WORK1, bndy_clinic, &
!                                field_loc_center, field_type_scalar)
         gsum = global_sum_prod(WORK1 , TAREA, distrb_clinic, &
                                 field_loc_center, RCALCT)*m2percm2
         if (my_task == master_task) then
            call cpl_fields_getField(label,k,cpl_fields_c2o_fields)
            write(stdout,1100)'ocn','recv', label ,gsum
         endif
      enddo
      if (my_task == master_task) call shr_sys_flush(stdout)
   endif

   deallocate(sbuf)

1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

!-----------------------------------------------------------------------
!EOC

 end subroutine recv_from_coupler

!***********************************************************************

!BOP
! !IROUTINE: send_to_coupler
! !INTERFACE:

 subroutine send_to_coupler

! !DESCRIPTION:
!  This routine packs fields into a message buffer and sends the
!  message to the flux coupler
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

   character (char_len)    :: label
 
   integer (int_kind) ::  &
      i,j,k,n,iblock,     &
      nsend

   real (r8), dimension(nx_block,ny_block) ::   &
      WORK1, WORK2,      &! local work space
      WORK3, WORK4

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
        WORKA               ! local work space with full block dimension

   real (r8) ::   &
      m2percm2,   &
      gsum

   type (block) :: this_block ! local block info

!-----------------------------------------------------------------------
!
!  initialize control buffer
!
!-----------------------------------------------------------------------

      nsend = cpl_interface_contractNumatt(contractS)
      iblock = 1
      this_block = get_block(blocks_clinic(iblock),iblock)

      allocate(sbuf((this_block%ie-this_block%ib+1)*(this_block%je-this_block%jb+1)  &
      ,       nsend))

      isbuf = 0

      if (check_time_flag(cpl_stop_now)) then
        isbuf(cpl_fields_ibuf_stopnow) = 1
      endif


      isbuf(cpl_fields_ibuf_cdate) = iyear*10000 + imonth*100 + iday
      isbuf(cpl_fields_ibuf_sec) =   &
         ihour*seconds_in_hour + iminute*seconds_in_minute + isecond

      if ( lsend_precip_fact )  &    ! send real as integer
        isbuf(cpl_fields_ibuf_precadj) = precip_fact * 1.0e6_r8  
 
!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points and rotate on T grid
!
!-----------------------------------------------------------------------

      n = 0
   do iblock = 1, nblocks_clinic

      call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2c_So_u),iblock)
      call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2c_So_v),iblock)

      WORK1 = (WORK3*cos(ANGLET(:,:,iblock))+WORK4*sin(-ANGLET(:,:,iblock)))  &
             * mpercm/tlast_coupled
      WORK2 = (WORK4*cos(ANGLET(:,:,iblock))-WORK3*sin(-ANGLET(:,:,iblock)))  &
             * mpercm/tlast_coupled

      this_block = get_block(blocks_clinic(iblock),iblock)

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         sbuf(n,index_o2c_So_u) = WORK1(i,j)
         sbuf(n,index_o2c_So_v) = WORK2(i,j)
      enddo
      enddo

   enddo
!-----------------------------------------------------------------------
!
!     convert and pack surface temperature
!
!-----------------------------------------------------------------------

      n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         sbuf(n,index_o2c_So_t) =   &
             SBUFF_SUM(i,j,iblock,index_o2c_So_t)/tlast_coupled + T0_Kelvin
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     convert and pack salinity
!
!-----------------------------------------------------------------------

      n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         sbuf(n,index_o2c_So_s) =   &
             SBUFF_SUM(i,j,iblock,index_o2c_So_s)*salt_to_ppt/tlast_coupled
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points, then rotate on T grid
!
!-----------------------------------------------------------------------

      n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2c_So_dhdx),iblock)
      call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2c_So_dhdy),iblock)
 
      WORK1 = (WORK3*cos(ANGLET(:,:,iblock)) + WORK4*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled
      WORK2 = (WORK4*cos(ANGLET(:,:,iblock)) - WORK3*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         sbuf(n,index_o2c_So_dhdx) = WORK1(i,j)
         sbuf(n,index_o2c_So_dhdy) = WORK2(i,j)
      enddo
      enddo

   enddo
!-----------------------------------------------------------------------
!
!     pack heat flux due to freezing/melting (W/m^2)
!     QFLUX computation and units conversion occurs in ice.F
!
!-----------------------------------------------------------------------

      n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         sbuf(n,index_o2c_Fioo_q) = QFLUX(i,j,iblock)
      enddo
      enddo
   enddo

   tlast_ice = c0
   AQICE     = c0
   QICE      = c0


!-----------------------------------------------------------------------
!
!  send fields to coupler
!
!-----------------------------------------------------------------------

      call cpl_interface_contractSend(cpl_fields_cplname,contractS,isbuf,sbuf)
 
!-----------------------------------------------------------------------
!
!     diagnostics
!
!-----------------------------------------------------------------------

      if (ldiag_cpl) then
        if (my_task == master_task) then
          call int_to_char (4,iyear   ,cyear  )
          call int_to_char (2,imonth  ,cmonth )
          call int_to_char (2,iday    ,cday   )
          call int_to_char (2,ihour   ,chour  )
          call int_to_char (2,iminute ,cminute)
          call int_to_char (2,isecond ,csecond)
          write(stdout,*) ' Global averages of fluxes sent to cpl at '  &
      ,                   ' ', cyear,'/',cmonth, '/',cday    &
      ,                   ' ', chour,':',cminute,':',csecond
          call shr_sys_flush(stdout)
        endif
 
         m2percm2  = mpercm*mpercm
         do k = 1,nsend
            n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               n = n + 1
               WORKA(i,j,iblock) = sbuf(n,k)
            enddo
            enddo
   enddo
      call update_ghost_cells(WORKA, bndy_clinic, &
                              field_loc_center, field_type_scalar)
      gsum = global_sum_prod(WORKA , TAREA, distrb_clinic, &
                                 field_loc_center, RCALCT)*m2percm2
            if (my_task == master_task) then
               call cpl_fields_getField(label,k,cpl_fields_o2c_fields)
               write(stdout,1100)'ocn','send', label ,gsum
            endif
         enddo
         if (my_task == master_task) call shr_sys_flush(stdout)
      endif

1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

      tlast_coupled = c0

      deallocate(sbuf)

!-----------------------------------------------------------------------
!EOC

 end subroutine send_to_coupler

!***********************************************************************

!BOP
! !IROUTINE: sum_buffer
! !INTERFACE:

 subroutine sum_buffer

! !DESCRIPTION:
!  This routine accumulates sums for averaging fields to
!  be sent to the coupler
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

   real (r8) ::   &
      delt                ! time interval since last step

   integer (int_kind) :: iblock

!-----------------------------------------------------------------------
!
!     zero buffer if this is the first time after a coupling interval
!
!-----------------------------------------------------------------------

      if (tlast_coupled == c0) SBUFF_SUM = c0

!-----------------------------------------------------------------------
!
!     update time since last coupling
!
!-----------------------------------------------------------------------

      if (avg_ts .or. back_to_back) then
        delt = p5*dtt
      else
        delt =    dtt
      endif
      tlast_coupled = tlast_coupled + delt

!-----------------------------------------------------------------------
!
!     accumulate sums of U,V,T,S and GRADP
!     ice formation flux is handled separately in ice routine
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
   SBUFF_SUM(:,:,iblock,index_o2c_So_u) =   &
      SBUFF_SUM(:,:,iblock,index_o2c_So_u) + delt*  &
                                   UVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2c_So_v) =   &
      SBUFF_SUM(:,:,iblock,index_o2c_So_v) + delt*  &
                                   VVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2c_So_t ) =   &
      SBUFF_SUM(:,:,iblock,index_o2c_So_t ) + delt*  &
                                   TRACER(:,:,1,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2c_So_s ) =   &
      SBUFF_SUM(:,:,iblock,index_o2c_So_s ) + delt*  &
                                   TRACER(:,:,1,2,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2c_So_dhdx) =   &
      SBUFF_SUM(:,:,iblock,index_o2c_So_dhdx) + delt*  &
                                   GRADPX(:,:,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2c_So_dhdy) =   &
      SBUFF_SUM(:,:,iblock,index_o2c_So_dhdy) + delt*  &
                                   GRADPY(:,:,curtime,iblock)
   enddo
   !$OMP END PARALLEL DO

 end subroutine sum_buffer


!-----------------------------------------------------------------------
!EOC
 
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

   !$OMP PARALLEL DO PRIVATE(iblock,this_block,WORK)

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

      if (tavg_requested(tavg_EVAP_F)) then
         call accumulate_tavg_field(EVAP_F(:,:,iblock), &
                                    tavg_EVAP_F,iblock,1)
      endif

      if (tavg_requested(tavg_PREC_F)) then
         call accumulate_tavg_field(PREC_F(:,:,iblock), &
                                    tavg_PREC_F,iblock,1)
      endif

      if (tavg_requested(tavg_SNOW_F)) then
         call accumulate_tavg_field(SNOW_F(:,:,iblock), &
                                    tavg_SNOW_F,iblock,1)
      endif

      if (tavg_requested(tavg_MELT_F)) then
         call accumulate_tavg_field(MELT_F(:,:,iblock), &
                                    tavg_MELT_F,iblock,1)
      endif

      if (tavg_requested(tavg_ROFF_F)) then
         call accumulate_tavg_field(ROFF_F(:,:,iblock), &
                                    tavg_ROFF_F,iblock,1)
      endif

      if (tavg_requested(tavg_SALT_F)) then
         call accumulate_tavg_field(SALT_F(:,:,iblock), &
                                    tavg_SALT_F,iblock,1)
      endif

      if (tavg_requested(tavg_SENH_F)) then
         call accumulate_tavg_field(SENH_F(:,:,iblock), &
                                    tavg_SENH_F,iblock,1)
      endif

      if (tavg_requested(tavg_LWUP_F)) then
         call accumulate_tavg_field(LWUP_F(:,:,iblock), &
                                    tavg_LWUP_F,iblock,1)
      endif

      if (tavg_requested(tavg_LWDN_F)) then
         call accumulate_tavg_field(LWDN_F(:,:,iblock), &
                                    tavg_LWDN_F,iblock,1)
      endif

      if (tavg_requested(tavg_MELTH_F)) then
         call accumulate_tavg_field(MELTH_F(:,:,iblock), &
                                    tavg_MELTH_F,iblock,1)
      endif



   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_coupled_forcing
#endif
 
!***********************************************************************

 end module forcing_coupled

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

