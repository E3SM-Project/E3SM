!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_CouplingMod

!BOP
! !MODULE: POP_CouplingMod
! !DESCRIPTION:
!  This module contains 
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod

#ifdef CCSMCOUPLED
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use constants
   use blocks
   use domain
   use exit_mod
   use forcing_shf, only: SHF_QSW
   use forcing_sfwf, only: lsend_precip_fact, precip_fact
   use forcing_fields
   use forcing_coupled, only : ncouple_per_day, timer_recv_from_cpl,  &
       update_ghost_cells_coupler_fluxes, rotate_wind_stress
   use ice, only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, tlast_ice
   use grid, only: TLAT, TLON, REGION_MASK
   use global_reductions, only: global_sum_prod
   use io_tools, only: document
   use named_field_mod, only: named_field_register, named_field_get_index, &
       named_field_set, named_field_get
   use prognostic
   use timers, only: get_timer
   use time_management

   use cpl_contract_mod
   use cpl_interface_mod
   use cpl_fields_mod
   use forcing_coupled
#endif

   implicit none
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: pop_coupling,                   &
             pop_send_to_coupler,            &
             pop_init_coupler_comm,          &
             pop_recv_from_coupler,          &
             pop_unpack_fluxes_from_coupler, &
             pop_sum_buffer,                 &
             pop_prepare_send_to_coupler,    &
             pop_allocate_sbuf_send,         &
             pop_allocate_sbuf_recv,         &
             pop_deallocate_sbuf

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
   real (POP_r8),   &
      dimension(:,:,:,:), allocatable ::  &
      SBUFF_SUM           ! accumulated sum of send buffer quantities
                          ! for averaging before being sent

   integer (POP_i4), private ::   &
      lsize_total,         &! aggregate total size of physical domain 
                            !   over all nblocks_clinic
      lisize_total,        &! aggregate total i points in physical domain
                            !   over all nblocks_clinic
      ljsize_total          ! aggregate total j points in physical domain
                            !   over all nblocks_clinic

   integer (POP_i4), dimension(cpl_fields_ibuf_total) ::  &
      isbuf,               &! integer control buffer for sends
      irbuf                 ! integer control buffer for receives
 
   type(cpl_contract) ::  &
      contractS,          &! contract for sends to coupler
      contractR            ! contract for receives from coupler
 
   real (POP_r8), dimension(:,:), allocatable ::  &
      sbuf                 ! temporary send/recv buffer

   real (POP_r8) ::  &
      tlast_coupled

   integer (POP_i4), private ::   &
      nsend, nrecv

!-----------------------------------------------------------------------
!  The following variables are used in the exchange of 
!  information between cpl6 and the ocean code.
!
!  ocn --> cpl6
!  ============
!    cpl_fields_ibuf_total -- length of integer ocean "send buffer" vector (isbuf)
!    nsend --  total number of 2D fields sent to cpl6 from ocn
!
!    integer send buffer indices (isbuf in pop_init_coupler_comm):  
!
!     o  cpl_fields_ibuf_cdate   -- oceans character date string (yyyymmdd)
!     o  cpl_fields_ibuf_sec     -- oceans character time string (seconds)
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
!    real send buffer indices (sbuf in pop_init_coupler_comm):
!
!     o  cpl_fields_grid_lon   -- radian*TLON(i,j)
!     o  cpl_fields_grid_lat   -- radian*TLAT(i,j)
!     o  cpl_fields_grid_area  -- TAREA(i,j)/(radius*radius)
!     o  cpl_fields_grid_mask  -- float(REGION_MASK(i,j))
!     o  cpl_fields_grid_index -- (j_global(j)-1)*(nx_global)+i_global(i)
!     o  cpl_fields_grid_frac  -- "cell fraction" information
!
!    real send buffer indices 
!
!      o  index_o2c_So_u     -- surface u velocity
!      o  index_o2c_So_v     -- surface v velocity
!      o  index_o2c_So_t     -- surface temperature
!      o  index_o2c_So_s     -- surface salinity
!      o  index_o2c_So_dhdx  -- e,w surface slope
!      o  index_o2c_So_dhdy  -- n,s surface slope
!      o  index_o2c_Fioo_q   -- qflux
!      o  index_o2c_Faoo_fco2-- co2 flux
!
!
!    cpl6 --> ocn  
!    ============
!
!    cpl_fields_ibuf_total -- length of integer ocean "receive buffer" vector (irbuf)
!
!    integer receive buffer indices (irbuf in pop_recv_from_coupler):
!
!     o  cpl_fields_ibuf_stopnow  -- stop ocean integration now
!     o  cpl_fields_ibuf_infobug  -- write ocean/coupler diagnostics now  
!     o  cpl_fields_ibuf_resteod  -- write ocean restart files at end of day
!     o  cpl_fields_ibuf_histeod  -- write ocean history files at end of day
!     o  cpl_fields_ibuf_histtavg -- write ocean "tavg"  files at end of day
!     o  cpl_fields_ibuf_diageod  -- write ocean diagnostics   at end of day
!
!    real receive buffer indices (sbuf in pop_recv_from_coupler):
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
!     o  index_c2o_Foxx_rofl   -- river runoff flux
!     o  index_c2o_Sa_co2prog  -- bottom atm level prognostic co2
!
!-----------------------------------------------------------------------

   integer(kind=POP_i4) :: index_o2c_So_t        ! temperature
   integer(kind=POP_i4) :: index_o2c_So_u        ! velocity, zonal
   integer(kind=POP_i4) :: index_o2c_So_v        ! velocity, meridional
   integer(kind=POP_i4) :: index_o2c_So_s        ! salinity
   integer(kind=POP_i4) :: index_o2c_So_dhdx     ! surface slope, zonal
   integer(kind=POP_i4) :: index_o2c_So_dhdy     ! surface slope, meridional
   integer(kind=POP_i4) :: index_o2c_Fioo_q      ! heat of fusion (q>0) melt pot (q<0)
   integer(kind=POP_i4) :: index_o2c_Faoo_fco2   ! co2 flux

   integer(kind=POP_i4) :: index_c2o_Si_ifrac    ! state: ice fraction
   integer(kind=POP_i4) :: index_c2o_Sa_pslv     ! state: sea level pressure
   integer(kind=POP_i4) :: index_c2o_Faoc_duu10n ! state: 10m wind speed squared
   integer(kind=POP_i4) :: index_c2o_Foxx_taux   ! wind stress: zonal
   integer(kind=POP_i4) :: index_c2o_Foxx_tauy   ! wind stress: meridional
   integer(kind=POP_i4) :: index_c2o_Foxx_swnet  ! heat flux: shortwave net
   integer(kind=POP_i4) :: index_c2o_Foxx_lat    ! heat flux: latent
   integer(kind=POP_i4) :: index_c2o_Foxx_sen    ! heat flux: sensible
   integer(kind=POP_i4) :: index_c2o_Foxx_lwup   ! heat flux: long-wave up
   integer(kind=POP_i4) :: index_c2o_Foxx_lwdn   ! heat flux: long-wave dow
   integer(kind=POP_i4) :: index_c2o_Foxx_melth  ! heat flux: melt
   integer(kind=POP_i4) :: index_c2o_Foxx_salt   ! salt flux
   integer(kind=POP_i4) :: index_c2o_Foxx_prec   ! water flux: rain+snow
   integer(kind=POP_i4) :: index_c2o_Foxx_snow   ! water flux: snow
   integer(kind=POP_i4) :: index_c2o_Foxx_rain   ! water flux: rain
   integer(kind=POP_i4) :: index_c2o_Foxx_evap   ! water flux: evap
   integer(kind=POP_i4) :: index_c2o_Foxx_meltw  ! water flux: melt
   integer(kind=POP_i4) :: index_c2o_Foxx_rofl   ! water flux: runoff
   integer(kind=POP_i4) :: index_c2o_Sa_co2prog  ! bottom atm level prognostic co2

   logical (POP_logical) ::   &
      ldiag_cpl = .false.

   integer (POP_i4), private ::   &
      cpl_stop_now,        &! flag id for stop_now flag
      cpl_write_restart,   &! flag id for write restart
      cpl_write_history,   &! flag id for write history
      cpl_write_tavg,      &! flag id for write tavg      
      cpl_diag_global,     &! flag id for computing diagnostics
      cpl_diag_transp       ! flag id for computing diagnostics
#endif


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: pop_coupling
! !INTERFACE:

 subroutine pop_coupling(lcoupled_ts, errorCode)

! !DESCRIPTION:
!  This routine call coupler communication routines to set
!  surface forcing data
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical(POP_logical), intent(in) :: &
      lcoupled_ts      ! flag indicating coupled timestep status

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode        ! returned error code

!EOP
!BOC

#if CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_logical), save ::    &
      send = .false.               ! flag for controlling pop_send_to_coupler
     
#endif
!-----------------------------------------------------------------------
!
!  call necessary send and recv routines
!-----------------------------------------------------------------------

   errorCode = POP_Success

#if CCSMCOUPLED
   if (send) then
      call pop_send_to_coupler(lcoupled_ts, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'pop_coupling: error in send_to_coupler')
         return
      endif
   endif

   if (lcoupled_ts .or. nsteps_run == 0) then
      call pop_recv_from_coupler(errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'pop_coupling: error in recv_from_coupler')
         return
      endif

      call pop_set_coupled_forcing
   endif

   send = .true.

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_coupling


!***********************************************************************

!BOP
! !IROUTINE: pop_send_to_coupler
! !INTERFACE:

 subroutine pop_send_to_coupler(lcoupled_ts, errorCode)


! !DESCRIPTION:
!  This routine calls the routines necessary to send pop fields to
!  the CCSM flux coupler
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical(POP_logical), intent(in) :: &
      lcoupled_ts      ! flag indicating coupled timestep status

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode        ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     send state variables to coupler
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

#if CCSMCOUPLED
   call pop_sum_buffer

   if (lcoupled_ts) then
     call timer_stop  (timer_recv_to_send)
   endif

   if (lcoupled_ts .or. nsteps_run == 0 ) then

     call timer_start (timer_send_to_cpl)
     call pop_allocate_sbuf_send

     call pop_prepare_send_to_coupler(errorCode)
     if (errorCode /= POP_Success) then
        call POP_ErrorSet(errorCode, &
           'pop_send_to_coupler: error in prepare_send')
        return
     endif

     call cpl_interface_contractSend(cpl_fields_cplname,contractS,isbuf,sbuf)
     call pop_deallocate_sbuf
     call timer_stop  (timer_send_to_cpl)

     call timer_start (timer_send_to_recv)
     call timer_stop  (timer_send_to_recv) 
 
     call timer_start (timer_recv_to_send)
   endif

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_send_to_coupler

!***********************************************************************

!BOP
! !IROUTINE: pop_init_coupler_comm 
! !INTERFACE:

 subroutine pop_init_coupler_comm 

! !DESCRIPTION:
! 
!  This routine initializes everything necessary for pop to couple
!  concurrently with the CCSM3 flux coupler, version 6 (cpl6). 
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::   &
      iblock, n,           &
      i,j,k,               &
      ib,ie,jb,je

   type (block) ::       &
      this_block          ! block information for current block


!-----------------------------------------------------------------------
!
!   ensure that this routine is called after pop_init_coupled is called;
!
!-----------------------------------------------------------------------

   if (.not. registry_match('pop_init_coupled')) then
     call exit_POP(sigAbort, &
        'ERROR: pop_init_coupler_comm must be called after pop_init_coupled')
   endif

   tlast_coupled = c0

!-----------------------------------------------------------------------
!
!   initialize the send buffer (but do not send to coupler)
!
!-----------------------------------------------------------------------

   isbuf = 0

   isbuf(cpl_fields_ibuf_cdate  ) = iyear*10000 + imonth*100 + iday
   isbuf(cpl_fields_ibuf_sec    ) =   &
   ihour*seconds_in_hour + iminute*seconds_in_minute + isecond

   !***  determine total size of buffers over all blocks
   lsize_total  = 0
   lisize_total = 0
   ljsize_total = 0

   do iblock = 1, nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)
     do i=this_block%ib,this_block%ie
       do j=this_block%jb,this_block%je
         lsize_total = lsize_total + 1
       enddo ! j
     enddo ! i
   enddo

   isbuf(cpl_fields_ibuf_lsize  ) = lsize_total
   isbuf(cpl_fields_ibuf_lisize ) = lisize_total
   isbuf(cpl_fields_ibuf_ljsize ) = ljsize_total
   isbuf(cpl_fields_ibuf_gsize  ) = nx_global*ny_global
   isbuf(cpl_fields_ibuf_gisize ) = nx_global
   isbuf(cpl_fields_ibuf_gjsize ) = ny_global
   isbuf(cpl_fields_ibuf_ncpl   ) = ncouple_per_day
   isbuf(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
   isbuf(cpl_fields_ibuf_dead   ) = 0           ! not a dead model

   allocate(sbuf(lsize_total,cpl_fields_grid_total))
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
         sbuf(n,cpl_fields_grid_frac ) = 1.0_POP_r8
      enddo
      enddo
   enddo


!-----------------------------------------------------------------------
!
!  initialize the contracts
!
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
   index_o2c_Faoo_fco2= cpl_interface_contractIndex(contractS,'Faoo_fco2',perrWith='quiet')

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
   index_c2o_Foxx_rofl   = cpl_interface_contractIndex(contractR,'Foxx_rofl')
   index_c2o_Sa_co2prog  = cpl_interface_contractIndex(contractR,'Sa_co2prog',  &
                           perrWith='quiet')

!-----------------------------------------------------------------------
!
!  register non-standard incoming fields
!
!-----------------------------------------------------------------------

   if (index_c2o_Sa_co2prog > 0) then
     call named_field_register('ATM_CO2', ATM_CO2_nf_ind)
   endif

!-----------------------------------------------------------------------
!
!  Initialize flags and shortwave absorption profile
!  Note that the cpl_write_xxx flags have _no_ default value;
!  therefore, they must be explicitly set .true. and .false.
!  at the appropriate times
!
!-----------------------------------------------------------------------

   call init_time_flag('stop_now',         cpl_stop_now, default=.false.)
   call init_time_flag('cpl_write_restart',cpl_write_restart)
   call init_time_flag('cpl_write_history',cpl_write_history)
   call init_time_flag('cpl_write_tavg'   ,cpl_write_tavg)
   call init_time_flag('cpl_diag_global'  ,cpl_diag_global)
   call init_time_flag('cpl_diag_transp'  ,cpl_diag_transp)


!-----------------------------------------------------------------------
!
!  register subroutine call
!
!-----------------------------------------------------------------------

   call register_string('pop_init_coupler_comm')

   call flushm (stdout)

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_init_coupler_comm

 
!***********************************************************************

!BOP
! !IROUTINE: pop_recv_from_coupler
! !INTERFACE:

 subroutine pop_recv_from_coupler(errorCode)

! !DESCRIPTION:
!  This routine receives surface flux data from coupler
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

#ifdef CCSMCOUPLED

      call timer_start (timer_recv_from_cpl)

      call pop_allocate_sbuf_recv

      call cpl_interface_contractRecv(cpl_fields_cplname,contractR, &
         irbuf,sbuf)

      call pop_unpack_fluxes_from_coupler(errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'pop_recv_from_coupler: error in unpack_fluxes')
         return
      endif

      call pop_deallocate_sbuf

      call timer_stop  (timer_recv_from_cpl)


#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_recv_from_coupler
 

!***********************************************************************

!BOP
! !IROUTINE: pop_unpack_fluxes_from_coupler
! !INTERFACE:

 subroutine pop_unpack_fluxes_from_coupler(errorCode)

! !DESCRIPTION:
!  This routine receives message from coupler with surface flux data
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len)    :: label, message
 
   integer (POP_i4) ::  &
      i,j,k,n,iblock

   real (POP_r8), dimension(nx_block,ny_block) ::  &
      WORKB

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
      WORK1, WORK2        ! local work space

   real (POP_r8) ::  &
      m2percm2,  &
      gsum

   type (block) :: this_block ! local block info
#endif

!-----------------------------------------------------------------------
!
!  zero out padded cells 
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

#ifdef CCSMCOUPLED
   WORK1 = c0
   WORK2 = c0

!-----------------------------------------------------------------------
!
!  set character date and time information (cyear, cmonth, etc)
!
!-----------------------------------------------------------------------
  
   call ccsm_char_date_and_time

!-----------------------------------------------------------------------
!
!  check coupler flags and respond appropriately
!
!-----------------------------------------------------------------------

   if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
     call override_time_flag(cpl_stop_now,value=.true.)
     write(message,'(6a,1x,5a)')  'cpl requests termination now: ', &
        cyear,'/',cmonth,'/',cday,   chour,':',cminute,':',csecond
     call document ('pop_recv_from_coupler', message)
     RETURN
   endif

   if (irbuf(cpl_fields_ibuf_infobug) >= 2) then
     ldiag_cpl = .true. 
   else
     ldiag_cpl = .false.
   endif

   if (irbuf(cpl_fields_ibuf_resteod) == 1) then
     call override_time_flag(cpl_write_restart,value=.true.)
     write(message,'(6a)') 'cpl requests restart file at eod  ',  &
                            cyear,'/',cmonth,'/',cday
     call document ('pop_recv_from_coupler', message)
   endif
  
   if (irbuf(cpl_fields_ibuf_diageod) == 1) then
     call override_time_flag(cpl_diag_global,value=.true.)
     call override_time_flag(cpl_diag_transp,value=.true.)
     write(message,'(6a)') ' cpl requests diagnostics at eod ',  &
                             cyear,'/',cmonth,'/',cday
     call document ('pop_recv_from_coupler', message)
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
   !*** do halo updates now to ensure correct T->U grid
   !***

   call POP_HaloUpdate(WORK1, POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'pop_unpack_fluxes: error updating taux halo')
      return
   endif

   call POP_HaloUpdate(WORK2, POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'pop_unpack_fluxes: error updating tauy halo')
      return
   endif

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

!-----------------------------------------------------------------------
!
!  rotate true zonal/meridional wind stress into local coordinates,
!  convert to dyne/cm**2, and shift SMFT to U grid
!
!-----------------------------------------------------------------------

      call rotate_wind_stress(WORK1, WORK2,iblock)

!-----------------------------------------------------------------------
!
!  unpack and distribute fresh water flux and salt flux
!
!  NOTE: if there are code changes associated with changing the names or
!        the number of fluxes received from the coupler, then subroutine
!        update_ghost_cells_coupler_fluxes will need to be modified also
!
!-----------------------------------------------------------------------


      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         SNOW_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_snow)
         WORKB (i,j       ) = sbuf(n,index_c2o_Foxx_rain)
         EVAP_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_evap)
         MELT_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_meltw)
         ROFF_F(i,j,iblock) = sbuf(n,index_c2o_Foxx_rofl)
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

!-----------------------------------------------------------------------
!
!  update ghost cells for fluxes received from the coupler
!
!-----------------------------------------------------------------------

   call update_ghost_cells_coupler_fluxes(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
        'pop_unpack_fluxes: error in update_ghost_cells_coupler_fluxes')
      return
   endif

!-----------------------------------------------------------------------
!
!  unpack atmospheric CO2
!
!-----------------------------------------------------------------------

   if (index_c2o_Sa_co2prog > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = sbuf(n,index_c2o_Sa_co2prog)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1, POP_haloClinic,          &
                          POP_gridHorzLocCenter,          &
                          POP_fieldKindScalar, errorCode, &
                          fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'pop_unpack_fluxes: error updating co2 halo')
         return
      endif

      call named_field_set(ATM_CO2_nf_ind, WORK1)
   endif
 
!-----------------------------------------------------------------------
!
!  diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then

     write(message,'(6a,1x,5a)')  &
         ' Global averages of fluxes received from cpl at ',  &
           cyear,'/',cmonth ,'/',cday,  chour,':',cminute,':',csecond
     call document ('pop_recv_from_coupler', message)
 
     m2percm2  = mpercm*mpercm
     nrecv = cpl_interface_contractNumatt(contractR)
     do k = 1,nrecv

         n = 0
         do iblock = 1, nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)

            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               n = n + 1
               WORK1(i,j,iblock) = sbuf(n,k)  ! mult. by TAREA in global_sum_prod
            enddo
            enddo
         enddo

         gsum = global_sum_prod(WORK1 , TAREA, distrb_clinic, &
                                 field_loc_center, RCALCT)*m2percm2
         if (my_task == master_task) then
            call cpl_fields_getField(label,k,cpl_fields_c2o_fields)
            write(stdout,1100)'ocn','recv', label ,gsum
            call shr_sys_flush(stdout)
         endif
      enddo
   endif


1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_unpack_fluxes_from_coupler
 
!***********************************************************************

!BOP
! !IROUTINE: pop_prepare_send_to_coupler
! !INTERFACE:

 subroutine pop_prepare_send_to_coupler(errorCode)

! !DESCRIPTION:
!  This routine packs fields into a message buffer in preparation
!  for sending to the coupler
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len)    :: label, message
 
   integer (POP_i4) ::  &
      i,j,k,n,iblock

   real (POP_r8), dimension(nx_block,ny_block) ::   &
      WORK1, WORK2,      &! local work space
      WORK3, WORK4

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
        WORKA               ! local work space with full block dimension

   real (POP_r8) ::   &
      m2percm2,   &
      gsum

   type (block) :: this_block ! local block info

#endif
!-----------------------------------------------------------------------
!
!  initialize control buffer
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

#ifdef CCSMCOUPLED
   if (.not. allocated(sbuf)) call exit_POP(sigAbort, &
       'ERROR: sbuf not allocated in subroutine pop_prepare_send_to_coupler')

   isbuf = 0

   if (check_time_flag(cpl_stop_now)) then
     isbuf(cpl_fields_ibuf_stopnow) = 1
   endif


   isbuf(cpl_fields_ibuf_cdate) = iyear*10000 + imonth*100 + iday
   isbuf(cpl_fields_ibuf_sec) =   &
      ihour*seconds_in_hour + iminute*seconds_in_minute + isecond

   if ( lsend_precip_fact )  &    ! send real as integer
     isbuf(cpl_fields_ibuf_precadj) = precip_fact * 1.0e6_POP_r8  
 
!-----------------------------------------------------------------------
!
!  interpolate onto T-grid points and rotate on T grid
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)

     call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2c_So_u),iblock)
     call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2c_So_v),iblock)

     WORK1 = (WORK3*cos(ANGLET(:,:,iblock))+WORK4*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled
     WORK2 = (WORK4*cos(ANGLET(:,:,iblock))-WORK3*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled

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
!  convert and pack surface temperature
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
!  convert and pack salinity
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
!  interpolate onto T-grid points, then rotate on T grid
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
!  pack heat flux due to freezing/melting (W/m^2)
!  QFLUX computation and units conversion occurs in ice.F
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
!  pack co2 flux, if requested (kg CO2/m^2/s)
!  units conversion occurs where co2 flux is computed
!
!-----------------------------------------------------------------------

   if (index_o2c_Faoo_fco2 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            sbuf(n,index_o2c_Faoo_fco2) = &
               SBUFF_SUM(i,j,iblock,index_o2c_Faoo_fco2)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!
!  diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then
      call ccsm_char_date_and_time
      write(message,'(6a,1x,5a)')' Global averages of fluxes sent to cpl at ', &
           cyear,'/',cmonth, '/',cday,  chour,':',cminute,':',csecond
      call document ('pop_recv_from_coupler', message)

      m2percm2  = mpercm*mpercm
      nsend = cpl_interface_contractNumatt(contractS)
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

        call POP_HaloUpdate(WORKA, POP_haloClinic, &
                            POP_gridHorzLocCenter, &
                            POP_fieldKindScalar,   &
                            errorCode, fillValue=0.0_POP_r8)

        if (errorCode /= POP_Success) then
           call POP_ErrorSet(errorCode, &
              'pop_prepare_send: error updating halo for state')
           return
        endif

        gsum = global_sum_prod(WORKA , TAREA, distrb_clinic, &
                                   field_loc_center, RCALCT)*m2percm2
        if (my_task == master_task) then
           call cpl_fields_getField(label,k,cpl_fields_o2c_fields)
           write(stdout,1100)'ocn','send', label ,gsum
         endif
      enddo ! k
      if (my_task == master_task) call shr_sys_flush(stdout)
   endif

1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

   tlast_coupled = c0

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_prepare_send_to_coupler

!***********************************************************************

!BOP
! !IROUTINE: pop_allocate_sbuf_send
! !INTERFACE:

 subroutine pop_allocate_sbuf_send

! !DESCRIPTION:
!  This routine allocates sbuf prior to sending information to the coupler
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

 
   integer (POP_i4) ::  &
      i,j,k,n,iblock

   type (block) :: this_block ! local block info


   nsend = cpl_interface_contractNumatt(contractS)
   allocate(sbuf(lsize_total, nsend))

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_allocate_sbuf_send

!***********************************************************************

!BOP
! !IROUTINE: pop_allocate_sbuf_recv
! !INTERFACE:

 subroutine pop_allocate_sbuf_recv

! !DESCRIPTION:
!  This routine allocates sbuf prior to receiving information to the coupler
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::  &
      iblock

   type (block) :: this_block ! local block info


   nrecv = cpl_interface_contractNumatt(contractR)
   allocate(sbuf((lsize_total), nrecv))


#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_allocate_sbuf_recv


!***********************************************************************

!BOP
! !IROUTINE: pop_deallocate_sbuf
! !INTERFACE:

 subroutine pop_deallocate_sbuf

! !DESCRIPTION:
!  This routine deallocates sbuf
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
#ifdef CCSMCOUPLED

      deallocate(sbuf)

#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine pop_deallocate_sbuf

!***********************************************************************

!BOP
! !IROUTINE: pop_sum_buffer
! !INTERFACE:

 subroutine pop_sum_buffer

! !DESCRIPTION:
!  This routine accumulates sums for averaging fields to
!  be sent to the coupler
!
! !REVISION HISTORY:
!  same as module
!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK                ! local work arrays

   real (POP_r8) ::   &
      delt,             & ! time interval since last step
      delt_last           ! time interval for previous step

   integer (POP_i4) :: &
      iblock,           & ! block index
      sflux_co2_nf_ind = 0! named field index of fco2

   logical (POP_logical) :: &
      first = .true.      ! only true for first call

!-----------------------------------------------------------------------
!
!  zero buffer if this is the first time after a coupling interval
!
!-----------------------------------------------------------------------

   if (tlast_coupled == c0) SBUFF_SUM = c0

!-----------------------------------------------------------------------
!
!  update time since last coupling
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
!  allow for fco2 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep
!
!-----------------------------------------------------------------------

   if (index_o2c_Faoo_fco2 > 0) then
      if (sflux_co2_nf_ind == 0) then
         call named_field_get_index('SFLUX_CO2', sflux_co2_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!-----------------------------------------------------------------------
!
!  accumulate sums of U,V,T,S and GRADP
!  accumulate sum of co2 flux, if requested
!     implicitly use zero flux if fco2 field not registered yet
!  ice formation flux is handled separately in ice routine
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

   if (index_o2c_Faoo_fco2 > 0 .and. sflux_co2_nf_ind > 0) then
      call named_field_get(sflux_co2_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2c_Faoo_fco2) = &
         SBUFF_SUM(:,:,iblock,index_o2c_Faoo_fco2) + delt_last*WORK(:,:,iblock)
   endif

   enddo
   !$OMP END PARALLEL DO

   first = .false.

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_sum_buffer
 
!***********************************************************************

 end module POP_CouplingMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
