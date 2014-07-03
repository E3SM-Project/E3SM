!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing

!BOP
! !MODULE: forcing
!
! !DESCRIPTION:
!  This is the main driver module for all surface and interior
!  forcing.  It contains necessary forcing fields as well as
!  necessary routines for call proper initialization and
!  update routines for those fields.
!
! !REVISION HISTORY:
!  SVN:$Id: forcing.F90 38321 2012-06-29 23:22:21Z mlevy@ucar.edu $
!
! !USES:

   use constants
   use blocks
   use distribution
   use domain
   use grid
   use ice, only: salice, tfreez, FW_FREEZE
   use forcing_ws
   use forcing_shf
   use forcing_sfwf
   use forcing_pt_interior
   use forcing_s_interior
   use forcing_ap
   use forcing_coupled, only: set_combined_forcing, tavg_coupled_forcing,  &
       liceform, qsw_12hr_factor, qsw_distrb_iopt, qsw_distrb_iopt_cosz, &
       tday00_interval_beg, interval_cum_dayfrac, QSW_COSZ_WGHT_NORM, &
       QSW_COSZ_WGHT, compute_cosz
   use forcing_tools
   use passive_tracers, only: set_sflux_passive_tracers
   use prognostic
   use tavg
   use movie, only: define_movie_field, movie_requested, update_movie_field
   use time_management
   use exit_mod
#ifdef CCSMCOUPLED
   use shr_sys_mod, only: shr_sys_abort
#endif

   !*** ccsm
   use sw_absorption, only: set_chl
   use registry
   use forcing_fields

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_forcing,           &
             set_surface_forcing,    &
             tavg_forcing,           &
             movie_forcing

!EOP
!BOC

   integer (int_kind) :: &
      tavg_SHF,          &! tavg_id for surface heat flux
      tavg_SHF_QSW,      &! tavg_id for short-wave solar heat flux
      tavg_SFWF,         &! tavg_id for surface freshwater flux
      tavg_SFWF_WRST,    &! tavg_id for weak restoring freshwater flux
      tavg_TAUX,         &! tavg_id for wind stress in X direction
      tavg_TAUX2,        &! tavg_id for wind stress**2 in X direction
      tavg_TAUY,         &! tavg_id for wind stress in Y direction
      tavg_TAUY2,        &! tavg_id for wind stress**2 in Y direction
      tavg_FW,           &! tavg_id for freshwater flux
      tavg_TFW_T,        &! tavg_id for T flux due to freshwater flux
      tavg_TFW_S,        &! tavg_id for S flux due to freshwater flux
      tavg_U10_SQR        ! tavg_id for U10_SQR 10m wind speed squared from cpl

!-----------------------------------------------------------------------
!
!  movie ids
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      movie_SHF,          &! movie id for surface heat flux
      movie_SFWF,         &! movie id for surface freshwater flux
      movie_TAUX,         &! movie id for wind stress in X direction
      movie_TAUY           ! movie id for wind stress in Y direction

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_forcing
! !INTERFACE:

   subroutine init_forcing

! !DESCRIPTION:
!  Initializes forcing by calling a separate routines for
!  wind stress, heat flux, fresh water flux, passive tracer flux,
!  interior restoring, and atmospheric pressure.
!
! !REVISION HISTORY:
!  same as module

!-----------------------------------------------------------------------
!
!  write out header for forcing options to stdout.
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a15)') 'Forcing options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

!-----------------------------------------------------------------------
!
!  initialize forcing arrays
!
!-----------------------------------------------------------------------

   ATM_PRESS = c0
   FW        = c0
   FW_OLD    = c0
   SMF       = c0
   SMFT      = c0
   STF       = c0
   TFW       = c0

!-----------------------------------------------------------------------
!
!  call individual initialization routines
!
!-----------------------------------------------------------------------

   call init_ws(SMF,SMFT,lsmft_avail)

   !*** NOTE: with bulk NCEP forcing init_shf must be called before
   !***       init_sfwf

   call init_shf (STF)
   call init_sfwf(STF)
   call init_pt_interior
   call init_s_interior
   call init_ap(ATM_PRESS)

!-----------------------------------------------------------------------
!
!  define tavg diagnostic fields
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_SHF, 'SHF', 2,                                &
                          long_name='Total Surface Heat Flux, Including SW', &
                          units='watt/m^2', grid_loc='2110',                 &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SHF_QSW, 'SHF_QSW', 2,                        &
                          long_name='Solar Short-Wave Heat Flux',            &
                          units='watt/m^2', grid_loc='2110',                 &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SFWF,'SFWF',2,                                   &
                          long_name='Virtual Salt Flux in FW Flux formulation', &
                          units='kg/m^2/s', grid_loc='2110',                    &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SFWF_WRST,'SFWF_WRST',2,                         &
                          long_name='Virtual Salt Flux due to weak restoring', &
                          units='kg/m^2/s', grid_loc='2110',                    &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_TAUX,'TAUX',2,                         &
                          long_name='Windstress in grid-x direction', &
                          units='dyne/centimeter^2', grid_loc='2220', &
                          coordinates='ULONG ULAT time')

   call define_tavg_field(tavg_TAUX2,'TAUX2',2,                         &
                          long_name='Windstress**2 in grid-x direction', &
                          units='dyne^2/centimeter^4', grid_loc='2220', &
                          coordinates='ULONG ULAT time')

   call define_tavg_field(tavg_TAUY,'TAUY',2,                         &
                          long_name='Windstress in grid-y direction', &
                          units='dyne/centimeter^2', grid_loc='2220', &
                          coordinates='ULONG ULAT time')

   call define_tavg_field(tavg_TAUY2,'TAUY2',2,                         &
                          long_name='Windstress**2 in grid-y direction', &
                          units='dyne^2/centimeter^4', grid_loc='2220', &
                          coordinates='ULONG ULAT time')

   call define_tavg_field(tavg_FW,'FW',2,                        &
                          long_name='Freshwater Flux',           &
                          units='centimeter/s', grid_loc='2110', &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_TFW_T,'TFW_T',2,                      &
                          long_name='T flux due to freshwater flux', &
                          units='watt/m^2', grid_loc='2110',         &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_TFW_S,'TFW_S',2,                      &
                          long_name='S flux due to freshwater flux (kg of salt/m^2/s)', &
                          units='kg/m^2/s', grid_loc='2110',         &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_U10_SQR,'U10_SQR',2,                  &
                          long_name='10m wind speed squared',      &
                          units='cm^2/^s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!
!  define movie diagnostic fields
!
!-----------------------------------------------------------------------

   call define_movie_field(movie_SHF,'SHF',0,                                &
                          long_name='Total Surface Heat Flux, Including SW', &
                          units='watt/m^2', grid_loc='2110')

   call define_movie_field(movie_SFWF,'SFWF',0,                                 &
                          long_name='Virtual Salt Flux in FW Flux formulation', &
                          units='kg/m^2/s', grid_loc='2110')

   call define_movie_field(movie_TAUX,'TAUX',0,                       &
                          long_name='Windstress in grid-x direction', &
                          units='dyne/centimeter^2', grid_loc='2220')

   call define_movie_field(movie_TAUY,'TAUY',0,                       &
                          long_name='Windstress in grid-y direction', &
                          units='dyne/centimeter^2', grid_loc='2220')

!-----------------------------------------------------------------------
!EOC

 end subroutine init_forcing


!***********************************************************************
!BOP
! !IROUTINE: set_surface_forcing
! !INTERFACE:

 subroutine set_surface_forcing

! !DESCRIPTION:
!  Calls surface forcing routines if necessary.
!  If forcing does not depend on the ocean state, then update
!     forcing if current time is greater than the appropriate
!     interpolation time or if it is the first step.
!  If forcing DOES depend on the ocean state, then call every
!     timestep.  interpolation check will be done within the set\_*
!     routine.
!  Interior restoring is assumed to take place every
!     timestep and is set in subroutine tracer\_update, but
!     updating the data fields must occur here outside
!     any block loops.
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

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      TFRZ               
   integer (int_kind) :: index_qsw, iblock
   real (r8) ::  &
      cosz_day,  &
      qsw_eps


#ifdef _HIRES
   qsw_eps = -1.e-8_r8
#else
   qsw_eps = c0
#endif

!-----------------------------------------------------------------------
!
!  Get any interior restoring data and interpolate if necessary.
!
!-----------------------------------------------------------------------

   call get_pt_interior_data
   call get_s_interior_data

!-----------------------------------------------------------------------
!
!  Call individual forcing update routines.
!
!-----------------------------------------------------------------------

   if (lsmft_avail) then
      call set_ws(SMF,SMFT=SMFT)
   else
      call set_ws(SMF)
   endif

   !*** NOTE: with bulk NCEP and partially-coupled forcing 
   !***       set_shf must be called before set_sfwf

   call set_shf(STF)
   call set_sfwf(STF,FW,TFW)

   if ( shf_formulation  == 'partially-coupled' .or.  &
        sfwf_formulation == 'partially-coupled' ) then
      call set_combined_forcing(STF,FW,TFW)
   endif


!-----------------------------------------------------------------------
!
! apply qsw 12hr if chosen
!
!-----------------------------------------------------------------------

      index_qsw = mod(nsteps_this_interval,nsteps_per_interval) + 1

      if ( qsw_distrb_iopt == qsw_distrb_iopt_cosz ) then
         cosz_day = tday00_interval_beg + interval_cum_dayfrac(index_qsw-1) &
            - interval_cum_dayfrac(nsteps_per_interval)

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock = 1, nblocks_clinic

            call compute_cosz(cosz_day, iblock, QSW_COSZ_WGHT(:,:,iblock))

            where (QSW_COSZ_WGHT_NORM(:,:,iblock) > c0)
               QSW_COSZ_WGHT(:,:,iblock) = QSW_COSZ_WGHT(:,:,iblock) &
                  * QSW_COSZ_WGHT_NORM(:,:,iblock)
            elsewhere
               QSW_COSZ_WGHT(:,:,iblock) = c1
            endwhere

            SHF_QSW(:,:,iblock) = QSW_COSZ_WGHT(:,:,iblock) &
               * SHF_COMP(:,:,iblock,shf_comp_qsw)

         enddo
         !$OMP END PARALLEL DO

      else

         if (registry_match('lcoupled')) then
            SHF_QSW = qsw_12hr_factor(index_qsw)*SHF_COMP(:,:,:,shf_comp_qsw)
         endif

      endif

      if ( registry_match('lcoupled') &
           .and. sfwf_formulation /= 'partially-coupled'  &
           .and. sfc_layer_type == sfc_layer_varthick .and. &
           .not. lfw_as_salt_flx .and. liceform ) then
        FW  = SFWF_COMP(:,:,:,  sfwf_comp_cpl)
        TFW =  TFW_COMP(:,:,:,:, tfw_comp_cpl)
      endif

      if ( sfc_layer_type == sfc_layer_varthick .and.   &
           .not. lfw_as_salt_flx .and. liceform ) then
        FW = FW + FW_FREEZE

        call tfreez(TFRZ,TRACER(:,:,1,2,curtime,:))

        TFW(:,:,1,:) = TFW(:,:,1,:) + FW_FREEZE(:,:,:)*TFRZ(:,:,:)
        TFW(:,:,2,:) = TFW(:,:,2,:) + FW_FREEZE(:,:,:)*salice
      endif


   call set_ap(ATM_PRESS)

   if (nt > 2)  &
      call set_sflux_passive_tracers(U10_SQR,IFRAC,ATM_PRESS,STF)

   call set_chl   

#ifdef CCSMCOUPLED
   if (ANY(SHF_QSW < qsw_eps)) then
      call shr_sys_abort('(set_surface_forcing) ERROR: SHF_QSW < qsw_eps in set_surface_forcing')
   endif
#endif


!-----------------------------------------------------------------------
!EOC

 end subroutine set_surface_forcing

!***********************************************************************
!BOP
! !IROUTINE: tavg_forcing
! !INTERFACE:

 subroutine tavg_forcing

! !DESCRIPTION:
!  This routine accumulates tavg diagnostics related to surface
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

      if (accumulate_tavg_now(tavg_SHF)) then
         where (KMT(:,:,iblock) > 0)
            WORK = (STF(:,:,1,iblock)+SHF_QSW(:,:,iblock))/ &
                   hflux_factor ! W/m^2
         elsewhere
            WORK = c0
         end where

         call accumulate_tavg_field(WORK,tavg_SHF,iblock,1)
      endif

      if (accumulate_tavg_now(tavg_SHF_QSW)) then
         where (KMT(:,:,iblock) > 0)
            WORK = SHF_QSW(:,:,iblock)/hflux_factor ! W/m^2
         elsewhere
            WORK = c0
         end where

         call accumulate_tavg_field(WORK,tavg_SHF_QSW,iblock,1)
      endif

      if (accumulate_tavg_now(tavg_SFWF)) then
         if (sfc_layer_type == sfc_layer_varthick .and. &
             .not. lfw_as_salt_flx) then
            where (KMT(:,:,iblock) > 0)
               WORK = FW(:,:,iblock)*seconds_in_year*mpercm ! m/yr
            elsewhere
               WORK = c0
            end where
         else
            where (KMT(:,:,iblock) > 0) ! convert to kg(freshwater)/m^2/s
               WORK = STF(:,:,2,iblock)/salinity_factor
            elsewhere
               WORK = c0
            end where
         endif

         call accumulate_tavg_field(WORK,tavg_SFWF,iblock,1)
      endif

      if (accumulate_tavg_now(tavg_SFWF_WRST)) then
         if ( sfwf_formulation == 'partially-coupled' ) then
            where (KMT(:,:,iblock) > 0) ! convert to kg(freshwater)/m^2/s
               WORK = SFWF_COMP(:,:,iblock,sfwf_comp_wrest)/salinity_factor
            elsewhere
               WORK = c0
            end where
         else
            WORK = c0
         endif
         call accumulate_tavg_field(WORK,tavg_SFWF_WRST,iblock,1)
      endif

      call accumulate_tavg_field(SMF(:,:,1,iblock), tavg_TAUX,iblock,1)
      call accumulate_tavg_field(SMF(:,:,1,iblock)**2, tavg_TAUX2,iblock,1)
      call accumulate_tavg_field(SMF(:,:,2,iblock), tavg_TAUY,iblock,1)
      call accumulate_tavg_field(SMF(:,:,2,iblock)**2, tavg_TAUY2,iblock,1)
      call accumulate_tavg_field(FW (:,:,iblock), tavg_FW,iblock,1)
      call accumulate_tavg_field(TFW(:,:,1,iblock)/hflux_factor, tavg_TFW_T,iblock,1)
      call accumulate_tavg_field(TFW(:,:,2,iblock)*rho_sw*c10, tavg_TFW_S,iblock,1)
      call accumulate_tavg_field(U10_SQR(:,:,iblock), tavg_U10_SQR,iblock,1)


   end do

   !$OMP END PARALLEL DO

   if (registry_match('lcoupled')) call tavg_coupled_forcing

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_forcing


!***********************************************************************
!BOP
! !IROUTINE: movie_forcing
! !INTERFACE:

 subroutine movie_forcing

! !DESCRIPTION:
!  This routine accumulates movie diagnostics related to surface
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
      WORK                ! local temp space for movie diagnostics

!-----------------------------------------------------------------------
!
!  compute and dump movie forcing diagnostics
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block,WORK)

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

!-----------------------------------------------------------------------
!
!     dump movie diagnostics if requested
!
!-----------------------------------------------------------------------

      if (movie_requested(movie_SHF) ) then
         where (KMT(:,:,iblock) > 0)
            WORK = (STF(:,:,1,iblock)+SHF_QSW(:,:,iblock))/ &
                   hflux_factor ! W/m^2
         elsewhere
            WORK = c0
         end where
         call update_movie_field(WORK, movie_SHF, iblock, 1)
      endif

      if (movie_requested(movie_SFWF) ) then
         if (sfc_layer_type == sfc_layer_varthick .and. &
             .not. lfw_as_salt_flx) then
            where (KMT(:,:,iblock) > 0)
               WORK = FW(:,:,iblock)*seconds_in_year*mpercm ! m/yr
            elsewhere
               WORK = c0
            end where
         else
            where (KMT(:,:,iblock) > 0) ! convert to kg(freshwater)/m^2/s
               WORK = STF(:,:,2,iblock)/salinity_factor
            elsewhere
               WORK = c0
            end where
         endif
         call update_movie_field(WORK, movie_SFWF, iblock, 1)
      endif

      if (movie_requested(movie_TAUX) ) then
         call update_movie_field(SMF(:,:,1,iblock), &
                                    movie_TAUX,iblock,1)
      endif

      if (movie_requested(movie_TAUY) ) then
         call update_movie_field(SMF(:,:,2,iblock), &
                                    movie_TAUY,iblock,1)
      endif


   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine movie_forcing


!***********************************************************************

 end module forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
