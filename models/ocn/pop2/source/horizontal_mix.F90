!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module horizontal_mix

!BOP
! !MODULE: horizontal_mix
!
! !DESCRIPTION:
!  This module contains driver routines for managing the individual
!  horizontal tracer and momentum mixing modules.
!
! !REVISION HISTORY:
!  SVN:$Id: horizontal_mix.F90 47361 2013-05-21 20:54:30Z mlevy@ucar.edu $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_ConstantsMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block
   use distribution, only: 
   use domain_size
   use domain, only: nblocks_clinic, distrb_clinic
   use constants, only: c0, blank_fmt, delim_fmt, ndelim_fmt
   use communicate, only: my_task, master_task
   use time_management, only: km, nt, mix_pass
   use broadcast, only: broadcast_scalar
   use grid, only: KMT, dz, partial_bottom_cells, DZT, dzr, dzwr
   use io_types, only: nml_in, nml_filename, stdout
   use hmix_del2, only: init_del2u, init_del2t, hdiffu_del2, hdifft_del2
   use hmix_del4, only: init_del4u, init_del4t, hdiffu_del4, hdifft_del4
   use hmix_gm, only: init_gm, hdifft_gm
   use hmix_aniso, only: init_aniso, hdiffu_aniso
   use topostress, only: ltopostress
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now, &
      tavg_in_which_stream, ltavg_on
   use timers, only: timer_start, timer_stop, get_timer
   use exit_mod, only: sigAbort, exit_pop, flushm
   use mix_submeso, only: init_submeso, submeso_flux, submeso_sf
   use hmix_gm_submeso_share, only: init_meso_mixing, tracer_diffs_and_isopyc_slopes
   use prognostic
   use vertical_mix

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_horizontal_mix, &
             hdiffu, hdifft, &
             iso_impvmixt_tavg

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  horizontal mixing choices
!
!-----------------------------------------------------------------------

   integer (POP_i4), parameter :: &! available choices for mixing type
      hmix_momentum_type_del2 = 1,  &
      hmix_momentum_type_del4 = 2,  &
      hmix_momentum_type_anis = 3,  &
      hmix_tracer_type_del2 = 1,    &
      hmix_tracer_type_del4 = 2,    &
      hmix_tracer_type_gm   = 3

   integer (POP_i4) ::            &
      hmix_momentum_itype,          &! users choice for type of mixing
      hmix_tracer_itype,            &! users choice for type of mixing
      tavg_HDIFT,                   &! tavg id for horizontal diffusion
      tavg_HDIFS                     ! tavg id for horizontal diffusion

   integer (POP_i4), dimension(nt) :: &
      tavg_HDIFE_TRACER,            &! tavg id for east face diffusive flux of tracer
      tavg_HDIFN_TRACER,            &! tavg id for north face diffusive flux of tracer
      tavg_HDIFB_TRACER              ! tavg id for bottom face diffusive flux of tracer

   logical (log_kind) ::            &
      lsubmesoscale_mixing           ! if true, submesoscale mixing is on


!-----------------------------------------------------------------------
!
!  timers
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      timer_hdiffu,      &! timer for horizontal momentum mixing
      timer_hdifft,      &! timer for horizontal tracer   mixing
      timer_submeso       ! timer for all submeso code

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_horizontal_mix
! !INTERFACE:

 subroutine init_horizontal_mix(errorCode)

! !DESCRIPTION:
!  Initializes choice of mixing method based on namelist input.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,        &! dummy loop index
      nml_error  ! error flag for namelist

   character (POP_charLength) ::  &! character choice for type of mixing
      hmix_momentum_choice, &
      hmix_tracer_choice

   namelist /hmix_nml/ hmix_momentum_choice, hmix_tracer_choice,  &
                       lsubmesoscale_mixing

!-----------------------------------------------------------------------
!
!  read namelist input
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   hmix_momentum_choice = 'unknown_hmix_momentum_choice'
   hmix_tracer_choice = 'unknown_hmix_tracer_choice'
   lsubmesoscale_mixing = .false.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=hmix_nml, iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading hmix_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a25)') 'Horizontal mixing options'
      write(stdout,blank_fmt)
      write(stdout,*) ' hmix_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,hmix_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)

      select case (hmix_momentum_choice(1:4))
      case ('del2')
         hmix_momentum_itype = hmix_momentum_type_del2
         write(stdout,'(a42)') &
           'Laplacian horizontal momentum mixing used.'
      case ('del4')
         hmix_momentum_itype = hmix_momentum_type_del4
         write(stdout,'(a43)') &
           'Biharmonic horizontal momentum mixing used.'
      case ('anis')
         hmix_momentum_itype = hmix_momentum_type_anis
         write(stdout,'(a44)') &
           'Anisotropic horizontal momentum mixing used.'
      case ('gent')
         hmix_momentum_itype = -1000
      case default
         hmix_momentum_itype = -2000
      end select

      select case (hmix_tracer_choice(1:4))
      case ('del2')
         hmix_tracer_itype = hmix_tracer_type_del2
         write(stdout,'(a44)') &
           'Laplacian horizontal tracer   mixing chosen.'
      case ('del4')
         hmix_tracer_itype = hmix_tracer_type_del4
         write(stdout,'(a43)') &
           'Biharmonic horizontal tracer   mixing used.'
      case ('gent')
         hmix_tracer_itype = hmix_tracer_type_gm
         write(stdout,'(a35)') &
          'Gent-McWilliams tracer mixing used.'
      case default
         hmix_tracer_itype = -1000
      end select

      if ( lsubmesoscale_mixing ) then
        write (stdout,blank_fmt)
        write (stdout, '(a48)') &
         'Submesoscale mixed layer parameterization is on.'
      endif

   endif

   call broadcast_scalar(hmix_momentum_itype, master_task)
   call broadcast_scalar(hmix_tracer_itype,   master_task)
   call broadcast_scalar(lsubmesoscale_mixing,master_task)

   if (hmix_momentum_itype == -1000) then
      call exit_POP(sigAbort, &
              'Gent-McWilliams can only be used for tracer mixing')
   else if (hmix_momentum_itype == -2000) then
      call exit_POP(sigAbort, &
                    'Unknown type for horizontal momentum mixing')
   endif

   if (hmix_tracer_itype == -1000) then
      call exit_POP(sigAbort, &
                    'Unknown type for horizontal tracer mixing')
   endif

!-----------------------------------------------------------------------
!
!  calculate additional coefficients based on mixing parameterization
!  initialize timers
!
!-----------------------------------------------------------------------

   select case (hmix_momentum_itype)
   case(hmix_momentum_type_del2)
      call init_del2u(errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_hmix: error initializing del2u')
         return
      endif

      call get_timer(timer_hdiffu,'HMIX_MOMENTUM_DEL2', &
                                  nblocks_clinic, distrb_clinic%nprocs)

   case(hmix_momentum_type_del4)
      call init_del4u(errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_hmix: error initializing del4u')
         return
      endif

      call get_timer(timer_hdiffu,'HMIX_MOMENTUM_DEL4', &
                                  nblocks_clinic, distrb_clinic%nprocs)

   case(hmix_momentum_type_anis)
      call init_aniso 
      call get_timer(timer_hdiffu,'HMIX_MOMENTUM_ANISO', &
                                  nblocks_clinic, distrb_clinic%nprocs)
   end select

   select case (hmix_tracer_itype)
   case(hmix_tracer_type_del2)
      call init_del2t(errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_hmix: error initializing del2t')
         return
      endif

      call get_timer(timer_hdifft,'HMIX_TRACER_DEL2', &
                                  nblocks_clinic, distrb_clinic%nprocs)

   case(hmix_tracer_type_del4)

      call init_del4t(errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_hmix: error initializing del4t')
         return
      endif

      call get_timer(timer_hdifft,'HMIX_TRACER_DEL4', &
                                  nblocks_clinic, distrb_clinic%nprocs)

   case(hmix_tracer_type_gm)
      call init_meso_mixing(hmix_tracer_itype,hmix_tracer_type_gm)
      call init_gm  ! variables used by GM parameterization
      call get_timer(timer_hdifft,'HMIX_TRACER_GM', &
                                  nblocks_clinic, distrb_clinic%nprocs)


   end select

!-----------------------------------------------------------------------
!
!  initialize submesoscale mixing
!
!-----------------------------------------------------------------------

   if ( lsubmesoscale_mixing )  then
        call get_timer(timer_submeso,'SUBMESO', &
                                  nblocks_clinic, distrb_clinic%nprocs)
   	call init_submeso
   	if ( .not. hmix_tracer_itype == hmix_tracer_type_gm ) then
     	 call init_meso_mixing(hmix_tracer_itype,hmix_tracer_type_gm)
	endif
   endif

  
!-----------------------------------------------------------------------
!
!  check for compatibility with topostress
!
!-----------------------------------------------------------------------

   if (ltopostress .and. &
       hmix_momentum_itype /= hmix_momentum_type_del2) then
      if (my_task == master_task) write(stdout,'(a59)') &
         'WARNING: TOPOSTRESS HAS NO EFFECT IF DEL2 MIXING NOT CHOSEN'
   endif

!-----------------------------------------------------------------------
!
!  define tavg field for tavg diagnostics
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_HDIFT,'HDIFT',2,                            &
                    long_name='Vertically Integrated Horz Mix T tendency', &
                          coordinates='TLONG TLAT time',                   &
                          units='centimeter degC/s', grid_loc='2110')

   call define_tavg_field(tavg_HDIFS,'HDIFS',2,                             &
                    long_name='Vertically Integrated Horz Diff S tendency', &
                          coordinates='TLONG TLAT time',                    &
                          scale_factor=1000.0_rtavg,                        &
                          units='centimeter gram/kilogram/s', grid_loc='2110')

   do n = 1,nt
      call define_tavg_field(tavg_HDIFE_TRACER(n),                        &
                             'HDIFE_' /&
                                       &/ trim(tracer_d(n)%short_name),3, &
                                long_name=trim(tracer_d(n)%short_name)   /&
                      &/ ' Horizontal Diffusive Flux in grid-x direction',&
                                units=trim(tracer_d(n)%tend_units),       &
                                scale_factor=tracer_d(n)%scale_factor,    &
                                grid_loc='3211',                          &
                                coordinates='ULONG TLAT z_t time' )

      call define_tavg_field(tavg_HDIFN_TRACER(n),                        &
                             'HDIFN_' /&
                                       &/ trim(tracer_d(n)%short_name),3, &
                                long_name=trim(tracer_d(n)%short_name)   /&
                      &/ ' Horizontal Diffusive Flux in grid-y direction',&
                                units=trim(tracer_d(n)%tend_units),       &
                                scale_factor=tracer_d(n)%scale_factor,    &
                                grid_loc='3121',                          &
                                coordinates='TLONG ULAT z_t time' )

      call define_tavg_field(tavg_HDIFB_TRACER(n),                        &
                             'HDIFB_' /&
                                       &/ trim(tracer_d(n)%short_name),3, &
                                long_name=trim(tracer_d(n)%short_name)   /&
                       &/ ' Horizontal Diffusive Flux across Bottom Face',&
                                units=trim(tracer_d(n)%tend_units),       &
                                scale_factor=tracer_d(n)%scale_factor,    &
                                grid_loc='3113',                          &
                                coordinates='TLONG TLAT z_w_bot time' )
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine init_horizontal_mix

!***********************************************************************
!BOP
! !IROUTINE: hdiffu
! !INTERFACE:

 subroutine hdiffu(k,HDUK,HDVK,UMIXK,VMIXK,this_block)

! !DESCRIPTION:
!  This routine returns tendencies for horizontal diffusion of
!  momentum.  It is a driver routine which simply branches to the
!  proper horizontal mix routine based on the user choice of mixing
!  method.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: k   ! depth level index

   real (POP_r8), dimension(nx_block,ny_block), intent(in) :: &
      UMIXK, VMIXK         ! U,V at level k and mix time level

   type (block), intent(in) :: &
      this_block           ! block information for this subblock

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block), intent(out) :: &
      HDUK,                   &! returned as Hdiff(U) at level k
      HDVK                     ! returned as Hdiff(V) at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  branch to the proper mix routine
!
!-----------------------------------------------------------------------

   call timer_start(timer_hdiffu, block_id=this_block%local_id)

   select case (hmix_momentum_itype)

   case (hmix_momentum_type_del2)
      call hdiffu_del2(k, HDUK, HDVK, UMIXK, VMIXK, this_block)
   case (hmix_momentum_type_del4)
      call hdiffu_del4(k, HDUK, HDVK, UMIXK, VMIXK, this_block)
   case (hmix_momentum_type_anis)
      call hdiffu_aniso(k, HDUK, HDVK, UMIXK, VMIXK, this_block)
   end select

   call timer_stop(timer_hdiffu, block_id=this_block%local_id)

!-----------------------------------------------------------------------
!EOC

 end subroutine hdiffu

!***********************************************************************
!BOP
! !IROUTINE: hdifft
! !INTERFACE:

 subroutine hdifft(k, HDTK, TMIX, UMIX, VMIX, this_block)

! !DESCRIPTION:
!  This routine returns tendencies for horizontal diffusion of
!  tracers.  It is a driver routine which simply branches to the
!  proper horizontal mix routine based on the user choice of mixing
!  method.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: k   ! depth level index

   real (POP_r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TMIX     ! tracers at mix time level

   real (POP_r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UMIX, VMIX   ! U,V velocities at mix time level

   type (block), intent(in) :: &
      this_block           ! block information for this subblock

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block,nt), intent(out) :: &
      HDTK                ! Hdiff(T) for nth tracer at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      bid                 ! local block id

   real (POP_r8), dimension(nx_block,ny_block) :: &
     WORK                 ! temporary to hold tavg field
   real (POP_r8), dimension(nx_block,ny_block,nt) :: &
      TDTK                ! Hdiff(T) for nth tracer at level k from submeso_flux code

!-----------------------------------------------------------------------
!
!  branch to the proper mix routine
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   call timer_start(timer_hdifft, block_id=bid)

   HDTK = c0

   select case (hmix_tracer_itype)
   case (hmix_tracer_type_del2)
      call hdifft_del2(k, HDTK, TMIX, tavg_HDIFE_TRACER, tavg_HDIFN_TRACER, this_block)
   case (hmix_tracer_type_del4)
      call hdifft_del4(k, HDTK, TMIX, tavg_HDIFE_TRACER, tavg_HDIFN_TRACER, this_block)
   case (hmix_tracer_type_gm)
      if (k == 1) then
	call tracer_diffs_and_isopyc_slopes(TMIX, this_block)
      endif
      call hdifft_gm(k, HDTK, TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
                     tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)
   end select
   
   call timer_stop(timer_hdifft, block_id=bid)
  
	
   if ( lsubmesoscale_mixing ) then
       call timer_start(timer_submeso, block_id=this_block%local_id)
        if (.not. hmix_tracer_itype == hmix_tracer_type_gm) then
	 if (k == 1) then
	  call tracer_diffs_and_isopyc_slopes(TMIX, this_block)
	 endif
	endif
        if (k == 1) then
	 call submeso_sf(TMIX, this_block)
	endif
        call submeso_flux(k, TDTK, TMIX, tavg_HDIFE_TRACER, &
                     tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)
	HDTK=HDTK+TDTK
       call timer_stop(timer_submeso, block_id=this_block%local_id)
   endif
   
  
   
   
!-----------------------------------------------------------------------
!
!  compute tavg diagnostic if requested
!
!-----------------------------------------------------------------------

   if (accumulate_tavg_now(tavg_HDIFT)) then
     WORK = c0
     if (partial_bottom_cells) then
        where (k <= KMT(:,:,bid)) WORK = DZT(:,:,k,bid)*HDTK(:,:,1)
     else
        where (k <= KMT(:,:,bid)) WORK = dz(k)*HDTK(:,:,1)
     endif
     call accumulate_tavg_field(WORK,tavg_HDIFT,bid,k)
   endif

   if (accumulate_tavg_now(tavg_HDIFS)) then
     WORK = c0
     if (partial_bottom_cells) then
        where (k <= KMT(:,:,bid)) WORK = DZT(:,:,k,bid)*HDTK(:,:,2)
     else
        where (k <= KMT(:,:,bid)) WORK = dz(k)*HDTK(:,:,2)
     endif
     call accumulate_tavg_field(WORK,tavg_HDIFS,bid,k)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine hdifft

!***********************************************************************

 subroutine iso_impvmixt_tavg(TNEW, bid)

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TNEW         ! on input, contains tracer to update from
                   ! on output, contains updated tracers at new time

   integer (int_kind), intent(in) :: &
      bid          ! local block address

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block) :: & 
      WORK1, WORK2

   integer (int_kind) ::  &
      k,n                  ! dummy loop indices

!-----------------------------------------------------------------------

   if (hmix_tracer_itype /= hmix_tracer_type_gm) return

   do n = 1,nt
      if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
         do k=1,km-1
            WORK1 = VDC_GM(:,:,k,bid)
            if (partial_bottom_cells) then
               WORK2 = merge(WORK1*(TNEW(:,:,k,n) - TNEW(:,:,k+1,n))/        &
                             (p5*(DZT(:,:,k,bid) + DZT(:,:,k+1,bid)))        &
                             ,c0, k < KMT(:,:,bid))*dzr(k)
            else
               WORK2 = merge(WORK1*(TNEW(:,:,k,n) - TNEW(:,:,k+1,n))*dzwr(k) &
                             ,c0, k < KMT(:,:,bid))*dzr(k)
            endif
            call accumulate_tavg_field(WORK2,tavg_HDIFB_TRACER(n),bid,k)
         end do
      endif
   end do

 end subroutine iso_impvmixt_tavg

!***********************************************************************

 end module horizontal_mix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
