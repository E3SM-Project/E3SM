!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice

!BOP
! !MODULE: ice
!
! !DESCRIPTION:
!  This module currently contains routines for computing ice 
!  formation and the heat flux associated with ice formation.
!  This heat flux is sent to the ice model via the flux coupler.
!  In the future, this module could contain the driver for a
!  subroutinized ice model.
!
! !REVISION HISTORY:
!  SVN:$Id: ice.F90 17759 2009-08-12 20:20:36Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod

   use kinds_mod, only: int_kind, log_kind, char_len, r8
   use blocks, only: nx_block, ny_block, block
   use domain_size
   use constants, only: cp_sw, latent_heat_fusion, delim_fmt, blank_fmt,     &
       sea_ice_salinity, ppt_to_salt, ocn_ref_salinity, c0, p5, hflux_factor,&
       grav, eps2, ndelim_fmt
   use broadcast, only: broadcast_scalar
   use communicate, only: my_task, master_task
   use io_types, only: nml_in, nml_filename, stdout
   use time_management, only: freq_opt_never, freq_opt_nyear, freq_opt_nday, &
       freq_opt_nhour, freq_opt_nsecond, freq_opt_nstep, init_time_flag,     &
       max_blocks_clinic, km, nt, avg_ts, back_to_back, dtt, check_time_flag,&
       partial_bottom_cells, KMT, DZT, DZ, freq_opt_nmonth, dt, tmix_matsuno,&
       tmix_iopt, ice_ts, access_time_flag
   use exit_mod, only: sigAbort, exit_pop, flushm
   use prognostic
   use passive_tracers, only: tracer_ref_val
   use grid, only: sfc_layer_varthick, sfc_layer_type

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_ice,           &
             increment_tlast_ice,&
             ice_formation,      &
             tfreez,             &
             ice_flx_to_coupler, &
             tmelt
             

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      liceform            ! flag to turn on/off ice formation

   logical (log_kind), public :: &
      lactive_ice         ! T ==> ocn is coupled to an active ice model
                          ! F ==> ocn is coupled to a dummy ice model

   integer (int_kind), public :: &
      ice_cpl_flag        ! time flag id for coupled timestep

   real (r8), public ::  &
      tlast_ice,         &! time since last ice flux computed
      time_weight,       &
      salice              ! sea ice salinity in msu

   real (r8), dimension(:,:,:), allocatable, public :: &
      AQICE,             &! sum of accumulated ice heat flux since tlast
      QFLUX,             &! internal ocn heat flux due to ice formation
      QICE                ! tot column cooling from ice form (in C*cm)

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), public ::  &
      FW_FREEZE  ! water flux at T points due to frazil ice formation

   real (r8), public ::  &
      cp_over_lhfusion    ! cp_sw/latent_heat_fusion

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ice_flag,          &! time flag id for ice formation
      kmxice              ! lowest level from which to integrate 
                          ! ice formation

   real (r8) ::          &
      salref              ! ocean ref salinity in msu


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:

 subroutine init_ice

! !DESCRIPTION:
!  This routine initializes ice variables.  It must be called
!  before initializing restarts because this module add the accumulated
!  ice heat flux to the restart file.
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

   integer (int_kind) ::   &
      nml_error             ! namelist i/o error flag

   character (char_len) :: &
      ice_freq_opt          ! option for frequency of computing ice

   integer (int_kind)   :: &
      ice_freq_iopt,       &! int option for freq units
      ice_freq              ! freq for computing ice in above units


   namelist /ice_nml/ kmxice, ice_freq_opt, ice_freq, lactive_ice

!-----------------------------------------------------------------------
!
!  read input namelists
!
!-----------------------------------------------------------------------

   cp_over_lhfusion = rho_sw*cp_sw/(latent_heat_fusion*rho_fw)

   kmxice           = 1
   ice_freq_opt     = 'never'
   ice_freq         = 100000
   lactive_ice      = .false.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=ice_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading ice_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) 'Ice:'
      write(stdout,blank_fmt)
      write(stdout,*) ' ice_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,ice_nml)
      write(stdout,delim_fmt)

      !***
      !*** define salice and salref in msu
      !***
 
      salice = sea_ice_salinity*ppt_to_salt
      salref = ocn_ref_salinity*ppt_to_salt

      liceform = .true.

      select case (ice_freq_opt)
      case ('never')
         write(stdout,'(a22)') 'Ice formation disabled'
         liceform = .false.
         ice_freq_iopt = freq_opt_never
      case ('coupled')
         write(stdout,'(a44)') &
           'Ice formation computed on coupler time steps'
         ice_freq_iopt = freq_opt_never ! check coupler flag instead
      case ('nyear')
         write(stdout,'(a29,i4,a9)') 'Ice formation computed every ', &
                                      ice_freq,' years.  '
         ice_freq_iopt = freq_opt_nyear
      case ('nmonth')
         write(stdout,'(a29,i4,a9)') 'Ice formation computed every ', &
                                      ice_freq,' months. '
         ice_freq_iopt = freq_opt_nmonth
      case ('nday')
         write(stdout,'(a29,i4,a9)') 'Ice formation computed every ', &
                                      ice_freq,' days.   '
         ice_freq_iopt = freq_opt_nday
      case ('nhour')
         write(stdout,'(a29,i4,a9)') 'Ice formation computed every ', &
                                      ice_freq,' hours.  '
         ice_freq_iopt = freq_opt_nhour
      case ('nsecond')
         write(stdout,'(a29,i4,a9)') 'Ice formation computed every ', &
                                      ice_freq,' seconds.'
         ice_freq_iopt = freq_opt_nsecond
      case ('nstep')
         write(stdout,'(a29,i4,a9)') 'Ice formation computed every ', &
                                      ice_freq,' steps.  '
         ice_freq_iopt = freq_opt_nstep
      case default
         ice_freq_iopt = -1000
      end select

      if (liceform) then
         write(stdout,'(a20,1pe10.3)') 'Ice salinity(ppt) = ', &
                                        sea_ice_salinity
         write(stdout,'(a30,i3,a13)') 'Ice formation computed in top ', &
                                       kmxice, ' levels only.'
      endif
      call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)

   endif

   call broadcast_scalar(ice_freq_iopt, master_task)
   call broadcast_scalar(lactive_ice,   master_task)

   if (ice_freq_iopt == -1000) then
      call exit_POP(sigAbort,'unknown restart frequency option')
   endif

   call broadcast_scalar(liceform,      master_task)

!-----------------------------------------------------------------------
!
!  if ice turned on, broadcast remaining vars and allocate memory
!
!-----------------------------------------------------------------------

   if (liceform) then
      call broadcast_scalar(ice_freq, master_task)
      call broadcast_scalar(kmxice,   master_task)
      call broadcast_scalar(salice,   master_task)
      call broadcast_scalar(salref,   master_task)
 
      !***
      !*** set up ice time flag and get coupled_ts id for local use
      !***

      call init_time_flag('ice',ice_flag, default=.false.,  &
                           owner    = 'init_ice',           &
                           freq_opt = ice_freq_iopt,        &
                           freq     = ice_freq)

      call access_time_flag('coupled_ts', ice_cpl_flag)

      !***
      !*** must keep track of time since last ice flux computed
      !***

      tlast_ice = c0

      !***
      !*** allocate and initialize ice flux arrays
      !***

      allocate( QICE(nx_block,ny_block,max_blocks_clinic), & 
               AQICE(nx_block,ny_block,max_blocks_clinic), &
               QFLUX(nx_block,ny_block,max_blocks_clinic))

      QICE  = c0
      AQICE = c0
      QFLUX = c0

   endif

   FW_FREEZE = c0


!-----------------------------------------------------------------------
!EOC

 call flushm (stdout)

 end subroutine init_ice

!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:
   subroutine increment_tlast_ice

! !DESCRIPTION:
!  This subroutine increments tlast_ice in a nonthreaded region.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  increment time since last evaluation
!
!-----------------------------------------------------------------------

   if (avg_ts .or. back_to_back) then
     time_weight = p5 
   else
     time_weight = c1 
   endif

   tlast_ice = tlast_ice + dtt*time_weight

!-----------------------------------------------------------------------
!EOC

 end subroutine increment_tlast_ice

!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:
   subroutine ice_formation(TNEW, SHF_IN, iblock,this_block,lfw_as_salt_flx)

! !DESCRIPTION:
!  This subroutine computes ocean heat flux to the sea-ice. it forms
!  the necessary ice in the ocean and adjusts the potential 
!  temperature and salinity fields accordingly. the logic of this 
!  subroutine is based on William Large''s 1-d model and is based
!  on a version from the NCOM model written by Gokhan Danabasoglu. 


! !REVISION HISTORY:
!  same as module
 
! !INPUT/OUTPUT PARAMETERS:
 
   real (r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
      TNEW                ! tracers at new time level
 
! !INPUT PARAMETERS:
 
   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SHF_IN              ! surface heat flux

   integer (int_kind), intent(in) :: &
      iblock              ! block information for current block

   type (block), intent(in) :: &
      this_block          ! block information for current block
 

   logical (log_kind), intent(in) :: &
      lfw_as_salt_flx         


!EOP
!BOC


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,                 &! vertical level index
      bid,               &! local block index
      n                   ! tracer index


   real (r8), dimension(nx_block,ny_block) :: &
     POTICE,           & ! potential amt of ice formation
     WORK1,            & ! work array
     TFRZ                ! freezing temp of water


   real (r8) ::        &
     ref_val             ! tracer reference value


!-----------------------------------------------------------------------
!
!  presently, the check_time_flag does not produce the same results
!  as ice_ts.  for now, ice_ts will be set in time_management and
!  used here
!-----------------------------------------------------------------------
!=====if (check_time_flag(ice_flag) .or.
!====&    check_time_flag(ice_cpl_flag)) then
 
!========================
   if ( ice_ts ) then
!========================

     bid = this_block%local_id

!-----------------------------------------------------------------------
!
!  initialize flux to zero
!
!-----------------------------------------------------------------------
 
     QICE(:,:,bid) = c0
     POTICE = c0

 
!-----------------------------------------------------------------------
!
!  compute frazil ice formation for sub-surface layers. if ice
!  forms in lower layers but layers above are warm - the heat is
!  used to melt the ice. the ice formation occurs at salinity, Si.
!  this volume is replaced with an equal volume at the salinity of
!  the layer above. the total ice heat flux is accumulated. 
!
!  WARNING: unless a monotone advection scheme is in place, 
!  advective errors could lead to temps that are far below freezing
!  in some locations and this scheme will form lots of ice.
!  ice formation should be limited to the top layer (kmxice=1)
!  if the advection scheme is not monotone.
!
!-----------------------------------------------------------------------
 
     do k=kmxice,2,-1

     !***
     !*** potice is the potential amount of ice formation 
     !*** (potice>0) or melting (potice<0) in layer k
     !***

       call tfreez(TFRZ,TNEW(:,:,k,2))
       if (partial_bottom_cells) then
          where (k <= KMT(:,:,bid)) &
             POTICE = (TFRZ - TNEW(:,:,k,1))*DZT(:,:,k,bid)
       else
          where (k <= KMT(:,:,bid)) &
             POTICE = (TFRZ - TNEW(:,:,k,1))*dz(k)
       endif

     !***
     !*** if potice < 0, use the heat to melt any ice
     !*** from lower layers
     !*** if potice > 0, keep on freezing (QICE < 0)
     !***

       POTICE = max(POTICE,QICE(:,:,bid))

     !***
     !*** adjust tracer values based on freeze/melt
     !***

       if (partial_bottom_cells) then
         TNEW(:,:,k,1) = TNEW(:,:,k,1) + POTICE/DZT(:,:,k,bid)
       else
         TNEW(:,:,k,1) = TNEW(:,:,k,1) + POTICE/dz(k)
       endif


       if (sfc_layer_type == sfc_layer_varthick .and. .not. lfw_as_salt_flx) then
         if (partial_bottom_cells) then
            TNEW(:,:,k,2) = ( TNEW(:,:,k,2)                               &
            * (DZT(:,:,k,bid) + cp_over_lhfusion * QICE(:,:,bid))         &
            + cp_over_lhfusion * (TNEW(:,:,k-1,2)                         &
            * (POTICE - QICE(:,:,bid)) - salice * POTICE) )/DZT(:,:,k,bid)
         else
            TNEW(:,:,k,2) = ( TNEW(:,:,k,2)                               &
            * (dz(k) + cp_over_lhfusion * QICE(:,:,bid))                  &
            + cp_over_lhfusion * (TNEW(:,:,k-1,2)                         &
            * (POTICE - QICE(:,:,bid)) - salice * POTICE) )/dz(k)
         endif
       else
          do n=2,nt
             ref_val = salref - salice
             if (n > 2)  &
                ref_val = ref_val * (tracer_ref_val(n) / salref)
             if (ref_val /= c0)  then
                if (partial_bottom_cells) then
                  TNEW(:,:,k,n) = TNEW(:,:,k,n) &
                  + ref_val*POTICE*cp_over_lhfusion/DZT(:,:,k,bid)
                else
                  TNEW(:,:,k,n) = TNEW(:,:,k,n) &
                  + ref_val*POTICE*cp_over_lhfusion/dz(k)
                endif
             endif
          end do
       endif

       !*** accumulate freezing potential
       QICE(:,:,bid) = QICE(:,:,bid) - POTICE

     enddo ! k loop

!-----------------------------------------------------------------------
!
!  now repeat the above algorithm for the surface layer. when fresh
!  water flux formulation is used, the surface layer does not get
!  any salt from other layers. instead, its volume changes. 
!
!-----------------------------------------------------------------------

     k = 1
      
     call tfreez(TFRZ,TNEW(:,:,k,2))

     if (partial_bottom_cells) then
       WORK1 = DZT(:,:,k,bid)
     else
       WORK1 = dz(k)
     endif

     if (sfc_layer_type == sfc_layer_varthick)  &
       WORK1 = WORK1 + PSURF(:,:,newtime,bid)/grav + eps2
     where ( k <= KMT(:,:,bid) )
       POTICE = (TFRZ - TNEW(:,:,k,1))*WORK1
     endwhere 

     POTICE = max(POTICE, QICE(:,:,bid))

     TNEW(:,:,k,1) = TNEW(:,:,k,1) + POTICE/WORK1
     if (sfc_layer_type == sfc_layer_varthick .and. .not. lfw_as_salt_flx) then
       TNEW(:,:,k,2) =  &
          (TNEW(:,:,k,2)*(WORK1 + cp_over_lhfusion*QICE(:,:,bid)) - &
           salice*QICE(:,:,bid)*cp_over_lhfusion )/WORK1
     else
        do n=2,nt 
           ref_val = salref - salice
           if (n > 2) ref_val = ref_val * (tracer_ref_val(n) / salref)
           if (ref_val /= c0)  &
              TNEW(:,:,k,n) = TNEW(:,:,k,n)  &
              + ref_val * POTICE * cp_over_lhfusion / WORK1
        end do
     endif

     QICE(:,:,bid) = QICE(:,:,bid) - POTICE

!-----------------------------------------------------------------------
!
!  let any residual heat in the upper layer melt previously formed ice
!
!-----------------------------------------------------------------------
 
     AQICE(:,:,bid) = AQICE(:,:,bid) + time_weight*QICE(:,:,bid)

!-----------------------------------------------------------------------
!
!  recalculate freezing potential based on adjusted T.
!  only interested in melt potential now (POTICE < 0) - use this 
!  melt to offset any accumulated freezing (AQICE < 0) and
!  adjust T and S to reflect this melting. when freshwater flux
!  formulation, compute the associated freshwater flux instead of
!  adjusting S.
!
!-----------------------------------------------------------------------

     call tfreez(TFRZ,TNEW(:,:,k,2))

     where (k <= KMT(:,:,bid))
       POTICE = (TFRZ - TNEW(:,:,k,1)) * WORK1 
     endwhere

     POTICE = max(POTICE, AQICE(:,:,bid))

     TNEW(:,:,k,1) = TNEW(:,:,k,1) + POTICE/WORK1

     if ( sfc_layer_type == sfc_layer_varthick .and. .not. lfw_as_salt_flx ) then
       FW_FREEZE(:,:,bid) = time_weight * min(POTICE(:,:),QICE(:,:,bid)) &
                  * cp_over_lhfusion / dt(k)
     else
       do n=2,nt
          ref_val = salref - salice
          if (n > 2) ref_val = ref_val * (tracer_ref_val(n)/salref)
          if (ref_val /= c0) &
             TNEW(:,:,k,n) = TNEW(:,:,k,n)  &
             + ref_val*POTICE*cp_over_lhfusion/WORK1 
       end do
     endif

     AQICE(:,:,bid) = AQICE(:,:,bid) - time_weight*POTICE

   endif ! time to do ice

!-----------------------------------------------------------------------
!EOC

   end subroutine ice_formation

!***********************************************************************

   subroutine ice_flx_to_coupler(TCUR,bid)

!-----------------------------------------------------------------------
!
!  This subroutine sets up the ice formation / melting potential
!  heat fluxes to be sent to the coupler. ice formation heat flux
!  is accumulated for time averaging.
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TCUR                ! tracers at new time level

   integer (int_kind), intent(in) :: &
      bid                 ! local block index

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      k                    ! vertical level index


   real (r8), dimension(nx_block,ny_block) :: & 
      WORK1, WORK2         ! work arrays

   real (r8), dimension(nx_block,ny_block) :: &
      TFRZ                 ! freezing temp of water

!-----------------------------------------------------------------------
!
!  compute the first layer thickness
!
!-----------------------------------------------------------------------

   k = 1 
   call tfreez(TFRZ,TCUR(:,:,k,2))

   if (partial_bottom_cells) then
     WORK1 = DZT(:,:,k,bid)
   else
     WORK1 = dz(k)
   endif 

   if ( sfc_layer_type == sfc_layer_varthick ) &
     WORK1 = WORK1 + PSURF(:,:,curtime,bid)/grav + eps2

!-----------------------------------------------------------------------
!
!  first compute the melt potential
!
!-----------------------------------------------------------------------

   WORK2 = c0
   where ( k <= KMT(:,:,bid) )
     WORK2 = (TFRZ - TCUR(:,:,k,1)) * WORK1
   endwhere

!-----------------------------------------------------------------------
!
!  adjust ice formation amount when mixing step is tavg
!
!-----------------------------------------------------------------------

   if ( tmix_iopt /= tmix_matsuno ) then
     AQICE(:,:,bid) = p5 * AQICE(:,:,bid)
   endif

!-----------------------------------------------------------------------
!
!  merge the ice formation and melt potential fluxes and compute
!
!-----------------------------------------------------------------------

   where ( AQICE(:,:,bid) < c0 ) 
     QICE(:,:,bid) = -AQICE(:,:,bid)
   elsewhere
     QICE(:,:,bid) = WORK2
   endwhere

   if (tlast_ice == c0) then
     QFLUX(:,:,bid) = c0
   else
     QFLUX(:,:,bid) = QICE(:,:,bid)/tlast_ice/hflux_factor
   endif

   end subroutine ice_flx_to_coupler


!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:

 subroutine tfreez(TFRZ,SALT)

! !DESCRIPTION:
!  This function computes the freezing point of salt water.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SALT                ! salinity in model units (g/g)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      TFRZ                ! freezing temperature of water in deg C

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  use only the first salinity term in the expansion
!
!-----------------------------------------------------------------------

   !TFRZ = -0.0544_r8*SALT*salt_to_ppt
   TFRZ = -1.8_r8

!-----------------------------------------------------------------------
!EOC

 end subroutine tfreez


!BOP
! !IROUTINE:
! !INTERFACE:

   subroutine tmelt (TMLT,SALT)

! !DESCRIPTION:
!  This subroutine sets the melting point temperature of ice.
!  For now, TMLT is a separate routine than TFRZ.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) ::  &
     SALT                ! salinity in model units (g/g)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: &
     TMLT                ! melting temperature in deg C
!EOP
!BOC

   if ( lactive_ice ) then
     TMLT = c0
   else
     call tfreez(TMLT,SALT)
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine tmelt


!***********************************************************************

 end module ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
