!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  module vertical_mix

!BOP
! !MODULE: vertical_mix
! !DESCRIPTION:
!  This module contains the routines for computing vertical mixing 
!  tendencies, implicit vertical mixing and convection.  Routines for 
!  computing the diffusion coefficients themselves are in separate 
!  modules for the specific parameterizations.
!  Currently, the options for vertical mixing include:
!  \begin{itemize}
!    \item constant: a constant (in space and time) vertical diffusion
!    \item Rich no:  vertical diffusion based on local Richardson number
!                    (e.g. Pacanowski and Philander)
!    \item KPP:      k-profile vertical mixing (from Large et al.)
!  \end{itemize}
!
! !REVISION HISTORY:
!  SVN:$Id: vertical_mix.F90 47279 2013-05-18 00:14:34Z mlevy@ucar.edu $

! !USES:

   use kinds_mod
   use blocks
   use distribution
   use domain
   use constants
   use grid
   use state_mod
   use broadcast
   use io
   use timers
   use tavg
   use time_management
   use diagnostics
   use vmix_const
   use vmix_rich
   use vmix_kpp
   use exit_mod
   use prognostic

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_vertical_mix,                    &
             vmix_coeffs,                          &
             vdifft, vdiffu,                       &
             impvmixt, impvmixt_correct, impvmixu, &
             convad, impvmixt_tavg

! !PUBLIC DATA MEMBERS:

   real (r8), dimension(:,:,:,:,:), allocatable, public, target :: &
      VDC                 ! tracer diffusivity - public to allow
                          ! possible modification by Gent-McWilliams
                          ! horizontal mixing parameterization

   real (r8), dimension(:,:,:,:), allocatable, public, target :: &
      VDC_GM              ! Gent-McWilliams contribution to VDC

   integer (int_kind), parameter, public :: &
      vmix_type_const = 1,  & ! integer identifiers for desired
      vmix_type_rich  = 2,  & ! mixing parameterization
      vmix_type_kpp   = 3

   integer (int_kind), public :: &
      vmix_itype              ! users choice for vmix parameterization

   logical (log_kind), public :: &
      implicit_vertical_mix   ! flag for computing vertical mixing
                              ! implicitly in time

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  vertical mixing quantities
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable, target :: &
      VVC                 ! momentum viscosity

!-----------------------------------------------------------------------
!
!  convection variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      convection_itype,  &! user choice of method to use for convection
      nconvad             ! number of times to convectively adjust

   real (r8) ::          &
      convect_diff,      &! tracer diffusion to use for convection
      convect_visc        ! viscosity to use for convection (momentum)

   integer (int_kind), parameter :: &
      convect_type_adjust = 1,      &! ids for convection choice
      convect_type_diff   = 2

!-----------------------------------------------------------------------
!
!  vectors needed for implicit vertical mixing
!
!-----------------------------------------------------------------------

   real (r8) ::          &
      aidif               ! time-centering parameter for implicit vmix

   real (r8), dimension(:), allocatable :: & 
      afac_u, afac_t

!-----------------------------------------------------------------------
!
!  these arrays store flux information to avoid re-computing fluxes
!  for each k level
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic) :: & 
      VTF        ! vertical tracer flux at top of T-box
 
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: & 
      VUF,VVF    ! vertical momentum fluxes at top of box

!-----------------------------------------------------------------------
!
!  bottom drag and bottom (geothermal) heat flux
!
!-----------------------------------------------------------------------

   logical (log_kind) ::    &
      lbottom_heat_flx       ! true if bottom heat flux /= 0

   real (r8) ::             &
      bottom_drag,          &! drag coefficient for bottom drag
      bottom_heat_flx,      &! bottom (geothermal) heat flux (W/m2)
      bottom_heat_flx_depth  ! depth below which bottom heat applied

!-----------------------------------------------------------------------
!
!  timer ids and tavg ids for time-averaged vertical mix diagnostics
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      timer_vmix_coeffs, &! timer number for vertical mix coeffs
      timer_vdiffu,      &! timer number for explicit momentum vmix
      timer_vdifft,      &! timer number for explicit tracer   vmix
      timer_impvmixt,    &! timer number for implicit tracer   vmix
      timer_impvmixu      ! timer number for implicit momentum vmix

   integer (int_kind) :: &
      tavg_VDC_T,        &! tavg id for total vertical TEMP diffusivity
      tavg_VDC_S,        &! tavg id for total vertical SALT diffusivity
      tavg_VVC,          &! tavg id for total vertical momentum viscosity
      tavg_VUF,          &! tavg id for vertical flux of U momentum
      tavg_VVF,          &! tavg id for vertical flux of V momentum
      tavg_PEC,          &! tavg id for pot energy release convection
      tavg_NCNV           ! tavg id for number of convective adjustments

   integer (int_kind), dimension(nt) :: &
      tavg_DIA_IMPVF_TRACER  ! tavg id for diabatic implicit vertical flux of tracer

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_vertical_mix
! !INTERFACE:

 subroutine init_vertical_mix

! !DESCRIPTION:
!  Initializes various mixing quantities and calls initialization
!  routines for specific parameterizations.
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

   integer (int_kind) ::  &
      k,                  &! vertical level index
      n,                  &! dummy loop index
      nu,                 &! i/o unit
      nml_error            ! namelist i/o error flag

   character (char_len) :: &
      vmix_choice          ! input choice for desired parameterization

   character (char_len) :: &
      convection_type      ! input choice for method for convection

   namelist /vertical_mix_nml/ vmix_choice, aidif, bottom_drag,        &
                               bottom_heat_flx, bottom_heat_flx_depth, &
                               implicit_vertical_mix,                  &
                               convection_type, nconvad,               &
                               convect_diff, convect_visc

!-----------------------------------------------------------------------
!
!  read input namelist and set mixing type
!
!-----------------------------------------------------------------------

   aidif = c1
   bottom_drag = 1.0e-3_r8
   bottom_heat_flx = c0
   bottom_heat_flx_depth = 1000.00e2_r8
   implicit_vertical_mix = .true.
   convection_type = 'diffusion'
   nconvad = 2
   convect_diff = 1000.0_r8
   convect_visc = 1000.0_r8

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=vertical_mix_nml, iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort, &
                    'ERROR reading vertical_mix_nml') 
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a27)') 'Vertical mixing options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)

      if (implicit_vertical_mix) then
         write(stdout,'(a33)') ' Implicit vertical mixing enabled'
      else
         write(stdout,'(a34)') ' Implicit vertical mixing disabled'
      endif

      select case (vmix_choice)

      case ('const')
         vmix_itype = vmix_type_const
         write(stdout,'(a32)') '  Using constant vertical mixing'

      case ('rich')
         vmix_itype = vmix_type_rich
         write(stdout,'(a32)') '  Using Richardson number mixing'

      case ('kpp')
         vmix_itype = vmix_type_kpp
         write(stdout,'(a32)') '  Using KPP vertical mixing     '

      case default
         vmix_itype = -1000
      end select

      write(stdout,'(a17,2x,1pe12.5)') '  aidif         =',aidif
      write(stdout,'(a17,2x,1pe12.5)') '  bottom_drag   =',bottom_drag

      select case (convection_type)
      case ('adjustment')
         convection_itype = convect_type_adjust
         write(stdout,'(a30)') ' Convective adjustment enabled'
         write(stdout,'(a7,i4,a7)') ' Using ',nconvad,' passes'
      case ('diffusion')
         convection_itype = convect_type_diff
         write(stdout,'(a41)') ' Convection simulated by strong diffusion'
         write(stdout,'(a18,1pe10.3)') '   convect_diff = ',convect_diff
         write(stdout,'(a18,1pe10.3)') '   convect_visc = ',convect_visc
      case default
         convection_itype = -1000
      end select

      if (bottom_heat_flx /= c0) then
         write(stdout,'(a36)') ' Bottom geothermal heat flux enabled'
         write(stdout,'(a17,1pe10.3)') '   Flux (W/m2) = ', &
                                       bottom_heat_flx
         write(stdout,'(a17,1pe10.3)') '   Depth (cm)  = ', &
                                       bottom_heat_flx_depth
      else
         write(stdout,'(a37)') ' Bottom geothermal heat flux disabled'
      endif
   endif

   call broadcast_scalar(vmix_itype,            master_task)
   call broadcast_scalar(aidif,                 master_task)
   call broadcast_scalar(bottom_drag,           master_task)
   call broadcast_scalar(bottom_heat_flx,       master_task)
   call broadcast_scalar(bottom_heat_flx_depth, master_task)
   call broadcast_scalar(implicit_vertical_mix, master_task)
   call broadcast_scalar(convection_itype,      master_task)
   call broadcast_scalar(nconvad,               master_task)
   call broadcast_scalar(convect_diff,          master_task)
   call broadcast_scalar(convect_visc,          master_task)
  
   bottom_heat_flx = bottom_heat_flx * hflux_factor

   if (bottom_heat_flx /= c0) then
      lbottom_heat_flx = .true.
   else
      lbottom_heat_flx = .false.
   endif

!-----------------------------------------------------------------------
!
!  do some error and consistency checks
!
!-----------------------------------------------------------------------

   if (vmix_itype == -1000) then
      call exit_POP(sigAbort,'Unknown vertical mixing choice')
   endif

   if (convection_itype == -1000) then
      call exit_POP(sigAbort, &
                    'ERROR: Unknown option for convection type')
   endif

   if (implicit_vertical_mix .and. &
       convection_itype == convect_type_adjust) then
      call exit_POP(sigAbort, &
          'Convective adjustment should not be used with implicit vmix')
   endif

   if (.not. implicit_vertical_mix .and. &
       convection_itype == convect_type_diff) then
      call exit_POP(sigAbort, &
           'Implicit vertical mixing required for diffusive convection')
   endif

   if (.not. implicit_vertical_mix .and. &
       vmix_itype == vmix_type_kpp) then
      call exit_POP(sigAbort, &
                    'Implicit vertical mixing required for KPP')
   endif

!-----------------------------------------------------------------------
!
!  allocate VDC, VVC arrays and define options for chosen 
!  parameterization
!
!-----------------------------------------------------------------------

   select case (vmix_itype)

   case(vmix_type_const)
      if (implicit_vertical_mix) then
         allocate (VDC(nx_block,ny_block,km,1,nblocks_clinic),  &
                   VVC(nx_block,ny_block,km,  nblocks_clinic))
      else
         allocate (VDC(nx_block,ny_block,1,1,nblocks_clinic),  &
                   VVC(nx_block,ny_block,1  ,nblocks_clinic))
      endif
      call init_vmix_const(VDC,VVC)
      call get_timer(timer_vmix_coeffs,'VMIX_COEFFICIENTS_CONSTANT', &
                                  nblocks_clinic, distrb_clinic%nprocs)

   case(vmix_type_rich)
      if (implicit_vertical_mix) then
         allocate (VDC(nx_block,ny_block,km,1,nblocks_clinic), &
                   VVC(nx_block,ny_block,km  ,nblocks_clinic))
      else
         !*** note that VVC must be 3-d because it is used
         !*** in a different k-loop from where it is defined
         allocate (VDC(nx_block,ny_block,1,1,nblocks_clinic), &
                   VVC(nx_block,ny_block,km ,nblocks_clinic))
      endif
      call init_vmix_rich(VDC,VVC)
      call get_timer(timer_vmix_coeffs,'VMIX_COEFFICIENTS_RICH', &
                                  nblocks_clinic, distrb_clinic%nprocs)

   case(vmix_type_kpp)
      allocate (VDC(nx_block,ny_block,0:km+1,2,nblocks_clinic), &
                VVC(nx_block,ny_block,km,      nblocks_clinic))
      call init_vmix_kpp(VDC,VVC)
      call get_timer(timer_vmix_coeffs,'VMIX_COEFFICIENTS_KPP', &
                                  nblocks_clinic, distrb_clinic%nprocs)

   end select

   call get_timer(timer_vdifft,'VMIX_EXPLICIT_TRACER', &
                               nblocks_clinic, distrb_clinic%nprocs)
   call get_timer(timer_vdiffu,'VMIX_EXPLICIT_MOMENTUM', &
                               nblocks_clinic, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!
!  set up implicit coefficients if implicit vertical mixing chosen
!
!-----------------------------------------------------------------------

   if (implicit_vertical_mix) then

      call get_timer(timer_impvmixt,'VMIX_IMPLICIT_TRACER', &
                                   nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_impvmixu,'VMIX_IMPLICIT_MOMENTUM', &
                                   nblocks_clinic, distrb_clinic%nprocs)

      allocate(afac_u(km), afac_t(km))
      do k=1,km
         afac_u(k) = aidif*dzwr(k)
         afac_t(k) = aidif*dzwr(k)
      enddo

   endif

!-----------------------------------------------------------------------
!
!  define fields for accumulating tavg diagnostics
!
!-----------------------------------------------------------------------

   if (vmix_itype == vmix_type_kpp) then
       call define_tavg_field(tavg_VDC_T,'VDC_T',3,             &
          long_name='total diabatic vertical TEMP diffusivity', &
          units='cm^2/s', grid_loc='3113',                      &
          coordinates='TLONG TLAT z_w_bot time')

       call define_tavg_field(tavg_VDC_S,'VDC_S',3,             &
          long_name='total diabatic vertical SALT diffusivity', &
          units='cm^2/s', grid_loc='3113',                      &
          coordinates='TLONG TLAT z_w_bot time')

       call define_tavg_field(tavg_VVC,'VVC',3,                 &
          long_name='total vertical momentum viscosity',        &
          units='cm^2/s', grid_loc='3113',                      &
          coordinates='TLONG TLAT z_w_bot time')
   endif

   call define_tavg_field(tavg_VUF,'VUF',3,                     &
                          long_name='Zonal viscous stress',     &
                          units='    ', grid_loc='3222')

   call define_tavg_field(tavg_VVF,'VVF',3,                      &
                          long_name='Meridional viscous stress', &
                          units='    ', grid_loc='3222')

   if (convection_itype == convect_type_adjust) then
     call define_tavg_field(tavg_PEC,'PEC',3,                          &
                      long_name='Potential energy release convection', &
                            units='g/cm^3', grid_loc='3111')

     call define_tavg_field(tavg_NCNV,'NCNV',3,                        &
                        long_name='Convective adjustments per second', &
                            units='adj/s', grid_loc='3111')
   endif

   do n = 1,nt
     call define_tavg_field(tavg_DIA_IMPVF_TRACER(n), 'DIA_IMPVF_'     /&
                                     &/ trim(tracer_d(n)%short_name),3, &
        long_name=trim(tracer_d(n)%short_name)                         /&
  &/ ' Flux Across Bottom Face from Diabatic Implicit Vertical Mixing', &
        units=trim(tracer_d(n)%flux_units), grid_loc='3113',            &
        scale_factor=tracer_d(n)%scale_factor,                          &
        coordinates='TLONG TLAT z_w_bot time')
   enddo

!-----------------------------------------------------------------------
!EOC

 call flushm(stdout)

 end subroutine init_vertical_mix

!***********************************************************************
!BOP
! !IROUTINE: vmix_coeffs
! !INTERFACE:

 subroutine vmix_coeffs(k, TMIX, UMIX, VMIX, UCUR, VCUR, RHOMIX, &
                           STF, SHF_QSW,                         &
                           this_block,                           &
                           SMF, SMFT)


! !DESCRIPTION:
!  This is a driver routine which calls the appropriate
!  parameterization to compute the mixing coefficients VDC, VVC.
!  This routine must be called successively for $k=1,2,3$...
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                      ! vertical level index

   real (r8), dimension(nx_block,ny_block,km,nt)  :: TMIX

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UMIX, VMIX,           &! U,V     at mix time level
      RHOMIX                 ! density at mix time level

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UCUR, VCUR             ! U,V     at cur time level

   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      STF                    ! surface forcing for all tracers

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SHF_QSW                ! short-wave forcing

   real (r8), dimension(nx_block,ny_block,2), intent(in), optional :: &
      SMF,                  &! surface momentum forcing at U points
      SMFT                   ! surface momentum forcing at T points
                            ! *** must pass either one or the other

   type (block), intent(in) :: &
      this_block             ! block info for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid,               &! local block address
      kk                  ! vertical level index

!-----------------------------------------------------------------------
!
!  determine which block we are working on
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   call timer_start(timer_vmix_coeffs, block_id=bid)

!-----------------------------------------------------------------------
!
!  call proper routine based on vmix choice
!
!-----------------------------------------------------------------------

   select case(vmix_itype)

   case(vmix_type_const)
      if (convection_itype == convect_type_diff) then
         call vmix_coeffs_const(k,VDC(:,:,:,:,bid),          &
                                  VVC(:,:,:,  bid),TMIX,     &
                                 this_block, convect_diff, convect_visc)
      else
         call vmix_coeffs_const(k,VDC(:,:,:,:,bid),          &
                                  VVC(:,:,:,  bid),TMIX,     &
                                  this_block)
      endif

   case(vmix_type_rich)
      if (convection_itype == convect_type_diff) then
         call vmix_coeffs_rich(k,VDC(:,:,:,:,bid),           &
                                 VVC(:,:,:,  bid),           &
                                 TMIX,UMIX,VMIX,RHOMIX,      &
                                 this_block, convect_diff, convect_visc)
      else
         call vmix_coeffs_rich(k,VDC(:,:,:,:,bid),           &
                                 VVC(:,:,:,  bid),           &
                                 TMIX,UMIX,VMIX,RHOMIX,      &
                                 this_block)
      endif

   case(vmix_type_kpp)
      if (k == 1) then     ! for KPP, compute coeffs for all levels
         if (.not. present(SMFT) .and. .not. present(SMF)) &
            call exit_POP(sigAbort, &
                    'vmix_coeffs: must supply either SMF,SMFT')

         if (present(SMFT)) then
            call vmix_coeffs_kpp(VDC(:,:,:,:,bid),           &
                                 VVC(:,:,:,  bid),           &
                                 TMIX,UMIX,VMIX,UCUR,VCUR,RHOMIX, &
                                 STF,SHF_QSW,                &
                                 this_block,                 &
                                 convect_diff, convect_visc, &
                                 SMFT=SMFT)
         else
            call vmix_coeffs_kpp(VDC(:,:,:,:,bid),           &
                                 VVC(:,:,:,  bid),           &
                                 TMIX,UMIX,VMIX,UCUR,VCUR,RHOMIX, &
                                 STF,SHF_QSW,                &
                                 this_block,                 &
                                 convect_diff, convect_visc, &
                                 SMF=SMF)
         endif

         if (accumulate_tavg_now(tavg_VDC_T) ) then
            do kk=1,km
              ! kk index is no longer shifted because z_w_bot is now an output coordinate
              ! VDC is at cell bottom
              call accumulate_tavg_field(VDC(:,:,kk,1,bid),tavg_VDC_T,bid,kk)
            end do
         endif

         if (accumulate_tavg_now(tavg_VDC_S) ) then
            do kk=1,km
              ! kk index is no longer shifted because z_w_bot is now an output coordinate
              ! VDC is at cell bottom
              call accumulate_tavg_field(VDC(:,:,kk,2,bid),tavg_VDC_S,bid,kk)
            end do
         endif

         if (accumulate_tavg_now(tavg_VVC) ) then
            do kk=1,km
              ! kk index is no longer shifted because z_w_bot is now an output coordinate
              ! VVC is at cell bottom
              call accumulate_tavg_field(VVC(:,:,kk,bid),tavg_VVC,bid,kk)
            end do
         endif

      endif

   end select
     
   call timer_stop(timer_vmix_coeffs, block_id=bid)

!-----------------------------------------------------------------------
!
!  call cfl diagnostics if not implicit mixing
!
!-----------------------------------------------------------------------

   if (.not. implicit_vertical_mix) then
      call cfl_vdiff(k, VDC(:,:,:,:,bid), VVC(:,:,:,bid), this_block)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine vmix_coeffs

!***********************************************************************
!BOP
! !IROUTINE: vdifft
! !INTERFACE:

 subroutine vdifft(k, VDTK, TOLD, STF, this_block)

! !DESCRIPTION:
!  Computes vertical diffusion of tracers
!  \begin{equation}
!   D_V(\varphi) = \delta_z(\kappa\delta_z\varphi)
!  \end{equation}
!  \begin{eqnarray}
!   \kappa\delta_z\varphi = Q_\varphi  &{\rm at }& z=0 \\
!   \kappa\delta_z\varphi = 0       &{\rm at }& z=-H_T
!  \end{eqnarray}
!
!  This routine must be called successively with $k = 1,2,3,$...
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k   ! vertical level

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TOLD                ! tracers at old time level

   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      STF                 ! surface forcing for all tracers

   type (block), intent(in) :: &
      this_block          ! block info for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(out) :: & 
      VDTK                ! returns VDIFF(Tracer(:,:,k,n))

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                  &! dummy loop counter for tracer number
      bid,                &! local block index
      kp1,                &! k+1
      kvdc,               &! k index into mixing coefficient array
      mt2                  ! index for when temperature VDC is different

   real (r8), dimension(nx_block,ny_block) :: &
      VTFB  ! vertical tracer flux across bottom of T-box

!-----------------------------------------------------------------------
!
!  set up vertical indices
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   call timer_start(timer_vdifft, block_id=bid)

   if (k  <  km) then
      kp1 = k + 1
   else
      kp1 = km
   endif

   kvdc = min(k,size(VDC,DIM=3))  !*** reduce to 1 if VDC 2-d array 

!-----------------------------------------------------------------------
!
!  start loop over tracers 
!  set mt2 to be either 1 (if VDC same for all tracers) or n
!    if coefficients are different for tracers
!
!-----------------------------------------------------------------------

   do n = 1,nt
      mt2 = min(n,size(VDC,DIM=4))

!-----------------------------------------------------------------------
!
!     surface tracer fluxes
!
!-----------------------------------------------------------------------

      if (k == 1) then
         VTF(:,:,n,bid) = merge(STF(:,:,n), c0, KMT(:,:,bid) >= 1)
      endif

!-----------------------------------------------------------------------
!
!     vertical tracer flux vdc*Dz(T) at bottom  of T box
!     calculate Dz(vdc*Dz(T))
!
!-----------------------------------------------------------------------

      if (partial_bottom_cells) then
!CDIR COLLAPSE
         VTFB = merge(VDC(:,:,kvdc,mt2,bid)*                    &
                      (TOLD(:,:,k  ,n) - TOLD(:,:,kp1,n))/      &
                      (p5*(DZT(:,:,k,bid) + DZT(:,:,kp1,bid)))  &
                      ,c0, KMT(:,:,bid) > k)

         if (lbottom_heat_flx .and. n == 1) then
!CDIR COLLAPSE
            VTFB = merge( -bottom_heat_flx, VTFB,       &
                         k == KMT(:,:,bid) .and.        &
                         (zt(k) + p5*DZT(:,:,k,bid)) >= &
                         bottom_heat_flx_depth)
         endif

!CDIR COLLAPSE
         VDTK(:,:,n) = merge((VTF(:,:,n,bid) - VTFB)/DZT(:,:,k,bid), &
                             c0, k <= KMT(:,:,bid))

      else
!CDIR COLLAPSE
         VTFB = merge(VDC(:,:,kvdc,mt2,bid)*                      &
                      (TOLD(:,:,k  ,n) - TOLD(:,:,kp1,n))*dzwr(k) &
                      ,c0, KMT(:,:,bid) > k)

         if (lbottom_heat_flx .and. n == 1) then
!CDIR COLLAPSE
            VTFB = merge( -bottom_heat_flx, VTFB,      &
                         k == KMT(:,:,bid) .and.       &
                         zw(k) >= bottom_heat_flx_depth)
         endif

!CDIR COLLAPSE
         VDTK(:,:,n) = merge((VTF(:,:,n,bid) - VTFB)*dzr(k), &
                             c0, k <= KMT(:,:,bid))
      endif

!-----------------------------------------------------------------------
!
!     set top value of VTF to bottom value for next pass at level k+1
!
!-----------------------------------------------------------------------

!CDIR COLLAPSE
      VTF(:,:,n,bid) = VTFB  

   enddo
      
   call timer_stop(timer_vdifft, block_id=bid)

!-----------------------------------------------------------------------
!EOC

 end subroutine vdifft

!***********************************************************************
!BOP
! !IROUTINE: vdiffu
! !INTERFACE:

 subroutine vdiffu(k, VDUK, VDVK, UOLD, VOLD, SMF, this_block)

! !DESCRIPTION:
!  Computes vertical diffusion of momentum:
!  \begin{equation}
!   D_V(u,v) = \delta_z(\mu\delta_z(u,v))
!  \end{equation}
!  \begin{eqnarray}
!   \mu\delta_z(u,v) = (\tau_x,\tau_y)    &{\rm at }& z=0 \\
!   \mu\delta_z(u,v) = c|{\bf\rm u}|(u,v) &{\rm at }& z=-H_U
!  \end{eqnarray}
!
!  This routine must be called successively with $k = 1,2,3,$...
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! vertical level index

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UOLD,                 &! U velocity at old time level
      VOLD                   ! V velocity at old time level

   real (r8), dimension(nx_block,ny_block,2), intent(in) :: &
      SMF                    ! surface momentum fluxes

   type (block), intent(in) :: &
      this_block             ! block info for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
       VDUK,             &! returns Vdiff(U)
       VDVK               ! returns Vdiff(V)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid,               &! local block index
      i,j,               &! loop indices
      kp1,               &! k+1
      kvvc                ! index into viscosity array

   real (r8) :: & 
      vmag                ! temp for bottom drag

   real (r8), dimension(nx_block,ny_block) :: & 
      VUFB,VVFB,        &! vertical momentum fluxes at bottom of box
      WORK               ! local work array

!-----------------------------------------------------------------------
!
!  set vertical level indices
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   call timer_start(timer_vdiffu, block_id=bid)

   if (k  <  km) then
      kp1 = k + 1
   else
      kp1 = km
   endif

   kvvc = min(k,size(VVC,DIM=3))  !*** reduce to 1 if VVC 2-d array 

!-----------------------------------------------------------------------
!
!  surface momentum fluxes (wind stress)
!
!-----------------------------------------------------------------------

   if (k == 1) then
      VUF(:,:,bid) = merge(SMF(:,:,1), c0, KMU(:,:,bid) >= 1)
      VVF(:,:,bid) = merge(SMF(:,:,2), c0, KMU(:,:,bid) >= 1)
   endif

!-----------------------------------------------------------------------
!
!  vertical momentum flux vvc*Dz{(Ub,Vb)} at bottom of U box
!
!-----------------------------------------------------------------------

   if (partial_bottom_cells) then

      if (k < km) then
         WORK = p5*(DZU(:,:,k,bid)+DZU(:,:,kp1,bid))
      else
         WORK = p5*DZU(:,:,kp1,bid)
      endif

!CDIR COLLAPSE
      VUFB = VVC(:,:,kvvc,bid)*(UOLD(:,:,k  ) -  &
                                UOLD(:,:,kp1))/WORK
!CDIR COLLAPSE
      VVFB = VVC(:,:,kvvc,bid)*(VOLD(:,:,k  ) -  &
                                VOLD(:,:,kp1))/WORK
   else
!CDIR COLLAPSE
      VUFB = VVC(:,:,kvvc,bid)*(UOLD(:,:,k  ) -  &
                                UOLD(:,:,kp1))*dzwr(k)
!CDIR COLLAPSE
      VVFB = VVC(:,:,kvvc,bid)*(VOLD(:,:,k  ) -  &
                                VOLD(:,:,kp1))*dzwr(k)
   endif

!-----------------------------------------------------------------------
!
!  bottom momentum fluxes from quadratic drag law
!
!-----------------------------------------------------------------------

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      if (k == KMU(i,j,bid)) then 
         vmag = bottom_drag*sqrt(UOLD(i,j,k)**2 + VOLD(i,j,k)**2)
         VUFB(i,j) = vmag*UOLD(i,j,k)
         VVFB(i,j) = vmag*VOLD(i,j,k)
      endif
   end do
   end do

!-----------------------------------------------------------------------
!
!  calculate Dz(vdc*Dz{(U,V)})
!
!-----------------------------------------------------------------------

   if (partial_bottom_cells) then
      VDUK = merge((VUF(:,:,bid) - VUFB)/DZU(:,:,k,bid), &
                   c0, k <= KMU(:,:,bid))
      VDVK = merge((VVF(:,:,bid) - VVFB)/DZU(:,:,k,bid), &
                   c0, k <= KMU(:,:,bid))
   else
      VDUK = merge((VUF(:,:,bid) - VUFB)*dzr(k), &
                   c0, k <= KMU(:,:,bid))
      VDVK = merge((VVF(:,:,bid) - VVFB)*dzr(k), &
                   c0, k <= KMU(:,:,bid))
   endif

!-----------------------------------------------------------------------
!
!  set top value of (VUF,VVF) to bottom value for next pass (level k+1)
!
!-----------------------------------------------------------------------

   VUF(:,:,bid) = VUFB  
   VVF(:,:,bid) = VVFB  

   call timer_stop(timer_vdiffu, block_id=bid)

!-----------------------------------------------------------------------
!
!  accumulate time-average of vertical diffusion fluxes if reqd
!
!-----------------------------------------------------------------------

   call accumulate_tavg_field(VUF(:,:,bid),tavg_VUF,bid,k)
   call accumulate_tavg_field(VVF(:,:,bid),tavg_VVF,bid,k)

!-----------------------------------------------------------------------
!EOC

 end subroutine vdiffu

!***********************************************************************
!BOP
! !IROUTINE: impvmixt
! !INTERFACE:

 subroutine impvmixt(TNEW, TOLD, PSFC, nfirst, nlast, this_block)

! !DESCRIPTION:
!  Computes the implicit vertical mixing of tracers
!  using tridiagonal solver in each vertical column.
!  Since the top boundary condition (given tracer flux) does not depend
!  on the tracer value, its influence is accounted for in the explicit
!  part (subroutine vdifft).  The bottom b.c. is zero flux.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
      TNEW         ! on input, contains right hand side
                   ! on output, contains updated tracers at new time

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TOLD         ! old tracer to update with del(Tracer)

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      PSFC         ! surface pressure for use in determining
                   ! variable thickness surface layer

   integer (int_kind), intent(in) :: & 
      nfirst, nlast        ! only update tracers from nfirst to nlast

   type (block), intent(in) :: &
      this_block             ! block info for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,k,n,           &! dummy loop indices
      ib,ie,jb,je,       &! beg,end indices for physical domain
      mt2,               &! index for selecting tracer coefficient
      bid                 ! local block address

   real (r8), dimension(nx_block,ny_block) :: &
      A,B,C,D             ! various temporaries

   real (r8), dimension(nx_block,ny_block,km) :: & 
      E,F                 ! various work arrays

   real (r8), dimension(km) :: & 
      hfac_t

   real (r8), dimension(nx_block,ny_block) :: &
      H1            ! factor containing full thickness of sfc layer

!-----------------------------------------------------------------------
!
!  check and initialize some necessary quantities
!
!-----------------------------------------------------------------------

   if (nfirst > nlast .or. nfirst > nt) return

   bid = this_block%local_id
   ib  = this_block%ib
   ie  = this_block%ie
   jb  = this_block%jb
   je  = this_block%je

   call timer_start(timer_impvmixt,block_id=bid)

   !*** note that this array is overwritten for partial bottom cells
   do k=1,km
      hfac_t(k) = dz(k)/c2dtt(k)
   end do

   !*** correct hfac for full surface layer thickness

   if (sfc_layer_type == sfc_layer_varthick) then
      H1 = hfac_t(1) + PSFC/(grav*c2dtt(1))
   else
      H1 = hfac_t(1)
   endif

!-----------------------------------------------------------------------
!
!  loop over tracers
!
!-----------------------------------------------------------------------

   do n = nfirst,nlast

      mt2 = min(n,size(VDC,DIM=4))

!CDIR COLLAPSE
      do j=jb,je
      do i=ib,ie

         !*** perform tridiagonal solve in the vertical for every
         !*** horizontal grid point

         A(i,j) = afac_t(1)*VDC(i,j,1,mt2,bid)
         D(i,j) = H1(i,j) + A(i,j)
         E(i,j,1) = A(i,j)/D(i,j)
         B(i,j) = H1(i,j)*E(i,j,1)
         F(i,j,1) = hfac_t(1)*TNEW(i,j,1,n)/D(i,j)
      end do
      end do

      do k=2,km

       if (partial_bottom_cells) then
        do j=jb,je
        do i=ib,ie
           C(i,j) = A(i,j)
           A(i,j) = aidif*VDC(i,j,k,mt2,bid)/ &
                    (p5*(DZT(i,j,k  ,bid) + &
                         DZT(i,j,k+1,bid)))
           hfac_t(k) = DZT(i,j,k,bid)/c2dtt(k)

           ! Note: this code is duplicated below for the case where
           !       partial_bottom_cells is .false.
           if (k > KMT(i,j,bid)) then
              F(i,j,k) = c0
           else
              if (k == KMT(i,j,bid)) then
                 D(i,j) = hfac_t(k)+B(i,j)
              else
                 D(i,j) = hfac_t(k)+A(i,j)+B(i,j)
              endif

              E(i,j,k) = A(i,j)/D(i,j)
              B(i,j) = (hfac_t(k) + B(i,j))*E(i,j,k)
   
              F(i,j,k) = (hfac_t(k)*TNEW(i,j,k,n) + &
                          C(i,j)*F(i,j,k-1))/D(i,j)
           endif

        end do
        end do
       else ! no partial_bottom_cells
        do j=jb,je
        do i=ib,ie
           C(i,j) = A(i,j)
           A(i,j) = afac_t(k)*VDC(i,j,k,mt2,bid)

           ! Note: this code is duplicated above for the case where
           !       partial_bottom_cells is .true.
           if (k > KMT(i,j,bid)) then
              F(i,j,k) = c0
           else
              if (k == KMT(i,j,bid)) then
                 D(i,j) = hfac_t(k)+B(i,j)
              else
                 D(i,j) = hfac_t(k)+A(i,j)+B(i,j)
              endif

              E(i,j,k) = A(i,j)/D(i,j)
              B(i,j) = (hfac_t(k) + B(i,j))*E(i,j,k)
   
              F(i,j,k) = (hfac_t(k)*TNEW(i,j,k,n) + &
                          C(i,j)*F(i,j,k-1))/D(i,j)
           endif

        end do
        end do
       endif ! partial_bottom_cells
      end do  ! k


      do k=km-1,1,-1
      do j=jb,je
      do i=ib,ie

         !*** back substitution

         if (k < KMT(i,j,bid)) &
            F(i,j,k) = F(i,j,k) + E(i,j,k)*F(i,j,k+1)

      end do
      end do
      end do

!-----------------------------------------------------------------------
!
!     The solution of the system (UZ) is DeltaTracer -- it has 
!     already been multiplied by dtt so we can advance the tracers 
!     directly
!
!-----------------------------------------------------------------------

!CDIR COLLAPSE
      do k=1,km
      do j=jb,je
      do i=ib,ie

         TNEW(i,j,k,n) = TOLD(i,j,k,n) + F(i,j,k)

      end do
      end do
      end do

!-----------------------------------------------------------------------
!
!  end tracer loop and stop timer
!
!-----------------------------------------------------------------------

   end do
   call timer_stop(timer_impvmixt,block_id=bid)

!-----------------------------------------------------------------------
!EOC

 end subroutine impvmixt

!***********************************************************************

 subroutine impvmixt_tavg(TNEW, bid)

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TNEW         ! on input, contains tracer to update from
                   ! on output, contains updated tracers at new time

   integer (int_kind), intent(in) ::  &
      bid                  ! local block address

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block) :: & 
      WORK1, WORK2

   integer (int_kind) ::  &
      k,n,                &! dummy loop indices
      mt2                  ! index for selecting tracer coefficient

!-----------------------------------------------------------------------

   do n = 1,nt
      mt2 = min(n,size(VDC,DIM=4))
      if (accumulate_tavg_now(tavg_DIA_IMPVF_TRACER(n))) then
         do k=1,km-1
            if (allocated(VDC_GM)) then
               WORK1 = VDC(:,:,k,mt2,bid) - VDC_GM(:,:,k,bid)
            else
               WORK1 = VDC(:,:,k,mt2,bid)
            endif
            if (partial_bottom_cells) then
               WORK2 = merge(WORK1*(TNEW(:,:,k,n) - TNEW(:,:,k+1,n))/        &
                             (p5*(DZT(:,:,k,bid) + DZT(:,:,k+1,bid)))        &
                             ,c0, k < KMT(:,:,bid))
               if (lbottom_heat_flx .and. n == 1) then
                  WORK2 = merge( -bottom_heat_flx, WORK2,      &
                                k == KMT(:,:,bid) .and.        &
                                (zt(k) + p5*DZT(:,:,k,bid)) >= &
                                bottom_heat_flx_depth)
               endif
            else
               WORK2 = merge(WORK1*(TNEW(:,:,k,n) - TNEW(:,:,k+1,n))*dzwr(k) &
                             ,c0, k < KMT(:,:,bid))
               if (lbottom_heat_flx .and. n == 1) then
                  WORK2 = merge( -bottom_heat_flx, WORK2,      &
                                k == KMT(:,:,bid) .and.        &
                                zw(k) >= bottom_heat_flx_depth)
               endif
            endif
            call accumulate_tavg_field(WORK2,tavg_DIA_IMPVF_TRACER(n),bid,k)
         end do
      endif
   end do

 end subroutine impvmixt_tavg

!***********************************************************************
!BOP
! !IROUTINE: impvmixt_correct
! !INTERFACE:

 subroutine impvmixt_correct(TNEW, PSFC, RHS, nfirst, nlast, this_block)

! !DESCRIPTION:
!
!  Computes implicit vertical mixing of tracers for a corrector step
!  using tridiagonal solver in each vertical column.
!  Since the top boundary condition (given tracer flux) does not depend
!  on the tracer value, its influence is accounted for in the explicit
!  part (subroutine vdifft).  The bottom b.c. is zero flux.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
      TNEW         ! on input, contains tracer to update from
                   ! on output, contains updated tracers at new time

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      PSFC         ! surface pressure for use in determining
                   ! variable thickness surface layer

   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      RHS          ! right hand side of linear system

   integer (int_kind), intent(in) :: & 
      nfirst, nlast        ! only update tracers from nfirst to nlast

   type (block), intent(in) :: &
      this_block             ! block info for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      i,j,k,n,            &! dummy loop indices
      ib,ie,jb,je,        &! beg,end indices for physical domain
      mt2,                &! index for selecting tracer coefficient
      bid                  ! local block address

   real (r8), dimension(nx_block,ny_block) ::           &
      A,B,C,D              ! various temporaries

   real (r8), dimension(nx_block,ny_block,km) :: &
      E,F                  ! various work arrays

   real (r8), dimension(km) :: &
      hfac_t

   real (r8), dimension(nx_block,ny_block) :: &
      H1            ! factor containing full thickness of sfc layer

!-----------------------------------------------------------------------
!
!  check and initialize some necessary quantities
!
!-----------------------------------------------------------------------

   if (nfirst > nlast .or. nfirst > nt) return

   bid = this_block%local_id
   ib  = this_block%ib
   ie  = this_block%ie
   jb  = this_block%jb
   je  = this_block%je

   call timer_start(timer_impvmixt,block_id=bid)

   !*** note that this array is overwritten for partial bottom cells
   do k=1,km
      hfac_t(k) = dz(k)/c2dtt(k)
   end do

   !*** correct hfac for full surface layer thickness

   if (sfc_layer_type == sfc_layer_varthick) then
      H1 = hfac_t(1) + PSFC/(grav*c2dtt(1))
   else
      H1 = hfac_t(1)
   endif

!-----------------------------------------------------------------------
!
!  loop over tracers
!
!-----------------------------------------------------------------------

   do n = nfirst,nlast

      mt2 = min(n,size(VDC,DIM=4))

      !*** perform tridiagonal solve in the vertical for every
      !*** horizontal grid point

!CDIR COLLAPSE
      do j=jb,je
      do i=ib,ie

         A(i,j) = afac_t(1)*VDC(i,j,1,mt2,bid)
         D(i,j) = H1(i,j) + A(i,j)
         E(i,j,1) = A(i,j)/D(i,j)
         B(i,j) = H1(i,j)*E(i,j,1)
         F(i,j,1) = hfac_t(1)*RHS(i,j,n)/D(i,j)

      end do
      end do

      do k=2,km

        if (partial_bottom_cells) then
          do j=jb,je
          do i=ib,ie
             C(i,j) = A(i,j)
             A(i,j) = aidif*VDC(i,j,k,mt2,bid)/ &
                      (p5*(DZT(i,j,k  ,bid) + DZT(i,j,k+1,bid)))
             hfac_t(k) = DZT(i,j,k,bid)/c2dtt(k)

             ! Note: this code is duplicated below for the case where
             !       partial_bottom_cells is .false.
             if (k > KMT(i,j,bid)) then
                F(i,j,k) = c0
             else
                if (k == KMT(i,j,bid)) then
                   D(i,j) = hfac_t(k)+B(i,j)
                else
                   D(i,j) = hfac_t(k)+A(i,j)+B(i,j)
                endif

                E(i,j,k) = A(i,j)/D(i,j)
                B(i,j) = (hfac_t(k) + B(i,j))*E(i,j,k)

                F(i,j,k) = C(i,j)*F(i,j,k-1)/D(i,j)
             endif

          end do
          end do
        else ! no partial_bottom_cells
          do j=jb,je
          do i=ib,ie
             C(i,j) = A(i,j)
             A(i,j) = afac_t(k)*VDC(i,j,k,mt2,bid)

             ! Note: this code is duplicated above for the case where
             !       partial_bottom_cells is .true.
             if (k > KMT(i,j,bid)) then
                F(i,j,k) = c0
             else
                if (k == KMT(i,j,bid)) then
                   D(i,j) = hfac_t(k)+B(i,j)
                else
                   D(i,j) = hfac_t(k)+A(i,j)+B(i,j)
                endif

                E(i,j,k) = A(i,j)/D(i,j)
                B(i,j) = (hfac_t(k) + B(i,j))*E(i,j,k)

                F(i,j,k) = C(i,j)*F(i,j,k-1)/D(i,j)
             endif

          end do
          end do
        endif ! partial_bottom_cells
      end do ! k

      !*** back substitution

      do k=km-1,1,-1
      do j=jb,je
      do i=ib,ie
         if (k < KMT(i,j,bid)) &
            F(i,j,k) = F(i,j,k) + E(i,j,k)*F(i,j,k+1)
      end do
      end do
      end do

!-----------------------------------------------------------------------
!
!     The solution of the system (UZ) is DeltaTracer -- it has 
!     already been multiplied by dtt so we can advance the tracers 
!     directly
!
!-----------------------------------------------------------------------

!CDIR COLLAPSE
      do k=1,km
      do j=jb,je
      do i=ib,ie
         TNEW(i,j,k,n) = TNEW(i,j,k,n) + F(i,j,k)
      end do
      end do
      end do

!-----------------------------------------------------------------------
!
!  end tracer loop and stop timer
!
!-----------------------------------------------------------------------

   end do
   call timer_stop(timer_impvmixt,block_id=bid)

!-----------------------------------------------------------------------
!EOC

 end subroutine impvmixt_correct

!***********************************************************************
!BOP
! !IROUTINE: impvmixu
! !INTERFACE:

 subroutine impvmixu(UNEW, VNEW, this_block)

! !DESCRIPTION:
!  Computes the implicit vertical mixing of momentum
!  using a tridiagonal solver.
!  Since the top boundary condition (given wind stress) does not 
!  depend on the velocity, its influence is accounted for in the 
!  explicit part (subroutine vdiffu).  On the other hand, the bottom 
!  b.c. does depend on velocity ($\tau = c|{\bf\rm u}|(u,v)$).,
!  but both the speed and velocity in this term are evaluated
!  at the old time, so the bottom drag is also accounted for
!  in the explicit part.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km), intent(inout) :: &
      UNEW, VNEW             ! U,V velocities at new time level

! !INPUT PARAMETERS:

   type (block), intent(in) :: &
      this_block             ! block info for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,k,             &! dummy loop indices
      ib,ie,jb,je,       &! beg,end indices for physical domain
      bid                 ! local block address

   real (r8), dimension(nx_block,ny_block) :: &
      A,B,C,D             ! various temp variables

   real (r8), dimension(nx_block,ny_block,km) :: & 
      E,F1,F2             ! various work arrays

   real (r8), dimension(km) :: & 
      hfac_u

!-----------------------------------------------------------------------
!
!  initialize
!
!-----------------------------------------------------------------------

   bid = this_block%local_id
   ib  = this_block%ib
   ie  = this_block%ie
   jb  = this_block%jb
   je  = this_block%je

   A=c0
   B=c0
   C=c0
   D=c0
   E=c0
   F1=c0
   F2=c0

   call timer_start(timer_impvmixu,block_id=bid)

   !*** note that this array is over-written for partial bottom cells
   do k=1,km
      hfac_u(k) = dz(k)/c2dtu
   end do


!-----------------------------------------------------------------------
!
!  solve tridiagonal system at each grid point
!
!-----------------------------------------------------------------------
 
!CDIR COLLAPSE
   do j=jb,je
   do i=ib,ie

      A(i,j)    = afac_u(1)*VVC(i,j,1,bid)
      D(i,j)    = hfac_u(1) + A(i,j)
      E(i,j,1)  = A(i,j)/D(i,j)
      B(i,j)    = hfac_u(1)*E(i,j,1)
      F1(i,j,1) = hfac_u(1)*UNEW(i,j,1)/D(i,j)
      F2(i,j,1) = hfac_u(1)*VNEW(i,j,1)/D(i,j)

   end do
   end do

   do k=2,km

      if (partial_bottom_cells) then
!CDIR COLLAPSE
        do j=jb,je
        do i=ib,ie
           C(i,j) = A(i,j)
           hfac_u(k) = DZU(i,j,k,bid)/c2dtu
           A(i,j) = aidif*VVC(i,j,k,bid)/(p5*(DZU(i,j,k,bid) + &
                                              DZU(i,j,k+1,bid)))

           ! Note: this code is duplicated below for the case where
           !       partial_bottom_cells is .false.
           if (k < KMU(i,j,bid)) then

              D(i,j)    = hfac_u(k) + A(i,j) + B(i,j)
              E(i,j,k)  = A(i,j)/D(i,j)
              B(i,j)    = (hfac_u(k) + B(i,j))*E(i,j,k)
              F1(i,j,k) = (hfac_u(k)*UNEW(i,j,k) + &
                           C(i,j)*F1(i,j,k-1))/D(i,j)
              F2(i,j,k) = (hfac_u(k)*VNEW(i,j,k) + &
                           C(i,j)*F2(i,j,k-1))/D(i,j)

           else if (k == KMU(i,j,bid)) then

              D(i,j)    = hfac_u(k) + B(i,j)
              E(i,j,k)  = A(i,j)/D(i,j)
              B(i,j)    = (hfac_u(k) + B(i,j))*E(i,j,k)
              F1(i,j,k) = (hfac_u(k)*UNEW(i,j,k) + &
                           C(i,j)*F1(i,j,k-1))/D(i,j)
              F2(i,j,k) = (hfac_u(k)*VNEW(i,j,k) + &
                           C(i,j)*F2(i,j,k-1))/D(i,j)
           else
              F1(i,j,k) = c0
              F2(i,j,k) = c0
         endif
        end do
        end do
      else ! no partial_bottom_cell
        do j=jb,je
        do i=ib,ie
           C(i,j) = A(i,j)
           A(i,j) = afac_u(k)*VVC(i,j,k,bid)

           ! Note: this code is duplicated above for the case where
           !       partial_bottom_cells is .true.
           if (k < KMU(i,j,bid)) then

              D(i,j)    = hfac_u(k) + A(i,j) + B(i,j)
              E(i,j,k)  = A(i,j)/D(i,j)
              B(i,j)    = (hfac_u(k) + B(i,j))*E(i,j,k)
              F1(i,j,k) = (hfac_u(k)*UNEW(i,j,k) + &
                           C(i,j)*F1(i,j,k-1))/D(i,j)
              F2(i,j,k) = (hfac_u(k)*VNEW(i,j,k) + &
                           C(i,j)*F2(i,j,k-1))/D(i,j)

           else if (k == KMU(i,j,bid)) then

              D(i,j)    = hfac_u(k) + B(i,j)
              E(i,j,k)  = A(i,j)/D(i,j)
              B(i,j)    = (hfac_u(k) + B(i,j))*E(i,j,k)
              F1(i,j,k) = (hfac_u(k)*UNEW(i,j,k) + &
                           C(i,j)*F1(i,j,k-1))/D(i,j)
              F2(i,j,k) = (hfac_u(k)*VNEW(i,j,k) + &
                           C(i,j)*F2(i,j,k-1))/D(i,j)
           else
              F1(i,j,k) = c0
              F2(i,j,k) = c0
         endif
        end do
        end do
      endif ! partial_bottom_cells
   end do ! k

   do k=km-1,1,-1
   do j=jb,je
   do i=ib,ie
      if (k < KMU(i,j,bid)) then
         F1(i,j,k) = F1(i,j,k) + E(i,j,k)*F1(i,j,k+1)
         F2(i,j,k) = F2(i,j,k) + E(i,j,k)*F2(i,j,k+1)
      endif
   end do
   end do
   end do

!CDIR COLLAPSE
   do k=1,km
   do j=jb,je
   do i=ib,ie
      UNEW(i,j,k) = F1(i,j,k)
      VNEW(i,j,k) = F2(i,j,k)
   end do
   end do
   end do

!-----------------------------------------------------------------------
!
!  UVEL,VVEL(newtime) now hold a modified rhs that has already been
!  multiplied by dtu to advance UVEL.
!
!-----------------------------------------------------------------------
   call timer_stop(timer_impvmixu,block_id=bid)

!-----------------------------------------------------------------------
!EOC

 end subroutine impvmixu

!***********************************************************************
!BOP
! !IROUTINE: convad
! !INTERFACE:

 subroutine convad(TNEW, RHONEW, this_block, bid)

! !DESCRIPTION:
!  Convectively adjusts tracers at adjacent depth
!  levels if they are gravitationally unstable.
!  This routine convectively adjusts tracers at the new
!  time (before they are updated at the end of step).
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
      TNEW                   ! tracers at new time level

   real (r8), dimension(nx_block,ny_block,km), intent(inout) :: &
      RHONEW                 ! density at new time level

! !INPUT PARAMETERS:

   type (block), intent(in) :: &
      this_block             ! block info for current block

   integer (int_kind), intent(in) :: &
      bid                    ! local block index

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nc,                &! dummy index for number of adjustment passes
      n,                 &! dummy index for tracer number
      k, ks               ! dummy vertical level indices

   real (r8), dimension(nx_block,ny_block) :: &
      RHOK,       &! density of water after adiab. displacement to k+1
      RHOKP,      &! actual density of water at level k+1
      WORK1,WORK2  ! local work space

   logical (log_kind) :: &
      lpec, lncnv         ! logical flags for computing tavg diagnostics

!-----------------------------------------------------------------------
!
!  do nothing if convective adjustment not chosen
!
!-----------------------------------------------------------------------

   if (convection_itype /= convect_type_adjust) return

!-----------------------------------------------------------------------
!
!  perform convective adjustment
!
!-----------------------------------------------------------------------

   do nc = 1,nconvad  ! repeat nconvad times

      do ks = 1,2
         do k = ks,km-1,2

            call state(k  ,k+1,TNEW(:,:,k  ,1), TNEW(:,:,k  ,2), &
                               this_block, RHOOUT=RHOK )
            call state(k+1,k+1,TNEW(:,:,k+1,1), TNEW(:,:,k+1,2), &
                               this_block, RHOOUT=RHOKP)
              
            if (partial_bottom_cells) then
               do n = 1,nt
                  where ((RHOK > RHOKP) .and. (k < KMT(:,:,bid)))
                     TNEW(:,:,k,n) = c1/                          &
                        (DZT(:,:,k  ,bid)/dttxcel(k  ) +             &
                         DZT(:,:,k+1,bid)/dttxcel(k+1))              &
                       *(DZT(:,:,k  ,bid)/dttxcel(k  )*TNEW(:,:,k,n) &
                       + DZT(:,:,k+1,bid)/dttxcel(k+1)*TNEW(:,:,k+1,n))
                     TNEW(:,:,k+1,n) = TNEW(:,:,k,n)
                  end where
               enddo
            else
               do n = 1,nt
                  where ((RHOK > RHOKP) .and. (k < KMT(:,:,bid)))
                     TNEW(:,:,k,n) = dzwxcel(k)*                   &
                                    (dztxcel(k  )*TNEW(:,:,k  ,n)  &
                                   + dztxcel(k+1)*TNEW(:,:,k+1,n))
                     TNEW(:,:,k+1,n) = TNEW(:,:,k,n)
                  endwhere
               enddo
            endif

         enddo
      enddo

   enddo

!-----------------------------------------------------------------------
!
!  compute new density based on new tracers and compute adjustment
!  diagnostic
!
!-----------------------------------------------------------------------

   lpec  = accumulate_tavg_now(tavg_PEC)
   lncnv = accumulate_tavg_now(tavg_NCNV)

   do k = 1,km
      if (lpec .or. lncnv) then
         WORK1 = RHONEW(:,:,k)
      endif

      call state(k,k,TNEW(:,:,k,1), TNEW(:,:,k,2),  &
                     this_block, RHOOUT=RHONEW(:,:,k))

      if (lpec) then
         WORK2 = RHONEW(:,:,k) - WORK1
         call accumulate_tavg_field(WORK2,tavg_PEC,bid,k)
      endif

      if (lncnv) then
         if (.not. lpec) WORK2 = RHONEW(:,:,k) - WORK1

         where (abs(WORK2) > 1.e-8_r8)
            WORK1 = c1
         elsewhere
            WORK1 = c0
         end where

         call accumulate_tavg_field(WORK1,tavg_NCNV,bid,k)
      endif

   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine convad

!***********************************************************************

 end module vertical_mix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
