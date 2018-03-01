!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      module hmix_gm

!BOP
! !MODULE: hmix_gm

! !DESCRIPTION:
!  This module contains routines for computing horizontal mixing
!  using the Gent-McWilliams eddy transport parameterization
!  and isopycnal diffusion.

! !REVISION HISTORY:
!  SVN:$Id: hmix_gm.F90 48834 2013-07-09 19:37:42Z mlevy@ucar.edu $

! !USES:

      use kinds_mod
      use blocks
      use distribution
      use domain
      use constants
      use broadcast
      use grid
      use io
      use vertical_mix
      use vmix_kpp
      use state_mod
      use time_management
      use tavg
      use diagnostics
      use exit_mod
      use registry
      use hmix_gm_submeso_share

#ifdef CCSMCOUPLED
   use shr_sys_mod
#endif
      use timers, only: timer_start, timer_stop, get_timer


      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: init_gm,   &
                hdifft_gm

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     variables to save from one call to next
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ktp = 1, kbt = 2     ! refer to the top and bottom halves of a
                              !  grid cell, respectively

      real (r8), dimension(:), allocatable :: &
         kappa_depth          ! depth dependence for KAPPA 

      real (r8), dimension(:,:,:), allocatable :: &
         HYXW, HXYS, &        ! west and south-shifted values of above
         RB,         &        ! Rossby radius
         RBR,        &        ! inverse of Rossby radius
         BTP,        &        ! beta plane approximation
         BL_DEPTH,   &        ! boundary layer depth
         UIT, VIT,   &        ! work arrays for isopycnal mixing velocities
         WTOP_ISOP, WBOT_ISOP ! vertical component of isopycnal velocities

      real (r8), dimension(:,:,:,:,:,:), allocatable :: &
         SF_SLX, SF_SLY       ! components of the merged streamfunction

      real (r8), dimension(:,:,:,:,:), allocatable :: &
         SLA_SAVE             ! isopycnal slopes

      real (r8), dimension(:,:,:,:), allocatable :: &
         FZTOP                ! vertical flux

      logical (log_kind), dimension(:), allocatable :: &
         compute_kappa        ! compute spatially varying coefficients
                              !  this time step?

      logical (log_kind) ::   &
         diff_tapering,       &   ! different tapering for two diffusivities
         cancellation_occurs, &   ! specified choices for the isopycnal and
                                  !  thickness diffusion coefficients result in 
                                  !  cancellation of some tensor elements
         read_n2_data             ! if .true., use climatological N^2 data 
                                  !  instead of model dependent N^2

!-----------------------------------------------------------------------
!
!     KAPPA_LATERAL and KAPPA_VERTICAL are 2D and 3D arrays, respectively,
!     containing the spatial variations of the isopycnal (KAPPA_ISOP)
!     and thickness (KAPPA_THIC) diffusion coefficients. Except in kappa_type_eg,
!     KAPPA_LATERAL has the actual diffusion coefficient values in cm^2/s,
!     whereas KAPPA_VERTICAL is unitless. So, the total coefficients are
!     constructed using
!
!      KAPPA_ISOP or KAPPA_THIC (:,:,:,k,bid) ~ KAPPA_LATERAL(:,:,bid)
!                                             * KAPPA_VERTICAL(:,:,k,bid)
!
!     When kappa_type_eg, KAPPA_VERTICAL contains the isopycnal diffusivity
!     coefficients in cm^2/s and KAPPA_LATERAL is not used!
!
!-----------------------------------------------------------------------

      real (r8), dimension(:,:,:,:,:), allocatable :: &
         KAPPA_ISOP, &      ! 3D isopycnal diffusion coefficient
                            !  for top and bottom half of a grid cell
         KAPPA_THIC, &      ! 3D thickness diffusion coefficient
                            !  for top and bottom half of a grid cell
         HOR_DIFF           ! 3D horizontal diffusion coefficient
                            !  for top and bottom half of a grid cell

      real (r8), dimension(:,:,:), allocatable :: &
         KAPPA_LATERAL      ! horizontal variation of KAPPA in cm^2/s

      real (r8), dimension(:,:,:,:), allocatable :: &
         KAPPA_VERTICAL     ! vertical variation of KAPPA (unitless),
                            !  e.g. normalized buoyancy frequency dependent 
                            !  profiles at the tracer grid points
                            !  ( = N^2 / N_ref^2 ) OR a time-independent
                            !  user-specified function

      real (r8), dimension(:,:,:,:), allocatable :: &
         BUOY_FREQ_SQ,    & ! N^2 defined at level interfaces
         SIGMA_TOPO_MASK    ! bottom topography mask used with kappa_type_eg

!-----------------------------------------------------------------------
!
!     GM specific options
!
!     kappa_freq = how often spatial variations of the diffusion 
!                  coefficients are computed. Same frequency is 
!                  used for both coefficients.
!     slope_control = tanh function (Danabasoglu and McWilliams 1995) or
!                     DM95 with replacement function to tanh or
!                     slope clipping or
!                     method of Gerdes et al (1991)
!     diag_gm_bolus = .true. for diagnostic bolus velocity computation.
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter ::   &
         kappa_type_const         = 1,   &
         kappa_type_depth         = 2,   &
         kappa_type_vmhs          = 3,   &
         kappa_type_hdgr          = 4,   &
         kappa_type_dradius       = 5,   &
         kappa_type_bfreq         = 6,   &
         kappa_type_bfreq_vmhs    = 7,   &
         kappa_type_bfreq_hdgr    = 8,   &
         kappa_type_bfreq_dradius = 9,   &
         kappa_type_eg            = 10,  &
         slope_control_tanh   = 1,       &
         slope_control_notanh = 2,       &
         slope_control_clip   = 3,       &
         slope_control_Gerd   = 4,       &
         kappa_freq_never           = 1, &
         kappa_freq_every_time_step = 2, &
         kappa_freq_once_a_day      = 3  

      integer (int_kind) :: &
         kappa_isop_type,   &   ! choice of KAPPA_ISOP
         kappa_thic_type,   &   ! choice of KAPPA_THIC
         kappa_freq,        &   ! frequency of KAPPA computations
         slope_control          ! choice for slope control

      logical (log_kind) ::  &
         diag_gm_bolus          ! true for diagnostic bolus velocity computation 

!-----------------------------------------------------------------------
!
!     if use_const_ah_bkg_srfbl = .true., then the specified constant
!     value of ah_bkg_srfbl is used as the "maximum" background horizontal 
!     diffusivity within the surface boundary layer. Otherwise,
!     KAPPA_ISOP is utilized as this "maximum".
!
!-----------------------------------------------------------------------
  
      logical (log_kind) ::      &
         use_const_ah_bkg_srfbl, & ! see above 
         transition_layer_on       ! control for transition layer parameterization
                               
      real (r8) ::      &
         ah,            &       ! isopycnal diffusivity
         ah_bolus,      &       ! thickness (GM bolus) diffusivity
         ah_bkg_bottom, &       ! backgroud horizontal diffusivity at k = KMT
         ah_bkg_srfbl,  &       ! backgroud horizontal diffusivity within the
                                !  surface boundary layer
         slm_r,         &       ! max. slope allowed for redi diffusion
         slm_b                  ! max. slope allowed for bolus transport

!-----------------------------------------------------------------------
!
!     the following set of variables are used in Eden and Greatbatch 
!     (2008) KAPPA formulation. They are in the input namelist.
!
!-----------------------------------------------------------------------

      real (r8) ::      &
         const_eg,      &  ! tuning parameter (unitless)
         gamma_eg,      &  ! (> 0) effective upper limit for inverse eddy 
                           !  time scale (unitless)
         kappa_min_eg,  &  ! minimum KAPPA (cm^2/s)
         kappa_max_eg      ! maximum KAPPA (cm^2/s)

!-----------------------------------------------------------------------
!
!     transition layer type variables
!
!----------------------------------------------------------------------- 

      type tlt_info
        real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
           DIABATIC_DEPTH,  &   ! depth of the diabatic region at the
                                !  surface, i.e. mean mixed or boundary layer
                                !  depth
           THICKNESS,       &   ! transition layer thickness
           INTERIOR_DEPTH       ! depth at which the interior, adiabatic
                                !  region starts, i.e.
                                !   = TLT%DIABATIC_DEPTH + TLT%THICKNESS
        integer (int_kind), &
              dimension(nx_block,ny_block,max_blocks_clinic) :: &
           K_LEVEL,  &          ! k level at or below which the interior,
                                !  adiabatic region starts
           ZTW                  ! designates if the interior region
                                !  starts below depth zt or zw.
                                !  ( = 1 for zt, = 2 for zw )
      end type tlt_info

      type (tlt_info) :: &
         TLT                    ! transition layer thickness related fields 

!-----------------------------------------------------------------------
!
!     tavg ids for tavg diagnostics related to diffusivities and
!     isopycnal velocities. Zonal and meridional refer here to logical 
!     space only.
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         tavg_UISOP,        &   ! zonal      isopycnal velocity
         tavg_VISOP,        &   ! meridional isopycnal velocity
         tavg_WISOP,        &   ! vertical   isopycnal velocity
         tavg_KAPPA_ISOP,   &   ! isopycnal  diffusion coefficient 
         tavg_KAPPA_THIC,   &   ! thickness  diffusion coefficient
         tavg_HOR_DIFF,     &   ! horizontal diffusion coefficient
         tavg_DIA_DEPTH,    &   ! depth of the diabatic region at the surface
         tavg_TLT,          &   ! transition layer thickness
         tavg_INT_DEPTH,    &   ! depth at which the interior region starts
         tavg_ADVT_ISOP,    &   ! vertically-integrated T eddy-induced
                                !  advection tendency
         tavg_ADVS_ISOP,    &   ! vertically-integrated S eddy-induced
                                !  advection tendency
         tavg_VNT_ISOP,     &   ! heat flux tendency in grid-y direction
                                !  due to eddy-induced velocity
         tavg_VNS_ISOP          ! salt flux tendency in grid-y direction
                                !  due to eddy-induced velocity

!-----------------------------------------------------------------------
!
!  timers
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      timer_nloop         ! main n loop

!EOC
!***********************************************************************

      contains

!***********************************************************************
!BOP
! !IROUTINE: init_gm
! !INTERFACE:

      subroutine init_gm

! !DESCRIPTION:
!  Initializes various choices and allocates necessary space.
!  Also computes some time-independent arrays.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   type (block) :: &
      this_block         ! block information for current block

   character (char_len) ::    &
      buoyancy_freq_filename, &! file name for the time-independent
                               !  buoyancy frequency squared
      buoyancy_freq_fmt,      &! format (bin or netcdf) of buoyancy file
      message                  ! string to hold error message

   type (datafile) :: &
      buoyancy_freq_data_file  ! data file descriptor for buoyancy freq data

   type (io_field_desc) :: &
      buoyancy_freq_data_in  ! io field descriptor for input buoyancy freq data

   type (io_dim) :: &
      i_dim, j_dim, &  ! dimension descriptors for horiz dims
      k_dim            ! dimension descriptor  for depth

!-----------------------------------------------------------------------
!
!     input namelist for setting various GM options
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error,         &  ! error flag for namelist
      k,                 &  ! vertical level index
      iblock,            &  ! block index
      i, j                  ! lateral indices

   real (r8) ::          &
      kappa_depth_1,     &  ! parameters for variation of KAPPA 
      kappa_depth_2,     &  !  with depth used with the kappa_type_depth
      kappa_depth_scale     !  option

   character (char_len) :: &
      kappa_choice,        &  ! (generic) choice for KAPPA
      kappa_freq_choice,   &  ! frequency choice for KAPPA computation
      kappa_isop_choice,   &  ! choice for KAPPA_ISOP
      kappa_thic_choice,   &  ! choice for KAPPA_THIC
      slope_control_choice    ! choice for slope control

   namelist /hmix_gm_nml/ kappa_isop_choice,                     &
                          kappa_thic_choice,                     &
                          kappa_freq_choice,                     &
                          slope_control_choice,                  &
                          kappa_depth_1, kappa_depth_2,          &
                          kappa_depth_scale, ah, ah_bolus,       &
                          use_const_ah_bkg_srfbl, ah_bkg_srfbl,  &
                          ah_bkg_bottom,                         &
                          slm_r, slm_b, diag_gm_bolus,           &
                          transition_layer_on, read_n2_data,     &
                          buoyancy_freq_filename,                &
                          buoyancy_freq_fmt,                     &
                          const_eg, gamma_eg,                    &
                          kappa_min_eg, kappa_max_eg

!-----------------------------------------------------------------------
!
!     read input namelist for additional GM options
!
!     DEFAULT SETUP: 
!     KAPPA_ISOP      : constant (= ah)
!     KAPPA_THIC      : constant (= ah_bolus)
!     kappa_freq      : never
!     slope control   : method by DM95 with replacing tanh by polynomial
!                        (= slope_control_notanh)
!
!     variation of kappa with depth used with kappa_type_depth is
!       kappa_depth_1 + kappa_depth_2*exp(-z/kappa_depth_scale)
!
!      with kappa_depth_1 = 1, kappa_depth_2 = 0, and 
!      kappa_depth_scale = 150000 cm, i.e. no depth variation
!
!     ah_bolus        : thickness diffusion coefficient (= 0.8e07 cm^2/s)
!     ah_bkg_bottom   : background horizontal diffusion at k=KMT (= 0)
!     use_const_ah_bkg_srfbl : use ah_bkg_srfbl value
!     ah_bkg_srfbl    : background horizontal diffusion within the
!                        surface boundary layer (= 0.8e07 cm^2/s)
!     slm_r           : max. slope allowed for isopycnal (redi) diffusion (= 0.3)
!     slm_b           : max. slope allowed for thickness (bolus) diffusion (= 0.3)
!     diag_gm_bolus   : .true.
!     transition_layer_on: .false.
!     read_n2_data    : .false.
!     const_eg        : 1.
!     gamma_eg        : 300.
!     kappa_min_eg    : 0.35e07 cm^2/s
!     kappa_max_eg    : 5.00e07 cm^2/s
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     register init_gm
!
!-----------------------------------------------------------------------

   call register_string('init_gm')
 
   kappa_isop_type = kappa_type_const
   kappa_thic_type = kappa_type_const
   kappa_freq    = kappa_freq_never
   slope_control = slope_control_notanh
   kappa_depth_1 = c1
   kappa_depth_2 = c0
   kappa_depth_scale = 150000.0_r8
   ah            = 0.8e7_r8
   ah_bolus      = 0.8e7_r8
   ah_bkg_bottom = c0
   ah_bkg_srfbl  = 0.8e7_r8
   slm_r = 0.3_r8
   slm_b = 0.3_r8
   diag_gm_bolus          = .true.
   use_const_ah_bkg_srfbl = .true.
   transition_layer_on    = .false.
   read_n2_data           = .false.
   buoyancy_freq_filename = 'unknown-buoyancy'
   buoyancy_freq_fmt      = 'nc'
   const_eg               = c1 
   gamma_eg               = 300.0_r8 
   kappa_min_eg           = 0.35e7_r8
   kappa_max_eg           = 5.0e7_r8

   if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     if (nml_error /= 0) then
       nml_error = -1
     else
       nml_error =  1
     endif
     do while (nml_error > 0)
       read(nml_in, nml=hmix_gm_nml,iostat=nml_error)
     end do
     if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR reading hmix_gm_nml')
   endif

   if (my_task == master_task) then

     write(stdout,*) ' '
     write(stdout,*) ' Document Namelist Parameters:'
     write(stdout,*) ' ============================ '
     write(stdout,*) ' '
     write(stdout,  hmix_gm_nml)
     write(stdout,*) ' '

     write(stdout,*) '  Gent-McWilliams options:'
     write(stdout,*) '    kappa_isop choice is ',  &
                     trim(kappa_isop_choice)
     write(stdout,*) '    kappa_thic choice is ',  &
                     trim(kappa_thic_choice)
     write(stdout,*) '    kappa computation frequency choice is ',  &
                     trim(kappa_freq_choice)
     write(stdout,*) '    slope control choice is ',  &
                                    trim(slope_control_choice)
     write(stdout,'(a28,1pe13.6)') ' isopycnal diffusion      = ',  &
                                     ah
     write(stdout,'(a28,1pe13.6)') ' thickness diffusion      = ',  &
                                     ah_bolus
     write(stdout,'(a28,1pe13.6)') ' backgroud diff. (bottom) = ',  &
                                     ah_bkg_bottom
     write(stdout,'(a28,1pe13.6)') ' backgroud diff. (srfbl)  = ',  &
                                     ah_bkg_srfbl
     write(stdout,'(a28,1pe13.6)') ' max. slope for redi      = ',  &
                                     slm_r
     write(stdout,'(a28,1pe13.6)') ' max. slope for bolus     = ',  &
                                     slm_b
     if ( diag_gm_bolus ) then
       write(stdout,'(a47)')  &
           '    diagnostic bolus velocity computation is on'
     endif

     if ( use_const_ah_bkg_srfbl ) then
       write(stdout,1001)
     else
       write(stdout,1002)
     endif
 1001   format(/,' Maximum horizontal background diffusion ', &
                 'coefficient', /,                            &
           '  within the boundary layer is constant at ah_bkg_srfbl.')
 1002   format(/,' Maximum horizontal background diffusion ', &
                 'coefficient', /,                            &
           '  within the boundary layer depends on KAPPA_ISOP.')

     if ( transition_layer_on ) then
       write(stdout,'(a33)') '  transition layer scheme is on. '
     endif

     select case (kappa_isop_choice(1:4))
     case ('cons')
       kappa_isop_type = kappa_type_const
     case ('dept')
       kappa_isop_type = kappa_type_depth
     case ('vmhs')
       kappa_isop_type = kappa_type_vmhs
     case ('hdgr')
       kappa_isop_type = kappa_type_hdgr
     case ('drad')
       kappa_isop_type = kappa_type_dradius 
     case ('bfre')
       kappa_isop_type = kappa_type_bfreq 
     case ('bfvm')
       kappa_isop_type = kappa_type_bfreq_vmhs
     case ('bfhd')
       kappa_isop_type = kappa_type_bfreq_hdgr
     case ('bfdr')
       kappa_isop_type = kappa_type_bfreq_dradius
     case ('edgr')
       kappa_isop_type = kappa_type_eg
     case default
       kappa_isop_type = -1000
     end select

     select case (kappa_thic_choice(1:4))
     case ('cons')
       kappa_thic_type = kappa_type_const
     case ('dept')
       kappa_thic_type = kappa_type_depth
     case ('vmhs')
       kappa_thic_type = kappa_type_vmhs
     case ('hdgr')
       kappa_thic_type = kappa_type_hdgr
     case ('drad')
       kappa_thic_type = kappa_type_dradius
     case ('bfre')
       kappa_thic_type = kappa_type_bfreq
     case ('bfvm')
       kappa_thic_type = kappa_type_bfreq_vmhs
     case ('bfhd')
       kappa_thic_type = kappa_type_bfreq_hdgr
     case ('bfdr')
       kappa_thic_type = kappa_type_bfreq_dradius
     case ('edgr')
       kappa_thic_type = kappa_type_eg
     case default
       kappa_thic_type = -1000
     end select

     select case (kappa_freq_choice(1:3))
     case ('nev')
       kappa_freq = kappa_freq_never
     case ('eve')
       kappa_freq = kappa_freq_every_time_step
     case ('onc')
       kappa_freq = kappa_freq_once_a_day
     case default
       kappa_freq = -1000
     end select

     select case (slope_control_choice(1:4))
     case ('tanh')
       slope_control = slope_control_tanh
     case ('nota')
       slope_control = slope_control_notanh
     case ('clip')
       slope_control = slope_control_clip
     case ('Gerd')
       slope_control = slope_control_Gerd
     case default
       slope_control = -1000
     end select

     if ( kappa_isop_type == kappa_type_bfreq          .or.  &
          kappa_isop_type == kappa_type_bfreq_vmhs     .or.  &
          kappa_isop_type == kappa_type_bfreq_hdgr     .or.  &
          kappa_isop_type == kappa_type_bfreq_dradius  .or.  &
          kappa_thic_type == kappa_type_bfreq          .or.  &
          kappa_thic_type == kappa_type_bfreq_vmhs     .or.  &
          kappa_thic_type == kappa_type_bfreq_hdgr     .or.  &
          kappa_thic_type == kappa_type_bfreq_dradius ) then
       if ( read_n2_data ) then
         write(stdout,'(a32)') ' using climatological N^2 data. '
       else
         write(stdout,'(a34)') ' using model data to compute N^2. '
       endif
      endif

     if ( kappa_isop_type == kappa_type_depth  .or.  &
          kappa_thic_type == kappa_type_depth ) then

       if ( kappa_isop_type == kappa_type_depth )  &
         write(stdout,'(a30)') ' KAPPA_ISOP varies with depth.'

       if ( kappa_thic_type == kappa_type_depth )  &
         write(stdout,'(a30)') ' KAPPA_THIC varies with depth.'

       write(stdout,'(a20,1pe13.6)') '    kappa_depth_1 = ',  &
                                          kappa_depth_1
       write(stdout,'(a20,1pe13.6)') '    kappa_depth_2 = ',  &
                                          kappa_depth_2
       write(stdout,'(a24,1pe13.6)') '    kappa_depth_scale = ',  &
                                          kappa_depth_scale

     endif

     if ( kappa_isop_type == kappa_type_eg  .or.  &
          kappa_thic_type == kappa_type_eg ) then

       write(stdout,'(a16,1pe13.6)') '     const_eg = ',  &
                                           const_eg 
       write(stdout,'(a16,1pe13.6)') '     gamma_eg = ',  &
                                           gamma_eg
       write(stdout,'(a16,1pe13.6)') ' kappa_min_eg = ',  &
                                       kappa_min_eg
       write(stdout,'(a16,1pe13.6)') ' kappa_max_eg = ',  &
                                       kappa_max_eg

     endif

   endif

   call broadcast_scalar(kappa_isop_type,        master_task)
   call broadcast_scalar(kappa_thic_type,        master_task)
   call broadcast_scalar(kappa_freq,             master_task)
   call broadcast_scalar(slope_control,          master_task)
   call broadcast_scalar(kappa_depth_1,          master_task)
   call broadcast_scalar(kappa_depth_2,          master_task)
   call broadcast_scalar(kappa_depth_scale,      master_task)
   call broadcast_scalar(ah,                     master_task)
   call broadcast_scalar(ah_bolus,               master_task)
   call broadcast_scalar(ah_bkg_bottom,          master_task)
   call broadcast_scalar(ah_bkg_srfbl,           master_task)
   call broadcast_scalar(slm_r,                  master_task)
   call broadcast_scalar(slm_b,                  master_task)
   call broadcast_scalar(diag_gm_bolus,          master_task)
   call broadcast_scalar(use_const_ah_bkg_srfbl, master_task)
   call broadcast_scalar(transition_layer_on,    master_task)
   call broadcast_scalar(read_n2_data,           master_task)
   call broadcast_scalar(buoyancy_freq_filename, master_task)
   call broadcast_scalar(buoyancy_freq_fmt,      master_task)
   call broadcast_scalar(const_eg,               master_task)
   call broadcast_scalar(gamma_eg,               master_task)
   call broadcast_scalar(kappa_min_eg,           master_task)
   call broadcast_scalar(kappa_max_eg,           master_task)

!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------

   if ( kappa_isop_type == -1000 ) then
     call exit_POP(sigAbort,  &
                   'unknown type for KAPPA_ISOP in GM setup')
   endif
   if ( kappa_thic_type == -1000 ) then
    call exit_POP(sigAbort,  &
                  'unknown type for KAPPA_THIC in GM setup')
   endif
   if ( kappa_freq == -1000 ) then
     call exit_POP(sigAbort,                                       &
                  'unknown type for KAPPA computation frequency in GM setup')
   endif
   if ( slope_control == -1000 ) then
     call exit_POP(sigAbort,  &
                   'unknown slope control method in GM setup')
   endif

   if ( kappa_isop_type == kappa_type_depth  .and.      &
        ( kappa_thic_type /= kappa_type_const  .and.     &
          kappa_thic_type /= kappa_type_depth ) ) then   
     message = 'when kappa_isop_type is kappa_type_depth, kappa_thic_type '/&
           &/  'should be either kappa_type_const or kappa_type_depth.'
     call exit_POP(sigAbort, message)
   endif

   if ( kappa_thic_type == kappa_type_depth  .and.      &
        ( kappa_isop_type /= kappa_type_const  .and.     &
          kappa_isop_type /= kappa_type_depth ) ) then   
     message = 'when kappa_thic_type is kappa_type_depth, kappa_isop_type ' /&
            &/ 'should be either kappa_type_const or kappa_type_depth.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( ( kappa_isop_type /= kappa_type_const   .and.     &
            kappa_isop_type /= kappa_type_depth ) .and.     &
          ( kappa_thic_type /= kappa_type_const   .and.     &
            kappa_thic_type /= kappa_type_depth  ) ) .and.  &
          ( kappa_isop_type /= kappa_thic_type ) ) then
     message = 'kappa_isop_type and kappa_thic_type should be the same ' /&
           &/  'if they are BOTH model fields dependent.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( ( kappa_isop_type == kappa_type_vmhs           .and.     &
            kappa_thic_type == kappa_type_vmhs )          .or.     &
          ( kappa_isop_type == kappa_type_bfreq_vmhs     .and.     &
            kappa_thic_type == kappa_type_bfreq_vmhs )    .or.     &
          ( kappa_isop_type == kappa_type_dradius        .and.     &
            kappa_thic_type == kappa_type_dradius )       .or.     &
          ( kappa_isop_type == kappa_type_bfreq_dradius  .and.     &
            kappa_thic_type == kappa_type_bfreq_dradius ) ) .and.  &
         ah /= ah_bolus ) then
     message = 'ah and ah_bolus should be equal when vmhs or dradius '  /&
          &/ ' related dependencies are used.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( kappa_isop_type == kappa_type_depth   .or.  &
          kappa_thic_type == kappa_type_depth ) .and. & 
        kappa_depth_2 == c0 ) then
     message = 'kappa_depth_2 should be non-zero if a time-independent '  /&
          &/   'depth variation is requested.'
     call exit_POP(sigAbort, message)
      endif

   if ( ( kappa_isop_type == kappa_type_vmhs         .or.   &
          kappa_thic_type == kappa_type_vmhs         .or.   &
          kappa_isop_type == kappa_type_bfreq_vmhs   .or.   &
          kappa_thic_type == kappa_type_bfreq_vmhs   .or.   &
          kappa_isop_type == kappa_type_hdgr         .or.   &
          kappa_thic_type == kappa_type_hdgr         .or.   &
          kappa_isop_type == kappa_type_bfreq_hdgr   .or.   &
          kappa_thic_type == kappa_type_bfreq_hdgr ) .and.  &
                    zt(km) <= 2.0e5_r8 ) then
     message = 'the max bottom depth should be > 2000.0 m' /&
           &/ ' when vmhs or hdgr dependency is chosen.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( kappa_isop_type == kappa_type_eg    .or.  &
          kappa_thic_type == kappa_type_eg )  .and. &
          use_const_ah_bkg_srfbl ) then
     message = 'use_const_ah_bkg_srfbl should be false when ' /&
            &/ 'kappa_type_eg is used.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( kappa_isop_type == kappa_type_vmhs           .or.    &
          kappa_thic_type == kappa_type_vmhs           .or.    &
          kappa_isop_type == kappa_type_hdgr           .or.    &
          kappa_thic_type == kappa_type_hdgr           .or.    &
          kappa_isop_type == kappa_type_dradius        .or.    &
          kappa_thic_type == kappa_type_dradius        .or.    &
          kappa_isop_type == kappa_type_bfreq_vmhs     .or.    &
          kappa_thic_type == kappa_type_bfreq_vmhs     .or.    &
          kappa_isop_type == kappa_type_bfreq_hdgr     .or.    &
          kappa_thic_type == kappa_type_bfreq_hdgr     .or.    &
          kappa_isop_type == kappa_type_bfreq_dradius  .or.    &
          kappa_thic_type == kappa_type_bfreq_dradius  .or.    &
          kappa_isop_type == kappa_type_eg             .or.    &
          kappa_thic_type == kappa_type_eg             .or.    &
          ( ( kappa_isop_type == kappa_type_bfreq  .or.        &
              kappa_thic_type == kappa_type_bfreq ) .and.      &
            .not. read_n2_data ) )                    .and.    &
          kappa_freq == kappa_freq_never ) then

     message = 'kappa_freq should not be set to never when model ' /&
            &/ 'fields dependent kappa types are chosen.'
     call exit_POP(sigAbort, message)
   endif

   if ( partial_bottom_cells ) then
     call exit_POP(sigAbort, &
       'hmix_gm currently incompatible with partial bottom cells')
   endif

   if (read_n2_data .and. trim(buoyancy_freq_filename) == 'unknown-buoyancy' ) then
     call exit_POP(sigAbort,  &
                   'Must define buoyancy_freq_filename if read_n2_data is .true.')
   endif

!-----------------------------------------------------------------------
!
!  allocate GM arrays
!
!-----------------------------------------------------------------------

    allocate (HYXW(nx_block,ny_block,nblocks_clinic),    &
             HXYS(nx_block,ny_block,nblocks_clinic),    &
             RBR (nx_block,ny_block,nblocks_clinic),    &
             BTP (nx_block,ny_block,nblocks_clinic),    &
             BL_DEPTH(nx_block,ny_block,nblocks_clinic))

    allocate (SF_SLX(nx_block,ny_block,2,2,km,nblocks_clinic),  &
             SF_SLY(nx_block,ny_block,2,2,km,nblocks_clinic))
    
    allocate (FZTOP(nx_block,ny_block,nt,nblocks_clinic))

    allocate (kappa_depth(km))

    allocate (KAPPA_ISOP(nx_block,ny_block,2,km,nblocks_clinic),  &
             KAPPA_THIC(nx_block,ny_block,2,km,nblocks_clinic),  &
             HOR_DIFF  (nx_block,ny_block,2,km,nblocks_clinic))

    allocate (KAPPA_LATERAL (nx_block,ny_block,nblocks_clinic),  &
             KAPPA_VERTICAL(nx_block,ny_block,km,nblocks_clinic))

    allocate (BUOY_FREQ_SQ(nx_block,ny_block,km,nblocks_clinic))

    allocate (VDC_GM(nx_block,ny_block,km,nblocks_clinic))

    allocate (compute_kappa(nblocks_clinic))

   HYXW     = c0
   HXYS     = c0
   RBR      = c0
   BTP      = c0
   BL_DEPTH = c0
   SF_SLX   = c0
   SF_SLY   = c0
   FZTOP    = c0
   VDC_GM   = c0

   if ( transition_layer_on ) then
     allocate (SLA_SAVE(nx_block,ny_block,2,km,nblocks_clinic))
     allocate (RB(nx_block,ny_block,nblocks_clinic))
     SLA_SAVE = c0
     RB = c0
   endif

!-----------------------------------------------------------------------
!
!  initialize various time-independent arrays
!
!-----------------------------------------------------------------------

   do k=1,km
     kappa_depth(k) = kappa_depth_1  &
                    + kappa_depth_2  &
                     *exp(-zt(k)/kappa_depth_scale)
   enddo

   do iblock = 1,nblocks_clinic

     this_block = get_block(blocks_clinic(iblock),iblock)

     KAPPA_LATERAL(:,:,iblock)    = ah
     KAPPA_VERTICAL(:,:,:,iblock) = c1

     KAPPA_ISOP(:,:,:,:,iblock) = ah
     KAPPA_THIC(:,:,:,:,iblock) = ah_bolus
     HOR_DIFF(:,:,:,:,iblock)   = ah

     if ( kappa_isop_type == kappa_type_depth  .or.  &
          kappa_thic_type == kappa_type_depth ) then

       do k=1,km 
         KAPPA_VERTICAL(:,:,k,iblock) = kappa_depth(k)
       enddo

     endif


     HYXW(:,:,iblock) = eoshift(HYX(:,:,iblock), dim=1, shift=-1)
     HXYS(:,:,iblock) = eoshift(HXY(:,:,iblock), dim=2, shift=-1)

!-----------------------------------------------------------------------
!  compute the Rossby radius which will be used to
!  contol KAPPA [Large et al (1997), JPO, 27,
!  pp 2418-2447]. Rossby radius = c/(2*omg*sin(latitude))
!  where c=200cm/s is the first baroclinic wave speed.
!  15km < Rossby radius < 100km
!-----------------------------------------------------------------------

     !*** Inverse of Rossby radius

     RBR(:,:,iblock) = abs(FCORT(:,:,iblock))  &
                       / 200.0_r8             ! |f|/Cg, Cg = 2 m/s = 200 cm/s
     RBR(:,:,iblock) = min(RBR(:,:,iblock),    &
                           c1/1.5e+6_r8)      ! Cg/|f| .ge. 15 km = 1.5e+6 cm
     RBR(:,:,iblock) = max(RBR(:,:,iblock),    &
                           1.e-7_r8)          ! Cg/|f| .le. 100 km = 1.e+7 cm

     if ( transition_layer_on ) then
       RB(:,:,iblock) = c1 / RBR(:,:,iblock)
     endif

     !*** beta at t-points

     call ugrid_to_tgrid(BTP(:,:,iblock),ULAT(:,:,iblock),iblock)
      
     BTP(:,:,iblock) = c2*omega*cos(BTP(:,:,iblock))/radius

   enddo

!-----------------------------------------------------------------------
!  HYXW, HXYS only needed in physical domain and should
!  be defined correctly there.  BTP invalid on south
!  and westernmost ghost cells, but not needed there
!  as long as number of ghost cells is >1
!-----------------------------------------------------------------------
 
!  call update_ghost_cells(HYXW, bndy_clinic, field_loc_t,     &
!                                             field_type_scalar)
!  call update_ghost_cells(HXYS, bndy_clinic, field_loc_t,     &
!                                             field_type_scalar)
!  call update_ghost_cells(BTP , bndy_clinic, field_loc_t,     &
!                                             field_type_scalar)

   BUOY_FREQ_SQ = c0

   if ( read_n2_data ) then

     buoyancy_freq_filename = '/home/bluesky/gokhan/buoyancy_freq.nc'

     buoyancy_freq_data_file = construct_file ( 'nc',         &
                          full_name=trim(buoyancy_freq_filename) )

     call data_set(buoyancy_freq_data_file, 'open_read')

     i_dim = construct_io_dim ( 'i', nx_global )
     j_dim = construct_io_dim ( 'j', ny_global )
     k_dim = construct_io_dim ( 'k', km )

     buoyancy_freq_data_in = construct_io_field                &
                         ('BUOYANCY_FREQUENCY',                &
                          dim1=i_dim, dim2=j_dim, dim3=k_dim,  &
                          field_loc = field_loc_center,        &
                          field_type = field_type_scalar,      &
                          d3d_array = BUOY_FREQ_SQ)

     call data_set (buoyancy_freq_data_file, 'define', &
                    buoyancy_freq_data_in)
     call data_set (buoyancy_freq_data_file, 'read',   &
                    buoyancy_freq_data_in)
     call data_set (buoyancy_freq_data_file, 'close')
     call destroy_io_field (buoyancy_freq_data_in)
     call destroy_file (buoyancy_freq_data_file)

     if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,'(a43,a)')                                 &
               ' Buoyancy frequency climatology file read: ',  &
                              trim(buoyancy_freq_filename)
     endif

   endif

   compute_kappa = .false.

   if (slm_r /= slm_b) then
     diff_tapering = .true.
   else
     diff_tapering = .false.
   endif

   cancellation_occurs = .true.

   if ( ( kappa_isop_type /= kappa_thic_type )  .or.     &
       ( diff_tapering )  .or.                           &
       ( ( kappa_isop_type == kappa_type_const  .and.    &
           kappa_thic_type == kappa_type_const )  .and.  &
         ah /= ah_bolus )  .or.                          &
       ( ( kappa_isop_type == kappa_type_depth  .and.    &
           kappa_thic_type == kappa_type_depth )  .and.  &
         ah /= ah_bolus )  .or.                          &
       ( ( kappa_isop_type == kappa_type_bfreq  .and.    &
           kappa_thic_type == kappa_type_bfreq )  .and.  &
         ah /= ah_bolus ) )                              &
    cancellation_occurs = .false. 

    !*** for transition layer cases, the following will always be true!!

   if ( transition_layer_on )  cancellation_occurs = .false.

!-----------------------------------------------------------------------
!
!  initialize topography mask used with kappa_type_eg. This mask eliminates
!  excessively large values of KAPPA near sloping topography.
!
!-----------------------------------------------------------------------

   if ( kappa_isop_type == kappa_type_eg  .or.  &
        kappa_thic_type == kappa_type_eg ) then 

     allocate (SIGMA_TOPO_MASK(nx_block,ny_block,km,nblocks_clinic))

     do iblock=1,nblocks_clinic

       do k=1,km
         where ( k < KMT(:,:,iblock) ) 
           SIGMA_TOPO_MASK(:,:,k,iblock) = c1
         elsewhere
           SIGMA_TOPO_MASK(:,:,k,iblock) = c0
         endwhere
       enddo

       do k=1,km-1
         do j=2,ny_block-1
           do i=2,nx_block-1 
             if ( k < KMT(i,j,iblock) ) then
               if ( k == KMT(i-1,j+1,iblock)  .or.  &
                    k == KMT(i  ,j+1,iblock)  .or.  &
                    k == KMT(i+1,j+1,iblock)  .or.  &
                    k == KMT(i-1,j  ,iblock)  .or.  &
                    k == KMT(i+1,j  ,iblock)  .or.  &
                    k == KMT(i-1,j-1,iblock)  .or.  &
                    k == KMT(i  ,j-1,iblock)  .or.  &
                    k == KMT(i+1,j-1,iblock) )      &
                 SIGMA_TOPO_MASK(i,j,k,iblock) = c0 
             endif
           enddo
         enddo
       enddo

     enddo

   endif

!-----------------------------------------------------------------------
!
!  define tavg fields related to diffusivities 
!
!-----------------------------------------------------------------------

   call define_tavg_field (tavg_KAPPA_ISOP, 'KAPPA_ISOP', 3,     &
                   long_name='Isopycnal diffusion coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_KAPPA_THIC, 'KAPPA_THIC', 3,     &
                   long_name='Thickness diffusion coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_HOR_DIFF, 'HOR_DIFF', 3,         &
                   long_name='Horizontal diffusion coefficient', &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!
!  define tavg fields related to the near-surface eddy flux
!  parameterization (i.e. transition layer is on). 
!
!-----------------------------------------------------------------------

   if ( transition_layer_on ) then

     call define_tavg_field (tavg_DIA_DEPTH, 'DIA_DEPTH', 2,       &
         long_name='Depth of the Diabatic Region at the Surface',  &
                     units='cm', grid_loc='2110',                  &
                     coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_TLT, 'TLT', 2,                   &
         long_name='Transition Layer Thickness',                   &
                     units='cm', grid_loc='2110',                  &
                     coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_INT_DEPTH, 'INT_DEPTH', 2,       &
         long_name='Depth at which the Interior Region Starts',    &
                     units='cm', grid_loc='2110',                  &
                     coordinates='TLONG TLAT time')

   endif

!-----------------------------------------------------------------------
!
!  allocate and initialize bolus velocity arrays if requested
!  define tavg fields related to bolus velocity
!
!-----------------------------------------------------------------------

   if ( diag_gm_bolus ) then

     call register_string ('diag_gm_bolus')

     allocate(WTOP_ISOP(nx_block,ny_block,nblocks_clinic), &
              WBOT_ISOP(nx_block,ny_block,nblocks_clinic), &
                    UIT(nx_block,ny_block,nblocks_clinic), &
                    VIT(nx_block,ny_block,nblocks_clinic))

     WTOP_ISOP = c0
     WBOT_ISOP = c0
     UIT       = c0
     VIT       = c0

     call define_tavg_field (tavg_UISOP, 'UISOP', 3,                &
      long_name='Bolus Velocity in grid-x direction (diagnostic)',  &
                   units='cm/s', grid_loc='3211',                   &
                   coordinates='ULONG TLAT z_t time')

     call define_tavg_field (tavg_VISOP, 'VISOP', 3,                &
      long_name='Bolus Velocity in grid-y direction (diagnostic)',  &
                   units='cm/s', grid_loc='3121',                   &
                   coordinates='TLONG ULAT z_t time')

     call define_tavg_field (tavg_WISOP, 'WISOP', 3,                &
      long_name='Vertical Bolus Velocity (diagnostic)',             &
                   units='cm/s', grid_loc='3112',                   &
                   coordinates='TLONG TLAT z_w time')

     call define_tavg_field (tavg_ADVT_ISOP, 'ADVT_ISOP', 2,                            &
      long_name='Vertically-Integrated T Eddy-Induced Advection Tendency (diagnostic)', &
                   units='cm degC/s', grid_loc='2110',                                  &
                   coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_ADVS_ISOP, 'ADVS_ISOP', 2,                            &
      long_name='Vertically-Integrated S Eddy-Induced Advection Tendency (diagnostic)', &
                   scale_factor=1000.0_rtavg,                                           &
                   units='cm gram/kilogram/s', grid_loc='2110',                         &
                   coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_VNT_ISOP, 'VNT_ISOP', 3,                               &
      long_name='Heat Flux Tendency in grid-y Dir due to Eddy-Induced Vel (diagnostic)', &
                   units='degC/s', grid_loc='3121',                                      &
                   coordinates='TLONG ULAT z_t time')

     call define_tavg_field (tavg_VNS_ISOP, 'VNS_ISOP', 3,                               &
      long_name='Salt Flux Tendency in grid-y Dir due to Eddy-Induced Vel (diagnostic)', &
                   scale_factor=1000.0_rtavg,                                            &
                   units='gram/kilogram/s', grid_loc='3121',                             &
                   coordinates='TLONG ULAT z_t time')

   endif

      call get_timer(timer_nloop,'HMIX_TRACER_GM_NLOOP', &
                                  nblocks_clinic, distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!EOC

   end subroutine init_gm

!***********************************************************************
!BOP
! !IROUTINE: hdifft_gm
! !INTERFACE:

      subroutine hdifft_gm (k, GTK, TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
                            tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

! !DESCRIPTION:
!  Gent-McWilliams eddy transport parameterization
!  and isopycnal diffusion.
!
!  This routine must be called successively with k = 1,2,3,...
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      integer (int_kind), intent(in) :: k  ! depth level index

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX            ! U,V  at all vertical levels
                               !   at mixing time level

      integer (int_kind), dimension(nt), intent(in) :: &
         tavg_HDIFE_TRACER, &! tavg id for east face diffusive flux of tracer
         tavg_HDIFN_TRACER, &! tavg id for north face diffusive flux of tracer
         tavg_HDIFB_TRACER   ! tavg id for bottom face diffusive flux of tracer

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

! !OUTPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,nt), intent(out) :: &
         GTK     ! diffusion+bolus advection for nth tracer at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ieast  = 1, iwest  = 2,       &
         jnorth = 1, jsouth = 2

      integer (int_kind) :: &
         i,j,n,kk,          &! dummy loop counters
         kid,ktmp,          &! array indices
         kk_sub, kp1,       & 
         bid                 ! local block address for this sub block

      real (r8) :: &
         fz, dz_bottom, factor

      real (r8), dimension(nx_block,ny_block) :: &
         CX, CY,                  &
         RZ,                      &! Dz(rho)
         SLA,                     &! absolute value of slope
         WORK1, WORK2,            &! local work space
         WORK3, WORK4,            &! local work space
         KMASK,                   &! ocean mask
         TAPER1, TAPER2, TAPER3,  &! tapering factors
         UIB, VIB,                &! work arrays for isopycnal mixing velocities
         U_ISOP, V_ISOP            ! horizontal components of isopycnal velocities

      real (r8), dimension(nx_block,ny_block,nt) :: &
         FX, FY                     ! fluxes across east, north faces

      real (r8), dimension(2) :: &
         reference_depth

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      U_ISOP = c0
      V_ISOP = c0
      WORK1  = c0
      WORK2  = c0
      WORK3  = c0
      WORK4  = c0

      if ( .not. implicit_vertical_mix )  &
        call exit_POP (sigAbort, &
         'implicit vertical mixing must be used with GM horiz mixing')

      if ( k == 1 ) then

        if ( diag_gm_bolus ) then
          UIB = c0
          VIB = c0
          UIT(:,:,bid) = c0
          VIT(:,:,bid) = c0
          WBOT_ISOP(:,:,bid) = c0
        endif

        HOR_DIFF(:,:,ktp,k,bid) = ah_bkg_srfbl

        BL_DEPTH(:,:,bid) = zw(k)
        if ( vmix_itype == vmix_type_kpp )  &
                    BL_DEPTH(:,:,bid) = KPP_HBLT(:,:,bid)

      endif

      if (diag_gm_bolus)  WTOP_ISOP(:,:,bid) = WBOT_ISOP(:,:,bid)

      CX = merge(HYX(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTE(:,:,bid)))
      CY = merge(HXY(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTN(:,:,bid)))

      if ( k == 1 ) then

        if ( transition_layer_on ) then

          if ( vmix_itype == vmix_type_kpp ) then

            call smooth_hblt ( .false., .true., bid,  &
                               SMOOTH_OUT=TLT%DIABATIC_DEPTH(:,:,bid) )
          else
            TLT%DIABATIC_DEPTH(:,:,bid) = zw(k)
          endif

          do kk=1,km
            do kk_sub = ktp,kbt
              kid = kk + kk_sub - 2
              SLA_SAVE(:,:,kk_sub,kk,bid) = dzw(kid)*sqrt(p5*(       &
                     (SLX(:,:,1,kk_sub,kk,bid)**2                    &
                    + SLX(:,:,2,kk_sub,kk,bid)**2)/DXT(:,:,bid)**2   &
                   + (SLY(:,:,1,kk_sub,kk,bid)**2                    &
                    + SLY(:,:,2,kk_sub,kk,bid)**2)/DYT(:,:,bid)**2)) &
                   + eps
            enddo
          enddo

          call transition_layer ( this_block )

        endif

	

!-----------------------------------------------------------------------
!
!     compute isopycnal and thickness diffusion coefficients if
!     they depend on the model fields
!
!-----------------------------------------------------------------------
        
	if ( ( kappa_isop_type == kappa_type_vmhs           .or.    &
               kappa_thic_type == kappa_type_vmhs           .or.    &
               kappa_isop_type == kappa_type_hdgr           .or.    &
               kappa_thic_type == kappa_type_hdgr           .or.    &
               kappa_isop_type == kappa_type_dradius        .or.    &
               kappa_thic_type == kappa_type_dradius        .or.    &
               kappa_isop_type == kappa_type_bfreq          .or.    &
               kappa_thic_type == kappa_type_bfreq          .or.    &
               kappa_isop_type == kappa_type_bfreq_vmhs     .or.    &
               kappa_thic_type == kappa_type_bfreq_vmhs     .or.    &
               kappa_isop_type == kappa_type_bfreq_hdgr     .or.    &
               kappa_thic_type == kappa_type_bfreq_hdgr     .or.    &
               kappa_isop_type == kappa_type_bfreq_dradius  .or.    &
               kappa_thic_type == kappa_type_bfreq_dradius  .or.    &
               kappa_isop_type == kappa_type_eg             .or.    &
               kappa_thic_type == kappa_type_eg )            .and.  &
           ( ( kappa_freq == kappa_freq_every_time_step )           &
        .or. ( kappa_freq == kappa_freq_once_a_day .and. eod_last ) &
        .or. ( nsteps_total == 1 ) ) )  compute_kappa(bid) = .true.

        if ( compute_kappa(bid) ) then

          if ( kappa_isop_type == kappa_type_vmhs        .or.  &
               kappa_thic_type == kappa_type_vmhs        .or.  &
               kappa_isop_type == kappa_type_bfreq_vmhs  .or.  & 
               kappa_thic_type == kappa_type_bfreq_vmhs ) then

            if ( nsteps_total == 1 ) then
              KAPPA_LATERAL(:,:,bid) = ah
              if ( kappa_isop_type == kappa_type_const )  &
                KAPPA_LATERAL(:,:,bid) = ah_bolus
            else
              call kappa_lon_lat_vmhs (TMIX, UMIX, VMIX, this_block)
            endif

          endif

          if ( kappa_isop_type == kappa_type_hdgr        .or.  &
               kappa_thic_type == kappa_type_hdgr        .or.  &
               kappa_isop_type == kappa_type_bfreq_hdgr  .or.  &
               kappa_thic_type == kappa_type_bfreq_hdgr )      &
            call kappa_lon_lat_hdgr (TMIX, this_block)

          if ( kappa_isop_type == kappa_type_dradius        .or.  &
               kappa_thic_type == kappa_type_dradius        .or.  &
               kappa_isop_type == kappa_type_bfreq_dradius  .or.  &
               kappa_thic_type == kappa_type_bfreq_dradius )      &
            call kappa_lon_lat_dradius (this_block)

          if ( kappa_isop_type == kappa_type_bfreq          .or.  &
               kappa_thic_type == kappa_type_bfreq          .or.  &
               kappa_isop_type == kappa_type_bfreq_vmhs     .or.  &
               kappa_thic_type == kappa_type_bfreq_vmhs     .or.  &
               kappa_isop_type == kappa_type_bfreq_hdgr     .or.  &
               kappa_thic_type == kappa_type_bfreq_hdgr     .or.  &
               kappa_isop_type == kappa_type_bfreq_dradius  .or.  &
               kappa_thic_type == kappa_type_bfreq_dradius )      &
            call buoyancy_frequency_dependent_profile (TMIX, this_block)

          if ( kappa_isop_type == kappa_type_eg  .or.  &
               kappa_thic_type == kappa_type_eg ) &
            call kappa_eg (TMIX, UMIX, VMIX, this_block) 

          compute_kappa(bid) = .false.

        endif  ! end of ( compute_kappa ) if statement


!-----------------------------------------------------------------------
!
!     reinitialize the diffusivity coefficients 
!
!-----------------------------------------------------------------------
	
        if ( kappa_isop_type == kappa_type_const ) then
          KAPPA_ISOP(:,:,:,:,bid) = ah
        elseif ( kappa_isop_type == kappa_type_eg ) then
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_ISOP(:,:,kk_sub,kk,bid) = KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo
        else
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_ISOP(:,:,kk_sub,kk,bid) =  KAPPA_LATERAL(:,:,bid)  &
                                         * KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo
        endif

        if ( .not. use_const_ah_bkg_srfbl )  &
          HOR_DIFF(:,:,ktp,k,bid) = KAPPA_ISOP(:,:,ktp,k,bid) 

        if ( kappa_thic_type == kappa_type_const ) then
          KAPPA_THIC(:,:,:,:,bid) = ah_bolus
        else if ( kappa_thic_type == kappa_type_depth  .or.  &
                  kappa_thic_type == kappa_type_bfreq ) then
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_THIC(:,:,kk_sub,kk,bid) =  ah_bolus  &
                                        * KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo 
        else if ( kappa_thic_type == kappa_type_eg ) then
          KAPPA_THIC(:,:,:,:,bid) = KAPPA_ISOP(:,:,:,:,bid)
        else
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_THIC(:,:,kk_sub,kk,bid) = KAPPA_LATERAL(:,:,bid)  &
                                        * KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo
        endif


!-----------------------------------------------------------------------
!
!     control slope of isopycnal surfaces or KAPPA
!
!-----------------------------------------------------------------------
	
        do kk=1,km

          kp1 = min(kk+1,km)
          reference_depth(ktp) = zt(kp1)
          reference_depth(kbt) = zw(kp1)
          if ( kk == km )  reference_depth(ktp) = zw(kp1)

          do kk_sub = ktp,kbt 

            kid = kk + kk_sub - 2

!-----------------------------------------------------------------------
!
!     control KAPPA to reduce the isopycnal mixing near the
!     ocean surface Large et al (1997), JPO, 27, pp 2418-2447.
!     WORK1 = ratio between the depth of water parcel and
!     the vertical displacement of isopycnal surfaces
!     where the vertical displacement =
!     Rossby radius * slope of isopycnal surfaces
!
!-----------------------------------------------------------------------

            if ( transition_layer_on ) then
              SLA = SLA_SAVE(:,:,kk_sub,kk,bid)
            else
              SLA = dzw(kid)*sqrt(p5*(                               &
                     (SLX(:,:,1,kk_sub,kk,bid)**2                    & 
                    + SLX(:,:,2,kk_sub,kk,bid)**2)/DXT(:,:,bid)**2   &
                   + (SLY(:,:,1,kk_sub,kk,bid)**2                    &
                    + SLY(:,:,2,kk_sub,kk,bid)**2)/DYT(:,:,bid)**2)) &
                   + eps
            endif

            TAPER1 = c1 
            if ( .not. transition_layer_on ) then

              if ( kk == 1 ) then
                dz_bottom = c0
              else
                dz_bottom = zt(kk-1)
              endif

              if (slope_control == slope_control_tanh) then

                WORK1 = min(c1,zt(kk)*RBR(:,:,bid)/SLA)
                TAPER1 = p5*(c1+sin(pi*(WORK1-p5)))

!     use the Rossby deformation radius tapering
!     only within the boundary layer

                TAPER1 = merge(TAPER1, c1,  &
                               dz_bottom <= BL_DEPTH(:,:,bid))

              else

!     sine function is replaced by
!     function = 4.*x*(1.-abs(x)) for |x|<0.5

                WORK1 = min(c1,zt(kk)*RBR(:,:,bid)/SLA)
                TAPER1 = (p5+c2*(WORK1-p5)*(c1-abs(WORK1-p5)))

                TAPER1 = merge(TAPER1, c1,  &
                               dz_bottom <= BL_DEPTH(:,:,bid))

              endif

            endif

!-----------------------------------------------------------------------
!
!     control KAPPA for numerical stability
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     methods to control slope
!
!-----------------------------------------------------------------------

            TAPER2 = c1 
            TAPER3 = c1 

            select case (slope_control)
            case (slope_control_tanh)

!     method by Danabasoglu & Mcwilliams (1995)

              TAPER2 = merge(p5*  &
                 (c1-tanh(c10*SLA/slm_r-c4)), c0, SLA < slm_r)

              if ( diff_tapering ) then
                TAPER3 = merge(p5*  &
                 (c1-tanh(c10*SLA/slm_b-c4)), c0, SLA < slm_b)
              else
                TAPER3 = TAPER2
              endif

            case (slope_control_notanh)

!     similar to DM95 except replacing tanh by
!     function = x*(1.-0.25*abs(x)) for |x|<2
!              = sign(x)            for |x|>2
!     (faster than DM95)

              do j=1,ny_block
                do i=1,nx_block
                  if (SLA(i,j) > 0.2_r8*slm_r .and. &
                      SLA(i,j) < 0.6_r8*slm_r) then
                    TAPER2(i,j) = &
                       p5*(c1-(2.5_r8*SLA(i,j)/slm_r-c1)*  &
                          (c4-abs(c10*SLA(i,j)/slm_r-c4)))
                  else if (SLA(i,j) >= 0.6_r8*slm_r) then
                    TAPER2(i,j) = c0
                  endif
                enddo
              enddo

              if ( diff_tapering ) then
                do j=1,ny_block
                  do i=1,nx_block
                    if (SLA(i,j) > 0.2_r8*slm_b .and. &
                        SLA(i,j) < 0.6_r8*slm_b) then
                      TAPER3(i,j) = &
                         p5*(c1-(2.5_r8*SLA(i,j)/slm_b-c1)* &
                            (c4-abs(c10*SLA(i,j)/slm_b-c4)))
                    else if (SLA(i,j) >= 0.6_r8*slm_b) then
                      TAPER3(i,j) = c0
                    endif
                  end do
                end do
              else
                TAPER3 = TAPER2
              endif

            case (slope_control_clip)

!     slope clipping

              do n=1,2
                do j=1,ny_block
                  do i=1,nx_block
                    if (abs(SLX(i,j,n,kk_sub,kk,bid)  &
                          * dzw(kid) / HUS(i,j,bid)) > slm_r) then
                      SLX(i,j,n,kk_sub,kk,bid) =             &
                                  sign(slm_r * HUS(i,j,bid)  &
                                     * dzwr(kid),            &
                                  SLX(i,j,n,kk_sub,kk,bid))
                    endif
                  enddo
                enddo
              enddo

              do n=1,2
                do j=1,ny_block
                  do i=1,nx_block
                    if (abs(SLY(i,j,n,kk_sub,kk,bid)  &
                          * dzw(kid) / HUW(i,j,bid)) > slm_r) then  
                      SLY(i,j,n,kk_sub,kk,bid) =             &
                                  sign(slm_r * HUW(i,j,bid)  &
                                     * dzwr(kid),            &
                                  SLY(i,j,n,kk_sub,kk,bid))
                    endif
                  enddo
                enddo
              enddo

            case (slope_control_Gerd)

!     method by Gerdes et al (1991)

              do j=1,ny_block
                do i=1,nx_block
                  if (SLA(i,j) > slm_r)  &
                    TAPER2(i,j) = (slm_r/SLA(i,j))**2
                enddo
              enddo

              if (diff_tapering) then
                do j=1,ny_block
                  do i=1,nx_block
                    if (SLA(i,j) > slm_b)  &
                      TAPER3(i,j) = (slm_b/SLA(i,j))**2
                  enddo
                enddo
              else
                TAPER3 = TAPER2
              endif

            end select

            if ( transition_layer_on ) then
              TAPER2 = merge(c1, TAPER2, reference_depth(kk_sub) &
                                         <= TLT%DIABATIC_DEPTH(:,:,bid))
              TAPER3 = merge(c1, TAPER3, reference_depth(kk_sub) &
                                         <= TLT%DIABATIC_DEPTH(:,:,bid))
            endif

            if ( transition_layer_on  .and.  use_const_ah_bkg_srfbl ) then

              HOR_DIFF(:,:,kk_sub,kk,bid) = ah_bkg_srfbl

            else if ( transition_layer_on .and.               &
                    ( .not. use_const_ah_bkg_srfbl      .or.  &
                      kappa_isop_type == kappa_type_eg  .or.  &
                      kappa_thic_type == kappa_type_eg ) ) then

              HOR_DIFF(:,:,kk_sub,kk,bid) = KAPPA_ISOP(:,:,kk_sub,kk,bid) 

            else

              if ( .not. ( kk == 1 .and. kk_sub == ktp ) ) then 

                if ( use_const_ah_bkg_srfbl ) then 
                  HOR_DIFF(:,:,kk_sub,kk,bid) =                       &
                       merge( ah_bkg_srfbl * (c1 - TAPER1 * TAPER2)   &
                              * KAPPA_VERTICAL(:,:,kk,bid),           &
                              c0, dz_bottom <= BL_DEPTH(:,:,bid) )
                else
                  HOR_DIFF(:,:,kk_sub,kk,bid) =                &
                         merge( KAPPA_ISOP(:,:,kk_sub,kk,bid)  &
                               * (c1 - TAPER1 * TAPER2),       &
                              c0, dz_bottom <= BL_DEPTH(:,:,bid) )
                endif

              endif

            endif

            KAPPA_ISOP(:,:,kk_sub,kk,bid) =  &
                  TAPER1 * TAPER2 * KAPPA_ISOP(:,:,kk_sub,kk,bid)
            KAPPA_THIC(:,:,kk_sub,kk,bid) =  &
                  TAPER1 * TAPER3 * KAPPA_THIC(:,:,kk_sub,kk,bid)

          end do  ! end of kk_sub loop


!-----------------------------------------------------------------------
!
!     impose the boundary conditions by setting KAPPA=0
!     in the quarter cells adjacent to rigid boundaries.
!     bottom B.C.s are considered within the kk-loop.
!
!-----------------------------------------------------------------------

!     B.C. at the bottom

          where (kk == KMT(:,:,bid))
            KAPPA_ISOP(:,:,kbt,kk,bid) = c0
            KAPPA_THIC(:,:,kbt,kk,bid) = c0
          end where

        enddo              ! end of kk-loop

!     B.C. at the top

        KAPPA_ISOP(:,:,ktp,1,bid) = c0
        KAPPA_THIC(:,:,ktp,1,bid) = c0

        FZTOP(:,:,:,bid) = c0        ! zero flux B.C. at the surface

        if ( transition_layer_on ) then

          call merged_streamfunction ( this_block )

          call apply_vertical_profile_to_isop_hor_diff ( this_block ) 

        else

          TLT%DIABATIC_DEPTH(:,:,bid) = c0
          TLT%THICKNESS(:,:,bid)      = c0
          TLT%INTERIOR_DEPTH(:,:,bid) = c0

          do kk=1,km

            do kk_sub=ktp,kbt

              where ( kk <= KMT(:,:,bid) ) 

                SF_SLX(:,:,1,kk_sub,kk,bid) =          &
                        KAPPA_THIC(:,:,kk_sub,kk,bid)  &
                      * SLX(:,:,1,kk_sub,kk,bid) * dz(kk)

                SF_SLX(:,:,2,kk_sub,kk,bid) =          &
                        KAPPA_THIC(:,:,kk_sub,kk,bid)  &
                      * SLX(:,:,2,kk_sub,kk,bid) * dz(kk)

                SF_SLY(:,:,1,kk_sub,kk,bid) =          &
                        KAPPA_THIC(:,:,kk_sub,kk,bid)  &
                      * SLY(:,:,1,kk_sub,kk,bid) * dz(kk)

                SF_SLY(:,:,2,kk_sub,kk,bid) =          &
                        KAPPA_THIC(:,:,kk_sub,kk,bid)  &
                      * SLY(:,:,2,kk_sub,kk,bid) * dz(kk)

              endwhere

            enddo  ! end of kk_sub-loop

          enddo    ! end of kk-loop

        endif

      endif  ! end of k==1 if statement

      KMASK = merge(c1, c0, k < KMT(:,:,bid))
      

!-----------------------------------------------------------------------
!
!     calculate effective vertical diffusion coefficient
!     NOTE: it is assumed that VDC has been set before this
!           in vmix_coeffs or something similar.
!
!     Dz(VDC * Dz(T)) where D is derivative rather than difference
!     VDC = (Az(dz*Ax(KAPPA*HYX*SLX**2)) + Az(dz*Ay(KAPPA*HXY*SLY**2)))*
!           dzw/TAREA
!
!-----------------------------------------------------------------------

      
      if ( k < km ) then

        WORK1 = dzw(k)*KMASK*TAREA_R(:,:,bid)*                &
                (dz(k  )*p25*KAPPA_ISOP(:,:,kbt,k,  bid)*     &
              (HYX (:,:,bid)*SLX(:,:,ieast, kbt,k,  bid)**2   &
             + HYXW(:,:,bid)*SLX(:,:,iwest, kbt,k,  bid)**2   &
             + HXY (:,:,bid)*SLY(:,:,jnorth,kbt,k,  bid)**2   &
             + HXYS(:,:,bid)*SLY(:,:,jsouth,kbt,k,  bid)**2)  &
                +dz(k+1)*p25*KAPPA_ISOP(:,:,ktp,k+1,bid)*     &
              (HYX (:,:,bid)*SLX(:,:,ieast, ktp,k+1,bid)**2   &
             + HYXW(:,:,bid)*SLX(:,:,iwest, ktp,k+1,bid)**2   &
             + HXY (:,:,bid)*SLY(:,:,jnorth,ktp,k+1,bid)**2   &
             + HXYS(:,:,bid)*SLY(:,:,jsouth,ktp,k+1,bid)**2))

        do n=1,size(VDC,DIM=4)
          VDC_GM(:,:,k,bid) = WORK1
          VDC(:,:,k,n,bid) = VDC(:,:,k,n,bid) + WORK1
        end do

      end if

!-----------------------------------------------------------------------
!
!     check if some horizontal diffusion needs to be added to the
!     bottom half of the bottom cell
!
!-----------------------------------------------------------------------

      if ( ah_bkg_bottom /= c0 ) then
        where ( k == KMT(:,:,bid) ) 
          HOR_DIFF(:,:,kbt,k,bid) = ah_bkg_bottom
        endwhere
      endif

!-----------------------------------------------------------------------
!
!     combine isopycnal and horizontal diffusion coefficients
!
!-----------------------------------------------------------------------

      do j=1,ny_block
        do i=1,nx_block-1
          WORK3(i,j) = KAPPA_ISOP(i,  j,ktp,k,bid)  &
                     + HOR_DIFF  (i,  j,ktp,k,bid)  &
                     + KAPPA_ISOP(i,  j,kbt,k,bid)  &
                     + HOR_DIFF  (i,  j,kbt,k,bid)  &
                     + KAPPA_ISOP(i+1,j,ktp,k,bid)  &
                     + HOR_DIFF  (i+1,j,ktp,k,bid)  &
                     + KAPPA_ISOP(i+1,j,kbt,k,bid)  &
                     + HOR_DIFF  (i+1,j,kbt,k,bid)
        enddo
      enddo
      
      do j=1,ny_block-1
        do i=1,nx_block
          WORK4(i,j) = KAPPA_ISOP(i,j,  ktp,k,bid)  &
                     + HOR_DIFF  (i,j,  ktp,k,bid)  & 
                     + KAPPA_ISOP(i,j,  kbt,k,bid)  &
                     + HOR_DIFF  (i,j,  kbt,k,bid)  &
                     + KAPPA_ISOP(i,j+1,ktp,k,bid)  &
                     + HOR_DIFF  (i,j+1,ktp,k,bid)  &
                     + KAPPA_ISOP(i,j+1,kbt,k,bid)  &
                     + HOR_DIFF  (i,j+1,kbt,k,bid)
        enddo
      enddo
      

!-----------------------------------------------------------------------
!
!     start loop over tracers
!
!-----------------------------------------------------------------------

      kp1 = k + 1
      if ( k == km )  kp1 = k

      if ( k < km ) then
        dz_bottom = dz(kp1)
        factor    = c1
      else
        dz_bottom = c0
        factor    = c0
      endif

      call timer_start(timer_nloop, block_id=this_block%local_id)

      do n = 1,nt

!-----------------------------------------------------------------------
!
!     calculate horizontal fluxes thru vertical faces of T-cell
!     FX = dz*HYX*Ax(Az(KAPPA))*Dx(T) : flux in x-direction
!     FY = dz*HXY*Ay(Az(KAPPA))*Dy(T) : flux in y-direction
!
!-----------------------------------------------------------------------

        FX(:,:,n) = dz(k) * CX * TX(:,:,k,n,bid) * WORK3
        FY(:,:,n) = dz(k) * CY * TY(:,:,k,n,bid) * WORK4 

      end do

      if ( .not. cancellation_occurs ) then

        do j=1,ny_block
          do i=1,nx_block-1
            WORK1(i,j) = KAPPA_ISOP(i,j,ktp,k,bid)                     &
                         * SLX(i,j,ieast,ktp,k,bid) * dz(k)            &
                         - SF_SLX(i,j,ieast,ktp,k,bid)
            WORK2(i,j) = KAPPA_ISOP(i,j,kbt,k,bid)                     &
                         * SLX(i,j,ieast,kbt,k,bid) * dz(k)            &
                         - SF_SLX(i,j,ieast,kbt,k,bid)
            WORK3(i,j) = KAPPA_ISOP(i+1,j,ktp,k,bid)                   &
                         * SLX(i+1,j,iwest,ktp,k,bid) * dz(k)          &
                         - SF_SLX(i+1,j,iwest,ktp,k,bid)
            WORK4(i,j) = KAPPA_ISOP(i+1,j,kbt,k,bid)                   &
                         * SLX(i+1,j,iwest,kbt,k,bid) * dz(k)          &
                         - SF_SLX(i+1,j,iwest,kbt,k,bid)
          enddo
        enddo
	
	do n = 1,nt

          if (n > 2 .and. k < km)  &
            TZ(:,:,k+1,n,bid) = TMIX(:,:,k  ,n) - TMIX(:,:,k+1,n)

          do j=1,ny_block
            do i=1,nx_block-1
              FX(i,j,n) = FX(i,j,n) - CX(i,j)                          &
               * ( WORK1(i,j) * TZ(i,j,k,n,bid)                        &
                   + WORK2(i,j) * TZ(i,j,kp1,n,bid)                    &
                   + WORK3(i,j) * TZ(i+1,j,k,n,bid)                    &
                   + WORK4(i,j) * TZ(i+1,j,kp1,n,bid) )
            enddo
          enddo

        end do


        do j=1,ny_block-1
          do i=1,nx_block
            WORK1(i,j) = KAPPA_ISOP(i,j,ktp,k,bid)                     &
                         * SLY(i,j,jnorth,ktp,k,bid) * dz(k)           &
                         - SF_SLY(i,j,jnorth,ktp,k,bid)
            WORK2(i,j) = KAPPA_ISOP(i,j,kbt,k,bid)                     &
                         * SLY(i,j,jnorth,kbt,k,bid) * dz(k)           &
                         - SF_SLY(i,j,jnorth,kbt,k,bid)
            WORK3(i,j) = KAPPA_ISOP(i,j+1,ktp,k,bid)                   &
                         * SLY(i,j+1,jsouth,ktp,k,bid) * dz(k)         &
                         - SF_SLY(i,j+1,jsouth,ktp,k,bid)
            WORK4(i,j) = KAPPA_ISOP(i,j+1,kbt,k,bid)                   &
                         * SLY(i,j+1,jsouth,kbt,k,bid) * dz(k)         &
                         - SF_SLY(i,j+1,jsouth,kbt,k,bid)
          enddo
        enddo

        do n = 1,nt

          do j=1,ny_block-1
            do i=1,nx_block
              FY(i,j,n) = FY(i,j,n) - CY(i,j)                          &
               * ( WORK1(i,j) * TZ(i,j,k,n,bid)                        &
                   + WORK2(i,j) * TZ(i,j,kp1,n,bid)                    &
                   + WORK3(i,j) * TZ(i,j+1,k,n,bid)                    &
                   + WORK4(i,j) * TZ(i,j+1,kp1,n,bid) )
            enddo
          enddo

        end do

      endif

      do n = 1,nt

!-----------------------------------------------------------------------
!
!     calculate vertical fluxes thru horizontal faces of T-cell
!     - Az(dz*Ax(HYX*KAPPA*SLX*TX)) - Az(dz*Ay(HXY*KAPPA*SLY*TY))
!     calculate isopycnal diffusion from flux differences
!     DTK = (Dx(FX)+Dy(FY)+Dz(FZ)) / volume
!
!-----------------------------------------------------------------------

        GTK(:,:,n) = c0

        if ( k < km ) then
          
          WORK3 = c0

          if ( .not. cancellation_occurs ) then

!pw loop split to improve performance  -- 2

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                WORK3(i,j) = WORK3(i,j)                               &
                    + ( dz(k) * KAPPA_ISOP(i,j,kbt,k,bid)             &
                    * ( SLX(i,j,ieast ,kbt,k  ,bid)                   &
                       * HYX(i  ,j  ,bid) * TX(i  ,j  ,k  ,n,bid)     &
                      + SLY(i,j,jnorth,kbt,k  ,bid)                   &
                       * HXY(i  ,j  ,bid) * TY(i  ,j  ,k  ,n,bid)     &
                      + SLX(i,j,iwest ,kbt,k  ,bid)                   &
                       * HYX(i-1,j  ,bid) * TX(i-1,j  ,k  ,n,bid)     &
                      + SLY(i,j,jsouth,kbt,k  ,bid)                   &
                       * HXY(i  ,j-1,bid) * TY(i  ,j-1,k  ,n,bid) ) )

               enddo
             enddo

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                WORK3(i,j) = WORK3(i,j)                               &
                    + ( SF_SLX(i  ,j  ,ieast ,kbt,k  ,bid)            &
                       * HYX(i  ,j  ,bid) * TX(i  ,j  ,k  ,n,bid)     &
                      + SF_SLY(i  ,j  ,jnorth,kbt,k  ,bid)            &
                       * HXY(i  ,j  ,bid) * TY(i  ,j  ,k  ,n,bid)     &
                      + SF_SLX(i  ,j  ,iwest ,kbt,k  ,bid)            &
                       * HYX(i-1,j  ,bid) * TX(i-1,j  ,k  ,n,bid)     &
                      + SF_SLY(i  ,j  ,jsouth,kbt,k  ,bid)            &
                       * HXY(i  ,j-1,bid) * TY(i  ,j-1,k  ,n,bid) )

               enddo
             enddo

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                WORK3(i,j) = WORK3(i,j)                               &
                    + ( dz_bottom * KAPPA_ISOP(i,j,ktp,kp1,bid)       &
                    * ( SLX(i  ,j  ,ieast ,ktp,kp1,bid)               &
                       * HYX(i  ,j  ,bid) * TX(i  ,j  ,kp1,n,bid)     &
                      + SLY(i  ,j  ,jnorth,ktp,kp1,bid)               &
                       * HXY(i  ,j  ,bid) * TY(i  ,j  ,kp1,n,bid)     &
                      + SLX(i  ,j  ,iwest ,ktp,kp1,bid)               &
                       * HYX(i-1,j  ,bid) * TX(i-1,j  ,kp1,n,bid)     &
                      + SLY(i  ,j  ,jsouth,ktp,kp1,bid)               &
                       * HXY(i  ,j-1,bid) * TY(i  ,j-1,kp1,n,bid) ) )

               enddo
             enddo

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                 WORK3(i,j) = WORK3(i,j)                              &
                    + ( factor                                        &
                    * ( SF_SLX(i  ,j  ,ieast ,ktp,kp1,bid)            &
                       * HYX(i  ,j  ,bid) * TX(i  ,j  ,kp1,n,bid)     &
                      + SF_SLY(i  ,j  ,jnorth,ktp,kp1,bid)            &
                       * HXY(i  ,j  ,bid) * TY(i  ,j  ,kp1,n,bid)     &
                      + SF_SLX(i  ,j  ,iwest ,ktp,kp1,bid)            &
                       * HYX(i-1,j  ,bid) * TX(i-1,j  ,kp1,n,bid)     &
                      + SF_SLY(i  ,j  ,jsouth,ktp,kp1,bid)            &
                       * HXY(i  ,j-1,bid) * TY(i  ,j-1,kp1,n,bid) ) )

               enddo
             enddo

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                fz = -KMASK(i,j) * p25 * WORK3(i,j)

                GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                             + FY(i,j,n) - FY(i,j-1,n)  &
                      + FZTOP(i,j,n,bid) - fz )*dzr(k)*TAREA_R(i,j,bid)

                FZTOP(i,j,n,bid) = fz

              enddo
            enddo

          else

!pw loop split to improve performance

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

               WORK3(i,j) =                                           &
                   ( dz(k) * KAPPA_ISOP(i,j,kbt,k,bid)                &
                   * (  SLX(i,j,ieast ,kbt,k  ,bid)                   &
                       * HYX(i  ,j  ,bid) * TX(i  ,j  ,k  ,n,bid)     &
                      + SLY(i,j,jnorth,kbt,k  ,bid)                   &
                       * HXY(i  ,j  ,bid) * TY(i  ,j  ,k  ,n,bid)     &
                      + SLX(i,j,iwest ,kbt,k  ,bid)                   &
                       * HYX(i-1,j  ,bid) * TX(i-1,j  ,k  ,n,bid)     &
                      + SLY(i,j,jsouth,kbt,k  ,bid)                   &
                       * HXY(i  ,j-1,bid) * TY(i  ,j-1,k  ,n,bid) ) )

              enddo
            enddo

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                WORK3(i,j) = WORK3(i,j)                               &
                    + ( dz_bottom * KAPPA_ISOP(i,j,ktp,kp1,bid)       &
                    * ( SLX(i  ,j  ,ieast ,ktp,kp1,bid)               &
                       * HYX(i  ,j  ,bid) * TX(i  ,j  ,kp1,n,bid)     &
                      + SLY(i  ,j  ,jnorth,ktp,kp1,bid)               &
                       * HXY(i  ,j  ,bid) * TY(i  ,j  ,kp1,n,bid)     &
                      + SLX(i  ,j  ,iwest ,ktp,kp1,bid)               &
                       * HYX(i-1,j  ,bid) * TX(i-1,j  ,kp1,n,bid)     &
                      + SLY(i  ,j  ,jsouth,ktp,kp1,bid)               &
                       * HXY(i  ,j-1,bid) * TY(i  ,j-1,kp1,n,bid) ) )

              enddo
            enddo

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                fz = -KMASK(i,j) * p5 * WORK3(i,j)

                GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                             + FY(i,j,n) - FY(i,j-1,n)  &
                      + FZTOP(i,j,n,bid) - fz )*dzr(k)*TAREA_R(i,j,bid)

                FZTOP(i,j,n,bid) = fz

              enddo
            enddo

          endif

        else                 ! k = km

          do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie

              GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                           + FY(i,j,n) - FY(i,j-1,n)  &
                    + FZTOP(i,j,n,bid) )*dzr(k)*TAREA_R(i,j,bid)

              FZTOP(i,j,n,bid) = c0 

            enddo
          enddo

        endif

!-----------------------------------------------------------------------
!
!     end of tracer loop
!
!-----------------------------------------------------------------------

      enddo

      call timer_stop(timer_nloop, block_id=this_block%local_id)

!-----------------------------------------------------------------------
!
!     diagnostic computation of the bolus velocities 
!
!-----------------------------------------------------------------------
      
      if ( diag_gm_bolus ) then

        do j=1,ny_block-1
          do i=1,nx_block-1

            WORK1(i,j) = (   SF_SLX(i  ,j,ieast,kbt,k,  bid)    &
                  + factor * SF_SLX(i  ,j,ieast,ktp,kp1,bid)    &
                           + SF_SLX(i+1,j,iwest,kbt,k,  bid)    &
                  + factor * SF_SLX(i+1,j,iwest,ktp,kp1,bid))   &
                   * p25 * HYX(i,j,bid)

            WORK2(i,j) = (   SF_SLY(i,j  ,jnorth,kbt,k,  bid)   &
                  + factor * SF_SLY(i,j  ,jnorth,ktp,kp1,bid)   &
                           + SF_SLY(i,j+1,jsouth,kbt,k,  bid)   &
                  + factor * SF_SLY(i,j+1,jsouth,ktp,kp1,bid))  &
                   * p25 * HXY(i,j,bid)

          enddo
        enddo

        UIB = merge( WORK1, c0, k < KMT(:,:,bid)  &
                          .and. k < KMTE(:,:,bid) )
        VIB = merge( WORK2, c0, k < KMT(:,:,bid)  &
                          .and. k < KMTN(:,:,bid) )

        WORK1 = merge( UIT(:,:,bid) - UIB, c0, k <= KMT(:,:,bid)  &
                                         .and. k <= KMTE(:,:,bid) )
        WORK2 = merge( VIT(:,:,bid) - VIB, c0, k <= KMT(:,:,bid)  &
                                         .and. k <= KMTN(:,:,bid) )

        U_ISOP = WORK1 * dzr(k) / HTE(:,:,bid)
        V_ISOP = WORK2 * dzr(k) / HTN(:,:,bid)

        if ( linertial .and. k == 2 ) then
          BOLUS_SP(:,:,bid) = 50.0_r8 * sqrt(U_ISOP**c2 + V_ISOP**c2)  
        endif

        do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie

            if ( k < KMT(i,j,bid) ) then
              WBOT_ISOP(i,j,bid) = WTOP_ISOP(i,j,bid)          &
                                   + TAREA_R(i,j,bid)          &
                               * ( WORK1(i,j) - WORK1(i-1,j)   &
                                 + WORK2(i,j) - WORK2(i,j-1) )
            else
              WBOT_ISOP(i,j,bid) = c0
            endif

          enddo
        enddo

      endif


!-----------------------------------------------------------------------
!
!     update remaining bottom-face fields to top-face fields for next
!     pass
!
!-----------------------------------------------------------------------

      
      if ( diag_gm_bolus ) then
        UIT(:,:,bid) = UIB
        VIT(:,:,bid) = VIB
      endif

!-----------------------------------------------------------------------
!
!     compute isopycnal diffusion cfl diagnostics if required
!
!-----------------------------------------------------------------------

      if (ldiag_cfl) then

        WORK2 = p5 * (KAPPA_ISOP(:,:,ktp,k,bid)  &
                    + KAPPA_ISOP(:,:,kbt,k,bid)) 

        WORK1 = merge(c4*WORK2*(DXTR(:,:,bid)**2 + DYTR(:,:,bid)**2),  &
                      c0, KMT(:,:,bid) > k)
        WORK2 = abs(WORK1)
        call cfl_hdiff (k,bid,WORK2,1,this_block)

      endif

!-----------------------------------------------------------------------
!
!     accumulate time average if necessary; testing is internal to
!       accumulate_tavg_field
!
!-----------------------------------------------------------------------

      if ( mix_pass /= 1 ) then

          call accumulate_tavg_field                      &
                   (p5 * (KAPPA_ISOP(:,:,ktp,k,bid)    &
                       +  KAPPA_ISOP(:,:,kbt,k,bid)),  &
                          tavg_KAPPA_ISOP, bid, k)

          call accumulate_tavg_field                      &
                   (p5 * (KAPPA_THIC(:,:,ktp,k,bid)    &
                       +  KAPPA_THIC(:,:,kbt,k,bid)),  &
                          tavg_KAPPA_THIC, bid, k)

          call accumulate_tavg_field                      &
                   (p5 * (HOR_DIFF(:,:,ktp,k,bid)      &
                       +  HOR_DIFF(:,:,kbt,k,bid)),    &
                          tavg_HOR_DIFF, bid, k)

        if ( transition_layer_on  .and.  k == 1 ) then

            call accumulate_tavg_field (TLT%DIABATIC_DEPTH(:,:,bid),  &
                                        tavg_DIA_DEPTH, bid, 1)  

            call accumulate_tavg_field (TLT%THICKNESS(:,:,bid),       &
                                        tavg_TLT, bid, 1)  

            call accumulate_tavg_field (TLT%INTERIOR_DEPTH(:,:,bid),  &
                                        tavg_INT_DEPTH, bid, 1)  

        endif

        if ( diag_gm_bolus ) then

            call accumulate_tavg_field (U_ISOP, tavg_UISOP, bid, k) 
            call accumulate_tavg_field (V_ISOP, tavg_VISOP, bid, k) 
            call accumulate_tavg_field (WTOP_ISOP(:,:,bid), tavg_WISOP,bid, k)

          if (accumulate_tavg_now(tavg_ADVT_ISOP)) then

            WORK1 = p5 * HTE(:,:,bid) * U_ISOP * ( TMIX(:,:,k,1)  &
                      + eoshift(TMIX(:,:,k,1), dim=1, shift=1) )
            WORK2 = eoshift(WORK1, dim=1, shift=-1)
            WORK3 = WORK1 - WORK2

            WORK1 = p5 * HTN(:,:,bid) * V_ISOP * ( TMIX(:,:,k,1)  &
                      + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )  
            WORK2 = eoshift(WORK1, dim=2, shift=-1)
            WORK3 = WORK3 + WORK1 - WORK2

            WORK1 = c0
            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
            enddo

            call accumulate_tavg_field (WORK1, tavg_ADVT_ISOP, bid, k)

          endif

           if (accumulate_tavg_now(tavg_ADVS_ISOP)) then

            WORK1 = p5 * HTE(:,:,bid) * U_ISOP * ( TMIX(:,:,k,2)  &
                      + eoshift(TMIX(:,:,k,2), dim=1, shift=1) )
            WORK2 = eoshift(WORK1, dim=1, shift=-1)
            WORK3 = WORK1 - WORK2

            WORK1 = p5 * HTN(:,:,bid) * V_ISOP * ( TMIX(:,:,k,2)  &
                      + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
            WORK2 = eoshift(WORK1, dim=2, shift=-1)
            WORK3 = WORK3 + WORK1 - WORK2

            WORK1 = c0
            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
            enddo

            call accumulate_tavg_field (WORK1, tavg_ADVS_ISOP, bid, k)

          endif

          if ( accumulate_tavg_now(tavg_VNT_ISOP)  .or.  &
               accumulate_tavg_now(tavg_VNS_ISOP) ) then

            WORK1 = p5 * V_ISOP * HTN(:,:,bid) * TAREA_R(:,:,bid) 

            if (accumulate_tavg_now(tavg_VNT_ISOP)) then
              WORK2 = WORK1 * (    TMIX(:,:,k,1)  &
                         + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )
              call accumulate_tavg_field (WORK2, tavg_VNT_ISOP, bid, k) 
            endif

            if (accumulate_tavg_now(tavg_VNS_ISOP)) then
              WORK2 = WORK1 * (    TMIX(:,:,k,2)  &
                         + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
              call accumulate_tavg_field (WORK2, tavg_VNS_ISOP, bid, k)
            endif

          endif

        endif ! bolus velocity option on

        do n = 1,nt
          if (accumulate_tavg_now(tavg_HDIFE_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FX(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            call accumulate_tavg_field(WORK1,tavg_HDIFE_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFN_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FY(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            call accumulate_tavg_field(WORK1,tavg_HDIFN_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FZTOP(i,j,n,bid)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            call accumulate_tavg_field(WORK1,tavg_HDIFB_TRACER(n),bid,k)
          endif
        enddo

      endif   ! mix_pass ne 1


!-----------------------------------------------------------------------
!EOC

      end subroutine hdifft_gm

!***********************************************************************
!BOP
! !IROUTINE: kappa_lon_lat_vmhs 
! !INTERFACE:

      subroutine kappa_lon_lat_vmhs (TMIX, UMIX, VMIX, this_block)

! !DESCRIPTION:
!  Variable kappa parameterization by Visbeck et al (1997):
!  \begin{equation}
!     KAPPA_LATERAL = C {{l^2}\over{T}},
!  \end{equation}
!  where $C$ is a constant, $T$ is the (Eady) time scale
!  $f/\surd\overline{Ri}$ and $l$ is the Rossby radius.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX               ! tracers at all vertical levels
                            !  at mixing time level

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX         ! U,V  at all vertical levels
                            !  at mixing time level

      type (block), intent(in) :: &
         this_block         ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         kp1,i,j,kk,k,      &
         k1,k2,             &
         bid                ! local block address for this sub block

      real (r8) :: &
         const, zmin1, zmin2

      real (r8), dimension(nx_block,ny_block) :: &
         UKT, VKT,                               &
         UKPT, VKPT, RNUM,                       &
         WORK1, WORK2, WORK3, WORK4,             &
         TK, TKP,                                & 
         RHOT, RHOS,                             &
         GRATE, LSC                              

!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------

      bid = this_block%local_id

      GRATE = c0
      LSC = c0
      k1 = 0
      k2 = 0

      vert_loop: do k=1,km

!     -2000m < z < -100m

        if ( zt(k) >= 1.0e4_r8 .and. zt(k) <= 2.0e5_r8 ) then

!-----------------------------------------------------------------------
!
!     start computing the Richardson number and the baroclinic wave speed 
!     average velocities at T points of level k
!
!-----------------------------------------------------------------------

          if ( k1 == 0 ) then

            k1 = k              ! k1 is the upper limit for integration

            call ugrid_to_tgrid(UKT,UMIX(:,:,k),bid)
            call ugrid_to_tgrid(VKT,VMIX(:,:,k),bid)

            TK = max(-c2, TMIX(:,:,k,1))

          endif

!-----------------------------------------------------------------------
!
!     compute RZ=Dz(rho) at k=k
!
!-----------------------------------------------------------------------

          if ( k < km ) then
            kp1 = k + 1
          else
            kp1 = km
          endif

          call state (k,kp1,TMIX(:,:,k,1),TMIX(:,:,k,2), &
                      this_block, DRHODT=RHOT, DRHODS=RHOS)

          TKP = max(-c2, TMIX(:,:,kp1,1))

          WORK1 = TK - TKP
          WORK2 = TMIX(:,:,k,2) - TMIX(:,:,kp1,2)

          WORK3 = RHOT * WORK1 + RHOS * WORK2
          WORK3 = min(WORK3,-eps2)

!-----------------------------------------------------------------------
!
!     average velocities at T points of level k+1
!
!-----------------------------------------------------------------------

          call ugrid_to_tgrid(UKPT,UMIX(:,:,kp1),bid)
          call ugrid_to_tgrid(VKPT,VMIX(:,:,kp1),bid)

!-----------------------------------------------------------------------
!
!     local Richardson number at the bottom of T box, level k
!     = -g*Dz(RHO)/RHO_0/(Dz(u)**2+Dz(v)**2)
!
!-----------------------------------------------------------------------

          if (partial_bottom_cells) then

            where (k < KMT(:,:,bid))

!     RHO_0 = 1 in denominator (approximately) in cgs unit.

              WORK4 = p5*(dz(k) + DZT(:,:,k+1,bid))  ! dzw equivalent

              RNUM = -WORK4/((UKT - UKPT)**2 + (VKT - VKPT)**2 + eps)
              GRATE = GRATE + grav*RNUM*WORK4*WORK3
              LSC = LSC - grav*WORK3

            end where

          else ! no partial bottom cells

            where (k < KMT(:,:,bid))

!     RHO_0 = 1 in denominator (approximately) in cgs unit.

              RNUM = -dzw(k)/((UKT - UKPT)**2 + (VKT - VKPT)**2 + eps)
              GRATE = GRATE + grav*RNUM*dzw(k)*WORK3
              LSC = LSC - grav*WORK3

            end where

          endif ! partial bottom cells

!-----------------------------------------------------------------------
!
!     bottom values become top values for next pass
!
!-----------------------------------------------------------------------

          UKT = UKPT
          VKT = VKPT
          TK  = TKP

        else

          if ( k1 /= 0 .and. k2 == 0 ) then
            k2 = k                ! k2 is the lower limit for integration
            exit vert_loop
          endif

        endif                     ! if(zt(k).ge.1.0e4.and.zt(k).lt.2.0e5)

      end do vert_loop

      do j=1,ny_block
        do i=1,nx_block
          kk = KMT(i,j,bid)
          if ( kk > 0 ) then
            zmin1 = min(zt(k1),zt(kk))
            zmin2 = min(zt(k2),zt(kk))
            GRATE(i,j) = GRATE(i,j)/(zmin2-zmin1+eps) ! Ri
            LSC(i,j) = LSC(i,j)*(zmin2-zmin1)         ! c_g^2=N^2*H^2
          else
            GRATE(i,j) = c0
            LSC(i,j)   = c0
          endif
        end do
      end do

!     equatorial inverse of time scale and square of length scale

      WORK1 = sqrt(c2*sqrt(LSC)*BTP(:,:,bid)) ! sqrt(c_g*2*beta)
      WORK2 = sqrt(LSC)/(c2*BTP(:,:,bid))     ! c_g/(2 beta)

!     inverse of time scale

      WORK1 = max(abs(FCORT(:,:,bid)),WORK1)
      GRATE = WORK1/sqrt(GRATE+eps)           ! 1/T = f/sqrt(Ri)

!     square of length scale

      LSC = LSC/(FCORT(:,:,bid)+eps)**2       ! L^2 = c_g^2/f^2
      LSC = min(LSC,WORK2)

!     grid size = lower limit of length scale

      WORK1 = min(DXT(:,:,bid)**2,DYT(:,:,bid)**2)
      LSC = max(LSC,WORK1)

!-----------------------------------------------------------------------
!
!     compute KAPPA_LATERAL
!
!-----------------------------------------------------------------------

!     const = 0.015_r8  ! constant taken from Visbeck et al (1997)
      const = 0.13_r8   
      WORK1 = const*GRATE*LSC                 ! const*L**2/T

!     KAPPA_LATERAL is bounded by 3.0e6 < KAPPA_LATERAL < 4.0e7

      KAPPA_LATERAL(:,:,bid) = min(4.0e7_r8,WORK1)
      KAPPA_LATERAL(:,:,bid) = max(3.0e6_r8,KAPPA_LATERAL(:,:,bid))

      where (KMT(:,:,bid) <= k1)     ! KAPPA_LATERAL=3.0e6 when depth < 100m
        KAPPA_LATERAL(:,:,bid) = 3.0e6_r8
      end where

!-----------------------------------------------------------------------
!EOC

      end subroutine kappa_lon_lat_vmhs 

!***********************************************************************
!BOP
! !IROUTINE: kappa_eg
! !INTERFACE:

      subroutine kappa_eg (TMIX, UMIX, VMIX, this_block)

! !DESCRIPTION:
!  Variable kappa parameterization by Eden and Greatbatch 
!  (2008, Ocean Modelling, v. 20, 223-239):
!  \begin{equation}
!     KAPPA_ISOP = KAPPA_THIC = c {{L^2}*{\sigma}},
!  \end{equation}
!  where $c$ is a tuning parameter (=const_eg), $\sigma$ is the inverse eddy 
!  time scale, and $L$ is the minimum of the Rossby deformation radius and
!  Rhines scale.
!
!  This subroutine returns KAPPA_ISOP in KAPPA_VERTICAL. Here, KAPPA_VERTICAL
!  serves as a temporary array because this subroutine may be called less
!  frequently than every time step and KAPPA_ISOP is modified in each time step.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX               ! tracers at all vertical levels
                            !  at mixing time level

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX         ! U,V at all vertical levels
                            !  at mixing time level

      type (block), intent(in) :: &
         this_block         ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kp1,            &  ! vertical level index
         bid                   ! local block address for this sub block

      real (r8), dimension(nx_block,ny_block) :: &
         TK, TKP,                                &  ! level k and k+1 TMIX
         WORK1, WORK2, WORK3,                    &  ! work arrays
         C_ROSSBY,                               &  ! first baroclinic Rossby wave speed
         L_ROSSBY,                               &  ! Rossby deformation radius
         RHOT, RHOS                                 ! dRHOdT and dRHOdS


      real (r8), dimension(nx_block,ny_block,km) :: &
         SIGMA,         &  ! inverse eddy time scale
         RI_NO             ! local Richardson number 

!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------

      bid = this_block%local_id

      SIGMA = c0
      RI_NO = c0

      BUOY_FREQ_SQ(:,:,:,bid) = c0

!-----------------------------------------------------------------------
!
!     compute buoyancy frequency and Richardson number at the bottom of 
!     T box:
!             Ri = -g*Dz(RHO)/RHO_0/(Dz(u)**2+Dz(v)**2)
!
!     RHO_0 ~ 1 in cgs units.
!
!-----------------------------------------------------------------------

      do k=1,km-1

        if ( k == 1 ) then
          TK = max(-c2, TMIX(:,:,k,1))
        endif

        kp1 = k+1

        call state (k, kp1, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                    this_block, DRHODT=RHOT, DRHODS=RHOS)

        TKP = max(-c2, TMIX(:,:,kp1,1))

        WORK1 = TK - TKP
        WORK2 = TMIX(:,:,k,2) - TMIX(:,:,kp1,2)
        WORK3 = RHOT * WORK1 + RHOS * WORK2
        WORK3 = min(WORK3,-eps2)

        WORK1 =  (UMIX(:,:,k) - UMIX(:,:,kp1))**2  &
               + (VMIX(:,:,k) - VMIX(:,:,kp1))**2

        call ugrid_to_tgrid (WORK2, WORK1, bid)

        where (k < KMT(:,:,bid))
          BUOY_FREQ_SQ(:,:,k,bid) = - grav * WORK3 * dzwr(k)
          WORK1 = dzw(k)**2 / ( WORK2 + eps2 )
          RI_NO(:,:,k) = WORK1 * BUOY_FREQ_SQ(:,:,k,bid) 
        end where

        TK  = TKP

      end do

!-----------------------------------------------------------------------
!
!     compute the first baroclinic gravity-wave phase speed.
!     Computation of Rossby deformation radius follows Chelton et al.(1998)
!
!-----------------------------------------------------------------------

      C_ROSSBY = c0

      k = 1
      where ( k < KMT(:,:,bid) )
        C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k-1)
      endwhere

      do k=1,km
        where ( k < KMT(:,:,bid) )
          C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k)
        endwhere
        where ( k > 1  .and.  k == KMT(:,:,bid) )
          C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k-1,bid)) * dzw(k)
        endwhere
      enddo

      C_ROSSBY = C_ROSSBY / pi

      L_ROSSBY = min( C_ROSSBY / (abs(FCORT(:,:,bid))+eps), &
                      sqrt( C_ROSSBY / (c2*BTP(:,:,bid)) ) )

!-----------------------------------------------------------------------
!
!     compute the inverse time scale
!
!-----------------------------------------------------------------------

      do k=1,km-1
        WORK1 = max( abs( FCORT(:,:,bid) ), sqrt(C_ROSSBY * c2 * BTP(:,:,bid)) )
        where (k < KMT(:,:,bid))
          SIGMA(:,:,k) = SIGMA_TOPO_MASK(:,:,k,bid) * WORK1  &
                        / sqrt( RI_NO(:,:,k) + gamma_eg )
        end where
      enddo

!----------------------------------------------------------------------
!
!     compute KAPPA_ISOP = c_kappa * sigma * min( Rossby_radius, Rhines scales)^2.
!     note that KAPPA_ISOP is stored in KAPPA_VERTICAL. KAPPA_VERTICAL is
!     located at the T-grid points. in the following, KAPPA_ISOP computed at
!     the lower interface of a T-grid box is assigned to the T-grid point 
!     just above this interface for simplicity.
!
!-----------------------------------------------------------------------

      do k=1,km

        WORK1 = min(L_ROSSBY, SIGMA(:,:,k)/BTP(:,:,bid))

        KAPPA_VERTICAL(:,:,k,bid) = const_eg * SIGMA(:,:,k) * WORK1**2 

      enddo

!----------------------------------------------------------------------
!
!     use below-diabatic-layer values of KAPPA_ISOP within the surface 
!     diabatic layer
!
!----------------------------------------------------------------------

      WORK1 = BL_DEPTH(:,:,bid)
      if ( transition_layer_on )  &
        WORK1 = TLT%DIABATIC_DEPTH(:,:,bid)

      do k=km-1,1,-1
        where ( zw(k) <= WORK1 )
          KAPPA_VERTICAL(:,:,k,bid) = KAPPA_VERTICAL(:,:,k+1,bid)
        endwhere
      enddo

!----------------------------------------------------------------------
!
!     impose lower and upper limits on KAPPA
!
!----------------------------------------------------------------------

      KAPPA_VERTICAL(:,:,:,bid) = max(kappa_min_eg, KAPPA_VERTICAL(:,:,:,bid))
      KAPPA_VERTICAL(:,:,:,bid) = min(kappa_max_eg, KAPPA_VERTICAL(:,:,:,bid))

      end subroutine kappa_eg

!***********************************************************************
!BOP
! !IROUTINE: kappa_lon_lat_hdgr 
! !INTERFACE:

      subroutine kappa_lon_lat_hdgr (TMIX, this_block)

! !DESCRIPTION:
!  Variable kappa parameterization based on the vertically averaged
!  horizontal density gradients as a measure of baroclinicity. This
!  approach is from GFDL. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k_min, k_max,  & ! min and max vertical indices for the vertical average
         k,             & ! vertical loop index
         i, j,          & ! horizontal loop indices
         kk,            & ! temporary value of KMT
         bid              ! local block address for this sub block

      real (r8), dimension(nx_block,ny_block) :: &
         KMASKE, KMASKN, &  ! ocean masks for the east and north faces of a T-cell
         WORK,           &  ! work array
         RHOT,           &  ! dRHO/dT
         RHOS,           &  ! dRHO/dS
         RHO_X,          &  ! grid x-direction density difference
         RHO_Y,          &  ! grid y-direction density difference
         RHO_XY_AVG         ! vertically averaged horizontal density gradient

      real (r8), dimension(nx_block,ny_block,2) :: &
         TRACER_X,       &  ! grid x-direction T and S differences
         TRACER_Y           ! grid y-direction T and S differences

      real (r8) :: &
         depth_min, depth_max, & ! actual min and max depths (reflecting model
                                 !  discretization) for computing vertical
                                 !  averages at every horizontal grid point
         depth_w_min,    &       ! protection against k_min = k_max = 0 in
         depth_w_max,    &       !  array zw
         scaling_coefficient     ! multiplicative factor

!-----------------------------------------------------------------------
!
!     parameter values for the formulation 
!
!-----------------------------------------------------------------------

      real (r8), parameter :: &
         tuning_const      = 2.0_r8,    & ! dimensionless tuning constant 
         length_scale_ref  = 50.0e5_r8, & ! constant length scale in cm
         buoyancy_freq_ref = 4.0e-3_r8, & ! constant buoyancy frequency in 1/s
         vert_average_min  = 1.0e4_r8,  & ! min and max depth range for the
         vert_average_max  = 2.0e5_r8,  & !  vertical average in cm
         kappa_min         = 3.0e6_r8,  & ! min KAPPA_LATERAL value in cm^2/s
         kappa_max         = 4.0e7_r8     ! max KAPPA_LATERAL value in cm^2/s

!-----------------------------------------------------------------------
!
!     initialization
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      k_min = -1 
      k_max =  0

      TRACER_X   = c0
      TRACER_Y   = c0
      RHO_X      = c0
      RHO_Y      = c0
      RHO_XY_AVG = c0

      scaling_coefficient = tuning_const * (length_scale_ref**2) * grav &
                           / ( rho_sw * buoyancy_freq_ref ) 

      vert_average_loop: do k=1,km

        if ( zt(k) >= vert_average_min  .and. &
             zt(k) <= vert_average_max ) then

!-----------------------------------------------------------------------
!
!     start computing the vertical average of the horizontal density 
!     gradients 
!
!-----------------------------------------------------------------------

          if ( k_min == -1 )  k_min = k-1 ! k_min is the upper limit for averaging 

          KMASKE = merge(c1, c0, k <= KMT(:,:,bid) .and. &
                         k <= KMTE(:,:,bid))
          KMASKN = merge(c1, c0, k <= KMT(:,:,bid) .and. &
                         k <= KMTN(:,:,bid))

          WORK = max(-c2, TMIX(:,:,k,1))

!-----------------------------------------------------------------------
!
!     compute tracer horizontal differences 
!
!-----------------------------------------------------------------------

          do j=1,ny_block
            do i=1,nx_block-1
              TRACER_X(i,j,1) = KMASKE(i,j) * ( WORK(i+1,j) &
                                              - WORK(i,j) )
              TRACER_X(i,j,2) = KMASKE(i,j)                 &
                          * ( TMIX(i+1,j,k,2) - TMIX(i,j,k,2) )
            enddo
          enddo

          do j=1,ny_block-1
            do i=1,nx_block
              TRACER_Y(i,j,1) = KMASKN(i,j) * ( WORK(i,j+1) &
                                              - WORK(i,j) )
              TRACER_Y(i,j,2) = KMASKN(i,j)                 &
                          * ( TMIX(i,j+1,k,2) - TMIX(i,j,k,2) )
            enddo
          enddo

!-----------------------------------------------------------------------
!
!     now compute horizontal density differences
!
!-----------------------------------------------------------------------

          call state( k, k, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                      this_block, DRHODT=RHOT, DRHODS=RHOS )

          RHO_X = RHOT * TRACER_X(:,:,1) + RHOS * TRACER_X(:,:,2)
          RHO_Y = RHOT * TRACER_Y(:,:,1) + RHOS * TRACER_Y(:,:,2)

!-----------------------------------------------------------------------
!
!     compute horizontal density gradients and perform their vertical 
!     integral 
!
!-----------------------------------------------------------------------

          do j=2,ny_block-1
            do i=2,nx_block-1
              RHO_XY_AVG(i,j) = RHO_XY_AVG(i,j) + dz(k) * sqrt( p5 &
                       * ( ( RHO_X(i,j)**2 + RHO_X(i-1,j)**2 )     &
                          * DXTR(i,j,bid)**2 )                         & 
                       + ( ( RHO_Y(i,j)**2 + RHO_Y(i,j-1)**2 )     &  
                          * DYTR(i,j,bid)**2 ) ) 
            enddo
          enddo

        else

          if ( k_min >= 0  .and.  k_max == 0 ) then
            k_max = k-1            ! k_max is the lower limit for averaging 
            exit vert_average_loop 
          endif

        endif 

      enddo vert_average_loop 

!-----------------------------------------------------------------------
!
!     compute the vertical average of horizontal density gradients for 
!     the specified vertical range
!
!-----------------------------------------------------------------------

      depth_w_min = c0
      depth_w_max = c0
      if ( k_min  > 0 )  depth_w_min = zw(k_min) 
      if ( k_max /= 0 )  depth_w_max = zw(k_max)

      do j=2,ny_block-1
        do i=2,nx_block-1
          kk = KMT(i,j,bid)
          if ( kk > 0 ) then
            depth_min = min( depth_w_min, zw(kk) )
            depth_max = min( depth_w_max, zw(kk) )
            RHO_XY_AVG(i,j) = RHO_XY_AVG(i,j)        &
                             / ( depth_max - depth_min + eps )
          else
            RHO_XY_AVG(i,j) = c0
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!
!     set KAPPA_LATERAL and enforce the limits 
!
!-----------------------------------------------------------------------

      KAPPA_LATERAL(:,:,bid) = scaling_coefficient * RHO_XY_AVG

      KAPPA_LATERAL(:,:,bid) = min( kappa_max, KAPPA_LATERAL(:,:,bid) )
      KAPPA_LATERAL(:,:,bid) = max( kappa_min, KAPPA_LATERAL(:,:,bid) )

      where ( KMT(:,:,bid) <= k_min )        ! KAPPA_LATERAL = kappa_min 
        KAPPA_LATERAL(:,:,bid) = kappa_min   !  when depth < vert_average_min
      endwhere

!-----------------------------------------------------------------------
!EOC

      end subroutine kappa_lon_lat_hdgr

!***********************************************************************
!BOP
! !IROUTINE: kappa_lon_lat_dradius 
! !INTERFACE:

      subroutine kappa_lon_lat_dradius (this_block)

! !DESCRIPTION:
!  Variable kappa parameterization based on
!  \begin{equation}
!      KAPPA_LATERAL = KAPPA_REF {LENGTH_SCALE \over LENGTH_SCALE_REF}
!  \end{equation}
!  where $LENGTH_SCALE=$ min (Rossby deformation radius, sqrt(TAREA)).
!  Computation of Rossby deformation radius follows Chelton et al.(1998)
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k,            &    ! vertical loop index
         bid                ! local block address for this sub block

      real (r8), dimension(nx_block,ny_block) :: &
         WORK1,        &    ! work array
         LENGTH_SCALE       ! mostly used as a temporary array,
                            !  but the final content is a length scale

      real (r8), parameter :: &
         length_scale_ref = 15.0e+5_r8, &  ! reference length scale in cm
         kappa_lateral_min = 0.1e+7_r8, &  ! minimum kappa in cm^2/s
         kappa_lateral_max = 8.0e+7_r8     ! maximum kappa in cm^2/s

!-----------------------------------------------------------------------
!
!     compute the first baroclinic gravity-wave phase speed, considering
!     the full-depth integral of N (consider pi later)
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      LENGTH_SCALE = c0

      k = 1
      where ( k < KMT(:,:,bid) )
        LENGTH_SCALE = LENGTH_SCALE  &
                      + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k-1)
      endwhere

      do k=1,km
        where ( k < KMT(:,:,bid) )
          LENGTH_SCALE = LENGTH_SCALE  &
                        + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k)
        endwhere
        where ( k == KMT(:,:,bid) )
          LENGTH_SCALE = LENGTH_SCALE  &
                        + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * p5 * dz(k)
        endwhere
      enddo

!-----------------------------------------------------------------------
!
!     now compute the Rossby radius of deformation
!
!-----------------------------------------------------------------------

      where ( abs(TLAT(:,:,bid)) > c5 / radian )
        LENGTH_SCALE = LENGTH_SCALE / ( pi * abs(FCORT(:,:,bid)) )
      elsewhere
        LENGTH_SCALE = sqrt( LENGTH_SCALE / ( c2 * pi * BTP(:,:,bid) ) )
      endwhere

!-----------------------------------------------------------------------
!
!     consider grid size as a possible lower limit of length scale.
!     note that if the vertical integral of BUOY_FREQ_SQ is zero, then
!     LENGTH_SCALE = 0. we choose to enforce limits on KAPPA later, rather
!     then on LENGTH_SCALE below.
!
!-----------------------------------------------------------------------

      WORK1 = min( DXT(:,:,bid), DYT(:,:,bid) )
      LENGTH_SCALE = min( LENGTH_SCALE, WORK1 )

!-----------------------------------------------------------------------
!
!     now compute KAPPA_LATERAL
!
!-----------------------------------------------------------------------

      KAPPA_LATERAL(:,:,bid) = ah * LENGTH_SCALE / length_scale_ref

      if ( kappa_isop_type == kappa_type_const )          &
        KAPPA_LATERAL(:,:,bid) = ah_bolus * LENGTH_SCALE  &
                                / length_scale_ref

      KAPPA_LATERAL(:,:,bid) = min(KAPPA_LATERAL(:,:,bid), &
                                   kappa_lateral_max)
      KAPPA_LATERAL(:,:,bid) = max(KAPPA_LATERAL(:,:,bid), &
                                   kappa_lateral_min)

!-----------------------------------------------------------------------
!EOC

      end subroutine kappa_lon_lat_dradius

!***********************************************************************
!BOP
! !IROUTINE: buoyancy_frequency_dependent_profile 
! !INTERFACE:

      subroutine buoyancy_frequency_dependent_profile (TMIX, this_block)

! !DESCRIPTION:
!  Computes normalized buoyancy frequency ($N$) dependent profiles that 
!  are used to represent the vertical variation of the isopycnal and/or
!  thickness diffusion coefficients. We use 
!  \begin{equation}
!   KAPPA_VERTICAL = {{N^2}\over {N_ref^2}}, 
!  \end{equation}
!  where $N_ref$ is the value of $N$ "just" below the surface diabatic
!  layer (SDL) provided that $N^2 > 0$ there. Otherwise, $N_ref$ is the
!  first stable $N^2$ below SDL. KAPPA_VERTICAL is set to unity for shallower
!  depths. Also, 0.1 <= KAPPA_VERTICAL <= 1.0 is enforced. Note that
!  this implies KAPPA_VERTICAL = 0.1 if the local $N^2 \le 0$.
!  KAPPA_VERTICAL is located at the model tracer points.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k,        &          ! vertical loop index
         bid                  ! local block address for this sub block

      integer (int_kind), dimension(nx_block,ny_block) :: &
         K_MIN                ! k index below SDL 

      real (r8), dimension(nx_block,ny_block) :: &
         TEMP_K,           &  ! temperature at level k
         TEMP_KP1,         &  ! temperature at level k+1
         RHOT,             &  ! dRHO/dT
         RHOS,             &  ! dRHO/dS
         BUOY_FREQ_SQ_REF, &  ! reference (normalization) value
                              !  of N^2
         SDL                  ! surface diabatic layer (see below)

      real (r8), dimension(nx_block,ny_block,km) :: &
         BUOY_FREQ_SQ_NORM    ! normalized N^2 defined at level interfaces

!-----------------------------------------------------------------------
!
!     initialize variables
!
!     note that "km+1" value in K_MIN is used as a flag, i.e. it indicates
!     that a minimum value has not been found yet or that it does not
!     exist.
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      BUOY_FREQ_SQ_NORM         = c0
      BUOY_FREQ_SQ_REF          = c0
      KAPPA_VERTICAL(:,:,:,bid) = c1

      K_MIN = merge( km+1, 0, KMT(:,:,bid) /= 0 )  

      SDL = zw(1)
      if ( vmix_itype == vmix_type_kpp )  SDL = KPP_HBLT(:,:,bid) 
      if ( transition_layer_on )    SDL = TLT%INTERIOR_DEPTH(:,:,bid)

      if ( .not. read_n2_data ) then

!-----------------------------------------------------------------------
!
!     compute N^2
!
!-----------------------------------------------------------------------

        do k=1,km-1

          if ( k == 1 ) &
            TEMP_K = max( -c2, TMIX(:,:,k,1) )

          TEMP_KP1 = max( -c2, TMIX(:,:,k+1,1) )

          call state( k, k+1, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                      this_block, DRHODT=RHOT, DRHODS=RHOS )

          where ( k < KMT(:,:,bid) ) 
            BUOY_FREQ_SQ(:,:,k,bid) = max( c0, - grav * dzwr(k) &
                * ( RHOT * ( TEMP_K - TEMP_KP1 )                &
                  + RHOS * ( TMIX(:,:,k,  2) - TMIX(:,:,k+1,2) ) ) )
          endwhere

          TEMP_K = TEMP_KP1

        enddo

      endif

!-----------------------------------------------------------------------
!
!     determine the reference buoyancy frequency and the associated 
!     level index (most likely just below SDL)
!
!-----------------------------------------------------------------------

      do k=1,km-1
        where ( ( K_MIN == km+1 ) .and. ( zw(k) > SDL ) .and.  &
                ( k <= KMT(:,:,bid) )  .and.                   &
                ( BUOY_FREQ_SQ(:,:,k,bid) > c0 ) ) 
          BUOY_FREQ_SQ_REF = BUOY_FREQ_SQ(:,:,k,bid)
          K_MIN = k
        endwhere
      enddo
        
!-----------------------------------------------------------------------
!
!     now compute the normalized profiles at the interfaces 
!
!-----------------------------------------------------------------------

      do k=1,km-1
        where ( ( k >= K_MIN ) .and. ( k < KMT(:,:,bid) ) .and. &
               ( BUOY_FREQ_SQ_REF /= c0 ) )
          BUOY_FREQ_SQ_NORM(:,:,k) =  &
              max( BUOY_FREQ_SQ(:,:,k,bid) / BUOY_FREQ_SQ_REF, 0.1_r8 )
          BUOY_FREQ_SQ_NORM(:,:,k) =  &
              min( BUOY_FREQ_SQ_NORM(:,:,k), c1 ) 
        elsewhere
          BUOY_FREQ_SQ_NORM(:,:,k) = c1
        endwhere
      enddo

      do k=1,km-1
        where ( k == KMT(:,:,bid)-1 ) 
          BUOY_FREQ_SQ_NORM(:,:,k+1) = BUOY_FREQ_SQ_NORM(:,:,k)
        endwhere
      enddo

!-----------------------------------------------------------------------
!
!     finally, do not average interface values of BUOY_FREQ_SQ_NORM to
!     the tracer points, but instead copy from above to preserve the 
!     extrema. 
!
!-----------------------------------------------------------------------

      do k=2,km
        where ( ( k > K_MIN ) .and. ( k <= KMT(:,:,bid) ) )
          KAPPA_VERTICAL(:,:,k,bid) = BUOY_FREQ_SQ_NORM(:,:,k-1)
        endwhere
      enddo

!-----------------------------------------------------------------------
!EOC

      end subroutine buoyancy_frequency_dependent_profile

!***********************************************************************
!BOP
! !IROUTINE: transition_layer 
! !INTERFACE:

      subroutine transition_layer ( this_block )

! !DESCRIPTION:
!  Compute transition layer related fields. the related algorithms
!  should work even with zero transition layer depth.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid                 ! local block address for this sub block

      integer (int_kind), dimension(nx_block,ny_block) :: &
         K_START,   &        ! work arrays for TLT%K_LEVEL and 
         K_SUB               !  TLT%ZTW, respectively

      logical (log_kind), dimension(nx_block,ny_block) :: &
         COMPUTE_TLT         ! flag

      real (r8), dimension(nx_block,ny_block) :: &
         WORK                ! work space for TLT%THICKNESS

      real (r8), dimension(2) :: &
         reference_depth     ! zt or zw

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      K_START = 0
      K_SUB   = 0

      TLT%THICKNESS(:,:,bid)      = c0
      TLT%INTERIOR_DEPTH(:,:,bid) = c0
      TLT%K_LEVEL(:,:,bid)        = 0
      TLT%ZTW(:,:,bid)            = 0

      COMPUTE_TLT = merge(.true., .false., KMT(:,:,bid) /= 0)

!-----------------------------------------------------------------------
!
!     initial pass to determine the minimum transition layer thickness.
!     the vertical level at or below which the interior region starts
!     is also computed.
!
!-----------------------------------------------------------------------

      do k=1,km

        where ( COMPUTE_TLT  .and.  &
                TLT%DIABATIC_DEPTH(:,:,bid) < zw(k) )

          K_START = k+1
          K_SUB   = ktp

          TLT%THICKNESS(:,:,bid) = zw(k) - TLT%DIABATIC_DEPTH(:,:,bid)
          TLT%K_LEVEL(:,:,bid)   = k
          TLT%ZTW(:,:,bid)       = 2

          COMPUTE_TLT = .false.

        endwhere

        where ( k /= 1  .and.  K_START == k+1  .and. &
                TLT%DIABATIC_DEPTH(:,:,bid) < zt(k) )

          K_START = k
          K_SUB   = kbt

          TLT%THICKNESS(:,:,bid) = zt(k) - TLT%DIABATIC_DEPTH(:,:,bid)
          TLT%K_LEVEL(:,:,bid)   = k
          TLT%ZTW(:,:,bid)       = 1

        endwhere

      enddo

#ifdef CCSMCOUPLED
#ifndef _HIRES
      if ( any(COMPUTE_TLT) ) then
        call shr_sys_abort ('Incorrect DIABATIC_DEPTH value in TLT'  &
                        //  ' computation')
      endif
#endif
#endif

!-----------------------------------------------------------------------
!
!     using R * |S| as the vertical displacement length scale associated
!     with a horizontal displacement equivalent to the Rossby deformation
!     radius, R, determine if the transition layer thickness needs to be
!     modified. |S| is the absolute value of the slope of the isopycnal
!     surfaces. First, start with columns where K_SUB = kbt.
!
!-----------------------------------------------------------------------

      where ( KMT(:,:,bid) == 0  .or.  K_START > KMT(:,:,bid)  .or.  &
              ( K_START == KMT(:,:,bid)  .and.  K_SUB == kbt ) )
        COMPUTE_TLT = .false.
      elsewhere
        COMPUTE_TLT = .true.
      endwhere

      do k=1,km-1

        WORK = c0

        where ( COMPUTE_TLT  .and.  K_SUB == kbt  .and.  &
                K_START < KMT(:,:,bid)  .and.  K_START == k )
          WORK = max(SLA_SAVE(:,:,kbt,k,bid), &
                     SLA_SAVE(:,:,ktp,k+1,bid)) * RB(:,:,bid)
        endwhere

        where ( WORK /= c0  .and.  &
                TLT%DIABATIC_DEPTH(:,:,bid) <  (zw(k) - WORK) )
          COMPUTE_TLT = .false.
        endwhere

        where ( WORK /= c0  .and.  &
                TLT%DIABATIC_DEPTH(:,:,bid) >= (zw(k) - WORK) )

          K_START = K_START + 1
          K_SUB   = ktp

          TLT%THICKNESS(:,:,bid) = zw(k) - TLT%DIABATIC_DEPTH(:,:,bid)
          TLT%K_LEVEL(:,:,bid)   = k
          TLT%ZTW(:,:,bid)       = 2

        endwhere

      enddo

!-----------------------------------------------------------------------
!
!     now consider the deeper levels
!
!-----------------------------------------------------------------------

      do k=2,km

        reference_depth(ktp) = zt(k)
        reference_depth(kbt) = zw(k)

        do kk=ktp,kbt

          WORK = c0

          if (kk == ktp) then
            where ( COMPUTE_TLT  .and.  K_START <= KMT(:,:,bid)  .and. &
                    K_START == k )
              WORK = max(SLA_SAVE(:,:,ktp,k,bid), &
                         SLA_SAVE(:,:,kbt,k,bid)) * RB(:,:,bid)
            endwhere
          else
            ! Checking k < km guarantees that k+1 is not out of bounds
            if (k.lt.km) then
              where ( COMPUTE_TLT  .and.  K_START < KMT(:,:,bid)  .and. &
                      K_START == k )
                WORK = max(SLA_SAVE(:,:,kbt,k,bid), &
                           SLA_SAVE(:,:,ktp,k+1,bid)) * RB(:,:,bid)
              endwhere
            end if
            where ( COMPUTE_TLT  .and.  K_START == KMT(:,:,bid)  .and. &
                    K_START == k )
              WORK = SLA_SAVE(:,:,kbt,k,bid) * RB(:,:,bid)
            endwhere
          endif

          where ( WORK /= c0  .and.  &
           TLT%DIABATIC_DEPTH(:,:,bid) <  (reference_depth(kk) - WORK) )
            COMPUTE_TLT = .false.
          endwhere

          where ( WORK /= c0  .and.  &
           TLT%DIABATIC_DEPTH(:,:,bid) >= (reference_depth(kk) - WORK) )
            TLT%THICKNESS(:,:,bid) = reference_depth(kk)  &
                                    - TLT%DIABATIC_DEPTH(:,:,bid)
            TLT%K_LEVEL(:,:,bid)   = k
            TLT%ZTW(:,:,bid)       = kk
          endwhere

        enddo

        where ( COMPUTE_TLT  .and.  K_START == k )
          K_START = K_START + 1
        endwhere

      enddo

#ifdef CCSMCOUPLED
#ifndef _HIRES
      if ( any(COMPUTE_TLT) ) then
        call shr_sys_abort ('Incorrect TLT computations')
      endif
#endif
#endif

!-----------------------------------------------------------------------
!
!     compute the depth at which the interior, adiabatic region starts
!
!-----------------------------------------------------------------------

      do k=1,km
        where ( TLT%K_LEVEL(:,:,bid) == k  .and.  &
                TLT%ZTW(:,:,bid) == 1 )
          TLT%INTERIOR_DEPTH(:,:,bid) = zt(k)
        endwhere
        where ( TLT%K_LEVEL(:,:,bid) == k  .and.  &
                TLT%ZTW(:,:,bid) == 2 )
          TLT%INTERIOR_DEPTH(:,:,bid) = zw(k)
        endwhere
      enddo

      COMPUTE_TLT = .false.
      where ( KMT(:,:,bid) /= 0  .and.  &
              TLT%INTERIOR_DEPTH(:,:,bid) == c0 )  &
        COMPUTE_TLT = .true.
      where ( KMT(:,:,bid) == 0  .and.  &
              TLT%INTERIOR_DEPTH(:,:,bid) /= c0 )  &
        COMPUTE_TLT = .true.

#ifdef CCSMCOUPLED
#ifndef _HIRES
      if ( any(COMPUTE_TLT) ) then
        call shr_sys_abort ('Incorrect TLT%INTERIOR_DEPTH computation')
      endif
#endif
#endif

!-----------------------------------------------------------------------
!EOC

      end subroutine transition_layer

!***********************************************************************
!BOP
! !IROUTINE: merged_streamfunction 
! !INTERFACE:

      subroutine merged_streamfunction ( this_block )

! !DESCRIPTION:
!  Construct a merged streamfunction that has the appropriate
!  behavior in the surface diabatic region, transition layer, and
!  adiabatic interior
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid                 ! local block address for this sub block
 
      real (r8), dimension(nx_block,ny_block,2) :: &
         WORK1, WORK2, WORK3, WORK4   ! work arrays

      real (r8), dimension(nx_block,ny_block) :: &
         WORK2_NEXT, WORK4_NEXT       ! WORK2 or WORK4 at next level

      real (r8), dimension(nx_block,ny_block) :: &
         WORK5, WORK6, WORK7          ! more work arrays

      logical (log_kind), dimension(nx_block,ny_block) :: &
         LMASK               ! flag

      real (r8), dimension(2) :: &
         reference_depth              ! zt or zw

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      SF_SLX(:,:,:,:,:,bid) = c0
      SF_SLY(:,:,:,:,:,bid) = c0

      WORK1 = c0
      WORK2 = c0
      WORK3 = c0
      WORK4 = c0
      WORK5 = c0
      WORK6 = c0
      WORK7 = c0
      WORK2_NEXT = c0
      WORK4_NEXT = c0

!-----------------------------------------------------------------------
!
!     compute the interior streamfunction and its first derivative at the
!     INTERIOR_DEPTH level. WORK1 and WORK2 contain the streamfunction
!     and its first derivative, respectively, for the zonal component
!     for the east and west sides of a grid cell. WORK3 and WORK4 are
!     the corresponding fields for the meridional component for the
!     north and south sides of a grid cell. Note that these definitions
!     include a "dz". Also, the first derivative computations assume 
!     that the streamfunctions are located in the middle of the top or
!     bottom half of a grid cell, hence a factor of two in WORK2 and
!     WORK4 calculations. 
!
!-----------------------------------------------------------------------

      do k=1,km-1

        do kk=1,2

          LMASK = TLT%K_LEVEL(:,:,bid) == k  .and.            &
                  TLT%K_LEVEL(:,:,bid) < KMT(:,:,bid)  .and.  &
                  TLT%ZTW(:,:,bid) == 1

          where ( LMASK ) 

            WORK1(:,:,kk) =  KAPPA_THIC(:,:,kbt,k,bid)  &
                           * SLX(:,:,kk,kbt,k,bid) * dz(k)
            WORK2(:,:,kk) = c2 * dzwr(k) * ( WORK1(:,:,kk)            &
              - KAPPA_THIC(:,:,ktp,k+1,bid) * SLX(:,:,kk,ktp,k+1,bid) &
                                            * dz(k+1) )

            WORK2_NEXT = c2 * ( &
              KAPPA_THIC(:,:,ktp,k+1,bid) * SLX(:,:,kk,ktp,k+1,bid) - &
              KAPPA_THIC(:,:,kbt,k+1,bid) * SLX(:,:,kk,kbt,k+1,bid) )

            WORK3(:,:,kk) =  KAPPA_THIC(:,:,kbt,k,bid)  &
                           * SLY(:,:,kk,kbt,k,bid) * dz(k)
            WORK4(:,:,kk) = c2 * dzwr(k) * ( WORK3(:,:,kk)            &
              - KAPPA_THIC(:,:,ktp,k+1,bid) * SLY(:,:,kk,ktp,k+1,bid) &
                                            * dz(k+1) )

            WORK4_NEXT = c2 * ( &
              KAPPA_THIC(:,:,ktp,k+1,bid) * SLY(:,:,kk,ktp,k+1,bid) - &
              KAPPA_THIC(:,:,kbt,k+1,bid) * SLY(:,:,kk,kbt,k+1,bid) )

          endwhere

          where ( LMASK .and. abs(WORK2_NEXT) < abs(WORK2(:,:,kk)) ) &
            WORK2(:,:,kk) = WORK2_NEXT

          where ( LMASK .and. abs(WORK4_NEXT) < abs(WORK4(:,:,kk)) ) &
            WORK4(:,:,kk) = WORK4_NEXT

          LMASK = TLT%K_LEVEL(:,:,bid) == k  .and.           &
                  TLT%K_LEVEL(:,:,bid) < KMT(:,:,bid)  .and. &
                  TLT%ZTW(:,:,bid) == 2

          where ( LMASK )

            WORK1(:,:,kk) =  KAPPA_THIC(:,:,ktp,k+1,bid)     & 
                           * SLX(:,:,kk,ktp,k+1,bid)
            WORK2(:,:,kk) =  c2 * ( WORK1(:,:,kk)                 &
                           - ( KAPPA_THIC(:,:,kbt,k+1,bid)        &
                              * SLX(:,:,kk,kbt,k+1,bid) ) )
            WORK1(:,:,kk) = WORK1(:,:,kk) * dz(k+1)

            WORK3(:,:,kk) =  KAPPA_THIC(:,:,ktp,k+1,bid)     &
                           * SLY(:,:,kk,ktp,k+1,bid)
            WORK4(:,:,kk) =  c2 * ( WORK3(:,:,kk)                 &
                           - ( KAPPA_THIC(:,:,kbt,k+1,bid)        &
                              * SLY(:,:,kk,kbt,k+1,bid) ) )
            WORK3(:,:,kk) = WORK3(:,:,kk) * dz(k+1)

          endwhere

          LMASK = LMASK .and. TLT%K_LEVEL(:,:,bid) + 1 < KMT(:,:,bid)

          if (k.lt.km-1) then ! added to avoid out of bounds access
            where ( LMASK )

              WORK2_NEXT = c2 * dzwr(k+1) * ( &
                KAPPA_THIC(:,:,kbt,k+1,bid) * SLX(:,:,kk,kbt,k+1,bid) * dz(k+1) - &
                KAPPA_THIC(:,:,ktp,k+2,bid) * SLX(:,:,kk,ktp,k+2,bid) * dz(k+2))

              WORK4_NEXT = c2 * dzwr(k+1) * ( &
                KAPPA_THIC(:,:,kbt,k+1,bid) * SLY(:,:,kk,kbt,k+1,bid) * dz(k+1) - &
                KAPPA_THIC(:,:,ktp,k+2,bid) * SLY(:,:,kk,ktp,k+2,bid) * dz(k+2))

            endwhere
          end if

          where ( LMASK .and. abs(WORK2_NEXT) < abs(WORK2(:,:,kk)) ) &
            WORK2(:,:,kk) = WORK2_NEXT

          where ( LMASK .and. abs(WORK4_NEXT) < abs(WORK4(:,:,kk)) ) &
            WORK4(:,:,kk) = WORK4_NEXT

        enddo

      enddo

!-----------------------------------------------------------------------
!
!     compute the depth independent interpolation factors used in the 
!     linear and quadratic interpolations within the diabatic and 
!     transition regions, respectively.
!
!-----------------------------------------------------------------------

      WORK5 = c0
      where (KMT(:,:,bid) /= 0)
        WORK5(:,:) = c1 / ( c2 * TLT%DIABATIC_DEPTH(:,:,bid) &
                   + TLT%THICKNESS(:,:,bid) )
      endwhere

      WORK6 = c0
      where ((KMT(:,:,bid) /= 0) .AND. (TLT%THICKNESS(:,:,bid) > eps))
        WORK6(:,:) = WORK5(:,:) / TLT%THICKNESS(:,:,bid)
      endwhere

!-----------------------------------------------------------------------
!
!     start of interpolation to construct the merged streamfunction
!
!-----------------------------------------------------------------------

      do k=1,km

        reference_depth(ktp) = zt(k) - p25 * dz(k)
        reference_depth(kbt) = zt(k) + p25 * dz(k)

        do kk=ktp,kbt

!-----------------------------------------------------------------------
!
!     diabatic region: use linear interpolation (in streamfunction) 
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) <= TLT%DIABATIC_DEPTH(:,:,bid)  &
                  .and.  k <= KMT(:,:,bid) ) 

            SF_SLX(:,:,1,kk,k,bid) = reference_depth(kk) * WORK5  &
                  * ( c2 * WORK1(:,:,1) + TLT%THICKNESS(:,:,bid)  &
                     * WORK2(:,:,1) )

            SF_SLX(:,:,2,kk,k,bid) = reference_depth(kk) * WORK5  &
                  * ( c2 * WORK1(:,:,2) + TLT%THICKNESS(:,:,bid)  &
                     * WORK2(:,:,2) )

            SF_SLY(:,:,1,kk,k,bid) = reference_depth(kk) * WORK5  &
                  * ( c2 * WORK3(:,:,1) + TLT%THICKNESS(:,:,bid)  &
                     * WORK4(:,:,1) )

            SF_SLY(:,:,2,kk,k,bid) = reference_depth(kk) * WORK5  &
                  * ( c2 * WORK3(:,:,2) + TLT%THICKNESS(:,:,bid)  &
                     * WORK4(:,:,2) )

          endwhere

!-----------------------------------------------------------------------
!
!     transition layer: use quadratic interpolation (in streamfunction) 
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) > TLT%DIABATIC_DEPTH(:,:,bid)   &
            .and.  reference_depth(kk) <= TLT%INTERIOR_DEPTH(:,:,bid) &
            .and.  k <= KMT(:,:,bid) )

            WORK7 = (TLT%DIABATIC_DEPTH(:,:,bid)  &
                     - reference_depth(kk))**2

            SF_SLX(:,:,1,kk,k,bid) = - WORK7 * WORK6            &
                * ( WORK1(:,:,1) + TLT%INTERIOR_DEPTH(:,:,bid)  &
                   * WORK2(:,:,1) )                             &
               + reference_depth(kk) * WORK5                    &
                * ( c2 * WORK1(:,:,1) + TLT%THICKNESS(:,:,bid)  &
                   * WORK2(:,:,1) )

            SF_SLX(:,:,2,kk,k,bid) = - WORK7 * WORK6            &
                * ( WORK1(:,:,2) + TLT%INTERIOR_DEPTH(:,:,bid)  &
                   * WORK2(:,:,2) )                             &
               + reference_depth(kk) * WORK5                    &
                * ( c2 * WORK1(:,:,2) + TLT%THICKNESS(:,:,bid)  &
                   * WORK2(:,:,2) )

            SF_SLY(:,:,1,kk,k,bid) = - WORK7 * WORK6            &
                * ( WORK3(:,:,1) + TLT%INTERIOR_DEPTH(:,:,bid)  &
                   * WORK4(:,:,1) )                             &
               + reference_depth(kk) * WORK5                    &
                * ( c2 * WORK3(:,:,1) + TLT%THICKNESS(:,:,bid)  &
                   * WORK4(:,:,1) )

            SF_SLY(:,:,2,kk,k,bid) = - WORK7 * WORK6            &
                * ( WORK3(:,:,2) + TLT%INTERIOR_DEPTH(:,:,bid)  &
                   * WORK4(:,:,2) )                             &
               + reference_depth(kk) * WORK5                    &
                * ( c2 * WORK3(:,:,2) + TLT%THICKNESS(:,:,bid)  &
                   * WORK4(:,:,2) )

          endwhere

!-----------------------------------------------------------------------
!
!     interior, adiabatic region: no interpolation is needed. note that
!     "dzw" is introduced here, too, for consistency. 
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) > TLT%INTERIOR_DEPTH(:,:,bid)  & 
                  .and.  k <= KMT(:,:,bid) )

            SF_SLX(:,:,1,kk,k,bid) =  KAPPA_THIC(:,:,kk,k,bid)  &
                              * SLX(:,:,1,kk,k,bid) * dz(k)

            SF_SLX(:,:,2,kk,k,bid) =  KAPPA_THIC(:,:,kk,k,bid)  &
                              * SLX(:,:,2,kk,k,bid) * dz(k)

            SF_SLY(:,:,1,kk,k,bid) =  KAPPA_THIC(:,:,kk,k,bid)  &
                              * SLY(:,:,1,kk,k,bid) * dz(k)

            SF_SLY(:,:,2,kk,k,bid) =  KAPPA_THIC(:,:,kk,k,bid)  &
                              * SLY(:,:,2,kk,k,bid) * dz(k)

          endwhere

        enddo  ! end of kk-loop

      enddo    ! end of k-loop

!-----------------------------------------------------------------------
!EOC

      end subroutine merged_streamfunction

!***********************************************************************
!BOP
! !IROUTINE: apply_vertical_profile_to_isop_hor_diff 
! !INTERFACE:

      subroutine apply_vertical_profile_to_isop_hor_diff ( this_block ) 

! !DESCRIPTION:
!  Apply vertical tapers to KAPPA_ISOP and HOR_DIFF based on their
!  vertical location with respect to the diabatic, transition, and
!  adiabatic regions.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid                 ! local block address for this sub block

      real (r8), dimension(2) :: &
         reference_depth

      bid = this_block%local_id

!-----------------------------------------------------------------------
!
!     start of tapering
!
!-----------------------------------------------------------------------

      do k=1,km

        reference_depth(ktp) = zt(k) - p25 * dz(k)
        reference_depth(kbt) = zt(k) + p25 * dz(k)

        do kk=ktp,kbt

!-----------------------------------------------------------------------
!
!     diabatic region: no isopycnal diffusion 
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) <= TLT%DIABATIC_DEPTH(:,:,bid)  &
                  .and.  k <= KMT(:,:,bid) ) 
            KAPPA_ISOP(:,:,kk,k,bid) = c0
          endwhere

!-----------------------------------------------------------------------
!
!      transition layer: a linear combination of isopcynal and horizontal
!      diffusion coefficients 
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) > TLT%DIABATIC_DEPTH(:,:,bid)   &
            .and.  reference_depth(kk) <= TLT%INTERIOR_DEPTH(:,:,bid) &
            .and.  k <= KMT(:,:,bid)  .and.                           &
                   TLT%THICKNESS(:,:,bid) > eps )

            HOR_DIFF(:,:,kk,k,bid) = ( TLT%INTERIOR_DEPTH(:,:,bid)  &
                  - reference_depth(kk) ) * HOR_DIFF(:,:,kk,k,bid)  &
                  / TLT%THICKNESS(:,:,bid)
            KAPPA_ISOP(:,:,kk,k,bid) = ( reference_depth(kk)        &
                                   - TLT%DIABATIC_DEPTH(:,:,bid) )  &
                 * KAPPA_ISOP(:,:,kk,k,bid) / TLT%THICKNESS(:,:,bid)

          endwhere

!-----------------------------------------------------------------------
!
!     interior region: no horizontal diffusion
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) > TLT%INTERIOR_DEPTH(:,:,bid)  &
                  .and.  k <= KMT(:,:,bid) )
            HOR_DIFF(:,:,kk,k,bid) = c0
          endwhere

        enddo  ! end of kk-loop

      enddo    ! end of k-loop

!-----------------------------------------------------------------------
!EOC

      end subroutine apply_vertical_profile_to_isop_hor_diff

!***********************************************************************

      end module hmix_gm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
