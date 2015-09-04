!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module niw_mixing

!BOP
! !MODULE: niw_mixing
!
! !DESCRIPTION:
! This module computes the near inertial wave (niw) mixing coefficients. 
! The form of the niw mixing is very similar to that of Jayne et al. and 
! Simmons et al., as in:
!
!    Jayne, S. R., and L. C. St. Laurent, 2001: Parameterizing tidal
!      dissipation over rough topography. Geophys. Res. Lett., 
!      v28, 811-814.
!
!    Simmons, H. L., S. R. Jayne, L. C. St. Laurent, and
!      A. J. Weaver, 2004: Tidally driven mixing in a numerical
!      model of the ocean general circulation. Ocean Modelling,
!      vol 6, 245-263.
!      
! The implementation allows both externally determined but time
! invariant as well as internally model computed energy source
! for niw mixing.
!
! In brief, the niw diffusivity (kn) is set by the expression: 
!
!   kn = (1-f_bl) q Gam En Fn(z)  /  rho N**2  
!
! where it is assumed that surface wind forcing generates a boundary
! layer energy source En (either external time invariant but spatially
! varying as En(x,y), or internal model computed and time varying En(x,y,t)).
! A fraction of this wave energy is absorbed directly in the boundary layer
! (f_bl) and does not affect mixing; the other fraction (1-f_bl) does. Of
! this, a portion is dissipated in the present (x,y) column (q), the rest
! is radiated away as near-inertial waves. Finally, of all the energy
! dissipated in the column, only a fraction Gam actually influences 
! mixing, the rest is dissipated thermally. The mixing is assumed to be
! proportional to En, and distributed vertically by a normalized vertical
! structure function Fn(z) that is a simple exponential decaying from the
! bottom of the boundary layer to ocean bottom with a user specified e-folding
! scale. The denominator ocean density rho and buoyancy frequency N**2 is part 
! of the parameterization as described in the above papers.
!
! This module evaluates the product of all the coefficients just mentioned,
! and then uses whatever En (external or internal) input to evaluate kn for
! a given model state (i.e. N).
!
!   Bruce P. Briegleb   April 2011 
!      

! !REVISION HISTORY:
! SVN:$Id$

! !USES

   use kinds_mod
   use domain_size
   use domain
   use blocks
   use io
   use io_tools
   use io_types
   use constants
   use exit_mod
   use grid
   use communicate
   use global_reductions
   use broadcast
   use tavg
   use time_management

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
 
   public :: init_niw_mixing, niw_mix

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public ::  &
      lniw_mixing             ! namelist variable; if true, niw mixing is on  

   character (char_len), public ::   &
      niw_energy_type         ! namelist variable; type (internal or external) for niw energy source

   real (r8), dimension(:,:,:), allocatable, public ::  &
      NIW_COEF,              &! time-independent part of the niw mixing
      NIW_COS_FACTOR          ! time-independent cos factor

   real (r8), dimension(:,:,:), allocatable, public ::  &
      NIW_ENERGY_FLUX         ! externally specified niw energy flux at T-grid points (W/m^2)

   real (r8),public ::       &
      niw_vert_decay_scale,  &! namelist variable; vertical decay scale for turbulence (cm)
      niw_mix_max             ! namelist variable; maximum diffusivity for niw (cm**2/sec)

   integer (int_kind) ::     &! diagnostic niw history fields
      tavg_KVNIW,            &! tavg id for near-inertial wave tracer diffusivity
      tavg_KVNIW_M,          &! tavg id for near-inertial wave vertical momentum viscosity
      tavg_N2,               &! tavg id for bouyancy frequency squared
      tavg_BFNIW              ! tavg id for bouyancy flux of near-inertial waves

!EOP
!BOC

!EOC
!***********************************************************************

   contains

!***********************************************************************
!BOP
! !IROUTINE: init_niw_mixing
! !INTERFACE:

 subroutine init_niw_mixing

! !DESCRIPTION:
!  Initializes parameters for near inertial wave (niw) mixing and computes 
!  the time independent part of the mixing coefficients. The latter accounts 
!  for various fractions of wave generation energy that is actually used to 
!  set the magnitude of the mixing coefficients.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     input namelist variables (for other public namelist variables, see above)
!
!-----------------------------------------------------------------------

   real (r8) ::                      &
      niw_boundary_layer_absorption, &! fraction of niw energy absorbed in boundary layer
      niw_local_mixing_fraction,     &! niw fraction of energy available for
                                      ! mixing local to the generation region
      niw_mixing_efficiency,         &! niw mixing efficiency (i.e. that portion
                                      ! producing mixing rather than thermal heating)
      niw_obs2model_ratio             ! ratio between observed and modelled NIW strength

   character (char_len) ::           &
      niw_energy_file,               &! input file for reading niw energy flux
      niw_energy_file_fmt             ! format (bin or nc) for input file

   namelist /niw_nml/ lniw_mixing, niw_boundary_layer_absorption,           &
                      niw_local_mixing_fraction, niw_mixing_efficiency,     &
                      niw_obs2model_ratio,                                  &
                      niw_vert_decay_scale,  niw_mix_max,                   &
                      niw_energy_type, niw_energy_file, niw_energy_file_fmt

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

    integer (int_kind) ::   &
      i,j,                  &! indices
      nml_error,            &! namelist i/o error flag
      nu                    ! i/o unit number

   character (char_len) :: &
      warning_string        ! for warning message

   character (char_len) :: &  
      string,              &! for defining history fields
      exit_string           ! exit error message string

   integer (int_kind) ::   &
      iblock                ! block index

   type (block) ::         &
      this_block            ! block information for current block

   type (datafile) ::      &
      niw_mixing_file_in    ! io file descriptor

   type (io_dim) ::        &
      i_dim,               &! dimension descriptors for horiz dims
      j_dim         

   type (io_field_desc) :: &
      NIW_ENERGY_FLUX_D    ! field descriptors for input field

!-----------------------------------------------------------------------
!
!     set defaults for niw parameters, then read them from namelist
!
!-----------------------------------------------------------------------

   lniw_mixing                    = .false.
   niw_boundary_layer_absorption  = 0.70_r8
   niw_local_mixing_fraction      = 0.50_r8
   niw_mixing_efficiency          = 0.20_r8
   niw_obs2model_ratio            = 2.00_r8
   niw_vert_decay_scale           = 2000.0e02_r8
   niw_mix_max                    = 100.0_r8
   niw_energy_type                = 'unknown_niw_energy_type' 
   niw_energy_file                = 'unknown_niw_energy_file' 
   niw_energy_file_fmt            = 'bin'
 
!-----------------------------------------------------------------------
!
!  read namelist input and broadcast variables
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=niw_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar (nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP (SigAbort, 'ERROR reading niw_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' NIW mixing information'
      write(stdout,blank_fmt)
      write(stdout,*) ' niw_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,niw_nml)
      write(stdout,blank_fmt)
      write(stdout,1010) ' lniw_mixing                   = ',  lniw_mixing
      write(stdout,1020) ' niw_boundary layer absorption = ',  niw_boundary_layer_absorption
      write(stdout,1020) ' niw_local_mixing_fraction     = ',  niw_local_mixing_fraction
      write(stdout,1020) ' niw_mixing_efficiency         = ',  niw_mixing_efficiency
      write(stdout,1020) ' niw_obs2model_ratio           = ',  niw_obs2model_ratio
      write(stdout,1020) ' niw_vert_decay_scale          = ',  niw_vert_decay_scale
      write(stdout,1020) ' niw_mix_max                   = ',  niw_mix_max
      write(stdout,1030) ' niw_energy_type               = ',  niw_energy_type
      write(stdout,1030) ' niw_energy_file               = ',  niw_energy_file
      write(stdout,1030) ' niw_energy_file_fmt           = ',  niw_energy_file_fmt
 
1010  format (a26,2x,l7)
1020  format (a26,2x,1pe12.5)
1030  format (a26,2x,a80)
   endif


   call broadcast_scalar (lniw_mixing,                   master_task)
   call broadcast_scalar (niw_local_mixing_fraction,     master_task)
   call broadcast_scalar (niw_mixing_efficiency,         master_task)
   call broadcast_scalar (niw_obs2model_ratio,           master_task)
   call broadcast_scalar (niw_boundary_layer_absorption, master_task)
   call broadcast_scalar (niw_vert_decay_scale,          master_task)
   call broadcast_scalar (niw_mix_max,                   master_task) 
   call broadcast_scalar (niw_energy_type,               master_task) 
   call broadcast_scalar (niw_energy_file,               master_task) 
   call broadcast_scalar (niw_energy_file_fmt,           master_task) 

!-----------------------------------------------------------------------
!
!  exit if niw mixing is not enabled
!
!-----------------------------------------------------------------------

   if (.not. lniw_mixing) return

!-----------------------------------------------------------------------
!
!  test for the existence of a valid niw_energy_file. if not, exit
!
!-----------------------------------------------------------------------
   if (niw_energy_file(1:7) == 'unknown') then
    exit_string = 'FATAL ERROR: lniw_mixing option is active, but niw_energy_file = '//trim(niw_energy_file)
    call document ('init_niw_mixing', exit_string)
    call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  warn user that if niw_energy_type=blke and time_mix_opt=matsuno, the
!  results are ill-defined.
!
!-----------------------------------------------------------------------

   if ( niw_energy_type .eq. 'blke' ) then
     warning_string = &
     'WARNING MESSAGE: niw mixing code with blke energy type ill-defined with matsuno'
     write (stdout,blank_fmt)
     write (stdout,*) warning_string
   endif

!-----------------------------------------------------------------------
!
!  use external energy flux file for niw
!
!-----------------------------------------------------------------------

   if( niw_energy_type .eq. 'external' ) then

     allocate ( NIW_ENERGY_FLUX(nx_block,ny_block,nblocks_clinic), &
                NIW_COEF       (nx_block,ny_block,nblocks_clinic))

!-----------------------------------------------------------------------
!
!     read the 2D energy flux array
!
!-----------------------------------------------------------------------

!***  first create input file 

     niw_mixing_file_in = construct_file(niw_energy_file_fmt,             &
                                       full_name=trim(niw_energy_file), &
                                       record_length=rec_type_dbl,        &
                                       recl_words=nx_global*ny_global)
     !*** open file and read attributes
     call data_set(niw_mixing_file_in, 'open_read')

     !*** define dimensions
     i_dim = construct_io_dim('i', nx_global)
     j_dim = construct_io_dim('j', ny_global)

     !*** define field to be read
     NIW_ENERGY_FLUX_D =    &
       construct_io_field('NIW_ENERGY_FLUX', dim1=i_dim, dim2=j_dim,  &
                   long_name='Input NIW Energy Flux at T-grid Points',&
                   units    ='W/m^2', grid_loc ='2110',               &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d2d_array =NIW_ENERGY_FLUX(:,:,:))
     call data_set (niw_mixing_file_in, 'define', NIW_ENERGY_FLUX_D)

     !***  read field
     call data_set (niw_mixing_file_in, 'read', NIW_ENERGY_FLUX_D)

     !*** after reading, get rid of io field descriptors and close file
     call destroy_io_field (NIW_ENERGY_FLUX_D)
     call data_set (niw_mixing_file_in, 'close')

     if (my_task == master_task) then
       write (stdout,blank_fmt)
       write (stdout,*) ' file read: ', trim(niw_energy_file)
     endif

     do iblock = 1,nblocks_clinic

!-----------------------------------------------------------------------
!
!  convert NIW_ENERGY_FLUX from W/m^2 to g/s^3.
!
!-----------------------------------------------------------------------

       NIW_ENERGY_FLUX(:,:,iblock) = c1000 * NIW_ENERGY_FLUX(:,:,iblock)

!-----------------------------------------------------------------------
!
!  compute the time independent part of the niw energy flux
!
!-----------------------------------------------------------------------

       NIW_COEF(:,:,iblock) = niw_local_mixing_fraction * niw_mixing_efficiency * &
                              niw_obs2model_ratio *                               &
                              (c1 - niw_boundary_layer_absorption) / rho_fw
     enddo ! iblock

   else if ( niw_energy_type .eq. 'blke' ) then

!-----------------------------------------------------------------------
!
!  use internal energy flux file for niw (computed prognostically)
!
!-----------------------------------------------------------------------

     allocate ( NIW_COEF(nx_block,ny_block,nblocks_clinic) )

     do iblock = 1,nblocks_clinic

       NIW_COEF(:,:,iblock) = niw_local_mixing_fraction * niw_mixing_efficiency * &
                              niw_obs2model_ratio *                               &
                              (c1 - niw_boundary_layer_absorption) / rho_fw
     enddo

   else

     call exit_POP (SigAbort, &
     'ERROR: unknown niw_energy_type; only external or blke allowed')

   endif

!-----------------------------------------------------------------------
!
!  pre-compute expensive coefficients
!
!-----------------------------------------------------------------------

     allocate ( NIW_COS_FACTOR(nx_block,ny_block,nblocks_clinic) )
     do iblock = 1,nblocks_clinic
       do j=1,ny_block
         do i=1,nx_block
            NIW_COS_FACTOR(i,j,iblock) =  (cos(pi2*TLATD(i,j,iblock)/10.0_r8) + c1)/c2
         enddo
       enddo
     enddo


!-----------------------------------------------------------------------
!
!  diagnostic niw history fields
!
!-----------------------------------------------------------------------

   string = 'Vertical diabatic diffusivity due to Near Inertial Wave Mixing with Background'
   call define_tavg_field(tavg_KVNIW,'KVNIW',3,               &
                          long_name=trim(string),             &
                          units='centimeter^2/s',             &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )

   string = 'Vertical viscosity due to Near Inertial Wave Mixing with Background'
   call define_tavg_field(tavg_KVNIW_M,'KVNIW_M',3,           &
                          long_name=trim(string),             &
                          units='centimeter^2/s',             &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )

   string = 'Bouyancy frequency squared'
   call define_tavg_field(tavg_N2,'N2',3,                   &
                          long_name=trim(string),           &
                          units='s^-2',                     &
                          grid_loc='3113',                  &
                          coordinates  ='TLONG TLAT z_w_bot time' )

   string = 'Bouyancy flux for near inertial waves'
   call define_tavg_field(tavg_BFNIW,'BFNIW',3,             &
                          long_name=trim(string),           &
                          units='cm^2 s^-3',                &
                          grid_loc='3113',                  &
                          coordinates  ='TLONG TLAT z_w_bot time' )

!-----------------------------------------------------------------------
!EOC

 end subroutine init_niw_mixing

!***********************************************************************

! !IROUTINE: niw_mix
! !INTERFACE:

 subroutine niw_mix(DBLOC, HBLT, KBL, En, VISC, VDC, this_block)

! !DESCRIPTION:
!  Specify the niw mixing coefficients. First, evaluate vertical structure 
!  function below the boundary layer depth, then compute normalization in 
!  the column, and finally evaluate the niw diffusivity contribution. It is 
!  assumed that the viscosity and diffusivity have been already initialized 
!  to background values. This routine must be called after bldepth, so that 
!  the boundary layer depth and first k index below boundary layer are available, 
!  and also after the background initialization.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      DBLOC             ! buoyancy difference between adjacent levels

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      HBLT              ! boundary layer depth

   integer (int_kind), dimension(nx_block,ny_block), intent(in) ::  &
      KBL               ! index of first level below hbl

   real (r8), dimension(nx_block,ny_block), intent(in) ::  &
      En                ! niw energy source flux

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,0:km+1), intent(inout) :: &
      VISC              ! viscosity

   real (r8), dimension(nx_block,ny_block,0:km+1,2), intent(inout) :: &
      VDC               ! diffusivity for tracer diffusion

   type (block), intent(in) :: &
      this_block        ! block information for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block)    ::  &
      WORK1,                &! WORK array
      WORK2                  ! WORK array
   real (r8), dimension(nx_block,ny_block,km) ::  &
      WORK3                  ! WORK array   
   real (r8), dimension(nx_block,ny_block)    ::  &
      WORK4                  ! WORK array
   real (r8), dimension(nx_block,ny_block)    ::  &
      KVNIW,                &! near-inertial wave diffusivity
      KVNIW_M                ! near-inertial wave viscosity
   integer (int_kind) ::    &
      k,                    &! vertical level index
      i,                    &! longitudinal index
      j,                    &! latitudinal index
      bid                    ! local block index
   real (r8) ::             &
      Prandtl                ! Prandtl number
   real (r8), dimension(0:km+1) :: &
      zgrid                  ! depth at cell interfaces

!-----------------------------------------------------------------------
!
!  set constants
!
!-----------------------------------------------------------------------

      Prandtl                = 10.0_r8
      zgrid(0) = eps
      do k=1,km
        zgrid(k) = -zt(k)
      enddo
      zgrid(km+1) = -zw(km)

!-----------------------------------------------------------------------
!
!  check if niw is active
!
!-----------------------------------------------------------------------

   if( lniw_mixing) then

     bid = this_block%local_id

!-----------------------------------------------------------------------
!
!  compute the time dependent niw mixing coefficient
!
!-----------------------------------------------------------------------

     WORK1(:,:)   = c0
     WORK2(:,:)   = c0
     WORK3(:,:,:) = c0
     WORK4(:,:)   = c0

!-----------------------------------------------------------------------
!
!  compute vertical structure function normalization.
!
!-----------------------------------------------------------------------

     do k=1,km
       where ( k >= KBL .and. k < KMT(:,:,bid) )
         WORK1 = WORK1 + &
             exp(-(zw(k) - HBLT)/niw_vert_decay_scale) * dzw(k)
       endwhere
     end do ! k

!-----------------------------------------------------------------------
!
!  compute bouyancy frequency N2 (in WORK2) and combine with niw diffusivity coefficient
!
!-----------------------------------------------------------------------

     do k=1,km
       WORK2 = DBLOC(:,:,k)/(zgrid(k) - zgrid(k+1))
       where (WORK2 > c0)
         WORK3(:,:,k) = En(:,:) /WORK2
       endwhere

      ! k index shifted because N2 is at cell bottom
      ! while output axis is at cell top
      call accumulate_tavg_field(WORK2,tavg_N2,bid,k)
     end do ! k

!-----------------------------------------------------------------------
!
!  combine niw diffusivity coefficient with the evaluated and normalized
!  vertical structure function, then combine with assumed input background 
!  diffusivity; save niw diffusivity just below bottom of mixed layer.
!
!-----------------------------------------------------------------------

     do k=1,km
       where ( k >= KBL .and. k < KMT(:,:,bid) )
         KVNIW = c0
         where ( WORK1 > c0 )
           KVNIW = WORK3(:,:,k) *                               &
              exp(-(zw(k) - HBLT)/niw_vert_decay_scale) / WORK1
         endwhere
       endwhere

       where ( k >= KBL .and. k < KMT(:,:,bid) )

!-----------------------------------------------------------------------
!
!  to avoid double-counting, max value of background and niw mixing taken
!
!-----------------------------------------------------------------------

         KVNIW = max(VDC(:,:,k,1),KVNIW)

!-----------------------------------------------------------------------
!
!  now limit diffusivity to specified input value
!
!-----------------------------------------------------------------------

         KVNIW = min(KVNIW,niw_mix_max)

!-----------------------------------------------------------------------
!
!  now set the viscosity and diffusivities
!
!-----------------------------------------------------------------------

         VISC(:,:,k)  = Prandtl * KVNIW
         VDC(:,:,k,1) = KVNIW
         VDC(:,:,k,2) = KVNIW
         where ( k == KBL ) 
            WORK4 = KVNIW
         endwhere

       endwhere

     enddo ! k

!-----------------------------------------------------------------------
!
!  set diffusivity in boundary layer to be top value of niw mixing; for
!  blmix code when matching interior with boundary layer diffusivities
!
!-----------------------------------------------------------------------

     do k=1,km
       where ( k < KBL )
         VISC(:,:,k)  = Prandtl * WORK4
         VDC(:,:,k,1) = WORK4
         VDC(:,:,k,2) = WORK4
       endwhere

!-----------------------------------------------------------------------
!
!  ensure maximum diffusivity and viscosity in the
!  column are the boundary layer values. 
!
!-----------------------------------------------------------------------

       VDC(:,:,k,1) = min(VDC(:,:,k,1),WORK4)
       VDC(:,:,k,2) = min(VDC(:,:,k,2),WORK4)
       VISC(:,:,k)  = min(VISC(:,:,k),Prandtl*WORK4)

       KVNIW   = VDC(:,:,k,1)
       KVNIW_M = VISC(:,:,k)

!-----------------------------------------------------------------------
!
!  save diffusivity and viscosity to history file
!  save bouyancy flux and fequency also
!
!-----------------------------------------------------------------------

       ! k index shifted because KVNIW and KVNIW_M are at cell bottom
       ! while output axis is at cell top
       call accumulate_tavg_field(KVNIW,tavg_KVNIW,bid,k)
       call accumulate_tavg_field(KVNIW_M,tavg_KVNIW_M,bid,k)

!-----------------------------------------------------------------------
!
!  compute bouyancy flux for near inertial waves and save to history file
!
!-----------------------------------------------------------------------

       if (accumulate_tavg_now(tavg_BFNIW)) then
          ! k index shifted because BFNIW is at cell bottom
          ! while output axis is at cell top
          WORK1 = c0
          WORK2 = DBLOC(:,:,k)/(zgrid(k) - zgrid(k+1))
          where (WORK2 > c0)
            WORK1 = KVNIW*WORK2
          endwhere
          call accumulate_tavg_field(WORK1,tavg_BFNIW,bid,k)
       endif

     enddo ! k

   endif ! if( lniw_mixing )

!-----------------------------------------------------------------------
!EOC

 end subroutine niw_mix

!***********************************************************************

end module niw_mixing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
