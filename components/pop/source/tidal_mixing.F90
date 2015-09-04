!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module tidal_mixing

!BOP
! !MODULE: tidal_mixing
!
! !DESCRIPTION:
! This module computes the time-independent part of the tidally 
! driven mixing coefficients. The formulas and the input tidal 
! energy flux are from
!
!    Jayne, S. R., and L. C. St. Laurent, 2001: Parameterizing tidal
!      dissipation over rough topography. Geophys. Res. Lett., 
!      v28, 811-814.
!
!    Simmons, H. L., S. R. Jayne, L. C. St. Laurent, and
!      A. J. Weaver, 2004: Tidally driven mixing in a numerical
!      model of the ocean general circulation. Ocean Modelling,
!      vol 6, 245-263.

! !REVISION HISTORY:
! SVN:$Id$

! !USES

   use kinds_mod
   use domain_size
   use domain
   use blocks
   use io
   use io_types
   use constants
   use exit_mod
   use grid
   use communicate
   use global_reductions
   use broadcast
   use tavg

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
 
   public :: init_tidal_mixing

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public ::  &
      ltidal_mixing        ! if true, tidal mixing is on  

   real (r8), dimension(:,:,:,:), allocatable, public ::  &
      TIDAL_COEF           ! time-independent part of the 
                           ! tidal mixing coefficients

   real (r8),public ::  &
     tidal_mix_max

!EOP
!BOC

!EOC
!***********************************************************************

   contains

!***********************************************************************
!BOP
! !IROUTINE: init_tidal_mixing
! !INTERFACE:

 subroutine init_tidal_mixing

! !DESCRIPTION:
!  Initializes parameters for tidally driven mixing and computes
!  the time independent part of the mixing coefficients
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     input namelist
!
!-----------------------------------------------------------------------

   real (r8) ::              &
      local_mixing_fraction,  &! fraction of energy available for
                               ! mixing local to the generation region
      mixing_efficiency,      &! mixing efficiency
      vertical_decay_scale     ! vertical decay scale for turbulence (cm)

   integer (int_kind) ::     &
      iblock                  ! block index

   character (char_len) ::   &
      tidal_energy_file,     &! input file for reading tidal energy flux
      tidal_energy_file_fmt   ! format (bin or nc) for input file

   namelist /tidal_nml/ ltidal_mixing, local_mixing_fraction,      &
                        mixing_efficiency, vertical_decay_scale,   &
                        tidal_energy_file, &
                        tidal_mix_max, tidal_energy_file_fmt

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      k,                  &! vertical level index
      nml_error,          &! namelist i/o error flag
      nu                   ! i/o unit number

   real (r8) ::  &
      coef 

   real (r8), dimension(:,:,:), allocatable ::  &
      TIDAL_ENERGY_FLUX,    &! input tidal energy flux at T-grid points 
                             ! (W/m^2)
      VERTICAL_FUNC,        &! vertical redistribution function (1/cm)
      WORK1                  ! WORK array

   type (block) ::          &
      this_block             ! block information for current block

   type (datafile) ::       &
      tidal_mixing_file_in   ! io file descriptor

   type (io_dim) :: &
      i_dim,        &! dimension descriptors for horiz dims
      j_dim         

   type (io_field_desc) ::   &
      TIDAL_ENERGY_FLUX_D     ! field descriptors for input field

!-----------------------------------------------------------------------
!
!     set defaults for tidal parameters, then read them from namelist
!
!-----------------------------------------------------------------------

   ltidal_mixing          = .false.
   local_mixing_fraction  = 0.33_r8
   mixing_efficiency      = 0.20_r8
   vertical_decay_scale   = 500.0e02_r8
   tidal_mix_max          = 100.0_r8
   tidal_energy_file      = 'unknown_tidal_energy_file' 
   tidal_energy_file_fmt  = 'bin'
 
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
         read(nml_in, nml=tidal_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar (nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP (SigAbort, 'ERROR reading tidal_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Tidal mixing information'
      write(stdout,blank_fmt)
      write(stdout,*) ' tidal_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,tidal_nml)
      write(stdout,blank_fmt)
      write(stdout,1010) ' ltidal_mixing          = ',  ltidal_mixing
      write(stdout,1020) ' local_mixing_fraction  = ',  local_mixing_fraction
      write(stdout,1020) ' mixing_efficiency      = ',  mixing_efficiency
      write(stdout,1020) ' vertical_decay_scale   = ',  vertical_decay_scale
      write(stdout,1030) ' tidal_energy_file      = ',  tidal_energy_file
 
1010  format (a26,2x,l7)
1020  format (a26,2x,1pe12.5)
1030  format (a26,2x,a80)
   endif


   call broadcast_scalar (ltidal_mixing,          master_task)
   call broadcast_scalar (local_mixing_fraction,  master_task)
   call broadcast_scalar (mixing_efficiency,      master_task)
   call broadcast_scalar (vertical_decay_scale,   master_task)
   call broadcast_scalar (tidal_mix_max,          master_task) 
   call broadcast_scalar (tidal_energy_file,      master_task) 
   call broadcast_scalar (tidal_energy_file_fmt,  master_task) 

!-----------------------------------------------------------------------
!
!  exit if tidal mixing is not enabled
!
!-----------------------------------------------------------------------

   if (.not. ltidal_mixing) return

!-----------------------------------------------------------------------
!
!  if tidal mixing and partial_bottom_cells are both enabled,
!  exit (until pbc option is added to tidal mixing code)
!
!-----------------------------------------------------------------------

   if ( ltidal_mixing .and. partial_bottom_cells)  &
      call exit_POP(SigAbort,  &
          'ERROR: partial bottom cells not implemented with tidal_mixing option')

   allocate ( TIDAL_ENERGY_FLUX(nx_block,ny_block,nblocks_clinic), &
              VERTICAL_FUNC    (nx_block,ny_block,nblocks_clinic), &
              WORK1            (nx_block,ny_block,nblocks_clinic), &
              TIDAL_COEF       (nx_block,ny_block,km,nblocks_clinic))

!-----------------------------------------------------------------------
!
!     read the 2D energy flux array
!
!-----------------------------------------------------------------------

   !***  first create input file 

   tidal_mixing_file_in = construct_file(tidal_energy_file_fmt,             &
                                         full_name=trim(tidal_energy_file), &
                                         record_length=rec_type_dbl,        &
                                         recl_words=nx_global*ny_global)

   !*** open file and read attributes
   call data_set(tidal_mixing_file_in, 'open_read')

   !*** define dimensions
   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)

   !*** define field to be read
   TIDAL_ENERGY_FLUX_D =    &
     construct_io_field('TIDAL_ENERGY_FLUX', dim1=i_dim, dim2=j_dim,  &
                 long_name='Input Tidal Energy Flux at T-grid Points',&
                 units    ='W/m^2', grid_loc ='2110',                 &
                 field_loc = field_loc_center,                        &
                 field_type = field_type_scalar,                      &
                 d2d_array =TIDAL_ENERGY_FLUX(:,:,:))
   call data_set (tidal_mixing_file_in, 'define', TIDAL_ENERGY_FLUX_D)

   !***  read field
   call data_set (tidal_mixing_file_in, 'read', TIDAL_ENERGY_FLUX_D)

   !*** after reading, get rid of io field descriptors and close file
   call destroy_io_field (TIDAL_ENERGY_FLUX_D)
   call data_set (tidal_mixing_file_in, 'close')

   if (my_task == master_task) then
     write (stdout,blank_fmt)
     write (stdout,*) ' file read: ', trim(tidal_energy_file)
   endif

   do iblock = 1,nblocks_clinic

!-----------------------------------------------------------------------
!
!  convert TIDAL_ENERGY_FLUX from W/m^2 to gr/s^3.
!
!-----------------------------------------------------------------------

     TIDAL_ENERGY_FLUX(:,:,iblock) = c1000 * TIDAL_ENERGY_FLUX(:,:,iblock)

!-----------------------------------------------------------------------
!
!  compute the time independent part of the tidal mixing coefficients 
!
!-----------------------------------------------------------------------

     WORK1(:,:,iblock) = c0

     do k=1,km
       where ( k < KMT(:,:,iblock) )
         WORK1(:,:,iblock) = WORK1(:,:,iblock) +    &
             exp(-(HT(:,:,iblock) - zw(k))/vertical_decay_scale) * dzw(k)
       endwhere
     end do ! k

     TIDAL_COEF(:,:,:,iblock) = c0

     coef = local_mixing_fraction * mixing_efficiency / rho_fw

     do k=1,km
       where ( k <= KMT(:,:,iblock) )
         VERTICAL_FUNC(:,:,iblock) =                     &
              exp(-(HT(:,:,iblock) - zw(k))              &
              /vertical_decay_scale) / WORK1(:,:,iblock)

         TIDAL_COEF(:,:,k,iblock) = coef *  &
              TIDAL_ENERGY_FLUX(:,:,iblock) * VERTICAL_FUNC(:,:,iblock)
       elsewhere
         TIDAL_COEF(:,:,k,iblock) = c0
       endwhere
     enddo ! k

   enddo ! iblock

   deallocate ( TIDAL_ENERGY_FLUX,VERTICAL_FUNC,WORK1)



!-----------------------------------------------------------------------
!EOC

 end subroutine init_tidal_mixing

!***********************************************************************

end module tidal_mixing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

