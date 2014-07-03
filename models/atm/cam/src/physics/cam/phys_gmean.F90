module phys_gmean
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Perform mixed layer global calculations for energy conservation checks.
!
! Methods: 
! Reproducible (nonscalable): 
!    Gather to a master processor who does all the work.
! Reproducible (scalable): 
!    Convert to fixed point (integer representation) to enable
!    reproducibility when using MPI collectives. Results compared with
!    a nonreproducible (but scalable) algorithm using floating point
!    and MPI_Allreduce to verify the results are good enough.
!
! Author: Byron Boville from SOM code by Jim Rosinski/Bruce Briegleb
! Modified: P. Worley to aggregate calculations (4/04)
! Modified: J. White/P. Worley to introduce scalable algorithms;
!           B. Eaton to remove dycore-specific dependencies and to 
!           introduce gmean_mass (10/07)
! Modified: P. Worley to replace in-place implementation with call
!           to repro_sum.
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use physconst,     only: pi
   use spmd_utils,    only: masterproc
   use ppgrid,        only: pcols, begchunk, endchunk
   use shr_reprosum_mod, only: shr_reprosum_calc, shr_reprosum_tolExceeded, &
                            shr_reprosum_reldiffmax, shr_reprosum_recompute
#if ( defined SPMD )
   use mpishorthand
#endif
   use perf_mod
   use cam_logfile,   only: iulog

   implicit none
   private
   save

   public :: &
      gmean,       &! compute global mean of 2D fields on physics decomposition
      gmean_mass    ! compute global mean mass of constituent fields on physics decomposition

   interface gmean
      module procedure gmean_arr
      module procedure gmean_scl
   endinterface

   CONTAINS

!
!========================================================================
!

   subroutine gmean_arr (arr, arr_gmean, nflds)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of each field in "arr" in the physics 
! chunked decomposition
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in)  :: nflds                 ! number of fields
      real(r8), intent(in) :: arr(pcols,begchunk:endchunk,nflds) 
                                                    ! Input array, chunked
      real(r8), intent(out):: arr_gmean(nflds)      ! global means
!
! Local workspace
!
      real(r8) :: rel_diff(2,nflds)                 ! relative differences between 
                                                    !  'fast' reproducible and 
                                                    !  nonreproducible means
      integer  :: ifld                              ! field index
      logical  :: write_warning
!
!-----------------------------------------------------------------------
!
      call t_startf ('gmean_fixed_repro')
      call gmean_fixed_repro(arr, arr_gmean, rel_diff, nflds)
      call t_stopf ('gmean_fixed_repro')

      ! check that "fast" reproducible sum is accurate enough. If not, calculate
      ! using old method
      write_warning = masterproc
      if ( shr_reprosum_tolExceeded('gmean', nflds, write_warning, &
                                  iulog, rel_diff) ) then
         if ( shr_reprosum_recompute ) then
            do ifld=1,nflds
               if ( rel_diff(1,ifld) > shr_reprosum_reldiffmax ) then
                  call t_startf ('gmean_float_repro')
                  call gmean_float_repro(arr(:,:,ifld), arr_gmean(ifld), 1)
                  call t_stopf ('gmean_float_repro')
               endif
            enddo
         endif
      endif

      return
   end subroutine gmean_arr 

!
!========================================================================
!

   subroutine gmean_scl (arr, gmean)
      use phys_grid, only : get_ncols_p

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of each field in "arr" in the physics 
! chunked decomposition
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      real(r8), intent(in) :: arr(pcols,begchunk:endchunk) 
                                                    ! Input array, chunked
      real(r8), intent(out):: gmean      ! global means
!
! Local workspace
!
      integer, parameter :: nflds = 1
      real(r8) :: gmean_array(nflds)
      real(r8) :: array(pcols,begchunk:endchunk,nflds) 
      integer  :: ncols, lchnk

      do lchnk=begchunk,endchunk
         ncols = get_ncols_p(lchnk)
         array(:ncols,lchnk,1) = arr(:ncols,lchnk)
      enddo
      call gmean_arr(array,gmean_array,nflds)
      gmean = gmean_array(1)

   end subroutine gmean_scl 

!
!========================================================================
!

   subroutine gmean_mass(title, state)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the 
! physics decomposition. Prints diagnostics to log file.
!
! Author: B. Eaton (based on gavglook)
!
!-----------------------------------------------------------------------
      use ppgrid,         only: pver
      use physconst,      only: gravit
      use phys_grid,      only: get_ncols_p
      use physics_types,  only: physics_state
      use constituents,   only: pcnst, cnst_name
!
! Arguments
!
      character(len=*),    intent(in) :: title    ! location of this call
      type(physics_state), intent(in) :: state(begchunk:endchunk)
!
! Local workspace
!
      character(len=*), parameter :: sub_name='gmean_mass: '

      integer :: c, i, k, m
      integer :: ierr
      integer :: ncols

      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
      real(r8) :: mass_wet_mean(pcnst)     ! global mean constituent masses assuming moist mmr
      real(r8) :: mass_dry_mean(pcnst)     ! global mean constituent masses assuming dry mmr
      real(r8) :: mmr_max(pcnst)           ! maximum constituent mmr in this process
      real(r8) :: mmr_min(pcnst)           ! minimum constituent mmr in this process
      real(r8) :: mmr_max_glob(pcnst)      ! global maximum constituent mmr
      real(r8) :: mmr_min_glob(pcnst)      ! global minimum constituent mmr
!
!-----------------------------------------------------------------------
!
      allocate(mass_wet(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      mmr_max(:) = -1.e36_r8
      mmr_min(:) =  1.e36_r8
      do m = 1, pcnst
         do c = begchunk, endchunk
            ncols = get_ncols_p(c)
            do i = 1, ncols

               ! Compute column masses assuming both dry and wet mixing ratios

               mass_wet(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
                                    state(c)%pdel(i,k)*state(c)%q(i,k,m)
                  mmr_max(m) = max(mmr_max(m), state(c)%q(i,k,m))
                  mmr_min(m) = min(mmr_min(m), state(c)%q(i,k,m))
               end do
               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

               mass_dry(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
                                    state(c)%pdeldry(i,k)*state(c)%q(i,k,m)
               end do
               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

            end do
         end do
      end do

      ! compute global mean mass
      call gmean(mass_wet, mass_wet_mean, pcnst)
      call gmean(mass_dry, mass_dry_mean, pcnst)

      ! global min/max mmr
#ifdef SPMD
      call mpi_reduce(mmr_max, mmr_max_glob, pcnst, mpir8, MPI_MAX, 0, mpicom, ierr)
      call mpi_reduce(mmr_min, mmr_min_glob, pcnst, mpir8, MPI_MIN, 0, mpicom, ierr)
#else
      mmr_max_glob(:) = mmr_max(:)
      mmr_min_glob(:) = mmr_min(:)
#endif

      ! report to log file
      if (masterproc) then

         do m = 1, pcnst
               write (6,66) trim(title)//' m=',m, &
                  'name='//trim(cnst_name(m))//' gavg dry, wet, min, max ', &
                  mass_dry_mean(m), mass_wet_mean(m), mmr_min_glob(m), mmr_max_glob(m)
66             format (a24,i2,a36,1p,4e25.13)
         end do

      endif

      deallocate(mass_wet)
      deallocate(mass_dry)

   end subroutine gmean_mass

!
!========================================================================
!

   subroutine gmean_float_repro (arr, arr_gmean, nflds)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of each field in "arr" in the physics 
! chunked decomposition - all work is done on the masterproc to avoid
! order of operations differences and assure bfb reproducibility.
! 
!-----------------------------------------------------------------------

      use rgrid,        only: nlon
      use dycore,       only: dycore_is
      use phys_grid,    only: gather_chunk_to_field
      use dyn_grid,     only: get_horiz_grid_dim_d, get_horiz_grid_d, get_dyn_grid_parm_real1d
!
! Arguments
!
      integer, intent(in)  :: nflds               ! number of fields
      real(r8), intent(in) :: &
         arr(pcols,begchunk:endchunk,nflds)       ! Input array, chunked
      real(r8), intent(out):: arr_gmean(nflds)    ! global means
!
! Local workspace
!
      real(r8), pointer :: w(:)
      real(r8) :: zmean                       ! zonal mean value
      real(r8) :: tmean                       ! temp global mean value
      integer :: i, j, ifld, n                ! longitude, latitude, field, 
                                              !  and global column indices
      integer :: hdim1, hdim2                 ! dimensions of rectangular horizontal 
                                              !  grid data structure, If 1D data 
                                              !  structure, then hdim2_d == 1.
      integer :: ngcols                       ! global column count (all)

      ! rectangular version of arr
      real(r8), allocatable :: arr_field(:,:,:)   

      ! column integration weight (from dynamics)
      real(r8), dimension(:), allocatable :: wght_d 

!
!-----------------------------------------------------------------------
!
      call get_horiz_grid_dim_d(hdim1, hdim2)
      allocate(arr_field(hdim1,hdim2,nflds))

      arr_field(:,:,:) = 0.0_r8
      call gather_chunk_to_field (1, 1, nflds, hdim1, arr, arr_field)

      if (masterproc) then

         if (dycore_is('UNSTRUCTURED')) then

            ngcols = hdim1*hdim2
            allocate ( wght_d(1:ngcols) )

            wght_d = 0.0_r8
            call get_horiz_grid_d(ngcols, wght_d_out=wght_d)

            do ifld=1,nflds
               arr_gmean(ifld) = 0._r8
               do j=1,hdim2
                  do i=1,hdim1
                     n = (j-1)*hdim1 + i
                     arr_gmean(ifld) = arr_gmean(ifld) + &
                                       arr_field(i,j,ifld)*wght_d(n)
                  end do
               end do
               arr_gmean(ifld) = arr_gmean(ifld) / (4.0_r8 * pi)
            end do

            deallocate ( wght_d )

         else
            w => get_dyn_grid_parm_real1d('w')
            do ifld=1,nflds
               tmean = 0._r8
               do j=1,hdim2
                  zmean = 0._r8
                  do i=1,hdim1
                     zmean = zmean + arr_field(i,j,ifld)
                  end do
                  tmean = tmean + zmean * 0.5_r8*w(j)/nlon(j)
               end do
               arr_gmean(ifld) = tmean
            end do

         end if

      end if

#if ( defined SPMD )
      call mpibcast (arr_gmean, nflds, mpir8, 0, mpicom)
#endif
      deallocate(arr_field)

      return

   end subroutine gmean_float_repro

!
!========================================================================
!
   subroutine gmean_fixed_repro (arr, arr_gmean, rel_diff, nflds)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of each field in "arr" in the physics 
! chunked decomposition with a reproducible yet scalable implementation
! based on a fixed-point algorithm.
!
!-----------------------------------------------------------------------
      use phys_grid, only    : get_ncols_p, get_wght_all_p, ngcols_p, &
                               get_nlcols_p
!
! Arguments
!
      integer, intent(in)  :: nflds             ! number of fields
      real(r8), intent(in) :: &
         arr(pcols,begchunk:endchunk,nflds)     ! Input array, chunked
      real(r8), intent(out):: arr_gmean(nflds)  ! global means
      real(r8), intent(out):: rel_diff(2,nflds) ! relative and absolute
                                                !  differences between 
                                                !  reproducible and nonreproducible
                                                !  means
!
! Local workspace
!
      integer :: lchnk, i, ifld                 ! chunk, column, field indices
      integer :: ncols                          ! number of columns in current chunk
      integer :: count                          ! summand count
      integer :: ierr                           ! MPI error return
#if (! defined SPMD)
      integer  :: mpicom = 0
#endif
      
      real(r8) :: wght(pcols)                   ! column for integration weights
      real(r8), allocatable :: xfld(:,:)        ! weighted summands
      integer :: nlcols
!
!-----------------------------------------------------------------------
!
      nlcols = get_nlcols_p()
      allocate(xfld(nlcols, nflds))

! pre-weight summands
      do ifld=1,nflds
         count = 0
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            call get_wght_all_p(lchnk, ncols, wght)
            do i=1,ncols
               count = count + 1
               xfld(count,ifld) = arr(i,lchnk,ifld)*wght(i)
            end do
         end do
      end do

! call fixed-point algorithm
      call shr_reprosum_calc (xfld, arr_gmean, count, nlcols, nflds, &
                      gbl_count=ngcols_p, commid=mpicom, rel_diff=rel_diff) 

      deallocate(xfld)
! final normalization
      arr_gmean(:) = arr_gmean(:) / (4.0_r8 * pi)

      return

   end subroutine gmean_fixed_repro

end module phys_gmean
