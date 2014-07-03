
module pfixer

!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: pfixer
!
! !DESCRIPTION
! Corrects (or fixes) mass fluxes and edge pressures to be consistent
! with offline surface pressure at next model time step.
!
! !USES
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use abortutils,    only: endrun
  use dynamics_vars, only: T_FVDYCORE_GRID
  use cam_logfile,   only: iulog
  use metdata,       only: met_rlx

! !PUBLIC MEMBER FUNCTIONS
  public :: adjust_press

! !REVISION HISTORY:
!   04.01.30  Stacy Walters    Creation
!   04.02.15  F Vitt  Fixed bug in edge pressure corrections
!   04.08.27  F Vitt  Added ability to handle 2D decomposition
!   05.07.06  Sawyer  Simplified interfaces with grid
!
! EOP
!----------------------------------------------------------------------- 

  private
  real(r8), parameter ::  D0_0                    =   0.0_r8
  real(r8), parameter ::  D0_5                    =   0.5_r8
  real(r8), parameter ::  D1_0                    =   1.0_r8
  real(r8), parameter ::  D100_0                  = 100.0_r8
contains

!-----------------------------------------------------------------------
!       ... adjust mass fluxes and pressures for lin-rood transport
!-----------------------------------------------------------------------

  subroutine adjust_press( grid, ps_mod,  ps_obs, mfx, mfy, pexy )

#if defined( SPMD )
    use mod_comm, only : mp_send3d, mp_recv3d, &
                         mp_sendirr, mp_recvirr
    use parutilitiesmodule, only: parcollective, SUMOP
#endif
    use time_manager, only : get_nstep
!!$    use history,        only: outfld

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    type (T_FVDYCORE_GRID), intent(in) :: grid
    real(r8), intent(in)    :: ps_obs(grid%im,grid%jfirst:grid%jlast)    ! surface pressure at t(n) (Pa)
    real(r8), intent(in)    :: ps_mod(grid%im,grid%jfirst:grid%jlast)    ! surface pressure at t(n) (Pa)
    real(r8), intent(inout) :: mfx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)     ! zonal mass flux
    real(r8), intent(inout) :: mfy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)   ! meridional mass flux
    real(r8), intent(inout) :: pexy(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer, parameter :: nstep0 = -1

    integer  :: i, j, k, km1
    integer  :: nstep
#if defined( SPMD )
    integer  :: dest, src
#endif
    integer  :: ndx(2)
    real(r8) :: dpi(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
    real(r8) :: dpixy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
    real(r8) :: dpi_in(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
    real(r8) :: dpi_inxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,1:grid%km)
    real(r8) :: dpi_c(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
    real(r8) :: dps   (grid%im,grid%jfirst:grid%jlast)
    real(r8) :: dps_in(grid%im,grid%jfirst:grid%jlast)
    real(r8) :: dps_inxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
    real(r8) :: ps_diffxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)

    real(r8) :: dmfx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)     ! original zonnal mass flux
    real(r8) :: dmfy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)   ! original meridional mass flux 
    real(r8) :: emfx(grid%im,grid%jfirst:grid%jlast)     ! zonal mass flux error
    real(r8) :: emfy(grid%im,grid%jfirst:grid%jlast+1)   ! meridional mass flux error
 
    real(r8) :: tmp2d(grid%im,grid%kfirst:grid%klast)
    logical :: debug = .false.
    logical :: method1 = .true.

    integer  :: im, jm, km, ifirstxy, ilastxy, jfirstxy, jlastxy
    integer  :: jfirst, jlast, kfirst, klast

    im       = grid%im
    jm       = grid%jm
    km       = grid%km

    ifirstxy = grid%ifirstxy
    ilastxy  = grid%ilastxy
    jfirstxy = grid%jfirstxy
    jlastxy  = grid%jlastxy

    jfirst   = grid%jfirst
    jlast    = grid%jlast
    kfirst   = grid%kfirst
    klast    = grid%klast
    
#if defined( SPMD )
    ! Send one latitude of mfy to the south
    if( mod(grid%iam,grid%npr_y) /= 0 ) then
       dest = grid%iam-1
    else
       dest = -1
    end if
    if( mod(grid%iam+1,grid%npr_y) /= 0 ) then
       src  = grid%iam+1
    else
       src = -1
    end if
    call mp_send3d( grid%commxy, dest, src, im, jm, km,                  &
                    1, im, jfirst, jlast+1, kfirst, klast,              &
                    1, im, jfirst, jfirst, kfirst, klast, mfy)
#endif

    do j = jfirst,jlast
       dps(:,j) = ps_obs(:,j) - ps_mod(:,j)
    end do

    ! ghost mfy
#if defined( SPMD )
    call mp_recv3d( grid%commxy, src, im, jm, km,                &
                    1, im, jfirst, jlast+1, kfirst, klast,      &
                    1, im, jlast+1, jlast+1, kfirst, klast, mfy)
#endif

    nstep = get_nstep()
!-----------------------------------------------------------------------
!       ... store incoming mass fluxes
!-----------------------------------------------------------------------
    if (debug) then
       do k = kfirst,klast
          do j = jfirst,jlast
             dmfx(:,j,k) = mfx(:,j,k)
             dmfy(:,j,k) = mfy(:,j,k)
          end do
          if( jlast /= jm ) then
             dmfy(:,jlast+1,k) = mfy(:,jlast+1,k)
          end if
       end do
    endif

!-----------------------------------------------------------------------
!       ... incoming mass flux divergence
!-----------------------------------------------------------------------
    call calc_divergence( grid, mfx, mfy, dpi_in )

!-----------------------------------------------------------------------
!       ... surface pressure from mass flux divergence
!-----------------------------------------------------------------------
!  Two different methods to compute change in ps give differnt 
!  results if 2D decomp is used (round off error).  Method 1 gives
!  identical 1D vs 2D decomposition results.
!-----------------------------------------------------------------------
    if (method1) then 

       ! xfer dpi_in to dpi_inxy
       if (grid%twod_decomp .eq. 1) then
#if defined (SPMD)
          call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                           grid%ijk_yz_to_xy%RecvDesc, dpi_in, dpi_inxy,             &
                           modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                           grid%ijk_yz_to_xy%RecvDesc, dpi_in, dpi_inxy,             &
                           modc=grid%modc_dynrun )
#endif
       else            
          dpi_inxy(:,:,:) = dpi_in(:,:,:)
       endif

       ! vertical sum
       do j = jfirstxy,jlastxy
          do i = ifirstxy,ilastxy
             dps_inxy(i,j) = sum( dpi_inxy(i,j,1:km) )
          end do
       end do

       ! xfer dps_inxy to dps_in
       ! Embed in 3D array since transpose machinery cannot handle 2D arrays
       if (grid%twod_decomp .eq. 1) then
#if defined (SPMD)
          do k = 1,km
             do j = jfirstxy,jlastxy
                do i = ifirstxy,ilastxy
                   dpixy(i,j,k) = dps_inxy(i,j)
                enddo
             enddo
          enddo

          call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                  &
                           grid%ijk_xy_to_yz%RecvDesc, dpixy, dpi,                   &
                           modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                  &
                           grid%ijk_xy_to_yz%RecvDesc, dpixy, dpi,                   &
                           modc=grid%modc_dynrun )

          do j = jfirst,jlast
             do i = 1,im
                dps_in(i,j) = dpi(i,j,kfirst)
             enddo
          enddo
#endif
       else            
          dps_in(:,:) = dps_inxy(:,:)
       endif

    else ! method1

       ! this method does not give identical results as the above method
       ! when two dimensional decomposition is used

       do j = jfirst,jlast
          do i = 1,im
             dps_in(i,j) = sum( dpi_in(i,j,kfirst:klast) )
          end do
       end do

#if ( defined SPMD )
       if (grid%twod_decomp .eq. 1) then
          call parcollective( grid%comm_z, SUMOP, im, jlast-jfirst+1,  dps_in )
       endif
#endif

    endif ! method1

!-----------------------------------------------------------------------
!       ... modify (fix) mass fluxes
!-----------------------------------------------------------------------
    call do_press_fix_llnl( grid, dps, dps_in, mfx, mfy )

!-----------------------------------------------------------------------
!       ... modified mass flux divergence
!-----------------------------------------------------------------------
    call calc_divergence( grid, mfx, mfy, dpi_c )
      
!-----------------------------------------------------------------------
!       ... differential mass flux divergence
!-----------------------------------------------------------------------
    do k = kfirst,klast
       do j = jfirst,jlast
          dpi(:,j,k) = dpi_c(:,j,k) - dpi_in(:,j,k)
       end do
    end do

    if (grid%twod_decomp .eq. 1) then
#if defined (SPMD)
       call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                        grid%ijk_yz_to_xy%RecvDesc, dpi, dpixy,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                        grid%ijk_yz_to_xy%RecvDesc, dpi, dpixy,                   &
                        modc=grid%modc_dynrun )
#endif
    else            
       dpixy(:,:,:) = dpi(:,:,:)
    endif


!-----------------------------------------------------------------------
!       ... modify pe
!-----------------------------------------------------------------------

    if (debug) then
       write(iulog,*) ' '
       write(iulog,*) 'adjust_press: max pe diff %  @ nstep,ifirstxy,ilastxy,jfirstxy,jlastxy = ',&
            nstep,ifirstxy,ilastxy,jfirstxy,jlastxy
    endif

    do k = 1+1,km+1
       km1 = k - 1

       if (debug) then
          do j = jfirstxy,jlastxy
             do i = ifirstxy,ilastxy
                ps_diffxy(i,j) = sum( dpixy(i,j,1:km1) )/ pexy(i,k,j )
             end do
          end do
       endif

       if( nstep > nstep0 ) then
          do j = jfirstxy,jlastxy
             do i = ifirstxy,ilastxy
                pexy(i,k,j) = pexy(i,k,j) + sum( dpixy(i,j,1:km1) ) 
             end do
          end do
       end if
       if (debug) then

          ndx(:)       = maxloc( abs( ps_diffxy(:,:) ) )

          ndx(1)       = ndx(1) + ifirstxy - 1
          ndx(2)       = ndx(2) + jfirstxy - 1

          write(iulog,'("pfixer press change error (% error,press adjmnt,new pe)",1x,3i5,1p,3g15.7)') &
               k,ndx(:),D100_0*abs( ps_diffxy(ndx(1),ndx(2)) ), &
               dpixy(ndx(1),ndx(2),km1),pexy(ndx(1),k,ndx(2))

       endif
    end do

    if (debug) then 
       write(iulog,*) ' '
       write(iulog,*) 'adjust_press: max  mass flux error  @ nstep,jfirst,jlast,kfirst,klast = ',&
            nstep,jfirst,jlast,kfirst,klast

       do k = kfirst,klast

          do j=jfirst,jlast
             do i=1,im
                emfx(i,j) = ( mfx(i,j,k)-dmfx(i,j,k) ) 
             enddo
          enddo

          ndx(:)       = maxloc( abs( emfx(:,:) ) )
          ndx(2)       = ndx(2) + jfirst - 1

          write(iulog,'("pfixer max x flux error (diff,fixed,orig) ",1x,3i5,1p,3g15.7)') &
               k,ndx(:), emfx( ndx(1),ndx(2) ) , &
               mfx(ndx(1),ndx(2),k),  dmfx(ndx(1),ndx(2),k)

          do j=jfirst,jlast+1
             do i=1,im
                emfy(i,j) = ( mfy(i,j,k)-dmfy(i,j,k) ) 
             enddo
          enddo

          ndx(:)       = maxloc( abs( emfy(:,:) ) )
          ndx(2)       = ndx(2) + jfirst - 1

          write(iulog,'("pfixer max y flux error (diff,fixed,orig) ",1x,3i5,1p,3g15.7)') &
               k,ndx(:), emfy( ndx(1),ndx(2) ) , &
               mfy(ndx(1),ndx(2),k),  dmfy(ndx(1),ndx(2),k)

       enddo
    endif

  end subroutine adjust_press

!-----------------------------------------------------------------------
!       ... calculate horizontal mass flux divergence
!-----------------------------------------------------------------------
  subroutine calc_divergence( grid, mfx, mfy, dpi )

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    type (T_FVDYCORE_GRID), intent(in) :: grid
    real(r8), intent(in)    :: mfx(grid%im,grid%jfirst:grid%jlast,      &
                                   grid%kfirst:grid%klast)          ! zonal mass flux
    real(r8), intent(in)    :: mfy(grid%im,grid%jfirst:grid%jlast+1,    &
                                   grid%kfirst:grid%klast)        ! meridional mass flux
    real(r8), intent(inout) :: dpi(grid%im,grid%jfirst:grid%jlast,      &
                                   grid%kfirst:grid%klast)          ! horizontal mass flux divergence

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer  :: i, j, k, js2g0, jn2g0
    real(r8) :: sum1
    integer  :: im, jm, km, jfirst, jlast, kfirst, klast

    im    = grid%im
    jm    = grid%jm
    km    = grid%km
    jfirst= grid%jfirst
    jlast = grid%jlast
    kfirst= grid%kfirst
    klast = grid%klast

    js2g0 = max( 2,jfirst )
    jn2g0 = min( jm-1,jlast )

!$omp parallel do private( j, k, sum1 )
    do k = kfirst,klast
!-----------------------------------------------------------------------
!       ... north-south component
!-----------------------------------------------------------------------
       do j = js2g0,jn2g0
          dpi(:,j,k) = (mfy(:,j,k) - mfy(:,j+1,k)) * grid%acosp(j)
       end do
!-----------------------------------------------------------------------
!       ... east-west component
!-----------------------------------------------------------------------
       do j = js2g0,jn2g0
          dpi(:im-1,j,k) = dpi(:im-1,j,k) + mfx(:im-1,j,k) - mfx(2:im,j,k)
          dpi(im,j,k)    = dpi(im,j,k) + mfx(im,j,k) - mfx(1,j,k)
       end do
!-----------------------------------------------------------------------
!       ... poles
!-----------------------------------------------------------------------
       if( jfirst == 1 ) then
          sum1 = -sum( mfy(:,2,k) )*grid%rcap
          dpi(:,1,k) = sum1
       end if
       if( jlast == jm ) then
          sum1 = sum( mfy(:,jm,k) ) * grid%rcap
          dpi(:,jm,k) = sum1
       end if
    end do
!$omp end parallel do

  end subroutine calc_divergence

!-----------------------------------------------------------------------
!       ... fix the mass fluxes to match the met field pressure tendency
!     See: http://asd.llnl.gov/pfix/index.html
!-----------------------------------------------------------------------
  subroutine do_press_fix_llnl( grid, dps, dps_ctm, mfx, mfy )

    use commap,        only : gw => w
#ifdef SPMD
    use mpishorthand,  only : mpicom, mpi_double_precision, mpi_success
    use spmd_dyn,      only : compute_gsfactors
#endif
    use spmd_utils,    only : npes
    use hycoef,        only : hybd
    implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
    type (T_FVDYCORE_GRID), intent(in) :: grid
! surface pressure change from met fields
    real(r8), intent(in)    :: dps(grid%im,grid%jfirst:grid%jlast)
! vert. sum of dpi from original mass fluxes
    real(r8), intent(in)    :: dps_ctm(grid%im,grid%jfirst:grid%jlast)
! zonal mass flux
    real(r8), intent(inout) :: mfx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
! meridional mass flux
    real(r8), intent(inout) :: mfy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer     :: i, j, jglob, k, astat, ierr
    integer     :: jn2g0, js2g0, jn2g1
    integer     :: cnt
#ifdef SPMD
    integer     :: numrecv(0:npes-1)
    integer     :: displs(0:npes-1)
#endif
    real(r8)    :: dpress_g                       ! global pressure error
    real(r8)    :: fxmean, factor
    real(r8)    :: ddps(grid%im,grid%jfirst:grid%jlast)        ! surface pressure change error
    real(r8)    :: dpresslat(grid%jm)
    real(r8)    :: mmfd(grid%jm)
    real(r8)    :: mmf(grid%jm+1)
    real(r8)    :: fxintegral(grid%im+1)
    real(r8)    :: xcolmass_fix(grid%im,grid%jfirst:grid%jlast)

    integer     :: im, jm, km, jfirst, jlast, kfirst, klast

    im    = grid%im
    jm    = grid%jm
    km    = grid%km
    jfirst= grid%jfirst
    jlast = grid%jlast
    kfirst= grid%kfirst
    klast = grid%klast

    js2g0 = max( 2,jfirst )
    jn2g0 = min( jm-1,jlast )
    jn2g1 = min( jm-1,jlast+1 )

    do j = jfirst,jlast
       ddps(:,j) = dps(:,j) - dps_ctm(:,j)
    end do
    factor = D0_5/im
    do j = jfirst,jlast
       dpresslat(j) = sum( ddps(:,j) ) * gw(j) * factor
    end do

#ifdef SPMD
    call compute_gsfactors( 1, cnt, numrecv, displs )
    call mpi_allgatherv( dpresslat(jfirst:jlast), cnt, mpi_double_precision, &
                         dpresslat, numrecv,   displs, mpi_double_precision, mpicom, ierr )
    if( ierr /= mpi_success ) then
       write(iulog,*) 'do_press_fix_llnl: mpi_allgatherv failed; error code = ',ierr
       call endrun
    end if
#endif

    dpress_g = sum( dpresslat(:) )
    if( grid%iam == 0 ) then
       write(iulog,*) 'do_press_fix_llnl: dpress_g = ',dpress_g
    end if

!-----------------------------------------------------------------------
!     calculate mean meridional flux divergence (df/dy).
!     note that mmfd is actually the zonal mean pressure change,
!     which is related to df/dy by geometrical factors.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     	... handle non-pole regions.
!-----------------------------------------------------------------------
    factor = D1_0/im
    do j = jfirst,jlast
       mmfd(j) = dpress_g - sum( ddps(:,j) ) * factor
    end do

#ifdef SPMD
    cnt = jlast - jfirst + 1
    call mpi_allgatherv( mmfd(jfirst:jlast), cnt, mpi_double_precision, &
         mmfd, numrecv, displs, mpi_double_precision, mpicom, ierr )
    if( ierr /= mpi_success ) then
       write(iulog,*) 'do_press_fix_llnl: mpi_allgatherv failed; error code = ',ierr
       call endrun
    end if
#endif

!-----------------------------------------------------------------------
!     calculate mean meridional fluxes (cosp*fy).
!     nb: this calculation is being done for global lats, i.e., (1,jm)
!-----------------------------------------------------------------------
    mmf(2) = mmfd(1) / (grid%rcap*im)
    do j = 2,jm-1
       mmf(j+1) = mmf(j) + mmfd(j) * grid%cosp(j)
    end do

!-----------------------------------------------------------------------
!     fix latitude bands.
!     note that we do not need to worry about geometry here because
!     all boxes in a latitude band are identical.
!     note also that fxintegral(im+1) should equal fxintegral(1),
!     i.e., zero.
!-----------------------------------------------------------------------
!$omp parallel do private( i, j, k, fxmean, fxintegral )
    do j = js2g0,jn2g0
       fxintegral(1) = D0_0
       do i = 1,im
          fxintegral(i+1)  = fxintegral(i) - (ddps(i,j) - dpress_g) - mmfd(j)
       end do
       fxintegral(1)       = fxintegral(im+1)
       fxmean              = sum( fxintegral(:im) ) * factor
       xcolmass_fix(:im,j) = fxintegral(:im) - fxmean
    end do
!$omp end parallel do

!-----------------------------------------------------------------------
!     	... distribute colmass_fix in vertical
!-----------------------------------------------------------------------
!$omp parallel do private( j, k )
    do k = kfirst,klast
       do j = js2g0,jn2g0
          mfx(:,j,k) = mfx(:,j,k) +  met_rlx(k) * xcolmass_fix(:,j) * hybd(k)
       end do
       do j = js2g0,jn2g1
          mfy(:,j,k) = mfy(:,j,k) +  met_rlx(k) * mmf(j) * hybd(k)
       end do
       if( jlast == jm ) then
          mfy(:,jm,k) = mfy(:,jm,k) +  met_rlx(k) * mmf(jm) * hybd(k)
       end if
    end do
!$omp end parallel do

  end subroutine do_press_fix_llnl

end module pfixer
