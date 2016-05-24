!-----------------------------------------------------------------------
!BOP
! !ROUTINE: dryairm --- Check dry air mass; set to a predefined value if
!                       nlres is false (initialization run)
!
! !INTERFACE:

subroutine dryairm( grid,  moun,  ps,   tracer,  delp,                   &
                    pe,    nlres_loc )

! !USES:
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use dynamics_vars,      only: T_FVDYCORE_GRID
#if defined( SPMD )
#define CPP_PRT_PREFIX  if( grid%iam == 0 )
#else
#define CPP_PRT_PREFIX
#endif

!fvitt
 use constituents,        only: cnst_type
 use mean_module,         only: gmeanxy
 use cam_control_mod,     only: aqua_planet
 use cam_logfile,         only: iulog
 implicit   none

 type (T_FVDYCORE_GRID), intent(in) :: grid
 logical, intent(in):: nlres_loc
 logical, intent(in):: moun

 real(r8), intent(inout) :: tracer(grid%ifirstxy:grid%ilastxy,                      &
                               grid%jfirstxy:grid%jlastxy,grid%km,grid%ntotq) ! Tracers
 real(r8), intent(inout) :: ps(grid%ifirstxy:grid%ilastxy,                          &
                               grid%jfirstxy:grid%jlastxy)   ! surface pressure
 real(r8), intent(inout) :: delp(grid%ifirstxy:grid%ilastxy,                        &
                                 grid%jfirstxy:grid%jlastxy,grid%km) ! press. thickness
 real(r8), intent(inout) :: pe(grid%ifirstxy:grid%ilastxy,grid%km+1,                &
                               grid%jfirstxy:grid%jlastxy)   ! edge pressure

! !DESCRIPTION:
!  Perform adjustment of the total dry-air-mass while preserving total
!  tracer mass
!  Developer: S.-J. Lin, Aug 2000
!
! !REVISION HISTORY:
!   AAM   01.06.27       Assure agreement thru roundoff for 2D decomp.
!   WS    05.07.06       Simplified interface with grid argument
!   WS    05.08.26       Modified for XY decomposition
!   WS    06.02.21       OMP bug fix (2nd to last DO), removed YZ ver.
!   WS    06.07.01       Transitioned tracers q to T_TRACERS
!
!EOP
!---------------------------------------------------------------------
!BOC

! Use work arrays psdk/psdkg to assure identical answers through roundoff
!    for different z decompositions

      real(r8), allocatable :: psdk(:,:,:)     ! local work array
      real(r8), allocatable :: psdkg(:,:,:)    ! global work array
! dry surface pressure
      real(r8)    psd(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
      real(r8)   drym,drym_loc            ! global mean dry air mass in pascals

      integer :: im, jm, km                            ! Dimensions
      integer :: ifirstxy, ilastxy, jfirstxy, jlastxy  ! XY slice
      integer :: nq                            ! Number of advective tracers         
      real(r8):: ptop

#if defined ( NAVY10 )
      parameter (drym = 98222.0_r8)           ! For US NAVY 10-min terrain
#else
      parameter (drym = 98288.0_r8)           ! For USGS terrain
#endif
      real(r8), parameter ::  D245_0        = 245._r8
      real(r8), parameter ::  D101325_0     = 101325._r8

      integer  i, j, k, ic
      real(r8) psm0, psm1
      real(r8) psdry
      real(r8) dpd

    im       = grid%im
    jm       = grid%jm
    km       = grid%km

    ifirstxy   = grid%ifirstxy
    ilastxy    = grid%ilastxy
    jfirstxy   = grid%jfirstxy
    jlastxy    = grid%jlastxy
    nq         = grid%nq
    ptop       = grid%ptop

    drym_loc = drym
    if (aqua_planet) then
       drym_loc = D101325_0 - D245_0
    end if

! Check global maximum/minimum

    call gmeanxy( grid, ps, psm0 )

    allocate (psdk(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
    allocate (psdkg(ifirstxy:ilastxy,jfirstxy:jlastxy,km))

!$omp  parallel do private(i,j,k)
    do k=1,km
       do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
             psdk(i,j,k) = 0._r8
          enddo
       enddo
    enddo

!$omp  parallel do private(i,j,k)
    do k=1,km
       do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
             psdkg(i,j,k) = 0._r8
          enddo
       enddo
    enddo

!$omp  parallel do private(i,j)
       do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
             psdk(i,j,1) = ptop
          enddo
       enddo

    if( nq .ne. 0 ) then
!$omp  parallel do private(i,j,k)
       do k=1,km
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                psdk(i,j,k) = psdk(i,j,k) +    &
                (1._r8-tracer(i,j,k,1))*(pe(i,k+1,j)-pe(i,k,j))
             enddo
          enddo
       enddo
    else

!$omp  parallel do private(i,j,k)
       do k=1,km
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                psdk(i,j,k) = psdk(i,j,k) +  pe(i,k+1,j) - pe(i,k,j)
             enddo
          enddo
       enddo

    endif

!$omp  parallel do private(i,j,k)
    do k=1,km
       do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
             psdkg(i,j,k) = psdk(i,j,k)
          enddo
       enddo
    enddo

!$omp  parallel do private(i,j)
    do j=jfirstxy,jlastxy
       do i=ifirstxy,ilastxy
          psd(i,j) = 0._r8
       enddo
    enddo

 !$omp  parallel do private(i,j,k)
    do j=jfirstxy,jlastxy
       do k=1,km
          do i=ifirstxy,ilastxy
             psd(i,j) = psd(i,j) + psdkg(i,j,k)
          enddo
       enddo
    enddo

    call gmeanxy( grid, psd, psdry )
 
 CPP_PRT_PREFIX write(iulog,*) 'Total Mass=', 0.01_r8*psm0, '(mb), Dry Mass=', 0.01_r8*psdry, '(mb)'
 CPP_PRT_PREFIX write(iulog,*) 'Total Precipitable Water =', (psm0-psdry)/9.80616_r8, '(kg/m**2)'

    deallocate (psdk)
    deallocate (psdkg)

    if( nlres_loc ) return

    if(moun) then
       dpd = drym_loc - psdry
    else
       dpd = 1000._r8*100._r8 - psdry
    endif
 CPP_PRT_PREFIX write(iulog,*) 'dry mass to be added =', 0.01_r8*dpd

!$omp  parallel do private(i, j, ic)

       do j=jfirstxy,jlastxy

          do ic=1,nq
             do i=ifirstxy,ilastxy
                ! fvitt
                ! don't want to change the initial dry mixing ratios of tracers
                if (cnst_type(ic).ne.'dry') tracer(i,j,km,ic) =        &
                   tracer(i,j,km,ic)*delp(i,j,km)/(delp(i,j,km)+dpd)
             enddo
          enddo

! Adjust the lowest Lagrangian layer
          do i=ifirstxy,ilastxy
             delp(i,j,km) = delp(i,j,km) + dpd
             pe(i,km+1,j) = pe(i,km,j) + delp(i,j,km)
             ps(i,j) = pe(i,km+1,j)
          enddo
       enddo

    call gmeanxy( grid, ps, psm1 )

 CPP_PRT_PREFIX write(iulog,*) 'Total moist surface pressure after adjustment (mb) = ',0.01_r8*psm1 

 return

!EOC
end subroutine dryairm
!---------------------------------------------------------------------
