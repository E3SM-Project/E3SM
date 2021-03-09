!-----------------------------------------------------------------------
!BOP
! !IROUTINE: p_d_adjust --- complete full physics update
!
! !INTERFACE: 
  subroutine p_d_adjust(grid, tracer, pelnxy, pkxy, pkzxy, zvir, &
                        cappa, delpxy, ptxy, pexy, psxy, ptop)

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use dynamics_vars, only : T_FVDYCORE_GRID
#if defined( SPMD )
    use parutilitiesmodule, only: parcollective, sumop
#endif
    use shr_reprosum_mod, only : shr_reprosum_calc, shr_reprosum_tolExceeded, &
                              shr_reprosum_reldiffmax, &
                              shr_reprosum_recompute
    use cam_logfile,   only : iulog
    use perf_mod
!-----------------------------------------------------------------------
    implicit none

! !INPUT PARAMETERS:
    type (T_FVDYCORE_GRID), intent(in) :: grid
    real(r8), intent(in) :: zvir
    real(r8), intent(in) :: ptop
    real(r8), intent(in) :: cappa

! !INPUT/OUTPUT PARAMETERS:
    real(r8), intent(inout) :: tracer(grid%ifirstxy:grid%ilastxy,                     &
                                      grid%jfirstxy:grid%jlastxy,grid%km,grid%ntotq) ! constituents
    real(r8), intent(inout) :: delpxy(grid%ifirstxy:grid%ilastxy,                     &
                                      grid%jfirstxy:grid%jlastxy,grid%km) ! pressure difference
    real(r8), intent(inout) :: ptxy (grid%ifirstxy:grid%ilastxy,                      &
                                     grid%jfirstxy:grid%jlastxy, grid%km) ! Virtual pot temp
    real(r8), intent(inout) :: pexy(grid%ifirstxy:grid%ilastxy,                       &
                                    grid%km+1,grid%jfirstxy:grid%jlastxy)

! !OUTPUT PARAMETERS
    real(r8), intent(out) :: psxy(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy)  ! surf. press
    real(r8), intent(out) :: pelnxy(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy) ! interface pres
    real(r8), intent(out) :: pkzxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy, grid%km)    ! Layer-mean value of PK
    real(r8), intent(out) :: pkxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy, grid%km+1)  ! PE**cappa

! !DESCRIPTION:
!
!   Complete adjustment of quantities after physics update
!
! !REVISION HISTORY:
!   00.06.01   Grant?     Creation
!   01.06.08   AAM        Created from p_d_coupling
!   02.04.24   WS         New mod_comm interface
!   02.05.01   WS         Fix of S.-J. and Phil to peln, pk update
!   03.03.31   BAB        dry mass adjustment moved to dme_adjust, just finish up here
!   05.07.06   WS         Use grid argument to get all grid-related data
!   05.09.23   WS         Transitioned to XY variables only
!   06.07.01   WS         Transitioned tracers q3 to T_TRACERS
!   08.06.25   PW         Added call to fixed point reproducible sum
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
    real(r8), parameter ::  D0_0                    =  0.0_r8
    real(r8), parameter ::  D1_0                    =  1.0_r8

    real(r8) :: pole(grid%ifirstxy:grid%ilastxy,grid%km,grid%ntotq+2)
                                               !  Array containing local pole values
    real(r8) :: pole_sum(grid%km,grid%ntotq+2) !  Array containing average of all pole values
    real(r8) :: rel_diff(2,grid%km,grid%ntotq+2)
    real(r8),allocatable :: pole_tmp(:)

    integer :: i, k, m, j ! indices
    integer :: im, jm, km, ntotq, lim
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy

    logical  :: write_warning

!---------------------------End Local workspace-------------------------

!
! ----------------------------------------------------
! Complete update of dynamics variables
! ----------------------------------------------------
!

    im       = grid%im
    jm       = grid%jm
    km       = grid%km
    ntotq    = grid%ntotq

    ifirstxy = grid%ifirstxy
    ilastxy  = grid%ilastxy
    jfirstxy = grid%jfirstxy
    jlastxy  = grid%jlastxy

    lim      = (ilastxy-ifirstxy) + 1

    ! Average the pole values (WS 2006/02/16, bug fix)

    if (jfirstxy==1) then

       !$omp parallel do private(i,k,m)
       do k = 1, km
          do i = ifirstxy, ilastxy
             pole(i,k,1) = delpxy(i,1,k)
             pole(i,k,2) =   ptxy(i,1,k)
          enddo
          do m = 1, ntotq
             do i = ifirstxy, ilastxy
                pole(i,k,m+2) = tracer(i,1,k,m)
             enddo
          enddo
       enddo

       call t_startf("pdadj_reprosum")
       call shr_reprosum_calc(pole, pole_sum, lim, lim, km*(ntotq+2), gbl_count=im, &
                      commid=grid%commxy_x, rel_diff=rel_diff) ! South pole
       call t_stopf("pdadj_reprosum")

       ! check that "fast" reproducible sum is accurate enough. If not, calculate
       ! using old method
       write_warning = .false.
       if (grid%myidxy_x == 0) write_warning = .true.
       if ( shr_reprosum_tolExceeded('p_d_adjust/South Pole', km*(ntotq+2), &
                                   write_warning, iulog, rel_diff) ) then
          if ( shr_reprosum_recompute ) then
             call t_startf("pdadj_sumfix")
             allocate( pole_tmp(im) )
             do m = 1, ntotq+2
                do k = 1, km
                   if (rel_diff(1,k,m) > shr_reprosum_reldiffmax) then
                      pole_tmp(:) = D0_0
                      do i = ifirstxy, ilastxy
                         pole_tmp(i) = pole(i,k,m)
                      enddo
#if defined(SPMD)
                      call parcollective(grid%commxy_x,sumop,im,pole_tmp)
#endif
                      pole_sum(k,m) = D0_0
                      do i = 1, im
                         pole_sum(k,m) = pole_sum(k,m) + pole_tmp(i)
                      enddo
                   endif
                enddo
             enddo
             deallocate( pole_tmp )
             call t_stopf("pdadj_sumfix")
          endif
       endif

       ! save results
       !$omp parallel do private(i,k,m)
       do k = 1, km
          ! normalize first
          do m = 1, ntotq+2
             pole_sum(k,m) = pole_sum(k,m)/im
          enddo
          do i = ifirstxy,ilastxy
             delpxy(i,1,k) = pole_sum(k,1)
             ptxy(i,1,k)   = pole_sum(k,2)
          enddo
          do m = 1, ntotq
             do i = ifirstxy,ilastxy
                tracer(i,1,k,m) = pole_sum(k,m+2)
             enddo
          enddo
       enddo

    endif ! jfirstxy==1

    if (jlastxy==jm) then

       !$omp parallel do private(i,k,m)
       do k = 1, km
          do i = ifirstxy, ilastxy
             pole(i,k,1) = delpxy(i,jm,k)
             pole(i,k,2) =   ptxy(i,jm,k)
          enddo
          do m = 1, ntotq
             do i = ifirstxy, ilastxy
                pole(i,k,m+2) = tracer(i,jm,k,m)
             enddo
          enddo
       enddo

       call t_startf("pdadj_reprosum")
       call shr_reprosum_calc(pole, pole_sum, lim, lim, km*(ntotq+2), gbl_count=im, &
                      commid=grid%commxy_x, rel_diff=rel_diff) ! North pole
       call t_stopf("pdadj_reprosum")

       ! check that "fast" reproducible sum is accurate enough. If not, calculate
       ! using old method
       write_warning = .false.
       if (grid%myidxy_x == 0) write_warning = .true.
       if ( shr_reprosum_tolExceeded('p_d_adjust/Nouth Pole', km*(ntotq+2), &
                                   write_warning, iulog, rel_diff) ) then
          if ( shr_reprosum_recompute ) then
             call t_startf("pdadj_sumfix")
             allocate( pole_tmp(im) )
             do m = 1, ntotq+2
                do k = 1, km
                   if (rel_diff(1,k,m) > shr_reprosum_reldiffmax) then
                      pole_tmp(:) = D0_0
                      do i = ifirstxy, ilastxy
                         pole_tmp(i) = pole(i,k,m)
                      enddo
#if defined(SPMD)
                      call parcollective(grid%commxy_x,sumop,im,pole_tmp)
#endif
                      pole_sum(k,m) = D0_0
                      do i = 1, im
                         pole_sum(k,m) = pole_sum(k,m) + pole_tmp(i)
                      enddo
                   endif
                enddo
             enddo
             deallocate( pole_tmp )
             call t_stopf("pdadj_sumfix")
          endif
       endif

       ! save results
       !$omp parallel do private(i,k,m)
       do k = 1, km
          ! normalize first
          do m = 1, ntotq+2
             pole_sum(k,m) = pole_sum(k,m)/im
          enddo
          do i = ifirstxy,ilastxy
             delpxy(i,jm,k) = pole_sum(k,1)
             ptxy(i,jm,k)   = pole_sum(k,2)
          enddo
          do m = 1, ntotq
             do i = ifirstxy,ilastxy
                tracer(i,jm,k,m) = pole_sum(k,m+2)
             enddo
          enddo
       enddo

    endif ! jlastxy==jm

    !
    ! Fix pe,ps if nontrivial z decomposition
    ! Transpose pe - change to better method (16-byte?) later on
    !

    !
    ! Compute pexy
    !
    !$omp parallel do private(i, j)
    do j = jfirstxy,jlastxy
       do i = ifirstxy, ilastxy
          pexy(i,1,j) = ptop
       enddo
    enddo

    !$omp parallel do private(i, j, k)
    do j = jfirstxy,jlastxy
       do k = 1, km
          do i = ifirstxy, ilastxy
             pexy(i,k+1,j) = pexy(i,k,j) + delpxy(i,j,k)
          enddo
       enddo
    enddo

    do j=jfirstxy,jlastxy
       do i=ifirstxy,ilastxy
          psxy(i,j) = pexy(i,km+1,j)
       enddo
    enddo

    !$omp parallel do private(i, j, k)
    do j=jfirstxy,jlastxy
       
       !
       ! Update pelnxy and pkxy
       !
       do k = 1, km+1
          do i = ifirstxy, ilastxy
             pelnxy(i,k,j) = log( pexy(i,k,j) )
             pkxy  (i,j,k) = pexy(i,k,j) ** cappa
          enddo
       enddo

       !
       ! Update pkzxy
       !
       do k = 1,km
          do i = ifirstxy,ilastxy
             pkzxy(i,j,k) = (pkxy(i,j,k+1)-pkxy(i,j,k))/(cappa*(pelnxy(i,k+1,j)-pelnxy(i,k,j)))
          enddo
       enddo
    enddo     ! jfirstxy:jlastxy loop

    !
    ! Calculate virtual potential temperature
    !

    !$omp parallel do private(i, j, k)
    do j=jfirstxy,jlastxy
       do k = 1,km
          do i = ifirstxy,ilastxy
             ptxy(i,j,k) = ptxy(i,j,k)*                  &
                (D1_0+zvir*tracer(i,j,k,1)) &
                /pkzxy(i,j,k)
          enddo
       enddo
    enddo     ! jfirstxy:jlastxy loop

    !EOC
 end subroutine p_d_adjust
!-----------------------------------------------------------------------
