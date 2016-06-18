 module mpas_add_anal_incr

 use mpas_derived_types
 use mpas_pool_routines
 use mpas_kind_types
 use mpas_dmpar
 use mpas_constants   

 !public :: atm_compute_iau_coef, atm_add_tend_anal_incr  
      
 contains

!==================================================================================================
 subroutine atm_compute_iau_coef (iopt, itimestep, dt, time_window_sec, wgt_iau)
!==================================================================================================
! Compute the coefficient (or weight) for the IAU forcing at itimestep.
 
      implicit none

      integer,           intent(in)  :: iopt, itimestep
      real (kind=RKIND), intent(in)  :: dt, time_window_sec
      real (kind=RKIND), intent(out) :: wgt_iau
      
      integer :: nsteps_iau       ! Total number of time steps where the IAU forcing is applied.
      
!      type (mpas_pool_type), intent(in) :: configs
!      type (MPAS_Time_type) :: startTime, stopTime      ! for the entire model integration period
!      type (MPAS_Time_type) :: time_begin, time_end     ! for the IAU window
!      type (MPAS_TimeInterval_type) :: runDuration
!      integer :: local_err       
!      character(len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time
!      real (kind=RKIND),      pointer :: time_window_sec, runtime_window


! Initialization
      wgt_iau = 0.
      

!     For config_IAU_option >= 1, we compute a weighting function here based on the time info in namelist.atmosphere.
!     The default option (config_IAU_option = 1) defines a constant forcing with the same weight
!     (= config_IAU_window_length_s/config_dt + 1) during the IAU time window.      
!     The model is assumed to be further advanced after the forcing (or the filtering) applied (as a free run),
!     we need to fill up the weighting function with zeros for the period from the end of the IAU window
!     all the way to config_stop_time (or for the rest of config_run_duration).      
!      call mpas_pool_get_config(configs, 'config_start_time',   config_start_time)
!      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
!      call mpas_pool_get_config(configs, 'config_stop_time',    config_stop_time)
!      call mpas_pool_get_config(configs, 'config_dt',           config_dt)
!      call mpas_pool_get_config(configs, 'config_IAU_window_length_s', time_window_sec)

      if(iopt == 1) then   ! HA: We should be able to expand this for more options later on.
            
         nsteps_iau = time_window_sec / dt + 1    ! FIXME: Do I need to apply an intrinsic function (INT) to convert real to integer? 
      
         if(itimestep <= nsteps_iau) then
            wgt_iau = 1./nsteps_iau                ! FIXME: Do I need to do FLOAT(nsteps_iau) here?
         end if
      
      end if  


!     FIXME: For now, the IAU forcing is applied from the initial time. If the start of the IAU window might be different,
!     I may need to introduce a new namelist parameter as config_IAU_start_time or something and figure out the time info.
!      
!      call mpas_set_time(startTime,  dateTimeString=config_start_time, ierr=local_err)
!      call mpas_set_time(time_begin, dateTimeString=config_start_time, ierr=local_err)
!      
!      if (trim(config_run_duration) /= "none") then
!         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=local_err)
!
!         if (trim(config_stop_time) /= "none") then
!            call mpas_set_time(stopTime, dateTimeString=config_stop_time, ierr=local_err)
!            if(startTime + runduration /= stopTime) then
!               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
!            end if
!         end if
!      else if (trim(config_stop_time) /= "none") then
!         call mpas_set_time(stopTime, dateTimeString=config_stop_time, ierr=local_err)
!         runDuration = stopTime - startTime
!      else
!          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
!          ierr = 1
!      end if
!       
!      call mpas_get_timeInterval(runDuration, dt=runtime_window)

      end subroutine atm_compute_iau_coef
      
!==================================================================================================
 subroutine atm_add_tend_anal_incr (tend, tend_iau, state, diag, mesh, wgt_iau)
!==================================================================================================
      
      implicit none

      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(inout) :: tend_iau
      type (mpas_pool_type), intent(in) :: state
      type (mpas_pool_type), intent(in) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      real (kind=RKIND),     intent(in) :: wgt_iau
      
      integer :: iEdge, iCell, i, j, k, n
      integer, pointer :: nCells, nEdges, nCellsSolve, nEdgesSolve, nVertLevels
      integer, pointer:: index_qv, index_qc, index_qr, index_qi, index_qs, index_qg
      
      real (kind=RKIND), dimension(:,:), pointer :: rho_edge, rho_zz, theta_m, theta, u, qv, zz, &
                                                    tend_u, tend_rho, tend_theta, tend_th, &
						    tend_qv, tend_qc, tend_qr, tend_qi, tend_qs, tend_qg
      real(kind=RKIND),dimension(:,:,:), pointer :: scalars, tend_scalars									   
      real(kind=RKIND),dimension(:,:), pointer:: u_amb, theta_amb, rho_amb, &
                                                 qv_amb, qc_amb, qr_amb, qi_amb, qs_amb, qg_amb

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      
      call mpas_pool_get_array(mesh, 'zz', zz)
      
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz',   rho_zz, 2)   ! FIXME: what does 2 mean here?
      call mpas_pool_get_array(diag , 'rho_edge', rho_edge)
      
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_dimension(state, 'index_qc', index_qc)
      call mpas_pool_get_dimension(state, 'index_qr', index_qr)
      call mpas_pool_get_dimension(state, 'index_qi', index_qi)
      call mpas_pool_get_dimension(state, 'index_qs', index_qs)
      call mpas_pool_get_dimension(state, 'index_qg', index_qg)
      qv => scalars(index_qv,:,:)

      call mpas_pool_get_array(tend, 'u', tend_u)
      call mpas_pool_get_array(tend, 'rho_zz',  tend_rho)
      call mpas_pool_get_array(tend, 'theta_m', tend_theta)
      call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)
      tend_qv => tend_scalars(index_qv,:,:)
      tend_qc => tend_scalars(index_qc,:,:)
      tend_qr => tend_scalars(index_qr,:,:)
      tend_qi => tend_scalars(index_qi,:,:)
      tend_qs => tend_scalars(index_qs,:,:)
      tend_qg => tend_scalars(index_qg,:,:)

      call mpas_pool_get_array(tend_iau, 'theta', theta_amb)
      call mpas_pool_get_array(tend_iau, 'rho',     rho_amb)
      call mpas_pool_get_array(tend_iau, 'u',         u_amb)
      call mpas_pool_get_array(tend_iau, 'qv',       qv_amb)
      call mpas_pool_get_array(tend_iau, 'qc',       qc_amb)
      call mpas_pool_get_array(tend_iau, 'qr',       qr_amb)
      call mpas_pool_get_array(tend_iau, 'qi',       qi_amb)
      call mpas_pool_get_array(tend_iau, 'qs',       qs_amb)
      call mpas_pool_get_array(tend_iau, 'qg',       qg_amb)
      
      allocate(theta(nVertLevels,nCellsSolve)  )
      allocate(tend_th(nVertLevels,nCellsSolve))

!     initialize the tendency for potential temperature
      tend_th = 0._RKIND
      
!     add coupled tendencies for u on edges
      do i = 1, nEdgesSolve
         do k = 1, nVertLevels
            tend_u(k,i) = tend_u(k,i) + wgt_iau * u_amb(k,i) * rho_edge(k,i)
         enddo
      enddo

!     add tendencies for rho_zz (instead of rho)
      do i = 1, nCellsSolve
         do k = 1, nVertLevels
            tend_rho(k,i) = tend_rho(k,i) + wgt_iau * rho_amb(k,i)/zz(k,i)
         enddo
      enddo
      
!     add coupled tendencies for other state variables on cell centers
      do i = 1, nCellsSolve
         do k = 1, nVertLevels
            tend_th(k,i)  = tend_th(k,i)  + wgt_iau * theta_amb(k,i)*rho_zz(k,i)
            tend_qv(k,i)  = tend_qv(k,i)  + wgt_iau *    qv_amb(k,i)*rho_zz(k,i)
            tend_qc(k,i)  = tend_qc(k,i)  + wgt_iau *    qc_amb(k,i)*rho_zz(k,i)
            tend_qr(k,i)  = tend_qr(k,i)  + wgt_iau *    qr_amb(k,i)*rho_zz(k,i)
            tend_qi(k,i)  = tend_qi(k,i)  + wgt_iau *    qi_amb(k,i)*rho_zz(k,i)
            tend_qs(k,i)  = tend_qs(k,i)  + wgt_iau *    qs_amb(k,i)*rho_zz(k,i)
            tend_qg(k,i)  = tend_qg(k,i)  + wgt_iau *    qg_amb(k,i)*rho_zz(k,i)
         enddo
      enddo

      !if non-hydrostatic core, convert the tendency for the potential temperature to a
      !tendency for the modified potential temperature
      do i = 1, nCellsSolve
         do k = 1, nVertLevels
            theta(k,i) = theta_m(k,i) / (1._RKIND + rvord * qv(k,i))
            tend_th(k,i) = (1. + rvord * qv(k,i))   * tend_th(k,i) &
                               + rvord * theta(k,i) * tend_qv(k,i)
            tend_theta(k,i) = tend_theta(k,i) + tend_th(k,i)
         enddo
      enddo
      
      deallocate(theta)
      deallocate(tend_th)

      
 end subroutine atm_add_tend_anal_incr
 !==================================================================================================

end module mpas_add_anal_incr
