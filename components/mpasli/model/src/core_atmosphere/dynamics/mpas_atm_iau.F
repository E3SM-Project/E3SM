! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_atm_iau

 use mpas_derived_types
 use mpas_pool_routines
 use mpas_kind_types
 use mpas_dmpar
 use mpas_constants   
 use mpas_log, only : mpas_log_write

 !public :: atm_compute_iau_coef, atm_add_tend_anal_incr  
      
 contains

!==================================================================================================
 real (kind=RKIND) function atm_iau_coef(configs, itimestep, dt) result(wgt_iau)
!==================================================================================================
! Compute the coefficient (or weight) for the IAU forcing at itimestep.
 
      implicit none

      type (mpas_pool_type), intent(in) :: configs
      integer,           intent(in)  :: itimestep
      real (kind=RKIND), intent(in)  :: dt

      integer :: nsteps_iau       ! Total number of time steps where the IAU forcing is applied.
      logical, parameter :: debug = .false.

      character(len=StrKIND), pointer :: iau_opt
      real (kind=RKIND), pointer :: time_window_sec
      
!      type (mpas_pool_type), intent(in) :: configs
!      type (MPAS_Time_type) :: startTime, stopTime      ! for the entire model integration period
!      type (MPAS_Time_type) :: time_begin, time_end     ! for the IAU window
!      type (MPAS_TimeInterval_type) :: runDuration
!      integer :: local_err       
!      character(len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time
!      real (kind=RKIND),      pointer :: time_window_sec, runtime_window

      call mpas_pool_get_config(configs, 'config_IAU_option', iau_opt)
      call mpas_pool_get_config(configs, 'config_IAU_window_length_s', time_window_sec)

! Initialization
      wgt_iau = 0.
      

!     For config_IAU_option /= 'off', we compute a weighting function here based on the time info in namelist.atmosphere.
!     The default option (config_IAU_option = 'on') defines a constant forcing with the same weight
!     (= config_IAU_window_length_s/config_dt + 1) during the IAU time window.      
!     The model is assumed to be further advanced after the forcing (or the filtering) applied (as a free run),
!     we need to fill up the weighting function with zeros for the period from the end of the IAU window
!     all the way to config_stop_time (or for the rest of config_run_duration).      
!      call mpas_pool_get_config(configs, 'config_start_time',   config_start_time)
!      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
!      call mpas_pool_get_config(configs, 'config_stop_time',    config_stop_time)
!      call mpas_pool_get_config(configs, 'config_dt',           config_dt)
!      call mpas_pool_get_config(configs, 'config_IAU_window_length_s', time_window_sec)

      if(trim(iau_opt) == 'on') then   ! HA: We should be able to expand this for more options later on.
            
         nsteps_iau = nint(time_window_sec / dt)
         !if(debug) call mpas_log_write('atm_compute_iau_coef: nsteps_iau = $i', intArgs=(/nsteps_iau/))
      
         if(itimestep <= nsteps_iau) then
            !wgt_iau = 1./nsteps_iau     
            wgt_iau = 1.0_RKIND / time_window_sec
            if(debug) call mpas_log_write('atm_compute_iau_coef: wgt_iau = $r', realArgs=(/wgt_iau/))
         end if
      
      end if  

 end function atm_iau_coef
      
!==================================================================================================
 subroutine atm_add_tend_anal_incr (configs, structs, itimestep, dt, tend_ru, tend_rtheta, tend_rho)
!==================================================================================================
      
      implicit none

      type (mpas_pool_type), intent(in) :: configs
      type (mpas_pool_type), intent(inout) :: structs
      integer,           intent(in)  :: itimestep
      real (kind=RKIND), intent(in)  :: dt
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_ru
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rtheta
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rho

      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: tend_iau
      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: mesh
      
      integer :: iEdge, iCell, i, j, k, n
      integer, pointer :: nCells, nEdges, nCellsSolve, nEdgesSolve, nVertLevels
      integer, pointer:: index_qv
      integer, pointer:: moist_start, moist_end
      
      real (kind=RKIND), dimension(:,:), pointer :: rho_edge, rho_zz, theta_m, theta, u, zz, &
                                                    tend_th, tend_w
!                                                    tend_u, tend_rho, tend_theta, tend_th, tend_w
      real(kind=RKIND),dimension(:,:,:), pointer :: scalars, tend_scalars, scalars_amb
      real(kind=RKIND),dimension(:,:), pointer:: u_amb, theta_amb, rho_amb

      real (kind=RKIND) :: wgt_iau

      !
      ! Compute weight for IAU forcing in this timestep, and return if weight
      ! is essentially zero.
      !
      wgt_iau = atm_iau_coef(configs, itimestep, dt)
      if (wgt_iau <= 1.0e-12_RKIND) then
         return
      end if

      call mpas_pool_get_subpool(structs, 'tend', tend)
      call mpas_pool_get_subpool(structs, 'tend_iau', tend_iau)
      call mpas_pool_get_subpool(structs, 'state', state)
      call mpas_pool_get_subpool(structs, 'diag', diag)
      call mpas_pool_get_subpool(structs, 'mesh', mesh)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      
      call mpas_pool_get_array(mesh, 'zz', zz)
      
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz',   rho_zz, 2) 
      call mpas_pool_get_array(diag , 'rho_edge', rho_edge)
      
      call mpas_pool_get_dimension(state, 'moist_start', moist_start)
      call mpas_pool_get_dimension(state, 'moist_end', moist_end)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      ! Joe did not recommend to add w tendecies in IAU.
      !call mpas_pool_get_array(tend, 'w', tend_w)       

!      call mpas_pool_get_array(tend, 'u', tend_u)
!      call mpas_pool_get_array(tend, 'rho_zz',  tend_rho)
!      call mpas_pool_get_array(tend, 'theta_m', tend_theta)
      call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)

      call mpas_pool_get_array(tend_iau, 'theta',     theta_amb)
      call mpas_pool_get_array(tend_iau, 'rho',         rho_amb)
      call mpas_pool_get_array(tend_iau, 'u',             u_amb)
      call mpas_pool_get_array(tend_iau, 'scalars', scalars_amb)
      !call mpas_pool_get_array(tend_iau, 'w',             w_amb)
      
      allocate(theta(nVertLevels,nCellsSolve)  )
      allocate(tend_th(nVertLevels,nCellsSolve))

!     initialize the tendency for potential temperature
      tend_th = 0._RKIND
      
!      call mpas_log_write('atm_add_tend_anal_incr: wgt_iau = $r', realArgs=(/wgt_iau/))

!     add coupled tendencies for u on edges
      do i = 1, nEdgesSolve
         do k = 1, nVertLevels
            tend_ru(k,i) = tend_ru(k,i) + wgt_iau * rho_edge(k,i) * u_amb(k,i)
         enddo
      enddo

!     add tendencies for rho_zz (instead of rho)
      do i = 1, nCellsSolve
         do k = 1, nVertLevels
            tend_rho(k,i) = tend_rho(k,i) + wgt_iau * rho_amb(k,i)/zz(k,i)
         enddo
      enddo
      
!     add tendencies for w (tend_w = 0 at k=1 and k=nVertLevelsP1) - Not tested yet
!      do i = 1, nCellsSolve
!         do k = 2, nVertLevels
!            tend_w(k,i) = tend_w(k,i) + wgt_iau * w_amb(k,i)*rho_zz(k,i)
!         enddo
!      enddo
      
      do i = 1, nCellsSolve
         do k = 1, nVertLevels
            theta(k,i) = theta_m(k,i) / (1._RKIND + rvord * scalars(index_qv,k,i))
         enddo
      enddo

!     add coupled tendencies for other state variables on cell centers
      do i = 1, nCellsSolve
         do k = 1, nVertLevels
            tend_th(k,i) = tend_th(k,i) + wgt_iau * (theta_amb(k,i)*rho_zz(k,i) + theta(k,i)*rho_amb(k,i)/zz(k,i))
            tend_scalars(moist_start:moist_end,k,i) = tend_scalars(moist_start:moist_end,k,i) &
                + wgt_iau * (scalars_amb(moist_start:moist_end,k,i)*rho_zz(k,i) + scalars(moist_start:moist_end,k,i)*rho_amb(k,i)/zz(k,i))
         enddo
      enddo

      !if non-hydrostatic core, convert the tendency for the potential temperature to a
      !tendency for the modified potential temperature
      do i = 1, nCellsSolve
         do k = 1, nVertLevels
            tend_th(k,i) = (1. + rvord * scalars(index_qv,k,i))   * tend_th(k,i) &
                               + rvord * theta(k,i) * tend_scalars(index_qv,k,i)
            tend_rtheta(k,i) = tend_rtheta(k,i) + tend_th(k,i)
         enddo
      enddo
      
      deallocate(theta)
      deallocate(tend_th)

      
 end subroutine atm_add_tend_anal_incr
 !==================================================================================================

end module mpas_atm_iau
