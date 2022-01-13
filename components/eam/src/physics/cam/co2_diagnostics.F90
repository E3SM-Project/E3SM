module co2_diagnostics

!-------------------------------------------------------------------------------
! Purpose:
! Gather CO2 mass, surface fluxes, and aircraft emissions,
! then check global mean carbon mass conservation for different periods
! Conservation check frequency controlled by following namelist flags:
!   co2_print_diags_timestep - time step level checking (default off)
!   co2_print_diags_monthly  - monthly level checking   (default on)
!   co2_print_diags_total    - full run checking        (default on)
!
! Author: Bryce Harrop
!-------------------------------------------------------------------------------

use shr_kind_mod   , only: r8 => shr_kind_r8
use camsrfexch     , only: cam_in_t
use co2_cycle      , only: c_i, co2_transport, co2_print_diags_timestep, &
                           co2_print_diags_monthly, co2_print_diags_total, &
                           co2_conserv_error_tol
use ppgrid         , only: pver, pcols, begchunk, endchunk
use physics_types  , only: physics_state, physics_tend, physics_ptend, &
                           physics_ptend_init
use constituents   , only: pcnst, cnst_name
use cam_logfile    , only: iulog
use spmd_utils     , only: masterproc
use cam_abortutils , only: endrun
use time_manager   , only: is_first_step, is_last_step, get_prev_date, &
                           get_curr_date, is_end_curr_month, & 
                           is_first_restart_step

implicit none
private
save

public co2_diags_init
public get_total_carbon
public get_carbon_sfc_fluxes
public get_carbon_air_fluxes
public print_global_carbon_diags

! Number of CO2 tracers
integer, parameter :: ncnst = 4      ! number of CO2 constituents

character(len=7), dimension(ncnst), parameter :: & ! constituent names
     c_names = (/'CO2_OCN', 'CO2_FFF', 'CO2_LND', 'CO2    '/)

integer :: co2_ocn_glo_ind ! global index of 'CO2_OCN'
integer :: co2_fff_glo_ind ! global index of 'CO2_FFF'
integer :: co2_lnd_glo_ind ! global index of 'CO2_LND'
integer :: co2_glo_ind     ! global index of 'CO2'

!----- formats -----
character(*),parameter :: C_FA0   = "('    ',12x,(42x,a10,2x),' | ',(3x,a10,2x))"
character(*),parameter :: C_FF    = "('    ',a41,e25.17,' | ',e25.17)"
character(*),parameter :: C_FS_2  = "('    ',6x,e25.17,5x,e25.17,5x,' | ',e25.17)"
character(*),parameter :: C_FS2_2 = "('    ',a12,11x,e25.17,18x,' | ',e25.17)"
character(*),parameter :: C_SA0_2 = "('    ',8x,2(11x,a3,15x),' |',(8x,a12,2x))"
character(*),parameter :: C_RER   = "('    ',a15,8x,e25.17,21x)"


!-------------------------------------------------------------------------------
contains


   subroutine co2_diags_init(state)
      !-------------------------------------------------
      ! Purpose: initialize state co2 fields to zero
      !-------------------------------------------------

      !arguments
      type(physics_state), intent(inout) :: state(begchunk:endchunk)

      ! local variables
      integer :: ncol, lchnk, i

      ! zero out all of the state co2 variables
      do lchnk = begchunk, endchunk
         ncol = state(lchnk)%ncol
         do i = 1, ncol
            ! carbon at current, initial, month start, and previous time steps
            state(lchnk)%tc_curr(i) = 0._r8
            state(lchnk)%tc_init(i) = 0._r8
            state(lchnk)%tc_mnst(i) = 0._r8
            state(lchnk)%tc_prev(i) = 0._r8
            ! carbon emissions and fluxes at current timestep
            state(lchnk)%c_flux_sfc(i) = 0._r8
            state(lchnk)%c_flux_air(i) = 0._r8
            ! monthly accumulated carbon emissions and fluxes
            state(lchnk)%c_mflx_sfc(i) = 0._r8
            state(lchnk)%c_mflx_air(i) = 0._r8
            state(lchnk)%c_mflx_sff(i) = 0._r8
            state(lchnk)%c_mflx_lnd(i) = 0._r8
            state(lchnk)%c_mflx_ocn(i) = 0._r8
            ! total time integrated carbon emissions and fluxes
            state(lchnk)%c_iflx_sfc(i) = 0._r8
            state(lchnk)%c_iflx_air(i) = 0._r8
            state(lchnk)%c_iflx_sff(i) = 0._r8
            state(lchnk)%c_iflx_lnd(i) = 0._r8
            state(lchnk)%c_iflx_ocn(i) = 0._r8
         end do
      end do

   end subroutine co2_diags_init

!-------------------------------------------------------------------------------

   subroutine get_total_carbon(state, wet_or_dry)
      !-------------------------------------------------
      ! Purpose: sum column carbon and store in state
      ! Called by: phys_run2
      !-------------------------------------------------
      use physconst,      only: rga

      type(physics_state), intent(inout) :: state
      character(len=3),    intent(in   ) :: wet_or_dry ! is co2 mmr wet or dry

      ! local variables
      real(r8) :: tc(state%ncol)            ! vertical integral of total carbon
      integer ncol                          ! number of atmospheric columns
      integer i, k, m                       ! column, level, constant indices
      !------------------------------------------------------------------------

      if ( .not. co2_transport() ) return

      ! Set CO2 global index
      do m = 1, ncnst
         select case (trim(c_names(m)))
         case ('CO2')
            co2_glo_ind = c_i(m)
         end select
      end do

      ! initialize array
      ncol = state%ncol
      do i = 1, ncol
         tc(i) = 0._r8
      end do
      
      ! sum column co2 mass
      select case (trim(wet_or_dry))
      case ('wet')
         do k = 1, pver
            do i = 1, ncol
               tc(i) = tc(i) + state%q(i,k,co2_glo_ind) * state%pdel(i,k)
            end do
         end do
      case ('dry')
         do k = 1, pver
            do i = 1, ncol
               tc(i) = tc(i) + state%q(i,k,co2_glo_ind) * state%pdeldry(i,k)
            end do
         end do
      end select

      do i = 1, ncol
         tc(i) = tc(i) * rga
      end do

      ! copy new value to state
      do i = 1, ncol
         state%tc_curr(i) = tc(i)
      end do

   end subroutine get_total_carbon

!-------------------------------------------------------------------------------

   subroutine get_carbon_sfc_fluxes(state, cam_in, dtime)
      !-------------------------------------------------
      ! Purpose: store surface carbon exchange in state
      ! Called by: tphysac
      !-------------------------------------------------

      type(physics_state), intent(inout) :: state
      type(cam_in_t),      intent(in   ) :: cam_in
      real(r8), intent(in)               :: dtime        ! physics time step

      ! local variables
      real(r8) :: tc(state%ncol)            ! vertical integral of total carbon
      integer ncol                          ! number of atmospheric columns
      integer i, m                          ! column, constant indices
      real(r8) :: sfc_flux(pcols)           ! surface carbon flux
      real(r8) :: sfc_flux_ocn(pcols)       ! surface ocean carbon flux
      real(r8) :: sfc_flux_fff(pcols)       ! surface fossil fuel carbon flux
      real(r8) :: sfc_flux_lnd(pcols)       ! surface land carbon flux
      !------------------------------------------------------------------------

      if ( .not. co2_transport() ) return

      ! Set CO2 global indices
      do m = 1, ncnst
         select case (trim(c_names(m)))
         case ('CO2_OCN')
            co2_ocn_glo_ind = c_i(m)
         case ('CO2_FFF')
            co2_fff_glo_ind = c_i(m)
         case ('CO2_LND')
            co2_lnd_glo_ind = c_i(m)
         case ('CO2')
            co2_glo_ind     = c_i(m)
         end select
      end do

      ! initialize arrays
      ncol  = state%ncol
      do i = 1, ncol
         sfc_flux(i)     = 0._r8
         sfc_flux_ocn(i) = 0._r8
         sfc_flux_fff(i) = 0._r8
         sfc_flux_lnd(i) = 0._r8
      end do

      ! gather surface fluxes
      do i = 1, ncol
         sfc_flux(i)     = sfc_flux(i)     + cam_in%cflx(i,co2_glo_ind)
         sfc_flux_ocn(i) = sfc_flux_ocn(i) + cam_in%cflx(i,co2_ocn_glo_ind)
         sfc_flux_fff(i) = sfc_flux_fff(i) + cam_in%cflx(i,co2_fff_glo_ind)
         sfc_flux_lnd(i) = sfc_flux_lnd(i) + cam_in%cflx(i,co2_lnd_glo_ind)
      end do

      ! put in state
      do i = 1, ncol
         state%c_flux_sfc(i) = sfc_flux(i)
      end do

      ! zero out monthly fluxes at start of each month
      if ( is_start_curr_month() ) then
         do i = 1, ncol
            state%c_mflx_sfc(i) = 0._r8
            state%c_mflx_ocn(i) = 0._r8
            state%c_mflx_sff(i) = 0._r8
            state%c_mflx_lnd(i) = 0._r8
         end do
      end if

      if ( .not. is_first_step() .and. .not. is_first_restart_step() ) then
         do i = 1, ncol
            state%c_iflx_sfc(i) = state%c_iflx_sfc(i) + (sfc_flux(i)     * dtime)
            state%c_iflx_ocn(i) = state%c_iflx_ocn(i) + (sfc_flux_ocn(i) * dtime)
            state%c_iflx_sff(i) = state%c_iflx_sff(i) + (sfc_flux_fff(i) * dtime)
            state%c_iflx_lnd(i) = state%c_iflx_lnd(i) + (sfc_flux_lnd(i) * dtime)
            state%c_mflx_sfc(i) = state%c_mflx_sfc(i) + (sfc_flux(i)     * dtime)
            state%c_mflx_ocn(i) = state%c_mflx_ocn(i) + (sfc_flux_ocn(i) * dtime)
            state%c_mflx_sff(i) = state%c_mflx_sff(i) + (sfc_flux_fff(i) * dtime)
            state%c_mflx_lnd(i) = state%c_mflx_lnd(i) + (sfc_flux_lnd(i) * dtime)
         end do
      end if

   end subroutine get_carbon_sfc_fluxes

!-------------------------------------------------------------------------------

   subroutine get_carbon_air_fluxes(state, pbuf, dtime)
      !-------------------------------------------------
      ! Purpose: store aircraft CO2 emissions in state
      ! Called by: tphysac
      !-------------------------------------------------
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

      type(physics_state), intent(inout) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      real(r8), intent(in)               :: dtime        ! physics time step

      ! local variables
      real(r8) :: tc(state%ncol)            ! vertical integral of total carbon
      integer ncol                          ! number of atmospheric columns
      integer i, k, m                       ! column, level, constant indices
      integer index_ac_CO2                  ! pbuf index for aircraft emissions
      real(r8), pointer :: ac_CO2(:,:)      ! aircraft emissions in pbuf
      real(r8) :: air_flux(pcols)           ! aircraft carbon flux
      !------------------------------------------------------------------------

      if ( .not. co2_transport() ) return

      ! Set CO2 global index
      do m = 1, ncnst
         select case (trim(c_names(m)))
         case ('CO2')
            co2_glo_ind = c_i(m)
         end select
      end do

      ! acquire aircraft fluxes from physics buffer
      index_ac_CO2 = pbuf_get_index('ac_CO2')   
      call pbuf_get_field(pbuf, index_ac_CO2, ac_CO2)

      ! initialize arrays
      ncol  = state%ncol
      do i = 1, ncol
         air_flux(i) = 0._r8
      end do

      ! gather aircraft fluxes
      do k = 1, pver
         do i = 1, ncol
            air_flux(i) = air_flux(i) + ac_CO2(i,k)
         end do
      end do

      ! put in state
      do i = 1, ncol
         state%c_flux_air(i) = air_flux(i)
      end do

      ! zero out monthly fluxes at start of each month
      if ( is_start_curr_month() ) then
         do i = 1, ncol
            state%c_mflx_air(i) = 0._r8
         end do
      end if

      if ( .not. is_first_step() .and. .not. is_first_restart_step() ) then
         do i = 1, ncol
            state%c_iflx_air(i) = state%c_iflx_air(i) + (air_flux(i) * dtime)
            state%c_mflx_air(i) = state%c_mflx_air(i) + (air_flux(i) * dtime)
         end do
      end if

   end subroutine get_carbon_air_fluxes

!-------------------------------------------------------------------------------

   subroutine print_global_carbon_diags(state, dtime, nstep)
      !-------------------------------------------------
      ! Purpose: Write out conservation checks
      ! Called by: phys_run2
      !-------------------------------------------------
      use phys_gmean,     only: gmean
      use phys_grid,      only: get_ncols_p

      type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
      real(r8), intent(in) :: dtime        ! physics time step
      integer , intent(in) :: nstep        ! current timestep number

      ! local variables
      integer :: ncol, lchnk, i
      integer :: ierr
      integer :: cdate, year, mon, day, sec
      integer, parameter :: c_num_var     = 4
      integer, parameter :: f_ts_num_var  = 2
      integer, parameter :: f_mon_num_var = 5
      integer, parameter :: f_run_num_var = 5
      character(len=*), parameter :: sub_name='print_global_carbon_diags: '
      real(r8) :: time_integrated_flux, state_net_change
      real(r8) :: tc_glob(c_num_var)
      real(r8) :: flux_ts_glob(f_ts_num_var)
      real(r8) :: flux_mon_glob(f_mon_num_var)
      real(r8) :: flux_run_glob(f_run_num_var)
      real(r8) :: gtc_curr, gtc_init, gtc_mnst, gtc_prev, gtc_delta
      real(r8) :: gtc_flux_sfc, gtc_flux_air
      real(r8) :: gtc_mflx_sfc, gtc_mflx_air, gtc_mflx_sff, gtc_mflx_lnd, gtc_mflx_ocn
      real(r8) :: gtc_iflx_sfc, gtc_iflx_air, gtc_iflx_sff, gtc_iflx_lnd, gtc_iflx_ocn
      real(r8) :: gtc_flux_tot, gtc_mflx_tot, gtc_iflx_tot
      real(r8) :: rel_error, expected_tc
      real(r8) :: rel_error_mon, expected_tc_mon
      real(r8) :: rel_error_run, expected_tc_run

      real(r8) :: tc(      pcols,begchunk:endchunk,c_num_var)     ! array for holding carbon variables
      real(r8) :: flux_ts( pcols,begchunk:endchunk,f_ts_num_var)  ! array for holding timestep fluxes
      real(r8) :: flux_mon(pcols,begchunk:endchunk,f_mon_num_var) ! array for holding monthly fluxes
      real(r8) :: flux_run(pcols,begchunk:endchunk,f_run_num_var) ! array for holding full run fluxes
      !------------------------------------------------------------------------

      if ( .not. co2_transport() ) return

      do lchnk = begchunk, endchunk
         ncol = get_ncols_p(lchnk)
         do i = 1, ncol
            ! total carbon mass at different time points
            tc(i,lchnk,1) = state(lchnk)%tc_curr(i)
            tc(i,lchnk,2) = state(lchnk)%tc_init(i)
            tc(i,lchnk,3) = state(lchnk)%tc_mnst(i)
            tc(i,lchnk,4) = state(lchnk)%tc_prev(i)
            ! carbon emissions and fluxes at current time step
            flux_ts(i,lchnk,1) = state(lchnk)%c_flux_sfc(i)
            flux_ts(i,lchnk,2) = state(lchnk)%c_flux_air(i)
            ! monthly accumulated carbon emissions and fluxes
            flux_mon(i,lchnk,1) = state(lchnk)%c_mflx_sfc(i)
            flux_mon(i,lchnk,2) = state(lchnk)%c_mflx_air(i)
            flux_mon(i,lchnk,3) = state(lchnk)%c_mflx_sff(i)
            flux_mon(i,lchnk,4) = state(lchnk)%c_mflx_lnd(i)
            flux_mon(i,lchnk,5) = state(lchnk)%c_mflx_ocn(i)
            ! total time integrated carbon emissions and fluxes
            flux_run(i,lchnk,1) = state(lchnk)%c_iflx_sfc(i)
            flux_run(i,lchnk,2) = state(lchnk)%c_iflx_air(i)
            flux_run(i,lchnk,3) = state(lchnk)%c_iflx_sff(i)
            flux_run(i,lchnk,4) = state(lchnk)%c_iflx_lnd(i)
            flux_run(i,lchnk,5) = state(lchnk)%c_iflx_ocn(i)
         end do
      end do



      ! Compute global means of carbon variables
      if ( ( co2_print_diags_timestep ) .or. &
           ( co2_print_diags_monthly .and. is_end_curr_month() ) .or. & 
           ( co2_print_diags_total .and. is_last_step() ) ) then
         call gmean(tc, tc_glob, c_num_var)
      end if

      if ( co2_print_diags_timestep) then
         call gmean(flux_ts,  flux_ts_glob,  f_ts_num_var)
      end if
      if ( co2_print_diags_monthly .and. is_end_curr_month() ) then
         call gmean(flux_mon, flux_mon_glob, f_mon_num_var)
      end if
      if ( co2_print_diags_total .and. is_last_step() ) then
         call gmean(flux_run, flux_run_glob, f_run_num_var)
      end if

      ! assign global means to readable variables
      gtc_curr     = tc_glob(1)
      gtc_init     = tc_glob(2)
      gtc_mnst     = tc_glob(3)
      gtc_prev     = tc_glob(4)

      gtc_flux_sfc = flux_ts_glob(1)
      gtc_flux_air = flux_ts_glob(2)

      gtc_mflx_sfc = flux_mon_glob(1)
      gtc_mflx_air = flux_mon_glob(2)
      gtc_mflx_sff = flux_mon_glob(3)
      gtc_mflx_lnd = flux_mon_glob(4)
      gtc_mflx_ocn = flux_mon_glob(5)

      gtc_iflx_sfc = flux_run_glob(1)
      gtc_iflx_air = flux_run_glob(2)
      gtc_iflx_sff = flux_run_glob(3)
      gtc_iflx_lnd = flux_run_glob(4)
      gtc_iflx_ocn = flux_run_glob(5)

      ! Compute important terms
      gtc_flux_tot    = gtc_flux_sfc + gtc_flux_air
      gtc_mflx_tot    = gtc_mflx_sfc + gtc_mflx_air
      gtc_iflx_tot    = gtc_iflx_sfc + gtc_iflx_air

      expected_tc     = gtc_prev + (gtc_flux_tot * dtime)
      rel_error       = ( expected_tc - gtc_curr ) / gtc_curr

      expected_tc_mon = gtc_mnst + gtc_mflx_tot ! dtime factor already included
      rel_error_mon   = (expected_tc_mon - gtc_curr) / gtc_curr

      expected_tc_run = gtc_init + gtc_iflx_tot ! dtime factor already included
      rel_error_run   = (expected_tc_run - gtc_curr) / gtc_curr

      gtc_delta    = gtc_curr - gtc_prev

      ! get the date
      call get_curr_date(year, mon, day, sec)
      cdate = year*10000 + mon*100 + day

      ! Time step level write outs----------------------------------------------
      if (masterproc .and. co2_print_diags_timestep .and. (.not. is_first_restart_step())) then
         write(iulog,*   )   ''
         write(iulog,*   )   'NET CO2 FLUXES : period = timestep : date = ',cdate,sec
         write(iulog,C_FA0 ) '  Time  ',   '  Time    '
         write(iulog,C_FA0 ) 'averaged',   'integrated'
         write(iulog,C_FA0 ) 'kgCO2/m2/s', 'kgCO2/m2'

         write(iulog, '(71("-"),"|",20("-"))')

         write(iulog,C_FF) 'Surface  Emissions', gtc_flux_sfc, gtc_flux_sfc * dtime
         write(iulog,C_FF) 'Aircraft Emissions', gtc_flux_air, gtc_flux_air * dtime

         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,C_FF) '   *SUM*', &
              gtc_flux_tot, gtc_flux_tot * dtime

         time_integrated_flux = gtc_flux_tot * dtime

         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,*)''
         write(iulog,*)'CO2 MASS (kgCO2/m2) : period = timestep : date = ',cdate,sec

         write(iulog,*)''
         write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
         write(iulog,C_FS_2) gtc_prev, gtc_curr, gtc_delta


         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,C_FS2_2)'       *SUM*', &
              (gtc_curr - gtc_prev), &
              (gtc_curr - gtc_prev)

         state_net_change = (gtc_curr - gtc_prev)

         write(iulog,C_RER)'Relative Error:', rel_error

         if (nstep > 0) then
            if (abs(rel_error) > co2_conserv_error_tol) then
               write(iulog,*) 'time integrated flux = ', time_integrated_flux
               write(iulog,*) 'net change in state  = ', state_net_change
               write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
               call endrun(trim(sub_name) // 'CO2 conservation failure detected')
            end if
         end if

         write(iulog, '(71("-"),"|",23("-"))')
      end if

      ! Whole run write outs----------------------------------------------------
      if ( is_last_step() .and. co2_print_diags_total ) then
         if (masterproc) then
            write(iulog,*   )   ''
            write(iulog,*   )   'NET CO2 FLUXES : period = full run : date = ',cdate,sec
            write(iulog,C_FA0 ) '  Time  ',   '  Time    '
            write(iulog,C_FA0 ) 'averaged',   'integrated'
            write(iulog,C_FA0 ) 'kgCO2/m2/s', 'kgCO2/m2'

            write(iulog, '(71("-"),"|",20("-"))')

            write(iulog,C_FF) 'Accumulated Surface Flux      ', gtc_iflx_sfc / dtime, gtc_iflx_sfc
            write(iulog,C_FF) 'Accumulated Aircraft Emissions', gtc_iflx_air / dtime, gtc_iflx_air

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FF) '   *SUM*', &
                 gtc_iflx_tot / dtime, gtc_iflx_tot

            time_integrated_flux = gtc_iflx_tot

            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) 'Accumulated Sfc Fssl Fuel Flux', gtc_iflx_sff / dtime, gtc_iflx_sff
            write(iulog,C_FF) 'Accumulated Land  Surface Flux ', gtc_iflx_lnd / dtime, gtc_iflx_lnd
            write(iulog,C_FF) 'Accumulated Ocean Surface Flux', gtc_iflx_ocn / dtime, gtc_iflx_ocn
            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) '   *SUM*', &
                 (gtc_iflx_sff + gtc_iflx_lnd + gtc_iflx_ocn) / dtime, &
                 (gtc_iflx_sff + gtc_iflx_lnd + gtc_iflx_ocn)

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,*)''
            write(iulog,*)'CO2 MASS (kgCO2/m2) : period = full run : date = ',cdate,sec

            write(iulog,*)''
            write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
            write(iulog,C_FS_2) gtc_init, gtc_curr, (gtc_curr - gtc_init)


            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FS2_2)'       *SUM*', &
                 (gtc_curr - gtc_init), &
                 (gtc_curr - gtc_init)

            state_net_change = (gtc_curr - gtc_init)

            write(iulog,C_RER)'Relative Error:', rel_error_run

            if (nstep > 0) then
               if (abs(rel_error_run) > co2_conserv_error_tol) then
                  write(iulog,*) 'time integrated flux = ', time_integrated_flux
                  write(iulog,*) 'net change in state  = ', state_net_change
                  write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
                  write(iulog,*) 'No point in erroring out now, but long-term carbon conservation is bad'
               end if
            end if

            write(iulog, '(71("-"),"|",23("-"))')
         end if ! (masterproc)
      end if ! ( is_last_step() )

      ! Monthly write outs------------------------------------------------------
      if ( is_end_curr_month() .and. co2_print_diags_monthly) then
         if (masterproc) then
            write(iulog,*   )   ''
            write(iulog,*   )   'NET CO2 FLUXES : period = monthly,: date = ',cdate,sec
            write(iulog,C_FA0 ) '  Time  ',   '  Time    '
            write(iulog,C_FA0 ) 'averaged',   'integrated'
            write(iulog,C_FA0 ) 'kgCO2/m2/s', 'kgCO2/m2'

            write(iulog, '(71("-"),"|",20("-"))')

            write(iulog,C_FF) 'Accumulated Surface Flux      ', gtc_mflx_sfc / dtime, gtc_mflx_sfc
            write(iulog,C_FF) 'Accumulated Aircraft Emissions', gtc_mflx_air / dtime, gtc_mflx_air

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FF) '   *SUM*', &
                 gtc_mflx_tot / dtime, gtc_mflx_tot

            time_integrated_flux = gtc_mflx_tot

            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) 'Accumulated Sfc Fssl Fuel Flux', gtc_mflx_sff / dtime, gtc_mflx_sff
            write(iulog,C_FF) 'Accumulated Land  Surface Flux ', gtc_mflx_lnd / dtime, gtc_mflx_lnd
            write(iulog,C_FF) 'Accumulated Ocean Surface Flux', gtc_mflx_ocn / dtime, gtc_mflx_ocn
            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) '   *SUM*', &
                 (gtc_mflx_sff + gtc_mflx_lnd + gtc_mflx_ocn) / dtime, &
                 (gtc_mflx_sff + gtc_mflx_lnd + gtc_mflx_ocn)

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,*)''
            write(iulog,*)'CO2 MASS (kgCO2/m2) : period = monthly,: date = ',cdate,sec

            write(iulog,*)''
            write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
            write(iulog,C_FS_2) gtc_mnst, gtc_curr, (gtc_curr - gtc_mnst)


            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FS2_2)'       *SUM*', &
                 (gtc_curr - gtc_mnst), &
                 (gtc_curr - gtc_mnst)

            state_net_change = (gtc_curr - gtc_mnst)

            write(iulog,C_RER)'Relative Error:', rel_error_mon

            if (nstep > 0) then
               if (abs(rel_error_mon) > co2_conserv_error_tol) then
                  write(iulog,*) 'time integrated flux = ', time_integrated_flux
                  write(iulog,*) 'net change in state  = ', state_net_change
                  write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
                  call endrun(trim(sub_name) // 'Monthly CO2 conservation failure detected')
               end if
            end if

            write(iulog, '(71("-"),"|",23("-"))')
         end if ! (masterproc)
      end if ! ( is_last_step() )

   end subroutine print_global_carbon_diags

!-------------------------------------------------------------------------------

   logical function is_start_curr_month()
   ! Return true if current timestep is first of the current month
   ! Based on is_end_curr_month in time_manager.F90

   ! Local variables
     integer :: &
        yr,   &! year
        mon,  &! month
        day,  &! day of month
        tod    ! time of day (seconds past 00Z)

     call get_prev_date(yr, mon, day, tod)
     is_start_curr_month = (day == 1  .and.  tod == 0)
   end function is_start_curr_month

!-------------------------------------------------------------------------------

end module co2_diagnostics
