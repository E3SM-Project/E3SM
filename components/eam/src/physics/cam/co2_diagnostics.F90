module co2_diagnostics

!------------------------------------------------------------------------------------------------
 
! Purpose:
! Write out mass of CO2 in each tracer (total, fossil fuel, land, and ocean)
! I wanted this to go in co2_cyle.F90, but it created a circular dependency
! with camsrfexch.F90, so this became its own module.
!
! Author: Bryce Harrop
! 
!                                              
!------------------------------------------------------------------------------------------------

use shr_kind_mod   , only: r8 => shr_kind_r8, cxx =>SHR_KIND_CXX, cl =>SHR_KIND_CL
use camsrfexch     , only: cam_in_t
use co2_cycle      , only: c_i, co2_transport
use ppgrid         , only: pver, pcols, begchunk, endchunk
use physics_types  , only: physics_state, physics_tend, physics_ptend, physics_ptend_init
use constituents   , only: pcnst, cnst_name
use cam_logfile    , only: iulog
use spmd_utils     , only: masterproc
use cam_abortutils , only: endrun
use time_manager   , only: is_first_step, is_last_step

implicit none
private
save

public co2_diags_register            ! setup co2 pbuf fields
public get_total_carbon
public get_carbon_emissions
public print_global_carbon_diags
public print_global_carbon_diags_scl

public check_co2_change_pr2
public co2_gmean_check_wflux         ! printout co2 global means
public co2_gmean_check2_wflux        ! higher level co2 checker

! Number of CO2 tracers
integer, parameter :: ncnst = 4                      ! number of constituents implemented

character(len=7), dimension(ncnst), parameter :: & ! constituent names
     c_names = (/'CO2_OCN', 'CO2_FFF', 'CO2_LND', 'CO2    '/)

integer :: co2_ocn_glo_ind ! global index of 'CO2_OCN'
integer :: co2_fff_glo_ind ! global index of 'CO2_FFF'
integer :: co2_lnd_glo_ind ! global index of 'CO2_LND'
integer :: co2_glo_ind     ! global index of 'CO2'

!----- formats -----
!character(*),parameter :: C_STRD  = "(1x,a30,1x,i8,3(1x,e25.17))"

character(*),parameter :: C_FA0   = "('    ',12x,(42x,a10,2x),' | ',(3x,a10,2x))" ! 71|
!character(*),parameter :: C_FF    = "('    ',a51,f15.8,' | ',f18.2)"
character(*),parameter :: C_FF    = "('    ',a41,e25.17,' | ',e25.17)" ! 71|

character(*),parameter :: FF2     = "('    ',a12,a15,' | ',f18.2)" ! 32|
character(*),parameter :: C_FS    = "('    ',a12,7(f18.2),26x,' | ',(f18.2))" ! 168|
character(*),parameter :: C_FS0   = "('    ',12x,8(a19),' | ',(a19))" ! 169|
character(*),parameter :: C_FS2   = "('    ',a12,67x,f18.2,67x,' | ',f18.2)" ! 169|
character(*),parameter :: C_FS3   = "('    ',a12,8(f18.2),8x,' | ',(f18.2))" ! 161|
!character(*),parameter :: C_FS_2  = "('    ',a25,f15.2,5x,f15.2,5x,' | ',f18.2)"
character(*),parameter :: C_FS_2  = "('    ',6x,e25.17,5x,e25.17,5x,' | ',e25.17)" ! 71|

character(*),parameter :: C_SA0   = "('    ',34x,2(5x,a3,8x),' | ',(8x,a12,2x))" ! 71|
!character(*),parameter :: C_FS2_2 = "('    ',a12,17x,f18.2,18x,' | ',f18.2)"
character(*),parameter :: C_FS2_2 = "('    ',a12,11x,e25.17,18x,' | ',e25.17)" ! 71|
!character(*),parameter :: C_SA0_2 = "('    ',31x,2(5x,a3,9x),' |',(8x,a12,2x))"
character(*),parameter :: C_SA0_2 = "('    ',8x,2(11x,a3,15x),' |',(8x,a12,2x))" ! 71|

!character(*),parameter :: C_FS3_3 = "('    ',a12,53x,' | ',(f18.2))"
character(*),parameter :: C_FS3_3 = "('    ',a12,54x,' | ',(f18.2))" ! 71|

contains

   subroutine co2_diags_register()
      !-------------------------------------------------
      ! Purpose: register co2 fields into pbuf
      !-------------------------------------------------
      use physics_buffer, only: pbuf_add_field, dtype_r8

      integer :: m,idx

      if (co2_transport()) then
         do m = 1,4
!            call outfld(trim(cnst_name(c_i(m)))//'_BOT', state%q(1,pver,c_i(m)), pcols, lchnk)
!            call outfld(sflxnam(c_i(m)), cam_in%cflx(:,c_i(m)), pcols, lchnk)
!            call pbuf_add_field(rflx_pname(i), 'physpkg', dtype_r8, (/pcols/), idx)
!
            call pbuf_add_field(trim(cnst_name(c_i(m)))//'_cur', 'physpkg', dtype_r8, (/pcols/), idx)
         enddo
      endif

   end subroutine co2_diags_register

   subroutine get_total_carbon(state, wet_or_dry)
      use physconst,      only: rga

      type(physics_state), intent(inout) :: state
      character(len=3),    intent(in   ) :: wet_or_dry ! is co2 mmr wet or dry

      ! local variables
      real(r8) :: tc(state%ncol)            ! vertical integral of total carbon
      integer ncol                          ! number of atmospheric columns
      integer i, k, m                       ! column, level, constant indices
      !------------------------------------------------------------------------

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

   subroutine get_carbon_emissions(state, cam_in, pbuf, dtime)
      use physconst,      only: rga
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

      type(physics_state), intent(inout) :: state
      type(cam_in_t),      intent(in   ) :: cam_in
      type(physics_buffer_desc), pointer :: pbuf(:)
      real(r8), intent(in)               :: dtime        ! physics time step

      ! local variables
      real(r8) :: tc(state%ncol)            ! vertical integral of total carbon
      integer ncol                          ! number of atmospheric columns
      integer i, k, m                       ! column, level, constant indices
      integer index_ac_CO2                  ! pbuf index for aircraft emissions
      real(r8), pointer :: ac_CO2(:,:)      ! aircraft emissions in pbuf
      real(r8) :: sfc_flux(pcols)           ! constituent surface flux
      real(r8) :: air_flux(pcols)           ! aircraft emission flux
      !------------------------------------------------------------------------

      ! Set CO2 global index
      do m = 1, ncnst
         select case (trim(c_names(m)))
         case ('CO2')
            co2_glo_ind = c_i(m)
         end select
      end do

      index_ac_CO2 = pbuf_get_index('ac_CO2')   
      call pbuf_get_field(pbuf, index_ac_CO2, ac_CO2)

      ! initialize arrays
      ncol  = state%ncol
      do i = 1, ncol
         sfc_flux(i) = 0._r8
         air_flux(i) = 0._r8
      end do

      ! gather surface and aircraft fluxes
      do i = 1, ncol
         sfc_flux(i) = sfc_flux(i) + cam_in%cflx(i,co2_glo_ind)
      end do

      do k = 1, pver
         do i = 1, ncol
            air_flux(i) = air_flux(i) + ac_CO2(i,k)
         end do
      end do

      ! Note: once future versions include chemical sources and sinks
      !       those will need to be included here too

      ! sum value and put in state
      do i = 1, ncol
         state%c_emis_sfc(i) = sfc_flux(i)
         state%c_emis_air(i) = air_flux(i)
      end do

      if ( .not. is_first_step() ) then
         do i = 1, ncol
            state%c_iflx_sfc(i) = state%c_iflx_sfc(i) + (sfc_flux(i) * dtime)
            state%c_iflx_air(i) = state%c_iflx_air(i) + (air_flux(i) * dtime)
         end do
      end if

   end subroutine get_carbon_emissions

   subroutine print_global_carbon_diags(state, dtime, nstep)
      use phys_gmean,     only: gmean

      type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
      real(r8), intent(in) :: dtime        ! physics time step
      integer , intent(in) :: nstep        ! current timestep number

      ! local variables
      integer ncol, lchnk
      real(r8) :: tc(pcols,begchunk:endchunk,7)
      real(r8) :: tc_glob(7)
!      real(r8) :: tc_start_glob, tc_end_glob, delta_tc_glob, c_emis_glob, c_air_glob
      real(r8) :: gtc_curr, gtc_init, gtc_prev, gtc_delta
      real(r8) :: gtc_emis_sfc, gtc_emis_air, gtc_iflx_sfc, gtc_iflx_air
      real(r8) :: rel_error, expected_tc
      !------------------------------------------------------------------------

      do lchnk = begchunk, endchunk
         ncol = state(lchnk)%ncol
         ! carbon at current time step
         tc(:ncol,lchnk,1) = state(lchnk)%tc_curr(:ncol)
         ! carbon at first time step
         tc(:ncol,lchnk,2) = state(lchnk)%tc_init(:ncol)
         ! carbon for previous time step
         tc(:ncol,lchnk,3) = state(lchnk)%tc_prev(:ncol)
         ! carbon emissions
         tc(:ncol,lchnk,4) = state(lchnk)%c_emis_sfc(:ncol)
         tc(:ncol,lchnk,5) = state(lchnk)%c_emis_air(:ncol)
         ! time integrated carbon fluxes
         tc(:ncol,lchnk,6) = state(lchnk)%c_iflx_sfc(:ncol)
         tc(:ncol,lchnk,7) = state(lchnk)%c_iflx_air(:ncol)
      end do

      ! Compute global means of input and output energies and of
      ! surface pressure for heating rate (assume uniform ptop)
!      if (nstep > 0) then
      call gmean(tc, tc_glob, 7)
!      end if

      gtc_curr     = tc_glob(1)
      gtc_init     = tc_glob(2)
      gtc_prev     = tc_glob(3)
      gtc_emis_sfc = tc_glob(4)
      gtc_emis_air = tc_glob(5)
      gtc_iflx_sfc = tc_glob(6)
      gtc_iflx_air = tc_glob(7)

!!      rel_error     = ( (delta_tc_glob / dtime) - c_emis_glob ) / c_emis_glob
!      expected_tc   = tc_start_glob + (c_emis_glob * dtime) + (c_air_glob * dtime)
!      rel_error     = ( expected_tc - tc_end_glob  ) / tc_end_glob
      expected_tc  = gtc_prev + ( (gtc_emis_sfc + gtc_emis_air) * dtime )
      rel_error    = ( expected_tc - gtc_curr ) / gtc_curr
      gtc_delta    = gtc_curr - gtc_prev

      if (masterproc) then
         write(iulog,'(1x,a30,1x,i8,3(1x,e25.17))') "nstep, tc start, tc end, diff ", nstep, gtc_prev, gtc_curr, gtc_delta
         write(iulog,'(1x,a30,1x,i8,3(1x,e25.17))') "nstep, d(tc)/dt, emis, rel err", nstep, gtc_delta/dtime, gtc_emis_sfc + gtc_emis_air, rel_error
         write(iulog,'(1x,a30,1x,i8,2(1x,e25.17))') "nstep, sfc_emis, aircraft_emis", nstep, gtc_emis_sfc, gtc_emis_air
      end if
      
   end subroutine print_global_carbon_diags

   subroutine print_global_carbon_diags_scl(state, dtime, nstep)
      use phys_gmean,     only: gmean

      type(physics_state), intent(in   ) :: state
      real(r8), intent(in) :: dtime        ! physics time step
      integer , intent(in) :: nstep        ! current timestep number

      ! local variables
      integer ncol, lchnk
      integer ierr
      character(len=*), parameter :: sub_name='print_global_carbon_diags_scl: '
      real(r8), parameter :: error_tol = 0.01_r8
      real(r8) :: time_integrated_flux, state_net_change
      real(r8) :: tc_glob(7)
!      real(r8) :: tc_start_glob, tc_end_glob, delta_tc_glob, c_emis_glob, c_air_glob
      real(r8) :: gtc_curr, gtc_init, gtc_prev, gtc_delta
      real(r8) :: gtc_emis_sfc, gtc_emis_air, gtc_iflx_sfc, gtc_iflx_air
      real(r8) :: rel_error, expected_tc

      real(r8), pointer :: tc(:,:,:) ! array for holding carbon variables
      !------------------------------------------------------------------------

      if ( .not. co2_transport() ) return

      allocate(tc(pcols,begchunk:endchunk,7), stat=ierr)
      if (ierr /= 0) write(iulog,*) trim(sub_name) // 'FAIL to allocate tc'

      ncol  = state%ncol
      lchnk = state%lchnk
      ! carbon at current time step
      tc(:ncol,lchnk,1) = state%tc_curr(:ncol)
      ! carbon at first time step
      tc(:ncol,lchnk,2) = state%tc_init(:ncol)
      ! carbon for previous time step
      tc(:ncol,lchnk,3) = state%tc_prev(:ncol)
      ! carbon emissions
      tc(:ncol,lchnk,4) = state%c_emis_sfc(:ncol)
      tc(:ncol,lchnk,5) = state%c_emis_air(:ncol)
      ! time integrated carbon fluxes
      tc(:ncol,lchnk,6) = state%c_iflx_sfc(:ncol)
      tc(:ncol,lchnk,7) = state%c_iflx_air(:ncol)



      ! Compute global means of input and output energies and of
      ! surface pressure for heating rate (assume uniform ptop)
!      if (nstep > 0) then
      call gmean(tc, tc_glob, 7)
!      end if

      gtc_curr     = tc_glob(1)
      gtc_init     = tc_glob(2)
      gtc_prev     = tc_glob(3)
      gtc_emis_sfc = tc_glob(4)
      gtc_emis_air = tc_glob(5)
      gtc_iflx_sfc = tc_glob(6)
      gtc_iflx_air = tc_glob(7)

!!      rel_error     = ( (delta_tc_glob / dtime) - c_emis_glob ) / c_emis_glob
!      expected_tc   = tc_start_glob + (c_emis_glob * dtime) + (c_air_glob * dtime)
!      rel_error     = ( expected_tc - tc_end_glob  ) / tc_end_glob
      expected_tc  = gtc_prev + ( (gtc_emis_sfc + gtc_emis_air) * dtime )
      rel_error    = ( expected_tc - gtc_curr ) / gtc_curr
      gtc_delta    = gtc_curr - gtc_prev

      if (masterproc) then
         write(iulog,'(1x,a30,1x,i8,3(1x,e25.17))') "nstep, tc start, tc end, diff ", nstep, gtc_prev, gtc_curr, gtc_delta
         write(iulog,'(1x,a30,1x,i8,3(1x,e25.17))') "nstep, d(tc)/dt, emis, rel err", nstep, gtc_delta/dtime, gtc_emis_sfc + gtc_emis_air, rel_error
         write(iulog,'(1x,a30,1x,i8,2(1x,e25.17))') "nstep, sfc_emis, aircraft_emis", nstep, gtc_emis_sfc, gtc_emis_air

         write(iulog,'(1x,a37,1x,i8,4(1x,e25.17))') "nstep, tc beg, tc end, d(tc), flux*dt", nstep, gtc_prev, gtc_curr, gtc_delta, dtime*(gtc_emis_sfc+gtc_emis_air)
      end if

      if (masterproc) then
         write(iulog,*   )   ''
         write(iulog,*   )   'NET CO2 FLUXES'! : period ',trim(pname(ip)),': date = ',cdate,sec
         write(iulog,C_FA0 ) '  Time  ',   '  Time    '
         write(iulog,C_FA0 ) 'averaged',   'integrated'
         write(iulog,C_FA0 ) 'kgCO2/m2/s', 'kgCO2/m2'

         write(iulog, '(71("-"),"|",20("-"))')

         write(iulog,C_FF) 'Surface  Emissions', gtc_emis_sfc, gtc_emis_sfc * dtime
         write(iulog,C_FF) 'Aircraft Emissions', gtc_emis_air, gtc_emis_air * dtime
!         do f = 1, f_size
!            write(iulog,C_FF)c_f_name(f),budg_fluxGpr(f,ip),budg_fluxG(f,ip)*unit_conversion*get_step_size()
!         end do

         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,C_FF) '   *SUM*', &
              (gtc_emis_sfc + gtc_emis_air), (gtc_emis_sfc + gtc_emis_air)*dtime

         time_integrated_flux = (gtc_emis_sfc + gtc_emis_air) * dtime

         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,*)''
         write(iulog,*)'CO2 MASS (kgCO2/m2)'!: period ',trim(pname(ip)),': date = ',cdate,sec

         write(iulog,*)''
         write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
         write(iulog,C_FS_2) gtc_prev, gtc_curr, gtc_delta
!         do s = 1,c_s_name_size-1
!            s_beg = s_totc_beg + s
!            s_end = s_totc_end + s
!            write(iulog,C_FS_2)c_s_name(s),&
!                 budg_stateG(s_beg,ip)*unit_conversion, &
!                 budg_stateG(s_end,ip)*unit_conversion, &
!                 (budg_stateG(s_end,ip) - budg_stateG(s_beg,ip))*unit_conversion
!         end do
!         write(iulog,C_FS_2)'Grid-level Err', 0._r8, rel_error , 0._r8


         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,C_FS2_2)'       *SUM*', &
              (gtc_curr - gtc_prev), &
              (gtc_curr - gtc_prev)

         state_net_change = (gtc_curr - gtc_prev)

         if (nstep > 0) then
            if (abs(time_integrated_flux - state_net_change) > error_tol) then
               write(iulog,*) 'time integrated flux = ', time_integrated_flux
               write(iulog,*) 'net change in state  = ', state_net_change
               write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
               call endrun(trim(sub_name) // 'CO2 conservation failure detected')
            end if
         end if

         write(iulog, '(71("-"),"|",23("-"))')
      end if

      if ( is_last_step() ) then
         if (masterproc) then
            write(iulog,*   )   ''
            write(iulog,*   )   'NET CO2 FLUXES'! : period ',trim(pname(ip)),': date = ',cdate,sec
            write(iulog,C_FA0 ) '  Time  ',   '  Time    '
            write(iulog,C_FA0 ) 'averaged',   'integrated'
            write(iulog,C_FA0 ) 'kgCO2/m2/s', 'kgCO2/m2'

            write(iulog, '(71("-"),"|",20("-"))')

            write(iulog,C_FF) 'Accumulated Surface  Emissions', gtc_iflx_sfc / dtime, gtc_iflx_sfc
            write(iulog,C_FF) 'Accumulated Aircraft Emissions', gtc_iflx_air / dtime, gtc_iflx_air

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FF) '   *SUM*', &
                 (gtc_iflx_sfc + gtc_iflx_air) / dtime, (gtc_iflx_sfc + gtc_iflx_air)

            time_integrated_flux = (gtc_iflx_sfc + gtc_iflx_air)

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,*)''
            write(iulog,*)'CO2 MASS (kgCO2/m2)'!: period ',trim(pname(ip)),': date = ',cdate,sec

            write(iulog,*)''
            write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
            write(iulog,C_FS_2) gtc_init, gtc_curr, (gtc_curr - gtc_init)


            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FS2_2)'       *SUM*', &
                 (gtc_curr - gtc_init), &
                 (gtc_curr - gtc_init)

            state_net_change = (gtc_curr - gtc_init)

            if (nstep > 0) then
               if (abs(time_integrated_flux - state_net_change) > error_tol) then
                  write(iulog,*) 'time integrated flux = ', time_integrated_flux
                  write(iulog,*) 'net change in state  = ', state_net_change
                  write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
                  call endrun(trim(sub_name) // 'Long CO2 conservation failure detected')
               end if
            end if

            write(iulog, '(71("-"),"|",23("-"))')
         end if ! (masterproc)
      end if ! ( is_last_step() )

      !if rel_error > rel_error_threshold then
         !error out with message about conservation
      !end if
      
      deallocate(tc)

   end subroutine print_global_carbon_diags_scl



   subroutine check_co2_change_pr2(state, tend, pbuf2d, cam_in, wet_or_dry)
      use physconst,      only: gravit
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, pbuf_set_field 
      use phys_control,   only: ieflx_opt
      use phys_gmean,     only: gmean

!      integer , intent(in) :: nstep        ! current timestep number
      type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
      type(physics_tend ), intent(in   ), dimension(begchunk:endchunk) :: tend
      type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
      type(physics_buffer_desc),          pointer :: pbuf2d(:,:)
      character(len=3),    intent(in   ) :: wet_or_dry    ! is co2 wet or dry at this point
 
      integer :: ncol                      ! number of active columns
      integer :: lchnk                     ! chunk index

      integer :: i, k, m

      character(len=*), parameter :: title='CO2 end of phys run2:)'
      character(len=*), parameter :: sub_name='co2_gmean_check: '


      real(r8) :: mass_wet_mean(ncnst)
      real(r8) :: mass_dry_mean(ncnst)
      real(r8) :: sfc_flux_mean(ncnst)
      real(r8) :: air_flux_mean(ncnst)

      real(r8) :: mass_wet(pcols, begchunk:endchunk, ncnst)
      real(r8) :: mass_dry(pcols, begchunk:endchunk, ncnst)
      real(r8) :: sfc_flux(pcols, begchunk:endchunk, ncnst)
      real(r8) :: air_flux(pcols, begchunk:endchunk, ncnst)
      

      if ( .not. co2_transport() ) return

      do lchnk = begchunk, endchunk
!         qflx(:ncol,lchnk) = cam_in(lchnk)%cflx(:ncol,1)
         mass_wet(:ncol, lchnk, :) = 0._r8
         mass_dry(:ncol, lchnk, :) = 0._r8
         sfc_flux(:ncol, lchnk, :) = 0._r8
         air_flux(:ncol, lchnk, :) = 0._r8
         do i = 1, ncol
            do k = 1, pver
               do m = 1, ncnst
                  mass_wet(i, lchnk, m) = mass_wet(i, lchnk, m) + &
                       state(lchnk)%pdel(i, k)*state(lchnk)%q(i, k, c_i(m))
                  mass_dry(i, lchnk, m) = mass_dry(i, lchnk, m) + &
                       state(lchnk)%pdel(i, k)*state(lchnk)%q(i, k, c_i(m))
               end do ! m = 1, ncnst
!               air_flux(i, lchnk) = air_flux(i, lchnk) + ac_CO2(i, k)
            end do ! k = 1, pver
            do m = 1, ncnst
               mass_wet(i, lchnk, m) = mass_wet(i, lchnk, m) / gravit
               mass_dry(i, lchnk, m) = mass_dry(i, lchnk, m) / gravit
               sfc_flux(i, lchnk, m) = sfc_flux(i, lchnk, m) + &
                    cam_in(lchnk)%cflx(i, c_i(m))
               air_flux_mean(m) = 0._r8  ! DONT KEEP THIS!!!
            end do ! m = 1, ncnst
         end do ! i = 1, ncol
      end do ! lchnk = begchunk, endchunk

      call gmean(mass_wet, mass_wet_mean, ncnst)
      call gmean(mass_dry, mass_dry_mean, ncnst)
      call gmean(sfc_flux, sfc_flux_mean, ncnst)
!      call gmean(air_flux, air_flux_mean, ncnst)

      if (begchunk .le. endchunk) then
         if (masterproc) then
            do m = 1, ncnst
               write (6,'(a32,i2,a36,1p,4e25.17)') trim(title)//' m=',c_i(m), &
                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc, air ', &
                  mass_dry_mean(m), mass_wet_mean(m), &
                  sfc_flux_mean(m), air_flux_mean(m)
            end do
         end if
      end if
      
   end subroutine check_co2_change_pr2

!   subroutine check_co2_change(wet_or_dry, state, ptend, cam_in, pbuf)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the 
! physics decomposition. Checks for conservation from previous timestep
!
! Author: Bryce Harrop
!
!-----------------------------------------------------------------------
!      use physconst,      only: gravit
!      use phys_grid,      only: get_ncols_p
!      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
!      use phys_gmean,     only: gmean
!
! Arguments
!
!      character(len=3),    intent(in) :: wet_or_dry    ! is co2 wet or dry at this point
!      type(physics_state), intent(in) :: state
!      type(physics_buffer_desc), pointer :: pbuf(:)
!      type(physics_ptend), intent(in) :: ptend
!      type(cam_in_t),      intent(in) :: cam_in
      !integer, dimension(:), intent(in) :: cnst_ind_arr
!
! Local workspace
!
!      character(len=*), parameter :: sub_name='check_co2_change: '

!      integer :: c, i, k, m
!      integer :: ierr, ncol, idx
!      integer :: index_co2_ocn, index_co2_fff, index_co2_lnd, index_co2
!      real(r8)          :: ac_CO2_tot(pcols)
!      real(r8), pointer :: ac_CO2(:,:)

!      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
!      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
!      real(r8) :: mass_wet_mean(ncnst)     ! global mean constituent masses assuming moist mmr
!      real(r8) :: mass_dry_mean(ncnst)     ! global mean constituent masses assuming dry mmr
!      real(r8), pointer :: sfc_flux(:,:,:) ! constituent surface flux
!      real(r8) :: sfc_flux_mean(ncnst)     ! global mean constituent surface flux
!      real(r8), pointer :: air_flux(:,:)   ! aircraft emission flux
!      real(r8) :: air_flux_mean            ! global mean aircraft emission flux

!      real(r8), pointer :: mass_cur(:,:,:) ! previous timestep mass
!      real(r8) :: mass_mean_cur(ncnst)

!      real(r8), pointer :: co2_cur(:,:,:)  ! previous co2 mixing ratios
!
!-----------------------------------------------------------------------
!
!      if ( .not. co2_transport() ) return

!      allocate(mass_wet(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

!      allocate(mass_dry(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

!      allocate(sfc_flux(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate sfc_flux'

!      allocate(air_flux(pcols,begchunk:endchunk), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate air_flux'

!      allocate(mass_cur(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate air_flux'

!      index_co2_ocn = pbuf_get_index(trim(cnst_name(c_i(1)))//'_cur')
!      index_co2_fff = pbuf_get_index(trim(cnst_name(c_i(2)))//'_cur')
!      index_co2_lnd = pbuf_get_index(trim(cnst_name(c_i(3)))//'_cur')
!      index_co2     = pbuf_get_index(trim(cnst_name(c_i(4)))//'_cur')
!      index_acco2   = pbuf_get_index('ac_CO2')  

!      call pbuf_get_field(pbuf, index_co2_ocn, co2_ocn_cur)
!      call pbuf_get_field(pbuf, index_co2_fff, co2_fff_cur)
!      call pbuf_get_field(pbuf, index_co2_lnd, co2_lnd_cur)
!      call pbuf_get_field(pbuf, index_co2    , co2_cur)
!      call pbuf_get_field(pbuf, index_acco2  , ac_CO2)

      ! Set CO2 global indices
!      do idx = 1, ncnst
!         select case (trim(c_names(idx)))
!         case ('CO2_OCN')
!            co2_ocn_glo_ind = c_i(idx)
!         case ('CO2_FFF')
!            co2_fff_glo_ind = c_i(idx)
!         case ('CO2_LND')
!            co2_lnd_glo_ind = c_i(idx)
!         case ('CO2')
!            co2_glo_ind     = c_i(idx)
!         end select
!      end do

!      ncol = state%ncol
!      c    = state%lchnk
!      !ncnst = size(cnst_ind_arr)
!      do m = 1, ncnst
!!         do c = begchunk, endchunk
!!            ncols = get_ncols_p(c)
!            do i = 1, ncol

!               ! Compute column masses assuming both dry and wet mixing ratios

!               mass_wet(i,c,m) = 0.0_r8
!               do k = 1, pver
!                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
!                                    state%pdel(i,k)*state%q(i,k,c_i(m))
!               end do
!               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

!               mass_dry(i,c,m) = 0.0_r8
!               do k = 1, pver
!                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
!                                    state%pdeldry(i,k)*state%q(i,k,c_i(m))
!               end do
!               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

!               sfc_flux(i,c,m) = 0.0_r8
!               sfc_flux(i,c,m) = sfc_flux(i,c,m) + &
!                                 cam_in%cflx(i,c_i(m))
               
!!               ac_CO2_tot(i) = 0.0_r8
!!               do k = 1, pver
!!                  ac_CO2_tot(i) = ac_CO2_tot(i) + ac_CO2(i,k)
!!               end do

!            end do
!!         end do
!!            if (m .eq. co2_glo_ind .or. m .eq. co2_fff_glo_ind) then
!!               air_flux(:ncol,c,m) = ac_CO2_tot(:ncol)
!!            else
!!               air_flux(:ncol,c,m) = 0.0_r8
!!            end if
!      end do

!      do i = 1, ncol
!         ac_CO2_tot(i) = 0.0_r8
!         do k = 1, pver
!            ac_CO2_tot(i) = ac_CO2_tot(i) + ac_CO2(i,k)
!         end do
!         air_flux(:ncol,c) = ac_CO2_tot(:ncol)
!      end do

!      do i = 1, ncol
!         mass_cur(i,c,:) = 0._r8
!         if (wet_or_dry == 'wet') then
!            mass_cur(i,c,1) = mass_wet(i,c,1) + &
!                                 state%pdel(i,k)*co2_ocn_cur(i,k)
!            mass_cur(i,c,2) = mass_wet(i,c,2) + &
!                                 state%pdel(i,k)*co2_fff_cur(i,k)
!            mass_cur(i,c,3) = mass_wet(i,c,3) + &
!                                 state%pdel(i,k)*co2_lnd_cur(i,k)
!            mass_cur(i,c,4) = mass_wet(i,c,4) + &
!                                 state%pdel(i,k)*co2_cur(i,k)
!         end if
!      end do

!      ! compute global mean mass
!      call gmean(mass_wet, mass_wet_mean, ncnst)
!      call gmean(mass_dry, mass_dry_mean, ncnst)
!      call gmean(sfc_flux, sfc_flux_mean, ncnst)
!      call gmean(air_flux, air_flux_mean)

!      call gmean(mass_cur, mass_mean_cur, ncnst)

!      ! report to log file
!      if (masterproc) then

!         do m = 1, ncnst
!               write (6,66) trim(title)//' m=',c_i(m), &
!                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc, air ', &
!                  mass_dry_mean(c_i(m)), mass_wet_mean(c_i(m)), &
!                  sfc_flux_mean(c_i(m)), air_flux_mean(c_i(m))
!66             format (a32,i2,a36,1p,4e25.15)
!!old version66             format (a32,i2,a36,1p,4e25.13)
!         end do

!      endif

!      deallocate(air_flux)
!      deallocate(sfc_flux)
!      deallocate(mass_dry)
!      deallocate(mass_wet)

!   end subroutine check_co2_change

!   subroutine check_co2_change_nope(state, pbuf2d)
      !-----------------------------------------------
      ! Purpose: check timestep level co2 conservation
      !
      ! Author: Bryce Harrop
      !-----------------------------------------------
!      use physics_types,  only: physics_state
!      use ppgrid,         only: begchunk, endchunk, pver, pcols
!      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
!                                pbuf_get_chunk, pbuf_get_index
!      use physconst,      only: rga
!      use co2_cycle,      only: c_i, co2_transport

!      implicit none
      
!      type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
!      type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)

!      !local vars
!      type(physics_buffer_desc), pointer     :: pbuf_chnk(:)

!      real(r8), pointer :: ac_CO2(:,:)

!      real(r8) :: ftem      (pcols,pver) ! tmp space
!      real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
!      real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor

!      real(r8) :: ftem(pcols,pver)         ! tmp space
!      real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
!      real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor
!      real(r8) :: CIDiff(pcols)            ! Difference in vertically integrated static energy

!      integer  :: index_co2_ocn, index_co2_fff, index_co2_lnd, index_co2
!      integer  :: ichnk, ncol, i, k, m, lchnk

!      if ( .not. co2_transport() ) return

!      do ichnk = begchunk, endchunk
!         ncol       =  state(ichnk)%ncol
!         pbuf_chunk => pbuf_get_chunk(pbuf2d, ichnk)

!         index_co2_ocn = pbuf_get_index(trim(cnst_name(c_i(1)))//'_cur')
!         index_co2_fff = pbuf_get_index(trim(cnst_name(c_i(2)))//'_cur')
!         index_co2_lnd = pbuf_get_index(trim(cnst_name(c_i(3)))//'_cur')
!         index_co2     = pbuf_get_index(trim(cnst_name(c_i(4)))//'_cur')
!         index_acco2   = pbuf_get_index('ac_CO2')   

!         call pbuf_get_field(pbuf_chnk, index_co2_ocn, co2_ocn_cur)
!         call pbuf_get_field(pbuf_chnk, index_co2_fff, co2_fff_cur)
!         call pbuf_get_field(pbuf_chnk, index_co2_lnd, co2_lnd_cur)
!         call pbuf_get_field(pbuf_chnk, index_co2    , co2_cur)
!         call pbuf_get_field(pbuf_chnk, index_acco2  , ac_CO2)

!         do i = 1, ncol
!            flux_ocn(i) = cam_in(ichnk)%cflx(i, c_i(1))
!            flux_fff(i) = cam_in(ichnk)%cflx(i, c_i(2)) + 
!            flux_lnd(i) = cam_in(ichnk)%cflx(i, c_i(3))
!            flux(i)     = flux_ocn(i) + flux_fff(i) + flux_lnd(i)
!         end do ! i = 1, ncol

!      end do ! ichnk = begchunk, endchunk


!      static_ener_ac_idx = pbuf_get_index('static_ener_ac')
!      call pbuf_get_field(pbuf, static_ener_ac_idx, static_ener_ac_2d )
!      water_vap_ac_idx   = pbuf_get_index('water_vap_ac')
!      call pbuf_get_field(pbuf, water_vap_ac_idx, water_vap_ac_2d )

      !Integrate column static energy
!      ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
!    do k=2,pver
!       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
!    end do
!    static_ener_ac_2d(:ncol) = ftem(:ncol,1)

!    static_ener_ac_idx = pbuf_get_index('static_ener_ac')
!    call pbuf_get_field(pbuf, static_ener_ac_idx, static_ener_ac_2d )
!    water_vap_ac_idx   = pbuf_get_index('water_vap_ac')
!    call pbuf_get_field(pbuf, water_vap_ac_idx, water_vap_ac_2d )

    ! Integrate and compute the difference
    ! CIDiff = difference of column integrated values
!    if( nstep == 0 ) then
!       CIDiff(:ncol) = 0.0_r8
!       call outfld('DTENDTH', CIDiff, pcols, lchnk )
!       call outfld('DTENDTQ', CIDiff, pcols, lchnk )
!    else
!       ! MSE first
!       ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
!       do k=2,pver
!          ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
!       end do
!       CIDiff(:ncol) = (ftem(:ncol,1) - static_ener_ac_2d(:ncol))*rtdt!

!       call outfld('DTENDTH', CIDiff, pcols, lchnk )
!       ! Water vapor second
!       ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
!       do k=2,pver
!          ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
!       end do
!       CIDiff(:ncol) = (ftem(:ncol,1) - water_vap_ac_2d(:ncol))*rtdt

!       call outfld('DTENDTQ', CIDiff, pcols, lchnk )
!    end if
      

!   end subroutine check_co2_change_nope

   subroutine co2_gmean_check_wflux(title, state, cam_in)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the 
! physics decomposition. Prints diagnostics to log file.
!
! Authors: B. Eaton (based on gavglook) & Bryce Harrop
!
!-----------------------------------------------------------------------
      use ppgrid,         only: pver, pcols, begchunk, endchunk
      use physconst,      only: gravit
      use phys_grid,      only: get_ncols_p
      use physics_types,  only: physics_state, physics_ptend
      use constituents,   only: pcnst, cnst_name
      use cam_logfile,    only: iulog
      use phys_gmean,     only: gmean
!
! Arguments
!
      character(len=*),    intent(in) :: title    ! location of this call
      type(physics_state), intent(in) :: state(begchunk:endchunk)
      type(cam_in_t),      intent(in) :: cam_in(begchunk:endchunk)
!
! Local workspace
!
      character(len=*), parameter :: sub_name='co2_gmean_check: '

      integer :: c, i, k, m
      integer :: ierr
      integer :: ncols

      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
      real(r8) :: mass_wet_mean(pcnst)     ! global mean constituent masses assuming moist mmr
      real(r8) :: mass_dry_mean(pcnst)     ! global mean constituent masses assuming dry mmr
      real(r8), pointer :: sfc_flux(:,:,:) ! constituent surface flux
      real(r8) :: sfc_flux_mean(pcnst)     ! global mean constituent surface flux

!
!-----------------------------------------------------------------------
!
      allocate(mass_wet(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_dry'

      allocate(sfc_flux(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate sfc_flux'

      do m = 1, pcnst
         do c = begchunk, endchunk
            ncols = get_ncols_p(c)
            do i = 1, ncols

               ! Compute column masses assuming both dry and wet mixing ratios

               mass_wet(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
                                    state(c)%pdel(i,k) * state(c)%q(i,k,m)
               end do
               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

               mass_dry(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
                                    state(c)%pdeldry(i,k) * state(c)%q(i,k,m)
               end do
               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

               sfc_flux(i,c,m) = 0.0_r8
               sfc_flux(i,c,m) = sfc_flux(i,c,m) + &
                                 cam_in(c)%cflx(i,m)

            end do
         end do
      end do

      ! compute global mean mass
      call gmean(mass_wet, mass_wet_mean, pcnst)
      call gmean(mass_dry, mass_dry_mean, pcnst)
      call gmean(sfc_flux, sfc_flux_mean, pcnst)

      ! report to log file
      if (masterproc) then

         do m = 1, ncnst
               write (6,66) trim(title)//' m=',c_i(m), &
                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc ', &
                  mass_dry_mean(c_i(m)), mass_wet_mean(c_i(m)), &
                  sfc_flux_mean(c_i(m))
66             format (a32,i2,a36,1p,3e25.17)
         end do

      endif

      deallocate(sfc_flux)
      deallocate(mass_dry)
      deallocate(mass_wet)

   end subroutine co2_gmean_check_wflux

   subroutine co2_gmean_check2_wflux(title, state, ptend, cam_in, pbuf)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the 
! physics decomposition. Prints diagnostics to log file.
!
! Authors: B. Eaton (based on gavglook) & Bryce Harrop
!
!-----------------------------------------------------------------------
      use ppgrid,         only: pver, pcols, begchunk, endchunk
      use physconst,      only: gravit
      use phys_grid,      only: get_ncols_p
      use physics_types,  only: physics_state, physics_ptend
      use constituents,   only: pcnst, cnst_name
      use cam_logfile,    only: iulog
      use phys_gmean,     only: gmean
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
!
! Arguments
!
      character(len=*),    intent(in) :: title    ! location of this call
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(physics_ptend), intent(in) :: ptend
      type(cam_in_t),      intent(in) :: cam_in
      !integer, dimension(:), intent(in) :: cnst_ind_arr
!
! Local workspace
!
      character(len=*), parameter :: sub_name='co2_gmean_check2: '

      integer :: c, i, k, m
      integer :: ierr
      integer :: ncol
      integer :: ifld, idx
      real(r8)          :: ac_CO2_tot(pcols)
      real(r8), pointer :: ac_CO2(:,:)

      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
      real(r8) :: mass_wet_mean(pcnst)     ! global mean constituent masses assuming moist mmr
      real(r8) :: mass_dry_mean(pcnst)     ! global mean constituent masses assuming dry mmr
      real(r8), pointer :: sfc_flux(:,:,:) ! constituent surface flux
      real(r8) :: sfc_flux_mean(pcnst)     ! global mean constituent surface flux
      real(r8), pointer :: air_flux(:,:,:) ! aircraft emission flux
      real(r8) :: air_flux_mean(pcnst)     ! global mean aircraft emission flux

!
!-----------------------------------------------------------------------
!
      allocate(mass_wet(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(sfc_flux(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate sfc_flux'

      allocate(air_flux(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate air_flux'

      ifld = pbuf_get_index('ac_CO2')   
      call pbuf_get_field(pbuf, ifld, ac_CO2)

      ! Set CO2 global indices
      do idx = 1, ncnst
         select case (trim(c_names(idx)))
         case ('CO2_OCN')
            co2_ocn_glo_ind = c_i(idx)
         case ('CO2_FFF')
            co2_fff_glo_ind = c_i(idx)
         case ('CO2_LND')
            co2_lnd_glo_ind = c_i(idx)
         case ('CO2')
            co2_glo_ind     = c_i(idx)
         end select
      end do

      ncol = state%ncol
      c    = state%lchnk
      !ncnst = size(cnst_ind_arr)
      do m = 1, pcnst
!         do c = begchunk, endchunk
!            ncols = get_ncols_p(c)
            do i = 1, ncol

               ! Compute column masses assuming both dry and wet mixing ratios

               mass_wet(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
                                    state%pdel(i,k)*state%q(i,k,m)
               end do
               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

               mass_dry(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
                                    state%pdeldry(i,k)*state%q(i,k,m)
               end do
               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

               sfc_flux(i,c,m) = 0.0_r8
               sfc_flux(i,c,m) = sfc_flux(i,c,m) + &
                                 cam_in%cflx(i,m)
               
               ac_CO2_tot(i) = 0.0_r8
               do k = 1, pver
                  ac_CO2_tot(i) = ac_CO2_tot(i) + ac_CO2(i,k)
               end do

            end do
!         end do
            if (m .eq. co2_glo_ind .or. m .eq. co2_fff_glo_ind) then
               air_flux(:ncol,c,m) = ac_CO2_tot(:ncol)
            else
               air_flux(:ncol,c,m) = 0.0_r8
            end if
      end do

      ! compute global mean mass
      call gmean(mass_wet, mass_wet_mean, pcnst)
      call gmean(mass_dry, mass_dry_mean, pcnst)
      call gmean(sfc_flux, sfc_flux_mean, pcnst)
      call gmean(air_flux, air_flux_mean, pcnst)

      ! report to log file
      if (masterproc) then

         do m = 1, ncnst
               write (6,66) trim(title)//' m=',c_i(m), &
                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc, air ', &
                  mass_dry_mean(c_i(m)), mass_wet_mean(c_i(m)), &
                  sfc_flux_mean(c_i(m)), air_flux_mean(c_i(m))
66             format (a32,i2,a36,1p,4e25.17)
!old version66             format (a32,i2,a36,1p,4e25.13)
         end do

      endif

      deallocate(air_flux)
      deallocate(sfc_flux)
      deallocate(mass_dry)
      deallocate(mass_wet)

   end subroutine co2_gmean_check2_wflux

end module co2_diagnostics
