module atm_import_export

  use shr_kind_mod  , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use cam_logfile      , only: iulog
  implicit none

#ifdef HAVE_MOAB
  ! to store all fields to be set in moab
  integer                :: mblsize, totalmbls, nsend, totalmbls_r, nrecv
  real(r8) , allocatable :: a2x_am(:,:) ! atm to coupler, on atm mesh, on atm component pes
  real(r8) , allocatable :: x2a_am(:,:) ! coupler to atm, on atm mesh, on atm component pes
#endif

contains

  subroutine atm_import( x2a, cam_in, restart_init)

    !-----------------------------------------------------------------------
    use cam_cpl_indices
    use camsrfexch,     only: cam_in_t
    use phys_grid ,     only: get_ncols_p
    use ppgrid    ,     only: begchunk, endchunk
    use shr_const_mod,  only: shr_const_stebol
    use seq_drydep_mod, only: n_drydep
    use co2_cycle     , only: c_i, co2_readFlux_ocn, co2_readFlux_fuel
    use co2_cycle     , only: co2_transport, co2_time_interp_ocn, co2_time_interp_fuel
    use co2_cycle     , only: data_flux_ocn, data_flux_fuel
    use iac_coupled_fields, only: iac_vertical_emiss
    use phys_control  , only: iac_present
    use physconst     , only: mwco2
    use time_manager  , only: is_first_step, get_curr_date
    use constituents  , only: pcnst
    use cam_abortutils, only: endrun
#ifdef HAVE_MOAB
    use seq_comm_mct,   only: mphaid
    use shr_kind_mod,   only: CXX=>shr_kind_cxx
    use iMOAB,          only: iMOAB_GetDoubleTagStorage
    use iso_c_binding,  only: C_NULL_CHAR
    use seq_flds_mod,   only: seq_flds_x2a_fields
#endif
    !
    ! Arguments
    !
    real(r8)      , intent(in)    :: x2a(:,:)
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)
    logical, optional, intent(in) :: restart_init
    !
    ! Local variables
    !
    integer            :: i,lat,n,c,ig  ! indices
    integer            :: ncols         ! number of columns
    logical, save      :: first_time = .true.
    integer, parameter :: ndst = 2
    integer, target    :: spc_ndx(ndst)
    integer, pointer   :: dst_a5_ndx, dst_a7_ndx
    integer, pointer   :: dst_a1_ndx, dst_a3_ndx
    integer :: icnst, mon_idx
    logical :: overwrite_flds
    integer :: idx_megan_end, idx_ddvel_end  ! end indices for array slices
#ifdef HAVE_MOAB
    character(CXX) :: tagname
    integer :: ent_type, ierr
#endif
    !-----------------------------------------------------------------------
    overwrite_flds = .true.
    ! don't overwrite fields if invoked during the initialization phase
    ! of a 'continue' or 'branch' run type with data from .rs file
    if (present(restart_init)) overwrite_flds = .not. restart_init

    if (iac_present) then
      mon_idx = get_month_index()
    endif

#ifdef HAVE_MOAB
    tagname=trim(seq_flds_x2a_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_GetDoubleTagStorage ( mphaid, tagname, totalmbls_r , ent_type, x2a_am )
    if ( ierr > 0) then
      call endrun('Error: fail to get  seq_flds_a2x_fields for atm physgrid moab mesh')
    endif
#endif

    ! E3SM sign convention is that fluxes are positive downward

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)

       ! initialize constituent surface fluxes to zero
       ! NOTE:overwrite_flds is .FALSE. for the first restart
       ! time step making cflx(:,1)=0.0 for the first restart time step.
       ! cflx(:,1) should not be zeroed out, start the second index of cflx from 2.

       ! +++ Update from 2022-09 +++
       ! For some of the new process coupling options in EAM, some of the constituents'
       ! cam_in%cflx are used not in tphysac but in the tphysbc call of the next time step.
       ! This means for an exact restart, we also need to write out cam_in%cflx(:,2:)
       ! and then read them back in. Because the present subroutine is called after
       ! the cflx variables are read in in the subroutine read_restart_physics
       ! in physics/cam/restart_physics.F90, we need to move the following line
       ! to that read_restart_physics to avoid incorrectly zeroing out the needed values.
       !
       !cam_in(c)%cflx(:,2:) = 0._r8
       !
       ! === Update from 2022-09 ===

       do i =1,ncols
          if (overwrite_flds) then
             ! Prior to this change, "overwrite_flds" was always .true. therefore wsx and wsy were always updated.
             ! Now, overwrite_flds is .false. for the first time step of the restart run. Move wsx and wsy out of
             ! this if-condition so that they are still updated everytime irrespective of the value of overwrite_flds.

             ! Move lhf to this if-block so that it is not overwritten to ensure BFB restarts when qneg4 correction
             ! occurs at the restart time step
             ! Modified by Wuyin Lin
             cam_in(c)%shf(i)    = -x2a_get(ig, index_x2a_Faxx_sen)
             cam_in(c)%cflx(i,1) = -x2a_get(ig, index_x2a_Faxx_evap)
             cam_in(c)%lhf(i)    = -x2a_get(ig, index_x2a_Faxx_lat)
          endif

          if (index_x2a_Faoo_h2otemp /= 0) then
             cam_in(c)%h2otemp(i) = -x2a_get(ig, index_x2a_Faoo_h2otemp)
          end if

          cam_in(c)%wsx(i)    = -x2a_get(ig, index_x2a_Faxx_taux)
          cam_in(c)%wsy(i)    = -x2a_get(ig, index_x2a_Faxx_tauy)
          cam_in(c)%lwup(i)      = -x2a_get(ig, index_x2a_Faxx_lwup)
          cam_in(c)%asdir(i)     =  x2a_get(ig, index_x2a_Sx_avsdr)
          cam_in(c)%aldir(i)     =  x2a_get(ig, index_x2a_Sx_anidr)
          cam_in(c)%asdif(i)     =  x2a_get(ig, index_x2a_Sx_avsdf)
          cam_in(c)%aldif(i)     =  x2a_get(ig, index_x2a_Sx_anidf)
          cam_in(c)%ts(i)        =  x2a_get(ig, index_x2a_Sx_t)
          cam_in(c)%sst(i)       =  x2a_get(ig, index_x2a_So_t)
          cam_in(c)%snowhland(i) =  x2a_get(ig, index_x2a_Sl_snowh)
          cam_in(c)%snowhice(i)  =  x2a_get(ig, index_x2a_Si_snowh)
          cam_in(c)%tref(i)      =  x2a_get(ig, index_x2a_Sx_tref)
          cam_in(c)%qref(i)      =  x2a_get(ig, index_x2a_Sx_qref)
          cam_in(c)%u10(i)       =  x2a_get(ig, index_x2a_Sx_u10)
          cam_in(c)%u10withgusts(i) = x2a_get(ig, index_x2a_Sx_u10withgusts)
          cam_in(c)%icefrac(i)   =  x2a_get(ig, index_x2a_Sf_ifrac)
          cam_in(c)%ocnfrac(i)   =  x2a_get(ig, index_x2a_Sf_ofrac)
          cam_in(c)%landfrac(i)  =  x2a_get(ig, index_x2a_Sf_lfrac)
          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  x2a_get(ig, index_x2a_Sl_ram1)
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  x2a_get(ig, index_x2a_Sl_fv)
          if ( associated(cam_in(c)%soilw) ) &
               cam_in(c)%soilw(i) =  x2a_get(ig, index_x2a_Sl_soilw)
          if ( associated(cam_in(c)%dstflx) ) then
             cam_in(c)%dstflx(i,1) = x2a_get(ig, index_x2a_Fall_flxdst1)
             cam_in(c)%dstflx(i,2) = x2a_get(ig, index_x2a_Fall_flxdst2)
             cam_in(c)%dstflx(i,3) = x2a_get(ig, index_x2a_Fall_flxdst3)
             cam_in(c)%dstflx(i,4) = x2a_get(ig, index_x2a_Fall_flxdst4)
          endif
          if ( associated(cam_in(c)%meganflx) ) then
             idx_megan_end = index_x2a_Fall_flxvoc+shr_megan_mechcomps_n-1
             cam_in(c)%meganflx(i,1:shr_megan_mechcomps_n) = &
                  x2a_slice(ig, index_x2a_Fall_flxvoc, idx_megan_end)
          endif

          ! dry dep velocities
          if ( index_x2a_Sl_ddvel/=0 .and. n_drydep>0 ) then
             idx_ddvel_end = index_x2a_Sl_ddvel+n_drydep-1
             cam_in(c)%depvel(i,:n_drydep) = &
                  x2a_slice(ig, index_x2a_Sl_ddvel, idx_ddvel_end)
          endif
          !
          ! fields needed to calculate water isotopes to ocean evaporation processes
          !
          cam_in(c)%ustar(i) = x2a_get(ig, index_x2a_So_ustar)
          cam_in(c)%re(i)    = x2a_get(ig, index_x2a_So_re)
          cam_in(c)%ssq(i)   = x2a_get(ig, index_x2a_So_ssq)
          !
          ! bgc scenarios
          !
          if (index_x2a_Fall_fco2_lnd /= 0) then
             cam_in(c)%fco2_lnd(i) = -x2a_get(ig, index_x2a_Fall_fco2_lnd)
          end if

          !------------------------------------------------------------------------------------------
          ! EHC fields do not need any interpolation: annual emissions were split assuming
          ! the monthly flux values were applied to the seconds in each month
          ! This is true for CEDS data as well. The last time step of the year is labelled as 
          ! nxty0101-00000, and needs month 12
          !------------------------------------------------------------------------------------------
          if (iac_present) then
             ! if surface emissions from EHC exist for this month, get them from coupler var
             if (index_x2a_Fazz_co2sfc_iac(mon_idx) /= 0) then
                cam_in(c)%fco2_surface_iac(i) = -x2a_get(ig, index_x2a_Fazz_co2sfc_iac(mon_idx))
             end if
             ! if aircraft lo emissions from EHC exist for this month, get them from coupler var
             if (index_x2a_Fazz_co2airlo_iac(mon_idx) /= 0) then
                iac_vertical_emiss(c)%fco2_low_height(i) = -x2a_get(ig, index_x2a_Fazz_co2airlo_iac(mon_idx))
             end if
             ! if aircraft lo emissions from EHC exist for this month, get them from coupler var
             if (index_x2a_Fazz_co2airhi_iac(mon_idx) /= 0) then
                iac_vertical_emiss(c)%fco2_high_height(i) = -x2a_get(ig, index_x2a_Fazz_co2airhi_iac(mon_idx))
             endif
          endif ! if (iac_present)

          if (index_x2a_Faoo_fco2_ocn /= 0) then
             cam_in(c)%fco2_ocn(i) = -x2a_get(ig, index_x2a_Faoo_fco2_ocn)
          end if
          if (index_x2a_Faoo_fdms_ocn /= 0) then
             cam_in(c)%fdms(i)     = -x2a_get(ig, index_x2a_Faoo_fdms_ocn)
          end if

          ig=ig+1

       end do
    end do

    ! Get total co2 flux from components,
    ! Note - co2_transport determines if cam_in(c)%cflx(i,c_i(1:4)) is allocated

    if (co2_transport().and.overwrite_flds) then

       ! Interpolate in time for flux data read in
       if (co2_readFlux_ocn) then
          call co2_time_interp_ocn
       end if
       if (co2_readFlux_fuel) then
          call co2_time_interp_fuel
       end if

       ! from ocn : data read in or from coupler or zero
       ! from fuel: data read in or zero
       ! from lnd : through coupler or zero
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)
          do i=1,ncols

             ! all co2 fluxes in unit kgCO2/m2/s ! co2 flux from ocn
             if (index_x2a_Faoo_fco2_ocn /= 0) then
                !FIXMEB: Instead of using hardwired numbers, 1,2 etc, can't we use integer parameters for c_i indices?
                cam_in(c)%cflx(i,c_i(1)) = cam_in(c)%fco2_ocn(i)
             else if (co2_readFlux_ocn) then
                ! convert from molesCO2/m2/s to kgCO2/m2/s
! The below section involves a temporary workaround for fluxes from data (read in from a file)
! There is an issue with infld that does not allow time-varying 2D files to be read correctly.
! The work around involves adding a singleton 3rd dimension offline and reading the files as
! 3D fields.  Once this issue is corrected, the old implementation can be reinstated.
! This is the case for both data_flux_ocn and data_flux_fuel
!++BEH  vvv old implementation vvv
!                cam_in(c)%cflx(i,c_i(1)) = &
!                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i)) &
!                     *mwco2*1.0e-3_r8
!       ^^^ old implementation ^^^   ///    vvv new implementation vvv
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,1,c)*(1._r8- cam_in(c)%landfrac(i)) &
                     *mwco2*1.0e-3_r8
!--BEH  ^^^ new implementation ^^^
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if

             ! co2 flux from fossil fuel
             if ( iac_present ) then
               if( index_x2a_Fazz_co2sfc_iac(mon_idx) /= 0) then
                  cam_in(c)%cflx(i,c_i(2)) = cam_in(c)%fco2_surface_iac(i)
               end if
             else if (co2_readFlux_fuel) then
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,1,c)
!--BEH  ^^^ new implementation ^^^
             else
                cam_in(c)%cflx(i,c_i(2)) = 0._r8
             end if

             ! co2 flux from land (cpl already multiplies flux by land fraction)
             if (index_x2a_Fall_fco2_lnd /= 0) then
                cam_in(c)%cflx(i,c_i(3)) = cam_in(c)%fco2_lnd(i)
             else
                cam_in(c)%cflx(i,c_i(3)) = 0._r8
             end if

             ! merged co2 flux
             cam_in(c)%cflx(i,c_i(4)) = cam_in(c)%cflx(i,c_i(1)) + &
                                        cam_in(c)%cflx(i,c_i(2)) + &
                                        cam_in(c)%cflx(i,c_i(3))
          end do
       end do
    end if
    !
    ! if first step, determine longwave up flux from the surface temperature
    !
    if (first_time) then
       if (is_first_step()) then
          do c=begchunk, endchunk
             ncols = get_ncols_p(c)
             do i=1,ncols
                cam_in(c)%lwup(i) = shr_const_stebol*(cam_in(c)%ts(i)**4)
             end do
          end do
       end if
       first_time = .false.
    end if

  contains

    pure function x2a_get(ig, index) result(val)
      integer, intent(in) :: ig
      integer, intent(in) :: index
      real(r8) :: val
#ifdef HAVE_MOAB
      val = x2a_am(ig, index)
#else
      val = x2a(index, ig)
#endif
    end function x2a_get

    pure function x2a_slice(ig, istart, iend) result(val)
      integer, intent(in) :: ig
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      real(r8) :: val(iend - istart + 1)
#ifdef HAVE_MOAB
      val = x2a_am(ig, istart:iend)
#else
      val = x2a(istart:iend, ig)
#endif
    end function x2a_slice

  end subroutine atm_import

  !===============================================================================

   function get_month_index() result(mon_idx)
   
      !-----------------------------------------------------------------------
      ! Determine which month's EHC data to use based on atm time manager date
      ! This function assumes that the atm time manager is always in sync with
      ! the EAM tm clock during atm_import calls
      !-----------------------------------------------------------------------
      use time_manager  , only: get_curr_date, is_first_step
      use shr_log_mod   , only: errMsg => shr_log_errMsg
      use cam_abortutils, only: endrun
      implicit none

      integer :: mon_idx !return value

      !local variables
      integer, parameter :: FIRST_MONTH = 1
      integer, parameter :: LAST_MONTH = 12
      integer :: yr, mon, day, tod
      character(len=256) :: errstr
      !---------------------------------------------------------------------
      ! Get the current model date
      !---------------------------------------------------------------------
      call get_curr_date( yr, mon, day, tod )
    
      !Sanity check for month (this shouldn't happen)
      if (mon < FIRST_MONTH .or. mon > LAST_MONTH) then
         write(errstr,*) 'ERROR! Month is out of bounds [',FIRST_MONTH,',', LAST_MONTH,'], current month is: ', mon,'. '
         call endrun(trim(errstr)//errMsg(__FILE__, __LINE__))
      end if
      !-----------------------------------------------------------------------
      ! Determine month index for emissions lookup
      !-----------------------------------------------------------------------
      ! EAM tm clock should always match atm Eclock/sync clock during
      ! atm_import calls
      if (mon == FIRST_MONTH .and. day == 1 .and. tod == 0) then
         if (is_first_step()) then
            ! Get this year's data because this timestep is run at the start of the year
            ! for the initial model start, the actual month 12 data for the previous
            ! data are not available
            mon_idx = FIRST_MONTH
         else
            ! Timestep usually run at the end of the year
            mon_idx = LAST_MONTH
         end if
      else if (day == 1 .and. tod == 0) then
         ! Last timestep of the previous month
         mon_idx = mon - 1
      else
         mon_idx = mon
      end if ! if (mon == 1 .and. day == 1 .and. tod == 0)

      !Sanity check for mon_idx
      if (mon_idx < FIRST_MONTH .or. mon_idx > LAST_MONTH) then
         write(errstr,*) 'ERROR! mon_idx is out of bounds [',FIRST_MONTH,',', LAST_MONTH,'], mon_idx is: ', mon_idx,'. '
         call endrun(trim(errstr)//errMsg(__FILE__, __LINE__))
      end if
   end function get_month_index

   !===============================================================================

  subroutine atm_export( cam_out, a2x )

    !-------------------------------------------------------------------
    use camsrfexch, only: cam_out_t
    use phys_grid , only: get_ncols_p
    use ppgrid    , only: begchunk, endchunk
    use cam_cpl_indices
    use phys_control, only: phys_getopts
    use lnd_infodata, only: precip_downscaling_method
    use cam_abortutils, only: endrun
#ifdef HAVE_MOAB
    use shr_kind_mod,  only: CXX=>shr_kind_cxx
    use seq_comm_mct,  only: mphaid
    use iMOAB,         only: iMOAB_WriteMesh, iMOAB_SetDoubleTagStorage
    use iso_c_binding, only: C_NULL_CHAR
    use seq_flds_mod,  only: seq_flds_a2x_fields
#endif

#ifdef MOABDEBUG
    character*100 outfile, wopts, lnum
    integer, save :: local_count = 0
    character*100 lnum2
#endif
    !
    ! Arguments
    !
    type(cam_out_t), intent(in)    :: cam_out(begchunk:endchunk)
    real(r8)       , intent(inout) :: a2x(:,:)
    !
    ! Local variables
    !
    integer :: i,c,n,ig         ! indices
    integer :: ncols            ! Number of columns
    logical :: linearize_pbl_winds
#ifdef HAVE_MOAB
    character(CXX) :: tagname
    integer :: ent_type, ierr
#endif
    !-----------------------------------------------------------------------

    call phys_getopts(linearize_pbl_winds_out=linearize_pbl_winds)

#ifdef HAVE_MOAB
    ! initialize/reset all export data to zero
    a2x_am = 0.0D0
#endif

    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector
    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          call a2x_set(ig, index_a2x_Sa_pslv, cam_out(c)%psl(i))
          call a2x_set(ig, index_a2x_Sa_z, cam_out(c)%zbot(i))
          call a2x_set(ig, index_a2x_Sa_u, cam_out(c)%ubot(i))
          call a2x_set(ig, index_a2x_Sa_v, cam_out(c)%vbot(i))
          if (linearize_pbl_winds) then
             call a2x_set(ig, index_a2x_Sa_wsresp, cam_out(c)%wsresp(i))
             call a2x_set(ig, index_a2x_Sa_tau_est, cam_out(c)%tau_est(i))
          end if
          ! This check is only for SCREAMv0; otherwise gustiness should always
          ! be exported.
          if (index_a2x_Sa_ugust /= 0) then
             call a2x_set(ig, index_a2x_Sa_ugust, cam_out(c)%ugust(i))
          end if
          call a2x_set(ig, index_a2x_Sa_tbot, cam_out(c)%tbot(i))
          call a2x_set(ig, index_a2x_Sa_ptem, cam_out(c)%thbot(i))
          call a2x_set(ig, index_a2x_Sa_pbot, cam_out(c)%pbot(i))
          call a2x_set(ig, index_a2x_Sa_shum, cam_out(c)%qbot(i,1))
          call a2x_set(ig, index_a2x_Sa_dens, cam_out(c)%rho(i))

          if (trim(adjustl(precip_downscaling_method)) == "FNM") then
             !if the land model's precip downscaling method is FNM, export uovern to the coupler
             call a2x_set(ig, index_a2x_Sa_uovern, cam_out(c)%uovern(i))
          end if
          call a2x_set(ig, index_a2x_Faxa_swnet, cam_out(c)%netsw(i))
          call a2x_set(ig, index_a2x_Faxa_lwdn, cam_out(c)%flwds(i))
          call a2x_set(ig, index_a2x_Faxa_rainc, (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8)
          call a2x_set(ig, index_a2x_Faxa_rainl, (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8)
          call a2x_set(ig, index_a2x_Faxa_snowc, cam_out(c)%precsc(i)*1000._r8)
          call a2x_set(ig, index_a2x_Faxa_snowl, cam_out(c)%precsl(i)*1000._r8)
          call a2x_set(ig, index_a2x_Faxa_swndr, cam_out(c)%soll(i))
          call a2x_set(ig, index_a2x_Faxa_swvdr, cam_out(c)%sols(i))
          call a2x_set(ig, index_a2x_Faxa_swndf, cam_out(c)%solld(i))
          call a2x_set(ig, index_a2x_Faxa_swvdf, cam_out(c)%solsd(i))

          ! aerosol deposition fluxes
          call a2x_set(ig, index_a2x_Faxa_bcphidry, cam_out(c)%bcphidry(i))
          call a2x_set(ig, index_a2x_Faxa_bcphodry, cam_out(c)%bcphodry(i))
          call a2x_set(ig, index_a2x_Faxa_bcphiwet, cam_out(c)%bcphiwet(i))
          call a2x_set(ig, index_a2x_Faxa_ocphidry, cam_out(c)%ocphidry(i))
          call a2x_set(ig, index_a2x_Faxa_ocphodry, cam_out(c)%ocphodry(i))
          call a2x_set(ig, index_a2x_Faxa_ocphiwet, cam_out(c)%ocphiwet(i))
          call a2x_set(ig, index_a2x_Faxa_dstwet1, cam_out(c)%dstwet1(i))
          call a2x_set(ig, index_a2x_Faxa_dstdry1, cam_out(c)%dstdry1(i))
          call a2x_set(ig, index_a2x_Faxa_dstwet2, cam_out(c)%dstwet2(i))
          call a2x_set(ig, index_a2x_Faxa_dstdry2, cam_out(c)%dstdry2(i))
          call a2x_set(ig, index_a2x_Faxa_dstwet3, cam_out(c)%dstwet3(i))
          call a2x_set(ig, index_a2x_Faxa_dstdry3, cam_out(c)%dstdry3(i))
          call a2x_set(ig, index_a2x_Faxa_dstwet4, cam_out(c)%dstwet4(i))
          call a2x_set(ig, index_a2x_Faxa_dstdry4, cam_out(c)%dstdry4(i))

          if (index_a2x_Sa_co2prog /= 0) then
             call a2x_set(ig, index_a2x_Sa_co2prog, cam_out(c)%co2prog(i)) ! atm prognostic co2
          end if
          if (index_a2x_Sa_co2diag /= 0) then
             call a2x_set(ig, index_a2x_Sa_co2diag, cam_out(c)%co2diag(i)) ! atm diagnostic co2
          end if

          ig=ig+1
       end do
    end do
#ifdef HAVE_MOAB
    tagname=trim(seq_flds_a2x_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, totalmbls , ent_type, a2x_am )
    if ( ierr > 0) then
      call endrun('Error: fail to set  seq_flds_a2x_fields for atm physgrid moab mesh')
    endif
#endif
#ifdef MOABDEBUG
    write(lnum,"(I0.2)")cur_atm_stepno
    local_count = local_count + 1
    write(lnum2,"(I0.2)")local_count
    outfile = 'atm_export_'//trim(lnum)//'_'//trim(lnum2)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mphaid, outfile, wopts)
    if (ierr > 0 )  &
      call endrun('Error: fail to write the atm phys mesh file with data')
#endif

  contains

    subroutine a2x_set(ig, index, val)
      integer, intent(in) :: ig
      integer, intent(in) :: index
      real(r8), intent(in) :: val
#ifdef HAVE_MOAB
      a2x_am(ig, index) = val
#else
      a2x(index, ig) = val
#endif
    end subroutine a2x_set

  end subroutine atm_export

end module atm_import_export
