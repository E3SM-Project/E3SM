module atm_import_export

  use shr_kind_mod  , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  implicit none

contains

  subroutine atm_import( x2a, cam_in, restart_init )

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
    use physconst     , only: mwco2
    use time_manager  , only: is_first_step
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
    logical :: overwrite_flds
    !-----------------------------------------------------------------------
    overwrite_flds = .true.
    ! don't overwrite fields if invoked during the initialization phase 
    ! of a 'continue' or 'branch' run type with data from .rs file
    if (present(restart_init)) overwrite_flds = .not. restart_init

    ! ccsm sign convention is that fluxes are positive downward

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c) 

       ! initialize constituent surface fluxes to zero
       ! NOTE:overwrite_flds is .FALSE. for the first restart
       ! time step making cflx(:,1)=0.0 for the first restart time step.
       ! cflx(:,1) should not be zeroed out, start the second index of cflx from 2.
       cam_in(c)%cflx(:,2:) = 0._r8 
                                               
       do i =1,ncols                                                               
          if (overwrite_flds) then
             ! Prior to this change, "overwrite_flds" was always .true. therefore wsx and wsy were always updated.
             ! Now, overwrite_flds is .false. for the first time step of the restart run. Move wsx and wsy out of 
             ! this if-condition so that they are still updated everytime irrespective of the value of overwrite_flds.

             ! Move lhf to this if-block so that it is not overwritten to ensure BFB restarts when qneg4 correction 
             ! occurs at the restart time step
             ! Modified by Wuyin Lin
             cam_in(c)%shf(i)    = -x2a(index_x2a_Faxx_sen, ig)     
             cam_in(c)%cflx(i,1) = -x2a(index_x2a_Faxx_evap,ig)                
             cam_in(c)%lhf(i)    = -x2a(index_x2a_Faxx_lat, ig)     
          endif
          cam_in(c)%wsx(i)    = -x2a(index_x2a_Faxx_taux,ig)     
          cam_in(c)%wsy(i)    = -x2a(index_x2a_Faxx_tauy,ig)     
          cam_in(c)%lwup(i)      = -x2a(index_x2a_Faxx_lwup,ig)    
          cam_in(c)%asdir(i)     =  x2a(index_x2a_Sx_avsdr, ig)  
          cam_in(c)%aldir(i)     =  x2a(index_x2a_Sx_anidr, ig)  
          cam_in(c)%asdif(i)     =  x2a(index_x2a_Sx_avsdf, ig)  
          cam_in(c)%aldif(i)     =  x2a(index_x2a_Sx_anidf, ig)
          cam_in(c)%ts(i)        =  x2a(index_x2a_Sx_t,     ig)  
          cam_in(c)%sst(i)       =  x2a(index_x2a_So_t,     ig)             
          cam_in(c)%snowhland(i) =  x2a(index_x2a_Sl_snowh, ig)  
          cam_in(c)%snowhice(i)  =  x2a(index_x2a_Si_snowh, ig)  
          cam_in(c)%tref(i)      =  x2a(index_x2a_Sx_tref,  ig)  
          cam_in(c)%qref(i)      =  x2a(index_x2a_Sx_qref,  ig)
          cam_in(c)%u10(i)       =  x2a(index_x2a_Sx_u10,   ig)
          cam_in(c)%icefrac(i)   =  x2a(index_x2a_Sf_ifrac, ig)  
          cam_in(c)%ocnfrac(i)   =  x2a(index_x2a_Sf_ofrac, ig)
	  cam_in(c)%landfrac(i)  =  x2a(index_x2a_Sf_lfrac, ig)
          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  x2a(index_x2a_Sl_ram1 , ig)
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  x2a(index_x2a_Sl_fv   , ig)
          if ( associated(cam_in(c)%soilw) ) &
               cam_in(c)%soilw(i) =  x2a(index_x2a_Sl_soilw, ig)
          if ( associated(cam_in(c)%dstflx) ) then
             cam_in(c)%dstflx(i,1) = x2a(index_x2a_Fall_flxdst1, ig)
             cam_in(c)%dstflx(i,2) = x2a(index_x2a_Fall_flxdst2, ig)
             cam_in(c)%dstflx(i,3) = x2a(index_x2a_Fall_flxdst3, ig)
             cam_in(c)%dstflx(i,4) = x2a(index_x2a_Fall_flxdst4, ig)
          endif
          if ( associated(cam_in(c)%meganflx) ) then
             cam_in(c)%meganflx(i,1:shr_megan_mechcomps_n) = &
                  x2a(index_x2a_Fall_flxvoc:index_x2a_Fall_flxvoc+shr_megan_mechcomps_n-1, ig)
          endif

          ! dry dep velocities
          if ( index_x2a_Sl_ddvel/=0 .and. n_drydep>0 ) then
             cam_in(c)%depvel(i,:n_drydep) = &
                  x2a(index_x2a_Sl_ddvel:index_x2a_Sl_ddvel+n_drydep-1, ig)
          endif
          !
          ! fields needed to calculate water isotopes to ocean evaporation processes
          !
          cam_in(c)%ustar(i) = x2a(index_x2a_So_ustar,ig)
          cam_in(c)%re(i)    = x2a(index_x2a_So_re   ,ig)
          cam_in(c)%ssq(i)   = x2a(index_x2a_So_ssq  ,ig)
          !
          ! bgc scenarios
          !
          if (index_x2a_Fall_fco2_lnd /= 0) then
             cam_in(c)%fco2_lnd(i) = -x2a(index_x2a_Fall_fco2_lnd,ig)
          end if
          if (index_x2a_Faoo_fco2_ocn /= 0) then
             cam_in(c)%fco2_ocn(i) = -x2a(index_x2a_Faoo_fco2_ocn,ig)
          end if
          if (index_x2a_Faoo_fdms_ocn /= 0) then
             cam_in(c)%fdms(i)     = -x2a(index_x2a_Faoo_fdms_ocn,ig)
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
                cam_in(c)%cflx(i,c_i(1)) = cam_in(c)%fco2_ocn(i)
             else if (co2_readFlux_ocn) then 
                ! convert from molesCO2/m2/s to kgCO2/m2/s
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i)) &
                     *mwco2*1.0e-3_r8
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if
             
             ! co2 flux from fossil fuel
             if (co2_readFlux_fuel) then
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,c)
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

  end subroutine atm_import

  !===============================================================================

  subroutine atm_export( cam_out, a2x )

    !-------------------------------------------------------------------
    use camsrfexch, only: cam_out_t
    use phys_grid , only: get_ncols_p
    use ppgrid    , only: begchunk, endchunk       
    use cam_cpl_indices
    !
    ! Arguments
    !
    type(cam_out_t), intent(in)    :: cam_out(begchunk:endchunk) 
    real(r8)       , intent(inout) :: a2x(:,:)
    !
    ! Local variables
    !
    integer :: avsize, avnat
    integer :: i,m,c,n,ig       ! indices
    integer :: ncols            ! Number of columns
    !-----------------------------------------------------------------------

    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          a2x(index_a2x_Sa_pslv   ,ig) = cam_out(c)%psl(i)
          a2x(index_a2x_Sa_z      ,ig) = cam_out(c)%zbot(i)   
          a2x(index_a2x_Sa_u      ,ig) = cam_out(c)%ubot(i)   
          a2x(index_a2x_Sa_v      ,ig) = cam_out(c)%vbot(i)   
          a2x(index_a2x_Sa_tbot   ,ig) = cam_out(c)%tbot(i)   
          a2x(index_a2x_Sa_ptem   ,ig) = cam_out(c)%thbot(i)  
          a2x(index_a2x_Sa_pbot   ,ig) = cam_out(c)%pbot(i)   
          a2x(index_a2x_Sa_shum   ,ig) = cam_out(c)%qbot(i,1) 
	  a2x(index_a2x_Sa_dens   ,ig) = cam_out(c)%rho(i)
          a2x(index_a2x_Faxa_swnet,ig) = cam_out(c)%netsw(i)      
          a2x(index_a2x_Faxa_lwdn ,ig) = cam_out(c)%flwds(i)  
          a2x(index_a2x_Faxa_rainc,ig) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
          a2x(index_a2x_Faxa_rainl,ig) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
          a2x(index_a2x_Faxa_snowc,ig) = cam_out(c)%precsc(i)*1000._r8
          a2x(index_a2x_Faxa_snowl,ig) = cam_out(c)%precsl(i)*1000._r8
          a2x(index_a2x_Faxa_swndr,ig) = cam_out(c)%soll(i)   
          a2x(index_a2x_Faxa_swvdr,ig) = cam_out(c)%sols(i)   
          a2x(index_a2x_Faxa_swndf,ig) = cam_out(c)%solld(i)  
          a2x(index_a2x_Faxa_swvdf,ig) = cam_out(c)%solsd(i)  

          ! aerosol deposition fluxes
          a2x(index_a2x_Faxa_bcphidry,ig) = cam_out(c)%bcphidry(i)
          a2x(index_a2x_Faxa_bcphodry,ig) = cam_out(c)%bcphodry(i)
          a2x(index_a2x_Faxa_bcphiwet,ig) = cam_out(c)%bcphiwet(i)
          a2x(index_a2x_Faxa_ocphidry,ig) = cam_out(c)%ocphidry(i)
          a2x(index_a2x_Faxa_ocphodry,ig) = cam_out(c)%ocphodry(i)
          a2x(index_a2x_Faxa_ocphiwet,ig) = cam_out(c)%ocphiwet(i)
          a2x(index_a2x_Faxa_dstwet1,ig)  = cam_out(c)%dstwet1(i)
          a2x(index_a2x_Faxa_dstdry1,ig)  = cam_out(c)%dstdry1(i)
          a2x(index_a2x_Faxa_dstwet2,ig)  = cam_out(c)%dstwet2(i)
          a2x(index_a2x_Faxa_dstdry2,ig)  = cam_out(c)%dstdry2(i)
          a2x(index_a2x_Faxa_dstwet3,ig)  = cam_out(c)%dstwet3(i)
          a2x(index_a2x_Faxa_dstdry3,ig)  = cam_out(c)%dstdry3(i)
          a2x(index_a2x_Faxa_dstwet4,ig)  = cam_out(c)%dstwet4(i)
          a2x(index_a2x_Faxa_dstdry4,ig)  = cam_out(c)%dstdry4(i)

          if (index_a2x_Sa_co2prog /= 0) then
             a2x(index_a2x_Sa_co2prog,ig) = cam_out(c)%co2prog(i) ! atm prognostic co2
          end if
          if (index_a2x_Sa_co2diag /= 0) then
             a2x(index_a2x_Sa_co2diag,ig) = cam_out(c)%co2diag(i) ! atm diagnostic co2
          end if

          ig=ig+1
       end do
    end do
    
  end subroutine atm_export 

end module atm_import_export
