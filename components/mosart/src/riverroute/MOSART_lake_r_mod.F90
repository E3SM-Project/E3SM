!
MODULE mosart_lake_r_mod
! Description: core code of lake module in MOSART framework
! A simple wier formula (Kindsvater and Carter (1959))is used to estimate outflow from lakes 
! Developed by Wondie Yigzaw, July 2020. 
! REVISION HISTORY:
! 
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_PI
    use shr_sys_mod   , only : shr_sys_abort
    use RtmTimeManager
    use RunoffMod , only : Tctl, TUnit, TRunoff, THeat, TUnit_lake_r, TLake_r, TPara, rtmCTL
    use rof_cpl_indices, only : nt_nliq, nt_nice
    use RtmVar         , only : iulog, ngeom, nlayers, rstraflag, lakeflag
    use MOSART_heat_mod
    use MOSART_lake_hydro_mod
    use MOSART_lake_heat_mod
    use MOSART_lake_geometry_mod
    implicit none
    
    real(r8), parameter :: TINYVALUE  = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits

    public mosart_lake_r
    
! !PUBLIC MEMBER FUNCTIONS:
    contains  
        
    subroutine mosart_lake_r(iunit,nt,localDeltaT,temp_Tr)
    ! !DESCRIPTION: calculate the water temperature of reservoir.
    
        use shr_sys_mod , only : shr_sys_flush
    
        implicit none
        integer,  intent(in) :: iunit, nt
        real(r8), intent(in) :: localDeltaT,temp_Tr
        character(len=*),parameter :: subname = '(mosart_lake_r)'
        integer :: j,i,w,k,ii,l,ww,m,n,nn,mm,jmax,jmin                          ! indices

        real(r8) :: c_w = 4.188e3                                               ! Water specific heat capacity  (J/kg/k)
        real(r8) :: F = 0.8_r8                                                  ! dimensionless factor for wind sheltering by riparian vegetation, Wu et al, 2012
        real(r8) :: grav   = 9.8062_r8                                          ! gravity constant (m/s2)
        real(r8) :: le= 2.501e6                                                 ! latent heat of vaporaization (J/kg)
        real(r8) :: rho_w = 1.e3                                                ! Water density  (kg/m3)
        real(r8) :: rho_a = 1.177_r8                                            ! Air density  (kg/m3)
        real(r8) :: st_bl = 5.67e-8                                             ! Stefan-Boltzmann constant ~ W/m^2/K^4
        real(r8) :: t_frz   = 273.15_r8                                         ! freezing temperature (K)
        real(r8) :: sar       = 1.0_r8                                          ! Surface area ratio
        integer  :: s_dtime                                                     ! number of sub hourly time step
        integer  :: d_n_n                                                       ! Adjusted layer number
        integer  :: outlet                                                      ! Outlet layer location
        integer  :: mix1,mix2                                                   ! used in convective mixing
        real(r8) :: dtime                                                       ! time step (sec)
        real(r8) :: es,ea                                                       ! Satuated, atmospheric vapor pressure(Pa)
        real(r8) :: evap                                                        ! evaporation rate (m/s)
        real(r8) :: kl                                                          ! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
        real(r8) :: d_res_sub                                                   ! Reservoir depth, taken as 0.95*(dam height) in GRanD database, (m)
        real(r8) :: df_eff(nlayers+1)                                           ! Effective diffusivity (molecular + eddy) [m^2/s]
        real(r8) :: rho_z(nlayers),rho_r                                        ! Depth based water density, density of inflow water   (kg/m^3)
        real(r8) :: drhodz(nlayers)                                             ! d [rhow] /dz (kg/m**4)
        real(r8) :: m1(nlayers),m2(nlayers),m3(nlayers)                         ! used in tridiagonal matrix
        real(r8) :: Fx(nlayers)                                                 ! used in tridiagonal matrix
        real(r8) :: a(nlayers),b(nlayers),c(nlayers)                            ! columns of tridiagonal matrix
        real(r8) :: r(nlayers)                                                  ! right  diagonal tridiagonal matrix
        real(r8) :: ri                                                          ! Richardson number
        real(r8) :: q_in,t_in,q_ou                                              ! inflow (m^3/s),inflow temperature (k), and outflow (m^3/s)
        real(r8) :: dv_in(nlayers),dv_ou(nlayers)                               ! volume increment/decrease at layer due to inflow/outflow(m^3)
        real(r8) :: d_evap,v_evap                                               ! Evaporated depth (m) and volume (m^3)
        real(r8) :: d_prcp,v_prcp                                               ! Evaporated depth (m) and volume (m^3)
        real(r8) :: v_spill                                                     ! additional flow volume when storage capacity exceeded (m^3)
        real(r8) :: delta_z                                                     ! depth change to calculate corresponding area/volume(m)
        real(r8) :: top_d                                                       ! top layer depth after mixing
        real(r8) :: v_zt0(nlayers+1)                                            ! Total reservoir volume at depth z from surface(m^3)
        real(r8) :: v_mix                                                       ! Total volume of mixed layer(m^3)
        real(r8) :: t_s,t_out                                                   ! reservoir surface and outflow temperature (k)
        real(r8) :: t_z(nlayers)                                                ! lake layer temperature
        real(r8) :: t_z_old(nlayers)                                            ! previous time lake layer temperature
        real(r8) :: bv_f                                                        ! brunt-vaisala frequency (/s**2)
        real(r8) :: eta                                                         ! light extinction coefficient
        real(r8) :: sh_net                                                      ! net short wave radiation
        real(r8) :: sh_mix                                                      ! net short wave radiation in mixed layer
        real(r8) :: beta                                                        ! shortwave absorbtion factor
        real(r8) :: lw_abr                                                      ! net atmospheric longwave absorbtion (W/m^2)
        real(r8) :: phi_o                                                       ! net surface radiation (W/m^2)
        real(r8) :: phi_z(nlayers)                                              ! radiation absorbed by layer (W)
        real(r8) :: phi_x(nlayers+1)                                            ! radiation absorbed by mixed layer (W/m^2)
        real(r8) :: lw_ems,lt_heat,sn_heat                                      ! longwave emission,latent heat,sensible heat (W/m^2)
        real(r8) :: alb_s                                                       ! surface albedo shortwave radiation
        real(r8) :: k_m                                                         ! molecular diffusivity
        real(r8) :: enr_in(nlayers),enr_ou(nlayers)                             ! Layer Energy from inflow,outflow (J/s,W)
        real(r8) :: enr_0(nlayers),enr_1(nlayers),enr_2(nlayers)                ! Initial inner energy,Inner energy after advection,Inner energy after stratification(J/s,W)
        real(r8) :: enr_phi(nlayers)                                            ! Inner energy from solar/atmospheric radiation(J/s,W)
        real(r8) :: enr_err1,enr_err2                                           ! Energy error (w) before stratification and after triadiagonal solution
        real(r8) :: th_en(nlayers)                                              ! Layer thermal energy (j/s)
        real(r8) :: e_a,e_b,e_ab                                                ! used in layer merging/split, factors
        real(r8) :: k_ew,k_ad(nlayers)                                          ! Effective,wind,convection,advective kinetic energy (kg.m^2/s^2)
        real(r8) :: c_d                                                         ! Drag coefficient
        real(r8) :: cfa,cfw                                                     ! used in diffusion calculation
        real(r8) :: Fr(nlayers)                                                 ! Froude number squared and inverted for diffusion coeff. calculation
        real(r8) :: dis_ad(nlayers)                                             ! rate of dissipation-inflow/outflow
        real(r8) :: dis_w                                                       ! rate of dissipation-wind
        real(r8) :: l_vel, s_vel                                                ! used in Froude number calculation
        real(r8) :: tau                                                         ! shear stress
        real(r8) :: q_adv(nlayers)                                              ! Layer flow rate (m^3/s)
        real(r8) :: denmix,tmix,tsum                                            ! used in convective mixing
        real(r8) :: mixvol1,mixvol2,sumvol,num,dem                              ! used in convective mixing
        real(r8) :: ta,tb,tab,dd_za,dd_zb,d_va,dv_oua,delta_a                   ! used in layer merging/split
        real(r8) :: dv_inb,dv_ina,dv_oub,dv_ouab,dv_inab,dd_zab,d_vab,d_vb      ! used in layer merging/split    
        real(r8) :: ddz_min,ddz_max                                             ! (nd) Minimum and Maximum non-top layer thickness limit for merge/split
        real(r8) :: tmp_evap, tmp_outflow, tmp_tout                             ! dummy variables for evaporation rate (m3/s), outflow(m3/s), outflow temperature (K) 
        integer  :: dflag                                                       ! number of sub hourly time step
        real(r8) :: ddz_top_min,ddz_top_max                                     ! (nd) Minimum and Maximum top layer thickness limit for merge/split (m)
        real(r8) :: temp_prev                                                   ! lake surface temperature before advective mixing (Kelvin)
        real(r8) :: temp_equi                                                   ! lake surface temperature calculated using the euilibrium approach (Kelvin)
        real(r8) :: mytmp(14)                                                   ! inputs for calculating lake surface temperature calculated using the euilibrium approach
        integer  :: yr, mon, day, ymd, tod                                      ! time information
        
!**************************************************************************************************************************************************
        
        if (nt == nt_nliq) then    ! Lake module active if there is natural lake
            ! Initialize for all sub-timesteps                                    
            s_dtime = 10
            dtime   = localDeltaT/s_dtime        !    Sub-timestep for numerical purpose
            ddz_min = 5._r8
            ddz_max = max(3.0_r8*TLake_r%ddz_local(iunit), 3.0_r8*ddz_min)
            ddz_top_min = 5._r8 !1.5_r8 !ddz_min !1.5_r8
            ddz_top_max = 3._r8 * ddz_top_min !4.0_r8
            t_in    = temp_Tr 
            tmp_outflow = 0._r8
            tmp_evap = 0._r8
            tmp_tout = 0._r8
            !TLake_r%lake_rain(iunit) = THeat%forc_rain(iunit) * TUnit_lake_r%A_max(iunit) * 1e6  ! note the unit for A_max is km2
            !TLake_r%lake_snow(iunit) = THeat%forc_snow(iunit) * TUnit_lake_r%A_max(iunit) * 1e6
            TLake_r%lake_inflow(iunit) = -TRunoff%erout(iunit,nt) !inflow does not change during sub timesteps                
            TLake_r%lake_Tsur(iunit) = THeat%Tr(iunit)
            do ww = 1,s_dtime    ! Start calculation for each sub-timestep    ************************************************            
                
                !!======================== water balance processes   
    !if(iunit == 198841) then 
    !    write(unit=21002,fmt="(i4, i4, 6(e18.10))") 1, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), t_in, THeat%forc_t(iunit), TLake_r%temp_lake(iunit,1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit))
    !end if

    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 0, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if

                ! initialize for current sub time step
                q_in = 0._r8
                q_ou = 0._r8
                !TLake_r%lake_inflow(iunit)  = 0._r8
                TLake_r%lake_outflow(iunit) = 0._r8
                TLake_r%lake_spillflow(iunit) = 0._r8
				v_spill = 0._r8
                d_prcp = 0._r8
                v_prcp = 0._r8
                evap = 0._r8
                d_evap = 0._r8
                v_evap = 0._r8
                d_n_n          = TLake_r%d_ns(iunit)                    ! carry number of layers if there is merging/split
                d_res_sub      = TLake_r%d_lake(iunit)
                do j = 1, nlayers+1    
                    TLake_r%v_zt0(iunit,j) = TLake_r%v_zt(iunit,j)
                    TLake_r%a_d0(iunit,j) = TLake_r%a_d(iunit,j)
                    TLake_r%d_z0(iunit,j) = TLake_r%d_z(iunit,j)
                end do
                do j = 1,nlayers   
                    TLake_r%dv_nt(iunit,j)    = 0._r8
                    dv_in(j)    = 0._r8
                    dv_ou(j)    = 0._r8
                end do
                                           
        !if(iunit == 186781) then
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 0, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 0, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
            !write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 0, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18005,fmt="(i4, i4, i4, 3(e14.6))") 0, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%lake_outflow(iunit)*dtime,  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  0, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  0, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if

                !! lake precipitation and evaporation
                !d_prcp = THeat%forc_rain(iunit) * dtime  //TODO: not enable lake prec. now to avoid double counting
                v_prcp = d_prcp * TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)
                t_s        = TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)) 
                call QSat_str(t_s, THeat%forc_pbot(iunit),es)                
                evap =  lake_evap_r(iunit,dtime,TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)),F,es,sar)
                d_evap = evap * dtime                
                v_evap = d_evap * TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)
                if (v_evap - 0.95_r8 * TLake_r%d_v(iunit, TLake_r%d_ns(iunit)) > TINYVALUE) then ! too little lake storage for lake evap. to occur fully
                    v_evap = TLake_r%d_v(iunit, TLake_r%d_ns(iunit)) * 0.95_r8
                end if                
                if(v_evap < 0._r8) then
                    v_evap = 0._r8
                end if
                d_evap = v_evap / TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)
                
                ! advective inflow/outflow only, don't consider prec. or evap. fluxes here
                !call mosart_lake_hydro_sub_channel(iunit,nt,dtime)
                q_in = TLake_r%lake_inflow(iunit) ! + TLake_r%lake_rain(iunit) + TLake_r%lake_snow(iunit)
                !q_ou = -TLake_r%lake_outflow(iunit)
                if(q_in*dtime + TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1) - TLake_r%v_zti(iunit,ngeom+1) > TINYVALUE) then ! too much inflow so the lake will be more than full
                    dv_in(TLake_r%d_ns(iunit)) = q_in * dtime
                    TLake_r%lake_spillflow(iunit) = - TLake_r%lake_inflow(iunit)
                    v_spill = -TLake_r%lake_spillflow(iunit) * dtime
					dv_ou(TLake_r%d_ns(iunit)) = 0._r8
                    do j = 1,TLake_r%d_ns(iunit)-1   
                        dv_in(j)    = 0._r8
                        dv_ou(j)    = 0._r8
                    end do
                else
                    ! Calculation of flow contibution due to total inflow/outflow                 
                    call av_in_r(TLake_r%d_ns(iunit),q_in,t_in,dtime,iunit,dv_in)
                end if
        !if(iunit == 186781) then
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 1, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 1, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 1, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
            !write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 1, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
            !write(unit=18005,fmt="(i4, i4, i4, 3(e14.6))") 1, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%lake_outflow(iunit)*dtime,  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  1, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  1, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if

                !! Resize layer thickness and numbers based on inflow/outflow contribution                            
                TLake_r%dv_nt(iunit,:)    = 0._r8
                TLake_r%dv_nt(iunit,TLake_r%d_ns(iunit)) = dv_in(TLake_r%d_ns(iunit)) - v_spill + v_prcp - v_evap
                !TLake_r%dv_nt(iunit,TLake_r%d_ns(iunit)) = dv_in(TLake_r%d_ns(iunit)) + v_prcp - v_evap
                do j = 1, TLake_r%d_ns(iunit)-1
                    TLake_r%dv_nt(iunit,j)    = dv_in(j)                       
                end do               
                !     Calculate layer volume (m3)
                do j = 1, nlayers    
                    if (j<=TLake_r%d_ns(iunit))  then   
                        TLake_r%v_zn(iunit,j) = TLake_r%d_v(iunit,j) + TLake_r%dv_nt(iunit,j)
                    else
                        TLake_r%v_zn(iunit,j) = 0._r8
                    end if
                end do

        !if(iunit == 186781) then
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 3, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 2, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 2, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 4(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  2, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  2, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if
                ! Calculate layer thickness, depth, area/volume after net advection
                if(TLake_r%d_ns(iunit) > 1 .or. (TLake_r%d_ns(iunit)==1 .and. TLake_r%d_v(iunit,1) >= 0.1_r8*TLake_r%a_d(iunit,2))) then 
                !if(TLake_r%d_ns(iunit) > 1) then 
        !if(iunit == 186781) then
        !    write(unit=18016,fmt="(i4, i4,  6(e14.6))") 1, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%v_zn(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%dv_nt(iunit,1),TLake_r%v_zt(iunit,i+1), TLake_r%v_zti(iunit,ngeom+1)
        !end if
                    call lakegeom_update_r(iunit)
        !if(iunit == 186781) then
        !    write(unit=18016,fmt="(i4, i4,  6(e14.6))") 2, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%v_zn(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%dv_nt(iunit,1),TLake_r%v_zt(iunit,i+1), TLake_r%v_zti(iunit,ngeom+1)
        !end if
                else
        !if(iunit == 186781) then
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))") 3, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%v_zn(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1)
        !end if
                    TLake_r%d_v(iunit,1)    = TLake_r%v_zn(iunit,1)
                    TLake_r%v_zt(iunit,2)   = TLake_r%v_zn(iunit,1) + TLake_r%v_zt(iunit,1)
        !if(iunit == 186781) then
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))") 4, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%v_zn(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1)
        !end if
                end if                
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 2, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if

        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 3, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 3, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18005,fmt="(i4, i4, i4, 3(e14.6))") 3, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%lake_outflow(iunit)*dtime,  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 4(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 6(e14.6))")  0, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2)
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  3, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  3, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if

                call mosart_lake_hydro_main_channel(iunit,nt,dtime)
                q_ou = -TLake_r%lake_outflow(iunit)
                if(abs(q_ou) > TINYVALUE) then
                    call av_out_r(TLake_r%d_ns(iunit),q_ou,dtime,iunit,dv_ou,t_out)
                    tmp_tout = tmp_tout + t_out
                end if                
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 3, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
                
                !if(isnan(t_out)) then
                !    write(unit=1098,fmt="(i10, i4, (e14.6))") iunit, ww, q_ou
                !end if
        !if(iunit == 186781) then
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 2, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 4, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 4, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18005,fmt="(i4, i4, i4, 3(e14.6))") 4, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%lake_outflow(iunit)*dtime,  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  4, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  4, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if
                !! Resize layer thickness and numbers based on inflow/outflow contribution                            
                TLake_r%dv_nt(iunit,:)    = 0._r8
                !TLake_r%dv_nt(iunit,TLake_r%d_ns(iunit)) = dv_in(TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)) + v_prcp - v_evap
                TLake_r%dv_nt(iunit,TLake_r%d_ns(iunit)) = - dv_ou(TLake_r%d_ns(iunit))
                do j = 1, TLake_r%d_ns(iunit)-1
                    !TLake_r%dv_nt(iunit,j)    = dv_in(j)-dv_ou(j)                        
                    TLake_r%dv_nt(iunit,j)    = -dv_ou(j)                        
                end do               
                !     Calculate layer volume (m3)
                do j = 1, nlayers    
                    if (j<=TLake_r%d_ns(iunit))  then   
                        TLake_r%v_zn(iunit,j) = TLake_r%d_v(iunit,j) + TLake_r%dv_nt(iunit,j)
                    else
                        TLake_r%v_zn(iunit,j) = 0._r8
                    end if
                end do
                !do j = 1, TLake_r%d_ns(iunit)
                !    TLake_r%d_v(iunit,j)    = TLake_r%v_zn(iunit,j)                       
                !end do               
                
        !if(iunit == 186781) then
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 3, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 5, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 5, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 4(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 6(e14.6))")  0, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2)
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  5, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  5, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if
                ! Calculate layer thickness, depth, area/volume after net advection                                 
                if(TLake_r%d_ns(iunit) > 1 .or. (TLake_r%d_ns(iunit)==1 .and. TLake_r%d_v(iunit,1) >= 0.1_r8*TLake_r%a_d(iunit,2))) then 
                !if(TLake_r%d_ns(iunit) > 1) then 
                    call lakegeom_update_r(iunit)
                else
                    TLake_r%d_v(iunit,1)    = TLake_r%v_zn(iunit,1)
                    TLake_r%v_zt(iunit,2)   = TLake_r%v_zn(iunit,1) + TLake_r%v_zt(iunit,1)
                end if                
                
        !if(iunit == 186781) then 
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 4, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 4, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 6, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 6, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !   write(unit=18005,fmt="(i4, i4, i4, 3(e14.6))") 6, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%lake_outflow(iunit)*dtime,  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  6, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  6, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if
                if ((TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1) - TLake_r%v_zti(iunit,ngeom+1))/TLake_r%v_zti(iunit,ngeom+1) > 0.05_r8) then ! storage is too large
                    call over_flow_r(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)                        
                end if
                if ((TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1) - TLake_r%v_zti(iunit,1))/TLake_r%v_zti(iunit,ngeom+1) < -TINYVALUE) then ! storage dries up
                    call under_flow_r(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)                        
                end if

                ! Adjust layers thickness for numerical stability                    
                !    Check if any layers are too thin
                if (TLake_r%d_ns(iunit) == 2) then
                    if (TLake_r%dd_z(iunit,2) <= ddz_top_min .or. TLake_r%dd_z(iunit,1) <= ddz_min) then
        !if(iunit == 186781) then
        !    write(unit=1009,fmt="(i4, i4, 3(e14.6))") 11, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
                        call toplayer_merge_r(iunit)
                    end if                                
        !if(iunit == 186781) then
        !    write(unit=1009,fmt="(i4, i4, 3(e14.6))") 12, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
                elseif(TLake_r%d_ns(iunit) > 2) then
                    if (TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)) <= ddz_top_min) then
        !if(iunit == 186781) then
        !    write(unit=1009,fmt="(i4, i4, 3(e14.6))") 21, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
                        call toplayer_merge_r(iunit)
                    end if   
        !if(iunit == 186781) then
        !    write(unit=1009,fmt="(i4, i4, 3(e14.6))") 22, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
                    
                    do i=1,TLake_r%d_ns(iunit)-1
                        if (TLake_r%dd_z(iunit,i) <= ddz_min) then
                            call layer_merge_r(iunit,i)
                            exit ! only allow merging occur once per each time step
                        end if                                
                    end do                        
                    ! repeat one more time in case any layer is still too thin after merging once (and layer number changed)
                    do i=1,TLake_r%d_ns(iunit)-1
                        if (TLake_r%dd_z(iunit,i) <= ddz_min) then
                            call layer_merge_r(iunit,i)
                            exit ! only allow merging occur once per each time step
                        end if                                
                    end do                        
        !if(iunit == 186781) then
        !    write(unit=1009,fmt="(i4, i4, 3(e14.6))") 23, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if

                end if                        
        !if(iunit == 186781) then
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 5, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 5, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 7, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 7, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 3(e14.6))") 7, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%lake_outflow(iunit)*dtime,  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  7, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  7, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if

                !    Check if any layers are too thick
                if (TLake_r%d_ns(iunit) > 1) then
                    if(TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)) > ddz_top_max) then
                       call layer_split_r(iunit,TLake_r%d_ns(iunit))
                    end if 
                    do i=1,TLake_r%d_ns(iunit)-1
                        if(TLake_r%dd_z(iunit,i) >= ddz_max) then 
                            call layer_split_r(iunit,i)
                            exit
                        end if                                
                    end do
                    ! repeat one more time in case any layer is still too thick after splitting once (and layer number changed)
                    do i=1,TLake_r%d_ns(iunit)-1
                        if(TLake_r%dd_z(iunit,i) >= ddz_max) then 
                            call layer_split_r(iunit,i)
                            exit
                        end if                                
                    end do                    
        !if(iunit == 186781) then
        !    write(unit=1009,fmt="(i4, i4, 3(e14.6))") 24, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
                end if

    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 6, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
        !if(iunit == 186781) then
        !    write(unit=1008,fmt="(i4, i4, 3(e14.6))") 6, TLake_r%d_ns(iunit),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,1), TLake_r%dd_z(iunit,2)
        !end if
        !if(iunit == 186781) then 
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 8, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18004,fmt="(i4, i4, i4, 9(e14.6))") 8, ww, TLake_r%d_ns(iunit),  sum(TLake_r%v_zn(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(TLake_r%d_v(iunit,:)) + sum(dv_in) - sum(dv_ou) + v_prcp - v_evap, sum(dv_in), sum(dv_ou), v_prcp - v_evap, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), dv_ou(TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)) - dv_ou(TLake_r%d_ns(iunit)),  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !    write(unit=18005,fmt="(i4, i4, i4, 8(e14.6))")  8, ww, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,1), TLake_r%dv_nt(iunit,1), TLake_r%d_v(iunit,2), TLake_r%dv_nt(iunit,2), TLake_r%d_v(iunit,3), TLake_r%dv_nt(iunit,3)
        !    write(unit=18006,fmt="(i4, i4,  3(e14.6))")  8, ww, TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1), TLake_r%v_zti(iunit,ngeom+1), TLake_r%v_zn(iunit,TLake_r%d_ns(iunit))
        !end if

                tmp_outflow = tmp_outflow + TLake_r%lake_outflow(iunit) + TLake_r%lake_spillflow(iunit)
                tmp_evap = tmp_evap + v_evap/dtime ! unit m3/s here, note the volume of evap. is calculated only over effective lake areas
                TLake_r%dV_str(iunit) = TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit) + TLake_r%lake_spillflow(iunit) + (v_prcp-v_evap)/dtime
                TLake_r%V_str(iunit) = TLake_r%V_str(iunit) + TLake_r%dV_str(iunit) * dtime ! dV_str = inflow - outflow
        !if(iunit == 186781) then 
        !    write(unit=18001,fmt="(i4, i4, i4, 7(e14.6))") 9, ww, TLake_r%d_ns(iunit),  sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)- TLake_r%v_zt(iunit,1), TLake_r%V_str(iunit) - TLake_r%v_zt(iunit,1), TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, (v_prcp-v_evap), (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit))*dtime + (v_prcp-v_evap)
        !    write(unit=18005,fmt="(i4, i4, i4, 3(e14.6))") 9, ww, TLake_r%d_ns(iunit),  TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%lake_outflow(iunit)*dtime,  TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  0, ww, 0, sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  0, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !end if
                !!! Please DO NOT delete. necessary for diagnosis & debugging 
                !Call lake_r_waterbalance_check(iunit)


                !!======================== energy balance processes   
                ! initialize for current sub time step
                phi_o          = 0._r8 !
                sh_net         = 0._r8            
                do j = 1,nlayers   
                    enr_0(j)    = 0._r8
                    enr_in(j)    = 0._r8
                    enr_ou(j)    = 0._r8
                    enr_1(j)    = 0._r8
                    enr_2(j)    = 0._r8
                    enr_phi(j)    = 0._r8
                    phi_z(j)    = 0._r8
                    rho_z(j)    = 0._r8
                    a(j)        = 0._r8
                    b(j)        = 0._r8
                    c(j)        = 0._r8
                    r(j)        = 0._r8
                    Fx(j)        = 0._r8
                end do
                           
                !     Intitialize/reassign layer storage
                !    Assign layer temperature and calculate reservoir density at depth z and incoming flow     
                !     Allocate layer temperature as old for assigning counter during averaging sub-timestep result        
                rho_r = den(t_in)                
                do j = 1, nlayers    
                    if (j<=TLake_r%d_ns(iunit)) then
                        !TLake_r%temp_lake(iunit,j) = TLake_r%temp_lake(iunit,j)
                        rho_z(j) = den(TLake_r%temp_lake(iunit,j))
                        t_z_old(j) = TLake_r%temp_lake(iunit,j)
                    else 
                        TLake_r%temp_lake(iunit,j) = 0._r8
                        rho_z(j) = 0._r8
                        t_z_old(j) = 0._r8
                    end if
                end do    
                                                               
                ! Calculate total mass (kg), layer energy (w)
                do j = 1, nlayers    
                    if (j<=TLake_r%d_ns(iunit))  then
                        TLake_r%v_zo(iunit,j) = TLake_r%d_v(iunit,j)
                        enr_0(j) = TLake_r%temp_lake(iunit,j)*TLake_r%d_v(iunit,j)*rho_z(j)*c_w/dtime
                    else
                        TLake_r%v_zo(iunit,j) = 0._r8
                        enr_0(j) = 0._r8
                    end if
                end do


                ! Lake is too empty/shallow, no need to invoke the energy balance method
                if(TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1) <= TINYVALUE) then    
                    do j = 1, nlayers    
                        if (j<=TLake_r%d_ns(iunit)) then
                            TLake_r%temp_lake(iunit,j) = THeat%Tr(iunit)
                        else 
                            TLake_r%temp_lake(iunit,j) = 0._r8
                        end if
                    end do
                    TLake_r%lake_Tsur(iunit) = THeat%Tr(iunit)
                    TLake_r%lake_Tout(iunit) = THeat%Tr(iunit)
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 7, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
                    cycle
                end if 
                                
                ! Calculate layer internal energy (w) due to inflow/outflow        
                do j = 1, nlayers    
                    if (j<=TLake_r%d_ns(iunit))  then   
                        enr_in(j) = (dv_in(j)*t_in*rho_r*c_w/dtime)        !Energy from inflow
                        enr_ou(j) = (dv_ou(j)*TLake_r%temp_lake(iunit,j)*rho_z(j)*c_w/dtime)    !Energy loss due to outflow
                    else
                        enr_in(j) = 0._r8
                        enr_ou(j) = 0._r8
                    end if
                end do
                
                ! Calculation of Surface fluxes and heat source 
                !    Net shortwave radiation (w/m^2)
                if (THeat%coszen(iunit)>TINYVALUE) then 
                    alb_s = 0.05_r8/(THeat%coszen(iunit) + 0.15_r8)!THeat%albedo(iunit)!
                else
                    alb_s = 0.06_r8    !0.06_r8
                end if
                ! alb_s = 0.03_r8    !0.06_r8    
                sh_net     = max(THeat%forc_solar(iunit)*(1._r8 - alb_s),0._r8)!            
                lw_abr     = (1._r8 - 0.03_r8)*THeat%forc_lwrad(iunit) !
                lw_ems     = 0.97_r8*st_bl*t_s**4_r8    
                sn_heat    = 1.5701_r8*THeat%forc_wind(iunit)* (t_s - THeat%forc_t(iunit))
                le         = 1000._r8*(2499.64_r8 - 2.51_r8 * (t_s-273.15_r8)) 

            !if(iunit == 186781) then 
            !    write(unit=88001,fmt="(i4, i4, i4, i4, 7(e18.10))") 1, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !end if

!        if(iunit == 186781) then
!            write(unit=12002,fmt="(i4, i4, i4, 6(e18.10))") 2, ww, dflag,  TLake_r%V_str(iunit), TLake_r%dV_str(iunit) * dtime, TLake_r%lake_inflow(iunit)*dtime, TLake_r%lake_outflow(iunit)*dtime, -v_evap, (TLake_r%lake_rain(iunit) + TLake_r%lake_snow(iunit))*dtime
!        end if
        !if(iunit == 186781) then
        !    write(unit=18001,fmt="(i4, i4, i4, 10(e14.6))") 8, ww, TLake_r%d_ns(iunit), sum(dv_in), sum(dv_ou), v_evap, TLake_r%lake_outflow(iunit)*dtime, sum(TLake_r%dv_nt(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit), TLake_r%v_zti(iunit,1), TLake_r%v_zti(iunit,ngeom+1)
        !    write(unit=18002,fmt="(i4, i4, i4, 6(e18.10))")  8, ww, TLake_r%d_ns(iunit), sum(TLake_r%d_v(iunit,:)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), sum(TLake_r%d_v(iunit,:))-TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)+TLake_r%v_zt(iunit,1), sum(TLake_r%d_v(iunit,:))-TLake_r%V_str(iunit)+TLake_r%v_zt(iunit,1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)-TLake_r%V_str(iunit) 
        !    write(unit=18003,fmt="(i4, i4, i4, 2(e14.6))")  8, ww, TLake_r%d_ns(iunit), sum(TLake_r%dv_nt(iunit,:)), q_in*dtime-q_ou*dtime - v_evap
        !end if

                ! latent heat needs to be consistent with evaporation    
                if (t_s > t_frz) then
                    lt_heat = rho_w * evap* le    !latent heat (w/m^2)
                else 
                    lt_heat = 0.0_r8
                end if

                !     Net surface heat flux  
                phi_o = sh_net + lw_abr - (lt_heat + lw_ems + sn_heat)!

!        if(iunit == 186781) then
!            write(unit=61001,fmt="(i4, i4, i4, 9(e18.10))") 0, ww, dflag, phi_o, sh_net, lw_abr, lt_heat, lw_ems, sn_heat, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), THeat%forc_t(iunit), t_in
!            write(unit=32000,fmt="(10(e14.6))") evap, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), THeat%forc_t(iunit), phi_o, sh_net, phi_o - sh_net, lw_abr, lt_heat, lw_ems, sn_heat
!        end if

    ! *********************************************************************************************************************************
            !     Calculate layer energy (w)
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 8, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
                
                if(TUnit_lake_r%one_layer(iunit) == 1 .or. TLake_r%d_lake(iunit) < 20._r8) then !if high-latitude or shallow lakes, empirical approach
                    if(THeat%forc_vp(iunit) < TINYVALUE) then
                        do j = 1, nlayers    
                            if (j<=TLake_r%d_ns(iunit))  then   
                                TLake_r%temp_lake(iunit,j) = THeat%Tr(iunit)
                            else
                                TLake_r%temp_lake(iunit,j) = 0._r8
                            end if
                        end do
                        TLake_r%lake_Tsur(iunit) = THeat%Tr(iunit)
                        TLake_r%lake_Tout(iunit) = THeat%Tr(iunit)
                        !write(iulog,*) 'Warning vapor pressure less than zero ! ', iunit, THeat%forc_vp(iunit)
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 9, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
    !do i=1, TLake_r%d_ns(iunit)
    !if(.not.(isnan(t_in)) .and. isnan(TLake_r%temp_lake(iunit,i))) then
    !    write(unit=1899,fmt="(i10, i4, i4, i4, 2(e14.6))") iunit, 1, i, TLake_r%d_ns(iunit), t_in, TLake_r%temp_lake(iunit,i)
    !    call shr_sys_abort('Goofy temperature -- NAN values #13')
        
    !end if
    !end do


    
                        cycle
                    else                    
                        temp_equi = TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)) - 273.15
                        mytmp(1) = QSat_huang(THeat%forc_t(iunit)-273.15_r8) * 0.001_r8 ! Pa to kPa
                        mytmp(2) = THeat%forc_vp(iunit) * 0.001_r8 ! Pa to kPa
                        mytmp(3) = THeat%forc_t(iunit) - 273.15 !kelvin to celcus
                        mytmp(4) = THeat%forc_wind(iunit) !
                        mytmp(5) = THeat%forc_pbot(iunit) * 0.001_r8 ! Pa to kPa
                        mytmp(6) = sqrt(TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)) ! here simplify Lf calculation without accounting for the wind direction and lake shape, 
                                                                                  ! lake temperature calculation seems not very sensentive to Lf
                        mytmp(7) = TLake_r%d_lake(iunit) !
                        mytmp(8) = rtmCTL%latc(iunit) * SHR_CONST_PI/180._r8 ! degree to radians
                        mytmp(9) = THeat%forc_solar(iunit) * 1e-6 * 86400 ! J·m-2·s-1 (w/m^2) to MJ·m-2·d-1
                        mytmp(10)= TUnit_lake_r%elev(iunit) !
                        mytmp(11)= alb_s !
                        mytmp(12)= TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)) - 273.15 !11.63
                        mytmp(13)= doy()*1._r8 ! return day of the year corresponding to the current time step
                        mytmp(14)= dtime/86400 ! seconds to day
                        temp_equi = temp_shallow_lake(mytmp) + 273.15 ! celcus to Kelvin
                    
                !if(isnan(temp_equi)) then
                !    write(unit=1099,fmt="(i10, 13(e14.6))") iunit, mytmp(1),mytmp(2),mytmp(3),mytmp(4),mytmp(5),mytmp(6),mytmp(7),mytmp(8),mytmp(9),mytmp(10),mytmp(11),mytmp(12),mytmp(13)
                    !call shr_sys_abort('Goofy temperature -- NAN values #1')
                    
                !end if
                        do j = 1, nlayers
                            if (j<=TLake_r%d_ns(iunit))  then   
                                TLake_r%temp_lake(iunit,j) = temp_equi
                            else
                                TLake_r%temp_lake(iunit,j) = 0._r8
                            end if
                        end do
                    end if
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 10, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if

            !if(iunit == 186781) then
            !    write(unit=71001,fmt="(i4, i4, i4, i4, 7(e18.10))") 0, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !    write(unit=71000,fmt="(i4, i4, i4, i4, 7(e18.10))") 0, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=11001,fmt="(14(e14.6))") mytmp(1), mytmp(2), mytmp(3), mytmp(4), mytmp(5), mytmp(6), mytmp(7), mytmp(8), mytmp(9), mytmp(10), mytmp(11), mytmp(12), mytmp(13), mytmp(14)
            !end if

                else
                                    
            !if(iunit == 186781) then
            !    write(unit=11001,fmt="(i4, 8(e18.10))") TLake_r%d_ns(iunit), enr_1(4), enr_0(4), enr_1(3), enr_0(3), enr_1(2), enr_0(2), enr_1(1), enr_0(1)
            !    write(unit=11002,fmt="(i4, 7(e18.10))") rho_r*c_w/dtime, t_in, dv_in(4), TLake_r%temp_lake(iunit,4), dv_ou(4), t_in, dv_in(3), TLake_r%temp_lake(iunit,3), dv_ou(3)
            !    write(unit=81002,fmt="(i4, i4, i4, i4, 7(e18.10))") 2, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !    write(unit=81000,fmt="(i4, i4, i4, i4, 7(e18.10))") 2, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !end if
                !    Check energy balance (w) after advective mixing    
                    !enr_err1 = (sum(enr_1) - (sum(enr_0)+ sum(enr_in) - sum(enr_ou)))
                    
                !******************************************************************************
                !     Calculate solar energy absorbed at each layer        
                    eta  = 1.4256_r8*(TLake_r%d_lake(iunit))**(-0.424_r8) ! as used in Subin et al, 2011 (citing Hakanson, 1995) but modified for actual reservoir depth, corrected following disucssion with Xing Fang 
                    beta = 0.175_r8 !
            
                    do j=1,nlayers+1
                        phi_x(j) = 0._r8
                    end do
                    
                    if (sh_net > TINYVALUE .and. TLake_r%d_ns(iunit) >1) then
                        k=0
                        top_d=TLake_r%d_lake(iunit)-TLake_r%d_z(iunit,TLake_r%d_ns(iunit)-k)
                        if(top_d < 0.61_r8) then
                            k=k+1
                        else
                            k=k+2
                        end if    
                        v_mix=(TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1)-TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)-k))    
                        !    Solar radiation energy absorbed at the mixed zone
                        sh_mix=sh_net*sar*(TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)-(1._r8-beta)*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)-k))
                        j=TLake_r%d_ns(iunit)-k
                        do i=j,TLake_r%d_ns(iunit)
                           phi_z(i)=sh_mix*TLake_r%d_v(iunit,i)/v_mix
                        end do
                        phi_x(TLake_r%d_ns(iunit)+1)=sh_net
                        phi_x(TLake_r%d_ns(iunit)-k)=(1._r8-beta)*sh_net
                        if(k>0) then
                            do i=1,k
                               ii=TLake_r%d_ns(iunit)-i+1
                               phi_x(ii)=(TLake_r%a_d(iunit,ii+1)*phi_x(ii+1)-phi_z(ii))/(TLake_r%a_d(iunit,ii))
                            end do
                        end if
                        
                        !    Solar radiation energy absorbed at sub layers
                        l=TLake_r%d_ns(iunit)-k-1
                        do j=1,l
                           i=(TLake_r%d_ns(iunit)-k)-j+1
                           phi_x(i-1)=phi_x(i)*exp(-eta*(TLake_r%d_z(iunit,i)-TLake_r%d_z(iunit,i-1)))
                        end do
                    
                        !    Solar radiation energy absorbed in each layer
                        j=TLake_r%d_ns(iunit)-k-1
                        do i=1,j
                           phi_z(i)=(sar*TLake_r%a_d(iunit,i+1)*phi_x(i+1)-sar*TLake_r%a_d(iunit,i)*phi_x(i))
                        end do
                    elseif (sh_net > TINYVALUE .and. TLake_r%d_ns(iunit)==1) then
                        phi_z(TLake_r%d_ns(iunit))=sh_net*sar*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)
                    else
                        do j=1,TLake_r%d_ns(iunit)
                            phi_z(j) = 0._r8
                        end do    
                    end if            
                    
                ! *********************************************************************************************************************************
                !     Calculation of effective diffusion coefficient     , Herb and Stefan (2004) citing Wu (1971)
                    if (THeat%forc_wind(iunit) >= 15._r8) then
                        c_d = 2.6e-3
                    else
                        c_d = 5.e-4*sqrt(THeat%forc_wind(iunit))
                    end if
                    tau = rho_a*c_d*THeat%forc_wind(iunit)**2._r8 ! Shear stress at surface
                    s_vel = sqrt(tau/rho_w) ! Shear velocity at surface
                    k_ew=tau*s_vel*sar*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)*dtime    ! Wind driven kinetic energy at surface
                    dis_w = k_ew/(rho_w*TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1)*dtime)  ! rate of dissipation-wind
                    cfw = 1.e-02!
                    cfa = 1.e-05!
                    k_m = 0.57_r8/(c_w*rho_w)                           !molecular diffusivity
                    df_eff(1)= 0._r8            !bottom interface
                    df_eff(TLake_r%d_ns(iunit)+1)= 0._r8        !air interface
                    do j = 2,TLake_r%d_ns(iunit) 
                        q_adv(j) = max((dv_in(j)+dv_ou(j)),0._r8)
                        k_ad(j)=0.5_r8*rho_w*q_adv(j)*dtime*(q_adv(j)/(TUnit_lake_r%Width(iunit)*TLake_r%dd_z(iunit,j)))**2._r8 ! Advection driven kinetic energy 
                        dis_ad(j)= k_ad(j)/(rho_w*TLake_r%v_zt(iunit,j)*dtime)    ! rate of dissipation-inflow/outflow                    
                    ! Calculate Richardson number
                        drhodz(j) = (rho_z(j-1)-rho_z(j))/0.5_r8*(TLake_r%dd_z(iunit,j)+TLake_r%dd_z(iunit,j-1))
                        bv_f = max((grav/rho_w)*drhodz(j),0._r8)
                        if (s_vel <= TINYVALUE) ri = 0._r8
                        ri = bv_f/((s_vel/(0.4_r8*TLake_r%d_z(iunit,j)))**2._r8)                
                    ! Calculate Froude number
                        l_vel = q_adv(j)*TUnit_lake_r%Length(iunit)/(sar*TLake_r%a_d(iunit,j)*TLake_r%dd_z(iunit,j))
                        if (q_adv(j) <= -TINYVALUE .or. drhodz(j) <= -TINYVALUE) then
                            Fr(j) = 0._r8
                        else    
                            Fr(j)= (grav*TLake_r%dd_z(iunit,j)*drhodz(j)/rho_w)/l_vel**2._r8
                        end if                        
                    ! Calculate diffusion coefficients                
                        df_eff(j)=min(max(dtime**2._r8*((cfw*dis_w/(1+ri))+(0.5_r8*cfa*(dis_ad(j)+dis_ad(j-1))/(1+Fr(j)))),k_m),5.56e-03) !5.56e-03!
                    end do
                
                !*****************************************************************
                ! Calculate matrix elements
                    do j = 1,TLake_r%d_ns(iunit)
                        if (j == 1 .and. TLake_r%d_ns(iunit)>1) then
                            m1(j) = 2._r8*dtime/(0.5_r8*(sar*TLake_r%a_d(iunit,j)+sar*TLake_r%a_d(iunit,j+1))*TLake_r%dd_z(iunit,j))
                            m2(j) = m1(j)*sar*TLake_r%a_d(iunit,j+1)*df_eff(j+1)/(TLake_r%dd_z(iunit,j)+TLake_r%dd_z(iunit,j+1))
                            m3(j) = 0._r8                        
                            Fx(j) = dtime*phi_z(j)/(TLake_r%d_v(iunit,j)*c_w*rho_z(j))
                            a(j) = - (m2(j))
                            b(j) = 1._r8 + (m2(j) + m3(j)) 
                            c(j) = 0._r8 
                            r(j) = TLake_r%temp_lake(iunit,j) + Fx(j) ! bottom boundary condition 
                        elseif (j <= TLake_r%d_ns(iunit)-1 .and. TLake_r%d_ns(iunit)>2) then
                            m1(j) = 2._r8*dtime/(0.5_r8*(sar*TLake_r%a_d(iunit,j)+sar*TLake_r%a_d(iunit,j+1))*TLake_r%dd_z(iunit,j))
                            m2(j) = m1(j)*sar*TLake_r%a_d(iunit,j+1)*df_eff(j+1)/(TLake_r%dd_z(iunit,j)+TLake_r%dd_z(iunit,j+1))
                            m3(j) = m1(j)*sar*TLake_r%a_d(iunit,j)*df_eff(j)/(TLake_r%dd_z(iunit,j)+TLake_r%dd_z(iunit,j-1))                        
                            Fx(j) = dtime*phi_z(j)/(TLake_r%d_v(iunit,j)*c_w*rho_z(j))
                            a(j) = - m2(j)
                            b(j) = 1._r8 + m2(j) + m3(j) 
                            c(j) = - m3(j)
                            r(j) = TLake_r%temp_lake(iunit,j) + Fx(j)
                        elseif (j == TLake_r%d_ns(iunit)) then!top layer
                            m1(j) = 2._r8*dtime/(0.5_r8*(sar*TLake_r%a_d(iunit,j)+sar*TLake_r%a_d(iunit,j+1))*TLake_r%dd_z(iunit,j))
                            m2(j) = 0._r8
                            m3(j) = m1(j)*sar*TLake_r%a_d(iunit,j)*df_eff(j)/(TLake_r%dd_z(iunit,j)+TLake_r%dd_z(iunit,j-1))
                            Fx(j) = dtime*((phi_o-sh_net)*sar*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)+phi_z(j))/(TLake_r%d_v(iunit,TLake_r%d_ns(iunit))*c_w*rho_z(j)) ! 
                            a(j) = 0._r8
                            b(j) = 1._r8 + (m2(j) + m3(j)) 
                            c(j) = - (m3(j))
                            r(j) = TLake_r%temp_lake(iunit,j) + Fx(j) ! top boundary condition                         
                        end if
                    end do    
            !if(iunit == 186781) then
            !    write(unit=81003,fmt="(i4, i4, i4, i4, 7(e18.10))") 3, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !    write(unit=81000,fmt="(i4, i4, i4, i4, 7(e18.10))") 3, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !end if
!          if(iunit == 186781) then 
!            write(unit=61002,fmt="(10(e18.10))") Fx(TLake_r%d_ns(iunit)), dtime, (phi_o-sh_net)*sar*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)+phi_z(TLake_r%d_ns(iunit)), c_w, rho_z(TLake_r%d_ns(iunit)), (phi_o-sh_net)*sar*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1), phi_z(TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), THeat%forc_t(iunit), t_in
!            write(unit=32000,fmt="(10(e14.6))") evap, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), THeat%forc_t(iunit), phi_o, sh_net, phi_o - sh_net, lw_abr, lt_heat, lw_ems, sn_heat
!        end if
                  
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 11, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
                    
                !    Solve for temperature
                    call solve_r(a,b,c,r,iunit,TLake_r%d_ns(iunit))
            !if(iunit == 186781) then
            !    write(unit=81004,fmt="(i4, i4, i4, i4, 7(e18.10))") 4, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !    write(unit=81000,fmt="(i4, i4, i4, i4, 7(e18.10))") 4, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !end if
!            if(iunit == 186781) then
!                write(unit=98000,fmt="(i4, i4, i4, i4, 7(e18.10))") 3, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=88003,fmt="(i4, i4, i4, i4, 7(e18.10))") 3, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=78000,fmt="(i4, 4(e18.10))") 3, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2),THeat%forc_t(iunit)
!            end if
    !if(iunit == 186781) then 
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 12, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !end if
                    
            !if(iunit == 186781) then
            !    write(unit=12002,fmt="(i4, i4, i4, (e14.6))") 4, ww, dflag, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit))
            !end if
                
                    ! Calculate layer intermediate internal energy (w)
                    if(0>1) then
                    do j = 1, nlayers    
                        if (j<TLake_r%d_ns(iunit))  then   
                            enr_2(j)   = TLake_r%temp_lake(iunit,j)*TLake_r%d_v(iunit,j)*rho_z(j)*c_w/dtime
                            enr_phi(j) = phi_z(j)
                        elseif (j==TLake_r%d_ns(iunit))  then
                            enr_2(j)   = TLake_r%temp_lake(iunit,j)*TLake_r%d_v(iunit,j)*rho_z(j)*c_w/dtime
                            enr_phi(j) = (phi_o-sh_net)*sar*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)        !Enegry from net surface flux (w)
                        else
                            enr_2(j) = 0.
                        end if
                    end do
                    end if                    
                
            !    Check energy balance error (w) after triadiagonal matrix solution    
                    !enr_err2 = (sum(enr_2) - (sum(enr_1) + sum(enr_phi)))
                    
            !***********************************************************************************************************************        
            ! !    Solve convective mixing
                    temp_prev = TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit))
                ! Recalculate layer density 
                    do j = 1,TLake_r%d_ns(iunit)   
                        rho_z(j) = den(TLake_r%temp_lake(iunit,j))
                    end do
                    
            !if(iunit == 186781) then
            !    write(unit=21001,fmt="(i4, i4, i4, i4, 6(e18.10))") 1, ww, dflag, TLake_r%d_ns(iunit), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-3),t_in, THeat%forc_t(iunit)
            !end if
                ! Check if instability exists
                    do j = 1,TLake_r%d_ns(iunit)-1 
                        if(rho_z(j) < rho_z(j+1)) then
                            ! Start mixing layers
                            mix1=j
                            mix2=mix1
                            sumvol=TLake_r%d_v(iunit,mix2)*1._r8
                            tsum=TLake_r%temp_lake(iunit,mix2)*TLake_r%d_v(iunit,mix2)    
                            
                            mix2=mix2+1
                            mixvol2=TLake_r%d_v(iunit,mix2)*1._r8
                            sumvol=sumvol + mixvol2
                            tsum=tsum+TLake_r%temp_lake(iunit,mix2)*mixvol2
                            tmix=tsum/sumvol
                            
                            ! Calculate density of mixed layer    
                            denmix = den(tmix)
                            
                            ! Check if instability exists below mixed layer    and mix layers    
                            if(rho_z(mix1-1) < denmix .and. mix1 >= 2) then
                                mix1=mix1-1    
                            
                                ! Calculate temperature of mixed layer    
                                mixvol1=TLake_r%d_v(iunit,mix1)*1._r8
                                sumvol=sumvol + mixvol1
                                tsum=tsum+TLake_r%temp_lake(iunit,mix1)*mixvol1
                                tmix= tsum/sumvol
                            end if
                            
                            ! Calculate density of mixed layer    
                            denmix = den(tmix)                    
                      
                            ! Set new layer temperature and density 
                            do i=mix1,mix2
                                rho_z(i)=denmix
                                TLake_r%temp_lake(iunit,i)=tmix
                            end do
            !if(iunit == 508 .and. j == (TLake_r%d_ns(iunit)-1)) then
            !    write(unit=81000,fmt="(i4, i4, i4, 3(e14.6))") ww, mix1, mix2, tmix, TLake_r%temp_lake(iunit,mix1), TLake_r%temp_lake(iunit,mix2) 
            !end if
    !if(iunit == 186781) then
    !    write(unit=78000,fmt="(i4, i4, 4(e18.10))") 10, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2),THeat%forc_t(iunit)
    !    if(isnan(TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)))) then
    !        call shr_sys_abort('Goofy temperature -- NAN values #8')
    !    end if
    !end if
!            if(iunit == 186781) then
!                write(unit=11001,fmt="(i4, i4, i4, 5(e14.6))") 2, ww, dflag, TLake_r%temp_lake(iunit,mix2), tmix, tsum, sumvol, denmix 
!            end if

!            if(iunit == 186781) then
!                write(unit=98000,fmt="(i4, i4, i4, i4, 7(e18.10))") 4, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=88004,fmt="(i4, i4, i4, i4, 7(e18.10))") 4, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!            end if
                            ! TODO: toplayer temp. does not go above freezing point if it is previously frozen and air temp. is also below
                            if(TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)) > 273.15_r8 .and. THeat%forc_t(iunit) <= 273.15) then 
                            if(temp_prev < 273.15_r8 - TINYVALUE) then 
                                TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)) = 273.15_r8
                            end if
                            end if

                        end if
                    end do
    !if(iunit == 186781) then
    !    write(unit=78000,fmt="(i4, i4, 5(e18.10))") 13, ww, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), t_in, THeat%forc_t(iunit)
    !    if(isnan(TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)))) then
    !        call shr_sys_abort('Goofy temperature -- NAN values #8')
    !    end if
    !end if

            !if(iunit == 186781) then
            !    write(unit=21001,fmt="(i4, i4, i4, i4, 6(e18.10))") 2, ww, dflag, TLake_r%d_ns(iunit), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-3),t_in, THeat%forc_t(iunit)
            !end if
            !if(iunit == 186781) then
            !    write(unit=81005,fmt="(i4, i4, i4, i4, 7(e18.10))") 5, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !    write(unit=81000,fmt="(i4, i4, i4, i4, 7(e18.10))") 5, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
            !end if
!            if(iunit == 186781) then
!                write(unit=98000,fmt="(i4, i4, i4, i4, 7(e18.10))") 5, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=88005,fmt="(i4, i4, i4, i4, 7(e18.10))") 5, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=78000,fmt="(i4, 4(e18.10))") 5, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2),THeat%forc_t(iunit)
!            end if
                
                end if
!            if(iunit == 186781) then
!                write(unit=98000,fmt="(i4, i4, i4, i4, 7(e18.10))") 6, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=88006,fmt="(i4, i4, i4, i4, 7(e18.10))") 6, ww, dflag, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1),t_in, THeat%forc_t(iunit)
!                write(unit=78000,fmt="(i4, 4(e18.10))") 6, TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-2),THeat%forc_t(iunit)
!            end if
                
                
                !! for the purpose of lake water balance check
                !if(iunit == 186781) then
                !    write(unit=12210,fmt="(i10, 6(e18.10))") iunit,TLake_r%V_str(iunit)/dtime, TLake_r%dV_str(iunit), TLake_r%lake_inflow(iunit), TLake_r%lake_outflow(iunit), -v_evap/dtime, TLake_r%lake_rain(iunit)+TLake_r%lake_snow(iunit)
                !    write(unit=12211,fmt="(i10, 6(e18.10))") iunit,TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1)/dtime, q_in, q_ou, TLake_r%lake_inflow(iunit)+TLake_r%lake_rain(iunit)+TLake_r%lake_snow(iunit), TLake_r%lake_outflow(iunit), -v_evap/dtime 
                !    write(unit=12212,fmt="(i10, 6(e18.10))") iunit,TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%V_str(iunit), q_in, (TLake_r%lake_inflow(iunit)+TLake_r%lake_rain(iunit)+TLake_r%lake_snow(iunit)), q_ou, -TLake_r%lake_outflow(iunit) 
                !end if
        !    ! gulu
        !    v_end1 = TLake_r%V_str(iunit)
        !    v_end2 = TLake_r%v_zt(iunit,TLake_r%d_ns(iunit))
        !    delta_v = (TLake_r%lake_inflow(iunit) + TLake_r%lake_outflow(iunit) + TLake_r%lake_evap(iunit) + TLake_r%lake_rain(iunit) + TLake_r%lake_snow(iunit))*localDeltaT
        !    if(iunit==499) then
        !        write(unit=1119,fmt="(i10, i4, 7(e14.6))") iunit, ww, v_end1 - v_strt1-delta_v, (v_end1 - v_strt1)/localDeltaT, delta_v/localDeltaT, TLake_r%lake_inflow(iunit), TLake_r%lake_outflow(iunit), TLake_r%lake_evap(iunit), TLake_r%lake_rain(iunit) + TLake_r%lake_snow(iunit)
        !    end if
        !    if(iunit==499) then
        !        write(unit=2119,fmt="(i10, i4, 7(e14.6))") iunit, ww, v_end2 - v_strt2-delta_v, (v_end2 - v_strt2)/localDeltaT, delta_v/localDeltaT, TLake_r%lake_inflow(iunit), TLake_r%lake_outflow(iunit), TLake_r%lake_evap(iunit), TLake_r%lake_rain(iunit) + TLake_r%lake_snow(iunit)
        !    end if
        !    ! gulu

    !do i=1, TLake_r%d_ns(iunit)
    !if(.not.(isnan(t_in)) .and. isnan(TLake_r%temp_lake(iunit,i))) then
    !    write(unit=1899,fmt="(i10, i4, i4, i4, 2(e14.6))") iunit, 2, i, TLake_r%d_ns(iunit), t_in, TLake_r%temp_lake(iunit,i)
    !    call shr_sys_abort('Goofy temperature -- NAN values #23')
        
    !end if
    !end do

    !if(iunit == 198841) then 
    !!    write(unit=21002,fmt="(i4, i4, 6(e18.10))") 2, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), t_in, THeat%forc_t(iunit), TLake_r%temp_lake(iunit,1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit))
    !end if
            end do ! sub-timestep

            tmp_outflow = tmp_outflow/s_dtime
            TRunoff%erout(iunit,nt_nliq) = tmp_outflow
            tmp_tout = tmp_tout/s_dtime
            TLake_r%lake_Tout(iunit) = tmp_tout
            THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nt_nliq)), TLake_r%lake_Tout(iunit))
            tmp_evap = tmp_evap/s_dtime
            TLake_r%lake_evap(iunit) = -tmp_evap
            
            TLake_r%lake_Tsur(iunit) = TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit))
            if(TLake_r%lake_Tsur(iunit) > 400._r8) then
                !write(unit=10009,fmt="(i10, 5(e14.6))") iunit, TLake_r%d_lake(iunit), t_in, THeat%forc_t(iunit), TLake_r%lake_Tsur(iunit), TLake_r%temp_lake(iunit,1)
                write(iulog,*) 'Lake temperature too high in r-zone! ', iunit, TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), t_in, THeat%forc_t(iunit), TLake_r%lake_Tsur(iunit), TLake_r%temp_lake(iunit,1)
                !call shr_sys_abort('mosart: water temperature: '//subname)
            end if
            
            ! the portion of local tributary flow that does not enter the lake will enter back the tributary channel
            !TLake_r%lake_Tout(iunit) = t_out
            !THeat%ha_lateral(iunit) = cr_advectheat(abs(TRunoff%erlateral(iunit,nt_nliq)+TRunoff%erlateral(iunit,nt_nice)), THeat%Tr(iunit))
            !if(iunit == 20975) then  
    !        if(iunit == 186781) then 
            !if(iunit == 138309) then
                !write(unit=51000,fmt="(i4, (e12.4))") TLake_r%d_ns(iunit), ddz_max
            !    write(unit=11000,fmt="(i4, 5(e12.4))") TLake_r%d_ns(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%d_lake(iunit), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), THeat%forc_t(iunit)
    !            write(unit=11001,fmt="(i4, i4, i4, i4, 8(e18.10))") 102, ww, dflag, TLake_r%d_ns(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1), TLake_r%d_lake(iunit), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%d_v(iunit,TLake_r%d_ns(iunit)-1), TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)), t_in, THeat%forc_t(iunit)
    !            write(unit=21001,fmt="(i4, i4, 9(e18.10))") dflag, TLake_r%d_ns(iunit), TRunoff%yt(iunit,1), TLake_r%d_lake(iunit), TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1), TLake_r%lake_inflow(iunit), TLake_r%lake_rain(iunit) + TLake_r%lake_snow(iunit), TRunoff%erlateral(iunit,nt_nliq), TLake_r%lake_outflow(iunit), TLake_r%lake_evap(iunit), TUnit_lake_r%F_local(iunit)
    !        end if
            !if(iunit == 186781) then
            
            !if(TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)) - THeat%forc_t(iunit) < 10._r8) then
            !    if(TLake_r%d_ns(iunit) >=2) then
            !        write(unit=8000,fmt="(i4, 3(e18.10))") TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1)
            !    else
            !        write(unit=8001,fmt="(i4, 2(e18.10))") TLake_r%d_ns(iunit), TLake_r%d_lake(iunit), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
            !    end if
            !end if
            !if(iunit == 62863) then 
            !    write(unit=21001,fmt="(i4, 7(e18.10))") TLake_r%d_ns(iunit), TRunoff%yr(iunit,1), TLake_r%d_lake(iunit), TLake_r%lake_Tout(iunit), THeat%Tr(iunit), TLake_t%lake_Tout(iunit), THeat%Tt(iunit), THeat%forc_t(iunit)
            !end if
            
        end if ! stratification
    end subroutine mosart_lake_r


    function lake_evap_r(iunit, dtime, t_sur, F, es, sar) result(evap)
    ! calculate lake surface evaporation
        implicit none
        integer,  intent(in) :: iunit
        real(r8), intent(in) :: dtime  ! local time step (s)
        real(r8), intent(in) :: t_sur  ! lake surface temperature (K)
        real(r8), intent(in) :: F      ! dimensionless factor for wind sheltering by riparian vegetation, Wu et al, 2012
        real(r8), intent(in) :: es     ! Satuated vapor pressure(Pa)
        real(r8), intent(in) :: sar    ! Surface area ratio
        real(r8)             :: evap   ! evaporation rate (m/s)
        
        real(r8) :: kl                 ! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
        real(r8) :: d_evap,v_evap      ! Evaporated depth (m) and volume (m^3)
        
        kl         = 0.211_r8 + 0.103_r8 * THeat%forc_wind(iunit) * F               
        evap = 0._r8
        d_evap = 0._r8
        v_evap = 0._r8
        if(TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1) - TLake_r%v_zti(iunit,1) < TINYVALUE) then
            evap = 0._r8
            d_evap = 0._r8
            v_evap = 0._r8
        else
            evap     = max(kl*(es - THeat%forc_vp(iunit))/100._r8,0._r8)  ! in mm/d
            evap     = evap/(86.4e6)  ! mm/d --> m/s
            
            d_evap   = evap*dtime
            v_evap   = max(d_evap*sar*TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1),0._r8)
            ! TODO. For now, the layer adjustment should ensure that evaporation occurs most of time
            ! in a extreme dry case, we don't let the lake to completely dry out for numerical stability
            if (v_evap - 0.95_r8 * TLake_r%d_v(iunit, TLake_r%d_ns(iunit)) > TINYVALUE) then ! too little lake storage for lake evap. to occur fully
                !write(unit=12202,fmt="(i10, 2(e14.6))") iunit, v_evap, TLake_r%V_str(iunit)
                v_evap = TLake_r%d_v(iunit, TLake_r%d_ns(iunit)) * 0.95_r8
            end if
            
            if(v_evap < 0._r8) then
                v_evap = 0._r8
            end if
            !if(v_evap - (TLake_r%v_zt(iunit, TLake_r%d_ns(iunit)+1) - TLake_r%v_zti(iunit, 1)) > TINYVALUE) then
            !    v_evap = TLake_r%v_zt(iunit, TLake_r%d_ns(iunit)+1) - TLake_r%v_zti(iunit, 1)
            !end if
                                
            d_evap = v_evap/TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)            
            evap = d_evap/dtime
        end if
        ! no evap if freezing surface
        if(t_sur <= 273.15_r8) then
            evap = 0._r8
        end if
        
    end function lake_evap_r

    SUBROUTINE PDEPTH_r(iunit,RHOIN,Q,S,FT,WIDTH,QENIN,IHP,HP)
    !
    !  COMPUTE DEPTH AT PLUNGE POINT AND INITIAL ENTRAINMENT AT
    !  THE PLUNGE POINT (GAMAIN). (BY AKVAMA)
    !
    
        implicit none
        integer,  intent(in) :: iunit   ! Index of local grid
        real(r8), intent(in) :: RHOIN   ! INITIAL DENSITY OF INFLOW (KG/M3)
        real(r8), intent(in) :: Q       ! Advective inflow (M3/S)
        real(r8), intent(in) :: S       ! SLOPE OF INFLOW CHANNEL (-)
        real(r8), intent(in) :: FT      ! MANNINGS FRICTION FACTOR OF INFLOW CHANNEL (-) *** I assume its actual unit is based on the SI unit system, not US customary units?
        real(r8), intent(in) :: WIDTH   ! WIDTH OF INFLOW CHANNEL (M) 
        real(r8), intent(out):: QENIN   ! INITIAL total ENTRAINMENT AT PLUNGING (M3/s)
        integer,  intent(out):: IHP     ! Index of deepest layer of the plunging point
        real(r8), intent(out):: HP      ! DEPTH AT PLUNGE POINT (M) 
    
        real(r8) :: RHOMIX  ! DENSITY IN THE MIXED ZONE, so the same as top layer (KG/M3)
        real(r8) :: SUMZ    ! Dummy variable (M)
        real(r8) :: EPSIN   ! relative density difference (-)
        real(r8) :: GAMAIN  ! 
        real(r8) :: X       ! dummy variable
        integer  :: I       ! local index
    
        
        RHOMIX = den(TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)))    
        EPSIN=(RHOIN-RHOMIX)/RHOMIX
        IF(EPSIN.LE.0._r8) then ! incoming flow lighter than any layers, no entrainment
            IHP   = TLake_r%d_ns(iunit)
            HP    = 0._r8
            QENIN = 0._r8
 !if(iunit == 186781) then
    ! write(unit=3001,fmt="(i10, i4, 4(e14.6))") iunit, TLake_r%d_ns(iunit), RHOIN, RHOMIX, EPSIN
 !end if
            RETURN
        end if
        
        X = Q/WIDTH
        IF(S.LE.6.66667E-3) THEN
           HP = 1.1*(FT*X*X/(S*EPSIN*9.81))**(.3333)
           GAMAIN=0.15
        ELSE
           HP=1.6*(X*X/(EPSIN*9.81))**(.3333)
           GAMAIN= 1.8
        ENDIF
        QENIN=GAMAIN*Q
    
    ! write(unit=3002,fmt="(i4, i10, i4, i4, 5(e14.6))") 1, iunit, TLake_r%d_ns(iunit), IHP, GAMAIN, QENIN, Q, SUMZ, HP
        SUMZ=0.
        IHP=TLake_r%d_ns(iunit)+1
        DO I=TLake_r%d_ns(iunit),1,-1 ! loop from top to bottom layer
           IHP = IHP-1
           SUMZ=SUMZ+TLake_r%dd_z(iunit,I)
           IF(SUMZ.GE.HP) then 
              exit
           end if
        end do  
        HP=SUMZ

    ! write(unit=3002,fmt="(i4, i10, i4, i4, 5(e14.6))") 2, iunit, TLake_r%d_ns(iunit), IHP, GAMAIN, QENIN, Q, SUMZ, HP

        RETURN
        
    END SUBROUTINE PDEPTH_r

    FUNCTION ENTRAIN(I,DCF,RDC,RHOAMB,DELZ,S,SUMZ,Q,IHP,WIDTH,FT,iunit) result(entrain_)
    !
    ! COMPUTE THE ENTRAINMENT FROM A LAYER INTO THE
    ! DENSITY CURRENT (FROM AKIYAMA)
    !
       implicit none
       integer,  intent(in) :: I,iunit       ! Index of the layer affected by entrainment
       real(r8), intent(in) :: DCF     ! Incoming flow IN THE DENSITY CURRENT (M3/s)
       real(r8), intent(in) :: RDC     ! DENSITY OF incoming DENSITY CURRENT (KG/M3)
 
       real(r8), intent(in) :: RHOAMB  ! AMBIENT (current, existing) DENSITY IN A LAYER (KG/M3)
       real(r8), intent(in) :: DELZ    ! Layer thickness (M)
       real(r8), intent(in) :: S       ! SLOPE OF INFLOW CHANNEL (-)
       real(r8), intent(in) :: SUMZ    ! DEPTH AT PLUNGE POINT (M)
       real(r8), intent(in) :: Q       ! Advective inflow (M3/s)
       integer,  intent(in) :: IHP     ! Index of deepest layer of the plunging point
       real(r8), intent(in) :: WIDTH   ! WIDTH OF INFLOW CHANNEL (M) 
       real(r8), intent(in) :: FT      ! MANNINGS FRICTION FACTOR OF INFLOW CHANNEL (-)
       real(r8)             :: entrain_! entrainment of current layer [m3/s]
 
       
       real(r8)             :: EPSI    ! relative density difference [-]
       real(r8)             :: FD      ! *** physical meaning and unit?
       real(r8)             :: F43     ! *** physical meaning and unit?
       real(r8)             :: X       ! local dummy variable
       real(r8)             :: GAMAI   ! *** physical meaning and unit?        
       
       if(SUMZ < TINYVALUE .or. S < TINYVALUE) then
           entrain_ = 0._r8
           return
       end if
       IF(I.GT.IHP) THEN
         EPSI=(RDC-RHOAMB)/RHOAMB
         if(EPSI < 0._r8) then
             entrain_ = 0._r8
             return
         end if
         
         !IF(I.EQ.(IHP+1)) THEN  !*** not sure whether I've translated the GO statement correctly, perhaps not
             FD=1.875E-4+FT
             F43=( (FD+SQRT(FD*FD+0.0045*S))/( 1.5*S) )**( 1.3333)
         !ENDIF
         X=DCF/WIDTH
         GAMAI=0.0015*DELZ*(9.81*EPSI/(X*X))**(.3333)/(F43*S)
         entrain_ = GAMAI*DCF
       ELSE
         entrain_ = Q*DELZ/SUMZ
       ENDIF
       
       RETURN
    END FUNCTION ENTRAIN
 
 
    subroutine av_in_r(d_n,in_f,in_t,dtime,iunit,dv_in)
!*******************************************************************************************************
!     calculte the layer volume and temperature changes due to advective inflow and its entrainment effects
!*******************************************************************************************************    
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        implicit none
        integer, intent(in)   :: d_n             ! number of effective layers
        integer, intent(in)   :: iunit           ! index of current grid
        real(r8), intent(in)  :: dtime           ! local time step (s)
        real(r8), intent(in)  :: in_f            ! initial inflow rate (without accounting for entrainment yet) (m3/s) 
        real(r8), intent(in)  :: in_t            ! initial inflow temperature (K)
        real(r8),dimension(nlayers), intent(out) :: dv_in    ! layer inflow volume (m3) 
        real(r8),dimension(nlayers):: rho_layers ! water density at each layer (kg/m3)
        real(r8),dimension(nlayers):: QE         ! entrainment at each layer (m3/s)
        real(r8)  :: rho_in                      ! initial inflow density (kg/m3)
        real(r8)  :: in_v                        ! total inflow volume (m3)
        integer   :: j, k                        ! indices
        integer   :: iplayer                     ! index of the layer receiving planging inflow
        real(r8)  :: QENIN                       ! INITIAL total ENTRAINMENT AT PLUNGING (M3/s)
        integer   :: IHP                         ! Index of deepest layer of the plunging point
        real(r8)  :: HP                          ! DEPTH AT PLUNGE POINT (M) 
        real(r8)  :: in_f_new                    ! updated inflow rate after accounting for entrainment (m3/s) 
        real(r8)  :: in_t_new                    ! updated inflow temperature (K)
        real(r8)  :: RDCF                        ! dummy variable

    !    Initialize
        rho_in = den(in_t)
        dv_in = 0._r8
        rho_layers = 0._r8
        do j=1,d_n
            rho_layers(j) = den(TLake_r%temp_lake(iunit,j))
        end do
        QE = 0._r8
        in_f_new = in_f
        in_t_new = in_t
        
 !if(iunit == 186781) then
 !    write(unit=2001,fmt="(i4, i4, 3(e14.6))") 0, d_n, sum(QE), sum(TLake_r%d_v), QENIN
 !end if
 !if(iunit == 186781) then
     !write(unit=2001,fmt="(i4, i4, 3(e14.6))") 1, d_n, sum(QE), sum(TLake_r%d_v), QENIN
 !    write(unit=2011,fmt="(i4, i4, i4, 7(e14.6))") 1, iplayer, d_n, sum(TLake_r%d_v(iunit,:)), in_f* dtime, in_f_new * dtime, QENIN * dtime, in_t, rho_in, in_t_new 
 !end if
        
        if(d_n>1) then
            ! determine IHP and QENIN
            call PDEPTH_r(iunit,rho_in,in_f,TUnit%rslp(iunit),TUnit%nr(iunit),TUnit%rwidth(iunit),QENIN,IHP,HP)
            in_v = in_f*dtime
            
            ! calculate entrainment due to plunging            
            if(QENIN > TINYVALUE) then
            do j=d_n,1,-1
                QE(j) =  ENTRAIN(j,in_v,rho_in,rho_layers(j),TLake_r%dd_z(iunit,j),TUnit%rslp(iunit),HP,QENIN,IHP,TUnit%rwidth(iunit),TUnit%nr(iunit),iunit)
                if(QE(j)*dtime > 0.95_r8*TLake_r%d_v(iunit,j)) QE(j) = 0.95_r8 * TLake_r%d_v(iunit,j) / dtime
                RDCF = 1/(in_f + QE(j))
 !if(iunit == 186781) then
 !    write(unit=2023,fmt="(i4, i4, 5(e14.6))") 1, j, in_t, in_f, TLake_r%temp_lake(iunit,j), QE(j), RDCF 
 !end if
                in_t_new = (in_t * in_f + TLake_r%temp_lake(iunit,j) * QE(j)) * RDCF
                in_f_new = in_f_new + QE(j)
 !if(iunit == 186781) then
 !    write(unit=2023,fmt="(i4, i4, 5(e14.6))") 2, j, in_t, in_f, TLake_r%temp_lake(iunit,j), QE(j), RDCF 
 !end if
                
                TLake_r%d_v(iunit,j) = TLake_r%d_v(iunit,j) - QE(j) * dtime
 !if(iunit == 508 .and. QE(j) > TINYVALUE) then
 !    write(unit=9001,fmt="(i4, i4, i4, 4(e14.6))") j, d_n, IHP, TLake_r%d_v(iunit,j), in_f*dtime, QE(j)*dtime, QENIN*dtime
 !end if
            end do
            end if            
 !    write(unit=9002,fmt="(i10, i4, (e14.6))") iunit, d_n, QENIN
        end if
 !if(iunit == 186781) then
 !    write(unit=2001,fmt="(i4, i4, 5(e14.6))") 1, d_n, QE(1), QE(2), TLake_r%d_v(iunit,1), TLake_r%d_v(iunit,2), QENIN
 !end if
        
!   Layer inflow based on density, assuming inflow itself has an uniform density will only enter into one layer instead of multiple
        iplayer = d_n ! initialize
        if(rho_in >= rho_layers(1)) then     ! inflow density larger than any layers, enter into the bottom layer
            iplayer = 1
        elseif(rho_in <= rho_layers(d_n)) then! inflow density less than any layers, enter into the top layer
            iplayer = d_n
        else
            do k=d_n-1, 1, -1
                if(rho_in >= rho_layers(k+1) .and. rho_in < rho_layers(k)) then
                    iplayer = k+1
                    exit
                end if                    
            end do
        end if           
        if(iplayer > d_n) iplayer = d_n
 !if(iunit == 186781) then
 !    write(unit=2013,fmt="(i4, i4, 6(e14.6))") 1, iplayer, RDCF, TLake_r%temp_lake(iunit,iplayer), TLake_r%d_v(iunit,iplayer), in_t_new, in_f_new * dtime, in_t 
 !end if
        
        RDCF = 1/(TLake_r%d_v(iunit,iplayer) + in_f_new * dtime)        
        TLake_r%temp_lake(iunit,iplayer) = (TLake_r%temp_lake(iunit,iplayer) * TLake_r%d_v(iunit,iplayer) + in_t_new * in_f_new * dtime) * RDCF
        TLake_r%d_v(iunit,iplayer) = TLake_r%d_v(iunit,iplayer) + in_f_new * dtime
        dv_in = 0._r8  ! TODO. Set it as zero, since layer volume change has been taken care of here already
 !if(iunit == 186781) then
     !write(unit=2001,fmt="(i4, i4, 3(e14.6))") 3, d_n, sum(QE), sum(TLake_r%d_v), QENIN
     !write(unit=2011,fmt="(i4, i4, i4, 4(e14.6))") 3, iplayer, d_n, sum(TLake_r%d_v(iunit,:)), in_f* dtime, in_f_new * dtime, QENIN * dtime 
 !    write(unit=2013,fmt="(i4, i4, 6(e14.6))") 2, iplayer, RDCF, TLake_r%temp_lake(iunit,iplayer), TLake_r%d_v(iunit,iplayer), in_t_new, in_f_new * dtime, in_t 
 !end if

                
    end subroutine av_in_r

    subroutine av_out_r(d_n,ou_f,dtime,iunit,dv_ou,t_out)
!*******************************************************************************************************
!     advective outflow
!*******************************************************************************************************    
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        implicit none
        integer, intent(in)  :: d_n             ! number of effective layers
        integer, intent(in)  :: iunit           ! index of current grid
        real(r8), intent(in) :: dtime           ! local time step (s)
        real(r8), intent(in) :: ou_f            ! outflow rate (m3/s)
        real(r8), dimension(nlayers), intent(out) :: dv_ou    ! layer outflow volume (m3) 
        real(r8), intent(out) :: t_out           ! average outflow temperature (K) 
        real(r8),dimension(nlayers):: rho_layers
        real(r8) :: ou_v                        ! total outflow volume (m3)
        real(r8) :: denom                       ! dummy variable
        integer :: j, k                          ! indices
        integer :: jmin, jmax                   ! range of indices for outflowing layers
    
        ou_v = ou_f*dtime

 !if(iunit == 186781) then
 !    write(unit=2001,fmt="(i4, i4, 3(e14.6))") 1, d_n, TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), sum(dv_ou), ou_v
 !end if
        
    !    Initialize
        do j=1,nlayers
            dv_ou(j)=0._r8
            rho_layers(j) = 0._r8
        end do
        do j=1,d_n
            rho_layers(j) = den(TLake_r%temp_lake(iunit,j))
        end do
                
        jmin=TLake_r%J_Min(iunit)
        !if(d_n>3)jmax=d_n-2 
        !if(d_n<=3)jmax=d_n-1
        jmax = d_n
        
        if(d_n==1) then ! outflow from single layer
            dv_ou(d_n) = ou_v
        else            ! outflow from multiple layers, using volume instead of thickness as weight to avoid any water balance error
            denom=sum(TLake_r%d_v(iunit,jmin:jmax))
            do j=jmin,jmax
                dv_ou(j) = ou_v*(TLake_r%d_v(iunit,j)/(denom))
                
 !   if(dv_ou(j) >= TLake_r%d_v(iunit,j)) then
 !       write(unit=2002,fmt="(i10, 4(i4), 4(e14.6))") iunit, j, jmax, jmin, TLake_r%d_ns(iunit), dv_ou(j), TLake_r%d_v(iunit,j), ou_v, denom
 !   end if
            end do
        end if

 !if(iunit == 186781) then
 !    write(unit=2002,fmt="(i4, i4, 2(e14.6))") jmax,jmin,dv_ou(1),dv_ou(2)
 !end if 

        t_out = 0._r8
        denom = 0._r8
        do j=jmin,jmax
            denom = denom + dv_ou(j) * TLake_r%temp_lake(iunit,j)
        end do
        t_out = denom / ou_v

  !if(isnan(t_out)) then
  !    write(unit=1088,fmt="(i10, 3(i4), 2(e14.6))") iunit, d_n, jmin,jmax, denom, ou_v
  !end if
        
    end subroutine av_out_r
       
    subroutine under_flow_r(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
!*******************************************************************************************************
!     lake is emptying or overflowing
!*******************************************************************************************************    
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        
        implicit none
        integer, intent(in) :: iunit,d_n_n              
        real(r8),intent(in) :: d_res_sub              
        real(r8),dimension(nlayers), intent(in) :: t_z_old 
        real(r8),dimension(nlayers), intent(out) :: dv_in,dv_ou    ! layer inflow,outflow(m3/s)
        character(len=*),parameter :: subname = '(over_under_flow_r)'
        integer :: j                            ! indices
        
        write(iulog,*) 'Warning under_flow_r occured ! ', iunit, TUnit_lake_r%V_max(iunit)*1e6, TLake_r%d_ns(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%v_zti(iunit,1), TLake_r%v_zti(iunit,ngeom+1)
        !call shr_sys_abort('mosart: imbalance in MOSART-lake: '//subname)

!if(0>1) then
        TLake_r%d_ns(iunit)      = d_n_n                
        TLake_r%d_lake(iunit)    = d_res_sub 
        do j = 1, nlayers+1    
            if (j<=nlayers)  then   
                TLake_r%d_v(iunit,j)         = TLake_r%v_zo(iunit,j)
                TLake_r%temp_lake(iunit,j) = t_z_old(j)
                TLake_r%v_zt(iunit,j)     = TLake_r%v_zt0(iunit,j)
                TLake_r%a_d(iunit,j)     = TLake_r%a_d0(iunit,j)
                TLake_r%d_z(iunit,j)     = TLake_r%d_z0(iunit,j)
                TLake_r%dd_z(iunit,j)    = TLake_r%d_z0(iunit,j+1)-TLake_r%d_z0(iunit,j)
                dv_in(j)    = 0._r8
                dv_ou(j)    = 0._r8
            else
                TLake_r%v_zt(iunit,j)     = TLake_r%v_zt0(iunit,j)
                TLake_r%a_d(iunit,j)     = TLake_r%a_d0(iunit,j)
                TLake_r%d_z(iunit,j)     = TLake_r%d_z0(iunit,j)
            end if
        end do
!end if

    end subroutine under_flow_r

    subroutine over_flow_r(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
!*******************************************************************************************************
!     lake is emptying or overflowing
!*******************************************************************************************************    
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        
        implicit none
        integer, intent(in) :: iunit,d_n_n              
        real(r8),intent(in) :: d_res_sub              
        real(r8),dimension(nlayers), intent(in) :: t_z_old 
        real(r8),dimension(nlayers), intent(out) :: dv_in,dv_ou    ! layer inflow,outflow(m3/s)
        character(len=*),parameter :: subname = '(over_under_flow_r)'
        integer :: j                            ! indices
        
        write(iulog,*) 'Warning over_flow_r occured ! ', iunit, TUnit_lake_r%V_max(iunit)*1e6, TLake_r%d_ns(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%v_zti(iunit,1), TLake_r%v_zti(iunit,ngeom+1)
        !call shr_sys_abort('mosart: imbalance in MOSART-lake: '//subname)

!if(0>1) then
        TLake_r%d_ns(iunit)      = d_n_n                
        TLake_r%d_lake(iunit)    = d_res_sub 
        do j = 1, nlayers+1    
            if (j<=nlayers)  then   
                TLake_r%d_v(iunit,j)         = TLake_r%v_zo(iunit,j)
                TLake_r%temp_lake(iunit,j) = t_z_old(j)
                TLake_r%v_zt(iunit,j)     = TLake_r%v_zt0(iunit,j)
                TLake_r%a_d(iunit,j)     = TLake_r%a_d0(iunit,j)
                TLake_r%d_z(iunit,j)     = TLake_r%d_z0(iunit,j)
                TLake_r%dd_z(iunit,j)    = TLake_r%d_z0(iunit,j+1)-TLake_r%d_z0(iunit,j)
                dv_in(j)    = 0._r8
                dv_ou(j)    = 0._r8
            else
                TLake_r%v_zt(iunit,j)     = TLake_r%v_zt0(iunit,j)
                TLake_r%a_d(iunit,j)     = TLake_r%a_d0(iunit,j)
                TLake_r%d_z(iunit,j)     = TLake_r%d_z0(iunit,j)
            end if
        end do
!end if

    end subroutine over_flow_r
     
    subroutine lake_r_waterbalance_check(iunit) 
    ! !DESCRIPTION: check the lake water balance
        implicit none    
        integer, intent(in) :: iunit
        character(len=*),parameter :: subname = '(lake water balance)'
      
        integer  :: i
        real(r8) :: v_sum ! total lake volume [m^3]
        real(r8) :: V_min ! local threshold for non-negligible lake storage
        real(r8) :: myTINYVALUE ! 
        
        V_min = 1000000._r8  ! relative error for very small lake volumes could be much higher than large lake volumes
        myTINYVALUE = 1e-5_r8
        
        v_sum = TLake_r%v_zt(iunit, 1)
        do i=1,TLake_r%d_ns(iunit)
           if (TLake_r%d_v(iunit,i) < -myTINYVALUE) then
              write(iulog,*) 'Type 1 error in r-zone lake water balance check ! ', iunit, i, TLake_r%d_ns(iunit), TLake_r%d_v(iunit,i)
              call shr_sys_abort('mosart: negative storage in MOSART-lake: '//subname)           
           end if
           v_sum = v_sum + TLake_r%d_v(iunit,i)
        end do
        if (TLake_r%V_str(iunit) < -myTINYVALUE) then
            write(iulog,*) 'Type 2 error in r-zone lake water balance check ! ', iunit, TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) +1)
            call shr_sys_abort('mosart: negative storage in MOSART-lake: '//subname)
        end if        

        if (v_sum > V_min .and. abs(TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1) - v_sum)/v_sum > myTINYVALUE) then
            write(iulog,*) 'Type 3 error in r-zone lake water balance check ! ', iunit, TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit)), v_sum, v_sum + TLake_r%v_zt(iunit,1)
            call shr_sys_abort('mosart: imbalance in MOSART-lake: '//subname)
        end if        

        if (v_sum > V_min .and. abs(TLake_r%V_str(iunit) - v_sum)/v_sum > myTINYVALUE) then
            write(iulog,*) 'Type 4 error in r-zone lake water balance check ! ', iunit, TLake_r%V_str(iunit), v_sum
            call shr_sys_abort('mosart: imbalance in MOSART-lake: '//subname)
        end if                       

        if (TLake_r%V_str(iunit) > V_min .and. abs(TLake_r%v_zt(iunit, TLake_r%d_ns(iunit)+1) - TLake_r%V_str(iunit))/TLake_r%V_str(iunit) > myTINYVALUE) then
            write(iulog,*) 'Type 5 error in r-zone lake water balance check ! ', iunit, TLake_r%v_zt(iunit, TLake_r%d_ns(iunit)+1), TLake_r%V_str(iunit), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit)+1) - TLake_r%V_str(iunit), TLake_r%lake_inflow(iunit)
            call shr_sys_abort('mosart: imbalance in MOSART-lake: '//subname)
        end if                       

    end subroutine lake_r_waterbalance_check

    subroutine nlayer2ngeom_r(iunit)
!*******************************************************************************************************
!     convert the nlayer-structure of lake layers into ngeom-structure
!*******************************************************************************************************    
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        
        implicit none
        integer, intent(in) :: iunit              
        character(len=*),parameter :: subname = '(nlayer2ngeom_r)'
        integer :: n_ngeom  ! number of effective layers in the ngeom frame
        real(r8), pointer :: v_local(:)  ! discretized (cumulative) layer volume in the ngeom frame
        real(r8) :: dd_in
        real(r8) :: ztop_ng, zbot_ng    ! top and bottom elevation of a layer in the ngeom frame (m)
        real(r8) :: ztop_nl1, ztop_nl2, zbot_nl    ! top and bottom elevation of a layer in the ngeom frame (m)
        real(r8) :: dv_in(nlayers),dv_out(nlayers)  ! volume increment/decrease at layer due to inflow/outflow(m^3)
        real(r8) :: enr_nl(nlayers)  ! inner energy (J)
        real(r8) :: c_w = 4.188e3                                                ! Water specific heat capacity  (J/kg/k)
        real(r8) :: rho_w = 1.e3                                                ! Water density  (kg/m3)
        real(r8) :: num,dem,delta_z,delta_a    !
        real(r8) :: delta_v1, delta_v2, delta_e1, delta_e2
        !real(r8) :: tab,e_ab,dd_zab,d_vab,dv_ouab,dv_inab,delta_z    
        integer :: ing, inl,j,k,mm                            ! indices
        
        allocate(v_local(ngeom))
        v_local = 0._r8
        dv_in = 0._r8 
        dv_out = 0._r8
        !rtmCTL%H_lake(iunit,:) = 0._r8
        !rtmCTL%T_lake(iunit,:) = 0._r8
        
        !rtmCTL%H_lake(iunit,1) = TLake_r%d_z(iunit,1)
        !rtmCTL%T_lake(iunit,1) = TLake_r%temp_lake(iunit,1)
        
        if(TLake_r%d_ns(iunit) <=1) then
            write(iulog,*) 'Warning -- no active lake layer, '//subname, iunit, TLake_r%d_ns(iunit)
            return
        end if
        
        ! first determine the effective depth of each layer in the ngeom frame
        n_ngeom = 2
        do ing=2, ngeom+1
            if((TLake_r%d_lake(iunit) - TLake_r%d_zi(iunit, ing-1) > TINYVALUE) .and. (TLake_r%d_lake(iunit) - TLake_r%d_zi(iunit, ing) <= -TINYVALUE)) then
                n_ngeom = ing
                exit
            end if
        end do
        dd_in = TLake_r%d_zi(iunit,2) - TLake_r%d_zi(iunit,1)        
        do ing=2, n_ngeom-1
            !rtmCTL%H_lake(iunit,ing) = rtmCTL%H_lake(iunit,ing-1) + dd_in
        end do
        !rtmCTL%H_lake(iunit,n_ngeom) = TLake_r%d_lake(iunit) - TLake_r%d_zi(iunit, n_ngeom-1) ! it is ok for this top layer depth to be small since it is not used for calculation
        
        ! now check the layer depth in the nlayer frame, to make each layer no less than dd_in
        do inl=2, TLake_r%d_ns(iunit)+1
            if(TLake_r%d_z(iunit, inl) - dd_in < - TINYVALUE) then
                !call layer_merge_r()
            end if
        end do
        
        ! now determine the effective volume of each layer (cumulative) and temperature in the ngeom frame
        v_local(1) = TLake_r%v_zti(iunit,1)
        do ing=2, n_ngeom
    !        ztop_ng = rtmCTL%H_lake(iunit,ing)
    !        zbot_ng = rtmCTL%H_lake(iunit,ing-1)
            do inl=2, TLake_r%d_ns(iunit)
                ztop_nl2 = TLake_r%d_z(iunit, inl+1)
                ztop_nl1 = TLake_r%d_z(iunit, inl)
                zbot_nl = TLake_r%d_z(iunit, inl-1)
                if((ztop_ng - ztop_nl1 < -TINYVALUE) .and. (zbot_ng - zbot_nl > TINYVALUE)) then ! the ngeom layer is complete within the middle of a nlayer layer
                    v_local(ing) = v_local(ing-1) + (TLake_r%v_zti(iunit, ing) - TLake_r%v_zti(iunit, ing-1))
    !                rtmCTL%T_lake(iunit,ing) = TLake_r%temp_lake(iunit, inl)
                    
                    exit
                end if
                if((ztop_ng - ztop_nl2 < -TINYVALUE) .and. (ztop_ng - ztop_nl1 > TINYVALUE) .and. (zbot_ng - zbot_nl > TINYVALUE)) then ! the ngeom layer is in the interface of two nlayer layers
                    delta_z = (ztop_ng-ztop_nl1)/(ztop_nl2-ztop_nl1)
                    delta_v2 = delta_z * TLake_r%v_zt(iunit,inl+1)
                    delta_e2 = delta_v2 * TLake_r%temp_lake(iunit, inl+1)
                    delta_z = (ztop_nl1-zbot_ng)/(ztop_nl1-zbot_nl)
                    delta_v1 = delta_z * TLake_r%v_zt(iunit,inl)
                    delta_e1 = delta_v1 * TLake_r%temp_lake(iunit, inl)
                    v_local(ing) = v_local(ing-1) + delta_v2 + delta_v1
    !                rtmCTL%T_lake(iunit, ing) = (delta_e2+delta_e1)/(delta_v2+delta_v1)
                    exit
                end if
                
                if((ztop_ng - ztop_nl1 > TINYVALUE) .and. (zbot_ng - zbot_nl < -TINYVALUE)) then ! the ngeom layer is completely including a nlayer layer
                    
                    
                    
                    exit
                end if
                
            
            end do
        end do
       
    end subroutine nlayer2ngeom_r


    subroutine solve_r(a,b,c,r,iunit,n) !
    ! subroutine to solve the tri-diagnal matrix. 
    ! Note that, here the lake layer numbers increase from bottom to top, similar to CE-QUAL-R1
        use shr_sys_mod , only : shr_sys_flush
    
        implicit none
!     a:left coefficient , b:middle coefficient, c:right coefficient, r:right hand side (known), n - index of top layer, aka, effective number of layers
        integer,intent(in) :: n,iunit
        real(r8),dimension(n),intent(in) :: a,b,c,r
        real(r8),dimension(n) :: bp,rp
        real(r8) :: m,tt 
        integer i        
        
!     initialize c-prime and d-prime
        bp(n)=b(n)
        rp(n)=r(n)
        do i=n-1,1,-1 
            tt=a(i)/b(i+1)
            bp(i) = b(i)-c(i+1)*tt
            rp(i) = r(i)-rp(i+1)*tt
        end do
       
!     Back substitution
        if(abs(bp(1)) > TINYVALUE) then
            TLake_r%temp_lake(iunit,1) = rp(1)/bp(1)
        else
            TLake_r%temp_lake(iunit,1) = 273.15_r8
        end if
        if(TLake_r%temp_lake(iunit,1) < 273.15_r8) then ! to be modified after including ice/snow effect
            TLake_r%temp_lake(iunit,1) = 273.15_r8
        end if
!     Back substitution
        do i = 2,n
            TLake_r%temp_lake(iunit,i) = (rp(i)-c(i)*TLake_r%temp_lake(iunit,i-1))/bp(i)
            if(TLake_r%temp_lake(iunit,i) < 273.15_r8) then ! to be modified after including ice/snow effect
                TLake_r%temp_lake(iunit,i) = 273.15_r8
            end if
            ! freezing do not start from the bottom
            !TLake_r%temp_lake(iunit,i) = (rp(i)-c(i)*TLake_r%temp_lake(iunit,i+1))/bp(i)
            if(TLake_r%temp_lake(iunit,i) <= 273.15_r8 .and. TLake_r%temp_lake(iunit,i+1) > 273.15_r8) then
                TLake_r%temp_lake(iunit,i) = TLake_r%temp_lake(iunit,i+1)
            end if
        end do
        
    !  Check numerical instability from reservoir geometry data     
        do i = 1,TLake_r%d_ns(iunit)
            if (isnan(TLake_r%temp_lake(iunit,i)) .or. TLake_r%temp_lake(iunit,i)>=huge(1._r8)) then
                !TLake_r%temp_lake(iunit,i) = THeat%Tr(iunit)
                write(iulog,*) 'error in lake temperature calculation! ', iunit, i, TLake_r%d_ns(iunit), TLake_r%temp_lake(iunit,i), TLake_r%d_lake(iunit)
                call shr_sys_abort('mosart-lake: unreasonable temperature value')
            end if    
        end do
        
    end subroutine solve_r 
 
  
end MODULE mosart_lake_r_mod