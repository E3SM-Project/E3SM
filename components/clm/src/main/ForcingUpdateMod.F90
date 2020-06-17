module ForcingUpdateMod

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use atm2lndType      , only : atm2lnd_type
  use TopounitDataType , only : top_af, top_as
  use GridcellType     , only : grc_pp
  use shr_const_mod    , only : SHR_CONST_TKFRZ, SHR_CONST_STEBOL
  use clm_varctl       , only: const_climate_hist, add_temperature, add_co2, use_cn
  use clm_varctl       , only: startdate_add_temperature, startdate_add_co2
  use clm_varctl       , only: co2_type, co2_ppmv,  use_c13, create_glacier_mec_landunit

  use clm_varcon       , only: rair, o2_molar_const, c13ratio
  use decompMod        , only : bounds_type
  use domainMod        , only : ldomain
  ! Constants to compute vapor pressure
  real(r8),parameter :: a0=6.107799961_r8, a1=4.436518521e-01_r8, &
       a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
       a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
       a6=6.136820929e-11_r8

  real(r8), parameter :: b0=6.109177956_r8, b1=5.034698970e-01_r8, &
       b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
       b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
       b6=1.838826904e-10_r8

  integer, dimension(13) :: caldaym= (/ 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 /)
  !$acc declare copyin(caldaym)
  !$acc declare copyin(a0,a1,a2,a3,a4,a5,a6)
  !$acc declare copyin(b0,b1,b2,b3,b4,b5,b6)

  public :: update_forcings_CPLBYPASS
  private :: szenith
contains

  real(r8) function tdc(t)
    !$acc routine seq
    implicit none
    real(r8), intent(in) :: t
    tdc = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
  end function tdc
  real(r8) function esatw(t)
    !$acc routine seq
    implicit none
    real(r8), intent(in) :: t
    esatw = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
  end function esatw
  real(r8) function esati(t)
    !$acc routine seq
    implicit none
    real(r8), intent(in) :: t
    esati= 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
  end function esati

  subroutine update_forcings_CPLBYPASS(bounds, atm2lnd_vars, dtime, thiscalday, &
    tod, yr, mon, nstep)
    !$acc routine seq
    implicit none
    type(bounds_type)  , intent(in)   :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    integer , value   , intent(in)    :: dtime
    real(r8), value   , intent(in)    :: thiscalday
    integer , value   , intent(in)    :: tod, yr, mon, nstep
    integer :: v, i, g, av, topo
    integer, parameter :: met_nvars = 7
    integer  :: swrad_period_len, swrad_period_start, thishr, thismin
    integer ::  aindex(2), starti(3), counti(3), tm,nindex(2)
    real(r8) :: wt1(14), wt2(14), tbot, t, qsat
    real(r8) :: swndf, swndr, swvdf, swvdr, ratio_rvrf, frac, q
    real(r8) :: e, ea, vp  ! vapor pressure (Pa)
    real(r8) :: avgcosz, thiscosz
    integer  :: sdate_addt, sy_addt, sm_addt, sd_addt
    integer  :: sdate_addco2, sy_addco2, sm_addco2, sd_addco2
    real(r8) :: forc_rainc    ! rainxy Atm flux mm/s
    real(r8) :: forc_t        ! atmospheric temperature (Kelvin)
    real(r8) :: forc_q        ! atmospheric specific humidity (kg/kg)
    real(r8) :: forc_pbot     ! atmospheric pressure (Pa)
    real(r8) :: forc_rainl    ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc    ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl    ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag ! temporary
    real(r8) :: co2_ppmv_prog ! temporary
    real(r8) :: co2_ppmv_val  ! temporary
    integer  :: co2_type_idx  ! integer flag for co2_type

    associate(  &
      tindex => atm2lnd_vars%tindex  &
      )
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       !!
       do v=1,met_nvars
         if (atm2lnd_vars%npf(v) - 1._r8 .gt. 1e-3) then
           if (v .eq. 4 .or. v .eq. 5 .or. (v .ge. 8 .and. v .le. 13)) then    !rad/Precipitation
             if (mod(tod/int(dtime),nint(atm2lnd_vars%npf(v))) == 1 .and. nstep .gt. 3) then
               tindex(g,v,1) = tindex(g,v,1)+1
               tindex(g,v,2) = tindex(g,v,2)+1
             end if
           else
             if (mod(tod/int(dtime)-1,nint(atm2lnd_vars%npf(v))) <= atm2lnd_vars%npf(v)/2._r8 .and. &
                 mod(tod/int(dtime),nint(atm2lnd_vars%npf(v))) > atm2lnd_vars%npf(v)/2._r8) then
               tindex(g,v,1) = tindex(g,v,1)+1
               tindex(g,v,2) = tindex(g,v,2)+1
             end if
           end if
         else
           tindex(g,v,1) = tindex(g,v,1)+nint(1/atm2lnd_vars%npf(v))
           tindex(g,v,2) = tindex(g,v,2)+nint(1/atm2lnd_vars%npf(v))
         end if

         if (const_climate_hist .or. yr .le. atm2lnd_vars%startyear_met) then
           if (tindex(g,v,1) .gt. atm2lnd_vars%timelen_spinup(v)) tindex(g,v,1) = 1
           if (tindex(g,v,2) .gt. atm2lnd_vars%timelen_spinup(v)) tindex(g,v,2) = 1
         else if (yr .gt. atm2lnd_vars%endyear_met_trans) then
           if (tindex(g,v,1) .gt. atm2lnd_vars%timelen(v)) then
              tindex(g,v,1) = atm2lnd_vars%timelen(v)-atm2lnd_vars%timelen_spinup(v)+1
           end if
           if (tindex(g,v,2) .gt. atm2lnd_vars%timelen(v)) then
              tindex(g,v,2) = atm2lnd_vars%timelen(v)-atm2lnd_vars%timelen_spinup(v)+1
           end if
         end if
       end do
       !!
       !get weights for linear interpolation
       do v=1,met_nvars
         if (atm2lnd_vars%npf(v) - 1._r8 .gt. 1e-3) then
            wt1(v) = 1._r8 - (mod((tod+86400)/dtime-atm2lnd_vars%npf(v)/2._r8, &
                atm2lnd_vars%npf(v))*1._r8)/atm2lnd_vars%npf(v)
                wt2(v) = 1._r8 - wt1(v)
          else
            wt1(v) = 0._r8
            wt2(v) = 1._r8
          end if
        end do

        !Air temperature
        atm2lnd_vars%forc_t_not_downscaled_grc(g)  = min(((atm2lnd_vars%atm_input(1,g,1,tindex(g,1,1))*atm2lnd_vars%scale_factors(1)+ &
                                                  atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,g,1,tindex(g,1,2))* &
                                                  atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                  atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)

        atm2lnd_vars%forc_th_not_downscaled_grc(g) = min(((atm2lnd_vars%atm_input(1,g,1,tindex(g,1,1))*atm2lnd_vars%scale_factors(1)+ &
                                                  atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,g,1,tindex(g,1,2))* &
                                                  atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                  atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)

        tbot = atm2lnd_vars%forc_t_not_downscaled_grc(g)

        !Air pressure
        atm2lnd_vars%forc_pbot_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(2,g,1,tindex(g,2,1))*atm2lnd_vars%scale_factors(2)+ &
                                                    atm2lnd_vars%add_offsets(2))*wt1(2) + (atm2lnd_vars%atm_input(2,g,1,tindex(g,2,2)) &
                                                    *atm2lnd_vars%scale_factors(2)+atm2lnd_vars%add_offsets(2))*wt2(2)) * &
                                                    atm2lnd_vars%var_mult(2,g,mon) + atm2lnd_vars%var_offset(2,g,mon), 4e4_r8)

        !Specific humidity
        atm2lnd_vars%forc_q_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(3,g,1,tindex(g,3,1))*atm2lnd_vars%scale_factors(3)+ &
                                                 atm2lnd_vars%add_offsets(3))*wt1(3) + (atm2lnd_vars%atm_input(3,g,1,tindex(g,3,2)) &
                                                 *atm2lnd_vars%scale_factors(3)+atm2lnd_vars%add_offsets(3))*wt2(3)) * &
                                                 atm2lnd_vars%var_mult(3,g,mon) + atm2lnd_vars%var_offset(3,g,mon), 1e-9_r8)
        if (atm2lnd_vars%metsource == 2) then  !convert RH to qbot
           if (tbot > SHR_CONST_TKFRZ) then
             e = esatw(tdc(tbot))
           else
             e = esati(tdc(tbot))
           end if
           qsat           = 0.622_r8*e / (atm2lnd_vars%forc_pbot_not_downscaled_grc(g) - 0.378_r8*e)
           atm2lnd_vars%forc_q_not_downscaled_grc(g) = qsat * atm2lnd_vars%forc_q_not_downscaled_grc(g) / 100.0_r8
        end if

        !use longwave from file if provided
        atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ((atm2lnd_vars%atm_input(7,g,1,tindex(g,7,1))*atm2lnd_vars%scale_factors(7)+ &
                                                    atm2lnd_vars%add_offsets(7))*wt1(7) + (atm2lnd_vars%atm_input(7,g,1,tindex(g,7,2)) &
                                                    *atm2lnd_vars%scale_factors(7)+atm2lnd_vars%add_offsets(7))*wt2(7)) * &
                                                    atm2lnd_vars%var_mult(7,g,mon) + atm2lnd_vars%var_offset(7,g,mon)
        !=======================================================================================================!!
         if (atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .le. 50 .or. atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .ge. 600) then
         !Longwave radiation (calculated from air temperature, humidity)
            e =  atm2lnd_vars%forc_pbot_not_downscaled_grc(g) * atm2lnd_vars%forc_q_not_downscaled_grc(g) / &
             (0.622_R8 + 0.378_R8 * atm2lnd_vars%forc_q_not_downscaled_grc(g) )
            ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
            atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ea * SHR_CONST_STEBOL * tbot**4
          end if
    !Shortwave radiation (cosine zenith angle interpolation)
    thishr = (tod-dtime/2)/3600
    if (thishr < 0) thishr=thishr+24
    thismin = mod((tod-dtime/2)/60, 60)
    thiscosz = max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr,thismin,0)* &
                    3.14159265358979/180.0d0), 0.001d0)
    avgcosz = 0d0
    if (atm2lnd_vars%npf(4) - 1._r8 .gt. 1e-3) then
      swrad_period_len   = dtime*nint(atm2lnd_vars%npf(4))
      swrad_period_start = ((tod-dtime/2)/swrad_period_len) * swrad_period_len
      !set to last period if first model timestep of the day
      if (tod-dtime/2 < 0) swrad_period_start = ((86400-dtime/2)/swrad_period_len) * swrad_period_len
      do tm=1,nint(atm2lnd_vars%npf(4))
        !Get the average cosine zenith angle over the time resolution of the input data
        thishr  = (swrad_period_start+(tm-1)*dtime+dtime/2)/3600
        if (thishr > 23) thishr=thishr-24
        thismin = mod((swrad_period_start+(tm-1)*dtime+dtime/2)/60, 60)
        avgcosz  = avgcosz + max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr, thismin, 0) &
                   *3.14159265358979/180.0d0), 0.001d0)/atm2lnd_vars%npf(4)
      end do

    else
      avgcosz = thiscosz
    end if
    if (thiscosz > 0.001d0) then
        wt2(4) = min(thiscosz/avgcosz, 10.0_r8)
    else
        wt2(4) = 0d0
    end if

    swndr = max(((atm2lnd_vars%atm_input(4,g,1,tindex(g,4,2))*atm2lnd_vars%scale_factors(4)+ &
                             atm2lnd_vars%add_offsets(4))*wt2(4)) * 0.50_R8, 0.0_r8)

    swndf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(g,4,2))*atm2lnd_vars%scale_factors(4)+ &
                            atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)

    swvdr = max(((atm2lnd_vars%atm_input(4,g,1,tindex(g,4,2))*atm2lnd_vars%scale_factors(4)+ &
                            atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)

    swvdf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(g,4,2))*atm2lnd_vars%scale_factors(4)+ &
                            atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)

    ratio_rvrf =   min(0.99_R8,max(0.29548_R8 + 0.00504_R8*swndr &
                   -1.4957e-05_R8*swndr**2 + 1.4881e-08_R8*swndr**3,0.01_R8))

    atm2lnd_vars%forc_solad_grc(g,2) = ratio_rvrf*swndr
    atm2lnd_vars%forc_solai_grc(g,2) = (1._R8 - ratio_rvrf)*swndf
    ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                       -9.0039e-06_R8*swvdr**2 +8.1351e-09_R8*swvdr**3,0.01_R8))
    atm2lnd_vars%forc_solad_grc(g,1) = ratio_rvrf*swvdr
    atm2lnd_vars%forc_solai_grc(g,1) = (1._R8 - ratio_rvrf)*swvdf

    frac = (atm2lnd_vars%forc_t_not_downscaled_grc(g) - SHR_CONST_TKFRZ)*0.5_R8       ! ramp near freezing
    frac = min(1.0_R8,max(0.0_R8,frac))           ! bound in [0,1]
    
    !Don't interpolate rainfall data
    forc_rainc = 0.1_R8 * frac * max((((atm2lnd_vars%atm_input(5,g,1,tindex(g,5,2))*atm2lnd_vars%scale_factors(5)+ &
                                  atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                  atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
    !!!!
    forc_rainl = 0.9_R8 * frac * max((((atm2lnd_vars%atm_input(5,g,1,tindex(g,5,2))*atm2lnd_vars%scale_factors(5)+ &
                                   atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                   atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
    !!!!
    forc_snowc = 0.1_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,g,1,tindex(g,5,2))*atm2lnd_vars%scale_factors(5)+ &
            atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
    !!!!
    forc_snowl = 0.9_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,g,1,tindex(g,5,2))*atm2lnd_vars%scale_factors(5)+ &
            atm2lnd_vars%add_offsets(5))) * atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
   
   
    !Wind
    atm2lnd_vars%forc_u_grc(g) = (atm2lnd_vars%atm_input(6,g,1,tindex(g,6,1))*atm2lnd_vars%scale_factors(6)+ &
                                 atm2lnd_vars%add_offsets(6))*wt1(6) + (atm2lnd_vars%atm_input(6,g,1,tindex(g,6,2))* &
                                 atm2lnd_vars%scale_factors(6)+atm2lnd_vars%add_offsets(6))*wt2(6)
    if (atm2lnd_vars%metsource == 5) then 
      atm2lnd_vars%forc_v_grc(g) = (atm2lnd_vars%atm_input(14,g,1,tindex(g,14,1))*atm2lnd_vars%scale_factors(14)+ &
                                 atm2lnd_vars%add_offsets(14))*wt1(14) + (atm2lnd_vars%atm_input(14,g,1,tindex(g,14,2))* &
                                 atm2lnd_vars%scale_factors(14)+atm2lnd_vars%add_offsets(14))*wt2(14)
    else
        atm2lnd_vars%forc_v_grc(g) = 0.0_R8 
    end if
    atm2lnd_vars%forc_hgt_grc(g) = 30.0_R8 !(atm2lnd_vars%atm_input(8,g,1,tindex(1))*wt1 + &
                                        !atm2lnd_vars%atm_input(8,g,1,tindex(2))*wt2)    ! zgcmxy  Atm state, default=30m

   
    !!!!
    !------------------------------------Fire data -------------------------------------------------------
    !get weights for interpolation
    wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
    wt2(1) = 1._r8 - wt1(1)
    atm2lnd_vars%forc_hdm(g) = atm2lnd_vars%hdm1(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt1(1) + &
                               atm2lnd_vars%hdm2(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt2(1)
    !------------------------------------------------------!
    atm2lnd_vars%forc_lnfm(g) = atm2lnd_vars%lnfm(g, ((int(thiscalday)-1)*8+tod/(3600*3))+1)

    !DMR note - ndep will NOT be correct if more than 1850 years of model
    !spinup (model year > 1850)
    nindex(1) = min(max(yr-1848,2), 168)
    nindex(2) = min(nindex(1)+1, 168)

    !get weights for interpolation
    wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
    wt2(1) = 1._r8 - wt1(1)

    atm2lnd_vars%forc_ndep_grc(g)    = (atm2lnd_vars%ndep1(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt1(1) + &
                                        atm2lnd_vars%ndep2(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt2(1)) / (365._r8 * 86400._r8)
    !!!!!!!!!!!!!!!!!!
    !------------------------------------Aerosol forcing--------------------------------------------------

    !get weights for interpolation (note this method doesn't get the month boundaries quite right..)
    aindex(1) = mon+1
    if (thiscalday .le. (caldaym(mon+1)+caldaym(mon))/2._r8) then
       wt1(1) = 0.5_r8 + (thiscalday-caldaym(mon))/(caldaym(mon+1)-caldaym(mon))
       aindex(2) = aindex(1)-1
    else
       wt1(1) = 1.0_r8 - (thiscalday-(caldaym(mon+1)+caldaym(mon))/2._r8)/   &
                      (caldaym(mon+1)-caldaym(mon))
       aindex(2) = aindex(1)+1
    end if
    wt2(1) = 1._r8 - wt1(1)

     do av = 1,14
       atm2lnd_vars%forc_aer_grc(g,av)  =  atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
         atm2lnd_vars%ndepind(g,2),aindex(1))*wt1(1)+atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
         atm2lnd_vars%ndepind(g,2),aindex(2))*wt2(1)
     end do

     !set the topounit-level atmospheric state and flux forcings (bypass mode)
     do topo = grc_pp%topi(g), grc_pp%topf(g)
       ! first, all the state forcings
       top_as%tbot(topo)    = atm2lnd_vars%forc_t_not_downscaled_grc(g)      ! forc_txy  Atm state K
       top_as%thbot(topo)   = atm2lnd_vars%forc_th_not_downscaled_grc(g)     ! forc_thxy Atm state K
       top_as%pbot(topo)    = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)   ! ptcmxy    Atm state Pa
       top_as%qbot(topo)    = atm2lnd_vars%forc_q_not_downscaled_grc(g)      ! forc_qxy  Atm state kg/kg
       top_as%ubot(topo)    = atm2lnd_vars%forc_u_grc(g)                     ! forc_uxy  Atm state m/s
       top_as%vbot(topo)    = atm2lnd_vars%forc_v_grc(g)                     ! forc_vxy  Atm state m/s
       top_as%zbot(topo)    = atm2lnd_vars%forc_hgt_grc(g)                   ! zgcmxy    Atm state m


       ! assign the state forcing fields derived from other inputs
       ! Horizontal windspeed (m/s)
       top_as%windbot(topo) = sqrt(top_as%ubot(topo)**2 + top_as%vbot(topo)**2)
       ! Relative humidity (percent)
       if (top_as%tbot(topo) > SHR_CONST_TKFRZ) then
          e = esatw(tdc(top_as%tbot(topo)))
       else
          e = esati(tdc(top_as%tbot(topo)))
       end if
       qsat           = 0.622_r8*e / (top_as%pbot(topo) - 0.378_r8*e)
       top_as%rhbot(topo) = 100.0_r8*(top_as%qbot(topo) / qsat)
       ! partial pressure of oxygen (Pa)
       top_as%po2bot(topo) = o2_molar_const * top_as%pbot(topo)
       ! air density (kg/m**3) - uses a temporary calculation of water vapor pressure (Pa)
       vp = top_as%qbot(topo) * top_as%pbot(topo)  / (0.622_r8 + 0.378_r8 * top_as%qbot(topo))
       top_as%rhobot(topo) = (top_as%pbot(topo) - 0.378_r8 * vp) / (rair * top_as%tbot(topo))

       ! second, all the flux forcings
       top_af%rain(topo)    = forc_rainc + forc_rainl            ! sum of convective and large-scale rain
       top_af%snow(topo)    = forc_snowc + forc_snowl            ! sum of convective and large-scale snow
       top_af%solad(topo,2) = atm2lnd_vars%forc_solad_grc(g,2)   ! forc_sollxy  Atm flux  W/m^2
       top_af%solad(topo,1) = atm2lnd_vars%forc_solad_grc(g,1)   ! forc_solsxy  Atm flux  W/m^2
       top_af%solai(topo,2) = atm2lnd_vars%forc_solai_grc(g,2)   ! forc_solldxy Atm flux  W/m^2
       top_af%solai(topo,1) = atm2lnd_vars%forc_solai_grc(g,1)   ! forc_solsdxy Atm flux  W/m^2
       top_af%lwrad(topo)   = atm2lnd_vars%forc_lwrad_not_downscaled_grc(g)     ! flwdsxy Atm flux  W/m^2
       ! derived flux forcings
       top_af%solar(topo) = top_af%solad(topo,2) + top_af%solad(topo,1) + &
                            top_af%solai(topo,2) + top_af%solai(topo,1)
     end do
    !!This is for co2_type_idx = 0
     co2_ppmv_val = co2_ppmv
     !!!!!========================================= !!!!!
     do topo = grc_pp%topi(g), grc_pp%topf(g)
       top_as%pco2bot(topo) = co2_ppmv_val * 1.e-6_r8 * top_as%pbot(topo)
       if (use_c13) then
          top_as%pc13o2bot(topo) = top_as%pco2bot(topo) * c13ratio;
       end if
     end do
     !
     co2_ppmv_prog = co2_ppmv
     co2_ppmv_diag = co2_ppmv
     ! Determine derived quantities for required fields
     forc_t = atm2lnd_vars%forc_t_not_downscaled_grc(g)
     forc_q = atm2lnd_vars%forc_q_not_downscaled_grc(g)
     forc_pbot = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)

     atm2lnd_vars%forc_hgt_u_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of wind [m]
     atm2lnd_vars%forc_hgt_t_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of temperature [m]
     atm2lnd_vars%forc_hgt_q_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of humidity [m]
     atm2lnd_vars%forc_vp_grc(g)    = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)
     atm2lnd_vars%forc_rho_not_downscaled_grc(g) = &
          (forc_pbot - 0.378_r8 * atm2lnd_vars%forc_vp_grc(g)) / (rair * forc_t)
     atm2lnd_vars%forc_po2_grc(g)   = o2_molar_const * forc_pbot
     atm2lnd_vars%forc_wind_grc(g)  = sqrt(atm2lnd_vars%forc_u_grc(g)**2 + atm2lnd_vars%forc_v_grc(g)**2)
     atm2lnd_vars%forc_solar_grc(g) = atm2lnd_vars%forc_solad_grc(g,1) + atm2lnd_vars%forc_solai_grc(g,1) + &
                                      atm2lnd_vars%forc_solad_grc(g,2) + atm2lnd_vars%forc_solai_grc(g,2)

     atm2lnd_vars%forc_rain_not_downscaled_grc(g)  = forc_rainc + forc_rainl
     atm2lnd_vars%forc_snow_not_downscaled_grc(g)  = forc_snowc + forc_snowl
     if (forc_t > SHR_CONST_TKFRZ) then
        e = esatw(tdc(forc_t))
     else
        e = esati(tdc(forc_t))
     end if
     qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)
     atm2lnd_vars%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
     ! Make sure relative humidity is properly bounded
     ! atm2lnd_vars%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_vars%forc_rh_grc(g) )
     ! atm2lnd_vars%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_vars%forc_rh_grc(g) )

     ! Determine derived quantities for optional fields
     ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
     ! Note that forc_pbot is in Pa

     co2_type_idx = 2
     !get weights/indices for interpolation (assume values represent annual averages)
     nindex(1) = min(max(yr,1850),2100)-1764
     if (thiscalday .le. 182.5) then
       nindex(2) = nindex(1)-1
     else
       nindex(2) = nindex(1)+1
     end if
     wt1(1) = 1._r8 - abs((182.5 - (thiscalday -1._r8))/365._r8)
     wt2(1) = 1._r8 - wt1(1)

     co2_ppmv_val = atm2lnd_vars%co2_input(1,1,nindex(1))*wt1(1) + atm2lnd_vars%co2_input(1,1,nindex(2))*wt2(1)
     !if (startdate_add_co2 .ne. '') then
     !  if ((yr == sy_addco2 .and. mon == sm_addco2 .and. day >= sd_addco2) .or. &
     !      (yr == sy_addco2 .and. mon > sm_addco2) .or. (yr > sy_addco2)) then
     !    co2_ppmv_val=co2_ppmv_val + add_co2
     !  end if
     !end if

     if (use_c13) then
       atm2lnd_vars%forc_pc13o2_grc(g) = (atm2lnd_vars%c13o2_input(1,1,nindex(1))*wt1(1) + &
            atm2lnd_vars%c13o2_input(1,1,nindex(2))*wt2(1)) * 1.e-6_r8 * forc_pbot
     end if

     co2_type_idx = 1
     atm2lnd_vars%forc_pco2_grc(g)   = co2_ppmv_val * 1.e-6_r8 * forc_pbot
     do topo = grc_pp%topi(g), grc_pp%topf(g)
       top_as%pco2bot(topo) = atm2lnd_vars%forc_pco2_grc(g)
       if (use_c13) then
          top_as%pc13o2bot(topo) = atm2lnd_vars%forc_pc13o2_grc(g)
       end if
     end do


   end do

    end associate

  end subroutine update_forcings_CPLBYPASS

  real(r8) function szenith(xcoor, ycoor, ltm, jday, hr, min, offset)
    !Function to calcualte solar zenith angle
    !Used in coupler bypass mode to compute inerpolation for incoming solar
    !$acc routine seq
    use shr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    !inputs
    real(r8) xcoor, ycoor, offset_min
    integer jday, hr, min, ltm, offset
    !working variables
    real(r8) d2r, r2d, lsn, latrad, decrad, decdeg, ha
    real(r8) hangle, harad, saltrad, saltdeg, sazirad, sazideg
    real(r8) szendeg,szenrad
    real,parameter :: pi = 3.14159265358979
    offset_min = offset/60d0   !note assumes 1hr or smaller timestep
    min = min - offset_min

    !adjust time for offsets
    if (min < 0) then
      hr = hr - 1
      min = min+60
    end if
    if (min >= 60) then
      hr = hr+1
      min = min-60
    end if
    if (hr < 0) then
      hr = hr+24
      jday = jday-1
    end if
    if (hr >= 24) then
      hr = hr-24
      jday = jday+1
    end if

    if (jday < 1) jday = 1
    if (xcoor > 180d0) xcoor = xcoor-360d0

    d2r     = pi/180d0
    r2d     = 1/d2r
    lsn     = 12.0d0+((ltm-xcoor)/15.0d0)
    latrad  = ycoor*d2r
    decrad  = 23.45*d2r*sin(d2r*360d0*(284d0+jday)/365d0)
    decdeg  = decrad*r2d
    ha      = hr+min/60.0d0
    hangle  = (lsn-ha)*60.0d0
    harad   = hangle*0.0043633d0

    saltrad = asin((sin(latrad)*sin(decrad))+(cos(latrad)*cos(decrad) &
         *cos(harad)))
    saltdeg = saltrad * r2d
    sazirad = asin(cos(decrad)*sin(harad)/cos(saltrad))
    sazideg = sazirad * r2d

    IF (saltdeg.LT.0.0d0 .OR. saltrad.GT.180.0d0) THEN  ! sun is below horizon
       saltdeg = 0.0d0
       saltrad = 0.0d0
       szendeg = 90.0d0
       szenrad = 90.0d0*d2r
       !mass    = 1229d0**.5d0             ! if solaralt=0 -> sin(0)=0
    ELSE
       szendeg = 90d0-saltdeg
       szenrad = szendeg*d2r
    ENDIF
    szenith = szendeg

  end function szenith

end module ForcingUpdateMod
