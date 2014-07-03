Module DryDepVelocity                                              

  !-----------------------------------------------------------------------  
  !  
  ! Purpose:  
  ! Deposition velocity (m/s) 
  !  
  ! Method:  
  ! This code simulates dry deposition velocities using the Wesely scheme.   
  ! Details of this method can be found in:  
  ! 
  ! M.L Wesely. Parameterization of surface resistances to gaseous dry deposition 
  ! in regional-scale numericl models. 1989. Atmospheric Environment vol.23 No.6  
  ! pp. 1293-1304. 
  ! 
  ! In Wesely (1998) "the magnitude of the dry deposition velocity can be found 
  ! as: 
  ! 
  !  |vd|=(ra+rb+rc)^-1 
  ! 
  ! where ra is the aerodynamic resistance (common to all gases) between a  
  ! specific height and the surface, rb is the quasilaminar sublayer resistance 
  ! (whose only dependence on the porperties of the gas of interest is its 
  ! molecular diffusivity in air), and rc is the bulk surface resistance". 
  ! 
  ! In this subroutine both ra and rb are calculated elsewhere in CLM.  Thus ra  
  ! and rb were "globalized" in order to gain access to them for the calculation. 
  ! "ram1" is the CLM variable used for ra.  ram1 was globalized in the following 
  ! subroutines; BareGroundFluxes.F90, Biogeophysics_lake.F90, CanopyFluxes.F90, 
  ! and clmtype.F90.  
  ! 
  ! "rb" is the CLM variable used for rb in the Wesely equation above.  rb was  
  ! globalized in the following subroutines; clmtype.F90     
  ! 
  ! In Wesely (1989) rc is estimated for five seasonal categories and 11 landuse 
  ! types.  For each season and landuse type, Wesely compiled data into a  
  ! look-up-table for several parameters used to calculate rc. In this subroutine 
  ! the same values are used as found in wesely's look-up-tables, the only  
  ! difference is that this subroutine uses a CLM generated LAI to select values 
  ! from the look-up-table instead of seasonality.  Inaddition, Wesely(1989)  
  ! land use types are "mapped" into CLM plant function types (PFT). 
  ! 
  ! Subroutine written to operate at the patch level. 
  ! 
  ! Output: 
  ! 
  ! vd(n_species) !Dry deposition velocity [m s-1] for each molecule or species 
  !  
  ! Author: Beth Holland and  James Sulzman 
  ! 
  ! Modified: Francis Vitt -- 30 Mar 2007
  !----------------------------------------------------------------------- 

  use shr_kind_mod,         only : r8 => shr_kind_r8 
  use clmtype 
  use abortutils,           only : endrun
  use clm_time_manager,     only : get_nstep, get_curr_date, get_curr_time 
  use clm_atmlnd,           only : clm_a2l
  use spmdMod,              only : masterproc
  use seq_drydep_mod,       only : n_drydep, drydep_list
  use seq_drydep_mod,       only : drydep_method, DD_XLND
  use seq_drydep_mod,       only : index_o3=>o3_ndx, index_o3a=>o3a_ndx, index_so2=>so2_ndx, index_h2=>h2_ndx, &
                                   index_co=>co_ndx, index_ch4=>ch4_ndx, index_pan=>pan_ndx, &
                                   index_xpan=>xpan_ndx
  implicit none 
  save 

  private

  public :: depvel_compute

CONTAINS 

  !----------------------------------------------------------------------- 
  ! computes the dry deposition velocity of tracers
  !----------------------------------------------------------------------- 
  subroutine depvel_compute( lbp , ubp ) 
    use shr_const_mod     , only :  tmelt => shr_const_tkfrz
    use seq_drydep_mod    , only :  seq_drydep_setHCoeff, mapping, drat, foxd, &
                                    rcls, h2_a, h2_b, h2_c, ri, rac, rclo, rlu, &
                                    rgss, rgso
    use clm_varcon        , only : istsoil, istice, istice_mec, istslak, istdlak, istwet, isturb
    use clm_varctl        , only : iulog
    use pftvarcon         , only : noveg, ndllf_evr_tmp_tree, ndllf_evr_brl_tree,   &
                                   ndllf_dcd_brl_tree,        nbrdlf_evr_trp_tree,  &
                                   nbrdlf_evr_tmp_tree,       nbrdlf_dcd_trp_tree,  &
                                   nbrdlf_dcd_tmp_tree,       nbrdlf_dcd_brl_tree,  &
                                   nbrdlf_evr_shrub,          nbrdlf_dcd_tmp_shrub, &
                                   nbrdlf_dcd_brl_shrub,      nc3_arctic_grass,     &
                                   nc3_nonarctic_grass,       nc4_grass, nc3crop,   &
                                   nc3irrig,       npcropmin, npcropmax

    implicit none 

    !----Arguments-----------------------------------------------------

    integer, intent(in) :: lbp, ubp                    ! pft bounds

    ! ------------------------ local variables ------------------------ 
    ! local pointers to implicit in arguments 
    logical , pointer :: pactive(:)       ! true=>do computations on this pft (see reweightMod for details)
    integer , pointer :: plandunit(:)     !pft's landunit index
    integer , pointer :: ivt(:)           !landunit type
    integer , pointer :: pgridcell(:)     !pft's gridcell index
    real(r8), pointer :: elai(:)          !one-sided leaf area index with burying by snow 
    real(r8), pointer :: forc_t(:)        !atmospheric temperature (Kelvin) 
    real(r8), pointer :: forc_q(:)        !atmospheric specific humidity (kg/kg) 
    real(r8), pointer :: forc_psrf(:)     !surface pressure (Pa) 
    real(r8), pointer :: latdeg(:)        !latitude (degrees) 
    real(r8), pointer :: londeg(:)        !longitude (degrees) 
    real(r8), pointer :: forc_rain(:)     !rain rate [mm/s] 
    real(r8), pointer :: forc_solad(:,:)  !direct beam radiation (visible only) 
    real(r8), pointer :: forc_solai(:,:)  !direct beam radiation (visible only) 
    real(r8), pointer :: ram1(:)          !aerodynamical resistance 
    real(r8), pointer :: vds(:)           !aerodynamical resistance 
    real(r8), pointer :: rssun(:)         !stomatal resistance 
    real(r8), pointer :: rssha(:)         !shaded stomatal resistance (s/m) 
    real(r8), pointer :: fsun(:)          !sunlit fraction of canopy 
    real(r8), pointer :: rb1(:)           !leaf boundary layer resistance [s/m]
    real(r8), pointer :: annlai(:,:)      !12 months of monthly lai from input data set 
    real(r8), pointer :: mlaidiff(:)      !difference in lai between month one and month two 
    real(r8), pointer :: velocity(:,:)
    real(r8), pointer :: snow_depth(:)        ! snow height (m)

    integer, pointer :: pcolumn(:)        ! column index associated with each pft
    integer :: c
    integer , pointer :: itypelun(:) 	   ! landunit type

    real(r8), pointer :: h2osoi_vol(:,:)    ! volumetric soil water (0<=h2osoi_vol<=watsat)
    real(r8) :: soilw, var_soilw, fact_h2, dv_soil_h2

    ! new local variables  
    integer :: pi,g, l
    integer :: ispec 
    integer :: length 
    integer :: wesveg              !wesely vegegation index  
    integer :: clmveg              !clm veg index from ivegtype 
    integer :: i 
    integer :: index_season        !seasonal index based on LAI.  This indexs wesely data tables 
    integer :: nstep               !current step 
    integer :: indexp 

    real(r8) :: pg          ! surface pressure 
    real(r8) :: tc          ! temperature in celsius  
    real(r8) :: rs          ! constant for calculating rsmx  
    real(r8) :: es          ! saturation vapor pressur 
    real(r8) :: ws          ! saturation mixing ratio 
    real(r8) :: rmx         ! resistance by vegetation 
    real(r8) :: qs          ! saturation specific humidity 
    real(r8) :: dewm        ! multiplier for rs when dew occurs 
    real(r8) :: crs         ! multiplier to calculate crs 
    real(r8) :: rdc         ! part of lower canopy resistance 
    real(r8) :: rain        ! rain fall 
    real(r8) :: spec_hum    ! specific humidity 
    real(r8) :: solar_flux  ! solar radiation(direct beam) W/m2 
    real(r8) :: lat         ! latitude in degrees 
    real(r8) :: lon         ! longitude in degrees 
    real(r8) :: sfc_temp    ! surface temp 
    real(r8) :: minlai      ! minimum of monthly lai 
    real(r8) :: maxlai      ! maximum of monthly lai 
    real(r8) :: rds

    logical  :: has_dew
    logical  :: has_rain
    real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

    ! local arrays: dependent on species only 
    ! 

    real(r8), dimension(n_drydep) :: rsmx !vegetative resistance (plant mesophyll) 
    real(r8), dimension(n_drydep) :: rclx  !lower canopy resistance  
    real(r8), dimension(n_drydep) :: rlux  !vegetative resistance (upper canopy) 
    real(r8), dimension(n_drydep) :: rgsx  !gournd resistance 
    real(r8), dimension(n_drydep) :: heff   
    real(r8) :: rc    !combined surface resistance 
    real(r8) :: cts   !correction to flu rcl and rgs for frost 
    real(r8) :: rlux_o3

    ! constants 
    real(r8), parameter :: slope = 0._r8      ! Used to calculate  rdc in (lower canopy resistance) 
    integer, parameter :: wveg_unset = -1     ! Unset Wesley vegetation type

    character(len=32), parameter :: subname = "depvel_compute"

    !-------------------------------------------------------------------------------------
    ! jfl : mods for PAN
    !-------------------------------------------------------------------------------------
    real(r8)                  :: dv_pan
    real(r8) :: c0_pan(11) = (/ 0.000_r8, 0.006_r8, 0.002_r8, 0.009_r8, 0.015_r8, &
                                0.006_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.002_r8, 0.002_r8 /)
    real(r8) :: k_pan (11) = (/ 0.000_r8, 0.010_r8, 0.005_r8, 0.004_r8, 0.003_r8, &
                                0.005_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.075_r8, 0.002_r8 /)
    !----------------------------------------------------------------------- 
    if ( n_drydep == 0 .or. drydep_method /= DD_XLND ) return

    ! local pointers to original implicit out arrays 

    ! Assign local pointers to derived subtypes components (column-level) 
    forc_t     => clm_a2l%forc_t
    forc_q     => clm_a2l%forc_q
    forc_psrf  => clm_a2l%forc_pbot
    forc_rain  => clm_a2l%forc_rain 

    latdeg     => clm3%g%latdeg
    londeg     => clm3%g%londeg
    pactive    => clm3%g%l%c%p%active
    ivt        => clm3%g%l%c%p%itype
    elai       => clm3%g%l%c%p%pps%elai 
    ram1       => clm3%g%l%c%p%pps%ram1 
    vds        => clm3%g%l%c%p%pps%vds
    fsun       => clm3%g%l%c%p%pps%fsun 
    rssun      => clm3%g%l%c%p%pps%rssun 
    rssha      => clm3%g%l%c%p%pps%rssha 
    rb1        => clm3%g%l%c%p%pps%rb1 
    mlaidiff   => clm3%g%l%c%p%pps%mlaidiff 
    annlai     => clm3%g%l%c%p%pps%annlai    

    forc_solai => clm_a2l%forc_solai 
    forc_solad => clm_a2l%forc_solad

    pgridcell  => clm3%g%l%c%p%gridcell
    plandunit  => clm3%g%l%c%p%landunit

    pcolumn    => clm3%g%l%c%p%column
    itypelun   => clm3%g%l%itype

    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol

    velocity   => clm3%g%l%c%p%pdd%drydepvel ! cm/sec

    snow_depth        => clm3%g%l%c%cps%snow_depth

    ! Assign local pointers to original implicit out arrays 
    !_________________________________________________________________ 
    ! 
    ! Begin loop through pfts 
    pft_loop: do pi = lbp,ubp
       l = plandunit(pi)

       active: if (pactive(pi)) then

          c = pcolumn(pi)
          g = pgridcell(pi)
          !solar_flux = forc_lwrad  !rename CLM variables to fit with Dry Dep variables 

          pg         = forc_psrf(g)  
          spec_hum   = forc_q(g)
          rain       = forc_rain(g) 
          sfc_temp   = forc_t(g) 
          lat        = latdeg(g) 
          lon        = londeg(g) 
          solar_flux = forc_solad(g,1) 
          clmveg     = ivt(pi) 
          soilw = h2osoi_vol(c,1)

          !  print *,'bb',pi,cps%npfts,lat,lon,clmveg 
          !map CLM veg type into Wesely veg type  
          wesveg = wveg_unset 
          if (clmveg == noveg                               ) wesveg = 8 
          if (clmveg == ndllf_evr_tmp_tree                  ) wesveg = 5 
          if (clmveg == ndllf_evr_brl_tree                  ) wesveg = 5 
          if (clmveg == ndllf_dcd_brl_tree                  ) wesveg = 5 
          if (clmveg == nbrdlf_evr_trp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_evr_tmp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_dcd_trp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_dcd_tmp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_dcd_brl_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_evr_shrub                    ) wesveg = 11 
          if (clmveg == nbrdlf_dcd_tmp_shrub                ) wesveg = 11 
          if (clmveg == nbrdlf_dcd_brl_shrub                ) wesveg = 11 
          if (clmveg == nc3_arctic_grass                    ) wesveg = 3 
          if (clmveg == nc3_nonarctic_grass                 ) wesveg = 3 
          if (clmveg == nc4_grass                           ) wesveg = 3 
          if (clmveg == nc3crop                             ) wesveg = 2 
          if (clmveg == nc3irrig                            ) wesveg = 2 
          if (clmveg >= npcropmin .and. clmveg <= npcropmax ) wesveg = 2 
          if (wesveg == wveg_unset )then
             write(iulog,*) 'clmveg = ', clmveg, 'itypelun = ', itypelun(l)
             call endrun( subname//': Not able to determine Wesley vegetation type')
          end if

          ! creat seasonality index used to index wesely data tables from LAI,  Bascially 
          !if elai is between max lai from input data and half that max the index_season=1 


          !mail1j and mlai2j are the two monthly lai values pulled from a CLM input data set 
          !/fs/cgd/csm/inputdata/lnd/clm2/rawdata/mksrf_lai.nc.  lai for dates in the middle  
          !of the month are interpolated using using these values and stored in the variable  
          !elai (done elsewhere).  If the difference between mlai1j and mlai2j is greater 
          !than zero it is assumed to be fall and less than zero it is assumed to be spring. 

          !wesely seasonal "index_season" 
          ! 1 - midsummer with lush vegetation 
          ! 2 - Autumn with unharvested cropland 
          ! 3 - Late autumn after frost, no snow 
          ! 4 - Winter, snow on ground and subfreezing 
          ! 5 - Transitional spring with partially green short annuals 


          !mlaidiff=jan-feb 
          minlai=minval(annlai(:,pi)) 
          maxlai=maxval(annlai(:,pi)) 

          index_season = -1

          if ( itypelun(l) /= istsoil )then
             if ( itypelun(l) == istice .or. itypelun(l) == istice_mec ) then
                wesveg       = 8
                index_season = 4
             elseif ( itypelun(l) == istdlak .or. itypelun(l) == istslak ) then
                wesveg       = 7
                index_season = 4
             elseif ( itypelun(l) == istwet ) then
                wesveg       = 9
                index_season = 2
             elseif ( itypelun(l) == isturb ) then
                wesveg       = 1
                index_season = 2
             end if
          else if ( snow_depth(c) > 0 ) then
             index_season = 4
          else if(elai(pi).gt.0.5_r8*maxlai) then  
             index_season = 1  
          endif

          if (index_season<0) then 
             if (elai(pi).lt.(minlai+0.05*(maxlai-minlai))) then  
                index_season = 3 
             endif
          endif

          if (index_season<0) then 
             if (mlaidiff(pi).gt.0.0_r8) then 
                index_season = 2 
             elseif (mlaidiff(pi).lt.0.0_r8) then 
                index_season = 5 
             elseif (mlaidiff(pi).eq.0.0_r8) then 
                index_season = 3 
             endif
          endif

          if (index_season<0) then 
             call endrun( subname//': not able to determine season')
          endif

          ! saturation specific humidity 
          ! 
          es = 611_r8*exp(5414.77_r8*((1._r8/tmelt)-(1._r8/sfc_temp))) 
          ws = .622_r8*es/(pg-es) 
          qs = ws/(1._r8+ws) 

          has_dew = .false.
          if( qs <= spec_hum ) then
             has_dew = .true.
          end if
          if( sfc_temp < tmelt ) then
             has_dew = .false.
          end if

          has_rain = rain > rain_threshold

          if ( has_dew .or. has_rain ) then
             dewm = 3._r8
          else
             dewm = 1._r8
          end if

          ! 
          ! constant in determining rs 
          ! 
          crs = 1.e36_r8

          tc = sfc_temp - tmelt 
          if(sfc_temp.gt.tmelt.and.sfc_temp.lt.313.15_r8) then 
             crs = (1._r8+(200._r8/(solar_flux+.1_r8))**2) * (400._r8/(tc*(40._r8-tc))) 
          endif
          ! 
          ! rdc (lower canopy res) 
          ! 
          rdc=100._r8*(1._r8+1000._r8/(solar_flux+10._r8))/(1._r8+1000._r8*slope) 

          ! surface resistance : depends on both land type and species 
          ! land types are computed seperately, then resistance is computed as average of values 
          ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1 
          ! 
          ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet 
          ! 

          !******************************************************* 
          call seq_drydep_setHCoeff( sfc_temp, heff(:n_drydep) )
          !********************************************************* 

          species_loop1: do ispec=1, n_drydep
             if(mapping(ispec).le.0) cycle 

             if(ispec.eq.index_o3.or.ispec.eq.index_o3a.or.ispec.eq.index_so2) then 
                rmx=0._r8
             else 
                rmx=1._r8/((heff(ispec)/3000._r8)+(100._r8*foxd(ispec))) 
             endif

             ! correction for frost 
             cts = 1000._r8*exp( -tc - 4._r8 )                 ! correction for frost
             rgsx(ispec) = cts + 1._r8/((heff(ispec)/(1.e5_r8*rgss(index_season,wesveg))) + & 
                                        (foxd(ispec)/rgso(index_season,wesveg))) 

             !-------------------------------------------------------------------------------------
             ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
             !-------------------------------------------------------------------------------------
             if( ispec == index_h2 .or. ispec == index_co .or. ispec == index_ch4 ) then

                if( ispec == index_co ) then
                   fact_h2 = 1.0_r8
                elseif ( ispec == index_h2 ) then
                   fact_h2 = 0.5_r8
                elseif ( ispec == index_ch4 ) then
                   fact_h2 = 50.0_r8
                end if
                !-------------------------------------------------------------------------------------
                ! no deposition on snow, ice, desert, and water
                !-------------------------------------------------------------------------------------
                if( wesveg == 1 .or. wesveg == 7 .or. wesveg == 8 .or. index_season == 4 ) then
                   rgsx(ispec) = 1.e36_r8
                else
                   var_soilw = max( .1_r8,min( soilw,.3_r8 ) )
                   if( wesveg == 3 ) then
                      var_soilw = log( var_soilw )
                   end if
                   dv_soil_h2 = h2_c(wesveg) + var_soilw*(h2_b(wesveg) + var_soilw*h2_a(wesveg))
                   if( dv_soil_h2 > 0._r8 ) then
                      rgsx(ispec) = fact_h2/(dv_soil_h2*1.e-4_r8)
                   end if
                end if
             end if

             rs=ri(index_season,wesveg)*crs 
             if(wesveg.eq.7) then ! over water
                rclx(ispec)=1.e36_r8
                rsmx(ispec)=1.e36_r8
                rlux(ispec)=1.e36_r8
             else 
                ! ??? fvitt   rs=(fsun(pi)*rssun(pi))+(rssha(pi)*(1.-fsun(pi))) -- made the same as mo_drydep
                rsmx(ispec) = (dewm*rs*drat(ispec)+rmx) 
             endif

             !-------------------------------------------------------------------------------------
             ! jfl : special case for PAN
             !-------------------------------------------------------------------------------------
             if( ispec == index_pan .or. ispec == index_xpan ) then
                dv_pan =  c0_pan(wesveg) * (1._r8 - exp( -k_pan(wesveg)*(dewm*rs*drat(ispec))*1.e-2_r8 ))
                if( dv_pan > 0._r8 .and. index_season /= 4 ) then
                   rsmx(ispec) = ( 1._r8/dv_pan )
                end if
             end if

             rclx(ispec) = cts + 1._r8/((heff(ispec)/(1.e5_r8*rcls(index_season,wesveg))) + & 
                           (foxd(ispec)/rclo(index_season,wesveg))) 
             rlux(ispec) = cts + rlu(index_season,wesveg)/(1.e-5_r8*heff(ispec)+foxd(ispec)) 

          end do species_loop1
          
          ! 
          ! no effect over water
          ! 
          no_water: if(wesveg.ne.7) then 
             ! 
             ! no effect if sfc_temp < O C 
             ! 
             non_freezing: if(sfc_temp.gt.tmelt) then

                if( has_dew ) then 
                   rlux_o3 = 1._r8/((1._r8/3000._r8)+(1._r8/(3._r8*rlu(index_season,wesveg)))) 
                   if (index_o3 > 0) then
                      rlux(index_o3) = rlux_o3
                   endif
                   if (index_o3a > 0) then
                      rlux(index_o3a) = rlux_o3
                   endif
                endif

                if(has_rain) then 
                   rlux_o3 = 1._r8/((1._r8/1000._r8)+(1._r8/(3._r8*rlu(index_season,wesveg)))) 
                   if (index_o3 > 0) then
                      rlux(index_o3) = rlux_o3
                   endif
                   if (index_o3a > 0) then
                      rlux(index_o3a) = rlux_o3
                   endif
                endif

                if ( index_o3 > 0 ) then
                   rclx(index_o3) = cts + rclo(index_season,wesveg)
                   rlux(index_o3) = cts + rlux(index_o3)
                end if
                if ( index_o3a > 0 ) then
                   rclx(index_o3a) = cts + rclo(index_season,wesveg)
                   rlux(index_o3a) = cts + rlux(index_o3a)
                end if

                species_loop2: do ispec=1,n_drydep 
                   if(mapping(ispec).le.0) cycle 
                   if(ispec.ne.index_o3.and.ispec.ne.index_o3a.and.ispec.ne.index_so2) then 

                      if( has_dew ) then
                         rlux(ispec)=1._r8/((1._r8/(3._r8*rlux(ispec)))+ & 
                              (1.e-7_r8*heff(ispec))+(foxd(ispec)/rlux_o3)) 
                      endif

                   elseif(ispec.eq.index_so2) then 

                      if( has_dew ) then
                         rlux(ispec) = 100._r8
                      endif

                      if(has_rain) then 
                         rlux(ispec) = 1._r8/((1._r8/5000._r8)+(1._r8/(3._r8*rlu(index_season,wesveg)))) 
                      endif

                      rclx(ispec) = cts + rcls(index_season,wesveg)
                      rlux(ispec) = cts + rlux(ispec)

                      if( has_dew .or. has_rain ) then
                         rlux(ispec)=50._r8
                      endif

                   endif
                end do species_loop2

             endif non_freezing

          endif no_water

          rds = 1._r8/vds(pi)

          species_loop3: do ispec=1,n_drydep 
             if(mapping(ispec).le.0) cycle 

             ! 
             ! compute rc 
             ! 
             rc = 1._r8/((1._r8/rsmx(ispec))+(1._r8/rlux(ispec)) + & 
                        (1._r8/(rdc+rclx(ispec)))+(1._r8/(rac(index_season,wesveg)+rgsx(ispec)))) 
             rc = max( 10._r8, rc)
             !
             ! assume no surface resistance for SO2 over water
             !
             if ( drydep_list(ispec) == 'SO2' .and. wesveg == 7 ) then
               rc = 0._r8
             end if

             select case( drydep_list(ispec) )
             case ( 'SO4' )
                velocity(pi,ispec) = (1._r8/(ram1(pi)+rds))*100._r8
             case ( 'NH4','NH4NO3','XNH4NO3' )
                velocity(pi,ispec) = (1._r8/(ram1(pi)+0.5_r8*rds))*100._r8
             case ( 'Pb' )
                velocity(pi,ispec) = 0.2_r8
             case ( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB', 'SOAX' )
                velocity(pi,ispec) = 0.10_r8
             case default
                velocity(pi,ispec) = (1._r8/(ram1(pi)+rb1(pi)+rc))*100._r8
             end select
          end do species_loop3

       endif active

    end do pft_loop

  end subroutine depvel_compute

end module DryDepVelocity                    
