module SoilTemperatureMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SoilTemperatureMod
!
! !DESCRIPTION:
! Calculates snow and soil temperatures including phase change
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilTemperature  ! Snow and soil temperatures including phase change
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilThermProp   ! Set therm conductivities and heat cap of snow/soil layers
  private :: PhaseChange     ! Calculation of the phase change within snow and soil layers
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilTemperature
!
! !INTERFACE:
  subroutine SoilTemperature(lbl, ubl, lbc, ubc, num_urbanl, filter_urbanl, &
                             num_nolakec, filter_nolakec, xmf, fact)
!
! !DESCRIPTION:
! Snow and soil temperatures including phase change
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), and the
!   conductivity of snow is from the formulation used in
!   SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction
!   in 10 soil layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two
!   neighboring layers (j, j+1) are derived from an assumption that
!   the flux across the interface is equal to that from the node j
!   to the interface and the flux from the interface to the node j+1.
!   The equation is solved using the Crank-Nicholson method and
!   results in a tridiagonal system equation.
!
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd    , only : clm_a2l
    use clm_time_manager  , only : get_step_size
    use clm_varctl    , only : iulog
    use clm_varcon    , only : sb, capr, cnfac, hvap, istice_mec, isturb, &
                               icol_roof, icol_sunwall, icol_shadewall, &
                               icol_road_perv, icol_road_imperv, istwet
    use clm_varpar    , only : nlevsno, nlevgrnd, max_pft_per_col, nlevurb
    use TridiagonalMod, only : Tridiagonal

!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                    ! column bounds
    integer , intent(in)  :: num_nolakec                 ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer , intent(in)  :: lbl, ubl                    ! landunit-index bounds
    integer , intent(in)  :: num_urbanl                  ! number of urban landunits in clump
    integer , intent(in)  :: filter_urbanl(ubl-lbl+1)    ! urban landunit filter
    real(r8), intent(out) :: xmf(lbc:ubc)                ! total latent heat of phase change of ground water
    real(r8), intent(out) :: fact(lbc:ubc, -nlevsno+1:nlevgrnd) ! used in computing tridiagonal matrix
!
! !CALLED FROM:
! subroutine Biogeophysics2 in module Biogeophysics2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! Changed references for tg to t_grnd, for consistency with the
! rest of the code (tg eliminated as redundant)
! 2/14/02, Peter Thornton: Migrated to new data structures. Added pft loop
! in calculation of net ground heat flux.
! 3/18/08, David Lawrence: Change nlevsoi to nlevgrnd for deep soil
! 03/28/08, Mark Flanner: Changes to allow solar radiative absorption in all snow layers and top soil layer
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: plandunit(:)       ! pft's landunit index
    integer , pointer :: clandunit(:)       ! column's landunit
    integer , pointer :: ltype(:)           ! landunit type
    integer , pointer :: ctype(:)           ! column type
    integer , pointer :: npfts(:)           ! column's number of pfts 
    integer , pointer :: pfti(:)            ! column's beginning pft index 
    real(r8), pointer :: pwtcol(:)          ! weight of pft relative to column
    real(r8), pointer :: pwtgcell(:)        ! weight of pft relative to corresponding gridcell
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    integer , pointer :: snl(:)             ! number of snow layers
    real(r8), pointer :: htvp(:)            ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: emg(:)             ! ground emissivity
    real(r8), pointer :: cgrnd(:)           ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: dlrad(:)           ! downward longwave radiation blow the canopy [W/m2]
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    integer , pointer :: frac_veg_nosno(:)  ! fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
    real(r8), pointer :: eflx_sh_grnd(:)    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:)   ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: zi(:,:)            ! interface level below a "z" level (m)
    real(r8), pointer :: dz(:,:)            ! layer depth (m)
    real(r8), pointer :: z(:,:)             ! layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)      ! soil temperature (Kelvin)
    real(r8), pointer :: eflx_lwrad_net(:)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: tssbef(:,:)        ! temperature at previous time step [K]
    real(r8), pointer :: t_building(:)      ! internal building temperature (K)
    real(r8), pointer :: t_building_max(:)  ! maximum internal building temperature (K)
    real(r8), pointer :: t_building_min(:)  ! minimum internal building temperature (K)
    real(r8), pointer :: hc_soi(:)          ! soil heat content (MJ/m2)
    real(r8), pointer :: hc_soisno(:)       ! soil plus snow plus lake heat content (MJ/m2)
    real(r8), pointer :: eflx_fgr12(:)      ! heat flux between soil layer 1 and 2 (W/m2)
    real(r8), pointer :: eflx_traffic(:)    ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_wasteheat(:)  ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_wasteheat_pft(:)  ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac(:) !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac_pft(:) !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: eflx_traffic_pft(:)    ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_anthro(:)         ! total anthropogenic heat flux (W/m**2)
    real(r8), pointer :: canyon_hwr(:)      ! urban canyon height to width ratio
    real(r8), pointer :: wtlunit_roof(:)    ! weight of roof with respect to landunit
    real(r8), pointer :: eflx_bot(:)        ! heat flux from beneath column (W/m**2) [+ = upward]
! 
! local pointers to  original implicit inout arguments
!
    real(r8), pointer :: t_grnd(:)          ! ground surface temperature [K]
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: eflx_gnet(:)          ! net ground heat flux into the surface (W/m**2)
    real(r8), pointer :: dgnetdT(:)            ! temperature derivative of ground net heat flux
    real(r8), pointer :: eflx_building_heat(:) ! heat flux from urban building interior to walls, roof (W/m**2)

! variables needed for SNICAR
    real(r8), pointer :: sabg_lyr(:,:)      ! absorbed solar radiation (pft,lyr) [W/m2]
    real(r8), pointer :: h2osno(:)          ! total snow water (col) [kg/m2]
    real(r8), pointer :: h2osoi_liq(:,:)    ! liquid water (col,lyr) [kg/m2]
    real(r8), pointer :: h2osoi_ice(:,:)    ! ice content (col,lyr) [kg/m2]

! Urban building HAC fluxes
    real(r8), pointer :: eflx_urban_ac(:)      ! urban air conditioning flux (W/m**2)
    real(r8), pointer :: eflx_urban_heat(:)    ! urban heating flux (W/m**2)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: j,c,p,l,g,pi                       !  indices
    integer  :: fc                                 ! lake filtered column indices
    integer  :: fl                                 ! urban filtered landunit indices
    integer  :: jtop(lbc:ubc)                      ! top level at each column
    real(r8) :: dtime                              ! land model time step (sec)
    real(r8) :: at (lbc:ubc,-nlevsno+1:nlevgrnd)   ! "a" vector for tridiagonal matrix
    real(r8) :: bt (lbc:ubc,-nlevsno+1:nlevgrnd)   ! "b" vector for tridiagonal matrix
    real(r8) :: ct (lbc:ubc,-nlevsno+1:nlevgrnd)   ! "c" vector for tridiagonal matrix
    real(r8) :: rt (lbc:ubc,-nlevsno+1:nlevgrnd)   ! "r" vector for tridiagonal solution
    real(r8) :: cv (lbc:ubc,-nlevsno+1:nlevgrnd)   ! heat capacity [J/(m2 K)]
    real(r8) :: tk (lbc:ubc,-nlevsno+1:nlevgrnd)   ! thermal conductivity [W/(m K)]
    real(r8) :: fn (lbc:ubc,-nlevsno+1:nlevgrnd)   ! heat diffusion through the layer interface [W/m2]
    real(r8) :: fn1(lbc:ubc,-nlevsno+1:nlevgrnd)   ! heat diffusion through the layer interface [W/m2]
    real(r8) :: brr(lbc:ubc,-nlevsno+1:nlevgrnd)   ! temporary
    real(r8) :: dzm                                ! used in computing tridiagonal matrix
    real(r8) :: dzp                                ! used in computing tridiagonal matrix
    real(r8) :: hs(lbc:ubc)                        ! net energy flux into the surface (w/m2)
    real(r8) :: dhsdT(lbc:ubc)                     ! d(hs)/dT
    real(r8) :: lwrad_emit(lbc:ubc)                ! emitted longwave radiation
    real(r8) :: dlwrad_emit(lbc:ubc)               ! time derivative of emitted longwave radiation
    integer  :: lyr_top                            ! index of top layer of snowpack (-4 to 0) [idx]
    real(r8) :: sabg_lyr_col(lbc:ubc,-nlevsno+1:1) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) :: eflx_gnet_top                      ! net energy flux into surface layer, pft-level [W/m2]
    real(r8) :: hs_top(lbc:ubc)                    ! net energy flux into surface layer (col) [W/m2]
    logical  :: cool_on(lbl:ubl)                   ! is urban air conditioning on?
    logical  :: heat_on(lbl:ubl)                   ! is urban heating on?
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_lwrad     => clm_a2l%forc_lwrad

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype          => lun%itype
    t_building     => lps%t_building
    t_building_max => lps%t_building_max
    t_building_min => lps%t_building_min
    eflx_traffic   => lef%eflx_traffic
    canyon_hwr     => lun%canyon_hwr
    eflx_wasteheat => lef%eflx_wasteheat
    eflx_heat_from_ac => lef%eflx_heat_from_ac
    wtlunit_roof   => lun%wtlunit_roof

    ! Assign local pointers to derived subtypes components (column-level)

    ctype          => col%itype
    clandunit      => col%landunit
    npfts          => col%npfts
    pfti           => col%pfti
    snl            => cps%snl
    htvp           => cps%htvp
    emg            => cps%emg
    t_grnd         => ces%t_grnd
    hc_soi         => ces%hc_soi
    hc_soisno      => ces%hc_soisno
    eflx_fgr12     => cef%eflx_fgr12
    zi             => cps%zi
    dz             => cps%dz
    z              => cps%z
    t_soisno       => ces%t_soisno
    eflx_building_heat => cef%eflx_building_heat
    tssbef             => ces%tssbef
    eflx_urban_ac      => cef%eflx_urban_ac
    eflx_urban_heat    => cef%eflx_urban_heat
    eflx_bot           => cef%eflx_bot

    ! Assign local pointers to derived subtypes components (pft-level)

    pgridcell      => pft%gridcell
    plandunit      => pft%landunit
    pwtcol         => pft%wtcol
    pwtgcell       => pft%wtgcell  
    frac_veg_nosno => pps%frac_veg_nosno
    cgrnd          => pef%cgrnd
    dlrad          => pef%dlrad
    sabg           => pef%sabg
    eflx_sh_grnd   => pef%eflx_sh_grnd
    qflx_evap_soi  => pwf%qflx_evap_soi
    qflx_tran_veg  => pwf%qflx_tran_veg
    eflx_gnet      => pef%eflx_gnet
    dgnetdT        => pef%dgnetdT
    eflx_lwrad_net => pef%eflx_lwrad_net
    eflx_wasteheat_pft => pef%eflx_wasteheat_pft
    eflx_heat_from_ac_pft => pef%eflx_heat_from_ac_pft
    eflx_traffic_pft => pef%eflx_traffic_pft
    eflx_anthro => pef%eflx_anthro

    sabg_lyr       => pef%sabg_lyr
    h2osno         => cws%h2osno
    h2osoi_liq     => cws%h2osoi_liq
    h2osoi_ice     => cws%h2osoi_ice

    ! Get step size

    dtime = get_step_size()

    ! Compute ground surface and soil temperatures

    ! Thermal conductivity and Heat capacity

    call SoilThermProp(lbc, ubc, num_nolakec, filter_nolakec, tk, cv)

    ! Net ground heat flux into the surface and its temperature derivative
    ! Added a pfts loop here to get the average of hs and dhsdT over 
    ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       lwrad_emit(c)  =    emg(c) * sb * t_grnd(c)**4
       dlwrad_emit(c) = 4._r8*emg(c) * sb * t_grnd(c)**3
    end do

    hs(lbc:ubc) = 0._r8
    dhsdT(lbc:ubc) = 0._r8
    do pi = 1,max_pft_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             l = plandunit(p)
             g = pgridcell(p)

             ! Note: Some glacier_mec pfts may have zero weight
             if (pwtgcell(p)>0._r8 .or. ltype(l)==istice_mec) then
                if (ltype(l) /= isturb) then
                   eflx_gnet(p) = sabg(p) + dlrad(p) &
                                  + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - lwrad_emit(c) &
                                  - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
                else
                   ! For urban columns we use the net longwave radiation (eflx_lwrad_net) because of 
                   ! interactions between urban columns.
                   
                   ! All wasteheat and traffic flux goes into canyon floor
                   if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                      eflx_wasteheat_pft(p) = eflx_wasteheat(l)/(1._r8-wtlunit_roof(l))
                      eflx_heat_from_ac_pft(p) = eflx_heat_from_ac(l)/(1._r8-wtlunit_roof(l))
                      eflx_traffic_pft(p) = eflx_traffic(l)/(1._r8-wtlunit_roof(l))
                   else
                      eflx_wasteheat_pft(p) = 0._r8
                      eflx_heat_from_ac_pft(p) = 0._r8
                      eflx_traffic_pft(p) = 0._r8
                   end if
                   ! Include transpiration term because needed for previous road
                   ! and include wasteheat and traffic flux
                   eflx_gnet(p) = sabg(p) + dlrad(p)  &
                                  - eflx_lwrad_net(p) &
                                  - (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + qflx_tran_veg(p)*hvap) &
                                  + eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
                   eflx_anthro(p) = eflx_wasteheat_pft(p) + eflx_traffic_pft(p)
                end if
                dgnetdT(p) = - cgrnd(p) - dlwrad_emit(c)
                hs(c) = hs(c) + eflx_gnet(p) * pwtcol(p)
                dhsdT(c) = dhsdT(c) + dgnetdT(p) * pwtcol(p)
             end if

          end if
       end do
    end do

    !       Additional calculations with SNICAR: 
    !       Set up tridiagonal matrix in a new manner. There is now 
    !       absorbed solar radiation in each snow layer, instead of 
    !       only the surface. Following the current implementation, 
    !       absorbed solar flux should be: S + ((delS/delT)*dT), 
    !       where S is absorbed radiation, and T is temperature. Now, 
    !       assume delS/delT is zero, then it is OK to just add S 
    !       to each layer

    ! Initialize:
    sabg_lyr_col(lbc:ubc,-nlevsno+1:1) = 0._r8
    hs_top(lbc:ubc) = 0._r8

    do pi = 1,max_pft_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          lyr_top = snl(c) + 1
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             l = plandunit(p)
             if (pwtgcell(p)>0._r8 .or. ltype(l)==istice_mec) then
                g = pgridcell(p)
                if (ltype(l) /= isturb )then

                   eflx_gnet_top = sabg_lyr(p,lyr_top) + dlrad(p) + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) &
                        - lwrad_emit(c) - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

                   hs_top(c) = hs_top(c) + eflx_gnet_top*pwtcol(p)
             
                   do j = lyr_top,1,1
                      sabg_lyr_col(c,j) = sabg_lyr_col(c,j) + sabg_lyr(p,j) * pwtcol(p)
                   enddo
                else

                   hs_top(c) = hs_top(c) + eflx_gnet(p)*pwtcol(p)

                   sabg_lyr_col(c,lyr_top) = sabg_lyr_col(c,lyr_top) + sabg(p) * pwtcol(p)
             
                endif
             endif

          endif
       enddo
    enddo

    ! Restrict internal building temperature to between min and max
    ! and determine if heating or air conditioning is on
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (ltype(l) == isturb) then
          cool_on(l) = .false. 
          heat_on(l) = .false. 
          if (t_building(l) > t_building_max(l)) then
            t_building(l) = t_building_max(l)
            cool_on(l) = .true.
            heat_on(l) = .false.
          else if (t_building(l) < t_building_min(l)) then
            t_building(l) = t_building_min(l)
            cool_on(l) = .false.
            heat_on(l) = .true.
          end if
       end if
    end do

    ! Determine heat diffusion through the layer interface and factor used in computing
    ! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
    ! matrix and solve system

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (j >= snl(c)+1) then
             if (j == snl(c)+1) then
                if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall .or. ctype(c)==icol_roof) then
                  fact(c,j) = dtime/cv(c,j)
                else
                  fact(c,j) = dtime/cv(c,j) * dz(c,j) / (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
                end if
                fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
             else if (j <= nlevgrnd-1) then
                fact(c,j) = dtime/cv(c,j)
                fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                dzm     = (z(c,j)-z(c,j-1))
             else if (j == nlevgrnd) then
                fact(c,j) = dtime/cv(c,j)

                ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                ! building temperature. (See Oleson urban notes of 6/18/03).
                if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall .or. ctype(c)==icol_roof) then
                   fn(c,j) = tk(c,j) * (t_building(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
                else
                   fn(c,j) = eflx_bot(c)
                end if
             end if
          end if
       enddo
    end do

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (j >= snl(c)+1) then
             if (j == snl(c)+1) then
                dzp     = z(c,j+1)-z(c,j)
                at(c,j) = 0._r8
                bt(c,j) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                ct(c,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                ! changed hs to hs_top
                rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs_top(c) - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
             else if (j <= nlevgrnd-1) then
                dzm     = (z(c,j)-z(c,j-1))
                dzp     = (z(c,j+1)-z(c,j))
                at(c,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                bt(c,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                ct(c,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp

                ! if this is a snow layer or the top soil layer,
                ! add absorbed solar flux to factor 'rt'
                if (j <= 1) then
                   rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) ) + (fact(c,j)*sabg_lyr_col(c,j))
                else
                   rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                endif

             else if (j == nlevgrnd) then

                ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                ! building temperature. (See Oleson urban notes of 6/18/03).
                if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall .or. ctype(c)==icol_roof) then
                   dzm     = ( z(c,j)-z(c,j-1))
                   dzp     = (zi(c,j)-z(c,j))
                   at(c,j) =   - (1._r8-cnfac)*fact(c,j)*(tk(c,j-1)/dzm)
                   bt(c,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j-1)/dzm + tk(c,j)/dzp)
                   ct(c,j) = 0._r8
                   rt(c,j) = t_soisno(c,j) + fact(c,j)*( fn(c,j) - cnfac*fn(c,j-1) )
                else
                   dzm     = (z(c,j)-z(c,j-1))
                   at(c,j) =   - (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                   bt(c,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                   ct(c,j) = 0._r8
                   rt(c,j) = t_soisno(c,j) - cnfac*fact(c,j)*fn(c,j-1) + fact(c,j)*fn(c,j)
                end if
             end if

          end if
       enddo
    end do

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       jtop(c) = snl(c) + 1
    end do
    call Tridiagonal(lbc, ubc, -nlevsno+1, nlevgrnd, jtop, num_nolakec, filter_nolakec, &
                     at, bt, ct, rt, t_soisno(lbc:ubc,-nlevsno+1:nlevgrnd))

    ! Melting or Freezing

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (j >= snl(c)+1) then
             if (j <= nlevgrnd-1) then
                fn1(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
             else if (j == nlevgrnd) then

                ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                ! building temperature. (See Oleson urban notes of 6/18/03).
                ! Note new formulation for fn, this will be used below in brr computation
                if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall .or. ctype(c)==icol_roof) then
                   fn1(c,j) = tk(c,j) * (t_building(l) - t_soisno(c,j))/(zi(c,j) - z(c,j))
                   fn(c,j)  = tk(c,j) * (t_building(l) - tssbef(c,j))/(zi(c,j) - z(c,j))
                else
                   fn1(c,j) = 0._r8
                end if
             end if
          end if
       end do
    end do

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (ltype(l) == isturb) then
         eflx_building_heat(c) = cnfac*fn(c,nlevurb) + (1-cnfac)*fn1(c,nlevurb)
         if (cool_on(l)) then
           eflx_urban_ac(c) = abs(eflx_building_heat(c))
           eflx_urban_heat(c) = 0._r8
         else if (heat_on(l)) then
           eflx_urban_ac(c) = 0._r8
           eflx_urban_heat(c) = abs(eflx_building_heat(c))
         else
           eflx_urban_ac(c) = 0._r8
           eflx_urban_heat(c) = 0._r8
         end if
       end if
    end do

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (j >= snl(c)+1) then
             if (j == snl(c)+1) then
                brr(c,j) = cnfac*fn(c,j) + (1._r8-cnfac)*fn1(c,j)
             else
                brr(c,j) = cnfac*(fn(c,j)-fn(c,j-1)) + (1._r8-cnfac)*(fn1(c,j)-fn1(c,j-1))
             end if
          end if
       end do
    end do

    call PhaseChange (lbc, ubc, num_nolakec, filter_nolakec, fact, brr, hs, dhsdT, xmf, hs_top, sabg_lyr_col)

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       t_grnd(c) = t_soisno(c,snl(c)+1)
    end do


! Initialize soil heat content
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (ltype(l) /= isturb) then
         hc_soisno(c) = 0._r8
         hc_soi(c)    = 0._r8
       end if
       eflx_fgr12(c)= 0._r8
    end do

! Calculate soil heat content and soil plus snow heat content
    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          eflx_fgr12(c) = -cnfac*fn(c,1) - (1._r8-cnfac)*fn1(c,1)
          if (ltype(l) /= isturb) then
            if (j >= snl(c)+1) then
               hc_soisno(c) = hc_soisno(c) + cv(c,j)*t_soisno(c,j) / 1.e6_r8
            endif
            if (j >= 1) then
               hc_soi(c) = hc_soi(c) + cv(c,j)*t_soisno(c,j) / 1.e6_r8
            end if
          end if
       end do
    end do

  end subroutine SoilTemperature

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilThermProp
!
! !INTERFACE:
  subroutine SoilThermProp (lbc, ubc,  num_nolakec, filter_nolakec, tk, cv)
!
! !DESCRIPTION:
! Calculation of thermal conductivities and heat capacities of
! snow/soil layers
! (1) The volumetric heat capacity is calculated as a linear combination
!     in terms of the volumetric fraction of the constituent phases.
!
! (2) The thermal conductivity of soil is computed from the algorithm of
!     Johansen (as reported by Farouki 1981), and of snow is from the
!     formulation used in SNTHERM (Jordan 1991).
! The thermal conductivities at the interfaces between two neighboring
! layers (j, j+1) are derived from an assumption that the flux across
! the interface is equal to that from the node j to the interface and the
! flux from the interface to the node j+1.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
    use clm_varcon  , only : denh2o, denice, tfrz, tkwat, tkice, tkair, &
                             cpice,  cpliq,  istice, istice_mec, istwet, &
                             icol_roof, icol_sunwall, icol_shadewall, &
                             icol_road_perv, icol_road_imperv
    use clm_varpar  , only : nlevsno, nlevgrnd, nlevurb, nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                       ! column bounds
    integer , intent(in)  :: num_nolakec                    ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(ubc-lbc+1)      ! column filter for non-lake points
    real(r8), intent(out) :: cv(lbc:ubc,-nlevsno+1:nlevgrnd)! heat capacity [J/(m2 K)]
    real(r8), intent(out) :: tk(lbc:ubc,-nlevsno+1:nlevgrnd)! thermal conductivity [W/(m K)]
!
! !CALLED FROM:
! subroutine SoilTemperature in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/13/02, Peter Thornton: migrated to new data structures
! 7/01/03, Mariana Vertenstein: migrated to vector code
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer , pointer :: ctype(:)         ! column type
    integer , pointer :: clandunit(:)     ! column's landunit
    integer , pointer :: ltype(:)       ! landunit type
    integer , pointer :: snl(:)           ! number of snow layers
    real(r8), pointer :: h2osno(:)        ! snow water (mm H2O)
!
! local pointers to original implicit in arrays
!
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: tksatu(:,:)      ! thermal conductivity, saturated soil [W/m-K]
    real(r8), pointer :: tkmg(:,:)        ! thermal conductivity, soil minerals  [W/m-K]
    real(r8), pointer :: tkdry(:,:)       ! thermal conductivity, dry soil (W/m/Kelvin)
    real(r8), pointer :: csol(:,:)        ! heat capacity, soil solids (J/m**3/Kelvin)
    real(r8), pointer :: dz(:,:)          ! layer depth (m)
    real(r8), pointer :: zi(:,:)          ! interface level below a "z" level (m)
    real(r8), pointer :: z(:,:)           ! layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(r8), pointer :: tk_wall(:,:)     ! thermal conductivity of urban wall
    real(r8), pointer :: tk_roof(:,:)     ! thermal conductivity of urban roof
    real(r8), pointer :: tk_improad(:,:)  ! thermal conductivity of urban impervious road
    real(r8), pointer :: cv_wall(:,:)     ! thermal conductivity of urban wall
    real(r8), pointer :: cv_roof(:,:)     ! thermal conductivity of urban roof
    real(r8), pointer :: cv_improad(:,:)  ! thermal conductivity of urban impervious road
    integer,  pointer :: nlev_improad(:)  ! number of impervious road layers

!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: l,c,j                     ! indices
    integer  :: fc                        ! lake filtered column indices
    real(r8) :: bw                        ! partial density of water (ice + liquid)
    real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
    real(r8) :: dke                       ! kersten number
    real(r8) :: fl                        ! fraction of liquid or unfrozen water to total water
    real(r8) :: satw                      ! relative total water content of soil.
    real(r8) :: thk(lbc:ubc,-nlevsno+1:nlevgrnd) ! thermal conductivity of layer
    real(r8) :: thk_bedrock = 3.0_r8      ! thermal conductivity of 'typical' saturated granitic rock 
                                          ! (Clauser and Huenges, 1995)(W/m/K)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype    => lun%itype

    ! Assign local pointers to derived subtypes components (column-level)

    ctype      => col%itype
    clandunit  => col%landunit
    snl        => cps%snl
    h2osno     => cws%h2osno
    watsat     => cps%watsat
    tksatu     => cps%tksatu
    tkmg       => cps%tkmg
    tkdry      => cps%tkdry
    csol       => cps%csol
    dz         => cps%dz
    zi         => cps%zi
    z          => cps%z
    t_soisno   => ces%t_soisno
    h2osoi_liq => cws%h2osoi_liq
    h2osoi_ice => cws%h2osoi_ice
    tk_wall    => lps%tk_wall
    tk_roof    => lps%tk_roof
    tk_improad => lps%tk_improad
    cv_wall    => lps%cv_wall
    cv_roof    => lps%cv_roof
    cv_improad => lps%cv_improad
    nlev_improad => lps%nlev_improad

    ! Thermal conductivity of soil from Farouki (1981)
    ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
    ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
    ! 41, 1011-1026.

    do j = -nlevsno+1,nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)

          ! Only examine levels from 1->nlevgrnd
          if (j >= 1) then    
             l = clandunit(c)
             if (ctype(c) == icol_sunwall .OR. ctype(c) == icol_shadewall) then
                thk(c,j) = tk_wall(l,j)
             else if (ctype(c) == icol_roof) then
                thk(c,j) = tk_roof(l,j)
             else if (ctype(c) == icol_road_imperv .and. j >= 1 .and. j <= nlev_improad(l)) then
                thk(c,j) = tk_improad(l,j)
             elseif (ltype(l) /= istwet .AND. ltype(l) /= istice  &
                                        .AND. ltype(l) /= istice_mec) then
                satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
                satw = min(1._r8, satw)
                if (satw > .1e-6_r8) then
                   fl = h2osoi_liq(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j))
                   if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                      dke = max(0._r8, log10(satw) + 1.0_r8)
                      dksat = tksatu(c,j)
                   else                               ! Frozen soil
                      dke = satw
                      dksat = tkmg(c,j)*0.249_r8**(fl*watsat(c,j))*2.29_r8**watsat(c,j)
                   endif
                   thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
                else
                   thk(c,j) = tkdry(c,j)
                endif
                if (j > nlevsoi) thk(c,j) = thk_bedrock
             else if (ltype(l) == istice .OR. ltype(l) == istice_mec) then
                thk(c,j) = tkwat
                if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
             else if (ltype(l) == istwet) then                         
                if (j > nlevsoi) then 
                   thk(c,j) = thk_bedrock
                else
                   thk(c,j) = tkwat
                   if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
                endif
             endif
          endif

          ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
          ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
          if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then  
             bw = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/dz(c,j)
             thk(c,j) = tkair + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkair)
          end if

       end do
    end do

    ! Thermal conductivity at the layer interface

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1 .AND. j <= nlevgrnd-1) then
             tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                       /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
          else if (j == nlevgrnd) then

             ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
             ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
             ! building temperature. (See Oleson urban notes of 6/18/03).
             if (ctype(c)==icol_sunwall .OR. ctype(c)==icol_shadewall .OR. ctype(c)==icol_roof) then
                tk(c,j) = thk(c,j)
             else
                tk(c,j) = 0._r8
             end if
          end if
       end do
    end do

    ! Soil heat capacity, from de Vires (1963)
    ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
    ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
    ! 41, 1011-1026.

    do j = 1, nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ctype(c)==icol_sunwall .OR. ctype(c)==icol_shadewall) then
             cv(c,j) = cv_wall(l,j) * dz(c,j)
          else if (ctype(c) == icol_roof) then
             cv(c,j) = cv_roof(l,j) * dz(c,j)
          else if (ctype(c) == icol_road_imperv .and. j >= 1 .and. j <= nlev_improad(l)) then
             cv(c,j) = cv_improad(l,j) * dz(c,j)
          elseif (ltype(l) /= istwet .AND. ltype(l) /= istice  &
                                     .AND. ltype(l) /= istice_mec) then
             cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) + (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          else if (ltype(l) == istwet) then 
             cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
             if (j > nlevsoi) cv(c,j) = csol(c,j)*dz(c,j)
          else if (ltype(l) == istice .OR. ltype(l) == istice_mec) then 
             cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          endif
          if (j == 1) then
             if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
                cv(c,j) = cv(c,j) + cpice*h2osno(c)
             end if
          end if
       enddo
    end do

    ! Snow heat capacity

    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
             cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
          end if
       end do
    end do

  end subroutine SoilThermProp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: PhaseChange
!
! !INTERFACE:
  subroutine PhaseChange (lbc, ubc, num_nolakec, filter_nolakec, fact, &
                          brr, hs, dhsdT, xmf, hs_top, sabg_lyr_col)
!
! !DESCRIPTION:
! Calculation of the phase change within snow and soil layers:
! (1) Check the conditions for which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e. melting),
!     or the layer temperature is less than the freezing point
!     and the liquid water mass is greater than the allowable supercooled 
!     liquid water calculated from freezing point depression (i.e. freezing).
! (2) Assess the rate of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point.
! (3) Re-adjust the ice and liquid mass, and the layer temperature
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_time_manager, only : get_step_size
    use clm_varcon  , only : tfrz, hfus, grav, istsoil, istice_mec, isturb, icol_road_perv
    use clm_varcon  , only : istcrop
    use clm_varpar  , only : nlevsno, nlevgrnd

!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                             ! column bounds
    integer , intent(in) :: num_nolakec                          ! number of column non-lake points in column filter
    integer , intent(in) :: filter_nolakec(ubc-lbc+1)            ! column filter for non-lake points
    real(r8), intent(in) :: brr   (lbc:ubc, -nlevsno+1:nlevgrnd) ! temporary
    real(r8), intent(in) :: fact  (lbc:ubc, -nlevsno+1:nlevgrnd) ! temporary
    real(r8), intent(in) :: hs    (lbc:ubc)                      ! net ground heat flux into the surface
    real(r8), intent(in) :: dhsdT (lbc:ubc)                      ! temperature derivative of "hs"
    real(r8), intent(out):: xmf   (lbc:ubc)                      ! total latent heat of phase change
    real(r8), intent(in) :: hs_top(lbc:ubc)                      ! net heat flux into the top snow layer [W/m2]
    real(r8), intent(in) :: sabg_lyr_col(lbc:ubc,-nlevsno+1:1)   ! absorbed solar radiation (col,lyr) [W/m2]

!
! !CALLED FROM:
! subroutine SoilTemperature in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/14/02, Peter Thornton: Migrated to new data structures.
! 7/01/03, Mariana Vertenstein: Migrated to vector code
! 04/25/07 Keith Oleson: CLM3.5 Hydrology
! 03/28/08 Mark Flanner: accept new arguments and calculate freezing rate of h2o in snow
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    real(r8), pointer :: qflx_snow_melt(:)    ! net snow melt
    integer , pointer :: snl(:)           !number of snow layers
    real(r8), pointer :: h2osno(:)        !snow water (mm H2O)
    integer , pointer :: ltype(:)         !landunit type
    integer , pointer :: clandunit(:)     !column's landunit
    integer , pointer :: ctype(:)         !column type
!
! local pointers to original implicit inout scalars
!
    real(r8), pointer :: snowdp(:)        !snow height (m)
!
! local pointers to original implicit out scalars
!
    real(r8), pointer :: qflx_snomelt(:)  !snow melt (mm H2O /s)
    real(r8), pointer :: eflx_snomelt(:)  !snow melt heat flux (W/m**2)
    real(r8), pointer :: eflx_snomelt_u(:)!urban snow melt heat flux (W/m**2)
    real(r8), pointer :: eflx_snomelt_r(:)!rural snow melt heat flux (W/m**2)
    real(r8), pointer :: qflx_snofrz_lyr(:,:)  !snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
    real(r8), pointer :: qflx_snofrz_col(:)  !column-integrated snow freezing rate (positive definite) [kg m-2 s-1]
    real(r8), pointer :: qflx_glcice(:)   !flux of new glacier ice (mm H2O/s) [+ = ice grows]
    real(r8), pointer :: qflx_glcice_melt(:)  !ice melt (positive definite) (mm H2O/s)
!
! local pointers to original implicit in arrays
!
    real(r8), pointer :: h2osoi_liq(:,:)  !liquid water (kg/m2) (new)
    real(r8), pointer :: h2osoi_ice(:,:)  !ice lens (kg/m2) (new)
    real(r8), pointer :: tssbef(:,:)      !temperature at previous time step [K]
    real(r8), pointer :: sucsat(:,:)      !minimum soil suction (mm)
    real(r8), pointer :: watsat(:,:)      !volumetric soil water at saturation (porosity)
    real(r8), pointer :: bsw(:,:)         !Clapp and Hornberger "b"
    real(r8), pointer :: dz(:,:)          !layer thickness (m)
!
! local pointers to original implicit inout arrays
!
    real(r8), pointer :: t_soisno(:,:)    !soil temperature (Kelvin)
!
! local pointers to original implicit out arrays
!
    integer, pointer :: imelt(:,:)        !flag for melting (=1), freezing (=2), Not=0 (new)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: j,c,g,l                            !do loop index
    integer  :: fc                                 !lake filtered column indices
    real(r8) :: dtime                              !land model time step (sec)
    real(r8) :: heatr                              !energy residual or loss after melting or freezing
    real(r8) :: temp1                              !temporary variables [kg/m2]
    real(r8) :: hm(lbc:ubc,-nlevsno+1:nlevgrnd)    !energy residual [W/m2]
    real(r8) :: xm(lbc:ubc,-nlevsno+1:nlevgrnd)    !melting or freezing within a time step [kg/m2]
    real(r8) :: wmass0(lbc:ubc,-nlevsno+1:nlevgrnd)!initial mass of ice and liquid (kg/m2)
    real(r8) :: wice0 (lbc:ubc,-nlevsno+1:nlevgrnd)!initial mass of ice (kg/m2)
    real(r8) :: wliq0 (lbc:ubc,-nlevsno+1:nlevgrnd)!initial mass of liquid (kg/m2)
    real(r8) :: supercool(lbc:ubc,nlevgrnd)        !supercooled water in soil (kg/m2) 
    real(r8) :: propor                             !proportionality constant (-)
    real(r8) :: tinc                               !t(n+1)-t(n) (K)
    real(r8) :: smp                                !frozen water potential (mm)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    snl          => cps%snl
    h2osno       => cws%h2osno
    snowdp       => cps%snowdp
    qflx_snow_melt => cwf%qflx_snow_melt
    qflx_snomelt => cwf%qflx_snomelt
    eflx_snomelt => cef%eflx_snomelt
    eflx_snomelt_u => cef%eflx_snomelt_u
    eflx_snomelt_r => cef%eflx_snomelt_r
    h2osoi_liq   => cws%h2osoi_liq
    h2osoi_ice   => cws%h2osoi_ice
    imelt        => cps%imelt
    t_soisno     => ces%t_soisno
    tssbef       => ces%tssbef
    bsw          => cps%bsw
    sucsat       => cps%sucsat
    watsat       => cps%watsat
    dz           => cps%dz
    ctype        => col%itype
    clandunit    => col%landunit
    ltype        => lun%itype
    qflx_snofrz_lyr => cwf%qflx_snofrz_lyr
    qflx_snofrz_col => cwf%qflx_snofrz_col
    qflx_glcice  => cwf%qflx_glcice
    qflx_glcice_melt  => cwf%qflx_glcice_melt

    ! Get step size

    dtime = get_step_size()

    ! Initialization

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)

       qflx_snomelt(c) = 0._r8
       qflx_snow_melt(c) = 0._r8
       xmf(c) = 0._r8
       qflx_snofrz_lyr(c,-nlevsno+1:0) = 0._r8
       qflx_snofrz_col(c) = 0._r8
       if (ltype(l)==istice_mec) then
          ! only need to initialize qflx_glcice_melt over ice_mec landunits, because
          ! those are the only places where it is computed
          qflx_glcice_melt(c) = 0._r8
       end if
    end do

    do j = -nlevsno+1,nlevgrnd       ! all layers
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Initialization
             imelt(c,j) = 0
             hm(c,j) = 0._r8
             xm(c,j) = 0._r8
             wice0(c,j) = h2osoi_ice(c,j)
             wliq0(c,j) = h2osoi_liq(c,j)
             wmass0(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

    do j = -nlevsno+1,0             ! snow layers 
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Melting identification
             ! If ice exists above melt point, melt some to liquid.
             if (h2osoi_ice(c,j) > 0._r8 .AND. t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
                t_soisno(c,j) = tfrz
             endif

             ! Freezing identification
             ! If liquid exists below melt point, freeze some to ice.
             if (h2osoi_liq(c,j) > 0._r8 .AND. t_soisno(c,j) < tfrz) then
                imelt(c,j) = 2
                t_soisno(c,j) = tfrz
             endif
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

    do j = 1,nlevgrnd             ! soil layers 
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (h2osoi_ice(c,j) > 0. .AND. t_soisno(c,j) > tfrz) then
             imelt(c,j) = 1
             t_soisno(c,j) = tfrz
          endif

          ! from Zhao (1997) and Koren (1999)
          supercool(c,j) = 0.0_r8
          if (ltype(l) == istsoil .or. ltype(l) == istcrop .or. ctype(c) == icol_road_perv) then
             if(t_soisno(c,j) < tfrz) then
                smp = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                supercool(c,j) = watsat(c,j)*(smp/sucsat(c,j))**(-1._r8/bsw(c,j))
                supercool(c,j) = supercool(c,j)*dz(c,j)*1000._r8       ! (mm)
             endif
          endif

          if (h2osoi_liq(c,j) > supercool(c,j) .AND. t_soisno(c,j) < tfrz) then
             imelt(c,j) = 2
             t_soisno(c,j) = tfrz
          endif

          ! If snow exists, but its thickness is less than the critical value (0.01 m)
          if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. j == 1) then
             if (t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
                t_soisno(c,j) = tfrz
             endif
          endif
       end do
    enddo

    do j = -nlevsno+1,nlevgrnd       ! all layers
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)

          if (j >= snl(c)+1) then

             ! Calculate the energy surplus and loss for melting and freezing
             if (imelt(c,j) > 0) then
                tinc = t_soisno(c,j)-tssbef(c,j)
                
                ! added unique cases for this calculation,
                ! to account for absorbed solar radiation in each layer
                if (j == snl(c)+1) then
                   ! top layer
                   hm(c,j) = hs_top(c) + dhsdT(c)*tinc + brr(c,j) - tinc/fact(c,j)
                elseif (j <= 1) then
                   ! snow layer or top soil layer (where sabg_lyr_col is defined)
                   hm(c,j) = brr(c,j) - tinc/fact(c,j) + sabg_lyr_col(c,j)
                else
                   ! soil layer
                   hm(c,j) = brr(c,j) - tinc/fact(c,j)
                endif

             endif

             ! These two errors were checked carefully (Y. Dai).  They result from the
             ! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
             if (imelt(c,j) == 1 .AND. hm(c,j) < 0._r8) then
                hm(c,j) = 0._r8
                imelt(c,j) = 0
             endif
             if (imelt(c,j) == 2 .AND. hm(c,j) > 0._r8) then
                hm(c,j) = 0._r8
                imelt(c,j) = 0
             endif

             ! The rate of melting and freezing

             if (imelt(c,j) > 0 .and. abs(hm(c,j)) > 0._r8) then
                xm(c,j) = hm(c,j)*dtime/hfus                           ! kg/m2

                ! If snow exists, but its thickness is less than the critical value
                ! (1 cm). Note: more work is needed to determine how to tune the
                ! snow depth for this case
                if (j == 1) then
                   if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. xm(c,j) > 0._r8) then
                      temp1 = h2osno(c)                           ! kg/m2
                      h2osno(c) = max(0._r8,temp1-xm(c,j))
                      propor = h2osno(c)/temp1
                      snowdp(c) = propor * snowdp(c)
                      heatr = hm(c,j) - hfus*(temp1-h2osno(c))/dtime   ! W/m2
                      if (heatr > 0._r8) then
                         xm(c,j) = heatr*dtime/hfus                    ! kg/m2
                         hm(c,j) = heatr                               ! W/m2
                      else
                         xm(c,j) = 0._r8
                         hm(c,j) = 0._r8
                      endif
                      qflx_snomelt(c) = max(0._r8,(temp1-h2osno(c)))/dtime   ! kg/(m2 s)
                      xmf(c) = hfus*qflx_snomelt(c)
                      qflx_snow_melt(c) = qflx_snomelt(c)
                   endif
                endif

                heatr = 0._r8
                if (xm(c,j) > 0._r8) then
                   h2osoi_ice(c,j) = max(0._r8, wice0(c,j)-xm(c,j))
                   heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                else if (xm(c,j) < 0._r8) then
                   if (j <= 0) then
                      h2osoi_ice(c,j) = min(wmass0(c,j), wice0(c,j)-xm(c,j))  ! snow
                   else
                      if (wmass0(c,j) < supercool(c,j)) then
                         h2osoi_ice(c,j) = 0._r8
                      else
                         h2osoi_ice(c,j) = min(wmass0(c,j) - supercool(c,j),wice0(c,j)-xm(c,j))
                      endif
                   endif
                   heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                endif

                h2osoi_liq(c,j) = max(0._r8,wmass0(c,j)-h2osoi_ice(c,j))

                if (abs(heatr) > 0._r8) then
                   if (j > snl(c)+1) then
                      t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr
                   else
                      t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr/(1._r8-fact(c,j)*dhsdT(c))
                   endif
                   if (j <= 0) then    ! snow
                      if (h2osoi_liq(c,j)*h2osoi_ice(c,j)>0._r8) t_soisno(c,j) = tfrz
                   end if
                endif

                xmf(c) = xmf(c) + hfus * (wice0(c,j)-h2osoi_ice(c,j))/dtime

                if (imelt(c,j) == 1 .AND. j < 1) then
                   qflx_snomelt(c) = qflx_snomelt(c) + max(0._r8,(wice0(c,j)-h2osoi_ice(c,j)))/dtime
                endif

                ! layer freezing mass flux (positive):
                if (imelt(c,j) == 2 .AND. j < 1) then
                   qflx_snofrz_lyr(c,j) = max(0._r8,(h2osoi_ice(c,j)-wice0(c,j)))/dtime
                endif

             endif
          endif   ! end of snow layer if-block

       ! For glacier_mec columns, compute negative ice flux from melted ice.
       ! Note that qflx_glcice can also include a positive component from excess snow,
       !  as computed in Hydrology2Mod.F90.

          l = clandunit(c)
          if (ltype(l)==istice_mec) then

             if (j>=1 .and. h2osoi_liq(c,j) > 0._r8) then   ! ice layer with meltwater
                ! melting corresponds to a negative ice flux
                qflx_glcice_melt(c) = qflx_glcice_melt(c) + h2osoi_liq(c,j)/dtime
                qflx_glcice(c) = qflx_glcice(c) - h2osoi_liq(c,j)/dtime

                ! convert layer back to pure ice by "borrowing" ice from below the column
                h2osoi_ice(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
                h2osoi_liq(c,j) = 0._r8

             endif  ! liquid water is present
          endif     ! istice_mec

       end do   ! end of column-loop
    enddo   ! end of level-loop

    ! Needed for history file output

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       eflx_snomelt(c) = qflx_snomelt(c) * hfus
       l = clandunit(c)
       if (ltype(l) == isturb) then
         eflx_snomelt_u(c) = eflx_snomelt(c)
       else if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         eflx_snomelt_r(c) = eflx_snomelt(c)
       end if
    end do

    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          qflx_snofrz_col(c) = qflx_snofrz_col(c) + qflx_snofrz_lyr(c,j)
       end do
    end do

  end subroutine PhaseChange


end module SoilTemperatureMod
