module ch4Mod
#ifdef LCH4

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ch4Mod
!
! !DESCRIPTION:
! Module holding routines to calculate methane fluxes
! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in ch4_tran.
!
! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varpar   , only: nlevsoi, ngases, nlevsno
  use clm_varcon   , only: denh2o, denice, tfrz, grav, spval, rgas
  use clm_varcon   , only: catomw, s_con, d_con_w, d_con_g, c_h_inv, kh_theta, kh_tbase
  use clm_atmlnd   , only : clm_a2l, clm_l2a
  use clm_time_manager, only : get_step_size, get_nstep
  use clm_varctl   , only : iulog
  use abortutils   , only : endrun

  implicit none
  save
  private
  real(r8) :: f_sat = 0.95_r8 ! volumetric soil water defining top of water table or where production is allowed

  ! Non-tunable constants
  real(r8) :: rgasm  ! J/mol.K; rgas / 1000; will be set below
  real(r8), parameter :: rgasLatm = 0.0821_r8 ! L.atm/mol.K

! !PUBLIC MEMBER FUNCTIONS:
  public :: ch4
  private :: ch4_prod
  private :: ch4_oxid
  private :: ch4_aere
  private :: ch4_ebul
  private :: ch4_tran
  private :: ch4annualupdate
  private :: get_jwt
!
! !REVISION HISTORY:
! Bill Riley & Zack Subin, Created 2008-2011
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4
!
! !INTERFACE:
subroutine ch4 (lbg, ubg, lbl, ubl, lbc, ubc, lbp, ubp, num_soilc, filter_soilc, num_lakec, filter_lakec, &
                num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Driver for the methane emissions model
!
! !USES:
   use clmtype
   use subgridAveMod , only : p2c, c2g
   use clm_varpar, only : nlevgrnd, nlevdecomp
   use pftvarcon,  only : noveg
   use ch4varcon,  only : replenishlakec, redoxlag, allowlakeprod, redoxlag_vertical
   use clm_varcon, only : secspday
   use ch4varcon,  only : ch4offline, fin_use_fsat, atmch4
   use nanmod   ,  only : nan, bigint
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbg, ubg           ! grid-index bounds
   integer, intent(in) :: lbl, ubl           ! landunit-level index bounds
   integer, intent(in) :: lbc, ubc           ! column-index bounds
   integer, intent(in) :: lbp, ubp           ! pft-level index bounds
   integer, intent(in) :: num_soilc          ! number of column soil points in column filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1)    ! column filter for soil points
   integer, intent(in) :: num_lakec          ! number of column lake points in column filter
   integer, intent(in) :: filter_lakec(ubc-lbc+1)    ! column filter for lake points
   integer, intent(in) :: num_soilp          ! number of soil points in pft filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1)    ! pft filter for soil points

!
! !CALLED FROM:
! driver.F90
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
! local pointers to implicit in variables
   integer , pointer :: cgridcell(:)   ! gridcell of corresponding column
   real(r8), pointer :: forc_t(:)      ! atmospheric temperature (Kelvin)
   real(r8), pointer :: forc_pbot(:)   ! atmospheric pressure (Pa)
   real(r8), pointer :: forc_pco2(:)   ! CO2 partial pressure (Pa)
   real(r8), pointer :: forc_pch4(:)   ! CH4 partial pressure (Pa)
   real(r8), pointer :: forc_po2(:)    ! O2 partial pressure (Pa)
   real(r8), pointer :: dz(:,:)        ! layer thickness (m)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: zi(:,:)        ! interface level below a "z" level (m)
   real(r8), pointer :: z(:,:)         ! layer depth (m) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: rootfr(:,:)    ! fraction of roots in each soil layer  (nlevgrnd)
   real(r8), pointer :: grnd_ch4_cond(:) ! tracer conductance for boundary layer [m/s]
   real(r8), pointer :: pwtc(:)        ! weight (relative to column) 
   integer , pointer :: ivt(:)         ! pft vegetation type
   integer , pointer :: pcolumn(:)     ! index into column level quantities
   real(r8), pointer :: zwt0(:)        ! decay factor for finundated (m)
   real(r8), pointer :: f0(:)          ! maximum gridcell fractional inundated area
   real(r8), pointer :: p3(:)          ! coefficient for qflx_surf_lag for finunated (s/mm)
   real(r8), pointer :: zwt(:)         ! water table depth (m) (from SoilHydrology)
   real(r8), pointer :: zwt_perched(:) ! perched water table depth (m)
   real(r8), pointer :: conc_o2_sat(:,:) ! O2 conc  in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: qflx_surf(:)   ! surface runoff (mm H2O /s)
   real(r8), pointer :: latdeg(:)      ! latitude (degrees)
   real(r8), pointer :: frac_h2osfc(:) ! fraction of ground covered by surface water (0 to 1)
   logical , pointer :: cactive(:)     ! true=>do computations on this column (see reweightMod for details)



! local pointers to implicit in/out variables
   real(r8), pointer :: fsat_bef(:)    ! finundated from previous timestep
   real(r8), pointer :: c_atm(:,:)     ! CH4, O2, CO2 atmospheric conc  (mol/m3)
   real(r8), pointer :: flux_ch4(:)    ! gridcell CH4 flux to atm. (kg C/m**2/s)
   real(r8), pointer :: ch4_surf_diff_sat(:)    ! CH4 surface flux (mol/m2/s)
   real(r8), pointer :: ch4_surf_diff_unsat(:)  ! CH4 surface flux (mol/m2/s)
   real(r8), pointer :: ch4_surf_diff_lake(:)   ! CH4 surface flux (mol/m2/s)
   real(r8), pointer :: ch4_surf_aere_sat(:)    ! Total column CH4 aerenchyma (mol/m2/s)
   real(r8), pointer :: ch4_surf_aere_unsat(:)  ! Total column CH4 aerenchyma (mol/m2/s)
   real(r8), pointer :: ch4_surf_ebul_sat(:)    ! CH4 ebullition to atmosphere (mol/m2/s)
   real(r8), pointer :: ch4_surf_ebul_unsat(:)  ! CH4 ebullition to atmosphere (mol/m2/s)
   real(r8), pointer :: ch4_surf_ebul_lake(:)   ! CH4 ebullition to atmosphere (mol/m2/s)
   real(r8), pointer :: ch4_oxid_depth_sat(:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_oxid_depth_unsat(:,:)! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_oxid_depth_lake(:,:)! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_prod_depth_sat(:,:) ! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: ch4_prod_depth_unsat(:,:)   ! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: ch4_prod_depth_lake(:,:)! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: ch4co2f(:)              ! gridcell CO2 production from CH4 oxidation (g C/m**2/s)
                                                ! needed to adjust nee for CN
   real(r8), pointer :: nem(:)                  ! gridcell average net methane correction to CO2 flux (g C/m^2/s)
   real(r8), pointer :: lake_soilc(:,:)         ! total soil organic matter found in level (g C / m^3) (nlevsoi)
                                                ! needed for lakes not using CN. initialized to total organic
                                                ! content read in from soil properties
   real(r8), pointer :: conc_ch4_sat(:,:)       ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_ch4_unsat(:,:)     ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_ch4_lake(:,:)      ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_o2_lake(:,:)       ! O2 conc  in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: ch4_dfsat_flux(:)       ! CH4 flux to atm due to decreasing finundated (kg C/m^2/s) [+]
   real(r8), pointer :: ch4prodg(:)             ! gridcell average CH4 production (g C/m^2/s)
   real(r8), pointer :: zwt_ch4_unsat(:)        ! depth of water table for unsaturated fraction (m)
   real(r8), pointer :: rootfr_col(:,:)         ! fraction of roots in each soil layer  (nlevgrnd)
   real(r8), pointer :: grnd_ch4_cond_col(:)    ! tracer conductance for boundary layer [m/s]
   real(r8), pointer :: totcolch4(:)            ! total methane in soil column (g C / m^2)
   real(r8), pointer :: finundated(:)           ! fractional inundated area in soil column (excluding dedicated wetland columns)
   real(r8), pointer :: qflx_surf_lag(:)        ! time-lagged surface runoff (mm H2O /s)
   real(r8), pointer :: finundated_lag(:)       ! time-lagged fractional inundated area
   real(r8), pointer :: layer_sat_lag(:,:)      ! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)


! !OTHER LOCAL VARIABLES:
   integer :: sat ! 0 = unsatured, 1 = saturated
   logical :: lake ! lake or not lake
   integer :: fc,c,g,j,fp,p               ! indices
   real(r8) :: dtime                      ! land model time step (sec)
   real(r8) :: dtime_ch4                  ! ch4 model time step (sec)
   integer  :: nstep
   integer  :: jwt(lbc:ubc)               !index of the soil layer right above the water table (-)
   real(r8) :: ch4_surf_flux_tot(lbc:ubc) !CH4 surface flux for column (kg C/m**2/s)
   real(r8) :: ch4_prod_tot(lbc:ubc)      !CH4 production for column (g C/m**2/s)
   real(r8) :: ch4_oxid_tot(lbc:ubc)      !CH4 oxidation for column (g C/m**2/s)
   real(r8) :: nem_col(lbc:ubc)           !net adjustment to atm. C flux from methane production (g C/m**2/s)
   real(r8) :: totalsat
   real(r8) :: totalunsat
   real(r8) :: dfsat
   real(r8) :: rootfraction(lbp:ubp, 1:nlevgrnd) 
   real(r8) :: totcolch4_bef(lbc:ubc) ! g C / m^2
   real(r8) :: errch4                 ! g C / m^2
   real(r8) :: zwt_actual
   real(r8) :: qflxlags               ! Time to lag qflx_surf_lag (s)
   real(r8) :: redoxlags              ! Redox time lag in s
   real(r8) :: redoxlags_vertical     ! Vertical redox lag time in s
   real(r8), parameter :: qflxlagd = 30._r8    !   days to lag qflx_surf_lag in the tropics (days)
   real(r8), parameter :: highlatfact = 2._r8  ! multiple of qflxlagd for high latitudes
   integer  :: dummyfilter(1)                  ! empty filter
   character(len=32) :: subname='ch4'          ! subroutine name

!EOP
!-----------------------------------------------------------------------
   ! Gridcell level pointers
   forc_t    => clm_a2l%forc_t
   forc_pbot => clm_a2l%forc_pbot
   c_atm     => clm3%g%gch4%c_atm
   forc_po2  => clm_a2l%forc_po2
   forc_pco2 => clm_a2l%forc_pco2
   forc_pch4 => clm_a2l%forc_pch4
   flux_ch4  => clm_l2a%flux_ch4
   ch4co2f   => clm3%g%gch4%ch4co2f
   ch4prodg  => clm3%g%gch4%ch4prodg
   nem       => clm3%g%gch4%nem
   latdeg    => clm3%g%latdeg
   ! Column level pointers
   cgridcell           => clm3%g%l%c%gridcell
   ch4_surf_diff_sat   => clm3%g%l%c%cch4%ch4_surf_diff_sat
   ch4_surf_diff_unsat => clm3%g%l%c%cch4%ch4_surf_diff_unsat
   ch4_surf_diff_lake  => clm3%g%l%c%cch4%ch4_surf_diff_lake
   ch4_surf_ebul_sat   => clm3%g%l%c%cch4%ch4_surf_ebul_sat
   ch4_surf_ebul_unsat => clm3%g%l%c%cch4%ch4_surf_ebul_unsat
   ch4_surf_ebul_lake  => clm3%g%l%c%cch4%ch4_surf_ebul_lake
   ch4_surf_aere_sat   => clm3%g%l%c%cch4%ch4_surf_aere_sat
   ch4_surf_aere_unsat => clm3%g%l%c%cch4%ch4_surf_aere_unsat
   fsat_bef            => clm3%g%l%c%cch4%fsat_bef
   ch4_oxid_depth_sat  => clm3%g%l%c%cch4%ch4_oxid_depth_sat
   ch4_oxid_depth_unsat=> clm3%g%l%c%cch4%ch4_oxid_depth_unsat
   ch4_oxid_depth_lake=> clm3%g%l%c%cch4%ch4_oxid_depth_lake
   ch4_prod_depth_sat  => clm3%g%l%c%cch4%ch4_prod_depth_sat
   ch4_prod_depth_unsat=> clm3%g%l%c%cch4%ch4_prod_depth_unsat
   ch4_prod_depth_lake => clm3%g%l%c%cch4%ch4_prod_depth_lake
   lake_soilc          => clm3%g%l%c%cch4%lake_soilc
   conc_ch4_sat        => clm3%g%l%c%cch4%conc_ch4_sat
   conc_ch4_unsat      => clm3%g%l%c%cch4%conc_ch4_unsat
   conc_ch4_lake       => clm3%g%l%c%cch4%conc_ch4_lake
   conc_o2_lake        => clm3%g%l%c%cch4%conc_o2_lake
   conc_o2_sat         => clm3%g%l%c%cch4%conc_o2_sat
   ch4_dfsat_flux      => clm3%g%l%c%cch4%ch4_dfsat_flux
   zwt_ch4_unsat       => clm3%g%l%c%cch4%zwt_ch4_unsat
   dz                  => clm3%g%l%c%cps%dz
   zi                  => clm3%g%l%c%cps%zi
   z                   => clm3%g%l%c%cps%z
   rootfr_col          => clm3%g%l%c%cps%pps_a%rootfr
   grnd_ch4_cond_col   => clm3%g%l%c%cps%pps_a%grnd_ch4_cond
   totcolch4           => clm3%g%l%c%cch4%totcolch4
   zwt0                => clm3%g%l%c%cps%zwt0
   f0                  => clm3%g%l%c%cps%f0
   p3                  => clm3%g%l%c%cps%p3
   finundated          => clm3%g%l%c%cws%finundated
   zwt                 => clm3%g%l%c%cws%zwt
   zwt_perched         => clm3%g%l%c%cws%zwt_perched
   qflx_surf           => clm3%g%l%c%cwf%qflx_surf
   qflx_surf_lag       => clm3%g%l%c%cch4%qflx_surf_lag
   finundated_lag      => clm3%g%l%c%cch4%finundated_lag
   layer_sat_lag       => clm3%g%l%c%cch4%layer_sat_lag
   frac_h2osfc         => clm3%g%l%c%cps%frac_h2osfc
   cactive             => clm3%g%l%c%active


   ! Pft level pointers
   rootfr              => clm3%g%l%c%p%pps%rootfr
   grnd_ch4_cond       => clm3%g%l%c%p%pps%grnd_ch4_cond
   pwtc                => clm3%g%l%c%p%wtcol
   ivt                 => clm3%g%l%c%p%itype 
   pcolumn             => clm3%g%l%c%p%column


   dtime = get_step_size()
   nstep = get_nstep()
   dtime_ch4 = dtime
   redoxlags = redoxlag*secspday ! days --> s
   redoxlags_vertical = redoxlag_vertical*secspday ! days --> s
   rgasm = rgas / 1000._r8

   jwt(:)            = bigint
   totcolch4_bef(:)  = nan

   ! Initialize local fluxes to zero: necessary for columns outside the filters because averaging up to gridcell will be done
   ch4_surf_flux_tot(:) = 0._r8
   ch4_prod_tot(:)      = 0._r8
   ch4_oxid_tot(:)      = 0._r8
   rootfraction(:,:)    = spval
   ! Adjustment to NEE for methane production - oxidation
   nem_col(:)           = 0._r8

   do g=lbg,ubg

      if (ch4offline) then
         forc_pch4(g) = atmch4*forc_pbot(g)
      else
         if (forc_pch4(g) == 0._r8) then
            write(iulog,*)'not using ch4offline, but methane concentration not passed from the atmosphere', &
                          'to land model! CLM Model is stopping.'
            call endrun( trim(subname)//' ERROR: Methane not being passed to atmosphere')
         end if
      end if

      c_atm(g,1) =  forc_pch4(g) / rgasm / forc_t(g) ! [mol/m3 air]
      c_atm(g,2) =  forc_po2(g) / rgasm / forc_t(g)  ! [mol/m3 air]
      c_atm(g,3) =  forc_pco2(g) / rgasm / forc_t(g) ! [mol/m3 air]
   end do

   ! Initialize CH4 balance and calculate finundated
   do fc = 1, num_soilc
      c = filter_soilc(fc)
      g = cgridcell(c)

      totcolch4_bef(c) = totcolch4(c)
      totcolch4(c) = 0._r8

      ! Update lagged surface runoff

      if (latdeg(g) < 45._r8) then
         qflxlags = qflxlagd * secspday ! 30 days
      else
         qflxlags = qflxlagd * secspday * highlatfact ! 60 days
      end if
      qflx_surf_lag(c) = qflx_surf_lag(c) * exp(-dtime/qflxlags) &
                       + qflx_surf(c) * (1._r8 - exp(-dtime/qflxlags))

!There may be ways to improve this for irrigated crop columns...
      if (fin_use_fsat) then
         finundated(c) = frac_h2osfc(c)
      else
         if (zwt0(c) > 0._r8) then
            if (zwt_perched(c) < z(c,nlevsoi)-1.e-5_r8 .and. zwt_perched(c) < zwt(c)) then
               zwt_actual = zwt_perched(c)
            else
               zwt_actual = zwt(c)
            end if
            finundated(c) = f0(c) * exp(-zwt_actual/zwt0(c)) + p3(c)*qflx_surf_lag(c)
         else
            finundated(c) = p3(c)*qflx_surf_lag(c)
         end if
      end if
      finundated(c) = max( min(finundated(c),1._r8), 0._r8)

      ! Update lagged finundated for redox calculation
      if (redoxlags > 0._r8) then
         finundated_lag(c) = finundated_lag(c) * exp(-dtime/redoxlags) &
                           + finundated(c) * (1._r8 - exp(-dtime/redoxlags))
      else
         finundated_lag(c) = finundated(c)
      end if

   end do

   do fc = 1, num_lakec
      c = filter_lakec(fc)

      totcolch4_bef(c) = totcolch4(c)
      totcolch4(c) = 0._r8
   end do
 
   ! Check to see if finundated changed since the last timestep.  If it increased, then reduce conc_ch4_sat
   ! proportionally.  If it decreased, then add flux to atm.

   do j=1,nlevsoi
      do fc = 1, num_soilc
         c = filter_soilc(fc)

         if (j==1) then
            ch4_dfsat_flux(c) = 0._r8
         end if

         if (fsat_bef(c) /= spval .and. finundated(c) > fsat_bef(c)) then !Reduce conc_ch4_sat
            dfsat = finundated(c) - fsat_bef(c)
            conc_ch4_sat(c,j) = (fsat_bef(c)*conc_ch4_sat(c,j) + dfsat*conc_ch4_unsat(c,j)) / finundated(c)
         else if (fsat_bef(c) /= spval .and. finundated(c) < fsat_bef(c)) then
            ch4_dfsat_flux(c) = ch4_dfsat_flux(c) + (fsat_bef(c) - finundated(c))*(conc_ch4_sat(c,j) - conc_ch4_unsat(c,j)) * &
                                dz(c,j) / dtime * catomw / 1000._r8 ! mol --> kg
         end if
      end do
   end do
   
   !!!! Begin biochemistry

   ! First for soil
   lake = .false.

   ! Do CH4 Annual Averages
   call ch4annualupdate(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, num_soilp, filter_soilp)

   if (nlevdecomp == 1) then
      ! Set rootfr = spval for non-veg points, unless pwtc > 0.99, in which case set it equal to uniform dist.
      do j=1, nlevsoi
         do fp = 1, num_soilp
            p = filter_soilp(fp)
            c = pcolumn(p)
   
            if (ivt(p) /= noveg) then
               rootfraction(p,j) = rootfr(p,j)
            else if (pwtc(p) < 0.99_r8) then
               rootfraction(p,j) = spval
            else
               rootfraction(p,j) = dz(c,j) / zi(c,nlevsoi)   ! Set equal to uniform distribution
            end if
         end do
      end do
   end if

   ! Call pft2col for rootfr & grnd_ch4_cond.
   if (nlevdecomp == 1) call p2c (lbp, ubp, lbc, ubc, nlevgrnd, rootfraction, rootfr_col, 'unity')
   ! Needed to use non-filter form above so that spval would be treated properly.
   call p2c (num_soilc, filter_soilc, grnd_ch4_cond, grnd_ch4_cond_col)

   if (nlevdecomp == 1) then
      ! Check for inactive columns
      do j=1, nlevsoi
         do fc = 1, num_soilc
            c = filter_soilc(fc)
   
            if (.not. cactive(c)) rootfr_col(c,j) = dz(c,j) / zi(c,nlevsoi)
         end do
      end do
   end if

   do fc = 1, num_soilc
      c = filter_soilc(fc)
      g = cgridcell(c)
      ! Set the atmospheric CH4 and O2 concentrations
      c_atm(g,1) =  forc_pch4(g) / rgasm / forc_t(g) ! [mol/m3 air]
      c_atm(g,2) =  forc_po2(g) / rgasm / forc_t(g) ! [mol/m3 air]
      !c_atm(g,3) =  forc_pco2(g) / rgasm / forc_t(g) ! [mol/m3 air]
      ! Not currently used
   enddo

   do sat = 0, 1 ! 0 == unsaturated; 1 = saturated

      ! Get index of water table
      if (sat == 0) then ! unsaturated
         call get_jwt (lbc, ubc, num_soilc, filter_soilc, jwt)
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            zwt_ch4_unsat(c) = zi(c,jwt(c))

         end do

         ! Update lagged saturation status of layer
         do j=1,nlevsoi
            do fc = 1, num_soilc
               c = filter_soilc(fc)

               if (j > jwt(c) .and. redoxlags_vertical > 0._r8) then ! saturated currently
                  layer_sat_lag(c,j) = layer_sat_lag(c,j) * exp(-dtime/redoxlags_vertical) &
                                     + (1._r8 - exp(-dtime/redoxlags_vertical))
               else if (redoxlags_vertical > 0._r8) then
                  layer_sat_lag(c,j) = layer_sat_lag(c,j) * exp(-dtime/redoxlags_vertical)
               else if (j > jwt(c)) then  ! redoxlags_vertical = 0
                  layer_sat_lag(c,j) = 1._r8
               else
                  layer_sat_lag(c,j) = 0._r8
               end if
            end do
         end do
            
      else ! saturated
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            jwt(c) = 0
         end do
      endif

      ! calculate CH4 production in each soil layer
      call ch4_prod (lbc, ubc, lbp, ubp, num_soilc, filter_soilc, num_soilp, filter_soilp, jwt, sat, lake)

      ! calculate CH4 oxidation in each soil layer
      call ch4_oxid (lbc, ubc, num_soilc, filter_soilc, jwt, sat, lake)

      ! calculate CH4 aerenchyma losses in each soil layer
      call ch4_aere (lbc, ubc, lbp, ubp, num_soilc, filter_soilc, num_soilp, filter_soilp, jwt, sat, lake)

      ! calculate CH4 ebullition losses in each soil layer
      call ch4_ebul (lbc, ubc, num_soilc, filter_soilc, jwt, sat, lake)

      ! Solve CH4 reaction/diffusion equation 
      call ch4_tran (lbc, ubc, num_soilc, filter_soilc, jwt, dtime_ch4, sat, lake)
      ! Competition for oxygen will occur here.

   enddo ! sat/unsat

   ! Now do over lakes
   if (allowlakeprod) then
      lake = .true.
      sat = 1
      do fc = 1, num_lakec
         c = filter_lakec(fc)
         jwt(c) = 0
      end do

      ! calculate CH4 production in each lake layer
      call ch4_prod (lbc, ubc, lbp, ubp, num_lakec, filter_lakec, 0, dummyfilter, jwt, sat, lake)

      ! calculate CH4 oxidation in each lake layer
      call ch4_oxid (lbc, ubc, num_lakec, filter_lakec, jwt, sat, lake)

      ! calculate CH4 aerenchyma losses in each lake layer
      call ch4_aere (lbc, ubc, lbp, ubp, num_lakec, filter_lakec, 0, dummyfilter, jwt, sat, lake)
      ! The p filter will not be used here; the relevant column vars will just be set to 0.

      ! calculate CH4 ebullition losses in each lake layer
      call ch4_ebul (lbc, ubc, num_lakec, filter_lakec, jwt, sat, lake)

      ! Solve CH4 reaction/diffusion equation 
      call ch4_tran (lbc, ubc, num_lakec, filter_lakec, jwt, dtime_ch4, sat, lake)
      ! Competition for oxygen will occur here.

   end if

   ! Average up to gridcell flux and column oxidation and production rate.
   ! First weight the soil columns by finundated.
   do j=1,nlevsoi
      do fc = 1, num_soilc
         c = filter_soilc(fc)

         if (j == 1) then
            totalsat = ch4_surf_diff_sat(c) + ch4_surf_aere_sat(c) + ch4_surf_ebul_sat(c)
            totalunsat = ch4_surf_diff_unsat(c) + ch4_surf_aere_unsat(c) + ch4_surf_ebul_unsat(c)
            ch4_surf_flux_tot(c) = (finundated(c)*totalsat + (1._r8 - finundated(c))*totalunsat) * &
                                catomw / 1000._r8
                                !Convert from mol to kg C
             ! ch4_oxid_tot and ch4_prod_tot are initialized to zero above
         end if

         ch4_oxid_tot(c) = ch4_oxid_tot(c) + (finundated(c)*ch4_oxid_depth_sat(c,j) + &
                               (1._r8 - finundated(c))*ch4_oxid_depth_unsat(c,j))*dz(c,j) * catomw
                             !Convert from mol to g C
         ch4_prod_tot(c) = ch4_prod_tot(c) + (finundated(c)*ch4_prod_depth_sat(c,j) + &
                            (1._r8 - finundated(c))*ch4_prod_depth_unsat(c,j))*dz(c,j) * catomw
                             !Convert from mol to g C
         if (j == nlevsoi) then
            ! Adjustment to NEE flux to atm. for methane production
            nem_col(c) = nem_col(c) - ch4_prod_tot(c)
            ! Adjustment to NEE flux to atm. for methane oxidation
            nem_col(c) = nem_col(c) + ch4_oxid_tot(c)
         end if
      end do
   end do

   !Correct for discrepancies in CH4 concentration from changing finundated
   do fc = 1, num_soilc
      c = filter_soilc(fc)

      if (fsat_bef(c) /= spval) then ! not first timestep
         ch4_surf_flux_tot(c) = ch4_surf_flux_tot(c) + ch4_dfsat_flux(c)
      end if
      fsat_bef(c) = finundated(c)
   end do

   if (allowlakeprod) then
      do j=1,nlevsoi
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            if (j == 1) then
                ! ch4_oxid_tot and ch4_prod_tot are initialized to zero above
               totalsat = ch4_surf_diff_sat(c) + ch4_surf_aere_sat(c) + ch4_surf_ebul_sat(c)
               ch4_surf_flux_tot(c) = totalsat*catomw / 1000._r8
            end if

            ch4_oxid_tot(c) = ch4_oxid_tot(c) + ch4_oxid_depth_sat(c,j)*dz(c,j)*catomw
            ch4_prod_tot(c) = ch4_prod_tot(c) + ch4_prod_depth_sat(c,j)*dz(c,j)*catomw

            if (.not. replenishlakec) then
               !Adjust lake_soilc for production.
               lake_soilc(c,j) = lake_soilc(c,j) - 2._r8*ch4_prod_depth_sat(c,j)*dtime*catomw
               ! Factor of 2 is for CO2 that comes off with CH4 because of stoichiometry
            end if

            if (j == nlevsoi) then
               ! Adjustment to NEE flux to atm. for methane production
               if (.not. replenishlakec) then
                  nem_col(c) = nem_col(c) + ch4_prod_tot(c)
                  ! Here this is positive because it is actually the CO2 that comes off with the methane
                  ! NOTE THIS MODE ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
                  ! IN THIS MODE.
               else ! replenishlakec
                  nem_col(c) = nem_col(c) - ch4_prod_tot(c)
                  ! Keep total C constant, just shift from CO2 to methane
               end if

               ! Adjustment to NEE flux to atm. for methane oxidation
               nem_col(c) = nem_col(c) + ch4_oxid_tot(c)

            end if


            !Set lake diagnostic output variables
            ch4_prod_depth_lake(c,j) = ch4_prod_depth_sat(c,j)
            conc_ch4_lake(c,j)       = conc_ch4_sat(c,j)
            conc_o2_lake(c,j)        = conc_o2_sat(c,j)
            ch4_oxid_depth_lake(c,j) = ch4_oxid_depth_sat(c,j)
            if (j == 1) then
               ch4_surf_diff_lake(c) = ch4_surf_diff_sat(c)
               ch4_surf_ebul_lake(c) = ch4_surf_ebul_sat(c)
            end if

         end do
      end do
   end if  ! ch4_surf_flux_tot, ch4_oxid_tot, and ch4_prod_tot should be initialized to 0 above if .not. allowlakeprod

   ! Finalize CH4 balance and check for errors

   do j=1,nlevsoi
      do fc = 1, num_soilc
         c = filter_soilc(fc)

         totcolch4(c) = totcolch4(c) + &
                (finundated(c)*conc_ch4_sat(c,j) + (1._r8-finundated(c))*conc_ch4_unsat(c,j))*dz(c,j)*catomw
                                                                                     ! mol CH4 --> g C

         if (j == nlevsoi .and. totcolch4_bef(c) /= spval) then ! not first timestep
         ! Check balance
            errch4 = totcolch4(c) - totcolch4_bef(c) - dtime*(ch4_prod_tot(c) - ch4_oxid_tot(c) &
                                  - ch4_surf_flux_tot(c)*1000._r8) ! kg C --> g C
            if (abs(errch4) > 1.e-7_r8) then ! g C / m^2 / timestep
               write(iulog,*)'CH4 Conservation Error in CH4Mod driver, nstep, c, errch4 (gC /m^2.timestep)', &
                             nstep,c,errch4
               g = cgridcell(c)
               write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
               call endrun( trim(subname)//' ERROR: Methane conservation error')
            end if 
         end if

      end do
      if (allowlakeprod) then
         do fc = 1, num_lakec
            c = filter_lakec(fc)
   
            totcolch4(c) = totcolch4(c) + conc_ch4_sat(c,j)*dz(c,j)*catomw ! mol CH4 --> g C
   
            if (j == nlevsoi .and. totcolch4_bef(c) /= spval) then ! not first timestep
            ! Check balance
               errch4 = totcolch4(c) - totcolch4_bef(c) - dtime*(ch4_prod_tot(c) - ch4_oxid_tot(c) &
                                     - ch4_surf_flux_tot(c)*1000._r8) ! kg C --> g C
               if (abs(errch4) > 1.e-7_r8) then ! g C / m^2 / timestep
                  write(iulog,*)'CH4 Conservation Error in CH4Mod driver for lake column, nstep, c, errch4 (gC/m^2.timestep)', &
                                nstep,c,errch4
                  g = cgridcell(c)
                  write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
                  call endrun( trim(subname)//' ERROR: Methane conservation error, allowlakeprod')
               end if
            end if
   
         end do
      end if
   end do


   ! Now average up to gridcell for fluxes
   call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                      &
              ch4_surf_flux_tot(lbc:ubc), flux_ch4(lbg:ubg), &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )

   call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                    &
              ch4_oxid_tot(lbc:ubc), ch4co2f(lbg:ubg),        &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )

   call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                    &
              ch4_prod_tot(lbc:ubc), ch4prodg(lbg:ubg),       &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )

   call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                    &
              nem_col(lbc:ubc), nem(lbg:ubg),       &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )


end subroutine ch4
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_prod
!
! !INTERFACE:
subroutine ch4_prod (lbc, ubc, lbp, ubp, num_methc, filter_methc, num_methp, &
                     filter_methp, jwt, sat, lake)
!
! !DESCRIPTION:
! Production is done below the water table, based on CN heterotrophic respiration.
! O2 is consumed by roots & by heterotrophic aerobes.
! Production is done separately for sat & unsat, and is adjusted for temperature, seasonal inundation,
! pH (optional), & redox lag factor.

! !USES:
   use clmtype
   use ch4varcon,        only : q10ch4base, q10ch4, rootlitfrac, cnscalefactor, mino2lim, & 
                                f_ch4, lake_decomp_fact, usephfact, anoxicmicrosites, ch4rmcnlim
   use clm_varctl,       only : anoxia
   use clm_varpar, only : nlevdecomp
#ifdef CENTURY_DECOMP
   use CNDecompCascadeMod_CENTURY, only : nlev_soildecomp_standard
#else
   use CNDecompCascadeMod_BGC, only : nlev_soildecomp_standard
#endif
   use pftvarcon, only  : noveg
   use nanmod,    only  : nan
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: lbp, ubp        ! pft-index bounds
   integer, intent(in) :: num_methc       ! number of column soil points in column filter
   integer, intent(in) :: filter_methc(ubc-lbc+1) ! column filter for soil points
   integer, intent(in) :: num_methp       ! number of soil points in pft filter
   integer, intent(in) :: filter_methp(ubp-lbp+1) ! pft filter for soil points
   integer, intent(in) :: jwt(lbc:ubc)    ! index of the soil layer right above the water table (-)
   integer, intent(in) :: sat             ! 0 = unsaturated; 1 = saturated
   logical, intent(in) :: lake            ! function called with lake filter
!
! !CALLED FROM:
! subroutine ch4()
!
! !REVISION HISTORY:
! 9/16/08: Created by William J. Riley
! 8/27/09: Modified for CLM 4 by Zack Subin
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   integer , pointer :: pcolumn(:)        ! index into column level quantities
   real(r8), pointer :: wtcol(:)          ! weight (relative to column) 
   integer , pointer :: ivt(:)            ! pft vegetation type
   real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity)
   real(r8), pointer :: rootfr_col(:,:)   ! fraction of roots in each soil layer  (nlevsoi)
   real(r8), pointer :: rootfr(:,:)       ! fraction of roots in each soil layer  (nlevsoi)

   real(r8), pointer :: dz(:,:)           ! layer thickness (m)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: z(:,:)            ! layer depth (m) (-nlevsno+1:nlevsoi) 
   real(r8), pointer :: zi(:,:)           ! interface level below a "z" level (m)
   real(r8), pointer :: somhr(:)          ! (gC/m2/s) soil organic matter heterotrophic respiration
   real(r8), pointer :: conc_o2(:,:)      ! O2 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: latdeg(:)         ! latitude (degrees)
   integer , pointer :: cgridcell(:)      ! gridcell of corresponding column
   real(r8), pointer :: annavg_finrw(:)   ! respiration-weighted annual average of finundated
   real(r8), pointer :: finundated(:)     !fractional inundated area in soil column
   !real(r8), pointer :: col_rr(:)        ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: lake_soilc(:,:)   ! total soil organic matter found in level (g C / m^3) (nlevsoi)
                                          ! needed for lakes not using CN. initialized to total organic
                                          ! content read in from soil properties
   real(r8), pointer :: pH(:)             ! soil water pH
   real(r8), pointer :: lithr(:)          ! (gC/m2/s) litter heterotrophic respiration 
   real(r8), pointer :: finundated_lag(:) ! time-lagged fractional inundated area
   real(r8), pointer :: hr_vr(:,:)        ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: o_scalar(:,:)     ! fraction by which decomposition is limited by anoxia
   real(r8), pointer :: fphr(:,:)         ! fraction of potential heterotrophic respiration
                                          ! excluding moisture and N limitation
   real(r8), pointer :: layer_sat_lag(:,:)! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
   real(r8), pointer :: pot_f_nit_vr(:,:) ! (gN/m3/s) potential soil nitrification flux


! local pointers to implicit in/out arrays
   real(r8), pointer :: ch4_prod_depth(:,:)    ! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: o2_decomp_depth(:,:)   ! O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: co2_decomp_depth(:,:)  ! CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: sif(:)                 ! (unitless) ratio applied to sat. prod. to account for seasonal inundation

!
! !OTHER LOCAL VARIABLES:
   integer :: p,c,j,g            ! indices
   integer :: fc                 ! column index
   integer :: fp                 ! PFT index
   real(r8) :: dtime
   real(r8):: base_decomp        ! base rate (mol/m2/s)
   real(r8) :: q10lake           ! For now, take to be the same as q10ch4 * 1.5.
   real(r8), parameter :: q10lakebase = 298._r8 ! (K) base temperature for lake CH4 production
   real(r8) :: partition_z

   ! added by Lei Meng to account for pH influence of CH4 production 
   real(r8), parameter :: pHmax = 9_r8
   real(r8), parameter :: pHmin = 2.2_r8    !
   real(r8), parameter :: pHopt = 6.5_r8    !optimal pH for methane production
   real(r8)            :: pH_fact_ch4       ! pH factor in methane production

   ! Factors for methanogen temperature dependence being greater than soil aerobes
   real(r8)            :: f_ch4_adj             ! Adjusted f_ch4
   real(r8)            :: t_fact_ch4            ! Temperature factor calculated using additional Q10
   ! O2 limitation on decomposition and methanogenesis
   real(r8)            :: seasonalfin           ! finundated in excess of respiration-weighted annual average
   real(r8)            :: oxinhib = 400._r8     ! inhibition of methane production by oxygen (m^3/mol)

   ! For calculating column average (rootfrac(p,j)*rr(p,j))
   real(r8)            :: rr_vr(lbc:ubc, 1:nlevsoi) ! vertically resolved column-mean root respiration (g C/m^2/s)
    

   character(len=32) :: subname='ch4_prod' ! subroutine name

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   ! Gridcell level
   latdeg          => clm3%g%latdeg

   ! Column level
   cgridcell       => clm3%g%l%c%gridcell
   t_soisno        => clm3%g%l%c%ces%t_soisno
   h2osoi_vol      => clm3%g%l%c%cws%h2osoi_vol
   watsat          => clm3%g%l%c%cps%watsat
   rootfr_col      => clm3%g%l%c%cps%pps_a%rootfr
   dz              => clm3%g%l%c%cps%dz
   z               => clm3%g%l%c%cps%z
   zi              => clm3%g%l%c%cps%zi
   somhr            => clm3%g%l%c%ccf%somhr
   lake_soilc       => clm3%g%l%c%cch4%lake_soilc
   fphr             => clm3%g%l%c%cch4%fphr
   annavg_finrw     => clm3%g%l%c%cch4%annavg_finrw
   sif              => clm3%g%l%c%cch4%sif
   finundated       => clm3%g%l%c%cws%finundated
   !col_rr           => clm3%g%l%c%ccf%pcf_a%rr
   pH               => clm3%g%l%c%cps%pH
   lithr            => clm3%g%l%c%ccf%lithr
   finundated_lag   => clm3%g%l%c%cch4%finundated_lag
   hr_vr            => clm3%g%l%c%ccf%hr_vr
   o_scalar         => clm3%g%l%c%ccf%o_scalar
   layer_sat_lag    => clm3%g%l%c%cch4%layer_sat_lag
#if (defined NITRIF_DENITRIF)
   pot_f_nit_vr     => clm3%g%l%c%cnf%pot_f_nit_vr
#endif

   ! PFT level
   pcolumn          => clm3%g%l%c%p%column
   wtcol            => clm3%g%l%c%p%wtcol
   ivt              => clm3%g%l%c%p%itype
   rootfr           => clm3%g%l%c%p%pps%rootfr
   rr               => clm3%g%l%c%p%pcf%rr


   if (sat == 0) then ! unsaturated
      ch4_prod_depth  => clm3%g%l%c%cch4%ch4_prod_depth_unsat
      o2_decomp_depth => clm3%g%l%c%cch4%o2_decomp_depth_unsat
      co2_decomp_depth=> clm3%g%l%c%cch4%co2_decomp_depth_unsat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_unsat
   else ! saturated
      ch4_prod_depth  => clm3%g%l%c%cch4%ch4_prod_depth_sat
      o2_decomp_depth => clm3%g%l%c%cch4%o2_decomp_depth_sat
      co2_decomp_depth=> clm3%g%l%c%cch4%co2_decomp_depth_sat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_sat
   endif

   dtime = get_step_size()
   q10lake = q10ch4 * 1.5_r8

   ! PFT loop to calculate vertically resolved column-averaged root respiration
   if (.not. lake) then
      rr_vr(:,:) = nan

      do fp = 1, num_methc
         c = filter_methc(fp)
         rr_vr(c,:) = 0.0_r8
      end do
      do j=1,nlevsoi
         do fp = 1, num_methp
            p = filter_methp(fp)
            c = pcolumn(p)

            if (wtcol(p) > 0._r8 .and. ivt(p) /= noveg) then
               rr_vr(c,j) = rr_vr(c,j) + rr(p)*rootfr(p,j)*wtcol(p)
            end if
         end do
      end do
   end if

   partition_z = 1._r8
   base_decomp = 0.0_r8

   ! column loop to partition decomposition_rate into each soil layer
   do j=1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)
         g = cgridcell(c)
         
         if (.not. lake) then
#if (defined CN)
            ! Use soil heterotrophic respiration (based on Wania)
            base_decomp = (somhr(c)+lithr(c)) / catomw
            ! Convert from gC to molC
            ! Multiply base_decomp by factor accounting for lower carbon stock in seasonally inundated areas than
            ! if it were inundated all year.
            ! This is to reduce emissions in seasonally inundated zones, because the eq.
            ! C-flux will be less than predicted by a non-O2-lim model
            if (sat == 1) then
               sif(c) = 1._r8
               if (.not. anoxia) then
                  if (annavg_finrw(c) /= spval) then
                     seasonalfin = max(finundated(c)-annavg_finrw(c), 0._r8)
                     if (seasonalfin > 0._r8) then
                        sif(c) = (annavg_finrw(c) + mino2lim*seasonalfin) / finundated(c)
                        base_decomp = base_decomp * sif(c)
                     end if
                  end if
               end if ! anoxia
            end if
#else
            call endrun( trim(subname)//' ERROR: No source for decomp rate in CH4Prod. CH4 model currently requires CN.' )
#endif
! defined CN
         ! For sensitivity studies
            base_decomp = base_decomp * cnscalefactor

         else !lake
            base_decomp = lake_decomp_fact * lake_soilc(c,j) * dz(c,j) * &
                                     q10lake**( (t_soisno(c,j)-q10lakebase)/10._r8) / catomw
                                     ! convert from g C to mol C
         end if

         ! For all landunits, prevent production or oxygen consumption when soil is at or below freezing.
         ! If using VERTSOILC, it is OK to use base_decomp as given because liquid water stress will limit decomp.
         if (t_soisno(c,j) <= tfrz .and. (nlevdecomp == 1 .or. lake)) base_decomp = 0._r8

         ! depth dependence of production either from rootfr or decomp model
         if (.not. lake) then ! use default rootfr, averaged to the column level in the ch4 driver, or vert HR
            if (nlevdecomp == 1) then ! not VERTSOILC
               if (j <= nlev_soildecomp_standard) then  ! Top 5 levels are also used in the CLM code for establishing temperature
                                 ! and moisture constraints on SOM activity
                  partition_z = rootfr_col(c,j)*rootlitfrac + (1._r8 - rootlitfrac)*dz(c,j)/zi(c,nlev_soildecomp_standard)
               else
                  partition_z = rootfr_col(c,j)*rootlitfrac
               end if
            else
               if ( (somhr(c) + lithr(c)) > 0._r8) then
                  partition_z = hr_vr(c,j) * dz(c,j) / (somhr(c) + lithr(c))
               else
                  partition_z = 1._r8
               end if
            end if
         else ! lake
            partition_z = 1._r8
         endif

         ! Adjust f_ch4 to account for the fact that methanogens may have a higher Q10 than aerobic decomposers.
         ! Note this is crude and should ideally be applied to all anaerobic decomposition rather than just the
         ! f_ch4.
         f_ch4_adj = 1.0_r8
         if (.not. lake) then
            t_fact_ch4 = q10ch4**((t_soisno(c,j) - q10ch4base)/10._r8)
            ! Adjust f_ch4 by the ratio
            f_ch4_adj = f_ch4 * t_fact_ch4

            ! Remove CN nitrogen limitation, as methanogenesis is not N limited.
            ! Also remove (low) moisture limitation
            if (ch4rmcnlim) then
               if (j > nlevdecomp) then
                  if (fphr(c,1) > 0._r8) then
                     f_ch4_adj = f_ch4_adj / fphr(c,1)
                  end if
               else ! j == 1 or VERTSOILC
                  if (fphr(c,j) > 0._r8) then
                     f_ch4_adj = f_ch4_adj / fphr(c,j)
                  end if
               end if
            end if

         else ! lake
            f_ch4_adj = 0.5_r8 ! For lakes assume no redox limitation. Production only depends on temp, soil C, and
                               ! lifetime parameter.
         end if

         ! If switched on, use pH factor for production based on spatial pH data defined in surface data.
         if (.not. lake .and. usephfact .and. pH(c).gt. pHmin .and.pH(c).lt. pHmax) then
            pH_fact_ch4 = 10._r8**(-0.2235_r8*pH(c)*pH(c) + 2.7727_r8*pH(c) - 8.6_r8)
            ! fitted function using data from Dunfield et al. 1993  
            ! Strictly less than one, with optimum at 6.5
            ! From Lei Meng
            f_ch4_adj = f_ch4_adj * pH_fact_ch4
         else
            ! if no data, then no pH effects
         end if

         ! Redox factor
         if ( (.not. lake) .and. sat == 1 .and. finundated_lag(c) < finundated(c)) then
            f_ch4_adj = f_ch4_adj * finundated_lag(c) / finundated(c)
         else if (sat == 0 .and. j > jwt(c)) then ! Assume lag in decay of alternative electron acceptors vertically
            f_ch4_adj = f_ch4_adj * layer_sat_lag(c,j)
         end if
         ! Alternative electron acceptors will be consumed first after soil is inundated.

         f_ch4_adj = min(f_ch4_adj, 0.5_r8)
         ! Must be less than 0.5 because otherwise the actual implied aerobic respiration would be negative.
         ! The total of aer. respiration + methanogenesis must remain equal to the SOMHR calculated in CN,
         ! so that the NEE is sensible. Even perfectly anaerobic conditions with no alternative
         ! electron acceptors would predict no more than 0.5 b/c some oxygen is present in organic matter.
         ! e.g. 2CH2O --> CH4 + CO2.


       ! Decomposition uses 1 mol O2 per mol CO2 produced (happens below WT also, to deplete O2 below WT)
       ! o2_decomp_depth is the demand in the absense of O2 supply limitation, in addition to autotrophic respiration.
                     ! Competition will be done in ch4_oxid

         o2_decomp_depth(c,j) = base_decomp * partition_z / dz (c,j)
         if (anoxia) then
         ! Divide off o_scalar to use potential O2-unlimited HR to represent aerobe demand for oxygen competition
            if (.not. lake .and. j > nlevdecomp) then
               if (o_scalar(c,1) > 0._r8) then
                  o2_decomp_depth(c,j) = o2_decomp_depth(c,j) / o_scalar(c,1)
               end if
            else if (.not. lake) then ! j == 1 or VERTSOILC
               if (o_scalar(c,j) > 0._r8) then
                  o2_decomp_depth(c,j) = o2_decomp_depth(c,j) / o_scalar(c,j)
               end if
            end if
         end if ! anoxia

         ! Add root respiration
         if (.not. lake) then
            !o2_decomp_depth(c,j) = o2_decomp_depth(c,j) + col_rr(c)*rootfr(c,j)/catomw/dz(c,j) ! mol/m^3/s
            o2_decomp_depth(c,j) = o2_decomp_depth(c,j) + rr_vr(c,j)/catomw/dz(c,j) ! mol/m^3/s
                                                       ! g C/m2/s ! gC/mol O2 ! m
         end if

         ! Add oxygen demand for nitrification
#if (defined NITRIF_DENITRIF)
         if (.not. lake) then
            o2_decomp_depth(c,j) = o2_decomp_depth(c,j) + pot_f_nit_vr(c,j) * 2.0_r8/14.0_r8
                                                        ! g N/m^3/s           mol O2 / g N
         end if
#endif

         if (j .gt. jwt(c)) then ! Below the water table so anaerobic CH4 production can occur
            ! partition decomposition to layer
            ! turn into per volume-total by dz
            ch4_prod_depth(c,j) = f_ch4_adj * base_decomp * partition_z / dz (c,j)! [mol/m3-total/s]
         else ! Above the WT
            if (anoxicmicrosites) then
               ch4_prod_depth(c,j) = f_ch4_adj * base_decomp * partition_z / dz (c,j) &
                                   / (1._r8 + oxinhib*conc_o2(c,j))
            else
               ch4_prod_depth(c,j) = 0._r8 ! [mol/m3 total/s]
            endif ! anoxicmicrosites
         endif ! WT

      end do ! fc
   end do ! nlevsoi

end subroutine ch4_prod

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_oxid
!
! !INTERFACE:
subroutine ch4_oxid (lbc, ubc, num_methc, filter_methc, jwt, sat, lake)
!
! !DESCRIPTION:
! Oxidation is based on double Michaelis-Mentin kinetics, and is adjusted for low soil moisture.
! Oxidation will be limited by available oxygen in ch4_tran.

! !USES:
   use clmtype
   use clm_time_manager, only : get_step_size
   use ch4varcon, only : vmax_ch4_oxid, k_m, k_m_o2, q10_ch4oxid, smp_crit, &
                         k_m_unsat, vmax_oxid_unsat

!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: num_methc       ! number of column soil points in column filter
   integer, intent(in) :: filter_methc(ubc-lbc+1) ! column filter for soil points
   integer, intent(in)  :: jwt(lbc:ubc)   ! index of the soil layer right above the water table (-)
   integer, intent(in)  :: sat            ! 0 = unsaturated; 1 = saturated
   logical, intent(in) :: lake            ! function called with lake filter
!
! !CALLED FROM:
! subroutine ch4()
!
! !REVISION HISTORY:
! 9/16/08: Created by William J. Riley
! 8/27/09: Modified by Zack Subin for CLM 4
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   real(r8), pointer :: h2osoi_vol(:,:)     ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   real(r8), pointer :: watsat(:,:)         ! volumetric soil water at saturation (porosity)
   real(r8), pointer :: t_soisno(:,:)       ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: smp_l(:,:)          ! soil matrix potential [mm]
   real(r8), pointer :: conc_o2(:,:)        ! O2 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_ch4(:,:)       ! CH4 conc in each soil layer (mol/m3) (nlevsoi)

! local pointers to implicit in/out arrays
   real(r8), pointer :: ch4_oxid_depth(:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_oxid_depth(:,:)  ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: co2_oxid_depth(:,:) ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_decomp_depth(:,:)! O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,j             ! indices
   integer :: fc              ! column index
   real(r8) :: dtime          ! land model time step (sec)
   real(r8):: t0              ! Base temperature for Q10
   real(r8):: porevol         ! air-filled volume ratio to total soil volume
   real(r8):: h2osoi_vol_min  ! h2osoi_vol restricted to be below watsat
   real(r8):: conc_ch4_rel    ! concentration with respect to water volume (mol/m^3 water)
   real(r8):: conc_o2_rel     ! concentration with respect to water volume (mol/m^3 water)
   real(r8):: oxid_a          ! Oxidation predicted by method A (temperature & enzyme limited) (mol CH4/m3/s)
   real(r8):: smp_fact        ! factor for reduction based on soil moisture (unitless)
   real(r8):: porewatfrac     ! fraction of soil pore space that is filled with water
   real(r8):: k_h_cc, k_h_inv ! see functions below for description
   real(r8):: k_m_eff         ! effective k_m
   real(r8):: vmax_eff        ! effective vmax
   
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   
   ! Column level
   h2osoi_vol     => clm3%g%l%c%cws%h2osoi_vol
   watsat         => clm3%g%l%c%cps%watsat
   t_soisno       => clm3%g%l%c%ces%t_soisno
   smp_l             => clm3%g%l%c%cws%smp_l

   if (sat == 0) then ! unsaturated
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_unsat
      o2_oxid_depth  => clm3%g%l%c%cch4%o2_oxid_depth_unsat
      co2_oxid_depth => clm3%g%l%c%cch4%co2_oxid_depth_unsat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_unsat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_unsat
      o2_decomp_depth => clm3%g%l%c%cch4%o2_decomp_depth_unsat
   else ! saturated
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_sat
      o2_oxid_depth  => clm3%g%l%c%cch4%o2_oxid_depth_sat
      co2_oxid_depth => clm3%g%l%c%cch4%co2_oxid_depth_sat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_sat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_sat
      o2_decomp_depth => clm3%g%l%c%cch4%o2_decomp_depth_sat
   endif

   ! Get land model time step
   dtime = get_step_size()

   t0 = tfrz + 12._r8 ! Walter, for Michigan site where the 45 M/h comes from

   ! Loop to determine oxidation in each layer
   do j=1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc(fc)

         if (sat == 1 .or. j > jwt(c)) then
         ! Literature (e.g. Bender & Conrad, 1992) suggests lower k_m and vmax for high-CH4-affinity methanotrophs in
         ! upland soils consuming ambient methane.
            k_m_eff = k_m
            vmax_eff = vmax_ch4_oxid
         else
            k_m_eff = k_m_unsat
            vmax_eff = vmax_oxid_unsat
         end if

         porevol = max(watsat(c,j) - h2osoi_vol(c,j), 0._r8)
         h2osoi_vol_min = min(watsat(c,j), h2osoi_vol(c,j))
         if (j <= jwt(c) .and. smp_l(c,j) < 0._r8) then
            smp_fact = exp(-smp_l(c,j)/smp_crit)
            ! Schnell & King, 1996, Figure 3
         else
            smp_fact = 1._r8
         end if

         if (j .le. jwt(c)) then ! Above the water table
            k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1)))
            k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
            conc_ch4_rel = conc_ch4(c,j) / (h2osoi_vol_min + porevol/k_h_cc)

            k_h_inv = exp(-c_h_inv(2) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(2)))
            k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
            conc_o2_rel  = conc_o2(c,j) / (h2osoi_vol_min + porevol/k_h_cc)
         else
            conc_ch4_rel = conc_ch4(c,j) / watsat(c,j)
            conc_o2_rel  = conc_o2(c,j) / watsat(c,j)
         endif

         oxid_a              = vmax_eff     * h2osoi_vol_min* conc_ch4_rel / (k_m_eff + conc_ch4_rel) &
         ![mol/m3-t/s]         [mol/m3-w/s]    [m3-w/m3-t]     [mol/m3-w]    [mol/m3-w]  [mol/m3-w]
                                                             * conc_o2_rel / (k_m_o2 + conc_o2_rel) &
                               * q10_ch4oxid ** ((t_soisno(c,j) - t0) / 10._r8) * smp_fact
        
         ! For all landunits / levels, prevent oxidation if at or below freezing
         if (t_soisno(c,j) <= tfrz) oxid_a = 0._r8

         ch4_oxid_depth(c,j) = oxid_a
         o2_oxid_depth(c,j) = ch4_oxid_depth(c,j) * 2._r8

      end do
   end do

end subroutine ch4_oxid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_aere
!
! !INTERFACE:
subroutine ch4_aere (lbc, ubc, lbp, ubp, num_methc, filter_methc, num_methp, &
                     filter_methp, jwt, sat, lake)
!
! !DESCRIPTION:
! Arctic c3 grass (which is often present in fens) and all vegetation in inundated areas is assumed to have
! some root porosity. Currently, root porosity is allowed to be different for grasses & non-grasses.
! CH4 diffuses out and O2 diffuses into the soil.  CH4 is also lossed via transpiration, which is both
! included in the "aere" variables and output separately.  In practice this value is small.
! By default upland veg. has small 5% porosity but this can be switched to be equal to inundated porosity.

! !USES:
   use clmtype
   use clm_varcon      , only : spval, rpi
   use clm_time_manager, only : get_step_size
   use pftvarcon,        only : nc3_arctic_grass, crop, nc3_nonarctic_grass, nc4_grass, noveg
   use ch4varcon,        only : scale_factor_aere, transpirationloss, nongrassporosratio, unsat_aere_ratio, &
                                usefrootc, aereoxid
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: lbp, ubp        ! pft-index bounds
   integer, intent(in) :: num_methc       ! number of column soil points in column filter
   integer, intent(in) :: filter_methc(ubc-lbc+1) ! column filter for soil points
   integer, intent(in) :: num_methp       ! number of soil points in pft filter
   integer, intent(in) :: filter_methp(ubp-lbp+1) ! pft filter for soil points
   integer, intent(in)  :: jwt(lbc:ubc)   ! index of the soil layer right above the water table (-)
   integer, intent(in) :: sat             ! 0 = unsaturated; 1 = saturated
   logical, intent(in) :: lake            ! function called with lake filter
!
! !CALLED FROM:
! subroutine ch4()
!
! !REVISION HISTORY:
! 9/16/08: Created by William J. Riley
! 8/27/09: Modified for CLM 4 by Zack Subin
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   integer , pointer :: ivt(:)          ! pft vegetation type
   real(r8), pointer :: wtcol(:)        ! weight (relative to column) 
   integer , pointer :: cgridcell(:)    ! gridcell of corresponding column
   integer , pointer :: pcolumn(:)      ! index into column level quantities
   real(r8), pointer :: rootfr(:,:)     ! fraction of roots in each soil layer  (nlevsoi)
   real(r8), pointer :: rootr(:,:)      ! effective fraction of roots in each soil layer  (nlevgrnd)
   real(r8), pointer :: c_atm(:,:)      ! CH4, O2, CO2 atmospheric conc  (mol/m3)
   real(r8), pointer :: elai(:)         ! one-sided leaf area index with burying by snow
   real(r8), pointer :: h2osoi_vol(:,:) ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   real(r8), pointer :: watsat(:,:)     !volumetric soil water at saturation (porosity)
   real(r8), pointer :: t_soisno(:,:)   ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: ch4_oxid_depth(:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: qflx_tran_veg(:)! vegetation transpiration (mm H2O/s) (+ = to atm)
   real(r8), pointer :: canopy_cond(:)  ! tracer conductance for canopy [m/s]
   real(r8), pointer :: conc_o2(:,:)    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_ch4(:,:)   ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: z(:,:)          ! layer depth (m) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: dz(:,:)         ! layer thickness (m)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: annsum_npp(:)   ! annual sum NPP (gC/m2/yr)
   real(r8), pointer :: annavg_agnpp(:) ! (gC/m2/s) annual average aboveground NPP
   real(r8), pointer :: annavg_bgnpp(:) ! (gC/m2/s) annual average belowground NPP
   real(r8), pointer :: ch4_prod_depth(:,:)   ! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: grnd_ch4_cond(:)! tracer conductance for boundary layer [m/s] (pft-level!)
   real(r8), pointer :: frootc(:)       ! (gC/m2) fine root C

! local pointers to implicit in/out arrays
   real(r8), pointer :: ch4_aere_depth(:,:) ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: co2_aere_depth(:,:) ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_tran_depth(:,:) ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_aere_depth(:,:)  ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
!
! !OTHER LOCAL VARIABLES:
   integer :: p,c,g,j       ! indices
   integer :: fc,fp         ! soil filter column index
   real(r8) :: f_oxid       ! fraction of CH4 oxidized in oxic zone around roots
   real(r8) :: diffus_aere  ! gas diffusivity through aerenchyma (m^2/s)
   real(r8) :: m_tiller 
   real(r8) :: n_tiller 
   real(r8) :: poros_tiller 
   real(r8) :: rob          ! root obliquity, e.g. csc of root angle relative to vertical
                            ! (ratio of root total length to depth)
   real(r8) :: area_tiller  ! cross-sectional area of tillers (m^2/m^2)
   real(r8) :: tranloss     ! loss due to transpiration (mol / m3 /s)
   real(r8):: aere, aeretran, oxaere ! (mol / m3 /s)
   real(r8):: k_h_cc, k_h_inv, dtime, oxdiffus, anpp, nppratio, h2osoi_vol_min, conc_ch4_wat
   real(r8):: aerecond      ! aerenchyma conductance (m/s)
   real(r8), parameter :: smallnumber = 1.e-12_r8
   real(r8), parameter :: porosmin = 0.05_r8 ! minimum aerenchyma porosity (unitless)
     
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   !gridcell level
   c_atm          => clm3%g%gch4%c_atm
   !column level
   cgridcell => clm3%g%l%c%gridcell
   z         => clm3%g%l%c%cps%z
   dz        => clm3%g%l%c%cps%dz
   h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
   watsat    => clm3%g%l%c%cps%watsat
   t_soisno  => clm3%g%l%c%ces%t_soisno
   !pft level
   wtcol      => clm3%g%l%c%p%wtcol
   ivt       => clm3%g%l%c%p%itype
   pcolumn   => clm3%g%l%c%p%column
   elai       => clm3%g%l%c%p%pps%elai
   frootc     => clm3%g%l%c%p%pcs%frootc
   rootfr    => clm3%g%l%c%p%pps%rootfr
   rootr     => clm3%g%l%c%p%pps%rootr
   qflx_tran_veg  => clm3%g%l%c%p%pwf%qflx_tran_veg
   canopy_cond => clm3%g%l%c%p%pps%canopy_cond
   annsum_npp            => clm3%g%l%c%p%pepv%annsum_npp
   annavg_agnpp          => clm3%g%l%c%p%pcf%annavg_agnpp
   annavg_bgnpp          => clm3%g%l%c%p%pcf%annavg_bgnpp
   grnd_ch4_cond       => clm3%g%l%c%p%pps%grnd_ch4_cond


   if (sat == 0) then ! unsaturated
      ch4_aere_depth => clm3%g%l%c%cch4%ch4_aere_depth_unsat
      ch4_tran_depth => clm3%g%l%c%cch4%ch4_tran_depth_unsat
      o2_aere_depth => clm3%g%l%c%cch4%o2_aere_depth_unsat
      co2_aere_depth => clm3%g%l%c%cch4%co2_aere_depth_unsat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_unsat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_unsat
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_unsat
      ch4_prod_depth => clm3%g%l%c%cch4%ch4_prod_depth_unsat
   else ! saturated
      ch4_aere_depth => clm3%g%l%c%cch4%ch4_aere_depth_sat
      ch4_tran_depth => clm3%g%l%c%cch4%ch4_tran_depth_sat
      o2_aere_depth => clm3%g%l%c%cch4%o2_aere_depth_sat
      co2_aere_depth => clm3%g%l%c%cch4%co2_aere_depth_sat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_sat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_sat
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_sat
      ch4_prod_depth => clm3%g%l%c%cch4%ch4_prod_depth_sat
   endif

   dtime = get_step_size()

   ! Initialize ch4_aere_depth
   do j=1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)
         ch4_aere_depth(c,j) = 0._r8
         ch4_tran_depth(c,j) = 0._r8
         o2_aere_depth(c,j) = 0._r8
      end do
   end do

   diffus_aere = d_con_g(1,1)*1.e-4_r8  ! for CH4: m^2/s
   rob = 3._r8 ! ratio of root length to vertical depth ("obliquity")
   ! This parameter is poorly constrained and should be done on a PFT-specific basis...

   ! point loop to partition aerenchyma flux into each soil layer
   if (.not. lake) then
      do j=1,nlevsoi
         do fp = 1, num_methp
            p = filter_methp (fp)
            c = pcolumn(p)
            g = cgridcell(c)

            ! Calculate transpiration loss
            if (transpirationloss .and. ivt(p) /= noveg) then !allow tloss above WT ! .and. j > jwt(c)) then
               ! Calculate water concentration
               h2osoi_vol_min = min(watsat(c,j), h2osoi_vol(c,j))
               k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1)))
               k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm
               conc_ch4_wat = conc_ch4(c,j) / ( (watsat(c,j)-h2osoi_vol_min)/k_h_cc + h2osoi_vol_min)

               tranloss = conc_ch4_wat *             rootr(p,j)*qflx_tran_veg(p) / dz(c,j) / 1000._r8
               ! mol/m3/s    mol/m3                                   mm / s         m           mm/m
               ! Use rootr here for effective per-layer transpiration, which may not be the same as rootfr
               tranloss = max(tranloss, 0._r8) ! in case transpiration is pathological
            else
               tranloss = 0._r8
            end if

            ! Calculate aerenchyma diffusion
            if (j > jwt(c) .and. t_soisno(c,j) > tfrz .and. ivt(p) /= noveg) then
            ! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
            ! the new components are in; if there is any tuning to be done to get a realistic global flux,
            ! this would probably be the place.  We will have to document clearly in the Tech Note
            ! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)

               anpp = annsum_npp(p) ! g C / m^2/yr
               anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools

               if (annavg_agnpp(p) /= spval .and. annavg_bgnpp(p) /= spval .and. &
                   annavg_agnpp(p) > 0._r8 .and. annavg_bgnpp(p) > 0._r8) then
                  nppratio = annavg_bgnpp(p) / (annavg_agnpp(p) + annavg_bgnpp(p))
               else
                  nppratio = 0.5_r8
               end if

               ! Estimate area of tillers (see Wania thesis)
               !m_tiller = anpp * r_leaf_root * lai ! (4.17 Wania)
               !m_tiller = 600._r8 * 0.5_r8 * 2._r8  ! used to be 300
               ! Note: this calculation is based on Arctic graminoids, and should be refined for woody plants, if not
               ! done on a PFT-specific basis.

               if (usefrootc) then
                  m_tiller = frootc(p) ! This will yield much smaller aere area.
               else
                  m_tiller = anpp * nppratio * elai(p)
               end if

               n_tiller = m_tiller / 0.22_r8

               if (ivt(p) == nc3_arctic_grass .or. crop(ivt(p)) == 1 .or. &
                     ivt(p) == nc3_nonarctic_grass .or. ivt(p) == nc4_grass) then
                  poros_tiller = 0.3_r8  ! Colmer 2003
               else
                  poros_tiller = 0.3_r8 * nongrassporosratio
               end if

               if (sat == 0) then
                  poros_tiller = poros_tiller * unsat_aere_ratio
               end if

               poros_tiller = max(poros_tiller, porosmin)

               area_tiller = scale_factor_aere * n_tiller * poros_tiller * rpi * 2.9e-3_r8**2._r8 ! (m2/m2)

               k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1))) ! (4.12) Wania (L atm/mol)
               k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
               aerecond = area_tiller * rootfr(p,j) * diffus_aere / (z(c,j)*rob)
               ! Add in boundary layer resistance
               aerecond = 1._r8 / (1._r8/(aerecond+smallnumber) + 1._r8/(grnd_ch4_cond(p)+smallnumber))
               
               aere = aerecond * (conc_ch4(c,j)/watsat(c,j)/k_h_cc - c_atm(g,1)) / dz(c,j) ![mol/m3-total/s]
               !ZS: Added watsat & Henry's const.
               aere = max(aere, 0._r8) ! prevent backwards diffusion

               ! Do oxygen diffusion into layer
               k_h_inv = exp(-c_h_inv(2) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(2)))
               k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
               oxdiffus = diffus_aere * d_con_g(2,1) / d_con_g(1,1) ! adjust for O2:CH4 molecular diffusion
               aerecond = area_tiller * rootfr(p,j) * oxdiffus / (z(c,j)*rob)
               aerecond = 1._r8 / (1._r8/(aerecond+smallnumber) + 1._r8/(grnd_ch4_cond(p)+smallnumber))
               oxaere = -aerecond *(conc_o2(c,j)/watsat(c,j)/k_h_cc - c_atm(g,2)) / dz(c,j) ![mol/m3-total/s]
               oxaere = max(oxaere, 0._r8)
               ! Diffusion in is positive; prevent backwards diffusion
               if (aereoxid >= 0._r8) then ! fixed aere oxid proportion; will be done in ch4_tran
                  oxaere = 0._r8
               end if
            else
               aere = 0._r8
               oxaere = 0._r8
            end if ! veg type, below water table, & above freezing

            ! Impose limitation based on available methane during timestep
            ! By imposing the limitation here, don't allow aerenchyma access to methane from other PFTs.
            aeretran = min(aere+tranloss, conc_ch4(c,j)/dtime + ch4_prod_depth(c,j))
            ch4_aere_depth (c, j) = ch4_aere_depth(c,j) + aeretran*wtcol(p) !pft weight in col.
            ch4_tran_depth (c, j) = ch4_tran_depth(c,j) + min(tranloss, aeretran)*wtcol(p)
            o2_aere_depth  (c, j) = o2_aere_depth (c,j) + oxaere*wtcol(p)
         end do ! p filter
      end do ! over levels
   end if ! not lake

end subroutine ch4_aere

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_ebul 
!
! !INTERFACE:
subroutine ch4_ebul (lbc, ubc, num_methc, filter_methc, jwt, sat, lake)
!
! !DESCRIPTION:
! Bubbling is based on temperature & pressure dependent solubility (k_h_cc), with assumed proportion of bubbles
! which are CH4, and assumed early nucleation at vgc_max sat (Wania).
! Bubbles are released to the water table surface in ch4_tran.

! !USES:
   use clmtype
   use clm_time_manager, only : get_step_size
   use ch4varcon       , only : vgc_max

!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: num_methc       ! number of column soil points in column filter
   integer, intent(in) :: filter_methc(ubc-lbc+1) ! column filter for soil points
   integer, intent(in)  :: jwt(lbc:ubc)   ! index of the soil layer right above the water table (-)
   integer, intent(in)  :: sat            ! 0 = unsaturated; 1 = saturated
   logical, intent(in) :: lake            ! function called with lake filter
!
! !CALLED FROM:
! subroutine ch4()
!
! !REVISION HISTORY:
! 9/16/08: Created by William J. Riley
! 8/27/09: Modified for CLM 4 by Zack Subin
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   integer , pointer :: cgridcell(:)      ! gridcell of corresponding column
   real(r8), pointer :: z(:,:)            ! soil layer depth (m)
   real(r8), pointer :: dz(:,:)           ! layer thickness (m)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: zi(:,:)           ! interface level below a "z" level (m)
   real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity)
   real(r8), pointer :: forc_pbot(:)      ! atmospheric pressure (Pa)
   real(r8), pointer :: lakedepth(:)      ! column lake depth (m)
   real(r8), pointer :: lake_icefrac(:,:) ! mass fraction of lake layer that is frozen
   real(r8), pointer :: ch4_aere_depth(:,:) ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_oxid_depth(:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: frac_h2osfc(:)    ! fraction of ground covered by surface water (0 to 1)
   real(r8), pointer :: h2osfc(:)         ! surface water (mm)

! local pointers to implicit in/out arrays
   real(r8), pointer :: ch4_ebul_depth(:,:) ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_ebul_total(:)   ! Total column CH4 ebullition (mol/m2/s)
   real(r8), pointer :: conc_ch4(:,:)       ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,j,g    ! indices
   integer :: fc       ! soil filter column index
   integer :: fp       ! soil filter pft index
   real(r8) :: dtime   ! land model time step (sec)
   real(r8) :: vgc     ! volumetric CH4 content (m3 CH4/m3 pore air)
   real(r8) :: vgc_min ! minimum aqueous CH4 content when ebullition ceases
   real(r8) :: k_h_inv ! 
   real(r8) :: k_h     ! 
   real(r8) :: k_h_cc  ! 
   real(r8) :: pressure! sum atmospheric and hydrostatic pressure
   real(r8) :: bubble_f! CH4 content in gas bubbles (Kellner et al. 2006)
   real(r8) :: ebul_timescale

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   ! Gridcell level
   forc_pbot => clm_a2l%forc_pbot
   ! Column level
   cgridcell => clm3%g%l%c%gridcell
   z         => clm3%g%l%c%cps%z
   dz        => clm3%g%l%c%cps%dz
   zi        => clm3%g%l%c%cps%zi
   t_soisno  => clm3%g%l%c%ces%t_soisno
   h2osoi_vol=> clm3%g%l%c%cws%h2osoi_vol
   watsat    => clm3%g%l%c%cps%watsat
   lakedepth      => clm3%g%l%c%cps%lakedepth
   lake_icefrac   => clm3%g%l%c%cws%lake_icefrac
   h2osfc         => clm3%g%l%c%cws%h2osfc
   frac_h2osfc    => clm3%g%l%c%cps%frac_h2osfc

   if (sat == 0) then ! unsaturated
      ch4_ebul_depth => clm3%g%l%c%cch4%ch4_ebul_depth_unsat
      ch4_ebul_total => clm3%g%l%c%cch4%ch4_ebul_total_unsat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_unsat
      ch4_aere_depth => clm3%g%l%c%cch4%ch4_aere_depth_unsat
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_unsat
   else ! saturated
      ch4_ebul_depth => clm3%g%l%c%cch4%ch4_ebul_depth_sat
      ch4_ebul_total => clm3%g%l%c%cch4%ch4_ebul_total_sat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_sat
      ch4_aere_depth => clm3%g%l%c%cch4%ch4_aere_depth_sat
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_sat
   endif

   ! Get land model time step
   dtime = get_step_size()

   bubble_f = 0.57_r8 ! CH4 content in gas bubbles (Kellner et al. 2006)
   vgc_min = vgc_max
   ebul_timescale = dtime ! Allow fast bubbling

   ! column loop to estimate ebullition CH4 flux from each soil layer
   do j=1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)
         g = cgridcell(c)

         if (j .gt. jwt(c) .and. t_soisno(c,j) > tfrz) then ! Ebullition occurs only below the water table

            k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1))) ! (4.12 Wania) (atm.L/mol)
            k_h = 1._r8 / k_h_inv ! (mol/L.atm)
            k_h_cc = t_soisno(c,j) * k_h * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)] 

            if (.not. lake) then
                pressure = forc_pbot(g) + denh2o * grav * (z(c,j)-zi(c,jwt(c))) ! (Pa)
                if (sat == 1 .and. frac_h2osfc(c) > 0._r8) then ! Add ponding pressure head
                   pressure = pressure + denh2o * grav * h2osfc(c)/1000._r8/frac_h2osfc(c)
                                                          ! mm     / mm/m
                end if
            else
                pressure = forc_pbot(g) + denh2o * grav * (z(c,j) + lakedepth(c))
            end if

            ! Compare partial pressure to ambient pressure.
            vgc = conc_ch4(c,j) / watsat(c,j) / k_h_cc * rgasm * t_soisno(c,j) / pressure
                  ! [mol/m3t]      [m3w/m3t]   [m3g/m3w]  [Pa/(mol/m3g)]          [Pa]

            if (vgc > vgc_max * bubble_f) then ! If greater than max value, remove amount down to vgc_min
               ch4_ebul_depth (c,j) = (vgc - vgc_min * bubble_f) * conc_ch4(c,j) / ebul_timescale
               ! [mol/m3t/s]                                       [mol/m3t]         [s]
            else
               ch4_ebul_depth (c,j) = 0._r8
            endif

         else ! above the water table or freezing
            ch4_ebul_depth (c,j) = 0._r8
         endif ! below the water table and not freezing

         ! Prevent ebullition from reaching the surface for frozen lakes
         if (lake .and. lake_icefrac(c,1) > 0.1_r8) ch4_ebul_depth(c,j) = 0._r8

      end do ! fc
   end do ! j

end subroutine ch4_ebul 

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_tran
!
! !INTERFACE:
subroutine ch4_tran (lbc, ubc, num_methc, filter_methc, jwt, dtime_ch4, sat, lake)
!
! !DESCRIPTION:
! Solves the reaction & diffusion equation for the timestep.  First "competition" between processes for
! CH4 & O2 demand is done.  Then concentrations are apportioned into gas & liquid fractions; only the gas
! fraction is considered for diffusion in unsat.  Snow and lake water resistance to diffusion is added as
! a bulk term in the ground conductance (which is really a surface layer conductance), but concentrations
! are not tracked and oxidation is not allowed inside snow and lake water.
! Diffusivity is set based on soil texture and organic matter fraction. A Crank-Nicholson solution is used.
! Then CH4 diffusive flux is calculated and consistency is checked.

! !USES:
   use clmtype
   use clm_time_manager, only : get_step_size, get_nstep
   use TridiagonalMod  , only : Tridiagonal
   use ch4varcon,        only : organic_max, satpow, aereoxid, scale_factor_gasdiff, scale_factor_liqdiff, ch4frzout

!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc           ! column-index bounds
   integer, intent(in) :: num_methc          ! number of column soil points in column filter
   integer, intent(in) :: filter_methc(ubc-lbc+1)    ! column filter for soil points
   integer, intent(in)  :: jwt(lbc:ubc)      ! index of the soil layer right above the water table (-)
   integer, intent(in)  :: sat               ! 0 = unsaturated; 1 = saturated
   logical, intent(in) :: lake               ! function called with lake filter
   real(r8), intent(in) :: dtime_ch4         ! time step for ch4 calculations
!
! !CALLED FROM:
! subroutine ch4()
!
! !REVISION HISTORY:
! 9/16/08: Created by William J. Riley
! 9/15/09: Adapted for CLM 4 by Zack Subin
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   integer , pointer :: cgridcell(:)      ! gridcell of corresponding column
   real(r8), pointer :: z(:,:)            ! soil layer depth (m)
   real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: forc_pbot(:)      ! atmospheric pressure (Pa)
   real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity)
   real(r8), pointer :: dz(:,:)           ! layer thickness (m)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: zi(:,:)           ! interface level below a "z" level (m)
   real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2) [for snow & soil layers]
   real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2) [for snow & soil layers]
   real(r8), pointer :: h2osno(:)         ! snow water (mm H2O)
   real(r8), pointer :: snow_depth(:)         ! snow height (m)
   real(r8), pointer :: lake_icefrac(:,:) ! mass fraction of lake layer that is frozen
   real(r8), pointer :: bsw(:,:)          ! Clapp and Hornberger "b" (nlevgrnd)  
   real(r8), pointer :: cellorg(:,:)      ! column 3D org (kg/m^3 organic matter) (nlevgrnd)
   real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)
   integer , pointer :: snl(:)            ! negative of number of snow layers
   real(r8), pointer :: frac_h2osfc(:)    ! fraction of ground covered by surface water (0 to 1)
   real(r8), pointer :: h2osfc(:)         ! surface water (mm)
   real(r8), pointer :: t_h2osfc(:)       ! surface water temperature




! local pointers to implicit in/out arrays
   real(r8), pointer :: ch4_prod_depth(:,:) ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_oxid_depth(:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_aere_depth(:,:) ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_surf_aere(:)    ! Total column CH4 aerenchyma (mol/m2/s)
   real(r8), pointer :: ch4_ebul_depth(:,:) ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_ebul_total(:)   ! Total column CH4 ebullition (mol/m2/s)
   real(r8), pointer :: ch4_surf_ebul(:)    ! CH4 ebullition to atmosphere (mol/m2/s)
   real(r8), pointer :: o2_oxid_depth(:,:)  ! O2 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_decomp_depth(:,:)! O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: co2_decomp_depth(:,:)! CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: conc_ch4(:,:)       ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: c_atm(:,:)          ! CH4, O2, CO2 atmospheric conc  (mol/m3)
   real(r8), pointer :: ch4_surf_diff(:)    ! CH4 surface flux (mol/m2/s)
   real(r8), pointer :: conc_o2(:,:)        ! O2 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: grnd_ch4_cond(:)    ! tracer conductance for boundary layer [m/s]
   real(r8), pointer :: o2_aere_depth(:,:)  ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2stress(:,:)       ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   real(r8), pointer :: ch4stress(:,:)      ! Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,j,g,p,s,i,ll ! indices
   integer :: fc             ! soil filter column index
   integer :: fp             ! soil filter pft index
   integer  :: jtop(lbc:ubc) ! top level at each column
   integer :: iter           ! iteration counter when dtime_ch4 < dtime
   real(r8) :: dtime         ! land model time step (sec)
   real(r8) :: at (lbc:ubc,0:nlevsoi)  ! "a" vector for tridiagonal matrix
   real(r8) :: bt (lbc:ubc,0:nlevsoi)  ! "b" vector for tridiagonal matrix
   real(r8) :: ct (lbc:ubc,0:nlevsoi)  ! "c" vector for tridiagonal matrix
   real(r8) :: rt (lbc:ubc,0:nlevsoi)  ! "r" vector for tridiagonal solution
   real(r8) :: f_a                     ! air-filled fraction of available pore space
   real(r8) :: diffus (lbc:ubc,0:nlevsoi) !diffusivity (m2/s)
   real(r8) :: k_h_inv                    ! 1/Henry's Law Constant in Latm/mol
   real(r8) :: k_h_cc(lbc:ubc,0:nlevsoi,ngases)   ! ratio of mol/m3 in liquid to mol/m3 in gas
   real(r8) :: dzj                         ! 
   real(r8) :: dp1_zp1 (lbc:ubc,0:nlevsoi) ! diffusivity/delta_z for next j
   real(r8) :: dm1_zm1 (lbc:ubc,0:nlevsoi) ! diffusivity/delta_z for previous j
   real(r8) :: t_soisno_c                  ! soil temperature (C)  (-nlevsno+1:nlevsoi)
   real(r8) :: eps                         ! either epsilon_a or epsilon_w, depending on where in soil, wrt WT
   real(r8) :: deficit                     ! mol CH4 /m^2 that must be subtracted from diffusive flux to atm. to make up
                                           ! for keeping concentrations always above zero
   real(r8) :: conc_ch4_bef(lbc:ubc,1:nlevsoi) ! concentration at the beginning of the timestep
   real(r8) :: errch4(lbc:ubc)             ! Error (Mol CH4 /m^2) [+ = too much CH4]
   real(r8) :: conc_ch4_rel(lbc:ubc,0:nlevsoi) ! Concentration per volume of air or water
   real(r8) :: conc_o2_rel(lbc:ubc,0:nlevsoi)  ! Concentration per volume of air or water
   real(r8) :: conc_ch4_rel_old(lbc:ubc,0:nlevsoi) ! Concentration during last Crank-Nich. loop
   real(r8) :: h2osoi_vol_min(lbc:ubc,1:nlevsoi)   ! h2osoi_vol restricted to be <= watsat
   real(r8), parameter :: smallnumber = 1.e-12_r8
   real(r8) :: snowdiff                            ! snow diffusivity (m^2/s)
   real(r8) :: snowres(lbc:ubc)                    ! Cumulative Snow resistance (s/m). Also includes
   real(r8) :: pondres                             ! Additional resistance from ponding, up to pondmx water on top of top soil layer (s/m)
   real(r8) :: pondz                               ! Depth of ponding (m)
   real(r8) :: ponddiff                            ! Pondwater diffusivity (m^2/s)
   real(r8) :: spec_grnd_cond(lbc:ubc,1:ngases)    ! species grnd conductance (s/m)
   real(r8) :: airfrac                             ! air fraction in snow
   real(r8) :: waterfrac                           ! water fraction in snow
   real(r8) :: icefrac                             ! ice fraction in snow
   real(r8) :: epsilon_t (lbc:ubc,1:nlevsoi,1:ngases) !
   real(r8) :: epsilon_t_old (lbc:ubc,1:nlevsoi,1:ngases) !epsilon_t from last time step !Currently deprecated
   real(r8) :: source (lbc:ubc,1:nlevsoi,1:ngases)      !source
   real(r8) :: source_old (lbc:ubc,1:nlevsoi,1:ngases)  !source from last time step !Currently deprecated
   real(r8) :: om_frac                             ! organic matter fraction
   real(r8) :: o2demand, ch4demand                 ! mol/m^3/s
   real(r8) :: liqfrac(lbc:ubc, 1:nlevsoi)
   real(r8) :: capthick = 100._r8                  ! (mm) min thickness before assuming h2osfc is impermeable


   integer  :: nstep                       ! time step number
   character(len=32) :: subname='ch4_tran' ! subroutine name

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   !pft level
   !column level
   z         => clm3%g%l%c%cps%z
   dz        => clm3%g%l%c%cps%dz
   zi        => clm3%g%l%c%cps%zi
   t_soisno  => clm3%g%l%c%ces%t_soisno
   cgridcell => clm3%g%l%c%gridcell
   h2osoi_vol     => clm3%g%l%c%cws%h2osoi_vol
   h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
   h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
   watsat         => clm3%g%l%c%cps%watsat
   grnd_ch4_cond  => clm3%g%l%c%cps%pps_a%grnd_ch4_cond
   h2osno    => clm3%g%l%c%cws%h2osno
   snow_depth       => clm3%g%l%c%cps%snow_depth
   lake_icefrac => clm3%g%l%c%cws%lake_icefrac
   bsw       => clm3%g%l%c%cps%bsw
   cellorg   => clm3%g%l%c%cps%cellorg
   t_grnd    => clm3%g%l%c%ces%t_grnd
   snl       => clm3%g%l%c%cps%snl
   h2osfc         => clm3%g%l%c%cws%h2osfc
   frac_h2osfc    => clm3%g%l%c%cps%frac_h2osfc
   t_h2osfc       => clm3%g%l%c%ces%t_h2osfc


   !gridcell level
   c_atm          => clm3%g%gch4%c_atm
   forc_pbot      => clm_a2l%forc_pbot

   if (sat == 0) then ! unsaturated
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_unsat
      o2_oxid_depth  => clm3%g%l%c%cch4%o2_oxid_depth_unsat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_unsat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_unsat
      ch4_prod_depth => clm3%g%l%c%cch4%ch4_prod_depth_unsat
      ch4_aere_depth => clm3%g%l%c%cch4%ch4_aere_depth_unsat
      ch4_surf_aere  => clm3%g%l%c%cch4%ch4_surf_aere_unsat
      ch4_ebul_depth => clm3%g%l%c%cch4%ch4_ebul_depth_unsat
      ch4_ebul_total => clm3%g%l%c%cch4%ch4_ebul_total_unsat
      ch4_surf_ebul  => clm3%g%l%c%cch4%ch4_surf_ebul_unsat
      ch4_surf_diff  => clm3%g%l%c%cch4%ch4_surf_diff_unsat
      o2_decomp_depth=> clm3%g%l%c%cch4%o2_decomp_depth_unsat
      o2_aere_depth  => clm3%g%l%c%cch4%o2_aere_depth_unsat
      co2_decomp_depth=> clm3%g%l%c%cch4%co2_decomp_depth_unsat
      o2stress       => clm3%g%l%c%cch4%o2stress_unsat
      ch4stress      => clm3%g%l%c%cch4%ch4stress_unsat
   else ! saturated
      ch4_oxid_depth => clm3%g%l%c%cch4%ch4_oxid_depth_sat
      o2_oxid_depth  => clm3%g%l%c%cch4%o2_oxid_depth_sat
      conc_ch4       => clm3%g%l%c%cch4%conc_ch4_sat
      conc_o2        => clm3%g%l%c%cch4%conc_o2_sat
      ch4_prod_depth => clm3%g%l%c%cch4%ch4_prod_depth_sat
      ch4_aere_depth => clm3%g%l%c%cch4%ch4_aere_depth_sat
      ch4_surf_aere  => clm3%g%l%c%cch4%ch4_surf_aere_sat
      ch4_ebul_depth => clm3%g%l%c%cch4%ch4_ebul_depth_sat
      ch4_ebul_total => clm3%g%l%c%cch4%ch4_ebul_total_sat
      ch4_surf_ebul  => clm3%g%l%c%cch4%ch4_surf_ebul_sat
      ch4_surf_diff  => clm3%g%l%c%cch4%ch4_surf_diff_sat
      o2_decomp_depth=> clm3%g%l%c%cch4%o2_decomp_depth_sat
      o2_aere_depth  => clm3%g%l%c%cch4%o2_aere_depth_sat
      co2_decomp_depth=> clm3%g%l%c%cch4%co2_decomp_depth_sat
      o2stress       => clm3%g%l%c%cch4%o2stress_sat
      ch4stress      => clm3%g%l%c%cch4%ch4stress_sat
   endif

   ! Get land model time step
   dtime = get_step_size()
   nstep = get_nstep()

   ! Perform competition for oxygen and methane in each soil layer if demands over the course of the timestep
   ! exceed that available. Assign to each process in proportion to the quantity demanded in the absense of
   ! the limitation.
   do j = 1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)

         o2demand = o2_decomp_depth(c,j) + o2_oxid_depth(c,j) ! o2_decomp_depth includes autotrophic root respiration
         if (o2demand > 0._r8) then
            o2stress(c,j) = min((conc_o2(c,j) / dtime + o2_aere_depth(c,j)) / o2demand, 1._r8)
         else
            o2stress(c,j) = 1._r8
         end if

         ch4demand = ch4_oxid_depth(c,j) + ch4_aere_depth(c,j) + ch4_ebul_depth(c,j)
         if (ch4demand > 0._r8) then
            ch4stress(c,j) = min((conc_ch4(c,j) / dtime + ch4_prod_depth(c,j)) / ch4demand, 1._r8)
         else
            ch4stress(c,j) = 1._r8
         end if

         ! Resolve methane oxidation
         if (o2stress(c,j) < 1._r8 .or. ch4stress(c,j) < 1._r8) then
            if (ch4stress(c,j) <= o2stress(c,j)) then ! methane limited
               if (o2stress(c,j) < 1._r8) then
                  ! Recalculate oxygen limitation
                  o2demand = o2_decomp_depth(c,j)
                  if (o2demand > 0._r8) then
                     o2stress(c,j) = min( (conc_o2(c,j) / dtime + o2_aere_depth(c,j) - ch4stress(c,j)*o2_oxid_depth(c,j) ) &
                                       / o2demand, 1._r8)
                  else
                     o2stress(c,j) = 1._r8
                  end if
               end if
               ! Reset oxidation
               ch4_oxid_depth(c,j) = ch4_oxid_depth(c,j) * ch4stress(c,j)
               o2_oxid_depth(c,j) = o2_oxid_depth(c,j) * ch4stress(c,j)
            else                                      ! oxygen limited
               if (ch4stress(c,j) < 1._r8) then
                  ! Recalculate methane limitation
                  ch4demand = ch4_aere_depth(c,j) + ch4_ebul_depth(c,j)
                  if (ch4demand > 0._r8) then
                     ch4stress(c,j) = min( (conc_ch4(c,j) / dtime + ch4_prod_depth(c,j) - &
                                            o2stress(c,j)*ch4_oxid_depth(c,j)) / ch4demand, 1._r8)
                  else
                     ch4stress(c,j) = 1._r8
                  end if
               end if
               ! Reset oxidation
               ch4_oxid_depth(c,j) = ch4_oxid_depth(c,j) * o2stress(c,j)
               o2_oxid_depth(c,j) = o2_oxid_depth(c,j) * o2stress(c,j)
            end if
         end if

         ! Reset non-methanotroph demands
         ch4_aere_depth(c,j) = ch4_aere_depth(c,j) * ch4stress(c,j)
         ch4_ebul_depth(c,j) = ch4_ebul_depth(c,j) * ch4stress(c,j)
         o2_decomp_depth(c,j) = o2_decomp_depth(c,j) * o2stress(c,j)

      end do !c
   end do !j


   ! Accumulate ebullition to place in first layer above water table, or directly to atmosphere
   do j = 1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)
         if (j == 1) ch4_ebul_total(c) = 0._r8
         ch4_ebul_total(c) = ch4_ebul_total(c) + ch4_ebul_depth(c,j) * dz(c,j)
      enddo
   enddo


   ! Set the Henry's Law coefficients
   do j = 0,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)
         
         do s=1,2         
            if (j == 0) then
               k_h_inv = exp(-c_h_inv(s) * (1._r8 / t_grnd(c) - 1._r8 / kh_tbase) + log (kh_theta(s)))
               ! (4.12) Wania (L atm/mol)
               k_h_cc(c,j,s) = t_grnd(c) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
            else
               k_h_inv = exp(-c_h_inv(s) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(s)))
               ! (4.12) Wania (L atm/mol)
               k_h_cc(c,j,s) = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
            end if
         end do
      end do
   end do
         

   ! Set the source term for each species (no need to do j=0, since epsilon_t and source not used there)
   ! Note that because of the semi-implicit diffusion and the 30 min timestep combined with explicit
   ! sources, occasionally negative concentration will result. In this case it is brought to zero and the
   ! surface flux is adjusted to conserve. This results in some inaccuracy as compared to a shorter timestep
   ! or iterative solution.
   do j = 1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)

         if (aereoxid >= 0._r8) then
         ! First remove the CH4 oxidation that occurs at the base of root tissues (aere), and add to oxidation
            ch4_oxid_depth(c,j) = ch4_oxid_depth(c,j) + aereoxid * ch4_aere_depth(c,j)
            ch4_aere_depth(c,j) = ch4_aere_depth(c,j) - aereoxid * ch4_aere_depth(c,j)
         end if ! else oxygen is allowed to diffuse in via aerenchyma

         source(c,j,1) = ch4_prod_depth(c,j) - ch4_oxid_depth(c,j) - &
                         ch4_aere_depth(c,j) - ch4_ebul_depth(c,j) ! [mol/m3-total/s]
                         ! aerenchyma added to surface flux below
                         ! ebul added to soil depth just above WT
         if (source(c,j,1) + conc_ch4(c,j) / dtime < -1.e-12_r8) then
            write(iulog,*) 'Methane demands exceed methane available. Error in methane competition (mol/m^3/s), c,j:', &
                            source(c,j,1) + conc_ch4(c,j) / dtime, c, j
            g = cgridcell(c)
            write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
            call endrun( trim(subname)//' ERROR: Methane demands exceed methane available.' )
         else if (ch4stress(c,j) < 1._r8 .and. source(c,j,1) + conc_ch4(c,j) / dtime > 1.e-12_r8) then
            write(iulog,*) 'Methane limited, yet some left over. Error in methane competition (mol/m^3/s), c,j:', &
                            source(c,j,1) + conc_ch4(c,j) / dtime, c, j
            g = cgridcell(c)
            write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
            call endrun( trim(subname)//' ERROR: Methane limited, yet some left over.' )
         end if

         source(c,j,2) = -o2_oxid_depth(c,j) - o2_decomp_depth(c,j) + o2_aere_depth(c,j) ! O2 [mol/m3/s]
         if (source(c,j,2) + conc_o2(c,j) / dtime < -1.e-12_r8) then
            write(iulog,*) 'Oxygen demands exceed oxygen available. Error in oxygen competition (mol/m^3/s), c,j:', &
                            source(c,j,2) + conc_o2(c,j) / dtime, c, j
            g = cgridcell(c)
            write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
            call endrun( trim(subname)//' ERROR: Oxygen demands exceed oxygen available.' )
         else if (o2stress(c,j) < 1._r8 .and. source(c,j,2) + conc_o2(c,j) / dtime > 1.e-12_r8) then
            write(iulog,*) 'Oxygen limited, yet some left over. Error in oxygen competition (mol/m^3/s), c,j:', &
                             source(c,j,2) + conc_o2(c,j) / dtime, c, j
            g = cgridcell(c)
            write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
            call endrun( trim(subname)//' ERROR: Oxygen limited, yet some left over.' )
         end if

         conc_ch4_bef(c,j) = conc_ch4(c,j) !For Balance Check
      enddo ! fc
   enddo ! j

   ! Accumulate aerenchyma to add directly to atmospheric flux
   do j = 1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)
         if (j==1) ch4_surf_aere(c) = 0._r8
         ch4_surf_aere(c) = ch4_surf_aere(c) + ch4_aere_depth(c,j) * dz(c,j)
      enddo
   enddo

   ! Add in ebullition to source at depth just above WT
   do fc = 1, num_methc
      c = filter_methc(fc)
      if (jwt(c) /= 0) then
         source(c,jwt(c),1) = source(c,jwt(c),1) + ch4_ebul_total(c)/dz(c,jwt(c))
      endif
   enddo ! fc

   ! Calculate concentration relative to m^3 of air or water: needed for the diffusion
   do j = 0,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)
         g = cgridcell(c)

         if (j == 0) then
            conc_ch4_rel(c,j) = c_atm(g,1)
            conc_o2_rel(c,j)  = c_atm(g,2)
         else
            h2osoi_vol_min(c,j) = min(watsat(c,j), h2osoi_vol(c,j))
            if (ch4frzout) then
               liqfrac(c,j) = max(0.05_r8, (h2osoi_liq(c,j)/denh2o+smallnumber)/ &
                                    (h2osoi_liq(c,j)/denh2o+h2osoi_ice(c,j)/denice+smallnumber))
            else
               liqfrac(c,j) = 1._r8
            end if
            if (j <= jwt(c)) then  ! Above the WT
               do s=1,2
                  epsilon_t(c,j,s) = watsat(c,j)- (1._r8-k_h_cc(c,j,s))*h2osoi_vol_min(c,j)*liqfrac(c,j)
               end do
               ! Partition between the liquid and gas phases. The gas phase will drive the diffusion.
            else ! Below the WT
               do s=1,2
                  epsilon_t(c,j,s) = watsat(c,j)*liqfrac(c,j)
               end do
            end if
            conc_ch4_rel(c,j) = conc_ch4(c,j)/epsilon_t(c,j,1)
            conc_o2_rel(c,j)  = conc_o2(c,j) /epsilon_t(c,j,2)
         end if
      end do
   end do


   ! Loop over species
   do s = 1, 2 ! 1=CH4; 2=O2; 3=CO2


      ! Adjust the grnd_ch4_cond to keep it positive, and add the snow resistance & pond resistance
      do j = -nlevsno + 1,0
         do fc = 1, num_methc
            c = filter_methc (fc)
      
            if (j == -nlevsno + 1) then
               if (grnd_ch4_cond(c) < smallnumber .and. s==1) grnd_ch4_cond(c) = smallnumber
               ! Needed to prevent overflow when ground is frozen, e.g. for lakes
               snowres(c) = 0._r8
            end if

            ! Add snow resistance
            if (j >= snl(c) + 1) then
               t_soisno_c = t_soisno(c,j) - tfrz
               icefrac = h2osoi_ice(c,j)/denice/dz(c,j)
               waterfrac = h2osoi_liq(c,j)/denh2o/dz(c,j)
               airfrac = max(1._r8 - icefrac - waterfrac, 0._r8)
               ! Calculate snow diffusivity
               if (airfrac > 0.05_r8) then
                  f_a = airfrac / (airfrac + waterfrac)
                  eps = airfrac ! Air-filled fraction of total snow volume
                  ! Use Millington-Quirk Expression, as hydraulic properties (bsw) not available
                  snowdiff = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
                             f_a**(10._r8/3._r8) / (airfrac+waterfrac)**2 &
                           * scale_factor_gasdiff
               else !solute diffusion in water only
                  eps = waterfrac  ! Water-filled fraction of total soil volume
                  snowdiff = eps**satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                           * scale_factor_liqdiff
               end if
               snowdiff = max(snowdiff, smallnumber)
               snowres(c) = snowres(c) + dz(c,j)/snowdiff
            end if

            if (j == 0) then ! final loop
               ! Add pond resistance
               pondres = 0._r8

               ! First old pond formulation up to pondmx
               if (.not. lake .and. snl(c) == 0 .and. h2osoi_vol(c,1) > watsat(c,1)) then
                  t_soisno_c = t_soisno(c,1) - tfrz
                  if (t_soisno(c,1) <= tfrz) then
                     ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                              * (h2osoi_liq(c,1)/denh2o+smallnumber)/ &
                                (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice+smallnumber) &
                              * scale_factor_liqdiff
                  else ! Unfrozen
                     ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                              * scale_factor_liqdiff
                  end if
                  pondz = dz(c,1) * (h2osoi_vol(c,1) - watsat(c,1))
                  pondres = pondz / ponddiff
               end if            

               ! Now add new h2osfc form
               if (.not. lake .and. sat == 1 .and. frac_h2osfc(c) > 0._r8 .and. t_h2osfc(c) >= tfrz) then
                  t_soisno_c = t_h2osfc(c) - tfrz
                  ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                              * scale_factor_liqdiff
                  pondz = h2osfc(c) / 1000._r8 / frac_h2osfc(c) ! Assume all h2osfc corresponds to sat area
                          ! mm      /  mm/m
                  pondres = pondres + pondz / ponddiff
               else if (.not. lake .and. sat == 1 .and. frac_h2osfc(c) > 0._r8 .and. &
                        h2osfc(c)/frac_h2osfc(c) > capthick) then ! Assuming short-circuit logic will avoid FPE here.
                       ! assume surface ice is impermeable
                  pondres = 1/smallnumber
               end if

               spec_grnd_cond(c,s) = 1._r8/(1._r8/grnd_ch4_cond(c) + snowres(c) + pondres)
            end if

         end do ! fc
      end do ! j

      ! Determine gas diffusion and fraction of open pore (f_a)
      do j = 1,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)
            g = cgridcell(c)

            t_soisno_c = t_soisno(c,j) - tfrz

            if (j <= jwt(c)) then  ! Above the WT
               f_a = 1._r8 - h2osoi_vol_min(c,j) / watsat(c,j)
               ! Provisionally calculate diffusivity as linear combination of the Millington-Quirk 
               ! expression in Wania (for peat) & Moldrup (for mineral soil)
               eps =  watsat(c,j)-h2osoi_vol_min(c,j) ! Air-filled fraction of total soil volume
               if (organic_max > 0._r8) then
                  om_frac = min(cellorg(c,j)/organic_max, 1._r8)
                  ! Use first power, not square as in iniTimeConst
               else
                  om_frac = 1._r8
               end if
               diffus (c,j) = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
                              (om_frac * f_a**(10._r8/3._r8) / watsat(c,j)**2._r8 + &
                               (1._r8-om_frac) * eps**2._r8 * f_a**(3._r8 / bsw(c,j)) ) &
                            * scale_factor_gasdiff
            else ! Below the WT use saturated diffusivity and only water in epsilon_t
               ! Note the following is not currently corrected for the effect on diffusivity of excess ice in soil under
               ! lakes (which is currently experimental only).
               eps = watsat(c,j)  ! Water-filled fraction of total soil volume
               diffus (c,j) = eps**satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                            * scale_factor_liqdiff
               if (t_soisno(c,j)<=tfrz) then
                  diffus(c,j) = diffus(c,j)*(h2osoi_liq(c,j)/denh2o+smallnumber)/ &
                                (h2osoi_liq(c,j)/denh2o+h2osoi_ice(c,j)/denice+smallnumber)
               end if
            endif ! Above/below the WT
            diffus(c,j) = max(diffus(c,j), smallnumber) ! Prevent overflow

         enddo ! fp
      enddo ! j

      do j = 1,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)

            ! Set up coefficients for tridiagonal solver.
            if (j == 1 .and. j /= jwt(c) .and. j /= jwt(c)+1) then
               dm1_zm1(c,j) = 1._r8/(1._r8/spec_grnd_cond(c,s)+dz(c,j)/(diffus(c,j)*2._r8))
               ! replace Diffusivity / Delta_z by conductance (grnd_ch4_cond) for top layer
               dp1_zp1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j+1)/diffus(c,j+1))
            else if (j == 1 .and. j == jwt(c)) then
               dm1_zm1(c,j) = 1._r8/(1._r8/spec_grnd_cond(c,s)+dz(c,j)/(diffus(c,j)*2._r8))
               ! layer resistance mult. by k_h_cc for dp1_zp1 term
               dp1_zp1(c,j) = 2._r8/(dz(c,j)*k_h_cc(c,j,s)/diffus(c,j)+dz(c,j+1)/diffus(c,j+1))
            else if (j == 1) then ! water table at surface: multiply ground resistance by k_h_cc
               dm1_zm1(c,j) = 1._r8/(k_h_cc(c,j-1,s)/spec_grnd_cond(c,s)+dz(c,j)/(diffus(c,j)*2._r8))
               ! air concentration will be mult. by k_h_cc below
               dp1_zp1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j+1)/diffus(c,j+1))
            else if (j <= nlevsoi-1 .and. j /= jwt(c) .and. j /= jwt(c)+1) then
               dm1_zm1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j-1)/diffus(c,j-1))
               dp1_zp1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j+1)/diffus(c,j+1))
            else if (j <= nlevsoi-1 .and. j == jwt(c)) then ! layer resistance mult. by k_h_cc for dp1_zp1 term
               dm1_zm1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j-1)/diffus(c,j-1))
               dp1_zp1(c,j) = 2._r8/(dz(c,j)*k_h_cc(c,j,s)/diffus(c,j)+dz(c,j+1)/diffus(c,j+1))
               ! Concentration in layer will be mult. by k_h_cc below
            else if (j <= nlevsoi-1) then ! j==jwt+1: layer above resistance mult. by k_h_cc for dm1_zm1 term
               dm1_zm1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j-1)*k_h_cc(c,j-1,s)/diffus(c,j-1))
               ! Concentration in layer above will be mult. by k_h_cc below
               dp1_zp1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j+1)/diffus(c,j+1))
            else if (j /= jwt(c)+1) then ! j ==nlevsoi
               dm1_zm1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j-1)/diffus(c,j-1))
            else                    ! jwt == nlevsoi-1: layer above resistance mult. by k_h_cc for dm1_zm1 term
               dm1_zm1(c,j) = 2._r8/(dz(c,j)/diffus(c,j)+dz(c,j-1)*k_h_cc(c,j-1,s)/diffus(c,j-1))
            end if
         enddo ! fp; pft
      end do ! j; nlevsoi

      ! Perform a second loop for the tridiagonal coefficients since need dp1_zp1 and dm1_z1 at each depth
      do j = 0,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)
            g = cgridcell(c)

            conc_ch4_rel_old(c,j) = conc_ch4_rel(c,j)

            if (j > 0) dzj = dz(c,j)
            if (j == 0) then ! top layer (atmosphere) doesn't change regardless of where WT is
               at(c,j) = 0._r8
               bt(c,j) = 1._r8
               ct(c,j) = 0._r8
               rt(c,j) = c_atm(g,s) ! 0th level stays at constant atmospheric conc
            elseif (j < nlevsoi .and. j == jwt(c)) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
               at(c,j) = -0.5_r8 / dzj * dm1_zm1(c,j)
               bt(c,j) = epsilon_t(c,j,s) / dtime_ch4 + 0.5_r8 / dzj * (dp1_zp1(c,j)*k_h_cc(c,j,s) + dm1_zm1(c,j))
               ct(c,j) = -0.5_r8 / dzj * dp1_zp1(c,j)
            elseif (j < nlevsoi .and. j == jwt(c)+1) then
               ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
               at(c,j) = -0.5_r8 / dzj * dm1_zm1(c,j) * k_h_cc(c,j-1,s)
               bt(c,j) = epsilon_t(c,j,s) / dtime_ch4 + 0.5_r8 / dzj * (dp1_zp1(c,j) + dm1_zm1(c,j))
               ct(c,j) = -0.5_r8 / dzj * dp1_zp1(c,j)
            elseif (j < nlevsoi) then
               at(c,j) = -0.5_r8 / dzj * dm1_zm1(c,j)
               bt(c,j) = epsilon_t(c,j,s) / dtime_ch4 + 0.5_r8 / dzj * (dp1_zp1(c,j) + dm1_zm1(c,j))
               ct(c,j) = -0.5_r8 / dzj * dp1_zp1(c,j)
            else if (j == nlevsoi .and. j== jwt(c)+1) then
               ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
               at(c,j) = -0.5_r8 / dzj * dm1_zm1(c,j) * k_h_cc(c,j-1,s)
               bt(c,j) = epsilon_t(c,j,s) / dtime_ch4 + 0.5_r8 / dzj * dm1_zm1(c,j)
               ct(c,j) = 0._r8
            else ! j==nlevsoi and jwt<nlevsoi-1 or jwt==nlevsoi: 0 flux at bottom
               at(c,j) = -0.5_r8 / dzj * dm1_zm1(c,j)
               bt(c,j) = epsilon_t(c,j,s) / dtime_ch4 + 0.5_r8 / dzj * dm1_zm1(c,j)
               ct(c,j) = 0._r8
            endif
         enddo ! fp; pft
      enddo ! j; nlevsoi

      do fc = 1, num_methc
         c = filter_methc (fc)
         jtop(c) = 0
      end do

      if (s == 1) then  ! CH4

         ! Set rt, since it depends on conc
         do j = 1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc (fc)

               ! For correct balance, deprecate source_old.
               ! The source terms are effectively constant over the timestep.
               source_old(c,j,s) = source(c,j,s)
               ! source_old could be removed later
               epsilon_t_old(c,j,s) = epsilon_t(c,j,s)
               ! epsilon_t acts like source also
               dzj = dz(c,j)
               if (j < nlevsoi .and. j == jwt(c)) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_ch4_rel(c,j) +           &
                        0.5_r8 / dzj * (dp1_zp1(c,j) * (conc_ch4_rel(c,j+1)-conc_ch4_rel(c,j)*k_h_cc(c,j,s)) - &
                        dm1_zm1(c,j) * (conc_ch4_rel(c,j)  -conc_ch4_rel(c,j-1))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               elseif (j < nlevsoi .and. j == jwt(c)+1) then
                  ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_ch4_rel(c,j) +           &
                        0.5_r8 / dzj * (dp1_zp1(c,j) * (conc_ch4_rel(c,j+1)-conc_ch4_rel(c,j)) - &
                        dm1_zm1(c,j) * (conc_ch4_rel(c,j) -conc_ch4_rel(c,j-1)*k_h_cc(c,j-1,s))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               elseif (j < nlevsoi) then
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_ch4_rel(c,j) +           &
                        0.5_r8 / dzj * (dp1_zp1(c,j) * (conc_ch4_rel(c,j+1)-conc_ch4_rel(c,j)) - &
                        dm1_zm1(c,j) * (conc_ch4_rel(c,j)  -conc_ch4_rel(c,j-1))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               else if (j == nlevsoi .and. j== jwt(c)+1) then
                  ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_ch4_rel(c,j) +           &
                        0.5_r8 / dzj * ( - dm1_zm1(c,j) * (conc_ch4_rel(c,j) -conc_ch4_rel(c,j-1)*k_h_cc(c,j-1,s))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               else  !j==nlevsoi
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_ch4_rel(c,j) +           &
                        0.5_r8 / dzj * ( - dm1_zm1(c,j) * (conc_ch4_rel(c,j)  -conc_ch4_rel(c,j-1))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               endif
               epsilon_t_old(c,j,s) = epsilon_t(c,j,s)
               source_old(c,j,s) = source(c,j,s)

            enddo ! fc; column
         enddo ! j; nlevsoi

         call Tridiagonal(lbc, ubc, 0, nlevsoi, jtop, num_methc, filter_methc, &
                          at, bt, ct, rt, conc_ch4_rel(lbc:ubc,0:nlevsoi))

         ! Calculate net ch4 flux to the atmosphere from the surface (+ to atm)
         do fc = 1, num_methc
            c = filter_methc (fc)
            g = cgridcell(c)
            if (jwt(c) /= 0) then ! WT not at the surface
               ch4_surf_diff(c) = dm1_zm1(c,1) * ( (conc_ch4_rel(c,1)+conc_ch4_rel_old(c,1))/2._r8 &
                                                        - c_atm(g,s)) ! [mol/m2/s]
               ch4_surf_ebul(c) = 0._r8 ! all the ebullition has already come out in the soil column (added to source)
               ! Try adding directly to atm. to prevent destabilization of diffusion
               !ch4_surf_ebul(c) = ch4_ebul_total(c) ! [mol/m2/s]
            else ! WT at the surface; i.e., jwt(c)==0
               ch4_surf_diff(c) = dm1_zm1(c,1) * ( (conc_ch4_rel(c,1)+conc_ch4_rel_old(c,1))/2._r8 &
                                                         - c_atm(g,s)*k_h_cc(c,0,s)) ! [mol/m2/s]
                                                         ! atmospheric concentration gets mult. by k_h_cc as above
               ch4_surf_ebul(c) = ch4_ebul_total(c) ! [mol/m2/s]
            endif
         enddo

         ! Ensure that concentrations stay above 0
         ! This should be done after the flux, so that the flux calculation is consistent.
         do j = 1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc (fc)

               if (conc_ch4_rel(c,j) < 0._r8) then
                  deficit = - conc_ch4_rel(c,j)*epsilon_t(c,j,1)*dz(c,j)  ! Mol/m^2 added
                  if (deficit > 1.e-3_r8 * scale_factor_gasdiff) then
                     if (deficit > 1.e-2_r8) then
                        write(iulog,*)'Note: sink > source in ch4_tran, sources are changing quickly relative to diffusion timestep, and/or diffusion is rapid.'
                        g = cgridcell(c)
                        write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
                        write(iulog,*)'This typically occurs when there is a larger than normal diffusive flux.'
                        write(iulog,*)'If this occurs frequently, consider reducing land model (or methane model) timestep, or reducing the max. sink per timestep in the methane model.'
                     end if
                     write(iulog,*) 'Negative conc. in ch4tran. c,j,deficit (mol):',c,j,deficit
                  end if
                  conc_ch4_rel(c,j) = 0._r8
                  ! Subtract deficit
                  ch4_surf_diff(c) = ch4_surf_diff(c) - deficit/dtime_ch4
               end if
            enddo
         enddo


      elseif (s == 2) then  ! O2

         ! Set rt, since it depends on conc
         do j = 1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc (fc)

               ! For correct balance, deprecate source_old.
               source_old(c,j,s) = source(c,j,s)
               ! source_old could be removed later
               epsilon_t_old(c,j,s) = epsilon_t(c,j,s)
               ! epsilon_t acts like source also
               dzj     = dz(c,j)
               if (j < nlevsoi .and. j == jwt(c)) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_o2_rel(c,j) +           &
                        0.5_r8 / dzj * (dp1_zp1(c,j) * (conc_o2_rel(c,j+1)-conc_o2_rel(c,j)*k_h_cc(c,j,s)) - &
                        dm1_zm1(c,j) * (conc_o2_rel(c,j)  -conc_o2_rel(c,j-1))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               elseif (j < nlevsoi .and. j == jwt(c)+1) then
                  ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_o2_rel(c,j) +           &
                        0.5_r8 / dzj * (dp1_zp1(c,j) * (conc_o2_rel(c,j+1)-conc_o2_rel(c,j)) - &
                        dm1_zm1(c,j) * (conc_o2_rel(c,j) -conc_o2_rel(c,j-1)*k_h_cc(c,j-1,s))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               elseif (j < nlevsoi) then
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_o2_rel(c,j) +           &
                        0.5_r8 / dzj * (dp1_zp1(c,j) * (conc_o2_rel(c,j+1)-conc_o2_rel(c,j)) - &
                        dm1_zm1(c,j) * (conc_o2_rel(c,j)  -conc_o2_rel(c,j-1))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               else if (j == nlevsoi .and. j== jwt(c)+1) then
                  ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_o2_rel(c,j) +           &
                        0.5_r8 / dzj * ( - dm1_zm1(c,j) * (conc_o2_rel(c,j) -conc_o2_rel(c,j-1)*k_h_cc(c,j-1,s))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               else  !j==nlevsoi
                  rt(c,j) = epsilon_t_old(c,j,s) / dtime_ch4 * conc_o2_rel(c,j) +           &
                        0.5_r8 / dzj * ( - dm1_zm1(c,j) * (conc_o2_rel(c,j)  -conc_o2_rel(c,j-1))) + &
                        0.5_r8 * (source(c,j,s) + source_old(c,j,s))
               endif
               epsilon_t_old(c,j,s) = epsilon_t(c,j,s)
               source_old(c,j,s) = source(c,j,s)

            enddo ! fc; column
         enddo ! j; nlevsoi

         call Tridiagonal(lbc, ubc, 0, nlevsoi, jtop, num_methc, filter_methc, &
                          at, bt, ct, rt, conc_o2_rel(lbc:ubc,0:nlevsoi))
         ! Ensure that concentrations stay above 0
         do j = 1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc (fc)
               g = cgridcell(c)
               conc_o2_rel(c,j) = max (conc_o2_rel(c,j), 1.e-12_r8)
               ! In case of pathologically large aerenchyma conductance. Should be OK in general but
               ! this will maintain stability even if a PFT with very small weight somehow has an absurd NPP or LAI.
               ! Also, oxygen above ambient will probably bubble.
               conc_o2_rel(c,j) = min (conc_o2_rel(c,j), c_atm(g,2)/epsilon_t(c,j,2))
            enddo
         enddo
       
      endif  ! species

   enddo  ! species

   ! Update absolute concentrations per unit volume
   do j = 1,nlevsoi ! No need to update the atm. level concentrations
      do fc = 1, num_methc
         c = filter_methc (fc)

         conc_ch4(c,j) = conc_ch4_rel(c,j)*epsilon_t(c,j,1)
         conc_o2(c,j)  = conc_o2_rel(c,j) *epsilon_t(c,j,2)
      end do
   end do

   ! Do Balance Check and absorb small
   !    discrepancy into surface flux.
   do j = 1,nlevsoi
      do fc = 1, num_methc
         c = filter_methc (fc)

         if (j == 1) errch4(c) = 0._r8
         errch4(c) = errch4(c) + (conc_ch4(c,j) - conc_ch4_bef(c,j))*dz(c,j)
         errch4(c) = errch4(c) - ch4_prod_depth(c,j)*dz(c,j)*dtime
         errch4(c) = errch4(c) + ch4_oxid_depth(c,j)*dz(c,j)*dtime
      end do
   end do

   do fc = 1, num_methc
      c = filter_methc (fc)

      ! For history make sure that grnd_ch4_cond includes snow, for methane diffusivity
      grnd_ch4_cond(c) = spec_grnd_cond(c,1)

      errch4(c) = errch4(c) + (ch4_surf_aere(c) + ch4_surf_ebul(c) + ch4_surf_diff(c))*dtime

      if (abs(errch4(c)) < 1.e-8_r8) then
         ch4_surf_diff(c) = ch4_surf_diff(c) - errch4(c)/dtime
      else ! errch4 > 1e-8 mol / m^2 / timestep
         write(iulog,*)'CH4 Conservation Error in CH4Mod during diffusion, nstep, c, errch4 (mol /m^2.timestep)', &
                       nstep,c,errch4(c)
         g = cgridcell(c)
         write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(g),clm3%g%londeg(g)
         call endrun( trim(subname)//' ERROR: CH4 Conservation Error in CH4Mod during diffusion' )
      end if
   end do


end subroutine ch4_tran


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_jwt 
!
! !INTERFACE:
subroutine get_jwt (lbc, ubc, num_methc, filter_methc, jwt)
!
! !DESCRIPTION:
! Finds the first unsaturated layer going up. Also allows a perched water table over ice.

! !USES:
   use clmtype

!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc           ! column-index bounds
   integer, intent(in) :: num_methc          ! number of column soil points in column filter
   integer, intent(in) :: filter_methc(ubc-lbc+1)    ! column filter for soil points
   integer, intent(out)  :: jwt(lbc:ubc)     ! index of the soil layer right above the water table (-)
!
! !CALLED FROM:
! subroutine ch4()
!
! !REVISION HISTORY:
! 9/16/08: Created by William J. Riley
! 8/27/09: Modified by Zack Subin for CLM 4
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   real(r8), pointer :: h2osoi_vol(:,:) ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   real(r8), pointer :: watsat(:,:)     !volumetric soil water at saturation (porosity)
   real(r8), pointer :: t_soisno(:,:)   ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,j,perch! indices
   integer :: fc       ! filter column index

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   h2osoi_vol     => clm3%g%l%c%cws%h2osoi_vol
   watsat         => clm3%g%l%c%cps%watsat
   t_soisno       => clm3%g%l%c%ces%t_soisno

   ! The layer index of the first unsaturated layer, i.e., the layer right above
   ! the water table.
   ! ZS: Loop is currently not vectorized.
    do fc = 1, num_methc
       c = filter_methc(fc)

       ! Check to see if any soil layers are frozen and saturated.  If so, start looking at the first layer above the top
       ! such layer.  This is potentially important for perched water tables in the Tundra.

       perch = nlevsoi
       do j = nlevsoi, 1, -1
          if (t_soisno(c,j) < tfrz .and. h2osoi_vol(c,j) > f_sat * watsat(c,j)) then
             ! strictly less than freezing because it could be permeable otherwise
             perch = j-1
          end if
       end do
       jwt(c) = perch

       do j = perch, 2, -1
          if(h2osoi_vol(c,j) > f_sat * watsat(c,j) .and. h2osoi_vol(c,j-1) < f_sat * watsat(c,j-1)) then
             jwt(c) = j-1
             exit
          end if
       enddo
       if (jwt(c) == perch .and. h2osoi_vol(c,1) > f_sat * watsat(c,1)) then ! missed that the top layer is saturated
          jwt(c) = 0
       endif
    end do

end subroutine get_jwt

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4annualupdate
!
! !INTERFACE:
subroutine ch4annualupdate(lbc, ubc, lbp, ubp, num_methc, filter_methc, num_methp, filter_methp)
!
! !DESCRIPTION: Annual mean fields.
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size, get_days_per_year, get_nstep
   use clm_varcon      , only: secspday
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc                ! column bounds
   integer, intent(in) :: lbp, ubp                ! pft bounds
   integer, intent(in) :: num_methc               ! number of soil columns in filter
   integer, intent(in) :: filter_methc(ubc-lbc+1) ! filter for soil columns
   integer, intent(in) :: num_methp               ! number of soil points in pft filter
   integer, intent(in) :: filter_methp(ubp-lbp+1) ! pft filter for soil points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 11/16/09: Zack Subin
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   real(r8), pointer :: somhr(:)          ! (gC/m2/s) soil organic matter heterotrophic respiration
   real(r8), pointer :: finundated(:)     ! fractional inundated area in soil column
   real(r8), pointer :: agnpp(:)          ! (gC/m2/s) aboveground NPP
   real(r8), pointer :: bgnpp(:)          ! (gC/m2/s) belowground NPP
   integer , pointer :: pcolumn(:)        ! index into column level quantities
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: annsum_counter(:)  ! seconds since last annual accumulator turnover
   ! This will point to a CH4 version, not the CN version.
   real(r8), pointer :: tempavg_somhr(:)   ! temporary average SOM heterotrophic resp. (gC/m2/s)
   real(r8), pointer :: annavg_somhr(:)    ! annual average SOM heterotrophic resp. (gC/m2/s)
   real(r8), pointer :: tempavg_finrw(:)   ! respiration-weighted annual average of finundated
   real(r8), pointer :: annavg_finrw(:)    ! respiration-weighted annual average of finundated
   real(r8), pointer :: tempavg_agnpp(:)   ! temporary average above-ground NPP (gC/m2/s)
   real(r8), pointer :: annavg_agnpp(:)    ! annual average above-ground NPP (gC/m2/s)
   real(r8), pointer :: tempavg_bgnpp(:)   ! temporary average below-ground NPP (gC/m2/s)
   real(r8), pointer :: annavg_bgnpp(:)    ! annual average below-ground NPP (gC/m2/s)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p       ! indices
   integer :: fc        ! soil column filter indices
   integer :: fp        ! soil pft filter indices
   real(r8):: dt        ! time step (seconds)
   real(r8):: secsperyear
   logical :: newrun

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays
   annsum_counter    => clm3%g%l%c%cch4%annsum_counter
   tempavg_somhr     => clm3%g%l%c%cch4%tempavg_somhr
   annavg_somhr      => clm3%g%l%c%cch4%annavg_somhr
   tempavg_finrw     => clm3%g%l%c%cch4%tempavg_finrw
   annavg_finrw      => clm3%g%l%c%cch4%annavg_finrw
   pcolumn           => clm3%g%l%c%p%column
   agnpp             => clm3%g%l%c%p%pcf%agnpp
   bgnpp             => clm3%g%l%c%p%pcf%bgnpp
   tempavg_agnpp     => clm3%g%l%c%p%pcf%tempavg_agnpp
   annavg_agnpp      => clm3%g%l%c%p%pcf%annavg_agnpp
   tempavg_bgnpp     => clm3%g%l%c%p%pcf%tempavg_bgnpp
   annavg_bgnpp      => clm3%g%l%c%p%pcf%annavg_bgnpp
   somhr             => clm3%g%l%c%ccf%somhr
   finundated        => clm3%g%l%c%cws%finundated

   ! set time steps
   dt = real(get_step_size(), r8)
   secsperyear = real( get_days_per_year() * secspday, r8)

   newrun = .false.

   ! column loop
   do fc = 1,num_methc
      c = filter_methc(fc)

      if (annsum_counter(c) == spval) then
      ! These variables are now in restart files for completeness, but might not be in inicFile and are not.
      ! set for arbinit.
         newrun = .true.
         annsum_counter(c)    = 0._r8
         tempavg_somhr(c)     = 0._r8
         tempavg_finrw(c)     = 0._r8
      end if

      annsum_counter(c) = annsum_counter(c) + dt
   end do

   ! pft loop
   do fp = 1,num_methp
      p = filter_methp(fp)

      if (newrun .or. tempavg_agnpp(p) == spval) then ! Extra check needed because for back-compatibility
         tempavg_agnpp(p) = 0._r8
         tempavg_bgnpp(p) = 0._r8
      end if
   end do

   do fc = 1,num_methc
      c = filter_methc(fc)
      if (annsum_counter(c) >= secsperyear) then

         ! update annual average somhr
         annavg_somhr(c)      =  tempavg_somhr(c)
         tempavg_somhr(c)     = 0._r8

         ! update annual average finrw
         if (annavg_somhr(c) > 0._r8) then
            annavg_finrw(c)      =  tempavg_finrw(c) / annavg_somhr(c)
         else
            annavg_finrw(c)      = 0._r8
         end if
         tempavg_finrw(c)     = 0._r8
      else
         tempavg_somhr(c)     = tempavg_somhr(c) + dt/secsperyear * somhr(c)
         tempavg_finrw(c)     = tempavg_finrw(c) + dt/secsperyear * finundated(c) * somhr(c)
      end if
   end do

   do fp = 1,num_methp
      p = filter_methp(fp)
      c = pcolumn(p)
      if (annsum_counter(c) >= secsperyear) then

         annavg_agnpp(p) = tempavg_agnpp(p)
         tempavg_agnpp(p) = 0._r8

         annavg_bgnpp(p) = tempavg_bgnpp(p)
         tempavg_bgnpp(p) = 0._r8

      else
         tempavg_agnpp(p) = tempavg_agnpp(p) + dt/secsperyear * agnpp(p)
         tempavg_bgnpp(p) = tempavg_bgnpp(p) + dt/secsperyear * bgnpp(p)
      end if
   end do

   ! column loop
   do fc = 1,num_methc
      c = filter_methc(fc)
      if (annsum_counter(c) >= secsperyear) annsum_counter(c) = 0._r8
   end do

end subroutine ch4annualupdate

#endif
!defined LCH4

end module ch4Mod
