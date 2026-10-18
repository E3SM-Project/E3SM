module SoilTemperatureMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates snow and soil temperatures including phase change
  !
  ! Local work arrays are compressed onto the column filter: they are indexed by
  ! the filter index fc (1:num_nolakec) rather than by the column index c.
  ! Column-level state (col_pp, col_es, ...) is still indexed by c = filter(fc).
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_infnan_mod    , only : nan => shr_infnan_nan
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use elm_varctl        , only : iulog, use_finetop_rad
  use elm_varcon        , only : spval, cnfac
  use elm_varpar        , only : nlevsno, nlevgrnd, nlevurb
  use column_varcon     , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon     , only : icol_road_perv, icol_road_imperv
  use UrbanParamsType   , only : urbanparams_type
  use CanopyStateType   , only : canopystate_type
  use SolarAbsorbedType , only : solarabs_type
  use SoilStateType     , only : soilstate_type
  use TopounitDataType  , only : top_af
  use LandunitType      , only : lun_pp
  use LandunitDataType  , only : lun_es, lun_ef
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_es, col_ef, col_ws, col_wf
  use VegetationType    , only : veg_pp
  use VegetationDataType, only : veg_ef, veg_wf
  use GridcellType      , only : grc_pp
  use timeinfoMod
  use perf_mod          , only : t_startf, t_stopf
  use shr_const_mod     , only : SHR_CONST_PI
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilTemperature
  !
  !    -> SoilTemperature:      soil/snow and ground temperatures
  !          -> SoilTermProp    thermal conductivities and heat capacities
  !          -> BandDiagonal    band diagonal matrix solution
  !          -> PhaseChange     phase change of liquid/ice contents
  !
  ! (1) Snow and soil temperatures
  !     o The volumetric heat capacity is calculated as a linear combination
  !       in terms of the volumetric fraction of the constituent phases.
  !     o The thermal conductivity of soil is computed from
  !       the algorithm of Johansen (as reported by Farouki 1981), and the
  !       conductivity of snow is from the formulation used in
  !       SNTHERM (Jordan 1991).
  !     o Boundary conditions:
  !       F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
  !     o Soil / snow temperature is predicted from heat conduction
  !       in 10 soil layers and up to 5 snow layers.
  !       The thermal conductivities at the interfaces between two
  !       neighboring layers (j, j+1) are derived from an assumption that
  !       the flux across the interface is equal to that from the node j
  !       to the interface and the flux from the interface to the node j+1.
  !       The equation is solved using the Crank-Nicholson method and
  !       results in a tridiagonal system equation.
  ! (2) Phase change
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilThermProp                  ! Set therm conduct. and heat cap of snow/soil layers
  private :: PhaseChangeH2osfc              ! When surface water freezes move ice to bottom snow layer
  private :: PhaseChange_beta               ! Calculation of the phase change within snow and soil layers
  private :: ComputeGroundHeatFluxAndDeriv  ! Computes G and dG/dT on surface of standing water, snow and soil
  private :: ComputeHeatDiffFluxAndFactor   ! Heat diffusion at layer interface and factor used in setting up of banded matrix
  private :: SolveTemperature               ! Assembles and solves the banded penta-diagonal system of equations
  private :: SetRHSVec                      ! Sets up the RHS vector for the numerical solution of temperature for snow/standing-water/soil
  private :: SetRHSVec_Snow                 ! Sets up the RHS vector corresponding to snow layers
  private :: SetRHSVec_StandingSurfaceWater ! Sets up the RHS vector corresponding to standing water layers
  private :: SetRHSVec_Soil                 ! Sets up the RHS vector corresponding to soil layers
  private :: SetMatrix                      ! Sets up the matrix for the numerical solution of temperature for snow/standing-water/soil
  private :: SetMatrix_Snow                 ! Set up the matrix entries corresponding to snow layers and snow-soil interaction
  private :: SetMatrix_Soil                 ! Set up the matrix entries corresponding to soil layers and soil-snow interaction
  private :: SetMatrix_StandingSurfaceWater ! Set up the matrix entries corresponding to standing surface water and its interaction with soil

  real(r8), private, parameter :: thin_sfclayer = 1.0e-6_r8   ! Threshold for thin surface layer
  integer , private, parameter :: nband = 5                   ! number of bands of the banded matrix
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilTemperature(bounds, num_urbanl, filter_urbanl, num_nolakec, filter_nolakec, &
       urbanparams_vars, canopystate_vars, solarabs_vars, soilstate_vars)
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
    use elm_varcon               , only : cpliq
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_urbanl                         ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:)                   ! urban landunit filter
    integer                , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                       ! indices
    integer  :: fc                                                          ! lake filtered column indices
    integer  :: fl                                                          ! urban filtered landunit indices
    real(r8) :: dtime                                                       ! land model time step (sec)
    integer  :: jtop(num_nolakec)                                           ! top level at each column
    integer  :: jbot(num_nolakec)                                           ! bottom level at each column
    real(r8) :: cv (num_nolakec,-nlevsno+1:nlevgrnd)                        ! heat capacity [J/(m2 K)]
    real(r8) :: tk (num_nolakec,-nlevsno+1:nlevgrnd)                        ! thermal conductivity [W/(m K)]
    real(r8) :: fn (num_nolakec,-nlevsno+1:nlevgrnd)                        ! heat diffusion through the layer interface [W/m2]
    real(r8) :: fn1(num_nolakec,-nlevsno+1:nlevgrnd)                        ! heat diffusion through the layer interface [W/m2]
    real(r8) :: sabg_lyr_col(num_nolakec,-nlevsno+1:1)                      ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) :: tvector(num_nolakec,-nlevsno:nlevgrnd)                      ! initial temperature solution [Kelvin]
    real(r8) :: hs_top(num_nolakec)                                         ! net energy flux into surface layer (col) [W/m2]
    real(r8) :: dz_h2osfc(num_nolakec)                                      ! height of standing surface water [m]
    real(r8) :: tk_h2osfc(num_nolakec)                                      ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8) :: dhsdT(num_nolakec)                                          ! temperature derivative of "hs" [col]
    real(r8) :: hs_soil(num_nolakec)                                        ! heat flux on soil [W/m2]
    real(r8) :: hs_top_snow(num_nolakec)                                    ! heat flux on top snow layer [W/m2]
    real(r8) :: hs_h2osfc(num_nolakec)                                      ! heat flux on standing water [W/m2]
    logical  :: cool_on(bounds%begl:bounds%endl)                            ! is urban air conditioning on?
    logical  :: heat_on(bounds%begl:bounds%endl)                            ! is urban heating on?
    !-----------------------------------------------------------------------

    associate(                                                                   &
         snl                     => col_pp%snl                                 , & ! Input:  [integer  (:)   ]  number of snow layers
         zi                      => col_pp%zi                                  , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         z                       => col_pp%z                                   , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)

         t_building_max          => urbanparams_vars%t_building_max         , & ! Input:  [real(r8) (:)   ]  maximum internal building temperature (K)
         t_building_min          => urbanparams_vars%t_building_min         , & ! Input:  [real(r8) (:)   ]  minimum internal building temperature (K)

         frac_sno_eff            => col_ws%frac_sno_eff        , & ! Input:  [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)
         h2osfc                  => col_ws%h2osfc              , & ! Input:  [real(r8) (:)   ]  surface water (mm)
         frac_h2osfc             => col_ws%frac_h2osfc         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)

         eflx_fgr12              => col_ef%eflx_fgr12          , & ! Output: [real(r8) (:)   ]  heat flux between soil layer 1 and 2 (W/m2)
         eflx_fgr                => col_ef%eflx_fgr            , & ! Output: [real(r8) (:,:) ]  (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
         eflx_building_heat      => col_ef%eflx_building_heat  , & ! Output: [real(r8) (:)   ]  heat flux from urban building interior to walls, roof (W/m**2)
         eflx_urban_ac           => col_ef%eflx_urban_ac       , & ! Output: [real(r8) (:)   ]  urban air conditioning flux (W/m**2)
         eflx_urban_heat         => col_ef%eflx_urban_heat     , & ! Output: [real(r8) (:)   ]  urban heating flux (W/m**2)

         hc_soi                  => col_es%hc_soi                           , & ! Output: [real(r8) (:)   ]  soil heat content (MJ/m2)
         hc_soisno               => col_es%hc_soisno                        , & ! Output: [real(r8) (:)   ]  soil plus snow plus lake heat content (MJ/m2)
         hc_lake                 => col_es%hc_lake                          , & ! Output: [real(r8) (:)   ]  lake-water heat content (MJ/m2; zero for non-lake columns)
         tssbef                  => col_es%t_ssbef                          , & ! Input:  [real(r8) (:,:) ]  temperature at previous time step [K]
         t_h2osfc                => col_es%t_h2osfc                         , & ! Output: [real(r8) (:)   ]  surface water temperature
         t_soisno                => col_es%t_soisno                         , & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_grnd                  => col_es%t_grnd                           , & ! Output: [real(r8) (:)   ]  ground surface temperature [K]
         t_building              => lun_es%t_building                       , & ! Output: [real(r8) (:)   ]  internal building temperature (K)
         xmf_h2osfc              => col_ef%xmf_h2osfc                       , & ! Output: [real(r8) (:)   ] latent heat of phase change of surface water [col]
         c_h2osfc                => col_es%c_h2osfc                           & ! Output: [real(r8) (:)   ] heat capacity of surface water [col]
         )

      ! Get step size

      dtime = dtime_mod !get_step_size()

      !$acc data create(jtop, jbot, cv, tk, fn, fn1, sabg_lyr_col, tvector, hs_top, &
      !$acc    dz_h2osfc, tk_h2osfc, dhsdT, hs_soil, hs_top_snow, hs_h2osfc, cool_on, heat_on)

      ! Restrict internal building temperature to between min and max
      ! and determine if heating or air conditioning is on
      !$acc parallel loop gang vector default(present) private(l)
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         if (lun_pp%urbpoi(l)) then
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

      ! set up compact matrix for band diagonal solver, requires additional
      !     sub/super diagonals (1 each), and one additional row for t_h2osfc
      !$acc parallel loop gang vector default(present) private(c)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         jtop(fc) = snl(c)
         ! compute jbot
         if ((col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
              .or. col_pp%itype(c) == icol_roof) ) then
            jbot(fc) = nlevurb
         else
            jbot(fc) = nlevgrnd
         endif
      end do

      !------------------------------------------------------
      ! Compute ground surface and soil temperatures
      !------------------------------------------------------

      ! Thermal conductivity and Heat capacity

      call SoilThermProp(num_nolakec, filter_nolakec, &
           tk, cv, tk_h2osfc, &
           urbanparams_vars, soilstate_vars)

      ! Net ground heat flux into the surface and its temperature derivative
      ! Added a patches loop here to get the average of hs and dhsdT over
      ! all Patches on the column. Precalculate the terms that do not depend on PFT.

      call t_startf('SoilTempGrndFlux')
      call ComputeGroundHeatFluxAndDeriv(num_nolakec, filter_nolakec, &
           hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, &
           canopystate_vars, solarabs_vars)
      call t_stopf('SoilTempGrndFlux')

      ! Determine heat diffusion through the layer interface and factor used in computing
      ! banded diagonal matrix and set up vector r and vectors a, b, c that define banded
      ! diagonal matrix and solve system

      call t_startf('SoilTempHeatDiff')
      call ComputeHeatDiffFluxAndFactor(num_nolakec, filter_nolakec, dtime, &
           tk, cv, fn)
      call t_stopf('SoilTempHeatDiff')

      ! compute thermal properties of h2osfc

      !$acc parallel loop gang vector default(present) private(c)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         if ( (h2osfc(c) > thin_sfclayer) .and. (frac_h2osfc(c) > thin_sfclayer) ) then
            c_h2osfc(c)   = max(thin_sfclayer, cpliq*h2osfc(c)/frac_h2osfc(c)  )
            dz_h2osfc(fc) = max(thin_sfclayer, 1.0e-3*h2osfc(c)/frac_h2osfc(c) )
         else
            c_h2osfc(c)   = thin_sfclayer
            dz_h2osfc(fc) = thin_sfclayer
         endif
      enddo

      ! initialize initial temperature vector

      call t_startf('SoilTempInit')
      !$acc parallel loop gang vector default(present) private(c)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         !$acc loop seq
         do j = -nlevsno, nlevgrnd
            tvector(fc,j) = spval
         end do
         !$acc loop seq
         do j = snl(c)+1, 0
            tvector(fc,j-1) = t_soisno(c,j)
         end do
         ! surface water layer has two coefficients
         tvector(fc,0) = t_h2osfc(c)

         ! soil layers; top layer will have one offset and one extra coefficient
         !$acc loop seq
         do j = 1, nlevgrnd
            tvector(fc,j) = t_soisno(c,j)
         end do
      enddo
      call t_stopf('SoilTempInit')

      ! Solve temperature

      call SolveTemperature(num_nolakec, filter_nolakec, dtime, &
           hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, &
           tk, tk_h2osfc, fn, dz_h2osfc, jtop, jbot, tvector)

      ! return temperatures to original array

      !$acc parallel loop gang vector default(present) private(c)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         !$acc loop seq
         do j = snl(c)+1, 0
            t_soisno(c,j) = tvector(fc,j-1)        !snow layers
         end do
         !$acc loop seq
         do j = 1, nlevgrnd
            t_soisno(c,j) = tvector(fc,j)          !soil layers
         end do

         if (frac_h2osfc(c) == 0._r8) then
            t_h2osfc(c)   = t_soisno(c,1)
         else
            t_h2osfc(c)   = tvector(fc,0)          !surface water
         endif
      enddo

      call t_startf('SoilTempFn1')
      ! Melting or Freezing

      !$acc parallel loop gang vector collapse(2) default(present) private(c,l)
      do j = -nlevsno+1,nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col_pp%landunit(c)
            if ((col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) .and. j <= nlevurb) then
               if (j >= snl(c)+1) then
                  if (j <= nlevurb-1) then
                     fn1(fc,j) = tk(fc,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j == nlevurb) then
                     ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                     ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                     ! building temperature. (See Oleson urban notes of 6/18/03).
                     ! Note new formulation for fn, this will be used below in net energey flux computations
                     fn1(fc,j) = tk(fc,j) * (t_building(l) - t_soisno(c,j))/(zi(c,j) - z(c,j))
                     fn(fc,j)  = tk(fc,j) * (t_building(l) - tssbef(c,j))/(zi(c,j) - z(c,j))
                  end if
               end if
            else if (col_pp%itype(c) /= icol_sunwall .and. col_pp%itype(c) /= icol_shadewall &
                 .and. col_pp%itype(c) /= icol_roof) then
               if (j >= snl(c)+1) then
                  if (j <= nlevgrnd-1) then
                     fn1(fc,j) = tk(fc,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j == nlevgrnd) then
                     fn1(fc,j) = 0._r8
                  end if
               end if
            end if
         end do
      end do

      !$acc parallel loop gang vector default(present) private(c,l)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col_pp%landunit(c)
         if (lun_pp%urbpoi(l)) then
            if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall .or. col_pp%itype(c) == icol_roof) then
               eflx_building_heat(c) = cnfac*fn(fc,nlevurb) + (1._r8-cnfac)*fn1(fc,nlevurb)
            else
               eflx_building_heat(c) = 0._r8
            end if
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

      call t_stopf('SoilTempFn1')

      ! compute phase change of h2osfc

      !$acc parallel loop gang vector default(present) private(c)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         xmf_h2osfc(c) = 0.
      end do

      call t_startf('PhaseChangeH2osfc')
      call PhaseChangeH2osfc (num_nolakec, filter_nolakec, dhsdT, dtime)
      call t_stopf('PhaseChangeH2osfc')

      call Phasechange_beta (num_nolakec, filter_nolakec, dhsdT, soilstate_vars, dtime)

      !$acc parallel loop gang vector default(present) private(c)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         ! this expression will (should) work whether there is snow or not
         if (snl(c) < 0) then
            if(frac_h2osfc(c) /= 0._r8) then
               t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                    + (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) &
                    + frac_h2osfc(c) * t_h2osfc(c)
            else
               t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                    + (1.0_r8 - frac_sno_eff(c)) * t_soisno(c,1)
            end if

         else
            if(frac_h2osfc(c) /= 0._r8) then
               t_grnd(c) = (1._r8 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
            else
               t_grnd(c) = t_soisno(c,1)
            end if
         endif
      end do

      ! Initialize soil heat content, then calculate soil heat content and soil
      ! plus snow heat content. The layer sums run sequentially within each
      ! column so they accumulate in the same order as a serial build.

      !$acc parallel loop gang vector default(present) private(c,l)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col_pp%landunit(c)
         if (.not. lun_pp%urbpoi(l)) then
            hc_soisno(c) = 0._r8
            hc_soi(c)    = 0._r8
            hc_lake(c)   = 0._r8
         end if
         eflx_fgr12(c)= 0._r8

         !$acc loop seq
         do j = -nlevsno+1,nlevgrnd
            if (j == 1) then ! this only needs to be done once
               eflx_fgr12(c) = -cnfac*fn(fc,1) - (1._r8-cnfac)*fn1(fc,1)
            end if
            if (j > 0 .and. j < nlevgrnd .and. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) then
               eflx_fgr(c,j) = -cnfac*fn(fc,j) - (1._r8-cnfac)*fn1(fc,j)
            else if (j == nlevgrnd .and. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) then
               eflx_fgr(c,j) = 0._r8
            end if

            if (.not. lun_pp%urbpoi(l)) then
               if (j >= snl(c)+1) then
                  hc_soisno(c) = hc_soisno(c) + cv(fc,j)*t_soisno(c,j) / 1.e6_r8
               endif
               if (j >= 1) then
                  hc_soi(c) = hc_soi(c) + cv(fc,j)*t_soisno(c,j) / 1.e6_r8
               end if
            end if
         end do
      end do

      !$acc end data

    end associate

  end subroutine SoilTemperature

  !-----------------------------------------------------------------------
  subroutine SolveTemperature(num_filter, filter, dtime, &
       hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, tk, &
       tk_h2osfc, fn, dz_h2osfc, jtop, jbot, tvector)
    !
    ! !DESCRIPTION:
    !  Assembles and solves the banded penta-diagonal system of equations
    !
    ! !USES:
    use BandDiagonalMod  , only : BandDiagonal
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                               ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                ! column filter
    real(r8) , intent(in)    :: dtime                                    ! land model time step (sec)
    real(r8) , intent(in)    :: hs_h2osfc(num_filter)                    ! heat flux on standing water [W/m2]
    real(r8) , intent(in)    :: hs_top_snow(num_filter)                  ! heat flux on top snow layer [W/m2]
    real(r8) , intent(in)    :: hs_soil(num_filter)                      ! heat flux on soil [W/m2]
    real(r8) , intent(in)    :: hs_top(num_filter)                       ! net energy flux into surface layer (col) [W/m2]
    real(r8) , intent(in)    :: dhsdT(num_filter)                        ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: sabg_lyr_col(num_filter,-nlevsno+1:1)    ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) , intent(in)    :: tk(num_filter,-nlevsno+1:nlevgrnd)       ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)    :: tk_h2osfc(num_filter)                    ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8) , intent(in)    :: fn(num_filter,-nlevsno+1:nlevgrnd)       ! heat diffusion through the layer interface [W/m2]
    real(r8) , intent(in)    :: dz_h2osfc(num_filter)                    ! Thickness of standing water [m]
    integer  , intent(in)    :: jtop(num_filter)                         ! top level at each column
    integer  , intent(in)    :: jbot(num_filter)                         ! bottom level at each column
    real(r8) , intent(inout) :: tvector(num_filter,-nlevsno:nlevgrnd)    ! Numerical solution of temperature
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, j, k
    real(r8) :: bmatrix(num_filter,nband,-nlevsno:nlevgrnd)              ! banded matrix for numerical solution of temperature
    real(r8) :: rvector(num_filter,-nlevsno:nlevgrnd)                    ! RHS vector for numerical solution of temperature
    !-----------------------------------------------------------------------

    !$acc data create(bmatrix, rvector)

    !$acc parallel loop gang vector collapse(2) default(present)
    do j = -nlevsno, nlevgrnd
       do fc = 1, num_filter
          rvector(fc,j) = spval
       end do
    end do

    !$acc parallel loop gang vector collapse(3) default(present)
    do j = -nlevsno, nlevgrnd
       do k = 1, nband
          do fc = 1, num_filter
             bmatrix(fc,k,j) = 0._r8
          end do
       end do
    end do

    ! Set up the RHS vector

    call t_startf('SoilTempRHS')
    call SetRHSVec(num_filter, filter, dtime, &
         hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, &
         tk_h2osfc, fn, dz_h2osfc, rvector)
    call t_stopf('SoilTempRHS')

    ! Set up the banded diagonal matrix

    call t_startf('SoilTempMatrix')
    call SetMatrix(num_filter, filter, dtime, &
         dhsdT, tk, tk_h2osfc, dz_h2osfc, bmatrix)
    call t_stopf('SoilTempMatrix')

    ! Solve the system
    call t_startf('SoilTempBandDiag')
    call BandDiagonal(num_filter, -nlevsno, nlevgrnd, jtop, jbot, &
         bmatrix, rvector, tvector)
    call t_stopf('SoilTempBandDiag')

    !$acc end data

  end subroutine SolveTemperature

  !-----------------------------------------------------------------------
  subroutine SoilThermProp (num_nolakec, filter_nolakec, &
       tk, cv, tk_h2osfc, &
       urbanparams_vars, soilstate_vars)
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
    use elm_varcon      , only : denh2o, denice, tfrz, tkwat, tkice, tkair, cpice,  cpliq, thk_bedrock
    use landunit_varcon , only : istice, istice_mec, istwet
    use elm_varctl      , only : use_T_rho_dependent_snowthk
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_nolakec                           ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec(:)                     ! column filter for non-lake points
    real(r8)               , intent(out)   :: tk(num_nolakec,-nlevsno+1:nlevgrnd)   ! thermal conductivity at the layer interface [W/(m K) ]
    real(r8)               , intent(out)   :: cv(num_nolakec,-nlevsno+1:nlevgrnd)   ! heat capacity [J/(m2 K)]
    real(r8)               , intent(out)   :: tk_h2osfc(num_nolakec)                ! thermal conductivity of h2osfc [W/(m K)]
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: l,c,j,i                   ! indices
    integer  :: nlevbed                   ! # levels to bedrock
    integer  :: fc                        ! lake filtered column indices
    real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
    real(r8) :: dke                       ! kersten number
    real(r8) :: fl                        ! volume fraction of liquid or unfrozen water to total water
    real(r8) :: satw                      ! relative total water content of soil.
    real(r8) :: zh2osfc
    real(r8), parameter :: rho_ice     = 917._r8
    real(r8), parameter :: k_snw_tmps(5) = [ 223.0_r8, 248.0_r8, 263.0_r8, 268.0_r8, 273.0_r8 ]
    real(r8), parameter :: k_snw_coe1(5) = [ 2.564_r8, 2.172_r8, 1.985_r8, 1.883_r8, 1.776_r8 ]
    real(r8), parameter :: k_snw_coe2(5) = [ -0.059_r8, 0.015_r8, 0.073_r8, 0.107_r8, 0.147_r8 ]
    real(r8), parameter :: k_snw_coe3(5) = [ 0.0205_r8, 0.0252_r8, 0.0336_r8, 0.0386_r8, 0.0455_r8 ]
    real(r8) :: k_snw_vals(5)
    !-----------------------------------------------------------------------
    call t_startf('SoilThermProp')

    associate(                                                 &
         snl          =>    col_pp%snl                       , & ! Input:  [integer  (:)   ]  number of snow layers
         dz           =>    col_pp%dz                        , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         zi           =>    col_pp%zi                        , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         z            =>    col_pp%z                         , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         nlev2bed     =>    col_pp%nlevbed                   , & ! Input:  [integer  (:)   ]  number of layers to bedrock

         nlev_improad =>    urbanparams_vars%nlev_improad    , & ! Input:  [integer  (:)   ]  number of impervious road layers
         tk_wall      =>    urbanparams_vars%tk_wall         , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban wall
         tk_roof      =>    urbanparams_vars%tk_roof         , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban roof
         tk_improad   =>    urbanparams_vars%tk_improad      , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban impervious road
         cv_wall      =>    urbanparams_vars%cv_wall         , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban wall
         cv_roof      =>    urbanparams_vars%cv_roof         , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban roof
         cv_improad   =>    urbanparams_vars%cv_improad      , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban impervious road

         t_soisno     =>    col_es%t_soisno                  , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)

         frac_sno     =>    col_ws%frac_sno_eff              , & ! Input:  [real(r8) (:)   ]  fractional snow covered area
         h2osfc       =>    col_ws%h2osfc                    , & ! Input:  [real(r8) (:)   ]  surface (mm H2O)
         h2osno       =>    col_ws%h2osno                    , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
         h2osoi_liq   =>    col_ws%h2osoi_liq                , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice   =>    col_ws%h2osoi_ice                , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         bw           =>    col_ws%bw                        , & ! Output: [real(r8) (:,:) ]  partial density of water in the snow pack (ice + liquid) [kg/m3]

         tkmg         =>    soilstate_vars%tkmg_col          , & ! Input:  [real(r8) (:,:) ]  thermal conductivity, soil minerals  [W/m-K]
         tkdry        =>    soilstate_vars%tkdry_col         , & ! Input:  [real(r8) (:,:) ]  thermal conductivity, dry soil (W/m/Kelvin)
         csol         =>    soilstate_vars%csol_col          , & ! Input:  [real(r8) (:,:) ]  heat capacity, soil solids (J/m**3/Kelvin)
         watsat       =>    soilstate_vars%watsat_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         thk          =>    soilstate_vars%thk_col             & ! Output: [real(r8) (:,:) ]  thermal conductivity of each layer  [W/m-K]
         )

      ! Thermal conductivity of soil from Farouki (1981)
      ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
      ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
      ! 41, 1011-1026.

      !$acc parallel loop gang vector collapse(2) default(present) &
      !$acc    private(c, l, nlevbed, satw, dke, fl, dksat)
      do j = 1,nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            nlevbed = nlev2bed(c)
            l = col_pp%landunit(c)
            if ((col_pp%itype(c) == icol_sunwall .OR. col_pp%itype(c) == icol_shadewall) .and. j <= nlevurb) then
               thk(c,j) = tk_wall(l,j)
            else if (col_pp%itype(c) == icol_roof .and. j <= nlevurb) then
               thk(c,j) = tk_roof(l,j)
            else if (col_pp%itype(c) == icol_road_imperv .and. j >= 1 .and. j <= nlev_improad(l)) then
               thk(c,j) = tk_improad(l,j)
            else if (lun_pp%itype(l) /= istwet .AND. lun_pp%itype(l) /= istice .AND. lun_pp%itype(l) /= istice_mec &
                 .AND. col_pp%itype(c) /= icol_sunwall .AND. col_pp%itype(c) /= icol_shadewall .AND. &
                 col_pp%itype(c) /= icol_roof) then

               satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
               satw = min(1._r8, satw)
               if (satw > .1e-6_r8) then
                  if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                     dke = max(0._r8, log10(satw) + 1.0_r8)
                  else                               ! Frozen soil
                     dke = satw
                  end if
                  fl = (h2osoi_liq(c,j)/(denh2o*dz(c,j))) / (h2osoi_liq(c,j)/(denh2o*dz(c,j)) + &
                       h2osoi_ice(c,j)/(denice*dz(c,j)))
                  dksat = tkmg(c,j)*tkwat**(fl*watsat(c,j))*tkice**((1._r8-fl)*watsat(c,j))
                  thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
               else
                  thk(c,j) = tkdry(c,j)
               endif
               if (j > nlevbed) thk(c,j) = thk_bedrock
            else if (lun_pp%itype(l) == istice .OR. lun_pp%itype(l) == istice_mec) then
               thk(c,j) = tkwat
               if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
            else if (lun_pp%itype(l) == istwet) then
               if (j > nlevbed) then
                  thk(c,j) = thk_bedrock
               else
                  thk(c,j) = tkwat
                  if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
               endif
            endif
         end do
      end do

      ! Thermal conductivity of snow
      if (use_T_rho_dependent_snowthk) then ! choose which snow thermal conductivity to use
         !$acc parallel loop gang vector collapse(2) default(present) private(c, i, k_snw_vals)
         do j = -nlevsno+1,0
            do fc = 1, num_nolakec
               c = filter_nolakec(fc)
               if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then
                  bw(c,j) = (h2osoi_ice(c,j) + h2osoi_liq(c,j)) / (frac_sno(c) * dz(c,j))

                  !$acc loop seq
                  do i = 1, 5
                     k_snw_vals(i) = k_snw_coe1(i) * (bw(c,j) / rho_ice)**2 - k_snw_coe2(i) * (bw(c,j) / rho_ice) + k_snw_coe3(i)
                  end do

                  !$acc loop seq
                  do i = 1, size(k_snw_tmps) - 1
                     if (k_snw_tmps(i) <= t_soisno(c,j) .and. t_soisno(c,j) <= k_snw_tmps(i + 1)) then
                        thk(c,j) = k_snw_vals(i) + (t_soisno(c,j) - k_snw_tmps(i)) * (k_snw_vals(i + 1)-k_snw_vals(i)) / (k_snw_tmps(i + 1) - k_snw_tmps(i))
                     end if
                  end do

                  ! Handle edge cases if t_soisno(c,j) is outside the given range
                  if (t_soisno(c,j) < k_snw_tmps(1)) then
                     thk(c,j) = k_snw_vals(1)
                  else if (t_soisno(c,j) > k_snw_tmps(size(k_snw_tmps))) then
                     thk(c,j) = k_snw_vals(size(k_snw_tmps))
                  end if
               end if
            end do
         end do
      else
         !$acc parallel loop gang vector collapse(2) default(present) private(c)
         do j = -nlevsno+1,0
            do fc = 1, num_nolakec
               c = filter_nolakec(fc)
               ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
               ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
               if (snl(c) + 1 < 1 .AND. (j >= snl(c) + 1) .AND. (j <= 0)) then
                  bw(c,j) = (h2osoi_ice(c,j) + h2osoi_liq(c,j)) / (frac_sno(c) * dz(c,j))
                  thk(c,j) = tkair + (7.75e-5_r8 * bw(c,j) + 1.105e-6_r8 * bw(c,j) * bw(c,j)) * (tkice - tkair)
               end if
            end do
         end do
      end if

      ! Thermal conductivity at the layer interface

      !$acc parallel loop gang vector collapse(2) default(present) private(c)
      do j = -nlevsno+1,nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if ((col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) .and. j <= nlevurb) then
               if (j >= snl(c)+1 .AND. j <= nlevurb-1) then
                  tk(fc,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                       /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
               else if (j == nlevurb) then

                  ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                  ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                  ! building temperature. (See Oleson urban notes of 6/18/03).
                  tk(fc,j) = thk(c,j)
               end if
            else if (col_pp%itype(c) /= icol_sunwall .and. col_pp%itype(c) /= icol_shadewall &
                 .and. col_pp%itype(c) /= icol_roof) then
               if (j >= snl(c)+1 .AND. j <= nlevgrnd-1) then
                  tk(fc,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                       /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
               else if (j == nlevgrnd) then
                  tk(fc,j) = 0._r8
               end if
            end if
         end do
      end do

      ! calculate thermal conductivity of h2osfc

      !$acc parallel loop gang vector default(present) private(c, zh2osfc)
      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         zh2osfc=1.0e-3*(0.5*h2osfc(c)) !convert to [m] from [mm]
         tk_h2osfc(fc)= tkwat*thk(c,1)*(z(c,1)+zh2osfc) &
              /(tkwat*z(c,1)+thk(c,1)*zh2osfc)
      enddo

      ! Soil heat capacity, from de Vires (1963)
      ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
      ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
      ! 41, 1011-1026.

      !$acc parallel loop gang vector collapse(2) default(present) private(c, l, nlevbed)
      do j = 1, nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col_pp%landunit(c)
            nlevbed = nlev2bed(c)
            if ((col_pp%itype(c) == icol_sunwall .OR. col_pp%itype(c) == icol_shadewall) .and. j <= nlevurb) then
               cv(fc,j) = cv_wall(l,j) * dz(c,j)
            else if (col_pp%itype(c) == icol_roof .and. j <= nlevurb) then
               cv(fc,j) = cv_roof(l,j) * dz(c,j)
            else if (col_pp%itype(c) == icol_road_imperv .and. j >= 1 .and. j <= nlev_improad(l)) then
               cv(fc,j) = cv_improad(l,j) * dz(c,j)
            else if (lun_pp%itype(l) /= istwet .AND. lun_pp%itype(l) /= istice .AND. lun_pp%itype(l) /= istice_mec &
                 .AND. col_pp%itype(c) /= icol_sunwall .AND. col_pp%itype(c) /= icol_shadewall .AND. &
                 col_pp%itype(c) /= icol_roof) then
               cv(fc,j) = csol(c,j)*(1._r8-watsat(c,j))*dz(c,j) + (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
               ! Port of CLM bedrock heat-capacity fix: do not reduce solid-rock
               ! heat capacity by soil porosity below the soil-bedrock boundary.
               ! csol is already set to bedrock mineral heat capacity for these layers.
               if (j > nlevbed) cv(fc,j) = csol(c,j)*dz(c,j)
            else if (lun_pp%itype(l) == istwet) then
               cv(fc,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
               if (j > nlevbed) cv(fc,j) = csol(c,j)*dz(c,j)
            else if (lun_pp%itype(l) == istice .OR. lun_pp%itype(l) == istice_mec) then
               cv(fc,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
            endif
            if (j == 1) then
               if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
                  cv(fc,j) = cv(fc,j) + cpice*h2osno(c)
               end if
            end if
         enddo
      end do

      ! Snow heat capacity

      !$acc parallel loop gang vector collapse(2) default(present) private(c)
      do j = -nlevsno+1,0
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
               if (frac_sno(c) > 0._r8) then
                  cv(fc,j) = max(thin_sfclayer, (cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j))/frac_sno(c))
               else
                  cv(fc,j) = thin_sfclayer
               end if
            end if
         end do
      end do

    end associate

    call t_stopf('SoilThermProp')

  end subroutine SoilThermProp

  !-----------------------------------------------------------------------
  subroutine PhaseChangeH2osfc (num_nolakec, filter_nolakec, dhsdT, dtime)
    !
    ! !DESCRIPTION:
    ! Only freezing is considered.  When water freezes, move ice to bottom snow layer.
    !
    ! !USES:
    use elm_varcon       , only : tfrz, hfus, denice, cpliq
    !
    ! !ARGUMENTS:
    integer  , intent(in) :: num_nolakec                   ! number of column non-lake points in column filter
    integer  , intent(in) :: filter_nolakec(:)             ! column filter for non-lake points
    real(r8) , intent(in) :: dhsdT(num_nolakec)            ! temperature derivative of "hs" [col]
    real(r8) , intent(in) :: dtime                         ! land model time step (sec)
    !
    ! !LOCAL VARIABLES:
    integer  :: c         !do loop index
    integer  :: fc        !lake filtered column indices
    real(r8) :: temp1     !temporary variables [kg/m2]
    real(r8) :: hm        !energy residual [W/m2]
    real(r8) :: xm        !melting or freezing within a time step [kg/m2]
    real(r8) :: tinc      !t(n+1)-t(n) (K)
    real(r8) :: rho_avg
    real(r8) :: z_avg
    real(r8) :: c1
    real(r8) :: c2
    !-----------------------------------------------------------------------
    call t_startf('PhaseChangeH2osfc')

    associate(                                                                   &
         snl                       =>    col_pp%snl                , & ! Input:  [integer  (:)   ] number of snow layers

         frac_sno                  =>    col_ws%frac_sno_eff       , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         frac_h2osfc               =>    col_ws%frac_h2osfc        , & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         h2osno                    =>    col_ws%h2osno             , & ! Output: [real(r8) (:)   ] snow water (mm H2O)
         h2osoi_ice                =>    col_ws%h2osoi_ice         , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)
         h2osfc                    =>    col_ws%h2osfc             , & ! Output: [real(r8) (:)   ] surface water (mm)
         int_snow                  =>    col_ws%int_snow           , & ! Output: [real(r8) (:)   ] integrated snowfall [mm]
         snow_depth                =>    col_ws%snow_depth         , & ! Output: [real(r8) (:)   ] snow height (m)

         qflx_h2osfc_to_ice        =>    col_wf%qflx_h2osfc_to_ice , & ! Output: [real(r8) (:)   ] conversion of h2osfc to ice

         eflx_h2osfc_to_snow_col   =>    col_ef%eflx_h2osfc_to_snow, & ! Output: [real(r8) (:)   ] col snow melt to h2osfc heat flux (W/m**2)

         fact                      =>    col_es%fact               , & ! Input:  [real(r8) (:,:) ] used in computing tridiagonal matrix
         c_h2osfc                  =>    col_es%c_h2osfc           , & ! Input:  [real(r8) (:)   ] heat capacity of surface water
         xmf_h2osfc                =>    col_ef%xmf_h2osfc         , & ! Output: [real(r8) (:)   ] latent heat of phase change of surface water
         t_soisno                  =>    col_es%t_soisno           , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)
         t_h2osfc                  =>    col_es%t_h2osfc             & ! Output: [real(r8) (:)   ] surface water temperature
         )

      ! Freezing identification

      !$acc parallel loop gang vector default(present) &
      !$acc    private(c, temp1, hm, xm, tinc, rho_avg, z_avg, c1, c2)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)

         ! Initialization
         xmf_h2osfc(c)              = 0._r8
         hm                         = 0._r8
         xm                         = 0._r8
         qflx_h2osfc_to_ice(c)      = 0._r8
         eflx_h2osfc_to_snow_col(c) = 0._r8

         ! If liquid exists below melt point, freeze some to ice.
         if ( frac_h2osfc(c) > 0._r8 .AND. t_h2osfc(c) <= tfrz) then
            tinc = tfrz - t_h2osfc(c)
            t_h2osfc(c) = tfrz

            ! energy absorbed beyond freezing temperature
            hm = frac_h2osfc(c)*(dhsdT(fc)*tinc - tinc*c_h2osfc(c)/dtime)

            ! mass of water converted from liquid to ice
            xm = hm*dtime/hfus
            temp1 = h2osfc(c) + xm

            z_avg=frac_sno(c)*snow_depth(c)
            if (z_avg > 0._r8) then
               rho_avg=min(800._r8,h2osno(c)/z_avg)
            else
               rho_avg=200._r8
            endif

            !=====================  xm < h2osfc  ====================================
            if(temp1 >= 0._r8) then ! add some frozen water to snow column

               ! add ice to snow column
               h2osno(c) = h2osno(c) - xm
               int_snow(c) = int_snow(c) - xm

               if(snl(c) < 0) h2osoi_ice(c,0) = h2osoi_ice(c,0) - xm

               ! remove ice from h2osfc
               h2osfc(c) = h2osfc(c) + xm

               xmf_h2osfc(c) = hm
               qflx_h2osfc_to_ice(c) = -xm/dtime

               ! update snow depth
               if (frac_sno(c) > 0 .and. snl(c) < 0) then
                  snow_depth(c)=h2osno(c)/(rho_avg*frac_sno(c))
               else
                  snow_depth(c)=h2osno(c)/denice
               endif

               ! adjust temperature of lowest snow layer to account for addition of ice
               if (snl(c) == 0) then
                  !initialize for next time step
                  t_soisno(c,0) = t_h2osfc(c)
                  eflx_h2osfc_to_snow_col(c) = 0._r8
               else
                  if (snl(c) == -1)then
                     c1=frac_sno(c)*(dtime/fact(c,0) - dhsdT(fc)*dtime)
                  else
                     c1=frac_sno(c)/fact(c,0)*dtime
                  end if
                  if ( frac_h2osfc(c) /= 0.0_r8 )then
                     c2=(-cpliq*xm - frac_h2osfc(c)*dhsdT(fc)*dtime)
                  else
                     c2=0.0_r8
                  end if
                  t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc(c)) &
                       /(c1 + c2)
                  eflx_h2osfc_to_snow_col(c) =(t_h2osfc(c)-t_soisno(c,0))*c2/dtime
               endif

               !=========================  xm > h2osfc  =============================
            else !all h2osfc converted to ice

               rho_avg=(h2osno(c)*rho_avg + h2osfc(c)*denice)/(h2osno(c) + h2osfc(c))
               h2osno(c) = h2osno(c) + h2osfc(c)
               int_snow(c) = int_snow(c) + h2osfc(c)

               qflx_h2osfc_to_ice(c) = h2osfc(c)/dtime

               ! excess energy is used to cool ice layer
               if(snl(c) < 0) h2osoi_ice(c,0) = h2osoi_ice(c,0) + h2osfc(c)

               ! NOTE: should compute and then use the heat capacity of frozen h2osfc layer
               !       rather than using heat capacity of the liquid layer. But this causes
               !       balance check errors as it doesn't know about it.
               ! compute heat capacity of frozen h2osfc layer

               ! cool frozen h2osfc layer with extra heat
               t_h2osfc(c) = t_h2osfc(c) - temp1*hfus/(dtime*dhsdT(fc) - c_h2osfc(c))

               xmf_h2osfc(c) = (hm - frac_h2osfc(c)*temp1*hfus/dtime)

               ! next, determine equilibrium temperature of combined ice/snow layer
               if (snl(c) == 0) then
                  !initialize for next time step
                  t_soisno(c,0) = t_h2osfc(c)
               else if (snl(c) == -1) then
                  c1=frac_sno(c)*(dtime/fact(c,0) - dhsdT(fc)*dtime)

                  if ( frac_h2osfc(c) /= 0.0_r8 )then
                     c2=frac_h2osfc(c)*(c_h2osfc(c) - dtime*dhsdT(fc))
                  else
                     c2=0.0_r8
                  end if
                  ! account for the change in t_soisno(c,0) via xmf_h2osfc(c)
                  t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc(c)) &
                       /(c1 + c2)
                  t_h2osfc(c) = t_soisno(c,0)

               else
                  c1=frac_sno(c)/fact(c,0)*dtime
                  if ( frac_h2osfc(c) /= 0.0_r8 )then
                     c2=frac_h2osfc(c)*(c_h2osfc(c) - dtime*dhsdT(fc))
                  else
                     c2=0.0_r8
                  end if
                  t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc(c)) &
                       /(c1 + c2)
                  t_h2osfc(c) = t_soisno(c,0)
               endif

               ! set h2osfc to zero (all liquid converted to ice)
               h2osfc(c) = 0._r8

               ! update snow depth
               if (frac_sno(c) > 0 .and. snl(c) < 0) then
                  snow_depth(c)=h2osno(c)/(rho_avg*frac_sno(c))
               else
                  snow_depth(c)=h2osno(c)/denice
               endif

            endif
         endif
      enddo

    end associate

    call t_stopf('PhaseChangeH2osfc')

  end subroutine PhaseChangeH2osfc

  !-----------------------------------------------------------------------
  subroutine Phasechange_beta (num_nolakec, filter_nolakec, dhsdT, &
       soilstate_vars, dtime)
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
    use elm_varcon       , only : tfrz, hfus, grav
    use landunit_varcon  , only : istice_mec, istice
    !
    ! !ARGUMENTS:
    integer              , intent(in) :: num_nolakec                  ! number of column non-lake points in column filter
    integer              , intent(in) :: filter_nolakec(:)            ! column filter for non-lake points
    real(r8)             , intent(in) :: dhsdT(num_nolakec)           ! temperature derivative of "hs" [col]
    type(soilstate_type) , intent(in) :: soilstate_vars
    real(r8)             , intent(in) :: dtime                        ! land model time step (sec)
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                               !do loop index
    integer  :: fc                                                  !lake filtered column indices
    real(r8) :: heatr                                               !energy residual or loss after melting or freezing
    real(r8) :: temp1                                               !temporary variables [kg/m2]
    real(r8) :: hm                                                  !energy residual [W/m2]
    real(r8) :: xm                                                  !melting or freezing within a time step [kg/m2]
    real(r8) :: wmass0                                              !initial mass of ice and liquid (kg/m2)
    real(r8) :: wice0                                               !initial mass of ice (kg/m2)
    real(r8) :: propor                                              !proportionality constant (-)
    real(r8) :: smp                                                 !frozen water potential (mm)
    real(r8) :: supercool(num_nolakec,nlevgrnd)                     !supercooled water in soil (kg/m2)
    real(r8) :: tinc(num_nolakec,-nlevsno+1:nlevgrnd)               !t(n+1)-t(n) (K)
    !-----------------------------------------------------------------------
    call t_startf('PhaseChangebeta')

    associate(                                                        &
         snl              =>    col_pp%snl                , & ! Input:  [integer  (:)   ] number of snow layers
         dz               =>    col_pp%dz                 , & ! Input:  [real(r8) (:,:) ] layer thickness (m)

         bsw              =>    soilstate_vars%bsw_col    , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
         sucsat           =>    soilstate_vars%sucsat_col , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
         watsat           =>    soilstate_vars%watsat_col , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)

         frac_sno_eff     =>    col_ws%frac_sno_eff       , & ! Input:  [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         frac_h2osfc      =>    col_ws%frac_h2osfc        , & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         snow_depth       =>    col_ws%snow_depth         , & ! Output: [real(r8) (:)   ] snow height (m)
         h2osno           =>    col_ws%h2osno             , & ! Output: [real(r8) (:)   ] snow water (mm H2O)
         h2osoi_liq       =>    col_ws%h2osoi_liq         , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)
         h2osoi_ice       =>    col_ws%h2osoi_ice         , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)

         qflx_snow_melt   =>    col_wf%qflx_snow_melt     , & ! Output: [real(r8) (:)   ] net snow melt
         qflx_snofrz_lyr  =>    col_wf%qflx_snofrz_lyr    , & ! Output: [real(r8) (:,:) ] snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
         qflx_snofrz_col  =>    col_wf%qflx_snofrz        , & ! Output: [real(r8) (:)   ] column-integrated snow freezing rate (positive definite) [kg m-2 s-1]
         qflx_glcice      =>    col_wf%qflx_glcice        , & ! Output: [real(r8) (:)   ] flux of new glacier ice (mm H2O/s) [+ = ice grows]
         qflx_glcice_melt =>    col_wf%qflx_glcice_melt   , & ! Output: [real(r8) (:)   ] ice melt (positive definite) (mm H2O/s)
         qflx_glcice_diag      =>    col_wf%qflx_glcice_diag      , & ! Output: [real(r8) (:)   ] flux of new glacier ice (mm H2O/s) [+ = ice grows]
         qflx_glcice_melt_diag =>    col_wf%qflx_glcice_melt_diag , & ! Output: [real(r8) (:)   ] ice melt (positive definite) (mm H2O/s)
         qflx_snomelt     =>    col_wf%qflx_snomelt       , & ! Output: [real(r8) (:)   ] snow melt (mm H2O /s)
         qflx_snomelt_lyr =>    col_wf%qflx_snomelt_lyr   , & ! Output: [real(r8) (:,:) ] snow melt in each layer (mm H2O /s)

         eflx_snomelt     =>    col_ef%eflx_snomelt       , & ! Output: [real(r8) (:)   ] snow melt heat flux (W/m**2)
         eflx_snomelt_r   =>    col_ef%eflx_snomelt_r     , & ! Output: [real(r8) (:)   ] rural snow melt heat flux (W/m**2)
         eflx_snomelt_u   =>    col_ef%eflx_snomelt_u     , & ! Output: [real(r8) (:)   ] urban snow melt heat flux (W/m**2)

         xmf              =>    col_ef%xmf                , & ! Output: [real(r8) (:)   ] total latent heat of phase change of ground water
         fact             =>    col_es%fact               , & ! Input:  [real(r8) (:,:) ] used in computing tridiagonal matrix

         imelt            =>    col_ef%imelt              , & ! Output: [integer  (:,:) ] flag for melting (=1), freezing (=2), Not=0 (new)
         t_soisno         =>    col_es%t_soisno             & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)
         )

      !$acc data create(supercool, tinc)

      ! Initialization

      !$acc parallel loop gang vector default(present) private(c)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)

         qflx_snomelt(c) = 0._r8
         !$acc loop seq
         do j = -nlevsno+1,0
            qflx_snomelt_lyr(c,j) = 0._r8
            qflx_snofrz_lyr(c,j)  = 0._r8
         end do
         xmf(c) = 0._r8
         qflx_snofrz_col(c) = 0._r8
         qflx_glcice_melt(c) = 0._r8
         qflx_glcice_melt_diag(c) = 0._r8
         qflx_snow_melt(c) = 0._r8
      end do

      !$acc parallel loop gang vector collapse(2) default(present) private(c)
      do j = -nlevsno+1,nlevgrnd       ! all layers
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (j >= snl(c)+1) then
               imelt(c,j) = 0
            endif   ! end of snow layer if-block
         end do   ! end of column-loop
      enddo   ! end of level-loop

      !--  snow layers  ---------------------------------------------------

      !$acc parallel loop gang vector collapse(2) default(present) private(c)
      do j = -nlevsno+1,0
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (j >= snl(c)+1) then

               ! Melting identification
               ! If ice exists above melt point, melt some to liquid.
               if (h2osoi_ice(c,j) > 0._r8 .AND. t_soisno(c,j) > tfrz) then
                  imelt(c,j) = 1
                  !                tinc(c,j) = t_soisno(c,j) - tfrz
                  tinc(fc,j) = tfrz - t_soisno(c,j)
                  t_soisno(c,j) = tfrz
               endif

               ! Freezing identification
               ! If liquid exists below melt point, freeze some to ice.
               if (h2osoi_liq(c,j) > 0._r8 .AND. t_soisno(c,j) < tfrz) then
                  imelt(c,j) = 2
                  !                tinc(c,j) = t_soisno(c,j) - tfrz
                  tinc(fc,j) = tfrz - t_soisno(c,j)
                  t_soisno(c,j) = tfrz
               endif
            endif   ! end of snow layer if-block
         end do   ! end of column-loop
      enddo   ! end of level-loop

      !-- soil layers   ---------------------------------------------------

      !$acc parallel loop gang vector collapse(2) default(present) private(c, smp)
      do j = 1,nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            supercool(fc,j) = 0.0_r8
            ! add in urban condition if-block
            if ((col_pp%itype(c) /= icol_sunwall .and. col_pp%itype(c) /= icol_shadewall &
                 .and. col_pp%itype(c) /= icol_roof) .or. ( j <= nlevurb)) then

               if (h2osoi_ice(c,j) > 0. .AND. t_soisno(c,j) > tfrz) then
                  imelt(c,j) = 1
                  !             tinc(c,j) = t_soisno(c,j) - tfrz
                  tinc(fc,j) = tfrz - t_soisno(c,j)
                  t_soisno(c,j) = tfrz
               endif

               ! from Zhao (1997) and Koren (1999)
               supercool(fc,j) = 0.0_r8
               if (col_pp%is_soil(c) .or. col_pp%is_crop(c) .or. col_pp%itype(c) == icol_road_perv) then
                  if(t_soisno(c,j) < tfrz) then
                     smp = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                     supercool(fc,j) = watsat(c,j)*(smp/sucsat(c,j))**(-1._r8/bsw(c,j))
                     supercool(fc,j) = supercool(fc,j)*dz(c,j)*1000._r8       ! (mm)
                  endif
               endif

               if (h2osoi_liq(c,j) > supercool(fc,j) .AND. t_soisno(c,j) < tfrz) then
                  imelt(c,j) = 2
                  !             tinc(c,j) = t_soisno(c,j) - tfrz
                  tinc(fc,j) = tfrz - t_soisno(c,j)
                  t_soisno(c,j) = tfrz
               endif

               ! If snow exists, but its thickness is less than the critical value (0.01 m)
               if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. j == 1) then
                  if (t_soisno(c,j) > tfrz) then
                     imelt(c,j) = 1
                     !                tincc,j) = t_soisno(c,j) - tfrz
                     tinc(fc,j) = tfrz - t_soisno(c,j)
                     t_soisno(c,j) = tfrz
                  endif
               endif

            endif

         end do
      enddo

      ! Melting/freezing. Each column walks its layers sequentially so the
      ! column fluxes (xmf, qflx_snomelt, qflx_glcice, ...) accumulate in layer order.

      !$acc parallel loop gang vector default(present) &
      !$acc    private(c, l, j, hm, xm, wmass0, wice0, temp1, propor, heatr)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col_pp%landunit(c)

         !$acc loop seq
         do j = -nlevsno+1,nlevgrnd       ! all layers

            if ((col_pp%itype(c) /= icol_sunwall .and. col_pp%itype(c) /= icol_shadewall &
                 .and. col_pp%itype(c) /= icol_roof) .or. ( j <= nlevurb)) then

               if (j >= snl(c)+1) then

                  ! Initialization
                  hm     = 0._r8
                  xm     = 0._r8
                  wice0  = h2osoi_ice(c,j)
                  wmass0 = h2osoi_ice(c,j) + h2osoi_liq(c,j)

                  ! Calculate the energy surplus and loss for melting and freezing
                  if (imelt(c,j) > 0) then

                     ! added unique cases for this calculation,
                     ! to account for absorbed solar radiation in each layer

                     !==================================================================
                     if (j == snl(c)+1) then ! top layer
                        if(j > 0) then
                           hm = dhsdT(fc)*tinc(fc,j) - tinc(fc,j)/fact(c,j)
                        else
                           hm = frac_sno_eff(c)*(dhsdT(fc)*tinc(fc,j) - tinc(fc,j)/fact(c,j))
                        endif

                        if ( j==1 .and. frac_h2osfc(c) /= 0.0_r8 ) then
                           hm = hm - frac_h2osfc(c)*(dhsdT(fc)*tinc(fc,j))
                        end if
                     else if (j == 1) then
                        hm = (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) &
                             *dhsdT(fc)*tinc(fc,j) - tinc(fc,j)/fact(c,j)
                     else ! non-interfacial snow/soil layers
                        if(j < 1) then
                           hm = - frac_sno_eff(c)*(tinc(fc,j)/fact(c,j))
                        else
                           hm = - tinc(fc,j)/fact(c,j)
                        endif
                     endif
                  endif

                  ! These two errors were checked carefully (Y. Dai).  They result from the
                  ! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
                  if (imelt(c,j) == 1 .AND. hm < 0._r8) then
                     hm = 0._r8
                     imelt(c,j) = 0
                  endif
                  if (imelt(c,j) == 2 .AND. hm > 0._r8) then
                     hm = 0._r8
                     imelt(c,j) = 0
                  endif

                  ! The rate of melting and freezing

                  if (imelt(c,j) > 0 .and. abs(hm) > 0._r8) then
                     xm = hm*dtime/hfus                           ! kg/m2

                     ! If snow exists, but its thickness is less than the critical value
                     ! (1 cm). Note: more work is needed to determine how to tune the
                     ! snow depth for this case
                     if (j == 1) then
                        if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. xm > 0._r8) then
                           temp1 = h2osno(c)                           ! kg/m2
                           h2osno(c) = max(0._r8,temp1-xm)
                           propor = h2osno(c)/temp1
                           snow_depth(c) = propor * snow_depth(c)
                           heatr = hm - hfus*(temp1-h2osno(c))/dtime   ! W/m2
                           if (heatr > 0._r8) then
                              xm = heatr*dtime/hfus                    ! kg/m2
                              hm = heatr                               ! W/m2
                           else
                              xm = 0._r8
                              hm = 0._r8
                           endif
                           qflx_snomelt(c) = max(0._r8,(temp1-h2osno(c)))/dtime   ! kg/(m2 s)
                           xmf(c) = hfus*qflx_snomelt(c)
                           qflx_snow_melt(c) = qflx_snomelt(c)
                        endif
                     endif

                     heatr = 0._r8
                     if (xm > 0._r8) then
                        h2osoi_ice(c,j) = max(0._r8, wice0-xm)
                        heatr = hm - hfus*(wice0-h2osoi_ice(c,j))/dtime
                     else if (xm < 0._r8) then
                        if (j <= 0) then
                           h2osoi_ice(c,j) = min(wmass0, wice0-xm)  ! snow
                        else
                           if (wmass0 < supercool(fc,j)) then
                              h2osoi_ice(c,j) = 0._r8
                           else
                              h2osoi_ice(c,j) = min(wmass0 - supercool(fc,j),wice0-xm)
                           endif
                        endif
                        heatr = hm - hfus*(wice0-h2osoi_ice(c,j))/dtime
                     endif

                     h2osoi_liq(c,j) = max(0._r8,wmass0-h2osoi_ice(c,j))

                     if (abs(heatr) > 0._r8) then
                        if (j == snl(c)+1) then

                           if(j==1) then
                              t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                   /(1._r8-(1.0_r8 - frac_h2osfc(c))*fact(c,j)*dhsdT(fc))
                           else
                              t_soisno(c,j) = t_soisno(c,j) + (fact(c,j)/frac_sno_eff(c))*heatr &
                                   /(1._r8-fact(c,j)*dhsdT(fc))
                           endif

                        else if (j == 1) then

                           t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                /(1._r8-(1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c))*fact(c,j)*dhsdT(fc))
                        else
                           if(j > 0) then
                              t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr
                           else
                              if(frac_sno_eff(c) > 0._r8) t_soisno(c,j) = t_soisno(c,j) + (fact(c,j)/frac_sno_eff(c))*heatr
                           endif
                        endif

                        if (j <= 0) then    ! snow
                           if (h2osoi_liq(c,j)*h2osoi_ice(c,j)>0._r8) t_soisno(c,j) = tfrz
                        end if
                     endif  ! end of heatr > 0 if-block

                     xmf(c) = xmf(c) + hfus*(wice0-h2osoi_ice(c,j))/dtime

                     if (imelt(c,j) == 1 .AND. j < 1) then
                        qflx_snomelt(c) = qflx_snomelt(c) + max(0._r8,(wice0-h2osoi_ice(c,j)))/dtime
                        qflx_snomelt_lyr(c,j) = max(0._r8,(wice0-h2osoi_ice(c,j)))/dtime
                     endif

                     ! layer freezing mass flux (positive):
                     if (imelt(c,j) == 2 .AND. j < 1) then
                        qflx_snofrz_lyr(c,j) = max(0._r8,(h2osoi_ice(c,j)-wice0))/dtime
                     endif

                  endif

               endif   ! end of snow layer if-block

            endif

            ! For glacier_mec columns, compute negative ice flux from melted ice.
            ! Note that qflx_glcice can also include a positive component from excess snow,
            ! as computed in HydrologyDrainageMod.F90.

            if ( lun_pp%itype(l) == istice_mec) then
               if (j>=1 .and. h2osoi_liq(c,j) > 0._r8) then   ! ice layer with meltwater
                  ! melting corresponds to a negative ice flux
                  qflx_glcice_melt(c) = qflx_glcice_melt(c) + h2osoi_liq(c,j)/dtime
                  qflx_glcice(c) = qflx_glcice(c) - h2osoi_liq(c,j)/dtime

                  ! convert layer back to pure ice by "borrowing" ice from below the column
                  h2osoi_ice(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
                  h2osoi_liq(c,j) = 0._r8

               endif  ! liquid water is present
            endif     ! istice_mec
            ! for diagnostic QICE SMB output only -
            ! these are to calculate SMB even without MECs
            if ( lun_pp%itype(l) == istice) then
               if (j>=1 .and. h2osoi_liq(c,j) > 0._r8) then   ! ice layer with meltwater
                  ! melting corresponds to a negative ice flux
                  qflx_glcice_melt_diag(c) = qflx_glcice_melt_diag(c) + h2osoi_liq(c,j)/dtime
                  qflx_glcice_diag(c) = qflx_glcice_diag(c) - h2osoi_liq(c,j)/dtime
               endif  ! liquid water is present
            endif     ! istice

         end do   ! end of level-loop
      enddo   ! end of column-loop

      ! Needed for history file output

      !$acc parallel loop gang vector default(present) private(c, l)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         eflx_snomelt(c) = qflx_snomelt(c) * hfus
         l = col_pp%landunit(c)
         if (lun_pp%urbpoi(l)) then
            eflx_snomelt_u(c) = eflx_snomelt(c)
         else if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then
            eflx_snomelt_r(c) = eflx_snomelt(c)
         end if

         !$acc loop seq
         do j = -nlevsno+1,0
            qflx_snofrz_col(c) = qflx_snofrz_col(c) + qflx_snofrz_lyr(c,j)
         end do
      end do

      !$acc end data

    end associate

    call t_stopf('PhaseChangebeta')

  end subroutine Phasechange_beta

  !-----------------------------------------------------------------------
  subroutine ComputeGroundHeatFluxAndDeriv(num_nolakec, filter_nolakec, &
       hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, &
       canopystate_vars, solarabs_vars)
    !
    ! !DESCRIPTION:
    ! Computes ground heat flux on:
    ! (1) The surface of standing water,
    ! (2) The surface of snow,
    ! (3) The surface of soil, and
    ! (4) Net energy flux into soil surface.
    ! Additionally, derivative of ground heat flux w.r.t to temeprature
    !
    ! !USES:
    use elm_varcon     , only : sb, hvap
    !
    ! !ARGUMENTS:
    implicit none
    integer                , intent(in)    :: num_nolakec                           ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec( : )                   ! column filter for non-lake points
    real(r8)               , intent(out)   :: hs_h2osfc(num_nolakec)                ! heat flux on standing water [W/m2]
    real(r8)               , intent(out)   :: hs_top_snow(num_nolakec)              ! heat flux on top snow layer [W/m2]
    real(r8)               , intent(out)   :: hs_soil(num_nolakec)                  ! heat flux on soil [W/m2]
    real(r8)               , intent(out)   :: hs_top(num_nolakec)                   ! net energy flux into surface layer (col) [W/m2]
    real(r8)               , intent(out)   :: dhsdT(num_nolakec)                    ! temperature derivative of "hs" [col]
    real(r8)               , intent(out)   :: sabg_lyr_col(num_nolakec,-nlevsno+1:1)! absorbed solar radiation (col,lyr) [W/m2]
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,p,l,t,g                                            ! indices
    integer  :: fc                                                     ! lake filtered column indices
    integer  :: lyr_top                                                ! index of top layer of snowpack (-4 to 0) [idx]
    real(r8) :: lwrad_emit                                             ! emitted longwave radiation
    real(r8) :: dlwrad_emit                                            ! time derivative of emitted longwave radiation
    real(r8) :: lwrad_emit_snow                                        !
    real(r8) :: lwrad_emit_soil                                        !
    real(r8) :: lwrad_emit_h2osfc                                      !
    real(r8) :: eflx_gnet_top                                          ! net energy flux into surface layer, pft-level [W/m2]
    real(r8) :: eflx_gnet_snow                                         !
    real(r8) :: eflx_gnet_soil                                         !
    real(r8) :: eflx_gnet_h2osfc                                       !
    real(r8) :: slope_rad
    real(r8), parameter :: deg2rad = SHR_CONST_PI/180._r8
    !-----------------------------------------------------------------------

    associate(                                                                &
         snl                     => col_pp%snl                              , & ! Input:  [integer (:)    ]  number of snow layers

         forc_lwrad              => top_af%lwrad                            , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)

         frac_veg_nosno          => canopystate_vars%frac_veg_nosno_patch   , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]

         frac_sno_eff            => col_ws%frac_sno_eff        , & ! Input:  [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)

         qflx_ev_snow            => veg_wf%qflx_ev_snow       , & ! Input:  [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]
         qflx_ev_soil            => veg_wf%qflx_ev_soil       , & ! Input:  [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]
         qflx_ev_h2osfc          => veg_wf%qflx_ev_h2osfc     , & ! Input:  [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]
         qflx_evap_soi           => veg_wf%qflx_evap_soi      , & ! Input:  [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_tran_veg           => veg_wf%qflx_tran_veg      , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)

         emg                     => col_es%emg                , & ! Input:  [real(r8) (:)   ]  ground emissivity
         t_h2osfc                => col_es%t_h2osfc           , & ! Input:  [real(r8) (:)   ]  surface water temperature
         t_grnd                  => col_es%t_grnd             , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]
         t_soisno                => col_es%t_soisno           , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)

         htvp                    => col_ef%htvp              , & ! Input:  [real(r8) (:)   ]  latent heat of vapor of water (or sublimation) [j/kg]
         cgrnd                   => veg_ef%cgrnd             , & ! Input:  [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]
         dlrad                   => veg_ef%dlrad             , & ! Input:  [real(r8) (:)   ]  downward longwave radiation blow the canopy [W/m2]
         eflx_traffic            => lun_ef%eflx_traffic      , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)
         eflx_wasteheat          => lun_ef%eflx_wasteheat    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_heat_from_ac       => lun_ef%eflx_heat_from_ac , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         eflx_sh_snow            => veg_ef%eflx_sh_snow      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_soil            => veg_ef%eflx_sh_soil      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_h2osfc          => veg_ef%eflx_sh_h2osfc    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from surface water (W/m**2) [+ to atm]
         eflx_sh_grnd            => veg_ef%eflx_sh_grnd      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_lwrad_net          => veg_ef%eflx_lwrad_net    , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
         eflx_wasteheat_patch    => veg_ef%eflx_wasteheat    , & ! Output: [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_heat_from_ac_patch => veg_ef%eflx_heat_from_ac , & ! Output: [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         eflx_traffic_patch      => veg_ef%eflx_traffic      , & ! Output: [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)
         eflx_anthro             => veg_ef%eflx_anthro       , & ! Output: [real(r8) (:)   ]  total anthropogenic heat flux (W/m**2)
         eflx_gnet               => veg_ef%eflx_gnet         , & ! Output: [real(r8) (:)   ]  net ground heat flux into the surface (W/m**2)
         dgnetdT                 => veg_ef%dgnetdT           , & ! Output: [real(r8) (:)   ]  temperature derivative of ground net heat flux

         sabg                    => solarabs_vars%sabg_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
         sabg_soil               => solarabs_vars%sabg_soil_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
         sabg_snow               => solarabs_vars%sabg_snow_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
         sabg_chk                => solarabs_vars%sabg_chk_patch            , & ! Output: [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
         sabg_lyr                => solarabs_vars%sabg_lyr_patch            , & ! Input:  [real(r8) (:,:) ]  absorbed solar radiation (pft,lyr) [W/m2]

         slope_deg               => grc_pp%slope_deg                          &
         )

      ! Each column sums over its patches sequentially, in the same order as a
      ! serial build, so the column fluxes are reproducible across thread counts.

      !$acc parallel loop gang vector default(present) &
      !$acc    private(c, p, l, t, g, j, lyr_top, lwrad_emit, dlwrad_emit, lwrad_emit_snow, &
      !$acc            lwrad_emit_soil, lwrad_emit_h2osfc, slope_rad, eflx_gnet_top, &
      !$acc            eflx_gnet_snow, eflx_gnet_soil, eflx_gnet_h2osfc)
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         g = col_pp%gridcell(c)
         l = col_pp%landunit(c)

         ! Net ground heat flux into the surface and its temperature derivative
         ! Added a pfts loop here to get the average of hs and dhsdT over
         ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

         lwrad_emit  =    emg(c) * sb * t_grnd(c)**4
         dlwrad_emit = 4._r8*emg(c) * sb * t_grnd(c)**3

         ! fractionate lwrad_emit; balanced in CanopyFluxes & Biogeophysics2
         lwrad_emit_snow    =    emg(c) * sb * t_soisno(c,snl(c)+1)**4
         lwrad_emit_soil    =    emg(c) * sb * t_soisno(c,1)**4
         lwrad_emit_h2osfc  =    emg(c) * sb * t_h2osfc(c)**4

         if (use_finetop_rad .and. (.not. lun_pp%urbpoi(l))) then
            slope_rad = slope_deg(g) * deg2rad
            lwrad_emit = lwrad_emit / cos(slope_rad)
            dlwrad_emit = dlwrad_emit / cos(slope_rad)
            lwrad_emit_snow = lwrad_emit_snow / cos(slope_rad)
            lwrad_emit_soil = lwrad_emit_soil / cos(slope_rad)
            lwrad_emit_h2osfc = lwrad_emit_h2osfc / cos(slope_rad)
         endif

         hs_soil(fc)   = 0._r8
         hs_h2osfc(fc) = 0._r8
         dhsdT(fc)     = 0._r8

         !$acc loop seq
         do p = col_pp%pfti(c), col_pp%pfti(c) + col_pp%npfts(c) - 1
            l = veg_pp%landunit(p)
            t = veg_pp%topounit(p)

            if (veg_pp%active(p)) then
               if (.not. lun_pp%urbpoi(l)) then
                  eflx_gnet(p) = sabg(p) + dlrad(p) &
                       + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(t) - lwrad_emit &
                       - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
                  ! save sabg for balancecheck, in case frac_sno is set to zero later
                  sabg_chk(p) = frac_sno_eff(c) * sabg_snow(p) + (1._r8 - frac_sno_eff(c) ) * sabg_soil(p)

                  eflx_gnet_snow = sabg_snow(p) + dlrad(p) &
                       + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(t) - lwrad_emit_snow &
                       - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

                  eflx_gnet_soil = sabg_soil(p) + dlrad(p) &
                       + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(t) - lwrad_emit_soil &
                       - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

                  eflx_gnet_h2osfc = sabg_soil(p) + dlrad(p) &
                       + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(t) - lwrad_emit_h2osfc &
                       - (eflx_sh_h2osfc(p)+qflx_ev_h2osfc(p)*htvp(c))
               else
                  ! For urban columns we use the net longwave radiation (eflx_lwrad_net) because of
                  ! interactions between urban columns.

                  ! All wasteheat and traffic flux goes into canyon floor
                  if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                     eflx_wasteheat_patch(p) = eflx_wasteheat(l)/(1._r8-lun_pp%wtlunit_roof(l))
                     eflx_heat_from_ac_patch(p) = eflx_heat_from_ac(l)/(1._r8-lun_pp%wtlunit_roof(l))
                     eflx_traffic_patch(p) = eflx_traffic(l)/(1._r8-lun_pp%wtlunit_roof(l))
                  else
                     eflx_wasteheat_patch(p) = 0._r8
                     eflx_heat_from_ac_patch(p) = 0._r8
                     eflx_traffic_patch(p) = 0._r8
                  end if
                  ! Include transpiration term because needed for previous road
                  ! and include wasteheat and traffic flux
                  eflx_gnet(p) = sabg(p) + dlrad(p)  &
                       - eflx_lwrad_net(p) &
                       - (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + qflx_tran_veg(p)*hvap) &
                       + eflx_wasteheat_patch(p) + eflx_heat_from_ac_patch(p) + eflx_traffic_patch(p)
                  eflx_anthro(p)   = eflx_wasteheat_patch(p) + eflx_traffic_patch(p)
                  eflx_gnet_snow   = eflx_gnet(p)
                  eflx_gnet_soil   = eflx_gnet(p)
                  eflx_gnet_h2osfc = eflx_gnet(p)
               end if
               dgnetdT(p) = - cgrnd(p) - dlwrad_emit
               dhsdT(fc) = dhsdT(fc) + dgnetdT(p) * veg_pp%wtcol(p)
               ! separate surface fluxes for soil/snow
               hs_soil(fc) = hs_soil(fc) + eflx_gnet_soil * veg_pp%wtcol(p)
               hs_h2osfc(fc) = hs_h2osfc(fc) + eflx_gnet_h2osfc * veg_pp%wtcol(p)
            end if
         end do

         ! Additional calculations with SNICAR:
         ! Set up tridiagonal matrix in a new manner. There is now
         ! absorbed solar radiation in each snow layer, instead of
         ! only the surface. Following the current implementation,
         ! absorbed solar flux should be: S + ((delS/delT)*dT),
         ! where S is absorbed radiation, and T is temperature. Now,
         ! assume delS/delT is zero, then it is OK to just add S
         ! to each layer

         ! Initialize:
         !$acc loop seq
         do j = -nlevsno+1,1
            sabg_lyr_col(fc,j) = 0._r8
         end do
         hs_top(fc)      = 0._r8
         hs_top_snow(fc) = 0._r8

         lyr_top = snl(c) + 1
         !$acc loop seq
         do p = col_pp%pfti(c), col_pp%pfti(c) + col_pp%npfts(c) - 1
            if (veg_pp%active(p)) then
               t = veg_pp%topounit(p)
               l = veg_pp%landunit(p)
               if (.not. lun_pp%urbpoi(l)) then

                  eflx_gnet_top = sabg_lyr(p,lyr_top) + dlrad(p) + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(t) &
                       - lwrad_emit - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

                  hs_top(fc) = hs_top(fc) + eflx_gnet_top*veg_pp%wtcol(p)

                  eflx_gnet_snow = sabg_lyr(p,lyr_top) + dlrad(p) + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(t) &
                       - lwrad_emit_snow - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

                  hs_top_snow(fc) = hs_top_snow(fc) + eflx_gnet_snow*veg_pp%wtcol(p)

                  !$acc loop seq
                  do j = lyr_top,1,1
                     sabg_lyr_col(fc,j) = sabg_lyr_col(fc,j) + sabg_lyr(p,j) * veg_pp%wtcol(p)
                  enddo
               else

                  hs_top(fc)      = hs_top(fc) + eflx_gnet(p)*veg_pp%wtcol(p)
                  hs_top_snow(fc) = hs_top_snow(fc) + eflx_gnet(p)*veg_pp%wtcol(p)
                  sabg_lyr_col(fc,lyr_top) = sabg_lyr_col(fc,lyr_top) + sabg(p) * veg_pp%wtcol(p)

               endif
            endif
         enddo
      enddo

    end associate

  end subroutine ComputeGroundHeatFluxAndDeriv

  !-----------------------------------------------------------------------
  subroutine ComputeHeatDiffFluxAndFactor(num_nolakec, filter_nolakec, dtime, &
       tk, cv, fn)
    !
    ! !DESCRIPTION:
    ! Computes:
    ! (1) Heat diffusion at the interface of layers.
    ! (2) Factor used in computing tridiagonal matrix (col_es%fact)
    !
    ! !USES:
    use elm_varcon     , only : capr
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)  :: num_nolakec                         ! number of column non-lake points in column filter
    integer  , intent(in)  :: filter_nolakec(:)                   ! column filter for non-lake points
    real(r8) , intent(in)  :: dtime                               ! land model time step (sec)
    real(r8) , intent(in)  :: tk(num_nolakec,-nlevsno+1:nlevgrnd) ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)  :: cv(num_nolakec,-nlevsno+1:nlevgrnd) ! heat capacity [J/(m2 K)]
    real(r8) , intent(out) :: fn(num_nolakec,-nlevsno+1:nlevgrnd) ! heat diffusion through the layer interface [W/m2]
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                             ! indices
    integer  :: fc                                                ! lake filtered column indices
    !-----------------------------------------------------------------------

    associate(&
         zi         => col_pp%zi         , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)
         dz         => col_pp%dz         , & ! Input:  [real(r8) (:,:) ] layer depth (m)
         z          => col_pp%z          , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         t_building => lun_es%t_building , & ! Input:  [real(r8) (:)   ] internal building temperature (K)
         t_soisno   => col_es%t_soisno   , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         fact       => col_es%fact       , & ! Output: [real(r8) (:,:) ] used in computing tridiagonal matrix
         eflx_bot   => col_ef%eflx_bot     & ! Input:  [real(r8) (:)   ] heat flux from beneath column (W/m**2) [+ = upward]
         )

      ! Determine heat diffusion through the layer interface and factor used in computing
      ! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
      ! matrix and solve system

      !$acc parallel loop gang vector collapse(2) default(present) private(c,l)
      do j = -nlevsno+1,nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col_pp%landunit(c)
            if ((col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) .and. j <= nlevurb) then
               if (j >= col_pp%snl(c)+1) then
                  if (j == col_pp%snl(c)+1) then
                     fact(c,j) = dtime/cv(fc,j)
                     fn(fc,j) = tk(fc,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j <= nlevurb-1) then
                     fact(c,j) = dtime/cv(fc,j)
                     fn(fc,j) = tk(fc,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j == nlevurb) then
                     fact(c,j) = dtime/cv(fc,j)
                     ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                     ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                     ! building temperature. (See Oleson urban notes of 6/18/03).
                     fn(fc,j) = tk(fc,j) * (t_building(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
                  end if
               end if
            else if (col_pp%itype(c) /= icol_sunwall .and. col_pp%itype(c) /= icol_shadewall &
                 .and. col_pp%itype(c) /= icol_roof) then
               if (j >= col_pp%snl(c)+1) then
                  if (j == col_pp%snl(c)+1) then
                     fact(c,j) = dtime/cv(fc,j) * dz(c,j) / (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
                     fn(fc,j) = tk(fc,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j <= nlevgrnd-1) then
                     fact(c,j) = dtime/cv(fc,j)
                     fn(fc,j) = tk(fc,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j == nlevgrnd) then
                     fact(c,j) = dtime/cv(fc,j)
                     fn(fc,j) = eflx_bot(c)
                  end if
               end if
            end if
         end do
      end do

    end associate

  end subroutine ComputeHeatDiffFluxAndFactor

  !-----------------------------------------------------------------------
  subroutine SetRHSVec(num_filter, filter, dtime, &
       hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, &
       tk_h2osfc, fn, dz_h2osfc, rvector)
    !
    ! !DESCRIPTION:
    ! Setup the RHS-Vector for the numerical solution of temperature for snow,
    ! standing surface water and soil layers.
    !
    !           |===========|
    !           |   Snow    |
    !           !===========|
    ! rvector = |    SSW    |
    !           !===========|
    !           !   Soil    |
    !           !===========|
    !
    ! rvector must be initialized by the caller; the entries for the inactive
    ! snow layers are not set.
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                               ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                ! column filter
    real(r8) , intent(in)    :: dtime                                    ! land model time step (sec)
    real(r8) , intent(in)    :: hs_h2osfc(num_filter)                    ! heat flux on standing water [W/m2]
    real(r8) , intent(in)    :: hs_top_snow(num_filter)                  ! heat flux on top snow layer [W/m2]
    real(r8) , intent(in)    :: hs_soil(num_filter)                      ! heat flux on soil [W/m2]
    real(r8) , intent(in)    :: hs_top(num_filter)                       ! net energy flux into surface layer (col) [W/m2]
    real(r8) , intent(in)    :: dhsdT(num_filter)                        ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: sabg_lyr_col(num_filter,-nlevsno+1:1)    ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) , intent(in)    :: tk_h2osfc(num_filter)                    ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8) , intent(in)    :: fn(num_filter,-nlevsno+1:nlevgrnd)       ! heat diffusion through the layer interface [W/m2]
    real(r8) , intent(in)    :: dz_h2osfc(num_filter)                    ! Thickness of standing water [m]
    real(r8) , intent(inout) :: rvector(num_filter,-nlevsno:nlevgrnd)    ! RHS vector used in numerical solution of temperature
    !-----------------------------------------------------------------------

    ! Set entries in RHS vector for snow layers
    call SetRHSVec_Snow(num_filter, filter, &
         hs_top_snow, hs_top, dhsdT, sabg_lyr_col, fn, rvector)

    ! Set entries in RHS vector for soil layers
    call SetRHSVec_Soil(num_filter, filter, &
         hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fn, rvector)

    ! Set entries in RHS vector for surface water layer, and add the standing
    ! water contribution to the top soil layer. Must follow SetRHSVec_Soil.
    call SetRHSVec_StandingSurfaceWater(num_filter, filter, dtime, &
         hs_h2osfc, hs_soil, dhsdT, tk_h2osfc, dz_h2osfc, rvector)

  end subroutine SetRHSVec

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_Snow(num_filter, filter, &
       hs_top_snow, hs_top, dhsdT, sabg_lyr_col, fn, rt)
    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to snow layers for all columns.
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                               ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                ! column filter
    real(r8) , intent(in)    :: hs_top_snow(num_filter)                  ! heat flux on top snow layer [W/m2]
    real(r8) , intent(in)    :: hs_top(num_filter)                       ! net energy flux into surface layer (col) [W/m2]
    real(r8) , intent(in)    :: dhsdT(num_filter)                        ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: sabg_lyr_col(num_filter,-nlevsno+1:1)    ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) , intent(in)    :: fn(num_filter,-nlevsno+1:nlevgrnd)       ! heat diffusion through the layer interface [W/m2]
    real(r8) , intent(inout) :: rt(num_filter,-nlevsno:nlevgrnd)         ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                      ! indices
    integer  :: fc                                                       ! lake filtered column indices
    real(r8) :: hs_top_lev                                               ! heat flux on top layer [W/m2]
    !-----------------------------------------------------------------------

    associate(                            &
         t_soisno => col_es%t_soisno    , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         fact     => col_es%fact          & ! Input:  [real(r8) (:,:) ] used in computing tridiagonal matrix
         )

      !$acc parallel loop gang vector collapse(2) default(present) private(c, hs_top_lev)
      do j = -nlevsno+1,0
         do fc = 1,num_filter
            c = filter(fc)
            if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) then
               ! urban sunwall, shadewall, and roof columns
               hs_top_lev = hs_top(fc)
            else
               ! urban road and non-urban columns
               hs_top_lev = hs_top_snow(fc)
            end if

            if (j == col_pp%snl(c)+1) then
               rt(fc,j-1) = t_soisno(c,j) +  fact(c,j)*( hs_top_lev &
                    - dhsdT(fc)*t_soisno(c,j) + cnfac*fn(fc,j) )
            else if (j > col_pp%snl(c)+1) then
               rt(fc,j-1) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(fc,j) - fn(fc,j-1) )
               rt(fc,j-1) = rt(fc,j-1) + fact(c,j)*sabg_lyr_col(fc,j)
            end if
         enddo
      end do

    end associate

  end subroutine SetRHSVec_Snow

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_StandingSurfaceWater(num_filter, filter, dtime, &
       hs_h2osfc, hs_soil, dhsdT, tk_h2osfc, dz_h2osfc, rt)
    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to standing surface water, and adds the
    ! standing surface water contribution to the top soil layer.
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                         ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                          ! column filter
    real(r8) , intent(in)    :: dtime                              ! land model time step (sec)
    real(r8) , intent(in)    :: hs_h2osfc(num_filter)              ! heat flux on standing water [W/m2]
    real(r8) , intent(in)    :: hs_soil(num_filter)                ! heat flux on soil [W/m2]
    real(r8) , intent(in)    :: dhsdT(num_filter)                  ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: tk_h2osfc(num_filter)              ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8) , intent(in)    :: dz_h2osfc(num_filter)              ! Thickness of standing water [m]
    real(r8) , intent(inout) :: rt(num_filter,-nlevsno:nlevgrnd)   ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: c                                                  ! indices
    integer  :: fc                                                 ! lake filtered column indices
    real(r8) :: dzm                                                ! used in computing tridiagonal matrix
    real(r8) :: fn_h2osfc                                          ! heat diffusion through standing-water/soil interface [W/m2]
    !-----------------------------------------------------------------------

    associate(                              &
         t_soisno    => col_es%t_soisno   , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         t_h2osfc    => col_es%t_h2osfc   , & ! Input:  [real(r8) (:)   ] surface water temperature
         c_h2osfc    => col_es%c_h2osfc   , & ! Input:  [real(r8) (:)   ] heat capacity of surface water
         fact        => col_es%fact       , & ! Input:  [real(r8) (:,:) ] used in computing tridiagonal matrix
         frac_h2osfc => col_ws%frac_h2osfc  & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         )

      !$acc parallel loop gang vector default(present) private(c, dzm, fn_h2osfc)
      do fc = 1,num_filter
         c = filter(fc)

         ! surface water layer has two coefficients
         dzm=(0.5*dz_h2osfc(fc)+col_pp%z(c,1))

         fn_h2osfc=tk_h2osfc(fc)*(t_soisno(c,1)-t_h2osfc(c))/dzm
         rt(fc,0)= t_h2osfc(c) +  (dtime/c_h2osfc(c)) &
              *( hs_h2osfc(fc) - dhsdT(fc)*t_h2osfc(c) + cnfac*fn_h2osfc )!rhs for h2osfc

         if ( frac_h2osfc(c) /= 0.0_r8 )then
            rt(fc,1)=rt(fc,1) &
                 -frac_h2osfc(c)*fact(c,1)*((hs_soil(fc) - dhsdT(fc)*t_soisno(c,1)) &
                 +cnfac*fn_h2osfc)
         end if
      enddo

    end associate

  end subroutine SetRHSVec_StandingSurfaceWater

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_Soil(num_filter, filter, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fn, rt)
    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to soil layers for all columns.
    ! The standing surface water contribution is added by
    ! SetRHSVec_StandingSurfaceWater.
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                               ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                ! column filter
    real(r8) , intent(in)    :: hs_top_snow(num_filter)                  ! heat flux on top snow layer [W/m2]
    real(r8) , intent(in)    :: hs_soil(num_filter)                      ! heat flux on soil [W/m2]
    real(r8) , intent(in)    :: hs_top(num_filter)                       ! net energy flux into surface layer (col) [W/m2]
    real(r8) , intent(in)    :: dhsdT(num_filter)                        ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: sabg_lyr_col(num_filter,-nlevsno+1:1)    ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) , intent(in)    :: fn(num_filter,-nlevsno+1:nlevgrnd)       ! heat diffusion through the layer interface [W/m2]
    real(r8) , intent(inout) :: rt(num_filter,-nlevsno:nlevgrnd)         ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                    ! indices
    integer  :: fc                                                       ! lake filtered column indices
    !-----------------------------------------------------------------------

    associate(                                &
         t_soisno     => col_es%t_soisno    , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         fact         => col_es%fact        , & ! Input:  [real(r8) (:,:) ] used in computing tridiagonal matrix
         frac_sno_eff => col_ws%frac_sno_eff  & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         )

      !
      ! urban sunwall, shadewall, and roof columns -------------------------------------
      !
      !$acc parallel loop gang vector collapse(2) default(present) private(c)
      do j = 1,nlevurb
         do fc = 1,num_filter
            c = filter(fc)
            if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) then
               if (j == col_pp%snl(c)+1) then
                  rt(fc,j) = t_soisno(c,j) +  fact(c,j)*( hs_top(fc) - dhsdT(fc)*t_soisno(c,j) + cnfac*fn(fc,j) )
               else if (j <= nlevurb-1) then
                  ! if this is a snow layer or the top soil layer,
                  ! add absorbed solar flux to factor 'rt'
                  if (j == 1) then
                     rt(fc,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(fc,j) - fn(fc,j-1) )
                     rt(fc,j) = rt(fc,j) + (fact(c,j)*sabg_lyr_col(fc,j))
                  else
                     rt(fc,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(fc,j) - fn(fc,j-1) )
                  endif
               else if (j == nlevurb) then
                  ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                  ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                  ! building temperature. (See Oleson urban notes of 6/18/03).
                  rt(fc,j) = t_soisno(c,j) + fact(c,j)*( fn(fc,j) - cnfac*fn(fc,j-1) )
               end if
            end if
         enddo
      end do

      !
      ! non-urban columns and urban road columns ---------------------------------------
      !
      !$acc parallel loop gang vector collapse(2) default(present) private(c,l)
      do j = 1,nlevgrnd
         do fc = 1,num_filter
            c = filter(fc)
            l = col_pp%landunit(c)
            if ((.not. lun_pp%urbpoi(l)) .or. &
                 col_pp%itype(c) == icol_road_imperv .or. col_pp%itype(c) == icol_road_perv) then
               if (j == col_pp%snl(c)+1) then
                  rt(fc,j) = t_soisno(c,j) +  fact(c,j)*( hs_top_snow(fc) &
                       - dhsdT(fc)*t_soisno(c,j) + cnfac*fn(fc,j) )
               else if (j == 1) then
                  ! this is the snow/soil interface layer
                  rt(fc,j) = t_soisno(c,j) + fact(c,j) &
                       *((1._r8-frac_sno_eff(c))*(hs_soil(fc) - dhsdT(fc)*t_soisno(c,j)) &
                       + cnfac*(fn(fc,j) - frac_sno_eff(c) * fn(fc,j-1)))

                  rt(fc,j) = rt(fc,j) +  frac_sno_eff(c)*fact(c,j)*sabg_lyr_col(fc,j)

               else if (j <= nlevgrnd-1) then
                  rt(fc,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(fc,j) - fn(fc,j-1) )

               else if (j == nlevgrnd) then
                  rt(fc,j) = t_soisno(c,j) - cnfac*fact(c,j)*fn(fc,j-1) + fact(c,j)*fn(fc,j)
               end if
            end if
         enddo
      end do

    end associate

  end subroutine SetRHSVec_Soil

  !-----------------------------------------------------------------------
  subroutine SetMatrix(num_filter, filter, dtime, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, bmatrix)
    !
    ! !DESCRIPTION:
    ! Setup the matrix for the numerical solution of temperature for snow,
    ! standing surface water and soil layers. bmatrix must be zeroed by the
    ! caller; each block below only sets its own nonzero entries.
    !
    !           |===========|===========|===========|
    !           |   Snow    |           | Snow-Soil |
    !           !===========|===========|===========|
    ! bmatrix = |           |    SSW    |  SSW-Soil |
    !           !===========|===========|===========|
    !           ! Soil-Snow | Soil-SSW  |    Soil   |
    !           !===========|===========|===========|
    !
    ! Non-zero pattern of bmatrix:
    !
    !        SNOW-LAYERS
    !            |
    !            |  STANDING-SURFACE-WATER
    !            |         |
    !            |         |              SOIL-LAYERS
    !            |         |                  |
    !            v         v                  v
    !
    !      -5 -4 -3 -2 -1| 0| 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    !      ==============================================================
    !  -5 | x  x         |  |                                            |
    !  -4 | x  x  x      |  |                                            |
    !  -3 |    x  x  x   |  |                                            |
    !  -2 |       x  x  x|  |                                            |
    !  -1 |          x  x|  | x                                          |
    !      ==============================================================
    !   0 |              | x| x                                          |
    !      ==============================================================
    !   1 |             x| x| x  x                                       |
    !   2 |              |  | x  x  x                                    |
    !   3 |              |  |    x  x  x                                 |
    !   4 |              |  |       x  x  x                              |
    !   5 |              |  |          x  x  x                           |
    !   6 |              |  |             x  x  x                        |
    !   7 |              |  |                x  x  x                     |
    !   8 |              |  |                   x  x  x                  |
    !   9 |              |  |                      x  x  x               |
    !  10 |              |  |                         x  x  x            |
    !  11 |              |  |                            x  x  x         |
    !  12 |              |  |                               x  x  x      |
    !  13 |              |  |                                  x  x  x   |
    !  14 |              |  |                                     x  x  x|
    !  15 |              |  |                                        x  x|
    !      ==============================================================
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                                  ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                   ! column filter
    real(r8) , intent(in)    :: dtime                                       ! land model time step (sec)
    real(r8) , intent(in)    :: dhsdT(num_filter)                           ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: tk(num_filter,-nlevsno+1:nlevgrnd)          ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)    :: tk_h2osfc(num_filter)                       ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)    :: dz_h2osfc(num_filter)                       ! Thickness of standing water [m]
    real(r8) , intent(inout) :: bmatrix(num_filter,nband,-nlevsno:nlevgrnd) ! matrix for numerical solution of temperature
    !-----------------------------------------------------------------------

    call SetMatrix_Snow(num_filter, filter, &
         dhsdT, tk, bmatrix)

    call SetMatrix_Soil(num_filter, filter, &
         dhsdT, tk, tk_h2osfc, dz_h2osfc, bmatrix)

    call SetMatrix_StandingSurfaceWater(num_filter, filter, dtime, &
         dhsdT, tk_h2osfc, dz_h2osfc, bmatrix)

  end subroutine SetMatrix

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow(num_filter, filter, &
       dhsdT, tk, bmatrix)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries corresponding to internal snow layers and
    ! snow-soil interaction for all columns
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                                  ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                   ! column filter
    real(r8) , intent(in)    :: dhsdT(num_filter)                           ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: tk(num_filter,-nlevsno+1:nlevgrnd)          ! thermal conductivity [W/(m K)]
    real(r8) , intent(inout) :: bmatrix(num_filter,nband,-nlevsno:nlevgrnd) ! matrix entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                         ! indices
    integer  :: fc                                                          ! lake filtered column indices
    real(r8) :: dzm                                                         ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                         ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    associate(                     &
         z    => col_pp%z        , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
         fact => col_es%fact       & ! Input:  [real(r8) (:,:)]  used in computing tridiagonal matrix
         )

      !$acc parallel loop gang vector collapse(2) default(present) private(c, dzm, dzp)
      do j = -nlevsno+1,0
         do fc = 1,num_filter
            c = filter(fc)
            if (j >= col_pp%snl(c)+1) then
               dzp     = z(c,j+1)-z(c,j)
               if (j == col_pp%snl(c)+1) then
                  bmatrix(fc,4,j-1) = 0._r8
                  bmatrix(fc,3,j-1) = 1._r8+(1._r8-cnfac)*fact(c,j)*tk(fc,j)/dzp-fact(c,j)*dhsdT(fc)
               else
                  dzm     = (z(c,j)-z(c,j-1))
                  bmatrix(fc,4,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j-1)/dzm
                  bmatrix(fc,3,j-1) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(fc,j)/dzp + tk(fc,j-1)/dzm)
               end if
               ! band 2 of the bottom snow layer (j=0) is the snow-soil coupling;
               ! band 1 at row -1 couples to the top soil layer
               if (j /= 0) then
                  bmatrix(fc,2,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j)/dzp
               else
                  bmatrix(fc,1,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j)/dzp
               end if
            end if
         enddo
      end do

    end associate

  end subroutine SetMatrix_Snow

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil(num_filter, filter, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, bmatrix)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries corresponding to internal soil layers and
    ! soil-snow interaction for all columns, including the diagonal
    ! correction for the presence of standing surface water
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                                  ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                   ! column filter
    real(r8) , intent(in)    :: dhsdT(num_filter)                           ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: tk(num_filter,-nlevsno+1:nlevgrnd)          ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)    :: tk_h2osfc(num_filter)                       ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)    :: dz_h2osfc(num_filter)                       ! Thickness of standing water [m]
    real(r8) , intent(inout) :: bmatrix(num_filter,nband,-nlevsno:nlevgrnd) ! matrix entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                       ! indices
    integer  :: fc                                                          ! lake filtered column indices
    real(r8) :: dzm                                                         ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                         ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    associate(                                 &
         zi           => col_pp%zi           , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
         z            => col_pp%z            , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
         fact         => col_es%fact         , & ! Input:  [real(r8) (:,:)]  used in computing tridiagonal matrix
         frac_h2osfc  => col_ws%frac_h2osfc  , & ! Input:  [real(r8) (:)  ]  fraction of ground covered by surface water (0 to 1)
         frac_sno_eff => col_ws%frac_sno_eff   & ! Input:  [real(r8) (:)  ]  fraction of ground covered by snow (0 to 1)
         )

      ! Band 5 of the top soil layer (j=1) is the soil-snow coupling.

      !
      ! urban sunwall, shadewall, and roof columns -------------------------------------
      !
      !$acc parallel loop gang vector collapse(2) default(present) private(c, dzm, dzp)
      do j = 1,nlevurb
         do fc = 1,num_filter
            c = filter(fc)
            if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) then
               if (j == col_pp%snl(c)+1) then
                  dzp     = z(c,j+1)-z(c,j)
                  if (j /= 1) then
                     bmatrix(fc,4,j) = 0._r8
                  else
                     bmatrix(fc,5,j) = 0._r8
                  end if
                  bmatrix(fc,3,j) = 1._r8+(1._r8-cnfac)*fact(c,j)*tk(fc,j)/dzp-fact(c,j)*dhsdT(fc)
                  bmatrix(fc,2,j) =  -(1._r8-cnfac)*fact(c,j)*tk(fc,j)/dzp
               else if (j <= nlevurb-1) then
                  dzm     = (z(c,j)-z(c,j-1))
                  dzp     = (z(c,j+1)-z(c,j))
                  if (j /= 1) then
                     bmatrix(fc,4,j) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j-1)/dzm
                  else
                     bmatrix(fc,5,j) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j-1)/dzm
                  end if
                  bmatrix(fc,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(fc,j)/dzp + tk(fc,j-1)/dzm)
                  bmatrix(fc,2,j) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j)/dzp
               else if (j == nlevurb) then
                  ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                  ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                  ! building temperature. (See Oleson urban notes of 6/18/03).
                  dzm     = ( z(c,j)-z(c,j-1))
                  dzp     = (zi(c,j)-z(c,j))
                  bmatrix(fc,4,j) =   - (1._r8-cnfac)*fact(c,j)*(tk(fc,j-1)/dzm)
                  bmatrix(fc,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(fc,j-1)/dzm + tk(fc,j)/dzp)
                  bmatrix(fc,2,j) = 0._r8
               end if
            end if
         enddo
      end do

      !
      ! non-urban columns and urban road columns ---------------------------------------
      !
      !$acc parallel loop gang vector collapse(2) default(present) private(c, l, dzm, dzp)
      do j = 1,nlevgrnd
         do fc = 1,num_filter
            c = filter(fc)
            l = col_pp%landunit(c)
            if ((.not. lun_pp%urbpoi(l)) .or. &
                 col_pp%itype(c) == icol_road_imperv .or. col_pp%itype(c) == icol_road_perv) then
               if (j == col_pp%snl(c)+1) then
                  dzp     = z(c,j+1)-z(c,j)
                  if (j /= 1) then
                     bmatrix(fc,4,j) = 0._r8
                  else
                     bmatrix(fc,5,j) = 0._r8
                  end if
                  bmatrix(fc,3,j) = 1._r8+(1._r8-cnfac)*fact(c,j)*tk(fc,j)/dzp-fact(c,j)*dhsdT(fc)
                  bmatrix(fc,2,j) =  -(1._r8-cnfac)*fact(c,j)*tk(fc,j)/dzp
               else if (j == 1) then
                  ! this is the snow/soil interface layer
                  dzm     = (z(c,j)-z(c,j-1))
                  dzp     = (z(c,j+1)-z(c,j))
                  bmatrix(fc,3,j) = 1._r8 + (1._r8-cnfac)*fact(c,j)*(tk(fc,j)/dzp &
                       + frac_sno_eff(c) * tk(fc,j-1)/dzm) &
                       - (1._r8 - frac_sno_eff(c))*fact(c,j)*dhsdT(fc)
                  bmatrix(fc,2,j) = - (1._r8-cnfac)*fact(c,j)*tk(fc,j)/dzp
                  bmatrix(fc,5,j) =   - frac_sno_eff(c) * (1._r8-cnfac) * fact(c,j) &
                       * tk(fc,j-1)/dzm
               else if (j <= nlevgrnd-1) then
                  dzm     = (z(c,j)-z(c,j-1))
                  dzp     = (z(c,j+1)-z(c,j))
                  bmatrix(fc,4,j) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j-1)/dzm
                  bmatrix(fc,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(fc,j)/dzp + tk(fc,j-1)/dzm)
                  bmatrix(fc,2,j) =   - (1._r8-cnfac)*fact(c,j)* tk(fc,j)/dzp
               else if (j == nlevgrnd) then
                  dzm     = (z(c,j)-z(c,j-1))
                  bmatrix(fc,4,j) =   - (1._r8-cnfac)*fact(c,j)*tk(fc,j-1)/dzm
                  bmatrix(fc,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*tk(fc,j-1)/dzm
                  bmatrix(fc,2,j) = 0._r8
               end if
            end if
         enddo
      end do

      !
      ! diagonal element correction for presence of h2osfc -----------------------------
      !
      !$acc parallel loop gang vector default(present) private(c, dzm)
      do fc = 1,num_filter
         c = filter(fc)

         ! surface water layer has two coefficients
         dzm=(0.5*dz_h2osfc(fc)+z(c,1))

         if ( frac_h2osfc(c) /= 0.0_r8 ) then
            bmatrix(fc,3,1)=bmatrix(fc,3,1)+ frac_h2osfc(c) &
                 *((1._r8-cnfac)*fact(c,1)*tk_h2osfc(fc)/dzm + fact(c,1)*dhsdT(fc))
         end if

      enddo

    end associate

  end subroutine SetMatrix_Soil

  !-----------------------------------------------------------------------
  subroutine SetMatrix_StandingSurfaceWater(num_filter, filter, dtime, &
       dhsdT, tk_h2osfc, dz_h2osfc, bmatrix)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries corresponding to standing surface water,
    ! standing surface water-soil interaction and soil-standing surface
    ! water interaction
    !
    ! !ARGUMENTS:
    implicit none
    integer  , intent(in)    :: num_filter                                  ! number of columns in the filter
    integer  , intent(in)    :: filter(:)                                   ! column filter
    real(r8) , intent(in)    :: dtime                                       ! land model time step (sec)
    real(r8) , intent(in)    :: dhsdT(num_filter)                           ! temperature derivative of "hs" [col]
    real(r8) , intent(in)    :: tk_h2osfc(num_filter)                       ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)    :: dz_h2osfc(num_filter)                       ! Thickness of standing water [m]
    real(r8) , intent(inout) :: bmatrix(num_filter,nband,-nlevsno:nlevgrnd) ! matrix entries
    !
    ! !LOCAL VARIABLES:
    integer  :: c                                                           ! indices
    integer  :: fc                                                          ! lake filtered column indices
    real(r8) :: dzm                                                         ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    associate(                               &
         c_h2osfc    => col_es%c_h2osfc    , & ! Input:  [real(r8) (:)  ]  heat capacity of surface water
         fact        => col_es%fact        , & ! Input:  [real(r8) (:,:)]  used in computing tridiagonal matrix
         frac_h2osfc => col_ws%frac_h2osfc   & ! Input:  [real(r8) (:)  ]  fraction of ground covered by surface water (0 to 1)
         )

      !$acc parallel loop gang vector default(present) private(c, dzm)
      do fc = 1,num_filter
         c = filter(fc)

         ! surface water layer has two coefficients
         dzm=(0.5*dz_h2osfc(fc)+col_pp%z(c,1))

         bmatrix(fc,3,0)= 1._r8+(1._r8-cnfac)*(dtime/c_h2osfc(c)) &
              *tk_h2osfc(fc)/dzm -(dtime/c_h2osfc(c))*dhsdT(fc) !interaction from atm

         bmatrix(fc,2,0)= -(1._r8-cnfac)*(dtime/c_h2osfc(c))*tk_h2osfc(fc)/dzm !flux to top soil layer

         ! top soil layer has sub coef shifted to 2nd super diagonal
         if ( frac_h2osfc(c) /= 0.0_r8 )then
            bmatrix(fc,4,1)=  - frac_h2osfc(c) * (1._r8-cnfac) * fact(c,1) &
                 * tk_h2osfc(fc)/dzm !flux from h2osfc
         end if
      enddo

    end associate

  end subroutine SetMatrix_StandingSurfaceWater

end module SoilTemperatureMod
