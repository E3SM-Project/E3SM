module SoilLateralFlowMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Lateral subsurface flow between grid cells, following Qiu et al. (2024),
  ! "Development of inter-grid-cell lateral unsaturated and saturated flow
  ! model in the E3SM Land Model (v2.0)", GMD 17, 143-167.
  !
  ! Lateral fluxes are computed explicitly from the state at the start of the
  ! soil-water solve:
  !
  ! - Unsaturated flow (Eq. 15): for each soil layer that lies entirely above
  !   the water table in both columns, a Darcy flux driven by the difference
  !   in matric potential plus elevation between the two columns. It enters
  !   the vertical (Zeng and Decker, 2009) Richards equation as a layer
  !   source term.
  ! - Saturated flow (Eq. 17): a Darcy flux driven by the difference in water
  !   table head (surface elevation minus water table depth), using a
  !   transmissivity based on the mean saturated thickness of the two columns.
  !   Each column's net saturated inflow is applied after the Richards solve
  !   at the water table (ApplySaturatedLateralFlux), as ELM does for
  !   baseflow: outflow drains the layers at and below the water table in
  !   turn, limited by the specific yield, and inflow fills the layers at and
  !   above the water table in turn. The layers below the water table stay
  !   saturated; the water table moves instead. (Zeng and Decker (2009) has
  !   no positive pressure below the water table, so a sink applied there
  !   inside the solve would desaturate the deep layers.)
  !
  ! Fluxes are first computed per connection and then summed for each column
  ! in the order of the neighbors' natural ids, so results do not depend on
  ! the number of MPI ranks.
  !
  ! Only naturally vegetated columns are connected (one per grid cell).
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use elm_varpar        , only : nlevgrnd
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_ws, col_wf
  use SoilStateType     , only : soilstate_type
  use SoilHydrologyType , only : soilhydrology_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ComputeLateralFlux         ! lateral fluxes into owned columns
  public :: ApplySaturatedLateralFlux  ! apply the saturated lateral flux at the water table
  public :: ThetaBasedWaterTable    ! diagnose the water table from soil moisture
  !
  ! !PRIVATE DATA:
  real(r8), parameter :: sat_lev = 0.96_r8   ! saturation level that defines the theta-based water table
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ComputeLateralFlux(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, soilstate_vars, qflx_lat_layer, qflx_sat, zwt_lat)
    !
    ! !DESCRIPTION:
    ! Computes, for the columns in filter_hydrologyc, the unsaturated lateral
    ! flux into each soil layer and the net saturated lateral flux into the
    ! column [mm H2O/s]. Sets col_wf%qflx_lateral (outflow positive, both
    ! parts), qflx_lateral_unsat, qflx_lateral_sat and qflx_lat_layer (the
    ! unsaturated part; ApplySaturatedLateralFlux adds the saturated part).
    !
    ! The fluxes use a water table diagnosed from the current soil moisture
    ! (zwt_lat), not zwt itself: after the previous time step's water-table
    ! diagnosis, Drainage moves zwt (by the specific yield of the drained
    ! water, and caps it at 80 m), so zwt can show a saturated zone that the
    ! soil moisture does not have (e.g. zwt = 80 m in an unsaturated 82 m deep
    ! column).
    !
    ! Must be called on all MPI ranks (it performs a halo exchange), with one
    ! clump per MPI rank: 'bounds' (clump bounds) must cover all owned columns.
    ! Ghost columns are addressed through the processor bounds.
    !
#ifdef MOAB_LATERAL
    use decompMod               , only : get_proc_bounds
    use elm_varctl              , only : lateral_unsat_flow, lateral_hk_anisotropy, lateral_theta_watertable
    use elm_varcon              , only : denh2o, denice
    use elm_varcon              , only : e_ice
    use ColumnConnectionSetType , only : c2c_connections
    use domainLateralMod        , only : NatVegColumnRealDataHaloExchange
#endif
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    real(r8)                 , intent(out)   :: qflx_lat_layer(bounds%begc:,1:)  ! unsaturated lateral flux into each soil layer [mm H2O/s]
    real(r8)                 , intent(out)   :: qflx_sat(bounds%begc:)           ! net saturated lateral flux into each column [mm H2O/s]
    real(r8)                 , intent(out)   :: zwt_lat(bounds%begc:)            ! water table depth used for the lateral fluxes [m]
    !
#ifdef MOAB_LATERAL
    ! !LOCAL VARIABLES:
    type(bounds_type)     :: bp              ! processor bounds (include ghost columns)
    integer               :: fc, c, j, k, iconn, cu, cd, nbu, nbd, nb, nconn, ncomp
    real(r8)              :: s1, bswl, hkl, impedl, smp_grad, face_area, area
    real(r8)              :: hsat_up, hsat_dn, trans, sgn
    real(r8), allocatable :: data(:,:)
    integer , allocatable :: jwt(:)
    real(r8), allocatable :: zw(:)           ! water table depth used for the fluxes, owned and ghost columns [m]
    real(r8)              :: vol(nlevgrnd)   ! volumetric water content [m3/m3]
    real(r8), allocatable :: qv_unsat(:,:)   ! [nconn, nlevgrnd] up-to-down volumetric flux per layer [mm H2O/s * m^2]
    real(r8), allocatable :: qv_sat(:)       ! [nconn] up-to-down volumetric saturated flux [mm H2O/s * m^2]
    real(r8)              :: unsat_in(nlevgrnd), sat_in
    !-----------------------------------------------------------------------

    associate(                                                  &
         zi           => col_pp%zi                            , & ! Input:  [real(r8) (:,:) ] interface depth (m)
         dz           => col_pp%dz                            , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         nlevbed      => col_pp%nlevbed                       , & ! Input:  [integer  (:)   ] number of hydrologically active layers
         h2osoi_vol   => col_ws%h2osoi_vol                    , & ! Input:  [real(r8) (:,:) ] volumetric soil water (m3/m3)
         smp_l        => soilstate_vars%smp_l_col             , & ! Input:  [real(r8) (:,:) ] soil matric potential (mm)
         watsat       => soilstate_vars%watsat_col            , & ! Input:  [real(r8) (:,:) ] porosity
         hksat        => soilstate_vars%hksat_col             , & ! Input:  [real(r8) (:,:) ] saturated hydraulic conductivity (mm/s)
         bsw          => soilstate_vars%bsw_col               , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
         origflag     => soilhydrology_vars%origflag          , & ! Input:  constant
         fracice      => soilhydrology_vars%fracice_col       , & ! Input:  [real(r8) (:,:) ] fractional impermeability
         icefrac      => soilhydrology_vars%icefrac_col       , & ! Input:  [real(r8) (:,:) ] fraction of ice
         zwt          => soilhydrology_vars%zwt_col           , & ! Input:  [real(r8) (:)   ] water table depth (m)
         conn         => c2c_connections                        &
         )

      call get_proc_bounds(bp)
      if (bounds%begc /= bp%begc .or. bounds%endc /= bp%endc) then
         call endrun(msg='ComputeLateralFlux: lateral flow requires one clump per MPI rank'//errMsg(__FILE__, __LINE__))
      end if

      ! --- Water table of owned columns, diagnosed from the current soil moisture
      allocate(zw(bp%begc_all:bp%endc_all))
      zw(:) = 0._r8
      do c = bp%begc, bp%endc
         zw(c) = zwt(c)
      end do
      do fc = 1, num_hydrologyc
         c  = filter_hydrologyc(fc)
         nb = nlevbed(c)
         if (lateral_theta_watertable) then
            do j = 1, nb
               vol(j) = col_ws%h2osoi_liq(c,j)/(dz(c,j)*denh2o) + col_ws%h2osoi_ice(c,j)/(dz(c,j)*denice)
            end do
            zw(c) = ThetaWaterTableDepth(nb, vol(1:nb), watsat(c,1:nb), col_pp%z(c,1:nb), zi(c,0:nb))
         else
            zw(c) = min(zwt(c), zi(c,nb))
         end if
      end do

      ! --- Copy the state of ghost columns from the owning ranks (one MPI round)
      ! Components: h2osoi_vol, smp_l, fracice, icefrac (each 1:nlevgrnd), water table
      ncomp = 4*nlevgrnd + 1
      allocate(data(bp%begc_all:bp%endc_all, ncomp))
      data(:,:) = 0._r8
      do c = bp%begc, bp%endc
         do j = 1, nlevgrnd
            data(c, 0*nlevgrnd + j) = h2osoi_vol(c,j)
            data(c, 1*nlevgrnd + j) = smp_l(c,j)
            data(c, 2*nlevgrnd + j) = fracice(c,j)
            data(c, 3*nlevgrnd + j) = icefrac(c,j)
         end do
         data(c, ncomp) = zw(c)
      end do

      call NatVegColumnRealDataHaloExchange(bp, 'lateral_flow_state', data)

      do c = bp%endc + 1, bp%endc_all
         do j = 1, nlevgrnd
            h2osoi_vol(c,j) = data(c, 0*nlevgrnd + j)
            smp_l(c,j)      = data(c, 1*nlevgrnd + j)
            fracice(c,j)    = data(c, 2*nlevgrnd + j)
            icefrac(c,j)    = data(c, 3*nlevgrnd + j)
         end do
         zw(c) = data(c, ncomp)
      end do
      deallocate(data)

      ! --- Index of the deepest layer that lies entirely above the water table
      allocate(jwt(bp%begc_all:bp%endc_all))
      do c = bp%begc_all, bp%endc_all
         jwt(c) = 0
         if (nlevbed(c) <= 0 .or. nlevbed(c) > nlevgrnd) cycle
         jwt(c) = nlevbed(c)
         do j = 1, nlevbed(c)
            if (zw(c) <= zi(c,j)) then
               jwt(c) = j - 1
               exit
            end if
         end do
      end do

      ! --- Fluxes per connection (positive from the up to the down column)
      nconn = conn%nconn
      allocate(qv_unsat(nconn, nlevgrnd)) ; qv_unsat(:,:) = 0._r8
      allocate(qv_sat(nconn))             ; qv_sat(:)     = 0._r8

      do iconn = 1, nconn
         cu  = conn%col_id_up(iconn)
         cd  = conn%col_id_dn(iconn)
         nbu = nlevbed(cu)
         nbd = nlevbed(cd)

         ! unsaturated flux (Eq. 15) in layers above the water table of both columns
         if (lateral_unsat_flow) then
            do j = 1, min(jwt(cu), jwt(cd))
               s1 = (h2osoi_vol(cu,j) + h2osoi_vol(cd,j)) / (watsat(cu,j) + watsat(cd,j))
               s1 = min(1._r8, max(s1, 0._r8))

               bswl = 0.5_r8*(bsw(cu,j) + bsw(cd,j))

               if (origflag == 1) then
                  impedl = 1._r8 - 0.5_r8*(fracice(cu,j) + fracice(cd,j))
               else
                  impedl = 10._r8**(-e_ice*(0.5_r8*(icefrac(cu,j) + icefrac(cd,j))))
               end if

               ! horizontal conductivity: geometric mean of the two columns
               hkl = lateral_hk_anisotropy * impedl * sqrt(hksat(cu,j)*hksat(cd,j)) &
                    * s1**(2._r8*bswl + 3._r8)

               ! gradient of (matric potential + elevation) along the land surface [-]
               smp_grad = (smp_l(cd,j) - smp_l(cu,j) + conn%dzg(iconn)*1.e3_r8) / (conn%dist(iconn)*1.e3_r8)

               face_area = conn%face_length(iconn) * 0.5_r8*(dz(cu,j) + dz(cd,j))   ! [m^2]

               qv_unsat(iconn, j) = -hkl * smp_grad * face_area
            end do
         end if

         ! saturated flux (Eq. 17) driven by the difference in water table head
         hsat_up = max(zi(cu,nbu) - zw(cu), 0._r8)   ! saturated thickness [m]
         hsat_dn = max(zi(cd,nbd) - zw(cd), 0._r8)

         trans = lateral_hk_anisotropy * sqrt(hksat(cu,nbu)*hksat(cd,nbd)) &
              * 0.5_r8*(hsat_up + hsat_dn)                                      ! [mm/s * m]

         ! head difference [m] = (elev_dn - zwt_dn) - (elev_up - zwt_up)
         qv_sat(iconn) = -trans * (conn%dzg(iconn) - zw(cd) + zw(cu)) / conn%dist(iconn) &
              * conn%face_length(iconn)
      end do

      ! --- Sum the fluxes for each owned column, in a decomposition-independent order
      qflx_lat_layer(bounds%begc:bounds%endc, :) = 0._r8
      qflx_sat(bounds%begc:bounds%endc)          = 0._r8
      zwt_lat(bounds%begc:bounds%endc)           = zw(bounds%begc:bounds%endc)

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         col_wf%qflx_lateral(c)       = 0._r8
         col_wf%qflx_lateral_unsat(c) = 0._r8
         col_wf%qflx_lateral_sat(c)   = 0._r8
         col_wf%qflx_lat_layer(c,:)   = 0._r8

         if (conn%col_nconn(c) == 0) cycle

         nb          = nlevbed(c)
         unsat_in(:) = 0._r8
         sat_in      = 0._r8

         do k = conn%col_conn_beg(c), conn%col_conn_beg(c) + conn%col_nconn(c) - 1
            iconn = conn%col_conn_id(k)
            sgn   = conn%col_conn_sign(k)    ! +1: c is the down column (up-to-down flux enters c)
            if (sgn > 0._r8) then
               area = conn%downarea(iconn)
            else
               area = conn%uparea(iconn)
            end if

            do j = 1, nb
               unsat_in(j) = unsat_in(j) + sgn*qv_unsat(iconn, j)/area
            end do
            sat_in = sat_in + sgn*qv_sat(iconn)/area
         end do

         do j = 1, nb
            qflx_lat_layer(c,j)          = unsat_in(j)
            col_wf%qflx_lat_layer(c,j)   = unsat_in(j)
            col_wf%qflx_lateral_unsat(c) = col_wf%qflx_lateral_unsat(c) + unsat_in(j)
         end do
         qflx_sat(c)                = sat_in
         col_wf%qflx_lateral_sat(c) = sat_in
         col_wf%qflx_lateral(c)     = -(col_wf%qflx_lateral_unsat(c) + sat_in)
      end do

      deallocate(jwt, zw, qv_unsat, qv_sat)

    end associate
#else
    call endrun(msg='ComputeLateralFlux: requires ELM built with -DMOAB_LATERAL'//errMsg(__FILE__, __LINE__))
#endif

  end subroutine ComputeLateralFlux

  !-----------------------------------------------------------------------
  subroutine ApplySaturatedLateralFlux(bounds, num_hydrologyc, filter_hydrologyc, dtime, &
       soilhydrology_vars, soilstate_vars, qflx_sat, zwt_lat)
    !
    ! !DESCRIPTION:
    ! Applies each column's net saturated lateral flux at the water table,
    ! after the Richards solve, following ELM's baseflow treatment (Drainage):
    !
    ! - Outflow is removed from the layer that contains the water table and
    !   then from the layers below it in turn; each layer gives at most its
    !   specific yield times its saturated thickness, and the water table
    !   drops accordingly. If the saturated zone cannot supply the outflow,
    !   the rest is taken from the unsaturated layers above, bottom up.
    ! - Inflow fills the layer that contains the water table and then the
    !   layers above it in turn, each up to its effective porosity, and the
    !   water table rises accordingly. Any inflow left when the column is full
    !   is added to the top layer, whose excess goes to surface runoff in
    !   Drainage.
    !
    ! The water table starts from zwt_lat, the one the fluxes were computed
    ! with. The applied flux is added to col_wf%qflx_lat_layer.
    !
    ! !USES:
    use elm_varcon , only : denice, watmin
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    real(r8)                 , intent(in)    :: dtime                ! time step [s]
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    real(r8)                 , intent(in)    :: qflx_sat(bounds%begc:) ! net saturated lateral flux into each column [mm H2O/s]
    real(r8)                 , intent(in)    :: zwt_lat(bounds%begc:)  ! water table depth used for the lateral fluxes [m]
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, j, jwt, nb
    real(r8) :: tot        ! flux still to apply [mm H2O] (positive into the column)
    real(r8) :: dw         ! water added to a layer [mm H2O]
    real(r8) :: s_y        ! specific yield [-]
    real(r8) :: avail      ! water a layer can give [mm H2O]
    real(r8) :: cap        ! water a layer can take [mm H2O]
    real(r8), parameter :: tol = 1.e-8_r8   ! [mm H2O]
    !-----------------------------------------------------------------------

    associate(                                         &
         zi         => col_pp%zi                     , & ! Input:  [real(r8) (:,:) ] interface depth (m)
         dz         => col_pp%dz                     , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         nlevbed    => col_pp%nlevbed                , & ! Input:  [integer  (:)   ] number of hydrologically active layers
         h2osoi_liq => col_ws%h2osoi_liq             , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice             , & ! Input:  [real(r8) (:,:) ] ice (kg/m2)
         watsat     => soilstate_vars%watsat_col     , & ! Input:  [real(r8) (:,:) ] porosity
         sucsat     => soilstate_vars%sucsat_col     , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
         bsw        => soilstate_vars%bsw_col        , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
         zwt        => soilhydrology_vars%zwt_col      & ! Output: [real(r8) (:)   ] water table depth (m)
         )

      do fc = 1, num_hydrologyc
         c   = filter_hydrologyc(fc)
         nb  = nlevbed(c)
         tot = qflx_sat(c)*dtime
         if (tot == 0._r8) cycle
         zwt(c) = zwt_lat(c)

         ! index of the deepest layer that lies entirely above the water table
         jwt = nb
         do j = 1, nb
            if (zwt(c) <= zi(c,j)) then
               jwt = j - 1
               exit
            end if
         end do

         if (tot < 0._r8) then
            ! outflow: drain the layers at and below the water table in turn
            do j = jwt+1, nb
               s_y = watsat(c,j) * (1._r8 - (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,j))**(-1._r8/bsw(c,j)))
               s_y = max(s_y, 0.02_r8)
               avail = min(s_y*(zi(c,j) - zwt(c))*1.e3_r8, max(h2osoi_liq(c,j) - watmin, 0._r8))
               dw = min(max(tot, -avail), 0._r8)
               h2osoi_liq(c,j) = h2osoi_liq(c,j) + dw
               col_wf%qflx_lat_layer(c,j) = col_wf%qflx_lat_layer(c,j) + dw/dtime
               tot = tot - dw
               if (tot >= 0._r8) then
                  zwt(c) = zwt(c) - dw/s_y/1.e3_r8
                  exit
               else
                  zwt(c) = zi(c,j)
               end if
            end do
            ! saturated zone exhausted: take the rest from the layers above
            do j = min(jwt, nb), 1, -1
               if (tot >= 0._r8) exit
               dw = min(max(tot, -max(h2osoi_liq(c,j) - watmin, 0._r8)), 0._r8)
               h2osoi_liq(c,j) = h2osoi_liq(c,j) + dw
               col_wf%qflx_lat_layer(c,j) = col_wf%qflx_lat_layer(c,j) + dw/dtime
               tot = tot - dw
            end do
            if (tot < -tol) then
               call endrun(msg='ApplySaturatedLateralFlux: lateral outflow exceeds the water in the column'// &
                    errMsg(__FILE__, __LINE__))
            end if
         else
            ! inflow: fill the layers at and above the water table in turn
            do j = min(jwt+1, nb), 1, -1
               cap = max(watsat(c,j) - h2osoi_ice(c,j)/(dz(c,j)*denice), 0.01_r8)*dz(c,j)*1.e3_r8
               cap = max(cap - h2osoi_liq(c,j), 0._r8)
               dw  = min(tot, cap)
               h2osoi_liq(c,j) = h2osoi_liq(c,j) + dw
               col_wf%qflx_lat_layer(c,j) = col_wf%qflx_lat_layer(c,j) + dw/dtime
               tot = tot - dw
               s_y = watsat(c,j) * (1._r8 - (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,j))**(-1._r8/bsw(c,j)))
               s_y = max(s_y, 0.02_r8)
               zwt(c) = max(zi(c,j-1), zwt(c) - dw/s_y/1.e3_r8)
               if (tot <= 0._r8) exit
            end do
            if (tot > 0._r8) then
               h2osoi_liq(c,1) = h2osoi_liq(c,1) + tot
               col_wf%qflx_lat_layer(c,1) = col_wf%qflx_lat_layer(c,1) + tot/dtime
            end if
         end if
      end do

    end associate

  end subroutine ApplySaturatedLateralFlux

  !-----------------------------------------------------------------------
  subroutine ThetaBasedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, soilhydrology_vars, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Diagnose the water table depth from soil moisture (as in the CLM5
    ! theta-based method): find the deepest layer, searching up from the
    ! bottom of the hydrologically active column, whose saturation is at or
    ! below sat_lev, and interpolate between it and the layer below.
    !
    ! !USES:
    use elm_varcon , only : denh2o, denice
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer             :: c, fc, k, nb
    !-----------------------------------------------------------------------

    associate(                                         &
         dz         => col_pp%dz                     , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         z          => col_pp%z                      , & ! Input:  [real(r8) (:,:) ] layer depth (m)
         zi         => col_pp%zi                     , & ! Input:  [real(r8) (:,:) ] interface depth (m)
         nlevbed    => col_pp%nlevbed                , & ! Input:  [integer  (:)   ] number of hydrologically active layers
         h2osoi_liq => col_ws%h2osoi_liq             , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice             , & ! Input:  [real(r8) (:,:) ] ice (kg/m2)
         h2osoi_vol => col_ws%h2osoi_vol             , & ! Output: [real(r8) (:,:) ] volumetric soil water (m3/m3)
         watsat     => soilstate_vars%watsat_col     , & ! Input:  [real(r8) (:,:) ] porosity
         zwt        => soilhydrology_vars%zwt_col      & ! Output: [real(r8) (:)   ] water table depth (m)
         )

      do fc = 1, num_hydrologyc
         c  = filter_hydrologyc(fc)
         nb = nlevbed(c)

         ! update h2osoi_vol from the bottom up to the first unsaturated layer
         do k = nb, 1, -1
            h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) + h2osoi_ice(c,k)/(dz(c,k)*denice)
            if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) exit
         end do

         zwt(c) = ThetaWaterTableDepth(nb, h2osoi_vol(c,1:nb), watsat(c,1:nb), z(c,1:nb), zi(c,0:nb))
      end do

    end associate

  end subroutine ThetaBasedWaterTable

  !-----------------------------------------------------------------------
  pure function ThetaWaterTableDepth(nb, vol, watsat, z, zi) result(zwt)
    !
    ! !DESCRIPTION:
    ! Water table depth diagnosed from soil moisture (as in the CLM5
    ! theta-based method): find the deepest layer, searching up from the
    ! bottom of the hydrologically active column, whose saturation is at or
    ! below sat_lev, and interpolate between it and the layer below. If the
    ! bottom layer is unsaturated, the water table is at the bottom of the
    ! column (no saturated zone).
    !
    ! !ARGUMENTS:
    integer  , intent(in) :: nb          ! number of hydrologically active layers
    real(r8) , intent(in) :: vol(:)      ! volumetric water content (1:nb) [m3/m3]
    real(r8) , intent(in) :: watsat(:)   ! porosity (1:nb)
    real(r8) , intent(in) :: z(:)        ! layer depth (1:nb) [m]
    real(r8) , intent(in) :: zi(0:)      ! interface depth (0:nb) [m]
    real(r8)              :: zwt         ! water table depth [m]
    !
    ! !LOCAL VARIABLES:
    integer  :: k, k_zwt
    logical  :: all_saturated
    real(r8) :: s1, s2, m, b
    !-----------------------------------------------------------------------

    k_zwt         = nb
    all_saturated = .true.
    do k = nb, 1, -1
       if (vol(k)/watsat(k) <= sat_lev) then
          k_zwt         = k
          all_saturated = .false.
          exit
       end if
    end do
    if (all_saturated) k_zwt = 1

    if (k_zwt == 1) then
       ! all layers below the first are saturated: water table at the bottom of layer 1
       zwt = zi(1)
    else if (k_zwt < nb) then
       ! interpolate between k_zwt and k_zwt+1
       s1  = vol(k_zwt  )/watsat(k_zwt  )
       s2  = vol(k_zwt+1)/watsat(k_zwt+1)
       m   = (z(k_zwt+1) - z(k_zwt))/(s2 - s1)
       b   = z(k_zwt+1) - m*s2
       zwt = max(0._r8, m*sat_lev + b)
    else
       ! bottom layer is unsaturated: water table at the bottom of the column
       zwt = zi(nb)
    end if

  end function ThetaWaterTableDepth

end module SoilLateralFlowMod
