module SoilLateralFlowMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Lateral subsurface flow between grid cells, following Qiu et al. (2024),
  ! "Development of inter-grid-cell lateral unsaturated and saturated flow
  ! model in the E3SM Land Model (v2.0)", GMD 17, 143-167.
  !
  ! Lateral fluxes are computed explicitly from the state at the start of the
  ! soil-water solve and enter the vertical (Zeng and Decker, 2009) Richards
  ! equation as layer source terms:
  !
  ! - Unsaturated flow (Eq. 15): for each soil layer that lies entirely above
  !   the water table in both columns, a Darcy flux driven by the difference
  !   in matric potential plus elevation between the two columns.
  ! - Saturated flow (Eq. 17): a Darcy flux driven by the difference in water
  !   table head (surface elevation minus water table depth), using a
  !   transmissivity based on the mean saturated thickness of the two columns.
  !   Each column's saturated inflow is spread over its saturated layers in
  !   proportion to their saturated thickness.
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
  public :: ComputeLateralFlux      ! lateral flux into each soil layer of owned columns
  public :: ThetaBasedWaterTable    ! diagnose the water table from soil moisture
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ComputeLateralFlux(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, soilstate_vars, qflx_lat_layer)
    !
    ! !DESCRIPTION:
    ! Computes the lateral subsurface flux into each soil layer [mm H2O/s] of
    ! the columns in filter_hydrologyc, and sets col_wf%qflx_lateral (outflow
    ! positive), qflx_lat_layer, qflx_lateral_unsat and qflx_lateral_sat.
    !
    ! Must be called on all MPI ranks (it performs a halo exchange), with one
    ! clump per MPI rank: 'bounds' (clump bounds) must cover all owned columns.
    ! Ghost columns are addressed through the processor bounds.
    !
#ifdef MOAB_LATERAL
    use decompMod               , only : get_proc_bounds
    use elm_varctl              , only : lateral_unsat_flow, lateral_hk_anisotropy
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
    real(r8)                 , intent(out)   :: qflx_lat_layer(bounds%begc:,1:)  ! lateral flux into each soil layer [mm H2O/s]
    !
#ifdef MOAB_LATERAL
    ! !LOCAL VARIABLES:
    type(bounds_type)     :: bp              ! processor bounds (include ghost columns)
    integer               :: fc, c, j, k, iconn, cu, cd, nbu, nbd, nb, nconn, ncomp
    real(r8)              :: s1, bswl, hkl, impedl, smp_grad, face_area, area
    real(r8)              :: hsat_up, hsat_dn, trans, sgn, hsat_c, sat_part
    real(r8), allocatable :: data(:,:)
    integer , allocatable :: jwt(:)
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

      ! --- Copy the state of ghost columns from the owning ranks (one MPI round)
      ! Components: h2osoi_vol, smp_l, fracice, icefrac (each 1:nlevgrnd), zwt
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
         data(c, ncomp) = zwt(c)
      end do

      call NatVegColumnRealDataHaloExchange(bp, 'lateral_flow_state', data)

      do c = bp%endc + 1, bp%endc_all
         do j = 1, nlevgrnd
            h2osoi_vol(c,j) = data(c, 0*nlevgrnd + j)
            smp_l(c,j)      = data(c, 1*nlevgrnd + j)
            fracice(c,j)    = data(c, 2*nlevgrnd + j)
            icefrac(c,j)    = data(c, 3*nlevgrnd + j)
         end do
         zwt(c) = data(c, ncomp)
      end do
      deallocate(data)

      ! --- Index of the deepest layer that lies entirely above the water table
      allocate(jwt(bp%begc_all:bp%endc_all))
      do c = bp%begc_all, bp%endc_all
         jwt(c) = 0
         if (nlevbed(c) <= 0 .or. nlevbed(c) > nlevgrnd) cycle
         jwt(c) = nlevbed(c)
         do j = 1, nlevbed(c)
            if (zwt(c) <= zi(c,j)) then
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
         hsat_up = max(zi(cu,nbu) - zwt(cu), 0._r8)   ! saturated thickness [m]
         hsat_dn = max(zi(cd,nbd) - zwt(cd), 0._r8)

         trans = lateral_hk_anisotropy * sqrt(hksat(cu,nbu)*hksat(cd,nbd)) &
              * 0.5_r8*(hsat_up + hsat_dn)                                      ! [mm/s * m]

         ! head difference [m] = (elev_dn - zwt_dn) - (elev_up - zwt_up)
         qv_sat(iconn) = -trans * (conn%dzg(iconn) - zwt(cd) + zwt(cu)) / conn%dist(iconn) &
              * conn%face_length(iconn)
      end do

      ! --- Sum the fluxes for each owned column, in a decomposition-independent order
      qflx_lat_layer(bounds%begc:bounds%endc, :) = 0._r8

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

         ! spread the saturated inflow over the saturated layers of the column
         hsat_c = zi(c,nb) - zwt(c)
         if (hsat_c > 0._r8) then
            do j = 1, nb
               sat_part = max(0._r8, zi(c,j) - max(zi(c,j-1), zwt(c)))
               qflx_lat_layer(c,j) = qflx_lat_layer(c,j) + sat_in * sat_part/hsat_c
            end do
         else
            qflx_lat_layer(c,nb) = qflx_lat_layer(c,nb) + sat_in
         end if

         do j = 1, nb
            qflx_lat_layer(c,j) = qflx_lat_layer(c,j) + unsat_in(j)
            col_wf%qflx_lateral_unsat(c) = col_wf%qflx_lateral_unsat(c) + unsat_in(j)
         end do
         col_wf%qflx_lateral_sat(c) = sat_in

         do j = 1, nb
            col_wf%qflx_lat_layer(c,j) = qflx_lat_layer(c,j)
            col_wf%qflx_lateral(c)     = col_wf%qflx_lateral(c) - qflx_lat_layer(c,j)
         end do
      end do

      deallocate(jwt, qv_unsat, qv_sat)

    end associate
#else
    call endrun(msg='ComputeLateralFlux: requires ELM built with -DMOAB_LATERAL'//errMsg(__FILE__, __LINE__))
#endif

  end subroutine ComputeLateralFlux

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
    integer             :: c, fc, k, k_zwt, nb
    logical             :: all_saturated
    real(r8)            :: s1, s2, m, b
    real(r8), parameter :: sat_lev = 0.96_r8   ! saturation level that defines the water table
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

         k_zwt         = nb
         all_saturated = .true.
         do k = nb, 1, -1
            h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) + h2osoi_ice(c,k)/(dz(c,k)*denice)
            if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then
               k_zwt         = k
               all_saturated = .false.
               exit
            end if
         end do
         if (all_saturated) k_zwt = 1

         if (k_zwt == 1) then
            ! all layers below the first are saturated: water table at the bottom of layer 1
            zwt(c) = zi(c,1)
         else if (k_zwt < nb) then
            ! interpolate between k_zwt and k_zwt+1
            s1 = h2osoi_vol(c,k_zwt  )/watsat(c,k_zwt  )
            s2 = h2osoi_vol(c,k_zwt+1)/watsat(c,k_zwt+1)
            m  = (z(c,k_zwt+1) - z(c,k_zwt))/(s2 - s1)
            b  = z(c,k_zwt+1) - m*s2
            zwt(c) = max(0._r8, m*sat_lev + b)
         else
            ! bottom layer is unsaturated: water table at the bottom of the column
            zwt(c) = zi(c,nb)
         end if
      end do

    end associate

  end subroutine ThetaBasedWaterTable

end module SoilLateralFlowMod
