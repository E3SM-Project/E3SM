module SoilWaterMovementMod

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  ! module contains different subroutines to couple soil and root water interactions
  !
  ! created by Jinyun Tang, Mar 12, 2014
  ! added variable DTB option for Zeng-Decker, Michael A. Brunke, Aug. 25, 2016
  !
  use ColumnDataType    , only : col_es, col_ws, col_wf
  use VegetationDataType, only : veg_wf
  !
  implicit none
  save 
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilWater            ! Calculate soil hydrology   
  public :: init_soilwater_movement
  public :: Compute_EffecRootFrac_And_VertTranSink
!  public :: Compute_EffecRootFrac_And_VertTranSink_HydStress
  !
  ! !PUBLIC DATA MEMBERS:
  logical, public :: zengdecker_2009_with_var_soil_thick
  !
  ! !PRIVATE DATA MEMBERS:
  integer, parameter, public :: zengdecker_2009 = 0
  integer, parameter, public :: vsfm = 1
  integer, public :: soilroot_water_method     !0: use the Zeng and deck method, this will be readin from namelist in the future

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine init_soilwater_movement()
    !
    !DESCRIPTION
    !specify method for doing soil&root water interactions
    !
    use elm_varctl, only : use_vsfm, use_var_soil_thick, use_hydrstress
    use spmdMod,    only : mpicom, MPI_LOGICAL
    use shr_sys_mod,only : shr_sys_abort
    ! !ARGUMENTS:
    implicit none
    integer :: ier ! error status
    !------------------------------------------------------------------------------

    soilroot_water_method = zengdecker_2009
    zengdecker_2009_with_var_soil_thick = .false.

    ! GB-FIX-ME: The call to control_spmd() [in subroutine control_init()] before
    !            call to init_hydrology() would avoid the mpi broadcast
    call mpi_bcast (use_vsfm, 1, MPI_LOGICAL, 0, mpicom, ier)
    if (use_vsfm) soilroot_water_method = vsfm

    call mpi_bcast (use_var_soil_thick, 1, MPI_LOGICAL, 0, mpicom, ier)
    if (use_var_soil_thick .and. soilroot_water_method .eq. zengdecker_2009) then
       zengdecker_2009_with_var_soil_thick = .true.
    end if
    if (use_var_soil_thick .and. soilroot_water_method .ne. zengdecker_2009) then
       call shr_sys_abort('ERROR: use_var_soil_thick not supported with anything but zengdecker_2009 at this time.')
    end if

  end subroutine init_soilwater_movement

  !-----------------------------------------------------------------------
  subroutine SoilWater(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
       waterflux_vars, waterstate_vars, temperature_vars, soil_water_retention_curve)
    !
    ! DESCRIPTION
    ! select one subroutine to do the soil and root water coupling
    !
    !USES
    use elm_varctl                 , only : use_betr
    use elm_varctl                 , only : use_var_soil_thick
    use shr_kind_mod               , only : r8 => shr_kind_r8
    use elm_varpar                 , only : nlevsoi    
    use decompMod                  , only : bounds_type   
    use abortutils                 , only : endrun   
    use SoilHydrologyType          , only : soilhydrology_type
    use SoilStateType              , only : soilstate_type
    use TemperatureType            , only : temperature_type
    use WaterFluxType              , only : waterflux_type
    use WaterStateType             , only : waterstate_type
    use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
    use elm_varcon                 , only : denh2o, denice, watmin
    use ColumnType                 , only : col_pp
    use ExternalModelConstants     , only : EM_VSFM_SOIL_HYDRO_STAGE
    use ExternalModelConstants     , only : EM_ID_VSFM
    use ExternalModelInterfaceMod  , only : EMI_Driver
    use clm_time_manager           , only : get_step_size, get_nstep
    !
    ! !ARGUMENTS:
    implicit none     
    type(bounds_type)        , intent(in)    :: bounds                ! bounds
    integer                  , intent(in)    :: num_hydrologyc        ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:)  ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc            ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)      ! column filter for urban points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    character(len=32)                        :: subname = 'SoilWater'       ! subroutine name
    real(r8)                                 :: xs(bounds%begc:bounds%endc) !excess soil water above urban ponding limit
    integer                                  :: nlevbed                     ! number of layers to bedrock

    integer  :: fc, c, j
    
    
    !------------------------------------------------------------------------------
    associate(                                                         &
      wa                 =>    soilhydrology_vars%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)
      dz                 =>    col_pp%dz                                , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)    
      zwt                =>    soilhydrology_vars%zwt_col            , & ! Input:  [real(r8) (:)   ]  water table depth (m)
      nlev2bed           =>    col_pp%nlevbed                           , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     
      h2osoi_ice         =>    col_ws%h2osoi_ice        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
      h2osoi_vol         =>    col_ws%h2osoi_vol        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
      h2osoi_liq         =>    col_ws%h2osoi_liq          & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
    )

    select case(soilroot_water_method)

    case (zengdecker_2009)

       call soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc, &
            num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
            waterflux_vars, waterstate_vars, temperature_vars, soil_water_retention_curve)

    case (vsfm)
#ifdef USE_PETSC_LIB

       call Prepare_Data_for_EM_VSFM_Driver(bounds, num_hydrologyc, filter_hydrologyc, &
            soilhydrology_vars, soilstate_vars, &
            waterflux_vars, waterstate_vars, temperature_vars)

       call EMI_Driver(EM_ID_VSFM, EM_VSFM_SOIL_HYDRO_STAGE, dt = get_step_size()*1.0_r8, &
            number_step = get_nstep(), &
            clump_rank  = bounds%clump_index, &
            num_hydrologyc=num_hydrologyc, filter_hydrologyc=filter_hydrologyc, &
            soilhydrology_vars=soilhydrology_vars, soilstate_vars=soilstate_vars, &
            waterflux_vars=waterflux_vars, waterstate_vars=waterstate_vars, &
            temperature_vars=temperature_vars)
#endif
    case default

       call endrun(subname // ':: a SoilWater implementation must be specified!')          

    end select

    if(use_betr)then
    !a work around of the negative liquid water embarrassment, which is
    !critical for a meaningufl tracer transport in betr. Jinyun Tang, Jan 14, 2015

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       nlevbed = nlev2bed(c)
       do j = 1, nlevbed-1
          if (h2osoi_liq(c,j) < 0._r8) then
             xs(c) = watmin - h2osoi_liq(c,j)
          else
             xs(c) = 0._r8
          end if
          h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
          h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
       end do
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       j = nlev2bed(c)
       if (h2osoi_liq(c,j) < watmin) then
          xs(c) = watmin-h2osoi_liq(c,j)
        else
          xs(c) = 0._r8
       end if
       h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
       wa(c) = wa(c) - xs(c)
    end do
    
    !update volumetric soil moisture for bgc calculation
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       nlevbed = nlev2bed(c)
       do j = 1, nlevbed
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) &
                            + h2osoi_ice(c,j)/(dz(c,j)*denice)
       enddo
    enddo
    endif
  end associate

  end subroutine SoilWater

  !-----------------------------------------------------------------------   
  subroutine soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
       waterflux_vars, waterstate_vars, temperature_vars, soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Soil hydrology
    ! Soil moisture is predicted from a 10-layer model (as with soil
    ! temperature), in which the vertical soil moisture transport is governed
    ! by infiltration, runoff, gradient diffusion, gravity, and root
    ! extraction through canopy transpiration.  The net water applied to the
    ! surface layer is the snowmelt plus precipitation plus the throughfall
    ! of canopy dew minus surface runoff and evaporation.
    ! CLM3.5 uses a zero-flow bottom boundary condition.
    !
    ! The vertical water flow in an unsaturated porous media is described by
    ! Darcy's law, and the hydraulic conductivity and the soil negative
    ! potential vary with soil water content and soil texture based on the work
    ! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
    ! integrated over the layer thickness, in which the time rate of change in
    ! water mass must equal the net flow across the bounding interface, plus the
    ! rate of internal source or sink. The terms of water flow across the layer
    ! interfaces are linearly expanded by using first-order Taylor expansion.
    ! The equations result in a tridiagonal system equation.
    !
    ! Note: length units here are all millimeter
    ! (in temperature subroutine uses same soil layer
    ! structure required but lengths are m)
    !
    ! Richards equation:
    !
    ! d wat      d     d wat d psi
    ! ----- = - -- [ k(----- ----- - 1) ] + S
    !   dt      dz       dz  d wat
    !
    ! where: wat = volume of water per volume of soil (mm**3/mm**3)
    ! psi = soil matrix potential (mm)
    ! dt  = time step (s)
    ! z   = depth (mm)
    ! dz  = thickness (mm)
    ! qin = inflow at top (mm h2o /s)
    ! qout= outflow at bottom (mm h2o /s)
    ! s   = source/sink flux (mm h2o /s)
    ! k   = hydraulic conductivity (mm h2o /s)
    !
    !                       d qin                  d qin
    ! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
    !                       d wat(j-1)             d wat(j)
    !                ==================|=================
    !                                  < qin
    !
    !                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
    !
    !                                  > qout
    !                ==================|=================
    !                        d qout               d qout
    ! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
    !                        d wat(j)             d wat(j+1)
    !
    !
    ! Solution: linearize k and psi about d wat and use tridiagonal
    ! system of equations to solve for d wat,
    ! where for layer j
    !
    !
    ! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
    !
    ! !USES:
    use elm_varctl           , only : use_var_soil_thick
    use shr_kind_mod         , only : r8 => shr_kind_r8     
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use decompMod            , only : bounds_type        
    use elm_varcon           , only : wimp,grav,hfus,tfrz
    use elm_varcon           , only : e_ice,denh2o, denice
    use elm_varpar           , only : nlevsoi, max_patch_per_col, nlevgrnd
    use clm_time_manager     , only : get_step_size
    use column_varcon        , only : icol_roof, icol_road_imperv
    use TridiagonalMod       , only : Tridiagonal
    use abortutils           , only : endrun     
    use SoilStateType        , only : soilstate_type
    use SoilHydrologyType    , only : soilhydrology_type
    use TemperatureType      , only : temperature_type
    use WaterFluxType        , only : waterflux_type
    use WaterStateType       , only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use VegetationType       , only : veg_pp
    use ColumnType           , only : col_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    type(soilstate_type)    , intent(inout) :: soilstate_vars
    type(waterflux_type)    , intent(inout) :: waterflux_vars
    type(waterstate_type)   , intent(inout) :: waterstate_vars
    type(temperature_type)  , intent(in)    :: temperature_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,fc,j                                     ! do loop indices
    integer  :: nlevbed                                      ! number of layers to bedrock
    integer  :: jtop(bounds%begc:bounds%endc)                ! top level at each column
    integer  :: jbot(bounds%begc:bounds%endc)                ! bottom level at each column
    real(r8) :: dtime                                        ! land model time step (sec)
    real(r8) :: hk(bounds%begc:bounds%endc,1:nlevgrnd)        ! hydraulic conductivity [mm h2o/s]
    real(r8) :: dhkdw(bounds%begc:bounds%endc,1:nlevgrnd)     ! d(hk)/d(vol_liq)
    real(r8) :: amx(bounds%begc:bounds%endc,1:nlevgrnd+1)     ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(bounds%begc:bounds%endc,1:nlevgrnd+1)     ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(bounds%begc:bounds%endc,1:nlevgrnd+1)     ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(bounds%begc:bounds%endc,1:nlevgrnd+1)     ! "r" forcing term of tridiagonal matrix
    real(r8) :: zmm(bounds%begc:bounds%endc,1:nlevgrnd+1)     ! layer depth [mm]
    real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevgrnd+1)    ! layer thickness [mm]
    real(r8) :: den                                          ! used in calculating qin, qout
    real(r8) :: dqidw0(bounds%begc:bounds%endc,1:nlevgrnd+1)  ! d(qin)/d(vol_liq(i-1))
    real(r8) :: dqidw1(bounds%begc:bounds%endc,1:nlevgrnd+1)  ! d(qin)/d(vol_liq(i))
    real(r8) :: dqodw1(bounds%begc:bounds%endc,1:nlevgrnd+1)  ! d(qout)/d(vol_liq(i))
    real(r8) :: dqodw2(bounds%begc:bounds%endc,1:nlevgrnd+1)  ! d(qout)/d(vol_liq(i+1))
    real(r8) :: dsmpdw(bounds%begc:bounds%endc,1:nlevgrnd+1)  ! d(smp)/d(vol_liq)
    real(r8) :: num                                          ! used in calculating qin, qout
    real(r8) :: qin(bounds%begc:bounds%endc,1:nlevgrnd+1)     ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout(bounds%begc:bounds%endc,1:nlevgrnd+1)    ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: s_node                                       ! soil wetness
    real(r8) :: s1                                           ! "s" at interface of layer
    real(r8) :: s2                                           ! k*s**(2b+2)
    real(r8) :: smp(bounds%begc:bounds%endc,1:nlevgrnd)       ! soil matrix potential [mm]
    real(r8) :: sdamp                                        ! extrapolates soiwat dependence of evaporation
    integer  :: pi                                           ! pft index
    real(r8) :: temp(bounds%begc:bounds%endc)                ! accumulator for rootr weighting
    integer  :: jwt(bounds%begc:bounds%endc)                 ! index of the soil layer right above the water table (-)
    real(r8) :: smp1,dsmpdw1,wh,wh_zwt,ka
    real(r8) :: dwat2(bounds%begc:bounds%endc,1:nlevgrnd+1)
    real(r8) :: dzq                                          ! used in calculating qin, qout (difference in equilbirium matric potential)
    real(r8) :: zimm(bounds%begc:bounds%endc,0:nlevgrnd)      ! layer interface depth [mm]
    real(r8) :: zq(bounds%begc:bounds%endc,1:nlevgrnd+1)      ! equilibrium matric potential for each layer [mm]
    real(r8) :: vol_eq(bounds%begc:bounds%endc,1:nlevgrnd+1)  ! equilibrium volumetric water content
    real(r8) :: tempi                                        ! temp variable for calculating vol_eq
    real(r8) :: temp0                                        ! temp variable for calculating vol_eq
    real(r8) :: voleq1                                       ! temp variable for calculating vol_eq
    real(r8) :: zwtmm(bounds%begc:bounds%endc)               ! water table depth [mm]
    real(r8) :: imped(bounds%begc:bounds%endc,1:nlevgrnd)             
    real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevgrnd)
    real(r8) :: z_mid
    real(r8) :: vwc_zwt(bounds%begc:bounds%endc)
    real(r8) :: vwc_liq(bounds%begc:bounds%endc,1:nlevgrnd+1) ! liquid volumetric water content
    real(r8) :: smp_grad(bounds%begc:bounds%endc,1:nlevgrnd+1)
    real(r8) :: dsmpds                                       !temporary variable
    real(r8) :: dhkds                                        !temporary variable
    real(r8) :: hktmp                                        !temporary variable
    !-----------------------------------------------------------------------

    associate(& 
         z                 =>    col_pp%z                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         zi                =>    col_pp%zi                             , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
         dz                =>    col_pp%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                             
         nlev2bed          =>    col_pp%nlevbed                        , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     

         origflag          =>    soilhydrology_vars%origflag        , & ! Input:  constant
         qcharge           =>    soilhydrology_vars%qcharge_col     , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)                      
         zwt               =>    soilhydrology_vars%zwt_col         , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
         fracice           =>    soilhydrology_vars%fracice_col     , & ! Input:  [real(r8) (:,:) ]  fractional impermeability (-)                   
         icefrac           =>    soilhydrology_vars%icefrac_col     , & ! Input:  [real(r8) (:,:) ]  fraction of ice                                 
         hkdepth           =>    soilhydrology_vars%hkdepth_col     , & ! Input:  [real(r8) (:)   ]  decay factor (m)                                  

         smpmin            =>    soilstate_vars%smpmin_col          , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)        
         watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
         hksat             =>    soilstate_vars%hksat_col           , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
         sucsat            =>    soilstate_vars%sucsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         eff_porosity      =>    soilstate_vars%eff_porosity_col    , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         
         smp_l             =>    soilstate_vars%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]                      
         hk_l              =>    soilstate_vars%hk_l_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s)                   
         rootr_pft         =>    soilstate_vars%rootr_patch         , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer  

         h2osoi_ice        =>    col_ws%h2osoi_ice     , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_liq        =>    col_ws%h2osoi_liq     , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_vol        =>    col_ws%h2osoi_vol     , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

         qflx_deficit      =>    col_wf%qflx_deficit    , & ! Input:  [real(r8) (:)   ]  water deficit to keep non-negative liquid water content
         qflx_infl         =>    col_wf%qflx_infl       , & ! Input:  [real(r8) (:)   ]  infiltration (mm H2O /s)                          
         qflx_rootsoi_col  =>    col_wf%qflx_rootsoi    , & ! Input: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
         t_soisno          =>    col_es%t_soisno        & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       
         )

      ! Get time step
      
      dtime = get_step_size()

      ! Because the depths in this routine are in mm, use local
      ! variable arrays instead of pointers

      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        nlevbed = nlev2bed(c)
        do j = 1, nlevbed
            zmm(c,j) = z(c,j)*1.e3_r8
            dzmm(c,j) = dz(c,j)*1.e3_r8
            zimm(c,j) = zi(c,j)*1.e3_r8

            ! calculate icefrac up here
            vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
            icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
            vwc_liq(c,j) = max(h2osoi_liq(c,j),1.0e-6_r8)/(dz(c,j)*denh2o)
         end do
      end do

      do fc = 1, num_hydrologyc 
         c = filter_hydrologyc(fc)
         zimm(c,0) = 0.0_r8
         zwtmm(c)  = zwt(c)*1.e3_r8
      end do



      !compute jwt index
      ! The layer index of the first unsaturated layer, i.e., the layer right above
      ! the water table

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         jwt(c) = nlevbed
         ! allow jwt to equal zero when zwt is in top layer
         do j = 1,nlevbed
            if (use_var_soil_thick) then
               if (zwt(c) <= zi(c,j) .and. zwt(c) < zi(c,nlevbed)) then
                  jwt(c) = j-1
                  exit
               end if
            else
               if (zwt(c) <= zi(c,j)) then
                  jwt(c) = j-1
                  exit
               end if
            end if
         enddo

         ! compute vwc at water table depth (mainly for case when t < tfrz)
         !     this will only be used when zwt is below the soil column
         vwc_zwt(c) = watsat(c,nlevbed)
         if(t_soisno(c,jwt(c)+1) < tfrz) then
            vwc_zwt(c) = vwc_liq(c,nlevbed)
            do j = nlevbed,nlevgrnd
               if(zwt(c) <= zi(c,j)) then
                  smp1 = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                  !smp1 = max(0._r8,smp1)
                  smp1 = max(sucsat(c,nlevsoi),smp1)
                  vwc_zwt(c) = watsat(c,nlevsoi)*(smp1/sucsat(c,nlevbed))**(-1._r8/bsw(c,nlevsoi))
                  ! for temperatures close to tfrz, limit vwc to total water content 
                  vwc_zwt(c) = min(vwc_zwt(c), 0.5*(watsat(c,nlevbed) + h2osoi_vol(c,nlevbed)) )
                  exit
               endif
            enddo
         endif
      end do

      ! calculate the equilibrium water content based on the water table depth

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         do j = 1, nlevbed
            if ((zwtmm(c) <= zimm(c,j-1))) then 
               vol_eq(c,j) = watsat(c,j)

               ! use the weighted average from the saturated part (depth > wtd) and the equilibrium solution for the
               ! rest of the layer, the equilibrium solution is based on Clapp-Hornberg parameterization
               ! and no extension to full range swrc is needed

            else if ((zwtmm(c) .lt. zimm(c,j)) .and. (zwtmm(c) .gt. zimm(c,j-1))) then
               tempi = 1.0_r8
               temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               voleq1 = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j-1))*(tempi-temp0)
               vol_eq(c,j) = (voleq1*(zwtmm(c)-zimm(c,j-1)) + watsat(c,j)*(zimm(c,j)-zwtmm(c)))/(zimm(c,j)-zimm(c,j-1))
               vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
               vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
            else
               tempi = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               vol_eq(c,j) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zimm(c,j)-zimm(c,j-1))*(tempi-temp0)
               vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
               vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
            endif
            zq(c,j) = -sucsat(c,j)*(max(vol_eq(c,j)/watsat(c,j),0.01_r8))**(-bsw(c,j))
            zq(c,j) = max(smpmin(c), zq(c,j))
         end do
      end do

      ! If water table is below soil column calculate zq for the 11th layer
      do fc=1, num_hydrologyc
         c = filter_hydrologyc(fc)
         j = nlev2bed(c)
         if(jwt(c) == nlevbed) then 
            tempi = 1._r8
            temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
            vol_eq(c,j+1) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j))*(tempi-temp0)
            vol_eq(c,j+1) = max(vol_eq(c,j+1),0.0_r8)
            vol_eq(c,j+1) = min(watsat(c,j),vol_eq(c,j+1))
            zq(c,j+1) = -sucsat(c,j)*(max(vol_eq(c,j+1)/watsat(c,j),0.01_r8))**(-bsw(c,j))
            zq(c,j+1) = max(smpmin(c), zq(c,j+1))
         end if
      end do

      ! Hydraulic conductivity and soil matric potential and their derivatives

      sdamp = 0._r8
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         do j = 1, nlevbed
            ! compute hydraulic conductivity based on liquid water content only

            if (origflag == 1) then
               s1 = 0.5_r8*(h2osoi_vol(c,j) + h2osoi_vol(c,min(nlevsoi, j+1))) / &
                    (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
            else
               s1 = 0.5_r8*(vwc_liq(c,j) + vwc_liq(c,min(nlevsoi, j+1))) / &
                    (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
            endif
            s1 = min(1._r8, s1)
            s2 = hksat(c,j)*s1**(2._r8*bsw(c,j)+2._r8)

            ! replace fracice with impedance factor, as in zhao 97,99
            if (origflag == 1) then
               imped(c,j)=(1._r8-0.5_r8*(fracice(c,j)+fracice(c,min(nlevsoi, j+1))))
            else
               imped(c,j)=10._r8**(-e_ice*(0.5_r8*(icefrac(c,j)+icefrac(c,min(nlevsoi, j+1)))))
            endif
            hk(c,j) = imped(c,j)*s1*s2
            dhkdw(c,j) = imped(c,j)*(2._r8*bsw(c,j)+3._r8)*s2* &
                 (1._r8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))

            !compute un-restricted hydraulic conductivity
            !call soil_water_retention_curve%soil_hk(hksat(c,j), imped(c,j), s1, bsw(c,j), hktmp, dhkds)
            !if(hktmp/=hk(c,j))write(10,*)'diff',hktmp,hk(c,j)
            !    call endrun('bad in hk')
            !endif    
            !apply ice impedance
            !hk(c,j) = imped(c,j)*hk(c,j)          
            !dhkdw(c,j) = imped(c,j) * dhkds * (1._r8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))


            ! compute matric potential and derivative based on liquid water content only
            if (origflag == 1) then
               s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
            else
               s_node = max(vwc_liq(c,j)/watsat(c,j), 0.01_r8)
            endif
            s_node = min(1.0_r8, s_node)

            !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp(c,j), dsmpds)

            smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
            smp(c,j) = max(smpmin(c), smp(c,j))
            !do not turn on the line below, which will cause bit to bit error, jyt, 2014 Mar 6
            !dsmpdw(c,j) = dsmpds/watsat(c,j)

            if (origflag == 1) then             
               dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))
            else
               dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/vwc_liq(c,j)
            endif

            smp_l(c,j) = smp(c,j)
            hk_l(c,j) = hk(c,j)

         end do
      end do

      ! aquifer (11th) layer
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         zmm(c,nlevbed+1) = 0.5*(1.e3_r8*zwt(c) + zmm(c,nlevbed))
         if(jwt(c) < nlevbed) then
            dzmm(c,nlevbed+1) = dzmm(c,nlevbed)
         else
            dzmm(c,nlevbed+1) = (1.e3_r8*zwt(c) - zmm(c,nlevbed))
         end if
      end do

      ! Set up r, a, b, and c vectors for tridiagonal solution

      ! Node j=1 (top)

      j = 1
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qin(c,j)    = qflx_infl(c)
         den    = (zmm(c,j+1)-zmm(c,j))
         dzq    = (zq(c,j+1)-zq(c,j))
         num    = (smp(c,j+1)-smp(c,j)) - dzq
         qout(c,j)   = -hk(c,j)*num/den
         dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
         dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
         rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)
         amx(c,j) =  0._r8
         bmx(c,j) =  dzmm(c,j)*(sdamp+1._r8/dtime) + dqodw1(c,j)
         cmx(c,j) =  dqodw2(c,j)
      end do

      ! Nodes j=2 to j=nlevsoi-1

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         do j = 2, nlevbed - 1
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            den    = (zmm(c,j+1)-zmm(c,j))
            dzq    = (zq(c,j+1)-zq(c,j))
            num    = (smp(c,j+1)-smp(c,j)) - dzq
            qout(c,j)   = -hk(c,j)*num/den
            dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
            dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
            rmx(c,j)    =  qin(c,j) - qout(c,j) -  qflx_rootsoi_col(c,j)
            amx(c,j)    = -dqidw0(c,j)
            bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j)    =  dqodw2(c,j)
         end do
      end do

      ! Node j=nlevsoi (bottom)

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         j = nlevbed
         if(j > jwt(c)) then !water table is in soil column
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            qout(c,j)   =  0._r8
            dqodw1(c,j) =  0._r8
            rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)
            amx(c,j)    = -dqidw0(c,j)
            bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j)    =  0._r8

            ! next set up aquifer layer; hydrologically inactive
            rmx(c,j+1) = 0._r8
            amx(c,j+1) = 0._r8
            bmx(c,j+1) = dzmm(c,j+1)/dtime
            cmx(c,j+1) = 0._r8
         else ! water table is below soil column

            ! compute aquifer soil moisture as average of layer 10 and saturation
            if(origflag == 1) then
               s_node = max(0.5*(1.0_r8+h2osoi_vol(c,j)/watsat(c,j)), 0.01_r8)
            else
               s_node = max(0.5*((vwc_zwt(c)+vwc_liq(c,j))/watsat(c,j)), 0.01_r8)
            endif
            s_node = min(1.0_r8, s_node)

            ! compute smp for aquifer layer
            !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp1, dsmpds)
            smp1 = -sucsat(c,j)*s_node**(-bsw(c,j))
            smp1 = max(smpmin(c), smp1)

            ! compute dsmpdw for aquifer layer
            !dsmpdw1 = dsmpds/watsat(c,j)
            dsmpdw1 = -bsw(c,j)*smp1/(s_node*watsat(c,j))

            ! first set up bottom layer of soil column
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            den    = (zmm(c,j+1)-zmm(c,j))
            dzq    = (zq(c,j+1)-zq(c,j))
            num    = (smp1-smp(c,j)) - dzq
            if (use_var_soil_thick) then
               qout(c,j) = 0._r8
               dqodw1(c,j) = 0._r8
               dqodw2(c,j) = 0._r8
            else
               qout(c,j)   = -hk(c,j)*num/den
               dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
               dqodw2(c,j) = -( hk(c,j)*dsmpdw1 + num*dhkdw(c,j))/den
            end if

            rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)
            amx(c,j) = -dqidw0(c,j)
            bmx(c,j) =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j) =  dqodw2(c,j)

            ! next set up aquifer layer; den/num unchanged, qin=qout
            qin(c,j+1)    = qout(c,j)
            dqidw0(c,j+1) = -(-hk(c,j)*dsmpdw(c,j) + num*dhkdw(c,j))/den
            dqidw1(c,j+1) = -( hk(c,j)*dsmpdw1   + num*dhkdw(c,j))/den
            qout(c,j+1)   =  0._r8  ! zero-flow bottom boundary condition
            dqodw1(c,j+1) =  0._r8  ! zero-flow bottom boundary condition
            if (use_var_soil_thick) then
               rmx(c,j+1) = 0._r8
               amx(c,j+1) = 0._r8
               bmx(c,j+1) = dzmm(c,j+1)/dtime
               cmx(c,j+1) = 0._r8
            else
               rmx(c,j+1) =  qin(c,j+1) - qout(c,j+1)
               amx(c,j+1) = -dqidw0(c,j+1)
               bmx(c,j+1) =  dzmm(c,j+1)/dtime - dqidw1(c,j+1) + dqodw1(c,j+1)
               cmx(c,j+1) =  0._r8
            end if
         endif
      end do

      ! Solve for dwat

      jtop(bounds%begc : bounds%endc) = 1
      ! Determination of how many layers (nlev2bed) to do for the tridiagonal
      ! at each column
      if (use_var_soil_thick) then
      	 do fc = 1,num_hydrologyc
            c = filter_hydrologyc(fc)
            jbot(c) = nlev2bed(c)
      	 end do
         call Tridiagonal(bounds, 1, nlevgrnd+1, &
              jtop(bounds%begc:bounds%endc),     &
              jbot(bounds%begc:bounds%endc),     &
              num_hydrologyc, filter_hydrologyc, &
              amx(bounds%begc:bounds%endc, :),   &
              bmx(bounds%begc:bounds%endc, :),   &
              cmx(bounds%begc:bounds%endc, :),   &
              rmx(bounds%begc:bounds%endc, :),   &
              dwat2(bounds%begc:bounds%endc, :) )
      else
         call Tridiagonal(bounds, 1, nlevsoi+1, &
              jtop(bounds%begc:bounds%endc), &
              num_hydrologyc, filter_hydrologyc, &
              amx(bounds%begc:bounds%endc, :), &
              bmx(bounds%begc:bounds%endc, :), &
              cmx(bounds%begc:bounds%endc, :), &
              rmx(bounds%begc:bounds%endc, :), &
              dwat2(bounds%begc:bounds%endc, :) )
      end if

      ! Renew the mass of liquid water
      ! also compute qcharge from dwat in aquifer layer
      ! update in drainage for case jwt < nlevsoi

      do fc = 1,num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         do j = 1, nlevbed
            h2osoi_liq(c,j) = h2osoi_liq(c,j) + dwat2(c,j)*dzmm(c,j)
         end do

         ! calculate qcharge for case jwt < nlevsoi
         if (use_var_soil_thick) then
            if (jwt(c) < nlevbed) then
               wh_zwt = 0._r8   !since wh_zwt = -sucsat - zq_zwt, where zq_zwt = -sucsat

               ! Recharge rate qcharge to groundwater (positive to aquifer)
               s_node = max(h2osoi_vol(c,jwt(c)+1)/watsat(c,jwt(c)+1), 0.01_r8)
               s1 = min(1._r8, s_node)

               !scs: this is the expression for unsaturated hk
               ka = imped(c,jwt(c)+1)*hksat(c,jwt(c)+1) &
                 *s1**(2._r8*bsw(c,jwt(c)+1)+3._r8)

               !compute unsaturated hk, this shall be tested later, because it
               !is not bit for bit
               !call soil_water_retention_curve%soil_hk(hksat(c,jwt(c)+1), s1, bsw(c,jwt(c)+1), ka)
               !apply ice impedance
               !ka = imped(c,jwt(c)+1) * ka 
               ! Recharge rate qcharge to groundwater (positive to aquifer)
               smp1 = max(smpmin(c), smp(c,max(1,jwt(c))))
               wh      = smp1 - zq(c,max(1,jwt(c)))

               !scs: original formulation
               if (jwt(c) == 0) then
                  qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)+1.e-3)*1000._r8)
               else
                  !             qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8)
                  !scs: 1/2, assuming flux is at zwt interface, saturation deeper than zwt
                  qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8*2.0)
               endif

               ! To limit qcharge  (for the first several timesteps)
               qcharge(c) = max(-10.0_r8/dtime,qcharge(c))
               qcharge(c) = min( 10.0_r8/dtime,qcharge(c))
            else
               ! if water table is below soil column, compute qcharge from dwat2(11)
               qcharge(c) = 0._r8
            endif
          else
            if (jwt(c) < nlevbed) then
               wh_zwt = 0._r8   !since wh_zwt = -sucsat - zq_zwt, where zq_zwt = -sucsat

               ! Recharge rate qcharge to groundwater (positive to aquifer)
               s_node = max(h2osoi_vol(c,jwt(c)+1)/watsat(c,jwt(c)+1), 0.01_r8)
               s1 = min(1._r8, s_node)

               !scs: this is the expression for unsaturated hk
               ka = imped(c,jwt(c)+1)*hksat(c,jwt(c)+1) &
                 *s1**(2._r8*bsw(c,jwt(c)+1)+3._r8)

               !compute unsaturated hk, this shall be tested later, because it
               !is not bit for bit
               !call soil_water_retention_curve%soil_hk(hksat(c,jwt(c)+1), s1, bsw(c,jwt(c)+1), ka)
               !apply ice impedance
               !ka = imped(c,jwt(c)+1) * ka 
               ! Recharge rate qcharge to groundwater (positive to aquifer)
               smp1 = max(smpmin(c), smp(c,max(1,jwt(c))))
               wh      = smp1 - zq(c,max(1,jwt(c)))

               !scs: original formulation
               if (jwt(c) == 0) then
                  qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)+1.e-3)*1000._r8)
               else
                  !             qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8)
                  !scs: 1/2, assuming flux is at zwt interface, saturation deeper than zwt
                  qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8*2.0)
               endif

               ! To limit qcharge  (for the first several timesteps)
               qcharge(c) = max(-10.0_r8/dtime,qcharge(c))
               qcharge(c) = min( 10.0_r8/dtime,qcharge(c))
            else
            ! if water table is below soil column, compute qcharge from dwat2(11)
               qcharge(c) = dwat2(c,nlevsoi+1)*dzmm(c,nlevsoi+1)/dtime
            endif
         endif
      end do

      ! compute the water deficit and reset negative liquid water content
      !  Jinyun Tang
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         qflx_deficit(c) = 0._r8
         do j = 1, nlevbed
            if(h2osoi_liq(c,j)<0._r8)then
               qflx_deficit(c) = qflx_deficit(c) - h2osoi_liq(c,j)
            endif
         enddo
      enddo

    end associate 

  end subroutine soilwater_zengdecker2009

  !-----------------------------------------------------------------------
  subroutine Prepare_Data_for_EM_VSFM_Driver(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, soilstate_vars, &
       waterflux_vars, waterstate_vars, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Prepare dataset that will be transfered from ALM to VSFM
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use decompMod                 , only : bounds_type
    use elm_varcon                , only : denh2o
    use elm_varpar                , only : nlevsoi, max_patch_per_col, nlevgrnd
    use clm_time_manager          , only : get_step_size
    use SoilStateType             , only : soilstate_type
    use SoilHydrologyType         , only : soilhydrology_type
    use TemperatureType           , only : temperature_type
    use WaterFluxType             , only : waterflux_type
    use WaterStateType            , only : waterstate_type
    use VegetationType            , only : veg_pp
    use ColumnType                , only : col_pp
    use elm_varcon                , only : watmin
    use LandunitType              , only : lun_pp
    use landunit_varcon           , only : istsoil, istcrop
    use elm_varctl                , only : lateral_connectivity
    use domainLateralMod          , only : ldomain_lateral
    !
    ! !ARGUMENTS:
    implicit none
    !
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    type(soilstate_type)    , intent(inout) :: soilstate_vars
    type(waterflux_type)    , intent(inout) :: waterflux_vars
    type(waterstate_type)   , intent(inout) :: waterstate_vars
    type(temperature_type)  , intent(in)    :: temperature_vars
    !
    ! !LOCAL VARIABLES:
    integer              :: p,c,fc,j,g                    ! do loop indices
    real(r8)             :: dtime                         ! land model time step (sec)
    real(r8)             :: temp(bounds%begc:bounds%endc) ! accumulator for rootr weighting
    integer              :: pi                            ! pft index
    real(r8)             :: dzsum                         ! summation of dzmm of layers below water table (mm)
    integer              :: jwt                           ! index of first unsaturated soil layer
    real(r8)             :: flux_unit_conversion          ! [mm/s] ---> [kg/s]
    real(r8)             :: area                          ! [m^2]
    real(r8)             :: z_up, z_dn                    ! [m]
    real(r8)             :: qflx_drain_layer              ! Drainage flux from a soil layer (mm H2O/s)
    real(r8)             :: qflx_drain_tot                ! Cummulative drainage flux from soil layers within a column (mm H2O/s)

    !-----------------------------------------------------------------------

    associate( &
         zi                        =>    col_pp%zi                                     , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         dz                        =>    col_pp%dz                                     , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         snl                       =>    col_pp%snl                                    , & ! Input:  [integer  (:)   ]  minus number of snow layers

         zwt                       =>    soilhydrology_vars%zwt_col                 , & ! Input:  [real(r8) (:)   ]  water table depth (m)

         watsat                    =>    soilstate_vars%watsat_col                  , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         rootr_col                 =>    soilstate_vars%rootr_col                   , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         rootr_pft                 =>    soilstate_vars%rootr_patch                 , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer

         h2osoi_liq                =>    col_ws%h2osoi_liq             , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         frac_h2osfc               =>    col_ws%frac_h2osfc            , & ! Input:  [real(r8) (:)   ]
         frac_sno                  =>    col_ws%frac_sno_eff           , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         vsfm_fliq_col_1d          =>    col_ws%vsfm_fliq_col_1d       , & ! Output: [real(r8) (:)   ]  1D fraction of liquid saturation for VSFM [-]
         vsfm_sat_col_1d           =>    col_ws%vsfm_sat_col_1d        , & ! Output: [real(r8) (:)   ]  1D liquid saturation from VSFM [-]
         vsfm_mass_col_1d          =>    col_ws%vsfm_mass_col_1d       , & ! Output: [real(r8) (:)   ]  1D liquid mass per unit area from VSFM [kg H2O/m^2]
         vsfm_smpl_col_1d          =>    col_ws%vsfm_smpl_col_1d       , & ! Output: [real(r8) (:)   ]  1D soil matrix potential liquid from VSFM [m]
         vsfm_soilp_col_1d         =>    col_ws%vsfm_soilp_col_1d      , & ! Output: [real(r8) (:)   ]  1D soil water pressure from VSFM [Pa]
         soilp_col                 =>    col_ws%soilp                  , & ! Output: [real(r8) (:,:) ]  soil water pressure (Pa)
         qflx_rootsoi_col          =>    col_wf%qflx_rootsoi            , & ! Input:  [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
         qflx_deficit              =>    col_wf%qflx_deficit            , & ! Input:  [real(r8) (:)   ]  water deficit to keep non-negative liquid water content
         qflx_infl                 =>    col_wf%qflx_infl               , & ! Input:  [real(r8) (:)   ]  infiltration (mm H2O /s)
         qflx_dew_snow             =>    col_wf%qflx_dew_snow           , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
         qflx_dew_grnd             =>    col_wf%qflx_dew_grnd           , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_sub_snow             =>    col_wf%qflx_sub_snow           , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_drain                =>    col_wf%qflx_drain              , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)

         mflx_infl_col             =>    col_wf%mflx_infl               , & ! Output: [real(r8) (:)   ]  infiltration source in top soil control volume (kg H2O /s)
         mflx_dew_col              =>    col_wf%mflx_dew                , & ! Output: [real(r8) (:)   ]  (liquid+snow) dew source in top soil control volume (kg H2O /s)
         mflx_snowlyr_disp_col     =>    col_wf%mflx_snowlyr_disp       , & ! Output: [real(r8) (:)   ]  mass flux to top soil layer due to disappearance of snow (kg H2O /s)
         mflx_sub_snow_col         =>    col_wf%mflx_sub_snow           , & ! Output: [real(r8) (:)   ]  mass flux from top soil layer due to sublimation of snow (kg H2O /s)
         mflx_et_col               =>    col_wf%mflx_et                 , & ! Output: [real(r8) (:)   ]  evapotranspiration sink from all soil coontrol volumes (kg H2O /s) (+ = to atm)
         mflx_drain_col            =>    col_wf%mflx_drain              , & ! Output: [real(r8) (:)   ]  drainage from groundwater and perched water table (kg H2O /s)
         mflx_snowlyr_col          =>    col_wf%mflx_snowlyr            , & ! Output: [real(r8) (:)   ]  mass flux to top soil layer due to disappearance of snow (kg H2O /s)
         mflx_neg_snow_col_1d      =>    col_wf%mflx_neg_snow_1d        , & ! Input:  [real(r8) (:)   ]  mass flux from top soil layer due to negative water content in snow layers (kg H2O /s)

         t_soisno                  =>    col_es%t_soisno                & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         )

      ! Get time step

      dtime = get_step_size()

      mflx_infl_col(:)              = 0.d0
      mflx_dew_col(:)               = 0.d0
      mflx_snowlyr_disp_col(:)      = 0.d0
      mflx_sub_snow_col(:)          = 0.d0
      mflx_et_col(:,:)              = 0.d0
      mflx_drain_col(:,:)           = 0.d0

      area = 1.d0 ! [m^2]

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

#ifdef USE_PETSC_LIB
         if (lateral_connectivity) then
            g    = col_pp%gridCell(c)
            area = ldomain_lateral%ugrid%areaGrid_ghosted(g)
         endif
#endif

         ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
         flux_unit_conversion     = area * denh2o * 1.0d-3

         do j = 1, nlevsoi
            ! ET sink
            mflx_et_col(c,j) = -qflx_rootsoi_col(c,j)*flux_unit_conversion
         end do

         ! Infiltration source term
         mflx_infl_col(c) = qflx_infl(c)*flux_unit_conversion

         ! Dew and snow sublimation source/sink term
         if (snl(c) >= 0) then
            mflx_dew_col(c)       = (qflx_dew_snow(c) + qflx_dew_grnd(c))* &
                                    (1._r8 - frac_h2osfc(c))*              &
                                     flux_unit_conversion

            mflx_sub_snow_col(c)  = -qflx_sub_snow(c)*          &
                                     (1._r8 - frac_h2osfc(c))*  &
                                      flux_unit_conversion
         end if

         if (qflx_drain(c) > 0.d0) then

            ! Find soil layer just above water table
            jwt = nlevgrnd

            ! allow jwt to equal zero when zwt is in top layer
            do j = 1,nlevgrnd
               if (zwt(c) <= zi(c,j)) then
                  jwt = j-1
                  exit
               end if
            enddo

            ! Now ensure the soil layer index corresponding to the water table
            ! is greater than or equal to the first soil layer.
            jwt = max(jwt,1)

            dzsum = 0.d0
            do j = jwt, nlevgrnd
               dzsum = dzsum + dz(c,j)
            end do

            qflx_drain_tot = 0.d0
            do j = jwt, nlevgrnd
               qflx_drain_layer = qflx_drain(c) * dz(c,j)/dzsum

               ! if the amount of water being drained from a given layer
               ! exceeds the allowable water, limit the drainage
               if (qflx_drain_layer*dtime > (h2osoi_liq(c,j)-watmin)) then
                  qflx_drain_layer = (h2osoi_liq(c,j)-watmin)/dtime
               endif
               qflx_drain_tot = qflx_drain_tot + qflx_drain_layer

               mflx_drain_col(c,j) = -qflx_drain_layer*flux_unit_conversion

           end do
           qflx_drain(c) = qflx_drain_tot

         endif

         ! The mass flux associated with disapperance of snow layer over the
         ! last time step.
         mflx_snowlyr_disp_col(c) = mflx_snowlyr_col(c)*area + &
                                    mflx_neg_snow_col_1d(c-bounds%begc+1)*area
         mflx_snowlyr_col(c) = 0._r8

      end do

      ! compute the water deficit and reset negative liquid water content
      !  Jinyun Tang
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qflx_deficit(c) = 0._r8
      enddo

    end associate

  end subroutine Prepare_Data_for_EM_VSFM_Driver

   ! ====================================================================================
   
   subroutine Compute_EffecRootFrac_And_VertTranSink(bounds, num_hydrologyc, &
         filter_hydrologyc, soilstate_inst, canopystate_inst, waterflux_inst, energyflux_inst)
      
      ! ---------------------------------------------------------------------------------
      ! This is a wrapper for calculating the effective root fraction and soil
      ! water sink due to plant transpiration. 
      ! Calculate Soil Water Sink to Roots over different types
      ! of columns and for different process modules
      ! The super-set of all columns that should have a root water sink
      ! is filter_hydrologyc
      ! There are three two of columns:
      ! 1) impervious roads, 2) non-natural vegetation and  natural vegetation
      ! There are several methods available.
      ! 1) the default version, 2) hydstress version and 3) fates boundary conditions
      !
      ! There are only two quantities that are the result of this routine, and its
      ! children:
      !   waterflux_inst%qflx_rootsoi_col(c,j)
      !   soilstate_inst%rootr_col(c,j)
      !
      !
      ! ---------------------------------------------------------------------------------

      use SoilStateType       , only : soilstate_type
      use WaterFluxType       , only : waterflux_type
      use CanopyStateType     , only : canopystate_type
      use EnergyFluxType      , only : energyflux_type
      use ColumnType          , only : col_pp 
      use LandunitType        , only : lun_pp
      use decompMod           , only : bounds_type   
      use column_varcon       , only : icol_road_perv
      use elm_varctl          , only : use_hydrstress, iulog
      use shr_log_mod         , only : errMsg => shr_log_errMsg
      use abortutils          , only : endrun

      ! Arguments
      type(bounds_type)       , intent(in)    :: bounds               ! bounds
      integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
      integer                 , intent(in)    :: filter_hydrologyc(num_hydrologyc) ! column filter for soil points
      type(soilstate_type)    , intent(inout) :: soilstate_inst
      type(waterflux_type)    , intent(inout) :: waterflux_inst
      type(canopystate_type)  , intent(in)    :: canopystate_inst
      type(energyflux_type)   , intent(in)    :: energyflux_inst

      ! Local Variables
      integer  :: filterc(bounds%endc-bounds%begc+1)           !column filter
      integer  :: num_filterc
      integer  :: num_filterc_tot
      integer  :: fc
      integer  :: c
      integer  :: l

      num_filterc_tot = 0

      ! 1) pervious roads
      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (col_pp%itype(c) == icol_road_perv) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do
      num_filterc_tot = num_filterc_tot+num_filterc
      call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
               num_filterc,filterc, soilstate_inst, waterflux_inst)


      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         l = col_pp%landunit(c)
         if ( (col_pp%itype(c) /= icol_road_perv) ) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do
      num_filterc_tot = num_filterc_tot+num_filterc
      if(use_hydrstress) then
         call Compute_EffecRootFrac_And_VertTranSink_HydStress(bounds, &
               num_filterc, filterc, waterflux_inst, soilstate_inst, &
               canopystate_inst, energyflux_inst)
      else
         call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
               num_filterc,filterc, soilstate_inst, waterflux_inst)
      end if

      if (num_hydrologyc /= num_filterc_tot) then
          write(iulog,*) 'The total number of columns flagged to root water uptake'
          write(iulog,*) 'did not match the total number calculated'
          write(iulog,*) 'This is likely a problem with the interpretation of column/lu filters.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
      end if


      return
   end subroutine Compute_EffecRootFrac_And_VertTranSink

   subroutine Compute_EffecRootFrac_And_VertTranSink_Default(bounds, num_filterc, &
         filterc, soilstate_vars, waterflux_vars)

    !
    ! Generic routine to apply transpiration as a sink condition that
    ! is vertically distributed over the soil column. Should be
    ! applicable to any Richards solver that is not coupled to plant
    ! hydraulics.
    !
    !USES:
    use decompMod        , only : bounds_type
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use elm_varpar       , only : nlevsoi, max_patch_per_col
    use SoilStateType    , only : soilstate_type
    use WaterFluxType    , only : waterflux_type
    use VegetationType   , only : veg_pp
    use ColumnType       , only : col_pp
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds                          ! bounds
    integer              , intent(in)    :: num_filterc                     ! number of column soil points in column filter
    integer              , intent(in)    :: filterc(num_filterc)            ! column filter for soil points
    type(waterflux_type) , intent(inout) :: waterflux_vars
    type(soilstate_type) , intent(inout) :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,fc,j                                              ! do loop indices
    integer  :: pi                                                    ! patch index
    integer  :: nlevbed                                               ! number of layers to bedrock
    real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting
    associate(& 
          nlev2bed            =>    col_pp%nlevbed                     , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     
          qflx_rootsoi_col    => col_wf%qflx_rootsoi    , & ! Output: [real(r8) (:,:) ]  
                                                                        ! vegetation/soil water exchange (m H2O/s) (+ = to atm)
          qflx_tran_veg_patch => veg_wf%qflx_tran_veg , & ! Input:  [real(r8) (:)   ]  
                                                                        ! vegetation transpiration (mm H2O/s) (+ = to atm) 
          qflx_tran_veg_col   => col_wf%qflx_tran_veg   , & ! Input:  [real(r8) (:)   ]  
                                                                        ! vegetation transpiration (mm H2O/s) (+ = to atm)
          qflx_rootsoi_frac_patch  =>    veg_wf%qflx_rootsoi_frac    , & ! Output: [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/s) (+ = to atm)
          rootr_patch         => soilstate_vars%rootr_patch         , & ! Input: [real(r8) (:,:) ]
                                                                        ! effective fraction of roots in each soil layer  
          rootr_col           => soilstate_vars%rootr_col             & ! Output: [real(r8) (:,:) ]  
                                                                        ! effective fraction of roots in each soil layer  
          )
      
      ! First step is to calculate the column-level effective rooting
      ! fraction in each soil layer. This is done outside the usual
      ! PATCH-to-column averaging routines because it is not a simple
      ! weighted average of the PATCH level rootr arrays. Instead, the
      ! weighting depends on both the per-unit-area transpiration
      ! of the PATCH and the PATCHEs area relative to all PATCHES.
      
      temp(bounds%begc : bounds%endc) = 0._r8
      
      do fc = 1, num_filterc
         c = filterc(fc)
         nlevbed = nlev2bed(c)
         do j = 1, nlevbed
            rootr_col(c,j) = 0._r8
         end do
      end do
      
      do pi = 1,max_patch_per_col
         do fc = 1, num_filterc
            c = filterc(fc)
            nlevbed = nlev2bed(c)
            do j = 1,nlevbed
               if (pi <= col_pp%npfts(c)) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then
                     rootr_col(c,j) = rootr_col(c,j) + rootr_patch(p,j) * &
                          qflx_tran_veg_patch(p) * veg_pp%wtcol(p)
                     qflx_rootsoi_frac_patch(p,j) = rootr_patch(p,j) * qflx_tran_veg_patch(p) * veg_pp%wtcol(p)
                  end if
               end if
            end do
         end do
         do fc = 1, num_filterc
            c = filterc(fc)
            if (pi <= col_pp%npfts(c)) then
               p = col_pp%pfti(c) + pi - 1
               if (veg_pp%active(p)) then
                  temp(c) = temp(c) + qflx_tran_veg_patch(p) * veg_pp%wtcol(p)
               end if
            end if
         end do
      end do
      
      do fc = 1, num_filterc
         c = filterc(fc)
         nlevbed = nlev2bed(c)
         do j = 1, nlevbed
            if (temp(c) /= 0._r8) then
               rootr_col(c,j) = rootr_col(c,j)/temp(c)
            end if
            qflx_rootsoi_col(c,j) = rootr_col(c,j)*qflx_tran_veg_col(c)

         end do
      end do

      do pi = 1,max_patch_per_col
         do j = 1,nlevsoi
            do fc = 1, num_filterc
               c = filterc(fc)
               if (pi <= col_pp%npfts(c)) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then
                    if(rootr_col(c,j)==0._r8)then
                      qflx_rootsoi_frac_patch(p,j) = 0._r8
                    else
                      qflx_rootsoi_frac_patch(p,j) = qflx_rootsoi_frac_patch(p,j)/(temp(c)*rootr_col(c,j))
                    endif
                  end if
               end if
            end do
         end do
       enddo
    end associate
    return
 end subroutine Compute_EffecRootFrac_And_VertTranSink_Default

   ! ==================================================================================
   
   subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress( bounds, &
           num_filterc, filterc, waterflux_vars, soilstate_vars, &
           canopystate_vars, energyflux_vars)


        !
        !USES:
        use decompMod        , only : bounds_type
        use elm_varpar       , only : nlevsoi
        use elm_varpar       , only : max_patch_per_col
        use SoilStateType    , only : soilstate_type
        use WaterFluxType    , only : waterflux_type
        use CanopyStateType  , only : canopystate_type
        use VegetationType   , only : veg_pp
        use ColumnType       , only : col_pp
        use elm_varctl       , only : iulog
        use PhotosynthesisMod, only : plc, params_inst
        use column_varcon    , only : icol_road_perv
        use shr_infnan_mod   , only : isnan => shr_infnan_isnan
        use EnergyFluxType   , only : energyflux_type
        use shr_kind_mod     , only : r8 => shr_kind_r8
        !
        ! !ARGUMENTS:
        type(bounds_type)    , intent(in)    :: bounds          ! bounds
        integer              , intent(in)    :: num_filterc     ! number of column soil points in column filter
        integer              , intent(in)    :: filterc(:)      ! column filter for soil points
        type(waterflux_type) , intent(inout) :: waterflux_vars
        type(soilstate_type) , intent(inout) :: soilstate_vars
        type(canopystate_type) , intent(in)  :: canopystate_vars
        type(energyflux_type), intent(in)    :: energyflux_vars
        !
        ! !LOCAL VARIABLES:
        integer  :: p,c,fc,j                                              ! do loop indices
        integer  :: pi                                                    ! patch index
        real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting
        real(r8) :: grav2                 ! soil layer gravitational potential relative to surface (mm H2O)
        integer , parameter :: soil=1,root=4  ! index values
        !-----------------------------------------------------------------------   
        
        associate(&
              k_soil_root         => soilstate_vars%k_soil_root_patch   , & ! Input:  [real(r8) (:,:) ]  
                                                                            ! soil-root interface conductance (mm/s)
              qflx_phs_neg_col    => waterflux_vars%qflx_phs_neg_col    , & ! Input:  [real(r8) (:)   ]  n
                                                                            ! net neg hydraulic redistribution flux(mm H2O/s)
              qflx_tran_veg_col   => col_wf%qflx_tran_veg   , & ! Input:  [real(r8) (:)   ]  
                                                                            ! vegetation transpiration (mm H2O/s) (+ = to atm)
              qflx_tran_veg_patch => veg_wf%qflx_tran_veg , & ! Input:  [real(r8) (:)   ]  
                                                                            ! vegetation transpiration (mm H2O/s) (+ = to atm)
              qflx_rootsoi_col    => col_wf%qflx_rootsoi    , & ! Output: [real(r8) (:)   ]
                                                                            ! col root and soil water 
                                                                            ! exchange [mm H2O/s] [+ into root]
              rootr_col           => soilstate_vars%rootr_col           , & ! Input:  [real(r8) (:,:) ]
                                                                            ! effective fraction of roots in each soil layer
              rootr_patch         => soilstate_vars%rootr_patch         , & ! Input:  [real(r8) (:,:) ]  
                                                                            ! effective fraction of roots in each soil layer
              smp                 => soilstate_vars%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix pot. [mm]
              frac_veg_nosno      => canopystate_vars%frac_veg_nosno_patch , & ! Input:  [integer  (:)  ] 
                                                                            ! fraction of vegetation not 
                                                                            ! covered by snow (0 OR 1) [-]  
              z                   => col_pp%z                              , & ! Input: [real(r8) (:,:) ]  layer node depth (m)
              vegwp               => canopystate_vars%vegwp_patch         & ! Input: [real(r8) (:,:) ]  vegetation water 
                                                                            ! matric potential (mm)
              )
          
          do fc = 1, num_filterc
             c = filterc(fc)
             qflx_phs_neg_col(c) = 0._r8
             
             do j = 1, nlevsoi
                grav2 = z(c,j) * 1000._r8
                temp(c) = 0._r8
                do pi = 1,max_patch_per_col
                   if (pi <= col_pp%npfts(c)) then
                      p = col_pp%pfti(c) + pi - 1
                      if (veg_pp%active(p).and.frac_veg_nosno(p)>0) then 
                         if (veg_pp%wtcol(p) > 0._r8) then
                            temp(c) = temp(c) + k_soil_root(p,j) &
                                  * (smp(c,j) - vegwp(p,4) - grav2)* veg_pp%wtcol(p)
                         endif
                      end if
                   end if
                end do
                qflx_rootsoi_col(c,j)= temp(c)
                
                if (temp(c) < 0._r8) qflx_phs_neg_col(c) = qflx_phs_neg_col(c) + temp(c)
             end do
             
             ! Back out the effective root density
             if( sum(qflx_rootsoi_col(c,1:nlevsoi))>0.0_r8 ) then
                do j = 1, nlevsoi
                   rootr_col(c,j) = qflx_rootsoi_col(c,j)/sum( qflx_rootsoi_col(c,1:nlevsoi))
                end do
             else
                rootr_col(c,:) = 0.0_r8
             end if
          end do
          
        end associate

        return
     end subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress


end module SoilWaterMovementMod
