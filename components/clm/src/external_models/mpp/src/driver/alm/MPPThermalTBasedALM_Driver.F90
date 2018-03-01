
module MPPThermalTBasedALM_Driver

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use decompMod        , only : bounds_type, get_proc_bounds 
  use abortutils       , only : endrun
  use mpp_varctl       , only : iulog
  ! 
  !-----------------------------------------
  ! Definition of component types
  !-----------------------------------------
  use LandunitType           , only : lun                
  use ColumnType             , only : col                

  implicit none
  
  public :: MPPThermalTBasedALM_Solve

  !------------------------------------------------------------------------
contains

  subroutine MPPThermalTBasedALM_Solve(bounds, num_filter, filter, &
       dtime,  sabg_lyr_col, dhsdT, hs_soil, hs_top_snow, hs_h2osfc, &
       waterstate_vars,  temperature_vars, tvector)

    !
    ! !DESCRIPTION:
    !
#include <petsc/finclude/petsc.h>
    ! !USES:
    use WaterstateType            , only : waterstate_type
    use TemperatureType           , only : temperature_type
    use clm_time_manager          , only : get_step_size, get_nstep
    use mpp_varpar                , only : nlevsno, nlevgrnd
    use mpp_varcon                , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
    use mpp_varcon                , only : istice, istice_mec, istsoil, istcrop
    use mpp_varcon                , only : istwet, istice, istice_mec, istsoil, istcrop
    use mpp_varcon                , only : capr
    use MultiPhysicsProbThermal   , only : thermal_mpp
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_AREAL_DEN
    use MultiPhysicsProbConstants , only : VAR_ICE_AREAL_DEN
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_SNOW_WATER
    use MultiPhysicsProbConstants , only : VAR_NUM_SNOW_LYR
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    use MultiPhysicsProbConstants , only : VAR_TUNING_FACTOR
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_DHS_DT
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use petscsys
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)  :: bounds
    integer                , intent(in)  :: num_filter                                         ! number of columns in the filter
    integer                , intent(in)  :: filter(:)                                          ! column filter
    real(r8)               , intent(in)  :: dtime                                              ! land model time step (sec)
    real(r8)               , intent(in)  :: sabg_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:1) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8)               , intent(in)  :: dhsdT(bounds%begc:bounds%endc)                     ! temperature derivative of "hs" [col]
    real(r8)               , intent(in)  :: hs_soil(bounds%begc:bounds%endc)                   ! heat flux on soil [W/m2]
    real(r8)               , intent(in)  :: hs_top_snow(bounds%begc:bounds%endc)               ! heat flux on top snow layer [W/m2]
    real(r8)               , intent(in)  :: hs_h2osfc(bounds%begc:bounds%endc)                 ! heat flux on standing water [W/m2]
    type(waterstate_type)  , intent(in)  :: waterstate_vars
    type(temperature_type) , intent(in)  :: temperature_vars
    real(r8)               , intent(out) :: tvector( bounds%begc: , -nlevsno: )                ! Numerical solution of temperature
    !
    ! !LOCAL VARIABLES:
    integer                              :: j,c,l,idx                                                                  !  indices
    integer                              :: offset
    integer                              :: fc                                                                         ! lake filtered column indices

    integer                              :: soe_auxvar_id

    ! internal auxvars
    real(r8) , pointer                   :: temperature_1d         (:)
    real(r8) , pointer                   :: liq_areal_den_1d       (:)
    real(r8) , pointer                   :: ice_areal_den_1d       (:)
    real(r8) , pointer                   :: snow_water_1d          (:)
    real(r8) , pointer                   :: dz_1d                  (:)
    real(r8) , pointer                   :: dist_up_1d             (:)
    real(r8) , pointer                   :: dist_dn_1d             (:)
    real(r8) , pointer                   :: frac_1d                (:)
    integer  , pointer                   :: num_snow_layer_1d      (:)
    logical  , pointer                   :: is_active_1d           (:)

    ! boundary auxvars
    real(r8) , pointer                   :: hs_snow_1d             (:)
    real(r8) , pointer                   :: hs_sh2o_1d             (:)
    real(r8) , pointer                   :: hs_soil_1d             (:)
    real(r8) , pointer                   :: dhsdT_snow_1d          (:)
    real(r8) , pointer                   :: dhsdT_sh2o_1d          (:)
    real(r8) , pointer                   :: dhsdT_soil_1d          (:)
    real(r8) , pointer                   :: frac_soil_1d           (:)

    ! source sink
    real(r8) , pointer                   :: sabg_snow_1d           (:)
    real(r8) , pointer                   :: sabg_soil_1d           (:)

    real(r8) , pointer                   :: tsurf_tuning_factor_1d (:)


    PetscErrorCode                       :: ierr                                                                       ! PETSc return error code
    PetscBool                            :: converged                                                                  ! Did thermal solver converge to a solution?
    PetscInt                             :: converged_reason
    !-----------------------------------------------------------------------
    associate(                                                                & 
         snl                     => col%snl                                 , & ! Input:  [integer  (:)   ]  number of snow layers                    
         zi                      => col%zi                                  , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m) 
         dz                      => col%dz                                  , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                       
         z                       => col%z                                   , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                   
         
         frac_sno_eff            => waterstate_vars%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)
         h2osno                  => waterstate_vars%h2osno_col              , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
         h2osfc                  => waterstate_vars%h2osfc_col              , & ! Input:  [real(r8) (:)   ]  surface water (mm)                      
         frac_h2osfc             => waterstate_vars%frac_h2osfc_col         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         h2osoi_liq              => waterstate_vars%h2osoi_liq_col          , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                  
         h2osoi_ice              => waterstate_vars%h2osoi_ice_col          , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                      
         
         t_h2osfc                => temperature_vars%t_h2osfc_col           , & ! Output: [real(r8) (:)   ]  surface water temperature               
         t_soisno                => temperature_vars%t_soisno_col           , & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)             
         
         begc                    => bounds%begc                             , &
         endc                    => bounds%endc                               &
         )

      ! Allocate memory
      allocate(temperature_1d         ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(liq_areal_den_1d       ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(ice_areal_den_1d       ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(snow_water_1d          ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(dz_1d                  ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(dist_up_1d             ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(dist_dn_1d             ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(frac_1d                ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))

      allocate(num_snow_layer_1d      ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      allocate(is_active_1d           ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))

      allocate(hs_snow_1d             ((bounds%endc-bounds%begc+1                      )))
      allocate(hs_sh2o_1d             ((bounds%endc-bounds%begc+1                      )))
      allocate(hs_soil_1d             ((bounds%endc-bounds%begc+1                      )))
      allocate(dhsdT_snow_1d          ((bounds%endc-bounds%begc+1                      )))
      allocate(dhsdT_sh2o_1d          ((bounds%endc-bounds%begc+1                      )))
      allocate(dhsdT_soil_1d          ((bounds%endc-bounds%begc+1                      )))
      allocate(frac_soil_1d           ((bounds%endc-bounds%begc+1                      )))

      allocate(sabg_snow_1d           ((bounds%endc-bounds%begc+1)*nlevsno ))
      allocate(sabg_soil_1d           ((bounds%endc-bounds%begc+1)*nlevgrnd))

      allocate(tsurf_tuning_factor_1d ((bounds%endc-bounds%begc+1)*(nlevgrnd+nlevsno+1 )))
      ! Initialize
      temperature_1d(:)         = 273.15_r8 ! temperature_vars%t_soisno_col + temperature_vars%t_h2osfc_col
      liq_areal_den_1d(:)       = 0._r8     ! waterstate_vars%h2osoi_liq_col
      ice_areal_den_1d(:)       = 0._r8     ! waterstate_vars%h2osoi_ice_col
      frac_soil_1d(:)           = 1._r8     ! 1.0 - waterstate_vars%frac_sno_eff_col - waterstate_vars%frac_h2osfc_col
      frac_1d(:)                = 1._r8     !
      num_snow_layer_1d(:)      = 0         ! col%snl
      is_active_1d(:)           = .false.   ! col%snl
      hs_snow_1d(:)             = 0._r8     ! LOCAL VARIABLE: hs_top_snow
      hs_sh2o_1d(:)             = 0._r8     ! LOCAL VARIABLE: hs_h2osfc
      hs_soil_1d(:)             = 0._r8     ! LOCAL VARIABLE: hs_soil
      dhsdT_snow_1d(:)          = 0._r8     ! LOCAL VARIABLE: hs_top_snow
      dhsdT_sh2o_1d(:)          = 0._r8     ! LOCAL VARIABLE: hs_h2osfc
      dhsdT_soil_1d(:)          = 0._r8     ! LOCAL VARIABLE: hs_soil
      sabg_snow_1d(:)           = 0._r8     ! LOCAL VARIABLE: sabg_lyr_col
      sabg_soil_1d(:)           = 0._r8     ! LOCAL VARIABLE: sabg_lyr_col
      snow_water_1d(:)          = 0._r8
      tsurf_tuning_factor_1d(:) = 1._r8
      dz_1d(:)                  = 0._r8
      dist_up_1d(:)             = 0._r8
      dist_dn_1d(:)             = 0._r8

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Save data for snow
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      offset = 0

      do fc = 1, num_filter
         c = filter(fc)
         do j = -nlevsno+1, 0

            l = col%landunit(c)

            ! Is this a soil column on which PETSc based thermal solver works?
            if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
            
               if (j >= snl(c)+1) then
                  
                  ! Index for internal SoE auxvars
                  idx = (c-bounds%begc)*nlevsno + j + nlevsno + offset

                  ! Save data for internal SoE auxvars
                  temperature_1d(idx)    = t_soisno(c,j)
                  dz_1d(idx)             = dz(c,j)
                  liq_areal_den_1d(idx)  = h2osoi_liq(c,j)
                  ice_areal_den_1d(idx)  = h2osoi_ice(c,j)
                  num_snow_layer_1d(idx) = -snl(c)
                  is_active_1d(idx)      = .true.
                  dist_up_1d(idx)        = col%zi(c,j  ) - col%z(c,j)
                  dist_dn_1d(idx)        = col%z(c,j) - col%zi(c,j-1)
                  frac_1d(idx)           = frac_sno_eff(c)

                  ! If not the top snow layer, save amount of absorbed solar
                  ! radiation
                  if (j /= snl(c) +  1) then
                     sabg_snow_1d(idx)      = sabg_lyr_col(c,j)
                  endif

                  ! Save follow data only for the top snow layer
                  if (j == snl(c)+1) then

                     ! Save tuning_factor for internal SoE auxvars
                     tsurf_tuning_factor_1d(idx) = dz(c,j) / &
                        (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))

                     ! Index for boundary SoE auxvars
                     idx = (c-bounds%begc)+1

                     ! Save data for boundary SoE auxvars
                     hs_snow_1d(idx)    = hs_top_snow(c)
                     dhsdT_snow_1d(idx) = dhsdT(c)
                     frac_soil_1d(idx)  = frac_soil_1d(idx) - frac_sno_eff(c)
                  endif

               endif
            endif
         enddo
      enddo

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Save data for h2osfc
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      offset = (bounds%endc - bounds%begc + 1)*nlevsno ! Number of data for snow

      do fc = 1, num_filter
         c = filter(fc)
         l = col%landunit(c)

         if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
            
            if (frac_h2osfc(c) > 0._r8) then
                  
               ! Index for internal SoE auxvars
               idx = (c-bounds%begc) + 1 + offset

               ! Save data for internal SoE auxvars
               temperature_1d(idx) = t_h2osfc(c)
               dz_1d(idx)          = 1.0e-3*h2osfc(c)
               is_active_1d(idx)   = .true.
               frac_1d(idx)        = frac_h2osfc(c)
               dist_up_1d(idx)     = dz_1d(idx)/2.d0
               dist_dn_1d(idx)     = dz_1d(idx)/2.d0

               ! Index for boundary SoE auxvars
               idx                 = (c-bounds%begc) + 1

               ! Save data for boundary SoE auxvars
               frac_soil_1d(idx)   = frac_soil_1d(idx) - frac_h2osfc(c)
               dhsdT_sh2o_1d(idx)  = dhsdT(c)               
               hs_sh2o_1d(idx)     = hs_h2osfc(c)
            endif
         endif
      enddo

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Save data for soil
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      offset = (bounds%endc - bounds%begc + 1)*(nlevsno + 1) ! Number of data for snow + sh2o

      do fc = 1, num_filter
         c = filter(fc)
         l = col%landunit(c)

         do j = 1,nlevgrnd

            ! Index for internal SoE auxvars
            idx = (c-bounds%begc)*nlevgrnd + j + offset

            if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then            
                  
               ! Save data for internal SoE auxvars
               temperature_1d(idx)    = t_soisno(c,j)
               dz_1d(idx)             = dz(c,j)
               is_active_1d(idx)      = .true.
               liq_areal_den_1d(idx)  = h2osoi_liq(c,j)
               ice_areal_den_1d(idx)  = h2osoi_ice(c,j)
               frac_1d(idx)           = 1.0_r8
               dist_up_1d(idx)        = col%zi(c,j)   - col%z(c,j)
               dist_dn_1d(idx)        = col%zi(c,j)   - col%z(c,j)

               ! Is this the top soil layer?
               if (j == 1) then

                  dz_1d(idx)             = z(c,j)*2.d0
                  num_snow_layer_1d(idx) = -snl(c)

                  if (snl(c) /= 0) then
                     ! Save data for boundary SoE auxvars
                     sabg_soil_1d((c-bounds%begc)*nlevgrnd + j) = frac_sno_eff(c)*sabg_lyr_col(c,j)
                     ! Save data for internal SoE auxvars
                     snow_water_1d(idx)                         = h2osno(c)
                  else
                     ! Save data for internal SoE auxvars
                     tsurf_tuning_factor_1d(idx) = dz(c,j) / &
                        (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
                  endif

                  ! Save data for boundary SoE auxvars
                  hs_soil_1d(   (c-bounds%begc) + 1) = hs_soil(c)
                  dhsdT_soil_1d((c-bounds%begc) + 1) = dhsdT(c)

               endif
            endif
         enddo
      enddo

      ! Set temperature
      call thermal_mpp%sysofeqns%SetSolnPrevCLM(temperature_1d)

      ! Set h2soi_liq
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_LIQ_AREAL_DEN, soe_auxvar_id, liq_areal_den_1d)

      ! Set h2osi_ice
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_ICE_AREAL_DEN, soe_auxvar_id, ice_areal_den_1d)

      ! Set snow water
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_SNOW_WATER, soe_auxvar_id, snow_water_1d)

      ! Set dz
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_DZ, soe_auxvar_id, dz_1d)

      ! Set dist_up
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_DIST_UP, soe_auxvar_id, dist_up_1d)

      ! Set dist_dn
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_DIST_DN, soe_auxvar_id, dist_dn_1d)

      ! Set number of snow layers
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetIDataFromCLM(AUXVAR_INTERNAL, &
           VAR_NUM_SNOW_LYR, soe_auxvar_id, num_snow_layer_1d)

      ! Set if cell is active
      call thermal_mpp%sysofeqns%SetBDataFromCLM(AUXVAR_INTERNAL, &
           VAR_ACTIVE, is_active_1d)

      ! Set tuning factor
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_TUNING_FACTOR, soe_auxvar_id, tsurf_tuning_factor_1d)

      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_INTERNAL, &
           VAR_FRAC, soe_auxvar_id, frac_1d)

      !
      ! Set heat flux for:
      !

      ! 1) top snow layer
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_BC, &
           VAR_BC_SS_CONDITION, soe_auxvar_id, hs_snow_1d)

      ! 2) top standing water layer
      soe_auxvar_id = 2;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_BC, &
           VAR_BC_SS_CONDITION, soe_auxvar_id, hs_sh2o_1d)
      
      ! 3) soil
      soe_auxvar_id = 3;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_BC, &
           VAR_BC_SS_CONDITION, soe_auxvar_id, hs_soil_1d)

      !
      ! Set derivative of heat flux w.r.t temperature for:
      !

      ! 1) top snow layer
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_BC, &
           VAR_DHS_DT, soe_auxvar_id, dhsdT_snow_1d)

      ! 2) top standing water layer
      soe_auxvar_id = 2;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_BC, &
           VAR_DHS_DT, soe_auxvar_id, dhsdT_sh2o_1d)
      
      ! 3) soil
      soe_auxvar_id = 3;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_BC, &
           VAR_DHS_DT, soe_auxvar_id, dhsdT_soil_1d)

      !
      ! Set fraction of soil not covered by snow and standing water
      !
      soe_auxvar_id = 3;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_BC, &
           VAR_FRAC, soe_auxvar_id, frac_soil_1d)

      
      ! Set absorbed solar radiation
      soe_auxvar_id = 1;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_SS, &
           VAR_BC_SS_CONDITION, soe_auxvar_id, sabg_snow_1d)

      ! Set absorbed solar radiation
      soe_auxvar_id = 2;
      call thermal_mpp%sysofeqns%SetRDataFromCLM(AUXVAR_SS, &
           VAR_BC_SS_CONDITION, soe_auxvar_id, sabg_soil_1d)
      
      
      ! Preform Pre-StepDT operations
      call thermal_mpp%sysofeqns%PreStepDT()

      ! Solve
      call thermal_mpp%sysofeqns%StepDT(dtime, get_nstep(), &
           converged, converged_reason, ierr); CHKERRQ(ierr)

      ! Did the model converge
      if (.not. converged) then
         call endrun(msg=' ERROR: PETSc thermal model failed to converge '//&
               errMsg(__FILE__, __LINE__))
      endif

      ! Get the updated soil tempreature
      call thermal_mpp%sysofeqns%GetSoln(temperature_1d)

      ! Put temperature back in CLM structure for snow
      offset = 0

      do fc = 1, num_filter
         c = filter(fc)

         do j = -nlevsno+1, 0
            idx = (c-bounds%begc)*nlevsno + j + nlevsno + offset
            l = col%landunit(c)
            if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
               if (j >= snl(c)+1) then                  
                  tvector(c,j-1) = temperature_1d(idx)
               endif
            endif
         enddo
      enddo

      ! Put temperature back in CLM structure for standing water
      offset = (bounds%endc - bounds%begc + 1)*nlevsno
      do fc = 1, num_filter
         c = filter(fc)
         idx = (c-bounds%begc) + 1 + offset
         l = col%landunit(c)
         if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
            if (frac_h2osfc(c) > 0._r8) then                  
               tvector(c,0) = temperature_1d(idx)
            endif
         endif
      enddo

      ! Put temperature back in CLM structure for soil
      offset = (bounds%endc - bounds%begc + 1)*(nlevsno + 1)
      do fc = 1, num_filter
         c = filter(fc)
         do j = 1,nlevgrnd
            idx = (c-bounds%begc)*nlevgrnd + j + offset
            l = col%landunit(c)
            if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
               tvector(c,j) = temperature_1d(idx)
            endif
         enddo
      enddo

      ! Free up memory
      deallocate (temperature_1d         )
      deallocate (liq_areal_den_1d       )
      deallocate (ice_areal_den_1d       )
      deallocate (snow_water_1d          )
      deallocate (dz_1d                  )
      deallocate (dist_up_1d             )
      deallocate (dist_dn_1d             )
      deallocate (frac_1d                )
      deallocate (num_snow_layer_1d      )
      deallocate (is_active_1d           )
      deallocate (hs_snow_1d             )
      deallocate (hs_sh2o_1d             )
      deallocate (hs_soil_1d             )
      deallocate (dhsdT_snow_1d          )
      deallocate (dhsdT_sh2o_1d          )
      deallocate (dhsdT_soil_1d          )
      deallocate (frac_soil_1d           )
      deallocate (sabg_snow_1d           )
      deallocate (sabg_soil_1d           )
      deallocate (tsurf_tuning_factor_1d )

    end associate

  end subroutine MPPThermalTBasedALM_Solve

#endif

end module MPPThermalTBasedALM_Driver
