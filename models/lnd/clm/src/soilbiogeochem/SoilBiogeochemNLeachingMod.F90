module SoilBiogeochemNLeachingMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for mineral nitrogen dynamics (deposition, fixation, leaching)
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use decompMod                       , only : bounds_type
  use clm_varcon                      , only : dzsoi_decomp, zisoi
  use clm_varctl                      , only : use_nitrif_denitrif, use_vertsoilc
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use WaterStateType                  , only : waterstate_type
  use WaterFluxType                   , only : waterflux_type
  use ColumnType                      , only : col                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: SoilBiogeochemNLeaching
  !
  ! !PRIVATE DATA:
  type, private :: params_type
     real(r8):: sf        ! soluble fraction of mineral N (unitless)
     real(r8):: sf_no3    ! soluble fraction of NO3 (unitless)
  end type params_type
  
  type(params_type), private ::  params_inst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use ncdio_pio   , only : file_desc_t,ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNNDynamicsParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    
    tString='sf_minn'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%sf=tempr

    tString='sf_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%sf_no3=tempr
   
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemNLeaching(bounds, num_soilc, filter_soilc, &
       waterstate_inst, waterflux_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen leaching rate
    ! as a function of soluble mineral N and total soil water outflow.
    !
    ! !USES:
    use clm_varpar       , only : nlevdecomp, nlevsoi
    use clm_time_manager , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(waterflux_type)                    , intent(in)    :: waterflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,fc                                 ! indices
    real(r8) :: dt                                     ! radiation time step (seconds)
    real(r8) :: sf                                     ! soluble fraction of mineral N (unitless)
    real(r8) :: sf_no3                                 ! soluble fraction of NO3 (unitless)
    real(r8) :: disn_conc                              ! dissolved mineral N concentration (gN/kg water)
    real(r8) :: tot_water(bounds%begc:bounds%endc)     ! total column liquid water (kg water/m2)
    real(r8) :: surface_water(bounds%begc:bounds%endc) ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot(bounds%begc:bounds%endc)     ! total drainage flux (mm H2O /s)
    real(r8), parameter :: depth_runoff_Nloss = 0.05   ! (m) depth over which runoff mixes with soil water for N loss to runoff
    !-----------------------------------------------------------------------

    associate(                                                                             & 
         h2osoi_liq          => waterstate_inst%h2osoi_liq_col                           , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

         qflx_drain          => waterflux_inst%qflx_drain_col                            , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf           => waterflux_inst%qflx_surf_col                             , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
         
         sminn_vr            => soilbiogeochem_nitrogenstate_inst%sminn_vr_col           , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                          
         smin_no3_vr         => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col        , & ! Input:  [real(r8) (:,:) ]                                                  

         sminn_leached_vr    => soilbiogeochem_nitrogenflux_inst%sminn_leached_vr_col    , & ! Output: [real(r8) (:,:) ]  rate of mineral N leaching (gN/m3/s)            
         smin_no3_leached_vr => soilbiogeochem_nitrogenflux_inst%smin_no3_leached_vr_col , & ! Output: [real(r8) (:,:) ]  rate of mineral NO3 leaching (gN/m3/s)          
         smin_no3_runoff_vr  => soilbiogeochem_nitrogenflux_inst%smin_no3_runoff_vr_col    & ! Output: [real(r8) (:,:) ]  rate of mineral NO3 loss with runoff (gN/m3/s)  
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      if (.not. use_nitrif_denitrif) then
         ! set constant sf 
         sf = params_inst%sf
      else
         ! Assume that 100% of the soil NO3 is in a soluble form
         sf_no3 =  params_inst%sf_no3 
      end if

      ! calculate the total soil water
      tot_water(bounds%begc:bounds%endc) = 0._r8
      do j = 1,nlevsoi
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            tot_water(c) = tot_water(c) + h2osoi_liq(c,j)
         end do
      end do

      ! for runoff calculation; calculate total water to a given depth
      surface_water(bounds%begc:bounds%endc) = 0._r8
      do j = 1,nlevsoi
         if ( zisoi(j) <= depth_runoff_Nloss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j)
            end do
         elseif ( zisoi(j-1) < depth_runoff_Nloss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * ( (depth_runoff_Nloss - zisoi(j-1)) / col%dz(c,j))
            end do
         endif
      end do

      ! Loop through columns
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         drain_tot(c) = qflx_drain(c)
      end do


      if (.not. use_nitrif_denitrif) then

         !----------------------------------------
         ! --------- NITRIF_NITRIF OFF------------
         !----------------------------------------

         do j = 1,nlevdecomp
            ! Loop through columns
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (.not. use_vertsoilc) then
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (tot_water(c) > 0._r8) then
                     disn_conc = (sf * sminn_vr(c,j) ) / tot_water(c)
                  end if

                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminn_leached_vr(c,j) = disn_conc * drain_tot(c)
               else
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (h2osoi_liq(c,j) > 0._r8) then
                     disn_conc = (sf * sminn_vr(c,j) * col%dz(c,j) )/(h2osoi_liq(c,j) )
                  end if

                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminn_leached_vr(c,j) = disn_conc * drain_tot(c) * h2osoi_liq(c,j) / ( tot_water(c) * col%dz(c,j) )

               end if

               ! limit the flux based on current sminn state
               ! only let at most the assumed soluble fraction
               ! of sminn be leached on any given timestep
               sminn_leached_vr(c,j) = min(sminn_leached_vr(c,j), (sf * sminn_vr(c,j))/dt)

               ! limit the flux to a positive value
               sminn_leached_vr(c,j) = max(sminn_leached_vr(c,j), 0._r8)

            end do
         end do

      else     

         !----------------------------------------
         ! --------- NITRIF_NITRIF ON-------------
         !----------------------------------------

         do j = 1,nlevdecomp
            ! Loop through columns
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (.not. use_vertsoilc) then
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (tot_water(c) > 0._r8) then
                     disn_conc = (sf_no3 * smin_no3_vr(c,j) )/tot_water(c)
                  end if

                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  smin_no3_leached_vr(c,j) = disn_conc * drain_tot(c)
               else
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (h2osoi_liq(c,j) > 0._r8) then
                     disn_conc = (sf_no3 * smin_no3_vr(c,j) * col%dz(c,j) )/(h2osoi_liq(c,j) )
                  end if
                  !
                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  smin_no3_leached_vr(c,j) = disn_conc * drain_tot(c) * h2osoi_liq(c,j) / ( tot_water(c) * col%dz(c,j) )
                  !
                  ! ensure that leaching rate isn't larger than soil N pool
                  smin_no3_leached_vr(c,j) = min(smin_no3_leached_vr(c,j), smin_no3_vr(c,j) / dt )
                  !
                  ! limit the leaching flux to a positive value
                  smin_no3_leached_vr(c,j) = max(smin_no3_leached_vr(c,j), 0._r8)
                  !
                  !
                  ! calculate the N loss from surface runoff, assuming a shallow mixing of surface waters into soil and removal based on runoff
                  if ( zisoi(j) <= depth_runoff_Nloss )  then
                     smin_no3_runoff_vr(c,j) = disn_conc * qflx_surf(c) * &
                          h2osoi_liq(c,j) / ( surface_water(c) * col%dz(c,j) )
                  elseif ( zisoi(j-1) < depth_runoff_Nloss )  then
                     smin_no3_runoff_vr(c,j) = disn_conc * qflx_surf(c) * &
                          h2osoi_liq(c,j) * ((depth_runoff_Nloss - zisoi(j-1)) / &
                          col%dz(c,j)) / ( surface_water(c) * (depth_runoff_Nloss-zisoi(j-1) ))
                  else
                     smin_no3_runoff_vr(c,j) = 0._r8
                  endif
                  !
                  ! ensure that runoff rate isn't larger than soil N pool
                  smin_no3_runoff_vr(c,j) = min(smin_no3_runoff_vr(c,j), smin_no3_vr(c,j) / dt - smin_no3_leached_vr(c,j))
                  !
                  ! limit the flux to a positive value
                  smin_no3_runoff_vr(c,j) = max(smin_no3_runoff_vr(c,j), 0._r8)


               endif
               ! limit the flux based on current smin_no3 state
               ! only let at most the assumed soluble fraction
               ! of smin_no3 be leached on any given timestep
               smin_no3_leached_vr(c,j) = min(smin_no3_leached_vr(c,j), (sf_no3 * smin_no3_vr(c,j))/dt)

               ! limit the flux to a positive value
               smin_no3_leached_vr(c,j) = max(smin_no3_leached_vr(c,j), 0._r8)

            end do
         end do
      endif

    end associate

  end subroutine SoilBiogeochemNLeaching

end module SoilBiogeochemNLeachingMod
