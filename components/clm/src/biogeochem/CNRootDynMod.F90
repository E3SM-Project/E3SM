module CNRootDynMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used for determining fine root distribution for all pfts.
  ! Includes dynamic root depth for crops
  !
  ! !USES:
  use shr_kind_mod        , only: r8 => shr_kind_r8
  use clm_time_manager    , only: get_step_size
  use clm_varpar          , only: nlevsoi, nlevgrnd
  use clm_varctl          , only: use_vertsoilc
  use decompMod           , only: bounds_type
  use pftvarcon           , only: noveg, npcropmin, roota_par, rootb_par,root_dmx
  use ColumnType          , only: col
  use PatchType           , only: pft
  use CNStateType         , only: cnstate_type
  use CNCarbonStateType   , only: carbonstate_type
  use CNCarbonFluxType    , only: carbonflux_type
  use CNNitrogenStateType , only: nitrogenstate_type
  use SoilStateType       , only: soilstate_type
  use CropType            , only: crop_type
  use SimpleMathMod       , only: array_normalization

  ! !PUBLIC TYPES:
  implicit none
  save
  private
  public :: CNRootDyn
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  !
  subroutine CNRootDyn(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonstate_vars, nitrogenstate_vars, carbonflux_vars,                    &
       cnstate_vars, crop_vars,  soilstate_vars)
    !
    ! !DESCRIPTION:
    ! This routine determine the fine root distribution
    ! Needs to be called after the photosynthesis calculation
    ! May need to update other subroutines that use the fixed root profile for calculations
    ! i.e. CNVerticalProfileMod
    !
    ! !USES:


    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds             ! bounds
    integer                  , intent(in)    :: num_soilc
    integer                  , intent(in)    :: filter_soilc(:)
    integer                  , intent(in)    :: num_soilp          ! number of soil pfts in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil pfts
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonflux_type)    , intent(in)    :: carbonflux_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars !
    type(crop_type)          , intent(in)    :: crop_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars

    !
    ! !LOCAL VARIABLES:

    integer  :: f,c,p,lev,j                                    ! indices
    real(r8) :: dt                                             ! radiation time step delta t (seconds)
    real(r8) :: w_limit(bounds%begp:bounds%endp)               ! soil water weighting factor
    real(r8) :: rswa(bounds%begp:bounds%endp,1:nlevgrnd)       ! soil water availability in each soil layer
    real(r8) :: rsmn(bounds%begp:bounds%endp,1:nlevgrnd)       ! soil nitrogen availability in each soil layer
    real(r8) :: sumrswa(bounds%begp:bounds%endp)               ! scaling soil water availability in each soil layer
    real(r8) :: sumrsmn(bounds%begp:bounds%endp)               ! scaling  soil mineral N availability in each soil layer
    real(r8) :: frootc_dz(bounds%begp:bounds%endp, 1:nlevgrnd) ! root carbon in each soil layer (gC)
    real(r8) :: sumfrootc(bounds%begp:bounds%endp)             ! fine root carbon total before turnover in each step
    real(r8) :: psi                                            ! soil moisture potential
    real(r8) :: maxpsi                                         ! maximum soil moisture potential
    real(r8) :: new_growth                                     ! new carbon allocated to roots this timestep
    real(r8), parameter :: minpsi = -1.5_r8                    ! minimum soil moisture potential - permanent wilting point (MPa)

    !-----------------------------------------------------------------------
    ! Assign local pointers to derived type arrays (in)
    associate(&
         ivt                    => pft%itype                                   , & ! Input  :  [integer (:)]  pft vegetation type
         pcolumn                => pft%column                                  , & ! Input  :  [integer (:)]  pft's column index
         croplive               => cnstate_vars%croplive_patch                 , & ! Input  :  [logical (:)]  flag, true if planted, not harvested
         cpool_to_frootc        => carbonflux_vars%cpool_to_frootc_patch       , & ! Input  :  [real(r8) (:)] allocation to fine root C (gC/m2/s)
         frootc_xfer_to_frootc  => carbonflux_vars%frootc_xfer_to_frootc_patch , & ! Input  :  [real(r8) (:)] fine root C growth from storage (gC/m2/s)
         dormant_flag           => cnstate_vars%dormant_flag_patch             , & ! Input  :  [real(r8) (:)]  dormancy flag
         root_depth             => soilstate_vars%root_depth_patch             , & ! InOut  :  [real(r8) (:)] current root depth
         dz                     => col%dz                                      , & ! Input  :  layer thickness (m)  (-nlevsno+1:nlevgrnd)
         zi                     => col%zi                                      , & ! Input  :  interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
         rootfr                 => soilstate_vars%rootfr_patch                 , & ! Output :  [real(r8) (:,:)]  fraction of roots in each soil layer
         sucsat                 => soilstate_vars%sucsat_col                   , & ! Input  :  minimum soil suction (mm)
         soilpsi                => soilstate_vars%soilpsi_col                  , & ! Input  :  soil water potential in each soil layer (MPa)
         sminn_vr               => nitrogenstate_vars%sminn_vr_col             , & ! Iniput :  [real(r8) (:,:)]  (gN/m3) soil mineral N
         frootc                 => carbonstate_vars%frootc_patch               , & ! Input  :  [real(r8) (:)]  (gC/m2) fine root C
         hui                    => crop_vars%gddplant_patch                    , & ! Input  :  [real(r8) (:)]  =gdd since planting (gddplant)
         huigrain               => cnstate_vars%huigrain_patch                   & ! Input  :  [real(r8) (:)]  same to reach vegetative maturity
         )

      ! set time steps
      dt = get_step_size()

      !initialize to 0
      w_limit(bounds%begp:bounds%endp)              = 0._r8
      sumrswa(bounds%begp:bounds%endp)              = 0._r8
      sumrsmn(bounds%begp:bounds%endp)              = 0._r8
      sumfrootc(bounds%begp:bounds%endp)            = 0._r8
      rswa(bounds%begp:bounds%endp,1:nlevgrnd)      = 0._r8
      rsmn(bounds%begp:bounds%endp,1:nlevgrnd)      = 0._r8
      frootc_dz(bounds%begp:bounds%endp,1:nlevgrnd) = 0._r8
      root_depth(bounds%begp:bounds%endp)           = 0._r8

      !---------------------------------------------------------------
      ! Set root depth, dynamic for crops, fixed for other vegetation
      !---------------------------------------------------------------

      do f = 1, num_soilp
         p = filter_soilp(f)
         c = pcolumn(p)
         if (ivt(p) /= noveg) then
            if ((ivt(p)) >= npcropmin) then !skip generic crop types
               if (huigrain(p) > 0._r8) then
                  root_depth(p) = max(zi(c,2), min(hui(p)/huigrain(p)* root_dmx(ivt(p)), root_dmx(ivt(p))))
               end if
            else
               ! this can be changed to any depth (i.e. the maximum soil depth)
               root_depth(p) = zi(c,nlevsoi)
            end if
         else
            root_depth(p) = 0._r8
         end if
      end do

      !----------------------------------------------------------------
      !  ! calculate a weighting function by soil depth that depends on the
      ! fine root distribution per pft and depth and the pft weight on the column.
      ! This will be used to weight the temperature and water potential scalars
      ! for decomposition control.

      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
      ! a comparison of models. Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

      do j = 1,nlevsoi
         do f = 1,num_soilp
            p = filter_soilp(f)
            c = pcolumn(p)
            maxpsi = sucsat(c,j) * (-9.8e-6_r8)
            psi = min(soilpsi(c,j),maxpsi)
            if (psi > minpsi) then

               ! First calculate water in the root zone
               if (root_depth(p) >  zi(c,3) .and. (zi(c,j) <= root_depth(p) .or. &
                    (zi(c,j-1) < root_depth(p) .and. zi(c,j) > root_depth(p)))) then
                  w_limit(p) = w_limit(p) + max(0._r8,log(minpsi/psi)/log(minpsi/maxpsi))*rootfr(p,j)
               end if
               ! Calculate the water in each soil layer
               if (root_depth(p) >= zi(c,j) .or. &
                    (zi(c,j-1) < root_depth(p) .and. zi(c,j) > root_depth(p))) then
                  rswa(p,j) = max(0._r8, (log(minpsi/psi)/log(minpsi/maxpsi)))
               end if
            end if
            sumrswa(p) = sumrswa(p) + rswa(p,j)

            ! Calculate the nitrogen profile in each layer
            ! For now, the profile for each PFT is equivilent to the
            ! column profile, in the future, this could be changed to a weighted profile
            if (use_vertsoilc) then !for vertical soil profile
               rsmn(p,j) = sminn_vr(c,j)
            else ! need to calculate a profile, top three layers are constant, and decrease linearly
               if (j <= 3) then
                  rsmn(p,j) = dz(c,j)
               end if
               if (j > 3) then
                  rsmn(p,j) = dz(c,j) * (zi(c,10) - zi(c,j)) / (zi(c,10) - zi(c,3))
               end if
            end if
            if (root_depth(p) >= zi(c,j).or. &
                 (zi(c,j-1) < root_depth(p) .and. zi(c,j) > root_depth(p))) then
               sumrsmn(p) = sumrsmn(p) + rsmn(p,j)
            end if
         end do
      end do


      !--------------------------------------------------------------------
      ! Now calculate the density of roots in each soil layer for each pft
      ! based on this timesteps growth
      !--------------------------------------------------------------------
      do lev = 1, nlevgrnd

         do f = 1, num_soilp
            p = filter_soilp(f) 
            c = pcolumn(p)

            new_growth = (cpool_to_frootc(p) + frootc_xfer_to_frootc(p))*dt
            if (zi(c,lev) <= root_depth(p) .or. &
                 (zi(c,lev-1) < root_depth(p) .and. zi(c,lev) > root_depth(p))) then
               if (sumrswa(p) <= 0._r8 .or. sumrsmn(p) <= 0._r8) then
                  ! when sumrswa or sumrsmn are less than or equal to 0 rootfr will not be updated
               else
                  frootc_dz(p,lev) = (frootc(p))*rootfr(p,lev) &
                       + new_growth * ((1._r8 - w_limit(p)) * rswa(p,lev) / sumrswa(p) &
                       + w_limit(p) * rsmn(p,lev) / sumrsmn(p))
               end if
            else
               frootc_dz(p,lev) = 0._r8
            end if

            sumfrootc(p) = sumfrootc(p) + frootc_dz(p,lev)

         end do
      end do

      !----------------------------------
      !Calculate root fraction
      !----------------------------------

      ! normalize the root fraction for each pft

      do lev = 1, nlevgrnd
         do f = 1, num_soilp
            p = filter_soilp(f)
            c = pcolumn(p)
            if (sumfrootc(p) > 0._r8) then
               rootfr(p,lev) = frootc_dz(p,lev)/sumfrootc(p)
            end if
            if (ivt(p) >= npcropmin .and. .not. croplive(p)) then
               ! CROPS are dormant, there are no roots!
               ! but, need an initial frootr so crops can start root production
               ! default is just the exponential over the top two soil layers
               if (lev <  2) then
                  rootfr(p,lev) = .5_r8*( exp(-roota_par(pft%itype(p)) * zi(c,lev-1))  &
                       + exp(-rootb_par(pft%itype(p)) * zi(c,lev-1))  &
                       - exp(-roota_par(pft%itype(p)) * zi(c,lev  ))  &
                       - exp(-rootb_par(pft%itype(p)) * zi(c,lev  )) )
               elseif (lev == 2) then
                  rootfr(p,lev) = .5_r8*( exp(-roota_par(pft%itype(p)) * zi(c,lev-1))  &
                       + exp(-rootb_par(pft%itype(p)) * zi(c,lev-1)) )
               else
                  rootfr(p,lev) =  0.0_r8
               end if

            end if
         end do
      end do


    end associate

  end subroutine CNRootDyn

end module CNRootDynMod
