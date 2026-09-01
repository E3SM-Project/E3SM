module SoilMoistStressMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates soil moisture stress for plant gpp and transpiration
  !
  ! After discussion with other developers, I have now removed all functions that
  ! return array, and decalared all variables that will be modified as intent(inout).
  ! The initialization will be done whenever the variable is initialized. This avoids
  ! code crash when initialization is not done appropriately, and make the code safer
  ! during the long-term maintenance
  ! Created by Jinyun Tang, Feb., 2014
  !
  use ColumnDataType   , only : col_es, col_ws
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calc_root_moist_stress
  public :: calc_effective_soilporosity
  public :: calc_effective_snowporosity
  public :: calc_volumetric_h2oliq
  public :: set_perchroot_opt
  public :: init_root_moist_stress
  private :: normalize_test
  !
  ! !PRIVATE DATA MEMBERS:
  integer ::   root_moist_stress_method
  integer, parameter :: moist_stress_clm_default  = 0  !default method for calculating root moisture stress
  logical,  private :: perchroot     = .false.  ! true => btran is based only on unfrozen soil levels
  logical,  private :: perchroot_alt = .false.  ! true => btran is based on active layer (defined over two years);
  !$acc declare create(perchroot)
  !$acc declare create(perchroot_alt)

  !--------------------------------------------------------------------------------

contains

subroutine normalize_test( lbj2, ubj2, numf, filter, arr2d_inout)
   !DESCRIPTIONS
   !do normalization with filter for the input array along dimension 2
   !
   !USES
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none

   integer,  intent(in) :: lbj2         !right bound of dim 1
   integer,  intent(in) :: ubj2         !right bound of dim 2
   integer,  intent(in) :: numf         !filter size
   integer,  intent(in) :: filter(:)    !filter
   real(r8), intent(inout) :: arr2d_inout(: , : )   !input 2d array

   !local variables
   integer  :: sz1, sz2     !array size
   integer  :: j2           !indices
   integer  :: f, p         !indices
   real(r8) :: arr_sum(numf), sum1

   !$acc enter data create(arr_sum(1:numf))
   !$acc parallel loop independent gang worker default(present) private(sum1)
   do f = 1, numf
     sum1 = 0._r8
     !$acc loop vector reduction(+:sum1)
     do j2 = lbj2, ubj2
         !obtain the total
         sum1=sum1+arr2d_inout(f,j2)
     enddo
     arr_sum(f) = sum1
   enddo

     !normalize with the total if arr_sum is non-zero
   !$acc parallel loop independent gang default(present)
   do j2 = lbj2, ubj2
     !$acc loop vector independent
     do f = 1, numf
      !I found I have to ensure >0._r8 because of some unknown reason, jyt May 23, 2014
      !I will test this later with arr_sum(p)/=0._r8
      if(arr_sum(f)>0._r8 .or. arr_sum(f)<0._r8)then
         arr2d_inout(f,j2) = arr2d_inout(f,j2)/arr_sum(f)
      endif
     enddo
   enddo
   !$acc exit data delete(arr_sum(:))
end subroutine normalize_test

  !--------------------------------------------------------------------------------
  subroutine init_root_moist_stress()
    !
    !DESCRIPTION
    !specify the method to compute root soil moisture stress
    !
    implicit none

    root_moist_stress_method = moist_stress_clm_default
  end subroutine init_root_moist_stress

  !--------------------------------------------------------------------------------
  subroutine set_perchroot_opt(perchroot_global, perchroot_alt_global)
    !
    !DESCRIPTIONS
    !set up local perchroot logical switches, in the future, this wil be
    !read in as namelist
    !
    ! !ARGUMENTS:
    !$acc routine seq
    implicit none
    logical, intent(in) :: perchroot_global
    logical, intent(in) :: perchroot_alt_global
    !------------------------------------------------------------------------------

    perchroot = perchroot_global
    perchroot_alt = perchroot_alt_global

    !$acc update device(perchroot, perchroot_alt)
  end subroutine set_perchroot_opt

  !--------------------------------------------------------------------------------
  subroutine calc_effective_soilporosity(bounds, ubj, numf, filter, &
       watsat, h2osoi_ice, denice, eff_por)
    !
    ! !DESCRIPTIONS
    ! compute the effective soil porosity
    !
    ! !USES
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use decompMod      , only : bounds_type
    use ColumnType     , only : col_pp
    use VegetationType , only : veg_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                          ! bounds
    integer           , intent(in)    :: ubj                             ! lbinning level indices
    integer           , intent(in)    :: numf                            ! filter dimension
    integer           , intent(in)    :: filter(:)                       ! filter
    real(r8)          , intent(in)    :: watsat( bounds%begc:, 1: )     ! soil porosity
    real(r8)          , intent(in)    :: h2osoi_ice( bounds%begc:,1 : ) ! ice water content, kg H2o/m2
    real(r8)          , intent(in)    :: denice                          ! ice density, kg/m3
    real(r8)          , intent(inout) :: eff_por( bounds%begc: ,1: )     ! effective porosity
    !
    ! !LOCAL VARIABLES:
    integer :: c, j, fp,p      !indices
    real(r8):: vol_ice    !volumetric ice
    !------------------------------------------------------------------------------
    !main calculation loop
    !it assumes the soil layers start from 1
    !$acc parallel loop independent gang default(present)
    do j = 1, ubj
      !$acc loop vector independent private(p,c,vol_ice)
       do fp = 1, numf
          p = filter(fp)
          c = veg_pp%column(p)
          !compute the volumetric ice content
          vol_ice=min(watsat(c,j), h2osoi_ice(c,j)/(denice*col_pp%dz(c,j)))
          !compute the maximum soil space to fill liquid water and air
          eff_por(c,j) = watsat(c,j) - vol_ice
       enddo
    enddo

  end subroutine calc_effective_soilporosity

  !--------------------------------------------------------------------------------
  subroutine calc_effective_snowporosity(bounds, lbj, jtop, numf, filter, &
       h2osoi_ice, denice, eff_por)
    !
    ! !DESCRIPTIONS
    ! compute the effective porosity snow
    !
    ! !USES
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use decompMod      , only : bounds_type
    use ColumnType     , only : col_pp
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds                            !bounds
    integer           , intent(in)    :: lbj                               !ubing level indices
    integer           , intent(in)    :: jtop( bounds%begc: )              !top level for each column [col]
    integer           , intent(in)    :: numf                              !filter dimension
    integer           , intent(in)    :: filter(:)                         !filter
    real(r8)          , intent(in)    :: h2osoi_ice( bounds%begc: , lbj: ) !ice water content, kg H2o/m2
    real(r8)          , intent(in)    :: denice                            !ice density, kg/m3
    real(r8)          , intent(inout) :: eff_por( bounds%begc: ,lbj: )     !returning effective porosity
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j, fc    !indices
    integer  :: ubj
    real(r8) :: vol_ice     !volumetric ice
    !------------------------------------------------------------------------------

    ubj = 0

    !main calculation loop

    !it assumes snow layer ends at 0
    do j = lbj,0
       do fc = 1, numf
          c = filter(fc)
          if (j>=jtop(c)) then
             !compute the volumetric ice content
             vol_ice=min(1._r8, h2osoi_ice(c,j)/(denice*col_pp%dz(c,j)))

             !compute the maximum snow void space to fill liquid water and air
             eff_por(c,j) = 1._r8 - vol_ice
          endif
       enddo
    enddo

  end subroutine calc_effective_snowporosity

  !--------------------------------------------------------------------------------
  subroutine calc_volumetric_h2oliq(bounds, lbj, ubj, numf, filter,&
       eff_porosity, h2osoi_liq, denh2o, vol_liq)
    !
    ! !DESCRIPTIONS
    ! compute the volumetric liquid water content
    !
    !
    ! !USES
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use decompMod      , only : bounds_type
    use ColumnType     , only : col_pp
    use VegetationType  , only : veg_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                             ! bounds
    integer           , intent(in)    :: lbj, ubj                           ! lbinning and ubing level indices
    integer           , intent(in)    :: numf                               ! filter dimension
    integer           , intent(in)    :: filter(:)                          ! filter
    real(r8)          , intent(in)    :: eff_porosity(bounds%begc: , lbj: ) ! effective soil porosity
    real(r8)          , intent(in)    :: h2osoi_liq(bounds%begc: , lbj: )   ! liquid water content [kg H2o/m2]
    real(r8)          , intent(in)    :: denh2o                             ! water density [kg/m3]
    real(r8)          , intent(inout) :: vol_liq(bounds%begc: , lbj: )      ! volumetric liquid water content
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, j, fp ! indices
    integer, parameter :: jtop = 1
    !------------------------------------------------------------------------------
    !main calculation loop
    !$acc parallel loop independent gang default(present)
    do j = lbj, ubj
      !$acc loop vector independent private(p,c)
       do fp = 1, numf
          p = filter(fp)
          c = veg_pp%column(p)
          if(j>=jtop)then
             !volume of liquid is no greater than effective void space
             vol_liq(c,j) = min(eff_porosity(c,j), h2osoi_liq(c,j)/(col_pp%dz(c,j)*denh2o))
          endif
       enddo
    enddo

  end subroutine calc_volumetric_h2oliq

  !--------------------------------------------------------------------------------
  subroutine normalize_unfrozen_rootfr(bounds, ubj, fn, filterp, &
       canopystate_vars, soilstate_vars, rootfr_unf)
    !
    ! !DESCRIPTIONS
    ! normalize root fraction for total unfrozen depth
    !
    ! !USES
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use elm_varcon      , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment
    use decompMod       , only : bounds_type
    use CanopyStateType , only : canopystate_type
    use EnergyFluxType  , only : energyflux_type
    use SoilStateType   , only : soilstate_type
    use SimpleMathMod   , only : array_normalization
    use VegetationType  , only : veg_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds                                     !bounds
    integer                , intent(in)    :: ubj                                        !ubinning level indices
    integer                , intent(in)    :: fn                                         !filter dimension
    integer                , intent(in)    :: filterp(:)                                 !filter
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    real(r8)               , intent(inout) :: rootfr_unf(:,:) !normalized root fraction in unfrozen layers
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, j, f  !indices
    real(r8) :: arr_sum(fn), sum1
    !------------------------------------------------------------------------------

    associate(     &
         rootfr               => soilstate_vars%rootfr_patch , & ! Input:  [real(r8)  (:,:) ]  fraction of roots in each soil layer
         t_soisno             => col_es%t_soisno             , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

         altmax_lastyear_indx => canopystate_vars%altmax_lastyear_indx_col , & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth of thaw
         altmax_indx          => canopystate_vars%altmax_indx_col            & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw
         )

      ! main calculation loop
      ! Initialize rootfr_unf to zero.
      ! I found it necessary to ensure the pgi compiler not
      ! to complain with float point exception. However, it raises a question how
      ! to make sure those values that are initialized with nan or spval are not reset
      ! to zero within similar coding style. Jinyun Tang, May 23, 2014.

      ! Define rootfraction for unfrozen soil only
      ! if (perchroot .or. perchroot_alt) then
         if (perchroot_alt) then
            ! use total active layer (defined ass max thaw depth for current and prior year)
            !$acc parallel loop independent gang default(present)
            do j = 1, ubj
               !$acc loop vector independent private(p,c)
               do f = 1, fn
                  p = filterp(f)
                  c = veg_pp%column(p)

                  if ( j <= max(altmax_lastyear_indx(c), altmax_indx(c), 1) )then
                     rootfr_unf(f,j) = rootfr(p,j)
                  else
                     rootfr_unf(f,j) = 0._r8
                  end if
               end do
            end do
         else
            ! use instantaneous temperature
            !$acc parallel loop independent gang default(present)
            do j = 1, ubj
               !$acc loop vector independent private(p,c)
               do f = 1, fn
                  p = filterp(f)
                  c = veg_pp%column(p)
                  if (t_soisno(c,j) >= tfrz) then
                     rootfr_unf(f,j) = rootfr(p,j)
                  else
                     rootfr_unf(f,j) = 0._r8
                  end if
               end do
            end do

         end if ! perchroot_alt
      ! end if ! perchroot

      !normalize the root fraction for each pft
      ! call normalize_test( 1, ubj, &
      !      fn, filterp, rootfr_unf(:, :))
      ! !!!
      !$acc enter data create(arr_sum(1:fn))
      !$acc parallel loop independent gang worker default(present) private(sum1)
      do f = 1, fn
        sum1 = 0._r8
        !$acc loop vector reduction(+:sum1)
        do j = 1, ubj
            !obtain the total
            sum1=sum1+rootfr_unf(f,j)
        enddo
        arr_sum(f) = sum1
      enddo

        !normalize with the total if arr_sum is non-zero
      !$acc parallel loop independent gang default(present)
      do j = 1, ubj
        !$acc loop vector independent
        do f = 1, fn
         !I found I have to ensure >0._r8 because of some unknown reason, jyt May 23, 2014
         !I will test this later with arr_sum(p)/=0._r8
         if(arr_sum(f)>0._r8 .or. arr_sum(f)<0._r8)then
            rootfr_unf(f,j) = rootfr_unf(f,j)/arr_sum(f)
         endif
        enddo
      enddo
      !$acc exit data delete(arr_sum(:))
    end associate

  end subroutine normalize_unfrozen_rootfr

  !--------------------------------------------------------------------------------
  subroutine calc_root_moist_stress_clm45default( &
                        nlevgrnd, fn, filterp, rootfr_unf, &
                        soilstate_vars, energyflux_vars)
    !
    ! DESCRIPTIONS
    ! compute the root water stress using the default clm45 approach
    !
    ! USES
    use shr_kind_mod         , only : r8 => shr_kind_r8
    use elm_varcon           , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment
    use VegetationPropertiesType     , only : veg_vp
    use SoilStateType        , only : soilstate_type
    use EnergyFluxType       , only : energyflux_type
    use VegetationType            , only : veg_pp
    use elm_varctl       , only : use_hydrstress
    !
    ! !ARGUMENTS:
    implicit none
    integer                , intent(in)    :: nlevgrnd                       !number of vertical layers
    integer                , intent(in)    :: fn                             !number of filters
    integer                , intent(in)    :: filterp(:)                     !filter array
    real(r8)               , intent(in)    :: rootfr_unf(: , : )
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8) :: smp_node, s_node  !temporary variables
    real(r8) :: smp_node_lf,sum1, sum2      !temporary variable
    integer :: p, f, j, c, l      !indices
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes

    associate(                                                &
         smpso         => veg_vp%smpso                  , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal opening (mm)
         smpsc         => veg_vp%smpsc                  , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal closure (mm)
         tc_stress     => veg_vp%tc_stress              , & ! Input:  [real(r8)       ]  critical soil temperature for soil water stress (C)
         t_soisno      => col_es%t_soisno     , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

         watsat        => soilstate_vars%watsat_col         , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)
         sucsat        => soilstate_vars%sucsat_col         , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                        (constant)
         bsw           => soilstate_vars%bsw_col            , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                         (constant)
         eff_porosity  => soilstate_vars%eff_porosity_col   , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice
         rootfr        => soilstate_vars%rootfr_patch       , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer
         rootr         => soilstate_vars%rootr_patch        , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer

         btran         => energyflux_vars%btran_patch       , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1) (integrated soil water stress)
         btran2        => energyflux_vars%btran2_patch      , & ! Output: [real(r8) (:)   ]  integrated soil water stress square
         rresis        => energyflux_vars%rresis_patch      , & ! Output: [real(r8) (:,:) ]  root soil water stress (resistance) by layer (0-1)  (nlevgrnd)

         h2osoi_vol    => col_ws%h2osoi_vol    , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osoi_liqvol => col_ws%h2osoi_liqvol   & ! Output: [real(r8) (:,:) ]  liquid volumetric moisture, will be used for BeTR
         )

      !$acc parallel loop independent gang default(present)
      do j = 1,nlevgrnd
         !$acc loop vector independent private(p,c,l,s_node,smp_node )
         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            l = veg_pp%landunit(p)

            ! Root resistance factors
            ! rootr effectively defines the active root fraction in each layer
            if (h2osoi_liqvol(c,j) .le. 0._r8 .or. t_soisno(c,j) .le. tfrz + tc_stress) then
                    rootr(p,j) = 0._r8
            else
               s_node = max(h2osoi_liqvol(c,j)/eff_porosity(c,j),0.01_r8)
               !smp_node = max(smpsc(veg_pp%itype(p)), -sucsat(c,j)*s_node**(-bsw(c,j)))
               !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_node)
               smp_node = -sucsat(c,j)*s_node**( -bsw(c,j) )

               smp_node = max(smpsc(veg_pp%itype(p)), smp_node)

               rresis(p,j) = min( (eff_porosity(c,j)/watsat(c,j))* &
                    (smp_node - smpsc(veg_pp%itype(p))) / (smpso(veg_pp%itype(p)) - smpsc(veg_pp%itype(p))), 1._r8)
               if (.not. (perchroot .or. perchroot_alt) ) then
                  rootr(p,j) = rootfr(p,j)*rresis(p,j)
               else
                  rootr(p,j) = rootfr_unf(f,j)*rresis(p,j)
               end if
            endif

         end do
      end do

      !calculate btran and btran2
      if(.not. use_hydrstress) then
         !$acc parallel loop independent gang worker default(present) private(p,c,sum1)
         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            sum1 = btran(p)
               !$acc loop reduction(+:sum1) private(s_node,smp_node_lf)
               do j = 1, nlevgrnd
                  if (.not. (h2osoi_liqvol(c,j) .le. 0._r8 .or. t_soisno(c,j) .le. tfrz + tc_stress)) then
                     sum1 = sum1 + max(rootr(p,j),0._r8)
                  end if
               end do
               btran(p) = sum1
         end do
      end if

      !$acc parallel loop independent gang worker default(present) private(p,c,sum2)
      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         sum2 = btran2(p)
            !$acc loop reduction(+:sum2) private(s_node,smp_node_lf)
            do j = 1,nlevgrnd
               if (.not. (h2osoi_liqvol(c,j) .le. 0._r8 .or. t_soisno(c,j) .le. tfrz + tc_stress)) then
                  !smp_node_lf = max(smpsc(veg_pp%itype(p)), -sucsat(c,j)*(h2osoi_vol(c,j)/watsat(c,j))**(-bsw(c,j)))
                  s_node = h2osoi_vol(c,j)/watsat(c,j)
                  !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_node_lf)
                  smp_node_lf = -sucsat(c,j)*s_node**( -bsw(c,j) )
                  !smp_node_lf =  -sucsat(c,j)*(h2osoi_vol(c,j)/watsat(c,j))**(-bsw(c,j))
                  smp_node_lf = max(smpsc(veg_pp%itype(p)), smp_node_lf)
                  sum2  = sum2 +rootfr(p,j)*min((smp_node_lf - smpsc(veg_pp%itype(p))) / &
                           (smpso(veg_pp%itype(p)) - smpsc(veg_pp%itype(p))), 1._r8)
               end if
            end do
            btran2(p) = sum2
      end do
      !$acc wait

      ! Normalize root resistances to get layer contribution to ET
      !$acc parallel loop independent gang default(present)
      do j = 1,nlevgrnd
         !$acc loop vector independent private(p)
         do f = 1, fn
            p = filterp(f)
            if (btran(p) > btran0) then
               rootr(p,j) = rootr(p,j)/btran(p)
            else
               rootr(p,j) = 0._r8
            end if
         end do
      end do

    end associate

  end subroutine calc_root_moist_stress_clm45default

  !--------------------------------------------------------------------------------
  subroutine calc_root_moist_stress(bounds, nlevgrnd, fn, filterp, &
       canopystate_vars, energyflux_vars,  soilstate_vars)
    !
    ! DESCRIPTIONS
    ! compute the root water stress using different approaches
    !
    ! USES
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use elm_varcon      , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment
    use decompMod       , only : bounds_type
    use CanopyStateType , only : canopystate_type
    use EnergyFluxType  , only : energyflux_type
    use SoilStateType   , only : soilstate_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds   !bounds
    integer                , intent(in)    :: nlevgrnd
    integer                , intent(in)    :: fn
    integer                , intent(in)    :: filterp(:)
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p, f, j, c, l                ! indices
    real(r8) :: smp_node, s_node            ! temporary variables
    real(r8) :: rootfr_unf(1:fn,1:nlevgrnd) ! Rootfraction defined for unfrozen layers only.
    !------------------------------------------------------------------------------

    !define normalized rootfraction for unfrozen soil
    rootfr_unf(1:fn,1:nlevgrnd) = 0._r8
    !$acc enter data copyin(rootfr_unf(1:fn,1:nlevgrnd))

    call normalize_unfrozen_rootfr(bounds,  &
         ubj = nlevgrnd,                    &
         fn = fn,                           &
         filterp = filterp,                 &
         canopystate_vars=canopystate_vars, &
         soilstate_vars=soilstate_vars,     &
         rootfr_unf=rootfr_unf(1:fn,1:nlevgrnd))

    !suppose h2osoi_liq, eff_porosity are already computed somewhere else

    select case (root_moist_stress_method)
       !add other methods later
    case (moist_stress_clm_default)

       call calc_root_moist_stress_clm45default( &
            nlevgrnd = nlevgrnd,                        &
            fn = fn,                                    &
            filterp = filterp,                          &
            energyflux_vars=energyflux_vars,            &
            soilstate_vars=soilstate_vars,              &
            rootfr_unf=rootfr_unf(1:fn,1:nlevgrnd))

    case default
    end select
    !$acc exit data delete(rootfr_unf(:,:))

  end subroutine calc_root_moist_stress

end module SoilMoistStressMod
