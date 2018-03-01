module SurfaceResistanceMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceResistanceMod
!
! !DESCRIPTION:
! Module holding routines for calculation of surface resistances of the different tracers
! transported with BeTR. The surface here refers to water and soil, not including canopy
!
! this will be merged with CLM's surfaceResistanceMod to form the full integration of betr into clm

  use shr_kind_mod       , only: r8 => shr_kind_r8
implicit none
  private
  public :: calc_snow_pondw_resistance


  real(r8), parameter :: smallnumber = 1.e-12_r8
contains


  subroutine calc_snow_pondw_resistance(bounds, numf, filter, lbj, ubj, jtop, col, bulkdif_sno, pondwdiff, t_soisno, islake, &
      num_tracers, waterstate_vars, soilstate_vars, grnd_res)
  !
  ! DESCRIPTION
  ! After struggling with the partition of soil surface into bare soil, snow covered soil and
  ! water ponded soil, I finally decided to lump the snow layer togeter and calculate
  ! the so called snow layer resistance using the approach by  Zack and Bill, but with different
  ! parameters. However, when a mechanistic treatment of snow layer is ready,
  ! this subroutine will no longer be needed.
  ! Jinyun Tang, June 5, 2014.

  !
  ! USES
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use clm_varcon            , only : tfrz, denice, denh2o
  use decompMod             , only: bounds_type
  use SoilStatetype         , only : soilstate_type
  use WaterStateType        , only : Waterstate_Type
  use ColumnType            , only : column_type
  implicit none
  type(bounds_type)         , intent(in) :: bounds                               ! bounds
  integer                   , intent(in) :: numf                                 ! number of columns in column filter
  integer                   , intent(in) :: filter(:)                            ! column filter
  integer                   , intent(in) :: lbj, ubj                             ! lower and upper bounds
  integer                   , intent(in) :: jtop(:)                              !indices
  type(column_type)         , intent(in) :: col                                  !column type
  real(r8)                  , intent(in) :: t_soisno(bounds%begc:bounds%endc,lbj:ubj)  !soil temperature
  integer                   , intent(in) :: num_tracers        ! betr configuration information
  real(r8)                  , intent(in) :: bulkdif_sno(bounds%begc:bounds%endc, lbj:ubj, 1:num_tracers)
  real(r8)                  , intent(in) :: pondwdiff(bounds%begc:bounds%endc,1:num_tracers) ! tracer diffusivity in ponded water
  logical                   , intent(in) :: islake  ! logical indicate
  type(Waterstate_Type)     , intent(in) :: waterstate_vars        ! water state variables
  type(soilstate_type)      , intent(in) :: soilstate_vars        ! column physics variable
  real(r8)                  , intent(inout):: grnd_res(bounds%begc:bounds%endc, 1:num_tracers)

  !local varibles
  real(r8) :: pondz                            !ponding depth, m
  real(r8) :: pondres                          !ponded layer resistance
  real(r8) :: ponddiff                         !diffusivity in ponded water, (s/m2)
  real(r8) :: snowres(bounds%begc:bounds%endc)   !snow layer resistance, lumped (s/m)
  integer :: kk, j, fc, c   !indices

  associate(&
   h2osoi_vol    =>   waterstate_vars%h2osoi_vol_col              , & !Input [real(r8)(:,:)] volumetric soil water content
   h2osoi_liq    =>   waterstate_vars%h2osoi_liq_col              , & !Input [real(r8)(:,:)] volumetric liquid water content
   h2osoi_ice    =>   waterstate_vars%h2osoi_ice_col              , & !Input [real(r8)(:,:)] volumetric ice content
   watsat        =>   soilstate_vars%watsat_col                   , & !Input [real(r8)(:,:)] saturated water content
   dz            =>   col%dz                                        & !Input [real(r8)(:,:)] layer thickness
  )

  do kk = 1, num_tracers
    snowres(bounds%begc:bounds%endc)= 0._r8
    do j = lbj, ubj
      do fc = 1, numf
        c = filter(fc)
        ! Add snow resistance
        if (j >= jtop(c)) then
          snowres(c) = snowres(c) + dz(c,j)/bulkdif_sno(c,j,kk)
        end if

        if (j == 0) then ! final loop
          ! Add pond resistance
          pondres = 0._r8

          ! First old pond formulation up to pondmx
          if (.not. islake .and. jtop(c) == 1 .and. h2osoi_vol(c,1) > watsat(c,1)) then

            if (t_soisno(c,1) <= tfrz) then
              ponddiff = pondwdiff(c,kk) * (h2osoi_liq(c,1)/denh2o+smallnumber)/ &
                          (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice+smallnumber)
            else ! Unfrozen
              ponddiff = pondwdiff(c,kk)
            end if
            pondz = dz(c,1) * (h2osoi_vol(c,1) - watsat(c,1))
            pondres = pondz / ponddiff
          end if

          !this version does not turn on wetland or inundated fraction, Jinyun Tang, June 9, 2014
          ! Now add new h2osfc form
          !if (.not. lake .and. sat == 1 .and. frac_h2osfc(c) > 0._r8 .and. t_h2osfc(c) >= tfrz) then
          !        t_soisno_c = t_h2osfc(c) - tfrz
          !        ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
          !             * scale_factor_liqdiff
          !        pondz = h2osfc(c) / 1000._r8 / frac_h2osfc(c) ! Assume all h2osfc corresponds to sat area
          !        ! mm      /  mm/m
          !        pondres = pondres + pondz / ponddiff
          !else if (.not. lake .and. sat == 1 .and. frac_h2osfc(c) > 0._r8 .and. &
          !          h2osfc(c)/frac_h2osfc(c) > capthick) then ! Assuming short-circuit logic will avoid FPE here.
          !        ! assume surface ice is impermeable
          !        pondres = 1/smallnumber
          !end if

          grnd_res(c,kk) =  snowres(c) + pondres
        endif !  end if

      end do ! fc
    end do ! j
  enddo !kk
  end associate
  end subroutine calc_snow_pondw_resistance


end module SurfaceResistanceMod
