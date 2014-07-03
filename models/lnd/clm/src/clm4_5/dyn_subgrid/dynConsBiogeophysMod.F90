module dynConsBiogeophysMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeophysical quantities (water & energy) with dynamic land
  ! cover.
  !
  ! !USES:
  use clmtype
  use decompMod      , only : bounds_type
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dyn_hwcontent_init            ! compute grid-level heat and water content, before land cover change
  public :: dyn_hwcontent_final           ! compute grid-level heat and water content, after land cover change; also compute dynbal fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS
  private :: dyn_hwcontent                ! do the actual computation of grid-level heat and water content
  !---------------------------------------------------------------------------

contains

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_init(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize variables used for dyn_hwcontent, and compute grid cell-level heat
    ! and water content before land cover change
    !
    ! Should be called BEFORE any subgrid weight updates this time step
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g   ! grid cell index
    !-------------------------------------------------------------------------------
    
    ! initialize heat and water content and dynamic balance fields to zero
    do g = bounds%begg, bounds%endg
       gwf%qflx_liq_dynbal(g) = 0._r8
       gws%gc_liq2(g)         = 0._r8
       gws%gc_liq1(g)         = 0._r8
       gwf%qflx_ice_dynbal(g) = 0._r8
       gws%gc_ice2(g)         = 0._r8 
       gws%gc_ice1(g)         = 0._r8
       gef%eflx_dynbal(g)     = 0._r8
       ges%gc_heat2(g)        = 0._r8
       ges%gc_heat1(g)        = 0._r8
    enddo

    call dyn_hwcontent( bounds, &
         gws%gc_liq1(bounds%begg:bounds%endg), &
         gws%gc_ice1(bounds%begg:bounds%endg), &
         ges%gc_heat1(bounds%begg:bounds%endg) )

  end subroutine dyn_hwcontent_init

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_final(bounds)
    !
    ! !DESCRIPTION:
    ! Compute grid cell-level heat and water content after land cover change, and compute
    ! the dynbal fluxes
    !
    ! Should be called AFTER all subgrid weight updates this time step
    !
    ! !USES:
    use clm_time_manager    , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: g     ! grid cell index
    real(r8) :: dtime ! land model time step (sec)
    !---------------------------------------------------------------------------

    call dyn_hwcontent( bounds, &
         gws%gc_liq2(bounds%begg:bounds%endg), &
         gws%gc_ice2(bounds%begg:bounds%endg), &
         ges%gc_heat2(bounds%begg:bounds%endg) )

    dtime = get_step_size()
    do g = bounds%begg, bounds%endg
       gwf%qflx_liq_dynbal(g) = (gws%gc_liq2 (g) - gws%gc_liq1 (g))/dtime
       gwf%qflx_ice_dynbal(g) = (gws%gc_ice2 (g) - gws%gc_ice1 (g))/dtime
       gef%eflx_dynbal    (g) = (ges%gc_heat2(g) - ges%gc_heat1(g))/dtime
    end do

  end subroutine dyn_hwcontent_final

  
  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent(bounds, gcell_liq, gcell_ice, gcell_heat)

    ! !DESCRIPTION:
    ! Compute grid-level heat and water content to track conservation with respect to
    ! dynamic land cover.

    ! !USES:
    use clm_varcon, only : istsoil,istice,istwet,istdlak,istice_mec,istcrop,&
                           icol_road_perv,icol_road_imperv,icol_roof,icol_sunwall,icol_shadewall,&
                           cpice,  cpliq, denh2o
    use clm_varpar, only : nlevsno, nlevgrnd, nlevurb, nlevlak
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    real(r8), intent(out) :: gcell_liq ( bounds%begg: )   ! [gridcell]
    real(r8), intent(out) :: gcell_ice ( bounds%begg: )   ! [gridcell]
    real(r8), intent(out) :: gcell_heat( bounds%begg: )   ! [gridcell]
    !
    ! !LOCAL VARIABLES:
    integer  :: li,lf         ! loop initial/final indicies
    integer  :: ci,cf         ! loop initial/final indicies
    integer  :: pi,pf         ! loop initial/final indicies

    integer  :: g,l,c,p,k     ! loop indicies (grid,lunit,column,pft,vertical level)

    real(r8) :: wtgcell       ! weight relative to grid cell
    real(r8) :: wtcol         ! weight relative to column
    real(r8) :: liq           ! sum of liquid water at column level
    real(r8) :: ice           ! sum of frozen water at column level
    real(r8) :: heat          ! sum of heat content at column level
    real(r8) :: cv            ! heat capacity [J/(m^2 K)]

    integer ,pointer :: nlev_improad(:)  ! number of impervious road layers
    real(r8),pointer :: cv_wall(:,:)     ! thermal conductivity of urban wall
    real(r8),pointer :: cv_roof(:,:)     ! thermal conductivity of urban roof
    real(r8),pointer :: cv_improad(:,:)  ! thermal conductivity of urban impervious road
    integer ,pointer :: snl(:)           ! number of snow layers
    real(r8),pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(r8),pointer :: h2osno(:)        ! snow water (mm H2O)
    real(r8),pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8),pointer :: h2osoi_ice(:,:)  ! frozen water (kg/m2)
    real(r8),pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(r8),pointer :: csol(:,:)        ! heat capacity, soil solids (J/m**3/Kelvin)
    real(r8),pointer :: dz(:,:)          ! layer depth (m)
    !-------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(gcell_liq)  == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(gcell_ice)  == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(gcell_heat) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))

    nlev_improad => lps%nlev_improad
    cv_wall      => lps%cv_wall
    cv_roof      => lps%cv_roof
    cv_improad   => lps%cv_improad
    snl          => cps%snl
    watsat       => cps%watsat
    csol         => cps%csol
    dz           => cps%dz
    t_soisno     => ces%t_soisno
    h2osoi_liq   => cws%h2osoi_liq
    h2osoi_ice   => cws%h2osoi_ice
    h2osno       => cws%h2osno

    ! Get relevant sizes

    do g = bounds%begg,bounds%endg ! loop over grid cells
       gcell_liq  (g) = 0.0_r8   ! sum for one grid cell
       gcell_ice  (g) = 0.0_r8   ! sum for one grid cell
       gcell_heat (g) = 0.0_r8   ! sum for one grid cell
    end do

    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       ci = lun%coli(l)
       cf = lun%colf(l)
       do c = ci,cf   ! loop over columns

          liq   = 0.0_r8 ! sum for one column
          ice   = 0.0_r8
          heat  = 0.0_r8

          !--- water & ice, above ground only ---
          if ( (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop         )  &
               .or. (lun%itype(l) == istwet                                   )  &
               .or. (lun%itype(l) == istice                                   )  &
               .or. (lun%itype(l) == istice_mec                               )  &           
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_roof       )  &
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_imperv)  &
               .or. (lun%itype(l) == istdlak                                  )  &
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_perv  )) then

             if ( snl(c) < 0 ) then
                do k = snl(c)+1,0 ! loop over snow layers
                   liq   = liq   + cws%h2osoi_liq(c,k)
                   ice   = ice   + cws%h2osoi_ice(c,k)
                end do
             else                 ! no snow layers exist
                ice = ice + cws%h2osno(c)
             end if
          end if

          !--- water & ice, below ground only ---
          if ( (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop         )  &
               .or. (lun%itype(l) == istwet                                   )  &
               .or. (lun%itype(l) == istice                                   )  &
               .or. (lun%itype(l) == istdlak                                  )  &
               .or. (lun%itype(l) == istice_mec                               )  &           
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_perv  )) then
             do k = 1,nlevgrnd
                liq   = liq   + cws%h2osoi_liq(c,k)
                ice   = ice   + cws%h2osoi_ice(c,k)
             end do
          end if

          !--- water & ice, below ground, for lakes ---
          if ( lun%itype(l) == istdlak ) then
             do k = 1,nlevlak
                liq   = liq   + (1 - cws%lake_icefrac(c,k))*cps%dz_lake(c,k)*denh2o
                ice   = ice   + cws%lake_icefrac(c,k)*cps%dz_lake(c,k)*denh2o
                ! lake layers do not change thickness when freezing, so denh2o should be used
                ! (thermal properties are appropriately adjusted; see SLakeTemperatureMod)
             end do
          end if

          !--- water in aquifer ---
          if ( (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop         )  &
               .or. (lun%itype(l) == istwet                                   )  &
               .or. (lun%itype(l) == istice                                   )  &
               .or. (lun%itype(l) == istice_mec                               )  &           
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_perv  )) then
             liq = liq + cws%wa(c)
          end if

          !--- water in canopy (at pft level) ---
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then   ! note: soil specified at LU level
             do p = col%pfti(c),col%pftf(c) ! loop over pfts
                if (pft%active(p)) then
                   liq = liq + pws%h2ocan(p) * pft%wtcol(p)
                end if
             end do
          end if

          !--- heat content, below ground only ---
          if (nlevurb > 0) then
             do k = 1,nlevurb
                if (col%itype(c)==icol_sunwall .OR. col%itype(c)==icol_shadewall) then
                   cv = cv_wall(l,k) * dz(c,k)
                   heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
                else if (col%itype(c) == icol_roof) then
                   cv = cv_roof(l,k) * dz(c,k)
                   heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
                end if
             end do
          end if
          do k = 1,nlevgrnd
             if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                  .and. col%itype(c) /= icol_roof) then
                if (col%itype(c) == icol_road_imperv .and. k >= 1 .and. k <= nlev_improad(l)) then
                   cv = cv_improad(l,k) * dz(c,k)
                else if (lun%itype(l) /= istwet .AND. lun%itype(l) /= istice .AND. lun%itype(l) /= istice_mec) then
                   cv = csol(c,k)*(1-watsat(c,k))*dz(c,k) + (h2osoi_ice(c,k)*cpice + h2osoi_liq(c,k)*cpliq)
                else
                   cv = (h2osoi_ice(c,k)*cpice + h2osoi_liq(c,k)*cpliq)
                endif
                heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
             end if
          end do

          !--- heat content, below ground in lake water, for lakes ---
          do k = 1,nlevlak
             if (lun%itype(l) == istdlak) then
                cv = denh2o*cps%dz_lake(c,k)*( cws%lake_icefrac(c,k)*cpice + &
                     (1 - cws%lake_icefrac(c,k))*cpliq )
                heat = heat + cv*ces%t_lake(c,k) / 1.e6_r8
             end if
          end do

          !--- heat content, above ground only ---
          if ( snl(c) < 0 ) then
             do k = snl(c)+1,0 ! loop over snow layers
                cv = cpliq*h2osoi_liq(c,k) + cpice*h2osoi_ice(c,k)
                heat = heat + cv*t_soisno(c,k) / 1.e6_r8
             end do
          else if ( h2osno(c) > 0.0_r8 .and. lun%itype(l) /= istdlak) then
             ! the heat capacity (not latent heat) of snow without snow layers
             ! is currently ignored in SLakeTemperature, so it should be ignored here
             k = 1
             cv = cpice*h2osno(c)
             heat = heat + cv*t_soisno(c,k) / 1.e6_r8
          end if

          !--- scale x/m^2 column-level values into x/m^2 gridcell-level values ---
          gcell_liq  (g) = gcell_liq  (g) + liq   * col%wtgcell(c)
          gcell_ice  (g) = gcell_ice  (g) + ice   * col%wtgcell(c)
          gcell_heat (g) = gcell_heat (g) + heat  * col%wtgcell(c)

       end do ! column loop      
    end do ! landunit loop

  end subroutine dyn_hwcontent

end module dynConsBiogeophysMod
