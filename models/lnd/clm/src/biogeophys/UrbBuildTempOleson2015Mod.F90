module UrbBuildTempOleson2015Mod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates internal building air temperature
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use perf_mod          , only : t_startf, t_stopf
  use clm_varctl        , only : iulog
  use UrbanParamsType   , only : urbanparams_type  
  use EnergyFluxType    , only : energyflux_type
  use TemperatureType   , only : temperature_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: BuildingTemperature ! Calculation of interior building air temperature, inner 
                                 ! surface temperatures of walls and roof, and floor temperature
  !-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildingTemperature
!
! !INTERFACE:
  subroutine BuildingTemperature (bounds, num_urbanl, filter_urbanl, num_nolakec, &
                                  filter_nolakec, tk, urbanparams_inst, temperature_inst, &
                                  energyflux_inst )
!
! !DESCRIPTION:
! Solve for t_building, inner surface temperatures of roof, sunw, shdw, and floor temperature
! Five equations, five unknowns (t_roof_inner,t_sunw_inner,t_shdw_inner,t_floor,t_building at n+1)
! Derived from energy balance equations at each surface and building air
! rd (radiation), cd (conduction), cv (convection)
! qrd_roof + qcd_roof + qcv_roof = 0
! qrd_sunw + qcd_sunw + qcv_sunw = 0
! qrd_shdw + qcd_shdw + qcv_shdw = 0
! qrd_floor + qcd_floor + qcv_floor = 0
! Vbld*rho_dair*cpair*(dt_building/dt) = sum(Asfc*hcv_sfc*(t_sfc - t_building) 
!                                        + Vvent*rho_dair*cpair*(taf - t_building)
!   where Vlbd is volume of building air,
!         rho_dair is density of dry air at t_building (kg m-3),
!         cpair is specific heat of dry air (J kg-1 K-1),
!         dt_building is change in interior building temperature (K),
!         dt is timestep (s),
!         Asfc is surface area of roof, sunw, shdw, floor (m2)
!         hcv_sfc is convective heat transfer coefficient for roof, sunw, shdw, floor (W m-2 K-1)
!         t_sfc is inner surface temperature of roof, sunw, shdw, floor (K)
!         t_building is interior building temperature (K)
!         Vvent is ventilation airflow rate (m3 s-1)
!         taf is urban canyon air temperature (K)
!
! This methodology was introduced as part of CLM5.0. 
!
! Conduction fluxes are obtained from terms of soil temperature equations
! Radiation fluxes are obtained from linearizing the longwave radiation equations taking into
! account view factors for each surface.

! qrd is positive away from the surface toward room air, so qrd = emitted - absorbed,
!   so positive qrd will result in a decrease in temperature
! qcd_floor is positive away from surface toward room air, so positive
!   qcd will result in a decrease in temperature 
! qcv is positive toward room air, so positive qcv (t_surface > t_room) will 
!   result in a decrease in temperature

! The LAPACK routine DGESV is used to compute the solution to the real system of linear equations
!     a * x = b,
!  where a is an n-by-n matrix and x and b are n-by-nrhs matrices.
! 
!  The LU decomposition with partial pivoting and row interchanges is
!  used to factor a as
!     a = P * L * U,
!  where P is a permutation matrix, L is unit lower triangular, and U is
!  upper triangular.  The factored form of a is then used to solve the
!  system of equations a * x = b.

! The following is from LAPACK documentation
! DGESV computes the solution to system of linear equations A * X = B for GE matrices
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Download DGESV + dependencies 
!  <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesv.f"> 
!  [TGZ]</a> 
!  <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesv.f"> 
!  [ZIP]</a> 
!  <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesv.f"> 
!  [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!       ..
!  
!
!  =============
! 
! 
!  DGESV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
! 
!  The LU decomposition with partial pivoting and row interchanges is
!  used to factor A as
!     A = P * L * U,
!  where P is a permutation matrix, L is unit lower triangular, and U is
!  upper triangular.  The factored form of A is then used to solve the
!  system of equations A * X = B.
!
!  Arguments:
!  ==========
!
!  \param[in] N
!           N is INTEGER
!           The number of linear equations, i.e., the order of the
!           matrix A.  N >= 0.
! 
!  \param[in] NRHS
!           NRHS is INTEGER
!           The number of right hand sides, i.e., the number of columns
!           of the matrix B.  NRHS >= 0.
! 
!  \param[in,out] A
!           A is DOUBLE PRECISION array, dimension (LDA,N)
!           On entry, the N-by-N coefficient matrix A.
!           On exit, the factors L and U from the factorization
!           A = P*L*U; the unit diagonal elements of L are not stored.
! 
!  \param[in] LDA
!           LDA is INTEGER
!           The leading dimension of the array A.  LDA >= max(1,N).
! 
!  \param[out] IPIV
!           IPIV is INTEGER array, dimension (N)
!           The pivot indices that define the permutation matrix P;
!           row i of the matrix was interchanged with row IPIV(i).
! 
!  \param[in,out] B
!           B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!           On entry, the N-by-NRHS matrix of right hand side matrix B.
!           On exit, if INFO = 0, the N-by-NRHS solution matrix X.
! 
!  \param[in] LDB
!           LDB is INTEGER
!           The leading dimension of the array B.  LDB >= max(1,N).
! 
!  \param[out] INFO
!           INFO is INTEGER
!           = 0:  successful exit
!           < 0:  if INFO = -i, the i-th argument had an illegal value
!           > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                 has been completed, but the factor U is exactly
!                 singular, so the solution could not be computed.
!
!  Authors:
!  ========
!
!  \author Univ. of Tennessee 
!  \author Univ. of California Berkeley 
!  \author Univ. of Colorado Denver 
!  \author NAG Ltd. 
!
!  \date November 2011
!
!  \ingroup doubleGEsolve

! !CALLED FROM:
! subroutine SoilTemperature in this module
!
! !REVISION HISTORY:
! 08/17/12 Keith Oleson: Initial code

!
! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_time_manager, only : get_step_size
    use clm_varcon      , only : rair, pstd, cpair, sb, hcv_roof, hcv_roof_enhanced, &
                                 hcv_floor, hcv_floor_enhanced, hcv_sunw, hcv_shdw, &
                                 em_roof_int, em_floor_int, em_sunw_int, em_shdw_int, &
                                 dz_floor, dens_floor, cp_floor, vent_ach
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varctl      , only : iulog
    use abortutils      , only : endrun
    use clm_varpar      , only : nlevurb, nlevsno, nlevgrnd
    use UrbanParamsType , only : urban_hac, urban_hac_off, urban_hac_on, urban_wasteheat_on
!
! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                   ! bounds
    integer , intent(in)  :: num_nolakec                      ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                ! column filter for non-lake points
    integer , intent(in)  :: num_urbanl                       ! number of urban landunits in clump
    integer , intent(in)  :: filter_urbanl(:)                 ! urban landunit filter
    real(r8), intent(in)  :: tk(bounds%begc: , -nlevsno+1: )  ! thermal conductivity (W m-1 K-1) [col, j]
    type(urbanparams_type), intent(in)    :: urbanparams_inst ! urban parameters
    type(temperature_type), intent(inout) :: temperature_inst ! temperature variables
    type(energyflux_type) , intent(inout) :: energyflux_inst  ! energy flux variables
!
! !LOCAL VARIABLES:
    integer, parameter :: neq = 5          ! number of equation/unknowns
    integer  :: fc,fl,c,l                  ! indices
    real(r8) :: dtime                      ! land model time step (s)
    real(r8) :: t_roof_inner_bef(bounds%begl:bounds%endl)  ! roof inside surface temperature at previous time step (K)              
    real(r8) :: t_sunw_inner_bef(bounds%begl:bounds%endl)  ! sunwall inside surface temperature at previous time step (K)              
    real(r8) :: t_shdw_inner_bef(bounds%begl:bounds%endl)  ! shadewall inside surface temperature at previous time step (K)              
    real(r8) :: t_floor_bef(bounds%begl:bounds%endl)       ! floor temperature at previous time step (K)              
    real(r8) :: t_building_bef(bounds%begl:bounds%endl)    ! internal building air temperature at previous time step [K]
    real(r8) :: t_building_bef_hac(bounds%begl:bounds%endl)! internal building air temperature before applying HAC [K]
    real(r8) :: hcv_roofi(bounds%begl:bounds%endl)         ! roof convective heat transfer coefficient (W m-2 K-1)
    real(r8) :: hcv_sunwi(bounds%begl:bounds%endl)         ! sunwall convective heat transfer coefficient (W m-2 K-1)
    real(r8) :: hcv_shdwi(bounds%begl:bounds%endl)         ! shadewall convective heat transfer coefficient (W m-2 K-1)
    real(r8) :: hcv_floori(bounds%begl:bounds%endl)        ! floor convective heat transfer coefficient (W m-2 K-1)
    real(r8) :: em_roofi(bounds%begl:bounds%endl)          ! roof inside surface emissivity (-)
    real(r8) :: em_sunwi(bounds%begl:bounds%endl)          ! sunwall inside surface emissivity (-)
    real(r8) :: em_shdwi(bounds%begl:bounds%endl)          ! shadewall inside surface emissivity (-)
    real(r8) :: em_floori(bounds%begl:bounds%endl)         ! floor inside surface emissivity (-)
    real(r8) :: dz_floori(bounds%begl:bounds%endl)         ! concrete floor thickness (m)
    real(r8) :: cp_floori(bounds%begl:bounds%endl)         ! concrete floor volumetric heat capacity (J m-3 K-1)
    real(r8) :: cv_floori(bounds%begl:bounds%endl)         ! intermediate calculation for concrete floor (W m-2 K-1)
    real(r8) :: rho_dair(bounds%begl:bounds%endl)          ! density of dry air at standard pressure and t_building (kg m-3)
    real(r8) :: vf_rf(bounds%begl:bounds%endl)             ! view factor of roof for floor (-)
    real(r8) :: vf_fr(bounds%begl:bounds%endl)             ! view factor of floor for roof (-)
    real(r8) :: vf_wf(bounds%begl:bounds%endl)             ! view factor of wall for floor (-)
    real(r8) :: vf_fw(bounds%begl:bounds%endl)             ! view factor of floor for wall (-)
    real(r8) :: vf_rw(bounds%begl:bounds%endl)             ! view factor of roof for wall (-)
    real(r8) :: vf_wr(bounds%begl:bounds%endl)             ! view factor of wall for roof (-)
    real(r8) :: vf_ww(bounds%begl:bounds%endl)             ! view factor of wall for wall (-)
    real(r8) :: zi_roof_innerl(bounds%begl:bounds%endl)    ! interface depth of nlevurb roof (m)
    real(r8) :: z_roof_innerl(bounds%begl:bounds%endl)     ! node depth of nlevurb roof (m)
    real(r8) :: zi_sunw_innerl(bounds%begl:bounds%endl)    ! interface depth of nlevurb sunwall (m)
    real(r8) :: z_sunw_innerl(bounds%begl:bounds%endl)     ! node depth of nlevurb sunwall (m)
    real(r8) :: zi_shdw_innerl(bounds%begl:bounds%endl)    ! interface depth of nlevurb shadewall (m)
    real(r8) :: z_shdw_innerl(bounds%begl:bounds%endl)     ! node depth of nlevurb shadewall (m)
    real(r8) :: t_roof_innerl_bef(bounds%begl:bounds%endl) ! roof temperature at nlevurb node depth at previous time step (K)
    real(r8) :: t_sunw_innerl_bef(bounds%begl:bounds%endl) ! sunwall temperature at nlevurb node depth at previous time step (K)
    real(r8) :: t_shdw_innerl_bef(bounds%begl:bounds%endl) ! shadewall temperature at nlevurb node depth at previous time step (K)
    real(r8) :: t_roof_innerl(bounds%begl:bounds%endl)     ! roof temperature at nlevurb node depth (K)
    real(r8) :: t_sunw_innerl(bounds%begl:bounds%endl)     ! sunwall temperature at nlevurb node depth (K)
    real(r8) :: t_shdw_innerl(bounds%begl:bounds%endl)     ! shadewall temperature at nlevurb node depth (K)
    real(r8) :: tk_roof_innerl(bounds%begl:bounds%endl)    ! roof thermal conductivity at nlevurb interface depth (W m-1 K-1)
    real(r8) :: tk_sunw_innerl(bounds%begl:bounds%endl)    ! sunwall thermal conductivity at nlevurb interface depth (W m-1 K-1)
    real(r8) :: tk_shdw_innerl(bounds%begl:bounds%endl)    ! shadewall thermal conductivity at nlevurb interface depth (W m-1 K-1)
    real(r8) :: qrd_roof(bounds%begl:bounds%endl)          ! roof inside net longwave for energy balance check (W m-2)
    real(r8) :: qrd_sunw(bounds%begl:bounds%endl)          ! sunwall inside net longwave for energy balance check (W m-2)
    real(r8) :: qrd_shdw(bounds%begl:bounds%endl)          ! shadewall inside net longwave for energy balance check (W m-2)
    real(r8) :: qrd_floor(bounds%begl:bounds%endl)         ! floor inside net longwave for energy balance check (W m-2)
    real(r8) :: qrd_building(bounds%begl:bounds%endl)      ! building inside net longwave for energy balance check (W m-2)
    real(r8) :: qcv_roof(bounds%begl:bounds%endl)          ! roof inside convection flux for energy balance check (W m-2)
    real(r8) :: qcv_sunw(bounds%begl:bounds%endl)          ! sunwall inside convection flux for energy balance check (W m-2)
    real(r8) :: qcv_shdw(bounds%begl:bounds%endl)          ! shadewall inside convection flux for energy balance check (W m-2)
    real(r8) :: qcv_floor(bounds%begl:bounds%endl)         ! floor inside convection flux for energy balance check (W m-2)
    real(r8) :: qcd_roof(bounds%begl:bounds%endl)          ! roof inside conduction flux for energy balance check (W m-2)
    real(r8) :: qcd_sunw(bounds%begl:bounds%endl)          ! sunwall inside conduction flux for energy balance check (W m-2)
    real(r8) :: qcd_shdw(bounds%begl:bounds%endl)          ! shadewall inside conduction flux for energy balance check (W m-2)
    real(r8) :: qcd_floor(bounds%begl:bounds%endl)         ! floor inside conduction flux for energy balance check (W m-2)
    real(r8) :: enrgy_bal_roof(bounds%begl:bounds%endl)    ! roof inside energy balance (W m-2)
    real(r8) :: enrgy_bal_sunw(bounds%begl:bounds%endl)    ! sunwall inside energy balance (W m-2)
    real(r8) :: enrgy_bal_shdw(bounds%begl:bounds%endl)    ! shadewall inside energy balance (W m-2)
    real(r8) :: enrgy_bal_floor(bounds%begl:bounds%endl)   ! floor inside energy balance (W m-2)
    real(r8) :: enrgy_bal_buildair(bounds%begl:bounds%endl)! building air energy balance (W m-2)
    real(r8) :: sum                        ! sum of view factors for floor, wall, roof
    integer  :: n                          ! number of linear equations (= neq)
    integer  :: nrhs                       ! number of right hand sides (= 1) 
    real(r8) :: a(neq,neq)                 ! n-by-n coefficient matrix a
    integer  :: lda                        ! leading dimension of the matrix a
    integer  :: ldb                        ! leading dimension of the matrix b
    real(r8) :: result(neq)                ! on entry, the right hand side of matrix b
                                           ! on exit, if info = 0, the n-by-nrhs solution matrix x
    integer  :: info                       ! exit information for LAPACK routine dgesv
    integer  :: ipiv(neq)                  ! the pivot indices that define the permutation matrix P
!EOP
!-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    associate(&
    clandunit         => col%landunit                      , & ! Input:  [integer (:)]  column's landunit
    ctype             => col%itype                         , & ! Input:  [integer (:)]  column type
    zi                => col%zi                            , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
    z                 => col%z                             , & ! Input:  [real(r8) (:,:)]  layer thickness (m)

    ht_roof           => lun%ht_roof                       , & ! Input:  [real(r8) (:)]  height of urban roof (m) 
    canyon_hwr        => lun%canyon_hwr                    , & ! Input:  [real(r8) (:)]  ratio of building height to street hwidth (-)
    wtlunit_roof      => lun%wtlunit_roof                  , & ! Input:  [real(r8) (:)]  weight of roof with respect to landunit
    urbpoi            => lun%urbpoi                        , & ! Input:  [logical (:)]  true => landunit is an urban point

    taf               => temperature_inst%taf_lun          , & ! Input:  [real(r8) (:)]  urban canopy air temperature (K)
    tssbef            => temperature_inst%t_ssbef_col      , & ! Input:  [real(r8) (:,:)]  temperature at previous time step (K)
    t_soisno          => temperature_inst%t_soisno_col     , & ! Input:  [real(r8) (:,:)]  soil temperature (K)
    t_roof_inner      => temperature_inst%t_roof_inner_lun , & ! InOut:  [real(r8) (:)]  roof inside surface temperature (K)
    t_sunw_inner      => temperature_inst%t_sunw_inner_lun , & ! InOut:  [real(r8) (:)]  sunwall inside surface temperature (K)
    t_shdw_inner      => temperature_inst%t_shdw_inner_lun , & ! InOut:  [real(r8) (:)]  shadewall inside surface temperature (K)
    t_floor           => temperature_inst%t_floor_lun      , & ! InOut:  [real(r8) (:)]  floor temperature (K)
    t_building        => temperature_inst%t_building_lun   , & ! InOut:  [real(r8) (:)]  internal building air temperature (K)

    t_building_max    => urbanparams_inst%t_building_max   , & ! Input:  [real(r8) (:)]  maximum internal building air temperature (K)
    t_building_min    => urbanparams_inst%t_building_min   , & ! Input:  [real(r8) (:)]  minimum internal building air temperature (K)

    eflx_building     => energyflux_inst%eflx_building_lun , & ! Output:  [real(r8) (:)]  building heat flux from change in interior building air temperature (W/m**2)
    eflx_urban_ac     => energyflux_inst%eflx_urban_ac_lun , & ! Output:  [real(r8) (:)]  urban air conditioning flux (W/m**2)
    eflx_urban_heat   => energyflux_inst%eflx_urban_heat_lun & ! Output:  [real(r8) (:)]  urban heating flux (W/m**2)
    )

    ! Get step size

    dtime = get_step_size()

    ! 1. Save t_* at previous time step
    ! 2. Set convective heat transfer coefficients (Bueno et al. 2012, GMD).
    !    An alternative is Salamanca et al. 2010, TAC, where they are all set to 8 W m-2 K-1.
    !    See clm_varcon.F90
    ! 3. Set inner surface emissivities (Bueno et al. 2012, GMD).
    ! 4. Set concrete floor properties (Salamanca et al. 2010, TAC).
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then
         t_roof_inner_bef(l)  = t_roof_inner(l)
         t_sunw_inner_bef(l)  = t_sunw_inner(l)
         t_shdw_inner_bef(l)  = t_shdw_inner(l)
         t_floor_bef(l)       = t_floor(l)
         t_building_bef(l)    = t_building(l)
         if (t_roof_inner_bef(l) .le. t_building_bef(l)) then
           hcv_roofi(l) = hcv_roof_enhanced
         else
           hcv_roofi(l) = hcv_roof
         end if
         if (t_floor_bef(l) .ge. t_building_bef(l)) then
           hcv_floori(l) = hcv_floor_enhanced
         else
           hcv_floori(l) = hcv_floor
         end if
         hcv_sunwi(l) = hcv_sunw
         hcv_shdwi(l) = hcv_shdw
         em_roofi(l)  = em_roof_int
         em_sunwi(l)  = em_sunw_int
         em_shdwi(l)  = em_shdw_int
         em_floori(l) = em_floor_int
         ! Concrete floor thickness (m)
         dz_floori(l) = dz_floor
         ! Concrete floor volumetric heat capacity (J m-3 K-1)
         cp_floori(l) = cp_floor
         ! Intermediate calculation for concrete floor (W m-2 K-1)
         cv_floori(l) = (dz_floori(l) * cp_floori(l)) / dtime
         ! Density of dry air at standard pressure and t_building (kg m-3)
         rho_dair(l) = pstd / (rair*t_building_bef(l))
       end if
    end do

    ! Get terms from soil temperature equations to compute conduction flux
    ! Negative is toward surface - heat added
    ! Note that the conduction flux here is in W m-2 wall area but for purposes of solving the set of
    ! simultaneous equations this must be converted to W m-2 ground area. This is done below when 
    ! setting up the equation coefficients.

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (urbpoi(l)) then
         if (ctype(c) == icol_roof) then
           zi_roof_innerl(l) = zi(c,nlevurb)
           z_roof_innerl(l) = z(c,nlevurb)
           t_roof_innerl_bef(l) = tssbef(c,nlevurb)
           t_roof_innerl(l) = t_soisno(c,nlevurb)
           tk_roof_innerl(l) = tk(c,nlevurb)
         else if (ctype(c) == icol_sunwall) then
           zi_sunw_innerl(l) = zi(c,nlevurb)
           z_sunw_innerl(l) = z(c,nlevurb)
           t_sunw_innerl_bef(l) = tssbef(c,nlevurb)
           t_sunw_innerl(l) = t_soisno(c,nlevurb)
           tk_sunw_innerl(l) = tk(c,nlevurb)
         else if (ctype(c) == icol_shadewall) then
           zi_shdw_innerl(l) = zi(c,nlevurb)
           z_shdw_innerl(l) = z(c,nlevurb)
           t_shdw_innerl_bef(l) = tssbef(c,nlevurb)
           t_shdw_innerl(l) = t_soisno(c,nlevurb)
           tk_shdw_innerl(l) = tk(c,nlevurb)
         end if
       end if
    end do

    ! Calculate view factors
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then

         vf_rf(l) = sqrt(1._r8 + canyon_hwr(l)**2._r8) - canyon_hwr(l)
         vf_fr(l) = vf_rf(l)

         ! This view factor implicitly converts from per unit wall area to per unit floor area
         vf_wf(l)  = 0.5_r8*(1._r8 - vf_rf(l))

         ! This view factor implicitly converts from per unit floor area to per unit wall area
         vf_fw(l) = vf_wf(l) / canyon_hwr(l)

         ! This view factor implicitly converts from per unit roof area to per unit wall area
         vf_rw(l)  = vf_fw(l)

         ! This view factor implicitly converts from per unit wall area to per unit roof area
         vf_wr(l)  = vf_wf(l)

         vf_ww(l)  = 1._r8 - vf_rw(l) - vf_fw(l)

      end if
    end do

    ! error check -- make sure view factor sums to one for floor, wall, and roof

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then

         sum = vf_rf(l) + 2._r8*vf_wf(l)
         if (abs(sum-1._r8) > 1.e-06_r8 ) then
            write (iulog,*) 'urban floor view factor error',sum
            write (iulog,*) 'clm model is stopping'
            call endrun()
         endif
         sum = vf_rw(l) + vf_fw(l) + vf_ww(l)
         if (abs(sum-1._r8) > 1.e-06_r8 ) then
            write (iulog,*) 'urban wall view factor error',sum
            write (iulog,*) 'clm model is stopping'
            call endrun()
         endif
         sum = vf_fr(l) + vf_wr(l) + vf_wr(l)
         if (abs(sum-1._r8) > 1.e-06_r8 ) then
            write (iulog,*) 'urban roof view factor error',sum
            write (iulog,*) 'clm model is stopping'
            call endrun()
         endif

       endif
    end do

    n = neq
    nrhs = 1
    lda = neq
    ldb = neq

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then

         ! ROOF
         a(1,1) =   0.5_r8*hcv_roofi(l) &
                  + 0.5_r8*tk_roof_innerl(l)/(zi_roof_innerl(l) - z_roof_innerl(l)) &
                  + 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8 &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rf(l)*(1._r8-em_floori(l))*vf_fr(l)

         a(1,2) = - 4._r8*em_roofi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wr(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l)

         a(1,3) = - 4._r8*em_roofi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wr(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l)

         a(1,4) = - 4._r8*em_roofi(l)*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fr(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wr(l)

         a(1,5) = - 0.5_r8*hcv_roofi(l)

         result(1) =   0.5_r8*tk_roof_innerl(l)*t_roof_innerl(l)/(zi_roof_innerl(l) - z_roof_innerl(l)) &
                     - 0.5_r8*tk_roof_innerl(l)*(t_roof_inner_bef(l)-t_roof_innerl_bef(l))/(zi_roof_innerl(l) &
                     - z_roof_innerl(l)) &
                     - 3._r8*em_roofi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wr(l) &
                     - 3._r8*em_roofi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wr(l) &
                     - 3._r8*em_roofi(l)*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fr(l) &
                     + 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8 &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rf(l)*(1._r8-em_floori(l))*vf_fr(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                     - 0.5_r8*hcv_roofi(l)*(t_roof_inner_bef(l) - t_building_bef(l))

         ! SUNWALL
         a(2,1) = - 4._r8*em_sunwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l)

         a(2,2) =   0.5_r8*hcv_sunwi(l)*canyon_hwr(l) &
                  + 0.5_r8*tk_sunw_innerl(l)/(zi_sunw_innerl(l) - z_sunw_innerl(l))*canyon_hwr(l) &
                  + 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8 &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l)

         a(2,3) = - 4._r8*em_sunwi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_ww(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l)

         a(2,4) = - 4._r8*em_sunwi(l)*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_shdwi(l))*vf_ww(l)
         a(2,5) = - 0.5_r8*hcv_sunwi(l)*canyon_hwr(l)

         result(2) =   0.5_r8*tk_sunw_innerl(l)*t_sunw_innerl(l)/(zi_sunw_innerl(l) - z_sunw_innerl(l))*canyon_hwr(l) &
                     - 0.5_r8*tk_sunw_innerl(l)*(t_sunw_inner_bef(l)-t_sunw_innerl_bef(l))/(zi_sunw_innerl(l) &
                     - z_sunw_innerl(l))*canyon_hwr(l) &
                     - 3._r8*em_sunwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l) &
                     - 3._r8*em_sunwi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_ww(l) &
                     - 3._r8*em_sunwi(l)*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l) &
                     + 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8 &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_ww(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                     - 0.5_r8*hcv_sunwi(l)*(t_sunw_inner_bef(l) - t_building_bef(l))*canyon_hwr(l)

         ! SHADEWALL
         a(3,1) = - 4._r8*em_shdwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l)

         a(3,2) = - 4._r8*em_shdwi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_ww(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l)

         a(3,3) =   0.5_r8*hcv_shdwi(l)*canyon_hwr(l) &
                  + 0.5_r8*tk_shdw_innerl(l)/(zi_shdw_innerl(l) - z_shdw_innerl(l))*canyon_hwr(l) &
                  + 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8 &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l)

         a(3,4) = - 4._r8*em_shdwi(l)*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_sunwi(l))*vf_ww(l)

         a(3,5) = - 0.5_r8*hcv_shdwi(l)*canyon_hwr(l)

         result(3) =   0.5_r8*tk_shdw_innerl(l)*t_shdw_innerl(l)/(zi_shdw_innerl(l) - z_shdw_innerl(l))*canyon_hwr(l) &
                     - 0.5_r8*tk_shdw_innerl(l)*(t_shdw_inner_bef(l)-t_shdw_innerl_bef(l))/(zi_shdw_innerl(l) &
                     - z_shdw_innerl(l))*canyon_hwr(l) &
                     - 3._r8*em_shdwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l) &
                     - 3._r8*em_shdwi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_ww(l) &
                     - 3._r8*em_shdwi(l)*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l) &
                     + 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8 &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_ww(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                     - 0.5_r8*hcv_shdwi(l)*(t_shdw_inner_bef(l) - t_building_bef(l))*canyon_hwr(l)

         ! FLOOR
         a(4,1) = - 4._r8*em_floori(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rf(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                  - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wf(l)

         a(4,2) = - 4._r8*em_floori(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wf(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wf(l) &
                  - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l)

         a(4,3) = - 4._r8*em_floori(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wf(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                  - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wf(l)

         a(4,4) =   (cv_floori(l) + 0.5_r8*hcv_floori(l)) &
                  + 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8 &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                  - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wf(l)

         a(4,5) = - 0.5_r8*hcv_floori(l)

         result(4) =   cv_floori(l)*t_floor_bef(l) &
                     - 3._r8*em_floori(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rf(l) &
                     - 3._r8*em_floori(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wf(l) &
                     - 3._r8*em_floori(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wf(l) &
                     + 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8 &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                     - 3._r8*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wf(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wf(l) &
                     - 3._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                     - 3._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                     - 3._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wf(l) &
                     - 0.5_r8*hcv_floori(l)*(t_floor_bef(l) - t_building_bef(l))

         ! Building air temperature
         a(5,1) = - 0.5_r8*hcv_roofi(l)
         a(5,2) = - 0.5_r8*hcv_sunwi(l)*canyon_hwr(l)

         a(5,3) = - 0.5_r8*hcv_shdwi(l)*canyon_hwr(l)

         a(5,4) = - 0.5_r8*hcv_floori(l)

         a(5,5) =  ((ht_roof(l)*rho_dair(l)*cpair)/dtime) + &
                   ((ht_roof(l)*vent_ach)/3600._r8)*rho_dair(l)*cpair + &
                   0.5_r8*hcv_roofi(l) + &
                   0.5_r8*hcv_sunwi(l)*canyon_hwr(l) + &
                   0.5_r8*hcv_shdwi(l)*canyon_hwr(l) + &
                   0.5_r8*hcv_floori(l)

         result(5) = (ht_roof(l)*rho_dair(l)*cpair/dtime)*t_building_bef(l) &
                      + ((ht_roof(l)*vent_ach)/3600._r8)*rho_dair(l)*cpair*taf(l) &
                      + 0.5_r8*hcv_roofi(l)*(t_roof_inner_bef(l) - t_building_bef(l)) &
                      + 0.5_r8*hcv_sunwi(l)*(t_sunw_inner_bef(l) - t_building_bef(l))*canyon_hwr(l) &
                      + 0.5_r8*hcv_shdwi(l)*(t_shdw_inner_bef(l) - t_building_bef(l))*canyon_hwr(l) &
                      + 0.5_r8*hcv_floori(l)*(t_floor_bef(l) - t_building_bef(l))

         ! Solve equations
         call dgesv(n, nrhs, a, lda, ipiv, result, ldb, info)

         ! If dgesv fails, abort 
         if (info /= 0) then
           write(iulog,*)'fl: ',fl
           write(iulog,*)'l: ',l
           write(iulog,*)'dgesv info: ',info
           write (iulog,*) 'dgesv error'
           write (iulog,*) 'clm model is stopping'
           call endrun()
         end if
         ! Assign new temperatures
         t_roof_inner(l)  = result(1)
         t_sunw_inner(l)  = result(2)
         t_shdw_inner(l)  = result(3)
         t_floor(l) = result(4)
         t_building(l)    = result(5)
       end if
    end do

    ! Energy balance checks
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then
         qrd_roof(l) = - em_roofi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wr(l) &
                       - 4._r8*em_roofi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wr(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - em_roofi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wr(l) &
                       - 4._r8*em_roofi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wr(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - em_roofi(l)*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fr(l) &
                       - 4._r8*em_roofi(l)*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fr(l)*(t_floor(l) - t_floor_bef(l)) &
                       - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                       - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wr(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                       - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wr(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rf(l)*(1._r8-em_floori(l))*vf_fr(l) &
                       - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rf(l)*(1._r8-em_floori(l))*vf_fr(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                       - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wr(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l) &
                       - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                       - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wr(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l) &
                       - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wf(l)*(1._r8-em_floori(l))*vf_fr(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wr(l) &
                       - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wr(l)*(t_floor(l) &
                       - t_floor_bef(l)) &
                       - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wr(l) &
                       - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wr(l)*(t_floor(l) &
                       - t_floor_bef(l)) &
                       + em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8 &
                       + 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*(t_roof_inner(l) - t_roof_inner_bef(l)) 

         qrd_sunw(l) = - em_sunwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l) &
                       - 4._r8*em_sunwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - em_sunwi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_ww(l)  &
                       - 4._r8*em_sunwi(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_ww(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - em_sunwi(l)*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l) &
                       - 4._r8*em_sunwi(l)*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(t_floor(l) - t_floor_bef(l)) &
                       - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                       - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3.*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_ww(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                       - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3.*vf_ww(l)*(1._r8-em_shdwi(l))*vf_ww(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                       - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3.*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                       - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3.*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                       - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3.*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rw(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                       - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3.*vf_rw(l)*(1._r8-em_shdwi(l))*vf_ww(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                       - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3.*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                       - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l)*(t_floor(l) &
                       - t_floor_bef(l)) &
                       - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fw(l)*(1._r8-em_shdwi(l))*vf_ww(l) &
                       - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*vf_fw(l)*(1._r8-em_shdwi(l))*vf_ww(l)*(t_floor(l) &
                       - t_floor_bef(l)) &
                       + em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8 &
                       + 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*(t_sunw_inner(l) - t_sunw_inner_bef(l))

         qrd_shdw(l) = - em_shdwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rw(l) &
                       - 4._r8*em_shdwi(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rw(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - em_shdwi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_ww(l) &
                       - 4._r8*em_shdwi(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_ww(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - em_shdwi(l)*em_floori(l)*sb*t_floor_bef(l)**4._r8*vf_fw(l) &
                       - 4._r8*em_shdwi(l)*em_floori(l)*sb*t_floor_bef(l)**3._r8*vf_fw(l)*(t_floor(l) - t_floor_bef(l)) &
                       - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                       - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3.*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_ww(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                       - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3.*vf_ww(l)*(1._r8-em_sunwi(l))*vf_ww(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                       - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3.*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l)*(t_shdw_inner(l) &
                       - t_shdw_inner_bef(l)) &
                       - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                       - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3.*vf_wf(l)*(1._r8-em_floori(l))*vf_fw(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                       - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3.*vf_wr(l)*(1._r8-em_roofi(l))*vf_rw(l)*(t_sunw_inner(l) &
                       - t_sunw_inner_bef(l)) &
                       - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rw(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                       - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3.*vf_rw(l)*(1._r8-em_sunwi(l))*vf_ww(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l) &
                       - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3.*vf_rf(l)*(1._r8-em_floori(l))*vf_fw(l)*(t_roof_inner(l) &
                       - t_roof_inner_bef(l)) &
                       - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l) &
                       - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*vf_fr(l)*(1._r8-em_roofi(l))*vf_rw(l)*(t_floor(l) &
                       - t_floor_bef(l)) &
                       - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fw(l)*(1._r8-em_sunwi(l))*vf_ww(l) &
                       - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*vf_fw(l)*(1._r8-em_sunwi(l))*vf_ww(l)*(t_floor(l) &
                       - t_floor_bef(l)) &
                       + em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8 &
                       + 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*(t_shdw_inner(l) - t_shdw_inner_bef(l))

         qrd_floor(l) = - em_floori(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8*vf_rf(l) &
                        - 4._r8*em_floori(l)*em_roofi(l)*sb*t_roof_inner_bef(l)**3._r8*vf_rf(l)*(t_roof_inner(l) &
                        - t_roof_inner_bef(l)) &
                        - em_floori(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8*vf_wf(l) &
                        - 4._r8*em_floori(l)*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3._r8*vf_wf(l)*(t_sunw_inner(l) &
                        - t_sunw_inner_bef(l)) &
                        - em_floori(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8*vf_wf(l) &
                        - 4._r8*em_floori(l)*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3._r8*vf_wf(l)*(t_shdw_inner(l) &
                        - t_shdw_inner_bef(l)) &
                        - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                        - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*vf_fr(l)*(1._r8-em_roofi(l))*vf_rf(l)*(t_floor(l) &
                        - t_floor_bef(l)) &
                        - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                        - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*vf_fw(l)*(1._r8-em_sunwi(l))*vf_wf(l)*(t_floor(l) &
                        - t_floor_bef(l)) &
                        - (em_floori(l)*sb*t_floor_bef(l)**4._r8)*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wf(l) &
                        - 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*vf_fw(l)*(1._r8-em_shdwi(l))*vf_wf(l)*(t_floor(l) &
                        - t_floor_bef(l)) &
                        - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wf(l) &
                        - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3.*vf_ww(l)*(1._r8-em_shdwi(l))*vf_wf(l)*(t_sunw_inner(l) &
                        - t_sunw_inner_bef(l)) &
                        - (em_sunwi(l)*sb*t_sunw_inner_bef(l)**4._r8)*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                        - 4._r8*em_sunwi(l)*sb*t_sunw_inner_bef(l)**3.*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l)*(t_sunw_inner(l) &
                        - t_sunw_inner_bef(l)) &
                        - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l) &
                        - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3.*vf_wr(l)*(1._r8-em_roofi(l))*vf_rf(l)*(t_shdw_inner(l) &
                        - t_shdw_inner_bef(l)) &
                        - (em_shdwi(l)*sb*t_shdw_inner_bef(l)**4._r8)*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                        - 4._r8*em_shdwi(l)*sb*t_shdw_inner_bef(l)**3.*vf_ww(l)*(1._r8-em_sunwi(l))*vf_wf(l)*(t_shdw_inner(l) &
                        - t_shdw_inner_bef(l)) &
                        - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wf(l) &
                        - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3.*vf_rw(l)*(1._r8-em_sunwi(l))*vf_wf(l)*(t_roof_inner(l) &
                        - t_roof_inner_bef(l)) &
                        - (em_roofi(l)*sb*t_roof_inner_bef(l)**4._r8)*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wf(l) &
                        - 4._r8*em_roofi(l)*sb*t_roof_inner_bef(l)**3.*vf_rw(l)*(1._r8-em_shdwi(l))*vf_wf(l)*(t_roof_inner(l) &
                        - t_roof_inner_bef(l)) &
                        + em_floori(l)*sb*t_floor_bef(l)**4._r8 &
                        + 4._r8*em_floori(l)*sb*t_floor_bef(l)**3.*(t_floor(l) - t_floor_bef(l))

         qrd_building(l) = qrd_roof(l) + canyon_hwr(l)*(qrd_sunw(l) + qrd_shdw(l)) + qrd_floor(l)

         if (abs(qrd_building(l)) > .10_r8 ) then
           write (iulog,*) 'urban inside building net longwave radiation balance error ',qrd_building(l)
           write (iulog,*) 'clm model is stopping'
           call endrun()
         end if

         qcv_roof(l) = 0.5_r8*hcv_roofi(l)*(t_roof_inner(l) - t_building(l)) + 0.5_r8*hcv_roofi(l)*(t_roof_inner_bef(l) &
                       - t_building_bef(l))
         qcd_roof(l) = 0.5_r8*tk_roof_innerl(l)*(t_roof_inner(l) - t_roof_innerl(l))/(zi_roof_innerl(l) - z_roof_innerl(l))  &
                       + 0.5_r8*tk_roof_innerl(l)*(t_roof_inner_bef(l) - t_roof_innerl_bef(l))/(zi_roof_innerl(l) &
                       - z_roof_innerl(l))
         enrgy_bal_roof(l) = qrd_roof(l) + qcv_roof(l) + qcd_roof(l)
         if (abs(enrgy_bal_roof(l)) > .10_r8 ) then
           write (iulog,*) 'urban inside roof energy balance error ',enrgy_bal_roof(l)
           write (iulog,*) 'clm model is stopping'
           call endrun()
         end if

         qcv_sunw(l) = 0.5_r8*hcv_sunwi(l)*(t_sunw_inner(l) - t_building(l)) + 0.5_r8*hcv_sunwi(l)*(t_sunw_inner_bef(l) &
                       - t_building_bef(l))
         qcd_sunw(l) = 0.5_r8*tk_sunw_innerl(l)*(t_sunw_inner(l) - t_sunw_innerl(l))/(zi_sunw_innerl(l) - z_sunw_innerl(l))  &
                       + 0.5_r8*tk_sunw_innerl(l)*(t_sunw_inner_bef(l) - t_sunw_innerl_bef(l))/(zi_sunw_innerl(l) &
                       - z_sunw_innerl(l))
         enrgy_bal_sunw(l) = qrd_sunw(l) + qcv_sunw(l)*canyon_hwr(l) + qcd_sunw(l)*canyon_hwr(l)
         if (abs(enrgy_bal_sunw(l)) > .10_r8 ) then
           write (iulog,*) 'urban inside sunwall energy balance error ',enrgy_bal_sunw(l)
           write (iulog,*) 'clm model is stopping'
           call endrun()
         end if

         qcv_shdw(l) = 0.5_r8*hcv_shdwi(l)*(t_shdw_inner(l) - t_building(l)) + 0.5_r8*hcv_shdwi(l)*(t_shdw_inner_bef(l) &
                       - t_building_bef(l))
         qcd_shdw(l) = 0.5_r8*tk_shdw_innerl(l)*(t_shdw_inner(l) - t_shdw_innerl(l))/(zi_shdw_innerl(l) - z_shdw_innerl(l))  &
                       + 0.5_r8*tk_shdw_innerl(l)*(t_shdw_inner_bef(l) - t_shdw_innerl_bef(l))/(zi_shdw_innerl(l) &
                       - z_shdw_innerl(l))
         enrgy_bal_shdw(l) = qrd_shdw(l) + qcv_shdw(l)*canyon_hwr(l) + qcd_shdw(l)*canyon_hwr(l)
         if (abs(enrgy_bal_shdw(l)) > .10_r8 ) then
           write (iulog,*) 'urban inside shadewall energy balance error ',enrgy_bal_shdw(l)
           write (iulog,*) 'clm model is stopping'
           call endrun()
         end if

         qcv_floor(l) = 0.5_r8*hcv_floori(l)*(t_floor(l) - t_building(l)) + 0.5_r8*hcv_floori(l)*(t_floor_bef(l) &
                        - t_building_bef(l))
         qcd_floor(l) = cv_floori(l)*(t_floor(l) - t_floor_bef(l))
         enrgy_bal_floor(l) = qrd_floor(l) + qcv_floor(l) + qcd_floor(l)
         if (abs(enrgy_bal_floor(l)) > .10_r8 ) then
           write (iulog,*) 'urban inside floor energy balance error ',enrgy_bal_floor(l)
           write (iulog,*) 'clm model is stopping'
           call endrun()
         end if

         enrgy_bal_buildair(l) = (ht_roof(l)*rho_dair(l)*cpair/dtime)*(t_building(l) - t_building_bef(l)) &
                                 - ht_roof(l)*(vent_ach/3600._r8)*rho_dair(l)*cpair*(taf(l) - t_building(l)) &
                                 - 0.5_r8*hcv_roofi(l)*(t_roof_inner(l) - t_building(l)) &
                                 - 0.5_r8*hcv_roofi(l)*(t_roof_inner_bef(l) - t_building_bef(l)) &
                                 - 0.5_r8*hcv_sunwi(l)*(t_sunw_inner(l) - t_building(l))*canyon_hwr(l) &
                                 - 0.5_r8*hcv_sunwi(l)*(t_sunw_inner_bef(l) - t_building_bef(l))*canyon_hwr(l) &
                                 - 0.5_r8*hcv_shdwi(l)*(t_shdw_inner(l) - t_building(l))*canyon_hwr(l) &
                                 - 0.5_r8*hcv_shdwi(l)*(t_shdw_inner_bef(l) - t_building_bef(l))*canyon_hwr(l) &
                                 - 0.5_r8*hcv_floori(l)*(t_floor(l) - t_building(l)) &
                                 - 0.5_r8*hcv_floori(l)*(t_floor_bef(l) - t_building_bef(l))
         if (abs(enrgy_bal_buildair(l)) > .10_r8 ) then
           write (iulog,*) 'urban building air energy balance error ',enrgy_bal_buildair(l)
           write (iulog,*) 'clm model is stopping'
           call endrun()
         end if
       end if
    end do

    ! Restrict internal building air temperature to between min and max
    ! Calculate heating or air conditioning flux from energy required to change
    ! internal building air temperature to t_building_min or t_building_max. 

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then
          if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
            t_building_bef_hac(l) = t_building(l)
!           rho_dair(l) = pstd / (rair*t_building(l))

            if (t_building_bef_hac(l) > t_building_max(l)) then
              t_building(l) = t_building_max(l)
              eflx_urban_ac(l) = wtlunit_roof(l) * abs( (ht_roof(l) * rho_dair(l) * cpair / dtime) * t_building(l) &
                                 - (ht_roof(l) * rho_dair(l) * cpair / dtime) * t_building_bef_hac(l) )
            else if (t_building_bef_hac(l) < t_building_min(l)) then
              t_building(l) = t_building_min(l)
              eflx_urban_heat(l) = wtlunit_roof(l) * abs( (ht_roof(l) * rho_dair(l) * cpair / dtime) * t_building(l) &
                                   - (ht_roof(l) * rho_dair(l) * cpair / dtime) * t_building_bef_hac(l) )
            else
              eflx_urban_ac(l) = 0._r8
              eflx_urban_heat(l) = 0._r8
            end if
          else
            eflx_urban_ac(l) = 0._r8
            eflx_urban_heat(l) = 0._r8
          end if
          eflx_building(l) = wtlunit_roof(l) * (ht_roof(l) * rho_dair(l)*cpair/dtime) * (t_building(l) - t_building_bef(l))
       end if
    end do

    end associate 
  end subroutine BuildingTemperature

  !-----------------------------------------------------------------------

end module UrbBuildTempOleson2015Mod
