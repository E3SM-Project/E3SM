	  MODULE rrtm_vars
!-----------------------------------------------------------------
!  This module specifies input and output variables for the
!  RRTMG radiation code.
!-----------------------------------------------------------------
	  USE parkind,   only: kind_rb, kind_im
	  USE rrtm_grid, only: nx, ny, nz, nzm

 	  IMPLICIT NONE
      SAVE

! RRTMG input variables -------------------------------------------------------

      REAL (KIND=kind_rb), allocatable, DIMENSION(:,:,:) :: &
          tabs, &    ! layer temperature (K)
          qv, &      ! layer vapor mixing ratio (g/g)
          qcl, &     ! layer cloud water mixing ratio (g/g)
          qci        ! layer cloud ice mixing ratio (g/g)
          
      REAL (KIND=kind_rb), allocatable, DIMENSION(:,:,:) :: &       ! JUNG_LocalP
          pres_loc,  &   ! model layer local pressure (hPa)
          presi_loc      ! model interface local pressure (hPa)

      REAL (KIND=kind_rb), allocatable, DIMENSION(:,:) :: &
          sstxy    ! sea surface temperature (K)

      REAL (KIND=kind_rb), DIMENSION(nzm) :: pres    ! model layer pressure (hPa)
      REAL (KIND=kind_rb), DIMENSION(nz)  :: presi   ! model interface pressure (hPa)

! RRTMG output variables ------------------------------------------------------

! Domain-averaged diagnostic fields
      REAL (KIND=kind_rb), DIMENSION(nz) :: &
          radlwup, &    ! LW up flux (average of lwUp)
          radlwdn, &    ! LW down flux (average of lwDown)
          radswup, &    ! SW up flux (average of swUp)
          radswdn       ! SW down flux (average of swDown)

! Domain-averaged diagnostic fields
      REAL (KIND=kind_rb), DIMENSION(nzm) :: &
          radqrlw, &    ! LW heating rate (K/day, average of lwHeatingRate)
          radqrsw, &    ! SW heating rate (K/day, average of swHeatingRate)
          radqrclw, &   ! LW clear sky heating rate (K/day)
          radqrcsw      ! SW clear sky heating rate (K/day)

! 2D diagnostics
      REAL (KIND=kind_rb), allocatable, DIMENSION(:,:) :: &
          lwns_xy, &                ! net surface LW flux (W/m2)
          lwnt_xy, &                ! net TOA LW flux (W/m2)
          swns_xy, &                ! net surface SW flux (W/m2)
          swnt_xy, &                ! net TOA SW flux (W/m2)
          solin_xy, &               ! TOA insolation (W/m2)
          lwnsc_xy, &               ! net surface LW flux, clear sky (W/m2)
          lwntc_xy, &               ! net TOA LW flux, clear sky (W/m2)
          swnsc_xy, &               ! net surface SW flux, clear sky (W/m2)
          swntc_xy, &               ! net TOA SW flux, clear sky (W/m2)
!          lwDownSurface, &          ! LW downward flux at surface (W/m2)
          lwDownSurfaceClearSky, &  ! LW downward flux at surface, clear sky (W/m2)
          lwUpSurface, &            ! LW upward flux at surface (W/m2)
          lwUpSurfaceClearSky, &    ! LW upward flux at surface, clear sky (W/m2)
          lwUpToa, &                ! LW upward flux at TOA (W/m2)
          lwUpToaClearSky, &        ! LW upward flux at TOA, clear sky (W/m2)
!          swDownSurface, &          ! SW downward flux at surface (W/m2)
          swDownSurfaceClearSky, &  ! SW downward flux at surface, clear sky (W/m2)
          swUpSurface, &            ! SW upward flux at surface (W/m2)
          swUpSurfaceClearSky, &    ! SW upward flux at surface, clear sky (W/m2)
          swDownToa, &              ! SW downward flux at TOA (W/m2)
          swUpToa, &                ! SW upward flux at TOA (W/m2)
          swUpToaClearSky           ! SW upward flux at TOA, clear sky (W/m2)
!          insolation_TOA            ! TOA insolation (W/m2)

! 1D diagnostics
      REAL (KIND=kind_rb) :: &
          s_flns, &     ! accumulated sum of lwns_xy
          s_flnt, &     ! accumulated sum of lwnt_xy
          s_fsns, &     ! accumulated sum of swns_xy
          s_fsnt, &     ! accumulated sum of swnt_xy
          s_flnsc, &    ! accumulated sum of lwnsc_xy
          s_flntc, &    ! accumulated sum of lwntc_xy
          s_fsnsc, &    ! accumulated sum of swnsc_xy
          s_fsntc, &    ! accumulated sum of swntc_xy
          s_solin, &    ! accumulated sum of solin_xy
          s_flds, &     ! accumulated sum of downwelling LW flux at surface
          s_fsds        ! accumulated sum of downwelling SW flux at surface

! Public data
!      REAL, (KIND=kind_rb), DIMENSION(nx, ny, nzm) :: &
!          qrad          ! Radiative heating rate (K/s), used for writing to restart files

!      REAL, (KIND=kind_rb), DIMENSION(nx, ny) :: &
!          lwnsxy, &   ! Net longwave flux (up + down) at surface (W/m2)
!          swnsxy      ! Net shortwave flux (up + down) at surface (W/m2)

!      REAL (KIND=kind_rb), DIMENSION(nx, nzm) :: &
!          lwHeatingRate, &          ! Longwave heating rate (K/day)
!          swHeatingRate, &          ! Shortwave heating rate (K/day)
!          swHeatingRateClearSky, &  ! Shortwave heating rate, clear sky (K/day)
!          lwHeatingRateClearSky     ! Longwave heating rate, clear sky (K/day)

!      REAL (KIND=kind_rb), DIMENSION(nx, nz) :: &
!          lwUp,   &          ! Longwave upward flux (W/m2)
!          lwDown, &          ! Longwave downward flux (W/m2)
!          swUp,   &          ! Shortwave upward flux (W/m2)
!          swDown, &          ! Shortwave downward flux (W/m2)
!          swUpClearSky, &    ! Shortwave upward flux, clear sky (W/m2)
!          swDownClearSky, &  ! Shortwave downward flux, clear sky (W/m2)
!          lwUpClearSky, &    ! Longwave upward flux, clear sky (W/m2)
!          lwDownClearSky     ! Longwave downward flux, clear sky (W/m2)

      public :: allocate_rrtm_vars

CONTAINS

  subroutine allocate_rrtm_vars()
  
    allocate(tabs(nx,ny,nzm))
    allocate(qv(nx,ny,nzm))
    allocate(qcl(nx,ny,nzm))
    allocate(qci(nx,ny,nzm))
    
    allocate(pres_loc(nx,ny,nzm))   ! LocalP
    allocate(presi_loc(nx,ny,nz))   ! LocalP
    
    allocate(sstxy(nx,ny))
    
    allocate(lwns_xy(nx,ny))
    allocate(lwnt_xy(nx,ny))
    allocate(swns_xy(nx,ny))
    allocate(swnt_xy(nx,ny))
    allocate(solin_xy(nx,ny))
    allocate(lwnsc_xy(nx,ny))
    allocate(lwntc_xy(nx,ny))
    allocate(swnsc_xy(nx,ny))
    allocate(swntc_xy(nx,ny))
    
    allocate(lwDownSurfaceClearSky(nx,ny))
    allocate(lwUpSurface(nx,ny))
    allocate(lwUpSurfaceClearSky(nx,ny))
    allocate(lwUpToa(nx,ny))
    allocate(lwUpToaClearSky(nx,ny))
    
    allocate(swDownSurfaceClearSky(nx,ny))
    allocate(swUpSurface(nx,ny))
    allocate(swUpSurfaceClearSky(nx,ny))
    allocate(swDownToa(nx,ny))
    allocate(swUpToa(nx,ny))
    allocate(swUpToaClearSky(nx,ny))
    
  end subroutine allocate_rrtm_vars

!=======================================================================
 	  END MODULE rrtm_vars
!=======================================================================
