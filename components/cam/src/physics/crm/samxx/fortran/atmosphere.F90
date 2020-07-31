
module atmosphere_mod
	implicit none

contains

  SUBROUTINE Atmosphere(alt, sigma, delta, theta)
    !   -------------------------------------------------------------------------
    ! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
    ! AUTHOR - Ralph Carmichael, Public Domain Aeronautical Software
    ! NOTE - If alt > 86, the values returned will not be correct, but they will
    !   not be too far removed from the correct values for density.
    !   The reference document does not use the terms pressure and temperature
    !   above 86 km.
    use params, only: crm_rknd
    IMPLICIT NONE
    !============================================================================
    !     A R G U M E N T S                                                     |
    !============================================================================
    REAL(crm_rknd),INTENT(IN)::  alt  ! geometric altitude, km.
    REAL(crm_rknd),INTENT(OUT):: sigma! density/sea-level standard density
    REAL(crm_rknd),INTENT(OUT):: delta! pressure/sea-level standard pressure
    REAL(crm_rknd),INTENT(OUT):: theta! temperature/sea-level standard temperature
    !============================================================================
    !     L O C A L   C O N S T A N T S                                         |
    !============================================================================
    REAL(crm_rknd),PARAMETER:: REARTH = 6369.0 ! radius of the Earth (km)
    REAL(crm_rknd),PARAMETER:: GMR = 34.163195 ! gas constant
    INTEGER,PARAMETER:: NTAB=8! number of entries in the defining tables
    !============================================================================
    !     L O C A L   V A R I A B L E S                                         |
    !============================================================================
    INTEGER:: i,j,k    ! counters
    REAL(crm_rknd):: h           ! geopotential altitude (km)
    REAL(crm_rknd):: tgrad, tbase! temperature gradient and base temp of this layer
    REAL(crm_rknd):: tlocal      ! local temperature
    REAL(crm_rknd):: deltah      ! height above base of this layer
    !============================================================================
    !     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
    !============================================================================
    REAL(crm_rknd),DIMENSION(NTAB),PARAMETER:: htab= (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0,84.852/)
    REAL(crm_rknd),DIMENSION(NTAB),PARAMETER:: ttab=  (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
    REAL(crm_rknd),DIMENSION(NTAB),PARAMETER:: ptab=  (/1.0, 2.233611e-1, &
    5.403295e-2, 8.5666784e-3, 1.0945601e-3, 6.6063531e-4, 3.9046834e-5, 3.68501e-6/)
    REAL(crm_rknd),DIMENSION(NTAB),PARAMETER:: gtab=  (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
    !----------------------------------------------------------------------------
    h=alt*REARTH/(alt+REARTH)! convert geometric to geopotential altitude

    i=1
    j=NTAB                     ! setting up for=binary search
    DO
      k=(i+j)/2
      IF (h < htab(k)) THEN
        j=k
      ELSE
        i=k
      END IF
      IF (j <= i+1) EXIT
    END DO

    tgrad=gtab(i)              ! i will be in 1...NTAB-1
    tbase=ttab(i)
    deltah=h-htab(i)
    tlocal=tbase+tgrad*deltah
    theta=tlocal/ttab(1)            ! temperature ratio

    IF (tgrad == 0.0) THEN        ! pressure ratio
      delta=ptab(i)*EXP(-GMR*deltah/tbase)
    ELSE
      delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
    END IF

    sigma=delta/theta            ! density ratio
    RETURN
  END Subroutine Atmosphere

end module atmosphere_mod
