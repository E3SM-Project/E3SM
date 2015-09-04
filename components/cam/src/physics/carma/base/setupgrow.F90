! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine defines time-independent parameters used to calculate
!! condensational growth/evaporation.
!!
!! The parameters defined for each gas are 
!1>
!!   gwtmol:   molecular weight [g/mol]
!!   diffus:   diffusivity      [cm^2/s]
!!   rlhe  :   latent heat of evaporation [cm^2/s^2]
!!   rlhm  :   latent heat of melting [cm^2/s^2]
!!<
!! Time-independent parameters that depend on particle radius are
!! defined in setupgkern.f.
!!
!! This routine requires that vertical profiles of temperature <T>,
!! and pressure <p> are defined.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine setupgrow(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc       !! return code, negative indicates failure
  
  ! Local Variable
  integer                        :: ielem    !! element index
  integer                        :: k        !! z index
  integer                        :: i
  real(kind=f)                   :: rhoa_cgs, aden
  ! Define formats
  1 format(a,':  ',12i6)
  2 format(a,':  ',i6)
  3 format(/' id  gwtmol   gasname',(/,i3,3x,f5.1,3x,a))
  5 format(/,'Particle growth mapping arrays (setupgrow):')


  !-----Check that values are valid------------------------------------------
  do ielem = 1, NELEM
    if( igrowgas(ielem) .gt. NGAS )then
      if (do_print) write(LUNOPRT,*) 'setupgrow::ERROR - component of igrowgas > NGAS'
      rc = -1
      return
    endif
  enddo

  ! Define parameters with weak time-dependence to be used in
  ! growth equation.
  do k = 1, NZ

    ! Diffusivity of water vapor in air from Pruppacher & Klett (eq. 13-3);
    ! units are [cm^2/s].
    if (igash2o /= 0) then 
      diffus(k, igash2o) = 0.211_f * (1.01325e+6_f / p(k)) * (t(k) / 273.15_f )**1.94_f
  
      ! Latent heat of evaporation for water; units are [cm^2/s^2]
      if (do_cnst_rlh) then
        rlhe(k, igash2o) = RLHE_CNST
      else
        ! from Stull
        rlhe(k, igash2o) = (2.5_f - .00239_f * (t(k) - 273.16_f)) * 1.e10_f      
      end if
  
      ! Latent heat of ice melting; units are [cm^2/s^2]
      if (do_cnst_rlh) then
        rlhm(k, igash2o) = RLHM_CNST
      else
  
        ! from Pruppacher & Klett (eq. 4-85b)
        !
        ! NOTE: This expression yields negative values for rlmh at mesospheric
        ! temperatures.
        rlhm(k, igash2o) = (79.7_f + 0.485_f * (t(k) - 273.16_f) - 2.5e-3_f * &
          ((t(k) - 273.16_f)**2)) * 4.186e7_f
      end if
    end if

    ! Properties for H2SO4
    if (igash2so4 /= 0) then
      ! Diffusivity
      rhoa_cgs = rhoa(k) / (xmet(k) * ymet(k) * zmet(k))
      aden     = rhoa_cgs * AVG / WTMOL_AIR
      diffus(k,igash2so4) = 1.76575e+17_f * sqrt(t(k)) / aden
  
      ! HACK: make H2SO4 latent heats same as water
      rlhe(k,igash2so4) = rlhe(k, igash2o)
      rlhm(k,igash2so4) = rlhe(k, igash2o)
    end if
    
  enddo

#ifdef DEBUG
  ! Report some initialization values
  if (do_print_init) then
    write(LUNOPRT,5)
    write(LUNOPRT,2) 'NGAS    ',NGAS
    write(LUNOPRT,1) 'igrowgas',(igrowgas(i),i=1,NELEM)
    write(LUNOPRT,3) (i,gwtmol(i),gasname(i),i=1,NGAS)
  endif
#endif

  ! Return to caller with particle growth mapping arrays and time-dependent
  ! parameters initialized.
  return
end
