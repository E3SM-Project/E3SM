! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

module wetr

contains

  !! This routine calculates the wet radius for hydrophilic particles that are
  !! assumed to grow in size based upon the realtive humidity.
  !!
  !! Parameterizations based upon Fitzgerald [1975] and Gerber [1985] are support and the
  !! particles are assumed to be spherical.
  !!
  !! @author  Chuck Bardeen, Pete Colarco
  !! @version May-2009 from Nov-2000 
  subroutine getwetr(carma, igroup, rh, rdry, rwet, rhopdry, rhopwet, rc, h2o_mass, h2o_vp, temp)
  
    ! types
    use carma_precision_mod
    use carma_enums_mod
    use carma_constants_mod
    use carma_types_mod
    use carmastate_mod
    use carma_mod
    use sulfate_utils
  
    implicit none
  
    type(carma_type), intent(in)         :: carma   !! the carma object
    integer, intent(in)                  :: igroup  !! group index
    real(kind=f), intent(in)             :: rh      !! relative humidity
    real(kind=f), intent(in)             :: rdry    !! dry radius [cm]
    real(kind=f), intent(out)            :: rwet    !! wet radius [cm]
    real(kind=f), intent(in)             :: rhopdry !! dry radius [cm]
    real(kind=f), intent(out)            :: rhopwet !! wet radius [cm]
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    real(kind=f), intent(in), optional   :: h2o_mass!! water vapor mass concentration (g/cm3)
    real(kind=f), intent(in), optional   :: h2o_vp  !! water eq. vaper pressure (dynes/cm2)    
    real(kind=f), intent(in), optional   :: temp    !! temperature [K]
  
    ! Local declarations
    real(kind=f)            :: humidity
    real(kind=f)            :: r_ratio
    real(kind=f)            :: wtpkelv, den1, den2, drho_dwt
    real(kind=f)            :: sigkelv, sig1, sig2, dsigma_dwt
    real(kind=f)            :: rkelvinH2O_a, rkelvinH2O_b, rkelvinH2O, h2o_kelv
        
    ! The following parameters relate to the swelling of seasalt like particles
    ! following Fitzgerald, Journal of Applied Meteorology, [1975].
    !
    ! Question - Should epsilon be 1._f? It means alpharat is 1 by definition.
    real(kind=f), parameter :: epsilon_  = 1._f     ! soluble fraction of deliquescing particle
    real(kind=f)            :: alphaComp
    real(kind=f)            :: alpha
    real(kind=f)            :: alpha1
    real(kind=f)            :: alpharat
    real(kind=f)            :: beta
    real(kind=f)            :: theta
    real(kind=f)            :: f1
    real(kind=f)            :: f2
  
    ! Parameters from Gerber [1985]
    real(kind=f)            :: c1
    real(kind=f)            :: c2
    real(kind=f)            :: c3
    real(kind=f)            :: c4
                               
    ! Define formats
    1 format(/,'Non-spherical particles specified for group ',i3, &
        ' (ishape=',i3,') but spheres assumed in wetr.f.'/)
  
    ! If humidty affects the particle, then determine the equilbirium
    ! radius and density based upon the relative humidity.
    if (irhswell(igroup) == I_NO_SWELLING) then
     
      ! No swelling, just use the dry values.
      rwet     = rdry
      rhopwet  = rhopdry
    else
    
      !  Warning message for non-spherical particles!
      if( ishape(igroup) .ne. I_SPHERE )then
        if (do_print) write(LUNOPRT,1) igroup, ishape(igroup)
        rc = RC_ERROR
        return
      endif
      
      ! The Parameterizations don't handly relative humidities of 0, and
      ! behave poorly when RH > 0.995, so cap the relative humidity
      ! used to these values.
      humidity = min(max(rh,tiny(1.0_f)), 0.995_f)
      
      ! Fitzgerald Parameterization
      if (irhswell(igroup) == I_FITZGERALD) then
        
        ! Calculate the alpha and beta parameters for the wet particle
        !            relative to amonium sulfate
        beta = exp((0.00077_f * humidity) / (1.009_f - humidity))
        if (humidity .le. 0.97_f) then
          theta = 1.058_f
        else
          theta = 1.058_f - (0.0155_f * (humidity - 0.97_f)) / (1.02_f - humidity**1.4_f)
        endif
           
        alpha1 = 1.2_f * exp((0.066_f * humidity) / (theta - humidity))
        f1 = 10.2_f - 23.7_f * humidity + 14.5_f * humidity**2
        f2 = -6.7_f + 15.5_f * humidity - 9.2_f * humidity**2
        alpharat = 1._f - f1 * (1._f - epsilon_) - f2 * (1._f - epsilon_**2)
  
        ! Scale the size based on the composition of the particle.
        select case(irhswcomp(igroup))
  
          case (I_SWF_NH42SO4)
            alphaComp = 1.00_f
          
          case(I_SWF_NH4NO3)
            alphaComp = 1.06_f
          
          case(I_SWF_NANO3)
            alphaComp = 1.17_f
          
          case(I_SWF_NH4CL)
            alphaComp = 1.23_f
          
          case(I_SWF_CACL2)
            alphaComp = 1.29_f
          
          case(I_SWF_NABR)
            alphaComp = 1.32_f
          
          case(I_SWF_NACL)
            alphaComp = 1.35_f
          
          case(I_SWF_MGCL2)
            alphaComp = 1.41_f
          
          case(I_SWF_LICL)
            alphaComp = 1.54_f
          
          case default
            if (do_print) write(LUNOPRT,*) "wetr:: ERROR - Unknown composition type  (", irhswcomp(igroup), &
              ") for Fitzgerald."
            rc = RC_ERROR
            return
        end select
  
        alpha = alphaComp * (alpha1 * alpharat)
  
        ! Determine the wet radius.
        !
        ! NOTE: Fitgerald's equations assume r in [um], so scale the cgs units
        ! appropriately.
        rwet = (alpha * (rdry * 1e4_f)**beta) * (1e-4_f)
      
        ! Determine the wet density from the wet radius.
        r_ratio  = (rdry / rwet)**3
        rhopwet = r_ratio * rhopdry + (1._f - r_ratio) * RHO_W
      end if
      
      
      ! Gerber Paremeterization
      if (irhswell(igroup) == I_GERBER) then
  
        ! Scale the size based on the composition of the particle.
        select case(irhswcomp(igroup))
  
          case (I_SWG_NH42SO4)
            c1 = 0.4809_f
            c2 = 3.082_f
            c3 = 3.110e-11_f
            c4 = -1.428_f
          
          case(I_SWG_URBAN)
            c1 = 0.3926_f
            c2 = 3.101_f
            c3 = 4.190e-11_f
            c4 = -1.404_f
          
          case(I_SWG_RURAL)
            c1 = 0.2789_f
            c2 = 3.115_f
            c3 = 5.415e-11_f
            c4 = -1.399_f
          
          case(I_SWG_SEA_SALT)
            c1 = 0.7674_f
            c2 = 3.079_f
            c3 = 2.572e-11_f
            c4 = -1.424_f
                  
          case default
            if (do_print) write(LUNOPRT,*) "wetr:: ERROR - Unknown composition type  (", irhswcomp(igroup), &
              ") for Gerber."
            rc = RC_ERROR
            return
        end select
        
        rwet  = ((c1 * rdry**c2 / (c3 * rdry**c4 - log10(humidity))) + rdry**3)**(1._f / 3._f)
      
        ! Determine the wet density from the wet radius.
        r_ratio  = (rdry / rwet)**3
        rhopwet = r_ratio * rhopdry + (1._f - r_ratio) * RHO_W
      end if
    end if
    
    
    ! Sulfate Aerosol, using weight percent.
    if (irhswell(igroup) == I_WTPCT_H2SO4) then
    
      ! Adjust calculation for the Kelvin effect of H2O:
      wtpkelv = 80._f ! start with assumption of 80 wt % H2SO4 
      den1 = 2.00151_f - 0.000974043_f * temp ! density at 79 wt %
      den2 = 2.01703_f - 0.000988264_f * temp ! density at 80 wt %
      drho_dwt = den2-den1 ! change in density for change in 1 wt %
      
      sig1 = 79.3556_f - 0.0267212_f * temp ! surface tension at 79.432 wt %
      sig2 = 75.608_f  - 0.0269204_f * temp ! surface tension at 85.9195 wt %      
      dsigma_dwt = (sig2-sig1) / (85.9195_f - 79.432_f) ! change in density for change in 1 wt %
      sigkelv = sig1 + dsigma_dwt * (80.0_f - 79.432_f)
      
      rwet = rdry * (100._f * rhopdry / wtpkelv / den2)**(1._f / 3._f)

      rkelvinH2O_b = 1._f + wtpkelv * drho_dwt / den2 - 3._f * wtpkelv &
          * dsigma_dwt / (2._f*sigkelv)

      rkelvinH2O_a = 2._f * gwtmol(igash2so4) * sigkelv / (den1 * RGAS * temp * rwet)     

      rkelvinH2O = exp (rkelvinH2O_a*rkelvinH2O_b)
            
      h2o_kelv = h2o_mass / rkelvinH2O
      wtpkelv = wtpct_tabaz(carma, temp, h2o_kelv, h2o_vp, rc)
      rhopwet   = sulfate_density(carma, wtpkelv, temp, rc)
      rwet      = rdry * (100._f * rhopdry / wtpkelv / rhopwet)**(1._f / 3._f)   
    end if
  
    ! Return to caller with wet radius evaluated.
    return
  end subroutine
end module
