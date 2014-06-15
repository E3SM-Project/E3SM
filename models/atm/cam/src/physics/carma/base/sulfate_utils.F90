! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

module sulfate_utils
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
 
  implicit none
  
  ! Declare the public methods.
  public wtpct_tabaz
  public sulfate_density
  public sulfate_surf_tens
  
  real(kind=f), public:: dnwtp(46), dnc0(46), dnc1(46)
   
  data dnwtp / 0._f, 1._f, 5._f, 10._f, 20._f, 25._f, 30._f, 35._f, 40._f, &
     41._f, 45._f, 50._f, 53._f, 55._f, 56._f, 60._f, 65._f, 66._f, 70._f, &
     72._f, 73._f, 74._f, 75._f, 76._f, 78._f, 79._f, 80._f, 81._f, 82._f, &
     83._f, 84._f, 85._f, 86._f, 87._f, 88._f, 89._f, 90._f, 91._f, 92._f, &
     93._f, 94._f, 95._f, 96._f, 97._f, 98._f, 100._f /
     
   data dnc0 / 1._f, 1.13185_f, 1.17171_f, 1.22164_f, 1.3219_f, 1.37209_f,       &
     1.42185_f, 1.4705_f, 1.51767_f, 1.52731_f, 1.56584_f, 1.61834_f, 1.65191_f, &
     1.6752_f, 1.68708_f, 1.7356_f, 1.7997_f, 1.81271_f, 1.86696_f, 1.89491_f,   &
     1.9092_f, 1.92395_f, 1.93904_f, 1.95438_f, 1.98574_f, 2.00151_f, 2.01703_f, &
     2.03234_f, 2.04716_f, 2.06082_f, 2.07363_f, 2.08461_f, 2.09386_f, 2.10143_f,&
     2.10764_f, 2.11283_f, 2.11671_f, 2.11938_f, 2.12125_f, 2.1219_f, 2.12723_f, &
     2.12654_f, 2.12621_f, 2.12561_f, 2.12494_f, 2.12093_f /
     
   data dnc1 / 0._f,  -0.000435022_f, -0.000479481_f, -0.000531558_f, -0.000622448_f,&
     -0.000660866_f, -0.000693492_f, -0.000718251_f, -0.000732869_f, -0.000735755_f, &
     -0.000744294_f, -0.000761493_f, -0.000774238_f, -0.00078392_f, -0.000788939_f,  &
     -0.00080946_f, -0.000839848_f, -0.000845825_f, -0.000874337_f, -0.000890074_f,  &
     -0.00089873_f, -0.000908778_f, -0.000920012_f, -0.000932184_f, -0.000959514_f,  &
     -0.000974043_f, -0.000988264_f, -0.00100258_f, -0.00101634_f, -0.00102762_f,    &
     -0.00103757_f, -0.00104337_f, -0.00104563_f, -0.00104458_f, -0.00104144_f,      &
     -0.00103719_f, -0.00103089_f, -0.00102262_f, -0.00101355_f, -0.00100249_f,      &
     -0.00100934_f, -0.000998299_f, -0.000990961_f, -0.000985845_f, -0.000984529_f,  &
     -0.000989315_f /  
contains

  !!  This function calculates the weight % H2SO4 composition of 
  !!  sulfate aerosol, using Tabazadeh et. al. (GRL, 1931, 1997).
  !!  Rated for T=185-260K, activity=0.01-1.0
  !!
  !!  Argument list input:   
  !!    temp = temperature (K)
  !!    h2o_mass = water vapor mass concentration (g/cm3)
  !!    h2o_vp = water eq. vaper pressure (dynes/cm2)
  !!
  !!  Output:
  !!    wtpct_tabaz = weight % H2SO4 in H2O/H2SO4 particle (0-100)
  !!
  !!  Include global constants and variables (BK=Boltzman constant,
  !!   AVG=Avogadro's constant)
  !!
  !! @author Jason English
  !! @ version Apr-2010
  function wtpct_tabaz(carma, temp, h2o_mass, h2o_vp, rc)
   
    real(kind=f)                         :: wtpct_tabaz
    type(carma_type), intent(in)         :: carma     !! the carma object
    real(kind=f), intent(in)             :: temp      !! temperature [K]
    real(kind=f), intent(in)             :: h2o_mass  !! water vapor mass concentration (g/cm3)
    real(kind=f), intent(in)             :: h2o_vp    !! water eq. vaper pressure (dynes/cm2) 
    integer, intent(inout)               :: rc        !! return code, negative indicates failure
      
    !  Declare variables for this routine only
    real(kind=f)     :: atab1,btab1,ctab1,dtab1,atab2,btab2,ctab2,dtab2
    real(kind=f)     :: h2o_num, p_h2o, vp_h2o
    real(kind=f)     :: contl, conth, contt, conwtp
    real(kind=f)     :: activ 
         
    ! Get number density of water (/cm3) from mass concentration (g/cm3)
    h2o_num=h2o_mass*AVG/gwtmol(1)

    !  Get partial pressure of water (dynes/cm2) from concentration (/cm3)
    ! Ideal gas law: P=nkT
    p_h2o=h2o_num*bk*temp

    !  Convert from dynes/cm2 to mb (hPa)
    p_h2o=p_h2o/1000.0_f     ! partial pressure
    vp_h2o=h2o_vp/1000.0_f   ! eq. vp

    !  Prevent a NaN calculation  
    !  In the upper thermosphere p_h2o can be very low and vp_h2o can be very high
    if (p_h2o.lt.1.e-10_f .and. vp_h2o.gt.0._f) p_h2o=1.e-10_f   
   
    !  Activity = water pp in mb / water eq. vp over pure water in mb
    activ = p_h2o/vp_h2o
 
    if (activ.lt.0.05_f) then
      activ = max(activ,1.e-6_f)    ! restrict minimum activity
      atab1 	= 12.37208932_f	
      btab1 	= -0.16125516114_f
      ctab1 	= -30.490657554_f
      dtab1 	= -2.1133114241_f
      atab2 	= 13.455394705_f	
      btab2 	= -0.1921312255_f
      ctab2 	= -34.285174607_f
      dtab2 	= -1.7620073078_f
    elseif (activ.ge.0.05_f.and.activ.le.0.85_f) then
      atab1 	= 11.820654354_f
      btab1 	= -0.20786404244_f
      ctab1 	= -4.807306373_f
      dtab1 	= -5.1727540348_f
      atab2 	= 12.891938068_f	
      btab2 	= -0.23233847708_f
      ctab2 	= -6.4261237757_f
      dtab2 	= -4.9005471319_f
    elseif (activ.gt.0.85_f) then
      activ = min(activ,1._f)      ! restrict maximum activity
      atab1 	= -180.06541028_f
      btab1 	= -0.38601102592_f
      ctab1 	= -93.317846778_f
      dtab1 	= 273.88132245_f
      atab2 	= -176.95814097_f
      btab2 	= -0.36257048154_f
      ctab2 	= -90.469744201_f
      dtab2 	= 267.45509988_f
    else
      if (do_print) write(LUNOPRT,*) 'invalid activity: activity,pp,vp=',activ, p_h2o
      rc = RC_ERROR
      return
    endif

    contl = atab1*(activ**btab1)+ctab1*activ+dtab1
    conth = atab2*(activ**btab2)+ctab2*activ+dtab2
      
    contt = contl + (conth-contl) * ((temp -190._f)/70._f)
    conwtp = (contt*98._f) + 1000._f

    wtpct_tabaz = (100._f*contt*98._f)/conwtp
    wtpct_tabaz = min(max(wtpct_tabaz,1._f),100._f) ! restrict between 1 and 100 %
      
    !  Note: restricting activity to 1.e-6 minimum allows for a maximum of
    !  98.5 wtpct at T=650K, 95.8 wtpct at T=300K, and 90.9 wtpct at 180K.
  
    return
  end function wtpct_tabaz
           
  !! Calculates specific gravity (g/cm3) of sulfate of 
  !! different compositions as a linear function of temperature,
  !! based of measurements of H2SO4/H2O solution densities made 
  !! at 0 to 100C tabulated in the International Critical Tables 
  !! (Washburn, ed., NRC, 1928). Measurements have confirmed that 
  !! this data may be linearly extrapolated to stratospheric 
  !! temperatures (180-380K) with excellent accuracy 
  !! (Beyer, Ravishankara, & Lovejoy, JGR, 1996).
  !!
  !! Argument list input:
  !!    wtp = aerosol composition in weight % H2SO4 (0-100)
  !!    temp = temperature in Kelvin
  !!
  !! Output:
  !!    sulfate_density (g/cm3) [function name]
  !!
  !! This function requires setup_sulfate_density to be run
  !! first to read in the density coefficients DNC0 and DNC1
  !! and the tabulated weight percents DNWTP.
  !!
  !! @author Mike Mills
  !! @version Mar-2013
  function sulfate_density(carma, wtp, temp, rc)

  !! Include global constants and variables

    real(kind=f)                         :: sulfate_density
    type(carma_type), intent(in)         :: carma   !! the carma object
    real(kind=f), intent(in)             :: wtp     !! weight percent
    real(kind=f), intent(in)             :: temp    !! temperature 
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    
    ! Local declarations
    integer           :: i
    real(kind=f)      :: den1, den2
    real(kind=f)      :: frac, temp_loc

    if (wtp .lt. 0.0_f .or. wtp .gt. 100.0_f) then
      if (do_print) write(LUNOPRT,*)'sulfate_density: Illegal value for wtp:',wtp
      rc = RC_ERROR
      return
    endif

    ! limit temperature to bounds of extrapolation
    temp_loc=min(temp, 380.0_f)
    temp_loc=max(temp_loc, 180.0_f)

    i=1

    do while (wtp .gt. dnwtp(i))
     i=i+1
    end do

    den2=dnc0(i)+dnc1(i)*temp_loc

    if (i.eq.1 .or. wtp.eq.dnwtp(i)) then
      sulfate_density=den2
      return
    endif

    den1=dnc0(i-1)+dnc1(i-1)*temp_loc
    frac=(dnwtp(i)-wtp)/(dnwtp(i)-dnwtp(i-1))
    sulfate_density=den1*frac+den2*(1.0_f-frac)

    return
  end function sulfate_density

  !!  Calculates surface tension (erg/cm2 = dyne/cm) of sulfate of 
  !!  different compositions as a linear function of temperature,
  !!  as described in Mills (Ph.D. Thesis, 1996), derived from
  !!  the measurements of Sabinina and Terpugow (1935).
  !!
  !!  Argument list input:
  !!     WTP = aerosol composition in weight % H2SO4 (0-100)
  !!     TEMP = temperature in Kelvin
  !!
  !!  Output:
  !!     sulfate_surf_tens (erg/cm2) [function name]
  !!
  !!  This function requires setup_sulfate_density to be run
  !!  first to read in the density coefficients DNC0 and DNC1
  !!  and the tabulated weight percents DNWTP.
  !!
  !! @author Mike Mills
  !! @version Mar-2013
  function sulfate_surf_tens(carma, wtp, temp, rc)  
  
    real(kind=f)                         :: sulfate_surf_tens
    type(carma_type), intent(in)         :: carma   !! the carma object
    real(kind=f), intent(in)             :: wtp     !! weight percent
    real(kind=f), intent(in)             :: temp    !! temperature 
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    
    ! Local declarations
    integer           :: i  
    real(kind=f)      :: sig1, sig2
    real(kind=f)      :: frac, temp_loc
    real(kind=f)      :: stwtp(15), stc0(15), stc1(15)
    
    data stwtp/0._f, 23.8141_f, 38.0279_f, 40.6856_f, 45.335_f, 52.9305_f, 56.2735_f, &
       & 59.8557_f, 66.2364_f, 73.103_f, 79.432_f, 85.9195_f, 91.7444_f, 97.6687_f, 100._f/
    
    data stc0/117.564_f, 103.303_f, 101.796_f, 100.42_f, 98.4993_f, 91.8866_f,     &
       & 88.3033_f, 86.5546_f, 84.471_f, 81.2939_f, 79.3556_f, 75.608_f, 70.0777_f,  &
       & 63.7412_f, 61.4591_f /
  
    data stc1/-0.153641_f, -0.0982007_f, -0.0872379_f, -0.0818509_f,           &
       & -0.0746702_f, -0.0522399_f, -0.0407773_f, -0.0357946_f, -0.0317062_f,   &
       & -0.025825_f, -0.0267212_f, -0.0269204_f, -0.0276187_f, -0.0302094_f,    &
       & -0.0303081_f /

    ! limit temperature to reasonable bounds of extrapolation
    temp_loc=min(temp, 380.0_f)
    temp_loc=max(temp_loc, 180.0_f)
       
    if (wtp .lt. 0.0_f .OR. wtp .gt. 100.0_f) then
      if (do_print) write(LUNOPRT,*)'sulfate_surf_tens: Illegal value for wtp:',wtp
      if (do_print) write(LUNOPRT,*)'sulfate_surf_tens: temp=',temp
      rc = RC_ERROR
      return
    endif
  
    i=1
  
    do while (wtp.gt.stwtp(i))
      i=i+1
    end do
  
    sig2=stc0(i)+stc1(i)*temp_loc
  
    if (i.eq.1 .or. wtp.eq.stwtp(i)) then
      sulfate_surf_tens=sig2
      return
    end if
  
    sig1=stc0(i-1)+stc1(i-1)*temp_loc
    frac=(stwtp(i)-wtp)/(stwtp(i)-stwtp(i-1))
    sulfate_surf_tens=sig1*frac+sig2*(1.0_f-frac)
  
    return
  end function sulfate_surf_tens
            
end module sulfate_utils
