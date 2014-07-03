!! This module is part of the CARMA module and contains enumerations that are part of
!! the CARMA and CARMASTATE objects.
!!
!! @author Chuck Bardeen
!! @ version July-2009
module carma_enums_mod

  !--
  ! Index values of CARMA's flags.  In a given list, begin with 1
  ! (instead of 0) so that undefined flags will produce an error. 
  !
  ! For example:
  ! if( itype(ielem) .eq. I_INVOLATILE )then
  !
  ! If itype(ielem) hasn't been defined (and is still 0), we do not want
  ! to execute the statements that follow.
  
  !  Define values of flag used for vertical transport
  !  boundary conditions (ixxxbnd_pc)
  integer, public, parameter :: I_FIXED_CONC = 1    !! Fixed Concentration
  integer, public, parameter :: I_FLUX_SPEC  = 2    !! Flux Specification
  
  !  Define values of flag used for particle element
  !  type specification (itype).
  integer, public, parameter :: I_INVOLATILE = 1    !! Involatile particle
  integer, public, parameter :: I_VOLATILE   = 2    !! Volatile particle
  integer, public, parameter :: I_COREMASS   = 3    !! Core Mass
  integer, public, parameter :: I_VOLCORE    = 4    !! Voltile Core
  integer, public, parameter :: I_CORE2MOM   = 5    !! Core Mass - 2 Moments
    
  !!  Define values of flag used for nucleation process
  !!  specification (inucproc).
  !!
  !!  NOTE: Some of these can be used in combination, so for aerosol freezing this is treated
  !!  as a bit mask. When setting for one (or more) of the Aerosol freezing methods, use:
  !!    IAERFREEZE + I_AF_xxx + I_AF_yyy + ...
  integer, public, parameter :: I_AF_TABAZADEH_2000 = 1     !! Aerosol Freezing, Tabazadeh[2000]
  integer, public, parameter :: I_AF_KOOP_2000      = 2     !! Aerosol Freezing, Koop[2000]
  integer, public, parameter :: I_AF_MOHLER_2010    = 4     !! Aerosol Freezing, Mohler[2010]
  integer, public, parameter :: I_AF_MURRAY_2010    = 8     !! Glassy Aerosol Freezing, Murray[2010]
  integer, public, parameter :: I_DROPACT           = 256   !! Droplet Activation
  integer, public, parameter :: I_AERFREEZE         = 512   !! Aerosol Freezing
  integer, public, parameter :: I_DROPFREEZE        = 1024  !! Droplet Freezing
  integer, public, parameter :: I_ICEMELT           = 2048  !! Ice Melting
  integer, public, parameter :: I_HETNUC            = 4096  !! Heterogeneous Nucleation
  integer, public, parameter :: I_HOMNUC            = 8192  !! Binary homogeneous gas-to-particle nucleation
  
  !  Define values of flag used for collection process (icollec)
  integer, public, parameter :: I_COLLEC_CONST = 1   !! Constant Collection Efficiency
  integer, public, parameter :: I_COLLEC_FUCHS = 2   !! Binwise Maxima of Fuchs' and Langmuir's Efficiencies
  integer, public, parameter :: I_COLLEC_DATA  = 3   !! Input Data
  
  !  Define values of flag used for coagulation operation (icoagop)
  integer, public, parameter :: I_COAGOP_CONST = 1   !! Constant Coagulation Kernel
  integer, public, parameter :: I_COAGOP_CALC  = 2   !! Calculate Coagulation Kernel

  !  Define values of flag used for particle shape (ishape)
  integer, public, parameter :: I_SPHERE   = 1   !! spherical
  integer, public, parameter :: I_HEXAGON  = 2   !! hexagonal prisms or plates
  integer, public, parameter :: I_CYLINDER = 3   !! circular disks, cylinders, or spheroids
  
  !  Define values of flag used for particle swelling parameterization (irhswell)
  integer, public, parameter :: I_NO_SWELLING  = 0   !! No swelling
  integer, public, parameter :: I_FITZGERALD   = 1   !! Fitzgerald
  integer, public, parameter :: I_GERBER       = 2   !! Gerber
  integer, public, parameter :: I_WTPCT_H2SO4  = 3   !! The weight percent method for sulfate aerosol
  
  !  Define vallues of flag used for particle swelling composition (Fiztgerald)
  integer, public, parameter :: I_SWF_NH42SO4   = 1   !! (NH4)2SO4
  integer, public, parameter :: I_SWF_NH4NO3    = 2   !! NH4NO3
  integer, public, parameter :: I_SWF_NANO3     = 3   !! NaNO3
  integer, public, parameter :: I_SWF_NH4CL     = 4   !! NH4Cl
  integer, public, parameter :: I_SWF_CACL2     = 5   !! CaCl2
  integer, public, parameter :: I_SWF_NABR      = 6   !! NaBr
  integer, public, parameter :: I_SWF_NACL      = 7   !! NaCl
  integer, public, parameter :: I_SWF_MGCL2     = 8   !! MgCl2
  integer, public, parameter :: I_SWF_LICL      = 9   !! LiCl

  !  Define vallues of flag used for particle swelling composition (Gerber)
  integer, public, parameter :: I_SWG_NH42SO4   = 11  !! (NH4)2SO4
  integer, public, parameter :: I_SWG_SEA_SALT  = 12  !! Sea Salt
  integer, public, parameter :: I_SWG_URBAN     = 13  !! Urban
  integer, public, parameter :: I_SWG_RURAL     = 14  !! Rural
  
  ! Routines to calculate gas vapor pressures
  integer, public, parameter :: I_VAPRTN_H2O_BUCK1981      = 1   !! H2O, Buck[1981]
  integer, public, parameter :: I_VAPRTN_H2O_MURPHY2005    = 2   !! H2O, Murphy & Koop [2005]
  integer, public, parameter :: I_VAPRTN_H2O_GOFF1946      = 3   !! H2O, Goff & Gratch [1946], used in CAM
  integer, public, parameter :: I_VAPRTN_H2SO4_AYERS1980   = 4   !! H2SO4, Ayers [1980] & Kumala [1990]

  ! Routines to calculate fall velocities
  integer, public, parameter :: I_FALLRTN_STD              = 1   !! Standard CARMA 2.3 routine (spherical only)
  integer, public, parameter :: I_FALLRTN_STD_SHAPE        = 2   !! Optional CARMA 2.3 routine (supports shapes)
  integer, public, parameter :: I_FALLRTN_HEYMSFIELD2010   = 3   !! Heymsfield & Westbrook [2010] (ice only)

  ! Routines to calculate mie optical properties
  integer, public, parameter :: I_MIERTN_TOON1981      = 1   !! Shell/Core, Toon & Ackerman [1981]
  integer, public, parameter :: I_MIERTN_BOHREN1983    = 2   !! Homogeneous Sphere, Bohren & Huffman [1983]
 
  ! Gas Composition  
  integer, public, parameter :: I_GCOMP_H2O             = 1   !! Water Vapor
  integer, public, parameter :: I_GCOMP_H2SO4           = 2   !! Sulphuric Acid
  integer, public, parameter :: I_GCOMP_SO2             = 3   !! Sulfer Dioxide
  
  ! How is the CARMA group represented in the parent model
  integer, public, parameter :: I_CNSTTYPE_PROGNOSTIC   = 1   !! Prognostic, advected constituent for each bin
  integer, public, parameter :: I_CNSTTYPE_DIAGNOSTIC   = 2   !! Diagnostic, bins diagonosed from model state
  
  ! Return Codes
  !
  ! NOTE: Also see error handling macros in globaer.h.
  integer, public, parameter :: RC_OK             = 0   !! Success
  integer, public, parameter :: RC_ERROR          = -1  !! Failure
  integer, public, parameter :: RC_WARNING        = 1   !! Warning
  integer, public, parameter :: RC_WARNING_RETRY  = 2   !! Warning, Retry Suggested


  !  Define values of symbols used to specify horizontal & vertical grid type.
  !   Grid selection is made by defining each of the variables
  !   <igridv> and <igridh> to one of the grid types known to the model.
  !
  !   Possible values for igridv:
  !       I_CART    cartesian
  !       I_SIG     sigma
  !       I_HYBRID  hybrid
  !
  !    Possible values for igridh:
  !       I_CART   cartesian
  !       I_LL     longitude_latitude
  !       I_LC     lambert_conformal
  !       I_PS     polar_stereographic
  !       I_ME     mercator
  integer, public, parameter   :: I_CART       = 1   !! Cartesian
  integer, public, parameter   :: I_SIG        = 2   !! Sigma
  integer, public, parameter   :: I_LL         = 3   !! Longitude & Latitude
  integer, public, parameter   :: I_LC         = 4   !! Lambert Conformal
  integer, public, parameter   :: I_PS         = 5   !! Polar Sterographic
  integer, public, parameter   :: I_ME         = 6   !! Mercator
  integer, public, parameter   :: I_HYBRID     = 7   !! Hybrid
end module

