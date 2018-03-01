module ThermalKSPTemperatureSoilAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !
  ! !USES:
  use ThermalKSPTemperatureBaseAuxType
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none

  private

  type, public, extends(therm_ksp_temp_base_auxvar_type)  :: therm_ksp_temp_soil_auxvar_type

     PetscReal :: liq_areal_den         ! [kg/m^2] = h2osoi_liq
     PetscReal :: ice_areal_den         ! [kg/m^2] = h2osoi_ice
     PetscReal :: snow_water            ! snow water (mm H2O)
     PetscInt  :: num_snow_layer        ! number of snow layer
     PetscReal :: tuning_factor         !
     PetscReal :: dz                    !

     ! parameters
     PetscReal :: por                   ! [m^3/m^3]
     PetscReal :: therm_cond_minerals   ! [W/(m K)]
     PetscReal :: therm_cond_dry        ! [W/(m K)]
     PetscReal :: heat_cap_minerals_puv ! [J/(m3 K)]
     PetscBool :: is_soil_shallow
     PetscInt  :: itype

   contains
     procedure, public :: Init                  => ThermKSPTempSoilAuxVarInit
     procedure, public :: AuxVarCompute         => ThermKSPTempSoilAuxVarCompute
  end type therm_ksp_temp_soil_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_soil_auxvar_type)   :: this

    call this%BaseInit()

    this%liq_areal_den         = 0.d0
    this%ice_areal_den         = 0.d0
    this%snow_water            = 0.d0
    this%num_snow_layer        = 0
    this%tuning_factor         = 1.d0
    this%dz                    = 0.d0

    this%por                   = 0.d0
    this%therm_cond_minerals   = 0.d0
    this%therm_cond_dry        = 0.d0
    this%heat_cap_minerals_puv = 0.d0

    this%is_soil_shallow       = PETSC_FALSE
    this%itype                 = -1

  end subroutine ThermKSPTempSoilAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilAuxVarCompute(this, dz, vol)
    !
    ! !DESCRIPTION:
    !
    use mpp_varcon , only : denh2o, denice, tfrz, tkwat, tkice, cpice,  cpliq, thk_bedrock
    use mpp_varcon , only : istice, istice_mec, istwet, istsoil, istcrop
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_soil_auxvar_type)   :: this
    PetscReal                                :: dz
    PetscReal                                :: vol
    !
    ! LOCAL VARIABLES
    PetscReal :: satw
    PetscReal :: fl
    PetscReal :: dke
    PetscReal :: dksat

    !select case(this%itype)
    !case (istsoil, istcrop)
    if (this%itype == istsoil .or. this%itype == istcrop) then
       if (this%is_soil_shallow) then

          satw = (this%liq_areal_den/denh2o + this%ice_areal_den/denice)/(dz*this%por)
          satw = min(1.d0, satw)

          if (satw > .1e-6) then
             if (this%temperature >= tfrz) then       ! Unfrozen soil
                dke = max(0.d0, log10(satw) + 1.0d0)
             else                                      ! Frozen soil
                dke = satw
             end if

             fl = (this%liq_areal_den/(denh2o*dz)) / (this%liq_areal_den/(denh2o*dz) + &
                                                      this%ice_areal_den/(denice*dz))

             dksat = this%therm_cond_minerals*    &
                  tkwat**(fl*this%por)*        &
                  tkice**((1.d0-fl)*this%por)

             this%therm_cond = dke*dksat + (1.d0-dke)*this%therm_cond_dry

          else
             this%therm_cond     = this%therm_cond_dry
          endif

          this%heat_cap_pva   = this%heat_cap_minerals_puv*(1.d0 - this%por)*dz + &
                                this%ice_areal_den*cpice                         + &
                                this%liq_areal_den*cpliq

          if (this%num_snow_layer == 0) then
             this%heat_cap_pva = this%heat_cap_pva + this%snow_water*cpice
          endif

       else

          this%therm_cond     = thk_bedrock
          this%heat_cap_pva   = this%heat_cap_minerals_puv*(1.d0 - this%por)*dz + &
                                this%ice_areal_den*cpice                         + &
                                this%liq_areal_den*cpliq
       endif

       this%heat_cap_pva = this%heat_cap_pva/dz

    else if (this%itype == istwet) then

       if (this%is_soil_shallow) then
          if (this%temperature < tfrz) then
             this%therm_cond   = tkice
          else
             this%therm_cond = tkwat
          endif
          this%heat_cap_pva = this%ice_areal_den*cpice + this%liq_areal_den*cpliq
          if (this%num_snow_layer == 0) then
             this%heat_cap_pva = this%heat_cap_pva + this%snow_water*cpice
          endif

          this%heat_cap_pva = this%heat_cap_pva/dz
       else
          this%therm_cond   = thk_bedrock
          this%heat_cap_pva = this%heat_cap_minerals_puv
       endif

    else if (this%itype == istice .or. this%itype == istice_mec) then
          if (this%temperature < tfrz) then
             this%therm_cond   = tkice
          else
             this%therm_cond = tkwat
          endif
          this%heat_cap_pva = this%ice_areal_den*cpice + this%liq_areal_den*cpliq
          if (this%num_snow_layer == 0) then
             this%heat_cap_pva = this%heat_cap_pva + this%snow_water*cpice
          endif

          this%heat_cap_pva = this%heat_cap_pva/dz
    end if
       

  end subroutine ThermKSPTempSoilAuxVarCompute

#endif
  
end module ThermalKSPTemperatureSoilAuxType
