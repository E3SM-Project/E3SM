#ifndef OMEGA_SURFACE_FLUX_H
#define OMEGA_SURFACE_FLUX_H
//===-- ocn/SurfaceFlux.h - Surface Flux Helpers ----------*- C++ -*-===//
//
/// \file
/// \brief Shared surface heat and freshwater flux calculations
//
//===----------------------------------------------------------------------===//

#include "Eos.h"
#include "GlobalConstants.h"

namespace OMEGA {

/// Surface heat flux used by KPP buoyancy forcing, excluding mass-flux
/// enthalpy carried by the tracer tendency.
KOKKOS_INLINE_FUNCTION Real sfcHeatFluxWithoutMassEnthalpy(
    I4 ICell, Real SaTop, Real PTopDb, EosType EosChoice,
    const Array1DReal &LongWaveHeatFluxUp,
    const Array1DReal &LongWaveHeatFluxDown,
    const Array1DReal &ShortWaveHeatFlux, const Array1DReal &SensibleHeatFlux,
    const Array1DReal &SeaIceHeatFlux, const Array1DReal &SeaIceFreshWaterFlux,
    const Array1DReal &SeaIceSaltFlux, const Array1DReal &LatentHeatFluxEvap,
    const Array1DReal &SnowFlux, const Array1DReal &IceRunoffFlux) {

   const Real SeaIceFreshWater = SeaIceFreshWaterFlux(ICell);
   const Real SeaIceSalt       = SeaIceSaltFlux(ICell);
   const Real SaSeaIce = SeaIceFreshWater > 0.0_Real && SeaIceSalt > 0.0_Real
                             ? SeaIceSalt / SeaIceFreshWater * Salt2PPt
                             : SaTop;

   const Real SeaIceLiqEnthalpyEstimate =
       SeaIceFreshWaterFlux(ICell) * Cp0Sw *
       Eos::calcCtFreezing(EosChoice, SaSeaIce, PTopDb, 0.0_Real);

   return LongWaveHeatFluxUp(ICell) + LongWaveHeatFluxDown(ICell) +
          ShortWaveHeatFlux(ICell) + SensibleHeatFlux(ICell) +
          SeaIceHeatFlux(ICell) - SeaIceLiqEnthalpyEstimate +
          LatentHeatFluxEvap(ICell) +
          (SnowFlux(ICell) + IceRunoffFlux(ICell)) * -LatIce;
}

/// Net surface freshwater mass flux (kg/m^2/s).
KOKKOS_INLINE_FUNCTION Real sfcFreshWaterFlux(
    I4 ICell, const Array1DReal &SnowFlux, const Array1DReal &RainFlux,
    const Array1DReal &EvaporationFlux, const Array1DReal &SeaIceFreshWaterFlux,
    const Array1DReal &IceRunoffFlux, const Array1DReal &RiverRunoffFlux) {
   return SnowFlux(ICell) + RainFlux(ICell) + EvaporationFlux(ICell) +
          SeaIceFreshWaterFlux(ICell) + IceRunoffFlux(ICell) +
          RiverRunoffFlux(ICell);
}

} // namespace OMEGA

#endif // OMEGA_SURFACE_FLUX_H
