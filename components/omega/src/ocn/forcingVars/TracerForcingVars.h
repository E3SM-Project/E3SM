#ifndef OMEGA_TRACER_FORCING_H
#define OMEGA_TRACER_FORCING_H

#include "DataTypes.h"
#include "HorzMesh.h"

#include <string>

namespace OMEGA {

// Forward declarations. Full definitions not needed in this header since only
// pointers are used.
class VertCoord;
class Eos;

class TracerForcingVars {
 public:
   Array1DReal SnowFluxCell;
   Array1DReal RainFluxCell;
   Array1DReal EvaporationFluxCell;
   Array1DReal SeaIceFreshWaterFluxCell;
   Array1DReal IceRunoffFluxCell;
   Array1DReal RiverRunoffFluxCell;

   Array1DReal LatentHeatFluxEvapCell;
   Array1DReal SensibleHeatFluxCell;
   Array1DReal LongWaveHeatFluxUpCell;
   Array1DReal LongWaveHeatFluxDownCell;
   Array1DReal SeaIceHeatFluxCell;
   Array1DReal ShortWaveHeatFluxCell;

   Array1DReal SeaIceSaltFluxCell;
   Array1DReal ExtraEnergyCell;

   /// KPP non-local surface tracer fluxes [NTracers, NCellsSize]. External
   /// forcing supplies non-temperature/salinity tracers; KPP derives and
   /// overwrites temperature and salinity from the coupled heat and mass
   /// fluxes.
   Array2DReal SurfaceTracerFluxCell;

   TracerForcingVars(const std::string &Suffix, const HorzMesh *Mesh);

   void registerFields(const std::string &MeshName) const;
   void unregisterFields() const;
};

} // namespace OMEGA

#endif
